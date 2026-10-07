// Oracle application for the electron code-to-code comparison of lindhard
// (docs/validation.md, "Electron oracles: Nebula and Geant4 MicroElec").
//
// Our own code, MIT-licensed like the rest of this tree. It is compiled
// against a Geant4 installation built outside this tree and contains no
// Geant4 source. Written from the Geant4 Book For Application Developers,
// release 11.4 (https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/):
// "How to Define the main() Program" (run manager, the three user
// initialization classes), "How to Define a Detector Geometry" and "How to
// Specify Materials in the Detector" (NIST materials such as G4_Si),
// "Optional User Actions" (run, event and stepping actions) and "Physics
// List" / "Electromagnetic physics" (G4EmStandardPhysics_option4, and
// G4EmParameters::AddMicroElec(region), which the guide lists among the
// G4EmParameters methods). The Physics Reference Manual 11.4, section "The
// MicroElec extension for microelectronics applications", states that the
// MicroElec models are valid for silicon only (G4_Si), that electrons are
// tracked step by step down to the eV scale, and that electrons below 16.7 eV
// are killed and deposit their energy locally.
//
// Geometry: a slab of the chosen NIST material occupying 0 <= z <= slab
// thickness and |x|, |y| <= half width, in a G4_Galactic world. The slab is
// its own region ("Target") with MicroElec active in it. Electrons start at
// z = -1 nm, travelling along +z, at the beam energy.
//
// Tallies (per worker thread, written at the end of the run into --out):
//   events_<thread>.txt   one line per event:
//       event fast slow primary_stopped primary_depth_nm
//     fast / slow: electrons leaving through the front face (z = 0) with
//     kinetic energy >= / < the split (50 eV, the convention of
//     lindhard::tally::SE_BSE_SPLIT_EV); primary_stopped is 1 if the primary
//     (parent ID 0) came to rest inside the slab, and primary_depth_nm its
//     depth z there (otherwise nan).
//   deposit_<thread>.txt  one line per non-empty cell:
//       batch bin energy_ev
//     energy deposited in the slab (G4Step::GetTotalEnergyDeposit, placed at
//     the post-step point), by radial distance from the beam axis in
//     --rbins linear bins over [0, --rmax-nm) plus one overflow bin
//     (bin = --rbins), and by batch = event ID modulo --batches.
//   run_<thread>.txt      events tallied and energy deposited by this thread.
//
// Usage:
//   lindhard_microelec --energy-ev 5000 --material G4_Si --events 10000 \
//     --threads 8 --seed 1 --out DIR --rmax-nm 5000 --rbins 5000 \
//     --batches 20 [--slab-nm 100000] [--half-width-nm 100000]

#include "G4Box.hh"
#include "G4EmParameters.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Region.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {

struct Options {
  double energy_ev = 0.0;
  std::string material = "G4_Si";
  long events = 0;
  int threads = 1;
  long seed = 1;
  std::string out;
  double rmax_nm = 0.0;
  int rbins = 0;
  int batches = 1;
  double slab_nm = 100000.0;
  double half_width_nm = 100000.0;
};

// SE/BSE split, eV: the default of lindhard::tally::SE_BSE_SPLIT_EV (Chen et
// al., Sci. Rep. 12, 18201 (2022)); E >= split is backscattered.
constexpr double kSplitEv = 50.0;

Options parse(int argc, char** argv) {
  Options o;
  for (int i = 1; i + 1 < argc; i += 2) {
    const std::string k = argv[i];
    const std::string v = argv[i + 1];
    if (k == "--energy-ev") o.energy_ev = std::stod(v);
    else if (k == "--material") o.material = v;
    else if (k == "--events") o.events = std::stol(v);
    else if (k == "--threads") o.threads = std::stoi(v);
    else if (k == "--seed") o.seed = std::stol(v);
    else if (k == "--out") o.out = v;
    else if (k == "--rmax-nm") o.rmax_nm = std::stod(v);
    else if (k == "--rbins") o.rbins = std::stoi(v);
    else if (k == "--batches") o.batches = std::stoi(v);
    else if (k == "--slab-nm") o.slab_nm = std::stod(v);
    else if (k == "--half-width-nm") o.half_width_nm = std::stod(v);
    else {
      std::cerr << "unknown option " << k << "\n";
      std::exit(2);
    }
  }
  if (o.energy_ev <= 0 || o.events <= 0 || o.out.empty() || o.rmax_nm <= 0 || o.rbins <= 0 ||
      o.batches <= 0 || o.threads <= 0) {
    std::cerr << "missing or invalid options; see the header of lindhard_microelec.cc\n";
    std::exit(2);
  }
  return o;
}

class Detector : public G4VUserDetectorConstruction {
 public:
  explicit Detector(const Options& o) : o_(o) {}

  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    G4Material* target = nist->FindOrBuildMaterial(o_.material);
    if (target == nullptr) {
      std::cerr << "unknown NIST material " << o_.material << "\n";
      std::exit(2);
    }
    const double hw = o_.half_width_nm * nm;
    const double t = o_.slab_nm * nm;
    // World: z from -t to t; the slab fills 0 <= z <= t.
    auto* world_s = new G4Box("World", 1.01 * hw, 1.01 * hw, t);
    auto* world_l = new G4LogicalVolume(world_s, vacuum, "World");
    auto* world_p = new G4PVPlacement(nullptr, G4ThreeVector(), world_l, "World", nullptr, false, 0);
    auto* slab_s = new G4Box("Target", hw, hw, 0.5 * t);
    auto* slab_l = new G4LogicalVolume(slab_s, target, "Target");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0.5 * t), slab_l, "Target", world_l, false, 0);
    auto* region = new G4Region("Target");
    region->AddRootLogicalVolume(slab_l);
    return world_p;
  }

 private:
  Options o_;
};

class Physics : public G4VModularPhysicsList {
 public:
  Physics() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    // The MicroElec models are discrete; production cuts do not apply to
    // them. A short cut keeps the standard processes outside the region
    // (vacuum only here) from mattering.
    SetDefaultCutValue(1 * nm);
  }
};

// Per-thread accumulators, owned by the run action.
struct Tally {
  explicit Tally(const Options& o)
      : o(o), deposit(static_cast<size_t>(o.batches) * (o.rbins + 1), 0.0) {}
  const Options& o;
  std::vector<double> deposit;
  std::vector<std::string> lines;
  long events = 0;
  double deposited_ev = 0.0;
  // Current event.
  long fast = 0;
  long slow = 0;
  bool primary_stopped = false;
  double primary_depth_nm = std::numeric_limits<double>::quiet_NaN();
  int batch = 0;
};

class RunAction : public G4UserRunAction {
 public:
  explicit RunAction(Tally* t) : t_(t) {}
  void EndOfRunAction(const G4Run*) override {
    if (IsMaster()) return;
    const int id = G4Threading::G4GetThreadId();
    const std::string sfx = std::to_string(id < 0 ? 0 : id) + ".txt";
    std::ofstream ev(t_->o.out + "/events_" + sfx);
    ev << "# event fast slow primary_stopped primary_depth_nm\n";
    for (const auto& l : t_->lines) ev << l << "\n";
    std::ofstream dep(t_->o.out + "/deposit_" + sfx);
    dep.precision(17);
    dep << "# batch bin energy_ev (bins: " << t_->o.rbins << " over [0, " << t_->o.rmax_nm
        << ") nm, then overflow)\n";
    const int nb = t_->o.rbins + 1;
    for (int b = 0; b < t_->o.batches; ++b)
      for (int i = 0; i < nb; ++i) {
        const double e = t_->deposit[static_cast<size_t>(b) * nb + i];
        if (e != 0.0) dep << b << " " << i << " " << e << "\n";
      }
    std::ofstream run(t_->o.out + "/run_" + sfx);
    run.precision(17);
    run << "events " << t_->events << "\ndeposited_ev " << t_->deposited_ev << "\n";
  }

 private:
  Tally* t_;
};

class EventAction : public G4UserEventAction {
 public:
  explicit EventAction(Tally* t) : t_(t) {}
  void BeginOfEventAction(const G4Event* e) override {
    t_->fast = 0;
    t_->slow = 0;
    t_->primary_stopped = false;
    t_->primary_depth_nm = std::numeric_limits<double>::quiet_NaN();
    t_->batch = static_cast<int>(e->GetEventID() % t_->o.batches);
  }
  void EndOfEventAction(const G4Event* e) override {
    t_->events += 1;
    char buf[160];
    std::snprintf(buf, sizeof buf, "%d %ld %ld %d %.9g", e->GetEventID(), t_->fast, t_->slow,
                  t_->primary_stopped ? 1 : 0, t_->primary_depth_nm);
    t_->lines.emplace_back(buf);
  }

 private:
  Tally* t_;
};

class SteppingAction : public G4UserSteppingAction {
 public:
  explicit SteppingAction(Tally* t) : t_(t) {}
  void UserSteppingAction(const G4Step* step) override {
    const G4StepPoint* pre = step->GetPreStepPoint();
    const G4StepPoint* post = step->GetPostStepPoint();
    const G4VPhysicalVolume* pre_v = pre->GetPhysicalVolume();
    if (pre_v == nullptr || pre_v->GetName() != "Target") return;
    const G4ThreeVector p = post->GetPosition();
    const double edep = step->GetTotalEnergyDeposit();
    if (edep > 0) {
      const double e_ev = edep / eV;
      const double r_nm = std::hypot(p.x(), p.y()) / nm;
      int bin = static_cast<int>(r_nm / t_->o.rmax_nm * t_->o.rbins);
      if (bin < 0 || bin >= t_->o.rbins) bin = t_->o.rbins;
      t_->deposit[static_cast<size_t>(t_->batch) * (t_->o.rbins + 1) + bin] += e_ev;
      t_->deposited_ev += e_ev;
    }
    G4Track* track = step->GetTrack();
    const G4VPhysicalVolume* post_v = post->GetPhysicalVolume();
    if (post_v != pre_v && p.z() <= 0.5 * nm) {
      // Left through the front face (z = 0): count it and stop following it.
      if (track->GetDefinition()->GetParticleName() == "e-") {
        const double ke_ev = post->GetKineticEnergy() / eV;
        if (ke_ev >= kSplitEv) t_->fast += 1;
        else t_->slow += 1;
      }
      track->SetTrackStatus(fStopAndKill);
      return;
    }
    // The primary came to rest in the slab: either its kinetic energy reached
    // zero, or a process stopped it (the MicroElec kill below its tracking
    // limit, Physics Reference Manual 11.4, sets the track to stop with its
    // remaining energy deposited locally). The track status is one of the
    // G4Track quantities listed in the Book For Application Developers 11.4,
    // "Tracking", "Access to Track and Step Information".
    const G4TrackStatus status = track->GetTrackStatus();
    const bool at_rest = post->GetKineticEnergy() <= 0.0 || status == fStopAndKill ||
                         status == fStopButAlive;
    if (track->GetParentID() == 0 && post_v == pre_v && at_rest) {
      t_->primary_stopped = true;
      t_->primary_depth_nm = p.z() / nm;
    }
  }

 private:
  Tally* t_;
};

class Primary : public G4VUserPrimaryGeneratorAction {
 public:
  explicit Primary(const Options& o) : gun_(1) {
    gun_.SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("e-"));
    gun_.SetParticleEnergy(o.energy_ev * eV);
    gun_.SetParticlePosition(G4ThreeVector(0, 0, -1 * nm));
    gun_.SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
  }
  void GeneratePrimaries(G4Event* e) override { gun_.GeneratePrimaryVertex(e); }

 private:
  G4ParticleGun gun_;
};

class Actions : public G4VUserActionInitialization {
 public:
  explicit Actions(const Options& o) : o_(o) {}
  void Build() const override {
    auto* t = new Tally(o_);
    SetUserAction(new Primary(o_));
    SetUserAction(new RunAction(t));
    SetUserAction(new EventAction(t));
    SetUserAction(new SteppingAction(t));
  }

 private:
  Options o_;
};

}  // namespace

int main(int argc, char** argv) {
  const Options o = parse(argc, argv);
  G4Random::setTheSeed(o.seed);
  auto* rm = G4RunManagerFactory::CreateRunManager();
  rm->SetNumberOfThreads(o.threads);
  rm->SetUserInitialization(new Detector(o));
  // MicroElec in the slab's region (Book For Application Developers,
  // "Electromagnetic physics", G4EmParameters::AddMicroElec(region)).
  G4EmParameters::Instance()->AddMicroElec("Target");
  rm->SetUserInitialization(new Physics());
  rm->SetUserInitialization(new Actions(o));
  rm->Initialize();
  rm->BeamOn(static_cast<G4int>(o.events));
  delete rm;
  return 0;
}
