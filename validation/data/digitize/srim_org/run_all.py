"""Digitize every hydrogen and helium stopping plot of a srim.org mirror and
write the compilation CSV.

Usage: python3 -I run_all.py <mirror-site-dir> <work-dir> <out.csv> [ion_z ...]

<mirror-site-dir> holds StopPlots/STOPzzZZ.gif and Citations/Cit zz-ZZ.pdf.
One JSON per plot is written to <work-dir> (kept for review; not committed),
then every attributable point goes to <out.csv>. Default ions: 1 and 2.
"""
import csv, glob, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import citpdf
import digitize

METHOD = "srim_org-digitize/1"
SYM = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
       "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb "
       "Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U").split()
COLUMNS = ["ion_z", "ion", "target_z", "target", "energy_kev_per_amu", "stopping_ev_cm2_per_1e15_atoms",
           "u_energy_rel", "u_stopping", "citation_keys", "old_refs", "attribution", "label", "plot", "colour",
           "pixel_x", "pixel_y", "origin", "method"]


def run(site, work, ions):
    os.makedirs(work, exist_ok=True)
    for f in sorted(os.listdir(os.path.join(site, "StopPlots"))):
        if not re.match(r"STOP\d{4}\.gif$", f) or int(f[4:6]) not in ions:
            continue
        out = os.path.join(work, f[:-4] + ".json")
        if os.path.exists(out):
            continue
        pdf = os.path.join(site, "Citations", f"Cit {f[4:6]}-{f[6:8]}.pdf")
        try:
            o = digitize.digitize(os.path.join(site, "StopPlots", f), pdf)
        except Exception as e:  # report and carry on; one bad plot must not stop the run
            print(f"{f}: {type(e).__name__}: {e}", file=sys.stderr)
            continue
        json.dump(o, open(out, "w"), indent=1)
        print(f"{f}: {o['n_points']} points; legend {o['legend_points']}, plot says {o['plot_points']}",
              file=sys.stderr)


def attribution(o, p):
    """unique: one legend row has the colour and the legend lists every data
    set; shared-colour: several rows share it; partial-legend: the plot lists
    only some of its data sets, so an unlisted one may own the point."""
    if not o.get("legend_complete"):
        return "partial-legend"
    return "shared-colour" if len(p["candidates"]) > 1 else "unique"


# Words in a citation file's comment line that name the ion of a plot.
ION_WORDS = {1: r"\b(H|D|T|P|p|H\+|H2|H3|Protons?|Deuterons?|Tritons?|Hydrogen|Light [Ii]ons)\b",
             2: r"(\bHe\b|He3|3He|4He|[Aa]lpha|Helium|Light [Ii]ons)"}


def key_for(row, z1, ref_index):
    """The citation key of a legend row. The plot's own citation file is
    trusted. Old citation numbers are not unique across the compilation, so a
    number found only in another system's file is used just when that entry's
    comment names this plot's ion; otherwise the row keeps "ref<N>"."""
    if row["key"]:
        return row["key"]
    e = ref_index.get(row["ref"])
    if e and re.search(ION_WORDS[z1], " ".join(e["lines"])):
        return e["key"]
    return "ref%d" % row["ref"]


def export(work, path, ref_index):
    w = csv.writer(open(path, "w", newline=""), lineterminator="\n")
    w.writerow(COLUMNS)
    n = 0
    for f in sorted(glob.glob(os.path.join(work, "STOP*.json"))):
        o = json.load(open(f))
        z1, z2 = int(o["plot"][4:6]), int(o["plot"][6:8])
        refs, keys = {}, {}
        for row in o["legend"]:
            refs.setdefault(row["colour"], []).append(str(row["ref"]))
            keys.setdefault(row["colour"], []).append(key_for(row, z1, ref_index))
        for p in sorted(o["points"], key=lambda p: (p["colour"], p["energy_kev_per_amu"], p["px"][1])):
            if not p["candidates"]:
                continue  # a colour with no legend row: unattributable
            lab = p.get("label", "?")
            w.writerow([z1, "D" if lab == "D" else SYM[z1 - 1], z2, SYM[z2 - 1], p["energy_kev_per_amu"],
                        p["stopping_ev_cm2_per_1e15"], p["u_energy_rel"], p["u_stopping"],
                        "|".join(keys[p["colour"]]), "|".join(refs[p["colour"]]), attribution(o, p), lab,
                        o["plot"], p["colour"], p["px"][0], p["px"][1], "digitized-compilation", METHOD])
            n += 1
    print(f"{n} points -> {path}", file=sys.stderr)


if __name__ == "__main__":
    ions = [int(z) for z in sys.argv[4:]] or [1, 2]
    ref_index = citpdf.index(sys.argv[1])
    run(sys.argv[1], sys.argv[2], ions)
    export(sys.argv[2], sys.argv[3], ref_index)
