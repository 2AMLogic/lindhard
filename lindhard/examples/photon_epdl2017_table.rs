//! Regenerate `lindhard/src/photon/epdl2017.toml` from a local copy of
//! `EPDL2017.ALL` (https://www.nndc.bnl.gov/endf/epics/ENDF2017/EPDL2017.ALL).
//!
//! ```text
//! cargo run -p lindhard --release --example photon_epdl2017_table -- \
//!     /path/to/EPDL2017.ALL > lindhard/src/photon/epdl2017.toml
//! ```
//!
//! The file is the credit header below plus `PhotonTable::to_toml_string`, so
//! no value is typed by hand. The SHA-256 of the source is printed to stderr.

use lindhard::photon::PhotonTable;
use sha2::{Digest, Sha256};

fn main() {
    let path = std::env::args().nth(1).expect("usage: <EPDL2017.ALL>");
    let bytes = std::fs::read(&path).expect("read the source file");
    let sha: String = Sha256::digest(&bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    eprintln!("sha256 {sha}");
    let provenance = format!(
        "EPDL2017 (EPICS2017), D. E. Cullen, IAEA-NDS-225 Rev. 1 (February 2018), ENDF-6 file EPDL2017.ALL (header \"Eval-Dec17\", \"Dist-Dec17\"; SHA-256 {sha}) retrieved 2026-11-10 from https://www.nndc.bnl.gov/endf/epics/ENDF2017/EPDL2017.ALL; Z=1..92, 1 keV to 1 MeV, thinned to 0.1 %"
    );
    let table = PhotonTable::from_epdl_endf6_str(
        &String::from_utf8(bytes).expect("ASCII source"),
        &provenance,
    )
    .expect("read the EPDL2017 file");
    println!(
        "# Photon interaction data, Z = 1..92, photon energy 1 keV to 1 MeV.
# Source: the Evaluated Photon Data Library (EPDL97) as distributed in EPICS2017,
# file EPDL2017.ALL (ENDF-6 format, header \"Eval-Dec17\", \"Dist-Dec17\"), retrieved
# 2026-11-10 from https://www.nndc.bnl.gov/endf/epics/ENDF2017/EPDL2017.ALL.
# SHA-256 of the source file: {sha}
# Sections read: MF=23 MT=522 (photoelectric), MT=504 (incoherent), MT=502
# (coherent), in barns against eV; MF=27 MT=502 (form factor) and MT=504
# (scattering function) against x = sin(theta/2)/lambda in 1/Angstrom.
# Credit: D. E. Cullen, author of the data library, \"A Survey of Photon Cross
# Section Data for use in EPICS2017\", IAEA-NDS-225 Rev. 1 (February 2018); issued
# by the IAEA Nuclear Data Section / National Nuclear Data Center (NNDC), BNL
# (see also IAEA-NDS-218).
# Notice: the source carries a warranty disclaimer only; no warranty is made for
# these values. Generated mechanically by lindhard's EPDL reader (window, then
# thinning to 0.1 % relative error); no value typed by hand. Committed under the
# operator ruling of 2026-11-10 recorded in docs/data-provenance.md.
"
    );
    print!("{}", table.to_toml_string().expect("serialize"));
}
