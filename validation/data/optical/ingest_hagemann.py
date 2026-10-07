#!/usr/bin/env python3
"""Convert a refractiveindex.info `tabulated nk` file of Hagemann, Gudat and
Kunz, J. Opt. Soc. Am. 65, 742 (1975) into an `OpticalElf` TOML.

Usage: ingest_hagemann.py INPUT.yml MATERIAL OUTPUT.toml

Only the transcribed (wavelength, n, k) rows are read; nothing is smoothed,
extrapolated or joined. Per row:

    E [eV]  = h c / lambda, with h c = 1.239841984 eV um (exact SI h, c, e)
    eps1    = n^2 - k^2
    eps2    = 2 n k
    ELF     = Im[-1/eps] = eps2 / (eps1^2 + eps2^2)

A repeated wavelength is an error (no row is ever dropped silently).
The input is untrusted data: it is parsed as plain text, never executed.
"""
import sys

H, C, E = 6.62607015e-34, 299792458.0, 1.602176634e-19
HC_EV_UM = H * C / E * 1e6

PROVENANCE = {
    "Al": (
        "Measured optical constants n, k of Al, far infrared to x-ray (about 1 meV "
        "to 120 keV): H.-J. Hagemann, W. Gudat and C. Kunz, J. Opt. Soc. Am. 65, "
        "742 (1975), doi:10.1364/JOSA.65.000742; full tables in DESY report SR-74/7 "
        "(1974), Table 2 (PDF pp. 40-42 of the open scan linked from the transcription), "
        "seen 2026-10-07; the JOSA paper itself (closed access) was not seen. The "
        "authors' joins (report Section 4, pp. 14-15): absorption from Drude parameters "
        "to 0.155 eV, interpolated to 0.8 eV, ellipsometry to 2.0 eV, Drude estimate to "
        "13.5 eV, k from angle-resolved reflectance and transmission to 36 eV, their own "
        "UHV transmission to 74 eV (fitted at 74 eV), literature transmission to 150 eV "
        "and to 300 eV (shape from two sets, scale from a third), Hubbell's compiled "
        "cross sections above 300 eV; n then by Kramers-Kronig. Rows read "
        "from the transcription Al/nk/Hagemann.yml of the CC0 refractiveindex.info "
        "database (https://github.com/polyanskiy/refractiveindex.info-database, "
        "commit c5c2f18), fetched 2026-10-07; n, k at 17 energies checked against the scan (secondread_hagemann1975.json). Converted here: E = hc/lambda, ELF = Im[-1/eps] = "
        "eps2/(eps1^2+eps2^2) with eps1 = n^2-k^2, eps2 = 2nk. No smoothing, "
        "joining or extrapolation by lindhard. Operator ruling #34: committed as "
        "cited facts. Issue #98."
    ),
    "Cu": (
        "Measured optical constants n, k of Cu, far infrared to x-ray (about 5 meV "
        "to 50 keV): H.-J. Hagemann, W. Gudat and C. Kunz, J. Opt. Soc. Am. 65, "
        "742 (1975), doi:10.1364/JOSA.65.000742; full tables in DESY report SR-74/7 "
        "(1974), Table 3 (PDF pp. 43-45 of the open scan linked from the transcription), "
        "seen 2026-10-07; the JOSA paper itself (closed access) was not seen. The "
        "authors' joins (report Section 4, pp. 15-16): absorption from Drude parameters "
        "to 0.5 eV, k from reflectance and transmission to 6.5 eV, an interpolated "
        "curve between two reflectance analyses, scaled to one of them, to 13 eV, their "
        "own transmission to 150 eV (with an additive correction, report p. 13), "
        "literature transmission to 450 eV, Hubbell's compiled cross sections above "
        "450 eV; n then by Kramers-Kronig. "
        "The report's 2.00E-03 eV row (no n printed) is not transcribed. Rows read "
        "from the transcription Cu/nk/Hagemann.yml of the CC0 refractiveindex.info "
        "database (https://github.com/polyanskiy/refractiveindex.info-database, "
        "commit c5c2f18), fetched 2026-10-07; n, k at 12 energies checked against the scan (secondread_hagemann1975.json). Converted here: E = hc/lambda, ELF = Im[-1/eps] = "
        "eps2/(eps1^2+eps2^2) with eps1 = n^2-k^2, eps2 = 2nk. No smoothing, "
        "joining or extrapolation by lindhard. Operator ruling #34: committed as "
        "cited facts. Issue #98."
    ),
}


def main(inp, material, out):
    rows, in_data = [], False
    for line in open(inp, encoding="utf-8"):
        s = line.strip()
        if s.startswith("data: |"):
            in_data = True
            continue
        if in_data:
            parts = s.split()
            if len(parts) != 3:
                break
            rows.append(tuple(float(x) for x in parts))
    pts = {}
    for lam, n, k in rows:
        e = HC_EV_UM / lam
        if e in pts:
            sys.exit(f"{inp}: duplicate wavelength {lam} um; refusing to drop a row")
        eps1, eps2 = n * n - k * k, 2 * n * k
        pts[e] = eps2 / (eps1 * eps1 + eps2 * eps2)
    es = sorted(pts)
    assert len(es) == len(rows), "every row read must be written"
    with open(out, "w", encoding="utf-8") as f:
        f.write(f'material = "{material}"\nprovenance = "{PROVENANCE[material]}"\n')
        f.write("energy_ev = [\n" + "".join(f"  {e:.6e},\n" for e in es) + "]\n")
        f.write("elf = [\n" + "".join(f"  {pts[e]:.6e},\n" for e in es) + "]\n")
    print(f"{material}: {len(rows)} rows read, {len(es)} written, {es[0]:.4g} to {es[-1]:.4g} eV")


if __name__ == "__main__":
    main(*sys.argv[1:4])
