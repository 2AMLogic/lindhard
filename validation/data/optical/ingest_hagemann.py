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

The input is untrusted data: it is parsed as plain text, never executed.
"""
import sys

H, C, E = 6.62607015e-34, 299792458.0, 1.602176634e-19
HC_EV_UM = H * C / E * 1e6

PROVENANCE = {
    "Al": (
        "Measured optical constants n, k of Al, far infrared to x-ray (about 1 meV "
        "to 120 keV): H.-J. Hagemann, W. Gudat and C. Kunz, J. Opt. Soc. Am. 65, "
        "742 (1975), doi:10.1364/JOSA.65.000742 (DESY report SR-74/7 (1974)); the "
        "journal table was not seen (closed access). Rows read from the "
        "transcription Al/nk/Hagemann.yml of the CC0 refractiveindex.info "
        "database (https://github.com/polyanskiy/refractiveindex.info-database, "
        "commit c5c2f18), fetched 2026-10-07. Converted here: E = hc/lambda, ELF = Im[-1/eps] = "
        "eps2/(eps1^2+eps2^2) with eps1 = n^2-k^2, eps2 = 2nk. No smoothing, "
        "joining or extrapolation by lindhard. Operator ruling #34: committed as "
        "cited facts. Issue #98."
    ),
    "Cu": (
        "Measured optical constants n, k of Cu, far infrared to x-ray (about 5 meV "
        "to 50 keV): H.-J. Hagemann, W. Gudat and C. Kunz, J. Opt. Soc. Am. 65, "
        "742 (1975), doi:10.1364/JOSA.65.000742 (DESY report SR-74/7 (1974)); the "
        "journal table was not seen (closed access). Rows read from the "
        "transcription Cu/nk/Hagemann.yml of the CC0 refractiveindex.info "
        "database (https://github.com/polyanskiy/refractiveindex.info-database, "
        "commit c5c2f18), fetched 2026-10-07. Converted here: E = hc/lambda, ELF = Im[-1/eps] = "
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
        eps1, eps2 = n * n - k * k, 2 * n * k
        pts[e] = eps2 / (eps1 * eps1 + eps2 * eps2)
    es = sorted(pts)
    with open(out, "w", encoding="utf-8") as f:
        f.write(f'material = "{material}"\nprovenance = "{PROVENANCE[material]}"\n')
        f.write("energy_ev = [\n" + "".join(f"  {e:.6e},\n" for e in es) + "]\n")
        f.write("elf = [\n" + "".join(f"  {pts[e]:.6e},\n" for e in es) + "]\n")
    print(f"{material}: {len(rows)} rows read, {len(es)} written, {es[0]:.4g} to {es[-1]:.4g} eV")


if __name__ == "__main__":
    main(*sys.argv[1:4])
