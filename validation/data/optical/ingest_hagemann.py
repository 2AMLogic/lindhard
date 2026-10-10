#!/usr/bin/env python3
"""Convert a refractiveindex.info `tabulated nk` file of Hagemann, Gudat and
Kunz, J. Opt. Soc. Am. 65, 742 (1975) into an `OpticalElf` TOML (Al, Cu, C,
Au).

Usage: ingest_hagemann.py INPUT.yml MATERIAL OUTPUT.toml
       ingest_hagemann.py --patch BASE.yml PATCH.yml KEY OUTPUT.toml

Only the transcribed (wavelength, n, k) rows are read; nothing is smoothed or
extrapolated. The first form converts one table and joins nothing. The
`--patch` form (issue #314, KEY in PATCHED) writes every PATCH row and the
BASE rows outside the PATCH energy span, with no smoothing at the two seams:
the Au Table 6 of the report covers 1.5 to 350 eV only. Per row:

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
    "C": (
        "Measured optical constants n, k of glassy carbon, far infrared to x-ray "
        "(about 10 meV to 30 keV): H.-J. Hagemann, W. Gudat and C. Kunz, J. Opt. Soc. "
        "Am. 65, 742 (1975), doi:10.1364/JOSA.65.000742; full tables in DESY report "
        "SR-74/7 (1974), Table 8 (PDF pp. 58-59 of the open scan linked from the "
        "transcription), seen 2026-10-08; the JOSA paper itself (closed access) was not "
        "seen. The authors' joins (report Section 4, 'Glassy carbon', p. 20): "
        "reflectance by Kramers-Kronig analysis to 0.5 eV, k from angle-resolved "
        "reflectance on glassy carbon to 80 eV, literature transmission to 700 eV "
        "(scaled by 0.75 as a density correction to give 2 effective K electrons), "
        "Hubbell's compiled cross sections above 700 eV; n then by Kramers-Kronig. "
        "The values are given on the basis of a density of 1.5 g/cm3 (report p. 20). "
        "Rows read from the transcription C/nk/Hagemann.yml of the CC0 "
        "refractiveindex.info database (https://github.com/polyanskiy/refractiveindex.info-database, "
        "commit c5c2f18), fetched 2026-10-08; n, k at 17 energies checked against the "
        "scan (secondread_hagemann1975.json). Converted here: E = hc/lambda, ELF = "
        "Im[-1/eps] = eps2/(eps1^2+eps2^2) with eps1 = n^2-k^2, eps2 = 2nk. No "
        "smoothing, joining or extrapolation by lindhard. Operator ruling #34: "
        "committed as cited facts. Issue #148."
    ),
    "Au": (
        "Measured optical constants n, k of Au, far infrared to x-ray (about 5 meV "
        "to 150 keV): H.-J. Hagemann, W. Gudat and C. Kunz, J. Opt. Soc. Am. 65, "
        "742 (1975), doi:10.1364/JOSA.65.000742; full tables in DESY report SR-74/7 "
        "(1974), Table 5 (PDF pp. 49-51 of the open scan linked from the transcription), "
        "seen 2026-10-08; the JOSA paper itself (closed access) was not seen. Table 5 "
        "is the first of the report's two Au versions, the one fitted to the authors' "
        "transmission results, which gives 79 effective electrons (report p. 18); "
        "Table 6, extrapolated from reflectance data, is not used. The authors' joins "
        "(report Section 4, 'Gold', p. 17): absorption from Drude parameters to 0.6 eV, "
        "k from reflectance and transmission to 2.4 eV, k from reflectance by "
        "Kramers-Kronig analysis (UHV) to 10.5 eV, a fitted segment shaped after "
        "reflectance measurements to 20 eV, their own transmission to 117 eV (with the "
        "additive correction of report p. 13), literature transmission to 300 eV, "
        "above that a shape interpolated to 500 eV and Hubbell's compiled cross "
        "sections, fitted for a smooth connection at 300 and 500 eV; n then by "
        "Kramers-Kronig. The report's 1.00E-03 eV row (no n printed) is not "
        "transcribed. Rows read from the transcription Au/nk/Hagemann.yml of the CC0 "
        "refractiveindex.info database (https://github.com/polyanskiy/refractiveindex.info-database, "
        "commit c5c2f18), fetched 2026-10-08; n, k at 19 energies checked against the "
        "scan (secondread_hagemann1975.json). Converted here: E = hc/lambda, ELF = "
        "Im[-1/eps] = eps2/(eps1^2+eps2^2) with eps1 = n^2-k^2, eps2 = 2nk. No "
        "smoothing, joining or extrapolation by lindhard. Operator ruling #34: "
        "committed as cited facts. Issue #148."
    ),
}

# `--patch` outputs: KEY -> (material, provenance).
PATCHED = {
    "Au-T6": (
        "Au",
        "Measured optical constants n, k of Au, far infrared to x-ray (about 5 meV "
        "to 150 keV), composite of the two Au versions of H.-J. Hagemann, W. Gudat and "
        "C. Kunz, J. Opt. Soc. Am. 65, 742 (1975), doi:10.1364/JOSA.65.000742; full "
        "tables in DESY report SR-74/7 (1974), open scan of 96 PDF pages (SHA-256 "
        "f2643a91539506dd6794dde3b4810b91c3c7a24a57f98cf8b8b810e8e0931ed9), seen "
        "2026-10-09; the JOSA paper itself (closed access) was not seen. From 1.5 to "
        "350 eV: Table 6, 'Optical constants of Au from reflectivity-data fitted to "
        "Canfield et al. (10)' (PDF pp. 52-54, 124 rows, 1.50 to 350 eV, the whole of "
        "that table). Below 1.5 eV and above 350 eV: Table 5 (PDF pp. 49-51), the same "
        "25 rows as au_elf_hagemann1975.toml. The report (p. 18) gives the two tables "
        "as alternative versions: Table 5 from fitting the adjacent values to the "
        "authors' transmission results, which gives 79 effective electrons; Table 6 "
        "from extrapolating the measured reflectance data, which 'reproduces the "
        "energy-loss-spectra more closely'. The report prints Table 6 for 1.5 to 350 eV "
        "only and does not say how to continue it; the join is lindhard's (issue #314), "
        "not the authors': each Table 6 row replaces the Table 5 row at the same energy "
        "(the two tables share all 124 energies), and nothing is smoothed, so the ELF "
        "steps at both seams: Table 6 / Table 5 = 0.470 at 1.5 eV and 0.229 at 350 eV, "
        "and the next Table 5 knot (500 eV) is 2.79 times the Table 6 value at 350 eV. "
        "Rows read from the transcriptions Au/nk/Hagemann-2.yml (Table 6) and "
        "Au/nk/Hagemann.yml (Table 5) of the CC0 refractiveindex.info database "
        "(https://github.com/polyanskiy/refractiveindex.info-database, commit "
        "c5c2f18), fetched 2026-10-09; Table 6 n, k at 27 energies from 1.5 to 350 eV "
        "checked against the scan (secondread_hagemann1975.json), all equal. Converted "
        "here: E = hc/lambda, ELF = Im[-1/eps] = eps2/(eps1^2+eps2^2) with eps1 = "
        "n^2-k^2, eps2 = 2nk. No smoothing or extrapolation by lindhard. Operator ruling "
        "#34: committed as cited facts. Issue #314.",
    ),
}



def read_elf(inp):
    """{E [eV]: ELF} from the transcribed (wavelength, n, k) rows of INPUT.yml."""
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
    assert len(pts) == len(rows), "every row read must be written"
    return pts


def write(out, material, provenance, pts):
    es = sorted(pts)
    with open(out, "w", encoding="utf-8") as f:
        f.write(f'material = "{material}"\nprovenance = "{provenance}"\n')
        f.write("energy_ev = [\n" + "".join(f"  {e:.6e},\n" for e in es) + "]\n")
        f.write("elf = [\n" + "".join(f"  {pts[e]:.6e},\n" for e in es) + "]\n")
    return es


def main(inp, material, out):
    pts = read_elf(inp)
    es = write(out, material, PROVENANCE[material], pts)
    print(f"{material}: {len(pts)} rows read, {len(es)} written, {es[0]:.4g} to {es[-1]:.4g} eV")


def main_patch(base, patch, key, out):
    """PATCH rows over their own energy span [first, last], BASE rows outside it.

    Every PATCH row is written; a BASE row is written only if its energy lies
    outside the PATCH span (a BASE row at a PATCH end point is replaced by the
    PATCH row there). Nothing is smoothed at the two seams.
    """
    material, provenance = PATCHED[key]
    b, p = read_elf(base), read_elf(patch)
    lo, hi = min(p), max(p)
    kept = {e: y for e, y in b.items() if e < lo * (1 - 1e-6) or e > hi * (1 + 1e-6)}
    pts = {**kept, **p}
    assert len(pts) == len(kept) + len(p), "a kept BASE row coincides with a PATCH row"
    es = write(out, material, provenance, pts)
    print(
        f"{key}: {len(p)} rows from {patch} ({lo:.4g} to {hi:.4g} eV), {len(kept)} of "
        f"{len(b)} from {base}, {len(es)} written, {es[0]:.4g} to {es[-1]:.4g} eV"
    )


if __name__ == "__main__":
    if sys.argv[1] == "--patch":
        main_patch(*sys.argv[2:6])
    else:
        main(*sys.argv[1:4])
