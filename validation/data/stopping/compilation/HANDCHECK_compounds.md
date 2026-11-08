# Hand check: compound stopping figures

Ten points of `srim_org_compounds.csv`, five on H in water (`SRIM/Compounds/H into Water.gif`)
and five on He in polystyrene (`SRIM/Compounds/He into Polystyrene.gif`), chosen at equal steps
through each figure's energy range, read by eye on 2026-10-30 at 8x
magnification from the operator's mirror of 2026-10-28 (attended session on
the mirror host; issue #498).

Method as in `../../ranges/compilation/HANDCHECK.md`: the centre of the
marker label (a letter on H in water, `He` on He in polystyrene) read by eye
and converted with the figure's axis references from `figures.json` (H in
water: x log, 10 keV/amu at px 104 and 1000 at px 870; y linear, 10 at px 452
and 30 at px 140 eV/(10^15 molecules/cm2). He in polystyrene: x log, 40 at px
242 and 400 at px 704; y linear, 10 at px 471 and 30 at px 198
eV/(10^15 atoms/cm2)). Checked on the tick labels: H in water `10`, `20`, `30`
centred at px 450, 297, 143; He in polystyrene `10` and `30` at px 473 and
200. E = energy (keV/amu), S = stopping in the figure's unit.

## H in water

 0  px=(102,438)  colour #0000FF  recorded E=10.21 S=10.51  by-eye centre=(106.6,443.9) -> E=10.16 (-0.5 %) S=10.52 (+0.1 %)
 1  px=(217,401)  colour #0000FF  recorded E=20.33 S=12.88  by-eye centre=(221.6,406.9) -> E=20.28 (-0.2 %) S=12.89 (+0.1 %)
 2  px=(333,280)  colour #0000FF  recorded E=40.73 S=20.64  by-eye centre=(337.4,285.9) -> E=40.68 (-0.1 %) S=20.65 (+0.0 %)  [B/D/E cluster on the curve; the E of the recorded point is the lower one]
 3  px=(401,228)  colour #0000FF  recorded E=61.2 S=23.97  by-eye centre=(405.0,231.8) -> E=61.08 (-0.2 %) S=24.12 (+0.6 %)  [the B glyph is taller than the others; centre read to about 1 px]
 4  px=(554,225)  colour #FF0000  recorded E=153 S=24.23  by-eye centre=(558.6,230.9) -> E=153.8 (+0.5 %) S=24.17 (-0.2 %)

## He in polystyrene

 0  px=(157,379)  colour #0000FF  recorded E=27.32 S=16.34  by-eye centre=(166.0,384.2) -> E=27.39 (+0.2 %) S=16.36 (+0.1 %)
 1  px=(305,278)  colour #0000FF  recorded E=57.12 S=23.74  by-eye centre=(314.0,283.6) -> E=57.27 (+0.2 %) S=23.73 (-0.0 %)
 2  px=(390,235)  colour #0000FF  recorded E=87.26 S=26.89  by-eye centre=(399.0,240.9) -> E=87.48 (+0.2 %) S=26.86 (-0.1 %)
 3  px=(556,240)  colour #FF00FF  recorded E=199.6 S=26.52  by-eye centre=(564.8,245.6) -> E=199.8 (+0.1 %) S=26.51 (-0.0 %)
 4  px=(695,316)  colour #AA4900  recorded E=399 S=20.95  by-eye centre=(703.8,321.2) -> E=399.5 (+0.1 %) S=20.97 (+0.1 %)

Result: every reading lands on the recorded point within 0.7 % (worst
0.6 %, under a pixel). No curve is read; the units and the
per-molecule / per-atom basis are those printed on each figure's y axis and
recorded in `value_unit`.
