# Hand check: B in Si ranges (srim.org figure `SRIM/SRIMPICS/RANGES/B into Si - Ranges.gif`)

Ten points of `srim_org_ranges.csv` (B in Si), chosen at equal steps through the
figure's energy range, read by eye on 2026-10-30 at 8x magnification from the
operator's mirror of 2026-10-28 (attended session on the mirror host; issue #495).

Method: for each point a 61 x 61 px window around the recorded pixel was
magnified 8x with a crosshair at the recorded pixel (the template's reference
corner); the centre of the `*B` label of the recorded colour was read by eye
and converted with the figure's axis references (`figures.json`: x log, 10 keV
at px 102 and 1000 keV at px 737; y log, 100 nm at px 460 and 1000 nm at px
249). The references were checked against the tick labels: on the figure the
`10`, `100` and `1000` keV labels are centred at px 104, 419 and 734, and the
`100` and `1000` nm labels at px 463 and 250 (label centres, grid lines a pixel
or two away). E = energy (keV), R = target depth (nm).

 0  px=(176,502)  colour #8092AA  recorded E=17.95 R=59.55  by-eye centre=(183.5,507.9) -> E=18.06 (+0.6 %) R=59.31 (-0.4 %)
 1  px=(382,361)  colour #55B600  recorded E=79.84 R=278.9  by-eye centre=(389.5,366.0) -> E=80.45 (+0.8 %) R=278.9 (+0.0 %)  [crowded cluster; the olive-green *B of the recorded colour is at the expected place; centre read to about 2 px]
 2  px=(412,354)  colour #AA4955  recorded E=99.64 R=299.4  by-eye centre=(418.9,359.9) -> E=99.55 (-0.1 %) R=298.2 (-0.4 %)
 3  px=(437,334)  colour #000000  recorded E=119.5 R=374.5  by-eye centre=(443.5,338.8) -> E=119 (-0.4 %) R=375.5 (+0.3 %)
 4  px=(468,304)  colour #55FF55  recorded E=149.7 R=519.6  by-eye centre=(474.5,307.8) -> E=149 (-0.4 %) R=526.7 (+1.4 %)  [by-eye depth 1.3 % above the recorded value (about 1.2 px)]
 5  px=(508,291)  colour #5592FF  recorded E=199.4 R=595.5  by-eye centre=(514.4,296.0) -> E=199 (-0.2 %) R=598.8 (+0.5 %)
 6  px=(508,293)  colour #55FFFF  recorded E=199.4 R=582.6  by-eye centre=(514.4,297.0) -> E=199 (-0.2 %) R=592.3 (+1.7 %)  [by-eye depth 1.6 % above the recorded value (about 1.5 px); cyan *B partly over the light-blue one of point 5]
 7  px=(563,272)  colour #009200  recorded E=298.4 R=732.7  by-eye centre=(569.5,276.8) -> E=296.8 (-0.5 %) R=738.7 (+0.8 %)
 8  px=(603,253)  colour #009200  recorded E=399 R=901.5  by-eye centre=(609.8,258.6) -> E=397.4 (-0.4 %) R=900.3 (-0.1 %)
 9  px=(698,206)  colour #009200  recorded E=792.6 R=1514  by-eye centre=(705.5,211.6) -> E=795.8 (+0.4 %) R=1504 (-0.7 %)  [asterisk hidden under the SRIM curve; the B is clear]

Result: every by-eye reading lands on the recorded point within 1.6 % in depth
and 0.7 % in energy (worst 1.7 %, about 1.5 px); the sign of the
differences is not systematic. The digitizer's stated uncertainty is half a
pixel plus the anchor scatter (0.7 px). Nothing here reads the SRIM curve.
