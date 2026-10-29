# History and literature of ion stopping

This chapter is a reading guide to how the theory of ion stopping and range
was built, from the first alpha-particle experiments to the codes of the
1980s. It names the papers that matter, says in a sentence what each
contributed, and links every model that lindhard implements to the manual
page where it lives. It is written for a reader who wants to know where a
model came from before trusting it. It is not a substitute for the papers.

## Where this chapter comes from

The shape of the story, and many of the references, follow the historical
review that James F. Ziegler wrote for *The Stopping and Range of Ions in
Solids* [Z85] and later put on the SRIM web site, srim.org [S03]. That
review, and the bibliography of measured stopping data that Ziegler compiled
over decades, are the reason a field this large has a shared literature at
all. The text below is our own summary, with every statement tied to a
primary paper. The srim.org site is no longer maintained, so each page
we cite is linked twice, to the live site and to a copy in the Internet
Archive:

| Page | srim.org | Internet Archive |
|---|---|---|
| Historical review | [HISTORY.htm](http://www.srim.org/SRIM/History/HISTORY.htm) | [2024-07-12 copy](https://web.archive.org/web/20240712200111/http://srim.org/SRIM/History/HISTORY.htm) |
| Stopping in compounds | [Compounds.htm](http://www.srim.org/SRIM/Compounds.htm) | [2025-05-24 copy](https://web.archive.org/web/20250524084924/http://www.srim.org/SRIM/Compounds.htm) |
| Ranges | [Ranges.htm](http://www.srim.org/SRIM/Ranges.htm) | [2024-12-12 copy](https://web.archive.org/web/20241212180843/http://www.srim.org/SRIM/Ranges.htm) |
| High-energy (H, He) stopping | [HESTOPPING.htm](http://www.srim.org/SRIM/SRIMPICS/HESTOPPING.htm) | [2024-10-06 copy](https://web.archive.org/web/20241006205259/http://www.srim.org/SRIM/SRIMPICS/HESTOPPING.htm) |

Prefer the archived copies: pages on the live site have since changed, and
some content on them is unrelated to the subject. lindhard builds on this
published literature as an independent implementation. It is not derived
from SRIM, uses none of its tables, and is not endorsed by its authors (see
[Conventions and the clean-room position](models/overview.md)). This chapter
contains no stopping data, curves or figures from srim.org.

The Monte Carlo side of the lineage, from the 1947 neutron histories to the
binary-collision codes, is told separately in
[`docs/history.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/history.md).
It cites some of the same papers; this chapter is about the stopping physics.

## Vocabulary

A fast ion in matter loses energy two ways. **Electronic stopping** is the
energy given to the target's electrons. **Nuclear stopping** is the energy
given to target atoms as a whole in elastic collisions. The rate of loss per
unit path, \\( dE/dx \\), is the *stopping power* (older papers) or
*stopping force* (newer ones); divided by the atomic density it is the
*stopping cross section*. Stopping depends on the ion's velocity much more
than on its energy, which is why the history below is organised by velocity.

## From radioactivity to the nuclear atom (1900–1915)

Soon after radioactivity was discovered, Marie Curie proposed that alpha
rays are material projectiles that lose speed in matter [C00]. Bragg and
Kleeman showed in 1905 that the range of alpha particles in a compound is
roughly the sum of contributions from its atoms [BK05]. That is the origin
of what lindhard calls [Bragg additivity](models/stopping-bragg.md).

Geiger and Marsden found that a small fraction of alpha particles bounce
back from a thin foil [GM09]. Rutherford showed that this needs a single
close collision with a compact, heavy, charged core [R11], which is the
nuclear atom and the scattering problem that
[the scattering integral](models/scattering.md) still solves.

Thomson [T12] and Bohr [B13, B15] then wrote the first theories of the
energy loss itself. Bohr treated target electrons as bound oscillators with
the frequencies of optical spectra, separated nuclear from electronic
losses, and saw that the unknown charge state of the moving ion limited his
accounting. The problems he listed (energy transfer to a quantised medium,
polarisation of the medium, the ion's charge state, screened nuclear
scattering, relativity) organise most of what follows.

## Quantum mechanics and fast particles (1930–1934)

Bethe [Be30, Be32] and Bloch [Bl33a, Bl33b] redid the problem with quantum
mechanics in the Born approximation. The result, the Bethe-Bloch formula,
is the standard description of a fast, fully stripped ion, and the base of
all high-energy stopping tables. Bethe and Heitler [BH34] extended the
treatment to radiative losses of fast electrons. Fermi's 1940 paper [F40]
added the polarisation of the medium, the origin of the *density effect*.
lindhard implements the formula, the effective charge options and the
density effect in
[Bethe-Bloch and the density effect](models/stopping-bethe-bloch.md).
For a survey of the relativistic regime see Ahlen [A80].

## Fission fragments and the effective charge (1939–1941)

The discovery of fission [HS39] produced heavy, partially stripped ions
with several MeV of energy, and a new question: how to apply a stopping
theory built for bare charges to an ion carrying its own electrons. Bohr
[B40, B41] proposed that an ion keeps the electrons that orbit more slowly
than it moves, and used the Thomas-Fermi atom [Th27] to turn this into an
effective charge that scales with velocity. Lamb [L40] argued from binding
energy instead, and Knipp and Teller [KT41] used the idea to scale proton
stopping to other ions. lindhard's effective-charge options are on the
[Bethe-Bloch](models/stopping-bethe-bloch.md) page.

## Screened nuclear collisions

Nuclear stopping needs the potential between two atoms. Thomas-Fermi
statistics [Th27] gave a universal screening function, and analytic
approximations to it followed: Lenz and Jensen in 1932 [Le32, J32] and
Molière in 1947 [M47]. Bohr [B48] argued for a screening length that
combines the two nuclear charges, and Firsov [Fi58] obtained a different
combination by fitting numerically computed potentials of colliding
Thomas-Fermi atoms. Later potentials replaced the statistical atom by
Hartree-Fock calculations: the Kr-C potential of Wilson, Haggmark and
Biersack [WHB77], and the universal potential of Ziegler, Biersack and
Littmark, fitted to many such calculations [BZ82, Z85]. lindhard offers
each of these in [Interatomic potentials](models/potentials.md), with the
screening lengths in [Screening lengths](models/screening-lengths.md).

## The slow-ion regime and the free electron gas (1947–1960)

For ions slower than the target's electrons the loss is proportional to the
ion's velocity. Fermi and Teller derived this for a degenerate electron gas
[FT47]. Lindhard's 1954 dielectric treatment [Li54], with the independent
work of Neufeld and Ritchie [NR55], put the electron gas on a full
quantum footing; Lindhard and Scharff used it in 1961 to give the
velocity-proportional stopping of keV ions [LS61]. That model is
[Lindhard-Scharff](models/stopping-lindhard-scharff.md). Northcliffe [N60]
collected measurements for many ions and targets into a ratio to proton
stopping at the same velocity, from which an empirical effective charge
follows. His review [N63] is a good entry point to the experimental side.

## LSS theory (1961–1963)

Lindhard, Scharff and Schiøtt [LSS63] joined the nuclear and electronic
parts into one theory with a universal reduced energy and reduced range, so
that range distributions could be computed for any ion-target pair from a
single curve. The reduced units, and the velocity-proportional electronic
loss, are the starting point of lindhard's
[Lindhard-Scharff](models/stopping-lindhard-scharff.md) model and its
[screening-length](models/screening-lengths.md) conventions. LSS theory
rests on the statistical atom, which is why it shows no shell structure.
Later work kept its framework and replaced its ingredients [Z85].

## Realistic atoms and computers (1970–1985)

Computers made it possible to replace statistical atoms by Hartree-Fock
ones. Rousseau, Chu and Powers did this for electronic stopping [RCP71],
and Wilson, Haggmark and Biersack for nuclear stopping [WHB77]. Brandt and
Kitagawa [BK82] revised Bohr's criterion for how many electrons an ion
keeps: they compared the ion's electron velocities with the Fermi velocity
of the solid, not with the ion velocity alone. In the same period Oen and
Robinson [OR76] gave a local, impact-parameter-dependent electronic loss,
which lindhard offers, mixed with the Lindhard-Scharff loss, in
[Oen-Robinson and the equipartition mix](models/stopping-oen-robinson.md).

On the transport side, Biersack and Haggmark's TRIM [BH80] made
binary-collision Monte Carlo practical for amorphous targets and introduced
the analytic scattering-angle formula that
[the scattering integral](models/scattering.md) keeps only as a
cross-check. [BCA transport and free-flight conventions](models/bca.md)
describes the transport loop. Ziegler, Biersack and Littmark's 1985 book
[Z85] and the SRIM program [ZZB10] that grew from it are the best-known
products of this period. For the earlier computer studies of ions in
crystals, and the theory of channeling, see
[`docs/history.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/history.md).

## Straggling and compounds

The fluctuation of energy loss about its mean, *straggling*, was already
treated by Bohr [B48] for free electrons; Fano's review [Fa63] covers the
corrections. lindhard's version is on the
[Energy-loss straggling](models/straggling.md) page. Stopping in compounds
and mixtures starts from Bragg's additivity rule [BK05] and corrects it for
chemical binding; the srim.org page on compounds (linked above) discusses the
corrections, and lindhard's treatment is
[Bragg additivity](models/stopping-bragg.md). Ziegler [Z85] and srim.org
describe a body of range and stopping data for many ion-target pairs. Data
and their provenance are tracked separately, in
[`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md).

## Where to read next

- Bohr's 1948 review [B48] is long but readable and covers most of what was
  known after 50 years.
- Northcliffe [N63] and Fano [Fa63] are the standard reviews of heavy-ion and
  light-ion stopping around the time of LSS.
- Whaling [W58] reviews the experiments of the 1950s. Ahlen [A80] covers
  relativistic ions.
- Sigmund's monograph *Particle Penetration and Radiation Effects* [Si06] is a
  modern textbook treatment.
- The International Commission on Radiation Units and Measurements publishes
  evaluated tables and their rationale: Report 49 for protons and alpha
  particles (1993) [ICRU49], and Report 37 for electrons and positrons
  (1984) [ICRU37]. lindhard does not use their tables; the reports are listed
  here for the physics in their text.
- A tribute to the man who compiled much of this literature:
  [James F. Ziegler (1935–2024)](history-ziegler.md).

## Corrections to the srim.org reference list

The list at the end of the srim.org review has inconsistencies. We checked
each reference below against the journal or publisher record (Crossref
unless noted) and used the corrected form in this chapter.

| srim.org gives | Record says | Note |
|---|---|---|
| Lenz, *Z. Phys.* 77, 713 (1932) and, in another entry, *Z. Phys.* 77, 722 (1932) | Lenz: 77, 713–721 (doi:10.1007/BF01342150). Jensen: 77, 722–745 (doi:10.1007/BF01342151) | 722 is the first page of Jensen's paper, which follows Lenz's in the same issue. |
| Rutherford, *Phil. Mag.* 21, 212 and 21, 699 (1911) | 21, 669–688 (doi:10.1080/14786440508637080) | One paper, not two. |
| Bloch, *Ann. Physik* 16, 287 (1933) | 16 (= vol. 408 in the continuous numbering), 285–320 (doi:10.1002/andp.19334080303) | First page is 285. |
| Rousseau, Chu and Powers, *Phys. Rev. A* 4, 1066 (**1970**) | 4, 1066–1070, published **1971** (doi:10.1103/PhysRevA.4.1066) | Year. |
| Lindhard, Nielsen, Scharff and Thomsen, *Mat. Fys. Medd.* 33, no. 10 (**1968**) | The same entry appears with the year 1963 elsewhere in the list; volume 33 is the 1963 volume that also holds LSS (no. 14) | Year inconsistent within the list; 1963 is consistent with the volume. No online journal record was found, so this is checked by volume only. |
| Robins, *New Principles of Gunnery*, dated 1804 and 1807 | First published 1742; the Internet Archive holds an 1805 edition with an account of his life (item `bub_gb_zokXxSnkinIC`) | The dates on srim.org belong to later editions. |
| Thomson, *Phil. Mag.* 6-23, 449 (1912) | Series 6, vol. 23, 449–457, "Ionization by moving electrified particles" (doi:10.1080/14786440408637241) | Same paper, series notation expanded. |

## References

Entries with a DOI were checked against Crossref. Entries marked *(no DOI)*
are not in Crossref; their citation is given as in the sources named, and
where we did not see the original we say so.

- [A80] S. P. Ahlen, "Theoretical and experimental aspects of the energy loss of relativistic heavily ionizing particles," *Rev. Mod. Phys.* 52, 121 (1980), doi:10.1103/RevModPhys.52.121.
- [B13] N. Bohr, "On the theory of the decrease of velocity of moving electrified particles on passing through matter," *Phil. Mag.* 25, 10 (1913), doi:10.1080/14786440108634305.
- [B15] N. Bohr, "On the decrease of velocity of swiftly moving electrified particles in passing through matter," *Phil. Mag.* 30, 581 (1915), doi:10.1080/14786441008635432.
- [B40] N. Bohr, "Scattering and stopping of fission fragments," *Phys. Rev.* 58, 654 (1940), doi:10.1103/PhysRev.58.654.
- [B41] N. Bohr, "Velocity-range relation for fission fragments," *Phys. Rev.* 59, 270 (1941), doi:10.1103/PhysRev.59.270.
- [B48] N. Bohr, "The penetration of atomic particles through matter," *Mat. Fys. Medd. Dan. Vid. Selsk.* 18, no. 8 (1948). *(no DOI; reprinted in* Niels Bohr Collected Works *vol. 8, doi:10.1016/S1876-0503(08)70172-5.)*
- [Be30] H. Bethe, "Zur Theorie des Durchgangs schneller Korpuskularstrahlen durch Materie," *Ann. Phys.* 397, 325 (1930), doi:10.1002/andp.19303970303.
- [Be32] H. Bethe, "Bremsformel für Elektronen relativistischer Geschwindigkeit," *Z. Phys.* 76, 293 (1932), doi:10.1007/BF01342532.
- [BH34] H. Bethe and W. Heitler, "On the stopping of fast particles and on the creation of positive electrons," *Proc. R. Soc. A* 146, 83 (1934), doi:10.1098/rspa.1934.0140.
- [BH80] J. P. Biersack and L. G. Haggmark, "A Monte Carlo computer program for the transport of energetic ions in amorphous targets," *Nucl. Instrum. Methods* 174, 257 (1980), doi:10.1016/0029-554X(80)90440-1.
- [Bl33a] F. Bloch, "Zur Bremsung rasch bewegter Teilchen beim Durchgang durch Materie," *Ann. Phys.* 408, 285 (1933), doi:10.1002/andp.19334080303.
- [Bl33b] F. Bloch, "Bremsvermögen von Atomen mit mehreren Elektronen," *Z. Phys.* 81, 363 (1933), doi:10.1007/BF01344553.
- [BK05] W. H. Bragg and R. Kleeman, "On the α particles of radium, and their loss of range in passing through various atoms and molecules," *Phil. Mag.* 10, 318 (1905), doi:10.1080/14786440509463378.
- [BK82] W. Brandt and M. Kitagawa, "Effective stopping-power charges of swift ions in condensed matter," *Phys. Rev. B* 25, 5631 (1982), doi:10.1103/PhysRevB.25.5631.
- [BZ82] J. P. Biersack and J. F. Ziegler, "Refined universal potentials in atomic collisions," *Nucl. Instrum. Methods* 194, 93 (1982), doi:10.1016/0029-554X(82)90496-7.
- [C00] M. Curie, *Comptes Rendus* 130, 76 (1900). *(no DOI; citation as given on srim.org, original not seen.)*
- [F40] E. Fermi, "The ionization loss of energy in gases and in condensed materials," *Phys. Rev.* 57, 485 (1940), doi:10.1103/PhysRev.57.485.
- [Fa63] U. Fano, "Penetration of protons, alpha particles, and mesons," *Ann. Rev. Nucl. Sci.* 13, 1 (1963), doi:10.1146/annurev.ns.13.120163.000245.
- [Fi58] O. B. Firsov, *Sov. Phys. JETP* 6, 534 (1958). *(no DOI; as cited on the* [screening-lengths](models/screening-lengths.md) *page.)*
- [FT47] E. Fermi and E. Teller, "The capture of negative mesotrons in matter," *Phys. Rev.* 72, 399 (1947), doi:10.1103/PhysRev.72.399.
- [GM09] H. Geiger and E. Marsden, "On a diffuse reflection of the α-particles," *Proc. R. Soc. A* 82, 495 (1909), doi:10.1098/rspa.1909.0054.
- [HS39] O. Hahn and F. Strassmann, "Über den Nachweis und das Verhalten der bei der Bestrahlung des Urans mittels Neutronen entstehenden Erdalkalimetalle," *Naturwissenschaften* 27, 11 (1939), doi:10.1007/BF01488241.
- [ICRU37] ICRU, *Stopping Powers for Electrons and Positrons*, ICRU Report 37 (1984). *(no DOI; title confirmed from the review by D. K. Brice, Nucl. Instrum. Methods B 12, 187 (1985), doi:10.1016/0168-583X(85)90718-9.)*
- [ICRU49] ICRU, *Stopping Powers and Ranges for Protons and Alpha Particles*, ICRU Report 49 (1993). *(no DOI; title confirmed from the review by J. O. Deasy, Med. Phys. 21, 709 (1994), doi:10.1118/1.597176.)*
- [J32] H. Jensen, "Die Ladungsverteilung in Ionen und die Gitterkonstante des Rubidiumbromids nach der statistischen Methode," *Z. Phys.* 77, 722 (1932), doi:10.1007/BF01342151.
- [KT41] J. Knipp and E. Teller, "On the energy loss of heavy ions," *Phys. Rev.* 59, 659 (1941), doi:10.1103/PhysRev.59.659.
- [L40] W. E. Lamb, "Passage of uranium fission fragments through matter," *Phys. Rev.* 58, 696 (1940), doi:10.1103/PhysRev.58.696.
- [Le32] W. Lenz, "Über die Anwendbarkeit der statistischen Methode auf Ionengitter," *Z. Phys.* 77, 713 (1932), doi:10.1007/BF01342150.
- [Li54] J. Lindhard, "On the properties of a gas of charged particles," *Mat. Fys. Medd. Dan. Vid. Selsk.* 28, no. 8 (1954). *(no DOI.)*
- [LS61] J. Lindhard and M. Scharff, "Energy dissipation by ions in the keV region," *Phys. Rev.* 124, 128 (1961), doi:10.1103/PhysRev.124.128.
- [LSS63] J. Lindhard, M. Scharff, H. E. Schiøtt, "Range concepts and heavy ion ranges (Notes on atomic collisions, II)," *Mat. Fys. Medd. Dan. Vid. Selsk.* 33, no. 14 (1963). *(no DOI.)*
- [M47] G. Molière, "Theorie der Streuung schneller geladener Teilchen I. Einzelstreuung am abgeschirmten Coulomb-Feld," *Z. Naturforsch. A* 2, 133 (1947), doi:10.1515/zna-1947-0302.
- [N60] L. C. Northcliffe, "Energy loss and effective charge of heavy ions in aluminum," *Phys. Rev.* 120, 1744 (1960), doi:10.1103/PhysRev.120.1744.
- [N63] L. C. Northcliffe, "Passage of heavy ions through matter," *Ann. Rev. Nucl. Sci.* 13, 67 (1963), doi:10.1146/annurev.ns.13.120163.000435.
- [NR55] J. Neufeld and R. H. Ritchie, "Passage of charged particles through plasma," *Phys. Rev.* 98, 1632 (1955), doi:10.1103/PhysRev.98.1632.
- [OR76] O. S. Oen and M. T. Robinson, "Computer studies of the reflection of light ions from solids," *Nucl. Instrum. Methods* 132, 647 (1976), doi:10.1016/0029-554X(76)90806-5.
- [R11] E. Rutherford, "The scattering of α and β particles by matter and the structure of the atom," *Phil. Mag.* 21, 669 (1911), doi:10.1080/14786440508637080.
- [RCP71] C. C. Rousseau, W. K. Chu, D. Powers, "Calculations of stopping cross sections for 0.8- to 2.0-MeV alpha particles," *Phys. Rev. A* 4, 1066 (1971), doi:10.1103/PhysRevA.4.1066.
- [S03] J. F. Ziegler, "Historical review," srim.org, 2003 (the date of the chapter): [live page](http://www.srim.org/SRIM/History/HISTORY.htm), [Internet Archive copy](https://web.archive.org/web/20240712200111/http://srim.org/SRIM/History/HISTORY.htm). The page says it is taken from [Z85].
- [Si06] P. Sigmund, *Particle Penetration and Radiation Effects*, Springer Series in Solid-State Sciences 151 (Springer, 2006), doi:10.1007/3-540-31718-X.
- [T12] J. J. Thomson, "Ionization by moving electrified particles," *Phil. Mag.* 23, 449 (1912), doi:10.1080/14786440408637241.
- [Th27] L. H. Thomas, "The calculation of atomic fields," *Proc. Camb. Phil. Soc.* 23, 542 (1927), doi:10.1017/S0305004100011683.
- [W58] W. Whaling, "The energy loss of charged particles in matter," in *Handbuch der Physik* 34, 193 (Springer, 1958), doi:10.1007/978-3-642-45898-9_5.
- [WHB77] W. D. Wilson, L. G. Haggmark, J. P. Biersack, "Calculations of nuclear stopping, ranges, and straggling in the low-energy region," *Phys. Rev. B* 15, 2458 (1977), doi:10.1103/PhysRevB.15.2458.
- [Z85] J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in Solids*, vol. 1 of the series *The Stopping and Ranges of Ions in Matter* (Pergamon, New York, 1985). *(no DOI; record checked in Open Library, which gives the title and publisher. The series title is written both with "Range" and "Ranges" in the sources we saw.)*
- [ZZB10] J. F. Ziegler, M. D. Ziegler, J. P. Biersack, "SRIM – The stopping and range of ions in matter (2010)," *Nucl. Instrum. Methods B* 268, 1818 (2010), doi:10.1016/j.nimb.2010.02.091.
