# Bethe-Bloch and the density effect

Code: `lindhard/src/ion/stopping/bethe.rs` (`BetheBloch`, `EffectiveCharge`,
`density_effect_single_oscillator`).

## Model

At high velocity the electronic stopping cross section per target atom is
(Bethe 1930, 1932; Fano 1963; the same equation as the Particle Data Group
review of the passage of particles through matter):

\\[ S = \frac{4\pi (e^2)^2 z^2 Z_2}{m c^2 \beta^2}
\left[ \frac{1}{2}\ln\frac{2 m c^2 \beta^2\gamma^2 W_\mathrm{max}}{I^2}
- \beta^2 - \frac{\delta}{2} - \frac{C}{Z_2} + L_1 + L_2 \right], \\]

\\[ W_\mathrm{max} = \frac{2 m c^2 \beta^2\gamma^2}{1 + 2\gamma m/M + (m/M)^2}, \\]

with \\( m \\) the electron mass, \\( M \\) the ion mass, \\( z \\) the
projectile effective charge, \\( I \\) the mean excitation energy,
\\( \delta \\) the density-effect correction and \\( C \\) the shell
correction. \\( \beta^2 \\) is computed as
\\( \tau(\tau + 2)/(1 + \tau)^2 \\) with \\( \tau = E/Mc^2 \\), which stays
accurate at low energy.

- **Bloch term** (Bloch 1933), on by default:
  \\( L_2 = -y^2 \sum_{n \ge 1} \frac{1}{n(n^2 + y^2)} \\),
  \\( y = z\alpha/\beta \\).
- **Mean excitation energy**: by default the Bloch rule
  \\( I = 10\ \mathrm{eV} \cdot Z_2 \\) (Bloch 1933), a rough estimate
  (about 20 % in \\( I \\), about 2 % in \\( S \\)). The library takes a
  measured, cited value through `BetheBloch::with_mean_excitation_ev`.
- **Shell correction** \\( C/Z_2 \\) and **density effect** \\( \delta \\):
  0 by default; the library accepts constants from a cited source
  (`shell_over_z`, `density_delta`).
- **Barkas term** \\( L_1 \\): not implemented (see
  [Validity ranges and declined terms](stopping-validity.md)).

### Effective charge

| Model | Rust | \\( z \\) |
|---|---|---|
| Bare (default) | `EffectiveCharge::Bare` | \\( z = Z_1 \\), fully stripped; right for protons and alphas at high energy |
| Barkas empirical | `EffectiveCharge::BarkasEmpirical` | \\( z = Z_1 \left[1 - \exp\left(-125\beta / Z_1^{2/3}\right)\right] \\) (Barkas 1963, quoted by Northcliffe 1963) |

The CLI uses the bare charge; the Barkas form is available in the library and
is not applied by default.

### Single-oscillator density effect (opt-in)

The dielectric formulation of the density effect (Fermi 1940; Sternheimer
1952), in the form given by Fano (1963), is

\\[ \delta = \sum_i f_i \ln\\!\left(1 + \frac{L^2}{\omega_i^2}\right) - \frac{L^2 (1 - \beta^2)}{\beta^2\omega_p^2},
\qquad \sum_i \frac{f_i\\, \omega_p^2}{\omega_i^2 + L^2} = \frac{1}{\beta^2} - 1. \\]

Specialised to one oscillator (\\( f = 1 \\), \\( \omega_0 = I/\hbar \\)) the
constraint solves in closed form, \\( L^2 = \omega_p^2(\beta\gamma)^2 - \omega_0^2 \\),
and

\\[ \delta = \ln\\!\left[(\beta\gamma)^2 \left(\frac{\hbar\omega_p}{I}\right)^2\right] - 1
+ \frac{(I/\hbar\omega_p)^2}{(\beta\gamma)^2} \quad \text{for } \beta\gamma > \frac{I}{\hbar\omega_p}, \\]

and 0 below, with \\( \hbar\omega_p = \hbar\sqrt{n e^2/\varepsilon_0 m} \\) the
free-electron plasma energy. The one-oscillator reduction is ours, not quoted
from a paper. It is enabled with `BetheBloch::with_density_effect` (library
only).

## Selecting it

`stopping = "bethe-bloch"` in `[physics]`, `StoppingChoice::BetheBloch`:
bare charge, Bloch term on, Bloch-rule \\( I \\), no shell or density
correction, all loss nonlocal.

**Use it only for problems that stay above its range.** Where the bracket is
not positive the model returns `NotApplicable`, and the transport stops with
that error rather than silently switching model. For any ion that slows down
in the target (a semi-infinite substrate always) the run fails once the
energy falls below the bracket's zero: for example, 300 keV H into Si with
`stopping = "bethe-bloch"` stops with "bethe-bloch is not applicable at
85306 eV". No low-to-high energy interpolation is implemented, for the
reasons given in [Validity ranges and declined terms](stopping-validity.md).

## Assumptions

- First Born approximation plus the Bloch correction; the projectile charge
  is fixed by the effective-charge model.
- Without shell and Barkas corrections, expect errors of a few percent at 1
  to 10 MeV/u.

## Validity

\\( v \ge 3 v_0 Z_1^{2/3} \\) up to 1 GeV/u (the upper end because the
density effect is off by default). The single-oscillator density effect
switches on only at \\( \beta\gamma > I/\hbar\omega_p \\), about 5 to 6 for Si
(a proton of about 5 GeV), so it is zero over the whole advisory range. Real
materials switch the effect on much earlier (around \\( \beta\gamma \approx 1.5 \\)
for Si), because their oscillator spectrum is spread: the one-oscillator
form underestimates \\( \delta \\) through the transition and is only the
correct asymptote for ultra-relativistic ions.

## Verification status

Hand-calculated proton-in-Si values (1, 10 and 100 MeV, \\( I = 140 \\) eV,
CODATA constants) are unit tests; they are our own arithmetic, not taken from
any stopping table. The density-effect form is tested through its limits
(\\( \delta \to 0 \\) continuously at threshold, and
\\( \delta \to 2\ln(\beta\gamma\\, \hbar\omega_p / I) - 1 \\) at large
\\( \beta\gamma \\)). The Bloch rule, the Bloch series and the Barkas
effective charge have not been checked against the original papers
([`docs/data-provenance.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/data-provenance.md)).

## References

- H. Bethe, Ann. Phys. 397, 325 (1930); H. Bethe, Z. Phys. 76, 293 (1932).
- F. Bloch, Ann. Phys. 408, 285 (1933).
- U. Fano, Ann. Rev. Nucl. Sci. 13, 1 (1963).
- Particle Data Group, "Passage of particles through matter", in the Review
  of Particle Physics.
- W. H. Barkas, *Nuclear Research Emulsions* I (Academic Press, 1963).
- L. C. Northcliffe, Ann. Rev. Nucl. Sci. 13, 67 (1963).
- E. Fermi, Phys. Rev. 57, 485 (1940).
- R. M. Sternheimer, Phys. Rev. 88, 851 (1952).
