# A historical note: from the 1947 neutron histories to lindhard

lindhard's transport loop is older than its physics. Its structure, a
random free flight, a collision, then a stack of secondaries to follow, is
the one John von Neumann sketched for fission neutrons at Los Alamos in
1947. This page traces that lineage and shows where lindhard's choices sit
in it. References are at the end. The facts were checked against them on
2026-10-04.

## Solitaire, a letter, and ENIAC (1946–1949)

In 1946, while recovering from an illness, Stanislaw Ulam wondered what
the odds were of winning a game of Canfield solitaire. Working the
combinatorics out seemed hopeless; simply playing many games and counting
did not. He saw that neutron diffusion could be treated the same way [1].

On **11 March 1947** von Neumann wrote to Robert Richtmyer, head of the Los
Alamos Theoretical Division, setting out how. The letter followed
individual neutron histories through a fissioning assembly: draw a free
path, decide what kind of collision occurs, draw the new direction and
energy, and track any fission neutrons as new histories. He proposed
running it on ENIAC [1, 2].

The calculations ran on ENIAC in three campaigns: **April–May 1948**,
**22 October–7 November 1948**, and **May–June 1949** [3, 21]. Klára Dán von
Neumann was the principal coder, and she and Nicholas Metropolis set up and
ran the machine. The first-run program, completed in December 1947 and run
in spring 1948 after ENIAC was converted to a new control mode, was the
first program written in the modern stored-program style ever executed
[3].

While ENIAC was unavailable, Enrico Fermi devised an analog shortcut: the
**FERMIAC**, a brass-and-acrylic trolley that a person rolled across a
drawing of the assembly to trace neutron paths by hand. Percy King built it
in 1947, and it is now in the Bradbury Science Museum in Los Alamos [4].

The method was described publicly in 1949 in Metropolis and Ulam's "The
Monte Carlo Method" [5]. MANIAC I, Los Alamos's own computer, came online
in 1952 [1, 6]. It is remembered mainly for the Metropolis algorithm
(1953) [7] and for the Fermi–Pasta–Ulam study (1955), which Mary Tsingou
coded [8].

## The loop, then and now

| Step | 1947 neutron history | lindhard ion history |
|---|---|---|
| Free flight | Exponential path length from the total cross section | Amorphous targets: a free-path convention set by the atomic density ([#5](https://github.com/2AMLogic/lindhard/issues/5)). Crystals: a deterministic search for the next lattice atom ([#20](https://github.com/2AMLogic/lindhard/issues/20)) |
| Collision | Scatter, absorb or fission, chosen by tabulated cross sections | Impact parameter → scattering angle from the screened-potential scattering integral ([#3](https://github.com/2AMLogic/lindhard/issues/3)) |
| Between collisions | Nothing: a neutron is neutral | Continuous electronic energy loss ([#4](https://github.com/2AMLogic/lindhard/issues/4)): Bohr and Bethe at high energy, Lindhard at low energy |
| Branching | Fission neutrons pushed onto a stack | Recoils above the displacement energy pushed onto a stack |
| Fate of the tree | Can multiply: a chain reaction | Can only divide the incident energy, so every cascade dies out |

The middle row is where the physics changes. A neutron loses energy only
at collisions. A charged particle also loses energy continuously to the
target's electrons. That continuous term goes back to Bohr (1913) [22] and
Bethe (1930) [23], whose stopping formulas hold at high velocity. What
Lindhard and Scharff's 1961 paper [9] and the LSS range theory of 1963 [10]
added is the low-velocity regime, where the loss is proportional to the
particle's velocity, and a unified theory of range.

## Charged particles: three branches

**Electrons.** Electrons undergo so many small collisions that following
each one was out of reach for the computers of the time. In 1963 Martin
Berger introduced *condensed history*, which groups many collisions into a
single step [11]. lindhard's electron engine
([#11](https://github.com/2AMLogic/lindhard/issues/11)) goes back to
simulating events one by one. At the energies it targets (about 10 eV to
50 keV) this is affordable today, and condensed history is least
trustworthy there. Its inelastic model starts from Lindhard's 1954
dielectric function [12], extended to finite momentum transfer by Penn
(1987) [13].

**Ions in crystals.** In 1965 Lindhard showed how a crystal lattice
steers energetic ions along its open channels [14]. The simulation came
first: at Oak Ridge in 1963, Robinson and Oen ran binary-collision computer
studies of ions slowing down in crystals and found channeling [24], two
years before Lindhard's theory. Robinson and Torrens's MARLOWE code (1974)
[15] was the mature, general form of that line of work. lindhard's M2
milestone ([#12](https://github.com/2AMLogic/lindhard/issues/12))
addresses the same problem class, built only from the published papers
(Lindhard 1965 [14]; Robinson and Oen 1963 [24]; Robinson and Torrens 1974
[15]).

**Ions in amorphous targets.** Biersack and Haggmark's TRIM (1980) [16]
made the free-flight-and-collision loop fast and practical for amorphous
targets, using an analytic "magic formula" for the scattering angle. With
Ziegler's stopping work [17] it grew into SRIM, the code whose problem class lindhard re-implements from
the published physics.

## Random numbers, then and now

Von Neumann generated his random numbers arithmetically with the
*middle-square* method, and he was frank about its weaknesses [18]. In
1955 RAND published *A Million Random Digits with 100,000 Normal Deviates*
as a trusted source [19]. lindhard's answer to the same worry is a
counter-based stream for each particle, keyed on (seed, particle index)
([#2](https://github.com/2AMLogic/lindhard/issues/2)). A run is
reproducible bit for bit at any thread count. In 1948 reproducing a run
meant keeping the punched cards.

## Jens Lindhard (1922–1997)

Jens Lindhard was born on 26 February 1922 and died on 15 October 1997. He
was professor of theoretical physics at Aarhus University and president
of the Royal Danish Academy of Sciences and Letters from 1981 to 1988
[20]. Three of his results run through this project: the dielectric
function (1954) [12], the theory of energy loss and range (1961, 1963)
[9, 10], and the theory of channeling (1965) [14]. Hence the name.

## References

1. R. Eckhardt, "Stan Ulam, John von Neumann, and the Monte Carlo Method," *Los Alamos Science* 15 (1987) 131–137. Reproduces the 1947 letter and Ulam's account.
2. A. Sood et al., "Neutronics Calculation Advances at Los Alamos: Manhattan Project to Monte Carlo," arXiv:2103.06260.
3. T. Haigh, M. Priestley, C. Rope, "Los Alamos Bets on ENIAC: Nuclear Monte Carlo Simulations, 1947–1948," *IEEE Annals of the History of Computing* 36(3), 42–63 (2014), doi:10.1109/MAHC.2014.40.
4. F. Coccetti, "The Fermiac or Fermi's Trolley," *Il Nuovo Cimento C* 39, 296 (2016), doi:10.1393/ncc/i2016-16296-7.
5. N. Metropolis, S. Ulam, "The Monte Carlo Method," *J. Am. Stat. Assoc.* 44(247), 335–341 (1949), doi:10.1080/01621459.1949.10483310.
6. N. Metropolis, "The Beginning of the Monte Carlo Method," *Los Alamos Science* 15 (1987) 125–130.
7. N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, E. Teller, "Equation of State Calculations by Fast Computing Machines," *J. Chem. Phys.* 21, 1087–1092 (1953), doi:10.1063/1.1699114.
8. E. Fermi, J. Pasta, S. Ulam, "Studies of Nonlinear Problems I," Los Alamos report LA-1940 (1955), OSTI 4376203.
9. J. Lindhard, M. Scharff, "Energy Dissipation by Ions in the keV Region," *Phys. Rev.* 124, 128 (1961), doi:10.1103/PhysRev.124.128.
10. J. Lindhard, M. Scharff, H. E. Schiøtt, "Range Concepts and Heavy Ion Ranges (Notes on Atomic Collisions, II)," *Mat. Fys. Medd. Dan. Vid. Selsk.* 33, no. 14 (1963).
11. M. J. Berger, "Monte Carlo Calculation of the Penetration and Diffusion of Fast Charged Particles," *Methods in Computational Physics* 1, 135–215 (Academic Press, 1963).
12. J. Lindhard, "On the Properties of a Gas of Charged Particles," *Mat. Fys. Medd. Dan. Vid. Selsk.* 28, no. 8 (1954).
13. D. R. Penn, "Electron Mean-Free-Path Calculations Using a Model Dielectric Function," *Phys. Rev. B* 35, 482 (1987), doi:10.1103/PhysRevB.35.482.
14. J. Lindhard, "Influence of Crystal Lattice on Motion of Energetic Charged Particles," *Mat. Fys. Medd. Dan. Vid. Selsk.* 34, no. 14 (1965).
15. M. T. Robinson, I. M. Torrens, "Computer Simulation of Atomic-Displacement Cascades in Solids in the Binary-Collision Approximation," *Phys. Rev. B* 9, 5008 (1974), doi:10.1103/PhysRevB.9.5008.
16. J. P. Biersack, L. G. Haggmark, "A Monte Carlo Computer Program for the Transport of Energetic Ions in Amorphous Targets," *Nucl. Instrum. Methods* 174, 257–269 (1980), doi:10.1016/0029-554X(80)90440-1.
17. J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in Solids* (Pergamon, 1985).
18. J. von Neumann, "Various Techniques Used in Connection with Random Digits," in *Monte Carlo Method*, NBS Applied Mathematics Series 12, 36–38 (1951).
19. RAND Corporation, *A Million Random Digits with 100,000 Normal Deviates* (Free Press, 1955).
20. J. U. Andersen, P. Sigmund, "Jens Lindhard" (obituary), *Physics Today* 51(9), 89–90 (1998), doi:10.1063/1.882460.
21. T. Haigh, M. Priestley, C. Rope, *ENIAC in Action: Making and Remaking the Modern Computer* (MIT Press, 2016).
22. N. Bohr, "On the Theory of the Decrease of Velocity of Moving Electrified Particles on Passing through Matter," *Phil. Mag.* 25, 10–31 (1913), doi:10.1080/14786440108634305.
23. H. Bethe, "Zur Theorie des Durchgangs schneller Korpuskularstrahlen durch Materie," *Ann. Phys.* 397, 325–400 (1930), doi:10.1002/andp.19303970303.
24. M. T. Robinson, O. S. Oen, "Computer Studies of the Slowing Down of Energetic Atoms in Crystals," *Phys. Rev.* 132, 2385 (1963), doi:10.1103/PhysRev.132.2385.
