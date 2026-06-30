# NOTES

The companion to the README: what the package computes, how the paper used it, what it found, and the choices behind it, for coauthors and for the next person to open the repository.

## What the package did for the paper

The package simulates one experiment. A chain of spins starts polarised, the paramagnetic ground state of the transverse field; the spin at the centre is flipped, and the disturbance spreads outward as a light cone bounded by the Lieb-Robinson limit. At the front of the cone the magnetisation focuses into a caustic, the same focusing that prints bright lines on the floor of a sunlit pool, and the spacing of its fringes shrinks with distance by a fixed power. The paper used the package to measure that power across the phase diagram, and the package produced the calculations reported there.

The chain is evolved with TEBD on a matrix product state, which reaches several hundred sites, so the cone runs long enough to resolve the scaling. Each run records two observables, the raw magnetisation <Z_j(t)> and its no-flip subtraction dZ of Eq. 22, and these are the inputs the scaling analysis reads.

## What it found

In the paramagnetic phase, where the field dominates and Jxx < hz, the fringe spacing scales with the 2/3 exponent of an Airy fold, the universal value for the simplest caustic. The fit converges to about 0.68 against that 2/3, and the origin of the residual gap, whether finite N, finite bond dimension, or a true many-body effect, is not settled. Approaching the transition at Jxx = hz the exponent departs from 2/3. The caustic and the 2/3 scaling survive a weak integrability-breaking coupling Jzz, and the caustic is lost as Jzz grows. These four regimes, the two phases with weak and strong integrability breaking, are the committed problem instances in problems/.

## How it works, and the conventions behind it

The central flip is injected as a single down spin on the polarised background, matching Eq. 2, and a no-flip companion runs on the same grid and the same gates, so the two subtract step for step and the boundary cancels. The Hamiltonian is recorded in the paper's Pauli operators and signs, then converted to the ITensors spin-half operators the evolver uses, where S_a = sigma_a / 2 carries a factor 2 on the field and 4 on the bonds. That conversion and the global sign are where a slip would otherwise be invisible, so the test suite pins them against the exact engine on every change. The field is applied once per site as its own layer, the bonds as a brickwall of even and odd pairs. The step dt = 0.1 is converged, and the truncation controls, the bond-dimension cap and the cutoff, live in parameters.jl. The environment is pinned by a committed Manifest, so the pipeline reproduces from a clean checkout.

One Hamiltonian carries the whole paper: Jzz = 0 is Eq. 1, Jzz > 0 is Eq. 23. The dual rotation in the paper's boundary section is a conceptual duality, not a separate model, so the package carries no separate type for it, and the evolver consumes the coupling and field terms without branching on the parameters.

## How we know it is right

The MPS evolution is checked against the exact engine on every change. At small N the test evolves the same chain on a full statevector under the same brickwall, with no truncation, and checks the MPS against it, which pins the operator factor and the sign. exact_caustic runs that statevector brickwall and reaches about 28 spins on a laptop; exact_propagator applies the true propagator exp(-i H dt) on a dense state at small N, which separates the Trotter error from the truncation; compare_caustic reports the difference as a single number. A coauthor reproduced the paper's results independently in Python, the cross-check cited in the manuscript.

A finite chain is the matched method here, not an infinite MPS. The single central flip breaks the translational invariance an infinite MPS assumes, and capturing the spreading excitation would need a unit cell as large as the light cone, which removes the advantage. A chain long enough that the cone never reaches the boundary within the run is sufficient, and the finite-size sweep in Appendix A is the check that it is long enough.
