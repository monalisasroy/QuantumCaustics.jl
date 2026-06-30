# Computing the caustic scaling exponent
A single spin flipped at the centre of a magnetic chain spreads outward as a wave, and at the front of that wave the magnetisation focuses into a sharp, repeating pattern: a caustic, the same focusing that draws the bright lines on the floor of a sunlit swimming pool. The spacing of that pattern shrinks with distance by a fixed power, and that power is a universal number, the same across a whole phase of the material and changing once the material crosses a phase transition. This plan describes how the chain is evolved and the caustic produced, and how that exponent is read from it.

The evolution is done in Julia with ITensor, using a matrix product state. A full exact calculation instead holds the $2^N$ amplitudes of the state and runs out of memory near $N = 28$ on a laptop, too few sites for the light cone to run long enough to resolve the scaling; the matrix product state reaches chains several times longer, at the cost that the answer is now limited by tensor-network truncation, which is largest where the chain is most entangled, late in the evolution and far from the field-dominated limit. The same first-order brickwall that drives the matrix-product-state run is reproduced exactly on small chains by the exact statevector engine, with a dense propagator $e^{-iH\,dt}$ as the true-dynamics reference, so the central result is checked from the inside against more than one method. This version provides the evolution, the measurements, and the caustic; the order parameter and the scaling-exponent fit are left to downstream analysis, computed from the dZ data the run writes, and the sections below describe how. The plan below is written for a reader meeting matrix product states for the first time.

## Scope

Included in this version: evolving the central spin flip on a single open chain with ITensor, measuring the magnetisation light cone and the half-chain entanglement, and producing the caustic colour map, with the exact-engine cross-checks.

Can be easily computed: the peak-delay order parameter $O_{1,2}$ in the paper and the scaling exponent $gamma$ across the paramagnetic phase, at sizes where the method is shown to converge. The package produces the dZ data; the sections below describe how the order parameter and the exponent are computed from it, the analysis a reader runs on the output.

Out of scope: Running the brickwall as a quantum circuit on a simulator or on hardware.

Assumed knowledge: quantum states of spin-1/2 chains, Pauli matrices, and light-cone spreading. Not assumed: matrix product states, ITensor, TEBD, or the Airy caustic, all introduced below.

## 1. What we are computing

The protocol is this: start the chain in the polarised state along the field axis, flip the single spin at the centre, and let it evolve under the Hamiltonian. The flipped spin spreads outward as a magnon, filling a light cone, and at the edge of that cone the amplitude focuses into a caustic, a fold of the wavefront where the magnetisation oscillates in a fixed, universal pattern. We want the exponent that governs how the spacing of those oscillations shrinks with distance.

The chain is $N$ spin-1/2 sites in a row, with Hamiltonian

$$H = -J^{xx} \sum_{i=1}^{N-1} X_i X_{i+1} - h^z \sum_{i=1}^{N} Z_i$$

and open boundaries. Here $X$ and $Z$ are Pauli matrices, $J^{xx}$ sets the strength of the spin-spin coupling and is the swept parameter, and $h^z = 1$ is the transverse field, held fixed. Site 1 is the leftmost spin. A second coupling $-J^{zz} \sum Z_i Z_{i+1}$ can be added to break integrability; it is zero in the pure TFIM and small when we test that the scaling survives integrability breaking.

The single flip is a superposition of momentum modes with dispersion $\varepsilon_k = 2h^z - 2J^{xx}\cos k$, so it spreads at the group velocity $v = 2J^{xx}$ and the light cone is $|j - j_0| = 2J^{xx} t$. Near the cone edge the local magnetisation takes the Airy form

$$\langle Z_j(t)\rangle - 1 \;\propto\; (J^{xx} t)^{-2/3}\, \mathrm{Ai}^2\!\left(\frac{|j - j_0| - 2J^{xx} t}{(J^{xx} t)^{1/3}}\right),$$

and the project rests on reproducing this fold and reading its scaling from the numerics.

We record four quantities along the run:

- Local magnetisation $\langle Z_j(t)\rangle$ and its no-flip subtraction $dZ_j(t)$: the light cone itself (Eq. 22 of the paper). The subtraction removes the evolution of the unflipped background and cancels the open boundary.
- Peak-delay order parameter $O_{12} = (t_2 - t_1)/t_1$ at each site: the normalised delay between the first two magnetisation peaks, which reads the spacing of the first two Airy fringes as the front passes that site.
- Scaling exponent $\gamma$: the slope magnitude of $\log O_{12}$ against $\log|j - j_0|$, the universal number, with theoretical value $2/3$.
- Half-chain entanglement entropy $S$: not a feature of the caustic itself but the quantity the truncation has to resolve, and so the diagnostic that tells us whether a given size is trustworthy.

The magnetisation is what this version produces. The order parameter and the exponent are the downstream analysis, computed from the dZ data. The entanglement is the convergence control, and Section 10 explains why it decides where the method can be used.

## 2. About the Computational State

The result needs no ensemble average and no thermal sampling. The state is pure at every instant, a single deterministic time evolution of one initial condition, so we evolve once and read the observables off the snapshot at each step. The only manipulation is the subtraction: we run a second, identical chain with no flip, and define

$$dZ_j(t) = \langle Z_j(t)\rangle_{\text{flip}} - \langle Z_j(t)\rangle_{\text{no flip}}.$$

On open boundaries the unflipped run carries the same boundary drift as the flipped run, so the difference isolates the contribution of the flipped spin and leaves the caustic clean. Two cheap evolutions, no Monte Carlo.

The chain is also finite by design, not by compromise. A referee may ask why an infinite-MPS method was not used. The quench is a single spin flip at one site, so the light cone it produces is inhomogeneous and breaks the translational invariance that an infinite, translationally invariant MPS assumes. Capturing a spreading single-site excitation in that framework needs a unit cell as large as the cone, which removes the advantage. A finite chain, long enough that the cone never reaches the boundary within the evolution time, is the matched method, and the finite-size sweep is the check that the chain is long enough.

## 3. How ITensor is used in the package

A matrix product state (MPS) writes the state of the whole chain as a row of small tensors, one per site, linked by internal indices called bonds. The size of the largest bond, the bond dimension $\chi$, sets both the storage cost and how much entanglement the state can hold: the polarised product state has $\chi = 1$, and $\chi$ grows as the spins become entangled across a cut. The exact statevector stores all $2^N$ amplitudes and stops near $N = 28$ on a laptop. An MPS stores about $N\chi^2$ numbers and reaches larger $N$ as long as $\chi$ stays bounded. The price is truncation: when the true state needs a larger $\chi$ than we allow, ITensor discards the smallest Schmidt weights at each gate, and the result carries an error set by the caps `maxdim` and `cutoff`.

To read the magnetisation we take $\langle Z_j\rangle$ directly from the MPS, an expectation that costs $O(N\chi^2)$. To read the half-chain entropy we take the Schmidt values at the central cut, which the MPS already holds in its bond, and form $S = -\sum_\alpha s_\alpha^2 \log s_\alpha^2$. Both are nearly free; the cost is in the evolution that grows $\chi$.

## 4. The model in our convention

ITensor has no built-in transverse-field Ising chain in our axes, and the axes matter here: the whole analysis rests on which Pauli operators carry the caustic, so we build the Hamiltonian ourselves. We use `siteinds("S=1/2", N; conserve_qns=false)` for the sites, and an `OpSum` with the $X_i X_{i+1}$ coupling and the $Z_i$ field, in Pauli normalisation, with the sign and factor of the Hamiltonian above. The operator factor and sign are not asserted but pinned by the test suite against the exact engine, so a convention slip cannot pass silently.

We set `conserve_qns=false` to begin. The $X X$ coupling does not conserve the magnetisation along z, so there is no abelian charge to exploit, and the single-flip excitation lives in its own sector without any symmetry bookkeeping. The model does keep a $\mathbb{Z}_2$ parity, the product of all $Z_i$, which the flip changes; leaving the quantum numbers off keeps the flipped and unflipped runs correct without tracking it.

The spin convention is fixed once, here: `"Up"` is spin up along z and `"Dn"` is spin down, and site 1 is the leftmost spin. Open boundaries only in v0.1, matching the paper's numerics, where the $dZ$ subtraction cancels the boundary.

## 5. Building the initial state

The polarised reference is a product state with $\chi = 1$: build it with `MPS(sites, fill("Up", N))`, every spin up along the field axis. This is the high-field ground state to leading order and the vacuum on which the single magnon sits.

The flip run is the same product state with the centre site set to `"Dn"`. With $N$ odd the centre is $j_0 = (N+1)/2$, a single site, so the flip is unambiguous. This state is still $\chi = 1$ and `from_product_state` builds it with no error.

The no-flip companion is the polarised state with no change, evolved alongside the flip run so that the two can be subtracted step for step. Keeping them on the same grid and the same gates is what makes the cancellation exact.

## 6. Evolving one step in time

The Hamiltonian couples only neighbouring sites, so TEBD (Time Evolving Block Decimation) applies. TEBD splits $U = e^{-iH\,dt}$ into a product of two-site gates using a Trotter decomposition: the bond terms $X_i X_{i+1}$ become two-site gates, the field term $Z_i$ folds into the layer, and the gates are ordered as a brickwall, even bonds then odd. For the first-order split the error per step is of order $dt^2\,\lVert[H_{\text{even}}, H_{\text{odd}}]\rVert$, and $dt = 0.1$ is the converged default. Each step is `apply(gates, psi; cutoff, maxdim)`, which sweeps the gates across and truncates as it goes.

The exact engine in `exact.jl` is a runnable part of the package, not hidden in the tests, so the MPS run can be compared directly against exact dynamics on small chains. `exact_caustic` applies the same brickwall to a full state vector with no truncation, so it isolates the MPS truncation, and holding $2^N$ amplitudes it reaches about 28 spins on a laptop. `exact_propagator` applies the true propagator $e^{-iH\,dt}$ on a dense state at small N, the reference that separates the Trotter error from the truncation. `compare_caustic` returns the maximum and mean differences in magnetisation and in entanglement, so the comparison is a number, not a glance at two plots.

## 7. From the time series to the exponent

At each step, read $\langle Z_j(t)\rangle$ for every site off both MPS runs and subtract to get $dZ_j(t)$. At a fixed site, the trace $dZ_j(t)$ stays near zero until the front arrives, then rises and oscillates: the first two local maxima are the times $t_1$ and $t_2$. A peak finder over the dZ trace locates them per site.

Form $O_{12} = (t_2 - t_1)/t_1$ at each site and collect it against the distance $|j - j_0|$. Because $t_n \approx |j - j_0|/(2J^{xx})$ with a correction of order $|j - j_0|^{1/3}$, the delay scales as $O_{12} \propto |j - j_0|^{-2/3}$. Fit $\log O_{12} = \text{const} - \gamma \log|j - j_0|$; the slope magnitude is the exponent $\gamma$, reported with its coefficient of determination $R^2$. The main-text figure keeps only the sites with $R^2 > 0.9$ (Fig. 3(b)); the full unfiltered set is kept for the appendix (Fig. 7(a)), and the relationship between the two is stated rather than left for the reader to infer.

Two conventions to hold throughout. Fix the logarithm base of the entropy once and match it to the stored convergence arrays, so a base mismatch does not masquerade as a discrepancy. And fix the read-off windows once: the half-chain entropy at the central cut, the magnetisation on every site, the peak search starting from the arrival of the front at each site.

## 8. Parameters

All settings live in one place, edited here and nowhere else.

```julia
# Parameters -- edit here
N        = 79         # chain length, odd so the flip sits on the centre site
JXX      = 0.5        # transverse Ising coupling (swept across the transition)
HZ       = 1.0        # transverse field, held fixed at 1
JZZ      = 0.0        # optional ZZ coupling; 0 is the pure TFIM, > 0 breaks integrability
DT       = 0.1        # Trotter sub-step, the converged default
TTOTAL   = 100.0      # total evolution time (see the boundary condition below)
MAXDIM   = 256        # bond-dimension cap (truncation control)
CUTOFF   = 1e-12      # smallest retained Schmidt weight (truncation control)
INIT     = :flip      # :flip evolves the central-flip state; the no-flip companion runs for dZ
BOUNDARY = OPEN       # open boundaries only in v0.1
```

Keep $2 \cdot \text{JXX} \cdot \text{TTOTAL}$ below the distance from the flip to the nearest boundary, $(N-1)/2$, so the light cone never reaches the edge within the run; this is the single constraint that ties $N$ to TTOTAL. `MAXDIM` and `CUTOFF` are the two truncation controls, and the result must be checked for convergence against both (Section 11). The example `JXX = 0.5` with `HZ = 1.0` sits in the paramagnetic phase, $J^{xx} < h^z$, where the $2/3$ scaling holds and is a sensible first operating point.

## 9. The pipeline 

1. Read the parameters.
2. Build the model with the $X X$ coupling and $Z$ field on `S=1/2` sites, in Pauli normalisation.
3. Build the polarised product state; flip the centre site for the flip run and leave it for the companion.
4. Evolve both runs by $dt$ with the brickwall, recording $\langle Z_j\rangle$ at each step, until $t$ reaches TTOTAL.
5. Subtract the two runs to get $dZ_j(t)$.
6. At each site find the first two peaks $t_1$ and $t_2$, and form $O_{12}$.
7. Fit the power law of $O_{12}$ against distance; read off $\gamma$ and its $R^2$.
8. Record the exponent for the current $(J^{xx}, h^z, J^{zz})$ cell, then move to the next cell.

## 10. Discussion

The accuracy is set by how entangled the chain becomes, because that is what the bond-dimension cap limits. Two things drive the entanglement. The first is time: the longer the run, the more the magnon front has dressed the chain, so a longer chain that needs a longer run carries more entanglement at the cut. The second is the coupling: deep in the paramagnetic phase, $J^{xx} \ll h^z$, the dynamics is nearly free and $\chi$ stays small; as $J^{xx}$ approaches $h^z$, the chain entangles faster and the cap is reached sooner.

This sorts the diagnostics by difficulty. The exponent fit reads the front, which stays comparatively weakly entangled, so small $\chi$ suffices and the MPS can run the long chains the fit wants, well past the exact ceiling. The half-chain entropy, by contrast, grows with both $N$ and time and is the quantity that forces $\chi$ up; it is reported mainly as the convergence diagnostic rather than as a target. The exact engine has no truncation and no $\chi$ to converge, which is why it is the matched method for the small chains where it fits, and why the test suite pins the MPS against it there. The MPS earns its place precisely where the exact engine cannot reach, on the long, weakly entangled fronts that the scaling fit needs.

One caveat carries through. The fitted exponent converges to about $0.68$ against the theoretical $2/3$. The origin of the residual gap is not settled: it may be finite $N$, finite bond dimension, or a true many-body effect. The convergence data are reported with the gap visible.

## 11. On Results' Accuracy

Before trusting any new size, reproduce the exact values at a few small chains: the MPS pipeline should match the exact-engine magnetisation and entropy to the tolerance the test suite sets, which `compare_caustic` reports as a single maximum and mean difference rather than a visual overlay.

As running checks while the pipeline executes:

- Confirm the polarised and flipped states are normalised and start at $\chi = 1$, and that $\chi$ grows smoothly rather than jumping to the cap, which would signal that `maxdim` is binding too early.
- Confirm the light cone stays inside the chain: the front at $2J^{xx} t$ must not reach the boundary within TTOTAL, or the $dZ$ subtraction no longer cancels it.
- Confirm the statevector brickwall and the true-dynamics propagator agree at one operating point, through `compare_caustic` and `exact_propagator`, which separates Trotter error from truncation error.
- Confirm convergence in both `maxdim` and `cutoff` before any large-$N$ exponent is reported, the sweep that the appendix figures document.

On a first run on a fresh checkout, three ITensor touch points are worth watching against the exact oracle: the gate construction in the brickwall, the `apply` call that truncates, and the expectation-value call that reads $\langle Z_j\rangle$. If those three agree with the exact engine on a small chain, the rest of the pipeline is arithmetic on numbers the MPS has already got right.
