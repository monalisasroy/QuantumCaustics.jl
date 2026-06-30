# Spin, basis, and phase-transition conventions, stated once for the whole package.
#
# Operators.
#   Gates are built from ITensors spin-1/2 operators, for which S_a = sigma_a / 2 with
#   eigenvalues plus or minus 1/2. The coupling Jxx, the integrability-breaking coupling
#   Jzz, and the field hz are given in the Pauli normalisation of Singh Roy et al. (2026),
#   where
#       H = - Jxx sum_j X_j X_{j+1} - hz sum_j Z_j                          (Eq. 1)
#   and the integrability-broken case adds - Jzz sum_j Z_j Z_{j+1}         (Eq. 23).
#   The spec layer (hamiltonians.jl) records terms with Pauli operator names "X","Z".
#   The evolver converts each Pauli-normalised term to its spin-operator gate coefficient.
#   The conversion factor is pinned against the small-N exact engine in the test suite.
#
# Basis.
#   |0> is spin up, Z|0> = +|0>; |1> is spin down, Z|1> = -|1>, matching the paper.
#   The quench of Eq. 2 is realised by preparing the uniform up background, the paramagnetic
#   ground state of the field term, and flipping the centre site down. The background
#   polarisation and the initial condition are recorded with each run, so the colormap can
#   be checked against Fig. 2.
#
# Phase transition.
#   In this normalisation the transition sits at Jxx = hz. Couplings are reported in the same
#   normalisation throughout, so an exponent that moves at Jxx/hz = 1 is read straight off
#   the input parameters with no rescaling.

const PAULI_OPS = ("X", "Y", "Z")
