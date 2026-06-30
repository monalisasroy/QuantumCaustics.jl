# Exact reference for the MPS evolution. The same brickwall the MPS applies is applied here to
# a full statevector with no truncation, so the difference between the two is the tensor-network
# truncation alone. The statevector holds 2^N complex amplitudes, so this engine reaches about
# 28 spins on a laptop before it runs out of memory, where the MPS reaches several hundred. A
# dense propagator, exp(-i H dt), supplies the true-dynamics reference at small N, against which
# the Trotter error of the brickwall is measured. The gates, the initial state, and the
# observables match the MPS path; only the engine differs.

using LinearAlgebra

const _PAULI_MAT = Dict(
    "X" => ComplexF64[0 1; 1 0],
    "Y" => ComplexF64[0 -im; im 0],
    "Z" => ComplexF64[1 0; 0 -1],
)
const _ID2 = ComplexF64[1 0; 0 1]

# Bond and field gates in the Pauli normalisation, identical to the spin-operator gates the
# MPS builds: a Pauli term with coefficient c contributes exp(-i dt c P). The conversion is
# pinned against the MPS in the test suite.
function _gate_pair(spec::HamiltonianSpec, dt::Float64)
    g4 = zeros(ComplexF64, 4, 4)
    for (A, B, c) in coupling_terms(spec)
        g4 .+= c .* kron(_PAULI_MAT[A], _PAULI_MAT[B])
    end
    g2 = zeros(ComplexF64, 2, 2)
    for (A, c) in field_terms(spec)
        g2 .+= c .* _PAULI_MAT[A]
    end
    return exp(-im * dt * g4), exp(-im * dt * g2)
end

# Statevector kernel. Site j is 1-indexed from the left and sits at bit position N - j of the
# basis index, so site 1 is the most significant bit, matching the Kronecker order of the dense
# reference. Gates are applied in place, so only one 2^N vector is held.

# One-site gate on site j: pair each basis state with its partner under flipping bit N - j.
function _apply1!(psi::Vector{ComplexF64}, g, j::Int, N::Int)
    p = N - j
    stride = 1 << p
    step = 1 << (p + 1)
    for base in 0:step:((1 << N) - 1)
        for off in 0:(stride - 1)
            k = base + off
            kp = k + stride
            a = psi[k + 1]
            b = psi[kp + 1]
            psi[k + 1] = g[1, 1] * a + g[1, 2] * b
            psi[kp + 1] = g[2, 1] * a + g[2, 2] * b
        end
    end
end

# Two-site gate on adjacent sites (j, j+1): site j is the higher bit (N - j), site j+1 the lower
# (N - j - 1). The 4-by-4 gate is built as kron(site j, site j+1), so its local index is
# 2 * bit_j + bit_{j+1}.
function _apply2!(psi::Vector{ComplexF64}, g, j::Int, N::Int)
    p2 = N - j - 1
    s2 = 1 << p2
    s1 = s2 << 1
    period = s1 << 1
    for base in 0:period:((1 << N) - 1)
        for off in 0:(s2 - 1)
            k = base + off
            k00 = k; k01 = k + s2; k10 = k + s1; k11 = k + s1 + s2
            v0 = psi[k00 + 1]; v1 = psi[k01 + 1]; v2 = psi[k10 + 1]; v3 = psi[k11 + 1]
            psi[k00 + 1] = g[1, 1] * v0 + g[1, 2] * v1 + g[1, 3] * v2 + g[1, 4] * v3
            psi[k01 + 1] = g[2, 1] * v0 + g[2, 2] * v1 + g[2, 3] * v2 + g[2, 4] * v3
            psi[k10 + 1] = g[3, 1] * v0 + g[3, 2] * v1 + g[3, 3] * v2 + g[3, 4] * v3
            psi[k11 + 1] = g[4, 1] * v0 + g[4, 2] * v1 + g[4, 3] * v2 + g[4, 4] * v3
        end
    end
end

# <Z_j> for every site, read from the amplitudes by the parity of bit N - j: bit 0 is spin up
# with Z = +1, bit 1 is spin down with Z = -1.
function _z_expect(psi::Vector{ComplexF64}, N::Int)
    Z = zeros(Float64, N)
    for k in 0:((1 << N) - 1)
        p = abs2(psi[k + 1])
        p == 0.0 && continue
        for j in 1:N
            bit = (k >> (N - j)) & 1
            Z[j] += p * (1 - 2 * bit)
        end
    end
    return Z
end

# Half-chain von Neumann entropy in bits across the cut to the left of site j0, from the
# singular values of the statevector reshaped at that cut.
function _exact_entanglement(psi::AbstractVector, j0::Int, N::Int)
    m = j0 - 1
    A = reshape(psi, 2^(N - m), 2^m)
    s = 0.0
    for lam in svdvals(A)
        q = lam^2
        q > 0 && (s -= q * log2(q))
    end
    return s
end

"""
    exact_caustic(spec; N, dt, ttotal, init=:flip)

Evolve the same chain on a full statevector under the brickwall the MPS applies, with no
truncation, and return a CausticRun. With init = :flip the centre site starts down, realising
the quench of Eq. 2; with init = :noflip it starts up. The statevector holds 2^N complex
amplitudes, so the engine reaches about 28 spins on a laptop before it runs out of memory.
"""
function exact_caustic(spec::HamiltonianSpec; N::Int, dt::Float64 = DEFAULT_DT,
                       ttotal::Float64 = DEFAULT_TTOTAL, init::Symbol = :flip)
    isodd(N) || error("N should be odd so the chain has a centre site")
    N <= 30 ||
        error("exact statevector evolution allocates 2^N complex amplitudes (2^28 is about 4 GB); N = $N is beyond a laptop's memory (try N <= 28)")
    boundary(spec) == OPEN || error("the exact path implements open boundaries, matching the MPS path")
    (init == :flip || init == :noflip) || error("init should be :flip or :noflip (got $init)")
    g_bond, g_field = _gate_pair(spec, dt)
    j0 = div(N + 1, 2)
    psi = zeros(ComplexF64, 1 << N)
    psi[(init == :flip ? (1 << (N - j0)) : 0) + 1] = 1.0
    nsteps = round(Int, ttotal / dt)
    times = collect(0.0:dt:(nsteps * dt))
    Z = Matrix{Float64}(undef, N, nsteps + 1)
    S = Vector{Float64}(undef, nsteps + 1)
    Z[:, 1] .= _z_expect(psi, N)
    S[1] = _exact_entanglement(psi, j0, N)
    for step in 1:nsteps
        for j in 1:2:(N - 1)
            _apply2!(psi, g_bond, j, N)
        end
        for j in 2:2:(N - 1)
            _apply2!(psi, g_bond, j, N)
        end
        for j in 1:N
            _apply1!(psi, g_field, j, N)
        end
        Z[:, step + 1] .= _z_expect(psi, N)
        S[step + 1] = _exact_entanglement(psi, j0, N)
    end
    return CausticRun(spec, N, j0, dt, times, Z, S, init)
end

# --- Dense true-dynamics reference, for the Trotter-error check at small N -----------------

_embed1(mat, site, N) = reduce(kron, [k == site ? mat : _ID2 for k in 1:N])

function _hamiltonian(spec::HamiltonianSpec, N::Int)
    H = zeros(ComplexF64, 2^N, 2^N)
    for (A, B, c) in coupling_terms(spec), j in 1:(N - 1)
        H .+= c .* (_embed1(_PAULI_MAT[A], j, N) * _embed1(_PAULI_MAT[B], j + 1, N))
    end
    for (A, c) in field_terms(spec), j in 1:N
        H .+= c .* _embed1(_PAULI_MAT[A], j, N)
    end
    return H
end

"""
    exact_propagator(spec; N, dt, ttotal, init=:flip)

Evolve the same chain under the true propagator exp(-i H dt), built as a dense 2^N by 2^N
operator, the dynamics the brickwall approximates. This is the reference against which the Trotter error is
measured: the brickwall agrees with it as dt decreases. Dense, so feasible for N up to about 12.
"""
function exact_propagator(spec::HamiltonianSpec; N::Int, dt::Float64 = DEFAULT_DT,
                          ttotal::Float64 = DEFAULT_TTOTAL, init::Symbol = :flip)
    isodd(N) || error("N should be odd so the chain has a centre site")
    N <= 14 || error("the dense propagator builds a 2^N by 2^N operator; N = $N is too large (try N <= 12)")
    boundary(spec) == OPEN || error("the exact path implements open boundaries, matching the MPS path")
    (init == :flip || init == :noflip) || error("init should be :flip or :noflip (got $init)")
    j0 = div(N + 1, 2)
    up = ComplexF64[1, 0]; dn = ComplexF64[0, 1]
    states = [up for _ in 1:N]
    init == :flip && (states[j0] = dn)
    psi = reduce(kron, states)
    U = exp(-im * dt * _hamiltonian(spec, N))
    nsteps = round(Int, ttotal / dt)
    times = collect(0.0:dt:(nsteps * dt))
    Z = Matrix{Float64}(undef, N, nsteps + 1)
    S = Vector{Float64}(undef, nsteps + 1)
    Zops = [_embed1(_PAULI_MAT["Z"], j, N) for j in 1:N]
    function record!(col, ket)
        for j in 1:N
            Z[j, col] = real(ket' * (Zops[j] * ket))
        end
        S[col] = _exact_entanglement(ket, j0, N)
    end
    record!(1, psi)
    for step in 1:nsteps
        psi = U * psi
        record!(step + 1, psi)
    end
    return CausticRun(spec, N, j0, dt, times, Z, S, init)
end

# --- Difference and comparison ------------------------------------------------------------

"""
    exact_difference(spec; N, dt, ttotal)

Exact analogue of caustic_difference: the local magnetisation difference dZ of Eq. 22 from the
statevector brickwall, the no-flip run minus the flip run. This dZ and the MPS dZ agree to
truncation. The test suite runs this check on the subtracted observable.
"""
function exact_difference(spec::HamiltonianSpec; N::Int, dt::Float64 = DEFAULT_DT,
                          ttotal::Float64 = DEFAULT_TTOTAL)
    flip = exact_caustic(spec; N = N, dt = dt, ttotal = ttotal, init = :flip)
    noflip = exact_caustic(spec; N = N, dt = dt, ttotal = ttotal, init = :noflip)
    dZ = noflip.Z .- flip.Z
    return (run = flip, dZ = dZ, times = flip.times)
end

"""
    compare_caustic(spec; N, dt, ttotal, init=:flip, maxdim, cutoff)

Run the MPS evolver and the exact statevector engine on the same chain and return both runs
with the agreement: the maximum and mean absolute difference in <Z_j(t)>, and the maximum
difference in the half-chain entanglement. Both apply the same brickwall, so the difference is
the MPS truncation, which at small N is none.
"""
function compare_caustic(spec::HamiltonianSpec; N::Int, dt::Float64 = DEFAULT_DT,
                         ttotal::Float64 = DEFAULT_TTOTAL, init::Symbol = :flip,
                         maxdim::Int = DEFAULT_MAXDIM, cutoff::Float64 = DEFAULT_CUTOFF)
    mps = evolve_caustic(spec; N = N, dt = dt, ttotal = ttotal, maxdim = maxdim, cutoff = cutoff, init = init)
    ex = exact_caustic(spec; N = N, dt = dt, ttotal = ttotal, init = init)
    dZ = abs.(mps.Z .- ex.Z)
    dS = abs.(mps.S .- ex.S)
    return (mps = mps, exact = ex,
            dZmax = maximum(dZ), dZmean = sum(dZ) / length(dZ), dSmax = maximum(dS))
end
