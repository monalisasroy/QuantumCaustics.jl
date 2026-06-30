# Time evolution: one first-order brickwall Trotter evolver over the shared interface,
# plus the exact statevector and dense-propagator references in exact.jl, which pin the
# operator convention against the evolver. The evolver consumes coupling_terms and field_terms, so it does not
# branch on the model.

using ITensors
using ITensorMPS

# Pauli operator name to ITensors spin-1/2 operator. Since S_a = sigma_a / 2, a one-site
# Pauli term carries a factor 2 and a two-site Pauli term a factor 4 when written with
# spin operators. These factors and the Hamiltonian sign are checked against the dense
# oracle in the test suite.
const _PAULI_TO_SPIN = Dict("X" => "Sx", "Y" => "Sy", "Z" => "Sz")

function _bond_gate(terms, sa, sb, dt)
    h = ITensor()
    first = true
    for (A, B, c) in terms
        contrib = (4.0 * c) * op(_PAULI_TO_SPIN[A], sa) * op(_PAULI_TO_SPIN[B], sb)
        h = first ? contrib : h + contrib
        first = false
    end
    return exp(-im * dt * h)
end

function _field_gate(terms, sa, dt)
    h = ITensor()
    first = true
    for (A, c) in terms
        contrib = (2.0 * c) * op(_PAULI_TO_SPIN[A], sa)
        h = first ? contrib : h + contrib
        first = false
    end
    return exp(-im * dt * h)
end

"""
    build_gates(spec, sites, dt)

Return the three brickwall layers as named tuples of ITensor gates: coupling on odd
bonds, coupling on even bonds, and the on-site field applied once per site. Applying the
field once, separate from the bond parity, corrects a double-application of the field
present in the exploratory scripts.
"""
function build_gates(spec::HamiltonianSpec, sites, dt)
    boundary(spec) == OPEN ||
        error("Only open boundaries are implemented in v0.1; this matches the paper's numerics. Periodic wrap is reserved.")
    N = length(sites)
    couplings = coupling_terms(spec)
    fields = field_terms(spec)
    odd = ITensor[]
    even = ITensor[]
    field = ITensor[]
    for j in 1:2:(N - 1)
        push!(odd, _bond_gate(couplings, sites[j], sites[j + 1], dt))
    end
    for j in 2:2:(N - 1)
        push!(even, _bond_gate(couplings, sites[j], sites[j + 1], dt))
    end
    for j in 1:N
        push!(field, _field_gate(fields, sites[j], dt))
    end
    return (odd = odd, even = even, field = field)
end

"A single evolved run: the magnetisation colormap, the entanglement trace, and the setup."
struct CausticRun
    spec::HamiltonianSpec
    N::Int
    j0::Int
    dt::Float64
    times::Vector{Float64}
    Z::Matrix{Float64}     # N rows (sites) by length(times) columns
    S::Vector{Float64}     # half-chain entanglement vs time
    init::Symbol
end

"""
    evolve_caustic(spec; N, dt, ttotal, maxdim, cutoff, init)

Evolve an N-site chain under `spec` and record <Z_j(t)> and the half-chain entanglement.
The background is polarised up, the paramagnetic ground state of the field term. With
init = :flip a single central down excitation is injected, realising the quench of Eq. 2.
N is odd so the chain has a centre. The evolution is a first-order Trotter brickwall; the
bond and field layers do not commute, so the error per step is of order dt^2, and dt = 0.1
is the converged choice of the paper.
"""
function evolve_caustic(spec::HamiltonianSpec; N::Int,
                        dt::Float64 = DEFAULT_DT, ttotal::Float64 = DEFAULT_TTOTAL,
                        maxdim::Int = DEFAULT_MAXDIM, cutoff::Float64 = DEFAULT_CUTOFF,
                        init::Symbol = :flip)
    isodd(N) || error("N should be odd so the chain has a centre site")
    (init == :flip || init == :noflip) || error("init should be :flip or :noflip (got $init)")
    sites = siteinds("S=1/2", N; conserve_qns = false)
    j0 = div(N + 1, 2)
    # The polarised reference is a chi = 1 product state along the field axis. The flip run
    # is the same product state with the centre site set to Dn, built directly so the
    # initial state carries no truncation. The no-flip companion leaves the centre Up.
    states = fill("Up", N)
    init == :flip && (states[j0] = "Dn")
    psi = MPS(sites, states)
    gates = build_gates(spec, sites, dt)
    nsteps = round(Int, ttotal / dt)
    times = collect(0.0:dt:(nsteps * dt))
    Z = Matrix{Float64}(undef, N, nsteps + 1)
    S = Vector{Float64}(undef, nsteps + 1)
    Z[:, 1] .= magnetisation(psi)
    S[1] = entanglement_entropy(psi, j0)
    for step in 1:nsteps
        psi = apply(gates.odd, psi; cutoff = cutoff, maxdim = maxdim)
        psi = apply(gates.even, psi; cutoff = cutoff, maxdim = maxdim)
        psi = apply(gates.field, psi; cutoff = cutoff, maxdim = maxdim)
        normalize!(psi)
        Z[:, step + 1] .= magnetisation(psi)
        S[step + 1] = entanglement_entropy(psi, j0)
    end
    return CausticRun(spec, N, j0, dt, times, Z, S, init)
end

"""
    caustic_difference(spec; kwargs...)

Run the flipped and unflipped evolutions and return the local magnetisation difference
dZ of Eq. 22, in which any boundary contribution cancels. The kwargs are passed to
evolve_caustic.
"""
function caustic_difference(spec::HamiltonianSpec; kwargs...)
    flip = evolve_caustic(spec; init = :flip, kwargs...)
    noflip = evolve_caustic(spec; init = :noflip, kwargs...)
    dZ = noflip.Z .- flip.Z
    return (run = flip, dZ = dZ, times = flip.times)
end
