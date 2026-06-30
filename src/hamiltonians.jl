# Hamiltonian specification: one type behind one interface.
#
# The paper is one Hamiltonian. Equation 1 is the transverse-field Ising model; Equation 23
# adds a Z Z coupling to that same Hamiltonian to break integrability. So one type carries
# Jxx, hz, and an optional Jzz: Jzz = 0 is Eq. 1, Jzz > 0 is Eq. 23. The evolver and the
# exact reference both consume the same two functions, coupling_terms and field_terms, so
# neither branches on the parameters. Terms are returned in the Pauli
# normalisation of conventions.jl.

"Boundary condition for the chain."
@enum Boundary OPEN PERIODIC

"Supertype for a Hamiltonian specification."
abstract type HamiltonianSpec end

"""
    TransverseFieldIsing(; Jxx, hz, Jzz=0.0, boundary=OPEN)

Transverse-field Ising model with an optional integrability-breaking Z Z coupling:

    H = - Jxx sum_j X_j X_{j+1} - Jzz sum_j Z_j Z_{j+1} - hz sum_j Z_j.

Jzz = 0 is the integrable TFIM of Eq. 1, driven across the quantum phase transition at
Jxx = hz. Jzz > 0 is the integrability-broken model of Eq. 23, in which, for the z-polarised
initial state, Jxx still drives the transition and Jzz acts as the perturbation.
"""
struct TransverseFieldIsing <: HamiltonianSpec
    Jxx::Float64
    Jzz::Float64
    hz::Float64
    boundary::Boundary
end
TransverseFieldIsing(; Jxx, hz, Jzz = 0.0, boundary::Boundary = OPEN) =
    TransverseFieldIsing(Jxx, Jzz, hz, boundary)

# --- Shared interface -------------------------------------------------------------
# coupling_terms: two-site terms as (op_a, op_b, coefficient), Pauli-normalised.
# field_terms:    one-site terms as (op, coefficient), Pauli-normalised.
# The Z Z term is present only when Jzz is nonzero, so the integrable case carries no
# zero-coefficient gate.

function coupling_terms(s::TransverseFieldIsing)
    terms = [("X", "X", -s.Jxx)]
    s.Jzz != 0.0 && push!(terms, ("Z", "Z", -s.Jzz))
    return terms
end

field_terms(s::TransverseFieldIsing) = [("Z", -s.hz)]

"Boundary condition of a spec."
boundary(s::HamiltonianSpec) = s.boundary
