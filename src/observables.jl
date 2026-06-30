# Observables read off an evolved MPS: magnetisation and half-chain entanglement.
# These are the quantities the figures and the scaling analysis consume.

using ITensors
using ITensorMPS

"""
    magnetisation(psi)

Return the Pauli magnetisation <Z_j> for every site as a vector in [-1, 1].
ITensors reports the spin operator <Sz> = <Z>/2, so the value is doubled here, once,
so the rest of the package works in the Pauli normalisation of the paper.
"""
magnetisation(psi::MPS) = 2.0 .* expect(psi, "Sz")

"""
    entanglement_entropy(psi, b)

Half-chain von Neumann entropy in bits across the bond to the left of site `b`, from the
singular values of the orthogonality centre. Used for the convergence study of
Appendix A. The chain is odd-length with its centre at b > 1, so the left link exists.
"""
function entanglement_entropy(psi::MPS, b::Int)
    psic = orthogonalize!(copy(psi), b)
    lb = linkind(psic, b - 1)
    sb = siteind(psic, b)
    _, S, _ = svd(psic[b], (lb, sb))
    SvN = 0.0
    for n in 1:dim(S, 1)
        p = S[n, n]^2
        p > 0 && (SvN -= p * log2(p))
    end
    return SvN
end
