module QuantumCaustics

# Files are included in dependency order, so the module file is the dependency graph.
# parameters and hamiltonians use only Base. observables and evolution add
# the MPS machinery; exact is the dense reference used to validate the evolution in the test
# suite; utils is the writing, plotting, and instance-loading layer built on a completed run.

include("parameters.jl")
include("hamiltonians.jl")
include("observables.jl")
include("evolution.jl")
include("exact.jl")
include("utils.jl")

export Boundary, OPEN, PERIODIC,
    HamiltonianSpec, TransverseFieldIsing,
    coupling_terms, field_terms, boundary,
    DEFAULT_N, DEFAULT_BOUNDARY, DEFAULT_JXX, DEFAULT_HZ, DEFAULT_JZZ,
    DEFAULT_DT, DEFAULT_TTOTAL, DEFAULT_CUTOFF, DEFAULT_MAXDIM, DEFAULT_TROTTER_ORDER,
    DEFAULT_OUTPUT_DIR,
    magnetisation, entanglement_entropy,
    CausticRun, build_gates, evolve_caustic, caustic_difference,
    exact_caustic, exact_propagator, exact_difference, compare_caustic,
    run_tag, colormap, write_caustic, load_problem

end # module QuantumCaustics
