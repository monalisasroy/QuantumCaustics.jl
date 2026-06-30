using QuantumCaustics
using Test

@testset "QuantumCaustics" begin

    @testset "Hamiltonian spec and shared interface" begin
        s = TransverseFieldIsing(; Jxx = 0.5, hz = 1.0)
        @test boundary(s) == OPEN
        @test field_terms(s) == [("Z", -1.0)]
        # Jzz = 0 is the pure TFIM: no Z Z term, no zero-coefficient gate.
        @test coupling_terms(s) == [("X", "X", -0.5)]

        sb = TransverseFieldIsing(; Jxx = 0.5, Jzz = 0.2, hz = 1.0, boundary = PERIODIC)
        @test boundary(sb) == PERIODIC
        @test coupling_terms(sb) == [("X", "X", -0.5), ("Z", "Z", -0.2)]
    end

    @testset "MPS evolution against the exact engine" begin
        # The cross-check that the CLI does not need: validate the MPS evolution against the
        # exact engine at small N. This pins the operator factor and the Hamiltonian sign.
        spec = TransverseFieldIsing(; Jxx = 0.5, hz = 1.0)
        N, dt, T = 5, 0.05, 1.0
        mps = evolve_caustic(spec; N = N, dt = dt, ttotal = T, maxdim = 64, cutoff = 1e-12)

        # Against the exact statevector brickwall, the same gates with no truncation: agreement
        # to truncation only, which at this size is machine precision.
        brick = exact_caustic(spec; N = N, dt = dt, ttotal = T)
        @test maximum(abs.(mps.Z .- brick.Z)) < 1e-8

        # Against the true dynamics exp(-iHdt) from the dense propagator: agreement to within
        # the first-order Trotter error.
        prop = exact_propagator(spec; N = N, dt = dt, ttotal = T)
        @test maximum(abs.(mps.Z .- prop.Z)) < 5e-2

        # compare_caustic returns the same agreement, in Z and in the entanglement.
        cmp = compare_caustic(spec; N = N, dt = dt, ttotal = T, maxdim = 64, cutoff = 1e-12)
        @test cmp.dZmax < 1e-8
        @test cmp.dSmax < 1e-8
    end

    @testset "Caustic difference, MPS against exact" begin
        # The subtracted dZ of Eq. 22 from the MPS path matches the exact brickwall dZ to
        # truncation, both for the integrable case and with a nonzero Z Z coupling.
        N, dt, T = 5, 0.05, 1.0
        for Jzz in (0.0, 0.2)
            spec = TransverseFieldIsing(; Jxx = 0.5, Jzz = Jzz, hz = 1.0)
            mps = caustic_difference(spec; N = N, dt = dt, ttotal = T, maxdim = 64, cutoff = 1e-12)
            ex = exact_difference(spec; N = N, dt = dt, ttotal = T)
            @test maximum(abs.(mps.dZ .- ex.dZ)) < 1e-8
        end
    end

    @testset "Problem instances load" begin
        dir = joinpath(@__DIR__, "..", "problems")
        para = load_problem(joinpath(dir, "tfim_paramagnetic.toml"))
        @test (para.Jxx, para.hz, para.Jzz) == (0.5, 1.0, 0.0)
        @test para.boundary == "open"
        @test !isempty(para.description)

        ferro = load_problem(joinpath(dir, "tfim_ferromagnetic.toml"))
        @test ferro.Jxx == 1.5

        weak = load_problem(joinpath(dir, "tfim_paramagnetic_jzz_weak.toml"))
        @test weak.Jzz == 0.2

        # A misspelled key is rejected rather than silently defaulted.
        bad = tempname() * ".toml"
        write(bad, "N = 9\nJx = 1.5\n")
        @test_throws ErrorException load_problem(bad)
        rm(bad; force = true)

        # A missing optional key still falls back to its default.
        minimal = tempname() * ".toml"
        write(minimal, "N = 9\nJxx = 0.5\nhz = 1.0\n")
        m = load_problem(minimal)
        @test m.Jzz == DEFAULT_JZZ
        @test m.maxdim == DEFAULT_MAXDIM
        rm(minimal; force = true)
    end

    @testset "Defaults" begin
        @test DEFAULT_DT == 0.1
        @test DEFAULT_JZZ == 0.0
        @test DEFAULT_TROTTER_ORDER == 1
    end

end
