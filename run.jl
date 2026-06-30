#!/usr/bin/env julia
# CLI entry point over the QuantumCaustics library. The library holds the logic; this file
# resolves parameters and dispatches. Two ways to set parameters: a named problem instance
# read from problems/, or parameter flags. Each run writes both observables, the raw
# magnetisation Z and the subtracted caustic dZ. Run with --help for the full option list.

const USAGE = """
QuantumCaustics run.jl: evolve the central-flip quench and write the caustic and its data.

Usage:
  julia --project=. run.jl --problem NAME [NAME ...] [--output-dir PATH]
  julia --project=. run.jl [parameter flags] [--output-dir PATH]

Two ways to choose parameters:

1. A named problem instance, read from problems/<NAME>.toml. The instances are the regimes of
   the paper; --problem all runs every one.

     julia --project=. run.jl --problem tfim_paramagnetic
     julia --project=. run.jl --problem tfim_paramagnetic tfim_ferromagnetic
     julia --project=. run.jl --problem all

2. Parameter flags. Every flag accepts one or more values, and the run repeats for each
   combination, so a sweep is one call.

     --N INT [...]          number of spins, odd                 (default 79)
     --Jxx FLOAT [...]      Ising coupling                       (default 0.5)
     --hz FLOAT [...]       transverse field                     (default 1.0)
     --Jzz FLOAT [...]      integrability-breaking ZZ coupling   (default 0.0, the pure TFIM)
     --dt FLOAT [...]       Trotter step                         (default 0.1)
     --T FLOAT [...]        total evolution time                 (default 40.0)
     --maxdim INT [...]     MPS bond-dimension cap               (default 256)
     --cutoff FLOAT [...]   MPS truncation cutoff                (default 1e-12)
     --boundary NAME [...]  open | periodic                      (default open)
     --output-dir PATH      directory for the run files          (default results)
     -h, --help             show this message and exit

     julia --project=. run.jl --N 79 --Jxx 0.5 --hz 1.0 --T 40 --output-dir results/
     julia --project=. run.jl --N 79 --Jxx 0.4 0.6 0.8 1.0 1.2 --hz 1.0 --output-dir results/sweep

Each run writes two observables as a .txt in the output directory, results/ by default, with a
matching .png in its figures/ subdirectory, results/figures/ by default: the raw magnetisation
"Z_<tag>", the measured light cone in which the boundary's reversed cone is visible, and the
subtracted caustic "dZ_<tag>", in which the boundary cancels and the 2/3 scaling lives. Every flag default above is defined in
src/parameters.jl; edit that file to change a default for every run. Each run logs the
parameters it resolved, with the centre site, the step count, and the phase, as a double-check
before it computes. Validation of the
MPS evolution against the exact engine is in the test suite.
"""

if any(a -> a == "--help" || a == "-h", ARGS)
    print(USAGE)
    exit(0)
end

using QuantumCaustics
using Printf

# Each flag collects the tokens up to the next flag, so any flag may take several values.
function parse_args(args)
    opts = Dict{String,Vector{String}}()
    i = 1
    while i <= length(args)
        tok = args[i]
        if startswith(tok, "--")
            key = tok[3:end]
            vals = String[]
            i += 1
            while i <= length(args) && !startswith(args[i], "--")
                push!(vals, args[i]); i += 1
            end
            opts[key] = isempty(vals) ? ["true"] : vals
        else
            i += 1
        end
    end
    return opts
end

# Each parameter resolves to a list: the values given, or the single default from parameters.jl.
list_float(opts, k, d) = haskey(opts, k) ? parse.(Float64, opts[k]) : [Float64(d)]
list_int(opts, k, d)   = haskey(opts, k) ? parse.(Int, opts[k]) : [Int(d)]
list_str(opts, k, d)   = haskey(opts, k) ? opts[k] : [String(d)]
gets(opts, k, d)       = haskey(opts, k) ? first(opts[k]) : d

# Run one parameter set on the MPS engine and write the two observables. The resolved
# parameters and a few derived quantities are logged first, as a double-check on what is run.
function run_one(; N, Jxx, hz, Jzz, dt, ttotal, maxdim, cutoff, boundary, outdir)
    bc = boundary == "open" ? OPEN : PERIODIC
    spec = TransverseFieldIsing(; Jxx = Jxx, hz = hz, Jzz = Jzz, boundary = bc)
    phase = Jxx < hz ? "paramagnetic" : Jxx > hz ? "ferromagnetic" : "critical"
    j0 = div(N + 1, 2)
    nsteps = round(Int, ttotal / dt)
    @info "running" N Jxx hz Jzz dt T=ttotal maxdim cutoff boundary phase j0 nsteps
    res = caustic_difference(spec; N = N, dt = dt, ttotal = ttotal, maxdim = maxdim, cutoff = cutoff)
    tag = write_caustic(res, outdir)
    @info "wrote" tag directory=outdir
    return tag
end

# Resolve --problem names to instance file paths under problems/.
function problem_paths(names, dir)
    paths = String[]
    for name in names
        if name == "all"
            for f in sort(readdir(dir))
                endswith(f, ".toml") && push!(paths, joinpath(dir, f))
            end
        else
            p = endswith(name, ".toml") ? name : joinpath(dir, name * ".toml")
            isfile(p) || error("problem instance not found: $p")
            push!(paths, p)
        end
    end
    return paths
end

function main()
    opts = parse_args(ARGS)
    outdir = gets(opts, "output-dir", DEFAULT_OUTPUT_DIR)

    if haskey(opts, "problem")
        dir = joinpath(@__DIR__, "problems")
        for path in problem_paths(opts["problem"], dir)
            pr = load_problem(path)
            @info "instance" name=basename(path)
            run_one(; N = pr.N, Jxx = pr.Jxx, hz = pr.hz, Jzz = pr.Jzz, dt = pr.dt,
                    ttotal = pr.ttotal, maxdim = pr.maxdim, cutoff = pr.cutoff,
                    boundary = pr.boundary, outdir = outdir)
        end
        return
    end

    Ns       = list_int(opts, "N", DEFAULT_N)
    Jxxs     = list_float(opts, "Jxx", DEFAULT_JXX)
    hzs      = list_float(opts, "hz", DEFAULT_HZ)
    Jzzs     = list_float(opts, "Jzz", DEFAULT_JZZ)
    dts      = list_float(opts, "dt", DEFAULT_DT)
    Ts       = list_float(opts, "T", DEFAULT_TTOTAL)
    maxdims  = list_int(opts, "maxdim", DEFAULT_MAXDIM)
    cutoffs  = list_float(opts, "cutoff", DEFAULT_CUTOFF)
    boundaries = list_str(opts, "boundary", DEFAULT_BOUNDARY)

    for N in Ns, bc_s in boundaries, Jxx in Jxxs, hz in hzs, Jzz in Jzzs,
        dt in dts, T in Ts, maxdim in maxdims, cutoff in cutoffs

        run_one(; N = N, Jxx = Jxx, hz = hz, Jzz = Jzz, dt = dt, ttotal = T,
                maxdim = maxdim, cutoff = cutoff, boundary = bc_s, outdir = outdir)
    end
end

main()
