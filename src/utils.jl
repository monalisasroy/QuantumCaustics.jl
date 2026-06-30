# Utilities a finished run needs: the run tag, the disk writers, the colour map, and the
# loader that reads a problem instance from a file. Each run writes two observables, the raw
# magnetisation Z and the subtracted caustic dZ. The data goes to the output directory as
# "<label>_<tag>.txt", and the figure to its figures/ subdirectory as "<label>_<tag>.png", so
# data and figure share one stem. Kept apart from the physics core; only this layer needs
# Printf, Plots, and TOML.

using Printf
using TOML
import Plots

function run_tag(run::CausticRun)
    s = run.spec
    parts = ["tfim", @sprintf("N%d", run.N), @sprintf("Jxx%.3f", s.Jxx)]
    s.Jzz != 0.0 && push!(parts, @sprintf("Jzz%.3f", s.Jzz))
    push!(parts, @sprintf("hz%.3f", s.hz), @sprintf("dt%.3f", run.dt), String(run.init))
    return join(parts, "_")
end

_write_matrix(path, M) = open(path, "w") do f
    for i in 1:size(M, 1)
        println(f, join(M[i, :], " "))
    end
end

"""
    colormap(run; data=run.Z, path=nothing, clabel="Z")

Heatmap of a data matrix, sites against time. Saves to `path` when given and returns the plot
object. This is the small version of Fig. 2.
"""
function colormap(run::CausticRun; data::AbstractMatrix = run.Z, path = nothing, clabel = "Z")
    p = Plots.heatmap(run.times, 1:run.N, data;
        xlabel = "time", ylabel = "site j", colorbar_title = clabel,
        title = @sprintf("TFIM, N=%d", run.N))
    path !== nothing && Plots.savefig(p, path)
    return p
end

"""
    write_caustic(res, dir) -> tag

Write both observables of one run into `dir` and their figures into `dir/figures/`: the raw
magnetisation "Z_<tag>", the measured light cone in which the boundary's reversed cone is
visible, and the subtracted caustic "dZ_<tag>", in which the boundary cancels and the 2/3
scaling lives. The data is a .txt in `dir`, the figure a matching .png in `dir/figures/`. The
shared time axis and a parameter summary are written once. `res` is the result of
caustic_difference. Returns the tag.
"""
function write_caustic(res, dir::AbstractString)
    run = res.run
    figdir = joinpath(dir, "figures")
    mkpath(figdir)
    tag = run_tag(run)
    _write_matrix(joinpath(dir, "Z_$(tag).txt"), run.Z)
    colormap(run; data = run.Z, path = joinpath(figdir, "Z_$(tag).png"), clabel = "Z")
    _write_matrix(joinpath(dir, "dZ_$(tag).txt"), res.dZ)
    colormap(run; data = res.dZ, path = joinpath(figdir, "dZ_$(tag).png"), clabel = "dZ")
    open(joinpath(dir, "times_$(tag).txt"), "w") do f
        println(f, join(run.times, " "))
    end
    open(joinpath(dir, "summary_$(tag).txt"), "w") do f
        s = run.spec
        println(f, "# QuantumCaustics run summary")
        println(f, "model tfim")
        println(f, "N ", run.N)
        println(f, "j0 ", run.j0)
        println(f, "dt ", run.dt)
        println(f, "ttotal ", last(run.times))
        println(f, "Jxx ", s.Jxx)
        println(f, "Jzz ", s.Jzz)
        println(f, "hz ", s.hz)
    end
    return tag
end

"""
    load_problem(path) -> NamedTuple

Read a problem instance from a file: the Hamiltonian and run parameters, with any missing key
filled from the defaults in parameters.jl. The file is key = value lines, the format of the
instances in problems/. Returns a named tuple with N, Jxx, hz, Jzz, dt, ttotal, maxdim,
cutoff, boundary, and description.
"""
function load_problem(path::AbstractString)
    t = TOML.parsefile(path)
    return (N = Int(get(t, "N", DEFAULT_N)),
            Jxx = Float64(get(t, "Jxx", DEFAULT_JXX)),
            hz = Float64(get(t, "hz", DEFAULT_HZ)),
            Jzz = Float64(get(t, "Jzz", DEFAULT_JZZ)),
            dt = Float64(get(t, "dt", DEFAULT_DT)),
            ttotal = Float64(get(t, "T", DEFAULT_TTOTAL)),
            maxdim = Int(get(t, "maxdim", DEFAULT_MAXDIM)),
            cutoff = Float64(get(t, "cutoff", DEFAULT_CUTOFF)),
            boundary = String(get(t, "boundary", DEFAULT_BOUNDARY)),
            description = String(get(t, "description", "")))
end
