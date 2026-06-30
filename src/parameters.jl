# Default parameters, defined once for the whole package. This is the single place to read or
# change a default.
#
# These are defaults, not fixed quantities. Every entry point takes the same parameters as
# keyword arguments, and run.jl takes them as command-line flags, so each run sets its own
# values: --N, --boundary, --Jxx, --hz, --Jzz, --dt, --T, --maxdim, --cutoff, and
# --output-dir each override the value here for that run. The const keyword fixes only the
# default binding, in one place, so a value is never repeated as a literal across the code;
# the value a run uses is a fresh local computed from the flag or the keyword argument. The
# numerical defaults are the converged choices of Singh Roy et al. (2026), for example
# dt = 0.1, not the dt = 1.0 of the exploratory scripts.

# Geometry.
const DEFAULT_N = 79                    # number of spins, odd so the chain has a centre
const DEFAULT_BOUNDARY = "open"         # open | periodic

# Couplings and field, in the Pauli normalisation of the paper.
const DEFAULT_JXX = 0.5                 # Ising coupling
const DEFAULT_HZ = 1.0                  # transverse field along z
const DEFAULT_JZZ = 0.0                 # integrability-breaking Z Z coupling; 0 is the pure TFIM

# Evolution and truncation.
const DEFAULT_DT = 0.1                  # Trotter step
const DEFAULT_TTOTAL = 40.0             # total evolution time in the same units as dt
const DEFAULT_CUTOFF = 1e-12            # MPS truncation cutoff (singular-value weight)
const DEFAULT_MAXDIM = 256              # bond-dimension cap; the paper sweeps chi up to 181
const DEFAULT_TROTTER_ORDER = 1         # first-order brickwall, matching the paper

# Run output.
const DEFAULT_OUTPUT_DIR = "results"    # directory for the run files
