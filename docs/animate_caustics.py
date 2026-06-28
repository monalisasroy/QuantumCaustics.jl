"""
animate_caustics.py
===================
Generate caustics_animation.gif and caustics_spacetime.png for the
QuantumCaustics.jl README landing page.

Physics
-------
Exact single-excitation formula for the quenched TFIM (Eqs. 9-10 of
Singh Roy et al., PRA 2026, arXiv:2410.06803).

    ΔZ_j(t) = 2 |1/N Σ_k exp(i [k(j−j₀) + 2 J^xx cos(k) t])|²

ΔZ_j = <Z_j>^{no-flip} − <Z_j>^{flip}.  The initial spin flip at j₀
gives ΔZ_{j₀}(0) = 2 exactly; values lie in [0, 2].

Colour scale
------------
magma colourmap, log scale (LogNorm, vmin = 0.005, vmax = 2).
The log stretch makes the caustic fringes (amplitude ~0.1–0.25 at late
times) legible against the full range up to 2.  Zero values (outside
the light cone) are masked and appear as the dark background.

Parameters reproduce Fig. 2 of the paper: N = 79, J^xx = 0.4, h^z = 1.0.

Output
------
    caustics_animation.gif      animated GIF for the README
    caustics_spacetime.png      static spacetime colourmap
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm
import matplotlib.ticker as mticker

# ── Output directory ─────────────────────────────────────────────────────────
# Change to "docs" when running from the repository root.
OUT_DIR = "."

# ── Physics ──────────────────────────────────────────────────────────────────
N   = 79
j0  = N // 2        # centre spin, 0-indexed (= 39)
Jxx = 0.4           # paramagnetic: J^xx < h^z = 1

k_grid = 2.0 * np.pi * np.arange(N) / N   # N-point DFT k-grid
j_grid = np.arange(N, dtype=float)

def delta_Z(t: float) -> np.ndarray:
    """ΔZ_j(t) from Eqs. 9-10 of the paper."""
    phase = k_grid * (j_grid[:, None] - j0) + 2.0 * Jxx * np.cos(k_grid) * t
    psi   = np.sum(np.exp(1j * phase), axis=1) / N
    return 2.0 * np.abs(psi) ** 2

# ── Precompute spacetime data ────────────────────────────────────────────────
T_MAX  = 33.0
N_T    = 400
t_all  = np.linspace(0.02, T_MAX, N_T)

print("Precomputing ΔZ_j(t) …")
dZ_all = np.stack([delta_Z(t) for t in t_all])   # (N_T, N)
print(f"  shape {dZ_all.shape},  max ΔZ = {dZ_all.max():.4f}  (should be ≈ 2)")

# ── Colour scale ─────────────────────────────────────────────────────────────
# magma colourmap; log scale stretches the low fringe amplitudes
# (~0.1–0.25 at late times) across the full colour range up to vmax = 2.
# Zero values (outside the light cone) are masked by LogNorm.
CMAP  = "magma"
NORM  = LogNorm(vmin=0.005, vmax=2.0)
CBAR_LABEL = r"$\Delta Z_j$  (log scale)"
CBAR_TICKS = [0.01, 0.1, 1.0, 2.0]   # decimal labels, log-spaced

# ── Figure background ─────────────────────────────────────────────────────────
BG = "#0e0e0e"

# ── Static spacetime figure ──────────────────────────────────────────────────
fig_s, ax_s = plt.subplots(figsize=(7.0, 5.0), facecolor=BG)
ax_s.set_facecolor(BG)
pm_s = ax_s.pcolormesh(
    t_all, j_grid, dZ_all.T,
    cmap=CMAP, norm=NORM, shading="nearest",
)
cbar_s = fig_s.colorbar(pm_s, ax=ax_s, pad=0.01, fraction=0.025)
cbar_s.set_ticks(CBAR_TICKS)
cbar_s.ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
cbar_s.ax.tick_params(colors="#888888")
cbar_s.set_label(CBAR_LABEL, color="#cccccc", fontsize=8)
ax_s.set_xlim(0, T_MAX)
ax_s.set_ylim(-0.5, N - 0.5)
ax_s.set_xlabel(r"Time, $t$", color="#cccccc", fontsize=12)
ax_s.set_ylabel(r"Position, $j$", color="#cccccc", fontsize=12)
ax_s.tick_params(colors="#888888")
for sp in ax_s.spines.values():
    sp.set_color("#555555")
ax_s.set_title(
    r"$\Delta Z_j(t)$   —   TFIM,  $J^{xx} = 0.4$,  $N = 79$",
    color="white", fontsize=12, pad=6,
)
ax_s.text(
    0.98, 0.03,
    "Singh Roy et al., PRA 2026  ·  arXiv:2410.06803",
    transform=ax_s.transAxes, color="#666666", fontsize=7.5,
    ha="right", va="bottom",
)
fig_s.tight_layout()
static_path = os.path.join(OUT_DIR, "caustics_spacetime.png")
fig_s.savefig(static_path, dpi=160, facecolor=BG)
plt.close(fig_s)
print(f"Saved static: {static_path}")

# ── Animated GIF ─────────────────────────────────────────────────────────────
N_FRAMES  = 100    # main sweep frames
HOLD_END  = 18     # extra frames pausing on the complete pattern
FPS       = 12

fi_main = np.round(np.linspace(0, N_T - 1, N_FRAMES)).astype(int)
fi_end  = np.full(HOLD_END, N_T - 1, dtype=int)
fi_seq  = np.concatenate([fi_main, fi_end])
TOTAL   = len(fi_seq)

fig = plt.figure(figsize=(8.6, 5.0), facecolor=BG)
gs  = GridSpec(
    2, 1, figure=fig,
    height_ratios=[4, 1],
    hspace=0.05,
    left=0.075, right=0.925, top=0.91, bottom=0.09,
)
ax_st = fig.add_subplot(gs[0])   # spacetime panel
ax_1d = fig.add_subplot(gs[1])   # 1-D snapshot strip

for ax in (ax_st, ax_1d):
    ax.set_facecolor(BG)

# ── Spacetime panel: pre-render everything, reveal with a moving cover ──────
pm = ax_st.pcolormesh(
    t_all, j_grid, dZ_all.T,
    cmap=CMAP, norm=NORM, shading="nearest",
)

cover = Rectangle(
    (0.0, -0.5), T_MAX + 1.0, N + 1.0,
    facecolor=BG, edgecolor="none", zorder=3,
)
ax_st.add_patch(cover)

ax_st.set_xlim(0, T_MAX)
ax_st.set_ylim(-0.5, N - 0.5)
ax_st.set_ylabel(r"Position, $j$", color="#cccccc", fontsize=10)
ax_st.tick_params(axis="both", colors="#666666", labelbottom=False)
for sp in ax_st.spines.values():
    sp.set_color("#444444")

ax_st.set_title(
    r"QuantumCaustics.jl   ·   "
    r"$\Delta Z_j(t) = \langle Z_j\rangle^{\rm no\text{-}flip}"
    r" - \langle Z_j\rangle^{\rm flip}$",
    color="white", fontsize=10.5, pad=6,
)

t_txt = ax_st.text(
    0.985, 0.955, "",
    transform=ax_st.transAxes,
    color="#ffcc88", fontsize=10, ha="right", va="top", zorder=5,
)

ax_st.text(
    0.015, 0.03,
    r"$J^{xx} = 0.4$,  $h^z = 1$,  $N = 79$",
    transform=ax_st.transAxes,
    color="#888888", fontsize=8.5, ha="left", va="bottom", zorder=5,
)

cbar = fig.colorbar(pm, ax=ax_st, pad=0.01, fraction=0.024)
cbar.set_ticks(CBAR_TICKS)
cbar.ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
cbar.ax.tick_params(colors="#666666")
cbar.set_label(CBAR_LABEL, color="#bbbbbb", fontsize=8.5)

# ── 1-D snapshot strip ───────────────────────────────────────────────────────
snap = ax_1d.imshow(
    np.zeros((1, N)),
    aspect="auto", cmap=CMAP, norm=NORM,
    origin="lower", extent=[-0.5, N - 0.5, 0, 1],
)
ax_1d.axvline(x=j0, color="white", lw=0.8, alpha=0.35, zorder=3)
ax_1d.set_xlim(-0.5, N - 0.5)
ax_1d.set_ylim(0, 1)
ax_1d.set_xlabel(r"Site, $j$", color="#cccccc", fontsize=10)
ax_1d.tick_params(axis="both", colors="#666666", left=False, labelleft=False)
for sp in ax_1d.spines.values():
    sp.set_color("#444444")

# ── Animation ────────────────────────────────────────────────────────────────
def init():
    cover.set_x(0.0)
    cover.set_width(T_MAX + 1.0)
    snap.set_data(np.zeros((1, N)))
    t_txt.set_text("")
    return []

def update(frame: int):
    ti    = fi_seq[frame]
    t_now = t_all[ti]
    cover.set_x(t_now)
    cover.set_width(max(0.0, T_MAX - t_now + 0.6))
    snap.set_data(dZ_all[ti : ti + 1, :])
    t_txt.set_text(rf"$t = {t_now:.1f}$")
    return []

ani = animation.FuncAnimation(
    fig, update, frames=TOTAL, init_func=init,
    interval=1000 // FPS, blit=False,
)

gif_path = os.path.join(OUT_DIR, "caustics_animation.gif")
print(f"Rendering {TOTAL} frames @ {FPS} fps …")
ani.save(gif_path, writer="pillow", fps=FPS, dpi=95)
plt.close(fig)
print(f"Saved animated GIF: {gif_path}")
