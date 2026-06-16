"""Figures for the numerical illustration. Reuses the model in numerics.py."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from numerics import (Params, minimize_obj, minimize_carbon, build_archetypes)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

p = Params()

# --------------------------------------------------------------------------
# Figure 1: the four-regime plane, with archetypes located by their ratios
# --------------------------------------------------------------------------

from matplotlib.patches import Rectangle
fig, ax = plt.subplots(figsize=(5.8, 4.6))
lo, hi = 0.2, 8.0
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)

# highlight the novel regime (top-right: faster AND more dispersed)
ax.add_patch(Rectangle((1, 1), hi - 1, hi - 1, color="0.92", zorder=0))
ax.axvline(1, color="0.4", lw=1.0, ls="--")
ax.axhline(1, color="0.4", lw=1.0, ls="--")

# quadrant labels in the corners (direction relative to cost optimum)
ax.text(0.26, 0.24, "Slow &\nconsolidate\n(conventional)",
        ha="center", va="center", fontsize=8.5, color="0.35")
ax.text(5.6, 0.24, "Slow,\nstock near", ha="center", va="center",
        fontsize=8.5, color="0.35")
ax.text(0.26, 5.6, "Faster,\npool upstream", ha="center", va="center",
        fontsize=8.5, color="0.35")
ax.text(2.0, 6.2, "Faster & more dispersed",
        ha="center", va="center", fontsize=9, fontweight="bold", color="0.0")

# archetypes (R_density on x, R_speed on y), labels placed with leader lines
arch = build_archetypes(p)
pts = {
    "Cold-chain pharma":   (arch["Cold-chain pharma (refrig., diesel last mile)"], (1.75, 3.3)),
    "E-grocery perishable":(arch["E-grocery perishables (refrig., van LM)"],       (6.4, 2.7)),
    "Bulky low-value":     (arch["Bulky low-value (cheap storage, diesel LM)"],    (4.33, 0.70)),
    "Ambient dry goods":   (arch["Ambient dry goods (stable, rail/efficient LM)"], (0.40, 0.72)),
}
for name, (pp, xytext) in pts.items():
    x, y = pp.R_density, pp.R_speed
    ax.scatter([x], [y], s=44, color="black", zorder=5)
    ax.annotate(name, (x, y), xytext, fontsize=8.5, ha="center", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.5",
                                shrinkA=2, shrinkB=4))

# baseline marker
ax.scatter([p.R_density], [p.R_speed], s=44, marker="s",
           facecolor="white", edgecolor="black", zorder=5)
ax.annotate("baseline", (p.R_density, p.R_speed),
            (1.25, 0.32), fontsize=8.5, ha="center",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="0.5",
                            shrinkA=2, shrinkB=4))

ax.set_xlabel(r"$R_{\mathrm{density}} = (\varepsilon_L/t_L)\,/\,(\varphi/F)$")
ax.set_ylabel(r"$R_{\mathrm{speed}} = (\eta/h)\,/\,(\varepsilon_0/t_1)$")
ax.set_xticks([0.3, 0.5, 1, 2, 4, 7]); ax.set_yticks([0.3, 0.5, 1, 2, 4, 7])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
fig.savefig("fig_regimes.pdf")
print("wrote fig_regimes.pdf")

# --------------------------------------------------------------------------
# Figure 2: cost-carbon Pareto frontier (baseline), swept over tau
# --------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(5.6, 4.0))
taus = np.concatenate([[0.0], np.geomspace(1, 5000, 60), [1e9]])
Ks, Cs = [], []
for tau in taus:
    _, _, K, C = minimize_obj(p, tau)
    Ks.append(K / 1e6); Cs.append(C / 1e3)
Ks, Cs = np.array(Ks), np.array(Cs)

ax.margins(x=0.12, y=0.10)
ax.plot(Ks, Cs, "-", color="0.2", lw=1.4)
ax.scatter([Ks[0]], [Cs[0]], s=45, color="black", zorder=5)
ax.annotate(r"cost-min ($\tau=0$)", (Ks[0], Cs[0]),
            (Ks[0] + 0.0010, Cs[0] - 0.2), fontsize=9, ha="left", va="center")
ax.scatter([Ks[-1]], [Cs[-1]], s=45, marker="D", color="black", zorder=5)
ax.annotate(r"carbon-min ($\tau\to\infty$)", (Ks[-1], Cs[-1]),
            (Ks[-1] - 0.0008, Cs[-1] + 1.1), fontsize=9, ha="right", va="center")

ax.set_xlabel(r"cost $K$ (\$ millions / week)")
ax.set_ylabel(r"carbon $C$ (tonnes / week)")
ax.grid(True, color="0.9", lw=0.6)
fig.savefig("fig_frontier.pdf")
print("wrote fig_frontier.pdf")
