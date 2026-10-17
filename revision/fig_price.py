"""Figure: how a carbon price changes emissions, the number of distribution centers, and line-haul time
in the reference network, for three action sets. Reads application_v22_price.json."""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
here = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 9, "axes.linewidth": 0.6})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e2e1dd"
COL = {"both": "#2a78d6", "network": "#eb6834", "replen": "#1baf7a"}
LAB = {"both": "Both decisions", "network": "Network only", "replen": "Replenishment only"}
LS = {"both": "-", "network": (0, (5, 2)), "replen": (0, (1.5, 1.5))}
rows = json.load(open(os.path.join(here, "application_v22_price.json")))
x = np.array([r["price"] for r in rows])
fig, axs = plt.subplots(3, 1, figsize=(6.0, 6.2), sharex=True, gridspec_kw=dict(hspace=0.3))
for ax in axs:
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(MUTED); ax.spines["bottom"].set_color(MUTED)
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.axvline(100, color=GRID, lw=1.0, zorder=0)
for mode in ("network", "both", "replen"):
    y = [np.array([r[mode][k] for r in rows]) for k in ("cut", "n", "days")]
    kw = dict(color=COL[mode], ls=LS[mode], lw=2.2 if mode == "replen" else 1.8, label=LAB[mode], zorder=4 if mode == "replen" else 3)
    axs[0].plot(x, y[0], **kw); axs[1].plot(x, y[1], drawstyle="steps-post", **kw); axs[2].plot(x, y[2], **kw)
axs[2].set_xscale("log"); axs[2].set_xlim(10, 5000)
ticks = [10, 30, 100, 300, 1000, 3000]
axs[2].xaxis.set_major_locator(FixedLocator(ticks)); axs[2].xaxis.set_minor_locator(NullLocator())
axs[2].xaxis.set_major_formatter(FixedFormatter([f"{t:,}" for t in ticks]))
axs[2].set_xlabel(r"Carbon price (US\$ per tonne CO$_2$e)")
axs[0].set_ylabel("Reduction in annual\n" r"emissions (\%)"); axs[1].set_ylabel("Distribution\ncenters")
axs[2].set_ylabel("Line-haul\ntransit (days)")
axs[0].set_ylim(0, 50); axs[1].set_ylim(0, 10); axs[1].set_yticks(range(0, 11, 2)); axs[2].set_ylim(2.0, 4.1)
for ax, t in zip(axs, ("(a) Emissions", "(b) Network", "(c) Replenishment")):
    ax.text(0.0, 1.02, t, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.5, color=INK)
axs[0].text(100, 48, r"\$100", ha="center", va="top", fontsize=7.5, color=MUTED, bbox=dict(fc="white", ec="none", pad=0.5))
h, l = axs[0].get_legend_handles_labels(); order = [l.index(k) for k in ("Replenishment only", "Network only", "Both decisions")]
axs[0].legend([h[i] for i in order], [l[i] for i in order], loc="upper left", frameon=False, fontsize=8, handlelength=2.6)
fig.savefig(os.path.join(here, "..", "fig_price.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(here, "fig_price.png"), dpi=200, bbox_inches="tight")
print("ok")
