"""Figure: meeting a carbon target with replenishment, the network, or both (reference network that finances its
stock in transit). Reads targets_v22.json."""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 9, "axes.linewidth": 0.6})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e2e1dd"
COL = {"both": "#2a78d6", "network": "#eb6834", "replen": "#1baf7a"}
LAB = {"both": "Both decisions", "network": "Network only", "replen": "Replenishment only"}
LS = {"both": "-", "network": (0, (5, 2)), "replen": (0, (1.5, 1.5))}
R = json.load(open(os.path.join(here, "targets_v22.json"))); ref = R["ref"]; rows = R["rows"]
x = np.array([r["x"] for r in rows])*100
def series(mode, key): return np.array([np.nan if r[mode] is None else r[mode][key] for r in rows], dtype=float)
fig, axs = plt.subplots(3, 1, figsize=(6.0, 6.4), sharex=True, gridspec_kw=dict(hspace=0.28))
for ax in axs:
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(MUTED); ax.spines["bottom"].set_color(MUTED)
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
for mode in ("network", "both", "replen"):
    K = series(mode, "K"); n = series(mode, "n"); d = series(mode, "days")
    kw = dict(color=COL[mode], ls=LS[mode], lw=2.2 if mode == "replen" else 1.8, label=LAB[mode], zorder=4 if mode == "replen" else 3)
    axs[0].plot(x, 100*(K/ref["K"]-1), **kw); axs[1].plot(x, n, drawstyle="steps-post", **kw); axs[2].plot(x, d, **kw)
    last = np.where(~np.isnan(K))[0].max()
    for ax, y in zip(axs, (100*(K[last]/ref["K"]-1), n[last], d[last])):
        ax.plot(x[last], y, marker="o", ms=4.5, color=COL[mode], mec="white", mew=0.8, zorder=5)
Kc = np.array([np.nan if r["both_cont"] is None else r["both_cont"]["K"] for r in rows], dtype=float)
nc = np.array([np.nan if r["both_cont"] is None else r["both_cont"]["n"] for r in rows], dtype=float)
dc = np.array([np.nan if r["both_cont"] is None else r["both_cont"]["days"] for r in rows], dtype=float)
axs[0].plot(x, 100*(Kc/ref["K"]-1), color=MUTED, lw=0.8, ls=(0, (3, 2)), label=r"Both, continuous number of centers")
axs[1].plot(x, nc, color=MUTED, lw=0.8, ls=(0, (3, 2))); axs[2].plot(x, dc, color=MUTED, lw=0.8, ls=(0, (3, 2)))
axs[0].set_ylabel("Increase in annual\n" r"logistics cost (\%)"); axs[1].set_ylabel("Distribution\ncenters")
axs[2].set_ylabel("Line-haul\ntransit (days)")
axs[2].set_xlabel(r"Required reduction in annual emissions from the cost-minimizing design (\%)")
axs[2].set_xlim(0, 50); axs[0].set_ylim(bottom=0); axs[1].set_ylim(0, max(12, np.nanmax(series("both", "n"))+1)); axs[2].set_ylim(2.0, 4.1)
for ax, t in zip(axs, ("(a) Cost", "(b) Network", "(c) Replenishment")):
    ax.text(0.0, 1.02, t, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.5, color=INK)
h, l = axs[0].get_legend_handles_labels(); order = [l.index(k) for k in ("Replenishment only", "Network only", "Both decisions", r"Both, continuous number of centers")]
axs[0].legend([h[i] for i in order], [l[i] for i in order], loc="upper left", bbox_to_anchor=(0.36, 1.0), frameon=False, fontsize=8, handlelength=2.6)
fig.savefig(os.path.join(here, "..", "fig_target.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(here, "fig_target.png"), dpi=200, bbox_inches="tight")
print("ok")
