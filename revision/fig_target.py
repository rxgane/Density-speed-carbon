"""Figure 3: least-cost design for each required cut with replenishment, the network, or both (reference network
that finances its stock in transit); each curve ends at the exact largest cut of its action set. Reads
targets_v22.json and levers_data.json."""
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
KC = json.load(open(os.path.join(here, "levers_data.json")))
END = {"replen": KC["carbonmin_speed_at_costmin_n"], "network": KC["carbonmin_n_at_costmin_speed"], "both": KC["carbonmin"]}
assert abs(KC["costmin"]["K"]*1e6 - ref["K"]) < 1 and abs(KC["costmin"]["C"]*1e3 - ref["C"]) < 1
x0 = x
for mode in ("network", "both", "replen"):
    K = series(mode, "K"); n = series(mode, "n"); d = series(mode, "days"); ok = ~np.isnan(K)
    e = END[mode]; xe = 100*(1 - e["C"]*1e3/ref["C"]); assert xe > x0[ok].max()   # exact largest cut of this action set
    x = np.append(x0[ok], xe); K = np.append(K[ok], e["K"]*1e6); n = np.append(n[ok], e["n"]); d = np.append(d[ok], e["days"])
    kw = dict(color=COL[mode], ls=LS[mode], lw=2.2 if mode == "replen" else 1.8, label=LAB[mode], zorder=4 if mode == "replen" else 3)
    axs[0].plot(x, 100*(K/ref["K"]-1), **kw); axs[1].plot(x, n, **kw); axs[2].plot(x, d, **kw)
    last = np.where(~np.isnan(K))[0].max()
    for ax, y in zip(axs, (100*(K[last]/ref["K"]-1), n[last], d[last])):
        ax.plot(x[last], y, marker="o", ms=4.5, color=COL[mode], mec="white", mew=0.8, zorder=5)
axs[0].set_ylabel("Increase in annual\n" r"logistics cost (\%)"); axs[1].set_ylabel("Distribution\ncenters")
axs[2].set_ylabel("Line-haul\ntransit (days)")
axs[2].set_xlabel(r"Required reduction in annual emissions from the cost-minimizing design (\%)")
axs[2].set_xlim(0, 50); axs[0].set_ylim(bottom=0); axs[1].set_ylim(0, 8); axs[1].set_yticks([0, 2, 4, 6, 8]); axs[2].set_ylim(2.0, 6.15)
for ax, t in zip(axs, ("(a) Cost", "(b) Network", "(c) Replenishment")):
    ax.text(0.0, 1.02, t, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.5, color=INK)
h, l = axs[0].get_legend_handles_labels(); order = [l.index(k) for k in ("Replenishment only", "Network only", "Both decisions")]
axs[0].legend([h[i] for i in order], [l[i] for i in order], loc="upper left", bbox_to_anchor=(0.36, 1.0), frameon=False, fontsize=8, handlelength=2.6)
fig.savefig(os.path.join(here, "..", "fig_target.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(here, "fig_target.png"), dpi=200, bbox_inches="tight")
print("ok")
