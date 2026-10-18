"""Figure 3: increase in annual logistics cost against the required cut in emissions, for replenishment alone, the
network alone, and both decisions (reference network); each curve ends at the exact largest cut of its action set.
Reads targets_v22.json and levers_data.json."""
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
KC = json.load(open(os.path.join(here, "levers_data.json")))
END = {"replen": KC["carbonmin_speed_at_costmin_n"], "network": KC["carbonmin_n_at_costmin_speed"], "both": KC["carbonmin"]}
assert abs(KC["costmin"]["K"]*1e6 - ref["K"]) < 1 and abs(KC["costmin"]["C"]*1e3 - ref["C"]) < 1
x0 = np.array([r["x"] for r in rows])*100
def series(mode, key): return np.array([np.nan if r[mode] is None else r[mode][key] for r in rows], dtype=float)
fig, ax = plt.subplots(figsize=(6.0, 3.0))
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.spines["left"].set_color(MUTED); ax.spines["bottom"].set_color(MUTED)
ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
xmax = 0
for mode in ("network", "both", "replen"):
    K = series(mode, "K"); ok = ~np.isnan(K)
    e = END[mode]; xe = 100*(1 - e["C"]*1e3/ref["C"]); assert xe > x0[ok].max()
    x = np.append(x0[ok], xe); y = np.append(100*(K[ok]/ref["K"]-1), 100*(e["K"]*1e6/ref["K"]-1)); xmax = max(xmax, xe)
    ax.plot(x, y, color=COL[mode], ls=LS[mode], lw=2.2 if mode == "replen" else 1.8, label=LAB[mode], zorder=4 if mode == "replen" else 3)
    ax.plot(x[-1], y[-1], marker="o", ms=4.5, color=COL[mode], mec="white", mew=0.8, zorder=5)
ax.set_xlim(0, 5*np.ceil(xmax/5)); ax.set_ylim(bottom=0)
ax.set_xlabel(r"Required reduction in annual emissions from the cost-minimizing design (\%)")
ax.set_ylabel("Increase in annual\n" r"logistics cost (\%)")
h, l = ax.get_legend_handles_labels(); order = [l.index(k) for k in ("Replenishment only", "Network only", "Both decisions")]
ax.legend([h[i] for i in order], [l[i] for i in order], loc="upper left", frameon=False, fontsize=8, handlelength=2.6)
fig.savefig(os.path.join(here, "..", "fig_target.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(here, "fig_target.png"), dpi=200, bbox_inches="tight")
print("ok")
