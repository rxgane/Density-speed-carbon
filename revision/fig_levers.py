"""Figure 2: annual cost and carbon of the reference network as each decision moves on its own from the
cost-minimizing design (the other held there), in percent of that design. Reads levers_data.json."""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 9, "axes.linewidth": 0.6})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e2e1dd"
D = json.load(open(os.path.join(here, "levers_data.json")))
K0, C0 = D["costmin"]["K"], D["costmin"]["C"]
def pct(rows, key): return np.array([100*(r[key]/(K0 if key == "K" else C0) - 1) for r in rows])
NMAX = 15.0
A = [r for r in sorted(D["vary_n_fixed_v"], key=lambda r: r["n"]) if r["n"] <= NMAX + 1e-9]
assert abs(A[-1]["n"] - NMAX) < 1e-9
B = sorted(D["vary_v_fixed_n"], key=lambda r: r["days"])
xa, xb = np.array([r["n"] for r in A]), np.array([r["days"] for r in B])
cm, cb = D["costmin"], D["carbonmin_speed_at_costmin_n"]; ca = A[-1]   # (a) ends at 15 centers; (b) at intermodal rail
COST = dict(color=INK, lw=1.7, zorder=3); CARB = dict(color=INK, lw=1.7, ls=(0, (4.5, 2)), zorder=3)
L = dict(fontsize=8, color=INK); M = dict(fontsize=7.5, color=MUTED)

fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.85), sharey=True, gridspec_kw=dict(wspace=0.10))
for ax in axs:
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(MUTED); ax.spines["bottom"].set_color(MUTED)
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)

def panel(ax, x, rows, xmin, xstar, smin):
    k, c = pct(rows, "K"), pct(rows, "C")
    ax.plot(x, k, **COST); ax.plot(x, c, **CARB)
    ax.plot(xmin, 0, marker="o", ms=5.5, mfc="white", mec=INK, mew=1.0, zorder=6)
    ks, cs = 100*(smin["K"]/K0 - 1), 100*(smin["C"]/C0 - 1)
    ax.plot([xstar, xstar], [cs, ks], color=MUTED, lw=0.7, ls=(0, (1, 1.5)), zorder=2)
    for y in (cs, ks): ax.plot(xstar, y, marker="o", ms=4.3, color=INK, mec="white", mew=0.8, zorder=5)
    return k, c, ks, cs

# (a) network: n moves, line haul held at the cost-minimizing speed
ax = axs[0]
k, c, ks, cs = panel(ax, xa, A, cm["n"], ca["n"], ca)
ax.text(ca["n"] - 0.25, ks + 1.5, rf"${ks:+.1f}\%$", ha="right", va="bottom", **L)
ax.text(ca["n"] - 0.2, cs + 1.6, rf"${cs:+.1f}\%$", ha="right", va="bottom", **L)
ax.text(11.0, float(np.interp(11.0, xa, k)) + 2.5, "Cost", ha="right", va="bottom", fontsize=8.5, color=INK)
ax.text(11.0, float(np.interp(11.0, xa, c)) + 2.0, "Carbon", ha="center", va="bottom", fontsize=8.5, color=INK)
ax.annotate("cost minimum", xy=(cm["n"], 0), xytext=(5.2, 20.0), ha="left", va="center", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=1, shrinkB=4), **M)
ax.set_xlim(1, NMAX); ax.set_xticks([1, 5, 10, 15])
ax.set_xlabel("Distribution centers", labelpad=12)
ax.set_ylabel(r"Change from the cost-minimizing" "\n" r"design (\%)")

# (b) replenishment: speed moves, n held at the cost-minimizing number of centers
ax = axs[1]
k, c, ks, cs = panel(ax, xb, B, cm["days"], cb["days"], cb)
assert abs(cb["days"] - 6) < 1e-6
ax.text(cb["days"] - 0.08, ks + 1.5, rf"${ks:+.1f}\%$", ha="right", va="bottom", **L)
ax.text(cb["days"] - 0.06, cs - 2.5, rf"${cs:.1f}\%$", ha="right", va="top", **L)
ax.text(4.0, float(np.interp(4.0, xb, k)) + 3.0, "Cost", ha="center", va="bottom", fontsize=8.5, color=INK)
ax.text(3.3, float(np.interp(3.3, xb, c)) - 3.5, "Carbon", ha="right", va="top", fontsize=8.5, color=INK)
ax.annotate("cost minimum", xy=(cm["days"], 0), xytext=(2.45, 20.0), ha="left", va="center", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=1, shrinkB=4), **M)
ax.set_xlim(2.0, 6.0); ax.set_xticks([2, 3, 4, 5, 6])
ax.set_xlabel("Line-haul transit (days)", labelpad=12)
ax.text(2.0, -0.115, "truckload", transform=ax.get_xaxis_transform(), ha="left", va="top", **M)
ax.text(4.0, -0.115, "intermodal", transform=ax.get_xaxis_transform(), ha="center", va="top", **M)
ax.text(6.0, -0.115, "carload", transform=ax.get_xaxis_transform(), ha="right", va="top", **M)

axs[0].set_ylim(-53, 36); axs[0].set_yticks([-40, -20, 0, 20]); axs[0].set_yticklabels(["$-40$", "$-20$", "$0$", "$+20$"])
for ax, t, s in ((axs[0], "(a) Network", "line haul at 2.15 days"), (axs[1], "(b) Replenishment", "3.9 centers")):
    ax.text(0.0, 1.03, t, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.5, color=INK)
    ax.text(1.0, 1.03, s, transform=ax.transAxes, ha="right", va="bottom", **M)
fig.savefig(os.path.join(here, "..", "fig_levers.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(here, "fig_levers.png"), dpi=200, bbox_inches="tight")
print("ok", round(cm["n"], 2), round(cm["days"], 3), round(K0, 2), round(C0))
