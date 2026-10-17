"""Figure 4: number of distribution centers against the carbon intensity of faster inbound service, r_t, with
delivery by diesel truck (x = 1.375 kg CO2e per loaded vehicle-mile). (a) carbon minimum, perishables that spoil;
(b) least-cost design for a 25% cut, two versions of the reference network. Curves from inbound_fine.json; labeled
numbers from inbound_data.json, checked against the fine grid."""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); FIG3 = HERE
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 9, "axes.linewidth": 0.6})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e2e1dd"
BLUE, ORANGE = "#2a78d6", "#eb6834"
X = 1.375
D = json.load(open(os.path.join(FIG3, "inbound_data.json"))); F = json.load(open(os.path.join(HERE, "inbound_fine.json")))

def rec(prod, pool, rt):
    r = [p for p in D["points"] if p["product"] == prod and p["pooling"] == pool and abs(p["x"]-X) < 1e-9 and abs(p["rt"]-rt) < 1e-9]
    assert r, (prod, pool, rt); return r[0]

def curve(prod, pool, crit):
    pts = sorted([p for p in F["points"] if p["product"] == prod and p["pooling"] == pool], key=lambda p: p["rt"])
    if crit == "cut": pts = [p for p in pts if p["cut"] is not None]
    return np.array([p["rt"] for p in pts]), np.array([p[crit]["n"] for p in pts])

def check(prod, pool, rt, crit):
    d = rec(prod, pool, rt)[crit]["n"]
    fr = [p for p in F["points"] if p["product"] == prod and p["pooling"] == pool and abs(p["rt"]-rt) < 1e-9]
    if fr: assert abs(fr[0][crit]["n"] - d) < 1e-9
    return d

OM = {(p, l): rec(p, l, 1.78)["omega"] for p in ("perishables", "reference") for l in ("lower", "higher")}
STY = {"higher": dict(color=BLUE, ls="-", lw=1.8, zorder=3), "lower": dict(color=ORANGE, ls=(0, (5, 2)), lw=1.8, zorder=2)}
T = dict(fontsize=7.8, color=INK, bbox=dict(fc="white", ec="none", pad=0.4)); M = dict(fontsize=7.4, color=MUTED)

fig, axs = plt.subplots(2, 1, figsize=(6.4, 4.6), sharex=True, gridspec_kw=dict(hspace=0.32))
for ax in axs:
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(MUTED); ax.spines["bottom"].set_color(MUTED)
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_xscale("log"); ax.set_xlim(0.01, 3.0); ax.set_ylabel("Distribution centers")

# (a) carbon minimum, perishables that spoil; their own service is team drivers
ax = axs[0]
for pool in ("higher", "lower"):
    r, n = curve("perishables", pool, "carbonmin"); ax.plot(r, n, **STY[pool])
hi03, hi178 = check("perishables", "higher", 0.03, "carbonmin"), check("perishables", "higher", 1.78, "carbonmin")
dfast = min(p["carbonmin"]["days"] for p in F["points"] if p["product"] == "perishables" and p["pooling"] == "higher")
hifast = max(p["carbonmin"]["n"] for p in F["points"] if p["product"] == "perishables" and p["pooling"] == "higher")
assert abs(dfast-1) < 1e-6 and abs(rec("perishables", "higher", 1.78)["carbonmin"]["days"]-2) < 1e-6
lo = [check("perishables", "lower", rt, "carbonmin") for rt in (0.03, 1.78)]; LO = round(np.mean(lo)); assert all(abs(v-LO) < 0.3 for v in lo)
ax.set_ylim(0, 1.3*hifast); ax.set_yticks([0, 5, 10, 15])
ax.axvline(0.03, color=GRID, lw=0.9, zorder=0); ax.text(0.03*1.07, 0.09*1.3*hifast, "team drivers ($r_t=0.03$)", ha="left", va="center", **M)
ax.text(0.0108, hifast+0.5, rf"{hifast:.1f} centers; line haul at 1 day, the fastest allowed", ha="left", va="bottom", **T)
ax.text(2.9, hi178-0.6, rf"{hi178:.1f} centers; line haul at 2 days, the slowest allowed", ha="right", va="top", **T)
ax.text(2.9, LO+0.5, rf"about {LO:.0f} centers at every $r_t$", ha="right", va="bottom", **T)
ax.text(0.0108, hifast-0.6, r"higher pooling weight ($\omega=%.2f$)" % OM[("perishables", "higher")], ha="left", va="top", fontsize=7.8, color=BLUE)
ax.text(0.0108, LO-0.6, r"lower pooling weight ($\omega=%.2f$)" % OM[("perishables", "lower")], ha="left", va="top", fontsize=7.8, color=ORANGE)
ax.set_title(r"(a) Carbon minimum, perishables that spoil", loc="left", fontsize=9.5, color=INK, pad=6)

# (b) least-cost design for a 25% cut, reference network; its own service choice is intermodal or truckload
ax = axs[1]
for pool in ("higher", "lower"):
    r, n = curve("reference", pool, "cut"); ax.plot(r, n, **STY[pool])
    ax.plot(r[0], n[0], marker="o", ms=5, mfc="white", mec=STY[pool]["color"], mew=1.3, zorder=6)
rH, nH = curve("reference", "higher", "cut"); rL, nL = curve("reference", "lower", "cut")
ymax = 1.25*max(nH.max(), nL.max()); ax.set_ylim(0, ymax); ax.set_yticks([0, 2, 4, 6, 8])
ax.text(rH[0]/1.10, nH[0], "25\\% cut cannot be met\nat lower $r_t$", ha="right", va="center", linespacing=1.0, **M)
pts = sorted([p for p in F["points"] if p["product"] == "reference" and p["pooling"] == "higher" and p["cut"]], key=lambda p: abs(p["rt"]-1.0))
c1 = pts[0]; r1 = c1["rt"]; c1 = c1["cut"]
ax.plot([r1], [c1["n"]], marker="o", ms=3.4, color=BLUE, zorder=5)
ax.annotate(r"%.1f centers; slower line haul" "\n" r"delivers %.0f\%% of the cut" % (c1["n"], 100*c1["share_inbound"]),
            xy=(r1, c1["n"]-0.05), xytext=(0.45, 3.9), ha="right", va="top", linespacing=1.0,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=1, shrinkB=2), **T)
cH, cL = rec("reference", "higher", 1.78)["cut"], rec("reference", "lower", 1.78)["cut"]
sh = [100*cH["share_inbound"], 100*cL["share_inbound"]]
ax.plot([1.78], [cH["n"]], marker="o", ms=3.4, color=BLUE, zorder=5); ax.plot([1.78], [cL["n"]], marker="o", ms=3.4, color=ORANGE, zorder=5)
ax.annotate(r"%.1f and %.1f centers; slower line haul" "\n" r"delivers %.0f\%% and %.0f\%% of the cut" % (cH["n"], cL["n"], sh[0], sh[1]),
            xy=(1.78, cH["n"]-0.12), xytext=(2.95, 2.95), ha="right", va="top", linespacing=1.0,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=1, shrinkB=2), **T)
ax.annotate("", xy=(1.78, cL["n"]-0.12), xytext=(2.2, 2.95), arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=0, shrinkB=2))
ax.axvline(1.78, color=GRID, lw=0.9, zorder=0); ax.text(1.78/1.07, 0.07*ymax, "intermodal to truckload ($r_t=1.78$)", ha="right", va="center", **M)
ax.text(rH[0]/1.10, nH[0]-1.25, r"higher pooling weight ($\omega=%.2f$)" % OM[("reference", "higher")], ha="right", va="center", fontsize=7.8, color=BLUE)
ax.text(rL[0]*1.12, nL[0]+0.2, r"lower pooling weight ($\omega=%.2f$)" % OM[("reference", "lower")], ha="left", va="center", fontsize=7.8, color=ORANGE)
ax.set_title(r"(b) Least-cost design for a 25\% cut, two versions of the reference network", loc="left", fontsize=9.5, color=INK, pad=6)
ax.set_xticks([0.01, 0.03, 0.1, 0.3, 1, 3]); ax.set_xticklabels(["0.01", "0.03", "0.1", "0.3", "1", "3"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel(r"Carbon intensity of faster inbound service, $r_t$ (kg CO$_2$e per dollar spent on speed)")

fig.savefig(os.path.join(HERE, "..", "fig_inbound.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(HERE, "fig_inbound.png"), dpi=200, bbox_inches="tight")
print("ok", round(hifast, 2), round(hi178, 2), [round(v, 2) for v in lo], round(r1, 3), round(c1["n"], 2), round(100*c1["share_inbound"]), round(cH["n"], 2), round(cL["n"], 2), [round(x) for x in sh], "thresholds", F["thresholds"])
