"""Figure: does the carbon intensity of faster inbound service change the number of distribution centers?
Reads inbound_v22.json (from targets_v22.py)."""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 9, "axes.linewidth": 0.6})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e2e1dd"
R = json.load(open(os.path.join(here, "inbound_v22.json"))); x = np.array(R["x"])
COL = {"0.03": "#2a78d6", "0.45": "#1baf7a", "1.78": "#eb6834"}
LS = {"0.03": "-", "0.45": (0, (1.5, 1.5)), "1.78": (0, (5, 2))}
LAB = {"0.03": r"$r_t=0.03$ (speed bought mostly with labor)", "0.45": r"$r_t=0.45$", "1.78": r"$r_t=1.78$ (speed bought mostly with fuel)"}
cut = int(round(100*R["bottom"]["lower"]["cut"]))
TITLES = {("top", "lower"): "(a) Carbon minimum, lower pooling", ("top", "higher"): "(b) Carbon minimum, higher pooling",
          ("bottom", "lower"): rf"(c) {cut}\% cut, lower pooling", ("bottom", "higher"): rf"(d) {cut}\% cut, higher pooling"}
fig, axs = plt.subplots(2, 2, figsize=(6.4, 5.3), sharex=True, gridspec_kw=dict(hspace=0.62, wspace=0.22))
for i, row in enumerate(("top", "bottom")):
    for j, env in enumerate(("lower", "higher")):
        ax = axs[i, j]; d = R[row][env]
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        ax.spines["left"].set_color(MUTED); ax.spines["bottom"].set_color(MUTED); ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
        for xv in (0.44, 1.375): ax.axvline(xv, color=GRID, lw=0.8, zorder=0)
        for rt in ("1.78", "0.45", "0.03"):
            y = np.array([np.nan if v is None else v for v in d["curves"][rt]], dtype=float)
            ax.plot(x, y, drawstyle="steps-mid", color=COL[rt], ls=LS[rt], lw=2.0 if rt == "0.45" else 1.7, label=LAB[rt], zorder=3 if rt == "0.45" else 2)
        ax.axhline(d["nK"], color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
        ax.text(x[-1], d["nK"], "cost\nminimum", ha="right", va="bottom", fontsize=7, color=MUTED, linespacing=0.9)
        ax.text(0.0, 1.03, TITLES[(row, env)] + r" ($\omega=%.2f$)" % d["omega"], transform=ax.transAxes, ha="left", va="bottom", fontsize=8.3)
        ax.set_xlim(0.1, 2.0)
        top = np.nanmax([np.nanmax(np.array([np.nan if v is None else v for v in d["curves"][k]], dtype=float)) for k in d["curves"]])
        ax.set_ylim(0, max(10, top*1.1))
axs[0, 0].text(0.44, axs[0, 0].get_ylim()[1]*0.98, "electric\ntruck", ha="center", va="top", fontsize=6.5, color=MUTED, linespacing=0.9)
axs[0, 0].text(1.375, axs[0, 0].get_ylim()[1]*0.98, "diesel\ntruck", ha="center", va="top", fontsize=6.5, color=MUTED, linespacing=0.9)
for ax in axs[1]: ax.set_xlabel("Delivery emissions per loaded\nvehicle-mile (kg CO$_2$e)")
axs[0, 0].set_ylabel("Distribution centers\nat the carbon minimum"); axs[1, 0].set_ylabel("Distribution centers that\nmeet the cut at least cost")
fig.text(0.06, 0.935, r"\textit{Spoiling perishables; solo or team drivers}", ha="left", va="bottom", fontsize=8.5)
fig.text(0.06, 0.475, r"\textit{Reference network: high-value goods; rail or truckload}", ha="left", va="bottom", fontsize=8.5)
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h[::-1], l[::-1], loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=3, frameon=False, fontsize=7.6, handlelength=2.4)
fig.savefig(os.path.join(here, "..", "fig_inbound.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(here, "fig_inbound.png"), dpi=200, bbox_inches="tight")
print("ok")
