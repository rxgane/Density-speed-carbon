"""Computational regime-boundary figure: actual rho=0.30 flip boundaries vs the
no-pooling thresholds R_speed=1, R_density=1, with archetype points."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dataclasses import replace
import numerics as NM
from numerics import Params, minimize_obj, build_archetypes

plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 0.8,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
p = replace(Params(), rho=0.30)
nK, vK, _, _ = minimize_obj(p, 0.0)              # cost optimum (carbon-rate invariant)
r_in, r_fac = p.eps0 / p.t1, p.phi / p.F

def carbmin(pp, nres=240, vres=220):
    ng = np.linspace(1, pp.N, nres); vg = np.linspace(pp.vmin, pp.vmax, vres)
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    C = NM.carbon_C(NN, VV, pp)
    i, j = np.unravel_index(np.argmin(C), C.shape)
    return ng[i], vg[j]

Rs = np.linspace(0.3, 6.0, 48)   # R_speed (y): eta/h = Rs * (eps0/t1)
Rd = np.linspace(0.3, 6.0, 48)   # R_density (x): epsL/tL = Rd * (phi/F)
dN = np.zeros((len(Rs), len(Rd))); dV = np.zeros_like(dN)
for i, rs in enumerate(Rs):
    for j, rd in enumerate(Rd):
        pm = replace(p, eta=rs * r_in * p.h, epsL=rd * r_fac * p.tL)
        nC, vC = carbmin(pm)
        dN[i, j] = nC - nK; dV[i, j] = vC - vK

X, Y = np.meshgrid(Rd, Rs)
fig, ax = plt.subplots(figsize=(5.8, 4.6))
# shade faster-and-more-dispersed region (computed)
ax.contourf(X, Y, np.where((dN > 0) & (dV > 0), 1.0, 0.0), levels=[0.5, 1.5], colors=["0.92"])
# computational rho=0.30 boundaries
cb_d = ax.contour(X, Y, dN, levels=[0], colors="black", linewidths=1.6)
cb_v = ax.contour(X, Y, dV, levels=[0], colors="0.45", linewidths=1.6, linestyles="solid")
# no-pooling reference thresholds
ax.axvline(1, color="0.5", lw=0.9, ls="--"); ax.axhline(1, color="0.5", lw=0.9, ls="--")
# proxy handles for legend
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],color="black",lw=1.6,label="density boundary ($n_C=n_K$), $\\rho{=}0.30$"),
                   Line2D([0],[0],color="0.45",lw=1.6,label="speed boundary ($v_C=v_K$), $\\rho{=}0.30$"),
                   Line2D([0],[0],color="0.5",lw=0.9,ls="--",label="no-pooling thresholds $R{=}1$")],
          loc="upper left", fontsize=8, frameon=False)
# points
arch = build_archetypes(p)
pts = {"Cold-chain\npharma": (4.00, 1.90, (2.6, 2.4)),
       "E-grocery": (5.00, 1.55, (5.2, 2.3)),
       "Bulky\nlow-value": (4.33, 0.44, (4.7, 0.75)),
       "Ambient\ndry goods": (0.62, 0.37, (0.95, 0.55))}
for name, (x, y, xytext) in pts.items():
    ax.scatter([x], [y], s=40, color="black", zorder=5)
    ax.annotate(name, (x, y), xytext, fontsize=8, ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.5", shrinkA=2, shrinkB=4))
ax.scatter([0.75], [0.38], s=40, marker="s", facecolor="white", edgecolor="black", zorder=5)
ax.annotate("baseline", (0.75, 0.38), (1.5, 0.34), fontsize=8, ha="center",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="0.5", shrinkA=2, shrinkB=4))
ax.set_xlim(0.3, 6); ax.set_ylim(0.3, 6)
ax.set_xlabel(r"$R_{\mathrm{density}}=(\varepsilon_L/t_L)/(\varphi/F)$")
ax.set_ylabel(r"$R_{\mathrm{speed}}=(\eta/h)/(\varepsilon_0/t_1)$")
fig.savefig("fig_regimes.pdf")
print("wrote fig_regimes.pdf ; nK=%.2f vK=%.2f ; dN range [%.1f,%.1f] dV range [%.2f,%.2f]"
      % (nK, vK, dN.min(), dN.max(), dV.min(), dV.max()))
