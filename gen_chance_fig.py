import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statistics import NormalDist
from dataclasses import replace
from numerics import Params, minimize_obj, minimize_carbon

plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 0.8,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
p = Params()
nK, vK, _, _ = minimize_obj(p, 0.0)
CVQ = (p.sigma / p.mu) * np.sqrt((1 + (p.N - 1) * p.rho) / p.N)
confs = np.linspace(0.50, 0.99, 40)
vC, nC = [], []
for c in confs:
    m = 1 + NormalDist().inv_cdf(float(c)) * CVQ
    pm = replace(p, eps0=m * p.eps0, epsL=m * p.epsL)
    n, v, _, _ = minimize_carbon(pm)
    vC.append(v); nC.append(n)

fig, ax1 = plt.subplots(figsize=(5.8, 4.0))
ax2 = ax1.twinx()
l1, = ax1.plot(confs * 100, vC, "-", color="0.15", lw=1.7, label=r"speed $v_C$")
l2, = ax2.plot(confs * 100, nC, "--", color="0.45", lw=1.7, label=r"density $n_C$")
ax1.axhline(vK, color="0.15", lw=0.8, ls=":")
ax2.axhline(nK, color="0.45", lw=0.8, ls=":")
ax1.text(50.5, vK, r" cost-optimal $v_K$", fontsize=8, color="0.15", va="center")
ax2.text(98.5, nK, r"cost-optimal $n_K$ ", fontsize=8, color="0.45", va="bottom", ha="right")
ax1.set_xlabel(r"confidence level $1-\beta$ (\%)")
ax1.set_ylabel(r"carbon-minimizing speed $v_C$")
ax2.set_ylabel(r"carbon-minimizing density $n_C$ (DCs)")
ax1.legend(handles=[l1, l2], loc="center left", fontsize=9, frameon=False)
fig.savefig("fig_chance.pdf")
print("wrote fig_chance.pdf  vK=%.2f nK=%.2f  vC[%.2f..%.2f] nC[%.2f..%.2f]"
      % (vK, nK, vC[0], vC[-1], nC[0], nC[-1]))
