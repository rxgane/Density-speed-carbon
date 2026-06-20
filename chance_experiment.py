"""Chance-constrained carbon cap illustration: scale the transport component by
m = q_{1-beta}/(N mu) and recompute the carbon-minimizing design."""
import numpy as np
from statistics import NormalDist
from dataclasses import replace
from numerics import Params, minimize_obj, minimize_carbon, cost_K, carbon_C

p = Params()
nK, vK, KK, CK = minimize_obj(p, 0.0)                  # cost optimum (carbon-rate independent)
CVQ = (p.sigma / p.mu) * np.sqrt((1 + (p.N - 1) * p.rho) / p.N)

print(f"CV_Q = {CVQ:.4f}")
print(f"cost optimum: n_K={nK:.2f}, v_K={vK:.2f}, K={KK:,.0f}, C={CK:,.0f}")
print(f"base ratios:  eta/h={p.eta/p.h:.3f}  eps0/t1={p.eps0/p.t1:.3f}  "
      f"epsL/tL={p.epsL/p.tL:.3f}  phi/F={p.phi/p.F:.3f}")
print(f"{'rule':10} {'m':>6} {'m*e0/t1':>8} {'m*eL/tL':>8} {'v_C':>6} {'n_C':>6} "
      f"{'dir':>22} {'costprem':>9} {'Cred':>7}")
print("-" * 92)
for label, conf in [("Expected", 0.50), ("90%", 0.90), ("95%", 0.95), ("99%", 0.99)]:
    z = NormalDist().inv_cdf(conf)
    m = 1.0 + z * CVQ
    pm = replace(p, eps0=m * p.eps0, epsL=m * p.epsL)
    nC, vC, _, _ = minimize_carbon(pm)
    Kc = cost_K(nC, vC, p)               # cost premium (cost rates unscaled)
    Cc = carbon_C(nC, vC, p)             # physical carbon at that design (unscaled, for reduction)
    sp_th, de_th = m * p.eps0 / p.t1, m * p.epsL / p.tL
    dv = "faster" if vC > vK + 1e-3 else ("slower" if vC < vK - 1e-3 else "same")
    dn = "more dispersed" if nC > nK + 1e-2 else ("more consolidated" if nC < nK - 1e-2 else "same")
    print(f"{label:10} {m:6.3f} {sp_th:8.3f} {de_th:8.3f} {vC:6.2f} {nC:6.2f} "
          f"{dv+', '+dn:>22} {100*(Kc/KK-1):8.2f}% {100*(1-Cc/CK):6.1f}%")
