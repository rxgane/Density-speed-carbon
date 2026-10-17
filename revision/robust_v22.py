"""Robustness of the least-cost design for a 30% cut in the reference network that finances its stock in transit."""
import os
import numpy as np
from scipy.optimize import minimize_scalar, brentq
import application_v22 as B
here = os.path.dirname(os.path.abspath(__file__)); Y = 365.0; X = 0.30
out = []
def say(s=""): print(s); out.append(s)
def best_v(n, p, cap, vset, o):
    f = lambda v: B.KC(n, v, p, **o)
    if vset is not None:
        feas = [(f(v)[0], v) for v in vset if cap is None or f(v)[1] <= cap*(1+1e-12)]
        return min(feas)[1] if feas else None
    lo, hi = p["vmin"], p["vmax"]
    vk = minimize_scalar(lambda v: f(v)[0], bounds=(lo, hi), method="bounded", options={"xatol": 1e-7}).x
    if cap is None or f(vk)[1] <= cap: return vk
    vc = minimize_scalar(lambda v: f(v)[1], bounds=(lo, hi), method="bounded", options={"xatol": 1e-7}).x
    if f(vc)[1] > cap: return None
    g = lambda v: f(v)[1] - cap
    a = lo if g(lo) <= 0 else brentq(g, lo, vc); b = hi if g(hi) <= 0 else brentq(g, vc, hi)
    return min(max(vk, a), b)
def solve(p, cap=None, vset=None, integer=True, **o):
    if integer:
        best = None
        for n in range(1, p["N"]+1):
            v = best_v(n, p, cap, vset, o)
            if v is None: continue
            k = B.KC(n, v, p, **o)[0]
            if best is None or k < best[0]: best = (k, n, v)
        return None if best is None else best[1:]
    def red(n):
        v = best_v(n, p, cap, vset, o); return 1e18 if v is None else B.KC(n, v, p, **o)[0]
    grid = np.exp(np.linspace(0, np.log(p["N"]), 400)); vals = [red(n) for n in grid]; i = int(np.argmin(vals))
    if vals[i] >= 1e18: return None
    n = minimize_scalar(red, bounds=(grid[max(i-1, 0)], grid[min(i+1, 399)]), method="bounded", options={"xatol": 1e-6}).x
    if red(n) > vals[i]: n = grid[i]          # refinement left the feasible set: keep the grid point
    return n, best_v(n, p, cap, vset, o)
def maxcut(p, base_C, vset=None, integer=False, **o):
    lo, hi = 0.0, 0.99
    for _ in range(40):
        mid = 0.5*(lo+hi)
        if solve(p, (1-mid)*base_C, vset, integer, **o) is None: hi = mid
        else: lo = mid
    return lo
def row(label, p, vset=None, integer=False, base_o=None, cap_o=None, budget=None):
    base_o = base_o or {}; cap_o = cap_o or {}
    nb, vb = solve(p, None, vset, integer, **base_o); Kb, Cb = B.KC(nb, vb, p, **base_o)
    cap = (1-X)*Cb if budget is None else budget
    sres = solve(p, cap, vset, integer, **{**base_o, **cap_o})
    if sres is None:
        mc = maxcut(p, Cb, vset, integer, **{**base_o, **cap_o}) if budget is None else None
        say(f"  {label:<44} cost min {nb:.1f} DCs {1000/vb:.2f} d | not attainable" + (f" (largest cut {100*mc:.1f}%)" if mc is not None else ""))
        return
    n, v = sres; K, C = B.KC(n, v, p, **base_o)
    say(f"  {label:<44} cost min {nb:.1f} DCs {1000/vb:.2f} d | {n:.1f} DCs {1000/v:.2f} d, cut {100*(1-C/Cb):.1f}%, cost +{100*(K/Kb-1):.1f}%")
if __name__ == "__main__":
    p = B.reference(); n0, v0 = solve(p, integer=False); C0 = B.KC(n0, v0, p)[1]
    say(f"ROBUSTNESS (continuous number of centers unless stated): least-cost design for a {100*X:.0f}% cut from each row's own cost minimum (realized caps: same expected budget as the reference)")
    row("Reference", p)
    row("Rail, intermodal, or truckload only", p, vset=(1000.0/6, 250.0, 500.0))
    row("Whole number of centers", p, integer=True)
    row("Warehouse carbon from site energy", B.reference(wh="site"))
    row("No financing of stock in transit", B.reference(financed=False))
    from statistics import NormalDist
    CVQ = (p["sigma"]*np.sqrt(7)/(p["mu"]*7))*np.sqrt((1+(p["N"]-1)*p["rho"])/p["N"]); z95 = NormalDist().inv_cdf(0.95)
    m_week, m_year = 1 + z95*CVQ, 1 + z95*CVQ/np.sqrt(52)
    say(f"  multipliers from primitives: CV_Q={CVQ:.5f}, weekly m={m_week:.5f}, annual m={m_year:.5f}")
    row("Annual cap on realized carbon, 95%", p, cap_o=dict(m=m_year), budget=(1-X)*C0)
    row("Weekly cap on realized carbon, 95%", p, cap_o=dict(m=m_week), budget=(1-X)*C0)
    for lam in (0.1, 0.25, 0.5):
        row(f"Load-sensitive delivery, lambda={lam}", p, base_o=dict(lam=lam, nref=4.0))
    with open(os.path.join(here, "robust_v22_output.txt"), "w") as fh: fh.write("\n".join(out) + "\n")
