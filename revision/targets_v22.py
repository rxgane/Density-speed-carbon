"""
Reference network that finances its stock in transit (application_v22.reference), solved as the model specifies:
the number of distribution centers n and the speed v are both continuous. Carbon-target exercise for three
action sets. Writes targets_v22.json and
targets_v22_output.txt.
"""
import os, json
import numpy as np
from scipy.optimize import minimize_scalar, brentq
import application_v21 as A
import application_v22 as B

here = os.path.dirname(os.path.abspath(__file__))
Y = 365.0
out = []
def say(s=""): print(s); out.append(s)
KC = B.KC

def interval(f, cap, lo, hi):
    """feasible interval {x in [lo,hi]: f(x) <= cap} of a quasiconvex f, or None"""
    xm = minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-9}).x
    if f(xm) > cap: return None
    g = lambda x: f(x) - cap
    a = lo if g(lo) <= 0 else brentq(g, lo, xm); b = hi if g(hi) <= 0 else brentq(g, xm, hi)
    return a, b

def best_v(n, p, cap=None):
    lo, hi = p["vmin"], p["vmax"]
    vk = minimize_scalar(lambda v: KC(n, v, p)[0], bounds=(lo, hi), method="bounded", options={"xatol": 1e-7}).x
    if cap is None or KC(n, vk, p)[1] <= cap: return vk
    iv = interval(lambda v: KC(n, v, p)[1], cap, lo, hi)
    return None if iv is None else min(max(vk, iv[0]), iv[1])

def best_n(v, p, cap=None):
    lo, hi = 1.0, float(p["N"])
    nk = minimize_scalar(lambda n: KC(n, v, p)[0], bounds=(lo, hi), method="bounded", options={"xatol": 1e-7}).x
    if cap is None or KC(nk, v, p)[1] <= cap: return nk
    iv = interval(lambda n: KC(n, v, p)[1], cap, lo, hi)
    return None if iv is None else min(max(nk, iv[0]), iv[1])

def solve(p, cap=None, mode="both", nfix=None, vfix=None):
    """least-cost (n, v) with carbon at most cap; n and v continuous"""
    if mode == "replen":
        v = best_v(nfix, p, cap); return None if v is None else (nfix, v)
    if mode == "network":
        n = best_n(vfix, p, cap); return None if n is None else (n, vfix)
    def red(n):
        v = best_v(n, p, cap); return 1e18 if v is None else KC(n, v, p)[0]
    grid = np.exp(np.linspace(0, np.log(p["N"]), 400)); vals = [red(n) for n in grid]; i = int(np.argmin(vals))
    if vals[i] >= 1e18: return None
    n = minimize_scalar(red, bounds=(grid[max(i-1, 0)], grid[min(i+1, 399)]), method="bounded", options={"xatol": 1e-7}).x
    if red(n) > vals[i]: n = grid[i]          # refinement left the feasible set (thin feasible interval): keep the grid point
    return n, best_v(n, p, cap)

def cmin(q):
    """carbon-minimizing (n, v), continuous"""
    def red(n):
        return minimize_scalar(lambda v: KC(n, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded", options={"xatol": 1e-7}).fun
    grid = np.exp(np.linspace(0, np.log(q["N"]), 300)); vals = [red(n) for n in grid]; i = int(np.argmin(vals))
    return minimize_scalar(red, bounds=(grid[max(i-1, 0)], grid[min(i+1, 299)]), method="bounded", options={"xatol": 1e-7}).x

def set_rt(p, rt):
    va = p["modes"][0]; e_min = p["e0"] + p["eps0"]*va**2; p["eps0"] = rt*p["t1"]; p["e0"] = e_min - p["eps0"]*va**2   # carbon of intermodal held fixed

if __name__ == "__main__":
    p = B.reference()
    n0, v0 = solve(p); K0, C0 = KC(n0, v0, p)
    say(f"REFERENCE (financed, continuous n): {n0:.2f} DCs, {1000/v0:.2f} days, cost ${K0*Y/1e6:.1f}M, carbon {C0*Y/1e3:,.0f} t")
    xs = np.round(np.arange(0, 0.6001, 0.0025), 4); rows = []
    for x in xs:
        cap = (1-x)*C0 if x > 0 else None; r = {"x": float(x)}
        for mode in ("replen", "network", "both"):
            sres = solve(p, cap, mode, nfix=n0, vfix=v0)
            if sres is None: r[mode] = None; continue
            K, C = KC(*sres, p)
            parts_K = None
            r[mode] = dict(n=float(sres[0]), days=1000/sres[1], K=K*Y, C=C*Y)
        rows.append(r)
    json.dump({"ref": dict(n=n0, days=1000/v0, K=K0*Y, C=C0*Y), "rows": rows}, open(os.path.join(here, "targets_v22.json"), "w"))
    for mode in ("replen", "network", "both"):
        last = max(r["x"] for r in rows if r[mode]); say(f"  largest attainable cut, {mode}: {100*last:.2f}% (grid)")
    for x in (0.05, 0.10, 0.20, 0.25, 0.30, 0.35, 0.40):
        r = next(q for q in rows if abs(q["x"]-x) < 1e-9); line = []
        for mode in ("replen", "network", "both"):
            q = r[mode]
            line.append(f"{mode}: n.a." if q is None else f"{mode}: {q['n']:.2f} DCs {q['days']:.2f} d cost +{100*(q['K']/(K0*Y)-1):.2f}% avg ${(q['K']-K0*Y)/(C0*Y-q['C'])*1e3:,.0f}/t")
        say(f"  cut {100*x:.0f}%: " + " | ".join(line))
    # share of the cut delivered by inbound freight and safety stock (speed) vs facilities and last mile (density), both decisions
    say("  sources of the cut, both decisions (change in t/yr: facilities, holding, inbound, last mile):")
    for x in (0.10, 0.20, 0.30):
        r = next(q for q in rows if abs(q["x"]-x) < 1e-9)["both"]
        if r is None: say(f"    cut {100*x:.0f}%: not attainable"); continue
        v = 1000/r["days"]; n = r["n"]
        S0, S1 = A.SS(n0, v0, p), A.SS(n, 1000/r["days"]*1, p)
        def comp(n_, v_):
            S = A.SS(n_, v_, p); t = p["N"]*p["mu"]
            return [p["phi"]*n_, p["eta"]*S, (p["e0"]+p["eps0"]*v_**2)*t*p["D"], p["rL"]*p["tL"]*t*A.d(n_, p)]
        d = [(b-a)*Y/1000 for a, b in zip(comp(n0, v0), comp(n, 1000/r["days"]))]
        say(f"    cut {100*x:.0f}%: " + ", ".join(f"{k} {v:+.0f}" for k, v in zip(("facilities", "holding", "inbound", "last mile"), d)))
    say("  price path (continuous):")
    for P in (50, 100, 200, 300, 500, 1000):
        tau = P/1000.0
        def red(n):
            return minimize_scalar(lambda v: KC(n, v, p)[0]+tau*KC(n, v, p)[1], bounds=(p["vmin"], p["vmax"]), method="bounded").fun
        n = minimize_scalar(red, bounds=(1, p["N"]), method="bounded", options={"xatol": 1e-7}).x
        v = minimize_scalar(lambda v: KC(n, v, p)[0]+tau*KC(n, v, p)[1], bounds=(p["vmin"], p["vmax"]), method="bounded").x
        K, C = KC(n, v, p)
        say(f"    ${P}/t: {n:.2f} DCs {1000/v:.2f} d cut {100*(1-C/C0):.1f}% cost +{100*(K/K0-1):.2f}% avg ${(K-K0)/(C0-C)*1000:,.0f}/t")

    with open(os.path.join(here, "targets_v22_output.txt"), "w") as fh: fh.write("\n".join(out) + "\n")
