"""
Reference network that finances its stock in transit (application_v22.reference): carbon-target exercise for
three action sets, and the inbound-transport link to network density. Writes targets_v22.json, inbound_v22.json
and targets_v22_output.txt.
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

def KC(n, v, p): return B.KC(n, v, p)

def best_v(n, p, cap=None):
    """least-cost speed at density n with carbon at most cap (cost and carbon are convex in v)"""
    lo, hi = p["vmin"], p["vmax"]
    vk = minimize_scalar(lambda v: KC(n, v, p)[0], bounds=(lo, hi), method="bounded", options={"xatol": 1e-7}).x
    if cap is None or KC(n, vk, p)[1] <= cap: return vk
    vc = minimize_scalar(lambda v: KC(n, v, p)[1], bounds=(lo, hi), method="bounded", options={"xatol": 1e-7}).x
    if KC(n, vc, p)[1] > cap: return None
    g = lambda v: KC(n, v, p)[1] - cap
    a = lo if g(lo) <= 0 else brentq(g, lo, vc); b = hi if g(hi) <= 0 else brentq(g, vc, hi)
    return min(max(vk, a), b)

def solve(p, cap=None, mode="both", nfix=None, vfix=None, integer=True):
    if mode == "replen":
        v = best_v(nfix, p, cap); return None if v is None else (nfix, v)
    if integer:
        best = None
        for n in range(1, p["N"]+1):
            if mode == "network":
                v = vfix
                if cap is not None and KC(n, v, p)[1] > cap: continue
            else:
                v = best_v(n, p, cap)
                if v is None: continue
            k = KC(n, v, p)[0]
            if best is None or k < best[0]: best = (k, n, v)
        return None if best is None else best[1:]
    def red(n):
        v = best_v(n, p, cap)
        return 1e18 if v is None else KC(n, v, p)[0]
    grid = np.exp(np.linspace(0, np.log(p["N"]), 400)); vals = [red(n) for n in grid]; i = int(np.argmin(vals))
    if vals[i] >= 1e18: return None
    n = minimize_scalar(red, bounds=(grid[max(i-1, 0)], grid[min(i+1, 399)]), method="bounded", options={"xatol": 1e-6}).x
    return n, best_v(n, p, cap)

def set_rt(p, rt):
    e_min = p["e0"] + p["eps0"]*p["vmin"]**2; p["eps0"] = rt*p["t1"]; p["e0"] = e_min - p["eps0"]*p["vmin"]**2

if __name__ == "__main__":
    p = B.reference()
    n0, v0 = solve(p); K0, C0 = KC(n0, v0, p)
    say(f"REFERENCE (financed): {n0} DCs, {1000/v0:.2f} days, cost ${K0*Y/1e6:.1f}M, carbon {C0*Y/1e3:,.0f} t")
    xs = np.round(np.arange(0, 0.5001, 0.0025), 4); rows = []
    for x in xs:
        cap = (1-x)*C0 if x > 0 else None; r = {"x": float(x)}
        for mode in ("replen", "network", "both"):
            sres = solve(p, cap, mode, nfix=n0, vfix=v0)
            r[mode] = None if sres is None else dict(n=int(sres[0]), days=1000/sres[1], K=KC(*sres, p)[0]*Y, C=KC(*sres, p)[1]*Y)
        sres = solve(p, cap, "both", integer=False)
        r["both_cont"] = None if sres is None else dict(n=float(sres[0]), days=1000/sres[1], K=KC(*sres, p)[0]*Y)
        rows.append(r)
    json.dump({"ref": dict(n=n0, days=1000/v0, K=K0*Y, C=C0*Y), "rows": rows}, open(os.path.join(here, "targets_v22.json"), "w"))
    for mode in ("replen", "network", "both"):
        last = max(r["x"] for r in rows if r[mode]); say(f"  largest attainable cut, {mode}: {100*last:.2f}% (grid)")
    first = next(r["x"] for r in rows if r["both"] and r["both"]["n"] > n0)
    say(f"  both decisions first add a DC at a cut of {100*first:.2f}% (grid)")
    for x in (0.10, 0.20, 0.25, 0.30, 0.35, 0.40):
        r = next(q for q in rows if abs(q["x"]-x) < 1e-9); line = []
        for mode in ("replen", "network", "both"):
            q = r[mode]
            line.append(f"{mode}: n.a." if q is None else f"{mode}: {q['n']} DCs {q['days']:.2f} d cost +{100*(q['K']/(K0*Y)-1):.1f}% avg ${(q['K']-K0*Y)/(C0*Y-q['C'])*1e3:,.0f}/t")
        say(f"  cut {100*x:.0f}%: " + " | ".join(line))

    say("\nINBOUND-TRANSPORT LINK (financed), carbon intensity of faster service varied with its cost schedule fixed")
    XS = np.round(np.arange(0.10, 2.0001, 0.05), 3)           # kg CO2e per loaded vehicle-mile, 20 pallets a load
    def cmin(q):
        best = None
        for n in range(1, q["N"]+1):
            v = minimize_scalar(lambda v: KC(n, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded", options={"xatol": 1e-6}).x
            c = KC(n, v, q)[1]
            if best is None or c < best[0]: best = (c, n, v)
        return best[1]
    res = {"x": XS.tolist(), "top": {}, "bottom": {}}
    TOP = {"lower": dict(F_year=6.0e6, over={}), "higher": dict(F_year=1.0e6, over=dict(sigma=150/np.sqrt(7), rho=0.1, z=2.326))}
    for env, e in TOP.items():
        base = B.reference(storage="refrigerated", value=1250.0, spoil_week=0.005, embodied=0.803, lane_kind="labor", F_year=e["F_year"])
        base.update(e["over"]); nb = solve(base)[0]
        res["top"][env] = {"omega": A.omega(base), "nK": int(nb), "curves": {}}
        for rt in (0.03, 0.45, 1.78):
            ys = []
            for xv in XS:
                q = B.reference(storage="refrigerated", value=1250.0, spoil_week=0.005, embodied=0.803, lane_kind="labor",
                                F_year=e["F_year"], rL=xv/20/0.1375)
                q.update(e["over"]); set_rt(q, rt); ys.append(int(cmin(q)))
            res["top"][env]["curves"][str(rt)] = ys
        say(f"  perishables, {env} pooling (omega {A.omega(base):.3f}, cost-min {nb} DCs): carbon-min DCs at diesel truck "
            + ", ".join(f"r_t={rt}: {res['top'][env]['curves'][str(rt)][list(XS).index(1.35)]}" for rt in (0.03, 0.45, 1.78)))
    CUT = 0.25
    BOT = {"lower": dict(sigma=50/np.sqrt(7), rho=0.6), "higher": dict(sigma=150/np.sqrt(7), rho=0.1, z=2.326)}
    for env, over in BOT.items():
        base = B.reference(); base.update(over); nb = solve(base)[0]
        res["bottom"][env] = {"omega": A.omega(base), "nK": int(nb), "curves": {}, "cut": CUT}
        for rt in (0.03, 0.45, 1.78):
            ys = []
            for xv in XS:
                q = B.reference(rL=xv/20/0.1375); q.update(over); set_rt(q, rt)
                nq, vq = solve(q); Cq = KC(nq, vq, q)[1]; sres = solve(q, cap=(1-CUT)*Cq)
                ys.append(None if sres is None else int(sres[0]))
            res["bottom"][env]["curves"][str(rt)] = ys
        say(f"  reference, {env} pooling (omega {A.omega(base):.3f}, cost-min {nb} DCs), {100*CUT:.0f}% cut, DCs at diesel truck "
            + ", ".join(f"r_t={rt}: {res['bottom'][env]['curves'][str(rt)][list(XS).index(1.35)]}" for rt in (0.03, 0.45, 1.78)))
    json.dump(res, open(os.path.join(here, "inbound_v22.json"), "w"))
    with open(os.path.join(here, "targets_v22_output.txt"), "w") as fh: fh.write("\n".join(out) + "\n")
