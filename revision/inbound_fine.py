"""Figure 4: fine r_t grid and the exact r_t at which the 25% cut becomes attainable. Writes inbound_fine.json."""
import sys, os, json, numpy as np
from multiprocessing import Pool
R = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, R)
import inbound_data as BD
X = 1.375

def maxcut(prod, pool, rt):
    import application_v22 as B, targets_v22 as T
    from scipy.optimize import minimize_scalar
    q = BD.make(prod, pool, X, rt); nK, vK = T.solve(q); CK = B.KC(nK, vK, q)[1]
    nC = float(T.cmin(q)); vC = minimize_scalar(lambda v: B.KC(nC, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded", options={"xatol": 1e-7}).x
    return 1 - B.KC(nC, vC, q)[1]/CK

def threshold(args):
    from scipy.optimize import brentq
    prod, pool, lo, hi = args
    return float(np.exp(brentq(lambda lr: maxcut(prod, pool, float(np.exp(lr))) - BD.CUT, np.log(lo), np.log(hi), xtol=1e-6)))

if __name__ == "__main__":
    D = json.load(open(os.path.join(R, "inbound_data.json")))
    br = {}
    for pool in ("lower", "higher"):
        pts = sorted([p for p in D["points"] if p["product"] == "reference" and p["pooling"] == pool and p["x"] == X], key=lambda p: p["rt"])
        ok = [p["cut"] is not None for p in pts]; i = ok.index(True)
        br[pool] = (pts[i-1]["rt"], pts[i]["rt"]) if i > 0 else None
    print("brackets", br)
    with Pool(max(2, os.cpu_count()-2)) as P:
        thr = {pool: (P.apply(threshold, (("reference", pool) + br[pool],)) if br[pool] else 0.01) for pool in ("lower", "higher")}
        print("thresholds", thr)
        RT = np.exp(np.linspace(np.log(0.01), np.log(3.0), 121)).tolist()
        jobs = [("perishables", pool, X, rt) for pool in ("lower", "higher") for rt in RT + [0.03, 1.78]]
        for pool in ("lower", "higher"):
            t = thr[pool]; near = (t*np.exp(np.linspace(1e-5, np.log(1.6), 40))).tolist()
            jobs += [("reference", pool, X, rt) for rt in sorted(set([r for r in RT if r > t] + near + [0.17, 1.78]))]
        res = P.map(BD.point, jobs, chunksize=4)
    json.dump(dict(x=X, thresholds=thr, points=res), open(os.path.join(R, "inbound_fine.json"), "w")); print(len(res), "points")
