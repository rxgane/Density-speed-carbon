"""Figure 4 data: carbon-minimizing density (perishables) and the least-cost design for a 25% cut (two versions
of the reference network) as the carbon intensity of faster inbound service r_t varies. Writes inbound_data.json."""
import sys, json, os, itertools
REV = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REV)
import numpy as np
from multiprocessing import Pool

ENV = {  # exactly the scenarios behind the current Figure 3
    ("perishables", "lower"): dict(kw=dict(storage="refrigerated", value=1250.0, spoil_week=0.005, embodied=0.803, lane_kind="labor", F_year=6.0e6), over={}),
    ("perishables", "higher"): dict(kw=dict(storage="refrigerated", value=1250.0, spoil_week=0.005, embodied=0.803, lane_kind="labor", F_year=1.0e6), over=dict(sigma=150/np.sqrt(7), rho=0.1, z=2.326)),
    ("reference", "lower"): dict(kw={}, over=dict(sigma=50/np.sqrt(7), rho=0.6)),
    ("reference", "higher"): dict(kw={}, over=dict(sigma=150/np.sqrt(7), rho=0.1, z=2.326)),
}
CUT = 0.25

def make(prod, pool, x=None, rt=None):
    import application_v22 as B, targets_v22 as T
    e = ENV[(prod, pool)]; kw = dict(e["kw"])
    if x is not None: kw["rL"] = x/20/0.1375
    q = B.reference(**kw); q.update(e["over"])
    if rt is not None: T.set_rt(q, rt)
    return q

def parts(n, v, q):
    import application_v21 as A
    S = A.SS(n, v, q); t = q["N"]*q["mu"]
    return dict(facilities=q["phi"]*n, holding=q["eta"]*S, inbound=(q["e0"]+q["eps0"]*v**2)*t*q["D"], lastmile=q["rL"]*q["tL"]*t*A.d(n, q))

def point(args):
    prod, pool, x, rt = args
    import application_v21 as A, application_v22 as B, targets_v22 as T
    from scipy.optimize import minimize_scalar
    q = make(prod, pool, x, rt)
    nK, vK = T.solve(q); KK, CK = B.KC(nK, vK, q)
    # carbon minimum (n, v)
    nC = float(T.cmin(q))
    vC = minimize_scalar(lambda v: B.KC(nC, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded", options={"xatol": 1e-7}).x
    KC_, CC = B.KC(nC, vC, q)
    r = dict(product=prod, pooling=pool, x=x, rt=rt, omega=float(A.omega(q)),
             costmin=dict(n=float(nK), days=1000/vK, K=KK*365, C=CK*365),
             carbonmin=dict(n=nC, days=1000/vC, K=KC_*365, C=CC*365, max_cut=1-CC/CK, at_slowest=bool(vC < q["vmin"]*1.001)),
             cut=None)
    s = T.solve(q, cap=(1-CUT)*CK)
    if s is not None:
        n, v = s; K, C = B.KC(n, v, q); p0, p1 = parts(nK, vK, q), parts(n, v, q)
        dC = {k: (p1[k]-p0[k])*365/1000 for k in p0}  # t/yr
        r["cut"] = dict(n=float(n), days=1000/v, K=K*365, C=C*365, cost_increase=K/KK-1, change_t_per_yr=dC,
                        share_inbound=-dC["inbound"]/((CK-C)*365/1000))
    return r

if __name__ == "__main__":
    XS = sorted(set(np.round(np.arange(0.10, 2.0001, 0.05), 3).tolist() + [0.44, 1.375]))
    RTS = sorted(set(np.round(np.exp(np.linspace(np.log(0.01), np.log(3.0), 41)), 4).tolist() + [0.03, 0.45, 1.78]))
    jobs = [(p, l, x, rt) for (p, l) in ENV for x in XS for rt in (0.03, 0.45, 1.78)]
    jobs += [(p, l, x, rt) for (p, l) in ENV for x in (0.44, 1.375) for rt in RTS]
    jobs = list(dict.fromkeys(jobs))
    with Pool(max(2, os.cpu_count()-2)) as pool: res = pool.map(point, jobs, chunksize=4)
    base = {f"{p}/{l}": point((p, l, None, None)) for (p, l) in ENV}  # each scenario at its own calibrated r_L and r_t
    out = dict(notes=dict(
        x="delivery (last-mile) emissions per loaded vehicle-mile, kg CO2e; model r_L = x/(20*0.1375); electric truck 0.44, diesel truck 1.375 (calibrated default for both products)",
        rt="r_t = carbon intensity of faster inbound service, kg CO2e per extra dollar spent on speed (epsilon_0/t_1); set_rt keeps cost schedule and emissions at the slowest speed fixed. 0.03 = calibrated perishables (team drivers, speed bought with labor); 1.78 = calibrated reference network (rail vs truckload, speed bought with fuel)",
        cut="least-cost (n, v) whose emissions are 25% below that scenario's own cost minimum; null = cut not attainable",
        carbonmin="design minimizing emissions; at_slowest = line haul at slowest allowed speed (4 days); max_cut = deepest attainable cut",
        scenarios="perishables = refrigerated dairy with spoilage, solo/team drivers (current Fig 3 top row); reference = $50k/pallet goods, rail-or-truckload, financed in transit (bottom row). lower/higher = pooling at stake (omega)"),
        x_grid=XS, rt_grid=RTS, base=base, points=res)
    json.dump(out, open(os.path.join(REV, "inbound_data.json"), "w"))
    print(len(res), "points")
    for k, b in base.items(): print(k, "omega %.3f" % b["omega"], "costmin %.2f DCs %.2f d" % (b["costmin"]["n"], b["costmin"]["days"]),
                                     "carbonmin %.2f DCs %.2f d maxcut %.0f%%" % (b["carbonmin"]["n"], b["carbonmin"]["days"], 100*b["carbonmin"]["max_cut"]),
                                     "cut25", None if b["cut"] is None else "%.2f DCs %.2f d +%.1f%% inbound share %.0f%%" % (b["cut"]["n"], b["cut"]["days"], 100*b["cut"]["cost_increase"], 100*b["cut"]["share_inbound"]))
