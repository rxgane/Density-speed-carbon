"""Table 4 rows: sensitivity of the least-cost design for a 20% cut to the holding, line-haul, warehouse, and store-delivery inputs.
Writes sensitivity_v24.json; the printed table is saved as sensitivity_v24_output.txt."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize_scalar, brentq
import application_v21 as A, application_v22 as B, targets_v22 as T
CUT = 0.20
Y = 365.0; TLc, TLk = 2.26/26, 0.48/3.66*10.21/26

def lane_carbon(q, k_im, k_tl):     # refit eps(v) through intermodal (250 mi/day) and truckload (500)
    e1 = (k_tl-k_im)/(500.0**2-250.0**2); q["eps0"] = e1; q["e0"] = k_im - e1*250.0**2

V = {  # name: (kwargs for B.reference, truckload kg per truck-mile with rail fixed, intermodal carbon as share of truckload,
       #        warehouse carbon factor, store-delivery r_L)
    "Reference": ({}, None, None, 1.0, None),
    "Capital charge 10%": (dict(capital=0.10), None, None, 1.0, None),
    "Capital charge 25%": (dict(capital=0.25), None, None, 1.0, None),
    "Pallet value $10,000": (dict(value=10_000.0), None, None, 1.0, None),
    "Pallet value $50,000": (dict(value=50_000.0), None, None, 1.0, None),
    "Efficient diesel truck (1.34 kg per mile)": ({}, 1.339, None, 1.0, None),
    "Electric truck (0.74 kg per mile)": ({}, 2.0*0.371, None, 1.0, None),
    "Intermodal 45% below truckload": ({}, None, 0.55, 1.0, None),
    "Intermodal 70% below truckload": ({}, None, 0.30, 1.0, None),
    "Warehouse carbon halved": ({}, None, None, 0.5, None),
    "Warehouse carbon doubled": ({}, None, None, 2.0, None),
    "Refrigerated storage": (dict(storage="refrigerated"), None, None, 1.0, None),
    "Electric store delivery (r_L = 0.16)": ({}, None, None, 1.0, 0.44/2.75),
}
def variant(name):
    kw, tl, im, wh, rL = V[name]
    q = B.reference(switch=1.0, **kw)                  # full line-haul curves; the switchable share is applied last
    k_tl, k_im = q["e0"]+q["eps0"]*500.0**2, q["e0"]+q["eps0"]*250.0**2
    if tl is not None: lane_carbon(q, k_im, tl/26.0)              # truck end changes, intermodal end fixed
    if im is not None: lane_carbon(q, im*k_tl, k_tl)              # intermodal end changes, truck end fixed
    q["phi"] *= wh; q["eta"] *= wh
    if rL is not None: q["rL"] = rL
    return B.switchable(q, B.SWITCH)

def price_for(q, C_target):
    def C_at(tau):
        n, v = B.design(q, tau, integer=False); return B.KC(n, v, q)[1]
    return brentq(lambda t: C_at(t) - C_target, 0.0, 50.0, xtol=1e-5)

def run(name):
    q = variant(name); nK, vK = T.solve(q); K0, C0 = B.KC(nK, vK, q)
    S = A.SS(nK, vK, q); tput = q["N"]*q["mu"]; lm = q["rL"]*q["tL"]*tput*A.d(nK, q)
    wh = q["phi"]*nK + q["eta"]*S; tr = C0 - wh
    s = T.solve(q, cap=(1-CUT)*C0)
    if s is None:
        nC = float(T.cmin(q)); vC = minimize_scalar(lambda v: B.KC(nC, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
        return dict(name=name, lm_share=lm/C0, nK=nK, dK=q["D"]/vK, C0=C0*Y/1e3, wh_share=wh/C0, n30=float("nan"), d30=float("nan"),
                    cost30=float("nan"), per_t=float("nan"), price30=float("nan"), slow_max=float("nan"), both_max=1-B.KC(nC, vC, q)[1]/C0)
    n30, v30 = s; K30, C30 = B.KC(n30, v30, q)
    vR = minimize_scalar(lambda v: B.KC(nK, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
    nC = float(T.cmin(q)); vC = minimize_scalar(lambda v: B.KC(nC, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
    tau = price_for(q, (1-CUT)*C0)
    return dict(name=name, lm_share=lm/C0, nK=nK, dK=q["D"]/vK, C0=C0*Y/1e3, wh_share=wh/C0, n30=n30, d30=q["D"]/v30, cost30=K30/K0-1,
                per_t=(K30-K0)*Y/((C0-C30)*Y/1e3), price30=tau*1000, slow_max=1-B.KC(nK, vR, q)[1]/C0, both_max=1-B.KC(nC, vC, q)[1]/C0)

NAMES = list(V)
if __name__ == "__main__":
    with Pool(min(len(NAMES), max(2, os.cpu_count()-2))) as P: R = P.map(run, NAMES)
    import json
    json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sensitivity_v24.json'), 'w'))
    print(f"(cut {int(100*CUT)}%)")
    print(f"{'Variant':<52}{'cost min':>16}{'carbon t':>10}{'wareh.':>8}{'store':>7} | {f'{int(CUT*100)}% cut design':>18}{'cost':>7}{'$/t':>7}{'price':>8} | {'slow max':>9}{'both max':>9}")
    for r in R:
        print(f"{r['name']:<52}{r['nK']:6.1f} DCs {r['dK']:.2f}d{r['C0']:>10,.0f}{100*r['wh_share']:7.1f}%{100*r['lm_share']:7.1f}% | {r['n30']:5.1f} DCs {r['d30']:.2f}d +{100*r['cost30']:4.1f}%{r['per_t']:>7,.0f}{r['price30']:>8,.0f} | {100*r['slow_max']:8.1f}%{100*r['both_max']:8.1f}%")
