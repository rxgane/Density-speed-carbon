"""Table 5: sensitivity of the least-cost design for a 30% cut to the carbon inputs and the storage environment.
Writes sensitivity_v24.json; the printed table is saved as sensitivity_v24_output.txt."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize_scalar, brentq
import application_v21 as A, application_v22 as B, targets_v22 as T
Y = 365.0; TLc, TLk = 2.26/26, 0.48/3.66*10.21/26

def lane_carbon(q, k_im, k_tl):     # refit eps(v) through intermodal (250 mi/day) and truckload (500)
    e1 = (k_tl-k_im)/(500.0**2-250.0**2); q["eps0"] = e1; q["e0"] = k_im - e1*250.0**2

def variant(name):
    q = B.reference(storage="refrigerated") if name.startswith("Refrigerated") else B.reference()
    if name == "Warehouse carbon from site energy": q = B.reference(wh="site")
    if name == "Warehouse carbon at Lineage's 0.28 kg/$": q["phi"] *= 0.28/0.244; q["eta"] *= 0.28/0.244
    if name == "Truckload emission factor -20%": lane_carbon(q, 0.8*0.4*TLk, 0.8*TLk)
    if name == "Truckload emission factor +20%": lane_carbon(q, 1.2*0.4*TLk, 1.2*TLk)
    if name == "Intermodal 45% less carbon than truckload": lane_carbon(q, 0.55*TLk, TLk)
    if name == "Intermodal 70% less carbon than truckload": lane_carbon(q, 0.30*TLk, TLk)
    if name == "Store delivery r_L = 0.30 (low end)": q["rL"] = 0.30
    if name == "Store delivery r_L = 0.60 (high end)": q["rL"] = 0.60
    if name == "Store delivery by electric truck (r_L = 0.16)": q["rL"] = 0.44/2.75
    return q

def price_for(q, C_target):
    def C_at(tau):
        n, v = B.design(q, tau, integer=False); return B.KC(n, v, q)[1]
    return brentq(lambda t: C_at(t) - C_target, 0.0, 50.0, xtol=1e-5)

def run(name):
    q = variant(name); nK, vK = T.solve(q); K0, C0 = B.KC(nK, vK, q)
    S = A.SS(nK, vK, q); tput = q["N"]*q["mu"]; lm = q["rL"]*q["tL"]*tput*A.d(nK, q)
    wh = q["phi"]*nK + q["eta"]*S; tr = C0 - wh
    s = T.solve(q, cap=0.7*C0); n30, v30 = s; K30, C30 = B.KC(n30, v30, q)
    vR = minimize_scalar(lambda v: B.KC(nK, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
    nC = float(T.cmin(q)); vC = minimize_scalar(lambda v: B.KC(nC, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
    tau = price_for(q, 0.7*C0)
    return dict(name=name, lm_share=lm/C0, nK=nK, dK=q["D"]/vK, C0=C0*Y/1e3, wh_share=wh/C0, n30=n30, d30=q["D"]/v30, cost30=K30/K0-1,
                per_t=(K30-K0)*Y/((C0-C30)*Y/1e3), price30=tau*1000, slow_max=1-B.KC(nK, vR, q)[1]/C0, both_max=1-B.KC(nC, vC, q)[1]/C0)

NAMES = ["Reference", "Warehouse carbon from site energy", "Warehouse carbon at Lineage's 0.28 kg/$",
         "Truckload emission factor -20%", "Truckload emission factor +20%",
         "Intermodal 45% less carbon than truckload", "Intermodal 70% less carbon than truckload",
         "Store delivery r_L = 0.30 (low end)", "Store delivery r_L = 0.60 (high end)", "Store delivery by electric truck (r_L = 0.16)",
         "Refrigerated storage (same goods)"]
if __name__ == "__main__":
    with Pool(min(len(NAMES), max(2, os.cpu_count()-2))) as P: R = P.map(run, NAMES)
    import json
    json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sensitivity_v24.json'), 'w'))
    print(f"{'Variant':<52}{'cost min':>16}{'carbon t':>10}{'wareh.':>8}{'store':>7} | {'30% cut design':>18}{'cost':>7}{'$/t':>7}{'price':>8} | {'slow max':>9}{'both max':>9}")
    for r in R:
        print(f"{r['name']:<52}{r['nK']:6.1f} DCs {r['dK']:.2f}d{r['C0']:>10,.0f}{100*r['wh_share']:7.1f}%{100*r['lm_share']:7.1f}% | {r['n30']:5.1f} DCs {r['d30']:.2f}d +{100*r['cost30']:4.1f}%{r['per_t']:>7,.0f}{r['price30']:>8,.0f} | {100*r['slow_max']:8.1f}%{100*r['both_max']:8.1f}%")
