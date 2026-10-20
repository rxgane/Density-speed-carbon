"""Carbon-parameter and operating sensitivities for the reference network (Section 5.3).
Writes sensitivity_v24.json and sensitivity_v24_output.txt, including the single-lever comparisons,
emissions-component changes, and the last-mile density-path illustration."""
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

def carbon_parts(n, v, q):
    """Daily emissions by component; their sum is the model's carbon objective."""
    tput = q["N"]*q["mu"]
    return dict(facilities=q["phi"]*n, holding=q["eta"]*A.SS(n, v, q),
                inbound=(q["e0"]+q["eps0"]*v**2)*tput*q["D"],
                last_mile=q["rL"]*q["tL"]*tput*A.d(n, q))

def run(name):
    q = variant(name); nK, vK = T.solve(q); K0, C0 = B.KC(nK, vK, q)
    S = A.SS(nK, vK, q); tput = q["N"]*q["mu"]; lm = q["rL"]*q["tL"]*tput*A.d(nK, q)
    wh = q["phi"]*nK + q["eta"]*S; tr = C0 - wh
    s = T.solve(q, cap=(1-CUT)*C0)
    if s is None:
        nC = float(T.cmin(q)); vC = minimize_scalar(lambda v: B.KC(nC, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
        return dict(name=name, lm_share=lm/C0, nK=nK, dK=q["D"]/vK, C0=C0*Y/1e3, wh_share=wh/C0, n30=float("nan"), d30=float("nan"),
                    cost30=float("nan"), per_t=float("nan"), price30=float("nan"), slow_max=float("nan"), both_max=1-B.KC(nC, vC, q)[1]/C0,
                    replen_cost=float("nan"), carbon_reductions=None)
    n30, v30 = s; K30, C30 = B.KC(n30, v30, q)
    vR = minimize_scalar(lambda v: B.KC(nK, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
    nC = float(T.cmin(q)); vC = minimize_scalar(lambda v: B.KC(nC, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
    tau = price_for(q, (1-CUT)*C0)
    replen = T.solve(q, cap=(1-CUT)*C0, mode="replen", nfix=nK)
    replen_cost = float("nan") if replen is None else B.KC(*replen, q)[0]/K0-1
    before, after = carbon_parts(nK, vK, q), carbon_parts(n30, v30, q)
    assert np.isclose(sum(before.values()), C0) and np.isclose(sum(after.values()), C30)
    reductions = {key: (before[key]-after[key])*Y/1e3 for key in before}
    return dict(name=name, lm_share=lm/C0, nK=nK, dK=q["D"]/vK, C0=C0*Y/1e3, wh_share=wh/C0, n30=n30, d30=q["D"]/v30, cost30=K30/K0-1,
                per_t=(K30-K0)*Y/((C0-C30)*Y/1e3), price30=tau*1000, slow_max=1-B.KC(nK, vR, q)[1]/C0, both_max=1-B.KC(nC, vC, q)[1]/C0,
                replen_cost=replen_cost, carbon_reductions=reductions)

def last_mile_path(rL=0.09, price=100.0):
    """Change last-mile carbon only; compare the initial response and carbon minimum."""
    q = B.reference(rL=rL)
    nK, vK = T.solve(q)
    nP, vP = B.design(q, price/1000, integer=False)
    nC = float(T.cmin(q))
    vhat = minimize_scalar(lambda v: B.KC(nK, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
    G = A.SS(nK, vK, q)*(1-q["rho"])/(2*A.Psi(nK, q))
    omega = q["h"]*G/(q["F"]+q["h"]*G)
    aK = B.share(nK, vK, q)
    rh, rt, rF = q["eta"]/q["h"], q["eps0"]/q["t1"], q["phi"]/q["F"]
    endpoint = (1-omega)*rF+omega*rh*np.sqrt(vK/vhat)
    initial = (1-omega)*rF+omega*(rh-(aK*rh-rt)/(6-aK))
    niK, _ = B.design(q, 0.0, integer=True)
    niP, _ = B.design(q, price/1000, integer=True)
    return dict(rL=rL, price=price, nK=nK, nP=nP, nC=nC, dK=q["D"]/vK, dP=q["D"]/vP,
                endpoint=endpoint, initial=initial, niK=niK, niP=niP)

NAMES = list(V)
if __name__ == "__main__":
    with Pool(min(len(NAMES), max(2, os.cpu_count()-2))) as P: R = P.map(run, NAMES)
    import json
    json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sensitivity_v24.json'), 'w'))
    output = []
    def say(line=""):
        print(line)
        output.append(line)
    say(f"(cut {int(100*CUT)}% from each row's own cost-minimizing emissions)")
    say(f"{'Variant':<52}{'cost min':>16}{'carbon t':>10}{'wareh.':>8}{'store':>7} | {f'{int(CUT*100)}% cut design':>18}{'cost':>7}{'$/t':>7}{'price':>8} | {'slow max':>9}{'both max':>9}")
    for r in R:
        say(f"{r['name']:<52}{r['nK']:6.1f} DCs {r['dK']:.2f}d{r['C0']:>10,.0f}{100*r['wh_share']:7.1f}%{100*r['lm_share']:7.1f}% | {r['n30']:5.1f} DCs {r['d30']:.2f}d +{100*r['cost30']:4.1f}%{r['per_t']:>7,.0f}{r['price30']:>8,.0f} | {100*r['slow_max']:8.1f}%{100*r['both_max']:8.1f}%")
    say("\nStrategy comparisons at each scenario's 20% target")
    for r in R:
        if r["carbon_reductions"] is None:
            continue
        changes = r["carbon_reductions"]
        net = sum(changes.values())
        shares = ", ".join(f"{key} {100*value/net:+.2f}%" for key, value in changes.items())
        say(f"{r['name']}: both +{100*r['cost30']:.3f}%, replenishment alone +{100*r['replen_cost']:.3f}%; {shares}")
    ref = next(r for r in R if r["name"] == "Reference")
    electric = next(r for r in R if r["name"] == "Electric truck (0.74 kg per mile)")
    say(f"\nElectric line haul: cost-minimizing emissions {electric['C0']:,.2f} t/year, "
        f"{100*(1-electric['C0']/ref['C0']):.2f}% below the reference. "
        f"Original absolute 20% target met before redesign: {electric['C0'] <= (1-CUT)*ref['C0']}. "
        "Freight costs are held fixed; vehicle acquisition is not evaluated.")
    path = last_mile_path()
    assert path["endpoint"] < path["rL"] < path["initial"]
    assert path["nP"] < path["nK"] < path["nC"]
    say(f"\nLast-mile sensitivity r_L={path['rL']}: endpoint threshold {path['endpoint']:.6f}, "
        f"initial threshold {path['initial']:.6f}; n_K={path['nK']:.6f}, "
        f"n at ${path['price']:.0f}/t={path['nP']:.6f}, n_C={path['nC']:.6f}; "
        f"whole-DC designs {path['niK']} and {path['niP']} at zero and ${path['price']:.0f}/t.")
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sensitivity_v24_output.txt'), 'w') as fh:
        fh.write("\n".join(output)+"\n")
