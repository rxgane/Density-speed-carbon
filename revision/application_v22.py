"""
Version 2.2 calibrated application: the reference network finances its stock in transit.

Reference: ambient goods worth $50,000 a pallet, rail-or-truckload service on a 1,000-mile inbound lane,
diesel truck store delivery, distribution centers at $3.0M a year; the firm owns goods in transit and
carries their capital cost (h_p). All inputs come from application_v21.py (sources in calibration_v2.py).

Outputs: reference inputs and ratios; the response of speed and density to a carbon price for three action
sets (replenishment only, network only, both); the line-haul link; robustness at one carbon price; rounding.
Units: period = day; distance = miles; quantity = pallets; money = US$; carbon = kg CO2e. Prices in $/t.
"""
import os, json
import numpy as np
from scipy.optimize import minimize_scalar
import application_v21 as A

out = []
def say(s=""): print(s); out.append(s)
Y = 365.0
VALUE = 50_000.0

def reference(financed=True, **kw):
    q = dict(storage="ambient", value=VALUE, rL=0.50, lane_kind="fuel"); q.update(kw)
    p = A.network(**q)
    p["pipe_h"] = A.CAPITAL*q["value"]/Y if financed else 0.0
    return p

def KC(n, v, p, lam=0.0, nref=1.0, m=1.0):
    """daily cost and carbon; lam = load sensitivity of last-mile rates, m = transport-carbon factor"""
    fill = (n/nref)**lam if lam else 1.0
    S = A.SS(n, v, p); tput = p["N"]*p["mu"]
    K = (p["F"]*n + p["h"]*S + p["pipe_h"]*tput*p["D"]/v + (p["t0"]+p["t1"]*v**2)*tput*p["D"]
         + p["tL"]*fill*tput*A.d(n, p))
    C = (p["phi"]*n + p["eta"]*S + m*(p["e0"]+p["eps0"]*v**2)*tput*p["D"] + m*p["rL"]*p["tL"]*fill*tput*A.d(n, p))
    return K, C

def vbest(n, p, tau, vset=None, **o):
    f = lambda v: KC(n, v, p, **o)[0] + tau*KC(n, v, p, **o)[1]
    if vset is not None: return min(vset, key=f)
    return minimize_scalar(f, bounds=(p["vmin"], p["vmax"]), method="bounded", options={"xatol": 1e-7}).x

def design(p, tau, mode="both", nfix=None, vfix=None, integer=True, vset=None, **o):
    J = lambda n, v: KC(n, v, p, **o)[0] + tau*KC(n, v, p, **o)[1]
    if mode == "replen": return nfix, vbest(nfix, p, tau, vset, **o)
    if mode == "network":
        return min(range(1, p["N"]+1), key=lambda n: J(n, vfix)), vfix
    if integer:
        n = min(range(1, p["N"]+1), key=lambda n: J(n, vbest(n, p, tau, vset, **o)))
        return n, vbest(n, p, tau, vset, **o)
    g = lambda n: J(n, vbest(n, p, tau, vset, **o))
    ns = np.exp(np.linspace(0, np.log(p["N"]), 300)); i = int(np.argmin([g(n) for n in ns]))
    n = minimize_scalar(g, bounds=(ns[max(i-1, 0)], ns[min(i+1, 299)]), method="bounded", options={"xatol": 1e-7}).x
    return n, vbest(n, p, tau, vset, **o)

def capped(p, B, **o):
    """least-cost integer design with carbon (scaled by m in o) at most B"""
    best = None
    for n in range(1, p["N"]+1):
        vs = np.linspace(p["vmin"], p["vmax"], 801)
        feas = [(KC(n, v, p, **o)[0], v) for v in vs if KC(n, v, p, **o)[1] <= B]
        if feas:
            k, v = min(feas)
            if best is None or k < best[0]: best = (k, n, v)
    return None if best is None else (best[1], best[2])

def share(n, v, p):
    S = A.SS(n, v, p); P = p["N"]*p["mu"]*p["D"]/v
    return p["h"]*S/(p["h"]*S + 2*p["pipe_h"]*P)

if __name__ == "__main__":
    p = reference()
    n0, v0 = design(p, 0.0); K0, C0 = KC(n0, v0, p)
    S0 = A.SS(n0, v0, p); P0 = p["N"]*p["mu"]*p["D"]/v0
    say("REFERENCE NETWORK (financed stock in transit)")
    say(f"  inputs: N={p['N']}, mu={p['mu']*7:.0f} pallets/week, CV={p['sigma']*np.sqrt(7)/(p['mu']*7):.2f}, rho={p['rho']}, z={p['z']:.3f}, "
        f"D={p['D']:.0f} mi, delta={p['delta']:.0f} mi, v in [{p['vmin']:.0f},{p['vmax']:.0f}] mi/day")
    say(f"  t(v): ${p['t0']+p['t1']*p['vmin']**2:.4f} -> ${p['t0']+p['t1']*p['vmax']**2:.4f} per pallet-mile; "
        f"eps(v): {p['e0']+p['eps0']*p['vmin']**2:.4f} -> {p['e0']+p['eps0']*p['vmax']**2:.4f} kg per pallet-mile")
    say(f"  t_L=${p['tL']:.4f}/pallet-mile, r_L={p['rL']}; F=${p['F']*Y/1e6:.1f}M/yr, phi={p['phi']*Y/1e3:.0f} t/yr; "
        f"h=${p['h']*Y:,.0f}/pallet-yr, eta={p['eta']*Y:.0f} kg/pallet-yr, h_p=${p['pipe_h']*Y:,.0f}/pallet-yr")
    say(f"  ratios: r_h={p['eta']/p['h']:.4f} r_t={p['eps0']/p['t1']:.3f} r_L={p['rL']:.3f} r_F={p['phi']/p['F']:.4f}")
    say(f"  cost minimum: n={n0}, line haul {1000/v0:.2f} days, safety stock {S0:,.0f} pallets, in transit {P0:,.0f} pallets, "
        f"a_K={share(n0, v0, p):.3f}, omega={A.omega(p):.3f}, cost ${K0*Y/1e6:.1f}M/yr, carbon {C0*Y/1e3:,.0f} t/yr")

    say("\nCARBON PRICE RESPONSE (integer DCs): emissions cut and cost increase vs the cost minimum")
    grid = np.unique(np.round(np.logspace(np.log10(10), np.log10(5000), 61), 1))
    rows = []
    for P in grid:
        tau = P/1000.0; r = {"price": float(P)}
        for mode in ("replen", "network", "both"):
            n, v = design(p, tau, mode, nfix=n0, vfix=v0); K, C = KC(n, v, p)
            r[mode] = dict(n=int(n), days=1000/v, cut=100*(1-C/C0), cost=100*(K/K0-1))
        rows.append(r)
    for P in (25, 50, 100, 200, 300, 500, 700, 1000, 2000, 5000):
        tau = P/1000.0; line = []
        for mode in ("replen", "network", "both"):
            n, v = design(p, tau, mode, nfix=n0, vfix=v0); K, C = KC(n, v, p)
            line.append(f"{mode}: n={n} {1000/v:.2f}d cut {100*(1-C/C0):.1f}% cost +{100*(K/K0-1):.2f}%")
        say(f"  ${P:>5}/t  " + " | ".join(line))
    first5 = next(r["price"] for r in rows if r["both"]["n"] > n0)
    say(f"  lowest grid price at which the both-decision design adds a DC: ${first5:,.0f}/t")
    json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "application_v22_price.json"), "w"))

    say("\nLINE-HAUL LINK: carbon-minimizing DC count as the carbon intensity of faster service varies (cost schedule fixed)")
    def set_rt(q, rt):
        e_min = q["e0"] + q["eps0"]*q["vmin"]**2; q["eps0"] = rt*q["t1"]; q["e0"] = e_min - q["eps0"]*q["vmin"]**2
    def cmin(q):
        best = None
        for n in range(1, q["N"]+1):
            v = minimize_scalar(lambda v: KC(n, v, q)[1], bounds=(q["vmin"], q["vmax"]), method="bounded").x
            c = KC(n, v, q)[1]
            if best is None or c < best[0]: best = (c, n, v)
        return best[1], 1000/best[2]
    for rt in (0.03, 0.45, 1.78):
        q = reference(); set_rt(q, rt); n, d = cmin(q)
        say(f"  reference, r_t={rt}: carbon-minimizing DCs {n}, line haul {d:.2f} days")
    for rt in (0.03, 1.78):
        q = A.network(storage="refrigerated", value=1250.0, spoil_week=0.005, embodied=0.803, lane_kind="labor", F_year=1.0e6)
        q.update(dict(sigma=150/np.sqrt(7), rho=0.1, z=2.326)); q["pipe_h"] = 0.0; set_rt(q, rt); n, d = cmin(q)
        say(f"  spoiling dairy, solo or team, CV 0.60, rho 0.1, 99% service, F $1.0M (omega={A.omega(q):.3f}), r_t={rt}: "
            f"carbon-minimizing DCs {n}, line haul {d:.2f} days")

    PR = 1000.0; tau = PR/1000.0
    say(f"\nROBUSTNESS at ${PR:,.0f}/t (both decisions; cost and cut relative to each row's own cost minimum unless stated)")
    def row(label, q, base_o=None, o=None, vset=None, integer=True, budget=None):
        base_o = base_o or {}; o = o or {}
        nb, vb = design(q, 0.0, vset=vset, integer=integer, **base_o); Kb, Cb = KC(nb, vb, q, **base_o)
        if budget is None:
            n, v = design(q, tau, vset=vset, integer=integer, **base_o); K, C = KC(n, v, q, **base_o)
        else:
            s = capped(q, budget, **{**base_o, **o})
            if s is None:
                say(f"  {label:<46} cost min {nb:.1f} DCs {1000/vb:.2f} d | not attainable"); return
            n, v = s; K, C = KC(n, v, q, **base_o)
        say(f"  {label:<46} cost min {nb:.1f} DCs {1000/vb:.2f} d | {n:.1f} DCs {1000/v:.2f} d, cut {100*(1-C/Cb):.1f}%, cost +{100*(K/Kb-1):.1f}%")
    row("Reference", p)
    row("Implementation: rail or truckload only", p, vset=(250.0, 500.0))
    row("Implementation: continuous number of DCs", p, integer=False)
    for beta in (0.25, 0.75):
        row(f"Allocation: building energy share beta={beta}", reference(beta=beta))
    nP, vP = design(p, tau); CP = KC(nP, vP, p)[1]
    for m, lab in ((1.05, "Accounting: annual cap on realized carbon, 95%"), (1.37, "Accounting: weekly cap on realized carbon, 95%")):
        row(lab, p, o=dict(m=m), budget=CP)
    row("Extension: no financing of stock in transit", reference(financed=False))
    for lam in (0.1, 0.25, 0.5):
        row(f"Extension: load-sensitive delivery, lambda={lam}", p, base_o=dict(lam=lam, nref=float(n0)))

    say("\nROUNDING: largest penalty from rounding the DC count, cost and carbon objectives and the price grid")
    worst = 0.0
    for P in (0, 50, 100, 200, 500, 1000, 2000, 5000):
        t = P/1000.0
        nc, vc = design(p, t, integer=False); Jc = KC(nc, vc, p)[0] + t*KC(nc, vc, p)[1]
        ni, vi = design(p, t); Ji = KC(ni, vi, p)[0] + t*KC(ni, vi, p)[1]
        worst = max(worst, Ji/Jc - 1)
    say(f"  largest penalty: {100*worst:.2f}%")

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "application_v22_output.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")
