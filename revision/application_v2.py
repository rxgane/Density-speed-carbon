"""
v2 calibrated application in real units: a US store-replenishment network for ambient packaged goods.

Model equations are those of the paper (eq. K, C, J):
    K(n,v) = F n + h SS(n,v) + t(v) N mu D + t_L N mu d(n)
    C(n,v) = phi n + eta SS(n,v) + eps(v) N mu D + eps_L N mu d(n)
    SS(n,v) = z sigma sqrt(D/v) sqrt(N Psi(n)),  Psi(n) = N rho + n (1-rho),  d(n) = delta n^{-1/2}
with t(v) = t0 + t1 v^2 and eps(v) = e0 + eps0 v^2 fitted through two real speed options on each lane.
The additive constants t0 and e0 do not depend on n or v (throughput is fixed), so they change reported
totals but not the optimal design.

Units: period = day; distance = miles; quantity = pallets; money = US$; carbon = kg CO2e.
Every carbon rate is (carbon intensity per dollar from calibration_v2.py) x (cost rate), so the
application inherits the source-traced ratios.

Inputs and sources (see calibration_v2.py for URLs):
  network      N = 40 metropolitan markets, 250 pallets/week each (10,000 pallets/week nationally),
               weekly CV 0.40, cross-market correlation 0.30, cycle service 95%
  geography    inbound distance D = 1,000 miles; DC-to-market distance d(n) = 700 n^{-1/2} miles
               (continuous approximation: about 700 miles from one national DC, about 250 miles from eight)
  facility     F = $3.0M per DC-year (rent about 40% at $10.12/ft2 on 120,000 ft2 [CW]); phi = 0.098 kg/$ x F
  holding      pallet value $1,500; capital charge 10%/yr; storage $240/pallet-yr [PAL][CARRY];
               storage energy 30 kg CO2e/pallet-yr [LIN, ESTAR]
  line haul    truckload $2.26/mile [ATRI], 26 pallets/truck; truck 1.34 kg/mile [ATRI, EPA-HUB]
               fuel lane:  intermodal 250 mi/day, 20% cheaper, 60% less carbon [IMX]; truck 500 mi/day
               labor lane: solo truck 500 mi/day; team 1,000 mi/day at +$0.35/mile, same fuel [TEAM]
  last mile    store delivery by truck, 20 pallets/load, $2.75/mile; carbon per dollar 0.50 (diesel)
               or 0.16 (electric) from calibration_v2.py
"""
import os
import numpy as np
from scipy.optimize import minimize_scalar
from statistics import NormalDist

out = []
def say(s=""): print(s); out.append(s)

YEAR = 365.0
base = dict(N=40, mu=250/7, sigma=100/np.sqrt(7), rho=0.30, z=NormalDist().inv_cdf(0.95),
            D=1000.0, delta=700.0,
            F=3.0e6/YEAR, rF=0.098,
            h=(0.10*1500 + 240)/YEAR, eta=30/YEAR,
            tL=2.75/20, rL=0.50)

def lane(kind):
    TLc, TLk = 2.26/26, 0.48/3.66*10.21/26          # $/pallet-mile and kg/pallet-mile, truckload
    if kind == "fuel":                               # intermodal (250 mi/day) <-> truckload (500 mi/day)
        v1, v2 = 250.0, 500.0
        c1, c2 = 0.8*TLc, TLc
        k1, k2 = 0.4*TLk, TLk
        vb = (200.0, 550.0)
    else:                                            # solo (500) <-> team (1,000 mi/day), same tractor
        v1, v2 = 500.0, 1000.0
        c1, c2 = TLc, TLc + 0.35/26
        k1, k2 = TLk, TLk + 0.03*0.35/26             # r_t = 0.03: negligible extra fuel
        vb = (450.0, 1100.0)
    t1 = (c2-c1)/(v2**2-v1**2); t0 = c1 - t1*v1**2
    e1 = (k2-k1)/(v2**2-v1**2); e0 = k1 - e1*v1**2
    return dict(t0=t0, t1=t1, e0=e0, eps0=e1, vmin=vb[0], vmax=vb[1])

def Psi(n, p): return p["N"]*p["rho"] + n*(1-p["rho"])
def SS(n, v, p): return p["z"]*p["sigma"]*np.sqrt(p["D"]/v)*np.sqrt(p["N"]*Psi(n, p))
def d(n, p): return p["delta"]*n**-0.5
def K(n, v, p):
    return p["F"]*n + p["h"]*SS(n, v, p) + (p["t0"]+p["t1"]*v**2)*p["N"]*p["mu"]*p["D"] + p["tL"]*p["N"]*p["mu"]*d(n, p)
def C(n, v, p):
    return (p["rF"]*p["F"])*n + p["eta"]*SS(n, v, p) + (p["e0"]+p["eps0"]*v**2)*p["N"]*p["mu"]*p["D"] \
           + (p["rL"]*p["tL"])*p["N"]*p["mu"]*d(n, p)
def J(n, v, p, tau): return K(n, v, p) + tau*C(n, v, p) if np.isfinite(tau) else C(n, v, p)

def vopt(n, p, tau):
    return minimize_scalar(lambda v: J(n, v, p, tau), bounds=(p["vmin"], p["vmax"]), method="bounded",
                           options={"xatol": 1e-9}).x
def opt(p, tau):
    ns = np.exp(np.linspace(0, np.log(p["N"]), 300)); vals = [J(n, vopt(n, p, tau), p, tau) for n in ns]
    i = int(np.argmin(vals))
    r = minimize_scalar(lambda n: J(n, vopt(n, p, tau), p, tau), bounds=(ns[max(i-1, 0)], ns[min(i+1, 299)]),
                        method="bounded", options={"xatol": 1e-9})
    return r.x, vopt(r.x, p, tau)

def describe(p, n, v):
    return dict(n=n, v=v, L=p["D"]/v, SS=SS(n, v, p), K=K(n, v, p)*YEAR/1e6, C=C(n, v, p)*YEAR/1e3)

prices = [("cost minimum", 0.0), ("$50/t", 0.05), ("$100/t", 0.10), ("$200/t", 0.20), ("carbon minimum", np.inf)]

def cost_side(p):
    """pooling weight and speed bands of Proposition prop:classify at the cost-minimizing design"""
    nK, vK = opt(p, 0.0)
    G = p["z"]*p["sigma"]*np.sqrt(p["D"]/vK)*np.sqrt(p["N"])*(1-p["rho"])/(2*np.sqrt(Psi(nK, p)))
    return nK, vK, p["h"]*G/(p["F"]+p["h"]*G), np.sqrt(Psi(nK, p)/p["N"]), np.sqrt(Psi(nK, p)/Psi(1.0, p))

if __name__ == "__main__":
    nK, vK, om, lo, hi = cost_side(dict(base, **lane("fuel")))
    say(f"COST SIDE (ambient goods, fuel lane): n_K = {nK:.2f}, v_K = {vK:.0f} mi/day, omega = {om:.3f}, "
        f"speed bands {lo:.3f} / {hi:.3f}")

    say("\nTABLE: cost- and carbon-minimizing designs (energy-only holding ratio 0.077 throughout)")
    say(f"  {'holding $/pal-yr':>16} {'speed bought with':<18}{'store delivery':<15}{'DCs: cost -> carbon':>21}"
        f"{'mi/day: cost -> carbon':>24}{'carbon at min':>15}{'cost at min':>13}{'carbon at $100/t':>18}")
    for hy in (390, 3000):
        for lk, lbl in (("fuel", "fuel (IM/truck)"), ("labor", "labor (team)")):
            for rL, lm in ((0.50, "diesel"), (0.16, "electric")):
                p = dict(base, rL=rL, h=hy/YEAR, eta=0.077*hy/YEAR, **lane(lk))
                nK, vK = opt(p, 0.0); nC, vC = opt(p, np.inf); n1, v1 = opt(p, 0.10)
                say(f"  {hy:>16} {lbl:<18}{lm:<15}{f'{nK:.1f} -> {nC:.1f}':>21}{f'{vK:.0f} -> {vC:.0f}':>24}"
                    f"{100*(C(nC, vC, p)/C(nK, vK, p)-1):14.1f}%{100*(K(nC, vC, p)/K(nK, vK, p)-1):12.1f}%"
                    f"{100*(C(n1, v1, p)/C(nK, vK, p)-1):17.2f}%")

    say("\nPRICE PATH (ambient goods, fuel lane, diesel store delivery)")
    p = dict(base, **lane("fuel")); nK, vK = opt(p, 0.0)
    for nm_, tau in prices:
        n, v = (nK, vK) if tau == 0 else opt(p, tau); r = describe(p, n, v)
        say(f"  {nm_:<16} DCs {r['n']:5.2f}  speed {r['v']:4.0f} mi/day ({r['L']:.1f} days)  "
            f"cost ${r['K']:.2f}M/yr  carbon {r['C']:,.0f} t/yr")

    say("\nROUNDING the DC count to the best integer")
    for lk in ("fuel", "labor"):
        for rL in (0.50, 0.16):
            q = dict(base, rL=rL, **lane(lk))
            for lbl, tau, f in (("cost", 0.0, K), ("carbon", np.inf, C)):
                n, v = opt(q, tau); best = f(n, v, q)
                bv, bk = min((f(k, vopt(k, q, tau), q), k) for k in (int(np.floor(n)), int(np.ceil(n))) if 1 <= k <= q["N"])
                say(f"  {lk:<5} r_L={rL:.2f} {lbl:<6}: n* = {n:5.2f} -> {bk:2d}, penalty {100*(bv/best-1):.3f}%")

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "application_v2_output.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")
