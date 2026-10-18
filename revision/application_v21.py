"""
Version 2.1 calibrated application (responds to the external review):
  * storage emissions are physical (per pallet position), not scaled with product value;
  * building energy and storage cost are split explicitly between the facility (fixed) and holding
    (varies with inventory held) by stated shares beta (energy) and gamma (storage cost), which are tested;
  * speed directions are computed with the actual speed bounds of each lane and, as a check, with the two
    real modes only; 'unchanged' is reported where speed stays at a bound;
  * a compact pooling-weight comparison shows how line-haul carbon intensity shifts the density threshold;
  * a pipeline check charges in-transit stock its capital cost.

Units: period = day; distance = miles; quantity = pallets; money = US$; carbon = kg CO2e.
Sources for every input are listed in calibration_v2.py.
"""
import os
import numpy as np
from scipy.optimize import minimize_scalar
from statistics import NormalDist

out = []
def say(s=""): print(s); out.append(s)
YEAR = 365.0

# ------------------------------------------------------------------ physical and cost inputs
GRID = 0.371                                   # kg CO2e per kWh (EPA)
EUI = {"ambient": 22.7, "refrigerated": 84.1}  # kBtu per ft2 per year, site energy (ENERGY STAR); all-electric assumption
DC_FT2 = 120_000                               # ft2 per DC
BUILDING_KG = {k: v/3.412*GRID*DC_FT2 for k, v in EUI.items()}       # kg CO2e per DC-year
POSITION_KG = {"ambient": 30.0, "refrigerated": 109.0}              # building carbon per pallet position-year
STORAGE_USD = {"ambient": 240.0, "refrigerated": 480.0}             # 3PL storage rate per pallet-year
F_YEAR = {"ambient": 3.0e6, "refrigerated": 6.0e6}   # operating cost per DC-year; cold-storage rent about 2x ambient
FIXED_SHARE = 0.6                              # share of a center's operating cost fixed per site (about $1.8M of $3.0M); the rest is handling
TL_KG_SMARTWAY = 1.52                          # truckload CO2 per truck-mile, EPA SmartWay 2023 dry-van median
CAPITAL = 0.10                                 # capital charge on inventory value per year
BETA, GAMMA = 0.5, 0.5                         # central shares of building energy / storage cost that vary with inventory
EEIO_WH = 0.244                                # kg CO2e per $ of warehousing, USEEIO v1.3 (NAICS 493110 and 493120) [EEIO]

def network(storage="ambient", value=1500.0, spoil_week=0.0, embodied=0.0, beta=BETA, gamma=GAMMA,
            F_year=None, rL=0.50, lane_kind="fuel", wh="site"):
    """parameter dictionary of the calibrated network for one product, lane and last-mile vehicle.
    wh="site": site energy of a median building times the grid factor (reference); wh="eeio": the EPA input-output
    rate per dollar of warehousing applied to the facility cost and the storage charge that varies with stock."""
    F_year = F_YEAR[storage] if F_year is None else F_year
    h_year = CAPITAL*value + gamma*STORAGE_USD[storage] + 52*spoil_week*value
    if wh == "eeio":
        eta_year = EEIO_WH*gamma*STORAGE_USD[storage] + 52*spoil_week*value*embodied
        phi_year = EEIO_WH*F_year
    else:
        eta_year = beta*POSITION_KG[storage] + 52*spoil_week*value*embodied
        phi_year = (1-beta)*BUILDING_KG[storage]*(F_year/F_YEAR[storage])   # facility carbon scales with facility size
    p = dict(N=40, mu=250/7, sigma=100/np.sqrt(7), rho=0.30, z=NormalDist().inv_cdf(0.95), D=1000.0, delta=700.0,
             F=F_year/YEAR, phi=phi_year/YEAR, h=h_year/YEAR, eta=eta_year/YEAR,
             tL=2.75/20, rL=rL, pipe_h=0.0)
    p.update(lane(lane_kind))
    return p

def lane(kind):
    TLc, TLk = 2.26/26, 0.48/3.66*10.21/26          # $ and kg per pallet-mile, truckload
    if kind == "fuel":
        v1, v2, c1, c2, k1, k2 = 250.0, 500.0, 0.8*TLc, TLc, 0.4*TLk, TLk
    else:
        v1, v2, c1, c2, k1, k2 = 500.0, 1000.0, TLc, TLc + 0.35/26, TLk, TLk + 0.03*0.35/26
    t1 = (c2-c1)/(v2**2-v1**2); e1 = (k2-k1)/(v2**2-v1**2)
    # the fuel lane's curves are fitted through intermodal (250 mi/day) and truckload (500); the firm may slow
    # further, to a 6-day transit (carload rail), along the same curves
    vlo = 1000.0/6 if kind == "fuel" else v1
    return dict(t0=c1-t1*v1**2, t1=t1, e0=k1-e1*v1**2, eps0=e1, vmin=vlo, vmax=v2, modes=(v1, v2))

def Psi(n, p): return p["N"]*p["rho"] + n*(1-p["rho"])
def SS(n, v, p): return p["z"]*p["sigma"]*np.sqrt(p["D"]/v)*np.sqrt(p["N"]*Psi(n, p))
def d(n, p): return p["delta"]*n**-0.5
def K(n, v, p):
    return (p["F"]*n + p["h"]*SS(n, v, p) + (p["t0"]+p["t1"]*v**2)*p["N"]*p["mu"]*p["D"]
            + p["tL"]*p["N"]*p["mu"]*d(n, p) + p["pipe_h"]*p["N"]*p["mu"]*p["D"]/v)
def C(n, v, p):
    return (p["phi"]*n + p["eta"]*SS(n, v, p) + (p["e0"]+p["eps0"]*v**2)*p["N"]*p["mu"]*p["D"]
            + p["rL"]*p["tL"]*p["N"]*p["mu"]*d(n, p))
def J(n, v, p, tau): return K(n, v, p) + tau*C(n, v, p) if np.isfinite(tau) else C(n, v, p)

def vopt(n, p, tau, discrete=False):
    if discrete:
        return min(p["modes"], key=lambda v: J(n, v, p, tau))
    return minimize_scalar(lambda v: J(n, v, p, tau), bounds=(p["vmin"], p["vmax"]), method="bounded",
                           options={"xatol": 1e-9}).x
def opt(p, tau, discrete=False):
    ns = np.exp(np.linspace(0, np.log(p["N"]), 300)); vals = [J(n, vopt(n, p, tau, discrete), p, tau) for n in ns]
    i = int(np.argmin(vals))
    r = minimize_scalar(lambda n: J(n, vopt(n, p, tau, discrete), p, tau),
                        bounds=(ns[max(i-1, 0)], ns[min(i+1, 299)]), method="bounded", options={"xatol": 1e-9})
    return r.x, vopt(r.x, p, tau, discrete)

def omega(p):
    nK, vK = opt(p, 0.0)
    G = p["z"]*p["sigma"]*np.sqrt(p["D"]/vK)*np.sqrt(p["N"])*(1-p["rho"])/(2*np.sqrt(Psi(nK, p)))
    return p["h"]*G/(p["F"]+p["h"]*G)

def direction(a, b, tol):
    return "faster" if b > a + tol else ("slower" if b < a - tol else "unchanged")

PRODUCTS = [("Ambient packaged goods", dict(storage="ambient", value=1500.0)),
            ("High-value goods",       dict(storage="ambient", value=20_000.0)),
            ("Refrigerated dairy",     dict(storage="refrigerated", value=1250.0)),
            ("Refrigerated dairy, spoilage 0.5%/wk", dict(storage="refrigerated", value=1250.0, spoil_week=0.005, embodied=0.803))]

if __name__ == "__main__":
    say("CALIBRATED NETWORK, v2.1: physical storage emissions; shares beta=%.2f (energy) gamma=%.2f (storage cost)" % (BETA, GAMMA))
    for nm_, kw in PRODUCTS:
        p = network(**kw)
        say(f"  {nm_:<38} h=${p['h']*YEAR:7.0f}/pallet-yr  eta={p['eta']*YEAR:6.1f} kg/pallet-yr  r_h={p['eta']/p['h']:.3f}  "
            f"r_F={p['phi']/p['F']:.3f}  omega={omega(p):.3f}")

    say("\nTABLE: cost- and carbon-minimizing designs (continuous speed between the two modes of each lane)")
    say(f"  {'product':<38}{'lane':<7}{'days: cost->carbon':>20}{'speed':>11}{'DCs diesel':>14}{'DCs electric':>15}"
        f"{'C@min diesel':>14}{'K@min':>8}{'C @$100/t':>11}")
    for nm_, kw in PRODUCTS:
        for lk in ("fuel", "labor"):
            row = []
            for rL in (0.50, 0.16):
                p = network(rL=rL, lane_kind=lk, **kw)
                nK, vK = opt(p, 0.0); nC, vC = opt(p, np.inf); n1, v1 = opt(p, 0.10)
                row.append((nK, vK, nC, vC, C(nC, vC, p)/C(nK, vK, p)-1, K(nC, vC, p)/K(nK, vK, p)-1, C(n1, v1, p)/C(nK, vK, p)-1))
            nK, vK, nC, vC, dC, dK, d100 = row[0]
            say(f"  {nm_:<38}{lk:<7}{f'{1000/vK:.2f} -> {1000/vC:.2f}':>20}{direction(vK, vC, 0.5):>11}"
                f"{f'{nK:.1f} -> {nC:.1f}':>14}{f'{row[1][0]:.1f} -> {row[1][2]:.1f}':>15}"
                f"{100*dC:13.1f}%{100*dK:7.1f}%{100*d100:10.2f}%")

    say("\nCHECK: two real modes only (discrete speed), diesel delivery")
    for nm_, kw in PRODUCTS:
        for lk in ("fuel", "labor"):
            p = network(lane_kind=lk, **kw)
            nK, vK = opt(p, 0.0, discrete=True); nC, vC = opt(p, np.inf, discrete=True)
            say(f"  {nm_:<38}{lk:<7} mode at cost minimum {vK:5.0f} mi/day, at carbon minimum {vC:5.0f} mi/day ({direction(vK, vC, 0.5)})")

    say("\nCHECK: allocation shares beta (energy) and gamma (storage cost); speed and density direction at the carbon minimum, diesel delivery")
    for nm_, kw in PRODUCTS:
        for lk in ("fuel", "labor"):
            res = {}
            for beta in (0.0, 0.25, 0.5, 0.75, 1.0):
                for gamma in (0.25, 0.5, 1.0):
                    p = network(lane_kind=lk, beta=beta, gamma=gamma, **kw)
                    nK, vK = opt(p, 0.0); nC, vC = opt(p, np.inf)
                    key = (direction(vK, vC, 0.5), "more DCs" if nC > nK*1.005 else ("fewer DCs" if nC < nK*0.995 else "same DCs"))
                    res[key] = res.get(key, 0) + 1
            tot = sum(res.values())
            say(f"  {nm_:<38}{lk:<7} " + "; ".join(f"{a}, {b}: {c}/{tot}" for (a, b), c in sorted(res.items(), key=lambda x: -x[1])))

    say("\nMAGNITUDE RANGE over the energy share beta in [0.25, 0.75] (gamma 0.5), diesel delivery, fuel lane")
    for nm_, kw in PRODUCTS:
        vals = []
        for beta in (0.25, 0.5, 0.75):
            p = network(lane_kind="fuel", beta=beta, **kw); nK, vK = opt(p, 0.0); nC, vC = opt(p, np.inf); n1, v1 = opt(p, 0.10)
            vals.append((nC, 100*(C(nC, vC, p)/C(nK, vK, p)-1), 100*(K(nC, vC, p)/K(nK, vK, p)-1), 100*(C(n1, v1, p)/C(nK, vK, p)-1), nK))
        say(f"  {nm_:<38} cost-min DCs {vals[1][4]:.1f}; carbon-min DCs {min(v[0] for v in vals):.1f}-{max(v[0] for v in vals):.1f}; "
            f"carbon at min {min(v[1] for v in vals):.1f}% to {max(v[1] for v in vals):.1f}%; cost at min {min(v[2] for v in vals):.1f}% to {max(v[2] for v in vals):.1f}%; "
            f"carbon at $100/t {min(v[3] for v in vals):.2f}% to {max(v[3] for v in vals):.2f}%")

    say("\nPOOLING WEIGHT AND THE LINE-HAUL LINK: density threshold (1-w) r_F + w r_h^(4/5) r_t^(1/5) and carbon-minimizing DCs")
    say("  (refrigerated dairy with spoilage 0.5%/wk, electric store delivery; DC fixed cost varied to change the pooling weight)")
    say(f"  {'DC fixed cost $M/yr':>20}{'omega':>8}{'r_F':>7}{'thr labor':>11}{'thr fuel':>10}{'DCs cost':>10}{'DCs carbon labor':>18}{'DCs carbon fuel':>17}")
    for Fy in (6.0e6, 2.0e6, 0.8e6, 0.3e6):
        pl = network(F_year=Fy, rL=0.16, lane_kind="labor", **PRODUCTS[3][1]); pf = network(F_year=Fy, rL=0.16, lane_kind="fuel", **PRODUCTS[3][1])
        w = omega(pl); rh, rF = pl["eta"]/pl["h"], pl["phi"]/pl["F"]
        th = [(1-w)*rF + w*rh**0.8*(p_["eps0"]/p_["t1"])**0.2 for p_ in (pl, pf)]
        nK, _ = opt(pl, 0.0); nCl, _ = opt(pl, np.inf); nCf, _ = opt(pf, np.inf)
        say(f"  {Fy/1e6:20.2f}{w:8.3f}{rF:7.3f}{th[0]:11.3f}{th[1]:10.3f}{nK:10.1f}{nCl:18.1f}{nCf:17.1f}")

    say("\nPIPELINE CHECK: in-transit stock charged its capital cost (no warehouse carbon), fuel lane, diesel delivery")
    for nm_, kw in PRODUCTS[:2]:
        for pipe in (False, True):
            p = network(lane_kind="fuel", **kw)
            if pipe: p["pipe_h"] = CAPITAL*kw["value"]/YEAR
            nK, vK = opt(p, 0.0); nC, vC = opt(p, np.inf)
            say(f"  {nm_:<38} pipeline {'on ' if pipe else 'off'}: days cost {1000/vK:.2f} -> carbon {1000/vC:.2f}")

    say("\nBENCHMARK RATIOS under the allocation (central beta=gamma=0.5; range over beta in [0.25,0.75], gamma in [0.25,1],")
    say("capital charge 8-15%, pallet value range, storage energy range, storage rate range)")
    import itertools
    VAL = {"Ambient packaged goods": (1000, 1500, 2500), "High-value goods": (10_000, 20_000, 50_000),
           "Refrigerated dairy": (1000, 1250, 1500), "Refrigerated dairy, spoilage 0.5%/wk": (1000, 1250, 1500)}
    EN = {"ambient": (20, 30, 50), "refrigerated": (80, 109, 140)}; ST = {"ambient": (144, 240, 300), "refrigerated": (360, 480, 720)}
    for nm_, kw in PRODUCTS:
        st = kw["storage"]; sp = kw.get("spoil_week", 0.0); em = kw.get("embodied", 0.0)
        def rh(V, cap, beta, gamma, e, srate):
            return (beta*e + 52*sp*V*em) / (cap*V + gamma*srate + 52*sp*V)
        caps = (0.10, 0.15, 0.25) if nm_ == "High-value goods" else (0.08, 0.10, 0.15)   # capital charge, low/central/high
        central = rh(VAL[nm_][1], caps[1], 0.5, 0.5, EN[st][1], ST[st][1])
        vals = [rh(V, c, b, g, e, r) for V, c, b, g, e, r in itertools.product(VAL[nm_], caps, (0.25, 0.5, 0.75),
                                                                              (0.25, 0.5, 1.0), EN[st], ST[st])]
        say(f"  holding  {nm_:<38} {central:.3f}  ({min(vals):.3f}-{max(vals):.3f})")
    for st in ("ambient", "refrigerated"):
        kwh = EUI[st]/3.412
        def rf(beta, eui_scale, rent_share):
            building = kwh*eui_scale*GRID*DC_FT2
            F = (10.12*(2 if st == "refrigerated" else 1))*DC_FT2/rent_share
            return (1-beta)*building/F
        central = (1-0.5)*BUILDING_KG[st]/(FIXED_SHARE*F_YEAR[st])          # carbon per dollar of the fixed charge
        vals = [rf(b, e, r)/fs for b, e, r, fs in itertools.product((0.25, 0.5, 0.75), (0.8, 1.0, 1.2), (0.3, 0.4, 0.5), (0.5, 0.6, 0.7))]
        say(f"  facility {st:<38} {central:.3f}  ({min(vals):.3f}-{max(vals):.3f})")

    say("\nROUNDING the DC count to the best integer (all rows of the design table, cost and carbon objectives)")
    worst = 0.0
    for nm_, kw in PRODUCTS:
        for lk in ("fuel", "labor"):
            for rL in (0.50, 0.16):
                q = network(rL=rL, lane_kind=lk, **kw)
                for tau, f in ((0.0, K), (np.inf, C)):
                    n, v = opt(q, tau); best = f(n, v, q)
                    bv = min(f(k, vopt(k, q, tau), q) for k in (int(np.floor(n)), int(np.ceil(n))) if 1 <= k <= q["N"])
                    worst = max(worst, bv/best - 1)
    say(f"  largest penalty: {100*worst:.3f}%")

    say("\nSPOILAGE BREAK-EVEN on the labor lane (refrigerated dairy, diesel delivery): smallest weekly spoilage at which")
    say("the carbon-minimizing chain replenishes faster than the cost-minimizing chain")
    lo, hi = 0.0, 0.01
    for _ in range(30):
        mid = 0.5*(lo+hi)
        p = network(storage="refrigerated", value=1250.0, spoil_week=mid, embodied=0.803, lane_kind="labor")
        nK, vK = opt(p, 0.0); nC, vC = opt(p, np.inf)
        if vC > vK + 0.5: hi = mid
        else: lo = mid
    say(f"  break-even spoilage: {100*hi:.2f}% of stock value per week")

    say("\nCARBON SHARE OF PRICED COST at central ratios, s = tau r/(1 + tau r); and the price at which s = 1/2")
    pa = network(**PRODUCTS[0][1])
    rows = [("holding, ambient packaged goods", pa["eta"]/pa["h"]), ("facilities, ambient DC", pa["phi"]/(FIXED_SHARE*pa["F"])),
            ("last mile, diesel local truck", 0.50), ("last mile, diesel van", 0.30),
            ("line haul, truck instead of intermodal", lane("fuel")["eps0"]/lane("fuel")["t1"]*TL_KG_SMARTWAY/1.339),
            ("line haul, faster vessel speed", 3114/580)]
    for nm_, r in rows:
        sh = [100*t*r/(1+t*r) for t in (0.05, 0.10, 0.20)]
        say(f"  {nm_:<40}{sh[0]:7.1f}%{sh[1]:7.1f}%{sh[2]:7.1f}%  s=1/2 at ${1000/r:,.0f}/t")

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "application_v21_output.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")
