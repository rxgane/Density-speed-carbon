"""Figure 2 data: annual cost and carbon of the reference network as each decision moves from the cost minimum,
with the exact largest cut of each action set for Figure 3. Writes levers_data.json.
Annual units: cost in $ million, carbon in metric tons CO2e."""
import sys, os, json
REV = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REV)
import numpy as np
from scipy.optimize import minimize_scalar
import application_v21 as A, application_v22 as B, targets_v22 as T

Y = 365.0
p = B.reference()
tput = p["N"]*p["mu"]; Dm = p["D"]

def parts(n, v):
    S = A.SS(n, v, p); dn = A.d(n, p)
    cost = dict(facilities=p["F"]*n + B.handling(n, p), safety_stock=p["h"]*S, stock_in_transit=p["pipe_h"]*tput*Dm/v,
                inbound_freight=(p["t0"]+p["t1"]*v**2)*tput*Dm, local_delivery=p["tL"]*tput*dn)
    carbon = dict(facilities=p["phi"]*n, safety_stock=p["eta"]*S, inbound_freight=(p["e0"]+p["eps0"]*v**2)*tput*Dm,
                  local_delivery=p["rL"]*p["tL"]*tput*dn)
    cost = {k: v_*Y/1e6 for k, v_ in cost.items()}; carbon = {k: v_*Y/1e3 for k, v_ in carbon.items()}
    K, C = B.KC(n, v, p)
    assert abs(sum(cost.values()) - K*Y/1e6) < 1e-6 and abs(sum(carbon.values()) - C*Y/1e3) < 1e-6
    return dict(n=float(n), days=float(Dm/v), cost=cost, carbon=carbon, K=K*Y/1e6, C=C*Y/1e3, safety_stock_pallets=float(S),
                transit_pallets=float(tput*Dm/v))

nK, vK = T.solve(p)
nC = float(T.cmin(p)); vC = minimize_scalar(lambda v: B.KC(nC, v, p)[1], bounds=(p["vmin"], p["vmax"]), method="bounded", options={"xatol": 1e-9}).x
vCK = minimize_scalar(lambda v: B.KC(nK, v, p)[1], bounds=(p["vmin"], p["vmax"]), method="bounded", options={"xatol": 1e-9}).x  # carbon-min speed at nK
nCK = minimize_scalar(lambda n: B.KC(n, vK, p)[1], bounds=(1, p["N"]), method="bounded", options={"xatol": 1e-9}).x  # carbon-min n at vK

NS = np.round(np.unique(np.concatenate([np.linspace(1, 30, 291), [nK, nC, nCK]])), 6)
VS = np.unique(np.concatenate([np.linspace(p["vmin"], p["vmax"], 201), [vK, vC, vCK]]))
out = dict(
    notes=dict(units="cost: $ million per year; carbon: t CO2e per year; days = line-haul transit D/v (2 = truckload, 4 = intermodal rail)",
               decisions="n = number of distribution centers (continuous, as in the model); days = line-haul transit time, the speed decision",
               components="cost: facilities, safety_stock (holding h*S), stock_in_transit (financing h_p*P), inbound_freight t(v)*N*mu*D, local_delivery t_L*N*mu*d(n); carbon: same without stock_in_transit",
               slices="vary_n_fixed_v: n moves, line haul fixed at the cost-minimizing speed; vary_v_fixed_n: speed moves, n fixed at the cost minimum; vary_n_best_v / vary_v_best_n: the other decision re-optimized for COST at each value",
               bounds="speed between 250 (rail) and 500 (truckload) miles a day, so days in [2, 4]; n from 1 to 30 shown (model allows up to N=40)"),
    costmin=parts(nK, vK), carbonmin=parts(nC, vC),
    carbonmin_n_at_costmin_speed=parts(nCK, vK), carbonmin_speed_at_costmin_n=parts(nK, vCK),
    vary_n_fixed_v=[parts(n, vK) for n in NS],
    vary_v_fixed_n=[parts(nK, v) for v in VS],
    vary_n_best_v=[parts(n, T.best_v(n, p)) for n in NS],
    vary_v_best_n=[parts(T.best_n(v, p), v) for v in VS],
    grid=dict(n=np.linspace(1, 25, 97).tolist(), days=(Dm/np.linspace(p["vmin"], p["vmax"], 41)).tolist(),
              K=[[B.KC(n, v, p)[0]*Y/1e6 for n in np.linspace(1, 25, 97)] for v in np.linspace(p["vmin"], p["vmax"], 41)],
              C=[[B.KC(n, v, p)[1]*Y/1e3 for n in np.linspace(1, 25, 97)] for v in np.linspace(p["vmin"], p["vmax"], 41)]),
)
R = json.load(open(os.path.join(REV, "targets_v22.json")))
out["least_cost_path_both"] = [dict(cut=r["x"], n=r["both"]["n"], days=r["both"]["days"]) for r in R["rows"] if r["both"] is not None]
json.dump(out, open(os.path.join(REV, "levers_data.json"), "w"))
for k in ("costmin", "carbonmin", "carbonmin_n_at_costmin_speed", "carbonmin_speed_at_costmin_n"):
    d = out[k]; print(f"{k}: n {d['n']:.2f} days {d['days']:.2f} K ${d['K']:.2f}M C {d['C']:,.0f} t | cost " +
                      ", ".join(f"{a} {b:.1f}" for a, b in d["cost"].items()) + " | carbon " + ", ".join(f"{a} {b:,.0f}" for a, b in d["carbon"].items()))
