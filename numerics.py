"""
Numerical illustration for
   "Cost, Carbon, and the Speed of Replenishment in Distribution Network Design"

Everything here evaluates the paper's own closed forms:
  cost            K(n,v)   eq.(K)
  carbon          C(n,v)   eq.(C)
  priced objective J = K + tau*C   eq.(J)
  safety stock    SS(n,v) = z*sigma*sqrt(D/v)*sqrt(N*Psi(n))   eq.(ss)
  pooling         Psi(n)  = N*rho + n*(1-rho)                  eq.(psi)
  last mile       d(n)    = delta * n^{-1/2}                   eq.(dn)
  inbound rates   t(v)=t0+t1 v^2 ,  eps(v)=eps0 v^2            (quadratic instantiation)

The four carbon-intensity-per-dollar ratios (Lemma 1) are the diagnostics:
  speed:    eta/h   vs  eps0/t1   ->  R_speed   = (eta/h)/(eps0/t1)
  density:  epsL/tL vs  phi/F     ->  R_density = (epsL/tL)/(phi/F)

Calibration philosophy: the absolute dollar magnitudes are chosen so the cost-
minimizing network is interior (n* ~ 10 DCs, v* ~ 1); the four kg-per-dollar
ratios are anchored in real-world ballparks (EEIO / GLEC / EPA-style intensities,
in kg CO2e per dollar of activity), since those ratios -- not the absolute
levels -- drive every result.
"""

import numpy as np
from statistics import NormalDist
from dataclasses import dataclass, replace

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

@dataclass
class Params:
    # demand / geography
    N: int = 40            # markets
    mu: float = 20000.0    # mean demand per market per week (units)
    sigma: float = 8000.0  # sd per market per week (CV = 0.40)
    rho: float = 0.30      # cross-market demand correlation
    D: float = 1.0         # mean inbound distance (scaled: lead time L = D/v weeks)
    delta: float = 120.0   # last-mile distance constant (miles): d(n)=delta*n^-0.5
    alpha: float = 0.95    # cycle service level -> z = Phi^{-1}(alpha)
    vmin: float = 0.30
    vmax: float = 3.00
    # ---- cost rates ($) ----
    F: float = 7500.0      # fixed cost per DC per week
    h: float = 0.50        # holding cost $/unit/week
    t0: float = 1.00       # inbound base cost $/unit
    t1: float = 0.05       # inbound speed-cost coefficient ($, on v^2)
    tL: float = 0.005      # last-mile cost $/unit-mile
    # ---- carbon rates (kg CO2e) ----
    phi: float = 3000.0    # facility carbon kg/DC/week
    eta: float = 0.15      # holding carbon kg/unit/week
    eps0: float = 0.04     # inbound speed-carbon coefficient (kg, on v^2)
    epsL: float = 0.0015   # last-mile carbon kg/unit-mile

    @property
    def z(self):
        return NormalDist().inv_cdf(self.alpha)

    # --- the four carbon-intensity-per-dollar ratios ---
    @property
    def r_hold(self):  return self.eta / self.h          # eta/h
    @property
    def r_in(self):    return self.eps0 / self.t1         # eps0/t1
    @property
    def r_last(self):  return self.epsL / self.tL         # epsL/tL
    @property
    def r_fac(self):   return self.phi / self.F           # phi/F
    @property
    def R_speed(self):   return self.r_hold / self.r_in   # >1 => carbon ships faster
    @property
    def R_density(self): return self.r_last / self.r_fac  # >1 => carbon disperses (rho=1 limit)


# --------------------------------------------------------------------------
# Model primitives
# --------------------------------------------------------------------------

def Psi(n, p):           return p.N * p.rho + n * (1.0 - p.rho)
def d_last(n, p):        return p.delta * n ** (-0.5)
def SS(n, v, p):         return p.z * p.sigma * np.sqrt(p.D / v) * np.sqrt(p.N * Psi(n, p))

def cost_K(n, v, p):
    fixed   = p.F * n
    holding = p.h * SS(n, v, p)
    inbound = (p.t0 + p.t1 * v**2) * p.N * p.mu * p.D
    lastm   = p.tL * p.N * p.mu * d_last(n, p)
    return fixed + holding + inbound + lastm

def carbon_C(n, v, p):
    fixed   = p.phi * n
    holding = p.eta * SS(n, v, p)
    inbound = (p.eps0 * v**2) * p.N * p.mu * p.D
    lastm   = p.epsL * p.N * p.mu * d_last(n, p)
    return fixed + holding + inbound + lastm

def buckets(n, v, p):
    """Return (transport_kg, holding_kg) carbon split of eq.(C)."""
    transport = (p.eps0 * v**2) * p.N * p.mu * p.D + p.epsL * p.N * p.mu * d_last(n, p)
    holding   = p.phi * n + p.eta * SS(n, v, p)
    return transport, holding

def J(n, v, p, tau):
    return cost_K(n, v, p) + tau * carbon_C(n, v, p)


# --------------------------------------------------------------------------
# Optimizers (grid search over both n and v; closed-form v used as a starting point)
# --------------------------------------------------------------------------

def _v_grid(p, m=1201):
    return np.linspace(p.vmin, p.vmax, m)

def _n_grid(p, m=1961):
    return np.linspace(1.0, p.N, m)

def minimize_obj(p, tau):
    """Minimize J = K + tau*C over the (n,v) box. Returns (n*, v*, K, C)."""
    ng = _n_grid(p); vg = _v_grid(p)
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    Jv = J(NN, VV, p, tau)
    i, j = np.unravel_index(np.argmin(Jv), Jv.shape)
    n, v = ng[i], vg[j]
    return n, v, cost_K(n, v, p), carbon_C(n, v, p)

def minimize_carbon(p):
    """Pure carbon minimizer (tau -> infinity)."""
    ng = _n_grid(p); vg = _v_grid(p)
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    Cv = carbon_C(NN, VV, p)
    i, j = np.unravel_index(np.argmin(Cv), Cv.shape)
    n, v = ng[i], vg[j]
    return n, v, cost_K(n, v, p), carbon_C(n, v, p)

def minimize_cost_with_cap(p, B, free="both", n_fix=None, v_fix=None):
    """Minimize K s.t. C <= B. `free` in {both, speed, density}."""
    ng = _n_grid(p); vg = _v_grid(p)
    if free == "speed":      ng = np.array([n_fix])
    elif free == "density":  vg = np.array([v_fix])
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    Kv = cost_K(NN, VV, p); Cv = carbon_C(NN, VV, p)
    feas = Cv <= B + 1e-6
    if not feas.any():
        return None
    Kv = np.where(feas, Kv, np.inf)
    i, j = np.unravel_index(np.argmin(Kv), Kv.shape)
    n, v = ng[i], vg[j]
    return n, v, cost_K(n, v, p), carbon_C(n, v, p)

def carbon_floor(p, free="both", n_fix=None, v_fix=None):
    """Minimum achievable carbon under a policy. Returns (n, v, C)."""
    ng = _n_grid(p); vg = _v_grid(p)
    if free == "speed":      ng = np.array([n_fix])    # only v free
    elif free == "density":  vg = np.array([v_fix])    # only n free
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    Cv = carbon_C(NN, VV, p)
    i, j = np.unravel_index(np.argmin(Cv), Cv.shape)
    return ng[i], vg[j], Cv[i, j]


# --------------------------------------------------------------------------
# Illustration A: baseline calibration
# --------------------------------------------------------------------------

def expA_baseline(p):
    print("=" * 74)
    print("A. BASELINE CALIBRATION  (representative product)")
    print("=" * 74)
    print(f"  N={p.N}, mu={p.mu:.0f}, sigma={p.sigma:.0f} (CV={p.sigma/p.mu:.2f}), "
          f"rho={p.rho}, alpha={p.alpha} (z={p.z:.3f})")
    print(f"  delta={p.delta} mi, D={p.D}, v in [{p.vmin},{p.vmax}]")
    print(f"  cost rates:   F={p.F}  h={p.h}  t0={p.t0}  t1={p.t1}  tL={p.tL}")
    print(f"  carbon rates: phi={p.phi}  eta={p.eta}  eps0={p.eps0}  epsL={p.epsL}")
    print(f"  carbon intensity per $ :  eta/h={p.r_hold:.2f}  eps0/t1={p.r_in:.2f}  "
          f"epsL/tL={p.r_last:.2f}  phi/F={p.r_fac:.2f}")
    print(f"  diagnostics:  R_speed={p.R_speed:.2f}   R_density={p.R_density:.2f}")
    nC, vC, KC, CC = minimize_obj(p, 0.0)        # cost min
    nE, vE, KE, CE = minimize_carbon(p)          # carbon min
    tA, hA = buckets(nC, vC, p)
    print("-" * 74)
    print(f"  cost-min   : n*={nC:5.2f}  v*={vC:5.2f}   K=${KC:,.0f}   C={CC:,.0f} kg")
    print(f"               carbon split: transport={tA:,.0f} kg  holding={hA:,.0f} kg")
    print(f"  carbon-min : n*={nE:5.2f}  v*={vE:5.2f}   K=${KE:,.0f}   C={CE:,.0f} kg")
    dirn_v = "FASTER" if vE > vC + 1e-3 else ("slower" if vE < vC - 1e-3 else "same")
    dirn_n = "DISPERSE" if nE > nC + 1e-2 else ("consolidate" if nE < nC - 1e-2 else "same")
    print(f"  direction carbon-min vs cost-min :  speed {dirn_v},  density {dirn_n}")
    print(f"  cost of carbon-min config : +{100*(KE/KC-1):.1f}% over cost-min  "
          f"for {100*(1-CE/CC):.1f}% less carbon")
    return dict(nC=nC, vC=vC, KC=KC, CC=CC, nE=nE, vE=vE, KE=KE, CE=CE)


# --------------------------------------------------------------------------
# Illustration B: four-regime map
# --------------------------------------------------------------------------

def _regime_map(p, rho):
    pp0 = replace(p, rho=rho)
    nC, vC, _, _ = minimize_obj(pp0, 0.0)   # cost-min invariant to carbon rates
    Rs_vals = [0.3, 0.6, 0.9, 1.0, 1.1, 1.5, 2.5]
    Rd_vals = [0.3, 0.6, 0.9, 1.0, 1.1, 1.5, 2.5]
    print(f"\n  rho = {rho}   cost-min reference: n={nC:.2f}, v={vC:.2f}")
    print("  rows = R_speed (eta/h), cols = R_density (epsL/tL)")
    header = "  R_s \\ R_d |" + "".join(f"{rd:>7.1f}" for rd in Rd_vals)
    print(header); print("  " + "-" * (len(header) - 2))
    for rs in Rs_vals:
        cells = []
        for rd in Rd_vals:
            pp = replace(pp0, eta=rs * pp0.r_in * pp0.h,     # eta/h = rs*(eps0/t1)
                              epsL=rd * pp0.r_fac * pp0.tL)  # epsL/tL = rd*(phi/F)
            nE, vE, _, _ = minimize_carbon(pp)
            s = "F" if vE > vC + 1e-3 else ("S" if vE < vC - 1e-3 else "=")
            d = "D" if nE > nC + 1e-2 else ("C" if nE < nC - 1e-2 else "=")
            cells.append(f"{s}{d}")
        print(f"  {rs:>7.1f}  |" + "".join(f"{c:>7}" for c in cells))

def expB_regime_map(p):
    print("\n" + "=" * 74)
    print("B. FOUR-REGIME MAP  (carbon-min direction vs cost-min, joint optima)")
    print("=" * 74)
    print("  Sweep R_speed via eta/h and R_density via epsL/tL; classify the")
    print("  joint carbon-min (n,v) against the joint cost-min (n,v).")
    print("  legend: 1st letter speed {F=faster, S=slower}; "
          "2nd density {D=disperse, C=consolidate}")
    print("\n  -- No pooling (rho=1): thresholds separate, R_speed=1 and R_density=1 --")
    _regime_map(p, 1.0 - 1e-9)
    print("\n  -- With pooling (rho=%.2f): density boundary tilts (pooling-adjusted threshold; holding is the inward force) --"
          % p.rho)
    _regime_map(p, p.rho)


# --------------------------------------------------------------------------
# Illustration C: value of joint redesign
# --------------------------------------------------------------------------

def build_archetypes(p):
    """Four products defined by their carbon-intensity-per-dollar ratios,
    anchored in real-world ballparks; returns {name: Params}."""
    def make(r_hold, r_in, r_last, r_fac):
        return replace(p, eta=r_hold * p.h, eps0=r_in * p.t1,
                          epsL=r_last * p.tL, phi=r_fac * p.F)
    return {
        "Cold-chain pharma (refrig., diesel last mile)": make(0.95, 0.50, 1.40, 0.35),
        "Ambient dry goods (stable, rail/efficient LM)": make(0.30, 0.80, 0.25, 0.40),
        "Bulky low-value (cheap storage, diesel LM)":    make(0.20, 0.45, 1.30, 0.30),
        "E-grocery perishables (refrig., van LM)":       make(0.85, 0.55, 1.50, 0.30),
    }


def expC_value_joint(p, name):
    print("\n" + "=" * 74)
    print("C. VALUE OF JOINT REDESIGN  (single-lever reach vs. joint, then common target)")
    print("=" * 74)
    p = build_archetypes(p)[name]
    print(f"  product: {name}")
    print(f"  R_speed={p.R_speed:.2f}, R_density={p.R_density:.2f} "
          f"(carbon-min wants faster + disperse -> levers reinforce)")
    nC, vC, KC, CC = minimize_obj(p, 0.0)
    print(f"  cost-min: n={nC:.2f}, v={vC:.2f}, K=${KC:,.0f}, C={CC:,.0f} kg")
    # maximum carbon reduction each policy can reach
    policies = [("joint (n,v)",            dict(free="both")),
                ("speed only (n fixed)",   dict(free="speed",   n_fix=nC)),
                ("density only (v fixed)", dict(free="density", v_fix=vC))]
    floors = {}
    print("\n  Maximum achievable carbon reduction by policy:")
    for label, kw in policies:
        _, _, Cf = carbon_floor(p, **kw)
        cut = 1 - Cf / CC
        floors[label] = (Cf, cut)
        print(f"    {label:24s}: floor C={Cf:,.0f} kg   max cut = {100*cut:.2f}%")
    best_single = max(floors["speed only (n fixed)"][1],
                      floors["density only (v fixed)"][1])
    joint_cut = floors["joint (n,v)"][1]
    # moderate target: 90% of the best single-lever reach (singles can still reach it)
    # deep target: 95% of the joint reach (beyond either single lever)
    for tag, target_cut in [("moderate (90% of best single lever)", 0.90 * best_single),
                            ("deep (95% of joint reach)",          0.95 * joint_cut)]:
        B = (1 - target_cut) * CC
        print(f"\n  Target {tag}: cut {100*target_cut:.2f}%  (B = {B:,.0f} kg). "
              f"Cost to comply:")
        for label, kw in policies:
            out = minimize_cost_with_cap(p, B, **kw)
            if out is None:
                print(f"    {label:24s}: INFEASIBLE (lever cannot reach this cut)")
                continue
            n, v, K, C = out
            print(f"    {label:24s}: n={n:5.2f} v={v:5.2f}  "
                  f"K=${K:,.0f}  (+{100*(K/KC-1):.2f}%)")
    print("\n  Reading: each single lever has a bounded reach; the deep target lies")
    print("  beyond both single levers and is met only by adjusting speed and density")
    print("  together -- the complementarity result.")


# --------------------------------------------------------------------------
# Illustration D: carbon-price path
# --------------------------------------------------------------------------

def expD_price_path(p):
    print("\n" + "=" * 74)
    print("D. CARBON-PRICE PATH  (sweep tau; trace n*, v*, K, C)")
    print("=" * 74)
    print(f"  R_speed={p.R_speed:.2f}, R_density={p.R_density:.2f}")
    taus = [0, 5, 10, 20, 35, 50, 100, 300, 800, 1e9]
    print(f"  {'tau':>10} {'n*':>7} {'v*':>7} {'K ($)':>14} {'C (kg)':>14}")
    print("  " + "-" * 56)
    for tau in taus:
        n, v, K, C = minimize_obj(p, tau)
        tlab = "inf" if tau >= 1e8 else f"{tau:.0f}"
        print(f"  {tlab:>10} {n:>7.2f} {v:>7.2f} {K:>14,.0f} {C:>14,.0f}")


# --------------------------------------------------------------------------
# Illustration E: product archetypes
# --------------------------------------------------------------------------

def expE_archetypes(p):
    print("\n" + "=" * 74)
    print("E. PRODUCT ARCHETYPES  (carbon-min direction + cost at 80% of own max cut)")
    print("=" * 74)
    for name, pp in build_archetypes(p).items():
        nC, vC, KC, CC = minimize_obj(pp, 0.0)
        nE, vE, KE, CE = minimize_carbon(pp)
        max_cut = 1 - CE / CC
        target_cut = 0.80 * max_cut
        B = (1 - target_cut) * CC
        capout = minimize_cost_with_cap(pp, B, free="both")
        dv = "faster" if vE > vC + 1e-3 else ("slower" if vE < vC - 1e-3 else "same")
        dn = "disperse" if nE > nC + 1e-2 else ("consolidate" if nE < nC - 1e-2 else "same")
        capcost = f"+{100*(capout[2]/KC-1):.2f}%" if capout else "infeasible"
        print(f"\n  {name}")
        print(f"    ratios: eta/h={pp.r_hold:.2f} eps0/t1={pp.r_in:.2f} "
              f"epsL/tL={pp.r_last:.2f} phi/F={pp.r_fac:.2f}  "
              f"-> R_speed={pp.R_speed:.2f}, R_density={pp.R_density:.2f}")
        print(f"    cost-min (n={nC:.1f},v={vC:.2f}) -> carbon-min (n={nE:.1f},v={vE:.2f}): "
              f"{dv}, {dn}   (max cut {100*max_cut:.1f}%)")
        print(f"    cost to cut {100*target_cut:.1f}% carbon (joint): {capcost}")


# --------------------------------------------------------------------------

if __name__ == "__main__":
    p = Params()
    expA_baseline(p)
    expB_regime_map(p)
    expC_value_joint(p, "E-grocery perishables (refrig., van LM)")
    expD_price_path(p)
    expE_archetypes(p)
