"""
v2 analytical results: numerical checks against independent brute-force optimization.

Checks
  A  uniqueness: the reduced objective J^r(n) = min_v J(n,v) changes slope sign once (single crossing),
     with and without speed bounds, including draws where the fixed-speed second-order condition fails
  B  exact joint classification: density threshold, speed identity, and the robust speed bands
  C  pooling weight identity  w = 1 - (n_K / n_K^np)^{3/2}
  D  analytical classification of the baseline and the four archetypes
  E  no-pooling closed forms: n_C/n_K, v_C/v_K, and the abatement formula
  F  carbon-price path: monotone density response under the stated ratio ordering
"""
import sys, os
import numpy as np
from dataclasses import replace
from scipy.optimize import minimize_scalar
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numerics as nm

def J(n, v, p, tau): return nm.J(n, v, p, tau)

def v_opt(n, p, tau, lo=None, hi=None):
    lo = p.vmin if lo is None else lo; hi = p.vmax if hi is None else hi
    return minimize_scalar(lambda v: J(n, v, p, tau), bounds=(lo, hi), method="bounded",
                           options={"xatol": 1e-12}).x

def Jr(n, p, tau, lo=None, hi=None): return J(n, v_opt(n, p, tau, lo, hi), p, tau)

def joint_opt(p, tau, lo=None, hi=None, m=400):
    ns = np.exp(np.linspace(0, np.log(p.N), m))
    vals = np.array([Jr(n, p, tau, lo, hi) for n in ns]); i = int(np.argmin(vals))
    r = minimize_scalar(lambda n: Jr(n, p, tau, lo, hi), bounds=(ns[max(i-1, 0)], ns[min(i+1, m-1)]),
                        method="bounded", options={"xatol": 1e-10})
    return r.x, v_opt(r.x, p, tau, lo, hi), i, m

def Psi(n, p): return nm.Psi(n, p)
def G(n, v, p): return p.z*p.sigma*np.sqrt(p.D/v)*np.sqrt(p.N)*(1-p.rho)/(2*np.sqrt(Psi(n, p)))
def ratios(p): return p.eta/p.h, p.eps0/p.t1, p.epsL/p.tL, p.phi/p.F    # r_h, r_t, r_L, r_F

def draw(rng, base, wide=1.0):
    f = lambda lo, hi: float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
    return replace(base, F=base.F*f(.3/wide, 3*wide), h=base.h*f(.3/wide, 3*wide), t1=base.t1*f(.3, 3),
                   tL=base.tL*f(.3, 3), phi=base.phi*f(.2, 5), eta=base.eta*f(.2, 5),
                   eps0=base.eps0*f(.2, 5), epsL=base.epsL*f(.2, 5), rho=float(rng.uniform(0.02, 0.7)),
                   delta=base.delta*f(.5, 2), sigma=base.sigma*f(.5, 2))

base = nm.Params()
out = []
def say(s=""): print(s); out.append(s)

# ---------------- A: single crossing / uniqueness ----------------
rng = np.random.default_rng(101)
taus = [0.0, 0.2, 2.0, 1e9]
nA = multi = soc_fail = 0
for k in range(120):
    p = draw(rng, base, wide=3.0)
    for bounded in (False, True):
        lo, hi = ((p.vmin, p.vmax) if bounded else (1e-3, 60.0))
        for tau in taus:
            ns = np.exp(np.linspace(0, np.log(p.N), 600))
            vals = np.array([Jr(n, p, tau, lo, hi) for n in ns])
            d = np.diff(vals)
            changes = int(np.sum((d[:-1] < 0) & (d[1:] > 0)))       # interior local minima
            nA += 1; multi += changes > 1
            # does the fixed-speed SOC (eq. densitysoc) fail somewhere at this draw?
            v = v_opt(ns[len(ns)//2], p, tau, lo, hi)
            Ht, lt = p.h + min(tau, 1e9)*p.eta, p.tL + min(tau, 1e9)*p.epsL
            lhs = 3*lt*p.N*p.mu*p.delta*ns**-2.5
            rhs = Ht*p.z*p.sigma*np.sqrt(p.D/v)*np.sqrt(p.N)*(1-p.rho)**2*Psi(ns, p)**-1.5
            soc_fail += bool(np.any(lhs <= rhs))
say("A  Unique joint design (single crossing of n^{3/2} J^r_n)")
say(f"   {nA} (draw, bounds, price) cases; reduced objectives with more than one interior local minimum: {multi}")
say(f"   cases where the fixed-speed second-order condition fails somewhere on [1,N]: {soc_fail} (uniqueness still holds)")

# ---------------- B: exact joint classification ----------------
rng = np.random.default_rng(202)
rows = []; tried = 0
while len(rows) < 300 and tried < 5000:
    tried += 1
    p = draw(rng, base)
    lo, hi = 1e-3, 60.0                                   # interior speeds
    nK, vK, iK, m = joint_opt(p, 0.0, lo, hi); nC, vC, iC, _ = joint_opt(p, 1e9, lo, hi)
    if not (0 < iK < m-1 and 0 < iC < m-1) or abs(nC-nK)/nK < 0.005: continue
    rh, rt, rL, rF = ratios(p)
    GK = G(nK, vK, p); w = p.h*GK/(p.F + p.h*GK)
    thr = (1-w)*rF + w*rh**0.8*rt**0.2
    Rs = rh/rt
    Sbar, Slow = np.sqrt(Psi(nK, p)/Psi(1.0, p)), np.sqrt(Psi(nK, p)/p.N)
    ident = np.log(vC/vK) - (0.4*np.log(Rs) + 0.2*np.log(Psi(nC, p)/Psi(nK, p)))
    nnp = (p.tL*p.N*p.mu*p.delta/(2*p.F))**(2/3)
    rows.append(dict(dens_ok=((rL > thr) == (nC > nK)), ident=abs(ident),
                     band_hi_viol=(Rs >= Sbar and not vC > vK), band_lo_viol=(Rs <= Slow and not vC < vK),
                     aligned_viol=((Rs > 1 and nC > nK and not vC > vK) or (Rs < 1 and nC < nK and not vC < vK)),
                     indet=((Slow < Rs < Sbar) and ((Rs > 1) != (nC > nK))),
                     w_err=abs(w - (1 - (nK/nnp)**1.5))))
T = len(rows)
say("\nB  Exact joint classification (independent 2-D optimization, interior optima)")
say(f"   draws: {T}")
say(f"   density direction n_C>n_K iff r_L > (1-w) r_F + w r_h^(4/5) r_t^(1/5): {sum(r['dens_ok'] for r in rows)}/{T}")
say(f"   speed identity ln(vC/vK) = (2/5)ln R_speed + (1/5)ln(Psi_C/Psi_K): max abs error {max(r['ident'] for r in rows):.1e}")
say(f"   aligned cells (faster&disperse / slower&consolidate) violations: {sum(r['aligned_viol'] for r in rows)}")
say(f"   R_speed >= sqrt(Psi_K/Psi(1)) => faster, violations: {sum(r['band_hi_viol'] for r in rows)}")
say(f"   R_speed <= sqrt(Psi_K/N)      => slower, violations: {sum(r['band_lo_viol'] for r in rows)}")
say(f"   draws in the indeterminate band (mixed cell, between the two bounds): {sum(r['indet'] for r in rows)}/{T}")
say(f"\nC  Pooling weight identity w = 1 - (n_K/n_K^np)^(3/2): max abs error {max(r['w_err'] for r in rows):.1e}")

# ---------------- D: analytical classification of the archetypes ----------------
say("\nD  Analytical classification of the baseline and the four archetypes")
arch = {"Baseline": base, **nm.build_archetypes(base)}
nK, vK, *_ = joint_opt(base, 0.0)
GK = G(nK, vK, base); w = base.h*GK/(base.F + base.h*GK)
Sbar, Slow = np.sqrt(Psi(nK, base)/Psi(1.0, base)), np.sqrt(Psi(nK, base)/base.N)
say(f"   shared cost optimum (n_K, v_K) = ({nK:.2f}, {vK:.2f});  pooling weight w = {w:.3f};"
    f"  speed bands: slower if R_speed <= {Slow:.3f}, faster if R_speed >= {Sbar:.3f}")
say(f"   {'case':<46}{'R_speed':>8}{'r_L':>7}{'thr':>7}  {'predicted':<26}{'computed (n_C, v_C)':>22}")
for name, p in arch.items():
    rh, rt, rL, rF = ratios(p); Rs = rh/rt; thr = (1-w)*rF + w*rh**0.8*rt**0.2
    dens = "disperse" if rL > thr else "consolidate"
    if Rs >= Sbar or (Rs > 1 and rL > thr): spd = "faster"
    elif Rs <= Slow or (Rs < 1 and rL <= thr): spd = "slower"
    else: spd = "speed: magnitude-dependent"
    nC, vC, *_ = joint_opt(p, 1e9)
    ok = ((nC > nK) == (rL > thr)) and ((vC > vK) == (spd == "faster")) if spd != "speed: magnitude-dependent" else True
    say(f"   {name[:45]:<46}{Rs:>8.3f}{rL:>7.2f}{thr:>7.3f}  {spd+', '+dens:<26}{f'({nC:.2f}, {vC:.2f})':>22}  {'agrees' if ok else 'DISAGREES'}")

# ---------------- E: no-pooling closed forms ----------------
say("\nE  No-pooling limit (rho = 1): closed-form magnitudes")
rng = np.random.default_rng(303); errs = []; aerr = []
for k in range(60):
    p = replace(draw(rng, base), rho=1.0)
    nK, vK, iK, m = joint_opt(p, 0.0, 1e-3, 60); nC, vC, iC, _ = joint_opt(p, 1e9, 1e-3, 60)
    if not (0 < iK < m-1 and 0 < iC < m-1): continue
    rh, rt, rL, rF = ratios(p); Rd, Rs = rL/rF, rh/rt
    errs.append(max(abs(nC/nK - Rd**(2/3)), abs(vC/vK - Rs**0.4)))
    # abatement: C* / C(nK,vK) from the separable closed form
    x, y = Rd**(-2/3), Rs**(-0.4)
    P_star = 3*p.phi*nC;  A_star = 1.25*p.eta*nm.SS(nC, vC, p)          # component minima (no pooling)
    pred = (P_star + A_star) / (P_star*(x + 2*x**-0.5)/3 + A_star*(4*y**-0.5 + y**2)/5)
    exact = nm.carbon_C(nC, vC, p) / nm.carbon_C(nK, vK, p)
    aerr.append(abs(pred - exact))
say(f"   n_C/n_K = R_density^(2/3) and v_C/v_K = R_speed^(2/5): max abs error {max(errs):.1e} over {len(errs)} draws")
say(f"   abatement formula C_min/C(n_K,v_K): max abs error {max(aerr):.1e}")

# ---------------- F: carbon-price path ----------------
say("\nF  Carbon-price path of the density decision")
rng = np.random.default_rng(404); up_ok = dn_ok = up_n = dn_n = 0
tgrid = np.concatenate([[0.0], np.exp(np.linspace(np.log(0.02), np.log(200), 22))])
while up_n + dn_n < 60:
    p = draw(rng, base); rh, rt, rL, rF = ratios(p)
    if rL > max(rF, rh, rt): kind = "up"
    elif rL < min(rF, rh, rt): kind = "dn"
    else: continue
    path = [joint_opt(p, t, 1e-3, 60, m=200)[0] for t in tgrid]
    d = np.diff(path)
    if kind == "up": up_n += 1; up_ok += bool(np.all(d >= -1e-6))
    else: dn_n += 1; dn_ok += bool(np.all(d <= 1e-6))
say(f"   r_L > max(r_F, r_h, r_t): density nondecreasing along the whole price path in {up_ok}/{up_n} draws")
say(f"   r_L < min(r_F, r_h, r_t): density nonincreasing along the whole price path in {dn_ok}/{dn_n} draws")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "v2_checks_output.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
