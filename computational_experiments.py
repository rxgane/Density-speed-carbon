"""Numerical illustrations: archetype endpoint premiums, chance x demand-CV,
robustness (integer n, discrete modes, rho, load factor), and a Monte Carlo
diagnostic of the chance constraint."""
import numpy as np
from statistics import NormalDist
from dataclasses import replace
from numerics import (Params, minimize_obj, minimize_carbon, minimize_cost_with_cap,
                      cost_K, carbon_C, SS, d_last)

p = Params()

def design(pp, carbon=False, integer_n=False, vmodes=None):
    ng = np.arange(1, pp.N + 1, dtype=float) if integer_n else np.linspace(1, pp.N, 1961)
    vg = np.array(vmodes, float) if vmodes is not None else np.linspace(pp.vmin, pp.vmax, 1201)
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    obj = carbon_C(NN, VV, pp) if carbon else cost_K(NN, VV, pp)
    i, j = np.unravel_index(np.argmin(obj), obj.shape)
    return ng[i], vg[j]

print("="*70, "\n1. ARCHETYPE PREMIUMS (cost at 80% target vs at carbon-min endpoint)")
from numerics import build_archetypes
for name, pp in build_archetypes(p).items():
    nK, vK, KK, CK = minimize_obj(pp, 0.0)
    nE, vE, KE, CE = minimize_carbon(pp)
    prem_end = cost_K(nE, vE, pp) / KK - 1
    maxcut = 1 - CE / CK
    out = minimize_cost_with_cap(pp, (1 - 0.8 * maxcut) * CK, free="both")
    prem80 = out[2] / KK - 1
    print(f"  {name[:34]:34} maxcut={100*maxcut:5.1f}%  prem@80%={100*prem80:5.2f}%  "
          f"prem@carbonmin={100*prem_end:5.2f}%  (nE={nE:.1f},vE={vE:.2f})")

print("="*70, "\n2. CHANCE x DEMAND CV (carbon-min design at 95% confidence)")
for cv in [0.20, 0.40, 0.60]:
    pcv = replace(p, sigma=cv * p.mu)
    CVQ = cv * np.sqrt((1 + (p.N - 1) * p.rho) / p.N)
    nK, vK, _, _ = minimize_obj(pcv, 0.0)
    m95 = 1 + NormalDist().inv_cdf(0.95) * CVQ
    pm = replace(pcv, eps0=m95 * pcv.eps0, epsL=m95 * pcv.epsL)
    nC, vC, _, _ = minimize_carbon(pm)
    dn = "disperse" if nC > nK + 0.01 else "consolidate"
    print(f"  CV={cv:.2f}  CV_Q={CVQ:.3f}  cost-opt(n={nK:.2f},v={vK:.2f})  "
          f"m@95%={m95:.3f}  carbon-min@95%(n={nC:.2f},v={vC:.2f})  {dn}")

print("="*70, "\n3. ROBUSTNESS")
nK, vK, _, _ = minimize_obj(p, 0.0); nE, vE, _, _ = minimize_carbon(p)
print(f"  continuous : cost(n={nK:.2f},v={vK:.2f})  carbon(n={nE:.2f},v={vE:.2f})")
nKi, vKi = design(p, carbon=False, integer_n=True); nEi, vEi = design(p, carbon=True, integer_n=True)
print(f"  integer n  : cost(n={nKi:.0f},v={vKi:.2f})  carbon(n={nEi:.0f},v={vEi:.2f})")
for modes in [[0.5, 1.0, 2.0]]:
    nKm, vKm = design(p, carbon=False, vmodes=modes); nEm, vEm = design(p, carbon=True, vmodes=modes)
    print(f"  modes {modes}: cost(n={nKm:.2f},v={vKm:.1f})  carbon(n={nEm:.2f},v={vEm:.1f})")
print("  rho sweep (cost-min n, carbon-min n):")
for rho in [0.0, 0.30, 0.60, 0.90]:
    pr = replace(p, rho=rho)
    nKr, _, _, _ = minimize_obj(pr, 0.0); nEr, _, _, _ = minimize_carbon(pr)
    print(f"    rho={rho:.2f}: cost n={nKr:.2f}, carbon n={nEr:.2f}")
print("  load factor on last mile (nref=8, 50% load-sensitive at reference density):")
def carbonmin_lf(pp, nref=8.0):
    ng = np.linspace(1, pp.N, 1961); vg = np.linspace(pp.vmin, pp.vmax, 1201)
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    gfac = 0.5 + 0.5 * (NN / nref)
    C = (pp.phi*NN + pp.eta*SS(NN,VV,pp) + (pp.eps0*VV**2)*pp.N*pp.mu*pp.D
         + gfac*pp.epsL*pp.N*pp.mu*d_last(NN,pp))
    i, j = np.unravel_index(np.argmin(C), C.shape)
    return ng[i], vg[j]
def costmin_lf(pp, nref=8.0):
    # common-multiplier case: the same fill multiplier applies to the last-mile COST rate
    ng = np.linspace(1, pp.N, 1961); vg = np.linspace(pp.vmin, pp.vmax, 1201)
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    gfac = 0.5 + 0.5 * (NN / nref)
    K = (pp.F*NN + pp.h*SS(NN,VV,pp) + (pp.t0 + pp.t1*VV**2)*pp.N*pp.mu*pp.D
         + gfac*pp.tL*pp.N*pp.mu*d_last(NN,pp))
    i, j = np.unravel_index(np.argmin(K), K.shape)
    return ng[i], vg[j]
for nm, pp in [("baseline", p), ("bulky", build_archetypes(p)["Bulky low-value (cheap storage, diesel LM)"])]:
    nK0, _, _, _ = minimize_obj(pp, 0.0); nE0, _, _, _ = minimize_carbon(pp)
    nLF, vLF = carbonmin_lf(pp); nKLF, vKLF = costmin_lf(pp)
    d_co = "disperse" if nLF > nK0 + 0.01 else "consolidate"
    d_sym = "disperse" if nLF > nKLF + 0.01 else "consolidate"
    print(f"    {nm:9}: constant-rate nK={nK0:.2f} nC={nE0:.2f} | carbon-only ell: nC={nLF:.2f} vs nK={nK0:.2f} ({d_co})"
          f" | common ell: nC={nLF:.2f} vs nK={nKLF:.2f} ({d_sym})")

print("="*70, "\n4. MONTE CARLO DIAGNOSTIC (design=cost optimum; cap at normal det.-equiv.)")
rng = np.random.default_rng(7)
M = 500_000
for cv, lbl in [(0.40, "baseline"), (0.80, "high-CV")]:
    pcv = replace(p, sigma=cv * p.mu)
    muQ = pcv.N * pcv.mu
    sigmaQ = pcv.sigma * np.sqrt(pcv.N * (1 + (pcv.N - 1) * pcv.rho))
    Qn = rng.normal(muQ, sigmaQ, M)
    k = muQ**2 / sigmaQ**2; theta = sigmaQ**2 / muQ
    Qg = rng.gamma(k, theta, M)
    print(f"  CV={cv:.2f} (CV_Q={sigmaQ/muQ:.3f}):")
    for conf in [0.90, 0.95, 0.99]:
        z = NormalDist().inv_cdf(conf); qstar = muQ + z * sigmaQ
        vn, vg = np.mean(Qn > qstar), np.mean(Qg > qstar)
        print(f"    nominal viol={1-conf:.3f}  normal-demand={vn:.4f}  gamma-demand={vg:.4f}")

print("="*70, "\n5. ACCOUNTING-BOUNDARY CHECK AND CAP HORIZON")
def transport_carbon(n, v, pp):
    return (pp.eps0*v**2)*pp.N*pp.mu*pp.D + pp.epsL*pp.N*pp.mu*d_last(n, pp)
ngT = np.linspace(1, p.N, 1961); vgT = np.linspace(p.vmin, p.vmax, 1201)
NNT, VVT = np.meshgrid(ngT, vgT, indexing="ij")
nK5, vK5, KK5, CK5 = minimize_obj(p, 0.0); CTK = transport_carbon(nK5, vK5, p)
for tau in [5.0]:
    J = cost_K(NNT, VVT, p) + tau*transport_carbon(NNT, VVT, p)
    i, j = np.unravel_index(np.argmin(J), J.shape); nT, vT = ngT[i], vgT[j]
    Ctot = carbon_C(nT, vT, p); Ct = transport_carbon(nT, vT, p)
    print(f"  tau={tau:.0f}: transport-only opt (n={nT:.2f},v={vT:.2f}): transport {100*(Ct/CTK-1):+.1f}%, "
          f"holding+facility {100*((Ctot-Ct)/(CK5-CTK)-1):+.1f}%, TOTAL {100*(Ctot/CK5-1):+.1f}%")
CVQh = (p.sigma/p.mu)*np.sqrt((1+(p.N-1)*p.rho)/p.N)
z95h = NormalDist().inv_cdf(0.95)
print(f"  cap horizon: m@95 weekly={1+z95h*CVQh:.3f};  annual cap, T=52 independent weeks: m@95={1+z95h*CVQh/np.sqrt(52):.3f}")
