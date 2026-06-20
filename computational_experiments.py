"""Numerical illustrations: archetype endpoint premiums, chance x demand-CV,
robustness (integer n, discrete modes, rho, load factor), and a Monte-Carlo
check of the chance constraint."""
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
print("  load factor on last mile (carbon-min n; nref=8, 50% load-sensitive):")
def carbonmin_lf(pp, nref=8.0):
    ng = np.linspace(1, pp.N, 1961); vg = np.linspace(pp.vmin, pp.vmax, 1201)
    NN, VV = np.meshgrid(ng, vg, indexing="ij")
    gfac = 0.5 + 0.5 * (NN / nref)
    C = (pp.phi*NN + pp.eta*SS(NN,VV,pp) + (pp.eps0*VV**2)*pp.N*pp.mu*pp.D
         + gfac*pp.epsL*pp.N*pp.mu*d_last(NN,pp))
    i, j = np.unravel_index(np.argmin(C), C.shape)
    return ng[i], vg[j]
for nm, pp in [("baseline", p), ("bulky", build_archetypes(p)["Bulky low-value (cheap storage, diesel LM)"])]:
    nE0, _, _, _ = minimize_carbon(pp); nLF, vLF = carbonmin_lf(pp)
    print(f"    {nm:9}: carbon-min n  no-LF={nE0:.2f}  with-LF={nLF:.2f}")

print("="*70, "\n4. MONTE-CARLO CHECK (design=cost optimum; cap at normal det.-equiv.)")
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
