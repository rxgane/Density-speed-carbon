import numpy as np
from dataclasses import replace
from numerics import Params, minimize_obj, minimize_carbon, SS, Psi, d_last, build_archetypes

p = Params()

def cond13(pp, n, v, carbon):
    """3 e n^-5/2 > b(v)(1-rho)^2 Psi^-3/2 ; carbon=True uses carbon coeffs, else cost."""
    e = (pp.epsL if carbon else pp.tL) * pp.N * pp.mu * pp.delta
    bcoef = (pp.eta if carbon else pp.h) * pp.z * pp.sigma * np.sqrt(pp.N) * np.sqrt(pp.D / v)
    lhs = 3 * e * n**(-2.5)
    rhs = bcoef * (1 - pp.rho)**2 * Psi(n, pp)**(-1.5)
    return lhs, rhs, lhs > rhs

def dCdn(pp, n, v):
    eps = 1e-4
    from numerics import carbon_C
    return (carbon_C(n + eps, v, pp) - carbon_C(n - eps, v, pp)) / (2 * eps)

print("=== Condition (13) strict convexity in n  [3e n^-5/2 > b(1-rho)^2 Psi^-3/2] ===")
nK, vK, _, _ = minimize_obj(p, 0.0); nC, vC, _, _ = minimize_carbon(p)
for tag, n, v, carb in [("baseline cost-min", nK, vK, False), ("baseline carbon-min", nC, vC, True)]:
    l, r, ok = cond13(p, n, v, carb)
    print(f"  {tag:22} n={n:5.2f} v={v:4.2f}: LHS={l:,.0f} {'>' if ok else '<='} RHS={r:,.0f}  {'HOLDS' if ok else 'FAILS'}")
for name, pp in build_archetypes(p).items():
    nCa, vCa, _, _ = minimize_carbon(pp)
    l, r, ok = cond13(pp, nCa, vCa, True)
    print(f"  {name[:30]:30} n_C={nCa:5.2f} v_C={vCa:4.2f}: LHS={l:,.0f} {'>' if ok else '<='} RHS={r:,.0f}  {'HOLDS' if ok else 'FAILS'}")

print("\n=== Interior optimum in n  [dC/dn(1) < 0 < dC/dn(N)] at carbon-min v ===")
for name, pp in build_archetypes(p).items():
    nCa, vCa, _, _ = minimize_carbon(pp)
    d1, dN = dCdn(pp, 1.0, vCa), dCdn(pp, float(pp.N), vCa)
    interior = (1.0 < nCa < pp.N) and (d1 < 0 < dN)
    print(f"  {name[:30]:30} n_C={nCa:5.2f}  dC/dn(1)={d1:+,.0f}  dC/dn(N)={dN:+,.0f}  "
          f"{'INTERIOR' if interior else 'BOUNDARY/other'}")
