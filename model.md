# Network Design with Carbon — Model

Working title: *Cost, Carbon, and Availability in Distribution Network Design:
Joint Choice of Density, Speed, and Safety Stock*

Motivation instance: the Wyllo supply-chain simulator (~/Documents/wyllo_deploy).
Wyllo is motivation + numerical validation only, NOT the contribution.

---

## 1. Setting

Single source `s` → `n` distribution centers (DCs) → `N` markets.
Single product, stationary stochastic demand, periodic review, base-stock.
Symmetric markets and DCs (tractability; relaxed in simulation).

## 2. Parameters (given)

- N         number of markets
- mu, sigma mean, std of demand per market per period
- rho       cross-market demand correlation (equicorrelation, in [0,1])
- D         inbound distance (source -> DC), constant in n
- delta     last-mile geography constant
- F, phi    fixed COST, CARBON per DC per period
- h, eta    holding COST, CARBON per unit per period
- t(v), eps(v)  inbound transport COST, CARBON per unit-mile (increasing, convex in v)
- t_L, eps_L    last-mile transport COST, CARBON per unit-mile (fixed local mode)
- alpha     required service level (fill rate)  -> z = Phi^{-1}(alpha)   [PARAMETER for now]
- tau / B   carbon price / carbon budget

Canonical convex forms (slow-steaming physics): eps(v) = eps_0 * v^2, t(v) = t_0 + t_1 * v^2.

## 3. Decision variables

- n   network density (# DCs),     1 <= n <= N        [continuous in core, rounded]
- v   inbound replenishment speed, v in [v_min, v_max]
- (alpha endogenous only in the "service currency" extension)

## 4. Derived quantities

- g = N / n                  markets per DC (pooling group size)
- L = D / v                  replenishment lead time (faster v -> shorter L)
- d(n) = delta / sqrt(n)     average last-mile distance (Daganzo spatial law)
- Psi(n) = N*rho + n*(1-rho) pooling expression (linear in n)
- z = Phi^{-1}(alpha)        safety factor

Total safety stock (built in model_derivation notes):
  SS(n,v) = z * sigma * sqrt(D/v) * sqrt( N * Psi(n) )
          = z * sigma * sqrt(L)   * sqrt( n*N*(1 + (N/n - 1)*rho) )
  Checks: rho=0,n=1 -> z*sigma*sqrt(L)*sqrt(N);  rho=0,n=N -> z*sigma*sqrt(L)*N;
          rho=1 any n -> z*sigma*sqrt(L)*N (constant).  At rho=0: SS ~ sqrt(n) (square-root law).

## 5. Objective functions (same 4-piece skeleton, $ vs kg rates)

COST:
  K(n,v) =  F*n                              (fixed facility)
          + h   * SS(n,v)                     (holding)
          + t(v)   * N*mu*D                   (inbound transport)
          + t_L * N*mu * delta/sqrt(n)        (last-mile transport)

CARBON:
  C(n,v) =  phi*n
          + eta * SS(n,v)
          + eps(v) * N*mu*D
          + eps_L * N*mu * delta/sqrt(n)

Two carbon buckets:
  TRANSPORT  = eps(v)*N*mu*D + eps_L*N*mu*delta/sqrt(n)   (up in v, down in n)
  HOLDING    = phi*n + eta*SS                              (down in v, up in n)
Both transport-cutting moves (slow down, spread out) inflate the holding bucket.

## 6. The program (three equivalent lenses)

(P_cap)    min_{n,v} K(n,v)  s.t.  C(n,v) <= B,  fill >= alpha,  boxes
(P_price)  min_{n,v} J(n,v;tau) = K + tau*C,  tau >= 0
(P_Pareto) efficient set traced by sweeping tau (or B)
Duality: shadow price of (P_cap) equals tau.

## 7. The workhorse J (blended coefficients)

  J(n,v) = a*n + b(v)*sqrt(Psi(n)) + c(v) + e * n^{-1/2}

  a    = F   + tau*phi
  b(v) = (h  + tau*eta) * z*sigma*sqrt(N) * sqrt(D/v)      = B0 * v^{-1/2},  B0=(h+tau*eta)*z*sigma*sqrt(N*D)
  c(v) = (t(v) + tau*eps(v)) * N*mu*D
  e    = (t_L + tau*eps_L) * N*mu*delta

  tau=0 -> cost optimum;  tau->inf -> carbon optimum;  sweep tau -> Pareto frontier.

## 8. First-order conditions

FOC in n:
  dJ/dn = a + b(v)*(1-rho)/(2*sqrt(Psi(n))) - (e/2)*n^{-3/2} = 0
  Balance:  [marginal cost of adding a DC] = [marginal last-mile saving]
            a + b(v)(1-rho)/(2 sqrt(Psi(n)))  =  (e/2) n^{-3/2}
  Note: holding pull carries factor (1-rho); vanishes at rho=1.
  SOC caveat: holding term makes J concave in n (sqrt of linear); interior optimum
  is a min only if (3e/4)n^{-5/2} > b(1-rho)^2/(4 Psi^{3/2}); else CORNER (n=1 or n=N).

FOC in v:
  dJ/dv = -(b(v)/(2v))*sqrt(Psi(n)) + (t'(v)+tau*eps'(v))*N*mu*D = 0
  Balance:  [marginal transport $ + carbon of speed] = [marginal holding saving]
            (t'(v)+tau*eps'(v))*N*mu*D  =  b(v)*sqrt(Psi(n))/(2v)

## 9. Closed-form optimal speed (canonical quadratic forms)

With eps(v)=eps_0 v^2, t(v)=t_0+t_1 v^2:

  v*(n) = [ (h+tau*eta)*z*sigma*sqrt(N*D)*sqrt(Psi(n)) / ( 4*(t_1+tau*eps_0)*N*mu*D ) ]^{2/5}

Comparative statics (exponent 2/5 > 0):
  v* increases in (h+tau*eta), z(alpha), sigma         -> heavier/volatile/high-service inventory => ship FASTER
  v* increases in sqrt(Psi(n))  (i.e., in n)           -> COMPLEMENTARITY: decentralize => ship faster
  v* decreases in (t_1+tau*eps_0)                      -> dirtier/costlier transport => slower

SPEED THRESHOLD (Theorem 2):
  v*_carbon > v*_cost   <=>   eta/h > eps_0/t_1   (holding "dirtier per dollar" than inbound transport)

## 10. Density threshold (Theorem 3, density part)

rho=1 (no pooling) closed form:  n* = (e / (2a))^{2/3}
DENSITY THRESHOLD:
  n*_carbon > n*_cost  <=>  eps_L/t_L > phi/F   (last-mile "dirtier per dollar" than facilities)
  For rho<1: holding (eta/h) re-enters with a CENTRALIZING pull.
  General rule: carbon moves n toward whichever bucket is dirtiest per dollar:
    facilities (phi/F) -> centralize;  holding (eta/h) -> centralize;  last-mile (eps_L/t_L) -> decentralize.

## 11. Propositions

P1 (interior optima + complementarity). Under SOCs a unique interior (n*,v*) solves the
   coupled FOCs; dv*/dn > 0 and dn*/dv > 0 (levers are complements / J supermodular after sign flip).
P2 (speed flip / Theorem 2). v*_carbon > v*_cost  iff  eta/h > eps_0/t_1. v* up in sigma, alpha, n.
P3 (density flip / Theorem 3). n*_carbon > n*_cost  iff  eps_L/t_L > phi/F (+ centralizing pull from eta/h if rho<1).
P4 (counterintuitive corner). When P2 and P3 both bind, carbon-optimal chain is FASTER AND MORE
   DECENTRALIZED, contradicting "slow + consolidate = green"; complementarity amplifies the shift.

## 12. Managerial map (how carbon reconfigures the network)

Columns set by density threshold (last-mile vs facilities); rows by speed threshold (holding vs inbound).

                 | Centralize (eps_L/t_L < phi/F) | Decentralize (eps_L/t_L > phi/F)
  Ship fast      | Pool & speed up                | FASTER & MORE SPREAD OUT  *** counterintuitive
  (eta/h>eps0/t1)|                                | (cold-chain, volatile demand, dirty last-mile)
  Ship slow      | Slow & consolidate (orthodoxy) | Slow the haul, spread out
  (eta/h<eps0/t1)|                                |

Key managerial insights:
  - No single "green network": depends on dirtiest-per-dollar bucket (the 2x2 is the decision rule).
  - "Consolidate to cut emissions" is half-wrong: centralizing LENGTHENS last-mile (raises transport
    carbon); its real green benefit is inventory POOLING + facility energy.
  - Higher volatility sigma / service alpha -> holding-heavier -> carbon pushes faster + more central.
  - Higher correlation rho -> pooling buys less -> tilts toward decentralization.
  - Higher carbon price tau -> amplifies the gap from cost-optimal design (does not change direction).
  - Levers are complements -> re-optimize network AND speed together; piecemeal carbon programs backfire.

## 13. Standing assumptions

1. Demand Normal, stationary, equicorrelated; symmetric markets/DCs.
2. t(v), eps(v) increasing convex; canonical eps(v)=eps_0 v^2.
3. D = MEAN source-to-DC distance (heterogeneous D_j enter only through their mean once throughput is
   conserved at N*mu); ~constant in n by region geometry, under POINT-TO-POINT replenishment. [If routed,
   inbound ~ sqrt(n) a la Cachon tour -> reinforces centralization in Thm 3, leaves Thm 2 unchanged.]
   Using mean D in the sqrt(L) safety-stock term is a mild Jensen approximation. Last-mile fixed local
   mode; last-mile time not in availability.
4. n relaxed to continuous; cycle stock constant (fixed review), set aside.

## 14. Extensions (deferred / sequel)

- HUB / vertical inventory placement: add a central hub -> genuine Graves-Willems-style safety-stock
  placement (which echelon decouples). Trades pooling (holding carbon) vs responsiveness/expediting.
- LOAD-FACTOR endogenous emission factor (2b): per-unit transport emission e(q)=a/q+b (fuller trucks
  cleaner) -> a new force making centralization cleaner per unit. Grounded (Bektas-Laporte, Bouchery).
- TIME-VARYING / stochastic carbon intensity ("carbon-aware logistics"): exogenous eps(t); needs a
  dynamic timing model -> a separate paper.
- Endogenous alpha -> "service currency" (buy availability with stock vs speed); the carbon-cheaper
  currency flips on the SAME eta/h vs eps/t threshold (Thm 2 via the availability lens).
