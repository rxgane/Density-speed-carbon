# v2 analysis: exact joint design and analytical product classification

Target: Transportation Research Part E. Branch `v2-tre`; v1 (OR Spectrum state) is tag `v1-orspectrum-ready`.

Notation as in `network_design_v2.tex`: $SS(n,v)=k\sqrt{\Psi(n)}\,v^{-1/2}$ with $k=z\sigma\sqrt{ND}$,
$\Psi(n)=N\rho+n(1-\rho)$, $d(n)=\delta n^{-1/2}$, blended rates $F_\tau,H_\tau,T_\tau(v),\ell_\tau$,
$G(n,v)=\partial_n SS=k(1-\rho)\,v^{-1/2}/(2\sqrt{\Psi(n)})$. Intensity ratios
$r_h=\eta/h$, $r_t=\varepsilon_0/t_1$, $r_L=\varepsilon_L/t_L$, $r_F=\varphi/F$;
$R_{\text{speed}}=r_h/r_t$, $R_{\text{density}}=r_L/r_F$.

All results below are checked against independent brute-force optimization in `v2_checks.py`
(output in `v2_checks_output.txt`).

---

## A. Unique joint design (any carbon price, no second-order condition)

**Proposition A.** Let $t,\varepsilon$ be increasing, differentiable and convex. For every $\tau\ge0$ the
priced problem $\min J(n,v;\tau)$ on $[1,N]\times[v_{\min},v_{\max}]$ has a unique minimizer. The reduced
objective $J^{r}(n)=\min_v J(n,v)$ is strictly quasiconvex on $[1,N]$, and so is $J(\cdot,v)$ for every fixed $v$.

*Proof.* For fixed $n$, $J$ is strictly convex in $v$, so the speed response $v^*(n)$ is unique and continuous;
it depends on $n$ only through $\Psi(n)$ and is nondecreasing in $n$. By the envelope theorem
$J^{r\prime}(n)=F_\tau+H_\tau G(n,v^*(n))-\tfrac12\ell_\tau N\mu\delta\,n^{-3/2}$. Multiply by $n^{3/2}>0$:
$$\phi(n)=F_\tau n^{3/2}+H_\tau\,n^{3/2}G(n,v^*(n))-\tfrac12\ell_\tau N\mu\delta .$$
The first term is strictly increasing. For the second, where $v^*$ is interior the speed condition
$H_\tau k\Psi^{1/2}v^{-3/2}=2N\mu D\,T_\tau'(v)$ gives
$n^{3/2}G(n,v^*)=\tfrac{(1-\rho)N\mu D}{H_\tau}\,n^{3/2}\Psi(n)^{-1}\,T_\tau'(v^*)v^*$; here
$n^{3/2}\Psi^{-1}$ has log-derivative $\tfrac32-n(1-\rho)/\Psi\ge\tfrac12>0$ and $T_\tau'(v)v$ is
nondecreasing in $v$ (convexity), so the product is strictly increasing. Where a speed bound binds,
$n^{3/2}G\propto n^{3/2}\Psi^{-1/2}$, also strictly increasing. Hence $\phi$ is strictly increasing and
$J^{r\prime}$ changes sign at most once, from negative to positive. The fixed-speed case is the same argument
with $v$ constant. $\square$

**Consequences for v1 results.**
- The density second-order condition (eq. `densitysoc`) is not needed for uniqueness or for the threshold of
  Proposition 3: $n_C(v)>n_K(v)$ iff $C_n(n_K,v)<0$ by strict quasiconvexity.
- Proposition 2 (service) and Proposition 4 (complementarity) hold as monotone comparative statics by single
  crossing ($J^r_{ns}>0$, $J_{nv}<0$), without $J_{nn}>0$.
- Corollary 1 (best-response iteration, Topkis) and the remark that a global statement "would require
  convexity over the whole feasible box" are superseded: the joint minimizer is unique although $J$ is not
  jointly convex.
- Checked: 960 (draw, bound, price) cases, zero reduced objectives with a second interior local minimum,
  including 284 cases where the fixed-speed second-order condition fails.

## B. Exact joint classification (replaces the numerical regime map)

Quadratic instantiation $t=t_0+t_1v^2$, $\varepsilon=\varepsilon_0v^2$; speeds interior; $n_K$ interior.

**Pooling weight.** $w=\dfrac{hG(n_K,v_K)}{F+hG(n_K,v_K)}=1-\Big(\dfrac{n_K}{n_K^{0}}\Big)^{3/2}$,
$\;n_K^{0}=\Big(\dfrac{t_LN\mu\delta}{2F}\Big)^{2/3}$ (cost-minimizing density without pooling). The share of
the marginal cost of a DC that is lost pooling; observable from the cost optimum alone.

**Theorem B.**
1. *Density.* $n_C>n_K$ iff
   $$\frac{\varepsilon_L}{t_L}\;>\;(1-w)\,\frac{\varphi}{F}\;+\;w\Big(\frac{\eta}{h}\Big)^{4/5}\Big(\frac{\varepsilon_0}{t_1}\Big)^{1/5},$$
   with $n_C=n_K$ at equality. Once speed re-optimizes, the holding intensity in the density threshold is
   replaced by a geometric blend of holding and inbound-transport intensities.
2. *Speed.* $\ln\dfrac{v_C}{v_K}=\tfrac25\ln R_{\text{speed}}+\tfrac15\ln\dfrac{\Psi(n_C)}{\Psi(n_K)}$
   (intensity effect plus density effect); hence $v_C>v_K$ iff
   $R_{\text{speed}}>\big(\Psi(n_K)/\Psi(n_C)\big)^{1/2}$.
3. *Classification.* With $\underline S=\sqrt{\Psi(n_K)/N}<1<\overline S=\sqrt{\Psi(n_K)/\Psi(1)}$:
   - density condition holds and $R_{\text{speed}}\ge1$: faster and more dispersed;
   - density condition fails and $R_{\text{speed}}\le1$: slower and more consolidated;
   - $R_{\text{speed}}\ge\overline S$: faster whatever the density (faster and more consolidated when the
     density condition fails);
   - $R_{\text{speed}}\le\underline S$: slower whatever the density (slower and more dispersed when the
     density condition holds);
   - only in a mixed cell with $\underline S<R_{\text{speed}}<\overline S$ does the speed direction depend on
     the size of the density move, resolved by part 2 with $n_C$ from the one-dimensional carbon condition.
4. *Power-law speed costs* ($v^m$): exponents become $\tfrac{2m}{2m+1},\tfrac{1}{2m+1}$ in part 1 and
   $\tfrac{2}{2m+1},\tfrac1{2m+1}$ in part 2; the speed condition in part 2 is independent of $m$.

*Proof.* Part 1: by Proposition A for carbon, $n_C>n_K$ iff $C^{r\prime}(n_K)<0$, where
$C^{r\prime}(n_K)=\varphi+\eta\,G(n_K,v_C(n_K))-\tfrac12\varepsilon_LN\mu\delta n_K^{-3/2}$. The cost condition
gives $\tfrac12 t_LN\mu\delta n_K^{-3/2}=F+hG_K$. Since $G\propto v^{-1/2}$ and, by eq. `vstar`,
$v_C(n)/v_K(n)=R_{\text{speed}}^{2/5}$ at a common $n$, $\eta G(n_K,v_C(n_K))=hG_K\,r_h^{4/5}r_t^{1/5}$.
Divide by $F+hG_K$. The identity for $w$ follows from the cost condition and $F=\tfrac12t_LN\mu\delta
(n_K^0)^{-3/2}$. Part 2: eq. `vstar` gives $v^*\propto(H_\tau/T_\tau)^{2/5}\Psi^{1/5}$. Part 3: $\Psi$ is
increasing with $\Psi(1)\le\Psi(n_C)\le\Psi(N)=N$. $\square$

**Checked** (300 interior draws, independent 2-D optimization): part 1 300/300; part 2 identity to
$7.5\times10^{-8}$; zero violations of the aligned cells and of both bands; 21/300 draws fall in the mixed
cell between the bands. Pooling-weight identity to $1.2\times10^{-7}$.

**Baseline and archetypes** (shared cost optimum $(7.84,1.03)$, $w=0.313$, bands $0.661$ and $1.174$):

| case | $R_{\text{speed}}$ | $r_L$ | threshold | Theorem B | computed $(n_C,v_C)$ |
|---|---|---|---|---|---|
| Baseline | 0.375 | 0.30 | 0.389 | slower, consolidate | (6.55, 0.69) |
| Ambient dry goods | 0.375 | 0.25 | 0.389 | slower, consolidate | (5.78, 0.69) |
| Cold-chain pharma | 1.900 | 1.40 | 0.502 | faster, disperse | (16.49, 1.42) |
| Bulky low-value | 0.444 | 1.30 | 0.280 | slower (below lower band), disperse | (22.84, 0.82) |
| E-grocery perishables | 1.545 | 1.50 | 0.450 | faster, disperse | (18.87, 1.32) |
| *Candidate:* refrigerated, electric last mile $(0.95,0.50,0.30,0.35)$ | 1.900 | 0.30 | 0.502 | faster (above upper band), consolidate | (5.45, 1.31) |

The v1 archetypes leave the faster-and-consolidated cell empty; the candidate fills it.

## C. No-pooling corollary: closed-form magnitudes

As $\rho\to1$ ($w=0$, $\Psi\equiv N$): $n_C/n_K=R_{\text{density}}^{2/3}$, $v_C/v_K=R_{\text{speed}}^{2/5}$.
Carbon separates into a proximity part $P(n)=\varphi n+\varepsilon_LN\mu\delta n^{-1/2}$ and an availability
part $A(v)=\eta k\sqrt N v^{-1/2}+\varepsilon_0N\mu Dv^2$, and
$$\frac{C(n_K,v_K)}{C_{\min}}=\theta_P\,\frac{x+2x^{-1/2}}{3}+\theta_A\,\frac{4y^{-1/2}+y^2}{5},
\quad x=R_{\text{density}}^{-2/3},\ y=R_{\text{speed}}^{-2/5},$$
with $\theta_P,\theta_A$ the shares of $P,A$ at the carbon minimum. To second order the excess carbon of the
cost-optimal design is $\theta_P(\ln R_{\text{density}})^2/9+\theta_A\,\tfrac{2}{25}(\ln R_{\text{speed}})^2$:
abatement from redesign grows with the squared log of each ratio and is small for products near parity.
Checked exact to $10^{-7}$ at $\rho=1$. With pooling ($\rho=0.3$) it is an approximation (archetypes:
exact 7.7/6.1/17.4/9.6% vs formula 10.0/5.9/14.9/12.2%; same ordering) and is not claimed there.

## D. Carbon-price path

(a) Fixed density: $\dfrac{d\ln v^*}{d\ln\tau}=\tfrac25(s_h-s_t)$, $s_x=\dfrac{\tau r_x}{1+\tau r_x}$, the
carbon share of the blended rate. (b) No pooling: $\dfrac{d\ln n^*}{d\ln\tau}=\tfrac23(s_L-s_F)$.
(c) With pooling: if $r_L>\max(r_F,r_h,r_t)$, $n^*(\tau)$ is nondecreasing along the whole price path; if
$r_L<\min(r_F,r_h,r_t)$, nonincreasing.

*Proof of (c).* Write $H_\tau G(n,v^*(n;\tau))=B_\tau g(n)$ with $B_\tau\propto H_\tau^{4/5}T_\tau^{1/5}$.
Dividing $\phi$ by $\ell_\tau$, its root moves right when $F_\tau/\ell_\tau$ and $B_\tau/\ell_\tau$ fall in
$\tau$. Each log-derivative is a combination of terms $1/(r_x^{-1}+\tau)$, increasing in $r_x$, so both are
nonpositive for all $\tau$ when $r_L$ exceeds the other ratios. $\square$

Checked: (a), (b) against numerical derivatives to $10^{-5}$; (c) 18/18 and 42/42 price paths monotone.

---

## Proposed v2 structure (for discussion; no manuscript edits yet)

- **Section 4.** Keep Lemma 1 and the two marginal thresholds (Propositions 1, 3), stated without
  `densitysoc`. Add Proposition A. Replace Corollary 1 with Theorem B. Keep complementarity and service as
  monotone comparative statics. Corollary C and Proposition D as compact results. Chance constraint unchanged.
- **Section 5.** Archetypes become instances of Theorem B: the table reports the theorem's classification
  (inputs: four ratios, $w$, the two bands) next to the computed magnitudes. Regime figure drawn from the
  closed-form boundaries. Add the faster-and-consolidated archetype so all four cells are populated.
- **Section 6.** Managerial implications rebuilt on Theorem B, Corollary C and Proposition D (outline below).
- **Back matter for TR-E** (Elsevier): highlights and Elsevier declarations return; author-year references.

## Managerial implications: outline (each point rests on a stated result)

1. *What a firm needs.* Four carbon-intensity-per-dollar ratios and one cost share $w$ (the lost-pooling share
   of a DC's marginal cost). Theorem B then gives the carbon-minimizing direction of both decisions.
2. *Speed rule with a margin of safety.* Replenish faster whatever happens to the network when holding exceeds
   inbound freight in carbon per dollar by the factor $\overline S$ (1.17 at the baseline); slower whenever it
   falls below $\underline S$ (0.66). Between the bands the network decision must be settled first.
3. *Inbound freight affects how many DCs to run.* Through speed re-optimization, cleaner inbound freight lowers
   the dispersal threshold. A firm electrifying line haul should revisit network density, not only speed.
4. *Technology changes move products across cells* (signs from Theorem B): cleaner facility energy and cleaner
   inbound freight favor dispersal; an electric last mile favors consolidation; more efficient refrigeration
   lowers both the speed index and the dispersal threshold.
5. *Service, volatility and correlation act through $w$ and the bands.* Higher service targets or volatility
   lower $n_K$ and raise $w$ (Proposition 2), tilting density toward the holding-transport blend. Higher
   correlation raises $n_K$ and lowers $w$: $H_\tau G(n,v^*(n))\propto(1-\rho)\Psi^{-3/5}$ falls in $\rho$
   because $\Psi=N\rho+n(1-\rho)$ rises in $\rho$ for $n<N$, so the root of $\phi$ moves right (single
   crossing, Proposition A). Since $w=1-(n_K/n_K^0)^{3/2}$ and both bands depend on $n_K$, every
   comparative static on $n_K$ moves the classification inputs in a known direction.
6. *Where redesign pays.* Without pooling, excess carbon of the cost-optimal design grows with the squared log
   ratios (Corollary C): products near parity gain little; prioritize products with ratios far from one.
7. *How the design should move with the carbon price.* Speed moves with the gap in carbon shares of holding and
   inbound freight (Proposition D); the network disperses monotonically along the whole price path when the
   last mile is the most carbon-intensive activity per dollar.
8. *Policy.* Transport-only instruments (unchanged from v1), now with the joint theorem behind the substitution.
