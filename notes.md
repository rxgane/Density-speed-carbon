# Notes — Positioning & Literature Scan

## Contribution (reweighted after lit scan)

LEAD with the cleanest white space:
1. SPEED REVERSAL (Thm 2): a carbon constraint can make the optimal chain FASTER, not slower,
   because faster shipping kills carbon-dirty safety stock. Threshold eta/h > eps_0/t_1. Overturns
   the slow-steaming orthodoxy via the safety-stock-carbon channel (a NEW channel).
2. JOINT (n,v) COUPLING / complementarity: carbon pushes faster AND more decentralized together.
   The joint closed-form is unoccupied (density papers lack speed; speed papers lack density).
3. DENSITY (Thm 3): NOT a standalone contribution. After reading Cachon (2014) in full, his model
   ALREADY has: U-shaped total cost/carbon in n (phi_c*c_c*n^-1/2 down vs (phi_t*c_t+phi_s*c_s)*n^1/2 up),
   sqrt(n) RISK POOLING in the inventory/space term (base-stock z*sigma_j, sigma_j ~ sqrt(1/n)),
   closed-form n* = phi_c*c_c/(phi_t*c_t+phi_s*c_s), emission-min vs cost-min density, and the
   "carbon constraint can raise emissions" direction. [The earlier scan agent's claim that Cachon
   lacks pooling was WRONG.] Our density result is essentially Cachon RELABELED for b2b firm
   distribution (firm last-mile instead of consumer driving). Use it as a BUILDING BLOCK / replication,
   not a claimed contribution.
4. SERVICE CURRENCY: buy availability with stock or speed; carbon-cheaper currency flips on the
   SAME eta/h vs eps/t threshold. Use fill rate (type-2), NOT service-time (that's GSM, occupied).

DO NOT claim: "first to jointly optimize" (Li et al. 2017 did it numerically); "first 3-way
cost/carbon/service" (Wang et al. 2025 has cost/carbon/service-time); "carbon can decentralize is new"
(Cachon 2014, Turken 2017, Harris 2011, McKinnon all have the direction).

Methodological hook: carbon-aware CONTINUOUS-APPROXIMATION (Daganzo) network design WITH the
inventory-pooling tradeoff is genuinely open (CA+carbon exists only for routing w/ constant EF).
Lean on Benjaafar-Li-Daskin "insights from simple/stylized models" license.

## Lit scan verdict table

| Result                         | Status                              | Closest threat                          | Differentiator                                  |
|--------------------------------|-------------------------------------|-----------------------------------------|-------------------------------------------------|
| Joint (n,v)+inventory closed-form | OPEN (done numerically)          | Li et al. 2017 (PLoS, MINLP/LINGO)      | closed-form/structural; continuous speed        |
| Thm 2 speed reversal           | OPEN (strongest)                    | Schaefer-Konur 2015 (TR-E); Hoen 2014   | continuous speed, eps~v^2, U-shape, eta/h>eps/t |
| Thm 3 density (decentralize)   | DIRECTION known; mechanism open     | Cachon 2014 (MnSc); Turken 2017; Harris | sqrt(n) pooling + Daganzo + closed-form thresh  |
| 3-way cost/carbon/service      | narrowly OPEN                       | Wang 2025 (TR-E), Ni-Shu 2015 (svc-TIME)| fill-rate axis + stock-vs-speed flip            |
| Method: CA + carbon + pooling  | OPEN                                | Figliozzi (CA+carbon, constant EF)      | CA network design WITH pooling tradeoff         |
| Ext: load-factor emission fac. | mechanism OCCUPIED, open in CA      | Bektas-Laporte 2011; Bouchery 2016      | embed in CA network model                       |

## Per-threat differentiation

READ & CONFIRMED (full text in directory):
- Li et al. (2017, PLoS ONE) [file.pdf]: CONFIRMED not pre-empting. Single-objective MINLP solved in
  LINGO; carbon enters as a TAX inside cost (NOT a cap/constraint/Pareto); transport = DISCRETE mode
  choice (6 modes); only numerical sensitivity (DC locations shift with tax rate). No closed-form, no
  thresholds, no continuous speed, no speed-reversal, no centralization-direction theorem. Their stated
  novelty = two-echelon inventory interaction. Cite as "joint problem solved numerically."
- Drent, Moradi & Arts (2023, EJOR 311(3):925-941) [..723004538]: CONFIRMED not pre-empting. Computational
  MILP (Dantzig-Wolfe); TWO discrete modes; BACKLOG cost (not fill-rate); single-echelon supply-mode
  selection across an assortment; NO network density, NO continuous speed, NO closed-form reversal
  threshold. "Slow can be cheaper+greener" = dynamic portfolio-switching effect, a DIFFERENT mechanism
  than our safety-stock-carbon reversal. (NB: appendix uses NTM load-dependent ship emission factor ->
  citation support for our 2b load-factor extension.)
- Ni & Shu (2015, IJPR) [EBSCO..]: service-TIME + carbon (GSM); not fill rate; no transport speed lever.
- Bouchery et al. (2012, EJOR) [..712003505]: sustainable order quantity (EOQ), multi-objective Pareto,
  single-location; notes frequent small batches raise carbon if distances large, and refrigerated storage
  can offset full-truck-load advantage -> supports our two-bucket tension; no network/n.

- Cachon (2014, Management Science) [Cachon_2014...pdf]: READ. CLOSER THAN EXPECTED — this is our true
  anchor. He HAS: Voronoi tessellation, consumer round-trip d_c=phi_c*n^-1/2 (down in n), retailer TSP
  d_t=phi_t*n^1/2 (up in n), base-stock z*sigma with sqrt(1/n) pooling, retail-space cost ~ n^1/2,
  U-shaped C(n)=phi_c*c_c*n^-1/2+(phi_t*c_t+phi_s*c_s)*n^1/2, closed-form n*, carbon price folded into
  each cost coefficient, emission-min vs cost-min density, "carbon constraint can raise emissions."
  WHAT HE DOES NOT HAVE = our wedge: NO transport SPEED/lead-time/mode decision (he assumes ~zero lead
  time). So Theorem 2 (speed reversal) and the speed x density complementarity are genuinely outside
  Cachon. REFRAME the paper as "Cachon (2014) extended with an endogenous transport-speed/lead-time lever."
  His scaffold (Voronoi, phi_c n^-1/2, phi_t n^1/2, base-stock z*sigma_j, closed-form n*) = our method template.
- Tsao et al. (2012, EJOR) [Tsao_2012...pdf]: READ — CA integrated facility-inventory, NO carbon, NO speed.
  Confirms CA-network-design + carbon + speed is open; another method template.

- Schaefer & Konur (2015, TR-E) [preview read]: CONFIRMED not pre-empting — and a STRONG FOIL.
  Bi-objective (cost, emissions) continuous-review (Q,r), SINGLE location, stochastic demand; transport =
  DISCRETE carrier (LTL vs TL); Pareto-frontier approximation (numerical). Lead time is EXOGENOUS
  (sensitivity), not a decision; they conclude "reducing lead time is a green investment" — a ONE-SIDED
  result because they never price the transport-carbon COST of going faster. No network density, no
  continuous-speed U-shape, no reversal threshold. => We COMPLETE them: add eps(v)~v^2 (faster = dirtier
  transport) -> genuine U-shape in speed + the eta/h>eps/t reversal. Cite as the paper whose "faster=green"
  message we qualify.
- Turken, Carrillo & Verter (2017, IJPE) [preview read]: CONFIRMED not pre-empting (density). Facility
  location + capacity; centralize-vs-decentralize driven by PRODUCTION economies of scale vs transport
  emissions; carbon tax + command-and-control; MIP + special-case analytics; auto/GM dataset. NO inventory/
  safety stock/pooling, NO speed. Decentralize iff "low production pollution, low demand, or high transport
  cost." => density direction is now DOUBLY occupied (Cachon = inventory pooling; Turken = production scale).

STILL TO PULL: Wang/Shu/Su/Zhen (2025, TR-E) — confirm service=service-time, no stock-vs-speed substitution.

## REFINED CONTRIBUTION STATEMENT (after reading Cachon, Schaefer-Konur, Turken)
The paper UNIFIES and COMPLETES two literatures:
 - Cachon (2014, MnSc): carbon-optimal network DENSITY via inventory pooling, but ZERO lead time / no speed.
 - Schaefer-Konur (2015, TR-E): inventory + transport + emissions at a SINGLE location, but exogenous lead
   time and a one-sided "faster = green."
Our wedge = an ENDOGENOUS CONTINUOUS TRANSPORT SPEED (faster = dirtier transport eps(v)~v^2 AND less
safety stock) embedded in a Cachon-style network. This yields (a) Thm 2: non-monotone carbon-in-speed with
the eta/h>eps/t reversal (completes Schaefer-Konur), and (b) the speed x density complementarity (neither
paper has both levers). Density alone is NOT claimed novel (Cachon + Turken own it).

## Terminology to adopt

location-inventory model; transport-mode/speed selection; fill rate (type-2 service);
carbon intensity per dollar (eta/h, eps/t); insights from stylized models; total logistics carbon
curve (McKinnon, qualitative -> we formalize); emission-minimizing vs cost-minimizing density (Cachon).

## Priority full-text downloads (Ram)

1. Cachon (2014, Management Science)
2. Schaefer & Konur (2015, TR-E)
3. Li et al. (2017, PLoS ONE)
4. Wang, Shu, Su & Zhen (2025, TR-E)
Secondary: Hoen et al. (2014); Arikan & Jammernegg (2014); Turken et al. (2017).

## Open items / next steps

- [ ] Pull the 4 priority PDFs; confirm differentiators above.
- [ ] Solve joint (n*,v*) fixed point explicitly; verify SOC region vs corner solutions.
- [ ] Decide target journal (M&SOM / POM for structural; TR-E / Omega for richer) — TR-E is where the
      nearest competitors sit (Schaefer-Konur, Wang 2025), so framing must distinguish hard.
- [ ] Numerical study: calibrate rate vectors from Wyllo (eps mode factors, warehouse factor, grid),
      sweep tau, map product archetypes onto the 2x2; second testbed = Willems' 38 networks.
- [x] Load-factor (2b): DEFERRED to a future-research paragraph in model.tex (not core).

## Model revisions (post review-instruction assessment) — model.tex compiles, 9pp
Applied 7 fixes after assessing an external review-instruction set:
1. alpha relabeled CYCLE SERVICE LEVEL / in-stock probability (Type-1, per Cachon) -- was mislabeled "fill rate."
2. "U-shaped in n" now under an explicit SOC (last-mile curvature > concave pooling term); else CORNER (n=1 or N).
3. eps(v), t(v) assumed convex increasing in GENERAL; quadratic eps_0 v^2 is the canonical slow-steaming
   instantiation; v read as a general replenishment-speed/mode index (not only vessel speed).
4. Speed threshold eta/h > eps_0/t_1 stated CONDITIONAL on common n; added general-convex marginal form
   eps'(v)/t'(v); added joint-optimum note (complementarity reinforces when carbon also decentralizes).
5. Facility carbon phi (baseline, inventory-independent) vs holding carbon eta (per-unit-stored) defined as
   DISTINCT -> not double-counted.
6. Cap-price duality lambda=tau flagged as within-regime only (J non-convex in n via pooling; n integer).
7. Load-factor e(q)=a/q+b added as a deferred future-research paragraph (restores consolidation's transport
   benefit our core omits; expected to soften, not reverse, the density threshold).
