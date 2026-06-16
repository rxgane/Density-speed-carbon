# Literature Review — scaffold (grows as full texts arrive)

Gap statement: pairwise combinations are published; the THREE-lever (network density + transport
speed + safety stock) joint problem under cost AND carbon AND availability, solved in CLOSED FORM with
structural thresholds, is open. Each pairwise strand below is "what's done"; our wedge is the join + the
analytics.

## Strand 1 — Location-inventory & risk pooling (our base model)
- Shen, Coullard & Daskin (2003); Daskin, Coullard & Shen (2002) — founding location-inventory, sqrt(n) pooling. [PDFs in dir]
- Snyder, Daskin & Teo (2007) — stochastic location with risk pooling.
- Ozsen, Coullard & Daskin (2008) — capacitated warehouse + risk pooling.
- Park, Lee & Sung (2010) — three-level network design with risk pooling + lead times.
- Tsao, Mangotra, Lu & Dong (2012) — CONTINUOUS APPROXIMATION facility-inventory (our method, NO carbon).
- Maister (1976) — square-root law (centralization).
GAP: none add transport speed/mode or carbon.

## Strand 2 — Carbon-aware inventory / EOQ (2-axis cost vs carbon)
- Hua, Cheng & Wang (2011) — carbon EOQ under cap-and-trade.
- Benjaafar, Li & Daskin (2013) — "insights from simple models" (our methodological license). [PDF in dir]
- Bouchery, Ghaffari, Jemai & Dallery (2012) — sustainable order quantity, Pareto; refrigerated storage
  offsets full-truck-load advantage (supports two-bucket tension). [PDF in dir]
- Arikan & Jammernegg (2014) — dual-sourcing newsvendor + carbon footprint (stock-vs-speed cousin).
GAP: single location; no network density; no continuous speed.

## Strand 3 — Carbon + transport mode / speed (-> Theorem 2)
- Hoen, Tan, Fransoo & van Houtum (2014, FSM + Transp Sci) — mode choice + inventory + carbon (discrete). [PDF in dir]
- Chen & Wang (2016, TR-E) — mode-shifting thresholds under carbon policies.
- Konur (2014, TR-E); Schaefer & Konur (2015, TR-E) — continuous-review inventory + transport + emissions [PULL].
- Drent, Moradi & Arts (2023, EJOR) — dynamic dual-mode selection, MILP; portfolio-switching reversal. [PDF in dir]
- Slow steaming: Psaraftis & Kontovas (2013); Cariou (2011) — the "slow=green" orthodoxy we overturn.
GAP: discrete modes; no continuous-speed U-shape; no eta/h>eps/t reversal threshold; no network density.

## Strand 4 — Carbon + network design / centralization (-> Theorem 3)
- Cachon (2014, Management Science) — density vs emissions (consumer driving; no inventory pooling). [PULL — heavyweight]
- Turken, Carrillo & Verter (2017, IJPE) — centralize vs decentralize under carbon (production scale; MIP).
- Harris et al. (2011, IJPE) — numerical "carbon decentralizes".
- Kohn & Brodin (2008) — "consolidate=green" foil (needs modal shift).
- McKinnon & Piecyk (2012) — qualitative total-logistics-CO2 curve (we formalize).
- Li, Li, Zhang & Zhou (2017, PLoS ONE) — JOINT 3-lever, NUMERICAL MINLP/LINGO, carbon tax. [PDF in dir]
- Alhaj, Svetinovic & Diabat (2016) — carbon-sensitive two-echelon inventory.
GAP: direction known; our wedge = sqrt(n) pooling + Daganzo last-mile + closed-form threshold + "pooling is the real green lever".

## Strand 5 — Carbon + safety-stock placement (GSM)
- Ni & Shu (2015, IJPR) — service-time vs carbon, GSM. [PDF in dir]
- Wang, Shu, Su & Zhen (2025, TR-E) — GSM cap-and-trade. [PULL]
GAP: lever is service TIME, not fill rate; transport speed/mode fixed.

## Strand 6 — Continuous-approximation methods (our toolkit)
- Daganzo (2005) — Logistics Systems Analysis (distance ~ 1/sqrt(n)).
- Tsao et al. (2012) — CA facility-inventory.
- Figliozzi (2017, 2020) — CA + carbon, last-mile (CONSTANT emission factor).
- Saberi & Verbas (2012) — CA emissions-minimizing routing.
GAP: CA + carbon exists only for routing w/ constant EF; CA network design is cost-only.

## Strand 7 — Load-factor / consolidation emissions (-> 2b extension)
- Bektas & Laporte (2011, TR-B) — pollution-routing; load-and-speed-dependent emissions (canonical form).
- Bouchery et al. (2012/2016) — fewer/fuller shipments cut emissions.
- NTM / GLEC accounting; Drent et al. (2023) appendix uses NTM load-dependent ship factor.
- Adopt per-unit form e(q) = a/q + b.

---

## Analysis techniques to borrow (log as PDFs arrive)
- Piecewise-linear approximation of nonlinear carbon/safety-stock terms (Ni & Shu 2015).
- Dantzig-Wolfe / column generation for multi-item emission-capped control (Drent et al. 2023).
- Lagrangian relaxation for location-inventory (Daskin-Coullard-Shen 2002).
- Continuous approximation (Daganzo; Tsao et al. 2012) — our core method.
- Copulas for correlated demand (Silbermayr et al. 2017) — for the rho / correlation analysis.
- Multi-objective Pareto + interactive selection (Bouchery et al. 2012) — for the cost/carbon/service frontier.
- [add per incoming paper: technique -> where it could help our analysis]
