# Attribution Map — which paper gets credit for each model element

Principle: we never overstate. Every borrowed building block is cited to its originator (and to the
paper we directly follow). Our claimed contribution is ONLY the speed layer + the speed x density
coupling on top of these. Everything else is reported and synthesized, not claimed.

| Model element (ours) | Exact form | Credit — originator | We directly follow / also |
|---|---|---|---|
| Continuous-approximation network modeling | region area a, n facilities | Daganzo (2005) Logistics Systems Analysis; Newell | Tsao et al. (2012) CA facility-inventory; Smilowitz-Daganzo |
| Last-mile distance falls as 1/sqrt(n) | d(n)=delta*n^(-1/2) | Beardwood, Halton & Hammersley (1959, TSP ~ sqrt(n)); Fisher & Hochbaum (1980, k-median) | Cachon (2014) uses phi_c*n^(-1/2) |
| Route / outbound distance grows as sqrt(n) | ~ phi_t*n^(1/2) | same TSP asymptotics | Cachon (2014) phi_t*n^(1/2) |
| Risk pooling / square-root law (total SS ~ sqrt(n)) | SS_total ~ z*sigma*sqrt(n) at rho=0 | Eppen (1979); Maister (1976, "square root law") | Cachon (2014) base-stock z*sigma_j, sigma_j~sqrt(1/n); location-inventory: Daskin-Coullard-Shen (2002), Shen-Coullard-Daskin (2003), Snyder-Daskin-Teo (2007) |
| Base-stock / safety stock over lead time | SS = z*sigma*sqrt(L), z=Phi^{-1}(alpha) | classical inventory theory (Zipkin); Simpson (1958) for placement | standard |
| Pooling with correlated demand | Psi(n)=N*rho+n*(1-rho) | Eppen (1979, independent); correlated extension standard | Silbermayr et al. (2017) copula dependence |
| Joint location-inventory model | open n DCs + inventory | Daskin, Coullard & Shen (2002); Shen, Coullard & Daskin (2003) | Ozsen-Coullard-Daskin (2008, capacitated) |
| Carbon as operational lever (vs tech investment); inventory carbon | two buckets: transport + holding | Hua, Cheng & Wang (2011); Benjaafar, Li & Daskin (2013) | Bouchery et al. (2012) |
| Carbon-optimal vs cost-optimal DENSITY; centralize/decentralize | U-shaped C(n), n* differs from cost | Cachon (2014, MnSc) [inventory pooling]; Turken et al. (2017, IJPE) [production scale] | McKinnon & Piecyk (2012, qualitative) |
| Transport mode/carrier + inventory + emissions | discrete carriers | Hoen et al. (2014); Konur (2014); Schaefer & Konur (2015) | Chen & Wang (2016) |
| Speed–emission relation eps(v) ~ v^2 (slow steaming, fuel ~ v^3) | convex increasing | Psaraftis & Kontovas (2013); Cariou (2011); Corbett-Wang-Winebrake | load+speed form: Bektas & Laporte (2011) / Barth |
| Lead time as a green lever | shorter L -> less SS | Arikan & Jammernegg (2014); Schaefer & Konur (2015) | -- |
| Carbon policy regimes (cap / tax / price / cap-and-trade); duality | min K s.t. C<=B == min K+tau*C | Hua-Cheng-Wang (2011); Benjaafar-Li-Daskin (2013) | cap-and-trade in placement: Wang, Shu, Su & Zhen (2025) |
| Guaranteed-service safety-stock PLACEMENT (hub/echelon extension only) | service-time decision | Simpson (1958); Graves & Willems (2000) | Ni & Shu (2015); Wang et al. (2025) |
| Load-factor endogenous emission factor (2b extension) | e(q)=a/q+b | Bektas & Laporte (2011); Bouchery et al. (2012) | NTM/GLEC accounting; Drent et al. (2023) appendix |
| Multi-objective cost-vs-carbon Pareto (method) | Pareto frontier | Bouchery et al. (2012) | Schaefer & Konur (2015) |
| "Insights from simple/stylized models" (methodological stance) | -- | Benjaafar, Li & Daskin (2013) | -- |
| Numerical validation testbed | real multi-echelon chains | Willems (2008) 38-chain dataset | (motivation sim: Wyllo) |

## OURS to claim (nothing borrowed here)
- Endogenous CONTINUOUS transport speed v with eps(v)~v^2 embedded in a Cachon-style network.
- Theorem 2: non-monotone total carbon in speed -> carbon-optimal speed > cost-optimal speed iff
  eta/h > eps_0/t_1 (completes Schaefer-Konur's one-sided "faster=green").
- Speed x density complementarity (joint (n,v)); neither Cachon (no speed) nor Schaefer-Konur (no density) has both.
- (Optional) service currency: stock vs speed for a fill-rate target (distinct from GSM service-time).

## Citations needing verification before use (classics)
- Beardwood, Halton & Hammersley (1959), "The shortest path through many points," Proc. Cambridge Phil. Soc. 55:299-327. [VERIFY pages]
- Fisher & Hochbaum (1980), "Probabilistic analysis of the planar k-median problem," Math. of OR 5(1):27-34. [VERIFY]
- Eppen (1979), "Effects of centralization on expected costs in a multi-location newsboy problem," Management Science 25(5):498-501. [VERIFY]
- Maister (1976), "Centralisation of inventories and the square root law," Int. J. Physical Distribution 6(3):124-134. [VERIFY]
- Simpson (1958), "In-process inventories," Operations Research 6(6):863-873. [VERIFY]
- Graves & Willems (2000), "Optimizing strategic safety stock placement in supply chains," M&SOM 2(1):68-83. [VERIFY]
