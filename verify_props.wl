(* Symbolic verification of the propositions in model.tex.
   Reserved names avoided: Nn=N, Dd=D, e0=eps0, eL=epsL, del=delta, sig=sigma,
   Kc=cost K, Cc=carbon C, ee=e coefficient, zb=z_beta. *)
$Assumptions = {F>0, phi>0, h>0, eta>0, t0>0, t1>0, e0>0, eL>0, tL>0,
  Nn>1, mu>0, sig>0, z>0, Dd>0, del>0, tau>0, 0<rho<1, n>1, n<Nn, v>0};

Psi[nn_] := Nn*rho + nn*(1-rho);
dd[nn_]  := del*nn^(-1/2);
ss[nn_,vv_] := z*sig*Sqrt[Dd/vv]*Sqrt[Nn*Psi[nn]];
eps[vv_] := e0*vv^2;
tt[vv_]  := t0 + t1*vv^2;
Kc[nn_,vv_] := F*nn + h*ss[nn,vv] + tt[vv]*Nn*mu*Dd + tL*Nn*mu*dd[nn];
Cc[nn_,vv_] := phi*nn + eta*ss[nn,vv] + eps[vv]*Nn*mu*Dd + eL*Nn*mu*dd[nn];
J[nn_,vv_,t_] := Kc[nn,vv] + t*Cc[nn,vv];

a = F + tau*phi;
b[vv_] := (h+tau*eta)*z*sig*Sqrt[Nn]*Sqrt[Dd/vv];
cco[vv_] := (tt[vv] + tau*eps[vv])*Nn*mu*Dd;
ee = (tL + tau*eL)*Nn*mu*del;
Jform[nn_,vv_] := a*nn + b[vv]*Sqrt[Psi[nn]] + cco[vv] + ee*nn^(-1/2);

pass[lbl_, cond_] := Print[lbl, ": ", If[TrueQ[cond], "PASS", Row[{"CHECK -> ", cond}]]];

Print["===== LEMMA 1: J = K + tau C equals the four-piece blended form ====="];
pass["J == a n + b(v) Sqrt[Psi] + c(v) + e n^-1/2",
     FullSimplify[J[n,v,tau] - Jform[n,v]] === 0];

Print["===== PROP 2: SPEED REVERSAL ====="];
pass["J strictly convex in v  (d2J/dv2 > 0)",
     Simplify[D[J[n,v,tau],{v,2}] > 0]];
vstar = ((h+tau*eta)*z*sig*Sqrt[Nn*Dd]*Sqrt[Psi[n]]/(4*(t1+tau*e0)*Nn*mu*Dd))^(2/5);
pass["v*(n) satisfies dJ/dv = 0",
     FullSimplify[(D[J[n,v,tau],v] /. v->vstar)] === 0];
ratio = (h+tau*eta)/(t1+tau*e0);
pass["d/dtau [(h+tau eta)/(t1+tau e0)] numerator == eta t1 - h e0",
     FullSimplify[Numerator[Together[D[ratio,tau]]] - (eta*t1 - h*e0)] === 0];
pass["v*|tau->inf > v*|tau=0  <=>  eta/h > e0/t1",
     Simplify[(Limit[ratio,tau->Infinity] > (ratio /. tau->0)) ==
              (eta/h > e0/t1)]];

Print["===== PROP 3: SERVICE & VOLATILITY ARE CARBON PARAMETERS ====="];
pass["dv*/dz > 0", Simplify[D[vstar,z] > 0]];
pass["dv*/dsig > 0", Simplify[D[vstar,sig] > 0]];
Jn = D[J[n,v,tau],n]; Jnn = D[J[n,v,tau],{n,2}];
holdLHS = b[v]*(1-rho)/(2*Sqrt[Psi[n]]);
pass["holding term in FOCn increasing in z", Simplify[D[holdLHS,z] > 0]];
pass["holding term in FOCn increasing in sig", Simplify[D[holdLHS,sig] > 0]];
(* IFT: dn*/dz = -Jnz/Jnn ; sign negative since Jnz>0, Jnn>0 under SOC *)
pass["Jnz > 0 (so dn*/dz < 0 under SOC: consolidation pull)",
     Simplify[D[Jn,z] > 0]];

Print["===== PROP 4: DENSITY & DISPERSION ====="];
pass["d2J/dn2 == (3/4) e n^-5/2 - (1/4) b(v)(1-rho)^2 Psi^-3/2",
     FullSimplify[D[J[n,v,tau],{n,2}]
        - (3/4*ee*n^(-5/2) - 1/4*b[v]*(1-rho)^2*Psi[n]^(-3/2))] === 0];
(* rho -> 1 : Psi = N, FOCn reduces to a = (e/2) n^-3/2 *)
Jn1 = D[J[nn,v,tau],nn] /. nn->n /. rho->1 // Simplify;
solN = Solve[(a - ee/2*n^(-3/2)) == 0 && n>0, n];
pass["rho=1: FOCn root n* = (e/2a)^(2/3)",
     FullSimplify[(n /. solN[[1]]) - (ee/(2 a))^(2/3)] === 0];
sratio = (tL+tau*eL)/(F+tau*phi);
pass["d/dtau [(tL+tau eL)/(F+tau phi)] numerator == eL F - tL phi",
     FullSimplify[Numerator[Together[D[sratio,tau]]] - (eL*F - tL*phi)] === 0];
pass["n*_carbon > n*_cost  <=>  eL/tL > phi/F",
     Simplify[(Limit[sratio,tau->Infinity] > (sratio /. tau->0)) ==
              (eL/tL > phi/F)]];

Print["===== PROP 5: COMPLEMENTARITY ====="];
pass["dv*/dn >= 0  (v* prop Psi^1/5, Psi' = 1-rho >= 0)",
     Simplify[D[vstar,n] >= 0]];
Jnv = D[J[n,v,tau],n,v];
pass["J_nv == -b(v)(1-rho)/(4 v Sqrt[Psi])  and < 0",
     FullSimplify[Jnv + b[v]*(1-rho)/(4*v*Sqrt[Psi[n]])] === 0 &&
     Simplify[Jnv < 0]];
pass["dn*/dv = -J_nv/J_nn > 0  under SOC (J_nn>0)",
     Simplify[Refine[-Jnv/Jnn > 0, Jnn > 0]]];

Print["===== PROP 7: DETERMINISTIC EQUIVALENT OF CHANCE CONSTRAINT ====="];
g = eps[v]*Dd + eL*dd[n];
sigQ = sig*Sqrt[Nn*(1+(Nn-1)*rho)];
muQ  = Nn*mu;
CVQ  = sigQ/muQ;
mfac = 1 + zb*CVQ;
pass["g*(muQ + zb sigQ) == m*(eps(v) N mu D + eL N mu d(n))",
     FullSimplify[g*(muQ + zb*sigQ) - mfac*(eps[v]*Nn*mu*Dd + eL*Nn*mu*dd[n])] === 0];
pass["CV_Q == (sig/mu) Sqrt[(1+(N-1)rho)/N]",
     FullSimplify[CVQ - (sig/mu)*Sqrt[(1+(Nn-1)*rho)/Nn]] === 0];
pass["Normal quantile: Quantile[N(muQ,sigQ),1-b] == muQ + sigQ*zb, zb=InvCDF",
     FullSimplify[Quantile[NormalDistribution[muQ,sigQ], 1-bb]
        - (muQ + sigQ*Quantile[NormalDistribution[],1-bb])] === 0];
pass["m >= 1  iff  beta <= 1/2  (zb = Phi^-1(1-beta) >= 0)",
     (Quantile[NormalDistribution[],1-1/2] === 0) &&
     Simplify[Quantile[NormalDistribution[],1-bb] >= 0, 0<bb<=1/2]];

Print["===== COR 8: UNCERTAINTY TIGHTENS THRESHOLDS ====="];
pass["m increasing in sig/mu  (dCVQ/dsig > 0)", Simplify[D[CVQ,sig] > 0]];
pass["m increasing in rho for N>1  (dCVQ/drho > 0)",
     Simplify[Refine[D[CVQ,rho] > 0, Nn>1]]];
pass["zb increasing in confidence 1-b", Simplify[D[Quantile[NormalDistribution[],cl],cl] > 0]];

Print["===== DONE ====="];
