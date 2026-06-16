$Assumptions={F>0,phi>0,h>0,eta>0,t0>0,t1>0,e0>0,eL>0,tL>0,Nn>1,mu>0,sig>0,
  z>0,Dd>0,del>0,tau>0,0<rho<1,n>1,n<Nn,v>0};
Psi[nn_]:=Nn*rho+nn*(1-rho); dd[nn_]:=del*nn^(-1/2);
ss[nn_,vv_]:=z*sig*Sqrt[Dd/vv]*Sqrt[Nn*Psi[nn]];
eps[vv_]:=e0*vv^2; tt[vv_]:=t0+t1*vv^2;
Kc[nn_,vv_]:=F*nn+h*ss[nn,vv]+tt[vv]*Nn*mu*Dd+tL*Nn*mu*dd[nn];
Cc[nn_,vv_]:=phi*nn+eta*ss[nn,vv]+eps[vv]*Nn*mu*Dd+eL*Nn*mu*dd[nn];
pass[l_,c_]:=Print[l,": ",If[TrueQ[c],"PASS",Row[{"-> ",c}]]];

Print["== PIPELINE: add (h+tau eta) N mu D / v (in-transit at holding rates) =="];
Jp = Kc[n,v] + tau*Cc[n,v] + (h+tau*eta)*Nn*mu*Dd/v;
foc = D[Jp,v];
A = z*sig*Sqrt[Nn*Dd]*Sqrt[Psi[n]]; P = Nn*mu*Dd;
(* claim: foc==0  <=>  (h+tau eta)/(t1+tau e0) = RHS(v), RHS strictly increasing *)
RHS = 2*Nn*mu*Dd*v/((1/2)*A*v^(-3/2) + P*v^(-2));
pass["foc factors as (t1+tau e0)2NmuD v - (h+tau eta)[(1/2)A v^-3/2 + P v^-2]",
  FullSimplify[foc - ((t1+tau*e0)*2*Nn*mu*Dd*v - (h+tau*eta)*((1/2)*A*v^(-3/2)+P*v^(-2)))]===0];
pass["RHS(v) strictly increasing in v  => v* monotone in (h+tau eta)/(t1+tau e0)",
  Simplify[D[RHS,v] > 0]];
pass["ratio (h+tau eta)/(t1+tau e0) increasing in tau iff eta/h > e0/t1 (threshold UNCHANGED)",
  Simplify[Sign[D[(h+tau*eta)/(t1+tau*e0),tau]] == Sign[eta*t1 - h*e0]]];

Print["== PIPELINE (cost only, no warehouse carbon): added to K, not C =="];
Jp2 = Kc[n,v] + h*Nn*mu*Dd/v + tau*Cc[n,v];
focCost = (D[Jp2,v]/.tau->0);   (* cost-min v FOC with pipeline *)
focCostNoPipe = (D[Kc[n,v],v]); (* cost-min v FOC without pipeline *)
pass["pipeline-cost adds a strictly negative term to dK/dv (pushes cost-min v UP)",
  Simplify[(focCost - focCostNoPipe) < 0]];

Print["== PROP 3: b(v*(n)) proportional to (z sigma)^(4/5) =="];
vstar=((h+tau*eta)*z*sig*Sqrt[Nn*Dd]*Sqrt[Psi[n]]/(4*(t1+tau*e0)*Nn*mu*Dd))^(2/5);
bvstar=(h+tau*eta)*z*sig*Sqrt[Nn]*Sqrt[Dd/vstar];
pass["b(v*) / (z sig)^(4/5) is independent of z",  FullSimplify[D[bvstar/(z*sig)^(4/5),z]]===0];
pass["b(v*) / (z sig)^(4/5) is independent of sig",FullSimplify[D[bvstar/(z*sig)^(4/5),sig]]===0];
Print["DONE"];
