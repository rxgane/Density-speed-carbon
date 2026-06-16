$Assumptions = {F>0,phi>0,h>0,eta>0,t0>0,t1>0,e0>0,eL>0,tL>0,Nn>1,mu>0,sig>0,
  z>0,Dd>0,del>0,0<rho<1,n>1,v>0,G>0};
Psi[nn_]:=Nn*rho+nn*(1-rho);
SS[nn_,vv_]:=z*sig*Sqrt[Dd/vv]*Sqrt[Nn*Psi[nn]];
Kc[nn_,vv_]:=F*nn+h*SS[nn,vv]+(t0+t1*vv^2)*Nn*mu*Dd+tL*Nn*mu*del*nn^(-1/2);
Cc[nn_,vv_]:=phi*nn+eta*SS[nn,vv]+e0*vv^2*Nn*mu*Dd+eL*Nn*mu*del*nn^(-1/2);
pass[l_,c_]:=Print[l,": ",If[TrueQ[c],"PASS",Row[{"-> ",c}]]];

Print["===== POOLING-ADJUSTED DENSITY THRESHOLD (Prop 2 / new core) ====="];
Gexpr = D[SS[n,v],n];
pass["G = dSS/dn = z sig Sqrt[D/v] Sqrt[N] (1-rho)/(2 Sqrt[Psi])",
  FullSimplify[Gexpr - z*sig*Sqrt[Dd/v]*Sqrt[Nn]*(1-rho)/(2*Sqrt[Psi[n]])]===0];
dKdn = D[Kc[n,v],n];  dCdn = D[Cc[n,v],n];
(* cost FOC at n_K : dKdn=0  =>  tL Nn mu del /(2 n^{3/2}) = F + h G  *)
focRel = Solve[dKdn==0, del][[1]];      (* express del from cost FOC *)
(* substitute that del into dCdn and into G; then check sign condition *)
dCdn_atNK = FullSimplify[dCdn /. focRel];
(* claim: dCdn_atNK == phi + eta*G(n,v) - eL*(F + h*G(n,v))/tL  *)
Gn = Gexpr;
claim = phi + eta*Gn - eL*(F + h*Gn)/tL;
pass["dC/dn at n_K  ==  phi + eta G - eL (F + h G)/tL",
  FullSimplify[dCdn_atNK - claim]===0];
pass["=> n_C > n_K  iff  eL/tL > (phi + eta G)/(F + h G)   [via sign(dC/dn|nK)<0]",
  Simplify[(claim < 0) == (eL/tL > (phi + eta*Gn)/(F + h*Gn))]];

Print["===== MEDIANT PROPERTIES of the pooling-adjusted threshold R(G)=(phi+eta G)/(F+h G) ====="];
R[g_]:=(phi+eta*g)/(F+h*g);
pass["R(0) = phi/F (no pooling, rho=1, G=0)", Simplify[R[0]-phi/F]===0];
pass["R(G)->eta/h as G->Infinity", Simplify[Limit[R[g],g->Infinity]-eta/h]===0];
pass["dR/dG sign = sign(eta/h - phi/F) (threshold moves toward eta/h)",
  Simplify[Sign[D[R[g],g]] == Sign[eta*F - h*phi]]];

Print["===== GENERAL-CONVEX SPEED THRESHOLD (Prop 1) ====="];
A=z*sig*Sqrt[Nn*Dd]*Sqrt[Psi[n]]; T=Nn*mu*Dd;
Kv[vv_]:=h*A*vv^(-1/2)+T*tt[vv];  Cv[vv_]:=eta*A*vv^(-1/2)+T*ee[vv];
(* cost FOC: dKv/dv=0 => T tt'(vK) = (1/2) h A vK^{-3/2} *)
dKv=D[Kv[v],v]; dCv=D[Cv[v],v];
costFOC = (T*tt'[v] == (1/2)*h*A*v^(-3/2));     (* equivalent statement *)
pass["cost FOC: dKv/dv=0  <=>  T t'(v) = (1/2) h A v^{-3/2}",
  FullSimplify[(dKv==0) == costFOC]];
(* substitute (1/2)A v^{-3/2} = T t'(v)/h into dCv *)
dCv_sub = dCv /. ( (1/2)*A*v^(-3/2) -> T*tt'[v]/h );
pass["dC_v/dv at v_K == T ( eps'(vK) - (eta/h) t'(vK) )",
  FullSimplify[dCv_sub - T*(ee'[v] - (eta/h)*tt'[v])]===0];
pass["=> v_C > v_K  iff  eta/h > eps'(vK)/t'(vK)   [assuming t'>0]",
  Simplify[(T*(ee'[v]-(eta/h)*tt'[v]) < 0) == (eta/h > ee'[v]/tt'[v]), tt'[v]>0]];
(* quadratic specialization *)
pass["quadratic: eps'(v)/t'(v) = e0/t1 with eps=e0 v^2, t=t0+t1 v^2",
  Simplify[(2*e0*v)/(2*t1*v) - e0/t1]===0];
Print["===== DONE ====="];
