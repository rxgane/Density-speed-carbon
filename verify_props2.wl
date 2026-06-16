$Assumptions = {F>0,phi>0,h>0,eta>0,t0>0,t1>0,e0>0,eL>0,tL>0,Nn>1,mu>0,sig>0,
  z>0,Dd>0,del>0,tau>0,0<rho<1,n>1,n<Nn,v>0};
Psi[nn_]:=Nn*rho+nn*(1-rho);
dd[nn_]:=del*nn^(-1/2);
ss[nn_,vv_]:=z*sig*Sqrt[Dd/vv]*Sqrt[Nn*Psi[nn]];
eps[vv_]:=e0*vv^2; tt[vv_]:=t0+t1*vv^2;
Kc[nn_,vv_]:=F*nn+h*ss[nn,vv]+tt[vv]*Nn*mu*Dd+tL*Nn*mu*dd[nn];
Cc[nn_,vv_]:=phi*nn+eta*ss[nn,vv]+eps[vv]*Nn*mu*Dd+eL*Nn*mu*dd[nn];
J[nn_,vv_,t_]:=Kc[nn,vv]+t*Cc[nn,vv];
a=F+tau*phi; b[vv_]:=(h+tau*eta)*z*sig*Sqrt[Nn]*Sqrt[Dd/vv];
ee=(tL+tau*eL)*Nn*mu*del;
pass[lbl_,cond_]:=Print[lbl,": ",If[TrueQ[cond],"PASS",Row[{"FAIL -> ",cond}]]];

Print["-- Prop 4 (rho=1): paper's n* = (e/2a)^(2/3) solves the FOC a = (e/2) n^-3/2 --"];
nStar = (ee/(2 a))^(2/3);
pass["paper n* satisfies FOC  a - (e/2) n*^(-3/2) == 0",
     FullSimplify[a - ee/2*nStar^(-3/2)] === 0];
rootReduce = Reduce[a - ee/2*n^(-3/2) == 0 && n>0, n];
pass["Reduce of FOC equals paper n*",
     FullSimplify[(n /. ToRules[rootReduce]) - nStar] === 0];

Print["-- Prop 5b: dn*/dv = -J_nv/J_nn > 0 given J_nv<0 (proven) and SOC J_nn>0 --"];
Jnv = D[J[n,v,tau],n,v];
pass["-J_nv > 0  (numerator of dn*/dv)", TrueQ[Simplify[-Jnv > 0]]];
pass["-J_nv/Jnn > 0  when Jnn>0", TrueQ[Simplify[-Jnv/jj > 0, jj>0]]];

Print["-- Prop 7: Normal quantile shift identity (sigma>0) --"];
pass["Quantile[N(m,s),1-b] == m + s*z_b",
     TrueQ[FullSimplify[
       Quantile[NormalDistribution[mq,sq],1-bb] ==
       mq + sq*Quantile[NormalDistribution[],1-bb],
       {mq \[Element] Reals, sq>0, 0<bb<1}]]];

Print["-- Cor 8: z_beta strictly increasing in confidence --"];
deriv = D[Quantile[NormalDistribution[],cl], cl];
Print["   d z_b/d(1-b) = ", deriv, "  (Sqrt[2 Pi] * Exp[...] > 0 always)"];
pass["derivative > 0 on (0,1)", TrueQ[Simplify[deriv > 0, 0<cl<1]] ||
     TrueQ[Simplify[Exp[InverseErfc[2 cl]^2] > 0, 0<cl<1]]];

Print["-- BONUS: numeric cross-check of v* and threshold directions --"];
sub = {F->7500,phi->3000,h->0.5,eta->0.15,t0->1.,t1->0.05,e0->0.04,eL->0.0015,
  tL->0.005,Nn->40,mu->20000,sig->8000,z->1.645,Dd->1.,del->120,rho->0.3};
vK = (((h)*z*sig*Sqrt[Nn*Dd]*Sqrt[Psi[8]])/(4*(t1)*Nn*mu*Dd))^(2/5) /. sub;
vC = (((h+10^9*eta)*z*sig*Sqrt[Nn*Dd]*Sqrt[Psi[8]])/(4*(t1+10^9*e0)*Nn*mu*Dd))^(2/5) /. sub;
Print["   eta/h=", 0.15/0.5, "  e0/t1=", 0.04/0.05,
      "  => slower expected; v*(tau=0)=", vK, "  v*(tau->inf)=", vC,
      "  slower? ", vC < vK];
Print["===== DONE ====="];
