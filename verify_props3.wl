$Assumptions = {F>0,phi>0,h>0,eta>0,t0>0,t1>0,e0>0,eL>0,tL>0,Nn>1,mu>0,sig>0,
  z>0,Dd>0,del>0,tau>0,0<rho<1,n>1,n<Nn,v>0};
Psi[nn_]:=Nn*rho+nn*(1-rho); dd[nn_]:=del*nn^(-1/2);
ss[nn_,vv_]:=z*sig*Sqrt[Dd/vv]*Sqrt[Nn*Psi[nn]];
eps[vv_]:=e0*vv^2; tt[vv_]:=t0+t1*vv^2;
Kc[nn_,vv_]:=F*nn+h*ss[nn,vv]+tt[vv]*Nn*mu*Dd+tL*Nn*mu*dd[nn];
Cc[nn_,vv_]:=phi*nn+eta*ss[nn,vv]+eps[vv]*Nn*mu*Dd+eL*Nn*mu*dd[nn];
J[nn_,vv_,t_]:=Kc[nn,vv]+t*Cc[nn,vv];
a=F+tau*phi; ee=(tL+tau*eL)*Nn*mu*del;
pass[lbl_,cond_]:=Print[lbl,": ",If[TrueQ[cond],"PASS",Row[{"FAIL -> ",cond}]]];

pass["Prop4 rho=1: paper n*=(e/2a)^(2/3) solves FOC",
  FullSimplify[a-ee/2*((ee/(2 a))^(2/3))^(-3/2)]===0];
Jnv=D[J[n,v,tau],n,v];
pass["Prop5b: -J_nv>0 (numerator of dn*/dv)", TrueQ[Simplify[-Jnv>0]]];
pass["Prop5b: -J_nv/Jnn>0 when Jnn>0", TrueQ[Simplify[-Jnv/jj>0, jj>0]]];
pass["Prop7: Quantile[N(m,s),1-b]==m+s z_b (s>0)",
  TrueQ[FullSimplify[Quantile[NormalDistribution[mq,sq],1-bb]==
    mq+sq*Quantile[NormalDistribution[],1-bb],{mq\[Element]Reals,sq>0,0<bb<1}]]];
deriv=D[Quantile[NormalDistribution[],cl],cl];
pass["Cor8: dz_b/dconf>0", TrueQ[Simplify[Exp[InverseErfc[2 cl]^2]>0,0<cl<1]]];
Print["   (dz_b/dconf = ", Simplify[deriv], ")"];
Print["===== DONE ====="];
