$Assumptions = {F>0,phi>0,h>0,eta>0,tL>0,eL>0,G>0,L>0,M>0,T>0,tp>0,ep>0,g>0};
pass[l_,c_]:=Print[l,": ",If[TrueQ[c],"PASS",Row[{"-> ",c}]]];

Print["== DENSITY: marginal balance. dK/dn = F + hG - L ; dC/dn = phi + eta G - (eL/tL)L =="];
Print["   (G=dSS/dn>0 pooling penalty; L=last-mile marginal magnitude; cost FOC: L=F+hG)"];
dCdn = phi + eta*G - (eL/tL)*L;
pass["dC/dn at n_K equals phi + eta G - eL(F+hG)/tL",
  FullSimplify[(dCdn /. L->F+h*G) - (phi+eta*G - eL*(F+h*G)/tL)]===0];
pass["n_C>n_K  iff  eL/tL > (phi+eta G)/(F+hG)",
  Simplify[((phi+eta*G - eL*(F+h*G)/tL) < 0) == (eL/tL > (phi+eta*G)/(F+h*G))]];
R[x_]:=(phi+eta*x)/(F+h*x);
pass["R(0)=phi/F", Simplify[R[0]==phi/F]];
pass["R(G)->eta/h as G->Inf", Simplify[Limit[R[x],x->Infinity]==eta/h]];
pass["dR/dG sign = sign(eta/h - phi/F)  (boundary tilts toward eta/h)",
  Simplify[Sign[D[R[x],x]]==Sign[eta/h-phi/F], F+h*x>0]];
pass["R(G) lies between phi/F and eta/h",
  Simplify[(R[g]-phi/F)*(R[g]-eta/h) <= 0]];

Print["== SPEED: marginal balance. dKv = T t'(v) - hM ; dCv = T eps'(v) - eta M =="];
Print["   (M=(1/2)A v^{-3/2}>0 holding-speed marginal; cost FOC: hM = T t'(vK))"];
dCv = T*ep - eta*M;
pass["dC_v/dv at v_K equals T(eps'(vK) - (eta/h) t'(vK))",
  FullSimplify[(dCv /. M->T*tp/h) - T*(ep-(eta/h)*tp)]===0];
pass["v_C>v_K  iff  eta/h > eps'(vK)/t'(vK)",
  Simplify[(T*(ep-(eta/h)*tp) < 0) == (eta/h > ep/tp), tp>0]];
Print["== DONE =="];
