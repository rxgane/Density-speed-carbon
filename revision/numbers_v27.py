import numpy as np, json, application_v21 as A, application_v22 as B, targets_v22 as T
from scipy.optimize import minimize_scalar
Y=365.0; p=B.reference(); nK,vK=T.solve(p); K0,C0=B.KC(nK,vK,p)
S=A.SS(nK,vK,p); P=p["N"]*p["mu"]*p["D"]/vK
G=p["z"]*p["sigma"]*np.sqrt(p["D"]/vK)*np.sqrt(p["N"])*(1-p["rho"])/(2*np.sqrt(A.Psi(nK,p)))
w=p["h"]*G/(p["F"]+p["h"]*G); aK=B.share(nK,vK,p); rh,rt,rF,rL=p["eta"]/p["h"],p["eps0"]/p["t1"],p["phi"]/p["F"],p["rL"]
vC=minimize_scalar(lambda v: B.KC(nK,v,p)[1],bounds=(p["vmin"],p["vmax"]),method="bounded").x
endthr=(1-w)*rF+w*rh*np.sqrt(vK/vC); locthr=(1-w)*rF+w*(rh-(aK*rh-rt)/(6-aK))
print(f"costmin n {nK:.3f} days {1000/vK:.3f} K {K0*Y/1e6:.2f} C {C0*Y/1e3:,.0f} SS {S:,.0f} transit {P:,.0f} aK {aK:.3f} omega {w:.3f}")
print(f"ratios r_h {rh:.4f} r_t {rt:.3f} r_F {rF:.4f} r_L {rL}; carbon-min speed at nK {1000/vC:.2f} d; endpoint threshold {endthr:.4f}; low-price threshold {locthr:.4f}")
for days in (4.0,):
    v=1000/days; K,C=B.KC(nK,v,p); print(f"at {days} days with nK: cost {100*(K/K0-1):+.1f}% carbon {100*(C/C0-1):+.1f}%")
print("line haul cost/carbon per pallet-mile 6d->2d: $%.4f->$%.4f, %.4f->%.4f kg"%(p['t0']+p['t1']*p['vmin']**2,p['t0']+p['t1']*p['vmax']**2,p['e0']+p['eps0']*p['vmin']**2,p['e0']+p['eps0']*p['vmax']**2))
print("intermodal 4d per pallet-mile: $%.4f, %.4f kg"%(p['t0']+p['t1']*250**2,p['e0']+p['eps0']*250**2))
# whole-center price jumps relative to the whole-center cost minimum
ni,vi=B.design(p,0.0,integer=True); Ki,Ci=B.KC(ni,vi,p); prev=None
for Pp in np.arange(100,1500,5):
    n,v=B.design(p,Pp/1000,integer=True); K,C=B.KC(n,v,p)
    if prev and n!=prev[0]: print(f"whole-center jump ${prev[2]:.0f}->{Pp:.0f}/t: {prev[0]} DCs cut {prev[1]:.1f}% -> {n} DCs cut {100*(1-C/Ci):.1f}%")
    prev=(n,100*(1-C/Ci),Pp)
# rounding penalty over the variants and prices
def vopt(n,q,tau):
    f=(lambda v: B.KC(n,v,q)[0]+tau*B.KC(n,v,q)[1]) if np.isfinite(tau) else (lambda v: B.KC(n,v,q)[1])
    return minimize_scalar(f,bounds=(q["vmin"],q["vmax"]),method="bounded",options={"xatol":1e-9}).x
def obj(n,q,tau):
    v=vopt(n,q,tau); K,C=B.KC(n,v,q); return K+tau*C if np.isfinite(tau) else C
def nopt(q,tau):
    g=lambda n: obj(n,q,tau); ns=np.exp(np.linspace(0,np.log(q["N"]),400)); i=int(np.argmin([g(n) for n in ns]))
    return minimize_scalar(g,bounds=(ns[max(i-1,0)],ns[min(i+1,399)]),method="bounded",options={"xatol":1e-9}).x
worst=(0,None); worst_fin=(0,None)
for nm,kw in {"reference":{},"capital 10%":dict(capital=0.10),"capital 25%":dict(capital=0.25),"$10K":dict(value=10_000.0),"$50K":dict(value=50_000.0),"refrigerated":dict(storage="refrigerated"),"electric delivery":dict(rL=0.16)}.items():
    q=B.reference(**kw)
    for tau in (0.0,0.1,0.5,1.0,2.0,np.inf):
        n=nopt(q,tau); best=obj(n,q,tau); bi=min(obj(k,q,tau) for k in (int(np.floor(n)),int(np.ceil(n))) if 1<=k<=q["N"]); pen=bi/best-1
        if pen>worst[0]: worst=(pen,(nm,tau,round(n,2)))
        if np.isfinite(tau) and pen>worst_fin[0]: worst_fin=(pen,(nm,tau,round(n,2)))
print("rounding: largest penalty %.3f%% %s; at finite prices %.3f%% %s"%(100*worst[0],worst[1],100*worst_fin[0],worst_fin[1]))
