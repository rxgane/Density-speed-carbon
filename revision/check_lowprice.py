import numpy as np
from scipy.optimize import minimize
rng=np.random.default_rng(7)
def make(P):
    N,rho,F,h,eta,hp,t1,e1,tL,eL,phi,kap,LN,D=P
    Psi=lambda n: N*rho+n*(1-rho)
    def K(x,y): n,v=np.exp(x),np.exp(y); S=kap*np.sqrt(N*Psi(n)/v); return F*n+h*S+hp*D/v+t1*v*v+tL*LN*n**-0.5
    def C(x,y): n,v=np.exp(x),np.exp(y); S=kap*np.sqrt(N*Psi(n)/v); return phi*n+eta*S+e1*v*v+eL*LN*n**-0.5
    def opt(tau,z0):
        r=minimize(lambda z: K(*z)+tau*C(*z),z0,method="Nelder-Mead",options=dict(xatol=1e-12,fatol=1e-14,maxiter=20000)); return r.x
    def stats(z):
        x,y=z; n,v=np.exp(x),np.exp(y); S=kap*np.sqrt(N*Psi(n)/v); gx=0.5*(1-rho)*n/Psi(n)
        M=n*F+h*S*gx; w=h*S*gx/M; a=h*S/(h*S+2*hp*D/v)
        return w,a
    return K,C,opt,stats
agree=0; tot=0; flips=0
for it in range(400):
    N=40; rho=rng.uniform(0,0.9); F=rng.uniform(0.2,5); h=rng.uniform(0.5,3); rh=rng.uniform(0.001,1.0); eta=rh*h
    hp=rng.choice([0.0,rng.uniform(0.1,3)]); t1=rng.uniform(0.2,3); rt=rng.uniform(0.005,3); e1=rt*t1
    tL=1.0; rL=rng.uniform(0.01,1.5); eL=rL*tL; rF=rng.uniform(0.01,1.0); phi=rF*F; kap=rng.uniform(0.05,0.5); LN=rng.uniform(5,60); D=rng.uniform(1,10)
    K,C,opt,stats=make((N,rho,F,h,eta,hp,t1,e1,tL,eL,phi,kap,LN,D))
    z0=opt(0.0,np.array([1.0,0.0])); n0=np.exp(z0[0])
    if not (1.05<n0<N-1): continue
    w,a=stats(z0)
    thr=(1-w)*rF+w*(rh-(a*rh-rt)/(6-a))
    if abs(rL-thr)<0.02*max(thr,1e-3): continue      # skip near-ties
    dt=1e-4; z1=opt(dt,z0); dn=np.exp(z1[0])-n0
    if abs(dn)<1e-9: continue
    tot+=1; agree+= (np.sign(dn)==np.sign(rL-thr))
    endthr=(1-w)*rF+w*rh**0.8*rt**0.2*a**-0.2
    if np.sign(rL-thr)!=np.sign(rL-endthr): flips+=1
print(f"local condition matches the numerical sign in {agree}/{tot} interior instances; {flips} instances where local and endpoint thresholds disagree")
