"""Fermeture rapide 2: (a) caractere T3 des modes d'ancrage; (b) anticroisement kz=2.5-2.8."""
import numpy as np, scipy.sparse.linalg as spl
from track import buildU, core_mask
ports=[(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
pidx={P:i for i,P in enumerate(ports)}
Nx=15; c=Nx//2; core=[(c,c),(c+1,c),(c,c+1)]
def gi(ix,iy,P): return 12*((ix%Nx)*Nx+(iy%Nx))+pidx[P]
# operateur T3: permute les colonnes du coeur c0->c1->c2->c0 et applique R(2pi/3)
# aux paires de polarisation des ports z du coeur (structure meme de la tresse)
dim=12*Nx*Nx
T=np.zeros((dim,dim))
th=2*np.pi/3; R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
for ix in range(Nx):
    for iy in range(Nx):
        if (ix,iy) in core:
            j=core.index((ix,iy)); tx,ty=core[(j+1)%3]
            for P in ports:
                d,s,p=P
                if d==2:
                    for a in range(2):
                        T[gi(tx,ty,(2,s,a)), gi(ix,iy,(2,s,p))]=R[a,p]
                else:
                    T[gi(tx,ty,P), gi(ix,iy,P)]=1.0
        else:
            for P in ports: T[gi(ix,iy,P), gi(ix,iy,P)]=1.0
kz=2.0
U=buildU(Nx,kz,core); mask=core_mask(Nx,core)
print("(a) caractere T3 des modes lies a kz=2.0 (racines cubiques attendues: 1, e^{+-2pi i/3}):")
for sig in (0.3325, 0.3923):
    w,V=spl.eigs(U,k=4,sigma=np.exp(1j*sig),which='LM')
    om=np.angle(w)
    P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True); wt=P2[mask].sum(0)
    i=np.argmax(wt)
    psi=V[:,i]/np.linalg.norm(V[:,i])
    ch=psi.conj()@(T@psi)
    print(f"  om={om[i]:.5f} poids={wt[i]:.3f}  <T3>={np.abs(ch):.3f} * e^(i {np.angle(ch):+.3f})   [2pi/3={2*np.pi/3:.3f}]")
print("\n(b) anticroisement: deux etats suivis, kz=2.40..2.90")
for kzv in np.arange(2.40,2.91,0.10):
    w,V=spl.eigs(buildU(Nx,kzv,core),k=8,sigma=np.exp(1j*0.29),which='LM')
    om=np.angle(w); P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True); wt=P2[mask].sum(0)
    o=np.argsort(-wt)[:2]
    print(f"  kz={kzv:.2f}: om1={om[o[0]]:.5f} (p={wt[o[0]]:.2f})  om2={om[o[1]]:.5f} (p={wt[o[1]]:.2f})  ecart={abs(om[o[0]]-om[o[1]]):.5f}")
