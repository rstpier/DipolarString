"""Variantes discretes du defaut D3 a kz=2.0: la frequence liee depend-elle d'un choix continu?
 v1: tresse + R(2pi/3) [reference, derive de la geometrie du triplet]
 v2: tresse pure (R=I)  [teste le role de la rotation de polarisation]
 v3: tresse + R(2pi/3) + phase scalaire e^{i 2pi/3} [autre branche discrete permise par la fermeture]
 v4: tresse + R(2pi/3) + phase scalaire e^{i 0.7} [PHASE ARBITRAIRE - test du bouton continu:
     si le spectre bouge continument avec elle ET qu'elle n'est pas fixee, condition de rejet]
"""
import numpy as np
from scipy.linalg import eig
ports = [(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
pidx = {P:i for i,P in enumerate(ports)}
S1 = np.load("S_johns.npy")
def run(kz, Nx, theta, extra):
    Nn=Nx*Nx; dim=12*Nn
    def gi(ix,iy,P): return 12*((ix%Nx)*Nx+(iy%Nx))+pidx[P]
    Sb=np.zeros((dim,dim))
    for n in range(Nn): Sb[12*n:12*n+12,12*n:12*n+12]=S1
    C=np.zeros((dim,dim),complex)
    R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])*np.exp(1j*extra)
    core=[(Nx//2,Nx//2),(Nx//2+1,Nx//2),(Nx//2,Nx//2+1)]; cs=set(core)
    for ix in range(Nx):
        for iy in range(Nx):
            for p in (1,2):
                C[gi(ix+1,iy,(0,-1,p)),gi(ix,iy,(0,+1,p))]=1
                C[gi(ix-1,iy,(0,+1,p)),gi(ix,iy,(0,-1,p))]=1
            for p in (0,2):
                C[gi(ix,iy+1,(1,-1,p)),gi(ix,iy,(1,+1,p))]=1
                C[gi(ix,iy-1,(1,+1,p)),gi(ix,iy,(1,-1,p))]=1
            if (ix,iy) in cs:
                j=core.index((ix,iy)); up=core[(j+1)%3]; dn=core[(j-1)%3]
                for a in range(2):
                    for b in range(2):
                        C[gi(up[0],up[1],(2,-1,a)),gi(ix,iy,(2,+1,b))]=R[a,b]*np.exp(1j*kz)
                        C[gi(dn[0],dn[1],(2,+1,a)),gi(ix,iy,(2,-1,b))]=np.conj(R.T)[a,b]*np.exp(-1j*kz)
            else:
                for p in (0,1):
                    C[gi(ix,iy,(2,-1,p)),gi(ix,iy,(2,+1,p))]=np.exp(1j*kz)
                    C[gi(ix,iy,(2,+1,p)),gi(ix,iy,(2,-1,p))]=np.exp(-1j*kz)
    U=C@Sb
    w,V=eig(U); om=np.angle(w)
    mask=np.zeros(dim,bool)
    for (cx,cy) in core: mask[12*((cx%Nx)*Nx+(cy%Nx)):12*((cx%Nx)*Nx+(cy%Nx))+12]=True
    P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True); wt=P2[mask].sum(0)
    sel=(om>0.02)&(om<kz/2-1e-6)
    o=np.argsort(-(wt*sel))[:2]
    return [(om[i],wt[i]) for i in o if sel[i]], np.abs(np.abs(w)-1).max()
kz=2.0; Nx=11
for name,th,ex in (("v1 R(2pi/3)",2*np.pi/3,0.0),("v2 R=I",0.0,0.0),
                   ("v3 +e^{i2pi/3}",2*np.pi/3,2*np.pi/3),("v4 +e^{i0.7} ARBITRAIRE",2*np.pi/3,0.7)):
    res,uerr = run(kz,Nx,th,ex)
    print(f"{name:26} unit={uerr:.0e}  " + "; ".join(f"om={o:.4f} poids={w:.3f}" for o,w in res))
