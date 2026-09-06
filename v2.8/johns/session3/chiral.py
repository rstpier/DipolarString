"""Test de coherence chirale: cycle inverse + R(-2pi/3) doit redonner le spectre de reference
(anti-corde = miroir complet). Cycle inverse + R(+2pi/3) = objet depareille (deja vu: 0.467)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
ports=[(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
pidx={P:i for i,P in enumerate(ports)}
S1=np.load("S_johns.npy")
def buildU(Nx,kz,core,theta):
    Nn=Nx*Nx
    def gi(ix,iy,P): return 12*((ix%Nx)*Nx+(iy%Nx))+pidx[P]
    rows=[];cols=[];vals=[]
    R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    cs=set(core)
    A=lambda r,c,v:(rows.append(r),cols.append(c),vals.append(v))
    for ix in range(Nx):
        for iy in range(Nx):
            for p in (1,2):
                A(gi(ix+1,iy,(0,-1,p)),gi(ix,iy,(0,+1,p)),1); A(gi(ix-1,iy,(0,+1,p)),gi(ix,iy,(0,-1,p)),1)
            for p in (0,2):
                A(gi(ix,iy+1,(1,-1,p)),gi(ix,iy,(1,+1,p)),1); A(gi(ix,iy-1,(1,+1,p)),gi(ix,iy,(1,-1,p)),1)
            if (ix,iy) in cs:
                j=core.index((ix,iy)); up=core[(j+1)%3]; dn=core[(j-1)%3]
                for a in range(2):
                    for b in range(2):
                        A(gi(up[0],up[1],(2,-1,a)),gi(ix,iy,(2,+1,b)),R[a,b]*np.exp(1j*kz))
                        A(gi(dn[0],dn[1],(2,+1,a)),gi(ix,iy,(2,-1,b)),R.T[a,b]*np.exp(-1j*kz))
            else:
                for p in (0,1):
                    A(gi(ix,iy,(2,-1,p)),gi(ix,iy,(2,+1,p)),np.exp(1j*kz))
                    A(gi(ix,iy,(2,+1,p)),gi(ix,iy,(2,-1,p)),np.exp(-1j*kz))
    C=sp.csr_matrix((vals,(rows,cols)),shape=(12*Nn,12*Nn))
    return (C@sp.block_diag([S1]*Nn,format='csr')).tocsc()
Nx=11; kz=2.0; c=Nx//2
ref=[(c,c),(c+1,c),(c,c+1)]; inv=[(c,c),(c,c+1),(c+1,c)]
for name,core,th in (("ref cycle+R(+)",ref,2*np.pi/3),
                     ("cycle inverse + R(-)",inv,-2*np.pi/3),
                     ("cycle inverse + R(+) [depareille]",inv,2*np.pi/3)):
    U=buildU(Nx,kz,core,th)
    from track import core_mask
    mask=core_mask(Nx,core); found={}
    for s0 in (0.15,0.35,0.55):
        w,V=spl.eigs(U,k=6,sigma=np.exp(1j*s0),which='LM')
        om=np.angle(w); P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True); wt=P2[mask].sum(0)
        for o,t in zip(om,wt):
            if 0.02<o<kz/2-1e-6 and t>0.15: found[round(o,5)]=max(found.get(round(o,5),0),t)
    top=sorted(found.items(), key=lambda x:-x[1])[:2]
    print(f"{name:36}: " + "; ".join(f"om={o:.5f} p={t:.3f}" for o,t in top))
