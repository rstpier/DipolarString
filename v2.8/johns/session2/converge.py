"""Convergence en taille (shift-invert creux autour de om*=0.3325) + profil de decroissance transverse."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
ports = [(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
pidx = {P:i for i,P in enumerate(ports)}
S1 = np.load("S_johns.npy")
def buildU(Nx,kz):
    Nn=Nx*Nx; dim=12*Nn
    def gi(ix,iy,P): return 12*((ix%Nx)*Nx+(iy%Nx))+pidx[P]
    rows=[];cols=[];vals=[]
    th=2*np.pi/3; R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
    core=[(Nx//2,Nx//2),(Nx//2+1,Nx//2),(Nx//2,Nx//2+1)]; cs=set(core)
    def C_add(r,c,v): rows.append(r);cols.append(c);vals.append(v)
    for ix in range(Nx):
        for iy in range(Nx):
            for p in (1,2):
                C_add(gi(ix+1,iy,(0,-1,p)),gi(ix,iy,(0,+1,p)),1)
                C_add(gi(ix-1,iy,(0,+1,p)),gi(ix,iy,(0,-1,p)),1)
            for p in (0,2):
                C_add(gi(ix,iy+1,(1,-1,p)),gi(ix,iy,(1,+1,p)),1)
                C_add(gi(ix,iy-1,(1,+1,p)),gi(ix,iy,(1,-1,p)),1)
            if (ix,iy) in cs:
                j=core.index((ix,iy)); up=core[(j+1)%3]; dn=core[(j-1)%3]
                for a in range(2):
                    for b in range(2):
                        C_add(gi(up[0],up[1],(2,-1,a)),gi(ix,iy,(2,+1,b)),R[a,b]*np.exp(1j*kz))
                        C_add(gi(dn[0],dn[1],(2,+1,a)),gi(ix,iy,(2,-1,b)),R.T[a,b]*np.exp(-1j*kz))
            else:
                for p in (0,1):
                    C_add(gi(ix,iy,(2,-1,p)),gi(ix,iy,(2,+1,p)),np.exp(1j*kz))
                    C_add(gi(ix,iy,(2,+1,p)),gi(ix,iy,(2,-1,p)),np.exp(-1j*kz))
    C=sp.csr_matrix((vals,(rows,cols)),shape=(dim,dim))
    Sb=sp.block_diag([S1]*Nn,format='csr')
    return (C@Sb).tocsc(), core
kz=2.0
for Nx in (11,15,21):
    U,core=buildU(Nx,kz)
    w,V=spl.eigs(U,k=6,sigma=np.exp(1j*0.3325),which='LM')
    om=np.angle(w)
    # poids coeur
    dim=12*Nx*Nx; mask=np.zeros(dim,bool)
    for (cx,cy) in core:
        b=12*((cx%Nx)*Nx+(cy%Nx)); mask[b:b+12]=True
    P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True); wt=P2[mask].sum(0)
    i=np.argmax(wt)
    print(f"Nx={Nx:2d}: om_lie={om[i]:.5f}  poids-coeur={wt[i]:.3f}")
    if Nx==21:
        # profil radial du mode
        psi=np.abs(V[:,i])**2
        prof={}
        for ix in range(Nx):
            for iy in range(Nx):
                r=np.hypot(ix-Nx//2, iy-Nx//2)
                prof.setdefault(round(r),0)
                prof[round(r)] += psi[12*(ix*Nx+iy):12*(ix*Nx+iy)+12].sum()
        rs=sorted(prof)[1:8]
        lp=[np.log(prof[r]) for r in rs]
        kap=-np.polyfit(rs,lp,1)[0]
        print("   profil ln P(r):", {r:round(np.log(prof[r]),2) for r in sorted(prof)[:8]})
        print(f"   kappa mesure = {kap:.3f}   (continuum attendu sqrt(kz^2-4om^2) = {np.sqrt(kz**2-4*om[i]**2):.3f})")
