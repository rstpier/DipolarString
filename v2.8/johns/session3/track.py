"""Session 3 (1)-(3): suivi de branche par continuite de recouvrement + v_g reel.
Ancre: kz=2.0, om*=0.33246 (converge). Pas Δkz=0.05, appariement par |<psi_prev|psi_new>|,
sous-espace k=8 autour de sigma=e^{i om_prev}. Nx=15. Arret si le mode rejoint le continuum."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
ports=[(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
pidx={P:i for i,P in enumerate(ports)}
S1=np.load("S_johns.npy")

def buildU(Nx,kz,core):
    Nn=Nx*Nx
    def gi(ix,iy,P): return 12*((ix%Nx)*Nx+(iy%Nx))+pidx[P]
    rows=[];cols=[];vals=[]
    th=2*np.pi/3; R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
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

def core_mask(Nx,core):
    m=np.zeros(12*Nx*Nx,bool)
    for cx,cy in core:
        b=12*((cx%Nx)*Nx+(cy%Nx)); m[b:b+12]=True
    return m

def track(Nx=15, kz0=2.0, om0=0.33246, dks=None):
    core=[(Nx//2,Nx//2),(Nx//2+1,Nx//2),(Nx//2,Nx//2+1)]
    mask=core_mask(Nx,core)
    out=[]
    for sign in (+1,-1):
        kz, om = kz0, om0; psi=None
        while 0.15 < kz < np.pi:
            U=buildU(Nx,kz,core)
            w,V=spl.eigs(U,k=8,sigma=np.exp(1j*om),which='LM')
            if psi is None:
                # depart: max poids coeur parmi les 8
                P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True)
                i=np.argmax(P2[mask].sum(0))
            else:
                i=np.argmax(np.abs(psi.conj()@V))
            psi=V[:,i]/np.linalg.norm(V[:,i])
            om=np.angle(w[i])
            wt=(np.abs(psi)**2)[mask].sum()/ (np.abs(psi)**2).sum()
            out.append((kz,om,wt))
            if om > kz/2-0.01 or wt < 0.05: break
            kz += sign*0.1
        psi=None
    out=sorted(set((round(k,3),round(o,5),round(t,3)) for k,o,t in out))
    return out

if __name__=='__main__':
    res=track()
    print(f"{'kz':>6} {'omega':>9} {'poids':>7} {'om/(kz/2)':>9}")
    for k,o,t in res: print(f"{k:6.2f} {o:9.5f} {t:7.3f} {o/(k/2):9.3f}")
    np.save("branch.npy", np.array(res))
    # v_g par differences centrees sur la partie liee
    arr=np.array([r for r in res if r[1] < r[0]/2-0.01 and r[2]>0.2])
    if len(arr)>4:
        k,o=arr[:,0],arr[:,1]
        vg=np.gradient(o,k)
        print("\n v_g = domega/dkz (branche suivie):")
        for i in range(len(k)): print(f"  kz={k[i]:.2f}  om={o[i]:.5f}  v_g={vg[i]:+.4f}  v_g/c0={vg[i]/0.5:+.3f}")
