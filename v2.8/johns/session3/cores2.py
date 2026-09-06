import numpy as np, scipy.sparse.linalg as spl
from track import buildU, core_mask
Nx=11; kz=2.0; c=Nx//2
geoms={
 "colineaire (0,0)(1,0)(2,0)": [(c,c),(c+1,c),(c+2,c)],
 "grand triangle (0,0)(2,0)(0,2)": [(c,c),(c+2,c),(c,c+2)],
}
for name,core in geoms.items():
    U=buildU(Nx,kz,core); mask=core_mask(Nx,core); found={}
    for s0 in (0.15,0.35,0.55,0.75,0.95):
        w,V=spl.eigs(U,k=6,sigma=np.exp(1j*s0),which='LM')
        om=np.angle(w); P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True); wt=P2[mask].sum(0)
        for o,t in zip(om,wt):
            if 0.02<o<kz/2-1e-6 and t>0.15: found[round(o,5)]=max(found.get(round(o,5),0),t)
    top=sorted(found.items(), key=lambda x:-x[1])[:3]
    print(f"{name:32}: " + ("; ".join(f"om={o:.5f} p={t:.3f}" for o,t in top) if top else "aucun"))
