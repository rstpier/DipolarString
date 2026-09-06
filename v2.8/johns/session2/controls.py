import numpy as np
from defect import build
from scipy.linalg import eig
kz=2.0; Nx=11
def analyse(model, mod_R=None, gamma=1.0):
    import defect
    U,core = build(Nx,Nx,kz,model,gamma=gamma)
    w,V = eig(U); om=np.angle(w)
    Nn=Nx*Nx; mask=np.zeros(12*Nn,bool)
    for (cx,cy) in core: mask[12*((cx%Nx)*Nx+(cy%Nx)):12*((cx%Nx)*Nx+(cy%Nx))+12]=True
    P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True); wt=P2[mask].sum(0)
    sel=(om>0.02)&(om<kz/2-1e-6)
    if not sel.any(): return "aucun mode sous cone"
    o=np.argsort(-(wt*sel))[:3]
    return "; ".join(f"om={om[i]:.4f} poids={wt[i]:.3f}" for i in o if sel[i])
print("Dcut Gamma=+1 :", analyse("Dcut", gamma=+1.0))
print("Dcut Gamma=-1 :", analyse("Dcut", gamma=-1.0))
