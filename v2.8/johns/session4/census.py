"""Session 4A/B — recensement spectral COMPLET (0<omega<pi), pas seulement sous le cone.
Motivation: en unites physiques hbar/Dt=0.487 MeV, omega_lattice=pi <-> E1=1.53 MeV et
omega_e=0.511 MeV <-> omega_lattice=1.0493. Un mode lie du defaut AU-DESSUS de la bande
du bulk (ou dans un gap interne) serait aussi legitime que sous le cone.
Dense 11x11 a kz donne: spectre D0 (continuum + bords de gaps) vs etats localises D3.
+ identification: nombre d'enroulement de phase autour du coeur pour les modes cles."""
import numpy as np
from scipy.linalg import eig
from track import buildU, core_mask
import scipy.sparse as sp
Nx=11; c=Nx//2
core=[(c,c),(c+1,c),(c,c+1)]
mask=core_mask(Nx,core)

def census(kz):
    res={}
    for model in ("D0","D3"):
        if model=="D0":
            # bulk: meme constructeur sans defaut -> utiliser core hors reseau? plus simple: defect.py build
            from defect import build
            U,_=build(Nx,Nx,kz,"D0")
        else:
            U=buildU(Nx,kz,core).toarray()
        w,V=eig(U); om=np.angle(w)
        P2=np.abs(V)**2; P2/=P2.sum(0,keepdims=True); wt=P2[mask].sum(0)
        res[model]=(om,wt,V)
    return res

def winding(V, i):
    """enroulement de phase d'une composante du champ sur un anneau r~3 autour du coeur."""
    ring=[]
    cx,cy=c+0.5,c+0.5
    import math
    pts=[]
    for ix in range(Nx):
        for iy in range(Nx):
            r=math.hypot(ix-cx,iy-cy)
            if 2.4<r<3.6: pts.append((math.atan2(iy-cy,ix-cx),ix,iy))
    pts.sort()
    # composante: somme des ports z-polarises en x du noeud (indice port (2,+1,0) = ?)
    ports=[(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
    pz=ports.index((2,+1,0))
    ph=[]
    for a,ix,iy in pts:
        amp=V[12*(ix*Nx+iy)+pz, i]
        ph.append(np.angle(amp))
    ph=np.unwrap(np.array(ph))
    return (ph[-1]-ph[0]+ (pts[0][0]+2*np.pi-pts[-1][0])*0 )/(2*np.pi)

for kz in (2.0, 0.6):
    r=census(kz)
    om0,wt0,_=r["D0"]; om3,wt3,V3=r["D3"]
    # gaps du bulk: trier les omega>0 du D0, chercher les trous > 0.05
    o=np.sort(om0[om0>1e-6])
    gaps=[(o[i],o[i+1]) for i in range(len(o)-1) if o[i+1]-o[i]>0.06]
    print(f"== kz={kz} ==")
    print(" gaps du bulk (>0.06):", [(round(a,3),round(b,3)) for a,b in gaps])
    loc=np.where((om3>0.02)&(wt3>0.15))[0]
    loc=loc[np.argsort(om3[loc])]
    print(" etats localises D3 (poids>0.15) sur tout (0,pi):")
    for i in loc:
        w_=winding(V3,i)
        print(f"   om={om3[i]:.5f}  (hbar_om={om3[i]*0.487:.4f} MeV)  poids={wt3[i]:.3f}  enroulement~{w_:+.2f}")
