"""Verification EXACTE des gaps du bulk: bandes 12x12 de Bloch balayees sur k_perp
(la supercellule 11x11 discretise k_perp -> faux gaps de taille finie possibles)."""
import numpy as np
ports=[(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
idx={P:i for i,P in enumerate(ports)}; N=12
S=np.load("S_johns.npy")
def om_all(q):
    C=np.zeros((N,N),complex)
    for i,(d,s,p) in enumerate(ports):
        C[idx[(d,-s,p)],i]=np.exp(1j*q[d]*s)
    return np.angle(np.linalg.eigvals(C@S))
for kz in (0.6, 2.0):
    oms=[]
    g=np.linspace(-np.pi,np.pi,61)
    for kx in g:
        for ky in g:
            oms.append(om_all(np.array([kx,ky,kz])))
    o=np.sort(np.concatenate(oms)); o=o[o>1e-6]
    gaps=[]
    for i in range(len(o)-1):
        if o[i+1]-o[i]>0.03: gaps.append((o[i],o[i+1]))
    print(f"kz={kz}: bande couvre [{o.min():.4f}, {o.max():.4f}]")
    print("  vrais gaps (>0.03):", [(round(a,4),round(b,4)) for a,b in gaps] or "AUCUN")
    # position de omega_e=1.0493 et des etats localises trouves
    for target,name in ((1.0493,"omega_e"),(0.33246,"ancre"),(1.4959,"famille mediane"),(2.5853,"famille haute")):
        inband = ((o>target-0.015)&(o<target+0.015)).sum()
        print(f"  {name} ({target:.4f}): densite locale du bulk ~{inband} etats/0.03 -> {'DANS la bande' if inband>0 else 'HORS bande (gap/stopband)'}")
