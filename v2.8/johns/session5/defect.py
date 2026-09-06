"""Session 2 — ligne de defaut dans le reseau de Johns.
Supercellule Nx x Ny (periodique transverse), Bloch k_z le long de la corde.
3 modeles: D0 (bulk), Dcut (lien z coupe, Gamma=+-1 = controle), D3 (tresse cyclique
de 3 colonnes + rotation de polarisation 2pi/3 imposee par la geometrie du triplet).
Etat: 12 ports/noeud. U(kz) = C(kz) . blockdiag(S). Valeurs propres e^{i omega}.
Mode lie cherche SOUS le cone du bulk: omega < omega_min_bulk(kz) = kz/2, localise (IPR).
"""
import numpy as np
from scipy.linalg import eig

ports = [(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
pidx = {P:i for i,P in enumerate(ports)}
S1 = np.load("S_johns.npy")

def build(Nx, Ny, kz, model, core=None, gamma=1.0):
    Nn = Nx*Ny; dim = 12*Nn
    def gi(ix,iy,P): return 12*((ix%Nx)*Ny + (iy%Ny)) + pidx[P]
    Sb = np.zeros((dim,dim))
    for n in range(Nn):
        Sb[12*n:12*n+12, 12*n:12*n+12] = S1
    C = np.zeros((dim,dim), complex)
    th = 2*np.pi/3; R = np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
    if core is None: core = [(Nx//2,Ny//2),(Nx//2+1,Ny//2),(Nx//2,Ny//2+1)]
    coreset = set(core)
    for ix in range(Nx):
        for iy in range(Ny):
            # liens x et y (transverse, periodiques)
            for p in (1,2):
                C[gi(ix+1,iy,(0,-1,p)), gi(ix,iy,(0,+1,p))] = 1.0
                C[gi(ix-1,iy,(0,+1,p)), gi(ix,iy,(0,-1,p))] = 1.0
            for p in (0,2):
                C[gi(ix,iy+1,(1,-1,p)), gi(ix,iy,(1,+1,p))] = 1.0
                C[gi(ix,iy-1,(1,+1,p)), gi(ix,iy,(1,-1,p))] = 1.0
            # liens z (Bloch)
            here = (ix,iy)
            if model=="D3" and here in coreset:
                j = core.index(here); up = core[(j+1)%3]; dn = core[(j-1)%3]
                for a in range(2):        # pol d'arrivee
                    for b in range(2):    # pol de depart
                        C[gi(up[0],up[1],(2,-1,a)), gi(ix,iy,(2,+1,b))] = R[a,b]*np.exp(1j*kz)
                        C[gi(dn[0],dn[1],(2,+1,a)), gi(ix,iy,(2,-1,b))] = R.T[a,b]*np.exp(-1j*kz)
            elif model=="Dcut" and here==core[0]:
                for p in (0,1):           # lien z reflechi totalement (controle lacune)
                    C[gi(ix,iy,(2,-1,p)), gi(ix,iy,(2,+1,p))] = gamma
                    C[gi(ix,iy,(2,+1,p)), gi(ix,iy,(2,-1,p))] = gamma
            else:
                for p in (0,1):
                    C[gi(ix,iy,(2,-1,p)), gi(ix,iy,(2,+1,p))] = np.exp(1j*kz)
                    C[gi(ix,iy,(2,+1,p)), gi(ix,iy,(2,-1,p))] = np.exp(-1j*kz)
    return C@Sb, core

def spectrum(Nx,Ny,kz,model,gamma=1.0):
    U,core = build(Nx,Ny,kz,model,gamma=gamma)
    w,V = eig(U)
    om = np.angle(w)
    # localisation: poids sur les colonnes du coeur
    Nn = Nx*Ny
    wt = np.zeros(len(w))
    mask = np.zeros(12*Nn, bool)
    for (cx,cy) in core:
        mask[12*((cx%Nx)*Ny+(cy%Ny)):12*((cx%Nx)*Ny+(cy%Ny))+12] = True
    P2 = np.abs(V)**2; P2 /= P2.sum(axis=0, keepdims=True)
    wt = P2[mask].sum(axis=0)
    return om, wt, np.abs(w)

if __name__ == "__main__":
    Nx=Ny=11
    print(f"supercellule {Nx}x{Ny}; poids-de-coeur bulk attendu ~ {3/(Nx*Ny):.4f}")
    for kz in (0.3, 0.6, 1.0):
        for model in ("D0","D3"):
            om, wt, mod = spectrum(Nx,Ny,kz,model)
            # fenetre sous le cone: 0 < omega < kz/2 (exclut modes plats a ~0)
            sel = (om>0.02) & (om < kz/2 - 1e-6)
            print(f"kz={kz:.2f} {model}: modes sous cone (0.02<om<kz/2): {sel.sum():3d} | "
                  f"max poids-coeur sous cone: {wt[sel].max() if sel.any() else 0:.4f} | "
                  f"unitarite max| |w|-1 |: {np.abs(mod-1).max():.1e}")
            if sel.any():
                top = np.argsort(-wt*sel)[:4]
                for t in top:
                    if sel[t]: print(f"      om={om[t]:+.4f}  poids-coeur={wt[t]:.3f}")
