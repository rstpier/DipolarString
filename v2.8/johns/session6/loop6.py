"""Session 6 — la boucle fermee reelle par compactification exacte.
La tresse a periode 3 en z: une boite z-periodique de Lz=3 cellules AVEC le defaut EST la
boucle fermee (la corde se referme a travers la direction compacte; circonference 3*l1 =
2*pi*lambda_C_barre exactement). Decomposition de Bloch exacte: les kz permis sont
{0, +2pi/3, -2pi/3}, donc
    E_boucle(fermee) = 0.487 * [ ds(0) + ds(+2pi/3) + ds(-2pi/3) ]
sans aucune matrice 3D. Geometrie torique: E_courbure et E_coins ABSENTS par construction
(declare); on isole E_fermeture = E_boucle(fermee) - 3*mu_s.
Secteurs de flux phi in {0, 2pi/3, 4pi/3} (e^{3i phi}=1, permis pour la boucle FERMEE):
multiplicite spectrale candidate."""
import numpy as np
from scipy.linalg import eigvals
import scipy.sparse as sp
ports=[(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
pidx={P:i for i,P in enumerate(ports)}
S1=np.load("S_johns.npy")

def U_of(Nx,kz,model,phi=0.0):
    Nn=Nx*Nx
    def gi(ix,iy,P): return 12*((ix%Nx)*Nx+(iy%Nx))+pidx[P]
    rows=[];cols=[];vals=[]
    th=2*np.pi/3
    R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])*np.exp(1j*phi)
    c=Nx//2; core=[(c,c),(c+1,c),(c,c+1)]; cs=set(core)
    A=lambda r,cc,v:(rows.append(r),cols.append(cc),vals.append(v))
    for ix in range(Nx):
        for iy in range(Nx):
            for p in (1,2):
                A(gi(ix+1,iy,(0,-1,p)),gi(ix,iy,(0,+1,p)),1); A(gi(ix-1,iy,(0,+1,p)),gi(ix,iy,(0,-1,p)),1)
            for p in (0,2):
                A(gi(ix,iy+1,(1,-1,p)),gi(ix,iy,(1,+1,p)),1); A(gi(ix,iy-1,(1,+1,p)),gi(ix,iy,(1,-1,p)),1)
            if model=="D3" and (ix,iy) in cs:
                j=core.index((ix,iy)); up=core[(j+1)%3]; dn=core[(j-1)%3]
                for a in range(2):
                    for b in range(2):
                        A(gi(up[0],up[1],(2,-1,a)),gi(ix,iy,(2,+1,b)),R[a,b]*np.exp(1j*kz))
                        A(gi(dn[0],dn[1],(2,+1,a)),gi(ix,iy,(2,-1,b)),np.conj(R.T)[a,b]*np.exp(-1j*kz))
            else:
                for p in (0,1):
                    A(gi(ix,iy,(2,-1,p)),gi(ix,iy,(2,+1,p)),np.exp(1j*kz))
                    A(gi(ix,iy,(2,+1,p)),gi(ix,iy,(2,-1,p)),np.exp(-1j*kz))
    C=sp.csr_matrix((vals,(rows,cols)),shape=(12*Nn,12*Nn))
    return (C@sp.block_diag([S1]*Nn,format='csr')).toarray()

def ds(Nx,kz,phi=0.0):
    w3=np.abs(np.angle(eigvals(U_of(Nx,kz,"D3",phi))))
    w0=np.abs(np.angle(eigvals(U_of(Nx,kz,"D0"))))
    return 0.5*(np.sort(w3).sum()-np.sort(w0).sum())

Nx=9
mu_s=12.77  # keV/cellule (valeur durcie, bords inclus, taille corrigee)
print("secteur phi=0:")
d0  = ds(Nx,0.0);          print(f"  ds(kz=0)      = {d0:+.6f}")
dp  = ds(Nx, 2*np.pi/3);   print(f"  ds(kz=+2pi/3) = {dp:+.6f}")
dm  = ds(Nx,-2*np.pi/3);   print(f"  ds(kz=-2pi/3) = {dm:+.6f}")
E0 = 487*(d0+dp+dm)
print(f"  E_boucle(fermee, phi=0) = {E0:.2f} keV ; 3*mu_s = {3*mu_s:.2f} keV ; E_fermeture = {E0-3*mu_s:+.2f} keV")
for phi in (2*np.pi/3, 4*np.pi/3):
    d0p=ds(Nx,0.0,phi); dpp=ds(Nx,2*np.pi/3,phi); dmp=ds(Nx,-2*np.pi/3,phi)
    E=487*(d0p+dpp+dmp)
    print(f"secteur phi={phi:.4f}: ds = {d0p:+.6f}, {dpp:+.6f}, {dmp:+.6f} -> E = {E:.2f} keV")
