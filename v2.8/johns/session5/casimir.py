"""Session 5 — energie propre du defaut par decalage spectral (Casimir de connectivite).
Idee: la corde n'est pas une excitation du champ mais un changement de connectivite;
son cout energetique = (1/2)hbar * somme des decalages de frequences propres entre D3 et D0,
integre sur kz. C'est un INVARIANT SPECTRAL: la tension de ligne mu_s sort du spectre seul.
Declaration: l'interpretation suppose la quantification demi-quantum (1/2)hbar*omega du
reseau classique — pas d'autre ingredient.
Unites: hbar/Dt = 0.487 MeV ; longueur de cellule l1 = 810.4 fm.
DeltaE par cellule de longueur = 0.487 * [ 2*int_0^pi dkz/(2pi) * (1/2)*Sum(|om|_D3 - |om|_D0) ].
"""
import numpy as np
from scipy.linalg import eigvals
from track import buildU
from defect import build as build_dense
import scipy.sparse as sp

def dsum(Nx,kz):
    c=Nx//2; core=[(c,c),(c+1,c),(c,c+1)]
    U3=buildU(Nx,kz,core).toarray()
    U0,_=build_dense(Nx,Nx,kz,"D0")
    w3=np.abs(np.angle(eigvals(U3))); w0=np.abs(np.angle(eigvals(U0)))
    return 0.5*(np.sort(w3).sum()-np.sort(w0).sum())

Nx=9
ks=np.linspace(0.05,np.pi-0.05,9)
vals=[]
for kz in ks:
    d=dsum(Nx,kz); vals.append(d)
    print(f"kz={kz:.3f}  (1/2)Sum(D|om|) = {d:+.5f}")
I = np.trapezoid(vals,ks)/np.pi   # 2*int_0^pi dk/(2pi) = int/pi
E_cell = 0.487*I
print(f"\nintegrale (1/2)SumD|om| dk/pi = {I:+.5f}  ->  DeltaE par cellule = {E_cell*1000:+.2f} keV")
print(f"tension de ligne mu_s = {E_cell*1000:+.2f} keV / l1")
print(f"energie de la boucle fermee (3 cellules) = {3*E_cell*1000:+.2f} keV")
print(f"cibles: m_e = 511 keV ; fermeture E (10.7 keV) ; auto-energie v1 (2.19 keV)")
# controle taille a kz=pi/2
d9=dsum(9,np.pi/2); d11=dsum(11,np.pi/2)
print(f"\ncontrole taille a kz=pi/2: Nx=9 -> {d9:+.5f} ; Nx=11 -> {d11:+.5f} (ecart {abs(d11-d9):.1e})")
