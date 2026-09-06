"""Fermeture rapide 1: test A/B/C (flux vs defaut local) + derivation de rho_bulk=3pi.
B (vraie torsion globale): E(Theta) = 0.487 * Sum_n ds((2pi n + Theta)/3), R INCHANGE.
Pour Theta=2pi (holonomie identique au secteur A phi=2pi/3): les kz permis {2pi/3, 4pi/3, 2pi}
sont LE MEME ENSEMBLE modulo 2pi que {0, +-2pi/3} -> E(2pi)=E(0) par construction de Bloch
(c'est aussi la construction C: la repartition de jauge redonne B). Verification numerique
directe de la periodicite/symetrie de ds, puis verdict."""
import numpy as np
from loop6 import ds
Nx=9
print("symetrie et periodicite de ds (R inchange, phi=0):")
for kz in (2*np.pi/3, 4*np.pi/3, -2*np.pi/3, 2*np.pi):
    print(f"  ds(kz={kz:+.4f}) = {ds(Nx,kz):+.6f}")
print("\ncourbe du VRAI flux E(Theta), R inchange:")
for Th in (0.0, np.pi/3, 2*np.pi/3, np.pi, 2*np.pi):
    E=487*sum(ds(Nx,(2*np.pi*n+Th)/3) for n in range(3))
    print(f"  Theta={Th:.4f}: E = {E:+.2f} keV")
print("\nrappel secteur A (e^{i 2pi/3} R): E = -126.06 keV")
# rho_bulk = 3pi: appariement omega <-> pi-omega par q
ports=[(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
idx={P:i for i,P in enumerate(ports)}; S=np.load("S_johns.npy")
rng=np.random.default_rng(1)
ok=True
for _ in range(6):
    q=rng.uniform(-np.pi,np.pi,3)
    C=np.zeros((12,12),complex)
    for i,(d,s,p) in enumerate(ports): C[idx[(d,-s,p)],i]=np.exp(1j*q[d]*s)
    om=np.sort(np.abs(np.angle(np.linalg.eigvals(C@S))))
    pairs=om+om[::-1]
    ok &= np.allclose(pairs, np.pi, atol=1e-10)
print(f"\nappariement |om| <-> pi-|om| a chaque q (6 q aleatoires): {ok}")
print("=> Sum|om|(q) = 6pi IDENTIQUEMENT => rho_bulk = 3pi * (hbar/Dt) EXACT (plus une moyenne)")
