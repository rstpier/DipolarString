"""La phase scalaire est-elle une boucle de Wilson (seul e^{3i phi} physique)?
Test: phi=0.7 vs phi=0.7+2pi/3 vs phi=0.7+4pi/3 -> spectres identiques attendus si oui.
Et les 3 secteurs discrets permis par la fermeture: phi=0, 2pi/3, 4pi/3."""
import numpy as np
from variants import run
kz=2.0; Nx=11; th=2*np.pi/3
for ex in (0.7, 0.7+2*np.pi/3, 0.7+4*np.pi/3):
    res,_=run(kz,Nx,th,ex)
    print(f"phi={ex:6.3f}: " + "; ".join(f"om={o:.5f}" for o,_ in res))
print("--- les trois secteurs de fermeture (e^{3i phi}=1):")
for ex in (0.0, 2*np.pi/3, 4*np.pi/3):
    res,_=run(kz,Nx,th,ex)
    print(f"phi={ex:6.3f}: " + "; ".join(f"om={o:.5f} poids={w:.3f}" for o,w in res))
