"""Dispersion de Bloch du reseau SCN homogene pour les 2 candidats isotropes."""
import numpy as np
ports = [(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
idx = {P:i for i,P in enumerate(ports)}; N=12
cands = np.load("candidates.npy")
S0, S7 = cands[0], cands[7]

def M(q, S):
    C = np.zeros((N,N), complex)
    for i,(d,s,p) in enumerate(ports):
        C[idx[(d,-s,p)], i] = np.exp(1j*q[d]*s)   # b sort cote s, entre chez le voisin s*d_hat
    return C@S

def branches(q, S):
    lam = np.linalg.eigvals(M(q,S))
    om = np.angle(lam)   # omega*Dt, Dt=Dl/c_lien=1
    return np.sort(om)

for name,S in (("cand0",S0),("cand7",S7)):
    print(f"== {name} ==")
    for qx in (0.02, 0.05, 0.1, 0.2):
        om = branches(np.array([qx,0,0]), S)
        pos = om[om>1e-9]
        print(f" q=({qx:.2f},0,0)  omegas>0: {np.round(pos,5)}  -> v=omega/q: {np.round(pos/qx,4)}")
    # direction diagonale
    for qq in (0.05, 0.2):
        q = np.array([qq,qq,qq])/np.sqrt(3)
        om = branches(q,S); pos = om[om>1e-9]
        print(f" q diag |q|={qq:.2f}      omegas>0: {np.round(pos,5)}  -> v: {np.round(pos/qq,4)}")
