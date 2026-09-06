"""Reconstruction du noeud SCN de Johns depuis les contraintes (C1-C5 style),
puis dispersion de Bloch du reseau homogene.
Ports: (d,s,p) d=axe de bras 0..2, s=+-1 cote, p=polarisation (p!=d). 12 ports.
Regle structurelle (orbites): b(d,s,p) couple a
  +1/2 : (t,+,p) et (t,-,p)   [meme polarisation, bras selon t = 3e axe]
  +-1/2: (p,+,d) et (p,-,d)   [ports croises, un + et un -]
Le motif de signes est determine par S=S^T + unitarite (recherche exhaustive 2^12).
"""
import numpy as np, itertools

ports = [(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
idx = {P:i for i,P in enumerate(ports)}
N = len(ports)
def third(d,p): return 3-d-p

def build_S(minus_choice):
    """minus_choice[i] in {0,1}: 0 -> minus sur (p,+,d), 1 -> minus sur (p,-,d)"""
    S = np.zeros((N,N))
    for i,(d,s,p) in enumerate(ports):
        t = third(d,p)
        S[i, idx[(t,+1,p)]] = 0.5
        S[i, idx[(t,-1,p)]] = 0.5
        sgnP = -0.5 if minus_choice[i]==0 else 0.5
        sgnM = -0.5 if minus_choice[i]==1 else 0.5
        S[i, idx[(p,+1,d)]] = sgnP
        S[i, idx[(p,-1,d)]] = sgnM
    return S

sols = []
for mc in itertools.product((0,1), repeat=N):
    S = build_S(mc)
    if not np.allclose(S, S.T): continue
    if not np.allclose(S@S.T, np.eye(N)): continue
    sols.append((mc, S))
print(f"solutions symetriques+unitaires: {len(sols)}")

# inversion spatiale: (d,s,p) -> (d,-s,p) avec signe - sur la polarisation (E impair)
Pinv = np.zeros((N,N))
for i,(d,s,p) in enumerate(ports):
    Pinv[idx[(d,-s,p)], i] = -1.0
kept = []
for mc,S in sols:
    inv_ok = np.allclose(Pinv@S@Pinv.T, S)
    # conservation de charge par polarisation (C2): somme des b pol p = somme des a pol p
    ok_charge = True
    for p in range(3):
        u = np.array([1.0 if pp==p else 0.0 for (_,_,pp) in ports])
        if not np.allclose(S.T@u, u): ok_charge = False
    kept.append((mc,S,inv_ok,ok_charge))
    print("mc:", "".join(map(str,mc)), "| inversion-invariant:", inv_ok, "| charge:", ok_charge)
np.save("candidates.npy", np.array([S for _,S,_,_ in kept]))
