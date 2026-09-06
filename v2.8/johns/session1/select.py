"""Selection parmi les 8: equivariance sous rotations (isotropie C4) et conservation du flux (C3)."""
import numpy as np
ports = [(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
idx = {P:i for i,P in enumerate(ports)}; N=12
def third(d,p): return 3-d-p
cands = np.load("candidates.npy")

# rotation de 90 deg autour de z: x->y, y->-x, z->z. Action sur axes: perm+signes.
# axe vecteur: ex->ey, ey->-ex, ez->ez
def rot_z(P):
    d,s,p = P
    mp = {0:(1,+1), 1:(0,-1), 2:(2,+1)}  # axe -> (axe', signe)
    d2,sd = mp[d]; p2,sp = mp[p]
    return (d2, s*sd, p2), sp   # le cote suit le signe de l'axe du bras; la pol porte le signe du champ
Rz = np.zeros((N,N))
for i,P in enumerate(ports):
    Q,sg = rot_z(P); Rz[idx[Q],i] = sg
def rot_x(P):
    d,s,p = P
    mp = {0:(0,+1), 1:(2,+1), 2:(1,-1)}  # y->z, z->-y
    d2,sd = mp[d]; p2,sp = mp[p]
    return (d2, s*sd, p2), sp
Rx = np.zeros((N,N))
for i,P in enumerate(ports):
    Q,sg = rot_x(P); Rx[idx[Q],i] = sg

# C3: conservation du flux magnetique: pour chaque t, vecteur circulation w_t:
# ports (d,s,p) avec {d,p}={les deux axes perp a t}; signe = orientation de la circulation:
# contribution de H_t d'une onde (d,s,p): signe de (d_hat x p_hat).t_hat * ... convention: w = s * eps(d,p,t)
import itertools
def eps(a,b,c):
    perm = (a,b,c)
    return {(0,1,2):1,(1,2,0):1,(2,0,1):1,(0,2,1):-1,(2,1,0):-1,(1,0,2):-1}.get(perm,0)
W = []
for t in range(3):
    w = np.zeros(N)
    for i,(d,s,p) in enumerate(ports):
        if d!=t and p!=t:
            w[i] = s*eps(d,p,t)
    W.append(w)

for k,S in enumerate(cands):
    e_rz = np.allclose(Rz@S@Rz.T, S)
    e_rx = np.allclose(Rx@S@Rx.T, S)
    fl = all(np.allclose(S.T@w, w) or np.allclose(S.T@w, -w) for w in W)
    flux_sign = [ 'w' if np.allclose(S.T@w,w) else ('-w' if np.allclose(S.T@w,-w) else 'x') for w in W]
    print(f"cand {k}: Rz-equiv {e_rz}  Rx-equiv {e_rx}  flux {flux_sign}")
