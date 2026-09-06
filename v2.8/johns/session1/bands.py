"""Structure de bandes complete du reseau de Johns (cand7) + inventaire des branches."""
import numpy as np
ports = [(d,s,p) for d in range(3) for s in (+1,-1) for p in range(3) if p!=d]
idx = {P:i for i,P in enumerate(ports)}; N=12
S = np.load("candidates.npy")[7]
np.save("S_johns.npy", S)

def M(q):
    C = np.zeros((N,N), complex)
    for i,(d,s,p) in enumerate(ports):
        C[idx[(d,-s,p)], i] = np.exp(1j*q[d]*s)
    return C@S

def om_all(q):
    return np.sort(np.angle(np.linalg.eigvals(M(q))))

# inventaire a q petit selon x
q = np.array([0.05,0,0])
om = om_all(q)
print("q=0.05 x  :", np.round(om,5))
# classement: zeros plats, acoustiques +-q/2, pi-modes
for lab,qv in [("x",[1,0,0]),("diag",[1,1,1]/np.sqrt(3)),("xy",[1,1,0]/np.sqrt(2))]:
    print(f"-- direction {lab}: v des branches acoustiques (fit sur q=0.01..0.2)")
    qs = np.linspace(0.01,0.2,8)
    acc = []
    for qq in qs:
        om = om_all(np.array(qv)*qq)
        pos = om[(om>1e-8)&(om<1.0)]
        acc.append(sorted(pos))
    n = min(len(a) for a in acc)
    acc = np.array([a[:n] for a in acc])
    for b in range(n):
        v = np.polyfit(qs, acc[:,b], 1)[0]
        print(f"   branche {b}: v = {v:.5f}")
# combien de modes exactement plats a omega=0 et omega=pi, a q generique?
q = np.array([0.31,0.17,0.53])
om = om_all(q)
print("q generique:", np.round(om,4))
print("modes ~0 :", np.sum(np.abs(om)<1e-9), " modes ~pi:", np.sum(np.abs(np.abs(om)-np.pi)<1e-9))
