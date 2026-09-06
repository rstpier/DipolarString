#!/usr/bin/env python3
"""
P2 master — Emergent relativistic dynamics of autonomous sub-gap breathers
Reproduit tous les resultats de l'annexe. Usage:
    python3 p2_master.py rest        # stabilite au repos + E0, om_b0 vs maille h et dt
    python3 p2_master.py boost      # E, P, m_inv, om_b(v) pour v=0.2,0.4 ; h=1 et 0.5 ; 2 breathers
    python3 p2_master.py force      # m_F = J/D(gamma v), +/-lambda, temoin lambda=0, rayonnement
    python3 p2_master.py force25    # cas m_F a h=0.25 (3e point de convergence, long)
    python3 p2_master.py open       # CONTROLE OUVERT (section 5): soliton pompe a phase imposee
                                    # sur chaine gappee, X dynamique -> m_dress ~ 0.15 U0
    python3 p2_master.py all

MODELE. Chaine sine-Gordon, pas de maille h:
    d2q_i/dt2 = (q_{i+1} - 2 q_i + q_{i-1})/h^2 - sin(q_i)         (c = Omega0 = 1)
Breather continu exact (omega^2 + eta^2 = 1):
    q = 4 atan[(eta/omega) sin(omega tau)/cosh(eta xi)],  xi = g(x-x0-vt), tau = g(t-v(x-x0))
Eponges: gamma(x) quadratique sur n_sp mailles a chaque bord (absorbant).

DEFINITIONS (mesurees sur la fenetre core, hors eponges):
    E     = sum_i h [ p_i^2/2 + (1-cos q_i) ] + sum_i h [ (q_{i+1}-q_i)/h ]^2 / 2
    P     = - sum_i h p_i (q_{i+1}-q_{i-1})/(2h)          (pseudo-impulsion canonique)
    X_E   = centroide d'energie ; v = pente de X_E(t) (fit lineaire)
    m_inv = sqrt(E^2 - P^2)
    om_b  = frequence des zeros ascendants de q interpole AU CENTROIDE EXACT
            (l'echantillonnage au site entier injecte une phase parasite g*v*om)
    m_F   = J / Delta(gamma v), J = integrale de la force externe mesuree
            F_ext(t) = -lam(t) sum_i h (x_i - x_r)' ... = -lam(t) sum_i h (1-cos q_i)
    rayonnement = energie dans les zones |x - X_E| > W_rad (loin du coeur)
Toutes les moyennes sur fenetres explicites; erreurs = dispersion entre demi-fenetres.
"""
import sys
import numpy as np

def make(h, Ncells):
    x = np.arange(Ncells)*h
    nsp = int(300/h)
    gam = np.zeros(Ncells)
    gam[:nsp] = 0.3*(1-np.arange(nsp)/nsp)**2
    gam[-nsp:] = gam[:nsp][::-1]
    return x, gam, nsp

def init(x, v, x0, om):
    eta = np.sqrt(1-om*om); A = eta/om
    g = 1/np.sqrt(1-v*v)
    xi = np.clip(g*(x-x0), -300, 300); tau = -g*v*(x-x0)
    ch = np.cosh(eta*xi); u = A*np.sin(om*tau)/ch
    q = 4*np.arctan(u)
    dudt = A*om*np.cos(om*tau)*g/ch + A*np.sin(om*tau)*np.tanh(eta*xi)/ch*eta*g*v
    return q, 4*dudt/(1+u*u)

def lap(q, h):
    l = np.empty_like(q)
    l[1:-1] = (q[2:]-2*q[1:-1]+q[:-2])/(h*h)
    l[0] = (q[1]-2*q[0])/(h*h); l[-1] = (q[-2]-2*q[-1])/(h*h)
    return l

def run(h=1.0, Ncells=4000, dt=0.02, T=500.0, v=0.0, om=0.95, lam0=0.0,
        t0f=80.0, Tr=40.0, Tpl=110.0, W_rad=150.0):
    x, gam, nsp = make(h, Ncells)
    x0 = x[Ncells//4] if v >= 0 else x[3*Ncells//4]
    q, p = init(x, v, x0, om)
    xr = x0
    t1f = t0f+Tr+Tpl
    def s5(u): return 10*u**3-15*u**4+6*u**5
    def lamt(t):
        if lam0 == 0.0 or t < t0f: return 0.0
        if t < t0f+Tr: return lam0*s5((t-t0f)/Tr)
        if t < t1f: return lam0
        if t < t1f+Tr: return lam0*(1-s5((t-t1f)/Tr))
        return 0.0
    core = slice(nsp, Ncells-nsp)
    xc = x[core]
    nst = int(T/dt); t = 0.0; J = 0.0; rec = []
    for i in range(nst):
        lm = lamt(t); lm2 = lamt(t+dt/2); lm3 = lamt(t+dt)
        def acc(q_, p_, l_): return lap(q_, h) - np.sin(q_) - l_*(x-xr)*np.sin(q_) - gam*p_
        k1q = p;            k1p = acc(q, p, lm)
        k2q = p+dt/2*k1p;   k2p = acc(q+dt/2*k1q, p+dt/2*k1p, lm2)
        k3q = p+dt/2*k2p;   k3p = acc(q+dt/2*k2q, p+dt/2*k2p, lm2)
        k4q = p+dt*k3p;     k4p = acc(q+dt*k3q, p+dt*k3p, lm3)
        q += dt/6*(k1q+2*k2q+2*k3q+k4q); p += dt/6*(k1p+2*k2p+2*k3p+k4p)
        J += dt*(-lamt(t+dt/2))*np.sum(1-np.cos(q[core]))*h
        t += dt
        if i % max(1, int(0.1/dt)) == 0:
            qc = q[core]; pc = p[core]
            e = 0.5*pc*pc+(1-np.cos(qc)); e[:-1] += 0.5*((qc[1:]-qc[:-1])/h)**2
            e *= h
            E = e.sum(); Xc = (xc*e).sum()/E
            P = -np.sum(pc[1:-1]*(qc[2:]-qc[:-2]))/2
            far = np.abs(xc-Xc) > W_rad
            Erad = e[far].sum()
            j = (Xc-xc[0])/h; j0 = int(j); fr = j-j0
            qcen = qc[j0]*(1-fr)+qc[min(j0+1, len(qc)-1)]*fr
            rec.append((t, E, P, Xc, qcen, Erad))
    return np.array(rec), J

def freq(t, s):
    z = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    if len(z) < 3: return np.nan
    tz = t[z]+(t[z+1]-t[z])*(-s[z])/(s[z+1]-s[z])
    d = np.diff(tz)
    return 2*np.pi/d.mean(), 2*np.pi*d.std()/(d.mean()**2*np.sqrt(len(d)))

def halves(y):
    n = len(y)//2
    return y.mean(), abs(y[:n].mean()-y[n:].mean())/2

def suite_rest():
    print("== REST: stabilite, E0, om_b0 (theorie continue: E0=16*eta, om_b0=om) ==")
    print(f"{'h':>5} {'dt':>6} {'om':>5} {'E0_mes':>9} {'E0_th':>8} {'dE/E':>9} {'om_b0':>8} {'d_om/om':>9} {'derive':>9} {'E_rad/E':>9}")
    for om in (0.95, 0.90):
        eta = np.sqrt(1-om*om); E0t = 16*eta
        for h, dt in ((1.0, 0.02), (0.5, 0.01), (0.25, 0.005), (1.0, 0.01)):
            N = int(4000/h)
            r, _ = run(h=h, Ncells=N, dt=dt, T=300.0, v=0.0, om=om)
            t_, E_, P_, X_, qc_, Er_ = r.T
            m = t_ > 40
            E0, dE = halves(E_[m])
            ob, dob = freq(t_[m], qc_[m])
            v2 = np.polyfit(t_[m], X_[m], 1)[0]
            print(f"{h:5.2f} {dt:6.3f} {om:5.2f} {E0:9.4f} {E0t:8.4f} {dE/E0:9.1e} {ob:8.5f} {(ob-om)/om:9.2e} {v2:9.2e} {Er_[m].mean()/E0:9.1e}")

def suite_boost():
    print("== BOOST: E/gam, m_inv, om_b*gam/om_b0 pour v=0.2,0.4 ==")
    print(f"{'h':>5} {'om':>5} {'v':>5} {'E/g/E0':>8} {'m_inv/E0':>9} {'obg/ob0':>8} {'derive':>9} {'E_rad/E':>9}")
    for om in (0.95, 0.90):
        eta = np.sqrt(1-om*om)
        for h in (1.0, 0.5):
            N = int(4000/h); dt = 0.02*h
            r0, _ = run(h=h, Ncells=N, dt=dt, T=200.0, v=0.0, om=om)
            m0 = r0[:, 0] > 40
            E0 = r0[m0, 1].mean(); ob0, _ = freq(r0[m0, 0], r0[m0, 4])
            for v in (0.2, 0.4):
                r, _ = run(h=h, Ncells=N, dt=dt, T=400.0, v=v, om=om)
                t_, E_, P_, X_, qc_, Er_ = r.T
                m = (t_ > 60) & (X_ > 350) & (X_ < N*h-350)
                th = t_[m]; Xh = X_[m]; n2 = len(th)//2
                v1 = np.polyfit(th[:n2], Xh[:n2], 1)[0]; v2 = np.polyfit(th[n2:], Xh[n2:], 1)[0]
                vm = np.polyfit(th, Xh, 1)[0]            # vitesse MESUREE
                g = 1/np.sqrt(1-vm*vm)                    # gamma mesure partout
                E, dE = halves(E_[m]); P, dP = halves(P_[m])
                minv = np.sqrt(max(E*E-P*P, 0))
                dminv = (E*dE+P*dP)/max(minv, 1e-12)
                ob, dob = freq(t_[m], qc_[m])
                print(f"{h:5.2f} {om:5.2f} {vm:6.4f} {E/g/E0:8.4f}({dE/g/E0:.4f}) {minv/E0:9.4f}({dminv/E0:.4f}) {ob/(ob0*np.sqrt(1-vm*vm)) if ob0 else np.nan:8.4f} {(v2-v1)/v1:9.2e} {Er_[m].mean()/E:9.1e}")

def suite_force():
    print("== FORCE: m_F=J/D(gv), +/-lam, temoin lam=0, rayonnement ==")
    om = 0.95; eta = np.sqrt(1-om*om); E0 = 16*eta
    print(f"{'h':>5} {'v':>5} {'lam0':>9} {'J':>9} {'Dgv':>9} {'m_F/E0':>8} {'E_rad/E fin':>11}")
    for h in (1.0, 0.5):
        N = int(4000/h); dt = 0.02*h
        for v in (0.0, 0.2):
            res = {}
            for lam0 in (4e-4, -4e-4, 0.0):
                r, J = run(h=h, Ncells=N, dt=dt, T=400.0, v=v, om=om, lam0=lam0)
                t_, E_, P_, X_, qc_, Er_ = r.T
                m1 = (t_ > 15) & (t_ < 75); m2 = (t_ > 285) & (t_ < 395)
                v1 = np.polyfit(t_[m1], X_[m1], 1)[0]; v2 = np.polyfit(t_[m2], X_[m2], 1)[0]
                g1 = 1/np.sqrt(1-v1*v1); g2 = 1/np.sqrt(1-v2*v2)
                dgv = g2*v2-g1*v1
                res[lam0] = (J, dgv, Er_[m2].mean()/E_[m2].mean())
                if lam0 == 0.0:
                    print(f"{h:5.2f} {v:5.2f} {'temoin':>9} {'-':>9} {dgv:9.2e} {'-':>8} {res[lam0][2]:11.1e}")
                else:
                    print(f"{h:5.2f} {v:5.2f} {lam0:9.1e} {J:9.5f} {dgv:9.5f} {J/dgv/E0:8.4f} {res[lam0][2]:11.1e}")
            Jp, dp, _ = res[4e-4]; Jm, dm, _ = res[-4e-4]
            print(f"{'':5} {'':5} {'antisym':>9} {'':9} {'':9} {(Jp-Jm)/(dp-dm)/E0:8.4f}")

def suite_force25():
    print("== FORCE h=0.25, v=0 (3e point de convergence m_F) ==")
    om = 0.95; E0 = 16*np.sqrt(1-om*om); h = 0.25; N = int(4000/h); dt = 0.02*h
    res = {}
    for lam0 in (4e-4, -4e-4):
        r, J = run(h=h, Ncells=N, dt=dt, T=400.0, v=0.0, om=om, lam0=lam0)
        t_, E_, P_, X_, qc_, Er_ = r.T
        m1 = (t_ > 15) & (t_ < 75); m2 = (t_ > 285) & (t_ < 395)
        v1 = np.polyfit(t_[m1], X_[m1], 1)[0]; v2 = np.polyfit(t_[m2], X_[m2], 1)[0]
        g1 = 1/np.sqrt(1-v1*v1); g2 = 1/np.sqrt(1-v2*v2)
        res[lam0] = (J, g2*v2-g1*v1)
        print(f"lam0={lam0:+.1e}: J={J:.5f} Dgv={res[lam0][1]:.5f} m_F/E0={J/res[lam0][1]/E0:.4f}")
    Jp, dp = res[4e-4]; Jm, dm = res[-4e-4]
    print(f"antisym: m_F/E0 = {(Jp-Jm)/(dp-dm)/E0:.4f}")

def suite_open():
    """Controle ouvert (annexe section 5): chaine lineaire GAPPEE (Om0=1), soliton pompe
    (source A*exp(-(x-X)^2/2s^2)*cos(ws t), ws=0.5 sous le gap, phase imposee par l'horloge
    du labo), X DYNAMIQUE: Mb X'' = F_champ + F_ext. Montre m_dress = Mtot - Mb ~ 0.15 U0."""
    N = 9000; Om0 = 1.0; ws = 0.5; sig = 3.0; A = 0.02; dt = 0.05
    x = np.arange(N, dtype=float); nsp = 600
    gam = np.zeros(N)
    gam[:nsp] = 0.5*(1-np.arange(nsp)/nsp)**2; gam[-nsp:] = gam[:nsp][::-1]
    core = slice(nsp, N-nsp)
    def lp(q):
        l = np.empty_like(q); l[1:-1] = q[2:]-2*q[1:-1]+q[:-2]
        l[0] = q[1]-2*q[0]; l[-1] = q[-2]-2*q[-1]
        return l
    Tramp = 100.0; t1, t2 = 150.0, 250.0; Tend = 650.0
    for Mb in (0.001, 0.002):
        J = 0.3*(Mb+0.16e-3)
        q = np.zeros(N); p = np.zeros(N); X = 1500.0; V = 0.0
        t = 0.0; out = []
        def Fx(tt):
            if t1 < tt < t2:
                u = (tt-t1)/(t2-t1); return J/(t2-t1)*6*u*(1-u)
            return 0.0
        def rhs(q_, p_, X_, V_, tt):
            d = x-X_
            env = np.exp(-d*d/(2*sig*sig))
            rmp = 0.5*(1+np.tanh((tt-Tramp/2)/(Tramp/6)))
            s = A*env*np.cos(ws*tt)*rmp
            dS = A*env*(d/(sig*sig))*np.cos(ws*tt)*rmp
            F = np.sum(q_*dS)
            return p_, s-Om0*Om0*q_+lp(q_)-gam*p_, V_, (F+Fx(tt))/Mb
        for i in range(int(Tend/dt)):
            # RK4 explicite lisible:
            a1q, a1p, a1X, a1V = rhs(q, p, X, V, t)
            a2q, a2p, a2X, a2V = rhs(q+dt/2*a1q, p+dt/2*a1p, X+dt/2*a1X, V+dt/2*a1V, t+dt/2)
            a3q, a3p, a3X, a3V = rhs(q+dt/2*a2q, p+dt/2*a2p, X+dt/2*a2X, V+dt/2*a2V, t+dt/2)
            a4q, a4p, a4X, a4V = rhs(q+dt*a3q, p+dt*a3p, X+dt*a3X, V+dt*a3V, t+dt)
            q += dt/6*(a1q+2*a2q+2*a3q+a4q); p += dt/6*(a1p+2*a2p+2*a3p+a4p)
            X += dt/6*(a1X+2*a2X+2*a3X+a4X); V += dt/6*(a1V+2*a2V+2*a3V+a4V)
            t += dt
            if i % 4 == 0:
                dq = q[core]; dp_ = p[core]
                U = 0.5*np.sum(dp_*dp_)+0.5*Om0*Om0*np.sum(dq*dq)+0.5*np.sum(np.diff(dq)**2)
                out.append((t, V, U))
        o = np.array(out); t_, V_, U_ = o.T
        U0 = U_[(t_ > 120) & (t_ < 148)].mean()
        Vf = V_[(t_ > 480) & (t_ < 640)].mean()
        md = J/Vf-Mb
        print(f"Mb={Mb:.4f}: U0={U0:.5f} Vf={Vf:.4f} m_dress={md:.6f} m_dress/U0={md/U0:.3f}")

if __name__ == '__main__':
    arg = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if arg in ('rest', 'all'): suite_rest()
    if arg in ('boost', 'all'): suite_boost()
    if arg in ('force', 'all'): suite_force()
    if arg in ('force25',): suite_force25()
    if arg in ('open', 'all'): suite_open()
