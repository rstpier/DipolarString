"""
2PN comparison: isotropic exponential metric (DS) vs Schwarzschild in isotropic coordinates.
U = GM/(r c^2).  DS:  g00 = -exp(-2U),  gij = exp(+2U) d_ij
Schwarzschild isotropic: g00 = -((1-U/2)/(1+U/2))^2,  gij = (1+U/2)^4 d_ij
"""
import sympy as sp

U = sp.symbols('U', positive=True)

g00_DS  = -sp.exp(-2*U)
gij_DS  =  sp.exp(2*U)
g00_Sch = -((1-U/2)/(1+U/2))**2
gij_Sch =  (1+U/2)**4

s00_DS  = sp.series(g00_DS,  U, 0, 4).removeO()
s00_Sch = sp.series(g00_Sch, U, 0, 4).removeO()
sij_DS  = sp.series(gij_DS,  U, 0, 4).removeO()
sij_Sch = sp.series(gij_Sch, U, 0, 4).removeO()

print("g00 DS          :", sp.expand(s00_DS))
print("g00 Schwarzschild:", sp.expand(s00_Sch))
print("Delta g00        :", sp.simplify(sp.expand(s00_DS - s00_Sch)))
print()
print("gij DS          :", sp.expand(sij_DS))
print("gij Schwarzschild:", sp.expand(sij_Sch))
print("Delta gij        :", sp.simplify(sp.expand(sij_DS - sij_Sch)))
print()

# PPN reading at 1PN: g00 = -(1 - 2U + 2 beta U^2), gij = (1 + 2 gamma U)
beta_DS  = sp.expand(s00_DS).coeff(U, 2)/2 * (-1) * (-1)  # coefficient of U^2 in -g00 /2... do it cleanly:
c2 = sp.expand(-s00_DS).coeff(U, 2); print("coeff U^2 in -g00 (DS)   = 2*beta =", c2, "=> beta =", c2/2)
c2s = sp.expand(-s00_Sch).coeff(U, 2); print("coeff U^2 in -g00 (Schw) = 2*beta =", c2s, "=> beta =", c2s/2)
c1 = sp.expand(sij_DS).coeff(U, 1); print("coeff U   in gij (DS)    = 2*gamma =", c1, "=> gamma =", c1/2)

# First deviations:
d00 = sp.expand(s00_DS - s00_Sch)   # expect O(U^3)
dij = sp.expand(sij_DS - sij_Sch)   # expect O(U^2)
print("\nLeading deviation in g00:", d00)
print("Leading deviation in gij:", dij)

# Order of magnitude for the double pulsar PSR J0737-3039:
# U at periastron ~ G(m1+m2)/(a(1-e) c^2). m_tot ~ 2.587 Msun, a ~ 8.8e8 m, e ~ 0.088
G, c, Msun = 6.674e-11, 2.998e8, 1.989e30
m_tot = 2.587*Msun; a = 8.784e8; e = 0.0878
U_peri = G*m_tot/(a*(1-e)*c**2)
print(f"\nDouble pulsar: U at periastron ~ {U_peri:.2e}")
print(f"Relative 2PN deviation of periastron advance ~ O(U) ~ {U_peri:.1e}")
print(f"Measured periastron precision (J0737-3039) ~ 1e-6 ; GR 2PN term itself ~ 1e-5 of 1PN")
