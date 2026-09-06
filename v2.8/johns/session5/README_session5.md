# Session 5 (23/07) — les trois nombres spectraux de la corde
Methode: energie propre du defaut = decalage spectral demi-quantum (Casimir de connectivite)
  mu_s = (hbar/2) * [Sum|om|_D3 - Sum|om|_D0] integre sur kz  (casimir.py)
Declaration: seule hypothese ajoutee = quantification (1/2)hbar*omega du reseau classique.
Convergence en taille: 3.4e-8 entre Nx=9 et Nx=11 (defaut local, shift local).

Resultats (hbar/Dt = 0.487 MeV, l1 = 810.4 fm):
  mu_s        = 12.33 keV/cellule
  E_boucle    = 3 cellules x mu_s = 37.0 keV      [boucle = 3 cellules droites: la
                courbure/les coins NE SONT PAS comptes — approximation declaree]
  rho_bulk    = (1/2)<Sum|om|>_BZ = 3*pi exactement (9.4248) -> 4.590 MeV/cellule
  l_s         = mu_s/rho_bulk = 0.00269 cellules = 2.18 fm

Confrontations (grille pre-enregistree: ±20% piste serieuse; x2-3 microstructure
incomplete; decades rejet):
  R : ECHEC STRUCTUREL (session 4, omega_e dans le continuum) — en attente de la
      reverification independante de Raff.
  S : l_s = 2.18 fm vs cible 4.5 fm -> facteur 2.07, Delta_log = 0.315
  E : 37.0 keV vs cible 10.7 keV   -> facteur 3.46, Delta_log = 0.539
Verdict grille: bande mediane «microstructure incomplete» pour S; E legerement au-dela.
Le critere strict Delta_log<0.08 N'EST PAS atteint.

RESULTAT STRUCTUREL NEUF: S et E sont MUTUELLEMENT INCOMPATIBLES pour toute tension
unique — S exige mu = 25.5 keV/cellule, E exige mu = 3.57 keV/cellule (rapport 7.2,
Delta_log 0.86). Les deux cibles continues ne peuvent pas etre satisfaites simultanement
par une ligne de defaut, quelle que soit sa valeur — sauf si leurs observables se
traduisent differemment au niveau reseau.

Veille numerologique (notes, non reclamees): m_e/E_boucle = 13.81 ~ kR=13.8 de P3;
rho_bulk = 3*pi exact (analytique probable, a deriver).
