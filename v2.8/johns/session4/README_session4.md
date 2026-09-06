# Session 4 (23/07) — recensement spectral et bandes exactes
- census.py: spectre COMPLET (0,pi) de D3 a kz=2.0 et 0.6, supercellule 11x11, etats
  localises (poids>0.15) + tentative de nombre d'enroulement (OUTIL TROP GROSSIER:
  anneau r~3, aliasing — les valeurs d'enroulement NE SONT PAS interpretables en l'etat).
- bulkgaps.py: bandes du bulk EXACTES (Bloch 12x12, grille 61x61 en k_perp) — les "gaps"
  de la supercellule 11x11 etaient pour la plupart des artefacts de taille finie.
Resultats graves:
1. SEUL vrai gap au-dessus du cone: le stopband superieur (ex.: (2.5708, pi) a kz=2.0).
2. Ancre 0.33246 a kz=2: confirmee HORS bande exacte (densite nulle) — lie authentique.
3. omega_e = 1.0493 (0.511 MeV) est DANS le continuum du bulk aux deux kz testes
   -> aucun etat lie du defaut rectiligne minimal ne peut exister a omega_e;
   la contrainte R (resonance transverse a omega_e) est INSATISFAISABLE pour ce defaut.
4. Famille liee superieure AUTHENTIQUE dans le stopband a kz=2.0:
   omega = 2.585, 2.658, 2.674, 2.804, 2.925, 2.995 (poids 0.30-0.52)
   soit hbar*omega = 1.259 a 1.459 MeV, sous E1 = pi x 0.487 = 1.530 MeV.
5. La "famille mediane" (0.68-0.81 MeV) tombe DANS la bande exacte: resonances ou
   artefacts, PAS des etats lies — retiree de l'inventaire.
6. kz = 2pi/3 = 2.094 (quantification n=1 d'une boucle fermee de 3 cellules,
   circonference 2*pi*lambda_C_barre = 3.0 cellules): voir sortie ci-dessus.
Non fait (reporte): identification propre des modes (caractere T3 / harmoniques),
anticroisement 2.5-2.8, designation de l'observable c_int depuis R/S/E.
