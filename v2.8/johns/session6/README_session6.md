# Session 6 (23/07) — la boucle fermee par compactification exacte
Methode (loop6.py): la tresse a periode 3 en z, donc une boite z-periodique de Lz=3 AVEC
le defaut EST la boucle fermee (circonference 3*l1 = 2*pi*lambda_C_barre exactement).
Decomposition de Bloch exacte: kz permis = {0, +-2pi/3}, d'ou
  E_boucle = 0.487*[ds(0)+ds(+2pi/3)+ds(-2pi/3)]  — AUCUNE matrice 3D requise.
Geometrie torique declaree: E_courbure et E_coins ABSENTS par construction; on isole
E_fermeture seul. Unitarite 1e-14 dans tous les secteurs; convergence Nx=7 vs 9: 3e-4 relatif.

## Resultats
Secteur phi=0 (defaut minimal):
  ds(0)=+0.041977, ds(+-2pi/3)=+0.020070 (degeneres)
  E_boucle = 39.99 keV ; 3*mu_s = 38.31 keV ; E_FERMETURE = +1.68 keV (+4.4%)
=> la fermeture est une petite correction. AUCUN facteur ~13 n'emerge de la fermeture.
Verdict sur le rapprochement m_e/E_boucle ~ kR: sous la regle pre-enregistree («emerger
ou mourir»), il MEURT au niveau fermeture — seul refuge restant: courbure/coins d'un
anneau reellement plie (non calcule ici), mais une correction de type courbure est
attendue a l'echelle du keV, pas d'un facteur 13.

Secteurs de flux phi=+-2pi/3 (permis pour la boucle FERMEE, e^{3i phi}=1, discrets):
  E = -126.06 keV, DOUBLEMENT DEGENERES (+-), sans parametre continu.
DECOUVERTE ET PROBLEME: le decalage spectral de point zero est NEGATIF — la boucle a flux
abaisse le point zero du vide de 126 keV. Trois lectures possibles, NON tranchees:
 (1) instabilite (le vide nucleerait spontanement des boucles a flux) — suspect;
 (2) l'interpretation demi-quantum (1/2)hbar*omega cesse d'etre valable pour ces secteurs
     (la "tension spectrale candidate" n'est pas une energie totale);
 (3) les secteurs a flux sont les etats habilles physiques (modes lies profonds tirent le
     spectre vers le bas) et il manque des contributions positives non comptees.
AUCUNE interpretation n'est adoptee. Veille numerologique: m_e/|E_flux| = 4.05 (note,
non reclame).

## Bilan pour l'objectif V2.8
Le spectre fournit: mu_s=12.77 keV/cellule, E_fermeture=+1.68 keV, secteurs discrets
{+40.0, -126.1, -126.1} keV pour la boucle fermee. Ce qui manque toujours: un invariant
spectral de l'ordre de 511 keV, et la resolution du signe negatif des secteurs de flux.
