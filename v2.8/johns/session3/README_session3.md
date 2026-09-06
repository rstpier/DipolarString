# V2.8 — pont Johns, session 3 (2026-07-23)

## Dependance
Tous les scripts chargent `S_johns.npy` (matrice de diffusion du noeud, reconstruite et
validee en session 1) DANS LE MEME REPERTOIRE. Il est inclus dans cette archive.
Reproduction session 1: `johns_scn.py` puis `select.py` (archive v28_johns_session1.zip).

## Domaine de validite de la branche (branch.npy: colonnes kz, omega, poids-coeur)
- Branche liee suivie par recouvrement, Nx=15, defaut D3 minimal (tresse + R(2pi/3), phi=0).
- DOMAINE INTERPRETE: 0.9 <= kz <= 2.4 (poids-coeur 0.42-0.50, suivi propre).
- v_g = domega/dkz < 0 (branche retrograde); maximum stable |v_g|/c0 = 0.1557 a kz = 1.5.
- ZONE NON INTERPRETEE: 2.5 <= kz <= 2.8 (anticroisement probable, poids chute a 0.12;
  toute derivee dans cette zone, y compris le |v_g|/c0 ~ 0.30 apparent, est NON exploitable).
- kz >= 2.9: le suivi retrouve un etat localise (poids ~0.35) mais l'identite de branche
  n'est pas garantie a travers la zone melangee.

## Statut de la vitesse
"Seconde vitesse effective du defaut": vitesse de groupe DISPERSIVE, pas une constante;
0.1557·c0 est un maximum de branche. L'observable representant c_int dans la relation
transverse (max|v_g|, |v_g(k*)|, limite de seuil, ou coefficient de courbure de bande)
n'est PAS encore choisie — elle devra etre imposee par la physique des contraintes R/S/E,
sous peine de constituer un nouveau bouton discret.

## Conversion CD (premiere passe, arithmetique)
hbar/Delta_t = 0.487 MeV (c_lien = 2c0 derive du noeud; l1 = 810.4 fm).
hbar*omega*(kz=2.0) = 0.33246 x 0.487 = 0.1619 MeV; omega_e/omega* = 3.156.
Proximite avec pi (0.46%) NOTEE, NON utilisee. l_ev(ancre) = l1/kappa = 398 fm vs
lambda_C_barre = 386 fm (3.1%): coincidence a l'ancre, kappa depend de kz.
