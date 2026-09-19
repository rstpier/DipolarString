# Travaux « matière » — résumé (au 19 septembre 2026, à jour de R98)

Tout ce qui suit est calculé par des scripts (`matter/scripts/`, `matter/redraw/`), chacun avec
des vérifications PASS/FAIL ; le détail étape par étape est dans `BASE.md` (R1 à R98).
Verdicts : **dérivé** (sort des règles de la base sans nombre réglé), **conditionnel** (dépend
d'une lecture non tranchée), **coïncidence** (nombre juste, sans mécanisme), **exclu**.

## 1. Point de départ et méthode

Le manuscrit V2.9/2.10 a été abandonné comme référence (« il est faux ») ; on a redessiné la
matière à partir des énoncés de l'auteur (R1 à R16, puis les suivants), en calculant chaque fois ce
que l'énoncé impose et ce qu'il casse. Contraintes utilisées : neutralité de l'atome, absence de
rayonnement de l'électron, absence d'électron excité, LEP, masses et moments mesurés, réseau QCD
pour la tension de corde et l'écart proton–neutron.

## 2. L'image finale, en clair

- **Le milieu.** L'espace est fait de cordes dipolaires neutres (DQD), lignes bifilaires accordées
  à l'impédance du vide Z₀ ; l'accord fixe la forme de la section (carrée), pas sa taille.
- **Les brins.** Toute particule est faite de brins portant −1/3, 0 ou +1/3, trois par fermion :
  e = (−,−,−), ν = (0,0,0), u = (+,+,0), d = (−,0,0). Trois brins donnent exactement les charges
  observées et rien d'autre ; la couleur est la place du brin impair (3 arrangements pour u et d,
  un seul pour e et ν). Un brin neutre est un DQD complet.
- **Le fluide.** Ce n'est pas une onde mais un courant qui tourne (vortex) : un quantum de
  circulation par circuit fermé, sans harmonique, sans rayonnement.
- **L'électron.** Un anneau fermé où le fluide tourne dans un sens (moment magnétique μ_B, spin ½,
  pas de rayonnement, pas d'état excité) ; la moitié de sa masse est dans la circulation
  (πℏc/trajet), l'autre moitié immobile dans les trois jonctions de l'anneau. Cette moitié
  immobile retrouve l'écart propre du DQD (le « 37,1 » du manuscrit) à 3 %.
- **Les générations.** μ = e + 4 brins neutres (7), τ = e + 8 (11) ; les masses suivent
  m = m_e(n/3)^{2π} à 1 % (exposant épinglé par les données à 0,15 %, unique parmi les constantes
  simples), mais la structure exacte est celle de Koide : le vecteur des racines de masse fait
  45° avec la diagonale (moitié commune, moitié dans les différences) et tourne de 2/9 rad,
  exact à 10⁻⁵ ; la loi de puissance donne déjà 0,220 rad (1 %). Le 2/9 est aussi le produit
  des charges u et d, sans mécanisme qui en fasse un angle.
- **Les quarks et les nucléons.** Un quark est le même anneau à trois brins, à l'échelle
  ℓ₁(9) = 0,81 fm : 254 MeV chacun ; un nucléon en a trois (763 MeV immobiles, la seule partition
  des 9 brins qui marche) plus une circulation d'anneau qui porte le spin et le moment. Le neutron
  est un proton avec une boucle négative de type électron accrochée par un joint.
- **Les mésons.** Des anneaux fermés de 4 à 6 brins à r_e/2 : 140 MeV pour le pion.
- **Les forces.** La force nucléaire est la jonction pôle à pôle (2 MeV à 0,7 fm, l'échelle du
  deutéron) ; la tension forte est l'énergie d'un quantum de circulation par unité de trajet,
  πℏc/ℓ₁(9)² = 939 MeV/fm, soit √σ = 430 MeV, au centre de la valeur du réseau (420–440).

## 3. Ce qui est dérivé

Sans paramètre supplémentaire ajusté, conditionnellement aux entrées α, m_e, l'exposant 2π,
le compte nucléonique 9 et la section du ruban ; 2π et 9 sont des entrées structurelles non dérivées.

| résultat | valeur | mesuré / réseau | écart |
|---|---|---|---|
| spectre de charge des fermions | {0, ±1/3, ±2/3, ±1} | idem | exact |
| moment magnétique de l'électron | μ_B, g = 2 | 2,0023 | 0,1 % |
| masse du muon, du tau : forme C₃ (Koide) avec a² = 2\|b\|² et φ = 2/9 **posés** (R76) | 105,66 ; 1777,0 MeV | 105,66 ; 1776,9 | 10⁻⁵ ; 7·10⁻⁵ |
| (l'échelle (n/3)^{2π} n'est plus une loi de masse : approximation à 1 %, exclue par l'absence de 4e génération, R73/R76) | 104,8 ; 1794 | | −0,8 % ; +1,0 % |
| pion (anneau fermé) | 2m_ec²/α = 140,05 MeV | 139,57 | +0,3 % |
| tension forte √σ | 430 MeV | 420–440 | centre |
| écart neutron–proton, part forte | 2,47 MeV | 2,52 ± 0,29 | dans l'erreur |
| écart neutron–proton, total | 1,27 à 1,35 MeV | 1,293 | −2 à +4,5 % |
| rayon de charge du proton (charges aux bouts des bras) | 0,81 fm | 0,84 | −3 % |
| Δ − N (trois fois la circulation d'anneau) | 250 à 336 MeV | 294 | ±15 % |
| force nucléaire (tri-jonction) | 2 MeV | deutéron 2,2 | 10 % |
| écart du DQD depuis g = 2 | ƛ/38,1 | manuscrit ƛ/37,1 | 2,7 % |

## 4. Ce qui est conditionnel ou coïncidence

- La largeur du ruban (w = 6ƛ/π² depuis les trois jonctions) : l'écart p–n n'y est sensible qu'à
  un facteur 2 près (±11 %) ; la forme du pôle (bouchon) vient du mécanisme d'entassement pris
  tel quel.
- Le rayon de l'anneau du nucléon, 0,67 fm : égal à √3 fois le rayon de circuit à 0,1 %, et
  l'enveloppe de trois circuits tangents vaut r_p à 0,6 % ; sans mécanisme.
- m_p = 4ℏc/r_p à 0,05 % ; le muon à 0,79 fm et 1,09 fm pour les brins + et − ; μ(Δ⁺⁺) = 6,4 μ_N.
- Le neutron sort 1,5 à 7 % trop lourd ; le spin est compté à la fois par circuits et par anneau.
- « Dynamisme de l'espace » (R57) : des DQD du vide qui circulent donneraient à l'électron son
  quantum par héritage, au prix d'une énergie du vide 10³³ fois la cosmologique ; l'expansion du
  milieu et le confinement par le milieu sont exclus.

## 5. Ce qui a été exclu en route

Nucléons et quarks sur l'échelle des leptons ; masses par comptes de brins ; « un joint = un
nombre » ; l'enroulement comme mécanisme de force ; le tube de flux dans la base telle qu'écrite ;
une densité de charge fixe du fluide ; la charge de vortex 5,85 e comme charge qui circule ; les
pôles libres et l'onde stationnaire (pas de moment, rayonne, harmoniques) ; l'épingle fermée (30 %
du moment) ; le mode transverse à 4 MeV (artefact de la lecture onde) ; la lecture « déformation »
de l'exposant.

## 6. Ce qui reste à fournir

Passe de cohérence (R71–R75) : cinq points examinés, trois fermés, un devenu une décision, un ouvert.

- Fermés : la brique n'a pas de taille (Z et l'énergie par longueur d'une ligne bifilaire sont
  invariantes d'échelle, R71) ; le nucléon a un seul compte de spin (l'anneau) et un seul nombre
  pour le neutron, −1,8 % avec le pôle de plus basse énergie (R72) ; l'anneau stable et l'émission
  de photons ne se contredisent pas, l'orbite stationnaire venant de la même horloge que le tour
  (de Broglie, R75).
- Décision : les masses des leptons sont Koide (à 10⁻⁵), n est une étiquette, 2π l'approximation
  à 1 % ; aucune correction à un paramètre ne fait mieux (R73).

Une prédiction sans nombre libre, qui mord (R89) : si l'électron est l'anneau, sa taille doit
reproduire le terme de Darwin de l'hydrogène, r_rms = (√3/2)ƛ = 334 fm ; la charge de R84 donne
363 fm (+18 % en ⟨r²⟩, +31 GHz sur le 1S, mesuré au kHz). Il faut un calcul atomique propre à la
base ; il tranchera. La matière noire, elle, n'est pas prédite : un défaut neutre verrouillé a
E = C/R sans minimum (aucune taille pour un objet sans charge), et épinglé au vide il pèse
~1 MeV, exclu par le BBN (R90) ; le secteur sombre de la base, ce sont les neutrinos, dont elle
ne fixe pas la taille. L'intrication : la base n'a pas de loi de mesure, et toute loi locale à
variables définies donne |S| ≤ 2 contre 2,83 mesuré (R91) : exclue comme théorie locale, pas
comme modèle entier. Ce qu'elle a déjà de non local, le quantum partagé de la paire (contrainte
globale) et la sphère de Bloch, ne suffit pas : sans la règle de Born (unique réponse donnant
−cos θ) la contrainte surchoote à S = 4, et l'onde du milieu (c₀ = c) est trop lente pour être le
canal (> 10⁴ c requis) ; aucune correction nouvelle ne sort (R92). Pour un seul anneau, Born
se réduit à la dynamique de l'analyseur : si un 2-port linéaire sans perte réalise le mélange
SU(2), flux quadratique + quantum indivisible ⇒ P(+) = |A₊|² (R93, conditionnel ; Tsirelson y
est une vérification avec le produit tensoriel importé, pas une dérivation). Mais aucun
analyseur géométrique réel du doublet E ne donne le demi-angle (toutes les constructions
donnent 1/2 à θ = 0 ou le mauvais angle, R94a) : le cos(θ/2) exige la coordonnée complexe du
fibré. Le relèvement SU(2) lui-même est dérivé : le mode antipériodique donne U(2π) = −1 sur
l'axe, et la loi de groupe des rotations force −1 pour tout axe, avec la connexion de courbure
Ω/2 et c₁ = 1 retrouvé sans Wigner (R95). Deux composantes circulaires candidates de charge
axiale ±½ sont compatibles avec le quantum R87 (un quantum de la rotation de section à ω_circ/2
vaut exactement E_circ) et avec un spin ½ par la carte de Hopf, mais ce n'est qu'un pont
mathématique : le fluide n'a pas encore ce C² (R96, corrigé). Les séparateurs de modes naturels
construits par projection du champ, avec ou sans la phase de circulation, et même après
orthonormalisation des canaux, ne portent pas le demi-angle (écart 0,25 à 0,44 à cos²(θ/2),
jusqu'à 75 % de la norme hors des deux canaux) : le e^{±iθ/2} est une propriété de la
transformation de l'orientation, pas des recouvrements de champs (R97 ; pas un théorème
d'inexistence). La réorientation par le seul couplage de la base, μ_B dans un champ : le
couplage conserve S·â exactement (précession, g = 2), la déflexion classique continue et
l'alignement amorti sont exclus par l'expérience, et deux sorties ±â avec conservation en
moyenne de S·â donnent p = cos²(θ/2), Born, sans |ψ|² (R98). Non dérivé : le saut discret
lui-même. Puis l'état joint (R94b) ; et le verrou N = 1.

Ce qui reste vraiment :

1. **Le −1 sous un tour complet** (Pauli) : la base a l'attache au milieu (tour de ceinture) mais
   pas de mode antipériodique sur un anneau de trois branches simples (R74). Candidat : le triple
   (Z₃) n'a pas d'élément d'ordre 2, la section carrée (Z₄) en a un ; une demi-torsion de section
   par circuit rendrait un mode impair antipériodique, à coût 0,5 % (R78) ; ce mode est construit :
   le motif dipolaire (cœur de la circulation hors du centre de la section), irrep E de Z₄, −1
   exactement sous une demi-torsion, +i sous un quart (R82). L'anneau ne peut pas être une toupie
   rigide (un j = 3/2 suivrait à 3mc² = 1,5 MeV, exclu) ; sans inertie d'orientation, l'orbite de
   rotation avec le mode dipolaire est la sphère de Bloch de j = ½, et l'attache par les jonctions
   donne l'échange par le tour de ceinture (R83). La courbure déplace le cœur de la circulation de
   34 fm vers l'intérieur (22 % de la section, flux exclu) : l'électron porte le motif dipolaire
   sans choix (R84, sous une condition constitutive posée : courant de surface, flux interne
   contraint). Le nombre de Chern du mode sur la sphère des orientations vaut 2 pour le cœur
   déplacé seul (hélicité 1) et 1 seulement avec le demi-tour de section : c₁ mesure t, il ne le
   dérive pas (R85). Le quantum de circulation de la base, πℏc/trajet, est exactement une action
   de demi-tour, ∮p dl = πℏ (R86, corrigé : conditionnel, la lecture onde et l'appel à g = 2
   étaient circulaires). Ce quantum se dérive sans g comme la moitié du quantum entier d'une
   boucle neutre périodique brisée symétriquement en deux (R87) ; alors E_circ = m/2, S = ℏ/2 et
   g = 2 sont des sorties. La mère est transitoire parce qu'elle n'a aucun verrou (charge 0,
   Lk = 0) et qu'elle est un photon fermé (E = hc/L_m, p = h/L_m) ; les filles sont verrouillées par
   leur charge et par leur demi-quantum, sans forme libre ; l'annihilation redonne deux quanta
   entiers de λ = L_d exactement (R88). Restent posés : les trois entrées de R87 (n = 1 sur la
   mère, p′ conservé à la brisure, partage C-symétrique), que la phase physique soit e^{iS/ℏ},
   la géométrie de la fermeture du photon sur trois DQD, et l'absence d'inertie d'orientation.
2. **Le mécanisme de Koide** : les masses sont les carrés des trois valeurs propres d'un opérateur
   de racine de masse C₃ (circulant), et Koide est la règle a² = 2|b|² (le 45°) avec la phase
   φ = 2/9, toutes deux posées (R76). Le 45° équivaut exactement à « normes isotrope et
   déviatorique égales » = un quantum par canal irréductible (pas par degré de liberté, qui
   donnerait deux leptons sans masse) ; aucun extremum ne le sélectionne, il faut une
   conservation ; la statique des trois brins donne la forme mais un doublet dégénéré, la phase
   est une holonomie indépendante des normes (R79). Ce quantum par canal n'est pas celui de R47 :
   dans l'anneau dessiné (trois brins en série) le doublet n'a même pas de circuit, et dans la
   seule topologie où il en a un (trois boucles coaxiales, compatible avec R53/R54) ni le couplage
   statique (M/L = 0,3–0,6 contre 0,71) ni un quantum égal par circuit ne donnent le 45° (R80).
   Mais le 45° est le « moitié-moitié » de g = 2 lu dans la base des positions : si l'opérateur
   de racine de masse est la matrice d'amplitudes du réseau (jonction = nœud = poids sur site,
   circulation = lien = poids de saut, masse = poids total), alors g = 2 ⇒ a² = 2|b|² sans √2 posé,
   au niveau de la famille (R81). Le pas de génération
   est deux DQD (le spin interdit un seul, R70) et le quatrième pas est le premier (U³ = I).
   La liaison des paires est un verrou topologique (Lk = Tw + Wr conservé, forme gelée, R77),
   mais ce verrou tient aussi une troisième paire (Wr = 1 à r₀ = 0,376) : l'absence de quatrième
   génération doit venir du spectre (trois valeurs propres), pas de la liaison.
3. **Le mur des durées de vie** : M = (96π²)^{1/4}/√G_F = 1,625 TeV par identité (R68, corrigé
   R76) : la loi de fuite a la structure dimensionnelle du faible, mais G_F reste une entrée ;
   avec elle, τ(n) = T(n)(M/m)⁴/C(n) donne le muon à 0,4 %, le tau à 13 % (0,5 % avec la
   correction forte), l'électron stable (R69). La base ne fabrique ni G_F, ni M_W, ni l'angle.
4. **La largeur du ruban**, non dérivée depuis l'anneau fermé (R53), sauf par les trois jonctions
   (R54) ; et le partage 2/3, 2/3, −1/3 de la circulation du nucléon (R72), posé.
5. **L'amplitude d'un saut** (émission, désintégration) : la base a les états et les horloges, pas
   les amplitudes.
6. Trois coïncidences à 1 % sans mécanisme : l'anneau du nucléon à 0,67 fm (√3 × R_c), le pion à
   r_e/2, m_p = 4ℏc/r_p.

## 7. Fichiers

- `matter/redraw/BASE.md` : les énoncés R1–R98 et leurs conséquences calculées, avec trois bilans
  (R20, R31, R49) et la passe de cohérence R71–R75.
- `matter/redraw/*.py` : 93 scripts (tous PASS), un par étape ; `matter/scripts/` : la phase
  précédente (holonomie, Faddeev, impédance du proton, topologie).
- `matter/OPEN_PROBLEMS.md`, `matter/CONSTRAINTS.md` : l'état des problèmes ouverts et des
  contraintes.
