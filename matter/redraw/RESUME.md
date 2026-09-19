# Travaux « matière » — résumé (au 19 septembre 2026, à jour de R58)

Tout ce qui suit est calculé par des scripts (`matter/scripts/`, `matter/redraw/`), chacun avec
des vérifications PASS/FAIL ; le détail étape par étape est dans `BASE.md` (R1 à R56).
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

## 3. Ce qui est dérivé (entrées : α, m_e, l'exposant 2π, le compte 9, la section du ruban)

| résultat | valeur | mesuré / réseau | écart |
|---|---|---|---|
| spectre de charge des fermions | {0, ±1/3, ±2/3, ±1} | idem | exact |
| moment magnétique de l'électron | μ_B, g = 2 | 2,0023 | 0,1 % |
| masse du muon, du tau (échelle 3, 7, 11) | 104,8 ; 1794 MeV | 105,7 ; 1777 | −0,8 % ; +1,0 % |
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

1. La règle qui produit Koide : le 45° a un candidat (un quantum de circulation dans le canal
   commun, un dans le canal des différences, R62) ; la phase du doublet 2/9 n'a aucun mécanisme
   (les identités trouvées à R63 sont des rationnels par construction). Et pourquoi 3, 7, 11
   (les DQD neutres s'ajoutent par paires, jamais un seul).
2. Une dynamique : aucun processus n'est calculé. La durée de vie comme « respiration de Z » (R64)
   exige une désadaptation en m² sur une échelle fixe de ~1 TeV que la base n'a pas ; une
   résonance à Z ≠ Z₀ est exclue (R65 : elle meurt en moins de cent tours).
3. La largeur du ruban, non dérivée depuis l'anneau fermé (R53), sauf par les trois jonctions (R54).
4. Sur quel cercle tourne la circulation collective des trois circuits du nucléon (0,67 fm).
5. Ce qui épingle l'anneau du pion à r_e/2.
6. Le spin : un seul compte (circuits ou anneau), et le neutron trop lourd de quelques %.
7. Le −1 sous un tour complet (statistique d'échange), non abordé.

## 7. Fichiers

- `matter/redraw/BASE.md` : les énoncés R1–R56 et leurs conséquences calculées, avec deux bilans.
- `matter/redraw/*.py` : 60 scripts (tous PASS), un par étape ; `matter/scripts/` : la phase
  précédente (holonomie, Faddeev, impédance du proton, topologie).
- `matter/OPEN_PROBLEMS.md`, `matter/CONSTRAINTS.md` : l'état des problèmes ouverts et des
  contraintes.
