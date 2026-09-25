# Le secteur hadronique de la base — consolidation (R109, 25 septembre 2026)

Tout ce qui suit est recalculé d'un bloc par `hadron_sector.py` (6/6), à partir des entrées de la
base, avec les comparaisons à l'expérience et au réseau et les sensibilités. Les étapes d'origine
sont dans `BASE.md` (R16–R49, R72) ; ce chapitre les réordonne en une seule chaîne.

## 1. Ce que ce secteur est, et ce qu'il n'est pas

- **Des rapports, ancrés sur ƛ_e.** Aucune masse absolue ; tout est une longueur en unités de
  ƛ_e = ℏc/m_e ou une énergie en unités de m_e c², comme la base l'exige (`CONSTRAINTS.md` §1).
- **Indépendant de l'électron.** Rien n'utilise la position de la charge de l'électron, son
  moment, son spin ni son rayon, tous exclus ou suspendus en R100–R108. La seule grandeur de
  l'électron qui entre est m_e, par ƛ_e. L'ancre du secteur est une longueur,
  ℓ₁(9) = (2π/3)·ƛ_e·3^{−2π} = 0,8128 fm ; « le quark est l'électron à l'échelle 9 » se lit
  E_quark = (m_e/2)·3^{2π}, et sous la lecture champ de R104 le facteur 2π/3 est posé, pas dérivé
  de la géométrie de l'électron.
- **Un secteur d'objets étendus.** Le proton mesure 0,84 fm ; les brins chargés à 0,8 fm y sont
  légitimes, là où ils ne l'étaient pas pour l'électron (G_E = 1 à 10⁻³ fm, R100).

## 2. Les entrées

| entrée | valeur | statut | où |
|---|---|---|---|
| ℏc, α, m_e | CODATA | données | — |
| l'exposant p de l'échelle ℓ₁(n) = ℓ₁(3)(3/n)^p | 2π (épinglé par μ et τ à 0,15 % ; fenêtre 6,2758–6,2926) | reparamétré | R14, R55, R73 |
| le compte de brins du nucléon | 9 (p = uud, n = udd, trois brins par quark) | structurel, non dérivé | R28 |
| la section du ruban | carrée, w_e = 4ƛ_e/π², w/ℓ₁ = 6/π³ à toute échelle | conditionnel (g = 2 et trois jonctions, R46, représentation aujourd'hui fermée) | R17, R46, R71 |

## 3. Les règles

1. **Un quantum de circulation par circuit fermé**, d'énergie πℏc/trajet (R40, R47 ; dérivation
   conditionnelle depuis une boucle mère brisée, R87). Posé.
2. **Trois brins par fermion**, charges {−1/3, 0, +1/3} ; u = (+,+,0), d = (−,0,0) ; la couleur
   est la position du brin impair ; p et n ont neuf brins (R28). Dérivé du spectre de charge, brins
   posés.
3. **Un circuit par quark** : des 30 partitions des 9 brins en circuits, seule {3, 3, 3} donne
   Σ1/k = 1, la part statique que μ_p exige (R41). Dérivé sous la règle 1.
4. **Un seul anneau porte le spin et le moment** ; les circuits de quark sont statiques (R72).
5. **La charge d'un quark est sur ses brins chargés** (deux pour u, un pour d) ; la section du
   ruban est la coupure du terme propre (plaque : a = w/4) ; les bras sont à 120° avec les charges
   aux bouts (R43–R45).
6. **Le pôle est la forme de plus basse énergie** dans la section, le bouchon cubique (R5, R72).
7. **La circulation se partage comme la charge** : 2/3, 2/3, −1/3 (R72). Posé (c'est le partage
   de SU(6)).

## 4. La chaîne

1. ℓ₁(9) = (2π/3)ƛ_e·3^{−2π} = 0,8128 fm ; w₉ = w_e·3^{−2π} = 0,1573 fm.
2. Circuit de quark : E_q = πℏc/(3ℓ₁(9)) = 254,2 MeV. Trois circuits : 762,7 MeV = 81,3 % de m_p.
3. Anneau : S = R·E/c = ℏ/2. Trois lectures : masse + spin, E = m_p − 763 = 176 MeV, R = 0,562 fm ;
   μ_p avec e sur l'anneau, R = μ_p ƛ_p = 0,587 fm, E = 168 MeV ; Δ − N = 2E (trois circulations
   alignées, 3ℏ/2), E = 147 MeV, R = 0,672 fm, d'où N = 910 MeV (−3,1 %) et Δ = 1203 MeV (−2,3 %)
   avec un seul rayon. Le rayon n'est pas fixé par la base.
4. Écart n − p : le brin neutre de plus du neutron est un DQD complet, d'énergie bifilaire
   4πK/(9ℓ₁(9)) = (2α/3)·3^{2π}·m_e c² = 2,474 MeV, la partie forte ; le Coulomb du proton, mutuel
   K/(√3ℓ₁) = 1,023 MeV (nul pour p, −0,34 pour n) et propre sur les brins chargés, donne
   1,324 MeV (plaque) ou 1,269 MeV (bouchon).
5. Tension forte : πℏc/ℓ₁(9)² = 938 MeV/fm, √σ = (3/(2√π))·3^{2π}·m_e c² = 430,3 MeV.
6. Force nucléaire : une paire de pôles unitaires, K/s = 2 MeV à s = 0,720 fm.
7. Le pion n'est pas fait de circuits (R42 D) : c'est l'anneau fermé à r_e/2, 2m_e/α = 140,05 MeV.

## 5. Les résultats

| grandeur | base | mesuré / réseau | écart | statut | R |
|---|---|---|---|---|---|
| spectre de charge, couleur, 9 brins | {0, ±1/3, ±2/3, ±1} | idem | exact | dérivé (brins posés) | R28 |
| circuit de quark | 254 MeV | quark constituant ~313 | −19 % | conditionnel (2π, 9) | R41 |
| part statique du nucléon | 763 MeV | 770–813 (μ_p) | −1 à −6 % | dérivé sous la règle 1 | R41 |
| nucléon, total | 910–938 MeV | 938,3 | −3 à 0 % | conditionnel (rayon posé) | R72, R48 |
| Δ − N | 250–351 MeV ; 294 avec R = 0,672 | 294 | ±15 % | conditionnel (rayon) | R48 |
| m_n − m_p, partie forte | 2,474 MeV | 2,52 ± 0,29 (BMW) | 1σ | dérivé | R42 |
| m_n − m_p, total | 1,27 (bouchon), 1,32 (plaque) MeV | 1,293 | −1,9 %, +2,4 % | conditionnel (section, forme du pôle) | R44, R45, R72 |
| part QED implicite | −1,20 MeV | −1,00 ± 0,16 | 1,3σ | conséquence | R44 |
| μ_p/μ_n | −3/2 | −1,460 | +2,7 % | posé (partage) | R72 |
| √σ | 430,3 MeV | 420–440 | centre | conditionnel (un quantum par bras) | R39, R40 |
| force nucléaire | 2 MeV à 0,72 fm | deutéron 2,22 | échelle | calibré | R20, R22 |
| pion | 140,05 MeV | 139,57 | +0,34 % | identification (2m_e/α), coïncidence connue | R16 |
| r_p | ℓ₁(9) = 0,813 ; enveloppe 0,836 fm | 0,841 | −3,3 % ; −0,6 % | coïncidences | R43, R49 |
| m_p = 4ℏc/r_p | 3,998 | 4 | 0,05 % | coïncidence | R20 |
| μ(Δ⁺⁺) | 6,4 μ_N | 3,7–7,5 | dedans | enregistré | R48 |

## 6. Les sensibilités

- **L'exposant** dans la fenêtre des données (6,2758 à 6,2926) : part statique 757–771 MeV,
  √σ 427–435 MeV, partie forte de n − p 2,45–2,50 MeV, total 1,26–1,28 MeV : tout bouge de moins
  de 2 %. L'exposant n'est pas ce qui décide.
- **Le compte 9** : 8 donnerait 364 MeV et √σ = 205 MeV, 10 donnerait 1479 MeV et 834 MeV
  (facteur 3^{2π/…} par brin) : c'est l'entrée qui pèse. Elle vient du contenu uud/udd à trois
  brins par quark, pas d'un ajustement.
- **La largeur du ruban** à un facteur 2 : n − p de 1,13 à 1,41 MeV (−12 à +9 %) ; le logarithme
  se compense presque entre u et d. L'écart n − p ne teste w qu'à un facteur 2 près.
- **Le rayon de l'anneau**, 0,56 à 0,67 fm : part circulante 176 à 147 MeV, N de 938 à 910 MeV.
  C'est le seul nombre non fixé qui pèse sur m_p ; √3 fois le rayon de circuit (0,672 fm) et
  l'enveloppe de trois circuits tangents (0,836 fm = r_p à 0,6 %) sont deux coïncidences de la
  même étoile, sans mécanisme (R49).

## 7. Ce qui est posé, ce qui est coïncidence, ce qui est ouvert

- Posé : le quantum par circuit ; les charges des brins ; le partage 2/3, 2/3, −1/3 ; la forme du
  pôle ; la section carrée ; le compte 9 ; l'exposant (reparamétré).
- Coïncidences enregistrées, non revendiquées : r_p comme enveloppe de trois circuits ;
  m_p c² = 4ℏc/r_p ; l'anneau à √3 rayons de circuit ; le pion à 2m_e/α.
- Ouvert : ce qui fixe le rayon de l'anneau (et donc m_p à 3 % près) ; ce qui épingle l'anneau
  fermé du pion à r_e/2 ; pourquoi 9 et pourquoi 2π ; μ_p exige 0,87 e en circulation au rayon de
  0,67 fm, pas une charge propre (R49) ; les autres lectures de la part circulante donnent le
  neutron de +1,5 à +7 % (R41) ; la section w n'a plus de dérivation depuis que la représentation
  de l'électron à trois jonctions est fermée (R103) : elle reste une entrée à un facteur 2 près.
- Fermé en route : nucléons et quarks sur l'échelle des leptons ; masses par comptes de brins ;
  « un joint = un nombre » ; l'enroulement comme force ; le tube de flux ; la densité de charge
  fixe (R18–R31, R33–R39).

## 8. Ce qui tuerait le secteur

- Une détermination sur réseau de la partie QCD de m_n − m_p qui s'écarte de 2,47 MeV de plus que
  son ±0,29 actuel, ou de √σ hors de 420–440 MeV de plus que les 2 % que l'exposant permet : les
  formules (2α/3)·3^{2π}·m_e c² et (3/(2√π))·3^{2π}·m_e c² n'ont aucun nombre libre.
- Un écart Δ − N ou un rapport de moments que plus aucun rayon d'anneau unique ne reproduit à
  15 % et 3 %.
- Une règle qui fixerait le rayon de l'anneau et donnerait m_p à mieux que 3 % serait le premier
  vrai test prédictif ; jusque-là le secteur est un ensemble de rapports à 1–3 % avec deux entrées
  structurelles.

## 9. Fichiers

- `hadron_sector.py` : la chaîne entière en un calcul, 6/6.
- Étapes d'origine : `closed_ring_scale.py` (R16), `quark_star.py` (R20), `quark_y.py` (R22),
  `junction_valence.py` (R24), `nucleon_moments.py` (R25), `three_strand.py` (R28),
  `mode_rule.py` (R29), `pitch_search.py` (R39), `circuit_rule.py` (R40),
  `circuit_partition.py` (R41), `pn_splitting.py` (R42), `coulomb_geometry.py` (R43),
  `log_constant.py` (R44), `partner_screening.py` (R45), `ribbon_width.py` (R46),
  `delta_shape.py` (R48), `ring_radius.py` (R49), `pole_shape.py`, `nucleon_single_count.py` (R72).
