# La CD comme variante de la théorie des cordes — définition, dictionnaire, premier audit

**22 novembre 2026.** Décision de l'auteur : « On va transformer dipolar string comme une variante de la TdC », avec
la mention à porter : « 22/11/2026 — Dipolar String change de direction ». Ce document fixe ce que cette phrase peut
vouloir dire, donne le dictionnaire, les résultats du premier audit (R127, `ds_as_string_variant.py`, 6/6) et les
critères d'arrêt. Il complète la feuille de route (`../../FEUILLE_DE_ROUTE_2026-11-22.md`).

## 1. Trois sens de « variante », et lesquels sont ouverts

1. **Un dictionnaire.** Relire chaque objet de la CD comme un objet de la théorie des cordes. Ouvert, et fait au §3.
2. **Une corde effective.** Le filament comme corde longue de tension donnée : ce sont les tubes de flux de la QCD.
   Ouvert : c'est le chemin B de la feuille de route.
3. **Une nouvelle théorie des cordes fondamentale.** Il faudrait une théorie de feuillet d'univers cohérente :
   anomalie conforme annulée, invariance modulaire, pas de tachyon. Les théories cohérentes connues sont classées
   (la corde bosonique à 26 dimensions, les cinq supercordes à 10). Non ouvert : hors de portée, et R127 ferme la
   version naïve à quatre dimensions.

## 2. La définition retenue

Le filament dipolaire devient une corde relativiste :

- sa loi élastique est l'action de Nambu–Goto, de tension σ, avec un terme de rigidité optionnel (Polyakov, Kleinert,
  1986) ; c'était l'objet réservé sous le numéro R127 ;
- il porte des charges +q et −q à ses extrémités : c'est une corde ouverte à charges opposées, vue de l'extérieur un
  dipôle ;
- il est quantifié, ce qui rend l'axiome d'amplitudes (A9) automatique.

Ce qu'on abandonne de la CD d'origine, déjà exclu par le dossier : le réseau au pas ℓ₁ (R123), le repère du substrat
et la gravité scalaire (R124), les objets sombres (R125), l'électron étendu (R100–R103).

## 3. Le dictionnaire

| Dans la CD | En théorie des cordes | Statut | Zwiebach |
|---|---|---|---|
| Filament dipolaire, pôles + et − | Corde ouverte à charges opposées aux extrémités ; un dipôle | Identité de structure | 15, 16 |
| Signal (V, I) le long du filament | Champs X(τ, σ) sur le feuillet : chaque coordonnée de la corde est une ligne de transmission | Identité (équation d'onde à 1 + 1 dimensions) | 4, 6 à 9 |
| Bout ouvert, bout court-circuité | Conditions de Neumann, de Dirichlet (D-brane) | Identité | 6, 15 |
| Échange L ↔ C, Z → Z_Q²/Z | Dualité T, R → α′/R | Identité conditionnelle (phase compacte) | 17, 18 |
| α = Z₀/(2R_K) | Rayon de compactification R²/α′ = 1/(2α) = 68,5 pour q = e | Identité ; pas une dérivation de α | 17 |
| Tension σ de la base | T = 1/(2πα′) | Entrée | 6 à 8 |
| « Un quantum par circuit », J = ½ | Trajectoire de Regge J = α′M² + a | Identité (R120) | 7, 8, 12 |
| Terme de Lüscher, point zéro | Ordonnée a = (D − 2)/24 et dimension critique | Identité | 12 |
| « L'énergie change le vide » (R122) | Électrodynamique de Born–Infeld sur la brane, champ critique E_c = T/q | Vrai en cordes ; exclu pour le vide réel (R122) | 19, 20 |
| Charge = brin non apparié | Corde tendue entre deux branes différentes | Analogie | 15, 21 |
| Tiers de charge, trois brins | L'hypercharge combine les U(1) des piles de branes, dont celle de trois branes de couleur | Analogie | 21 |
| Gravité par un champ Z(x) scalaire | Graviton = mode sans masse de la corde fermée | Remplacement | 13 |
| Spin ½ absent | Fermions du feuillet (secteur de Ramond) | Remplacement, au prix de la supersymétrie | 14 |
| Repère du substrat, réseau à ℓ₁ | Aucun repère, aucun réseau ; ℓ_s est la seule longueur | Abandon | — |
| Générations, compte 9, exposant 2π | Rien d'équivalent | Non transposable | — |

## 4. Premier audit (R127)

- **A. Identités.** Avec la tension de la base, α′ = 0,859 GeV⁻² et ℓ_s = 0,183 fm. La part statique de la base est
  ℏc/(√2 ℓ_s) = 763 MeV ; la boucle fermée à un quantum a le rayon ℓ_s, l'anneau rigide le rayon ℓ_s/√2. Toutes les
  longueurs de la CD « corde » sont ℓ_s : aucun nombre nouveau.
- **B. Cohérence quantique.** À quatre dimensions l'ordonnée vaut 1/12 au lieu de 1 : fondamental tachyonique
  (311 MeV), premier niveau massif (1033 MeV) avec deux états là où un vecteur massif en exige trois. L'invariance de
  Lorentz est perdue. Un vecteur sans masse exige 26 dimensions, 10 avec fermions. **Une corde fondamentale à quatre
  dimensions est exclue.** Une corde effective, longue devant ℓ_s, est cohérente en toute dimension (Polchinski et
  Strominger, 1991) ; sa limite de validité est le rayon d'Arvis, (π/√3) ℓ_s = 0,33 fm.
- **C. Gravitation.** Le couplage du graviton est G ~ g_s²α′. Avec la tension de la base il faudrait
  g_s = ℓ_Planck/ℓ_s = 9·10⁻²⁰ ; à couplage d'ordre 1 la gravité serait 10³⁸ fois trop forte. **Une seule corde pour
  les hadrons et la gravité est exclue.** *Rectificatif (R128) : cet énoncé vaut en espace plat. En dix dimensions, une
  géométrie déformée réduit la tension effective de la corde, et la même corde peut faire les deux (§7).*
- **D. Le dipôle.** Champ critique E_c = σ/q : 9,4·10²³ V/m pour q = e, soit 7·10⁵ fois le champ de Schwinger. Le vide
  réel crée des paires bien avant.
- **E. Ligne de transmission.** Une ligne d'impédance Z et de quantum de charge q est un boson compact de rayon
  R²/α′ = h/(q²Z). Pour le vide adapté et q = e : 1/(2α) = 68,5. La dualité T est l'échange L ↔ C ; le point
  auto-dual est Z = R_K = 25,8 kΩ. C'est la lecture exacte de l'idée d'origine, « le vide comme ligne adaptée » : en
  cordes, la ligne est le feuillet, pas l'espace.
- **F. « 10 = 3 + 3 + 3 + 1 ».** Capacitif et inductif correspondent à enroulement et impulsion, deux descriptions
  duales des mêmes dimensions compactes. Le 3 + 3 dédouble trois dimensions, il n'en fait pas six : sept dimensions
  physiques, pas dix. Le compte ne se transpose pas (c'est la réduction de R121). *Rectificatif (R128) : sous dix dimensions il se
  transpose autrement, 3 + 3 = trois tores (§7).*

## 5. Critères d'arrêt, écrits le 22 novembre, et leur état

1. Corde fondamentale à quatre dimensions incohérente : **joué** (R127 B).
2. Une seule corde pour les hadrons et la gravité : joué en espace plat, 10³⁸ (R127 C) ; **levé en dix dimensions** par
   un facteur de déformation (R128 E), au prix d'un nombre choisi.
3. Corde effective : la CD n'y apporte qu'un nombre propre, σ = 0,185 GeV² tiré de 9 et 2π, déjà 2 à 7σ sous le
   réseau (R119). **En suspens** ; les calculs B1 à B3 de la feuille de route tranchent.
4. Si le dictionnaire ne produit aucun énoncé testable distinct de la théorie des cordes ou de la corde effective
   standard, la « variante » est une lecture, pas une théorie. **C'est l'état au 22 novembre.**

## 6. Ce que la mention peut dire

Elle peut dire : la CD d'origine, le vide comme réseau de lignes adaptées, est close par l'erratum V2.11 ; le projet se
réoriente vers la théorie des cordes ; le filament dipolaire y est relu comme corde ouverte à charges opposées ; le
dictionnaire est exact. Elle ne peut pas dire que la CD *est* une variante de la théorie des cordes au sens d'une
théorie nouvelle et cohérente, ni annoncer une prédiction. Depuis le choix de dix dimensions (§7), elle peut ajouter
que la CD devient une construction connue de la théorie des cordes, où les nombres sont choisis et non dérivés.

## 7. Dix dimensions (décision du 22 novembre, R128)

L'auteur, devant l'incohérence à quatre dimensions : « ok prend 10 dimensions ». Dix dimensions, c'est la supercorde :
des fermions sur le feuillet, pas de tachyon, un graviton sans masse. La construction la plus simple qui respecte le
découpage de l'auteur (un temps, trois dimensions « massiques », trois « capacitives », trois « inductives ») est la
**supercorde de type IIA sur trois tores, avec des D6-branes qui se coupent** : c'est le chapitre 21 de Zwiebach.
Audit : `ds_ten_dimensions.py`, 5/5.

| Motif de la CD | Sa place en dix dimensions | Statut |
|---|---|---|
| Le filament dipolaire | La corde ouverte tendue entre deux D6-branes, chargée + sur l'une et − sur l'autre | Identité de structure |
| 3 capacitives + 3 inductives | Les deux directions de chacun des trois tores : six dimensions réelles, trois complexes | Choix de compactification |
| Une brane | Trois dimensions visibles et une direction par tore | Définition d'une D6-brane |
| Les générations | Le nombre d'intersection de deux branes, produit de \|n₁m₂ − m₁n₂\| sur les trois tores | Entier choisi |
| Le spin ½ | 8 états fermioniques sans masse pour 8 bosoniques, puis 128 + 128 | Acquis, au prix de la supersymétrie |
| La tension des hadrons | La tension fondamentale réduite par un facteur de déformation e^{2A} | Nombre choisi |

- **Spectre.** Bosons et fermions en nombre égal à chaque niveau : 8, 128, 1152, 7680, 42112 (identité de Jacobi,
  vérifiée en série).
- **Générations.** La formule d'intersection est vérifiée par comptage direct. Trois intersections s'obtiennent par
  exemple avec les enroulements (1,0)(1,0)(1,0) et (1,3)(1,1)(1,1). Mais 3 n'est pas préféré : sa part parmi les
  paires de branes va de 0 % à 14 % selon la plage d'enroulements. Il se choisit.
- **Échelles.** Avec des tores de la taille de la corde, M_s = √π g_s M_Planck, soit 4·10¹⁷ à 4·10¹⁸ GeV. La borne
  basse vient du LHC, environ 8 TeV (valeur citée de mémoire). La longueur de corde est donc sous 2,5·10⁻²⁰ m : les
  longueurs de la base, 0,81 fm et 809 fm, ne sont pas des longueurs de corde.
- **Tension des hadrons.** Il faut un facteur de déformation e^{A} de 10⁻⁴ (M_s = 8 TeV) à 10⁻¹⁹ (M_s près de
  Planck). C'est le mécanisme de la QCD holographique, et c'est ce qui rectifie R127 C.
- **Non observé.** Aucun partenaire supersymétrique, loi en 1/r² vérifiée jusqu'à 52 micromètres, aucune résonance
  de corde. Rien n'est prédit aux énergies accessibles.

**Verdict.** En dix dimensions la CD devient une construction connue, où ses motifs ont une place exacte. Mais le
groupe de jauge, le 3 des générations et le facteur de déformation y sont choisis, pas dérivés : c'est une
paramétrisation, comme le catalogue de particules retiré en V2.10. Ce n'est pas une variante de la théorie des cordes,
c'en est un cas, parmi un très grand nombre.

Critères d'arrêt pour cette lecture, écrits le 22 novembre : (1) si chaque nombre mesuré exige un entier ou un facteur
choisi, c'est une paramétrisation : **c'est l'état actuel** ; (2) si rien ne distingue la construction des modèles de
branes déjà publiés, elle n'ajoute rien : **à établir par la lecture du chapitre 21** ; (3) règle 1 de la feuille de
route : l'auteur doit pouvoir refaire à la main le compte 8 + 8 et un nombre d'intersection.

## 8. Suite

Les chemins A et B de la feuille de route, avec sa règle 1 : ne rien revendiquer qu'on ne peut pas refaire soi-même.
Les deux premiers calculs à refaire à la main sont R127 A (Zwiebach, chapitres 7 et 8) et R127 B (chapitre 12).

*Références.* Vérifiées par recherche (titre, revue) : L. E. Ibáñez, F. Marchesano et R. Rabadán, JHEP 11 (2001) 002,
hep-th/0105155 ; J. G. Lee et al., Phys. Rev. Lett. 124, 101101 (2020) ; CMS, arXiv:1911.03947, paru dans JHEP
en 2020 (la valeur de la borne sur les résonances de cordes, 7,9 TeV de mémoire, n'a pas pu être confirmée).
*Données de mémoire, à vérifier à la source :* J. Polchinski et A. Strominger, Phys. Rev. Lett. 67, 1681
(1991) ; A. M. Polyakov, Nucl. Phys. B 268, 406 (1986) ; H. Kleinert, Phys. Lett. B 174, 335 (1986) ; J. F. Arvis,
Phys. Lett. B 127, 106 (1983) ; B. Zwiebach, *A First Course in String Theory*, 2ᵉ éd. (2009).
