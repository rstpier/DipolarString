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
que la CD devient une construction connue de la théorie des cordes, où les nombres sont choisis et non dérivés. Depuis le
choix du Calabi-Yau (§9), elle peut nommer l'espace caché ; la conclusion ne change pas.

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

## 8. Trois groupes de trois (précisions du 22 novembre, R129)

L'auteur : « 9 dim : 3 groupe de 3 dimensions », puis, à la question de savoir si c'est le même espace que les trois
tores du §7 : « tores, non probablement pas ». Audit : `ds_three_groups.py`, 9/9.

**Si chaque groupe caché est fermé et plat**, rien ne change : six dimensions regroupées par paires (les trois tores
du §7) ou par triplets (le groupe capacitif, le groupe inductif), c'est le même espace. Ce que le regroupement par
trois apporte :

- les objets qui remplissent des groupes entiers ont 3, 6 ou 9 dimensions d'espace (D3, D6, D9) ;
- retourner un groupe entier (dualité T sur ses trois directions) échange dans ce groupe impulsion et enroulement,
  l'« inductif » et le « capacitif », et fait passer du type IIA au type IIB ;
- deux objets alignés sur les groupes se croisent 0 ou 1 fois : une génération au plus. Trois générations demandent
  un objet en biais entre une direction capacitive et une inductive (enroulement 2 au moins), ou une symétrie
  d'ordre 3 ;
- cette symétrie d'ordre 3 a 3 points fixes par tore, donc 3, 9, 27 : c'est la place du « 3ⁿ » de l'auteur ;
- même groupe : paires non chirales, annihilation et masse permises ; groupes différents : fermion chiral.

**Si les groupes ne sont pas des tores**, la forme d'un groupe de trois dimensions n'est pas libre, parce qu'en trois
dimensions la courbure est fixée entièrement par le contenu (6 composantes de courbure, 6 équations). Trois cas :

| Forme d'un groupe caché | Ce qu'il faut | Ce que cela coûte | Statut |
|---|---|---|---|
| Plat et fermé | Rien | Rien. Un 3-tore ou l'un de ses neuf replis (ordres 2, 3, 4, 6) | Seul cas gratuit en produit |
| Sphère | Un flux qui la tient, rayon R² = kα′ | L'espace visible paie : couplage multiplié par e tous les R/√2, ou courbure visible de même rayon | Exclu |
| Infini | Rien | Gravité en 1/r⁸, pas d'orbite stable | Exclu |

La seule solution exacte connue à trois groupes de trois courbes est AdS₃ × S³ × S³ × S¹ : 3 + 3 + 3 + 1, mais son
groupe visible n'a que deux dimensions d'espace et il est courbe au rayon des sphères. Le compte de l'auteur avec
trois sphères et le temps n'a pas de solution statique.

**L'issue sans coût** est de renoncer au produit : un espace de six dimensions d'un bloc peut être courbe sans rien
coûter, c'est une variété de Calabi-Yau. Elle garde quelque chose des deux groupes : on peut la voir comme une famille
de 3-tores (un groupe) posée au-dessus d'une base de trois dimensions (l'autre groupe), la forme du tore variant d'un
point à l'autre de la base (Strominger, Yau et Zaslow 1996, conjecture) ; retourner le groupe torique y est la
symétrie miroir. Les deux groupes y sont entrelacés, et l'un reste torique. Le nombre de générations y est un nombre
topologique de l'espace : 100 pour la quintique, 3 pour le quotient d'ordre 3 de la variété de Tian et Yau (χ = −6).
Il se choisit avec l'espace.

**Verdict.** « Trois groupes de trois » est compatible avec le §7. « Pas des tores » n'est possible qu'en entrelaçant
les deux groupes cachés ; les sphères et les groupes infinis sont exclus. Dans tous les cas le 3 des générations reste
un choix : paramétrisation. Hors champ de cet audit : les flux de Ramond et les plans d'orientifold.

## 9. Calabi-Yau (décision du 22 novembre, R130)

L'auteur, devant le §8 : « ok prend Calabi-Yau ». Définition retenue : **supercorde de type IIA sur une variété de
Calabi-Yau, avec des D6-branes enroulées sur des cycles de dimension trois** ; la variété est vue comme une famille de
3-tores (le groupe torique) au-dessus d'une base de trois dimensions (le groupe de base). Exemple de travail : la
quintique, la plus simple. Audit : `ds_calabi_yau.py`, 7/7.

| Question | Réponse calculée | Statut |
|---|---|---|
| Que coûte la courbure ? | Rien : à six dimensions, 105 composantes de courbure pour 21 équations | Dérivé |
| Que reste-t-il de la supersymétrie ? | Un quart : 8 charges sur 32, 4 avec des branes ; le tore gardait tout | Dérivé |
| Que deviennent les deux groupes cachés ? | Une base et un tore dont la forme varie et se pince ; ils ne sont plus interchangeables | Identité de structure |
| Que vaut l'échange capacitif ↔ inductif ? | Retourner le groupe torique : c'est la symétrie miroir, qui échange tailles et formes | Identité |
| Que peut-on calculer ? | Des entiers. Quintique : 1 taille, 101 formes, 204 cycles de dimension trois, nombre d'Euler −200 | Calculé, à refaire à la main |
| D'où viennent les générations ? | Du produit d'intersection de deux branes dans le réseau des cycles : tout entier est disponible | Entier choisi |
| Et la règle du nombre d'Euler ? | Elle appartient à la corde hétérotique, sans corde ouverte, donc sans dipôle ; et −6 n'apparaît jamais dans un échantillon de 2362 espaces | Non disponible |
| Les masses ? | Elles dépendent de la forme exacte de l'espace, qui existe mais que personne ne sait écrire | Non calculable |
| Le prix ? | 52 à 410 champs sans masse selon l'espace, à fixer par un mécanisme de plus | Choix supplémentaires |

- **Les entiers.** Trois comptes indépendants concordent sur dix variétés classiques : les formes par comptage des
  déformations, le nombre d'Euler par les classes de Chern, les tailles par différence. Pour la quintique : 126 monômes
  de degré 5 moins 25 changements de coordonnées, soit 101 formes ; nombre d'Euler −200 ; une taille.
- **Les générations.** Le dipôle exige des cordes ouvertes, donc des branes ; le nombre de générations est alors le
  produit antisymétrique de deux branes dans le réseau des cycles de dimension trois. Sur le 6-tore c'est la formule
  du §7 (un déterminant, réseau de rang 20) ; sur la quintique le réseau a le rang 204. Tout entier s'y obtient.
- **Complément au §7.** Dans un espace fermé les lignes de force n'ont nulle part où aller : la charge totale des
  branes est nulle sur chaque cycle. Deux familles de branes seules ont alors un nombre d'intersection nul. L'exemple
  à trois générations du §7 ne ferme donc pas seul son bilan : il lui faut d'autres branes ou des plans d'orientifold.
- **Le prix.** 6-tore : 70 champs scalaires sans masse, dont 22 de taille, de forme et de couplage. Quintique : 410,
  dont 204. Quotient de Tian et Yau : 52, dont 25. Chacun porterait une force de portée infinie ; il faut leur donner
  une masse supérieure à 3,8 meV par des flux, dont les choix entiers font le « paysage ».
- **Le Calabi-Yau le plus proche du §7** est celui des trois tores repliés par leur symétrie d'ordre 3 : la courbure
  y est concentrée aux 27 = 3³ points fixes, il n'a aucun paramètre de forme et deux cycles de dimension trois
  seulement. Mais la rotation envoie la direction capacitive de chaque tore sur l'inductive : les deux groupes cachés
  de l'auteur n'y font plus qu'un. Trois générations y sont fréquentes aux petits enroulements (34 %, puis 16 %, puis
  2 %) : effet de plage, comme au §7.

**Verdict.** Le Calabi-Yau est un choix cohérent, sans coût pour l'espace visible, qui garde deux groupes cachés. Il ne
donne aucun nombre mesuré : l'espace se choisit, parmi au moins 30 108 types connus ; le 3 des générations est un entier
choisi ; les masses ne sont pas calculables ; et le choix ajoute des dizaines à des centaines de champs à fixer. C'est
une paramétrisation, plus large que celle du §7. Le critère d'arrêt (1) du §7 reste atteint.

## 10. La chaîne de maillons (décision du 22 novembre, R131)

L'auteur demandait si les cordes peuvent se lier en longueur, en fil. Oui, bout à bout ; et son filament de dipôles en
chaîne a un équivalent connu. L'auteur : « ok prend la chaîne de maillons ». Définition retenue : la corde du §9
écrite comme une **chaîne de M maillons**, chacun portant la même part d'impulsion et lié à ses voisins par un ressort.
C'est le modèle des « bits de corde » de Giles et Thorn (1977). Audit : `ds_chain_of_links.py`, 6/6.

| Motif de la CD | Sa place dans la chaîne | Statut |
|---|---|---|
| Le dipôle élémentaire | Un maillon : une charge sur un support, la charge opposée sur un autre | Identité |
| Le filament dipolaire | La chaîne : les charges voisines s'annulent, il reste un plus et un moins aux bouts | Identité |
| Deux fils qui se lient en longueur | La mise bout à bout de deux chaînes ; le nombre de maillons se conserve | Identité |
| Les cellules L et C de la ligne | La chaîne elle-même : même dispersion, même impédance | Identité |
| L'adaptation d'impédance | L'égalité des tensions : l'impédance de la chaîne est la tension | Identité |
| La pente de Regge | Les modes lents de la chaîne, équidistants | Retrouvé |
| Le 1/12 qui condamnait la corde à quatre dimensions | L'énergie de point zéro de la chaîne | Retrouvé |
| La taille du maillon | Aucune : un maillon de taille finie fixe un repère et un axe | Exclu |

- **La chaîne est la corde.** Ses fréquences sont (2T₀/m) sin(nπ/2M). Pour les modes lents elles sont équidistantes :
  c'est le spectre de la corde, avec la pente α′ = 1/(2πσ) = 0,86 GeV⁻².
- **Le point zéro.** L'énergie de point zéro de la chaîne vaut 4M²/π − M − π/12 − … en unités de la tension. Les deux
  premiers termes demandent deux contre-termes ajustés ; le troisième donne l'ordonnée (D − 2)/24. Avec la
  supersymétrie tout s'annule.
- **Lorentz.** Avec un nombre fini de maillons, deux niveaux qui doivent coïncider se séparent dans le rapport
  cos(π/2M). Cinq maillons par fermi, la taille des longueurs de la base : 5 % d'écart. Les tests d'isotropie en
  demandent plus de 10¹⁴ par nucléon. Et le nombre de maillons dépend du repère : un proton du LHC en porte 14 500
  fois plus que le même proton au repos. Le maillon n'est pas un constituant, c'est une tranche aussi fine qu'on veut.
- **Où vivent les maillons.** Dans huit dimensions, deux visibles et six cachées ; la troisième dimension visible est
  portée par leur nombre.
- **À refaire à la main.** Les fréquences des chaînes de 2, 3 et 4 maillons : √2 ; 1 et √3 ; 0,765, 1,414 et 1,848.

**Verdict.** La chaîne de maillons est une écriture connue de la corde, exacte quand les maillons sont infiniment fins.
Elle donne sa place exacte au filament dipolaire et aux cellules LC de la CD d'origine. Elle n'apporte aucun nombre
nouveau, et elle ferme une porte : le maillon ne peut pas avoir la taille d'une longueur de la base. C'est la même
cause que R123 et R124 : toute structure discrète de taille finie fixe un repère.

## 11. Le maillon infiniment fin (décision du 22 novembre, R132)

L'auteur, devant le §10 : « ok prend le maillon infiniment fin ». C'est la limite continue de la chaîne : le nombre de
maillons tend vers l'infini, la part de chacun vers zéro. Audit : `ds_fine_link_limit.py`, 6/6.

| Question | Réponse calculée | Statut |
|---|---|---|
| La limite existe-t-elle ? | Oui : l'écart de la chaîne à la corde tombe comme 1/M² | Dérivé |
| À quelle condition respecte-t-elle la relativité ? | Le premier niveau doit être sans masse : 26 dimensions, ou 10 avec des fermions | Dérivé, classique |
| Que devient le choix de dix dimensions du §7 ? | Une conséquence, pour une corde valable à toute longueur | Dérivé |
| Et à quatre dimensions ? | Premier niveau à 1033 MeV avec 2 états pour 3 ; seule une corde effective, longue, reste permise | Exclu comme corde fondamentale |
| Que coûte la limite ? | Deux contre-termes ajustés exactement, ou la supersymétrie, qui les rend inutiles | Dérivé |
| Que devient le dipôle élémentaire ? | Son moment tend vers zéro ; le fil est une ligne uniformément polarisée, charges aux deux bouts | Identité |
| Quelle est la taille du fil ? | Aucune taille propre : α′(ln M + 0,82) par dimension, elle croît avec la résolution | Dérivé |
| Quelle longueur reste ? | ℓ_s = ℓ₁(9)/(π√2) = 0,183 fm, amplitude des vibrations, pas taille d'un élément | Identité |

- **Dix dimensions, à refaire à la main.** Le premier niveau de vibration a D − 2 états. C'est le nombre d'états d'un
  vecteur sans masse ; un vecteur massif en aurait D − 1. Il doit donc être sans masse. Or sa masse au carré vaut
  1 − (D − 2)/24 en unités de 1/α′ : zéro pour D = 26. Avec un fermion par dimension elle vaut 1/2 − (D − 2)/16 : zéro
  pour D = 10.
- **La taille.** Chaque mode n ajoute α′/n au carré de la taille : la somme croît comme le logarithme du nombre de
  maillons. De mille à 10³⁰ maillons, 0,5 à 1,5 fm par dimension avec la tension de la base : l'ordre de grandeur d'un
  hadron, sans valeur définie. C'est la loi connue de l'élargissement du tube de flux et du rétrécissement de Regge ;
  le seul nombre qui y entre est la tension.
- **L'élément.** Dans cette lecture, la question d'origine du projet, caractériser l'élément du filament, reçoit sa
  réponse : il n'y a pas d'élément. Le rapport D/r = 2 cosh π de la CD d'origine n'y a pas de place.

**Verdict.** Le maillon infiniment fin est la corde continue. La limite n'est cohérente, pour une corde valable à toute
longueur, qu'à dix dimensions avec des fermions : le choix du §7 en découle. Le fil est fixé par sa tension seule ;
aucun nombre nouveau.

## 12. La corde effective (décision du 22 novembre, R133)

L'auteur, devant le §11 : « ok prend la corde effective ». C'est la corde à quatre dimensions, valable seulement quand
elle est longue : la description établie des tubes de flux entre un quark et un antiquark. C'est aussi le chemin B de
la feuille de route, dont R133 fait le premier calcul. Audit : `ds_effective_string.py`, 6/6.

| Question | Réponse calculée | Statut |
|---|---|---|
| Est-elle cohérente à quatre dimensions ? | Oui, pour une corde longue devant 0,33 fm | Acquis, classique |
| Que couvre-t-elle ? | Les hadrons, entre 0,33 fm et la brisure vers 1,25 fm ; ni leptons, ni photon, ni gravité | Limite du choix |
| Le ρ et le nucléon ? | Longs de 0,45 et 0,56 fm : au bord du domaine ; la « part statique » de la base n'est pas contrôlée | Hors domaine |
| Le calcul B1 reproduit-il la littérature ? | Oui : pente commune de six familles 0,883 ± 0,014 GeV⁻², publiée 0,884 | Fait, à refaire par l'auteur |
| La tension de la base est-elle confirmée par les trajectoires ? | Elle est dans leur enveloppe, 426 ± 6 MeV ; mais la pente varie de 20 % d'une famille à l'autre | Non concluant |
| Et par le réseau ? | 445 ± 7 MeV : la base, 430,3 MeV, est à deux écarts types en dessous | Non concluant |
| Les deux tensions sont-elles la même ? | Écart de 19 MeV, deux écarts types | Non tranché |

- **Le calcul.** Une corde qui tourne avec une masse à chaque bout : σℓ = mγ²β². On ajuste la pente, avec une
  ordonnée par famille, sur six familles de mésons et deux de baryons.
- **Ce que le calcul apprend.** La pente dite universelle ne l'est qu'à 20 % près : de 0,74 pour η/h à 0,92 pour ω/f.
  Selon les familles retenues, la base est compatible ou trop haute en tension de 2 à 4 écarts types. R120 est à
  nuancer dans ce sens.
- **Ce qui ne dépend pas de la CD.** Le terme de Lüscher, −51,7 MeV·fm/r, et la forme du potentiel. À 0,5 fm le
  potentiel mesuré sur le réseau n'est pas celui d'une corde seule : la tension qu'on en tire dépend de la fenêtre.
- **Le raccord avec le §7.** La corde effective est la partie horizontale d'une corde fondamentale au fond d'une
  géométrie déformée ; les masses d'extrémité en sont les segments verticaux.

**Verdict.** La corde effective est de la physique établie et testable, mais elle réduit le projet aux hadrons, et la
CD n'y apporte qu'un nombre, la tension, ni confirmée ni exclue. Restent à faire : B1 par l'auteur, B2 avec des tables
de réseau, la recherche bibliographique de B3.

## 13. Choix faits par délégation (22 novembre, R134 à R138)

L'auteur : « ok prend la tension du réseau », puis : « Puisque tu connais les meilleures réponses aux questions que tu
me poses, fais les choix à ma place. » Les choix de physique sont désormais faits par l'assistant, consignés avec leur
raison dans `../redraw/BASE.md`, et révocables. La publication, le dépôt, la mention publique et le calendrier restent
à l'auteur, ainsi que la règle 1 : refaire soi-même ce que l'on revendique.

| Choix | Décision | Raison |
|---|---|---|
| La tension de la corde effective | Celle du réseau, 445 ± 7 MeV | C'est la pente du potentiel à grande distance, que le réseau mesure directement |
| La formule de la base pour la tension | Retirée comme prédiction, gardée au registre | Aucun mécanisme ; accord non significatif (R126) ; 2π devrait valoir 6,314 ± 0,014 |
| Le rôle des deux lectures | Dix dimensions : étude. Quatre dimensions : tout calcul confronté aux mesures | Seule la corde effective est testable |
| La table des hadrons de la base | Pas de remise à l'échelle | Ce serait de la numérologie ; la « part statique » est hors domaine |
| La suite | Lire avant de calculer ; aucune nouvelle adoption | Chaque choix de plus serait une paramétrisation |
| La piste de la détente du fil | Calculée (R138), puis fermée | Elle déplace le désaccord sans le réduire : 2,9 à 3,0 écarts types, comme sans elle |

**Ce que la tension du réseau entraîne (R134, `ds_lattice_tension.py`, 6/6).** Le fil pèse 1004 MeV par fermi au lieu
de 938. Dans la corde effective il ne reste aucun nombre propre à la CD. Et l'écart change de côté : avec cette
tension la pente des trajectoires sort trop faible de 9 à 11 %, à 2,7 et 3,6 écarts types. Aucune valeur ne réconcilie
les deux mesures ; la meilleure valeur commune serait 434 ± 5 MeV.

**L'état de la question B3 (R135, `veille_b3.py`, 6/6).** Côté spectres la littérature donne σ = 0,180(2) GeV², soit
424 ± 2 MeV ; côté réseau 445(3)(6) MeV. Écart de 21 MeV, deux à trois écarts types selon le budget d'erreur. La
recherche n'a trouvé aucune source qui confronte directement les deux. Une piste publiée existe : le potentiel
s'aplatit à grande distance quand des paires se créent, ce qui abaisse les pentes. Elle est à lire avant tout calcul.
La recherche a aussi confirmé deux valeurs prises de mémoire en R133 : la brisure à 1,22 et 1,29 fm, et r₀ = 0,472 fm.

**Lecture à la source (22 novembre).** Tentée et bloquée : les sites des articles ne sont pas lisibles depuis la session
en cours. La pente publiée 0,884 GeV⁻², les masses d'extrémité et les masses des familles autres que ρ/a, N et Δ
restent donc citées de mémoire. Les résumés donnent un élément de plus : selon Badalian et Bakker, l'aplatissement du
potentiel abaisse les pentes d'environ 30 %, le sens et l'ordre de grandeur qu'il faut pour l'écart de 9 à 11 %.

**Lecture à la source, faite (R136, `ds_b1_reproduction.py`, 6/6).** L'auteur a autorisé le téléchargement depuis
arXiv. L'ajustement de l'article de 2014 est refait à l'identique : neuf familles, 38 états, pente 0,8848 contre
0,884 publié, même mesure d'écart, 13,13·10⁻⁴, et les 38 masses calculées à 4 MeV près. La tension des spectres vaut
424 ± 2,5 MeV. Celle de la base en est à 2,5 écarts types, celle du réseau à 2,9. Côté réseau, la tension dépend de sa
définition : 445 MeV dans le modèle qui traite la brisure, 467 à 482 MeV dans des ajustements du potentiel sous 1 fm.
Et la piste publiée n'est pas une réponse : elle vient d'un autre modèle, qui garde la tension des spectres et ne la
compare pas au réseau. La question B3 se pose donc des deux côtés.

**Masses du PDG 2026, lues à la source (R137, `ds_pdg2026_update.py`, 5/5).** L'auteur a autorisé le téléchargement
depuis le PDG. Sur trente et un états cités de mémoire, vingt-huit sont identiques aux tables ; une correction compte,
celle du Δ(2420), estimé à 2450 MeV, ce qui ramène la pente de sa famille de 0,92 à 0,89, celle des mésons. Avec les
masses de 2026 l'ajustement ne bouge pas : tension des spectres 423,7 ± 2,5 MeV, base à 2,6 écarts types, réseau
à 3,0. La dernière masse a été relue ensuite : D₂(2740)⁰, 2747 ± 6 MeV. Plus aucune valeur n'est citée de mémoire.

**La détente du fil appliquée à la corde tournante (R138, `ds_string_relaxation.py`, 8/8).** Le long d'un fil long, des
paires de quarks apparaissent et le fil tire moins fort : c'est la piste publiée pour l'écart entre les deux tensions.
Critère écrit avant le calcul : l'écart, 3,0 écarts types, doit tomber sous 2. Relue, la fonction publiée est faible
aux longueurs des mésons : 1,9 % à 1,18 fm, 4,4 % à 1,6 fm ; les 40 % ne sont atteints qu'au-delà de 3 fm. Les « 30 % »
des résumés valent dans le modèle de ses auteurs, pas pour la corde tournante.

| | Tension demandée par les spectres | Perte d'ajustement | Écart complet au réseau |
|---|---|---|---|
| Sans détente | 423,7 MeV | 0 | 2,97 écarts types |
| Détente publiée, force égale à la tension | 430,1 MeV | 1,95 écart type | 2,86 |
| Détente publiée, force égale à la dérivée de l'énergie | 434,2 MeV | 2,54 écarts types | 3,01 |

- **Ce que le calcul apprend.** La détente déplace le désaccord sans le réduire : ce qu'elle gagne sur la tension, 6 à
  10 MeV sur 21, elle le perd sur l'ajustement. Sur toute la fourchette publiée l'écart complet reste au-dessus de 2,8.
  La seconde version est fragile : l'a₆ et le f₆ y sont au bout de la branche, et son écart monte à 3,8 si le profil
  monte 1,8 % plus vite.
- **Les spectres ne demandent pas de détente.** L'ajustement se dégrade dès la plus petite valeur ; à deux écarts types
  ils tolèrent γ jusqu'à 0,36 à 0,41, la valeur publiée étant 0,40.
- **Un fait nouveau, qui n'est pas de la CD.** Si la force est la dérivée de l'énergie, un fil tournant n'a plus de
  rotation stable au-delà de 2,1 fm : la trajectoire des mésons légers s'arrête à J = 6, le plus haut spin observé.
  Cela dépend d'un profil que ses auteurs disent arbitraire.

**Verdict.** La chaîne des adoptions est close. Il reste de la physique établie, les tubes de flux, et une question
ouverte, B3, pour laquelle la CD n'apporte aucun nombre. La seule piste publiée pour B3, la détente du fil, ne tient
pas pour la corde tournante (R138). Reste B2, avec les niveaux publiés du réseau.

## 14. Suite

Les chemins A et B de la feuille de route, avec sa règle 1 : ne rien revendiquer qu'on ne peut pas refaire soi-même.
Les deux premiers calculs à refaire à la main sont R127 A (Zwiebach, chapitres 7 et 8) et R127 B (chapitre 12). Pour le
§9 : le compte 126 − 25 = 101 et le nombre d'Euler −200 de la quintique. Pour le §10 : les fréquences des chaînes de 2,
3 et 4 maillons. Pour le §11 : 1 − (D − 2)/24 = 0 et 1/2 − (D − 2)/16 = 0.

*Références.* Vérifiées par recherche (titre, revue) : L. E. Ibáñez, F. Marchesano et R. Rabadán, JHEP 11 (2001) 002,
hep-th/0105155 ; J. G. Lee et al., Phys. Rev. Lett. 124, 101101 (2020) ; CMS, arXiv:1911.03947, paru dans JHEP
en 2020 (la valeur de la borne sur les résonances de cordes, 7,9 TeV de mémoire, n'a pas pu être confirmée).
*Données de mémoire, à vérifier à la source :* J. Polchinski et A. Strominger, Phys. Rev. Lett. 67, 1681
(1991) ; A. M. Polyakov, Nucl. Phys. B 268, 406 (1986) ; H. Kleinert, Phys. Lett. B 174, 335 (1986) ; J. F. Arvis,
Phys. Lett. B 127, 106 (1983) ; B. Zwiebach, *A First Course in String Theory*, 2ᵉ éd. (2009).
*Pour le §8, de mémoire, à vérifier :* L. Bieberbach, Math. Ann. 70, 297 (1911) ; W. Hantzsche et H. Wendt, Math. Ann.
110, 593 (1935) ; C. Callan, J. Harvey et A. Strominger, Nucl. Phys. B 359, 611 (1991) ; S. Elitzur, O. Feinerman,
A. Giveon et D. Tsabar, Phys. Lett. B 449, 180 (1999) ; A. Strominger, S.-T. Yau et E. Zaslow, Nucl. Phys. B 479, 243
(1996) ; P. Candelas, G. Horowitz, A. Strominger et E. Witten, Nucl. Phys. B 258, 46 (1985) ; G. Tian et S.-T. Yau
(1987).
*Pour le §9, de mémoire, à vérifier :* S.-T. Yau, Proc. Natl. Acad. Sci. 74, 1798 (1977) ; P. Candelas, A. Dale,
C. Lütken et R. Schimmrigk, Nucl. Phys. B 298, 493 (1988) ; M. Kreuzer et H. Skarke, Adv. Theor. Math. Phys. 4, 1209
(2000) ; B. Greene, K. Kirklin, P. Miron et G. Ross, Nucl. Phys. B 278, 667 (1986) ; L. Dixon, J. Harvey, C. Vafa et
E. Witten, Nucl. Phys. B 261, 678 (1985) ; R. Blumenhagen, B. Körs, D. Lüst et T. Ott, Nucl. Phys. B 616, 3 (2001) ;
M. Gross, Invent. Math. 144, 75 (2001).
*Pour le §10, de mémoire, à vérifier :* R. Giles et C. B. Thorn, Phys. Rev. D 16, 366 (1977) ; I. Klebanov et
L. Susskind, Nucl. Phys. B 309, 175 (1988) ; C. B. Thorn, hep-th/9405069 ; O. Bergman et C. B. Thorn, Phys. Rev. D 52,
5980 (1995) ; M. Smiciklas et al., Phys. Rev. Lett. 107, 171604 (2011).
*Pour le §11, de mémoire, à vérifier :* P. Goddard, J. Goldstone, C. Rebbi et C. B. Thorn, Nucl. Phys. B 56, 109
(1973) ; L. Brink et H. B. Nielsen, Phys. Lett. B 45, 332 (1973) ; M. Karliner, I. Klebanov et L. Susskind, Int. J.
Mod. Phys. A 3, 1981 (1988) ; M. Lüscher, G. Münster et P. Weisz, Nucl. Phys. B 180, 1 (1981).
*Pour le §12, de mémoire, à vérifier :* J. Sonnenschein et D. Weissman, JHEP 08 (2014) 013 et JHEP 02 (2015) 147 ;
M. Lüscher et P. Weisz, JHEP 07 (2002) 049 ; R. Sommer, Nucl. Phys. B 411, 839 (1994) ; G. Bali et al., Phys. Rev. D
71, 114513 (2005) ; masses des états autres que ρ/a, N et Δ.
*Pour le §13, vérifiées par recherche le 22 novembre :* J. Bulava et al., Phys. Lett. B 854 (2024) 138754,
arXiv:2403.00754 ; A. M. Badalian et B. L. G. Bakker, arXiv:1901.10280 ; J. Sonnenschein et D. Weissman, JHEP 06 (2018)
148, arXiv:1801.00798 ; FLAG Review 2021, arXiv:2111.09849. *Lus à la source le 22 novembre (R136) :* arXiv:1402.5603,
arXiv:2403.00754, arXiv:1901.10280. *Pour R138 :* arXiv:1901.10280, éq. 51 à 54, relues à la source.
