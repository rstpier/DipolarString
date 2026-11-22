# Revue scientifique du programme « Cordes dipolaires » (CD / DS)

**Rapport d'arbitre interne — 19 novembre 2026**

**Pièces examinées.** Le manuscrit V2.10 (`v2.10/manuscript/DS_model_V2_10_Zenodo.tex`, gelé le
7 septembre 2026, vérification 40/40 + 9 contrôles de source), lu en entier ; `matter/CONSTRAINTS.md`
et `matter/OPEN_PROBLEMS.md` ; le registre du secteur matière (`matter/redraw/BASE.md`, R1–R118,
113 scripts ; `RESUME.md` ; `HADRONS.md`) ; le brouillon `matter/paper/DS_matter_sector.tex` et son
`README.md` ; les notes de `matter/notes/` (septembre–novembre 2026), dont les notes de principe
d'octobre–novembre enregistrées dans `BASE.md` et la note « tout est énergie ondulatoire ».

**Conflit d'intérêts, à lire d'abord.** Ce rapport est rédigé par le même système d'IA qui a produit
la couche formelle du manuscrit (déclaration de l'auteur, section « Acknowledgments and Methodological
Disclosure ») et qui a écrit les 113 scripts du secteur matière. Il ne vaut donc pas revue par un pair
indépendant, et il est exposé au biais de confirmation de ses propres calculs. Il est écrit comme un
rapporteur extérieur l'écrirait, sur pièces, en citant les emplacements ; mais sa conclusion la plus
solide est précisément qu'une relecture par un physicien des particules et par un relativiste humains
est nécessaire et n'a pas eu lieu.

---

## 0. L'axe d'évaluation

L'auteur rappelle (19 novembre 2026) l'objet du programme : **reconstruire les interactions
universelles et les particules élémentaires à partir d'un filament dipolaire hypothétique, dont
l'élément est à caractériser.** C'est l'axe de ce rapport : l'élément (§2), les quatre interactions
(§3), les particules (§4), la structure quantique sans laquelle « particule » n'a pas son sens moderne
(§5), la méthode (§6), ce qui est acquis (§7), les recommandations (§8) et le verdict (§9).

Les critères sont ceux du manuscrit lui-même : une grandeur est **dérivée** si elle sort des axiomes
sans nombre réglé ; **conditionnelle** si elle dépend d'une entrée posée ; **reparamétrée** ou
**coïncidence** si elle réécrit une donnée ou tombe juste sans mécanisme. Le registre matière ajoute
**exclu**. Les problèmes ouverts du manuscrit sont cités par leur rang dans la section « Open Problems
and Limitations » (n° 1 à n° 20).

---

## 1. Verdict en une page

| Secteur | Ce que le programme veut | Ce qui est atteint (V2.10 + registre R1–R118) | Verdict |
|---|---|---|---|
| L'élément (filament / DQD) | Un objet caractérisé par une dynamique | Une géométrie transverse (D/r = 2cosh π) et des règles posées ; pas de hamiltonien, pas d'équation d'état du fluide de charge ; raideur de phase et longueur de cohérence inconnues | **Partiellement caractérisé ; dynamique absente** |
| Électromagnétisme | Dériver Maxwell du filament | Unicité du nœud de Johns sous C1–C5 (Théorème 4) ; c'est Maxwell réécrit en variables de circuit sous des entrées posées, ce que V2.9.2 reconnaît | **Re-dérivation sous entrées, pas une dérivation** |
| Gravitation | Dériver la gravitation du milieu | Secteur scalaire PPN γ = β = 1 ; bootstrap de Deser affirmé, non effectué → extérieur Schwarzschild, aucune prédiction extérieure ; intérieur newtonien non raccordé ; α₁, α₂ non calculés | **RG héritée + un intérieur conjectural** |
| Interaction forte | Dériver le confinement et les hadrons | Rapports hadroniques à 1–3 % sur deux entiers (9, 2π) et une règle posée ; deux formules sans paramètre libre, l'une à 2,2σ, l'autre en tension sur la coupe QCD/QED | **Paramétrisation, pas une interaction** |
| Interaction faible | Dériver W, Z, la violation de P | Rien ; P agit comme C dans le modèle, le secteur faible est hors champ par construction | **Absente** |
| Particules | Reconstruire e, ν, quarks, bosons | Catalogue réduit à photon, DQD, électron ; l'électron étendu est exclu par les facteurs de forme (R100–R103) ; neutrinos bloqués en principe ; W, Z, H retirés ; spin ½ et statistique absents | **Aucune particule reconstruite** |
| Structure quantique | Donner un substrat à la MQ | Pas d'espace d'états, de règle de Born ni de carte de mesure (problème ouvert n° 3) ; équation d'enveloppe de Schrödinger sur un réseau substitut, côté cinématique seul | **Non construite ; l'axiome manquant (A3) est identifié** |

**En trois phrases.** Le programme est un modèle exploratoire du milieu, documenté avec une discipline
rare (grille des statuts, retraits, prédictions enregistrées, scripts PASS/FAIL), qui n'a pas atteint
son objet : aucune interaction n'est dérivée du filament au-delà de ce que les entrées posées
contiennent déjà, et aucune particule n'est construite comme état stable d'une dynamique. Son contenu
scientifique le plus solide est négatif : une série de résultats d'impossibilité internes (pas
d'échelle de masse, pas de soliton 3D stable, pas de charge par connectivité, pas de doublet, pas
d'inversion pour les neutrinos, électron étendu exclu) qui ferment proprement toute une classe de
modèles « vide = réseau de lignes de transmission ». La décision de ne pas publier sans découverte
notable (4 octobre 2026) est la bonne, et ce rapport la confirme.

---

## 2. L'élément : le filament dipolaire

### 2.1 Fiche d'état

| Propriété du filament | Valeur / énoncé | Statut (source) | Appréciation |
|---|---|---|---|
| Rapport transverse D/r | 2cosh π = 23,18 | Dérivé de l'adaptation Z_bif = Z₀ (manuscrit ; R116) | Solide **comme résultat de ligne bifilaire** ; ne fixe qu'un rapport transverse, ni la taille, ni le nombre de brins, ni l'exposant d'échelle (R116) |
| Coefficient d'interface Γ_pole | 1/3 | Rangé « Derived (Predictive) » (manuscrit) | Le coefficient de réflexion peut être dérivé ; son identification à la fraction de charge e/3 est une identification (`CONSTRAINTS.md` §2 : « Γ_pole = 1/3 is a reflection coefficient ») |
| Pas de cellule ℓ₁ | 2πƛ_e/3 = 808,7 fm | Posé (condition initiale cosmologique ; problème ouvert n° 8) | Toute la métrologie du modèle (coupure de Brillouin, E₁ = 1,53 MeV, M_crit ∝ ℓ₁²) est conditionnelle à cette entrée |
| Coupure de Brillouin | ħω_max = (3/π) m_e c² = 0,487 MeV | Dérivé sous ℓ₁ | Le manuscrit enregistre lui-même que les photons observés (2,5 PeV, LHAASO) dépassent la coupure d'un facteur 5·10⁹ et que la dispersion du réseau correspondrait à E_QG,2 ≈ 0,69 MeV contre des bornes ≳ 10¹⁰ GeV. La sortie (photon = soliton N = 1, non collectif) n'est accompagnée d'aucune relation de dispersion calculée pour ce régime : la borne n'est pas évaluée pour le mode qui porte réellement les photons observés |
| Section du ruban w | 4ƛ_e/π² | Entrée à un facteur 2 (R103 ferme la dérivation à trois jonctions) | Entre dans m_n − m_p à ±11 % |
| Fluide de charge e/3 par brin | Équation d'état **non spécifiée** | `CONSTRAINTS.md` §2 ; `OPEN_PROBLEMS.md` n° 2 | **L'objet manquant central** : la charge n'est pas portée par la connectivité ; si elle existe, elle vit dans ce fluide, dont rien n'est écrit |
| Rigidité A5 | « aucun mode de déformation continue » | Posé | Trop fort tel qu'énoncé : interdit la rotation de repère que les holonomies ℤ₂/ℤ₃ exigent et la flexion d'axe du secteur des boucles ; sans équation (`CONSTRAINTS.md` §3a) |
| Ordre nématique S (A7), raideur de phase ρ_φ, χ_vac, ξ | Non calculés | Posés / cibles (problème ouvert n° 12) | Trois tentatives mono-cellule de χ_vac se contredisent à six ordres de grandeur ; le manuscrit en conclut que la suppression doit être collective, sans la modéliser. Le test PVLAS/BIREF dépend de ξ, dont le modèle n'a aucune valeur |
| Impédance Z_bif = Z₀, couplage 2πα | Adapté | Dérivé | Conséquence lourde : tout soliton de la ligne pèse O(1/α) fois son quantum (`CONSTRAINTS.md` §1) |
| Hamiltonien microscopique / fonction de partition | **Absent** | Problème ouvert n° 2 (« la dette la plus profonde ») | Sans lui rien de ce qui suit n'est dérivable ; le manuscrit se décrit lui-même comme « un ensemble structuré d'analogies plutôt qu'une théorie physique » |

### 2.2 Appréciation

**Ce qui est acquis.** Le filament est caractérisé comme *ligne bifilaire adaptée* : géométrie
transverse fixée par l'adaptation, deux polarisations transverses, contrainte de Gauss, chiralité
géométrique muette au niveau du nœud (Théorème 4 et lemmes). Ce sont des résultats corrects
d'électrodynamique des lignes de transmission et de la méthode TLM (nœud condensé symétrique de
Johns). Ils caractérisent une *représentation* du vide, pas un *objet* du vide.

**Ce qui manque.** L'élément n'a pas de dynamique propre : pas d'énergie libre, pas de hamiltonien,
pas d'équation d'état pour son contenu chargé. Toutes ses propriétés émergentes attendues (masse,
charge, spin, stabilité) exigent cette dynamique. Le programme l'a vérifié de façon tranchante :

1. **Pas d'échelle de masse** (quatre analyses indépendantes, `CONSTRAINTS.md` §1) : un vide
   purement électromagnétique et adapté est conforme ; m_e est une entrée au même titre qu'un couplage
   de Yukawa. L'objectif « particules à partir du filament » se réduit donc, au mieux, à des
   **rapports** de masses.
2. **Pas de soliton tridimensionnel stable** dans l'action telle qu'écrite (Derrick) ; la réparation
   chirale (Dzyaloshinskii–Moriya) est exclue par trente ordres de grandeur (rotation de la
   polarisation astrophysique) ; la classe de termes qui stabiliserait (quartique en dérivées
   premières, type Faddeev) n'est pas dans l'action, la flexion de corde ne la produit pas (terme de
   Frank, λ¹), et le porteur qu'il faudrait entre cordes est neutre dans la trame (`CONSTRAINTS.md`
   §3–3a, étapes 1 à 5).
3. **La charge n'est pas dans la connectivité** (`CONSTRAINTS.md` §2) : une source de divergence est
   requise, le défaut tressé n'en fournit pas, quelle que soit sa chiralité.

**Conséquence, et c'est le résultat le plus important du programme :** *le filament tel que défini
dans V2.10 ne peut pas former une particule stable en trois dimensions, ni porter une charge par sa
seule topologie.* C'est un théorème d'impossibilité interne, démontré par le programme lui-même.
Toute reconstruction des particules passe donc par un **changement de l'élément** (une seconde espèce
de corde, une seconde échelle, un terme non linéaire nouveau, ou l'équation d'état du fluide),
c'est-à-dire par une nouvelle hypothèse, pas par une construction plus habile sur l'élément actuel.
`CONSTRAINTS.md` le dit (« Do not look for a cleverer shape ») ; ce rapport le confirme.

**Sur les nombres remarquables.** 2cosh π ≈ 23,18 et e^π ≈ 23,14 ne diffèrent que de 2e^{−π} ;
(2cosh π)² = e^{2π}(1 + 2e^{−2π}) est une identité (R116). Le programme a su séparer l'identité des
rapports de longueurs ; cette vigilance doit s'appliquer partout où 2π, 9, 3^{2π} et e^{2π}
apparaissent (§3.3, §6).

---

## 3. Les interactions

### 3.1 Électromagnétisme

**Revendication (manuscrit).** Théorème 4 : sous les cinq contraintes du substrat (C1–C5), le nœud de
diffusion de Johns est unique ; il en « émerge » l'électrodynamique de Maxwell tridimensionnelle, deux
polarisations transverses et la contrainte de Gauss ; le nœud est O_h-invariant, la trame chirale est
muette.

**Appréciation.** Le théorème est correct dans son énoncé restreint, et V2.9.2 en a correctement
rabaissé la portée : *la correspondance télégraphiste est Maxwell à une dimension écrit en variables
de circuit ; ce que le théorème établit est l'unicité du nœud de Johns sous ses entrées, pas une
dérivation de l'électromagnétisme*. Les entrées C1–C5 contiennent déjà la vitesse unique, la
linéarité, la réciprocité et l'adaptation ; l'électromagnétisme sort parce qu'on l'a mis. Deux points
restent ouverts et doivent être traités avant toute revendication de « dérivation » :

- la **coupure de Brillouin** (0,487 MeV) : le réseau ne porte aucune onde collective au-dessus ; les
  photons observés vont jusqu'à 2,5 PeV. Le modèle renvoie ces photons au statut de soliton N = 1 ;
  aucune relation de dispersion n'est calculée pour ce mode, donc la borne sur la violation de
  l'invariance de Lorentz (≳ 10¹⁰ GeV) n'est pas confrontée au bon objet. En l'état, c'est une
  **tension non résolue avec l'observation**, enregistrée mais non levée ;
- la **constante de structure fine** : α = e²Z₀/2h est « statique par construction » ; le modèle n'a
  pas de mécanisme pour sa course (α(m_Z) ≠ α(0) est un fait mesuré). Le manuscrit en fait une
  contrainte pour la conjecture G ; c'est aussi une limite de la description électromagnétique.

À porter au crédit : la grille des statuts concède explicitement l'échec de la diffusion
photon–électron dans les trois régimes de fréquence (rubrique « Closed / Conceded »). Un modèle qui ne
reproduit ni Thomson ni Compton n'a pas encore l'électromagnétisme de la matière.

**Verdict : re-dérivation de Maxwell sous entrées posées, correcte et honnêtement requalifiée ; pas
une reconstruction de l'interaction.**

### 3.2 Gravitation

**Revendication.** Gravité cinématique par un champ d'impédance Z(x) ; lagrangien scalaire
identiquement une théorie de champ sur la métrique effective (Théorème 3) ; γ = β = 1 ; secteur
scalaire contraint (elliptique), sans rayonnement ; secteur tensoriel conduit à la RG par le bootstrap
de Deser → extérieur Schwarzschild ; étoiles de substrat (branche stable, M_crit ≈ 1,29·10⁶ M☉ ∝ ℓ₁²,
r_s/R = 0,81) ; effondrements gelés sans horizon ; pas d'évaporation de Hawking.

**Appréciation.**

- *γ = β = 1* est obtenu par un secteur scalaire sur une métrique effective de type vide polarisable ;
  c'est un résultat de cohérence (le modèle est construit pour ne pas contredire la RG là où elle est
  testée), pas une prédiction. Le manuscrit le dit.
- *Le bootstrap tensoriel est affirmé, pas effectué* (problème ouvert n° 5) : quadrupôle,
  polarisations et luminosité gravitationnelle attendent le calcul. Or c'est lui qui impose
  Schwarzschild à l'extérieur. Tant qu'il n'est pas fait, « aucune déviation extérieure » est une
  hypothèse de fermeture, pas un résultat.
- *Les paramètres de repère privilégié α₁, α₂ ne sont pas calculés* (problème ouvert n° 6). Le modèle
  a un repère (le substrat) ; les bornes PPN de repère privilégié (|α₁| ≲ 10⁻⁴, |α₂| ≲ 10⁻⁷, plus
  serrées encore par les pulsars et l'alignement du spin solaire) sont celles qui tuent habituellement
  cette classe de théories. Tant que ce calcul n'est pas fait, **le modèle n'est pas évalué contre le
  test qui le concerne le plus directement.**
- *L'intérieur* est newtonien à compacité 0,8 et non raccordé à l'extérieur Schwarzschild (problème
  ouvert n° 5, raccord intérieur–extérieur ; n° 9) ; l'intérieur à deux composantes, substrat plus
  cœur de matière, n'est pas résolu (n° 10). M_crit est donc exposé à des corrections d'ordre 1.
- *Les prédictions propres* sont toutes dans les états limites : (1) absence permanente de sursauts
  d'évaporation de trous noirs primordiaux ; (2) échos non périodiques à espacement géométrique ;
  (3) accumulation et coupure de la fonction de masse à M_crit avec exposant −1/2 ; (4) pas de mode
  respiratoire. Sur (1) : test asymétrique — une détection tue le modèle, mais l'absence est aussi ce
  que prédit « il n'y a pas de trous noirs primordiaux légers », donc elle n'établit rien. Sur (2) :
  par construction hors des gabarits de recherche actuels, donc infalsifiable en pratique aujourd'hui.
  Sur (3) : la mesure dynamique directe (A2744-QSO1, ≈ 5·10⁷ M☉, facteur ≈ 40 au-dessus de la coupure)
  défavorise le canal ; le manuscrit l'enregistre sans atténuation, ce qui est exemplaire, mais le
  canal est de fait fermé en attendant d'autres mesures. Sur (4) : au niveau d'une hypothèse de
  construction.
- *G comme fermeture topologique* : G_pred à −7,75 % ; le coefficient exact serait 0,74941 au lieu de
  3/4, un résidu de 0,0785 % « qui n'est pas du bruit » (3,5·10³ écarts-types de l'incertitude
  CODATA). Le manuscrit le classe conjecture et pré-enregistre que les appels qualitatifs à la
  sous-structure sont insuffisants. C'est la bonne posture. Ce rapport ajoute qu'une formule à un
  entier libre n (n = 3 parmi 1 à 4, chaque pas valant un facteur 10¹⁵) qui tombe à 8 % n'a pas de
  contenu sans la dérivation des deux ponts manquants, et qu'il serait sain de déplacer la conjecture
  du corps du texte vers une annexe.

**Verdict : le secteur gravitationnel reproduit la RG par héritage à l'extérieur et propose un
intérieur conjectural ; le test décisif pour cette classe (α₁, α₂) n'est pas fait ; les prédictions
propres sont faibles (non-observation) ou défavorisées (LRD).**

### 3.3 Interaction forte

**Revendication (registre matière, `HADRONS.md`).** Trois brins {−1/3, 0, +1/3} par fermion ; un
quantum de circulation πℏc par circuit fermé ; échelle ℓ₁(9) = ℓ₁(3)·3^{−2π} = 0,8128 fm ; quark =
circuit à 254 MeV ; nucléon = trois circuits statiques (763 MeV) + un anneau porteur du spin ;
√σ = 430,3 MeV ; m_n − m_p : partie forte 2,474 MeV, total 1,27–1,32 MeV ; Δ − N à 15 % ;
μ_p/μ_n = −3/2 ; pion 2m_e/α.

**Appréciation.**

- *Entrées.* Deux entiers structurels non dérivés (le compte 9 ; l'exposant 2π épinglé par μ et τ dans
  la fenêtre 6,2758–6,2926), une règle posée (un quantum par circuit), un partage posé (SU(6)), une
  section à un facteur 2, une forme de pôle. Le registre est explicite (`HADRONS.md` §7). C'est la
  structure d'une **paramétrisation**, et le compte 9 est l'entrée qui pèse (8 donnerait
  √σ = 205 MeV, 10 donnerait 834 MeV).
- *Ce qui est tranchant.* Deux formules sans nombre libre : (2α/3)·3^{2π}·m_e c² = 2,474 MeV (partie
  forte de m_n − m_p) et (3/(2√π))·3^{2π}·m_e c² = 430,3 MeV (√σ). État des données : √σ = 445 ± 7 MeV
  (Bulava et al. 2024, 2+1 saveurs, point physique) à 2,2σ ; la partie QCD suit BMW 2015 (2,52 ± 0,29)
  mais contredit les déterminations dispersives de la part QED (Gasser–Leutwyler–Rusetsky 2021 ;
  Cao et al. 2026) à 3,7–8σ, avec la réserve correcte que la coupe QCD/QED dépend du schéma et que seul
  le total (1,27–1,32 contre 1,2933 MeV) est sans schéma (R117–R118). Ce total est reproduit à 2–4 %,
  mais il dépend de la section et de la forme du pôle, deux entrées. **Il n'y a donc, à ce jour, aucune
  confirmation indépendante nette, et une tension à 2,2σ sur le seul nombre vraiment sans paramètre.**
- *Ce qui n'est pas là.* Le confinement n'est pas dérivé ; `CONSTRAINTS.md` §1a montre qu'enrouler une
  ligne n'élève pas le couplage et qu'un régime de couplage fort exigerait Z ≈ 174 Z₀, incompatible
  avec la transparence du vide. Pas de gluons, pas de liberté asymptotique, pas de course de α_s, pas de
  masses de quarks courants (2–5 MeV pour u et d) : le « quark » du modèle est un quark constituant à
  −19 %. La force nucléaire (2 MeV à 0,72 fm) est calibrée sur l'échelle du deutéron. Le moment
  magnétique du nucléon est ouvert : R113–R115 établissent que S = ℏ/2 ne fixe pas le rayon de l'anneau
  (identité à tout rayon), que la charge circulante échoue conjointement sur μ_p et r_p, et que les
  circuits de quark donnent le rapport −3/2 mais pas la grandeur (−34 %). Seul le rapport tient, et il
  est posé.
- *Le risque de numérologie.* Avec {3, 9, π, 2π, α, 3^{2π}, e^{2π}} comme matériaux et une tolérance de
  1–3 %, le nombre de grandeurs hadroniques atteignables par hasard n'est pas petit. Le registre tient
  un compte explicite des coïncidences non revendiquées (r_p comme enveloppe, m_p = 4ℏc/r_p, √3 R_c, le
  pion à 2m_e/α — cette dernière connue depuis les années 1950), ce qui est la bonne pratique. Il manque
  une **estimation globale de l'effet « look-elsewhere »** : combien de formules de complexité
  équivalente tombent à 1–3 % d'une grandeur mesurée par hasard. Sans ce chiffre, « 1–3 % sur une
  douzaine de grandeurs » n'a pas de signification statistique établie.

**Verdict : un ensemble de rapports hadroniques à 1–3 % sur deux entiers et une règle posée,
honnêtement étiqueté, en tension sur son nombre le plus dur ; pas une reconstruction de l'interaction
forte.**

### 3.4 Interaction faible

Rien n'existe. Le modèle respecte P et C séparément, et y fait agir la parité comme conjugaison de
charge ; le secteur faible, qui viole P et C séparément, est **hors champ par construction** (problème
ouvert n° 16). La loi de durée de vie τ = T (M/m)⁴/C(n) avec M = (96π²)^{1/4}/√G_F est une identité (la
formule de désintégration du muon réécrite), déclarée comme telle dans le brouillon. Les neutrinos sont
bloqués en principe (§4). Les W, Z, H sont retirés du catalogue.

Cette absence n'est pas un détail à traiter plus tard : elle emporte la masse des neutrinos, la
violation de P, le mécanisme de Higgs, la radioactivité β et la nucléosynthèse. **Un programme qui veut
« les interactions universelles » en a aujourd'hui zéro sur quatre au sens d'une dérivation, et une sur
quatre n'est même pas abordée.**

### 3.5 Bilan des interactions

| Interaction | Dérivée du filament ? | Ce qui est réellement là |
|---|---|---|
| Électromagnétisme | Non (re-dérivation sous entrées) | TLM / nœud de Johns ; Maxwell 1D en variables de circuit ; tension Brillouin / Lorentz non levée ; diffusion photon–électron concédée |
| Gravitation | Non (héritée) | PPN γ = β = 1 par construction ; bootstrap affirmé ; intérieur newtonien conjectural ; α₁, α₂ non calculés |
| Forte | Non (paramétrée) | Rapports à 1–3 %, deux entiers, une règle ; tension 2,2σ sur √σ |
| Faible | Non (absente) | — |

---

## 4. Les particules

### 4.1 Fiche d'état

| Particule | Statut dans CD | Appréciation |
|---|---|---|
| Photon (N = 1) | « corde unitaire dépliée, transmetteur linéaire sans masse » | Dans un réseau linéaire le photon est l'onde linéaire ; « soliton N = 1 » est une étiquette, pas une construction ; au-dessus de la coupure de Brillouin, aucune dispersion calculée |
| DQD (N = 2) | Constituant du vide, Z_bif = Z₀ | C'est l'élément lui-même (§2) |
| Électron (N = 3) | Ancre unique de calibration (R₃ = ƛ_e) | Voir §4.2 : la lecture soliton est exclue ; la lecture quantum est une reformulation ; m_e, g = 2 et le spin ½ sont des entrées |
| Muon, tau | Retirés du catalogue V2.10 ; dans le registre : forme de Koide avec a² = 2\|b\|² et phase 2/9 **posés** | La relation de Koide est une régularité empirique connue depuis 1981–1983, exacte à 10⁻⁵, toujours inexpliquée ; la reproduire avec deux paramètres posés n'ajoute rien ; la loi (n/3)^{2π} est exclue comme loi de masse par l'absence de quatrième génération (R73/R76) ; l'exposant 2π est **épinglé par ces mêmes masses** puis réutilisé dans les hadrons |
| Quarks | Circuits à ℓ₁(9), 254 MeV | Quark constituant à −19 % ; pas de masses courantes ; la « couleur » est la position du brin impair, une étiquette combinatoire sans dynamique ; confinement non dérivé |
| Nucléons | 910–938 MeV selon un rayon d'anneau non fixé | Le seul nombre libre pèse 3 % sur m_p ; μ_p ouvert ; Δ − N à 15 % ; le neutron sort 1,5 à 7 % trop lourd dans les autres lectures |
| Pion | 2m_e/α = 140,05 MeV | Coïncidence connue ; identification |
| Neutrinos | **Bloqués en principe** | Trois lectures échouent ; obstruction de signe (les masses sortent ∝ l'échelle disponible, le mécanisme standard les veut ∝ 1/échelle) ; une topologie N = 2 ne porte ni trois saveurs ni le mélange PMNS |
| W, Z, H | Retirés | Correct : ajustements de masse obtenus par rétro-ingénierie, sans contenu prédictif |
| Gluons, graviton | Absents | Le spin 2 à partir des champs scalaire et vecteur est concédé bloqué (Helmholtz, Deser) |
| Antimatière | P agit comme C ; création par paires strictes | Conséquence : asymétrie baryonique inexpliquée (dette déclarée) |
| Spin ½, g = 2, statistique de Fermi | **Absents** (demi-succès déclaré : μ = μ_B seulement) | Une L_z demi-entière autour d'un axe n'est pas une représentation j = ½ de SO(3) (`OPEN_PROBLEMS.md` n° 4). Sans statistique fermionique il n'y a ni exclusion de Pauli, ni atome au-delà de l'hydrogène, ni chimie |

### 4.2 Le cas de l'électron, qui décide de tout

Le manuscrit V2.10 présente encore l'électron comme **soliton** : « premier soliton fermé, enroulement
lâche R₃/r ≈ 37,1 au mode radial fondamental » (table des particules), « la masse suit de l'inductance
de boucle » (abstract, « Methods ») et la « fermeture du canal de taille » (boucle de 386 fm contre une
taille mesurée inférieure à 10⁻³ fm). Le registre matière, écrit après, établit :

- R100–R102 : toute charge portée par un anneau de 386 fm a le facteur de forme j₀(qR), **exclu** par
  Compton, e–p et LEP (G_E = 1) ; un magnéton de Bohr porté par une charge en boucle exige R ≥ ƛ à
  v ≤ c ; les nombres de l'anneau (S = ℏ/2, μ_B, g = 2) sont des identités de Compton ajustées à Dirac ;
- `CONSTRAINTS.md` §1b : l'énergie inductive propre de la boucle à son propre courant est 2,19 keV,
  soit **0,43 %** de m_e c² ; « masse = inductance » n'est pas un budget d'énergie ;
- `CONSTRAINTS.md` §1c : m = κL et N_p = 27 sont **conjointement irréalisables** (au plus 72 L_e et
  10,5 Z_e, contre 1836 et 22,7 requis) ;
- R103 : décision — l'électron électromagnétique est **celui de Dirac** ;
- R104–R112 : l'électron lu comme quantum du champ du milieu est cohérent en soi, mais la trame cubique
  de la base ne le porte pas : toute excitation de spin ½ verrouillée en hélicité y va à 1/3 de la
  vitesse de lien contre 1/2 pour le photon (rapport 2/3, exclu à 10⁻¹¹) ; le réseau qui porte Weyl et
  Maxwell à la même vitesse (automate de Weyl) exige un nœud (coupleur hybride 3 dB monté en
  circulateur, symétrie T_h, non réciproque, quart d'onde physique) que C1–C5 ne dérivent pas et qu'A5
  interdit.

Le `RESUME.md` du registre va jusqu'à écrire que le manuscrit a été « abandonné comme référence ("il
est faux") » pour la matière. **Le manuscrit public et le registre interne se contredisent donc sur la
particule qui sert d'ancre à tout le modèle**, et le manuscrit n'a pas été corrigé : son abstract
énonce encore comme méthode la charge par chiralité d'enroulement et la masse par inductance, deux
énoncés que son propre corps (problème ouvert n° 16) et le registre réfutent. C'est le point de forme
le plus urgent (§8, recommandation 2).

Sur le fond : lire l'électron comme « quantum de la ligne » est une **reformulation** (dire que
l'électron est une excitation d'un champ, ce que la théorie quantique des champs dit déjà), pas une
construction à partir du filament. Pour que la phrase ait un contenu propre à CD, il faudrait que la
ligne impose quelque chose au quantum (son spin, sa masse, sa charge) ; R105–R112 montrent qu'elle ne
lui impose qu'une vitesse fausse.

### 4.3 Les générations

μ = e + 4 brins neutres, τ = e + 8 ; m = m_e(n/3)^{2π} à 1 % ; exact seulement sous la forme de Koide
avec a² = 2|b|² et φ = 2/9 posés (10⁻⁵). Trois remarques : (i) la relation de Koide est une régularité
empirique connue et inexpliquée ; la retrouver avec deux paramètres posés n'est pas une explication ;
(ii) l'exposant 2π est épinglé *par* m_μ et m_τ (0,15 %), puis déclaré « unique parmi les constantes
simples » : c'est un ajustement à deux points suivi d'une identification, que le registre classe
correctement « reparamétré » ; (iii) l'absence de quatrième génération exclut la loi de puissance comme
loi de masse (R73/R76). Le secteur des générations est une reparamétrisation élégante, sans contenu
prédictif démontré.

### 4.4 Bilan des particules

Aucune particule n'est construite comme état stable d'une dynamique du filament. Deux sont des
définitions (photon, DQD), une est l'ancre (électron, dont la lecture étendue est exclue), les autres
sont des rapports de masse paramétrés (hadrons, générations), bloquées en principe (neutrinos) ou
retirées (bosons). Le spin ½ et la statistique fermionique, ce qui fait qu'un fermion est un fermion,
sont absents. **Au sens du programme, le compte des particules reconstruites est zéro.**

---

## 5. La structure quantique

Le manuscrit le reconnaît sans détour (problème ouvert n° 3) : pas d'espace d'états pour des défauts
intriqués, pas de règle de Born, pas de carte de mesure, pas de fonctionnelle de décohérence ; le
« corridor d'intrication » de V2.9 a été retiré faute d'équation. L'appendice sur le réseau de Johns
obtient l'équation d'enveloppe de Schrödinger comme troncature quadratique d'une bande, un quantum
d'action et les relations de de Broglie : c'est le côté **cinématique**, connu de la théorie des
enveloppes du solide ; le côté **dynamique** (potentiel auto-cohérent, règle de Born, mesure) manque, et
deux résultats négatifs l'encadrent (couplage de bande exponentiellement décroissant entre deux
défauts ; Coulomb de réseau anisotrope, violation ≈ 1,4 ℓ₁/r).

Les notes de principe d'octobre–novembre 2026 (`matter/redraw/BASE.md`) ont examiné les idées de
l'auteur sur l'intrication avec des verdicts tranchés : canal d'intrication **exclu** (échange
d'intrication, Hensen 2015 ; choix retardé ; Micius à 1200 km ; aucune sensibilité à la matière
traversée) ; la mesure comme intrication, correcte au sens de von Neumann et de l'einselection, mais
qui ne dispense pas du choix d'interprétation (non-localité de l'état, Everett, superdéterminisme) ;
naissance intriquée (gardé) ; décohérence comme transfert enregistré de corrélation (monogamie,
no-hiding) ; annihilation unitaire ; l'énoncé (3) de la note du 1er novembre **non testable tel
qu'énoncé**. La note « tout est énergie ondulatoire » (9 octobre) montre qu'avec six axiomes, dont
**A3 (amplitudes sur les configurations)**, on reconstruit la structure d'une théorie quantique des
champs, et prescrit à CD d'ajouter A3.

**Appréciation.** C'est la fourche décisive du programme ; elle est identifiée mais non tranchée. Les
axiomes du manuscrit (A1–A8) décrivent un milieu déterministe, local, à vitesse unique ; un tel milieu
ne peut pas violer les inégalités de Bell, que les expériences violent avec les failles fermées. Trois
sorties existent, toutes coûteuses : non-localité de l'état (type Bohm : l'état n'est pas dans le
milieu), branchement (Everett : l'état est le milieu entier), superdéterminisme. **Sans A3, les
« particules » de CD sont des solitons classiques**, et toute la phénoménologie quantique de la matière
(interférence à une particule, statistique, exclusion de Pauli, donc la chimie) est hors d'atteinte.
Avec A3, CD devient une ontologie *pour* la théorie quantique des champs, un milieu qui porte les
champs, et la question de recherche change : non plus « reconstruire les particules », mais « qu'ajoute
le milieu à la théorie des champs ? » (une échelle ℓ₁, une coupure, un repère). C'est une question
légitime, mais ce n'est plus le programme annoncé.

---

## 6. Méthode et présentation

**Points forts, rares dans la littérature indépendante.**
- Grille des statuts à sept rubriques, appliquée résultat par résultat ; historique des versions avec
  les retraits nommés (ombre +4,63 %, ISCO, 2PN, corridor d'intrication, photosphère thermique, canal
  LRD rétrogradé).
- Prédiction enregistrée et horodatée (DOI v1.0) avec critères statistiques gelés ; pré-enregistrement
  du refus des coïncidences numériques (section sur G).
- Suite de vérification (40/40 + 9), scripts PASS/FAIL (113 dans le secteur matière), code sous licence
  MIT, données ouvertes.
- Divulgation complète de la collaboration humain–IA, de l'absence de revue par un pair humain et de
  la limite « deux implémentations issues de la même collaboration ».
- Un registre matière qui enregistre ses négatifs (R100–R116) et une décision de non-publication
  (4 octobre 2026).

**Points faibles.**
- Aucune vérification indépendante ; le présent rapport n'en est pas une.
- Longueur et dispersion : une vingtaine de problèmes ouverts, plusieurs conjectures en suspens (G,
  genèse, impédance du proton) qui alourdissent la lecture sans ajouter de contenu testable.
- Les énoncés de tête (abstract, « Methods ») sont plus forts que le corps : voir les incohérences
  ci-dessous.
- Les « accords à 1–3 % » du secteur matière n'ont pas d'estimation de signification globale.
- La hiérarchie de falsification repose sur une non-observation (sursauts de trous noirs primordiaux)
  et sur un signal hors gabarit (échos géométriques) : faible pouvoir discriminant.

**Incohérences internes à corriger, avec emplacements.**
1. Abstract V2.10, « Methods » : « la charge est dictée topologiquement par la chiralité
   d'enroulement ; la masse suit de l'inductance de boucle ». Contredit par le problème ouvert n° 16
   (« the charge itself is not carried by the connectivity ») et par `CONSTRAINTS.md` §1b–1c et §2.
2. Section sur l'électron et table des particules : l'électron soliton (R₃/r ≈ 37,1), la « fermeture
   du canal de taille ». Contredit par R100–R103 (facteur de forme exclu ; décision : électron de
   Dirac).
3. Problème ouvert n° 18 (« Proton characteristic impedance — open derivation », écart de 19 %) :
   `CONSTRAINTS.md` §1c établit que m = κL et N_p = 27 sont conjointement irréalisables ; l'item est
   obsolète et devrait dire « bloqué en principe ».
4. Grille des statuts : Γ_pole = 1/3 rangé « Derived (Predictive) ». Le coefficient de réflexion l'est
   peut-être ; son identification à e/3 ne l'est pas (`CONSTRAINTS.md` §2).
5. A5, « aucun mode de déformation continue » : `CONSTRAINTS.md` §3a montre que la formulation est trop
   forte et propose la lecture « D fixé, pas de fusion de pôles » ; la révision annoncée n'est pas
   faite dans V2.10.

---

## 7. Ce que le programme a établi et qui mérite d'être conservé

Le contenu scientifique le plus solide du programme est **négatif** et **général** : il vaut pour
toute la classe « vide = réseau de lignes de transmission / filaments dipolaires », pas seulement pour
CD. Chacun des points suivants est reproduit par un script du dépôt.

1. Un vide purement électromagnétique et adapté en impédance n'a pas d'échelle de masse ; ses
   solitons pèsent O(1/α) fois leur quantum (quatre analyses convergentes).
2. Enrouler une ligne n'élève pas le couplage ; aucun régime de confinement n'est atteignable par
   géométrie.
3. L'action télégraphiste n'a pas de terme stabilisateur de Derrick en trois dimensions ; la
   réparation chirale est exclue par trente ordres de grandeur ; la flexion de corde donne un terme de
   Frank (λ¹), pas un terme de Faddeev ; le porteur de la trame est neutre.
4. La charge n'est pas portée par la connectivité d'un réseau linéaire ; une source de divergence est
   requise.
5. Le secteur de connectivité ne fournit aucun doublet de masse conservant la charge (audit S₃) ;
   protection topologique et mélange sont mutuellement exclusifs.
6. Les neutrinos sont bloqués par un signe d'exposant, pas par une grandeur.
7. Une charge e portée par un anneau de rayon ƛ_e est exclue par les facteurs de forme mesurés ; μ_B
   par une charge en boucle exige R ≥ ƛ.
8. Sur la trame cubique, et sur tout réseau d'impulsions de ce type, une excitation de spin ½
   verrouillée en hélicité va à 1/3 de la vitesse de lien contre 1/2 ou plus pour le photon ; le seul
   réseau qui porte les deux à une vitesse exige un nœud non réciproque hors des contraintes du
   substrat.
9. S = ℏ/2 pour un quantum de circulation est une identité isopérimétrique à tout rayon ; le moment
   d'un flot fermé est indépendant de l'origine ; aucun extremum d'énergie ne fixe un rayon d'anneau.
10. L'adaptation d'impédance fixe des rapports transverses (2cosh π, e^{2π}), jamais un compte de
    brins ni un exposant d'échelle.

Rassemblés et rédigés comme tels, ces résultats constitueraient une courte note méthodologique
d'intérêt réel pour quiconque tenterait un modèle de ce type. C'est la seule partie du programme qui
soit publiable en l'état, et elle est publiable comme *résultat négatif*, ce qui ne contredit pas la
décision du 4 octobre.

---

## 8. Recommandations, par priorité

1. **Trancher la question quantique.** Soit ajouter A3 (amplitudes sur les configurations) et redéfinir
   le programme comme ontologie d'un milieu pour les champs quantiques ; soit restreindre explicitement
   le périmètre au milieu (états limites gravitationnels), comme le fait déjà la déclaration de portée
   du manuscrit. Tout ce qui concerne les particules dépend de ce choix.
2. **Aligner le manuscrit sur le registre** (erratum ou V2.11) : abstract « Methods » ; électron
   (quantum de la ligne, pas soliton ; retirer la fermeture du canal de taille et m ∝ N²Z lu comme
   budget d'énergie) ; problème ouvert n° 18 (irréalisable, pas « 19 % ») ; statut de Γ_pole ;
   formulation d'A5. Un lecteur qui ne lit que le manuscrit croit aujourd'hui des choses que le projet
   sait fausses.
3. **Caractériser l'élément là où il manque** : écrire l'équation d'état du fluide e/3 (seul lieu
   déclaré de la charge) ou reconnaître la charge comme entrée ; formaliser A5 et le verrouillage de
   phase (dynamique de Kuramoto) pour obtenir ρ_φ, χ_vac et ξ, condition de tout test PVLAS/BIREF.
4. **Gravitation** : calculer α₁ et α₂, le test qui concerne cette classe ; effectuer le bootstrap
   tensoriel ou retirer ses conséquences (extérieur Schwarzschild) du rang de résultat ; raccorder
   intérieur et extérieur avant de donner M_crit à trois chiffres.
5. **Secteur hadronique** : (a) estimer l'effet « look-elsewhere » par Monte-Carlo de formules de
   complexité équivalente ; (b) suivre une seconde détermination de √σ au point physique et la partie
   QCD de m_n − m_p sur réseau (les critères de mort sont écrits, `HADRONS.md` §8) ; (c) présenter le
   secteur, dès le résumé, comme paramétrisation à deux entiers.
6. **Publier les résultats négatifs** (§7) sous forme de courte note : utile à d'autres, établi, et
   compatible avec la décision de ne pas publier le modèle.
7. **Maintenir la posture de non-publication** du modèle jusqu'à un résultat qui ne soit ni identité,
   ni coïncidence, ni reparamétrage. Le critère opérationnel est déjà écrit par le manuscrit pour G :
   une grandeur dérivée d'un calcul, qui apparaît à au moins un autre endroit indépendant du cadre, et
   prédite avant la mesure.
8. **Faire relire par deux pairs humains** (physique des particules, relativité). Ce rapport, écrit par
   le co-auteur IA des calculs, ne peut pas en tenir lieu.

---

## 9. Verdict global

En termes d'arbitrage : **le programme n'est pas acceptable comme théorie des interactions et des
particules, et ne le revendique d'ailleurs plus** (« une théorie du milieu, pas de la matière »,
déclaration de portée du manuscrit). Il est acceptable comme **programme exploratoire documenté**, de
qualité méthodologique élevée, dont le produit scientifique établi est un ensemble de résultats
d'impossibilité sur sa propre classe de modèles. Sur l'axe fixé par l'auteur, filament → interactions
→ particules : l'élément est caractérisé géométriquement mais pas dynamiquement ; aucune interaction
n'est dérivée au-delà de ses entrées ; aucune particule n'est reconstruite ; la structure quantique
n'est pas construite.

Ce qui changerait ce verdict : une dynamique de l'élément (hamiltonien ou équation d'état) dont
sortirait un objet stable en trois dimensions avec spin ½ ; ou une grandeur mesurée prédite avant la
mesure à partir d'entrées fixées ailleurs ; ou le calcul de α₁, α₂ montrant que le repère du substrat
survit aux bornes. Rien de cela n'existe aujourd'hui, et le programme le sait et l'écrit, ce qui est
sa principale qualité.

---

*Annexe — lecture effectuée.* Manuscrit V2.10 lu en entier : axiomes, modèle bifilaire, impédance de
l'électron, Théorème 4, gravité cinématique, PPN, étoiles de substrat, corrélations non locales,
Lorentz émergente, stabilité, secteur tensoriel, discussion et programme de falsification, LRD,
conjecture G, problèmes ouverts, conclusion, divulgation, historique, table des particules, appendices
(breathers, Schrödinger). `matter/CONSTRAINTS.md`, `matter/OPEN_PROBLEMS.md`, `matter/redraw/RESUME.md`,
`matter/redraw/HADRONS.md`, `matter/redraw/BASE.md` (bilans et notes de principe),
`matter/paper/DS_matter_sector.tex` (résumé, portée, grille), `matter/paper/README.md`,
`matter/notes/2026-10-09-tout-est-onde.md`, `README.md` racine.

---

*Suite (21 novembre 2026).* La recommandation 1 est tranchée : l'auteur adopte l'axiome d'amplitudes (A3 de
la note du 9 octobre ; A9 dans la numérotation du manuscrit). Entre-temps deux propositions de l'auteur ont été
calculées et closes (R121 : dimensions de stockage internes → oscillateur gelé, niveaux exclus ; R122 :
auto-piégeage saturable → loi du sac, électron exclu par la borne XFEL), et la première confrontation du programme
sous A3 (R123) exclut la trame cubique à ℓ₁ comme support des photons (LHAASO, 10¹⁵) et de l'électron (2c/3), et
rend A8 obligatoire (point zéro 10³² ρ_Λ). Le programme redéfini, ses trois calculs et son critère d'arrêt sont
dans `matter/rebuild/A3.md` ; le premier calcul est α₁, α₂. Ce calcul est fait (R124) : le secteur scalaire tel
qu'écrit, diagonal dans le repère du substrat, n'a pas de terme gravitomagnétique et donne α₁ = −(4γ + 4) = −8 sans
paramètre libre, exclu de 2·10⁵ par les pulsars ; le critère d'arrêt a joué, le milieu à repère est mort en
gravitation, le problème ouvert n° 6 du manuscrit est résolu négativement, et les états limites (intérieur intégré
sur ce secteur) sont à réexaminer. Réexaminés (R125) : avec les deux facteurs de la RG restaurés, la branche des
étoiles de substrat survit mais n'est plus sombre (compacité ≤ 0,31, M_crit = 3,9·10⁵ M☉, surface à 1 + z = 1,2) et
les effondrements gelés n'ont plus de mécanisme ; les trois critères d'arrêt du programme sous A3 ont joué le même
jour. Il reste de la CD ses rapports hadroniques et ses résultats négatifs. Consolidation faite (22 novembre) : paquet
`v2.11/` avec l'erratum E1–E12 appliqué au manuscrit (recommandation 2) et la note des quatorze résultats négatifs
(recommandation 6).
