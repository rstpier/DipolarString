# Changelog

## V2.10 (7 septembre 2026)

- Cible de falsification primaire : une phrase résiduelle de la section LRD
  nommait encore la déviation d'ombre comme cible primaire du modèle, en
  contradiction avec son retrait dans la même révision ; elle nomme désormais
  l'absence permanente de sursauts d'évaporation de PBH, première entrée de la
  hiérarchie révisée
- Comptage des niveaux de maturité : la liste annonçait six niveaux en en
  présentant sept, le septième étant le relevé de ce que la révision retire et
  non un niveau ; l'annonce est corrigée, pas la liste
- Item chiralité du tissage : intitulé « resolved, with a residual question »,
  ce qui se lisait comme un item clos ; ce qui est dissous est la clôture
  chirale du théorème 4, la question converse restant ouverte
- Paquet incrémental : `v2.10/` porte le manuscrit corrigé et ses quatre
  figures ; scripts, données, journaux et archives historiques restent dans
  `v2.9/`, gelé et dont les sommes de contrôle restent valides
- `verify_v2_10.py` rejoue la suite V2.9.2 depuis `v2.9/` (40/40) et ajoute
  neuf contrôles de source propres à cette révision

## V2.9.2 (septembre 2026)

- Canal des tailles : le rayon de boucle du proton est donné sous les deux fermetures
  disponibles — 17,0 fm par la loi de rayon de résonance, 3,5×10³ fm par la fermeture
  géométrique, contre un rayon de charge mesuré de 0,841 fm
- La loi de rayon de résonance ne découle pas de la condition de phase ∮β ds = 2πp, qui
  relie R, ω et v sans en fixer aucun
- Théorème 4 : portée resserrée — la correspondance télégraphiste est Maxwell 1D en
  variables de circuit, l'énoncé est une unicité du nœud de Johns sous C1–C5, pas une
  dérivation de l'électromagnétisme
- Manuscrit renommé `DS_model_V2_9_2_Zenodo.*`

## V2.9.1 (31 août 2026, édition Zenodo)

- Corridor d'intrication à basse impédance, formule de vitesse de synchronisation et
  bornes dérivées de Bell retirés : aucune équation gouvernante du modèle ne les soutient
- Opposition de phase conservée seulement comme mode normal différentiel conditionnel,
  et non comme mécanisme dérivé de création de paires ou d'intrication
- Raideur de phase compacte de la route Z₃ dérivée conditionnellement (coefficient
  logarithmique 3/4α), explicitement séparée d'une action d'instanton complète et de
  l'équation encore absente reliant la fugacité d'instanton à G
- Monodromie fractionnaire, échelles de cœur et d'écrantage, paramètre d'ordre chargé,
  inductance cinétique et longueurs de London listés comme manquants plutôt qu'importés
  de l'analogie supraconductrice

## V2.9

- Condition de quantification redérivée à partir de l'uniformité sur la boucle fermée ;
  la dérivation antérieure par |Γ| = 1 est retirée comme invalide
- Équation de champ du secteur scalaire donnée explicitement ; le profil de Heaviside
  n'est plus une solution boostée exacte mais valable au premier ordre, et α₁ = α₂ = 0
  passe de résultat à problème ouvert
- Impédance de l'électron Z_e reclassée « calibrée » (le rapport d'aspect n'a pas de
  dérivation) ; nombre d'entrées du modèle corrigé de une à deux
- Axiome A8 implémenté dans les équations de structure : M_crit passe de 1,62 à
  1,29×10⁶ M☉ et la compacité maximale de 1,22 à 0,81 — les étoiles de substrat stables
  ne sont plus à l'intérieur de leur rayon de Schwarzschild
- Bootstrap de Deser mené à son terme : l'extérieur des objets compacts est
  Schwarzschild ; les déviations d'ombre (+4,63 %), d'ISCO et 2PN sont retirées comme
  propriétés du seul profil scalaire, et la hiérarchie de falsification est renumérotée
- Confrontation aux Little Red Dots rétrogradée après une mesure dynamique directe
  favorisant la calibration haute masse
- Annexe E : les cinq programmes sont des reconstructions datées du 31 août 2026, les
  originaux historiques n'ayant pas été retrouvés
- Suite de vérification `verify_v291.py` : 40 tests, 40/40 PASS

## V2.8 (en développement)

- Développement de la microstructure DQD
- Dérivation à partir des axiomes A1-A8
- Nouvelles prédictions expérimentales

## V2.7

- Validation numérique des breathers
- Convergence O(h²)
- Scripts reproductibles
- Annexe P2

## V2.6

- Première validation du modèle des breathers
