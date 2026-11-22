# Changelog

## V2.11 (21 novembre 2026) — édition erratum, paquet gelé le 22 novembre

- Le manuscrit V2.10 gelé reçoit un erratum appliqué sur place, douze points marqués
  « Erratum En » dans le texte et listés dans l'historique du document ; chacun réduit
  les revendications, et deux (E6, E7) effectuent des calculs que V2.10 déclarait non
  faits, tous deux contre le modèle. Détail : `v2.11/ERRATUM.md`
- E1 : charge par chiralité d'enroulement et masse par inductance de boucle retirées du
  résumé (la charge n'est pas portée par la connectivité ; l'énergie inductive de la
  boucle vaut 0,43 % de m_e c²)
- E2 : l'électron comme soliton étendu est retiré (facteur de forme j₀(qR₃), magnéton,
  quantum LC) ; l'électron électromagnétique est celui de Dirac
- E3 : l'impédance du proton est bloquée en principe (m = κL et N_p = 27 irréalisables
  ensemble) ; E4 : Γ_pole = 1/3 est un coefficient de réflexion, son identification à
  e/3 une identification ; E5 : A5 reformulé (« D fixé, pas de fusion de pôles »)
- E6 : paramètres de repère privilégié calculés : α₁ = −(4γ + 4) = −8 sans paramètre
  libre, exclu de 2·10⁵ par les pulsars ; problème ouvert n° 6 clos négativement
- E7 : les étoiles de substrat réintégrées avec les deux facteurs de la RG : M_crit =
  3,9·10⁵ M☉, compacité maximale stable 0,31 ; aucun objet sombre, plus de mécanisme de
  gel ; falsificateurs 1–3 et prédiction LRD enregistrée retirés
- E8 : sous l'axiome d'amplitudes le photon porte la dispersion du réseau ; LHAASO exclut
  ℓ₁ comme pas ; E9 : un spineur de la trame cubique va à 2c₀/3 ; E10 : point zéro à ℓ₁ =
  10³² ρ_Λ, A8 obligatoire ; E11 : grille des statuts reclassée ; E12 : adoption de
  l'axiome d'amplitudes A9 par l'auteur (21 novembre 2026) enregistrée
- Note d'accompagnement `v2.11/negative_results/` : quatorze résultats négatifs structurels
  valables pour toute la classe « vide = réseau de lignes adaptées »
- `verify_v2_11.py` : contrôles de source, note, relance des sept scripts de preuve
  (R119–R125), suite V2.9.2 déléguée si ses dépendances sont présentes

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
