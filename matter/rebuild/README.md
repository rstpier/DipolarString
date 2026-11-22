# Reconstruction de l'objet — `matter/rebuild/`

Décision de l'auteur (20 novembre 2026) : reconstruire l'élément (le filament dipolaire) **de zéro**, en
partant de son énergie, pas de sa forme. Ce dossier est un front de recherche, pas une version ; il ne modifie
ni le manuscrit V2.10 ni le registre `../redraw/` (qui garde la numérotation R des étapes).

- `SPEC.md` — le cahier des charges : les trois exigences que la revue du 19 novembre impose à tout nouvel
  objet, les bornes expérimentales qu'il doit respecter (R119 : u_s > 10³ u₀ ; R120 : la tension fixe un rayon),
  et pour chaque candidat le test qui le tue.
- `A3.md` — le programme sous l'axiome d'amplitudes (décision de l'auteur, 21 novembre 2026) : définition,
  première confrontation (R123), les trois calculs à faire dans l'ordre, et le critère d'arrêt écrit d'avance.
- `TDC.md` — la CD comme variante de la théorie des cordes (décision de l'auteur, 22 novembre 2026) : les trois sens
  de « variante », la définition retenue, le dictionnaire CD ↔ cordes, le premier audit (R127) et les critères d'arrêt.
- `veille_b3.py` — R135 : état de la littérature sur la question B3 (les deux tensions sont-elles la même ?) ;
  spectres 424 ± 2 MeV, réseau 445(3)(6) MeV, écart de deux à trois écarts types, aucune source ne tranche.
- `ds_lattice_tension.py` — R134, premier choix fait par délégation : la tension du réseau à la place de celle de la
  base ; exposant requis 6,314 ± 0,014 ; pentes de Regge trop faibles de 9 à 11 % ; aucun nombre propre à la CD ne
  reste dans la corde effective.
- `ds_effective_string.py` — R133, après le choix « ok prend la corde effective » : le calcul B1 de la feuille de
  route (corde tournante à extrémités massives, six familles de mésons, pente commune 0,883 ± 0,014 GeV⁻² contre
  0,884 publié) ; domaine de 0,33 à 1,25 fm, ρ et nucléon au bord ; tension de la base dans l'enveloppe des
  trajectoires (426 ± 6 MeV) et à deux écarts types sous le réseau (445 ± 7) : ni confirmée ni exclue.
- `ds_fine_link_limit.py` — R132, après le choix « ok prend le maillon infiniment fin » : la limite continue existe
  (écart en 1/M²) ; elle n'est invariante de Lorentz qu'à 26 dimensions, ou 10 avec des fermions (premier niveau sans
  masse) ; deux contre-termes ou la supersymétrie ; filament = ligne uniformément polarisée ; taille du fil
  α′(ln M + 0,82), sans valeur propre ; aucune taille d'élément.
- `ds_chain_of_links.py` — R131, après le choix « ok prend la chaîne de maillons » : la corde comme chaîne de maillons
  à ressorts (bits de corde) ; fréquences exactes et pente de Regge ; ordonnée (D − 2)/24 comme énergie de point zéro ;
  un nombre fini de maillons brise l'invariance de Lorentz (cos(π/2M)) ; maillon = dipôle, chaîne = produit de
  matrices, jonction = concaténation ; chaîne = échelle LC d'impédance égale à la tension.
- `ds_calabi_yau.py` — R130, après le choix « ok prend Calabi-Yau » : un quart des supersymétries ; les entiers de la
  quintique (1 taille, 101 formes, 204 cycles, χ = −200) retrouvés de trois façons sur dix variétés ; générations comme
  produit d'intersection dans le réseau des cycles (tout entier disponible) ; 52 à 410 champs sans masse à fixer ; les
  trois tores repliés (36, 0) ; χ = −6 absent d'un échantillon de 2362 configurations.
- `ds_three_groups.py` — R129, après les précisions « 9 dim : 3 groupes de 3 dimensions » et « tores, non probablement
  pas » : même espace que les trois tores si les groupes sont fermés et plats ; sinon un groupe vide est plat (à trois
  dimensions Riemann est une fonction de Ricci), une sphère tenue par un flux se paie dans l'espace visible (exclue),
  des groupes infinis donnent une gravité en 1/r⁸ ; reste le Calabi-Yau (χ = −200, −18, −6).
- `ds_ten_dimensions.py` — R128, après le choix de dix dimensions : spectre de la supercorde (8 + 8, 128 + 128), les
  trois tores, générations comme nombres d'intersection (3 non préféré), échelles, facteur de déformation ; deux
  rectificatifs de R127.
- `ds_as_string_variant.py` — R127 : identités (tout est ℓ_s), incohérence d'une corde fondamentale à D = 4, exclusion
  d'une corde unique hadrons + gravité, champ critique du dipôle, ligne de transmission comme boson compact.
- `limit_states_tov.py` — les états limites sous R124 (R125) : la branche des étoiles de substrat avec les deux
  facteurs de la RG restaurés (TOV, conventions A8) ; contrôle sur le système du manuscrit ; pas d'objet sombre.
- `alpha12_preferred_frame.py` — premier calcul du programme sous A3 (R124) : les paramètres de repère privilégié
  α₁, α₂, α₃ du secteur scalaire, substrat en mouvement (SymPy ; deux lectures PPN, contrôle RG) : α₁ = −8, exclu.
- `step1_saturable_selftrapping.py` — étape 1 (R122) : l'idée de base de l'auteur (« l'énergie monte Z là où elle
  est ; assez de quanta au même endroit se piègent ») écrite comme milieu de Kerr saturable et calculée : seuil
  N_min en 3D, échelle d'énergie E_s = (u_s (ℏc)³)^{1/4}, loi de masse M = N_min^{1/4} E_s N^{3/4} (le sac), bornes
  HIBEF et Schwinger. Verdict dans `SPEC.md` §4 et dans `../redraw/BASE.md` (R122).

Méthode inchangée : chaque étape est un script autonome avec vérifications PASS/FAIL et une ligne
`RESULT: N/N PASS; 0 FAIL`, enregistrée dans `../redraw/BASE.md` avec son verdict.
