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
