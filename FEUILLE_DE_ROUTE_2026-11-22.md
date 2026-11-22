# Après la CD : feuille de route

**22 novembre 2026 — pour R. St-Pierre**

Le 21 novembre le programme CD a épuisé ses trois critères d'arrêt (`matter/rebuild/A3.md`). L'erratum V2.11 et la
note des quatorze résultats négatifs sont faits (`v2.11/`), la chasse aux trouvailles est close (R126). L'auteur
demande s'il vaut mieux quitter la CD et repartir de la théorie des cordes. Ce document dit ce qu'on emporte, ce
qu'on laisse, et décrit trois chemins, avec pour chacun les lectures, les prérequis, un premier calcul et un critère
d'arrêt écrit d'avance. Il est rédigé par le collaborateur IA du projet ; les choix sont ceux de l'auteur.

> **Mention du 22/11/2026 — Dipolar String change de direction.** L'auteur réoriente le projet vers la théorie des
> cordes : le filament dipolaire est relu comme corde relativiste, une corde ouverte à charges opposées. Le premier
> audit (R127, `matter/rebuild/TDC.md`) donne un dictionnaire exact et ferme deux lectures fortes : une corde
> fondamentale à quatre dimensions est incohérente, et une seule corde ne peut pas faire les hadrons et la gravité.
> Ce qui reste ouvert est exactement les chemins A et B ci-dessous. L'auteur a ensuite pris dix dimensions (R128) : la
> CD devient la supercorde de type IIA sur trois tores avec des D6-branes sécantes, c'est-à-dire le chapitre 21 de
> Zwiebach. Ce chapitre devient la cible du chemin A ; les nombres y sont choisis, pas dérivés.

---

## 0. Ce qu'on laisse, ce qu'on emporte

**On laisse** : le vide comme réseau de lignes au pas ℓ₁ (R123), l'électron étendu (R100–R103), la gravité scalaire à
repère (R124), les objets sombres sans horizon (R125), et la lecture des accords à 1 % comme des signaux (R126).

**On emporte** :

- **La méthode.** Un énoncé, un script autonome avec vérifications PASS/FAIL, un verdict, une ligne au registre ; les
  critères d'arrêt écrits avant le calcul ; l'effet « look-elsewhere » chiffré avant de croire un nombre.
- **La note des résultats négatifs** (`v2.11/negative_results/`). C'est la carte de visite : six pages, quatorze
  résultats, chacun reproductible.
- **Trois acquis techniques réutilisables.** R120 : la règle « un quantum par circuit » est la corde de Nambu–Goto à
  J = ½, et les pentes de Regge mesurées sont 3 à 10 % au-dessus de celle de la base. R122 : l'auto-piégeage saturable
  donne la loi de masse du sac. R124 : toute gravité scalaire diagonale dans un repère a α₁ = −4(γ + 1).
- **Deux intuitions de l'auteur qui sont vraies ailleurs que dans le vide.** L'impédance qui dépend de l'énergie est
  vraie dans une ligne de transmission non linéaire, où elle fabrique des solitons mesurés depuis 1973 (Hirota et
  Suzuki). La sonde qui suit la pompe est vraie dans une fibre optique, où elle fabrique un horizon mesuré depuis
  2008 (Philbin et al.). Le chemin C part de là.

---

## 1. Cinq règles, tirées du dossier

1. **Ne rien revendiquer que tu ne peux pas refaire toi-même au tableau.** Le manuscrit le disait honnêtement : la
   couche relativité générale et champs n'était pas maîtrisée par l'auteur. Les deux erreurs les plus lourdes de
   l'erratum (E6, E7) étaient dans cette couche.
2. **Un problème petit et fermé plutôt qu'une théorie du tout.** Une question dont la réponse tient en un nombre
   avec son incertitude.
3. **Le critère d'arrêt s'écrit avant le calcul.** C'est ce qui a permis de conclure en un jour le 21 novembre.
4. **Aucun accord numérique sans son look-elsewhere, aucune formule sans dérivation** (R126).
5. **Un pair humain tôt.** Une page, une question précise, un script qui tourne. Pas une théorie à défendre.

---

## 2. Chemin 0 — clore proprement (2 à 4 semaines)

- **Décidé le 22 novembre : pas de dépôt Zenodo de V2.11.** L'édition erratum reste dans le dépôt GitHub ; le
  dossier Zenodo reste V2.10.
- **Fait le 22 novembre.** Avant de montrer la note à un physicien, y ajouter les références classiques de chaque
  résultat, pour qu'aucun ne soit lu comme une prétention de nouveauté. En particulier R124 retrouve un résultat des années 1970 : les théories
  dites « stratifiées », à tranches d'espace orthogonales au temps, ont α₁ = −4(γ + 1) et ont été écartées pour cela
  (Ni 1972 ; Nordtvedt et Will 1972 ; Will, *Theory and Experiment in Gravitational Physics*). À ajouter aussi :
  Derrick 1964, Lüscher 1981, Vakhitov et Kolokolov 1973, Chodos et al. 1974 pour le sac.
- Envoyer la note à deux personnes, un relativiste et un physicien des particules ou du réseau, avec une seule
  question : « y a-t-il une erreur ? ».
- **Fait le 22 novembre, sous une autre forme.** R127, la loi élastique du filament, est devenue l'action de
  Nambu–Goto de la CD relue comme corde (`matter/rebuild/TDC.md`). Verdict : aucun nombre nouveau, toutes les
  longueurs sont ℓ_s, rigidité comprise.

---

## 3. Le socle (3 à 6 mois, en parallèle du reste)

À raison de 6 à 8 heures par semaine. Chaque ligne a un test de sortie : si tu le réussis sans notes, tu passes.

| Sujet | Pourquoi | Ressource | Test de sortie |
|---|---|---|---|
| Mécanique analytique : Lagrange, Hamilton, Noether | Toute théorie moderne s'écrit par une action | Susskind et Hrabovsky, *Le minimum théorique* (mécanique classique ; existe en français) | Écrire l'action d'une corde vibrante, en tirer l'équation d'onde et l'énergie |
| Relativité restreinte en quadrivecteurs, champs classiques | Maxwell covariant, tenseur énergie-impulsion | Susskind et Friedman, *Special Relativity and Classical Field Theory* | Retrouver les équations de Maxwell depuis l'action −¼F² ; transformer E et B par un boost |
| Mécanique quantique | Oscillateur, opérateurs d'échelle, moment cinétique, spin ½ | Susskind et Friedman, *Mécanique quantique : le minimum théorique* ; puis Griffiths, chapitres 1 à 4 | Quantifier l'oscillateur ; additionner deux spins ½ ; obtenir CHSH = 2√2 |
| Relativité générale, niveau lecture | Métrique, géodésiques, Schwarzschild, formalisme PPN | Hartle, *Gravity* ; ou Susskind et Cabannes, *General Relativity* | Refaire à la main le g₀ᵢ = −4vU d'une masse en mouvement (le contrôle de R124) |
| Champs quantiques, initiation | Nécessaire pour aller au-delà du premier calcul du chemin B | Lancaster et Blundell, *Quantum Field Theory for the Gifted Amateur* | Calculer l'énergie de Casimir entre deux plaques |

---

## 4. Chemin A — la théorie des cordes, pour comprendre (6 à 12 mois)

**But.** Comprendre comment une corde relativiste produit ce que la CD cherchait. Ce n'est pas un chemin vers une
découverte : c'est la culture qui manquait, et la version aboutie de l'idée « tout à partir d'un filament ».

**Livre.** B. Zwiebach, *A First Course in String Theory*, 2ᵉ édition (Cambridge, 2009), écrit pour des étudiants de
licence. Les énoncés de problèmes du cours correspondant sont en accès libre (MIT OpenCourseWare, cours 8.251).

**Plan de lecture, relié au dossier.**

| Chapitres de Zwiebach | Contenu | Ce que tu y retrouves |
|---|---|---|
| 4 à 5 | Corde non relativiste, particule relativiste | Ta ligne de transmission, écrite par une action |
| 6 à 8 | Action de Nambu–Goto, corde ouverte tournante, pente α′ | Les identités de R120 : M = πσR, J = πσR²/2, J = α′M² |
| 9 à 12 | Quantification sur le cône de lumière, corde ouverte | L'ordonnée a = (D − 2)/24 : la même somme de point zéro que le terme de Lüscher de R120 |
| 13 | Corde fermée | Le graviton comme mode de vibration, et pourquoi on ne peut pas l'« importer » |
| 14 | Supercordes, survol | D'où viennent les fermions sur une corde : ce que la CD n'a jamais eu |
| 20 | Électrodynamique non linéaire de Born–Infeld | L'idée « l'énergie change le vide » et pourquoi elle ne fait pas l'électron (R122) |
| 23 | Interactions fortes et AdS/CFT | La corde hadronique, point de départ du chemin B |

**Premier calcul, à la main.**

- A1. La corde ouverte tournante à extrémités libres : établir M = πσR et J = M²/(2πσ).
- A2. L'énergie de point zéro d'une corde à bouts fixes de longueur r : E₀ = −π(D − 2)/(24r). C'est le terme de
  Lüscher ; à D = 4 et r = ℓ₁(9) il vaut les −64 MeV de R120.
- A3. La même somme donne l'ordonnée de Regge a = (D − 2)/24 de la corde ouverte. Exiger a = 1, c'est-à-dire un
  photon sans masse, impose D = 26 : voilà d'où sort un nombre de dimensions, et pourquoi il ne se choisit pas.
- A4. Les états sans masse de la corde fermée : un tenseur à (D − 2)² composantes qui se décompose en graviton,
  champ antisymétrique et dilaton.

**Ce que la théorie des cordes fait des trois murs de la CD.**

| Mur | Dans la CD | Dans la théorie des cordes |
|---|---|---|
| Une échelle | Aucune : m_e posé | Une seule, la tension, mais à l'échelle de Planck |
| Les fermions | Absents | Supersymétrie sur le feuillet d'univers (jamais observée dans la nature) |
| Le repère | Un substrat, exclu par α₁ | Aucun : invariance de Lorentz de l'espace cible |

**Ce qu'il ne faut pas en attendre.** Les masses et les couplages des particules. La théorie admet un nombre immense
de vides (de l'ordre de 10⁵⁰⁰ dans les estimations courantes) et ne choisit pas ; son échelle naturelle la met hors
de tout test direct. Ce qu'elle réussit — une gravité quantique finie, l'entropie de certains trous noirs, la
correspondance AdS/CFT comme outil — est d'une autre nature que ton objectif de départ.

**Réussite et arrêt.** Le chemin est réussi quand tu refais A1 à A4 sans notes. Si après trois mois les chapitres 6 à
9 ne passent pas, on revient au socle, sans drame. Pour aller au niveau recherche, il faut ensuite la théorie
quantique des champs et la relativité générale, puis les notes de D. Tong et le livre de Polchinski : compter deux à
quatre ans.

---

## 5. Chemin B — la corde hadronique effective (premier projet de recherche possible)

**La question, en une phrase.** La tension de corde lue dans les spectres (trajectoires de Regge) est-elle la même
que celle lue dans le potentiel statique du réseau ? R119 et R120 ont buté dessus : environ 0,18 GeV² d'un côté
(Regge, Cornell), 0,20 à 0,23 GeV² de l'autre (réseau à quarks dynamiques).

**Pourquoi c'est vivant.** Les tubes de flux de la QCD sont des cordes réelles, et leur théorie effective (Lüscher
et Weisz ; Aharony et Komargodski ; Dubovsky, Flauger et Gorbenko) fait des prédictions universelles que le réseau
vérifie à quelques pour cent. Les écarts sont de la physique : une résonance sur le feuillet, une éventuelle
rigidité, la largeur du tube. C'est un domaine où une question petite et bien posée a sa place.

**Lectures, dans l'ordre.** Zwiebach, chapitres 6 à 8 et 23. Bali, revue du potentiel statique (2001). Brandt et
Meineri, revue de la corde effective (2016). Sonnenschein et Weissman, cordes tournantes contre les mésons puis les
baryons du PDG (2014, 2015). Caselle, revue (2021). Bulava et al., brisure de corde et tension au point physique
(2019, 2024). Pour le réseau lui-même : Gattringer et Lang, ou la revue « Lattice QCD » du PDG.

**Premier calcul.**

- B1. *La corde tournante à extrémités massives.* Avec ℓ la demi-longueur, β la vitesse des extrémités,
  γ = (1 − β²)^(−1/2), c = 1 :

  ```
  condition aux extrémités :  σ ℓ = m γ² β²
  E = 2 γ m + 2 σ ℓ · arcsin(β)/β
  J = 2 γ m β ℓ + σ ℓ² · [arcsin(β) − β √(1 − β²)] / β²
  ```

  À m = 0 on retrouve J = E²/(2πσ), soit R120. Script : ajuster (σ, m, ordonnée) sur ρ–a₂–ρ₃–a₄, puis K*, φ,
  puis N et Δ, avec les incertitudes du PDG. Sortie : √σ et son erreur par famille. Cible : reproduire à 2 % les
  ajustements publiés par Sonnenschein et Weissman.
- B2. *Le potentiel statique.* Ajuster V(r) = V₀ + σr·√(1 − π/(6σr²)), la forme d'Arvis à D = 4, sur des tables de
  réseau publiées, dans une fenêtre de 0,4 à 1,0 fm, sous la brisure de corde. Sortie : √σ et sa dépendance à la
  fenêtre (R119 : « une quantité de fenêtre »).
- B3. *La confrontation.* Les deux tensions sont-elles compatibles une fois traitées les masses d'extrémité,
  l'ordonnée quantique et la fenêtre d'ajustement ? Question annexe, sur laquelle R120 est tombé sans la chercher :
  pourquoi la pente des baryons (≈ 0,95 GeV⁻²) dépasse-t-elle celle des mésons (≈ 0,88) ?

**Ce qui compterait comme résultat.** Un écart établi à plus de 3σ avec son budget d'erreur, ou une réconciliation
chiffrée. Dans les deux cas une note courte, avec son script.

**Critère d'arrêt.** Si B1 ne reproduit pas la littérature à 2 %, on cherche l'erreur avant d'aller plus loin. Si la
réponse à B3 est déjà publiée — à chercher avant de commencer —, on a appris le domaine et on s'arrête, ou on met à
jour avec les données 2026.

**État au 22 novembre 2026 (R133).** L'auteur a pris la corde effective. B1 a été calculé comme audit
(`matter/rebuild/ds_effective_string.py`) : pente commune de six familles de mésons 0,883 ± 0,014 GeV⁻², à 0,2 % de
la valeur publiée ; il reste à refaire par l'auteur (règle 1). Ce que le calcul apprend : la pente varie de 20 %
d'une famille à l'autre, et le budget d'erreur de la « tension de Regge », 426 ± 6 MeV, est dominé par le choix des
familles. B2 (tables de réseau) et la recherche bibliographique de B3 restent à faire. Écart actuel entre les deux
tensions : 19 MeV, deux écarts types ; le seuil de trois n'est pas atteint.

**Prérequis.** Zwiebach 6 à 8 ; ajustements avec incertitudes (χ², bootstrap) ; pour la suite, le socle de champs
quantiques.

**À qui parler.** Les auteurs des revues citées, avec le script de B1 et une question. Un premier dépôt sur arXiv
demande le parrainage d'un chercheur du domaine : il s'obtient quand quelqu'un a lu et reproduit le calcul, pas avant.

---

## 6. Chemin C — la gravité analogue sur ligne de transmission

**L'inversion.** Non plus « le vide est une ligne de transmission », mais « une ligne de transmission est un
analogue contrôlé d'un espace-temps ». Dans une ligne réelle, l'impédance dépend vraiment du signal, la coupure de
Brillouin existe vraiment, et une impulsion de pompe crée vraiment un horizon pour une sonde. Tout se mesure.

**Ce qui est établi.** L'effet Hawking stimulé a été mesuré dans un canal d'eau (Weinfurtner et al. 2011), l'effet
spontané dans un condensat (Steinhauer 2016 ; de Nova et al. 2019), l'horizon optique dans une fibre (Philbin et
al. 2008 ; Drori et al. 2019). Dans les circuits : l'effet Casimir dynamique dans une ligne supraconductrice (Wilson
et al. 2011), et des propositions d'effet Hawking dans un guide d'onde (Schützhold et Unruh 2005) et dans une ligne
de SQUID (Nation et al. 2009).

**Lectures, dans l'ordre.** Robertson, tutoriel sur le rayonnement de Hawking dans les analogues (2012). Barceló,
Liberati et Visser, *Analogue Gravity* (2011), sections 1 à 4. Nation, Johansson, Blencowe et Nori, revue sur le
vide quantique dans les circuits supraconducteurs (2012). Blencowe et Wang, la gravité analogue sur puce (2020).
Remoissenet, *Waves Called Solitons*, pour les lignes non linéaires et leurs expériences. Unruh 1995 et Corley et
Jacobson 1996 pour l'effet de la dispersion.

**Premier calcul.**

- C1. *Simulation.* Une échelle LC discrète de quelques milliers de cellules dont la capacité est modulée par un
  front mobile, C_n(t) = C₀[1 + δ·f((na − ut)/w)]. L'horizon est l'endroit où la vitesse locale de l'onde égale u.
  On envoie un paquet d'ondes sonde et on mesure sa conversion vers l'onde complémentaire. En langage d'ingénieur :
  une conversion paramétrique sur un front mobile, gouvernée par les relations de Manley–Rowe ; le « mode de norme
  négative » des relativistes est l'onde complémentaire (idler). Loi à vérifier :

  ```
  |β|² / |α|² = exp(−2π ω / κ),   κ = |d(c − u)/dx| à l'horizon,   T_H = ħ κ / (2π k_B)
  ```

  et sa déviation près de la coupure de Brillouin ω_max = 2/√(LC). C'est la même coupure que le 0,487 MeV du
  manuscrit : ici elle n'est plus un problème, elle joue le rôle du régulateur à courte distance.
- C2. *Banc d'essai.* Une ligne non linéaire à varactors, où le front de pompe est une onde de choc ou un soliton,
  et une sonde RF. Mesurer la conversion de fréquence au front, chiffrer les pertes et le bruit, comparer à C1.
  C'est l'effet stimulé, classique ; l'effet quantique spontané demande du cryogénique et n'est pas l'objectif.

**Avant de commencer.** Une vraie recherche bibliographique sur les analogues électriques discrets : il est possible
que C2 existe déjà, et il vaut mieux le savoir d'abord.

**Ce qui compterait comme résultat.** C1 reproduit la loi thermique à 10 % dans le régime κ ≪ ω_max et documente
proprement la sortie de ce régime. C2 : une mesure de conversion stimulée, avec sa loi exponentielle en fréquence,
sur un banc RF décrit de façon reproductible.

**Critère d'arrêt.** Si C1 ne retrouve pas la loi exponentielle pour un front lent, c'est la simulation qui est
fausse, pas la physique : on corrige avant d'interpréter. Si les pertes d'une ligne réelle noient l'effet d'un
facteur 10 dans le devis de C2, on s'arrête à la simulation.

**Prérequis.** Les ondes dans une échelle LC (acquis), vitesse de groupe, effet Doppler, normes de modes (Robertson).

**À qui parler.** La communauté de la gravité analogue est petite, expérimentale et ouverte aux ingénieurs ; son
numéro de référence est celui des *Philosophical Transactions* de 2020 (Jacquet, Weinfurtner et König).

---

## 7. Comparatif

| | A — cordes | B — corde hadronique | C — gravité analogue |
|---|---|---|---|
| Nature | Compréhension | Petit projet de phénoménologie | Simulation, puis banc d'essai |
| Prérequis neufs | Mécanique analytique, relativité, quantique | Zwiebach 6 à 8, statistique d'ajustement ; champs pour la suite | Ondes dispersives, normes de modes |
| Recouvrement avec ton métier | Faible | Moyen (méthode du dossier, R120) | Fort (lignes, RF, mesure) |
| Premier jalon | A1 à A4 à la main, vers 6 mois | B1 reproduit, 3 mois après le chapitre 7 | C1, environ 3 mois |
| Chance d'un résultat publiable | Nulle, ce n'est pas le but | Faible à modérée | Modérée pour C1 ; réelle si C2 est neuf |
| Risque principal | Décrochage mathématique | Déjà fait par des professionnels | Déjà fait ; pertes et bruit |
| Lien avec « tout à partir d'un filament » | Direct : c'est la théorie qui le fait | La corde réelle de la nature, à l'échelle du fermi | Aucun : c'est un analogue, pas une ontologie |

---

## 8. Calendrier proposé (12 mois)

| Mois | Travail | Jalon |
|---|---|---|
| 0 à 1 | Chemin 0 | Note relue par deux pairs (décision Zenodo prise : pas de dépôt) |
| 1 à 4 | Socle (mécanique analytique, relativité, quantique) ; Zwiebach 1 à 5 | Tests de sortie des trois premières lignes du socle |
| 4 à 6 | Zwiebach 6 à 9 ; script B1 dès la fin du chapitre 7 | A1 à la main ; B1 reproduit la littérature |
| 6 à 9 | Zwiebach 10 à 14 ; B2 | A2 à A4 à la main ; √σ du potentiel avec sa fenêtre |
| 9 à 12 | B3, ou C1 si le banc d'essai t'attire plus que les tables | Une note courte, ou une simulation validée |
| 12 | Point de décision | Voir ci-dessous |

**Au mois 12, trois questions, écrites maintenant.** Sais-tu refaire A1 à A4 sans notes ? B1 reproduit-il la
littérature à 2 % ? Un pair humain a-t-il lu au moins un de tes calculs ? Trois oui : on continue sur B3 ou C2.
Moins de deux : on reste au socle et à la lecture, qui se suffisent.

**Ma recommandation est celle donnée en conversation** : A pour la culture, B comme premier projet, parce qu'il part de ce
qui a survécu de la CD et qu'il traite des cordes que la nature fabrique réellement. C est le chemin où ton métier
pèse le plus ; si au bout du socle le banc d'essai t'attire davantage que les tables du PDG, c'est celui-là qu'il
faut prendre, et il n'y aura rien à regretter.

---

## 9. Le rôle de l'IA dans ce cadre

Vérifier tes calculs à la main, pas les remplacer. Écrire avec toi les scripts de B1 à B3 et de C1, tests compris.
Chercher et résumer la littérature, en disant ce qui n'a pas pu être vérifié à la source. Jouer l'arbitre. La
règle 1 s'applique aux deux : rien n'entre au dossier que tu ne puisses expliquer toi-même.

---

## Références

Elles sont données de mémoire ; plusieurs sites, dont arXiv, n'étaient pas accessibles depuis la session qui a
rédigé ce document. Chacune est à vérifier à la source avant d'être citée. Les « références classiques » du dernier
paragraphe, ainsi que Lüscher 1981 et Chiao, Garmire et Townes 1964, ont été vérifiées par recherche le 22 novembre
(titres, volumes, pages) et ajoutées à la note.

**Socle et cordes.** L. Susskind et G. Hrabovsky, *The Theoretical Minimum* (2013) ; L. Susskind et A. Friedman,
*Quantum Mechanics: The Theoretical Minimum* (2014) et *Special Relativity and Classical Field Theory* (2017) ;
L. Susskind et A. Cabannes, *General Relativity: The Theoretical Minimum* (2023). D. Griffiths, *Introduction to
Quantum Mechanics*. J. Hartle, *Gravity* (2003). T. Lancaster et S. Blundell, *Quantum Field Theory for the Gifted
Amateur* (2014). B. Zwiebach, *A First Course in String Theory*, 2ᵉ éd. (2009). D. Tong, *Lectures on String
Theory*, arXiv:0908.0333. J. Polchinski, *String Theory* (1998). J. Scherk et J. H. Schwarz, Nucl. Phys. B 81, 118
(1974).

**Corde hadronique.** M. Lüscher, Nucl. Phys. B 180, 317 (1981). M. Lüscher et P. Weisz, JHEP 07 (2002) 049.
J. F. Arvis, Phys. Lett. B 127, 106 (1983). O. Aharony et Z. Komargodski, JHEP 05 (2013) 118. S. Dubovsky,
R. Flauger et V. Gorbenko, JHEP 09 (2012) 044 ; Phys. Rev. Lett. 111, 062006 (2013). A. Athenodorou, B. Bringoltz et
M. Teper, JHEP 02 (2011) 030. B. B. Brandt et M. Meineri, Int. J. Mod. Phys. A 31, 1643001 (2016). M. Caselle,
Universe 7, 170 (2021) ; M. Caselle et al., JHEP 01 (2015) 105. J. Sonnenschein et D. Weissman, JHEP 08 (2014) 013 ;
JHEP 02 (2015) 147. S. Hellerman et I. Swanson, Phys. Rev. Lett. 114, 111601 (2015). G. S. Bali, Phys. Rept. 343, 1
(2001). J. Bulava et al., Phys. Lett. B 793, 493 (2019) ; Phys. Lett. B 854, 138754 (2024). C. Gattringer et
C. B. Lang, *Quantum Chromodynamics on the Lattice* (2010). A. Chodos et C. B. Thorn, Nucl. Phys. B 72, 509 (1974).

**Gravité analogue.** W. G. Unruh, Phys. Rev. Lett. 46, 1351 (1981) ; Phys. Rev. D 51, 2827 (1995). S. Corley et
T. Jacobson, Phys. Rev. D 54, 1568 (1996). C. Barceló, S. Liberati et M. Visser, Living Rev. Relativ. 14, 3 (2011).
S. J. Robertson, J. Phys. B 45, 163001 (2012). S. Weinfurtner et al., Phys. Rev. Lett. 106, 021302 (2011).
J. Steinhauer, Nature Phys. 12, 959 (2016). J. R. Muñoz de Nova et al., Nature 569, 688 (2019). T. G. Philbin et
al., Science 319, 1367 (2008). J. Drori et al., Phys. Rev. Lett. 122, 010404 (2019). R. Schützhold et W. G. Unruh,
Phys. Rev. Lett. 95, 031301 (2005). P. D. Nation et al., Phys. Rev. Lett. 103, 087004 (2009). C. M. Wilson et al.,
Nature 479, 376 (2011). P. D. Nation, J. R. Johansson, M. P. Blencowe et F. Nori, Rev. Mod. Phys. 84, 1 (2012).
M. P. Blencowe et H. Wang, Phil. Trans. R. Soc. A 378, 20190224 (2020). M. J. Jacquet, S. Weinfurtner et F. König,
Phil. Trans. R. Soc. A 378, 20190239 (2020). R. Hirota et K. Suzuki, Proc. IEEE 61, 1483 (1973). M. Remoissenet,
*Waves Called Solitons* (3ᵉ éd., 1999). G. E. Volovik, *The Universe in a Helium Droplet* (2003).

**Références classiques à ajouter à la note des résultats négatifs.** G. H. Derrick, J. Math. Phys. 5, 1252 (1964).
W.-T. Ni, Astrophys. J. 176, 769 (1972). K. Nordtvedt et C. M. Will, Astrophys. J. 177, 775 (1972). C. M. Will,
*Theory and Experiment in Gravitational Physics* (2ᵉ éd., 2018). N. G. Vakhitov et A. A. Kolokolov, Radiophys.
Quantum Electron. 16, 783 (1973). A. Chodos, R. L. Jaffe, K. Johnson, C. B. Thorn et V. F. Weisskopf, Phys. Rev. D
9, 3471 (1974).
