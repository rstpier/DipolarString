# Redraw — la base, énoncée par l'auteur

*(English: a fresh foundation for the matter sector, written from the author's statements of
15 September 2026, independent of the V2.10 manuscript. Each statement is his; under it, the one
precision needed to compute with it. Nothing here is derived yet.)*

**Règle du jeu.** Le manuscrit V2.10 n'est plus la référence. Ses hypothèses écartées ici :
corde = conducteur ; fluide = charge répartie sur la corde ; jonction = contact ponctuel sans
angle ; « Γ_pole = 1/3 fixe e/3 » ; N_p = 27. Ce qui est gardé de la session : les outils de
calcul (électrostatique de lignes, modes d'anneaux, holonomies) et les bornes mesurées
(neutralité des atomes, PVLAS, coupure de Brillouin si elle survit).

## Les énoncés (R1–R8)

**R1 — La corde et son fluide.** Une corde contient un fluide EM polarisé, avec un pôle + et un
pôle − à ses extrémités.
→ À préciser : le fluide a-t-il une charge nette (par corde), ou seulement ses deux pôles ?
Se déplace-t-il le long de la corde, à quelle vitesse, et qu'est-ce qui le retient aux pôles ?

**R2 — L'assemblage.** Les cordes s'assemblent sur la longueur par leurs pôles (+ contre −).
→ Calculé : bout à bout, une jonction lie dès que les pôles valent 0,3 de la charge nette ;
côte à côte seulement à charges opposées. À préciser : la jonction impose-t-elle un angle ?

**R3 — La mère.** Toute particule vient d'un assemblage de DQD de spin 1 (symétriques) ; la mère
est une chaîne ouverte de DQD, `=:=:=:=:=:=:=` ; trois DQD font `===`.
→ À préciser : `=` est-il deux cordes parallèles liées côte à côte, ou une corde repliée ?
`:` est-il une jonction pôle à pôle ?

**R4 — La brisure.** Une chaîne de 3 DQD se brise en deux par une action unique : un positon et un
électron, rotations opposées, `(+)———   ———(−)`.
→ Calculé : si la mère porte un mode d'enroulement 1, chaque fille hérite d'un demi-enroulement
(même fréquence exacte). À préciser : « rotation » = circulation du fluide dans la corde, ou
rotation de la corde entière dans l'espace ?

**R5 — La charge active.** Le fluide s'accumule dans le pôle suiveur et crée la charge active.
→ Contrainte mesurée : la charge ne peut dépendre du mouvement dans le substrat (neutralité des
atomes, 10⁻²¹). À préciser : « suiveur » par rapport à quoi (le mouvement de la corde, la
circulation du fluide, la rotation) ? et qu'est-ce qui empêche le fluide accumulé de repartir ?

**R6 — Le spin.** Le spin est le nombre d'axes de symétrie divisé par deux : `—` et `=` en ont
deux (spin 1), `———(−)` en a un (spin ½).
→ À préciser : axes de symétrie de la forme au repos, ou du mouvement ? Que vaut la règle pour un
anneau fermé et pour un triangle ?

**R7 — La courbure.** Chaque liaison de corde produit une masse de liaison qui engendre une
courbure de déplacement ; c'est une conséquence mécanique simple.
→ À préciser en une ligne : quelle masse (celle de la jonction ?), quel déplacement (la corde qui
avance ? le fluide ?), et quelle force courbe la trajectoire.

**R8 — L'équilibre de charge.** Ce qui force les quarks à tenir ensemble est l'interaction EM
universelle : tout soliton doit s'équilibrer à la charge q.
→ À préciser : q = e exactement ? et par quel mécanisme un anneau à ±1/3 est-il instable seul
(dans un milieu linéaire il ne l'est pas).

## Ce que la base doit produire, dans l'ordre

1. **La charge** : où elle est, pourquoi e/3 par corde, pourquoi elle ne dépend pas du mouvement.
2. **Le spin** : pourquoi ½ pour la fille et 1 pour la mère, et le signe −1 sous un tour complet.
3. **La masse** : ce qu'elle est (mode du fluide ? masse de liaison des jonctions ? les deux ?) et
   le rapport m_p/m_e.
4. **La stabilité** : pourquoi la fille ne se referme pas en anneau et ne se dissout pas.

## Méthode

Chaque énoncé précisé devient un calcul avec un script et un verdict (DÉRIVÉ / CONDITIONNEL /
POSTULÉ / EXCLU), comme pour le reste du dossier `matter/`. Rien n'est repris du manuscrit sans
être re-dérivé de R1–R8.

---

## Premier calcul sur la base : R3–R6 en mécanique (`spin_rod.py`, 8/8)

Lecture adoptée, à corriger si fausse : « rotation » = rotation de la fille dans l'espace autour
du point de brisure ; « pôle suiveur » = le bout éloigné de ce point, où la force centrifuge
chasse le fluide. Barre rigide de longueur λ_C (trois cordes), masse m_e, fluide de charge −e.

| ce que ta mécanique donne | valeur |
|---|---|
| rotation pour un spin ½ (pivot au point de brisure) | ω = 2,9×10¹⁹ rad/s, bout éloigné à 0,24 c |
| force centrifuge sur le fluide / sa répulsion de Coulomb | 49 : le fluide s'accumule bien au pôle éloigné (R5) |
| indépendance de la charge au mouvement | vraie si l'accumulation est complète (neutralité des atomes, 10⁻²¹) |
| axes de symétrie | 1 pour la fille, 2 pour la mère (R6) |
| **facteur g** | **3** pour une barre uniforme avec la charge au bout ; 1 si la masse est au pôle avec la charge ; **2,0023 seulement si la moitié de la masse est au pôle et l'autre au pivot** |
| rayonnement de Larmor | 25 W : l'énergie de repos part en 3×10⁻¹⁵ s |

Ce que la base doit encore fournir : pourquoi la charge en rotation ne rayonne pas (un mode
stationnaire réactif, pas une forme), et le −1 sous un tour complet, qu'une rotation classique ne
donne pas. La condition sur g est la première prédiction que la base peut être tenue de vérifier.

**R9 — La particule de l'espace.** Le DQD est la particule stable de l'espace.
→ À préciser : sa masse (nulle ?), sa rotation propre, et pourquoi il est stable là où la fille
ne l'est pas.

**R10 — La vitesse de rotation.** La vitesse de rotation est stabilisée par Z₀ de l'espace.
→ Lecture calculée (`rod_z0.py`, 6/6) : le fluide est EM, il va à c₀ (c₀ et Z₀ viennent des mêmes
L₀, C₀), donc le pôle éloigné de la fille tourne à c₀. Avec S = ℏ/2 et la répartition de masse
que g = 2 impose (moitié au pôle, moitié au pivot), la fille mesure exactement ℏ/(m_e c) = 386 fm
et son moment vaut exactement un magnéton de Bohr : trois nombres de l'électron (spin, taille de
Compton, μ_B) sortis d'une seule barre, sans le manuscrit. Dette inchangée et plus aiguë : une
charge qui tourne à c₀ rayonne 3×10⁵ W classiquement (3×10⁻¹⁹ s) ; « stabilisée par Z₀ » doit
vouloir dire « mode lié, non rayonnant », et c'est cela qu'il faut écrire.

**R11 — Le rayonnement.** La vitesse c₀ annule son propre rayonnement pour un observateur
extérieur, et s'éteint instantanément dans son référentiel.
→ Calculé (`no_radiation.py`, 6/6) : exact sous une condition. Une charge ponctuelle au pôle
tournant à c₀ rayonne sans limite (Larmor × γ⁴, 10⁻¹⁹ s) ; une onde de charge qui circule sur
l'anneau est une antenne boucle résonante et se vide en dix tours ; **une charge répartie
uniformément sur l'orbite avec un courant continu est une source stationnaire et ne rayonne
jamais**, en gardant μ = μ_B et S = ℏ/2 (moitié de la masse sur l'anneau à c₀). Donc R11 fixe la
lecture de R5 : à c₀ le fluide *est* toute l'orbite pour l'observateur extérieur, pas un point ;
et R6 doit être relu pour un anneau de courant (aucun axe miroir dans le plan, un axe de rotation).
Cohérent avec « s'éteint dans son référentiel » : un fluide guidé à c₀ n'a pas de temps propre,
il n'émet pas, il ne peut que fuir, et une source stationnaire ne fuit rien.

**R1, précisé — Le fluide est sans masse** (pour le moment ; problème ouvert).
→ Calculé (`massless_fluid.py`, 8/8) : un fluide sans masse guidé à c₀ porte S = R·E_circ/c et
μ = e·c·R/2, donc **g = 1 si toute l'énergie de repos circule avec la charge, et g = 2,0023 si et
seulement si la moitié circule et la moitié est statique**. Alors S = ℏ/2 fixe R = ℏ/(m_e c)
= 386 fm exactement, μ = μ_B exactement, et le quantum qui circule (m_e c²/2) est le mode
d'enroulement ½ : l'ancrage R₃ de l'ancien manuscrit et l'anneau antipériodique sont un seul
objet. La moitié statique n'est pas l'énergie de champ de l'anneau (0,4 à 0,9 %) ; si c'est la
masse de liaison des jonctions (R7), il faut 128 keV par jonction pour deux jonctions ou 85 keV
pour trois, soit des pôles de 1,2 à 1,5 e. Tension avec R11 : un anneau stationnaire (ω = 0) ne
tient que ½LI² ≈ 2 keV ; la masse est donc dans un mode (ω ≠ 0), qui ne rayonne pas seulement
si son champ est confiné dans un guide fermé. **La base doit dire où vit l'onde : dans une paire
(le DQD, R9) ou sur la fille ouverte.**

**R7, précisé par l'auteur — pour n = 3, deux liaisons portent la masse.**
→ Calculé (`junction_mass.py`, 7/7) : 511/2 = 255,5 keV par jonction (260,5 × 2 = 521 keV
dépasserait m_e c² de 2 %). Une liaison de 255 keV pôle à pôle demande des pôles de ~2,1 e. Mais
une masse **statique** dans les jonctions ne tourne pas (S = 0, μ = 0) : cette énergie est celle
du **mode** aux jonctions, en partie circulante. Pour un fluide sans masse, g = (1 + ρ)/(1 − ρ)
avec ρ le rapport des puissances arrière/avant créé par les réflexions aux pôles : **g = 2 ⇔ les
pôles réfléchissent un tiers de la puissance** (|Γ| = 0,58 en amplitude). L'ancien Γ_pole = 1/3
en amplitude (ρ = 1/9) donnerait g = 1,25. Le rayon reste ℏ/(m_e c). Prédiction tenue par la
base : réflectivité des pôles 1/3 en puissance.

**R12 — La fille est un 3/4 de spire, et il contribue.**
→ Calculé (`three_quarter_turn.py`, 7/7), image P2 : moitié statique de la masse dans les deux
jonctions (128 keV chacune, pôles ~1,5 e), moitié dans le mode qui vit sur l'arc et se réfléchit
aux deux pôles. S = ℏ/2 et la résonance du fondamental sur l'arc donnent **R = (4/3)·ℏ/(m_e c)
= 515 fm et une fraction circulante exactement 3/4** : pour un arc de fraction a, la fraction
circulante est f = a, le 3/4 de spire contribue 3/4. μ = μ_B et g = 2 tiennent (ils ne demandent
que f·R = ℏ/(m_e c) et la moitié statique). Les pôles doivent alors réfléchir **1/7 de la
puissance** (|Γ| = 0,38) ; l'ancien Γ_pole = 1/3 correspondrait à un 4/5 de spire. Deux images
restent en lice : P1 (toute la masse dans le mode aux jonctions, ρ = 1/3, R = ℏ/(m_e c)) et P2
(moitié statique, arc 3/4, ρ = 1/7, R = 515 fm). Elles diffèrent par le rayon et par la
réflectivité des pôles ; c'est là que la base doit choisir.
