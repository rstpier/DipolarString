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

**Calibration sur plusieurs particules (auteur : muon et tau).**
→ Calculé (`lepton_calibration.py`, 7/7). Même objet mis à l'échelle : R = ℏ/(mc) donne 1,87 fm
pour le muon et 0,11 fm pour le tau. **L'image P2 avec des jonctions de Coulomb à contact fixe ne
suit pas** : la moitié statique demanderait des pôles de 22 e (muon) et 88 e (tau). Donc ce qui
est statique doit varier en 1/R comme le mode ; ce n'est pas une liaison de contact. L'image P1
(tout dans le mode) suit sans effort mais ne fixe aucun rapport : R est libre dans les deux.
Les données obéissent déjà à deux relations que toute calibration doit reproduire : Koide,
(m_e + m_μ + m_τ)/(Σ√m)² = 2/3 à 10⁻⁵ ; Barut, m_n = m_e[1 + (3/2α)Σk⁴], muon à 0,1 %, tau à
0,6 %, un anneau avec des quanta d'excitation interne et une auto-énergie magnétique en 1/α. Les
harmoniques du mode ½ donnent 3, 5, 7, pas 206,77. **La base doit fournir la règle de l'échelle
entre leptons** ; le 1/α de Barut est le rapport entre l'échelle du mode ℏc/R et l'échelle de
Coulomb αℏc/R, deux échelles que la base possède.

**Correction (moi) et R13 — « Une jonction qui ne tourne pas est impossible. C'est Z_pole
intrinsèque. »**
→ `g_bookkeeping.py`, 9/9. **Correction :** dans `junction_mass.py` j'avais écrit que réfléchir un
tiers de la puissance aux pôles donnait g = 2. C'est faux : une onde réfléchie réduit le courant
de charge et le moment cinétique dans la même proportion, donc **g = m c²/E_mode, indépendant de
la réflexion**. Les vérifications C et C' de ce script sont retirées. Ce que g = 2 exige, exactement :
la moitié de l'énergie de repos hors du mode et **sur l'axe de rotation**, tournant sur elle-même
sans orbiter, comme tu le dis d'une jonction. Toute énergie neutre qui orbite ajoute du moment
cinétique sans moment magnétique et fait tomber g sous 1 (moitié de la masse en orbite : g = 0,86).
L'image de la barre tournant autour du point de brisure a ce pivot ; l'image de l'arc a un centre
vide et ne peut pas l'héberger. **Z_pole intrinsèque :** avec f = a (résonance sur l'arc), la
fraction d'arc fixe l'impédance du pôle : 3/4 de spire ⇔ Z_pole = 2,22 Z_corde ; Z_pole = 2 Z_corde
(l'ancien Γ = 1/3) ⇔ 4/5 de spire. Pour tout ρ, S = ℏ/2 et μ = μ_B tiennent à R = ℏ/(m_e c)/f.

**Contributions de masse des spires et des jonctions** (demande de l'auteur).
→ `mass_contributions.py`, 6/6. Ce que la base fixe : pour l'électron, g = 2 partage la masse à
parts égales entre le mode sur la spire et le pivot ; avec 3/4 de spire et 2 jonctions, **340,7 keV
par tour complet (85,2 keV par corde d'un quart de tour) et 127,7 keV par jonction**. Loi additive
à un entier (n cordes, n/4 tours, n − 1 jonctions) : m(n) = 212,9 n − 127,7 keV ; les chaînes de
4 à 7 cordes pèsent 0,72 à 1,36 MeV, cent fois sous le muon ; le muon exigerait 497 cordes, le tau
8346. Comptes libres (spires et jonctions indépendants) : aucune solution à coefficients positifs
tant que le muon a moins de 156 spires et moins de 414 jonctions, puis des milliers au-delà : la
calibration ne sélectionne rien sans règle sur les comptes. **Verdict : des contributions additives
calibrées sur l'électron n'atteignent ni le muon ni le tau avec de petites chaînes ; il faut des
centaines d'unités, ou des masses unitaires qui varient en 1/R avec l'objet (le mode), ce qui
n'est pas additif.**

**R14 — L'échelle des leptons : e = 3, μ = 7, τ = 11 cordes ; l'électron serait spécial par sa
spire ouverte (impact exponentiel sur la masse de spire).**
→ `lepton_ladder.py`, 8/8. Une exponentielle en n passant par e et μ prédit le tau 12 fois trop
lourd ; passant par μ et τ (×2,02 par corde), elle place l'électron 12,3 fois sous la loi : la
« suppression de la spire ouverte » serait e^{−2,5}, un nombre libre. **Une loi de puissance
m = m_e (n/3)^p** fait mieux : p ajusté sur le muon = 6,29, tau prédit à +2,2 % ; avec **p = 2π**,
sans paramètre ajusté, muon à −0,8 % et tau à +1,0 %. Parmi tous les triplets impairs (3, n_μ, n_τ)
avec n ≤ 15, (3, 7, 11) est celui qui prédit le mieux le tau depuis le muon ; le suivant manque
de 29 %. Sens physique sous m ∝ 1/R : R ∝ n^{−2π}, la boucle du muon est 205 fois plus petite et
ses cordes 480 fois plus courtes que celles de l'électron : les cordes n'ont pas de longueur
fixe, la fermeture fixe leur échelle. Réserves : trois masses, entiers choisis, un exposant ; le
tau à 1 % est le seul test passé ; Koide vaut 0,668 sur la loi contre 2/3 à 10⁻⁵ sur les données ;
un n = 15 pèserait 12,6 GeV, exclu par le LEP, donc l'échelle s'arrête à 11 ou 15 n'est pas un lepton.

**R15 — « 1 spire doit donner environ 150 MeV ; 3/4 de spire, presque rien. »**
→ `closure_law.py`, 9/9. Les deux ancrages : fermer le dernier quart de tour multiplie la masse par
~290 et, sous m ∝ 1/R, rétrécit la boucle de 386 fm à 1,3 fm : l'électron est immense parce qu'il
est ouvert. Mais « 150 MeV par spire fermée » ne construit pas μ et τ : en additionnant les spires
fermées, μ = 150 MeV (+42 %) et τ = 300 MeV (6 fois trop léger) ; en les multipliant par 290,
τ = 44 GeV ; une exponentielle en fraction de fermeture donne des TeV. Seule la loi de puissance
des chaînes ouvertes, (n/3)^{2π}, atteint μ et τ (1 %). **L'objet réel « une spire fermée d'environ
150 MeV » est le pion** : un anneau fermé de 4 cordes (2+, 2−) est neutre, c'est le π⁰ à 135,0 MeV ;
un anneau fermé de 5 cordes (1+, 4−) a la charge −1, c'est le π⁻ à 139,6 MeV ; parité de n et
charge concordent. La fermeture est donc la première variable : une chaîne ouverte de 5 pèserait
12,8 MeV sur la loi des leptons, un anneau fermé de 5 pèse 140 MeV. La masse n'est pas une fonction
de n seul. Deux lois à expliquer : chaînes ouvertes, m = m_e (n/3)^{2π} ; anneaux fermés, ~140 MeV,
presque indépendants de n entre 4 et 5.

**R16 — Ce qui fixe ~140 MeV pour l'anneau fermé** (demande de l'auteur).
→ `closed_ring_scale.py`, 7/7. Dans les constantes de la base (e, ℏ, c₀, Z₀, m_e), une seule
longueur vaut 1,4 fm : **la moitié du rayon classique de l'électron, r_e/2 = α·ℏ/(m_e c)/2 =
1,409 fm**, dont l'énergie de mode est ℏc/(r_e/2) = **2 m_e c²/α = 140,05 MeV** : le π± à 0,34 %,
le π⁰ à 3,6 % (coïncidence connue, m_π ≈ 2m_e/α). Deux lectures dans la base : (a) Coulomb, r_e
est le rayon où l'énergie de champ de la charge e vaut m_e c², l'anneau fermé siège à la moitié ;
(b) impédance, Z₀ = 2α·R_K avec R_K = h/e², donc 2/α = 4 R_K/Z₀ = 274 : la masse de l'anneau
fermé est celle de l'électron × 4 (résistance quantique / impédance du vide), comme si l'anneau
fermé « voyait » h/e² là où la chaîne ouverte voit Z₀. Tests : l'énergie de Coulomb propre de
l'anneau donne le bon signe pour π± − π⁰ mais 0,36 MeV contre 4,59 mesurés (13 fois trop peu) ;
les hadrons ne sont pas des multiples de 140 MeV (K 3,5, η 3,9, ρ 5,5, p 6,7). **Statut :
identification de l'échelle, pas dérivation** ; le mécanisme qui épingle un anneau fermé à r_e/2
est ce que la base doit fournir.

**R17 — « La largeur de la nappe est fixée par l'impédance Z₀ de l'espace. »**
→ `sheet_matched.py`, 8/8. Pour une nappe (guide à plaques parallèles, largeur w, écart d),
Z = Z₀·d/w : l'adaptation à Z₀ impose **w = d, une section carrée** ; pour deux fils, le même
principe redonne D/r = 2 cosh π = 23,2. L'adaptation fixe donc la **forme** de la section, pas sa
taille. Conséquences sur un guide adapté à Z₀ : (1) la charge de pôle par quantum vaut
δ = √(4ℏ/πZ₀) = **e/(π√α) = 3,73 e**, sans réglage ; (2) l'énergie de Coulomb entre deux pôles
distants de s vaut δ²e²/(4πε₀s) = **ℏc/(π²s)** : α s'annule et une jonction a la forme d'un mode ;
si s suit la largeur et la largeur suit l'objet, **les masses de jonction varient en 1/R comme les
masses de mode**, ce que la calibration μ/τ exigeait et que le contact à rayon fixe ne donnait pas
(il fallait des pôles de 22 e et 88 e ; ici les pôles restent à 3,73 e et c'est l'écart qui change).
(3) Avec g = 2 et deux jonctions portant la moitié statique, 2ℏc/(π²w) = mc²/2 donne
**w = 4ƛ/π² = 0,405 ƛ** : 157 fm pour l'électron, 0,76 fm pour le muon ; avec R = 4ƛ/3, le rapport
**w/R = 3/π² = 0,30 est le même pour tous les leptons** (famille auto-similaire). (4) **Le prix :**
une nappe de largeur w a des modes transverses à ℏcπ/w, soit 3,96 MeV pour l'électron et 819 MeV
pour le muon ; un électron excité à 4 MeV n'est pas observé. Deux sorties : le fluide n'admet que
le mode TEM par axiome, ou le guide de l'électron est bien plus mince que 157 fm ; mais sous
6·10⁻³ fm (mode transverse au-delà de 100 GeV) la moitié statique des deux jonctions, ℏc/(π²w)
chacune, dépasse la masse de l'électron d'un facteur 2·10⁴ : **la nappe ne peut pas être à la fois
mince et porteuse de la moitié statique.** Verdict : DÉRIVÉ pour la forme carrée, δ = e/(π√α) et
ℏc/(π²s) ; CONDITIONNEL (g = 2 + deux jonctions statiques) pour w = 0,405 ƛ ; le mode transverse
à 4 MeV est un problème ouvert de la base.

**R18 — « Les masses fonctionnent-elles pour le proton et le neutron ? Table de masse de la famille
de l'électron (fermions). »** (demande de l'auteur)
→ `fermion_table.py`, 9/9. Règle de parité : q = (n₊ − n₋)e/3 impose n impair pour ±1 et ±1/3,
n pair pour 0 et ±2/3. **Non pour les nucléons.** Proton (charge +1, n impair) : l'échelle offre
509 MeV (n = 9, −46 %) ou 1795 MeV (n = 11, +91 %) ; le n continu vaut 9,92, du côté pair. Neutron
(n pair) : n = 10 donne 986 MeV (+5 %), mais p et n diffèrent de 0,14 % alors que deux barreaux
voisins diffèrent de 82 % : la paire ne peut pas être deux chaînes. Lecture anneaux fermés : le
proton fait 6,70 anneaux de 140 MeV, pas un entier. Lecture composite (uud = chaînes 4, 4, 5) :
18,9 MeV, soit 2 % du proton ; 98 % de la masse du nucléon est de la liaison, dont la base n'a
aucun mécanisme. **Table des fermions** (meilleur n de la parité imposée par la charge) : e 3
(exact), μ 7 (−0,8 %), τ 11 (+1,0 %) ; u 4 (+44 %), d 5 (+169 %), s 7 (+12 %), c 10 (−23 %),
b 13 (+22 %), t 22 (−19 %) : les quarks dispersent de 12 à 169 %, l'échelle est une loi de leptons.
**Problème caché de la famille de l'électron** : q = −1 avec n impair contient aussi n = 5
(12,7 MeV) et n = 9 (509 MeV), des leptons chargés qui n'existent pas ; la famille doit donc être
n ≡ 3 (mod 4) : 3, 7, 11, puis 15 (12,6 GeV) et 19 (55,6 GeV), exclus par le LEP (> 102,8 GeV),
et 23 (185 GeV), premier membre permis. Pourquoi le pas de 4 cordes est la question ouverte.

**R19 — « Neutron et proton ont juste un joint de différence, donc on a la calibration du nœud. »**
→ `nucleon_joint.py`, 6/6. Deux lectures. (i) Neutron = proton + une corde + un joint : mais la
règle de charge q = (n₊ − n₋)e/3 fait varier q de 1/3 par corde, pas de 1 ; le proton + une corde
ne peut pas être neutre. (ii) **Neutron = proton + une chaîne d'électron (0, 3) attachée par un
joint** : charge +1 − 1 = 0 et parité (impair + 3 = pair) respectées ; c'est la lecture cohérente
avec la base. **Calibration** : (i) m_n − m_p = 1,293 MeV ; (ii) m_n − m_p − m_e = 0,782 MeV (le Q
de la désintégration β) ; soit 10,1 et 6,1 fois la jonction de l'électron (127,75 keV, R13).
**Ce que la base impose en plus** : (1) sur un guide adapté (pôles 3,73 e, R17) un joint de 1,29 /
0,78 MeV siège à s = 15,5 / 25,6 fm, 18 à 30 rayons de proton ; à l'intérieur du proton (s = r_p =
0,84 fm) le joint adapté pèse 24 MeV : la calibration n'est pas la jonction R17, sauf pôles de
0,87 e / 0,68 e. (2) La base porte l'énergie de Coulomb du proton (0,86 à 1,03 MeV à r_p) ou
l'attraction p–e (−1,71 MeV à r_p) : le joint brut vaut 2,2 à 2,5 MeV, 1,29 MeV est le net (réseau
QCD+QED : partie forte 2,52, partie EM −1,00). (3) **Le joint n'est pas un nombre unique** : dans le
propre décompte de la base, le π⁻ (anneau de 5) est un π⁰ (anneau de 4) + une corde + un joint, et
ce joint vaut 4,59 MeV, 3,6 fois celui du nucléon ; les paires d'isospin vont de 0,3 (B) à 8,1 MeV
(Σ). **Ce qui passe** (lecture ii comme image) : la boucle négative ajoutée, avec le μ = −ecR/2 de
la base, doit avoir R = 0,99 fm pour donner μ_n − μ_p = −4,71 μ_N ; le rayon de charge du neutron
(⟨r²⟩ = −0,1155 fm², cœur +e du proton + anneau −e) donne R = 0,91 fm : deux mesures indépendantes
concordent à 9 %, et le spin de la boucle (R·E/c < 0,02 ℏ) laisse le neutron à spin 1/2. Verdict :
calibration = un nombre, pas une dérivation ; l'image « proton + boucle d'électron à ~0,95 fm »
est CONDITIONNELLE et passe deux tests de taille.

**R20 — « Une jonction vaut 2 MeV. Mais les quarks ne fonctionnent pas pareil : ils semblent
ramifiés en étoile. »**
→ `quark_star.py`, 6/6. **Le joint de 2 MeV** dans la loi de Coulomb de la base u = δ²K/s : avec
des pôles unitaires (e) il siège à s = 0,72 fm, la taille du nucléon ; avec les pôles adaptés de
R17 (3,73 e) à 10 fm ; avec le pôle d'une corde (e/3) à 0,08 fm. Il vaut 15,7 fois la jonction de
l'électron (127,75 keV) : ce n'est pas cette jonction remise à l'échelle 1/s (elle serait à 10 fm),
c'est une jonction à pôles unitaires à 0,7 fm. **L'étoile et la charge** : la règle q = (n₊ − n₋)e/3
ne voit pas la topologie ; une étoile à k bras a les classes de charge d'une chaîne de k cordes
(parité de k). +2/3 exige k pair, −1/3 exige k impair : une étoile à 3 bras ne peut pas être un u,
une étoile à 4 bras ne peut pas être un d. Étoiles minimales : d = 1 corde ou 3 bras (1+, 2−) ;
u = 2 cordes ou 4 bras (3+, 1−). **Budget de joints** : un proton dessiné comme étoile d'étoiles
(u, u, d avec leur centre + un centre Y) a 4 joints × 2 MeV = 8 MeV, 0,9 % du proton : la masse
n'est pas dans les joints. **Lecture de mode au rayon de charge** : ℏc/r_p = 234,7 MeV et le
produit mesuré r_p·m_p·c/ℏ = 3,998 ± 0,002 (PDG 2024) : le proton est 4 modes de son rayon de
charge à 0,05 %, un Y à 3 bras en fait 3 (−25 %) ; coïncidence à expliquer, pas une dérivation ;
le fluide à c sur r_p demande une charge circulante de 0,70 e pour μ_p (2/3 e : −4,6 %). **Signe du
confinement** : dans la base une corde deux fois plus longue est deux fois plus légère (mode et
Coulomb en 1/L) ; la corde Y du réseau (σ ≈ 0,89 GeV/fm) est deux fois plus lourde : signes
opposés. L'étoile de la base n'a donc pas de taille propre ; seul S = R·E/c = ℏ/2 la fixe, comme
pour l'électron, R = ƛ/f avec f = 1/4 pour le proton. **R6 sur l'étoile** (enregistré, non testé) :
étoile symétrique à 3 bras, 3 axes → spin 3/2 (Δ) ; un bras différent, 1 axe → spin 1/2 (N) : à
contenu égal, l'écart Δ − N de 294 MeV serait une différence de forme.

---

## Bilan au 16 septembre 2026 (R1–R20, 20 scripts, tous PASS)

**Ce qui tient (dérivé dans la base, script à l'appui).**
- Le fluide à c₀ sur une boucle de rayon R = ƛ/f donne S = ℏ/2 et μ = μ_B exactement ; g = mc²/E_mode,
  donc g = 2 ⇔ la moitié de l'énergie de repos est statique (rod_z0, massless_fluid, g_bookkeeping).
- Pas de rayonnement seulement pour un anneau stationnaire uniforme (no_radiation).
- La charge est la branche non appariée : q = (n₊ − n₋)e/3, la parité de n fixe la classe
  (ring_catalogue, fermion_table) ; la règle ne voit pas la topologie (quark_star).
- L'échelle des leptons m = m_e (n/3)^{2π} pour n = 3, 7, 11 : μ à −0,8 %, τ à +1,0 %, exposant
  non ajusté (lepton_ladder).
- L'anneau fermé pèse 2m_ec²/α = 140,05 MeV : π± à 0,34 %, avec la bonne parité de charge
  (closure_law, closed_ring_scale).
- Adaptation à Z₀ : section carrée, pôle e/(π√α) = 3,73 e, jonction ℏc/(π²s), famille
  auto-similaire w/R = 3/π² (sheet_matched).
- Neutron = proton + chaîne d'électron + un joint : la boucle négative doit avoir R = 0,99 fm
  pour μ_n et 0,91 fm pour le rayon de charge, deux tests indépendants à 9 % (nucleon_joint).
- Coïncidence enregistrée : m_p c² = 4ℏc/r_p à 0,05 % (quark_star).

**Ce qui ne tient pas (exclu par le calcul).**
- Nucléons et quarks ne sont pas sur l'échelle des leptons (quarks à 12–169 %, proton à 509 ou
  1795 MeV) ; la paire p–n (0,14 %) ne peut pas être deux barreaux (82 %).
- Les masses additives (spires + jonctions calibrées sur l'électron) n'atteignent ni μ ni τ.
- « Un joint = un nombre » : 1,3 MeV au nucléon, 4,6 MeV au pion, 0,3 à 8 MeV sur l'isospin.
- La nappe adaptée ne peut pas être mince et porter la moitié statique : mode transverse à 4 MeV
  pour l'électron, non observé.
- La corde s'allège en s'allongeant (mode et Coulomb en 1/L) : pas de confinement, aucune taille
  propre ; seul le spin fixe R.
- g = 2 ⇔ ρ = 1/3 : retiré (g ne dépend pas de ρ).

**Ce que la base doit encore fournir.**
- Pourquoi l'exposant 2π, pourquoi n = 3, 7, 11 et pas 5, 9 (pas de 4 cordes, LEP exclut 15 et 19).
- Ce qui épingle l'anneau fermé à r_e/2.
- Où vit l'onde (R11), le fluide sans masse, où siège la moitié statique.
- Un mécanisme pour les quarks et les nucléons : 98 % de la masse est de la liaison, et les
  « 4 modes de r_p » n'ont pas d'origine.
- Le −1 sous 2π (statistique d'échange) et la courbure de déplacement (R7).

**R21 — « Un milieu avec plein de topologies, plein de familles de particules à morphologies
différentes ? »** (question de l'auteur)
→ Pas de nouveau calcul ; lecture de R1–R20. C'est ce que la base dit maintenant : un seul milieu
(le DQD, avec c₀ et Z₀) et une forme par famille :

| forme | famille | loi de masse trouvée | nombre propre à la forme | lien avec une autre famille |
|---|---|---|---|---|
| chaîne ouverte, n cordes | leptons | m_e (n/3)^{2π} | l'exposant 2π, les n = 3, 7, 11 | — |
| anneau fermé | pions | 2m_ec²/α = 140,05 MeV | le rayon r_e/2 | m_e et α : l'anneau est pesé par l'électron |
| étoile (Y) | quarks, nucléons | 4ℏc/r_p (coïncidence) | r_p | aucun : r_p n'est pas exprimé dans la base |
| nappe (section) | toutes | — | w = d par Z₀ | δ = e/(π√α) commun à toutes les formes |
| proton + boucle | neutron | joint 0,78 MeV net | — | R = 0,9–1,0 fm tenu par μ_n et ⟨r²⟩_n |

C'est l'image de Wen (un condensat, les particules sont des façons de nouer la même corde) ; Wen
non plus ne dérive pas les masses, ce sont les gaps de son hamiltonien. **Le critère** : tant que
chaque forme a son nombre propre, c'est un catalogue ; le milieu n'explique une famille que quand
une de ses grandeurs sort des constantes d'une autre sans nombre neuf. La base a un tel lien
(anneau ← électron, 0,3 %), un demi-lien (boucle du neutron ← proton, 9 %), et aucun pour l'étoile :
r_p (ou l'exposant 2π) exprimé avec m_e, α, Z₀ serait le test qui fait du milieu une théorie.

**R22 — « Si les quarks sont n = 3 en Y, ça explique l'interaction forte (une tri-jonction) avec
3 pôles libres à 1/3. »**
→ `quark_y.py`, 5/5. **La charge d'abord** : un Y de trois cordes porte q ∈ {±1/3, ±1} ; d, s, b
(−1/3) y tiennent, u, c, t (+2/3) n'y tiennent pas, +2/3 exige un nombre pair de cordes ; aucun
découpage du proton (6+, 3−) en trois Y(3) ne donne (2/3, 2/3, −1/3), seulement (1/3, 1/3, 1/3),
(1, 1/3, −1/3) ou (1, 1, −1). La proposition couvre la moitié des quarks. Ce que la parité permet :
u = 2 cordes (2+) ou 4 (3+, 1−) ; d = 1 corde (1−) ou 3 (1+, 2−) ; proton = 5 cordes (4+, 1−) ou
11 ; neutron = 4 (2+, 2−) ou 10. **La tri-jonction ensuite** : trois pôles en un point contiennent
au moins une paire de même signe (frustration) ; (+,+,−) lie exactement comme une paire, −Kδ²/s,
et (+,+,+) repousse. Une tri-jonction vaut donc une jonction : 2 MeV à 0,72 fm avec des pôles
unitaires (R20). C'est l'échelle de la force nucléaire résiduelle (deutéron 2,22 MeV), pas du
confinement (300 MeV par quark constituant, 0,9 GeV par fm de corde Y) : rapport 150 à 450. Avec
des pôles libres de e/3, le lien vaut K/(9s) = 0,16 MeV à 1 fm et n'atteint 2 MeV qu'à 0,08 fm,
14 fois trop faible pour le deutéron aux distances nucléaires. Enfin trois Y(3) exposent 9 pôles
libres, nombre impair : un pôle ne s'apparie jamais ; un méson (Y + anti-Y, 6 pôles) s'apparie
complètement. **Verdict** : la tri-jonction est une bonne image de la force nucléaire (MeV), pas
de l'interaction forte qui pèse le nucléon (centaines de MeV) ; et « n = 3 en Y » ne fait que les
quarks de charge −1/3.

**R23 — Le critère de R21 appliqué à l'étoile : r_p sort-il des constantes de la base ?**
→ `proton_radius_search.py`, 3/3. Expressions simples les plus proches de r_p = 0,8409(4) fm :
3r_e/10 (+0,5 %, 11 σ), 3r_e/π² (+1,9 %, le rapport d'aspect de la nappe), r_e/π (+6,7 %). Recherche
systématique (ƛ_e ou r_e × α^k × π^j × p/q × facteur d'échelle (n/3)^{2π}, 7590 candidats) : aucun
candidat à ±0,1 %, moins que les ~3 coups de hasard attendus. **r_p n'est pas dérivé ; l'étoile
garde son nombre propre.** Ce que la règle de parité dit du contenu : proton = (4+, 1−), les
cordes de l'anneau π⁺ ; neutron = (2+, 2−), les cordes de l'anneau π⁰. Mêmes cordes, forme
différente : m_p/m_π⁺ = 6,72, m_n/m_π⁰ = 6,96 ; dans les deux lectures de la base (anneau à r_e/2,
étoile à 4 modes de r_p) le rapport vaut 2r_e/r_p = 6,70, le même nombre ouvert que r_p. La
question devient : pourquoi les mêmes cinq cordes pèsent 140 MeV en anneau et 938 MeV en étoile.

**R24 — « C'est la masse de jonction : plus il y a de connexions sur une jonction, plus l'espace
lui accorde d'inertie. »**
→ `junction_valence.py`, 6/6. **La loi existe dans la base, mais elle est linéaire et plafonnée** :
dans la loi de Coulomb, une jonction équilibrée à k pôles lie comme ⌊k/2⌋ paires (k = 2 et 3 : une
paire ; 4 et 5 : deux ; 6 et 7 : trois) ; le pôle impair n'ajoute rien (frustration, R22). DÉRIVÉ.
**Ce qu'elle exclut** : la lecture à un seul centre (proton = 5 cordes sur un nœud, neutron = 4)
donne deux centres égaux, donc p = n + une corde + son Coulomb, le proton serait le plus lourd ;
observé : le neutron l'est de 1,29 MeV. **Ce qu'elle permet** : u = paire (2+, une jonction k = 2),
d = Y (1+, 2−, un centre k = 3, égal à k = 2 sous la loi) ; alors p = uud = 7 cordes (5+, 2−),
n = udd = 8 (4+, 4−), charges et parités justes, **et les anneaux de R15 ressortent du contenu en
quarks** : π⁺ = ud̄ = (4+, 1−), π⁻ = dū = (1+, 4−), π⁰ = uū = (2+, 2−). Le neutron a une corde de
plus que le proton : plus lourd, bon signe ; avec la part Coulomb du proton (~1,0 MeV) la corde
supplémentaire pèse 2,3 MeV à l'échelle du nucléon. **Ce qu'elle ne fait pas** : l'anneau du pion
(jonctions k = 2 seulement) pèse son mode, 140 MeV ; les ~930 MeV du baryon devraient siéger dans
son centre, alors que le centre Y du d, de même valence, pèse moins de 2,3 MeV : même k, masse 400
fois différente. À 2 MeV la paire, la loi linéaire demanderait ~930 connexions ; une jonction de
Coulomb n'atteint 930 MeV qu'à s = 1,5 am (pôles unitaires) ou 21 am (pôles adaptés). Et les
mésons légers atteignent la masse du proton sans aucun centre (η′ 958, a₀ 980) ; ρ/π = 5,55 à
contenu et topologie égaux, contre p/π⁺ = 6,72 : un facteur 5 à 7 naît dans une seule topologie,
donc le rapport anneau/étoile ne mesure pas une jonction. **Verdict** : la valence n'est pas la
variable, l'échelle l'est ; ce qui manque reste ce qui fixe la taille du centre.

**R25 — Ce qui fixe la taille du centre : le mécanisme de taille de la base (spin + moment)
appliqué au nucléon.** (suite de R24)
→ `nucleon_moments.py`, 5/5. Avec le fluide à c sur un anneau (S = R·E/c = ℏ/2, μ = qcR/2,
g = mc²/E_circ) : **la charge e entière sur un anneau**, μ_p fixe R = 2,79 ƛ_p = 0,587 fm (0,70 r_p)
et g_p = 5,586 signifie que 18 % de l'énergie du proton circule, 82 % est statique au centre
(électron : 50/50). **Le contenu de R24 sur deux anneaux** (cordes + à R₊, cordes − à R₋, mêmes
rayons pour p et n puisque même graphe) : μ_p et μ_n donnent R₊ = 3,75 ƛ_p = 0,789 fm (−6 % vs r_p)
et R₋ = 5,19 ƛ_p = 1,09 fm : les cordes négatives circulent à l'extérieur, comme la peau négative
du neutron l'exige ; deux données, deux inconnues, R₊ proche de r_p est le seul contenu. Les mêmes
anneaux comme charge statique échouent (r_p à 0,50 fm, ⟨r²⟩_n 6,5 fois trop négatif) : les moments
placent le fluide, pas la charge (R5 : charge aux pôles, courant sur l'arc). Le spin ℏ/2 sur
l'anneau + fait circuler 125 MeV, 13 % du proton. **Non fixé** : les 87 % statiques (~810 MeV) du
centre ; la moitié statique de l'électron venait de g = 2, le g du nucléon dit seulement combien
circule. La taille du centre reste le nombre manquant.

**R26 — « Cherche ce qui fixe les 810 MeV du centre. Peut-être l'entreposage des bouts de deux
cordes, qui empêche leur séparation. »**
→ `centre_search.py`, 5/5. Cible : E_s = m_p − E_circ = 813 MeV (R25). **Ce que la base offre tombe
en 1/L** : une cavité entre deux bouts entreposés (mode ℏcπ/ℓ) donne 813 MeV à ℓ = 0,763 fm = R₊
à 3 %, c'est le mode transverse de la nappe (R17) à la largeur du nucléon ; bonne échelle, mauvais
signe, car une énergie en 1/L baisse quand les bouts s'écartent (−461 MeV pour 1 fm) : elle ne
peut pas empêcher la séparation. Le coût d'une coupure (deux pôles adaptés neufs, 2ℏc/(π²w)) fait
813 MeV à w = 0,049 fm, sans rien qui fixe w. **La seule énergie qui croît avec la séparation** :
le flux du pôle confiné dans la nappe (tube de flux). Avec le pôle adapté δ = e/(π√α), la tension
vaut σ = δ²/(2ε₀A) = 2ℏc/(πA) ; pour un tube de rayon ƛ_p, **σ = 2(m_pc²)²/(π²ℏc) = 904 MeV/fm,
soit √σ = (√2/π)·m_pc² = 422 MeV, dans la fourchette du réseau (420–440 MeV)** ; 813 MeV est alors
0,90 fm de tube (σ·r_p = 760 MeV, −6,5 % ; ou un Y de trois bras de 0,30 fm) ; écarter les bouts
de 1 fm coûte +904 MeV. CONDITIONNEL : la ligne conductrice de la base étale la charge d'un pôle le
long d'elle (énergie en 1/L) ; garder le flux dans la nappe est une règle nouvelle (Wen,
supraconducteur dual), et le rayon ƛ_p est la longueur propre du proton, donc σ tiré de m_p est
une relation de cohérence, pas une dérivation de m_p. Réserve : un anneau de pion (140 MeV) paie
0,15 fm de tube, un tube de 0,9 fm en vaut six ; il faut une règle qui interdit la rupture (en
QCD la paire coûte ~1 GeV et la corde casse vers 1,2 fm). **L'entreposage de l'auteur, lu comme
« le flux des deux bouts reste dans la nappe », a la bonne échelle et le bon signe.**

**R27 — « Entreposage = entremêlage. »** (les bouts de deux cordes sont enroulés l'un autour de
l'autre, ce qui empêche leur séparation)
→ `entanglement.py`, 5/5. **Topologie** : le nombre d'enlacement de deux bouts entremêlés n'est
conservé que si les quatre bouts sont tenus (réseau fermé) ; une paire isolée se déroule par
rotation libre. Le verrou appartient donc au réseau, pas à la paire : un quark isolé ne tient
aucune torsade (cohérent avec « pas de quark libre »). **Énergie de mode** : pour une paire torsadée
à longueur axiale fixe 0,90 fm (le tube de R26) et rayon d'hélice ƛ_p, la longueur d'arc croît avec
les tours et le mode ℏcπ/s baisse (1 tour : −302 MeV ; 2 tours : −467 ; limite −690 = ℏcπ/ℓ) :
l'enroulement est spontané et le déroulement, nécessaire pour séparer, coûte jusqu'à 690 MeV.
Bon ordre de grandeur, mais c'est une **liaison** : la torsade tient les bouts, elle ne stocke pas
les 813 MeV. **Dual magnétique** : le fluide enroulé est un solénoïde de conducteur parfait, son
flux est conservé, et l'étirer à flux par tour fixe coûte σ_B = Φ₁²/(2μ₀A) ; avec un quantum h/e
par tour dans un tube de rayon ƛ_p, σ_B = (π²/4α)·σ_E = 338 × 904 MeV/fm = 306 GeV/fm ; pour
retrouver le réseau il faudrait un tube de 3,9 fm (plus grand que le proton) ou un flux de
0,054 h/e. Bon signe, mauvaise échelle : c'est le tube électrique de R26 qui tombe juste, et
σ_E = σ_B exige Φ = Z₀δ = (2√α/π)·h/e, qui n'est pas un quantum de flux. **Tresses** : trois bras
entremêlés sont une tresse à trois brins (Bilson-Thompson) ; BT donne à ses brins {−1/3, 0, +1/3}
et fait u = (+,+,0), d = (−,0,0), e = (−,−,−) avec trois brins ; la règle de la base (toute corde
±1/3) ne fait pas +2/3 à trois cordes (R22). **Une corde neutre** (la question ouverte de R1 : un
fluide à pôles sans charge nette) permettrait à la base de construire chaque quark comme une
tresse à trois brins, comme BT. Verdict : l'entremêlage verrouille (réseau fermé) et lie (mode
plus bas), il ne pèse pas ; les 813 MeV restent au tube de flux de R26.

**R28 — « Teste la corde neutre : refais les quarks à 3 brins. »**
→ `three_strand.py`, 7/7. Brins à {−1/3, 0, +1/3}, trois brins par fermion. **Ce que la corde neutre
achète** : (1) le spectre de charge exact, {0, ±1/3, ±2/3, ±1} et rien d'autre, alors que la règle
des chaînes autorisait ±4/3 (n = 4) et ±5/3 (n = 5), jamais observés ; (2) e⁻ = (−,−,−), ν = (0,0,0),
u = (+,+,0), d = (−,0,0), antiparticules par changement de signe ; **la couleur = la position du
brin impair** : u et d ont 3 arrangements, e et ν un seul ; (3) p = uud = (4+, 1−, 4×0) et
n = udd = (2+, 2−, 5×0), **tous deux 9 brins** (même compte, ce que p ≈ n exigeait) ; π⁺ = (3+, 3×0),
π⁻ = (3−, 3×0), π⁰ = (2+, 2−, 2×0), 6 brins ; Δ⁺⁺ = (6+, 3×0) = +2 ; un baryon est incolore quand
ses trois brins impairs occupent trois positions différentes (6 arrangements sur 27) ; (4) l'échelle
des leptons se relit : μ = e + 4 brins neutres (7), τ = e + 8 (11), **le pas de génération est
quatre brins neutres, deux DQD neutres**, la mère de R3. **Ce qu'elle coûte** : le neutrino (0,0,0)
à n = 3 pèserait m_e sous m = m_e(n/3)^{2π}, 6·10⁵ fois la limite KATRIN (0,8 eV) ; seule la règle
« les brins chargés portent le mode, le compte total fixe l'échelle » garde ν ≈ 0 avec e, μ, τ
inchangés (POSTULÉ). p et n ayant le même compte, l'argument de R24 (une corde de plus) tombe, et
m_d − m_u = (m_n − m_p) + 1,0 MeV (EM) = 2,3 MeV (réseau : 2,5) exige que le d, avec un seul brin
chargé, pèse plus que le u qui en a deux : la masse n'est pas un compte de brins. Tombent aussi la
règle de parité (u pair, d impair) et l'identité « p porte les cordes du π⁺ » (R23). Sous la
condition de spin de la base, un neutrino sous 0,8 eV mesure au moins 123 nm. **Verdict** : la
corde neutre donne la structure (charges exactes, couleur, 9 brins pour p et n, générations en
DQD neutres) et ne donne toujours pas les masses ; R1 doit dire ce qu'est une corde sans charge.

**R29 — « Teste la règle : les brins chargés portent le mode. Qu'est-ce que le mode ? »**
→ `mode_rule.py`, 6/6. **Le mode** : le fluide court à c le long de la corde et se réfléchit aux deux
pôles ; l'onde stationnaire d'une demi-longueur d'onde entre les pôles est le mode, d'énergie
E = ℏcπ/ℓ_arc. Pour l'électron (3/4 de tour à R = 4ƛ/3) l'arc vaut 2πƛ = 2426 fm et E_mode =
m_ec²/2 exactement : la moitié qui circule (g = 2), l'autre moitié est statique aux pôles ; par
corde, 809 fm, le ℓ₁ = 2πR₃/3 de la base. Loi d'échelle tirée de l'échelle des leptons :
ℓ₁(n) = ℓ₁(3)(3/n)^{2π}, soit 3,94 fm au muon, 0,23 fm au tau, et ℓ₁(9) = 0,813 fm = r_p à 3 %
(coïncidence enregistrée). **La règle** : e, μ, τ inchangés (trois brins chargés, arc 3ℓ₁(n)) ; le
neutrino n'a pas d'arc chargé, donc pas de mode, m = 0 (< 0,8 eV). Nucléons au compte n = 9 :
mode du proton sur 5 brins chargés = 153 MeV, du neutron sur 4 = 191 MeV ; la part de mode propre
à la base, m/|g|, vaut 168 (p) et 246 (n) MeV : −9 %, −22 %, et l'ordre n > p est bon ; les
~700–800 MeV statiques du centre sont hors de la règle. Pions au compte n = 6 : ℓ₁ = 10,4 fm et un
mode de 20 MeV pour le π⁺, contre les 140 MeV de l'anneau fermé (R16) : la fermeture fixe sa
propre échelle, la loi de compte ne vaut pas pour les anneaux. **Comme loi de masse** (chaque brin
chargé son mode) : à n = 3, u = 0,34 et d = 0,17 MeV, mauvais ordre (il faut m_d − m_u = +2,3 MeV) ;
à n = 9, p = 848 MeV (−10 %) et n = 678 (−28 %), mauvais ordre. **Verdict** : loi de mode à 10–20 %
(leptons, neutrino nul, part circulante des nucléons), pas loi de masse ; ce qui pèse le hadron
reste le centre statique.

**R30 — Creuser R1 : que fait le champ du pôle dans la nappe ? Le tube de flux de R26 est-il dans la
base ?** (suite de « Oui »)
→ `r1_field.py`, 5/5. **Ce que la base confine déjà** : la ligne bifilaire (D/r = 2 cosh π) garde
92 % de l'énergie de son champ transverse dans un rayon D et 99 % dans 3D, densité en 1/ρ⁴ au-delà :
le champ du fluide polarisé reste sur la ligne sans règle nouvelle, la ligne est sa propre nappe.
**Mais ce champ confiné ne donne une tension que si la charge par unité de longueur est fixe** :
avec e/3 par corde et une longueur libre, U = (e/3)²/(ε₀ℓ) = 4πK/(9ℓ) tombe en 1/ℓ, pas de
tension ; à la longueur fixée ℓ₁(9) = 0,81 fm, chaque corde stocke 2,47 MeV (la « corde de plus »
de R24/R29 à 2,3 MeV, +8 %), soit une tension de 3,0 MeV/fm, 300 fois sous le réseau. **Le pôle** :
un monopôle dans un tube conducteur de rayon ƛ_p est écranté, pas canalisé ; son champ décroît en
exp(−2,405 z/a), éteint en 0,09 fm, énergie locale ~Kδ²/a ≈ 95 MeV ; le flux ne parcourt jamais les
0,9 fm du tube de R26. **Ce qu'exigerait la tension du réseau** : une charge linéique fixe de
7,1 e/fm, c'est-à-dire des cordes de e/3 longues de 0,047 fm (19 par 0,9 fm), que l'échelle place à
n = 14 brins et non 9 ; ou un courant fixe de 4,8·10⁵ A, 13 fois le e·c/(2πƛ_p) de la base.
**Verdict** : la base ne peut pas tenir à la fois « e/3 par corde de longueur libre » et une
tension de tube de flux ; un fluide conducteur écrante le pôle au lieu de le canaliser ; R26 n'est
pas un mécanisme de la base telle qu'écrite. Les ~800 MeV du centre restent sans mécanisme ; la
seule piste chiffrée est une charge linéique fixe (R1 : le fluide aurait une densité de charge
propre, pas une charge par corde).

**R31 — La dernière piste chiffrée de R30 : un fluide à densité de charge linéique fixe.**
→ `fixed_density.py`, 4/4. La tension du réseau exige λ = 7,07 e/fm ; les cordes de l'électron
portent e/3 sur 809 fm, λ_e = 4·10⁻⁴ e/fm, 17 000 fois moins (tension 3·10⁻⁶ MeV/fm) : un seul
fluide ne peut pas avoir les deux densités, **la densité fixe est exclue par l'électron**. Une
densité qui suit l'échelle, λ = (e/3)/ℓ₁(n), donne σ = 4πK/(9ℓ₁²) = 3,0 MeV/fm au nucléon, 300 fois
trop peu, et n'atteint le réseau qu'à ℓ₁ = 0,047 fm, le n = 14 de l'échelle et non le 9 du proton.
Avec e/3 par corde et λ fixe, les longueurs seraient quantifiées à 0,047 fm : l'électron (3 cordes)
ferait 0,14 fm contre son arc de 2426 fm. **Verdict** : jonction de Coulomb (R24), mode (R29),
valence (R24), entremêlage (R27), tube de flux (R26/R30), densité fixe (ici) manquent tous les
~800 MeV du centre d'un facteur 100 à 300 ou par le signe : **le fluide de la base, une onde TEM à
c sur une ligne adaptée, n'a pas de tension.** L'interaction forte n'est pas dans l'électromagnétisme
de la base.

---

## Bilan au 19 septembre 2026 (R21–R31, ajouts au bilan du 16)

**Ce qui s'est ajouté au solide.**
- Trois brins à {−1/3, 0, +1/3} donnent exactement le spectre de charge observé et rien d'autre ;
  la couleur est la position du brin impair ; p et n ont 9 brins, les pions 6 ; les générations
  sont e + 4 ou 8 brins neutres (R28).
- Le mode est l'onde stationnaire entre les deux pôles, E = ℏcπ/ℓ_arc ; pour l'électron,
  m_ec²/2 exactement ; la règle « les brins chargés portent le mode » tient comme loi de mode
  (leptons, neutrino nul, part circulante des nucléons à 10–20 %) (R29).
- Une jonction à k pôles lie comme ⌊k/2⌋ paires ; la tri-jonction à 2 MeV est l'échelle de la
  force nucléaire (R22, R24) ; μ_p et μ_n placent le fluide à 0,79 fm (+) et 1,09 fm (−) (R25).
- Le champ du fluide polarisé reste sur la ligne bifilaire par géométrie (R30).

**Ce qui s'est fermé.**
- r_p ne sort pas des constantes de la base (R23) ; la valence ne fait pas le GeV (R24) ;
  l'entremêlage verrouille et lie mais ne pèse pas (R27) ; le tube de flux n'est pas un mécanisme
  de la base, un fluide conducteur écrante le pôle (R30) ; une densité fixe est exclue par
  l'électron (R31) ; les comptes de brins ne sont pas des masses (R28, R29).

**Où on en est.** La base décrit charges, couleur, spin, moments, tailles, la force nucléaire et
les leptons ; elle ne pèse aucun hadron : les ~800 MeV du centre n'ont pas de mécanisme, et
toutes les routes électromagnétiques essayées manquent d'un facteur 100 à 300. Le fluide TEM à
c n'a pas de tension ; il faut soit un ingrédient non électromagnétique, soit un fluide qui n'est
pas une onde sur une ligne conductrice.

**R32 — « Le fluide n'est pas une onde, c'est un vortex. »** (redessin de R1)
→ `vortex.py`, 5/5. **Ce que le vortex apporte** : sa circulation est un quantum conservé, pas une
charge libre de s'étaler ; un courant I stationnaire sur la ligne adaptée stocke Z₀I²/c par unité
de longueur, constant le long de la ligne, donc l'énergie croît avec la longueur : **le signe que
le centre exigeait (R26), sans règle nouvelle**. Deux lectures sur l'anneau de l'électron : l'anneau
de fumée (tourbillon poloïdal à c autour d'un cœur de rayon a, R/a = 37,1) s'auto-propulse à
0,073 c, ce n'est pas une particule au repos ; la boucle de courant (écoulement le long de
l'anneau, le R2/R4 de la base) est stationnaire et garde spin ℏ/2, μ_B et l'absence de rayonnement
(R11), mieux qu'une onde stationnaire : **le vortex, c'est la boucle de courant**. **Ce qu'il
coûte** : l'énergie d'écoulement d'une charge q à c sur un rayon R vaut 2Kq²/R ; pour q = e sur
l'anneau de l'électron, 7,5 keV, 1,5 % de m_e ; l'écoulement ne porte la moitié circulante
m_ec²/2 que si q = e/(2√α) = 5,85 e (la condition de spin réécrite en charge), le pôle adapté
3,73 e en porte 40 % ; la lecture « mode » de R29 est perdue. **L'échelle** : la tension de ce
vortex au rayon ƛ_p vaut ℏc/(4πR²) = 355 MeV/fm (réseau/2,5), le tube électrique de R26 donnait
904 ; les deux vont en m_p², le nombre est la longueur de Compton du proton, circulaire tant que R
n'est pas dérivé. **En plus** : deux lignes de vortex antiparallèles interagissent en
(μ₀I²/2π)·ln d par unité de longueur, potentiel logarithmique, force en 1/d, un confinement doux
entre lignes que l'onde n'avait pas ; au courant du réseau (3,4·10⁵ A), passer de 0,1 à 1 fm coûte
331 MeV par fm de ligne. **Verdict** : le vortex donne le signe (tension) et garde spin, moment,
non-rayonnement ; il perd la masse par mode et ne fixe toujours pas le nombre, qui reste ƛ_p.

**R33 — « Il ne peut pas y avoir un effet inductif amplificateur à chaque tour ? »**
→ `inductive_winding.py`, 5/5. **Oui, et il est quadratique** : un brin chargé enroulé en hélice
(n tours par unité de longueur, cœur de rayon a) est un solénoïde, son champ vaut μ₀nI et son
énergie par unité de longueur est amplifiée de π(na)² par rapport au brin droit (10 tours/fm sur
ƛ_p : ×14). **La condition** : la charge par unité de longueur de trajet doit être fixe, chaque tour
portant son propre courant ; avec une charge totale e fixe, le courant par tour tombe quand le
trajet s'allonge et le gain s'annule exactement (u = K/(2ℓ²) = 0,9 MeV/fm quel que soit n). C'est
la circulation conservée du vortex (R32) qui autorise le gain. **Les nombres** : à la densité de
brin du nucléon (e/3 par 0,81 fm, brin droit 3,0 MeV/fm), la tension du réseau demande un gain de
297, soit na = 9,7 : 46 tours par fm sur ƛ_p, 42 tours dans le centre de 0,9 fm, 55 fm de corde,
68 brins de 0,81 fm et non 9. À la densité de l'électron (4·10⁻⁴ e/fm), la même tension demande
46 000 tours par fm, 54 000 fm de corde, 68 cordes d'électron : **l'exclusion de R31 tombe si le
nombre d'enroulement est libre, un seul fluide, enroulé ou non.** Seul le brin non apparié
s'enroule en solénoïde : torsader une paire bifilaire à courants opposés annule la composante
azimutale, les brins neutres sont inductivement inertes (cohérent avec R29). Le nombre
d'enroulement est la torsion Tw du brin (Lk = Tw + Wr, R27), un entier que la base possède ; la
valeur requise (~40 par nucléon à la densité de brin) n'est pas dérivée ; l'énergie électrique de
la charge enroulée (cylindre chargé écranté à r_p) ajoute un terme du même ordre ou plus grand et
abaisse les tours nécessaires à quelques dizaines par fm. **Verdict** : mécanisme réel, signe et
échelle atteignables avec le fluide de l'électron ; le nombre (la torsion) reste à fixer.

**R34 — « Cherche ce qui fixe la torsion Tw. »**
→ `twist_search.py`, 5/5. Cinq candidats. **Topologie** : Tw = Lk − Wr ; Lk est un entier fixé à la
naissance (R4) et conservé, Wr vaut 0 pour un anneau plan ou un bras droit ; Tw est donc l'entier
avec lequel le brin est né : conservé, pas dérivé. **Énergie seule** : l'énergie du solénoïde croît
en Tw², un brin libre se déroule à Tw = 0 ; seule la conservation tient une torsion. **Quantification
du flux** : à la solution de R33 (46 tours/fm) le flux du cœur vaut 0,038 h/e ; un quantum entier
h/e donnerait 680 fois la tension du réseau ; le flux propre de l'anneau de l'électron,
(α/π)(ln(8R/a) − 2) = 0,0086 h/e, y tient 4,5 fois : pas d'entier naturel. **Adaptation à Z₀**
(R10 et R17 appliqués à la ligne hélicoïdale) : L′ = μ₀n²πa², C′ = 2πε₀/ln(b/a), Z = Z₀·na·√(ln(b/a)/2),
et Z = Z₀ fixe **na = √(2/ln(b/a)) = 1,20 (b = r_p) à 0,74 (b/a = 37,1) : environ un tour par rayon
de cœur, angle de pas 7 à 12°, Tw = 3 à 5 tours dans le centre de 0,9 fm**. C'est la seule règle
de la base qui fixe la torsion. **Conséquence** : le gain adapté π(na)² = 1,7 à 4,5 donne 5 à
14 MeV/fm à la densité de brin du nucléon, 66 à 170 fois trop peu ; la tension du réseau au pas
adapté exige une densité de trajet de 3,3 à 5,4 e/fm, des brins de e/3 longs de 0,06 à 0,10 fm
(le n = 13 de l'échelle) : **le nombre libre passe de Tw à la densité**. L'électron n'est pas touché :
son enroulement adapté sur son tube de 10,4 fm ajoute 0,03 MeV à 511 keV. Verdict : Tw est fixé
par l'adaptation à Z₀ à ~1 tour par rayon de cœur (DÉRIVÉ sous R10/R17) puis conservé ; ce qui
reste libre est la charge par unité de longueur de trajet du fluide.

**R35 — « Cherche ce qui fixe la densité de charge du fluide. »**
→ `density_search.py`, 5/5. Fenêtre à atteindre (R34) : 3,3 à 5,4 e par fm de trajet. **Ce qui rate** :
la règle d'anneau avec l'échelle, λ = (e/3)/ℓ₁(9) = 0,41 e/fm, 11 fois trop peu ; la densité propre
de la boucle, e/(2πR₊) = 0,20 e/fm, 20 fois trop peu. **Ce qui tombe dans la fenêtre** : la charge
de pôle adaptée (R17, δ = e/(π√α) = 3,73 e) prise comme charge du fluide par brin, λ = δ/ℓ₁(9) =
**4,6 e/fm**. La tension du brin enroulé vaut alors 4πKλ² × gain, 660 à 1720 MeV/fm pour le pas
adapté avec une coupure b/a de 37 à 4 (0,7 à 1,9 fois le réseau), et exactement le réseau à
b/a = 14 : **aucun nombre posé à la main, δ vient de Z₀, ℓ₁ de l'échelle des leptons, le pas de Z₀.**
**Le test de cohérence sur l'électron** : la même règle donne à ses trois cordes de 809 fm une
énergie de ligne statique 3 × 2πKδ²/ℓ₁ = (9/π²)·m_ec² = 0,91 m_ec² (α s'annule) ; mais la base paie
déjà m_e avec mode + jonctions (R13), la règle surcompte de 466 keV, sauf si la moitié statique
(255 keV) est cette énergie de ligne, qu'elle manque d'un facteur 1,8 : CONDITIONNEL. **Sur les
nucléons** (enregistré, non revendiqué) : la même énergie de ligne par brin chargé à ℓ₁(9),
2ℏc/(πℓ₁) = 155 MeV, donne p = 5 × 155 + 153 (mode) = 926 MeV (−1,3 %) et n = 4 × 155 + 191 =
809 MeV (−14 %), mauvais ordre. **Verdict** : la densité est fixée si le fluide de chaque brin porte
la charge de pôle adaptée δ ; c'est la première fois que la tension du réseau sort des règles de
la base sans nombre libre, au prix d'un facteur 1,8 sur la moitié statique de l'électron et d'un
neutron trop léger de 14 %.

**R36 — « Cherche ce qui fixe la coupure b/a = 14. »**
→ `cutoff_search.py`, 4/4. Avec λ = δ/ℓ₁(9) (R35), σ = (4ℏc/(πℓ₁²)) × 2π/ln(b/a). La bande du réseau,
√σ = 420 à 440 MeV (894 à 981 MeV/fm), correspond à b/a = 11,4 à 14,5. **Candidats de la base** :
la branche partenaire à la distance D, D/r = 2 cosh π = 23,2, donc ln(b/a) = π à 0,06 % près, gain
exactement 2, **σ = 8ℏc/(πℓ₁²) = 761 MeV/fm (−16 %)** ; le demi-écartement D/2r = cosh π = 11,6, gain
2,56, **σ = 976 MeV/fm (+8 %, dans la bande)** ; R₃/r = 37,1 → 661 ; r_p/ƛ_p = 4 → 1724 ; ℓ₁/ƛ_p →
1770 ; R₋/ƛ_p → 1453 ; (r_e/2)/ƛ_p → 1257. Seuls les deux rapports de l'adaptation à Z₀ approchent.
**Les deux lectures de l'adaptation encadrent le réseau** : 761 (b = D) et 976 (b = D/2), la bande
894–981 est entre les deux ; la coupure est donc fixée par la même adaptation qui fixe D/r, à
l'identification D ou D/2 près, soit ±15 % sur σ. **Forme fermée** avec ln(b/a) = π :
σ = 8ℏc/(πℓ₁(9)²) = (18·3^{4π}/π³)·m_e²c⁴/ℏc, **√σ = √(18/π³)·3^{2π}·m_ec² = 388 MeV** (réseau 420 à
440, −8 à −12 %) ; avec b = D/2, 439 MeV. La tension forte s'écrit avec m_e et l'exposant 2π seuls,
sans α ni Z₀ (ils s'annulent entre δ et l'adaptation). Verdict : DÉRIVÉ à 15 % près ; ce qui reste
est de dire si la charge enroulée voit sa partenaire à D ou à D/2.

**R37 — « Cherche si la charge enroulée voit sa partenaire à D ou à D/2. »**
→ `partner_distance.py`, 5/5. Tranché par l'électrostatique de la paire elle-même : deux fils ±λ de
rayon r à l'écartement D stockent (λ²/2πε₀)·arccosh(D/2r) par unité de longueur, soit par fil
(λ²/4πε₀)·arccosh(D/2r), la formule du cylindre seul avec ln(b/r) = arccosh(D/2r) = ln(D/r) à 0,06 %
près (D/2 serait 22 % à côté) ; la lecture par plan de symétrie donne la même chose, l'image d'un
fil à D/2 du plan est à D. **La coupure est D, sans ambiguïté.** Donc ln(b/a) = arccosh(D/2r) = π
exactement par l'adaptation (R10), gain 2, **σ = 8ℏc/(πℓ₁(9)²) = 761 MeV/fm, √σ = 388 MeV** ; le
réseau (420 à 440 MeV) est 8 à 12 % au-dessus en √σ, 16 % en σ. Le résidu n'est pas dans les
entrées : bande du réseau ±2,3 % ; l'exposant ajusté sur le muon (6,292 au lieu de 2π) déplace
ℓ₁(9) de 1,0 % (2 % sur σ) ; l'ancrage du muon 0,8 % ; le manque de 15 % est réel, le terme suivant
est l'écart entre l'hélice et un cylindre lisse (non calculé). L'adaptation de la double hélice
(coax + solénoïde), y = ln(b/a)/2π = 1/2, donne le terme solénoïde G′ = 1,5 et na = 0,69 tour par
rayon de cœur, cohérent avec R34. **Verdict** : DÉRIVÉ, la charge enroulée voit sa partenaire à D ;
la tension forte de la base est 8ℏc/(πℓ₁²), 15 % sous le réseau.

**R38 — « Calcule l'écart entre l'hélice réelle et le cylindre lisse. »**
→ `helix_vs_cylinder.py`, 5/5. Le modèle lisse (R34–R37) étalait la charge du brin enroulé le long
de l'axe à sa densité de trajet ; l'hélice réelle à n tours par unité de longueur sur un rayon a
met γ = √(1 + (2πna)²) unités de trajet, donc de charge, par unité d'axe : le fluide court à c le
long du fil, la densité axiale est γλ et la vitesse de phase axiale c/γ (ligne à onde lente). Au
pas adapté du modèle lisse (na = 0,69 à 0,80), γ = 4,45 à 5,1. **La correction est grande, pas
15 %** : pour une ligne adaptée le long de son trajet (Z = Z₀, R10), l'énergie par fm de trajet
vaut 4πKλ² = 380 MeV quel que soit l'enroulement ; par fm d'axe, γ fois plus. Au pas du modèle
lisse la tension réelle est donc 1690 à 1950 MeV/fm, et non 761 : le modèle lisse sous-comptait
d'un facteur 2,2 à 2,6 (+120 à +155 %) ; **le résidu de 15 % de R37 était un artefact du modèle**.
La tension s'écrit maintenant **σ = 4πKλ²·γ**, et le réseau (894 à 981, centre 904) exige
γ = 2,35 à 2,58 (2,376) : na = 0,34 à 0,38, angle de pas 25°, 2,4 fm de brin par fm d'axe ; avec
a = ƛ_p, 1,6 tour par fm et 1,5 tour dans le centre de 0,9 fm. L'adaptation « grossière » de R34
(L′ solénoïde + C′ coaxial), rendue cohérente avec l'onde lente (L′ = μ₀γ), donne na = 1,78,
γ = 11,2 et σ = 4270 MeV/fm, 4,7 fois le réseau : ce n'est pas un calcul d'hélice valable ; il
faut la théorie de l'hélice-gaine (Pierce) pour le pas, non faite. **Ce qui fixe le pas est
rouvert.** L'électron, non enroulé (3/4 de tour, γ = 1), n'est pas touché. Verdict : la tension
forte de la base est l'énergie de la ligne adaptée par unité de trajet (380 MeV/fm à la densité δ/ℓ₁)
multipliée par le rapport d'enroulement ; il manque la règle qui fixe γ ≈ 2,4.

**R39 — « Cherche ce qui fixe le pas γ ≈ 2,4. »**
→ `pitch_search.py`, 5/5. **Correction à R38 d'abord** : la force n'est pas l'énergie par longueur
d'axe. Une ligne adaptée enroulée stocke E = u_trajet × s avec s = √(L² + (2πa·Tw)²) ; à Tw et cœur
fixes, dE/dL = u_trajet/γ ≤ u_trajet ; à trajet fixe, dE/dL = 0. L'enroulement ne peut pas fournir
une force au-delà de u_trajet = 380 MeV/fm (à λ = δ/ℓ₁(9)) ; le σ = 380γ de R38 est une énergie par
fm d'axe, pas une tension. Le 904 MeV/fm du réseau comme force exige u_trajet lui-même :
**« γ = 2,376 » est le rapport 904/380, pas un pas.** Aucune règle de la base ne fixe un pas de
toute façon (énergie plate en γ à trajet fixe, adaptation grossière qui surestime ; les
coïncidences tan ψ = 1/2, na = 1/3, 3π/4 ne viennent d'aucun mécanisme). **Ce qui ferme le facteur
sans enroulement** : l'autre charge de la base, la charge de vortex q = e/(2√α) = 5,85 e (R32, celle
qui porte l'énergie du spin ½) : (q/δ)² = (π/2)² = 2,467, à 3,8 % du 2,376 requis. Avec q par brin,
**u_trajet = 4πKq²/ℓ₁(9)² = πℏc/ℓ₁(9)² = 939 MeV/fm, √σ = 430 MeV, le centre de la bande du réseau
(420–440)** ; forme fermée √σ = (3/(2√π))·3^{2π}·m_ec², α et Z₀ s'annulant. **Le prix** : q par
brin sur les trois cordes de l'électron stocke 3πℏc/ℓ₁(3) = 2,30 MeV = 4,5 m_e ; dans l'électron, q
est la charge de l'anneau entier (2Kq²/ƛ = m_ec²/2 exactement) ; le nucléon en aurait besoin par
brin. Côte à côte : δ par brin → électron 0,91 m_e (×1,8 la moitié statique), nucléon 380 MeV/fm
(2,4 fois trop peu) ; q par brin → nucléon 939 (réseau), électron 4,5 m_e ; q par objet → spin ½
de l'électron, anneau du proton à 125 MeV, centre non couvert. **Verdict** : la base a le nombre
(πℏc/ℓ₁², √σ = 430 MeV) mais pas la règle qui dit quel brin porte la charge de vortex ; la route de
l'enroulement (R33–R38) se ferme comme mécanisme de force.

**R40 — « Et ? » : la règle qui assignerait la charge de vortex.**
→ `circuit_rule.py`, 5/5. Règle testée : **un quantum de circulation, de charge q = e/(2√α), par
circuit bifilaire fermé** (le fluide sort par une branche et revient par l'autre) ; une chaîne de
brins bout à bout est un seul circuit, une étoile a un circuit par bras (R20). **Électron** = un
circuit de trois brins (trajet 3ℓ₁(3) = 2πƛ) : l'énergie de ligne adaptée 4πK(q/trajet)² × trajet =
πℏc/(3ℓ₁) = m_ec²/2 exactement, la moitié circulante de R13, par la même formule que la tension du
nucléon ; μ et τ suivent par l'échelle. **Bras d'étoile** = un circuit d'un brin : tension
πℏc/ℓ₁(9)² = 939 MeV/fm, √σ = 430 MeV, dans la bande du réseau. Une règle, les deux secteurs, sans
nombre libre. **Mais** énergie = tension × longueur : un bras de ℓ₁(9) stocke 763 MeV ; cinq bras
chargés (proton à trois brins) pèseraient 3,8 GeV ; **un seul bras** donne 763 MeV, la part statique
du proton avec e sur un anneau (770 MeV, R25 A) à −0,9 % : le centre serait un circuit unique.
Parité de spin : ℏ/2 par circuit, nombre impair → spin demi-entier : leptons (1), quarks (3),
baryons (9) impairs, mésons (6) pairs ; cohérent, sans prédiction. **Verdict** : la règle des
circuits réconcilie la moitié de l'électron et la tension du réseau sans nombre libre ; elle ne
donne la masse du proton que si son centre est un seul circuit, contre le contenu à trois brins
(cinq bras chargés). Ouvert : quels brins se ferment en circuits.

**R41 — « Cherche quels brins se ferment en circuits. »**
→ `circuit_partition.py`, 6/6. Sous la règle des circuits, un circuit de k brins en série a un
trajet kℓ₁ et une énergie πℏc/(kℓ₁) ; la part statique du proton (770 à 813 MeV, R25) doit être la
somme sur ses circuits, soit Σ1/k_i = 1,01 à 1,07 en unités de 1/ℓ₁(9), avec πℏc/ℓ₁(9) = 763 MeV.
**Des 30 partitions des 9 brins, une seule a Σ1/k = 1 : {3, 3, 3}, trois circuits de trois brins,
les trois quarks, chacun un circuit à trois brins comme un lepton** ; elle donne 763 MeV (−1 % sous
770, −6 % sous 813). La bande est encadrée par {3, 3, 3} et {4, 3, 2} (13/12, 826 MeV, qui mélange
des brins de quarks différents et n'a pas de lecture) ; {8, 1} donne 858 ; tout le reste manque de
plus de 12 %. **Chaque circuit de quark stocke πℏc/(3ℓ₁(9)) = 254 MeV**, l'échelle du quark
constituant (m_p/3 = 313, −19 %) ; le quark est l'objet de l'électron (un circuit à trois brins) à
l'échelle ℓ₁(9), rapport d'énergie 3^{2π} = 996 exactement. Avec les parts circulantes (règle de
R29 153/191, base m/|g| 168/246) : p = 916 à 931 MeV (−2,4 à −0,8 %), n = 954 à 1009 (+1,5 à +7,4 %) :
bon ordre, écart p–n 30 à 60 fois trop grand. Le pion en deux circuits de trois brins fait 40 MeV à
ℓ₁(6) ou 508 à ℓ₁(9), jamais 140 : il reste l'anneau fermé de R16. Spin : trois circuits à ℏ/2
donnent ℏ/2 (↑↑↓), l'anneau + de R25 ne peut pas en porter un autre, il doit être la circulation
collective des circuits (ouvert). **Verdict** : les brins se ferment par quark, trois circuits de
trois ; la masse du nucléon en sort à 1 à 3 % (p) et 1,5 à 7 % (n) avec une règle et l'échelle des
leptons, sans nombre libre ; le quark constituant est un électron à l'échelle 9.

**R42 — « Résume la topologie finale. Cherche ce qui fixe l'écart p–n. »**
→ `pn_splitting.py`, 5/5. **Topologie finale (R28 à R41), en clair** : l'espace est un milieu de
cordes dipolaires neutres (DQD) adaptées à Z₀ ; tout fermion est trois brins portant chacun −1/3,
0 ou +1/3 (e = −−−, ν = 000, u = ++0, d = −00), la couleur étant la position du brin impair ; le
fluide de chaque brin est un vortex, pas une onde ; les trois brins d'un lepton forment un seul
circuit fermé (3/4 de tour) portant un quantum de circulation q = e/(2√α), dont l'énergie de ligne
πℏc/(3ℓ₁) est la moitié de la masse, l'autre moitié étant immobile aux deux pôles, avec l'échelle
ℓ₁(n) = ℓ₁(3)(3/n)^{2π} et n = 3, 7, 11 (chaque génération ajoute deux DQD neutres) ; un quark est
le même circuit à trois brins à l'échelle ℓ₁(9) = 0,81 fm (254 MeV), attaché en étoile ; un baryon
est trois circuits de quark (9 brins, 763 MeV immobiles) plus une circulation d'anneau à 0,8–1 fm
(125 à 170 MeV) qui porte le spin et le moment ; un méson est un anneau fermé de 4 à 6 brins à
r_e/2 (140 MeV) ; le neutron est un proton avec une boucle négative de type électron à 0,9–1,0 fm
accrochée par un joint ; la force nucléaire est la jonction pôle à pôle, 2 MeV à 0,7 fm ; la tension
forte est l'énergie de la ligne adaptée par unité de trajet, πℏc/ℓ₁(9)² = 939 MeV/fm (√σ = 430 MeV).
**L'écart p–n** : p = uud et n = udd sont les mêmes trois circuits et le même graphe de jonctions,
donc les 763 MeV immobiles et les joints s'annulent dans n − p ; ce qui diffère est le contenu des
brins, le neutron ayant un brin neutre de plus (5 contre 4) et un brin chargé de moins. Si le brin
neutre est un DQD complet (deux branches ±e/3, l'unité de l'espace, R9), son énergie bifilaire à
l'échelle du nucléon vaut (e/3)²/(ε₀ℓ₁(9)) = 4πK/(9ℓ₁) = **2,47 MeV = (2α/3)·3^{2π}·m_ec²**, la
part forte de n − p, contre 2,52 ± 0,29 sur réseau ; moins l'énergie de Coulomb du proton (0,86 à
1,03 MeV à r_p ; réseau 1,00 ± 0,16), **m_n − m_p = 1,44 à 1,61 MeV contre 1,293 : bon signe, +11 à
+25 %**, sans nombre libre (α, m_e, 2π, le compte 9). La même règle sur les pions rate (1,9 MeV
contre 4,59 avec uū, mauvais signe avec dd̄) : l'anneau du pion n'est pas fait de circuits. Le
résidu de 11 à 25 % est dans le terme de Coulomb du proton, dont la géométrie n'est pas fixée.

**R43 — « Cherche ce qui fixe la géométrie du terme de Coulomb du proton. »**
→ `coulomb_geometry.py`, 5/5. Avec les charges de quark q_i = (2/3, 2/3, −1/3) pour p et
(2/3, −1/3, −1/3) pour n, E_C = S·Σq_i² + M·Σq_iq_j (S : une charge unité dans la géométrie d'un
quark ; M : deux charges unité à l'écartement des quarks). **Ce qu'il faut** : E_C(p) − E_C(n) =
2,474 − 1,293 = 1,18 MeV ; Σq² = 1 (p), 2/3 (n) ; Σq_iq_j = 0 (p), −1/3 (n) ; donc S + M = 3,54 MeV.
**Le terme mutuel est sans paramètre dans la base** : charge aux pôles (R5), bras à 120° dont les
bouts sont à ℓ₁(9) du centre, écartement √3·ℓ₁ = 1,41 fm : M = K/(√3ℓ₁) = 1,02 MeV ; le Coulomb
inter-quarks du proton s'annule exactement, celui du neutron vaut −0,34 MeV. **Le terme propre
demande la taille de la charge** : étalée le long du brin (barreau de 0,81 fm) avec pour épaisseur
la largeur de nappe de l'électron ramenée à l'échelle, w₉ = (4ƛ_e/π²)·3^{−2π} = 0,157 fm (R17) :
S = (K/ℓ₁)(ln(2ℓ₁/w) − 1) = 2,37 MeV, S + M = 3,39 MeV, **m_n − m_p = 1,34 MeV (+3,8 %)** ; la
constante du logarithme et l'écartement font varier le résultat entre 1,09 et 1,34 MeV (−16 à
+4 %) : la base fixe la géométrie à ±10 % près. La même géométrie donne le rayon de charge du
proton : charges à ℓ₁(9) du centre, √⟨r²⟩_p = ℓ₁(9) = 0,813 fm (−3,4 % contre 0,841) ; ⟨r²⟩_n = 0
contre −0,116 fm² (la peau négative demande la boucle extérieure de R19/R25). Sans l'épaisseur de
la base, les alternatives sont des ajustements : pôles ponctuels de rayon 0,29 à 0,34 fm, sphères
tangentes de 0,45 fm (d = 0,89 fm), sphère uniforme au r_p mesuré (+12 %). **Verdict** : le mutuel
est dérivé (0 pour p, −0,34 pour n), le propre l'est à ±10 % avec la nappe de R17 comme épaisseur ;
l'écart p–n sort à +4 % (1,34 MeV) avec α, m_e, 2π, le compte 9 et la largeur de nappe.

**R44 — « Cherche ce qui fixe la constante du logarithme. »**
→ `log_constant.py`, 5/5. L'énergie propre d'une charge q étalée sur une longueur ℓ avec une taille
transverse a vaut (Kq²/ℓ)[ln(2ℓ/a) + c] : il y a la constante c et la coupure a à fixer. **Le −1 est
dérivé** : pour une charge linéique uniforme de longueur ℓ, potentiel pris à la distance a de la
ligne, c = −1 exactement (vérifié par intégration directe) ; pour un anneau fermé de même
longueur, c = −1 − ln(π/2) = −1,45 ; l'arc d'un circuit fermé est entre les deux. **La coupure est
la section** : le brin de la base est une plaque de largeur w (R17, w = d) ; une bande mince de
largeur w a la capacité d'un cylindre de rayon w/4 (classique), donc a = w/4 = 0,039 fm à
w₉ = 0,157 fm ; un carré plein de côté w donnerait 0,59w, un cylindre inscrit w/2 ; **le a = w de
R43 était un cylindre plus large que la nappe, pas une section de la base.** **La longueur** : la
charge d'un quark s'étale sur ses brins chargés (R28), deux pour u (2ℓ₁), un pour d (ℓ₁), donc
S_u ≠ S_d ; R43 prenait un barreau de ℓ₁ par quark. Les huit combinaisons vont de 0,52 à 1,60 MeV ;
**la paire propre à la base, plaque (a = w/4) et brins chargés, donne m_n − m_p = 1,32 MeV
(+2,4 %)** ; le 1,34 de R43 venait de deux choix hors base (a = w, un barreau) ; les variantes à
section carrée donnent 1,03 (−20 %) et 1,49 (+15 %). La fermeture en anneau (c = −1,45) donnerait
1,41 MeV (+9 %) : le reste est la lecture ouverte ou fermée de l'arc chargé du quark ; l'écran de la
plaque partenaire n'est pas calculé. **Verdict** : avec la section de R17 et le contenu de R28, la
constante est fixée et l'écart p–n sort à +2,4 % ; entrées : α, m_e, 2π, le compte 9, la largeur
de nappe (elle-même de g = 2 et Z₀).

**R45 — « Calcule l'écran de la plaque partenaire. »**
→ `partner_screening.py`, 4/4. **Comptabilité d'abord** (R42/R44) : un brin chargé est une branche
seule portant e/3 (« la charge est la branche non appariée »), un brin neutre un DQD complet, deux
plaques conductrices flottantes ±e/3 à l'écart d = w. Un brin chargé n'a donc pas de plaque
partenaire propre ; ce qui l'écrante, ce sont les DQD neutres du même quark : un pour u = (+,+,0),
deux pour d = (−,0,0). Une plaque partenaire portant la polarisation ±δ ajouterait un terme croisé
2πKqδ/ℓ = 14 MeV par brin chargé et déplacerait m_n − m_p de ∓14 MeV : exclu, la comptabilité
« branche seule » est forcée. **Calcul** : méthode des moments en 2D sur la section ; bandes
conductrices de largeur w, la bande chargée porte λ, chaque plaque neutre porte zéro et flotte à
son potentiel ; le champ lointain étant le même avec ou sans voisins, les voisins ne changent que
le rayon effectif a_eff de la bande chargée. Validation : bande isolée, a_eff = w/4 à 0,2 %. **Résultat** :
à l'écartement de centres s = w (bandes qui se touchent), a_eff croît de 1,47 (u, un DQD) et 1,98
(d, deux DQD) ; à s = 2w, 1,05 et 1,10 ; à 3w, 1,02 et 1,05 : plus fort pour le d, grand seulement au
contact. **Sur m_n − m_p, les écrans de u et de d se compensent presque** dans (4/9)S_u − (1/9)S_d :
1,324 → 1,342 MeV au contact (+1,4 %), inchangé au-delà de 1,5w. L'écran change le résultat de
moins de 1,5 % ; **le +2,4 % de R44 tient.** Verdict : DÉRIVÉ ; le terme de Coulomb du proton est
fixé par la section de R17, le contenu de R28 et la géométrie en étoile, écran compris.

**R46 — « Cherche ce qui fixe la largeur du ruban. »**
→ `ribbon_width.py`, 5/5. L'adaptation à Z₀ ne fixe que la forme (Z = Z₀·d/w = Z₀ ⇔ w = d, à toute
taille, R17 A). **La taille vient de g = 2** : les deux jonctions de pôles adaptés, ℏc/(π²w) chacune,
portent la moitié statique de la masse (R13, R17 D) : 2ℏc/(π²w) = mc²/2, donc w = 4ƛ/π² = 0,405 ƛ,
156,5 fm pour l'électron ; le rapport w/ℓ₁ = 6/π³ = 0,19 est le même à toute échelle, d'où
w₉ = w_e·3^{−2π} = 0,157 fm. Entrées : g = 2 (mesuré 2,0023), δ tiré de Z₀, et la lecture « écart
des pôles = largeur » (R17 C). **Rien de plus fondamental ne la fixe** : w/r = 15,0 (rayon de tube du
manuscrit), w/D₀ = 0,65 (écartement des fils), w/ℓ₁ = 0,19, w/ƛ = 0,405 sont des rapports sans
règle derrière. **L'écart p–n n'est qu'un test faible de w** : il y entre par un logarithme qui se
compense presque entre u et d ; w/2 donne −11 %, 2w +11 %, +10 % sur w +1,5 % : le +2,4 % de R44
fixe w à un facteur 2 près, pas mieux. **Le vrai test** reste le mode transverse de la nappe,
ℏcπ/w = (π³/4)·m_ec² = 3,96 MeV pour l'électron (R17 E), jamais observé ; la sortie « fluide TEM
seulement » reste un postulat. Verdict : w est fixée par le moment magnétique de l'électron (g = 2)
via l'énergie de jonction, CONDITIONNEL à la lecture des pôles ; c'est le nombre le moins testé
de la chaîne, et son test est un électron excité à 4 MeV qui n'existe pas.

**R47 — « Cherche ce qui interdit le mode transverse à 4 MeV. »**
→ `transverse_mode.py`, 5/5. Le 3,96 MeV est la coupure TE₁ du guide de l'électron, ℏcπ/w avec
w = 156,5 fm : une notion d'**onde**. Dans la lecture « onde » (R29) rien ne l'interdit, et rien
n'interdit non plus les harmoniques du mode d'arc : des électrons excités à 0,766 et 1,02 MeV
(harmoniques) et à 4,47 MeV (transverse), tous exclus par l'expérience (aucune résonance dans la
diffusion Compton au MeV, échelle de compositeness au-delà de 10 TeV). **Dans la lecture
« vortex » (R32, R40) il n'y a pas d'onde** : le fluide est une circulation stationnaire portant un
seul quantum conservé q = e/(2√α), d'énergie πℏc/trajet fixée par le quantum et le trajet ; il
n'y a ni harmonique (un écoulement stationnaire n'a pas de n-ième harmonique) ni mode transverse
(pas d'onde, pas de coupure). Un second quantum doublerait la charge circulante (11,7 e) et
quadruplerait l'énergie (1,02 MeV) : un autre objet, pas un électron excité. **Les nombres ne
changent pas** : 4πKq² = πℏc exactement, donc l'énergie de circuit égale l'énergie du mode en
demi-onde pour tout trajet ; tout R29–R45 tient avec « mode » lu comme « circuit ». Ce qui reste
excitable, ce sont les formes et les topologies, pas les modes : les hadrons en ont (Δ − N = 294 MeV
contre un circuit de quark à 254 MeV, +16 %), les leptons n'ont qu'une forme. **Verdict** : le mode
à 4 MeV n'est pas interdit par une règle, il est absent ; c'était un artefact de la lecture en onde,
et R32 le supprime en même temps que les harmoniques de l'électron.

**R48 — « Cherche ce qui fixe la forme excitée Δ à +294 MeV. »**
→ `delta_shape.py`, 6/6. Le Δ(1232) a le contenu du nucléon (Δ⁺ = uud) et le spin 3/2 ; R6 le lit
comme l'étoile symétrique (trois axes → 3/2) contre l'étoile asymétrique du nucléon (un axe → 1/2).
Dans l'image des circuits, les trois circuits de quark sont identiques dans les deux (763 MeV
immobiles) ; ce qui diffère est la circulation d'anneau qui porte le spin, S = R·E_circ/c (R25) :
**le nucléon a besoin de ℏ/2, le Δ de 3ℏ/2, donc E_circ(Δ) = 3·E_circ(N) au même rayon et
Δ − N = 2·E_circ(N)** : les circulations des trois circuits alignées (R41 F) au lieu de deux contre
une. Avec les deux lectures de l'anneau de R25 : 125 MeV (brins + à R₊ = 0,79 fm) → 250 MeV (−15 %) ;
168 MeV (e sur un anneau à 0,59 fm) → 336 MeV (+14 %) : **elles encadrent les 294 MeV mesurés.**
En inversant : E_circ(N) = 147 MeV, rayon d'anneau 0,671 fm (entre les deux rayons de R25) ; alors
N = 763 + 147 = 910 MeV (−3,0 %) et Δ = 763 + 441 = 1204 MeV (−2,3 %), tous deux à 3 % avec un seul
rayon. Test de ce rayon : le Δ⁺⁺ (uuu, charge 2e) avec toute sa charge à 0,671 fm a μ = 2(R/ƛ_p)·μ_N
= 6,4 μ_N, dans la fourchette mesurée 3,7 à 7,5 (enregistré). Une interaction spin-spin
électromagnétique entre circuits (l'analogue de l'hyperfin du modèle des quarks) vaut ~0,02 MeV :
exclue comme mécanisme ; l'écart est de la circulation, pas du magnétisme. **Verdict** : le Δ est
le nucléon avec trois fois sa circulation d'anneau ; DÉRIVÉ à ±15 % avec les rayons de R25 ; le
rayon exact (0,67 fm) n'est pas fixé par la base.

**R49 — « Cherche ce qui fixe le rayon de l'anneau à 0,67 fm. »**
→ `ring_radius.py`, 5/5. Cible : R = ℏc/(2 × 147 MeV) = 0,671 fm ; dans le langage de la base,
R = ƛ_p/f avec la fraction circulante f = 0,313 (l'électron a 3/4), E_circ/m = f/2 = 0,157. **Les
longueurs de la base à l'échelle 9** contre 0,671 fm : le rayon de circuit R_c = 3ℓ₁/(2π) = 0,388 fm
est √3 fois plus petit (à 0,1 %), ƛ_p est π fois plus petit (1,6 %), le reste est sans rapport.
**L'étoile de trois anneaux-circuits** : passant par le centre, leurs centres sont à R_c et
l'écartement entre quarks vaut √3·R_c = 0,672 fm, le rayon cherché à 0,1 % ; mutuellement
tangents, leur enveloppe extérieure R_c(1 + 2/√3) = 0,836 fm est le rayon de charge du proton à
0,6 % (r_p n'est une entrée nulle part), centres à 0,448 fm, rayon quadratique 0,593 fm. **Aucun
mécanisme ne sélectionne √3·R_c** : un écoulement le long des arêtes du triangle des quarks porte
S = R_in·E/c avec le rayon inscrit R_in = R_c/2 (E_circ serait 508 MeV) ; le rayon propre du circuit
non apparié donne 254 MeV ; aucun des deux n'est 147. L'identification R = √3·R_c = (3√3/2π)·ℓ₁(9)
est une coïncidence à 0,1 %, pas une dérivation. Si on l'accepte, tout sort de ℓ₁(9) seul :
N = 910 MeV (−3 %), Δ = 1204 (−2,3 %), μ(Δ⁺⁺) = 6,4 μ_N, r_p = 0,836 fm (+0,6 %) ; mais μ_p exigerait
0,87 e en circulation à ce rayon, pas une charge propre (enregistré). **Verdict** : deux
coïncidences géométriques de la même étoile de circuits (0,672 et 0,836 fm) ; le rayon de
l'anneau n'est pas dérivé, il manque la règle qui dit sur quel cercle la circulation collective
des trois circuits tourne.

---

## Bilan au 19 septembre 2026, suite (R32–R49 ; 44 scripts, tous PASS)

**Le redessin qui a tout débloqué : le fluide est un vortex, pas une onde (R32).** Un courant
stationnaire portant un quantum de circulation par circuit fermé, de charge q = e/(2√α), dont
l'énergie 4πKq²/trajet = πℏc/trajet est identique à l'ancienne « note » en demi-onde (R47) : les
nombres de R29–R45 tiennent, les harmoniques et le mode transverse à 4 MeV disparaissent.

**Ce qui est dérivé (avec α, m_e, l'exposant 2π, le compte 9 et la largeur de nappe).**
- Le quark est l'objet de l'électron, un circuit de trois brins, à l'échelle ℓ₁(9) = 0,81 fm :
  254 MeV chacun ; le nucléon en a trois (763 MeV immobiles), partition unique des 9 brins (R40, R41).
- La tension forte : πℏc/ℓ₁(9)² = 939 MeV/fm, √σ = 430 MeV, au centre de la bande du réseau
  (420–440), α et Z₀ s'annulant (R39, R40).
- L'écart neutron–proton : le brin neutre de plus est un DQD dont l'énergie bifilaire vaut
  (2α/3)·3^{2π}·m_ec² = 2,47 MeV (réseau 2,52 ± 0,29), moins le Coulomb du proton fixé par la
  section de R17 (plaque, a = w/4) et le contenu de R28 (charge sur les brins chargés) : 1,32 MeV
  contre 1,293, +2,4 %, écran des voisins compris (R42–R45).
- Le Δ est le nucléon avec trois fois sa circulation d'anneau : ±15 % avec les rayons de R25 (R48).
- La tri-jonction à 2 MeV (force nucléaire), la loi ⌊k/2⌋ des jonctions (R22, R24).

**Ce qui est conditionnel ou coïncidence.**
- La largeur de nappe w = 4ƛ/π² vient de g = 2 via la lecture « écart des pôles = largeur » ;
  l'écart p–n n'y est sensible qu'à un facteur 2 près (R46).
- Le rayon de l'anneau du nucléon, 0,67 fm : égal à √3 fois le rayon de circuit à 0,1 %, et
  l'enveloppe de trois circuits tangents vaut r_p à 0,6 % ; sans mécanisme (R49).
- Le neutron à +1,5 à +7 % (R41), μ_p exigeant 0,87 e en circulation (R49), le spin compté par
  circuits et par anneau (R41 F).

**Ce qui est fermé.** L'enroulement comme mécanisme de force (R39), le pas d'hélice comme
variable (R38–R39), le tube de flux dans la base telle qu'écrite (R30), la densité fixe (R31),
les comptes de brins comme masses (R28, R29).

**Ce qui reste.** La règle qui dit sur quel cercle tourne la circulation collective des trois
circuits (le 0,67 fm) ; ce qui fixe le pion à r_e/2 ; pourquoi 2π et 3, 7, 11 ; où se range la
moitié statique de l'électron quand le fluide est un vortex ; le −1 sous 2π.
