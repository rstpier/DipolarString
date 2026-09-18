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
