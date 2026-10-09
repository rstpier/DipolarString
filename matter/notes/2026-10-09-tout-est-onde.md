# Tout est énergie ondulatoire — une théorie minimale des champs et de la matière

*Note de méthode, 9 octobre 2026. Demande de l'auteur : « bâtir une théorie des champs et de la
matière en présumant que tout est énergie ondulatoire ». Construite comme le redessin : des
axiomes, ce qui en découle, ce qui reste posé, ce qui tuerait l'ensemble. Pas de calcul nouveau ;
chaque étape « dérivée » est un théorème connu, cité. Le résultat est annoncé d'avance : cette
théorie existe, elle s'appelle la théorie quantique des champs, et l'exercice sert à voir
exactement où la prémisse « tout est onde » porte, et où elle s'arrête.*

## 0. La règle du jeu

Il n'y a qu'une sorte de chose : un milieu qui remplit l'espace, et de l'énergie qui y ondule.
« Particule » n'est pas un objet premier ; c'est le nom d'un quantum d'une onde du milieu. On
s'interdit tout point matériel, toute bille, tout nœud rigide. On note **D** ce qui est dérivé des
axiomes, **P** ce qui est posé.

## 1. Six axiomes

**A1 — Le milieu.** Un milieu occupe l'espace ; son état classique est une configuration de champs
φ(x) (une ou plusieurs composantes par point) ; ses ondes ont une vitesse finie c. *(La base DS en
donne une réalisation : des lignes à Z₀, c₀ = 1/√(LC).)*

**A2 — Ondes linéaires et locales.** Une petite perturbation obéit à une équation d'onde linéaire
où chaque point ne dépend que de ses voisins. Pour un mode, la plus générale est
□φ + (ω₀/c)² φ = 0 : une vitesse c commune et, par mode, une fréquence de coupure ω₀ qui peut
être nulle.

**A3 — Amplitudes sur les configurations.** L'état du milieu n'est pas une configuration φ(x),
c'est une amplitude Ψ[φ] sur toutes les configurations. Les amplitudes s'additionnent, les
probabilités sont leurs carrés, l'évolution est linéaire et conserve la norme. *(C'est le seul
axiome qui n'est pas classique ; c'est celui que la base DS n'a pas encore.)*

**A4 — Symétrie.** Le milieu paraît le même sous translation, rotation, et changement de vitesse
uniforme : une seule vitesse c pour tous ses modes, et les équations gardent leur forme (Lorentz).
*(Pour un milieu réel, éther ou réseau, c'est une symétrie émergente ; le repère privilégié doit
rester caché sous 10⁻²⁰, voir §5.)*

**A5 — Couplages locaux.** Les modes se couplent entre eux, localement, avec une force (une
charge) par couple de modes ; le couplage respecte A4. Pour un mode vectoriel, respecter A4
impose l'invariance de jauge, donc la conservation de la charge.

**A6 — Modes à demi-tour.** Certains modes changent de signe sous une rotation de 2π (spineurs).
Les autres n'en changent pas.

## 2. Ce qui en découle

Chaque ligne est un théorème ; la référence dit où il est démontré.

| conséquence | de | statut |
|---|---|---|
| Une onde de fréquence ω porte son énergie par paquets ℏω ; une « particule » est un quantum d'un mode ; E = ℏω, p = ℏk (de Broglie) | A3 appliqué à un oscillateur | D |
| **Masse = fréquence de coupure** : mc² = ℏω₀ ; dispersion E² = (pc)² + (mc²)² ; les modes à ω₀ = 0 vont à c (photon), les autres moins vite ; portée statique ℏ/mc (Yukawa 1935) | A2 + A4 | D |
| Le quantum au repos a une horloge interne mc²/ℏ ; pour un spineur, un va-et-vient entre ses deux moitiés à 2mc²/ℏ (Zitterbewegung) | A2 + A6 | D |
| Spin ½ : le mode à demi-tour a deux composantes, le moment cinétique ℏ/2, et g = 2 au premier ordre du couplage (Lévy-Leblond 1967) | A5 + A6 | D |
| **Antiparticules** : les solutions de fréquence négative de A2 doivent être relues (A3) comme des quanta de charge opposée ; création de paires dès que E ≥ 2mc² (Dirac 1930, Feynman–Stückelberg). C'est ici que « la matière est de l'énergie ondulatoire » devient littéral : un photon de 1,022 MeV devient un électron et un positron | A2 + A3 + A4 | D |
| **Forces** : échanger un quantum d'un mode entre deux autres donne une force ; 1/r pour un mode vectoriel sans masse (Coulomb), e^{−r/λ}/r pour un mode massif ; le signe dépend du spin du mode échangé : spin 1 fait se repousser les charges de même signe, spin 0 et 2 les attirent (théorème) | A5 | D |
| **Statistique** : les quanta des modes à demi-tour ne peuvent pas occuper deux fois le même état (Pauli) ; ceux des modes entiers s'entassent (Pauli 1940, spin-statistique). D'où : la matière est faite des premiers (elle a des couches, une chimie, une rigidité), les forces des seconds (elles se superposent en champs classiques) | A3 + A4 + A6 | D |
| **Intrication** : après un couplage, l'amplitude conjointe de deux modes ne se factorise plus ; violation de Bell ; aucun transport, puisque A3 vit sur les configurations et non dans l'espace | A3 + A5 | D |
| Lois de conservation : énergie, impulsion, moment cinétique, charge, chacune d'une symétrie (Noether 1918) | A4 + A5 | D |
| Le vide a un demi-quantum par mode : force de Casimir (mesurée au pourcent), déplacement de Lamb, et une intrication entre régions voisines proportionnelle à leur surface | A3 | D |
| Rayonnement thermique : loi de Planck pour le mode sans masse | A3 | D |
| Corrections au deuxième ordre du couplage : moment anormal α/2π, vérifié à 10⁻¹² | A5 en perturbation | D |
| Rigidité de la matière : aucun son plus rapide que c ; un solide est un état intriqué stable, pas instantané | A2 + Pauli | D |

Tout cela sort de six axiomes et d'aucun nombre. C'est le contenu de la théorie quantique des
champs, et il est vérifié partout où on l'a regardé.

## 3. Ce qui reste posé

| entrée | nombre | statut |
|---|---|---|
| la liste des modes : quels champs existent, avec quelles symétries internes (3 copies de 4 modes de matière, 4 modes vectoriels, 1 scalaire) | structure | P |
| les forces de couplage : α, α_s, le couplage faible | 3 | P |
| les fréquences de coupure (masses) : 12 fermions, W, Z, Higgs ; ou, de façon équivalente, le fond uniforme v qui fixe le taux de conversion gauche↔droite et 12 couplages à ce fond | 15 | P |
| les angles de mélange entre copies (CKM, PMNS) | 8 | P |
| pourquoi trois copies | structure | P |

La prémisse « tout est onde » explique **le mécanisme de tout** et **la valeur de rien**. Elle
dit ce qu'est une masse (une coupure), pas pourquoi 0,511 MeV. Elle est exactement aussi muette
que le Modèle standard sur ces nombres, parce qu'elle **est** le Modèle standard, écrit depuis
l'autre bout.

## 4. Le noyau mathématique (une page)

Pour fixer les mots, voici les équations qui portent chaque axiome.

```
A1–A2  mode scalaire :    (1/c²) ∂²φ/∂t² − ∇²φ + (ω₀/c)² φ = 0
       mode à demi-tour : i ℏ ∂ψ/∂t = −i ℏ c σ·∇ ψ          (paire de Weyl, sans masse)
       masse = conversion : + m c² (ψ_L ↔ ψ_R)               (taux m c²/ℏ)
A3     [φ(x), π(y)] = i ℏ δ(x − y)   ;   Ψ[φ] évolue par i ℏ ∂Ψ/∂t = H Ψ
       énergie d'un mode : (n + ½) ℏ ω   ;   n = nombre de quanta
A4     une seule vitesse c dans toutes les équations ; forme invariante sous Lorentz
A5     ∂ → ∂ − i q A   ;   charge q par mode ; invariance de jauge ⇒ ∂·j = 0
A6     rotation de 2π : ψ → −ψ   ;   spin ½   ;   antisymétrie à l'échange (Pauli)
```

Rien ici n'est nouveau. Ce qui est utile, c'est de voir que la ligne « masse = conversion » est la
seule où entre un nombre par espèce : c'est là, et seulement là, que la théorie cesse d'être
dérivée.

## 5. Ce qui tuerait cette théorie, et ce qui a été cherché

| prédiction de « tout est onde » | test | état |
|---|---|---|
| une seule vitesse limite pour tous les modes | électron : |v_max − c|/c < 10⁻¹¹ (Hohensee 2009) ; neutrinos : < 2·10⁻⁹ (SN 1987A) ; photons : pas de dispersion jusqu'à l'énergie de Planck (Fermi, sursauts gamma) | tient |
| aucun quantum n'a de taille propre | électron < 10⁻¹⁸ m (LEP) | tient |
| superposition exactement linéaire (A3) | non-linéarités bornées sous 10⁻²⁰ (spectroscopie, Weinberg 1989) ; superposition d'une molécule de 25 000 u (2019), d'un cristal de 16 μg (2023) | tient |
| E = ℏω pour tout mode | photons, électrons, neutrons, atomes, molécules | tient |
| seuil de création de paires exactement 2mc² | 1,022 MeV | tient |
| Pauli exact | violation bornée sous 10⁻²⁸ (VIP) | tient |
| énergie du vide réelle | Casimir au pourcent, Lamb, g − 2 à 10⁻¹² | tient |
| aucun transport plus vite que c, y compris par l'intrication | Bell sans faille (2015), non-signalement à chaque test | tient |
| pour un milieu réel : un repère privilégié | aucune anisotropie de c sous 10⁻¹⁷ (résonateurs), aucune dépendance de la vitesse de la Terre | caché sous 10⁻¹⁷ à 10⁻²⁰ |

La dernière ligne est le prix spécifique d'un **milieu** : si A1 est pris au sérieux (un éther, un
réseau), A4 n'est qu'approché, et la différence doit se cacher sous les bornes existantes. Un
réseau avec un tic doit donc avoir son pas très au-dessous de tout ce qu'on sonde, ou être
invariant au premier ordre comme l'automate de Weyl (R111).

## 6. Ce qui manque encore à tout le monde

- **La gravitation.** Dans ce langage, c'est le milieu lui-même qui se déforme : sa vitesse de
  propagation dépend de la densité d'énergie, et ses ondes voient une géométrie courbe (c'est
  exactement ce qu'un fluide en mouvement fait à ses ondes sonores, Unruh 1981). Faire obéir cette
  déformation à A3, c'est la gravité quantique ; personne ne l'a.
- **Les nombres du §3.**
- **La mesure.** A3 appliqué à l'appareil donne une superposition d'aiguilles ; on en voit une.
  Everett garde A3 tel quel, les modèles d'effondrement le cassent au-delà d'une masse (testable,
  en cours), Bohm le complète. Indécidé.

## 7. Ce que ça change pour la base DS

La base est une réalisation de A1 (lignes, Z₀, c₀) qui essaie de **dériver** A4 (l'isotropie par le
nœud de Johns) et les nombres du §3 (comptes de brins, échelle 2π). Le redessin a établi :

- A1 concret marche pour le photon (R106 : le SCN donne Maxwell) ;
- il ne porte pas le spineur à c (R107–R110), sauf à changer le nœud (R111–R112) ;
- les nombres restent posés (R116) ;
- **A3 manque.** Les lignes de la base portent des ondes, V(x, t), pas des amplitudes Ψ[V]. C'est
  la superposition des ondes (linéarité classique, Bell ≤ 2), pas la superposition des
  configurations. Sans A3 : pas d'intrication, pas de Pauli, pas de création de paires, pas de
  Casimir, pas de g − 2 ; avec A3, le SCN devient un automate quantique et la base devient une
  théorie quantique des champs sur réseau avec repère privilégié, c'est-à-dire exactement ce que
  D'Ariano construit, avec les mêmes 15 nombres à poser.

La prescription tient en une ligne : **garder les lignes, remplacer l'onde par l'amplitude sur les
ondes.** C'est un axiome de plus, et c'est le seul qui transforme un modèle de milieu en théorie de
la matière. Il ne rapproche d'aucun des nombres ; il rend tout le reste possible.

## 8. Bilan

- Dérivé (six axiomes, zéro nombre) : la cinématique, la masse comme coupure, le spin, les
  antiparticules, les forces et leur signe, la statistique, l'intrication, les conservations, le
  vide et ses effets, la rigidité de la matière.
- Posé : la liste des modes, 3 couplages, 15 coupures, 8 angles, le nombre 3.
- Ouvert : la gravité quantique, la mesure, les nombres.
- Pour la base : un axiome à ajouter (A3), un nœud à changer (R112), des nombres à trouver
  ailleurs que dans la géométrie des lignes (R116).
