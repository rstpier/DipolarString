# Cahier des charges — reconstruire l'objet (21 novembre 2026)

Décision de l'auteur (20 novembre 2026) : l'élément de la CD, le filament dipolaire, est à reconstruire de zéro, en
partant de son **énergie**, pas de sa forme. La revue du 19 novembre (`../../REVUE_SCIENTIFIQUE_2026-11-19.md`, §2.2)
a établi que le filament de V2.10 ne peut former ni particule stable en trois dimensions ni charge par sa seule
topologie ; toute reconstruction change donc l'élément, c'est-à-dire pose une hypothèse nouvelle. Ce cahier des
charges dit ce que la nouvelle hypothèse doit contenir, ce qu'elle doit respecter, et ce qui a déjà été essayé.

## 1. Les trois exigences

1. **Une dynamique avant une géométrie.** Une énergie E(taille, …) avec au moins deux termes qui tirent en sens
   contraire (R113 : tout en 1/R ne fixe rien ; R120 : tension × longueur + quantum/R fixe un rayon). La géométrie
   (2cosh π, neuf brins) doit en sortir, pas y entrer.
2. **Ce que le vide adapté à Z₀ n'a pas, posé explicitement** (`../CONSTRAINTS.md` §1–§3) : une échelle (tension σ,
   ℓ₁ ou densité de saturation u_s, déclarée comme entrée) ; un terme qui résiste à l'effondrement en 3D (quartique
   en gradients, saturation, ou tension d'une boucle fermée) ; le lieu de la charge (l'équation d'état du fluide
   e/3, ou la charge comme entrée).
3. **Classique ou quantique, décidé avant le premier calcul** (revue §5). Sans amplitudes sur les configurations
   (axiome A3), tout lump est un objet classique : ni spin ½, ni statistique de Fermi, ni Bell. Avec A3, la CD
   devient une ontologie pour les champs quantiques et la question change (« qu'ajoute le milieu ? »).

## 2. Les bornes et acquis à respecter

| Contrainte | Valeur | Source |
|---|---|---|
| Électron ponctuel, G_E = 1 jusqu'à 10⁻³ fm ; μ_B et g = 2 par Dirac | aucun électron étendu | R100–R103 |
| Non-linéarité du vide à basse énergie | σ_γγ < 7·10⁴ × QED à 116 eV → u_s > 5,7·10²⁶ J/m³ = 1,2·10³ u₀ pour toute Z(u) | R119 (HIBEF, oct.–nov. 2026) |
| Le vide réel sature par création de paires | u_Schwinger = 7,7·10²⁴ J/m³, 73 fois sous la borne HIBEF | R122 |
| La tension fixe un rayon de boucle fermée | R* = √(κℏc/2πσ) ; 0,13 fm et 1,53 GeV pour σ_base, κ = ½ | R120 |
| Premier terme quantique d'une corde | Lüscher −πℏc/12r = −64 MeV à ℓ₁(9) | R120 |
| Pas de soliton 3D dans l'action de V2.10 ; réparation chirale exclue | — | `../CONSTRAINTS.md` §3 |
| Pas d'échelle de masse dans un vide EM adapté | les masses absolues sont des entrées | `../CONSTRAINTS.md` §1 |
| Aucun niveau excité de l'électron | e* > 100 GeV ; g − 2 à 10⁻¹² | R121 |

## 3. Les candidats et le test qui les tue

| Candidat | Ce qu'il donne | Ce qui le tue | État |
|---|---|---|---|
| Milieu de Kerr saturable, Z(u) (idée de base de l'auteur) | lumps stables au-dessus d'un seuil ; loi M = Ñ_min^{1/4} E_s N^{3/4} (le sac) | électron : E_s requis 0,14 MeV, 6·10⁴ sous la borne HIBEF en u_s ; toute non-linéarité qui piège agit au-delà de Schwinger | **étape 1 faite (R122) : exclu pour l'électron ; sac pour les hadrons avec E_s posé** |
| Boucle fermée de Nambu–Goto + quantum de circulation | rayon et masse fixés (0,13 fm, 1,5 GeV) ; tours de Regge | pas de spin ½ ; brisure de corde ; Lüscher dans le mauvais sens pour le nucléon | R120 : identité avec la règle posée, pas d'objet nouveau |
| Oscillateurs LC internes (« coût = énergie pour masse ») | masse = énergie interne (Zitterbewegung) | niveaux excités à 612, 713 keV… exclus | R121 : exclu |
| Soliton topologique (Hopfion, terme de Faddeev) | stabilité sans seuil | aucun porteur dans la trame de V2.10 | `../CONSTRAINTS.md` §3a : fermé |
| Ajouter A3 (amplitudes) | spin ½, statistique, Bell deviennent possibles | ce n'est plus la reconstruction d'un objet mais d'un cadre | ouvert : décision de l'auteur |

## 4. Résultat de l'étape 1 (R122, `step1_saturable_selftrapping.py`, 5/5)

L'idée « l'énergie monte Z là où elle est ; assez de quanta au même endroit se piègent » écrite comme milieu de
Kerr saturable donne, pour l'enveloppe d'une porteuse, l'équation de Schrödinger non linéaire saturable. En 1D tout
paquet se piège ; en 2D le seuil est celui de Townes (Ñ_c = 11,70, contrôle de la méthode) ; en 3D cubique rien
n'est stable (effondrement) ; en 3D saturable il existe un seuil universel **Ñ_min = 158/√x_s** (x_s = u_sat/u_s),
atteint à μ/x_s = 0,049, avec une branche stable au-dessus. Retour aux unités : N_q = Ñ (E_s/E_c)⁴ quanta,
M c² = Ñ_min^{1/4} E_s N_q^{3/4}, avec E_s = (u_s (ℏc)³)^{1/4} l'échelle d'énergie de la saturation : **la loi de masse
du sac** (MIT : M ∝ B^{1/4} N^{3/4}).

- Échelles : E_s(u₀) = 0,39 MeV ; E_s(borne HIBEF) = 2,3 MeV ; E_s(Schwinger) = 0,78 MeV ; E_s(Euler–Heisenberg) =
  9,2 MeV ; E_s(sac, B^{1/4} = 145 MeV) = 145 MeV.
- Le lump classique le plus léger (10 quanta) à la borne HIBEF pèse 46 MeV, 90 fois l'électron. L'électron comme lump
  d'un quantum exigerait E_s = 0,14 MeV, u_s = 9·10²¹ J/m³, 6·10⁴ fois sous la borne : **exclu**. Un nucléon de trois
  quanta exige E_s = 116 MeV, u_s = 4·10³³ J/m³ = 8·10⁹ u₀ : l'échelle du substrat n'entre pas ; c'est une seconde
  échelle, posée, celle du sac.
- La borne HIBEF met u_s à 73 fois la densité de Schwinger : toute non-linéarité assez forte pour piéger agit là où
  le vide réel crée déjà des paires ; la « saturation » du vide QED est la création de paires, c'est-à-dire
  l'électron lui-même. Et un lump de quanta sans masse libres de sortir est métastable (oscillon) : la stabilité
  exige une charge conservée (boucle fermée) ou un extérieur interdit (confinement).

**Verdict.** Le mécanisme est juste et c'est un acquis : il fournit la dynamique demandée (deux termes opposés, un
seuil, une loi de masse) et il explique pourquoi « deux quanta ne suffisent pas ». Il est **exclu comme origine de
l'électron** et il est, pour les hadrons, **le modèle du sac** avec une échelle E_s ≈ 120–200 MeV posée et sans lien
avec ℓ₁. L'étape 1 ne fournit pas d'objet nouveau ; elle fixe ce que l'objet doit avoir : une charge conservée ou un
extérieur interdit, et une échelle qui n'est pas celle du substrat.

## 5. Décision proposée (le physicien responsable)

La route « énergie → Z » est fermée pour l'électron et redondante pour les hadrons (le sac et la QCD existent). Il
reste deux voies cohérentes avec tout le dossier : **(a)** restreindre la CD au milieu (états limites gravitationnels),
comme la déclaration de portée du manuscrit le fait déjà ; **(b)** adopter A3 et redéfinir le programme comme
ontologie d'un milieu pour les champs quantiques, en demandant ce que le milieu ajoute (ℓ₁, une coupure, un repère).
Aucune étape 2 n'est lancée sans ce choix, qui appartient à l'auteur.
