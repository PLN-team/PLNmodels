# Espèces dégénérées dans PLNnetwork : exploration de trois remèdes

1er octobre 2026 — branche `explore-degenerate`. Fait suite aux issues #180 et #184 et aux
PR #191 à #193.

## Le problème

Dans `PLNnetwork()`, une espèce souvent absente mais abondante quand elle est présente
devient un hub artificiel : elle se retrouve reliée à presque toutes les autres. Sur
`oaks ~ 1`, `f_OTU_1011` (présente sur les 38 arbres « intermediate », absente ailleurs)
porte 103 des 103 premières arêtes du chemin. Le phénomène existe aussi sur `barents` et
`mollusk`, avec ou sans pénalité sur la diagonale.

Deux mécanismes se cumulent.

1. **Effet d'échelle.** La pénalité ℓ1 du graphical Lasso porte sur les entrées de `Ω`,
   qui ne sont pas invariantes d'échelle : si le latent de l'espèce `j` est multiplié par
   `c`, ses `Ω_jk` sont divisés par `c`. Une espèce à forte variance latente a donc des
   arêtes presque gratuites. Or les zéros d'une telle espèce sont ajustés par des moyennes
   latentes très négatives, ce qui gonfle sa variance. Cet effet agit même sans explosion :
   dans la simulation à n = 200, la variance reste autour de 20 (aucune espèce signalée
   par `$degenerate_species`) et les espèces contaminées portent pourtant 34 % des arêtes
   au lieu de 15 %.
2. **Direction plate à l'infini.** Quand le latent de l'espèce est bien prédit par les
   autres (p proche de n), `M → −∞` sur les zéros et `S² → ∞` laissent l'ELBO presque
   inchangée : le maximum est atteint à l'infini, comme dans une séparation en régression
   logistique. La variance latente part alors en milliers (jusqu'à 50 000 sur `oaks`).

## Les trois pistes

Prototypes non invasifs (`proto.R` remplace `PLNnetworkfit$optimize` dans la session).

- **A. Plancher sur les moyennes variationnelles** : `M_ij ≥ log ε − O_ij`, avec `S²`
  ré-optimisé sur les cellules bornées. C'est une restriction de la famille
  variationnelle : le modèle et le critère (l'ELBO) sont inchangés, et la borne reste une
  borne inférieure de la même vraisemblance.
- **B. Exclusion du réseau** : les arêtes des espèces signalées reçoivent un poids de
  pénalité prohibitif (l'espèce reste dans la vraisemblance, comme nœud isolé), et le
  chemin est réajusté jusqu'à ce que plus rien ne soit signalé.
- **C. Pénalité à l'échelle des corrélations** : `ρ_ij = λ √(S_ii S_jj)`, recalculée à
  chaque M-step. Cela revient à appliquer le graphical Lasso à la matrice de corrélation
  résiduelle (`Ω = D^-1/2 Θ D^-1/2`), comme il est d'usage pour un modèle graphique
  gaussien. `λ` devient sans dimension, entre 0 et 1.

Une variante de A bornant le comptage attendu (`O + M + S²/2 ≥ log ε`) a été écartée : elle
ne stoppe pas l'explosion, puisque `S²` part à l'infini avec `M`.

## Simulation avec vérité terrain (`sim_known_network.R`)

p = 40 espèces, réseau creux connu, 20 répétitions. Trois espèces sont rendues absentes de
65 % des échantillons (un groupe différent par espèce), et le modèle est ajusté sans la
covariable correspondante. Les arêtes sont évaluées entre espèces saines ; la part des
arêtes touchant une espèce contaminée devrait être de 15 %.

Avec 3 espèces contaminées :

| méthode | n = 200 : part contaminée | n = 200 : F1 à taille vraie | n = 50 : part contaminée | n = 50 : F1 à taille vraie | n = 50 : F1 au BIC | n = 50 : variance max |
|---|---|---|---|---|---|---|
| actuel | 34 % | 0.76 | 93 % | 0.13 | 0.06 | 190 |
| A, ε = 1e-2 | 32 % | 0.79 | 80 % | 0.31 | 0.28 | 17 |
| A, ε = 1e-3 | 33 % | 0.78 | 85 % | 0.24 | 0.20 | 23 |
| B | 34 % | 0.76 | 53 % | 0.46 | 0.41 | 47 |
| **C** | **0 %** | **0.85** | **2 %** | **0.70** | **0.62** | 26 |
| C + A, ε = 1e-3 | 0 % | 0.85 | 2 % | 0.70 | 0.62 | 22 |

Sans contamination, le réglage actuel donne un F1 à taille vraie de 0.85 (n = 200) et 0.63
(n = 50) ; C donne 0.87 et 0.70. C ne dégrade donc rien sur des données saines, et fait
mieux quand n est petit.

Lecture :
- **C supprime l'artefact dans les deux régimes** et ramène la qualité au niveau des
  données saines. En contrepartie, les espèces contaminées perdent leurs vraies arêtes
  (part de 0 à 2 % au lieu de 15 %).
- **A contient la variance mais pas les hubs** : les espèces contaminées portent encore
  80 à 85 % des arêtes à n = 50.
- **B n'agit que si une espèce est signalée** : il est sans effet à n = 200, et partiel à
  n = 50.

## Données réelles

Au modèle retenu par le BIC (`real_bic.R`) :

| | actuel | A, ε = 1e-3 | B | C | C + A |
|---|---|---|---|---|---|
| espèces dégénérées sur le chemin, `oaks ~1` | 3 | 0 | 1 | 8 (variance 4e16) | 0 |
| idem, `oaks ~tree` | 7 | 0 | 0 | 0 | 0 |
| idem, `barents` | 1 | 0 | 0 | 4 (variance 5e10) | 0 |
| idem, `mollusk` | 1 | 0 | 0 | 0 | 0 |
| espèces exclues par B | – | – | 30 sur 114 (`oaks ~1`, non convergé en 10 tours), 19 (`oaks ~tree`), 7 sur 30 (`barents`), 4 sur 32 (`mollusk`) | – | – |

À taille de réseau fixée, environ p arêtes (`real_fixed_size.R`), degré maximal :

| | actuel | A, ε = 1e-3 | C | C + A |
|---|---|---|---|---|
| `oaks ~1` (p = 114) | 107 | 22 | 12 | 10 |
| `oaks ~tree` | 8 (106 à 2p arêtes) | 18 | 12 | 12 |
| `mollusk` (p = 32) | 26 | 7 | 11 | 11 |
| `barents` (p = 30) | 28 | 9 | 7 | 7 |

Lecture :
- **B produit bien la cascade redoutée** : exclure une espèce en fait basculer d'autres.
- **C seul peut encore diverger** (`oaks ~1`, `barents`) : sans ancrage d'échelle, la
  direction plate reste ouverte. Combiné au plancher, plus aucune espèce n'est dégénérée
  sur les quatre jeux, et les hubs passent d'un degré proche de p à un degré de 7 à 20.
- **A est une vraie régularisation**, sensible à ε : sur `mollusk`, 12 à 59 % des cellules
  sont à la borne et la loglik perd 50 à 650 points ; le réseau retenu par le BIC change
  de façon non monotone avec ε (`oaks ~tree` : 170, 222 puis 48 arêtes pour ε = 1e-2,
  1e-3, 1e-4).
- Sur `mollusk`, C laisse les trois espèces à zéros structurels en tête des degrés (37 %
  des extrémités d'arêtes pour 12 % attendus), avec un degré modéré de 11.
- Le BIC n'est pas un juge fiable sur ces données : il récompense la dégénérescence.

## Conclusion provisoire

La pénalité à l'échelle des corrélations (C) traite la cause, là où A et B traitent des
symptômes. Elle a besoin d'un garde-fou contre la divergence de la variance, rôle que
remplit le plancher (A) avec un ε petit. La combinaison C + A, ε = 1e-3 est la seule à
n'avoir ni espèce dégénérée ni hub artificiel sur l'ensemble des cas testés.

## Implémentation

Les pistes C et A sont maintenant des options du paquet, désactivées par défaut :
`PLNnetwork_param(penalty_scale = "correlation", latent_floor = 1e-3)`. Les scripts
`sim_known_network_package.R` et `real_package.R` refont la simulation et la comparaison sur
données réelles avec ces options, et redonnent les résultats des prototypes. Le récit complet
est dans `inst/devlog/DEVLOG_2026-09-30_10-01.md`.

## Limites

- Le plancher est prototypé par projection après le pas de Newton, non par un pas de Newton
  sous contrainte : les temps de calcul et les ELBO rapportés sont indicatifs.
- La simulation ne couvre qu'un type de contamination (absence dans un groupe aléatoire),
  une taille de réseau et une proportion de zéros.
- `ZIPLNnetwork()` n'a pas été exploré.
- Avec C, le critère pénalisé dépend de `S` à travers les poids : ce n'est plus la
  maximisation d'un critère fixe, mais un point fixe de l'algorithme alterné. Les
  propriétés de la sélection par BIC, EBIC ou StARS sur ce chemin restent à vérifier.

## Fichiers

- `proto.R` : prototypes (options `proto.floor`, `proto.corr`), `fit_path()`,
  `fit_exclusion()`.
- `sim_known_network.R` : simulation, `Rscript sim_known_network.R 20 50` pour 20
  répétitions à n = 50 ; résultats dans `sim_known_network_n50.rds` et `_n200.rds`.
- `real_bic.R`, `real_fixed_size.R` : données réelles ; résultats dans les `.rds` de même
  nom.
- `sim_known_network_package.R`, `real_package.R` : les mêmes comparaisons avec les options
  du paquet (`penalty_scale`, `latent_floor`), sans les prototypes.
