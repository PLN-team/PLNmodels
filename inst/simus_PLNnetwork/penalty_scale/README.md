# Quelle échelle de pénalité par défaut pour les réseaux ?

2 octobre 2026 — version 1.3.2.9350, plancher ciblé actif (`latent_floor = 1e-3`).

Comparaison de `penalty_scale = "covariance"` (défaut actuel) et `"correlation"` sur des
données simulées à réseau connu. 1 960 jeux de données, chacun ajusté sur les deux échelles.

## Plan

`sim_penalty_scale.R`, p = 40 sauf mention contraire, chemins de 30 pénalités
(`min_ratio = 0.05`), 20 répétitions par configuration (10 pour StARS).

| bloc | configurations | ce qui varie |
|---|---|---|
| `clean` | 64 | graphe (aléatoire, à hubs, en communautés, en chaîne) × (n, p) ∈ {(50, 40), (100, 40), (200, 40), (100, 80)} × variances latentes (unité, ou écarts-types de 0.5 à 2) × abondances (fortes, ou mêlées) |
| `contaminated` | 12 | graphe (aléatoire, à hubs) × n ∈ {50, 100, 200} × contamination (3 espèces absentes au hasard, ou 8 espèces ne vivant que dans un groupe d'échantillons sur trois) |
| `covariate` | 12 | graphe × n × variances, avec une covariable dans les données et dans le modèle |
| `zi` | 8 | données avec 20 % de zéros ajoutés ; `ZIPLNnetwork()` et `PLNnetwork()` |
| `stars` | 4 | graphe × variances, n = 100, sélection par StARS |

Mesures, sur les arêtes entre espèces non contaminées : F1 au modèle de la taille du vrai
réseau, meilleur F1 du chemin, aire sous la courbe précision-rappel, F1 aux modèles retenus
par BIC, EBIC et StARS. Comparaisons appariées (même jeu de données, deux échelles).

## Résultats

**L'échelle des corrélations fait mieux dans les 100 configurations, sans exception**, sur le
F1 à taille vraie comme sur le meilleur F1 du chemin.

| bloc | F1 à taille vraie, covariance | corrélation | répétitions gagnées par la corrélation | gain moyen par configuration |
|---|---|---|---|---|
| `clean` | 0.576 | 0.740 | 1 209 sur 1 280 | de +0.008 à +0.389 |
| `contaminated` | 0.265 | 0.784 | 230 sur 240 | de +0.031 à +0.839 |
| `covariate` | 0.618 | 0.802 | 226 sur 240 | de +0.013 à +0.331 |
| `zi` | 0.230 | 0.281 | 136 sur 160 | de +0.036 à +0.066 |
| `stars` | 0.632 | 0.818 | 37 sur 40 | de +0.037 à +0.331 |

### Données saines (`clean`)

Le facteur décisif est l'hétérogénéité des variances latentes :

| variances | abondances | covariance | corrélation | répétitions gagnées |
|---|---|---|---|---|
| hétérogènes | fortes | 0.481 | 0.795 | 320 sur 320 |
| hétérogènes | mêlées | 0.420 | 0.664 | 320 sur 320 |
| unité | fortes | 0.759 | 0.811 | 284 sur 320 |
| unité | mêlées | 0.644 | 0.691 | 285 sur 320 |

Avec des variances hétérogènes, la pénalité sur l'échelle des covariances favorise les
arêtes des espèces à forte variance : le F1 perd 0.22 à 0.28 par rapport aux variances
unité, alors qu'il ne bouge presque pas sur l'échelle des corrélations (0.730 contre 0.751).
Même à variances égales, la corrélation gagne 0.05, parce que les variances estimées
diffèrent d'une espèce à l'autre.

Le gain tient pour les quatre graphes (de +0.128 pour les communautés à +0.191 pour la
chaîne) et pour toutes les tailles (de +0.148 à n = 200 à +0.186 pour p = 80).

### Données contaminées

| contamination | n | covariance | corrélation | part des arêtes sur les espèces contaminées (covariance → corrélation) |
|---|---|---|---|---|
| au hasard | 50 | 0.235 | 0.700 | 86 % → 1 % |
| au hasard | 100 | 0.535 | 0.813 | 64 % → 0 % |
| au hasard | 200 | 0.819 | 0.866 | 33 % → 0 % |
| par groupes | 50 | 0.002 | 0.709 | 100 % → 41 % |
| par groupes | 100 | 0.000 | 0.786 | 100 % → 48 % |
| par groupes | 200 | 0.000 | 0.832 | 100 % → 51 % |

Avec des absences par groupes, l'échelle des covariances ne retrouve rien, même à n = 200,
et son chemin n'atteint la taille du vrai réseau que dans 65 % des ajustements du bloc : la
grille est calée sur les covariances des espèces contaminées. (Pour les groupes, 40 % des
arêtes touchent une espèce contaminée dans le vrai réseau.)

### Sélection de modèle

| critère | bloc | covariance | corrélation | répétitions gagnées / perdues |
|---|---|---|---|---|
| BIC | `clean` | 0.466 | 0.559 | 996 / 278 |
| EBIC | `clean` | 0.304 | 0.340 | 588 / 452 |
| BIC | `contaminated` | 0.235 | 0.532 | 173 / 65 |
| StARS | `stars` | 0.594 | 0.713 | 36 / 4 |
| BIC | `stars` | 0.506 | 0.588 | 30 / 10 |

Le gain passe à la sélection, plus faiblement : le BIC retient environ deux fois trop
d'arêtes sur les deux échelles (une centaine pour 55 vraies), l'EBIC trop peu (une
trentaine). Sur le BIC, la corrélation fait mieux dans 55 configurations saines sur 64 ; à
variances unité les deux échelles sont à égalité (0.560 contre 0.547).

### Zéro-inflation

Avec 20 % de zéros ajoutés au hasard, la reconstruction est mauvaise sur les deux échelles
(F1 à taille vraie de 0.17 à 0.34), pour `ZIPLNnetwork()` comme pour `PLNnetwork()`, et BIC
comme EBIC retiennent le réseau vide. La corrélation garde un avantage de 0.05, mais ce bloc
ne dit rien de la sélection. Un taux de zéros plus faible serait plus instructif.

### Temps de calcul

Équivalent : 0.77 s par chemin contre 0.79 s sur données saines.

## Lecture

- Rien, dans ces simulations, ne plaide pour garder l'échelle des covariances.
- Le gain est le plus fort là où les données réelles se situent : variances latentes
  hétérogènes, espèces souvent absentes.
- Limites : une seule force d'arêtes, des graphes de degré moyen 2 à 3, p ≤ 80, pas de
  poids de pénalité, un seul taux de zéro-inflation, peu instructif.

## Décision

`penalty_scale = "correlation"` est le défaut depuis la version 1.3.2.9500, pour
`PLNnetwork()`, `ZIPLNnetwork()` et `ZIPLN()` en covariance creuse. Conséquences :

- **Les pénalités changent de sens.** Elles sont sans dimension, entre 0 et 1.
- **Retour au comportement précédent** : `PLNnetwork_param(penalty_scale = "covariance")`.
- **Pénalités fournies sans choisir l'échelle** : jusqu'à 1, elles sont prises telles
  quelles, avec un message par session ; au-dessus de 1, elles sont prises pour des
  pénalités sur l'échelle des covariances et divisées par le rapport entre la plus grande
  covariance résiduelle de l'inception et sa plus grande corrélation, avec un avertissement.
  Cette conversion fait correspondre les pénalités donnant le réseau vide sur les deux
  échelles ; elle est exacte à variances égales, et seulement indicative sinon (sur
  trichoptera, les pénalités 3, 1 et 0.3 donnaient 1, 19 et 35 arêtes ; converties, 0, 31
  et 56).
- Le critère n'est plus fixe : l'algorithme alterné cherche un point fixe, les poids de
  pénalité étant recalculés à chaque M-step.
- Avec `penalize_diagonal = TRUE`, l'échelle des corrélations gonfle les variances de façon
  multiplicative (`Σ_ii = S_ii (1 + λ)`) : combinaison à éviter.

## Fichiers

- `sim_penalty_scale.R` : `Rscript sim_penalty_scale.R <bloc> [répétitions] [cœurs]` ;
  résultats dans `sim_<bloc>.rds`. L'ensemble tourne en 6 minutes sur 16 cœurs.
- `summarize.R` : comparaisons appariées, `Rscript summarize.R [bloc ...]`.
