# Réglages d'optimisation de PLNnetwork : quel effet sur le réseau retrouvé ?

14 novembre 2026, PR #202 (version 1.3.2.9550 + défauts builtin de PLN).

Le script `sim_optim_settings.R` reprend les simulations et les mesures de `../penalty_scale/`
(échelle des corrélations, plancher actif). Il compare des réglages de la boucle alternée
(graphical Lasso, puis étape VE à Ω fixée) du backend builtin. On simule 10 répétitions par
configuration, avec les mêmes graines pour tous les réglages, et les comparaisons sont appariées.

| réglage | étape VE à Ω fixée | boucle alternée |
|---|---|---|
| `legacy` | défaut de la version installée (HEAD de #202 avant ces changements) : VE complète, puisque `maxit_ve = 1` était sans effet ; inception tronquée à 20 itérations EM | `ftol_em = 1e-5`, `maxit_em = 20` |
| `ve1_out5` (nouveau défaut) | 1 pas de Newton (`maxit_ve = 1`, désormais effectif) | `ftol_out = 1e-5`, `maxit_out = 20` |
| `ve1_out7` | 1 pas de Newton | `ftol_out = 1e-7`, `maxit_out = 1000` |
| `full_out5`, `full_out7` | VE complète, pas cellule par cellule | 1e-5 / 20, 1e-7 / 1000 |
| `cg5_out5`, `cg5_out7` | VE complète, pas de Newton-CG (`cg_max = 5`) | 1e-5 / 20, 1e-7 / 1000 |

Dans tous les nouveaux réglages, l'inception (PLN diagonal) a les défauts de `PLN()` : SQUAREM et
tests sur la borne avec ses constantes.

## Résultats

| bloc | réglage | F1 à taille vraie | AUPR | itérations externes (chemin) | temps (s) |
|---|---|---|---|---|---|
| `clean` (640 jeux) | legacy | 0.738 | 0.780 | 137 | 0.93 |
| | ve1_out5 | 0.738 | 0.780 | 135 | 0.97 |
| | full_out5 | 0.738 | 0.780 | 136 | 1.02 |
| | cg5_out5 | 0.738 | 0.780 | 135 | 1.08 |
| | ve1_out7 | 0.738 | 0.780 | 329 | 1.92 |
| | full_out7 | 0.738 | 0.780 | 311 | 2.02 |
| | cg5_out7 | 0.738 | 0.780 | 303 | 2.22 |
| `contaminated` (120 jeux) | legacy | 0.780 | 0.853 | 171 | 1.28 |
| | ve1_out5 | 0.777 | 0.853 | 206 | 1.03 |
| | full_out5 | 0.780 | 0.853 | 165 | 1.27 |
| | ve1_out7 | 0.779 | 0.853 | 1048 | 4.26 |

**Le réglage de l'optimisation ne change pas le réseau retrouvé.** Les écarts appariés de F1 avec
`legacy` restent sous 0.003, et ceux d'AUPR sous 0.0002. Faire converger la boucle alternée à 1e-7
double le temps sur `clean` et le triple sur `contaminated`, sans aucun gain.

Avec la pénalité sur l'échelle des corrélations, l'alternance cherche un point fixe : les poids
λ√(S_ii S_jj) dépendent de la covariance résiduelle S. Elle ne maximise donc pas un critère fixe. Sur
données réelles (`../../benchmark/pln_network/`), la poursuivre fait monter la pen_loglik
cumulée sur barents et la fait baisser sur mollusk, où la variance latente maximale dérive de 19 à 28.

On garde donc `maxit_ve = 1` (désormais effectif), `ftol_out = 1e-5` et `maxit_out = 20`.

```sh
Rscript sim_optim_settings.R clean 10 16 ve1_out5,ve1_out7,full_out5,full_out7,cg5_out5,cg5_out7
Rscript sim_optim_settings.R clean 10 16 legacy <bibliothèque où est installée l'ancienne version>
Rscript summarize.R
```
