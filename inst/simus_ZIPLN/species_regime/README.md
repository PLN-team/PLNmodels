# Régime des zéros espèce par espèce dans ZIPLN

Suite de `inst/simus_ZIPLN/init/` (PR #199). Les départs « lm » et « positive » sont les deux
extrêmes d'un choix fait espèce par espèce : les zéros d'une espèce sont expliqués soit par
la variance latente, soit par l'inflation. Un départ mixte, choisi par espèce sur des
ajustements diagonaux, battait les deux de 177 sur microcosm (inflation par espèce).

## Fichiers

- `cases.R` : les 13 cas réels du banc d'essai.
- `diagnose.R` : profil des variances latentes et des taux d'inflation de l'ajustement retenu.
- `flips.R` : bascule de régime après convergence, puis réajustement depuis l'état mixte.

## Résultats (13/11/2026)

- **Pas d'espèce dégénérée dans ZIPLN** : aucune variance latente au-dessus de 100 (au plus
  21), 0 à 3 espèces au-dessus de 5 fois la médiane. Le régime « variance » est diffus, pas
  concentré sur quelques espèces explosives : le lien avec `degenerate_species` ne tient pas
  pour ZIPLN.
- **Indicateur de bascule** : écart entre la part des zéros en excès sous un Poisson-lognormal
  ajusté aux comptages positifs, et la part attribuée à l'inflation par l'ajustement. Il
  signale presque toutes les espèces (écart médian 0,5 à 0,9) : ajuster sur les seuls
  comptages positifs ignore la troncature en zéro, ce qui sous-estime la probabilité d'un zéro
  et fait voir de l'excès partout.
- **Mais une bascule bien choisie paie beaucoup** : sur microcosm (inflation par espèce),
  basculer une seule espèce vers l'inflation gagne **+151** ; 9 espèces +101, 31 espèces +25.
  Sur barents (inflation par espèce) +8. Basculer en bloc les espèces signalées est en revanche
  désastreux sur oaks (−772) et mollusk (−93).

Conclusion d'étape : la bascule est un bon mouvement, à condition de la faire espèce par espèce
et de la refuser quand elle fait baisser la borne. Prochaine étape : recherche gloutonne.

## Recherche gloutonne (`greedy.R`, `control_noop.R`, `greedy_converged.R`)

Bascule d'une espèce à la fois (10 candidates vers l'inflation, 3 vers la variance, dans
l'ordre de l'indicateur), réajustement, bascule gardée seulement si la borne monte.

- **Piège de convergence** : avec le plafond par défaut `maxit_out = 200`, 6 cas sur 13 ne
  convergent pas. Un simple réajustement sans bascule gagne alors jusqu'à +167 (microcosm,
  inflation par espèce), ce qui gonflait les gains de la première recherche (+251). À
  convergence (jusqu'à 560 itérations), ce témoin ne gagne plus que 0,4 au plus.
- **Gains propres des bascules, à convergence** : microcosm (par espèce) +78, microcosm (site)
  +24, oaks (`tree | tree`) +10,6, barents (par espèce) +7,9, oaks (par espèce) +5,4, barents
  (unique) +3,8, et +0,1 à +2,3 ailleurs. Jamais négatifs par construction.
- **Effet sur les estimations** : quasi nul en moyenne (π moyen et variance latente médiane
  inchangés à 2 % près) ; les bascules ne touchent que quelques espèces.
- **Coût** : environ la moitié d'un ajustement sur microcosm (53 s contre 112 s), 1 à 3 fois
  l'ajustement sur les petits jeux. L'indicateur ordonne mal les candidates : 76 bascules
  gardées sur 154 essais, et des essais refusés qui coûtent jusqu'à −160.

Conclusions d'étape : la bascule trouve de meilleurs optima, mais les gains restent modestes
au regard de la borne et des estimations ; le plafond de convergence de ZIPLN pèse plus lourd.
