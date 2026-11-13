# Point de départ de ZIPLN

Banc d'essai du 13/11/2026 (branche `zipln-init`). La relecture de la famille variationnelle
de ZIPLN (branche `zipln-j3`, `inst/simus_ZIPLN/variational_J3/`) a montré que, sur données
simulées à forte inflation et espèces abondantes, `ZIPLN()` s'arrête sur un mauvais optimum
local depuis son point de départ.

## Fichiers

- `simulate.R` : simulations ZIPLN (n = 100, p = 20 ; ρ ∈ {0 ; 0,5 ; 0,8} ; inflation
  faible/forte ; abondance basse/élevée).
- `inits.R` : points de départ candidats (`current` = régression linéaire sur log(1 + Y) ;
  `positive` = B et variances latentes sur les comptages positifs, M = XB sur les zéros,
  R par un Poisson-lognormal calé sur le taux de zéros ; `zip` = ZIP univarié par EM ;
  `truth`), et ajustement de `ZIPLN()` avec un point de départ imposé.
- `check_oracle.R` : 6 réplicats à forte inflation et forte abondance, backends builtin et
  nlopt.
- `benchmark.R` : `sim` (12 configurations × 10 réplicats) et `real` (trichoptera,
  barents, mollusk, oaks, microcosm).

## Résultats

- **Forte abondance** : le départ actuel termine 100 à 200 sous l'optimum atteint depuis les
  vrais paramètres, avec builtin comme avec nlopt ; `positive` (et `zip`) l'atteint
  exactement. Erreur sur π divisée par 2,5 à 4,5, sur les variances par 2,4 à 2,8.
- **Faible abondance** : le départ actuel trouve une ELBO un peu plus haute, même que les vrais
  paramètres. L'ELBO y préfère expliquer les zéros par la variance latente ; aucune
  initialisation n'y change rien.
- **Données réelles** : `positive` gagne 705 sur microcosm (inflation selon le site) et 25 sur
  barents (inflation par espèce), mais perd 26 à 780 sur oaks et 37 à 93 sur mollusk. `zip`
  plante sur oaks et microcosm.

Aucun départ ne domine : `ZIPLN()` part désormais des deux (`init = "best"`) et garde la
meilleure borne. Le temps de calcul est environ doublé ; les temps de ce banc d'essai ont été
mesurés avec un build de débogage et ne valent qu'en relatif.

## Gagner du temps (13/11/2026)

`trajectories.R` enregistre la trajectoire de l'objectif depuis chaque départ, sur les 120
simulations et 13 jeux réels.

- **Course précoce** (`analyze_trajectories.R`) : ajuster les deux départs k itérations et
  garder le meilleur. Inutilisable : 40 décisions fausses sur 133 à k = 5, encore 17 à k = 20,
  avec des pertes jusqu'à 182 (microcosm) ; l'ordre s'inverse souvent tard.
- **Choix par espèce sur le modèle diagonal** (`species_start.R`) : l'ELBO diagonale est
  séparable par espèce. Comme raccourci, il se trompe sur oaks et mollusk (le modèle diagonal
  ignore les corrélations) et n'est pas moins cher sur de petits p. Comme troisième candidat,
  il gagne sur 6 jeux réels sur 13 (microcosm +177), mais pas sur les simulations : piste
  gardée pour plus tard.
- **Abandon du second départ** (`analyze_extrapolation.R`, puis `validate_abandon.R`) :
  ajuster `"lm"` jusqu'au bout, puis `"positive"` en l'abandonnant dès que l'extrapolation
  géométrique de sa trajectoire (fenêtre 10, marge ×3, après 30 itérations) ne peut plus
  atteindre la borne de `"lm"`. Sur les 133 cas : résultat identique au meilleur des deux
  ajustements complets, second départ abandonné dans 52 cas, environ un quart des itérations
  économisées (rien quand `"positive"` gagne). **Retenu.**
- **Ordonner les départs selon un indicateur du scénario** (`order_features.R`,
  `order_evaluate.R`) : mettre en premier le départ prédit gagnant, pour que l'autre soit
  abandonné. Meilleur indicateur, calculé sans ajustement : p0, la probabilité moyenne d'un
  zéro sous le Poisson-lognormal ajusté aux comptages positifs (« positive » en premier si
  p0 < 0,10). Il prédit le gagnant dans 85 % des cas (11 sur 13 sur données réelles), sans
  jamais changer le résultat. Mais le gain est minime : 25,5 % d'itérations économisées contre
  23,7 % avec « lm » toujours en premier, et 31 % au mieux avec un ordre parfait. La limite
  est ailleurs : le départ perdant tourne une trentaine d'itérations (le minimum de
  `cannot_beat()`) avant d'être jugé, et l'abaisser fait des erreurs. **Non retenu** ; p0
  reste un bon descripteur du scénario (zéros bien plus nombreux que ce qu'explique la
  partie Poisson-lognormal → « positive » gagne).
