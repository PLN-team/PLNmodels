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
