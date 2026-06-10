---
name: feedback-cpp-audit
description: Corrections C++ appliquées sur la branche code-enhancement en juin 2026 — 8 fixes, compilation propre, 72 tests passent
metadata:
  type: feedback
---

8 améliorations C++ appliquées intégralement le 10/06/2026. Toutes les corrections sont sur la branche `code-enhancement`.

**Fixes appliqués :**

1. **Guard P_X (d=0)** — 4 fichiers `nlopt_*.cpp` : `(X.n_cols > 0) ? arma::solve(...) : arma::mat(0, Y.n_rows)` évite le crash quand la formule est `~ 0`.

2. **O(p³)→O(np) dans `nlopt_fixed_cov.cpp`** : `trace(Omega * (...))` remplacé par `full_cov_obj_grad_impl` (accu elementwise).

3. **`DenseOmegaImpl` base struct** dans `CovarianceTraits.h` : `FullCovTraits` et `FixedCovTraits` héritent désormais de `DenseOmegaImpl` qui contient les 6 méthodes statiques identiques (`cov_diag`, `grad_hess_M`, `times_Omega`, `penalty_M`, `objective_cov`, `final_loglik`).

4. **`NewtonConfig` struct** dans `utils.h` : centralise les 4 `containsElementNamed` parsings.

5. **Adoption `NewtonConfig`** : `newton_{full,diag,spherical,fixed}_cov.cpp` + `nlopt_full_cov.cpp`.

6. **`nlopt_impl.h`** : nouveau header partagé avec les 3 helpers `inline` ; copies `static` supprimées des `.cpp`.

7. **Méthode `update()`** dans chaque `State` de `CovarianceTraits.h` : constructeur et `mstep` délèguent vers `s.update(M, S2, w, w_bar)`.

8. **Style `SphericalCovTraits::output_cov`** : `.fill(s.sigma2)` à la place de `= arma::ones<arma::vec>(p) * s.sigma2`.

**Why:** factorisation, efficacité (O(p³)→O(np)), et robustesse (guard d=0).

**How to apply:** avant toute modification C++ future, relire `CovarianceTraits.h` pour comprendre la hiérarchie `DenseOmegaImpl → FullCovTraits/FixedCovTraits` et vérifier que les nouvelles méthodes communes sont ajoutées à la base, pas dupliquées.
