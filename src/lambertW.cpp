/* lambertW.cpp

Copyright (C) 2015, Avraham Adler
All rights reserved.

SPDX-License-Identifier: BSD-2-Clause

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

References:

Corless, R. M.; Gonnet, G. H.; Hare, D. E.; Jeffrey, D. J. & Knuth, D. E.
 "On the Lambert W function", Advances in Computational Mathematics,
 Springer, 1996, 5, 329-359

Fritsch, F. N.; Shafer, R. E. & Crowley, W. P.
 "Solution of the transcendental equation (we^w = x)",
 Communications of the ACM, Association for Computing Machinery (ACM),
 1973, 16, 123-124
*/

#include "lambertW.h"

  /* Fritsch Iteration
  * W_{n+1} = W_n * (1 + e_n)
  * z_n = ln(x / W_n) - W_n
  * q_n = 2 * (1 + W_n) * (1 + W_n + 2 / 3 * z_n)
  * e_n = z_n / (1 + W_n) * (q_n - z_n) / (q_n - 2 * z_n)
  */

double FritschIter(double x, double w){
  int MaxEval = 5;
  bool CONVERGED = false;
  double k = 2.0 / 3.0;
  int i = 0;
  do {
    double z = std::log(x / w) - w;
    double w1 = w + 1.0;
    double q = 2.0 * w1 * (w1 + k * z);
    double qmz = q - z;
    double e = z / w1 * qmz / (qmz - z);
    CONVERGED = std::abs(e) <= EPS;
    w *= (1.0 + e);
    ++i;
  } while (!CONVERGED && i < MaxEval);
  return(w);
}

double lambertW0_CS(double x) {
  if (x == R_PosInf) {
    return(R_PosInf);
  } else if (x < -M_1_E) {
    return(R_NaN);
  } else if (std::abs(x + M_1_E) <= EPS) {
    return(-1.0);
  } else if (std::abs(x) <= 1e-16) {
    /* This close to 0 the W_0 branch is best estimated by its Taylor/Pade
     expansion whose first term is the value x and remaining terms are below
     machine double precision. See
     https://math.stackexchange.com/questions/1700919
     */
    return(x);
  } else {
    double w;
    if (std::abs(x) <= 6.4e-3) {
      /* When this close to 0 the Fritsch iteration may underflow. Instead,
       * function will use degree-6 minimax polynomial approximation of Halley
       * iteration-based values. Should be more accurate by three orders of
       * magnitude than Fritsch's equation (5) in this range.
       */

      // Minimax Approximation calculated using R package minimaxApprox 0.1.0
      return((((((-1.0805085529250425e1 * x + 5.2100070265741278) * x -
             2.6666665063383532) * x + 1.4999999657268301) * x -
             1.0000000000016802) * x + 1.0000000000001752) * x +
             2.6020852139652106e-18);

    } else if (x <= M_E) {
      /* Use expansion in Corliss 4.22 to create (2, 2) Pade approximant.
       * Equation with a few extra terms is:
       * -1 + p - 1/3p^2 + 11/72p^3 - 43/540p^4 + 689453/8398080p^4 - O(p^5)
       * This is just used to estimate a good starting point for the Fritsch
       * iteration process itself.
      */
      double p = std::sqrt(2.0 * (M_E * x + 1.0));
      double Numer = (0.2787037037037037 * p + 0.311111111111111) * p - 1.0;
      double Denom = (0.0768518518518518 * p + 0.688888888888889) * p + 1.0;
      w = Numer / Denom;
    } else {
      /* Use first five terms of Corliss et al. 4.19 */
      w = std::log(x);
      double L_2 = std::log(w);
      double L_3 = L_2 / w;
      double L_3_sq = L_3 * L_3;
      w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
        w + L_3_sq * L_3 / 3.0;
    }
    return(FritschIter(x, w));
  }
}

double lambertWm1_CS(double x){
  if (x == 0.0) {
    return(R_NegInf);
  } else if (x < -M_1_E || x > 0.0) {
    return(R_NaN);
  } else if (std::abs(x + M_1_E) <= EPS) {
    return(-1.0);
  } else {
    double w;
    /* Use first five terms of Corliss et al. 4.19 */
    w = std::log(-x);
    double L_2 = std::log(-w);
    double L_3 = L_2 / w;
    double L_3_sq = L_3 * L_3;
    w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
      w + L_3_sq * L_3 / 3.0;
    return(FritschIter(x, w));
  }
}

