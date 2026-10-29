"""Vérifications calculatoires du chapitre « ar-l2-series-entieres » (analyse réelle, L2) :
cartes, exercices et QCM de content/*/analyse-reelle/l2-ar-l2-series-entieres.json et content/qcm/ar-l2-series-entieres.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/ar_l2_series_entieres.py (doit finir sans erreur, en quelques secondes).
"""

import cmath
import itertools
import math
from fractions import Fraction as F

import mpmath as mp
from sympy import (Rational, binomial, cos, exp, factorial, limit, log, oo, series, sin, sqrt, symbols, atan, asin,
                   erf, pi, simplify, diff, expand, Poly)

mp.mp.dps = 30
n = symbols("n", positive=True, integer=True)
x = symbols("x")


def coefs(expr, N):
    """Coefficients 0..N-1 du développement en série de expr en 0."""
    s = series(expr, x, 0, N).removeO()
    return [s.coeff(x, k) for k in range(N)]


def proche(a, b, eps=1e-20):
    return abs(mp.mpf(a) - mp.mpf(b)) < eps if not isinstance(a, complex) else abs(a - b) < 1e-12


# ---------------------------------------------------------------------------
# ar-l2-series-entieres-abel (carte)
# ---------------------------------------------------------------------------
# Leurres et pièges : a_n = 1, z0 = 2, z = 3/2 : (a_n) bornée, |z| < |z0|, mais (3/2)^n non borné.
assert F(3, 2) ** 60 > 10**10
# Leurre « série divergente en z0 » / « a_n z0^n ≠ 0 » : a_n = n!, z0 = 1, z = 1/2 : n!/2^n -> +inf.
assert math.factorial(40) / 2**40 > 10**30
# Contre-exemple : sum z^n/n converge en -1 (vers -ln 2), diverge en 1.
assert abs(mp.nsum(lambda k: (-1) ** k / k, [1, mp.inf]) + mp.log(2)) < mp.mpf(10) ** -25
# Affirmation 4 : a_n = (-1)^n/n ; en z = -1, a_n (-1)^n = 1/n (série harmonique divergente).
assert all(F((-1) ** k, k) * (-1) ** k == F(1, k) for k in range(1, 50))

# ar-l2-series-entieres-abel-ex1 (+ v1, v2) : a_n = (p + (-1)^n)^n, R = 1/(p+1).
for p, R, leurres in [(2, F(1, 3), [F(1, 2), 1, 3]), (3, F(1, 4), [F(1, 3), F(1, 2), 4]), (4, F(1, 5), [F(1, 4), F(1, 3), 5])]:
    a = [(p + (-1) ** k) ** k for k in range(60)]
    assert all(a[k] * R**k <= 1 for k in range(60))           # bornée en r = R
    r = R * F(11, 10)
    assert a[58] * r**58 > 100                                 # non bornée au-delà (rangs pairs)
    assert R not in leurres
# ar-l2-series-entieres-abel-ex2 (+ v1, v2) : exemples réalisant les deux bornes.
# R = 2 : a_n = 1/(n^2 2^n), terme en -3 de module (3/2)^n/n^2 non borné ; R = 3 : a_n = (-1)^n/(n 3^n), série en -3 = sum 1/n.
assert F(3, 2) ** 200 / 200**2 > 10**20
assert all(F((-1) ** k, k * 3**k) * (-3) ** k == F(1, k) for k in range(1, 30))
# v1 : a_n = 1/n^2 (R = 1, terme en 4i de module 4^n/n^2) ; a_n = 1/(n (4i)^n) (R = 4, série en 4i = sum 1/n).
assert all(abs(1 / (k * (4j) ** k) * (4j) ** k - 1 / k) < 1e-12 for k in range(1, 30))
assert abs(0.5j) < 1 and abs(-5) > 4
# v2 : a_n = 1/(n^2 5^n) (R = 5) ; a_n = 1/(n 6^n) (R = 6) ; points 4i (|4i| = 4 < 5) et -7 (7 > 6).
assert F(6, 5) ** 400 / 400**2 > 10**20 and abs(4j) < 5 and 7 > 6
# ar-l2-series-entieres-abel-ex3 (+ v1, v2) : sin n et cos n ne tendent pas vers 0 (|.| >= 1/2 infiniment souvent).
assert sum(1 for k in range(10**4, 10**4 + 200) if abs(math.sin(k)) >= 0.5) > 50
assert sum(1 for k in range(10**4, 10**4 + 200) if abs(math.cos(k)) >= 0.5) > 50
assert math.sin(1) != 0
# ar-l2-series-entieres-abel-ex4 : a_n = 1, b_n = 2^-n - 1 : a_n + b_n = 2^-n (rayon 2), b_n ~ -1 (rayon 1).
assert all(1 + (F(1, 2**k) - 1) == F(1, 2**k) for k in range(30))

# ---------------------------------------------------------------------------
# ar-l2-series-entieres-calcul-rayon (carte)
# ---------------------------------------------------------------------------
assert limit((n + 1) ** 2 / (3 * n**2), n, oo) == Rational(1, 3)                       # exemple 1 : R = 3
assert limit((2 * n + 2) * (2 * n + 1) / (n + 1) ** 2, n, oo) == 4                     # exemple 2 : R = 1/4
assert limit(((n + 1) / (n**3 + 2)) / (1 / n**2), n, oo) == 1                          # exemple 3 : a_n ~ 1/n^2
# Contre-exemple : a_n = 2 + (-1)^n, quotient alternativement 1/3 et 3.
assert {F(2 + (-1) ** (k + 1), 2 + (-1) ** k) for k in range(10)} == {F(1, 3), F(3)}
# Piège 2 : a_n = 2^n (n pair), 1 (n impair) : a_{2n+1}/a_{2n} = 4^-n -> 0, mais R = 1/2.
assert all(F(1, 2 ** (2 * k)) == F(1, 4**k) for k in range(10))
# Pièges et leurres avec a_n = 2^-n (R = 2) : quotient 1/2 <= 1, |a_n/a_{n+1}| = 2, |a_{n+1}| - |a_n| -> 0.
assert F(1, 2) <= 1 and 1 / F(1, 2) == 2 != 1
# Leurre « |a_n| -> l » : a_n = 1/(n+1) -> 0 mais R = 1 (pas +inf).
assert limit((n + 1) / (n + 2), n, oo) == 1
# Affirmation 3 : 2^n -> R = 1/2, 4^n -> R = 1/4.
assert limit(2 ** (n + 1) / 2**n, n, oo) == 2 and limit(4 ** (n + 1) / 4**n, n, oo) == 4

# ar-l2-series-entieres-calcul-rayon-ex1 (+ v1, v2).
assert limit(((n + 1) ** 2 + 1) / (2 * (n**2 + 1)), n, oo) == Rational(1, 2)            # R = 2
assert limit((n + 1) ** 3 / (5 * n**3), n, oo) == Rational(1, 5)                        # v1 : R = 5
assert limit(3 * (n + 3) / (n + 2), n, oo) == 3                                         # v2 : R = 1/3
# ar-l2-series-entieres-calcul-rayon-ex2 (+ v1, v2).
assert simplify((n + 1) ** 2 / ((2 * n + 2) * (2 * n + 1)) - (n + 1) / (2 * (2 * n + 1))) == 0
assert limit((n + 1) / (2 * (2 * n + 1)), n, oo) == Rational(1, 4)                      # R = 4
assert limit((n / (n + 1)) ** n, n, oo) == exp(-1)                                     # v1 : R = e
k = symbols("k", positive=True, integer=True)
ratio = binomial(3 * (k + 1), k + 1) / binomial(3 * k, k)
assert all(ratio.subs(k, j) == F((3 * j + 3) * (3 * j + 2) * (3 * j + 1), (j + 1) * (2 * j + 2) * (2 * j + 1)) for j in range(1, 8))
assert limit((3 * n + 3) * (3 * n + 2) * (3 * n + 1) / ((n + 1) * (2 * n + 2) * (2 * n + 1)), n, oo) == Rational(27, 4)  # v2 : R = 4/27
# ar-l2-series-entieres-calcul-rayon-ex3 (+ v1, v2) : séries lacunaires.
assert limit(3 * (n + 1) / (n + 2), n, oo) == 3          # sum 3^n/(n+1) w^n : R_w = 1/3, R = 1/sqrt 3
assert limit((n + 1) / (2 * n), n, oo) == Rational(1, 2)  # v1 : R_w = 2, R = 2^(1/3)
assert limit(5 * (n + 1) ** 2 / n**2, n, oo) == 5         # v2 : R_w = 1/5, R = 1/sqrt 5
# ar-l2-series-entieres-calcul-rayon-ex4 (+ v1, v2) : développement de n^2 ln(1 + c/n).
u = symbols("u", positive=True)   # u = 1/n
assert series(log(1 + u) / u**2, u, 0, 1).removeO() == 1 / u - Rational(1, 2)          # n - 1/2
assert series(log(1 + 2 * u) / u**2, u, 0, 1).removeO() == 2 / u - 2                   # v1 : 2n - 2
assert series(log(1 - u) / u**2, u, 0, 1).removeO() == -1 / u - Rational(1, 2)         # v2 : -n - 1/2
# Terme général au bord : a_n R^n -> e^-1/2, e^-2, e^-1/2 (vérification numérique à n = 10^6).
N = 10**6
assert abs(mp.exp(N**2 * mp.log(1 + mp.mpf(1) / N) - N) - mp.exp(-0.5)) < 1e-5
assert abs(mp.exp(N**2 * mp.log(1 + mp.mpf(2) / N) - 2 * N) - mp.exp(-2)) < 1e-5
assert abs(mp.exp(N**2 * mp.log(1 - mp.mpf(1) / N) + N) - mp.exp(-0.5)) < 1e-5

# ---------------------------------------------------------------------------
# ar-l2-series-entieres-convergence-normale (carte et exercices)
# ---------------------------------------------------------------------------
assert abs(mp.nsum(lambda j: mp.mpf(3) ** j / mp.factorial(j), [0, mp.inf]) - mp.e**3) < 1e-25   # exemple : sum r^n/n! = e^r
# ex1 (+ v1, v2) : sum des sup = 1/(1 - a/b).
for b_, a_, val, leurres in [(3, 2, F(3), [F(3, 2), F(2), F(1, 3)]), (4, 3, F(4), [F(4, 3), F(3), F(1, 4)]),
                             (5, 1, F(5, 4), [F(5), F(1, 5), F(4, 5)])]:
    assert 1 / (1 - F(a_, b_)) == val and val not in leurres
# ex2 : reste de sum x^n : x^(n+1)/(1-x) ; en x = 1 - 10^-6, R_10 > 10^5.
assert (1 - 1e-6) ** 11 / 1e-6 > 1e5
# ex3 (+ v1, v2) : intégration terme à terme.
assert abs(mp.nsum(lambda j: 1 / (j * mp.mpf(2) ** j), [1, mp.inf]) - mp.log(2)) < 1e-25
assert abs(mp.nsum(lambda j: 1 / (j * mp.mpf(3) ** j), [1, mp.inf]) - mp.log(1.5)) < 1e-25
assert abs(mp.nsum(lambda j: (-1) ** j / ((2 * j + 1) * mp.mpf(2) ** (2 * j + 1)), [0, mp.inf]) - mp.atan(0.5)) < 1e-25
for faux in [1, mp.log(1.5), 2 * mp.log(2)]:
    assert abs(faux - mp.log(2)) > 0.1
for faux in [mp.log(mp.mpf(4) / 3), 0.5, mp.log(3)]:
    assert abs(faux - mp.log(1.5)) > 0.05
for faux in [mp.pi / 8, mp.log(1.25), 0.5]:
    assert abs(faux - mp.atan(0.5)) > 0.01
# ex4 : sum 1/n diverge (sommes partielles > 10 pour n = 10^5).
assert mp.harmonic(10**5) > 10

# ---------------------------------------------------------------------------
# ar-l2-series-entieres-derivation (carte et exercices)
# ---------------------------------------------------------------------------
assert abs(mp.nsum(lambda j: j / mp.mpf(2) ** j, [1, mp.inf]) - 2) < 1e-25                 # exemple
# Affirmation 2 : S'(x) = -ln(1-x)/x tend vers +inf en 1^-.
assert -mp.log(1e-12) / (1 - 1e-12) > 27
# ex1 (+ v1, v2) : sum n x^n = x/(1-x)^2.
for xx, val, leurres in [(F(1, 3), F(3, 4), [F(3, 2), F(1, 2), F(9, 4)]), (F(1, 5), F(5, 16), [F(5, 4), F(1, 4), F(25, 16)]),
                         (F(-1, 2), F(-2, 9), [F(2, 9), F(-2, 3), F(2)])]:
    assert xx / (1 - xx) ** 2 == val and val not in leurres
    assert abs(mp.nsum(lambda j: j * mp.mpf(xx.numerator) ** j / mp.mpf(xx.denominator) ** j, [1, mp.inf]) - mp.mpf(val.numerator) / val.denominator) < 1e-20
# ex2 (+ v1, v2) : sum n^2 x^n = x(1+x)/(1-x)^3 ; sum n(n-1) x^n = 2x^2/(1-x)^3.
X = symbols("X")
assert simplify(2 * X**2 / (1 - X) ** 3 + X / (1 - X) ** 2 - X * (1 + X) / (1 - X) ** 3) == 0
assert F(1, 2) * F(3, 2) / F(1, 2) ** 3 == 6 and F(1, 3) * F(4, 3) / F(2, 3) ** 3 == F(3, 2) and 2 * F(1, 4) / F(1, 2) ** 3 == 4
assert abs(mp.nsum(lambda j: j**2 / mp.mpf(2) ** j, [1, mp.inf]) - 6) < 1e-20
assert abs(mp.nsum(lambda j: j**2 / mp.mpf(3) ** j, [1, mp.inf]) - 1.5) < 1e-20
assert abs(mp.nsum(lambda j: j * (j - 1) / mp.mpf(2) ** j, [2, mp.inf]) - 4) < 1e-20
# ex3 (+ v1, v2) : S(x) = 1 + (1-x) ln(1-x)/x.
for xx, val in [(mp.mpf(1) / 2, 1 - mp.log(2)), (mp.mpf(1) / 3, 1 - 2 * mp.log(1.5)), (mp.mpf(1) / 4, 1 - 3 * mp.log(mp.mpf(4) / 3))]:
    s = mp.nsum(lambda j: xx**j / (j * (j + 1)), [1, mp.inf])
    assert abs(s - val) < 1e-25 and abs(s - (1 + (1 - xx) * mp.log(1 - xx) / xx)) < 1e-25
assert min(abs(f - (1 - mp.log(2))) for f in [mp.log(2), 2 * mp.log(2) - 1, 1 - mp.log(2) / 2]) > 0.05
assert min(abs(f - (1 - 2 * mp.log(1.5))) for f in [mp.log(1.5), 1 - mp.log(1.5), 3 * mp.log(1.5) - 1]) > 0.02
assert min(abs(f - (1 - 3 * mp.log(mp.mpf(4) / 3))) for f in [mp.log(mp.mpf(4) / 3), 1 - mp.log(mp.mpf(4) / 3), 4 * mp.log(mp.mpf(4) / 3) - 1]) > 0.01
# ex4 : dilogarithme ; f(x) + f(-x) = f(x^2)/2 en x = 0,7 ; sum (-1)^(n-1)/n^2 = pi^2/12.
dilog = lambda t: mp.nsum(lambda j: t**j / j**2, [1, mp.inf])
assert abs(dilog(mp.mpf("0.7")) + dilog(mp.mpf("-0.7")) - dilog(mp.mpf("0.49")) / 2) < 1e-25
assert abs(mp.nsum(lambda j: (-1) ** (j - 1) / j**2, [1, mp.inf]) - mp.pi**2 / 12) < 1e-25

# ---------------------------------------------------------------------------
# ar-l2-series-entieres-dse (carte et exercices)
# ---------------------------------------------------------------------------
# Démonstration et ex4 : f^(n)(x) = P_n(1/x) e^(-1/x^2) avec P_{n+1}(X) = 2X^3 P_n(X) - X^2 P_n'(X).
P = Poly(1, X)
for _ in range(5):
    g = P.as_expr().subs(X, 1 / x) * exp(-1 / x**2)
    P = Poly(2 * X**3 * P.as_expr() - X**2 * diff(P.as_expr(), X), X)
    assert simplify(diff(g, x) - P.as_expr().subs(X, 1 / x) * exp(-1 / x**2)) == 0
# Contre-exemple g = sum e^-k cos(k^2 x) : (2n)! <= (2n)^(2n) et e^-2n (2n)^4n / (2n)! >= (2n/e)^(2n).
assert all(math.factorial(2 * m) <= (2 * m) ** (2 * m) for m in range(1, 30))
# Affirmation 4 : 1/(1+x^2) = sum (-1)^n x^(2n).
assert coefs(1 / (1 + x**2), 10) == [1, 0, -1, 0, 1, 0, -1, 0, 1, 0]
# ex1 (+ v1, v2) : coefficients des fractions rationnelles.
for expr, form, leurres in [
    (1 / ((1 - x) * (2 - x)), lambda m: 1 - F(1, 2 ** (m + 1)), [lambda m: 1 - F(1, 2**m), lambda m: F(1, 2 ** (m + 1)) - 1]),
    (1 / ((1 - x) * (3 - x)), lambda m: F(1, 2) * (1 - F(1, 3 ** (m + 1))), [lambda m: 1 - F(1, 3 ** (m + 1)), lambda m: F(1, 2) * (1 - F(1, 3**m))]),
    (1 / ((1 + x) * (2 - x)), lambda m: F(1, 3) * ((-1) ** m + F(1, 2 ** (m + 1))), [lambda m: F(1, 3) * (1 + F(1, 2 ** (m + 1))), lambda m: F(1, 3) * ((-1) ** m - F(1, 2 ** (m + 1)))]),
]:
    c = coefs(expr, 10)
    assert all(c[m] == form(m) for m in range(10))
    assert all(any(c[m] != lf(m) for m in range(10)) for lf in leurres)
# ex2 (+ v1, v2) : dérivées en 0.
assert coefs(x**3 * cos(x), 10)[9] * factorial(9) == -504
assert coefs(x**2 * exp(x), 9)[8] * factorial(8) == 56
assert coefs(x * sin(x), 9)[8] * factorial(8) == -8
assert diff(x**3 * cos(x), x, 9).subs(x, 0) == -504 and diff(x * sin(x), x, 8).subs(x, 0) == -8
# ex3 : e^(x^2) int_0^x e^(-t^2) dt = sqrt(pi)/2 e^(x^2) erf(x) ; coefficients 4^k k!/(2k+1)!.
c = coefs(sqrt(pi) / 2 * exp(x**2) * erf(x), 12)
assert all(c[2 * j] == 0 for j in range(6))
assert all(c[2 * j + 1] == Rational(4**j * math.factorial(j), math.factorial(2 * j + 1)) for j in range(6))
for lf in [lambda j: Rational(2**j * math.factorial(j), math.factorial(2 * j + 1)), lambda j: Rational(4**j, math.factorial(2 * j + 1)),
           lambda j: Rational((-1) ** j, math.factorial(j) * (2 * j + 1))]:
    assert any(c[2 * j + 1] != lf(j) for j in range(6))   # leurres faux (4^k/(2k+1)! coïncide pour k <= 1, pas pour k = 2)

# ---------------------------------------------------------------------------
# ar-l2-series-entieres-usuels (carte et exercices)
# ---------------------------------------------------------------------------
assert coefs(log(1 + x), 8) == [0] + [Rational((-1) ** (m - 1), m) for m in range(1, 8)]
assert coefs(atan(x), 10) == [0 if m % 2 == 0 else Rational((-1) ** ((m - 1) // 2), m) for m in range(10)]
assert coefs(sqrt(1 + x), 4) == [1, Rational(1, 2), Rational(-1, 8), Rational(1, 16)]
al = symbols("alpha")
cn = lambda m: Rational(1) if m == 0 else expand(math.prod([(al - j) for j in range(m)]) / factorial(m))
assert all(expand((m + 1) * cn(m + 1) - (al - m) * cn(m)) == 0 for m in range(6))
# Leurres du trou de ln(1+x) et arctan, évalués en x = 0,3.
t3 = mp.mpf("0.3")
assert abs(mp.log(1 + t3) - mp.nsum(lambda j: (-1) ** (j - 1) * t3**j / j, [1, mp.inf])) < 1e-25
assert min(abs(mp.log(1 + t3) - v) for v in [-mp.log(1 + t3), -mp.log(1 - t3), mp.exp(-t3)]) > 0.05
assert min(abs(mp.atan(t3) - v) for v in [mp.sin(t3), mp.atanh(t3), 1 / (1 + t3**2)]) > 1e-3
# ex1 (+ v1, v2).
assert coefs(log(1 + 2 * x), 6)[1:] == [Rational((-1) ** (m - 1) * 2**m, m) for m in range(1, 6)]
assert coefs(log(1 - 3 * x), 6)[1:] == [Rational(-(3**m), m) for m in range(1, 6)]
assert coefs(1 / (1 + 3 * x) ** 2, 6) == [(-1) ** m * (m + 1) * 3**m for m in range(6)]
# ex2 (+ v1, v2) : ln(1 + x + ... + x^(m-1)) : a_n = 1/n, ou -(m-1)/n si m | n.
for mm in (3, 4, 5):
    c = coefs(log(sum(x**j for j in range(mm))), 16)
    assert all(c[j] == (Rational(-(mm - 1), j) if j % mm == 0 else Rational(1, j)) for j in range(1, 16))
# ex3 : arcsin.
c = coefs(asin(x), 12)
assert all(c[2 * j + 1] == binomial(2 * j, j) / (4**j * (2 * j + 1)) for j in range(6))
assert c[3] == Rational(1, 6) and c[5] == Rational(3, 40)
assert coefs(1 / sqrt(1 - x**2), 12)[::2] == [binomial(2 * j, j) / 4**j for j in range(6)]
# ex4 : série de Leibniz, majoration 1/(2n+3) et nombre de termes (n = 499).
assert abs(mp.nsum(lambda j: (-1) ** j / (2 * j + 1), [0, mp.inf]) - mp.pi / 4) < 1e-25
assert 1 / F(2 * 499 + 3) <= F(1, 1000) < 1 / F(2 * 498 + 3)
P499 = mp.fsum((-1) ** j / mp.mpf(2 * j + 1) for j in range(500))
assert abs(P499 - mp.pi / 4) <= mp.mpf(1) / 1001

# ---------------------------------------------------------------------------
# ar-l2-series-entieres-exp-complexe (carte et exercices)
# ---------------------------------------------------------------------------
# Démonstration : cos 2 <= 1 - 2 + 2/3 = -1/3 ; u_n = 4^n/(2n)! décroissante pour n >= 1.
assert all(F(4**(j + 1), math.factorial(2 * j + 2)) < F(4**j, math.factorial(2 * j)) for j in range(1, 20))
assert math.cos(2) <= -1 / 3
assert abs(math.cos(math.pi / 2)) < 1e-15 and 1 < math.pi / 2 < 2
# Exemples et contre-exemples.
assert abs(cmath.exp(1j * math.pi) + 1) < 1e-15
assert abs(cmath.exp(math.log(2) + 1j * math.pi / 2) - 2j) < 1e-14
assert abs(abs(cmath.exp(1 + 1j)) - math.e) < 1e-14
assert abs(cmath.exp(1 + 1j * math.pi / 2) - 1j * math.e) < 1e-14
# Leurres du trou 1 : exp(ab) vs exp(a)exp(b) (a = b = 1) ; exp(0+0) = 1 != 2.
assert abs(math.exp(1) - math.exp(1) ** 2) > 1
# ex1 (+ v1, v2).
for z, w, leurres in [(math.log(2) + 1j * math.pi / 3, 1 + 1j * math.sqrt(3), [math.sqrt(3) + 1j, 2 + 1j * math.sqrt(3), 0.5 + 1j * math.sqrt(3) / 2]),
                      (math.log(3) - 1j * math.pi / 2, -3j, [3j, -3, 3 - 1j]),
                      (-math.log(2) + 3j * math.pi / 4, (-math.sqrt(2) + 1j * math.sqrt(2)) / 4,
                       [(math.sqrt(2) + 1j * math.sqrt(2)) / 4, -math.sqrt(2) + 1j * math.sqrt(2), (-math.sqrt(2) - 1j * math.sqrt(2)) / 4])]:
    assert abs(cmath.exp(z) - w) < 1e-14 and all(abs(cmath.exp(z) - l) > 0.1 for l in leurres)
# ex2 (+ v1, v2) : solutions de exp(z) = w ; les leurres échouent pour k = 0 ou k = 1.
for w, sol, leurres in [
    (1 + 1j, lambda k: math.log(2) / 2 + 1j * (math.pi / 4 + 2 * k * math.pi),
     [lambda k: math.log(2) + 1j * (math.pi / 4 + 2 * k * math.pi), lambda k: math.log(2) / 2 + 1j * (math.pi / 4 + k * math.pi)]),
    (-2, lambda k: math.log(2) + 1j * (math.pi + 2 * k * math.pi),
     [lambda k: math.log(2) + 2j * k * math.pi, lambda k: -math.log(2) + 1j * (math.pi + 2 * k * math.pi), lambda k: math.log(2) + 1j * (math.pi + k * math.pi)]),
    (-1 + 1j * math.sqrt(3), lambda k: math.log(2) + 1j * (2 * math.pi / 3 + 2 * k * math.pi),
     [lambda k: math.log(2) + 1j * (math.pi / 3 + 2 * k * math.pi), lambda k: math.log(2) + 1j * (-2 * math.pi / 3 + 2 * k * math.pi),
      lambda k: math.log(2) / 2 + 1j * (2 * math.pi / 3 + 2 * k * math.pi)]),
]:
    assert all(abs(cmath.exp(sol(k)) - w) < 1e-12 for k in range(-3, 4))
    assert all(any(abs(cmath.exp(l(k)) - w) > 0.1 for k in (0, 1)) for l in leurres)
# ex3 : cos(pi/4) = sin(pi/4) = sqrt(2)/2.
assert abs(cmath.exp(1j * math.pi / 4) - math.sqrt(2) / 2 * (1 + 1j)) < 1e-15
# ex4 : encadrement de pi.
assert 1 - F(9, 8) + F(81, 384) - F(729, 46080) == F(3231, 46080) > 0
assert 1 - 2 + F(16, 24) == F(-1, 3)
assert math.cos(1.5) >= 3231 / 46080 and 3 < math.pi < 4

# ---------------------------------------------------------------------------
# ar-l2-series-entieres-generatrices (carte et exercices)
# ---------------------------------------------------------------------------
fib = [0, 1]
while len(fib) < 20:
    fib.append(fib[-1] + fib[-2])
assert coefs(x / (1 - x - x**2), 20) == fib                                            # exemple : Fibonacci
phi, psi = (1 + math.sqrt(5)) / 2, (1 - math.sqrt(5)) / 2
assert all(abs((phi**m - psi**m) / math.sqrt(5) - fib[m]) < 1e-6 for m in range(20))
assert all(fib[m] <= 2**m for m in range(20))
# Contre-exemple : |c_n| >= 2(n+1)/(n+2) >= 1 pour le produit de Cauchy de (-1)^n/sqrt(n+1) par elle-même.
for m in range(60):
    cm = sum(1 / math.sqrt((j + 1) * (m - j + 1)) for j in range(m + 1))
    assert cm >= 2 * (m + 1) / (m + 2) - 1e-12 >= 1 - 1e-12
# Loi P(X = n) = 1/(n(n+1)) : somme 1 (télescopage).
assert sum(F(1, j * (j + 1)) for j in range(1, 200)) == 1 - F(1, 200)
# Leurre du trou 1 « sum a_k b_k » : A = x, b_n = 2^n : x/(1-2x) a pour coefficient de x 1, pas 2.
assert coefs(x / (1 - 2 * x), 3)[1] == 1 != 2
# Affirmation 4 : nombre de couples (a, b) avec a + b = n : n + 1 = coefficient de 1/(1-x)^2.
assert coefs(1 / (1 - x) ** 2, 8) == [m + 1 for m in range(8)] and coefs(1 / (1 - x**2), 2)[1] == 0
# ex1 (+ v1, v2) : façons de payer.
paye = lambda p1, p2, s: sum(1 for a in range(s + 1) for b in range(s + 1) if p1 * a + p2 * b == s)
assert paye(1, 2, 10) == 6 and paye(1, 3, 10) == 4 and paye(2, 3, 12) == 3
assert all(paye(1, 2, s) == F(2 * s + 3 + (-1) ** s, 4) for s in range(30))
assert coefs(1 / ((1 - x) * (1 - x**2)), 11)[10] == 6 and coefs(1 / ((1 - x**2) * (1 - x**3)), 13)[12] == 3
# ex2 (+ v1, v2) : somme de variables de Poisson indépendantes (convolution exacte).
for l1, l2, kk, val, leurres in [
    (mp.mpf(1), mp.mpf(2), 2, mp.mpf(9) / 2 * mp.exp(-3), [9 * mp.exp(-3), mp.mpf(5) / 2 * mp.exp(-3), 3 * mp.exp(-3)]),
    (mp.mpf(1) / 2, mp.mpf(3) / 2, 3, mp.mpf(4) / 3 * mp.exp(-2), [8 * mp.exp(-2), 2 * mp.exp(-2), mp.mpf(7) / 12 * mp.exp(-2)]),
    (mp.mpf(1), mp.mpf(1), 1, 2 * mp.exp(-2), [mp.exp(-2), mp.exp(-1), 4 * mp.exp(-2)]),
]:
    pois = lambda lam, j: mp.exp(-lam) * lam**j / mp.factorial(j)
    conv = mp.fsum(pois(l1, j) * pois(l2, kk - j) for j in range(kk + 1))
    assert abs(conv - val) < 1e-25 and abs(conv - pois(l1 + l2, kk)) < 1e-25
    assert all(abs(l - val) > 1e-3 for l in leurres)
# ex3 (+ v1, v2) : récurrences a_{n+1} = c a_n + b^n, a_0 = 0.
for c_, b_, form, leurres in [
    (2, 3, lambda m: 3**m - 2**m, [lambda m: 3**m - 1, lambda m: m * F(3) ** (m - 1), lambda m: F(3**m - 1, 2)]),
    (4, 2, lambda m: F(4**m - 2**m, 2), [lambda m: 4**m - 2**m, lambda m: F(4**m - 1, 3), lambda m: m * F(4) ** (m - 1)]),
    (1, 2, lambda m: 2**m - 1, [lambda m: F(2) ** (m - 1), lambda m: m * F(2) ** (m - 1), lambda m: 2**m - m]),
]:
    a = [0]
    for m in range(20):
        a.append(c_ * a[-1] + b_**m)
    assert all(a[m] == form(m) for m in range(21))
    assert all(any(a[m] != lf(m) for m in range(21)) for lf in leurres)
assert all(math.prod([1]) * 0 <= 3**m - 2**m <= 3**m for m in range(20))
assert simplify(x / ((1 - 2 * x) * (1 - 3 * x)) - (1 / (1 - 3 * x) - 1 / (1 - 2 * x))) == 0
assert simplify(x / ((1 - 4 * x) * (1 - 2 * x)) - (1 / (1 - 4 * x) - 1 / (1 - 2 * x)) / 2) == 0
assert simplify(x / ((1 - x) * (1 - 2 * x)) - (1 / (1 - 2 * x) - 1 / (1 - x))) == 0
# ex4 : nombres de Catalan, par énumération des mots de Dyck (n <= 7).
def dyck(m):
    total = 0
    for mot in itertools.product((1, -1), repeat=2 * m):
        h = 0
        for pas in mot:
            h += pas
            if h < 0:
                break
        else:
            total += h == 0
    return total
cat = [dyck(m) for m in range(8)]
assert cat == [math.comb(2 * m, m) // (m + 1) for m in range(8)] and cat[5] == 42
assert all(cat[m + 1] == sum(cat[j] * cat[m - j] for j in range(m + 1)) for m in range(7))
assert coefs(sqrt(1 - 4 * x), 10)[1:] == [-Rational(2, m) * binomial(2 * m - 2, m - 1) for m in range(1, 10)]
assert coefs((1 - sqrt(1 - 4 * x)) / (2 * x), 10) == [binomial(2 * m, m) / (m + 1) for m in range(10)]
assert all(any(cat[m] != v(m) for m in range(8)) for v in [lambda m: math.comb(2 * m, m), lambda m: F(4**m, m + 1), lambda m: F(math.comb(2 * m, m), 2 * m + 1)])

# ---------------------------------------------------------------------------
# QCM ar-l2-series-entieres
# ---------------------------------------------------------------------------
assert limit((2 ** (n + 1) + 3 ** (n + 1)) / (2**n + 3**n), n, oo) == 3                # q2 : R = 1/3
assert 1 / (1 - F(1, 3)) ** 2 == F(9, 4)                                                # q4
assert Rational(1, 2) * Rational(-1, 2) * Rational(-3, 2) / 6 == Rational(1, 16)        # q6
assert abs(cmath.exp(1j * math.pi) + 1) < 1e-15 and abs(cmath.exp(0) - 1) < 1e-15     # q7 (leurre k = -1 : z = 0)
assert math.exp(math.pi) > 0 and abs(cmath.exp(1j * math.pi / 2) + 1) > 1              # q7 (autres leurres)
G = x / (2 - x)
assert diff(G, x).subs(x, 1) == 2 and coefs(G, 6) == [0] + [Rational(1, 2**m) for m in range(1, 6)]   # q8
print("ar_l2_series_entieres : toutes les vérifications passent.")
