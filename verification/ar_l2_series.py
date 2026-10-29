"""Vérifications calculatoires du chapitre « ar-l2-series » (analyse réelle, L2, séries numériques) :
cartes, exercices et QCM de content/*/analyse-reelle/l2-ar-l2-series.json et content/qcm/ar-l2-series.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/ar_l2_series.py (doit finir sans erreur, en quelques secondes).
"""

import cmath
import math
from fractions import Fraction as F

import mpmath as mp
from sympy import (Rational, binomial, cos, exp, integrate, limit, log, oo, pi,
                   simplify, sin, sqrt, symbols, zeta)

mp.mp.dps = 30
x, t = symbols("x t", positive=True)
n, k = symbols("n k", positive=True, integer=True)


def eq_lim(expr, cible):
    """expr / cible -> 1 quand n -> +oo."""
    return simplify(limit(expr / cible, n, oo) - 1) == 0


# ---------------------------------------------------------------------------
# ar-l2-series-positives (carte)
# ---------------------------------------------------------------------------
# Exemple : ln(1+1/n^2) ~ 1/n^2 et sin(1/n) ~ 1/n.
assert eq_lim(log(1 + 1 / n**2), 1 / n**2) and eq_lim(sin(1 / n), 1 / n)
# Contre-exemple 1 : ((-1)^n/sqrt n + 1/n) / ((-1)^n/sqrt n) = 1 + (-1)^n/sqrt n -> 1.
assert all(abs((((-1) ** m / math.sqrt(m) + 1 / m) / ((-1) ** m / math.sqrt(m))) - 1 - (-1) ** m / math.sqrt(m)) < 1e-12
           for m in range(1, 50))
# Contre-exemple 3 : H_n >= ln(n+1).
H = 0.0
for m in range(1, 2000):
    H += 1 / m
    assert H >= math.log(m + 1)
# Affirmation 3 : u_n = 1/n aux carrés, 1/n^2 sinon ; sommes partielles majorées par 2 zeta(2), n u_n = 1 aux carrés.
s = sum(F(1, m) if int(math.isqrt(m)) ** 2 == m else F(1, m * m) for m in range(1, 3000))
assert s < 2 * math.pi**2 / 6
# Trou 3, leurre u_n = O(v_n) et u_n = o(v_n) : u = 1/n^2, v = 1/n ; leurre u - v -> 0 : u = 1/n, v = 1/n^2.
assert limit((1 / n**2) / (1 / n), n, oo) == 0 and limit(1 / n - 1 / n**2, n, oo) == 0

# ar-l2-series-positives-ex1 (+ v1, v2) : équivalents rationnels.
assert eq_lim((2 * n + 1) / (n**3 + n + 1), 2 / n**2) and not eq_lim((2 * n + 1) / (n**3 + n + 1), 1 / n**2)
assert eq_lim((n**2 + 3) / (n**3 + 5 * n), 1 / n) and not eq_lim((n**2 + 3) / (n**3 + 5 * n), Rational(3, 5) / n)
assert eq_lim((sqrt(n) + 1) / (n**2 + 1), n ** Rational(-3, 2)) and not eq_lim((sqrt(n) + 1) / (n**2 + 1), 1 / n**2)
# ar-l2-series-positives-ex2 (+ v1, v2) : 1 <= 2 + cos n <= 3 ; ln n <= n ; ln 3 > 1 (minoration pour n >= 3).
assert all(1 <= 2 + math.cos(m) <= 3 for m in range(1, 1000))
assert all(math.log(m) <= m for m in range(1, 1000)) and math.log(3) > 1
# Le rapport (2 + cos n)/2 ne tend pas vers 1 (leurre « équivalent à 2/n^{3/2} »).
assert max(2 + math.cos(m) for m in range(1, 100)) > 2.9 and min(2 + math.cos(m) for m in range(1, 100)) < 1.1
# ar-l2-series-positives-ex3 (+ v1, v2) : développements limités.
assert eq_lim(exp(1 / n) - 1 - 1 / n, 1 / (2 * n**2))
assert eq_lim(1 / n - sin(1 / n), 1 / (6 * n**3))
assert eq_lim(cos(1 / sqrt(n)) - 1, -1 / (2 * n))
# ar-l2-series-positives-ex4 : Bertrand, primitives -1/ln t et ln ln t ; sqrt(u_n/n) = 1/(n ln n).
assert simplify(integrate(1 / (t * log(t) ** 2), t) + 1 / log(t)) == 0
assert simplify(integrate(1 / (t * log(t)), t) - log(log(t))) == 0
assert abs(math.sqrt((1 / (7 * math.log(7) ** 2)) / 7) - 1 / (7 * math.log(7))) < 1e-15

# ---------------------------------------------------------------------------
# ar-l2-series-riemann (carte)
# ---------------------------------------------------------------------------
# Exemple : sum ((1+i)/2)^n = 1 + i ; |(1+i)/2| = sqrt(2)/2.
q = (1 + 1j) / 2
assert abs(1 / (1 - q) - (1 + 1j)) < 1e-15 and abs(abs(q) - math.sqrt(2) / 2) < 1e-15
# Démonstration : S_n <= 1 + 1/(alpha-1) = alpha/(alpha-1) (affirmation 3), testé pour quelques alpha.
for a in [1.5, 2, 3]:
    assert sum(1 / m**a for m in range(1, 20000)) <= a / (a - 1)
# Affirmation 1 : n^{1/n} -> 1, donc 1/n^{1+1/n} ~ 1/n.
assert limit(n ** (1 / n), n, oo) == 1
# Affirmation 4 : 1/n, rapport n/(n+1) < 1.
# Trous : leurres 1/(1+q), q/(1-q), 1/(1-|q|) différents de 1/(1-q) pour q = 1/2 ou -1/2.
assert F(1, 1) / (1 - F(1, 2)) == 2 and 1 / (1 + F(1, 2)) != 2 and F(1, 2) / (1 - F(1, 2)) != 2
assert 1 / (1 - F(1, 2)) != 1 / (1 + F(1, 2))  # q = -1/2 : 1/(1-|q|) = 2 et 1/(1-q) = 2/3

# ar-l2-series-riemann-ex1 (+ v1, v2) : sommes géométriques et leurres (départs en 0, 1, premier terme).
def geom(a, r, n0):
    return a * r**n0 / (1 - r)


assert geom(3, F(1, 4), 2) == F(1, 4) and [geom(3, F(1, 4), 0), geom(3, F(1, 4), 1), 3 * F(1, 16)] == [4, 1, F(3, 16)]
assert geom(2, F(-1, 3), 1) == F(-1, 2) and geom(2, F(-1, 3), 0) == F(3, 2) and geom(2, F(1, 3), 1) == 1
assert geom(1, F(2, 3), 3) == F(8, 9) and [geom(1, F(2, 3), m) for m in (0, 1, 2)] == [3, 2, F(4, 3)]
# ar-l2-series-riemann-ex2 (+ v1, v2) : exposants.
assert eq_lim(n ** Rational(1, 3) / sqrt(n**3 + 1), n ** (Rational(1, 3) - Rational(3, 2)))
assert eq_lim(n ** 2 / (n**3 + n + 1), n ** (2 - 3))
assert eq_lim(n / (n**2 + n) ** Rational(1, 3), n ** (1 - Rational(2, 3)))
for a, beta in [(Rational(1, 2), Rational(3, 2) - Rational(1, 2)), (Rational(-1, 3), Rational(2, 3) + Rational(1, 3)),
                (1, 3 - 2)]:
    assert beta == 1  # frontière : a = 1/2, -1/3, 1 donne un exposant 1 (divergence)
# ar-l2-series-riemann-ex3 (+ v1, v2) : sum cos(n theta)/2^n = Re 1/(1 - e^{i theta}/2) et leurres.
for th, rep, leurres in [(math.pi / 3, 1, [2, math.sqrt(3) / 3, 4 / 3]),
                         (math.pi / 2, 4 / 5, [2, 2 / 5, 1]),
                         (2 * math.pi / 3, 5 / 7, [2, math.sqrt(3) / 7, 4 / 5])]:
    z = 1 / (1 - cmath.exp(1j * th) / 2)
    sdir = sum(math.cos(m * th) / 2**m for m in range(200))
    assert abs(z.real - rep) < 1e-12 and abs(sdir - rep) < 1e-12
    assert abs(z.imag - leurres[1]) < 1e-12 and abs(1 / (1 - math.cos(th) / 2) - leurres[2]) < 1e-12
    assert all(abs(l - rep) > 1e-3 for l in leurres)
# ar-l2-series-riemann-ex4 (+ v1, v2) : plus petit N avec 1/((alpha-1) N^{alpha-1}) <= eps.
def plus_petit(f, eps):
    m = 1
    while f(m) > eps:
        m += 1
    return m


assert plus_petit(lambda m: F(1, m), F(1, 1000)) == 1000
assert plus_petit(lambda m: F(1, 2 * m * m), F(1, 10**4)) == 71 and 70**2 == 4900 and 71**2 == 5041
assert plus_petit(lambda m: 2 / math.sqrt(m), 1e-2) == 40000
# Leurres : 999, 32, 500 ; 70, 100, 50 ; 20000, 200, 10000 ne conviennent pas (bien que 999 < 1000, etc.).
assert 1 / 999 > 1e-3 and 1 / 2 / 70**2 > 1e-4 and 2 / math.sqrt(20000) > 1e-2
# Majoration (question 1) : 1/n^alpha <= (1/(alpha-1)) (1/(n-1)^{alpha-1} - 1/n^{alpha-1}).
for a in [1.5, 2, 3]:
    assert all(1 / m**a <= (1 / (m - 1) ** (a - 1) - 1 / m ** (a - 1)) / (a - 1) for m in range(2, 500))

# ---------------------------------------------------------------------------
# ar-l2-series-sommation-comparaison (carte)
# ---------------------------------------------------------------------------
# Exemple : sum_{k>n} (1/(k-1) - 1/k) = 1/n ; 1/k^2 ~ 1/(k(k-1)).
assert eq_lim(1 / n**2, 1 / (n * (n - 1)))
# Cours : 2(sqrt(k+1) - sqrt k) ~ 1/sqrt k.
assert eq_lim(2 * (sqrt(n + 1) - sqrt(n)), 1 / sqrt(n))
# Contre-exemple 2 : sum 1/n^2 + 1/n^3 = zeta(2) + zeta(3) != zeta(2).
assert zeta(3) > 1
# Affirmation 2 : sum 1/((n+1)(n+2)) = 1 et sum 1/(n+1)^2 = pi^2/6.
assert sum(F(1, (m + 1) * (m + 2)) for m in range(100)) == 1 - F(1, 101)
# Affirmation 4 : ln ln(n+1) - ln ln n ~ 1/(n ln n).
assert eq_lim(log(log(n + 1)) - log(log(n)), 1 / (n * log(n)))

# ar-l2-series-sommation-comparaison-ex1 (+ v1, v2) : restes.
assert eq_lim(1 / (n - 1) ** 2 - 1 / n**2, 2 / n**3)
assert eq_lim(1 / (n - 1) ** 3 - 1 / n**3, 3 / n**4)
assert eq_lim(1 / sqrt(n - 1) - 1 / sqrt(n), 1 / (2 * n ** Rational(3, 2)))
for p, cible in [(3, lambda m: 1 / (2 * m * m)), (4, lambda m: 1 / (3 * m**3)), (1.5, lambda m: 2 / math.sqrt(m))]:
    m0 = 2000
    reste = float(mp.zeta(p, m0 + 1))  # zêta de Hurwitz : somme des 1/j^p pour j > m0
    assert abs(reste / cible(m0) - 1) < 1e-2
# ar-l2-series-sommation-comparaison-ex2 (+ v1, v2) : sommes partielles.
for f, cible in [(math.sqrt, lambda m: 2 / 3 * m**1.5), (lambda j: 1 / math.sqrt(j), lambda m: 2 * math.sqrt(m))]:
    m0 = 10**6
    sp = math.fsum(f(j) for j in range(1, m0 + 1))
    assert abs(sp / cible(m0) - 1) < 1e-2
assert eq_lim((n + 1) * log(n + 1) - n * log(n), log(n))
assert abs(math.lgamma(10**6 + 1) / (10**6 * math.log(10**6)) - 1) < 0.08  # ln n! ~ n ln n (convergence lente)
# ar-l2-series-sommation-comparaison-ex3 (+ v1, v2) : limites des différences.
assert limit(1 / sin(x) ** 2 - 1 / x**2, x, 0) == Rational(1, 3)
assert limit(1 / log(1 + x) - 1 / x, x, 0) == Rational(1, 2)
assert simplify((1 + x**2) ** 2 / x**2 - 1 / x**2 - (2 + x**2)) == 0
for f, cible in [(math.sin, lambda m: math.sqrt(3 / m)), (math.log1p, lambda m: 2 / m),
                 (lambda u: u / (1 + u * u), lambda m: 1 / math.sqrt(2 * m))]:
    u = 1.0
    for m in range(1, 200001):
        u = f(u)
    assert abs(u / cible(200000) - 1) < 2e-2
# ar-l2-series-sommation-comparaison-ex4 : a_{n+1} - a_n ~ -1/(2n^2) et H_n - ln n - gamma ~ 1/(2n).
assert eq_lim(1 / (n + 1) - log(1 + 1 / n), -1 / (2 * n**2))
m0 = 1000
Hm = mp.harmonic(m0)
assert abs((Hm - mp.log(m0) - mp.euler) * 2 * m0 - 1) < 1e-3

# ---------------------------------------------------------------------------
# ar-l2-series-integrale (carte)
# ---------------------------------------------------------------------------
# Contre-exemple 1 : integrale de sin^2(pi t) sur [n, n+1] = 1/2.
assert integrate(sin(pi * t) ** 2, (t, 3, 4)) == Rational(1, 2)
# Exemple : gamma ~ 0,5772.
assert abs(mp.euler - mp.mpf("0.5772")) < 1e-4
# Trou 1, leurre « croissante » (f(t) = t) : int_1^{n+1} t dt > sum_{k<=n} k.
assert all(F((m + 1) ** 2 - 1, 2) > F(m * (m + 1), 2) for m in range(1, 50))
# Trou 2, leurres (f = 1/t^2) : n = 1 et n = 2.
assert 1 > 0 and 1 > F(1, 2) and F(5, 4) > F(1, 4) + F(1, 2)
# Affirmation 1 : primitive ln ln t ; affirmation 2 : reste de 1/k^2 <= 1/n.
assert all(float(mp.nsum(lambda j: 1 / j**2, [m + 1, mp.inf])) <= 1 / m for m in (1, 2, 5, 50))

# ar-l2-series-integrale-ex1 (+ v1, v2) : limites ln 2, ln 3, ln(5/2).
for a, b, rep, leurres in [(1, 2, math.log(2), [0.5, 1]), (1, 3, math.log(3), [2 / 3, 2, math.log(2)]),
                           (2, 5, math.log(2.5), [math.log(3), 0.6, math.log(5)])]:
    m0 = 10**5
    s = math.fsum(1 / j for j in range(a * m0 + 1, b * m0 + 1))
    assert abs(s - rep) < 1e-4 and all(abs(l - rep) > 0.05 for l in leurres)
    assert math.log((b * m0 + 1) / (a * m0 + 1)) <= s <= math.log(b / a)
# ar-l2-series-integrale-ex2 : primitives de 1/(t (ln t)^beta).
assert simplify(integrate(1 / (t * log(t) ** 3), t) + 1 / (2 * log(t) ** 2)) == 0
# ar-l2-series-integrale-ex3 (+ v1, v2) : primitives et équivalents.
assert simplify(integrate(t ** Rational(-1, 3), t) - Rational(3, 2) * t ** Rational(2, 3)) == 0
assert simplify(integrate(log(t) / t, t) - log(t) ** 2 / 2) == 0
assert simplify((1 - log(t)) / t**2 - (log(t) / t).diff(t)) == 0
m0 = 10**6
S_b = math.fsum(1 / (j * math.log(j)) for j in range(2, m0 + 1))  # encadrement de la correction (convergence lente)
L2 = math.log(math.log(2))
assert math.log(math.log(m0 + 1)) - L2 <= S_b <= 1 / (2 * math.log(2)) + math.log(math.log(m0)) - L2
assert abs(math.fsum(j ** (-1 / 3) for j in range(1, m0 + 1)) / (1.5 * m0 ** (2 / 3)) - 1) < 1e-3
assert abs(math.fsum(math.log(j) / j for j in range(1, m0 + 1)) / (math.log(m0) ** 2 / 2) - 1) < 0.02
# ar-l2-series-integrale-ex4 (+ v1, v2) : primitives et encadrements des restes.
assert simplify(integrate(log(t) / t**2, t) + (log(t) + 1) / t) == 0
assert simplify((1 - 2 * log(t)) / t**3 - (log(t) / t**2).diff(t)) == 0 and 2 > math.sqrt(math.e)
# Encadrement sur une tranche finie [m0+1, M] : int_{m0+1}^{M+1} f <= somme <= int_{m0}^{M} f (F = primitive tendant vers 0).
for f, Fp in [(lambda u: 1 / (u * math.log(u) ** 2), lambda u: -1 / math.log(u)),
              (lambda u: math.log(u) / u**2, lambda u: -(math.log(u) + 1) / u),
              (lambda u: 1 / (u * math.log(u) ** 3), lambda u: -1 / (2 * math.log(u) ** 2))]:
    m0, M = 50, 10**5
    tranche = math.fsum(f(j) for j in range(m0 + 1, M + 1))
    assert Fp(M + 1) - Fp(m0 + 1) <= tranche <= Fp(M) - Fp(m0)

# ---------------------------------------------------------------------------
# ar-l2-series-cauchy-absolue (carte)
# ---------------------------------------------------------------------------
# Contre-exemple 2 : sum_{k=p+1}^{2p} 1/k >= 1/2.
assert all(sum(F(1, j) for j in range(p + 1, 2 * p + 1)) >= F(1, 2) for p in range(1, 60))
# Trou 3, leurres : u = (1, -1, 0, ...) (égalité et inégalité inverse fausses) ; u_n = 1/2^n (sup = 1 < 2).
assert abs(1 - 1) < abs(1) + abs(-1)
# Affirmation 4 : sum_{k=n+1}^{2n} 1/(k ln k) <= 1/ln n.
assert all(sum(1 / (j * math.log(j)) for j in range(m + 1, 2 * m + 1)) <= 1 / math.log(m) for m in range(2, 200))

# ar-l2-series-cauchy-absolue-ex2 (+ v1, v2) : |a + i|/2 < 1 <=> a < sqrt 3, etc.
for a2max, rayon, im in [(3, 2, 1), (5, 3, 2), (1, math.sqrt(2), 1)]:
    a = math.sqrt(a2max)
    assert abs(abs(complex(a, im)) / rayon - 1) < 1e-12  # frontière : module 1
    assert abs(complex(a * 0.999, im)) / rayon < 1 and abs(complex(a * 1.001, im)) / rayon > 1
# ar-l2-series-cauchy-absolue-ex3 : 1/(n + sqrt n) décroissante, n u_n -> 1.
assert limit(n / (n + sqrt(n)), n, oo) == 1

# ---------------------------------------------------------------------------
# ar-l2-series-alternees (carte)
# ---------------------------------------------------------------------------
# Cours et contre-exemple 1 : (-1)^n/(sqrt n + (-1)^n) = (-1)^n/sqrt n - 1/n + O(n^{-3/2}) (n pair et impair).
for sgn in (1, -1):
    d = sgn / (sqrt(x) + sgn) - (sgn / sqrt(x) - 1 / x)
    assert limit(d * x ** Rational(3, 2), x, oo) == sgn  # reste en n^{-3/2}
# Exemple : 1/(2n+3) <= 1/100 <=> n >= 49 (48,5).
assert plus_petit(lambda m: F(1, 2 * m + 3), F(1, 100)) == 49
# Trous : a = 1/(n+1) donne R_0 = ln 2 - 1 < 0 ; |R_0| = 1 - ln 2 < 1/2 (leurre |R_n| >= a_{n+1}).
assert math.log(2) - 1 < 0 and 1 - math.log(2) < 0.5
# Leurres |R_n| <= a_{n+2} et <= a_{n+1}/2 : a = (1, 1, 0, 0, ...), R_0 = -1.

# ar-l2-series-alternees-ex1 (+ v1, v2) : plus petits n.
assert plus_petit(lambda m: F(1, m + 2), F(1, 1000)) == 998
assert plus_petit(lambda m: F(1, (m + 2) ** 2), F(1, 10**4)) == 98
assert plus_petit(lambda m: F(1, math.factorial(m + 1)), F(1, 10**6)) == 9
assert math.factorial(9) == 362880 and math.factorial(10) == 3628800
# ar-l2-series-alternees-ex2 (+ v1, v2) : sommes partielles S_2, S_3 et leurres faux.
def S(a, m):
    return sum((-1) ** j * a(j) for j in range(m + 1))


for a, s2, s3, val, leurres in [
    (lambda j: F(1, 2 * j + 1), F(13, 15), F(76, 105), math.pi / 4, [(F(13, 15), 1), (F(2, 3), F(76, 105)), (0, F(2, 3))]),
    (lambda j: F(1, j + 1), F(5, 6), F(7, 12), math.log(2), [(F(5, 6), 1), (F(1, 2), F(7, 12)), (0, F(1, 2))]),
    (lambda j: F(1, (j + 1) ** 2), F(31, 36), F(115, 144), math.pi**2 / 12,
     [(F(31, 36), 1), (F(3, 4), F(115, 144)), (0, F(3, 4))]),
]:
    assert S(a, 2) == s2 and S(a, 3) == s3 and s3 <= val <= s2
    assert all(not (lo <= val <= hi) for lo, hi in leurres)
# ar-l2-series-alternees-ex3 (+ v1, v2) : développements.
for sgn in (1, -1):
    d1 = sgn / (x + sgn) - (sgn / x - 1 / x**2)
    assert limit(d1 * x**3, x, oo) == sgn
    d2 = log(1 + sgn / sqrt(x)) - (sgn / sqrt(x) - 1 / (2 * x))
    assert limit(d2 * x ** Rational(3, 2), x, oo) == Rational(sgn, 3)
# ar-l2-series-alternees-ex4 : r_n + r_{n+1} = 1/(n+1), R_n = (-1)^{n+1} r_n, r_n ~ 1/(2n).
for m in (0, 3, 10):
    r = lambda j: mp.quad(lambda u: u**j / (1 + u), [0, 1])
    assert abs(r(m) + r(m + 1) - mp.mpf(1) / (m + 1)) < 1e-20
    Rm = mp.nsum(lambda j: (-1) ** j / j, [m + 1, mp.inf])
    assert abs(Rm - (-1) ** (m + 1) * r(m)) < 1e-15
assert abs(mp.quad(lambda u: u**2000 / (1 + u), [0, 1]) * 2 * 2000 - 1) < 1e-3

# ---------------------------------------------------------------------------
# ar-l2-series-abel (carte)
# ---------------------------------------------------------------------------
# Formule d'Abel et leurres du trou 3, sur des suites aléatoires fixes.
import random  # noqa: E402

random.seed(1)
for _ in range(20):
    N = random.randint(1, 8)
    a = [random.randint(-5, 5) for _ in range(N + 1)]
    b = [random.randint(-5, 5) for _ in range(N + 1)]
    B = [sum(b[: j + 1]) for j in range(N + 1)]
    gauche = sum(a[j] * b[j] for j in range(N + 1))
    assert gauche == a[N] * B[N] + sum((a[j] - a[j + 1]) * B[j] for j in range(N))
assert 1 * 1 + 0 * 0 != 0 * 1 - (1 - 0) * 1  # leurre « - » : a = (1, 0), b = (1, 0)
assert 1 * 1 + 1 * 0 != 1 * 1 + 1 * 1  # leurre sum a_k B_k : a = (1, 1), b = (1, 0)
assert 1 * 1 + 0 * 0 != 1 * 1 + (0 - 1) * 1  # leurre a_0 B_n + ... : a = (1, 0), b = (1, 0)
# Majoration |sum e^{ik theta}| <= 1/|sin(theta/2)| et |1 - e^{i theta}| = 2|sin(theta/2)|.
for th in [0.3, 1, 2, 3, 5]:
    assert abs(abs(1 - cmath.exp(1j * th)) - 2 * abs(math.sin(th / 2))) < 1e-12
    s = 0
    for m in range(500):
        s += cmath.exp(1j * m * th)
        assert abs(s) <= 1 / abs(math.sin(th / 2)) + 1e-9
# Exemple : sum sin n / n = (pi - 1)/2.
# (somme de -log(1 - e^{i}) = sum e^{in}/n, partie imaginaire)
assert abs(float(mp.im(-mp.log(1 - mp.exp(1j)))) - (math.pi - 1) / 2) < 1e-15

# ar-l2-series-abel-ex1 (+ v1, v2) : plus petites constantes 2, sqrt 2, 1.
for th, rep, leurres in [(math.pi / 3, 2, [1, math.sqrt(3), 2 / math.sqrt(3)]),
                         (math.pi / 2, math.sqrt(2), [1, 2, math.sqrt(2) / 2]),
                         (2 * math.pi / 3, 1, [2 / math.sqrt(3), math.sqrt(3), 2])]:
    s, mx = 0, 0
    for m in range(60):
        s += cmath.exp(1j * m * th)
        mx = max(mx, abs(s))
    assert abs(mx - rep) < 1e-9 and all(abs(l - rep) > 0.1 for l in leurres)
# ar-l2-series-abel-ex3 (+ v1, v2) : plus petits n et leurres (sans facteur sin, avec facteur 2, n + 1).
for M, eps, p, rep, leurres in [(1 / math.sin(0.5), 1e-2, 1, 208, [209, 99, 417]),
                                (1 / math.sin(1), 1e-2, 1, 118, [119, 99, 237]),
                                (1 / math.sin(0.5), 1e-1, 0.5, 435, [436, 99, 1740])]:
    borne = lambda m: M / (m + 1) ** p
    assert plus_petit(borne, eps) == rep
    assert math.ceil((1 / eps) ** (1 / p)) - 1 == 99 and math.ceil((2 * M / eps) ** (1 / p)) - 1 == leurres[2]
# Vérification numérique de la majoration du reste de sum sin k / k pour quelques n.
for m in (5, 20, 100):
    Rm = (math.pi - 1) / 2 - math.fsum(math.sin(j) / j for j in range(1, m + 1))
    assert abs(Rm) <= 1 / ((m + 1) * math.sin(0.5))
# ar-l2-series-abel-ex2 : |sin n| >= sin^2 n.
assert all(abs(math.sin(m)) >= math.sin(m) ** 2 for m in range(1, 1000))

# ---------------------------------------------------------------------------
# ar-l2-series-produit-cauchy (carte)
# ---------------------------------------------------------------------------
def cauchy(a, b, N):
    return [sum(a(j) * b(m - j) for j in range(m + 1)) for m in range(N)]


# Exemple : sum (n+1)/2^n = 4.
assert abs(sum((m + 1) / 2**m for m in range(200)) - 4) < 1e-12
# Contre-exemple 1 : |c_n| >= 2(n+1)/(n+2) >= 1 pour a_n = b_n = (-1)^n/sqrt(n+1).
c = cauchy(lambda j: (-1) ** j / math.sqrt(j + 1), lambda j: (-1) ** j / math.sqrt(j + 1), 200)
assert all(abs(c[m]) >= 2 * (m + 1) / (m + 2) - 1e-12 >= 1 - 1e-12 for m in range(200))
# Contre-exemple 2 : sum 1/4^n = 4/3.
# Trou 1, leurre a_k b_{n+k} avec a = b = 1/2^n : somme 16/7 != 4.
s_leurre = sum(sum(F(1, 2**j) * F(1, 2 ** (m + j)) for j in range(m + 1)) for m in range(60))
assert abs(float(s_leurre) - 16 / 7) < 1e-12
# Trou 3, leurre sur les modules : a = (1, 1), b = (1, -1) : sum |c| = 2 != 4.
cc = cauchy(lambda j: [1, 1][j] if j < 2 else 0, lambda j: [1, -1][j] if j < 2 else 0, 3)
assert cc == [1, 0, -1] and sum(abs(v) for v in cc) == 2
# Affirmation 3 : a = (2, 2, 4, 8, ...), b = (-1, 1, 1, ...) : c_0 = -2, c_n = 0.
cc = cauchy(lambda j: 2 if j == 0 else 2**j, lambda j: -1 if j == 0 else 1, 40)
assert cc[0] == -2 and all(v == 0 for v in cc[1:])
# Affirmation 2 : binôme.
assert simplify(sum(binomial(5, j) * x**j * t ** (5 - j) for j in range(6)) - (x + t) ** 5) == 0

# ar-l2-series-produit-cauchy-ex1 (+ v1, v2) : 1/(1-q)^2 et leurres.
for qq, rep, leurres in [(F(1, 3), F(9, 4), [F(3, 2), F(3, 4), F(9, 2)]),
                         (F(-1, 2), F(4, 9), [F(2, 3), F(4), F(-4, 9)]),
                         (F(1, 4), F(16, 9), [F(4, 3), F(16, 3), F(4, 9)])]:
    assert 1 / (1 - qq) ** 2 == rep and rep not in leurres
    assert abs(sum((m + 1) * float(qq) ** m for m in range(300)) - float(rep)) < 1e-12
# ar-l2-series-produit-cauchy-ex2 (+ v1, v2) : expressions de c_n et sommes.
for r1, r2, expr, total, leurres in [
    (F(1, 2), F(1, 3), lambda m: 3 * F(1, 2) ** m - 2 * F(1, 3) ** m, 3,
     [lambda m: F(1, 6) ** m, lambda m: 2 * F(1, 2) ** m - F(1, 3) ** m, lambda m: (m + 1) * F(1, 6) ** m]),
    (F(1, 2), F(1, 4), lambda m: 2 * F(1, 2) ** m - F(1, 4) ** m, F(8, 3),
     [lambda m: F(1, 8) ** m, lambda m: 4 * F(1, 2) ** m - 3 * F(1, 4) ** m, lambda m: (m + 1) * F(1, 8) ** m]),
    (F(2, 3), F(1, 3), lambda m: 2 * F(2, 3) ** m - F(1, 3) ** m, F(9, 2),
     [lambda m: F(2, 9) ** m, lambda m: 3 * F(2, 3) ** m - 2 * F(1, 3) ** m, lambda m: (m + 1) * F(2, 9) ** m]),
]:
    cc = cauchy(lambda j: r1**j, lambda j: r2**j, 30)
    assert all(cc[m] == expr(m) for m in range(30))
    assert 1 / (1 - r1) * 1 / (1 - r2) == total
    assert all(any(cc[m] != l(m) for m in range(3)) for l in leurres)
# ar-l2-series-produit-cauchy-ex3 (+ v1, v2) : 2e, 3/(2e), 3e^2/2.
for a, b, rep in [(lambda j: 1 / math.factorial(j), lambda j: 1 / 2**j, 2 * math.e),
                  (lambda j: (-1) ** j / math.factorial(j), lambda j: 1 / 3**j, 3 / (2 * math.e)),
                  (lambda j: 2**j / math.factorial(j), lambda j: 1 / 3**j, 3 * math.e**2 / 2)]:
    assert abs(sum(cauchy(a, b, 120)) - rep) < 1e-10
assert all(abs(l - 2 * math.e) > 0.1 for l in [math.e, math.e / 2, math.e + 2])
assert all(abs(l - 3 / (2 * math.e)) > 0.1 for l in [3 * math.e / 2, 2 / (3 * math.e), 1 / math.e + 1.5])
assert all(abs(l - 3 * math.e**2 / 2) > 0.1 for l in [3 * math.e**2, 2 * math.e**2 / 3, math.e**2 + 1.5])

# ---------------------------------------------------------------------------
# QCM ar-l2-series-q1 … q8
# ---------------------------------------------------------------------------
# q1 : 1 - cos(1/n) ~ 1/(2n^2) ; les trois autres ~ 1/n.
assert eq_lim(1 - cos(1 / n), 1 / (2 * n**2))
assert all(eq_lim(e, 1 / n) for e in [sin(1 / n), log(1 + 1 / n), exp(1 / n) - 1])
# q3 : ln(1 + 1/sqrt k) ~ 1/sqrt k.
assert eq_lim(log(1 + 1 / sqrt(n)), 1 / sqrt(n))
# q4 : ln(n+1) <= H_n <= 1 + ln n.
Hs = 0
for m in range(1, 3000):
    Hs += 1 / m
    assert math.log(m + 1) <= Hs <= 1 + math.log(m) + 1e-12
# q6 : S - 7/8 = R_1 ~ 0,0265, entre 37/1728 et 1/27.
R1 = 3 / 4 * float(zeta(3)) - 7 / 8
assert F(37, 1728) <= R1 <= F(1, 27) and R1 > 1 / 64
# q8 : (1 - 1)^n / n!.
cc = cauchy(lambda j: F(1, math.factorial(j)), lambda j: F((-1) ** j, math.factorial(j)), 20)
assert cc[0] == 1 and all(v == 0 for v in cc[1:])

print("ar_l2_series : toutes les vérifications passent.")
