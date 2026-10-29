"""Vérifications calculatoires de l'enrichissement du chapitre « ar-l1-suites » (analyse réelle, L1) :
questions de content/qcm/ar-l1-suites.json et exercices ajoutés à content/exercices/analyse-reelle/l1-ar-l1-suites.json.

Chaque assert cite l'identifiant de la question ou de l'exercice (variante comprise : v1, v2).
Lancer : python3 verification/enrichi_ar_l1_suites.py (doit finir sans erreur).
"""

import math
from fractions import Fraction as F

from sympy import E, Rational, atan, cos, exp, limit, log, oo, series, simplify, sin, sqrt, symbols

n = symbols("n", positive=True, integer=True)
x = symbols("x", positive=True)


def suite2(u0, u1, a, b, c, N):
    """u_0..u_N pour u_{n+2} = a u_{n+1} + b u_n + c (calcul exact)."""
    u = [u0, u1]
    while len(u) <= N:
        u.append(a * u[-1] + b * u[-2] + c)
    return u


# ---------------------------------------------------------------------------
# ar-l1-suites-q1 : u0 = 0, u1 = 1, u_{n+2} = u_{n+1} + 2 u_n : u_n = (2^n - (-1)^n)/3.
# ---------------------------------------------------------------------------
u = suite2(0, 1, 1, 2, 0, 30)
assert all(u[k] == F(2**k - (-1) ** k, 3) for k in range(31))
leurres = [lambda k: F(2**k + (-1) ** k, 3), lambda k: F(2**k - (-1) ** k, 2), lambda k: F(2**k - 1)]
assert all(any(f(k) != u[k] for k in range(5)) for f in leurres)

# ---------------------------------------------------------------------------
# ar-l1-suites-q2 : (1 + 1/k)^k -> e, donc par Cesàro la moyenne tend vers e (leurres 1, 0, +oo faux).
# ---------------------------------------------------------------------------
assert limit((1 + 1 / n) ** n, n, oo) == E
moy = sum((1 + 1 / k) ** k for k in range(1, 200001)) / 200000
assert abs(moy - math.e) < 1e-3

# ---------------------------------------------------------------------------
# ar-l1-suites-q3 (vrai) : u_n = 1 - 1/n croissante, v_n = 1 + 1/n^2 décroissante, v_n - u_n -> 0.
# ---------------------------------------------------------------------------
assert simplify((1 - 1 / (n + 1)) - (1 - 1 / n) - 1 / (n * (n + 1))) == 0
assert all(F(1, (k + 1) ** 2) < F(1, k**2) for k in range(1, 100))
assert limit(1 / n**2 + 1 / n, n, oo) == 0

# ---------------------------------------------------------------------------
# ar-l1-suites-q5 (faux) : u_n = 1 si n est premier avec 6, 0 sinon ; u_{2n} = u_{3n} = 0 mais u_{6k+1} = 1.
# ---------------------------------------------------------------------------
w = [1 if math.gcd(k, 6) == 1 else 0 for k in range(0, 3000)]
assert all(w[2 * k] == 0 for k in range(1500)) and all(w[3 * k] == 0 for k in range(1000))
assert all(w[6 * k + 1] == 1 for k in range(400))

# ---------------------------------------------------------------------------
# ar-l1-suites-q6 (faux) : u_n = n (n pair), 0 (n impair) : non bornée, sous-suite des rangs impairs constante.
# ---------------------------------------------------------------------------
v = [k if k % 2 == 0 else 0 for k in range(100)]
assert max(v) == 98 and all(v[2 * k + 1] == 0 for k in range(50))

# ---------------------------------------------------------------------------
# ar-l1-suites-q7 : u_n = (n + (-1)^n n^2)/(n^2 + 1) : rangs pairs -> 1, impairs -> -1, suite bornée.
# ---------------------------------------------------------------------------
assert limit((2 * n + 4 * n**2) / (4 * n**2 + 1), n, oo) == 1
assert limit((2 * n + 1 - (2 * n + 1) ** 2) / ((2 * n + 1) ** 2 + 1), n, oo) == -1
assert all(abs(F(k + (-1) ** k * k * k, k * k + 1)) <= 2 for k in range(0, 500))  # leurre « non bornée » faux

# ---------------------------------------------------------------------------
# ar-l1-suites-q8 : S_4 = 7/12, S_5 = 47/60 encadrent ln 2 ; les leurres ne contiennent pas ln 2.
# ---------------------------------------------------------------------------
S = [sum(F((-1) ** (k + 1), k) for k in range(1, m + 1)) for m in range(0, 6)]
assert S[4] == F(7, 12) and S[5] == F(47, 60) and S[2] == F(1, 2) and S[3] == F(5, 6) and S[1] == 1
l2 = math.log(2)
assert 7 / 12 < l2 < 47 / 60
assert not (1 / 2 <= l2 <= 7 / 12) and not (47 / 60 <= l2 <= 5 / 6) and not (5 / 6 <= l2 <= 1)

# ---------------------------------------------------------------------------
# ar-l1-suites-cesaro-ex5 : u_{n+1} = sin u_n.
# ---------------------------------------------------------------------------
# 1/sin^2 x - 1/x^2 -> 1/3 ; développement admis de sin à l'ordre 3
assert limit(1 / sin(x) ** 2 - 1 / x**2, x, 0) == Rational(1, 3)
assert series(sin(x), x, 0, 4).removeO() == x - x**3 / 6
assert limit((x**2 - (x - x**3 / 6) ** 2) / x**4, x, 0) == Rational(1, 3)
t = 1.0
for k in range(1, 2_000_001):
    t = math.sin(t)
assert abs(math.sqrt(2_000_000) * t - math.sqrt(3)) < 5e-3  # sqrt(n) u_n -> sqrt(3) ; leurres 3, 1/sqrt3, 0
assert all(abs(val - math.sqrt(3)) > 0.5 for val in (3, 1 / math.sqrt(3), 0))
# v1 : u_{n+1} = ln(1 + u_n) : 1/ln(1+x) - 1/x -> 1/2, n u_n -> 2 ; leurres 1/2, 1, 0
assert limit(1 / log(1 + x) - 1 / x, x, 0) == Rational(1, 2)
assert limit((x - (x - x**2 / 2)) / x**2, x, 0) == Rational(1, 2)
t = 1.0
for k in range(1, 2_000_001):
    t = math.log1p(t)
assert abs(2_000_000 * t - 2) < 2e-2
# v2 : u_{n+1} = u_n/(1 + u_n^2) : 1/u_{n+1}^2 - 1/u_n^2 = 2 + u_n^2, sqrt(n) u_n -> 1/sqrt2 ; leurres sqrt2, 1/2, 2
assert simplify(1 / (x / (1 + x**2)) ** 2 - 1 / x**2 - (2 + x**2)) == 0
t = 1.0
for k in range(1, 2_000_001):
    t = t / (1 + t * t)
assert abs(math.sqrt(2_000_000) * t - math.sqrt(2) / 2) < 5e-3
assert all(abs(val - math.sqrt(2) / 2) > 0.2 for val in (math.sqrt(2), 0.5, 2))

# ---------------------------------------------------------------------------
# ar-l1-suites-recurrence-lineaire-ex5 et variantes : second membre constant.
# ---------------------------------------------------------------------------
# exercice : u_{n+2} = 5u_{n+1} - 6u_n + 4, u0 = 3, u1 = 5 -> u_n = 3^n + 2 ; c = 2
assert 2 == 5 * 2 - 6 * 2 + 4
u = suite2(3, 5, 5, -6, 4, 25)
assert all(u[k] == 3**k + 2 for k in range(26))
for f in (lambda k: 2**k + 3**k, lambda k: 2 * 3**k + 1, lambda k: 3**k + 2 * k + 2):
    assert any(f(k) != u[k] for k in range(4))
# v1 : u_{n+2} = u_{n+1} + 2u_n - 4, u0 = 2, u1 = 5 -> u_n = 2^n - (-1)^n + 2 ; c = 2
assert 2 == 2 + 2 * 2 - 4
u = suite2(2, 5, 1, 2, -4, 25)
assert all(u[k] == 2**k - (-1) ** k + 2 for k in range(26))
for f in (lambda k: 2**k + (-1) ** k, lambda k: 2**k - (-1) ** k + 2 * k + 2, lambda k: 2 ** (k + 1) - (-1) ** k + 1):
    assert any(f(k) != u[k] for k in range(4))
# v2 : u_{n+2} = 6u_{n+1} - 9u_n + 8, u0 = 3, u1 = 8 -> u_n = (n+1)3^n + 2 ; c = 2
assert 2 == 6 * 2 - 9 * 2 + 8
u = suite2(3, 8, 6, -9, 8, 25)
assert all(u[k] == (k + 1) * 3**k + 2 for k in range(26))
for f in (lambda k: 3**k + 2, lambda k: (2 * k + 1) * 3**k + 2, lambda k: (k + 1) * 3**k + 2 * k + 2):
    assert any(f(k) != u[k] for k in range(4))

# ---------------------------------------------------------------------------
# ar-l1-bolzano-weierstrass-ex4 : lemme des pics (essai numérique sur une suite finie).
# ---------------------------------------------------------------------------
import random

random.seed(2)
for _ in range(50):
    L = [random.randint(-20, 20) for _ in range(60)]
    pics = [i for i in range(len(L)) if all(L[i] >= L[j] for j in range(i + 1, len(L)))]
    assert all(L[pics[k]] >= L[pics[k + 1]] for k in range(len(pics) - 1))  # les pics forment une suite décroissante

# ===========================================================================
# Relecture indépendante (2026-10-29) : réponses retrouvées à l'aveugle, recalculées autrement.
# ===========================================================================


def _iter(u0, u1, a, b, c, N=30):
    u = [u0, u1]
    while len(u) < N:
        u.append(a * u[-1] + b * u[-2] + c)
    return u


# ar-l1-suites-recurrence-lineaire-ex5 : u_n = 3^n + 2 ; v1 : 2^n - (-1)^n + 2 ; v2 : (n+1)3^n + 2 (30 termes).
assert _iter(3, 5, 5, -6, 4) == [3**k + 2 for k in range(30)]
assert _iter(2, 5, 1, 2, -4) == [2**k - (-1) ** k + 2 for k in range(30)]
assert _iter(3, 8, 6, -9, 8) == [(k + 1) * 3**k + 2 for k in range(30)]
# leurres faux (premier terme où ils diffèrent)
assert _iter(3, 5, 5, -6, 4)[:2] != [2**k + 3**k for k in range(2)]
assert [2 * 3**k + 1 for k in range(2)] != [3, 5] and [3**k + 2 * k + 2 for k in range(2)] != [3, 5]
assert [2**k + (-1) ** k for k in range(2)] != [2, 5] and [2**k - (-1) ** k + 2 * k + 2 for k in range(2)] != [2, 5]
assert [2 ** (k + 1) - (-1) ** k + 1 for k in range(2)] != [2, 5]
assert [3**k + 2 for k in range(2)] != [3, 8] and [(2 * k + 1) * 3**k + 2 for k in range(2)] != [3, 8]
assert [(k + 1) * 3**k + 2 * k + 2 for k in range(2)] != [3, 8]

# ar-l1-suites-q1 : u_n = (2^n - (-1)^n)/3.
assert _iter(0, 1, 1, 2, 0) == [(2**k - (-1) ** k) // 3 for k in range(30)]

# ar-l1-suites-cesaro-ex5 : limites des quotients (SymPy) et contrôle numérique de sqrt(n) u_n.
_x = symbols("_x", positive=True)
assert limit(1 / sin(_x) ** 2 - 1 / _x**2, _x, 0) == Rational(1, 3)
assert limit(1 / log(1 + _x) - 1 / _x, _x, 0, "+") == Rational(1, 2)
assert simplify((1 + _x**2) ** 2 / _x**2 - 1 / _x**2 - (2 + _x**2)) == 0  # v2
_u, _v, _w = 1.0, 1.0, 1.0
for _k in range(10**6):
    _u, _v, _w = math.sin(_u), math.log1p(_v), _w / (1 + _w * _w)
_N = 10**6
assert abs(math.sqrt(_N) * _u - math.sqrt(3)) < 1e-2  # principal : sqrt(3)
assert abs(_N * _v - 2) < 1e-2  # v1 : 2
assert abs(math.sqrt(_N) * _w - math.sqrt(2) / 2) < 1e-2  # v2 : sqrt(2)/2

# ar-l1-suites-q5 : u_n = 1 si pgcd(n, 6) = 1 ; u_{2n} = u_{3n} = 0 mais u_{6k+1} = 1.
_c = lambda k: 1 if math.gcd(k, 6) == 1 else 0
assert all(_c(2 * k) == 0 and _c(3 * k) == 0 and _c(6 * k + 1) == 1 for k in range(500))

# ar-l1-suites-q8 : S_4 = 7/12, S_5 = 47/60, et ln 2 entre les deux.
_S = [sum(F((-1) ** (k + 1), k) for k in range(1, m + 1)) for m in (4, 5)]
assert _S == [F(7, 12), F(47, 60)] and F(7, 12) < math.log(2) < F(47, 60)

print("enrichi_ar_l1_suites : toutes les vérifications passent.")
