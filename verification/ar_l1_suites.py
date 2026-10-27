"""Vérifications calculatoires du chapitre « ar-l1-suites » (analyse réelle, L1) :
cartes et exercices de content/*/analyse-reelle/l1-ar-l1-suites.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise : v1, v2).
Lancer : python3 verification/ar_l1_suites.py (doit finir sans erreur).
"""

import math
from fractions import Fraction as F

from sympy import (E, I, Rational, cos, exp, expand, factor, factorial, limit, log, oo, pi, simplify, sin, sqrt,
                   symbols, zeta)

n = symbols("n", positive=True, integer=True)
X = symbols("x", real=True)


def suite2(u0, u1, a, b, N):
    """u_0..u_N pour u_{n+2} = a u_{n+1} + b u_n (calcul exact)."""
    u = [u0, u1]
    while len(u) <= N:
        u.append(a * u[-1] + b * u[-2])
    return u


# ---------------------------------------------------------------------------
# ar-l1-suites-adjacentes (carte)
# ---------------------------------------------------------------------------
# Exemple : approximations décimales de sqrt(2) à 10^-3.
assert F(1414, 1000) < math.sqrt(2) < F(1415, 1000)
# Contre-exemple 1 : u_n = -1/(n+1) croissante, v_n = 1 + 1/(n+1) décroissante, limites 0 et 1.
assert all(-F(1, k + 1) < -F(1, k + 2) and 1 + F(1, k + 1) > 1 + F(1, k + 2) for k in range(100))
# Contre-exemple 2 : v_n = n + 1/(n+1) croissante.
assert all((k + 1 + F(1, k + 2)) - (k + F(1, k + 1)) > 0 for k in range(100))
# Trou 1, leurre « majorée » : v_n = -n, u_n = -n - 1/(n+1) (u majorée par 0, v - u -> 0, divergence).
assert all(-k - F(1, k + 1) <= 0 for k in range(100))
# Trou 1, leurre « positive » : u_n = 2 + (1+(-1)^n)/(n+1) > 0, v_n = 2 + 2/(n+1) décroissante, v - u -> 0, limite 2, u_0 = 4 > 2.
uu = [2 + F(1 + (-1) ** k, k + 1) for k in range(200)]
vv = [2 + F(2, k + 1) for k in range(200)]
assert all(t > 0 for t in uu) and all(vv[k + 1] < vv[k] for k in range(199)) and uu[0] == 4 > 2
assert all(vv[k] - uu[k] in (0, F(2, k + 1)) for k in range(200))
# Trou 3, leurre |v_n - l| <= |u_n - l| : u_n = -1/(n+1), v_n = 2/(n+1), l = 0 ; n = 0 : 2 > 1.
assert abs(F(2, 1)) > abs(F(-1, 1))
# Affirmation 4 : (-1)^n/(n+1) et 2(-1)^n/(n+1) ne sont pas monotones.
w = [F((-1) ** k, k + 1) for k in range(4)]
assert w[1] < w[0] and w[2] > w[1]

# ar-l1-suites-adjacentes-ex1 (+ v1, v2) : encadrements u_n <= l <= v_n.
k = symbols("k", positive=True, integer=True)
assert simplify(1 / (n + 1) ** 2 + 1 / (n + 1) - 1 / n + 1 / (n * (n + 1) ** 2)) == 0
assert simplify(1 / (n + 1) ** 3 + 1 / (n + 1) ** 2 - 1 / n**2 + (n**2 + 3 * n + 1) / (n**2 * (n + 1) ** 3)) == 0
for p, corr, N, l, (bu, bv), leurres in [
    (2, lambda m: F(1, m), 4, math.pi**2 / 6, (F(205, 144), F(241, 144)),
     [(F(205, 144), F(107, 72)), (F(169, 144), F(205, 144)), (F(241, 144), F(277, 144))]),
    (3, lambda m: F(1, m * m), 3, float(zeta(3)), (F(251, 216), F(275, 216)),
     [(F(251, 216), F(259, 216)), (F(227, 216), F(251, 216)), (F(275, 216), F(299, 216))]),
    (2, lambda m: F(1, m), 5, math.pi**2 / 6, (F(5269, 3600), F(5989, 3600)),
     [(F(5269, 3600), F(5413, 3600)), (F(4549, 3600), F(5269, 3600)), (F(5989, 3600), F(6709, 3600))]),
]:
    u = lambda m: sum(F(1, j**p) for j in range(1, m + 1))
    assert u(N) == bu and u(N) + corr(N) == bv and bu <= l <= bv
    assert all(u(m + 1) > u(m) and u(m + 1) + corr(m + 1) < u(m) + corr(m) for m in range(1, 60))
    for a, b in leurres:  # chaque leurre est un encadrement FAUX de l
        assert not (a <= l <= b)

# ar-l1-suites-adjacentes-ex2 (+ v1, v2) : u_{n+1} = (u+2v)/3, v_{n+1} = (u+3v)/4 ; 3u + 8v invariant.
for u0, v0, l, leurres in [(1, 12, 9, [F(13, 2), F(25, 3), F(37, 4)]),
                           (2, 13, 10, [F(15, 2), F(28, 3), F(41, 4)]),
                           (-3, 8, 5, [F(5, 2), F(13, 3), F(21, 4)])]:
    u, v = F(u0), F(v0)
    assert F(3 * u0 + 8 * v0, 11) == l and v0 - u0 == 11
    assert leurres == [F(u0 + v0, 2), F(u0 + 2 * v0, 3), F(u0 + 3 * v0, 4)] and l not in leurres
    for m in range(40):
        un, vn = (u + 2 * v) / 3, (u + 3 * v) / 4
        assert vn - un == (v - u) / 12 and un > u and vn < v and 3 * un + 8 * vn == 3 * u0 + 8 * v0
        u, v = un, vn
    assert abs(u - l) < F(1, 10**30) and abs(v - l) < F(1, 10**30)

# ar-l1-suites-adjacentes-ex3 : b_{n+1} - a_{n+1} <= (b_n - a_n)/2 (moyenne arithmético-géométrique).
for a, b in [(1.0, 2.0), (0.1, 50.0), (3.0, 3.5)]:
    for m in range(4):
        a2, b2 = math.sqrt(a * b), (a + b) / 2
        assert a <= a2 + 1e-12 and a2 <= b2 + 1e-12 and b2 <= b + 1e-12 and b2 - a2 <= (b - a) / 2 + 1e-12
        a, b = a2, b2

# ar-l1-suites-adjacentes-ex4 : v_{n+1} - v_n = -1/(n(n+1)(n+1)!) ; u_n < e < v_n.
for m in range(1, 30):
    um = sum(F(1, math.factorial(j)) for j in range(m + 1))
    vm = um + F(1, m * math.factorial(m))
    um1 = um + F(1, math.factorial(m + 1))
    vm1 = um1 + F(1, (m + 1) * math.factorial(m + 1))
    assert vm1 - vm == -F(1, m * (m + 1) * math.factorial(m + 1))
    assert (m * (m + 1) + m - (m + 1) ** 2) == -1
    assert um < math.e < vm or m > 15  # au-delà, la précision flottante ne suffit plus
    assert (math.factorial(m) * um).denominator == 1  # q! u_q est entier

# ---------------------------------------------------------------------------
# ar-l1-suites-pairs-impairs (carte)
# ---------------------------------------------------------------------------
# Exemple : (-1)^n n/(n+1) : limites 1 et -1 des rangs pairs et impairs.
assert limit(2 * n / (2 * n + 1), n, oo) == 1 and limit(-(2 * n + 1) / (2 * n + 2), n, oo) == -1
# Contre-exemple : u_n = 0 si n pair ou multiple de 3, 1 sinon.
uc = lambda m: 0 if m % 2 == 0 or m % 3 == 0 else 1
assert all(uc(2 * m) == 0 and uc(3 * m) == 0 and uc(6 * m + 1) == 1 for m in range(200))
# Trou 1, leurres : (u_{2n}, u_{4n}) et (u_{2n+1}, u_{4n+1}) constants pour u_n = (-1)^n, qui diverge.
assert all((-1) ** (4 * m) == 1 and (-1) ** (4 * m + 1) == -1 for m in range(50))

# ar-l1-suites-pairs-impairs-ex1 (+ v1, v2) : limites des sous-suites.
cas = [
    ((3 * 2 * n + 2 * n) / (2 * n + 1), (3 * (2 * n + 1) - (2 * n + 1)) / (2 * n + 2), (4, 2)),
    ((4 * n + 1) / (2 * n + 3), (2 * (2 * n + 1) - 1) / (2 * n + 1 + 3), (2, 2)),
    ((2 * n) / (4 * n + 1), -(2 * n + 1) / (2 * (2 * n + 1) + 1), (Rational(1, 2), Rational(-1, 2))),
]
for pair, impair, (lp, li) in cas:
    assert limit(pair, n, oo) == lp and limit(impair, n, oo) == li
# Vérification directe des expressions des sous-suites à partir de u_n.
for m in range(1, 30):
    for N in (2 * m, 2 * m + 1):
        assert F(3 * N + (-1) ** N * N, N + 1) == (F(8 * m, 2 * m + 1) if N % 2 == 0 else F(2 * (2 * m + 1), 2 * m + 2))
        assert F(2 * N + (-1) ** N, N + 3) == (F(4 * m + 1, 2 * m + 3) if N % 2 == 0 else F(4 * m + 1, 2 * m + 4))
        assert F((-1) ** N * N, 2 * N + 1) == (F(2 * m, 4 * m + 1) if N % 2 == 0 else -F(2 * m + 1, 4 * m + 3))

# ar-l1-suites-pairs-impairs-ex2 Q3 : sqrt(n+1) - sqrt(n) -> 0.
assert limit(sqrt(n + 1) - sqrt(n), n, oo) == 0

# ar-l1-suites-pairs-impairs-ex4 (+ v1, v2) : u_{n+1} = 1 - u_n^2.
f = lambda t: 1 - t * t
assert expand(f(f(X)) - X - (-X * (X - 1) * (X**2 + X - 1))) == 0
alpha = (math.sqrt(5) - 1) / 2
assert abs(f(alpha) - alpha) < 1e-15
for u0, (lp, li) in [(F(1, 2), (0, 1)), (F(4, 5), (1, 0)), (F(1, 3), (0, 1))]:
    assert (u0 < alpha) == (lp == 0)
    u = float(u0)
    suite = [u]
    for m in range(400):
        u = f(u)
        suite.append(u)
    assert abs(suite[400] - lp) < 1e-12 and abs(suite[399] - li) < 1e-12
    # monotonie de (u_{2n}) : décroissante si u0 < alpha, croissante sinon
    pairs = suite[0::2][:30]
    assert all((pairs[j + 1] < pairs[j]) == (lp == 0) for j in range(len(pairs) - 1) if abs(pairs[j] - lp) > 1e-14)
# leurres : u_2 n'est pas u_0 (pas de suite stationnaire : 1/2 -> 3/4 -> 7/16)
assert f(F(1, 2)) == F(3, 4) and f(F(3, 4)) == F(7, 16) and f(F(4, 5)) == F(9, 25) and f(F(1, 3)) == F(8, 9)

# ---------------------------------------------------------------------------
# ar-l1-suites-cauchy (carte)
# ---------------------------------------------------------------------------
H = lambda m: sum(F(1, j) for j in range(1, m + 1))
# Contre-exemple : H_{2n} - H_n >= 1/2 et H_{n+p} - H_n <= p/(n+1).
assert all(H(2 * m) - H(m) >= F(1, 2) for m in range(1, 80))
assert all(H(m + p) - H(m) <= F(p, m + 1) for m in range(1, 40) for p in range(1, 6))
# Exemple : somme des 1/k^2 de q+1 à p <= 1/q - 1/p.
assert all(sum(F(1, j * j) for j in range(q + 1, p + 1)) <= F(1, q) - F(1, p) for q in range(1, 15) for p in range(q + 1, 40))
# Contre-exemple Q : |q_p - q_q| <= 10^-min(p,q) pour les approximations décimales de sqrt 2.
qd = [F(math.isqrt(2 * 10 ** (2 * j)), 10**j) for j in range(12)]
assert all(abs(qd[p] - qd[q]) <= F(1, 10 ** min(p, q)) for p in range(12) for q in range(12))
# Affirmation 3 : somme des 1/2^k pour k de q à p-1 < 2/2^q.
assert all(sum(F(1, 2**j) for j in range(q, p)) < F(2, 2**q) for q in range(10) for p in range(q, 30))

# ar-l1-suites-cauchy-ex1 (+ v1, v2) : plus petit N tel que |u_p - u_q| < eps pour p, q >= N.
# Pour une suite décroissante positive de limite 0, la condition équivaut à u_N <= eps.
for u, eps, N, debut, leurres in [(lambda m: F(1, m), F(1, 100), 100, 1, [99, 101, 50]),
                                   (lambda m: F(2, m + 1), F(1, 50), 99, 0, [100, 98, 49]),
                                   (lambda m: F(1, 3**m), F(1, 100), 5, 0, [4, 6, 100])]:
    plus_petit = next(m for m in range(debut, 1000) if u(m) <= eps)
    assert plus_petit == N and N not in leurres
    # N - 1 échoue : avec p = N - 1 et q grand, |u_p - u_q| dépasse eps
    assert any(u(N - 1) - u(q) >= eps for q in range(N, 20000))
assert 3**4 == 81 < 100 <= 243 == 3**5
assert F(1, 99) - F(1, 9900) == F(1, 100)  # cas limite cité dans la correction

# ar-l1-suites-cauchy-ex2 Q3 : u_{n+1} = cos(u_n)/2 converge vers l avec 2 l = cos l.
u = 0.0
for m in range(200):
    u = math.cos(u) / 2
assert abs(2 * u - math.cos(u)) < 1e-14

# ar-l1-suites-cauchy-ex4 Q2 : V_n <= 2 - 1/(n+1) ; Q3 : |u_{k+1} - u_k| = 1/(k+2) + 1/(k+1).
assert all(1 + sum(F(1, (j + 1) ** 2) for j in range(1, m + 1)) <= 2 - F(1, m + 1) for m in range(1, 60))
uq = lambda j: F((-1) ** j, j + 1)
assert all(abs(uq(j + 1) - uq(j)) == F(1, j + 2) + F(1, j + 1) for j in range(100))
u = 1.0
for m in range(5000):
    u = u + math.cos(m * u) / (m + 1) ** 2
assert abs(u) < 10  # suite bornée (convergente)

# ---------------------------------------------------------------------------
# ar-l1-suites-cesaro (carte)
# ---------------------------------------------------------------------------
# Contre-exemple « bornée » : u_n = 1 si 4^m <= n < 2.4^m, 0 sinon ; moyennes aux rangs 2.4^m - 1 et 4^(m+1) - 1.
ub = lambda j: 1 if any(4**m <= j < 2 * 4**m for m in range(12)) else 0
cumul, s = {}, 0
for j in range(1, 4**7):
    s += ub(j)
    cumul[j] = s
for m in range(6):
    assert cumul[2 * 4**m - 1] == cumul[4 ** (m + 1) - 1] == F(4 ** (m + 1) - 1, 3)
    assert F(cumul[4 ** (m + 1) - 1], 4 ** (m + 1) - 1) == F(1, 3)
assert abs(F(cumul[2 * 4**6 - 1], 2 * 4**6 - 1) - F(2, 3)) < F(1, 1000)
# Contre-exemple « réciproque » : sommes partielles de (-1)^n valent -1 ou 0.
assert all(sum((-1) ** j for j in range(1, m + 1)) in (-1, 0) for m in range(1, 50))
# Contre-exemple « écarts » : moyenne de sqrt k >= (1/2) sqrt(n/2).
for N in (10, 100, 1000):
    assert sum(math.sqrt(j) for j in range(1, N + 1)) / N >= 0.5 * math.sqrt(N / 2)
# Trou 2, leurre « monotone » : u_n = (-1)^n/n donne v_1 = -1, v_2 = -1/4, v_3 = -5/18.
v = [sum(F((-1) ** j, j) for j in range(1, m + 1)) / m for m in range(1, 4)]
assert v == [-1, F(-1, 4), F(-5, 18)] and v[1] > v[0] and v[2] < v[1]
# Exemple : u_{n+1} - u_n = 2 + 1/(n+1) -> 2.
assert limit(2 + 1 / (n + 1), n, oo) == 2

# ar-l1-suites-cesaro-ex1 (+ v1, v2) : limites des termes moyennés, et moyennes numériques.
for expr, l, leurres in [((2 * n + 1) / (n + 3), 2, [0, 1, oo]),
                         ((3 * n**2 + 1) / (n**2 + n), 3, [0, 1, oo]),
                         ((n - 1) / (2 * n + 5), Rational(1, 2), [0, 2, Rational(-1, 5)])]:
    assert limit(expr, n, oo) == l and l not in leurres
from math import fsum
assert abs(fsum((2 * j + 1) / (j + 3) for j in range(1, 10**5 + 1)) / 10**5 - 2) < 1e-3
assert abs(fsum((3 * j * j + 1) / (j * j + j) for j in range(1, 10**5 + 1)) / 10**5 - 3) < 1e-3
assert abs(fsum((j - 1) / (2 * j + 5) for j in range(1, 10**5 + 1)) / 10**5 - 0.5) < 1e-3

# ar-l1-suites-cesaro-ex2 (+ v1, v2) : u_n / n.
for incr, l, leurres in [((3 * n + 1) / (n + 2), 3, [0, Rational(1, 3), oo]),
                         ((n + 4) / (2 * n + 1), Rational(1, 2), [0, 2, oo]),
                         ((5 * n**2 - 1) / (n**2 + 1), 5, [0, Rational(1, 5), oo])]:
    assert limit(incr, n, oo) == l and l not in leurres
assert (5 * 0**2 - 1) / (0**2 + 1) == -1  # premier accroissement négatif (v2)
uu = 0.0
for j in range(10**5):
    uu += (3 * j + 1) / (j + 2)
assert abs(uu / 10**5 - 3) < 1e-3

# ar-l1-suites-cesaro-ex3 (+ v1, v2) : u_{n+1} = u_n + c/(d u_n), u_n/sqrt(n) -> sqrt(2c/d).
for c, d, l, leurres in [(1, 1, sqrt(2), [2, 1, 1 / sqrt(2)]), (2, 1, 2, [4, sqrt(2), 1]), (1, 2, 1, [Rational(1, 2), sqrt(2), 1 / sqrt(2)])]:
    assert simplify(expand((X + Rational(c, d) / X) ** 2 - X**2) - (2 * Rational(c, d) + Rational(c, d) ** 2 / X**2)) == 0
    assert sqrt(2 * Rational(c, d)) == l and all(simplify(l - t) != 0 for t in leurres)
    u = 1.0
    N = 10**6
    for j in range(N):
        u += c / (d * u)
    assert abs(u / math.sqrt(N) - float(l)) < 1e-2

# ar-l1-suites-cesaro-ex4 (+ v1, v2) : a_n^(1/n) -> limite du quotient a_{n+1}/a_n.
assert limit((1 + 1 / n) ** n, n, oo) == E
assert simplify((2 * n + 2) * (2 * n + 1) / (n + 1) ** 2 - 2 * (2 * n + 1) / (n + 1)) == 0
assert limit(2 * (2 * n + 1) / (n + 1), n, oo) == 4
for lna, l, leurres in [(lambda m: m * math.log(m) - math.lgamma(m + 1), math.e, [1, 1 / math.e]),
                        (lambda m: math.lgamma(2 * m + 1) - 2 * math.lgamma(m + 1), 4, [2, math.e, 1]),
                        (lambda m: m * math.log(m) - m * math.log(2) - math.lgamma(m + 1), math.e / 2, [2 * math.e, math.e, 0.5])]:
    m = 10**6
    approx = math.exp(lna(m) / m)
    assert abs(approx - l) < 1e-3 * l * 10
    assert all(abs(t - l) > 0.1 for t in leurres)
# Quotients exacts (petits n) : (n+1)^n/n^n, (2n+2)(2n+1)/(n+1)^2, (1/2)(1+1/n)^n.
for m in range(1, 20):
    a = lambda j: F(j**j, math.factorial(j))
    b = lambda j: F(math.factorial(2 * j), math.factorial(j) ** 2)
    c2 = lambda j: F(j**j, 2**j * math.factorial(j))
    assert a(m + 1) / a(m) == F((m + 1) ** m, m**m)
    assert b(m + 1) / b(m) == F((2 * m + 2) * (2 * m + 1), (m + 1) ** 2)
    assert c2(m + 1) / c2(m) == F((m + 1) ** m, 2 * m**m)

# ---------------------------------------------------------------------------
# ar-l1-suites-recurrence-lineaire (carte)
# ---------------------------------------------------------------------------
phi, psi = (1 + sqrt(5)) / 2, (1 - sqrt(5)) / 2
fib = suite2(0, 1, 1, 1, 30)
assert all(simplify((phi**j - psi**j) / sqrt(5) - fib[j]) == 0 for j in range(15))  # Binet
assert abs(float(psi)) < 1 and simplify(psi**2 - psi - 1) == 0  # affirmation 2 (fausse)
# Contre-exemple a = b = 0 : u = 0, 1, 0, 0, ... ; (l + m n) 0^n vaut 0 en n = 1.
assert suite2(0, 1, 0, 0, 5) == [0, 1, 0, 0, 0, 0]
# Racine double : n 2^n vérifie u_{n+2} = 4u_{n+1} - 4u_n.
assert all((j + 2) * 2 ** (j + 2) == 4 * (j + 1) * 2 ** (j + 1) - 4 * j * 2**j for j in range(30))
# Démonstration (cas double) : r0^n [ n (r0^2 - a r0 - b) + r0 (2 r0 - a) ] avec r0 = a/2, b = -a^2/4.
a_ = symbols("a", real=True)
r0 = a_ / 2
assert simplify(r0**2 - a_ * r0 + a_**2 / 4) == 0 and simplify(r0 * (2 * r0 - a_)) == 0
# Second membre : u_{n+2} = u_{n+1} + 2u_n + 1, u_0 = u_1 = 0 -> u_2 = 1, alors que la forme homogène impose 0.
assert 0 + 2 * 0 + 1 == 1
# Affirmation 4 : 2^n cos(n pi/3) vérifie u_{n+2} = 2u_{n+1} - 4u_n, et u_{6m} = 2^{6m}.
assert all(simplify(2 ** (j + 2) * cos((j + 2) * pi / 3) - 2 * 2 ** (j + 1) * cos((j + 1) * pi / 3) + 4 * 2**j * cos(j * pi / 3)) == 0 for j in range(12))
assert simplify(((1 + I * sqrt(3)) - 2 * exp(I * pi / 3)).rewrite(cos)) == 0
# Affirmation 3 : u_{n+2} = 2u_{n+1} - u_n : suite arithmétique.
assert suite2(F(3), F(5), 2, -1, 10) == [3 + 2 * j for j in range(11)]

# ar-l1-suites-recurrence-lineaire-ex1..ex3 (+ v1, v2) : formules contre la récurrence ; leurres faux.
j = symbols("j", integer=True, nonnegative=True)
cas = [
    # (u0, u1, a, b, formule, leurres)
    (1, 4, 5, -6, 2 * 3**j - 2**j, [3**j + 2**j, 3 * 2**j - 2 * 3**j, 4 * 3**j - 3 * 2**j]),
    (0, 1, 1, 6, (3**j - (-2) ** j) / 5, [(3**j + (-2) ** j) / 5, (3**j - 2**j) / 5, (2**j - (-3) ** j) / 5]),
    (2, 3, 3, -2, 2**j + 1, [2 ** (j + 1), j + 2, 3 * 2**j - 1]),
    (1, 4, 4, -4, (j + 1) * 2**j, [4**j, (2 * j + 1) * 2**j, (3 * j + 2) * 2 ** (j - 1)]),
    (2, 3, 6, -9, (2 - j) * 3**j, [(2 + j) * 3**j, (2 - 3 * j) * 3**j, 3**j + 1]),
    (3, -1, -2, -1, (3 - 2 * j) * (-1) ** j, [(3 + 2 * j) * (-1) ** j, (3 - 4 * j) * (-1) ** j, 3 * (-1) ** j + 2 * j]),
    (1, 0, 2, -2, 2 ** (j / 2) * (cos(j * pi / 4) - sin(j * pi / 4)),
     [cos(j * pi / 4) - sin(j * pi / 4), 2 ** (j / 2) * (cos(j * pi / 4) + sin(j * pi / 4)), 2 ** (j / 2) * cos(j * pi / 4)]),
    (2, 2, 2, -4, 2 ** (j + 1) * cos(j * pi / 3),
     [2 * cos(j * pi / 3), 2 ** (j + 1) * cos(j * pi / 6), 2**j * (2 * cos(j * pi / 3) + sin(j * pi / 3))]),
    (0, 3, -1, -1, 2 * sqrt(3) * sin(2 * j * pi / 3),
     [2 * sqrt(3) * sin(j * pi / 3), 3 * sin(2 * j * pi / 3), 2 * sqrt(3) * cos(2 * j * pi / 3)]),
]
for u0, u1, a, b, formule, leurres in cas:
    u = suite2(u0, u1, a, b, 20)
    assert all(simplify(formule.subs(j, t) - u[t]) == 0 for t in range(21))
    for l in leurres:
        assert any(simplify(l.subs(j, t) - u[t]) != 0 for t in range(4))

# ar-l1-suites-recurrence-lineaire-ex4 (+ v1, v2) : u_0 = 1, bornée ssi le coefficient de la grande racine est nul.
u1 = symbols("u1")
for a, b, (r1, r2), bonne, autres in [(1, 1, (phi, psi), psi, [phi, 0, -1, 1, Rational(1, 2)]),
                                      (2, 3, (3, -1), -1, [3, 1, 0]),
                                      (Rational(3, 2), 1, (2, Rational(-1, 2)), Rational(-1, 2), [2, Rational(1, 2), 0])]:
    assert simplify(r1**2 - a * r1 - b) == 0 and simplify(r2**2 - a * r2 - b) == 0
    assert abs(float(r1)) > 1 >= abs(float(r2))
    lam = (u1 - r2) / (r1 - r2)  # coefficient de r1^n quand u_0 = 1
    assert simplify(lam.subs(u1, bonne)) == 0
    for t in autres:  # pour ces valeurs de u_1 (leurres), lambda != 0 : suite non bornée
        assert simplify(lam.subs(u1, t)) != 0
    # numériquement : avec u_1 = bonne, la suite reste bornée ; sinon elle explose
    uu = suite2(1.0, float(bonne), float(a), float(b), 30)
    assert max(abs(t) for t in uu) <= 1.0 + 1e-6
    uu = suite2(1.0, float(autres[0]), float(a), float(b), 60)
    assert max(abs(t) for t in uu) > 1e5
# v1 : lambda = (1 + u1)/4 ; v2 : lambda = (2 u1 + 1)/5 ; principal : lambda = (u1 - psi)/sqrt 5.
assert simplify((u1 + 1) / 4 - (u1 - (-1)) / (3 - (-1))) == 0
assert simplify((2 * u1 + 1) / 5 - (u1 + Rational(1, 2)) / (2 + Rational(1, 2))) == 0
assert simplify(phi - psi - sqrt(5)) == 0

print("ar_l1_suites : toutes les vérifications passent.")
