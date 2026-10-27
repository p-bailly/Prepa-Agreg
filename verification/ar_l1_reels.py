"""Vérifications calculatoires du chapitre « ar-l1-reels » (analyse réelle, L1) :
cartes et exercices de content/*/analyse-reelle/l1-ar-l1-reels.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise : v1, v2).
Lancer : python3 verification/ar_l1_reels.py (doit finir sans erreur).
"""

import itertools
import math
from fractions import Fraction as F

from sympy import Interval, Rational, S, Symbol, Union, floor, pi, solve_univariate_inequality, sqrt

x = Symbol("x", real=True)


def solutions_abs(a, b, r):
    """Ensemble des x réels tels que |a x + b| <= r."""
    return solve_univariate_inequality(abs(a * x + b) <= r, x, relational=False)


def grille(a, b, pas=F(1, 37)):
    """Rationnels de ]a, b[ sur une grille."""
    t = a + pas
    while t < b:
        yield t
        t += pas


# ---------------------------------------------------------------------------
# ar-l1-reels-triangulaire (carte)
# ---------------------------------------------------------------------------
# Exemple : |x-2| <= 1, |y+1| <= 3 => -3 <= x+y <= 5 (bornes atteintes en x=1,y=-4 et x=3,y=2).
assert 1 + (-4) == -3 and 3 + 2 == 5
assert abs(-3 + (-4)) == 7 == abs(-3) + abs(-4) and abs(-3 + 4) == 1 < 7
# Contre-exemples : stricte (2, -5) ; égalité sans x = y ni |x| = |y| (1, 4) ; x+y >= 0 sans égalité (2, -1) ;
# on ne soustrait pas (0, 1).
assert abs(2 - 5) == 3 < abs(2) + abs(-5) == 7
assert abs(1 + 4) == abs(1) + abs(4) and abs(1) != abs(4)
assert 2 + (-1) >= 0 and abs(2 - 1) < abs(2) + abs(-1)
assert not (abs(0 - 1) <= abs(0) - abs(1))
# Cas d'égalité « |x+y| = |x|+|y| <=> xy >= 0 » vérifié sur une grille d'entiers.
for a, b in itertools.product(range(-6, 7), repeat=2):
    assert (abs(a + b) == abs(a) + abs(b)) == (a * b >= 0)
    assert abs(a + b) <= abs(a) + abs(b) and abs(abs(a) - abs(b)) <= abs(a - b)
# Pièges : « |x| = |y| » (x = 1, y = -1 : |0| < 2) ; « x, y > 0 » (x = y = -1 : égalité).
assert abs(1 - 1) < 2 and abs(-1 - 1) == abs(-1) + abs(-1)
# Trou 1 (leurres faux) : « >= » et « = » faux pour (1, -1) ; « < » faux pour (1, 1).
assert not (abs(1 - 1) >= 2) and not (abs(1 + 1) < 2)
# Trou 3 (leurres faux) : |x|-|y| (0, 1) ; |x-y|/2 (1, 0) ; min(|x|,|y|) (1, 0).
assert not (abs(0 - 1) <= abs(0) - abs(1))
assert not (abs(abs(1) - abs(0)) <= F(abs(1 - 0), 2))
assert not (abs(abs(1) - abs(0)) <= min(1, 0))
# Affirmations : 2 fausse (x = y = -1), 4 fausse (x = 1, y = -1 : 0 != 2).
assert abs(abs(1) - abs(-1)) == 0 != abs(1 - (-1))

# ar-l1-reels-triangulaire-ex1 (+ v1, v2) : |2x-3| <= 5 ; |3x+2| <= 7 ; |4x-1| <= 9.
for (a, b, r), bonne, leurres in [
    ((2, -3, 5), Interval(-1, 4), [Interval(-2, 8), Union(Interval(-S.Infinity, -1), Interval(4, S.Infinity)), Interval(-4, 1)]),
    ((3, 2, 7), Interval(-3, Rational(5, 3)), [Interval(-9, 5), Union(Interval(-S.Infinity, -3), Interval(Rational(5, 3), S.Infinity)), Interval(Rational(-5, 3), 3)]),
    ((4, -1, 9), Interval(-2, Rational(5, 2)), [Interval(-8, 10), Union(Interval(-S.Infinity, -2), Interval(Rational(5, 2), S.Infinity)), Interval(Rational(-5, 2), 2)]),
]:
    sol = solutions_abs(a, b, r)
    assert sol == bonne
    for l in leurres:
        assert l != sol

# ar-l1-reels-triangulaire-ex2 (+ v1, v2) : maximum d'une expression affine sur un pavé fermé (atteint en un sommet).
for (cx, rx, cy, ry, f), bonne, leurres in [
    ((2, F(1, 10), 5, F(1, 20), lambda X, Y: 2 * X - Y + 1), F(1, 4), [F(3, 20), F(3, 10), 0]),
    ((-1, F(1, 5), 3, F(1, 10), lambda X, Y: 3 * X + Y), F(7, 10), [F(3, 10), F(1, 2), 0]),
    ((1, F(1, 4), -2, F(1, 8), lambda X, Y: X - 2 * Y - 5), F(1, 2), [F(3, 8), F(1, 4), 0]),
]:
    assert f(cx, cy) == 0  # l'expression s'annule au centre (combinaison sans constante)
    sommets = [abs(f(cx + sx * rx, cy + sy * ry)) for sx in (-1, 1) for sy in (-1, 1)]
    assert max(sommets) == bonne
    assert all(l != bonne for l in leurres)
# Points d'égalité donnés dans les corrections.
assert 2 * F(21, 10) - F(99, 20) + 1 == F(1, 4)
assert 3 * F(-4, 5) + F(31, 10) == F(7, 10)
assert F(5, 4) - 2 * F(-17, 8) - 5 == F(1, 2)

# ar-l1-reels-triangulaire-ex3 : |x+y|/(1+|x+y|) <= |x|/(1+|x|) + |y|/(1+|y|) sur une grille.
g = lambda t: F(abs(t)) / (1 + abs(t))
for a, b in itertools.product([F(k, 3) for k in range(-12, 13)], repeat=2):
    assert g(a + b) <= g(a) + g(b)

# ar-l1-reels-triangulaire-ex4 (+ v1, v2) : minimum de |x-a|+|x-b|+|x-c| (convexe, affine par morceaux).
for (pts), (vmin, xmin), mauvais in [
    ((1, 4, 6), (5, 4), [(F(16, 3), F(11, 3)), (2, 4)]),
    ((-2, 1, 5), (7, 1), [(F(22, 3), F(4, 3)), (4, 1)]),
    ((0, 3, 9), (9, 3), [(10, 4), (6, 3)]),
]:
    f = lambda t, p=pts: sum(abs(t - q) for q in p)
    assert min(f(q) for q in pts) == vmin == f(xmin)
    # unicité : f > vmin ailleurs (sommets et points d'une grille fine)
    assert all(f(t) > vmin for t in grille(F(pts[0] - 3), F(pts[2] + 3), F(1, 7)) if t != xmin)
    # leurre « vmin atteint sur tout [a, c] » faux : f(a) != vmin
    assert f(pts[0]) != vmin
    # leurres (valeur, point) : la valeur en la moyenne est > vmin ; l'autre valeur n'est pas le minimum
    (vm, xm), (vb, xb) = mauvais
    assert xm == F(sum(pts), 3) and f(xm) == vm and vm > vmin
    assert vb != vmin

# ar-l1-reels-triangulaire-ex5 Q3 : somme des 1/2^k pour k < n = 2 - 1/2^(n-1) < 2.
for n in range(1, 40):
    s = sum(F(1, 2**k) for k in range(n))
    assert s == 2 - F(1, 2 ** (n - 1)) and s < 2

# ---------------------------------------------------------------------------
# ar-l1-reels-densite (carte)
# ---------------------------------------------------------------------------
assert floor(pi) == 3 and math.floor(-2.5) == -3 and math.floor(-3) == -3
assert [math.floor(10**k * math.sqrt(2)) for k in range(4)] == [1, 14, 141, 1414]
assert sqrt(2) < Rational(1415, 1000) < Rational(142, 100)
# Contre-exemples : x = -1/2 n'a pas de partie entière dans N ; aucun entier dans ]1/10, 1/5[.
assert math.floor(F(-1, 2)) == -1 and not any(F(n) <= F(-1, 2) for n in range(0, 10))
assert not any(F(1, 10) < n < F(1, 5) for n in range(-5, 5))
# Affirmations 1 et 2 fausses : x = y = 1/2.
assert math.floor(F(1, 2) + F(1, 2)) == 1 != math.floor(F(1, 2)) + math.floor(F(1, 2))
assert math.floor(-F(1, 2)) == -1 != -math.floor(F(1, 2))
# Affirmation 4 vraie : floor(x + k) = floor(x) + k (grille).
assert all(math.floor(F(p, 7) + k) == math.floor(F(p, 7)) + k for p in range(-30, 30) for k in range(-5, 6))
# Trou 2 (leurres faux en x = 1) : n < x <= n+1 donne n = 0 ; x-1 <= n < x donne n = 0 ; n <= x <= n+1 : n = 0 ou 1.
assert [n for n in range(-3, 4) if n < 1 <= n + 1] == [0]
assert [n for n in range(-3, 4) if 1 - 1 <= n < 1] == [0]
assert [n for n in range(-3, 4) if n <= 1 <= n + 1] == [0, 1]
# Démonstration de la densité : p = floor(q a) + 1 avec q (b - a) > 1 (test sur quelques couples).
for a, b in [(F(1, 3), F(17, 50)), (F(-7, 5), F(-139, 100)), (F(2), F(2001, 1000))]:
    q = math.floor(1 / (b - a)) + 1
    p = math.floor(q * a) + 1
    assert q * (b - a) > 1 and a < F(p, q) < b

# ar-l1-reels-densite-ex1 (+ v1, v2) : floor(2x+1) = 4 ; floor(3x-2) = -1 ; floor((x+3)/2) = -2.
eps = F(1, 10**6)
for f, k, (lo, hi), leurres in [
    (lambda t: 2 * t + 1, 4, (F(3, 2), F(2)), [(F(2), F(5, 2))]),
    (lambda t: 3 * t - 2, -1, (F(1, 3), F(2, 3)), [(F(-1, 3), F(0))]),
    (lambda t: (t + 3) / 2, -2, (F(-7), F(-5)), [(F(-4), F(-2))]),
]:
    sol = lambda t: math.floor(f(t)) == k
    # solution [lo, hi[ : lo inclus, hi exclu, intérieur inclus, extérieur proche exclu
    assert sol(lo) and not sol(hi) and sol(hi - eps) and not sol(lo - eps)
    assert all(sol(t) for t in grille(lo, hi, (hi - lo) / 50))
    # leurres « ]lo, hi] » et « [lo, hi] » : contiennent hi, qui n'est pas solution ; leurre décalé : son bord gauche n'est pas solution
    for a, b in leurres:
        assert not sol(a)

# ar-l1-reels-densite-ex2 (+ v1, v2) : somme des floor(sqrt k).
for N, bonne, leurres in [(20, 54, [50, 57, 70]), (30, 100, [95, 106, 70]), (12, 25, [24, 28, 34])]:
    s = sum(math.isqrt(k) for k in range(1, N + 1))
    assert s == bonne and bonne not in leurres
assert 1 * 3 + 2 * 5 + 3 * 7 + 4 * 5 == 54 and 3 + 10 + 21 + 36 + 30 == 100 and 3 + 10 + 12 == 25

# ar-l1-reels-densite-ex3 : floor(x) + floor(x + 1/2) = floor(2x) (grille de rationnels).
for p in range(-400, 400):
    t = F(p, 37)
    assert math.floor(t) + math.floor(t + F(1, 2)) == math.floor(2 * t)

# ar-l1-reels-densite-ex4 Q1 : r_n <= x < s_n = r_n + 10^-n.
for t in [math.pi, -math.e, 2**0.5, -3.0]:
    for n in range(8):
        r = math.floor(10**n * t) / 10**n
        assert r <= t < r + 10**-n + 1e-12

# ---------------------------------------------------------------------------
# ar-l1-reels-sup-sequentielle (carte)
# ---------------------------------------------------------------------------
# Exemple : A = {1/p + 1/q} : éléments > 0, 2/n -> 0, maximum 2 en p = q = 1.
elts = [F(1, p) + F(1, q) for p in range(1, 60) for q in range(1, 60)]
assert all(e > 0 for e in elts) and max(elts) == 2 and min(elts) == F(2, 59)
# Contre-exemples : A = [0,1], M = 1/2 (non majorant) ; M = 2 ; a_n = 1 - 1/(n+1) croissante dans A.
assert F(1, 2) < 1 and all(1 - F(1, n + 1) < 1 - F(1, n + 2) <= 1 for n in range(50))
# Affirmation 2 fausse : 1/(n+1) dans ]0,1], limite 0 hors de A.
assert all(0 < F(1, n + 1) <= 1 for n in range(100))
# Trous : leurre « diverge » faux pour A = {0, 1} (suite 0,1,0,1 divergente dans une partie bornée) ;
# leurre « n'est pas bornée » faux pour A = ]-inf, 0] (suite -n non bornée, A majorée).
assert {0, 1} == {n % 2 for n in range(10)}

# ar-l1-reels-sup-sequentielle-ex1 (+ v1, v2) : bornes d'ensembles à termes alternés (n de 1 à 4000).
NMAX = 4000
A0 = [(-1) ** n * (1 - F(1, n)) for n in range(1, NMAX)]
assert all(-1 < a < 1 for a in A0)
assert 1 - max(A0) < F(1, 1000) and min(A0) + 1 < F(1, 1000)
assert max(A0[:4]) == F(3, 4) > F(1, 2)  # leurre « sup = 1/2 atteint » faux (a_4 = 3/4)
A1 = [(-1) ** n + F(1, n) for n in range(1, NMAX)]
assert max(A1) == F(3, 2) == A1[1] and all(a > -1 for a in A1) and min(A1) + 1 < F(1, 1000)
assert F(3, 2) > 1  # leurre « sup = 1 » faux
A2 = [2 + F((-1) ** n, n) for n in range(1, NMAX)]
assert max(A2) == F(5, 2) == A2[1] and min(A2) == 1 == A2[0]
assert any(a > 2 for a in A2) and any(a < 2 for a in A2)  # leurre « sup = inf = 2 » faux

# ar-l1-reels-sup-sequentielle-ex3 (+ v1, v2) : A = {xy, x dans ]a,b[, y dans ]c,d[}.
for (a, b, c, d), (s, i), leurres in [
    ((-1, 2, -3, 1), (3, -6), [(2, -1), (2, -6)]),
    ((-2, 1, -1, 4), (4, -8), [(2, -1), (4, -1)]),
    ((-3, 2, -2, 3), (6, -9), [(6, -4), (9, -9)]),
]:
    coins = [p * q for p in (a, b) for q in (c, d)]
    assert max(coins) == s and min(coins) == i
    pas = F(1, 23)
    prods = [p * q for p in grille(F(a), F(b), pas) for q in grille(F(c), F(d), pas)]
    assert all(i < v < s for v in prods)  # bornes strictes : non atteintes
    assert s - max(prods) < F(1, 2) and min(prods) - i < F(1, 2)
    for ls, li in leurres:  # (sup, inf) est unique et vérifié ci-dessus : chaque leurre en diffère
        assert (ls, li) != (s, i)
        assert max(prods) > ls or min(prods) < li or ls > s  # et il est contredit directement

# ar-l1-reels-sup-sequentielle-ex2 Q2 : f = indicatrice de [1, +inf[, A = [0, 1[ : sup f(A) = 0 != f(1) = 1.
ind = lambda t: 1 if t >= 1 else 0
assert max(ind(t) for t in grille(F(0), F(1), F(1, 1000))) == 0 and ind(1) == 1

print("ar_l1_reels : toutes les vérifications passent.")
