"""Vérifications calculatoires de l'enrichissement du chapitre « ar-l1-reels » (analyse réelle, L1) :
questions de content/qcm/ar-l1-reels.json et exercices ajoutés à content/exercices/analyse-reelle/l1-ar-l1-reels.json.

Chaque assert cite l'identifiant de la question ou de l'exercice (variante comprise : v1, v2).
Lancer : python3 verification/enrichi_ar_l1_reels.py (doit finir sans erreur).
"""

import random
from fractions import Fraction as F

from sympy import (Abs, Eq, Interval, Rational, S, Union, expand, factorial, floor, limit, multiplicity, oo, simplify,
                   solveset, sqrt, symbols)

x = symbols("x", real=True)
n = symbols("n", positive=True, integer=True)

# ---------------------------------------------------------------------------
# ar-l1-reels-q1 : partie entière de -7/2 et de 7/2.
# ---------------------------------------------------------------------------
assert floor(Rational(-7, 2)) == -4 and floor(Rational(7, 2)) == 3
assert floor(Rational(-7, 2)) + floor(Rational(7, 2)) == -1  # bonne réponse -1 ; leurres 0, 1, -2 faux

# ---------------------------------------------------------------------------
# ar-l1-reels-q2 : |x-1| + |x+2| = 3 a pour ensemble de solutions [-2, 1].
# ---------------------------------------------------------------------------
sol = solveset(Eq(Abs(x - 1) + Abs(x + 2), 3), x, S.Reals)
assert sol == Interval(-2, 1)
# leurres : {-2,1} et {0} sont strictement contenus (pas l'ensemble entier) ; la réunion hors de ]-2,1[ est fausse
assert sol != S.EmptySet and Rational(1, 2) in sol and 2 not in sol

# ---------------------------------------------------------------------------
# ar-l1-reels-q3 : A = {(-1)^n + 1/n, n >= 1} : sup 3/2 atteint (n = 2), inf -1 non atteint.
# ---------------------------------------------------------------------------
valeurs = [F((-1) ** k) + F(1, k) for k in range(1, 2001)]
assert max(valeurs) == F(3, 2) and valeurs[1] == F(3, 2)
assert all(v > -1 for v in valeurs)  # -1 n'est jamais atteint
assert F(-1) + F(1, 1999) - (-1) < F(1, 1000)  # les rangs impairs s'approchent de -1
assert valeurs[0] == 0  # n = 1 donne 0 (et non 2) : leurre « sup A = 2 » faux
# leurre « inf A = 0 atteint » faux : -1 + 1/3 = -2/3 < 0
assert valeurs[2] == F(-2, 3)

# ---------------------------------------------------------------------------
# ar-l1-reels-q4 (faux) : sup(AB) = sup A sup B est faux en général : A = B = {-1, 0}.
# ---------------------------------------------------------------------------
A = [-1, 0]
AB = [a * b for a in A for b in A]
assert max(AB) == 1 and max(A) * max(A) == 0

# ---------------------------------------------------------------------------
# ar-l1-reels-q5 (vrai) : floor(x) + floor(-x) vaut 0 si x entier, -1 sinon.
# ---------------------------------------------------------------------------
for v in [Rational(7, 3), Rational(-5, 2), sqrt(2), -sqrt(3), S(4), S(-6), S(0)]:
    s = floor(v) + floor(-v)
    assert s == (0 if v.is_integer else -1)

# ---------------------------------------------------------------------------
# ar-l1-reels-q6 (faux) : sqrt(x^2) = x est faux pour x = -1 (sqrt(x^2) = |x|).
# ---------------------------------------------------------------------------
assert sqrt(S(-1) ** 2) == 1 != -1
assert simplify(sqrt(x**2) - Abs(x)) == 0

# ---------------------------------------------------------------------------
# ar-l1-reels-q7 : floor(x^2) = 2 <=> 2 <= x^2 < 3 <=> x dans ]-sqrt3, -sqrt2] u [sqrt2, sqrt3[.
# ---------------------------------------------------------------------------
sol7 = solveset(2 <= x**2, x, S.Reals).intersect(solveset(x**2 < 3, x, S.Reals))
attendu = Union(Interval.Lopen(-sqrt(3), -sqrt(2)), Interval.Ropen(sqrt(2), sqrt(3)))
assert simplify(sol7.symmetric_difference(attendu)) == S.EmptySet
# leurres : [sqrt2, sqrt3[ oublie -1,5 ; [sqrt2, sqrt3] contient sqrt3 (floor(3) = 3) ; ]-sqrt3,-sqrt2] u [sqrt2,sqrt3] idem
assert Rational(-3, 2) in attendu and floor(Rational(9, 4)) == 2
assert floor(sqrt(3) ** 2) == 3

# ---------------------------------------------------------------------------
# ar-l1-reels-q8 : A = {2^-n + 3^-m, n, m dans N} : sup 2 atteint en (0,0), inf 0 non atteint.
# ---------------------------------------------------------------------------
vals = [F(1, 2**a) + F(1, 3**b) for a in range(40) for b in range(40)]
assert max(vals) == 2 and min(vals) > 0
assert limit(2**-n + 3**-n, n, oo) == 0

# ---------------------------------------------------------------------------
# ar-l1-borne-sup-ex4 : sous-groupes de R ; Z + sqrt2 Z contient (sqrt2 - 1)^k -> 0, de la forme a + b sqrt2.
# ---------------------------------------------------------------------------
for k in range(1, 8):
    e = expand((sqrt(2) - 1) ** k)
    b = e.coeff(sqrt(2))
    a = expand(e - b * sqrt(2))
    assert a.is_integer and b.is_integer and 0 < e.evalf() < 1
assert (sqrt(2) - 1).evalf() ** 7 < 0.01
# question 3 : n = floor(x/g) + 1 vérifie x < n g <= x + g < y (essais numériques exacts)
random.seed(1)
for _ in range(200):
    xx = F(random.randint(-1000, 1000), random.randint(1, 50))
    yy = xx + F(random.randint(1, 100), random.randint(1, 100))
    g = (yy - xx) / random.randint(2, 9)
    m = (xx / g).__floor__() + 1
    assert xx < m * g < yy

# ---------------------------------------------------------------------------
# ar-l1-reels-densite-ex5 : formule de Legendre.
# ---------------------------------------------------------------------------
def legendre(N, p):
    s, q = 0, p
    while q <= N:
        s += N // q
        q *= p
    return s


# question 1 : nombre de multiples de d dans {1..n} = floor(n/d)
for N in range(1, 60):
    for d in range(1, 15):
        assert sum(1 for k in range(1, N + 1) if k % d == 0) == N // d
# exercice : zéros de 100! = v5(100!) = 20 + 4 = 24 (v2 >= v5) ; leurres 20, 25, 22
assert legendre(100, 5) == 24 == multiplicity(5, factorial(100)) and legendre(100, 2) >= 24
assert multiplicity(10, factorial(100)) == 24 and str(factorial(100)).endswith("0" * 24)
assert str(factorial(100))[-25] != "0"
assert all(v != 24 for v in (20, 25, 22))
# v1 : plus grand k tel que 3^k divise 100! : 33 + 11 + 3 + 1 = 48 ; leurres 33, 44, 47
assert legendre(100, 3) == 48 == multiplicity(3, factorial(100))
assert all(v != 48 for v in (33, 44, 47))
# v2 : zéros de 250! : 50 + 10 + 2 = 62 ; leurres 50, 60, 63
assert legendre(250, 5) == 62 == multiplicity(10, factorial(250))
assert all(v != 62 for v in (50, 60, 63))

# ---------------------------------------------------------------------------
# ar-l1-reels-sup-sequentielle-ex5 et variantes : bornes d'ensembles indexés par m, n >= 1.
# ---------------------------------------------------------------------------
m_, n_ = symbols("m n", positive=True)
# exercice : mn/(m^2+n^2) <= 1/2 avec égalité si m = n ; > 0 ; m = 1, n -> oo donne 0 : sup 1/2 atteint, inf 0 non atteint
assert simplify(Rational(1, 2) - m_ * n_ / (m_**2 + n_**2) - (m_ - n_) ** 2 / (2 * (m_**2 + n_**2))) == 0
assert limit(n_ / (1 + n_**2), n_, oo) == 0
E0 = [F(a * b, a * a + b * b) for a in range(1, 80) for b in range(1, 80)]
assert max(E0) == F(1, 2) and min(E0) > 0
# v1 : (m+n)^2/(m^2+n^2) = 1 + 2mn/(m^2+n^2) dans ]1, 2] ; sup 2 atteint, inf 1 non atteint
assert simplify((m_ + n_) ** 2 / (m_**2 + n_**2) - 1 - 2 * m_ * n_ / (m_**2 + n_**2)) == 0
assert limit((1 + n_) ** 2 / (1 + n_**2), n_, oo) == 1
E1 = [F((a + b) ** 2, a * a + b * b) for a in range(1, 80) for b in range(1, 80)]
assert max(E1) == 2 and min(E1) > 1
# v2 : (m-n)/(m+n) dans ]-1, 1[ ; sup 1 et inf -1 non atteints
assert limit((n_ - 1) / (n_ + 1), n_, oo) == 1 and limit((1 - n_) / (1 + n_), n_, oo) == -1
E2 = [F(a - b, a + b) for a in range(1, 80) for b in range(1, 80)]
assert max(E2) < 1 and min(E2) > -1

print("enrichi_ar_l1_reels : toutes les vérifications passent.")
