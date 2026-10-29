"""Vérifications calculatoires de l'enrichissement du chapitre « ar-l1-continuite » (analyse réelle, L1) :
questions de content/qcm/ar-l1-continuite.json et exercices ajoutés à
content/exercices/analyse-reelle/l1-ar-l1-continuite.json.

Chaque assert cite l'identifiant de la question ou de l'exercice (variante comprise : v1, v2).
Lancer : python3 verification/enrichi_ar_l1_continuite.py (doit finir sans erreur).
"""

import math

from sympy import Interval, Rational, S, atan, diff, exp, limit, oo, pi, simplify, sin, solve, sqrt, symbols
from sympy.calculus.util import function_range

x = symbols("x", real=True)
y = symbols("y", real=True)
n = symbols("n", positive=True, integer=True)

# ---------------------------------------------------------------------------
# ar-l1-continuite-q1 : image de ]-1, 2] par x^2 - 2x : [-1, 3[.
# ---------------------------------------------------------------------------
img = function_range(x**2 - 2 * x, x, Interval.Lopen(-1, 2))
assert img == Interval.Ropen(-1, 3)
for leurre in (Interval.Ropen(0, 3), Interval(-1, 3), Interval.open(-1, 3)):
    assert leurre != img

# ---------------------------------------------------------------------------
# ar-l1-continuite-q2 (faux) : indicatrice de {0} : f(u_n) bornée pour toute suite, f non continue en 0.
# ---------------------------------------------------------------------------
ind = lambda t: 1 if t == 0 else 0
assert all(ind(1 / k) == 0 for k in range(1, 100)) and ind(0) == 1

# ---------------------------------------------------------------------------
# ar-l1-continuite-q3 (vrai) : f(x) = f(2x) et f continue => f(x) = f(x/2^n) -> f(0).
# ---------------------------------------------------------------------------
assert limit(Rational(7) / 2**n, n, oo) == 0

# ---------------------------------------------------------------------------
# ar-l1-continuite-q4 (vrai) : un polynôme de degré 3 tend vers -oo et +oo (coefficient dominant > 0).
# ---------------------------------------------------------------------------
P = 2 * x**3 - 5 * x + 7
assert limit(P, x, oo) == oo and limit(P, x, -oo) == -oo

# ---------------------------------------------------------------------------
# ar-l1-continuite-q5 : f(x) = x/(1+|x|), f^-1(y) = y/(1-|y|) sur ]-1, 1[.
# ---------------------------------------------------------------------------
yp = symbols("yp", positive=True)
xp = symbols("xp", positive=True)
g = lambda t, s: t / (1 - s * t)  # s = signe de t : |t| = s t
f = lambda t, s: t / (1 + s * t)
assert simplify(f(g(yp, 1), 1) - yp) == 0 and simplify(f(g(-yp, -1), -1) + yp) == 0  # f(f^-1(y)) = y (signe conservé)
assert simplify(g(f(xp, 1), 1) - xp) == 0 and simplify(g(f(-xp, -1), -1) + xp) == 0
assert function_range(xp / (1 + xp), xp, Interval.open(0, oo)) == Interval.open(0, 1)
# leurres testés en y = 1/3 et y = -1/3 : f^-1(1/3) = 1/2 et f^-1(-1/3) = -1/2
for yy, attendu in ((Rational(1, 3), Rational(1, 2)), (Rational(-1, 3), Rational(-1, 2))):
    s_ = 1 if yy > 0 else -1
    assert f(attendu, s_) == yy and yy / (1 - abs(yy)) == attendu  # bonne réponse y/(1-|y|)
    assert yy / (1 + abs(yy)) != attendu  # leurre y/(1+|y|)
    assert (1 - abs(yy)) / yy != attendu  # leurre (1-|y|)/y
assert abs(Rational(-1, 3)) / (1 - Rational(-1, 3)) != Rational(-1, 2)  # leurre |y|/(1-y) (positif en y < 0)

# ---------------------------------------------------------------------------
# ar-l1-continuite-q6 : seule sqrt est uniformément continue (sur [0,+oo[) parmi les quatre.
# ---------------------------------------------------------------------------
# |sqrt x - sqrt y| <= sqrt|x-y| (voir heine-ex5 plus bas). x^2 : x_n = n + 1/n, y_n = n.
assert limit((n + 1 / n) ** 2 - n**2, n, oo) == 2
# 1/x sur ]0,1] : x_n = 1/n, y_n = 1/(n+1) ; écart -> 0, images à distance 1
assert limit(1 / n - 1 / (n + 1), n, oo) == 0
# sin(1/x) : x_n = 1/(2n pi + pi/2), y_n = 1/(2n pi) : sin vaut 1 et 0
assert limit(1 / (2 * n * pi) - 1 / (2 * n * pi + pi / 2), n, oo) == 0
assert sin(2 * pi + pi / 2) == 1 and sin(2 * pi) == 0

# ---------------------------------------------------------------------------
# ar-l1-continuite-q7 : arctan n'atteint pas sa borne sup pi/2 ; les trois autres l'atteignent.
# ---------------------------------------------------------------------------
assert solve(diff(x * (1 - x), x), x) == [Rational(1, 2)]
assert solve(diff(x * exp(-x), x), x) == [1] and limit(x * exp(-x), x, oo) == 0
assert sin(pi / 2) == 1
assert limit(atan(x), x, oo) == pi / 2 and all(math.atan(t) < math.pi / 2 for t in (1, 10, 1e6))

# ---------------------------------------------------------------------------
# ar-l1-continuite-q8 (vrai) : uniformément continue sur ]0,1[ => bornée : |f(x)| <= |f(1/2)| + N avec N eta >= 1/2.
# ---------------------------------------------------------------------------
eta = Rational(1, 7)
N = math.ceil(Rational(1, 2) / eta)
assert N * eta >= Rational(1, 2)

# ---------------------------------------------------------------------------
# ar-l1-continuite-image-intervalle-ex5 : cordes horizontales.
# ---------------------------------------------------------------------------
# question 2 : la somme télescopique des g(k/n) vaut f(1) - f(0) ; essai avec une fonction quelconque
F = lambda t: math.sin(7 * t) * (t - 1) * t + 3
for m in range(2, 9):
    s = sum(F((k + 1) / m) - F(k / m) for k in range(m))
    assert abs(s - (F(1) - F(0))) < 1e-12
# question 3 : f(x) = sin^2(5 pi x/2) - x : f(0) = f(1) = 0 et f(x + 2/5) - f(x) = -2/5 pour tout x
f3 = sin(5 * pi * x / 2) ** 2 - x
assert f3.subs(x, 0) == 0 and simplify(f3.subs(x, 1)) == 0
assert simplify(f3.subs(x, x + Rational(2, 5)) - f3) == Rational(-2, 5)
# forme générale : sin^2(pi x/a) - x sin^2(pi/a) ; f(x+a) - f(x) = -a sin^2(pi/a) (ici a = 2/5)
a = symbols("a", positive=True)
fa = sin(pi * x / a) ** 2 - x * sin(pi / a) ** 2
assert simplify(fa.subs(x, x + a) - fa + a * sin(pi / a) ** 2) == 0

# ---------------------------------------------------------------------------
# ar-l1-bornes-atteintes-ex4 : exemples 1/(1+x) (inf 0 non atteint) et -1/(1+x) (sup 0 non atteint).
# ---------------------------------------------------------------------------
assert limit(1 / (1 + x), x, oo) == 0 and limit(-1 / (1 + x), x, oo) == 0
assert all(1 / (1 + t) > 0 for t in range(0, 1000))

# ---------------------------------------------------------------------------
# ar-l1-continuite-heine-ex5 et variantes : plus grand eta = eps^2 pour sqrt sur [0, +oo[.
# ---------------------------------------------------------------------------
p, q = symbols("p q", nonnegative=True)
# |sqrt x - sqrt y| <= sqrt|x - y| : pour x = (p+q)^2 >= y = q^2, sqrt x - sqrt y = p et x - y = p^2 + 2pq >= p^2
assert simplify(((p + q) ** 2 - q**2) - p**2 - 2 * p * q) == 0
for eps, eta_ok, leurres in [
    (Rational(1, 10), Rational(1, 100), (Rational(1, 10), Rational(1, 20), Rational(1, 50))),  # exercice
    (Rational(1, 5), Rational(1, 25), (Rational(1, 5), Rational(2, 25), Rational(1, 10))),  # v1
    (Rational(3, 10), Rational(9, 100), (Rational(3, 10), Rational(3, 20), Rational(1, 10))),  # v2
]:
    assert eta_ok == eps**2 and sqrt(eta_ok) == eps  # eta = eps^2 convient (x = eta, y = 0 : égalité)
    for l in leurres:
        assert l > eps**2 and sqrt(l) > eps  # x = l, y = 0 met en défaut : leurre faux

print("enrichi_ar_l1_continuite : toutes les vérifications passent.")
