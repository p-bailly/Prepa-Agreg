"""Vérifications calculatoires du chapitre « ar-t-continuite » (analyse réelle, terminale) : cartes et
exercices de content/*/analyse-reelle/terminale-ar-t-continuite.json (les cartes ar-t-tvi et ar-t-point-fixe
sont vérifiées dans analyse_reelle.py).

Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise : v1, v2).
Lancer : python3 verification/ar_t_continuite.py (doit finir sans erreur).
"""

import math
from fractions import Fraction as F

from sympy import (E, Rational, cos, diff, exp, floor, limit, log, nsolve, oo, pi, real_roots, simplify, sin, solve,
                   sqrt, symbols)

x = symbols("x", real=True)
n = symbols("n", positive=True)

# ---------------------------------------------------------------------------
# ar-t-continuite-bijection (carte)
# ---------------------------------------------------------------------------
# Exemple : x^3 + x - 1, f' = 3x^2 + 1 > 0 ; f(0,6) = -0,184 ; f(0,7) = 0,043 ; racine unique dans ]0,6;0,7[.
f = x**3 + x - 1
assert diff(f, x) == 3 * x**2 + 1
assert f.subs(x, Rational(6, 10)) == Rational(-184, 1000) and f.subs(x, Rational(7, 10)) == Rational(43, 1000)
r = real_roots(f)
assert len(r) == 1 and 0.6 < float(r[0]) < 0.7
# Contre-exemple : x^3 - x sur [-2,2] : f(-2) = -6, f(2) = 6, racines -1, 0, 1.
assert (x**3 - x).subs(x, -2) == -6 and (x**3 - x).subs(x, 2) == 6 and sorted(solve(x**3 - x, x)) == [-1, 0, 1]
# Contre-exemple (continuité) : f(x) = x si x < 1, x + 1 si x >= 1 sur [0,2] : f(x) = 3/2 sans solution.
g = lambda t: t if t < 1 else t + 1
assert g(0) == 0 and g(2) == 3
assert not (F(3, 2) < 1) and not (F(1, 2) >= 1)  # x = 3/2 n'est pas < 1 ; x + 1 = 3/2 donne x = 1/2, pas >= 1
pts = [F(k, 1000) for k in range(0, 2001)]
assert all(g(a) < g(b) for a, b in zip(pts, pts[1:]))  # strictement croissante (contrôle sur une grille)
# Contre-exemple (monotonie large) : min(x, 1) = 1 sur tout [1, 2].
assert all(min(t, 1) == 1 for t in pts if t >= 1)
# Affirmation 1 : f(x) = x - 1 (x < 1), x (x >= 1) sur [0,2] : strictement croissante, f(0) = -1, f(2) = 2, sans zéro.
h = lambda t: t - 1 if t < 1 else t
assert h(0) == -1 and h(2) == 2 and all(h(t) != 0 for t in pts) and all(h(a) < h(b) for a, b in zip(pts, pts[1:]))
assert all(h(t) < 0 for t in pts if t < 1) and all(h(t) >= 1 for t in pts if t >= 1)
# Affirmation 3 : 1/(1-x) croissante sur ]0,1[, limite +oo en 1.
assert diff(1 / (1 - x), x) == 1 / (1 - x) ** 2 and limit(1 / (1 - x), x, 1, "-") == oo
assert limit(1 / (1 - x), x, 0, "+") == 1
# Contre-exemple : exp a pour limites 0 et +oo, et ne s'annule pas.
assert limit(exp(x), x, -oo) == 0 and limit(exp(x), x, oo) == oo and solve(exp(x), x) == []
# Trous : leurre « positive » (trou 1) : x + 1 (x < 1), x + 2 (x >= 1) sur [0,2] ; k = 5/2 sans solution.
p = lambda t: t + 1 if t < 1 else t + 2
assert p(0) == 1 and p(2) == 4 and all(p(t) != F(5, 2) for t in pts) and all(p(t) > 0 for t in pts)
# Trou 2, leurres « positive » / « non constante » : x^2 sur [-1,1], continue, f(-1) = f(1) = 1, deux solutions.
assert sorted(solve(x**2 - 1, x)) == [-1, 1]

# ---------------------------------------------------------------------------
# ar-t-continuite-bijection-ex1 (+ v1, v2) : f croît de -2 (en -3) à 5 (en 1), décroît jusqu'à 1 (en 4).
# Nombre de solutions de f(x) = k : sur [-3,1] si -2 <= k <= 5, sur [1,4] si 1 <= k <= 5, en retirant le double
# compte de x = 1 (seulement si k = 5).
# ---------------------------------------------------------------------------
def nb_solutions(k):
    gauche = -2 <= k <= 5
    droite = 1 <= k <= 5
    return gauche + droite - (k == 5)


assert nb_solutions(3) == 2  # ex1 : 2 (leurres 1, 0, 3)
assert nb_solutions(0) == 1  # v1 : 1 (leurres 2, 0, 3)
assert nb_solutions(1) == 2  # v2 : 2 (dont x = 4 ; leurres 1, 3, 0)
# Contrôle sur une fonction affine par morceaux qui réalise ce tableau.
fp = lambda t: -2 + F(7, 4) * (t + 3) if t <= 1 else 5 - F(4, 3) * (t - 1)
assert fp(-3) == -2 and fp(1) == 5 and fp(4) == 1
grille = [F(k, 600) for k in range(-1800, 2401)]
for k, attendu in ((3, 2), (0, 1), (1, 2)):
    # solutions exactes : antécédents sur chaque morceau, dans le bon intervalle
    sol = set()
    t1 = F(k + 2, 1) / F(7, 4) - 3
    if -3 <= t1 <= 1:
        sol.add(t1)
    t2 = 1 + (5 - F(k)) / F(4, 3)
    if 1 <= t2 <= 4:
        sol.add(t2)
    assert len(sol) == attendu and all(fp(s) == k for s in sol)

# ---------------------------------------------------------------------------
# ar-t-continuite-bijection-ex2 (+ v1, v2) : encadrements au dixième.
# ---------------------------------------------------------------------------
for pol, d, (a, b), vals, leurres in (
    (x**3 + 3 * x - 5, 3 * x**2 + 3, (Rational(11, 10), Rational(12, 10)), (-0.369, 0.328),
     [(1, 1.1), (1.2, 1.3), (1.5, 1.6)]),
    (x**3 + x - 3, 3 * x**2 + 1, (Rational(12, 10), Rational(13, 10)), (-0.072, 0.497),
     [(1.1, 1.2), (1.3, 1.4), (1.5, 1.6)]),
    (x**3 + 4 * x + 1, 3 * x**2 + 4, (Rational(-3, 10), Rational(-2, 10)), (-0.227, 0.192),
     [(-0.2, -0.1), (-0.4, -0.3), (0.2, 0.3)]),
):
    assert simplify(diff(pol, x) - d) == 0
    r = real_roots(pol)
    assert len(r) == 1
    al = float(r[0])
    assert float(a) < al < float(b)
    assert math.isclose(float(pol.subs(x, a)), vals[0], abs_tol=1e-12)
    assert math.isclose(float(pol.subs(x, b)), vals[1], abs_tol=1e-12)
    for lo, hi in leurres:
        assert not (lo < al < hi)
# Valeurs f(1), f(2) (ex2, v1) et f(-1), f(0) (v2) citées.
assert (x**3 + 3 * x - 5).subs(x, 1) == -1 and (x**3 + 3 * x - 5).subs(x, 2) == 9
assert (x**3 + x - 3).subs(x, 1) == -1 and (x**3 + x - 3).subs(x, 2) == 7
assert (x**3 + 4 * x + 1).subs(x, -1) == -4 and (x**3 + 4 * x + 1).subs(x, 0) == 1

# ---------------------------------------------------------------------------
# ar-t-continuite-bijection-ex3 (+ v1, v2) : nombre de racines réelles et de racines > 0.
# ---------------------------------------------------------------------------
for pol, total, positives, valeurs in (
    (x**3 - 3 * x + 1, 3, 2, {-1: 3, 1: -1, 0: 1}),
    (x**3 - 3 * x - 1, 3, 1, {-1: 1, 1: -3, 0: -1}),
    (x**3 - 12 * x - 20, 1, 1, {-2: -4, 2: -36, 4: -4, 5: 45}),
):
    r = [float(t) for t in real_roots(pol)]
    assert len(set(round(t, 9) for t in r)) == total and sum(t > 0 for t in r) == positives
    for t, vt in valeurs.items():
        assert pol.subs(x, t) == vt
# Les leurres (« 3 dont 1 », « 2 dont 1 », « 1 positive », …) diffèrent de (total, positives) ci-dessus ;
# x^3 - 3x + 1 : racine du milieu dans ]0,1[ ; x^3 - 3x - 1 : dans ]-1,0[ ; x^3 - 12x - 20 : dans ]4,5[.
assert any(0 < float(t) < 1 for t in real_roots(x**3 - 3 * x + 1))
assert any(-1 < float(t) < 0 for t in real_roots(x**3 - 3 * x - 1))
assert 4 < float(real_roots(x**3 - 12 * x - 20)[0]) < 5

# ---------------------------------------------------------------------------
# ar-t-continuite-bijection-ex4 (+ v1, v2) : x^n + c x - d = 0, unique racine u_n dans ]0,1[, limite.
# ---------------------------------------------------------------------------
def racine(c, d, k):
    """Racine dans ]0,1[ de x^k + c x - d (dichotomie en flottants, 200 itérations)."""
    lo, hi = 0.0, 1.0
    for _ in range(200):
        m = (lo + hi) / 2
        if m**k + c * m - d < 0:
            lo = m
        else:
            hi = m
    return (lo + hi) / 2


for c, d, u1, u2, lim, leurres in (
    (1, 1, F(1, 2), (sqrt(5) - 1) / 2, 1, [Rational(1, 2), (sqrt(5) - 1) / 2, 0]),
    (2, 1, F(1, 3), sqrt(2) - 1, Rational(1, 2), [1, Rational(1, 3), sqrt(2) - 1]),
    (3, 2, F(1, 2), (sqrt(17) - 3) / 2, Rational(2, 3), [1, Rational(1, 2), (sqrt(17) - 3) / 2]),
):
    assert (0 + c * 0 - d) < 0 < (1 + c - d)  # f_n(0) < 0 < f_n(1)
    assert u1 ** 1 + c * u1 - d == 0  # u_1
    assert simplify(u2**2 + c * u2 - d) == 0 and 0 < float(u2) < 1  # u_2
    us = [racine(c, d, k) for k in range(1, 400)]
    assert math.isclose(us[0], float(u1), abs_tol=1e-12) and math.isclose(us[1], float(u2), abs_tol=1e-12)
    assert all(a < b for a, b in zip(us[:40], us[1:40]))  # croissante (strict, tant que le flottant distingue)
    assert all(a <= b for a, b in zip(us, us[1:]))
    # f_{n+1}(u_n) = u_n^{n+1} - u_n^n < 0
    assert all(us[k - 1] ** (k + 1) + c * us[k - 1] - d < 0 for k in range(1, 40))
    if lim != 1:
        assert all(t < float(lim) for t in us[:40]) and all(t <= float(lim) for t in us)  # u_n < d/c
        assert all(F(d, c) ** k + c * F(d, c) - d == F(d, c) ** k > 0 for k in range(1, 50))  # f_n(d/c) = (d/c)^n
        assert Rational(d, c) == lim
    assert abs(us[-1] - float(lim)) < 0.02 and all(abs(us[-1] - float(le)) > 0.04 for le in leurres)
# Limite exacte en v1, v2 : u_n = (d - u_n^n)/c avec u_n^n <= (d/c)^n -> 0.
assert limit(Rational(1, 2) ** n, n, oo) == 0 and limit(Rational(2, 3) ** n, n, oo) == 0

# ---------------------------------------------------------------------------
# ar-t-continuite-image-suite (carte)
# ---------------------------------------------------------------------------
# Exemple : (4n+1)/n -> 4, sqrt -> 2.
assert limit((4 * n + 1) / n, n, oo) == 4 and limit(sqrt((4 * n + 1) / n), n, oo) == 2
# Contre-exemple : E(1 - 1/(n+2)) = 0, E(1) = 1 ; u_n dans ]0,1[.
assert all(0 < 1 - F(1, k + 2) < 1 and math.floor(1 - F(1, k + 2)) == 0 for k in range(200)) and math.floor(1) == 1
assert limit(1 - 1 / (n + 2), n, oo) == 1
# Contre-exemple : 1/x sur ]0,+oo[, u_n = 1/(n+1) -> 0, f(u_n) = n+1 -> +oo.
assert limit(1 / (n + 1), n, oo) == 0 and simplify(1 / (1 / (n + 1)) - (n + 1)) == 0
# Affirmation 1 : E(1 + 1/(n+2)) = 1 = E(1), mais E(1 - 1/(n+2)) = 0.
assert all(math.floor(1 + F(1, k + 2)) == 1 for k in range(200))
# Affirmation 2 : x^2 et (-1)^n.
assert all(((-1) ** k) ** 2 == 1 for k in range(20))
# Affirmation 3 : x/(1+x^2) en 2 vaut 2/5.
assert (x / (1 + x**2)).subs(x, 2) == Rational(2, 5)
# Trou 2, leurres : f(u_n) -> f(l) ≠ l en général (ex. f = exp, l = 0 : f(l) = 1).
assert exp(0) != 0

# ---------------------------------------------------------------------------
# ar-t-continuite-image-suite-ex1..ex3 (+ v1, v2)
# ---------------------------------------------------------------------------
def lim_et_leurres(expr, bonne, leurres):
    L = limit(expr, n, oo)
    assert (L == bonne) if bonne in (oo, -oo) else simplify(L - bonne) == 0
    for le in leurres:
        assert le == oo and L != oo or le != oo and (L in (oo, -oo) or abs(float(L) - float(le)) > 1e-9)


# ex1 : u_n = (2n+1)/(n+3) -> 2 ; sqrt(u_n^2 + 5) -> 3 ; leurres sqrt5, sqrt7 (= sqrt(2+5)), +oo.
u = (2 * n + 1) / (n + 3)
assert limit(u, n, oo) == 2
lim_et_leurres(sqrt(u**2 + 5), 3, [sqrt(5), sqrt(7), oo])
u = (3 * n - 1) / (n + 2)  # v1 : -> 3 ; sqrt(u^2+16) -> 5 ; leurres 4, sqrt19, 3
assert limit(u, n, oo) == 3
lim_et_leurres(sqrt(u**2 + 16), 5, [4, sqrt(19), 3])
u = (n + 4) / (2 * n)  # v2 : -> 1/2 ; sqrt(8u^2+7) -> 3 ; leurres sqrt7, sqrt11 (= sqrt(8/2+7)), 1/2
assert limit(u, n, oo) == Rational(1, 2) and simplify(u - (Rational(1, 2) + 2 / n)) == 0
lim_et_leurres(sqrt(8 * u**2 + 7), 3, [sqrt(7), sqrt(11), Rational(1, 2)])
# ex2 (+ v1, v2) : cos / sin d'une suite convergente.
lim_et_leurres(cos(pi * n / (3 * n + 1)), Rational(1, 2), [sqrt(3) / 2, 1])
lim_et_leurres(sin(pi * n**2 / (4 * n**2 + 1)), sqrt(2) / 2, [0, Rational(1, 2)])
lim_et_leurres(cos((2 * pi * n + 1) / (3 * n)), Rational(-1, 2), [Rational(1, 2), -sqrt(3) / 2])
assert limit(pi * n / (3 * n + 1), n, oo) == pi / 3 and limit((2 * pi * n + 1) / (3 * n), n, oo) == 2 * pi / 3
# ex3 (+ v1, v2) : forme indéterminée puis continuité.
assert limit(sqrt(n**2 + 2 * n) - n, n, oo) == 1
lim_et_leurres(exp(sqrt(n**2 + 2 * n) - n), E, [1, oo, E**2])
assert limit(sqrt(n**2 + 6 * n) - n, n, oo) == 3
lim_et_leurres(log(sqrt(n**2 + 6 * n) - n + 1), log(4), [log(7), 0, oo])
assert limit(n * (sqrt(n**2 + 1) - n), n, oo) == Rational(1, 2)
lim_et_leurres(sqrt(n * (sqrt(n**2 + 1) - n)), sqrt(2) / 2, [0, Rational(1, 2), 1])

# ---------------------------------------------------------------------------
# ar-t-continuite-image-suite-ex4 : f(2x) = f(x) ; contre-exemple sans continuité : f(0) = 0, f(x) = 1 sinon.
# ---------------------------------------------------------------------------
fc = lambda t: 0 if t == 0 else 1
assert all(fc(2 * t) == fc(t) for t in [F(k, 7) for k in range(-50, 51)])
assert all(fc(F(1, 2**k)) == 1 for k in range(60)) and fc(0) == 0

print("ar_t_continuite : toutes les vérifications passent.")
