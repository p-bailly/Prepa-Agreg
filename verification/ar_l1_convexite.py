"""Vérifications calculatoires du chapitre « ar-l1-convexite » (Jensen, croissance des pentes, convexité des fonctions
dérivables, inégalité arithmético-géométrique) : cartes et exercices de content/*/analyse-reelle/l1-ar-l1-convexite.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice concerné (v1, v2 = variantes 1 et 2).
Lancer : python3 verification/ar_l1_convexite.py (doit finir sans erreur).
"""

import itertools
import random

from sympy import (Max, Rational, cos, diff, exp, expand, factor, limit, log, nsimplify, oo, pi, root, simplify, sin,
                   sqrt, symbols, solve, E)

x, y, t, a, h = symbols("x y t a h", real=True)
random.seed(2026)


def egal(u, v):
    return simplify(u - v) == 0


def pente(f, u, v):
    return (f(v) - f(u)) / (v - u)


# ===========================================================================
# Carte ar-l1-convexite-jensen
# ===========================================================================
# Contre-exemples : poids (2, -1) avec x^2 ; poids (1/2, 1/4) avec f = 1 ; cos(pi x) convexe près de ±1.
sq = lambda u: u**2
assert sq(2 * 1 - 1 * 0) == 4 and 2 * sq(1) - sq(0) == 2 and 4 > 2
assert 1 > Rational(1, 2) + Rational(1, 4)
assert simplify(diff(cos(pi * x), x, 2).subs(x, 1)) == pi**2 and simplify(diff(cos(pi * x), x, 2).subs(x, -1)) == pi**2
assert cos(0) == 1 and (cos(-pi) + cos(pi)) / 2 == -1
# Affirmations : x * x^2 = x^3 non convexe (-1 > -4) ; -(x^2) non convexe (0 > -1).
assert (-1)**3 == -1 and ((-2)**3 + 0) / 2 == -4
assert -(0**2) == 0 and (-(1) - 1) / 2 == -1
# Trou 1, leurres : définition « ≥ » (concavité), « ≥ min » ; avec √ (qui vérifie ces deux propriétés), Jensen est faux.
assert sqrt(Rational(1, 2) * 0 + Rational(1, 2) * 4) > (sqrt(0) + sqrt(4)) / 2      # √2 > 1
assert sqrt(2) <= Max(sqrt(0), sqrt(4))                                              # √ vérifie le leurre « ≤ max »
# Trou 2 et 3, leurres : poids (2,-1) (non nul, rationnel) ; f(Σλx) ≤ Σ f(λx) faux pour x^2, λ = (1/2,1/2), x = (1,1).
assert sq(1) > sq(Rational(1, 2)) + sq(Rational(1, 2))

# ar-l1-convexite-jensen-ex1 et variantes : D = (1-t)f(x) + t f(y) - f((1-t)x + ty).
def D(f):
    return expand((1 - t) * f(x) + t * f(y) - f((1 - t) * x + t * y))


assert egal(D(lambda u: u**2), t * (1 - t) * (x - y)**2)                          # ex1
assert egal(D(lambda u: 3 * u**2 - 4 * u + 1), 3 * t * (1 - t) * (x - y)**2)      # ex1 v1
assert egal(D(lambda u: -u**2 + 2 * u), -t * (1 - t) * (x - y)**2)                # ex1 v2
assert D(lambda u: -u**2 + 2 * u).subs({x: 0, y: 2, t: Rational(1, 2)}) == -1     # ex1 v2 : non convexe
# Leurres : expressions différentes de D.
assert not egal(D(lambda u: u**2), t * (1 - t) * (x + y)**2) and not egal(D(lambda u: u**2), t**2 * (x - y)**2)
assert not egal(D(lambda u: 3 * u**2 - 4 * u + 1), t * (1 - t) * (x - y)**2)
assert not egal(D(lambda u: 3 * u**2 - 4 * u + 1), 3 * t * (1 - t) * (x - y)**2 - 4 * t * (1 - t) * (x - y))
assert not egal(D(lambda u: -u**2 + 2 * u), -t * (1 - t) * (x - y)**2 + 2 * t * (1 - t) * (x - y))


# ar-l1-convexite-jensen-ex2 et variantes : majoration par la corde, fonction g = max de deux affines qui l'atteint.
def verifie_corde(points, x0, s_court, s_long, M, borne_longue, g_aff):
    (x1, f1), (x2, f2), (x3, f3) = points            # x1 < x2 < x0 < x3
    tc = Rational(x0 - x2, x3 - x2)
    tl = Rational(x0 - x1, x3 - x1)
    assert tc == s_court and tl == s_long
    assert (1 - tc) * f2 + tc * f3 == M
    assert (1 - tl) * f1 + tl * f3 == borne_longue and M < borne_longue
    g = lambda u: Max(*[p * u + q for (p, q) in g_aff])
    assert all(g(xi) == fi for xi, fi in points) and g(x0) == M
    # g est convexe (max de deux affines) ; pentes croissantes de l'interpolation
    assert pente(g, x1, x2) <= pente(g, x2, x3)


verifie_corde([(0, 2), (2, -1), (6, 5)], 3, Rational(1, 4), Rational(1, 2), Rational(1, 2), Rational(7, 2),
              [(Rational(-3, 2), 2), (Rational(3, 2), -4)])                              # ex2 : M = 1/2
verifie_corde([(-1, 4), (1, 0), (4, 3)], 2, Rational(1, 3), Rational(3, 5), 1, Rational(17, 5),
              [(-2, 2), (1, -1)])                                                        # ex2 v1 : M = 1
verifie_corde([(0, 0), (1, -2), (5, 6)], 2, Rational(1, 4), Rational(2, 5), 0, Rational(12, 5),
              [(-2, 0), (2, -4)])                                                        # ex2 v2 : M = 0
# Leurres distincts de M : 7/2, 2, -5/2 (ex2) ; 17/5, 3/2, -2 (v1) ; 12/5, 2, -4 (v2).
assert Rational(1, 2) not in (Rational(7, 2), 2, Rational(-5, 2))
assert 1 not in (Rational(17, 5), Rational(3, 2), -2) and 0 not in (Rational(12, 5), 2, -4)

# ar-l1-convexite-jensen-ex3 et variantes : minimums 1/n, n, n^2 (contrôle exact pour n = 1..6 au point d'égalité,
# et échantillonnage aléatoire de l'inégalité).
for n in range(1, 7):
    assert n * Rational(1, n)**2 == Rational(1, n)                  # ex3 : au point x_i = 1/n
    assert n * 1**3 == n                                            # ex3 v1 : x_i = 1
    assert n * n == n**2                                            # ex3 v2 : x_i = 1/n, Σ 1/x_i = n·n
    for _ in range(300):
        w = [random.random() for _ in range(n)]
        s = sum(w)
        z = [wi / s for wi in w]                                    # somme 1, positifs
        assert sum(zi**2 for zi in z) >= 1 / n - 1e-12              # ex3
        assert sum((n * zi)**3 for zi in z) >= n - 1e-9             # ex3 v1 (somme n)
        assert sum(1 / zi for zi in z) >= n**2 - 1e-9               # ex3 v2
        r = [random.uniform(-3, 3) for _ in range(n)]               # ex3 : réels de signes quelconques
        sr = sum(r)
        if abs(sr) > 1e-3:
            zr = [ri / sr for ri in r]
            assert sum(zi**2 for zi in zr) >= 1 / n - 1e-9
# Leurres : 1, 1/n^2, 0 ; 1, n^3, n^2 ; n, 1/n, 1 diffèrent du minimum pour n = 2.
assert Rational(1, 2) not in (1, Rational(1, 4), 0) and 2 not in (1, 8, 4) and 4 not in (2, Rational(1, 2), 1)

# ar-l1-convexite-jensen-ex4 et variantes : maximums 3√3/2, 3√3/8, 1/8 au triangle équilatéral.
assert simplify(diff(sin(x), x, 2) + sin(x)) == 0                                    # sin'' = -sin ≤ 0 sur [0, π]
assert simplify(diff(log(sin(x)), x, 2) + 1 / sin(x)**2) == 0                        # (ln sin)'' = -1/sin² < 0
assert egal(3 * sin(pi / 3), 3 * sqrt(3) / 2)                                        # ex4
assert egal(sin(pi / 3)**3, 3 * sqrt(3) / 8)                                         # ex4 v1
assert egal(sin(pi / 6)**3, Rational(1, 8))                                          # ex4 v2
import math
for _ in range(20000):
    A_, B_ = random.uniform(0, math.pi), random.uniform(0, math.pi)
    C_ = math.pi - A_ - B_
    if C_ <= 0:
        continue
    assert math.sin(A_) + math.sin(B_) + math.sin(C_) <= 3 * math.sqrt(3) / 2 + 1e-12
    assert math.sin(A_) * math.sin(B_) * math.sin(C_) <= 3 * math.sqrt(3) / 8 + 1e-12
    assert math.sin(A_ / 2) * math.sin(B_ / 2) * math.sin(C_ / 2) <= 1 / 8 + 1e-12
# Leurres : 3, 3/2, 3√3/4 ; 1, 3√3/2, 1/8 ; 3√3/8, 1/2, 1/27 sont différents des maximums.
assert all(not egal(v, 3 * sqrt(3) / 2) for v in (3, Rational(3, 2), 3 * sqrt(3) / 4))
assert all(not egal(v, 3 * sqrt(3) / 8) for v in (1, 3 * sqrt(3) / 2, Rational(1, 8)))
assert all(not egal(v, Rational(1, 8)) for v in (3 * sqrt(3) / 8, Rational(1, 2), Rational(1, 27)))

# ar-l1-convexite-jensen-ex5 : t_k = floor(2^k t)/2^k est dyadique dans [0, 1] et 0 ≤ t - t_k < 2^-k.
for tt in (Rational(1, 3), Rational(5, 7), Rational(1)):
    for k in range(0, 12):
        tk = Rational(int(2**k * tt), 2**k)
        assert 0 <= tk <= 1 and 0 <= tt - tk < Rational(1, 2**k)

# ===========================================================================
# Carte ar-l1-convexite-pentes
# ===========================================================================
# Exemple : τ_a(x) = x + a pour x^2 ; pentes de exp en (0, 1, 2).
assert egal((x**2 - a**2) / (x - a), x + a)
p01, p02, p12 = float(E - 1), float((E**2 - 1) / 2), float(E**2 - E)
assert abs(p01 - 1.72) < 0.005 and abs(p02 - 3.19) < 0.005 and abs(p12 - 4.67) < 0.005 and p01 <= p02 <= p12
# Contre-exemples : √ en (1, 4, 9) ; -|t| en (-1, 0, 1) ; t^2 en (0, 2, 1).
assert pente(sqrt, 1, 4) == Rational(1, 3) and pente(sqrt, 1, 9) == Rational(1, 4)
mabs = lambda u: -abs(u)
assert (pente(mabs, -1, 0), pente(mabs, -1, 1), pente(mabs, 0, 1)) == (1, 0, -1)
assert pente(sq, 0, 2) == 2 and pente(sq, 0, 1) == 1
# Trous : leurres faux pour x^2 en (0, 1, 2) : pentes 1, 2, 3 ; τ positive faux pour f = -x ; stricte faux pour affine.
assert (pente(sq, 0, 1), pente(sq, 0, 2), pente(sq, 1, 2)) == (1, 2, 3)
assert not (2 <= 1) and not (3 <= 2)
assert pente(lambda u: -u, 0, 1) == -1
# Affirmation : f = 0 sur [0,1[, f(1) = 1 est convexe (vérification exhaustive sur une grille rationnelle).
f_saut = lambda u: 1 if u == 1 else 0
grille = [Rational(i, 8) for i in range(9)]
for u, v, s in itertools.product(grille, grille, grille):
    assert f_saut((1 - s) * u + s * v) <= (1 - s) * f_saut(u) + s * f_saut(v)

# ar-l1-convexite-pentes-ex1 et variantes : m = f(x2) + p(x1,x2)(x3 - x2), atteint par la fonction affine.
for (x1, f1), (x2, f2), x3, m, leurres in [((0, 1), (2, 3), 5, 6, (3, 4, Rational(15, 2))),          # ex1
                                          ((1, 2), (3, -2), 6, -8, (4, -2, -14)),                   # ex1 v1
                                          ((-2, 5), (0, 1), 3, -5, (7, 1, -3))]:                    # ex1 v2
    p = Rational(f2 - f1, x2 - x1)
    assert f2 + p * (x3 - x2) == m
    aff = lambda u: f1 + p * (u - x1)
    assert aff(x2) == f2 and aff(x3) == m
    assert m not in leurres

# ar-l1-convexite-pentes-ex2 et variantes : dérivées à gauche et à droite de max(x^2, affine).
for aff, x0, (dg, dd), (r1, r2) in [(x + 2, 2, (1, 4), (-1, 2)),          # ex2
                                    (2 * x + 3, -1, (-2, 2), (-1, 3)),    # ex2 v1
                                    (3 * x - 2, 1, (2, 3), (1, 2))]:      # ex2 v2
    assert set(solve(x**2 - aff, x)) == {r1, r2}
    # près de x0 : à gauche le morceau le plus grand en x0 - 1/100, à droite celui en x0 + 1/100
    gauche = x**2 if (x**2 - aff).subs(x, x0 - Rational(1, 100)) > 0 else aff
    droite = x**2 if (x**2 - aff).subs(x, x0 + Rational(1, 100)) > 0 else aff
    f0 = Max(x**2, aff).subs(x, x0)
    assert gauche.subs(x, x0) == f0 and droite.subs(x, x0) == f0
    assert limit((gauche - f0) / (x - x0), x, x0, "-") == dg
    assert limit((droite - f0) / (x - x0), x, x0, "+") == dd
    assert dg < dd
assert (1, 4) != (4, 1) and (1, 4) != (4, 4) and (1, 4) != (1, 2)

# ar-l1-convexite-pentes-ex3 : e^{-x} convexe, majorée par 1 sur ]0, +∞[ ; x^2 convexe minorée.
assert diff(exp(-x), x, 2) == exp(-x)

# ar-l1-convexite-pentes-ex4 : τ_1(x) = 1/(1 - x) → +∞ pour la fonction à saut.
assert egal((0 - 1) / (x - 1), 1 / (1 - x)) and limit(1 / (1 - x), x, 1, "-") == oo

# ===========================================================================
# Carte ar-l1-convexite-derivables
# ===========================================================================
# Exemple : x ln x, f' = ln x + 1, tangente en 1 : y = x - 1 ; f = x^2/2 (x ≥ 0), 0 sinon, dérivable en 0.
assert egal(diff(x * log(x), x), log(x) + 1)
assert limit((x**2 / 2) / x, x, 0, "+") == 0
# Contre-exemples : ln 2 > (ln 1 + ln 3)/2 ; x^3 en (-2, 0) ; fonction nulle sur [0,1[ valant -1 en 1.
assert log(2).evalf() > (log(3) / 2).evalf() and abs(float(log(2)) - 0.693) < 5e-4 and abs(float(log(3) / 2) - 0.549) < 5e-4
assert diff(x**3, x, 2).subs(x, 1) == 6 and (-1)**3 > ((-2)**3 + 0) / 2
assert 0 > (0 + (-1)) / Rational(2)
# Affirmations : x^4 convexe avec f''(0) = 0 ; x^2 avec f' croissante non monotone.
assert diff(x**4, x, 2) == 12 * x**2

# ar-l1-convexite-derivables-ex1 et variantes : signe de f''.
assert egal(diff(x**4 - 6 * x**2 + 1, x, 2), 12 * (x - 1) * (x + 1))                        # ex1
assert set(solve(diff(x**4 - 6 * x**2 + 1, x), x)) == {0, sqrt(3), -sqrt(3)}                    # ex1 (leurre ±√3)
assert egal(diff(x * exp(x), x, 2), (x + 2) * exp(x)) and egal(diff(x * exp(x), x), (x + 1) * exp(x))   # ex1 v1
assert egal(diff(1 / (1 + x**2), x, 2), (6 * x**2 - 2) / (1 + x**2)**3)                          # ex1 v2
assert set(solve(6 * x**2 - 2, x)) == {1 / sqrt(3), -1 / sqrt(3)}
assert diff(1 / (1 + x**2), x, 2).subs(x, Rational(1, 2)) < 0                                    # v2 : [0,+∞[ faux

# ar-l1-convexite-derivables-ex2 et variantes : tangentes et leurres faux.
assert egal(log(E) + (x - E) / E, x / E)                                                          # ex2
for u in (Rational(1, 2), 1, 2, 3, 10, 100):
    assert log(u) <= u / E
assert log(1) < 1 / E                       # leurre « ln x ≥ x/e » faux en 1
assert log(E) > E / E - 1                   # leurre « ln x ≤ x/e - 1 » faux en e
assert log(E) > E - E                       # leurre « ln x ≤ x - e » faux en e
assert egal(E + E * (x - 1), E * x)                                                               # ex2 v1
assert exp(0) > E * 0 and exp(1) < E * 1 + 1 and exp(1) < 1 + E       # leurres v1 faux (en 0, en 1, en 1)
assert egal(2 + (x - 4) / 4, 1 + x / 4) and diff(sqrt(x), x).subs(x, 4) == Rational(1, 4)          # ex2 v2
assert egal(1 + x / 4 - sqrt(x), (sqrt(x) - 2)**2 / 4)
assert sqrt(0) < 1 + 0 and sqrt(1) > Rational(1, 4)     # leurres v2 « ≥ » faux en 0, « ≤ x/4 » faux en 1
assert sqrt(16) > Rational(3, 2) + Rational(16, 8)      # leurre v2 « ≤ 3/2 + x/8 » faux en 16

# ar-l1-convexite-derivables-ex3 et variantes : point critique et minimum.
for f, x0, mini, f0, leurres in [
        (exp(x) + exp(-2 * x), log(2) / 3, 3 / root(4, 3), 2, (2, Rational(3, 2), root(2, 3))),            # ex3
        (exp(2 * x) + 4 * exp(-x), log(2) / 3, 3 * root(4, 3), 5, (5, 3 * root(2, 3), 4)),                 # ex3 v1
        (exp(3 * x) + 6 * exp(-x), log(2) / 4, 4 * root(8, 4), 7, (7, 8, 4 * root(2, 4)))]:               # ex3 v2
    fp = diff(f, x)
    assert egal(fp.subs(x, x0), 0)
    assert diff(fp, x).subs(x, x0).evalf() > 0 and simplify(diff(fp, x)).is_positive is not False
    assert egal(f.subs(x, x0), mini) and f.subs(x, 0) == f0 and mini.evalf() < f0
    assert all(not egal(l, mini) for l in leurres)
assert egal(exp(3 * x) + 6 * exp(-x), exp(3 * x) + 6 * exp(-x)) and egal(4 * root(8, 4), 2**Rational(11, 4))

# ar-l1-convexite-derivables-ex4 : f'' = α(α-1)(1+x)^(α-2) ; contre-exemple α = 1/2, x = 1 ; α = -1 et α = 2 vérifient.
al = symbols("alpha", real=True)
assert egal(diff((1 + h)**al, h, 2), al * (al - 1) * (1 + h)**(al - 2))
assert sqrt(2) < Rational(3, 2)
for alpha in (Rational(-1), Rational(-1, 2), Rational(0), Rational(1), Rational(2), Rational(5, 2)):
    for xv in (Rational(-9, 10), Rational(-1, 2), Rational(1, 3), Rational(2), Rational(10)):
        assert ((1 + xv)**alpha - 1 - alpha * xv).evalf() >= -1e-12
for alpha in (Rational(1, 3), Rational(1, 2), Rational(9, 10)):
    for xv in (Rational(-9, 10), Rational(-1, 2), Rational(1, 3), Rational(2), Rational(10)):
        assert ((1 + xv)**alpha - 1 - alpha * xv).evalf() < 0

# ===========================================================================
# Carte ar-l1-convexite-iag
# ===========================================================================
# Contre-exemples : (-2, -8) ; (-1, -1, 4) ; (0, 1) ; (1/2, 3/2).
assert sqrt((-2) * (-8)) == 4 and Rational(-2 - 8, 2) == -5
assert (-1) * (-1) * 4 == 4 and root(4, 3).evalf() > Rational(2, 3) and abs(float(root(4, 3)) - 1.59) < 5e-3
assert sqrt(0 * 1) < Rational(1, 2) and sqrt(Rational(3, 4)) < 1 and abs(float(sqrt(Rational(3, 4))) - 0.87) < 5e-3
# Affirmations : ((x+y)/2)^2 - xy = ((x-y)/2)^2 ; (-1, -1, 2) : xyz = 2 > 0.
assert egal(((x + y) / 2)**2 - x * y, ((x - y) / 2)**2)
assert (-1) * (-1) * 2 == 2 and ((-1 - 1 + 2) / 3)**3 == 0
# Trous : leurres faux pour (0, 1), (2, 2), (1/2, 1/2) ; cas d'égalité « sommes = n » faux pour (1/2, 3/2) ;
# « tous égaux à 1 » faux pour (2, 2) qui réalise l'égalité.
assert Rational(1, 2) > sqrt(0) and 2 * 2 > Rational(2 + 2, 2) and sqrt(Rational(1, 4)) > Rational(1, 4)
assert sqrt(2 * 2) == Rational(2 + 2, 2)
# Démonstration : somme des ln(x_i/m) ≤ somme (x_i/m - 1) = 0.
g = x - 1 - log(x)
assert g.subs(x, 1) == 0 and egal(diff(g, x), 1 - 1 / x)

# ar-l1-convexite-iag-ex1 et variantes : minimum de αx + β/x = 2√(αβ) en x = √(β/α).
for al_, be_, mini, x0, leurres in [(1, 4, 4, 2, [(2, 1), (4, 4), (5, 1)]),                                     # ex1
                                    (2, 9, 6 * sqrt(2), 3 * sqrt(2) / 2,
                                     [(6, Rational(3, 2)), (3 * sqrt(2), 3 * sqrt(2) / 2), (6 * sqrt(2), Rational(3, 2))]),  # v1
                                    (3, 12, 12, 2, [(6, 2), (12, 4), (15, 1)])]:                                 # ex1 v2
    f = al_ * x + be_ / x
    assert egal(f.subs(x, x0), mini) and egal(2 * sqrt(al_ * be_), mini)
    assert egal(diff(f, x).subs(x, x0), 0)
    for (v, pt) in leurres:      # chaque leurre est faux : mauvaise valeur, ou valeur non atteinte au point annoncé
        assert not (egal(v, mini) and egal(pt, x0))
        assert not egal(f.subs(x, pt), v) or not egal(v, mini)
assert (2 * Rational(3, 2) + 9 / Rational(3, 2)) == 9                       # v1 : f(3/2) = 9

# ar-l1-convexite-iag-ex2 et variantes : volumes maximaux 64, 216, 36.
assert 4**3 == 64 and 4 + 4 + 4 == 12                                                     # ex2
assert 6**2 * 6 == 216 and 2 * 6 + 6 == 18                                                # ex2 v1
assert 6 * 3 * 2 == 36 and 6 + 2 * 3 + 3 * 2 == 18 and Rational(6**3, 6) == 36            # ex2 v2
for _ in range(20000):
    u, v = random.uniform(0, 12), random.uniform(0, 12)
    w = 12 - u - v
    if w > 0:
        assert u * v * w <= 64 + 1e-9                                                     # ex2
    aa = random.uniform(0, 9)
    assert aa**2 * (18 - 2 * aa) <= 216 + 1e-9                                            # ex2 v1
    b_, c_ = random.uniform(0, 9), random.uniform(0, 6)
    a_ = 18 - 2 * b_ - 3 * c_
    if a_ > 0:
        assert a_ * b_ * c_ <= 36 + 1e-9                                                  # ex2 v2
assert 64 not in (48, 36, 216) and 216 not in (324, 108, 162) and 36 not in (216, 8, 24)

# ar-l1-convexite-iag-ex3 et variantes : minimums 8, 3, 9 en a = b = c ; échantillonnage.
assert (2 * 2 * 2) / 1 == 8 and 1 + 1 + 1 == 3 and 3 * 3 == 9
for _ in range(20000):
    p_, q_, r_ = (random.uniform(0.01, 10) for _ in range(3))
    assert (p_ + q_) * (q_ + r_) * (r_ + p_) / (p_ * q_ * r_) >= 8 - 1e-9                # ex3
    assert p_ / q_ + q_ / r_ + r_ / p_ >= 3 - 1e-9                                        # ex3 v1
    assert (p_ + q_ + r_) * (1 / p_ + 1 / q_ + 1 / r_) >= 9 - 1e-9                        # ex3 v2
assert 8 not in (6, 9, 2) and 3 not in (1, Rational(3, 2), 6) and 9 not in (3, 6, 27)

# ar-l1-convexite-iag-ex4 : (u_n) strictement croissante, (v_n) strictement décroissante, v_n - u_n = u_n / n, v_1 = 4.
U = lambda n: (1 + Rational(1, n))**n
V = lambda n: (1 + Rational(1, n))**(n + 1)
for n in range(1, 60):
    assert U(n) < U(n + 1) and V(n + 1) < V(n) and V(n) - U(n) == U(n) / n
    # sommes et produits des réels utilisés
    assert 1 + n * (1 + Rational(1, n)) == n + 2 and (1 + Rational(1, n))**n == U(n)
    assert 1 + (n + 1) * Rational(n, n + 1) == n + 1 and Rational(n, n + 1)**(n + 1) == 1 / V(n)
    assert Rational(n + 1, n + 2)**(n + 2) == 1 / V(n + 1)
assert V(1) == 4 and U(1) == 2 and abs(float(V(5)) - 2.99) < 5e-3 and float(V(5)) > float(E) > 2

# ar-l1-convexite-iag-ex5 : Young pour (p, q) = (3, 3/2) et (2, 2) sur une grille ; Hölder sur des exemples.
for pp in (Rational(3), Rational(2), Rational(5, 4)):
    qq = pp / (pp - 1)
    assert 1 / pp + 1 / qq == 1
    for av in (0, Rational(1, 3), 1, 2, 5):
        for bv in (0, Rational(1, 2), 1, 3):
            assert av * bv <= (Rational(av)**pp / pp + Rational(bv)**qq / qq).evalf() + 1e-12
    for _ in range(200):
        aa = [random.uniform(0, 3) for _ in range(4)]
        bb = [random.uniform(0, 3) for _ in range(4)]
        P, Q = float(pp), float(qq)
        assert sum(ai * bi for ai, bi in zip(aa, bb)) <= sum(ai**P for ai in aa)**(1 / P) * sum(bi**Q for bi in bb)**(1 / Q) + 1e-9

print("ar_l1_convexite : toutes les vérifications passent.")
