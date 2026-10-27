"""Vérifications calculatoires du chapitre « ar-l1-continuite » (caractérisation séquentielle, image d'un intervalle,
réciproque d'une fonction continue strictement monotone, théorème de Heine) : cartes et exercices de
content/*/analyse-reelle/l1-ar-l1-continuite.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice concerné (v1, v2 = variantes 1 et 2).
Lancer : python3 verification/ar_l1_continuite.py (doit finir sans erreur).
"""

from sympy import (Abs, E, Interval, Rational, S, atan, cos, exp, floor, limit, log, oo, pi, simplify, sin, solve,
                   sqrt, symbols)
from sympy.calculus.util import function_range

x, y, t = symbols("x y t", real=True)
n = symbols("n", integer=True, positive=True)
h = symbols("h", positive=True)


def egal(u, v):
    return simplify(u - v) == 0


# ===========================================================================
# Carte ar-l1-continuite-sequentielle
# ===========================================================================
# Exemple : sin(2nπ) = 0 et sin(2nπ + π/2) = 1.
assert sin(2 * n * pi) == 0 and sin(2 * n * pi + pi / 2) == 1
# Contre-exemple 1 / affirmation 1 : sin(π n) = 0, sin(π (4n+1)/2) = 1.
assert sin(pi * n) == 0 and sin(pi * (4 * n + 1) / 2) == 1
# Affirmation 1 : 2/(4n+1) est dans ]0,1] seulement pour n ≥ 1 (pour n = 0, il vaut 2).
assert all(0 < Rational(2, 4 * k + 1) <= 1 for k in range(1, 20)) and Rational(2, 1) > 1
# Exemple : r_n = floor(10^n x)/10^n vérifie x - 10^-n < r_n ≤ x (contrôle numérique pour x = sqrt 2, n ≤ 8).
for k in range(9):
    r = floor(10**k * sqrt(2)) / S(10)**k
    assert sqrt(2) - S(10)**(-k) < r <= sqrt(2)

# ar-l1-continuite-sequentielle-ex1 et variantes : valeurs exactes des suites (1, -1 ; 1, -1 ; 1, 0).
assert cos(2 * n * pi) == 1 and cos((2 * n + 1) * pi) == -1                                  # ex1
assert sin(2 * n * pi + pi / 2) == 1 and sin(2 * n * pi - pi / 2) == -1                      # ex1 v1
assert 2 * pi * 1 - pi / 2 > 0                                                               # ex1 v1 : v_n > 0
assert cos(n * pi)**2 == 1 and cos(n * pi + pi / 2)**2 == 0                                  # ex1 v2

# ar-l1-continuite-sequentielle-ex2 et variantes : f(√2) = 2-√2 ; f(√3) = √3/2 ; f(π) = π-1 ; leurres faux.
v = sqrt(2)**2 - sqrt(2)
assert egal(v, 2 - sqrt(2)) and not egal(v, 0) and not egal(v, 2 + sqrt(2))                 # ex2
v = 2 * sqrt(3) / (1 + 3)
assert egal(v, sqrt(3) / 2) and not egal(v, 0) and not egal(v, sqrt(3))                      # ex2 v1
v = Abs(pi - 1)
assert egal(v, pi - 1) and not egal(v, 0) and not egal(v, 1 - pi)                            # ex2 v2

# ar-l1-continuite-sequentielle-ex3 et variantes : points où les deux expressions coïncident.
assert set(solve(x**2 - (x + 2), x)) == {-1, 2}                                              # ex3
assert set(solve(x - (1 - x), x)) == {Rational(1, 2)}                                        # ex3 v1
assert set(solve(x**3 - x, x)) == {-1, 0, 1}                                                 # ex3 v2

# ar-l1-continuite-sequentielle-ex4 Q3 : g(1-1/n) = sin(nπ) = 0, g(1-2/(4n+1)) = 1, et v_n ∈ [0,1[.
assert sin(pi / (1 - (1 - S(1) / n))) == 0
assert simplify(sin(pi / (1 - (1 - S(2) / (4 * n + 1))))) == 1
assert 0 <= 1 - Rational(2, 5) < 1

# ===========================================================================
# Carte ar-l1-continuite-image-intervalle
# ===========================================================================
# Cours / affirmation 1 : sin(]0,2π[) = [-1,1] ; 1/x envoie ]0,1] sur [1,+∞[.
assert function_range(sin(x), x, Interval.open(0, 2 * pi)) == Interval(-1, 1)
assert function_range(1 / x, x, Interval.Lopen(0, 1)) == Interval(1, oo)
# Exemple : x/(1+|x|) a pour image ]-1,1[ (étude sur R+ puis imparité).
assert function_range(x / (1 + x), x, Interval(0, oo)) == Interval.Ropen(0, 1)
# Affirmation 4 : u + (v-u)/√2 irrationnel pour u ≠ v rationnels (√2 = (v-u)/(w-u) serait rationnel) ; exemple u=0, v=1.
assert not (1 / sqrt(2)).is_rational

# ar-l1-continuite-image-intervalle-ex1 et variantes : images ; leurres différents.
assert function_range(x**2, x, Interval.Ropen(-1, 2)) == Interval.Ropen(0, 4)                 # ex1
assert function_range(x**2, x, Interval.Lopen(-3, 1)) == Interval.Ropen(0, 9)                 # ex1 v1
assert function_range(cos(x), x, Interval.open(0, 3 * pi / 2)) == Interval.Ropen(-1, 1)       # ex1 v2
for bonne, leurres in [(Interval.Ropen(0, 4), [Interval.Ropen(1, 4), Interval(0, 4), Interval.open(1, 4)]),
                       (Interval.Ropen(0, 9), [Interval.Ropen(1, 9), Interval(0, 9), Interval.open(0, 9)]),
                       (Interval.Ropen(-1, 1), [Interval.open(0, 1), Interval.open(-1, 1), Interval(-1, 1)])]:
    assert all(l != bonne for l in leurres)

# ar-l1-continuite-image-intervalle-ex2 et variantes : nombre de valeurs possibles = nombre de fonctions constantes.
yy = symbols("yy")
assert len(solve(yy**2 - 1, yy)) == 2                                                         # ex2
assert len(solve(yy**3 - yy, yy)) == 3                                                        # ex2 v1
assert len(solve((yy**2 - 1) * (yy**2 - 4), yy)) == 4                                         # ex2 v2

# ar-l1-continuite-image-intervalle-ex3 (exemples des cas possibles).
f = (1 + sin(2 * pi * x)) / 2
assert f.subs(x, Rational(1, 4)) == 1 and f.subs(x, Rational(3, 4)) == 0                     # ex3 (a)
assert limit(sin(pi * (x - S(1) / 2)) / cos(pi * (x - S(1) / 2)), x, 0, "+") == -oo           # ex3 (c)
assert limit(sin(pi * (x - S(1) / 2)) / cos(pi * (x - S(1) / 2)), x, 1, "-") == oo
f = (1 + sin(x)) / 2
assert f.subs(x, pi / 2) == 1 and f.subs(x, -pi / 2) == 0                                     # ex3 (d)
assert limit(x**2 / (1 - x**2), x, 1, "-") == oo                                              # ex3 v1 (b)
assert function_range((1 - x**2) / (1 + x**2), x, Interval(0, oo)) == Interval.Lopen(-1, 1)   # ex3 v1 (c)
assert sin(pi / 2) == 1 and sin(3 * pi / 2) == -1                                             # ex3 v1 (d)
assert (10 * x - 5).subs(x, 0) == -5 and (10 * x - 5).subs(x, 1) == 5                         # ex3 v2 (b)
assert function_range(1 / x - 1, x, Interval.open(0, 1)) == Interval.open(0, oo)              # ex3 v2 (c)

# ===========================================================================
# Carte ar-l1-continuite-reciproque
# ===========================================================================
# Exemple : f(x) = x + sin(x)/2, f' = 1 + cos(x)/2 ≥ 1/2.
assert function_range(1 + cos(x) / 2, x, Interval(0, 2 * pi)) == Interval(Rational(1, 2), Rational(3, 2))
# Affirmation 3 : le taux de la racine cubique en 0 tend vers +∞.
assert limit(h**Rational(1, 3) / h, h, 0, "+") == oo

# ar-l1-continuite-reciproque-ex1 et variantes : f(bonne) = valeur, f(leurre) ≠ valeur.
for f, val, bonne, leurres in [(x**3 + x, 10, 2, [1010, Rational(1, 1010), 10**Rational(1, 3)]),                # ex1
                               (x**5 + 2 * x, 36, 2, [18, Rational(1, 36), 36**Rational(1, 5)]),                # ex1 v1
                               (exp(x) + x, 1, 0, [1, E + 1, -1])]:                                             # ex1 v2
    assert egal(f.subs(x, bonne), val)
    assert all(not egal(f.subs(x, l), val) for l in leurres)
assert egal((x**3 + x).subs(x, 10), 1010)

# ar-l1-continuite-reciproque-ex2 et variantes : f'(x), limites, f(g(y)) = y ; leurres de formule faux.
for f, g, a0, lim_a, lim_b, leurres_g in [
        ((2 * x + 1) / (x + 3), (3 * y - 1) / (2 - y), -3, -oo, 2, [(y + 3) / (2 * y + 1), (3 * y + 1) / (2 - y)]),     # ex2
        ((x - 1) / (x + 2), (2 * y + 1) / (1 - y), -2, -oo, 1, [(y + 2) / (y - 1), (2 * y - 1) / (1 - y)])]:           # ex2 v1
    num = simplify(f.diff(x) * (x - a0)**2)
    assert num.is_number and num > 0          # f'(x) = num/(x - a0)^2 > 0 (5 puis 3)
    assert limit(f, x, a0, "+") == lim_a and limit(f, x, oo) == lim_b
    assert egal(f.subs(x, g), y)
    assert all(not egal(f.subs(x, l), y) for l in leurres_g)
f = x**2 - 4 * x + 1                                                                                                  # ex2 v2
assert egal(f, (x - 2)**2 - 3) and f.subs(x, 2) == -3
assert egal(f.subs(x, 2 + sqrt(y + 3)), y)
assert not egal(f.subs(x, 2 + sqrt(y - 3)), y)
assert (2 - sqrt(1 + 3)) < 2      # 2 - √(y+3) n'est pas dans [2, +∞[ (ici y = 1)

# ar-l1-continuite-reciproque-ex3 et variantes : g(valeur), puis lim g(y)/y = lim x/f(x) (y = f(x), x → ∓∞).
assert (2 * x + exp(x)).subs(x, 0) == 1 and limit(x / (2 * x + exp(x)), x, -oo) == Rational(1, 2)                      # ex3
assert (3 * x - exp(-x)).subs(x, 0) == -1 and limit(x / (3 * x - exp(-x)), x, oo) == Rational(1, 3)                    # ex3 v1
assert egal((x + 2 * atan(x)).subs(x, 1), 1 + pi / 2) and limit(x / (x + 2 * atan(x)), x, oo) == 1                     # ex3 v2
for bonne, leurres in [(Rational(1, 2), [2, 0, 1]), (Rational(1, 3), [3, 0, 1]), (1, [0, Rational(1, 2), pi / 2])]:
    assert all(not egal(l, bonne) for l in leurres)

# ar-l1-continuite-reciproque-ex4 Q3 : x ↦ 5 - x envoie [2,3] sur [2,3], disjoint de [0,1].
assert (5 - x).subs(x, 2) == 3 and (5 - x).subs(x, 3) == 2

# ===========================================================================
# Carte ar-l1-continuite-heine
# ===========================================================================
# Exemple : (√x + √(y-x))^2 - y = 2√(x(y-x)) ≥ 0 ; taux de √ en 0 : 1/√x non borné.
xx, dd = symbols("xx dd", nonnegative=True)
assert egal((sqrt(xx) + sqrt(dd))**2 - (xx + dd), 2 * sqrt(xx * dd))
assert limit(sqrt(h) / h, h, 0, "+") == oo
# Contre-exemples : 1/n - 1/(n+1) = 1/(n(n+1)) ; (n+1/n)^2 - n^2 = 2 + 1/n^2 ; sin(1/x) aux deux suites.
assert egal(S(1) / n - S(1) / (n + 1), S(1) / (n * (n + 1)))
assert egal((n + S(1) / n)**2 - n**2, 2 + S(1) / n**2)
assert sin(2 * n * pi + pi / 2) - sin(2 * n * pi) == 1

# ar-l1-continuite-heine-ex1 et variantes : la bonne valeur de η convient, les leurres non (pire cas à l'extrémité).
e = symbols("e", positive=True)
# ex1 : x^2 sur [0,3], ε = 1/2 : |x^2-y^2| ≤ 6|x-y| ; pire écart 9 - (3-η)^2 = 6η - η^2.
assert 6 * Rational(1, 12) == Rational(1, 2)
assert egal(9 - (3 - e)**2, 6 * e - e**2)
for eta, val in [(Rational(1, 6), Rational(35, 36)), (Rational(1, 10), Rational(59, 100)), (Rational(1, 2), Rational(11, 4))]:
    assert (6 * eta - eta**2) == val and val > Rational(1, 2)
# ex1 v1 : 1/x sur [1,2], ε = 1/10 : |1/x-1/y| ≤ |x-y| ; pire écart 1 - 1/(1+η) = η/(1+η).
for eta, val in [(Rational(1, 5), Rational(1, 6)), (Rational(1, 4), Rational(1, 5)), (Rational(1, 2), Rational(1, 3))]:
    assert eta / (1 + eta) == val and val > Rational(1, 10)
assert Rational(1, 10) / (1 + Rational(1, 10)) <= Rational(1, 10)
# ex1 v2 : x^3 sur [0,2], ε = 1/10 : |x^3-y^3| ≤ 12|x-y| ; pire écart 8 - (2-η)^3 = 12η - 6η^2 + η^3.
assert egal(8 - (2 - e)**3, 12 * e - 6 * e**2 + e**3) and 12 * Rational(1, 120) == Rational(1, 10)
for eta in [Rational(1, 12), Rational(1, 10), Rational(1, 40)]:
    assert 12 * eta - 6 * eta**2 + eta**3 > Rational(1, 10)
assert abs(float(12 * Rational(1, 12) - 6 * Rational(1, 144) + Rational(1, 1728)) - 0.959) < 1e-3
assert abs(float(12 * Rational(1, 10) - 6 * Rational(1, 100) + Rational(1, 1000)) - 1.141) < 1e-3
assert abs(float(12 * Rational(1, 40) - 6 * Rational(1, 1600) + Rational(1, 64000)) - 0.296) < 1e-3

# ar-l1-continuite-heine-ex2 et variantes : limites de x_n - y_n et de f(x_n) - f(y_n).
assert egal((n + S(1) / n**2)**3 - n**3, 3 + S(3) / n**3 + S(1) / n**6)                      # ex2
assert limit((n + S(1) / n**2)**3 - n**3, n, oo) == 3
assert limit(S(1) / n**2, n, oo) == 0
assert limit(log(1 + S(1) / n), n, oo) == 0                                                  # ex2 v1
assert egal(exp(log(n + 1)) - exp(log(n)), 1)
assert limit((n + 1 / sqrt(n))**Rational(3, 2) - n**Rational(3, 2), n, oo) == Rational(3, 2)  # ex2 v2
assert limit(1 / sqrt(n), n, oo) == 0

# ar-l1-continuite-heine-ex3 Q2 : x_n^2 - y_n^2 = π/2, x_n - y_n → 0, écart des sinus = 1.
xn, yn = sqrt(2 * n * pi + pi / 2), sqrt(2 * n * pi)
assert egal(xn**2 - yn**2, pi / 2) and limit(xn - yn, n, oo) == 0
assert sin(xn**2) - sin(yn**2) == 1

# ar-l1-continuite-heine-ex4 Q4 : x sin(1/x) → 0 en 0 ; sin(1/x) vaut 0 et 1 sur deux suites.
assert limit(x * sin(1 / x), x, 0, "+") == 0
assert sin(2 * n * pi) == 0 and sin(2 * n * pi + pi / 2) == 1

print("ar_l1_continuite : toutes les vérifications passent.")
