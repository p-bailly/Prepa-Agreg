"""Vérifications calculatoires des questions de QCM (content/qcm/ar-t-fonctions-usuelles.json) et des exercices ajoutés au
chapitre ar-t-fonctions-usuelles (content/exercices/analyse-reelle/terminale-ar-t-fonctions-usuelles.json).
Voir CLAUDE.md, « Contrôle de l'exactitude »."""
from sympy import (E, Rational as R, cos, diff, exp, limit, log, nsolve, oo, pi, simplify, sin, solve, solveset, sqrt,
                   symbols, Interval, Union, S, Reals, asin, nsimplify)

x = symbols('x', real=True)
X = symbols('X', real=True)

# ar-t-fonctions-usuelles-q1 : e^{3x-1} >= e^{-x} <=> x >= 1/4 ; leurres faux.
assert solveset(3 * x - 1 + x >= 0, x, Reals) == Interval(R(1, 4), oo)
assert (exp(3 * x - 1) - 1 / exp(x)).subs(x, R(3, 10)).evalf() > 0   # 0,3 est solution : [1/2,+oo[ et [1/3,+oo[ incomplets
assert (exp(3 * x - 1) - 1 / exp(x)).subs(x, 0).evalf() < 0           # 0 n'est pas solution : ]-oo,1/4] faux

# ar-t-fonctions-usuelles-q2 : (e^{2x}-1)/(e^x+1) = e^x - 1 (affirmation vraie).
assert simplify((exp(2 * x) - 1) / (exp(x) + 1) - (exp(x) - 1)) == 0

# ar-t-fonctions-usuelles-q3 : ln(sqrt(e^3)) - ln(1/e) = 5/2 ; leurres 1/2, 3/2, e^{5/2} faux.
v = simplify(log(sqrt(exp(3))) - log(1 / E))
assert v == R(5, 2) and all(v != l for l in (R(1, 2), R(3, 2), exp(R(5, 2))))

# ar-t-fonctions-usuelles-q4 : x(x-2) = 3 a pour solutions -1 et 3 (toutes deux avec x(x-2) > 0) ; la seconde équation
# (x > 2) n'a que 3 (affirmation fausse).
assert sorted(solve(x * (x - 2) - 3, x)) == [-1, 3]
assert all((s * (s - 2)) > 0 for s in (-1, 3))
assert simplify(log(3) + log(3 - 2) - log(3)) == 0

# ar-t-fonctions-usuelles-q5 : limite 3 ; leurres 0, 1, +oo faux.
f = (x**2 * exp(-x) + 3 * x) / (x + log(x))
assert limit(f, x, oo) == 3

# ar-t-fonctions-usuelles-q6 : x (ln x)^2 -> 0 en 0+ ; identité x(ln x)^2 = 4 (sqrt x ln sqrt x)^2.
assert limit(x * log(x)**2, x, 0, '+') == 0
t = symbols('t', positive=True)
assert simplify(t**2 * (log(t**2))**2 - 4 * (t * log(t))**2) == 0

# ar-t-fonctions-usuelles-q7 : sin x / (sqrt(1+x) - 1) -> 2 ; leurres 1/2, 1, 0 faux.
assert limit(sin(x) / (sqrt(1 + x) - 1), x, 0) == 2
assert limit((sqrt(1 + x) - 1) / x, x, 0) == R(1, 2)

# ar-t-fonctions-usuelles-q8 : 2 sin x = 1 a 4 solutions dans [0, 4 pi] ; leurres 2, 1, 8 faux.
sols = solveset(2 * sin(x) - 1, x, Interval(0, 4 * pi))
assert sols == {pi / 6, 5 * pi / 6, 13 * pi / 6, 17 * pi / 6}
assert sin(4 * pi) == 0

# ar-t-fonctions-usuelles-croissances-ex5 : f = ln x / x.
f = log(x) / x
assert limit(f, x, 0, '+') == -oo and limit(f, x, oo) == 0
assert simplify(diff(f, x) - (1 - log(x)) / x**2) == 0 and f.subs(x, E) == 1 / E
# question 3 (exercice et variantes) : pi^e < e^pi, 3^e < e^3, 2^e < e^2.
assert (E * log(pi) < pi) and (pi**E < E**pi)
assert (E * log(3) < 3) and (3**E < E**3)
assert (E * log(2) < 2) and (2**E < E**2)
# question 4 : x^2 = 2^x a exactement deux solutions positives (2 et 4), une troisième négative ; leurres 1, 3, 0 faux.
assert 2**4 == 4**2 and simplify(f.subs(x, 4) - f.subs(x, 2)) == 0
assert 0 < (log(2) / 2).evalf() < (1 / E).evalf()
neg = nsolve(x**2 - 2**x, x, -0.8)
assert -1 < neg < 0
# vérification numérique du nombre de racines positives : changements de signe de g = 2 ln x - x ln 2 sur une grille fine
def nb_racines(g, a, b, pas):
    vals = [g.subs(x, a + k * pas).evalf() for k in range(int((b - a) / pas) + 1)]
    return sum(1 for u, w in zip(vals, vals[1:]) if u == 0 or u * w < 0)
g = 2 * log(x) - x * log(2)
assert simplify(g.subs(x, 2)) == 0 and simplify(g.subs(x, 4)) == 0
assert nb_racines(g, R(1, 10), R(201, 10), R(3, 10)) == 2  # changements de signe (ou zéros exacts) sur une grille fine
assert nb_racines(3 * log(x) - x * log(3), R(1, 10), R(201, 10), R(3, 10)) == 2
# variante 1 : x^3 = 3^x : solutions 3 et beta dans ]1, e[ (environ 2,48) ; ln 3 / 3 < 1/e.
assert 0 < (log(3) / 3).evalf() < (1 / E).evalf()
beta = nsolve(3 * log(x) - x * log(3), x, 2.5)
assert 1 < beta < E and abs(beta - 2.478) < 0.01
assert simplify((3 * log(x) - x * log(3)).subs(x, 3)) == 0
# x^3 = 3^x n'a pas de solution négative (x^3 < 0 < 3^x) : 2 solutions au total.
# variante 2 : x^e = e^x <=> e ln x = x <=> f(x) = 1/e : seule solution e (maximum strict).
h = E * log(x) - x
assert simplify(h.subs(x, E)) == 0 and simplify(diff(h, x).subs(x, E)) == 0
assert all(h.subs(x, p).evalf() < 0 for p in (R(1, 10), 1, 2, R(27, 10), R(28, 10), 3, 5, 50))

# ar-t-fonctions-usuelles-limites-zero-ex5 : (e^x - cos x)/x -> 1 ; (e^{3x} - cos x)/x -> 3 ; (ln(1+x) - sin 2x)/x -> -1.
assert limit((exp(x) - cos(x)) / x, x, 0) == 1
assert limit((exp(3 * x) - cos(x)) / x, x, 0) == 3
assert limit((log(1 + x) - sin(2 * x)) / x, x, 0) == -1
assert limit((cos(x) - 1) / x, x, 0) == 0
# leurres : 0, 2 / 1, 4, 0 / 1, 3, -2 sont distincts des limites.

# ar-t-fonctions-usuelles-ln-ex6 : inéquations en X = ln x.
assert solveset(X**2 - X - 2 > 0, X, Reals) == Union(Interval.open(-oo, -1), Interval.open(2, oo))
assert solveset(log(x)**2 - log(x) - 2 > 0, x, Interval.open(0, oo)) == Union(Interval.open(0, exp(-1)), Interval.open(exp(2), oo))
assert solveset(X**2 + X - 6 <= 0, X, Reals) == Interval(-3, 2)
assert solveset(log(x)**2 + log(x) - 6 <= 0, x, Interval.open(0, oo)) == Interval(exp(-3), exp(2))
assert solveset(2 * X**2 - 3 * X + 1 < 0, X, Reals) == Interval.open(R(1, 2), 1)
assert solveset(2 * log(x)**2 - 3 * log(x) + 1 < 0, x, Interval.open(0, oo)) == Interval.open(sqrt(E), E)
# leurres de la variante 2 : ]e/2, e[ contient 1,5 qui n'est pas solution.
assert (2 * log(x)**2 - 3 * log(x) + 1).subs(x, R(3, 2)).evalf() > 0 and E / 2 < R(3, 2) < sqrt(E)
# Relecture indépendante (2026-10-29) : leurres des exercices ajoutés.
from sympy import E as _E, Interval as _I, Union as _U, exp as _exp, oo as _oo, sqrt as _sqrt
# ar-t-fonctions-usuelles-limites-zero-ex5 (et variantes) : leurres distincts de la limite (1, 3, -1).
assert all(_l != 1 for _l in (0, 2)) and all(_l != 3 for _l in (1, 4, 0)) and all(_l != -1 for _l in (1, 3, -2))
# ar-t-fonctions-usuelles-ln-ex6 (et variantes) : chaque leurre diffère de l'ensemble des solutions.
_sol = (_U(_I.open(0, _exp(-1)), _I.open(_exp(2), _oo)), _I(_exp(-3), _exp(2)), _I.open(_sqrt(_E), _E))
_leurres = ((_U(_I.open(-_oo, -1), _I.open(2, _oo)), _I.open(_exp(-1), _exp(2)), _I.open(_exp(2), _oo)),
            (_I(-3, 2), _U(_I.Lopen(0, _exp(-3)), _I.Ropen(_exp(2), _oo)), _I.Lopen(0, _exp(2))),
            (_I.open(0.5, 1), _U(_I.open(0, _sqrt(_E)), _I.open(_E, _oo)), _I.open(_E / 2, _E)))
for _s, _ls in zip(_sol, _leurres):
    assert all(_l != _s for _l in _ls)
print('enrichi_ar_t_fonctions_usuelles : tout est vérifié')
