"""Vérifications calculatoires des questions de QCM (content/qcm/ar-t-continuite.json) et des exercices ajoutés au chapitre
ar-t-continuite (content/exercices/analyse-reelle/terminale-ar-t-continuite.json). Voir CLAUDE.md, « Contrôle de l'exactitude »."""
from sympy import (E, Rational as R, diff, exp, factor, limit, nroots, oo, pi, real_roots, simplify, sin, solve, sqrt,
                   symbols, Poly, S, Interval, solveset, Reals, nsolve)

x = symbols('x', real=True)
n = symbols('n', positive=True, integer=True)

# ar-t-continuite-q1 : g(x) = e^x + x - 3, g(0) = -2 < 0, g(1) = e - 2 > 0 (affirmation vraie).
g = exp(x) + x - 3
assert g.subs(x, 0) == -2 and (g.subs(x, 1)).evalf() > 0

# ar-t-continuite-q2 : 1/x = 1/2 seulement pour x = 2, hors de [-1, 1] ; f(0) = 0 != 1/2 (affirmation fausse).
assert solve(1 / x - R(1, 2), x) == [2]

# ar-t-continuite-q3 : x^4 - 4x + 1 a exactement 2 racines réelles, l'une dans ]0,1[, l'autre dans ]1,2[ ; leurres 0, 1, 4 faux.
P = x**4 - 4 * x + 1
racines = real_roots(Poly(P, x))
assert len(racines) == 2
r1, r2 = sorted(float(r.evalf()) for r in racines)
assert 0 < r1 < 1 < r2 < 2
assert P.subs(x, 0) == 1 and P.subs(x, 1) == -2 and P.subs(x, 2) == 9
assert simplify(diff(P, x) - 4 * (x - 1) * (x**2 + x + 1)) == 0

# ar-t-continuite-q4 : exemple f(x) = 2 - 1/x - 1/x^2... on vérifie le raisonnement sur un exemple : f(x) = 2 - 1/x
# (continue, strictement croissante, -oo en 0+, 2 en +oo) : f(x) = 3 n'a pas de solution, f(x) = 1 en a une.
f = 2 - 1 / x
assert limit(f, x, 0, '+') == -oo and limit(f, x, oo) == 2 and simplify(diff(f, x) - 1 / x**2) == 0
assert solveset(f - 3, x, Interval.open(0, oo)) == S.EmptySet
assert solveset(f - 1, x, Interval.open(0, oo)) == {1}

# ar-t-continuite-q5 : (1-2n)/(n+3) -> -2, limite e^{-2} ; leurres 0, e^{1/3} (valeur pour n = 0), 1 faux.
assert limit((1 - 2 * n) / (n + 3), n, oo) == -2
assert limit(exp((1 - 2 * n) / (n + 3)), n, oo) == exp(-2)
assert exp(-2) != exp(R(1, 3)) and exp(-2) != 1 and exp(R(1, 3)) == exp(S(1 - 0) / (0 + 3))

# ar-t-continuite-q6 : u_n = 1/(n+1) -> 0, f(u_n) = 1 pour tout n (affirmation fausse).
assert limit(1 / (n + 1), n, oo) == 0

# ar-t-continuite-q7 : points fixes de (x^2+6)/5 : 2 et 3 ; u_n reste dans [0, 2] ; limite 2 ; leurres 3, 7/5 (= u_1), 6/5 faux.
fq = (x**2 + 6) / 5
assert sorted(solve(fq - x, x)) == [2, 3]
u = S(1)
for _ in range(6):  # exact sur quelques rangs (les dénominateurs grossissent vite)
    assert 0 <= u <= 2
    u = fq.subs(x, u)
u = 1.0
for _ in range(300):
    assert 0 <= u <= 2
    u = (u * u + 6) / 5
assert abs(u - 2) < 1e-6  # convergence géométrique de rapport 4/5 près de 2
assert fq.subs(x, 1) == R(7, 5) and fq.subs(x, 0) == R(6, 5) and fq.subs(x, 2) == 2

# ar-t-continuite-q8 : f(x) = x/2 + 1/x, f(x) - sqrt 2 = (x - sqrt 2)^2/(2x), point fixe positif sqrt 2 (affirmation vraie).
fh = x / 2 + 1 / x
assert simplify(fh - sqrt(2) - (x - sqrt(2))**2 / (2 * x)) == 0
assert simplify(fh - x - (2 - x**2) / (2 * x)) == 0
u = S(2)
for _ in range(6):
    u = fh.subs(x, u)
assert abs(float(u) - float(sqrt(2))) < 1e-12

# ar-t-point-fixe-ex5 : f(x) = (3x+2)/(x+2), u_0 = 0 ; v_n = (u_n - 2)/(u_n + 1) géométrique de raison 1/4.
def verifie_homographique(a, b, c, u0, l, raison, formule, leurres, intervalle):
    f = (a * x + b) / (x + 2)
    assert simplify(diff(f, x) - (2 * a - b) / (x + 2)**2) == 0 and 2 * a - b > 0
    assert sorted(solve(f - x, x)) == sorted([l, -1])
    lo, hi = intervalle
    if hi is not None:
        assert lo <= f.subs(x, lo) and f.subs(x, hi) <= hi
    else:
        assert f.subs(x, lo) == lo
    # v_{n+1} = raison * v_n (identité en u)
    v = (x - l) / (x + 1)
    assert simplify(((f - l) / (f + 1)) - raison * v) == 0
    u = S(u0)
    termes = []
    for m in range(8):
        termes.append(u)
        assert u == formule(m)
        u = f.subs(x, u)
    # monotonie et encadrement
    for m in range(7):
        if u0 < l:
            assert 0 <= termes[m] <= termes[m + 1] <= l
        else:
            assert l <= termes[m + 1] <= termes[m]
    for lf in leurres:  # chaque leurre diffère de u_n pour au moins un des premiers rangs
        assert any(lf(m) != termes[m] for m in range(4))
    return termes

t = verifie_homographique(3, 2, None, 0, 2, R(1, 4), lambda m: R(2 * 4**m - 2, 4**m + 2),
                          [lambda m: 2 - 2 * R(1, 4)**m, lambda m: R(2 * 4**m - 2, 4**m + 1), lambda m: R(2 * 4**m + 2, 4**m - 2)], (0, 2))
assert t[1] == 1 and (S(0) - 2) / (0 + 1) == -2
# variante 1 : f(x) = (4x+3)/(x+2), u_0 = 0, limite 3, raison 1/5, v_0 = -3.
verifie_homographique(4, 3, None, 0, 3, R(1, 5), lambda m: R(3 * 5**m - 3, 5**m + 3),
                      [lambda m: 3 - 3 * R(1, 5)**m, lambda m: R(3 * 5**m - 3, 5**m + 1), lambda m: R(3 * 5**m + 3, 5**m - 3)], (0, 3))
# variante 2 : f(x) = (5x+4)/(x+2), u_0 = 6, limite 4, raison 1/6, v_0 = 2/7, u_1 = 17/4.
t = verifie_homographique(5, 4, None, 6, 4, R(1, 6), lambda m: R(28 * 6**m + 2, 7 * 6**m - 2),
                          [lambda m: 4 + 2 * R(1, 6)**m, lambda m: R(28 * 6**m - 2, 7 * 6**m + 2), lambda m: R(4 * 6**m + 2, 6**m - 2)], (4, None))
assert t[1] == R(17, 4) and (S(6) - 4) / (6 + 1) == R(2, 7)

# ar-t-tvi-ex5 : g(x) = e^x + x^2 - c ; signes aux bornes et encadrements d'amplitude 0,1.
for c, (a, b), (lo, hi) in ((2, (0, 1), (R(5, 10), R(6, 10))), (3, (0, 1), (R(8, 10), R(9, 10))), (4, (1, 2), (S(1), R(11, 10)))):
    g = exp(x) + x**2 - c
    assert g.subs(x, -2).evalf() > 0 and g.subs(x, -1).evalf() < 0
    assert g.subs(x, a).evalf() < 0 and g.subs(x, b).evalf() > 0
    assert g.subs(x, lo).evalf() < 0 and g.subs(x, hi).evalf() > 0
    # unicité dans ]a, b[ : g' = e^x + 2x > 0 sur [0, +oo[
    assert all((exp(x) + 2 * x).subs(x, R(j, 10) + a).evalf() > 0 for j in range(11))
    alpha = nsolve(g, x, (lo + hi) / 2)
    assert lo < alpha < hi
    # leurres : intervalles disjoints de ]lo, hi[, donc faux
# valeurs citées dans les corrections
assert abs(float(exp(R(1, 2)) - R(175, 100)) + 0.10) < 0.01 and abs(float(exp(R(6, 10)) - R(164, 100)) - 0.18) < 0.01
assert abs(float(exp(R(8, 10)) - R(236, 100)) + 0.13) < 0.01 and abs(float(exp(R(9, 10)) - R(219, 100)) - 0.27) < 0.01
assert abs(float(E - 3) + 0.28) < 0.01 and abs(float(exp(R(11, 10)) - R(279, 100)) - 0.21) < 0.01

# ar-t-continuite-image-suite-ex5 : f(a_n) = sin(n pi) = 0, f(b_n) = sin((4n+1)pi/2) = 1 ; |x sin(1/x)| <= |x|.
for m in range(1, 20):
    assert sin(1 / (1 / (m * pi))) == 0
    assert simplify(sin(1 / (2 / ((4 * m + 1) * pi)))) == 1
assert limit(1 / (n * pi), n, oo) == 0 and limit(2 / ((4 * n + 1) * pi), n, oo) == 0
assert limit(x * sin(1 / x), x, 0) == 0
# Relecture indépendante (2026-10-29) : contrôles complémentaires par itération.
from sympy import Rational as _R, sqrt as _sqrt
# ar-t-continuite-q7 : u_0 = 1, u_{n+1} = (u_n^2+6)/5 reste dans [6/5, 2] dès n = 1 et tend vers 2.
_u = 1.0
for _m in range(200):
    _u = (_u**2 + 6) / 5
    assert 1.2 - 1e-12 <= _u <= 2
assert abs(_u - 2) < 1e-6
# ar-t-continuite-q8 : u_0 = 2, u_{n+1} = u_n/2 + 1/u_n décroît vers sqrt 2.
_u = _R(2)
for _m in range(6):
    _v = _u / 2 + 1 / _u
    assert _sqrt(2) <= _v <= _u
    _u = _v
assert abs(float(_u) - 2**0.5) < 1e-12
# ar-t-point-fixe-ex5 (et variantes) : la formule explicite coïncide avec l'itération pour n <= 8.
for (_a, _b, _u0, _formule) in ((3, 2, _R(0), lambda m: _R(2 * 4**m - 2, 4**m + 2)),
                                (4, 3, _R(0), lambda m: _R(3 * 5**m - 3, 5**m + 3)),
                                (5, 4, _R(6), lambda m: _R(28 * 6**m + 2, 7 * 6**m - 2))):
    _u = _u0
    for _m in range(9):
        assert _u == _formule(_m)
        _u = (_a * _u + _b) / (_u + 2)
print('enrichi_ar_t_continuite : tout est vérifié')
