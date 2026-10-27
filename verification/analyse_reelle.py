"""Vérifications calculatoires du thème « analyse-reelle » (cartes et exercices terminale et l1).

Chaque assert cite l'identifiant de la carte ou de l'exercice concerné.
Lancer : python3 verification/analyse_reelle.py (doit finir sans erreur).
"""

import math

import sympy
from sympy import (Abs, Rational, atan, cos, diff, factor, limit, log, oo, pi, real_roots, simplify, sin, solve,
                   sqrt, symbols)

x, n, t = symbols("x n t", real=True)
np_ = symbols("n", positive=True, integer=True)


def proche(a, b, tol=1e-12):
    return abs(float(a) - float(b)) < tol


# ---------------------------------------------------------------------------
# ar-t-limite-monotone
# ---------------------------------------------------------------------------

# Exemple : u_n = somme des 1/k^2 reste < 2 (majoration télescopique 1 + 1 - 1/n).
s = 0.0
for k in range(1, 5001):
    s += 1 / k**2
    assert s <= 2 - 1 / k + 1e-12
# Contre-exemple / affirmation 2 : 1 - 1/(n+1) est majorée par 2 mais tend vers 1.
assert limit(1 - 1 / (np_ + 1), np_, oo) == 1
# Contre-exemple : sqrt(n+1) - sqrt(n) -> 0 alors que sqrt(n) -> +oo.
assert limit(sqrt(np_ + 1) - sqrt(np_), np_, oo) == 0 and limit(sqrt(np_), np_, oo) == oo

# ar-t-limite-monotone-ex1 : u_n = somme_{k=0}^n 1/(2^k+1) < 2 - 1/2^n - 1/2 < 3/2, limite ≈ 1,26 (donc ≠ 2, ≤ 2).
u = 0.0
for k in range(0, 80):
    u += 1 / (2**k + 1)
    assert u < 2 - 1 / 2**k - 0.5 + 1e-12 and u < 1.5
assert 1.26 < u < 1.27
# Le terme général est majoré par la série géométrique de somme 2.
assert sum(Rational(1, 2**k) for k in range(0, 30)) < 2

# ar-t-limite-monotone-ex2 : a_n = somme 1/sqrt(k) - 2 sqrt(n) décroissante, > -2, a_1 = -1, limite zeta(1/2) ≈ -1,46.
a_prec = None
S = 0.0
for k in range(1, 200001):
    S += 1 / math.sqrt(k)
    a = S - 2 * math.sqrt(k)
    assert a > -2
    if a_prec is not None:
        assert a < a_prec
    a_prec = a
    if k == 1:
        assert proche(a, -1)
zeta_demi = float(sympy.zeta(Rational(1, 2)))
assert abs(a_prec - zeta_demi) < 0.01 and proche(zeta_demi, -1.4603545088, 1e-8)
assert abs(zeta_demi - (-2)) > 0.5  # le leurre « de limite -2 » est faux
# Inégalité clé : 1/sqrt(k) > 2(sqrt(k+1) - sqrt(k)) (différence > 0 pour k = 1..1000).
for k in range(1, 1001):
    assert 1 / math.sqrt(k) - 2 * (math.sqrt(k + 1) - math.sqrt(k)) > 0

# ar-t-limite-monotone-ex3 : sqrt(n+1) - sqrt(n) < 1/2 pour n >= 1 ; 2 cos(1/2) > 1,7.
assert 1 / (sqrt(2) + 1) < Rational(1, 2)
for k in range(1, 1000):
    assert math.sqrt(k + 1) - math.sqrt(k) < 0.5
assert 2 * math.cos(0.5) > 1.7 and math.cos(0.5) > 0.87
# Identité cos a - cos b = -2 sin((a+b)/2) sin((a-b)/2).
a_, b_ = symbols("a b", real=True)
assert simplify(cos(a_) - cos(b_) + 2 * sin((a_ + b_) / 2) * sin((a_ - b_) / 2)) == 0

# ---------------------------------------------------------------------------
# ar-t-gendarmes
# ---------------------------------------------------------------------------

# Exemple : floor(n sqrt 2)/n -> sqrt 2 (encadrement numérique).
for k in (10, 1000, 10**6):
    q = math.floor(k * math.sqrt(2)) / k
    assert math.sqrt(2) - 1 / k < q <= math.sqrt(2)
# Affirmation 1 : 0 < 1/(n+1) et les deux limites valent 0.
assert limit(1 / (np_ + 1), np_, oo) == 0

# ar-t-gendarmes-ex1 : les bornes (3n-1)/(n+2) et (3n+1)/(n+2) tendent vers 3 ; leurres 3/2 et 4 faux.
assert limit((3 * np_ - 1) / (np_ + 2), np_, oo) == 3 and limit((3 * np_ + 1) / (np_ + 2), np_, oo) == 3
assert proche((3 * 10**8 + math.cos(10**8)) / (10**8 + 2), 3, 1e-6)

# ar-t-gendarmes-ex2 : (2^n + 3^n)^(1/n) -> 3 (et non 5, sqrt 6, 1).
assert limit(3 * 2 ** (1 / np_), np_, oo) == 3
for k in (1, 5, 50, 300):
    v = (2**k + 3**k) ** (1 / k)
    assert 3 <= v <= 3 * 2 ** (1 / k) + 1e-9
assert abs((2**300 + 3**300) ** (1 / 300) - 3) < 1e-9
assert all(abs(3 - c) > 0.4 for c in (5, math.sqrt(6), 1))

# ar-t-gendarmes-ex3 : (l - eps = l' + eps) avec eps = (l - l')/2.
l1, l2 = symbols("l1 l2", real=True)
assert simplify((l1 - (l1 - l2) / 2) - (l2 + (l1 - l2) / 2)) == 0

# ---------------------------------------------------------------------------
# ar-t-tvi
# ---------------------------------------------------------------------------

# Exemple : x^3 + x - 1 : f(0) = -1, f(1) = 1.
f = x**3 + x - 1
assert f.subs(x, 0) == -1 and f.subs(x, 1) == 1
# Contre-exemple / affirmation 1 : x^3 - x sur [-2, 2] : f(-2) = -6, f(2) = 6, racines -1, 0, 1.
f = x**3 - x
assert f.subs(x, -2) == -6 and f.subs(x, 2) == 6 and f.subs(x, -2) * f.subs(x, 2) == -36
assert sorted(solve(f, x)) == [-1, 0, 1]
# Contre-exemple : sin(1/x) vaut 1 en 1/(pi/2 + 2 n pi).
assert simplify(sin(1 / (1 / (pi / 2 + 2 * pi * np_)))) == 1

# ar-t-tvi-ex1 : x^5 + 3x - 1 : f(0) = -1, f(1/2) = 17/32, f' > 0, une seule racine réelle, dans ]0, 1/2[.
f = x**5 + 3 * x - 1
assert f.subs(x, 0) == -1 and f.subs(x, Rational(1, 2)) == Rational(17, 32)
racines = real_roots(f)
assert len(racines) == 1
r = float(racines[0])
assert 0 < r < 0.5  # les leurres [1/2,1], [1,2], [-1,0] ne contiennent pas la racine
# f'(x) = 5x^4 + 3 > 0.
assert diff(f, x) == 5 * x**4 + 3

# ar-t-tvi-ex2 : (1+x)/2 = x seulement pour x = 1 (pas de point fixe dans [0,1[).
assert solve((1 + x) / 2 - x, x) == [1]

# ar-t-tvi-ex3 : x^3 - 3x + 1 : valeurs en -2, -1, 0, 1, 2 et exactement 3 racines réelles.
P = x**3 - 3 * x + 1
assert [P.subs(x, v) for v in (-2, -1, 0, 1, 2)] == [-1, 3, 1, -1, 3]
rr = [float(z) for z in real_roots(P)]
assert len(rr) == 3
assert -2 < rr[0] < -1 and 0 < rr[1] < 1 and 1 < rr[2] < 2

# ar-t-tvi-ex4 : contre-exemple f(x) = x : f(c + 1/2) - f(c) = 1/2.
assert (x + Rational(1, 2)) - x == Rational(1, 2)
# Télescopage : somme des g(k/n) = f(1) - f(0) pour une fonction test.
ftest = sympy.Lambda(x, sin(7 * x) + x**2)
for m in (2, 3, 5):
    tele = sum(ftest(Rational(k + 1, m)) - ftest(Rational(k, m)) for k in range(m))
    assert simplify(tele - (ftest(1) - ftest(0))) == 0

# ---------------------------------------------------------------------------
# ar-t-point-fixe
# ---------------------------------------------------------------------------

# Exemple : sqrt(x+6) : point fixe 3, u_{n+1}^2 - u_n^2 = (3-u)(u+2), f(0) = sqrt 6, f(3) = 3.
assert solve(sqrt(x + 6) - x, x) == [3]
assert sympy.expand((3 - x) * (x + 2) - (x + 6 - x**2)) == 0
u = 1.0
for _ in range(60):
    u = math.sqrt(u + 6)
assert proche(u, 3)
# Contre-exemple : f(x) = (1+x)/2 si x < 1, 2 sinon ; u_0 = 0 donne u_n = 1 - 2^{-n} < 1, et f(1) = 2.
def f_saut(v):
    return (1 + v) / 2 if v < 1 else 2
u = 0.0
for k in range(1, 40):
    u = f_saut(u)
    assert proche(u, 1 - 2.0**-k) and u < 1
assert f_saut(1) == 2 and f_saut(2) == 2  # 2 est point fixe ; pas d'autre point fixe dans [0,1[
assert solve((1 + x) / 2 - x, x) == [1]
# Contre-exemple : f(x) = -x, u_0 = 1 : (-1)^n.
assert [(-1) ** k for k in range(4)] == [1, -1, 1, -1]

# ar-t-point-fixe-ex1 : (l + 4)/3 = l donne l = 2 (et pas 4/3, 4, 0).
assert solve((x + 4) / 3 - x, x) == [2]
u = 0.0
for _ in range(80):
    u = (u + 4) / 3
assert proche(u, 2)

# ar-t-point-fixe-ex2 : l^2 = 2l + 3 donne l = -1 ou 3 ; seule 3 est dans [1,3] ; u_1 = sqrt 5.
assert sorted(solve(x**2 - 2 * x - 3, x)) == [-1, 3]
assert sympy.expand((3 - x) * (x + 1) - (2 * x + 3 - x**2)) == 0
u = 1.0
assert proche(math.sqrt(2 * u + 3), math.sqrt(5))
for _ in range(80):
    u_suiv = math.sqrt(2 * u + 3)
    assert 1 <= u <= u_suiv <= 3
    u = u_suiv
assert proche(u, 3)

# ar-t-point-fixe-ex3 : x - x^2 n'a que 0 comme point fixe ; f(]0,1[) ⊂ ]0,1/4] ; v_0 = 2 donne v_1 = -2 puis -> -oo.
assert solve(x - x**2 - x, x) == [0]
assert sympy.maximum(x - x**2, x, sympy.Interval(0, 1)) == Rational(1, 4)
v = 2
vs = []
for k in range(1, 8):
    v = v - v**2
    vs.append(v)
    assert v <= -2 - 4 * (k - 1)
assert vs[0] == -2 and vs[1] == -6
u = 0.5
for _ in range(10**5):
    u = u - u * u
assert 0 < u < 1e-4

# ar-t-point-fixe-ex4 : f o f = (1+x)/(2+x), point fixe (sqrt5 - 1)/2, u_2 = 1/2, u_3 = 2/3.
F = 1 / (1 + x)
G = simplify(F.subs(x, F))
assert simplify(G - (1 + x) / (2 + x)) == 0
pts = solve(G - x, x)
assert set(pts) == {(-1 + sqrt(5)) / 2, (-1 - sqrt(5)) / 2}
assert [p for p in pts if 0 <= float(p) <= 1] == [(-1 + sqrt(5)) / 2]
u = [0.0]
for _ in range(60):
    u.append(1 / (1 + u[-1]))
assert u[1] == 1 and proche(u[2], 0.5) and proche(u[3], 2 / 3)
assert all(u[2 * k] <= u[2 * k + 2] for k in range(25)) and all(u[2 * k + 1] >= u[2 * k + 3] for k in range(25))
assert proche(u[-1], (math.sqrt(5) - 1) / 2)
assert not proche(u[-1], 0.5, 1e-3) and not proche(u[-1], (1 + math.sqrt(5)) / 2, 1e-3)

# ---------------------------------------------------------------------------
# ar-t-derivee-monotonie
# ---------------------------------------------------------------------------

# Contre-exemple : -1/x a pour dérivée 1/x^2 > 0 mais f(-1) = 1 > f(1) = -1.
assert simplify(diff(-1 / x, x) - 1 / x**2) == 0 and (-1 / x).subs(x, -1) == 1 and (-1 / x).subs(x, 1) == -1
# x^3 : dérivée nulle en 0.
assert diff(x**3, x).subs(x, 0) == 0
# Exemple : (x - sin x)' = 1 - cos x, nulle en 2k pi.
assert diff(x - sin(x), x) == 1 - cos(x)

# ar-t-derivee-monotonie-ex1 : f' = 3(x-1)(x+1) ; décroissante sur [-1,1] ; zéros de f : 0, ±sqrt 3.
f = x**3 - 3 * x
assert sympy.expand(diff(f, x) - 3 * (x - 1) * (x + 1)) == 0
assert set(solve(f, x)) == {0, sqrt(3), -sqrt(3)}
assert all(diff(f, x).subs(x, Rational(k, 10)) <= 0 for k in range(-10, 11))
assert diff(f, x).subs(x, Rational(3, 2)) > 0 and diff(f, x).subs(x, -Rational(3, 2)) > 0

# ar-t-derivee-monotonie-ex2 : g = x - ln(1+x), g' = x/(1+x), g(0) = 0.
g = x - log(1 + x)
assert simplify(diff(g, x) - x / (1 + x)) == 0 and g.subs(x, 0) == 0
assert all(float(g.subs(x, v)) > 0 for v in (-0.9, -0.5, 0.1, 3, 100))

# ar-t-derivee-monotonie-ex3 : f = x - 1/x ; f(-1) = 0, f(1/2) = -3/2 ; x - 1/x = 3 a deux solutions de signes opposés.
f = x - 1 / x
assert simplify(diff(f, x) - (1 + 1 / x**2)) == 0
assert f.subs(x, -1) == 0 and f.subs(x, Rational(1, 2)) == -Rational(3, 2)
sol = solve(x**2 - 3 * x - 1, x)
assert set(sol) == {(3 + sqrt(13)) / 2, (3 - sqrt(13)) / 2}
assert len(sol) == 2 and sorted(sympy.sign(s_) for s_ in sol) == [-1, 1]
assert all(simplify(f.subs(x, s_) - 3) == 0 for s_ in sol)

# ar-t-derivee-monotonie-ex4 (et affirmation 4) : f = x + 2x^2 sin(1/x) ; f'(0) = 1 ; f'(1/(2k pi)) = -1.
f = x + 2 * x**2 * sin(1 / x)
fp = diff(f, x)
assert simplify(fp - (1 + 4 * x * sin(1 / x) - 2 * cos(1 / x))) == 0
assert limit((f - 0) / x, x, 0) == 1
for k in range(1, 6):
    assert simplify(fp.subs(x, 1 / (2 * k * pi))) == -1
# Les leurres 1, 3, 0 sont faux.
assert all(simplify(fp.subs(x, 1 / (2 * pi))) != c for c in (1, 3, 0))

# ---------------------------------------------------------------------------
# ar-l1-borne-sup
# ---------------------------------------------------------------------------

# Exemple : 1 - 1/n > 1 - eps dès que n > 1/eps.
for eps in (Rational(1, 2), Rational(1, 100), Rational(1, 10**5), Rational(3, 7)):
    m = math.floor(1 / eps) + 1
    assert 1 - Rational(1, m) > 1 - eps

# ar-l1-borne-sup-ex1 : (2n+1)/(n+2) = 2 - 3/(n+2), croissante, a_0 = 1/2, < 2.
assert simplify((2 * n + 1) / (n + 2) - (2 - 3 / (n + 2))) == 0
termes = [Rational(2 * k + 1, k + 2) for k in range(0, 200)]
assert termes[0] == Rational(1, 2) and all(a < b for a, b in zip(termes, termes[1:])) and all(a < 2 for a in termes)
assert limit((2 * np_ + 1) / (np_ + 2), np_, oo) == 2

# ar-l1-borne-sup-ex3 : sqrt 2 est irrationnel.
assert not sqrt(2).is_rational

# ---------------------------------------------------------------------------
# ar-l1-bolzano-weierstrass
# ---------------------------------------------------------------------------

# ar-l1-bolzano-weierstrass-ex1 : cos(n pi/3) selon n mod 6 : 1, 1/2, -1/2, -1, -1/2, 1/2.
valeurs = [cos(k * pi / 3) for k in range(6)]
assert valeurs == [1, Rational(1, 2), -Rational(1, 2), -1, -Rational(1, 2), Rational(1, 2)]
assert set(valeurs) == {-1, -Rational(1, 2), Rational(1, 2), 1}
for k in range(6, 60):
    assert cos(k * pi / 3) == valeurs[k % 6]

# ---------------------------------------------------------------------------
# ar-l1-bornes-atteintes
# ---------------------------------------------------------------------------

# Contre-exemple : x/(1+x) < 1 sur [0, +oo[ et tend vers 1.
assert limit(x / (1 + x), x, oo) == 1 and simplify(1 - x / (1 + x) - 1 / (1 + x)) == 0
# Affirmation 1 : arctan tend vers ±pi/2.
assert limit(atan(x), x, oo) == pi / 2 and limit(atan(x), x, -oo) == -pi / 2

# ar-l1-bornes-atteintes-ex1 : x^3 - 3x sur [0,2] : f(0) = 0, f(1) = -2, f(2) = 2.
f = x**3 - 3 * x
assert [f.subs(x, v) for v in (0, 1, 2)] == [0, -2, 2]
assert sympy.maximum(f, x, sympy.Interval(0, 2)) == 2 and sympy.minimum(f, x, sympy.Interval(0, 2)) == -2

# ar-l1-bornes-atteintes-ex2 : contre-exemple sans continuité : x^2 (x≠0), 1 en 0 : inf 0 non atteint.
assert sympy.minimum(x**2, x, sympy.Interval.Lopen(0, 1)) == 0  # borne inférieure 0, approchée mais 0 exclu

# ar-l1-bornes-atteintes-ex3 : f = x - x^2 : f/x = 1 - x, de max 1 - a sur [a,1], de borne sup 1 sur ]0,1].
assert simplify((x - x**2) / x - (1 - x)) == 0
assert sympy.maximum(x - x**2, x, sympy.Interval(0, 1)) == Rational(1, 4)

# ---------------------------------------------------------------------------
# ar-l1-rolle
# ---------------------------------------------------------------------------

# Contre-exemple : e^{it} : |f'(t)| = 1.
assert simplify(Abs(diff(sympy.exp(sympy.I * t), t))) == 1
# Affirmation 3 : sin sur [0, 3pi] : zéros de cos dans ]0, 3pi[ : pi/2, 3pi/2, 5pi/2.
assert sin(3 * pi) == 0
assert sorted(sympy.solveset(cos(x), x, sympy.Interval.open(0, 3 * pi))) == [pi / 2, 3 * pi / 2, 5 * pi / 2]

# ar-l1-rolle-ex1 : x^2 - 4x + 1 sur [0,4] : f(0) = f(4) = 1, f'(2) = 0 ; racines de f : 2 ± sqrt 3.
f = x**2 - 4 * x + 1
assert f.subs(x, 0) == 1 and f.subs(x, 4) == 1 and solve(diff(f, x), x) == [2]
assert set(solve(f, x)) == {2 - sqrt(3), 2 + sqrt(3)}

# ar-l1-rolle-ex2 : x^5 - 5x + 1 : P(-2) = -21, P(-1) = 5, P(1) = -3, P(2) = 23 ; P' a 2 racines réelles ; P en a 3.
P = x**5 - 5 * x + 1
assert [P.subs(x, v) for v in (-2, -1, 1, 2)] == [-21, 5, -3, 23]
assert sorted(real_roots(diff(P, x))) == [-1, 1]
assert len(real_roots(P)) == 3

# ar-l1-rolle-ex3 : g = e^{it}, g' = -sin t + i cos t ; Re g'(pi) = 0, Im g'(pi/2) = 0.
gp = diff(cos(t) + sympy.I * sin(t), t)
assert sympy.re(gp.subs(t, pi)) == 0 and sympy.im(gp.subs(t, pi / 2)) == 0
assert gp.subs(t, pi) != 0 and gp.subs(t, pi / 2) != 0

# ---------------------------------------------------------------------------
# ar-l1-taf
# ---------------------------------------------------------------------------

# Contre-exemple : |x| sur [-1,2] : pente moyenne 1/3.
assert (Abs(2) - Abs(-1)) / Rational(3) == Rational(1, 3)
# Contre-exemple : sin(ln x) : dérivée cos(ln x)/x ; valeurs 0 et 1 en e^{2n pi} et e^{2n pi + pi/2}.
assert simplify(diff(sin(log(x)), x) - cos(log(x)) / x) == 0
assert sin(log(sympy.exp(4 * pi))) == 0 and sin(log(sympy.exp(4 * pi + pi / 2))) == 1
# Exemple : (y-x)/y < ln y - ln x < (y-x)/x pour quelques 0 < x < y.
for (a, b) in ((1, 2), (0.5, 0.7), (3, 10)):
    assert (b - a) / b < math.log(b) - math.log(a) < (b - a) / a
# Affirmation 4 : f = sin(x^2)/x, f'(sqrt(2n pi)) = 2 et f'(sqrt((2n+1) pi)) = -2.
f = sin(x**2) / x
fp = diff(f, x)
assert simplify(fp - (2 * cos(x**2) - sin(x**2) / x**2)) == 0
for k in range(1, 4):
    assert simplify(fp.subs(x, sqrt(2 * k * pi))) == 2 and simplify(fp.subs(x, sqrt((2 * k + 1) * pi))) == -2

# ar-l1-taf-ex1 : x^3 sur [0,3] : 27 = 9c^2, c = sqrt 3 (et non 3/2, 3).
assert set(solve(9 * x**2 - 27, x)) == {sqrt(3), -sqrt(3)}
assert 0 < math.sqrt(3) < 3 and 9 * (1.5) ** 2 != 27 and 9 * 3**2 != 27

# ar-l1-taf-ex2 : encadrement de arctan(4/3) : seule la première proposition est vraie.
val = math.atan(4 / 3)
p4 = math.pi / 4
assert proche(Rational(1, 3) * Rational(9, 25), Rational(3, 25)) and proche(Rational(1, 3) * Rational(1, 2), Rational(1, 6))
assert p4 + 3 / 25 < val < p4 + 1 / 6  # bonne réponse
assert not (p4 + 1 / 6 < val < p4 + 1 / 3)  # leurre 1
assert not (p4 + 1 / 25 < val < p4 + 3 / 25)  # leurre 2
assert not (p4 + 4 / 25 < val < p4 + 1 / 3)  # leurre 3

# ar-l1-taf-ex3 : x^2 sin(1/x) : taux x sin(1/x) -> 0, dérivée 2x sin(1/x) - cos(1/x).
assert limit(x * sin(1 / x), x, 0) == 0
assert simplify(diff(x**2 * sin(1 / x), x) - (2 * x * sin(1 / x) - cos(1 / x))) == 0

# ar-l1-taf-ex4 : sin(x^2) : dérivée 2x cos(x^2), qui vaut 2 sqrt(2n pi) en sqrt(2n pi).
assert diff(sin(x**2), x) == 2 * x * cos(x**2)
assert simplify(diff(sin(x**2), x).subs(x, sqrt(2 * 3 * pi)) - 2 * sqrt(6 * pi)) == 0

# ===========================================================================
# Variantes (enrichissement) et exercice ar-t-derivee-monotonie-ex5
# ===========================================================================


def faux_partout(valeur, leurres, tol=1e-6):
    """Chaque leurre numérique diffère de la bonne valeur."""
    return all(abs(float(valeur) - float(c)) > tol for c in leurres)


# ar-t-limite-monotone-ex1 variante 1 : somme 1/(3^k+2) croissante, < 3/2 (géométrique), et même <= 5/6 ; limite ≈ 0,68.
assert sum(Rational(1, 3**k) for k in range(0, 40)) < Rational(3, 2)
assert simplify((1 - Rational(1, 3) ** (np_ + 1)) / (1 - Rational(1, 3)) - Rational(3, 2) * (1 - Rational(1, 3) ** (np_ + 1))) == 0
u = Rational(0)
for k in range(0, 60):
    u += Rational(1, 3**k + 2)
    assert u < Rational(3, 2) and u <= Rational(5, 6)
assert Rational(1, 3) + Rational(1, 2) == Rational(5, 6) and Rational(3, 2) - Rational(2, 3) == Rational(5, 6)
assert 0.67 < float(u) < 0.69  # donc ≠ 3/2 (leurre) et pas >= 3/2 (leurre)
# ar-t-limite-monotone-ex1 variante 2 : somme 1/(4^k+1) < 4/3, <= 4/3 - 1/2 = 5/6 ; limite ≈ 0,78.
u = Rational(0)
for k in range(0, 60):
    u += Rational(1, 4**k + 1)
    assert u < Rational(4, 3) and u <= Rational(5, 6)
assert Rational(4, 3) - Rational(1, 2) == Rational(5, 6)
assert 0.77 < float(u) < 0.79

# ar-t-limite-monotone-ex2 variante 1 : a_n = somme 1/sqrt(k+1) - 2 sqrt n, décroissante, > -2 sqrt 2, limite zeta(1/2) - 1 ≈ -2,46.
S, a_prec = 0.0, None
for k in range(1, 200001):
    S += 1 / math.sqrt(k + 1)
    a = S - 2 * math.sqrt(k)
    assert a > -2 * math.sqrt(2)
    if a_prec is not None:
        assert a < a_prec
    a_prec = a
assert abs(a_prec - (zeta_demi - 1)) < 0.01 and abs(a_prec - (-2 * math.sqrt(2))) > 0.3  # leurre « limite -2 sqrt 2 » faux
for k in range(1, 1001):
    assert 1 / math.sqrt(k + 1) - 2 * (math.sqrt(k + 2) - math.sqrt(k + 1)) > 0
    assert math.sqrt(k + 1) + math.sqrt(k) < 2 * math.sqrt(k + 2)
# Télescopage : somme_{k=1}^n (sqrt(k+2) - sqrt(k+1)) = sqrt(n+2) - sqrt 2.
for m in (1, 2, 7, 30):
    assert simplify(sum(sqrt(k + 2) - sqrt(k + 1) for k in range(1, m + 1)) - (sqrt(m + 2) - sqrt(2))) == 0
# ar-t-limite-monotone-ex2 variante 2 : a_n = somme 1/sqrt k - 2 sqrt(n+1), croissante, < 0, limite zeta(1/2) ≈ -1,46.
S, a_prec = 0.0, None
for k in range(1, 200001):
    S += 1 / math.sqrt(k)
    a = S - 2 * math.sqrt(k + 1)
    assert a < 0
    if a_prec is not None:
        assert a > a_prec
    a_prec = a
assert abs(a_prec - zeta_demi) < 0.01 and abs(a_prec) > 1  # leurre « limite 0 » faux
for k in range(1, 1001):
    assert 1 / math.sqrt(k) < 2 * (math.sqrt(k) - math.sqrt(k - 1))
for m in (1, 2, 7, 30):
    assert simplify(sum(sqrt(k) - sqrt(k - 1) for k in range(1, m + 1)) - sqrt(m)) == 0

# ar-t-gendarmes-ex1 variante 1 : (5n ± 1)/(2n+3) -> 5/2 (leurres 5/3, 5 faux).
assert limit((5 * np_ - 1) / (2 * np_ + 3), np_, oo) == Rational(5, 2) == limit((5 * np_ + 1) / (2 * np_ + 3), np_, oo)
assert faux_partout(Rational(5, 2), (Rational(5, 3), 5))
# ar-t-gendarmes-ex1 variante 2 : (4n^2 ± n)/(n^2+1) -> 4 (leurres 5, 3 faux).
assert limit((4 * np_**2 - np_) / (np_**2 + 1), np_, oo) == 4 == limit((4 * np_**2 + np_) / (np_**2 + 1), np_, oo)
assert proche((4 * 10**16 + 10**8 * math.cos(10**8)) / (10**16 + 1), 4, 1e-6)

# ar-t-gendarmes-ex2 variante 1 : (4^n + 5^n)^(1/n) -> 5 ; leurres 9, 2 sqrt 5 = sqrt 20, 1.
assert limit(5 * 2 ** (1 / np_), np_, oo) == 5
for k in (1, 5, 50, 300):
    v = (4**k + 5**k) ** (1 / k)
    assert 5 <= v <= 5 * 2 ** (1 / k) + 1e-9
assert simplify(2 * sqrt(5) - sqrt(20)) == 0 and faux_partout(5, (9, 2 * math.sqrt(5), 1), 0.4)
# ar-t-gendarmes-ex2 variante 2 : (3^n + 4^n + 6^n)^(1/n) -> 6 ; leurres 13, 2·9^(1/3) = 72^(1/3), 1.
assert limit(6 * 3 ** (1 / np_), np_, oo) == 6
for k in (1, 5, 50, 300):
    v = (3**k + 4**k + 6**k) ** (1 / k)
    assert 6 <= v <= 6 * 3 ** (1 / k) + 1e-9
assert simplify(2 * sympy.cbrt(9) - sympy.cbrt(72)) == 0 and faux_partout(6, (13, 2 * 9 ** (1 / 3), 1), 0.4)

# ar-t-tvi-ex1 variante 1 : x^3 + 2x - 2 : f(1/2) = -7/8, f(1) = 1, une seule racine (≈ 0,77) dans ]1/2, 1[.
f = x**3 + 2 * x - 2
assert f.subs(x, Rational(1, 2)) == -Rational(7, 8) and f.subs(x, 1) == 1 and diff(f, x) == 3 * x**2 + 2
rr = real_roots(f)
assert len(rr) == 1 and 0.5 < float(rr[0]) < 1  # pas dans [0,1/2], [1,2], [-1,0]
# ar-t-tvi-ex1 variante 2 : x^5 + x + 1 : f(-1) = -1, f(-1/2) = 15/32, une seule racine (≈ -0,75) dans ]-1, -1/2[.
f = x**5 + x + 1
assert f.subs(x, -1) == -1 and f.subs(x, -Rational(1, 2)) == Rational(15, 32) and diff(f, x) == 5 * x**4 + 1
rr = real_roots(f)
assert len(rr) == 1 and -1 < float(rr[0]) < -0.5  # pas dans [-1/2,0], [0,1], [-2,-1]

# ar-t-tvi-ex3 variante 1 : x^3 - 4x + 2 : valeurs -13, 2, 2, -1, 2 en -3, -2, 0, 1, 2 ; 3 racines.
P = x**3 - 4 * x + 2
assert [P.subs(x, v) for v in (-3, -2, 0, 1, 2)] == [-13, 2, 2, -1, 2]
rr = [float(z) for z in real_roots(P)]
assert len(rr) == 3 and -3 < rr[0] < -2 and 0 < rr[1] < 1 and 1 < rr[2] < 2
# ar-t-tvi-ex3 variante 2 : x^3 - 3x + 3 : P(-1) = 5, P(1) = 1 > 0, P(-3) = -15, P(-2) = 1 ; une seule racine, dans ]-3,-2[.
P = x**3 - 3 * x + 3
assert sympy.expand(diff(P, x) - 3 * (x - 1) * (x + 1)) == 0
assert [P.subs(x, v) for v in (-3, -2, -1, 1)] == [-15, 1, 5, 1]
rr = [float(z) for z in real_roots(P)]
assert len(rr) == 1 and -3 < rr[0] < -2
assert sympy.minimum(P, x, sympy.Interval(-1, oo)) == 1

# ar-t-point-fixe-ex1 variante 1 : l/2 + 3 = l donne 6 (leurres 3, 2, 1 faux).
assert solve(x / 2 + 3 - x, x) == [6]
u = 1.0
for _ in range(80):
    u = u / 2 + 3
assert proche(u, 6)
# ar-t-point-fixe-ex1 variante 2 : (3l - 2)/4 = l donne -2 (leurres 2, -1/2, -2/7 faux).
assert solve((3 * x - 2) / 4 - x, x) == [-2]
u = 2.0
for _ in range(200):
    u = (3 * u - 2) / 4
assert proche(u, -2, 1e-9) and faux_partout(-2, (2, -0.5, -2 / 7))

# ar-t-point-fixe-ex2 variante 1 : sqrt(u+2), u_0 = 0 : [0,2] stable, (2-u)(u+1), limite 2 (pas -1, sqrt 2).
assert sorted(solve(x**2 - x - 2, x)) == [-1, 2]
assert sympy.expand((2 - x) * (x + 1) - (x + 2 - x**2)) == 0
u = 0.0
for _ in range(80):
    u_suiv = math.sqrt(u + 2)
    assert 0 <= u <= u_suiv <= 2
    u = u_suiv
assert proche(u, 2) and proche(math.sqrt(0 + 2), math.sqrt(2))
# ar-t-point-fixe-ex2 variante 2 : sqrt(3u+4), u_0 = 1 : [1,4] stable (7 <= 3u+4 <= 16), (4-u)(u+1), limite 4.
assert sorted(solve(x**2 - 3 * x - 4, x)) == [-1, 4]
assert sympy.expand((4 - x) * (x + 1) - (3 * x + 4 - x**2)) == 0
assert 3 * 1 + 4 == 7 and 3 * 4 + 4 == 16
u = 1.0
for _ in range(80):
    u_suiv = math.sqrt(3 * u + 4)
    assert 1 <= u <= u_suiv <= 4
    u = u_suiv
assert proche(u, 4)

# ar-t-point-fixe-ex3 variante 1 : f = x + x^2 : seul point fixe 0 ; f(]-1,0[) ⊂ [-1/4, 0[ ; v_0 = 1 : v_1 = 2, v_n >= 2 + 4(n-1).
assert solve(x + x**2 - x, x) == [0]
assert sympy.minimum(x + x**2, x, sympy.Interval(-1, 0)) == -Rational(1, 4)
assert sympy.maximum(x + x**2, x, sympy.Interval(-1, 0)) == 0
v = 1
for k in range(1, 7):
    v = v + v**2
    assert v >= 2 + 4 * (k - 1)
    if k == 1:
        assert v == 2
u = -0.5
for _ in range(10**5):
    u = u + u * u
assert -1e-4 < u < 0
# ar-t-point-fixe-ex3 variante 2 : f = x - 2x^2 : seul point fixe 0 ; f([0,1/2]) = [0,1/8] ; v_0 = 1 : v_1 = -1, v_n <= -1 - 2(n-1).
assert solve(x - 2 * x**2 - x, x) == [0]
assert sympy.maximum(x - 2 * x**2, x, sympy.Interval(0, Rational(1, 2))) == Rational(1, 8)
assert sympy.minimum(x - 2 * x**2, x, sympy.Interval(0, Rational(1, 2))) == 0
assert sympy.expand(Rational(1, 8) - 2 * (x - Rational(1, 4)) ** 2 - (x - 2 * x**2)) == 0
v = 1
for k in range(1, 7):
    v = v - 2 * v**2
    assert v <= -1 - 2 * (k - 1)
    if k == 1:
        assert v == -1
u = 0.25
for _ in range(10**5):
    u = u - 2 * u * u
assert 0 < u < 1e-4

# ar-t-point-fixe-ex4 variante 1 : f = 2/(1+x) sur [0,2] : f o f = 2(1+x)/(3+x), points fixes 1 et -2 ; u_1..u_3 = 2, 2/3, 6/5.
F = 2 / (1 + x)
G = simplify(F.subs(x, F))
assert simplify(G - 2 * (1 + x) / (3 + x)) == 0
assert set(solve(G - x, x)) == {1, -2}
assert F.subs(x, 0) == 2 and F.subs(x, 2) == Rational(2, 3)
assert sympy.diff(G, x).subs(x, 1) > 0 and simplify(sympy.diff(G, x) - 4 / (3 + x) ** 2) == 0
u = [Rational(0)]
for _ in range(3):
    u.append(2 / (1 + u[-1]))
assert u[1:] == [2, Rational(2, 3), Rational(6, 5)]
u = [0.0]
for _ in range(80):
    u.append(2 / (1 + u[-1]))
assert all(u[2 * k] <= u[2 * k + 2] for k in range(35)) and all(u[2 * k + 1] >= u[2 * k + 3] for k in range(35))
assert proche(u[-1], 1, 1e-9) and faux_partout(1, (-2, 2, 2 / 3))
# ar-t-point-fixe-ex4 variante 2 : f = 1/(2+x) sur [0,1/2] : f o f = (2+x)/(5+2x), point fixe sqrt 2 - 1 ; u_1..u_3 = 1/2, 2/5, 5/12.
F = 1 / (2 + x)
G = simplify(F.subs(x, F))
assert simplify(G - (2 + x) / (5 + 2 * x)) == 0
assert set(solve(G - x, x)) == {-1 + sqrt(2), -1 - sqrt(2)}
assert F.subs(x, 0) == Rational(1, 2) and F.subs(x, Rational(1, 2)) == Rational(2, 5)
assert simplify(sympy.diff(G, x) - 1 / (5 + 2 * x) ** 2) == 0
u = [Rational(0)]
for _ in range(3):
    u.append(1 / (2 + u[-1]))
assert u[1:] == [Rational(1, 2), Rational(2, 5), Rational(5, 12)]
u = [0.0]
for _ in range(80):
    u.append(1 / (2 + u[-1]))
assert all(u[2 * k] <= u[2 * k + 2] for k in range(35)) and all(u[2 * k + 1] >= u[2 * k + 3] for k in range(35))
assert proche(u[-1], math.sqrt(2) - 1, 1e-12) and faux_partout(math.sqrt(2) - 1, (math.sqrt(2) + 1, -1 - math.sqrt(2), 0.5), 1e-3)

# ar-t-derivee-monotonie-ex1 variante 1 : x^3 - 12x : f' = 3(x-2)(x+2), décroissante sur [-2,2] ; zéros de f : 0, ±2 sqrt 3.
f = x**3 - 12 * x
assert sympy.expand(diff(f, x) - 3 * (x - 2) * (x + 2)) == 0
assert set(solve(f, x)) == {0, 2 * sqrt(3), -2 * sqrt(3)}
assert sympy.solveset(diff(f, x) <= 0, x, sympy.S.Reals) == sympy.Interval(-2, 2)
# ar-t-derivee-monotonie-ex1 variante 2 : 2x^3 - 3x^2 - 12x : f' = 6(x-2)(x+1), décroissante sur [-1,2] (pas [-2,1]).
f = 2 * x**3 - 3 * x**2 - 12 * x
assert sympy.expand(diff(f, x) - 6 * (x - 2) * (x + 1)) == 0
assert sympy.solveset(diff(f, x) <= 0, x, sympy.S.Reals) == sympy.Interval(-1, 2)

# ar-t-derivee-monotonie-ex3 variante 1 : x - 4/x : f(-2) = 0 > f(1) = -3 ; f(x) = 1 : (1 ± sqrt 17)/2, de signes opposés.
f = x - 4 / x
assert simplify(diff(f, x) - (1 + 4 / x**2)) == 0
assert f.subs(x, -2) == 0 and f.subs(x, 1) == -3
sol = solve(x**2 - x - 4, x)
assert set(sol) == {(1 + sqrt(17)) / 2, (1 - sqrt(17)) / 2} and sorted(sympy.sign(s_) for s_ in sol) == [-1, 1]
assert all(simplify(f.subs(x, s_) - 1) == 0 for s_ in sol)
# ar-t-derivee-monotonie-ex3 variante 2 : 2x - 1/x : f(-1) = -1 > f(1/4) = -7/2 ; f(x) = -1 : x = -1 ou 1/2.
f = 2 * x - 1 / x
assert simplify(diff(f, x) - (2 + 1 / x**2)) == 0
assert f.subs(x, -1) == -1 and f.subs(x, Rational(1, 4)) == -Rational(7, 2)
assert set(solve(f + 1, x)) == {-1, Rational(1, 2)}

# ar-t-derivee-monotonie-ex4 variante 1 : x + 3x^2 sin(1/x) : f'(0) = 1, f'(1/(2k pi)) = -2 (leurres 1, 4, 0 faux).
f = x + 3 * x**2 * sin(1 / x)
fp = diff(f, x)
assert simplify(fp - (1 + 6 * x * sin(1 / x) - 3 * cos(1 / x))) == 0
assert limit(f / x, x, 0) == 1
for k in range(1, 6):
    assert simplify(fp.subs(x, 1 / (2 * k * pi))) == -2
# ar-t-derivee-monotonie-ex4 variante 2 : 2x + 3x^2 cos(1/x) : f'(0) = 2, f'(1/(2k pi - pi/2)) = -1 (leurres 2, 5, 0 faux).
f = 2 * x + 3 * x**2 * cos(1 / x)
fp = diff(f, x)
assert simplify(fp - (2 + 6 * x * cos(1 / x) + 3 * sin(1 / x))) == 0
assert limit(f / x, x, 0) == 2
for k in range(1, 6):
    xk = 1 / (2 * k * pi - pi / 2)
    assert xk.is_positive and simplify(fp.subs(x, xk)) == -1

# ar-t-derivee-monotonie-ex5 : e^x - 3x : minimum 3 - 3 ln 3 < 0 en ln 3, limites +oo ; exactement 2 solutions (≈ 0,62 et 1,51).
f = sympy.exp(x) - 3 * x
assert solve(diff(f, x), x) == [log(3)]
assert simplify(f.subs(x, log(3)) - (3 - 3 * log(3))) == 0 and float(3 - 3 * log(3)) < 0 and log(3) > 1
assert limit(f, x, oo) == oo and limit(f, x, -oo) == oo
r1, r2 = sympy.nsolve(f, x, 0.6), sympy.nsolve(f, x, 1.5)
assert 0.61 < r1 < 0.63 and 1.50 < r2 < 1.52
# Comptage numérique des changements de signe sur une grille fine (f a au plus 2 zéros car strictement monotone de part et d'autre de ln 3).
grille = [i / 1000 for i in range(-5000, 5001)]
changements = sum(1 for a, b in zip(grille, grille[1:]) if (math.exp(a) - 3 * a) * (math.exp(b) - 3 * b) < 0)
assert changements == 2
# ar-t-derivee-monotonie-ex5 variante 1 : e^x - 2x : minimum 2 - 2 ln 2 > 0 : aucune solution.
f = sympy.exp(x) - 2 * x
assert solve(diff(f, x), x) == [log(2)]
assert simplify(f.subs(x, log(2)) - (2 - 2 * log(2))) == 0 and float(2 - 2 * log(2)) > 0 and log(2) < 1
assert all(math.exp(a) - 2 * a > 0 for a in grille)
# ar-t-derivee-monotonie-ex5 variante 2 : e^x - e x : minimum f(1) = 0 : une seule solution, x = 1.
f = sympy.exp(x) - sympy.E * x
assert solve(diff(f, x), x) == [1] and f.subs(x, 1) == 0
assert all(math.exp(a) - math.e * a > 0 for a in grille if abs(a - 1) > 1e-9)

# ar-l1-borne-sup-ex1 variante 1 : (3n+2)/(n+1) = 3 - 1/(n+1), croissante, min 2, sup 3 non atteint.
assert simplify((3 * n + 2) / (n + 1) - (3 - 1 / (n + 1))) == 0
termes = [Rational(3 * k + 2, k + 1) for k in range(0, 200)]
assert termes[0] == 2 and all(a < b for a, b in zip(termes, termes[1:])) and all(a < 3 for a in termes)
assert limit((3 * np_ + 2) / (np_ + 1), np_, oo) == 3
# ar-l1-borne-sup-ex1 variante 2 : (n+4)/(2n+1) = 1/2 + 7/(2(2n+1)), décroissante, max 4, inf 1/2 non atteint.
assert simplify((n + 4) / (2 * n + 1) - (Rational(1, 2) + 7 / (2 * (2 * n + 1)))) == 0
termes = [Rational(k + 4, 2 * k + 1) for k in range(0, 200)]
assert termes[0] == 4 and all(a > b for a, b in zip(termes, termes[1:])) and all(a > Rational(1, 2) for a in termes)
assert limit((np_ + 4) / (2 * np_ + 1), np_, oo) == Rational(1, 2)

# ar-l1-bolzano-weierstrass-ex1 variante 1 : sin(n pi/2) selon n mod 4 : 0, 1, 0, -1 ; ensemble {-1, 0, 1}.
valeurs = [sin(k * pi / 2) for k in range(4)]
assert valeurs == [0, 1, 0, -1] and set(valeurs) == {-1, 0, 1}
assert all(sin(k * pi / 2) == valeurs[k % 4] for k in range(4, 60))
# ar-l1-bolzano-weierstrass-ex1 variante 2 : cos(2n pi/3) selon n mod 3 : 1, -1/2, -1/2 ; ensemble {-1/2, 1}.
valeurs = [cos(2 * k * pi / 3) for k in range(3)]
assert valeurs == [1, -Rational(1, 2), -Rational(1, 2)] and set(valeurs) == {-Rational(1, 2), 1}
assert all(cos(2 * k * pi / 3) == valeurs[k % 3] for k in range(3, 60))

# ar-l1-bornes-atteintes-ex1 variante 1 : x^3 - 12x sur [-1,3] : f(-1) = 11, f(2) = -16, f(3) = -9 ; max 11, min -16.
f = x**3 - 12 * x
assert [f.subs(x, v) for v in (-1, 2, 3, -2)] == [11, -16, -9, 16]
assert sympy.maximum(f, x, sympy.Interval(-1, 3)) == 11 and sympy.minimum(f, x, sympy.Interval(-1, 3)) == -16
# ar-l1-bornes-atteintes-ex1 variante 2 : x^4 - 2x^2 sur [-1/2, 2] : valeurs -7/16, 0, -1, 8 ; max 8, min -1.
f = x**4 - 2 * x**2
assert [f.subs(x, v) for v in (-Rational(1, 2), 0, 1, 2)] == [-Rational(7, 16), 0, -1, 8]
assert sympy.maximum(f, x, sympy.Interval(-Rational(1, 2), 2)) == 8
assert sympy.minimum(f, x, sympy.Interval(-Rational(1, 2), 2)) == -1

# ar-l1-rolle-ex1 variante 1 : x^3 - 3x sur [0, sqrt 3] : f(0) = f(sqrt 3) = 0, f' nul en ±1, seul 1 dans ]0, sqrt 3[.
f = x**3 - 3 * x
assert f.subs(x, 0) == 0 and simplify(f.subs(x, sqrt(3))) == 0
assert set(solve(diff(f, x), x)) == {-1, 1} and 0 < 1 < math.sqrt(3)
# ar-l1-rolle-ex1 variante 2 : x^3 - 6x^2 + 9x + 1 sur [0,3] : f(0) = f(3) = 1, f' = 3(x-1)(x-3), seul 1 dans ]0,3[.
f = x**3 - 6 * x**2 + 9 * x + 1
assert f.subs(x, 0) == 1 and f.subs(x, 3) == 1
assert sympy.expand(diff(f, x) - 3 * (x - 1) * (x - 3)) == 0 and set(solve(diff(f, x), x)) == {1, 3}

# ar-l1-rolle-ex2 variante 1 : x^4 - 4x + 1 : P' a une seule racine réelle (1) ; P(0) = 1, P(1) = -2, P(2) = 9 ; 2 racines.
P = x**4 - 4 * x + 1
assert real_roots(diff(P, x)) == [1] and sympy.expand(diff(P, x) - 4 * (x - 1) * (x**2 + x + 1)) == 0
assert [P.subs(x, v) for v in (0, 1, 2)] == [1, -2, 9]
rr = [float(z) for z in real_roots(P)]
assert len(rr) == 2 and 0 < rr[0] < 1 < rr[1] < 2
# ar-l1-rolle-ex2 variante 2 : x^5 - 10x + 3 : P' a 2 racines réelles (±2^(1/4)) ; P(-2) = -9, P(-1) = 12, P(1) = -6, P(2) = 15 ; 3 racines.
P = x**5 - 10 * x + 3
assert len(real_roots(diff(P, x))) == 2
assert sympy.expand(diff(P, x) - 5 * (x**2 - sqrt(2)) * (x**2 + sqrt(2))) == 0
assert set(solve(x**4 - 2, x)) & {r_ for r_ in solve(x**4 - 2, x) if r_.is_real} == {2 ** Rational(1, 4), -(2 ** Rational(1, 4))}
assert [P.subs(x, v) for v in (-2, -1, 1, 2)] == [-9, 12, -6, 15]
assert len(real_roots(P)) == 3

# ar-l1-taf-ex1 variante 1 : x^3 sur [0,2] : 8 = 6c^2, c = 2 sqrt 3/3 (≈ 1,15) ; leurres 1, 4/3 faux.
assert set(solve(6 * x**2 - 8, x)) == {2 * sqrt(3) / 3, -2 * sqrt(3) / 3}
assert 0 < float(2 * sqrt(3) / 3) < 2 and 6 * 1**2 != 8 and 6 * Rational(4, 3) ** 2 != 8
# ar-l1-taf-ex1 variante 2 : 1/x sur [1,4] : -3/4 = -3/c^2, c = 2 ; leurres 5/2, 4 faux.
assert set(solve(-3 / x**2 - (Rational(1, 4) - 1), x)) == {2, -2}
assert -3 / Rational(5, 2) ** 2 != -Rational(3, 4) and -3 / Rational(4) ** 2 != -Rational(3, 4)

# ar-l1-taf-ex2 variante 1 : arctan 2 : 1/5 < arctan 2 - pi/4 < 1/2, seule la première proposition est vraie.
val = math.atan(2) - p4
assert 1 / 5 < val < 1 / 2
assert not (1 / 2 < val < 1) and not (1 / 10 < val < 1 / 5) and not (2 / 5 < val < 1)
# ar-l1-taf-ex2 variante 2 : arctan(3/2) : (1/2)(4/13) = 2/13 < arctan(3/2) - pi/4 < 1/4.
assert Rational(1, 2) * Rational(4, 13) == Rational(2, 13) and 1 + Rational(3, 2) ** 2 == Rational(13, 4)
val = math.atan(1.5) - p4
assert 2 / 13 < val < 1 / 4
assert not (1 / 4 < val < 1 / 2) and not (1 / 13 < val < 2 / 13) and not (3 / 13 < val < 1 / 2)

print("analyse_reelle.py : toutes les vérifications passent.")
