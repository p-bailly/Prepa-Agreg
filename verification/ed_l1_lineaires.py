"""Vérifications calculatoires du chapitre « ed-l1-lineaires » (équations différentielles linéaires, L1).

Cartes et exercices de content/{cartes,exercices}/equations-differentielles/l1-ed-l1-lineaires.json.
Chaque assert cite l'identifiant de la carte ou de l'exercice concerné (v1, v2 = variantes).
Lancer : python3 verification/ed_l1_lineaires.py (doit finir sans erreur).
"""

from sympy import (E, I, Matrix, N, Rational, Sum, atan, cos, diff, exp, im, integrate, limit, log, oo, pi, re,
                   simplify, sin, sqrt, symbols, tan)

x, t, h = symbols("x t h", real=True)
lam, mu, C, k0 = symbols("lam mu C k0", real=True)
POINTS = [Rational(1, 3), Rational(7, 10), Rational(13, 10), Rational(-1, 2), Rational(-6, 5)]


def nul(e, points=POINTS):
    """e (expression en x) est identiquement nulle : simplification, sinon contrôle numérique en plusieurs points."""
    if simplify(e) == 0:
        return True
    return all(abs(N(e.subs(x, p), 40)) < 1e-30 for p in points)


def L1(y, a):
    return diff(y, x) + a * y


def L2(y, a, b, c):
    return a * diff(y, x, 2) + b * diff(y, x) + c * y


def sol1(y, a, b, x0=None, y0=None, points=POINTS):
    ok = nul(L1(y, a) - b, points)
    if x0 is not None:
        ok = ok and nul(y.subs(x, x0) - y0, [0])
    return ok


def sol2(y, a, b, c, f, cond=(), points=POINTS):
    """cond : liste de (ordre de dérivation, point, valeur)."""
    ok = nul(L2(y, a, b, c) - f, points)
    for o, p, v in cond:
        ok = ok and simplify(diff(y, x, o).subs(x, p) - v) == 0
    return ok


POS = [Rational(1, 3), Rational(7, 10), Rational(13, 10), 2, 3]  # points de ]0, +oo[
TRIG = [Rational(1, 3), Rational(7, 10), Rational(-1, 2), Rational(-6, 5), Rational(1, 10)]  # points de ]-pi/2, pi/2[

# =============================================================================
# ed-l1-lineaires-ordre1 (carte)
# =============================================================================
# Trou 3 (a(x) = 2x, A = x^2) : C e^{-A} solution ; leurres C e^{A}, C e^{-x a(x)}, e^{-A} + C non solutions (C = 1).
assert sol1(C * exp(-x**2), 2 * x, 0)
assert not sol1(exp(x**2), 2 * x, 0)
assert not sol1(exp(-x * 2 * x), 2 * x, 0)
assert not sol1(exp(-x**2) + 1, 2 * x, 0)
# Trou 2 / piège « bornées » : y' = b avec b = 1_{[0,+oo[} : y = x + c à droite, c à gauche ; taux en 0 : 1 à droite, 0 à gauche.
assert limit(((h + 5) - 5) / h, h, 0, "+") == 1 and limit((5 - 5) / h, h, 0, "-") == 0
# Démo : y_p = e^{-A} * int_{x0}^x b e^{A} est solution (a = 2x, b = 2x, x0 = 0).
G = integrate(2 * t * exp(t**2), (t, 0, x))
assert sol1(G * exp(-x**2), 2 * x, 2 * x)
# Exemple : 1 + 2e^{-x^2} solution de y' + 2xy = 2x, y(0) = 3 ; C x solution de y' - y/x = 0.
assert sol1(1 + 2 * exp(-x**2), 2 * x, 2 * x, 0, 3)
assert sol1(C * x, -1 / x, 0, points=POS)
# Contre-exemple (intervalle) : C/x solution de y' + y/x = 0 (C = 1, 2), sur R+* et R-*.
assert sol1(1 / x, 1 / x, 0) and sol1(2 / x, 1 / x, 0)
# Contre-exemple (coefficient de y' nul en 0) / affirmation 2 : x^2 et 0 vérifient x y' - 2y = 0 ; dérivée en 0 du recollement nulle.
assert nul(x * diff(x**2, x) - 2 * x**2)
assert limit(h**2 / h, h, 0, "+") == 0
# Affirmation 3 : 2 vérifie y' + y = 2, pas 1.
assert sol1(2 + 0 * x, 1, 2) and not sol1(2 + 0 * x, 1, 1)

# ed-l1-lineaires-ordre1-ex1 (et v1, v2) : bonne proposition, leurres faux.
for a, y0, bon, leurres in [
    (2 * x, 3, 3 * exp(-x**2), [3 * exp(x**2), 3 * exp(-2 * x), exp(-x**2) + 2]),
    (-3 * x**2, -1, -exp(x**3), [-exp(-x**3), -exp(3 * x**2), exp(x**3) - 2]),
    (cos(x), 2, 2 * exp(-sin(x)), [2 * exp(sin(x)), 2 * exp(-cos(x)), exp(-sin(x)) + 1]),
]:
    assert sol1(bon, a, 0, 0, y0)
    for l in leurres:
        assert not sol1(l, a, 0, 0, y0)

# ed-l1-lineaires-ordre1-ex2 (et v1 sur ]0,+oo[, v2 sur R) : équations normalisées.
for a, x0, y0, bon, leurres, pts in [
    (2 / x, 1, 5, 5 / x**2, [5 * x**2, 5 / x, 5 * exp(2 - 2 * x)], POS),
    (-3 / x, 1, 2, 2 * x**3, [2 / x**3, 2 * x**2, 2 * exp(3 * x - 3)], POS),
    (2 * x / (1 + x**2), 0, 4, 4 / (1 + x**2), [4 * (1 + x**2), 4 * exp(-x**2), 4 / (1 + x)**2], POS),
]:
    assert sol1(bon, a, 0, x0, y0, pts)
    for l in leurres:
        assert not sol1(l, a, 0, x0, y0, pts)
# ex2 : forme non normalisée x y' + 2y = 0 vérifiée directement ; v2 : (1+x^2) y' + 2x y = 0.
assert nul(x * diff(5 / x**2, x) + 2 * 5 / x**2)
assert nul((1 + x**2) * diff(4 / (1 + x**2), x) + 2 * x * 4 / (1 + x**2))

# ed-l1-lineaires-ordre1-ex3 (et v1, v2) : constante particulière, solution avec la condition.
for a, b, k, y0, bon, leurres in [
    (2 * x, 4 * x, 2, 0, 2 - 2 * exp(-x**2), [2 * exp(-x**2) - 2, 2 - 2 * exp(x**2), 2 * x**2 * exp(-x**2)]),
    (-3 * x**2, 6 * x**2, -2, 0, 2 * exp(x**3) - 2, [2 - 2 * exp(x**3), 2 * exp(-x**3) - 2, 2 * exp(-x**3) - 2 * exp(x**3)]),
    (cos(x), cos(x), 1, 3, 1 + 2 * exp(-sin(x)), [3 * exp(-sin(x)), 1 + 2 * exp(sin(x)), 4 * exp(-sin(x)) - 1]),
]:
    assert sol1(k + 0 * x, a, b)  # solution constante
    assert sol1(bon, a, b, 0, y0)
    for l in leurres:
        assert not sol1(l, a, b, 0, y0)

# ed-l1-lineaires-ordre1-ex4 : x y' - 2y = 0 ; recollement λx² / μx² dérivable en 0 (taux λh → 0).
assert nul(x * diff(lam * x**2, x) - 2 * lam * x**2)
assert limit(lam * h**2 / h, h, 0) == 0
# Leurre « C x² » : la fonction x² (x ≥ 0), 0 (x < 0) n'est pas de la forme C x² (valeurs en 1 et -1 : 1 et 0).
# v1 : x y' - y = 0 : C x solution ; recollement λx/μx : taux λ à droite, μ à gauche ; |x| non dérivable en 0.
assert nul(x * diff(C * x, x) - C * x)
assert limit(abs(h) / h, h, 0, "+") == 1 and limit(abs(h) / h, h, 0, "-") == -1
# v2 : x y' + y = 0 : C/x solution sur chaque demi-droite, non bornée en 0+ ; leurre C x non solution.
assert nul(x * diff(C / x, x) + C / x)
assert limit(1 / h, h, 0, "+") == oo
assert not nul(x * diff(x, x) + x)

# ed-l1-lineaires-ordre1-ex5 : sin : m = 0, solutions C e^{cos x} 2π-périodiques ; 1 + sin : m = 2π.
assert integrate(sin(t), (t, 0, 2 * pi)) == 0 and integrate(1 + sin(t), (t, 0, 2 * pi)) == 2 * pi
assert sol1(exp(cos(x)), sin(x), 0) and simplify(exp(cos(x + 2 * pi)) - exp(cos(x))) == 0
y5 = exp(-x + cos(x))
assert sol1(y5, 1 + sin(x), 0) and simplify(y5.subs(x, x + 2 * pi) - exp(-2 * pi) * y5) == 0

# =============================================================================
# ed-l1-lineaires-variation-constante (carte)
# =============================================================================
# Démo : y = C(x) e^{-A} donne y' + a y = C' e^{-A} (a = 2x, A = x², C quelconque).
Cf = symbols("Cf", cls=__import__("sympy").Function)
assert simplify(L1(Cf(x) * exp(-x**2), 2 * x) - diff(Cf(x), x) * exp(-x**2)) == 0
# Trou 1, leurres (a = 2x, b = 0, C = 1) : A = a' = 2, A = primitive de b = 0, A = x a(x) = 2x² ne donnent pas de solution.
for Afaux in [2 + 0 * x, 0 * x, 2 * x**2]:
    assert not sol1(exp(-Afaux), 2 * x, 0)
# Trou 2, leurres (a = 1, b = 1, A = x) : C = -e^{-x} vérifie C' = b e^{-A} mais y = -e^{-2x} n'est pas solution ;
# C = e^x + 1 donne la solution 1 + e^{-x} sans vérifier C = b e^{A} ; C = 1 vérifie C' + C = 1 mais e^{-x} n'est pas solution.
assert not sol1(-exp(-2 * x), 1, 1)
assert sol1((exp(x) + 1) * exp(-x), 1, 1) and not nul(exp(x) + 1 - exp(x))
assert not sol1(exp(-x), 1, 1)
# Exemple : x/2 + k/x solution de y' + y/x = 1 sur ]0,+oo[, et x/2 + 1/(2x) vérifie y(1) = 1.
assert sol1(x / 2 + k0 / x, 1 / x, 1, points=POS) and sol1(x / 2 + 1 / (2 * x), 1 / x, 1, 1, 1, POS)
# Contre-exemple (intervalle) : sur ]-oo,0[, x/2 + k/x est aussi solution.
assert sol1(x / 2 + k0 / x, 1 / x, 1, points=[-1, -2, Rational(-1, 3)])
# Contre-exemple / affirmation 3 : y' = 0, y(0) = 1, A = 1 : la formule donne e^{-1} ≠ 1.
assert exp(-1) * (1 + 0) != 1
# Affirmation 1 : C = e^x et C = e^x + 1 donnent les solutions 1 et 1 + e^{-x} de y' + y = 1.
assert sol1(1 + 0 * x, 1, 1) and sol1(1 + exp(-x), 1, 1)

# ed-l1-lineaires-variation-constante-ex1 (et v1, v2).
for a, b, y0, bon, leurres in [
    (1, exp(-x), 2, (x + 2) * exp(-x), [x * exp(-x) + 2, (2 - x) * exp(-x), (x + 2) * exp(x)]),
    (-2, exp(2 * x), -1, (x - 1) * exp(2 * x), [x * exp(2 * x) - 1, (1 - x) * exp(2 * x), (x - 1) * exp(-2 * x)]),
    (3, 2 * exp(-3 * x), 1, (2 * x + 1) * exp(-3 * x), [(x + 1) * exp(-3 * x), 2 * x * exp(-3 * x) + 1, (2 * x + 1) * exp(3 * x)]),
]:
    assert sol1(bon, a, b, 0, y0)
    for l in leurres:
        assert not sol1(l, a, b, 0, y0)

# ed-l1-lineaires-variation-constante-ex2 (et v1, v2) sur ]0,+oo[.
for a, b, x0, y0, bon, leurres in [
    (1 / x, 3 * x, 1, 2, x**2 + 1 / x, [x**2 + 1, 2 * x**2, x**2 + x]),
    (-1 / x, x, 1, 3, x**2 + 2 * x, [x**2 + 2, 3 * x**2, x**2 + 2 / x]),  # xy' - y = x² normalisée
    (2 / x, 1, 1, 1, x / 3 + 2 / (3 * x**2), [x / 3 + Rational(2, 3), x, x / 3 + 2 / (3 * x)]),
]:
    assert sol1(bon, a, b, x0, y0, POS)
    for l in leurres:
        assert not sol1(l, a, b, x0, y0, POS)
assert nul(x * diff(x**2 + 2 * x, x) - (x**2 + 2 * x) - x**2)  # v1 sous la forme de l'énoncé

# ed-l1-lineaires-variation-constante-ex3 (sur ]-π/2, π/2[, v1 aussi ; v2 sur R).
assert nul(diff(-log(cos(x)), x) - tan(x), TRIG)  # primitive de tan
for a, b, y0, bon, leurres, pts in [
    (tan(x), cos(x)**2, 1, (1 + sin(x)) * cos(x), [(1 + sin(x)) / cos(x), cos(x) + sin(x), cos(x)**2 + sin(x)], TRIG),
    (-tan(x), 1 / cos(x), 2, (x + 2) / cos(x), [(x + 2) * cos(x), (2 + sin(x)) / cos(x), x / cos(x) + 2], TRIG),
    (2 * x / (1 + x**2), 1 / (1 + x**2), -1, (x - 1) / (1 + x**2), [(x - 1) * (1 + x**2), (atan(x) - 1) / (1 + x**2), x - 1], POINTS),
]:
    assert sol1(bon, a, b, 0, y0, pts)
    for l in leurres:
        assert not sol1(l, a, b, 0, y0, pts)

# ed-l1-lineaires-variation-constante-ex4 : identité de la question 2 (ℓ = e^{-x}(ℓ + ∫_0^x ℓ e^t dt)) et majoration ε(1 - e^{X-x}).
ell, X, eps = symbols("ell X eps", real=True)
assert simplify(exp(-x) * (ell + integrate(ell * exp(t), (t, 0, x))) - ell) == 0
assert simplify(exp(-x) * integrate(eps * exp(t), (t, X, x)) - eps * (1 - exp(X - x))) == 0
# Exemple de contrôle : b = 1 + 1/(1+x)... on teste la formule de la question 1 sur b = 3 (y = 3 + (y0 - 3)e^{-x}).
y0s = symbols("y0s", real=True)
yq1 = exp(-x) * (y0s + integrate(3 * exp(t), (t, 0, x)))
assert sol1(yq1, 1, 3, 0, y0s) and limit(yq1, x, oo) == 3

# =============================================================================
# ed-l1-lineaires-superposition (carte)
# =============================================================================
# Trou 1, leurres : y' + y² (y1 = 1/x, b1 = 0 ; 2y1 non solution), y' + y + 1 (y1 = 0, b1 = 1 ; 2y1 donne 1 ≠ 2), y y' (y1 = x, 2y1 donne 4x ≠ 2x).
assert nul(diff(1 / x, x) + (1 / x)**2) and not nul(diff(2 / x, x) + (2 / x)**2)
assert nul(diff(0 * x, x) + 0 + 1 - 1) and not nul(diff(0 * x, x) + 0 + 1 - 2)
assert nul(x * diff(x, x) - x) and not nul(2 * x * diff(2 * x, x) - 2 * x)
# Trou 2, leurres (L(y) = y', y1 = x, b1 = 1, y2 = 0, b2 = 0) : λ = 2, μ = 0 contredit b1 + b2 ; λ = μ = 1 avec y2 = x contredit λμ b1 b2 ;
# λ = 1, μ = 0 contredit μ b1 + λ b2.
assert diff(2 * x, x) != 1 + 0
assert diff(x + x, x) != 1 * 1 * 1 * 1
assert diff(x, x) != 0 * 1 + 1 * 0
# Exemple 1 : x + e^x/2 solution de y'' + y = x + e^x.
assert sol2(x + exp(x) / 2, 1, 0, 1, x + exp(x))
# Exemple 2 : (cos + sin)/2 solution de y' + y = cos ; c = 1/(1+i) = (1-i)/2.
assert simplify(1 / (1 + I) - (1 - I) / 2) == 0
assert sol1((cos(x) + sin(x)) / 2, 1, cos(x))
# Contre-exemples : e^{-x} + e^{-2x} ; cos x et y' - i y.
yy = exp(-x) + exp(-2 * x)
assert nul(L1(yy, 1) + exp(-2 * x)) and nul(L1(yy, 2) - exp(-x))
assert nul(diff(exp(I * x), x) - I * exp(I * x)) and not nul(diff(cos(x), x) - I * cos(x))
# Affirmation 1 : 2x vérifie y'' + y = 2x ≠ x ; affirmation 3 : (y1 + y2)/2 (y1 = x, y2 = x + sin x) solution de y'' + y = x.
assert sol2(2 * x, 1, 0, 1, 2 * x) and not sol2(2 * x, 1, 0, 1, x)
assert sol2((x + x + sin(x)) / 2, 1, 0, 1, x)
# Affirmation 2 : z = e^{ix}/(χ(i)) avec χ(r) = r² + 3r + 2, Re z solution de y'' + 3y' + 2y = cos x.
z = exp(I * x) / (I**2 + 3 * I + 2)
assert nul(L2(z, 1, 3, 2) - exp(I * x))
assert sol2(re(z.expand(complex=True)), 1, 3, 2, cos(x))

# ed-l1-lineaires-superposition-ex1 (et v1, v2).
assert sol2(x, 1, 0, 1, x) and sol2(exp(x) / 2, 1, 0, 1, exp(x))
for l in [3 * x - 2 * exp(x), x - exp(x) / 2, 3 * x - exp(x) / 2]:
    assert not sol2(l, 1, 0, 1, 3 * x - 2 * exp(x))
assert sol2(3 * x - exp(x), 1, 0, 1, 3 * x - 2 * exp(x))
assert sol1(x - 1, 1, x) and sol1(exp(x) / 2, 1, exp(x))
assert sol1(2 * x - 2 + 2 * exp(x), 1, 2 * x + 4 * exp(x))
for l in [2 * x + 2 * exp(x), 2 * x - 2 + 4 * exp(x), x - 1 + exp(x) / 2]:
    assert not sol1(l, 1, 2 * x + 4 * exp(x))
assert sol2(-x, 1, 0, -1, x) and sol2(-cos(x) / 2, 1, 0, -1, cos(x))
assert sol2(-5 * x - cos(x), 1, 0, -1, 5 * x + 2 * cos(x))
for l in [-5 * x - 2 * cos(x), 5 * x + cos(x), -x - cos(x) / 2]:
    assert not sol2(l, 1, 0, -1, 5 * x + 2 * cos(x))

# ed-l1-lineaires-superposition-ex2 (et v1, v2) : c et la solution réelle.
assert sol1((cos(x) + sin(x)) / 2, 1, cos(x))
for l in [(cos(x) - sin(x)) / 2, cos(x) / 2, cos(x) + sin(x)]:
    assert not sol1(l, 1, cos(x))
assert simplify(1 / (I - 2) - (-2 - I) / 5) == 0
assert simplify(im(((-2 - I) / 5 * (cos(x) + I * sin(x))).expand()) + (2 * sin(x) + cos(x)) / 5) == 0
assert sol1(-(2 * sin(x) + cos(x)) / 5, -2, sin(x))
for l in [(2 * sin(x) + cos(x)) / 5, -(2 * cos(x) + sin(x)) / 5, -sin(x) / 2]:
    assert not sol1(l, -2, sin(x))
assert simplify(1 / (3 + 2 * I) - (3 - 2 * I) / 13) == 0
assert sol1((3 * cos(2 * x) + 2 * sin(2 * x)) / 13, 3, cos(2 * x))
for l in [(3 * cos(2 * x) - 2 * sin(2 * x)) / 13, (3 * cos(2 * x) + 2 * sin(2 * x)) / 5, cos(2 * x) / 3]:
    assert not sol1(l, 3, cos(2 * x))

# ed-l1-lineaires-superposition-ex3 (et v1, v2) : linéarisation, solution avec y(0) = 0.
assert nul(cos(x)**2 - (1 + cos(2 * x)) / 2) and nul(sin(x)**2 - (1 - cos(2 * x)) / 2)
assert sol1((cos(2 * x) + 2 * sin(2 * x)) / 5, 1, cos(2 * x))
assert sol1((-cos(2 * x) + 2 * sin(2 * x)) / 5, -1, cos(2 * x))
assert simplify(1 / (2 * I - 1) - (-1 - 2 * I) / 5) == 0
for a, b, bon, leurres in [
    (1, cos(x)**2, Rational(1, 2) + (cos(2 * x) + 2 * sin(2 * x)) / 10 - Rational(3, 5) * exp(-x),
     [Rational(1, 2) + (cos(2 * x) + 2 * sin(2 * x)) / 10 - exp(-x) / 2,
      Rational(1, 2) + (cos(2 * x) - 2 * sin(2 * x)) / 10 - Rational(3, 5) * exp(-x),
      Rational(1, 2) + (cos(2 * x) + 2 * sin(2 * x)) / 5 - Rational(7, 10) * exp(-x)]),
    (1, sin(x)**2, Rational(1, 2) - (cos(2 * x) + 2 * sin(2 * x)) / 10 - Rational(2, 5) * exp(-x),
     [Rational(1, 2) - (cos(2 * x) - 2 * sin(2 * x)) / 10 - Rational(2, 5) * exp(-x),
      Rational(1, 2) + (cos(2 * x) + 2 * sin(2 * x)) / 10 - Rational(3, 5) * exp(-x),
      Rational(1, 2) - (cos(2 * x) + 2 * sin(2 * x)) / 10 - exp(-x) / 2]),
    (-1, cos(x)**2, -Rational(1, 2) + (2 * sin(2 * x) - cos(2 * x)) / 10 + Rational(3, 5) * exp(x),
     [-Rational(1, 2) + (cos(2 * x) + 2 * sin(2 * x)) / 10 + Rational(2, 5) * exp(x),
      Rational(1, 2) + (2 * sin(2 * x) - cos(2 * x)) / 10 - Rational(2, 5) * exp(x),
      -Rational(1, 2) + (2 * sin(2 * x) - cos(2 * x)) / 10 + Rational(3, 5) * exp(-x)]),
]:
    assert sol1(bon, a, b, 0, 0)
    for l in leurres:
        assert not sol1(l, a, b, 0, 0)

# ed-l1-lineaires-superposition-ex4 (et v1, v2) : contrôle pour n = 1, 2, 3, et périodicité.
for n in [1, 2, 3]:
    ks = range(1, n + 1)
    s_cos = sum(cos(k * x) for k in ks)
    s_sin = sum(sin(k * x) for k in ks)
    for a, b, bon, leurres, hom in [
        (1, s_cos, sum((cos(k * x) + k * sin(k * x)) / (1 + k**2) for k in ks),
         [sum((cos(k * x) - k * sin(k * x)) / (1 + k**2) for k in ks), sum((cos(k * x) + k * sin(k * x)) / k**2 for k in ks)], exp(-x)),
        (2, s_cos, sum((2 * cos(k * x) + k * sin(k * x)) / (4 + k**2) for k in ks),
         [sum((2 * cos(k * x) - k * sin(k * x)) / (4 + k**2) for k in ks), sum((cos(k * x) + k * sin(k * x)) / (1 + k**2) for k in ks)], exp(-2 * x)),
        (-1, s_sin, -sum((sin(k * x) + k * cos(k * x)) / (1 + k**2) for k in ks),
         [sum((sin(k * x) + k * cos(k * x)) / (1 + k**2) for k in ks), -sum((cos(k * x) - k * sin(k * x)) / (1 + k**2) for k in ks)], exp(x)),
    ]:
        assert sol1(bon, a, b)
        assert simplify(bon.subs(x, x + 2 * pi) - bon) == 0  # 2π-périodique
        for l in leurres:
            assert not sol1(l, a, b)
        # Le leurre « homogène + y_p » est solution mais pas 2π-périodique.
        assert sol1(bon + hom, a, b) and simplify((bon + hom).subs(x, 2 * pi) - (bon + hom).subs(x, 0)) != 0

# =============================================================================
# ed-l1-lineaires-ordre2 (carte)
# =============================================================================
r = symbols("r")
# Démo : avec y = z e^{r1 x}, a y'' + b y' + c y = a (z'' + (r1 - r2) z') e^{r1 x} (r1, r2 racines, b = -a(r1+r2), c = a r1 r2).
a_, r1, r2 = symbols("a_ r1 r2")
zf = __import__("sympy").Function("z")(x)
b_, c_ = -a_ * (r1 + r2), a_ * r1 * r2
assert simplify(L2(zf * exp(r1 * x), a_, b_, c_) - a_ * (diff(zf, x, 2) + (r1 - r2) * diff(zf, x)) * exp(r1 * x)) == 0
# Démo : partie réelle (p+iq)e^{(α+iβ)x} + (s+it)e^{(α-iβ)x}.
p, q, s, tt, al, be = symbols("p q s tt al be", real=True)
expr = ((p + I * q) * exp((al + I * be) * x) + (s + I * tt) * exp((al - I * be) * x)).expand(complex=True)
assert simplify(re(expr) - exp(al * x) * ((p + s) * cos(be * x) + (tt - q) * sin(be * x))) == 0
# Démo : déterminants des trois systèmes de conditions initiales.
x0 = symbols("x0", real=True)
def wr(f, g):
    return simplify(Matrix([[f, g], [diff(f, x), diff(g, x)]]).det().subs(x, x0))
assert simplify(wr(exp(r1 * x), exp(r2 * x)) - (r2 - r1) * exp((r1 + r2) * x0)) == 0
assert simplify(wr(x * exp(r1 * x), exp(r1 * x)) + exp(2 * r1 * x0)) == 0
assert simplify(wr(exp(al * x) * cos(be * x), exp(al * x) * sin(be * x)) - be * exp(2 * al * x0)) == 0
# Trou 1 (racine double r = 1 : y'' - 2y' + y = 0) : x e^x solution (absente de λe^{rx}) ; e^{-x} et x² e^x non solutions.
assert sol2(x * exp(x), 1, -2, 1, 0) and not sol2(exp(-x), 1, -2, 1, 0) and not sol2(x**2 * exp(x), 1, -2, 1, 0)
# Trou 2 (y'' + 2y' + 5y = 0, α = -1, β = 2) : bonne forme ; leurres faux.
assert sol2(exp(-x) * cos(2 * x), 1, 2, 5, 0) and sol2(exp(-x) * sin(2 * x), 1, 2, 5, 0)
for l in [exp(2 * x) * cos(-x), exp(-x) + 0 * x, exp(-x) * cos(x)]:
    assert not sol2(l, 1, 2, 5, 0)
# Trou 3, leurres / piège / contre-exemple (y'' + y = 0) : sin x et 2 sin x ont même y(0) et même y(π) ; cos x et 2 cos x même y'(0).
assert sol2(sin(x), 1, 0, 1, 0, [(0, 0, 0), (0, pi, 0)]) and sol2(2 * sin(x), 1, 0, 1, 0, [(0, 0, 0), (0, pi, 0)])
assert sol2(cos(x), 1, 0, 1, 0, [(1, 0, 0)]) and sol2(2 * cos(x), 1, 0, 1, 0, [(1, 0, 0)])
# Exemples.
assert sol2(exp(2 * x) - exp(x), 1, -3, 2, 0, [(0, 0, 0), (1, 0, 1)])
assert sol2(exp(-x) * sin(2 * x), 1, 2, 5, 0)
w = symbols("w", positive=True)
assert nul(diff(sin(w * x), x, 2) + w**2 * sin(w * x))
# Contre-exemple a = 0 : C e^x vérifie y'(0) = y(0).
assert simplify(diff(C * exp(x), x).subs(x, 0) - C) == 0
# Affirmation 2 : cos solution de y'' + y = 0 sans limite (cos(2πn) = 1, cos(2πn + π) = -1).
assert sol2(cos(x), 1, 0, 1, 0)
# Affirmation 3 (b, c > 0) : contrôle sur des exemples (racines de partie réelle < 0).
from sympy import roots as racines
for bb, cc in [(1, 1), (3, 2), (2, 1), (Rational(1, 10), 5)]:
    assert all(re(rr) < 0 for rr in racines(r**2 + bb * r + cc, r))
# Affirmation 4 : λ cos x + μ sin x = A cos(x - φ).
phi, A_ = symbols("phi A_", real=True)
assert simplify((A_ * cos(x - phi)).expand(trig=True) - (A_ * cos(phi) * cos(x) + A_ * sin(phi) * sin(x))) == 0

# ed-l1-lineaires-ordre2-ex1 (et v1, v2) : bases de solutions ; leurres.
assert sol2(exp(2 * x), 1, -5, 6, 0) and sol2(exp(3 * x), 1, -5, 6, 0)
for l in [exp(-2 * x), exp(-3 * x), exp(Rational(5, 2) * x) * cos(x / 2)]:
    assert not sol2(l, 1, -5, 6, 0)
assert sol2(exp(-x), 1, 2, 1, 0) and sol2(x * exp(-x), 1, 2, 1, 0)
assert not sol2(x * exp(x), 1, 2, 1, 0) and not sol2(exp(x), 1, 2, 1, 0)  # leurres (λx+μ)e^x et λe^{-x}+μe^x ; λe^{-x} manque xe^{-x}
assert racines(r**2 + 4 * r + 13, r) == {-2 + 3 * I: 1, -2 - 3 * I: 1}
assert sol2(exp(-2 * x) * cos(3 * x), 1, 4, 13, 0) and sol2(exp(-2 * x) * sin(3 * x), 1, 4, 13, 0)
for l in [exp(-3 * x) * cos(2 * x), exp(2 * x) * cos(3 * x), exp(3 * x)]:
    assert not sol2(l, 1, 4, 13, 0)

# ed-l1-lineaires-ordre2-ex2 (et v1, v2) : problèmes de Cauchy.
for (a, b, c), cond, bon, leurres in [
    ((1, 0, 4), [(0, 0, 1), (1, 0, 2)], cos(2 * x) + sin(2 * x), [cos(2 * x) + 2 * sin(2 * x), cos(2 * x) - sin(2 * x), exp(2 * x)]),
    ((1, -3, 2), [(0, 0, 0), (1, 0, 1)], exp(2 * x) - exp(x), [exp(x) - exp(2 * x), x * exp(x), exp(-x) - exp(-2 * x)]),
    ((1, 2, 2), [(0, 0, 1), (1, 0, 0)], exp(-x) * (cos(x) + sin(x)), [exp(-x) * cos(x), exp(-x) * (cos(x) - sin(x)), cos(x) + sin(x)]),
]:
    assert sol2(bon, a, b, c, 0, cond)
    for l in leurres:
        assert not sol2(l, a, b, c, 0, cond)

# ed-l1-lineaires-ordre2-ex3 (et v1, v2) : conditions aux deux bouts.
for (a, b, c), cond, bon, leurres in [
    ((1, 2, 5), [(0, 0, 0), (0, pi / 4, 1)], exp(pi / 4 - x) * sin(2 * x),
     [exp(-x) * sin(2 * x), exp(x - pi / 4) * sin(2 * x), sqrt(2) * exp(pi / 4 - x) * sin(x)]),
    ((1, 0, -1), [(0, 0, 1), (0, 1, 0)], (exp(1 - x) - exp(x - 1)) / (E - exp(-1)),
     [1 - x, (exp(1 - x) + exp(x - 1)) / (E + exp(-1)), cos(x) - cos(1) / sin(1) * sin(x)]),
    ((1, -2, 2), [(0, 0, 0), (0, pi / 2, 1)], exp(x - pi / 2) * sin(x), [exp(x) * sin(x), exp(pi / 2 - x) * sin(x), sin(x)]),
]:
    assert sol2(bon, a, b, c, 0, cond)
    for l in leurres:
        assert not sol2(l, a, b, c, 0, cond)
# v1 : λ = 1/(1-e²) = -e^{-1}/(e - e^{-1}), μ = e²/(e²-1) = e/(e - e^{-1}).
assert simplify(1 / (1 - E**2) + exp(-1) / (E - exp(-1))) == 0 and simplify(E**2 / (E**2 - 1) - E / (E - exp(-1))) == 0
# Remarque : μ sin x vérifie y(0) = y(π) = 0 pour y'' + y = 0.
assert sol2(mu * sin(x), 1, 0, 1, 0, [(0, 0, 0), (0, pi, 0)])

# ed-l1-lineaires-ordre2-ex4 (et v1, v2) : valeurs propres.
om = symbols("om", positive=True)
# k < 0 : e^{ωπ} - e^{-ωπ} > 0 et e^{ω} - e^{-ω} > 0 et e^{ωπ} + e^{-ωπ} > 0 (contrôle sur quelques ω).
for o in [Rational(1, 5), 1, 3]:
    assert exp(o * pi) - exp(-o * pi) > 0 and exp(o) - exp(-o) > 0
# Bonnes réponses : pour n = 1, 2, 3, une solution non nulle existe.
for n in [1, 2, 3]:
    assert sol2(sin(n * x), 1, 0, n**2, 0, [(0, 0, 0), (0, pi, 0)])
    assert sol2(sin(n * pi * x), 1, 0, n**2 * pi**2, 0, [(0, 0, 0), (0, 1, 0)])
for n in [0, 1, 2]:
    w_ = n + Rational(1, 2)
    assert sol2(sin(w_ * x), 1, 0, w_**2, 0, [(0, 0, 0), (1, pi, 0)])
# Leurres : une valeur de k sans solution non nulle (après y(0) = 0, il reste μ sin(ωx), condition au bout non nulle).
assert sin(sqrt(2) * pi) != 0          # « k > 0 quelconque » : k = 2
assert sin(sqrt(2) * pi) != 0          # « k = n » : n = 2
assert exp(pi) - exp(-pi) != 0         # « k = -n² » : n = 1 (k < 0)
assert sin(1) != 0                     # v1 « k = n² » : n = 1, sin(1·1) ≠ 0
assert sin(sqrt(pi)) != 0              # v1 « k = nπ » : n = 1, ω = √π
assert sin(pi * 1) == 0 and sin(sqrt(pi) * 1) != 0  # v1 « k = n²π » : n = 1, ω = √π
assert cos(1 * pi) != 0                # v2 « k = n² » : n = 1, cos(π) ≠ 0
assert N(cos(sqrt(Rational(1, 2)) * pi)) != 0  # v2 « k = n + 1/2 » : n = 0, ω = 1/√2
assert cos(1 * pi) != 0                # v2 « k = (2n+1)² » : n = 0, ω = 1

# =============================================================================
# ed-l1-lineaires-exp-polynome (carte)
# =============================================================================
m_ = symbols("m_")
# Démo : (z e^{mx}) donne (a z'' + χ'(m) z' + χ(m) z) e^{mx}.
aa, bb2, cc2 = symbols("aa bb2 cc2")
chi = aa * m_**2 + bb2 * m_ + cc2
assert simplify(L2(zf * exp(m_ * x), aa, bb2, cc2) - (aa * diff(zf, x, 2) + (2 * aa * m_ + bb2) * diff(zf, x) + chi * zf) * exp(m_ * x)) == 0
# Exemples : -x e^x pour y'' - 3y' + 2y = e^x ; x sin x pour y'' + y = 2cos x ; -x² - 2 pour y'' - y = x².
assert sol2(-x * exp(x), 1, -3, 2, exp(x))
assert simplify(re((-I * x * exp(I * x)).expand(complex=True)) - x * sin(x)) == 0
assert sol2(x * sin(x), 1, 0, 1, 2 * cos(x))
assert sol2(-x**2 - 2, 1, 0, -1, x**2)
# Contre-exemples : α e^x donne 0 ; α cos x + β sin x donne 0 ; ln x + λx + μ solution de y'' = -1/x² ; ln x / x → 0.
al2, be2 = symbols("al2 be2")
assert nul(L2(al2 * exp(x), 1, -3, 2))
assert nul(L2(al2 * cos(x) + be2 * sin(x), 1, 0, 1))
assert nul(diff(log(x) + lam * x + mu, x, 2) + 1 / x**2, POS)
assert limit(log(x) / x, x, oo) == 0
# Trou 1 et 2, leurres (unicité de la solution polynôme × e^{mx} à un terme homogène près) :
#  y'' + y = x : la solution polynomiale est x (unique, m = 0 non racine) ⇒ « k = deg P », « Q de degré n + 1 », « n + 2 », « au plus n - 1 » faux ;
Pa, Pb, Pc = symbols("Pa Pb Pc")
from sympy import solve
sols = solve([c for c in __import__("sympy").Poly(L2(Pa * x**2 + Pb * x + Pc, 1, 0, 1) - x, x).all_coeffs()], [Pa, Pb, Pc], dict=True)
assert sols == [{Pa: 0, Pb: 1, Pc: 0}]
#  y'' - y = e^{2x} (deux racines réelles, m = 2 non racine) : seule e^{2x}/3 est de la forme c·e^{2x} ; c x² e^{2x} n'est pas solution ;
assert sol2(exp(2 * x) / 3, 1, 0, -1, exp(2 * x))
assert solve(__import__("sympy").Poly(simplify(L2(Pa * x**2 * exp(2 * x), 1, 0, -1) * exp(-2 * x) - 1), x).all_coeffs(), [Pa]) == []
#  y'' - y = e^x (m = 1 racine simple, n = 0) : x(αx + β)e^x solution ⇔ α = 0 (« Q de degré n + k = 1 » faux).
eqs = __import__("sympy").Poly(simplify(L2(x * (Pa * x + Pb) * exp(x), 1, 0, -1) * exp(-x) - 1), x).all_coeffs()
assert solve(eqs, [Pa, Pb], dict=True) == [{Pa: 0, Pb: Rational(1, 2)}]
#  y'' = 1 (m = 0 racine double, n = 0) : x²(αx + β) solution ⇔ α = 0, β = 1/2 (« Q de degré n + 2 » ou « n + 1 » faux).
eqs2 = __import__("sympy").Poly(diff(x**2 * (Pa * x + Pb), x, 2) - 1, x).all_coeffs()
assert solve(eqs2, [Pa, Pb], dict=True) == [{Pa: 0, Pb: Rational(1, 2)}]
#  « k = 0 toujours » : y'' - 3y' + 2y = e^x sans solution c e^x (ci-dessus).
# Affirmations : -x e^x (ci-dessus) ; α cos 2x + β sin 2x donne 0 pour y'' + 4y ; x sin(2x)/4 solution de y'' + 4y = cos 2x ;
# -x² - 2 (ci-dessus) ; x² e^{-x}/2 solution de y' + y = x e^{-x}.
assert nul(L2(al2 * cos(2 * x) + be2 * sin(2 * x), 1, 0, 4))
assert sol2(x * sin(2 * x) / 4, 1, 0, 4, cos(2 * x))
assert sol1(x**2 * exp(-x) / 2, 1, x * exp(-x))
assert nul(diff(-x * exp(x), x) + (x + 1) * exp(x)) and nul(diff(-x * exp(x), x, 2) + (x + 2) * exp(x))

# ed-l1-lineaires-exp-polynome-ex1 (et v1, v2).
for (a, b, c), f, bon, leurres in [
    ((1, -3, 2), 2 * x + 1, x + 2, [x + Rational(1, 2), x - 1, 2 * x + 1]),
    ((1, 1, -2), 4 * x, -2 * x - 1, [-2 * x, -2 * x + 1, 2 * x + 1]),
    ((1, 0, 4), 8 * x**2, 2 * x**2 - 1, [2 * x**2, 2 * x**2 + 1, 2 * x**2 - 4]),
]:
    assert sol2(bon, a, b, c, f)
    for l in leurres:
        assert not sol2(l, a, b, c, f)

# ed-l1-lineaires-exp-polynome-ex2 (et v1, v2) : racine simple, χ'(m).
for (a, b, c), m, f, bon, leurres in [
    ((1, -1, -2), 2, 3 * exp(2 * x), x * exp(2 * x), [-x * exp(2 * x), 3 * x * exp(2 * x), Rational(3, 2) * x * exp(2 * x)]),
    ((1, 1, -6), 2, 10 * exp(2 * x), 2 * x * exp(2 * x), [x * exp(2 * x), Rational(10, 3) * x * exp(2 * x), 10 * x * exp(2 * x)]),
    ((1, 0, -4), -2, 8 * exp(-2 * x), -2 * x * exp(-2 * x), [2 * x * exp(-2 * x), -x * exp(-2 * x), -8 * x * exp(-2 * x)]),
]:
    assert a * m**2 + b * m + c == 0 and 2 * a * m + b != 0  # racine simple
    assert sol2(bon, a, b, c, f)
    for l in leurres:
        assert not sol2(l, a, b, c, f)
assert 2 * 2 - 1 == 3 and 2 * 2 + 1 == 5 and 2 * (-2) == -4  # χ'(m) cités

# ed-l1-lineaires-exp-polynome-ex3 (et v1, v2) : résonance, y(0) = y'(0) = 0.
CI0 = [(0, 0, 0), (1, 0, 0)]
assert nul(L2(-I * x * exp(I * x), 1, 0, 1) - 2 * exp(I * x))
assert nul(L2(-I * x * exp(2 * I * x), 1, 0, 4) - 4 * exp(2 * I * x))
assert simplify(im((-I * x * exp(2 * I * x)).expand(complex=True)) + x * cos(2 * x)) == 0
for (a, b, c), f, bon, leurres in [
    ((1, 0, 1), 2 * cos(x), x * sin(x), [x * cos(x), x * sin(x) / 2, sin(x) - x * cos(x)]),
    ((1, 0, 4), 4 * sin(2 * x), sin(2 * x) / 2 - x * cos(2 * x), [-x * cos(2 * x), x * sin(2 * x), sin(2 * x) - x * cos(2 * x)]),
    ((1, -2, 1), 2 * exp(x), x**2 * exp(x), [x * exp(x), 2 * x**2 * exp(x), (x**2 - x) * exp(x)]),
]:
    assert sol2(bon, a, b, c, f, CI0)
    for l in leurres:
        assert not sol2(l, a, b, c, f, CI0)
assert sol2(-x * cos(2 * x), 1, 0, 4, 4 * sin(2 * x))  # y_p de v1

# ed-l1-lineaires-exp-polynome-ex4 (et v1, v2) : polynôme × exponentielle en résonance.
assert sol2(-(x**2 + 3 * x) * exp(x), 1, -3, 2, (2 * x + 1) * exp(x))
assert sol2(x**2 * exp(x), 1, 1, -2, (6 * x + 2) * exp(x))
assert sol2(x**3 * exp(2 * x), 1, -4, 4, 6 * x * exp(2 * x))
for (a, b, c), f, cond, bon, leurres in [
    ((1, -3, 2), (2 * x + 1) * exp(x), [(0, 0, 0), (1, 0, 0)], 3 * exp(2 * x) - (x**2 + 3 * x + 3) * exp(x),
     [-(x**2 + 3 * x) * exp(x), exp(2 * x) - (x**2 + x + 1) * exp(x), (x**2 + 3 * x + 3) * exp(x) - 3 * exp(2 * x)]),
    ((1, 1, -2), (6 * x + 2) * exp(x), [(0, 0, 0), (1, 0, 3)], (x**2 + 1) * exp(x) - exp(-2 * x),
     [x**2 * exp(x), (x**2 - 1) * exp(x) + exp(-2 * x), (x**2 + 3 * x) * exp(x)]),
    ((1, -4, 4), 6 * x * exp(2 * x), [(0, 0, 1), (1, 0, 0)], (x**3 - 2 * x + 1) * exp(2 * x),
     [(x**3 + 1) * exp(2 * x), (x**2 - 2 * x + 1) * exp(2 * x), (3 * x**2 - 2 * x + 1) * exp(2 * x)]),
]:
    assert sol2(bon, a, b, c, f, cond)
    for l in leurres:
        assert not sol2(l, a, b, c, f, cond)

print("ed_l1_lineaires : toutes les vérifications passent.")
