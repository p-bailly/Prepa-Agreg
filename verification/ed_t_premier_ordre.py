"""Vérifications calculatoires du chapitre « ed-t-premier-ordre » (équations différentielles, Terminale).

Cartes et exercices de content/{cartes,exercices}/equations-differentielles/terminale-ed-t-premier-ordre.json.
Chaque assert cite l'identifiant de la carte ou de l'exercice concerné (v1, v2 = variantes).
Lancer : python3 verification/ed_t_premier_ordre.py (doit finir sans erreur).
"""

from sympy import (E, Function, Rational, cos, diff, exp, limit, log, oo, simplify, sin, sqrt, symbols)

x, t, C, k = symbols("x t C k", real=True)
n = symbols("n", positive=True, integer=True)


def nul(e):
    return simplify(e) == 0


def solution(y, second_membre, x0, y0, var=x):
    """y est solution de y' = second_membre(var, y) avec y(x0) = y0."""
    return nul(diff(y, var) - second_membre(var, y)) and nul(y.subs(var, x0) - y0)


def primitive(F, f, x0, y0):
    return nul(diff(F, x) - f) and nul(F.subs(x, x0) - y0)


# ---------------------------------------------------------------------------
# ed-t-premier-ordre-ay (carte)
# ---------------------------------------------------------------------------
a = symbols("a", real=True)
# Démo : (C e^{ax})' = a C e^{ax} ; g = y e^{-ax} a une dérivée nulle quand y = C e^{ax}.
assert nul(diff(C * exp(a * x), x) - a * C * exp(a * x))
y0, x0 = symbols("y0 x0", real=True)
assert solution(y0 * exp(a * (x - x0)), lambda v, y: a * y, x0, y0)
# Trou 2 : leurres e^{ax}+C et C x e^{ax} ne sont pas solutions (ex. a = 1, C = 1).
assert not nul(diff(exp(x) + 1, x) - (exp(x) + 1))
assert not nul(diff(x * exp(x), x) - x * exp(x))
# Exemple : y = 5 e^{-2(x-1)} solution de y' = -2y, y(1) = 5 ; demi-vie T = ln2/λ.
assert solution(5 * exp(-2 * (x - 1)), lambda v, y: -2 * y, 1, 5)
lam = symbols("lam", positive=True)
assert nul(exp(-lam * log(2) / lam) - Rational(1, 2))
# Contre-exemple (intervalle) : e^x et 2e^x vérifient y' = y chacun sur leur demi-droite.
assert nul(diff(2 * exp(x), x) - 2 * exp(x))
# Contre-exemple (a constant) / affirmation 3 : e^{2x^2} n'est pas solution de y' = 2xy ; e^{x^2} l'est.
assert nul(diff(exp(2 * x**2), x) - 4 * x * exp(2 * x**2))
assert not nul((diff(exp(2 * x**2), x) - 2 * x * exp(2 * x**2)).subs(x, 1))
assert nul(diff(exp(x**2), x) - 2 * x * exp(x**2))
# Contre-exemple / affirmation 4 : -1/x solution de y' = y^2, mais pas -2/x.
assert nul(diff(-1 / x, x) - (-1 / x) ** 2)
assert not nul(diff(-2 / x, x) - (-2 / x) ** 2)
# Affirmation 2 : -e^{-x} solution de y' = -y et croissante.
assert solution(-exp(-x), lambda v, y: -y, 0, -1) and nul(diff(-exp(-x), x) - exp(-x))

# ed-t-premier-ordre-ay-ex1 (et v1, v2) : la bonne proposition, et les leurres faux.
cas = [
    (lambda v, y: -3 * y, 4, 4 * exp(-3 * x), [4 * exp(3 * x), exp(-3 * x) + 3, -12 * exp(-3 * x)]),
    (lambda v, y: 5 * y, -2, -2 * exp(5 * x), [-2 * exp(-5 * x), exp(5 * x) - 3, 5 * exp(-2 * x)]),
    (lambda v, y: y / 2, 3, 3 * exp(x / 2), [3 * exp(2 * x), exp(x / 2) + 2, exp(3 * x) / 2]),
]
for sm, val, bon, leurres in cas:
    assert solution(bon, sm, 0, val)
    for l in leurres:
        assert not solution(l, sm, 0, val)

# ed-t-premier-ordre-ay-ex2 (et v1, v2) : f(0).
cas = [
    (2, 1, 3, 3 * exp(-2), [3 * exp(2), 3, 6 * exp(-2)]),
    (-1, 2, 5, 5 * exp(2), [5 * exp(-2), 5, -5 * exp(2)]),
    (3, -1, 2, 2 * exp(3), [2 * exp(-3), 2, 6 * exp(3)]),
]
for aa, xx0, yy0, bon, leurres in cas:
    f = yy0 * exp(aa * (x - xx0))
    assert solution(f, lambda v, y: aa * y, xx0, yy0)
    assert nul(f.subs(x, 0) - bon)
    for l in leurres:
        assert not nul(f.subs(x, 0) - l)

# ed-t-premier-ordre-ay-ex3 (et v1, v2) : a à partir de deux mesures.
# Principal : f(0)=8, f(2)=2 -> a = -ln2, f(3) = 1 ; leurres -1 (droite), 1/8 (e^a = 1/4), 1/2.
f = 8 * exp(-log(2) * t)
assert nul(f.subs(t, 2) - 2) and nul(f.subs(t, 3) - 1)
assert 8 + 3 * (2 - 8) / 2 == -1 and nul(8 * Rational(1, 4) ** 3 - Rational(1, 8))
for l in [-1, Rational(1, 8), Rational(1, 2)]:
    assert not nul(f.subs(t, 3) - l)
# v1 : f(0)=3, f(2)=12 -> a = ln2, f(3) = 24 ; leurres 33/2 (droite), 192 (e^a = 4), 36.
f = 3 * exp(log(2) * t)
assert nul(f.subs(t, 2) - 12) and nul(f.subs(t, 3) - 24)
assert 3 + 3 * Rational(12 - 3, 2) == Rational(33, 2) and 3 * 4**3 == 192
for l in [Rational(33, 2), 192, 36]:
    assert not nul(f.subs(t, 3) - l)
# v2 : f(0)=5, f(3)=40 -> a = ln2, f(4) = 80 ; leurres 155/3 (droite), 20480 (e^a = 8), 45.
f = 5 * exp(log(2) * t)
assert nul(f.subs(t, 3) - 40) and nul(f.subs(t, 4) - 80)
assert 5 + 4 * Rational(40 - 5, 3) == Rational(155, 3) and 5 * 8**4 == 20480
for l in [Rational(155, 3), 20480, 45]:
    assert not nul(f.subs(t, 4) - l)

# ed-t-premier-ordre-ay-ex4 (et v1, v2) : inégalités différentielles.
# Principal : f' <= 2f ; g = f e^{-2x}, g' = (f' - 2f) e^{-2x}. Exemple f = e^x réfute les trois leurres en x = -1.
F = Function("F")
assert nul(diff(F(x) * exp(-2 * x), x) - (diff(F(x), x) - 2 * F(x)) * exp(-2 * x))
f = exp(x)
assert all((diff(f, x) - 2 * f).subs(x, v) <= 0 for v in [-3, -1, 0, 1, 3])
xm = -1
assert f.subs(x, xm) >= f.subs(x, 0) * exp(2 * xm)                       # bonne : f(x) >= f(0) e^{2x}
assert not f.subs(x, xm) <= f.subs(x, 0) * exp(2 * xm)                   # leurre <=
assert not nul(f.subs(x, xm) - f.subs(x, 0) * exp(2 * xm))               # leurre =
assert not f.subs(x, xm) >= f.subs(x, 0) * exp(-2 * xm)                  # leurre f(0) e^{-2x}
# v1 : f' >= -f, g = f e^{x} ; f = 1 réfute les leurres en x = -1.
assert nul(diff(F(x) * exp(x), x) - (diff(F(x), x) + F(x)) * exp(x))
f = 1 + 0 * x
assert diff(f, x) >= -f
assert f.subs(x, xm) <= f.subs(x, 0) * exp(-xm)                          # bonne
assert not f.subs(x, xm) >= f.subs(x, 0) * exp(-xm)
assert not nul(f.subs(x, xm) - f.subs(x, 0) * exp(-xm))
assert not f.subs(x, xm) <= f.subs(x, 0) * exp(xm)
# v2 : f' <= -3f, g = f e^{3x} ; f = -e^{-2x} réfute les leurres en x = -1.
assert nul(diff(F(x) * exp(3 * x), x) - (diff(F(x), x) + 3 * F(x)) * exp(3 * x))
f = -exp(-2 * x)
assert nul(diff(f, x) + 3 * f + exp(-2 * x))  # f' + 3f = -e^{-2x} <= 0
assert f.subs(x, xm) >= f.subs(x, 0) * exp(-3 * xm)                      # bonne
assert not f.subs(x, xm) <= f.subs(x, 0) * exp(-3 * xm)
assert not nul(f.subs(x, xm) - f.subs(x, 0) * exp(-3 * xm))
assert not f.subs(x, xm) >= f.subs(x, 0) * exp(3 * xm)

# ed-t-premier-ordre-ay-ex5 : f(x+y) = f(x) f(y). Les e^{ax} et 0 conviennent ; 2e^{ax} non (f(0) = 2, f(0)^2 = 4).
yy = symbols("yy", real=True)
assert nul(exp(a * (x + yy)) - exp(a * x) * exp(a * yy))
assert 2 * exp(0) != (2 * exp(0)) ** 2
# Leurre « e^x seulement » : e^{2x} convient aussi.
assert nul(exp(2 * (x + yy)) - exp(2 * x) * exp(2 * yy))

# ---------------------------------------------------------------------------
# ed-t-premier-ordre-ay-b (carte)
# ---------------------------------------------------------------------------
b = symbols("b", real=True)
an = symbols("an", nonzero=True, real=True)
# Énoncé / démo : C e^{ax} - b/a solution ; solution avec y(x0) = y0.
assert nul(diff(C * exp(an * x) - b / an, x) - (an * (C * exp(an * x) - b / an) + b))
assert solution((y0 + b / an) * exp(an * (x - x0)) - b / an, lambda v, y: an * y + b, x0, y0)
# Trou 2 : leurres b/a, -a/b, b non solutions constantes en général (a = 1, b = 1 : il faut 0 = a k + b).
for (aa, bb) in [(2, 3)]:
    assert aa * Rational(-bb, aa) + bb == 0                                # -b/a
    for kk in [Rational(bb, aa), Rational(-aa, bb), bb]:
        assert aa * kk + bb != 0
# Exemple : y' = -2y + 6, y(0) = 1 -> 3 - 2e^{-2x}, limite 3.
assert solution(3 - 2 * exp(-2 * x), lambda v, y: -2 * y + 6, 0, 1)
assert limit(3 - 2 * exp(-2 * x), x, oo) == 3
# Contre-exemples : y' = y + x, solution -x - 1 ; y' = 1, solution x ; e^x - 1 et 2e^x - 1 solutions de y' = y + 1.
assert nul(diff(-x - 1, x) - ((-x - 1) + x))
assert nul(diff(2 * exp(x) - 1, x) - ((2 * exp(x) - 1) + 1)) and nul(diff(exp(x) - 1, x) - (exp(x) - 1 + 1))
# Affirmation 2 : constante 1 solution de y' = y - 1 ; affirmation 3 : -1 solution de y' = y + 1, pas -2.
assert 1 - 1 == 0 and -1 + 1 == 0 and -2 + 1 != 0

# ed-t-premier-ordre-ay-b-ex1 (et v1, v2).
cas = [
    (lambda v, y: 3 * y - 6, 5, 3 * exp(3 * x) + 2, [7 * exp(3 * x) - 2, 5 * exp(3 * x) + 2, 3 * exp(-3 * x) + 2]),
    (lambda v, y: -2 * y + 8, 1, -3 * exp(-2 * x) + 4, [5 * exp(-2 * x) - 4, exp(-2 * x) + 4, -3 * exp(2 * x) + 4]),
    (lambda v, y: y / 2 + 1, 0, 2 * exp(x / 2) - 2, [-2 * exp(x / 2) + 2, 2 * exp(2 * x) - 2, exp(x / 2) - 2]),
]
for sm, val, bon, leurres in cas:
    assert solution(bon, sm, 0, val)
    for l in leurres:
        assert not solution(l, sm, 0, val)

# ed-t-premier-ordre-ay-b-ex2 (et v1, v2) : refroidissement, instant où theta = seuil.
cas = [  # (theta0, piece, k, seuil, bonne, leurres)
    (80, 20, Rational(1, 10), 50, 10 * log(2), [log(2) / 10, 10 * log(Rational(6, 5)), 10 * log(Rational(8, 5))]),
    (95, 25, Rational(1, 5), 60, 5 * log(2), [log(2) / 5, 5 * log(Rational(19, 12)), 5 * log(Rational(7, 6))]),
    (90, 18, Rational(1, 20), 42, 20 * log(3), [log(3) / 20, 20 * log(Rational(15, 7)), 20 * log(Rational(12, 7))]),
]
for th0, p, kk, s, bon, leurres in cas:
    th = p + (th0 - p) * exp(-kk * t)
    assert solution(th, lambda v, y: -kk * (y - p), 0, th0, var=t)
    assert nul(th.subs(t, bon) - s)
    for l in leurres:
        assert not nul(th.subs(t, l) - s)
# Leurres « oubli de la pièce » : (theta0 - piece) e^{-kt} = seuil, puis theta0 e^{-kt} = seuil.
assert nul(60 * exp(-10 * log(Rational(6, 5)) / 10) - 50) and nul(80 * exp(-log(Rational(8, 5))) - 50)
assert nul(70 * exp(-log(Rational(7, 6))) - 60) and nul(95 * exp(-log(Rational(19, 12))) - 60)
assert nul(72 * exp(-log(Rational(12, 7))) - 42) and nul(90 * exp(-log(Rational(15, 7))) - 42)
assert Rational(9, 10) == Rational(1, 20) * 18  # v2 : theta' = -0,05 theta + 0,9

# ed-t-premier-ordre-ay-b-ex3 (et v1, v2) : second membre affine.
cas = [  # (a, alpha, beta [second membre alpha x + beta], p, q, y(0), bonne, leurres)
    (-1, 2, 3, 2, 1, 4, 2 * x + 1 + 3 * exp(-x), [2 * x + 1 + 4 * exp(-x), 2 * x + 1 + 3 * exp(x), 2 * x + 3 + exp(-x)]),
    (2, -4, 0, 2, 1, 3, 2 * x + 1 + 2 * exp(2 * x), [2 * x + 1 + 3 * exp(2 * x), 2 * x + 1 + 2 * exp(-2 * x), 2 * x + 2 + exp(2 * x)]),
    (-3, 3, 4, 1, 1, 0, x + 1 - exp(-3 * x), [x + 1 - exp(3 * x), x + 1 + exp(-3 * x), x - 1 + exp(-3 * x)]),
]
for aa, al, be, p, q, val, bon, leurres in cas:
    sm = lambda v, y, aa=aa, al=al, be=be: aa * y + al * v + be
    assert nul(diff(p * x + q, x) - sm(x, p * x + q))
    assert solution(bon, sm, 0, val)
    for l in leurres:
        assert not solution(l, sm, 0, val)

# ed-t-premier-ordre-ay-b-ex4 (et v1, v2) : chute, v = (10/k)(1 - e^{-kt}).
kp = symbols("kp", positive=True)
vv = 10 / kp * (1 - exp(-kp * t))
assert solution(vv, lambda s, y: 10 - kp * y, 0, 0, var=t)
assert limit(vv, t, oo) == 10 / kp
cas = [  # (vitesse limite, cible, bonne, leurres)
    (50, 25, 5 * log(2), [log(2) / 5, Rational(5, 2), 50 * log(2)]),
    (40, 30, 4 * log(4), [log(4) / 4, 3, 40 * log(4)]),
    (20, 15, 2 * log(4), [log(4) / 2, Rational(3, 2), 20 * log(4)]),
]
for vl, cible, bon, leurres in cas:
    kk = Rational(10, vl)
    v = vv.subs(kp, kk)
    assert nul(v.subs(t, bon) - cible)
    for l in leurres:
        assert not nul(v.subs(t, l) - cible)
assert nul(4 * log(4) - 8 * log(2)) and nul(2 * log(4) - 4 * log(2))
# Remarque : si k < 0, (10/k)(1 - e^{-kt}) -> +oo.
assert limit(10 / (-1) * (1 - exp(t)), t, oo) == oo

# ed-t-premier-ordre-ay-b-ex5 (et v1, v2) : méthode d'Euler.
cas = [  # (a, b, y0, constante, u_n formule, limite, leurres)
    (-1, 1, 0, 1, 1 - (1 - 1 / n) ** n, 1 - exp(-1), [0, 1, exp(-1)]),
    (2, 2, 0, -1, (1 + 2 / n) ** n - 1, exp(2) - 1, [0, oo, exp(2)]),
    (-3, 3, 2, 1, 1 + (1 - 3 / n) ** n, 1 + exp(-3), [2, 1, 1 + exp(3)]),
]
for aa, bb, yi, cst, un, lim, leurres in cas:
    yex = (yi - cst) * exp(aa * x) + cst
    assert solution(yex, lambda v, y: aa * y + bb, 0, yi)
    assert nul(yex.subs(x, 1) - lim)
    # récurrence d'Euler calculée exactement pour quelques n, comparée à la formule fermée
    for N in [4, 5, 10]:
        h = Rational(1, N)
        u = Rational(yi)
        for _ in range(N):
            u = u + h * (aa * u + bb)
        assert nul(u - un.subs(n, N))
    assert simplify(limit(un, n, oo) - lim) == 0
    for l in leurres:
        assert l != lim and not (l != oo and nul(l - lim))

# ---------------------------------------------------------------------------
# ed-t-premier-ordre-primitives (carte)
# ---------------------------------------------------------------------------
# Exemple : x^2 + sin x + 3 ; e^{x^2} ; ln(x^2+1).
assert primitive(x**2 + sin(x) + 3, 2 * x + cos(x), 0, 3)
assert nul(diff(exp(x**2), x) - 2 * x * exp(x**2)) and nul(diff(log(x**2 + 1), x) - 2 * x / (x**2 + 1))
# Contre-exemple (intervalle) : -1/x et -1/x + 1 ont la même dérivée 1/x^2.
assert nul(diff(-1 / x, x) - 1 / x**2) and nul(diff(-1 / x + 1, x) - 1 / x**2)
# Contre-exemple (Heaviside) : taux d'accroissement à gauche nul, alors que H(0) = 1 (raisonnement, rien à calculer).
# Affirmation 3 : x primitive de 1, négative en -1.
assert diff(x, x) == 1 and (x).subs(x, -1) < 0

# ed-t-premier-ordre-primitives-ex1 (et v1, v2).
cas = [
    (6 * x**2 + 1 / x**2, 0, 2 * x**3 - 1 / x - 1, [2 * x**3 + 1 / x - 3, 12 * x - 2 / x**3 - 10, 2 * x**3 - 1 / x]),
    (4 * x**3 - 3 / x**2, 2, x**4 + 3 / x - 2, [x**4 - 3 / x + 4, x**4 + 3 / x, 12 * x**2 + 6 / x**3 - 16]),
    (2 * x + 1 / sqrt(x), 0, x**2 + 2 * sqrt(x) - 3,
     [x**2 + sqrt(x) / 2 - Rational(3, 2), x**2 - 2 / sqrt(x) + 1, x**2 + 2 * sqrt(x)]),
]
for f, val, bon, leurres in cas:
    assert primitive(bon, f, 1, val)
    for l in leurres:
        assert not primitive(l, f, 1, val)

# ed-t-premier-ordre-primitives-ex2 (et v1, v2).
cas = [
    (x * exp(-x**2), 1, Rational(3, 2) - exp(-x**2) / 2,
     [exp(-x**2) / 2 + Rational(1, 2), -2 * exp(-x**2) + 3, x**2 / 2 * exp(-x**2) + 1]),
    (2 * x / (x**2 + 1), 3, log(x**2 + 1) + 3, [1 / (x**2 + 1) + 2, 2 * log(x**2 + 1) + 3, log(x**2 + 1)]),
    (cos(x) * exp(sin(x)), 2, exp(sin(x)) + 1, [exp(sin(x)) + 2, sin(x) * exp(sin(x)) + 2, exp(cos(x)) + 2 - E]),
]
for f, val, bon, leurres in cas:
    assert primitive(bon, f, 0, val)
    for l in leurres:
        assert not primitive(l, f, 0, val)

# ed-t-premier-ordre-primitives-ex3 (et v1, v2) : identification (ax+b)e^{±x}.
cas = [
    (x * exp(-x), -(x + 1) * exp(-x), 1 - (x + 1) * exp(-x),
     [(x - 1) * exp(-x) + 1, x**2 / 2 * exp(-x), 1 - (x + 1) * exp(x)]),
    ((2 * x + 1) * exp(x), (2 * x - 1) * exp(x), (2 * x - 1) * exp(x) + 1,
     [(2 * x + 1) * exp(x) - 1, (x**2 + x) * exp(x), (2 * x - 1) * exp(x)]),
    ((x - 2) * exp(-x), (1 - x) * exp(-x), (1 - x) * exp(-x) - 1,
     [(x - 1) * exp(-x) + 1, (1 - x) * exp(-x), (x**2 / 2 - 2 * x) * exp(-x)]),
]
for f, Fp, bon, leurres in cas:
    assert nul(diff(Fp, x) - f)
    assert primitive(bon, f, 0, 0)
    for l in leurres:
        assert not primitive(l, f, 0, 0)
assert limit(1 - (x + 1) * exp(-x), x, oo) == 1

# ed-t-premier-ordre-primitives-ex4 (et v1, v2) : une constante par intervalle sur R*.
cas = [  # (f, primitive négative, primitive positive, F(-1), F(1), bonne, leurres numériques)
    (1 / x**2, -1 / x, -1 / x, 0, 0, 1, [-1, 0]),
    (2 / x**3, -1 / x**2, -1 / x**2, 0, 2, 2, [0, -2]),
    (1 / x, log(-x), log(x), 1, 0, -1, [0, 1]),
]
for f, Gm, Gp, vm, vp, bon, leurres in cas:
    assert nul(diff(Gm, x) - f) and nul(diff(Gp, x) - f)
    C1 = vm - Gm.subs(x, -1)
    C2 = vp - Gp.subs(x, 1)
    ecart = (Gp.subs(x, 2) + C2) - (Gm.subs(x, -2) + C1)
    assert nul(ecart - bon)
    for l in leurres:
        assert not nul(ecart - l)
    # Il existe une solution (la formule ci-dessus) ; avec une seule constante on aurait C1 = C2, ce qui est faux ici.
    assert not nul(C1 - C2)

# ed-t-premier-ordre-primitives-ex5 : parité. Exemple x + 1, primitive de 1 (paire), non impaire.
assert (x + 1).subs(x, -1) != -(x + 1).subs(x, 1)
# Primitive de x (impaire) : x^2/2 + C paire pour tout C.
assert nul((x**2 / 2 + C).subs(x, -x) - (x**2 / 2 + C))

# Indices des variantes de primitives-ex1 : -3/x primitive de 3/x^2, 2 sqrt(x) primitive de 1/sqrt(x).
assert nul(diff(-3 / x, x) - 3 / x**2) and nul(diff(2 * sqrt(x), x) - 1 / sqrt(x))
# Cours des primitives : u'u^n -> u^{n+1}/(n+1) (ex. u = x^2 + 1, n = -2, u ne s'annule pas).
assert nul(diff((x**2 + 1) ** (-1) / (-1), x) - 2 * x * (x**2 + 1) ** (-2))

print("ed_t_premier_ordre : toutes les vérifications passent.")
