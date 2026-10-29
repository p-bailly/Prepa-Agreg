"""Vérifications calculatoires de l'enrichissement du chapitre « ed-l1-lineaires » (L1).

Questions de content/qcm/ed-l1-lineaires.json et exercices ajoutés (ordre2-ex5, variation-constante-ex5,
exp-polynome-ex5) de content/exercices/equations-differentielles/l1-ed-l1-lineaires.json.
Chaque assert cite l'identifiant de la question ou de l'exercice (v1, v2 = variantes).
Lancer : python3 verification/enrichi_ed_l1_lineaires.py (doit finir sans erreur).
"""

from sympy import (Function, I, Rational, atan, cos, diff, dsolve, exp, expand, log, oo, limit, simplify, sin, solve,
                   symbols, roots, Eq)

x = symbols("x", positive=True)  # pour les équations d'Euler (x > 0)
X = symbols("X", real=True)
t = symbols("t", real=True)
al, be, lam, mu, C = symbols("alpha beta lambda mu C")


def nul(e):
    return simplify(e) == 0


# --- QCM ---------------------------------------------------------------------
# ed-l1-lineaires-q1 : y' + cos(x) y = 0 : C e^{−sin x} ; leurres faux (C = 1).
L1 = lambda y: diff(y, X) + cos(X) * y
assert nul(L1(C * exp(-sin(X))))
for y in (exp(sin(X)), exp(-cos(X)), cos(X)):
    assert not nul(L1(y))

# ed-l1-lineaires-q2 : y'' + 2y' + 2y = 0, racines −1 ± i.
L2 = lambda y: diff(y, X, 2) + 2 * diff(y, X) + 2 * y
r = symbols("r")
assert set(roots(r**2 + 2 * r + 2, r)) == {-1 + I, -1 - I}
assert nul(L2(exp(-X) * (lam * cos(X) + mu * sin(X))))
for y in (exp(X) * cos(X), exp(-X) * cos(2 * X), X * exp(-X)):  # chaque leurre contient une fonction non solution
    assert not nul(L2(y))

# ed-l1-lineaires-q3 (vrai) : y = λe^{2x} + μe^{−2x}, y(0) = 1, y'(0) = 2 ⇒ λ = 1, μ = 0.
y = lam * exp(2 * X) + mu * exp(-2 * X)
assert nul(diff(y, X, 2) - 4 * y)
assert solve([y.subs(X, 0) - 1, diff(y, X).subs(X, 0) - 2], [lam, mu]) == {lam: 1, mu: 0}

# ed-l1-lineaires-q4 (faux) : α x e^x n'est jamais solution de y'' + y = x e^x ; une solution est (x−1)e^x/2.
y = al * X * exp(X)
reste = expand(simplify((diff(y, X, 2) + y - X * exp(X)) * exp(-X)))
assert solve([reste.coeff(X, 1), reste.coeff(X, 0)], [al]) == []
assert nul(diff((X - 1) * exp(X) / 2, X, 2) + (X - 1) * exp(X) / 2 - X * exp(X))

# ed-l1-lineaires-q5 : y'' − 4y' + 4y = (x+1)e^{2x} : x²(αx+β)e^{2x} convient (α = 1/6, β = 1/2), pas les autres formes.
L5 = lambda y: diff(y, X, 2) - 4 * diff(y, X) + 4 * y
cible = (X + 1) * exp(2 * X)


def coeffs_possibles(forme):
    reste = expand(simplify((L5(forme) - cible) * exp(-2 * X)))
    eqs = [reste.coeff(X, k) for k in range(5)]
    return solve(eqs, [al, be], dict=True)


assert coeffs_possibles(X**2 * (al * X + be) * exp(2 * X)) == [{al: Rational(1, 6), be: Rational(1, 2)}]
for forme in ((al * X + be) * exp(2 * X), X * (al * X + be) * exp(2 * X), al * X**2 * exp(2 * X)):
    assert coeffs_possibles(forme) == []

# ed-l1-lineaires-q6 : y' + y = x + e^{−x} : x − 1 + x e^{−x} ; leurres faux.
L6 = lambda y: diff(y, X) + y - X - exp(-X)
assert nul(L6(X - 1 + X * exp(-X)))
for y in (X + 1 + X * exp(-X), X - 1 + exp(-X) / 2, X - 1 - X * exp(-X)):
    assert not nul(L6(y))

# ed-l1-lineaires-q7 (faux) : y = e^{−x} + (sin x − cos x)/2 est solution de y' + y = sin x et non bornée en −∞.
y = exp(-X) + (sin(X) - cos(X)) / 2
assert nul(diff(y, X) + y - sin(X))
assert limit(y, X, -oo) == oo

# ed-l1-lineaires-q8 : y' − 2y = e^{2x}/(1+x²), y(0) = 0 : e^{2x} arctan x ; leurres faux.
L8 = lambda y: diff(y, X) - 2 * y - exp(2 * X) / (1 + X**2)
assert nul(L8(exp(2 * X) * atan(X))) and atan(0) == 0
for y in (exp(-2 * X) * atan(X), exp(2 * X) * log(1 + X**2) / 2, X * exp(2 * X) / (1 + X**2)):
    assert y.subs(X, 0) == 0 and not nul(L8(y))

# --- ed-l1-lineaires-ordre2-ex5 : équation d'Euler x²y'' + p x y' + q y = 0 sur ]0,+∞[ ----------------------------
z = Function("z")
yf = Function("y")
# Q1-Q2 : si z(t) = y(e^t), alors z' = e^t y'(e^t), z'' = e^t y'(e^t) + e^{2t} y''(e^t), d'où
# x y'(x) = z'(t) et x² y''(x) = z''(t) − z'(t) en x = e^t.
zt = yf(exp(t))
z1, z2 = diff(zt, t), diff(zt, t, 2)
yp = diff(yf(x), x).subs(x, exp(t))
ypp = diff(yf(x), x, 2).subs(x, exp(t))
assert nul(z1 - exp(t) * yp) and nul(z2 - (exp(t) * yp + exp(2 * t) * ypp))
assert nul(exp(2 * t) * ypp - (z2 - z1))


def euler(p, q, general, cond_sol, leurres):
    E = lambda y: x**2 * diff(y, x, 2) + p * x * diff(y, x) + q * y
    assert nul(E(general))
    # (E') : z'' + (p − 1) z' + q z = 0 ; sa solution générale, composée avec ln, redonne la solution générale.
    sol = dsolve(z(t).diff(t, 2) + (p - 1) * z(t).diff(t) + q * z(t), z(t)).rhs
    assert nul(E(sol.subs(t, log(x))))
    assert nul(E(cond_sol)) and nul(cond_sol.subs(x, 1)) and nul(diff(cond_sol, x).subs(x, 1) - 1)
    for L in leurres:
        assert not (nul(E(L)) and nul(L.subs(x, 1)) and nul(diff(L, x).subs(x, 1) - 1))


# main : x²y'' − 2xy' + 2y = 0 ⇒ z'' − 3z' + 2z = 0 ⇒ y = λx + μx² ; y(1) = 0, y'(1) = 1 ⇒ y = x² − x.
euler(-2, 2, lam * x + mu * x**2, x**2 - x, [exp(2 * (x - 1)) - exp(x - 1), x * log(x), x - x**2])
# Le leurre e^{2(x−1)} − e^{x−1} est la solution de (E') (en t) avec ces conditions, pas de (E).
assert nul(diff(exp(2 * (t - 1)) - exp(t - 1), t, 2) - 3 * diff(exp(2 * (t - 1)) - exp(t - 1), t) + 2 * (exp(2 * (t - 1)) - exp(t - 1)))
# v1 : x²y'' + xy' − 4y = 0 ⇒ z'' − 4z = 0 ⇒ y = λx² + μx^{−2} ; y = x²/4 − 1/(4x²).
euler(1, -4, lam * x**2 + mu / x**2, x**2 / 4 - 1 / (4 * x**2),
      [(exp(2 * (x - 1)) - exp(-2 * (x - 1))) / 4, x * log(x), 1 / (4 * x**2) - x**2 / 4])
# v2 : x²y'' − xy' + y = 0 ⇒ z'' − 2z' + z = 0 ⇒ y = (λ + μ ln x) x ; y = x ln x.
euler(-1, 1, (lam + mu * log(x)) * x, x * log(x), [(x - 1) * exp(x - 1), x**2 - x, -x * log(x)])

# --- ed-l1-lineaires-variation-constante-ex5 : (1+x²)y' + 2xy = b(x) ------------------------------------------------
def vc(b, y0, attendu, leurres):
    E = lambda y: (1 + X**2) * diff(y, X) + 2 * X * y - b
    # e^{−A} = 1/(1+x²) avec A = ln(1+x²) primitive de 2x/(1+x²)
    assert nul(diff(log(1 + X**2), X) - 2 * X / (1 + X**2))
    assert nul(E(attendu)) and nul(attendu.subs(X, 0) - y0)
    for L in leurres:
        assert not (nul(E(L)) and nul(L.subs(X, 0) - y0))


# main : b = 3x² + 1, y(0) = 2 ⇒ C = x³ + x + 2, y = x + 2/(1+x²).
vc(3 * X**2 + 1, 2, X + 2 / (1 + X**2), [(X**3 + X) / (1 + X**2) + 2, X + 2 * exp(-X**2), (X**3 + X + 2) * (1 + X**2)])
assert nul((X**3 + X + 2) / (1 + X**2) - (X + 2 / (1 + X**2)))
# v1 : b = 4x³ + 2x, y(0) = −1 ⇒ y = x² − 1/(1+x²).
vc(4 * X**3 + 2 * X, -1, X**2 - 1 / (1 + X**2), [X**2 - 1, X**2 - exp(-X**2), (X**2 - 1) / (1 + X**2)])
assert nul((X**4 + X**2 - 1) / (1 + X**2) - (X**2 - 1 / (1 + X**2)))
# v2 : b = cos x, y(0) = 1 ⇒ y = (1 + sin x)/(1+x²).
vc(cos(X), 1, (1 + sin(X)) / (1 + X**2), [cos(X) / (1 + X**2), (1 + sin(X)) * exp(-X**2), (1 + sin(X)) / (1 + X**2) ** 2])

# --- ed-l1-lineaires-exp-polynome-ex5 : racine double ----------------------------------------------------------------
def double(m, c, y0, yp0, alpha, attendu, leurres):
    E = lambda y: diff(y, X, 2) - 2 * m * diff(y, X) + m**2 * y - c * exp(m * X)
    # y = z e^{mx} ⇒ y'' − 2my' + m²y = z'' e^{mx}
    zz = Function("zz")
    assert nul(diff(zz(X) * exp(m * X), X, 2) - 2 * m * diff(zz(X) * exp(m * X), X) + m**2 * zz(X) * exp(m * X)
               - diff(zz(X), X, 2) * exp(m * X))
    assert nul(E(alpha * X**2 * exp(m * X)))
    assert nul(E((lam + mu * X + alpha * X**2) * exp(m * X)))
    ok = lambda y: nul(E(y)) and nul(y.subs(X, 0) - y0) and nul(diff(y, X).subs(X, 0) - yp0)
    assert ok(attendu)
    for L in leurres:
        assert not ok(L)


# main : y'' − 2y' + y = 6e^x, α = 3, y(0) = 1, y'(0) = 0 ⇒ y = (3x² − x + 1)e^x.
double(1, 6, 1, 0, 3, (3 * X**2 - X + 1) * exp(X),
       [(3 * X**2 + X + 1) * exp(X), (3 * X**2 + 1) * exp(X), (X**2 - X + 1) * exp(X)])
# v1 : y'' + 4y' + 4y = 4e^{−2x}, α = 2, y(0) = 0, y'(0) = 1 ⇒ y = (2x² + x)e^{−2x}.
double(-2, 4, 0, 1, 2, (2 * X**2 + X) * exp(-2 * X),
       [(2 * X**2 - X) * exp(-2 * X), (X**2 + X) * exp(-2 * X), (2 * X**2 + X) * exp(2 * X)])
# v2 : y'' − 6y' + 9y = 2e^{3x}, α = 1, y(0) = 2, y'(0) = 5 ⇒ y = (x² − x + 2)e^{3x}.
double(3, 2, 2, 5, 1, (X**2 - X + 2) * exp(3 * X),
       [(X**2 + 5 * X + 2) * exp(3 * X), (X**2 - X + 2) * exp(-3 * X), (2 * X**2 - X + 2) * exp(3 * X)])

# --- Relecture indépendante (2026-10-29) : réponses retrouvées à l'aveugle, recalculées par dsolve -------------------
import sympy as sp
_x = sp.symbols("x", positive=True)
_y = sp.Function("y")
def _ok(expr, eq):
    return sp.simplify(eq(expr)) == 0
# ed-l1-lineaires-ordre2-ex5 (et v1, v2) : équations d'Euler, solution avec y(1) = 0, y'(1) = 1.
for p, q, att in ((-2, 2, _x**2 - _x), (1, -4, (_x**2 - _x**-2) / 4), (-1, 1, _x * sp.log(_x))):
    eq = lambda u, p=p, q=q: _x**2 * u.diff(_x, 2) + p * _x * u.diff(_x) + q * u
    assert _ok(att, eq) and att.subs(_x, 1) == 0 and sp.simplify(att.diff(_x).subs(_x, 1) - 1) == 0
# ed-l1-lineaires-variation-constante-ex5 (et v1, v2) : ((1+x²)y)' = b.
_X = sp.symbols("X", real=True)
for bb, y0, att in ((3 * _X**2 + 1, 2, _X + 2 / (1 + _X**2)), (4 * _X**3 + 2 * _X, -1, _X**2 - 1 / (1 + _X**2)),
                    (sp.cos(_X), 1, (1 + sp.sin(_X)) / (1 + _X**2))):
    assert sp.simplify((1 + _X**2) * att.diff(_X) + 2 * _X * att - bb) == 0 and att.subs(_X, 0) == y0
# ed-l1-lineaires-exp-polynome-ex5 (et v1, v2) : racine double, dsolve avec conditions initiales.
for m, c, y0, yp0, att in ((1, 6, 1, 0, (3 * _X**2 - _X + 1) * sp.exp(_X)), (-2, 4, 0, 1, (2 * _X**2 + _X) * sp.exp(-2 * _X)),
                           (3, 2, 2, 5, (_X**2 - _X + 2) * sp.exp(3 * _X))):
    eq = sp.Eq(_y(_X).diff(_X, 2) - 2 * m * _y(_X).diff(_X) + m**2 * _y(_X), c * sp.exp(m * _X))
    s_ = sp.dsolve(eq, _y(_X), ics={_y(0): y0, _y(_X).diff(_X).subs(_X, 0): yp0}).rhs
    assert sp.simplify(s_ - att) == 0

print("enrichi_ed_l1_lineaires : tout est vérifié.")
