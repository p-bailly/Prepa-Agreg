"""Vérifications calculatoires du chapitre « int-l2-parametre » (intégration, L2) :
cartes, exercices et QCM de content/*/integration/l2-int-l2-parametre.json et content/qcm/int-l2-parametre.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Les intégrales impropres sont calculées numériquement à 30 chiffres (mpmath) ; les identités en SymPy.
Lancer : python3 verification/int_l2_parametre.py (doit finir sans erreur).
"""

from fractions import Fraction as Fr

import mpmath as mp
from sympy import (Rational, atan, diff, exp, integrate, limit, log, oo, pi, simplify, sin, cos, sqrt, symbols,
                   apart, factorial, gamma, Symbol)

mp.mp.dps = 30
x, t, u, p, xi = symbols("x t u p xi", real=True)
xp = symbols("xp", positive=True)
EPS = mp.mpf(10) ** -20


def q(f, a=0, b=mp.inf):
    """Intégrale numérique de f sur [a, b] (ou sur la subdivision a, si a est une liste)."""
    return mp.quad(f, a if isinstance(a, list) else [a, b])


def proche(a, b, tol=EPS):
    return abs(mp.mpf(a) - mp.mpf(b)) < tol if not isinstance(a, mp.mpc) and not isinstance(b, mp.mpc) else abs(a - b) < tol


def distincts_et_faux(vraie, leurres, tol=mp.mpf(10) ** -6):
    """La bonne valeur diffère de chaque leurre (en valeur numérique)."""
    for l in leurres:
        assert abs(mp.mpc(vraie) - mp.mpc(l)) > tol, (vraie, l)


# ===========================================================================
# int-l2-parametre-segment (carte)
# ===========================================================================
# Contre-exemple 1 : f(x,t) = 2x²t/(x²+t²)², F(x) = 1/(1+x²) pour x ≠ 0, F(0) = 0 ; f(x,x) = 1/(2x).
X = symbols("X", positive=True)
assert simplify(integrate(2 * X**2 * t / (X**2 + t**2) ** 2, (t, 0, 1)) - 1 / (1 + X**2)) == 0
assert simplify((2 * X**2 * t / (X**2 + t**2) ** 2).subs(t, X) - 1 / (2 * X)) == 0
# Contre-exemple 2 : f(x,t) = t e^{-t/√x} (x > 0) ; F(x) = x(1 - (1 + 1/√x)e^{-1/√x}) et F(x)/x → 1 ; ∂f/∂x(0,t) = 0.
Fx = integrate(t * exp(-t / sqrt(X)), (t, 0, 1))
assert simplify(Fx - X * (1 - (1 + 1 / sqrt(X)) * exp(-1 / sqrt(X)))) == 0
assert limit(Fx / X, X, 0, "+") == 1
T = symbols("T", positive=True)
assert limit(T * exp(-T / sqrt(X)) / X, X, 0, "+") == 0
# ∂f/∂x pour x > 0 : t² e^{-t/√x} / (2 x^{3/2}), de limite 0 quand x → 0+ (t > 0 fixé) : continuité séparée de ∂f/∂x.
assert simplify(diff(T * exp(-T / sqrt(X)), X) - T**2 * exp(-T / sqrt(X)) / (2 * X ** Rational(3, 2))) == 0
assert limit(T**2 * exp(-T / sqrt(X)) / (2 * X ** Rational(3, 2)), X, 0, "+") == 0
assert max(s * mp.e ** -s for s in mp.linspace(0, 5, 501)) <= 1 / mp.e + EPS
# Exemple : ∫_0^1 t^x dt = 1/(x+1), ∫_0^1 t^x ln t dt = -1/(x+1)², ∫_0^1 t ln t = -1/4.
assert simplify(integrate(t**X, (t, 0, 1)) - 1 / (X + 1)) == 0
assert simplify(integrate(t**X * log(t), (t, 0, 1)) + 1 / (X + 1) ** 2) == 0
assert integrate(t * log(t), (t, 0, 1)) == Rational(-1, 4)
# Trou 3 : avec f = xt sur [0,1], F' = 1/2 ; les leurres donnent x, 1, 0.
assert integrate(x * t, (t, 0, 1)).diff(x) == Rational(1, 2)
assert integrate(diff(x * t, t), (t, 0, 1)) == x and (diff(x * t, x).subs(t, 1) - diff(x * t, x).subs(t, 0)) == 1
assert diff(x * t, x).subs(t, 0) == 0
# Trou 1, leurre « bornée » : f = 1 si x > 0, 0 sinon, F = (b-a)·1_{x>0} discontinue (évident).
# Affirmation 3 : G(x) = (e^{x²}-1)/x, G'(0) = 1 alors que ∫_0^0 t e^{0} = 0.
G = (exp(x**2) - 1) / x
assert limit(G / x, x, 0) == 1 and simplify(integrate(exp(X * t), (t, 0, X)) - (exp(X**2) - 1) / X) == 0

# int-l2-parametre-segment-ex1 (+ v1, v2) : F'(0) = -1/2, -1/3, 1/4.
for f, val, leurres in [(1 / (1 + x * t), Rational(-1, 2), [Rational(1, 2), -1, Rational(-1, 3)]),
                        (exp(-x * t**2), Rational(-1, 3), [Rational(1, 3), Rational(-1, 2), -1]),
                        (sqrt(1 + x * t), Rational(1, 4), [Rational(-1, 4), Rational(1, 2), Rational(1, 3)])]:
    assert integrate(diff(f, x).subs(x, 0), (t, 0, 1)) == val and val not in leurres
assert simplify(log(1 + X) / X - integrate(1 / (1 + X * t), (t, 0, 1))) == 0

# int-l2-parametre-segment-ex2 (+ v1, v2) : F'(x) = ln(1+x²)/(2x²) ; F'(1), F'(2), F'(1/2).
Fat = atan(X) - log(1 + X**2) / (2 * X)
assert simplify(diff(Fat, X) - log(1 + X**2) / (2 * X**2)) == 0
assert simplify(integrate(atan(X * t), (t, 0, 1)) - Fat) == 0
for xv, val, leurres in [(1, mp.log(2) / 2, [mp.log(2), mp.pi / 4, mp.log(2) / 4]),
                         (2, mp.log(5) / 8, [mp.log(5) / 4, mp.log(5) / 16, mp.mpf(1) / 5]),
                         (mp.mpf(1) / 2, 2 * mp.log(mp.mpf(5) / 4), [mp.log(mp.mpf(5) / 4), mp.log(mp.mpf(5) / 4) / 2, mp.mpf(4) / 5])]:
    assert proche(q(lambda s: s / (1 + xv**2 * s**2), 0, 1), val)
    distincts_et_faux(val, leurres)

# int-l2-parametre-segment-ex3 (+ v1, v2) : ∫_0^1 (t^a - 1)/ln t = ln(1+a).
Fr_ = lambda a: q(lambda s: (s**a - 1) / mp.log(s), 0, 1)
for a in [2, 3, 4, mp.mpf("0.5")]:
    assert proche(Fr_(a), mp.log(1 + a))
assert proche(q(lambda s: (s**3 - s) / mp.log(s), 0, 1), mp.log(2))  # v1
distincts_et_faux(mp.log(3), [mp.log(2), 2, 2 * mp.log(3)])
distincts_et_faux(mp.log(2), [mp.log(3), mp.log(4), mp.log(1.5)])
distincts_et_faux(mp.log(5), [mp.log(4), 4, mp.log(5) / 4])

# int-l2-parametre-segment-ex4 : y = ∫_0^x sin(x-t) g(t) dt vérifie y'' + y = g, y(0) = y'(0) = 0 (g = e^t, g = 1).
for g in [exp(t), 1 + 0 * t, t**2]:
    y = integrate(sin(x - t) * g, (t, 0, x))
    assert simplify(diff(y, x, 2) + y - g.subs(t, x)) == 0 and y.subs(x, 0) == 0 and diff(y, x).subs(x, 0) == 0
    assert simplify(diff(y, x) - integrate(cos(x - t) * g, (t, 0, x))) == 0
# Leurres avec g = 1 : y' = sin x ; g(x) + sin x, -sin x, sin x - sin x = 0 sont faux.
y1 = 1 - cos(x)
assert simplify(integrate(sin(x - t), (t, 0, x)) - y1) == 0
assert all(simplify(l - sin(x)) != 0 for l in [1 + sin(x), -sin(x), sin(x) - sin(x)])

# ===========================================================================
# int-l2-parametre-continuite (carte)
# ===========================================================================
# Contre-exemple : ∫_0^∞ x e^{-xt} dt = 1 (x > 0) ; sup_x x e^{-xt} = 1/(et), atteint en x = 1/t.
assert integrate(X * exp(-X * T), (T, 0, oo)) == 1
assert simplify((X * exp(-X * T)).subs(X, 1 / T) - 1 / (exp(1) * T)) == 0
assert diff(X * exp(-X * T), X).subs(X, 1 / T).simplify() == 0
# Trou 2 : -x e^{-xt} ≤ 0 ≤ e^{-t} et ∫ = -1 ; |1/(1+t)|² intégrable, 1/(1+t) non ; |∫ x e^{-xt}| = 1 = ∫ e^{-t}.
assert integrate(1 / (1 + T) ** 2, (T, 0, oo)) == 1 and integrate(1 / (1 + T), (T, 0, oo)) == oo
assert integrate(exp(-T), (T, 0, oo)) == 1
# Exemple : ∫_0^∞ cos(xt)/(1+t²) dt = (π/2) e^{-|x|} (x = 1 et x = 2, numériquement).
for xv in [1, 2]:
    assert proche(mp.quadosc(lambda s: mp.cos(xv * s) / (1 + s**2), [0, mp.inf], omega=xv), mp.pi / 2 * mp.e ** -xv, mp.mpf(10) ** -12)
# Affirmation 3 : ∫_0^∞ e^{-xt} dt = 1/x.
assert integrate(exp(-X * T), (T, 0, oo)) == 1 / X

# int-l2-parametre-continuite-ex1 (+ v1, v2) : π/2, 1, 2.
assert integrate(1 / (1 + T**2), (T, 0, oo)) == pi / 2
assert integrate(T**2 * exp(-T), (T, 0, oo)) == 2
distincts_et_faux(mp.pi / 2, [0, 1])
distincts_et_faux(1, [0, 0.5, mp.pi / 2])
distincts_et_faux(2, [1, 0, 6])

# int-l2-parametre-continuite-ex2 (+ v1, v2) : domaines ]1/2,∞[, ]1,∞[, ]1/3,∞[ ; F(1) = π/2, F(2) = 1 (v1).
assert limit(T / sqrt(1 + T**2), T, oo) == 1          # x = 1/2 : équivalent à 1/t, non intégrable
assert limit(T / (1 + T), T, oo) == 1                 # v1, x = 1
assert limit(T / (1 + T**3) ** Rational(1, 3), T, oo) == 1  # v2, x = 1/3
assert proche(q(lambda s: (1 + s**2) ** -mp.mpf(0.75), [0, 1, 10, 100, 10**4, mp.inf]), mp.sqrt(mp.pi) * mp.gamma(0.25) / (2 * mp.gamma(0.75)), 1e-8)  # converge pour x = 3/4
assert integrate(1 / (1 + T) ** 2, (T, 0, oo)) == 1 and simplify(integrate((1 + T) ** -3, (T, 0, oo)) - Rational(1, 2)) == 0

# int-l2-parametre-continuite-ex3 (+ v1, v2) : limites 1, 1/2, π/4.
for xv in [3, 5, 50]:
    assert proche(q(lambda s: 1 / (1 + s**xv), 0, 1) + q(lambda s: 1 / (1 + s**xv), 1, mp.inf), (mp.pi / xv) / mp.sin(mp.pi / xv), mp.mpf(10) ** -15)
    assert proche(q(lambda s: s / (1 + s**xv), 0, 1) + q(lambda s: s / (1 + s**xv), 1, mp.inf), (mp.pi / xv) / mp.sin(2 * mp.pi / xv), mp.mpf(10) ** -15)
assert limit((pi / X) / sin(pi / X), X, oo) == 1 and limit((pi / X) / sin(2 * pi / X), X, oo) == Rational(1, 2)
v2_400 = q(lambda s: 1 / ((1 + s**400) * (1 + s**2)), 0, 1) + q(lambda s: 1 / ((1 + s**400) * (1 + s**2)), 1, mp.inf)
assert abs(v2_400 - mp.pi / 4) < 0.01
distincts_et_faux(1, [0, 0.5, 2]); distincts_et_faux(0.5, [1, 0, mp.mpf(1) / 3]); distincts_et_faux(mp.pi / 4, [mp.pi / 2, mp.pi / 8, 1])
# Dominations : pour x ≥ 2 et t ≥ 1, 1/(1+t^x) ≤ 1/t² ; pour x ≥ 4, t/(1+t^x) ≤ 1/t³.
assert all(1 / (1 + s**xv) <= 1 / s**2 and s / (1 + s**(xv + 2)) <= 1 / s**3 for xv in [2, 3, 7] for s in mp.linspace(1, 50, 200))

# int-l2-parametre-continuite-ex4 (+ v1, v2) : F(x) = e^x E1(x) = -ln x - γ + o(1) ; v1 : √x F(x) → π/2 ; v2 : F + (ln x)/2 bornée.
for xv in [mp.mpf(10) ** -3, mp.mpf(10) ** -6, mp.mpf(10) ** -9]:
    F = q(lambda s: mp.e ** -s / (xv + s), 0, 1) + q(lambda s: mp.e ** -s / (xv + s), 1, mp.inf)
    assert proche(F, mp.e ** xv * mp.e1(xv), mp.mpf(10) ** -15)
    assert abs(F - mp.log(1 + 1 / xv)) <= 1 + 1 / mp.e
    F1 = q(lambda s: mp.e ** -s / (xv + s**2), [0, mp.sqrt(xv), 1, mp.inf])
    assert abs(F1 - mp.pi / (2 * mp.sqrt(xv))) <= mp.log(1 + 1 / xv) / 2 + 1
    y = mp.sqrt(xv)
    F2 = q(lambda s: mp.e ** -s / (y + s), [0, y, 1, mp.inf])
    assert abs(F2 - mp.log(1 + 1 / y)) <= 1 + 1 / mp.e
xs = mp.mpf(10) ** -12
assert abs(mp.sqrt(xs) * q(lambda s: mp.e ** -s / (xs + s**2), [0, mp.sqrt(xs), 1, mp.inf]) - mp.pi / 2) < 10**-3
# Les leurres ne sont pas équivalents : -ln x / (1/x) → 0, etc. (rapports ≠ 1 en 0+).
assert limit(-log(X) * X, X, 0, "+") == 0 and limit(-log(X) * sqrt(X), X, 0, "+") == 0
assert limit((-log(X)) / (-log(X) / 2), X, 0, "+") == 2 and limit((pi / (2 * sqrt(X))) / (pi / (2 * X)), X, 0, "+") == 0

# ===========================================================================
# int-l2-parametre-derivation (carte)
# ===========================================================================
# Exemple : ∫_0^∞ e^{-t} sin(xt)/t dt = arctan x (x = 2) ; ∫ e^{-t} cos(xt) = 1/(1+x²).
assert proche(q(lambda s: mp.e ** -s * mp.sin(2 * s) / s), mp.atan(2))
assert simplify(integrate(exp(-T) * cos(X * T), (T, 0, oo)) - 1 / (1 + X**2)) == 0
# Contre-exemple : ∫_0^∞ x²/(x²+t²) dt = π|x|/2 ; domination 1/(1+t²) pour |x| < 1 ; ∂f/∂x(t,t) = 1/(2t).
assert simplify(integrate(X**2 / (X**2 + T**2), (T, 0, oo)) - pi * X / 2) == 0
assert simplify((1 / (1 + T**2) - X**2 / (X**2 + T**2)) * (1 + T**2) * (X**2 + T**2) - T**2 * (1 - X**2)) == 0
assert simplify(diff(X**2 / (X**2 + T**2), X).subs(X, T) - 1 / (2 * T)) == 0
# Leurre |x ∂f/∂x| ≤ h : 2x²t²/(x²+t²)² ≤ 1/2, et pour t ≥ 1, |x| < 1 : croissante en u = x² sur [0,1] (max 2t²/(1+t²)²).
U = symbols("U", positive=True)
gU = 2 * U * T**2 / (U + T**2) ** 2
assert simplify(Rational(1, 2) - gU - (U - T**2) ** 2 / (2 * (U + T**2) ** 2)) == 0
assert simplify(diff(gU, U) - 2 * T**2 * (T**2 - U) / (U + T**2) ** 3) == 0
assert integrate(2 * T**2 / (1 + T**2) ** 2, (T, 1, oo)) < oo
# Leurre trou 3 : f = -x e^{-t} : F' = -1, ∫|∂f| = 1, ∫ ∂_t f = x, ∫ h = 2 avec h = 2e^{-t}.
assert integrate(diff(-X * exp(-T), T), (T, 0, oo)) == X and integrate(2 * exp(-T), (T, 0, oo)) == 2
# Affirmation 3 : ∫_0^∞ t e^{-xt} dt = 1/x².
assert integrate(T * exp(-X * T), (T, 0, oo)) == 1 / X**2

# int-l2-parametre-derivation-ex1 (+ v1, v2) : ∫ t^n e^{-at} = n!/a^{n+1} : 3/8, 2/27, 24.
for n_, a, val, leurres in [(3, 2, Rational(3, 8), [Rational(3, 16), Rational(1, 8), Rational(3, 4)]),
                            (2, 3, Rational(2, 27), [Rational(2, 9), Rational(1, 27), Rational(2, 81)]),
                            (4, 1, 24, [6, 120, 4])]:
    assert integrate(T**n_ * exp(-a * T), (T, 0, oo)) == val == factorial(n_) / Rational(a) ** (n_ + 1) and val not in leurres
    assert simplify(diff(1 / X, X, n_) - (-1) ** n_ * factorial(n_) / X ** (n_ + 1)) == 0

# int-l2-parametre-derivation-ex2 (+ v1, v2) : F(x) = ½ ln(1 + x²/a²) pour e^{-at}.
for a, xv, val, leurres in [(1, 1, mp.log(2) / 2, [mp.log(2), mp.pi / 4, 1 - 1 / mp.sqrt(2)]),
                            (2, 1, mp.log(1.25) / 2, [mp.log(2) / 2, mp.log(1.25), mp.atan(0.5)]),
                            (1, mp.sqrt(3), mp.log(2), [mp.log(2) / 2, mp.pi / 3, mp.log(4)])]:
    num = q(lambda s: (1 - mp.cos(xv * s)) / s * mp.e ** (-a * s), [0, 1, 10, mp.inf])
    assert proche(num, val, mp.mpf(10) ** -15) and proche(val, mp.log(1 + xv**2 / a**2) / 2)
    distincts_et_faux(val, leurres)
assert simplify(integrate(sin(X * T) * exp(-2 * T), (T, 0, oo)) - X / (4 + X**2)) == 0

# int-l2-parametre-derivation-ex3 : F'' + F = 1/x (x = 1 et x = 2, numériquement) ; leurres faux.
for xv in [1, 2]:
    F = q(lambda s: mp.e ** (-xv * s) / (1 + s**2))
    F2 = q(lambda s: s**2 * mp.e ** (-xv * s) / (1 + s**2))
    assert proche(F2 + F, mp.mpf(1) / xv)
    assert abs(F2 - F - mp.mpf(1) / xv) > 0.01 and abs(F2 + F) > 0.01
assert abs(q(lambda s: s**2 * mp.e ** (-2 * s) / (1 + s**2)) + q(lambda s: mp.e ** (-2 * s) / (1 + s**2)) - mp.mpf(1) / 4) > 0.01

# int-l2-parametre-derivation-ex4 (+ v1, v2) : F(x) = (π/2) ln(1+x).
assert simplify(1 / ((1 + X**2 * T**2) * (1 + T**2)) - (1 / (1 - X**2)) * (1 / (1 + T**2) - X**2 / (1 + X**2 * T**2))) == 0
assert simplify((1 / (1 - X**2)) * (pi / 2 - pi * X / 2) - pi / (2 * (1 + X))) == 0
for xv, val, leurres in [(1, mp.pi / 2 * mp.log(2), [mp.log(2), mp.pi / 4 * mp.log(2), mp.pi / 2]),
                         (2, mp.pi / 2 * mp.log(3), [mp.pi / 2 * mp.log(2), mp.pi * mp.log(3), mp.pi / 4 * mp.log(3)]),
                         (mp.mpf(1) / 2, mp.pi / 2 * mp.log(1.5), [mp.pi / 4 * mp.log(1.5), mp.pi / 2 * mp.log(2), mp.pi / 4 * mp.log(3)])]:
    assert proche(q(lambda s: mp.atan(xv * s) / (s * (1 + s**2)), [0, 1, mp.inf]), val)
    distincts_et_faux(val, leurres)

# ===========================================================================
# int-l2-parametre-gamma (carte et exercices)
# ===========================================================================
assert [gamma(k) for k in range(1, 6)] == [1, 1, 2, 6, 24]
assert integrate(T**3 * exp(-2 * T), (T, 0, oo)) == Rational(3, 8)          # exemple
assert proche(mp.diff(mp.gamma, 1), -mp.euler)                                # Γ'(1) = -γ
assert proche(q(lambda s: mp.e ** -s * mp.log(s), [0, 1, mp.inf]), -mp.euler)
assert all(mp.e ** -s / s >= mp.e ** -1 / s for s in mp.linspace(0.01, 1, 50))  # contre-exemple x = 0
# Trous : Γ(x) = xΓ(x+1) faux en x = 2 ; Γ(2x) = 2xΓ(x) faux en x = 1 ; Γ(2) = 1 ≠ 2! ; 1^1 e^{-1} ≠ 1.
assert gamma(2) != 2 * gamma(3) and gamma(2) != 2 * gamma(1) and gamma(2) != 2 * 1 and exp(-1) != 1
assert proche(q(lambda s: s ** -0.5 * mp.e ** -s), mp.sqrt(mp.pi))           # converge pour x = 1/2 (leurres x > 1, x ≥ 1)
# Affirmation 2 / ex4 : minimum en α ∈ ]1,2[, Γ'(1) < 0 < Γ'(2) = 1 + Γ'(1), Γ(α) ≈ 0,886.
alpha = mp.findroot(lambda s: mp.diff(mp.gamma, s), 1.46)
assert 1 < alpha < 2 and abs(alpha - mp.mpf("1.46163")) < 1e-4 and abs(mp.gamma(alpha) - mp.mpf("0.8856")) < 1e-3
assert proche(mp.diff(mp.gamma, 2), 1 + mp.diff(mp.gamma, 1), mp.mpf(10) ** -15) and mp.diff(mp.gamma, 1) < 0
assert limit(X * gamma(X), X, 0, "+") == 1

# int-l2-parametre-gamma-ex1 (+ v1, v2) : ∫_0^1 (ln u)^4 = 24, (ln u)^3 = -6, u (ln u)² = 1/4.
assert integrate(log(U) ** 4, (U, 0, 1)) == 24 and integrate(log(U) ** 3, (U, 0, 1)) == -6
assert integrate(U * log(U) ** 2, (U, 0, 1)) == Rational(1, 4)
assert 24 not in [-24, 4, 120] and -6 not in [6, -2, -24] and Rational(1, 4) not in [2, Rational(1, 2), Rational(1, 8)]

# int-l2-parametre-gamma-ex2 (+ v1, v2) : ∫_0^∞ u^n e^{-u²} = ½Γ((n+1)/2) : 1, 1/2, 3.
for n_, val, leurres in [(5, 1, [2, Rational(1, 2), 60]), (3, Rational(1, 2), [1, 3, Rational(1, 4)]), (7, 3, [6, 12, 24])]:
    assert integrate(U**n_ * exp(-U**2), (U, 0, oo)) == val == gamma(Rational(n_ + 1, 2)) / 2 and val not in leurres

# int-l2-parametre-gamma-ex3 (+ v1, v2) : Γ(n + 1/2) = (2n)! √π / (4^n n!).
for n_ in range(6):
    assert simplify(gamma(n_ + Rational(1, 2)) - factorial(2 * n_) * sqrt(pi) / (4**n_ * factorial(n_))) == 0
sp = sqrt(pi)
for arg, val, leurres in [(Rational(7, 2), 15 * sp / 8, [15 * sp / 4, 105 * sp / 16, 3 * sp / 4]),
                          (Rational(9, 2), 105 * sp / 16, [105 * sp / 8, 15 * sp / 8, 945 * sp / 32]),
                          (Rational(5, 2), 3 * sp / 4, [3 * sp / 2, sp / 2, 15 * sp / 8])]:
    assert simplify(gamma(arg) - val) == 0 and all(simplify(l - val) != 0 for l in leurres)

# ===========================================================================
# int-l2-parametre-laplace (carte et exercices)
# ===========================================================================
P_ = symbols("P", positive=True)
assert simplify(integrate(T**3 * exp(-P_ * T), (T, 0, oo)) - 6 / P_**4) == 0           # exemple L(t^n)
# Contre-exemple 1 : L(1/(1+t))(p) ≥ e^{-1} ln(1 + 1/p) ; ex. p = 10^-3.
pv = mp.mpf(10) ** -3
assert q(lambda s: mp.e ** (-pv * s) / (1 + s), [0, 1, 1 / pv, mp.inf]) >= mp.e ** -1 * mp.log(1 + 1 / pv)
# Contre-exemple 2 : ∫ t e^{-pt}/(1+t²) ≥ ln(1/p)/(2e) ; t/(1+t²) ≥ 1/(2t) pour t ≥ 1.
assert q(lambda s: s * mp.e ** (-pv * s) / (1 + s**2), [0, 1, 1 / pv, mp.inf]) >= mp.log(1 / pv) / (2 * mp.e)
# Trous : f = 1 : L1 = 1/p, (L1)' = -1/p² ; leurres 1/p², 0, -1 ; p·L1 = 1 ; L(e^{-t}) = 1/(p+1) → 0 ≠ 1.
assert simplify(integrate(exp(-P_ * T), (T, 0, oo)) - 1 / P_) == 0
assert simplify(integrate(exp(-T) * exp(-P_ * T), (T, 0, oo)) - 1 / (P_ + 1)) == 0
assert simplify(integrate(exp(T) * exp(-(P_ + 1) * T), (T, 0, oo)) - 1 / P_) == 0       # L(e^{at}) avec a = 1, p+1 > 1

pn = mp.mpf(3)   # valeur de p pour les contrôles numériques


def L(f, pv=pn):
    return q(lambda s: f(s) * mp.e ** (-pv * s), [0, 1, 10, mp.inf])


# int-l2-parametre-laplace-ex1 (+ v1, v2).
for f, val, leurres in [(lambda s: mp.sin(2 * s), 2 / (pn**2 + 4), [pn / (pn**2 + 4), 2 / (pn**2 - 4), 1 / (pn**2 + 4)]),
                        (lambda s: mp.cos(3 * s), pn / (pn**2 + 9), [3 / (pn**2 + 9), pn / (pn**2 - 9), pn / (pn**2 + 3)]),
                        (lambda s: mp.sinh(s), 1 / (pn**2 - 1), [pn / (pn**2 - 1), 1 / (pn**2 + 1), 1 / (pn - 1) ** 2])]:
    assert proche(L(f), val)
    distincts_et_faux(val, leurres)
# int-l2-parametre-laplace-ex2 (+ v1, v2).
for f, val, leurres in [(lambda s: s**2 * mp.e ** -s, 2 / (pn + 1) ** 3, [1 / (pn + 1) ** 3, -2 / (pn + 1) ** 3, 2 / (pn + 1) ** 2]),
                        (lambda s: s * mp.sin(s), 2 * pn / (pn**2 + 1) ** 2, [-2 * pn / (pn**2 + 1) ** 2, 1 / (pn**2 + 1) ** 2, 2 / (pn**2 + 1) ** 2]),
                        (lambda s: s * mp.cos(s), (pn**2 - 1) / (pn**2 + 1) ** 2, [(1 - pn**2) / (pn**2 + 1) ** 2, 2 * pn / (pn**2 + 1) ** 2, pn**2 / (pn**2 + 1) ** 2])]:
    assert proche(L(f), val)
    distincts_et_faux(val, leurres)
assert simplify(-diff(P_ / (P_**2 + 1), P_) - (P_**2 - 1) / (P_**2 + 1) ** 2) == 0
# int-l2-parametre-laplace-ex3 (+ v1, v2) : ∫ sin t / t e^{-pt} = π/2 - arctan p.
for pv, val, leurres in [(1, mp.pi / 4, [mp.pi / 2, mp.log(2), 0.5]),
                         (mp.sqrt(3), mp.pi / 6, [mp.pi / 3, mp.pi / 2, mp.sqrt(3) / 2]),
                         (1 / mp.sqrt(3), mp.pi / 3, [mp.pi / 6, mp.pi / 2, 0.5])]:
    assert proche(L(lambda s: mp.sin(s) / s if s else mp.mpf(1), pv), val) and proche(val, mp.pi / 2 - mp.atan(pv))
    distincts_et_faux(val, leurres)
# int-l2-parametre-laplace-ex4 (+ v1, v2) : p·Lf(p) → ℓ quand p → 0+, → f(0) quand p → +∞.
for f, ell, f0 in [(mp.atan, mp.pi / 2, 0), (lambda s: s / (1 + s), 1, 0), (lambda s: 2 - mp.e ** -s, 2, 1)]:
    small, big = mp.mpf(10) ** -5, mp.mpf(10) ** 5
    assert abs(small * q(lambda s: f(s) * mp.e ** (-small * s), [0, 1, 10, 100, 1000, 10**4, 10**5, 10**6, mp.inf]) - ell) < 1e-3
    assert abs(big * q(lambda s: f(s) * mp.e ** (-big * s), [0, 10**-5, 10**-3, 1, mp.inf]) - f0) < 1e-3
assert simplify(integrate((2 - exp(-T)) * exp(-P_ * T), (T, 0, oo)) - (2 / P_ - 1 / (P_ + 1))) == 0
# Leurres d'équivalents : rapports à la bonne réponse non tendant vers 1.
assert limit((1 / P_) / (pi / (2 * P_)), P_, 0, "+") == 2 / pi and limit((-log(P_)) / (1 / P_), P_, 0, "+") == 0

# ===========================================================================
# int-l2-parametre-fourier (carte et exercices)
# ===========================================================================


def TF(f, x_, a=-mp.inf, b=mp.inf, pts=None):
    pts = pts or [a, -1, 0, 1, b]
    return mp.quad(lambda s: f(s) * mp.exp(-1j * x_ * s), pts)


z = mp.mpf("1.3")
assert abs(TF(lambda s: mp.e ** -abs(s), z) - 2 / (1 + z**2)) < EPS                        # exemple
# Contre-exemple : ∫ (cos u - 1)/u² du = -π < 0 ; FT de 1/(1+t²) en ξ = 1 vaut π/e.
assert abs(2 * mp.quadosc(lambda s: (mp.cos(s) - 1) / s**2, [0, mp.inf], omega=1) + mp.pi) < 1e-12
assert abs(2 * mp.quadosc(lambda s: mp.cos(s) / (1 + s**2), [0, mp.inf], omega=1) - mp.pi / mp.e) < 1e-12
# Trou 2 : indicatrice de [-1,1] : FT(0) = 2 > 1 = sup|f| ; f impaire sign(t)e^{-|t|} : FT = -2iξ/(1+ξ²), FT(0) = 0.
assert abs(TF(lambda s: mp.sign(s) * mp.e ** -abs(s), z) - (-2j * z / (1 + z**2))) < EPS
# int-l2-parametre-fourier-ex1 (+ v1, v2).
assert abs(TF(lambda s: 1, mp.pi / 2, -1, 1, [-1, 1]) - 4 / mp.pi) < EPS
assert abs(TF(lambda s: 1, mp.pi / 4, -2, 2, [-2, 2]) - 8 / mp.pi) < EPS
assert abs(TF(lambda s: 1, mp.pi, 0, 1, [0, 1]) - (-2j / mp.pi)) < EPS
distincts_et_faux(4 / mp.pi, [2 / mp.pi, 4 / mp.pi**2, 2]); distincts_et_faux(8 / mp.pi, [4 / mp.pi, 4 * mp.sqrt(2) / mp.pi, 16 / mp.pi])
distincts_et_faux(-2j / mp.pi, [2j / mp.pi, 2 / mp.pi, 0])
# int-l2-parametre-fourier-ex2 (+ v1, v2) : 4/(4+ξ²), 6/(9+ξ²), 1/(1+iξ) (contrôle en ξ = 1,3).
assert abs(TF(lambda s: mp.e ** (-2 * abs(s)), z) - 4 / (4 + z**2)) < EPS
assert abs(TF(lambda s: mp.e ** (-3 * abs(s)), z) - 6 / (9 + z**2)) < EPS
assert abs(TF(lambda s: mp.e ** -s if s >= 0 else 0, z) - 1 / (1 + 1j * z)) < EPS
distincts_et_faux(4 / (4 + z**2), [2 / (4 + z**2), 4 / (2 + z**2), 1 / (4 + z**2)])
distincts_et_faux(6 / (9 + z**2), [3 / (9 + z**2), 6 / (3 + z**2), 2 / (9 + z**2)])
distincts_et_faux(1 / (1 + 1j * z), [1 / (1 - 1j * z), 2 / (1 + z**2), 1 / (1 + z**2)])
# int-l2-parametre-fourier-ex3 (+ v1, v2).
g2 = 2 / (1 + xi**2)
assert simplify(diff(g2, xi, 2) - (12 * xi**2 - 4) / (1 + xi**2) ** 3) == 0
for f, val, leurres in [(lambda s: s * mp.e ** -abs(s), -4j * z / (1 + z**2) ** 2, [4j * z / (1 + z**2) ** 2, -4 * z / (1 + z**2) ** 2, 2j / (1 + z**2) ** 2]),
                        (lambda s: s * mp.e ** (-2 * abs(s)), -8j * z / (4 + z**2) ** 2, [8j * z / (4 + z**2) ** 2, -8 * z / (4 + z**2) ** 2, -4j * z / (4 + z**2) ** 2]),
                        (lambda s: s**2 * mp.e ** -abs(s), (4 - 12 * z**2) / (1 + z**2) ** 3, [(12 * z**2 - 4) / (1 + z**2) ** 3, 4 / (1 + z**2) ** 3, -4 / (1 + z**2) ** 2])]:
    assert abs(TF(f, z) - val) < EPS
    distincts_et_faux(val, leurres)
# int-l2-parametre-fourier-ex4 (+ v1, v2) : FT de e^{-at²} = √(π/a) e^{-ξ²/(4a)} (ξ = 1,3).
for a, val, leurres in [(mp.mpf(1) / 2, mp.sqrt(2 * mp.pi) * mp.e ** (-z**2 / 2), [mp.sqrt(2 * mp.pi) * mp.e ** (-z**2), mp.sqrt(mp.pi) * mp.e ** (-z**2 / 2), mp.sqrt(2 * mp.pi) * mp.e ** (-2 * z**2)]),
                        (1, mp.sqrt(mp.pi) * mp.e ** (-z**2 / 4), [mp.sqrt(mp.pi) * mp.e ** (-z**2), mp.sqrt(2 * mp.pi) * mp.e ** (-z**2 / 2), mp.sqrt(mp.pi) * mp.e ** (-z**2 / 2)]),
                        (2, mp.sqrt(mp.pi / 2) * mp.e ** (-z**2 / 8), [mp.sqrt(mp.pi / 2) * mp.e ** (-z**2 / 2), mp.sqrt(2 * mp.pi) * mp.e ** (-z**2 / 8), mp.sqrt(mp.pi) * mp.e ** (-z**2 / 8)])]:
    assert abs(TF(lambda s: mp.e ** (-a * s**2), z) - val) < EPS
    distincts_et_faux(val, leurres)

# ===========================================================================
# int-l2-parametre-gauss (carte et exercices)
# ===========================================================================
A = symbols("A", positive=True)
assert integrate(exp(-t**2), (t, -oo, oo)) == sqrt(pi) and simplify(integrate(exp(-A * t**2), (t, -oo, oo)) - sqrt(pi / A)) == 0
assert gamma(Rational(1, 2)) == sqrt(pi)
# Démonstration : f(x) + g(x) = π/4 et 0 ≤ g(x) ≤ e^{-x²} (x = 0,7 et x = 2).
for xv in [mp.mpf("0.7"), mp.mpf(2)]:
    fv = q(lambda s: mp.e ** (-s**2), 0, xv) ** 2
    gv = q(lambda s: mp.e ** (-xv**2 * (1 + s**2)) / (1 + s**2), 0, 1)
    assert proche(fv + gv, mp.pi / 4) and 0 <= gv <= mp.e ** (-xv**2)
# Trous et affirmations : leurres faux en a = 1 ou a = 4 ; ∫ e^{-4t²} = √π/2 ; ∫ t² e^{-t²} = √π/2 ; ∫ e^{-t²/2} = √(2π).
assert integrate(exp(-4 * t**2), (t, -oo, oo)) == sqrt(pi) / 2 != sqrt(pi) / 4
assert integrate(t**2 * exp(-t**2), (t, -oo, oo)) == sqrt(pi) / 2
assert simplify(integrate(exp(-t**2 / 2), (t, -oo, oo)) - sqrt(2 * pi)) == 0
assert sqrt(pi / 2) != sqrt(pi) and sqrt(pi) / 4 != sqrt(pi) / 2 and pi != sqrt(pi)
# int-l2-parametre-gauss-ex1 (+ v1, v2).
assert simplify(integrate(exp(-3 * t**2), (t, -oo, oo)) - sqrt(pi / 3)) == 0
assert integrate(exp(-t**2 / 4), (t, 0, oo)) == sqrt(pi)
assert simplify(integrate(exp(-t**2 / 8), (t, -oo, oo)) - 2 * sqrt(2 * pi)) == 0
distincts_et_faux(mp.sqrt(mp.pi / 3), [mp.sqrt(mp.pi) / 3, mp.sqrt(3 * mp.pi), mp.pi / 3])
distincts_et_faux(mp.sqrt(mp.pi), [mp.sqrt(mp.pi) / 2, 2 * mp.sqrt(mp.pi), mp.sqrt(2 * mp.pi)])
distincts_et_faux(2 * mp.sqrt(2 * mp.pi), [mp.sqrt(2 * mp.pi), mp.sqrt(mp.pi / 8), 8 * mp.sqrt(mp.pi)])
# int-l2-parametre-gauss-ex2 (+ v1, v2) : compléter le carré.
for e_, val in [(-t**2 + 2 * t, exp(1) * sqrt(pi)), (-t**2 - t, exp(Rational(1, 4)) * sqrt(pi)), (-2 * t**2 + 4 * t, exp(2) * sqrt(pi / 2))]:
    assert simplify(integrate(exp(e_), (t, -oo, oo)) - val) == 0
distincts_et_faux(mp.e * mp.sqrt(mp.pi), [mp.sqrt(mp.pi), mp.sqrt(mp.pi) / mp.e, mp.e**2 * mp.sqrt(mp.pi)])
distincts_et_faux(mp.e ** 0.25 * mp.sqrt(mp.pi), [mp.e ** -0.25 * mp.sqrt(mp.pi), mp.sqrt(mp.pi), mp.e ** 0.5 * mp.sqrt(mp.pi)])
distincts_et_faux(mp.e**2 * mp.sqrt(mp.pi / 2), [mp.e**2 * mp.sqrt(mp.pi), mp.e**4 * mp.sqrt(mp.pi / 2), mp.e * mp.sqrt(mp.pi / 2)])
# int-l2-parametre-gauss-ex3 (+ v1, v2) : M_2 = 3√(2π), M_3 = 15√(2π), N_1 = √π/4.
assert simplify(integrate(t**4 * exp(-t**2 / 2), (t, -oo, oo)) - 3 * sqrt(2 * pi)) == 0
assert simplify(integrate(t**6 * exp(-t**2 / 2), (t, -oo, oo)) - 15 * sqrt(2 * pi)) == 0
assert integrate(t**2 * exp(-t**2), (t, 0, oo)) == sqrt(pi) / 4
# int-l2-parametre-gauss-ex4 (+ v1, v2) : I(a) = (√π/2) e^{-2a} ; ∫ e^{-(t-1/t)²} = √π/2.
for a in [1, 2, mp.mpf("0.3")]:
    assert proche(q(lambda s: mp.e ** (-s**2 - a**2 / s**2), [0, 1, mp.inf]), mp.sqrt(mp.pi) / 2 * mp.e ** (-2 * a))
assert proche(q(lambda s: mp.e ** (-(s - 1 / s) ** 2), [0, 1, mp.inf]), mp.sqrt(mp.pi) / 2)
distincts_et_faux(mp.sqrt(mp.pi) / 2 * mp.e**-2, [mp.sqrt(mp.pi) / 2 * mp.e**-1, mp.sqrt(mp.pi) * mp.e**-2, mp.sqrt(mp.pi) / 2 * mp.e**-0.5])
distincts_et_faux(mp.sqrt(mp.pi) / 2 * mp.e**-4, [mp.sqrt(mp.pi) / 2 * mp.e**-2, mp.sqrt(mp.pi) * mp.e**-4, mp.sqrt(mp.pi) / 2 * mp.e**-8])
distincts_et_faux(mp.sqrt(mp.pi) / 2, [mp.sqrt(mp.pi), mp.sqrt(mp.pi) / 2 * mp.e**-2, mp.sqrt(mp.pi) / 2 * mp.e**2])

# ===========================================================================
# QCM int-l2-parametre-q1 … q8
# ===========================================================================
assert integrate(t**2, (t, 0, 1)) == Rational(1, 3)                                                   # q1
assert integrate(X * exp(-X * T), (T, 0, oo)) == 1                                                    # q2
assert proche(q(lambda s: (mp.e ** (-3 * s) - mp.e ** -s) / s, [0, 1, mp.inf]), -mp.log(3))          # q3
assert gamma(Rational(3, 2)) == sqrt(pi) / 2                                                          # q4
assert abs(TF(lambda s: mp.e ** -abs(s), 0) - 2) < EPS                                                # q6
assert simplify(integrate(exp(-t**2 / 2), (t, -oo, oo)) - sqrt(2 * pi)) == 0                          # q7

print("int_l2_parametre : toutes les vérifications passent.")
