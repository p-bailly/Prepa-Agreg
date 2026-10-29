"""Vérifications calculatoires de l'enrichissement du chapitre « ed-t-premier-ordre » (Terminale).

Questions de content/qcm/ed-t-premier-ordre.json et exercices ajoutés (ay-ex6, ay-b-ex6, primitives-ex6) de
content/exercices/equations-differentielles/terminale-ed-t-premier-ordre.json.
Chaque assert cite l'identifiant de la question ou de l'exercice (v1, v2 = variantes).
Lancer : python3 verification/enrichi_ed_t_premier_ordre.py (doit finir sans erreur).
"""

from sympy import (E, Rational, cos, diff, exp, limit, log, oo, pi, simplify, sin, solve, sqrt, symbols, tan, discriminant)

x = symbols("x", real=True)
a, b, C = symbols("a b C", real=True)


def nul(e):
    return simplify(e) == 0


def sol(y, membre, x0=None, y0=None):
    """y vérifie y' = membre(y) (et y(x0) = y0 si donné)."""
    ok = nul(diff(y, x) - membre(y))
    if x0 is not None:
        ok = ok and nul(y.subs(x, x0) - y0)
    return ok


# --- QCM ---------------------------------------------------------------------
# ed-t-premier-ordre-q1 : f = K e^{ax}, f(1) = 2, f(3) = 18 ⇒ e^{2a} = 9, f(2) = 2 e^{a} = 6 ; leurres 10, 9, 4 faux.
K = symbols("K", real=True)
s = solve([K * exp(a) - 2, K * exp(3 * a) - 18], [K, a], dict=True)
s = [d for d in s if d[a].is_real]
assert len(s) == 1 and nul(s[0][K] * exp(2 * s[0][a]) - 6)
assert all(v != 6 for v in (10, 9, 4))

# ed-t-premier-ordre-q2 : 2y' + y = 4, y(0) = 1 ⇒ y = 4 − 3e^{−x/2} ; leurres faux.
eq = lambda y: (4 - y) / 2
assert sol(4 - 3 * exp(-x / 2), eq, 0, 1)
for y in (4 - 3 * exp(-2 * x), 2 - exp(-x / 2), 4 - 3 * exp(x / 2)):
    assert not sol(y, eq, 0, 1)

# ed-t-premier-ordre-q3 (faux) : y = 3 + e^{−2x} est solution de y' = −2y + 6 et y(0) = 4 > 3.
assert sol(3 + exp(-2 * x), lambda y: -2 * y + 6) and (3 + exp(-2 * x)).subs(x, 0) == 4

# ed-t-premier-ordre-q4 : primitive de 2x/(x²+1)² nulle en 0 : x²/(x²+1) ; leurres faux.
f = 2 * x / (x**2 + 1) ** 2
assert nul(diff(x**2 / (x**2 + 1), x) - f) and (x**2 / (x**2 + 1)).subs(x, 0) == 0
for F in (log(x**2 + 1), -1 / (x**2 + 1), x**2 / (x**2 + 1) ** 2):
    assert not (nul(diff(F, x) - f) and F.subs(x, 0) == 0)
assert (-1 / (x**2 + 1)).subs(x, 0) == -1  # bonne dérivée, mauvaise valeur en 0

# ed-t-premier-ordre-q5 (vrai) : (x ln x − x)' = ln x.
assert nul(diff(x * log(x) - x, x) - log(x))

# ed-t-premier-ordre-q6 : N(t+5) = 2N(t) ⟺ e^{5a} = 2 ⟺ a = ln2/5 ; leurres faux.
assert nul(exp(5 * log(2) / 5) - 2)
for v in (5 * log(2), Rational(2, 5), log(5) / 2):
    assert not nul(exp(5 * v) - 2)

# ed-t-premier-ordre-q7 (vrai) : e^{2x+3} solution de y' = 2y.
assert sol(exp(2 * x + 3), lambda y: 2 * y)

# ed-t-premier-ordre-q8 : y' = y/2 − 1 : y = (y(0) − 2) e^{x/2} + 2, bornée sur [0,+∞[ ssi y(0) = 2.
y0 = symbols("y0", real=True)
yy = (y0 - 2) * exp(x / 2) + 2
assert sol(yy, lambda y: y / 2 - 1, 0, y0)
for v in (0, -2, Rational(1, 2)):
    assert limit(yy.subs(y0, v), x, oo) == -oo  # non bornée
assert limit(yy.subs(y0, 2), x, oo) == 2

# --- ed-t-premier-ordre-ay-ex6 -----------------------------------------------
def probleme(k, c, f_attendu, leurres):
    g = x**2 * exp(-k * x)
    membre = lambda y: -k * y + 2 * x * exp(-k * x)
    assert sol(g, membre)                                    # Q1
    # Q3 : solutions (x² + C) e^{−kx}
    assert sol((x**2 + C) * exp(-k * x), membre)
    assert sol(f_attendu, membre, 0, c)                      # Q4
    for L in leurres:
        assert not sol(L, membre, 0, c)
    assert limit(f_attendu, x, oo) == 0                      # Q5 : limite
    return simplify(diff(f_attendu, x) * exp(k * x))


# ay-ex6 : y' + 2y = 2x e^{−2x}, f(0) = 1, f = (x²+1)e^{−2x}, f' = −2(x² − x + 1)e^{−2x} < 0.
d = probleme(2, 1, (x**2 + 1) * exp(-2 * x),
             [x**2 * exp(-2 * x) + 1, (x**2 + 1) * exp(2 * x), (x + 1) ** 2 * exp(-2 * x)])
assert nul(d + 2 * (x**2 - x + 1)) and discriminant(x**2 - x + 1, x) == -3
# ay-ex6 v1 : y' + y = 2x e^{−x}, f(0) = 1, f = (x²+1)e^{−x}, f' = −(x−1)² e^{−x}.
d = probleme(1, 1, (x**2 + 1) * exp(-x),
             [x**2 * exp(-x) + 1, (x**2 + 1) * exp(x), (x + 1) ** 2 * exp(-x)])
assert nul(d + (x - 1) ** 2)
# ay-ex6 v2 : y' + y = 2x e^{−x}, f(0) = −3, f = (x²−3)e^{−x}, f' = −(x+1)(x−3)e^{−x}.
d = probleme(1, -3, (x**2 - 3) * exp(-x),
             [x**2 * exp(-x) - 3, (x**2 - 3) * exp(x), (x**2 - 3) * exp(-2 * x)])
assert nul(d + (x + 1) * (x - 3))
assert nul(((x**2 - 3) * exp(-x)).subs(x, -1) + 2 * E)
assert nul(((x**2 - 3) * exp(-x)).subs(x, 3) - 6 * exp(-3))

# --- ed-t-premier-ordre-ay-b-ex6 ---------------------------------------------
def retrouver(f0, fp0, lim, a_att, b_att, f_att, leurres):
    # f' = a f + b en 0 ; limite −b/a (a < 0)
    s = solve([a * f0 + b - fp0, -b / a - lim], [a, b], dict=True)
    assert s == [{a: a_att, b: b_att}]
    assert a_att < 0
    assert sol(f_att, lambda y: a_att * y + b_att, 0, f0)
    assert nul(diff(f_att, x).subs(x, 0) - fp0) and limit(f_att, x, oo) == lim
    for (aa, bb) in leurres:
        assert not (aa * f0 + bb == fp0 and aa < 0 and -bb / aa == lim)


R = Rational
# ay-b-ex6 : f(0) = 5, f'(0) = −2, limite 1 ⇒ a = −1/2, b = 1/2, f = 4e^{−x/2} + 1.
retrouver(5, -2, 1, R(-1, 2), R(1, 2), 4 * exp(-x / 2) + 1, [(R(1, 2), R(-1, 2)), (-2, 2), (R(-1, 2), R(-1, 2))])
# v1 : f(0) = 3, f'(0) = −4, limite −1 ⇒ a = b = −1, f = 4e^{−x} − 1.
retrouver(3, -4, -1, -1, -1, 4 * exp(-x) - 1, [(1, 1), (-1, 1), (-4, -4)])
# v2 : f(0) = 2, f'(0) = 3, limite 5 ⇒ a = −1, b = 5, f = 5 − 3e^{−x}.
retrouver(2, 3, 5, -1, 5, 5 - 3 * exp(-x), [(1, -5), (-1, -5), (R(-3, 5), 3)])

# --- ed-t-premier-ordre-primitives-ex6 ---------------------------------------
def prim(F, f, x0, y0):
    return nul(diff(F, x) - f) and nul(F.subs(x, x0) - y0)


# main : 1/(x ln x) sur ]1,+∞[, F(e) = 1 ⇒ F = ln(ln x) + 1.
f = 1 / (x * log(x))
assert prim(log(log(x)) + 1, f, E, 1)
for F in (log(x) ** 2 / 2 + R(1, 2), log(x), 1 / log(x)):
    assert nul(F.subs(x, E) - 1) and not nul(diff(F, x) - f)
# v1 : tan sur ]−π/2, π/2[, F(0) = 2 ⇒ F = 2 − ln(cos x).
assert prim(2 - log(cos(x)), tan(x), 0, 2)
for F in (2 + log(cos(x)), 1 + 1 / cos(x) ** 2, 2 + sin(x) ** 2 / 2):
    assert nul(F.subs(x, 0) - 2) and not nul(diff(F, x) - tan(x))
# v2 : e^x/(e^x+1) sur R, F(0) = 0 ⇒ F = ln(e^x + 1) − ln 2.
f = exp(x) / (exp(x) + 1)
assert prim(log(exp(x) + 1) - log(2), f, 0, 0)
for F in (log(exp(x) + 1), f - R(1, 2), x - log(exp(x) + 1) + log(2)):
    assert not prim(F, f, 0, 0)

print("enrichi_ed_t_premier_ordre : tout est vérifié.")
