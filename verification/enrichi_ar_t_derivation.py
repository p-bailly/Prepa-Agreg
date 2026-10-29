"""Vérifications calculatoires des questions de QCM (content/qcm/ar-t-derivation.json) et des exercices ajoutés au chapitre
ar-t-derivation (content/exercices/analyse-reelle/terminale-ar-t-derivation.json). Voir CLAUDE.md, « Contrôle de l'exactitude »."""
from sympy import (E, Rational as R, cos, diff, exp, factor, limit, oo, simplify, sin, solve, sqrt, symbols, pi, expand,
                   Interval, S, solveset, maximum)

x = symbols('x', real=True)

def egal(a, b):
    return simplify(a - b) == 0

def different(a, b, points=(R(-3, 2), R(-1, 3), R(1, 7), R(2, 9))):
    """a et b diffèrent en au moins un point de test (leurre faux)."""
    return any(simplify((a - b).subs(x, p)) != 0 for p in points)

# ar-t-derivation-q1 : (e^{x^2-3x})' = (2x-3) e^{x^2-3x} ; leurres faux.
f = exp(x**2 - 3 * x)
assert egal(diff(f, x), (2 * x - 3) * exp(x**2 - 3 * x))
for l in (exp(2 * x - 3), (x**2 - 3 * x) * exp(x**2 - 3 * x - 1), (2 * x - 3) * exp(2 * x - 3)):
    assert different(diff(f, x), l)

# ar-t-derivation-q2 : (sqrt(1-2x))' = -1/sqrt(1-2x) sur ]-oo, 1/2[ ; leurres faux.
f = sqrt(1 - 2 * x)
pts = (R(-3), R(-1, 2), R(0), R(1, 4))
assert egal(diff(f, x), -1 / sqrt(1 - 2 * x))
for l in (1 / (2 * sqrt(1 - 2 * x)), -2 / sqrt(1 - 2 * x), -1 / (2 * sqrt(1 - 2 * x))):
    assert different(diff(f, x), l, pts)

# ar-t-derivation-q3 : f = x^3 - 3x, f'' = 6x < 0 sur [-1, 0[ ; f(-1/2) = 11/8 > (f(-1)+f(0))/2 = 1 (affirmation fausse).
f = x**3 - 3 * x
assert diff(f, x, 2) == 6 * x
assert f.subs(x, -1) == 2 and f.subs(x, 0) == 0 and f.subs(x, R(-1, 2)) == R(11, 8) > 1

# ar-t-derivation-q4 : f'' = 12(x-1)^2 >= 0 : aucun point d'inflexion ; leurres 1, 2, 3 faux.
f = x**4 - 4 * x**3 + 6 * x**2 + x
assert egal(diff(f, x, 2), 12 * (x - 1)**2)

# ar-t-derivation-q5 : 1/(1+x) - (1-x) = x^2/(1+x) >= 0 pour x >= 0 ; h'' = 2/(1+x)^3 ; tangente en 0 : y = 1 - x.
h = 1 / (1 + x)
assert egal(h - (1 - x), x**2 / (1 + x)) and egal(diff(h, x, 2), 2 / (1 + x)**3)
assert h.subs(x, 0) == 1 and diff(h, x).subs(x, 0) == -1

# ar-t-derivation-q6 : tangente en 1 : 2 + 3(x-1), vaut 8 en 3. Contre-exemples : 3x - 1 et 2 + 3(x-1) + (x-1)^2.
T = 2 + 3 * (x - 1)
assert T.subs(x, 3) == 8
for g in (3 * x - 1, 2 + 3 * (x - 1) + (x - 1)**2):
    assert g.subs(x, 1) == 2 and diff(g, x).subs(x, 1) == 3 and diff(g, x, 2).is_nonnegative
g1, g2 = 3 * x - 1, 2 + 3 * (x - 1) + (x - 1)**2
assert g1.subs(x, 0) == -1 < 0 and g2.subs(x, 3) == 12 > 8 and g2.subs(x, 2) == 6 != 5

# ar-t-derivation-q7 : (x + cos x)' = 1 - sin x >= 0 (affirmation vraie).
assert egal(diff(x + cos(x), x), 1 - sin(x))

# ar-t-derivation-q8 : f' = 3(x+1)(x-3) <= 0 exactement sur [-1, 3] ; leurres faux.
f = x**3 - 3 * x**2 - 9 * x + 1
assert egal(diff(f, x), 3 * (x + 1) * (x - 3))
assert solveset(diff(f, x) <= 0, x, S.Reals) == Interval(-1, 3)
assert diff(f, x).subs(x, -2) > 0 and diff(f, x).subs(x, 4) > 0  # [-3,1], ]-oo,-1], [3,+oo[ contiennent des points de croissance

# ar-t-derivation-inflexion-ex5 : f(x) = (x+a)e^{-x} (a = 2, 3) et (x-1)e^x.
def etude(f, xinfl, extremum_x, extremum_val, tangente, leurres):
    f1, f2 = diff(f, x), diff(f, x, 2)
    assert f1.subs(x, extremum_x) == 0 and egal(f.subs(x, extremum_x), extremum_val)
    assert f2.subs(x, xinfl) == 0
    assert (f2.subs(x, xinfl - R(1, 2)) * f2.subs(x, xinfl + R(1, 2))) < 0  # changement de signe
    T = f.subs(x, xinfl) + f1.subs(x, xinfl) * (x - xinfl)
    assert egal(T, tangente)
    for l in leurres:
        assert different(T, l)
    d = f - T
    for p in (xinfl - 3, xinfl - 1, xinfl - R(1, 10)):
        assert d.subs(x, p).evalf() < 0  # courbe sous la tangente à gauche
    for p in (xinfl + R(1, 10), xinfl + 1, xinfl + 3):
        assert d.subs(x, p).evalf() > 0  # au-dessus à droite
    return f1, f2

f1, f2 = etude((x + 2) * exp(-x), 0, -1, E, 2 - x, [2 + x, -x, 2 - 2 * x])
assert egal(f1, -(x + 1) * exp(-x)) and egal(f2, x * exp(-x))
assert limit((x + 2) * exp(-x), x, oo) == 0 and limit((x + 2) * exp(-x), x, -oo) == -oo
# variante 1
f1, f2 = etude((x + 3) * exp(-x), -1, -2, E**2, -E * x + E, [-E * x + 2 * E, E * x + 3 * E, -x + 2 * E])
assert egal(f1, -(x + 2) * exp(-x)) and egal(f2, (x + 1) * exp(-x)) and egal(((x + 3) * exp(-x)).subs(x, -1), 2 * E)
assert limit((x + 3) * exp(-x), x, oo) == 0 and limit((x + 3) * exp(-x), x, -oo) == -oo
# variante 2 : minimum f(0) = -1, inflexion en -1, tangente y = -(x+3)/e.
f1, f2 = etude((x - 1) * exp(x), -1, 0, -1, -(x + 3) / E, [-(x + 1) / E, (x - 1) / E, -(x + 2) / E])
assert egal(f1, x * exp(x)) and egal(f2, (x + 1) * exp(x)) and egal(((x - 1) * exp(x)).subs(x, -1), -2 / E)
assert limit((x - 1) * exp(x), x, -oo) == 0 and limit((x - 1) * exp(x), x, oo) == oo
assert egal(f1.subs(x, -1), -1 / E)

# ar-t-derivation-convexite-ex5 : f'' = (x+2)e^x, (x-1)e^x, (x-2)e^{-x} ; convexe sur [-2,+oo[, [1,+oo[, [2,+oo[.
for f, fp, fs, c, leurres_bornes in (
        (x * exp(x), (x + 1) * exp(x), (x + 2) * exp(x), -2, ('[-1', ']-oo,-2', '[0')),
        ((x - 3) * exp(x), (x - 2) * exp(x), (x - 1) * exp(x), 1, ('[2', ']-oo,1', '[3')),
        (x * exp(-x), (1 - x) * exp(-x), (x - 2) * exp(-x), 2, (']-oo,1', ']-oo,2', '[1'))):
    assert egal(diff(f, x), fp) and egal(diff(f, x, 2), fs)
    assert solveset(fs >= 0, x, S.Reals) == Interval(c, oo)
# leurres : [-1,+oo[ et [0,+oo[ ne sont pas le plus grand intervalle (f'' > 0 sur ]-2,-1[) ; ]-oo,-2] : f'' < 0 ; etc.
assert ((x + 2) * exp(x)).subs(x, R(-3, 2)) > 0 and ((x + 2) * exp(x)).subs(x, -3) < 0
assert ((x - 1) * exp(x)).subs(x, R(3, 2)) > 0 and ((x - 1) * exp(x)).subs(x, 0) < 0
assert ((x - 2) * exp(-x)).subs(x, 0) < 0 and ((x - 2) * exp(-x)).subs(x, R(3, 2)) < 0

# ar-t-derivee-monotonie-ex6 : V(x) = x(c-2x)^2 sur ]0, c/2[, V' = 12(c/2 - x)(c/6 - x), maximum en c/6.
for c, Vmax, leurres in ((12, 128, (144, 108, 216)), (18, 432, (486, 384, 729)), (30, 2000, (2250, 1800, 3375))):
    V = x * (c - 2 * x)**2
    assert egal(diff(V, x), (c - 2 * x) * (c - 6 * x)) and egal(diff(V, x), 12 * (R(c, 2) - x) * (R(c, 6) - x))
    assert V.subs(x, R(c, 6)) == Vmax
    assert maximum(V, x, Interval.open(0, R(c, 2))) == Vmax
    assert all(l != Vmax for l in leurres)
assert (12 * x * (12 - 2 * x)**2).subs(x, 3) / 12 == 108  # V(3) = 108 pour la plaque de 12 cm
print('enrichi_ar_t_derivation : tout est vérifié')
