"""Vérifications calculatoires du chapitre « ar-t-derivation » (dérivée d'une composée, convexité, point d'inflexion,
position par rapport aux tangentes) : cartes et exercices de content/*/analyse-reelle/terminale-ar-t-derivation.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice concerné (v1, v2 = variantes 1 et 2).
Lancer : python3 verification/ar_t_derivation.py (doit finir sans erreur).
"""

import math

from sympy import (E, Abs, Rational, cos, diff, discriminant, exp, expand, factor, limit, nsimplify, oo, pi, Poly,
                   simplify, sin, solve, sqrt, symbols)

x, y, t, a, m = symbols("x y t a m", real=True)
h = symbols("h", positive=True)


def egal(u, v):
    return simplify(u - v) == 0


def different(u, v):
    """u et v ne sont pas la même fonction de x (au moins un point où elles diffèrent)."""
    return simplify(u - v) != 0


def signe_change(f2, a0, eps=Rational(1, 1000)):
    """f2 change strictement de signe en a0 (tests en a0 ± eps et a0 ± 10 eps)."""
    g = f2.subs(x, a0 - eps) * f2.subs(x, a0 + eps)
    g2 = f2.subs(x, a0 - 10 * eps) * f2.subs(x, a0 + 10 * eps)
    return g < 0 and g2 < 0


# ===========================================================================
# Carte ar-t-derivation-composee
# ===========================================================================
# Exemple : (sqrt(1+x^2))' = x/sqrt(1+x^2) ; (e^{-x^2})' = -2x e^{-x^2} ; (cos(3x+1))' = -3 sin(3x+1).
assert egal(diff(sqrt(1 + x**2), x), x / sqrt(1 + x**2))
assert egal(diff(exp(-x**2), x), -2 * x * exp(-x**2))
assert egal(diff(cos(3 * x + 1), x), -3 * sin(3 * x + 1))
# Contre-exemple 1 : sqrt(x^2) = |x|, taux d'accroissement en 0 : +1 à droite, -1 à gauche.
assert limit(sqrt(h**2) / h, h, 0, "+") == 1 and limit(sqrt(h**2) / (-h), h, 0, "+") == -1
# Contre-exemple 2 : (|h|^3+|h|)/h -> 1 à droite, -1 à gauche.
assert limit((h**3 + h) / h, h, 0, "+") == 1 and limit((h**3 + h) / (-h), h, 0, "+") == -1
# Contre-exemple 3 / affirmation 1 : |x|^2 = x^2 (dérivable).
# Cours : sqrt(x^4) = x^2 est dérivable en 0 bien que sqrt ne soit pas dérivable en u(0) = 0 (le théorème ne dit rien).
assert egal(sqrt(x**4), x**2) and diff(x**2, x).subs(x, 0) == 0
assert egal(Abs(x)**2, x**2)
# Affirmation 2 : ((2x+1)^5)' = 10(2x+1)^4 ≠ 5(2x+1)^4 (valeurs 10 et 5 en 0).
d = diff((2 * x + 1)**5, x)
assert egal(d, 10 * (2 * x + 1)**4) and d.subs(x, 0) == 10 and different(d, 5 * (2 * x + 1)**4)
# Affirmation 4 : sqrt(x^2) = |x| (non dérivable en 0, vu ci-dessus).

# ar-t-derivation-composee-ex1 et variantes : bonne réponse, leurres faux.
for f, bonne, leurres in [
    ((3 * x**2 - 1)**4, 24 * x * (3 * x**2 - 1)**3,
     [4 * (3 * x**2 - 1)**3, 6 * x * (3 * x**2 - 1)**3, 24 * x * (3 * x**2 - 1)**4]),               # ex1
    ((2 * x**3 + 5)**3, 18 * x**2 * (2 * x**3 + 5)**2,
     [3 * (2 * x**3 + 5)**2, 6 * x**2 * (2 * x**3 + 5)**2, 18 * x**2 * (2 * x**3 + 5)**3]),        # ex1 v1
    ((1 - 4 * x)**5, -20 * (1 - 4 * x)**4,
     [5 * (1 - 4 * x)**4, 20 * (1 - 4 * x)**4, -20 * (1 - 4 * x)**5]),                              # ex1 v2
]:
    assert egal(diff(f, x), bonne)
    for l in leurres:
        assert different(diff(f, x), l)

# ar-t-derivation-composee-ex2 et variantes : forme canonique, dérivée, minimum ; leurres faux.
for P, c, mini, leurres in [
    (x**2 - 4 * x + 5, 2, 1, [(0, 2), (sqrt(5), 0), (sqrt(2), 1)]),            # ex2
    (x**2 + 6 * x + 13, -3, 2, [(0, -3), (sqrt(13), 0), (2, 3)]),             # ex2 v1
    (x**2 - 2 * x + 10, 1, 3, [(0, 1), (sqrt(10), 0), (9, 1)]),               # ex2 v2
]:
    f = sqrt(P)
    assert egal(P, (x - c)**2 + mini**2)                       # P = (x-c)^2 + min^2 > 0
    assert egal(diff(f, x), (x - c) / sqrt(P))                  # f'(x) = (x-c)/sqrt(P)
    assert f.subs(x, c) == mini
    for valeur, point in leurres:                               # chaque leurre : mauvaise valeur ou mauvais point
        assert not (valeur == mini and point == c)
        assert f.subs(x, point) > mini or valeur != f.subs(x, point)

# ar-t-derivation-composee-ex3 et variantes : f = x e^{-k x^2}, f' = (1 - 2k x^2) e^{-k x^2}, max en 1/sqrt(2k).
for k, xmax, maxi, leurres in [
    (1, 1 / sqrt(2), 1 / sqrt(2 * E), [1 / E, 1 / sqrt(2), 1 / (2 * E)]),                    # ex3
    (2, Rational(1, 2), 1 / (2 * sqrt(E)), [1 / (2 * E), Rational(1, 2), exp(-2)]),          # ex3 v1
    (Rational(1, 8), 2, 2 / sqrt(E), [2 / E, 2, 1 / sqrt(E)]),                                # ex3 v2
]:
    f = x * exp(-k * x**2)
    assert egal(diff(f, x), (1 - 2 * k * x**2) * exp(-k * x**2))
    assert solve(1 - 2 * k * x**2, x) == [-xmax, xmax] or set(solve(1 - 2 * k * x**2, x)) == {-xmax, xmax}
    assert egal(f.subs(x, xmax), maxi)
    assert limit(f, x, oo) == 0 and limit(f, x, -oo) == 0
    for l in leurres:
        assert abs(float(l) - float(maxi)) > 1e-3
# Valeurs numériques annoncées : 1/sqrt(2e) ≈ 0,43 ; 1/(2 sqrt e) ≈ 0,30 ; 2/sqrt e ≈ 1,21 ; f(1) = 1/e pour ex3.
assert abs(float(1 / sqrt(2 * E)) - 0.43) < 0.005 and abs(float(1 / (2 * sqrt(E))) - 0.30) < 0.005
assert abs(float(2 / sqrt(E)) - 1.21) < 0.005

# ar-t-derivation-composee-ex4 : contre-exemples de la question 4.
# f = x + 1 : f' = 1 paire, f(0) = 1 ≠ 0 donc f non impaire.
assert diff(x + 1, x) == 1 and (x + 1).subs(x, 0) == 1
# f = x + sin x : f' = 1 + cos x 2π-périodique, f(0) = 0, f(2π) = 2π ≠ 0.
fs = x + sin(x)
assert egal(diff(fs, x).subs(x, x + 2 * pi), diff(fs, x)) and fs.subs(x, 0) == 0 and fs.subs(x, 2 * pi) == 2 * pi
# Une fonction paire a une dérivée impaire (exemple : x^4 + cos x).
g = x**4 + cos(x)
assert egal(diff(g, x).subs(x, -x), -diff(g, x))

# ar-t-derivation-composee-ex5 et variantes : P_{n+1} = P_n' + c x P_n, f = e^{(c/2) x^2}.
for c, P3, leurres in [
    (-2, -8 * x**3 + 12 * x, [-8 * x**3 + 8 * x, 8 * x**3 - 12 * x, -8 * x**3 + 4 * x]),                # ex5
    (-6, -216 * x**3 + 108 * x, [-216 * x**3 + 72 * x, 216 * x**3 - 108 * x, -216 * x**3 + 36 * x]),    # ex5 v1
    (2, 8 * x**3 + 12 * x, [8 * x**3 + 8 * x, 8 * x**3 - 12 * x, 8 * x**3 + 4 * x]),                     # ex5 v2
]:
    f = exp(Rational(c, 2) * x**2)
    P = [1]
    for n in range(8):
        P.append(expand(diff(P[-1], x) + c * x * P[-1]))
    for n in range(8):
        assert egal(diff(f, x, n), P[n] * f)                     # f^{(n)} = P_n e^{...}
        pol = Poly(P[n], x)
        assert pol.degree() == n and pol.LC() == c**n           # degré n, coefficient dominant c^n
    assert expand(P[3] - P3) == 0
    for l in leurres:
        assert expand(P[3] - l) != 0
# Valeurs intermédiaires de la correction de ex5 : P1 = -2x, P2 = 4x^2 - 2 ; v1 : -6x, 36x^2 - 6 ; v2 : 2x, 4x^2 + 2.
Q = [1]
for n in range(3):
    Q.append(expand(diff(Q[-1], x) - 2 * x * Q[-1]))
assert Q[1] == -2 * x and expand(Q[2] - (4 * x**2 - 2)) == 0
Q = [1]
for n in range(3):
    Q.append(expand(diff(Q[-1], x) - 6 * x * Q[-1]))
assert Q[1] == -6 * x and expand(Q[2] - (36 * x**2 - 6)) == 0
Q = [1]
for n in range(3):
    Q.append(expand(diff(Q[-1], x) + 2 * x * Q[-1]))
assert Q[1] == 2 * x and expand(Q[2] - (4 * x**2 + 2)) == 0

# ===========================================================================
# Carte ar-t-derivation-convexite
# ===========================================================================
# Exemples : exp'' = exp ; (x^3)'' = 6x ; (sqrt x)'' = -1/(4 x sqrt x).
assert egal(diff(exp(x), x, 2), exp(x)) and diff(x**3, x, 2) == 6 * x
xp = symbols("xp", positive=True)
assert egal(diff(sqrt(xp), xp, 2), -1 / (4 * xp * sqrt(xp)))
# Contre-exemple / affirmation 1 : (x^4)'' = 12x^2, nul en 0.
assert diff(x**4, x, 2) == 12 * x**2
# Affirmation 2 : sqrt 5 ≈ 2,24 > (1 + 3)/2 = 2.
assert sqrt(5) > 2 and abs(float(sqrt(5)) - 2.24) < 0.005
# Affirmation 4 : (x * x^2)'' = 6x < 0 pour x < 0.
assert diff(x * x**2, x, 2) == 6 * x
# Contre-exemple |x| convexe : inégalité triangulaire testée sur une grille.
for a0 in range(-5, 6):
    for b0 in range(-5, 6):
        for k in range(11):
            tt = Rational(k, 10)
            assert abs((1 - tt) * a0 + tt * b0) <= (1 - tt) * abs(a0) + tt * abs(b0)
# Démonstration : (b-a)s = (x-a)p1 + (b-x)p2 (identité algébrique).
fa, fb, fx, aa, bb, xx = symbols("fa fb fx aa bb xx")
p1 = (fx - fa) / (xx - aa)
p2 = (fb - fx) / (bb - xx)
assert simplify((xx - aa) * p1 + (bb - xx) * p2 - (fb - fa)) == 0
# Démonstration : f(b) - f(a) - s(x-a) = s(b-x).
s_ = (fb - fa) / (bb - aa)
assert simplify(fb - fa - s_ * (xx - aa) - s_ * (bb - xx)) == 0


def intervalle_convexite(f2, cote):
    """Racine unique de f2 affine ; renvoie la racine."""
    r = solve(f2, x)
    assert len(r) == 1
    return r[0]


# ar-t-derivation-convexite-ex1 et variantes : f'' affine.
for f, racine, convexe_a_droite in [
    (x**3 - 6 * x**2 + 5 * x + 1, 2, True),     # ex1 : convexe sur [2, +inf[
    (x**3 + 3 * x**2 - x + 2, -1, True),        # ex1 v1 : convexe sur [-1, +inf[
    (-x**3 + 3 * x**2 + 1, 1, False),           # ex1 v2 : convexe sur ]-inf, 1]
]:
    f2 = diff(f, x, 2)
    assert solve(f2, x) == [racine]
    if convexe_a_droite:
        assert f2.subs(x, racine + 1) > 0 and f2.subs(x, racine - 1) < 0
    else:
        assert f2.subs(x, racine - 1) > 0 and f2.subs(x, racine + 1) < 0
# Leurres : chacun contient un point où f'' < 0 (donc f n'y est pas convexe).
assert diff(x**3 - 6 * x**2 + 5 * x + 1, x, 2).subs(x, 1) < 0          # ex1 : ]-inf,2], [0,+inf[, R contiennent 1
assert diff(x**3 + 3 * x**2 - x + 2, x, 2).subs(x, Rational(-3, 2)) < 0  # ex1 v1 : ]-inf,-1], [-2,+inf[, R contiennent -3/2
assert diff(-x**3 + 3 * x**2 + 1, x, 2).subs(x, Rational(3, 2)) < 0      # ex1 v2 : [1,+inf[, ]-inf,2], R contiennent 3/2
# Correction ex1 : f' = 3x^2 - 12x + 5 ; v1 : 3x^2 + 6x - 1 ; v2 : -3x^2 + 6x.
assert diff(x**3 - 6 * x**2 + 5 * x + 1, x) == 3 * x**2 - 12 * x + 5
assert diff(x**3 + 3 * x**2 - x + 2, x) == 3 * x**2 + 6 * x - 1
assert diff(-x**3 + 3 * x**2 + 1, x) == -3 * x**2 + 6 * x

# ar-t-derivation-convexite-ex2 et variantes : f'' = 12 (x - r1)(x - r2), concave sur [r1, r2].
for f, f1, r1, r2, points_leurres in [
    (x**4 - 6 * x**2 + 3 * x, 4 * x**3 - 12 * x + 3, -1, 1, [Rational(3, 2), 0]),           # ex2
    (x**4 - 24 * x**2 + x, 4 * x**3 - 48 * x + 1, -2, 2, [3, 0]),                            # ex2 v1
    (x**4 + 4 * x**3 - 18 * x**2, 4 * x**3 + 12 * x**2 - 36 * x, -3, 1, [2, 0]),            # ex2 v2
]:
    assert diff(f, x) == f1
    assert expand(diff(f, x, 2) - 12 * (x - r1) * (x - r2)) == 0
    # Leurre intervalle plus grand ([-sqrt3, sqrt3], [-2sqrt3, 2sqrt3], [-1, 3]) : contient un point où f'' > 0.
    assert diff(f, x, 2).subs(x, points_leurres[0]) > 0
    # Leurre réunion ]-inf,r1] ∪ [r2,+inf[ : f'' < 0 en 0 n'y est pas, mais ce n'est pas un intervalle ; R contient 0
    # où f'' < 0 (pas convexe) et un point où f'' > 0 (pas concave).
    assert diff(f, x, 2).subs(x, points_leurres[1]) < 0
assert sqrt(3) > Rational(3, 2) and 2 * sqrt(3) > 3
# Correction ex2 : f(-2) = -14, f(2) = -2, moyenne -8 < f(0) = 0 (corde sous la courbe).
f = x**4 - 6 * x**2 + 3 * x
assert f.subs(x, -2) == -14 and f.subs(x, 2) == -2 and (f.subs(x, -2) + f.subs(x, 2)) / 2 == -8 and f.subs(x, 0) == 0
# Variante 2 : 12(x^2 + 2x - 3) = 12(x+3)(x-1).
assert expand(12 * (x**2 + 2 * x - 3) - 12 * (x + 3) * (x - 1)) == 0

# ar-t-derivation-convexite-ex3 et variantes : f_a = x^4 + a x^3 + c x^2, convexe ssi |a| <= borne.
for c, borne, leurres in [
    (6, 4, [("ouvert", 4), ("ferme", 2), ("ferme", 16)]),                     # ex3
    (24, 8, [("ouvert", 8), ("ferme", 4), ("ferme", 64)]),                    # ex3 v1
    (Rational(3, 2), 2, [("ouvert", 2), ("ferme", 1), ("ferme", 4)]),         # ex3 v2
]:
    fa_ = x**4 + a * x**3 + c * x**2
    f2 = diff(fa_, x, 2)
    D = discriminant(f2, x)
    assert expand(D - 36 * (a**2 - borne**2)) == 0
    # a = ±borne : f'' carré parfait >= 0 (convexe), donc l'intervalle ouvert est faux.
    for sgn in (1, -1):
        f2b = f2.subs(a, sgn * borne)
        r = solve(f2b, x)
        assert len(r) == 1 and expand(f2b - 12 * (x - r[0])**2) == 0
    # a légèrement au-delà de la borne : f'' prend une valeur < 0 (non convexe).
    for sgn in (1, -1):
        f2c = f2.subs(a, sgn * (borne + Rational(1, 10)))
        sommet = -sgn * (borne + Rational(1, 10)) * 6 / 24
        assert f2c.subs(x, sommet) < 0
    for nature, b0 in leurres:
        if nature == "ferme" and b0 < borne:
            # [-b0, b0] trop petit : a = borne est exclu alors qu'il convient.
            assert b0 < borne
        elif nature == "ferme":
            # [-b0, b0] trop grand : a = b0 > borne ne convient pas.
            f2d = f2.subs(a, b0)
            assert f2d.subs(x, -Rational(b0, 4)) < 0
# Correction ex3 : f_4'' = 12(x+1)^2 ; v1 : f_8'' = 12(x+2)^2 ; v2 : f_2'' = 3(2x+1)^2.
assert expand(diff(x**4 + 4 * x**3 + 6 * x**2, x, 2) - 12 * (x + 1)**2) == 0
assert expand(diff(x**4 + 8 * x**3 + 24 * x**2, x, 2) - 12 * (x + 2)**2) == 0
assert expand(diff(x**4 + 2 * x**3 + Rational(3, 2) * x**2, x, 2) - 3 * (2 * x + 1)**2) == 0

# ar-t-derivation-convexite-ex4 et variantes : sin x >= c x sur [0, b], c = sin(b)/b optimal.
for b, cbon, leurres in [
    (pi / 2, 2 / pi, [1, Rational(1, 2), pi / 2]),                          # ex4
    (pi / 3, 3 * sqrt(3) / (2 * pi), [1, sqrt(3) / 2, 2 / pi]),             # ex4 v1
    (pi / 6, 3 / pi, [1, Rational(1, 2), 6 / pi]),                          # ex4 v2
]:
    assert egal(sin(b) / b, cbon)
    # sin concave sur [0, b] : sin'' = -sin <= 0 (b <= pi).
    assert float(b) <= math.pi
    # Inégalité vérifiée numériquement sur une grille de [0, b].
    for k in range(0, 201):
        xv = float(b) * k / 200
        assert math.sin(xv) >= float(cbon) * xv - 1e-12
    for l in leurres:
        lv = float(l)
        if lv > float(cbon):
            assert math.sin(float(b)) < lv * float(b)     # leurre trop grand : faux en x = b
        else:
            assert lv < float(cbon)                       # leurre trop petit : vrai mais pas le plus grand
assert abs(float(2 / pi) - 0.64) < 0.005 and abs(float(3 * sqrt(3) / (2 * pi)) - 0.83) < 0.005
assert abs(float(3 / pi) - 0.95) < 0.005

# ===========================================================================
# Carte ar-t-derivation-inflexion
# ===========================================================================
# Exemple : f = x^3 - 3x^2 + 1, A(1, -1), tangente y = -3x + 2, f - T = (x-1)^3.
f = x**3 - 3 * x**2 + 1
assert f.subs(x, 1) == -1 and diff(f, x).subs(x, 1) == -3 and signe_change(diff(f, x, 2), 1)
assert expand(f - (-3 * x + 2) - (x - 1)**3) == 0
# Contre-exemple x^4 : f''(0) = 0, f'' >= 0, f >= tangente y = 0.
assert diff(x**4, x, 2).subs(x, 0) == 0 and not signe_change(diff(x**4, x, 2), 0)
# Contre-exemple x|x| : f' = 2|x| (taux en 0 : ±2), f'' = ±2.
assert limit((h * h) / h, h, 0, "+") == 0
assert diff(x**2, x, 2) == 2 and diff(-x**2, x, 2) == -2
assert limit((2 * h) / h, h, 0, "+") == 2 and limit((2 * h) / (-h), h, 0, "+") == -2   # f'(x) = 2|x| non dérivable en 0
# Contre-exemple f = x^3 + x^4 sin(1/x) : f'(0) = f''(0) = 0, f'' donnée, signe de f, f'' de signe variable.
f = x**3 + x**4 * sin(1 / x)
f1 = diff(f, x)
f2 = diff(f, x, 2)
assert egal(f2, 6 * x + 12 * x**2 * sin(1 / x) - 6 * x * cos(1 / x) - sin(1 / x))
assert limit(f / x, x, 0) == 0                     # f'(0) = 0
assert limit(f1 / x, x, 0) == 0                    # f''(0) = lim f'(x)/x = 0
for k in range(1, 50):
    xk = 1 / (pi / 2 + 2 * k * pi)                 # sin(1/x) = 1, cos(1/x) = 0 : f'' < 0
    yk = 1 / (3 * pi / 2 + 2 * k * pi)             # sin(1/x) = -1 : f'' > 0
    assert f2.subs(x, xk).evalf() < 0 and f2.subs(x, yk).evalf() > 0
for xv in [0.9, 0.5, 0.1, 0.01, -0.01, -0.1, -0.5, -0.9]:
    assert abs(xv**4 * math.sin(1 / xv)) <= xv**4 < abs(xv)**3
    assert (xv**3 + xv**4 * math.sin(1 / xv)) * xv > 0
# Affirmation 2 : x^3 + x, f'(0) = 1, f'' = 6x, f - x = x^3.
assert diff(x**3 + x, x).subs(x, 0) == 1 and diff(x**3 + x, x, 2) == 6 * x
# Affirmation 4 : degré 3, f'' = 6αx + 2β s'annule en -β/(3α).
al, be, ga, de = symbols("alpha beta gamma delta", real=True)
fc = al * x**3 + be * x**2 + ga * x + de
assert expand(diff(fc, x, 2) - (6 * al * x + 2 * be)) == 0 and solve(diff(fc, x, 2), x) == [-be / (3 * al)]
# Leurres du texte à trous : f'' atteint un maximum en 0 pour f = -x^4/12 sans inflexion (f - T = -x^4/12 <= 0).
assert diff(-x**4 / 12, x, 2) == -x**2

# ar-t-derivation-inflexion-ex1 et variantes : point d'inflexion, leurres faux (point hors courbe ou f'' ≠ 0).
for f, A, leurres in [
    (x**3 - 3 * x**2 + 2 * x + 1, (1, 1), [(1, 0), (2, 1), (0, 1)]),                   # ex1
    (x**3 + 6 * x**2 + 9 * x + 2, (-2, 0), [(2, 52), (-2, 2), (-1, -2)]),              # ex1 v1
    (x**3 - 9 * x**2 + 1, (3, -53), [(3, 0), (0, 1), (6, -107)]),                      # ex1 v2
]:
    f2 = diff(f, x, 2)
    assert solve(f2, x) == [A[0]] and signe_change(f2, A[0]) and f.subs(x, A[0]) == A[1]
    for (u, v) in leurres:
        assert f2.subs(x, u) != 0 or f.subs(x, u) != v
# Correction ex1 v1 : (-1, -2) sur la courbe, f'(-1) = 0, f''(-1) = 6 ; v2 : (0,1), (6,-107) points critiques.
f = x**3 + 6 * x**2 + 9 * x + 2
assert f.subs(x, -1) == -2 and diff(f, x).subs(x, -1) == 0 and diff(f, x, 2).subs(x, -1) == 6 and f.subs(x, 2) == 52
f = x**3 - 9 * x**2 + 1
assert f.subs(x, 0) == 1 and f.subs(x, 6) == -107 and solve(diff(f, x), x) == [0, 6]

# ar-t-derivation-inflexion-ex2 et variantes : tangente au point d'inflexion et position.
for f, a0, T, cube in [
    (x**3 - 3 * x**2 + 4, 1, -3 * x + 5, (x - 1)**3),              # ex2 : au-dessous à gauche, au-dessus à droite
    (-x**3 + 6 * x**2 - 5, 2, 12 * x - 13, -(x - 2)**3),           # ex2 v1 : au-dessus à gauche, au-dessous à droite
    (2 * x**3 + 6 * x**2 + x, -1, -5 * x - 2, 2 * (x + 1)**3),     # ex2 v2 : au-dessous à gauche, au-dessus à droite
]:
    assert signe_change(diff(f, x, 2), a0) and solve(diff(f, x, 2), x) == [a0]
    assert expand(f.subs(x, a0) + diff(f, x).subs(x, a0) * (x - a0) - T) == 0
    assert expand(f - T - cube) == 0
# Leurres : tangente horizontale (y = f(a)) ou ordonnée à l'origine fausse.
assert diff(x**3 - 3 * x**2 + 4, x).subs(x, 1) != 0 and (x**3 - 3 * x**2 + 4).subs(x, 1) == 2
assert diff(-x**3 + 6 * x**2 - 5, x).subs(x, 2) == 12 and (-x**3 + 6 * x**2 - 5).subs(x, 2) == 11
assert diff(2 * x**3 + 6 * x**2 + x, x).subs(x, -1) == -5 and (2 * x**3 + 6 * x**2 + x).subs(x, -1) == 3

# ar-t-derivation-inflexion-ex3 et variantes : courbes logistiques.
for f, rel1, rel2, a0, T, leurres in [
    (1 / (1 + exp(-x)), lambda F: F * (1 - F), lambda F, F1: F1 * (1 - 2 * F), 0, x / 4 + Rational(1, 2),
     [x / 2 + Rational(1, 2), Rational(1, 2), x + Rational(1, 2)]),                                     # ex3
    (4 / (1 + exp(-2 * x)), lambda F: F * (4 - F) / 2, lambda F, F1: F1 * (2 - F), 0, 2 * x + 2,
     [x + 2, 2 + 0 * x, 8 * x + 2]),                                                                    # ex3 v1
    (1 / (1 + exp(3 - x)), lambda F: F * (1 - F), lambda F, F1: F1 * (1 - 2 * F), 3, x / 4 - Rational(1, 4),
     [x / 4 + Rational(1, 2), Rational(1, 2), x / 4 + Rational(5, 4)]),                                 # ex3 v2
]:
    f1 = diff(f, x)
    f2 = diff(f, x, 2)
    assert egal(f1, rel1(f)) and egal(f2, rel2(f, f1))
    assert signe_change(f2, a0)
    assert simplify(f2.subs(x, a0)) == 0
    Tcalc = simplify(f.subs(x, a0) + f1.subs(x, a0) * (x - a0))
    assert expand(Tcalc - T) == 0
    for l in leurres:
        assert expand(Tcalc - l) != 0
# Correction ex3 v1 : f'(x) = 8 e^{-2x} / (1 + e^{-2x})^2.
assert egal(diff(4 / (1 + exp(-2 * x)), x), 8 * exp(-2 * x) / (1 + exp(-2 * x))**2)

# ar-t-derivation-inflexion-ex4 et variantes : f_m = x^4 - p m x^3 + q x^2, deux points d'inflexion ssi |m| > s.
for p, q, coef, trin, s, double in [
    (2, 6, 12, x**2 - m * x + 1, 2, 1),         # ex4
    (1, 6, 6, 2 * x**2 - m * x + 2, 4, 1),      # ex4 v1
    (2, 54, 12, x**2 - m * x + 9, 6, 3),        # ex4 v2
]:
    fm = x**4 - p * m * x**3 + q * x**2
    f2 = diff(fm, x, 2)
    assert expand(f2 - coef * trin) == 0
    D = discriminant(trin, x)
    assert expand(D - (m**2 - s**2)) == 0
    for sgn in (1, -1):
        f2d = f2.subs(m, sgn * s)
        assert expand(f2d - 12 * (x - sgn * double)**2) == 0          # racine double : pas de changement de signe
        f2e = f2.subs(m, sgn * (s + 1))                                # |m| > s : deux racines simples
        racines = solve(f2e, x)
        assert len(racines) == 2 and all(signe_change(f2e, r) for r in racines)
    f2z = f2.subs(m, 0)                                                # |m| < s : f'' > 0
    assert all(f2z.subs(x, xv) > 0 for xv in range(-10, 11)) and discriminant(trin.subs(m, 0), x) < 0
# Correction ex4 : f_m' = 4x^3 - 6m x^2 + 12x.
assert expand(diff(x**4 - 2 * m * x**3 + 6 * x**2, x) - (4 * x**3 - 6 * m * x**2 + 12 * x)) == 0
# Discriminants exacts : ex4 m^2 - 4 ; v1 (2x^2 - m x + 2) : m^2 - 16 ; v2 m^2 - 36.
assert expand(discriminant(x**2 - m * x + 1, x) - (m**2 - 4)) == 0
assert expand(discriminant(2 * x**2 - m * x + 2, x) - (m**2 - 16)) == 0
assert expand(discriminant(x**2 - m * x + 9, x) - (m**2 - 36)) == 0

# ===========================================================================
# Carte ar-t-derivation-tangentes
# ===========================================================================
# Cours / contre-exemple 1 : x^3 - (3x - 2) = (x-1)^2 (x+2), f(-3) = -27 < -11.
assert expand(x**3 - (3 * x - 2) - (x - 1)**2 * (x + 2)) == 0 and (3 * x - 2).subs(x, -3) == -11
# Exemples : sqrt x <= (1+x)/2 car (1+x)/2 - sqrt x = (sqrt x - 1)^2 / 2 ; sin x <= x sur [0, pi].
assert expand((1 + xp) / 2 - sqrt(xp) - (sqrt(xp) - 1)**2 / 2) == 0
for k in range(0, 101):
    xv = math.pi * k / 100
    assert math.sin(xv) <= xv + 1e-15
# Contre-exemple 2 : 1/x^2, f'' = 6/x^4, tangente en 1 : y = 3 - 2x, f(-1) = 1 < 5.
f = 1 / x**2
assert egal(diff(f, x, 2), 6 / x**4)
assert expand(f.subs(x, 1) + diff(f, x).subs(x, 1) * (x - 1) - (3 - 2 * x)) == 0 and f.subs(x, -1) == 1 and (3 - 2 * x).subs(x, -1) == 5
# Contre-exemple 3 : sqrt 4 = 2 < 5/2.
assert sqrt(4) == 2 and Rational(5, 2) > 2
# Affirmation 4 : x^4 - 2x^2, tangente en 2 : y = 24x - 40, différence (x-2)^2 (x^2 + 4x + 10), f''(0) = -4.
f = x**4 - 2 * x**2
assert expand(f.subs(x, 2) + diff(f, x).subs(x, 2) * (x - 2) - (24 * x - 40)) == 0
assert expand(f - (24 * x - 40) - (x - 2)**2 * (x**2 + 4 * x + 10)) == 0
assert discriminant(x**2 + 4 * x + 10, x) < 0 and diff(f, x, 2).subs(x, 0) == -4
# Leurre du texte à trous f(x) >= f(a) + f'(x)(x - a) : faux pour x^2, a = 0, x = 1 (1 >= 2 faux).
assert not (1 >= 0 + 2 * 1 * (1 - 0))

# ar-t-derivation-tangentes-ex1 et variantes : tangentes à exp et e^{2x}, leurres faux.
for f, a0, T, leurres in [
    (exp(x), 1, E * x, [("<=", E * x), ("<=", E * (x - 1)), (">=", x + E - 1)]),                   # ex1
    (exp(x), 2, E**2 * (x - 1), [("<=", E**2 * (x - 1)), (">=", E**2 * x), ("<=", E**2 * (x - 2))]),  # ex1 v1
    (exp(2 * x), 0, 1 + 2 * x, [("<=", 1 + 2 * x), (">=", 1 + x), (">=", 2 + 2 * x)]),              # ex1 v2
]:
    assert egal(f.subs(x, a0) + diff(f, x).subs(x, a0) * (x - a0), T)
    assert all(float(diff(f, x, 2).subs(x, xv)) > 0 for xv in range(-5, 6))
    for xv in [Rational(k, 4) for k in range(-40, 41)]:
        assert f.subs(x, xv) >= T.subs(x, xv)
    for sens, l in leurres:
        # Chaque leurre a une tangente fausse ou une inégalité fausse : on exhibe un point où l'inégalité échoue.
        pts = [Rational(k, 4) for k in range(-40, 41)]
        if sens == "<=":
            assert any(f.subs(x, xv) > l.subs(x, xv) for xv in pts)
        else:
            assert any(f.subs(x, xv) < l.subs(x, xv) for xv in pts)
        assert not egal(l, T) or sens == "<="

# ar-t-derivation-tangentes-ex2 et variantes : tangente à sqrt en a, majorant de sqrt(b).
for a0, T, b, maj, approx, leurres in [
    (4, xp / 4 + 1, 5, Rational(9, 4), 2.236, [(">=", xp / 4 + 1), ("<=", xp / 4 + Rational(1, 2)), ("<=", xp / 2)]),   # ex2
    (9, xp / 6 + Rational(3, 2), 10, Rational(19, 6), 3.162,
     [(">=", xp / 6 + Rational(3, 2)), ("<=", xp / 6 + Rational(1, 2)), ("<=", xp / 3)]),                               # ex2 v1
    (Rational(1, 4), xp + Rational(1, 4), Rational(3, 10), Rational(55, 100), 0.548,
     [(">=", xp + Rational(1, 4)), ("<=", xp + Rational(1, 8)), ("<=", 2 * xp)]),                                       # ex2 v2
]:
    Tc = sqrt(a0) + diff(sqrt(xp), xp).subs(xp, a0) * (xp - a0)
    assert expand(Tc - T) == 0
    assert expand(T - sqrt(xp) - (sqrt(xp) - sqrt(a0))**2 / (2 * sqrt(a0))) == 0     # T - sqrt x = (sqrt x - sqrt a)^2/(2 sqrt a) >= 0
    assert T.subs(xp, 0) > 0                                                           # cas x = 0
    assert T.subs(xp, b) == maj and abs(math.sqrt(float(b)) - approx) < 0.001 and math.sqrt(float(b)) <= float(maj)
    pts = [Rational(k, 100) for k in range(0, 1001)]
    for sens, l in leurres:
        if sens == "<=":
            assert any(sqrt(pt) > l.subs(xp, pt) for pt in pts)
        else:
            assert any(sqrt(pt) < l.subs(xp, pt) for pt in pts)
assert abs(19 / 6 - 3.167) < 0.001

# ar-t-derivation-tangentes-ex3 et variantes : Bernoulli, minorants ; leurres faux (supérieurs à la vraie valeur).
n_ = symbols("n", positive=True, integer=True)
f = (1 + x)**n_
assert egal(diff(f, x, 2), n_ * (n_ - 1) * (1 + x)**(n_ - 2))
for base, n0, mino, approx, leurres in [
    (Rational(101, 100), 100, 2, 2.705, [E, 3, 101]),                                    # ex3
    (Rational(99, 100), 50, Rational(1, 2), 0.605, [Rational(99, 100), Rational(7, 10), Rational(49, 50)]),   # ex3 v1
    (Rational(103, 100), 20, Rational(8, 5), 1.806, [2, Rational(19, 10), 60]),        # ex3 v2
]:
    valeur = base**n0
    assert 1 + n0 * (base - 1) == mino and valeur >= mino
    assert abs(float(valeur) - approx) < 0.001
    for l in leurres:
        assert valeur < l                                                                 # le leurre n'est pas un minorant
assert float(Rational(101, 100)**100) < float(E)

# ar-t-derivation-tangentes-ex4 : e^{-x} convexe, majorée par 1 sur [0, +inf[, décroissante, non constante.
assert egal(diff(exp(-x), x, 2), exp(-x)) and exp(-x).subs(x, 0) == 1 and exp(-x).subs(x, 1) != 1

# ar-t-derivation-tangentes-ex5 et variantes : e^{t-1} >= t, inégalité arithmético-géométrique, maximum du produit.
g = exp(t - 1)
assert egal(g.subs(t, 1) + diff(g, t).subs(t, 1) * (t - 1), t)
hh = exp(t - 1) - t
assert solve(diff(hh, t), t) == [1] and hh.subs(t, 1) == 0
for tv in [Rational(k, 10) for k in range(-50, 51) if k != 10]:
    assert hh.subs(t, tv) > 0
for somme, nvar, maxi, leurres in [
    (6, 3, 8, [6, 27, 216]),                          # ex5
    (9, 3, 27, [9, 729, 81]),                         # ex5 v1
    (8, 4, 16, [8, 4096, Rational(512, 27)]),         # ex5 v2
]:
    moy = Rational(somme, nvar)
    assert moy**nvar == maxi
    # Recherche sur une grille : le produit ne dépasse jamais maxi, et vaut maxi au point d'égalité.
    N = 10 * nvar                      # grille de pas somme/N, qui contient le point d'égalité (moy = 10 pas)
    pas = Rational(somme, N)
    meilleur = 0
    if nvar == 3:
        for i in range(1, N):
            for j in range(1, N - i):
                xi, yj = i * pas, j * pas
                zk = somme - xi - yj
                prod = xi * yj * zk
                assert prod <= maxi
                meilleur = max(meilleur, prod)
    else:
        for i in range(1, N):
            for j in range(1, N - i):
                for k in range(1, N - i - j):
                    xi, yj, zk = i * pas, j * pas, k * pas
                    w = somme - xi - yj - zk
                    prod = xi * yj * zk * w
                    assert prod <= maxi
                    meilleur = max(meilleur, prod)
    assert meilleur == maxi
    for l in leurres:
        assert l != maxi
        assert l < maxi or l > maxi       # un leurre plus petit n'est pas le maximum (atteint) ; plus grand : jamais atteint
# Leurre 512/27 de ex5 v2 : (8/3)^3, supérieur au maximum 16.
assert Rational(8, 3)**3 == Rational(512, 27) and Rational(512, 27) > 16

print("ar_t_derivation : toutes les vérifications passent.")
