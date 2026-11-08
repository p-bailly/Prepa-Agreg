"""Vérifications calculatoires du chapitre « geo-l3-coniques » (géométrie, L3) : cartes, exercices et QCM de
content/*/geometrie/l3-geo-l3-coniques.json et content/qcm/geo-l3-coniques.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question vérifiée.
Calcul exact avec SymPy (rationnels, radicaux), sauf quelques contrôles numériques sur des points tirés au hasard.
Lancer : python3 verification/geo_l3_coniques.py.
"""

import random
from fractions import Fraction as Fr

import sympy as sp

x, y, z, t, u, v, th = sp.symbols('x y z t u v theta', real=True)
a, b, c, e, p, h, k, m = sp.symbols('a b c e p h k m', positive=True)
random.seed(1)


def egal(A, B):
    return sp.simplify(sp.nsimplify(A) - sp.nsimplify(B)) == 0


def dist(P, Q):
    return sp.sqrt(sum((pi - qi) ** 2 for pi, qi in zip(P, Q)))


# ---------------------------------------------------------------------------------------------------------------
# Carte geo-l3-coniques-definitions
# ---------------------------------------------------------------------------------------------------------------
# Identité (2) de la démonstration : MF^2 - (a - e x)^2 = (1 - e^2) x^2 + y^2 - (a^2 - c^2) avec c = e a.
MF2 = (x - e * a) ** 2 + y ** 2
assert sp.expand(MF2 - (a - e * x) ** 2 - ((1 - e ** 2) * x ** 2 + y ** 2 - (a ** 2 - (e * a) ** 2))) == 0
# (1) : MF'^2 - MF^2 = 4 c x
assert sp.expand(((x + c) ** 2 + y ** 2) - ((x - c) ** 2 + y ** 2) - 4 * c * x) == 0
# a^2/c = a/e quand e = c/a
assert sp.simplify(a ** 2 / c - a / (c / a)) == 0


def bifocal_ellipse_ok(cc, aa, n=200):
    """Contrôle numérique : sur l'ellipse x^2/a^2 + y^2/b^2 = 1, MF + MF' = 2a et MF = e d(M, D)."""
    bb = (aa ** 2 - cc ** 2) ** 0.5
    ee = cc / aa
    for _ in range(n):
        tt = random.uniform(0, 6.3)
        X, Y = aa * sp.cos(tt).evalf(), bb * sp.sin(tt).evalf()
        X, Y = float(X), float(Y)
        mf = ((X - cc) ** 2 + Y ** 2) ** 0.5
        mf2 = ((X + cc) ** 2 + Y ** 2) ** 0.5
        assert abs(mf + mf2 - 2 * aa) < 1e-9
        assert abs(mf - ee * abs(X - aa / ee)) < 1e-9
    return True


def bifocal_hyperbole_ok(cc, aa, n=200):
    bb = (cc ** 2 - aa ** 2) ** 0.5
    ee = cc / aa
    import math
    for _ in range(n):
        tt = random.uniform(-3, 3)
        for s in (1, -1):
            X, Y = s * aa * math.cosh(tt), bb * math.sinh(tt)
            mf = ((X - cc) ** 2 + Y ** 2) ** 0.5
            mf2 = ((X + cc) ** 2 + Y ** 2) ** 0.5
            assert abs(abs(mf - mf2) - 2 * aa) < 1e-7 * (1 + abs(X))
            assert abs(mf - ee * abs(X - aa / ee)) < 1e-7 * (1 + abs(X))
    return True


assert bifocal_ellipse_ok(3.0, 5.0) and bifocal_ellipse_ok(0.5, 2.0)
assert bifocal_hyperbole_ok(5.0, 4.0) and bifocal_hyperbole_ok(5.0, 3.0)
# exemple : foyers (±3,0), 2a = 10 : b^2 = 16, e = 3/5, directrices ±25/3 ; foyers (±5,0), a = 4 : b^2 = 9, e = 5/4, ±16/5
assert 25 - 9 == 16 and Fr(25, 3) == Fr(5 ** 2, 3) and 25 - 16 == 9 and Fr(16, 5) == Fr(4 ** 2, 5)
# réciproque : l'autre point à distance c de F convient seulement si e = 1
ee = sp.symbols('ee', positive=True)
aa = sp.Integer(1)
hh_ell = aa / ee - ee * aa  # e < 1 : h = a/e - c
assert set(sp.solve(sp.Eq(sp.Abs(hh_ell - ee * aa) ** 2, (aa / ee) ** 2), ee)) <= {1}
hh_hyp = ee * aa - aa / ee  # e > 1 : h = c - a/e
assert sp.solve(sp.Eq(hh_hyp + ee * aa, aa / ee), ee) == [1]
# a = e h / |1 - e^2| : pour e < 1, h = a/e - c donne bien a
assert sp.simplify(ee * hh_ell / (1 - ee ** 2) - aa) == 0
assert sp.simplify(ee * hh_hyp / (ee ** 2 - 1) - aa) == 0
# contre-exemples : (3,0) sur x^2/9 - y^2/16 = 1, MF - MF' = 2 - 8 = -6 avec F = (5,0)
assert dist((3, 0), (5, 0)) - dist((3, 0), (-5, 0)) == -6
# directrice du mauvais côté : sommet (5,0), MF = 2, (3/5)(5 + 25/3) = 8
assert dist((5, 0), (3, 0)) == 2 and Fr(3, 5) * (5 + Fr(25, 3)) == 8
# trou 1, leurre MF.MF' = a^2 : le point de l'axe vérifiant cette égalité est x = sqrt(a^2 + c^2) != a
sols_cassini = sp.solve(sp.Eq((x - c) * (x + c), a ** 2), x)
assert set(sols_cassini) == {sp.sqrt(a ** 2 + c ** 2), -sp.sqrt(a ** 2 + c ** 2)}
assert all(sp.simplify(s_ ** 2 - a ** 2) != 0 for s_ in sols_cassini)
# trou 1, leurre MF + MF' = a : le sommet (a,0) donne 2a != a ; leurre |MF - MF'| = 2a : vide car |MF - MF'| <= 2c < 2a
# trou 2 : avec e' = a/c, c^2/a^2 ou 2c/a, le sommet (a,0) n'est pas sur la conique de foyer F et de directrice D
for (A_, C_) in [(5, 3), (13, 12), (2, 1)]:
    A_, C_ = sp.Integer(A_), sp.Integer(C_)
    d_vertex = A_ ** 2 / C_ - A_  # distance du sommet (a,0) à la directrice x = a^2/c
    assert A_ - C_ == (C_ / A_) * d_vertex  # la bonne excentricité
    for eps in [A_ / C_, C_ ** 2 / A_ ** 2, 2 * C_ / A_]:
        assert A_ - C_ != eps * d_vertex  # (trou 2, leurres)
    # trou 3 : OH = c^2/a, OH = a, FH = a^2/c (OH = a^2/c + c) diffèrent de a^2/c
    for OH in [C_ ** 2 / A_, A_, A_ ** 2 / C_ + C_]:
        assert OH != A_ ** 2 / C_
    # pièges : directrice à c^2/a, ou du côté de F' : le sommet (a,0) ne vérifie plus MF = e d
    assert A_ - C_ != (C_ / A_) * abs(A_ - C_ ** 2 / A_)
    assert A_ - C_ != (C_ / A_) * (A_ + A_ ** 2 / C_)

# ---------------------------------------------------------------------------------------------------------------
# Carte geo-l3-coniques-reduites
# ---------------------------------------------------------------------------------------------------------------
# Repère d'origine F, D : x = h, p = e h : x^2 + y^2 - e^2 (x - h)^2 = (1 - e^2) x^2 + 2 e p x + y^2 - p^2
expr = x ** 2 + y ** 2 - e ** 2 * (x - h) ** 2
assert sp.expand(expr - ((1 - e ** 2) * x ** 2 + 2 * e * (e * h) * x + y ** 2 - (e * h) ** 2)) == 0
# complétion du carré : x0 = e p/(e^2 - 1), (1-e^2) X^2 + y^2 = p^2/(1-e^2)
x0 = e * p / (e ** 2 - 1)
lhs = (1 - e ** 2) * x ** 2 + 2 * e * p * x + y ** 2 - p ** 2
assert sp.simplify(lhs - ((1 - e ** 2) * (x - x0) ** 2 + y ** 2 - p ** 2 / (1 - e ** 2))) == 0
# e < 1 : a = p/(1-e^2), b^2 = p^2/(1-e^2), foyer X = e a, directrice X = h - x0 = a/e ; c^2 = a^2 - b^2 ; b^2/a = p
A_ = p / (1 - e ** 2)
B2 = p ** 2 / (1 - e ** 2)
assert sp.simplify(-x0 - e * A_) == 0 and sp.simplify(p / e - x0 - A_ / e) == 0
assert sp.simplify((e * A_) ** 2 - (A_ ** 2 - B2)) == 0 and sp.simplify(B2 / A_ - p) == 0
# e > 1 : a = p/(e^2-1), b^2 = p^2/(e^2-1), foyer X = -e a, directrice X = -a/e ; c^2 = a^2 + b^2
A_ = p / (e ** 2 - 1)
B2 = p ** 2 / (e ** 2 - 1)
assert sp.simplify(-x0 + e * A_) == 0 and sp.simplify(p / e - x0 + A_ / e) == 0
assert sp.simplify((e * A_) ** 2 - (A_ ** 2 + B2)) == 0 and sp.simplify(B2 / A_ - p) == 0
# réciproque : MF^2 - e^2 d^2 = b^2 (x^2/a^2 + y^2/b^2 - 1) (ellipse), -b^2(...) (hyperbole), y^2 - 2px (parabole)
for signe in (1, -1):
    bb2 = signe * (a ** 2 - c ** 2)  # b^2 = a^2 - c^2 ou c^2 - a^2
    ident = (x - c) ** 2 + y ** 2 - (c / a) ** 2 * (x - a ** 2 / c) ** 2
    cible = signe * bb2 * (x ** 2 / a ** 2 + signe * y ** 2 / bb2 - 1)
    assert sp.simplify(ident - cible) == 0
assert sp.expand((x - p / 2) ** 2 + y ** 2 - (x + p / 2) ** 2 - (y ** 2 - 2 * p * x)) == 0
# exemple : x^2/9 + y^2/5 : c = 2, e = 2/3, p = 5/3, h = 9/2 - 2 = 5/2, e h = p ; y^2 = 6x : p = 3 ; x^2/4 - y^2/12 : c = 4, e = 2, p = 6
assert Fr(2, 3) * (Fr(9, 2) - 2) == Fr(5, 3) == Fr(5, 3)
assert sp.sqrt(4 + 12) == 4 and Fr(12, 2) == 6 and Fr(2, 1) * (4 - Fr(4, 2)) == 4  # e h = 2 (4 - 1) = 6 ? voir ci-dessous
assert 2 * (4 - 1) == 6  # h = c - a/e = 4 - 1 = 3, p = e h = 6 (directrice x = a/e = 1)
# contre-exemples : x^2/4 + y^2/9 : foyers (0, ±sqrt5), e = sqrt5/3 ; x^2/16 - y^2/9 : c = 5
assert sp.sqrt(9 - 4) / 3 == sp.sqrt(5) / 3 and sp.sqrt(16 + 9) == 5
# trous : pour e dans ]0,1[, chaque leurre diffère de la bonne valeur
for E in [sp.Rational(1, 3), sp.Rational(1, 2), 1 / sp.sqrt(2), sp.Rational(3, 5), sp.Rational(9, 10)]:
    P_ = sp.Integer(1)
    A_ = P_ / (1 - E ** 2)
    B_ = P_ / sp.sqrt(1 - E ** 2)
    C_ = E * A_
    for le in [P_ / (1 + E ** 2), P_ / sp.sqrt(1 - E ** 2), P_ / (1 - E)]:
        assert not egal(le, A_)  # trou 1
    assert egal(C_ / A_, E)
    for le in [C_ / B_, A_ / C_, C_ ** 2 / A_ ** 2]:
        assert not egal(le, E)  # trou 2
    assert egal(B_ ** 2 / A_, P_)
    for le in [A_ ** 2 / B_, B_ ** 2 / C_, 2 * B_ ** 2 / A_]:
        assert not egal(le, P_)  # trou 3
# affirmations : x^2/4 + y^2 et x^2/16 + y^2/4 ont e = sqrt3/2 ; hyperbole équilatère e = sqrt2
assert sp.sqrt(4 - 1) / 2 == sp.sqrt(16 - 4) / 4 == sp.sqrt(3) / 2 and sp.sqrt(1 + 1) == sp.sqrt(2)

# ---------------------------------------------------------------------------------------------------------------
# Carte geo-l3-coniques-classification-affine
# ---------------------------------------------------------------------------------------------------------------


def invariants(P):
    """Rang et signature (non ordonnée) de la partie quadratique q et de la forme homogénéisée Q."""
    Ph = sp.expand(P.subs({x: x / t, y: y / t}) * t ** 2)
    Q = sp.hessian(Ph, (x, y, t)) / 2
    q = sp.hessian(sp.Poly(P, x, y).as_expr(), (x, y)) / 2
    q = sp.Matrix([[sp.Poly(P, x, y).coeff_monomial(x ** 2), sp.Poly(P, x, y).coeff_monomial(x * y) / 2],
                   [sp.Poly(P, x, y).coeff_monomial(x * y) / 2, sp.Poly(P, x, y).coeff_monomial(y ** 2)]])

    def sig(M):
        vp = [complex(sp.N(l)).real for l in M.eigenvals(multiple=True)]
        pos = sum(1 for l in vp if l > 1e-12)
        neg = sum(1 for l in vp if l < -1e-12)
        return tuple(sorted((pos, neg), reverse=True))
    return sig(q), sig(Q)


formes = [x**2 + y**2 - 1, x**2 - y**2 - 1, x**2 - y, x**2 + y**2 + 1, x**2 + y**2, x**2 - y**2, x**2 - 1, x**2, x**2 + 1]
attendus = [((2, 0), (2, 1)), ((1, 1), (2, 1)), ((1, 0), (2, 1)), ((2, 0), (3, 0)), ((2, 0), (2, 0)), ((1, 1), (1, 1)),
            ((1, 0), (1, 1)), ((1, 0), (1, 0)), ((1, 0), (2, 0))]
inv = [invariants(f) for f in formes]
assert inv == attendus and len(set(inv)) == 9  # démonstration : les neuf couples sont distincts
assert sp.expand(x**2 + (y - t)**2 / 4 - (y + t)**2 / 4 - (x**2 - y * t)) == 0
# exemples
assert sp.expand((x + y/2)**2 + sp.Rational(3, 4) * y**2 - 1 - (x**2 + x*y + y**2 - 1)) == 0
assert sp.expand(((x + y)/2)**2 - ((x - y)/2)**2 - x*y) == 0
assert sp.expand((x - y)**2 + x - (x**2 - 2*x*y + y**2 + x)) == 0
# contre-exemples et affirmations : x^2+y^2+1 et x^2+1 (rangs 3 et 2), x^2 - 2xy + y^2 - 1 = (x-y)^2 - 1
assert invariants(x**2 + y**2 + 1)[1] == (3, 0) and invariants(x**2 + 1)[1] == (2, 0)
assert sp.factor(x**2 - 2*x*y + y**2 - 1) == (x - y - 1) * (x - y + 1)
# trous : orthonormé, translation, rotation ne suffisent pas (x^2 + 2y^2 = 1 n'est pas équivalente à x^2 + y^2 = 1 par
# isométrie et facteur : le rapport des valeurs propres de q est conservé) ; xy - 1 n'a pas de partie quadratique
# proportionnelle à x^2 - y^2 (translation : q inchangée)
assert not any(sp.simplify(lam * (x**2 - y**2) - x*y) == 0 for lam in [1, -1, 2, sp.Rational(1, 2)])

# ---------------------------------------------------------------------------------------------------------------
# Carte geo-l3-coniques-polaire
# ---------------------------------------------------------------------------------------------------------------
# r = p/(1 + e cos θ) : r = p - e x, donc |r| = |p - e x| = e |p/e - x|
r = p / (1 + e * sp.cos(th))
assert sp.simplify(r - (p - e * r * sp.cos(th))) == 0
# exemple : r = 2/(1 + cos θ) : y^2 = 4 - 4x
for T in [sp.pi / 3, sp.Rational(1, 2), 2]:
    R = 2 / (1 + sp.cos(T))
    X, Y = R * sp.cos(T), R * sp.sin(T)
    assert sp.simplify(Y**2 - (4 - 4 * X)) == 0
# exemple : e = 1/2, p = 3 : r(0) = 2, r(π) = 6, a = 4, c = 2, b^2 = p a = 12
assert Fr(3) / Fr(3, 2) == 2 and Fr(3) / Fr(1, 2) == 6 and 4 * 4 - 2 * 2 == 12 == 3 * 4
# contre-exemples : (6,0) sur la directrice x = 6 ; e = 2, p = 4 : θ = π donne r = -4, point (4,0), OM = 4 = 2·(4-2)
assert Fr(3) / (1 - Fr(1, 2)) == 6 and Fr(4) / (1 - 2) == -4 and 4 == 2 * (4 - 2) and Fr(4, 2) < 4
# pôle au centre : 1 + e = 1 - e impose e = 0
assert sp.solve(sp.Eq(1 + e, 1 - e), e) == []
# trous : pour e < 1, les leurres ne donnent pas les bons points (θ = 0 donne p/(1+e))
for E in [sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(4, 5)]:
    P_ = 1
    vrai = P_ / (1 + E)
    assert P_ / (1 - E) != vrai and P_ / (1 + 0) != vrai and P_ * (1 + E) != vrai  # trou 1 (θ = 0)
    assert P_ / (1 + E * sp.cos(sp.pi)) == P_ / (1 - E)  # θ = π a cos θ < 0 : exclu par le leurre « cos θ > 0 » (trou 2)
    vraies = {P_ / (1 + E), P_ / (1 - E)}
    for le in [{P_ * (1 - E), P_ * (1 + E)}, {sp.Rational(P_, 2), P_}, {P_ / (1 + E**2), P_ / (1 - E**2)}]:
        assert le != vraies  # trou 3
# trou 3, leurre p/(1+e) et p/e : coïnciderait pour e = 1/2, c'est pourquoi il n'est pas utilisé
assert sp.solve(sp.Eq(1 / (1 - e), 1 / e), e) == [sp.Rational(1, 2)]
# affirmation : p/(1+e) + p/(1-e) = 2p/(1-e^2)
assert sp.simplify(p / (1 + e) + p / (1 - e) - 2 * p / (1 - e**2)) == 0

# ---------------------------------------------------------------------------------------------------------------
# Carte geo-l3-coniques-tangentes
# ---------------------------------------------------------------------------------------------------------------
# démonstration (ellipse) : d/dt (MF + MF') = 0 le long de (a cos t, b sin t), et M' orthogonal à u + u'
A_, B_ = 5, 3
C_ = 4
M = sp.Matrix([A_ * sp.cos(t), B_ * sp.sin(t)])
F, Fp = sp.Matrix([C_, 0]), sp.Matrix([-C_, 0])
s = (M - F).norm() + (M - Fp).norm()
for T in [0.3, 1.1, 2.5, 4.0]:
    assert abs(sp.N(s.subs(t, T)) - 10) < 1e-12
    Mv = M.subs(t, T)
    uu = (F - Mv) / (F - Mv).norm()
    uup = (Fp - Mv) / (Fp - Mv).norm()
    assert abs(sp.N(M.diff(t).subs(t, T).dot(uu + uup))) < 1e-12
# exemple : M = (3, 12/5), tangente 9x + 20y = 75, u + u' = -(30/481)(9, 20)
Mv = sp.Matrix([3, sp.Rational(12, 5)])
assert sp.Rational(9, 25) + sp.Rational(144, 25) / 9 == 1 and 9 * 3 + 20 * sp.Rational(12, 5) == 75
uu = (F - Mv) / (F - Mv).norm()
uup = (Fp - Mv) / (Fp - Mv).norm()
assert (F - Mv).norm() == sp.Rational(13, 5) and (Fp - Mv).norm() == sp.Rational(37, 5)
assert uu + uup == -sp.Rational(30, 481) * sp.Matrix([9, 20])
# exemple hyperbole x^2/4 - y^2/5 = 1, M = (3, 5/2) : tangente 3x - 2y = 4 dirigée par (2,3), u + u' = -(6/13)(2,3)
Mv = sp.Matrix([3, sp.Rational(5, 2)])
F, Fp = sp.Matrix([3, 0]), sp.Matrix([-3, 0])
uu = (F - Mv) / (F - Mv).norm()
uup = (Fp - Mv) / (Fp - Mv).norm()
assert uu + uup == -sp.Rational(6, 13) * sp.Matrix([2, 3]) and 3 * 3 - 2 * sp.Rational(5, 2) == 4
assert uu - uup == sp.Matrix([sp.Rational(12, 13), -sp.Rational(8, 13)])  # contre-exemple : bissectrice extérieure = normale
# contre-exemple (piège MF^2 + MF'^2 = 2a^2) : cercle unité, M = (3/5, 4/5), u - u' ∝ (3,-1), tangente ∝ (4,-3)
Mv = sp.Matrix([sp.Rational(3, 5), sp.Rational(4, 5)])
F, Fp = sp.Matrix([1, 0]), sp.Matrix([-1, 0])
assert (Mv - F).norm()**2 + (Mv - Fp).norm()**2 == 2 * 2  # 2 a^2 avec a = sqrt2
uu = (F - Mv) / (F - Mv).norm()
uup = (Fp - Mv) / (Fp - Mv).norm()
w = uu - uup
assert sp.simplify(w[0] * 1 + w[1] * 3) == 0  # w ∝ (3,-1)
assert sp.simplify(w[0] * (-3) - w[1] * 4) != 0  # non colinéaire à (4,-3)
# piège MF·MF' = a^2 (ovale de Cassini, F = (±1,0), a^2 = 2) : en un point générique, la tangente n'est pas la
# bissectrice extérieure
G = ((x - 1)**2 + y**2) * ((x + 1)**2 + y**2) - 4
Y0 = sp.nsolve(G.subs(x, sp.Rational(1, 2)), y, 1)
Mv = sp.Matrix([0.5, Y0])
grad = sp.Matrix([G.diff(x), G.diff(y)]).subs({x: 0.5, y: Y0})
tang = sp.Matrix([-grad[1], grad[0]])
F, Fp = sp.Matrix([1, 0]), sp.Matrix([-1, 0])
uu = (F - Mv) / (F - Mv).norm()
uup = (Fp - Mv) / (Fp - Mv).norm()
w = uu - uup
assert abs(tang[0] * w[1] - tang[1] * w[0]) > 1e-3
# contre-exemple parabole y^2 = 4x, M = (1,2) : tangente y = x + 1, u + v = (-1,-1)
assert 2 * 2 == 2 * (1 + 1)  # y y0 = 2 (x + x0) en (1,1) ... tangente 2y = 2(x+1)
# affirmation : tangente xx0/a^2 + yy0/b^2 = 1 orthogonale à M'(t)
assert sp.simplify(sp.Matrix([sp.cos(t) / a, sp.sin(t) / b]).dot(sp.Matrix([-a * sp.sin(t), b * sp.cos(t)]))) == 0
# démonstration (parabole) : u - v = HF/MF
M = sp.Matrix([t**2 / (2 * p), t])
F = sp.Matrix([p / 2, 0])
H = sp.Matrix([-p / 2, t])
MFn = t**2 / (2 * p) + p / 2
assert sp.simplify((M - F).dot(M - F) - MFn**2) == 0
assert sp.simplify(M.diff(t).dot((F - H))) == 0  # la tangente est orthogonale à FH

# ---------------------------------------------------------------------------------------------------------------
# Carte geo-l3-coniques-sections-cone
# ---------------------------------------------------------------------------------------------------------------
sq = sp.sqrt(1 + m**2)
pt = (u / sq, v, h + m * u / sq)
eq = sp.expand(pt[0]**2 + pt[1]**2 - k**2 * pt[2]**2)
A_ = (1 - k**2 * m**2) / (1 + m**2)
assert sp.simplify(eq - (A_ * u**2 - 2 * k**2 * m * h / sq * u + v**2 - k**2 * h**2)) == 0
# 1 - A = m^2 (1 + k^2)/(1 + m^2) = cos^2 β / cos^2 α avec cos^2 β = m^2/(1+m^2), cos^2 α = 1/(1+k^2)
assert sp.simplify(1 - A_ - (m**2 / (1 + m**2)) / (1 / (1 + k**2))) == 0
# C = k^2 h^2 + k^4 m^2 h^2/(s^2 A) = k^2 h^2/(s^2 A)
assert sp.simplify(k**2 * h**2 + k**4 * m**2 * h**2 / (sq**2 * A_) - k**2 * h**2 / (sq**2 * A_)) == 0
# plan x = d : y^2 + d^2 = k^2 z^2, e^2 = 1 + k^2 = 1/cos^2 α
# exemple : x^2 + y^2 = 3z^2, x = 1 : e = 2 ; plan z = x/sqrt3 + 1 : k m = 1 (parabole)
assert sp.sqrt(1 + 3) == 2 and sp.sqrt(3) * (1 / sp.sqrt(3)) == 1
# contre-exemples : x^2 + y^2 = z^2, z = 2x + 1 : A = -3/5 ; z = 2x (par le sommet) : deux droites v = ±sqrt(3/5) u
assert A_.subs({k: 1, m: 2}) == -sp.Rational(3, 5)
assert sp.N(sp.atan(2) * 180 / sp.pi) > 63.4 and sp.N(sp.atan(sp.Rational(1, 2)) * 180 / sp.pi) < 26.6
# affirmation : α = π/4, β = π/3 : ellipse, e = sqrt2/2
assert sp.simplify(sp.cos(sp.pi / 3) / sp.cos(sp.pi / 4) - sp.sqrt(2) / 2) == 0
# trous : β = 2α > α (ellipse si 2α <= π/2) ; leurres d'excentricité pour une ellipse (β > α, cercle compris) :
# cos α / cos β > 1, sin β / sin α > 1 et 1/cos α > 1, alors que cos β / cos α < 1
for (al, be) in [(sp.pi / 6, sp.pi / 3), (sp.pi / 8, sp.pi / 4), (sp.pi / 5, sp.pi / 2), (sp.pi / 3, 2 * sp.pi / 5)]:
    vrai = sp.cos(be) / sp.cos(al)
    assert sp.N(vrai) < 1 and sp.N(sp.sin(be) / sp.sin(al)) > 1 and sp.N(1 / sp.cos(al)) > 1
    if be != sp.pi / 2:
        assert sp.N(sp.cos(al) / sp.cos(be)) > 1

# ---------------------------------------------------------------------------------------------------------------
# Carte geo-l3-coniques-kepler
# ---------------------------------------------------------------------------------------------------------------
T_ = sp.symbols('T', real=True)
uf = sp.Function('u')
C, K = sp.symbols('C K', positive=True)
# Binet : avec r = 1/u(θ(t)), θ' = C u^2 : r' = -C u', r'' = -C^2 u^2 u'' et r'' - r θ'^2 = -C^2 u^2 (u + u'')
U = uf(th)
rp = sp.diff(1 / U, th) * C * U**2
assert sp.simplify(rp + C * sp.diff(U, th)) == 0
rpp = sp.diff(rp, th) * C * U**2
assert sp.simplify(rpp + C**2 * U**2 * sp.diff(U, th, 2)) == 0
assert sp.simplify(rpp - (1 / U) * (C * U**2)**2 + C**2 * U**2 * (U + sp.diff(U, th, 2))) == 0
assert sp.simplify(rp**2 + (1 / U)**2 * (C * U**2)**2 - C**2 * (U**2 + sp.diff(U, th)**2)) == 0
# cas newtonien : u = 1/p + A cos(θ - θ0) vérifie u'' + u = K/C^2 = 1/p
A1, t0 = sp.symbols('A1 t0', real=True)
usol = 1 / p + A1 * sp.cos(th - t0)
assert sp.simplify(sp.diff(usol, th, 2) + usol - 1 / p) == 0
# troisième loi : T^2 = 4π^2 a^2 b^2/C^2 avec b^2 = p a, C^2 = K p
assert sp.simplify(4 * sp.pi**2 * a**2 * (p * a) / (K * p) - 4 * sp.pi**2 * a**3 / K) == 0
# exemple : orbite circulaire v^2 = K/R ; Jupiter 5,2^(3/2) ≈ 11,9
assert abs(5.2**1.5 - 11.86) < 0.01
# contre-exemples : x = cos t, y = 2 sin t : a = -OM, C = 2, foyers (0, ±sqrt3) ; spirale r = 1/θ (u'' = 0 si K = C^2)
Mh = sp.Matrix([sp.cos(t), 2 * sp.sin(t)])
assert sp.simplify(Mh.diff(t, 2) + Mh) == sp.zeros(2, 1)
assert sp.simplify(Mh[0] * Mh.diff(t)[1] - Mh[1] * Mh.diff(t)[0]) == 2 and sp.sqrt(4 - 1) == sp.sqrt(3)
assert sp.diff(th, th, 2) == 0
# trous (kepler) : pour l'ellipse e = 1/2, p = 3/2, K = 1, la vitesse et <r,v> ne sont pas constantes ; v ∧ a non plus
pK, eK = sp.Rational(3, 2), sp.Rational(1, 2)
CK = sp.sqrt(pK)
uK = (1 + eK * sp.cos(th)) / pK
v2 = CK**2 * (uK**2 + sp.diff(uK, th)**2)
assert v2.subs(th, 0) != v2.subs(th, sp.pi)
rv = (1 / uK) * (-CK * sp.diff(uK, th))  # <r, v> = r r'
assert sp.simplify(rv.subs(th, sp.pi / 2)) != 0 and rv.subs(th, 0) == 0
va = K * CK * uK**3  # |v ∧ a| = (K/r^2) r θ' = K C u^3 (non constant)
assert va.subs({th: 0, K: 1}) != va.subs({th: sp.pi, K: 1})
# trou 2 : u'' + 2u = 1/p et u'' - u = 1/p sont incompatibles avec u'' + u = 1/p (u > 0) ; trou 3 : 2C^2/K, C^2/(2K), C^2/(4π^2 K) != C^2/K

# ---------------------------------------------------------------------------------------------------------------
# Exercices des cartes definitions et reduites
# ---------------------------------------------------------------------------------------------------------------


def ell_bifocale(cc, aa):
    return (aa**2, aa**2 - cc**2)


# definitions-ex1 (et v1, v2) : (a^2, b^2) ; leurres (a^2, c^2), (a^2, a^2+c^2), (4a^2, 4b^2)
for (cc, aa, bonne, leurres) in [(4, 5, (25, 9), [(25, 16), (25, 41), (100, 36)]),
                                 (12, 13, (169, 25), [(169, 144), (169, 313), (676, 100)]),
                                 (3, 5, (25, 16), [(25, 9), (25, 34), (100, 64)])]:
    assert ell_bifocale(cc, aa) == bonne and all(l != bonne for l in leurres)
    assert bifocal_ellipse_ok(float(cc), float(aa), 20)
# definitions-ex2 (et v1, v2) : hyperbole, b^2 = c^2 - a^2
for (cc, aa, bonne, dirs) in [(5, 3, (9, 16), Fr(9, 5)), (13, 5, (25, 144), Fr(25, 13)), (10, 6, (36, 64), Fr(18, 5))]:
    assert (aa**2, cc**2 - aa**2) == bonne and Fr(aa**2, cc) == dirs
    assert bifocal_hyperbole_ok(float(cc), float(aa), 20)
    assert (aa**2, cc**2 + aa**2) != bonne and (4 * aa**2, 4 * (cc**2 - aa**2)) != bonne
# definitions-ex3 (et v1, v2) : foyer O, directrice x = d, excentricité ε : équation, centre, second foyer


def conique_foyer_origine(d, eps):
    eqc = sp.expand(x**2 + y**2 - eps**2 * (x - d)**2)
    return eqc


for (d, eps, centre, aa, bb2, Fp) in [(3, sp.Rational(1, 2), -1, 2, 3, -2), (4, sp.Rational(1, 3), sp.Rational(-1, 2), sp.Rational(3, 2), 2, -1),
                                      (6, sp.Rational(1, 2), -2, 4, 12, -4)]:
    eqc = conique_foyer_origine(d, eps)
    red = (x - centre)**2 / aa**2 + y**2 / bb2 - 1
    rapport = sp.cancel(eqc / sp.expand(red))
    assert rapport.free_symbols == set() and rapport != 0  # même équation à un facteur près
    cc = sp.sqrt(aa**2 - bb2)
    assert centre + cc == 0 and centre - cc == Fp and cc / aa == eps and centre + aa / eps == d
# leurres de definitions-ex3 : centre, symétrique de la directrice, opposé du second foyer
assert len({(-2, 0), (-1, 0), (-3, 0), (2, 0)}) == 4

# definitions-ex4 : demi-grand axe R/2 (MF + MF' = R)

# reduites-ex1 (et v1, v2) : e = c/a, p = b^2/a ; leurres
for (A2, B2, eb, pb, leurres) in [(25, 16, Fr(3, 5), Fr(16, 5), [(Fr(3, 5), Fr(25, 3)), (Fr(4, 5), Fr(16, 5)), (Fr(3, 4), Fr(9, 4))]),
                                  (169, 144, Fr(5, 13), Fr(144, 13), [(Fr(5, 13), Fr(169, 5)), (Fr(12, 13), Fr(144, 13)), (Fr(5, 12), Fr(25, 12))])]:
    aa, bb = sp.sqrt(A2), sp.sqrt(B2)
    cc = sp.sqrt(A2 - B2)
    assert cc / aa == sp.Rational(eb.numerator, eb.denominator) and B2 / aa == sp.Rational(pb.numerator, pb.denominator)
    assert all(l != (eb, pb) for l in leurres)
    # contrôle : e × distance du foyer à la directrice = p
    assert (cc / aa) * (aa**2 / cc - cc) == B2 / aa
# reduites-ex1 v2 : x^2/4 + y^2/3 : e = 1/2, p = 3/2 ; leurres (1/2, 4), (sqrt3/2, 3/2), (1/sqrt3, 1/sqrt3)
assert sp.sqrt(4 - 3) / 2 == sp.Rational(1, 2) and sp.Rational(3, 2) == sp.Rational(3, 2)
assert {(sp.Rational(1, 2), 4), (sp.sqrt(3) / 2, sp.Rational(3, 2)), (1 / sp.sqrt(3), 1 / sp.sqrt(3))}.isdisjoint({(sp.Rational(1, 2), sp.Rational(3, 2))})
# reduites-ex2 (et v1, v2) : hyperbole, e = c/a, directrices ±a^2/c
for (A2, B2, eb, db) in [(16, 9, sp.Rational(5, 4), sp.Rational(16, 5)), (9, 16, sp.Rational(5, 3), sp.Rational(9, 5)), (4, 5, sp.Rational(3, 2), sp.Rational(4, 3))]:
    aa, cc = sp.sqrt(A2), sp.sqrt(A2 + B2)
    assert cc / aa == eb and A2 / cc == db
    # leurres : c = sqrt|a^2 - b^2| et directrices a^2/c' ; directrices c^2/a ; échange de a et b
    c2 = sp.sqrt(abs(A2 - B2))
    assert c2 / aa != eb and cc**2 / aa != db and sp.sqrt(B2 + A2) / sp.sqrt(B2) != eb
# reduites-ex3 (et v1, v2) : foyer O, directrice x = -d, excentricité ε : a = p/|1 - ε^2| ; leurres p/(1+ε^2), p, c
for (d, eps, abonne, leurres) in [(2, 2, sp.Rational(4, 3), [sp.Rational(4, 5), 4, sp.Rational(8, 3)]),
                                  (3, sp.Rational(1, 2), 2, [sp.Rational(6, 5), sp.Rational(3, 2), 1]),
                                  (1, 3, sp.Rational(3, 8), [sp.Rational(3, 10), 3, sp.Rational(9, 8)])]:
    eps = sp.nsimplify(eps)
    pp = eps * d
    assert pp / abs(1 - eps**2) == abonne and pp / (1 + eps**2) == leurres[0] and pp == leurres[1] and eps * abonne == leurres[2]
    eqc = sp.expand(x**2 + y**2 - eps**2 * (x + d)**2)
    xc = sp.solve(sp.diff(eqc, x), x)[0]
    const = -eqc.subs({x: xc, y: 0})  # coefficient constant après translation
    coefx = eqc.coeff(x, 2)
    assert sp.simplify(const / coefx - abonne**2) == 0  # a^2 = (second membre)/(coefficient de X^2)
    assert abonne not in leurres
# reduites-ex4 (et v1, v2) : excentricité par valeurs propres
for (S, k_, ebonne, leurres) in [(sp.Matrix([[5, 3], [3, 5]]), 8, sp.sqrt(3) / 2, [sp.Rational(1, 2), sp.Rational(3, 4), sp.sqrt(5) / 2]),
                                 (sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]]), 1, sp.sqrt(6) / 3, [sp.sqrt(3) / 3, sp.Rational(2, 3), 2 * sp.sqrt(3) / 3]),
                                 (sp.Matrix([[3, 5], [5, 3]]), 8, sp.sqrt(5), [sp.sqrt(3), 5, sp.sqrt(5) / 2])]:
    l1, l2 = sorted(S.eigenvals(multiple=True), reverse=True)
    if l1 * l2 > 0:
        demi = sorted([sp.sqrt(k_ / l1), sp.sqrt(k_ / l2)], reverse=True)
        ecc = sp.sqrt(demi[0]**2 - demi[1]**2) / demi[0]
    else:
        aa, bb = sp.sqrt(k_ / l1), sp.sqrt(-k_ / l2)
        ecc = sp.sqrt(aa**2 + bb**2) / aa
    assert sp.simplify(ecc - ebonne) == 0 and all(not egal(l, ebonne) for l in leurres)
# reduites-ex5 : rapport p'/p (homothétie : le paramètre est multiplié par k)

# ---------------------------------------------------------------------------------------------------------------
# Exercices de la carte classification-affine
# ---------------------------------------------------------------------------------------------------------------
# ex1 : parabole ; v1 : deux droites parallèles ; v2 : point (1,-1)
assert invariants(x**2 - 2*x*y + y**2 + x) == ((1, 0), (2, 1))
assert sp.factor(x**2 - 2*x*y + y**2 + x - y) == (x - y) * (x - y + 1) and invariants(x**2 - 2*x*y + y**2 + x - y) == ((1, 0), (1, 1))
assert sp.expand((x + y)**2 + (x - 1)**2 - (2*x**2 + 2*x*y + y**2 - 2*x + 1)) == 0
assert invariants(2*x**2 + 2*x*y + y**2 - 2*x + 1) == ((2, 0), (2, 0))
# ex2 (et v1, v2) : centre et nature
for (P_, cb, leurres) in [(x**2 + 4*x*y + 5*y**2 - 2*x - 2*y, (3, -1), [(1, -1), (-1, 3), (1, 1)]),
                          (x**2 + 2*x*y + 3*y**2 - 4*x - 8*y, (1, 1), [(-1, -1), (2, 4), (2, 0)]),
                          (x**2 - 2*x*y + 2*y**2 + 2*x - 6*y, (1, 2), [(-1, -2), (2, 1), (-1, 3)])]:
    sol = sp.solve([P_.diff(x), P_.diff(y)], [x, y])
    assert (sol[x], sol[y]) == cb and all(l != cb for l in leurres)
    assert invariants(P_) == ((2, 0), (2, 1))  # ellipse non vide
assert (x**2 + 4*x*y + 5*y**2 - 2*x - 2*y).subs({x: 3, y: -1}) == -2
assert (x**2 + 2*x*y + 3*y**2 - 4*x - 8*y).subs({x: 1, y: 1}) == -6
assert (x**2 - 2*x*y + 2*y**2 + 2*x - 6*y).subs({x: 1, y: 2}) == -5
# ex3 : rangs 3 et 2
assert invariants(x**2 + y**2 + 1)[1] == (3, 0) and invariants(x**2 + 1)[1] == (2, 0)
assert invariants(x**2 + y**2)[1] == (2, 0) and invariants(x**2)[1] == (1, 0)
# ex4 (et v1, v2) : discussion selon λ


def nature(P_):
    return {((2, 0), (2, 1)): 'ellipse', ((1, 1), (2, 1)): 'hyperbole', ((1, 0), (2, 1)): 'parabole',
            ((2, 0), (2, 0)): 'point', ((1, 1), (1, 1)): 'secantes', ((1, 0), (1, 1)): 'paralleles',
            ((2, 0), (3, 0)): 'vide'}.get(invariants(P_), 'autre')


lam = sp.Rational
fam0 = lambda L: x**2 + 2*L*x*y + y**2 + 2*x + 1
fam1 = lambda L: x**2 + 2*x*y + L*y**2 + 2*x
fam2 = lambda L: L*x**2 + y**2 + 2*x
tests = [lam(-3), -1, lam(-1, 2), 0, lam(1, 3), lam(1, 2), 1, lam(3, 2), 2, 5]
ell0 = {L for L in tests if nature(fam0(L)) == 'ellipse'}
assert ell0 == {L for L in tests if 0 < abs(L) < 1}  # ex4 : 0 < |λ| < 1
assert nature(fam0(0)) == 'point' and nature(fam0(1)) == 'parabole' and nature(fam0(-1)) == 'parabole' and nature(fam0(2)) == 'hyperbole'
assert sp.Matrix([[1, lam(1, 2), 1], [lam(1, 2), 1, 0], [1, 0, 1]]).det() == -lam(1, 4)
ell1 = {L for L in tests if nature(fam1(L)) == 'ellipse'}
assert ell1 == {L for L in tests if L > 1}  # v1 : λ > 1
assert nature(fam1(1)) == 'parabole' and nature(fam1(0)) == 'secantes' and nature(fam1(lam(1, 2))) == 'hyperbole'
ell2 = {L for L in tests if nature(fam2(L)) == 'ellipse'}
assert ell2 == {L for L in tests if L > 0}  # v2 : λ > 0
assert nature(fam2(0)) == 'parabole' and nature(fam2(-1)) == 'hyperbole'
# ex5 : contrôle numérique sur x^2/4 + y^2 = 1 : les milieux des cordes de direction (1,1) sont alignés avec le centre
mil = []
for c0 in [-1, -0.5, 0.3, 1.2]:
    sols = sp.solve([x**2 / 4 + y**2 - 1, y - x - c0], [x, y])
    (x1, y1), (x2, y2) = sols
    mil.append(((x1 + x2) / 2, (y1 + y2) / 2))
assert all(abs(sp.N(mx * mil[0][1] - my * mil[0][0])) < 1e-12 for (mx, my) in mil)

# ---------------------------------------------------------------------------------------------------------------
# Exercices de la carte polaire
# ---------------------------------------------------------------------------------------------------------------
# ex1 (et v1, v2) : forme normalisée
for (num, cst, coef, eb, pb) in [(12, 3, 1, lam(1, 3), 4), (10, 5, 3, lam(3, 5), 2), (6, 1, 3, 3, 6)]:
    assert lam(coef, cst) == eb and lam(num, cst) == pb
assert lam(4) / (1 + lam(1, 3)) == 3 and lam(4) / (1 - lam(1, 3)) == 6 and lam(4) / (1 - lam(1, 9)) == lam(9, 2)
assert lam(2) / (1 + lam(3, 5)) == lam(5, 4) and lam(2) / (1 - lam(3, 5)) == 5 and lam(2) / (1 - lam(9, 25)) == lam(25, 8)
# ex2 (et v1, v2) : e = (R - r)/(R + r), p = 2rR/(r + R), a = (r + R)/2 ; leurres a, r + R, r
for (rm, rM, eb, pb, ab) in [(2, 6, lam(1, 2), 3, 4), (1, 9, lam(4, 5), lam(9, 5), 5), (3, 5, lam(1, 4), lam(15, 4), 4)]:
    assert lam(rM - rm, rM + rm) == eb and lam(2 * rm * rM, rm + rM) == pb and lam(rm + rM, 2) == ab
    assert pb / (1 + eb) == rm and pb / (1 - eb) == rM and pb / (1 - eb**2) == ab
    assert pb not in (ab, rm + rM, rm)
# ex3 (et v1, v2) : 1/FM + 1/FN = 2/p
for (pp, FM, FNb, bornes, leurres) in [(lam(9, 5), 3, lam(9, 7), (1, 9), [lam(7, 9), 7, 3]),
                                       (3, 4, lam(12, 5), (2, 6), [lam(5, 12), 4, 2]),
                                       (2, 5, lam(5, 4), (1, sp.oo), [lam(4, 5), 5, 1])]:
    assert 1 / (2 / lam(pp) - lam(1, FM)) == FNb and bornes[0] <= FM <= bornes[1] and bornes[0] <= FNb <= bornes[1]
    assert FNb not in leurres
# vérification géométrique de v2 : M = (4,4), N = (1/4,-1) sur y^2 = 4x et alignés avec F = (1,0)
assert sp.Matrix([[4 - 1, 4], [lam(1, 4) - 1, -1]]).det() == 0 and (-1)**2 == 4 * lam(1, 4)
# vérification de ex3 sur l'ellipse x^2/25 + y^2/9 : une corde focale avec FM = 3 existe et FN = 9/7
thm = sp.acos((lam(9, 5) / 3 - 1) / lam(4, 5))
assert sp.simplify(lam(9, 5) / (1 - lam(4, 5) * sp.cos(thm)) - lam(9, 7)) == 0
# ex4 (et v1, v2) : r > 0 ssi cos θ > -1/e ; bornes 2π/3, 3π/4, 5π/6
for (ee_, borne) in [(2, 2 * sp.pi / 3), (sp.sqrt(2), 3 * sp.pi / 4), (2 / sp.sqrt(3), 5 * sp.pi / 6)]:
    assert sp.simplify(1 + ee_ * sp.cos(borne)) == 0
    for le in [sp.pi / 3, sp.pi / 2, 2 * sp.pi / 3, 3 * sp.pi / 4, 5 * sp.pi / 6, sp.pi / 4, sp.pi / 6]:
        if le != borne:
            assert sp.simplify(1 + ee_ * sp.cos(le)) != 0
# centres et asymptotes : (2,0), a = 1, b = sqrt3 ; (sqrt2,0), a = b = 1 ; (4 sqrt3, 0), a = 6, b = 2 sqrt3
for (pp, ee_, centre, aa, bb) in [(3, 2, 2, 1, sp.sqrt(3)), (1, sp.sqrt(2), sp.sqrt(2), 1, 1), (2, 2 / sp.sqrt(3), 4 * sp.sqrt(3), 6, 2 * sp.sqrt(3))]:
    s1 = sp.Integer(pp) / (1 + ee_)  # sommet θ = 0
    s2 = -sp.Integer(pp) / (1 - ee_)  # θ = π, r < 0 : le point d'abscisse pp/(ee_ - 1)
    nul = lambda w: sp.simplify(sp.radsimp(sp.sympify(w))) == 0
    assert nul((s1 + s2) / 2 - centre) and nul((s2 - s1) / 2 - aa)
    assert nul(centre - ee_ * aa) and nul(bb / aa - sp.sqrt(ee_**2 - 1))
assert sp.simplify(2 * sp.sqrt(3) / (sp.sqrt(3) + 2) - (4 * sp.sqrt(3) - 6)) == 0
# ex5 : L(θ) = 2p/(1 - e^2 cos^2 θ) et 1/L1 + 1/L2 = (2 - e^2)/(2p) ; leurres différents
L = p / (1 + e * sp.cos(th)) + p / (1 - e * sp.cos(th))
assert sp.simplify(L - 2 * p / (1 - e**2 * sp.cos(th)**2)) == 0
S2 = sp.simplify(1 / L + 1 / L.subs(th, th + sp.pi / 2))
assert sp.simplify(S2 - (2 - e**2) / (2 * p)) == 0
for le in [2 / p, (2 - e**2) / p, (1 - e**2) / p]:
    assert sp.simplify(le - (2 - e**2) / (2 * p)) != 0

# ---------------------------------------------------------------------------------------------------------------
# Exercices de la carte tangentes
# ---------------------------------------------------------------------------------------------------------------


def tangente(Pq, M):
    gx, gy = Pq.diff(x).subs({x: M[0], y: M[1]}), Pq.diff(y).subs({x: M[0], y: M[1]})
    return sp.expand(gx * (x - M[0]) + gy * (y - M[1]))


def meme_droite(D1, D2):
    """D1 = 0 et D2 = 0 sont la même droite (polynômes de degré 1 proportionnels)."""
    return sp.Matrix([[D1.coeff(x), D1.coeff(y), D1.subs({x: 0, y: 0})],
                      [D2.coeff(x), D2.coeff(y), D2.subs({x: 0, y: 0})]]).rank() == 1


for (Pq, M, bonne, leurres) in [(x**2 / 8 + y**2 / 2 - 1, (2, 1), x + 2*y - 4, [x + 4*y - 8, 2*x - y - 3, 2*x + y - 5]),
                                (x**2 / 4 + y**2 / 3 - 1, (-1, lam(3, 2)), x - 2*y + 4, [x + 2*y - 4, 4*x + 2*y + 1, x - 2*y + 2]),
                                (x**2 / 4 - y**2 / 5 - 1, (3, lam(5, 2)), 3*x - 2*y - 4, [3*x + 2*y - 4, 2*x + 3*y - lam(27, 2), 3*x - 2*y - 2])]:
    assert Pq.subs({x: M[0], y: M[1]}) == 0
    T_ = tangente(Pq, M)
    assert meme_droite(T_, bonne) and all(not meme_droite(T_, l) for l in leurres)
# normales citées : 2x - y = 3 (ex1), 4x + 2y = -1 (v1), 2x + 3y = 27/2 (v2) passent par M et sont orthogonales
for (Pq, M, nor) in [(x**2 / 8 + y**2 / 2 - 1, (2, 1), 2*x - y - 3), (x**2 / 4 + y**2 / 3 - 1, (-1, lam(3, 2)), 4*x + 2*y + 1),
                     (x**2 / 4 - y**2 / 5 - 1, (3, lam(5, 2)), 2*x + 3*y - lam(27, 2))]:
    T_ = tangente(Pq, M)
    assert nor.subs({x: M[0], y: M[1]}) == 0 and T_.coeff(x) * nor.coeff(x) + T_.coeff(y) * nor.coeff(y) == 0
# ex2 (et v1, v2) : pied de la normale N = e^2 x0, NF/NF' = MF/MF' ; leurres x0, a^2/x0, e x0
for (A2, B2, M, Nb, leurres) in [(25, 9, (4, lam(9, 5)), lam(64, 25), [4, lam(25, 4), lam(16, 5)]),
                                 (169, 144, (5, lam(144, 13)), lam(125, 169), [5, lam(169, 5), lam(25, 13)]),
                                 (4, 3, (1, lam(3, 2)), lam(1, 4), [1, 4, lam(1, 2)])]:
    Pq = x**2 / A2 + y**2 / B2 - 1
    assert Pq.subs({x: M[0], y: M[1]}) == 0
    g = (lam(2 * M[0], A2), 2 * M[1] / B2)
    tt = -M[1] / g[1]
    assert M[0] + tt * g[0] == Nb and Nb not in leurres
    cc = sp.sqrt(A2 - B2)
    MF, MFp = dist(M, (cc, 0)), dist(M, (-cc, 0))
    assert sp.simplify((cc - Nb) / (Nb + cc) - MF / MFp) == 0
# ex3 (et v1, v2) : T = (-x0, 0), FT = FM ; leurres
for (M, Tb, leurres) in [((4, 4), (-4, 0), [(4, 0), (-2, 0), (-1, 0)]), ((9, 6), (-9, 0), [(9, 0), (-lam(9, 2), 0), (-1, 0)]),
                         ((lam(1, 4), 1), (-lam(1, 4), 0), [(lam(1, 4), 0), (-lam(1, 8), 0), (-1, 0)])]:
    T_ = tangente(y**2 - 4*x, M)
    xT = sp.solve(T_.subs(y, 0), x)[0]
    assert (xT, 0) == Tb and Tb not in leurres
    assert dist((1, 0), Tb) == dist((1, 0), M)
# ex4 : réflexion par rapport à la bissectrice extérieure (contrôle numérique) et OP = a
A_, B_ = 5, 3
for T in [0.4, 1.3, 2.2, 3.7, 5.5]:
    Mv = sp.Matrix([A_ * sp.cos(T), B_ * sp.sin(T)]).evalf()
    tg = sp.Matrix([-A_ * sp.sin(T), B_ * sp.cos(T)]).evalf()
    tg = tg / tg.norm()
    Fv = sp.Matrix([4, 0])
    w = Fv - Mv
    Pp = Mv + tg * (w.dot(tg))  # projeté de F sur la tangente
    assert abs(Pp.norm() - 5) < 1e-10
# tangentes-ex4, remarque de la correction : tout point P du cercle principal est obtenu ; la droite passant par
# P = (a cos φ, a sin φ) et orthogonale à FP (F = (c,0)) vérifie la condition de tangence a^2 u^2 + b^2 v^2 = w^2
ph_, aa_, cc_ = sp.symbols('phi_ a_ c_', positive=True)
nu_, nv_ = aa_ * sp.cos(ph_) - cc_, aa_ * sp.sin(ph_)
nw_ = nu_ * aa_ * sp.cos(ph_) + nv_ * aa_ * sp.sin(ph_)
assert sp.simplify(aa_**2 * nu_**2 + (aa_**2 - cc_**2) * nv_**2 - nw_**2) == 0
# ex5 : condition de tangence a^2 u^2 + b^2 v^2 = w^2 et trace nulle sur le cercle x^2 + y^2 = a^2 + b^2
U_, V_, W_ = sp.symbols('U V W', real=True)
ligne_y = sp.solve(U_ * x + V_ * y - W_, y)[0]
disc = sp.discriminant(sp.expand((x**2 / a**2 + ligne_y**2 / b**2 - 1) * a**2 * b**2 * V_**2), x)
# tangentes-ex5, question 1 : pour V != 0, la droite coupe l'ellipse en un seul point ssi le discriminant est nul
assert sp.expand(disc - 4 * V_**2 * a**2 * b**2 * (a**2 * U_**2 + b**2 * V_**2 - W_**2)) == 0
x0_, y0_ = sp.symbols('x0 y0', real=True)
phi = sp.expand((U_ * x0_ + V_ * y0_)**2 - a**2 * U_**2 - b**2 * V_**2)
Smat = sp.Matrix([[phi.coeff(U_, 2), phi.coeff(U_, 1).coeff(V_, 1) / 2], [phi.coeff(U_, 1).coeff(V_, 1) / 2, phi.coeff(V_, 2)]])
assert sp.simplify(Smat.trace() - (x0_**2 + y0_**2 - a**2 - b**2)) == 0
# exemple numérique : a = 2, b = 1, P = (sqrt5, 0) sur le cercle orthoptique : les deux pentes des tangentes ont un produit -1
sl = sp.symbols('sl')
P0 = (sp.sqrt(5), 0)
cond = sp.expand((sl * P0[0] - P0[1])**2 - (4 * sl**2 + 1))  # droite y = sl (x - x0) + y0 tangente ssi c^2 = a^2 sl^2 + b^2
ra = sp.solve(cond, sl)
assert len(ra) == 2 and sp.simplify(ra[0] * ra[1] + 1) == 0
# leurres du rayon : a + b, sqrt(a^2 - b^2), sqrt2 a diffèrent de sqrt(a^2 + b^2) pour a > b > 0
for (A_, B_) in [(2, 1), (5, 3), (3, 2)]:
    r0 = sp.sqrt(A_**2 + B_**2)
    assert all(not egal(l, r0) for l in [A_ + B_, sp.sqrt(A_**2 - B_**2), sp.sqrt(2) * A_])

# ---------------------------------------------------------------------------------------------------------------
# Exercices de la carte sections-cone
# ---------------------------------------------------------------------------------------------------------------
# ex1 (et v1, v2) : e = cos β / cos α
for (al, be, nat, eb) in [(sp.pi / 3, sp.pi / 4, 'hyperbole', sp.sqrt(2)), (sp.pi / 6, sp.pi / 3, 'ellipse', sp.sqrt(3) / 3),
                          (sp.pi / 4, sp.pi / 6, 'hyperbole', sp.sqrt(6) / 2)]:
    assert sp.simplify(sp.cos(be) / sp.cos(al) - eb) == 0
    assert (nat == 'ellipse') == (be > al)
    assert not egal(sp.cos(al) / sp.cos(be), eb) and not egal(sp.sin(be) / sp.sin(al), eb)
# leurres explicites : ex1 √2/2 et √6/3 ; v1 √3 et 1/2 ; v2 √6/3 et √2 ne sont pas les bonnes excentricités
assert not egal(sp.sqrt(2) / 2, sp.sqrt(2)) and not egal(sp.Rational(1, 2), sp.sqrt(3) / 3) and not egal(sp.sqrt(2), sp.sqrt(6) / 2)


def section(mm, hh=1, kk=1):
    """Nature et excentricité de la section de x^2 + y^2 = k^2 z^2 par z = m x + h."""
    mm, kk = sp.nsimplify(mm), sp.nsimplify(kk)
    A_ = (1 - kk**2 * mm**2) / (1 + mm**2)
    if A_ == 0:
        return 'parabole', 1
    return ('ellipse' if A_ > 0 else 'hyperbole'), sp.sqrt(1 - A_)


# ex2 (et v1, v2)
assert section(lam(1, 2)) == ('ellipse', sp.sqrt(10) / 5)
assert section(2) == ('hyperbole', 2 * sp.sqrt(10) / 5)
assert section(1) == ('parabole', 1)
assert sp.expand((x + 1)**2 - x**2 - (2*x + 1)) == 0  # v2 : y^2 = 2x + 1
assert not egal(sp.sqrt(10) / 2, sp.sqrt(10) / 5) and not egal(2, 2 * sp.sqrt(10) / 5)
# ex3 (et v1, v2) : grand axe 2 sqrt(1+m^2)/(1-m^2), petit axe 2 sqrt(C), C = 1/((1+m^2) A)
for (mm, gb, leurres, petit) in [(lam(1, 2), 4 * sp.sqrt(5) / 3, [lam(8, 3), 2 * sp.sqrt(5) / 3, 4 * sp.sqrt(3) / 3], 4 * sp.sqrt(3) / 3),
                                 (lam(1, 3), 3 * sp.sqrt(10) / 4, [lam(9, 4), 3 * sp.sqrt(10) / 8, 3 * sp.sqrt(2) / 2], 3 * sp.sqrt(2) / 2),
                                 (lam(3, 4), lam(40, 7), [lam(32, 7), lam(20, 7), 8 * sp.sqrt(7) / 7], 8 * sp.sqrt(7) / 7)]:
    A_ = (1 - mm**2) / (1 + mm**2)
    C_ = 1 / ((1 + mm**2) * A_)
    assert sp.simplify(2 * sp.sqrt(C_ / A_) - gb) == 0 and sp.simplify(2 * sp.sqrt(C_) - petit) == 0
    P1 = (1 / (1 - mm), 0, 1 / (1 - mm))
    P2 = (-1 / (1 + mm), 0, 1 / (1 + mm))
    assert sp.simplify(dist(P1, P2) - gb) == 0 and P1[0] - P2[0] == leurres[0]
    assert all(not egal(l, gb) for l in leurres)
    # équation dans le plan : vérification de la complétion du carré
    s_ = sp.sqrt(1 + mm**2)
    eqp = sp.expand((u / s_)**2 + v**2 - (1 + mm * u / s_)**2)
    u0 = mm / (s_ * A_)
    assert sp.simplify(eqp - (A_ * (u - u0)**2 + v**2 - C_)) == 0
# ex4 : AA' = d/cos α (longueur d'un segment de génératrice entre deux plans perpendiculaires à l'axe)
# sections-cone-ex3 (et v1, v2), contrôle direct dans le plan y = 0 : extrémités du grand axe sur x = ±z
for (mm_, gb, pb_) in [(lam(1, 2), 4 * sp.sqrt(5) / 3, 4 * sp.sqrt(3) / 3), (lam(1, 3), 3 * sp.sqrt(10) / 4, 3 * sp.sqrt(2) / 2),
                       (lam(3, 4), lam(40, 7), 8 * sp.sqrt(7) / 7)]:
    xa_, xb_ = 1 / (1 - mm_), -1 / (1 + mm_)  # x = z et x = -z avec z = m x + 1
    long_ = sp.sqrt((xa_ - xb_)**2 * (1 + mm_**2))
    assert sp.simplify(long_ - gb) == 0 and sp.simplify(gb - pb_) > 0
    # e^2 = 1 - (petit/grand)^2 = cos^2 β / cos^2 α = 2 m^2/(1 + m^2)
    assert sp.simplify(1 - (pb_ / gb)**2 - 2 * mm_**2 / (1 + mm_**2)) == 0

# ---------------------------------------------------------------------------------------------------------------
# Exercices de la carte kepler
# ---------------------------------------------------------------------------------------------------------------
# ex1 (et v1, v2) : v_a = r_p v_p / r_a ; leurres v_p r_a/r_p, v_p r_p^2/r_a^2, v_p
for (rp_, vp_, ra_, vb, leurres) in [(2, 6, 3, 4, [9, lam(8, 3), 6]), (3, 10, 5, 6, [lam(50, 3), lam(18, 5), 10]),
                                     (1, 12, 4, 3, [48, lam(3, 4), 12])]:
    assert lam(rp_ * vp_, ra_) == vb
    assert leurres == [lam(vp_ * ra_, rp_), lam(vp_ * rp_**2, ra_**2), vp_] and vb not in leurres
# ex2 (et v1, v2) : T2 = (a2/a1)^(3/2)
for (q_, Tb, leurres) in [(4, 8, [4, 16, 64]), (9, 27, [9, 81, 729]), (lam(1, 4), lam(1, 8), [lam(1, 4), lam(1, 16), lam(1, 64)])]:
    assert sp.Rational(q_)**sp.Rational(3, 2) == Tb and Tb not in leurres
assert abs(365.25 / 8 - 46) < 0.5
# ex3 (et v1, v2) : u'' + u en fonction de u
uc = 1 / (2 * sp.cos(th))
assert sp.simplify(sp.diff(uc, th, 2) + uc - 8 * uc**3) == 0  # cercle : n = 5
assert sp.simplify(sp.diff(uc, th, 2) - (2 - sp.cos(th)**2) / (2 * sp.cos(th)**3)) == 0
ul = sp.exp(-th)
assert sp.simplify(sp.diff(ul, th, 2) + ul - 2 * ul) == 0  # spirale : n = 2 + 1 = 3
uk = 1 / (1 + sp.cos(th))
assert sp.simplify(sp.diff(uk, th, 2) + uk - 3 * uk**2) == 0  # cardioïde : n = 2 + 2 = 4
assert sp.simplify(sp.diff(uk, th, 2) - (sp.cos(th) * (1 + sp.cos(th)) + 2 * sp.sin(th)**2) / (1 + sp.cos(th))**3) == 0
# ex4 (et v1, v2) : C^2 = r0^2 v0^2, p = C^2/K, λ = 1/r0 - 1/p, e = |λ| p
for (v02, eb) in [(lam(3, 2), lam(1, 2)), (2, 1), (3, 2)]:
    pp = lam(v02)
    lam_ = 1 - 1 / pp
    assert abs(lam_) * pp == eb == abs(v02 - 1)
assert lam(3, 2) / (1 - lam(1, 4)) == 2 and lam(3, 2) / (1 + lam(1, 2)) == 1
# leurres d'excentricité : 3/4, 3/2 (ex4) ; 1/2, 2 (v1) ; 1/2, 3 (v2) différents des bonnes valeurs
# kepler-ex4 v2 : u = 1/3 + (2/3) cos θ, soit r = 3/(1 + 2 cos θ), branche cos θ > -1/2
assert sp.simplify(1 / (lam(1, 3) + lam(2, 3) * sp.cos(th)) - 3 / (1 + 2 * sp.cos(th))) == 0
# ex5 : E = K(e^2 - 1)/(2p) et a = -K/(2E)
ph = sp.symbols('phi', real=True)
uu_ = (1 + e * sp.cos(ph)) / p
E_ = K * p / 2 * (uu_**2 + sp.diff(uu_, ph)**2) - K * uu_
assert sp.simplify(E_ - K * (e**2 - 1) / (2 * p)) == 0
assert sp.simplify(-K / (2 * (K * (e**2 - 1) / (2 * p))) - p / (1 - e**2)) == 0
for le in [4 * sp.pi**2 * a**2 / K, 2 * sp.pi * a**3 / K]:
    assert sp.simplify(le - 4 * sp.pi**2 * a**3 / K) != 0

# ---------------------------------------------------------------------------------------------------------------
# QCM
# ---------------------------------------------------------------------------------------------------------------
assert sp.sqrt(5**2 - 3**2) == 4 and 4 not in (8, sp.sqrt(34), 5)  # q1
assert sp.sqrt(9 + 16) / 3 == lam(5, 3) and lam(5, 3) not in (lam(5, 4), sp.sqrt(7) / 3, lam(4, 3))  # q3
rq = 4 / (2 - sp.cos(th))
assert rq.subs(th, sp.pi) == lam(4, 3) and rq.subs(th, 0) == 4 and lam(2) / (1 - lam(1, 4)) == lam(8, 3)  # q5
assert lam(4, 3) not in (4, 2, lam(8, 3))
assert sp.Integer(27)**2 == 9**3 and 9 not in (27, 3, 81)  # q8

print("geo_l3_coniques : toutes les vérifications passent.")
