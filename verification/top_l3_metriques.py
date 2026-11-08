"""Vérifications calculatoires du chapitre « top-l3-metriques » (topologie, L3) : cartes, exercices et QCM de
content/*/topologie-calcul-diff/l3-top-l3-metriques.json et content/qcm/top-l3-metriques.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question. Calculs
exacts avec SymPy et les fractions, valeurs numériques avec mpmath (30 chiffres), inégalités testées sur des grilles
ou par tirages aléatoires. Lancer : python3 verification/top_l3_metriques.py.
"""

import math
import random
from fractions import Fraction as F
from math import ceil

import mpmath as mp
import sympy as sp

mp.mp.dps = 30
random.seed(1)
x, y, t, u, s, c, lam, eta, n = sp.symbols("x y t u s c lambda eta n", real=True)

# ---------------------------------------------------------------------------------------------------------------
# Carte top-l3-metriques-continuite
# ---------------------------------------------------------------------------------------------------------------
# exemple et ex1 : diamètre de R pour |arctan x - arctan y| : 2 arctan n -> pi, et la borne pi n'est pas atteinte.
assert sp.limit(2 * sp.atan(n), n, sp.oo) == sp.pi
assert all(abs(mp.atan(a) - mp.atan(b)) < mp.pi for a, b in [(1e6, -1e6), (1e12, -1e12)])
# ex1 : leurres pi/2, +oo, 2pi différents de pi.
assert len({sp.pi, sp.pi / 2, sp.oo, 2 * sp.pi}) == 4
# ex1 v1 : phi(x) = x/(1+|x|) strictement croissante (dérivée 1/(1+x)^2 > 0 sur [0,+oo[, impaire), diamètre 2.
assert sp.simplify(sp.diff(x / (1 + x), x) - 1 / (1 + x) ** 2) == 0
assert sp.limit(2 * n / (1 + n), n, sp.oo) == 2
assert len({2, 1, sp.oo, 4}) == 4
# ex1 v2 : |e^0 - e^n| = e^n - 1 -> +oo, diamètre infini ; leurres 1, e, e - 1 finis.
assert sp.limit(sp.exp(n) - 1, n, sp.oo) == sp.oo

# ex2 (et v1, v2) : distance de (0,c) à la parabole y = x^2 : min_u (u + (u-c)^2) = c - 1/4, atteint en u = c - 1/2.
for cc, rep, leurres in [
    (2, sp.sqrt(7) / 2, [2, sp.Rational(7, 4), sp.sqrt(2)]),  # ex2
    (3, sp.sqrt(11) / 2, [3, sp.Rational(11, 4), sp.sqrt(3)]),  # ex2 v1
    (4, sp.sqrt(15) / 2, [4, sp.Rational(15, 4), 2]),  # ex2 v2
]:
    q = u + (u - cc) ** 2
    assert sp.expand(q - ((u - (cc - sp.Rational(1, 2))) ** 2 + cc - sp.Rational(1, 4))) == 0
    assert cc - sp.Rational(1, 2) >= 0  # le minimum est atteint pour u = t^2 >= 0
    assert sp.sqrt(cc - sp.Rational(1, 4)) == rep
    # contrôle numérique sur une grille de t
    m = min(float(tt) ** 2 + (float(tt) ** 2 - cc) ** 2 for tt in [k / 1000 for k in range(-4000, 4001)])
    assert abs(m - (cc - 0.25)) < 1e-5
    # leurres : distance au sommet (c), carré de la distance (c - 1/4), distance au point (sqrt c, c) ; tous faux
    assert leurres[0] == cc and sp.sqrt(cc) == leurres[2]
    for L in leurres:
        assert sp.simplify(L - rep) != 0
    # le point (sqrt c, c) est sur la parabole et à distance sqrt c de (0,c)
    assert sp.sqrt(sp.sqrt(cc) ** 2 + (cc - cc) ** 2) == sp.sqrt(cc)

# ex4 (question 4) : ||(n, 1/n) - (n, 0)|| = 1/n -> 0, et (n, 1/n) est sur xy = 1.
assert sp.limit(1 / n, n, sp.oo) == 0

# ex5 : min(1,|s-t|) vérifie l'inégalité triangulaire (tirages) ; sum 2^-k = 2.
rho = lambda a, b: min(1.0, abs(a - b))
for _ in range(20000):
    a, b, cq = (random.uniform(-3, 3) for _ in range(3))
    assert rho(a, cq) <= rho(a, b) + rho(b, cq) + 1e-12
assert sp.summation(sp.Rational(1, 2) ** sp.Symbol("k", integer=True, nonnegative=True), (sp.Symbol("k", integer=True, nonnegative=True), 0, sp.oo)) == 2

# contre-exemple : la fonction égale à 1 sur [0,+oo[ et 0 ailleurs ; f^{-1}(]1/2,3/2[) = [0,+oo[ contient 0 mais pas -1e-9.
ind = lambda v: 1 if v >= 0 else 0
assert 0.5 < ind(0) < 1.5 and not (0.5 < ind(-1e-9) < 1.5)

# ---------------------------------------------------------------------------------------------------------------
# Carte top-l3-metriques-compacite
# ---------------------------------------------------------------------------------------------------------------
# exemple et ex1 (v1, v2) : max de a cos t + b sin t = sqrt(a^2 + b^2), atteint en (a, b)/sqrt(a^2+b^2).
for a_, b_, rep, leurres in [
    (1, 2, sp.sqrt(5), [3, 2, 5]),  # ex1 et exemple
    (3, 4, sp.Integer(5), [7, 4, 25]),  # ex1 v1
    (1, -3, sp.sqrt(10), [4, 3, 10]),  # ex1 v2
]:
    R = sp.sqrt(a_ ** 2 + b_ ** 2)
    assert R == rep
    g = a_ * sp.cos(t) + b_ * sp.sin(t)
    pt = (sp.Integer(a_) / R, sp.Integer(b_) / R)  # point où le maximum est atteint
    assert sp.simplify(pt[0] ** 2 + pt[1] ** 2) == 1 and sp.simplify(a_ * pt[0] + b_ * pt[1] - R) == 0
    m = max(a_ * math.cos(2 * math.pi * k / 100000) + b_ * math.sin(2 * math.pi * k / 100000) for k in range(100000))
    assert abs(m - math.sqrt(a_ ** 2 + b_ ** 2)) < 1e-8
    for L in leurres:
        assert sp.simplify(L - rep) != 0
    # le leurre « somme des coefficients en valeur absolue » est bien |a| + |b|
    assert leurres[0] == abs(a_) + abs(b_)

# ex2 (v1, v2) : distance du cercle unité à la droite ax + by = c (c > sqrt(a^2+b^2)) = c/sqrt(a^2+b^2) - 1.
for a_, b_, c_, rep, leurres in [
    (1, 1, 4, 2 * sp.sqrt(2) - 1, [2 * sp.sqrt(2), 4 - sp.sqrt(2), 3]),  # ex2
    (1, 1, 3, 3 * sp.sqrt(2) / 2 - 1, [3 * sp.sqrt(2) / 2, 3 - sp.sqrt(2), 2]),  # ex2 v1
    (1, -2, 5, sp.sqrt(5) - 1, [sp.sqrt(5), 5 - sp.sqrt(5), 4]),  # ex2 v2
]:
    nr = sp.sqrt(a_ ** 2 + b_ ** 2)
    assert c_ > nr  # la droite ne rencontre pas le cercle
    assert sp.simplify((c_ - nr) / nr - rep) == 0
    # contrôle numérique : min sur le cercle de la distance à la droite
    m = min(abs(a_ * math.cos(2 * math.pi * k / 200000) + b_ * math.sin(2 * math.pi * k / 200000) - c_) / math.sqrt(a_ ** 2 + b_ ** 2)
            for k in range(200000))
    assert abs(m - float(rep)) < 1e-8
    for L in leurres:
        assert abs(sp.N(L - rep)) > 0.1

# ex4 : |e^{-i/n} - 1| = 2 sin(1/(2n)) <= 1/n.
for k in range(1, 200):
    assert abs(mp.exp(-1j / k) - 1) <= mp.mpf(1) / k

# ex5 (question 5) et carte point-fixe (contre-exemple) : h(x) = x + 1/x sur [1,+oo[.
h = x + 1 / x
assert sp.simplify((h.subs(x, y) - h) - (y - x) * (1 - 1 / (x * y))) == 0
assert sp.solve(sp.Eq(h, x), x) == []  # pas de point fixe
assert sp.minimum(h, x, sp.Interval(1, sp.oo)) == 2  # h envoie [1,+oo[ dans [2,+oo[
for _ in range(20000):
    a, b = sorted(random.uniform(1, 50) for _ in range(2))
    if b > a:
        d = (b + 1 / b) - (a + 1 / a)
        assert 0 < d < b - a

# contre-exemples : hyperbole et axe à distance 1/n ; f(x)=x sur ]0,1], f(0)=1 : image ]0,1], inf 0 non atteint.
assert sp.limit(1 / n, n, sp.oo) == 0

# ---------------------------------------------------------------------------------------------------------------
# Carte top-l3-metriques-heine
# ---------------------------------------------------------------------------------------------------------------
# contre-exemples : (n + 1/n)^2 - n^2 = 2 + 1/n^2 ; sin(1/x) avec 1/(2n pi) et 1/(2n pi + pi/2) ; sin(x^2).
assert sp.expand((n + 1 / n) ** 2 - n ** 2) == 2 + 1 / n ** 2
assert sp.sin(2 * sp.pi * 7) == 0 and sp.sin(2 * sp.pi * 7 + sp.pi / 2) == 1
xn, yn = sp.sqrt(2 * n * sp.pi + sp.pi / 2), sp.sqrt(2 * n * sp.pi)
assert sp.limit(xn - yn, n, sp.oo) == 0
assert sp.simplify(sp.sin(xn ** 2).subs(n, 5) - sp.sin(yn ** 2).subs(n, 5)) == 1
# affirmation (sqrt|x| non lipschitzienne mais uniformément continue) : |sqrt u - sqrt v| <= sqrt|u - v|.
for _ in range(20000):
    a, b = random.uniform(-100, 100), random.uniform(-100, 100)
    assert abs(abs(a) ** 0.5 - abs(b) ** 0.5) <= abs(a - b) ** 0.5 + 1e-12
# affirmation « uniformément continue sur ]0,1] => bornée » : nombre de pas ceil(2/eta) de longueur eta/2 < eta.
for et in [0.3, 0.05, 1.0]:
    N = ceil(2 / et)
    assert N * (et / 2) >= 1  # on va de 1 à n'importe quel x de ]0,1] en N pas au plus

# ex1 (v1, v2) : plus grand eta tel que sup_{|x-y|<=eta} |x^2 - y^2| <= eps sur [0,a] ; ce sup vaut 2a eta - eta^2.
for a_, eps, rep, leurres in [
    (3, 5, 1, [F(5, 6), F(5, 3), None]),  # ex1 (leurre sqrt 5)
    (5, 16, 2, [F(8, 5), F(16, 5), 4]),  # ex1 v1
    (4, 15, 3, [F(15, 8), F(15, 4), None]),  # ex1 v2 (leurre sqrt 15)
]:
    assert sp.solve(sp.Eq(2 * a_ * eta - eta ** 2, eps), eta)[0] == rep  # plus petite racine
    assert rep == a_ - sp.sqrt(a_ ** 2 - eps)
    assert a_ ** 2 > eps  # pour eta >= a, le sup a^2 dépasse eps

    def sup_ecart(e, a=a_):
        # max de x^2 - y^2 pour 0 <= y <= x <= a, x - y <= e : grille fine
        best = 0.0
        for i in range(0, 2001):
            xx = a * i / 2000
            yy = max(0.0, xx - e)
            best = max(best, xx * xx - yy * yy)
        return best

    assert abs(sup_ecart(rep) - eps) < 1e-9  # en eta = rep, la borne eps est atteinte
    assert sup_ecart(rep + 1e-3) > eps  # au-delà, elle est dépassée : rep est bien le plus grand
    # leurres : eps/(2a) et eps/a (le premier convient mais est plus petit, le second ne convient pas), sqrt(eps)
    assert leurres[0] == F(eps, 2 * a_) and leurres[0] < rep
    assert leurres[1] == F(eps, a_) and sup_ecart(float(leurres[1])) > eps
    assert sup_ecart(eps ** 0.5) > eps
    if leurres[2] is not None:
        assert leurres[2] == eps ** 0.5
# ex1 : la majoration grossière |x^2 - y^2| <= 2a|x - y| donne eta = eps/(2a), strictement plus petit.

# ex2 : f2 = x sin(1/x) -> 0 ; contre-exemples f1, f3, f4.
assert sp.limit(x * sp.sin(1 / x), x, 0, "+") == 0
assert sp.limit(1 / (2 * n * sp.pi) - 1 / (2 * n * sp.pi + sp.pi / 2), n, sp.oo) == 0
assert sp.simplify(1 / (1 / (n + 1)) - 1 / (1 / n)) == 1  # f3
assert sp.simplify(sp.log(sp.exp(-n)) - sp.log(sp.exp(-n - 1))) == 1  # f4
assert sp.limit(sp.exp(-n) - sp.exp(-n - 1), n, sp.oo) == 0

# ex3 : x^2 non uniformément continue (croissance non affine), x + sin x 2-lipschitzienne.
assert sp.limit(x ** 2 / x, x, sp.oo) == sp.oo
assert sp.maximum(sp.diff(x + sp.sin(x), x), x, sp.Interval(-10, 10)) == 2

# exemple (intégrale à paramètre) : rien de calculatoire.

# QCM q5 : x sin x avec 2n pi + 1/n et 2n pi : écart des images -> 2 pi.
assert sp.limit((2 * n * sp.pi + 1 / n) * sp.sin(1 / n), n, sp.oo) == 2 * sp.pi
assert sp.limit(sp.sqrt(2 * n * sp.pi + sp.pi / 2) - sp.sqrt(2 * n * sp.pi), n, sp.oo) == 0

# ---------------------------------------------------------------------------------------------------------------
# Carte top-l3-metriques-connexite
# ---------------------------------------------------------------------------------------------------------------
# démonstration (sinus du topologiste) : les points d'abscisse 1/(theta + 2n pi) ont pour ordonnée sin(theta).
for th in [-sp.pi / 2, sp.Rational(1, 3), sp.pi / 2]:
    for k in [1, 5, 40]:
        assert sp.simplify(sp.sin(1 / (1 / (th + 2 * k * sp.pi))) - sp.sin(th)) == 0
        assert 0 < 1 / (th + 2 * k * sp.pi) <= 1
# les abscisses 1/(pi/2 + 2k pi) et 1/(-pi/2 + 2k pi) donnent les ordonnées 1 et -1 et tendent vers 0.
assert sp.sin(sp.pi / 2 + 2 * 3 * sp.pi) == 1 and sp.sin(-sp.pi / 2 + 2 * 3 * sp.pi) == -1

# ex1 (v1, v2) : nombre d'intervalles = nombre de racines réelles distinctes + 1.
for P, rep, leurres in [
    (x ** 3 - x, 4, [3, 2, 5]),  # ex1
    (x ** 4 - 5 * x ** 2 + 4, 5, [4, 3, 6]),  # ex1 v1
    (x ** 4 - 1, 3, [2, 5, 4]),  # ex1 v2
]:
    rac = sorted(set(sp.real_roots(sp.Poly(P, x))))
    assert len(rac) + 1 == rep
    assert rep not in leurres and len(set(leurres)) == 3
assert sorted(set(sp.real_roots(sp.Poly(x ** 4 - 5 * x ** 2 + 4, x)))) == [-2, -1, 1, 2]

# ex2 : det(tA) = t^2 det A en taille 2 ; det(diag(-1,1)) = -1 ; tr(-I2) = -2.
A = sp.Matrix(2, 2, sp.symbols("a0:4"))
assert sp.expand((t * A).det() - t ** 2 * A.det()) == 0
assert sp.diag(-1, 1).det() == -1 and sp.eye(2).det() == 1
assert (-sp.eye(2)).trace() == -2 and sp.eye(2).trace() == 2
# diag(-1,1) est orthogonale (O2 contient des éléments de déterminant -1 et 1).
assert sp.diag(-1, 1).T * sp.diag(-1, 1) == sp.eye(2)

# ex3 : x -> 1/2 + arctan(x)/pi est une bijection de R sur ]0,1[, de réciproque tan(pi (t - 1/2)).
assert sp.simplify(sp.tan(sp.pi * ((sp.Rational(1, 2) + sp.atan(x) / sp.pi) - sp.Rational(1, 2))) - x) == 0
assert sp.limit(sp.Rational(1, 2) + sp.atan(x) / sp.pi, x, sp.oo) == 1
assert sp.limit(sp.Rational(1, 2) + sp.atan(x) / sp.pi, x, -sp.oo) == 0

# ex4 (question 4) : det((1-t) I2 + t diag(-1,1)) = 1 - 2t, nul en 1/2.
M = (1 - t) * sp.eye(2) + t * sp.diag(-1, 1)
assert sp.expand(M.det()) == 1 - 2 * t
# ex4 (question 2) : P(z) = det((1-z)A + zB) est un polynôme de degré <= n, P(0) = det A, P(1) = det B.
z = sp.symbols("z")
Aa = sp.Matrix([[1, 2], [3, 4]])
Bb = sp.Matrix([[0, 1], [-1, 5]])
Pz = sp.expand(((1 - z) * Aa + z * Bb).det())
assert sp.degree(Pz, z) <= 2 and Pz.subs(z, 0) == Aa.det() and Pz.subs(z, 1) == Bb.det()

# contre-exemple : x^2 et [1,4] : image réciproque [-2,-1] u [1,2].
assert set(sp.solve(sp.Eq(x ** 2, 1), x)) == {-1, 1} and set(sp.solve(sp.Eq(x ** 2, 4), x)) == {-2, 2}

# ---------------------------------------------------------------------------------------------------------------
# Carte top-l3-metriques-point-fixe
# ---------------------------------------------------------------------------------------------------------------
# exemple : cos envoie [0,1] dans [cos 1, 1] subset [0,1], |sin| <= sin 1 < 1 sur [0,1], point fixe 0,739085.
assert 0 < mp.cos(1) < 1 and mp.sin(1) < 1
pf = mp.findroot(lambda v: mp.cos(v) - v, 0.7)
assert abs(pf - mp.mpf("0.739085133215160641655312087673873")) < 1e-25
# exemple : f = e^{x/2} vérifie f(x) = 1 + (1/2) int_0^x f.
assert sp.simplify(1 + sp.Rational(1, 2) * sp.integrate(sp.exp(s / 2), (s, 0, x)) - sp.exp(x / 2)) == 0

# trous, leurres de l'estimation a priori : f(x) = x/2 + 1, x0 = 0 ; x_n = 2 - 2^{1-n}, l = 2, k = 1/2, d(x0,x1) = 1.
k_, l_ = F(1, 2), 2
xs = [F(0)]
for _ in range(12):
    xs.append(xs[-1] / 2 + 1)
for m in range(1, 12):
    err = l_ - xs[m]
    assert err == F(2, 2 ** m)
    assert err <= k_ ** m / (1 - k_) * (xs[1] - xs[0])  # vraie estimation a priori
    assert err <= k_ / (1 - k_) * (xs[m] - xs[m - 1])  # vraie estimation a posteriori
    assert err > k_ ** m * (xs[1] - xs[0])  # leurre k^n d(x0,x1) : faux
    assert err > k_ ** (m + 1) / (1 - k_) * (xs[1] - xs[0])  # leurre k^{n+1}/(1-k) d(x0,x1) : faux
    assert err > k_ ** m / (1 + k_) * (xs[1] - xs[0])  # leurre k^n/(1+k) d(x0,x1) : faux
# contre-exemple : x/2 n'a pas de point fixe dans ]0,1] (le seul est 0).
assert sp.solve(sp.Eq(x / 2, x), x) == [0]
# contre-exemple : arctan n -> pi/2 (suite de Cauchy pour la distance transportée) mais arctan l < pi/2.
assert sp.limit(sp.atan(n), n, sp.oo) == sp.pi / 2

# ex1 (v1, v2) : plus petit n tel que k^n/(1-k) |x1 - x0| <= 10^-6.
for g, kq, rep, leurres in [
    (sp.cos(x) / 2, F(1, 2), 20, [19, 21, 6]),  # ex1
    (1 + sp.sin(x) / 3, F(1, 3), 13, [12, 14, 6]),  # ex1 v1
    (1 + sp.atan(x) / 4, F(1, 4), 11, [10, 12, 6]),  # ex1 v2
]:
    # constante de Lipschitz : sup |g'| = k (atteint ou approché)
    gp = sp.lambdify(x, sp.diff(g, x), "math")
    mg = max(abs(gp(-10 + 20 * i / 100000)) for i in range(100001))
    assert mg <= float(kq) + 1e-15 and abs(mg - float(kq)) < 1e-6
    x1 = F(str(sp.nsimplify(g.subs(x, 0))))
    borne = lambda m: kq ** m / (1 - kq) * abs(x1)
    nmin = next(m for m in range(1, 100) if borne(m) <= F(1, 10 ** 6))
    assert nmin == rep
    assert borne(rep - 1) > F(1, 10 ** 6)  # le leurre rep - 1 ne suffit pas
    assert rep + 1 in leurres and 6 in leurres and rep - 1 in leurres  # rep + 1 suffit mais n'est pas le plus petit
    # contrôle : l'erreur réelle après rep itérations est bien sous 10^-6
    gf = sp.lambdify(x, g, "mpmath")
    v = mp.mpf(0)
    for _ in range(200):
        v = gf(v)
    ell = v
    v = mp.mpf(0)
    for _ in range(rep):
        v = gf(v)
    assert abs(v - ell) <= mp.mpf("1e-6")
assert 2 ** 19 == 524288 and 2 ** 20 == 1048576 and 3 ** 12 == 531441 and 3 ** 13 == 1594323
assert 4 ** 10 == 1048576 and 4 ** 11 == 4194304

# ex2 (v1, v2) : estimation a posteriori.
for g, a_, b_, kq, kk, x0, attendus, rep, leurres, ell in [
    (lambda v: 1 / (2 + v), 0, 1, F(1, 4), F(1, 3), F(0), [F(1, 2), F(2, 5), F(5, 12)], F(1, 180),
     [F(1, 60), F(1, 240), F(1, 96)], sp.sqrt(2) - 1),  # ex2
    (lambda v: 1 / (3 + v), 0, 1, F(1, 9), F(1, 8), F(0), [F(1, 3), F(3, 10), F(10, 33)], F(1, 2640),
     [F(1, 330), F(1, 2970), F(1, 1944)], (sp.sqrt(13) - 3) / 2),  # ex2 v1
    (lambda v: 1 / (1 + v), F(1, 2), 1, F(4, 9), F(4, 5), F(1), [F(1, 2), F(2, 3), F(3, 5)], F(4, 75),
     [F(1, 15), F(4, 135), F(32, 405)], (sp.sqrt(5) - 1) / 2),  # ex2 v2
]:
    # stabilité du segment (g décroissante : image [g(b), g(a)]) et constante sup|g'| = |g'(a)| = k
    assert a_ <= g(F(b_)) <= g(F(a_)) <= b_
    vv = sp.Symbol("vv")
    expr = g(vv)
    assert sp.Abs(sp.diff(expr, vv)).subs(vv, a_) == sp.Rational(kq.numerator, kq.denominator)
    assert kq / (1 - kq) == kk
    seq = [x0]
    for _ in range(3):
        seq.append(g(seq[-1]))
    assert seq[1:] == attendus
    assert kk * abs(seq[3] - seq[2]) == rep
    # leurres : |x3 - x2| seul, k |x3 - x2|, estimation a priori k^3/(1-k) |x1 - x0|
    assert leurres == [abs(seq[3] - seq[2]), kq * abs(seq[3] - seq[2]), kq ** 3 / (1 - kq) * abs(seq[1] - seq[0])]
    assert rep not in leurres
    # point fixe et erreur réelle sous la borne
    assert sp.simplify(expr.subs(vv, ell) - ell) == 0 and a_ <= ell <= b_
    assert abs(sp.N(sp.Rational(seq[3].numerator, seq[3].denominator) - ell)) <= float(rep)

# ex3 (v1, v2) : T(f)(x) = x + lam int_0^1 x t f(t) dt ; f = c x avec c = 3/(3 - lam).
for lv, rep, leurres in [
    (1, F(3, 2), [1, F(2), F(4, 3)]),  # ex3
    (F(1, 2), F(6, 5), [1, F(5, 4), F(8, 7)]),  # ex3 v1
    (-1, F(3, 4), [1, F(3, 2), F(2, 3)]),  # ex3 v2
]:
    L = sp.Rational(lv.numerator, lv.denominator) if isinstance(lv, F) else sp.Integer(lv)
    cc = sp.solve(sp.Eq(c * x, x + L * x * sp.integrate(t * c * t, (t, 0, 1))), c)[0]
    assert cc == sp.Rational(rep.numerator, rep.denominator) == 3 / (3 - L)
    assert abs(L) / 2 < 1  # contraction
    for Lr in leurres:
        assert F(Lr) != rep
    assert sp.integrate(t, (t, 0, 1)) == sp.Rational(1, 2)

# ex4 : T(f)(x) = 1 + int_0^x f ; exp est point fixe ; leurres non points fixes ; T non contractante (f=1, h=0).
for cand, est_pf in [(sp.exp(x), True), (1 + x, False), (sp.exp(x) - 1, False), (sp.exp(-x), False)]:
    assert (sp.simplify(1 + sp.integrate(cand.subs(x, s), (s, 0, x)) - cand) == 0) == est_pf
assert sp.integrate(s ** 3 / sp.factorial(3), (s, 0, x)) == x ** 4 / sp.factorial(4)  # pas de la récurrence
assert sp.integrate(1, (s, 0, 1)) == 1  # ||T(1) - T(0)|| = sup x = 1

# ex5 : int_0^t e^{2Ls} ds = (e^{2Lt} - 1)/(2L), d'où la constante 1/2.
Ls, tt = sp.symbols("L t", positive=True)
assert sp.simplify(Ls * sp.integrate(sp.exp(2 * Ls * s), (s, 0, tt)) - (sp.exp(2 * Ls * tt) - 1) / 2) == 0

# ---------------------------------------------------------------------------------------------------------------
# Carte top-l3-metriques-cauchy-limites
# ---------------------------------------------------------------------------------------------------------------
# affirmation : sqrt(x+1) - sqrt(x) = 1/(sqrt(x+1) + sqrt(x)) -> 0 alors que sqrt x -> +oo.
assert sp.simplify(sp.sqrt(x + 1) - sp.sqrt(x) - 1 / (sp.sqrt(x + 1) + sp.sqrt(x))) == 0
assert sp.limit(sp.sqrt(x + 1) - sp.sqrt(x), x, sp.oo) == 0
# contre-exemple : sin(1/x) vaut 0 en 1/(n pi) et 1 en 1/(2n pi + pi/2).
assert sp.sin(5 * sp.pi) == 0 and sp.sin(2 * 5 * sp.pi + sp.pi / 2) == 1

# ex1 : f1 = x cos(1/x^2) -> 0 ; f3 = cos(ln x) vaut 1 en e^{-2n pi} et -1 en e^{-(2n+1) pi} ; f4 écart 1 + x_n.
assert sp.limit(x * sp.cos(1 / x ** 2), x, 0, "+") == 0
assert sp.cos(sp.log(sp.exp(-2 * 3 * sp.pi))) == 1 and sp.cos(sp.log(sp.exp(-(2 * 3 + 1) * sp.pi))) == -1
xn = 1 / (2 * 3 * sp.pi + sp.pi / 2)
assert sp.simplify((1 + xn) * sp.sin(1 / xn) - (1 + xn)) == 0

# ex2 (v1, v2) : identités d'intégration par parties et majorations, vérifiées numériquement ; seuils A.
cas = [
    (lambda v: mp.sin(v) / v, lambda v: mp.mpf(2) / v, F(1, 100), 200, [100, 50, 400]),  # ex2
    (lambda v: mp.cos(v) / v, lambda v: mp.mpf(2) / v, F(1, 1000), 2000, [1000, 500, 4000]),  # ex2 v1
    (lambda v: mp.sin(v) / mp.sqrt(v), lambda v: mp.mpf(2) / mp.sqrt(v), F(1, 100), 40000, [200, 20000, 10000]),  # ex2 v2
]
for f_, maj, eps_, A_, leurres in cas:
    for a, b in [(1, 7.5), (3, 40), (10, 11), (50, 300)]:
        I = mp.quad(f_, mp.linspace(a, b, 40))
        assert abs(I) <= maj(a) + mp.mpf("1e-20")
    # plus petit A tel que maj(x) <= eps pour tout x >= A (maj décroissante) : maj(A) = eps
    assert abs(maj(A_) - mp.mpf(eps_.numerator) / eps_.denominator) < 1e-25
    for Lr in leurres:
        assert Lr != A_
        assert maj(Lr) > mp.mpf(eps_.numerator) / eps_.denominator or Lr > A_  # trop petit, ou pas le plus petit
# formules exactes des intégrations par parties (dérivée de chaque membre).
assert sp.simplify(sp.diff(-sp.cos(t) / t, t) - sp.cos(t) / t ** 2 - sp.sin(t) / t) == 0
assert sp.simplify(sp.diff(sp.sin(t) / t, t) + sp.sin(t) / t ** 2 - sp.cos(t) / t) == 0
assert sp.simplify(sp.diff(-sp.cos(t) / sp.sqrt(t), t) - sp.cos(t) / (2 * t ** sp.Rational(3, 2)) - sp.sin(t) / sp.sqrt(t)) == 0
assert sp.simplify(sp.Rational(1, 2) * sp.integrate(s ** sp.Rational(-3, 2), (s, x, y)) - (1 / sp.sqrt(x) - 1 / sp.sqrt(y))) == 0

# ex5 (question 4) : aire des pics 1/(2n^2), somme pi^2/12 ; intervalles [n - 1/(2n^2), n + 1/(2n^2)] disjoints.
k = sp.Symbol("k", integer=True, positive=True)
assert sp.summation(1 / (2 * k ** 2), (k, 1, sp.oo)) == sp.pi ** 2 / 12
for m in range(1, 200):
    assert m + F(1, 2 * m * m) < (m + 1) - F(1, 2 * (m + 1) ** 2)
    assert m - F(1, 2 * m * m) >= 0

# exemple : f(x) = int_x^1 cos(1/t) dt est 1-lipschitzienne (|cos| <= 1) ; contrôle numérique d'une limite en 0.
fv = lambda a: mp.quad(lambda v: mp.cos(1 / v), [a, 1])
assert abs(fv(mp.mpf("0.001")) - fv(mp.mpf("0.002"))) <= mp.mpf("0.001")

# QCM q2 : exp(R) = ]0,+oo[, adhérence [0,+oo[ ; q7 : 1/n de Cauchy, limite 0 hors de ]0,1].
assert sp.limit(sp.exp(x), x, -sp.oo) == 0
# QCM q6, q1 : la réciproque de x^3 est la racine cubique, continue.
for v in [-8, -1, 0, 2, 27]:
    assert abs(float(sp.real_root(v ** 3, 3)) - v) < 1e-12

print("top_l3_metriques : toutes les vérifications passent.")
