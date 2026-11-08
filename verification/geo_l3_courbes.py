"""Vérifications calculatoires du chapitre « geo-l3-courbes » (géométrie, L3 : courbes paramétrées) : cartes, exercices
et QCM de content/*/geometrie/l3-geo-l3-courbes.json et content/qcm/geo-l3-courbes.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question. Les calculs
sont symboliques (SymPy), en fractions exactes, ou numériques à 30 chiffres (mpmath) pour les longueurs sans forme close.
Pour chaque réponse finale, on vérifie la bonne proposition et que chaque leurre est faux.
Lancer : python3 verification/geo_l3_courbes.py (doit finir sans erreur).
"""

from fractions import Fraction as F
from math import comb

import mpmath as mp
from sympy import (Abs, Function, Interval, Matrix, Rational, S, atan, binomial, cos, cosh, diff, exp, expand,
                   expand_trig, factor, limit, log, nsolve, oo, pi, simplify, sin, sinh, solveset, sqrt, symbols, tan,
                   trigsimp, asinh, integrate)

mp.mp.dps = 30
t, s, h, u, x, y, z, th = symbols("t s h u x y z theta", real=True)


# ---------------------------------------------------------------------------------------------------------------
# Outils
# ---------------------------------------------------------------------------------------------------------------
def nul(e):
    return simplify(e) == 0


def egal(a, b):
    return simplify(a - b) == 0


def d(f, n=1, var=t):
    return [diff(c, var, n) for c in f]


def det2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def det3(a, b, c):
    return Matrix([list(a), list(b), list(c)]).det()


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def colineaires(a, b):
    """a et b (vecteurs du plan ou de l'espace) sont colinéaires."""
    if len(a) == 2:
        return nul(det2(a, b))
    return all(nul(c) for c in cross(a, b))


def pq(f, t0, var=t, kmax=20):
    """Entiers p et q de l'étude locale en t0 (None si q n'existe pas jusqu'à kmax)."""
    der = [None] + [[simplify(diff(c, var, k).subs(var, t0)) for c in f] for k in range(1, kmax + 1)]
    p = next(k for k in range(1, kmax + 1) if any(c != 0 for c in der[k]))
    for k in range(p + 1, kmax + 1):
        if not colineaires(der[p], der[k]):
            return p, k
    return p, None


NAT = {(1, 0): "ordinaire", (1, 1): "inflexion", (0, 1): "rebroussement 1", (0, 0): "rebroussement 2"}


def nature(p, q):
    return NAT[(p % 2, q % 2)]


def courbure(f, var=t):
    f1, f2 = d(f, 1, var), d(f, 2, var)
    return det2(f1, f2) / (f1[0] ** 2 + f1[1] ** 2) ** Rational(3, 2)


def choix(bonne, leurres, valeur):
    """La bonne proposition vaut `valeur`, chaque leurre en diffère (comparaison exacte ou numérique)."""
    assert simplify(bonne - valeur) == 0, (bonne, valeur)
    for l in leurres:
        assert abs(complex(l - valeur)) > 1e-12, (l, valeur)
    assert len({str(simplify(v)) for v in [bonne, *leurres]}) == 1 + len(leurres)


def signe_change(expr, t0, var=t, eps=Rational(1, 10**6)):
    a, b = expr.subs(var, t0 - eps), expr.subs(var, t0 + eps)
    return float(a) * float(b) < 0


def zeros_numeriques(expr, a, b, var=t, n=200000):
    """Points de [a,b[ où expr change de signe (balayage fin)."""
    import numpy as np
    fn = __import__("sympy").lambdify(var, expr, "numpy")
    xs = np.linspace(float(a), float(b), n, endpoint=False)
    v = fn(xs)
    return [xs[i] for i in range(len(xs) - 1) if v[i] * v[i + 1] < 0]


def bern(k, n, var):
    return binomial(n, k) * var**k * (1 - var) ** (n - k)


def bezier(P, tt):
    n = len(P) - 1
    return tuple(sum(F(comb(n, k)) * tt**k * (1 - tt) ** (n - k) * P[k][i] for k in range(n + 1)) for i in range(len(P[0])))


def casteljau(P, tt):
    niv = [[tuple(F(c) for c in p) for p in P]]
    while len(niv[-1]) > 1:
        q = niv[-1]
        niv.append([tuple((1 - tt) * a + tt * b for a, b in zip(q[i], q[i + 1])) for i in range(len(q) - 1)])
    return niv


# ===============================================================================================================
# Carte geo-l3-courbes-etude-locale
# ===============================================================================================================
# exemple : (t^2, t^3) en 0 : p=2, q=3, rebroussement de première espèce
assert pq((t**2, t**3), 0) == (2, 3) and nature(2, 3) == "rebroussement 1"
# exemple : (t^2, t^4+t^5) : f''=(2,0), f'''=0, f''''=(0,24) ; p=2, q=4 ; demi-branches y = x^2 ± x^(5/2) au-dessus
f = (t**2, t**4 + t**5)
assert d(f, 2) and [c.subs(t, 0) for c in d(f, 2)] == [2, 0] and [c.subs(t, 0) for c in d(f, 3)] == [0, 0]
assert [c.subs(t, 0) for c in d(f, 4)] == [0, 24] and pq(f, 0) == (2, 4) and nature(2, 4) == "rebroussement 2"
X = symbols("X", positive=True)
assert egal((t**4 + t**5).subs(t, sqrt(X)), X**2 + X ** Rational(5, 2))
assert egal((t**4 + t**5).subs(t, -sqrt(X)), X**2 - X ** Rational(5, 2))
assert all((X**2 - X ** Rational(5, 2)).subs(X, v) > 0 for v in [Rational(1, 10), Rational(1, 100), Rational(1, 2)])
# exemple et affirmation (singuliers) : (cos t, sin 2t) en pi/2 : f'=(-1,-2), f''=0, f'''=(1,8), det=-6, p=1, q=3
f = (cos(t), sin(2 * t))
assert [c.subs(t, pi / 2) for c in d(f, 1)] == [-1, -2] and [c.subs(t, pi / 2) for c in d(f, 2)] == [0, 0]
assert [c.subs(t, pi / 2) for c in d(f, 3)] == [1, 8] and det2((-1, -2), (1, 8)) == -6 and pq(f, pi / 2) == (1, 3)
# contre-exemple et pièges : (t+t^2, t^3) ; f'(0)=(1,0), f''(0)=(2,0), f'''(0)=(0,6) ; vrai (p,q)=(1,3), inflexion
f = (t + t**2, t**3)
assert [c.subs(t, 0) for c in d(f, 1)] == [1, 0] and [c.subs(t, 0) for c in d(f, 2)] == [2, 0]
assert [c.subs(t, 0) for c in d(f, 3)] == [0, 6] and pq(f, 0) == (1, 3)
# dans la base (f'(0), f'''(0)) = ((1,0),(0,6)) : X = t+t^2, Y = t^3/6 (change de signe : inflexion)
assert egal(f[1] / 6, t**3 / 6) and signe_change(f[1] / 6, 0)
# piège « q non nul » : premier j>1 avec f^(j)(0) non nul = 2 -> conclurait « allure ordinaire » : faux (Y change de signe)
assert nature(1, 2) == "ordinaire" and nature(*pq(f, 0)) == "inflexion"
# piège « p pas le plus petit » : (f''(0), f'''(0)) libre -> « p=2, q=3 » = rebroussement, or x = t+t^2 change de signe
assert det2((2, 0), (0, 6)) != 0 and nature(2, 3) == "rebroussement 1" and signe_change(t + t**2, 0)
# cours : (t, e^{-1/t^2}) : les dérivées de e^{-1/t^2} (t != 0) tendent vers 0 en 0 (ordres 1 à 4)
for k in range(1, 5):
    assert limit(diff(exp(-1 / t**2), t, k), t, 0) == 0
# affirmation : (t^3, t^5) : p=3, q=5, inflexion
assert pq((t**3, t**5), 0) == (3, 5) and nature(3, 5) == "inflexion"
# QCM q1 et cours : (t^3, t^4) : p=3, q=4, allure ordinaire ; y = |x|^{4/3}
assert pq((t**3, t**4), 0) == (3, 4) and nature(3, 4) == "ordinaire"
assert egal((t**4).subs(t, X ** Rational(1, 3)), X ** Rational(4, 3))
# trou 1 : le leurre f'(t0) est nul dès que p >= 2 (ex. (t^2,t^3)) : il ne dirige rien
assert [c.subs(t, 0) for c in d((t**2, t**3), 1)] == [0, 0]

# démonstration, sur un exemple : coordonnées X(h), Y(h) dans (f^(p), f^(q)) équivalentes à h^p/p!, h^q/q!
f = (t**2 + t**3, 2 * t**2 + 2 * t**3 + t**5)  # p=2 (f''=(2,4)), f'''=(6,12) colinéaire, q=5
assert pq(f, 0) == (2, 5)
uvec, vvec = [c.subs(t, 0) for c in d(f, 2)], [c.subs(t, 0) for c in d(f, 5)]
M = Matrix([uvec, vvec]).T
XY = M.solve(Matrix(f))
assert limit(XY[0] / (t**2 / 2), t, 0) == 1 and limit(XY[1] / (t**5 / 120), t, 0) == 1

# ===============================================================================================================
# Exercices geo-l3-courbes-etude-locale-ex1 (et v1, v2)
# ===============================================================================================================
f = (t**2 - 2 * t, t**3 - 3 * t)
assert [c.subs(t, 1) for c in d(f, 1)] == [0, 0] and [c.subs(t, 1) for c in d(f, 2)] == [2, 6]
assert [c.subs(t, 1) for c in d(f, 3)] == [0, 6] and pq(f, 1) == (2, 3) and nature(2, 3) == "rebroussement 1"
assert [c.subs(t, 1) for c in f] == [-1, -2] and egal(3 * (-1) + 1, -2)  # tangente y = 3x+1 passe par (-1,-2), pente 6/2
f = (t**2, t**4 + t**5)  # v1
assert pq(f, 0) == (2, 4) and nature(2, 4) == "rebroussement 2" and det2((2, 0), (0, 24)) == 48
f = (t + t**2, t**3 + t**4)  # v2
assert [c.subs(t, 0) for c in d(f, 1)] == [1, 0] and [c.subs(t, 0) for c in d(f, 2)] == [2, 0]
assert [c.subs(t, 0) for c in d(f, 3)] == [0, 6] and pq(f, 0) == (1, 3) and nature(1, 3) == "inflexion"
# les quatre natures proposées sont distinctes : chaque leurre est une autre nature, donc faux
assert len(set(NAT.values())) == 4

# ===============================================================================================================
# geo-l3-courbes-etude-locale-ex2 (et v1, v2) : points d'inflexion de courbes de Lissajous
# ===============================================================================================================
c_, s_ = cos(t), sin(t)
for f, D_attendu, n_attendu in [
    ((cos(t), sin(2 * t)), 2 * cos(t) * (1 + 2 * sin(t) ** 2), 2),  # ex2
    ((cos(t), sin(3 * t)), -3 * (8 * cos(t) ** 4 - 12 * cos(t) ** 2 + 3), 4),  # v1
    ((sin(2 * t), sin(3 * t)), -6 * sin(t) * (8 * sin(t) ** 4 - 10 * sin(t) ** 2 + 5), 2),  # v2
]:
    D = det2(d(f), d(f, 2))
    assert nul(expand_trig(D - D_attendu))
    vit = d(f)[0] ** 2 + d(f)[1] ** 2
    import numpy as np
    xs = np.linspace(0, 2 * np.pi, 20001)
    assert min(__import__("sympy").lambdify(t, vit, "numpy")(xs)) > 0.2  # régulière
    zs = zeros_numeriques(D, -1e-3, 2 * np.pi - 1e-3)  # inclut t = 0 pour v2
    assert len(zs) == n_attendu, (f, zs)
# ex2 : zéros en pi/2 et 3pi/2, tous deux en l'origine
assert [c.subs(t, pi / 2) for c in (cos(t), sin(2 * t))] == [0, 0] == [c.subs(t, 3 * pi / 2) for c in (cos(t), sin(2 * t))]
# v1 : racines du trinôme 8X^2-12X+3, une seule dans [0,1]
r1, r2 = (3 - sqrt(3)) / 4, (3 + sqrt(3)) / 4
assert nul(8 * r1**2 - 12 * r1 + 3) and nul(8 * r2**2 - 12 * r2 + 3) and 0 < r1 < 1 < r2
# v1 : 9 sin t sin 3t + 3 cos t cos 3t = 27s^2-36s^4+12c^4-9c^2
assert nul(expand_trig(9 * s_ * sin(3 * t) + 3 * c_ * cos(3 * t) - (27 * s_**2 - 36 * s_**4 + 12 * c_**4 - 9 * c_**2)))
# v2 : linéarisations, sin 5t, discriminant négatif, f(0)=f(pi)=(0,0)
assert nul(expand_trig(12 * sin(2 * t) * cos(3 * t) - 6 * (sin(5 * t) - s_)))
assert nul(expand_trig(18 * cos(2 * t) * sin(3 * t) - 9 * (sin(5 * t) + s_)))
assert nul(expand_trig(sin(5 * t) - (16 * s_**5 - 20 * s_**3 + 5 * s_))) and 100 - 160 < 0
assert [c.subs(t, pi) for c in (sin(2 * t), sin(3 * t))] == [0, 0]
# réponses : ex2 -> 2 (leurres 0, 4, 1) ; v1 -> 4 (leurres 2, 6, 0) ; v2 -> 2 (leurres 4, 6, 0)
choix(2, [0, 4, 1], 2)
choix(4, [2, 6, 0], 4)
choix(2, [4, 6, 0], 2)

# ===============================================================================================================
# geo-l3-courbes-etude-locale-ex3 (et v1, v2) : asymptotes
# ===============================================================================================================
for k_, b_att in [(3, -1), (6, -2)]:  # ex3 (folium) et v1
    xt, yt = k_ * t / (1 + t**3), k_ * t**2 / (1 + t**3)
    assert limit(yt / xt, t, -1) == -1
    b = limit(yt + xt, t, -1)
    assert b == b_att
    assert nul(xt + yt - b - (k_ // 3) * (t + 1) ** 2 / (t**2 - t + 1))  # position : > 0 de part et d'autre
    assert limit(Abs(xt), t, -1, "+") == oo
# ex3 : x+y+1=0 ; leurres x+y-1, x+y+3, x+y : constantes 1, -1, 3, 0 -> seule la constante 1 convient
choix(1, [-1, 3, 0], -(-1))
choix(2, [-2, 6, 1], -(-2))  # v1 : x+y+2=0
# v2 : x = t^2/(t-1), y = t/(t-1), t -> 1 : a = 1, b = -1, position y-(x-1) = 1-t
xt, yt = t**2 / (t - 1), t / (t - 1)
assert limit(yt / xt, t, 1) == 1 and limit(yt - xt, t, 1) == -1 and nul(yt - (xt - 1) - (1 - t))
choix(-1, [1, 0, -2], -1)  # ordonnées à l'origine de y = x-1 (leurres y=x+1, y=x, y=x-2)

# ===============================================================================================================
# geo-l3-courbes-etude-locale-ex4 : graphe local, phi'' = det(f',f'')/x'^3
# ===============================================================================================================
xf, yf = Function("x")(t), Function("y")(t)
phi1 = diff(yf, t) / diff(xf, t)  # phi'(x(t))
phi2 = diff(phi1, t) / diff(xf, t)  # phi''(x(t))
Dg = diff(xf, t) * diff(yf, t, 2) - diff(xf, t, 2) * diff(yf, t)
assert nul(phi2 - Dg / diff(xf, t) ** 3)
# leurres faux sur un exemple : f(t) = (t + t^2, t^3) en t = 1
fe = (t + t**2, t**3)
D1 = det2(d(fe), d(fe, 2)).subs(t, 1)
xp, xs2, ys2 = d(fe)[0].subs(t, 1), d(fe, 2)[0].subs(t, 1), d(fe, 2)[1].subs(t, 1)
bonne = D1 / xp**3
choix(bonne, [D1 / xp**2, D1 / (xp**2 + d(fe)[1].subs(t, 1) ** 2) ** Rational(3, 2), ys2 / xs2], phi2.subs({xf: fe[0], yf: fe[1]}).doit().subs(t, 1))

# ===============================================================================================================
# Carte geo-l3-courbes-singuliers
# ===============================================================================================================
f = (t - sin(t), 1 - cos(t))  # cycloïde : f'=0 ssi t dans 2piZ ; en 0, f''=(0,1), f'''=(1,0), p=2, q=3
assert solveset(d(f)[0] ** 2 + d(f)[1] ** 2, t, Interval.Ropen(0, 2 * pi)) == {0}
assert [c.subs(t, 0) for c in d(f, 2)] == [0, 1] and [c.subs(t, 0) for c in d(f, 3)] == [1, 0] and pq(f, 0) == (2, 3)
# exemple : tan admissible, (tan s, tan^2 s) régulier
assert nul(diff(tan(s), s) - (1 + tan(s) ** 2))
# contre-exemple : (s^3, s^6) : p=3, q=6 ; (t,t^2) : p=1, q=2
assert pq((s**3, s**6), 0, var=s) == (3, 6) and pq((t, t**2), 0) == (1, 2)
# contre-exemple : theta(s)=s/2 sur [0,2pi] -> demi-cercle supérieur (sin >= 0)
assert all(sin(v / 2) >= 0 for v in [0, 1, 2, 3, 4, 5, 6])
# affirmation : t^3, t^4 allure ordinaire ; démonstration (lemme et entiers p, q) sur un exemple :
# f = (t^2, t^3), theta(s) = 2s + s^2 (theta(0)=0, theta'(0)=2) : g = f o theta a les mêmes (p,q) et g''(0) = 2^2 f''(0)
g = (((2 * s + s**2) ** 2), (2 * s + s**2) ** 3)
assert pq(g, 0, var=s) == (2, 3) == pq((t**2, t**3), 0)
assert [diff(c, s, 2).subs(s, 0) for c in g] == [4 * 2, 0]
# lemme : g^(j) = theta'^j f^(j)(theta) + combinaison des f^(m)(theta), m < j. Avec f = exp(a t), f^(m)(theta) =
# a^m e^(a theta) : g^(j) e^(-a theta) est un polynôme en a de degré j, sans terme constant, de coefficient dominant theta'^j
aa_, thf = symbols("alpha"), Function("theta")(s)
for j in range(1, 7):
    P_ = __import__("sympy").Poly(simplify(diff(exp(aa_ * thf), s, j) * exp(-aa_ * thf)), aa_)
    assert P_.degree() == j and P_.coeff_monomial(1) == 0 and nul(P_.coeff_monomial(aa_**j) - diff(thf, s) ** j)
# trou 1 et pièges : s^3 strictement croissante mais g'(0)=0 alors que f'(0) != 0 ; theta(s)=s^2 sur [0,1] (theta''=2)
assert diff(s**3, s).subs(s, 0) == 0 and diff(s**2, s, 2) == 2 and diff(s**2, s).subs(s, 0) == 0
# trou 2 : f(t) = (t, 0), theta(s) = 2s, en s = 1/3 : g'(s) = (2, 0) ; leurres f'(theta(s)) = (1, 0),
# theta(s) f'(theta(s)) = (2/3, 0), f'(theta'(s)) = (1, 0) : tous faux
fl, thl = (t, 0 * t), 2 * s
gl = [diff(c.subs(t, thl), s).subs(s, Rational(1, 3)) for c in fl]
fprime = [diff(c, t) for c in fl]
assert gl == [2, 0]
assert [c.subs(t, thl).subs(s, Rational(1, 3)) for c in fprime] != gl
assert [(thl * c.subs(t, thl)).subs(s, Rational(1, 3)) for c in fprime] != gl
assert [c.subs(t, diff(thl, s)) for c in fprime] != gl

# ===============================================================================================================
# geo-l3-courbes-singuliers-ex1 (et v1, v2) : points singuliers
# ===============================================================================================================
for f, n_att, t0, f2, f3 in [
    ((cos(t) ** 3, sin(t) ** 3), 4, 0, [-3, 0], [0, 6]),  # astroïde
    ((2 * cos(t) - cos(2 * t), 2 * sin(t) - sin(2 * t)), 1, 0, [2, 0], [0, 6]),  # v1 cardioïde
    ((2 * cos(t) + cos(2 * t), 2 * sin(t) - sin(2 * t)), 3, 0, [-6, 0], [0, 6]),  # v2 deltoïde
]:
    vit = trigsimp(d(f)[0] ** 2 + d(f)[1] ** 2)
    zs = solveset(vit, t, Interval.Ropen(0, 2 * pi))
    assert len(zs) == n_att, (f, zs)
    assert [c.subs(t, t0) for c in d(f, 2)] == f2 and [c.subs(t, t0) for c in d(f, 3)] == f3 and pq(f, t0) == (2, 3)
assert nul(expand_trig(d((2 * cos(t) - cos(2 * t), 0))[0] ** 2 + d((0, 2 * sin(t) - sin(2 * t)))[1] ** 2 - 8 * (1 - cos(t))))
assert nul(expand_trig((-2 * sin(t) - 2 * sin(2 * t)) ** 2 + (2 * cos(t) - 2 * cos(2 * t)) ** 2 - 8 * (1 - cos(3 * t))))
# DL : 2cos h - cos 2h = 1 + h^2 + O(h^4) ; 2 sin h - sin 2h = h^3 + O(h^5) ; 2cos h + cos 2h = 3 - 3h^2 + O(h^4)
assert (2 * cos(h) - cos(2 * h)).series(h, 0, 4).removeO() == 1 + h**2
assert (2 * sin(h) - sin(2 * h)).series(h, 0, 5).removeO() == h**3
assert (2 * cos(h) + cos(2 * h)).series(h, 0, 4).removeO() == 3 - 3 * h**2
assert (cos(h) ** 3).series(h, 0, 4).removeO() == 1 - 3 * h**2 / 2 and (sin(h) ** 3).series(h, 0, 4).removeO() == h**3
# symétries : astroïde f(t+pi/2) = rotation de f(t) ; deltoïde z(t+2pi/3) = e^{2i pi/3} z(t)
assert nul(expand_trig(cos(t + pi / 2) ** 3 + sin(t) ** 3)) and nul(expand_trig(sin(t + pi / 2) ** 3 - cos(t) ** 3))
from sympy import I as Ii
zt = 2 * exp(Ii * t) + exp(-2 * Ii * t)
assert nul(expand(zt.subs(t, t + 2 * pi / 3) - exp(2 * Ii * pi / 3) * zt))
choix(4, [1, 2, 3], 4)
choix(1, [2, 3, 4], 1)
choix(3, [1, 2, 4], 3)

# ===============================================================================================================
# geo-l3-courbes-singuliers-ex2 (et v1, v2) : changements de paramètre
# ===============================================================================================================
assert (3 * s**2 + 1).subs(s, 0) > 0 and diff(s**3, s).subs(s, 0) == 0
assert limit(atan(s), s, oo) == pi / 2 and diff(s + sin(s), s).subs(s, pi) == 0  # ex2
assert diff(log(s) ** 3, s).subs(s, 1) == 0 and diff(log(s), s) == 1 / s  # v1
assert diff(log(s) + sin(log(s)), s).subs(s, exp(pi)) == 0
assert simplify(diff(tan(s) ** 3, s).subs(s, 0)) == 0 and simplify(diff(tan(s) - s, s) - tan(s) ** 2) == 0  # v2
assert limit(tan(s) - s, s, pi / 2, "-") == oo

# ===============================================================================================================
# geo-l3-courbes-singuliers-ex3 (et v1, v2) : (p,q) de f et du graphe
# ===============================================================================================================
for f, pq_f, h_, pq_g, nat in [
    ((t**3, t**3 + t**6), (3, 6), x + x**2, (1, 2), "ordinaire"),
    ((t**5, t**5 + t**10), (5, 10), x + x**2, (1, 2), "ordinaire"),  # v1
    ((t**3, t**3 - t**9), (3, 9), x - x**3, (1, 3), "inflexion"),  # v2
]:
    assert pq(f, 0) == pq_f and nature(*pq_f) == nat
    assert pq((x, h_), 0, var=x) == pq_g and nature(*pq_g) == nat
    assert nul(f[1] - h_.subs(x, f[0]))  # le support est dans le graphe
# leurres : ex3 (1,2), (3,4), (3,5) ; v1 (1,2), (5,6), (5,15) ; v2 (1,3), (3,6), (3,4) : tous différents de (p,q)
for vrai, leurres in [((3, 6), [(1, 2), (3, 4), (3, 5)]), ((5, 10), [(1, 2), (5, 6), (5, 15)]), ((3, 9), [(1, 3), (3, 6), (3, 4)])]:
    assert vrai not in leurres and len(set(leurres)) == 3

# ===============================================================================================================
# geo-l3-courbes-singuliers-ex4 : la pointe de y^2 = x^3
# ===============================================================================================================
for v in [Rational(-3, 2), Rational(-1, 3), 0, Rational(2, 7), 2]:
    X_, Y_ = v**2, v**3
    assert X_ >= 0 and Y_**2 == X_**3
for X_, Y_ in [(Rational(4, 9), Rational(8, 27)), (Rational(4, 9), -Rational(8, 27)), (0, 0)]:
    tt = sqrt(X_) if Y_ >= 0 else -sqrt(X_)
    assert tt**2 == X_ and tt**3 == Y_
# un paramétrage C^1 passant par l'origine a une vitesse nulle : exemple g(s) = (s^2, s^3) et g(s) = (s^4, s^6)
assert [diff(c, s).subs(s, 0) for c in (s**4, s**6)] == [0, 0]

# ===============================================================================================================
# Carte geo-l3-courbes-polaires
# ===============================================================================================================
rho = Function("rho")(th)
uv, vv = (cos(th), sin(th)), (-sin(th), cos(th))
fp = [rho * c for c in uv]
f1, f2 = d(fp, 1, th), d(fp, 2, th)
assert all(nul(f1[i] - (diff(rho, th) * uv[i] + rho * vv[i])) for i in range(2))
assert nul(f1[0] ** 2 + f1[1] ** 2 - (rho**2 + diff(rho, th) ** 2))
assert all(nul(f2[i] - ((diff(rho, th, 2) - rho) * uv[i] + 2 * diff(rho, th) * vv[i])) for i in range(2))
assert nul(det2(f1, f2) - (rho**2 + 2 * diff(rho, th) ** 2 - rho * diff(rho, th, 2)))
# u = 1/rho : rho^2 + 2 rho'^2 - rho rho'' = rho^3 (u + u'') (cours, polaires-ex4)
uf = Function("u")(th)
r_ = 1 / uf
assert nul(r_**2 + 2 * diff(r_, th) ** 2 - r_ * diff(r_, th, 2) - r_**3 * (uf + diff(uf, th, 2)))
# cours : distance algébrique de rho u(theta) à la droite d'angle theta0 = rho sin(theta - theta0)
th0 = symbols("theta0", real=True)
assert nul(expand_trig(rho * cos(th) * (-sin(th0)) + rho * sin(th) * cos(th0) - rho * sin(th - th0)))
# exemple cardioïde : rho(pi) = rho'(pi) = 0 ; tan V = -1/tan(theta/2) ; V = theta/2 + pi/2 (vecteur f' colinéaire)
r = 1 + cos(th)
assert r.subs(th, pi) == 0 and diff(r, th).subs(th, pi) == 0
assert all(abs(float((r / diff(r, th) + 1 / tan(th / 2)).subs(th, v))) < 1e-12 for v in [0.3, 1.1, -2.0, 2.9])
assert nul((r / diff(r, th) + 1 / tan(th / 2)).rewrite(tan))
assert nul(trigsimp(det2((diff(r, th), r), (cos(th / 2 + pi / 2), sin(th / 2 + pi / 2)))))
# spirale : tan V = 1/a ; affirmation e^theta : V = pi/4
a_ = symbols("a", nonzero=True, real=True)
assert nul(exp(a_ * th) / diff(exp(a_ * th), th) - 1 / a_)
# contre-exemple rho = 1 + |theta| : demi-tangentes (1,1) et (-1,1) en M(0) = (1,0)
assert [1 * 1 + 1 * 0, 1 * 0 + 1 * 1] == [1, 1] and [-1 * 1 + 1 * 0, -1 * 0 + 1 * 1] == [-1, 1]
# affirmations : rho = cos theta en pi/2 : rho = 0, rho' = -1 ; M(theta + pi) = M(theta) ; rho = 1/cos : x = 1
assert cos(pi / 2) == 0 and diff(cos(th), th).subs(th, pi / 2) == -1
assert all(nul(cos(th + pi) * c.subs(th, th + pi) - cos(th) * c) for c in uv)
assert nul(1 / cos(th) * cos(th) - 1)
# trou 2 : leurres faux, par exemple pour rho = e^{2 theta} (rho/rho' = 1/2)
r = exp(2 * th)
choix(r / diff(r, th), [diff(r, th) / r, -r / diff(r, th), (r * diff(r, th)).subs(th, 1)], Rational(1, 2))

# ===============================================================================================================
# geo-l3-courbes-polaires-ex1 (et v1, v2) : tan V
# ===============================================================================================================
for r, t0, val, leurres in [
    (1 + cos(th), pi / 3, -sqrt(3), [-1 / sqrt(3), sqrt(3), -sqrt(3) / 2]),
    (2 + cos(th), pi / 2, -2, [-Rational(1, 2), 2, Rational(1, 2)]),  # v1
    (1 + 2 * cos(th), pi / 2, -Rational(1, 2), [-2, Rational(1, 2), 2]),  # v2
]:
    rv, rpv = r.subs(th, t0), diff(r, th).subs(th, t0)
    assert rv != 0 and rpv < 0 < rv  # cos V < 0, sin V > 0
    choix(val, leurres, rv / rpv)
assert abs(float(pi - atan(2)) * 180 / float(pi) - 116.565) < 1e-3
assert abs(float(pi - atan(Rational(1, 2))) * 180 / float(pi) - 153.435) < 1e-3

# ===============================================================================================================
# geo-l3-courbes-polaires-ex2 (et v1, v2) : tangentes horizontales
# ===============================================================================================================
C = symbols("C", real=True)
for r, poly, cval, xval, leurres in [
    (1 + cos(th), 2 * C**2 + C - 1, Rational(1, 2), Rational(3, 4), [Rational(3, 2), 3 * sqrt(3) / 4, Rational(1, 2)]),
    (2 + cos(th), 2 * C**2 + 2 * C - 1, (sqrt(3) - 1) / 2, sqrt(3) / 2, [(sqrt(3) - 1) / 2, (sqrt(3) + 3) / 2, Rational(1, 2)]),
    (3 + cos(th), 2 * C**2 + 3 * C - 1, (sqrt(17) - 3) / 4, (3 * sqrt(17) - 5) / 8,
     [(sqrt(17) - 3) / 4, (sqrt(17) + 9) / 4, (3 * sqrt(17) + 5) / 8]),
]:
    yp = diff(r * sin(th), th)
    assert nul(expand_trig(yp - poly.subs(C, cos(th)).subs(sin(th) ** 2, 1 - cos(th) ** 2)))
    assert nul(poly.subs(C, cval)) and -1 <= cval <= 1
    autres = [rt for rt in __import__("sympy").solve(poly, C) if simplify(rt - cval) != 0]
    assert all(abs(float(rt)) > 1 or simplify(rt + 1) == 0 for rt in autres)
    choix(xval, leurres, simplify((r.subs(cos(th), cval) * cval)))
# ex2 : points (3/4, ±3 sqrt3/4) ; au pôle (theta = pi) le point est singulier
assert nul((1 + cos(pi / 3)) * sin(pi / 3) - 3 * sqrt(3) / 4)

# ===============================================================================================================
# geo-l3-courbes-polaires-ex3 (et v1, v2) : angle V constant
# ===============================================================================================================
for r, V, leurres in [
    (exp(th / sqrt(3)), pi / 3, [pi / 6, 2 * pi / 3, pi / 4]),
    (exp(th), pi / 4, [3 * pi / 4, pi / 3, pi / 6]),  # v1
    (exp(-th), 3 * pi / 4, [pi / 4, 2 * pi / 3, pi / 2]),  # v2
]:
    rp = diff(r, th)
    for th_v in (0, 1):  # V ne dépend pas de theta
        rv, rpv = r.subs(th, th_v), rp.subs(th, th_v)
        choix(V, leurres, __import__("sympy").atan2(simplify(rv / sqrt(rv**2 + rpv**2)), simplify(rpv / sqrt(rv**2 + rpv**2))))
# question 2 : solutions de rho' = rho / tan V0
V0 = symbols("V0", positive=True)
r = exp(th / tan(V0))
assert nul(diff(r, th) - r / tan(V0))

# ===============================================================================================================
# geo-l3-courbes-polaires-ex4 : inflexions de rho = 1/(2 + cos 2 theta)
# ===============================================================================================================
uu = 2 + cos(2 * th)
assert nul(uu + diff(uu, th, 2) - (2 - 3 * cos(2 * th)))
r = 1 / uu
assert nul(simplify(r**2 + 2 * diff(r, th) ** 2 - r * diff(r, th, 2) - r**3 * (2 - 3 * cos(2 * th))))
zs = solveset(2 - 3 * cos(2 * th), th, Interval.Ropen(0, 2 * pi))
assert len(zs) == 4 and all(abs(float(sin(2 * z_))) > 0.5 for z_ in zs)
choix(4, [0, 2, 8], len(zs))
# question 2 : u = a cos + b sin vérifie u + u'' = 0
A_, B_ = symbols("A B", real=True)
uu = A_ * cos(th) + B_ * sin(th)
assert nul(uu + diff(uu, th, 2))

# ===============================================================================================================
# Carte geo-l3-courbes-espace
# ===============================================================================================================
a, b = symbols("a b", positive=True)
H = (a * cos(t), a * sin(t), b * t)
cr = cross(d(H), d(H, 2))
assert all(nul(cr[i] - [a * b * sin(t), -a * b * cos(t), a**2][i]) for i in range(3))
n2 = simplify(sum(c**2 for c in cr))
assert nul(n2 - a**2 * (a**2 + b**2))
kap = sqrt(n2) / sqrt(sum(c**2 for c in d(H))) ** 3
tau = det3(d(H), d(H, 2), d(H, 3)) / n2
assert nul(kap - a / (a**2 + b**2)) and nul(tau - b / (a**2 + b**2))
# torsion invariante par changement de paramètre (exemple : hélice, theta(s) = s^3 + s)
Hs = [c.subs(t, s**3 + s) for c in (cos(t), sin(t), t)]
crs = cross(d(Hs, 1, s), d(Hs, 2, s))
tau_s = det3(d(Hs, 1, s), d(Hs, 2, s), d(Hs, 3, s)) / sum(c**2 for c in crs)
assert nul(tau_s - Rational(1, 2))
# exemple hélice (cos t, sin t, t) : f' ^ f'' = (sin t, -cos t, 1) ; en 0 normal (0,-1,1) ; det = 1
H1 = (cos(t), sin(t), t)
assert [simplify(c) for c in cross(d(H1), d(H1, 2))] == [sin(t), -cos(t), 1]
assert [c.subs(t, 0) for c in cross(d(H1), d(H1, 2))] == [0, -1, 1] and nul(det3(d(H1), d(H1, 2), d(H1, 3)) - 1)
# exemple cubique : f'^f'' = (6t^2, -6t, 2), det = 12, plan osculateur 3t^2 x - 3t y + z = t^3
Cu = (t, t**2, t**3)
assert [expand(c) for c in cross(d(Cu), d(Cu, 2))] == [6 * t**2, -6 * t, 2] and det3(d(Cu), d(Cu, 2), d(Cu, 3)) == 12
assert nul(3 * t**2 * t - 3 * t * t**2 + t**3 - t**3)
# trou 1 : sur l'hélice en 0, f''' n'est pas dans le plan osculateur ; (f'',f''') engendre un autre plan
n0 = [0, -1, 1]
assert sum(n0[i] * [c.subs(t, 0) for c in d(H1, 3)][i] for i in range(3)) != 0
assert not colineaires(cross([c.subs(t, 0) for c in d(H1, 2)], [c.subs(t, 0) for c in d(H1, 3)]), n0)
# trou 2 : cercle (cos t, sin t, 1) plan mais det(f, f', f'') = 1 ; hélice : ||f'^f''|| constante, non plane
Ci = (cos(t), sin(t), 1)
assert nul(det3(Ci, d(Ci), d(Ci, 2)) - 1) and nul(det3(d(Ci), d(Ci, 2), d(Ci, 3)))
assert nul(sum(c**2 for c in cross(d(H1), d(H1, 2))) - 2)
# démonstration : changement de paramètre, déterminant de (g', g'') dans (f', f'') = theta'^3
thp, thpp = symbols("thp thpp")
assert Matrix([[thp, 0], [thpp, thp**2]]).det() == thp**3
# contre-exemple par morceaux : C^3, régulière, det nul, quatre points non coplanaires
gauche, droite = (t, t**4, 0), (t, 0, t**4)
for k in range(0, 4):
    assert [diff(c, t, k).subs(t, 0) for c in gauche] == [diff(c, t, k).subs(t, 0) for c in droite]
assert [diff(c, t, 4).subs(t, 0) for c in gauche] != [diff(c, t, 4).subs(t, 0) for c in droite]
assert det3(d(gauche), d(gauche, 2), d(gauche, 3)) == 0 and det3(d(droite), d(droite, 2), d(droite, 3)) == 0
assert det3((-1, 1, 0), (-2, 16, 0), (1, 0, 1)) == -14
# contre-exemple : deux arcs de cercle dans z=0 et z=1 : chacun de déterminant nul, et 4 points non coplanaires
assert det3((cos(-1.5) - cos(-1.2), sin(-1.5) - sin(-1.2), 0), (cos(-1.8) - cos(-1.2), sin(-1.8) - sin(-1.2), 0),
            (cos(1.5) - cos(-1.2), sin(1.5) - sin(-1.2), 1)) != 0

# ===============================================================================================================
# geo-l3-courbes-espace-ex1 (et v1, v2) : plans osculateurs
# ===============================================================================================================
def plan_osc(f, t0):
    P = [c.subs(t, t0) for c in f]
    n = [c.subs(t, t0) for c in cross(d(f), d(f, 2))]
    return P, n


def est_plan(P, n, coeffs):
    """coeffs = (A, B, C, D) pour Ax+By+Cz=D : plan passant par P, de normal colinéaire à n."""
    A, B, Cc, D = coeffs
    return colineaires([A, B, Cc], n) and A * P[0] + B * P[1] + Cc * P[2] == D


for f, t0, bonne, leurres in [
    ((cos(t), sin(t), t), 0, (0, 1, -1, 0), [(0, 1, 1, 0), (1, 0, 0, 1), (1, 1, -1, 1)]),
    ((2 * cos(t), 2 * sin(t), t), 0, (0, 1, -2, 0), [(0, 2, -1, 0), (0, 2, 1, 0), (1, 0, 0, 2)]),  # v1
    ((t, t**2, t**3), 1, (3, -3, 1, 1), [(1, 2, 3, 6), (3, 3, 1, 7), (3, -3, 1, 0)]),  # v2
]:
    P, n = plan_osc(f, t0)
    assert any(c != 0 for c in n) and est_plan(P, n, bonne)
    assert all(not est_plan(P, n, l) for l in leurres)
# v1 / ex1 : les leurres « plan normal » sont orthogonaux à la tangente
assert est_plan([1, 0, 0], [0, 1, 1], (0, 1, 1, 0)) and est_plan([2, 0, 0], [0, 2, 1], (0, 2, 1, 0))

# ===============================================================================================================
# geo-l3-courbes-espace-ex2 (et v1, v2) : torsion et courbure d'hélices
# ===============================================================================================================
for (av, bv), tau_v, kap_v, leurres in [
    ((1, 1), Rational(1, 2), Rational(1, 2), [1, 1 / sqrt(2), 2]),
    ((3, 4), Rational(4, 25), Rational(3, 25), [Rational(3, 25), Rational(4, 5), Rational(25, 4)]),  # v1
    ((2, 1), Rational(1, 5), Rational(2, 5), [Rational(2, 5), Rational(1, 2), 5]),  # v2
]:
    Hv = (av * cos(t), av * sin(t), bv * t)
    crv = cross(d(Hv), d(Hv, 2))
    tv = simplify(det3(d(Hv), d(Hv, 2), d(Hv, 3)) / sum(c**2 for c in crv))
    kv = simplify(sqrt(sum(c**2 for c in crv)) / sqrt(sum(c**2 for c in d(Hv))) ** 3)
    assert kv == kap_v
    choix(tau_v, leurres, tv)
assert simplify(sum(c**2 for c in cross(d((3 * cos(t), 3 * sin(t), 4 * t)), d((3 * cos(t), 3 * sin(t), 4 * t), 2)))) == 225

# ===============================================================================================================
# geo-l3-courbes-espace-ex3 (et v1, v2) : plans osculateurs de la cubique passant par A
# ===============================================================================================================
for A, n_att, leurres in [((0, -1, 0), 3, [1, 2, 0]), ((0, 1, 0), 1, [3, 2, 0]), ((1, 0, 0), 2, [3, 1, 0])]:
    poly = t**3 - 3 * A[0] * t**2 + 3 * A[1] * t - A[2]
    # A est sur le plan osculateur en M(t) ssi 3t^2 x - 3t y + z = t^3 en A
    assert nul(poly + (3 * t**2 * A[0] - 3 * t * A[1] + A[2] - t**3))
    racines = solveset(poly, t, S.Reals)
    choix(n_att, leurres, len(racines))

# ===============================================================================================================
# geo-l3-courbes-espace-ex4 : la courbe en deux morceaux n'est pas birégulière en 0 seulement
# ===============================================================================================================
for tv in [-2, -1, Rational(-1, 3), Rational(1, 3), 1, 2]:
    morceau = gauche if tv < 0 else droite
    assert any(c != 0 for c in cross([diff(c, t).subs(t, tv) for c in morceau], [diff(c, t, 2).subs(t, tv) for c in morceau]))
assert [diff(c, t, 2).subs(t, 0) for c in gauche] == [0, 0, 0]

# ===============================================================================================================
# geo-l3-courbes-espace-ex5 : dérivée de det(g, f', f'') = det(g, f', f''') (exemple de l'hélice, A quelconque)
# ===============================================================================================================
A3 = symbols("A1:4", real=True)
gvec = [H1[i] - A3[i] for i in range(3)]
assert nul(diff(det3(gvec, d(H1), d(H1, 2)), t) - det3(gvec, d(H1), d(H1, 3)))

# ===============================================================================================================
# Carte geo-l3-courbes-longueur
# ===============================================================================================================
cyc = (t - sin(t), 1 - cos(t))
assert nul(expand_trig(d(cyc)[0] ** 2 + d(cyc)[1] ** 2 - 4 * sin(t / 2) ** 2))
assert integrate(2 * sin(t / 2), (t, 0, 2 * pi)) == 8
R_ = symbols("R", positive=True)
gc = (R_ * cos(s / R_), R_ * sin(s / R_))
assert nul(d(gc, 1, s)[0] ** 2 + d(gc, 1, s)[1] ** 2 - 1)
# contre-exemple (t^2, t^3) : s(t) ~ t|t| (t > 0 : s(t)/t^2 -> 1), et g_1 ~ |sigma|
sfun = integrate(u * sqrt(4 + 9 * u**2), (u, 0, t))
assert limit(sfun / t**2, t, 0, "+") == 1 and limit(t**2 / sfun, t, 0, "+") == 1
# contre-exemple (t,|t|) : ||f'|| = sqrt 2 ; cercle sur [0,4pi] : 4 pi
assert integrate(1, (t, 0, 4 * pi)) == 4 * pi
# affirmation : spirale e^{-theta}, theta >= 0 : longueur sqrt 2
r = exp(-th)
assert nul(sqrt(r**2 + diff(r, th) ** 2) - sqrt(2) * exp(-th)) and integrate(sqrt(2) * exp(-th), (th, 0, oo)) == sqrt(2)
# trou 2 : (t - t0) ||f'(t)|| n'est pas l'abscisse curviligne pour f(t) = (t^2, 0), t0 = 1
assert integrate(2 * u, (u, 1, t)) != (t - 1) * 2 * t

# ===============================================================================================================
# geo-l3-courbes-longueur-ex1 (et v1, v2)
# ===============================================================================================================
choix(8, [4, 2 * pi, 16], integrate(2 * sin(t / 2), (t, 0, 2 * pi)))
ast = (cos(t) ** 3, sin(t) ** 3)
assert nul(trigsimp(d(ast)[0] ** 2 + d(ast)[1] ** 2 - 9 * sin(t) ** 2 * cos(t) ** 2))
L_ast = 4 * integrate(3 * sin(t) * cos(t), (t, 0, pi / 2))
choix(6, [Rational(3, 2), 0, 3], L_ast)
assert integrate(3 * sin(t) * cos(t), (t, 0, 2 * pi)) == 0
L_ch = integrate(cosh(t), (t, 0, log(2)))
choix(Rational(3, 4), [Rational(5, 4), log(2), Rational(3, 2)], simplify(L_ch.rewrite(exp)))
assert simplify(cosh(log(2)).rewrite(exp)) == Rational(5, 4)

# ===============================================================================================================
# geo-l3-courbes-longueur-ex2 (et v1, v2) : spirales logarithmiques
# ===============================================================================================================
aa = symbols("aa", negative=True)
for k_, Lfin, leurres, raison in [
    (1, sqrt(2), [1, 2], exp(-2 * pi)),
    (2, sqrt(5) / 2, [Rational(1, 2), sqrt(5)], exp(-4 * pi)),  # v1
    (Rational(1, 2), sqrt(5), [2, sqrt(5) / 2], exp(-pi)),  # v2
]:
    r = exp(k_ * th)
    vit = simplify(sqrt(r**2 + diff(r, th) ** 2))
    La = integrate(vit, (th, aa, 0))
    assert nul(limit(La, aa, -oo) - Lfin)
    choix(Lfin, leurres, limit(La, aa, -oo))  # le leurre « elle est infinie » est faux : limite finie
    tour = lambda b_: integrate(vit, (th, b_ - 2 * pi, b_))
    assert nul(simplify(tour(-2 * pi) / tour(0)) - raison)

# ===============================================================================================================
# geo-l3-courbes-longueur-ex3 (et v1, v2) : paramétrages normaux
# ===============================================================================================================
sig = symbols("sigma", real=True)
assert nul(integrate(cosh(u), (u, 0, t)) - sinh(t))
gch = (asinh(sig), sqrt(1 + sig**2))
assert nul(cosh(asinh(sig)) - sqrt(1 + sig**2)) and nul(d(gch, 1, sig)[0] ** 2 + d(gch, 1, sig)[1] ** 2 - 1)
for l in [cosh(sig), 1 + sig**2 / 2, sinh(sig)]:
    assert abs(float((l - sqrt(1 + sig**2)).subs(sig, Rational(3, 2)))) > 1e-6
# v1 : hélice, s = sqrt2 t, g = (cos(s/sqrt2), sin(s/sqrt2), s/sqrt2)
gh = (cos(sig / sqrt(2)), sin(sig / sqrt(2)), sig / sqrt(2))
assert nul(sum(c**2 for c in d(gh, 1, sig)) - 1)
for l in [sig, sqrt(2) * sig, sig / 2]:
    assert abs(float((l - sig / sqrt(2)).subs(sig, 1))) > 1e-6
# v2 : arc d'astroïde, s = 3/2 sin^2 t, g = ((1-2s/3)^{3/2}, (2s/3)^{3/2}) sur ]0, 3/2[
assert nul(integrate(3 * sin(u) * cos(u), (u, 0, t)) - Rational(3, 2) * sin(t) ** 2)
ga = ((1 - 2 * sig / 3) ** Rational(3, 2), (2 * sig / 3) ** Rational(3, 2))
assert nul(sum(c**2 for c in d(ga, 1, sig)) - 1)
assert nul(ga[1].subs(sig, Rational(3, 2) * sin(t) ** 2) - (sin(t) ** 2) ** Rational(3, 2))
for l in [(3 * sig / 2) ** Rational(3, 2), 2 * sig / 3, sin(sig) ** 3]:
    assert abs(float((l - ga[1]).subs(sig, Rational(1, 2)))) > 1e-6

# ===============================================================================================================
# geo-l3-courbes-longueur-ex4 : polygones inscrits dans le cercle
# ===============================================================================================================
n = symbols("n", positive=True)
ell = 2 * n * sin(pi / n)
assert limit((2 * pi - ell) / (pi**3 / (3 * n**2)), n, oo) == 1
for l in [pi**3 / (6 * n**2), pi**2 / n**2, pi**3 / (3 * n**3)]:
    assert limit((2 * pi - ell) / l, n, oo) != 1
assert abs(abs(mp.e ** (2j * mp.pi / 7) - 1) - 2 * mp.sin(mp.pi / 7)) < 1e-25
l96 = 192 * mp.sin(mp.pi / 96)
assert abs(l96 - mp.mpf("6.28206")) < 1e-5 and abs((2 * mp.pi - l96) - mp.mpf("1.1e-3")) < 0.05e-3

# ===============================================================================================================
# geo-l3-courbes-longueur-ex5 : périmètre de l'ellipse
# ===============================================================================================================
assert integrate(sin(t) ** 2, (t, 0, 2 * pi)) == pi and integrate(cos(t) ** 2, (t, 0, 2 * pi)) == pi
L21 = mp.quad(lambda v: mp.sqrt(4 * mp.sin(v) ** 2 + mp.cos(v) ** 2), [0, mp.pi / 2, mp.pi, 3 * mp.pi / 2, 2 * mp.pi])
assert abs(L21 - mp.mpf("9.68845")) < 1e-5
assert 3 * mp.pi < L21 < mp.pi * mp.sqrt(10)
assert abs(3 * mp.pi - mp.mpf("9.42477")) < 1e-5 and abs(mp.pi * mp.sqrt(10) - mp.mpf("9.93459")) < 1e-5
# leurres : encadrements faux (l'un des deux bords est violé)
assert not (L21 <= mp.pi * mp.sqrt(5))  # 3 pi <= L <= pi sqrt5
assert not (L21 <= 3 * mp.pi)  # 3pi/2 <= L <= 3 pi
assert not (mp.pi * mp.sqrt(10) <= L21)  # pi sqrt10 <= L <= 3 sqrt2 pi
# concavité : sqrt(lambda x + (1-lambda) y) >= lambda sqrt x + (1-lambda) sqrt y (contrôle numérique)
for v in [0.1, 0.7, 1.3, 2.9]:
    assert mp.sqrt(4 * mp.sin(v) ** 2 + mp.cos(v) ** 2) >= 2 * mp.sin(v) ** 2 + mp.cos(v) ** 2
# cercle (a = b) : égalité 2 pi a
assert 2 * mp.pi == mp.pi * (1 + 1) == mp.pi * mp.sqrt(2 * (1 + 1))

# ===============================================================================================================
# Carte geo-l3-courbes-frenet
# ===============================================================================================================
rr = symbols("r", positive=True)
assert nul(courbure((rr * cos(t), rr * sin(t))) - 1 / rr)
par = (t, t**2)
kpar = simplify(courbure(par))
assert nul(kpar - 2 / (1 + 4 * t**2) ** Rational(3, 2)) and kpar.subs(t, 0) == 2


def centre(f, t0):
    f1 = [c.subs(t, t0) for c in d(f)]
    k = courbure(f).subs(t, t0)
    v = sqrt(f1[0] ** 2 + f1[1] ** 2)
    N = (-f1[1] / v, f1[0] / v)
    return tuple(simplify(f[i].subs(t, t0) + N[i] / k) for i in range(2))


Cp = centre(par, t)
assert nul(Cp[0] + 4 * t**3) and nul(Cp[1] - Rational(1, 2) - 3 * t**2)
# contre-exemples : (t, t|t|) courbure ±2/(1+4t^2)^{3/2} ; (t^2,t^3) formule ; (cos t, -sin t) : kappa = -1
tp = symbols("tp", positive=True)
assert nul(courbure((tp, tp**2), tp).subs(tp, t) - 2 / (1 + 4 * t**2) ** Rational(3, 2))
assert nul(courbure((tp, -tp**2), tp).subs(tp, t) + 2 / (1 + 4 * t**2) ** Rational(3, 2))
assert nul(courbure((tp**2, tp**3), tp) - 6 / (tp * (4 + 9 * tp**2) ** Rational(3, 2)))
assert nul(courbure((cos(t), -sin(t))) + 1) and nul(courbure((cos(t), sin(t))) - 1)
# affirmations : ellipse en (a,0) : kappa = a/b^2 ; vitesse constante : |kappa| = ||f''|| / ||f'||^2
el = (a * cos(t), b * sin(t))
assert nul(courbure(el).subs(t, 0) - a / b**2)
vc = (rr * cos(t), rr * sin(t))
assert nul(Abs(courbure(vc)) - sqrt(sum(c**2 for c in d(vc, 2))) / sum(c**2 for c in d(vc)))
# trou 1 : leurres faux pour la parabole en t = 1 (||f''||/||f'||^2 = 2/5) et le cercle de rayon 2
k1 = courbure(par).subs(t, 1)
assert nul(k1 - 2 / (5 * sqrt(5)))
assert simplify(sqrt(sum(c**2 for c in d(par, 2))) / sum(c**2 for c in d(par))).subs(t, 1) != k1
c2 = (2 * cos(t), 2 * sin(t))
assert simplify(det2(d(c2), d(c2, 2)) / sum(c**2 for c in d(c2))) != Rational(1, 2) and det2(d(c2), d(c2, 2)).simplify() != Rational(1, 2)
# trou 2 : sur le cercle de rayon 2, dT/ds = kappa N (kappa = 1/2), donc ni -kappa N ni N/kappa
T2 = [c / 2 for c in d(c2)]
dTds = [diff(c, t) / 2 for c in T2]
N2 = [-T2[1], T2[0]]
assert all(nul(dTds[i] - Rational(1, 2) * N2[i]) for i in range(2))
assert not all(nul(dTds[i] + Rational(1, 2) * N2[i]) for i in range(2)) and not all(nul(dTds[i] - 2 * N2[i]) for i in range(2))
# démonstration : cercle osculateur c(h) : c(0) = M, c'(0) = T0, c''(0) = kappa0 N0, vitesse unité
k0, hh = symbols("k0 hh", nonzero=True, real=True)
Cx, Cy = symbols("Cx Cy", real=True)
T0, N0 = (1, 0), (0, 1)
cc = [[Cx, Cy][i] + (sin(k0 * hh) * T0[i] - cos(k0 * hh) * N0[i]) / k0 for i in range(2)]
assert [nul(cc[i].subs(hh, 0) - ([Cx, Cy][i] - N0[i] / k0)) for i in range(2)] == [True, True]
assert [simplify(diff(c, hh).subs(hh, 0)) for c in cc] == [1, 0] and [simplify(diff(c, hh, 2).subs(hh, 0)) for c in cc] == [0, k0]
assert nul(sum(diff(c, hh) ** 2 for c in cc) - 1)
# osculation sur la parabole en t = 1 (contrôle numérique) : g(h) - c(h) = O(h^3)
mp.mp.dps = 30
sabs = lambda tt: mp.quad(lambda w: mp.sqrt(1 + 4 * w**2), [1, tt])
def g_par(hv):
    tt = mp.findroot(lambda w: sabs(w) - hv, 1 + hv / mp.sqrt(5))
    return mp.matrix([tt, tt**2])
kk = 2 / mp.sqrt(5) ** 3
T0n, N0n = mp.matrix([1, 2]) / mp.sqrt(5), mp.matrix([-2, 1]) / mp.sqrt(5)
Cn = mp.matrix([1, 1]) + N0n / kk
def c_par(hv):
    return Cn + (mp.sin(kk * hv) * T0n - mp.cos(kk * hv) * N0n) / kk
e1, e2 = mp.norm(g_par(mp.mpf("1e-2")) - c_par(mp.mpf("1e-2"))), mp.norm(g_par(mp.mpf("1e-3")) - c_par(mp.mpf("1e-3")))
assert e1 / mp.mpf("1e-4") < 1e-1 and e2 / mp.mpf("1e-6") < 1e-2 and 500 < e1 / e2 < 2000  # décroissance en h^3

# ===============================================================================================================
# geo-l3-courbes-frenet-ex1 (et v1, v2)
# ===============================================================================================================
choix(2 / (5 * sqrt(5)), [Rational(2, 5), 5 * sqrt(5) / 2, 2 / sqrt(5)], courbure(par).subs(t, 1))
cub = (t, t**3)
choix(3 / (5 * sqrt(10)), [Rational(3, 5), 5 * sqrt(10) / 3, 6 / sqrt(10)], courbure(cub).subs(t, 1))
choix(2, [Rational(1, 2), Rational(1, 4), 4], courbure((2 * cos(t), sin(t))).subs(t, 0))
assert abs(float(2 / (5 * sqrt(5))) - 0.179) < 1e-3 and abs(float(5 * sqrt(5) / 2) - 5.59) < 1e-2
assert abs(float(3 / (5 * sqrt(10))) - 0.190) < 1e-3 and abs(float(5 * sqrt(10) / 3) - 5.27) < 1e-2

# ===============================================================================================================
# geo-l3-courbes-frenet-ex2 (et v1, v2) : centres de courbure
# ===============================================================================================================
def choix_point(bonne, leurres, valeur):
    assert tuple(simplify(c) for c in bonne) == tuple(simplify(c) for c in valeur)
    assert all(tuple(simplify(c) for c in l) != tuple(simplify(c) for c in valeur) for l in leurres)
    assert len({tuple(bonne)} | {tuple(l) for l in leurres}) == 4


R2 = Rational
choix_point((-4, R2(7, 2)), [(6, -R2(3, 2)), (4, R2(7, 2)), (-4, R2(5, 2))], centre(par, 1))
choix_point((-R2(1, 2), R2(5, 4)), [(R2(3, 2), -R2(3, 4)), (R2(1, 2), R2(5, 4)), (-R2(1, 2), R2(3, 4))], centre(par, R2(1, 2)))
choix_point((-4, R2(8, 3)), [(6, -R2(2, 3)), (4, R2(8, 3)), (-4, R2(5, 3))], centre(cub, 1))
# leurres « M - N/kappa » : symétriques du centre par rapport à M
assert (2 * 1 - (-4), 2 * 1 - R2(7, 2)) == (6, -R2(3, 2)) and (2 * 1 - (-4), 2 * 1 - R2(8, 3)) == (6, -R2(2, 3))
assert (2 * R2(1, 2) + R2(1, 2), 2 * R2(1, 4) - R2(5, 4)) == (R2(3, 2), -R2(3, 4))

# ===============================================================================================================
# geo-l3-courbes-frenet-ex3 (et v1, v2) : rayon de courbure en polaires
# ===============================================================================================================
def kpol(r):
    r1, r2 = diff(r, th), diff(r, th, 2)
    return (r**2 + 2 * r1**2 - r * r2) / (r**2 + r1**2) ** Rational(3, 2)


r = 1 + cos(th)
assert nul(expand_trig(r**2 + diff(r, th) ** 2 - 2 * (1 + cos(th))))
assert nul(expand_trig(r**2 + 2 * diff(r, th) ** 2 - r * diff(r, th, 2) - 3 * (1 + cos(th))))
# formule polaire = formule générale (cardioïde, en theta = 1)
fc = ((1 + cos(th)) * cos(th), (1 + cos(th)) * sin(th))
assert abs(float((courbure(fc, th) - kpol(r)).subs(th, 1))) < 1e-12
choix(Rational(4, 3), [Rational(3, 4), 2, Rational(3, 2)], 1 / kpol(r).subs(th, 0))
choix(Rational(9, 4), [Rational(4, 9), 3, Rational(27, 4)], 1 / kpol(2 + cos(th)).subs(th, 0))
choix(Rational(2, 3), [Rational(3, 2), Rational(1, 2), Rational(4, 3)], simplify(1 / kpol(r).subs(th, 2 * pi / 3)))

# ===============================================================================================================
# geo-l3-courbes-frenet-ex4 : C' = -kappa'/kappa^2 N ; développée de la parabole
# ===============================================================================================================
kp = courbure(par)
v_ = sqrt(1 + 4 * t**2)
Np = (-2 * t / v_, 1 / v_)
Cpt = [par[i] + Np[i] / kp for i in range(2)]
assert all(nul(diff(Cpt[i], t) + diff(kp, t) / kp**2 * Np[i]) for i in range(2))
XX, YY = -4 * t**3, Rational(1, 2) + 3 * t**2
assert nul(27 * XX**2 - 16 * (YY - Rational(1, 2)) ** 3)
for l in [16 * XX**2 - 27 * (YY - Rational(1, 2)) ** 3, 27 * XX**2 - 16 * YY**3, XX**2 - (YY - Rational(1, 2)) ** 3]:
    assert not nul(l)
assert [diff(XX, t).subs(t, 0), diff(YY, t).subs(t, 0)] == [0, 0] and diff(kp, t).subs(t, 0) == 0
# réciproque : (X, Y) avec Y >= 1/2 est atteint par t = -signe(X) sqrt((Y-1/2)/3)
for Xv in [R2(-4), R2(1, 2), R2(4)]:
    Yv = R2(1, 2) + 3 * (abs(Xv) / 4) ** R2(2, 3)
    tv = -__import__("sympy").sign(Xv) * sqrt((Yv - R2(1, 2)) / 3)
    assert nul(-4 * tv**3 - Xv) and nul(27 * Xv**2 - 16 * (Yv - R2(1, 2)) ** 3)

# ===============================================================================================================
# geo-l3-courbes-frenet-ex5 : angle total de la tangente pour kappa(s) = 1/(1+s^2)
# ===============================================================================================================
total = integrate(1 / (1 + s**2), (s, -oo, oo))
choix(pi, [pi / 2, 2 * pi], total)
assert integrate(1 / (1 + u**2), (u, 0, s)) == atan(s)
# question 3 : g' = e^{i phi} unitaire, g'' = i kappa g'
phi = Function("phi")(s)
gp = exp(Ii * phi)
assert nul(diff(gp, s) - Ii * diff(phi, s) * gp)

# ===============================================================================================================
# Carte geo-l3-courbes-bezier
# ===============================================================================================================
for nn in range(1, 8):
    assert nul(sum(bern(k, nn, t) for k in range(nn + 1)) - 1)
    # base de R_n[X] : matrice des coefficients inversible
    Mb = Matrix([[__import__("sympy").Poly(bern(k, nn, t), t).coeff_monomial(t**j) for j in range(nn + 1)] for k in range(nn + 1)])
    assert Mb.det() != 0
    # dérivée : B_{k,n}' = n (B_{k-1,n-1} - B_{k,n-1})
    for k in range(nn + 1):
        gauche_ = bern(k - 1, nn - 1, t) if k >= 1 else 0
        droite_ = bern(k, nn - 1, t) if k <= nn - 1 else 0
        assert nul(diff(bern(k, nn, t), t) - nn * (gauche_ - droite_))
# De Casteljau = Bernstein, et B'(t) = n (P_1^{(n-1)} - P_0^{(n-1)}) : points de contrôle symboliques, n <= 5
for nn in range(1, 6):
    Ps = symbols(f"p0:{nn + 1}")
    niv = [list(Ps)]
    while len(niv[-1]) > 1:
        q_ = niv[-1]
        niv.append([(1 - t) * q_[i] + t * q_[i + 1] for i in range(len(q_) - 1)])
    Bt = sum(bern(k, nn, t) * Ps[k] for k in range(nn + 1))
    assert nul(niv[-1][0] - Bt)
    assert nul(diff(Bt, t) - nn * (niv[-2][1] - niv[-2][0]))
    assert nul(diff(Bt, t).subs(t, 0) - nn * (Ps[1] - Ps[0]))
# exemple : P = (0,0),(1,2),(2,0) : B = (2t, 4t-4t^2), y = 2x - x^2, tangentes, B(1/2) = (1,1)
Bx = sum(bern(k, 2, t) * [0, 1, 2][k] for k in range(3))
By = sum(bern(k, 2, t) * [0, 2, 0][k] for k in range(3))
assert nul(Bx - 2 * t) and nul(By - (4 * t - 4 * t**2)) and nul(By - (2 * Bx - Bx**2))
assert bezier([(0, 0), (1, 2), (2, 0)], F(1, 2)) == (1, 1)
assert casteljau([(0, 0), (1, 2), (2, 0)], F(1, 2))[1] == [(F(1, 2), 1), (F(3, 2), 1)]
assert nul(diff(2 * Bx - Bx**2, t).subs(t, 1) / diff(Bx, t).subs(t, 1) + 2)  # pente -2 en P2
# contre-exemples : B(2) = (2,0) pour n = 1 ; P0 = P1 : B = t^2 (1,1) ; ordonnée 2t(1-t) <= 1/2
assert tuple(-1 * p + 2 * q for p, q in zip((0, 0), (1, 0))) == (2, 0)
assert nul(sum(bern(k, 2, t) * [0, 0, 1][k] for k in range(3)) - t**2)
assert __import__("sympy").maximum(2 * t * (1 - t), t, Interval(0, 1)) == Rational(1, 2)
# affirmation : tangente en P2 dirigée par P1P2 = (1,-2) et non par (P0P2) ; trou 2 : B'(0) = n P0P1
assert det2((1, -2), (2, 0)) != 0
# trou 3 : P_1^{(n-1)} et P_0^{(n-1)} diffèrent de B(t) en général (n = 2, t = 1/2)
niv = casteljau([(0, 0), (1, 2), (2, 0)], F(1, 2))
assert niv[1][1] != niv[2][0] and niv[1][0] != niv[2][0]

# ===============================================================================================================
# geo-l3-courbes-bezier-ex1 (et v1, v2)
# ===============================================================================================================
def choix_pt(bonne, leurres, val):
    assert tuple(F(c) for c in bonne) == tuple(val)
    assert all(tuple(F(c) for c in l) != tuple(val) for l in leurres) and len({bonne, *leurres}) == 4


choix_pt((1, 1), [(1, 2), (1, F(2, 3)), (1, F(1, 2))], bezier([(0, 0), (1, 2), (2, 0)], F(1, 2)))
choix_pt((F(3, 2), F(9, 4)), [(F(3, 2), 3), (F(3, 2), F(3, 2)), (F(3, 2), F(3, 4))], bezier([(0, 0), (0, 3), (3, 3), (3, 0)], F(1, 2)))
choix_pt((2, F(4, 3)), [(2, 2), (3, F(3, 2)), (1, 1)], bezier([(0, 0), (3, 3), (6, 0)], F(1, 3)))
assert bezier([(0, 0), (3, 3), (6, 0)], F(1, 2)) == (3, F(3, 2))
assert (sum(F(c) for c in (0, 1, 2)) / 3, sum(F(c) for c in (0, 2, 0)) / 3) == (1, F(2, 3))  # isobarycentre

# ===============================================================================================================
# geo-l3-courbes-bezier-ex2 (et v1, v2) : De Casteljau
# ===============================================================================================================
for P, tt, niveaux, sym, leurres in [
    ([(0, 0), (1, 2), (3, 2), (4, 0)], F(1, 4),
     [[(F(1, 4), F(1, 2)), (F(3, 2), 2), (F(13, 4), F(3, 2))], [(F(9, 16), F(7, 8)), (F(31, 16), F(15, 8))], [(F(29, 32), F(9, 8))]],
     (F(99, 32), F(9, 8)), [(F(99, 32), F(9, 8)), (F(9, 16), F(7, 8)), (2, F(3, 2))]),
    ([(0, 0), (0, 3), (3, 3), (3, 0)], F(1, 3),
     [[(0, 1), (1, 3), (3, 2)], [(F(1, 3), F(5, 3)), (F(5, 3), F(8, 3))], [(F(7, 9), 2)]],
     (F(20, 9), 2), [(F(20, 9), 2), (F(1, 3), F(5, 3)), (F(3, 2), F(9, 4))]),
    ([(0, 0), (2, 0), (2, 2), (0, 2)], F(1, 3),
     [[(F(2, 3), 0), (2, F(2, 3)), (F(4, 3), 2)], [(F(10, 9), F(2, 9)), (F(16, 9), F(10, 9))], [(F(4, 3), F(14, 27))]],
     (F(4, 3), F(40, 27)), [(F(4, 3), F(40, 27)), (F(10, 9), F(2, 9)), (F(3, 2), 1)]),
]:
    niv = casteljau(P, tt)
    assert [[tuple(F(c) for c in p) for p in l] for l in niveaux] == niv[1:]
    assert niv[-1][0] == bezier(P, tt) and bezier(P, 1 - tt) == sym
    choix_pt(niveaux[-1][0], leurres, bezier(P, tt))

# ===============================================================================================================
# geo-l3-courbes-bezier-ex3 (et v1, v2) : point d'inflexion
# ===============================================================================================================
for ys, ypoly, t0, pt, leurres in [
    ((0, 1, -1, 0), 3 * (2 * t**3 - 3 * t**2 + t), Rational(1, 2), (Rational(3, 2), 0),
     [(1, 1), (Rational(3, 2), Rational(1, 4)), (2, -1)]),
    ((0, 2, -1, 0), 3 * (3 * t**3 - 5 * t**2 + 2 * t), Rational(5, 9), (Rational(5, 3), Rational(20, 81)),
     [(Rational(3, 2), Rational(3, 8)), (Rational(5, 3), Rational(2, 9)), (Rational(4, 3), Rational(20, 81))]),
    ((0, -1, 2, 0), 3 * (-3 * t**3 + 4 * t**2 - t), Rational(4, 9), (Rational(4, 3), Rational(20, 81)),
     [(Rational(3, 2), Rational(3, 8)), (Rational(5, 3), Rational(20, 81)), (Rational(4, 3), Rational(2, 9))]),
]:
    Bx = sum(bern(k, 3, t) * k for k in range(4))
    By = sum(bern(k, 3, t) * ys[k] for k in range(4))
    assert nul(Bx - 3 * t) and nul(By - ypoly)
    D = det2((diff(Bx, t), diff(By, t)), (diff(Bx, t, 2), diff(By, t, 2)))
    assert solveset(D, t, S.Reals) == {t0} and signe_change(D, t0)
    val = (Bx.subs(t, t0), By.subs(t, t0))
    assert val == pt and all(l != val for l in leurres) and len({pt, *leurres}) == 4
# ex3 : y(1/2) = 3/8 pour v1 et v2 (leurre B(1/2)), centre de symétrie (3/2, 0) pour l'exercice
assert sum(bern(k, 3, t) * (0, 2, -1, 0)[k] for k in range(4)).subs(t, R2(1, 2)) == R2(3, 8)
assert sum(bern(k, 3, t) * (0, -1, 2, 0)[k] for k in range(4)).subs(t, R2(1, 2)) == R2(3, 8)
yb = sum(bern(k, 3, t) * (0, 1, -1, 0)[k] for k in range(4))
assert nul(yb.subs(t, 1 - t) + yb)

# ===============================================================================================================
# geo-l3-courbes-bezier-ex4 : subdivision
# ===============================================================================================================
for nn in range(1, 8):
    for i in range(nn + 1):
        for j in range(i, nn + 1):
            assert comb(nn, j) * comb(j, i) == comb(nn, i) * comb(nn - i, j - i)
lam = symbols("lam")
for nn in range(1, 6):
    Ps = symbols(f"p0:{nn + 1}")
    Q = [sum(bern(i, j, lam) * Ps[i] for i in range(j + 1)) for j in range(nn + 1)]
    gauche_ = sum(bern(j, nn, u) * Q[j] for j in range(nn + 1))
    droite_ = sum(bern(i, nn, lam * u) * Ps[i] for i in range(nn + 1))
    assert nul(expand(gauche_ - droite_))
niv = casteljau([(0, 0), (2, 4), (4, 4), (6, 0)], F(1, 2))
Qs = [niv[j][0] for j in range(4)]
assert Qs == [(0, 0), (1, 2), (2, 3), (3, 3)]
choix_pt((2, 3), [(1, 2), (3, 3), (2, 4)], Qs[2])

# ===============================================================================================================
# geo-l3-courbes-bezier-ex5 : polynômes de Bernstein d'une fonction
# ===============================================================================================================
for nn in range(1, 9):
    assert nul(sum(k * bern(k, nn, t) for k in range(nn + 1)) - nn * t)
    assert nul(sum(k * (k - 1) * bern(k, nn, t) for k in range(nn + 1)) - nn * (nn - 1) * t**2)
    assert nul(sum(Rational(k, nn) ** 2 * bern(k, nn, t) for k in range(nn + 1)) - (t**2 + t * (1 - t) / nn))
    assert nul(sum((Rational(k, nn) - t) ** 2 * bern(k, nn, t) for k in range(nn + 1)) - t * (1 - t) / nn)
nn = symbols("N", positive=True)
assert __import__("sympy").maximum(t * (1 - t) / nn, t, Interval(0, 1)) == 1 / (4 * nn)
for l in [1 / nn, 1 / (2 * nn), 1 / (4 * nn**2)]:
    assert not nul(l - 1 / (4 * nn))
# question 4 : contrôle numérique de |B_n(f) - f| <= K/(2 sqrt n) pour f(x) = |x - 1/3| (K = 1)
for nv in [1, 3, 10, 40]:
    for tv in [mp.mpf(k) / 37 for k in range(38)]:
        Bn = sum(abs(mp.mpf(k) / nv - mp.mpf(1) / 3) * mp.binomial(nv, k) * tv**k * (1 - tv) ** (nv - k) for k in range(nv + 1))
        assert abs(Bn - abs(tv - mp.mpf(1) / 3)) <= 1 / (2 * mp.sqrt(nv)) + mp.mpf(10) ** -25

# ===============================================================================================================
# QCM
# ===============================================================================================================
# q1 : voir plus haut ((t^3, t^4) : p=3, q=4, allure ordinaire)
assert nature(*pq((t**3, t**4), 0)) == "ordinaire"
# q2 : (s^3, s^6) a une vitesse nulle en 0
assert [diff(c, s).subs(s, 0) for c in (s**3, s**6)] == [0, 0]
# q3 : cardioïde singulière seulement en theta = pi (modulo 2 pi) : rho = rho' = 0
r = 1 + cos(th)
assert solveset(r**2 + diff(r, th) ** 2, th, Interval.Ropen(0, 2 * pi)) == {pi}
# q4 : la cubique gauche n'est pas plane
assert det3(d(Cu), d(Cu, 2), d(Cu, 3)) == 12
# q5 : longueur 4 pi
choix(4 * pi, [2 * pi, 0, 8 * pi], integrate(1, (t, 0, 4 * pi)))
# q6 : kappa = -1/3 pour (3 cos t, -3 sin t)
choix(-Rational(1, 3), [Rational(1, 3), -3, 3], simplify(courbure((3 * cos(t), -3 * sin(t)))))
# q7 : (t, t^4) : kappa(0) = 0 mais q = 4 (pas d'inflexion) ; kappa(t0) != 0 => (p,q) = (1,2)
assert courbure((t, t**4)).subs(t, 0) == 0 and pq((t, t**4), 0) == (1, 4) and nature(1, 4) == "ordinaire"
# q8 : B'(0) = 2 (P1 - P0) = (2,2) ; B(t) = (2t, 2t(1-t))
Bq = [sum(bern(k, 2, t) * P[i] for k, P in enumerate([(0, 0), (1, 1), (2, 0)])) for i in range(2)]
assert nul(Bq[0] - 2 * t) and nul(Bq[1] - 2 * t * (1 - t)) and [diff(c, t).subs(t, 0) for c in Bq] == [2, 2]
assert all(l != (2, 2) for l in [(1, 1), (2, 0), (1, 0)])

print("geo_l3_courbes : toutes les vérifications passent.")
