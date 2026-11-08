"""Vérifications calculatoires du chapitre « cd-l3-differentiabilite » (calcul différentiel, L3) : cartes, exercices et
QCM de content/*/topologie-calcul-diff/l3-cd-l3-differentiabilite.json et content/qcm/cd-l3-differentiabilite.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question vérifiée.
Convention : d2(f, x, y) = ∂/∂x(∂f/∂y). Calculs exacts en SymPy (limites, dérivées, simplifications), complétés par des
tirages aléatoires pour les majorations. Lancer : python3 verification/cd_l3_differentiabilite.py.
"""

import random

from sympy import (Abs, lambdify, Matrix, Rational, S, atan, cos, diff, exp, expand, factor, limit, ln, nsimplify, oo, pi, simplify,
                   sin, solve, sqrt, symbols, N, Function, eye, zeros)

random.seed(2)
x, y, z, t, h, k, s, u, v = symbols("x y z t h k s u v", real=True)
r, th = symbols("r theta", positive=True)
v1, v2 = symbols("v1 v2", real=True)


def en(f, a):
    """Valeur de f au point a = (a1, a2[, a3])."""
    return simplify(f.subs({x: a[0], y: a[1], **({z: a[2]} if len(a) > 2 else {})}, simultaneous=True))


def grad(f, vs=(x, y)):
    return [diff(f, w) for w in vs]


def taux(f, a, w):
    """(f(a + t w) - f(a)) / t pour une fonction de (x, y)."""
    return (f.subs({x: a[0] + t * w[0], y: a[1] + t * w[1]}, simultaneous=True) - f.subs({x: a[0], y: a[1]})) / t


def dir_deriv(f, a, w):
    return limit(simplify(taux(f, a, w)), t, 0)


def d2(f, p, q):
    """∂/∂p (∂f/∂q)."""
    return diff(diff(f, q), p)


def faux(valeur, *leurres):
    for l in leurres:
        assert simplify(S(l) - S(valeur)) != 0, (valeur, l)


def tirages(n=3000, borne=1.0):
    for _ in range(n):
        yield random.uniform(-borne, borne), random.uniform(-borne, borne)


# ===========================================================================
# Carte cd-l3-differentiabilite-differentielle
# ===========================================================================
# exemple : f = x^2 y en (1,2), reste 2h1^2 + 2h1h2 + h1^2 h2, majoré par 5||h||_inf^2 si ||h||_inf <= 1.
f = x**2 * y
reste = expand(f.subs({x: 1 + h, y: 2 + k}) - 2 - 4 * h - k)
assert reste == expand(2 * h**2 + 2 * h * k + h**2 * k)
for a1, a2 in tirages():
    assert abs(float(reste.subs({h: a1, k: a2}))) <= 5 * max(abs(a1), abs(a2)) ** 2 + 1e-12
assert 4 * 1 + 1 * 1 == 5  # D_(1,1) f(1,2)
# exemple : (M+H)^2 = M^2 + MH + HM + H^2, et MH + HM != 2MH pour des matrices qui ne commutent pas.
M = Matrix([[0, 1], [0, 0]])
H = Matrix([[0, 0], [1, 0]])
assert (M + H) ** 2 == M**2 + M * H + H * M + H**2 and M * H + H * M != 2 * M * H
# exemple : ||a+h||^2 = ||a||^2 + 2<a,h> + ||h||^2.
a1, a2, b1, b2 = symbols("a1 a2 b1 b2", real=True)
assert expand((a1 + b1) ** 2 + (a2 + b2) ** 2 - (a1**2 + a2**2) - 2 * (a1 * b1 + a2 * b2) - (b1**2 + b2**2)) == 0

# contre-exemple 1 (et affirmation 1, cd-l3-differentiabilite-ex3) : xy^2/(x^2+y^4).
f1 = x * y**2 / (x**2 + y**4)
q = simplify(f1.subs({x: t * v1, y: t * v2}) / t)
assert simplify(q - v1 * v2**2 / (v1**2 + t**2 * v2**4)) == 0
assert limit(q.subs({v1: 2, v2: 3}), t, 0) == Rational(9, 2)  # v2^2/v1
assert simplify(q.subs(v1, 0)) == 0
assert simplify(f1.subs(x, y**2)) == Rational(1, 2)

# contre-exemple 2 : x^3/(x^2+y^2) : |f| <= |x| (car |x|(x^2+y^2) - |x|^3 = |x| y^2 >= 0), f(tv) = t f(v).
f2 = x**3 / (x**2 + y**2)
assert simplify(f2.subs({x: t * v1, y: t * v2}) - t * f2.subs({x: v1, y: v2})) == 0
assert [f2.subs({x: 1, y: 0}), f2.subs({x: 0, y: 1}), f2.subs({x: 1, y: 1})] == [1, 0, Rational(1, 2)]
for a1_, a2_ in tirages():
    if (a1_, a2_) != (0, 0):
        assert abs(float(f2.subs({x: a1_, y: a2_}))) <= abs(a1_) + 1e-12

# contre-exemple 3 : xy/sqrt(x^2+y^2) : |f| <= sqrt(x^2+y^2)/2, f(t,t)/t = ±1/sqrt(2).
f3 = x * y / sqrt(x**2 + y**2)
for a1_, a2_ in tirages():
    if (a1_, a2_) != (0, 0):
        assert abs(float(f3.subs({x: a1_, y: a2_}))) <= 0.5 * (a1_**2 + a2_**2) ** 0.5 + 1e-12
assert limit(f3.subs({x: t, y: t}) / t, t, 0, "+") == 1 / sqrt(2)
assert limit(f3.subs({x: t, y: t}) / t, t, 0, "-") == -1 / sqrt(2)
assert f3.subs({x: t, y: 0}) == 0 and f3.subs({x: 0, y: t}) == 0

# affirmation 3 : x^4 y/(x^6+y^2) : f(tv)/t -> 0 pour tout v, |f| <= |x|/2, f(x, x^3) = x/2.
f4 = x**4 * y / (x**6 + y**2)
q4 = simplify(f4.subs({x: t * v1, y: t * v2}) / t)
assert simplify(q4 - t**2 * v1**4 * v2 / (t**4 * v1**6 + v2**2)) == 0
for a, b in [(1, 1), (2, -3), (0, 1), (1, 0), (-5, 7)]:
    assert limit(q4.subs({v1: a, v2: b}), t, 0) == 0
assert expand((x**3 - y) ** 2) == expand(x**6 + y**2 - 2 * x**3 * y)  # x^6 + y^2 >= 2|x^3 y|
assert simplify(f4.subs(y, x**3)) == x / 2

# trou 1 (leurres) : x^2 sin(1/x^2) est dérivable en 0 (|f(h)| <= h^2) mais f' est non bornée près de 0.
g = x**2 * sin(1 / x**2)
gp = diff(g, x)
for kk in [1, 10, 100]:
    pt = 1 / sqrt(2 * kk * pi)
    assert simplify(gp.subs(x, pt) + 2 * sqrt(2 * kk * pi)) == 0  # f'(pt) = -2 sqrt(2k pi), non borné
# trou 2 (leurres) : f(x) = -x sur R : D_1 f = -1 alors que |||df||| * 1 = 1 ; f = x, v = 2 : D_v f = 2, df(v/|v|) = 1.
assert -1 != 1 and 2 != 1

# ===========================================================================
# Exercices de la carte differentielle
# ===========================================================================
# cd-l3-differentiabilite-differentielle-ex1 (et v1, v2) : D_v f(a) par le taux d'accroissement.
for f_, a_, w_, dev, rep, leurres in [
    (x**2 * y - y**3, (1, 2), (3, -1), -6 + 23 * t + 6 * t**2 - 8 * t**3, 23, (1, 12, -23)),
    (x**3 + x * y**2, (2, 1), (1, 2), 10 + 21 * t + 18 * t**2 + 5 * t**3, 21, (13, 17, 42)),
    (x**2 * y + 2 * x * y**2, (1, -1), (2, 3), 1 - 9 * t + 2 * t**2 + 48 * t**3, -9, (9, -3, -6)),
]:
    assert expand(f_.subs({x: a_[0] + t * w_[0], y: a_[1] + t * w_[1]}, simultaneous=True)) == expand(dev)
    assert dir_deriv(f_, a_, w_) == rep
    gx, gy = [en(c, a_) for c in grad(f_)]
    assert gx * w_[0] + gy * w_[1] == rep
    faux(rep, *leurres)
# ex1 : leurre 1 = 4*3 + (-11)*1 (erreur de signe sur v2) ; v1 : 13 = dérivée partielle en x seule.
assert 4 * 3 + (-11) * 1 == 1 and 3 * 4 + 1 == 13

# cd-l3-differentiabilite-differentielle-ex2 (et v1, v2) : fonctions homogènes de degré 1.
for f_, w_, rep, leurres, lin in [
    (x**3 / (x**2 + y**2), (1, 1), Rational(1, 2), (1, 0, 2), 1 + 0),
    (x**2 * y / (x**2 + y**2), (1, 2), Rational(2, 5), (0, 2, Rational(1, 2)), 0),
    ((x**3 + y**3) / (x**2 + y**2), (2, 1), Rational(9, 5), (3, 0, 1), 2 * 1 + 1 * 1),
]:
    assert simplify(f_.subs({x: t * x, y: t * y}, simultaneous=True) - t * f_) == 0
    assert dir_deriv(f_, (0, 0), w_) == rep
    assert rep != lin  # non linéarité de v -> D_v f(0)
    faux(rep, *leurres)
# v2 : |x^3 + y^3| <= (|x|+|y|)(x^2+y^2).
for a1_, a2_ in tirages():
    assert abs(a1_**3 + a2_**3) <= (abs(a1_) + abs(a2_)) * (a1_**2 + a2_**2) + 1e-12

# cd-l3-differentiabilite-differentielle-ex3 (et v1, v2) : dérivées selon tout vecteur sans continuité.
fA = x * y**2 / (x**2 + y**4)
assert simplify(fA.subs(x, y**2)) == Rational(1, 2)
faux(Rational(1, 2), 0, 1)
fB = x * y**2 / (x**2 + 2 * y**4)  # v1
assert simplify(fB.subs({x: t * v1, y: t * v2}) / t - v1 * v2**2 / (v1**2 + 2 * t**2 * v2**4)) == 0
assert simplify(fB.subs(x, y**2)) == Rational(1, 3)
faux(Rational(1, 3), Rational(1, 2), 0, 1)
fC = x**2 * y / (x**4 + 3 * y**2)  # v2
qC = simplify(fC.subs({x: t * v1, y: t * v2}) / t)
assert simplify(qC - v1**2 * v2 / (t**2 * v1**4 + 3 * v2**2)) == 0
assert limit(qC.subs({v1: 3, v2: 2}), t, 0) == Rational(9, 6)  # v1^2/(3 v2)
assert simplify(fC.subs(y, x**2)) == Rational(1, 4)
faux(Rational(1, 4), Rational(1, 3), 0, Rational(1, 2))

# cd-l3-differentiabilite-differentielle-ex4 : g_alpha = (x^3 + alpha y^3)/(x^2+y^2), homogène de degré 1, jamais
# linéaire : l'additivité sur (1,1) impose alpha = -1, sur (1,-1) impose alpha = 1.
al = symbols("alpha", real=True)
ga = (x**3 + al * y**3) / (x**2 + y**2)
assert simplify(ga.subs({x: t * x, y: t * y}, simultaneous=True) - t * ga) == 0
assert solve(ga.subs({x: 1, y: 1}) - (ga.subs({x: 1, y: 0}) + ga.subs({x: 0, y: 1})), al) == [-1]
assert solve(ga.subs({x: 1, y: -1}) - (ga.subs({x: 1, y: 0}) - ga.subs({x: 0, y: 1})), al) == [1]
for a1_, a2_ in tirages(500):
    for alpha in [-3, -1, 0, 1, 2.5]:
        val = (a1_**3 + alpha * a2_**3) / (a1_**2 + a2_**2)
        assert abs(val) <= (1 + abs(alpha)) * (a1_**2 + a2_**2) ** 0.5 + 1e-9

# cd-l3-differentiabilite-differentielle-ex5 : normes.
# Q1 : identité du numérateur et majoration 2||h||^2/||a||.
na, nah, pah, nh2 = symbols("na nah pah nh2", real=True)
assert expand((2 * pah + nh2) * na - pah * (nah + na) - (pah * (na - nah) + nh2 * na)) == 0
for _ in range(2000):
    A = [random.uniform(-2, 2) for _ in range(3)]
    Hh = [random.uniform(-1, 1) * 10 ** random.uniform(-4, 0) for _ in range(3)]
    nA = sum(c * c for c in A) ** 0.5
    if nA < 1e-3:
        continue
    nAH = sum((a_ + b_) ** 2 for a_, b_ in zip(A, Hh)) ** 0.5
    nH = sum(c * c for c in Hh) ** 0.5
    ecart = nAH - nA - sum(a_ * b_ for a_, b_ in zip(A, Hh)) / nA
    assert abs(ecart) <= 2 * nH**2 / nA + 1e-12
# Q4 : ||.||_inf est différentiable en (2,-1) (|a1| != |a2|), pas en (1,1), (1,-1), (0,0) : en a = (c eps1, ±c),
# t -> ||a + t(eps1,0)||_inf vaut c + t pour t >= 0 et c pour t <= 0.
def ninf(p):
    return max(abs(p[0]), abs(p[1]))
for a_ in [(1, 1), (1, -1)]:
    droite = (ninf((a_[0] + 1e-6, a_[1])) - ninf(a_)) / 1e-6
    gauche = (ninf((a_[0] - 1e-6, a_[1])) - ninf(a_)) / -1e-6
    assert abs(droite - 1) < 1e-6 and abs(gauche) < 1e-6
assert abs(2) != abs(-1)

# ===========================================================================
# Carte cd-l3-differentiabilite-jacobienne
# ===========================================================================
F = Matrix([x**2 * y, x - y**2, x * y])
J = F.jacobian([x, y])
assert J.subs({x: 1, y: 2}) == Matrix([[4, 1], [1, -4], [2, 1]])  # exemple et ex1
assert J.subs({x: 1, y: 2}) * Matrix([1, 1]) == Matrix([5, -3, 3])
# exemple : coordonnées polaires, déterminant r ; xyz.
Phi = Matrix([r * cos(th), r * sin(th)])
assert simplify(Phi.jacobian([r, th]).det()) == r
assert Matrix([x * y * z]).jacobian([x, y, z]) == Matrix([[y * z, x * z, x * y]])
# contre-exemples : xy/(x^2+y^2) (dérivées partielles nulles, f(t,t) = 1/2) ; x^3/(x^2+y^2) (D_(1,1) = 1/2 != 1) ;
# (y,0) et sa transposée.
f5 = x * y / (x**2 + y**2)
assert f5.subs({x: t, y: 0}) == 0 and f5.subs({x: 0, y: t}) == 0 and simplify(f5.subs({x: t, y: t})) == Rational(1, 2)
assert dir_deriv(x**3 / (x**2 + y**2), (0, 0), (1, 1)) == Rational(1, 2)
assert Matrix([y, 0]).jacobian([x, y]) == Matrix([[0, 1], [0, 0]])
assert Matrix([0, x]).jacobian([x, y]) == Matrix([[0, 1], [0, 0]]).T
# affirmation 3 : jacobienne de x -> Ax égale à A, différente de sa transposée pour A = [[0,1],[0,0]].
A = Matrix([[0, 1], [0, 0]])
assert (A * Matrix([x, y])).jacobian([x, y]) == A != A.T
# trou 3 (leurres) : xy/(x^2+y^2) n'a pas de dérivée selon (1,1) en 0 : f(t,t)/t = 1/(2t).
assert simplify(f5.subs({x: t, y: t}) / t) == 1 / (2 * t)

# cd-l3-differentiabilite-jacobienne-ex1 (et v1, v2).
for Fv, a_, hv, rep, leurres in [
    (Matrix([x**2 * y, x - y**2, x * y]), (1, 2), (1, 1), (5, -3, 3), [(4, 1, 2), (5, 3, 3), (1, -4, 1)]),
    (Matrix([x * y**2, x**2 + y, x - y]), (2, 1), (1, -1), (-3, 3, 2), [(1, 4, 1), (5, 5, 0), (3, -3, 2)]),
    (Matrix([x + y**2, x * y, y**3]), (1, -1), (2, 1), (0, -1, 3), [(1, -1, 0), (-1, 0, 3), (0, 1, 3)]),
]:
    Ja = Fv.jacobian([x, y]).subs({x: a_[0], y: a_[1]})
    res = tuple(Ja * Matrix(hv))
    assert res == rep
    for l in leurres:
        assert l != rep
# v1 : la première colonne de J vaut (1,4,1) (leurre), la somme des colonnes (5,5,0).
Jv1 = Matrix([x * y**2, x**2 + y, x - y]).jacobian([x, y]).subs({x: 2, y: 1})
assert tuple(Jv1[:, 0]) == (1, 4, 1) and tuple(Jv1[:, 0] + Jv1[:, 1]) == (5, 5, 0)

# cd-l3-differentiabilite-jacobienne-ex2 (et v1, v2) : déterminants jacobiens.
assert simplify(Phi.jacobian([r, th]).det().subs({r: 3, th: pi / 4})) == 3
faux(3, 9, 1, 3 * sqrt(2) / 2)
P1 = Matrix([u**2 - v**2, 2 * u * v])
assert expand(P1.jacobian([u, v]).det()) == 4 * u**2 + 4 * v**2 and P1.jacobian([u, v]).det().subs({u: 1, v: 2}) == 20
assert 2 * 2 - 4 * 4 == -12  # leurre : erreur de signe dans le déterminant
faux(20, -12, 10, 5)
P2 = Matrix([exp(x) * cos(y), exp(x) * sin(y)])
assert simplify(P2.jacobian([x, y]).det()) == exp(2 * x)
assert simplify(exp(2 * x).subs(x, ln(2))) == 4
faux(4, 2, 1, 8)

# cd-l3-differentiabilite-jacobienne-ex3 (et v1, v2) : d(M^2)(H) = MH + HM ; matrice dans (E11, E12, E21, E22).
def mat_carre(d1, d2):
    Mm = Matrix([[d1, 0], [0, d2]])
    col = []
    for (i, j) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        E = zeros(2, 2)
        E[i, j] = 1
        img = Mm * E + E * Mm
        col.append([img[0, 0], img[0, 1], img[1, 0], img[1, 1]])
    return Matrix(col).T
for d1_, d2_, rep, leurres in [(1, 2, 72, (64, 24, 10)), (1, 3, 192, (144, 48, 16)), (2, -1, -8, (64, 8, 4))]:
    Mat = mat_carre(d1_, d2_)
    assert Mat == Matrix.diag(2 * d1_, d1_ + d2_, d1_ + d2_, 2 * d2_)
    assert Mat.det() == rep
    faux(rep, *leurres)
    # leurre « 2MH » : diag(2d1, 2d1, 2d2, 2d2).
    assert Matrix.diag(2 * d1_, 2 * d1_, 2 * d2_, 2 * d2_).det() == leurres[0]
# ex3 : (M+H)^2 - M^2 - MH - HM = H^2 (matrices génériques).
m = Matrix(2, 2, symbols("m0:4"))
hh = Matrix(2, 2, symbols("h0:4"))
assert expand((m + hh) ** 2 - m**2 - m * hh - hh * m - hh**2) == zeros(2, 2)

# cd-l3-differentiabilite-jacobienne-ex4 (et v1, v2) : d(det)(M)(H) = tr(com(M)^T H).
def ddet(Mm, Hm):
    return diff((Mm + t * Hm).det(), t).subs(t, 0)
mg = Matrix(3, 3, symbols("g0:9"))
hg = Matrix(3, 3, symbols("k0:9"))
assert expand(ddet(mg, hg) - (mg.adjugate() * hg).trace()) == 0  # adjugate = com(M)^T
for Mm, Hm, rep, leurres in [
    (Matrix([[1, 2], [3, 4]]), Matrix([[0, 1], [1, 0]]), -5, (-2, 0, 5)),
    (Matrix([[2, 1], [1, 3]]), eye(2), 5, (2, 10, -5)),
    (Matrix([[1, 1], [0, 2]]), Matrix([[0, 0], [1, 0]]), -1, (1, 0, 2)),
]:
    assert ddet(Mm, Hm) == rep == (Mm.adjugate() * Hm).trace()
    assert Mm.det() * (Mm.inv() * Hm).trace() == rep
    faux(rep, *leurres)
assert Matrix([[1, 2], [3, 4]]).adjugate().T == Matrix([[4, -3], [-2, 1]])  # com(M) de la question 4
assert Matrix([[2, 1], [1, 3]]).det() * eye(2).trace() == 10  # leurre v1 : det(M) tr(H)

# ===========================================================================
# Carte cd-l3-differentiabilite-classe-c1
# ===========================================================================
fe = x * exp(x * y)
assert simplify(diff(fe, x) - (1 + x * y) * exp(x * y)) == 0 and simplify(diff(fe, y) - x**2 * exp(x * y)) == 0
assert [en(c, (1, 0)) for c in grad(fe)] == [1, 1]
f6 = x**2 * y**2 / (x**2 + y**2)
assert simplify(diff(f6, x) - 2 * x * y**4 / (x**2 + y**2) ** 2) == 0
# contre-exemple : (x^2+y^2) sin(1/r) ; dérivée partielle en x sur le demi-axe x > 0.
xp = symbols("xp", positive=True)
fs = (x**2 + y**2) * sin(1 / sqrt(x**2 + y**2))
dxs = simplify(diff(fs, x).subs({x: xp, y: 0}))
assert simplify(dxs - (2 * xp * sin(1 / xp) - cos(1 / xp))) == 0
for kk in [1, 5, 50]:
    assert simplify(dxs.subs(xp, 1 / (2 * kk * pi))) == -1
    assert simplify(dxs.subs(xp, 1 / ((2 * kk + 1) * pi))) == 1
# contre-exemple : x^3/(x^2+y^2) a des dérivées partielles bornées (bornes 9/8 <= 3 et 3 sqrt(3)/8 <= 1).
dgx = simplify(diff(x**3 / (x**2 + y**2), x))
dgy = simplify(diff(x**3 / (x**2 + y**2), y))
assert simplify(dgx - x**2 * (x**2 + 3 * y**2) / (x**2 + y**2) ** 2) == 0
assert simplify(dgy + 2 * x**3 * y / (x**2 + y**2) ** 2) == 0
cc = symbols("c", nonnegative=True)
assert simplify(dgx.subs({x: cos(th), y: sin(th)}) - cos(th) ** 2 * (3 - 2 * cos(th) ** 2)) == 0
assert solve(diff(cc * (3 - 2 * cc), cc), cc) == [Rational(3, 4)]
assert Rational(3, 4) * (3 - 2 * Rational(3, 4)) == Rational(9, 8) <= 3
assert N(3 * sqrt(3) / 8) < 1
for a1_, a2_ in tirages():
    if (a1_, a2_) != (0, 0):
        assert abs(float(dgy.subs({x: a1_, y: a2_}))) <= 3 * 3**0.5 / 8 + 1e-12
        assert 0 <= float(dgx.subs({x: a1_, y: a2_})) <= 9 / 8 + 1e-12
# trou 2 (leurres) : f(x,y) = phi(x - y) vérifie df(e1 + e2) = 0 ; f(x,y) = phi(y) vérifie df(e1) = 0.
phi = Function("phi")
assert simplify(diff(phi(x - y), x) + diff(phi(x - y), y)) == 0 and diff(phi(y), x) == 0

# cd-l3-differentiabilite-classe-c1-ex1 (et v1, v2).
for f_, a_, hv, rep, leurres in [
    (x * exp(x * y), (1, 0), (2, 3), 5, (2, 3, 6)),
    (y * exp(x) + x * y**3, (0, 1), (1, 2), 4, (2, 3, 5)),
    (cos(x) + x * y**2, (0, 2), (3, -1), 12, (4, 13, -12)),
]:
    gx, gy = [en(c, a_) for c in grad(f_)]
    assert gx * hv[0] + gy * hv[1] == rep
    faux(rep, *leurres)

# cd-l3-differentiabilite-classe-c1-ex2 (et v1, v2) : valeurs des dérivées partielles oscillantes.
fc = (x**2 + y**2) * cos(1 / sqrt(x**2 + y**2))
dxc = simplify(diff(fc, x).subs({x: xp, y: 0}))
assert simplify(dxc - (2 * xp * cos(1 / xp) + sin(1 / xp))) == 0
for kk in [1, 3, 20]:
    assert simplify(dxs.subs(xp, 1 / (2 * kk * pi))) == -1  # ex
    assert simplify(dxc.subs(xp, 2 / ((4 * kk + 1) * pi))) == 1  # v1
    yp = 1 / ((2 * kk + 1) * pi)
    assert simplify(diff(fs, y).subs({x: 0, y: yp})) == 1  # v2
faux(-1, 0, 1)
faux(1, -1, 0, 2)

# cd-l3-differentiabilite-classe-c1-ex3 (et v1, v2) : plus petit k donnant la classe C^1.
kk_ = symbols("k", positive=True, integer=True)
# ex : f_k = x^k y/(x^2+y^2) ; dérivées partielles hors de l'origine.
for kv in [3, 4, 5]:
    fk = x**kv * y / (x**2 + y**2)
    assert simplify(diff(fk, x) - (kv * x ** (kv - 1) * y * (x**2 + y**2) - 2 * x ** (kv + 1) * y) / (x**2 + y**2) ** 2) == 0
    assert simplify(diff(fk, y) - x**kv * (x**2 - y**2) / (x**2 + y**2) ** 2) == 0
assert simplify((x * y / (x**2 + y**2)).subs({x: t, y: t})) == Rational(1, 2)  # f_1 discontinue
assert simplify((x**2 * y / (x**2 + y**2)).subs({x: t, y: t}) / t) == Rational(1, 2)  # f_2 non différentiable
faux(3, 1, 2, 4)
# v1 : f_k = x y^k/(x^4+y^4) ; (x^2+y^2)^2 <= 2(x^4+y^4).
assert expand(2 * (x**4 + y**4) - (x**2 + y**2) ** 2) == expand((x**2 - y**2) ** 2)
for kv in [5, 6]:
    fk = x * y**kv / (x**4 + y**4)
    assert simplify(diff(fk, x) - y**kv * (y**4 - 3 * x**4) / (x**4 + y**4) ** 2) == 0
    assert simplify(diff(fk, y) - x * y ** (kv - 1) * (kv * (x**4 + y**4) - 4 * y**4) / (x**4 + y**4) ** 2) == 0
for kv in [1, 2, 3]:
    assert simplify((x * y**kv / (x**4 + y**4)).subs({x: t, y: t}) - t ** (kv - 3) / 2) == 0
assert simplify((x * y**4 / (x**4 + y**4)).subs({x: t, y: t}) / t) == Rational(1, 2)  # f_4 non différentiable
faux(5, 4, 3, 6)
# v2 : f_k = x^k y^2/(x^2+y^2)^2.
for kv in [4, 5]:
    fk = x**kv * y**2 / (x**2 + y**2) ** 2
    assert simplify(diff(fk, x) - x ** (kv - 1) * y**2 * (kv * (x**2 + y**2) - 4 * x**2) / (x**2 + y**2) ** 3) == 0
    assert simplify(diff(fk, y) - 2 * x**kv * y * (x**2 - y**2) / (x**2 + y**2) ** 3) == 0
for kv in [1, 2]:
    assert simplify((x**kv * y**2 / (x**2 + y**2) ** 2).subs({x: t, y: t}) - t ** (kv - 2) / 4) == 0
assert simplify((x**3 * y**2 / (x**2 + y**2) ** 2).subs({x: t, y: t}) / t) == Rational(1, 4)  # f_3 non différentiable
faux(4, 3, 2, 5)
# Majorations utilisées (tirages) : |∂f_k| <= C r^(degré) pour les trois familles.
majorations = []
for kv in [3, 4]:
    ex_ = x**kv * y / (x**2 + y**2)
    majorations += [(diff(ex_, x), (kv + 2), kv - 2), (diff(ex_, y), 1, kv - 2)]
e1_ = x * y**5 / (x**4 + y**4)
majorations += [(diff(e1_, x), 16, 1), (diff(e1_, y), 4 * 5 + 16, 1)]
e2_ = x**4 * y**2 / (x**2 + y**2) ** 2
majorations += [(diff(e2_, x), 4 + 4, 1), (diff(e2_, y), 2, 1)]
majorations = [(lambdify((x, y), e), C, d) for e, C, d in majorations]
for a1_, a2_ in tirages(1500):
    if (a1_, a2_) == (0, 0):
        continue
    rr = (a1_**2 + a2_**2) ** 0.5
    for fn, C, d in majorations:
        assert abs(fn(a1_, a2_)) <= C * rr**d + 1e-12
    assert abs(a1_ * a2_**4) / (a1_**4 + a2_**4) <= 2 * rr + 1e-12  # |f_4| <= 2r (v1)

# cd-l3-differentiabilite-classe-c1-ex4 : borne supérieure 9/8 (calculée ci-dessus), leurres faux.
faux(Rational(9, 8), 1, 3, 3 * sqrt(3) / 8)

# ===========================================================================
# Carte cd-l3-differentiabilite-composee
# ===========================================================================
# exemple : polaires, G = r^2 cos 2θ.
gpol = x**2 - y**2
G = gpol.subs({x: r * cos(th), y: r * sin(th)})
assert simplify(G - r**2 * cos(2 * th)) == 0
formule = (-r * sin(th) * diff(gpol, x) + r * cos(th) * diff(gpol, y)).subs({x: r * cos(th), y: r * sin(th)})
assert simplify(formule - diff(G, th)) == 0 and simplify(diff(G, th) + 2 * r**2 * sin(2 * th)) == 0
# exemple : courbe, et déterminants de f(x,y) = (x+y, xy), g(u,v) = (u^2-2v, u).
fF = Matrix([x + y, x * y])
gG = Matrix([u**2 - 2 * v, u])
comp = gG.subs({u: fF[0], v: fF[1]}, simultaneous=True)
assert expand(comp - Matrix([x**2 + y**2, x + y])) == zeros(2, 1)
assert gG.jacobian([u, v]).subs({u: 3, v: 2}).det() == 2 and fF.jacobian([x, y]).subs({x: 2, y: 1}).det() == 1
assert comp.jacobian([x, y]).subs({x: 2, y: 1}) == Matrix([[4, 2], [1, 1]]) and Matrix([[4, 2], [1, 1]]).det() == 2
# contre-exemple : g = x^3/(x^2+y^2) le long de (t,t) : g(t,t) = t/2, alors que la formule donnerait 1.
assert simplify((x**3 / (x**2 + y**2)).subs({x: t, y: t})) == t / 2

# cd-l3-differentiabilite-composee-ex1 (et v1, v2) : dérivée le long d'une courbe.
for f_, gam, t0, rep, leurres in [
    (x**2 * y + y**3, (exp(t), 1 + 2 * t), 0, 10, (6, 2, 8)),
    (x * y - y**2, (cos(t), sin(t)), pi / 2, -1, (1, -2, 0)),
    (x * exp(y), (t**2, 3 * t), 1, 5 * exp(3), (exp(3), 6 * exp(3), 2 * exp(3))),
]:
    direct = simplify(diff(f_.subs({x: gam[0], y: gam[1]}, simultaneous=True), t).subs(t, t0))
    pt = (gam[0].subs(t, t0), gam[1].subs(t, t0))
    regle = sum(en(c, pt) * diff(gc, t).subs(t, t0) for c, gc in zip(grad(f_), gam))
    assert simplify(direct - rep) == 0 and simplify(regle - rep) == 0
    faux(rep, *leurres)
assert 2 + 4 == 6  # ex1 : leurre 6 = somme des dérivées partielles sans gamma'

# cd-l3-differentiabilite-composee-ex2 (et v1, v2) : polaires.
for f_, var_, pt, rep, leurres in [
    (x**2 - y**2, th, (2, pi / 6), -4 * sqrt(3), (4 * sqrt(3), -4, -2 * sqrt(3))),
    (x * y, th, (2, pi / 3), -2, (2, 4, -2 * sqrt(3))),
    (x**2 + 2 * y**2, r, (3, pi / 4), 9, (Rational(27, 2), 6, 3)),
]:
    Gf = f_.subs({x: r * cos(th), y: r * sin(th)}, simultaneous=True)
    if var_ == th:
        regle = (-r * sin(th) * diff(f_, x) + r * cos(th) * diff(f_, y))
    else:
        regle = (cos(th) * diff(f_, x) + sin(th) * diff(f_, y))
    regle = regle.subs({x: r * cos(th), y: r * sin(th)}, simultaneous=True)
    assert simplify(regle - diff(Gf, var_)) == 0
    assert simplify(diff(Gf, var_).subs({r: pt[0], th: pt[1]}) - rep) == 0
    faux(rep, *leurres)
assert simplify((x**2 + 2 * y**2).subs({x: 3 * cos(pi / 4), y: 3 * sin(pi / 4)})) == Rational(27, 2)  # leurre v2 : g

# cd-l3-differentiabilite-composee-ex3 (et v1, v2) : équations de transport a f_x + b f_y = 0.
def transport(a_, b_, f_):
    return simplify(a_ * diff(f_, x) + b_ * diff(f_, y))
assert transport(2, -3, (3 * x + 2 * y) ** 2) == 0
assert [factor(transport(2, -3, e)) for e in [(2 * x - 3 * y) ** 2, (2 * x + 3 * y) ** 2, (3 * x - 2 * y) ** 2]] == \
    [factor(26 * (2 * x - 3 * y)), factor(-10 * (2 * x + 3 * y)), factor(24 * (3 * x - 2 * y))]
assert transport(1, 2, sin(2 * x - y)) == 0
assert [transport(1, 2, e) for e in [sin(x - 2 * y), sin(x + 2 * y), sin(2 * x + y)]] == \
    [-3 * cos(x - 2 * y), 5 * cos(x + 2 * y), 4 * cos(2 * x + y)]
assert transport(3, 1, exp(x - 3 * y)) == 0
assert [transport(3, 1, e) for e in [exp(3 * x - y), exp(x + 3 * y), exp(3 * x + y)]] == \
    [8 * exp(3 * x - y), 6 * exp(x + 3 * y), 10 * exp(3 * x + y)]
# changements de variables : f(x,y) = g(3x+2y, y) avec g(u,v) = f((u-2v)/3, v), etc. ; dérivée en v.
fx_, fy_ = symbols("fx fy")
assert expand(Rational(-2, 3) * fx_ + fy_ + Rational(1, 3) * (2 * fx_ - 3 * fy_)) == 0
assert expand(Rational(1, 2) * fx_ + fy_ - Rational(1, 2) * (fx_ + 2 * fy_)) == 0
assert simplify((u - 2 * v) / 3 - x).subs({u: 3 * x + 2 * y, v: y}) == 0
assert simplify(((u + v) / 2).subs({u: 2 * x - y, v: y}) - x) == 0

# cd-l3-differentiabilite-composee-ex4 : différentielle de l'inversion.
Mi = Matrix([[2, 1], [1, 1]])
Hi = Matrix([[1, 0], [0, 0]])
dinv = diff((Mi + t * Hi).inv(), t).subs(t, 0)
assert dinv == -Mi.inv() * Hi * Mi.inv() == Matrix([[-1, 1], [1, -1]])
assert Mi.inv() == Matrix([[1, -1], [-1, 2]])
for l in [Matrix([[1, -1], [-1, 1]]), -Hi, -Mi.inv() ** 2 * Hi]:
    assert l != dinv
assert -Mi.inv() ** 2 * Hi == Matrix([[-2, 0], [3, 0]])
# identité de la question 1 (matrices génériques inversibles, vérifiée sur des exemples rationnels).
for _ in range(20):
    Mr = Matrix(3, 3, lambda i, j: Rational(random.randint(-5, 5), random.randint(1, 4)))
    Hr = Matrix(3, 3, lambda i, j: Rational(random.randint(-5, 5), random.randint(1, 4)))
    if Mr.det() == 0 or (Mr + Hr).det() == 0:
        continue
    gauche = (Mr + Hr).inv() - Mr.inv() + Mr.inv() * Hr * Mr.inv()
    assert gauche == (Mr + Hr).inv() * Hr * Mr.inv() * Hr * Mr.inv()

# cd-l3-differentiabilite-composee-ex5 : relation d'Euler, fonction de la question 4.
fE = x**2 + y**2
assert simplify(x * diff(fE, x) + y * diff(fE, y) - 2 * fE) == 0
for l in [sqrt(x**2 + y**2), S(1), (x**2 + y**2) ** 2]:
    assert simplify(x * diff(l, x) + y * diff(l, y) - 2 * l) != 0 or simplify(l.subs({x: 1, y: 0}) - 1) != 0
# Euler pour une fonction homogène de degré alpha générique : f(x,y) = r^alpha * psi(θ) ; ici exemples.
for fh, deg in [(x**3 / (x**2 + y**2), 1), ((x**2 + y**2) ** Rational(3, 4), Rational(3, 2)), (x * y / (x**2 + y**2), 0)]:
    assert simplify(x * diff(fh, x) + y * diff(fh, y) - deg * fh) == 0

# ===========================================================================
# Carte cd-l3-differentiabilite-accroissements-finis
# ===========================================================================
# exemple : ||∇(sin x cos y)||^2 <= 1 : 1 - (cos^2 x cos^2 y + sin^2 x sin^2 y) = cos^2 x sin^2 y + sin^2 x cos^2 y.
assert simplify(1 - (cos(x) ** 2 * cos(y) ** 2 + sin(x) ** 2 * sin(y) ** 2)
                - (cos(x) ** 2 * sin(y) ** 2 + sin(x) ** 2 * cos(y) ** 2)) == 0
# contre-exemple : fonction angle sur la couronne fendue.
R_ = sqrt(x**2 + y**2)
theta = 2 * atan(y / (x + R_))
assert simplify(diff(theta, x) + y / (x**2 + y**2)) == 0 and simplify(diff(theta, y) - x / (x**2 + y**2)) == 0
aA = (Rational(-3, 2), Rational(1, 2))
bB = (Rational(-3, 2), Rational(-1, 2))
assert 1 < aA[0] ** 2 + aA[1] ** 2 == Rational(5, 2) < 4
tA = theta.subs({x: aA[0], y: aA[1]})
tB = theta.subs({x: bB[0], y: bB[1]})
assert abs(N(tA - (pi - atan(Rational(1, 3))), 30)) < 1e-25 and abs(N(tB + (pi - atan(Rational(1, 3))), 30)) < 1e-25
diffAB = 2 * pi - 2 * atan(Rational(1, 3))
assert abs(N(tA - tB - diffAB, 30)) < 1e-25 and N(diffAB) > 5.63 and N(diffAB) > 1
# 2) θ(ρ cos φ, ρ sin φ) = φ : tirages.
for _ in range(500):
    rho = random.uniform(1.01, 1.99)
    ph = random.uniform(-3.1, 3.1)
    val = float(theta.subs({x: rho * cos(ph), y: rho * sin(ph)}))
    assert abs(val - ph) < 1e-9
# leurres de cd-l3-differentiabilite-accroissements-finis-ex4 : 2 arctan(1/3), π - arctan(1/3), 1.
faux(diffAB, 2 * atan(Rational(1, 3)), pi - atan(Rational(1, 3)), 1)
# contre-exemples : x - x^3/3 ; x + y pour la norme infinie.
fm = x - x**3 / 3
assert diff(fm, x).subs(x, 1) == 0 == diff(fm, x).subs(x, -1) and fm.subs(x, 1) - fm.subs(x, -1) == Rational(4, 3)
assert abs(2) > 1 * max(1, 1)
# affirmation 1 : f(t) = (cos 2πt, sin 2πt) : ||f'|| = 2π.
assert simplify(diff(cos(2 * pi * t), t) ** 2 + diff(sin(2 * pi * t), t) ** 2) == 4 * pi**2

# cd-l3-differentiabilite-accroissements-finis-ex1 (et v1, v2).
gr1 = grad(sin(x + 2 * y))
assert simplify(gr1[0] ** 2 + gr1[1] ** 2 - 5 * cos(x + 2 * y) ** 2) == 0
faux(sqrt(5), 3, 1, 5)
gr2 = grad(cos(3 * x - 4 * y))
assert simplify(gr2[0] ** 2 + gr2[1] ** 2 - 25 * sin(3 * x - 4 * y) ** 2) == 0
assert [c.subs({x: -pi / 6, y: 0}) for c in gr2] == [3, -4]
faux(5, 7, 1, 25)
assert limit(sin(2 * t) / t, t, 0) == 2  # v2 : f(t(1,1)) = sin 2t, norme infinie
faux(2, sqrt(2), 1, 4)

# cd-l3-differentiabilite-accroissements-finis-ex2 (et v1, v2) : nombre de c tels que f(1) - f(0) = f'(c).
def nb_c(fx_e, fy_e):
    dx_, dy_ = diff(fx_e, t), diff(fy_e, t)
    cible = (fx_e.subs(t, 1) - fx_e.subs(t, 0), fy_e.subs(t, 1) - fy_e.subs(t, 0))
    sols = [c for c in solve(dx_ - cible[0], t) if c.is_real and 0 < c < 1]
    return [c for c in sols if simplify(dy_.subs(t, c) - cible[1]) == 0]
assert nb_c(t**2, t**3) == []
assert nb_c(t**3, t) == [1 / sqrt(3)]
assert simplify(diff(cos(pi * t), t).subs(t, Rational(1, 2)) + pi) == 0 and N(-pi) != -2  # v2 : c = 1/2 donne -π
assert nb_c(sin(pi * t), cos(pi * t)) == []  # v2, lu sur la seconde coordonnée puis la première
assert 4 + 9 == 13 and 9 + 1 == 10 and 2 <= N(pi)  # bornes de la question 2

# cd-l3-differentiabilite-accroissements-finis-ex3 (et v1, v2) : normes subordonnées à ||.||_inf (sommes des lignes).
def sup_lignes(Fv, pts):
    Jf = lambdify((x, y), Fv.jacobian([x, y]))
    best = 0
    for p in pts:
        Jp = Jf(*p)
        best = max(best, max(abs(Jp[i][0]) + abs(Jp[i][1]) for i in range(2)))
    return best
grille = [(i * 0.13, j * 0.17) for i in range(-25, 26) for j in range(-25, 26)] + [(0, 0), (0, float(pi) / 2)]
for Fv, rep, leurres in [
    (Matrix([sin(x + y) / 3, cos(x - y) / 3]), Rational(2, 3), (Rational(1, 3), Rational(4, 3), sqrt(2) / 3)),
    (Matrix([sin(x + y) / 4, cos(x - y) / 4]), Rational(1, 2), (Rational(1, 4), 1, sqrt(2) / 4)),
    (Matrix([sin(x) / 2 + cos(y) / 4, sin(x - y) / 3]), Rational(3, 4), (Rational(2, 3), Rational(1, 2), Rational(17, 12))),
]:
    sup = sup_lignes(Fv, grille)
    assert abs(sup - float(rep)) < 1e-12  # atteinte en (0,0) ou (0, π/2), jamais dépassée sur la grille
    faux(rep, *leurres)
# borne théorique v2 : 1/2 + 1/4 = 3/4 et 2/3 < 3/4.
assert Rational(1, 2) + Rational(1, 4) == Rational(3, 4) > Rational(2, 3)

# cd-l3-differentiabilite-accroissements-finis-ex5 : g = (sin y / 2, cos x / 2), borne 1/2 atteinte en (0,0).
assert abs(sup_lignes(Matrix([sin(y) / 2, cos(x) / 2]), grille) - 0.5) < 1e-12
faux(Rational(1, 2), 1, Rational(1, 4), sqrt(2) / 2)
# itérations de la question 3 : convergence vers une solution de f(x) = y pour un exemple.
import math
cible_ = (1.3, -0.7)
X = (0.0, 0.0)
for _ in range(80):
    X = (cible_[0] - 0.5 * math.sin(X[1]), cible_[1] - 0.5 * math.cos(X[0]))
assert abs(X[0] + 0.5 * math.sin(X[1]) - cible_[0]) < 1e-12 and abs(X[1] + 0.5 * math.cos(X[0]) - cible_[1]) < 1e-12

# ===========================================================================
# Carte cd-l3-differentiabilite-schwarz
# ===========================================================================
fS = x**3 * y**2 + sin(x * y)
assert simplify(d2(fS, x, y) - d2(fS, y, x)) == 0
assert simplify(d2(fS, x, y) - (6 * x**2 * y + cos(x * y) - x * y * sin(x * y))) == 0
P_, Q_ = 2 * x * y + y, x**2 + x
assert diff(P_, y) == diff(Q_, x) and grad(x**2 * y + x * y) == [P_, Q_]
assert diff(y, y) != diff(-x, x)  # (y, -x) n'est pas un gradient
# contre-exemple : xy(x^2-y^2)/(x^2+y^2).
fX = x * y * (x**2 - y**2) / (x**2 + y**2)
y0 = symbols("y0", nonzero=True, real=True)
assert limit(fX.subs({x: h, y: y0}) / h, h, 0) == -y0
assert limit(fX.subs({x: y0, y: h}) / h, h, 0) == y0
assert fX.subs({x: h, y: 0}) == 0 and fX.subs({x: 0, y: h}) == 0
# dérivées secondes homogènes de degré 0 (f homogène de degré 2), dérivées premières de degré 1.
lam = symbols("lam", positive=True)
assert simplify(fX.subs({x: lam * x, y: lam * y}, simultaneous=True) - lam**2 * fX) == 0
for p_, q_ in [(x, x), (x, y), (y, x), (y, y)]:
    dd = simplify(d2(fX, p_, q_))
    assert simplify(dd.subs({x: lam * x, y: lam * y}, simultaneous=True) - dd) == 0
# ∂f/∂x(x, 0) = 0 et ∂f/∂y(0, y) = 0, d'où ∂²f/∂x² (0,0) = ∂²f/∂y² (0,0) = 0.
assert simplify(diff(fX, x).subs(y, 0)) == 0 and simplify(diff(fX, y).subs(x, 0)) == 0

# cd-l3-differentiabilite-schwarz-ex1 (et v1, v2).
for f_, a_, rep, leurres in [
    (x**3 * y**2 - 2 * x * y**3, (1, 2), -12, (12, -24, 24)),
    (x**2 * y**3 + 3 * x**3 * y, (1, 1), 15, (6, 9, -15)),
    (x * exp(y) + x**3 * y**2, (2, 0), 1, (0, 2, 13)),
]:
    assert en(d2(f_, x, y), a_) == en(d2(f_, y, x), a_) == rep
    faux(rep, *leurres)
assert simplify(d2(x**3 * y**2 - 2 * x * y**3, x, y) - (6 * x**2 * y - 6 * y**2)) == 0

# cd-l3-differentiabilite-schwarz-ex2 (et v1, v2).
for f_, rep_yx, rep_xy in [
    (x * y * (x**2 - y**2) / (x**2 + y**2), -1, 1),
    (x * y * (x**2 - 2 * y**2) / (x**2 + y**2), -2, 1),
    (x * y * (2 * x**2 - y**2) / (x**2 + y**2), -1, 2),
]:
    dx0 = limit(f_.subs({x: h, y: y0}) / h, h, 0)  # ∂f/∂x(0, y0), y0 != 0
    dy0 = limit(f_.subs({x: y0, y: h}) / h, h, 0)  # ∂f/∂y(x0, 0), x0 != 0
    assert f_.subs({x: h, y: 0}) == 0 and f_.subs({x: 0, y: h}) == 0
    assert simplify(dx0 / y0) == rep_yx and simplify(dy0 / y0) == rep_xy
faux(-1, 1, 0, 2)
faux(-2, 1, 0, 2)
faux(2, -1, 0, 1)

# cd-l3-differentiabilite-schwarz-ex3 (et v1, v2) : valeur de a et potentiel.
for P_, Q_, rep, pot, leurres in [
    (2 * x * y + y, x**2 + al * x, 1, x**2 * y + x * y, (0, 2, -1)),
    (3 * x**2 * y**2 + al * y, 2 * x**3 * y + 2 * x, 2, x**3 * y**2 + 2 * x * y, (0, 1, 3)),
    (al * y**2 + 2 * x, 6 * x * y, 3, 3 * x * y**2 + x**2, (6, 1, 2)),
]:
    assert solve(expand(diff(P_, y) - diff(Q_, x)).subs({x: 1, y: 1}), al) == [rep]
    assert simplify(diff(P_, y) - diff(Q_, x)).subs(al, rep) == 0
    assert grad(pot) == [P_.subs(al, rep), Q_.subs(al, rep)]
    faux(rep, *leurres)

# cd-l3-differentiabilite-schwarz-ex4 : équation des ondes.
cpos = symbols("c", positive=True)
X_, T_ = symbols("X T")
# vérification sur une fonction test polynomiale arbitraire F(X,T) = X^3 T + X T^2 + T^4.
Ft = X_**3 * T_ + X_ * T_**2 + T_**4
gt = Ft.subs({X_: (u + v) / 2, T_: (u - v) / (2 * cpos)}, simultaneous=True)
lhs = diff(gt, u, v)
rhs = (Rational(1, 4) * (diff(Ft, X_, 2) - diff(Ft, T_, 2) / cpos**2)).subs({X_: (u + v) / 2, T_: (u - v) / (2 * cpos)}, simultaneous=True)
assert simplify(lhs - rhs) == 0
def onde(f_):
    return simplify(diff(f_, t, 2) - diff(f_, x, 2))
assert onde(x**2 + t**2) == 0
assert [onde(e) for e in [x**2 - t**2, x * t**2, exp(x) * cos(t)]] == [-4, 2 * x, -2 * exp(x) * cos(t)]
assert expand(Rational(1, 2) * (x + t) ** 2 + Rational(1, 2) * (x - t) ** 2) == x**2 + t**2
assert onde(sin(x) * sin(t)) == 0  # (raison pour laquelle cette fonction n'est pas un leurre)

# cd-l3-differentiabilite-schwarz-ex5 : laplacien en polaires et fonctions radiales harmoniques.
for ftest in [x**3 * y - y**2 * x, exp(x) * sin(y), x**4 + 3 * x * y**2]:
    gt = ftest.subs({x: r * cos(th), y: r * sin(th)}, simultaneous=True)
    lap = (diff(ftest, x, 2) + diff(ftest, y, 2)).subs({x: r * cos(th), y: r * sin(th)}, simultaneous=True)
    assert simplify(diff(gt, r, 2) + diff(gt, r) / r + diff(gt, th, 2) / r**2 - lap) == 0
aa, bb = symbols("a b", real=True)
def lap_radial(ph):
    return simplify(diff(ph, r, 2) + diff(ph, r) / r)
assert lap_radial(aa * ln(r) + bb) == 0
for l in [aa / r + bb, aa * r**2 + bb, aa / r**2 + bb]:
    assert lap_radial(l) != 0
assert simplify(diff(r * diff(Function("p")(r), r), r) / r - (diff(Function("p")(r), r, 2) + diff(Function("p")(r), r) / r)) == 0

# ===========================================================================
# Carte cd-l3-differentiabilite-gradient
# ===========================================================================
f7 = x**2 + 3 * x * y
assert [en(c, (1, 1)) for c in grad(f7)] == [5, 3] and sqrt(5**2 + 3**2) == sqrt(34)
f8 = x**2 + x * y + y**3
assert en(f8, (1, 1)) == 3 and [en(c, (1, 1)) for c in grad(f8)] == [3, 4] and 3 * 4 + 4 * (-3) == 0
assert 3 * 1 + 4 * 1 == 7
A2 = Matrix(2, 2, symbols("p0:4"))
X2 = Matrix([x, y])
qf = expand((X2.T * A2 * X2)[0])
assert expand(Matrix(grad(qf)) - (A2 + A2.T) * X2) == zeros(2, 1)
assert grad(x**2 + y**2) == [2 * x, 2 * y]
# contre-exemples : bases non orthonormées ; produit scalaire 2u1v1 + u2v2.
e1, e2 = Matrix([1, 0]), Matrix([1, 1])
fy = lambda P: P[1]
assert fy(s * e1 + t * e2) == t  # dérivées partielles 0 et 1 dans la base (e1, e2)
assert 0 * e1 + 1 * e2 == Matrix([1, 1]) != Matrix([0, 1])
e1b, e2b = Matrix([2, 0]), Matrix([0, 1])
assert (s * e1b + t * e2b)[0] == 2 * s and 2 * e1b == Matrix([4, 0]) != Matrix([1, 0])
Sm = Matrix([[2, 0], [0, 1]])
assert (Matrix([Rational(1, 2), 0]).T * Sm)[0, 0] == 1 and (Matrix([Rational(1, 2), 0]).T * Sm)[0, 1] == 0
assert (x * y).subs({x: t, y: t}) == t**2 and (x * y).subs({x: t, y: -t}) == -t**2
# trou 2 (leurre) : base de vecteurs unitaires non orthogonaux, e1 = (1,0), e2 = (1,1)/sqrt 2, f = y.
e2u = Matrix([1, 1]) / sqrt(2)
assert simplify((s * e1 + t * e2u)[1] - t / sqrt(2)) == 0
assert simplify(0 * e1 + e2u / sqrt(2) - Matrix([Rational(1, 2), Rational(1, 2)])) == zeros(2, 1)
# affirmation : grad ||x|| = x/||x|| ; pour x = (2,0,0), (1,0,0) != x.
nrm = sqrt(x**2 + y**2 + z**2)
assert [simplify(en(c, (2, 0, 0))) for c in grad(nrm, (x, y, z))] == [1, 0, 0]

# cd-l3-differentiabilite-gradient-ex1 (et v1, v2).
for f_, a_, rep, leurres in [
    (x**2 * y + y * z**3, (1, 2, 1), (4, 2, 6), [(4, 1, 6), (2, 2, 6), (6, 2, 4)]),
    (x * y * z + x**2, (1, -1, 2), (0, 2, -1), [(-2, 2, -1), (2, 2, -1), (0, -1, 2)]),
    (exp(x * y) + z**2, (0, 3, 1), (3, 0, 2), [(0, 3, 2), (3, 0, 1), (1, 0, 2)]),
]:
    assert tuple(en(c, a_) for c in grad(f_, (x, y, z))) == rep
    for l in leurres:
        assert l != rep

# cd-l3-differentiabilite-gradient-ex2 (et v1, v2) : pente maximale = ||∇f(a)||.
for f_, a_, rep, leurres in [
    (x**2 + 3 * x * y, (1, 1), sqrt(34), (8, 34, 5)),
    (x * y**2, (2, 1), sqrt(17), (5, 17, 4)),
    (ln(x**2 + y**2), (1, 2), 2 * sqrt(5) / 5, (Rational(6, 5), Rational(4, 5), sqrt(5))),
]:
    gv = [en(c, a_) for c in grad(f_)]
    assert simplify(sqrt(gv[0] ** 2 + gv[1] ** 2) - rep) == 0
    faux(rep, *leurres)
    # maximum de <∇f, (cos φ, sin φ)> sur le cercle unité, par échantillonnage.
    mx = max(float(gv[0] * math.cos(p) + gv[1] * math.sin(p)) for p in [i * 2 * math.pi / 20000 for i in range(20000)])
    assert abs(mx - float(rep)) < 1e-6
assert [en(c, (1, 2)) for c in grad(ln(x**2 + y**2))] == [Rational(2, 5), Rational(4, 5)]

# cd-l3-differentiabilite-gradient-ex3 (et v1, v2) : tangente à une ligne de niveau.
for f_, a_, niveau, droite, leurres in [
    (x**2 + x * y + y**3, (1, 1), 3, 3 * x + 4 * y - 7, [4 * x - 3 * y - 1, 3 * x + 4 * y - 3, 2 * x + y - 3]),
    (x**3 + y**3 - x * y, (1, 1), 1, x + y - 2, [x - y, 2 * x + 2 * y - 1, x + y - 1]),
    (x * exp(y) + y, (2, 0), 2, x + 3 * y - 2, [3 * x - y - 6, x + 3 * y - 6, 3 * x + y - 6]),
]:
    assert en(f_, a_) == niveau
    gv = [en(c, a_) for c in grad(f_)]
    # la droite passe par a et a pour vecteur normal ∇f(a) (coefficients proportionnels).
    assert en(droite, a_) == 0
    cx, cy = diff(droite, x), diff(droite, y)
    assert simplify(cx * gv[1] - cy * gv[0]) == 0
    for l in leurres:
        lx, ly = diff(l, x), diff(l, y)
        assert en(l, a_) != 0 or simplify(lx * gv[1] - ly * gv[0]) != 0

# cd-l3-differentiabilite-gradient-ex4 (et v1, v2) : gradient pour <u,v>_S = u^T S v.
for Sm_, f_, a_, rep, leurres in [
    (Matrix([[2, 1], [1, 1]]), x**2 + y, (1, 0), (1, 0), [(2, 1), (5, 3), (1, 1)]),
    (Matrix([[1, 0], [0, 4]]), x * y, (2, 1), (1, Rational(1, 2)), [(1, 2), (1, 8), (2, 1)]),
    (Matrix([[3, 1], [1, 1]]), x + y**2, (0, 1), (Rational(-1, 2), Rational(5, 2)), [(1, 2), (5, 3), (Rational(1, 2), Rational(-5, 2))]),
]:
    assert Sm_.is_symmetric() and all(ev > 0 for ev in Sm_.eigenvals())
    gcan = Matrix([en(c, a_) for c in grad(f_)])
    gS = Sm_.inv() * gcan
    assert tuple(gS) == rep and Sm_ * gS == gcan
    for l in leurres:
        assert l != rep
    assert tuple(Sm_ * gcan) in [tuple(l) for l in leurres] or tuple(gcan) in [tuple(l) for l in leurres]

# cd-l3-differentiabilite-gradient-ex5 : méthode du gradient, A = diag(1,4) : convergence ssi 0 < ρ < 1/2.
rho = symbols("rho", positive=True)
Aq = Matrix([[1, 0], [0, 4]])
bq = Matrix([1, 1])
xq = Matrix([x, y])
fq = (Rational(1, 2) * xq.T * Aq * xq - bq.T * xq)[0]
assert Matrix(grad(fq)) == Aq * xq - bq
def converge(rv):
    return all(abs(1 - rv * lv) < 1 for lv in [1, 4])
assert converge(Rational(49, 100)) and not converge(Rational(1, 2)) and not converge(Rational(3, 4))
assert not converge(Rational(3, 2)) and converge(Rational(1, 5))
# leurres : 0<ρ<2, 0<ρ<1/4, 0<ρ<1 ne décrivent pas l'ensemble ]0, 1/2[ (ρ = 3/4 ou ρ = 0,3 les réfutent).
assert not converge(Rational(3, 4)) and converge(Rational(3, 10))

# ===========================================================================
# Carte cd-l3-differentiabilite-extremums
# ===========================================================================
def nature(f_, a_):
    rr_, ss_, tt_ = en(diff(f_, x, 2), a_), en(diff(f_, x, y), a_), en(diff(f_, y, 2), a_)
    assert [en(c, a_) for c in grad(f_)] == [0, 0], "point non critique"
    dlt = rr_ * tt_ - ss_**2
    if dlt > 0:
        return "min" if rr_ > 0 else "max"
    if dlt < 0:
        return "col"
    return "?"
fP = x**3 + y**3 - 3 * x * y
assert set(map(tuple, [[s_[x], s_[y]] for s_ in solve(grad(fP), [x, y], dict=True) if s_[x].is_real and s_[y].is_real])) == {(0, 0), (1, 1)}
assert nature(fP, (0, 0)) == "col" and nature(fP, (1, 1)) == "min" and en(fP, (1, 1)) == -1
# Identités de la démonstration.
rr, ss, tt = symbols("r s t", real=True)
QQ = rr * h**2 + 2 * ss * h * k + tt * k**2
assert expand(rr * QQ - ((rr * h + ss * k) ** 2 + (rr * tt - ss**2) * k**2)) == 0
assert expand(QQ.subs({h: -ss, k: rr}) - rr * (rr * tt - ss**2)) == 0
# contre-exemples : r, s, t des cas rt - s^2 = 0 ; x^2 + 3xy + y^2 ; x + x^2 + y^2 non critique.
for f_ in [x**2 + y**4, x**2 - y**4, x**2 + y**3]:
    assert nature(f_, (0, 0)) == "?" and [en(diff(f_, x, 2), (0, 0)), en(diff(f_, x, y), (0, 0)), en(diff(f_, y, 2), (0, 0))] == [2, 0, 0]
assert nature(x**2 + 3 * x * y + y**2, (0, 0)) == "col" and expand((x**2 + 3 * x * y + y**2).subs({x: t, y: -t})) == -t**2
assert [en(c, (0, 0)) for c in grad(x + x**2 + y**2)] == [1, 0]
assert (x + x**2 + y**2).subs({x: Rational(-1, 2), y: 0}) < 0
# trous (leurres) : -x^2 - y^2 (rt - s^2 > 0, r < 0 : maximum) ; y^2 + x^4 (r = 0, minimum).
assert nature(-x**2 - y**2, (0, 0)) == "max"
assert en(diff(y**2 + x**4, x, 2), (0, 0)) == 0
# affirmation : x^4 + y^4 en (0,0) : r = s = t = 0, minimum strict.
assert nature(x**4 + y**4, (0, 0)) == "?"

# cd-l3-differentiabilite-extremums-ex1 (et v1, v2) : bonnes réponses 0 (min), 2 (col), 1 (max).
assert nature(x**2 + x * y + y**2 + x**3, (0, 0)) == "min"
assert nature(x**2 + 4 * x * y + y**2 + y**3, (0, 0)) == "col"
assert expand((x**2 + 4 * x * y + y**2 + y**3).subs({x: t, y: -t})) == -2 * t**2 - t**3
assert nature(-2 * x**2 + 2 * x * y - y**2 + x**4, (0, 0)) == "max"

# cd-l3-differentiabilite-extremums-ex2 (et v1, v2).
def points_critiques(f_):
    sols = solve(grad(f_), [x, y], dict=True)
    return {(s_[x], s_[y]) for s_ in sols if s_[x].is_real and s_[y].is_real}
fQ = x**3 - 3 * x + y**2 - 2 * y
assert points_critiques(fQ) == {(1, 1), (-1, 1)}
assert nature(fQ, (1, 1)) == "min" and nature(fQ, (-1, 1)) == "col" and en(fQ, (1, 1)) == -3
fR = x**4 + y**4 - 4 * x * y
assert points_critiques(fR) == {(0, 0), (1, 1), (-1, -1)}
assert nature(fR, (0, 0)) == "col" and nature(fR, (1, 1)) == "min" == nature(fR, (-1, -1)) and en(fR, (1, 1)) == -2
fT = x**2 * y + x * y**2 - 3 * x * y
assert points_critiques(fT) == {(0, 0), (0, 3), (3, 0), (1, 1)}
assert [nature(fT, p) for p in [(0, 0), (0, 3), (3, 0)]] == ["col"] * 3 and nature(fT, (1, 1)) == "min"
assert en(fT, (1, 1)) == -1

# cd-l3-differentiabilite-extremums-ex3 : f_1 et f_4 seulement.
f4b = -x**4 - y**4
assert nature(f4b, (0, 0)) == "?" and (x**2 - y**4).subs({x: 0, y: t}) == -t**4

# cd-l3-differentiabilite-extremums-ex4 (et v1, v2) : extremums sur un compact (bord + intérieur + grille).
def extremums(f_, sommets, bords, interieur):
    vals = []
    for p in points_critiques(f_):
        if interieur(p):
            vals.append(en(f_, p))
    for (par, lo, hi) in bords:  # par : (x(s), y(s)) pour s dans [lo, hi]
        g_ = f_.subs({x: par[0], y: par[1]}, simultaneous=True)
        cand = [lo, hi] + [c for c in solve(diff(g_, s), s) if c.is_real and lo <= c <= hi]
        vals += [simplify(g_.subs(s, c)) for c in cand]
    return min(vals), max(vals), vals
fT1 = x**2 - x * y + y**2 - x
mn, mx_, vals = extremums(fT1, None, [((s, 0), 0, 2), ((0, s), 0, 2), ((s, 2 - s), 0, 2)],
                          lambda p: p[0] > 0 and p[1] > 0 and p[0] + p[1] < 2)
assert (mn, mx_) == (Rational(-1, 3), 4)
assert set(vals) >= {Rational(-1, 4), Rational(-1, 12), 0, 2, 4, Rational(-1, 3)}
assert expand(fT1.subs(y, 2 - x)) == 3 * x**2 - 7 * x + 4
fT2 = x**2 + y**2 - x * y - x - y
mn, mx_, vals = extremums(fT2, None, [((s, 0), 0, 3), ((0, s), 0, 3), ((s, 3 - s), 0, 3)],
                          lambda p: p[0] > 0 and p[1] > 0 and p[0] + p[1] < 3)
assert (mn, mx_) == (-1, 6) and Rational(-3, 4) in vals
assert expand(fT2.subs(y, 3 - x)) == 3 * x**2 - 9 * x + 6
fK = x**2 - 2 * x * y + 2 * y**2 - 2 * y
mn, mx_, vals = extremums(fK, None, [((0, s), 0, 2), ((2, s), 0, 2), ((s, 0), 0, 2), ((s, 2), 0, 2)],
                          lambda p: 0 < p[0] < 2 and 0 < p[1] < 2)
assert (mn, mx_) == (-1, 4) and Rational(-1, 2) in vals and en(fK, (0, 2)) == 4 == en(fK, (2, 0))
for (f_, test, mn_, mx2) in [(fT1, lambda a, b: a + b <= 2, Rational(-1, 3), 4), (fT2, lambda a, b: a + b <= 3, -1, 6),
                             (fK, lambda a, b: a <= 2 and b <= 2, -1, 4)]:
    for i in range(0, 61):
        for j in range(0, 61):
            a, b = i * 0.05, j * 0.05
            if test(a, b):
                val = float(f_.subs({x: a, y: b}))
                assert mn_ - 1e-12 <= val <= mx2 + 1e-12
# leurres : (-1/4, 4), (-1/3, 2), (-1/12, 4) ; (-3/4, 6), (-1, 3), (-1/4, 6) ; (-1/2, 4), (-1, 8), (0, 4).
for bon, leurres in [((Rational(-1, 3), 4), [(Rational(-1, 4), 4), (Rational(-1, 3), 2), (Rational(-1, 12), 4)]),
                     ((-1, 6), [(Rational(-3, 4), 6), (-1, 3), (Rational(-1, 4), 6)]),
                     ((-1, 4), [(Rational(-1, 2), 4), (-1, 8), (0, 4)])]:
    for l in leurres:
        assert l != bon

# cd-l3-differentiabilite-extremums-ex5 : exemple de Peano.
fPe = (y - x**2) * (y - 2 * x**2)
assert expand(fPe) == y**2 - 3 * x**2 * y + 2 * x**4
assert [en(c, (0, 0)) for c in grad(fPe)] == [0, 0]
assert [en(diff(fPe, x, 2), (0, 0)), en(diff(fPe, x, y), (0, 0)), en(diff(fPe, y, 2), (0, 0))] == [0, 0, 2]
mm = symbols("m", real=True)
assert expand(fPe.subs(y, mm * x)) == mm**2 * x**2 - 3 * mm * x**3 + 2 * x**4
assert fPe.subs({x: 0}) == y**2
assert simplify(fPe.subs(y, Rational(3, 2) * x**2)) == -x**4 / 4
faux(-x**4 / 4, x**4 / 4, -x**4 / 2, 0)

# ===========================================================================
# QCM cd-l3-differentiabilite
# ===========================================================================
# q1 : x^2 y/(x^2+y^2) : continue (|f| <= |y|), D_(1,1) = 1/2, dérivées partielles nulles.
fq1 = x**2 * y / (x**2 + y**2)
assert dir_deriv(fq1, (0, 0), (1, 1)) == Rational(1, 2) and dir_deriv(fq1, (0, 0), (1, 0)) == 0 == dir_deriv(fq1, (0, 0), (0, 1))
# q2 : déterminant jacobien des polaires = r (vérifié plus haut).
# q3 : x^3/(x^2+y^2) : dérivées partielles bornées, D_(1,1) = 1/2 != 1 (vérifié plus haut).
# q4 : g(t) = f(t^2, 3t), g'(1) = 2 f_x(1,3) + 3 f_y(1,3).
X1, Y1 = symbols("X1 Y1")
fpoly = X1**3 * Y1 + 2 * X1 * Y1**2
gpoly = fpoly.subs({X1: t**2, Y1: 3 * t})
assert diff(gpoly, t).subs(t, 1) == 2 * diff(fpoly, X1).subs({X1: 1, Y1: 3}) + 3 * diff(fpoly, Y1).subs({X1: 1, Y1: 3})
assert diff(gpoly, t).subs(t, 1) != diff(fpoly, X1).subs({X1: 1, Y1: 3}) + diff(fpoly, Y1).subs({X1: 1, Y1: 3})
# q5 : la ligne brisée p -> m -> q a une longueur qui tend vers ||p - q|| quand m tend vers le milieu.
pq_p, pq_q = (-1.0, 0.0), (2.0, 0.0)
for eps in [1e-1, 1e-3, 1e-6]:
    mpt = (0.5, eps)
    long = math.dist(pq_p, mpt) + math.dist(mpt, pq_q)
    assert long - 3.0 < eps
# q6 : e^{xy} sin x est C^2 et ses dérivées croisées sont égales ; |x| y n'a pas de dérivée en x en (0, y), y != 0.
fq6 = exp(x * y) * sin(x)
assert simplify(d2(fq6, x, y) - d2(fq6, y, x)) == 0
assert limit((Abs(h) * 1 - 0) / h, h, 0, "+") == 1 and limit((Abs(h) * 1 - 0) / h, h, 0, "-") == -1
# q7 : grad <Ax,x> = (A + A^T) x, différent de 2Ax pour A = [[0,1],[0,0]].
A7 = Matrix([[0, 1], [0, 0]])
q7 = expand((X2.T * A7 * X2)[0])
assert q7 == x * y and Matrix(grad(q7)) == Matrix([y, x]) != 2 * A7 * X2
# q8 : r = 1, s = 2, t = 5 : rt - s^2 = 1 > 0 et r > 0 : minimum local strict.
assert 1 * 5 - 2**2 == 1 > 0

print("cd_l3_differentiabilite : toutes les vérifications passent.")
