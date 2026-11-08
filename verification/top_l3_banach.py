"""Vérifications calculatoires du chapitre « top-l3-banach » (topologie et calcul différentiel, L3) : cartes, exercices
et QCM de content/*/topologie-calcul-diff/l3-top-l3-banach.json et content/qcm/top-l3-banach.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question. Calculs
exacts en SymPy et en fractions, contrôles numériques avec mpmath (exponentielles de matrices, normes subordonnées).
Lancer : python3 verification/top_l3_banach.py.
"""

from fractions import Fraction as Fr
import random

import mpmath as mp
import sympy as sp

t, x, s, n_ = sp.symbols("t x s n", positive=True)
X = sp.symbols("X")
n = sp.symbols("n", integer=True, positive=True)
k = sp.symbols("k", integer=True, nonnegative=True)
mp.mp.dps = 30


def egal(a, b):
    return sp.simplify(sp.nsimplify(a) - sp.nsimplify(b)) == 0


def differents(bonne, leurres, valeurs=None):
    """Les leurres sont distincts de la bonne réponse (et entre eux), pour chaque valeur de n proposée."""
    for v in (valeurs or [None]):
        sub = (lambda e: e) if v is None else (lambda e: sp.sympify(e).subs(n, v))
        b = sp.nsimplify(sub(bonne))
        vals = [sp.nsimplify(sub(l)) for l in leurres]
        assert all(sp.simplify(b - l) != 0 for l in vals), (bonne, leurres, v)


# ===========================================================================
# Carte top-l3-banach-banach
# ===========================================================================
# trous : leurre 1 (terme général -> 0, terme de norme décroissante, sommes partielles bornées... dans R complet) :
# la série harmonique diverge ; leurre « sommes partielles bornées » : sum (-1)^n diverge
assert sp.summation(1 / (k + 1), (k, 0, sp.oo)) == sp.oo
assert sp.limit(sp.Rational(1, 1) / (n + 1), n, sp.oo) == 0
# (-1)^n : sommes partielles 1, 0, 1, 0... bornées mais non convergentes
sommes = [sum((-1) ** j for j in range(m + 1)) for m in range(10)]
assert set(sommes) == {0, 1}
# trou 2 : x0 = 1, x1 = -1 : |somme| = 0 < 2 ; un seul terme non nul : égalité (donc « strictement inférieur » faux)
assert abs(1 + (-1)) == 0 < 2 and abs(3) == 3

# contre-exemple : R[X] et la norme des coefficients, sommes partielles de X^n/2^n ; coefficient de X^(d+1) = 2^-(d+1)
for d in range(6):
    for N in range(d + 1, d + 6):
        SN = sum(X ** j / sp.Integer(2) ** j for j in range(N + 1))
        assert sp.Poly(SN, X).coeff_monomial(X ** (d + 1)) == sp.Rational(1, 2 ** (d + 1))

# exemple : ||cos(nt)/n^2||_inf = 1/n^2 (atteint en t = 0) et la série des 1/n^2 converge
assert sp.summation(1 / k ** 2, (k, 1, sp.oo)) == sp.pi ** 2 / 6

# banach-ex1 (et v1, v2) : sum ||f_n|| et ||S||_inf
for r_, S_expr, sup_attendu, somme_attendue, leurres in [
    (sp.Rational(1, 3), 3 / (3 + t), 1, sp.Rational(3, 2),
     [(sp.Rational(3, 4), sp.Rational(3, 2)), (sp.Rational(3, 2), sp.Rational(3, 2)), (1, sp.Rational(3, 4))]),
    (sp.Rational(1, 2), 2 / (2 + t), 1, 2,
     [(sp.Rational(2, 3), 2), (2, 2), (1, sp.Rational(2, 3))]),
    (sp.Rational(1, 3), 3 / (4 - 2 * t), sp.Rational(3, 2), sp.Rational(3, 2),
     [(1, sp.Rational(3, 2)), (sp.Rational(3, 4), sp.Rational(3, 2)), (sp.Rational(3, 2), sp.Rational(3, 4))]),
]:
    assert sp.summation(r_ ** k, (k, 0, sp.oo)) == somme_attendue
    vals = [S_expr.subs(t, v) for v in (0, 1)]
    assert sp.diff(S_expr, t).subs(t, sp.Rational(1, 2)) != 0  # S monotone : extrêmes aux bords
    assert max(vals) == sup_attendu
    for l in leurres:
        assert l != (sup_attendu, somme_attendue)
# vérification des sommes géométriques en un point
assert sp.simplify(sp.summation((-t / 3) ** k, (k, 0, sp.oo)).rewrite(sp.Piecewise).args[0][0] - 3 / (3 + t)) == 0
assert sp.simplify(sp.summation((-t / 2) ** k, (k, 0, sp.oo)).rewrite(sp.Piecewise).args[0][0] - 2 / (2 + t)) == 0
u_ = sp.symbols("u")
assert sp.simplify(1 / (1 - (2 * t - 1) / 3) - 3 / (4 - 2 * t)) == 0

# banach-ex2 : ||P_q - P_p|| = 1/(p+2)^2 (plus grand coefficient de la différence)
for p in range(8):
    for q in range(p + 1, p + 6):
        diff = [sp.Rational(1, (j + 1) ** 2) for j in range(p + 1, q + 1)]
        assert max(diff) == sp.Rational(1, (p + 2) ** 2)

# banach-ex3 (v1, v2) : int_0^{1/n^a} (t^{-1/a'} - n) dt
for alpha, borne, attendu, leurres in [
    (sp.Rational(1, 2), 1 / n ** 2, 1 / n, [2 / n, 1 / n ** 2, 1 / (2 * n)]),
    (sp.Rational(1, 3), 1 / n ** 3, 1 / (2 * n ** 2), [sp.Rational(3, 2) / n ** 2, 1 / n ** 2, 1 / (2 * n ** 3)]),
    (sp.Rational(1, 4), 1 / n ** 4, 1 / (3 * n ** 3), [sp.Rational(4, 3) / n ** 3, 1 / n ** 3, 1 / (3 * n ** 4)]),
]:
    I_ = sp.integrate(t ** (-alpha) - n, (t, 0, borne))
    assert sp.simplify(I_ - attendu) == 0
    assert sp.simplify(borne ** (-alpha) - n) == 0  # point de raccord
    differents(attendu, leurres, valeurs=range(2, 30))
# ||f_m - f_n||_1 = 1/n - 1/m (main) et série télescopique de somme 1/2
assert sp.summation(1 / n - 1 / (n + 1), (n, 2, sp.oo)) == sp.Rational(1, 2)

# banach-ex4 : ||x^(p) - x||_inf = 1/(p+2) ; x = (1/(n+1)) n'est pas dans l^1
for p in range(10):
    assert max(Fr(1, j + 1) for j in range(p + 1, p + 50)) == Fr(1, p + 2)

# banach-ex5 : L(g_n) = sqrt(n)/2 = g_n'(0) ; leurres distincts pour n >= 2 ; |g_n - sqrt| <= 1/sqrt(n)
gn = sp.sqrt(x + 1 / n) - sp.sqrt(1 / n)
assert sp.simplify(sp.diff(gn, x).subs(x, 0) - sp.sqrt(n) / 2) == 0
assert sp.simplify(sp.diff(gn, x, 2)) != 0  # g_n' strictement décroissante (dérivée seconde < 0)
differents(sp.sqrt(n) / 2, [sp.sqrt(n), n / 2, sp.Rational(1, 2)], valeurs=range(2, 40))
for nn in (2, 5, 50):
    for xv in [i / 200 for i in range(201)]:
        assert abs(float(mp.sqrt(xv + 1 / nn) - mp.sqrt(1 / nn) - mp.sqrt(xv))) <= 1 / mp.sqrt(nn) + 1e-15

# ===========================================================================
# Carte top-l3-banach-uniforme
# ===========================================================================
# contre-exemple : sum x^n (1 - x) = 1 sur [0,1[ ; ||f_n|| = n^n/(n+1)^(n+1) >= 1/(e(n+1))
assert sp.simplify(sp.summation(x ** k * (1 - x), (k, 0, sp.oo)).rewrite(sp.Piecewise).args[0][0]) == 1
for m in range(1, 60):
    fx = x ** m * (1 - x)
    xm = sp.solve(sp.diff(fx, x), x)
    xm = [v for v in xm if 0 < v < 1][0]
    assert xm == sp.Rational(m, m + 1)
    sup = sp.Rational(m ** m, (m + 1) ** (m + 1))
    assert fx.subs(x, xm) == sup
    assert sup >= 1 / (sp.E * (m + 1))

# uniforme-ex1 (v1, v2) : sup ||S(t)|| = somme des ||f_n|| (atteint en 0) ; |S| = 1/|q^-1 - e^{it}| pour q = 1/2, 1/3
for q_, sup_att, leurres in [(sp.Rational(1, 2), 1, [sp.Rational(1, 2), sp.Rational(1, 3), 2]),
                             (sp.Rational(1, 3), sp.Rational(1, 2), [sp.Rational(1, 4), sp.Rational(1, 3), 1])]:
    assert sp.summation(q_ ** k, (k, 1, sp.oo)) == sup_att
    inv = 1 / q_
    vals = [1 / abs(inv - mp.e ** (1j * th)) for th in [2 * mp.pi * i / 400 for i in range(400)]]
    assert abs(max(vals) - sup_att) < 1e-20 and abs(min(vals) - 1 / (inv + 1)) < 1e-20
    assert all(l != sup_att for l in leurres)
assert sp.summation(1 / k ** 2, (k, 1, sp.oo)) == sp.pi ** 2 / 6
differents(sp.pi ** 2 / 6, [1, sp.pi ** 2 / 12, sp.pi / 6])

# uniforme-ex2 (v1, v2) : max de x^a e^{-nx} sur [0, +oo[ atteint en a/n
for a, attendu, leurres in [
    (2, 4 / (sp.E ** 2 * n ** 2), [2 / (sp.E * n ** 2), 4 / (sp.E * n ** 2), 1 / (sp.E ** 2 * n ** 2)]),
    (3, 27 / (sp.E ** 3 * n ** 3), [27 / (sp.E * n ** 3), 9 / (sp.E ** 3 * n ** 3), 1 / (sp.E ** 3 * n ** 3)]),
    (1, 1 / (sp.E * n), [1 / n, 2 / (sp.E * n), 1 / (sp.E ** 2 * n)]),
]:
    f = x ** a * sp.exp(-n * x)
    crit = sp.solve(sp.diff(f, x), x)
    assert sp.Rational(a, 1) / n in crit
    assert sp.simplify(f.subs(x, a / n) - attendu) == 0
    differents(attendu, leurres, valeurs=range(1, 20))
    # somme pour x > 0 : x^a / (e^x - 1)
    somme = sp.summation(x ** a * sp.exp(-x) ** k, (k, 1, sp.oo)).rewrite(sp.Piecewise).args[0][0]
    assert sp.simplify(somme - x ** a / (sp.exp(x) - 1)) == 0
assert sp.limit(x ** 2 / (sp.exp(x) - 1), x, 0) == 0 and sp.limit(x / (sp.exp(x) - 1), x, 0) == 1
assert sp.summation(4 / (sp.E ** 2 * n ** 2), (n, 1, sp.oo)).is_finite
assert sp.summation(1 / (sp.E * n), (n, 1, sp.oo)) == sp.oo

# uniforme-ex3 (v1, v2) : (Kf)(x) = lam * x * m(f), m(f) = int_0^1 t f(t) dt ; solution f = 1 + c x
c = sp.symbols("c")
m_iota = sp.integrate(t * t, (t, 0, 1))
m_un = sp.integrate(t, (t, 0, 1))
assert m_iota == sp.Rational(1, 3) and m_un == sp.Rational(1, 2)
for lam, c_att, leurres in [(1, sp.Rational(3, 4), [sp.Rational(1, 2), sp.Rational(3, 2), sp.Rational(2, 3)]),
                            (sp.Rational(3, 2), sp.Rational(3, 2), [sp.Rational(3, 4), 3, sp.Rational(1, 2)]),
                            (-1, -sp.Rational(3, 8), [-sp.Rational(1, 2), sp.Rational(3, 8), -sp.Rational(3, 4)])]:
    assert abs(sp.Rational(lam)) / 2 < 1  # |||K||| = |lam|/2 < 1 : la série de Neumann converge
    sol = sp.solve(sp.Eq(c, lam * (m_un + c * m_iota)), c)
    assert sol == [c_att]
    # série de Neumann : f = 1 + sum_{k>=1} lam^k m(1) 3^{-(k-1)} x
    neumann = sp.summation(lam ** (k + 1) * m_un / sp.Integer(3) ** k, (k, 0, sp.oo))
    assert sp.simplify(neumann - c_att) == 0
    for l in leurres:  # f = 1 + l x ne vérifie pas l'équation
        assert sp.simplify(l - lam * (m_un + l * m_iota)) != 0
# remarque : |||K^k||| = 1/(2*3^(k-1)) pour lam = 1 (sup de |m(f)| = 1/2, K^k f = m(f) x / 3^(k-1))

# uniforme-ex5 : fonction de Takagi aux points dyadiques (somme finie), r_k entier de parité k
def phi(y):
    y = Fr(y)
    fl = y.numerator // y.denominator
    return min(y - fl, fl + 1 - y)


def takagi(y):
    y = Fr(y)
    tot, j = Fr(0), 0
    while True:
        terme = phi(2 ** j * y) / 2 ** j
        if (2 ** j * y).denominator == 1:
            return tot
        tot += terme
        j += 1


assert takagi(Fr(1, 4)) == Fr(1, 2) and takagi(Fr(3, 8)) == Fr(5, 8)
xv = Fr(1, 3)
u3 = Fr((8 * xv.numerator) // xv.denominator, 8)
assert u3 == Fr(1, 4)
r3 = (takagi(u3 + Fr(1, 8)) - takagi(u3)) / Fr(1, 8)
assert r3 == 1 and all(r3 != l for l in (0, 2, 3))
random.seed(1)
for _ in range(200):
    xv = Fr(random.randint(-10 ** 6, 10 ** 6), random.randint(1, 10 ** 6))
    for kk in range(0, 14):
        uk = Fr((2 ** kk * xv.numerator) // xv.denominator, 2 ** kk)
        rk = (takagi(uk + Fr(1, 2 ** kk)) - takagi(uk)) * 2 ** kk
        assert rk.denominator == 1 and (rk.numerator - kk) % 2 == 0 and uk <= xv < uk + Fr(1, 2 ** kk)

# ===========================================================================
# Carte top-l3-banach-exponentielle
# ===========================================================================
A = sp.Matrix([[0, 1], [0, 0]])
B = sp.Matrix([[0, 0], [1, 0]])
I2 = sp.eye(2)
eA, eB, eAB = A.exp(), B.exp(), (A + B).exp()
assert eA == I2 + A and eB == I2 + B
assert eA * eB == sp.Matrix([[2, 1], [1, 1]]) and eB * eA == sp.Matrix([[1, 1], [1, 2]])
assert sp.simplify(eAB - sp.Matrix([[sp.cosh(1), sp.sinh(1)], [sp.sinh(1), sp.cosh(1)]])) == sp.zeros(2)
assert sp.simplify((eA * eB).trace() - eAB.trace() - (3 - 2 * sp.cosh(1))) == 0
assert (3 - 2 * sp.cosh(1)).evalf() < 0  # ≈ -0,086 : e^{A+B} ≠ e^A e^B (carte, affirmation, trou 2, ex3)
assert A * B - B * A == sp.Matrix([[1, 0], [0, -1]])
# trou 2, leurre « u et v inversibles » : u = I + A, v = I + B
U, V_ = I2 + A, I2 + B
assert U.det() != 0 and V_.det() != 0 and sp.simplify((U + V_).exp() - U.exp() * V_.exp()) != sp.zeros(2)
# leurre « u^2 = v^2 » : A^2 = B^2 = 0 ; leurre « |||u||| = |||v||| » : normes subordonnées égales (1 pour |.|_2)
assert A ** 2 == sp.zeros(2) and B ** 2 == sp.zeros(2)
assert max((A.T * A).eigenvals()) == 1 and max((B.T * B).eigenvals()) == 1
# trou 1 : E = R, u = -1 : e^{-1} ≠ e, e^{-1} < e ; u = 1 : e > 2 = 1 + |u|
assert sp.exp(-1) < sp.E and sp.E > 2
# trou 3 : u = 0 (dérivée 0, leurre e^{0} = 1) ; u = 1 en t = 0 (u e^{0} = 1 ≠ 0) ; u = -1 (|u| e^{-t} ≠ -e^{-t})
# rotation J : |||J||| = theta, e^J rotation de norme 1 ; e^{2 pi J1} = I
th = sp.symbols("theta", positive=True)
J = sp.Matrix([[0, -th], [th, 0]])
eJ = sp.simplify(J.exp())
assert sp.simplify(eJ - sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])) == sp.zeros(2)
assert sp.simplify(eJ.T * eJ) == I2
assert sp.simplify(eJ.subs(th, 2 * sp.pi)) == I2
# exemple : projecteur, e^{tp} = id + (e^t - 1) p ; multiplication par g(x) = -x : sup e^{-x} = 1 < e
P = sp.Matrix([[1, 1], [0, 0]])
assert sp.simplify((t * P).exp() - (I2 + (sp.exp(t) - 1) * P)) == sp.zeros(2)

# exponentielle-ex1 (v1, v2) : e^p pour trois projecteurs, leurres
for Pm, leurres in [
    (sp.Matrix([[1, 1], [0, 0]]), [[[sp.E, sp.E], [0, 1]], [[sp.E, sp.E - 1], [0, 0]], [[sp.E, 1], [0, 1]]]),
    (sp.Matrix([[0, 0], [2, 1]]), [[[1, 0], [2 * sp.E, sp.E]], [[0, 0], [2 * (sp.E - 1), sp.E]], [[1, 0], [2, sp.E]]]),
    (sp.Matrix([[1, 0], [3, 0]]), [[[sp.E, 0], [3 * sp.E, 1]], [[sp.E, 0], [3 * (sp.E - 1), 0]], [[sp.E, 0], [3, 1]]]),
]:
    assert Pm ** 2 == Pm
    eP = sp.simplify(Pm.exp())
    assert sp.simplify(eP - (I2 + (sp.E - 1) * Pm)) == sp.zeros(2)
    for l in leurres:
        assert sp.simplify(eP - sp.Matrix(l)) != sp.zeros(2)
assert sp.Matrix([[1, 1], [0, 0]]).exp().applyfunc(sp.simplify) == sp.Matrix([[sp.E, sp.E - 1], [0, 1]])

# exponentielle-ex2 (v1, v2) : |||M_g||| = sup|g|, |||e^{M_g}||| = sup e^g sur [0,1] (g affine : extrêmes aux bords)
for g, nM, nE, leurres in [(-2 * x, 2, 1, [(2, sp.E ** 2), (2, sp.exp(-2)), (1, 1)]),
                           (1 - 3 * x, 2, sp.E, [(2, sp.E ** 2), (1, sp.E), (2, sp.exp(-2))]),
                           (x - 2, 2, sp.exp(-1), [(2, sp.E ** 2), (1, sp.exp(-1)), (2, sp.exp(-2))])]:
    bords = [g.subs(x, 0), g.subs(x, 1)]
    assert max(abs(b) for b in bords) == nM
    assert max(sp.exp(b) for b in bords) == nE
    assert all(l != (nM, nE) for l in leurres)
    assert nE < sp.exp(nM)

# exponentielle-ex3 (v1, v2) : traces
for a, b, attendu, leurres in [(1, 1, 3 - 2 * sp.cosh(1), [0, 3 - 2 * sp.cos(1), 2 - 2 * sp.cosh(1)]),
                               (2, 2, 6 - 2 * sp.cosh(2), [0, 6 - 2 * sp.cos(2), 5 - 2 * sp.cosh(2)]),
                               (-1, 1, 1 - 2 * sp.cos(1), [0, 1 - 2 * sp.cosh(1), 2 - 2 * sp.cos(1)])]:
    Aa, Bb = sp.Matrix([[0, a], [0, 0]]), sp.Matrix([[0, 0], [b, 0]])
    val = sp.simplify((Aa.exp() * Bb.exp()).trace() - (Aa + Bb).exp().trace())
    assert sp.simplify((val - attendu).rewrite(sp.exp)) == 0
    for l in leurres:
        assert abs(sp.N(attendu - l)) > 1e-3
assert sp.Matrix([[0, 2], [0, 0]]).exp() * sp.Matrix([[0, 0], [2, 0]]).exp() == sp.Matrix([[5, 2], [2, 1]])
assert sp.Matrix([[0, -1], [0, 0]]).exp() * sp.Matrix([[0, 0], [1, 0]]).exp() == sp.Matrix([[0, -1], [1, 1]])
assert sp.Matrix([[0, 0], [1, 0]]).exp() * sp.Matrix([[0, -1], [0, 0]]).exp() == sp.Matrix([[1, -1], [1, 0]])


def norme1(M):  # norme subordonnée à |.|_1 (max des sommes de colonnes)
    return mp.mnorm(M, 1)


# exponentielle-ex4 : |||e^u - e^v||| <= |||u - v||| e^{max} sur des matrices aléatoires ; contre-exemples dans R
random.seed(2)
for _ in range(100):
    Ur = mp.matrix([[random.uniform(-2, 2) for _ in range(3)] for _ in range(3)])
    Vr = mp.matrix([[random.uniform(-2, 2) for _ in range(3)] for _ in range(3)])
    r = max(norme1(Ur), norme1(Vr))
    assert norme1(mp.expm(Ur) - mp.expm(Vr)) <= norme1(Ur - Vr) * mp.e ** r + mp.mpf(10) ** -20
assert mp.e - 1 > 1                                   # leurre |||u - v||| (u = 1, v = 0)
assert mp.e ** 10 - mp.e ** 9 > mp.e - 1              # leurre e^{|||u-v|||} - 1 (u = 10, v = 9)
assert mp.e - 1 > 1 * mp.e ** 0                       # leurre e^{min} (u = 1, v = 0)

# exponentielle-ex5 : formule de Trotter ; |||e^w - id - w||| <= |||w|||^2 e^{|||w|||} / 2 ; erreur <= C e^{a+b}/n
for _ in range(50):
    W = mp.matrix([[random.uniform(-2, 2) for _ in range(3)] for _ in range(3)])
    cw = norme1(W)
    assert norme1(mp.expm(W) - mp.eye(3) - W) <= cw ** 2 * mp.e ** cw / 2 + mp.mpf(10) ** -20
Am, Bm = mp.matrix([[0, 1], [0, 0]]), mp.matrix([[0, 0], [1, 0]])
a_, b_ = norme1(Am), norme1(Bm)
C = a_ * b_ + (1 + a_) * b_ ** 2 * mp.e ** b_ / 2 + a_ ** 2 * mp.e ** (a_ + b_) / 2 + (a_ + b_) ** 2 * mp.e ** (a_ + b_) / 2
eApB = mp.expm(Am + Bm)
for nn in (1, 2, 5, 10, 50, 200):
    Xn = mp.expm(Am / nn) * mp.expm(Bm / nn)
    Yn = mp.expm((Am + Bm) / nn)
    assert norme1(Xn - Yn) <= C / nn ** 2
    assert norme1(Xn ** nn - eApB) <= C * mp.e ** (a_ + b_) / nn
    # e^{A/n} e^{B/n} = [[1 + 1/n^2, 1/n], [1/n, 1]] et coefficient (2,2) de Y_n = ch(1/n) > 1
    assert norme1(Xn - mp.matrix([[1 + mp.mpf(1) / nn ** 2, mp.mpf(1) / nn], [mp.mpf(1) / nn, 1]])) < 1e-25
    assert Yn[1, 1] - 1 > 0 and abs(Yn[1, 1] - mp.cosh(mp.mpf(1) / nn)) < 1e-25
# la limite est e^{A+B} et non e^A e^B
assert norme1(mp.expm(Am) * mp.expm(Bm) - eApB) > 0.05
# matrices qui commutent : (e^{u/n} e^{v/n})^n = e^{u+v}
D1, D2 = mp.matrix([[1, 2], [0, 1]]), mp.matrix([[3, -1], [0, 3]])
assert norme1(D1 * D2 - D2 * D1) < 1e-25
for nn in (1, 3, 7):
    assert norme1((mp.expm(D1 / nn) * mp.expm(D2 / nn)) ** nn - mp.expm(D1 + D2)) < 1e-15

# ===========================================================================
# Carte top-l3-banach-bessel
# ===========================================================================
e0, e1 = sp.Integer(1), sp.sqrt(3) * (2 * t - 1)
ps = lambda f, g: sp.integrate(f * g, (t, 0, 1))
assert ps(e0, e0) == 1 and ps(e1, e1) == 1 and ps(e0, e1) == 0
# exemple : t^4 -> c0 = 1/5, c1 = 2 sqrt(3)/15, d^2 = 4/225
assert ps(e0, t ** 4) == sp.Rational(1, 5) and sp.simplify(ps(e1, t ** 4) - 2 * sp.sqrt(3) / 15) == 0
assert sp.simplify(ps(t ** 4, t ** 4) - ps(e0, t ** 4) ** 2 - ps(e1, t ** 4) ** 2) == sp.Rational(4, 225)
# exemple R^3 : e0 = (1,0,0), e1 = (0,1,1)/sqrt2, x = (0,1,0)
E0, E1 = sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 1]) / sp.sqrt(2)
xv3 = sp.Matrix([0, 1, 0])
assert (E0.dot(xv3)) ** 2 + (E1.dot(xv3)) ** 2 == sp.Rational(1, 2)
assert sp.simplify(E0.dot(xv3) * E0 + E1.dot(xv3) * E1 - sp.Matrix([0, 1, 1]) / 2) == sp.zeros(3, 1)
# contre-exemples (produit scalaire des coefficients) : 4 > 1 ; 1 + 2 = 3 > 2
assert (2 * 1) ** 2 > 1 and 1 ** 2 + (sp.Rational(2) / sp.sqrt(2)) ** 2 == 3 > 2
# affirmation : x = (1/2^n) dans l^2, sum 1/4^n = 4/3
assert sp.summation(sp.Rational(1, 4) ** k, (k, 0, sp.oo)) == sp.Rational(4, 3)

# bessel-ex1 (v1, v2) : distance au plan Vect(e0, e1)
for xv, d_att, leurres in [((1, 2, 3), sp.sqrt(2) / 2, [sp.Rational(1, 2), sp.sqrt(2), sp.sqrt(14) / 2]),
                           ((1, -1, 2), 3 * sp.sqrt(2) / 2, [sp.Rational(9, 2), sp.sqrt(6), sp.sqrt(2) / 2]),
                           ((2, 3, -1), 2 * sp.sqrt(2), [8, sp.sqrt(14), sp.sqrt(6)])]:
    xm = sp.Matrix(xv)
    proj = E0.dot(xm) * E0 + E1.dot(xm) * E1
    d = sp.sqrt((xm - proj).dot(xm - proj))
    assert sp.simplify(d - d_att) == 0
    assert E0.dot(xm) ** 2 + E1.dot(xm) ** 2 <= xm.dot(xm)
    differents(d_att, leurres)

# bessel-ex2 (v1, v2) : min de int (f - alpha - beta t)^2 par les équations normales
al, be = sp.symbols("alpha beta")
for f, mini, proj_att, leurres in [
    (t ** 2, sp.Rational(1, 180), t - sp.Rational(1, 6), [sp.Rational(1, 5), sp.Rational(7, 36), sp.Rational(1, 12)]),
    (t ** 3, sp.Rational(9, 700), sp.Rational(9, 10) * t - sp.Rational(1, 5), [sp.Rational(1, 7), sp.Rational(13, 100), sp.Rational(1, 28)]),
    (sp.sqrt(t), sp.Rational(1, 450), sp.Rational(4, 5) * t + sp.Rational(4, 15), [sp.Rational(1, 2), sp.Rational(112, 225), sp.Rational(1, 18)]),
]:
    J_ = sp.integrate((f - al - be * t) ** 2, (t, 0, 1))
    sol = sp.solve([sp.diff(J_, al), sp.diff(J_, be)], [al, be], dict=True)[0]
    assert sp.simplify(J_.subs(sol) - mini) == 0
    assert sp.simplify(sol[al] + sol[be] * t - proj_att) == 0
    c0, c1 = ps(e0, f), ps(e1, f)
    assert sp.simplify(ps(f, f) - c0 ** 2 - c1 ** 2 - mini) == 0
    differents(mini, leurres)
assert sp.Rational(1, 16) + sp.Rational(27, 400) == sp.Rational(13, 100)
assert sp.Rational(4, 9) + sp.Rational(4, 75) == sp.Rational(112, 225)

# bessel-ex3 (v1, v2) : famille sqrt2 sin(n pi t), coefficients et sommes de Parseval
en = sp.sqrt(2) * sp.sin(n * sp.pi * t)
m_ = sp.symbols("m", integer=True, positive=True)
for nn in range(1, 7):
    for mm in range(1, 7):
        val = sp.integrate(2 * sp.sin(nn * sp.pi * t) * sp.sin(mm * sp.pi * t), (t, 0, 1))
        assert val == (1 if nn == mm else 0)
for f, coef_att, somme_att, leurres in [
    (sp.Integer(1), sp.sqrt(2) * (1 - (-1) ** n) / (n * sp.pi), 1, [8 / sp.pi ** 2, sp.Rational(1, 2), sp.pi ** 2 / 8]),
    (t, sp.sqrt(2) * (-1) ** (n + 1) / (n * sp.pi), sp.Rational(1, 3), [2 / sp.pi ** 2, sp.Rational(1, 6), sp.Rational(1, 2)]),
    (1 + t, sp.sqrt(2) * (1 - 2 * (-1) ** n) / (n * sp.pi), sp.Rational(7, 3), [sp.Rational(9, 4), 18 / sp.pi ** 2, 2]),
]:
    for nn in range(1, 12):
        assert sp.simplify(sp.integrate(f * en.subs(n, nn), (t, 0, 1)) - coef_att.subs(n, nn)) == 0
    impairs = sp.summation((coef_att.subs(n, 2 * k + 1)) ** 2, (k, 0, sp.oo))
    pairs = sp.summation((coef_att.subs(n, 2 * k + 2)) ** 2, (k, 0, sp.oo))
    assert sp.simplify(impairs + pairs - somme_att) == 0
    assert sp.simplify(sp.integrate(f ** 2, (t, 0, 1)) - somme_att) == 0
    differents(somme_att, leurres)
assert sp.summation(1 / (2 * k + 1) ** 2, (k, 0, sp.oo)) == sp.pi ** 2 / 8
assert sp.summation(1 / (4 * k ** 2), (k, 1, sp.oo)) == sp.pi ** 2 / 24

# bessel-ex4 : au plus ||x||^2/eps^2 = 4 indices ; x = (1 + X + X^2 + X^3)/2 de norme 1 atteint la borne
assert Fr(1) / Fr(1, 4) == 4 and sum(Fr(1, 2) ** 2 for _ in range(4)) == 1

# bessel-ex5 : B = 1 + 1/sqrt2, meilleure constante = plus grande valeur propre de M
M = sp.Matrix([[sp.Rational(3, 2), sp.Rational(1, 2)], [sp.Rational(1, 2), sp.Rational(1, 2)]])
u1, u2 = sp.Matrix([1, 0]), sp.Matrix([1, 1]) / sp.sqrt(2)
assert sp.simplify(u1 * u1.T + u2 * u2.T - M) == sp.zeros(2)
vp = list(M.eigenvals())
assert sp.simplify(max(vp, key=lambda v: sp.N(v)) - (1 + sp.sqrt(2) / 2)) == 0
Bc = 1 + abs(u1.dot(u2))
assert sp.simplify(Bc - (1 + sp.sqrt(2) / 2)) == 0
differents(1 + sp.sqrt(2) / 2, [2, sp.Rational(3, 2), 1 - sp.sqrt(2) / 2])
# inégalité de la question 3 sur des familles aléatoires
for _ in range(100):
    N_ = random.randint(1, 5)
    us = [mp.matrix([random.uniform(-1, 1) for _ in range(4)]) for _ in range(N_)]
    xr = mp.matrix([random.uniform(-1, 1) for _ in range(4)])
    dot = lambda a, b: sum(a[i] * b[i] for i in range(4))
    Bv = max(sum(abs(dot(us[i], us[j])) for j in range(N_)) for i in range(N_))
    assert sum(dot(u, xr) ** 2 for u in us) <= Bv * dot(xr, xr) + mp.mpf(10) ** -20

# ===========================================================================
# Carte top-l3-banach-hilbert
# ===========================================================================
# exemple : <x, y> = sum 1/((n+1) 2^n) = 2 ln 2 ; ||x||^2 = pi^2/6 ; ||y||^2 = 4/3
assert sp.simplify(sp.summation(1 / ((k + 1) * 2 ** k), (k, 0, sp.oo)) - 2 * sp.log(2)) == 0
# trou 2 : a_n = (-1)^n / sqrt(n+1) : sum a_n converge (Leibniz), sum a_n^2 diverge
assert sp.summation(1 / (k + 1), (k, 0, sp.oo)) == sp.oo
assert abs(mp.nsum(lambda j: (-1) ** j / mp.sqrt(j + 1), [0, mp.inf])) < 1
# trou 3 : a = (1,1) : sqrt2 ≠ 2 ; a0 = 2 : 2 ≠ 4 ; inégalité stricte fausse (égalité)
assert sp.sqrt(2) != 2 and 2 != 4
# carte (C([0,1]) non complet) et hilbert-ex3 (v1, v2) : int_0^1 (f_n - h)^2
nn_ = sp.symbols("N", positive=True)
I_main = sp.integrate((1 - nn_ * s) ** 2, (s, 0, 1 / nn_))
assert sp.simplify(I_main - 1 / (3 * nn_)) == 0
I_v1 = sp.integrate((nn_ * s / 2) ** 2, (s, 0, 1 / nn_)) + sp.integrate(((1 - nn_ * s) / 2) ** 2, (s, 0, 1 / nn_))
assert sp.simplify(I_v1 - 1 / (6 * nn_)) == 0
I_v2 = sp.integrate((1 - nn_ ** 2 * s) ** 2, (s, 0, 1 / nn_ ** 2))
assert sp.simplify(I_v2 - 1 / (3 * nn_ ** 2)) == 0
differents(1 / (3 * n), [1 / n, 1 / (2 * n), 1 / (3 * n ** 2)], valeurs=range(2, 30))
differents(1 / (6 * n), [1 / (12 * n), 1 / (3 * n), 1 / (6 * n ** 2)], valeurs=range(2, 30))
differents(1 / (3 * n ** 2), [1 / n ** 2, 1 / (3 * n), 1 / (3 * n ** 4)], valeurs=range(2, 30))
# v1 : f_n(t) = 1/2 + (n/2)(t - 1/2) vaut 0 en 1/2 - 1/n et 1 en 1/2 + 1/n
fn = sp.Rational(1, 2) + nn_ / 2 * (t - sp.Rational(1, 2))
assert sp.simplify(fn.subs(t, sp.Rational(1, 2) - 1 / nn_)) == 0 and sp.simplify(fn.subs(t, sp.Rational(1, 2) + 1 / nn_)) == 1

# hilbert-ex1 (v1, v2)
for rx, ry, attendu, nx, ny, leurres in [
    (sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(6, 5), sp.Rational(4, 3), sp.Rational(9, 8),
     [sp.Rational(3, 2), sp.Rational(5, 6), sp.sqrt(6) / 2]),
    (sp.Rational(1, 2), -sp.Rational(1, 2), sp.Rational(4, 5), sp.Rational(4, 3), sp.Rational(4, 3),
     [sp.Rational(4, 3), 0, -sp.Rational(4, 5)]),
    (sp.Rational(1, 3), -sp.Rational(1, 2), sp.Rational(6, 7), sp.Rational(9, 8), sp.Rational(4, 3),
     [sp.Rational(6, 5), sp.Rational(7, 6), -sp.Rational(6, 7)]),
]:
    assert sp.summation((rx * ry) ** k, (k, 0, sp.oo)) == attendu
    assert sp.summation((rx ** 2) ** k, (k, 0, sp.oo)) == nx and sp.summation((ry ** 2) ** k, (k, 0, sp.oo)) == ny
    assert attendu ** 2 <= nx * ny
    differents(attendu, leurres)
assert sp.sqrt(sp.Rational(4, 3) * sp.Rational(9, 8)) == sp.sqrt(6) / 2

# hilbert-ex2 (v1, v2) : d(e0, v^perp) = |v_0| / ||v||
for nv2, d_att, leurres in [(sp.Rational(4, 3), sp.sqrt(3) / 2, [sp.Rational(3, 4), 2 / sp.sqrt(3), 1]),
                            (sp.Rational(9, 8), 2 * sp.sqrt(2) / 3, [sp.Rational(8, 9), 3 / (2 * sp.sqrt(2)), 1]),
                            (sp.pi ** 2 / 6, sp.sqrt(6) / sp.pi, [6 / sp.pi ** 2, sp.pi / sp.sqrt(6), 1])]:
    assert sp.simplify(1 / sp.sqrt(nv2) - d_att) == 0
    differents(d_att, leurres)
assert sp.summation(sp.Rational(1, 9) ** k, (k, 0, sp.oo)) == sp.Rational(9, 8)

# hilbert-ex4 : sum 4/(n+1)^2 et sum 1/4^n convergent ; e_k/(k+1) -> 0
assert sp.summation(4 / (k + 1) ** 2, (k, 0, sp.oo)) == 2 * sp.pi ** 2 / 3

# hilbert-ex5 : f_n(t) = min(1, n t) sur [0,1] : int_0^1 (1 - f_n)^2 = 1/(3n) ; d(1, F)^2 = 1
assert sp.simplify(sp.integrate((1 - nn_ * t) ** 2, (t, 0, 1 / nn_)) - 1 / (3 * nn_)) == 0
assert sp.integrate(1, (t, -1, 0)) == 1
assert all(l != 1 for l in (0, 2, sp.Rational(1, 2)))

# ===========================================================================
# Carte top-l3-banach-polynomes-orthogonaux
# ===========================================================================
def gram_schmidt(moment, N):
    """Polynômes unitaires orthogonaux P_0..P_N pour le produit scalaire défini par les moments m_k."""
    def prod(Pp, Qq):
        c = sp.Poly(sp.expand(Pp * Qq), X).all_coeffs()[::-1]
        return sp.simplify(sum(ci * moment(i) for i, ci in enumerate(c)))
    Ps = []
    for d in range(N + 1):
        Pd = X ** d - sum(prod(Pk, X ** d) / prod(Pk, Pk) * Pk for Pk in Ps)
        Ps.append(sp.expand(Pd))
    return Ps, prod


poids = {
    "legendre": (lambda i: sp.Rational(2, i + 1) if i % 2 == 0 else 0, (-1, 1)),
    "uniforme01": (lambda i: sp.Rational(1, i + 1), (0, 1)),
    "uniforme02": (lambda i: sp.Rational(2 ** (i + 1), i + 1), (0, 2)),
    "hermite": (lambda i: sp.gamma(sp.Rational(i + 1, 2)) if i % 2 == 0 else 0, (-sp.oo, sp.oo)),
    "laguerre": (lambda i: sp.factorial(i), (0, sp.oo)),
    "tchebychev": (lambda i: sp.pi * sp.binomial(i, i // 2) / 2 ** i if i % 2 == 0 else 0, (-1, 1)),
}
# moments de Tchebychev : int_{-1}^1 t^{2j}/sqrt(1-t^2) dt = int_0^pi cos^{2j} = pi C(2j, j)/4^j
for j in range(5):
    assert sp.simplify(sp.integrate(sp.cos(t) ** (2 * j), (t, 0, sp.pi)) - sp.pi * sp.binomial(2 * j, j) / 4 ** j) == 0
# moments d'Hermite : int t^{2j} e^{-t^2} = Gamma(j + 1/2) ; m0 = sqrt(pi), m2 = sqrt(pi)/2, m4 = 3 sqrt(pi)/4
tt = sp.symbols("tt", real=True)
for j in range(4):
    assert sp.simplify(sp.integrate(tt ** (2 * j) * sp.exp(-tt ** 2), (tt, -sp.oo, sp.oo)) - sp.gamma(sp.Rational(2 * j + 1, 2))) == 0
assert poids["hermite"][0](2) == sp.sqrt(sp.pi) / 2 and poids["hermite"][0](4) == 3 * sp.sqrt(sp.pi) / 4

GS = {}
for nom, (mom, (aa, bb)) in poids.items():
    Ps, prod = gram_schmidt(mom, 7)
    GS[nom] = (Ps, prod)
    for i in range(8):
        for j in range(i):
            assert prod(Ps[i], Ps[j]) == 0
    # récurrence à trois termes (carte)
    for d in range(0, 7):
        a_d = prod(X * Ps[d], Ps[d]) / prod(Ps[d], Ps[d])
        rec = (X - a_d) * Ps[d] - (prod(Ps[d], Ps[d]) / prod(Ps[d - 1], Ps[d - 1]) * Ps[d - 1] if d >= 1 else 0)
        assert sp.expand(rec - Ps[d + 1]) == 0
        if d >= 1:
            assert prod(Ps[d], Ps[d]) / prod(Ps[d - 1], Ps[d - 1]) > 0
    # racines réelles, simples, dans ]a, b[ (carte)
    for d in range(1, 8):
        racines = sp.Poly(Ps[d], X).nroots(n=30)
        assert all(abs(sp.im(r)) < 1e-20 for r in racines)
        reelles = sorted(sp.re(r) for r in racines)
        assert all(reelles[i + 1] - reelles[i] > 1e-10 for i in range(d - 1))
        assert all((aa == -sp.oo or r > aa) and (bb == sp.oo or r < bb) for r in reelles)

PL, prodL = GS["legendre"]
# carte : Legendre P_2, P_3, normes, b_1 = 1/3, b_2 = 4/15, b_n = n^2/(4n^2 - 1)
assert PL[2] == X ** 2 - sp.Rational(1, 3) and PL[3] == X ** 3 - sp.Rational(3, 5) * X
assert [prodL(PL[i], PL[i]) for i in range(3)] == [2, sp.Rational(2, 3), sp.Rational(8, 45)]
for d in range(1, 7):
    assert prodL(PL[d], PL[d]) / prodL(PL[d - 1], PL[d - 1]) == sp.Rational(d ** 2, 4 * d ** 2 - 1)
# Hermite : b_n = n/2, P_3 = X^3 - 3X/2 ; Laguerre : P_3 = X^3 - 9X^2 + 18X - 6, b_n = n^2, a_n = 2n + 1
PH, prodH = GS["hermite"]
for d in range(1, 7):
    assert sp.simplify(prodH(PH[d], PH[d]) / prodH(PH[d - 1], PH[d - 1]) - sp.Rational(d, 2)) == 0
PLa, prodLa = GS["laguerre"]
for d in range(0, 6):
    assert prodLa(X * PLa[d], PLa[d]) / prodLa(PLa[d], PLa[d]) == 2 * d + 1
    if d >= 1:
        assert prodLa(PLa[d], PLa[d]) / prodLa(PLa[d - 1], PLa[d - 1]) == d ** 2
# Tchebychev : P_n = 2^{1-n} T_n, racines cos((2k-1) pi / (2n))
PT, prodT = GS["tchebychev"]
T_ = [sp.Integer(1), X]
for d in range(1, 7):
    T_.append(sp.expand(2 * X * T_[d] - T_[d - 1]))
for d in range(1, 8):
    assert sp.expand(PT[d] - T_[d] / 2 ** (d - 1)) == 0
    for kk in range(1, d + 1):
        assert abs(sp.N(T_[d].subs(X, sp.cos((2 * kk - 1) * sp.pi / (2 * d))), 30)) < 1e-25
for d in range(6):
    assert sp.simplify(sp.expand_trig(T_[d].subs(X, sp.cos(t)) - sp.cos(d * t))) == 0

# trous : leurres de a_n (poids 1 sur ]0,1[, n = 1 : a_1 = 1/2) et de b_n (Legendre, n = 1 : b_1 = 1/3)
P01, prod01 = GS["uniforme01"]
a1 = prod01(X * P01[1], P01[1]) / prod01(P01[1], P01[1])
assert a1 == sp.Rational(1, 2)
assert prod01(X * P01[1], P01[1]) == sp.Rational(1, 24) != a1                       # leurre <X P_n, P_n>
assert prod01(X * P01[1], P01[0]) / prod01(P01[1], P01[1]) == 1 != a1               # leurre <XP_n,P_{n-1}>/||P_n||^2
assert prod01(P01[1], P01[1]) / prod01(X * P01[1], P01[1]) == 2 != a1               # leurre ||P_n||^2/<XP_n,P_n>
assert PL[1] == X and prodL(X * PL[1], PL[0]) / prodL(PL[1], PL[1]) == 1 != 0       # Legendre : a_1 = 0
b1 = prodL(PL[1], PL[1]) / prodL(PL[0], PL[0])
assert b1 == sp.Rational(1, 3)
for l in (1 / b1, sp.sqrt(b1), -b1):
    assert l != b1
# trou 1 : pour n = 2 (Legendre), deux racines simples dans ]-1,1[ : les trois leurres sont faux
assert sp.solve(PL[2], X) == [-sp.sqrt(3) / 3, sp.sqrt(3) / 3]
# contre-exemples de la carte
assert sp.integrate(t * (t ** 2 + 1), (t, -1, 1)) == 0 and sp.integrate(t ** 2 + 1, (t, -1, 1)) == sp.Rational(8, 3)
assert sp.integrate(t, (t, -2, -1)) + sp.integrate(t, (t, 1, 2)) == 0
assert sp.integrate(t ** 2 - 3, (t, -1, 1)) == -sp.Rational(16, 3)
# affirmation : P_2 et P_3 de Legendre n'ont pas de racine commune
assert sp.gcd(sp.Poly(PL[2], X), sp.Poly(PL[3], X)).degree() == 0

# polynomes-orthogonaux-ex1 (v1, v2) : P_2 et leurres (non orthogonaux à 1 ou à X)
for nom, P2, leurres in [
    ("legendre", X ** 2 - sp.Rational(1, 3), [X ** 2 - sp.Rational(1, 2), X ** 2 - sp.Rational(2, 3), X ** 2 + sp.Rational(1, 3)]),
    ("uniforme01", X ** 2 - X + sp.Rational(1, 6), [X ** 2 - X + sp.Rational(1, 4), X ** 2 - X / 2 + sp.Rational(1, 6), X ** 2 - X - sp.Rational(1, 6)]),
    ("uniforme02", X ** 2 - 2 * X + sp.Rational(2, 3), [X ** 2 - 2 * X + 1, X ** 2 - 2 * X + sp.Rational(1, 3), X ** 2 - X + sp.Rational(2, 3)]),
]:
    Ps, prod = GS[nom]
    assert sp.expand(Ps[2] - P2) == 0
    for l in leurres:
        assert prod(l, 1) != 0 or prod(l, X) != 0
assert GS["uniforme01"][0][1] == X - sp.Rational(1, 2) and GS["uniforme02"][0][1] == X - 1

# polynomes-orthogonaux-ex2 (v1, v2) : racines de T_3, T_4, T_5 et leurres
assert T_[3] == 4 * X ** 3 - 3 * X and T_[4] == 8 * X ** 4 - 8 * X ** 2 + 1 and T_[5] == 16 * X ** 5 - 20 * X ** 3 + 5 * X
for d, bonnes, leurres in [
    (3, [0, sp.sqrt(3) / 2, -sp.sqrt(3) / 2],
     [[0, sp.Rational(1, 2), -sp.Rational(1, 2)], [sp.sqrt(3) / 2, -sp.sqrt(3) / 2, 1], [0, sp.sqrt(2) / 2, -sp.sqrt(2) / 2]]),
    (4, [sp.cos(sp.pi / 8), -sp.cos(sp.pi / 8), sp.cos(3 * sp.pi / 8), -sp.cos(3 * sp.pi / 8)],
     [[sp.cos(sp.pi / 4), -sp.cos(sp.pi / 4), sp.cos(3 * sp.pi / 4), -sp.cos(3 * sp.pi / 4)],
      [sp.cos(sp.pi / 8), -sp.cos(sp.pi / 8), sp.cos(sp.pi / 4), -sp.cos(sp.pi / 4)],
      [0, sp.sqrt(2) / 2, -sp.sqrt(2) / 2, 1]]),
    (5, [0, sp.cos(sp.pi / 10), -sp.cos(sp.pi / 10), sp.cos(3 * sp.pi / 10), -sp.cos(3 * sp.pi / 10)],
     [[0, sp.cos(sp.pi / 5), -sp.cos(sp.pi / 5), sp.cos(2 * sp.pi / 5), -sp.cos(2 * sp.pi / 5)],
      [sp.cos(sp.pi / 10), -sp.cos(sp.pi / 10), sp.cos(3 * sp.pi / 10), -sp.cos(3 * sp.pi / 10), 1],
      [0, sp.cos(sp.pi / 10), -sp.cos(sp.pi / 10), sp.cos(sp.pi / 5), -sp.cos(sp.pi / 5)]]),
]:
    assert all(abs(sp.N(T_[d].subs(X, r), 30)) < 1e-25 for r in bonnes) and len({sp.N(r, 20) for r in bonnes}) == d
    for l in leurres:
        assert any(abs(sp.N(T_[d].subs(X, r), 30)) > 1e-5 for r in l) or len({sp.N(r, 20) for r in l}) < d
assert sp.simplify(sp.cos(sp.pi / 8) - sp.sqrt(2 + sp.sqrt(2)) / 2) == 0

# polynomes-orthogonaux-ex3 (v1, v2) : P_3 et coefficients intermédiaires
assert GS["hermite"][0][3] == X ** 3 - sp.Rational(3, 2) * X and GS["hermite"][0][2] == X ** 2 - sp.Rational(1, 2)
assert GS["laguerre"][0][2] == X ** 2 - 4 * X + 2 and GS["laguerre"][0][3] == X ** 3 - 9 * X ** 2 + 18 * X - 6
assert GS["uniforme01"][0][3] == X ** 3 - sp.Rational(3, 2) * X ** 2 + sp.Rational(3, 5) * X - sp.Rational(1, 20)
assert prod01(P01[2], P01[2]) == sp.Rational(1, 180) and prod01(X * P01[2], P01[2]) == sp.Rational(1, 360)
assert prodLa(PLa[2], PLa[2]) == 4 and prodLa(X * PLa[2], PLa[2]) == 20
for nom, leurres in [("hermite", [X ** 3 - 3 * X, X ** 3 - X / 2, X ** 3 - X]),
                     ("laguerre", [X ** 3 - 9 * X ** 2 + 18 * X + 6, X ** 3 - 6 * X ** 2 + 18 * X - 6, X ** 3 - 9 * X ** 2 + 12 * X - 6]),
                     ("uniforme01", [X ** 3 - sp.Rational(3, 2) * X ** 2 + sp.Rational(3, 5) * X + sp.Rational(1, 20),
                                     X ** 3 - sp.Rational(3, 2) * X ** 2 + X / 2, X ** 3 - sp.Rational(3, 5) * X])]:
    Ps, prod = GS[nom]
    for l in leurres:
        assert any(prod(l, X ** j) != 0 for j in range(3))
# Hermite : H_3 = 8X^3 - 12X proportionnel à P_3 ; X^3 - 3X est le polynôme pour le poids e^{-t^2/2}
assert sp.expand(sp.hermite(3, X) / 8 - GS["hermite"][0][3]) == 0
mom_demi = lambda i: sp.sqrt(2 * sp.pi) * sp.factorial2(i - 1) if i % 2 == 0 else 0
assert gram_schmidt(mom_demi, 3)[0][3] == X ** 3 - 3 * X

# polynomes-orthogonaux-ex4 (v1, v2) : quadrature de Gauss à deux points
for nom, noeuds, lam, exact, approx, leurres in [
    ("legendre", [-1 / sp.sqrt(3), 1 / sp.sqrt(3)], 1, sp.Rational(2, 5), sp.Rational(2, 9), [sp.Rational(2, 5), sp.Rational(4, 9), sp.Rational(1, 9)]),
    ("tchebychev", [-1 / sp.sqrt(2), 1 / sp.sqrt(2)], sp.pi / 2, 3 * sp.pi / 8, sp.pi / 4, [3 * sp.pi / 8, sp.pi / 2, sp.pi / 8]),
    ("hermite", [-1 / sp.sqrt(2), 1 / sp.sqrt(2)], sp.sqrt(sp.pi) / 2, 3 * sp.sqrt(sp.pi) / 4, sp.sqrt(sp.pi) / 4,
     [3 * sp.sqrt(sp.pi) / 4, sp.sqrt(sp.pi) / 2, sp.sqrt(sp.pi) / 8]),
]:
    Ps, prod = GS[nom]
    mom = poids[nom][0]
    assert all(sp.simplify(Ps[2].subs(X, r)) == 0 for r in noeuds)
    assert sp.simplify(2 * lam - mom(0)) == 0
    for d in range(4):  # exacte jusqu'au degré 3
        assert sp.simplify(lam * sum(r ** d for r in noeuds) - mom(d)) == 0
    assert sp.simplify(mom(4) - exact) == 0
    assert sp.simplify(lam * sum(r ** 4 for r in noeuds) - approx) == 0
    differents(approx, leurres)
    # poids : lambda_i = int L_i w (Lagrange)
    L1 = (X - noeuds[1]) / (noeuds[0] - noeuds[1])
    c1 = sp.Poly(sp.expand(L1), X).all_coeffs()[::-1]
    assert sp.simplify(sum(ci * mom(i) for i, ci in enumerate(c1)) - lam) == 0

# polynomes-orthogonaux-ex5 : Rodrigues, ||L_n||^2, I_n
for d in range(0, 8):
    Ln = sp.diff((X ** 2 - 1) ** d, X, d)
    assert sp.Poly(Ln, X).degree() == d and sp.Poly(Ln, X).LC() == sp.factorial(2 * d) / sp.factorial(d)
    assert sp.expand(sp.factorial(d) / sp.factorial(2 * d) * Ln - PL[d]) == 0
    nL = sp.integrate(Ln ** 2, (X, -1, 1))
    assert nL == 2 ** (2 * d + 1) * sp.factorial(d) ** 2 / (2 * d + 1)
    In = sp.integrate((1 - X ** 2) ** d, (X, -1, 1))
    assert In == 2 ** (2 * d + 1) * sp.factorial(d) ** 2 / sp.factorial(2 * d + 1)
    if d >= 1:
        assert In == sp.Rational(2 * d, 2 * d + 1) * sp.integrate((1 - X ** 2) ** (d - 1), (X, -1, 1))
        for l in (sp.Rational(2, 2 * d + 1), 2 ** (2 * d + 1) * sp.factorial(d) ** 2 / sp.factorial(2 * d + 1),
                  2 ** (d + 1) * sp.factorial(d) ** 2 / (2 * d + 1)):
            assert l != nL
    assert sp.integrate((Ln / (2 ** d * sp.factorial(d))) ** 2, (X, -1, 1)) == sp.Rational(2, 2 * d + 1)
    assert (Ln / (2 ** d * sp.factorial(d))).subs(X, 1) == 1
    for kk in range(d):  # dérivées de U_n d'ordre < n nulles en ±1
        Uk = sp.diff((X ** 2 - 1) ** d, X, kk)
        assert Uk.subs(X, 1) == 0 and Uk.subs(X, -1) == 0

# ===========================================================================
# QCM
# ===========================================================================
# q4 : u = -id : e^{tu} = e^{-t} ; les leurres diffèrent en t = 1
assert sp.simplify(sp.cosh(1) - sp.sinh(1) - sp.exp(-1)) == 0
for l in (sp.cos(1) - sp.sin(1), sp.E, 2 - sp.E):
    assert abs(sp.N(l - sp.exp(-1))) > 0.1
# q5 : e^{-1} < e ; q6 : sum 1/4^n = 4/3, 1 < 4/3 < 16/9 < 2
assert sp.Rational(1) < sp.Rational(4, 3) < sp.Rational(16, 9) < 2
# q8 : b_2 = 4/15 et seule cette valeur rend X P_2 - b X orthogonal à X
for b in (sp.Rational(4, 15), sp.Rational(1, 3), sp.Rational(15, 4), sp.Rational(2, 5)):
    orth = prodL(X * PL[2] - b * X, X) == 0
    assert orth == (b == sp.Rational(4, 15))
assert sp.integrate(t ** 4 - sp.Rational(3, 5) * t ** 2, (t, -1, 1)) == 0

print("top_l3_banach : toutes les vérifications passent.")
