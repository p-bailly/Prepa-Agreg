"""Vérifications calculatoires du chapitre « top-l2-evn » (espaces vectoriels normés, L2) :
cartes, exercices et QCM de content/*/topologie-calcul-diff/l2-top-l2-evn.json et content/qcm/top-l2-evn.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variantes v1, v2) ou de la question vérifiée.
Lancer : python3 verification/top_l2_evn.py (doit finir sans erreur, en quelques secondes).
"""

import itertools
import math
import random
from fractions import Fraction as F

from sympy import Matrix, Rational, Symbol, diff, expand, factor, limit, oo, solve, sqrt, symbols, E, exp, log, sin, cos, pi

random.seed(1)
x, y, t = symbols("x y t", real=True)
n = Symbol("n", positive=True, integer=True)


def n1(v):
    return sum(abs(c) for c in v)


def n2(v):
    return math.sqrt(sum(c * c for c in v))


def ninf(v):
    return max(abs(c) for c in v)


# ---------------------------------------------------------------------------
# top-l2-evn-normes (carte)
# ---------------------------------------------------------------------------
# Énoncé : ||x||_inf <= ||x||_2 <= ||x||_1 <= sqrt(n)||x||_2 <= n||x||_inf (tests aléatoires).
for _ in range(2000):
    d = random.randint(1, 6)
    v = [random.uniform(-5, 5) for _ in range(d)]
    assert ninf(v) <= n2(v) + 1e-12 <= n1(v) + 2e-12 <= math.sqrt(d) * n2(v) + 3e-12 <= d * ninf(v) + 4e-12
# Exemple : x = (1,-2,2) : 5, 3, 2 et 2 <= 3 <= 5 <= 3 sqrt 3 <= 6.
assert (n1([1, -2, 2]), n2([1, -2, 2]), ninf([1, -2, 2])) == (5, 3.0, 2) and 5 <= 3 * math.sqrt(3) <= 6
# Trous 2 et 3 (n >= 2) : chaque leurre est faux pour un vecteur bien choisi.
for d in range(2, 8):
    un = [1] * d
    assert n1(un) > n2(un) and n1(un) > ninf(un) and n1(un) > math.sqrt(d) * ninf(un)  # leurres trou 2
    e1 = [1] + [0] * (d - 1)
    assert ninf(e1) > n1(e1) / d and ninf(e1) > n2(e1) / math.sqrt(d) and ninf(e1) > min(abs(c) for c in e1)  # trou 3
# Contre-exemple : N_inf(t^n) = 1, N_1(t^n) = 1/(n+1).
from sympy import integrate
assert all(integrate(t**k, (t, 0, 1)) == Rational(1, k + 1) for k in range(1, 8))
# Affirmation min(N,N') : N = |x|+2|y|, N' = 2|x|+|y| : min = 1 en (1,0) et (0,1), 3 en (1,1).
mn = lambda a, b: min(abs(a) + 2 * abs(b), 2 * abs(a) + abs(b))
assert mn(1, 0) == 1 and mn(0, 1) == 1 and mn(1, 1) == 3 > 2

# top-l2-evn-normes-ex1 (+ v1, v2) : triplet (||x||_1, ||x||_2, ||x||_inf) et leurres.
for v, bon, leurres in [((3, -4, 12), (19, 13, 12), [(11, 13, 12), (19, 169, 12), (19, 13, 3)]),
                        ((2, -6, 9), (17, 11, 9), [(5, 11, 9), (17, 121, 9), (17, 11, 2)]),
                        ((-1, 4, 8), (13, 9, 8), [(11, 9, 8), (13, 81, 8), (13, 9, 1)])]:
    assert (n1(v), round(n2(v)), ninf(v)) == bon and n2(v) == bon[1]
    assert leurres[0][0] == sum(v) and leurres[1][1] == bon[1] ** 2 and leurres[2][2] == min(abs(c) for c in v)
    assert all(l != bon for l in leurres)
    assert bon[0] ** 2 <= 3 * bon[1] ** 2  # ||x||_1 <= sqrt 3 ||x||_2

# top-l2-evn-normes-ex2 (+ v1, v2) : contre-exemples des non-normes, injectivité des normes.
Na = lambda a, b: max(abs(a), abs(a + b))
assert Matrix([[1, 0], [1, 1]]).det() != 0  # N_a = ||(x, x+y)||_inf
assert abs(1 - 1) == 0  # N_b(1,1) = |1-1| = 0
assert math.sqrt(0 + 2**4) == 4 != 2 * math.sqrt(0 + 1)  # N_c(0,2) = 4, 2 N_c(0,1) = 2
Nd = lambda a, b: (math.sqrt(abs(a)) + math.sqrt(abs(b))) ** 2
assert Nd(1, 0) == Nd(0, 1) == 1 and Nd(1, 1) == 4 > 2
assert Matrix([[2, 1], [0, 1]]).det() != 0  # v1 : |2x+y|+|y|
assert abs(1 - 1) + abs(2 - 2) == 0 and 2**2 + 0 == 4 != 2 * 1  # v1 : N_a(1,-1) = 0, N_b(2,0) = 4
assert Matrix([[1, 0], [1, -1]]).det() != 0  # v2 : sqrt(x^2 + (x-y)^2)
assert 1 + 1 - abs(1 - (-1)) == 0 and math.sqrt(0 + 4) == 2 != 4 * math.sqrt(1) and max(abs(1 - 1), abs(2 - 2)) == 0  # v2
for _ in range(500):  # positivité de |x|+|y|-|x-y| (v2) et triangle pour N_a
    a, b, c, d = (random.uniform(-3, 3) for _ in range(4))
    assert abs(a) + abs(b) - abs(a - b) >= -1e-12
    assert Na(a + c, b + d) <= Na(a, b) + Na(c, d) + 1e-12

# top-l2-evn-normes-ex3 (+ v1, v2) : constantes sqrt n, sqrt n, n atteintes en (1,...,1) ; leurres distincts pour n >= 2.
for d in range(2, 30):
    un = [1] * d
    assert math.isclose(n1(un), math.sqrt(d) * n2(un)) and math.isclose(n2(un), math.sqrt(d) * ninf(un)) and n1(un) == d * ninf(un)
    assert len({math.sqrt(d), d, 1, d * d}) == 4

# top-l2-evn-normes-ex4 (+ v1, v2) : quotients n+1, (n+1)/sqrt(2n+1), sqrt(2n+1) et leurres distincts pour n >= 1.
for k in range(1, 40):
    I1, I2 = F(1, k + 1), F(1, 2 * k + 1)  # N_1(f_k), N_2(f_k)^2
    assert 1 / I1 == k + 1
    r1 = (k + 1) / math.sqrt(2 * k + 1)
    assert math.isclose(math.sqrt(I2) / I1, r1)
    assert len({k + 1, k, 1 / (k + 1), 2 * k + 1}) == 4
    assert len({round(r1, 9), round(math.sqrt(2 * k + 1), 9), round((k + 1) / (2 * k + 1), 9), round(1 / r1, 9)}) == 4
    assert len({round(math.sqrt(2 * k + 1), 9), 2 * k + 1, round(1 / math.sqrt(2 * k + 1), 9), round(math.sqrt(k + 1), 9)}) == 4
assert limit((n + 1) / sqrt(2 * n + 1), n, oo) == oo

# ---------------------------------------------------------------------------
# top-l2-evn-ouverts-fermes : contre-exemples et exercices (faits d'ensembles vérifiés sur des points tests)
# ---------------------------------------------------------------------------
# Contre-exemple : intersection des ]-1/k, 1/k[ réduite à {0} (tout r > 0 sort d'un intervalle).
assert all(any(not (-F(1, k) < F(1, m) < F(1, k)) for k in range(1, 3 * m)) for m in range(1, 50))
# Équivalence des normes : B_N(a, r/beta) ⊂ B_N'(a, r) avec N = ||.||_inf, N' = ||.||_1, beta = 2 dans R^2.
for _ in range(500):
    v = [random.uniform(-1, 1) for _ in range(2)]
    if ninf(v) < 0.5:
        assert n1(v) < 1
# ex3 : int([0,1] ∪ [1,2]) contient 1 ; v2 : int(]-oo,0]) ∪ int([0,+oo[) = R \ {0}.
assert all(0 <= 1 + h <= 2 for h in (-0.5, 0.5))
# ex4 : x_n = a + (1 - 1/n)(x - a) dans B(a, r) et y = a + (1+s)(x-a) hors de B'(a, r).
for r in (0.5, 1, 3):
    for k in range(1, 20):
        assert (1 - 1 / k) * r < r and (1 + 1 / (k + 1)) * r > r

# ---------------------------------------------------------------------------
# top-l2-evn-suites
# ---------------------------------------------------------------------------
# Exemple : H = {xy = 1} non bornée ; (k, 1/k) dans H.
assert all(F(k) * F(1, k) == 1 for k in range(1, 20))
# top-l2-evn-suites-ex1 (+ v1, v2) : limites des coordonnées.
h = Symbol("h", positive=True)
assert limit((1 + 1 / n) ** n, n, oo) == E and limit(n * sin(1 / n), n, oo) == 1
assert limit(n * (exp(1 / n) - 1), n, oo) == 1 and limit((1 - 2 / n) ** n, n, oo) == exp(-2)
assert limit(n * log(1 + 1 / n), n, oo) == 1 and limit(sqrt(n**2 + n) - n, n, oo) == Rational(1, 2)


# top-l2-evn-suites-ex2 (+ v1, v2) : points ajoutés absents de A, limites des sous-suites.
def elements(f, N=4000):
    return {f(k) for k in range(1, N)}


A0 = elements(lambda k: (-1) ** k + F(1, k))
assert -1 not in A0 and 1 not in A0 and 0 in A0  # 0 = -1 + 1/1 est dans A (d'où le leurre « A » seul)
sinv = lambda k: [0, 1, 0, -1][k % 4]
A1 = elements(lambda k: sinv(k) + F(1, k))
assert all(v not in A1 for v in (-1, 0, 1)) and 2 in A1
cosv = lambda k: F(1) if k % 3 == 0 else F(-1, 2)
A2 = elements(lambda k: cosv(k) + F(1, k))
assert F(-1, 2) not in A2 and 1 not in A2 and min(A2) > F(-1, 2)  # -1 n'est pas adhérent
assert F(1, 2) in A2 and 0 in A2  # pourquoi les leurres A ∪ {-1/2, 1/2, 1} etc. ont été évités

# top-l2-evn-suites-ex3 (+ v1, v2) : 0 hors de D, de D'(v1), 1 hors de D' (v2) ; les points ajoutés ne sont pas dans A.
D = {F(1, k) for k in range(1, 500)}
Dv1 = {F((-1) ** k, k) for k in range(1, 500)}
Dp = {1 - F(1, m) for m in range(1, 500)}
assert 0 not in D and 0 not in Dv1 and 1 not in Dp and 0 in Dp  # 0 = 1 - 1/1 appartient à D'
# Leurre « A ∪ (D×{0}) ∪ ({0}×D) » : il manque (0,0) ; v2 : il manque (0,1).
assert (0, 0) not in {(a, 0) for a in D} | {(0, b) for b in D}
assert (0, 1) not in {(a, 1) for a in D} | {(0, b) for b in Dp}

# ---------------------------------------------------------------------------
# top-l2-evn-continuite-lineaire
# ---------------------------------------------------------------------------
# Contre-exemples : N(X^k) = 1 et N(k X^(k-1)) = k sur [0,1] ; f_k = (1-t)^k : N_1 = 1/(k+1).
assert all(integrate((1 - t) ** k, (t, 0, 1)) == Rational(1, k + 1) for k in range(1, 6))
# Affirmation v∘u = 0 pour u(x,y) = (x,0), v(x,y) = (0,y).
U, V = Matrix([[1, 0], [0, 0]]), Matrix([[0, 0], [0, 1]])
assert V * U == Matrix.zeros(2, 2)
# Trou 3, leurre <= |||v||| + |||u||| : u = v = 3 id : 9 > 6.
assert 9 > 3 + 3


def norme_inf_mat(M):
    """Norme subordonnée à ||.||_inf : max des sommes de lignes, contrôlée par les vecteurs de signes."""
    p = M.shape[1]
    s = max(sum(abs(c) for c in M.row(i)) for i in range(M.shape[0]))
    best = max(max(abs(c) for c in M * Matrix(sg)) for sg in itertools.product([-1, 1], repeat=p))
    assert best == s
    return s


def norme_1_mat(M):
    return max(sum(abs(c) for c in M.col(j)) for j in range(M.shape[1]))


# top-l2-evn-continuite-lineaire-ex1 (+ v1, v2) : 7, 8, 4 ; leurres = colonnes, somme totale, max des coefficients.
for M, val, leurres, xe, Axe in [
    (Matrix([[1, -2], [3, 4]]), 7, [6, 10, 4], [1, 1], [-1, 7]),
    (Matrix([[2, -1, 0], [1, 3, -1], [0, -2, 6]]), 8, [7, 16, 6], [1, -1, 1], [3, -3, 8]),
    (Matrix([[-3, 1], [2, 2]]), 4, [5, 8, 3], [-1, 1], [4, 0]),
]:
    assert norme_inf_mat(M) == val
    assert leurres == [norme_1_mat(M), sum(abs(c) for c in M), max(abs(c) for c in M)] and val not in leurres
    assert list(M * Matrix(xe)) == Axe and max(abs(c) for c in Axe) == val
# top-l2-evn-q4 : norme subordonnée à ||.||_1 de la même matrice = 6, atteinte en (0,1) ; 7 (lignes) est un leurre.
assert norme_1_mat(Matrix([[1, -2], [3, 4]])) == 6 and list(Matrix([[1, -2], [3, 4]]) * Matrix([0, 1])) == [-2, 4]
for _ in range(300):
    v = Matrix([random.uniform(-1, 1), random.uniform(-1, 1)])
    assert n1(list(Matrix([[1, -2], [3, 4]]) * v)) <= 6 * n1(list(v)) + 1e-9


# top-l2-evn-continuite-lineaire-ex2 (+ v1, v2) : valeurs atteintes par les fonctions de norme 1.
def pl(points):
    """Fonction affine par morceaux passant par les points donnés (abscisses croissantes)."""
    def f(s):
        for (a, fa), (b, fb) in zip(points, points[1:]):
            if a < b and a <= s <= b:
                return fa + (fb - fa) * (s - a) / (b - a)
    return f


grille = [F(i, 120) for i in range(121)]
for phi, f, val, leurres in [
    (lambda g: g(0) - 2 * g(1), pl([(0, 1), (1, -1)]), 3, [1, 2, math.sqrt(5)]),
    (lambda g: 3 * g(0) + g(F(1, 2)) - g(1), pl([(0, 1), (F(1, 2), 1), (1, -1)]), 5, [3, 1, math.sqrt(11)]),
    (lambda g: g(0) - g(F(1, 2)) + g(1), pl([(0, 1), (F(1, 2), -1), (1, 1)]), 3, [1, 2, math.sqrt(3)]),
]:
    assert max(abs(f(s)) for s in grille) == 1 and phi(f) == val and val not in leurres
assert pl([(0, 1), (F(1, 2), 1), (1, -1)])(F(3, 4)) == 3 - 4 * F(3, 4)  # f = 3 - 4t sur [1/2, 1]
assert pl([(0, 1), (F(1, 2), -1), (1, 1)])(F(1, 5)) == abs(4 * F(1, 5) - 2) - 1  # f = |4t-2| - 1


# top-l2-evn-continuite-lineaire-ex3 (+ v1, v2) : phi(f_n) = 1 - 1/n, 2 - 1/n, 5/3 - 3/(2n) (intégrales exactes de
# fonctions affines par morceaux : méthode des trapèzes, exacte sur chaque morceau).
def integrale(points, a, b):
    f = pl(points)
    noeuds = sorted({a, b} | {p for p, _ in points if a < p < b})
    return sum((q - p) * (f(p) + f(q)) / 2 for p, q in zip(noeuds, noeuds[1:]))


for k in range(3, 30):
    fk = [(0, 1), (F(1, 2) - F(1, k), 1), (F(1, 2) + F(1, k), -1), (1, -1)]
    assert integrale(fk, 0, F(1, 2)) - integrale(fk, F(1, 2), 1) == 1 - F(1, k)
    gk = [(0, 1), (1 - F(1, k), 1), (1, -1)]
    assert integrale(gk, 0, 1) - (-1) == 2 - F(1, k)
    hk = [(0, 1), (F(1, 3) - F(1, k), 1), (F(1, 3) + F(1, k), -1), (1, -1)]
    assert integrale(hk, 0, F(1, 3)) == F(1, 3) - F(1, 2 * k) and integrale(hk, F(1, 3), 1) == -F(2, 3) + F(1, 2 * k)
    assert integrale(hk, 0, F(1, 3)) - 2 * integrale(hk, F(1, 3), 1) == F(5, 3) - F(3, 2 * k)
assert F(1, 3) + 2 * F(2, 3) == F(5, 3)  # majoration de v2

# top-l2-evn-continuite-lineaire-ex4 : N'(P') = sum_{k>=1} k! |a_k| <= N'(P), égalité pour P = X.
for _ in range(200):
    a = [random.randint(-9, 9) for _ in range(7)]
    Np = sum(math.factorial(k) * abs(c) for k, c in enumerate(a))
    der = [k * a[k] for k in range(1, 7)]
    assert sum(math.factorial(k) * abs(c) for k, c in enumerate(der)) == Np - abs(a[0]) <= Np

# ---------------------------------------------------------------------------
# top-l2-evn-equivalence-dim-finie
# ---------------------------------------------------------------------------
# Contre-exemple : P_n = 1 + X + ... + X^n : N_1 = n+1, N_inf = 1, P_n(1) = n+1.
assert all(n1([1] * (k + 1)) == k + 1 == sum([1] * (k + 1)) and ninf([1] * (k + 1)) == 1 for k in range(10))


# top-l2-evn-equivalence-dim-finie-ex1 (+ v1, v2) : (alpha, beta) = (1,3), (1,4), (1/2,4).
# N est convexe et affine par morceaux : extremums sur le bord du carré unité aux sommets et aux points de cassure ;
# contrôle sur une grille rationnelle fine du bord (qui contient ces points).
def bord(m=240):
    pts = set()
    for i in range(-m, m + 1):
        s = F(i, m)
        pts |= {(s, F(1)), (s, F(-1)), (F(1), s), (F(-1), s)}
    return pts


B = bord()
for N, alpha, beta, leurres in [
    (lambda a, b: abs(a) + abs(a - b), 1, 3, [(1, 2), (F(1, 2), 3), (0, 3)]),
    (lambda a, b: 2 * abs(a) + abs(a + b), 1, 4, [(1, 3), (2, 4), (F(1, 2), 4)]),
    (lambda a, b: abs(a + 2 * b) + abs(b), F(1, 2), 4, [(1, 4), (F(1, 2), 3), (F(1, 3), 4)]),
]:
    vals = [N(a, b) for a, b in B]
    assert min(vals) == alpha and max(vals) == beta and (alpha, beta) not in leurres
assert abs(0) + abs(0 - 1) == 1 and abs(1) + abs(1 + 1) == 3  # égalités du corrigé (ex1)
assert abs(1 + 2 * F(-1, 2)) + abs(F(-1, 2)) == F(1, 2)  # v2 : égalité en (1, -1/2)

# top-l2-evn-equivalence-dim-finie-ex2 : det(A - tI) de coefficient dominant (-1)^n.
for d in range(1, 5):
    M = Matrix(d, d, lambda i, j: random.randint(-3, 3))
    P = (M - t * Matrix.eye(d)).det()
    assert expand(P).coeff(t, d) == (-1) ** d

# top-l2-evn-equivalence-dim-finie-ex3 : Taylor-Lagrange, |e^t - P_n(t)| <= e/(n+1)! sur [0,1].
for k in range(0, 12):
    Pk = sum(F(1, math.factorial(j)) for j in range(k + 1))
    assert abs(math.e - float(Pk)) <= math.e / math.factorial(k + 1) + 1e-15

# top-l2-evn-equivalence-dim-finie-ex4 (+ v1, v2) : constantes optimales max_i sup |L_i|.
s = Symbol("s", real=True)


def sup_abs(p, a, b):
    crit = [c for c in solve(diff(p, s), s) if c.is_real and a <= c <= b]
    return max(abs(p.subs(s, c)) for c in crit + [a, b])


L0, L1, L2 = (s - 1) * (s - 2) / 2, -s * (s - 2), s * (s - 1) / 2
assert [sup_abs(L, 0, 2) for L in (L0, L1, L2)] == [1, 1, 1]  # ex4 : C = 1
assert [sup_abs(L, 0, 3) for L in (L0, L1, L2)] == [1, 3, 3]  # v1 : C = 3
assert [sup_abs(L, 0, 1) for L in (1 - s, s)] == [1, 1]  # v2 : C = 1
assert sup_abs(L0, 1, 2) == Rational(1, 8) and sup_abs(L2, 0, 1) == Rational(1, 8)
assert expand(L0 + L1 + L2) == 1 and [L1.subs(s, v) for v in (0, 1, 2)] == [0, 1, 0]
q3 = s * (s - 1) * (s - 2)
assert [q3.subs(s, v) for v in (0, 1, 2)] == [0, 0, 0]  # X(X-1)(X-2) : N n'est pas une norme pour n >= 3
# Le leurre C = 1/3 (ex4) : P = L0 donne sup|P| = 1 = N(P), donc C >= 1.
assert sum(abs(L0.subs(s, v)) for v in (0, 1, 2)) == 1

# top-l2-evn-q7 : A^k -> 0 ; |||A|||_inf = 3/4 ; valeurs propres 3/4 et 1/4.
A = Matrix([[Rational(1, 2), Rational(1, 4)], [Rational(1, 4), Rational(1, 2)]])
assert norme_inf_mat(A) == Rational(3, 4) and set(A.eigenvals()) == {Rational(3, 4), Rational(1, 4)}
assert max(abs(c) for c in A**40) < Rational(1, 10**4)
# top-l2-evn-q1 : x = (1,-1,1,-1) : ||x||_1 / ||x||_2 = 4/2 = 2.
assert n1([1, -1, 1, -1]) / n2([1, -1, 1, -1]) == 2
# top-l2-evn-q5 : X^k vaut 2^k en 2 pour une norme sup sur [0,1] égale à 1.
assert all(2**k >= k for k in range(30))

# ---------------------------------------------------------------------------
# top-l2-evn-compacts
# ---------------------------------------------------------------------------
# top-l2-evn-compacts-ex2 (+ v1, v2) : min de x^2 + (x^2 - h)^2 = h - 1/4 en x^2 = (2h-1)/2 ; d = sqrt(h - 1/4).
for hh, x2, fmin, dist, leurres in [
    (2, Rational(3, 2), Rational(7, 4), sqrt(7) / 2, [2, Rational(7, 4), sqrt(Rational(3, 2))]),
    (3, Rational(5, 2), Rational(11, 4), sqrt(11) / 2, [3, Rational(11, 4), sqrt(Rational(5, 2))]),
    (1, Rational(1, 2), Rational(3, 4), sqrt(3) / 2, [1, Rational(3, 4), sqrt(2) / 2]),
]:
    g = x**2 + (x**2 - hh) ** 2
    assert expand(diff(g, x) - 2 * x * (2 * x**2 - (2 * hh - 1))) == 0
    crit = solve(diff(g, x), x)
    assert min(g.subs(x, c) for c in crit) == fmin and g.subs(x, sqrt(x2)) == fmin and fmin < hh**2
    assert sqrt(fmin) == dist and dist < hh and all(sqrt(fmin) != l for l in leurres)
    assert len({float(dist)} | {float(l) for l in leurres}) == 4

# top-l2-evn-compacts-ex3 : tr(A^T A) = n pour A orthogonale ; diag(k, 1/k) dans SL_2 ; [[1,k],[0,0]] projecteur.
R = Matrix([[cos(1), -sin(1)], [sin(1), cos(1)]])
assert (R.T * R).applyfunc(lambda c: c.simplify()) == Matrix.eye(2) and (R.T * R).trace().simplify() == 2
for k in range(1, 6):
    assert Matrix([[k, 0], [0, Rational(1, k)]]).det() == 1
    Bk = Matrix([[1, k], [0, 0]])
    assert Bk**2 == Bk

# top-l2-evn-compacts-ex4 : h(x) = sqrt(x^2+1) strictement 1-lipschitzienne, sans point fixe.
hf = lambda u: math.sqrt(u * u + 1)
for _ in range(2000):
    a, b = random.uniform(-50, 50), random.uniform(-50, 50)
    if a != b:
        assert abs(hf(a) - hf(b)) < abs(a - b)
assert solve(sqrt(x**2 + 1) - x, x) == []

# top-l2-evn-q8 : x^4 + y^4 = 1 entraîne |x| <= 1 et |y| <= 1.
assert all(abs(c) <= 1 for c in (1, 0, 2 ** -0.25))

print("top_l2_evn : toutes les vérifications passent.")
