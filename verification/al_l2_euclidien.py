"""Vérifications calculatoires du chapitre « al-l2-euclidien » (algèbre linéaire, L2) : cartes, exercices et QCM
de content/*/algebre-lineaire/l2-al-l2-euclidien.json et content/qcm/al-l2-euclidien.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/al_l2_euclidien.py (doit finir sans erreur).
"""

from fractions import Fraction as F
from itertools import product

from sympy import (I, Matrix, Rational, cos, diag, eye, integrate, nsimplify, pi, simplify, solve, sqrt, symbols,
                   zeros)

R = Rational
t = symbols("t", real=True)


def ps(a, b):
    return (Matrix(a).T * Matrix(b))[0]


def nrm(a):
    return sqrt(ps(a, a))


def egal(a, b):
    return simplify(Matrix(a) - Matrix(b)) == zeros(*Matrix(a).shape)


def schmidt(vs, prod=ps):
    fs, gs = [], []
    for v in vs:
        v = Matrix(v)
        g = v - sum((prod(v, f) * f for f in fs), zeros(*v.shape))
        g = g.applyfunc(simplify)
        gs.append(g)
        fs.append((g / sqrt(prod(g, g))).applyfunc(simplify))
    return fs, gs


def refl(a):
    a = Matrix(a)
    return eye(len(a)) - 2 * a * a.T / ps(a, a)


# =========================================================================================================
# al-l2-euclidien-cauchy-schwarz (carte)
# =========================================================================================================
# Exemple : (n(n+1)/2)^2 <= n * n(n+1)(2n+1)/6, égalité seulement pour n = 1.
for n in range(1, 60):
    g, d = F(n * (n + 1), 2) ** 2, n * F(n * (n + 1) * (2 * n + 1), 6)
    assert g <= d and ((g == d) == (n == 1)) and ((3 * (n + 1) <= 2 * (2 * n + 1)) == True)
# Contre-exemple « sans symétrie » : b(x,y) = x1y1 + x2y2 + 2(x1y2 - x2y1), b(x,x) = |x|^2, b(e1,e2) = 2 > 1.
b = lambda x, y: x[0] * y[0] + x[1] * y[1] + 2 * (x[0] * y[1] - x[1] * y[0])
assert b((3, -7), (3, -7)) == 58 and b((1, 0), (0, 1)) == 2 and b((1, 0), (1, 0)) == 1 == b((0, 1), (0, 1))
# Contre-exemple « sans positivité » : b(x,y) = x1y1 - x2y2, x = (2,1), y = (1,0) : 2 > sqrt(3).
m = lambda x, y: x[0] * y[0] - x[1] * y[1]
assert m((2, 1), (2, 1)) == 3 and m((1, 0), (1, 0)) == 1 and m((2, 1), (1, 0)) == 2 and 2 > sqrt(3)
# Piège « non dégénérée » : la même forme est non dégénérée (matrice diag(1,-1) inversible).
assert diag(1, -1).det() != 0
# Contre-exemple « sans caractère défini » : b(x,y) = x1y1, x = (0,1), y = (1,0) : égalité 0 = 0 avec famille libre.
assert 0 * 1 == 0 and Matrix([[0, 1], [1, 0]]).det() != 0
# Contre-exemples des normes 1 et infini (parallélogramme).
assert (1 + 1) ** 2 + (1 + 1) ** 2 == 8 and 2 * 1 + 2 * 1 == 4          # norme 1
assert 1 ** 2 + 1 ** 2 == 2 and 2 * 1 + 2 * 1 == 4                      # norme infini (affirmation 3)
# Trou 1 : les leurres sont faux (x = y = 2e1 ; x = e1, y = e2).
x, y = Matrix([2, 0]), Matrix([2, 0])
assert not (abs(ps(x, y)) <= min(nrm(x), nrm(y)))                       # 4 <= 2 faux
assert not (ps(x, y) ** 2 <= ps(x, x) + ps(y, y))                       # 16 <= 8 faux
assert not (nrm((1, 0)) * nrm((0, 1)) <= abs(ps((1, 0), (0, 1))))       # 1 <= 0 faux
# Trou 2 : orthogonaux non nuls -> inégalité stricte ; x = e1, y = -e1 : égalité sans y = lambda x, lambda >= 0.
assert abs(ps((1, 0), (0, 1))) < 1 and abs(ps((1, 0), (-1, 0))) == 1
# Trou 3 et affirmation 2 : y = -x : famille liée mais |x + y| = 0 < 2|x| ; e1, e2 : |e1 + e2| = sqrt 2 < 2.
assert nrm((0, 0)) == 0 and nrm((1, 1)) < 2

# =========================================================================================================
# al-l2-euclidien-cauchy-schwarz-ex1 (+ v1, v2) : max de <u, X> sur la sphère = |u|.
for u, N, s in [((2, 3, 6), 7, 11), ((1, 4, 8), 9, 13), ((2, 6, 9), 11, 17)]:
    assert nrm(u) == N and sum(u) == s
    assert ps(u, Matrix(u) / N) == N                                     # atteint en u/|u|
    assert len({N, s, N ** 2, sqrt(N)}) == 4                              # propositions distinctes, leurres faux
# al-l2-euclidien-cauchy-schwarz-ex2 (+ v1, v2) : angles pi/4, pi/3, 2pi/3.
for A, B, th in [([[1, 0], [0, 1]], [[1, 1], [1, 1]], pi / 4), ([[1, 1], [0, 0]], [[0, 1], [1, 0]], pi / 3),
                 ([[1, 1], [0, 0]], [[0, -1], [1, 0]], 2 * pi / 3)]:
    A, B = Matrix(A), Matrix(B)
    assert (A.T * B).trace() == sum(A[k] * B[k] for k in range(4))
    c = (A.T * B).trace() / sqrt((A.T * A).trace() * (B.T * B).trace())
    assert simplify(c - cos(th)) == 0
    for autre in [pi / 6, pi / 4, pi / 3, 2 * pi / 3]:
        if autre != th:
            assert simplify(c - cos(autre)) != 0
# al-l2-euclidien-cauchy-schwarz-ex3 : min de sum 1/x_k = n^2 ; leurres n, 1/n, n(n+1)/2 différents pour n >= 2.
for n in range(2, 30):
    assert n * n not in (n, F(1, n), F(n * (n + 1), 2))
    xs = [F(1, n)] * n
    assert sum(xs) == 1 and sum(1 / v for v in xs) == n * n
    xs = [F(1, 2 * n)] * (n - 1) + [1 - F(n - 1, 2 * n)]                  # un autre point : valeur > n^2
    assert sum(xs) == 1 and sum(1 / v for v in xs) > n * n
# al-l2-euclidien-cauchy-schwarz-ex4 : I + N, N nilpotente non nulle : tr = n, tr(t(I+N)(I+N)) = n + |N|^2 > n.
for n in (2, 3):
    N = zeros(n, n)
    N[0, n - 1] = 1
    M = eye(n) + N
    assert M.trace() == n and (M.T * M).trace() == n + 1
    O = Matrix([[0, -1], [1, 0]]) if n == 2 else None                     # une orthogonale de trace != n
    if O is not None:
        assert O.T * O == eye(2) and O.trace() != 2
    S = diag(n, *([0] * (n - 1)))                                         # symétrique de trace n, tr(S^2) = n^2 != n
    assert S.trace() == n and (S.T * S).trace() != n
# Question 3 : projecteur orthogonal diag(1,0) : tr = 1 dans [0, 2].
P = diag(1, 0)
assert P * P == P and P.T == P and 0 <= P.trace() <= 2

# =========================================================================================================
# al-l2-euclidien-schmidt (carte)
# =========================================================================================================
(f1, f2, f3), (g1, g2, g3) = schmidt([(1, 1, 0), (1, 0, 1), (0, 1, 1)])
assert egal(g2, [R(1, 2), -R(1, 2), 1]) and egal(f2, Matrix([1, -1, 2]) / sqrt(6))
assert egal(g3, [-R(2, 3), R(2, 3), R(2, 3)]) and egal(f3, Matrix([-1, 1, 1]) / sqrt(3))
assert ps((0, 1, 1), f1) == 1 / sqrt(2) and simplify(ps((0, 1, 1), f2) - 1 / sqrt(6)) == 0
# Contre-exemples : famille liée (g2 = 0) ; ordre ; rotation d'angle t : <f1,e1> = <f2,e2> = cos t > 0.
assert schmidt([(1, 0), (1, 1)])[0] == [Matrix([1, 0]), Matrix([0, 1])]
assert egal(schmidt([(1, 1), (1, 0)])[0][1], Matrix([1, -1]) / sqrt(2))
assert Matrix([2, 0]) - 2 * Matrix([1, 0]) == zeros(2, 1)
# Affirmation 2 : e2 = (1,1)/sqrt2, e2 - <e2,e1>e1 = (0, 1/sqrt2) n'est pas normé.
assert nrm(Matrix([1, 1]) / sqrt(2) - Matrix([1, 0]) / sqrt(2)) == 1 / sqrt(2)
# Affirmation 4 : base (1,0),(1,1), x = y = e2 : coordonnées (0,1), somme 1, mais |e2|^2 = 2.
assert ps((1, 1), (1, 1)) == 2
# Trou 3, leurre <f_k, e_k> = |e_k| : impose f_k = e_k/|e_k| ; impossible si e1, e2 non orthogonaux.
assert ps((1, 0), (1, 1)) != 0
# Trou 2, leurre « de vecteurs unitaires » : (1,0) puis f2 = (cos s, sin s) avec sin s > 0 hors orthogonalité.
assert ps((1, 0), (R(3, 5), R(4, 5))) != 0 and ps((R(3, 5), R(4, 5)), (0, 1)) > 0

# al-l2-euclidien-schmidt-ex1 (+ v1, v2) : f3 et leurres (-f3, g3, normalisé de e3 - <e3,f1>f1).
for fam, f3_att in [([(1, 1, 1), (0, 1, 1), (0, 0, 1)], Matrix([0, -1, 1]) / sqrt(2)),
                    ([(1, 2, 2), (1, 0, 1), (0, 0, 1)], Matrix([-2, -1, 2]) / 3),
                    ([(0, 1, 1), (1, 0, 1), (1, 1, 1)], Matrix([1, 1, -1]) / sqrt(3))]:
    assert Matrix(fam).det() != 0
    (a1, a2, a3), (_, _, h3) = schmidt(fam)
    assert egal(a3, f3_att) and ps(a3, fam[2]) > 0
    h = Matrix(fam[2]) - ps(fam[2], a1) * a1
    h = h / nrm(h)
    for L in (-a3, h3, h):
        assert not egal(L, a3)
# v0 : g3 = (0,-1/2,1/2), leurre h = (-1,-1,2)/sqrt6.
assert egal(schmidt([(1, 1, 1), (0, 1, 1), (0, 0, 1)])[1][2], [0, -R(1, 2), R(1, 2)])

# al-l2-euclidien-schmidt-ex2 (+ v1, v2) : Schmidt sur (1, X, X^2) pour trois produits scalaires.
X = symbols("X")
P2 = lambda p, q, lo, hi: integrate(p * q, (X, lo, hi))
for prod, f3_att, g3_att, leurre_norm in [
    (lambda p, q: P2(p, q, -1, 1), sqrt(10) / 4 * (3 * X**2 - 1), X**2 - R(1, 3), sqrt(10) / 2 * X**2),
    (lambda p, q: P2(p, q, 0, 1), sqrt(5) * (6 * X**2 - 6 * X + 1), X**2 - X + R(1, 6), sqrt(5) * X**2),
    (lambda p, q: sum((p * q).subs(X, v) for v in (-1, 0, 1)), (3 * X**2 - 2) / sqrt(6), X**2 - R(2, 3), X**2 / sqrt(2)),
]:
    fs, gs_ = [], []
    for v in (1 + 0 * X, X, X**2):
        g = simplify(v - sum(prod(v, f) * f for f in fs))
        fs.append(simplify(g / sqrt(prod(g, g))))
        gs_.append(g)
    assert simplify(fs[2] - f3_att) == 0 and simplify(gs_[2] - g3_att) == 0
    assert simplify(prod(X**2, X**2) - (leurre_norm / X**2) ** -2) == 0       # X^2 normé = leurre 4
    for L in (g3_att, -f3_att, leurre_norm):
        assert simplify(L - f3_att) != 0
assert integrate((t**2 - R(1, 3)) ** 2, (t, -1, 1)) == R(8, 45)
assert integrate((t**2 - t + R(1, 6)) ** 2, (t, 0, 1)) == R(1, 180)
assert integrate((t - R(1, 2)) ** 2, (t, 0, 1)) == R(1, 12)
# v1 : leurre sqrt5 (6X^2 - 6X - 1) différent.
assert simplify(sqrt(5) * (6 * X**2 - 6 * X - 1) - sqrt(5) * (6 * X**2 - 6 * X + 1)) != 0

# al-l2-euclidien-schmidt-ex3 (+ v1, v2) : A = QR.
for A, Ratt in [([[3, 1], [4, 2]], [[5, R(11, 5)], [0, R(2, 5)]]), ([[1, 3], [1, 1]], [[sqrt(2), 2 * sqrt(2)], [0, sqrt(2)]]),
                ([[1, 7], [2, 4]], [[sqrt(5), 3 * sqrt(5)], [0, 2 * sqrt(5)]])]:
    A, Ratt = Matrix(A), Matrix(Ratt)
    Q = A * Ratt.inv()
    assert egal(Q.T * Q, eye(2)) and Ratt[0, 0] > 0 and Ratt[1, 1] > 0
    c2 = A[:, 1]
    for L in (Matrix([[Ratt[0, 0], -Ratt[0, 1]], [0, Ratt[1, 1]]]), Matrix([[Ratt[0, 0], Ratt[0, 1]], [0, nrm(c2)]]),
              Matrix([[Ratt[0, 0], 0], [Ratt[0, 1], Ratt[1, 1]]])):
        QL = A * L.inv()                                                    # avec un leurre, Q ne serait pas orthogonale
        assert not egal(QL.T * QL, eye(2))
# al-l2-euclidien-schmidt-ex4 : Hadamard, n = 2 et n = 4 ; leurres n!, n^n, n distincts de n^(n/2) pour n >= 2... (sauf
# n = 1, hors sujet) : pour n = 4, n! = 24 > 16 n'est pas la borne de Hadamard.
H2 = Matrix([[1, 1], [1, -1]])
H4 = Matrix([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]])
assert abs(H2.det()) == 2 and abs(H4.det()) == 16 and H4.T * H4 == 4 * eye(4)
# Leurres : pour n = 4, la borne de Hadamard 4^2 = 16 diffère de 4! = 24, 4^4 = 256 et 4.
assert 4 ** 2 == 16 and len({16, 24, 256, 4}) == 4
for n in (3, 4, 5):
    assert simplify(sqrt(n) ** n - n ** R(n, 2)) == 0 and len({n ** R(n, 2), sqrt(n) ** n}) == 1
# Test aléatoire de Hadamard sur des matrices 3x3 entières.
import random
random.seed(1)
for _ in range(200):
    M = Matrix(3, 3, lambda i, j: random.randint(-3, 3))
    assert abs(M.det()) ** 2 <= (M.T * M)[0, 0] * (M.T * M)[1, 1] * (M.T * M)[2, 2]

# =========================================================================================================
# al-l2-euclidien-projection (carte)
# =========================================================================================================
# Exemple : moindres carrés (0,1), (1,1), (2,3) : a = 1, c = 2/3, somme des carrés 2/3.
A = Matrix([[0, 1], [1, 1], [2, 1]])
B = Matrix([1, 1, 3])
assert A.T * A == Matrix([[5, 3], [3, 3]]) and A.T * B == Matrix([7, 5])
s = (A.T * A).inv() * A.T * B
assert s == Matrix([1, R(2, 3)])
r = A * s - B
assert r == Matrix([-R(1, 3), R(2, 3), -R(1, 3)]) and (r.T * r)[0] == R(2, 3) and A.T * r == zeros(2, 1)
# Contre-exemples : base non orthonormale ((2,1) != (1,0)) ; base orthogonale non normée (facteur 4).
assert ps((1, 0), (1, 0)) * Matrix([1, 0]) + ps((1, 0), (1, 1)) * Matrix([1, 1]) == Matrix([2, 1])
assert ps((1, 0), (2, 0)) * Matrix([2, 0]) == Matrix([4, 0])
# Moindres carrés sans rang plein : A = [[1,1],[1,1]], b = (1,0) : tA A X = tA b <=> s + t = 1/2.
A = Matrix([[1, 1], [1, 1]])
s_, t_ = symbols("s t")
assert solve(list(A.T * A * Matrix([s_, t_]) - A.T * Matrix([1, 0])), [s_, t_]) == {s_: R(1, 2) - t_}
# Projecteur non orthogonal p(x,y) = (x - y, 0).
Pm = Matrix([[1, -1], [0, 0]])
assert Pm * Pm == Pm and ps((1, 0), (1, 1)) != 0 and Pm * Matrix([0, 1]) == Matrix([-1, 0])
# Affirmation « det G(x1,x2,x3) = 0 » : (1,0),(0,1),(1,1).
V = Matrix([[1, 0, 1], [0, 1, 1]])
assert (V.T * V).det() == 0
# Gram : det G = det(M)^2 pour une famille libre (exemple), et formule de la distance par Gram.
M3 = Matrix([[1, 2, 0], [0, 1, 1], [1, 0, 3]])
assert (M3.T * M3).det() == M3.det() ** 2
xb = [Matrix([1, 1, 1, 1]), Matrix([1, -1, 1, -1])]
xx = Matrix([1, 2, 3, 4])
G2 = Matrix(2, 2, lambda i, j: ps(xb[i], xb[j]))
G3 = Matrix(3, 3, lambda i, j: ps((xb + [xx])[i], (xb + [xx])[j]))
assert G3.det() / G2.det() == 4                                             # = d(x,F)^2 de projection-ex1

# al-l2-euclidien-projection-ex1 (+ v1, v2).
U1, U2 = Matrix([1, 1, 1, 1]), Matrix([1, -1, 1, -1])
for x, p_att, d in [((1, 2, 3, 4), (2, 3, 2, 3), 2), ((3, 1, 0, 2), (R(3, 2),) * 4, sqrt(5)),
                    ((2, 0, 1, 5), (R(3, 2), R(5, 2), R(3, 2), R(5, 2)), sqrt(13))]:
    x = Matrix(x)
    p = ps(x, U1) / 4 * U1 + ps(x, U2) / 4 * U2
    assert p == Matrix(p_att) and ps(x - p, U1) == 0 == ps(x - p, U2) and nrm(x - p) == d
    assert len({d, d**2, nrm(p), nrm(x)}) == 4

# al-l2-euclidien-projection-ex2 (+ v1, v2).
a_, b_ = symbols("a b")
for f, lo, hi, ab, mn, leurres in [(t**2, 0, 1, (1, -R(1, 6)), R(1, 180), [R(1, 5), R(4, 45), sqrt(5) / 30]),
                                   (t**2, -1, 1, (0, R(1, 3)), R(8, 45), [R(2, 5), R(4, 45), 2 * sqrt(10) / 15]),
                                   (t**3, 0, 1, (R(9, 10), -R(1, 5)), R(9, 700), [R(1, 7), R(9, 112), 3 * sqrt(7) / 70])]:
    J = integrate((f - a_ * t - b_) ** 2, (t, lo, hi))
    sol = solve([J.diff(a_), J.diff(b_)], [a_, b_])
    assert (sol[a_], sol[b_]) == ab and simplify(J.subs(sol) - mn) == 0
    assert leurres[0] == integrate(f**2, (t, lo, hi))                       # |f|^2 (valeur en a = b = 0)
    assert all(simplify(L - mn) != 0 for L in leurres)
assert sqrt(R(1, 180)) == sqrt(5) / 30 and simplify(sqrt(R(8, 45)) - 2 * sqrt(10) / 15) == 0
assert simplify(sqrt(R(9, 700)) - 3 * sqrt(7) / 70) == 0
assert R(1, 5) - R(1, 3) ** 2 == R(4, 45) and R(1, 7) - R(1, 4) ** 2 == R(9, 112)   # approximation par constantes

# al-l2-euclidien-projection-ex3 (+ v1, v2) : droites des moindres carrés et leurres.
for pts, ab, leurres in [([(0, 0), (1, 1), (2, 3), (3, 4)], (R(7, 5), -R(1, 10)), [(-R(1, 10), R(7, 5)), (R(4, 3), 0), (R(7, 5), 2)]),
                         ([(1, 1), (2, 3), (3, 4), (4, 4)], (1, R(1, 2)), [(R(1, 2), 1), (1, 0), (1, 3)]),
                         ([(0, 2), (1, 1), (2, 1), (3, 0)], (-R(3, 5), R(19, 10)), [(R(19, 10), -R(3, 5)), (-R(2, 3), 2), (-R(3, 5), 1)])]:
    A = Matrix([[p[0], 1] for p in pts])
    Y = Matrix([p[1] for p in pts])
    sol = (A.T * A).inv() * A.T * Y
    assert tuple(sol) == ab
    E = lambda a, b: sum((a * p[0] + b - p[1]) ** 2 for p in pts)
    for L in leurres:
        assert E(*L) > E(*ab)

# al-l2-euclidien-projection-ex4 (+ v1, v2).
for M, dS, dH in [([[1, 2], [0, 3]], sqrt(2), 2 * sqrt(2)), ([[2, -1], [3, 0]], 2 * sqrt(2), sqrt(2)),
                  ([[4, 1], [-2, 1]], 3 * sqrt(2) / 2, 5 * sqrt(2) / 2)]:
    M = Matrix(M)
    An = (M - M.T) / 2
    assert simplify(sqrt((An.T * An).trace()) - dS) == 0 and simplify(abs(M.trace()) / sqrt(2) - dH) == 0
    Sym = (M + M.T) / 2
    assert ((M - Sym).T * Sym).trace() == 0                                 # M - p(M) orthogonal à p(M)
    paires = [(dS, dH), (dH, dS), (2 * dS, dH), (dS, abs(M.trace()))]
    assert len({(simplify(p), simplify(q)) for p, q in paires}) == 4

# =========================================================================================================
# al-l2-euclidien-adjoint (carte)
# =========================================================================================================
Mu = Matrix([[1, 2], [0, 3]])
assert Mu.T == Matrix([[1, 0], [2, 3]])                                     # u*(x,y) = (x, 2x + 3y)
# Exemple matriciel : <AM, N> = <M, tA N> (test sur des matrices entières).
random.seed(2)
for _ in range(20):
    A_, M_, N_ = (Matrix(2, 2, lambda i, j: random.randint(-3, 3)) for _ in range(3))
    assert ((A_ * M_).T * N_).trace() == (M_.T * (A_.T * N_)).trace()
# Contre-exemple : p(x,y) = (x,0) dans la base (1,0),(1,1) : matrice [[1,1],[0,0]].
Pb = Matrix([[1, 1], [0, 1]])
assert Pb.inv() * Matrix([[1, 0], [0, 0]]) * Pb == Matrix([[1, 1], [0, 0]])
# u(x,y) = (x + y, y) : F = Vect(e1) stable, u(e2) = (1,1) hors de F^perp.
assert Matrix([[1, 1], [0, 1]]) * Matrix([0, 1]) == Matrix([1, 1])
# Trou 2, leurre « orthogonale » : base (1,0),(0,2), G = diag(1,4), M = E12 : G^-1 tM G != tM.
G_, E12 = diag(1, 4), Matrix([[0, 1], [0, 0]])
assert G_.inv() * E12.T * G_ != E12.T
# Trou 3 et affirmation 2 : u = E12, u* = E21 : noyaux Vect(e1) et Vect(e2), Im u = Vect(e1).
assert E12 * Matrix([1, 0]) == zeros(2, 1) and E12.T * Matrix([0, 1]) == zeros(2, 1)
# Projection orthogonale sur une droite : Ker u* = Ker u != (Ker u)^perp.
# Affirmation 4 : (E12 E21)* = E11, E21 E12 = E22.
assert E12 * E12.T == diag(1, 0) and E12.T * E12 == diag(0, 1)

# al-l2-euclidien-adjoint-ex1 (+ v1, v2) : N = G^-1 tM G et leurres.
for G, M, N_att in [([[2, 0], [0, 1]], [[1, 2], [0, 1]], [[1, 0], [4, 1]]), ([[1, 0], [0, 3]], [[2, 1], [1, 0]], [[2, 3], [R(1, 3), 0]]),
                    ([[2, 1], [1, 1]], [[0, 1], [0, 0]], [[-2, -1], [4, 2]])]:
    G, M, N = Matrix(G), Matrix(M), Matrix(N_att)
    assert G.det() > 0 and G.trace() > 0 and G.inv() * M.T * G == N
    xs_, ys_ = symbols("x1 x2"), symbols("y1 y2")
    xv, yv = Matrix(xs_), Matrix(ys_)
    assert simplify(((M * xv).T * G * yv)[0] - (xv.T * G * (N * yv))[0]) == 0
    for L in (M.T, G * M.T * G.inv(), M.T * G):
        assert L != N
# al-l2-euclidien-adjoint-ex2 : phi*(N) = tA N tB ; les leurres échouent sur un exemple.
A_, B_ = Matrix([[1, 2], [0, 1]]), Matrix([[0, 1], [3, 1]])
M_, N_ = Matrix([[1, 0], [2, 1]]), Matrix([[2, 1], [1, 0]])
ip = lambda P, Q: (P.T * Q).trace()
g = ip(A_ * M_ * B_, N_)
assert g == ip(M_, A_.T * N_ * B_.T)
for L in (B_.T * N_ * A_.T, B_ * N_ * A_, B_.T * N_.T * A_.T):
    assert ip(M_, L) != g
assert ip(M_.T, N_) == ip(M_, N_.T)                                         # transposition autoadjointe
# al-l2-euclidien-adjoint-ex3 : A = (1, i)^T : tA A = 0, rang 1.
Ac = Matrix([1, I])
assert simplify((Ac.T * Ac)[0]) == 0 and Ac.rank() == 1
# al-l2-euclidien-adjoint-ex4 : rangs 2 et 4 atteints par des blocs J dans R^5 ; antisymétrique d'ordre impair singulière.
J = Matrix([[0, -1], [1, 0]])
K2 = diag(J, zeros(3, 3))
K4 = diag(J, J, zeros(1, 1))
assert K2.T == -K2 and K4.T == -K4 and K2.rank() == 2 and K4.rank() == 4
for _ in range(50):
    Z = Matrix(5, 5, lambda i, j: random.randint(-3, 3))
    Z = Z - Z.T
    assert Z.det() == 0 and Z.rank() % 2 == 0

# =========================================================================================================
# al-l2-euclidien-orthogonal (carte)
# =========================================================================================================
Sa = eye(3) - R(2, 3) * Matrix.ones(3, 3)
assert Sa == refl((1, 1, 1)) and Sa.T * Sa == eye(3) and Sa.det() == -1
assert Sa * Matrix([3, 0, 0]) == Matrix([1, -2, -2]) and nrm((1, -2, -2)) == 3
# Contre-exemples : u(x) = |x| e non linéaire ; 2 id ; transvection ; échange de (1,0) et (1,1).
T_ = Matrix([[1, 1], [0, 1]])
assert T_.det() == 1 and nrm(T_ * Matrix([0, 1])) == sqrt(2)
assert nrm((1, 1)) == sqrt(2)
# Trou 1 : s(x,y) = (x, x - y) vérifie s∘s = id mais n'est pas orthogonale.
Sx = Matrix([[1, 0], [1, -1]])
assert Sx * Sx == eye(2) and Sx.T * Sx != eye(2)
# Trou 2 : rotation d'angle pi/2 dans la base (1,0),(0,2) : matrice non orthogonale.
Pb = diag(1, 2)
Mr = Pb.inv() * J * Pb
assert Mr == Matrix([[0, -2], [R(1, 2), 0]]) and Mr.T * Mr != eye(2)
# Affirmations : réflexion de déterminant -1, composée de deux réflexions de déterminant 1.
assert refl((1, 2, 0)).det() == -1 and (refl((1, 2, 0)) * refl((0, 1, 1))).det() == 1

# al-l2-euclidien-orthogonal-ex1 (+ v1, v2).
for k, a, b, c in [(3, (1, 2, 2), (2, 1, -2), (-2, 2, -1)), (7, (2, 3, 6), (3, -6, 2), (6, 2, -3)), (9, (1, 4, 8), (4, 7, -4), (-8, 4, -1))]:
    c1, c2, c3 = Matrix(a) / k, Matrix(b) / k, Matrix(c) / k
    Mo = Matrix.hstack(c1, c2, c3)
    assert Mo.T * Mo == eye(3) and Mo.det() == 1
    for L in (-c3, Matrix([c3[1], c3[0], c3[2]]), Matrix([c3[0], c3[1], -c3[2]])):
        ML = Matrix.hstack(c1, c2, L)
        assert not (ML.T * ML == eye(3) and ML.det() == 1)
# al-l2-euclidien-orthogonal-ex2 (+ v1, v2).
for x, y, z, sz in [((3, 4, 0), (0, 0, 5), (1, 0, 0), (R(16, 25), -R(12, 25), R(3, 5))),
                    ((1, 2, 2), (3, 0, 0), (0, 1, 0), (R(2, 3), R(1, 3), -R(2, 3))),
                    ((2, -1, 2), (0, 3, 0), (0, 0, 1), (-R(1, 3), R(2, 3), R(2, 3)))]:
    a = Matrix(x) - Matrix(y)
    S_ = refl(a)
    assert S_ * Matrix(x) == Matrix(y) and S_ * Matrix(z) == Matrix(sz)
    z = Matrix(z)
    for L in (z - ps(z, a) / ps(a, a) * a, z + 2 * ps(z, a) / ps(a, a) * a, z - 2 * ps(z, a) * a):
        assert L != Matrix(sz)
# al-l2-euclidien-orthogonal-ex3 (+ v1, v2) : décompositions en réflexions et leurres.
e = [None, Matrix([1, 0, 0]), Matrix([0, 1, 0]), Matrix([0, 0, 1])]
for M, bon, leurres in [
    ([[0, 0, 1], [1, 0, 0], [0, 1, 0]], [e[2] - e[1], e[3] - e[2]], [[e[3] - e[2], e[2] - e[1]], [e[2] - e[1], e[3] - e[1]], [e[1] + e[2] + e[3]]]),
    ([[0, 1, 0], [0, 0, 1], [1, 0, 0]], [e[3] - e[1], e[3] - e[2]], [[e[3] - e[2], e[3] - e[1]], [e[3] - e[1], e[2] - e[1]], [e[1] - e[2] - e[3]]]),
    ([[0, -1, 0], [1, 0, 0], [0, 0, 1]], [e[2] - e[1], e[2]], [[e[2], e[2] - e[1]], [e[1] + e[2], e[2]], [e[3]]]),
]:
    M = Matrix(M)
    assert M.T * M == eye(3) and M.det() == 1
    prod_ = lambda vs: eval("*".join(["eye(3)"] + [f"refl(vs[{i}])" for i in range(len(vs))]))
    assert prod_(bon) == M
    for L in leurres:
        assert prod_(L) != M
# Remarque de la correction : la décomposition n'est pas unique (v0 : s_{e3-e1} ∘ s_{e2-e1} convient aussi).
assert refl(e[3] - e[1]) * refl(e[2] - e[1]) == Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])

# =========================================================================================================
# al-l2-euclidien-spectral (carte)
# =========================================================================================================
S = Matrix([[2, 1], [1, 2]])
Pp = Matrix([[1, 1], [1, -1]]) / sqrt(2)
assert egal(Pp * diag(3, 1) * Pp.T, S) and egal(Pp.T * Pp, eye(2))
Sc = Matrix([[1, I], [I, -1]])
assert (Sc * Sc).applyfunc(simplify) == zeros(2, 2) and Sc != zeros(2, 2)
Md = Matrix([[1, 1], [0, 2]])
assert Md * Matrix([1, 1]) == 2 * Matrix([1, 1]) and ps((1, 0), (1, 1)) != 0
assert Matrix([[1, 1], [0, 1]]).inv() * Md * Matrix([[1, 1], [0, 1]]) == diag(1, 2)
lam = symbols("lam")
assert (J - lam * eye(2)).det() == lam**2 + 1
# Trou 2 : S = [[0,1],[1,0]] n'est pas P D tP avec P triangulaire supérieure.
p11, p12, p22, d1, d2 = symbols("p11 p12 p22 d1 d2")
Pt = Matrix([[p11, p12], [0, p22]])
eqs = list(Pt * diag(d1, d2) * Pt.T - Matrix([[0, 1], [1, 0]]))
assert solve(eqs, [p11, p12, p22, d1, d2], dict=True) == []

# al-l2-euclidien-spectral-ex1 (+ v1, v2).
for S, l1, v1, l2, v2 in [([[5, 2], [2, 2]], 6, (2, 1), 1, (1, -2)), ([[-3, 4], [4, 3]], 5, (1, 2), -5, (2, -1)),
                          ([[7, -2], [-2, 4]], 8, (2, -1), 3, (1, 2))]:
    S = Matrix(S)
    P = Matrix.hstack(Matrix(v1), Matrix(v2)) / sqrt(5)
    D = diag(l1, l2)
    assert egal(P * D * P.T, S) and egal(P.T * P, eye(2))
    for Dl, Pl in [(diag(l2, l1), P), (D, Matrix.hstack(P[:, 1], P[:, 0])), (D, P * sqrt(5))]:
        assert not egal(Pl * Dl * Pl.T, S)
# al-l2-euclidien-spectral-ex2 (+ v1, v2) : S = alpha I + beta J.
for al, be, simple, double, leurres in [(1, 1, 4, 1, [(1, 4), (3, 0)]), (-1, 1, 2, -1, [(-1, 2), (3, 0)]), (2, -1, -1, 2, [(2, -1), (3, 0)])]:
    S = al * eye(3) + be * Matrix.ones(3, 3)
    assert S.eigenvals() == {simple: 1, double: 2}
    for (s_l, d_l) in leurres:
        assert {s_l: 1, d_l: 2} != S.eigenvals()
Pq = Matrix([[1 / sqrt(3), 1 / sqrt(2), 1 / sqrt(6)], [1 / sqrt(3), -1 / sqrt(2), 1 / sqrt(6)], [1 / sqrt(3), 0, -2 / sqrt(6)]])
assert egal(Pq.T * Pq, eye(3)) and egal(Pq * diag(4, 1, 1) * Pq.T, eye(3) + Matrix.ones(3, 3))
# leurres « trois valeurs propres simples » : (4,2,0), (2,1,-3), (-1,1,3) n'ont pas le bon produit (déterminant).
assert (eye(3) + Matrix.ones(3, 3)).det() == 4 != 4 * 2 * 0
assert (-eye(3) + Matrix.ones(3, 3)).det() == 2 != 2 * 1 * (-3)
assert (2 * eye(3) - Matrix.ones(3, 3)).det() == -4 != -1 * 1 * 3
# al-l2-euclidien-spectral-ex3 (+ v1, v2) : demi-axes.
for a, b2, c, k, d1_, d2_, leurres in [(5, 2, 8, 9, R(3, 2), 1, [(2, 3), (R(9, 4), 1), (R(2, 3), 1)]),
                                       (6, 2, 9, 10, sqrt(2), 1, [(sqrt(5), sqrt(10)), (2, 1), (sqrt(2) / 2, 1)]),
                                       (3, 1, 3, 16, 2, 2 * sqrt(2), [(2, sqrt(2)), (4, 8), (R(1, 2), sqrt(2) / 4)])]:
    S = Matrix([[a, b2], [b2, c]])
    ev = S.eigenvals()
    axes = sorted([sqrt(R(k, l)) for l in ev], key=float)
    assert sorted([d1_, d2_], key=float) == axes
    for L in leurres:
        assert sorted(L, key=float) != axes
# al-l2-euclidien-spectral-ex4 : rotation d'angle 2pi/3 : B^3 = I, B != I, non symétrique ; (-I)^3 = -I.
Brot = Matrix([[-R(1, 2), -sqrt(3) / 2], [sqrt(3) / 2, -R(1, 2)]])
assert egal(Brot**3, eye(2)) and Brot != Brot.T and (-eye(2)) ** 3 == -eye(2)

# =========================================================================================================
# al-l2-euclidien-positifs (carte)
# =========================================================================================================
Tq = Matrix([[sqrt(3) + 1, sqrt(3) - 1], [sqrt(3) - 1, sqrt(3) + 1]]) / 2
assert egal(Tq**2, Matrix([[2, 1], [1, 2]])) and egal(Pp * diag(sqrt(3), 1) * Pp.T, Tq)
Aa = Matrix([[1, 1], [0, 1]])
assert Aa.T * Aa == Matrix([[1, 1], [1, 2]])
assert set((Aa.T * Aa).eigenvals()) == {R(3, 2) + sqrt(5) / 2, R(3, 2) - sqrt(5) / 2}
assert simplify(sqrt(R(3, 2) + sqrt(5) / 2) - (1 + sqrt(5)) / 2) == 0
# Contre-exemples : [[1,2],[2,1]] en (1,-1) donne -2 ; -I2 ; diag(1,-1)^2 = I ; E12 de norme 1.
assert (Matrix([[1, -1]]) * Matrix([[1, 2], [2, 1]]) * Matrix([1, -1]))[0] == -2
assert diag(1, -1) ** 2 == eye(2) and (-eye(2)).det() == 1
assert max((E12.T * E12).eigenvals()) == 1
# Trou 1, leurre trace : diag(2,-1) ; trou 3 : u = I/2 donne 1/2 > 1/4 ; u = I, x = 2e : 4 > 2 ; lambda + |x|^2 : u = 2I, x = e1 : 2 <= 3 vrai,
# mais x = 3e1 : 18 > 2 + 9.
assert diag(2, -1).trace() > 0 and R(1, 2) > R(1, 4) and 4 > 2 and 2 * 9 > 2 + 9
# Trou 2 : I = diag(1,-1)^2 (autoadjointe non positive) ; 4I n'a pas de racine orthogonale (v orthogonale => v^2 orthogonale).
assert (4 * eye(2)).T * (4 * eye(2)) != eye(2)
# Affirmations : AB non symétrique.
assert Matrix([[1, 0], [0, 0]]) * Matrix([[1, 1], [1, 1]]) == Matrix([[1, 1], [0, 0]])

# al-l2-euclidien-positifs-ex1 (+ v1, v2) : contre-exemples des leurres.
aa = symbols("a", real=True)
S0 = Matrix([[2, aa, 0], [aa, 2, aa], [0, aa, 2]])
assert set(S0.eigenvals()) == {2, 2 - sqrt(2) * aa, 2 + sqrt(2) * aa}
Sv = lambda M, v: min(float(x) for x in M.subs(aa, nsimplify(v)).eigenvals())
assert Sv(S0, 1.3) > 0 and Sv(S0, 1.5) < 0                                  # |a| < sqrt2 ; leurres |a| < 2 et |a| <= sqrt2 faux
assert abs(Sv(S0, sqrt(2))) < 1e-12                                   # |a| = sqrt2 : valeur propre nulle
S1 = Matrix([[1, aa], [aa, 4]])
assert Sv(S1, 1.9) > 0 and abs(Sv(S1, 2)) < 1e-12 and Sv(S1, 3) < 0 and Sv(S1, -1) > 0   # leurres |a| <= 2, |a| < 4, 0 < a < 2
S2 = Matrix([[aa, 1, 1], [1, aa, 1], [1, 1, aa]])
assert set(S2.eigenvals()) == {aa + 2, aa - 1}
assert Sv(S2, 1.1) > 0 and Sv(S2, 0.5) < 0 and abs(Sv(S2, 1)) < 1e-12        # leurres a > -2, a >= 1, a > 0 faux
assert Sv(S0, 0.9) > 0                                                       # leurre |a| < 1 : manque a = 1.2 (définie positive)
assert Sv(S0, 1.2) > 0
# al-l2-euclidien-positifs-ex2 (+ v1, v2).
for A, nA, rho, fro, lmax in [([[3, 0], [4, 5]], 3 * sqrt(5), 5, 5 * sqrt(2), 45), ([[5, 3], [0, 4]], 2 * sqrt(10), 5, 5 * sqrt(2), 40),
                              ([[5, 0], [12, 13]], 5 * sqrt(13), 13, 13 * sqrt(2), 325)]:
    A = Matrix(A)
    assert max((A.T * A).eigenvals()) == lmax and simplify(sqrt(lmax) - nA) == 0
    assert max(abs(v) for v in A.eigenvals()) == rho and simplify(sqrt((A.T * A).trace()) - fro) == 0
    assert len({simplify(nA), rho, simplify(fro), lmax}) == 4
# al-l2-euclidien-positifs-ex3 (+ v1, v2).
for a, b, T, Lsw, Ld in [(5, 4, [[2, 1], [1, 2]], [[1, 2], [2, 1]], (3, 1)), (10, 6, [[3, 1], [1, 3]], [[1, 3], [3, 1]], (4, 2)),
                         (13, 12, [[3, 2], [2, 3]], [[2, 3], [3, 2]], (5, 1))]:
    S = Matrix([[a, b], [b, a]])
    T, Lsw = Matrix(T), Matrix(Lsw)
    assert T**2 == S and min(T.eigenvals()) > 0
    assert Lsw**2 == S and min(Lsw.eigenvals()) < 0                         # racine carrée non positive
    assert S.applyfunc(sqrt) ** 2 != S and diag(*Ld) ** 2 != S
# al-l2-euclidien-positifs-ex4 (+ v1, v2) : décomposition polaire.
for A, S, O in [([[1, 2], [2, 1]], [[2, 1], [1, 2]], [[0, 1], [1, 0]]), ([[-1, -3], [3, 1]], [[3, 1], [1, 3]], [[0, -1], [1, 0]]),
                ([[-3, -2], [2, 3]], [[3, 2], [2, 3]], [[-1, 0], [0, 1]])]:
    A, S, O = Matrix(A), Matrix(S), Matrix(O)
    assert O * S == A and O.T * O == eye(2) and S**2 == A.T * A and min(S.eigenvals()) > 0
    for L in (-O, eye(2), O.T, O * diag(1, -1), diag(1, -1) * O):
        assert L == O or L * S != A

# =========================================================================================================
# al-l2-euclidien-produit-vectoriel (carte)
# =========================================================================================================
x, y = Matrix([1, 2, 3]), Matrix([4, 5, 6])
w = x.cross(y)
assert w == Matrix([-3, 6, -3]) and ps(w, x) == 0 == ps(w, y)
assert ps(w, w) == 54 == ps(x, x) * ps(y, y) - ps(x, y) ** 2 and sqrt(54) == 3 * sqrt(6)
# Non-associativité ; base directe non orthonormale (e1,e2,2e3) : det(e1,e2,e3) = 1/2.
assert e[1].cross(e[1]).cross(e[2]) == zeros(3, 1) and e[1].cross(e[1].cross(e[2])) == -e[2]
assert Matrix.hstack(e[1], e[2], 2 * e[3]).inv() * e[3] == Matrix([0, 0, R(1, 2)])
assert Matrix.hstack(e[1], e[2], 2 * e[3]).det() == 2
# Double produit vectoriel (identité symbolique) et invariance par rotation.
xs3 = Matrix(symbols("x1:4"))
ys3 = Matrix(symbols("y1:4"))
zs3 = Matrix(symbols("z1:4"))
assert (xs3.cross(ys3.cross(zs3)) - (ps(xs3, zs3) * ys3 - ps(xs3, ys3) * zs3)).expand() == zeros(3, 1)
assert (ps(xs3.cross(ys3), xs3.cross(ys3)) + ps(xs3, ys3) ** 2 - ps(xs3, xs3) * ps(ys3, ys3)).expand() == 0
# Trou 3 : x = 2e1, y = e2 : |x ^ y| = 2 et non 4.
assert nrm((2 * e[1]).cross(e[2])) == 2
# al-l2-euclidien-produit-vectoriel-ex1 (+ v1, v2).
for A, B, C, aire in [((1, 0, 0), (0, 2, 0), (0, 0, 3), R(7, 2)), ((1, 1, 0), (2, 3, 1), (0, 2, 2), 3 * sqrt(3) / 2),
                      ((2, 0, 1), (4, 1, 3), (1, 2, 3), sqrt(65) / 2)]:
    A, B, C = Matrix(A), Matrix(B), Matrix(C)
    N = nrm((B - A).cross(C - A))
    assert simplify(N / 2 - aire) == 0 and len({simplify(N / 2), N, simplify(N**2 / 2), simplify(N / 6)}) == 4
# al-l2-euclidien-produit-vectoriel-ex2 (+ v1, v2).
for A, u, M, d in [((1, 0, 0), (1, 1, 1), (0, 1, 3), 2 * sqrt(2)), ((0, 1, 2), (2, 1, -2), (2, 0, 2), 2),
                   ((1, 1, 1), (0, 3, 4), (2, 0, 1), sqrt(41) / 5)]:
    A, u, M = Matrix(A), Matrix(u), Matrix(M)
    am = M - A
    lam_ = ps(am, u) / ps(u, u)
    H = A + lam_ * u                                                         # projeté orthogonal
    assert simplify(nrm(M - H) - d) == 0 and simplify(nrm(am.cross(u)) / nrm(u) - d) == 0
    vals = [simplify(d), simplify(nrm(am.cross(u))), simplify(abs(ps(am, u)) / nrm(u)), simplify(nrm(am))]
    assert len(set(vals)) == 4
# al-l2-euclidien-produit-vectoriel-ex3 (+ v1, v2).
for a, b, x0 in [((1, 1, 0), (1, -1, 2), (1, -1, -1)), ((1, 0, 2), (2, 1, -1), (-R(2, 5), 1, R(1, 5))),
                 ((0, 1, 1), (2, 1, -1), (-1, 1, -1))]:
    a, b, x0 = Matrix(a), Matrix(b), Matrix(x0)
    assert ps(a, b) == 0 and a.cross(b) / ps(a, a) == x0
    lm = symbols("lm")
    assert ((x0 + lm * a).cross(a) - b).expand() == zeros(3, 1)
    assert (-x0).cross(a) != b and (x0 + b).cross(a) != b and (x0 + a).cross(a) == b   # leurres faux
# al-l2-euclidien-produit-vectoriel-ex4 : u(x) ^ u(y) = det(u) u(x ^ y) pour u orthogonale (tests).
s3 = diag(1, 1, -1)
assert (s3 * e[1]).cross(s3 * e[3]) == e[2] == -s3 * e[1].cross(e[3])
assert s3 * e[1].cross(e[3]) == -e[2] and e[1].cross(e[3]) == -e[2]          # leurres u(x^y) et x^y faux ici
assert (-e[1]).cross(-e[2]) == e[1].cross(e[2])                               # u = -id : leurre -x^y faux
Rz = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
for Q in (Rz, refl((1, 2, 2)), Sa):
    assert Q.T * Q == eye(3)
    assert egal((Q * xs3).cross(Q * ys3), Q.det() * Q * xs3.cross(ys3))

# =========================================================================================================
# QCM al-l2-euclidien-q1 à q8
# =========================================================================================================
assert egal(schmidt([(1, 1), (1, 0)])[0][1], Matrix([1, -1]) / sqrt(2))                    # q2
assert schmidt([(1, 1), (1, 0)])[1][1] == Matrix([R(1, 2), -R(1, 2)])
pD = ps((1, 2, 3), (1, 1, 1)) / 3 * Matrix([1, 1, 1])                                        # q3
assert nrm(Matrix([1, 2, 3]) - pD) == sqrt(2) and nrm(pD) == 2 * sqrt(3)
Oq = Matrix([[3, -4], [4, 3]]) / 5                                                           # q5
assert Oq.T * Oq == eye(2)
for L in (Matrix([[1, 1], [-1, 1]]), Matrix([[1, 0], [1, 1]]), Matrix([[1, sqrt(3)], [sqrt(3), 1]]) / 2):
    assert not egal(L.T * L, eye(2))
Aq = Matrix([[0, 2], [0, 0]])                                                                # q7
assert max((Aq.T * Aq).eigenvals()) == 4
assert Matrix([1, 2, 0]).cross(Matrix([0, 1, 3])) == Matrix([6, -3, 1])                      # q8
# =========================================================================================================
# Compléments de la relecture indépendante
# =========================================================================================================
# al-l2-euclidien-adjoint, trou 3 : leurres (Ker u)^perp et Im u faux pour u = id (Ker u* = {0}, mais (Ker u)^perp
# = Im u = E) ; leurre Ker u faux pour u = E12 (Ker u = Vect(e1), Ker u* = Vect(e2)).
assert eye(2).T.nullspace() == [] and eye(2).rank() == 2
assert E12.nullspace() == [Matrix([1, 0])] and E12.T.nullspace() == [Matrix([0, 1])]
# al-l2-euclidien-schmidt, trou 3 : leurre « Vect(f_k) = Vect(e_k) » : pour e1 = (1,0), e2 = (1,1), f1 et f2 seraient
# colinéaires à e1 et e2, non orthogonaux : aucune famille orthonormale ne convient.
assert ps((1, 0), (1, 1)) != 0
# al-l2-euclidien-spectral, trou 2 : leurre « P nilpotente » : P D tP serait de rang < 2, donc différente de I2 ;
# trou 3, leurre « D nilpotente » : I2 = P N tP avec P orthogonale donnerait N = I2, non nilpotente.
Pn = Matrix([[0, 1], [0, 0]])
d1_, d2_ = symbols("d1_ d2_")
assert (Pn * diag(d1_, d2_) * Pn.T).rank() <= 1
# al-l2-euclidien-produit-vectoriel, trou 2 : le leurre [x,z,y] définirait y ^ x = -(x ^ y) (ici (1,2,3), (4,5,6)).
for z_ in (e[1], e[2], e[3]):
    assert Matrix.hstack(x, z_, y).det() == ps(-w, z_)
# al-l2-euclidien-positifs-ex2 (+ v1, v2) : rayon spectral <= norme subordonnée <= norme de Frobenius, inégalités
# strictes ; leurres (plus grand module des valeurs propres, Frobenius, lambda_max(tAA)) tous distincts.
for A_, nor in [(Matrix([[3, 0], [4, 5]]), 3 * sqrt(5)), (Matrix([[5, 3], [0, 4]]), 2 * sqrt(10)),
                (Matrix([[5, 0], [12, 13]]), 5 * sqrt(13))]:
    lmax = max((A_.T * A_).eigenvals())
    rho = max(abs(v) for v in A_.eigenvals())
    frob = sqrt(sum(c**2 for c in A_))
    assert sqrt(lmax) == nor and rho < nor < frob and len({nor, rho, frob, lmax}) == 4
# al-l2-euclidien-positifs (carte), énoncé : sans symétrie, tXSX >= 0 ne suffit pas (rotation d'angle pi/2).
xs2 = Matrix(symbols("xa xb"))
assert (xs2.T * Matrix([[0, -1], [1, 0]]) * xs2)[0] == 0
# al-l2-euclidien-q1 : f(t) = t donne (1/2)^2 = 1/4 <= 1/3, inégalité stricte pour f non constante.
assert integrate(t, (t, 0, 1)) ** 2 == R(1, 4) and integrate(t**2, (t, 0, 1)) == R(1, 3)
# al-l2-euclidien-cauchy-schwarz-ex4, question 2 : pour A orthogonale (rotation), tr(tAA) = n mais tr A = n
# seulement pour A = I_n ; matrice symétrique de trace 2, diag(2, 0), a tr(tAA) = 4 != 2.
assert (J.T * J).trace() == 2 and J.trace() == 0 and (diag(2, 0).T * diag(2, 0)).trace() == 4
# al-l2-euclidien-schmidt-ex4, question 3 : pour n = 3, la borne 3^(3/2) ~ 5,196 est atteinte par aucune matrice
# à coefficients ±1 (le maximum de |det| vaut 4), ce qui n'est pas demandé ; on vérifie seulement la borne.
assert max(abs(Matrix(3, 3, list(c)).det()) for c in product((-1, 1), repeat=9)) == 4 <= 3 * sqrt(3)

print("al_l2_euclidien : toutes les vérifications passent.")
