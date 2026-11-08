"""Vérifications calculatoires du chapitre « al-l3-formes-quadratiques » (algèbre linéaire, L3) : cartes, exercices et
QCM de content/*/algebre-lineaire/l3-al-l3-formes-quadratiques.json et content/qcm/al-l3-formes-quadratiques.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Chaque réduction de Gauss est vérifiée en développant, chaque signature par les signes des valeurs propres (règle de
Descartes sur le polynôme caractéristique, exacte puisque ses racines sont réelles).
Lancer : python3 verification/al_l3_formes_quadratiques.py (doit finir sans erreur).
"""

from itertools import product

from sympy import I, Matrix, Poly, Rational, conjugate, diag, expand, eye, factor, simplify, sqrt, symbols, zeros

R = Rational
x, y, z, lam, a = symbols("x y z lambda a")
x1, x2, x3, x4, x5, x6 = symbols("x1:7")
X3 = Matrix([x, y, z])


def mat(q, v=(x, y, z)):
    """Matrice (symétrique) de la forme quadratique q dans la base canonique."""
    n = len(v)
    return Matrix(n, n, lambda i, j: R(1, 2) * q.diff(v[i]).diff(v[j]))


def forme(M, v=(x, y, z)):
    X = Matrix(v)
    return expand((X.T * M * X)[0])


def signature(M):
    """(p, s) par la règle de Descartes appliquée au polynôme caractéristique (racines réelles : M symétrique)."""
    c = Poly(M.charpoly(lam).as_expr(), lam).all_coeffs()
    def changements(cs):
        cs = [simplify(t) for t in cs]
        cs = [t for t in cs if t != 0]
        return sum(1 for u, w in zip(cs, cs[1:]) if (u > 0) != (w > 0))
    p = changements(c)
    n = len(c) - 1
    cm = [t * (-1) ** (n - k) for k, t in enumerate(c)]
    s = changements(cm)
    return p, s


def egal(A, B):
    return simplify(Matrix(A) - Matrix(B)) == zeros(*Matrix(A).shape)


def adj(M):
    return M.H


# =========================================================================================================
# al-l3-formes-quadratiques-bilineaires (carte)
# =========================================================================================================
# Exemple : q = x1^2 - x2^2 + 4x1x2 - 6x2x3, matrice, déterminant -9, et tPMP dans la base ((1,0,0),(-2,1,0),(0,0,1)).
q = x**2 - y**2 + 4*x*y - 6*y*z
M = mat(q)
assert M == Matrix([[1, 2, 0], [2, -1, -3], [0, -3, 0]]) and M.det() == -9
P = Matrix([[1, -2, 0], [0, 1, 0], [0, 0, 1]])
assert P.T * M * P == Matrix([[1, 0, 0], [0, -5, -3], [0, -3, 0]]) and P.det() == 1
y1, y2, y3 = symbols("y1:4")
assert expand(q.subs({x: y1 - 2*y2, y: y2, z: y3})) == expand(y1**2 - 5*y2**2 - 6*y2*y3)
# Forme polaire de l'exemple.
xp, yp, zp = symbols("xp yp zp")
b = x*xp - y*yp + 2*(x*yp + y*xp) - 3*(y*zp + z*yp)
assert expand(b.subs({xp: x, yp: y, zp: z})) == expand(q)
# Trou 1 : les leurres P^{-1}MP, P M tP et t(P^{-1}) M P^{-1} diffèrent de tPMP sur des exemples.
Pt, Mt = Matrix([[1, 1], [0, 1]]), diag(1, 0)
assert Pt.T * Mt * Pt == Matrix([[1, 1], [1, 1]])
assert Pt * Mt * Pt.T != Pt.T * Mt * Pt                      # leurre P M tP faux
assert Pt.inv() * Mt * Pt != Pt.T * Mt * Pt                 # leurre P^{-1} M P faux
P2 = 2 * eye(2)
assert (P2.inv()).T * eye(2) * P2.inv() != P2.T * eye(2) * P2  # leurre t(P^{-1}) M P^{-1} faux
# Trou 2 : q(t) = t^2 sur R ; b(1,1) = 1 ; leurres q(x+y)-q(x)-q(y) = 2, (q(x+y)-q(x-y))/2 = 2 ; (q(1)+q(1))/4 = 1/2 != b(1,0) = 0.
qq = lambda t: t**2
assert qq(2) - 1 - 1 == 2 and R(qq(2) - qq(0), 2) == 2 and R(qq(1) + qq(1), 4) == R(1, 2)
# Trou 3 : x1^2 - x2^2 non dégénérée, isotrope (1,1), trace nulle ; 2x1x2 non dégénérée de diagonale nulle.
assert diag(1, -1).det() != 0 and (1 - 1) == 0 and diag(1, -1).trace() == 0
assert Matrix([[0, 1], [1, 0]]).det() != 0
# Contre-exemple Z/2Z : b(x,y) = x1y2 + x2y1 symétrique, b(x,x) = 2x1x2 = 0 mod 2, non nulle.
assert all((2 * u * v) % 2 == 0 for u, v in product(range(2), repeat=2)) and (1 * 1 + 0 * 0) % 2 == 1
# Affirmation 5 : forme positive, cône isotrope = noyau (exemple de bilineaires-ex2 : (x+y+z)^2 + (y+2z)^2).

# =========================================================================================================
# al-l3-formes-quadratiques-bilineaires-ex1 (+ v1, v2) : matrices, formes polaires, leurres.
for qe, Mb, leurres, bpol in [
    (x**2 + 3*y**2 - z**2 + 4*x*y - 6*y*z, [[1, 2, 0], [2, 3, -3], [0, -3, -1]],
     [[[1, 4, 0], [4, 3, -6], [0, -6, -1]], [[1, 2, -3], [2, 3, 0], [-3, 0, -1]], [[1, 2, 0], [2, 3, 3], [0, 3, -1]]],
     x*xp + 3*y*yp - z*zp + 2*(x*yp + y*xp) - 3*(y*zp + z*yp)),
    (2*x**2 - y**2 + z**2 + 6*x*y + 2*x*z, [[2, 3, 1], [3, -1, 0], [1, 0, 1]],
     [[[2, 6, 2], [6, -1, 0], [2, 0, 1]], [[2, 3, 0], [3, -1, 1], [0, 1, 1]], [[2, 3, 1], [3, 1, 0], [1, 0, 1]]],
     2*x*xp - y*yp + z*zp + 3*(x*yp + y*xp) + (x*zp + z*xp)),
    (-x**2 + 4*y**2 + 2*z**2 - 2*x*y + 8*y*z, [[-1, -1, 0], [-1, 4, 4], [0, 4, 2]],
     [[[-1, -2, 0], [-2, 4, 8], [0, 8, 2]], [[-1, -1, 4], [-1, 4, 0], [4, 0, 2]], [[-1, 1, 0], [1, 4, 4], [0, 4, 2]]],
     -x*xp + 4*y*yp + 2*z*zp - (x*yp + y*xp) + 4*(y*zp + z*yp)),
]:
    assert mat(qe) == Matrix(Mb)
    assert all(forme(Matrix(L)) != expand(qe) for L in leurres)
    assert expand(bpol.subs({xp: x, yp: y, zp: z})) == expand(qe)
    assert expand(bpol - (Matrix([[x, y, z]]) * Matrix(Mb) * Matrix([xp, yp, zp]))[0]) == 0

# al-l3-formes-quadratiques-bilineaires-ex2 (+ v1, v2) : déterminant nul, noyau, rang 2, Gauss, leurres.
for qe, Mb, k, k_faux, plan, gauss in [
    (x**2 + 2*y**2 + 5*z**2 + 2*x*y + 2*x*z + 6*y*z, [[1, 1, 1], [1, 2, 3], [1, 3, 5]], (1, -2, 1), (1, 2, 1),
     (1, -2, 1), (x + y + z)**2 + (y + 2*z)**2),
    (x**2 + 5*y**2 + 5*z**2 + 4*x*y + 2*x*z + 8*y*z, [[1, 2, 1], [2, 5, 4], [1, 4, 5]], (3, -2, 1), (3, 2, 1),
     (3, -2, 1), (x + 2*y + z)**2 + (y + 2*z)**2),
    (x**2 + 2*y**2 + 5*z**2 - 2*x*y + 4*x*z - 6*y*z, [[1, -1, 2], [-1, 2, -3], [2, -3, 5]], (-1, 1, 1), (1, 1, 1),
     (-1, 1, 1), (x - y + 2*z)**2 + (y - z)**2),
]:
    Mb = Matrix(Mb)
    assert mat(qe) == Mb and Mb.det() == 0 and Mb.rank() == 2
    assert Mb * Matrix(k) == zeros(3, 1) and Mb * Matrix(k_faux) != zeros(3, 1)
    assert expand(gauss) == expand(qe)
    # leurre « plan » : un vecteur du plan d'équation k.X = 0 hors du noyau
    v = Matrix([plan[1], -plan[0], 0])
    assert (Matrix(plan).T * v)[0] == 0 and Mb * v != zeros(3, 1)
    assert signature(Mb) == (2, 0)

# al-l3-formes-quadratiques-bilineaires-ex3 (+ v1, v2) : tPMP, P^{-1}MP, déterminants, leurres.
for Mb, Pb, N, faux, expr in [
    ([[1, 2], [2, 1]], [[1, 1], [1, -1]], [[6, 0], [0, -2]], [[[3, 0], [0, -1]], [[6, 0], [0, 2]], [[6, 2], [2, -2]]], (6, 0, -2)),
    ([[3, 1], [1, 3]], [[1, 1], [1, -1]], [[8, 0], [0, 4]], [[[4, 0], [0, 2]], [[8, 0], [0, -4]], [[8, 2], [2, 4]]], (8, 0, 4)),
    ([[1, 1], [1, -1]], [[1, -2], [2, 1]], [[1, -7], [-7, -1]], [[[1, -14], [-14, -1]], [[R(1, 5), R(-7, 5)], [R(-7, 5), R(-1, 5)]], [[1, 7], [7, -1]]], (1, -14, -1)),
]:
    Mb, Pb, N = Matrix(Mb), Matrix(Pb), Matrix(N)
    assert Pb.T * Mb * Pb == N and all(Matrix(f) != N for f in faux)
    assert N.det() == Pb.det()**2 * Mb.det() and (N.det() > 0) == (Mb.det() > 0)
    al, be = symbols("al be")
    w = Pb * Matrix([al, be])
    assert expand(forme(Mb, (x1, x2)).subs({x1: w[0], x2: w[1]})) == expand(expr[0]*al**2 + expr[1]*al*be + expr[2]*be**2)
assert Matrix([[1, 1], [1, -1]]).inv() * Matrix([[1, 2], [2, 1]]) * Matrix([[1, 1], [1, -1]]) == diag(3, -1)
assert Matrix([[1, 1], [1, -1]]).inv() * Matrix([[3, 1], [1, 3]]) * Matrix([[1, 1], [1, -1]]) == diag(4, 2)
assert Matrix([[1, 1], [1, -1]]) * Matrix([[-2], [1]]) == Matrix([[-1], [-3]])                        # v2 : M v

# al-l3-formes-quadratiques-bilineaires-ex4 : q(P) = P(0)P(1) sur R_2[X], base (1, X, X^2).
c0, c1, c2 = symbols("c0 c1 c2")
Pp = lambda t: c0 + c1*t + c2*t**2
qP = expand(Pp(0) * Pp(1))
MP = mat(qP, (c0, c1, c2))
assert MP == Matrix([[1, R(1, 2), R(1, 2)], [R(1, 2), 0, 0], [R(1, 2), 0, 0]]) and MP.rank() == 2
assert MP * Matrix([0, 1, -1]) == zeros(3, 1)                       # noyau Vect(X - X^2)
assert MP * Matrix([0, 1, 0]) != zeros(3, 1)                        # leurre Vect(X)
assert expand(qP.subs({c0: 1, c1: 0, c2: 0})) == 1                  # q(1) = 1 : X + (1 - X) non isotrope
assert all(expand(qP.subs({c0: 0, c1: s, c2: t})) == 0 for s, t in [(1, 0), (0, 1), (2, -3)])     # F0 isotrope
assert all(expand(qP.subs({c0: -s - t, c1: s, c2: t})) == 0 for s, t in [(1, 0), (0, 1), (2, -3)])  # F1 isotrope

# =========================================================================================================
# al-l3-formes-quadratiques-gauss (carte)
# =========================================================================================================
# Exemple : x^2 + 4xy + 2xz + 3y^2 + 2yz = (x+2y+z)^2 - (y+z)^2, noyau Vect((1,-1,1)).
q = x**2 + 4*x*y + 2*x*z + 3*y**2 + 2*y*z
assert expand((x + 2*y + z)**2 - (y + z)**2) == q and mat(q) * Matrix([1, -1, 1]) == zeros(3, 1) and mat(q).rank() == 2
assert expand(R(1, 2)*(x + y)**2 - R(1, 2)*(x - y)**2) == 2*x*y
# Contre-exemple et affirmation : (x-y)^2+(y-z)^2+(z-x)^2 = 2(x-(y+z)/2)^2 + 3/2 (y-z)^2, rang 2.
q = (x - y)**2 + (y - z)**2 + (z - x)**2
assert expand(2*(x - (y + z)/2)**2 + R(3, 2)*(y - z)**2 - q) == 0 and mat(q).rank() == 2
assert mat(q) * Matrix([1, 1, 1]) == zeros(3, 1)
# Affirmation 3 : x^2 + 2xy + 2y^2 = (x+y)^2 + y^2 = 2(y + x/2)^2 + x^2/2.
assert expand((x + y)**2 + y**2 - (x**2 + 2*x*y + 2*y**2)) == 0
assert expand(2*(y + x/2)**2 + x**2/2 - (x**2 + 2*x*y + 2*y**2)) == 0
# Trou 1 : 0.x^2 (r = 1, rang 0) réfute « quelconques », « positifs ou nuls », « entiers ».
# Trou 2 : x^2 + y^2 + (x+y)^2 sur R^2 (r = 3, rang 2) réfute « génératrice » et « deux à deux non proportionnelles ».
assert mat(x**2 + y**2 + (x + y)**2, (x, y)).rank() == 2
# Trou 3 : x^2 sur R (r = 1, noyau de dim 0) ; -x^2 (aucune valeur propre > 0) ; x^2 sur R^2 (dim 2).
assert mat(-x**2, (x,)) == Matrix([[-1]]) and mat(x**2, (x, y)).rank() == 1

# al-l3-formes-quadratiques-gauss-ex1 (+ v1, v2) : réductions et leurres.
for qe, bonne, leurres, sig, temoin in [
    (x**2 + 2*y**2 + 3*z**2 + 2*x*y - 2*x*z, (x + y - z)**2 + (y + z)**2 + z**2,
     [(x + y - z)**2 + (y + z)**2 - z**2, (x + y - z)**2 + (y - z)**2 + z**2, (x + y - z)**2 + y**2 + 2*z**2], (3, 0), None),
    (x**2 + 5*y**2 + z**2 + 4*x*y + 2*x*z + 2*y*z, (x + 2*y + z)**2 + (y - z)**2 - z**2,
     [(x + 2*y + z)**2 + (y - z)**2 + z**2, (x + 2*y + z)**2 + (y + z)**2 - z**2, (x + 2*y + z)**2 - (y - z)**2 + z**2], (2, 1), ((-3, 1, 1), -1)),
    (x**2 + 3*y**2 + 2*z**2 - 2*x*y + 4*x*z, (x - y + 2*z)**2 + 2*(y + z)**2 - 4*z**2,
     [(x - y + 2*z)**2 + 2*(y + z)**2 - 2*z**2, (x - y + 2*z)**2 + 2*(y - z)**2 - 4*z**2, (x - y + 2*z)**2 + (y + z)**2 - 4*z**2], (2, 1), ((-3, -1, 1), -4)),
]:
    assert expand(bonne - qe) == 0 and all(expand(L - qe) != 0 for L in leurres)
    assert signature(mat(qe)) == sig
    if temoin:
        assert qe.subs(dict(zip((x, y, z), temoin[0]))) == temoin[1]
assert expand(-(2*y + z)**2 + 5*y**2 + z**2 + 2*y*z - (y**2 - 2*y*z)) == 0              # v1, étape intermédiaire
assert expand(-(2*z - y)**2 + 3*y**2 + 2*z**2 - (2*y**2 + 4*y*z - 2*z**2)) == 0         # v2, étape intermédiaire

# al-l3-formes-quadratiques-gauss-ex2 (+ v1, v2) : carrés de formes liées.
for qe, gauss, rel, k, k_faux, sig in [
    ((x - y)**2 + (y - z)**2 + (z - x)**2, 2*(x - (y + z)/2)**2 + R(3, 2)*(y - z)**2, (x - y) + (y - z) + (z - x), (1, 1, 1), (1, -1, 0), (2, 0)),
    ((x + y)**2 + (y + z)**2 + (x - z)**2, 2*(x + (y - z)/2)**2 + R(3, 2)*(y + z)**2, (x + y) - (y + z) - (x - z), (1, -1, 1), (1, 1, -1), (2, 0)),
    ((x + y)**2 + (y + z)**2 - (x + 2*y + z)**2, R(1, 2)*(x - z)**2 - 2*(y + (x + z)/2)**2, (x + y) + (y + z) - (x + 2*y + z), (1, -1, 1), (1, 1, 1), (1, 1)),
]:
    Mq = mat(expand(qe))
    assert expand(gauss - qe) == 0 and expand(rel) == 0 and Mq.rank() == 2
    assert Mq * Matrix(k) == zeros(3, 1) and Mq * Matrix(k_faux) != zeros(3, 1) and signature(Mq) == sig
assert expand((x + y)**2 + (y + z)**2 - (x + 2*y + z)**2 + 2*(x + y)*(y + z)) == 0               # v2 : = -2(x+y)(y+z)

# al-l3-formes-quadratiques-gauss-ex3 (+ v1, v2) : formes sans carré.
for qe, prod, gauss, Lm, sig in [
    (x*y + y*z + z*x, (x + z)*(y + z) - z**2, R(1, 4)*(x + y + 2*z)**2 - R(1, 4)*(x - y)**2 - z**2, [[1, 1, 2], [1, -1, 0], [0, 0, 1]], (1, 2)),
    (x*y + 2*y*z + 3*z*x, (x + 2*z)*(y + 3*z) - 6*z**2, R(1, 4)*(x + y + 5*z)**2 - R(1, 4)*(x - y - z)**2 - 6*z**2, [[1, 1, 5], [1, -1, -1], [0, 0, 1]], (1, 2)),
    (x*y - y*z + z*x, (x - z)*(y + z) + z**2, R(1, 4)*(x + y)**2 - R(1, 4)*(x - y - 2*z)**2 + z**2, [[1, 1, 0], [1, -1, -2], [0, 0, 1]], (2, 1)),
]:
    assert expand(prod - qe) == 0 and expand(gauss - qe) == 0 and Matrix(Lm).det() == -2
    assert signature(mat(qe)) == sig
assert mat(x*y + y*z + z*x).eigenvals() == {1: 1, R(-1, 2): 2}
assert mat(x*y - y*z + z*x).eigenvals() == {-1: 1, R(1, 2): 2}
Mv1 = mat(x*y + 2*y*z + 3*z*x)
assert Mv1.trace() == 0 and Mv1.det() == R(3, 2)

# al-l3-formes-quadratiques-gauss-ex4 (+ v1, v2) : base orthogonale = colonnes de L^{-1}, leurres.
for qe, gauss, Lm, alpha, base, leurres in [
    (x**2 + 2*x*y + 3*y**2 + 4*y*z + z**2, (x + y)**2 + 2*(y + z)**2 - z**2, [[1, 1, 0], [0, 1, 1], [0, 0, 1]], (1, 2, -1),
     [(1, 0, 0), (-1, 1, 0), (1, -1, 1)], [[(1, 1, 0), (0, 1, 1), (0, 0, 1)], [(1, 0, 0), (1, 1, 0), (1, 1, 1)], [(1, 0, 0), (-1, 1, 0), (0, -1, 1)]]),
    (x**2 - 4*x*y + 2*x*z + 3*y**2 - 10*y*z - 6*z**2, (x - 2*y + z)**2 - (y + 3*z)**2 + 2*z**2, [[1, -2, 1], [0, 1, 3], [0, 0, 1]], (1, -1, 2),
     [(1, 0, 0), (2, 1, 0), (-7, -3, 1)], [[(1, -2, 1), (0, 1, 3), (0, 0, 1)], [(1, 0, 0), (2, 1, 0), (-7, 3, 1)], [(1, 0, 0), (2, 1, 0), (1, -3, 1)]]),
    (2*x**2 + 4*x*y - 4*x*z + 3*y**2 - 8*y*z + 3*z**2, 2*(x + y - z)**2 + (y - 2*z)**2 - 3*z**2, [[1, 1, -1], [0, 1, -2], [0, 0, 1]], (2, 1, -3),
     [(1, 0, 0), (-1, 1, 0), (-1, 2, 1)], [[(1, 1, -1), (0, 1, -2), (0, 0, 1)], [(1, 0, 0), (-1, 1, 0), (1, 2, 1)], [(1, 0, 0), (1, 1, 0), (-1, 2, 1)]]),
]:
    Mq = mat(qe)
    assert expand(gauss - qe) == 0
    Pb = Matrix(Lm).inv()
    assert Pb == Matrix.hstack(*[Matrix(v) for v in base])
    assert Pb.T * Mq * Pb == diag(*alpha)
    for L in leurres:
        G = Matrix.hstack(*[Matrix(v) for v in L])
        assert not (G.T * Mq * G).is_diagonal()
assert expand(x**2 - 4*x*y + 2*x*z + 3*y**2 - 10*y*z - 6*z**2 - (x - 2*y + z)**2 - (-y**2 - 6*y*z - 7*z**2)) == 0   # v1
assert expand(2*x**2 + 4*x*y - 4*x*z + 3*y**2 - 8*y*z + 3*z**2 - 2*(x + y - z)**2 - (y**2 - 4*y*z + z**2)) == 0     # v2

# al-l3-formes-quadratiques-gauss-ex5 : S = tTDT, D = diag(Δk/Δk-1), leurres.
S = Matrix([[2, 1, 0], [1, 2, 1], [0, 1, 2]])
Ds = [S[:k, :k].det() for k in range(1, 4)]
assert Ds == [2, 3, 4]
D = diag(2, R(3, 2), R(4, 3))
T = Matrix([[1, R(1, 2), 0], [0, 1, R(2, 3)], [0, 0, 1]])
assert T.T * D * T == S
assert expand(2*(x + y/2)**2 + R(3, 2)*(y + 2*z/3)**2 + R(4, 3)*z**2 - forme(S)) == 0
assert S.eigenvals() == {2: 1, 2 + sqrt(2): 1, 2 - sqrt(2): 1}
# Unicité de (T, D) (D déterminée par les mineurs) : les trois leurres diffèrent de D.
assert all(L != D for L in [diag(2, 3, 4), diag(R(1, 2), R(2, 3), R(3, 4)), diag(2 + sqrt(2), 2, 2 - sqrt(2))])
# Question 2 (récurrence) testée sur une matrice aléatoire à mineurs non nuls.
S4 = Matrix([[1, 2, 0, 1], [2, 3, 1, 0], [0, 1, 2, 2], [1, 0, 2, 1]])
L, Dl = S4.LDLdecomposition(hermitian=False)          # S4 = L Dl tL, avec T = tL triangulaire supérieure unipotente
mins = [S4[:k, :k].det() for k in range(1, 5)]
assert all(m != 0 for m in mins) and (L.T).T * Dl * L.T == S4 and all(L[k, k] == 1 for k in range(4))
assert [Dl[k, k] for k in range(4)] == [mins[0]] + [mins[k] / mins[k - 1] for k in range(1, 4)]

# =========================================================================================================
# al-l3-formes-quadratiques-sylvester (carte)
# =========================================================================================================
# Exemple : x^2 + y^2 + z^2 + 4xy = (x+2y)^2 - 3y^2 + z^2, valeurs propres 3, -1, 1.
q = x**2 + y**2 + z**2 + 4*x*y
assert expand((x + 2*y)**2 - 3*y**2 + z**2 - q) == 0 and mat(q).eigenvals() == {3: 1, -1: 1, 1: 1} and signature(mat(q)) == (2, 1)
# Contre-exemple : 2xy dans la base ((1,1),(1,-1)) : diag(2,-2).
Pb = Matrix([[1, 1], [1, -1]])
assert Pb.T * Matrix([[0, 1], [1, 0]]) * Pb == diag(2, -2)
# Trou 1 et affirmation 1 : 4I = t(2I) I (2I) ; I et -I mêmes rang et déterminant, signatures distinctes.
assert (2*eye(2)).T * eye(2) * (2*eye(2)) == 4*eye(2) and eye(2).det() == (-eye(2)).det() and signature(eye(2)) != signature(-eye(2))
# Trou 2 : q = 0 positive sur R^n (p = 0) ; -x^2 non dégénérée sur R (p = 0) ; x^2 nulle seulement sur {0} (p = 1).
assert signature(Matrix([[-1]])) == (0, 1) and signature(Matrix([[1]])) == (1, 0)
# Piège « sur C » : z^2 dans la base (i) vaut -1.
assert I**2 == -1

# al-l3-formes-quadratiques-sylvester-ex1 (+ v1, v2).
for qe, sig, rg in [((x + 2*y)**2 + (y - z)**2 - 3*z**2, (2, 1), 3), ((x - y + z)**2 - (y + 2*z)**2 - z**2, (1, 2), 3),
                    (2*(x - z)**2 - (y + z)**2, (1, 1), 2)]:
    Mq = mat(expand(qe))
    assert signature(Mq) == sig and Mq.rank() == rg
assert mat(expand(2*(x - z)**2 - (y + z)**2)) * Matrix([1, -1, 1]) == zeros(3, 1)

# al-l3-formes-quadratiques-sylvester-ex2 (+ v1, v2) : q = x^2+y^2+z^2 + 2a(xy+yz+zx), a = 2, -1, -1/2.
J = Matrix.ones(3, 3)
for av, vps, sig in [(2, {5: 1, -1: 2}, (1, 2)), (-1, {-1: 1, 2: 2}, (2, 1)), (R(-1, 2), {0: 1, R(3, 2): 2}, (2, 0))]:
    Mq = mat(x**2 + y**2 + z**2 + 2*av*(x*y + y*z + z*x))
    assert Mq == (1 - av)*eye(3) + av*J and Mq.eigenvals() == vps and signature(Mq) == sig
q = x**2 + y**2 + z**2 + 4*(x*y + y*z + z*x)
assert q.subs({x: 1, y: 1, z: 1}) == 15 and q.subs({x: 1, y: -1, z: 0}) == -2
q = x**2 + y**2 + z**2 - 2*(x*y + y*z + z*x)
assert q.subs({x: 1, y: 1, z: 1}) == -3 and q.subs({x: 1, y: -1, z: 0}) == 4
assert expand(x**2 + y**2 + z**2 - (x*y + y*z + z*x) - R(1, 2)*((x - y)**2 + (y - z)**2 + (z - x)**2)) == 0

# al-l3-formes-quadratiques-sylvester-ex3 (+ v1, v2) : définie positive selon a.
Ma = Matrix([[1, a, 0], [a, 1, a], [0, a, 1]])
assert factor(Ma.charpoly(lam).as_expr()) == factor((lam - 1)*((lam - 1)**2 - 2*a**2))
assert Ma * Matrix([1, 0, -1]) == Matrix([1, 0, -1]) and expand(Ma.det() - (1 - 2*a**2)) == 0
Mb_ = (1 - a)*eye(3) + a*J
assert factor(Mb_.det()) == factor((1 + 2*a)*(1 - a)**2) and Mb_[:2, :2].det() == expand(1 - a**2)
Mc = Matrix([[2, a, 0], [a, 2, a], [0, a, 2]])
assert expand(Mc.det() - 2*(4 - 2*a**2)) == 0 and Mc * Matrix([1, 0, -1]) == 2*Matrix([1, 0, -1])
for M_, bonne, leurre_vals in [
    (Ma, lambda t: abs(t) < 1/2**0.5, [R(6, 10), R(7, 10), R(-7, 10)]),       # |a| < 1 et |a| <= √2/2 et |a| < 1/2 : réfutés
    (Mb_, lambda t: -0.5 < t < 1, [R(-6, 10), R(3, 4), R(1, 1)]),
    (Mc, lambda t: abs(t) < 2**0.5, [R(3, 2), R(13, 10), R(-13, 10)]),
]:
    for t in [R(k, 20) for k in range(-50, 51)]:
        Mt = M_.subs(a, t)
        assert (signature(Mt) == (3, 0)) == bonne(float(t))
# Leurres : |a|<1 contient 0,75 (non défini positive pour v0) ; |a|<1/2 exclut 0,6 (défini positive pour v0).
assert signature(Ma.subs(a, R(3, 4))) != (3, 0) and signature(Ma.subs(a, R(6, 10))) == (3, 0)
assert signature(Ma.subs(a, sqrt(2)/2)) == (2, 0) and signature(Ma.subs(a, 1)) == (2, 1)
assert signature(Mb_.subs(a, R(-3, 4))) != (3, 0) and signature(Mb_.subs(a, R(3, 4))) == (3, 0)
assert signature(Mb_.subs(a, 1)) == (1, 0) and signature(Mb_.subs(a, 2)) == (1, 2) and signature(Mb_.subs(a, -1)) == (2, 1)
assert signature(Mb_.subs(a, R(-1, 2))) == (2, 0)
assert signature(Mc.subs(a, R(3, 2))) == (2, 1)        # √2 < 3/2 < 2 : leurre |a| < 2 réfuté
assert signature(Mc.subs(a, R(13, 10))) == (3, 0)      # 1 < 1,3 < √2 : leurre |a| < 1 réfuté
assert signature(Mc.subs(a, sqrt(2))) == (2, 0)        # leurre |a| <= √2 réfuté

# al-l3-formes-quadratiques-sylvester-ex4 : q(M) = tr(M^2), signature (n(n+1)/2, n(n-1)/2), n = 2, 3.
for n, sig in [(2, (3, 1)), (3, (6, 3))]:
    vs = symbols(f"m0:{n*n}")
    Mm = Matrix(n, n, vs)
    G = mat(expand((Mm * Mm).trace()), vs)
    assert signature(G) == sig and G.det() != 0 and sig == (n*(n + 1)//2, n*(n - 1)//2)
    assert sig != (n*n, 0) and sig != (n, n*n - n) and sig != (n*(n - 1)//2, n*(n + 1)//2)

# al-l3-formes-quadratiques-sylvester-ex5 : x1x2 + x3x4 + x5x6 : signature (3,3), totalement isotrope Vect(e1,e3,e5).
v6 = (x1, x2, x3, x4, x5, x6)
q6 = x1*x2 + x3*x4 + x5*x6
assert signature(mat(q6, v6)) == (3, 3)
assert q6.subs({x2: 0, x4: 0, x6: 0}) == 0
# Construction de la question 3 : e + f isotropes et orthogonaux deux à deux (b = diag(1,-1)).
B2 = diag(1, -1, 1, -1)
w1, w2 = Matrix([1, 1, 0, 0]), Matrix([0, 0, 1, 1])
assert (w1.T * B2 * w1)[0] == 0 and (w2.T * B2 * w2)[0] == 0 and (w1.T * B2 * w2)[0] == 0

# =========================================================================================================
# al-l3-formes-quadratiques-simultanee (carte)
# =========================================================================================================
A = Matrix([[1, 1], [1, 2]])
Bm = Matrix([[1, 1], [1, 3]])
P = Matrix([[1, -1], [0, 1]])
assert P.T * A * P == eye(2) and P.T * Bm * P == diag(1, 2)
assert factor((Bm - lam*A).det()) == factor((1 - lam)*(2 - lam))
# Contre-exemple : diag(1,-1) et [[0,1],[1,0]] : A^{-1}B a pour valeurs propres ±i.
assert (diag(1, -1).inv() * Matrix([[0, 1], [1, 0]])).eigenvals() == {I: 1, -I: 1}
# Contre-exemple et trou 1 : A = 2I, B = I : tPP = I/2, coefficients 1/2.
assert (diag(2, 2).inv() * eye(2)).eigenvals() == {R(1, 2): 2}
# Trou 1, leurre « diagonale » : A = [[2,1],[1,2]], P = diag(p,q) inversible : coefficient non diagonal pq ≠ 0.
p_, q_ = symbols("p q", nonzero=True)
assert (diag(p_, q_) * Matrix([[2, 1], [1, 2]]) * diag(p_, q_))[0, 1] == p_*q_
# Trou 2 : A = 2I, B = I : racines de det(B - λI) = 1, det(A - λB) = 2, det(AB - λI) = 2, au lieu de 1/2.
A2, B2_ = 2*eye(2), eye(2)
for poly_, rac in [((B2_ - lam*A2).det(), R(1, 2)), ((B2_ - lam*eye(2)).det(), 1), ((A2 - lam*B2_).det(), 2), ((A2*B2_ - lam*eye(2)).det(), 2)]:
    assert poly_.subs(lam, rac) == 0
# Affirmation 3 : A = diag(1,2), B = [[0,1],[1,0]] : A^{-1}B non symétrique.
C_ = diag(1, 2).inv() * Matrix([[0, 1], [1, 0]])
assert C_ == Matrix([[0, 1], [R(1, 2), 0]]) and C_ != C_.T

# al-l3-formes-quadratiques-simultanee-ex1 (+ v1, v2).
for A_, B_, rac, AB, vpB, leurres in [
    (Matrix([[1, 1], [1, 2]]), Matrix([[3, 1], [1, -1]]), {-1, 4}, Matrix([[5, 3], [-2, -2]]), {1 - sqrt(5), 1 + sqrt(5)},
     [{1 - sqrt(5), 1 + sqrt(5)}, {-4, 1}, {-1, R(1, 4)}]),
    (Matrix([[1, 1], [1, 2]]), Matrix([[0, 2], [2, 1]]), {-4, 1}, Matrix([[-2, 3], [2, -1]]), {R(1, 2) - sqrt(17)/2, R(1, 2) + sqrt(17)/2},
     [{R(1, 2) - sqrt(17)/2, R(1, 2) + sqrt(17)/2}, {-1, 4}, {R(-1, 4), 1}]),
    (Matrix([[2, 1], [1, 1]]), Matrix([[3, 1], [1, 1]]), {1, 2}, Matrix([[2, 0], [-1, 1]]), {2 - sqrt(2), 2 + sqrt(2)},
     [{2 - sqrt(2), 2 + sqrt(2)}, {R(1, 2), 1}, {-2, -1}]),
]:
    assert A_.trace() > 0 and A_.det() > 0
    assert set(Poly((B_ - lam*A_).det(), lam).all_roots()) == rac
    assert A_.inv() * B_ == AB and set(B_.eigenvals()) == vpB
    assert all(L != rac for L in leurres)

# al-l3-formes-quadratiques-simultanee-ex2 (+ v1, v2) : P = P1 O, leurres.
for A_, B_, P1, B1, O, Pb, D_, leurres in [
    (Matrix([[1, 1], [1, 2]]), Matrix([[2, 3], [3, 6]]), Matrix([[1, -1], [0, 1]]), Matrix([[2, 1], [1, 2]]),
     Matrix([[1, 1], [1, -1]])/sqrt(2), Matrix([[0, 2], [1, -1]])/sqrt(2), diag(3, 1),
     [Matrix([[1, 1], [1, -1]])/sqrt(2), Matrix([[0, 2], [1, -1]]), Matrix([[2, 0], [-1, 1]])/sqrt(2)]),
    (Matrix([[1, 1], [1, 2]]), Matrix([[1, 3], [3, 3]]), Matrix([[1, -1], [0, 1]]), Matrix([[1, 2], [2, -2]]),
     Matrix([[2, 1], [1, -2]])/sqrt(5), Matrix([[1, 3], [1, -2]])/sqrt(5), diag(2, -3),
     [Matrix([[2, 1], [1, -2]])/sqrt(5), Matrix([[1, 3], [1, -2]]), Matrix([[3, 1], [-2, 1]])/sqrt(5)]),
    (Matrix([[1, -1], [-1, 2]]), Matrix([[3, -2], [-2, 4]]), Matrix([[1, 1], [0, 1]]), Matrix([[3, 1], [1, 3]]),
     Matrix([[1, 1], [1, -1]])/sqrt(2), Matrix([[2, 0], [1, -1]])/sqrt(2), diag(4, 2),
     [Matrix([[1, 1], [1, -1]])/sqrt(2), Matrix([[2, 0], [1, -1]]), Matrix([[0, 2], [-1, 1]])/sqrt(2)]),
]:
    assert P1.T * A_ * P1 == eye(2) and P1.T * B_ * P1 == B1
    assert egal(O.T * O, eye(2)) and egal(O.T * B1 * O, D_) and egal(P1 * O, Pb)
    assert egal(Pb.T * A_ * Pb, eye(2)) and egal(Pb.T * B_ * Pb, D_)
    for L in leurres:
        assert not (egal(L.T * A_ * L, eye(2)) and egal(L.T * B_ * L, D_))
assert Matrix([[2, 3], [3, 6]]) * Matrix([[1, -1], [0, 1]]) == Matrix([[2, 1], [3, 3]])
assert Matrix([[1, 3], [3, 3]]) * Matrix([[1, -1], [0, 1]]) == Matrix([[1, 2], [3, 0]])
assert Matrix([[3, -2], [-2, 4]]) * Matrix([[1, 1], [0, 1]]) == Matrix([[3, 1], [-2, 2]])

# al-l3-formes-quadratiques-simultanee-ex3 (+ v1, v2) : extremums sur l'ellipse, contrôle numérique.
import math
for A_, B_, mx, mn, f, ell in [
    (diag(1, 2), Matrix([[0, 2], [2, 0]]), sqrt(2), -sqrt(2), lambda u, v: 4*u*v, (1, 2)),
    (diag(2, 1), Matrix([[0, 1], [1, 0]]), sqrt(2)/2, -sqrt(2)/2, lambda u, v: 2*u*v, (2, 1)),
    (diag(1, 4), Matrix([[1, 2], [2, 0]]), (1 + sqrt(5))/2, (1 - sqrt(5))/2, lambda u, v: u*u + 4*u*v, (1, 4)),
]:
    assert set(Poly((B_ - lam*A_).det(), lam).all_roots()) == {simplify(mx), simplify(mn)}
    vals = [f(math.cos(t)/math.sqrt(ell[0]), math.sin(t)/math.sqrt(ell[1])) for t in [k*2*math.pi/20000 for k in range(20000)]]
    assert abs(max(vals) - float(mx)) < 1e-6 and abs(min(vals) - float(mn)) < 1e-6
# Leurres : 2 (valeur propre de B), 2√2, 1 ; 1, √2, 1/2 ; (1+√17)/2, 1+√5, √5/2.
assert set(Matrix([[0, 2], [2, 0]]).eigenvals()) == {2, -2} and set(Matrix([[1, 2], [2, 0]]).eigenvals()) == {R(1, 2) - sqrt(17)/2, R(1, 2) + sqrt(17)/2}

# al-l3-formes-quadratiques-simultanee-ex4 : det(A+B) > det A + det B (n >= 2), contrôles numériques exacts.
for A_, B_ in [(Matrix([[2, 1], [1, 1]]), Matrix([[1, 0], [0, 3]])), (eye(3), diag(1, 2, 3)),
               (Matrix([[2, 1, 0], [1, 2, 1], [0, 1, 2]]), Matrix([[1, 1, 0], [1, 2, 0], [0, 0, 1]]))]:
    assert (A_ + B_).det() > A_.det() + B_.det()
assert (2*eye(2)).det() != 2 * eye(2).det()          # leurre « égalité lorsque B = A » réfuté pour n = 2

# al-l3-formes-quadratiques-simultanee-ex5 : Minkowski, égalité si et seulement si B = cA.
A_ = Matrix([[2, 1], [1, 1]])
for B_, egalite in [(3*A_, True), (A_, True), (Matrix([[1, 0], [0, 3]]), False), (diag(1, 2), False)]:
    g = float((A_ + B_).det())**0.5
    d = float(A_.det())**0.5 + float(B_.det())**0.5
    assert g >= d - 1e-12 and ((abs(g - d) < 1e-12) == egalite)
# Leurre « AB = BA » : A = I et B = diag(1,2) commutent sans égalité.
assert abs(float(diag(2, 3).det())**0.5 - (1 + float(diag(1, 2).det())**0.5)) > 1e-3
# Question 2 : prod(1+λ)^{1/n} >= 1 + prod(λ)^{1/n}, contrôle sur une grille.
for l1 in [0.1, 0.5, 1, 2, 7]:
    for l2 in [0.1, 0.5, 1, 2, 7]:
        assert ((1 + l1)*(1 + l2))**0.5 >= 1 + (l1*l2)**0.5 - 1e-12
        assert (abs(((1 + l1)*(1 + l2))**0.5 - 1 - (l1*l2)**0.5) < 1e-12) == (l1 == l2)

# =========================================================================================================
# al-l3-formes-quadratiques-coniques (carte)
# =========================================================================================================
S_ = Matrix([[1, R(1, 2)], [R(1, 2), 1]])
assert S_.eigenvals() == {R(3, 2): 1, R(1, 2): 1} and S_ * Matrix([1, 1]) == R(3, 2)*Matrix([1, 1])
assert Matrix([[0, R(1, 2)], [R(1, 2), 0]]).eigenvals() == {R(1, 2): 1, R(-1, 2): 1}
assert Matrix([[1, 1], [1, 2]]).det() == 1 and Matrix([[1, 2], [2, 2]]).det() == -2          # piège 2b
assert (x**2 + 2*x*y + 2*y**2).subs({x: 1, y: 0}) == 1
assert factor(x**2 + 2*x*y + y**2 + x + y) == (x + y)*(x + y + 1)                            # affirmation 1
assert factor(x**2 + x) == x*(x + 1)

# al-l3-formes-quadratiques-coniques-ex1 (+ v1, v2).
def centre(S, d, e):
    return S.inv() * Matrix([-R(d, 2), -R(e, 2)])
# ex : x^2 + 4xy + y^2 - 2x
S0 = Matrix([[1, 2], [2, 1]]); Om = centre(S0, -2, 0)
assert S0.det() == -3 and Om == Matrix([R(-1, 3), R(2, 3)])
F0 = (x**2 + 4*x*y + y**2 - 2*x)
assert F0.subs({x: Om[0], y: Om[1]}) == R(1, 3)
Xs, Ys = symbols("Xs Ys")
assert expand(F0.subs({x: Om[0] + Xs, y: Om[1] + Ys}) - (Xs**2 + 4*Xs*Ys + Ys**2 + R(1, 3))) == 0
assert S0.eigenvals() == {3: 1, -1: 1}
# v1 : x^2 + 2xy + 3y^2 + 2x
S1 = Matrix([[1, 1], [1, 3]]); Om = centre(S1, 2, 0)
F1 = x**2 + 2*x*y + 3*y**2 + 2*x
assert S1.det() == 2 and Om == Matrix([R(-3, 2), R(1, 2)]) and F1.subs({x: Om[0], y: Om[1]}) == R(-3, 2)
assert F1.subs({x: 0, y: 0}) == 0 and signature(S1) == (2, 0)
# v2 : x^2 - 2xy + y^2 + 2x : parabole (partie linéaire non orthogonale au noyau (1,1)).
S2 = Matrix([[1, -1], [-1, 1]])
assert S2.det() == 0 and S2 * Matrix([1, 1]) == zeros(2, 1) and Matrix([[2, 0]]) * Matrix([1, 1]) != zeros(1, 1)
U_, V_ = symbols("U V")
assert expand((x**2 - 2*x*y + y**2 + 2*x).subs({x: (U_ + V_)/sqrt(2), y: (V_ - U_)/sqrt(2)}) - (2*U_**2 + sqrt(2)*U_ + sqrt(2)*V_)) == 0

# al-l3-formes-quadratiques-coniques-ex2 (+ v1, v2) : demi-axes.
for S, k, vecs, demi in [
    (Matrix([[5, -3], [-3, 5]]), 8, [((1, 1), 2), ((1, -1), 8)], {2, 1}),
    (Matrix([[3, 1], [1, 3]]), 8, [((1, 1), 4), ((1, -1), 2)], {sqrt(2), 2}),
    (Matrix([[7, 2], [2, 4]]), 24, [((2, 1), 8), ((1, -2), 3)], {sqrt(3), 2*sqrt(2)}),
]:
    for v, l in vecs:
        assert S * Matrix(v) == l * Matrix(v)
    assert {sqrt(R(k, l)) for _, l in vecs} == demi
# Leurres (ensembles distincts de la bonne réponse).
assert {4, 1} != {2, 1} and {R(1, 2), 1} != {2, 1} and {sqrt(2), 2*sqrt(2)} != {2, 1}
assert {4, 2} != {sqrt(2), 2} and {sqrt(2)/2, R(1, 2)} != {sqrt(2), 2} and {sqrt(2), 1} != {sqrt(2), 2}
assert {8, 3} != {2*sqrt(2), sqrt(3)} and {sqrt(3)/3, sqrt(2)/4} != {2*sqrt(2), sqrt(3)} and {sqrt(3), sqrt(2)} != {2*sqrt(2), sqrt(3)}
# Contrôles des corrections.
assert (5*x**2 - 6*x*y + 5*y**2).subs({x: sqrt(2), y: sqrt(2)}) == 8
assert (3*x**2 + 2*x*y + 3*y**2).subs({x: sqrt(2), y: -sqrt(2)}) == 8
assert simplify((7*x**2 + 4*x*y + 4*y**2).subs({x: 2*sqrt(R(3, 5)), y: sqrt(R(3, 5))})) == 24

# al-l3-formes-quadratiques-coniques-ex3 (+ v1, v2) : paraboles, paramètre p.
Xp, Yp = symbols("Xp Yp")
# ex : X = (x+y)/√2, Y = (x-y)/√2.
e0 = x**2 + 2*x*y + y**2 - 2*sqrt(2)*x + 2*sqrt(2)*y
assert expand(e0.subs({x: (Xp + Yp)/sqrt(2), y: (Xp - Yp)/sqrt(2)}) - (2*Xp**2 - 4*Yp)) == 0          # X^2 = 2Y, p = 1
assert e0.subs({x: 2*sqrt(2), y: 0}) == 0
# v1 : X = (x-y)/√2, Y = -(x+y)/√2 : x = (X - Y)/√2, y = -(X + Y)/√2.
e1 = x**2 - 2*x*y + y**2 + 4*sqrt(2)*(x + y)
assert expand(e1.subs({x: (Xp - Yp)/sqrt(2), y: -(Xp + Yp)/sqrt(2)}) - (2*Xp**2 - 8*Yp)) == 0         # X^2 = 4Y, p = 2
assert expand(e1.subs({x: 0, y: -4*sqrt(2)})) == 0
# v2 : U = (x+y)/√2, V = (x-y)/√2, puis X = U - √2/2, Y = V : X^2 = 2√2 Y, p = √2, sommet (1/2, 1/2).
e2 = x**2 + 2*x*y + y**2 - 6*x + 2*y + 1
Uq = Xp + sqrt(2)/2
assert expand(e2.subs({x: (Uq + Yp)/sqrt(2), y: (Uq - Yp)/sqrt(2)}) - (2*Xp**2 - 4*sqrt(2)*Yp)) == 0
assert e2.subs({x: R(1, 2), y: R(1, 2)}) == 0
# p est intrinsèque : les leurres 2, 1/2, √2 (ex), 4, 1, 2√2 (v1), 2√2, 2, √2/2 (v2) diffèrent de la bonne valeur.

# al-l3-formes-quadratiques-coniques-ex4 (+ v1, v2) : xy + yz + zx = k.
Sq = mat(x*y + y*z + z*x)
assert Sq == (J - eye(3))/2 and Sq.eigenvals() == {1: 1, R(-1, 2): 2} and signature(Sq) == (1, 2)
O3 = Matrix.hstack(Matrix([1, 1, 1])/sqrt(3), Matrix([1, -1, 0])/sqrt(2), Matrix([1, 1, -2])/sqrt(6))
assert egal(O3.T * O3, eye(3)) and egal(O3.T * Sq * O3, diag(1, R(-1, 2), R(-1, 2)))
# Sections X = c : Y^2 + Z^2 = 2(c^2 - k) : k = 1 (vide pour |c| < 1 : deux nappes), k = -1 (toujours : une nappe), k = 0 (cône).

# al-l3-formes-quadratiques-coniques-ex5 : x^2 + 2mxy + y^2 + 2x + 2y = 0.
m = symbols("m")
Sm = Matrix([[1, m], [m, 1]])
Om = -Matrix([1, 1])/(1 + m)
assert simplify(Sm * Om + Matrix([1, 1])) == zeros(2, 1)
Fm = x**2 + 2*m*x*y + y**2 + 2*x + 2*y
assert simplify(Fm.subs({x: Om[0], y: Om[1]}) + 2/(1 + m)) == 0
assert simplify(expand(Fm.subs({x: Om[0] + Xs, y: Om[1] + Ys})) - (Xs**2 + 2*m*Xs*Ys + Ys**2 - 2/(1 + m))) == 0
assert factor(Fm.subs(m, 1)) == (x + y)*(x + y + 2)
assert expand(Fm.subs(m, -1).subs({x: (U_ + V_)/sqrt(2), y: (V_ - U_)/sqrt(2)}) - (2*U_**2 + 2*sqrt(2)*V_)) == 0
assert expand(Fm.subs(m, 0) - ((x + 1)**2 + (y + 1)**2 - 2)) == 0
for mv in [R(-1, 2), 0, R(1, 2)]:
    assert signature(Sm.subs(m, mv)) == (2, 0) and 2/(1 + mv) > 0                       # ellipses
for mv in [2, -2, 5, -5]:
    assert signature(Sm.subs(m, mv)) == (1, 1) and 2/(1 + mv) != 0                      # hyperboles

# =========================================================================================================
# al-l3-formes-quadratiques-hermitiens (carte)
# =========================================================================================================
def ph(u, v):
    return simplify((Matrix(u).H * Matrix(v))[0])

Uc = Matrix([[2, I*sqrt(3)], [-I*sqrt(3), 2]])
assert simplify(Uc.T * Uc) == eye(2) and simplify(Uc.H * Uc) != eye(2)
assert set(Uc.eigenvals()) == {2 + sqrt(3), 2 - sqrt(3)}
th = symbols("theta", real=True)
from sympy import cos, sin, exp
Rt = Matrix([[cos(th), -sin(th)], [sin(th), cos(th)]])
assert (Rt * Matrix([1, -I]) - exp(I*th)*Matrix([1, -I])).applyfunc(lambda t: simplify(t.rewrite(cos))) == zeros(2, 1)
assert (Rt * Matrix([1, I]) - exp(-I*th)*Matrix([1, I])).applyfunc(lambda t: simplify(t.rewrite(cos))) == zeros(2, 1)
assert ph((1, -I), (1, I)) == 0
assert ph((I,), (1,)) == -I                                              # affirmation 3 : <i.1, 1> = -i
# Base non orthonormale : u(e1) = 2e2, u(2e2) = e1 ; matrice [[0,1],[1,0]] dans (e1, 2e2).
Pb = Matrix([[1, 0], [0, 2]]); Mu = Pb * Matrix([[0, 1], [1, 0]]) * Pb.inv()
assert Mu * Matrix([1, 0]) == Matrix([0, 2])
# Polarisation hermitienne (cours) : <x,y> = 1/4 Σ i^{-k} ||x + i^k y||^2, testée numériquement.
xv, yv = Matrix([1 + 2*I, -I]), Matrix([3, 2 - I])
assert simplify(sum(I**(-k) * ph(xv + I**k*yv, xv + I**k*yv) for k in range(4))/4 - ph(xv, yv)) == 0
# Trou 2 : e^{i} n'est pas une racine de l'unité (π irrationnel) ; iI est unitaire de valeur propre i.
assert simplify((I*eye(2)).H * (I*eye(2))) == eye(2)

# al-l3-formes-quadratiques-hermitiens-ex1 (+ v1, v2).
for u, v, ps, nu, nv, leurres in [
    ((1, I), (2, 1 + I), 3 - I, 2, 6, [3 + I, 1 + I, 1 - I]),
    ((2, I), (1 - I, 3), 2 - 5*I, 5, 11, [2 + 5*I, 2 + I, 2 - I]),
    ((1, 1 + I), (I, 2), 2 - I, 3, 5, [2 + I, 2 + 3*I, 2 - 3*I]),
]:
    assert ph(u, v) == ps and ph(u, u) == nu and ph(v, v) == nv and ph(v, u) == conjugate(ps)
    assert simplify((Matrix(u).T * Matrix(v))[0]) == leurres[1]             # produit sans conjugaison
    assert ps not in leurres and abs(ps)**2 <= nu*nv

# al-l3-formes-quadratiques-hermitiens-ex2 (+ v1, v2).
for Um, vps, leurres in [
    (Matrix([[1 + I, 1 - I], [1 - I, 1 + I]])/2, {1, I}, [{1, -I}, {1 + I, 1 - I}, {(1 + I)/2, (1 - I)/2}]),
    (Matrix([[1, I], [I, 1]])/sqrt(2), {(1 + I)/sqrt(2), (1 - I)/sqrt(2)}, [{1 + I, 1 - I}, {1, I}, {(1 + I)/2, (1 - I)/2}]),
    (Matrix([[0, I], [I, 0]]), {I, -I}, [{1, -1}, {1, I}, {0, I}]),
]:
    assert simplify(Um.H * Um) == eye(2)
    ev = {simplify(t) for t in Um.eigenvals()}
    assert ev == {simplify(t) for t in vps} and all({simplify(t) for t in L} != ev for L in leurres)
    assert simplify(Um * Matrix([1, 1]) - Matrix([1, 1]) * simplify((Um * Matrix([1, 1]))[0])) == zeros(2, 1)

# al-l3-formes-quadratiques-hermitiens-ex3 (+ v1, v2) : Schmidt hermitien, leurres.
for v1, v2, f2, faux_, leurres in [
    ((1, I, 0), (0, 1, 1), Matrix([I, 1, 2])/sqrt(6), Matrix([-I, 3, 2])/sqrt(14), [Matrix([-I, 1, 2])/sqrt(6), Matrix([I/2, R(1, 2), 1])]),
    ((1, 0, I), (0, 1, 1), Matrix([I, 2, 1])/sqrt(6), Matrix([-I, 2, 3])/sqrt(14), [Matrix([-I, 2, 1])/sqrt(6), Matrix([I/2, 1, R(1, 2)])]),
    ((1, 1, 0), (I, 0, 1), Matrix([I, -I, 2])/sqrt(6), Matrix([3*I, I, 2])/sqrt(14), [Matrix([I, I, 2])/sqrt(6), Matrix([I/2, -I/2, 1])]),
]:
    V1, V2 = Matrix(v1), Matrix(v2)
    f1 = V1/sqrt(2)
    g2 = V2 - ph(f1, V2)*f1
    assert egal(g2/sqrt(ph(g2, g2)), f2) and ph(f2, f2) == 1 and ph(V1, f2) == 0
    gf = V2 - ph(V2, f1)*f1
    assert egal(gf/sqrt(ph(gf, gf)), faux_) and ph(V1, faux_) != 0
    for L in leurres:
        assert ph(V1, L) != 0 or ph(L, L) != 1

# al-l3-formes-quadratiques-hermitiens-ex4 : décomposition QR (contrôle sur un exemple), unicité.
Aq = Matrix([[1, I], [1, 2]])
cols = [Aq[:, 0], Aq[:, 1]]
q1 = cols[0]/sqrt(ph(cols[0], cols[0]))
g = cols[1] - ph(q1, cols[1])*q1
q2 = g/sqrt(ph(g, g))
Q = Matrix.hstack(q1, q2)
Rq = simplify(Q.H * Aq)
assert simplify(Q.H * Q) == eye(2) and Rq[1, 0] == 0 and Rq[0, 0].is_positive and simplify(Rq[1, 1]).is_positive
assert egal(Q * Rq, Aq)

# al-l3-formes-quadratiques-hermitiens-ex5 : Σ|λ|^2 <= tr(A*A), égalité si et seulement si A normale.
for Am, normale in [(Matrix([[1, 1], [0, 2]]), False), (I*eye(2), True), (2*eye(2), True), (Matrix([[1, 2], [-2, 1]]), True),
                    (Matrix([[0, 1], [0, 0]]), False)]:
    s_vp = sum(abs(t)**2 * mult for t, mult in Am.eigenvals().items())
    tr = simplify((Am.H * Am).trace())
    assert simplify(tr - s_vp) >= 0 and (simplify(tr - s_vp) == 0) == normale
    assert (simplify(Am.H * Am - Am * Am.H) == zeros(2, 2)) == normale
assert sum(abs(t)**2 for t in [1, 2]) == 5 and (Matrix([[1, 1], [0, 2]]).H * Matrix([[1, 1], [0, 2]])).trace() == 6

# =========================================================================================================
# al-l3-formes-quadratiques-endo-hermitiens (carte)
# =========================================================================================================
Hc = Matrix([[1, I], [-I, 1]])
assert Hc.H == Hc and Hc.eigenvals() == {0: 1, 2: 1}
assert Hc * Matrix([1, I]) == zeros(2, 1) and Hc * Matrix([1, -I]) == 2*Matrix([1, -I])
Uh = Matrix([[1, 1], [I, -I]])/sqrt(2)
assert simplify(Uh.H * Uh) == eye(2) and egal(Uh * diag(0, 2) * Uh.H, Hc)
Ns = Matrix([[1, I], [I, -1]])
assert Ns.T == Ns and Ns.H != Ns and Ns**2 == zeros(2, 2)
# Trou 1 : [[1,1],[0,0]] idempotente non hermitienne ; iI unitaire non hermitienne.
assert Matrix([[1, 1], [0, 0]])**2 == Matrix([[1, 1], [0, 0]]) and Matrix([[1, 1], [0, 0]]).H != Matrix([[1, 1], [0, 0]])
# Trou 2 : U triangulaire supérieure, D réelle diagonale : UDU* = [[0,1],[1,0]] impossible (coefficient (2,2) = |c|^2 d2).
aa, bb, cc, d1, d2 = symbols("aa bb cc d1 d2")
Ut = Matrix([[aa, bb], [0, cc]])
Pr = Ut * diag(d1, d2) * Ut.H
assert simplify(Pr[1, 1] - cc*conjugate(cc)*d2) == 0 and simplify(Pr[1, 0] - cc*d2*conjugate(bb)) == 0

# al-l3-formes-quadratiques-endo-hermitiens-ex1 (+ v1, v2).
for Hm, vps, Um, D_, leurres in [
    (Matrix([[3, 1 + I], [1 - I, 2]]), {1, 4}, Matrix([[(1 + I)/sqrt(3), (1 + I)/sqrt(6)], [1/sqrt(3), -2/sqrt(6)]]), diag(4, 1),
     [{2, 3}, {-1, 6}, {1 + I, 4 - I}]),
    (Matrix([[1, 1 - I], [1 + I, 2]]), {0, 3}, Matrix([[(1 - I)/sqrt(6), (1 - I)/sqrt(3)], [2/sqrt(6), -1/sqrt(3)]]), diag(3, 0),
     [{1, 2}, {-1, 4}, {1 - I, 2 + I}]),
    (Matrix([[3, 2*I], [-2*I, 0]]), {4, -1}, Matrix([[2*I, 1], [1, 2*I]])/sqrt(5), diag(4, -1),
     [{3, 0}, {5, -2}, {2*I, -2*I}]),
]:
    assert Hm.H == Hm and set(Hm.eigenvals()) == vps and all(L != vps for L in leurres)
    assert simplify(Um.H * Um) == eye(2) and egal(Um.H * Hm * Um, D_)

# al-l3-formes-quadratiques-endo-hermitiens-ex2 (+ v1, v2) : définie positive selon a.
for Hm, ok, vals in [
    (Matrix([[1, I*a], [-I*a, 1]]), lambda t: abs(t) < 1, [R(k, 4) for k in range(-8, 9)]),
    (Matrix([[2, a*(1 + I)], [a*(1 - I), 4]]), lambda t: abs(t) < 2, [R(k, 4) for k in range(-12, 13)]),
    (Matrix([[a, 1 + I], [1 - I, a]]), lambda t: t > 2**0.5, [R(k, 4) for k in range(-12, 13)]),
]:
    for t in vals:
        Ht = Hm.subs(a, t)
        assert Ht.H == Ht
        assert (all(simplify(e) > 0 for e in Ht.eigenvals())) == ok(float(t))
assert set(Matrix([[1, I*a], [-I*a, 1]]).eigenvals()) == {1 - a, 1 + a}
assert expand(Matrix([[2, a*(1 + I)], [a*(1 - I), 4]]).det()) == 8 - 2*a**2
assert set(Matrix([[a, 1 + I], [1 - I, a]]).eigenvals()) == {a - sqrt(2), a + sqrt(2)}
# Leurres réfutés : ex, a = 1 (|a| <= 1) et a = -1/2 (0 <= a < 1) ; v1, a = 3/2 (|a| < √2) et a = 2 (|a| <= 2) ;
# v2, a = 3/2 (a > 2), a = 5/4 (a > 1) et a = -2 (|a| > √2) : ces cas sont dans la grille ci-dessus.

# al-l3-formes-quadratiques-endo-hermitiens-ex3 : H hermitienne, H^3 = I => H = I ; jI unitaire avec (jI)^3 = I.
from sympy import pi
j = exp(2*I*pi/3)
assert simplify(j**3) == 1 and simplify(abs(j)) == 1 and Poly(lam**3 - 1, lam).count_roots() == 1   # une seule racine cubique réelle
assert (-eye(2))**3 == -eye(2)                                                                         # leurre -I_n

# al-l3-formes-quadratiques-endo-hermitiens-ex4 : transformée de Cayley.
for lv, mu in [(0, -1), (1, -I)]:
    assert simplify((lv - I)/(lv + I) - mu) == 0
Hm = Matrix([[R(1, 2), R(1, 2)], [R(1, 2), R(1, 2)]])                    # valeurs propres 0 et 1
assert Hm.eigenvals() == {0: 1, 1: 1}
Cm = (Hm - I*eye(2)) * (Hm + I*eye(2)).inv()
assert simplify(Cm.H * Cm) == eye(2) and {simplify(t) for t in Cm.eigenvals()} == {-1, -I}
assert egal(I*(eye(2) + Cm)*(eye(2) - Cm).inv(), Hm)

# al-l3-formes-quadratiques-endo-hermitiens-ex5 : Weyl, exemple H = diag(3,1), K = [[1,1],[1,1]].
Hm, Km = diag(3, 1), Matrix([[1, 1], [1, 1]])
assert set(Km.eigenvals()) == {0, 2} and set((Hm + Km).eigenvals()) == {3 + sqrt(2), 3 - sqrt(2)}
l_HK = sorted([float(t) for t in (Hm + Km).eigenvals()], reverse=True)
assert l_HK[0] >= 3 and l_HK[1] >= 1 and l_HK[0] <= 3 + 2 and l_HK[1] <= 1 + 2
# Leurre « λk(H+K) = λk(H) + λk(K) » : H = diag(1,0), K = diag(0,1).
assert (diag(1, 0) + diag(0, 1)).eigenvals() == {1: 2}

# =========================================================================================================
# QCM
# =========================================================================================================
# q1 : forme polaire de x1^2 - 6x1x2 + 2x2^2.
bq = x1*xp - 3*(x1*yp + x2*xp) + 2*x2*yp
assert expand(bq.subs({xp: x1, yp: x2})) == expand(x1**2 - 6*x1*x2 + 2*x2**2)
assert expand((x1*xp - 6*(x1*yp + x2*xp) + 2*x2*yp).subs({xp: x1, yp: x2})) == expand(x1**2 - 12*x1*x2 + 2*x2**2)
assert expand((x1*xp - 3*x1*yp + 2*x2*yp).subs({xp: x1, yp: x2})) == expand(x1**2 - 3*x1*x2 + 2*x2**2)
# q2 : (x-y)^2 + (y-z)^2 + (z-x)^2 nulle en (1,1,1).
assert ((x - y)**2 + (y - z)**2 + (z - x)**2).subs({x: 1, y: 1, z: 1}) == 0
# q3 : 2xy + z^2, signature (2,1).
assert signature(mat(2*x*y + z**2)) == (2, 1)
# q4 : voir le contre-exemple de la carte simultanee (valeurs propres ±i).
# q5 : x^2 + 6xy + y^2 = 1 : valeurs propres 4 et -2, hyperbole.
assert Matrix([[1, 3], [3, 1]]).eigenvals() == {4: 1, -2: 1}
# q6 : <(1,i),(1,-i)> = 0.
assert ph((1, I), (1, -I)) == 0 and simplify((Matrix([1, I]).T * Matrix([1, -I]))[0]) == 2
# q7 : [[0,-i],[i,0]] hermitienne, valeurs propres ±1.
Mq7 = Matrix([[0, -I], [I, 0]])
assert Mq7.H == Mq7 and Mq7.eigenvals() == {1: 1, -1: 1}
# q8 : H^2 hermitienne positive (exemple).
H8 = Matrix([[1, 2 + I], [2 - I, -3]])
assert H8.H == H8 and (H8**2).H == H8**2 and all(simplify(t) >= 0 for t in (H8**2).eigenvals())

print("al_l3_formes_quadratiques : toutes les vérifications passent.")
