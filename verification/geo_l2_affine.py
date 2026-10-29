"""Vérifications calculatoires du chapitre « geo-l2-affine » (géométrie affine, L2) :
cartes, exercices et QCM de content/*/geometrie/l2-geo-l2-affine.json et content/qcm/geo-l2-affine.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/geo_l2_affine.py (doit finir sans erreur).
"""

from fractions import Fraction as F
from itertools import combinations

from sympy import Matrix, Rational as R, eye, solve, symbols, zeros

m, s, t, a, b, c, d = symbols("m s t a b c d")


def det(*cols):
    """Déterminant de la matrice dont les colonnes sont les vecteurs donnés."""
    return Matrix([list(v) for v in cols]).T.det()


def cr(o, p, q):
    """det(OP, OQ) dans le plan."""
    return (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0])


# ===========================================================================
# geo-l2-affine-sous-espaces
# ===========================================================================
# Exemple : P : x+y+z=1, w=(1,0,0) hors de la direction, B=(0,2,3) -> point (-4,2,3).
assert 1 + 0 + 0 != 0 and -4 + 2 + 3 == 1
# Contre-exemple 2 : D1=(t,0,0), D2=(0,t,1) ; AB=(0,0,1) hors de Vect(e1,e2) : det non nul.
assert det((0, 0, 1), (1, 0, 0), (0, 1, 0)) != 0
# Trous : axes du plan, A=(1,0), B=(0,1), AB=(-1,1) hors des deux axes, et les axes se coupent.
assert (-1, 1)[0] != 0 and (-1, 1)[1] != 0
# Affirmation 2 : P1={(x,y,0,0)}, P2={(x,0,z,1)} : (0,0,0,1) hors de Vect(e1,e2,e3).
assert det((0, 0, 0, 1), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)) != 0

# geo-l2-affine-sous-espaces-ex1 (+ v1, v2) : point dans S, base de la direction, leurres.
def verif_sea(eqs, point, base, leurres):
    A = Matrix([e[:4] for e in eqs])
    rhs = Matrix([e[4] for e in eqs])
    noyau_dim = 4 - A.rank()
    assert A * Matrix(point) == rhs
    assert all(A * Matrix(v) == zeros(len(eqs), 1) for v in base)
    assert Matrix([list(v) for v in base]).rank() == len(base) == noyau_dim
    for p, bs in leurres:  # chaque leurre : point hors de S, ou un vecteur hors du noyau, ou mauvaise dimension
        faux = (A * Matrix(p) != rhs) or any(A * Matrix(v) != zeros(len(eqs), 1) for v in bs) or len(bs) != noyau_dim
        assert faux


verif_sea([(1, 1, 1, 1, 1), (1, -1, 0, 0, 2)], (2, 0, -1, 0), [(1, 1, -2, 0), (0, 0, 1, -1)],
          [((2, 0, -1, 0), [(1, 1, -2, 0), (1, -1, 0, 0)]), ((1, 1, -1, 0), [(1, 1, -2, 0), (0, 0, 1, -1)]),
           ((2, 0, -1, 0), [(1, 1, -2, 0)])])
verif_sea([(1, 2, -1, 1, 3), (0, 1, 1, -1, 1)], (1, 1, 0, 0), [(-3, 1, 0, 1), (0, 0, 1, 1)],
          [((1, 1, 0, 0), [(-3, 1, 0, 1), (1, 0, 0, 0)]), ((3, 0, 0, 0), [(-3, 1, 0, 1), (0, 0, 1, 1)]),
           ((1, 1, 0, 0), [(-3, 1, 0, 1)])])
verif_sea([(1, 1, 1, 1, 3), (1, 0, -1, 0, 1), (0, 1, 0, -1, 0)], (2, 0, 1, 0), [(1, -1, 1, -1)],
          [((2, 0, 1, 0), [(1, -1, 1, -1), (1, 0, 1, 0)]), ((1, 1, 1, 0), [(1, -1, 1, -1)]),
           ((2, 0, 1, 0), [(1, 1, 1, 1)])])


# geo-l2-affine-sous-espaces-ex2 (+ v1, v2) : valeur de m, point commun, leurres.
def verif_droites(A, u, B, v, m0, pt, leurre_pt):
    A, u, B, v = map(Matrix, (A, u, B, v))
    assert Matrix.hstack(u, v).rank() == 2  # non parallèles
    dd = Matrix.hstack(B - A, u, v).det()
    assert solve(dd, m) == [m0]
    Bm = B.subs(m, m0)
    sol = solve(list(A + t * u - Bm - s * v), [t, s], dict=True)
    assert len(sol) == 1 and list((A + t * u).subs(sol[0])) == list(pt)
    assert dd.subs(m, -m0) != 0  # leurre m = -m0
    # leurre : point qui n'est pas commun aux deux droites (pour m = m0)
    P = Matrix(leurre_pt)
    sur_A = solve(list(A + t * u - P), [t], dict=True)
    sur_B = solve(list(Bm + s * v - P), [s], dict=True)
    assert not (sur_A and sur_B)
    return dd


assert verif_droites((1, 0, 2), (1, 1, 0), (0, m, 1), (2, -1, 1), 2, (2, 1, 2), (3, 2, 2)).expand() == 2 - m
assert verif_droites((0, 1, -1), (1, 0, 2), (3, 2, m), (1, 1, 1), 4, (2, 1, 3), (3, 2, 4)).expand() == m - 4
assert verif_droites((1, 1, 0), (1, 0, -1), (0, m, 2), (1, 2, 0), 3, (-1, 1, 2), (1, 1, 0)).expand() == 3 - m

# geo-l2-affine-sous-espaces-ex3 : suites x_{n+1} = alpha x_n + beta ; bonne réponse et leurres (alpha=3, beta=2, lambda=0 ou 1).
al, be = F(3), F(2)
def solution(x0, N=6):
    xs = [x0]
    for _ in range(N):
        xs.append(al * xs[-1] + be)
    return xs
for lam in (F(0), F(1), F(-2)):
    assert [be / (1 - al) + lam * al**n for n in range(7)] == solution(be / (1 - al) + lam)
assert [be + 0 * al**n for n in range(7)] != solution(be)                       # leurre beta + lambda alpha^n
assert [be / (al - 1) + 0 * al**n for n in range(7)] != solution(be / (al - 1))  # leurre beta/(alpha-1)
assert [0 + be / (1 - al) * al**n for n in range(7)] != solution(be / (1 - al))  # leurre lambda + c alpha^n, lambda=0


# geo-l2-affine-sous-espaces-ex4 (+ v1, v2) : déterminant des directions, point commun, leurres hors de P ou de Q.
def verif_plans(A, u1, u2, B, v1, v2, detv, pt, leurres):
    A, u1, u2, B, v1, v2 = map(Matrix, (A, u1, u2, B, v1, v2))
    M = Matrix.hstack(u1, u2, v1, v2)
    assert M.det() == detv != 0
    x = M.solve(B - A)
    P = A + x[0] * u1 + x[1] * u2
    assert list(P) == list(pt) and P == B - x[2] * v1 - x[3] * v2
    for L in leurres:
        L = Matrix(L)
        dans_P = Matrix.hstack(u1, u2, L - A).rank() == 2
        dans_Q = Matrix.hstack(v1, v2, L - B).rank() == 2
        assert not (dans_P and dans_Q)


verif_plans((1, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 1), (0, 2, 1, 0), (0, 1, 0, 1), (1, 0, 0, 0), -1, (4, 3, 1, 1),
            [(3, 3, 1, 1), (4, 3, 1, 0)])
verif_plans((0, 1, 0, 0), (1, 0, 1, 0), (0, 1, 0, -1), (2, 0, 1, 1), (1, 1, 0, 0), (0, 0, 1, 1), -2,
            (R(3, 2), R(-1, 2), R(3, 2), R(3, 2)), [(R(-1, 2), R(3, 2), R(3, 2), R(3, 2)), (R(3, 2), R(-1, 2), R(3, 2), R(-3, 2))])
verif_plans((2, 0, 1, 0), (1, 1, 0, 0), (0, 0, 1, -1), (0, 1, 0, 3), (1, 0, 0, 1), (0, 1, 1, 0), 2,
            (R(1, 2), R(-3, 2), R(-5, 2), R(7, 2)), [(R(-1, 2), R(-3, 2), R(-5, 2), R(7, 2)), (R(1, 2), R(-3, 2), R(7, 2), R(-5, 2))])
# Déterminant de la correction principale (développement selon la dernière colonne).
assert Matrix([[1, 0, 1], [0, 1, 0], [0, 1, 1]]).det() == 1

# QCM q1, q2 : Grassmann (2+2-3 = 1 ; 2+2 = 4).
assert 2 + 2 - 3 == 1 and 2 + 2 == 4

# ===========================================================================
# geo-l2-affine-volume
# ===========================================================================
# Exemple : u=(3,1), v=(1,2) : 5, -5 ; dans B'=(2e1,e2) : 5/2.
assert det((3, 1), (1, 2)) == 5 and det((1, 2), (3, 1)) == -5
Bp = Matrix([[2, 0], [0, 1]])
assert (Bp.inv() * Matrix([[3, 1], [1, 2]])).det() == R(5, 2)
# Contre-exemples : (e1,-e2) donne -1, (2e1,e2) donne 1/2 ; affirmation 4 : matrice diag(2,1).
assert (Matrix([[1, 0], [0, -1]]).inv()).det() == -1 and (Bp.inv()).det() == R(1, 2) and Bp.det() == 2
# Affirmations 1 et 2 : (e2,e3,e1) de déterminant 1 ; (-e1,-e2) de déterminant 1 en dimension 2.
assert det((0, 1, 0), (0, 0, 1), (1, 0, 0)) == 1 and det((-1, 0), (0, -1)) == 1
# Trou 2, leurre det_B(B') : B'=(2e1,e2) donne 2 au lieu de 1/2 ; leurre |det| : B'=(e2,e1) donne -1.
assert det((0, 1), (1, 0)) == -1

# geo-l2-affine-volume-ex1 (+ v1, v2) et leurres.
assert det((1, 2, 0), (0, 1, 3), (2, 0, 1)) == 13
assert det((2, 1, 1), (1, 3, 0), (0, 1, 1)) == 6
assert det((1, 2, 1), (2, 1, 1), (1, 1, 2)) == -4
assert 11 != 13 and 4 != 6 and -6 != -4

# geo-l2-affine-volume-ex2 : déterminants des bases permutées ; cycle de longueur n : (-1)^(n-1).
assert [det(*B) for B in [((0, 1, 0), (1, 0, 0), (0, 0, 1)), ((0, 1, 0), (0, 0, 1), (1, 0, 0)),
                           ((-1, 0, 0), (0, -1, 0), (0, 0, -1)), ((1, 0, 0), (0, 1, 0), (0, 0, -1))]] == [-1, 1, -1, -1]
for n in range(2, 8):
    cycle = Matrix([[1 if i == (j + 1) % n else 0 for j in range(n)] for i in range(n)])  # colonnes e2,...,en,e1
    assert cycle.det() == (-1) ** (n - 1) and (-eye(n)).det() == (-1) ** n


# geo-l2-affine-volume-ex3 (+ v1, v2) : signes de part et d'autre de (AC), convexité, aires et leurres.
def verif_quad(A, B, C, D, aire, leurres):
    P = [A, B, C, D]
    tours = [cr(P[i], P[(i + 1) % 4], P[(i + 2) % 4]) for i in range(4)]
    assert all(x > 0 for x in tours) or all(x < 0 for x in tours)  # quadrilatère convexe
    d1, d2 = cr(A, C, B), cr(A, C, D)
    assert d1 * d2 < 0
    assert R(abs(d1) + abs(d2), 2) == aire
    assert aire not in leurres
    return d1, d2


assert verif_quad((0, 0), (5, 1), (4, 4), (1, 3), 12, [24, 8, 4]) == (-16, 8)
assert verif_quad((1, 0), (6, 2), (5, 5), (0, 3), 17, [34, R(17, 2), 15]) == (-17, 17)
assert verif_quad((0, 0), (6, 0), (5, 3), (1, 4), R(35, 2), [35, 9, R(1, 2)]) == (-18, 17)
assert R(abs(-16 + 8), 2) == 4 and R(abs(-18 + 17), 2) == R(1, 2)  # erreur typique citée

# geo-l2-affine-volume-ex4 (+ v1, v2) : det_B(B'), det_B(u,v,w), det_B'(u,v,w), contrôlé par les coordonnées dans B'.
for Bp_, uvw, dBp, duvw, rep, leurres in [
    ([(1, 1, 0), (0, 1, 1), (1, 0, 1)], [(1, 0, 1), (0, 2, 1), (1, 1, 0)], 2, -3, R(-3, 2), [-6, -3, R(3, 2)]),
    ([(1, 0, 0), (2, 1, 0), (1, 1, 3)], [(2, 1, 0), (1, 0, 1), (0, 1, 3)], 3, -5, R(-5, 3), [-15, -5, R(5, 3)]),
    ([(0, 1, 0), (1, 0, 0), (1, 1, 1)], [(1, 1, 1), (1, 2, 3), (0, 1, 4)], -1, 2, -2, [2, R(-1, 2), R(1, 2)]),
]:
    P = Matrix([list(v) for v in Bp_]).T
    U = Matrix([list(v) for v in uvw]).T
    assert P.det() == dBp and U.det() == duvw and (P.inv() * U).det() == rep == R(duvw, dBp)
    assert rep not in leurres


# geo-l2-affine-volume-ex5 (+ v1, v2) : test des segments, point commun, leurres.
def verif_seg(A, B, C, D, coupe, pt, leurres_pts):
    assert all(cr(*tr) != 0 for tr in combinations([A, B, C, D], 3))
    phiC, phiD, psiA, psiB = cr(A, B, C), cr(A, B, D), cr(C, D, A), cr(C, D, B)
    sol = solve([A[0] + s * (B[0] - A[0]) - C[0] - t * (D[0] - C[0]), A[1] + s * (B[1] - A[1]) - C[1] - t * (D[1] - C[1])], [s, t])
    vraiment = 0 <= sol[s] <= 1 and 0 <= sol[t] <= 1
    assert vraiment == coupe == (phiC * phiD < 0 and psiA * psiB < 0)
    X = (A[0] + sol[s] * (B[0] - A[0]), A[1] + sol[s] * (B[1] - A[1]))
    if coupe:
        assert X == pt and sol[t] == R(phiC, phiC - phiD)
    for L in leurres_pts:  # aucun leurre n'est dans [AB] ∩ [CD]
        assert not (coupe and L == X)
    return phiC, phiD, psiA, psiB, X


assert verif_seg((0, 0), (4, 2), (1, 3), (3, -1), True, (2, 1), [(1, R(1, 2)), (R(5, 2), R(5, 4))])[:4] == (10, -10, -10, 10)
r1 = verif_seg((0, 0), (4, 2), (1, 3), (2, 2), False, None, [(R(8, 3), R(4, 3)), (2, 1), (R(3, 2), R(5, 2))])
assert r1[:4] == (10, 4, -4, 2) and r1[4] == (R(8, 3), R(4, 3))  # intersection des droites, hors de [CD]
assert verif_seg((-1, 0), (3, 4), (0, 3), (2, -1), True, (R(2, 3), R(5, 3)), [(1, 1), (R(4, 3), R(7, 3))])[:4] == (8, -16, -10, 14)

# QCM q3 : (e1,-e2) orthonormée, déterminant -1.
assert det((1, 0), (0, -1)) == -1

# ===========================================================================
# geo-l2-affine-applications
# ===========================================================================
# Exemple : f(x,y) = (x+2y, 3x+y+4) : det -5, det(M-I) = -6, point fixe (-4/3, 0).
M = Matrix([[1, 2], [3, 1]])
assert M.det() == -5 and (M - eye(2)).det() == -6
assert M * Matrix([R(-4, 3), 0]) + Matrix([0, 4]) == Matrix([R(-4, 3), 0])
# Contre-exemple : f(x,y)=(x+1,2y) sans point fixe, det 2 ; affirmation 3 : h1∘h2 translation de vecteur O1O2 ;
# affirmation 4 : transvection de déterminant 1. Trou 3, leurres : homothétie de rapport 2 en dimension n (2^n ≠ 2^-n), symétrie (-1 ≠ 1).
assert Matrix([[1, 0], [0, 2]]).det() == 2 and Matrix([[1, 1], [0, 1]]).det() == 1
O1, O2, X0 = Matrix([0, 0]), Matrix([3, 5]), Matrix(symbols("x y"))
assert (O1 + 2 * (O2 + R(1, 2) * (X0 - O2) - O1) - X0).expand() == O2 - O1


# geo-l2-affine-applications-ex1 (+ v1, v2) : point fixe unique, leurres non fixes.
def verif_pf(M, C, pt, leurres):
    M, C = Matrix(M), Matrix(C)
    assert (M - eye(2)).det() != 0
    assert M * Matrix(pt) + C == Matrix(pt)
    for L in leurres:
        assert M * Matrix(L) + C != Matrix(L)
    return (M - eye(2)).det()


assert verif_pf([[2, 1], [1, 3]], [-1, 2], (4, -3), [(-4, 3), (3, -4)]) == 1
assert verif_pf([[3, -1], [1, 1]], [2, -1], (1, 4), [(4, 1), (-1, -4)]) == 1
assert verif_pf([[1, -2], [3, 2]], [4, -1], (R(-1, 3), 2), [(R(1, 3), 2), (2, R(-1, 3))]) == 6


# geo-l2-affine-applications-ex2 (+ v1, v2) : centre et rapport ; leurres (centre faux ou rapport faux).
def verif_hom(k, C, centre, faux_centres):
    C = Matrix(C)
    f = lambda P: k * Matrix(P) + C
    Om = Matrix(centre)
    assert f(Om) == Om
    P = Matrix([7, -2])
    assert f(P) == Om + k * (P - Om)
    for L in faux_centres:
        assert f(Matrix(L)) != Matrix(L)


verif_hom(-2, [3, 1], (1, R(1, 3)), [(R(1, 3), 1), (3, 1)])
verif_hom(3, [-4, 2], (2, -1), [(-2, 1), (-4, 2)])
verif_hom(R(1, 2), [1, -3], (2, -6), [(1, -3), (-2, 6)])


# geo-l2-affine-applications-ex3 (+ v1, v2) : composée de deux homothéties.
def hom(O, k):
    O = Matrix(O)
    return lambda P: O + k * (Matrix(P) - O)


def compo(O1, k1, O2, k2):
    h1, h2 = hom(O1, k1), hom(O2, k2)
    return lambda P: h1(h2(P))


g = compo((0, 0), 2, (3, 0), R(1, 4))
assert g((9, 0)) == Matrix([9, 0]) and g((1, 1)) == Matrix([9, 0]) + R(1, 2) * (Matrix([1, 1]) - Matrix([9, 0]))
assert g((1, 1)) - Matrix([1, 1]) != g((0, 0)) - Matrix([0, 0])  # ce n'est pas une translation
assert list(g((0, 0))) == [R(9, 2), 0]  # h1∘h2(x,y) = (x/2 + 9/2, y/2)
g = compo((1, 1), 3, (0, 2), -1)
assert g((R(-1, 2), R(5, 2))) == Matrix([R(-1, 2), R(5, 2)]) and list(g((0, 0))) == [-2, 10]
assert g((R(1, 2), R(-5, 2))) != Matrix([R(1, 2), R(-5, 2)])
g = compo((0, 0), 2, (1, 1), R(1, 2))
assert all(g(P) - Matrix(P) == Matrix([1, 1]) for P in [(0, 0), (3, -2), (R(1, 3), 5)])


# geo-l2-affine-applications-ex4 (+ v1, v2) : projection et symétrie ; leurres (symétrie de l'autre côté, orthogonale).
def verif_proj(n, dd, u, M0, p0, s0, s_faux, p_orth):
    n, u, M0 = Matrix(n), Matrix(u), Matrix(M0)
    assert n.dot(u) != 0
    lam = (dd - n.dot(M0)) / n.dot(u)
    p = M0 + lam * u
    assert n.dot(p) == dd and list(p) == list(p0) and list(M0 + 2 * lam * u) == list(s0)
    assert list(M0 + 2 * lam * u) != list(s_faux)
    po = M0 + (dd - n.dot(M0)) / n.dot(n) * n
    assert list(po) == list(p_orth[0]) and list(2 * po - M0) == list(p_orth[1]) and po != p
    return lam


assert verif_proj((1, 1, 1), 1, (1, 1, 0), (3, 1, 2), (R(1, 2), R(-3, 2), 2), (-2, -4, 2), (8, 6, 2),
                  ((R(4, 3), R(-2, 3), R(1, 3)), (R(-1, 3), R(-7, 3), R(-4, 3)))) == R(-5, 2)
assert verif_proj((1, -1, 2), 3, (0, 0, 1), (1, 2, 4), (1, 2, 2), (1, 2, 0), (1, 2, 6),
                  ((R(1, 3), R(8, 3), R(8, 3)), (R(-1, 3), R(10, 3), R(4, 3)))) == -2
assert verif_proj((2, 1, -1), 1, (1, 0, 1), (2, 1, 0), (-2, 1, -4), (-6, 1, -8), (10, 1, 8),
                  ((R(2, 3), R(1, 3), R(2, 3)), (R(-2, 3), R(-1, 3), R(4, 3)))) == -4
# Expressions générales de p et s (version principale et variantes).
x, y, z = symbols("x y z")
lam0 = (1 - x - y - z) / 2
assert [e.expand() for e in (x + lam0, y + lam0, x + 2 * lam0, y + 2 * lam0)] == [
    ((1 + x - y - z) / 2).expand(), ((1 - x + y - z) / 2).expand(), 1 - y - z, 1 - x - z]
lam1 = (3 - x + y - 2 * z) / 2
assert (z + lam1).expand() == ((3 - x + y) / 2).expand() and (z + 2 * lam1).expand() == 3 - x + y - z
lam2 = 1 - 2 * x - y + z
assert [(x + lam2).expand(), (z + lam2).expand(), (x + 2 * lam2).expand(), (z + 2 * lam2).expand()] == [
    1 - x - y + z, 1 - 2 * x - y + 2 * z, 2 - 3 * x - 2 * y + 2 * z, 2 - 4 * x - 2 * y + 3 * z]


# geo-l2-affine-applications-ex5 (+ v1, v2) : triangle des céviennes (Routh), points et rapport ; leurres.
def inter(P, Q, Rr, S):
    sol = solve(list(P + a * (Q - P) - Rr - b * (S - Rr)), [a, b])
    return P + sol[a] * (Q - P)


for tt, pts, rapport, leurres in [
    (R(1, 3), [(R(4, 7), R(2, 7)), (R(1, 7), R(4, 7)), (R(2, 7), R(1, 7))], R(1, 7), [R(1, 9), R(1, 6), R(1, 4)]),
    (R(1, 4), [(R(9, 13), R(3, 13)), (R(1, 13), R(9, 13)), (R(3, 13), R(1, 13))], R(4, 13), [R(1, 4), R(1, 13), R(3, 13)]),
    (R(2, 5), [(R(9, 19), R(6, 19)), (R(4, 19), R(9, 19)), (R(6, 19), R(4, 19))], R(1, 19), [R(1, 7), R(4, 25), R(1, 5)]),
]:
    A0, B0, C0 = Matrix([0, 0]), Matrix([1, 0]), Matrix([0, 1])
    Ap, Bp_, Cp = B0 + tt * (C0 - B0), C0 + tt * (A0 - C0), A0 + tt * (B0 - A0)
    P, Q, Rr = inter(A0, Ap, B0, Bp_), inter(B0, Bp_, C0, Cp), inter(C0, Cp, A0, Ap)
    assert [tuple(P), tuple(Q), tuple(Rr)] == pts
    assert abs(Matrix.hstack(Q - P, Rr - P).det()) / 2 / R(1, 2) == rapport == (1 - 2 * tt) ** 2 / (1 - tt + tt**2)
    assert rapport not in leurres
    # Invariance affine : même rapport pour un triangle quelconque (image par une bijection affine).
    Lm, Cv = Matrix([[2, 1], [-1, 3]]), Matrix([5, -2])
    f = lambda X: Lm * X + Cv
    A1, B1, C1 = f(A0), f(B0), f(C0)
    Ap1, Bp1, Cp1 = B1 + tt * (C1 - B1), C1 + tt * (A1 - C1), A1 + tt * (B1 - A1)
    P1, Q1, R1 = inter(A1, Ap1, B1, Bp1), inter(B1, Bp1, C1, Cp1), inter(C1, Cp1, A1, Ap1)
    assert abs(Matrix.hstack(Q1 - P1, R1 - P1).det()) / abs(Matrix.hstack(B1 - A1, C1 - A1).det()) == rapport
# Équations des droites (version principale) : (AA') x=2y, (BB') 2x+3y=2, (CC') 3x+y=1.
assert (R(2, 3), R(1, 3))[0] == 2 * R(1, 3) and 2 * 0 + 3 * R(2, 3) == 2 and 3 * R(1, 3) + 0 == 1

# QCM q4 : f(x,y) = (x+2, 3y) sans point fixe.
assert solve([x + 2 - x, 3 * y - y], [x, y], dict=True) == []

# ===========================================================================
# geo-l2-affine-barycentres
# ===========================================================================
# Affirmation 1 : barycentre de (0,2), (1,-1) = -1 ; affirmation 3 : 7/4 contre 3/2 ; contre-exemple x -> x^2.
assert F(2 * 0 - 1, 1) == -1 and F(0 + 1 + 6, 4) == F(7, 4) and F(0 * 1 + 2 * 3, 4) == F(3, 2)
assert ((0 + 2) / 2) ** 2 != (0**2 + 2**2) / 2


# geo-l2-affine-barycentres-ex1 (+ v1, v2) : barycentre et leurres.
def bary(pts, poids):
    L = sum(poids)
    return tuple(F(sum(w * p[i] for p, w in zip(pts, poids)), L) for i in range(2))


assert bary([(1, 2), (3, -1), (0, 4)], [2, -1, 3]) == (F(-1, 4), F(17, 4))
assert tuple(v * 4 / 6 for v in (F(-1, 4), F(17, 4))) == (F(-1, 6), F(17, 6))  # leurre : division par 6
assert (F(2 * 1 + 3, 4), F(2 * 2 - 1 + 12, 4)) == (F(5, 4), F(15, 4))           # leurre : signe du poids de B oublié
assert bary([(0, 1), (2, 3), (-1, 0)], [1, 3, -2]) == (4, 5)
assert (F(8, 6), F(10, 6)) == (F(4, 3), F(5, 3)) and (F(0 + 6 - 2, 2), F(1 + 9, 2)) == (2, 5)
assert bary([(2, 0), (-1, 1), (3, -2)], [-1, 2, 4]) == (F(8, 5), F(-6, 5))
assert (F(2 - 2 + 12, 5), F(-6, 5)) == (F(12, 5), F(-6, 5))


# geo-l2-affine-barycentres-ex2 (+ v1, v2) : coordonnées barycentriques et leurres.
def coords(A, B, C, M):
    sol = Matrix([[B[0] - A[0], C[0] - A[0]], [B[1] - A[1], C[1] - A[1]]]).solve(Matrix([M[0] - A[0], M[1] - A[1]]))
    return (1 - sol[0] - sol[1], sol[0], sol[1])


def point(A, B, C, w):
    return tuple(w[0] * A[i] + w[1] * B[i] + w[2] * C[i] for i in range(2))


for A, B, C, M, rep, dt, faux in [
    ((1, 1), (5, 2), (2, 6), (3, 3), (R(5, 19), R(8, 19), R(6, 19)), 19, [(R(8, 19), R(6, 19), R(5, 19)), (R(1, 3),) * 3]),
    ((0, 0), (4, 1), (1, 3), (4, 4), (R(-9, 11), R(8, 11), R(12, 11)), 11, [(R(-9, 11), R(12, 11), R(8, 11))]),
    ((2, 1), (6, 3), (0, 5), (3, 3), (R(3, 10), R(2, 5), R(3, 10)), 20, [(R(2, 5), R(3, 10), R(3, 10)), (R(1, 4), R(1, 2), R(1, 4))]),
]:
    assert cr(A, B, C) == dt and coords(A, B, C, M) == rep
    for w in faux:
        assert point(A, B, C, w) != tuple(M)
assert point((1, 1), (5, 2), (2, 6), (R(1, 3),) * 3) == (R(8, 3), 3)

# geo-l2-affine-barycentres-ex3 et q5 : AG = 3/4 AG_A ; milieux des arêtes opposées (tétraèdre quelconque).
T4 = [Matrix([0, 0, 0]), Matrix([3, 1, 0]), Matrix([1, 4, 2]), Matrix([-1, 2, 5])]
G = sum(T4, zeros(3, 1)) / 4
GA = (T4[1] + T4[2] + T4[3]) / 3
assert G - T4[0] == R(3, 4) * (GA - T4[0])
for (i, j), (k, l) in [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]:
    assert ((T4[i] + T4[j]) / 2 + (T4[k] + T4[l]) / 2) / 2 == G

# geo-l2-affine-barycentres-ex4 (+ v1, v2) : Ceva.
for r1_, r2_, poids, rep_rapport, faux_rapports in [
    (2, 3, (6, 1, 2), R(1, 6), [6]), (2, 5, (10, 1, 2), R(1, 10), [10]), (R(1, 2), 3, (3, 2, 1), R(2, 3), [R(3, 2), R(1, 6)]),
]:
    al_, be_, ga_ = map(R, poids)
    assert ga_ / be_ == r1_ and al_ / ga_ == r2_ and be_ / al_ == rep_rapport and r1_ * r2_ * rep_rapport == 1
    assert rep_rapport not in faux_rapports
    # Contrôle géométrique sur le triangle (0,0), (1,0), (0,1) : P est sur (AA') et (BB').
    A0, B0, C0 = Matrix([0, 0]), Matrix([1, 0]), Matrix([0, 1])
    Ap = B0 + R(r1_) / (1 + R(r1_)) * (C0 - B0)          # BA' = r1 A'C
    Bp_ = C0 + R(r2_) / (1 + R(r2_)) * (A0 - C0)         # CB' = r2 B'A
    P = inter(A0, Ap, B0, Bp_)
    S = al_ + be_ + ga_
    assert P == (al_ * A0 + be_ * B0 + ga_ * C0) / S
    Cp = inter(C0, P, A0, B0)
    k = (Cp - A0)[0] / (B0 - Cp)[0]                 # AC' = k C'B
    assert k == rep_rapport
assert (R(6, 9), R(1, 9), R(2, 9)) == (R(2, 3), R(1, 9), R(2, 9)) and (R(3, 6), R(2, 6), R(1, 6)) == (R(1, 2), R(1, 3), R(1, 6))

# QCM q6 : somme des poids 1 + 1 - 1 = 1.
assert 1 + 1 - 1 != 0

# ===========================================================================
# geo-l2-affine-convexite
# ===========================================================================
# Exemple : (1/2,1/4,1/4) positifs ; (-1/2,1,1/2) de somme 1 avec une coordonnée négative.
assert F(1, 2) + F(1, 4) + F(1, 4) == 1 and F(-1, 2) + 1 + F(1, 2) == 1
# Contre-exemples et affirmations : milieu (0,1) de (-1,1) et (1,1) hors de la parabole ; conv(A∩B) ≠ conv(A)∩conv(B).
assert 0**2 != 1
# Trou 1, leurres isobarycentres / milieux : Q est stable par moyennes, pas convexe (sqrt(2)/2 ∈ [0,1] irrationnel).
# geo-l2-affine-convexite-ex1 : contre-exemples de convexité.
assert 1 * 1 >= 1 and (-1) * (-1) >= 1 and not (0 * 0 >= 1)
assert not (1 <= 0**2 + 0**2 <= 4) and 1 <= 1 <= 4
assert -8 >= (-2) ** 3 and 0 >= 0**3 and not (-4 >= (-1) ** 3)
# Convexité de y >= x^2 : (1-t)x1^2 + t x2^2 - ((1-t)x1 + t x2)^2 = t(1-t)(x1-x2)^2 >= 0.
x1, x2 = symbols("x1 x2")
assert ((1 - t) * x1**2 + t * x2**2 - ((1 - t) * x1 + t * x2) ** 2 - t * (1 - t) * (x1 - x2) ** 2).expand() == 0


# geo-l2-affine-convexite-ex3 (+ v1, v2) : enveloppe convexe (chaîne monotone exacte), aire, leurres.
def enveloppe(P):
    P = sorted(set(P))
    def chaine(pts):
        h = []
        for p in pts:
            while len(h) >= 2 and cr(h[-2], h[-1], p) <= 0:
                h.pop()
            h.append(p)
        return h
    bas, haut = chaine(P), chaine(list(reversed(P)))
    return bas[:-1] + haut[:-1]


def aire(H):
    return R(sum(H[i][0] * H[(i + 1) % len(H)][1] - H[(i + 1) % len(H)][0] * H[i][1] for i in range(len(H))), 2)


for pts, ordre, ar, tours, interieurs in [
    ([(0, 0), (6, 0), (5, 4), (1, 5), (2, 1), (3, 6)], [(0, 0), (6, 0), (5, 4), (3, 6), (1, 5)], R(51, 2), [24, 6, 6, 9, 30],
     [((2, 1), ((0, 0), (6, 0), (5, 4)), (R(5, 8), R(1, 8), R(1, 4)))]),
    ([(0, 0), (4, 0), (6, 3), (2, 5), (3, 2), (0, 3)], [(0, 0), (4, 0), (6, 3), (2, 5), (0, 3)], 21, [12, 16, 12, 6, 12],
     [((3, 2), ((0, 0), (6, 3), (2, 5)), (R(5, 12), R(11, 24), R(1, 8)))]),
    ([(1, 0), (5, 1), (4, 5), (0, 4), (2, 2), (3, 3)], [(1, 0), (5, 1), (4, 5), (0, 4)], 17, [17, 17, 17, 17],
     [((2, 2), ((1, 0), (4, 5), (0, 4)), (R(10, 17), R(6, 17), R(1, 17))),
      ((3, 3), ((1, 0), (5, 1), (4, 5)), (R(6, 17), R(1, 17), R(10, 17)))]),
]:
    H = enveloppe(pts)
    assert set(H) == set(ordre) and aire(ordre) == ar
    n = len(ordre)
    assert [cr(ordre[i], ordre[(i + 1) % n], ordre[(i + 2) % n]) for i in range(n)] == tours
    for M, tri, co in interieurs:
        assert coords(*tri, M) == co and all(w > 0 for w in co)
# Découpages cités : 24 + 18 + 9 = 51 ; 12 + 24 + 6 = 42 ; leurre « 4 sommets » : aire sans P6 = 45/2, sans (0,3) = 18.
assert [cr((0, 0), (6, 0), (5, 4)), cr((0, 0), (5, 4), (3, 6)), cr((0, 0), (3, 6), (1, 5))] == [24, 18, 9]
assert [cr((0, 0), (4, 0), (6, 3)), cr((0, 0), (6, 3), (2, 5)), cr((0, 0), (2, 5), (0, 3))] == [12, 24, 6]
assert aire([(0, 0), (6, 0), (5, 4), (1, 5)]) == R(45, 2) and aire([(0, 0), (4, 0), (6, 3), (2, 5)]) == 18
assert det((4, 1), (-1, 4)) == 17

# geo-l2-affine-convexite-ex5 : Gauss-Lucas, contrôle numérique sur P = (X-1)(X-2i)(X+3) : racines de P' barycentres positifs.
from sympy import I, Poly, nroots  # noqa: E402
X = symbols("X")
racines = [1, 2 * I, -3]
Pp = Poly((X - 1) * (X - 2 * I) * (X + 3), X).diff(X)
for zr in nroots(Pp.as_expr()):
    w = [1 / abs(complex(zr) - complex(r)) ** 2 for r in racines]
    bar = sum(wi * complex(r) for wi, r in zip(w, racines)) / sum(w)
    assert abs(bar - complex(zr)) < 1e-10

# QCM q7 : 3/2 hors de [0,1] ∪ [2,3].
assert not (0 <= F(3, 2) <= 1 or 2 <= F(3, 2) <= 3)

# ===========================================================================
# geo-l2-affine-optimisation
# ===========================================================================
# Exemple : x + 2y aux sommets (0,0), (3,0), (2,2), (0,3) : 0, 3, 6, 6 ; quadrilatère convexe.
Qd = [(0, 0), (3, 0), (2, 2), (0, 3)]
assert [p[0] + 2 * p[1] for p in Qd] == [0, 3, 6, 6] and all(cr(Qd[i], Qd[(i + 1) % 4], Qd[(i + 2) % 4]) > 0 for i in range(4))
# Affirmation 3 : (x-1/3)^2+(y-1/3)^2 vaut 2/9, 5/9, 5/9 aux sommets.
assert [(F(p[0]) - F(1, 3)) ** 2 + (F(p[1]) - F(1, 3)) ** 2 for p in [(0, 0), (1, 0), (0, 1)]] == [F(2, 9), F(5, 9), F(5, 9)]


# geo-l2-affine-optimisation-ex1 et ex2 (+ variantes) : sommets admissibles par force brute, maximum, leurres.
def sommets(contraintes):
    """contraintes : liste de (p, q, r) pour p x + q y <= r ; renvoie les sommets admissibles."""
    S = set()
    for (p1, q1, r1), (p2, q2, r2) in combinations(contraintes, 2):
        D = p1 * q2 - p2 * q1
        if D == 0:
            continue
        P = (F(r1 * q2 - r2 * q1, D), F(p1 * r2 - p2 * r1, D))
        if all(p * P[0] + q * P[1] <= r for p, q, r in contraintes):
            S.add(P)
    return S


def verif_lp(contraintes, f, sommets_attendus, opt, point, leurres):
    S = sommets(contraintes)
    assert S == {tuple(map(F, v)) for v in sommets_attendus}
    vals = {P: f(*P) for P in S}
    assert max(vals.values()) == opt and vals[tuple(map(F, point))] == opt
    assert list(vals.values()).count(opt) == 1
    for (L, valeur) in leurres:  # chaque leurre : point non optimal admissible, ou point non admissible
        admissible = all(p * L[0] + q * L[1] <= r for p, q, r in contraintes)
        assert f(*L) == valeur and (not admissible or valeur < opt)


pos = [(-1, 0, 0), (0, -1, 0)]
verif_lp(pos + [(1, 1, 4), (1, 3, 6)], lambda x, y: 3 * x + 2 * y, [(0, 0), (4, 0), (3, 1), (0, 2)], 12, (4, 0),
         [((3, 1), 11), ((6, 0), 18), ((0, 4), 8)])
verif_lp(pos + [(1, 1, 5), (2, 1, 8), (0, 1, 3)], lambda x, y: 2 * x + 5 * y, [(0, 0), (4, 0), (3, 2), (2, 3), (0, 3)], 19, (2, 3),
         [((3, 2), 16), ((0, 5), 25), ((F(5, 2), 3), 20)])
verif_lp(pos + [(2, 1, 10), (1, 2, 8)], lambda x, y: x + 3 * y, [(0, 0), (5, 0), (4, 2), (0, 4)], 12, (0, 4),
         [((4, 2), 10), ((0, 10), 30), ((8, 0), 8)])
verif_lp(pos + [(2, 1, 100), (1, 2, 80)], lambda x, y: 30 * x + 40 * y, [(0, 0), (50, 0), (40, 20), (0, 40)], 2000, (40, 20),
         [((50, 0), 1500), ((0, 40), 1600), ((20, 40), 2200)])
verif_lp(pos + [(3, 1, 90), (1, 1, 50)], lambda x, y: 50 * x + 30 * y, [(0, 0), (30, 0), (20, 30), (0, 50)], 1900, (20, 30),
         [((30, 0), 1500), ((0, 50), 1500), ((30, 20), 2100)])
verif_lp(pos + [(1, 3, 90), (1, 1, 50)], lambda x, y: 20 * x + 50 * y, [(0, 0), (50, 0), (30, 20), (0, 30)], 1600, (30, 20),
         [((0, 30), 1500), ((50, 0), 1000), ((20, 30), 1900)])
assert 20 + 2 * 40 == 100 > 80 and 3 * 30 + 20 == 110 > 90 and 20 + 3 * 30 == 110 > 90  # leurres non réalisables cités

# geo-l2-affine-optimisation-ex3 (+ v1, v2) : max aux sommets, min par forme canonique (point dans le triangle).
for f, tri, mx, mn, pmin, leurres in [
    (lambda x, y: x**2 + y**2 - 2 * x, [(0, 0), (3, 0), (0, 2)], 4, -1, (1, 0), [(4, 0), (3, -1), (4, -2)]),
    (lambda x, y: x**2 + y**2 - 4 * y, [(0, 0), (2, 0), (0, 3)], 4, -4, (0, 2), [(4, -3), (0, -4), (13, -4)]),
    (lambda x, y: x**2 + y**2 - 2 * x - 2 * y, [(0, 0), (4, 0), (0, 4)], 8, -2, (1, 1), [(8, 0), (0, -2), (16, -2)]),
]:
    assert max(f(*p) for p in tri) == mx and f(*pmin) == mn
    co = coords(*tri, pmin)
    assert all(w >= 0 for w in co)                          # le point du minimum est dans le triangle
    # forme canonique : f - mn est un carré de distance au point pmin
    assert (f(x, y) - mn - ((x - pmin[0]) ** 2 + (y - pmin[1]) ** 2)).expand() == 0
    assert all((L[0], L[1]) != (mx, mn) for L in leurres)
assert [x**2 + y**2 - 4 * y for x, y in [(0, 0), (2, 0), (0, 3)]] == [0, 4, -3]


# geo-l2-affine-optimisation-ex4 (+ v1, v2) : fonction de Leibniz.
def leibniz(pts, poids, M):
    return sum(w * ((M[0] - p[0]) ** 2 + (M[1] - p[1]) ** 2) for p, w in zip(pts, poids))


for pts, poids, Gatt, mn, faux in [
    ([(0, 0), (4, 0), (0, 4)], [2, 1, 1], (1, 1), 24, [((F(4, 3), F(4, 3)), 24), ((1, 1), 12), ((1, 1), 0)]),
    ([(0, 0), (6, 0), (0, 4)], [1, 1, 2], (F(3, 2), 2), 43, [((2, F(4, 3)), 43), ((F(3, 2), 2), 86), ((F(3, 2), 2), 0)]),
    ([(1, 1), (6, 1), (1, 6)], [3, 1, 1], (2, 2), 40, [((F(8, 3), F(8, 3)), 40), ((2, 2), 20), ((2, 2), 0)]),
]:
    Gc = bary(pts, poids)
    assert Gc == tuple(map(F, Gatt)) and leibniz(pts, poids, Gc) == mn
    # identité de Leibniz en un point symbolique
    assert (leibniz(pts, poids, (x, y)) - sum(poids) * ((x - Gc[0]) ** 2 + (y - Gc[1]) ** 2) - mn).expand() == 0
    for (P, v) in faux:
        assert not (tuple(map(F, P)) == Gc and v == mn)
    # l'isobarycentre n'est pas le point de minimum
    iso = bary(pts, [1, 1, 1])
    assert leibniz(pts, poids, iso) > mn
assert leibniz([(0, 0), (4, 0), (0, 4)], [2, 1, 1], (F(4, 3), F(4, 3))) == F(224, 9)

# QCM q8 : 2x - y aux sommets du carré : 0, 2, 1, -1.
assert [2 * p[0] - p[1] for p in [(0, 0), (1, 0), (1, 1), (0, 1)]] == [0, 2, 1, -1]

# ===========================================================================
# Ajouts de la relecture indépendante
# ===========================================================================
# geo-l2-affine-barycentres-ex1 v1 : A(0,1), B(2,3), C(-1,0) sont alignés sur y = x+1, et G = (4,5) aussi
# (la correction ne dit plus « hors du triangle »).
assert cr((0, 1), (2, 3), (-1, 0)) == 0 and all(p[1] == p[0] + 1 for p in [(0, 1), (2, 3), (-1, 0), (4, 5)])
# geo-l2-affine-optimisation-ex1 v1 : l'intersection (0,8) de x=0 et 2x+y=8 viole y <= 3.
assert 2 * 0 + 8 == 8 and 8 > 3
# geo-l2-affine-convexite-ex3 (+ v1, v2) : déterminants de deux côtés consécutifs cités dans les corrections.
def tours(poly):
    n = len(poly)
    # det(P_{i-1}P_i, P_iP_{i+1}) = det(P_{i-1}P_i, P_{i-1}P_{i+1})
    return [cr(poly[i - 1], poly[i % n], poly[(i + 1) % n]) for i in range(1, n + 1)]
assert tours([(0, 0), (6, 0), (5, 4), (3, 6), (1, 5)]) == [24, 6, 6, 9, 30]
assert tours([(0, 0), (4, 0), (6, 3), (2, 5), (0, 3)]) == [12, 16, 12, 6, 12]
assert tours([(1, 0), (5, 1), (4, 5), (0, 4)]) == [17, 17, 17, 17]
# geo-l2-affine-applications-ex3 v2 : k1 k2 = 1, translation de vecteur (k1 - 1) O1O2 = (1,1).
assert (2 - 1) * (1 - 0) == 1 and 2 * R(1, 2) == 1
# geo-l2-affine-volume (cours) : deux bases de même orientation non orthonormées peuvent donner le même volume :
# B' = (2e1, e2/2) a pour déterminant 1 dans la base canonique.
assert det((2, 0), (0, R(1, 2))) == 1

print("geo_l2_affine : toutes les vérifications passent.")
