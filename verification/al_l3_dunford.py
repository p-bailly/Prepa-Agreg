"""Vérifications calculatoires du chapitre « al-l3-dunford » (algèbre linéaire, L3) : décomposition de Dunford et
exponentielle. Cartes, exercices et QCM de content/*/algebre-lineaire/l3-al-l3-dunford.json et content/qcm/al-l3-dunford.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/al_l3_dunford.py (doit finir sans erreur, en quelques secondes).
"""

import random
from fractions import Fraction as F

from sympy import (E, I as i_, Matrix, Rational, cos, diag, diff, dsolve, Eq, exp, expand, eye, factor, Function,
                   pi, simplify, sin, symbols, zeros, Poly, re, factorial)

X = symbols("X")
t, th = symbols("t theta", real=True)
I2, I3 = eye(2), eye(3)
E11 = Matrix([[1, 0], [0, 0]])
E12 = Matrix([[0, 1], [0, 0]])
J = Matrix([[1, 1], [0, 1]])
R = Matrix([[0, -1], [1, 0]])  # rotation d'angle pi/2 (notée R, C ou K selon les textes)


def chi(A):
    return A.charpoly(X).as_expr()


def evalp(p, A):
    """p(A) pour un polynôme p en X (Horner)."""
    n = A.shape[0]
    Res = zeros(n)
    for c in Poly(expand(p), X).all_coeffs():
        Res = Res * A + c * eye(n)
    return Res


def nilp(N):
    return N ** N.shape[0] == zeros(N.shape[0])


def dunford_ok(A, D, N):
    """(D, N) est la décomposition de Dunford de A."""
    return A == D + N and D.is_diagonalizable() and nilp(N) and D * N == N * D


def eqm(M1, M2):
    return simplify(M1 - M2) == zeros(*M1.shape)


def multiset_vp(A):
    out = []
    for v, m in A.eigenvals().items():
        out += [v] * m
    return sorted(out, key=lambda z: (re(z), str(z)))


def vp(*xs):
    return sorted(list(xs), key=lambda z: (re(z), str(z)))


# ===========================================================================
# al-l3-dunford-simultanee (carte)
# ===========================================================================
A = diag(1, 1, 2)
B = Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 5]])
# Exemple : A et B commutent ; e1 propre pour A, pas pour B ; base commune (e1+e2, e1-e2, e3).
assert A * B == B * A and B * Matrix([1, 0, 0]) == Matrix([0, 1, 0])
P = Matrix([[1, 1, 0], [1, -1, 0], [0, 0, 1]])
assert P.inv() * A * P == diag(1, 1, 2) and P.inv() * B * P == diag(1, -1, 5)
assert multiset_vp(A + B) == vp(0, 2, 7)
# Contre-exemple 1 (et trou 2, leurres « même polynôme caractéristique », « mêmes valeurs propres ») :
A1 = Matrix([[1, 1], [0, 0]])
assert A1.is_diagonalizable() and E11.is_diagonalizable() and chi(A1) == chi(E11)
assert A1 * E11 == E11 and E11 * A1 == A1 and A1 * E11 != E11 * A1
# Trou 2, leurre « somme diagonalisable » : A1 + E11 = [[2,1],[0,0]] est diagonalisable, sans base commune.
assert (A1 + E11).is_diagonalizable()
# Contre-exemple 2 : A1 + diag(-1,0) nilpotente non nulle.
assert A1 + diag(-1, 0) == E12 and not E12.is_diagonalizable()
# Contre-exemple 3, piège 1, trou 1 (leurres trigonalisables, inversibles) : I2 et J commutent, J non diagonalisable.
assert I2 * J == J * I2 and not J.is_diagonalizable() and J.det() != 0
# Trou 1, leurre « nilpotents » : 0 et E12 commutent, E12 non diagonalisable.
assert zeros(2) * E12 == E12 * zeros(2) and nilp(E12)
# Piège 2 : A1 et E11 ne commutent pas, mais tous les produits A1 E11, E11 A1, A1^2, E11^2 sont diagonalisables.
assert all(M.is_diagonalizable() for M in (A1 * E11, E11 * A1, A1 * A1, E11 * E11))
# Affirmation 2 : u = E12, v = -E12 commutent, u + v = 0.
assert E12 * (-E12) == (-E12) * E12 and E12 + (-E12) == zeros(2)
# Affirmation 3 : I2 commute avec J ; J e2 = e1 + e2.
assert J * Matrix([0, 1]) == Matrix([1, 1])

# --- al-l3-dunford-simultanee-ex1 (et v1, v2) : A, B, v1, v2, valeurs propres de AB, leurres.
for (A, B, v1, v2, a, b, AB, bonne, leurres) in [
    ([[1, 2], [2, 1]], [[3, -1], [-1, 3]], (1, 1), (1, -1), (3, -1), (2, 4), [[1, 5], [5, 1]], (6, -4), [(12, -2), (5, 3), (6, 4)]),
    ([[7, -3], [6, -2]], [[-7, 5], [-10, 8]], (1, 2), (1, 1), (1, 4), (3, -2), [[-19, 11], [-22, 14]], (3, -8), [(-2, 12), (4, 2), (3, 8)]),
    ([[7, -10], [5, -8]], [[9, -8], [4, -3]], (2, 1), (1, 1), (2, -3), (5, 1), [[23, -26], [13, -16]], (10, -3), [(2, -15), (7, -2), (10, 3)]),
]:
    A, B, v1, v2, AB = Matrix(A), Matrix(B), Matrix(v1), Matrix(v2), Matrix(AB)
    assert A * B == B * A == AB
    assert A * v1 == a[0] * v1 and A * v2 == a[1] * v2 and B * v1 == b[0] * v1 and B * v2 == b[1] * v2
    assert multiset_vp(A * B) == vp(*bonne)
    assert multiset_vp(A + B) == vp(a[0] + b[0], a[1] + b[1])
    for l in leurres:
        assert vp(*l) != multiset_vp(A * B)

# --- al-l3-dunford-simultanee-ex2 (et v1, v2)
for (A, B, base, dA, dB, vpsum, leurres, pasvp) in [
    ([[1, 0, 0], [0, 1, 0], [1, -1, 2]], [[-3, 3, 0], [-6, 6, 0], [0, 0, 3]], [(1, 1, 0), (1, 2, 1), (0, 0, 1)], (1, 1, 2), (0, 3, 3),
     (1, 4, 5), [(2, 4, 4), (-2, 7, 5), (1, 3, 5)], (-1, 0, 1)),
    ([[2, 3, -3], [0, 2, 0], [0, 3, -1]], [[-1, -2, 2], [-2, -1, 2], [-2, -2, 3]], [(0, 1, 1), (1, 1, 1), (1, 0, 1)], (2, 2, -1), (1, -1, 1),
     (0, 1, 3), [(-2, 3, 3), (1, 1, 2), (0, 2, 3)], (1, 0, 0)),
    ([[3, 3, -3], [3, 3, -3], [3, 3, -3]], [[1, 0, 0], [-1, 1, 1], [-1, 0, 2]], [(1, 0, 1), (0, 1, 1), (1, 1, 1)], (0, 0, 3), (1, 2, 1),
     (1, 2, 4), [(1, 1, 5), (-1, 4, 4), (2, 2, 3)], (1, -1, 0)),
]:
    A, B = Matrix(A), Matrix(B)
    assert A * B == B * A and A.is_diagonalizable() and B.is_diagonalizable()
    P = Matrix.hstack(*[Matrix(v) for v in base])
    assert P.det() != 0 and P.inv() * A * P == diag(*dA) and P.inv() * B * P == diag(*dB)
    assert multiset_vp(A + B) == vp(*vpsum)
    for l in leurres:
        assert vp(*l) != vp(*vpsum)
    # un vecteur du sous-espace propre double de A qui n'est pas propre pour B
    w = Matrix(pasvp)
    lam2 = [k for k, m in A.eigenvals().items() if m == 2][0]
    assert A * w == lam2 * w and Matrix.hstack(w, B * w).rank() == 2
# rangs et traces cités dans les corrections
assert (Matrix([[1, 0, 0], [0, 1, 0], [1, -1, 2]]) - I3).rank() == 1
assert (Matrix([[2, 3, -3], [0, 2, 0], [0, 3, -1]]) - 2 * I3).rank() == 1 and (Matrix([[-1, -2, 2], [-2, -1, 2], [-2, -2, 3]]) - I3).rank() == 1
assert (Matrix([[1, 0, 0], [-1, 1, 1], [-1, 0, 2]]) - I3).rank() == 1

# --- al-l3-dunford-simultanee-ex3 (et v1, v2) : nombre de X réelles avec X^2 = A.
a, b, c = symbols("a b c")
for A, nb, vps in [([[1, 1, 0], [0, 4, 1], [0, 0, 9]], 8, (1, 4, 9)), ([[4, 1, 1], [0, 0, 1], [0, 0, 1]], 4, (4, 0, 1)),
                   ([[1, 2, 0], [0, -1, 1], [0, 0, 4]], 0, (1, -1, 4))]:
    A = Matrix(A)
    assert multiset_vp(A) == vp(*vps)
    # X commute avec A (valeurs propres distinctes) : X est un polynôme de degré <= 2 en A.
    Xs = symbols("x0:9")
    Xm = Matrix(3, 3, Xs)
    from sympy import linsolve
    sol = linsolve(list(Xm * A - A * Xm), Xs)
    assert len(list(sol)[0].free_symbols) == 3  # commutant de dimension 3 = Vect(I, A, A^2)
    from sympy import solve
    M = a * I3 + b * A + c * A * A
    sols = solve(list(M * M - A), [a, b, c], dict=True)
    reelles = [s for s in sols if all(v.is_real for v in s.values())]
    assert len(reelles) == nb
    # les leurres « 8 », « 4 », « 0 », « une infinité » : un seul est juste
    assert nb in (0, 4, 8)

# --- al-l3-dunford-simultanee-ex4 : I2 et R commutent, R sans vecteur propre réel ; E11, E12.
assert I2 * R == R * I2 and factor(chi(R)) == X ** 2 + 1
assert E11 * E12 == E12 and E12 * E11 == zeros(2)  # réciproques fausses (vecteur propre commun e1)
assert E11 * Matrix([1, 0]) == Matrix([1, 0]) and E12 * Matrix([1, 0]) == zeros(2, 1)

# --- al-l3-dunford-simultanee-ex5 : le groupe des diag(±1) a 2^m éléments ; leurres jamais égaux à 2^m.
import itertools
from math import factorial
for m in range(1, 5):
    G = [diag(*s) for s in itertools.product([1, -1], repeat=m)]
    assert len(G) == 2 ** m and all(g * g == eye(m) for g in G)
    assert all(g * h in G for g in G for h in G)  # sous-groupe
for m in range(1, 30):
    assert factorial(m) != 2 ** m and 2 ** (m - 1) != 2 ** m

# ===========================================================================
# al-l3-dunford-dunford (carte)
# ===========================================================================
A = Matrix([[1, 1, 1], [0, 1, 1], [0, 0, 2]])
assert factor(chi(A)) == (X - 2) * (X - 1) ** 2
assert expand((X - 1) ** 2 - X * (X - 2)) == 1  # Bézout de l'exemple
p2 = (A - I3) ** 2
D = I3 + p2
N = A - D
assert D == Matrix([[1, 0, 2], [0, 1, 1], [0, 0, 2]]) and N == Matrix([[0, 1, -1], [0, 0, 0], [0, 0, 0]])
assert N ** 2 == zeros(3) and dunford_ok(A, D, N)
# la partie diagonale ne convient pas
Sd = A - diag(1, 1, 2)
assert diag(1, 1, 2) * Sd != Sd * diag(1, 1, 2)
# Contre-exemple « partie diagonale » et affirmation 2, q3 : [[1,1],[0,2]] diagonalisable.
A2 = Matrix([[1, 1], [0, 2]])
assert A2.is_diagonalizable() and diag(1, 2) * E12 == E12 and E12 * diag(1, 2) == 2 * E12
assert dunford_ok(A2, A2, zeros(2))
# Trou 3, leurres « d∘n nilpotent » et « tr(d∘n) = 0 » : deux écritures de A2 qui les vérifient.
for (d, n) in [(A2, zeros(2)), (diag(1, 2), E12)]:
    assert d + n == A2 and d.is_diagonalizable() and nilp(n) and nilp(d * n) and (d * n).trace() == 0
# Trou 3, leurre « d∘n = 0 » : pour J, n nilpotente non nulle de taille 2 vérifie n^2 = 0, donc d n = J n ≠ 0.
p_, q_, r_ = symbols("p q r")
n_gen = Matrix([[p_, q_], [r_, -p_]])
assert expand((n_gen ** 2 - (p_ ** 2 + q_ * r_) * I2)) == zeros(2)  # n^2 = -det(n) I, nul si n nilpotente
assert J.det() != 0  # J n = 0 impose n = 0, et alors d = J non diagonalisable
# Trou 1 : leurre « trigonalisable » (J = J + 0 = I + E12), « inversible » (u = 0), « nilpotent » (u = I).
assert dunford_ok(J, I2, E12) and J + zeros(2) == J
# Trou 2 : « diagonalisable » (I = I + 0 = 0 + I), « rang au plus 1 » (diag(1,0)), « trace nulle » (diag(1,-1)).
assert diag(1, 0).rank() == 1 and diag(1, -1).trace() == 0
# Contre-exemple « polynôme non scindé » : chi_R = X^2 + 1 sans racine réelle.
assert factor(chi(R)) == X ** 2 + 1
# Affirmation 4 : polynômes minimaux de J et de I2.
assert (J - I2) != zeros(2) and (J - I2) ** 2 == zeros(2)

# --- al-l3-dunford-dunford-ex1 (et v1, v2)
for A, lam, Nm, dg in [([[3, 1], [-1, 1]], 2, [[1, 1], [-1, -1]], (3, 1)), ([[5, -1], [4, 1]], 3, [[2, -1], [4, -2]], (5, 1)),
                       ([[-1, 4], [-1, 3]], 1, [[-2, 4], [-1, 2]], (-1, 3))]:
    A, Nm = Matrix(A), Matrix(Nm)
    assert factor(chi(A)) == (X - lam) ** 2
    assert dunford_ok(A, lam * I2, Nm) and not A.is_diagonalizable()
    R0 = A - diag(*dg)
    assert not nilp(R0) and R0 ** 2 == R0[0, 1] * R0[1, 0] * I2  # la partie diagonale ne convient pas (leurre)
    lower = Matrix([[A[0, 0], 0], [A[1, 0], A[1, 1]]])
    assert lower != lam * I2 and A != lam * I2
assert (Matrix([[0, 1], [-1, 0]])) ** 2 == -I2 and Matrix([[0, -1], [4, 0]]) ** 2 == -4 * I2 and Matrix([[0, 4], [-1, 0]]) ** 2 == -4 * I2

# --- al-l3-dunford-dunford-ex2 (et v1, v2) : triangulaire à diagonale distincte, N = 0.
for A, dg, N0, P1, P2 in [
    ([[1, 2, 0], [0, 3, 1], [0, 0, 5]], (1, 3, 5), [[0, 2, 0], [0, 0, 1], [0, 0, 0]], [[0, 2, 0], [0, 0, 3], [0, 0, 0]], [[0, 6, 0], [0, 0, 5], [0, 0, 0]]),
    ([[2, -1, 3], [0, -1, 1], [0, 0, 4]], (2, -1, 4), [[0, -1, 3], [0, 0, 1], [0, 0, 0]], [[0, -2, 6], [0, 0, -1], [0, 0, 0]], [[0, 1, 12], [0, 0, 4], [0, 0, 0]]),
    ([[0, 1, 1], [0, 1, 2], [0, 0, -1]], (0, 1, -1), [[0, 1, 1], [0, 0, 2], [0, 0, 0]], [[0, 0, 0], [0, 0, 2], [0, 0, 0]], [[0, 1, -1], [0, 0, -2], [0, 0, 0]]),
]:
    A, N0, P1, P2 = Matrix(A), Matrix(N0), Matrix(P1), Matrix(P2)
    assert A == diag(*dg) + N0 and nilp(N0)
    assert diag(*dg) * N0 == P1 and N0 * diag(*dg) == P2 and P1 != P2
    assert dunford_ok(A, A, zeros(3))

# --- al-l3-dunford-dunford-ex3 (et v1, v2) : D = Q(A), Q de degré <= 2 ; leurres Q(A) != D.
for A, lam, mu, rg, S, D, N, Q, leurres in [
    ([[-1, 0, 0], [1, 3, -1], [-2, 1, 1]], 2, -1, 2, [[9, 0, 0], [0, 0, 0], [9, 0, 0]], [[-1, 0, 0], [0, 2, 0], [-3, 0, 2]],
     [[0, 0, 0], [1, 1, -1], [1, 1, -1]], 2 - (X - 2) ** 2 / 3, [2 + (X - 2) ** 2 / 3, 2 - (X - 2) ** 2 / 9, 2 - (X - 2) ** 2]),
    ([[1, 0, 0], [-1, -1, 1], [0, -1, 1]], 0, 1, 2, [[1, 0, 0], [0, 0, 0], [1, 0, 0]], [[1, 0, 0], [0, 0, 0], [1, 0, 0]],
     [[0, 0, 0], [-1, -1, 1], [-1, -1, 1]], X ** 2, [X ** 2 / 2, -X ** 2, 2 * X ** 2]),  # leurres revus à la relecture
    ([[1, 0, 0], [-3, 3, 1], [-2, 0, 3]], 3, 1, 2, [[4, 0, 0], [4, 0, 0], [4, 0, 0]], [[1, 0, 0], [-2, 3, 0], [-2, 0, 3]],
     [[0, 0, 0], [-1, 0, 1], [0, 0, 0]], 3 - (X - 3) ** 2 / 2, [3 + (X - 3) ** 2 / 2, 3 - (X - 3) ** 2 / 4, 3 - (X - 3) ** 2]),
]:
    A, S, D, N = Matrix(A), Matrix(S), Matrix(D), Matrix(N)
    assert expand(chi(A) - (X - lam) ** 2 * (X - mu)) == 0
    assert (A - lam * I3).rank() == rg and not A.is_diagonalizable()
    assert (A - lam * I3) ** 2 == S and dunford_ok(A, D, N)
    assert Q.subs(X, mu) == mu and expand(Q - lam).subs(X, lam) == 0 and diff(Q, X).subs(X, lam) == 0
    assert evalp(Q, A) == D
    for Ql in leurres:
        assert evalp(Ql, A) != D
    # polynôme minimal = chi (A non diagonalisable) : (A - lam)(A - mu) != 0
    assert (A - lam * I3) * (A - mu * I3) != zeros(3)
# rangs cités : colonnes opposées (ex3 et v1), colonne nulle (v2)
M_ = Matrix([[-1, 0, 0], [1, 3, -1], [-2, 1, 1]]) - 2 * I3
assert M_[:, 1] == -M_[:, 2]
M_ = Matrix([[1, 0, 0], [-1, -1, 1], [0, -1, 1]])
assert M_[:, 1] == -M_[:, 2]
M_ = Matrix([[1, 0, 0], [-3, 3, 1], [-2, 0, 3]]) - 3 * I3
assert M_[:, 1] == zeros(3, 1)

# --- al-l3-dunford-dunford-ex4 : u = E12 : u^2 = 0 diagonalisable, u non diagonalisable ; leurres réfutés.
assert E12 ** 2 == zeros(2) and not E12.is_diagonalizable()
U = diag(1, 0)
assert (U ** 2).is_diagonalizable() and U.det() == 0 and not nilp(U)  # « nilpotent ou inversible » faux
# u^2 = d^2 + (2d + n) n pour une décomposition de Dunford (exemple de la carte)
A = Matrix([[1, 1, 1], [0, 1, 1], [0, 0, 2]])
D = Matrix([[1, 0, 2], [0, 1, 1], [0, 0, 2]]); N = A - D
assert A ** 2 == D ** 2 + (2 * D + N) * N and dunford_ok(A ** 2, D ** 2, (2 * D + N) * N)

# --- al-l3-dunford-dunford-ex5 : conjugaison par diag(1, eps, eps^2).
eps = symbols("epsilon", positive=True)
T = Matrix(3, 3, lambda i, j: symbols(f"t{i}{j}") if j >= i else 0)
De = diag(1, eps, eps ** 2)
C = De.inv() * T * De
assert all(simplify(C[i, j] - eps ** (j - i) * T[i, j]) == 0 for i in range(3) for j in range(3))
assert C.subs(eps, 0) == diag(T[0, 0], T[1, 1], T[2, 2])
Jn = E12
assert diag(1, eps).inv() * Jn * diag(1, eps) == Matrix([[0, eps], [0, 0]])

# ===========================================================================
# al-l3-dunford-exponentielle (carte)
# ===========================================================================
A = Matrix([[1, 1, 1], [0, 1, 1], [0, 0, 2]])
p2 = (A - I3) ** 2
p1 = I3 - p2
eD = E * p1 + E ** 2 * p2
assert eqm(eD, Matrix([[E, 0, 2 * E ** 2 - 2 * E], [0, E, E ** 2 - E], [0, 0, E ** 2]]))
eA = Matrix([[E, E, 2 * E ** 2 - 3 * E], [0, E, E ** 2 - E], [0, 0, E ** 2]])
assert eqm(eD * (I3 + (A - I3 - p2)), eA) and eqm(A.exp(), eA)
assert simplify(eA.det() - exp(A.trace())) == 0
# Contre-exemple « sans commutation » : A1 = diag(1,0) + E12, A1^2 = A1.
A1 = Matrix([[1, 1], [0, 0]])
assert A1 ** 2 == A1 and eqm(A1.exp(), Matrix([[E, E - 1], [0, 1]]))
assert eqm(diag(1, 0).exp() * E12.exp(), Matrix([[E, E], [0, 1]])) and not eqm(A1.exp(), diag(1, 0).exp() * E12.exp())
# Trou 1, pièges 1 : e^{d+n} != e^d e^n sans commutation (même exemple).
# Contre-exemple / affirmation 2 / pièges 2, 3 : B = 2 pi R, e^B = I2.
B2 = 2 * pi * R
assert eqm(B2.exp(), I2) and not B2.is_diagonalizable(reals_only=True) and B2.is_diagonalizable()
# piège 3 : u = diag(0) + 2 pi R sur R^3 : chi a une racine réelle, e^u = I3, u non diagonalisable sur R.
U3 = diag(Matrix([[0]]), B2)
assert eqm(U3.exp(), I3) and not U3.is_diagonalizable(reals_only=True)
assert factor(chi(U3)) == X * (X ** 2 + 4 * pi ** 2)
# piège 2 (K = R, chi de e^u scindé) : chi de e^B = (X - 1)^2 scindé sur R, chi_B = X^2 + 4 pi^2 ne l'est pas.
assert expand(chi(B2.exp().applyfunc(simplify)) - (X - 1) ** 2) == 0 and factor(chi(B2)) == X ** 2 + 4 * pi ** 2
# al-l3-dunford-dunford, piège 1 (K = R, chi_u a une racine réelle) : u = diag(1) + R sur R^3, chi = (X-1)(X^2+1).
# Si u = d + n (d diagonalisable sur R, n nilpotent, dn = nd), n laisse stables les sous-espaces propres de d, sur
# lesquels u = lambda id + n : chi_u serait scindé sur R. Donc u n'a pas de décomposition de Dunford réelle.
assert factor(chi(diag(Matrix([[1]]), R))) == (X - 1) * (X ** 2 + 1)
# Trou 1, leurres : e^0 + e^0 = 2I ; e^I = e I mais e^0 = I ; pour J : e^J = e(I+E12) != e I + E12.
assert eqm(J.exp(), E * (I2 + E12)) and not eqm(J.exp(), E * I2 + E12)
# Trou 2, leurres : u = I (e^1 != 1, e != 2) ; u = [[0,1],[1,0]] : diagonale nulle, valeurs propres de e^u : e, 1/e.
S2 = Matrix([[0, 1], [1, 0]])
assert set(S2.exp().eigenvals().keys()) == {E, exp(-1)}
# Trou 3, leurres : u = I (d = I != 0, e^u diagonalisable) ; u = J (n^2 = 0, J inversible, e^J non diagonalisable).
assert not J.exp().is_diagonalizable()
# Affirmation 3 : u nilpotente => e^u - I nilpotente.
Nn = Matrix([[0, 1, 2], [0, 0, 3], [0, 0, 0]])
assert nilp(Nn.exp() - I3)
# Affirmation 4 : det e^A = e^{tr A} > 0.
assert simplify(Matrix([[1, 2], [3, 4]]).exp().det() - exp(5)) == 0

# --- al-l3-dunford-exponentielle-ex1 (et v1, v2) : det e^A = e^{tr A} ; leurres.
for A, bonne, leurres in [([[1, 2], [3, 4]], exp(5), [exp(-2), exp(10), E + exp(4)]),
                          ([[0, 3], [1, -2]], exp(-2), [exp(-3), exp(2), 1 + exp(-2)]),
                          ([[2, 1, 0], [-1, 0, 3], [4, 1, -1]], E, [exp(5), exp(9), exp(2) + 1 + exp(-1)])]:
    A = Matrix(A)
    assert simplify(exp(A.trace()) - bonne) == 0
    assert exp(A.det()) in leurres  # leurre e^{det A}
    assert exp(sum(A)) in leurres  # leurre e^{somme des coefficients}
    assert simplify(sum(exp(A[k, k]) for k in range(A.shape[0])) - [l for l in leurres if l.is_Add][0]) == 0
    for l in leurres:
        assert simplify(l - bonne) != 0
assert simplify(Matrix([[1, 2], [3, 4]]).exp().det() - exp(5)) == 0

# --- al-l3-dunford-exponentielle-ex2 (et v1, v2) : coefficient (1,3) de e^{tA}, A = lambda I + N.
for lam, Nm, c13, leurres in [
    (-1, [[0, 1, 1], [0, 0, 2], [0, 0, 0]], (t + t ** 2) * exp(-t), [(t + 2 * t ** 2) * exp(-t), t * exp(-t), (1 + t + t ** 2) * exp(-t)]),
    (1, [[0, 2, -1], [0, 0, 1], [0, 0, 0]], (t ** 2 - t) * exp(t), [(2 * t ** 2 - t) * exp(t), -t * exp(t), (1 - t + t ** 2) * exp(t)]),
    (2, [[0, 3, 1], [0, 0, 2], [0, 0, 0]], (t + 3 * t ** 2) * exp(2 * t), [(t + 6 * t ** 2) * exp(2 * t), t * exp(2 * t), (1 + t + 3 * t ** 2) * exp(2 * t)]),
]:
    Nm = Matrix(Nm)
    A = lam * I3 + Nm
    assert dunford_ok(A, lam * I3, Nm)
    EtA = exp(lam * t) * (I3 + t * Nm + t ** 2 / 2 * Nm ** 2)
    assert eqm((t * A).exp(), EtA) and simplify(EtA[0, 2] - c13) == 0
    for l in leurres:
        assert simplify(l - c13) != 0

# --- al-l3-dunford-exponentielle-ex3 (et v1, v2) : e^A par Dunford ; leurres.
for A, lam, mu, Dfor, D, N, eDm, eAm, leurres in [
    ([[0, 0, 0], [-1, 2, -1], [0, 1, 0]], 1, 0, None, [[0, 0, 0], [0, 1, 0], [1, 0, 1]], [[0, 0, 0], [-1, 1, -1], [-1, 1, -1]],
     [[1, 0, 0], [0, E, 0], [E - 1, 0, E]], [[1, 0, 0], [-E, 2 * E, -E], [-1, E, 0]],
     [[[1, 1, 1], [exp(-1), exp(2), exp(-1)], [1, E, 1]], [[1, 0, 0], [-1, E + 1, -1], [E - 2, 1, E - 1]]]),
    ([[1, 0, 0], [1, 0, 1], [-3, -1, -2]], -1, 1, None, [[1, 0, 0], [0, -1, 0], [-2, 0, -1]], [[0, 0, 0], [1, 1, 1], [-1, -1, -1]],
     [[E, 0, 0], [0, exp(-1), 0], [exp(-1) - E, 0, exp(-1)]], [[E, 0, 0], [exp(-1), 2 * exp(-1), exp(-1)], [-E, -exp(-1), 0]],
     [[[E, 1, 1], [E, 1, E], [exp(-3), exp(-1), exp(-2)]], [[E, 0, 0], [1, exp(-1) + 1, 1], [exp(-1) - E - 1, -1, exp(-1) - 1]]]),
    ([[0, 0, 0], [-1, 1, -1], [3, 1, 3]], 2, 0, None, [[0, 0, 0], [0, 2, 0], [2, 0, 2]], [[0, 0, 0], [-1, -1, -1], [1, 1, 1]],
     [[1, 0, 0], [0, exp(2), 0], [exp(2) - 1, 0, exp(2)]], [[1, 0, 0], [-exp(2), 0, -exp(2)], [2 * exp(2) - 1, exp(2), 2 * exp(2)]],
     [[[1, 1, 1], [exp(-1), E, exp(-1)], [exp(3), E, exp(3)]], [[1, 0, 0], [-1, exp(2) - 1, -1], [exp(2), 1, exp(2) + 1]]]),
]:
    A, D, N, eDm, eAm = Matrix(A), Matrix(D), Matrix(N), Matrix(eDm), Matrix(eAm)
    assert expand(chi(A) - (X - lam) ** 2 * (X - mu)) == 0
    assert dunford_ok(A, D, N) and N ** 2 == zeros(3)
    assert eqm(D.exp(), eDm) and eqm(A.exp(), eAm) and eqm(eDm * (I3 + N), eAm)
    assert simplify(eAm.det() - exp(A.trace())) == 0
    # leurres : e^D (sans N), exponentielle coefficient par coefficient, e^D + N
    assert eqm(A.applyfunc(exp), Matrix(leurres[0])) and eqm(eDm + N, Matrix(leurres[1]))
    for L in [eDm, Matrix(leurres[0]), Matrix(leurres[1])]:
        assert not eqm(L, eAm)
# formules de D : I - (A-I)^2 ; -I + (A+I)^2/2 ; 2I - (A-2I)^2/2
assert Matrix([[0, 0, 0], [0, 1, 0], [1, 0, 1]]) == I3 - (Matrix([[0, 0, 0], [-1, 2, -1], [0, 1, 0]]) - I3) ** 2
assert Matrix([[1, 0, 0], [0, -1, 0], [-2, 0, -1]]) == -I3 + (Matrix([[1, 0, 0], [1, 0, 1], [-3, -1, -2]]) + I3) ** 2 / 2
assert Matrix([[0, 0, 0], [0, 2, 0], [2, 0, 2]]) == 2 * I3 - (Matrix([[0, 0, 0], [-1, 1, -1], [3, 1, 3]]) - 2 * I3) ** 2 / 2

# --- al-l3-dunford-exponentielle-ex4 : e^{theta K} ; theta = 2 pi : I2 ; theta = pi : -I2 ; leurres.
K = R
assert eqm((th * K).exp(), Matrix([[cos(th), -sin(th)], [sin(th), cos(th)]]))
assert eqm((2 * pi * K).exp(), I2) and eqm((pi * K).exp(), -I2) and not eqm((pi / 2 * K).exp(), I2)

# --- al-l3-dunford-exponentielle-ex5 : exp(L(1)) = I + N pour N nilpotente ; non-surjectivité réelle.
for Nm in [Matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]]), Matrix([[0, 2, -1], [0, 0, 3], [0, 0, 0]])]:
    m = Nm.shape[0]
    L = sum(((-1) ** (k + 1) * Nm ** k / k for k in range(1, m)), zeros(m))
    assert eqm(L.exp(), eye(m) + Nm)
    Lt = sum(((-1) ** (k + 1) * t ** k * Nm ** k / k for k in range(1, m)), zeros(m))
    assert eqm(Lt.diff(t) * (eye(m) + t * Nm), Nm)
    # correction revue : L(t)^m = 0, L'(t) L(t)^(m-1) = 0, d/dt e^{L(t)} = L'(t) e^{L(t)} (somme finie),
    # (I + tN)^{-1} = sum_{j<m} (-tN)^j et sa dérivée -N (I + tN)^{-2}, g(t) = (I + tN)^{-1} e^{L(t)} constante.
    assert Lt ** m == zeros(m) and (Lt.diff(t) * Lt ** (m - 1)).applyfunc(expand) == zeros(m)
    eLt = sum((Lt ** k / factorial(k) for k in range(m)), zeros(m))
    assert (eLt.diff(t) - Lt.diff(t) * eLt).applyfunc(expand) == zeros(m)
    inv_t = sum(((-t * Nm) ** j for j in range(m)), zeros(m))
    assert ((eye(m) + t * Nm) * inv_t).applyfunc(expand) == eye(m)
    assert (inv_t.diff(t) + Nm * inv_t * inv_t).applyfunc(expand) == zeros(m)
    assert (inv_t * eLt).applyfunc(expand) == eye(m)
assert Matrix([[-1, 0], [0, 2]]).det() < 0 and eqm((pi * K).exp(), -I2)
assert eqm(diag(Rational(1), 1).exp(), diag(E, E)) and eqm(E12.exp(), J)

# ===========================================================================
# al-l3-dunford-densite (carte)
# ===========================================================================
k = symbols("k", positive=True)
Jk = Matrix([[1 / k, 1], [0, 0]])
assert Jk.eigenvals() == {1 / k: 1, 0: 1} and (E12 - I2 / k).det() == 1 / k ** 2
# chi_AB = chi_BA sur des matrices aléatoires (dont non inversibles)
random.seed(1)
for _ in range(20):
    Am = Matrix(3, 3, lambda i, j: random.randint(-2, 2))
    Bm = Matrix(3, 3, lambda i, j: random.randint(-2, 2))
    assert expand(chi(Am * Bm) - chi(Bm * Am)) == 0
# Contre-exemple sur R : discriminant < 0 dans la boule de rayon 1/4 (échantillonnage exact) autour de R.
for _ in range(2000):
    aa, bb, cc, dd = [F(random.randint(-249, 249), 1000) for _ in range(4)]
    a_, b_, c_, d_ = aa, -1 + bb, 1 + cc, dd
    assert (a_ - d_) ** 2 + 4 * b_ * c_ < 0
# borne : 1/4 - 9/4 < 0
assert F(1, 4) - F(9, 4) < 0
# Trou 3 (m >= 2), piège 2 : R est diagonalisable sur C mais pas trigonalisable sur R.
assert R.is_diagonalizable() and not R.is_diagonalizable(reals_only=True) and factor(chi(R)) == X ** 2 + 1
# Démonstration : chi de diag(R, I_{m-2}).
assert factor(chi(diag(R, eye(2)))) == (X - 1) ** 2 * (X ** 2 + 1)
# Affirmation 5 : I/k inversible tend vers 0.

# --- al-l3-dunford-densite-ex1 (et v1, v2) : racines de det(A + t I).
for A, racines, leurres in [([[1, 2], [2, 4]], {0, -5}, [{0, 5}, {-5}, {0}]), ([[2, -1], [-4, 2]], {0, -4}, [{0, 4}, {-4}, {0}]),
                            ([[1, -3], [2, -6]], {0, 5}, [{0, -5}, {5}, {0}])]:
    A = Matrix(A)
    from sympy import solve as _solve
    s = set(_solve((A + t * I2).det(), t))
    assert s == racines and all(l != racines for l in leurres)
    assert all(1 / F(kk) not in racines for kk in range(1, 50))
assert expand((1 + t) * (4 + t) - 4 - t * (t + 5)) == 0 and expand((2 + t) ** 2 - 4 - t * (t + 4)) == 0
assert expand((1 + t) * (t - 6) + 6 - t * (t - 5)) == 0

# --- al-l3-dunford-densite-ex2 (et v1, v2) : rayon de la boule d'inversibles autour de diag(a, b), norme max.
r = symbols("r")
for (a0, b0, rho, leurres) in [(1, 2, F(2, 3), [1, F(1, 3), F(1, 2)]), (1, 1, F(1, 2), [1, F(1, 4), F(2, 3)]), (2, 3, F(6, 5), [2, F(3, 5), F(5, 6)])]:
    borne = (a0 - r) * (b0 - r) - r ** 2
    assert expand(borne - (a0 * b0 - (a0 + b0) * r)) == 0
    assert borne.subs(r, rho) == 0 and rho < min(a0, b0)
    Hs = Matrix([[-rho, rho], [rho, -rho]])
    assert (diag(a0, b0) + Hs).det() == 0
    # aucune matrice singulière strictement plus près (échantillonnage exact sur la boule de rayon < rho)
    for _ in range(400):
        h = [F(random.randint(-999, 999), 1000) * rho for _ in range(4)]
        assert (a0 + h[0]) * (b0 + h[3]) - h[1] * h[2] > 0
    assert all(F(l) != rho for l in leurres)
    # distance pour la norme subordonnée euclidienne (preuve ajoutée à la relecture) : ||A^{-1}|| = 1/a0
    # (norme subordonnée euclidienne d'une matrice diagonale = plus grand |coefficient|), et A - a0 E11 est singulière,
    # avec ||a0 E11|| = a0 : la distance vaut a0.
    assert min(a0, b0) == a0 and max(F(1, a0), F(1, b0)) == F(1, a0)
    assert (diag(a0, b0) - a0 * E11).det() == 0
    assert max(abs(x) for x in (a0 * E11).singular_values()) == a0

# --- al-l3-dunford-densite-ex3 : AB = E12, BA = 0.
assert E11 * E12 == E12 and E12 * E11 == zeros(2) and chi(E12) == chi(zeros(2)) == X ** 2

# --- al-l3-dunford-densite-ex4 (et v1, v2) : distance de C aux diagonalisables réelles (norme max).
for (c12, c21, rmax, rho, leurres) in [(-1, 1, 1, F(1, 2), [1, F(1, 4)]), (-1, 4, 1, F(4, 5), [1, 2, F(1, 2)]), (-2, 2, 2, F(1), [2, F(1, 2)])]:
    bb, cc = -c12, c21
    bound = 4 * r ** 2 - 4 * (bb - r) * (cc - r)
    assert bound.subs(r, rho) == 0 and rho < rmax
    # majorant de delta sur la boule fermée de rayon r (échantillonnage exact), pour quelques r < rho et r > rho
    for rr in [rho * F(9, 10), rho, (rho + rmax) / 2]:
        for _ in range(300):
            aa, dd, be, ga = [F(random.randint(-1000, 1000), 1000) * rr for _ in range(4)]
            delta = (aa - dd) ** 2 + 4 * (c12 + be) * (c21 + ga)
            assert delta <= bound.subs(r, rr)
    # B_r diagonalisable pour r dans ]rho, rmax[
    rr = (rho + rmax) / 2
    Br = Matrix([[rr, c12 + rr], [c21 - rr, -rr]])
    assert max(abs(x) for x in (Br - Matrix([[0, c12], [c21, 0]]))) == rr
    assert (Br.trace()) ** 2 - 4 * Br.det() > 0
    # en r = rho, la seule matrice extrémale (au signe près) est nilpotente non nulle
    Bc = Matrix([[rho, c12 + rho], [c21 - rho, -rho]])
    assert Bc.trace() == 0 and Bc.det() == 0 and Bc != zeros(2)
    assert all(F(l) != rho for l in leurres)
assert expand(4 * r ** 2 - 4 * (1 - r) ** 2 - 4 * (2 * r - 1)) == 0
assert expand(4 * r ** 2 - 4 * (1 - r) * (4 - r) - 4 * (5 * r - 4)) == 0
assert expand(4 * r ** 2 - 4 * (2 - r) ** 2 - 16 * (r - 1)) == 0

# --- al-l3-dunford-densite-ex5 : polynômes caractéristiques des matrices proposées.
assert expand(chi(J) - (X - 1) ** 2) == 0
assert expand(chi(Matrix([[1, -1], [1, 1]])) - (X ** 2 - 2 * X + 2)) == 0
assert expand(chi(Matrix([[0, 1], [-2, 0]])) - (X ** 2 + 2)) == 0
# |Q(z)| >= |Im z|^m pour Q scindé sur R (échantillonnage)
Q = (X - 1) * (X + 2) * (X - F(1, 3))
from sympy import Abs, im, I as II
for _ in range(50):
    z = random.uniform(-3, 3) + random.uniform(-3, 3) * 1j
    val = abs(complex(Q.subs(X, z).evalf()))
    assert val >= abs(z.imag) ** 3 - 1e-9

# ===========================================================================
# al-l3-dunford-systemes (carte)
# ===========================================================================
A = Matrix([[1, 1, 1], [0, 1, 1], [0, 0, 2]])
p2 = (A - I3) ** 2
p1 = I3 - p2
assert p1 == Matrix([[1, 0, -2], [0, 1, -1], [0, 0, 0]]) and (A - I3) * p1 == Matrix([[0, 1, -1], [0, 0, 0], [0, 0, 0]])
EtA = exp(t) * (p1 + t * (A - I3) * p1) + exp(2 * t) * p2
Ref = Matrix([[exp(t), t * exp(t), 2 * exp(2 * t) - (t + 2) * exp(t)], [0, exp(t), exp(2 * t) - exp(t)], [0, 0, exp(2 * t)]])
assert eqm(EtA, Ref) and eqm((t * A).exp(), Ref)
# Contre-exemples : E12 (solution (t, 1)) ; 0_2 ; -2 I.
assert eqm((t * E12).exp() * Matrix([0, 1]), Matrix([t, 1]))
# Trou 1, leurre « partie imaginaire nulle » : rotation, solutions bornées ; trou 2, leurre « sous-espace propre de dimension 1 » : J.
assert eqm((t * R).exp(), Matrix([[cos(t), -sin(t)], [sin(t), cos(t)]]))
# Affirmation 4 : A antisymétrique => d/dt |X|^2 = 0 (exemple de taille 3).
Aas = Matrix([[0, 1, -2], [-1, 0, 3], [2, -3, 0]])
Xv = Matrix(symbols("y1:4"))
assert expand((Xv.T * Aas * Xv)[0]) == 0

# --- al-l3-dunford-systemes-ex1 (et v1, v2)
for A, vps, cas in [([[-1, 2, 0], [0, -3, 1], [0, 0, -2]], (-1, -3, -2), 0), ([[0, 1, 0], [0, -1, 2], [0, 0, -2]], (0, -1, -2), 1),
                    ([[0, 1, 0], [0, 0, 3], [0, 0, -1]], (0, 0, -1), 2)]:
    A = Matrix(A)
    assert multiset_vp(A) == vp(*vps)
    if cas == 0:
        assert A * Matrix([1, 0, 0]) == -Matrix([1, 0, 0])
    if cas == 1:
        assert A * Matrix([1, 0, 0]) == zeros(3, 1) and A * Matrix([1, -1, 0]) == -Matrix([1, -1, 0])
    if cas == 2:
        assert A.rank() == 2 and eqm((t * A).exp() * Matrix([0, 1, 0]), Matrix([t, 1, 0]))

# --- al-l3-dunford-systemes-ex2 (et v1, v2) : solution issue de (1,0,0) ; leurres.
for A, lam, mu, zm, y0, ay0, sol, leurres in [
    ([[2, 0, 0], [-1, 0, 1], [0, -1, 2]], 1, 2, (1, 0, 1), (0, 0, -1), (0, -1, -1),
     [exp(2 * t), -t * exp(t), exp(2 * t) - (1 + t) * exp(t)],
     [[exp(2 * t), 0, exp(2 * t) - exp(t)], [exp(2 * t), t * exp(t), exp(2 * t) - (1 - t) * exp(t)], [exp(2 * t), -t * exp(t), exp(2 * t) - exp(t)]]),
    ([[0, 0, 0], [1, 0, -1], [2, 1, -2]], -1, 0, (1, 0, 1), (0, 0, -1), (0, 1, 1),
     [1, t * exp(-t), 1 - (1 - t) * exp(-t)],
     [[1, 0, 1 - exp(-t)], [1, -t * exp(-t), 1 - (1 + t) * exp(-t)], [1, t * exp(-t), 1 - exp(-t)]]),
    ([[-1, 0, 0], [-2, 0, 1], [-1, 0, 0]], 0, -1, (1, 1, 1), (0, -1, -1), (0, -1, 0),
     [exp(-t), exp(-t) - 1 - t, exp(-t) - 1],
     [[exp(-t), exp(-t) - 1, exp(-t) - 1], [exp(-t), exp(-t) - 1 + t, exp(-t) - 1], [exp(-t), exp(-t) - 1 - t, exp(-t) - 1 - t]]),
]:
    A = Matrix(A)
    assert expand(chi(A) - (X - lam) ** 2 * (X - mu)) == 0
    assert A * Matrix(zm) == mu * Matrix(zm) and ((A - lam * I3) ** 2) * Matrix(y0) == zeros(3, 1)
    assert ((A - lam * I3) ** 2).nullspace() and all(v[0] == 0 for v in ((A - lam * I3) ** 2).nullspace())
    assert Matrix(zm) + Matrix(y0) == Matrix([1, 0, 0]) and (A - lam * I3) * Matrix(y0) == Matrix(ay0)
    S = Matrix(sol)
    assert eqm((t * A).exp() * Matrix([1, 0, 0]), S)
    assert eqm(S.diff(t), A * S) and S.subs(t, 0) == Matrix([1, 0, 0])
    for L in leurres:
        L = Matrix(L)
        assert L.subs(t, 0) == Matrix([1, 0, 0]) and not eqm(L, S)

# --- al-l3-dunford-systemes-ex3 (et v1, v2) : résonance.
x = Function("x")
for w, y_expr, x_expr, leurres in [
    (1, sin(t), (sin(t) - t * cos(t)) / 2, [t * sin(t) / 2, (sin(t) + t * cos(t)) / 2, -t * cos(t) / 2]),
    (1, cos(t), t * sin(t) / 2, [(sin(t) - t * cos(t)) / 2, t * sin(t), t * cos(t) / 2]),
    (2, sin(2 * t), (sin(2 * t) - 2 * t * cos(2 * t)) / 8, [(sin(2 * t) - 2 * t * cos(2 * t)) / 4, t * sin(2 * t) / 4, -t * cos(2 * t) / 4]),
]:
    A = Matrix([[0, 1, 0, 0], [-w ** 2, 0, 1, 0], [0, 0, 0, 1], [0, 0, -w ** 2, 0]])
    assert expand(chi(A) - (X ** 2 + w ** 2) ** 2) == 0
    assert (A - w * i_ * eye(4)).rank() == 3 and not A.is_diagonalizable()
    assert simplify(diff(y_expr, t, 2) + w ** 2 * y_expr) == 0
    assert simplify(diff(x_expr, t, 2) + w ** 2 * x_expr - y_expr) == 0
    assert x_expr.subs(t, 0) == 0 and diff(x_expr, t).subs(t, 0) == 0
    for l in leurres:
        ok = (simplify(diff(l, t, 2) + w ** 2 * l - y_expr) == 0 and l.subs(t, 0) == 0 and diff(l, t).subs(t, 0) == 0)
        assert not ok
# conditions initiales des trois cas
assert sin(0) == 0 and diff(sin(t), t).subs(t, 0) == 1 and cos(0) == 1 and diff(sin(2 * t), t).subs(t, 0) == 2

# --- al-l3-dunford-systemes-ex4 (et v1, v2) : factorisations et bornitude.
def bornees(P):
    """Toutes les solutions de l'équation de polynôme caractéristique P sont bornées sur [0, +inf[."""
    from sympy import roots
    rts = roots(Poly(P, X))
    assert sum(rts.values()) == Poly(P, X).degree()
    for z, mlt in rts.items():
        rz = re(z)
        if rz > 0 or (rz == 0 and mlt >= 2):
            return False
    return True


for bonne, leurres_ in [
    (X ** 3 + X, [X ** 4 + 2 * X ** 2 + 1, X ** 3, X ** 2 - X]),
    (X ** 4 + 5 * X ** 2 + 4, [X ** 4 + 2 * X ** 2 + 1, X ** 4 - 1, X ** 3 + X ** 2]),
    (X ** 3 + 2 * X ** 2 + X, [X ** 4 + 4 * X ** 2, X ** 3 - X ** 2 + X - 1, X ** 3 - X]),
]:
    assert bornees(bonne) and not any(bornees(l) for l in leurres_)
assert expand(X * (X ** 2 + 1) - (X ** 3 + X)) == 0 and expand((X ** 2 + 1) * (X ** 2 + 4) - (X ** 4 + 5 * X ** 2 + 4)) == 0
assert expand((X - 1) * (X + 1) * (X ** 2 + 1) - (X ** 4 - 1)) == 0 and expand(X * (X + 1) ** 2 - (X ** 3 + 2 * X ** 2 + X)) == 0
assert expand((X - 1) * (X ** 2 + 1) - (X ** 3 - X ** 2 + X - 1)) == 0
# t e^{lambda t} solution si lambda racine double : (t cos t) pour x'''' + 2x'' + x = 0
assert simplify(diff(t * cos(t), t, 4) + 2 * diff(t * cos(t), t, 2) + t * cos(t)) == 0
# dérivées de t e^{lambda t}
lam_ = symbols("lambda")
for kk in range(1, 6):
    assert simplify(diff(t * exp(lam_ * t), t, kk) - (t * lam_ ** kk + kk * lam_ ** (kk - 1)) * exp(lam_ * t)) == 0

# --- al-l3-dunford-systemes-ex5 (et v1, v2) : plus petite période commune.
for P, sols, T0, leurres in [
    (X ** 4 + 10 * X ** 2 + 9, [cos(t), sin(t), cos(3 * t), sin(3 * t)], 2 * pi, [2 * pi / 3, pi, 6 * pi]),
    (X ** 4 + 20 * X ** 2 + 64, [cos(2 * t), sin(2 * t), cos(4 * t), sin(4 * t)], pi, [pi / 2, 2 * pi, pi / 4]),
    (X ** 3 + 4 * X, [1 + 0 * t, cos(2 * t), sin(2 * t)], pi, [2 * pi, pi / 2, 4 * pi]),
]:
    from sympy import roots
    rts = roots(Poly(P, X))
    assert all(m == 1 for m in rts.values()) and all(re(z) == 0 for z in rts)
    for s in sols:
        # s est solution de l'équation
        n = Poly(P, X).degree()
        coeffs = Poly(P, X).all_coeffs()
        assert simplify(sum(c * diff(s, t, n - j) for j, c in enumerate(coeffs))) == 0
        assert simplify(s.subs(t, t + T0) - s) == 0
    # aucune période strictement plus petite parmi les leurres : un leurre < T0 échoue pour une solution
    for l in leurres:
        if l < T0:
            assert any(simplify(s.subs(t, t + l) - s) != 0 for s in sols)
    # T0 = plus petit T > 0 tel que e^{T lambda} = 1 pour toute racine lambda
    for z in rts:
        assert simplify(exp(T0 * z) - 1) == 0
assert expand((X ** 2 + 1) * (X ** 2 + 9) - (X ** 4 + 10 * X ** 2 + 9)) == 0
assert expand((X ** 2 + 4) * (X ** 2 + 16) - (X ** 4 + 20 * X ** 2 + 64)) == 0

# ===========================================================================
# QCM al-l3-dunford-q1 … q8
# ===========================================================================
# q1 : commutant de diag(1,1,2) de dimension 5.
Xs = symbols("z0:9")
Xm = Matrix(3, 3, Xs)
from sympy import linsolve
solq1 = list(linsolve(list(Xm * diag(1, 1, 2) - diag(1, 1, 2) * Xm), Xs))[0]
assert len(solq1.free_symbols) == 5
# q2 : I2 et J.
assert I2 * J == J * I2 and not J.is_diagonalizable()
# q3 : [[1,1],[0,2]] diagonalisable, Dunford (A, 0).
assert dunford_ok(A2, A2, zeros(2)) and A2 != diag(1, 2) and A2 != I2 and A2 != Rational(3, 2) * I2
# q4 : D polynôme en A (exemple de la carte) : commute avec tout ce qui commute avec A.
A = Matrix([[1, 1, 1], [0, 1, 1], [0, 0, 2]])
assert evalp(1 + (X - 1) ** 2, A) == Matrix([[1, 0, 2], [0, 1, 1], [0, 0, 2]])
# q5 : e^{pi K} = -I2 ; leurres.
assert eqm((pi * K).exp(), -I2) and not eqm((pi * K).exp(), I2)
# q6 : epsilon E12 non diagonalisable.
assert not (Rational(1, 10) * E12).is_diagonalizable()
# q7 : seule [[1,-1],[1,1]] n'est pas trigonalisable sur R.
for M_, trig in [([[1, -1], [1, 1]], False), ([[0, 1], [0, 0]], True), ([[1, 1], [0, 1]], True), ([[2, 1], [1, 2]], True)]:
    M_ = Matrix(M_)
    disc = M_.trace() ** 2 - 4 * M_.det()
    assert (disc >= 0) == trig
# q8 : E12, solution (t, 1).
assert eqm((t * E12).exp() * Matrix([0, 1]), Matrix([t, 1]))

print("al_l3_dunford : toutes les vérifications passent.")
