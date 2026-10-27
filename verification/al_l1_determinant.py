"""Vérifications calculatoires du chapitre « al-l1-determinant » (thème algebre-lineaire, niveau l1).

Cartes content/cartes/algebre-lineaire/l1-al-l1-determinant.json et exercices
content/exercices/algebre-lineaire/l1-al-l1-determinant.json. Chaque assert cite l'identifiant de la carte ou de
l'exercice vérifié (variante comprise : « v1 », « v2 »). Les leurres des questions à choix et des textes à trous sont
vérifiés FAUX (sur un exemple, ou pour au moins une valeur du paramètre).
Lancer : python3 verification/al_l1_determinant.py (doit finir sans erreur).
"""

from itertools import permutations

from sympy import (Matrix, Poly, Rational, binomial, diff, eye, exp, expand, factor, factorial, ff, simplify, sqrt,
                   symbols, zeros)

a, b, c, m, t, x, X, lam = symbols("a b c m t x X lambda")


def colonnes(*vs):
    """Matrice dont les colonnes sont les vecteurs donnés."""
    return Matrix([list(v) for v in vs]).T


def sarrus(M):
    return (M[0, 0] * M[1, 1] * M[2, 2] + M[0, 1] * M[1, 2] * M[2, 0] + M[0, 2] * M[1, 0] * M[2, 1]
            - M[0, 2] * M[1, 1] * M[2, 0] - M[0, 0] * M[1, 2] * M[2, 1] - M[0, 1] * M[1, 0] * M[2, 2])


def signature(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def leibniz(M):
    n = M.rows
    total = 0
    for p in permutations(range(n)):
        terme = signature(p)
        for j in range(n):
            terme *= M[p[j], j]
        total += terme
    return expand(total)


def mineur(M, i, j):
    """Δ_{ij} (indices à partir de 1)."""
    return M.minor_submatrix(i - 1, j - 1).det()


def vdm(*xs):
    return Matrix([[xi ** j for j in range(len(xs))] for xi in xs]).det()


def produit_vdm(*xs):
    p = 1
    for i in range(len(xs)):
        for j in range(i + 1, len(xs)):
            p *= xs[j] - xs[i]
    return p


def coords(P):
    P = Poly(P, X)
    return [P.coeff_monomial(X ** k) for k in range(3)]


# ===========================================================================
# Cartes
# ===========================================================================
# al-l1-determinant-famille : la formule explicite (Leibniz) coïncide avec le déterminant, Sarrus pour n = 3.
for n in (2, 3, 4):
    S = Matrix(n, n, symbols(f"s0:{n * n}"))
    assert expand(leibniz(S) - S.det()) == 0  # al-l1-determinant-famille (formule explicite)
S3 = Matrix(3, 3, symbols("q0:9"))
assert expand(sarrus(S3) - S3.det()) == 0  # al-l1-determinant-famille (règle de Sarrus)
T = Matrix([[a, 1, 2], [0, b, 3], [0, 0, c]])
assert T.det() == a * b * c and T.T.det() == a * b * c  # al-l1-determinant-famille (triangulaire)

# al-l1-determinant-famille, exemple
x1, x2, x3 = (1, 0, 2), (2, 1, 0), (0, 3, 1)
A = colonnes(x1, x2, x3)
assert A == Matrix([[1, 2, 0], [0, 1, 3], [2, 0, 1]])
assert sarrus(A) == 13 and A.det() == 13  # al-l1-determinant-famille (exemple : 13)
assert colonnes(x2, x1, x3).det() == -13  # al-l1-determinant-famille (exemple : échange)
assert colonnes(x1, tuple(x2[i] + 5 * x1[i] for i in range(3)), x3).det() == 13  # al-l1-determinant-famille (exemple)
Bp = colonnes((2, 0, 0), (0, 1, 0), (0, 0, 1))  # B' = (2e1, e2, e3)
assert Bp.det() == 2 and (Bp.inv() * A).det() == Rational(13, 2)  # al-l1-determinant-famille (exemple, det_B'(x) = 13/2)

# al-l1-determinant-famille, contre-exemples (R^2, base canonique)
y1, y2, z1, z2, al, be = symbols("y1 y2 z1 z2 alpha beta")
f1 = lambda u, v: u[0] * v[1]
assert f1((1, 0), (0, 1)) == 1 and f1((0, 1), (1, 0)) == 0 != colonnes((0, 1), (1, 0)).det()  # al-l1-determinant-famille (sans alternance)
f2 = lambda u, v: u[0] * v[1] + u[1] * v[0]
assert f2((1, 0), (0, 1)) == 1 and f2((0, 1), (1, 0)) == 1  # al-l1-determinant-famille (symétrique, f ≠ det)
assert expand(f2((y1, y2), (z1, z2)) - f2((z1, z2), (y1, y2))) == 0  # al-l1-determinant-famille (f2 symétrique)
f3 = lambda u, v: (u[0] * v[1] - u[1] * v[0]) * (1 + (v[1] - 1) ** 2)
assert expand(f3((y1, y2), (y1, y2))) == 0  # al-l1-determinant-famille (f3 alternée)
assert expand(f3((al * y1 + be * z1, al * y2 + be * z2), (x, t))
              - al * f3((y1, y2), (x, t)) - be * f3((z1, z2), (x, t))) == 0  # al-l1-determinant-famille (f3 linéaire en x)
assert f3((1, 0), (0, 1)) == 1 and f3((1, 0), (0, 2)) == 4 and colonnes((1, 0), (0, 2)).det() == 2  # al-l1-determinant-famille
assert colonnes((1, 0), (2, 0)).det() == 0  # al-l1-determinant-famille (famille liée : aucune forme alternée n'y vaut 1)
# pièges : sans alternance / symétrique / linéaire en la première variable seulement : deux formes distinctes valent 1 en B
assert f1((1, 0), (0, 1)) == colonnes((1, 0), (0, 1)).det() and f1((0, 1), (1, 0)) != colonnes((0, 1), (1, 0)).det()

# al-l1-determinant-famille, affirmations
assert colonnes((2, 0), (0, 2)).det() == 4 != 2  # al-l1-determinant-famille (det(λx) = λ^n det x)
Bp2 = colonnes((2, 0), (0, 1))
assert Bp2.det() == 2 and (Bp2.inv() * eye(2)).det() == Rational(1, 2)  # al-l1-determinant-famille (dépend de la base)
assert colonnes((1, 0), (0, 1)).det() == 1 and colonnes((1, 0), (0, 0)).det() + colonnes((0, 0), (0, 1)).det() == 0  # al-l1-determinant-famille
# al-l1-determinant-famille, trous : leurres faux
assert colonnes(x2, x1, x3).det() == -A.det() and -A.det() not in (A.det(), 2 * A.det(), 0)  # trou 2
assert colonnes(x1, tuple(x2[i] + 5 * x1[i] for i in range(3)), x3).det() not in (-13, 0, 3 * 13)  # trou 3
g1 = lambda u, v: u[0] * v[1]
g2 = lambda u, v: u[0] * v[1] + u[0] * v[0]  # bilinéaire, vaut 1 en (e1, e2), différente de g1
assert g1((1, 0), (0, 1)) == g2((1, 0), (0, 1)) == 1 and g1((1, 0), (1, 0)) != g2((1, 0), (1, 0))  # trou 1 (« n-linéaire »)
h1 = lambda u, v: u[0] * v[1] + u[1] * v[0]
h2 = lambda u, v: u[0] * v[1] + u[1] * v[0] + u[0] * v[0]  # bilinéaires symétriques valant 1 en (e1, e2)
assert h1((1, 0), (0, 1)) == h2((1, 0), (0, 1)) == 1 and h1((1, 0), (1, 0)) != h2((1, 0), (1, 0))  # trou 1 (« symétrique »)
l1_ = lambda u, v: u[0]
l2_ = lambda u, v: v[1]  # formes linéaires sur E^2 valant 1 en (e1, e2)
assert l1_((1, 0), (0, 1)) == l2_((1, 0), (0, 1)) == 1 and l1_((1, 0), (0, 0)) != l2_((1, 0), (0, 0))  # trou 1 (« linéaire »)

# al-l1-determinant-produit
S = Matrix(3, 3, symbols("p0:9"))
R = Matrix(3, 3, symbols("r0:9"))
assert expand((S * R).det() - S.det() * R.det()) == 0  # al-l1-determinant-produit (formule du produit, n = 3)
A2 = Matrix([[1, 2], [3, 4]])
B2 = Matrix([[0, 1], [1, 1]])
assert A2.det() == -2 and B2.det() == -1  # al-l1-determinant-produit (exemple)
assert A2 * B2 == Matrix([[2, 3], [4, 7]]) and (A2 * B2).det() == 2  # al-l1-determinant-produit (exemple)
assert B2 * A2 == Matrix([[3, 4], [4, 6]]) and (B2 * A2).det() == 2 and A2 * B2 != B2 * A2  # al-l1-determinant-produit (exemple)
assert (A2 ** 5).det() == -32 and (3 * A2).det() == -18  # al-l1-determinant-produit (exemple)
L = Matrix([[1, 0]])
C = Matrix([[1], [0]])
assert (L * C).det() == 1 and (C * L).det() == 0  # al-l1-determinant-produit (contre-exemple : non carrées)
assert (eye(2) + eye(2)).det() == 4 and eye(2).det() + eye(2).det() == 2  # al-l1-determinant-produit (pas d'additivité)
assert (2 * eye(3)).det() == 8  # al-l1-determinant-produit (pas de linéarité)
N = Matrix([[0, 1], [0, 0]])
assert N * N == zeros(2, 2) and N != zeros(2, 2)  # al-l1-determinant-produit (affirmation AB = 0)
assert Matrix([[1, 1], [0, 1]]).det() == 1 == eye(2).det()  # al-l1-determinant-produit (même det, non semblables)
assert (-eye(3)).det() == -1  # al-l1-determinant-produit (A^2 = -I_3 impossible : (det A)^2 = -1)
# trous : leurres faux
d1, d2 = A2.det(), B2.det()
assert (A2 * B2).det() not in (d1 + d2, 2 * d1 * d2, (d1 * d2) ** 2)  # al-l1-determinant-produit (trou 1, n = 2)
assert (-1 * A2).det() == (-1) ** 2 * d1 and (-1 * A2).det() not in (-d1, 2 * (-1) * d1)  # trou 2 (λ = -1)
assert (Rational(-1) * eye(3)).det() != abs(-1) ** 3  # trou 2 (|λ|^n faux pour λ = -1, n = 3)
P = Matrix([[2, 1], [1, 2]])
D = P.inv() * A2 * P
assert P.det() == 3 and D.det() == d1 and D.det() not in (-d1, 1 / d1, P.det() * d1)  # al-l1-determinant-produit (trou 3)

# al-l1-determinant-inversibilite
Am = Matrix([[1, 1, 1], [1, m, 1], [1, 1, m]])
assert factor(Am.det()) == (m - 1) ** 2  # al-l1-determinant-inversibilite (exemple)
assert Am.subs(m, 2).inv().det() == 1  # al-l1-determinant-inversibilite (exemple, m = 2)
Ac = Matrix.hstack(Am[:, 0], Am[:, 1] - Am[:, 0], Am[:, 2] - Am[:, 0])
assert Ac == Matrix([[1, 0, 0], [1, m - 1, 0], [1, 0, m - 1]]) and expand(Ac.det() - Am.det()) == 0  # inversibilite (exemple, colonnes)
assert Matrix([[1, 2], [2, 4]]).det() == 0  # al-l1-determinant-inversibilite (coefficients non nuls)
assert colonnes((1, 0, 1), (0, 1, 1), (1, 1, 2)).det() == 0  # al-l1-determinant-inversibilite (deux à deux non colinéaires)
for u, v in (((1, 0, 1), (0, 1, 1)), ((1, 0, 1), (1, 1, 2)), ((0, 1, 1), (1, 1, 2))):
    assert colonnes(u, v).rank() == 2  # al-l1-determinant-inversibilite (non colinéaires deux à deux)
assert Matrix([[1, 1], [1, 1]]).det() == 0  # al-l1-determinant-inversibilite (diagonale non nulle)
J = Matrix([[0, 1], [1, 0]])
assert J * J == eye(2) and J.det() == -1  # al-l1-determinant-inversibilite (signe sans importance ; leurre « det A > 0 »)
assert (eye(2) + (-eye(2))) == zeros(2, 2)  # al-l1-determinant-inversibilite (affirmation A + B)
assert (2 * eye(2)).det() == 4  # al-l1-determinant-inversibilite (leurre « det A = 1 » : 2 I_2 inversible)
A4 = 2 * eye(2)
assert A4.inv().det() == Rational(1, 4) and Rational(1, 4) not in (A4.det(), -A4.det(), A4.det() ** (2 - 1))  # trou 2
assert colonnes((1, 0), (0, 2)).det() == 2  # al-l1-determinant-inversibilite (trou 3 : base de det 2 ≠ 1)
assert colonnes((1, 0), (2, 0)).det() == 0  # al-l1-determinant-inversibilite (trou 3 : vecteurs non nuls, famille liée)

# al-l1-determinant-developpement
for n in (2, 3, 4):
    S = Matrix(n, n, symbols(f"d0:{n * n}"))
    assert expand(S.T.det() - S.det()) == 0  # al-l1-determinant-developpement (transposée)
    for j in range(1, n + 1):
        assert expand(sum((-1) ** (i + j) * S[i - 1, j - 1] * mineur(S, i, j) for i in range(1, n + 1)) - S.det()) == 0  # colonne j
    for i in range(1, n + 1):
        assert expand(sum((-1) ** (i + j) * S[i - 1, j - 1] * mineur(S, i, j) for j in range(1, n + 1)) - S.det()) == 0  # ligne i
    Com = Matrix(n, n, lambda i, j: (-1) ** (i + j) * mineur(S, i + 1, j + 1))
    assert expand(S * Com.T - S.det() * eye(n)) == zeros(n, n)  # al-l1-determinant-developpement (comatrice)
    assert expand(Com.T * S - S.det() * eye(n)) == zeros(n, n)  # al-l1-determinant-developpement (comatrice)
A = Matrix([[2, 0, 1, 3], [1, 0, 0, 2], [0, 1, 4, 1], [3, 0, 2, 0]])
assert mineur(A, 3, 2) == 4 and A.det() == -4  # al-l1-determinant-developpement (exemple)
assert Matrix([[2, 1, 3], [1, 0, 2], [3, 2, 0]]).det() == 4  # al-l1-determinant-developpement (exemple, Δ32)
assert Matrix([[1, 3], [2, 0]]).det() == -6 and Matrix([[2, 1], [3, 2]]).det() == 1  # al-l1-determinant-developpement (exemple)
Bs = Matrix([[a, b], [c, symbols("d")]])
assert Bs.adjugate() == Matrix([[symbols("d"), -b], [-c, a]])  # al-l1-determinant-developpement (exemple, t(Com B))
# contre-exemples avec A2 = [[1, 2], [3, 4]]
assert 1 * 4 + 3 * 2 == 10 != A2.det() and -1 * 4 + 3 * 2 == 2 != A2.det()  # al-l1-determinant-developpement (signes)
Com2 = Matrix([[4, -3], [-2, 1]])
assert Com2 == A2.adjugate().T  # al-l1-determinant-developpement (Com(A))
fausse = Com2 / -2
assert fausse == Matrix([[-2, Rational(3, 2)], [1, Rational(-1, 2)]]) and (A2 * fausse)[0, 0] == 0  # comatrice non transposée
assert A2.inv() == Matrix([[-2, 1], [Rational(3, 2), Rational(-1, 2)]]) == Com2.T / -2  # al-l1-determinant-developpement
Bc = Matrix([[1, 2, 0], [0, 1, 3], [4, 0, 1]])
assert Bc.det() == 25  # al-l1-determinant-developpement (contre-exemple, indices échangés)
assert [mineur(Bc, 1, i) for i in (1, 2, 3)] == [1, -12, -4]
assert sum((-1) ** (i + 1) * Bc[i - 1, 0] * mineur(Bc, 1, i) for i in (1, 2, 3)) == -15  # piège Δ_ji : -15 ≠ 25
U = Matrix([[1, 1], [1, 1]])
assert U.adjugate().T == Matrix([[1, -1], [-1, 1]]) and U.det() == 0  # al-l1-determinant-developpement (cofacteurs non nuls)
# trous : leurres faux, sur A2 et Bc
assert sum((-1) ** (i * 1) * A2[i - 1, 0] * mineur(A2, i, 1) for i in (1, 2)) != A2.det()  # trou 1 : (-1)^{ij}
assert sum(Bc[i - 1, 1] * mineur(Bc, 1, 2) * (-1) ** (i + 2) for i in (1, 2, 3)) != 25  # trou 2 : première ligne
assert sum(Bc[i - 1, 0] * mineur(Bc, i, 3) * (-1) ** (i + 1) for i in (1, 2, 3)) != 25  # trou 2 : dernière colonne
for leurre in (Com2 / A2.det(), A2.det() * Com2.T, Com2.T):
    assert A2 * leurre != eye(2)  # al-l1-determinant-developpement (trou 3, leurres faux)
assert Matrix([[2, 0], [0, 1]]).inv() == Matrix([[Rational(1, 2), 0], [0, 1]])  # al-l1-determinant-developpement (affirmation)
p_, q_, r_ = symbols("p q r")
Anti = Matrix([[0, p_, q_], [-p_, 0, r_], [-q_, -r_, 0]])
assert expand(Anti.det()) == 0  # al-l1-determinant-developpement (antisymétrique d'ordre 3)

# al-l1-determinant-vandermonde
for n in (2, 3, 4, 5):
    xs = symbols(f"v1:{n + 1}")
    assert expand(vdm(*xs) - produit_vdm(*xs)) == 0  # al-l1-determinant-vandermonde (formule)
    assert expand(Matrix([[xi ** j for j in range(n)] for xi in xs]).T.det() - produit_vdm(*xs)) == 0  # transposée
    assert expand(Matrix([[xi ** (j + 1) for j in range(n)] for xi in xs]).det()
                  - produit_vdm(*xs) * __import__("sympy").prod(xs)) == 0  # puissances 1..n : x1⋯xn V
    pi = 1
    for i in range(n):
        for j in range(i + 1, n):
            pi *= xs[i] - xs[j]
    assert expand(pi - (-1) ** (n * (n - 1) // 2) * produit_vdm(*xs)) == 0  # signe (-1)^{n(n-1)/2}
assert vdm(1, 2, 3) == 2 and vdm(0, 1, -1) == 2  # al-l1-determinant-vandermonde (exemple)
Pi = X ** 2 - 2 * X + 3
assert [Pi.subs(X, k) for k in (1, 2, 3)] == [2, 3, 6]  # al-l1-determinant-vandermonde (exemple, interpolation)
assert vdm(1, 2, 1) == 0 and vdm(3, 3) == 0  # al-l1-determinant-vandermonde (contre-exemples)
assert Matrix([[2, 4], [3, 9]]).det() == 6 and vdm(2, 3) == 1  # al-l1-determinant-vandermonde (puissances 1..n)
assert vdm(1, 2) == 1 and vdm(2, 1) == -1  # al-l1-determinant-vandermonde (affirmation : permutation)
assert vdm(0, 1) != 0  # interpolation (0, 1) -> (1, 1) : c0 = 1, c1 = 0, polynôme constant
assert Matrix([[1, 0], [1, 1]]).solve(Matrix([1, 1])) == Matrix([1, 0])  # al-l1-determinant-vandermonde (affirmation)
# trous : leurres faux
assert (1 - 2) != vdm(1, 2)  # trou 1 : produit des x_i - x_j (n = 2)
assert (2 - 1) + (3 - 1) + (3 - 2) != vdm(1, 2, 3)  # trou 1 : somme
assert 1 * 2 * 3 * vdm(1, 2, 3) != vdm(1, 2, 3)  # trou 1 : x1 x2 x3 V
assert vdm(1, 1) == 0 and vdm(2, 1) != 0 and vdm(1, 2, 1) == 0  # trou 2 : non nuls, ordre croissant, consécutifs

# ===========================================================================
# Exercices
# ===========================================================================
# al-l1-determinant-famille-ex1 (et variantes) : valeur et leurres (−d, autres erreurs de Sarrus)
for (u, v, w), d, leurres in [(((1, 2, 0), (3, -1, 4), (2, 0, 1)), 9, (-9, 21, 15)),
                              (((2, 1, -1), (0, 3, 1), (1, -2, 2)), 20, (-20, 6, 13)),
                              (((1, -1, 3), (2, 0, 1), (-1, 4, 2)), 25, (-25, 17, 21))]:
    M = colonnes(u, v, w)
    assert M.det() == d and sarrus(M) == d and d not in leurres  # al-l1-determinant-famille-ex1 (v1, v2)
M = colonnes((1, 2, 0), (3, -1, 4), (2, 0, 1))
assert (M[0, 0] * M[1, 1] * M[2, 2], M[0, 2] * M[1, 0] * M[2, 1], M[0, 1] * M[1, 0] * M[2, 2]) == (-1, 16, 6)  # ex1
M = colonnes((2, 1, -1), (0, 3, 1), (1, -2, 2))
assert (M[0, 0] * M[1, 1] * M[2, 2], M[0, 2] * M[1, 0] * M[2, 1], M[0, 2] * M[1, 1] * M[2, 0], M[0, 0] * M[1, 2] * M[2, 1]) == (12, 1, -3, -4)  # ex1 v1
M = colonnes((1, -1, 3), (2, 0, 1), (-1, 4, 2))
assert (M[0, 1] * M[1, 2] * M[2, 0], M[0, 2] * M[1, 0] * M[2, 1], M[0, 0] * M[1, 2] * M[2, 1], M[0, 1] * M[1, 0] * M[2, 2]) == (24, 1, 4, -4)  # ex1 v2

# al-l1-determinant-famille-ex2 : on prend u, v, w génériques (coordonnées symboliques)
U3 = Matrix(3, 3, symbols("w0:9"))
u, v, w = U3[:, 0], U3[:, 1], U3[:, 2]
dd = U3.det()
assert expand(Matrix.hstack(u - 2 * v, 3 * v, w + u + v).det() - 3 * dd) == 0  # famille-ex2 : 3 det, soit 12
assert 3 * 4 == 12 and 12 not in (4, -12, 24)
assert expand(Matrix.hstack(2 * u + w, v - u, 5 * w).det() - 10 * dd) == 0  # famille-ex2 v1 : 10 det, soit -20
assert 10 * (-2) == -20 and -20 not in (20, -10, 10)
assert expand(Matrix.hstack(w, u + v, -2 * v).det() + 2 * dd) == 0  # famille-ex2 v2 : -2 det, soit -6
assert expand(Matrix.hstack(w, u, v).det() - dd) == 0 and -2 * 3 == -6 and -6 not in (6, -3, 3)

# al-l1-determinant-famille-ex3 (et variantes)
for Bp_, F, dBp, dF, rep, leurres in [([1 + X, 1 - X, X ** 2 + 2], [X ** 2, X, 1], -2, -1, Rational(1, 2), (2, Rational(-1, 2), -2)),
                                      ([1 + X, X + X ** 2, 1 + X ** 2], [1 + X ** 2, X, 1], 2, -1, Rational(-1, 2), (-2, Rational(1, 2), 2)),
                                      ([1 + X, 2 * X + X ** 2, 1 + X ** 2], [X ** 2, X, 1], 3, -1, Rational(-1, 3), (-3, Rational(1, 3), 3))]:
    MB = colonnes(*[coords(p) for p in Bp_])
    MF = colonnes(*[coords(p) for p in F])
    assert MB.det() == dBp and sarrus(MB) == dBp and MF.det() == dF  # al-l1-determinant-famille-ex3 (v1, v2)
    assert (MB.inv() * MF).det() == rep == Rational(dF, dBp) and rep not in leurres  # det_B'(F) = coordonnées dans B'

# al-l1-determinant-famille-ex4 (et variantes) : vérification pour n = 2..6
for n in range(2, 7):
    Jn = Matrix.ones(n, n)
    In = eye(n)
    d0 = (Jn - In).det()
    assert d0 == (-1) ** (n - 1) * (n - 1)  # famille-ex4 : x_j = s - e_j
    assert d0 != (-1) ** n * (n - 1) and d0 != (-1) ** (n - 1) * n  # famille-ex4 (leurres, pour tout n)
    d1_ = (Jn + In).det()
    assert d1_ == n + 1 and d1_ != n  # famille-ex4 v1 (et leurre n)
    d2_ = (2 * Jn + In).det()
    assert d2_ == 2 * n + 1 and d2_ not in (2 * n, 2 * n - 1, 3 ** n)  # famille-ex4 v2
assert (Matrix.ones(2, 2) - eye(2)).det() == -1 != 2 - 1  # famille-ex4 : leurre n - 1 faux pour n = 2
assert (Matrix.ones(2, 2) + eye(2)).det() == 3 not in ((-1) ** 1 * 3, 2 ** 2)  # famille-ex4 v1 : leurres faux pour n = 2

# al-l1-determinant-produit-ex1 (et variantes)
assert (-2) ** 3 * 3 * (-2) * 3 == 144 and 144 not in (-144, 36, -18)  # produit-ex1
assert 3 ** 2 * (-1) ** 2 * 4 == 36 and 36 not in (12, -36, 108)  # produit-ex1 v1
assert (-1) ** 4 * 2 * 3 ** 2 == 18 and 18 not in (-18, 12, 36)  # produit-ex1 v2
Ad, Bd = Matrix([[3, 0, 0], [0, 1, 0], [0, 0, 1]]), Matrix([[-2, 0, 0], [0, 1, 0], [0, 0, 1]])
assert (-2 * Ad * Bd * Ad).det() == 144  # produit-ex1 (exemple concret)

# al-l1-determinant-produit-ex2 (et variantes) : valeurs nécessaires et réalisations
dsym = symbols("delta")
from sympy import solve
assert set(solve(dsym ** 3 - 4 ** 3 * dsym, dsym)) == {-8, 0, 8}  # produit-ex2
for k in (0, 2, -2):
    Ak = k * eye(3)
    assert Ak ** 3 == 4 * Ak  # produit-ex2 (réalisations)
assert {Matrix(k * eye(3)).det() for k in (0, 2, -2)} == {-8, 0, 8}
assert set(solve(dsym ** 2 - 3 ** 2 * dsym, dsym)) == {0, 9}  # produit-ex2 v1
assert (3 * eye(2)) ** 2 == 3 * (3 * eye(2)) and (3 * eye(2)).det() == 9
assert set(solve(dsym ** 2 - 2 ** 3 * dsym, dsym)) == {0, 8}  # produit-ex2 v2
assert (2 * eye(3)) ** 2 == 2 * (2 * eye(3)) and (2 * eye(3)).det() == 8
# leurres : ensembles différents de la bonne réponse
assert {-2, 0, 2} != {-8, 0, 8} and {0, 8} != {-8, 0, 8} and {-64, 0, 64} != {-8, 0, 8}
assert {0, 3} != {0, 9} and {-9, 0, 9} != {0, 9} and {9} != {0, 9}
assert {0, 2} != {0, 8} and {-8, 0, 8} != {0, 8} and {8} != {0, 8}

# al-l1-determinant-produit-ex3 (et variantes) : matrice de φ dans (1, X, X^2)
for sub, Mat_, dphi, d2, leurres in [(2 * X + 1, Matrix([[1, 1, 1], [0, 2, 4], [0, 0, 4]]), 8, 64, (8, 16, 512)),
                                     (3 * X - 1, Matrix([[1, -1, 1], [0, 3, -6], [0, 0, 9]]), 27, 729, (27, 81, 54)),
                                     (2 - X, Matrix([[1, 2, 4], [0, -1, -4], [0, 0, 1]]), -1, 1, (-1, 2, -2))]:
    Mphi = colonnes(*[coords(expand((X ** k).subs(X, sub))) for k in range(3)])
    assert Mphi == Mat_ and Mphi.det() == dphi  # al-l1-determinant-produit-ex3 (v1, v2)
    assert (Mphi * Mphi).det() == d2 and d2 not in leurres  # det(φ∘φ)
    Pb = colonnes(*[coords(expand((X - 1) ** k)) for k in range(3)])  # base (1, X - 1, (X - 1)^2)
    assert Pb.det() == 1 and (Pb.inv() * Mphi * Pb).det() == dphi  # même déterminant dans l'autre base
assert expand((X ** 2).subs(X, 2 - X).subs(X, 2 - X)) == X ** 2  # produit-ex3 v2 : φ∘φ = id

# al-l1-determinant-produit-ex4
A_ = Matrix([[1, 0], [0, -1]])
B_ = Matrix([[0, 1], [1, 0]])
assert A_.det() == -1 and B_.det() == -1 and A_ * B_ == -B_ * A_  # produit-ex4 (question 2)
assert A_ * B_ == Matrix([[0, 1], [-1, 0]]) and B_ * A_ == Matrix([[0, -1], [1, 0]])
from sympy import diag
A4_, B4_ = diag(A_, A_), diag(B_, B_)
assert A4_ * B4_ == -B4_ * A4_ and A4_.det() != 0 and B4_.det() != 0  # produit-ex4 (question 3, n = 4)
assert eye(3) * zeros(3, 3) == -zeros(3, 3) * eye(3)  # produit-ex4 (question 4 : B = 0, n impair)

# al-l1-determinant-inversibilite-ex1 (et variantes)
A1 = Matrix([[1, m, 0], [0, 1, m], [m, 0, 1]])
assert expand(A1.det()) == 1 + m ** 3 and expand(sarrus(A1) - 1 - m ** 3) == 0  # inversibilite-ex1
assert [r for r in solve(1 + m ** 3, m) if r.is_real] == [-1]  # seule racine réelle : -1
assert A1.subs(m, 1).det() != 0 and A1.subs(m, 0).det() != 0  # leurres (m = 1, m = 0 ne sont pas exclus)
A1 = Matrix([[m, 1, 1], [1, m, 1], [1, 1, m]])
assert factor(A1.det()) == (m - 1) ** 2 * (m + 2)  # inversibilite-ex1 v1
assert A1.subs(m, -1).det() != 0 and A1.subs(m, 2).det() != 0 and A1.subs(m, -2).det() == 0  # leurres v1
A1 = Matrix([[2, m, 1], [m, 2, 0], [1, 0, 1]])
assert expand(A1.det()) == 2 - m ** 2  # inversibilite-ex1 v2
assert A1.subs(m, -sqrt(2)).det() == 0 and A1.subs(m, 2).det() != 0 and A1.subs(m, 1).det() != 0  # leurres v2

# al-l1-determinant-inversibilite-ex2 (et variantes)
for n, dA, lam_, rep, leurres in [(3, 2, 2, Rational(1, 16), (4, Rational(1, 4), 16)),
                                  (2, -3, 3, Rational(-1, 27), (-3, Rational(-1, 9), Rational(1, 27))),
                                  (4, 3, -2, Rational(1, 48), (Rational(-1, 48), Rational(-1, 6), Rational(16, 3)))]:
    Aex = diag(dA, *[1] * (n - 1))
    assert Aex.det() == dA and ((lam_ * Aex).inv()).det() == rep == Rational(1, lam_ ** n * dA)  # inversibilite-ex2
    assert rep not in leurres
assert (2 * diag(2, 1, 1).inv()).det() == 4  # leurre « det(2 A^{-1}) »
assert (3 * diag(-3, 1).inv()).det() == -3  # leurre v1
assert (-2 * diag(3, 1, 1, 1).inv()).det() == Rational(16, 3)  # leurre v2

# al-l1-determinant-inversibilite-ex3 (et variantes)
for mats, dets, bonne, inv in [
    ([[[1, 2, 0], [0, 1, 1], [1, 0, 1]], [[2, 1, 0], [1, 1, 0], [0, 3, 1]], [[1, 1, 1], [1, 2, 3], [1, 3, 5]], [[2, 0, 0], [1, 1, 0], [5, 3, 1]]],
     [3, 1, 0, 2], 1, [[1, -1, 0], [-1, 2, 0], [3, -6, 1]]),
    ([[[1, 1, 0], [2, 1, 1], [0, 1, 3]], [[1, 2, 1], [0, 1, 2], [1, 0, 1]], [[3, 1, 2], [1, 1, 1], [2, 0, 1]], [[1, 3, 1], [2, 5, 0], [1, 2, 0]]],
     [-4, 4, 0, -1], 3, [[0, -2, 5], [0, 1, -2], [1, -1, 1]]),
    ([[[2, 3, 1], [1, 2, 1], [1, 1, 1]], [[2, 1, 1], [1, 2, 1], [1, 1, 2]], [[1, 0, 2], [2, 1, 3], [0, 1, 2]], [[1, 2, 3], [2, 3, 4], [3, 4, 5]]],
     [1, 4, 3, 0], 0, [[1, -2, 1], [0, 1, -1], [-1, 1, 1]])]:
    Ms = [Matrix(M_) for M_ in mats]
    assert [sarrus(M_) for M_ in Ms] == dets  # al-l1-determinant-inversibilite-ex3 (v1, v2)
    assert [abs(d) == 1 for d in dets].count(True) == 1 and abs(dets[bonne]) == 1  # une seule candidate
    assert Ms[bonne].inv() == Matrix(inv) and Ms[bonne] * Matrix(inv) == eye(3)
    for k, M_ in enumerate(Ms):
        if k != bonne and M_.det() != 0:
            assert any(not e.is_integer for e in M_.inv())  # leurres : inverse non entier

# al-l1-determinant-inversibilite-ex4 (et variantes)
for Ase, ens, leurres in [(Matrix([[1, 2, 0], [0, 3, 0], [1, 1, 2]]), {1, 2, 3}, ({1, 3}, {2, 3}, {-3, -2, -1})),
                          (Matrix([[2, 1, 0], [1, 2, 0], [0, 0, 3]]), {1, 3}, ({2, 3}, {1, 2, 3}, {-1, 3})),
                          (Matrix([[0, 1, 0], [0, 0, 1], [2, -1, 2]]), {2}, ({-1, 1, 2}, {1, 2}, set()))]:
    chi = (Ase - t * eye(3)).det()
    assert expand(chi - sarrus(Ase - t * eye(3))) == 0
    assert Poly(chi, t).degree() == 3 and Poly(chi, t).LC() == -1  # degré n, coefficient dominant (-1)^n
    assert {r for r in solve(chi, t) if r.is_real} == ens and all(ens != l_ for l_ in leurres)  # inversibilite-ex4
assert expand((Matrix([[0, 1, 0], [0, 0, 1], [2, -1, 2]]) - t * eye(3)).det() + (t - 2) * (t ** 2 + 1)) == 0  # ex4 v2
assert expand((Matrix([[2, 1, 0], [1, 2, 0], [0, 0, 3]]) - t * eye(3)).det() + (t - 3) ** 2 * (t - 1)) == 0  # ex4 v1
for n in (1, 2, 3, 4):
    S = Matrix(n, n, symbols(f"k0:{n * n}"))
    P_ = Poly((S - t * eye(n)).det(), t)
    assert P_.degree() == n and P_.LC() == (-1) ** n  # inversibilite-ex4 (question 1)

# al-l1-determinant-developpement-ex1 (et variantes)
for Ae, (i0, j0), Dm, d, leurres in [(Matrix([[1, 0, 2, -1], [3, 0, 1, 2], [0, 2, 1, 1], [2, 0, 0, 1]]), (3, 2), 5, -10, (10, -5, 5)),
                                     (Matrix([[2, 1, 0, 3], [0, 0, 3, 0], [1, 2, 1, 1], [1, 0, 2, 2]]), (2, 3), 1, -3, (3, -1, 1)),
                                     (Matrix([[1, 2, 3, 0], [0, 1, 1, 0], [2, 1, 0, -1], [1, 0, 2, 0]]), (3, 4), 1, 1, (-1, 2, -2))]:
    assert mineur(Ae, i0, j0) == Dm and Ae.det() == d == (-1) ** (i0 + j0) * Ae[i0 - 1, j0 - 1] * Dm  # developpement-ex1
    assert d not in leurres
assert Matrix([[2, -1], [1, 2]]).det() == 5 and Matrix([[1, 2], [3, 1]]).det() == -5  # ex1 (calcul de Δ32)
assert Matrix([[1, 1], [1, 2]]).det() == 1 and Matrix([[2, 3], [1, 2]]).det() == 1  # ex1 v1
assert Matrix([[1, 3], [1, 2]]).det() == -1 and Matrix([[1, 2], [1, 0]]).det() == -2  # ex1 v2

# al-l1-determinant-developpement-ex2 (et variantes) : coefficient (i, j) de M^{-1} = C_ji / det
for Me, (i0, j0), d, rep, leurres in [(Matrix([[1, 2, 0], [0, 1, 3], [4, 0, 1]]), (1, 2), 25, Rational(-2, 25), (Rational(12, 25), Rational(2, 25), Rational(-12, 25))),
                                      (Matrix([[2, 1, 1], [0, 3, 1], [1, 0, 2]]), (1, 3), 10, Rational(-1, 5), (Rational(-3, 10), Rational(1, 5), Rational(3, 10))),
                                      (Matrix([[1, 0, 2], [3, 1, 0], [0, 2, 1]]), (2, 1), 13, Rational(-3, 13), (Rational(4, 13), Rational(3, 13), Rational(-4, 13)))]:
    assert Me.det() == d and Me.inv()[i0 - 1, j0 - 1] == rep  # developpement-ex2 (v1, v2)
    Cji = (-1) ** (i0 + j0) * mineur(Me, j0, i0)
    Cij = (-1) ** (i0 + j0) * mineur(Me, i0, j0)
    assert Rational(Cji, d) == rep and Rational(Cij, d) == leurres[0]  # leurre : comatrice non transposée
    assert rep not in leurres
assert mineur(Matrix([[1, 2, 0], [0, 1, 3], [4, 0, 1]]), 2, 1) == 2  # ex2 : Δ21 = 2
assert mineur(Matrix([[2, 1, 1], [0, 3, 1], [1, 0, 2]]), 3, 1) == -2  # ex2 v1 : Δ31 = -2
assert mineur(Matrix([[1, 0, 2], [3, 1, 0], [0, 2, 1]]), 1, 2) == 3  # ex2 v2 : Δ12 = 3

# al-l1-determinant-developpement-ex3 (et variantes) : déterminants tridiagonaux, n = 1..8


def tridiag(n, d_, h_, b_):
    return Matrix(n, n, lambda i, j: d_ if i == j else (h_ if j == i + 1 else (b_ if i == j + 1 else 0)))


for (d_, h_, b_), formule, leurres in [((2, -1, -1), lambda n: n + 1, (lambda n: n, lambda n: 2 ** n, lambda n: n + 2)),
                                      ((5, 2, 2), lambda n: Rational(4 ** (n + 1) - 1, 3), (lambda n: 4 ** n + 1, lambda n: 5 * 4 ** (n - 1), lambda n: Rational(4 ** n - 1, 3))),
                                      ((3, 1, 2), lambda n: 2 ** (n + 1) - 1, (lambda n: 2 ** n + 1, lambda n: 3 ** n, lambda n: 2 ** (n + 1) + 1))]:
    vals = [tridiag(n, d_, h_, b_).det() for n in range(1, 9)]
    assert vals == [formule(n) for n in range(1, 9)]  # developpement-ex3 (v1, v2)
    for n in range(3, 9):
        assert vals[n - 1] == d_ * vals[n - 2] - h_ * b_ * vals[n - 3]  # récurrence D_n = a D_{n-1} - bc D_{n-2}
    for l_ in leurres:
        assert any(l_(n) != formule(n) for n in (1, 2))  # leurres faux (dès n = 1 ou 2)
assert tridiag(2, 2, -1, -1).det() == 3 and tridiag(2, 5, 2, 2).det() == 21 and tridiag(2, 3, 1, 2).det() == 7  # D_2

# al-l1-determinant-developpement-ex4 (et variantes)
for Ma, detf, ens, leurres in [(Matrix([[a, 1, 0], [1, a, 1], [0, 1, a]]), a ** 3 - 2 * a, {-1, 1}, ({1}, {-1, 0, 1}, set())),
                               (Matrix([[1, a, 0], [a, 1, a], [0, a, 1]]), 1 - 2 * a ** 2, {-1, 0, 1}, ({0}, {-1, 1}, set())),
                               (Matrix([[a, 2, 1], [1, a, 0], [0, 1, 1]]), a ** 2 - 1, {0}, ({-1, 1}, {-1, 0, 1}, set()))]:
    assert expand(Ma.det() - detf) == 0  # developpement-ex4 (v1, v2)
    sol = {k for k in range(-50, 51) if detf.subs(a, k) in (1, -1)}
    assert sol == ens and all(ens != l_ for l_ in leurres)
    for k in ens:
        assert all(e.is_integer for e in Ma.subs(a, k).inv())  # inverse entier
assert all(abs(k ** 3 - 2 * k) >= 4 for k in range(-50, 51) if abs(k) >= 2)  # ex4 : |a| ≥ 2 exclu

# al-l1-determinant-developpement-ex5 (et variantes) : det(Com A) = (det A)^{n-1}
for n in (2, 3):
    S = Matrix(n, n, symbols(f"c0:{n * n}"))
    Com = S.adjugate().T
    assert expand(Com.det() - S.det() ** (n - 1)) == 0  # developpement-ex5 (formule générale, n = 2, 3)
Sing = Matrix([[1, 2, 3], [2, 4, 6], [1, 0, 1]])
assert Sing.det() == 0 and Sing.adjugate().det() == 0 and Sing * Sing.adjugate() == zeros(3, 3)  # ex5 (cas det A = 0)
for n, d, rep, leurres in [(4, -2, -8, (16, 8, -2)), (3, -3, 9, (-27, -9, 27)), (5, 2, 16, (32, 8, 4))]:
    Aex = diag(d, *[1] * (n - 1))
    assert Aex.adjugate().T.det() == rep == d ** (n - 1) and rep not in leurres  # developpement-ex5 (v1, v2)

# al-l1-determinant-vandermonde-ex1 (et variantes)
assert Matrix([[1, 1, 1], [1, 2, 4], [1, 4, 16]]).det() == 6 == vdm(1, 2, 4) and 6 not in (-6, 3, 12)  # vandermonde-ex1
assert Matrix([[1, -1, 1], [1, 0, 0], [1, 3, 9]]).det() == 12 == vdm(-1, 0, 3) and 12 not in (-12, 6, 36)  # ex1 v1
M4 = Matrix([[1, 2, 4, 8], [1, 3, 9, 27], [1, 5, 25, 125], [1, 7, 49, 343]])
assert M4.det() == 240 == vdm(2, 3, 5, 7) and 240 not in (-240, 120, 30)  # vandermonde-ex1 v2

# al-l1-determinant-vandermonde-ex2 (et variantes)
c0_, c1_, c2_ = symbols("c0 c1 c2")
for pts, Pv, x0, val, leurres in [([(0, 1), (1, 3), (2, 7)], X ** 2 + X + 1, 3, 13, (11, 12, 15)),
                                  ([(-1, 2), (0, 1), (1, 4)], 2 * X ** 2 + X + 1, 2, 11, (7, 9, 13)),
                                  ([(1, 0), (2, 1), (3, 4)], (X - 1) ** 2, 4, 9, (7, 8, 16))]:
    W = Matrix([[1, p, p ** 2] for p, _ in pts])
    assert W.det() == vdm(*[p for p, _ in pts]) == 2  # vandermonde-ex2 (v1, v2) : déterminant non nul
    coef = W.solve(Matrix([q for _, q in pts]))
    assert expand(coef[0] + coef[1] * X + coef[2] * X ** 2 - Pv) == 0
    assert Pv.subs(X, x0) == val and val not in leurres

# al-l1-determinant-vandermonde-ex3 (et variantes)
V3 = (b - a) * (c - a) * (c - b)
for ligne, rep, leurres in [(lambda z: [1, z, z ** 3], (a + b + c) * V3, (a * b * c * V3, (a * b + b * c + c * a) * V3, -(a + b + c) * V3)),
                            (lambda z: [1, z ** 2, z ** 3], (a * b + b * c + c * a) * V3, ((a + b + c) * V3, a * b * c * V3, -(a * b + b * c + c * a) * V3)),
                            (lambda z: [1, z, z ** 4], (a ** 2 + b ** 2 + c ** 2 + a * b + b * c + c * a) * V3,
                             ((a + b + c) ** 2 * V3, (a ** 2 + b ** 2 + c ** 2) * V3, (a + b + c) * (a * b + b * c + c * a) * V3))]:
    Dv = Matrix([ligne(a), ligne(b), ligne(c)]).det()
    assert expand(Dv - rep) == 0  # al-l1-determinant-vandermonde-ex3 (v1, v2)
    for l_ in leurres:
        assert expand(Dv - l_) != 0  # leurres faux
assert expand((b + a) * (c ** 2 + c * a + a ** 2) - (c + a) * (b ** 2 + a * b + a ** 2) - (c - b) * (a * b + b * c + c * a)) == 0  # ex3 v1
assert expand((c ** 3 - b ** 3) + a * (c ** 2 - b ** 2) + a ** 2 * (c - b)
              - (c - b) * (c ** 2 + b * c + b ** 2 + a * (b + c) + a ** 2)) == 0  # ex3 v2

# al-l1-determinant-vandermonde-ex4 (et variantes) : exposants deux à deux distincts
for expos, interdits, leurres in [((1, a, a ** 2), {-1, 0, 1}, ({0, 1}, {1}, {-1, 1})),
                                  ((1, a, a + 2), {-1, 1}, ({1}, {-2, -1, 1}, {0, 1})),
                                  ((a, a ** 2, 2), {-sqrt(2), 0, 1, sqrt(2), 2}, ({0, 1, 2}, {0, 1, sqrt(2), 2}, {-sqrt(2), sqrt(2)}))]:
    sols = set()
    for i in range(3):
        for j in range(i + 1, 3):
            sols |= set(solve(expos[j] - expos[i], a))
    assert sols == interdits and all(interdits != l_ for l_ in leurres)  # al-l1-determinant-vandermonde-ex4 (v1, v2)
# question 1 : la matrice des dérivées en 0 est la transposée de Vandermonde
l1s, l2s, l3s = symbols("l1 l2 l3")
fs = [exp(l1s * x), exp(l2s * x), exp(l3s * x)]
Mder = Matrix(3, 3, lambda k, i: diff(fs[i], x, k).subs(x, 0))
assert expand(Mder.det() - produit_vdm(l1s, l2s, l3s)) == 0  # vandermonde-ex4 (transposée de Vandermonde)

# al-l1-determinant-vandermonde-ex5 (et variantes)


def Bk(k):
    return 1 if k == 0 else ff(X, k) / factorial(k)


for k in range(0, 7):
    for mm in range(-10, 11):
        assert (Bk(k).subs(X, mm) if k else 1).is_integer  # vandermonde-ex5 (question 3)
        if mm < 0 and k >= 1:
            assert Bk(k).subs(X, mm) == (-1) ** k * binomial(-mm + k - 1, k)
xs4 = symbols("y1:5")
Mb = Matrix(4, 4, lambda i, j: Bk(j).subs(X, xs4[i]) if j else 1)
assert simplify(Mb.det() - produit_vdm(*xs4) / 12) == 0  # vandermonde-ex5 (question 2, n = 4 : 0!1!2!3! = 12)
assert factorial(0) * factorial(1) * factorial(2) * factorial(3) == 12
for xs_, V_, rep, leurres in [((0, 2, 3, 7), 840, 70, (840, 140, 35)), ((1, 4, 5, 9), 1920, 160, (1920, 320, 80)),
                              ((-1, 0, 2, 6), 1008, 84, (1008, 168, 42))]:
    Mn = Matrix(4, 4, lambda i, j: Bk(j).subs(X, xs_[i]) if j else 1)
    assert produit_vdm(*xs_) == V_ and Mn.det() == rep == V_ / 12 and rep not in leurres  # vandermonde-ex5 (v1, v2)
import random
random.seed(1)
for _ in range(200):
    n = random.randint(2, 5)
    ent = [random.randint(-20, 20) for _ in range(n)]
    sf = 1
    for k in range(n):
        sf *= factorial(k)
    assert produit_vdm(*ent) % sf == 0  # vandermonde-ex5 (question 4, entiers au hasard)

print("al_l1_determinant : toutes les vérifications passent.")
