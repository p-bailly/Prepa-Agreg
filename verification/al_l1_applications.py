"""Vérifications calculatoires du chapitre « al-l1-applications » (thème algebre-lineaire, niveau l1).

Cartes content/cartes/algebre-lineaire/l1-al-l1-applications.json et exercices
content/exercices/algebre-lineaire/l1-al-l1-applications.json. Chaque assert cite l'identifiant de la carte ou de
l'exercice vérifié (variante comprise : « v1 », « v2 »). Les leurres des questions à choix sont vérifiés FAUX.
Lancer : python3 verification/al_l1_applications.py (doit finir sans erreur).
"""

from sympy import (Abs, Matrix, Poly, Rational, diag, diff, eye, expand, factor, integrate, simplify, solve, symbols,
                   zeros)

x, y, z, w, X, m, t, lam = symbols("x y z w X m t lam")
PREF = "al-l1-applications-"


def mat(*lignes):
    return Matrix([list(l) for l in lignes])


def meme_vect(F, G):
    """Vect(F) = Vect(G) (familles de vecteurs, éventuellement vides) dans un même espace de dimension n."""
    F = [Matrix(v) for v in F]
    G = [Matrix(v) for v in G]
    rF = Matrix.hstack(*F).rank() if F else 0
    rG = Matrix.hstack(*G).rank() if G else 0
    if not F or not G:
        return rF == rG == 0
    return rF == rG == Matrix.hstack(*F, *G).rank()


def noyau(A):
    return [list(v) for v in A.nullspace()]


def image(A):
    return [list(v) for v in A.columnspace()]


def coeffs(P, deg):
    """Coordonnées de P dans (1, X, ..., X^deg)."""
    p = Poly(expand(P), X)
    return [p.coeff_monomial(X**k) for k in range(deg + 1)]


def mat_poly(f, deg):
    """Matrice de l'endomorphisme f de R_deg[X] dans la base (1, X, ..., X^deg)."""
    return Matrix.hstack(*[Matrix(coeffs(f(X**k), deg)) for k in range(deg + 1)])


def ker_poly(phi_coords, deg):
    """Noyau d'une application linéaire sur R_deg[X] donnée par la matrice des coordonnées de phi."""
    M = Matrix.hstack(*[Matrix(phi_coords(X**k)) for k in range(deg + 1)])
    return [sum(c * X**k for k, c in enumerate(v)) for v in M.nullspace()], M


def vect_poly_egal(L1, L2, deg):
    return meme_vect([coeffs(P, deg) for P in L1], [coeffs(P, deg) for P in L2])


# ===========================================================================
# Cartes
# ===========================================================================
# al-l1-applications-noyau, exemple : f(x, y, z) = (x - y, y - z) ; noyau Vect((1, 1, 1)), image R^2.
A = mat((1, -1, 0), (0, 1, -1))
assert meme_vect(noyau(A), [(1, 1, 1)])
assert A.rank() == 2
assert meme_vect([(1, 0), (0, -1)], [(1, 0), (0, 1)])  # (1, 0) et (0, -1) engendrent R^2
# al-l1-applications-noyau, contre-exemples : x -> x^2 (f(-1) = f(1)), zéros de x^2 - 1, D(1, X) = (0, 1).
assert (-1) ** 2 == 1 ** 2 and set(solve(x**2 - 1, x)) == {-1, 1}
assert [diff(1, X), diff(X, X)] == [0, 1]
# al-l1-applications-noyau, affirmations : g o f = id avec f(x) = (x, 0), g(x, y) = x ; g(0, 1) = 0.
Ff = mat((1,), (0,))
Gg = mat((1, 0))
assert Gg * Ff == eye(1) and Gg * Matrix([0, 1]) == zeros(1, 1)
assert mat((1, 0)).rank() == 1  # (x, y) -> x envoie la base canonique sur (1, 0) : famille liée de R (2 vecteurs, rang 1)

# al-l1-applications-rang, exemple : f(x, y, z, t) = (x + y, z + t) surjective, noyau de base (1, -1, 0, 0), (0, 0, 1, -1).
A = mat((1, 1, 0, 0), (0, 0, 1, 1))
assert A.rank() == 2 and len(noyau(A)) == 2 and meme_vect(noyau(A), [(1, -1, 0, 0), (0, 0, 1, -1)])
# al-l1-applications-rang, exemple : matrices de trace nulle de dimension n^2 - 1 (n = 2, 3).
for n in (2, 3):
    trace = Matrix([[1 if (i == j) else 0 for i in range(n) for j in range(n)]])  # tr comme forme linéaire sur R^(n^2)
    assert trace.rank() == 1 and len(trace.nullspace()) == n * n - 1
# al-l1-applications-rang, contre-exemple : f(x, y) = (y, 0), Ker f = Im f = Vect((1, 0)).
N = mat((0, 1), (0, 0))
assert meme_vect(noyau(N), [(1, 0)]) and meme_vect(image(N), [(1, 0)])
# al-l1-applications-rang, trou 2 : leurres faux sur des exemples.
# x -> (x, 0) de R dans R^2 : dim E = 1, dim Ker = 0, rg = 1, dim F = 2.
assert 2 - 0 != 1 and 0 + 2 != 1  # dim F - dim Ker, dim Ker + dim F
assert 0 - 1 != 1  # application nulle sur R : rg - dim Ker = -1 ≠ dim E = 1
# al-l1-applications-rang, trou 1 (texte corrigé) : leurres « F » et « Im f » faux avec l'évaluation en 0 de R[X] dans R
# (image R, noyau contenant X, X^2, ..., X^n libres pour tout n : dimension infinie) ; leurre « Ker f » faux avec id de R[X]
# (noyau nul, image contenant 1, X, ..., X^n libres pour tout n).
for n in (5, 12):
    ev0 = Matrix([[1] + [0] * n])  # P -> P(0) sur R_n[X], dans la base (1, X, ..., X^n)
    assert ev0.rank() == 1 and len(ev0.nullspace()) == n  # le noyau grandit avec n
    assert eye(n + 1).rank() == n + 1 and eye(n + 1).nullspace() == []
# al-l1-applications-rang, contre-exemple : restrictions de f(x, y) = (y, 0) à Vect((1, 1)) (isomorphisme sur Im f) et au noyau (nulle).
assert N * Matrix([1, 1]) == Matrix([1, 0]) and N * Matrix([1, 0]) == zeros(2, 1)
# al-l1-applications-rang-matrice, cours : le rang peut être conservé par une matrice non inversible (E·E = E).
Bm1 = mat((1, 0), (0, 0))
assert Bm1 * Bm1 == Bm1 and Bm1.rank() == 1 and Bm1.det() == 0
# al-l1-applications-rang, affirmations : (x, y, z) -> (x, y) s'annule en (0, 0, 1).
assert mat((1, 0, 0), (0, 1, 0)) * Matrix([0, 0, 1]) == zeros(2, 1)

# al-l1-applications-isomorphismes, exemple : interpolation en 0, 1, 2 (matrice de Vandermonde inversible).
V = Matrix([[a**k for k in range(3)] for a in (0, 1, 2)])
assert V.det() != 0
# al-l1-applications-isomorphismes, contre-exemples : D o I = id, (I o D)(1) = 0, X*P non surjective.
I_op = lambda P: integrate(P, (X, 0, X))  # noqa: E731
for k in range(6):
    assert expand(diff(I_op(X**k), X) - X**k) == 0
assert I_op(diff(1, X)) == 0
# al-l1-applications-isomorphismes, affirmation AB = I => BA = I : exemple numérique.
Aa = mat((1, 2), (3, 5))
assert Aa * Aa.inv() == eye(2) and Aa.inv() * Aa == eye(2)

# al-l1-applications-matrice, exemple : A = [[1, 1], [0, 2]], B' = ((1, 0), (1, 1)).
A = mat((1, 1), (0, 2))
P = mat((1, 1), (0, 1))
assert P.inv() == mat((1, -1), (0, 1))
assert P.inv() * A * P == diag(1, 2)
assert A * Matrix([1, 0]) == Matrix([1, 0]) and A * Matrix([1, 1]) == Matrix([2, 2])
# al-l1-applications-matrice, contre-exemple : P A P^{-1} = [[1, 2], [0, 2]] ≠ diag(1, 2).
assert P * A * P.inv() == mat((1, 2), (0, 2)) and P * A * P.inv() != diag(1, 2)
# al-l1-applications-matrice, exemple et contre-exemple : dérivation de R_2[X].
MD = mat_poly(lambda P: diff(P, X), 2)
assert MD == mat((0, 1, 0), (0, 0, 2), (0, 0, 0))
assert MD.T * Matrix([0, 1, 0]) == Matrix([0, 0, 2])  # la transposée donnerait D(X) = 2X^2
# al-l1-applications-matrice, contre-exemple : (1, 1) a deux systèmes de coordonnées dans ((1,0),(0,1),(1,1)).
G = mat((1, 0, 1), (0, 1, 1))
assert G * Matrix([1, 1, 0]) == G * Matrix([0, 0, 1]) == Matrix([1, 1])
# al-l1-applications-matrice, trou 2 : sur un exemple (A, P, Q), les leurres diffèrent de Q^{-1} A P.
A = mat((1, 2), (3, 4))
P = mat((1, 1), (0, 1))
Q = mat((2, 1), (1, 1))
bonne = Q.inv() * A * P
for leurre in (P.inv() * A * Q, Q * A * P.inv(), P * A * Q.inv()):
    assert leurre != bonne
# al-l1-applications-matrice, affirmations : I_2 et [[1, 1], [0, 1]] (même trace, même rang), matrice de passage.
T = mat((1, 1), (0, 1))
assert T.trace() == 2 and T.rank() == 2 and T != eye(2)
Pp = mat((1, 0), (1, 1))
assert Pp.inv() == mat((1, 0), (-1, 1)) and Pp.inv() != Pp

# al-l1-applications-rang-matrice, exemple : A = [[1,2,3],[2,4,6],[1,1,2]] de rang 2, C3 = C1 + C2, L2 = 2 L1.
A = mat((1, 2, 3), (2, 4, 6), (1, 1, 2))
assert A.rank() == 2 and A[:, 2] == A[:, 0] + A[:, 1] and A[1, :] == 2 * A[0, :]
# al-l1-applications-rang-matrice, contre-exemples et affirmations.
N = mat((0, 1), (0, 0))
B = mat((1, 0), (0, 0))
assert N.rank() == 1 and (N * N).rank() == 0
assert (N * B) == zeros(2, 2) and (B * N) == N and (B * N).rank() == 1
assert mat((1, 1), (1, 1)).rank() == 1
# al-l1-applications-rang-matrice, trous : leurres faux.
assert mat((1, 1)).rank() == 1  # 2 colonnes non nulles, 2 coefficients non nuls, 2 colonnes : rang 1
assert (N + N.T).rank() == 2 != N.rank()  # A + tA
assert mat((0, 1))[:, :1].rank() == 0 != mat((0, 1)).rank()  # A privée de sa dernière colonne
# Matrices non inversibles « non nulle, carrée, triangulaire » : multiplier N par B (non nulle, triangulaire) fait chuter le rang.
assert (B * N * B).rank() == 0
# al-l1-applications-rang-matrice : rg(tA) = rg A sur des matrices rectangulaires.
for M in (mat((1, 2, 0, 1), (2, 4, 1, 3), (3, 6, 1, 4)), mat((1, -1), (2, -2), (0, 1))):
    assert M.rank() == M.T.rank()

# al-l1-applications-projecteurs, exemple : p(x, y) = (x - y, 0), s = 2p - id = (x - 2y, -y).
Pm = mat((1, -1), (0, 0))
Sm = 2 * Pm - eye(2)
assert Pm * Pm == Pm and Sm == mat((1, -2), (0, -1)) and Sm * Sm == eye(2)
assert meme_vect(image(Pm), [(1, 0)]) and meme_vect(noyau(Pm), [(1, 1)])
# al-l1-applications-projecteurs, contre-exemples : |.|, N^2 = 0, deux projecteurs de même image, Z/2Z.
for v in (-3, -1, 0, 2, 5):
    assert Abs(Abs(v)) == Abs(v)
P1 = mat((1, 0), (0, 0))
assert P1 * P1 == P1 and meme_vect(image(P1), image(Pm)) and not meme_vect(noyau(P1), noyau(Pm))
assert P1 * Matrix([0, 1]) != Pm * Matrix([0, 1])
S2 = mat((1, 1), (0, 1))
assert (S2 * S2 - eye(2)).applyfunc(lambda c: c % 2) == zeros(2, 2)
# Modulo 2, S - I = S + I = [[0, 1], [0, 0]] : même noyau Vect((1, 0)).
assert (S2 - eye(2)).applyfunc(lambda c: c % 2) == (S2 + eye(2)).applyfunc(lambda c: c % 2) == mat((0, 1), (0, 0))
# al-l1-applications-projecteurs, trou 1 et affirmation : -id vérifie p o p = id mais n'est pas un projecteur ; 2 id non plus.
assert (-eye(2)) * (-eye(2)) == eye(2) and (-eye(2)) * (-eye(2)) != -eye(2)
assert (2 * eye(2)) ** 2 != 2 * eye(2)
# id - p projecteur (affirmation vraie), sur l'exemple.
assert (eye(2) - Pm) ** 2 == eye(2) - Pm

# al-l1-applications-pivot, exemple : solutions (1, 1, 0) + z(1, -2, 1).
A = mat((1, 1, 1), (1, 2, 3), (1, 3, 5))
Bv = Matrix([2, 3, 4])
assert A.rank() == 2 and A.row_join(Bv).rank() == 2
assert A * Matrix([1, 1, 0]) == Bv and meme_vect(noyau(A), [(1, -2, 1)])
# al-l1-applications-pivot, contre-exemples : x + y = 1, x - y = 0.
A = mat((1, 1), (1, -1))
assert A.solve(Matrix([1, 0])) == Matrix([Rational(1, 2), Rational(1, 2)])
assert mat((1, 1), (0, 0)).rank() == 1  # L2 <- 0 * L2 : une droite de solutions
# Opérations simultanées : nouvelles lignes (0, 2 | 1) et (0, -2 | -1), de rang 1.
L1, L2 = mat((1, 1, 1)), mat((1, -1, 0))
assert (L1 - L2) == mat((0, 2, 1)) and (L2 - L1) == mat((0, -2, -1)) and (L1 - L2).col_join(L2 - L1).rank() == 1
# Oublier le second membre : (1, 0) ne vérifie pas x - y = 0.
assert 1 - 0 != 0
# al-l1-applications-pivot, trou 3 : A = (1 0 0), n = 1, p = 3, r = 1 : dim Ker = 2 et les leurres n - r, r, n - p valent 0, 1, -2.
A = mat((1, 0, 0))
assert len(noyau(A)) == 2 and 2 not in (1 - 1, 1, 1 - 3)
# al-l1-applications-pivot, trou 2 : rg A = n n'est pas nécessaire (exemple), B = 0 non plus, rg A = p non plus.
A = mat((1, 1), (2, 2))
assert A.rank() == A.row_join(Matrix([1, 2])).rank() == 1  # compatible avec rg A = 1 < n = 2 = p, et B ≠ 0
# al-l1-applications-pivot, affirmations : x + y + z = 0 et = 1 incompatible ; x + y = 1, 2x + 2y = 2 ; x + y = 1, x + y = 2.
A = mat((1, 1, 1), (1, 1, 1))
assert A.rank() == 1 and A.row_join(Matrix([0, 1])).rank() == 2
assert mat((1, 1), (2, 2)).rank() == mat((1, 1, 1), (2, 2, 2)).rank() == 1
assert mat((1, 1), (1, 1)).rank() == 1 and mat((1, 1, 1), (1, 1, 2)).rank() == 2

# ===========================================================================
# Exercices
# ===========================================================================
# --- noyau-ex1 (et variantes) : noyaux d'applications de R^3 dans R^2.
for nom, A, bon, leurres in [
    ("noyau-ex1", mat((1, 1, -1), (2, 1, 1)), [(-2, 3, 1)], [[], [(2, -3, 1)], [(1, 1, -1), (2, 1, 1)]]),
    ("noyau-ex1 v1", mat((1, -1, 2), (1, 1, 0)), [(-1, 1, 1)], [[], [(1, -1, 1)], [(1, -1, 2), (1, 1, 0)]]),
    ("noyau-ex1 v2", mat((1, 2, 3), (2, 4, 6)), [(-2, 1, 0), (-3, 0, 1)], [[(1, 2, 3)], [(-2, 1, 0)], []]),
]:
    K = noyau(A)
    assert meme_vect(K, bon), nom
    for L in leurres:
        assert not (meme_vect(L, K) if L else len(K) == 0), (nom, L)

# --- noyau-ex2 (et variantes) : noyaux d'applications sur R_2[X].
cas = [
    ("noyau-ex2", lambda P: [P.subs(X, 1), diff(P, X).subs(X, 1)], (X - 1) ** 2, [X - 1, X**2 - 1]),
    ("noyau-ex2 v1", lambda P: [P.subs(X, 0), P.subs(X, 2)], X**2 - 2 * X, [X - 2]),
    ("noyau-ex2 v2", lambda P: coeffs(P - X * diff(P, X), 2), X, [1 + 0 * X, X**2]),
]
for nom, phi, bon, leurres in cas:
    K, M = ker_poly(lambda P: phi(expand(P)), 2)
    assert len(K) == 1 and vect_poly_egal(K, [bon], 2), nom
    for L in leurres:
        assert not vect_poly_egal(K, [L], 2), (nom, L)
    assert len(K) != 0  # leurre {0} faux
K, M = ker_poly(lambda P: [P.subs(X, 0), P.subs(X, 2)], 2)
assert not vect_poly_egal(K, [X, X - 2], 2)  # noyau-ex2 v1 : leurre Vect(X, X - 2)
assert expand(X**2 - 2 * X - X * (X - 2)) == 0  # noyau-ex2 v1 : X^2 - 2X = X(X - 2)
assert expand((X - 1) ** 2 - (1 - 2 * X + X**2)) == 0  # noyau-ex2


# --- noyau-ex3 (et variantes) : M -> AM sur M_2(R) ou M_{2,3}(R) ; dimension du noyau = (nb colonnes) x dim Ker A.
def dim_ker_mult_gauche(A, ncol):
    a = symbols(f"a0:{A.cols * ncol}")
    M = Matrix(A.cols, ncol, a)
    eqs = list(A * M)
    J = Matrix([[diff(e, v) for v in a] for e in eqs])
    return len(a) - J.rank()


assert dim_ker_mult_gauche(mat((1, 2), (2, 4)), 2) == 2  # noyau-ex3 : réponse 2 (leurres 1, 0, 4 faux)
assert mat((1, 2), (2, 4)) * Matrix([-2, 1]) == zeros(2, 1)  # noyau-ex3 : colonne (-2, 1)
assert mat((1, 2), (2, 4)) * mat((-2, 0), (1, 0)) == zeros(2, 2) and mat((1, 2), (2, 4)) * mat((0, -2), (0, 1)) == zeros(2, 2)
assert dim_ker_mult_gauche(mat((1, 2), (3, 4)), 2) == 0  # noyau-ex3 v1 : réponse 0
assert dim_ker_mult_gauche(mat((1, 2), (2, 4)), 3) == 3  # noyau-ex3 v2 : réponse 3 (leurres 1, 2, 6 faux)

# --- noyau-ex4 (remarque) : D sur R[X] ; Ker D = constantes, Ker D^2 = R_1[X], Im D = Im D^2 = R[X] (primitives).
assert diff(X, X, 2) == 0 and diff(X, X) != 0
for k in range(5):
    assert diff(X ** (k + 2) / ((k + 1) * (k + 2)), X, 2) == X**k  # X^k a un antécédent par D^2

# --- rang-ex1 (et variantes) : dim Ker = dim E - rg.
assert 5 - 3 == 2 and 2 not in (3, 5, 8)  # rang-ex1
assert 4 - 3 == 1 and 1 not in (3, 4, 7)  # rang-ex1 v1 : dim M_2(R) = 4
assert 5 - 2 == 3 and 3 not in (2, 1, 5)  # rang-ex1 v2 : dim R_4[X] = 5

# --- rang-ex2 (et variantes) : noyau, rang, équation de l'image.
for nom, A, k, n_bon, leurres in [
    ("rang-ex2", mat((1, -1, 0), (0, 1, -1), (-1, 0, 1)), (1, 1, 1), (1, 1, 1), [(1, -1, 0)]),
    ("rang-ex2 v1", mat((1, 2, 1), (2, 1, -1), (1, -1, -2)), (1, -1, 1), (1, -1, 1), [(1, 1, 1)]),
    ("rang-ex2 v2", mat((1, 1, 0), (0, 1, 1), (1, 2, 1)), (1, -1, 1), (1, 1, -1), [(1, -1, 1)]),
]:
    assert meme_vect(noyau(A), [k]) and A.rank() == 2, nom
    assert Matrix([n_bon]) * A == zeros(1, 3), nom  # l'équation est vérifiée par toutes les images
    for nl in leurres:
        assert Matrix([nl]) * A != zeros(1, 3), (nom, nl)  # plan leurre : ne contient pas l'image
    assert A.rank() != 3 and A.rank() != 1, nom  # leurres « R^3 » et « droite »
# rang-ex2 : la droite Vect((1, 1, 1)) n'est pas l'image (rang 2) ; v1, v2 : (1, -1, 1) n'est pas dans l'image.
assert mat((1, 2, 1), (2, 1, -1), (1, -1, -2)).row_join(Matrix([1, -1, 1])).rank() == 3
assert mat((1, 1, 0), (0, 1, 1), (1, 2, 1)).row_join(Matrix([1, -1, 1])).rank() == 3

# --- rang-ex3 (et variantes) : base du noyau d'une évaluation sur R_3[X].
for nom, phi, preimages, base, leurres in [
    ("rang-ex3", lambda P: [P.subs(X, 1), P.subs(X, -1)], [(1 + X) / 2, (1 - X) / 2], [X**2 - 1, X**3 - X],
     [[X - 1, X + 1], [X**2 - 1], [X**2 - 1, X**3 - 1]]),
    ("rang-ex3 v1", lambda P: [P.subs(X, 0), diff(P, X).subs(X, 0)], [1 + 0 * X, X], [X**2, X**3],
     [[1 + 0 * X, X], [X, X**2], [X**3]]),
    ("rang-ex3 v2", lambda P: [P.subs(X, 2), P.subs(X, 3)], [3 - X, X - 2], [X**2 - 5 * X + 6, X**3 - 5 * X**2 + 6 * X],
     [[X - 2, X - 3], [X**2 - 5 * X + 6], [X**2 - 5 * X + 6, X**3 - 6]]),
]:
    assert [phi(expand(P)) for P in preimages] == [[1, 0], [0, 1]], nom
    K, M = ker_poly(lambda P: phi(expand(P)), 3)
    assert len(K) == 2 and vect_poly_egal(K, base, 3), nom
    assert Matrix.hstack(*[Matrix(coeffs(P, 3)) for P in base]).rank() == 2, nom  # famille libre
    for L in leurres:
        est_base = len(L) == 2 and vect_poly_egal(K, L, 3) and Matrix.hstack(*[Matrix(coeffs(P, 3)) for P in L]).rank() == 2
        assert not est_base, (nom, L)

# --- rang-ex4 (et variantes) : rang maximal d'un endomorphisme de carré nul en dimension n : n // 2, atteint.
for nom, n, bon, leurres in [("rang-ex4", 7, 3, (4, 6, 2)), ("rang-ex4 v1", 6, 3, (2, 5, 6)), ("rang-ex4 v2", 9, 4, (5, 8, 3))]:
    r = n // 2
    assert r == bon, nom
    F = zeros(n, n)
    for i in range(r):
        F[i + r, i] = 1  # f(e_{i+1}) = e_{i+1+r}
    assert F * F == zeros(n, n) and F.rank() == r, nom
    for L in leurres:
        assert L != bon and (2 * L > n or L < bon), (nom, L)  # rang L impossible (2L > n) ou non maximal

# --- rang-ex5 : exemple d'égalité f = g = N sur R^2.
N = mat((0, 1), (0, 0))
assert N * N == zeros(2, 2) and N.rank() + N.rank() == 2
# Formule rg(g o f) = rg f - dim(Ker g ∩ Im f) testée sur des matrices.
for Fm, Gm in [(mat((1, 0), (0, 1), (1, 1)), mat((1, -1, 0), (0, 0, 1))), (mat((1, 2), (2, 4), (0, 0)), mat((2, -1, 0), (0, 0, 1)))]:
    dim_inter = Fm.rank() + len(noyau(Gm)) - Matrix.hstack(*[Matrix(v) for v in image(Fm) + noyau(Gm)]).rank()
    assert (Gm * Fm).rank() == Fm.rank() - dim_inter
    assert (Gm * Fm).rank() >= Fm.rank() + Gm.rank() - Gm.cols

# --- isomorphismes-ex1 (et variantes) : bijectivité par le noyau.
assert mat((1, 1, 0), (0, 1, 1), (1, 0, 1)).det() != 0  # isomorphismes-ex1 : bijectif (réponse « Oui »)
assert meme_vect(noyau(mat((1, -1, 0), (0, 1, -1), (-1, 0, 1))), [(1, 1, 1)])  # v1 : noyau une droite
A = mat((1, 2, 0), (0, 1, 2), (2, 0, 1))
assert A.det() == 9 and noyau(A) == []  # v2 : x = -8x, 9x = 0 ; bijectif

# --- isomorphismes-ex2 (et variantes) : dimensions.
assert (3 + 1, 2 * 2) == (4, 4)  # R_3[X] et M_2(R) : dimension 4
assert (5 + 1, 2 * 3) == (6, 6)  # v1 : R_5[X] et M_{2,3}(R) : dimension 6
assert (5, 2 * 2) == (5, 4)  # v2 : R^5 et M_2(R) non isomorphes

# --- isomorphismes-ex3 (et variantes) : interpolation.
a0, a1, a2 = symbols("a0 a1 a2")
for nom, pts, vals, bon, leurres in [
    ("isomorphismes-ex3", (0, 1, 2), (1, 2, 5), X**2 + 1, [2 * X**2 - X + 1, X**2 + 2, 3 * X + 1]),
    ("isomorphismes-ex3 v1", (0, 1, 2), (2, 1, 4), 2 * X**2 - 3 * X + 2, [2 * X**2 + 3 * X + 2, X**2 - 2 * X + 2, -X + 2]),
    ("isomorphismes-ex3 v2", (-1, 0, 1), (3, 1, 3), 2 * X**2 + 1, [X**2 + 2, 2 * X**2 + X + 1, 3 * X**2 + 1]),
]:
    P = a0 + a1 * X + a2 * X**2
    sol = solve([P.subs(X, p) - v for p, v in zip(pts, vals)], [a0, a1, a2])
    assert expand(P.subs(sol) - bon) == 0, nom
    for L in leurres:
        assert [L.subs(X, p) for p in pts] != list(vals), (nom, L)

# --- isomorphismes-ex4 : I(D(P)) = P - P(0) ; D o I = id (déjà vérifié pour la carte).
for P in (X**3 + 2 * X + 5, 7 + 0 * X, X**4):
    assert expand(I_op(diff(P, X)) - (P - P.subs(X, 0))) == 0

# --- matrice-ex1 (et variantes) : matrices dans (1, X, X^2).
for nom, f, bon, leurres in [
    ("matrice-ex1", lambda P: diff(P, X), mat((0, 1, 0), (0, 0, 2), (0, 0, 0)),
     [mat((0, 0, 0), (1, 0, 0), (0, 2, 0)), mat((0, 1, 0), (0, 0, 1), (0, 0, 0)), mat((1, 0, 0), (0, 2, 0), (0, 0, 0))]),
    ("matrice-ex1 v1", lambda P: P.subs(X, X + 1), mat((1, 1, 1), (0, 1, 2), (0, 0, 1)),
     [mat((1, 0, 0), (1, 1, 0), (1, 2, 1)), mat((1, 1, 0), (0, 1, 1), (0, 0, 1)), mat((1, 1, 1), (0, 1, 1), (0, 0, 1))]),
    ("matrice-ex1 v2", lambda P: X * diff(P, X), diag(0, 1, 2),
     [mat((0, 1, 0), (0, 0, 2), (0, 0, 0)), mat((0, 0, 0), (1, 0, 0), (0, 2, 0)), diag(1, 2, 3)]),
]:
    M = mat_poly(f, 2)
    assert M == bon, nom
    for L in leurres:
        assert L != M, (nom, L)

# --- matrice-ex2 (et variantes) : P^{-1} A P.
for nom, A, Pm, bon, leurres in [
    ("matrice-ex2", mat((4, -2), (1, 1)), mat((1, 2), (1, 1)), diag(2, 3), [diag(3, 2), mat((-6, 12), (-6, 11)), mat((4, 1), (-2, 1))]),
    ("matrice-ex2 v1", mat((5, -6), (3, -4)), mat((1, 2), (1, 1)), diag(-1, 2), [diag(2, -1), mat((-25, 36), (-18, 26)), mat((5, 3), (-6, -4))]),
    ("matrice-ex2 v2", mat((2, 1), (1, 2)), mat((1, 1), (1, 0)), mat((3, 1), (0, 1)), [mat((3, 0), (1, 1)), mat((2, 1), (1, 2)), mat((1, 0), (1, 3))]),
]:
    assert Pm.det() != 0 and Pm.inv() * A * Pm == bon, nom
    for L in leurres:
        assert L != bon, (nom, L)
assert mat((1, 2), (1, 1)).inv() == mat((-1, 2), (1, -1))  # matrice-ex2 et v1 : P^{-1}
assert mat((1, 1), (1, 0)).inv() == mat((0, 1), (1, -1))  # matrice-ex2 v2 : P^{-1}
assert mat((1, 2), (1, 1)) * mat((4, -2), (1, 1)) * mat((1, 2), (1, 1)).inv() == mat((-6, 12), (-6, 11))  # leurre P A P^{-1}
assert mat((1, 2), (1, 1)) * mat((5, -6), (3, -4)) * mat((1, 2), (1, 1)).inv() == mat((-25, 36), (-18, 26))  # v1
# matrice-ex2 : images f(u1), f(u2).
assert mat((4, -2), (1, 1)) * Matrix([1, 1]) == Matrix([2, 2]) and mat((4, -2), (1, 1)) * Matrix([2, 1]) == Matrix([6, 3])
assert mat((5, -6), (3, -4)) * Matrix([1, 1]) == Matrix([-1, -1]) and mat((5, -6), (3, -4)) * Matrix([2, 1]) == Matrix([4, 2])
assert mat((2, 1), (1, 2)) * Matrix([1, 1]) == Matrix([3, 3]) and mat((2, 1), (1, 2)) * Matrix([1, 0]) == Matrix([2, 1])

# --- matrice-ex3 (et variantes) : Q^{-1} A P et leurres (A P, Q A P, Q^{-1} A).
h = Rational(1, 2)
for nom, A, Pm, Q, Qi, bon, AP in [
    ("matrice-ex3", mat((1, 1, 1), (1, -1, 0)), mat((1, 1, 1), (0, 1, 1), (0, 0, 1)), mat((1, 1), (1, -1)),
     h * mat((1, 1), (1, -1)), mat((1, 1, 3 * h), (0, 1, 3 * h)), mat((1, 2, 3), (1, 0, 0))),
    ("matrice-ex3 v1", mat((1, 0, -1), (0, 1, 1)), mat((1, 0, 1), (1, 1, 0), (0, 1, 1)), mat((1, 1), (0, 1)),
     mat((1, -1), (0, 1)), mat((0, -3, -1), (1, 2, 1)), mat((1, -1, 0), (1, 2, 1))),
    ("matrice-ex3 v2", mat((2, 1, 0), (0, 1, -1)), mat((1, 1, 1), (0, 1, 1), (0, 0, 1)), mat((1, 0), (2, 1)),
     mat((1, 0), (-2, 1)), mat((2, 3, 3), (-4, -5, -6)), mat((2, 3, 3), (0, 1, 0))),
]:
    assert Q.inv() == Qi and A * Pm == AP and Pm.det() != 0, nom
    assert Q.inv() * A * Pm == bon, nom
    for L in (A * Pm, Q * A * Pm, Q.inv() * A):
        assert L != bon, nom
assert mat((1, 1), (1, -1)) * mat((1, 2, 3), (1, 0, 0)) == mat((2, 2, 3), (0, 2, 3))  # matrice-ex3 : leurre Q A P
assert h * mat((1, 1), (1, -1)) * mat((1, 1, 1), (1, -1, 0)) == mat((1, 0, h), (0, 1, h))  # matrice-ex3 : leurre Q^{-1} A
assert mat((1, 1), (0, 1)) * mat((1, -1, 0), (1, 2, 1)) == mat((2, 1, 1), (1, 2, 1))  # v1 : Q A P
assert mat((1, -1), (0, 1)) * mat((1, 0, -1), (0, 1, 1)) == mat((1, -1, -2), (0, 1, 1))  # v1 : Q^{-1} A
assert mat((1, 0), (2, 1)) * mat((2, 3, 3), (0, 1, 0)) == mat((2, 3, 3), (4, 7, 6))  # v2 : Q A P
assert mat((1, 0), (-2, 1)) * mat((2, 1, 0), (0, 1, -1)) == mat((2, 1, 0), (-4, -1, -1))  # v2 : Q^{-1} A
# Décompositions dans C' : (a, b) en coordonnées.
a, b = symbols("a b")
assert expand((a + b) / 2 * Matrix([1, 1]) + (a - b) / 2 * Matrix([1, -1]) - Matrix([a, b])) == zeros(2, 1)  # matrice-ex3
assert expand((a - b) * Matrix([1, 0]) + b * Matrix([1, 1]) - Matrix([a, b])) == zeros(2, 1)  # v1
assert expand(a * Matrix([1, 2]) + (b - 2 * a) * Matrix([0, 1]) - Matrix([a, b])) == zeros(2, 1)  # v2

# --- matrice-ex4 (et variantes) : M^2 = 0, P = (Mx | x) avec x = e1 ; leurres faux.
Nn = mat((0, 1), (0, 0))
for nom, M, bon, leurres in [
    ("matrice-ex4", mat((2, -4), (1, -2)), mat((2, 1), (1, 0)), [mat((1, 2), (0, 1)), mat((2, -4), (1, -2)), mat((1, 0), (2, 1))]),
    ("matrice-ex4 v1", mat((3, 9), (-1, -3)), mat((3, 1), (-1, 0)), [mat((1, 3), (0, -1)), mat((3, 9), (-1, -3)), mat((1, 0), (3, 1))]),
    ("matrice-ex4 v2", mat((-1, 1), (-1, 1)), mat((-1, 1), (-1, 0)), [mat((1, -1), (0, -1)), mat((-1, 1), (-1, 1)), mat((1, 0), (-1, 1))]),
]:
    assert M * M == zeros(2, 2) and M != zeros(2, 2), nom
    assert bon[:, 0] == M * Matrix([1, 0]) and bon.det() != 0 and bon.inv() * M * bon == Nn, nom
    for L in leurres:
        assert L.det() == 0 or L.inv() * M * L != Nn, (nom, L)

# --- rang-matrice-ex1 (et variantes).
for nom, A, r in [("rang-matrice-ex1", mat((1, 2, 0, 1), (2, 4, 1, 3), (3, 6, 1, 4)), 2),
                  ("rang-matrice-ex1 v1", mat((1, 1, 0, 2), (0, 1, 1, 0), (1, 0, 1, 1)), 3),
                  ("rang-matrice-ex1 v2", mat((1, -2, 3, 1), (-2, 4, -6, -2), (3, -6, 9, 3)), 1)]:
    assert A.rank() == r, nom  # leurres : les trois autres valeurs de {1, 2, 3, 4}
A = mat((1, 2, 0, 1), (2, 4, 1, 3), (3, 6, 1, 4))
assert A[2, :] == A[0, :] + A[1, :] and A[:, 1] == 2 * A[:, 0] and A[:, 3] == A[:, 0] + A[:, 2]  # rang-matrice-ex1
E3 = mat((1, 1, 0, 2), (0, 1, 1, 0), (1, 0, 1, 1))
E3[2, :] = E3[2, :] - E3[0, :] + E3[1, :]
assert E3[2, :] == mat((0, 0, 2, -1))  # rang-matrice-ex1 v1 : forme échelonnée

# --- rang-matrice-ex2 (et variantes) : dim Ker A = p - r, dim Ker tA = n - r (réponse), sur des matrices exemples.
for nom, n, p, r, props, bon in [("rang-matrice-ex2", 4, 3, 2, (2, 1, 3, 4), 2), ("rang-matrice-ex2 v1", 3, 5, 3, (0, 2, 3, 5), 0),
                                 ("rang-matrice-ex2 v2", 5, 2, 1, (4, 1, 3, 5), 4)]:
    A = Matrix(n, p, lambda i, j: 1 if (i == j and i < r) else 0)
    assert A.rank() == r and len(A.nullspace()) == p - r and len(A.T.nullspace()) == n - r == bon, nom
    assert props[0] == bon and len(set(props)) == 4, nom

# --- rang-matrice-ex3 (et variantes) : valeurs de m pour lesquelles le rang vaut 2.
for nom, A, bon, leurres in [
    ("rang-matrice-ex3", mat((1, 1, m), (1, m, 1), (m, 1, 1)), {-2}, [{1, -2}, {1}, set()]),
    ("rang-matrice-ex3 v1", mat((1, 2, 1), (2, m, 0), (1, 0, m)), {0, 5}, [{5}, {0}, set()]),
    ("rang-matrice-ex3 v2", mat((1, 1, 1), (1, m, m**2), (1, 2, 4)), {1, 2}, [{1}, {2}, set()]),
]:
    racines = set(solve(A.det(), m))
    rang2 = {r for r in racines if A.subs(m, r).rank() == 2}
    assert rang2 == bon, nom
    for L in leurres:
        assert L != bon, (nom, L)
    assert A.subs(m, 7).rank() == 3, nom  # hors des racines du déterminant, rang 3
assert mat((1, 1, m), (1, m, 1), (m, 1, 1)).subs(m, 1).rank() == 1  # rang-matrice-ex3 : m = 1 donne le rang 1
assert factor(mat((1, 1, m), (1, m, 1), (m, 1, 1)).det()) == factor(-(m - 1) ** 2 * (m + 2))
assert expand(-4 + (m - 4) * (m - 1) - (m**2 - 5 * m)) == 0  # rang-matrice-ex3 v1 : dernier coefficient
assert expand(m**2 - 1 - 3 * (m - 1) - (m - 1) * (m - 2)) == 0  # rang-matrice-ex3 v2 : dernier coefficient

# --- rang-matrice-ex4 (et variantes) : A^k = tr(A)^(k-1) A ; leurres faux (au moins pour k = 2).
k = symbols("k", integer=True, positive=True)
for nom, A, U, Vv, tr, leurres in [
    ("rang-matrice-ex4", mat((1, 2), (2, 4)), Matrix([1, 2]), Matrix([1, 2]), 5, [lambda j: 5**j, lambda j: 1, lambda j: 4 ** (j - 1)]),
    ("rang-matrice-ex4 v1", mat((1, -1, 2), (-1, 1, -2), (2, -2, 4)), Matrix([1, -1, 2]), Matrix([1, -1, 2]), 6,
     [lambda j: 6**j, lambda j: 1, lambda j: 4 ** (j - 1)]),
    ("rang-matrice-ex4 v2", mat((2, -6), (1, -3)), Matrix([2, 1]), Matrix([1, -3]), -1,
     [lambda j: 1, lambda j: (-1) ** j, lambda j: 5 ** (j - 1)]),
]:
    assert U * Vv.T == A and A.rank() == 1 and A.trace() == tr and (Vv.T * U)[0, 0] == tr, nom
    for j in range(1, 7):
        assert A**j == tr ** (j - 1) * A, (nom, j)
    for L in leurres:
        assert any(A**j != L(j) * A for j in range(1, 4)), nom

# --- projecteurs-ex1 (et variantes).
for nom, Pm, im, ke, leurres in [
    ("projecteurs-ex1", mat((1, 1), (0, 0)), (1, 0), (1, -1), [((1, -1), (1, 0)), ((1, 0), (0, 1))]),
    ("projecteurs-ex1 v1", mat((2, -1), (2, -1)), (1, 1), (1, 2), [((1, 2), (1, 1)), ((1, 1), (1, -1))]),
    ("projecteurs-ex1 v2", mat((0, 0), (1, 1)), (0, 1), (1, -1), [((1, -1), (0, 1)), ((0, 1), (1, 0))]),
]:
    assert Pm * Pm == Pm, nom  # leurre « n'est pas un projecteur » faux
    assert meme_vect(image(Pm), [im]) and meme_vect(noyau(Pm), [ke]), nom
    for (i2, k2) in leurres:
        assert not (meme_vect(image(Pm), [i2]) and meme_vect(noyau(Pm), [k2])), (nom, i2, k2)

# --- projecteurs-ex2 (et variantes) : projection sur un plan parallèlement à une droite.
for nom, nvec, d, v, bon, leurres in [
    ("projecteurs-ex2", (1, 1, 1), (1, 0, 0), (1, 2, 3), (-5, 2, 3), [(6, 0, 0), (-1, 0, 1), (-11, 2, 3)]),
    ("projecteurs-ex2 v1", (1, -1, 2), (0, 0, 1), (1, 3, 2), (1, 3, 1), [(0, 0, 1), (1, 3, 0), (Rational(2, 3), Rational(10, 3), Rational(4, 3))]),
    ("projecteurs-ex2 v2", (2, 1, -1), (1, 1, 1), (1, 0, 0), (0, -1, -1), [(1, 1, 1), (-1, -2, -2), (Rational(1, 3), Rational(-1, 3), Rational(1, 3))]),
]:
    nn, dd, vv = Matrix(nvec), Matrix(d), Matrix(v)
    assert (nn.T * dd)[0, 0] != 0, nom  # D non incluse dans P : somme directe
    tt = (nn.T * vv)[0, 0] / (nn.T * dd)[0, 0]
    pv = vv - tt * dd
    assert pv == Matrix(bon) and (nn.T * pv)[0, 0] == 0, nom
    for L in leurres:
        assert Matrix(L) != pv, (nom, L)
assert 2 * Matrix([-5, 2, 3]) - Matrix([1, 2, 3]) == Matrix([-11, 2, 3])  # projecteurs-ex2 : leurre s(v)
assert Matrix([1, 2, 3]) - 2 * Matrix([1, 1, 1]) == Matrix([-1, 0, 1])  # projecteurs-ex2 : leurre projection orthogonale
assert 2 * Matrix([1, 3, 1]) - Matrix([1, 3, 2]) == Matrix([1, 3, 0])  # v1 : leurre s(v)
assert Matrix([1, 3, 2]) - Rational(1, 3) * Matrix([1, -1, 2]) == Matrix([Rational(2, 3), Rational(10, 3), Rational(4, 3)])  # v1
assert 2 * Matrix([0, -1, -1]) - Matrix([1, 0, 0]) == Matrix([-1, -2, -2])  # v2 : leurre s(v)
assert Matrix([1, 0, 0]) - Rational(2, 6) * Matrix([2, 1, -1]) == Matrix([Rational(1, 3), Rational(-1, 3), Rational(1, 3)])  # v2
# Expressions générales de p.
assert simplify(Matrix([x, y, z]) - (x + y + z) * Matrix([1, 0, 0]) - Matrix([-y - z, y, z])) == zeros(3, 1)  # projecteurs-ex2
assert simplify(Matrix([x, y, z]) - (x - y + 2 * z) / 2 * Matrix([0, 0, 1]) - Matrix([x, y, (y - x) / 2])) == zeros(3, 1)  # v1

# --- projecteurs-ex3 (et variantes) : parties symétrique et antisymétrique.
for nom, M, S, Aa, leurres in [
    ("projecteurs-ex3", mat((1, 2), (4, 3)), mat((1, 3), (3, 3)), mat((0, -1), (1, 0)), [mat((0, -1), (1, 0)), mat((1, 2), (2, 3)), mat((2, 6), (6, 6))]),
    ("projecteurs-ex3 v1", mat((0, 5), (1, 2)), mat((0, 3), (3, 2)), mat((0, 2), (-2, 0)), [mat((0, 2), (-2, 0)), mat((0, 5), (5, 2)), mat((0, 6), (6, 4))]),
    ("projecteurs-ex3 v2", mat((1, 0, 2), (4, 1, 0), (0, 2, 3)), mat((1, 2, 1), (2, 1, 1), (1, 1, 3)), mat((0, -2, 1), (2, 0, -1), (-1, 1, 0)),
     [mat((0, -2, 1), (2, 0, -1), (-1, 1, 0)), mat((2, 4, 2), (4, 2, 2), (2, 2, 6)), mat((1, 0, 2), (0, 1, 0), (2, 0, 3))]),
]:
    assert (M + M.T) / 2 == S and (M - M.T) / 2 == Aa and S + Aa == M and S.T == S and Aa.T == -Aa, nom
    for L in leurres:
        assert L != S, (nom, L)

# --- projecteurs-ex4 : exemples numériques (p q = q p = 0 => p + q projecteur ; sinon non).
p1 = mat((1, 0, 0), (0, 0, 0), (0, 0, 0))
q1 = mat((0, 0, 0), (0, 1, 0), (0, 0, 0))
assert p1 * q1 == q1 * p1 == zeros(3, 3) and (p1 + q1) ** 2 == p1 + q1
q2 = mat((1, 1), (0, 0))
p2 = mat((1, 0), (0, 0))
assert (p2 + q2) ** 2 != p2 + q2 and p2 * q2 != zeros(2, 2)

# --- projecteurs-ex5 (et variantes) : P^2 = P et rg P = tr P.
for nom, Pm, r, props in [
    ("projecteurs-ex5", mat((2, -2, -4), (-1, 3, 4), (1, -2, -3)), 2, (2, 3, 1, 0)),
    ("projecteurs-ex5 v1", mat((0, -1, -2), (-2, -1, -4), (1, 1, 3)), 2, (2, 3, 1, 0)),
    ("projecteurs-ex5 v2", mat((1, 1, 0, -1), (1, 0, 1, 0), (0, -1, 1, 1), (1, 0, 1, 0)), 2, (2, 3, 1, 4)),
]:
    assert Pm * Pm == Pm and Pm.trace() == r == Pm.rank(), nom
    assert props[0] == r and len(set(props)) == 4, nom
assert eye(3) - Matrix([1, 2, -1]) * Matrix([[1, 1, 2]]) == mat((0, -1, -2), (-2, -1, -4), (1, 1, 3))  # v1 : I - u tv

# --- pivot-ex1 (et variantes) : systèmes de Cramer.
for nom, A, B, bon, leurres in [
    ("pivot-ex1", mat((1, 1, 1), (1, -1, 2), (2, 1, -1)), (6, 5, 1), (1, 2, 3), [(1, 3, 2), (2, 1, 3), (3, 2, 1)]),
    ("pivot-ex1 v1", mat((1, 2, -1), (2, -1, 1), (1, 1, 3)), (-1, 6, 4), (2, -1, 1), [(1, -1, 2), (2, 1, -1), (-1, 2, 1)]),
    ("pivot-ex1 v2", mat((1, 1, 1), (2, 3, 1), (1, -1, 2)), (1, 0, 3), (-1, 0, 2), [(1, 0, 0), (0, -1, 2), (2, 0, -1)]),
]:
    assert A.det() != 0 and A.solve(Matrix(B)) == Matrix(bon), nom
    for L in leurres:
        assert A * Matrix(L) != Matrix(B), (nom, L)

# --- pivot-ex2 (et variantes) : solutions X0 + Ker A.
for nom, A, B, X0, base, leurres in [
    ("pivot-ex2", mat((1, 1, 1), (1, 2, 3), (2, 3, 4)), (3, 6, 9), (0, 3, 0), [(1, -2, 1)], [((1, 1, 1), []), ((0, 3, 0), [(1, 2, 1)])]),
    ("pivot-ex2 v1", mat((1, -1, 1), (2, -1, 3), (1, 0, 2)), (1, 4, 3), (3, 2, 0), [(-2, -1, 1)], [((3, 2, 0), [(2, 1, 1)]), ((3, 2, 0), [])]),
    ("pivot-ex2 v2", mat((1, 1, 1, 1), (1, -1, 1, -1)), (2, 0), (1, 1, 0, 0), [(-1, 0, 1, 0), (0, -1, 0, 1)],
     [((1, 1, 0, 0), [(-1, 0, 1, 0)]), ((1, 1, 0, 0), [(1, 0, 1, 0), (0, 1, 0, 1)])]),
]:
    assert A * Matrix(X0) == Matrix(B) and meme_vect(noyau(A), base) and len(noyau(A)) == len(base), nom
    assert A.rank() == A.row_join(Matrix(B)).rank(), nom  # compatible : le leurre « ensemble vide » est faux
    for (x0, b2) in leurres:
        # L'ensemble leurre x0 + Vect(b2) est différent : soit x0 n'est pas solution, soit la direction diffère.
        assert A * Matrix(x0) != Matrix(B) or not (meme_vect(b2, noyau(A)) if b2 else False), (nom, x0, b2)
assert mat((1, 1, 1), (1, 2, 3), (2, 3, 4)) * Matrix([1, 1, 1]) == Matrix([3, 6, 9])  # pivot-ex2 : (1,1,1) solution, mais pas la seule

# --- pivot-ex3 (et variantes) : compatibilité.
for nom, A, B, bon, X0, d in [
    ("pivot-ex3", mat((1, 1, 1), (1, 2, 3), (2, 3, 4)), (1, 2, m), 3, (0, 1, 0), (1, -2, 1)),
    ("pivot-ex3 v1", mat((1, 1, -1), (1, -1, 2), (1, 3, -4)), (2, 0, m), 4, (1, 1, 0), (-1, 3, 2)),
    ("pivot-ex3 v2", mat((2, 1, 0), (1, 0, -1), (3, 1, -1)), (1, -1, m), 0, (0, 1, 1), (1, -2, 1)),
]:
    assert A.rank() == 2, nom
    Bm = Matrix(B)
    for val in (bon, bon + 1, bon - 2, 10):
        compatible = A.row_join(Bm.subs(m, val)).rank() == 2
        assert compatible == (val == bon), (nom, val)
    assert A * Matrix(X0) == Bm.subs(m, bon) and meme_vect(noyau(A), [d]), nom

# --- pivot-ex4 (et variantes) : équations de F = Vect(u1, u2).
for nom, u1, u2, bon, leurres in [
    ("pivot-ex4", (1, 1, 0, 1), (0, 1, 1, 2), [(1, -1, 1, 0), (1, 0, 2, -1)], [[(1, -1, 1, 0)], [(1, -1, 1, 0), (1, 0, 1, -1)], [(1, 1, -1, 0), (1, 0, 2, -1)]]),
    ("pivot-ex4 v1", (1, 0, 1, 0), (1, 1, 0, 1), [(1, -1, -1, 0), (0, 1, 0, -1)], [[(1, -1, -1, 0)], [(1, -1, -1, 0), (1, 0, 0, -1)], [(1, 1, -1, 0), (0, 1, 0, -1)]]),
    ("pivot-ex4 v2", (1, 2, 0, 1), (0, 1, 1, 1), [(2, -1, 1, 0), (1, 0, 1, -1)], [[(2, -1, 1, 0)], [(2, -1, 1, 0), (1, 0, -1, -1)], [(2, 1, -1, 0), (1, 0, 1, -1)]]),
]:
    F = [u1, u2]
    assert meme_vect(noyau(Matrix(bon)), F) and Matrix([u1, u2]).rank() == 2, nom
    for L in leurres:
        assert not meme_vect(noyau(Matrix(L)), F), (nom, L)

# --- pivot-ex5 (et variantes) : valeurs de m sans solution.
for nom, A, B, bon in [
    ("pivot-ex5", mat((1, 1, m), (1, m, 1), (m, 1, 1)), (1, 1, 1), {-2}),
    ("pivot-ex5 v1", mat((1, 1, 1), (1, m, 1), (1, 1, m)), (1, 1, m), set()),
    ("pivot-ex5 v2", mat((1, 1, 1), (1, 2, m), (1, m, 2)), (1, 2, 3), {0, 2}),
]:
    Bm = Matrix(B)
    sans = {r for r in solve(A.det(), m) if A.subs(m, r).rank() < A.subs(m, r).row_join(Bm.subs(m, r)).rank()}
    assert sans == bon, nom
# Leurres : pivot-ex5 {1}, {1, -2}, aucune ; v1 {1}, {-1}, {0, 1} ; v2 {2}, {0}, aucune.
assert {-2} not in ({1}, {1, -2}, set()) and set() not in ({1}, {-1}, {0, 1}) and {0, 2} not in ({2}, {0}, set())
# Solutions explicites.
A, Bm = mat((1, 1, m), (1, m, 1), (m, 1, 1)), Matrix([1, 1, 1])
assert simplify(A * Matrix([1, 1, 1]) / (m + 2) - Bm) == zeros(3, 1)  # pivot-ex5 : x = y = z = 1/(m + 2)
assert mat((1, 1, 1), (1, m, 1), (1, 1, m)) * Matrix([0, 0, 1]) == Matrix([1, 1, m])  # pivot-ex5 v1 : (0, 0, 1)
A, Bm = mat((1, 1, 1), (1, 2, m), (1, m, 2)), Matrix([1, 2, 3])
sol = Matrix([(m - 3) / m, (2 * m - 3) / (m * (m - 2)), (m - 3) / (m * (m - 2))])
assert simplify(A * sol - Bm) == zeros(3, 1)  # pivot-ex5 v2 : solution unique pour m ∉ {0, 2}
assert expand(1 - (m - 1) ** 2 - m * (2 - m)) == 0  # pivot-ex5 v2 : facteur m(2 - m)

print("al_l1_applications : toutes les vérifications passent.")
