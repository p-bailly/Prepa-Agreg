"""Vérifications calculatoires de l'enrichissement du chapitre « al-l1-determinant » (QCM et exercices ajoutés).

Fichiers content/qcm/al-l1-determinant.json et fin de content/exercices/algebre-lineaire/l1-al-l1-determinant.json
(al-l1-determinant-produit-ex5, al-l1-determinant-famille-ex5, al-l1-determinant-inversibilite-ex5). Chaque assert
cite l'identifiant de la question ou de l'exercice (variantes : « v1 », « v2 ») ; les leurres sont vérifiés FAUX.
Lancer : python3 verification/enrichi_al_l1_determinant.py (doit finir sans erreur).
"""

import random
from itertools import permutations

from sympy import Matrix, Rational, diag, eye, factor, solve, sqrt, symbols, zeros

t = symbols("t")

# --- QCM ---

# al-l1-determinant-q1 : det A = 2, A d'ordre 4, det(3A) = 3^4 * 2 = 162 ; leurres 6, 54, 81 faux.
A = diag(2, 1, 1, 1)
assert (3 * A).det() == 162
assert all(v != 162 for v in (6, 54, 81))

# al-l1-determinant-q2 : ordre 3, det A = -2, det B = 3 : det(2 A^{-1} B^2) = 8 * (-1/2) * 9 = -36.
A = diag(-2, 1, 1)
B = diag(3, 1, 1)
assert (2 * A.inv() * B**2).det() == -36
assert all(v != -36 for v in (36, -9, -18))
# Sur des matrices non diagonales aussi.
A2 = Matrix([[1, 2, 0], [0, 1, 3], [1, 0, 0]])  # det = 6 ? on normalise ci-dessous
assert A2.det() == 6
A2 = A2 * diag(Rational(-1, 3), 1, 1)  # det -2
B2 = Matrix([[1, 1, 0], [0, 3, 1], [0, 0, 1]])  # det 3
assert A2.det() == -2 and B2.det() == 3
assert (2 * A2.inv() * B2**2).det() == -36

# al-l1-determinant-q3 : déterminant 4x4 développé selon la deuxième colonne.
M = Matrix([[2, 0, 1, 3], [1, -1, 0, 2], [0, 0, 3, 0], [4, 0, 1, 1]])
assert M.det() == 30
assert M.minor_submatrix(1, 1) == Matrix([[2, 1, 3], [0, 3, 0], [4, 1, 1]])
assert Matrix([[2, 1, 3], [0, 3, 0], [4, 1, 1]]).det() == -30
assert Matrix([[2, 3], [4, 1]]).det() == -10
assert all(v != 30 for v in (-30, 10, -10))

# al-l1-determinant-q4 : Vandermonde V(1, 2, m) = (2-1)(m-1)(m-2) ; inversible ssi m ∉ {1, 2}.
m = symbols("m")
V = Matrix([[1, 1, 1], [1, 2, 4], [1, m, m**2]])
assert factor(V.det()) == factor((m - 1) * (m - 2))
assert set(solve(V.det(), m)) == {1, 2}
# Leurres : « m ∉ {1, 4} » (m = 2 non inversible mais autorisé), « m ≠ 2 » (m = 1 autorisé à tort),
# « m ∉ {0, 1, 2} » (m = 0 donne une matrice inversible, exclue à tort).
assert V.subs(m, 2).det() == 0 and V.subs(m, 4).det() != 0
assert V.subs(m, 1).det() == 0
assert V.subs(m, 0).det() == 2

# al-l1-determinant-q5 (faux) : A^2 = I_n mais det A = -1.
A = diag(-1, 1)
assert A**2 == eye(2) and A.det() == -1
A = -eye(3)
assert A**2 == eye(3) and A.det() == -1

# al-l1-determinant-q6 (vrai) : tA A = I_3 implique (det A)^2 = 1 ; exemples det 1 et det -1.
for A in (Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]), Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])):
    assert A.T * A == eye(3) and A.det() in (1, -1)
assert Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]]).det() == -1

# al-l1-determinant-q7 (faux) : renverser l'ordre des colonnes d'une matrice d'ordre 4 ne change pas le déterminant
# (deux échanges) ; contre-exemple I_4 -> matrice antidiagonale de déterminant +1. À l'ordre 3, c'est -det.
J4 = Matrix(4, 4, lambda i, j: 1 if i + j == 3 else 0)
assert J4.det() == 1
random.seed(1)
for _ in range(5):
    R = Matrix(4, 4, lambda i, j: random.randint(-5, 5))
    Rrev = R[:, ::-1]
    assert Rrev.det() == R.det()
    R3 = Matrix(3, 3, lambda i, j: random.randint(-5, 5))
    assert R3[:, ::-1].det() == -R3.det()

# al-l1-determinant-q8 : det(u, v, w) = 5 ; det(w, u, v + 2u) = det(w, u, v) = det(u, v, w) = 5 (permutation
# circulaire, paire). Leurres -5, 10, 15 faux.
u, v, w = Matrix([1, 2, 0]), Matrix([0, 1, 1]), Matrix([2, 0, 1])
P = Matrix.hstack(u, v, w)
assert P.det() != 0
facteur = Rational(5) / P.det()
u = u * facteur  # det(u, v, w) = 5
assert Matrix.hstack(u, v, w).det() == 5
assert Matrix.hstack(w, u, v + 2 * u).det() == 5
assert all(x != 5 for x in (-5, 10, 15))

# --- Exercices ---


def blocs(A, Bm, D):
    p, q = A.rows, D.rows
    M = zeros(p + q, p + q)
    M[:p, :p] = A
    M[:p, p:] = Bm
    M[p:, p:] = D
    return M


# al-l1-determinant-produit-ex5 : factorisation M = diag(I_p, D) * [[A, B], [0, I_q]] et det M = det A det D,
# testée sur des matrices aléatoires (p, q de 1 à 3), puis les applications numériques.
random.seed(2)
for p in range(1, 4):
    for q in range(1, 4):
        A = Matrix(p, p, lambda i, j: random.randint(-4, 4))
        Bm = Matrix(p, q, lambda i, j: random.randint(-4, 4))
        D = Matrix(q, q, lambda i, j: random.randint(-4, 4))
        G = blocs(eye(p), zeros(p, q), D)
        H = blocs(A, Bm, eye(q))
        assert G * H == blocs(A, Bm, D)
        assert G.det() == D.det() and H.det() == A.det()
        assert blocs(A, Bm, D).det() == A.det() * D.det()
        # Récurrence de la question 1 : développer selon la première colonne enlève un bloc identité.
        if p >= 2:
            assert G.det() == blocs(eye(p - 1), zeros(p - 1, q), D).det()
        # Question 2 : développer H selon sa dernière ligne enlève la dernière colonne de B.
        if q >= 2:
            assert H.det() == blocs(A, Bm[:, :q - 1], eye(q - 1)).det()

# Application (exercice) : det = (2*3 - 1*1)(4*3 - 1*2) = 5 * 10 = 50.
M = Matrix([[2, 1, 5, -3], [1, 3, 0, 7], [0, 0, 4, 1], [0, 0, 2, 3]])
assert M[:2, :2].det() == 5 and M[2:, 2:].det() == 10 and M.det() == 50
# Leurres : 5 + 10 = 15, -50, 4*... ; on vérifie seulement qu'ils diffèrent de la bonne valeur.
assert all(x != 50 for x in (15, -50, 10))
# v1 : det = 5 * (-7) = -35.
M = Matrix([[3, -1, 2, 4], [2, 1, -5, 1], [0, 0, 1, 2], [0, 0, 3, -1]])
assert M[:2, :2].det() == 5 and M[2:, 2:].det() == -7 and M.det() == -35
assert all(x != -35 for x in (35, -2, -7))
# v2 : det = (-5) * (-1) = 5.
M = Matrix([[1, 4, 0, 2], [2, 3, 6, -1], [0, 0, 5, 2], [0, 0, 3, 1]])
assert M[:2, :2].det() == -5 and M[2:, 2:].det() == -1 and M.det() == 5
assert all(x != 5 for x in (-5, -6, 6))
# Question 4 (exercice et variantes) : M inversible ssi A et D le sont ; un bloc singulier rend M singulière.
assert blocs(Matrix([[1, 2], [2, 4]]), Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])).det() == 0


# al-l1-determinant-famille-ex5 : diagonale impaire, hors diagonale paire => det impair (donc non nul).
def leibniz_termes(M):
    n = M.rows
    for s in permutations(range(n)):
        prod = 1
        for j in range(n):
            prod *= M[s[j], j]
        yield s, prod


random.seed(3)
for n in range(2, 6):
    for _ in range(20):
        M = Matrix(n, n, lambda i, j: 2 * random.randint(-3, 3) + 1 if i == j else 2 * random.randint(-3, 3))
        for s, prod in leibniz_termes(M):
            if list(s) != list(range(n)):
                assert prod % 2 == 0  # tout terme avec σ ≠ id contient un facteur hors diagonale
        assert M.det() % 2 == 1
# Leurres de la réponse : « pair » et « nul » faux (déterminant impair), « vaut ±1 » faux : diag(3, 1) donne 3.
assert diag(3, 1).det() == 3

# al-l1-determinant-inversibilite-ex5 : base ssi det != 0.
U = Matrix([[t, 1, 0], [1, t, 1], [0, 1, t]])
assert factor(U.det()) == t * (t**2 - 2)
assert set(solve(U.det(), t)) == {0, sqrt(2), -sqrt(2)}
# Leurres : t ∉ {-√2, √2} (oublie 0), t ≠ 0, t ∉ {-1, 0, 1} (t = 1 donne une base).
assert U.subs(t, 1).det() != 0 and U.subs(t, 0).det() == 0
# v1 :
U = Matrix([[t, 2, 0], [2, t, 2], [0, 2, t]])
assert factor(U.det()) == t * (t**2 - 8)
assert set(solve(U.det(), t)) == {0, 2 * sqrt(2), -2 * sqrt(2)}
assert U.subs(t, 2).det() != 0  # leurre t ∉ {-2, 0, 2}
# v2 :
U = Matrix([[1, t, 0], [t, 1, t], [0, t, 1]])
assert (U.det() - (1 - 2 * t**2)).expand() == 0
assert set(solve(U.det(), t)) == {sqrt(2) / 2, -sqrt(2) / 2}
assert U.subs(t, 1).det() != 0 and U.subs(t, 0).det() != 0  # leurres t ∉ {-1, 1}, t ≠ 0 (exclusions à tort)

print("enrichi_al_l1_determinant : tout est vérifié.")
