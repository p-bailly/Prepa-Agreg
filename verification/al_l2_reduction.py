"""Vérifications calculatoires du chapitre « al-l2-reduction » (algèbre linéaire, L2) :
cartes, exercices et QCM de content/*/algebre-lineaire/l2-al-l2-reduction.json et content/qcm/al-l2-reduction.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/al_l2_reduction.py (doit finir sans erreur, en quelques secondes).
"""

from fractions import Fraction as F
from itertools import product

from sympy import (Matrix, Rational, binomial, diag, diff, exp, eye, expand, factor, simplify, symbols, zeros, I as i_,
                   Poly, sqrt)

X, t, k = symbols("X t k")
I2, I3 = eye(2), eye(3)


def chi(A):
    return A.charpoly(X).as_expr()


def evalp(p, A):
    """p(A) pour un polynôme p en X (Horner)."""
    n = A.shape[0]
    R = zeros(n)
    for c in Poly(expand(p), X).all_coeffs():
        R = R * A + c * eye(n)
    return R


def dimE(A, lam):
    return len((A - lam * eye(A.shape[0])).nullspace())


def diagonalisable_R(A):
    return A.is_diagonalizable(reals_only=True)


R = Matrix([[0, -1], [1, 0]])          # rotation d'angle pi/2
J = Matrix([[1, 1], [0, 1]])           # bloc de Jordan
M3 = Matrix([[1, 0, 0], [0, 0, -1], [0, 1, 0]])  # (1) ⊕ rotation, non trigonalisable sur R

# ---------------------------------------------------------------------------
# al-l2-reduction-valeurs-propres (carte)
# ---------------------------------------------------------------------------
# Exemple : A = [[2,1],[1,2]], vecteurs propres (1,1) (3) et (1,-1) (1).
A = Matrix([[2, 1], [1, 2]])
assert A * Matrix([1, 1]) == 3 * Matrix([1, 1]) and A * Matrix([1, -1]) == Matrix([1, -1])
# Contre-exemple / affirmation 1 : la rotation n'a aucune valeur propre réelle.
assert factor(chi(R)) == X**2 + 1 and not any(r.is_real for r in R.eigenvals())
# Piège 2 / contre-exemple : diag(1,2), (e1, 2e1, e2) liée ; piège 3 : id, (e1, e2, e1+e2) liée.
assert Matrix([[1, 2, 0], [0, 0, 1]]).rank() == 2 and Matrix([[1, 0, 1], [0, 1, 1]]).rank() == 2
# Trou 2, leurre « orthogonale » : [[1,1],[0,2]] a des vecteurs propres (1,0), (1,1) non orthogonaux.
B = Matrix([[1, 1], [0, 2]])
assert B * Matrix([1, 0]) == Matrix([1, 0]) and B * Matrix([1, 1]) == 2 * Matrix([1, 1])
# Trou 3, leurre det u = lambda : diag(1,2), 1 est valeur propre mais det = 2.
assert diag(1, 2).det() == 2
# Affirmation 2 : u(e1+e2) = (1,2) non colinéaire à (1,1).
assert Matrix([[1, 1], [1, 2]]).det() != 0

# valeurs-propres-ex1 (+ v1, v2) : x=(1,2), y=(1,1), matrice dans (x,y).
x, y = Matrix([1, 2]), Matrix([1, 1])
for A, l1, l2, leurres in [
    (Matrix([[5, -2], [4, -1]]), 1, 3, [diag(3, 1), Matrix([[1, 1], [0, 3]]), diag(5, -1)]),
    (Matrix([[4, -1], [2, 1]]), 2, 3, [diag(3, 2), Matrix([[2, 1], [0, 3]]), diag(4, 1)]),
    (Matrix([[-1, 2], [-4, 5]]), 3, 1, [diag(1, 3), Matrix([[3, 1], [0, 1]]), diag(-1, 5)]),
]:
    assert A * x == l1 * x and A * y == l2 * y
    P = Matrix.hstack(x, y)
    D = P.inv() * A * P
    assert D == diag(l1, l2) and all(D != L for L in leurres)

# valeurs-propres-ex2 (+ v1, v2) : phi(X^k) = c_k X^k sur R_n[X].
for phi, n, vals in [(lambda p: X * diff(p, X), 3, [0, 1, 2, 3]),
                     (lambda p: 2 * X * diff(p, X) - p, 3, [-1, 1, 3, 5]),
                     (lambda p: X * diff(p, X) + p, 2, [1, 2, 3])]:
    for kk in range(n + 1):
        assert expand(phi(X**kk) - vals[kk] * X**kk) == 0
    assert len(set(vals)) == n + 1  # n+1 valeurs propres distinctes = dimension

# valeurs-propres-ex4 (+ v1, v2) : Vandermonde, nombre de sous-espaces stables.
assert Matrix([[1, 1, 1], [1, 2, 3], [1, 4, 9]]).det() == 2
assert Matrix(4, 4, lambda a, b: (b + 1) ** a).det() == 12
assert 2**3 == 8 and 2**4 == 16


def stables_diag(d):
    """Nombre de sous-espaces engendrés par des vecteurs de la base canonique (tous stables, vérifié ici) ;
    l'exercice montre que ce sont les seuls quand les coefficients diagonaux sont distincts."""
    n = len(d)
    A = diag(*d)
    nb = 0
    for parts in product([0, 1], repeat=n):
        cols = [eye(n)[:, j] for j in range(n) if parts[j]]
        if cols:
            F_ = Matrix.hstack(*cols)
            assert Matrix.hstack(F_, A * F_).rank() == F_.rank()  # stable
        nb += 1
    return nb


assert stables_diag([1, 2, 3]) == 8 and stables_diag([1, 2, 3, 4]) == 16
# v2 : diag(1,1,2), une infinité de droites stables Vect(e1 + t e2).
A = diag(1, 1, 2)
for tt in range(-3, 4):
    v = Matrix([1, tt, 0])
    assert A * v == v
# Une droite non engendrée par des e_i, non stable pour diag(1,2,3) : Vect(e1+e2).
assert Matrix.hstack(Matrix([1, 1, 0]), diag(1, 2, 3) * Matrix([1, 1, 0])).rank() == 2

# ---------------------------------------------------------------------------
# al-l2-reduction-polynome-caracteristique (carte)
# ---------------------------------------------------------------------------
assert expand(chi(Matrix([[4, -2], [1, 1]])) - (X - 2) * (X - 3)) == 0                      # exemple
B = Matrix([[1, 2, 3], [0, -1, 4], [0, 0, 1]])
assert expand(chi(B) - (X - 1) ** 2 * (X + 1)) == 0 and B.trace() == 1 and B.det() == -1  # exemple
assert expand(chi(Matrix([[0, 1], [1, 0]])) - (X**2 - 1)) == 0                             # contre-exemple
assert chi(I2) == chi(J) and J != I2                                                       # contre-exemple, affirmation 2
assert diag(4, 1).det() == 4 and 2 not in diag(4, 1).eigenvals()                           # piège det u = lambda^n
# Affirmation 4 : A, B nilpotentes, A + B de valeurs propres ±1.
Aa, Bb = Matrix([[0, 1], [0, 0]]), Matrix([[0, 0], [1, 0]])
assert set((Aa + Bb).eigenvals()) == {1, -1} and set(Aa.eigenvals()) == {0} and set(Bb.eigenvals()) == {0}
# Trous : leurres réfutés sur I_2 (coefficient -2, terme constant 1), sur (1) et (2) en taille 1.
c = Poly(chi(I2), X).all_coeffs()
assert c == [1, -2, 1] and c[1] not in (I2.trace(), I2.det(), -I2.det()) and c[2] not in (-I2.det(), I2.trace())
assert Poly(chi(Matrix([[2]])), X).all_coeffs() == [1, -2]          # coefficient dominant 1 ≠ det = 2, ≠ (-1)^1
assert Poly(chi(Matrix([[1]])), X).all_coeffs()[1] == -1            # terme constant -1 ≠ det = 1
# Formules générales : coefficients de chi en taille 3 symbolique.
a = symbols("a0:9")
A = Matrix(3, 3, a)
cc = Poly(chi(A), X).all_coeffs()
assert cc[0] == 1 and expand(cc[1] + A.trace()) == 0 and expand(cc[3] + A.det()) == 0

# polynome-caracteristique-ex1 (+ v1, v2).
for A, vals, v1, v2, leurres in [
    (Matrix([[1, 2], [3, 0]]), {3, -2}, (3, [1, 1]), (-2, [2, -3]), [{2, -3}, {1, 0}, {-1, 6}]),
    (Matrix([[3, 4], [1, 0]]), {4, -1}, (4, [4, 1]), (-1, [1, -1]), [{1, -4}, {3, 0}, {-2, 2}]),
    (Matrix([[0, 1], [-6, 5]]), {2, 3}, (2, [1, 2]), (3, [1, 3]), [{-2, -3}, {0, 5}, {1, 6}]),
]:
    assert set(A.eigenvals()) == vals and all(set(A.eigenvals()) != L for L in leurres)
    for lam, v in (v1, v2):
        assert A * Matrix(v) == lam * Matrix(v)

# polynome-caracteristique-ex2 (+ v1, v2).
for A, ch, dims, vecs, leurres in [
    (Matrix([[3, 1, -1], [1, 3, -1], [0, 0, 2]]), (X - 2) ** 2 * (X - 4), {2: 2, 4: 1}, {4: [1, 1, 0]},
     [(X - 2) * (X - 4) ** 2, (X - 2) * (X - 3) ** 2, (X - 1) * (X - 2) * (X - 4)]),
    (Matrix([[1, 2, 0], [2, 1, 0], [1, 1, 3]]), (X - 3) ** 2 * (X + 1), {3: 1, -1: 1}, {3: [0, 0, 1], -1: [1, -1, 0]},
     [(X - 3) * (X + 1) ** 2, (X + 3) ** 2 * (X - 1), (X - 1) * (X - 3) * (X + 1)]),
    (Matrix([[4, 1, 0], [2, 3, 0], [-1, 2, 5]]), (X - 5) ** 2 * (X - 2), {5: 1, 2: 1}, {5: [0, 0, 1], 2: [3, -6, 5]},
     [(X - 5) * (X - 2) ** 2, (X - 3) * (X - 4) * (X - 5), (X + 5) ** 2 * (X + 2)]),
]:
    assert expand(chi(A) - ch) == 0 and all(expand(chi(A) - L) != 0 for L in leurres)
    assert all(dimE(A, l) == d for l, d in dims.items())
    for l, v in vecs.items():
        assert A * Matrix(v) == l * Matrix(v)
# ex2 : plan E_2 d'équation x + y - z = 0.
assert Matrix([[3, 1, -1], [1, 3, -1], [0, 0, 2]]) * Matrix([1, 0, 1]) == 2 * Matrix([1, 0, 1])

# polynome-caracteristique-ex3 : chi_J = X^(n-1)(X-n) ; leurres faux ; chi(aI+bJ).
a_, b_ = symbols("a b")
for n in range(2, 7):
    Jn = Matrix.ones(n, n)
    assert expand(chi(Jn) - X ** (n - 1) * (X - n)) == 0
    for L in [X ** (n - 1) * (X - 1), X * (X - 1) ** (n - 1), X**n - n]:
        assert expand(chi(Jn) - L) != 0
    if n <= 4:
        assert expand((a_ * eye(n) + b_ * Jn).charpoly(X).as_expr() - (X - a_) ** (n - 1) * (X - a_ - n * b_)) == 0

# polynome-caracteristique-ex4 : AB et BA (question 4), et chi_AB = chi_BA sur des exemples non inversibles.
A_, B_ = Matrix([[1, 0], [0, 0]]), Matrix([[0, 1], [0, 0]])
assert A_ * B_ == Matrix([[0, 1], [0, 0]]) and B_ * A_ == zeros(2) and chi(A_ * B_) == chi(B_ * A_) == X**2
A_ = Matrix([[1, 2, 3], [2, 4, 6], [0, 1, 1]])
B_ = Matrix([[0, 1, 0], [5, -1, 2], [3, 3, 3]])
assert A_.det() == 0 and expand(chi(A_ * B_) - chi(B_ * A_)) == 0

# ---------------------------------------------------------------------------
# al-l2-reduction-diagonalisabilite (carte)
# ---------------------------------------------------------------------------
A = Matrix([[1, 0, 1], [0, 1, 0], [0, 0, 2]])
B = Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 2]])
assert chi(A) == chi(B) and expand(chi(A) - (X - 1) ** 2 * (X - 2)) == 0                 # exemple
P = Matrix([[1, 0, 1], [0, 1, 0], [0, 0, 1]])
assert P.inv() * A * P == diag(1, 1, 2) and dimE(B, 1) == 1 and not B.is_diagonalizable()
assert dimE(J, 1) == 1 and not J.is_diagonalizable()                                      # contre-exemple, pièges 2-3
assert expand(chi(M3) - (X - 1) * (X**2 + 1)) == 0 and dimE(M3, 1) == 1 and not diagonalisable_R(M3)  # piège 1
assert not R.is_diagonalizable(reals_only=True) and R.is_diagonalizable()                 # dépend du corps
N_ = Matrix([[0, 1], [0, 0]])
assert N_**2 == zeros(2) and not N_.is_diagonalizable()                                   # affirmation 4
# Trou 1, leurre tr u ≠ 0 / u inversible : la matrice nulle est diagonalisable.
assert zeros(2).is_diagonalizable()

# diagonalisabilite-ex1 (+ v1, v2) : P commune, D = diag(lam, mu, mu).
P = Matrix([[1, 1, 1], [1, -1, 0], [1, 0, -1]])
for c_, lam, mu in [(0, 2, -1), (1, 3, 0), (2, 4, 1)]:
    A = Matrix.ones(3, 3) + (c_ - 1) * I3
    assert expand(chi(A) - (X - lam) * (X - mu) ** 2) == 0
    assert P.inv() * A * P == diag(lam, mu, mu) and (A - mu * I3).rank() == 1
    s_ = X - c_  # coefficient diagonal de XI - A ; après opérations : X - mu sur la diagonale
    assert expand(s_ + 1 - (X - mu)) == 0
    M = Matrix([[1, -1, -1], [1, s_, -1], [1, -1, s_]])
    assert expand((X - lam) * M.det() - chi(A)) == 0

# diagonalisabilite-ex2 (+ v1, v2) : conditions sur les paramètres (balayage).
for a1, b1 in product(range(-3, 4), repeat=2):
    A = Matrix([[2, 1, a1], [0, 1, b1], [0, 0, 2]])
    assert A.is_diagonalizable() == (a1 + b1 == 0)
    A = Matrix([[2, 0, a1], [b1, 1, 0], [0, 0, 2]])
    assert expand(chi(A) - (X - 2) ** 2 * (X - 1)) == 0 and A.is_diagonalizable() == (a1 == 0)
for a1, b1, c1 in product(range(-1, 2), repeat=3):
    assert Matrix([[1, a1, b1], [0, 1, c1], [0, 0, 1]]).is_diagonalizable() == (a1 == b1 == c1 == 0)
# ex2 : leurres réfutés ((1,1) : a = b mais non diagonalisable ; (1,-1) diagonalisable).
assert not Matrix([[2, 1, 1], [0, 1, 1], [0, 0, 2]]).is_diagonalizable()
assert Matrix([[2, 1, 1], [0, 1, -1], [0, 0, 2]]).is_diagonalizable()
A = Matrix([[2, 1, 1], [0, 1, -1], [0, 0, 2]])
assert A * Matrix([0, 1, -1]) == 2 * Matrix([0, 1, -1]) and A * Matrix([1, -1, 0]) == Matrix([1, -1, 0])

# diagonalisabilite-ex3 : rang 1, diagonalisable ssi trace ≠ 0 (exemples et chi).
for M in [Matrix([[1, 1], [1, 1]]), Matrix([[1, 1], [-1, -1]]), Matrix([[1, 2, 3], [2, 4, 6], [-1, -2, -3]]),
          Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])]:
    assert M.rank() == 1
    n = M.shape[0]
    assert expand(chi(M) - X ** (n - 1) * (X - M.trace())) == 0
    assert M.is_diagonalizable() == (M.trace() != 0)
assert Matrix([[1, 1], [-1, -1]]) ** 2 == zeros(2)

# diagonalisabilite-ex4 : transposition, tr T = n, det T = (-1)^(n(n-1)/2), leurres faux.
for n in range(2, 5):
    T = zeros(n * n)
    for p_, q_ in product(range(n), repeat=2):
        T[q_ * n + p_, p_ * n + q_] = 1  # E_pq -> E_qp (coordonnées en ordre ligne)
    assert T * T == eye(n * n) and T.trace() == n and T.det() == (-1) ** (n * (n - 1) // 2)
    assert T.is_diagonalizable()
assert (-1) ** (2 * 1 // 2) != (-1) ** 2  # leurre (-1)^n faux pour n = 2

# ---------------------------------------------------------------------------
# al-l2-reduction-lemme-noyaux (carte)
# ---------------------------------------------------------------------------
A = Matrix([[1, 2], [0, -1]])
assert A**2 == I2 and (A + I2) / 2 == Matrix([[1, 1], [0, 0]])                          # exemple
assert (A - I2) * Matrix([1, 0]) == zeros(2, 1) and (A + I2) * Matrix([1, -1]) == zeros(2, 1)
assert expand(Rational(1, 2) * (X + 1) - Rational(1, 2) * (X - 1) - 1) == 0
assert (J - I2) ** 2 == zeros(2) and J != I2                                              # contre-exemple
assert R**2 + I2 == zeros(2)                                                              # contre-exemple
# Affirmation 4 : A = [[0,-1],[1,-1]], A^2, A^3 = I, chi = X^2 + X + 1 sans racine réelle.
A = Matrix([[0, -1], [1, -1]])
assert A**2 == Matrix([[-1, 1], [-1, 0]]) and A**3 == I2 and expand(chi(A) - (X**2 + X + 1)) == 0
assert not diagonalisable_R(A)

# lemme-noyaux-ex1 (+ v1, v2) : inverses, leurres faux (testés sur une matrice diagonale admissible).
for poly, r1, r2, inv, leurres in [
    (X**2 - 3 * X + 2, 1, 2, lambda M: (3 * eye(2) - M) / 2,
     [lambda M: (M - 3 * eye(2)) / 2, lambda M: 3 * eye(2) - M, lambda M: (2 * eye(2) - M) / 3]),
    (X**2 - X - 6, 3, -2, lambda M: (M - eye(2)) / 6,
     [lambda M: (eye(2) - M) / 6, lambda M: M - eye(2), lambda M: (M + eye(2)) / 6]),
    (X**2 + 2 * X - 3, 1, -3, lambda M: (M + 2 * eye(2)) / 3,
     [lambda M: (M - 2 * eye(2)) / 3, lambda M: M + 2 * eye(2), lambda M: (M + 3 * eye(2)) / 2]),
]:
    assert expand(poly - (X - r1) * (X - r2)) == 0
    for M in [diag(r1, r2), diag(r1, r1, r2), Matrix([[r1, 1], [0, r2]])]:
        assert evalp(poly, M) == zeros(M.shape[0])
        n = M.shape[0]
        Mi = inv(M) if n == 2 else None
        if n == 2:
            assert M * Mi == eye(2)
            assert all(M * L(M) != eye(2) for L in leurres)

# lemme-noyaux-ex2 (+ v1, v2) : polynômes minimaux.
A = Matrix([[3, -1, 1], [2, 0, 2], [1, -1, 3]])
N_ = A - 2 * I3
assert expand(chi(A) - (X - 2) ** 3) == 0 and N_ != zeros(3) and N_**2 == zeros(3)
A = Matrix([[2, 0, 1], [1, 1, 1], [0, 0, 1]])
assert expand(chi(A) - (X - 1) ** 2 * (X - 2)) == 0 and (A - I3) * (A - 2 * I3) == zeros(3) and A.is_diagonalizable()
A = Matrix([[2, 1, 1], [0, 1, 0], [0, -1, 2]])
M_ = (A - 2 * I3) * (A - I3)
assert expand(chi(A) - (X - 2) ** 2 * (X - 1)) == 0 and M_ == Matrix([[0, -1, 1], [0, 0, 0], [0, 0, 0]])
assert (A - 2 * I3) * M_ == zeros(3) and not A.is_diagonalizable()
assert (A - 2 * I3) ** 2 != zeros(3)  # leurre (X-2)^2 n'annule pas A

# lemme-noyaux-ex3 : projecteurs id - u^2, (u^2+u)/2, (u^2-u)/2 sur un exemple u^3 = u.
u = diag(0, 1, -1, 1)
for Q, cible in [(1 - X**2, 0), ((X**2 + X) / 2, 1), ((X**2 - X) / 2, -1)]:
    pq = evalp(Q, u)
    assert pq == diag(*[1 if v == cible else 0 for v in (0, 1, -1, 1)])
# leurres : u^2 et id - u ne sont pas le projecteur sur Ker u.
assert u**2 != diag(1, 0, 0, 0) and eye(4) - u != diag(1, 0, 0, 0)

# lemme-noyaux-ex4 : contre-exemple de la question 3.
A_, B_ = Matrix([[1, 1], [0, 0]]), Matrix([[-1, 0], [0, 0]])
assert A_.is_diagonalizable() and (A_ + B_) == Matrix([[0, 1], [0, 0]]) and not (A_ + B_).is_diagonalizable()
assert A_ * B_ == Matrix([[-1, 0], [0, 0]]) and B_ * A_ == Matrix([[-1, -1], [0, 0]])

# lemme-noyaux-ex5 (+ v1, v2) : exemples et rangs.
A = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
assert A**3 + A == zeros(3) and A.rank() == 2
A = Matrix([[0, -2, 0], [2, 0, 0], [0, 0, 0]])
assert A**2 == diag(-4, -4, 0) and A**3 + 4 * A == zeros(3) and A.rank() == 2
A = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
assert A**3 - A**2 + A - I3 == zeros(3) and (A - I3).rank() == 2
assert expand(X**3 - X**2 + X - 1 - (X - 1) * (X**2 + 1)) == 0
# Parité : (det B)^2 = det(-c I_d) = (-c)^d est négatif pour d impair (c = 1 ou 4).
assert all((-1) ** d < 0 and (-4) ** d < 0 for d in (1, 3))

# ---------------------------------------------------------------------------
# al-l2-reduction-cayley-hamilton (carte)
# ---------------------------------------------------------------------------
A = Matrix([[1, 2], [3, 4]])
assert A**2 - 5 * A - 2 * I2 == zeros(2) and A.inv() == (A - 5 * I2) / 2                 # exemple
assert A.inv() == Matrix([[-2, 1], [Rational(3, 2), Rational(-1, 2)]])
D_ = diag(1, 0)
assert 2 * D_ - D_.trace() * I2 == diag(1, -1)                                            # fausse preuve
assert evalp((X - 1) * (X - 2), I2) == zeros(2) and chi(I2) != chi(diag(1, 2))            # chi_A(B) = 0
# Trous : chi(u) = id faux pour u = 0 ; chi(0) ≠ 0 pour id ; chi(tr u) ≠ 0 pour id_2.
assert chi(Matrix([[0]])) == X and chi(I2).subs(X, 2) == 1 and chi(I2).subs(X, 0) == 1
# Cayley-Hamilton vérifié symboliquement en taille 3, et compagnons (carte, démonstration) jusqu'à p = 6.
A = Matrix(3, 3, a)
assert evalp(chi(A), A).applyfunc(expand) == zeros(3)
for p in range(1, 7):
    cs = symbols(f"c0:{p}")
    C = zeros(p)
    for j in range(p - 1):
        C[j + 1, j] = 1
    for j in range(p):
        C[j, p - 1] = -cs[j]
    Pc = X**p + sum(cs[j] * X**j for j in range(p))
    assert expand(chi(C) - Pc) == 0
    if p <= 4:
        assert evalp(Pc, C).applyfunc(expand) == zeros(p)
# Affirmation : A nilpotente complexe => A^n = 0 (exemple d'indice 3 en taille 3).
N_ = Matrix([[0, 1, 2], [0, 0, 3], [0, 0, 0]])
assert N_**3 == zeros(3) and chi(N_) == X**3

# cayley-hamilton-ex1 (+ v1, v2) : A^3 = alpha A + beta I.
for A, al, be, leurres in [
    (Matrix([[2, 1], [1, 3]]), 20, -25, [(25, -20), (20, 25), (5, -5)]),
    (Matrix([[1, 3], [2, 2]]), 13, 12, [(12, 13), (13, -12), (3, 4)]),
    (Matrix([[3, 1], [2, 2]]), 21, -20, [(20, -21), (21, 20), (5, -4)]),
]:
    assert A**3 == al * A + be * I2 and all(A**3 != p_ * A + q_ * I2 for p_, q_ in leurres)
assert Matrix([[2, 1], [1, 3]]) ** 2 == Matrix([[5, 5], [5, 10]]) and Matrix([[2, 1], [1, 3]]) ** 3 == Matrix([[15, 20], [20, 35]])
assert Matrix([[1, 3], [2, 2]]) ** 2 == Matrix([[7, 9], [6, 10]]) and Matrix([[3, 1], [2, 2]]) ** 2 == Matrix([[11, 5], [10, 6]])

# cayley-hamilton-ex2 (+ v1, v2) : inverse par Cayley-Hamilton.
for A, ch, cst, q, A2, invm in [
    (Matrix([[2, 1, 0], [0, 1, 1], [1, 0, 1]]), X**3 - 4 * X**2 + 5 * X - 3, 3, X**2 - 4 * X + 5,
     Matrix([[4, 3, 1], [1, 1, 2], [3, 1, 1]]), Matrix([[1, -1, 1], [1, 2, -2], [-1, 1, 2]]) / 3),
    (Matrix([[0, 1, 0], [0, 0, 1], [2, -5, 4]]), X**3 - 4 * X**2 + 5 * X - 2, 2, X**2 - 4 * X + 5,
     Matrix([[0, 0, 1], [2, -5, 4], [8, -18, 11]]), Matrix([[5, -4, 1], [2, 0, 0], [0, 2, 0]]) / 2),
    (Matrix([[2, 1, 1], [0, 1, 0], [0, -1, 2]]), X**3 - 5 * X**2 + 8 * X - 4, 4, X**2 - 5 * X + 8,
     Matrix([[4, 2, 4], [0, 1, 0], [0, -3, 4]]), Matrix([[2, -3, -1], [0, 4, 0], [0, 2, 2]]) / 4),
]:
    assert expand(chi(A) - ch) == 0 and A**2 == A2 and A.inv() == invm == evalp(q, A) / cst
    assert A.trace() == -Poly(ch, X).all_coeffs()[1] and A.det() == cst
    # leurres : (A^2 + c1 A - c0 I)/cst, cst q(A), q(A) ne sont pas l'inverse.
    qq = Poly(q, X).all_coeffs()
    for L in [evalp(X**2 - qq[1] * X - qq[2], A) / cst, cst * evalp(q, A), evalp(q, A)]:
        assert L != A.inv()
assert expand((X - 2) * (X - 1) ** 2 - 1 - (X**3 - 4 * X**2 + 5 * X - 3)) == 0

# cayley-hamilton-ex3 : A^2 = (tr A) A - (det A) I en taille 2 ; exemple ; leurres de la question 3.
A = Matrix(2, 2, a[:4])
assert (A**2 - A.trace() * A + A.det() * I2).applyfunc(expand) == zeros(2)
A = Matrix([[1, 2], [3, -1]])
assert A.det() == -7 and A**2 == 7 * I2
assert diag(1, 0) ** 2 != -diag(1, 0).det() * I2 and I2**2 != -I2.det() * I2

# ---------------------------------------------------------------------------
# al-l2-reduction-trigonalisation (carte)
# ---------------------------------------------------------------------------
A = Matrix([[3, 1], [-1, 1]])
P = Matrix([[1, 0], [-1, 1]])
assert expand(chi(A) - (X - 2) ** 2) == 0 and P.inv() * A * P == Matrix([[2, 1], [0, 2]])  # exemple
assert A * Matrix([0, 1]) == Matrix([1, 1])
# Trou 3 : leurres sur diag(2,3) et (2).
assert diag(2, 3).det() == 6 and 6 not in (5, 13) and Matrix([[2]]).det() == 2

# trigonalisation-ex1 (+ v1, v2).
for A, lam, v1, T, leurres in [
    (Matrix([[1, 1], [-1, 3]]), 2, [1, 1], Matrix([[2, 1], [0, 2]]), [Matrix([[2, 0], [0, 2]]), Matrix([[2, 2], [0, 2]]), Matrix([[1, 1], [0, 3]])]),
    (Matrix([[5, -1], [4, 1]]), 3, [1, 2], Matrix([[3, -1], [0, 3]]), [Matrix([[3, 1], [0, 3]]), Matrix([[3, 0], [0, 3]]), Matrix([[5, -1], [0, 1]])]),
    (Matrix([[-1, 1], [-4, 3]]), 1, [1, 2], Matrix([[1, 1], [0, 1]]), [Matrix([[1, -1], [0, 1]]), Matrix([[1, 2], [0, 1]]), Matrix([[-1, 1], [0, 3]])]),
]:
    assert expand(chi(A) - (X - lam) ** 2) == 0 and not A.is_diagonalizable()
    assert A * Matrix(v1) == lam * Matrix(v1)
    P = Matrix([[v1[0], 0], [v1[1], 1]])
    assert P.inv() * A * P == T and all(L != T for L in leurres)

# trigonalisation-ex2 (+ v1, v2) : matrices compagnons (diagonalisables ou non) des polynômes donnés.
for ch, Q, val, autres in [((X - 1) ** 2 * (X + 2), X**2 + X + 1, 27, (9, 3, 12)),
                           ((X - 2) * (X + 1) ** 2, X**2 + 1, 20, (10, 5, 8)),
                           (X * (X - 1) * (X - 3), X + 2, 30, (0, 8, 15))]:
    cs = Poly(expand(ch), X).all_coeffs()[::-1]
    C = zeros(3)
    C[1, 0] = C[2, 1] = 1
    for j in range(3):
        C[j, 2] = -cs[j]
    assert expand(chi(C) - ch) == 0 and evalp(Q, C).det() == val and val not in autres
# Valeurs de tr, det, tr(A^2) annoncées.
assert [Matrix([[0, 0, 2], [1, 0, -3], [0, 1, 0]]).trace()] == [0]
for ch, trv, detv, tr2 in [((X - 1) ** 2 * (X + 2), 0, -2, 6), ((X - 2) * (X + 1) ** 2, 0, 2, 6), (X * (X - 1) * (X - 3), 4, 0, 10)]:
    cs = Poly(expand(ch), X).all_coeffs()[::-1]
    C = zeros(3)
    C[1, 0] = C[2, 1] = 1
    for j in range(3):
        C[j, 2] = -cs[j]
    assert C.trace() == trv and C.det() == detv and (C**2).trace() == tr2

# trigonalisation-ex4 : J + eps diag(1,2) diagonalisable ; discriminant -4 en R.
eps = Rational(1, 10)
assert (J + eps * diag(1, 2)).is_diagonalizable()
assert R.trace() ** 2 - 4 * R.det() == -4

# ---------------------------------------------------------------------------
# al-l2-reduction-sous-espaces-caracteristiques (carte)
# ---------------------------------------------------------------------------
A = Matrix([[1, 1, 1], [0, 1, 1], [0, 0, 2]])
assert expand(chi(A) - (X - 1) ** 2 * (X - 2)) == 0 and (A - I3) ** 2 == Matrix([[0, 0, 2], [0, 0, 1], [0, 0, 1]])
assert dimE(A, 1) == 1 and A * Matrix([2, 1, 1]) == 2 * Matrix([2, 1, 1])
P = Matrix([[1, 0, 2], [0, 1, 1], [0, 0, 1]])
assert P.inv() * A * P == Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 2]])
# Trou 1, leurres pour J : Ker(J - I) droite, Im (J-I)^2 = 0, Ker(J^2 - I) droite.
assert dimE(J, 1) == 1 and (J - I2) ** 2 == zeros(2) and len((J**2 - I2).nullspace()) == 1

# sous-espaces-caracteristiques-ex1 (+ v1, v2) : bases des sous-espaces caractéristiques ; leurres.
for A, lam, base, leurres in [
    (Matrix([[2, 1, 0], [0, 2, 0], [1, 0, 3]]), 2, [[1, -1, 0], [1, 0, -1]], [[[1, 0, -1]], [[1, -1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0]]]),
    (Matrix([[3, 0, 1], [0, 3, 0], [0, 1, 1]]), 3, [[1, 0, 0], [0, 2, 1]], [[[1, 0, 0]], [[1, 0, 0], [1, 0, -2]], [[0, 1, 0], [0, 0, 1]]]),
    (Matrix([[2, 0, 1], [1, 1, 1], [0, 0, 1]]), 1, [[0, 1, 0], [1, 0, -1]], [[[0, 1, 0]], [[1, 1, 0], [0, 1, 0]], [[1, 0, 0], [0, 0, 1]]]),
]:
    Nl = (A - lam * I3) ** 2
    assert Nl.rank() == 1
    B_ = Matrix.hstack(*[Matrix(v) for v in base])
    assert B_.rank() == 2 and Nl * B_ == zeros(3, 2)
    for L in leurres:
        BL = Matrix.hstack(*[Matrix(v) for v in L])
        assert BL.rank() < 2 or Nl * BL != zeros(3, BL.shape[1])
assert (Matrix([[2, 1, 0], [0, 2, 0], [1, 0, 3]]) - 2 * I3) ** 2 == Matrix([[0, 0, 0], [0, 0, 0], [1, 1, 1]])
assert (Matrix([[3, 0, 1], [0, 3, 0], [0, 1, 1]]) - 3 * I3) ** 2 == Matrix([[0, 1, -2], [0, 0, 0], [0, -2, 4]])
assert (Matrix([[2, 0, 1], [1, 1, 1], [0, 0, 1]]) - I3) ** 2 == Matrix([[1, 0, 1], [1, 0, 1], [0, 0, 0]])
assert Matrix([[3, 0, 1], [0, 3, 0], [0, 1, 1]]) * Matrix([1, 0, -2]) == Matrix([1, 0, -2])

# sous-espaces-caracteristiques-ex2 (+ v1, v2) : Dunford.
for A, a1, b1, pbm, D, Nm, leurres in [
    (Matrix([[1, -1, 1], [0, 2, 0], [0, 1, 1]]), 1, 2, Matrix([[0, 0, 0], [0, 1, 0], [0, 1, 0]]),
     Matrix([[1, 0, 0], [0, 2, 0], [0, 1, 1]]), Matrix([[0, -1, 1], [0, 0, 0], [0, 0, 0]]),
     [diag(1, 2, 1), Matrix([[1, -1, 0], [0, 2, 0], [0, 1, 1]]), Matrix([[1, 0, 1], [0, 2, 0], [0, 0, 1]])]),
    (Matrix([[2, 1, 1], [0, 1, 0], [0, -1, 2]]), 2, 1, Matrix([[0, -2, 0], [0, 1, 0], [0, 1, 0]]),
     Matrix([[2, 2, 0], [0, 1, 0], [0, -1, 2]]), Matrix([[0, -1, 1], [0, 0, 0], [0, 0, 0]]),
     [diag(2, 1, 2), Matrix([[2, 1, 0], [0, 1, 0], [0, -1, 2]]), Matrix([[2, -2, 0], [0, 1, 0], [0, 1, 2]])]),
    (Matrix([[3, 1, 0], [0, 3, 0], [1, 0, 2]]), 3, 2, Matrix([[0, 0, 0], [0, 0, 0], [-1, 1, 1]]),
     Matrix([[3, 0, 0], [0, 3, 0], [1, -1, 2]]), Matrix([[0, 1, 0], [0, 0, 0], [0, 1, 0]]),
     [diag(3, 3, 2), Matrix([[3, 0, 0], [0, 3, 0], [1, 0, 2]]), Matrix([[3, 0, 0], [0, 3, 0], [-1, 1, 2]])]),
]:
    assert expand(chi(A) - (X - a1) ** 2 * (X - b1)) == 0
    assert expand((X - a1) ** 2 - (X - b1) * (X - 2 * a1 + b1) - 1) == 0          # relation de Bézout
    pb = (A - a1 * I3) ** 2
    assert pb == pbm and pb * pb == pb and (A - b1 * I3) * (A - (2 * a1 - b1) * I3) == pb - I3
    assert a1 * (I3 - pb) + b1 * pb == D and A - D == Nm
    assert Nm**2 == zeros(3) and D * Nm == Nm * D and D.is_diagonalizable()
    assert all(L != D for L in leurres)
    # « partie diagonale » : ne commute pas avec le reste.
    Dg = diag(*[A[j, j] for j in range(3)])
    assert Dg * (A - Dg) != (A - Dg) * Dg

# sous-espaces-caracteristiques-ex4 (+ v1, v2) : systèmes différentiels.
for A, lam, Y0, sol, leurres in [
    (Matrix([[2, 1], [-1, 4]]), 3, [1, 0], [(1 - t) * exp(3 * t), -t * exp(3 * t)],
     [[(1 + t) * exp(3 * t), t * exp(3 * t)], [exp(3 * t), 0], [(1 - t) * exp(2 * t), -t * exp(4 * t)]]),
    (Matrix([[4, -1], [1, 2]]), 3, [1, 0], [(1 + t) * exp(3 * t), t * exp(3 * t)],
     [[(1 - t) * exp(3 * t), -t * exp(3 * t)], [exp(3 * t), t * exp(3 * t)], [(1 + t) * exp(4 * t), t * exp(2 * t)]]),
    (Matrix([[0, 1], [-1, 2]]), 1, [0, 1], [t * exp(t), (1 + t) * exp(t)],
     [[-t * exp(t), (1 - t) * exp(t)], [t * exp(2 * t), (1 + t) * exp(2 * t)], [0, exp(t)]]),
]:
    assert expand(chi(A) - (X - lam) ** 2) == 0 and (A - lam * I2) ** 2 == zeros(2)
    Y = Matrix(sol)
    assert simplify(Y.diff(t) - A * Y) == zeros(2, 1) and Y.subs(t, 0) == Matrix(Y0)
    assert simplify(Y - exp(lam * t) * (I2 + t * (A - lam * I2)) * Matrix(Y0)) == zeros(2, 1)
    for L in leurres:
        YL = Matrix(L)
        assert simplify(YL.diff(t) - A * YL) != zeros(2, 1) or YL.subs(t, 0) != Matrix(Y0)

# ---------------------------------------------------------------------------
# al-l2-reduction-puissances (carte)
# ---------------------------------------------------------------------------
A = Matrix([[1, 2], [2, 1]])
for kk in range(6):
    assert A**kk == Matrix([[3**kk + (-1) ** kk, 3**kk - (-1) ** kk], [3**kk - (-1) ** kk, 3**kk + (-1) ** kk]]) / 2
u_ = [kk**2 + 1 for kk in range(10)]                                                         # exemple (X-1)^3
assert u_[:3] == [1, 2, 5] and all(u_[j + 3] == 3 * u_[j + 2] - 3 * u_[j + 1] + u_[j] for j in range(7))
Aa, Bb = Matrix([[0, 1], [0, 0]]), Matrix([[0, 0], [1, 0]])
assert (Aa + Bb) ** 2 == I2 and Aa**2 + 2 * Aa * Bb + Bb**2 == Matrix([[2, 0], [0, 0]])      # contre-exemple
assert all((kk + 1) * 2 ** (kk + 1) != 2 * kk * 2**kk for kk in range(10))                     # k 2^k non solution
# Trou 3 : leurres faux pour P = [[1,1],[0,1]], D = diag(1,2), k = 2.
P, D = Matrix([[1, 1], [0, 1]]), diag(1, 2)
A = P * D * P.inv()
assert A**2 == P * D**2 * P.inv() and A**2 not in (P**2 * D * P.inv() ** 2, P.inv() * D**2 * P, D**2)
# Trou 2, leurres : u_k = k n'est pas solution de u_{k+2} = 3u_{k+1} - 2u_k (« degré au plus p-1 ») ;
# la suite constante 1 est solution de u_{k+2} = 2u_{k+1} - u_k (« degré exactement m_j - 1 » l'exclurait).
assert expand((k + 2) - 3 * (k + 1) + 2 * k) != 0 and 1 == 2 * 1 - 1
# Affirmations : (alpha + beta k) 2^k solution de u_{k+2} = 4u_{k+1} - 4u_k ; (-1)^k solution de u_{k+2} = u_k.
al, be = symbols("al be")
uu = lambda m: (al + be * m) * 2**m
assert expand(uu(k + 2) - 4 * uu(k + 1) + 4 * uu(k)) == 0
lam_, N2 = symbols("lam"), Matrix([[0, 1], [0, 0]])
for kk in range(1, 6):
    assert ((lam_ * I2 + N2) ** kk - (lam_**kk * I2 + kk * lam_ ** (kk - 1) * N2)).applyfunc(expand) == zeros(2)

# puissances-ex1 (+ v1, v2) : coefficient (1,2) de A^k, leurres faux.
for A, coef, leurres in [
    (Matrix([[2, 1], [1, 2]]), lambda m: F(3**m - 1, 2), [lambda m: F(3**m + 1, 2), lambda m: 3**m - 1, lambda m: F(3**m - 1, 3)]),
    (Matrix([[3, 1], [1, 3]]), lambda m: F(4**m - 2**m, 2), [lambda m: F(4**m + 2**m, 2), lambda m: 4**m - 2**m, lambda m: F(4**m - 2**m, 4)]),
    (Matrix([[1, 3], [3, 1]]), lambda m: F(4**m - (-2) ** m, 2), [lambda m: F(4**m + (-2) ** m, 2), lambda m: F(4**m - 2**m, 2), lambda m: 4**m - (-2) ** m]),
]:
    for m in range(8):
        assert (A**m)[0, 1] == coef(m)
    assert all(any((A**m)[0, 1] != L(m) for m in range(8)) for L in leurres)


def suite(init, rel, n=30):
    u = list(init)
    while len(u) < n:
        u.append(rel(u))
    return u


# puissances-ex2, puissances-ex4 (+ v1, v2) : formules explicites et leurres.
for init, rel, ch, fact, formule, leurres in [
    ((3, 2, 6), lambda u: 2 * u[-1] + u[-2] - 2 * u[-3], X**3 - 2 * X**2 - X + 2, (X - 2) * (X - 1) * (X + 1),
     lambda m: 1 + (-1) ** m + 2**m, [lambda m: 1 - (-1) ** m + 2 ** (m + 1), lambda m: 2 + 2**m, lambda m: 3 * 2**m - m]),
    ((1, 6, 12), lambda u: 4 * u[-1] - u[-2] - 6 * u[-3], X**3 - 4 * X**2 + X + 6, (X + 1) * (X - 2) * (X - 3),
     lambda m: 2**m + 3**m - (-1) ** m, [lambda m: 2**m + 3**m + (-1) ** m, lambda m: 3**m - 2**m + (-1) ** m, lambda m: 2 * 3**m - 1]),
    ((1, 2, 6), lambda u: 6 * u[-1] - 11 * u[-2] + 6 * u[-3], X**3 - 6 * X**2 + 11 * X - 6, (X - 1) * (X - 2) * (X - 3),
     lambda m: 1 - 2**m + 3**m, [lambda m: 1 + 2**m - 3**m, lambda m: 2**m + 3**m - 1, lambda m: 3**m - m]),
    ((1, 1, 6), lambda u: 3 * u[-2] + 2 * u[-3], X**3 - 3 * X - 2, (X + 1) ** 2 * (X - 2),
     lambda m: m * (-1) ** m + 2**m, [lambda m: (-1) ** m + 2**m, lambda m: (1 + m) * (-1) ** m + 2**m, lambda m: m * (-1) ** (m + 1) + 2**m]),
    ((1, 3, 9), lambda u: 5 * u[-1] - 8 * u[-2] + 4 * u[-3], X**3 - 5 * X**2 + 8 * X - 4, (X - 1) * (X - 2) ** 2,
     lambda m: 1 + m * 2**m, [lambda m: 1 + 2**m, lambda m: 1 + m * 2 ** (m + 1), lambda m: 2 ** (m + 1) - 1 + m]),
    ((1, 0, 3), lambda u: -u[-1] + u[-2] + u[-3], X**3 + X**2 - X - 1, (X - 1) * (X + 1) ** 2,
     lambda m: 1 + m * (-1) ** m, [lambda m: 1 + (-1) ** m, lambda m: 1 - m * (-1) ** m, lambda m: (1 + m) * (-1) ** m]),
]:
    assert expand(ch - fact) == 0 and Poly(ch, X).all_coeffs()[-1] != 0
    u = suite(init, rel)
    assert all(u[m] == formule(m) for m in range(30))
    assert all(any(u[m] != L(m) for m in range(3)) for L in leurres)  # chaque leurre rate une condition initiale
# puissances-ex4 : sans le terme k(-1)^k, système incompatible (gamma = 2/3, alpha = 1/3, puis 3 ≠ 6).
assert F(1, 3) + 4 * F(2, 3) == 3 != 6

# puissances-ex3 (+ v1, v2) : coefficient (1,3) de A^k.
for A, lam, c13, leurres in [
    (Matrix([[2, 1, 2], [0, 2, 3], [0, 0, 2]]), 2, lambda m: m * 2**m + F(3 * m * (m - 1) * 2**m, 8),
     [lambda m: m * 2**m, lambda m: F(3 * m * (m - 1) * 2**m, 8), lambda m: F(m * 2**m, 2) + F(3 * m * (m - 1) * 2**m, 8)]),
    (Matrix([[3, 1, 0], [0, 3, 2], [0, 0, 3]]), 3, lambda m: F(m * (m - 1) * 3**m, 9),
     [lambda m: F(m * (m - 1) * 3**m, 18), lambda m: F(m * 3**m, 3), lambda m: F(2 * m * 3**m, 3)]),
    (Matrix([[1, 2, 1], [0, 1, 1], [0, 0, 1]]), 1, lambda m: m * m, [lambda m: m, lambda m: m * (m - 1), lambda m: F(m * (m + 1), 2)]),
]:
    N_ = A - lam * I3
    assert N_**3 == zeros(3)
    for m in range(2, 10):
        assert (A**m)[0, 2] == c13(m)
        assert A**m == lam**m * I3 + m * lam ** (m - 1) * N_ + binomial(m, 2) * lam ** (m - 2) * N_**2
    assert all(any((A**m)[0, 2] != L(m) for m in range(2, 6)) for L in leurres)
# ex3 : coefficients (1,2) et (2,3) annoncés pour lambda = 2 ; (2,3) pour lambda = 3 ; (1,2), (2,3) pour lambda = 1.
for m in range(2, 8):
    A2_ = Matrix([[2, 1, 2], [0, 2, 3], [0, 0, 2]]) ** m
    assert A2_[0, 1] == m * 2 ** (m - 1) and A2_[1, 2] == 3 * m * 2 ** (m - 1)
    assert (Matrix([[3, 1, 0], [0, 3, 2], [0, 0, 3]]) ** m)[1, 2] == 2 * m * 3 ** (m - 1)
    B2_ = Matrix([[1, 2, 1], [0, 1, 1], [0, 0, 1]]) ** m
    assert B2_[0, 1] == 2 * m and B2_[1, 2] == m

# puissances-ex5 (+ v1, v2) : moyennes itérées.
A = Matrix([[0, 1, 1], [1, 0, 1], [1, 1, 0]]) / 2
assert A.is_diagonalizable() and set(A.eigenvals()) == {1, Rational(-1, 2)} and A.eigenvals()[Rational(-1, 2)] == 2
for x0, m, Y0, leurres in [((6, 0, 3), 3, (3, -3, 0), (F(9, 2), F(3, 2), 2)),
                           ((1, 2, 9), 4, (-3, -2, 5), (5, F(11, 2), F(3, 2))),
                           ((2, 7, 9), 6, (-4, 1, 3), (F(11, 2), 8, F(9, 2)))]:
    v = [F(c_) for c_ in x0]
    assert sum(v) == 3 * m and tuple(c_ - m for c_ in x0) == Y0
    for kk in range(40):
        assert all(v[j] == m + Y0[j] * F(-1, 2) ** kk for j in range(3))
        v = [(v[1] + v[2]) / 2, (v[0] + v[2]) / 2, (v[0] + v[1]) / 2]
    assert m not in leurres

# ---------------------------------------------------------------------------
# QCM (content/qcm/al-l2-reduction.json)
# ---------------------------------------------------------------------------
assert expand(X**2 - 5 * X + 6 - (X - 2) * (X - 3)) == 0                                   # q1
assert not J.is_diagonalizable()                                                           # q2
Pq = diag(1, 1, 1, 0, 0)
assert Pq**2 == Pq and Pq.trace() == 3 and Pq.rank() == 3                                  # q3
assert chi(I2) != X - 1                                                                    # q4
assert diag(1, 2, -1) ** 2 + I3 == diag(2, 5, 2) and (diag(1, 2, -1) ** 2 + I3).det() == 20  # q5
assert diag(1, 2, -1).det() ** 2 + 1 == 5 and (diag(1, 2, -1) ** 2).trace() + 3 == 9       # q5 : leurres
for M in [diag(1, 1, 2), Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 2]])]:                        # q6
    assert len(((M - I3) ** 2).nullspace()) == 2
uq = suite((1, 4), lambda u: 4 * u[-1] - 4 * u[-2])                                         # q7
assert all(uq[m] == (1 + m) * 2**m for m in range(30))
assert uq[2] == 12 and 4**2 != 12 and (1 + 2) * 2 != 4 and 2**2 - 1 != 4

print("al_l2_reduction : toutes les vérifications passent.")
