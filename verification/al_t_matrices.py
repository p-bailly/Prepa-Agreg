# Vérification calculatoire du chapitre « al-t-matrices » (thème algebre-lineaire, niveau terminale) :
# cartes content/cartes/algebre-lineaire/terminale-al-t-matrices.json et exercices
# content/exercices/algebre-lineaire/terminale-al-t-matrices.json. Calcul exact (SymPy, rationnels), dénombrements
# de chemins par énumération complète. Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise :
# « v1 », « v2 »). Voir CLAUDE.md, « Contrôle de l'exactitude ».
from itertools import combinations, product
import sympy as sp
from sympy import Matrix as M, Rational as R, eye, zeros, ones, symbols

n_ = symbols("n", integer=True, nonnegative=True)


def puiss(A, k):
    """A^k, avec A^0 = I."""
    return A ** k if k > 0 else eye(A.shape[0])


def adj(p, aretes, oriente=False):
    A = zeros(p, p)
    for i, j in aretes:
        A[i - 1, j - 1] = 1
        if not oriente:
            A[j - 1, i - 1] = 1
    return A


def chemins(A, n, i, j, distincts=False):
    """Nombre de suites de sommets (s0..sn), s0 = i, sn = j, avec une arête de s_{k-1} vers s_k (énumération)."""
    p = A.shape[0]
    c = 0
    for milieu in product(range(p), repeat=n - 1):
        s = (i - 1,) + milieu + (j - 1,)
        if distincts and len(set(s)) < len(s):
            continue
        if all(A[s[k - 1], s[k]] for k in range(1, n + 1)):
            c += 1
    return c


def taille_produit(*tailles):
    """Taille du produit de matrices de tailles données (lignes, colonnes), ou None s'il n'est pas défini."""
    l, c = tailles[0]
    for (l2, c2) in tailles[1:]:
        if c != l2:
            return None
        c = c2
    return (l, c)


E11, E12, E21, E22 = (M([[1, 0], [0, 0]]), M([[0, 1], [0, 0]]), M([[0, 0], [1, 0]]), M([[0, 0], [0, 1]]))
a, b, c, d = symbols("a b c d")
A2 = M([[a, b], [c, d]])

# ==================================================================================================
# al-t-matrices-produit (carte)
# ==================================================================================================
# Trou 1 : (AB)_ij = sum_k a_ik b_kj ; les leurres (Hadamard, A B^T, A^T B) donnent d'autres matrices.
A, B = M([[1, 2], [3, 4]]), M([[0, 1], [5, 7]])
AB = A * B
assert AB == M([[sum(A[i, k] * B[k, j] for k in range(2)) for j in range(2)] for i in range(2)]), "al-t-matrices-produit : formule du produit"
assert A.multiply_elementwise(B) != AB, "al-t-matrices-produit : leurre a_ij b_ij"
assert A * B.T != AB, "al-t-matrices-produit : leurre sum a_ik b_jk"
assert A.T * B != AB, "al-t-matrices-produit : leurre sum a_ki b_kj"
# Trou 2 : (A+B)^2 = A^2 + AB + BA + B^2 (symbolique), égal à A^2 + 2AB + B^2 ssi AB = BA.
Xs = M(2, 2, symbols("x0:4"))
Ys = M(2, 2, symbols("y0:4"))
assert sp.expand((Xs + Ys) ** 2 - (Xs ** 2 + Xs * Ys + Ys * Xs + Ys ** 2)) == zeros(2, 2), "al-t-matrices-produit : carré d'une somme"
assert sp.expand((Xs + Ys) ** 2 - (Xs ** 2 + 2 * Xs * Ys + Ys ** 2) - (Ys * Xs - Xs * Ys)) == zeros(2, 2)
# Leurre « AB = 0 » : A = E11, B = E21 : AB = 0 mais l'identité échoue (BA = B).
assert E11 * E21 == zeros(2, 2) and E21 * E11 == E21, "al-t-matrices-produit : AB = 0, BA = B"
assert (E11 + E21) ** 2 - (E11 ** 2 + 2 * E11 * E21 + E21 ** 2) == E21 != zeros(2, 2), "al-t-matrices-produit : contre-exemple carré"
# Leurre « A ou B inversible » : A = (1 1; 0 1) inversible, B = E11 : l'identité échoue.
Ai = M([[1, 1], [0, 1]])
assert Ai.det() != 0 and (Ai + E11) ** 2 != Ai ** 2 + 2 * Ai * E11 + E11 ** 2, "al-t-matrices-produit : leurre inversible"
# Leurre « A = B » : A = I, B = 2I commutent, l'identité vaut, et A != B.
assert (eye(2) + 2 * eye(2)) ** 2 == eye(2) + 4 * eye(2) + 4 * eye(2), "al-t-matrices-produit : leurre A = B"
# Démonstration : E (1,1) et F (1,2) en ordre n = 2..5 : EF = F, FE = 0.
for n in range(2, 6):
    E = zeros(n, n); E[0, 0] = 1
    F = zeros(n, n); F[0, 1] = 1
    assert E * F == F and F * E == zeros(n, n), "al-t-matrices-produit : EF = F, FE = 0"
# Associativité et distributivité sur des matrices rectangulaires symboliques.
Ar = M(2, 3, symbols("p0:6")); Br = M(3, 4, symbols("q0:12")); Cr = M(4, 2, symbols("r0:8")); Br2 = M(3, 4, symbols("s0:12"))
assert sp.expand((Ar * Br) * Cr - Ar * (Br * Cr)) == zeros(2, 2), "al-t-matrices-produit : associativité"
assert sp.expand(Ar * (Br + Br2) - Ar * Br - Ar * Br2) == zeros(2, 4), "al-t-matrices-produit : distributivité"
assert eye(2) * Ar == Ar and Ar * eye(3) == Ar, "al-t-matrices-produit : neutre"
# Exemple.
A, B = M([[1, 2], [0, 1]]), M([[0, 1], [1, 0]])
assert A * B == M([[2, 1], [1, 0]]) and B * A == M([[0, 1], [1, 2]]), "al-t-matrices-produit : exemple"
assert M([[1, 2, 3]]) * M([4, 5, 6]) == M([[32]]) and (M([4, 5, 6]) * M([[1, 2, 3]])).shape == (3, 3), "al-t-matrices-produit : ligne x colonne"
# Contre-exemples.
assert E12 ** 2 == zeros(2, 2), "al-t-matrices-produit : N^2 = 0"
assert taille_produit((2, 3), (3, 2)) == (2, 2) and taille_produit((3, 2), (2, 3)) == (3, 3) and taille_produit((2, 3), (2, 3)) is None
x1, y1 = symbols("x1 y1")
assert M([[x1]]) * M([[y1]]) == M([[y1]]) * M([[x1]]), "al-t-matrices-produit : ordre 1"
# Affirmations : (AB)^2 != A^2B^2 ; AB = 0 avec A = B = N ; AB et BA définis => carrées.
assert (E12 * E21) == E11 and (E12 * E21) ** 2 == E11 and E12 ** 2 * E21 ** 2 == zeros(2, 2), "al-t-matrices-produit : affirmation 1"
assert E12 * E12 == zeros(2, 2), "al-t-matrices-produit : affirmation 2"
for (nn, pp) in [(2, 3), (3, 2), (1, 4)]:
    assert taille_produit((nn, pp), (pp, nn)) == (nn, nn) and taille_produit((pp, nn), (nn, pp)) == (pp, pp)

# ---------------- al-t-matrices-produit-ex1 (+ v1, v2) : taille de BCA
assert taille_produit((3, 4), (4, 2), (2, 3)) == (3, 3), "al-t-matrices-produit-ex1 : 3 x 3"
assert taille_produit((3, 4), (4, 2)) == (3, 2)
assert taille_produit((2, 5), (5, 3), (3, 2)) == (2, 2), "al-t-matrices-produit-ex1 v1 : 2 x 2"
assert taille_produit((3, 2), (3, 4)) is None and taille_produit((3, 4), (2, 3)) is None, "al-t-matrices-produit-ex1 v2 : non défini"
# Leurres : les autres tailles proposées ne sont pas la bonne (vérifié par différence avec la réponse calculée).
for bonne, leurres in [((3, 3), [(2, 2), (4, 4), None]), ((2, 2), [(3, 3), (5, 5), None]), (None, [(3, 3), (2, 2), (3, 4)])]:
    assert all(l != bonne for l in leurres), "al-t-matrices-produit-ex1 : leurres"

# ---------------- al-t-matrices-produit-ex2 (+ v1, v2) : AB, BA ; leurres BA, Hadamard, A B^T
for (A, B, ABv, BAv, had, ABt, nom) in [
    (M([[1, 2], [0, -1]]), M([[3, 0], [1, 2]]), M([[5, 4], [-1, -2]]), M([[3, 6], [1, 0]]), M([[3, 0], [0, -2]]), M([[3, 5], [0, -2]]), ""),
    (M([[2, -1], [1, 0]]), M([[1, 1], [0, 3]]), M([[2, -1], [1, 1]]), M([[3, -1], [3, 0]]), M([[2, -1], [0, 0]]), M([[1, -3], [1, 0]]), " v1"),
    (M([[0, 1], [2, 3]]), M([[1, -1], [1, 0]]), M([[1, 0], [5, -2]]), M([[-2, -2], [0, 1]]), M([[0, -1], [2, 0]]), M([[-1, 0], [-1, 2]]), " v2"),
]:
    assert A * B == ABv and B * A == BAv, "al-t-matrices-produit-ex2" + nom + " : AB et BA"
    assert A.multiply_elementwise(B) == had and A * B.T == ABt, "al-t-matrices-produit-ex2" + nom + " : leurres"
    assert len({tuple(X) for X in [ABv, BAv, had, ABt]}) == 4, "al-t-matrices-produit-ex2" + nom + " : propositions distinctes, une seule juste"

# ---------------- al-t-matrices-produit-ex3 (+ v1, v2) : commutant
x, y, z, t = symbols("x y z t")
Mg = M([[x, y], [z, t]])
FAM = {
    "haut": M([[x, y], [0, x]]), "bas": M([[x, 0], [z, x]]), "diag": M([[x, 0], [0, t]]),
    "tri": M([[x, y], [0, t]]), "trib": M([[x, 0], [z, t]]), "scal": M([[x, 0], [0, x]]),
}
def dim_commutant(A):
    """Dimension de l'espace des M qui commutent avec A (4 - rang du système linéaire AM - MA = 0)."""
    sysm = sp.linear_eq_to_matrix(list(A * Mg - Mg * A), [x, y, z, t])[0]
    return 4 - sysm.rank()
def commute_toujours(A, fam):
    return sp.expand(A * fam - fam * A) == zeros(2, 2)
def egale_commutant(A, fam):
    # Deux sous-espaces : la famille est contenue dans le commutant et a la même dimension (nombre de paramètres).
    return commute_toujours(A, fam) and len(fam.free_symbols) == dim_commutant(A)
for A, bonne, leurres, nom in [
    (M([[1, 1], [0, 1]]), "haut", ["diag", "tri", "scal"], ""),
    (M([[2, 0], [0, 3]]), "diag", ["haut", "tri", "scal"], " v1"),
    (M([[1, 0], [1, 1]]), "bas", ["haut", "trib", "diag"], " v2"),
]:
    assert egale_commutant(A, FAM[bonne]), "al-t-matrices-produit-ex3" + nom + " : bonne réponse"
    for l in leurres:
        assert not egale_commutant(A, FAM[l]), "al-t-matrices-produit-ex3" + nom + " : leurre " + l

# ---------------- al-t-matrices-produit-ex4 (+ v1, v2) : A^n = lam^n I + n lam^(n-1) N
for A, N, lam, formule, leurres, nom in [
    (M([[2, 1], [0, 2]]), M([[0, 1], [0, 0]]), 2, lambda k: M([[2 ** k, k * 2 ** (k - 1)], [0, 2 ** k]]),
     [lambda k: M([[2 ** k, k * 2 ** k], [0, 2 ** k]]), lambda k: M([[2 ** k, 1], [0, 2 ** k]]), lambda k: M([[2 ** k, R(2) ** (k - 1)], [0, 2 ** k]])], ""),
    (M([[3, 2], [0, 3]]), M([[0, 2], [0, 0]]), 3, lambda k: M([[3 ** k, 2 * k * 3 ** (k - 1)], [0, 3 ** k]]),
     [lambda k: M([[3 ** k, 2 * k * 3 ** k], [0, 3 ** k]]), lambda k: M([[3 ** k, 2 ** k], [0, 3 ** k]]), lambda k: M([[3 ** k, k * 3 ** (k - 1)], [0, 3 ** k]])], " v1"),
    (M([[1, 0], [4, 1]]), M([[0, 0], [4, 0]]), 1, lambda k: M([[1, 0], [4 * k, 1]]),
     [lambda k: M([[1, 0], [4 ** k, 1]]), lambda k: M([[1, 0], [k * 4 ** (k - 1), 1]]), lambda k: M([[k, 0], [4 * k, k]])], " v2"),
]:
    assert A - lam * eye(2) == N and N ** 2 == zeros(2, 2), "al-t-matrices-produit-ex4" + nom + " : N^2 = 0"
    for k in range(1, 13):
        assert A ** k == formule(k) == lam ** k * eye(2) + k * R(lam) ** (k - 1) * N, "al-t-matrices-produit-ex4" + nom + " : A^n"
    for f in leurres:
        assert any(f(k) != A ** k for k in range(1, 6)), "al-t-matrices-produit-ex4" + nom + " : leurre faux"

# ---------------- al-t-matrices-produit-ex5 : s(AB) = s(BA), pas de AB - BA = I
As = M(2, 2, symbols("a0:4")); Bs = M(2, 2, symbols("b0:4"))
assert sp.expand((As * Bs).trace() - (Bs * As).trace()) == 0, "al-t-matrices-produit-ex5 : s(AB) = s(BA)"
assert eye(2).trace() == 2, "al-t-matrices-produit-ex5 : s(I) = 2"
assert E12 * E21 == E11 and E21 * E12 == E22 and E12 * E21 - E21 * E12 == M([[1, 0], [0, -1]]), "al-t-matrices-produit-ex5 : question 4"

# ==================================================================================================
# al-t-matrices-inverse (carte)
# ==================================================================================================
adjA = M([[d, -b], [-c, a]])
assert sp.expand(A2 * adjA) == (a * d - b * c) * eye(2) and sp.expand(adjA * A2) == (a * d - b * c) * eye(2), "al-t-matrices-inverse : A adj(A) = det I"
# Trou 1 : leurres de critère.
assert M([[1, 1], [1, 1]]).det() == 0 and 1 * 1 + 1 * 1 != 0, "al-t-matrices-inverse : leurre ad + bc"
assert eye(2).det() != 0 and (1 * 0 - 0 * 1) == 0, "al-t-matrices-inverse : leurre ab - cd (I inversible, ab - cd = 0)"
assert M([[1, 1], [1, 1]]).det() == 0, "al-t-matrices-inverse : leurre a, d non nuls"
# Trou 2 : leurres de formule (faux pour A = (2 1; 1 3), det 5).
At = M([[2, 1], [1, 3]]); de = At.det()
assert de == 5 and At * (adjA.subs({a: 2, b: 1, c: 1, d: 3}) / de) == eye(2)
for L in [M([[2, -1], [-1, 3]]) / de, M([[3, 1], [1, 2]]) / de, de * M([[3, -1], [-1, 2]])]:
    assert At * L != eye(2), "al-t-matrices-inverse : leurre de la formule de l'inverse"
# Trou 3 et contre-exemple : (AB)^-1 = B^-1 A^-1 != A^-1 B^-1, A^-1 + B^-1, BA.
A, B = M([[1, 1], [0, 1]]), M([[1, 0], [1, 1]])
assert (A * B).inv() == B.inv() * A.inv() == M([[1, -1], [-1, 2]]), "al-t-matrices-inverse : (AB)^-1"
assert A.inv() * B.inv() == M([[2, -1], [-1, 1]]) != (A * B).inv(), "al-t-matrices-inverse : A^-1 B^-1"
assert A.inv() + B.inv() != (A * B).inv() and B * A != (A * B).inv(), "al-t-matrices-inverse : leurres"
# Exemple.
A = M([[2, 5], [1, 3]])
assert A.det() == 1 and A.inv() == M([[3, -5], [-1, 2]]) and A.inv() * M([1, 2]) == M([-7, 3]), "al-t-matrices-inverse : exemple"
# Contre-exemples.
assert M([[1, 2], [2, 4]]).det() == 0, "al-t-matrices-inverse : coefficients non nuls, non inversible"
L1, C1 = M([[1, 0]]), M([1, 0])
assert L1 * C1 == eye(1) and C1 * L1 == M([[1, 0], [0, 0]]) != eye(2), "al-t-matrices-inverse : non carrées"
# Affirmations.
assert (eye(2) + (-eye(2))).det() == 0, "al-t-matrices-inverse : affirmation 1"
assert E12 ** 2 == zeros(2, 2) and E12.det() == 0, "al-t-matrices-inverse : affirmation 2 (exemple)"
A = M([[2, 1], [1, 3]])
assert (A ** 2).inv() == A.inv() ** 2, "al-t-matrices-inverse : affirmation 4"

# ---------------- al-t-matrices-inverse-ex1 (+ v1, v2)
for A, inv, leurres, nom in [
    (M([[3, 5], [1, 2]]), M([[2, -5], [-1, 3]]), [M([[2, 5], [1, 3]]), M([[R(1, 3), R(1, 5)], [1, R(1, 2)]]), M([[3, -5], [-1, 2]])], ""),
    (M([[4, 7], [1, 2]]), M([[2, -7], [-1, 4]]), [M([[2, 7], [1, 4]]), M([[R(1, 4), R(1, 7)], [1, R(1, 2)]]), M([[4, -7], [-1, 2]])], " v1"),
    (M([[2, 3], [4, 5]]), M([[R(-5, 2), R(3, 2)], [2, -1]]), [M([[R(5, 2), R(-3, 2)], [-2, 1]]), M([[R(1, 2), R(1, 3)], [R(1, 4), R(1, 5)]]), M([[-5, 3], [4, -2]])], " v2"),
]:
    assert A.inv() == inv and A * inv == eye(2), "al-t-matrices-inverse-ex1" + nom
    for L in leurres:
        assert A * L != eye(2), "al-t-matrices-inverse-ex1" + nom + " : leurre"
assert M([[2, 3], [4, 5]]).det() == -2

# ---------------- al-t-matrices-inverse-ex2 (+ v1, v2) : systèmes
for A, Y, inv, sol, leurres, nom in [
    (M([[2, 3], [1, 2]]), M([7, 4]), M([[2, -3], [-1, 2]]), M([2, 1]), [(1, 2), (5, -1), (-1, 3)], ""),
    (M([[3, 1], [5, 2]]), M([5, 9]), M([[2, -1], [-5, 3]]), M([1, 2]), [(2, 1), (0, 5), (2, -1)], " v1"),
    (M([[1, -2], [3, 1]]), M([1, 10]), M([[1, 2], [-3, 1]]) / 7, M([3, 1]), [(1, 3), (5, 2), (1, 0)], " v2"),
]:
    assert A.inv() == inv and inv * Y == sol and A * sol == Y, "al-t-matrices-inverse-ex2" + nom
    for l in leurres:
        assert A * M(l) != Y, "al-t-matrices-inverse-ex2" + nom + " : leurre"

# ---------------- al-t-matrices-inverse-ex3 (+ v1, v2) : relation polynomiale
for A, A2v, rel, inv, leurres, nom in [
    (M([[2, 1, 1], [0, 1, 0], [0, 0, 1]]), M([[4, 3, 3], [0, 1, 0], [0, 0, 1]]), lambda X: X ** 2 - 3 * X + 2 * eye(3),
     M([[R(1, 2), R(-1, 2), R(-1, 2)], [0, 1, 0], [0, 0, 1]]),
     [M([[1, -1, -1], [0, 2, 0], [0, 0, 2]]), M([[R(-1, 2), R(1, 2), R(1, 2)], [0, -1, 0], [0, 0, -1]]), M([[R(1, 2), 1, 1], [0, 1, 0], [0, 0, 1]])], ""),
    (M([[3, 2, 0], [0, 1, 0], [0, 0, 3]]), M([[9, 8, 0], [0, 1, 0], [0, 0, 9]]), lambda X: X ** 2 - 4 * X + 3 * eye(3),
     M([[R(1, 3), R(-2, 3), 0], [0, 1, 0], [0, 0, R(1, 3)]]),
     [M([[1, -2, 0], [0, 3, 0], [0, 0, 1]]), M([[R(-1, 3), R(2, 3), 0], [0, -1, 0], [0, 0, R(-1, 3)]]), M([[R(1, 3), R(1, 2), 0], [0, 1, 0], [0, 0, R(1, 3)]])], " v1"),
    (M([[-1, 3, 0], [0, 2, 0], [0, 0, 2]]), M([[1, 3, 0], [0, 4, 0], [0, 0, 4]]), lambda X: X ** 2 - X - 2 * eye(3),
     M([[-1, R(3, 2), 0], [0, R(1, 2), 0], [0, 0, R(1, 2)]]),
     [M([[-2, 3, 0], [0, 1, 0], [0, 0, 1]]), M([[0, R(3, 2), 0], [0, R(3, 2), 0], [0, 0, R(3, 2)]]), M([[-1, R(1, 3), 0], [0, R(1, 2), 0], [0, 0, R(1, 2)]])], " v2"),
]:
    assert A ** 2 == A2v and rel(A) == zeros(3, 3), "al-t-matrices-inverse-ex3" + nom + " : relation"
    assert A.inv() == inv and A * inv == eye(3) == inv * A, "al-t-matrices-inverse-ex3" + nom + " : inverse"
    for L in leurres:
        assert A * L != eye(3), "al-t-matrices-inverse-ex3" + nom + " : leurre"
assert M([[2, 1, 1], [0, 1, 0], [0, 0, 1]]).inv() == (3 * eye(3) - M([[2, 1, 1], [0, 1, 0], [0, 0, 1]])) / 2
assert M([[3, 2, 0], [0, 1, 0], [0, 0, 3]]).inv() == (4 * eye(3) - M([[3, 2, 0], [0, 1, 0], [0, 0, 3]])) / 3
assert M([[-1, 3, 0], [0, 2, 0], [0, 0, 2]]).inv() == (M([[-1, 3, 0], [0, 2, 0], [0, 0, 2]]) - eye(3)) / 2

# ---------------- al-t-matrices-inverse-ex4 (+ v1, v2) : (I - N)^-1 = I + N + N^2
for N, N2, A, inv, X, leurres, nom in [
    (M([[0, 2, -1], [0, 0, 3], [0, 0, 0]]), M([[0, 0, 6], [0, 0, 0], [0, 0, 0]]), M([[1, -2, 1], [0, 1, -3], [0, 0, 1]]), M([[1, 2, 5], [0, 1, 3], [0, 0, 1]]), M([8, 4, 1]),
     [M([[1, 2, -1], [0, 1, 3], [0, 0, 1]]), M([[1, -2, 1], [0, 1, -3], [0, 0, 1]]), M([[1, -2, -5], [0, 1, -3], [0, 0, 1]])], ""),
    (M([[0, 1, 2], [0, 0, -1], [0, 0, 0]]), M([[0, 0, -1], [0, 0, 0], [0, 0, 0]]), M([[1, -1, -2], [0, 1, 1], [0, 0, 1]]), M([[1, 1, 1], [0, 1, -1], [0, 0, 1]]), M([3, 0, 1]),
     [M([[1, 1, 2], [0, 1, -1], [0, 0, 1]]), M([[1, -1, -1], [0, 1, 1], [0, 0, 1]]), M([[1, 1, 3], [0, 1, -1], [0, 0, 1]])], " v1"),
    (M([[0, -3, 1], [0, 0, 2], [0, 0, 0]]), M([[0, 0, -6], [0, 0, 0], [0, 0, 0]]), M([[1, 3, -1], [0, 1, -2], [0, 0, 1]]), M([[1, -3, -5], [0, 1, 2], [0, 0, 1]]), M([-7, 3, 1]),
     [M([[1, -3, 1], [0, 1, 2], [0, 0, 1]]), M([[1, -3, 7], [0, 1, 2], [0, 0, 1]]), M([[1, 3, 5], [0, 1, -2], [0, 0, 1]])], " v2"),
]:
    assert N ** 2 == N2 and N ** 3 == zeros(3, 3), "al-t-matrices-inverse-ex4" + nom + " : N^2, N^3"
    assert eye(3) - N == A and A.inv() == inv == eye(3) + N + N ** 2, "al-t-matrices-inverse-ex4" + nom + " : inverse"
    assert (eye(3) - N) * (eye(3) + N + N ** 2) == eye(3) - N ** 3 == (eye(3) + N + N ** 2) * (eye(3) - N)
    assert inv * M([1, 1, 1]) == X and A * X == M([1, 1, 1]), "al-t-matrices-inverse-ex4" + nom + " : solution"
    for L in leurres:
        assert A * L != eye(3), "al-t-matrices-inverse-ex4" + nom + " : leurre"

# ==================================================================================================
# al-t-matrices-suites (carte)
# ==================================================================================================
def suite(A, C, U0, k):
    U = U0
    for _ in range(k):
        U = A * U + C
    return U
# Formule générale U_n = A^n (U0 - L) + L sur quelques exemples (dont A non inversible).
for A, C, U0 in [(M([[R(1, 2), 1], [0, R(1, 3)]]), M([1, -2]), M([5, 7])), (M([[2, 1], [1, 1]]), M([0, 1]), M([1, 0])), (M([[0, 0], [0, R(1, 2)]]), M([3, 1]), M([1, 1]))]:
    L = (eye(2) - A).inv() * C
    assert A * L + C == L
    for k in range(8):
        assert suite(A, C, U0, k) == puiss(A, k) * (U0 - L) + L, "al-t-matrices-suites : U_n = A^n(U0 - L) + L"
# Trous : contre-exemples en taille 1 (A = 2, C = 1, U0 = 0, L = -1).
A1, C1, U01 = M([[2]]), M([[1]]), M([[0]])
assert suite(A1, C1, U01, 1) == M([[1]]) and suite(A1, C1, U01, 2) == M([[3]])
assert (eye(1) - A1).inv() * C1 == M([[-1]])
for Lf in [M([[0]]), M([[R(1, 2)]])]:  # L = AL (L = 0) ; AL = C et L = A^-1 C (L = 1/2)
    assert A1 * Lf + C1 != Lf and A1 * (U01 - Lf) + Lf != suite(A1, C1, U01, 1), "al-t-matrices-suites : leurres du trou 1"
assert A1 ** 2 * U01 + 2 * C1 != suite(A1, C1, U01, 2), "al-t-matrices-suites : leurre A^n U0 + nC"
assert A1 * (U01 + M([[-1]])) - M([[-1]]) != suite(A1, C1, U01, 1), "al-t-matrices-suites : leurre A^n(U0 + L) - L"
assert A1 * U01 + M([[-1]]) != suite(A1, C1, U01, 1), "al-t-matrices-suites : leurre A^n U0 + L"
# Trou 3 : A = I2 inversible, C = (1, 0) non nulle, coefficients positifs : pas de L.
xx, yy = symbols("xx yy")
assert sp.solve(list(eye(2) * M([xx, yy]) + M([1, 0]) - M([xx, yy])), [xx, yy]) == [], "al-t-matrices-suites : pas de colonne fixe"
# Exemple.
A, C = M([[R(1, 2), 0], [0, R(1, 3)]]), M([1, 2])
assert (eye(2) - A).inv() * C == M([2, 3]), "al-t-matrices-suites : exemple L"
ua, ub = symbols("ua ub")
for k in range(6):
    assert sp.simplify(suite(A, C, M([ua, ub]), k) - M([2 + (ua - 2) / R(2) ** k, 3 + (ub - 3) / R(3) ** k])) == zeros(2, 1), "al-t-matrices-suites : exemple U_n"
# Contre-exemples.
U0g = M([ua, ub])
assert all(sp.simplify(suite(eye(2), M([1, 0]), U0g, k) - (U0g + k * M([1, 0]))) == zeros(2, 1) for k in range(5))
for u0 in [0, 1, 5]:
    for k in range(6):
        assert suite(A1, C1, M([[u0]]), k) == M([[2 ** k * (u0 + 1) - 1]]), "al-t-matrices-suites : divergence A = 2"
assert suite(eye(2), zeros(2, 1), M([1, 1]), 5) == M([1, 1]), "al-t-matrices-suites : A = I, C = 0"
# Affirmations : U2 = 3 != A^2 U0 + 2C = 2 ; U_n = 2^n - 1 diverge.
assert suite(A1, C1, U01, 2) == M([[3]]) and A1 ** 2 * U01 + 2 * C1 == M([[2]]), "al-t-matrices-suites : affirmation 1"
assert all(suite(A1, C1, U01, k) == M([[2 ** k - 1]]) for k in range(10)), "al-t-matrices-suites : affirmation 3"

# ---------------- al-t-matrices-suites-ex1 (+ v1, v2)
for A, U0, U1, U2, U3, A3, leurres, nom in [
    (M([[1, 1], [0, 1]]), M([1, 2]), M([3, 2]), M([5, 2]), M([7, 2]), M([[1, 3], [0, 1]]), [(5, 2), (9, 2), (7, 8)], ""),
    (M([[2, 0], [1, 1]]), M([1, 0]), M([2, 1]), M([4, 3]), M([8, 7]), M([[8, 0], [7, 1]]), [(4, 3), (8, 3), (6, 7)], " v1"),
    (M([[0, 1], [1, 1]]), M([1, 1]), M([1, 2]), M([2, 3]), M([3, 5]), M([[1, 2], [2, 3]]), [(2, 3), (5, 8), (3, 4)], " v2"),
]:
    assert A * U0 == U1 and A * U1 == U2 and A * U2 == U3 and A ** 3 == A3 and A3 * U0 == U3, "al-t-matrices-suites-ex1" + nom
    assert all(M(l) != U3 for l in leurres), "al-t-matrices-suites-ex1" + nom + " : leurres"

# ---------------- al-t-matrices-suites-ex2 (+ v1, v2)
for A, C, L, IA, leurres, nom in [
    (M([[0, R(1, 2)], [R(1, 2), 0]]), M([1, 2]), M([R(8, 3), R(10, 3)]), M([[1, R(-1, 2)], [R(-1, 2), 1]]), [(1, 2), (2, R(5, 2)), (0, 2)], ""),
    (M([[R(1, 2), 0], [R(1, 4), R(1, 2)]]), M([2, 1]), M([4, 4]), M([[R(1, 2), 0], [R(-1, 4), R(1, 2)]]), [(2, 1), (4, 2), (1, 2)], " v1"),
    (M([[0, R(1, 3)], [R(2, 3), 0]]), M([1, 1]), M([R(12, 7), R(15, 7)]), M([[1, R(-1, 3)], [R(-2, 3), 1]]), [(1, 1), (R(15, 7), R(12, 7)), (R(3, 2), 2)], " v2"),
]:
    assert eye(2) - A == IA and IA.det() == {"": R(3, 4), " v1": R(1, 4), " v2": R(7, 9)}[nom], "al-t-matrices-suites-ex2" + nom + " : det(I - A) = 3/4, 1/4, 7/9"
    assert IA.inv() * C == L and A * L + C == L, "al-t-matrices-suites-ex2" + nom + " : L"
    for l in leurres:
        assert A * M(l) + C != M(l), "al-t-matrices-suites-ex2" + nom + " : leurre"
# Le leurre (0, 2) de l'exercice est (I + A)^-1 C (erreur de signe).
assert (eye(2) + M([[0, R(1, 2)], [R(1, 2), 0]])).inv() * M([1, 2]) == M([0, 2])

# ---------------- al-t-matrices-suites-ex3 (+ v1, v2) : matrice triangulaire
for q, C, L, u, v, leurres, nom in [
    (2, M([1, 1]), M([4, 2]), lambda k: 4 - R(2 * k + 4, 2 ** k), lambda k: 2 - R(2, 2 ** k),
     [lambda k: 4 - R(4, 2 ** k), lambda k: 4 - R(k + 4, 2 ** k), lambda k: 4 + R(2 * k + 4, 2 ** k)], ""),
    (2, M([0, 2]), M([4, 4]), lambda k: 4 - R(4 * k + 4, 2 ** k), lambda k: 4 - R(4, 2 ** k),
     [lambda k: 4 - R(4, 2 ** k), lambda k: 4 - R(2 * k + 4, 2 ** k), lambda k: 4 + R(4 * k + 4, 2 ** k)], " v1"),
    (3, M([2, 2]), M([R(9, 2), 3]), lambda k: R(9, 2) - R(6 * k + 9, 2 * 3 ** k), lambda k: 3 - R(3, 3 ** k),
     [lambda k: R(9, 2) - R(9, 2 * 3 ** k), lambda k: R(9, 2) - R(3 * k + 9, 2 * 3 ** k), lambda k: 3 - R(3, 3 ** k)], " v2"),
]:
    A = M([[R(1, q), R(1, q)], [0, R(1, q)]])
    assert (eye(2) - A).inv() * C == L and A * L + C == L, "al-t-matrices-suites-ex3" + nom + " : L"
    for k in range(12):
        assert puiss(A, k) == M([[R(1, q ** k), R(k, q ** k)], [0, R(1, q ** k)]]), "al-t-matrices-suites-ex3" + nom + " : A^n"
        Uk = suite(A, C, zeros(2, 1), k)
        assert Uk == M([u(k), v(k)]) == puiss(A, k) * (-L) + L, "al-t-matrices-suites-ex3" + nom + " : u_n, v_n"
    for f in leurres:
        assert any(f(k) != suite(A, C, zeros(2, 1), k)[0] for k in range(4)), "al-t-matrices-suites-ex3" + nom + " : leurre"
    assert sp.limit(sp.sympify(k_ := symbols("k")) / q ** k_, k_, sp.oo) == 0

# ---------------- al-t-matrices-suites-ex4 (+ v1, v2) : A = P D P^-1
A = M([[R(1, 2), R(1, 4)], [R(1, 4), R(1, 2)]])
P = M([[1, 1], [1, -1]]); D = M([[R(3, 4), 0], [0, R(1, 4)]])
assert P.det() == -2 and P.inv() == P / 2, "al-t-matrices-suites-ex4 : P^-1"
assert P * D == M([[R(3, 4), R(1, 4)], [R(3, 4), R(-1, 4)]]) and P * D * P.inv() == A, "al-t-matrices-suites-ex4 : A = PDP^-1"
assert (eye(2) - A).det() == R(3, 16), "al-t-matrices-suites-ex4 : det(I - A)"
for k in range(10):
    al, be = R(3, 4) ** k, R(1, 4) ** k
    assert puiss(A, k) == P * puiss(D, k) * P.inv() == M([[al + be, al - be], [al - be, al + be]]) / 2, "al-t-matrices-suites-ex4 : A^n"
for C, U0, L, W, u, vv, leurres, nom in [
    (M([1, 1]), M([0, 2]), M([4, 4]), M([-4, -2]), lambda k: 4 - 3 * R(3, 4) ** k - R(1, 4) ** k, lambda k: 4 - 3 * R(3, 4) ** k + R(1, 4) ** k,
     [lambda k: 4 - 3 * R(3, 4) ** k + R(1, 4) ** k, lambda k: 4 - 4 * R(3, 4) ** k, lambda k: -3 * R(3, 4) ** k - R(1, 4) ** k], ""),
    (M([1, 1]), M([2, 0]), M([4, 4]), M([-2, -4]), lambda k: 4 - 3 * R(3, 4) ** k + R(1, 4) ** k, lambda k: 4 - 3 * R(3, 4) ** k - R(1, 4) ** k,
     [lambda k: 4 - 3 * R(3, 4) ** k - R(1, 4) ** k, lambda k: 4 - 2 * R(3, 4) ** k, lambda k: -3 * R(3, 4) ** k + R(1, 4) ** k], " v1"),
    (M([2, 2]), M([4, 0]), M([8, 8]), M([-4, -8]), lambda k: 8 - 6 * R(3, 4) ** k + 2 * R(1, 4) ** k, lambda k: 8 - 6 * R(3, 4) ** k - 2 * R(1, 4) ** k,
     [lambda k: 8 - 6 * R(3, 4) ** k - 2 * R(1, 4) ** k, lambda k: 8 - 4 * R(3, 4) ** k, lambda k: 4 - 6 * R(3, 4) ** k + 2 * R(1, 4) ** k], " v2"),
]:
    assert (eye(2) - A).inv() * C == L and U0 - L == W, "al-t-matrices-suites-ex4" + nom + " : L"
    for k in range(10):
        assert suite(A, C, U0, k) == M([u(k), vv(k)]), "al-t-matrices-suites-ex4" + nom + " : u_n, v_n"
    for f in leurres:
        assert any(f(k) != u(k) for k in range(4)), "al-t-matrices-suites-ex4" + nom + " : leurre"

# ==================================================================================================
# al-t-matrices-adjacence (carte)
# ==================================================================================================
graphes_test = [adj(4, [(1, 2), (2, 3), (3, 4), (4, 1), (1, 3)]), adj(3, [(1, 2), (2, 3), (3, 1), (1, 3)], True),
                adj(4, [(1, 2), (2, 2), (2, 3), (4, 1)], True), adj(3, [(1, 2), (2, 3)])]
for G in graphes_test:
    p = G.shape[0]
    for n in range(1, 5):
        for i in range(1, p + 1):
            for j in range(1, p + 1):
                assert (G ** n)[i - 1, j - 1] == chemins(G, n, i, j), "al-t-matrices-adjacence : (M^n)_ij = nombre de chemins"
# Trous : leurres faux sur le graphe de l'exemple (chemin 1-2-3).
G = adj(3, [(1, 2), (2, 3)])
assert (2 * G)[1, 1] != chemins(G, 2, 2, 2) and (G ** 3)[1, 1] != chemins(G, 2, 2, 2) and G[1, 1] != chemins(G, 2, 2, 2), "al-t-matrices-adjacence : leurres nM, M^(n+1), M^(n-1)"
assert (G ** 2)[1, 1] == 2 and chemins(G, 2, 2, 2, distincts=True) == 0, "al-t-matrices-adjacence : répétitions comptées"
assert (G ** 2)[0, 1] == 0 and chemins(G, 1, 1, 2) == 1, "al-t-matrices-adjacence : leurre « au plus n » (1 chemin de longueur 1 de 1 à 2)"
assert (G ** 2)[0, 2] == 1 and chemins(G, 3, 1, 3) == 0, "al-t-matrices-adjacence : leurre « longueur n+1 »"
# Exemple.
assert G == M([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) and G ** 2 == M([[1, 0, 1], [0, 2, 0], [1, 0, 1]]), "al-t-matrices-adjacence : exemple"
# Contre-exemples : orientation ; deux sommets reliés.
Go = adj(2, [(1, 2)], True)
assert Go == M([[0, 1], [0, 0]]) and Go != Go.T
K2 = adj(2, [(1, 2)])
assert all(K2 ** k == (K2 if k % 2 else eye(2)) for k in range(1, 9)), "al-t-matrices-adjacence : puissances alternées"
# Affirmations : degré ; somme des puissances.
for G in [adj(4, [(1, 2), (2, 3), (3, 4), (4, 1), (1, 3)]), adj(5, [(1, 2), (1, 3), (1, 4), (4, 5)])]:
    p = G.shape[0]
    for i in range(p):
        assert (G ** 2)[i, i] == sum(G[i, :]), "al-t-matrices-adjacence : (M^2)_ii = degré"
    S = sum((G ** k for k in range(1, 4)), zeros(p, p))
    assert all(S[i, j] == sum(chemins(G, k, i + 1, j + 1) for k in range(1, 4)) for i in range(p) for j in range(p))
assert any(x == 0 for x in K2 ** 1), "al-t-matrices-adjacence : affirmation 2 (p = 2)"
assert all(x > 0 for x in eye(2) + K2)

# ---------------- al-t-matrices-adjacence-ex1 (+ v1, v2)
for p, ar, Mv, M2v, i, j, rep, leurres, nom in [
    (4, [(1, 2), (1, 3), (2, 3), (3, 4)], M([[0, 1, 1, 0], [1, 0, 1, 0], [1, 1, 0, 1], [0, 0, 1, 0]]), M([[2, 1, 1, 1], [1, 2, 1, 1], [1, 1, 3, 0], [1, 1, 0, 1]]), 1, 4, 1, [0, 2, 3], ""),
    (5, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (1, 3)], M([[0, 1, 1, 0, 1], [1, 0, 1, 0, 0], [1, 1, 0, 1, 0], [0, 0, 1, 0, 1], [1, 0, 0, 1, 0]]),
     M([[3, 1, 1, 2, 0], [1, 2, 1, 1, 1], [1, 1, 3, 0, 2], [2, 1, 0, 2, 0], [0, 1, 2, 0, 2]]), 1, 4, 2, [1, 0, 3], " v1"),
    (4, [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)], M([[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 0]]), M([[3, 2, 1, 1], [2, 3, 1, 1], [1, 1, 2, 2], [1, 1, 2, 2]]), 3, 4, 2, [0, 1, 4], " v2"),
]:
    G = adj(p, ar)
    assert G == Mv and G ** 2 == M2v, "al-t-matrices-adjacence-ex1" + nom + " : M et M^2"
    assert M2v[i - 1, j - 1] == rep == chemins(G, 2, i, j), "al-t-matrices-adjacence-ex1" + nom + " : réponse"
    assert rep not in leurres
# Listes de chaînes données dans les corrections.
assert [k for k in range(1, 5) if adj(4, [(1, 2), (1, 3), (2, 3), (3, 4)])[0, k - 1] and adj(4, [(1, 2), (1, 3), (2, 3), (3, 4)])[k - 1, 3]] == [3]
assert [k for k in range(1, 6) if adj(5, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (1, 3)])[0, k - 1] and adj(5, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (1, 3)])[k - 1, 3]] == [3, 5]
assert [k for k in range(1, 5) if adj(4, [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)])[2, k - 1] and adj(4, [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)])[k - 1, 3]] == [1, 2]

# ---------------- al-t-matrices-adjacence-ex2 (+ v1, v2) : graphe orienté
def liste_chemins(G, n, i, j):
    p = G.shape[0]
    return [(i,) + tuple(m + 1 for m in mi) + (j,) for mi in product(range(p), repeat=n - 1)
            if all(G[s - 1, t - 1] for s, t in zip((i,) + tuple(m + 1 for m in mi), tuple(m + 1 for m in mi) + (j,)))]
for p, ar, Mv, M2v, M3v, i, j, rep, liste, leurres, nom in [
    (3, [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1)], M([[0, 1, 1], [1, 0, 1], [1, 0, 0]]), M([[2, 0, 1], [1, 1, 1], [0, 1, 1]]), M([[1, 2, 2], [2, 1, 2], [2, 0, 1]]), 1, 3, 2,
     [(1, 2, 1, 3), (1, 3, 1, 3)], [1, 3, 0], ""),
    (3, [(1, 2), (1, 3), (2, 3), (3, 1)], M([[0, 1, 1], [0, 0, 1], [1, 0, 0]]), M([[1, 0, 1], [1, 0, 0], [0, 1, 1]]), M([[1, 1, 1], [0, 1, 1], [1, 0, 1]]), 1, 3, 1,
     [(1, 3, 1, 3)], [2, 3, 0], " v1"),
    (4, [(1, 2), (1, 4), (2, 3), (2, 4), (3, 4), (4, 1), (4, 3)], M([[0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1], [1, 0, 1, 0]]), M([[1, 0, 2, 1], [1, 0, 1, 1], [1, 0, 1, 0], [0, 1, 0, 2]]),
     M([[1, 1, 1, 3], [1, 1, 1, 2], [0, 1, 0, 2], [2, 0, 3, 1]]), 1, 4, 3, [(1, 2, 3, 4), (1, 4, 1, 4), (1, 4, 3, 4)], [2, 1, 4], " v2"),
]:
    G = adj(p, ar, True)
    assert G == Mv and G ** 2 == M2v and G ** 3 == M3v, "al-t-matrices-adjacence-ex2" + nom + " : M, M^2, M^3"
    assert M3v[i - 1, j - 1] == rep == chemins(G, 3, i, j) and sorted(liste_chemins(G, 3, i, j)) == sorted(liste), "al-t-matrices-adjacence-ex2" + nom + " : réponse"
    assert rep not in leurres

# ---------------- al-t-matrices-adjacence-ex3 (+ v1, v2) : graphe complet
for p, m2, m3, off, diag, leurres, nom in [(4, (2, 1), (7, -1), 7, 6, [6, 9, 27], ""), (5, (3, 1), (13, -1), 13, 12, [12, 16, 64], " v1"), (3, (1, 1), (3, -1), 3, 2, [2, 4, 8], " v2")]:
    J = ones(p, p); Kp = J - eye(p)
    assert Kp == adj(p, list(combinations(range(1, p + 1), 2))) and J ** 2 == p * J, "al-t-matrices-adjacence-ex3" + nom + " : J^2 = pJ"
    assert Kp ** 2 == m2[0] * J + m2[1] * eye(p) and Kp ** 3 == m3[0] * J + m3[1] * eye(p), "al-t-matrices-adjacence-ex3" + nom + " : M^2, M^3"
    assert chemins(Kp, 3, 1, 2) == off and chemins(Kp, 3, 1, 1) == diag and (p - 1) * off + diag == (p - 1) ** 3, "al-t-matrices-adjacence-ex3" + nom + " : dénombrement"
    assert off not in leurres

# ---------------- al-t-matrices-adjacence-ex4 (+ v1, v2) : triangles
def triangles(G):
    p = G.shape[0]
    return [s for s in combinations(range(1, p + 1), 3) if all(G[u - 1, v - 1] for u, v in combinations(s, 2))]
for G in [adj(5, [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5)]), adj(6, [(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 6), (6, 4), (1, 4)]), ones(4, 4) - eye(4)]:
    for i in range(1, G.shape[0] + 1):
        assert (G ** 3)[i - 1, i - 1] == 2 * sum(1 for s in triangles(G) if i in s), "al-t-matrices-adjacence-ex4 : (M^3)_ii = 2 t_i"
    assert (G ** 3).trace() == 6 * len(triangles(G)), "al-t-matrices-adjacence-ex4 : trace = 6 t"
for G, diag, somme, tri, leurres, nom in [
    (adj(5, [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5)]), [2, 4, 6, 4, 2], 18, [(1, 2, 3), (2, 3, 4), (3, 4, 5)], [6, 18, 9], ""),
    (adj(5, [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4), (4, 5)]), [4, 2, 4, 2, 0], 12, [(1, 2, 3), (1, 3, 4)], [4, 12, 6], " v1"),
    (ones(4, 4) - eye(4), [6, 6, 6, 6], 24, list(combinations(range(1, 5), 3)), [6, 24, 8], " v2"),
]:
    assert [(G ** 3)[i, i] for i in range(G.shape[0])] == diag and sum(diag) == somme, "al-t-matrices-adjacence-ex4" + nom + " : diagonale de M^3"
    assert triangles(G) == tri and somme == 6 * len(tri), "al-t-matrices-adjacence-ex4" + nom + " : triangles"
    assert len(tri) not in leurres

# ---------------- al-t-matrices-adjacence-ex5 (+ v1, v2) : carré
C4 = adj(4, [(1, 2), (2, 3), (3, 4), (4, 1)])
assert C4 == M([[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]) and C4 ** 2 == M([[2, 0, 2, 0], [0, 2, 0, 2], [2, 0, 2, 0], [0, 2, 0, 2]])
assert C4 ** 3 == 4 * C4, "al-t-matrices-adjacence-ex5 : M^3 = 4M"
for k in range(1, 7):
    assert C4 ** (2 * k) == 4 ** (k - 1) * C4 ** 2 and C4 ** (2 * k + 1) == 4 ** k * C4, "al-t-matrices-adjacence-ex5 : puissances"
for n, i, j, rep, leurres, nom in [(10, 1, 3, 512, [256, 1024, 0], ""), (11, 1, 4, 1024, [512, 2048, 0], " v1"), (8, 2, 2, 128, [64, 256, 0], " v2")]:
    assert (C4 ** n)[i - 1, j - 1] == rep and rep not in leurres, "al-t-matrices-adjacence-ex5" + nom
assert (C4 ** 10)[0, 2] == 4 ** 4 * 2 and (C4 ** 11)[0, 3] == 4 ** 5 and (C4 ** 8)[1, 1] == 4 ** 3 * 2
assert chemins(C4, 8, 2, 2) == 128, "al-t-matrices-adjacence-ex5 v2 : énumération"

# ==================================================================================================
# al-t-matrices-markov (carte)
# ==================================================================================================
def loi_enumeree(pi0, P, n):
    """Loi de X_n par énumération des trajectoires (probabilités totales, propriété de Markov)."""
    p = P.shape[0]
    loi = [0] * p
    for traj in product(range(p), repeat=n + 1):
        pr = pi0[traj[0]]
        for k in range(n):
            pr *= P[traj[k], traj[k + 1]]
        loi[traj[-1]] += pr
    return M([loi])
P3 = M([[R(1, 2), R(1, 2), 0], [R(1, 4), R(1, 2), R(1, 4)], [R(1, 3), 0, R(2, 3)]])
pi0 = M([[R(1, 5), R(3, 5), R(1, 5)]])
for n in range(5):
    assert loi_enumeree(pi0, P3, n) == pi0 * puiss(P3, n), "al-t-matrices-markov : pi_n = pi_0 P^n"
    for i in range(3):
        e = zeros(1, 3); e[i] = 1
        assert loi_enumeree(e, P3, n) == puiss(P3, n)[i, :], "al-t-matrices-markov : P(X_n = j | X_0 = i) = (P^n)_ij"
# Deux états : invariante et convergence (symbolique).
aa, bb, xx0 = symbols("aa bb xx0", positive=True)
P2s = M([[1 - aa, aa], [bb, 1 - bb]])
pis = M([[bb / (aa + bb), aa / (aa + bb)]])
assert sp.simplify(pis * P2s - pis) == zeros(1, 2), "al-t-matrices-markov : distribution invariante à deux états"
xs = symbols("xs")
assert sp.solve(sp.Eq(xs * (1 - aa) + (1 - xs) * bb, xs), xs) == [bb / (aa + bb)], "al-t-matrices-markov : unicité"
for k in range(6):
    piK = M([[xx0, 1 - xx0]]) * puiss(P2s, k)
    assert sp.simplify(piK[0] - (bb / (aa + bb) + (1 - aa - bb) ** k * (xx0 - bb / (aa + bb)))) == 0, "al-t-matrices-markov : x_n - x = (1-a-b)^n (x_0 - x)"
# Trous : leurres.
Pex = M([[R(9, 10), R(1, 10)], [R(6, 10), R(4, 10)]]); p0 = M([[1, 0]])
try:
    Pex * M([[1, 0], [0, 1]]) ** 0 * p0
    raise AssertionError("al-t-matrices-markov : P^n pi_0 devrait être non défini")
except sp.ShapeError:
    pass
assert p0 * Pex ** 3 != p0 * Pex ** 2 and 2 * p0 * Pex != p0 * Pex ** 2, "al-t-matrices-markov : leurres pi_0 P^(n+1), n pi_0 P"
pi_inv = M([[R(6, 7), R(1, 7)]])
assert pi_inv * Pex == pi_inv, "al-t-matrices-markov : exemple invariante 6/7, 1/7"
assert sum(pi_inv * Pex) == 1, "al-t-matrices-markov : leurre pi P = 0 impossible"
assert p0 * Pex != p0, "al-t-matrices-markov : leurre pi_0 P = pi (pi = pi_1 non invariante)"
assert (p0 * Pex) * Pex != p0 * Pex
a_, b_ = R(1, 10), R(6, 10)
for L in [M([[a_ / (a_ + b_), b_ / (a_ + b_)]]), M([[R(1, 2), R(1, 2)]]), M([[1 - a_, a_]])]:
    assert L * Pex != L, "al-t-matrices-markov : leurre du trou 3"
# Exemple : P(X_n = 1) = 6/7 + 0.3^n (x0 - 6/7).
for k in range(6):
    assert (M([[xx0, 1 - xx0]]) * puiss(Pex, k))[0].expand() == (R(6, 7) + R(3, 10) ** k * (xx0 - R(6, 7))).expand(), "al-t-matrices-markov : exemple"
# Contre-exemples.
assert all(sp.simplify(M([[xs, 1 - xs]]) * eye(2) - M([[xs, 1 - xs]])) == zeros(1, 2) for _ in [0]), "al-t-matrices-markov : P = I"
S = M([[0, 1], [1, 0]])
assert M([[R(1, 2), R(1, 2)]]) * S == M([[R(1, 2), R(1, 2)]])
assert sp.solve(list(M([[xs, 1 - xs]]) * S - M([[xs, 1 - xs]])), xs) == {xs: R(1, 2)} or sp.solve(list(M([[xs, 1 - xs]]) * S - M([[xs, 1 - xs]])), xs) == [R(1, 2)]
assert all(M([[1, 0]]) * S ** k == (M([[0, 1]]) if k % 2 else M([[1, 0]])) for k in range(8)), "al-t-matrices-markov : alternance"
assert Pex * M([1, 0]) == M([R(9, 10), R(6, 10)]) and sum(Pex * M([1, 0])) == R(3, 2), "al-t-matrices-markov : colonne"
assert M([[1, 0]]) * Pex == M([[R(9, 10), R(1, 10)]])
# Affirmations.
assert S[0, 0] == 0 and (S ** 2)[0, 0] == 1, "al-t-matrices-markov : (P^n)_ij != (p_ij)^n"

# ---------------- al-t-matrices-markov-ex1 (+ v1, v2)
for P, p0, p1, p2, leurres, nom in [
    (M([[R(8, 10), R(2, 10)], [R(3, 10), R(7, 10)]]), M([[1, 0]]), M([[R(8, 10), R(2, 10)]]), M([[R(7, 10), R(3, 10)]]),
     [(R(8, 10), R(2, 10)), (R(64, 100), R(4, 100)), (R(6, 10), R(4, 10))], ""),
    (M([[R(9, 10), R(1, 10)], [R(4, 10), R(6, 10)]]), M([[0, 1]]), M([[R(4, 10), R(6, 10)]]), M([[R(6, 10), R(4, 10)]]),
     [(R(4, 10), R(6, 10)), (R(16, 100), R(36, 100)), (R(8, 10), R(2, 10))], " v1"),
    (M([[R(1, 2), R(1, 2)], [R(2, 10), R(8, 10)]]), M([[R(1, 2), R(1, 2)]]), M([[R(35, 100), R(65, 100)]]), M([[R(305, 1000), R(695, 1000)]]),
     [(R(35, 100), R(65, 100)), (R(3, 10), R(7, 10)), (R(695, 1000), R(305, 1000))], " v2"),
]:
    assert p0 * P == p1 and p1 * P == p2 == p0 * P ** 2, "al-t-matrices-markov-ex1" + nom
    assert all(M([l]) != p2 for l in leurres), "al-t-matrices-markov-ex1" + nom + " : leurres"
assert R(8, 10) * R(8, 10) + R(2, 10) * R(3, 10) == R(7, 10)

# ---------------- al-t-matrices-markov-ex2 (+ v1, v2)
for reste, repare, pi, leurres, nom in [
    (R(7, 10), R(6, 10), M([[R(2, 3), R(1, 3)]]), [(R(1, 3), R(2, 3)), (R(1, 2), R(1, 2)), (R(7, 10), R(3, 10))], ""),
    (R(9, 10), R(4, 10), M([[R(4, 5), R(1, 5)]]), [(R(1, 5), R(4, 5)), (R(9, 10), R(1, 10)), (R(1, 2), R(1, 2))], " v1"),
    (R(5, 10), R(2, 10), M([[R(2, 7), R(5, 7)]]), [(R(5, 7), R(2, 7)), (R(1, 2), R(1, 2)), (R(1, 5), R(4, 5))], " v2"),
]:
    P = M([[reste, 1 - reste], [repare, 1 - repare]])
    assert pi * P == pi and sum(pi) == 1, "al-t-matrices-markov-ex2" + nom + " : invariante"
    assert sp.solve(sp.Eq(reste * xs + repare * (1 - xs), xs), xs) == [pi[0]], "al-t-matrices-markov-ex2" + nom + " : unicité"
    for l in leurres:
        assert M([l]) * P != M([l]), "al-t-matrices-markov-ex2" + nom + " : leurre"

# ---------------- al-t-matrices-markov-ex3 (+ v1, v2)
X3, Y3, Z3 = symbols("X3 Y3 Z3")
for P, pi, leurres, nom in [
    (M([[R(1, 2), R(1, 2), 0], [R(1, 4), R(1, 2), R(1, 4)], [0, R(1, 2), R(1, 2)]]), (R(1, 4), R(1, 2), R(1, 4)),
     [(R(1, 3), R(1, 3), R(1, 3)), (R(1, 4), R(1, 4), R(1, 2)), (R(1, 2), R(1, 2), 0)], ""),
    (M([[R(1, 2), R(1, 2), 0], [0, R(1, 2), R(1, 2)], [R(1, 2), 0, R(1, 2)]]), (R(1, 3), R(1, 3), R(1, 3)),
     [(R(1, 4), R(1, 2), R(1, 4)), (R(1, 2), R(1, 4), R(1, 4)), (R(1, 2), 0, R(1, 2))], " v1"),
    (M([[R(2, 10), R(8, 10), 0], [0, R(5, 10), R(5, 10)], [R(5, 10), 0, R(5, 10)]]), (R(5, 21), R(8, 21), R(8, 21)),
     [(R(1, 3), R(1, 3), R(1, 3)), (R(8, 21), R(5, 21), R(8, 21)), (R(1, 5), R(2, 5), R(2, 5))], " v2"),
]:
    assert all(sum(P[i, :]) == 1 for i in range(3))
    sol = sp.solve(list(M([[X3, Y3, Z3]]) * P - M([[X3, Y3, Z3]])) + [X3 + Y3 + Z3 - 1], [X3, Y3, Z3], dict=True)
    assert sol == [{X3: pi[0], Y3: pi[1], Z3: pi[2]}], "al-t-matrices-markov-ex3" + nom + " : unique invariante"
    for l in leurres:
        assert M([l]) * P != M([l]), "al-t-matrices-markov-ex3" + nom + " : leurre"

# ---------------- al-t-matrices-markov-ex4 (+ v1, v2)
for P, p0, q, bconst, ell, xn, leurres, nom in [
    (M([[R(8, 10), R(2, 10)], [R(3, 10), R(7, 10)]]), M([[1, 0]]), R(5, 10), R(3, 10), R(3, 5), lambda k: R(3, 5) + R(2, 5) * R(1, 2) ** k,
     [lambda k: R(6, 10) - R(4, 10) * R(1, 2) ** k, lambda k: R(6, 10) + R(4, 10) * R(8, 10) ** k, lambda k: R(8, 10) ** k], ""),
    (M([[R(6, 10), R(4, 10)], [R(2, 10), R(8, 10)]]), M([[0, 1]]), R(4, 10), R(2, 10), R(1, 3), lambda k: R(1, 3) - R(1, 3) * R(4, 10) ** k,
     [lambda k: R(1, 3) + R(1, 3) * R(4, 10) ** k, lambda k: R(1, 3) - R(1, 3) * R(6, 10) ** k, lambda k: R(2, 3) - R(2, 3) * R(4, 10) ** k], " v1"),
    (M([[R(1, 10), R(9, 10)], [R(3, 10), R(7, 10)]]), M([[1, 0]]), R(-2, 10), R(3, 10), R(1, 4), lambda k: R(1, 4) + R(3, 4) * R(-2, 10) ** k,
     [lambda k: R(1, 4) + R(3, 4) * R(2, 10) ** k, lambda k: R(3, 4) + R(1, 4) * R(-2, 10) ** k, lambda k: R(1, 4) - R(3, 4) * R(-2, 10) ** k], " v2"),
]:
    assert P[0, 0] - P[1, 0] == q and P[1, 0] == bconst, "al-t-matrices-markov-ex4" + nom + " : x_(n+1) = q x_n + b"
    assert ell == q * ell + bconst and M([[ell, 1 - ell]]) * P == M([[ell, 1 - ell]]), "al-t-matrices-markov-ex4" + nom + " : invariante"
    for k in range(12):
        assert (p0 * puiss(P, k))[0] == xn(k), "al-t-matrices-markov-ex4" + nom + " : x_n"
    for f in leurres:
        assert any(f(k) != xn(k) for k in range(4)), "al-t-matrices-markov-ex4" + nom + " : leurre"
assert R(6, 10) == R(3, 5)

# ---------------- al-t-matrices-markov-ex5 (+ v1, v2) : marche sur un graphe complet
for p, P, cc, diag, leurres, nom in [
    (3, (ones(3, 3) - eye(3)) / 2, R(-1, 2), lambda k: R(1, 3) + R(2, 3) * R(-1, 2) ** k,
     [lambda k: R(1, 3) - R(1, 3) * R(-1, 2) ** k, lambda k: R(1, 2) ** k, lambda k: R(1, 3) + R(2, 3) * R(1, 2) ** k], ""),
    (4, (ones(4, 4) - eye(4)) / 3, R(-1, 3), lambda k: R(1, 4) + R(3, 4) * R(-1, 3) ** k,
     [lambda k: R(1, 4) - R(1, 4) * R(-1, 3) ** k, lambda k: R(1, 3) ** k, lambda k: R(1, 4) + R(3, 4) * R(1, 3) ** k], " v1"),
    (3, (eye(3) + ones(3, 3)) / 4, R(1, 4), lambda k: R(1, 3) + R(2, 3) * R(1, 4) ** k,
     [lambda k: R(1, 3) - R(1, 3) * R(1, 4) ** k, lambda k: R(1, 2) ** k, lambda k: R(1, 3) + R(2, 3) * R(1, 2) ** k], " v2"),
]:
    J = ones(p, p); Q = J / p; I = eye(p)
    assert all(sum(P[i, :]) == 1 for i in range(p)) and P == P.T, "al-t-matrices-markov-ex5" + nom + " : matrice de transition symétrique"
    assert J ** 2 == p * J and Q ** 2 == Q and Q * (I - Q) == zeros(p, p) and (I - Q) ** 2 == I - Q, "al-t-matrices-markov-ex5" + nom + " : Q"
    assert P == Q + cc * (I - Q), "al-t-matrices-markov-ex5" + nom + " : P = Q + c(I - Q)"
    for k in range(10):
        assert puiss(P, k) == Q + cc ** k * (I - Q) and puiss(P, k)[0, 0] == diag(k), "al-t-matrices-markov-ex5" + nom + " : P^n"
    for f in leurres:
        assert any(f(k) != diag(k) for k in range(4)), "al-t-matrices-markov-ex5" + nom + " : leurre"
    pistar = ones(1, p) / p
    assert pistar * P == pistar and (P.T - I).nullspace() and len((P.T - I).nullspace()) == 1, "al-t-matrices-markov-ex5" + nom + " : invariante unique"
# Le leurre 1/3 - 1/3 (-1/2)^n est P(X_n = 2).
assert all(puiss((ones(3, 3) - eye(3)) / 2, k)[0, 1] == R(1, 3) - R(1, 3) * R(-1, 2) ** k for k in range(8))
# La variante v2 : P(X_1 = 1) = 1/2 (on reste sur place).
assert ((eye(3) + ones(3, 3)) / 4)[0, 0] == R(1, 2)

# al-t-matrices-suites (piège « C n'est pas la colonne nulle ») : avec C = 0 la formule vaut (L = 0, U_n = A^n U_0).
Az = M([[R(1, 2), 1], [0, R(1, 3)]])
assert all(suite(Az, zeros(2, 1), M([3, 5]), k) == puiss(Az, k) * M([3, 5]) for k in range(6)), "al-t-matrices-suites : C = 0 permis"
# al-t-matrices-markov (piège « pi_0 invariante ») : pi_0 = (1, 0) non invariante, et pi_n converge quand même vers (6/7, 1/7).
assert M([[1, 0]]) * Pex != M([[1, 0]]) and sp.limit((M([[1, 0]]) * puiss(Pex, 0))[0] * 0 + R(6, 7) + R(3, 10) ** n_ * (1 - R(6, 7)), n_, sp.oo) == R(6, 7)
assert all((M([[1, 0]]) * puiss(Pex, k))[0] == R(6, 7) + R(3, 10) ** k * R(1, 7) for k in range(8)), "al-t-matrices-markov : convergence sans pi_0 invariante"
# al-t-matrices-produit (cours) : pour n = 3 la formule du binôme peut valoir sans commutation (E12, E23 d'ordre 3).
E12_3 = zeros(3, 3); E12_3[0, 1] = 1
E23_3 = zeros(3, 3); E23_3[1, 2] = 1
assert E12_3 * E23_3 != E23_3 * E12_3 and (E12_3 + E23_3) ** 3 == E12_3 ** 3 + 3 * E12_3 ** 2 * E23_3 + 3 * E12_3 * E23_3 ** 2 + E23_3 ** 3, "al-t-matrices-produit : « garanties » seulement"

print("al_t_matrices : toutes les vérifications passent.")
