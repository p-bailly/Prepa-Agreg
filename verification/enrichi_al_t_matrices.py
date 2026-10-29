# Vérification calculatoire de l'enrichissement du chapitre « al-t-matrices » (thème algebre-lineaire, niveau
# terminale) : questions de content/qcm/al-t-matrices.json et exercices ajoutés à
# content/exercices/algebre-lineaire/terminale-al-t-matrices.json (inverse-ex5, adjacence-ex6, markov-ex6, variantes
# « v1 », « v2 » comprises). Calcul exact (SymPy, rationnels), chemins comptés par énumération. Chaque assert cite
# l'identifiant de la question ou de l'exercice. Voir CLAUDE.md, « Contrôle de l'exactitude ».
from itertools import product
import sympy as sp
from sympy import Matrix as M, Rational as R, eye, zeros, ones, symbols

m, n_, a, b, c, d = symbols("m n a b c d")


def adj(p, aretes, oriente=False):
    A = zeros(p, p)
    for i, j in aretes:
        A[i - 1, j - 1] = 1
        if not oriente:
            A[j - 1, i - 1] = 1
    return A


def chaines(A, n, i, j):
    """Liste des suites de sommets (s0..sn), s0 = i, sn = j, reliées par des arêtes (énumération)."""
    p = A.shape[0]
    res = []
    for milieu in product(range(p), repeat=n - 1):
        s = (i - 1,) + milieu + (j - 1,)
        if all(A[s[k - 1], s[k]] for k in range(1, n + 1)):
            res.append(tuple(x + 1 for x in s))
    return res


# al-t-matrices-q1 : AB échange les colonnes de A ; leurres = BA, produit terme à terme, A
A, B = M([[1, 2], [3, 4]]), M([[0, 1], [1, 0]])
assert A * B == M([[2, 1], [4, 3]])
assert B * A == M([[3, 4], [1, 2]]) and A.multiply_elementwise(B) == M([[0, 2], [3, 0]])
for l in (M([[3, 4], [1, 2]]), M([[0, 2], [3, 0]]), A):
    assert A * B != l

# al-t-matrices-q2 (faux) : contre-exemple (A+B)(A-B) ≠ A² - B²
A, B = M([[1, 0], [0, 0]]), M([[0, 1], [0, 0]])
assert A * B == M([[0, 1], [0, 0]]) and B * A == zeros(2, 2)
assert (A + B) * (A - B) - (A ** 2 - B ** 2) == B * A - A * B == M([[0, -1], [0, 0]])

# al-t-matrices-q3 : det = m(m-3) - 4 = (m+1)(m-4)
dt = M([[m, 2], [2, m - 3]]).det()
assert sp.expand(dt - (m + 1) * (m - 4)) == 0 and set(sp.solve(dt, m)) == {-1, 4}
# leurres : pour m = 1, -4, 0, 3 la matrice... m=1 et m=-4, m=0, m=3 sont inversibles (donc les ensembles proposés sont faux)
for v in (1, -4, 0, 3):
    assert dt.subs(m, v) != 0

# al-t-matrices-q4 (vrai) : exemple A² = A, non inversible
A = M([[1, 0], [0, 0]])
assert A ** 2 == A and A.det() == 0

# al-t-matrices-q5 : A^n = [[1, 2n],[0,1]], U_10 = (21, 1)
A, U0 = M([[1, 2], [0, 1]]), M([1, 1])
for k in range(0, 15):
    assert A ** k == M([[1, 2 * k], [0, 1]])
assert A ** 10 * U0 == M([21, 1])
for l in ([11, 1], [20, 1], [1024, 1]):
    assert A ** 10 * U0 != M(l)

# al-t-matrices-q6 (vrai) : M^n symétrique pour un graphe non orienté (contrôle sur des graphes au hasard)
import random
random.seed(3)
for _ in range(30):
    p = random.randint(2, 6)
    aretes = [(i, j) for i in range(1, p + 1) for j in range(i + 1, p + 1) if random.random() < 0.5]
    G = adj(p, aretes)
    for k in range(1, 5):
        assert (G ** k).T == G ** k

# al-t-matrices-q7 : graphe orienté, (M^3)_{11} = 2, chemins (1,2,4,1) et (1,3,4,1)
G = adj(4, [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 1)], oriente=True)
assert G ** 3 == M([[2, 0, 0, 1], [1, 1, 1, 0], [0, 1, 1, 0], [0, 0, 1, 2]])
assert sorted(chaines(G, 3, 1, 1)) == [(1, 2, 4, 1), (1, 3, 4, 1)]
assert len(chaines(G, 4, 1, 1)) >= 1 and (1, 2, 3, 4, 1) in chaines(G, 4, 1, 1)

# al-t-matrices-q8 : (P^2)_{11} = 0,61 ; leurres 0,49 ; 0,7 ; 0,42 faux
P = M([[R(7, 10), R(3, 10)], [R(4, 10), R(6, 10)]])
assert (P ** 2)[0, 0] == R(61, 100)
for l in (R(49, 100), R(7, 10), R(42, 100)):
    assert (P ** 2)[0, 0] != l

# al-t-matrices-inverse-ex5
J = ones(3, 3)
assert J ** 2 == 3 * J                                                           # Q1
assert sp.simplify((a * eye(3) + b * J) * (c * eye(3) + d * J) - (a * c * eye(3) + (a * d + b * c + 3 * b * d) * J)) == zeros(3, 3)  # Q2
assert sp.simplify((a * eye(3) + b * J) * (c * eye(3) + d * J) - (c * eye(3) + d * J) * (a * eye(3) + b * J)) == zeros(3, 3)
Ainv = eye(3) / a - b / (a * (a + 3 * b)) * J                                     # Q3
assert sp.simplify((a * eye(3) + b * J) * Ainv - eye(3)) == zeros(3, 3)
assert sp.simplify(Ainv * (a * eye(3) + b * J) - eye(3)) == zeros(3, 3)
# Q4 : a + 3b = 0 -> (1,1,1) dans le noyau ; a = 0 -> (1,-1,0) dans le noyau
assert ((-3 * b) * eye(3) + b * J) * M([1, 1, 1]) == zeros(3, 1)
assert (b * J) * M([1, -1, 0]) == zeros(3, 1)
# Q5 (exercice, v1, v2) : (a, b, Y, solution, leurres)
cas = [
    (2, 1, [4, 6, 10], [0, 1, 3], [[2, 3, 5], [0, 3, 1], [-2, -1, 1]]),
    (3, 1, [6, 3, 9], [1, 0, 2], [[2, 1, 3], [0, 1, 2], [1, 2, 0]]),
    (-1, 1, [1, 0, 3], [1, 2, -1], [[-1, 0, -3], [2, 1, -1], [1, -1, 2]]),
]
systemes = [M([[3, 1, 1], [1, 3, 1], [1, 1, 3]]), M([[4, 1, 1], [1, 4, 1], [1, 1, 4]]), M([[0, 1, 1], [1, 0, 1], [1, 1, 0]])]
for (av, bv, Y, X, leurres), S in zip(cas, systemes):
    Am = av * eye(3) + bv * J
    assert Am == S                                       # la matrice du système est bien aI + bJ
    assert av != 0 and av + 3 * bv != 0
    inv = Ainv.subs({a: av, b: bv})
    assert inv == Am.inv()
    assert inv * M(Y) == M(X) and Am * M(X) == M(Y)
    assert J * M(Y) == sum(Y) * M([1, 1, 1])
    for l in leurres:
        assert Am * M(l) != M(Y)
# valeurs de A^{-1} citées : 1/2 I - 1/10 J ; 1/3 I - 1/18 J ; -I + 1/2 J
assert Ainv.subs({a: 2, b: 1}) == eye(3) / 2 - J / 10
assert Ainv.subs({a: 3, b: 1}) == eye(3) / 3 - J / 18
assert Ainv.subs({a: -1, b: 1}) == -eye(3) + J / 2

# al-t-matrices-adjacence-ex6 (exercice, v1, v2) : lignes 1 de M, M², M³ ; distance ; chaînes minimales
cas = [
    ([(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)], [[0, 1, 1, 0, 0], [2, 0, 0, 2, 0], [0, 4, 4, 0, 2]], 3,
     [(1, 2, 4, 5), (1, 3, 4, 5)], [(3, 1), (2, 2), (4, 2)]),
    ([(1, 2), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)], [[0, 1, 1, 1, 0], [3, 0, 0, 0, 3], [0, 6, 6, 6, 0]], 2,
     [(1, 2, 5), (1, 3, 5), (1, 4, 5)], [(2, 1), (3, 6), (1, 3)]),
    ([(1, 2), (2, 3), (2, 4), (3, 4), (4, 5)], [[0, 1, 0, 0, 0], [1, 0, 1, 1, 0], [0, 3, 1, 1, 1]], 3,
     [(1, 2, 4, 5)], [(3, 3), (4, 1), (2, 1)]),
]
for aretes, lignes, dist, minimales, leurres in cas:
    G = adj(5, aretes)
    for k in (1, 2, 3):
        assert list((G ** k).row(0)) == lignes[k - 1]
    # distance = plus petit n tel que (M^n)_{15} > 0, et nombre de chaînes = ce coefficient
    nmin = min(k for k in range(1, 6) if (G ** k)[0, 4] > 0)
    assert nmin == dist and sorted(chaines(G, dist, 1, 5)) == sorted(minimales)
    assert (G ** dist)[0, 4] == len(minimales)
    for dl, nl in leurres:                               # couples (distance, nombre) proposés à tort
        assert (dl, nl) != (dist, len(minimales))

# al-t-matrices-markov-ex6 (exercice, v1, v2) : matrice complétée, π1, π2, probabilité demandée et leurres
cas = [
    ([[R(5, 10), R(3, 10)], [R(4, 10), R(4, 10)], [R(6, 10)]], 0, 2, R(16, 100), [R(2, 10), R(4, 100), R(36, 100)],
     [R(49, 100), R(35, 100), R(16, 100)]),
    ([[R(6, 10), R(1, 10)], [R(5, 10), R(2, 10)], [R(7, 10)]], 0, 1, R(17, 100), [R(1, 10), R(1, 100), R(21, 100)],
     [R(62, 100), R(17, 100), R(21, 100)]),
    ([[R(7, 10), R(2, 10)], [R(3, 10), R(5, 10)], [R(5, 10)]], 1, 0, R(46, 100), [R(3, 10), R(9, 100), R(41, 100)],
     [R(46, 100), R(41, 100), R(13, 100)]),
]
for donnees, depart, cible, bonne, leurres, pi2_cite in cas:
    (p11, p12), (p21, p22), (p31,) = donnees
    P = M([[p11, p12, 1 - p11 - p12], [p21, p22, 1 - p21 - p22], [p31, 1 - p31, 0]])
    assert all(x >= 0 for x in P) and all(sum(P.row(i)) == 1 for i in range(3))
    pi0 = zeros(1, 3)
    pi0[depart] = 1
    pi1, pi2 = pi0 * P, pi0 * P ** 2
    assert pi1 == P.row(depart) and list(pi2) == pi2_cite and sum(pi2) == 1
    assert pi2[cible] == bonne
    for l in leurres:
        assert l != bonne
# compléments cités dans les corrections : 0,2 ; 0,2 ; 0,4 (exercice) etc.
assert 1 - R(5, 10) - R(3, 10) == R(2, 10) and 1 - R(7, 10) - R(2, 10) == R(1, 10) and 1 - R(6, 10) - R(1, 10) == R(3, 10)

# --- Relecture indépendante (2026-10-29) : réponses retrouvées à l'aveugle, recalculées ici.
from sympy import Matrix as _M, Rational as _R, eye as _eye, ones as _ones
from itertools import product as _prod
# al-t-matrices-q7 : M³ du graphe orienté, coefficient (1,1) = 2
_G = _M([[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [1, 0, 0, 0]])
assert _G ** 3 == _M([[2, 0, 0, 1], [1, 1, 1, 0], [0, 1, 1, 0], [0, 0, 1, 2]])
# al-t-matrices-inverse-ex5 Q5 (exercice, v1, v2) : résolution directe des trois systèmes
for _a, _Y, _sol in [(2, [4, 6, 10], [0, 1, 3]), (3, [6, 3, 9], [1, 0, 2]), (-1, [1, 0, 3], [1, 2, -1])]:
    _A = _a * _eye(3) + _ones(3, 3)
    assert _A.solve(_M(_Y)) == _M(_sol)
    # A^{-1} = (1/a) I - b/(a(a+3b)) J avec b = 1
    assert _A.inv() == _R(1, _a) * _eye(3) - _R(1, _a * (_a + 3)) * _ones(3, 3)
# al-t-matrices-adjacence-ex6 (exercice, v1, v2) : distance et nombre de chaînes minimales par énumération des suites de sommets
def _chaines(aretes, n):
    adj = {frozenset(e) for e in aretes}
    return [c for c in _prod(range(1, 6), repeat=n - 1) if all(frozenset(p) in adj for p in zip((1,) + c, c + (5,)))]
for _ar, _d, _nb in [([(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)], 3, 2), ([(1, 2), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)], 2, 3),
                     ([(1, 2), (2, 3), (2, 4), (3, 4), (4, 5)], 3, 1)]:
    assert all(len(_chaines(_ar, _k)) == 0 for _k in range(1, _d)) and len(_chaines(_ar, _d)) == _nb
# al-t-matrices-markov-ex6 (exercice, v1, v2) : π2 = π0 P²
for _P, _pi0, _pi2 in [([[5, 3, 2], [4, 4, 2], [6, 4, 0]], [1, 0, 0], [49, 35, 16]), ([[6, 1, 3], [5, 2, 3], [7, 3, 0]], [1, 0, 0], [62, 17, 21]),
                       ([[7, 2, 1], [3, 5, 2], [5, 5, 0]], [0, 1, 0], [46, 41, 13])]:
    assert _M([_pi0]) * (_M(_P) / 10) ** 2 == _M([_pi2]) / 100

print("enrichi_al_t_matrices : tout est vérifié.")
