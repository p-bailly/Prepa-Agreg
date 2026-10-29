# Vérification calculatoire de l'enrichissement du chapitre « al-l1-applications » (thème algebre-lineaire, niveau
# l1) : questions de content/qcm/al-l1-applications.json et exercices ajoutés à
# content/exercices/algebre-lineaire/l1-al-l1-applications.json (rang-ex6, isomorphismes-ex5, matrice-ex5,
# variantes « v1 », « v2 » comprises). Calcul exact avec SymPy. Chaque assert cite l'identifiant de la question ou de
# l'exercice. Voir CLAUDE.md, « Contrôle de l'exactitude ».
import sympy as sp
from sympy import Matrix as M, symbols, Poly, expand, Rational as R, zeros, eye, factor

X, a, N = symbols("X a N")


def coords(p, n):
    """Coordonnées de p dans (1, X, ..., X^n)."""
    pp = Poly(expand(p), X)
    return [pp.coeff_monomial(X ** k) for k in range(n + 1)]


def matrice(f, n):
    """Matrice de l'endomorphisme f de R_n[X] dans la base (1, X, ..., X^n)."""
    return M([coords(f(X ** j), n) for j in range(n + 1)]).T


# al-l1-applications-q1 : Ker(P -> XP' - P) = Vect(X) dans R_2[X]
f = lambda p: X * sp.diff(p, X) - p
A = matrice(f, 2)
ns = A.nullspace()
assert len(ns) == 1 and list(ns[0]) == [0, 1, 0]
assert f(1) == -1 and expand(f(X ** 2)) == X ** 2 and f(X) == 0

# al-l1-applications-q2 (faux) : p(x, y) = (x, 0) : Ker ∩ Im = {0} mais p non injectif
P = M([[1, 0], [0, 0]])
assert P.nullspace() == [M([0, 1])] and P.columnspace() == [M([1, 0])]
assert M.hstack(M([0, 1]), M([1, 0])).rank() == 2          # intersection réduite à {0}

# al-l1-applications-q3 : rang de P -> P'' sur R_3[X] = 2
A = matrice(lambda p: sp.diff(p, X, 2), 3)
assert A.rank() == 2 and len(A.nullspace()) == 2

# al-l1-applications-q4 (faux) : dim Ker + dim Im = 3 impair ; exemple en dimension 2
assert all(2 * k != 3 for k in range(4))
A2 = M([[0, 1], [0, 0]])
assert A2.nullspace() == [M([1, 0])] and A2.columnspace() == [M([1, 0])]

# al-l1-applications-q5 : matrice de P -> P(X+1) sur R_2[X] ; leurres faux
A = matrice(lambda p: p.subs(X, X + 1) if p != 1 else sp.Integer(1), 2)
assert A == M([[1, 1, 1], [0, 1, 2], [0, 0, 1]])
for l in ([[1, 0, 0], [1, 1, 0], [1, 2, 1]], [[1, 1, 1], [0, 1, 1], [0, 0, 1]], [[1, 1, 0], [0, 1, 2], [0, 0, 1]]):
    assert M(l) != A

# al-l1-applications-q6 : rang 2, L3 = L1 + L2
A = M([[1, 1, 1, 1], [1, 2, 3, 4], [2, 3, 4, 5]])
assert A.rank() == 2 and A.row(2) == A.row(0) + A.row(1)

# al-l1-applications-q7 (vrai) : A^2 = A, image Vect((2,1)), noyau Vect((1,1)), trace = rang = 1
A = M([[2, -2], [1, -1]])
assert A ** 2 == A and A.trace() == 1 == A.rank()
assert A * M([1, 1]) == zeros(2, 1) and A * M([2, 1]) == M([2, 1])

# al-l1-applications-q8 : système incompatible ; leurres faux (aucune solution)
x, y, z = symbols("x y z")
assert sp.linsolve([x + 2 * y - z - 1, 2 * x + 4 * y - 2 * z - 3], [x, y, z]) == sp.EmptySet

# al-l1-applications-rang-ex6 (exercice, v1, v2)
for n in range(1, 7):
    D = matrice(lambda p: expand(p.subs(X, X + 1) - p) if p != 1 else sp.Integer(0), n)
    assert D.rank() == n                                          # Q3 : rang n
    assert D.nullspace() == [M([1] + [0] * n)]                    # Q2 : noyau = constantes
    assert all(D[n, j] == 0 for j in range(n + 1))                # image dans R_{n-1}[X]
    for k in range(1, n + 1):                                     # Q1 : degré k -> degré k-1, coefficient k a_k
        d = Poly(expand((X + 1) ** k - X ** k), X)
        assert d.degree() == k - 1 and d.LC() == k
cas = [
    (X ** 2, 3, X * (X - 1) * (2 * X - 1) / 6, N * (N + 1) * (2 * N + 1) / 6, lambda k: k ** 2,
     [X * (X + 1) * (2 * X + 1) / 6, X ** 3 / 3, X * (X - 1) / 2]),
    (X ** 3, 4, X ** 2 * (X - 1) ** 2 / 4, N ** 2 * (N + 1) ** 2 / 4, lambda k: k ** 3,
     [X ** 2 * (X + 1) ** 2 / 4, X ** 4 / 4, X * (X - 1) * (2 * X - 1) / 6]),
    (X, 2, X * (X - 1) / 2, N * (N + 1) / 2, lambda k: k,
     [X * (X + 1) / 2, X ** 2 / 2, (X ** 2 - 1) / 2]),
]
for cible, n, Pb, somme, terme, leurres in cas:
    Delta = lambda p: expand(p.subs(X, X + 1) - p)
    assert Delta(Pb) == expand(cible) and Pb.subs(X, 0) == 0 and Poly(expand(Pb), X).degree() <= n
    # unicité : les solutions de Δ(P) = cible, P(0) = 0, dans R_n[X], par identification des coefficients
    cs = symbols(f"c0:{n + 1}")
    Pg = sum(c * X ** k for k, c in enumerate(cs))
    sol = sp.solve(coords(Delta(Pg) - cible, n) + [Pg.subs(X, 0)], cs, dict=True)
    assert len(sol) == 1 and expand(Pg.subs(sol[0]) - Pb) == 0
    for l in leurres:                                             # chaque leurre ne convient pas
        assert Delta(l) != expand(cible) or l.subs(X, 0) != 0
    # Q5 : somme télescopique
    assert expand(Pb.subs(X, N + 1) - somme) == 0
    for Nv in range(0, 12):
        assert sum(terme(k) for k in range(Nv + 1)) == somme.subs(N, Nv)
# coefficients cités : a = 1/3, b = -1/2, c = 1/6 ; a = 1/4, b = -1/2, c = 1/4, d = 0 ; a = 1/2, b = -1/2
assert expand(X ** 3 / 3 - X ** 2 / 2 + X / 6 - X * (X - 1) * (2 * X - 1) / 6) == 0
assert expand(X ** 4 / 4 - X ** 3 / 2 + X ** 2 / 4 - X ** 2 * (X - 1) ** 2 / 4) == 0

# al-l1-applications-isomorphismes-ex5 (exercice, v1, v2) : valeurs de a où f_a n'est pas bijectif, vecteur du noyau
cas = [
    (M([[1, a, 0], [0, 1, a], [a, 0, 1]]), {-1}, {-1: [1, 1, 1]}),
    (M([[a, 1, 0], [1, a, 0], [0, 0, 1]]), {-1, 1}, {1: [1, -1, 0], -1: [1, 1, 0]}),
    (M([[1, 1, 1], [1, a, 0], [1, 0, a]]), {0, 2}, {0: [0, 1, -1], 2: [2, -1, -1]}),
]
for A, interdites, noyau in cas:
    reels = {r for r in sp.solve(A.det(), a) if r.is_real}
    assert reels == interdites
    for av, v in noyau.items():
        assert A.subs(a, av) * M(v) == zeros(3, 1) and M(v) != zeros(3, 1)
    for av in (-3, -2, R(1, 2), 3, 5):
        if av not in interdites:
            assert A.subs(a, av).rank() == 3
# leurres : chaque ensemble proposé à tort diffère de l'ensemble exact
assert {1} != {-1} and {-1, 1} != {-1} and set() != {-1}
assert {-1} != {-1, 1} and {0} != {-1, 1} and {1} != {-1, 1}
assert {2} != {0, 2} and {0, -2} != {0, 2} and {0} != {0, 2}
assert sp.factor(1 + a ** 3) == (a + 1) * (a ** 2 - a + 1)

# al-l1-applications-matrice-ex5 (exercice, v1, v2) : matrices de f, g, g∘f = BA, f∘g = AB, bijectivité
cas = [
    (lambda x, y: (x + y, x - y, 2 * x), lambda x, y, z: (x - z, y + z),
     [[1, 1], [1, -1], [2, 0]], [[1, 0, -1], [0, 1, 1]], [[-1, 1], [3, -1]], [[1, 1, 0], [1, -1, -2], [2, 0, -2]], True,
     [[[-1, 3], [1, -1]], [[1, -1], [3, -1]], [[-1, 1], [3, 1]]]),
    (lambda x, y: (x, x + y, y), lambda x, y, z: (x + y + z, x - z),
     [[1, 0], [1, 1], [0, 1]], [[1, 1, 1], [1, 0, -1]], [[2, 2], [1, -1]], [[1, 1, 1], [2, 1, 0], [1, 0, -1]], True,
     [[[2, 1], [2, -1]], [[2, 2], [1, 1]], [[1, 2], [1, -1]]]),
    (lambda x, y: (x - y, y, x + y), lambda x, y, z: (x + z, y - z),
     [[1, -1], [0, 1], [1, 1]], [[1, 0, 1], [0, 1, -1]], [[2, 0], [-1, 0]], [[1, -1, 2], [0, 1, -1], [1, 1, 0]], False,
     [[[2, -1], [0, 0]], [[2, 0], [1, 0]], [[0, 2], [0, -1]]]),
]
for f, g, A, B, BA, AB, bij, leurres in cas:
    A, B = M(A), M(B)
    assert A == M.hstack(M(f(1, 0)), M(f(0, 1)))
    assert B == M.hstack(M(g(1, 0, 0)), M(g(0, 1, 0)), M(g(0, 0, 1)))
    assert B * A == M(BA) and A * B == M(AB)
    # g∘f calculée directement sur la base canonique
    assert M.hstack(M(g(*f(1, 0))), M(g(*f(0, 1)))) == M(BA)
    assert (M(BA).det() != 0) == bij
    assert (A * B).rank() <= 2 and (A * B).det() == 0
    for l in leurres:
        assert M(l) != M(BA)
# déterminants cités : -2, -4 ; vecteur (0, 1) du noyau dans la variante 2
assert M([[-1, 1], [3, -1]]).det() == -2 and M([[2, 2], [1, -1]]).det() == -4
assert M([[2, 0], [-1, 0]]) * M([0, 1]) == zeros(2, 1)

# --- Relecture indépendante (2026-10-29) : réponses retrouvées à l'aveugle, recalculées ici.
from sympy import Matrix as _M, symbols as _sym, expand as _expand, factor as _factor, summation as _sum, Rational as _R
_X, _N, _k, _a = _sym('X N k a')
# al-l1-applications-rang-ex6 (exercice, v1, v2) : Δ(P) = X², X³, X ; P(0) = 0 ; sommes de puissances ; leurres faux
for _P, _cible, _leurres in [(_X * (_X - 1) * (2 * _X - 1) / 6, _X ** 2, [_X * (_X + 1) * (2 * _X + 1) / 6, _X ** 3 / 3, _X * (_X - 1) / 2]),
                             (_X ** 2 * (_X - 1) ** 2 / 4, _X ** 3, [_X ** 2 * (_X + 1) ** 2 / 4, _X ** 4 / 4, _X * (_X - 1) * (2 * _X - 1) / 6]),
                             (_X * (_X - 1) / 2, _X, [_X * (_X + 1) / 2, _X ** 2 / 2, (_X ** 2 - 1) / 2])]:
    _D = lambda Q: _expand(Q.subs(_X, _X + 1) - Q)
    assert _D(_P) == _cible and _P.subs(_X, 0) == 0
    assert all(_D(L) != _cible or L.subs(_X, 0) != 0 for L in _leurres)
    assert _factor(_sum(_cible.subs(_X, _k), (_k, 0, _N)) - _P.subs(_X, _N + 1)) == 0
# al-l1-applications-isomorphismes-ex5 (exercice, v1, v2) : déterminants 1 + a³, a² - 1, a(a - 2) et vecteurs des noyaux
_mats = [_M([[1, _a, 0], [0, 1, _a], [_a, 0, 1]]), _M([[_a, 1, 0], [1, _a, 0], [0, 0, 1]]), _M([[1, 1, 1], [1, _a, 0], [1, 0, _a]])]
assert [_factor(m.det()) for m in _mats] == [_factor(1 + _a ** 3), _factor(_a ** 2 - 1), _factor(_a * (_a - 2))]
for _m, _val, _v in [(0, -1, (1, 1, 1)), (1, 1, (1, -1, 0)), (1, -1, (1, 1, 0)), (2, 0, (0, 1, -1)), (2, 2, (-2, 1, 1))]:
    assert _mats[_m].subs(_a, _val) * _M(_v) == _M([0, 0, 0])
# al-l1-applications-matrice-ex5 (exercice, v1, v2) : BA et AB recalculés à partir des formules de f et g
for _A, _B, _BA, _AB in [([[1, 1], [1, -1], [2, 0]], [[1, 0, -1], [0, 1, 1]], [[-1, 1], [3, -1]], [[1, 1, 0], [1, -1, -2], [2, 0, -2]]),
                         ([[1, 0], [1, 1], [0, 1]], [[1, 1, 1], [1, 0, -1]], [[2, 2], [1, -1]], [[1, 1, 1], [2, 1, 0], [1, 0, -1]]),
                         ([[1, -1], [0, 1], [1, 1]], [[1, 0, 1], [0, 1, -1]], [[2, 0], [-1, 0]], [[1, -1, 2], [0, 1, -1], [1, 1, 0]])]:
    assert _M(_B) * _M(_A) == _M(_BA) and _M(_A) * _M(_B) == _M(_AB) and (_M(_AB)).rank() <= 2

print("enrichi_al_l1_applications : tout est vérifié.")
