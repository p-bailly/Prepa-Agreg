"""Vérifications calculatoires du chapitre « al-l3-dualite » (algèbre linéaire, L3 : dualité et quotients) :
cartes, exercices et QCM de content/*/algebre-lineaire/l3-al-l3-dualite.json et content/qcm/al-l3-dualite.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/al_l3_dualite.py (doit finir sans erreur, en quelques secondes).
"""

from itertools import product

from sympy import (Matrix, Rational as R, diff, eye, expand, integrate, prod, symbols, zeros, Poly, factor)

X, t, x, y, z = symbols("X t x y z")


def lagrange(pts):
    """Polynômes de Lagrange associés à des points distincts."""
    return [expand(prod((X - b) / (a - b) for b in pts if b != a)) for a in pts]


def forme(coeffs, v):
    return sum(c * w for c, w in zip(coeffs, v))


def proportionnelles(u, v):
    return Matrix([list(u), list(v)]).rank() < 2


# ---------------------------------------------------------------------------
# al-l3-dualite-hyperplans (carte)
# ---------------------------------------------------------------------------
phi = (1, 2, -1)
# Exemple : v = (1,1,1) = (1,1,3) - 2 (0,0,1), avec (1,1,3) dans x + 2y - z = 0.
assert forme(phi, (1, 1, 1)) == 2 and forme(phi, (0, 0, 1)) == -1 and forme(phi, (1, 1, 3)) == 0
assert Matrix([1, 1, 3]) - 2 * Matrix([0, 0, 1]) == Matrix([1, 1, 1])
# Exemple : P = (P - P(1)) + P(1) ; P - P(1) s'annule en 1.
P = 3 * X**2 - X + 5
assert (P - P.subs(X, 1)).subs(X, 1) == 0
# Pièges 2 et 3 / contre-exemples : x et y sur R^2 ou R^3 non proportionnelles.
assert not proportionnelles((1, 0, 0), (0, 1, 0))
# Contre-exemple : g(x,y) = (y,x) n'est pas un multiple de id.
G = Matrix([[0, 1], [1, 0]])
assert all(G != c * eye(2) for c in (-1, 0, 1)) and G * Matrix([1, 0]) == Matrix([0, 1])
# Polynômes pairs : phi(X) X^3 - phi(X^3) X est impair (aucun monôme pair), donc nul s'il est dans P.
a1, a3 = symbols("a1 a3")
Q = expand(a1 * X**3 - a3 * X)
assert all(Poly(Q, X).coeff_monomial(X**k) == 0 for k in (0, 2, 4))
# Trou 1 : leurre « au moins n-1 » faux (E lui-même, de dimension n, n'est pas un hyperplan) ; n-2 : dimension d'une
# intersection de deux hyperplans distincts, pas d'un hyperplan. Trou 2 : phi et 2 phi ont le même noyau sans être égales,
# ni de somme nulle, et ne sont pas indépendantes.
assert proportionnelles((1, 2, -1), (2, 4, -2)) and (1, 2, -1) != (2, 4, -2)
# Affirmation 1 : (1,0) et (0,1) dans x = 0 ou y = 0, pas (1,1).
assert (1 * 0 == 0) and not (1 == 0 or 1 == 0)
# Affirmation 3 : dim(H1 ∩ H2) = n - 2 (exemple dans R^4 : x1 = 0 et x2 = 0).
assert len(Matrix([[1, 0, 0, 0], [0, 1, 0, 0]]).nullspace()) == 2

# hyperplans-ex1 (+ v1, v2) : v = h + lambda a.
for (ph, a, v, lam, h, leurres) in [
    ((2, -1, 3), (1, 0, 0), (1, 1, 1), 2, (-1, 1, 1), [(R(1, 2), (R(1, 2), 1, 1)), (4, (-3, 1, 1)), (-2, (3, 1, 1))]),
    ((1, 1, -2), (0, 0, 1), (3, 1, 0), -2, (3, 1, 2), [(2, (3, 1, -2)), (R(-1, 2), (3, 1, R(1, 2))), (4, (3, 1, -4))]),
    ((3, -2, 1), (0, 1, 0), (1, 2, 3), -1, (1, 3, 3), [(1, (1, 1, 3)), (2, (1, 0, 3)), (-2, (1, 4, 3))]),
]:
    assert forme(ph, a) != 0
    assert R(forme(ph, v), forme(ph, a)) == lam and forme(ph, h) == 0
    assert Matrix(h) + lam * Matrix(a) == Matrix(v)
    for (l2, h2) in leurres:  # décompositions cohérentes v = h + lambda a, mais h hors de H
        assert Matrix(h2) + l2 * Matrix(a) == Matrix(v) and forme(ph, h2) != 0

# hyperplans-ex2 (+ v1, v2) : équation de Vect(u, v, w) dans R^4.
for (vecs, eq, leurres) in [
    ([(1, 1, 0, 0), (0, 1, 1, 0), (0, 0, 1, 1)], (1, -1, 1, -1), [(1, 1, 1, 1), (1, -1, -1, 1), (1, 0, -1, 0)]),
    ([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)], (1, 1, 1, -1), [(1, 1, 1, 1), (1, -1, 1, -1), (1, 1, 1, -3)]),
    ([(1, 2, 0, 0), (0, 1, -1, 0), (0, 0, 1, 3)], (6, -3, -3, 1), [(2, -1, -1, 1), (6, -3, 3, 1), (6, 3, 3, -1)]),
]:
    M = Matrix(vecs)
    assert M.rank() == 3
    assert all(forme(eq, w) == 0 for w in vecs)
    for L in leurres:
        assert any(forme(L, w) != 0 for w in vecs)  # pas une équation de H

# hyperplans-ex3 : phi(AB) = phi(BA) => phi = lambda tr ; dim T = 1 ; Ker tr engendré par des commutateurs (n = 3).
n = 3
def E(i, j):
    M = zeros(n)
    M[i, j] = 1
    return M
for i, j, k, l in product(range(n), repeat=4):
    assert E(i, j) * E(k, l) == (E(i, l) if j == k else zeros(n))
A = Matrix(n, n, symbols("a0:9"))
B = Matrix(n, n, symbols("b0:9"))
assert expand((A * B).trace() - (B * A).trace()) == 0
coms = [E(i, i) * E(i, j) - E(i, j) * E(i, i) for i in range(n) for j in range(n) if i != j]
coms += [E(i, 0) * E(0, i) - E(0, i) * E(i, 0) for i in range(1, n)]
assert Matrix([list(c) for c in coms]).rank() == n * n - 1 and all(c.trace() == 0 for c in coms)
# Système linéaire : les formes phi telles que phi(E_ij E_kl) = phi(E_kl E_ij) forment un espace de dimension 1.
cs = symbols("c0:9")
def phi_of(M):
    return sum(cs[3 * i + j] * M[i, j] for i in range(3) for j in range(3))
eqs = [phi_of(E(i, j) * E(k, l) - E(k, l) * E(i, j)) for i, j, k, l in product(range(n), repeat=4)]
Msys = Matrix([[e.coeff(c) for c in cs] for e in eqs])
assert 9 - Msys.rank() == 1  # dim T = 1 (leurres 0, n, n^2 - 1 faux pour n >= 2)

# hyperplans-ex4 : v(t) = e1 + t e2 + ... ; phi(v(t)) polynôme de coefficients phi(e_i).
c = symbols("p1:4")
assert expand(forme(c, (1, t, t**2))) == c[0] + c[1] * t + c[2] * t**2

# hyperplans-ex5 : tr(A E_ij) = a_ji ; P = E_12 + E_23 + E_31 inversible ; I - lambda E_ij de déterminant 1.
for i, j in product(range(n), repeat=2):
    assert (A * E(i, j)).trace() == A[j, i]
Pc = E(0, 1) + E(1, 2) + E(2, 0)
assert Pc.det() != 0 and all(Pc[i, i] == 0 for i in range(n))
lam = symbols("lam")
assert (eye(n) - lam * E(0, 2)).det() == 1 and (eye(n) - lam * E(2, 0)).det() == 1
# n = 2 : la permutation E_12 + E_21.
assert Matrix([[0, 1], [1, 0]]).det() == -1

# ---------------------------------------------------------------------------
# al-l3-dualite-base-duale (carte)
# ---------------------------------------------------------------------------
# Exemple : base ((1,1),(1,2)), duale 2x - y et -x + y.
Pm = Matrix([[1, 1], [1, 2]])
assert Pm.inv() == Matrix([[2, -1], [-1, 1]])
assert Pm.inv() * Pm == eye(2)
# Cours / contre-exemple : e1* = x pour e2 = (0,1), mais x - y pour e2 = (1,1).
assert Matrix([[1, 0], [0, 1]]).inv()[0, :] == Matrix([[1, 0]])
assert Matrix([[1, 1], [0, 1]]).inv()[0, :] == Matrix([[1, -1]])
# Exemple Lagrange en 0, 1, 2.
L0, L1, L2 = lagrange([0, 1, 2])
assert L0 == expand((X - 1) * (X - 2) / 2) and L1 == expand(-X * (X - 2)) and L2 == expand(X * (X - 1) / 2)
assert expand(L0 + 3 * L1 + 7 * L2) == X**2 + X + 1
# Affirmation 4 : somme des L_i = 1 ; trou 3 : leurres faux (n = 2).
assert expand(L0 + L1 + L2) == 1
pts = [0, 1, 2]
for i, a in enumerate(pts):
    sans_den = expand(prod(X - b for b in pts if b != a))
    inverse = expand(prod((X - a) / (b - a) for b in pts if b != a))
    assert sans_den != lagrange(pts)[i] and inverse.subs(X, a) == 0 and X**i != lagrange(pts)[i]
# Trou 1 : 1 - delta_ij pour n = 2 donne la base (e2*, e1*), pas la base duale.
# Contre-exemple : X(X-1)(X-2) nul en 0, 1, 2 ; X(X-1) nul en 0, 0, 1.
assert all(expand(X * (X - 1) * (X - 2)).subs(X, a) == 0 for a in (0, 1, 2))

# base-duale-ex1 (+ v1, v2) : e2* = deuxième ligne de P^{-1}.
for (cols, e2, leurres) in [
    ([(1, 0, 0), (1, 1, 0), (1, 1, 1)], (0, 1, -1), [(0, 1, 0), (1, -1, 0), (0, -1, 1)]),
    ([(1, 1, 0), (0, 1, 1), (1, 0, 1)], (R(-1, 2), R(1, 2), R(1, 2)), [(-1, 1, 1), (R(1, 2), R(1, 2), R(-1, 2)), (0, 1, 0)]),
    ([(1, 2, 0), (0, 1, 0), (0, 3, 1)], (-2, 1, -3), [(0, 1, 0), (-2, 1, 0), (2, 1, 3)]),
]:
    Pm = Matrix(cols).T
    assert Pm.det() != 0
    assert Pm.inv()[1, :] == Matrix([e2])
    for L in leurres:
        assert Matrix([L]) != Pm.inv()[1, :]
assert Matrix([(1, 1, 0), (0, 1, 1), (1, 0, 1)]).T.det() == 2
assert Matrix([(1, 1, 0), (0, 1, 1), (1, 0, 1)]).T.inv() == R(1, 2) * Matrix([[1, 1, -1], [-1, 1, 1], [1, -1, 1]])

# base-duale-ex2 (+ v1, v2) : coordonnées d'une forme dans une base d'évaluations, et test sur X^3.
for (pts, fphi, bon, leurres, exact3) in [
    ([-1, 0, 1], lambda p: integrate(p, (X, -1, 1)), (R(1, 3), R(4, 3), R(1, 3)),
     [(R(2, 3),) * 3, (R(1, 6), R(2, 3), R(1, 6)), (R(1, 2), 1, R(1, 2))], True),
    ([0, R(1, 2), 1], lambda p: integrate(p, (X, 0, 1)), (R(1, 6), R(2, 3), R(1, 6)),
     [(R(1, 3),) * 3, (R(1, 3), R(4, 3), R(1, 3)), (R(1, 4), R(1, 2), R(1, 4))], True),
    ([0, 1, 2], lambda p: diff(p, X).subs(X, 0), (R(-3, 2), 2, R(-1, 2)),
     [(-1, 1, 0), (R(3, 2), -2, R(1, 2)), (R(-3, 2), 2, R(1, 2))], False),
]:
    Ls = lagrange(pts)
    assert tuple(fphi(Li) for Li in Ls) == bon
    for L in leurres:
        assert any(fphi(X**k) != sum(c * (X**k).subs(X, a) for c, a in zip(L, pts)) for k in range(3))
    assert (fphi(X**3) == sum(c * a**3 for c, a in zip(bon, pts))) == exact3
assert lagrange([0, R(1, 2), 1]) == [expand((2 * X - 1) * (X - 1)), expand(-4 * X * (X - 1)), expand(X * (2 * X - 1))]

# base-duale-ex3 (+ v1, v2) : base antéduale = colonnes de M^{-1}.
for (fs, minv, P1, leurres) in [
    ([lambda p: p.subs(X, 0), lambda p: diff(p, X).subs(X, 0), lambda p: integrate(p, (X, 0, 1))],
     Matrix([[1, 0, 0], [0, 1, 0], [-3, R(-3, 2), 3]]), 1 - 3 * X**2, [1, 1 - X**2, 1 - 2 * X]),
    ([lambda p: p.subs(X, 1), lambda p: diff(p, X).subs(X, 1), lambda p: p.subs(X, -1)],
     R(1, 4) * Matrix([[3, -2, 1], [2, 0, -2], [-1, 2, 1]]), (1 + X) * (3 - X) / 4, [(1 + X) / 2, (1 + X)**2 / 4, (1 + X) * (3 - X) / 2]),
    ([lambda p: p.subs(X, 0), lambda p: p.subs(X, 1), lambda p: diff(p, X, 2).subs(X, 0)],
     Matrix([[1, 0, 0], [-1, 1, R(-1, 2)], [0, 0, R(1, 2)]]), 1 - X, [1, (1 - X)**2, 1 - X**2]),
]:
    M = Matrix([[f(X**k) for k in range(3)] for f in fs])
    assert M.det() != 0 and M.inv() == minv
    ante = [expand(sum(minv[k, jj] * X**k for k in range(3))) for jj in range(3)]
    assert ante[0] == expand(P1)
    assert all(fs[i](ante[jj]) == (1 if i == jj else 0) for i in range(3) for jj in range(3))
    for L in leurres:
        assert [f(expand(L)) for f in fs] != [1, 0, 0]
# Valeurs citées dans les corrections.
assert Matrix([[1, 0, 0], [0, 1, 0], [1, R(1, 2), R(1, 3)]]).det() == R(1, 3)
assert Matrix([[1, 1, 1], [0, 1, 2], [1, -1, 1]]).det() == 4 and Matrix([[1, 0, 0], [1, 1, 1], [0, 0, 2]]).det() == 2

# base-duale-ex4 : formule interpolatoire ; somme des poids = b - a ; trapèzes ; contre-exemple de degré 2n + 2.
aa, bb = symbols("a b")
Lt = lagrange([aa, bb])
assert [factor(integrate(Li, (X, aa, bb))) for Li in Lt] == [factor((bb - aa) / 2)] * 2
ws = [integrate(Li, (X, -1, 1)) for Li in lagrange([-1, R(1, 3), 1])]
assert sum(ws) == 2  # b - a ; leurres (b-a)/2, (b-a)/(n+1), 0 faux
Qp = prod((X - yy)**2 for yy in (-1, R(1, 3), 1))
assert integrate(Qp, (X, -1, 1)) > 0 and Poly(Qp, X).degree() == 6

# base-duale-ex5 : Hermite. P = X + X^2 - X^3 ; leurres faux ; K_i = (X - a_i) L_i^2 (points 0, 1, 3).
Ph = X + X**2 - X**3
cond = lambda p: (p.subs(X, 0), diff(p, X).subs(X, 0), p.subs(X, 1), diff(p, X).subs(X, 1))
assert cond(Ph) == (0, 1, 1, 0)
for L in (3 * X**2 - 2 * X**3, X, X + 2 * X**2 - 2 * X**3):
    assert cond(L) != (0, 1, 1, 0)
pts = [0, 1, 3]
Ls = lagrange(pts)
for i, a in enumerate(pts):
    K = expand((X - a) * Ls[i]**2)
    assert Poly(K, X).degree() <= 5
    assert all(K.subs(X, b) == 0 for b in pts)
    assert all(diff(K, X).subs(X, b) == (1 if b == a else 0) for b in pts)
# Q2 : P -> (P(a_i), P'(a_i)) bijective sur R_5[X] (points 0, 1, 3).
Mh = Matrix([[(X**k).subs(X, a) for k in range(6)] for a in pts] + [[diff(X**k, X).subs(X, a) for k in range(6)] for a in pts])
assert Mh.det() != 0

# ---------------------------------------------------------------------------
# al-l3-dualite-orthogonalite (carte)
# ---------------------------------------------------------------------------
# Exemple : F = Vect((1,1,0,0),(0,1,1,1)) ; F^perp = Vect(x1-x2+x3, x1-x2+x4) ; x3 - x4 = phi1 - phi2.
Fv = Matrix([[1, 1, 0, 0], [0, 1, 1, 1]])
p1, p2 = Matrix([1, -1, 1, 0]), Matrix([1, -1, 0, 1])
assert Fv * p1 == zeros(2, 1) and Fv * p2 == zeros(2, 1) and len(Fv.nullspace()) == 2
assert p1 - p2 == Matrix([0, 0, 1, -1])
# (F^perp)° = F : les solutions du système sont exactement F.
assert Matrix.hstack(*Matrix([p1.T, p2.T]).nullspace()).rank() == 2
assert Matrix.vstack(Matrix([p1.T, p2.T]).nullspace()[0].T, Fv).rank() == 2
# Contre-exemple : formes liées x + y, 2x + 2y sur R^3 : intersection de dimension 2.
assert len(Matrix([[1, 1, 0], [2, 2, 0]]).nullspace()) == 2
# Contre-exemple (famille infinie) : P = X (X - 2)(X + 1)(X - 5) s'annule en 0, 2, -1, 5 et P'(0) = (-2)(1)(-5) = 10.
Pz = X * (X - 2) * (X + 1) * (X - 5)
assert all(Pz.subs(X, a) == 0 for a in (0, 2, -1, 5)) and diff(Pz, X).subs(X, 0) == 10
# Affirmation 3 : 2 équations dans R^5 laissent au moins 3 dimensions.
assert all(5 - Matrix([r1, r2]).rank() >= 3 for r1 in [(1, 0, 0, 0, 0)] for r2 in [(0, 1, 0, 0, 0), (1, 0, 0, 0, 0)])
# Trou 3 : leurre n - p faux pour une famille liée (n = 3, p = 2, r = 1).
assert 3 - Matrix([[1, 1, 0], [2, 2, 0]]).rank() == 2 != 3 - 2

# orthogonalite-ex1 (+ v1, v2) : équations de F dans R^4.
def definit(Fvecs, eqs):
    M = Matrix(eqs)
    return all(forme(e, w) == 0 for e in eqs for w in Fvecs) and 4 - M.rank() == Matrix(Fvecs).rank()
for (Fvecs, bon, leurres) in [
    ([(1, 0, 1, 0), (0, 1, 1, 1)], [(1, 1, -1, 0), (0, 1, 0, -1)],
     [[(1, 1, -1, 0), (2, 2, -2, 0)], [(1, 0, -1, 0), (0, 1, 0, -1)], [(1, 1, -1, 0), (1, 0, 0, 1)]]),
    ([(1, 1, 0, 0), (0, 0, 1, 1)], [(1, -1, 0, 0), (0, 0, 1, -1)],
     [[(1, -1, 0, 0), (0, 0, 1, 1)], [(1, -1, 0, 0), (-2, 2, 0, 0)], [(1, 1, -1, -1), (1, -1, 0, 0)]]),
    ([(1, 2, 0, 1), (0, 1, 1, 0)], [(1, 0, 0, -1), (0, 1, -1, -2)],
     [[(1, 0, 0, -1), (2, 0, 0, -2)], [(1, 0, 0, -1), (0, 1, -1, 0)], [(1, 0, 0, 1), (0, 1, -1, -2)]]),
]:
    assert definit(Fvecs, bon)
    assert not any(definit(Fvecs, L) for L in leurres)

# orthogonalite-ex2 (+ v1, v2) : rang des formes, dim F, psi combinaison.
for (forms, psi, coeffs, d, Fbase) in [
    ([(1, 1, 1), (1, -1, 0), (3, 1, 2)], (1, 3, 2), (2, -1), 1, [(1, 1, -2)]),
    ([(1, -1, 0, 0), (0, 1, -1, 0), (1, 0, -1, 0)], (1, 1, -2, 0), (1, 2), 2, [(1, 1, 1, 0), (0, 0, 0, 1)]),
    ([(1, 1, 1, 1), (1, -1, 0, 0), (0, 0, 1, -1), (2, 0, 1, 1)], (3, -1, 2, 0), (1, 2, 1), 1, [(1, 1, -1, -1)]),
]:
    M = Matrix(forms)
    nvar = M.shape[1]
    assert nvar - M.rank() == d and len(M.nullspace()) == d
    assert all(M * Matrix(f) == zeros(M.shape[0], 1) for f in Fbase) and Matrix(Fbase).rank() == d
    assert sum((ci * Matrix([forms[k]]) for k, ci in enumerate(coeffs)), zeros(1, nvar)) == Matrix([psi])
    assert all(forme(psi, f) == 0 for f in Fbase)
assert Matrix([(3, 1, 2)]) == 2 * Matrix([(1, 1, 1)]) + Matrix([(1, -1, 0)])
assert Matrix([(2, 0, 1, 1)]) == Matrix([(1, 1, 1, 1)]) + Matrix([(1, -1, 0, 0)])

# orthogonalite-ex3 (+ v1, v2) : dim(F^perp + G^perp) = n - dim(F ∩ G) si F + G = E ; exemple concret.
for (nn, df, dg, ans) in [(5, 3, 3, 4), (6, 4, 3, 5), (7, 5, 5, 4)]:
    di = df + dg - nn
    assert nn - di == ans == (nn - df) + (nn - dg)
    I_ = eye(nn)
    Fm = Matrix.hstack(*[I_[:, k] for k in range(df)])            # e1..e_df
    Gm = Matrix.hstack(*[I_[:, k] for k in range(nn - dg, nn)])   # derniers vecteurs
    Fp, Gp = Fm.T.nullspace(), Gm.T.nullspace()                    # F^perp, G^perp (en coordonnées duales)
    assert Matrix.hstack(Fm, Gm).rank() == nn and Matrix.hstack(*(Fp + Gp)).rank() == ans

# orthogonalite-ex4 : G° = {0} et psi = P'(0) hors de G (vu dans la carte) ; un polynôme non nul a un nombre fini de racines.
Pq = X * (X - 3) * (X - R(1, 2))
assert diff(Pq, X).subs(X, 0) == R(3, 2) and all(Pq.subs(X, a) == 0 for a in (0, 3, R(1, 2)))

# ---------------------------------------------------------------------------
# al-l3-dualite-transposee (carte)
# ---------------------------------------------------------------------------
# Exemple : u(x,y,z) = (x+2y, y-z) ; tu(f1*) = e1* + 2 e2*, tu(f2*) = e2* - e3* ; Ker u = Vect((-2,1,1)).
A = Matrix([[1, 2, 0], [0, 1, -1]])
assert A.T * Matrix([1, 0]) == Matrix([1, 2, 0]) and A.T * Matrix([0, 1]) == Matrix([0, 1, -1])
assert A * Matrix([-2, 1, 1]) == zeros(2, 1) and len(A.nullspace()) == 1
assert forme((1, 2, 0), (-2, 1, 1)) == 0 and forme((0, 1, -1), (-2, 1, 1)) == 0
# Contre-exemple : t(BA) = tA tB != tB tA.
Am, Bm = Matrix([[0, 1], [0, 0]]), Matrix([[0, 0], [1, 0]])
assert (Bm * Am).T == Am.T * Bm.T == Matrix([[0, 0], [0, 1]]) and Bm.T * Am.T == Matrix([[1, 0], [0, 0]])
# Contre-exemple : x -> (x, 0) ; tu envoie (a,b) -> b sur 0.
assert Matrix([[1], [0]]).T * Matrix([0, 1]) == zeros(1, 1)
# Affirmation 3 : même polynôme caractéristique.
Mt = Matrix([[1, 2, 3], [0, 4, 5], [7, 0, 1]])
assert Mt.charpoly(X) == Mt.T.charpoly(X)
# Rang des lignes = rang des colonnes (exemple).
assert Mt.rank() == Mt.T.rank()

# transposee-ex1 (+ v1, v2) : ligne de psi o u = (ligne de psi) A.
for (A, psi, res, leurres) in [
    (Matrix([[1, 1, 0], [0, 1, -2]]), (2, -1), (2, 1, 2), [(2, 1, -2), (2, -1, 0), (2, 3, 2)]),
    (Matrix([[1, 0, -1], [0, 2, 1]]), (1, 3), (1, 6, 2), [(1, 6, -4), (1, 3, 0), (1, 6, 4)]),
    (Matrix([[1, 1], [1, -1], [2, 0]]), (1, -1, 1), (2, 2), [(2, 0), (0, 2), (4, 2)]),
]:
    assert Matrix([psi]) * A == Matrix([res]) and A.T * Matrix(psi) == Matrix(res)
    assert all(Matrix([L]) != Matrix([res]) for L in leurres)
u1 = lambda x_, y_, z_: (x_ + y_, y_ - 2 * z_)
assert expand(2 * u1(x, y, z)[0] - u1(x, y, z)[1]) == 2 * x + y + 2 * z

# transposee-ex2 (+ v1, v2) : tD(delta_a) dans la base (delta_-1, delta_0, delta_1).
Ls = lagrange([-1, 0, 1])
for (fphi, bon, leurres) in [
    (lambda p: diff(p, X).subs(X, 0), (R(-1, 2), 0, R(1, 2)), [(-1, 0, 1), (R(1, 2), 0, R(-1, 2)), (0, -1, 1)]),
    (lambda p: diff(p, X).subs(X, 1), (R(1, 2), -2, R(3, 2)), [(R(-1, 2), 0, R(1, 2)), (0, -1, 1), (R(3, 2), -2, R(1, 2))]),
    (lambda p: diff(p, X, 2).subs(X, 0), (1, -2, 1), [(R(1, 2), -1, R(1, 2)), (0, -2, 1), (1, 0, 1)]),
]:
    assert tuple(fphi(Li) for Li in Ls) == bon
    for L in leurres:
        assert any(fphi(X**k) != sum(c * (X**k).subs(X, a) for c, a in zip(L, [-1, 0, 1])) for k in range(3))
assert Ls == [expand(X * (X - 1) / 2), 1 - X**2, expand(X * (X + 1) / 2)]

# transposee-ex3 (+ v1, v2) : plans stables = noyaux des vecteurs propres de tA (valeurs propres distinctes).
def plans_stables(A):
    res = []
    for lam_, mult, vs in A.T.eigenvects():
        assert mult == 1 and len(vs) == 1
        res.append(vs[0])
    return res
def stable(A, L):
    # le plan d'équation L.X = 0 est stable par A ssi L A est proportionnelle à L
    return Matrix.vstack(Matrix([L]), Matrix([L]) * A).rank() == 1
for (A, bons, leurres) in [
    (Matrix([[1, 1, 0], [0, 2, 1], [0, 0, 3]]), [(0, 0, 1), (0, 1, -1), (2, -2, 1)],
     [[(1, 0, 0), (0, 1, 0), (0, 0, 1)], [(1, 1, 0), (0, 1, 1), (0, 0, 1)], [(0, 0, 1)]]),
    (Matrix([[1, 0, 0], [1, 2, 0], [1, 1, 3]]), [(1, 0, 0), (1, 1, 0), (1, 1, 1)],
     [[(1, 0, 0), (0, 1, 0), (0, 0, 1)], [(1, 0, 0)], [(1, 0, 0), (1, -1, 0), (1, 1, 1)]]),
    (Matrix([[3, 1, 1], [0, 2, 1], [0, 0, 1]]), [(0, 0, 1), (0, 1, 1), (1, 1, 1)],
     [[(1, 0, 0), (0, 1, 0), (0, 0, 1)], [(0, 0, 1)], [(0, 0, 1), (0, 1, -1), (1, 1, 1)]]),
]:
    vps = plans_stables(A)
    assert len(vps) == 3 and all(stable(A, L) for L in bons)
    assert all(any(proportionnelles(L, v) for v in vps) for L in bons)
    for liste in leurres:  # un plan non stable, ou des plans stables manquants
        assert len(liste) < 3 or not all(stable(A, L) for L in liste)

# transposee-ex4 (+ v1, v2) : AX = B compatible ssi Y0 B = 0.
b1, b2, b3 = symbols("b1 b2 b3")
for (A, y0, leurres) in [
    (Matrix([[1, 1, 0], [0, 1, 1], [1, 2, 1]]), (1, 1, -1), [(1, -1, 1), (1, 1, 1), (0, 0, 1)]),
    (Matrix([[1, 0, 2], [0, 1, 1], [2, 1, 5]]), (2, 1, -1), [(1, 2, -1), (2, -1, 1), (1, 1, 1)]),
    (Matrix([[1, -1, 0], [2, 1, 3], [0, 3, 3]]), (2, -1, 1), [(2, 1, -1), (1, -2, 1), (0, 1, -1)]),
]:
    assert A.rank() == 2 and Matrix([y0]) * A == zeros(1, 3)
    ns = A.T.nullspace()
    assert len(ns) == 1 and proportionnelles(ns[0], y0)
    for L in leurres:
        assert not proportionnelles(L, y0)
    # B = colonne de A : toujours compatible ; la condition doit l'accepter.
    for k in range(3):
        assert forme(y0, A[:, k]) == 0

# transposee-ex5 : J nilpotent, un seul plan stable z = 0.
J = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
ev = J.T.eigenvects()
assert len(ev) == 1 and ev[0][0] == 0 and len(ev[0][2]) == 1 and proportionnelles(ev[0][2][0], (0, 0, 1))
assert stable(J, (0, 0, 1)) and not stable(J, (1, 0, 0)) and not stable(J, (0, 1, 0))
assert J.charpoly(X).as_expr() == X**3

# ---------------------------------------------------------------------------
# al-l3-dualite-quotients (carte)
# ---------------------------------------------------------------------------
# Exemple : u(x,y,z) = (x-y, y-z), noyau Vect((1,1,1)), u(a+b, b, 0) = (a, b) ; (1,0,0) et (2,1,1) même image.
Uq = Matrix([[1, -1, 0], [0, 1, -1]])
assert Uq.nullspace()[0] == Matrix([1, 1, 1]) and Uq.rank() == 2
aa_, bb_ = symbols("aa bb")
assert Uq * Matrix([aa_ + bb_, bb_, 0]) == Matrix([aa_, bb_])
assert Uq * Matrix([1, 0, 0]) == Uq * Matrix([2, 1, 1])
# Contre-exemple : F = {(0,0),(1,0)} : relation non symétrique.
Fpart = {(0, 0), (1, 0)}
assert (1 - 0, 0 - 0) in Fpart and (0 - 1, 0 - 0) not in Fpart

# quotients-ex1 (+ v1, v2) : phi de noyau F ; classe de v.
for (f1, f2, ph, v, bon, autres) in [
    ((1, 1, 1), (1, 0, -1), (1, -2, 1), (1, 0, 0), (0, 0, 1), [(0, 1, 0), (1, 1, 0), (2, 0, 0)]),
    ((1, 0, 1), (0, 1, 1), (1, 1, -1), (1, 0, 0), (0, 1, 0), [(0, 0, 1), (1, 1, 0), (1, 0, 1)]),
    ((1, 2, 0), (0, 1, 1), (2, -1, 1), (0, 1, 0), (0, 0, -1), [(0, 0, 1), (-1, 0, 0), (1, 0, 0)]),
]:
    assert forme(ph, f1) == 0 and forme(ph, f2) == 0 and Matrix([f1, f2]).rank() == 2
    # même classe <=> v - w dans F <=> v - w combinaison de f1, f2
    dans_F = lambda w: Matrix([f1, f2, tuple(a - b for a, b in zip(v, w))]).rank() == 2
    assert dans_F(bon) and forme(ph, bon) == forme(ph, v)
    assert not any(dans_F(w) for w in autres)

# quotients-ex2 (+ v1, v2) : E/F isomorphe à R^k.
for (deg, fs, k) in [
    (5, [lambda p: p.subs(X, 1), lambda p: diff(p, X).subs(X, 1)], 2),
    (4, [lambda p: p.subs(X, -1), lambda p: p.subs(X, 0), lambda p: p.subs(X, 1)], 3),
    (6, [lambda p, m=m: diff(p, X, m).subs(X, 0) for m in range(4)], 4),
]:
    M = Matrix([[f(X**j) for j in range(deg + 1)] for f in fs])
    assert M.rank() == k  # u surjective ; dim E/F = k ; dim F = deg + 1 - k
    assert len(M.nullspace()) == deg + 1 - k
assert [f(1 + 2 * (X - 1)) for f in [lambda p: p.subs(X, 1), lambda p: diff(p, X).subs(X, 1)]] == [1, 2]
a_, b_, c_, d_ = symbols("a_ b_ c_ d_")
Pt = a_ + b_ * X + c_ / 2 * X**2 + d_ / 6 * X**3
assert [diff(Pt, X, m).subs(X, 0) for m in range(4)] == [a_, b_, c_, d_]

# quotients-ex3 (+ v1, v2) : dim (E/F)/(G/F) = dim E - dim G.
for (dE, dG, dF, ans, leurres) in [(7, 5, 2, 2, [3, 5, 0]), (9, 6, 4, 3, [2, 5, 0]), (10, 4, 1, 6, [3, 9, 0])]:
    assert (dE - dF) - (dG - dF) == dE - dG == ans and ans not in leurres
    assert leurres[:2] == [dG - dF, dE - dF]

# quotients-ex4 (+ v1, v2) : chi de l'endomorphisme induit sur E/Vect(e1).
for (M, chi, leurres) in [
    (Matrix([[2, 1, 3], [0, 1, 4], [0, -1, 5]]), (X - 3)**2, [(X - 2) * (X - 3)**2, (X - 1) * (X - 5), X**2 - 6 * X + 13]),
    (Matrix([[1, 2, 0], [0, 3, 1], [0, 1, 3]]), (X - 2) * (X - 4), [(X - 1) * (X - 2) * (X - 4), (X - 3)**2, X**2 - 6 * X + 10]),
    (Matrix([[0, 1, 1], [0, 1, -2], [0, 1, 4]]), (X - 2) * (X - 3), [X * (X - 2) * (X - 3), (X - 1) * (X - 4), (X - 2) * (X + 3)]),
]:
    assert M[1:, 0] == zeros(2, 1)  # F = Vect(e1) stable
    Cb = M[1:, 1:]
    assert expand(Cb.charpoly(X).as_expr() - chi) == 0
    assert expand(M.charpoly(X).as_expr() - (X - M[0, 0]) * chi) == 0
    assert all(expand(L - chi) != 0 for L in leurres)
# Nilpotence : u_F et u barre nilpotents => u nilpotent (exemple triangulaire strict).
Nn = Matrix([[0, 5, 7], [0, 0, 2], [0, 0, 0]])
assert Nn**3 == zeros(3)

# quotients-ex5 : delta_0, delta_1, intégrale sur [0,1] libres (rang 3) : codimension 3.
Mc = Matrix([[1, 0, 0], [1, 1, 1], [1, R(1, 2), R(1, 3)]])
assert Mc.det() == R(-1, 6) and Mc.rank() == 3

# ---------------------------------------------------------------------------
# QCM
# ---------------------------------------------------------------------------
# q1 : dans R^5, deux hyperplans distincts : intersection de dimension 3.
assert len(Matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]]).nullspace()) == 3
# q2 : évaluations en points distincts : matrice de Vandermonde inversible.
assert Matrix([[a**k for k in range(4)] for a in (-2, 0, 1, 5)]).det() != 0
# q3 : L1 + 2 L2 = X (points 0, 1, 2) ; leurres faux.
L0, L1, L2 = lagrange([0, 1, 2])
assert expand(L1 + 2 * L2) == X and all(expand(L1 + 2 * L2 - w) != 0 for w in (X**2, 1, 3 * X))
# q4 : deux droites distinctes de R^2.
assert Matrix([[1, 0], [0, 1]]).rank() == 2  # F + G = R^2, F ∩ G = {0}
# q5 : quatre formes de rang 3 dans un espace de dimension 7.
F4 = Matrix([[1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0, 0]])
assert F4.rank() == 3 and len(F4.nullspace()) == 4
# q6 : u surjective <=> tu injective (exemple : A de rang 2 dans M_{2,3}).
assert Matrix([[1, 0, 1], [0, 1, 1]]).T.nullspace() == []
# q7 : A dans M_{3,5} de rang 2 : dim Ker tu = 1.
A7 = Matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [1, 1, 0, 0, 0]])
assert A7.rank() == 2 and len(A7.T.nullspace()) == 1 and len(A7.nullspace()) == 3
# q8 : les X^(2k+1) sont libres et forment une famille infinie dans le supplémentaire impair.
assert Matrix([[1 if j == 2 * k + 1 else 0 for j in range(12)] for k in range(6)]).rank() == 6

# ---------------------------------------------------------------------------
# Relecture indépendante : leurres des textes à trous (énoncés quantifiés « pour tout »), pièges.
# ---------------------------------------------------------------------------
# al-l3-dualite-base-duale, trou 1 (n >= 2) : les formes (x+y, x+y) vérifient e_i*(e_i) = 1 pour la base canonique de R^2
# sans être la base duale (famille liée) ; e_i*(e_j) = 0 (i != j) est vérifié par (2x, 3y), qui n'est pas la base duale.
assert all(forme((1, 1), e) == 1 for e in [(1, 0), (0, 1)]) and Matrix([[1, 1], [1, 1]]).rank() == 1
assert forme((2, 0), (0, 1)) == 0 and forme((0, 3), (1, 0)) == 0 and forme((2, 0), (1, 0)) != 1
# al-l3-dualite-orthogonalite, trou 1 (n = 2) : F = {0} donne dim F^perp = 2 != dim F = 0 ; F = E donne 0 != n = 2
# et 0 != n + dim F = 4. Trou 3 : trois formes nulles sur R^2 (p = 3, r = 0) : noyaux d'intersection R^2, de dimension
# 2 = n - r, alors que n - p = -1, r = 0 et p - r = 3.
assert len(zeros(1, 2).nullspace()) == 2 and len(eye(2).nullspace()) == 0
Z3 = zeros(3, 2)
assert Z3.rank() == 0 and len(Z3.nullspace()) == 2 and 2 not in (2 - 3, 0, 3 - 0)
# al-l3-dualite-transposee, trous 2 et 3 : u = 0 de R^2 dans R^3 : Ker tu = (R^3)* (dimension 3, donc pas {0}), rg tu = 0,
# alors que dim E - rg u = 2, dim F - rg u = 3, dim F = 3.
U0 = zeros(3, 2)
assert len(U0.T.nullspace()) == 3 and U0.T.rank() == 0 and 0 not in (2 - 0, 3 - 0, 3)
# al-l3-dualite-transposee, piège 2 : B' = B, C' = (2 f_1, 2 f_2) ; Mat_{B,C'}(u) = A/2, et la matrice de tu dans
# (C'*, B*) est t(A/2), différente de tA dès que A != 0 (A de l'exemple de la carte).
A = Matrix([[1, 2, 0], [0, 1, -1]])
Pc2 = 2 * eye(2)                        # colonnes : coordonnées de C' dans C
assert (Pc2.inv() * A).T != A.T and (Pc2.inv() * A).T == A.T / 2
# Piège 1 (contre-exemple de la carte) : u = id de R, base chi = 2x de E*, psi = y de F* : tu(psi) = psi = chi / 2.
assert R(1, 2) * 2 == 1
# al-l3-dualite-quotients, trou 3 (E = R^2) : F = E donne dim E/F = 0, alors que dim E = 2 et dim E + dim F = 4 ;
# F = {0} donne dim E/F = 2 != dim F = 0.
assert (2 - 2, 2 - 0) == (0, 2) and 0 not in (2, 4) and 2 != 0
# al-l3-dualite-hyperplans-ex5 : pour n = 1, M_1(K) = K, son seul hyperplan est {0}, sans matrice inversible ;
# pour n = 2, l'hyperplan des matrices de trace nulle contient E_12 + E_21, inversible.
assert Matrix([[0, 1], [1, 0]]).trace() == 0 and Matrix([[0, 1], [1, 0]]).det() != 0
# hyperplans-ex5 Q4 (n = 3) : I - lambda E_ij est dans H = {tr(AM) = 0} pour lambda = tr(A) / a_ji (a_ji != 0).
A5 = Matrix([[1, 2, 0], [3, 4, 5], [6, 7, 8]])
lam5 = A5.trace() / A5[2, 0]            # i = 0, j = 2 : tau_A(E_02) = a_20
M5 = eye(3) - lam5 * E(0, 2)
assert (A5 * M5).trace() == 0 and M5.det() == 1
# base-duale-ex4 Q3 (n = 1, [0,1], points 0 et 1/3) : Q = X^2 (X - 1/3)^2 d'intégrale > 0, nulle aux deux points.
Q3 = X**2 * (X - R(1, 3))**2
assert integrate(Q3, (X, 0, 1)) > 0 and Q3.subs(X, 0) == 0 and Q3.subs(X, R(1, 3)) == 0
# orthogonalite-ex4 Q3, cas où 0 n'est pas parmi les a_k : P = X (X - 1)(X + 2), nul en 1 et -2, P'(0) = -2 != 0.
P4 = X * (X - 1) * (X + 2)
assert P4.subs(X, 1) == 0 and P4.subs(X, -2) == 0 and diff(P4, X).subs(X, 0) == -2
# quotients-ex5 Q4 : les trois formes restent libres sur R[X] (déjà sur R_2[X]) ; F contient X(X-1)(5X^2 - 5X + 1)
# (intégrale nulle), donc F est de dimension infinie mais de codimension 3.
P5 = X * (X - 1) * (5 * X**2 - 5 * X + 1)
assert P5.subs(X, 0) == 0 and P5.subs(X, 1) == 0 and integrate(P5, (X, 0, 1)) == 0

print("al_l3_dualite : toutes les vérifications passent.")
