"""Vérifications calculatoires du chapitre « al-l1-espaces » (thème algebre-lineaire, niveau l1).

Cartes content/cartes/algebre-lineaire/l1-al-l1-espaces.json et exercices
content/exercices/algebre-lineaire/l1-al-l1-espaces.json. Chaque assert cite l'identifiant de la carte ou de
l'exercice vérifié (variante comprise : « v1 », « v2 »). Les leurres des questions à choix sont vérifiés FAUX.
Lancer : python3 verification/al_l1_espaces.py (doit finir sans erreur).
"""

from sympy import (Matrix, Rational, binomial, cosh, diff, exp, expand, factor, simplify, sinh, solve, symbols,
                   zeros)

x, y, z, X, m, a, t, n = symbols("x y z X m a t n")


def rang(*vecteurs):
    return Matrix([list(v) for v in vecteurs]).rank()


def libre(*vecteurs):
    return rang(*vecteurs) == len(vecteurs)


def dans(v, *gen):
    """v appartient à Vect(gen)."""
    return rang(*gen) == rang(*gen, v)


def dim_inter(F, G):
    """dim(F ∩ G) par la formule de Grassmann, F et G donnés par des familles génératrices."""
    return rang(*F) + rang(*G) - rang(*F, *G)


def noyau(*equations):
    """Dimension de l'espace des solutions d'un système homogène (lignes = coefficients)."""
    M = Matrix([list(e) for e in equations])
    return M.cols - M.rank()


# ---------------------------------------------------------------------------
# Cartes
# ---------------------------------------------------------------------------
# al-l1-espaces-sev : contre-exemples (axes : (1,0)+(0,1) hors de xy = 0) ; exemple x + 2y - z = 0 linéaire.
assert 1 * 1 != 0
lam, x1, y1, z1, x2, y2, z2 = symbols("lam x1 y1 z1 x2 y2 z2")
assert expand((lam * x1 + x2) + 2 * (lam * y1 + y2) - (lam * z1 + z2) - (lam * (x1 + 2 * y1 - z1) + (x2 + 2 * y2 - z2))) == 0

# al-l1-espaces-familles : exemple, base ((1,1,0),(1,0,1),(0,1,1)) et coordonnées (0, 2, 3) de (2, 3, 5).
B = Matrix([[1, 1, 0], [1, 0, 1], [0, 1, 1]]).T
assert B.det() == -2
assert B.solve(Matrix([2, 3, 5])) == Matrix([0, 2, 3])
X_, Y_, Z_ = symbols("X_ Y_ Z_")
assert B.solve(Matrix([X_, Y_, Z_])) == Matrix([(X_ + Y_ - Z_) / 2, (X_ - Y_ + Z_) / 2, (-X_ + Y_ + Z_) / 2])
# familles : contre-exemples ; (1,0),(0,1),(1,1) liée ; (e1,e2,e3,x) liée avec e3 = e1 + e2 ; (u+v, u-v) libre.
assert not libre((1, 0), (0, 1), (1, 1))
assert not libre((1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1)) and not dans((0, 0, 1), (1, 0, 0), (0, 1, 0), (1, 1, 0))
assert Matrix([[1, 1], [1, -1]]).det() != 0
# familles, affirmation 4 : ((1,0),(2,0),(0,1)) génératrice liée, ((1,0),(2,0)) non génératrice.
assert rang((1, 0), (2, 0), (0, 1)) == 2 and rang((1, 0), (2, 0)) == 1

# al-l1-espaces-dimension : exemple (u, e1, e3) base, e2 dans Vect(u, e1).
assert libre((1, 1, 0), (1, 0, 0), (0, 0, 1)) and dans((0, 1, 0), (1, 1, 0), (1, 0, 0))
# dimension, affirmation 4 : 4 vecteurs du plan z = 0 deux à deux non colinéaires, rang 2.
V4 = [(1, 0, 0), (0, 1, 0), (1, 1, 0), (1, -1, 0)]
assert rang(*V4) == 2 and all(rang(V4[i], V4[j]) == 2 for i in range(4) for j in range(i + 1, 4))

# al-l1-espaces-somme-directe : exemple D = Vect((1,1,1)) et x + y - z = 0 ; contre-exemples.
assert 1 + 1 - 1 != 0
assert rang((1, 0), (0, 1), (1, 1)) == 2  # 2 ≠ 1 + 1 + 1
assert dim_inter([(1, 0, 0), (0, 1, 0)], [(0, 1, 0), (0, 0, 1)]) == 1

# al-l1-espaces-grassmann : exemple P : x + y + z = 0 et Vect((1,0,0)), Vect((0,1,0)) supplémentaires.
P = [(1, -1, 0), (1, 0, -1)]
assert rang(*P, (1, 0, 0)) == 3 and rang(*P, (0, 1, 0)) == 3
# grassmann : trois droites du plan, 2 ≠ 3 ; non-unicité ; F = G droite.
assert rang((1, 0), (0, 1)) == 2 and rang((1, 0), (1, 1)) == 2
assert dim_inter([(1, 0, 0)], [(0, 1, 0)]) == 0 and rang((1, 0, 0), (0, 1, 0)) == 2
# grassmann, démonstration : S_t = Vect(s1 + t f1, ...) supplémentaire, exemple numérique dans R^2.
for tt in range(-3, 4):
    assert libre((1, 0), (0 + tt, 1))

# al-l1-espaces-inclusion : exemple P : x + y + z = 0, u = (1,1,-2), v = (2,-1,-1).
for v in [(1, 1, -2), (2, -1, -1), (1, -1, 0), (1, 0, -1)]:
    assert sum(v) == 0
assert libre((1, 1, -2), (2, -1, -1)) and libre((1, -1, 0), (1, 0, -1))
assert dim_inter([(1, 0, 0), (0, 1, 0)], [(0, 1, 0), (0, 0, 1)]) == 1  # plans z = 0 et x = 0

# ---------------------------------------------------------------------------
# Exercices
# ---------------------------------------------------------------------------
# sev-ex1 : C = {xy = 0}, D = {x >= 0} ; v1 : x^2 = y^2, y = x^2, x + y dans Z ; v2 : 2X - 3 en 1 et 2.
assert (1 + 1) ** 2 != (1 - 1) ** 2  # v1 : (1,1)+(1,-1) = (2,0) hors de A
assert 2 != 2 ** 2  # v1 : 2(1,1) hors de B
assert Rational(1, 2) != int(Rational(1, 2))  # v1 : 1/2 non entier
Pp = 2 * X - 3
assert Pp.subs(X, 1) * Pp.subs(X, 2) == -1  # v2
Q = symbols("q0:4")
PQ = sum(Q[i] * X ** i for i in range(4))
assert expand(diff(PQ, X).subs(X, 1)) == expand(sum(i * Q[i] for i in range(4)))  # v2 : D linéaire

# sev-ex2 : a^2 - 4 = 0 ; v1 : a^2 - 3a = 0 ; v2 : a^3 - a = 0 ; leurres faux.
assert set(solve(a ** 2 - 4, a)) == {-2, 2}
assert set(solve(a ** 2 - 3 * a, a)) == {0, 3}
assert set(solve(a ** 3 - a, a)) == {-1, 0, 1}
assert {2} != {-2, 2} and {0} != {-2, 2} and {-3, 0} != {0, 3} and {-1, 1} != {-1, 0, 1}

# sev-ex3 : u = (0,1,1,...) croissante, v = (0,0,-1,...) décroissante, somme (0,1,0,...) non monotone.
u = [0] + [1] * 5
w = [0, 0] + [-1] * 4
s = [u[i] + w[i] for i in range(6)]
assert s[:3] == [0, 1, 0]
# sev-ex3 v2 : E12, E21 de carré nul, (E12 + E21)^2 = I.
E12 = Matrix([[0, 1], [0, 0]])
E21 = Matrix([[0, 0], [1, 0]])
assert E12 ** 2 == zeros(2) and E21 ** 2 == zeros(2) and (E12 + E21) ** 2 == Matrix.eye(2)

# familles-ex1 : u3 = u2 - u1 ; leurres faux. v1 : libre. v2 : u3 = 2u1 + 3u2.
u1, u2, u3 = Matrix([1, 2, 1]), Matrix([2, 1, 0]), Matrix([1, -1, -1])
assert u3 == u2 - u1 and u3 != u1 - u2 and u3 != u1 + u2 and rang(u1, u2, u3) == 2
u1, u2, u3 = Matrix([1, 0, 2]), Matrix([0, 1, 1]), Matrix([1, 1, 4])
assert libre(u1, u2, u3)
u1, u2, u3 = Matrix([1, 1, 0]), Matrix([0, 1, 1]), Matrix([2, 5, 3])
assert u3 == 2 * u1 + 3 * u2 and u3 != 3 * u1 + 2 * u2 and u3 != 2 * u1 - 3 * u2

# familles-ex2 : coordonnées et leurres.
def coords(base, v):
    return list(Matrix([list(b) for b in base]).T.solve(Matrix(v)))


assert coords([(1, 1, 1), (0, 1, 1), (0, 0, 1)], (3, 5, 4)) == [3, 2, -1]
assert coords([(1, 0, 0), (1, 1, 0), (1, 1, 1)], (2, -1, 3)) == [3, -4, 3]  # v1
# familles-ex2 v2 : P = X^2 + 3X + 1 dans (1, X-1, (X-1)^2) ; leurres faux.
Pv = X ** 2 + 3 * X + 1
base_pol = [1, X - 1, (X - 1) ** 2]
assert expand(5 + 5 * (X - 1) + (X - 1) ** 2 - Pv) == 0
for c in [(1, 3, 1), (1, 5, 5), (5, 3, 1)]:
    assert expand(sum(ci * bi for ci, bi in zip(c, base_pol)) - Pv) != 0
assert coords([(1, 1, 1), (0, 1, 1), (0, 0, 1)], (3, 5, 4)) not in ([3, 5, 4], [3, 2, 1], [-1, 2, 3])
assert coords([(1, 0, 0), (1, 1, 0), (1, 1, 1)], (2, -1, 3)) not in ([2, -1, 3], [3, -4, -3], [3, 4, 3])

# familles-ex3 : det = -(m-1)^2 (m+2) ; v1 : det = 1 + m^3 ; v2 : det = m - 5.
d0 = Matrix([[1, 1, m], [1, m, 1], [m, 1, 1]]).det()
assert factor(d0) == factor(-(m - 1) ** 2 * (m + 2))
assert set(solve(d0, m)) == {1, -2}
d1 = Matrix([[1, m, 0], [0, 1, m], [m, 0, 1]]).det()
assert expand(d1 - (1 + m ** 3)) == 0 and [r for r in solve(d1, m) if r.is_real] == [-1]
d2 = Matrix([[1, 1, 1], [1, 2, 3], [1, 3, m]]).det()
assert expand(d2 - (m - 5)) == 0
# relations dans les cas liés
assert Matrix([1, 1, -2]) + Matrix([1, -2, 1]) + Matrix([-2, 1, 1]) == zeros(3, 1)  # m = -2
assert Matrix([1, -1, 0]) + Matrix([0, 1, -1]) + Matrix([-1, 0, 1]) == zeros(3, 1)  # v1, m = -1
assert Matrix([1, 1, 1]) - 2 * Matrix([1, 2, 3]) + Matrix([1, 3, 5]) == zeros(3, 1)  # v2, m = 5
# leurres : pour ces valeurs le déterminant ne s'annule pas (m = -1, 2 pour ex3 ; m = 1, 0 pour v1 ; m = 4, 1, 3 pour v2)
assert all(d0.subs(m, val) != 0 for val in (-1, 2))
assert all(d1.subs(m, val) != 0 for val in (1, 0))
assert all(d2.subs(m, val) != 0 for val in (4, 1, 3))

# familles-ex4, question 3 : ch, sh, e^{2x} exprimés avec e^{-x}, e^x, e^{2x}.
lm, mu, nu = symbols("lm mu nu")
comb = lm * cosh(x) + mu * sinh(x) + nu * exp(2 * x)
cible = (lm - mu) / 2 * exp(-x) + (lm + mu) / 2 * exp(x) + nu * exp(2 * x)
assert simplify((comb - cible).rewrite(exp)) == 0

# dimension-ex1 : dim 2 ; v1 : dim 2 ; v2 : dim 1.
assert noyau((1, 1, 1, 1), (1, -1, 0, 0)) == 2
assert noyau((0, 1, 1)) == 2  # v1 : P = a + bX + cX^2, b + c = 0
assert noyau((1, -2, 0), (0, 2, -3)) == 1  # v2 : x = 2y, 2y = 3z
assert Matrix([6, 3, 2]).T * Matrix([[1, 0], [-2, 2], [0, -3]]) == Matrix([[0, 0]])
assert all(c == 0 for c in Matrix([[1, 1, 1, 1], [1, -1, 0, 0]]) * Matrix([1, 1, 0, -2])) and \
    all(c == 0 for c in Matrix([[1, 1, 1, 1], [1, -1, 0, 0]]) * Matrix([0, 0, 1, -1]))

# dimension-ex2 : e3 ∈ Vect(u1,u2), e1, e2 non ; v1 : e2 seul dedans ; v2 : aucun.
E3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
def completent(u1, u2):
    return [i + 1 for i, e in enumerate(E3) if libre(u1, u2, e)]


assert libre((1, -1, 2), (2, -2, 1)) and completent((1, -1, 2), (2, -2, 1)) == [1, 2]
assert Rational(2, 3) * Matrix([1, -1, 2]) - Rational(1, 3) * Matrix([2, -2, 1]) == Matrix([0, 0, 1])
assert libre((1, 2, -1), (-1, 1, 1)) and completent((1, 2, -1), (-1, 1, 1)) == [1, 3]  # v1
assert libre((1, 1, 1), (1, 2, 3)) and completent((1, 1, 1), (1, 2, 3)) == [1, 2, 3]  # v2

# dimension-ex3 : dimension du commutant (résolution du système AM = MA).
def dim_commutant(A):
    k = A.shape[0]
    ms = symbols(f"c0:{k * k}")
    M = Matrix(k, k, ms)
    eqs = list(A * M - M * A)
    J = Matrix([[e.coeff(v) for v in ms] for e in eqs])
    return k * k - J.rank()


assert dim_commutant(Matrix.diag(1, 1, 2)) == 5
assert dim_commutant(Matrix.diag(1, 2, 3)) == 3  # v1
assert dim_commutant(Matrix([[0, 1], [0, 0]])) == 2  # v2

# dimension-ex4 : base de Bernstein, coordonnées de 1 (n = 3, n = 4) et de X (n = 3) ; leurres faux.
def coord_bernstein(nn, cible):
    Ps = [X ** k * (1 - X) ** (nn - k) for k in range(nn + 1)]
    M = Matrix([[expand(p).coeff(X, j) for p in Ps] for j in range(nn + 1)])
    assert M.det() != 0  # famille libre, donc base de R_n[X]
    b = Matrix([expand(cible).coeff(X, j) for j in range(nn + 1)])
    return list(M.solve(b))


assert coord_bernstein(3, 1) == [1, 3, 3, 1]
assert coord_bernstein(4, 1) == [1, 4, 6, 4, 1]  # v1
assert coord_bernstein(3, X) == [0, 1, 2, 1]  # v2
assert [binomial(3, k) for k in range(4)] == [1, 3, 3, 1]
for c in ([1, 1, 1, 1], [1, -3, 3, -1], [0, 0, 0, 1]):
    assert c != coord_bernstein(3, 1)
for c in ([1, 1, 1, 1, 1], [1, -4, 6, -4, 1], [0, 0, 0, 0, 1]):
    assert c != coord_bernstein(4, 1)
for c in ([1, 2, 1, 0], [0, 1, 1, 1], [0, 1, -2, 1]):
    assert c != coord_bernstein(3, X)

# somme-directe-ex1 : t(1,a,2) ∈ G ⟺ t(a-1) = 0 ; v1 : 2a - 4 ; v2 : G : z = x + 2y, a + 1.
assert solve(1 + a - 2, a) == [1]
assert solve(2 * a - 1 - 3, a) == [2]  # v1
assert rang((1, 0, 1), (0, 1, 2), (1, -1, -1)) == 2 and rang((1, 0, 1), (0, 1, 2), (1, -1, 0)) == 3  # v2 : a = -1 seul
assert solve(a - (1 - 2), a) == [-1]

# somme-directe-ex2 : décompositions ; leurres faux.
def decomp(dirv, normale, v):
    tt = solve(sum(nc * (vc - t * dc) for nc, vc, dc in zip(normale, v, dirv)), t)[0]
    f = tuple(tt * dc for dc in dirv)
    g = tuple(vc - fc for vc, fc in zip(v, f))
    assert sum(nc * gc for nc, gc in zip(normale, g)) == 0
    assert sum(nc * dc for nc, dc in zip(normale, dirv)) != 0  # somme directe
    return f, g


assert decomp((1, 1, 1), (1, 1, -1), (3, 1, 2)) == ((2, 2, 2), (1, -1, 0))
assert decomp((1, 0, 1), (1, -1, 2), (5, 2, 0)) == ((1, 0, 1), (4, 2, -1))  # v1
assert decomp((1, 2, 0), (1, 1, 1), (4, 5, 0)) == ((3, 6, 0), (1, -1, 0))  # v2
for leurre in [(1, 1, 1), (1, -1, 0), (3, 3, 3)]:
    assert leurre != (2, 2, 2)
for leurre in [(4, 2, -1), (5, 0, 5), (2, 0, 2)]:
    assert leurre != (1, 0, 1)
for leurre in [(1, -1, 0), (1, 2, 0), (4, 8, 0)]:
    assert leurre != (3, 6, 0)

# somme-directe-ex3 : déterminants m - 2, 1 - m, m + 2 ; intersections deux à deux nulles pour tout m.
for (v1_, v2_, v3_), racine in [(((1, 0, 1), (0, 1, 1), (1, 1, m)), 2),
                                (((1, 1, 0), (0, 1, -1), (m, 0, 1)), 1),
                                (((1, 2, 0), (0, 1, 1), (1, 0, m)), -2)]:
    dd = Matrix([v1_, v2_, v3_]).det()
    assert solve(dd, m) == [racine]
    for val in range(-5, 6):
        vs = [Matrix(v).subs(m, val) for v in (v1_, v2_, v3_)]
        assert all(rang(vs[i], vs[j]) == 2 for i in range(3) for j in range(i + 1, 3))
    # leurres « m ≠ 0 », « m ≠ ±1 », « m ≠ 2 » … : la racine est différente de la valeur exclue par le leurre
    assert racine not in ({2: [0], 1: [-1], -2: [2]}[racine])
assert Matrix([1, 0, 1]) + Matrix([0, 1, 1]) - Matrix([1, 1, 2]) == zeros(3, 1)
assert -Matrix([1, 1, 0]) + Matrix([0, 1, -1]) + Matrix([1, 0, 1]) == zeros(3, 1)  # v1
assert -Matrix([1, 2, 0]) + 2 * Matrix([0, 1, 1]) + Matrix([1, 0, -2]) == zeros(3, 1)  # v2

# somme-directe-ex4 : parties impaires ; leurres faux.
def impaire(f):
    return simplify(((f - f.subs(x, -x)) / 2).rewrite(exp))


def egal(f, g):
    return simplify((f - g).rewrite(exp)) == 0


f0 = (x + 1) * exp(x)
assert egal(impaire(f0), x * cosh(x) + sinh(x))
for l in [x * sinh(x) + cosh(x), x * exp(x), sinh(x)]:
    assert not egal(impaire(f0), l)
f1 = x * exp(2 * x)  # v1
assert egal(impaire(f1), x * cosh(2 * x))
for l in [x * sinh(2 * x), sinh(2 * x), x * exp(-2 * x)]:
    assert not egal(impaire(f1), l)
f2 = (2 - x) * exp(x)  # v2
assert egal(impaire(f2), 2 * sinh(x) - x * cosh(x))
for l in [2 * cosh(x) - x * sinh(x), 2 * sinh(x) + x * cosh(x), -x * exp(x)]:
    assert not egal(impaire(f2), l)

# somme-directe-ex5, question 3 : R^4 = F ⊕ G ⊕ H.
assert libre((1, 1, 0, 0), (0, 1, 1, 0), (0, 0, 1, 0), (0, 0, 0, 1))

# grassmann-ex1 : max(0, p + q - n) ; exemples atteignant la borne ; leurres faux.
for nn, p, q, rep, leurres in [(7, 4, 5, 2, [0, 1, 3]), (10, 6, 7, 3, [0, 2, 4]), (6, 2, 3, 0, [-1, 1, 2])]:
    assert max(0, p + q - nn) == rep and rep not in leurres
    I = Matrix.eye(nn)
    Fb = [list(I.row(i)) for i in range(p)]
    Gb = [list(I.row(i)) for i in range(nn - q, nn)]
    assert rang(*Fb) == p and rang(*Gb) == q and dim_inter(Fb, Gb) == rep

# grassmann-ex2 : dim(F ∩ G) = 1 ; v1 : 0 ; v2 : 2.
assert dim_inter([(1, 0, 1, 0), (0, 1, 0, 1)], [(1, 1, 1, 1), (1, 0, 0, 1)]) == 1
assert rang((1, 0, 1, 0), (0, 1, 0, 1), (1, 1, 1, 1), (1, 0, 0, 1)) == 3
assert dim_inter([(1, 1, 0, 0), (0, 0, 1, 1)], [(1, 0, 1, 0), (0, 1, 0, -1)]) == 0  # v1
assert dim_inter([(1, 2, 0, 1), (0, 1, 1, 0)], [(1, 3, 1, 1), (1, 1, -1, 1)]) == 2  # v2
assert Matrix([1, 3, 1, 1]) == Matrix([1, 2, 0, 1]) + Matrix([0, 1, 1, 0])
assert Matrix([1, 1, -1, 1]) == Matrix([1, 2, 0, 1]) - Matrix([0, 1, 1, 0])

# grassmann-ex3 : F = noyau ; dim(F ∩ G), dim(F + G).
def base_noyau(*equations):
    return [list(v) for v in Matrix([list(e) for e in equations]).nullspace()]


F3 = base_noyau((1, 1, -1, 0), (0, 1, 0, 1))
assert len(F3) == 2
assert (dim_inter(F3, [(1, 0, 1, 0), (0, 1, 1, 1)]), rang(*F3, (1, 0, 1, 0), (0, 1, 1, 1))) == (1, 3)
assert (dim_inter(F3, [(1, 0, 0, 0), (0, 1, 0, 0)]), rang(*F3, (1, 0, 0, 0), (0, 1, 0, 0))) == (0, 4)  # v1
F3b = base_noyau((1, -1, 1, 0), (0, 0, 1, 1))
assert (dim_inter(F3b, [(1, 2, 1, -1), (1, 0, -1, 1)]), rang(*F3b, (1, 2, 1, -1), (1, 0, -1, 1))) == (2, 2)  # v2
assert Matrix([1, 2, 1, -1]) == Matrix([1, 1, 0, 0]) + Matrix([0, 1, 1, -1])

# grassmann-ex4 : composantes sur G = R_1[X] (interpolation en deux points) ; leurres faux.
def comp_G(P, p1, p2):
    al, be = symbols("al be")
    sol = solve([al + be * p1 - P.subs(X, p1), al + be * p2 - P.subs(X, p2)], [al, be])
    R = sol[al] + sol[be] * X
    reste = expand(P - R)
    assert reste.subs(X, p1) == 0 and reste.subs(X, p2) == 0
    return expand(R), factor(reste)


R0, F0 = comp_G(X ** 3 + 2 * X ** 2 + 3, 0, 1)
assert R0 == 3 * X + 3 and expand(F0 - X * (X - 1) * (X + 3)) == 0
assert all(expand(R0 - l) != 0 for l in [3, 6 * X + 3, 3 * X])
R1, F1 = comp_G(X ** 3 + X ** 2 + 1, 1, -1)  # v1
assert R1 == X + 2 and expand(F1 - (X ** 2 - 1) * (X + 1)) == 0
assert all(expand(R1 - l) != 0 for l in [1, X + 1, 2 * X + 1])
R2, F2 = comp_G(X ** 3 - X, 0, 2)  # v2
assert R2 == 3 * X and expand(F2 - X * (X - 2) * (X + 2)) == 0
assert all(expand(R2 - l) != 0 for l in [-X, 6 * X, 3 * X + 6])

# inclusion-ex1 : F : x - 2y + z = 0 ; G = Vect(...) ; v1, v2.
Fp = base_noyau((1, -2, 1))
assert len(Fp) == 2
assert all(v[0] - 2 * v[1] + v[2] == 0 for v in [(1, 1, 1), (2, 1, 0), (2, 1, 0), (-1, 0, 1)])
assert rang((1, 1, 1), (2, 1, 0)) == 2 and dim_inter(Fp, [(1, 1, 1), (2, 1, 0)]) == 2  # F = G
assert rang((1, 1, 1), (2, 2, 2)) == 1 and dim_inter(Fp, [(1, 1, 1), (2, 2, 2)]) == 1  # v1 : G droite ⊂ F
assert rang((1, 1, 1), (1, 0, 0)) == 2 and dim_inter(Fp, [(1, 1, 1), (1, 0, 0)]) == 1  # v2 : F ∩ G droite
assert 1 - 0 + 0 != 0  # v2 : (1,0,0) ∉ F

# inclusion-ex2 : dimensions des espaces de polynômes (conditions linéaires sur les coefficients).
def dim_pol(deg, conditions):
    cs = symbols(f"k0:{deg + 1}")
    Pk = sum(cs[i] * X ** i for i in range(deg + 1))
    eqs = [cond(Pk) for cond in conditions]
    J = Matrix([[expand(e).coeff(v) for v in cs] for e in eqs])
    return deg + 1 - J.rank()


assert dim_pol(3, [lambda P: P.subs(X, 1)]) == 3
assert dim_pol(3, [lambda P: P.subs(X, 1), lambda P: P.subs(X, -1)]) == 2  # v1
assert dim_pol(4, [lambda P: P.subs(X, 0), lambda P: diff(P, X).subs(X, 0)]) == 3  # v2
assert dim_pol(4, [lambda P: P.subs(X, 0)]) == 4  # v2 : dim H = 4
assert (X - 1).subs(X, -1) == -2  # v1 : X - 1 ∈ H \ F

# inclusion-ex3 : dim H = 3, dim(F ∩ H).
H = base_noyau((1, 1, 1, 1))
assert len(H) == 3
assert dim_inter([(1, 0, 0, 0), (0, 1, -1, 0), (0, 0, 1, -1)], H) == 2
assert rang((1, 0, 0, 0), (0, 1, -1, 0), (0, 0, 1, -1)) == 3
assert dim_inter([(1, 1, 0, 0), (0, 0, 1, 1)], H) == 1  # v1
assert dim_inter([(1, -1, 0, 0), (0, 1, -1, 0), (0, 0, 1, -1)], H) == 3  # v2

# inclusion-ex4 : racines, coefficients, suites vérifiant la relation ; leurres faux (en n = 0 ou n = 1).
r = symbols("r")
for (A1, A0), (u0, u1), racines, sol, leurres in [
    ((1, 2), (1, 5), {2, -1}, 2 ** (n + 1) - (-1) ** n,
     [2 ** (n + 1) + (-1) ** n, 3 * 2 ** n - 2 * (-1) ** n, 2 ** (n + 2) - 3 * (-1) ** n]),
    ((5, -6), (0, 1), {2, 3}, 3 ** n - 2 ** n,
     [2 ** n - 3 ** n, 3 ** n + 2 ** n, 2 * 3 ** n - 2 ** (n + 1)]),
    ((1, 6), (2, 1), {3, -2}, 3 ** n + (-2) ** n,
     [3 ** n - (-2) ** n, 3 ** n + 2 ** n, (-3) ** n + 2 ** n]),
]:
    assert set(solve(r ** 2 - A1 * r - A0, r)) == racines
    assert 0 - A1 * 0 - A0 * 1 != 0  # la suite (0^n) = (1, 0, 0, ...) n'est pas solution
    vals = [sol.subs(n, k) for k in range(8)]
    assert vals[0] == u0 and vals[1] == u1
    assert all(vals[k + 2] == A1 * vals[k + 1] + A0 * vals[k] for k in range(6))
    for l in leurres:
        assert (l.subs(n, 0), l.subs(n, 1)) != (u0, u1)

print("al_l1_espaces : toutes les vérifications passent.")
