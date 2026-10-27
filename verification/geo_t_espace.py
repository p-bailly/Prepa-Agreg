"""Vérifications calculatoires du chapitre « geo-t-espace » (géométrie dans l'espace, Terminale).

Fichiers vérifiés : content/cartes/geometrie/terminale-geo-t-espace.json et
content/exercices/geometrie/terminale-geo-t-espace.json. Chaque assert cite l'identifiant de la carte ou de
l'exercice (et la variante : v1, v2). Les leurres des propositions sont vérifiés FAUX.
Lancer : python3 verification/geo_t_espace.py (doit finir sans erreur).
"""

from sympy import Matrix, Rational as Q, nsimplify, simplify, solve, sqrt, symbols

a, b, c, t, s, m, k = symbols("a b c t s m k")


def V(*x):
    return Matrix([nsimplify(v) for v in x])


def colin(u, v):
    return Matrix.hstack(u, v).rank() <= 1


def coplan(u, v, w):
    return Matrix.hstack(u, v, w).det() == 0


def decomp(u, v, w):
    """Couples (a, b) tels que w = a u + b v (liste vide si aucun)."""
    return solve(list(a * u + b * v - w), [a, b], dict=True)


def sur_droite(P, u, M):
    return colin(M - P, u)


def meme_droite(P, u, Q_, v):
    return colin(u, v) and sur_droite(P, u, Q_)


def position(P, u, Q_, v):
    """Position relative de deux droites (P, u) et (Q, v)."""
    if colin(u, v):
        return "confondues" if sur_droite(P, u, Q_) else "paralleles"
    sol = solve(list(P + t * u - Q_ - s * v), [t, s], dict=True)
    if sol:
        return ("secantes", P + sol[0][t] * u)
    return "non coplanaires"


def plan(n, A):
    """Équation (a, b, c, d) du plan passant par A de vecteur normal n (repère orthonormé)."""
    return (n[0], n[1], n[2], -(n.dot(A)))


def val(p, M):
    return p[0] * M[0] + p[1] * M[1] + p[2] * M[2] + p[3]


def meme_plan(p, q):
    return Matrix([list(p), list(q)]).rank() == 1


def normale(p):
    return V(p[0], p[1], p[2])


def dist(M, p):
    return simplify(abs(val(p, M)) / normale(p).norm())


def projete(M, p):
    n = normale(p)
    t0 = -val(p, M) / n.dot(n)
    return M + t0 * n


def plan3(A, B, C):
    n = (B - A).cross(C - A)
    assert n != V(0, 0, 0)
    return plan(n, A)


# =============================================================================
# geo-t-espace-coplanaires (carte)
# =============================================================================
u, v, w = V(1, 0, 2), V(0, 1, -1), V(2, 3, 1)
assert not colin(u, v) and decomp(u, v, w) == [{a: 2, b: 3}]  # exemple : w = 2u + 3v
# contre-exemple 1 : u(1,0,0), v(2,0,0) colinéaires, w(0,1,0) : coplanaires mais pas de décomposition
u, v, w = V(1, 0, 0), V(2, 0, 0), V(0, 1, 0)
assert colin(u, v) and coplan(u, v, w) and decomp(u, v, w) == []
assert 2 * u - v + 0 * w == V(0, 0, 0)  # affirmation 2 : relation non triviale
# contre-exemple 2 / affirmation 1 : deux à deux non colinéaires mais coplanaires
u, v, w = V(1, 0, 0), V(0, 1, 0), V(1, 1, 0)
assert not colin(u, v) and not colin(u, w) and not colin(v, w) and coplan(u, v, w) and w == u + v
assert solve(list(a * u + b * v + c * w - V(0, 0, 1)), [a, b, c], dict=True) == []
# trous : leurres faux. « non nuls » / « de même norme » : u = v non nul, w hors de la droite
assert decomp(V(1, 0, 0), V(1, 0, 0), V(0, 1, 0)) == []
# « couple d'entiers » : w = u/2 ; « réels positifs » : w = -u ; « a + b = 1 » : w = 2u (u, v non colinéaires)
u, v = V(1, 0, 0), V(0, 1, 0)
assert decomp(u, v, u / 2) == [{a: Q(1, 2), b: 0}]
assert decomp(u, v, -u) == [{a: -1, b: 0}]
assert decomp(u, v, 2 * u) == [{a: 2, b: 0}] and 2 + 0 != 1

# geo-t-espace-coplanaires-ex1 : w = 2u + v ; leurres faux
u, v, w = V(1, 2, -1), V(0, 1, 3), V(2, 5, 1)
assert not colin(u, v) and decomp(u, v, w) == [{a: 2, b: 1}]
assert u + 2 * v != w and 2 * u - v != w and coplan(u, v, w)
# ex1 v1 : w = 3u - v
u, v, w = V(1, -1, 2), V(1, 0, 1), V(2, -3, 5)
assert not colin(u, v) and decomp(u, v, w) == [{a: 3, b: -1}]
assert -u + 3 * v != w and 3 * u + v != w
# ex1 v2 : non coplanaires (det = -1)
u, v, w = V(2, 1, 0), V(1, 0, -1), V(1, 2, 4)
assert not colin(u, v) and decomp(u, v, w) == [] and Matrix.hstack(u, v, w).det() == -1
assert 2 * u - 3 * v != w and 2 * u - 4 * v != w and u + v != w

# geo-t-espace-coplanaires-ex2 : m = 3 (unique), leurres 2, 1 faux, « aucune » faux
def m_coplan(u, v, w):
    return solve(Matrix.hstack(u, v, w).det(), m)


assert m_coplan(V(1, 1, 0), V(0, 1, 1), V(1, m, 2)) == [3]
assert decomp(V(1, 1, 0), V(0, 1, 1), V(1, 3, 2)) == [{a: 1, b: 2}]
# ex2 v1 : m = 4 (w = 3u + 2v) ; leurres 8, 2 faux
assert m_coplan(V(1, 0, 2), V(0, 1, -1), V(3, 2, m)) == [4]
assert decomp(V(1, 0, 2), V(0, 1, -1), V(3, 2, 4)) == [{a: 3, b: 2}]
# ex2 v2 : m = 8 (w = 3u + 2v) ; leurres 4, 7 faux
assert m_coplan(V(2, 1, 1), V(1, -1, 0), V(m, 1, 3)) == [8]
assert decomp(V(2, 1, 1), V(1, -1, 0), V(8, 1, 3)) == [{a: 3, b: 2}]

# geo-t-espace-coplanaires-ex3 : 2 IJ = AC + BD pour un tétraèdre quelconque (coordonnées symboliques)
xs = symbols("xa ya za xb yb zb xc yc zc xd yd zd")
A_, B_, C_, D_ = (Matrix(xs[i:i + 3]) for i in range(0, 12, 3))
I_, J_ = (A_ + B_) / 2, (C_ + D_) / 2
assert simplify(2 * (J_ - I_) - (C_ - A_) - (D_ - B_)) == Matrix([0, 0, 0])

# geo-t-espace-coplanaires-ex4 : repère (A; AB, AC, AD), A = 0, B = e1, C = e2, D = e3
A0, B0, C0, D0 = V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)


def t_coplan(I, M, K):
    L = C0 + t * (D0 - C0)
    assert not colin(M - I, K - I)
    return solve(Matrix.hstack(M - I, K - I, L - I).det(), t), L


# original : I milieu [AB], M milieu [BC], AK = AD/3 -> t = 1/3, a = 4/3, b = 1, KL = (2/3) AC
I, M, K = (A0 + B0) / 2, (B0 + C0) / 2, D0 / 3
sol, L = t_coplan(I, M, K)
assert sol == [Q(1, 3)]
L = L.subs(t, Q(1, 3))
assert decomp(M - I, K - I, L - I) == [{a: Q(4, 3), b: 1}]
assert L - K == Q(2, 3) * (C0 - A0) and M - I == Q(1, 2) * (C0 - A0)
# v1 : AI = AB/3, BM = (2/3) BC, AK = (2/5) AD -> t = 2/5, a = 9/10, KL = (3/5) AC
I, M, K = B0 / 3, B0 + Q(2, 3) * (C0 - B0), Q(2, 5) * D0
assert M == V(Q(1, 3), Q(2, 3), 0)
sol, L = t_coplan(I, M, K)
assert sol == [Q(2, 5)]
L = L.subs(t, Q(2, 5))
assert decomp(M - I, K - I, L - I) == [{a: Q(9, 10), b: 1}]
assert L - K == Q(3, 5) * C0 and M - I == Q(2, 3) * C0
# v2 : I milieu [AB], M milieu [BD], AK = AC/4 -> t = 3/4, a = 3/2, KL = (3/4) AD
I, M, K = (A0 + B0) / 2, (B0 + D0) / 2, C0 / 4
sol, L = t_coplan(I, M, K)
assert sol == [Q(3, 4)]
L = L.subs(t, Q(3, 4))
assert decomp(M - I, K - I, L - I) == [{a: Q(3, 2), b: 1}]
assert L - K == Q(3, 4) * D0 and M - I == Q(1, 2) * D0
# leurres 2/3, 1/2 (original), 3/5, 1/3 (v1), 1/4, 1/2 (v2) : différents de la valeur unique trouvée
assert {Q(2, 3), Q(1, 2)}.isdisjoint({Q(1, 3)}) and {Q(3, 5), Q(1, 3)}.isdisjoint({Q(2, 5)})
assert {Q(1, 4), Q(1, 2)}.isdisjoint({Q(3, 4)})

# =============================================================================
# geo-t-espace-droite-param (carte)
# =============================================================================
A, B = V(1, -2, 3), V(3, 0, -1)
assert B - A == V(2, 2, -4)
assert A + Q(-1, 2) * (B - A) == V(0, -3, 5)  # exemple : C de paramètre t = -1/2
assert B + (-3) * V(1, 1, -2) == V(0, -3, 5)  # exemple : s = -3 dans la seconde représentation
assert meme_droite(A, B - A, B, V(1, 1, -2))
# contre-exemple / affirmation 1 : même vecteur directeur, droites distinctes
assert not meme_droite(V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(1, 0, 0))
# affirmation 2 : x = t, y = t, z = 0 et x = 1 + 2s, y = 1 + 2s, z = 0 : même droite
assert meme_droite(V(0, 0, 0), V(1, 1, 0), V(1, 1, 0), V(2, 2, 0))
# affirmation 3 : (Ox) et la droite x = 0, y = s, z = 1 : non coplanaires
assert position(V(0, 0, 0), V(1, 0, 0), V(0, 0, 1), V(0, 1, 0)) == "non coplanaires"

# geo-t-espace-droite-param-ex1 : C de paramètre t = -2 ; leurres 2 et -1 faux
A, B, C = V(2, -1, 0), V(1, 1, 3), V(4, -5, -6)
assert B - A == V(-1, 2, 3) and A + (-2) * (B - A) == C
assert A + 2 * (B - A) != C and A + (-1) * (B - A) != C
# ex1 v1 : t = 3 ; leurres -3, 2 faux
A, B, C = V(0, 3, -1), V(2, 2, 1), V(6, 0, 5)
assert B - A == V(2, -1, 2) and A + 3 * (B - A) == C
assert A + (-3) * (B - A) != C and A + 2 * (B - A) != C
# ex1 v2 : C n'est pas sur (AB) (t = 2 convient pour x et y mais z = 0 != 1)
A, B, C = V(1, 0, 2), V(2, 2, 1), V(3, 4, 1)
assert not sur_droite(A, B - A, C) and A + 2 * (B - A) == V(3, 4, 0)

# geo-t-espace-droite-param-ex2 : seule (c) représente D
D0p, D0u = V(3, 1, 0), V(-2, 4, -1)
syst = {"a": (V(1, 5, 1), V(2, -4, 1)), "b": (V(3, 1, 0), V(2, 4, -1)), "c": (V(1, 5, -1), V(2, -4, 1)),
        "d": (V(-2, 4, -1), V(3, 1, 0))}
assert [x for x, (P, w_) in syst.items() if meme_droite(D0p, D0u, P, w_)] == ["c"]
assert D0p + 1 * D0u == V(1, 5, -1)
# ex2 v1 : seule (a)
D0p, D0u = V(1, -2, 4), V(1, 3, -2)
syst = {"a": (V(3, 4, 0), V(-1, -3, 2)), "b": (V(3, 4, 1), V(-1, -3, 2)), "c": (V(1, -2, 4), V(1, 3, 2)),
        "d": (V(1, 3, -2), V(1, -2, 4))}
assert [x for x, (P, w_) in syst.items() if meme_droite(D0p, D0u, P, w_)] == ["a"]
assert D0p + 2 * D0u == V(3, 4, 0)
# ex2 v2 : seule (d)
D0p, D0u = V(0, 2, 5), V(-1, 1, 3)
syst = {"a": (V(2, 0, 1), V(2, -2, -6)), "b": (V(0, 2, 5), V(2, 2, 6)), "c": (V(-1, 1, 3), V(0, 2, 5)),
        "d": (V(2, 0, -1), V(2, -2, -6))}
assert [x for x, (P, w_) in syst.items() if meme_droite(D0p, D0u, P, w_)] == ["d"]
assert D0p + (-2) * D0u == V(2, 0, -1)

# geo-t-espace-droite-param-ex3 : sécantes en (2,1,3) (t = 1, s = 2) ; leurre (3,0,6) = point t = 2 de D1
P1, u1, P2, u2 = V(1, 2, 0), V(1, -1, 3), V(0, -3, 5), V(1, 2, -1)
assert position(P1, u1, P2, u2) == ("secantes", V(2, 1, 3))
assert solve(list(P1 + t * u1 - P2 - s * u2), [t, s], dict=True) == [{t: 1, s: 2}]
assert P1 + 2 * u1 == V(3, 0, 6) and not sur_droite(P2, u2, V(3, 0, 6))
# ex3 v1 : non coplanaires ; les deux équations donnent s = -2, t = -3, points (-3,-2,3) et (-3,-2,0)
P1, u1, P2, u2 = V(0, 1, 0), V(1, 1, -1), V(1, 0, 2), V(2, 1, 1)
assert position(P1, u1, P2, u2) == "non coplanaires"
assert P1 + (-3) * u1 == V(-3, -2, 3) and P2 + (-2) * u2 == V(-3, -2, 0)
# ex3 v2 : parallèles distinctes ; leurre (1,3,3) = point t = 1 de D1, pas sur D2
P1, u1, P2, u2 = V(2, 1, 3), V(-1, 2, 0), V(1, 0, 3), V(2, -4, 0)
assert position(P1, u1, P2, u2) == "paralleles" and P1 + u1 == V(1, 3, 3)
assert not sur_droite(P2, u2, V(1, 3, 3))


# geo-t-espace-droite-param-ex4 : sécantes ssi m = 1, en (0,1,1) ; sinon non coplanaires
def valeurs_m(P1, u1, P2m, u2):
    assert not colin(u1, u2)
    sols = solve(list(P1 + t * u1 - P2m - s * u2), [t, s, m], dict=True)
    return sols


assert valeurs_m(V(1, 0, 2), V(1, -1, 1), V(0, 1, m), V(3, 1, -1)) == [{t: -1, s: 0, m: 1}]
assert V(1, 0, 2) + (-1) * V(1, -1, 1) == V(0, 1, 1)
assert position(V(1, 0, 2), V(1, -1, 1), V(0, 1, -1), V(3, 1, -1)) == "non coplanaires"  # leurre m = -1
# ex4 v1 : m = 0, point (2,2,0) ; leurre m = 2 : non coplanaires
assert valeurs_m(V(0, 1, -1), V(2, 1, 1), V(1, m, 1), V(1, 2, -1)) == [{t: 1, s: 1, m: 0}]
assert V(0, 1, -1) + V(2, 1, 1) == V(2, 2, 0)
assert position(V(0, 1, -1), V(2, 1, 1), V(1, 2, 1), V(1, 2, -1)) == "non coplanaires"
# ex4 v2 : m = 3, point (-1,2,0) ; leurre m = 1 : non coplanaires
assert valeurs_m(V(-1, 2, 0), V(1, 0, 3), V(0, m, 2), V(1, 1, 2)) == [{t: 0, s: -1, m: 3}]
assert position(V(-1, 2, 0), V(1, 0, 3), V(0, 1, 2), V(1, 1, 2)) == "non coplanaires"

# =============================================================================
# geo-t-espace-equation-plan (carte)
# =============================================================================
assert plan(V(2, -1, 3), V(1, 2, -1)) == (2, -1, 3, 3)  # exemple
assert val((1, 0, -4, 7), V(-7, 0, 0)) == 0
# contre-exemple « repère non orthonormé » : e1 = i, e2 = i + j, e3 = k ; plan x = 0 dirigé par e2, e3
e1, e2, e3 = V(1, 0, 0), V(1, 1, 0), V(0, 0, 1)
assert e1.dot(e2) == 1  # e1 n'est pas normal au plan
# leurre de trou « orthogonal non normé » : base 2i, j, k ; plan x' + y' = 0 ; n de coordonnées (1,1,0)
f1, f2, f3 = V(2, 0, 0), V(0, 1, 0), V(0, 0, 1)
direction = [f1 - f2, f3]  # points x' + y' = 0 : vecteurs (x', -x', z') -> x'(f1 - f2) + z' f3
n_coord = f1 + f2
assert any(n_coord.dot(d) != 0 for d in direction)
# leurre de trou « normé non orthogonal » : base i, (i + j)/sqrt2, k ; plan x' = 0 ; n de coordonnées (1,0,0)
g2 = V(1, 1, 0) / sqrt(2)
assert simplify(e1.dot(g2)) != 0
# contre-exemple « un seul vecteur » / affirmation 4 : plan z = 0, n = (1,0,0)
assert V(1, 0, 0).dot(V(0, 1, 0)) == 0 and V(1, 0, 0).dot(V(1, 0, 0)) == 1
# affirmation 1 : x + y - z + 1 = 0 et 2x + 2y - 2z + 1 = 0 distincts ; (0,0,1)
assert not meme_plan((1, 1, -1, 1), (2, 2, -2, 1))
assert val((1, 1, -1, 1), V(0, 0, 1)) == 0 and val((2, 2, -2, 1), V(0, 0, 1)) == -1
# affirmation 2 : x + 2y = 3 contient trois points non alignés
pts = [V(3, 0, 0), V(1, 1, 0), V(3, 0, 5)]
assert all(val((1, 2, 0, -3), p) == 0 for p in pts) and not colin(pts[1] - pts[0], pts[2] - pts[0])
# démonstration : (-b, a, 0) et (-c, 0, a) orthogonaux à (a, b, c), non colinéaires si a != 0
n_ = Matrix([a, b, c])
assert Matrix([-b, a, 0]).dot(n_) == 0 and Matrix([-c, 0, a]).dot(n_) == 0
assert Matrix.hstack(Matrix([-b, a, 0]), Matrix([-c, 0, a])).extract([1, 2], [0, 1]).det() == a**2

# geo-t-espace-equation-plan-ex1 : x + 3y - 2z + 9 = 0 ; leurres faux (ne passent pas par A, ou mauvaise normale)
A = V(2, -1, 4)
p = plan(V(1, 3, -2), A)
assert p == (1, 3, -2, 9)
for q in [(1, 3, -2, -9), (1, 3, -2, 0), (2, -1, 4, 9)]:
    assert not meme_plan(p, q)
# ex1 v1 : 4x - y + 2z + 7 = 0
p = plan(V(4, -1, 2), V(0, 5, -1))
assert p == (4, -1, 2, 7)
for q in [(4, -1, 2, -7), (4, -1, 2, 0), (0, 5, -1, 7)]:
    assert not meme_plan(p, q)
# ex1 v2 : 2x - 5z + 11 = 0
p = plan(V(2, 0, -5), V(-3, 2, 1))
assert p == (2, 0, -5, 11)
for q in [(2, 0, -5, -11), (2, 0, -5, 0), (2, -5, 0, 11)]:
    assert not meme_plan(p, q)

# geo-t-espace-equation-plan-ex2 : plan parallèle passant par A
for P_, A, attendu, leurres in [
    ((2, -1, 1, -7), V(1, 1, 1), (2, -1, 1, -2), [(2, -1, 1, -7), (2, -1, 1, 2), (1, 1, 1, -3)]),  # original
    ((1, 4, -2, 1), V(0, -2, 3), (1, 4, -2, 14), [(1, 4, -2, 1), (1, 4, -2, -14), (0, -2, 3, -13)]),  # v1
    ((3, 0, 2, -4), V(-1, 2, 5), (3, 0, 2, -7), [(3, 0, 2, -4), (3, 0, 2, 7), (3, 2, 0, -1)]),  # v2
]:
    p = plan(normale(P_), A)
    assert p == attendu and val(P_, A) != 0
    for q in leurres:
        assert not meme_plan(p, q)

# geo-t-espace-equation-plan-ex3 : plan par trois points (n vérifié orthogonal à AB et AC)
for A, B, C, n_att, attendu, leurres in [
    (V(1, 2, 0), V(2, 0, 1), V(0, 1, 3), V(5, 4, 3), (5, 4, 3, -13),
     [(5, 4, 3, 13), (1, -2, 1, 3), (5, -4, 3, 3)]),  # original
    (V(2, 0, 1), V(0, 1, 1), V(1, 1, 0), V(1, 2, 1), (1, 2, 1, -3),
     [(1, 2, 1, 3), (-2, 1, 0, 4), (1, 2, -1, -1)]),  # v1
    (V(1, -1, 2), V(3, 0, 1), V(0, 1, 4), V(4, -3, 5), (4, -3, 5, -17),
     [(4, -3, 5, 17), (2, 1, -1, 1), (4, 3, 5, -11)]),  # v2
]:
    assert not colin(B - A, C - A)
    assert n_att.dot(B - A) == 0 and n_att.dot(C - A) == 0
    p = plan(n_att, A)
    assert p == attendu and val(p, B) == 0 and val(p, C) == 0
    assert meme_plan(p, plan3(A, B, C))
    for q in leurres:
        assert not meme_plan(p, q)
# leurres plausibles du v1 : x + 2y - z - 1 = 0 passe par A et B mais pas par C
assert val((1, 2, -1, -1), V(2, 0, 1)) == 0 and val((1, 2, -1, -1), V(0, 1, 1)) == 0
assert val((1, 2, -1, -1), V(1, 1, 0)) != 0


# geo-t-espace-equation-plan-ex4 : famille de plans ; D incluse ssi m = m0, sinon sécante en A
def etude_famille(Pt, u, coeffs):
    """coeffs(m) -> (a, b, c, d) ; renvoie m0 tel que D soit incluse, après contrôle de A dans tous les P_m."""
    p = coeffs(m)
    assert simplify(val(p, Pt)) == 0
    expr = simplify(val(p, Pt + t * u))
    coef = simplify(expr.coeff(t, 1))
    assert simplify(expr - coef * t) == 0  # pas de terme constant : t = 0 (point A) toujours solution
    return solve(coef, m)


assert etude_famille(V(1, -1, 3), V(2, 1, -1), lambda m: (1, m, 1, m - 4)) == [-1]
assert normale((1, -1, 1, 0)).dot(V(2, 1, -1)) == 0
# Delta : x = k, y = -1, z = 4 - k dans tous les P_m
assert simplify(val((1, m, 1, m - 4), V(k, -1, 4 - k))) == 0
# leurres m = 1, m = 4 : D non incluse (coefficient 1 + m non nul)
assert all(1 + mm != 0 for mm in [1, 4])
# v1 : m = 3 ; Delta : x = k, y = 1, z = 2 - k ; leurres -3, 2
assert etude_famille(V(2, 1, 0), V(1, -1, 2), lambda m: (1, m, 1, -m - 2)) == [3]
assert simplify(val((1, m, 1, -m - 2), V(k, 1, 2 - k))) == 0
assert all(3 - mm != 0 for mm in [-3, 2])
# v2 : m = -4 ; Delta : x = k, y = 2 - k, z = -1 ; leurres 4, 2
assert etude_famille(V(0, 2, -1), V(1, 3, 1), lambda m: (1, 1, m, m - 2)) == [-4]
assert simplify(val((1, 1, m, m - 2), V(k, 2 - k, -1))) == 0
assert all(4 + mm != 0 for mm in [4, 2])

# =============================================================================
# geo-t-espace-orthogonalite (carte)
# =============================================================================
assert V(4, -2, 4) == 2 * V(2, -1, 2)  # exemple
# cube : AE orthogonal à AB et AD, et à BD
assert V(0, 0, 1).dot(V(1, 0, 0)) == 0 and V(0, 0, 1).dot(V(0, 1, 0)) == 0 and V(0, 0, 1).dot(V(-1, 1, 0)) == 0
# contre-exemples : u(0,1,1) orthogonal à (1,0,0) et (2,0,0) mais pas à (0,1,0)
uu = V(0, 1, 1)
assert uu.dot(V(1, 0, 0)) == 0 and uu.dot(V(2, 0, 0)) == 0 and uu.dot(V(0, 1, 0)) == 1
assert not colin(uu, V(0, 0, 1))
# (Oz) et la droite y = 1, z = 0 : orthogonales, non sécantes
assert V(0, 0, 1).dot(V(1, 0, 0)) == 0 and position(V(0, 0, 0), V(0, 0, 1), V(0, 1, 0), V(1, 0, 0)) == "non coplanaires"
# affirmation 2 : (Ox), (Oy) orthogonales à (Oz) et sécantes
assert position(V(0, 0, 0), V(1, 0, 0), V(0, 0, 0), V(0, 1, 0)) == ("secantes", V(0, 0, 0))


# geo-t-espace-orthogonalite-ex1 : réponses (0 : colinéaire ; 2 : parallèle ; 3 : sécante non orthogonale)
def nature(u, n, Pt, p):
    if colin(u, n):
        return 0
    if u.dot(n) == 0:
        assert val(p, Pt) != 0  # strictement parallèle
        return 2
    return 3


assert nature(V(-4, 2, -8), V(2, -1, 4), V(1, 0, 3), (2, -1, 4, 1)) == 0 and V(-4, 2, -8) == -2 * V(2, -1, 4)
assert V(-4, 2, -8).dot(V(2, -1, 4)) == -42
assert nature(V(1, 1, 1), V(1, 3, -1), V(2, 0, -1), (1, 3, -1, -2)) == 3 and V(1, 1, 1).dot(V(1, 3, -1)) == 3  # v1
assert nature(V(2, 1, 0), V(1, -2, 1), V(1, 3, 2), (1, -2, 1, 4)) == 2 and val((1, -2, 1, 4), V(1, 3, 2)) == 1  # v2

# geo-t-espace-orthogonalite-ex2 : plan passant par A de normale u
for A, u, attendu, leurres in [
    (V(2, 1, -1), V(-1, 2, 3), (-1, 2, 3, 3), [(-1, 2, 3, -3), (3, 1, 5, -2), (1, 2, 3, -1)]),
    (V(0, 3, 1), V(2, -1, 1), (2, -1, 1, 2), [(2, -1, 1, -2), (1, 0, 4, -4), (2, 1, 1, -4)]),
    (V(-1, 2, 2), V(5, -1, 2), (5, -1, 2, 3), [(5, -1, 2, -3), (0, 2, -3, 2), (5, 1, 2, -1)]),
]:
    p = plan(u, A)
    assert p == attendu
    for q in leurres:
        assert not meme_plan(p, q)
    # trois des leurres passent par A (erreur de vecteur normal), pas le premier
    assert val(leurres[0], A) != 0 and all(val(q, A) == 0 for q in leurres[1:])

# geo-t-espace-orthogonalite-ex3 : plan médiateur = {MA^2 = MB^2}
x, y, z = symbols("x y z")
Mxyz = Matrix([x, y, z])
for A, B, attendu, leurres in [
    (V(1, 2, -1), V(3, 0, 3), (1, -1, 2, -3), [(1, -1, 2, 3), (1, -1, 2, -9), (2, 1, 1, -6)]),
    (V(0, 1, 4), V(2, -3, 0), (1, -2, -2, 1), [(1, -2, -2, 10), (1, -2, -2, -8), (1, -1, 2, -6)]),
    (V(-2, 0, 1), V(4, 2, -1), (3, 1, -1, -4), [(3, 1, -1, 7), (3, 1, -1, -15), (1, 1, 0, -2)]),
]:
    diff = ((Mxyz - A).dot(Mxyz - A) - (Mxyz - B).dot(Mxyz - B)).expand()
    eq = val(attendu, Mxyz)
    assert simplify(diff / eq).is_number  # même équation à un facteur près
    I = (A + B) / 2
    assert val(attendu, I) == 0 and colin(B - A, normale(attendu))
    for q in leurres:
        assert not meme_plan(attendu, q)
# développements écrits dans les corrections (facteur 4 : 4x - 4y + 8z - 12 ; v1 : 4x - 8y - 8z + 4 ; v2 : 12x + 4y - 4z - 16)
assert (((Mxyz - V(1, 2, -1)).dot(Mxyz - V(1, 2, -1)) - (Mxyz - V(3, 0, 3)).dot(Mxyz - V(3, 0, 3))).expand()
        == 4 * x - 4 * y + 8 * z - 12)
assert (((Mxyz - V(0, 1, 4)).dot(Mxyz - V(0, 1, 4)) - (Mxyz - V(2, -3, 0)).dot(Mxyz - V(2, -3, 0))).expand()
        == 4 * x - 8 * y - 8 * z + 4)
assert (((Mxyz - V(-2, 0, 1)).dot(Mxyz - V(-2, 0, 1)) - (Mxyz - V(4, 2, -1)).dot(Mxyz - V(4, 2, -1))).expand()
        == 12 * x + 4 * y - 4 * z - 16)

# geo-t-espace-orthogonalite-ex4 : cube d'arête 1, A origine
cube = dict(A=V(0, 0, 0), B=V(1, 0, 0), C=V(1, 1, 0), D=V(0, 1, 0), E=V(0, 0, 1), F=V(1, 0, 1), G=V(1, 1, 1),
            H=V(0, 1, 1))
for (o, e), (p1, p2, p3), K_att, leurres in [
    (("A", "G"), ("B", "D", "E"), V(Q(1, 3), Q(1, 3), Q(1, 3)),
     [V(Q(1, 2), Q(1, 2), Q(1, 2)), V(Q(2, 3), Q(2, 3), Q(2, 3)), V(Q(1, 3), Q(1, 3), 0)]),
    (("A", "G"), ("C", "F", "H"), V(Q(2, 3), Q(2, 3), Q(2, 3)),
     [V(Q(1, 3), Q(1, 3), Q(1, 3)), V(Q(1, 2), Q(1, 2), Q(1, 2)), V(Q(2, 3), Q(2, 3), 0)]),
    (("B", "H"), ("A", "C", "F"), V(Q(2, 3), Q(1, 3), Q(1, 3)),
     [V(Q(1, 3), Q(2, 3), Q(2, 3)), V(Q(1, 2), Q(1, 2), Q(1, 2)), V(Q(1, 3), Q(1, 3), Q(1, 3))]),
]:
    O_, d = cube[o], cube[e] - cube[o]
    P1, P2, P3 = cube[p1], cube[p2], cube[p3]
    assert d.dot(P2 - P1) == 0 and d.dot(P3 - P1) == 0 and not colin(P2 - P1, P3 - P1)
    p = plan(d, P1)
    assert val(p, P2) == 0 and val(p, P3) == 0
    sol = solve(val(p, O_ + t * d), t)
    assert len(sol) == 1 and O_ + sol[0] * d == K_att
    assert all((K_att - P).dot(K_att - P) == Q(2, 3) for P in (P1, P2, P3))  # KB = KD = KE = sqrt(6)/3
    for L_ in leurres:
        assert val(p, L_) != 0 or not sur_droite(O_, d, L_)
assert plan(V(1, 1, 1), cube["B"]) == (1, 1, 1, -1) and plan(V(1, 1, 1), cube["C"]) == (1, 1, 1, -2)
assert plan(V(-1, 1, 1), cube["A"]) == (-1, 1, 1, 0)
assert dist(cube["A"], (1, 1, 1, -1)) == sqrt(3) / 3  # AK = distance de A à (BDE)

# =============================================================================
# geo-t-espace-projete (carte)
# =============================================================================
p = (1, 2, -2, 3)
assert projete(V(3, 3, -3), p) == V(1, -1, 1) and dist(V(3, 3, -3), p) == 6  # exemple
assert (V(3, 3, -3) - V(1, -1, 1)).norm() == 6
# démonstration : MN^2 = MH^2 + HN^2 pour N dans P (paramétré) ; minimalité
M_ = V(3, 3, -3)
H_ = projete(M_, p)
u1, u2 = V(2, -1, 0), V(2, 0, 1)
assert normale(p).dot(u1) == 0 and normale(p).dot(u2) == 0
N_ = H_ + a * u1 + b * u2
assert simplify((N_ - M_).dot(N_ - M_) - (H_ - M_).dot(H_ - M_) - (N_ - H_).dot(N_ - H_)) == 0
# contre-exemple repère non orthonormé : vraie distance 1/sqrt2, formule naïve 1
assert dist(V(1, 0, 0), (1, -1, 0, 0)) == 1 / sqrt(2) and 1 / sqrt(2) != 1
# contre-exemple / affirmation 2 : P : z = 0, M(0,0,1), N(0,5,0)
assert (V(0, 5, 0) - V(0, 0, 1)).dot(V(1, 0, 0)) == 0 and projete(V(0, 0, 1), (0, 0, 1, 0)) == V(0, 0, 0)
assert (V(0, 5, 0) - V(0, 0, 1)).norm() == sqrt(26)
# affirmation 1 : d(O, 2x + 2y + z - 6 = 0) = 2 (et non 6), H(4/3, 4/3, 2/3)
assert dist(V(0, 0, 0), (2, 2, 1, -6)) == 2 and projete(V(0, 0, 0), (2, 2, 1, -6)) == V(Q(4, 3), Q(4, 3), Q(2, 3))
# affirmation 4 : (0,0,1) et (1,0,1) à distance 1 de z = 0, projetés distincts
assert dist(V(0, 0, 1), (0, 0, 1, 0)) == dist(V(1, 0, 1), (0, 0, 1, 0)) == 1
assert projete(V(0, 0, 1), (0, 0, 1, 0)) != projete(V(1, 0, 1), (0, 0, 1, 0))

# geo-t-espace-projete-ex1 : distances ; leurres faux
assert dist(V(1, -2, 3), (2, -1, 2, -1)) == 3
assert {9, 1, Q(10, 3)}.isdisjoint({3}) and Q(2 + 2 + 6, 3) == Q(10, 3)  # leurre « oubli de d »
assert dist(V(1, 0, -2), (6, 2, -3, 1)) == Q(13, 7)  # v1
assert Q(12, 7) == Q(abs(6 + 0 + 6), 7) and {13, Q(13, 49), Q(12, 7)}.isdisjoint({Q(13, 7)})
assert dist(V(3, -1, 2), (1, -1, 1, -2)) == 4 * sqrt(3) / 3  # v2
assert all(simplify(x_ - 4 * sqrt(3) / 3) != 0 for x_ in [4, Q(4, 3), 2 * sqrt(3)])

# geo-t-espace-projete-ex2 : projetés ; leurres = M + kn (signe, symétrique, division par la norme), hors du plan
for M_, p, H_att, MH, leurres in [
    (V(2, -2, 3), (1, -2, 2, -3), V(1, 0, 1), 3, [V(3, -4, 5), V(0, 2, -1), V(-1, 4, -3)]),
    (V(4, 1, -2), (2, 1, -1, 1), V(0, -1, 0), 2 * sqrt(6), [V(8, 3, -4), V(-4, -3, 2), V(Q(1, 3), Q(-5, 6), Q(-1, 6))]),
    (V(6, 7, 3), (3, 0, 4, -5), V(3, 7, -1), 5, [V(9, 7, 7), V(0, 7, -5), V(-9, 7, -17)]),
]:
    assert projete(M_, p) == H_att and dist(M_, p) == MH and (M_ - H_att).norm() == MH
    for L_ in leurres:
        assert val(p, L_) != 0
# leurre v1 « oubli de d » : t = -11/6
assert V(4, 1, -2) + Q(-11, 6) * V(2, 1, -1) == V(Q(1, 3), Q(-5, 6), Q(-1, 6))

# geo-t-espace-projete-ex3 : cercle de centre H, rayon sqrt(R^2 - MH^2) ; leurres faux
for M_, p, H_att, MH, R_, r_ in [
    (V(1, 2, 3), (1, 1, 1, 0), V(-1, 0, 1), 2 * sqrt(3), 4, 2),
    (V(3, 0, 3), (2, -1, 2, -3), V(1, 1, 1), 3, 5, 4),
    (V(3, 4, 3), (1, 2, 2, 1), V(1, 0, -1), 6, 10, 8),
]:
    assert projete(M_, p) == H_att and dist(M_, p) == MH
    assert R_ > MH and sqrt(R_**2 - MH**2) == r_  # intersection non vide, rayon r
    assert r_ != R_ and sqrt(R_**2 + MH**2) != r_  # leurres « rayon R » et « rayon sqrt(R^2 + MH^2) »
assert sqrt(16 + 12) == 2 * sqrt(7) and sqrt(25 + 9) == sqrt(34) and sqrt(100 + 36) == 2 * sqrt(34)

# geo-t-espace-projete-ex4 : aire de ABC = 3V / delta, contrôlée par le produit vectoriel ; leurres faux
for (xa, yb, zc), n_att, d_att, delta, vol, aire, leurres in [
    ((2, 3, 6), V(3, 2, 1), -6, 3 * sqrt(14) / 7, 6, 3 * sqrt(14), [sqrt(14), 6 * sqrt(14), 3 * sqrt(14) / 7]),
    ((1, 2, 2), V(2, 1, 1), -2, sqrt(6) / 3, Q(2, 3), sqrt(6), [sqrt(6) / 3, 2 * sqrt(6), 3 * sqrt(6)]),
    ((4, 4, 2), V(1, 1, 2), -4, 2 * sqrt(6) / 3, Q(16, 3), 4 * sqrt(6), [4 * sqrt(6) / 3, 8 * sqrt(6), 2 * sqrt(6) / 3]),
]:
    A, B, C = V(xa, 0, 0), V(0, yb, 0), V(0, 0, zc)
    assert n_att.dot(B - A) == 0 and n_att.dot(C - A) == 0
    p = plan(n_att, A)
    assert p[3] == d_att and val(p, B) == 0 and val(p, C) == 0
    assert simplify(dist(V(0, 0, 0), p) - delta) == 0
    assert Q(1, 3) * Q(xa * yb, 2) * zc == vol
    assert simplify(3 * vol / delta - aire) == 0
    assert simplify((B - A).cross(C - A).norm() / 2 - aire) == 0
    assert simplify(vol / delta - leurres[0]) == 0  # leurre « V / delta » (oubli du facteur 3)
    for L_ in leurres:
        assert simplify(L_ - aire) != 0

print("geo_t_espace : toutes les vérifications passent.")
