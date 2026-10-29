"""Vérifications calculatoires de l'enrichissement du chapitre « geo-t-espace » (Terminale).

Questions de content/qcm/geo-t-espace.json et exercices ajoutés (droite-param-ex5, equation-plan-ex5, projete-ex5)
de content/exercices/geometrie/terminale-geo-t-espace.json. Repère orthonormé quand il le faut.
Chaque assert cite l'identifiant de la question ou de l'exercice (v1, v2 = variantes).
Lancer : python3 verification/enrichi_geo_t_espace.py (doit finir sans erreur).
"""

from sympy import Matrix, Rational, expand, solve, sqrt, symbols, Abs, linsolve

t, s, m = symbols("t s m", real=True)
V = lambda *c: Matrix(c)


def colineaires(u, v):
    return u.cross(v) == V(0, 0, 0)


def dans_plan(P, n, d):
    return n.dot(P) + d == 0


# --- QCM ---------------------------------------------------------------------
# geo-t-espace-q1 : plan (ABC), A(1,0,2), B(0,1,1), C(2,1,0) : x + 3y + 2z − 5 = 0 ; leurres faux.
A, B, C = V(1, 0, 2), V(0, 1, 1), V(2, 1, 0)
n = (B - A).cross(C - A)
assert colineaires(n, V(1, 3, 2))
assert all(dans_plan(P, V(1, 3, 2), -5) for P in (A, B, C))
for (nn, d) in ((V(1, 3, 2), 5), (V(1, -3, 2), -5), (V(3, 1, 2), -7)):
    assert not all(dans_plan(P, nn, d) for P in (A, B, C))
assert dans_plan(A, V(1, -3, 2), -5) and dans_plan(A, V(3, 1, 2), -7)  # leurres passant par A

# geo-t-espace-q2 : droite passant par A(1,−1,2) orthogonale à 2x − y + 3z − 1 = 0.
A, nP = V(1, -1, 2), V(2, -1, 3)


def droite_ok(P0, u):
    contient_A = colineaires(A - P0, u)
    return colineaires(u, nP) and contient_A


assert droite_ok(V(1, -1, 2), V(2, -1, 3))
assert not droite_ok(V(1, -1, 2), V(2, 1, 3))
assert not droite_ok(V(2, -1, 3), V(1, -1, 2))
assert not droite_ok(V(3, -2, 4), V(2, -1, 3))  # bonne direction, mais ne passe pas par A

# geo-t-espace-q3 (vrai) : A(1,0,0), B(0,1,0), C(0,0,1), D(1,1,−1) coplanaires (AD = AB − AC).
A, B, C, D = V(1, 0, 0), V(0, 1, 0), V(0, 0, 1), V(1, 1, -1)
assert D - A == (B - A) - (C - A)
assert Matrix.hstack(B - A, C - A, D - A).det() == 0

# geo-t-espace-q4 (faux) : D(t, 1+t, 2t), D'(1+s, s, 1−s) : ni sécantes ni parallèles.
assert solve([t - (1 + s), 1 + t - s, 2 * t - (1 - s)], [t, s]) == []
assert not colineaires(V(1, 1, 2), V(1, 1, -1))

# geo-t-espace-q5 : distance de A(2,−1,3) à x − 2y + 2z + 3 = 0 : 13/3 ; leurres faux.
d = Abs(2 + 2 + 6 + 3) / sqrt(1 + 4 + 4)
assert d == Rational(13, 3)
assert all(v != d for v in (13, Rational(13, 9), Rational(7, 3)))

# geo-t-espace-q6 : u(2,m,−1) colinéaire à n(4,−6,−2) ⟺ m = −3 ; 5/3 rend u orthogonal à n (droite parallèle au plan).
assert solve(list(V(2, m, -1).cross(V(4, -6, -2))), m) == {m: -3}
assert solve(V(2, m, -1).dot(V(4, -6, -2)), m) == [Rational(5, 3)]
for v in (3, Rational(5, 3), -6):
    assert not colineaires(V(2, v, -1), V(4, -6, -2))

# geo-t-espace-q7 (vrai) : x + y = 3 pour tout point (1+t, 2−t, 3t).
assert expand((1 + t) + (2 - t) - 3) == 0

# geo-t-espace-q8 (faux) : H(0,1,1) ∈ P : x+y+z−2 = 0, mais MH non colinéaire à n ; le projeté de M(2,1,2) est (1,0,1).
M, H, nP = V(2, 1, 2), V(0, 1, 1), V(1, 1, 1)
assert dans_plan(H, nP, -2) and not colineaires(H - M, nP)
Hv = V(1, 0, 1)
assert dans_plan(Hv, nP, -2) and colineaires(Hv - M, nP)

# --- geo-t-espace-droite-param-ex5 : point de D le plus proche de B ------------------------------------------------
def plus_proche(A0, u, B, t_min, H, dist2, d_plan, leurres):
    Mt = A0 + t * u
    assert solve(list(Mt - B), t) == []                       # Q1 : B ∉ D
    BM2 = expand((Mt - B).dot(Mt - B))
    a2 = u.dot(u)
    assert expand(BM2 - (a2 * (t - t_min) ** 2 + dist2)) == 0  # Q2-Q3 : forme canonique
    assert Mt.subs(t, t_min) == H                              # H
    assert (H - B).dot(u) == 0                                 # Q4
    assert dans_plan(B, u, d_plan) and solve(u.dot(Mt) + d_plan, t) == [t_min]  # Q5
    assert all(L != sqrt(dist2) for L in leurres)


# main : D(1+t, t, 2−t), B(2,4,1) ⇒ BM² = 3(t−2)² + 6, H(3,2,0), BH = √6 ; P : x + y − z − 5 = 0.
plus_proche(V(1, 0, 2), V(1, 1, -1), V(2, 4, 1), 2, V(3, 2, 0), 6, -5, [6, 3 * sqrt(2), sqrt(3)])
# v1 : D(2t, 1−t, −1+t), B(3,1,−1) ⇒ BM² = 6(t−1)² + 3, H(2,0,0), BH = √3 ; P : 2x − y + z − 4 = 0.
plus_proche(V(0, 1, -1), V(2, -1, 1), V(3, 1, -1), 1, V(2, 0, 0), 3, -4, [3, sqrt(6), sqrt(2)])
# v2 : D(2+t, −1+2t, 2t), B(3,−4,−2) ⇒ BM² = 9(t+1)² + 5, H(1,−3,−2), BH = √5 ; P : x + 2y + 2z + 9 = 0.
plus_proche(V(2, -1, 0), V(1, 2, 2), V(3, -4, -2), -1, V(1, -3, -2), 5, 9, [5, 3, sqrt(14)])

# --- geo-t-espace-equation-plan-ex5 : intersection droite-plan -------------------------------------------------------
def intersection(A0, u, n, d, I, leurres):
    assert u.dot(n) != 0                                       # sécants
    tt = solve(n.dot(A0 + t * u) + d, t)
    assert len(tt) == 1 and A0 + tt[0] * u == I
    for L in leurres:
        assert not dans_plan(L, n, d)


intersection(V(2, -1, 1), V(1, 2, -1), V(2, -1, 1), -9, V(-1, -7, 4), [V(5, 5, -2), V(-1, -7, -2), V(2, -1, 1)])
intersection(V(1, 0, 2), V(-1, 3, 1), V(1, 2, -1), -7, V(-1, 6, 4), [V(3, -6, 0), V(1, 0, 2), V(-1, 6, -4)])
intersection(V(0, 1, -2), V(3, -1, 2), V(2, 1, -1), -9, V(6, -1, 2), [V(-6, 3, -6), V(0, 1, -2), V(6, 1, 2)])

# --- geo-t-espace-projete-ex5 : symétrique par rapport à un plan ------------------------------------------------------
def symetrique(M, n, d, t0, H, MH2, Mp, leurres):
    tt = solve(n.dot(M + t * n) + d, t)
    assert tt == [t0] and M + t0 * n == H and dans_plan(H, n, d)
    assert (H - M).dot(H - M) == MH2
    assert 2 * H - M == Mp and (Mp + M) / 2 == H
    for L in leurres:
        assert L != Mp


# main : P : 2x − y + 2z − 3 = 0, M(1,−2,4) ⇒ H(−1,−1,2), MH = 3, M'(−3,0,0).
symetrique(V(1, -2, 4), V(2, -1, 2), -3, -1, V(-1, -1, 2), 9, V(-3, 0, 0), [V(-1, -1, 2), V(5, -4, 8), V(3, 0, 0)])
assert V(1, -2, 4) + 2 * V(2, -1, 2) == V(5, -4, 8)  # leurre : mauvais sens
# v1 : P : x + y + z − 3 = 0, M(3,4,5) ⇒ H(0,1,2), MH = 3√3, M'(−3,−2,−1).
symetrique(V(3, 4, 5), V(1, 1, 1), -3, -3, V(0, 1, 2), 27, V(-3, -2, -1), [V(0, 1, 2), V(9, 10, 11), V(3, 2, 1)])
# v2 : P : x − 2y + 2z + 4 = 0, M(1,−1,1) ⇒ H(0,1,−1), MH = 3, M'(−1,3,−3).
symetrique(V(1, -1, 1), V(1, -2, 2), 4, -1, V(0, 1, -1), 9, V(-1, 3, -3), [V(0, 1, -1), V(3, -5, 5), V(1, -3, 3)])

# --- Relecture indépendante (2026-10-29) : réponses retrouvées à l'aveugle, recalculées ------------------------------
import sympy as sp
_t = sp.symbols("t", real=True)
# geo-t-espace-droite-param-ex5 (et v1, v2) : minimum de BM_t², H, BH, plan orthogonal à D passant par B.
for A0, u, B, H, d2 in (((1, 0, 2), (1, 1, -1), (2, 4, 1), (3, 2, 0), 6), ((0, 1, -1), (2, -1, 1), (3, 1, -1), (2, 0, 0), 3),
                        ((2, -1, 0), (1, 2, 2), (3, -4, -2), (1, -3, -2), 5)):
    M = sp.Matrix(A0) + _t * sp.Matrix(u)
    D2 = sp.expand((M - sp.Matrix(B)).dot(M - sp.Matrix(B)))
    t0 = sp.solve(D2.diff(_t), _t)[0]
    assert tuple(M.subs(_t, t0)) == H and D2.subs(_t, t0) == d2
    assert all(M.subs(_t, t1) != sp.Matrix(B) for t1 in sp.solve(M[0] - B[0], _t))  # B ∉ D
    assert (sp.Matrix(H) - sp.Matrix(B)).dot(sp.Matrix(u)) == 0
# geo-t-espace-equation-plan-ex5 (et v1, v2) : point d'intersection.
for A0, u, n, d, I_ in (((2, -1, 1), (1, 2, -1), (2, -1, 1), -9, (-1, -7, 4)), ((1, 0, 2), (-1, 3, 1), (1, 2, -1), -7, (-1, 6, 4)),
                        ((0, 1, -2), (3, -1, 2), (2, 1, -1), -9, (6, -1, 2))):
    M = sp.Matrix(A0) + _t * sp.Matrix(u)
    assert sp.Matrix(n).dot(sp.Matrix(u)) != 0
    assert tuple(M.subs(_t, sp.solve(sp.Matrix(n).dot(M) + d, _t)[0])) == I_
# geo-t-espace-projete-ex5 (et v1, v2) : projeté H, distance, symétrique M'.
for Mp, n, d, H, dist, Ms in (((1, -2, 4), (2, -1, 2), -3, (-1, -1, 2), 3, (-3, 0, 0)),
                              ((3, 4, 5), (1, 1, 1), -3, (0, 1, 2), 3 * sp.sqrt(3), (-3, -2, -1)),
                              ((1, -1, 1), (1, -2, 2), 4, (0, 1, -1), 3, (-1, 3, -3))):
    M, N_ = sp.Matrix(Mp), sp.Matrix(n)
    k = -(N_.dot(M) + d) / N_.dot(N_)
    assert tuple(M + k * N_) == H and tuple(M + 2 * k * N_) == Ms
    assert sp.simplify(sp.Abs(N_.dot(M) + d) / sp.sqrt(N_.dot(N_)) - dist) == 0

print("enrichi_geo_t_espace : tout est vérifié.")
