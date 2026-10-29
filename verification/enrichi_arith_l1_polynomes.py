"""Vérifications calculatoires de l'enrichissement du chapitre « arith-l1-polynomes » : questions de QCM
(content/qcm/arith-l1-polynomes.json) et exercices ajoutés à content/exercices/arithmetique/l1-arith-l1-polynomes.json
(arith-l1-polynomes-division-ex5, arith-l1-polynomes-viete-ex5, arith-l1-polynomes-bezout-ex5, variantes comprises).

Chaque assert cite l'identifiant de la question ou de l'exercice concerné. Voir CLAUDE.md, « Contrôle de l'exactitude ».
Lancer : python3 verification/enrichi_arith_l1_polynomes.py (doit finir sans erreur).
"""

from sympy import (Poly, Rational, cos, diff, discriminant, expand, factor_list, gcd, pi, rem, div, roots, simplify,
                   sqrt, symbols, I, prod, real_roots, cbrt, together)

X = symbols("X")

# --- QCM -------------------------------------------------------------------------------------------------------
# arith-l1-polynomes-q1 : reste de X^100 par X² − 1 = 1 ; 100X − 99 est le reste par (X − 1)².
assert rem(X ** 100, X ** 2 - 1, X) == 1
assert rem(X ** 100, (X - 1) ** 2, X) == 100 * X - 99

# arith-l1-polynomes-q2 (vrai) : 1 racine de multiplicité exactement 2 de X³ − 3X + 2 = (X − 1)²(X + 2).
P = X ** 3 - 3 * X + 2
assert P.subs(X, 1) == 0 and diff(P, X).subs(X, 1) == 0 and diff(P, X, 2).subs(X, 1) == 6
assert expand((X - 1) ** 2 * (X + 2) - P) == 0

# arith-l1-polynomes-q3 : X³ + aX + 16 a une racine multiple (a réel) ssi a = −12 ; leurres 12, −16, −4 sans racine multiple.
a = symbols("a")
disc = discriminant(X ** 3 + a * X + 16, X)
assert expand(disc - (-4 * a ** 3 - 27 * 256)) == 0
assert [r for r in roots(Poly(disc, a)).keys() if r.is_real] == [-12]
assert expand((X - 2) ** 2 * (X + 4) - (X ** 3 - 12 * X + 16)) == 0
assert all(disc.subs(a, v) != 0 for v in (12, -16, -4))

# arith-l1-polynomes-q4 (faux) : X³ − 2 = (X − ∛2)(X² + ∛2 X + ∛4) avec un second facteur de discriminant < 0.
c2 = cbrt(2)
assert simplify(expand((X - c2) * (X ** 2 + c2 * X + c2 ** 2)) - (X ** 3 - 2)) == 0
assert simplify(c2 ** 2 - 4 * c2 ** 2) < 0
assert len(factor_list(X ** 3 - 2)[1]) == 1  # irréductible dans Q[X]

# arith-l1-polynomes-q5 : X^5 − 1 a 3 facteurs irréductibles unitaires dans R[X], 5 dans C[X].
F = (X - 1) * (X ** 2 - 2 * cos(2 * pi / 5) * X + 1) * (X ** 2 - 2 * cos(4 * pi / 5) * X + 1)
assert simplify(expand(F) - (X ** 5 - 1)) == 0
assert all(simplify(4 * cos(t) ** 2 - 4) < 0 for t in (2 * pi / 5, 4 * pi / 5))
assert len(real_roots(X ** 5 - 1)) == 1

# arith-l1-polynomes-q6 : σ3·σ1 = 6 pour X³ − 2X² + 5X − 3.
r = list(roots(X ** 3 - 2 * X ** 2 + 5 * X - 3, multiple=True))
assert len(r) == 3
s = simplify(r[0] ** 2 * r[1] * r[2] + r[0] * r[1] ** 2 * r[2] + r[0] * r[1] * r[2] ** 2)
assert abs(complex(s.evalf(30)) - 6) < 1e-20

# arith-l1-polynomes-q7 (faux) : PGCD de X³ + 1 et X² − 1 = X + 1.
assert gcd(X ** 3 + 1, X ** 2 - 1) == X + 1

# arith-l1-polynomes-q8 : P = 2X² − X + 1, P(2) = 7 ; leurres 5, 3, 11.
P = 2 * X ** 2 - X + 1
assert [P.subs(X, t) for t in (-1, 0, 1, 2)] == [4, 1, 2, 7]

# --- arith-l1-polynomes-division-ex5 (X^n modulo (X − 1)²(X + 1)) --------------------------------------------------
B = (X - 1) ** 2 * (X + 1)
assert diff(B, X).subs(X, 1) == 0 and B.subs(X, 1) == 0 and B.subs(X, -1) == 0


def reste_formule(m):
    if m % 2 == 0:
        return Rational(m, 2) * (X ** 2 - 1) + 1
    return Rational(m - 1, 2) * (X ** 2 - 1) + X


for m in range(0, 40):  # formule générale de la question 3
    R = rem(X ** m, B, X)
    assert expand(R - reste_formule(m)) == 0
    assert R.subs(X, 1) == 1 and diff(R, X).subs(X, 1) == m and R.subs(X, -1) == (-1) ** m
# Unicité (question 3) : le système a + b + c = 1, 2a + b = n, a − b + c = (−1)^n est de Cramer.
from sympy import Matrix
assert Matrix([[1, 1, 1], [2, 1, 0], [1, -1, 1]]).det() != 0
# exercice : n = 6 → 3X² − 2 ; leurres 6X − 5 (reste par (X − 1)²), 3X² + X − 3, X² faux.
assert expand(rem(X ** 6, B, X) - (3 * X ** 2 - 2)) == 0 and rem(X ** 6, (X - 1) ** 2, X) == 6 * X - 5
# variante 1 : n = 7 → 3X² + X − 3 ; leurres 7X − 6, 3X² − 2, X.
assert expand(rem(X ** 7, B, X) - (3 * X ** 2 + X - 3)) == 0 and rem(X ** 7, (X - 1) ** 2, X) == 7 * X - 6
# variante 2 : n = 9 → 4X² + X − 4 ; leurres 9X − 8, 4X² − 3, X.
assert expand(rem(X ** 9, B, X) - (4 * X ** 2 + X - 4)) == 0 and rem(X ** 9, (X - 1) ** 2, X) == 9 * X - 8
for m, R in ((6, 3 * X ** 2 - 2), (7, 3 * X ** 2 + X - 3), (9, 4 * X ** 2 + X - 4)):
    D = X ** m - R
    assert D.subs(X, 1) == 0 and diff(D, X).subs(X, 1) == 0 and D.subs(X, -1) == 0

# --- arith-l1-polynomes-viete-ex5 ((1 + x1)(1 + x2)(1 + x3)) -------------------------------------------------------


def produit(P):
    rs = roots(P, X, multiple=True)
    assert len(rs) == 3
    return complex(prod([1 + t for t in rs]).evalf(40))


# exercice : P = X³ − 4X² + 2X + 5 : 1 + 4 + 2 − 5 = 2, P(−1) = −2 ; leurres −2, 12, 4.
P = X ** 3 - 4 * X ** 2 + 2 * X + 5
assert abs(produit(P) - 2) < 1e-25 and P.subs(X, -1) == -2
# variante 1 : P = X³ + 3X² − X + 2 : 1 − 3 − 1 − 2 = −5, P(−1) = 5.
P = X ** 3 + 3 * X ** 2 - X + 2
assert abs(produit(P) + 5) < 1e-25 and P.subs(X, -1) == 5
# variante 2 : P = 2X³ − X² + 3X − 4 : σ1 = 1/2, σ2 = 3/2, σ3 = 2, produit 5, P(−1) = −10 = −2·5.
P = 2 * X ** 3 - X ** 2 + 3 * X - 4
assert abs(produit(P) - 5) < 1e-25 and P.subs(X, -1) == -10

# --- arith-l1-polynomes-bezout-ex5 (PGCD et racines communes) ------------------------------------------------------


def verifier_pgcd(A, B, D, q1, R1, q2, facA, facB, communes, leurres):
    assert div(A, B, X) == (q1, R1)  # première division
    assert div(B, R1, X) == (q2, 0)
    assert expand(gcd(A, B) - D) == 0  # PGCD unitaire
    assert expand(A - facA) == 0 and expand(B - facB) == 0
    assert set(roots(D, X).keys()) == set(communes)
    assert set(roots(A, X).keys()) & set(roots(B, X).keys()) == set(communes)
    for L in leurres:
        assert expand(L - D) != 0


# exercice : A = X^5 − X^4 + 2X³ − 2X² + X − 1, B = X^4 − 1 ; D = X³ − X² + X − 1 = (X − 1)(X² + 1).
verifier_pgcd(X ** 5 - X ** 4 + 2 * X ** 3 - 2 * X ** 2 + X - 1, X ** 4 - 1, X ** 3 - X ** 2 + X - 1,
              X - 1, 2 * X ** 3 - 2 * X ** 2 + 2 * X - 2, Rational(1, 2) * X + Rational(1, 2),
              (X - 1) * (X ** 2 + 1) ** 2, (X - 1) * (X + 1) * (X ** 2 + 1), [1, I, -I],
              [X - 1, X ** 2 + 1, 2 * X ** 3 - 2 * X ** 2 + 2 * X - 2])
assert div(X ** 4 - 1, X ** 3 - X ** 2 + X - 1, X) == (X + 1, 0)
# variante 1 : A = X^4 + 2X³ + 2X² + 2X + 1, B = X^4 − 1 ; D = X³ + X² + X + 1 = (X + 1)(X² + 1).
verifier_pgcd(X ** 4 + 2 * X ** 3 + 2 * X ** 2 + 2 * X + 1, X ** 4 - 1, X ** 3 + X ** 2 + X + 1,
              1, 2 * X ** 3 + 2 * X ** 2 + 2 * X + 2, Rational(1, 2) * X - Rational(1, 2),
              (X + 1) ** 2 * (X ** 2 + 1), (X - 1) * (X + 1) * (X ** 2 + 1), [-1, I, -I],
              [X + 1, X ** 2 + 1, 2 * X ** 3 + 2 * X ** 2 + 2 * X + 2])
assert div(X ** 4 - 1, X ** 3 + X ** 2 + X + 1, X) == (X - 1, 0)
# variante 2 : A = X^4 + X³ + 2X² + X + 1, B = X³ + X² + X + 1 ; D = X² + 1.
verifier_pgcd(X ** 4 + X ** 3 + 2 * X ** 2 + X + 1, X ** 3 + X ** 2 + X + 1, X ** 2 + 1,
              X, X ** 2 + 1, X + 1,
              (X ** 2 + 1) * (X ** 2 + X + 1), (X + 1) * (X ** 2 + 1), [I, -I],
              [X + 1, X ** 3 + X ** 2 + X + 1, 1])
assert discriminant(X ** 2 + X + 1, X) == -3

print("enrichi_arith_l1_polynomes : tout est vérifié.")
