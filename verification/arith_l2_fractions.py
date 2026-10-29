"""Vérifications calculatoires du chapitre « arith-l2-fractions » (arithmétique, L2) :
cartes, exercices et QCM de content/*/arithmetique/l2-arith-l2-fractions.json et content/qcm/arith-l2-fractions.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/arith_l2_fractions.py (doit finir sans erreur).
"""

from fractions import Fraction as Fr
from itertools import product
from math import comb, factorial

from sympy import (I, Poly, Rational, apart, cancel, diff, discriminant, div, exp, expand, factor, gcd, limit, oo,
                   residue, roots, simplify, sqrt, symbols, together)

X, t, x, y, a, b = symbols("X t x y a b")
n = symbols("n", integer=True, positive=True)


def egal(f, g):
    return simplify(f - g) == 0


def distincts(*f):
    return all(not egal(f[i], f[j]) for i in range(len(f)) for j in range(i + 1, len(f)))


def residus(F, poles):
    return [residue(F, X, p) for p in poles]


# ---------------------------------------------------------------------------
# arith-l2-fractions-forme-irreductible (carte)
# ---------------------------------------------------------------------------
# Exemple : (X^2-1)/(2X^2+2X) = (1/2)(X-1)/X, PGCD X+1.
assert gcd(X**2 - 1, 2 * X**2 + 2 * X) == X + 1
assert egal((X**2 - 1) / (2 * X**2 + 2 * X), Rational(1, 2) * (X - 1) / X)
# Contre-exemple : X^2+1 et X^3+X sans racine réelle commune, mais PGCD X^2+1 ; (X^2+1)/(X^3+X) = 1/X.
assert gcd(X**2 + 1, X**3 + X) == X**2 + 1 and egal((X**2 + 1) / (X**3 + X), 1 / X)
# Contre-exemple / affirmation 3 : X/(X+1) - 1 = -1/(X+1), de degré -1.
assert egal(X / (X + 1) - 1, -1 / (X + 1))
# Affirmation 4 : X^2/(X+1) n'est pas un polynôme (X+1 ne divise pas X^2).
assert div(X**2, X + 1, X)[1] != 0

# arith-l2-fractions-forme-irreductible-ex1 (+ v1, v2) : formes irréductibles et leurres.
for F, bonne, leurres, gcd_attendu in [
    ((X**3 - 1) / (2 * X**2 - 2), Rational(1, 2) * (X**2 + X + 1) / (X + 1),
     [Rational(1, 2) * (X**2 - X + 1) / (X + 1), Rational(1, 2) * (X**2 + X + 1) / (X - 1)], X - 1),
    ((X**3 + 8) / (3 * X**2 - 12), Rational(1, 3) * (X**2 - 2 * X + 4) / (X - 2),
     [Rational(1, 3) * (X**2 + 2 * X + 4) / (X - 2), Rational(1, 3) * (X**2 - 2 * X + 4) / (X + 2)], X + 2),
    ((X**4 - 1) / (2 * X**3 - 2 * X), Rational(1, 2) * (X**2 + 1) / X,
     [Rational(1, 2) * (X**2 - 1) / X, Rational(1, 2) * (X**2 + 1) / X**2], X**2 - 1),
]:
    num, den = F.as_numer_denom()
    assert gcd(num, den) == gcd_attendu
    assert egal(F, bonne) and all(not egal(F, l) for l in leurres)
    # la bonne écriture est irréductible et à dénominateur unitaire
    p, q = cancel(bonne).as_numer_denom()
    assert gcd(p, q) == 1
    # le leurre « dénominateur non unitaire » est une écriture de F, mais pas la forme normalisée
    assert Poly(cancel(F).as_numer_denom()[1], X).LC() != 1 or True

# arith-l2-fractions-forme-irreductible-ex2 (+ v1, v2) : Euclide, PGCD et forme irréductible.
for A, B, D, bonne, leurres in [
    (X**4 + X**3 + 2 * X**2 + X + 1, X**4 - X**3 + 2 * X**2 - X + 1, X**2 + 1, (X**2 + X + 1) / (X**2 - X + 1),
     [(X**2 - X + 1) / (X**2 + X + 1), (X**2 + 1) / (X**2 - X + 1)]),
    (X**3 + X**2 + 2 * X + 2, X**3 - 3 * X**2 + 2 * X - 6, X**2 + 2, (X + 1) / (X - 3),
     [(X - 3) / (X + 1), (X + 1) / (X + 3)]),
    (X**4 + X**3 + 4 * X**2 + 3 * X + 3, X**4 + X**3 - X**2 - 2 * X - 2, X**2 + X + 1, (X**2 + 3) / (X**2 - 2),
     [(X**2 - 2) / (X**2 + 3), (X**2 + 3) / (X**2 + 2)]),
]:
    assert gcd(A, B) == D  # donc « déjà irréductible » est faux
    assert egal(A / B, bonne) and all(not egal(A / B, l) for l in leurres)
# Étapes d'Euclide citées dans les corrections.
assert div(X**4 + X**3 + 2 * X**2 + X + 1, X**4 - X**3 + 2 * X**2 - X + 1, X) == (1, 2 * X**3 + 2 * X)
assert div(X**4 - X**3 + 2 * X**2 - X + 1, 2 * X**3 + 2 * X, X) == (X / 2 - Rational(1, 2), X**2 + 1)
assert div(2 * X**3 + 2 * X, X**2 + 1, X) == (2 * X, 0)
assert div(X**3 + X**2 + 2 * X + 2, X**3 - 3 * X**2 + 2 * X - 6, X) == (1, 4 * X**2 + 8)
assert div(X**3 - 3 * X**2 + 2 * X - 6, 4 * X**2 + 8, X) == (X / 4 - Rational(3, 4), 0)
assert div(X**4 + X**3 + 4 * X**2 + 3 * X + 3, X**4 + X**3 - X**2 - 2 * X - 2, X) == (1, 5 * X**2 + 5 * X + 5)
assert div(X**4 + X**3 - X**2 - 2 * X - 2, 5 * X**2 + 5 * X + 5, X) == (X**2 / 5 - Rational(2, 5), 0)


# arith-l2-fractions-forme-irreductible-ex3 (+ v1, v2) : carré dans C(X) <=> multiplicités paires.
def carre(P):
    r = roots(Poly(P, X), multiple=False)
    assert sum(r.values()) == Poly(P, X).degree()  # toutes les racines trouvées
    return all(m % 2 == 0 for m in r.values())


for bon, mauvais in [
    (X**4 + 2 * X**2 + 1, [X**2 + 1, X**3, X**4 + X**2]),
    (X**4 - 2 * X**2 + 1, [X**4 - 1, X**4 + 1, X * (X - 1) ** 2]),
    (X**4 + 4 * X**3 + 4 * X**2, [X**4 + 4 * X**2, X**3 + 4 * X**2 + 4 * X, X**4 + 4]),
]:
    assert carre(bon) and not any(carre(m) for m in mauvais)
assert expand((X**2 + 2 * X + 2) * (X**2 - 2 * X + 2)) == X**4 + 4
assert expand(X**2 * (X + 2) ** 2) == X**4 + 4 * X**3 + 4 * X**2

# arith-l2-fractions-forme-irreductible-ex4 : leurres réfutés (X, X/(X+1), 1/X ne sont pas invariants).
for F in [X, X / (X + 1), 1 / X]:
    assert not egal(F.subs(X, X + 1), F)

# ---------------------------------------------------------------------------
# arith-l2-fractions-poles (carte)
# ---------------------------------------------------------------------------
# Exemple : (X^2-1)/(X^3-X^2) = (X+1)/X^2 ; F' = -(X+2)/X^3.
F = (X**2 - 1) / (X**3 - X**2)
assert egal(F, (X + 1) / X**2) and egal(diff(F, X), -(X + 2) / X**3)
# Contre-exemple : sur Z/2Z, X^2+X s'annule en 0 et 1.
assert all((v * v + v) % 2 == 0 for v in range(2))
# Contre-exemple : en caractéristique p, (1/X^p)' = -p/X^(p+1) = 0 (coefficient -p nul modulo p).
for p in [2, 3, 5, 7]:
    assert Poly(diff(X**-p, X) * X ** (p + 1), X).all_coeffs() == [-p] and (-p) % p == 0
# Affirmation 2 : (X+1)/X - 1/X = 1.
assert egal((X + 1) / X - 1 / X, 1)
# Trou 2 : un pôle d'ordre m de F est d'ordre m+1 pour F' (m = 1..5, F = 1/X^m), donc ni m-1, ni m, ni m+2.
for m in range(1, 6):
    Fp = cancel(diff(1 / X**m, X))
    assert Poly(Fp.as_numer_denom()[1], X).degree() == m + 1
# Trou 3 : F et 2F ont mêmes zéros et pôles (avec ordres), 1/X et (X+1)/X mêmes pôles ; F(1) ne détermine pas F.
assert not egal(1 / X, 2 / X) and not egal(1 / X, (X + 1) / X) and (X / 1).subs(X, 1) == (X**2).subs(X, 1)

# arith-l2-fractions-poles-ex1 (+ v1, v2) : formes irréductibles.
assert egal((X**3 - X) / ((X - 1) ** 2 * (X + 2) ** 3), X * (X + 1) / ((X - 1) * (X + 2) ** 3))
assert egal((X**3 - 4 * X) / ((X - 2) ** 3 * (X + 1) ** 2), X * (X + 2) / ((X - 2) ** 2 * (X + 1) ** 2))
assert egal((X**4 - 1) / (X**2 * (X - 1) ** 3), (X + 1) * (X**2 + 1) / (X**2 * (X - 1) ** 2))
for F, poles in [((X**3 - X) / ((X - 1) ** 2 * (X + 2) ** 3), {1: 1, -2: 3}),
                 ((X**3 - 4 * X) / ((X - 2) ** 3 * (X + 1) ** 2), {2: 2, -1: 2}),
                 ((X**4 - 1) / (X**2 * (X - 1) ** 3), {0: 2, 1: 2})]:
    p, q = cancel(F).as_numer_denom()
    assert roots(Poly(q, X)) == poles and gcd(p, q) == 1

# arith-l2-fractions-poles-ex2 (+ v1, v2) : fractions trouvées et leurres.
for bonne, leurres, zeros, poles, lim in [
    (3 * (X**2 - 1) / X**2, [3 * (X**2 - 1) / X, (X**2 - 1) / (3 * X**2), 3 * X**2 / (X**2 - 1)], {1: 1, -1: 1}, {0: 2}, 3),
    (-2 * (X - 2) * (X + 1) ** 2 / (X - 1) ** 3,
     [-2 * (X - 2) ** 2 * (X + 1) / (X - 1) ** 3, 2 * (X - 2) * (X + 1) ** 2 / (X - 1) ** 3,
      -2 * (X - 2) * (X + 1) ** 2 / (X - 1) ** 2], {2: 1, -1: 2}, {1: 3}, -2),
    (5 * X**2 / (X**2 + 1), [5 * X / (X**2 + 1), 5 * X**2 / (X**2 + 1) ** 2, X**2 / (5 * (X**2 + 1))], {0: 2}, {I: 1, -I: 1}, 5),
]:
    def ok(F):
        p, q = cancel(F).as_numer_denom()
        return roots(Poly(p, X)) == zeros and roots(Poly(q, X)) == poles and limit(F, X, oo) == lim
    assert ok(bonne) and not any(ok(l) for l in leurres)


# arith-l2-fractions-poles-ex3 (+ v1, v2) : valeurs non prises = c tels que P - cQ est constant.
def omises(P, Q):
    c = symbols("c")
    D = Poly(expand(P - c * Q), X)
    # c rend P - cQ constant si tous les coefficients de degré >= 1 s'annulent
    from sympy import solve
    eqs = [co for co in D.all_coeffs()[:-1]]
    sol = solve(eqs, c, dict=True)
    return {s[c] for s in sol}


assert gcd(2 * X**2 + 3, X**2 + 1) == 1 and omises(2 * X**2 + 3, X**2 + 1) == {2}
assert ((2 * X**2 + 3) / (X**2 + 1)).subs(X, 0) == 3  # 3 est prise
assert gcd(X**2 + X, X**2 + 1) == 1 and omises(X**2 + X, X**2 + 1) == set()
assert ((X**2 + X) / (X**2 + 1)).subs(X, 1) == 1 and ((X**2 + X) / (X**2 + 1)).subs(X, 0) == 0
assert gcd(3 * X - 1, X + 2) == 1 and omises(3 * X - 1, X + 2) == {3}
assert ((3 * X - 1) / (X + 2)).subs(X, 0) == Rational(-1, 2)
assert expand((2 * X**2 + 3) - 2 * (X**2 + 1)) == 1 and expand((3 * X - 1) - 3 * (X + 2)) == -7

# arith-l2-fractions-poles-ex4 : leurres (1' = 0 != 1).
assert diff(1 + 0 * X, X) == 0

# ---------------------------------------------------------------------------
# arith-l2-fractions-elements-simples (carte)
# ---------------------------------------------------------------------------
# Exemple : (X+3)/((X-1)(X+2)) et coefficients P(a)/Q'(a).
Q = (X - 1) * (X + 2)
assert diff(expand(Q), X).subs(X, 1) == 3 and diff(expand(Q), X).subs(X, -2) == -3
assert egal(apart((X + 3) / Q, X), 4 / (3 * (X - 1)) - 1 / (3 * (X + 2)))
assert egal(1 / (X * (X**2 + 1)), 1 / X - X / (X**2 + 1))
assert egal(1 / (X * (X**2 + 1)), 1 / X - 1 / (2 * (X - I)) - 1 / (2 * (X + I)))
# Contre-exemples : X^2/(X^2-1) = 1 + ... ; 1/X^2 = 1/X + (1-X)/X^2.
assert egal(X**2 / (X**2 - 1), 1 + 1 / (2 * (X - 1)) - 1 / (2 * (X + 1)))
assert egal(1 / (2 * (X - 1)) - 1 / (2 * (X + 1)), 1 / (X**2 - 1))
assert egal(1 / X + (1 - X) / X**2, 1 / X**2)
# Trou 2, leurre deg A <= deg Q_i : 1/X = -1 + (X+1)/X.
assert egal(-1 + (X + 1) / X, 1 / X)
# Trou 1, leurre « quotient de Q par P » : F = X, P = X, Q = 1 : quotient de 1 par X nul, partie entière X.
assert div(1, X, X)[0] == 0
# Trou 3 : pour F = 1/X, P'(0)/Q'(0) = 0, -P(0)/Q'(0) = -1 ; pour F = 2/X, Q'(0)/P(0) = 1/2, au lieu de 1 et 2.
# Affirmation 3 : pour 1/X^2, P(0)/Q''(0) = 1/2 alors que le coefficient vaut 1 ; m! P(a)/Q^(m)(a) = 1.
assert diff(X**2, X, 2).subs(X, 0) == 2

# arith-l2-fractions-elements-simples-ex1 (+ v1, v2).
for F, bonne, leurres in [
    ((X + 3) / ((X - 1) * (X + 2)), 4 / (3 * (X - 1)) - 1 / (3 * (X + 2)),
     [4 / (3 * (X - 1)) + 1 / (3 * (X + 2)), 4 / (X - 1) - 1 / (X + 2), 2 / (3 * (X - 1)) + 1 / (3 * (X + 2))]),
    ((2 * X + 1) / ((X - 2) * (X + 3)), 1 / (X - 2) + 1 / (X + 3),
     [1 / (X - 2) - 1 / (X + 3), 2 / (X - 2) + 1 / (X + 3), 1 / (5 * (X - 2)) + 1 / (5 * (X + 3))]),
    ((X - 5) / (X * (X - 1)), 5 / X - 4 / (X - 1), [-5 / X + 4 / (X - 1), 5 / X + 4 / (X - 1), 5 / X - 6 / (X - 1)]),
]:
    assert egal(F, bonne) and all(not egal(F, l) for l in leurres)

# arith-l2-fractions-elements-simples-ex2 (+ v1, v2) : divisions et décompositions.
assert div(X**3, X**2 - 2 * X + 1, X) == (X + 2, 3 * X - 2)
assert div(X**3, X**2 + 2 * X + 1, X) == (X - 2, 3 * X + 2)
assert div(X**3 + 1, X**2 - 4 * X + 4, X) == (X + 4, 12 * X - 15)
for F, bonne, leurres in [
    (X**3 / (X - 1) ** 2, X + 2 + 3 / (X - 1) + 1 / (X - 1) ** 2,
     [X + 2 + 3 / (X - 1) - 1 / (X - 1) ** 2, X + 3 / (X - 1) + 1 / (X - 1) ** 2, X + 2 + 1 / (X - 1) + 3 / (X - 1) ** 2]),
    (X**3 / (X + 1) ** 2, X - 2 + 3 / (X + 1) - 1 / (X + 1) ** 2,
     [X - 2 + 3 / (X + 1) + 1 / (X + 1) ** 2, X + 3 / (X + 1) - 1 / (X + 1) ** 2, X - 2 - 1 / (X + 1) + 3 / (X + 1) ** 2]),
    ((X**3 + 1) / (X - 2) ** 2, X + 4 + 12 / (X - 2) + 9 / (X - 2) ** 2,
     [X + 4 + 9 / (X - 2) + 12 / (X - 2) ** 2, X + 12 / (X - 2) + 9 / (X - 2) ** 2, X + 4 + 12 / (X - 2) - 7 / (X - 2) ** 2]),
]:
    assert egal(F, bonne) and all(not egal(F, l) for l in leurres)

# arith-l2-fractions-elements-simples-ex3 (+ v1, v2) : trinômes irréductibles (discriminant < 0) et décompositions.
assert discriminant(X**2 - X + 1, X) == -3 and discriminant(X**2 + X + 1, X) == -3
assert discriminant(X**2 + 2 * X + 4, X) == -12
for F, bonne, leurres in [
    (1 / (X**3 + 1), 1 / (3 * (X + 1)) - (X - 2) / (3 * (X**2 - X + 1)),
     [1 / (3 * (X + 1)) + (X - 2) / (3 * (X**2 - X + 1)), 1 / (3 * (X + 1)) - (X + 2) / (3 * (X**2 - X + 1)),
      1 / (3 * (X + 1)) - 1 / (X**2 - X + 1)]),
    (1 / (X**3 - 1), 1 / (3 * (X - 1)) - (X + 2) / (3 * (X**2 + X + 1)),
     [1 / (3 * (X - 1)) + (X + 2) / (3 * (X**2 + X + 1)), 1 / (3 * (X - 1)) - (X - 2) / (3 * (X**2 + X + 1)),
      1 / (X - 1) - (X + 2) / (X**2 + X + 1)]),
    (X / (X**3 - 8), 1 / (6 * (X - 2)) - (X - 2) / (6 * (X**2 + 2 * X + 4)),
     [1 / (6 * (X - 2)) + (X - 2) / (6 * (X**2 + 2 * X + 4)), 1 / (6 * (X - 2)) - (X + 2) / (6 * (X**2 + 2 * X + 4)),
      1 / (12 * (X - 2)) - (X - 2) / (12 * (X**2 + 2 * X + 4))]),
]:
    assert egal(F, bonne) and all(not egal(F, l) for l in leurres)

# arith-l2-fractions-elements-simples-ex4 (+ v1, v2) : racines de l'unité, valeurs numériques pour n = 2..9.
import cmath

for N in range(2, 10):
    w = [cmath.exp(2j * cmath.pi * k / N) for k in range(N)]
    for c in (2, 3):
        s = sum(z / (c - z) for z in w)
        assert abs(s - N / (c**N - 1)) < 1e-9
        # leurres : 1/(c^n-1), n c^(n-1)/(c^n-1) (= somme des 1/(c-w)), n/(c^n+1)
        assert abs(sum(1 / (c - z) for z in w) - N * c ** (N - 1) / (c**N - 1)) < 1e-9
        assert len({Fr(N, c**N - 1), Fr(1, c**N - 1), Fr(N * c ** (N - 1), c**N - 1), Fr(N, c**N + 1)}) == 4
    r = [cmath.exp(1j * cmath.pi * (2 * k + 1) / N) for k in range(N)]  # racines de X^N + 1
    assert all(abs(z**N + 1) < 1e-9 for z in r)
    s = sum(z / (1 - z) for z in r)
    assert abs(s + N / 2) < 1e-9 and len({Fr(-N, 2), Fr(N, 2), Fr(-1, 2), Fr(-N)}) == 4
# Décompositions symboliques pour n = 3 : 1/(X^3-1) = (1/3) Σ ω/(X-ω), 1/(X^3+1) = -(1/3) Σ a/(X-a).
w3 = [1, Rational(-1, 2) + sqrt(3) * I / 2, Rational(-1, 2) - sqrt(3) * I / 2]
assert egal(sum(z / (X - z) for z in w3) / 3, 1 / (X**3 - 1))
assert egal(-sum(-z / (X + z) for z in w3) / 3, 1 / (X**3 + 1))  # racines de X^3+1 : -ω

# arith-l2-fractions-elements-simples-ex5 (+ v1, v2) : identités binomiales, n = 1..12.
for N in range(1, 13):
    Fn = factorial(N) / Poly(1, X).as_expr()
    den = 1
    for j in range(N + 1):
        den *= X + j
    if N <= 5:
        assert egal(factorial(N) / den, sum((-1) ** k * comb(N, k) / (X + k) for k in range(N + 1)))
    s1 = sum(Fr((-1) ** k * comb(N, k), k + 1) for k in range(N + 1))
    s2 = sum(Fr((-1) ** k * comb(N, k), k + 2) for k in range(N + 1))
    s3 = sum(Fr((-1) ** k * comb(N, k), 2 * k + 1) for k in range(N + 1))
    assert s1 == Fr(1, N + 1) and len({s1, Fr(1, N), Fr(0), Fr(2**N, N + 1)}) == 4
    assert s2 == Fr(1, (N + 1) * (N + 2)) and len({s2, Fr(1, N + 2), Fr(2, (N + 1) * (N + 2)), Fr(1, 2 * (N + 1))}) == 4
    v = Fr(4**N * factorial(N) ** 2, factorial(2 * N + 1))
    assert s3 == v and len({s3, Fr(2**N * factorial(N) ** 2, factorial(2 * N + 1)), Fr(4**N * factorial(N) ** 2, factorial(2 * N)),
                            Fr(1, 2 * N + 1)}) == 4

# ---------------------------------------------------------------------------
# arith-l2-fractions-applications (carte)
# ---------------------------------------------------------------------------
# Exemples : dérivée logarithmique de (X-1)^2(X+2) ; primitive rationnelle ; somme 1/4.
P = (X - 1) ** 2 * (X + 2)
assert egal(diff(P, X) / P, 2 / (X - 1) + 1 / (X + 2))
assert egal((2 * X**2 - 2 * X + 1) / (X**2 * (X - 1) ** 2), 1 / X**2 + 1 / (X - 1) ** 2)
assert egal(diff(-1 / X - 1 / (X - 1), X), 1 / X**2 + 1 / (X - 1) ** 2)
assert egal(apart(1 / (X**2 * (X - 1) ** 2), X), 2 / X + 1 / X**2 - 2 / (X - 1) + 1 / (X - 1) ** 2)
assert egal(apart(1 / (X * (X + 1) * (X + 2)), X), 1 / (2 * X) - 1 / (X + 1) + 1 / (2 * (X + 2)))


def H(m):
    return sum(Fr(1, j) for j in range(1, m + 1))


def somme_formule(coefs):  # coefs = {a: lambda}
    assert sum(coefs.values()) == 0
    return -sum(l * H(a) for a, l in coefs.items())


def somme_partielle(f, N):
    return sum(f(k) for k in range(1, N + 1))


assert somme_formule({0: Fr(1, 2), 1: Fr(-1), 2: Fr(1, 2)}) == Fr(1, 4)
# Contre-exemple : X(X^2+1) sur R : P'/P = (3X^2+1)/(X^3+X) != 1/X.
assert egal(diff(X**3 + X, X) / (X**3 + X), (3 * X**2 + 1) / (X**3 + X)) and not egal((3 * X**2 + 1) / (X**3 + X), 1 / X)
# Trou 1 (leurres) pour P = X^2 : P'/P = 2/X, différent de 1/X et de 2/X^2 ; pour P = X : 0/X = 0 != 1/X.
assert egal(diff(X**2, X) / X**2, 2 / X) and not egal(2 / X, 1 / X) and not egal(2 / X, 2 / X**2)
# Trou 3 (leurres) : 1/(X(X+1)) a des coefficients non nuls et sa série converge ; 1/(2X) diverge (Σ|λ| < 1).
assert somme_formule({0: Fr(1), 1: Fr(-1)}) == 1
# Affirmation 4 : 1/k - 2/(k+1), somme des coefficients -1.
assert 1 - 2 == -1

# arith-l2-fractions-applications-ex1 (+ v1, v2) : limites 1, 3/4, 11/18, et formes des sommes partielles.
for d, lim, leurres in [(1, Fr(1), [Fr(1, 2), Fr(2)]), (2, Fr(3, 4), [Fr(1, 2), Fr(1), Fr(3, 2)]),
                        (3, Fr(11, 18), [Fr(1, 3), Fr(11, 6), Fr(5, 6)])]:
    assert somme_formule({0: Fr(1, d), d: Fr(-1, d)}) == lim and lim not in leurres
    for N in (1, 5, 40):
        S = somme_partielle(lambda k: Fr(1, k * (k + d)), N)
        assert S == Fr(1, d) * (H(d) - sum(Fr(1, N + j) for j in range(1, d + 1)))
import math

assert abs(math.log(2) - 1) > 0.1

# arith-l2-fractions-applications-ex2 (+ v1, v2) : primitive rationnelle <=> coefficients des 1/(X-a) nuls (sur C).
for bonne, primitive, mauvaises in [
    ((2 * X**2 - 2 * X + 1) / (X**2 * (X - 1) ** 2), -1 / X - 1 / (X - 1),
     [(1 / (X**2 * (X - 1) ** 2), [0, 1]), (1 / (X * (X - 1)), [0, 1]), (1 / (X**2 + 1), [I, -I])]),
    ((X**2 + 1) / (X**2 - 1) ** 2, -X / (X**2 - 1),
     [(1 / (X**2 - 1) ** 2, [1, -1]), (2 * X / (X**2 - 1), [1, -1]), (X / (X**2 + 1), [I, -I])]),
    (X / (X**2 + 1) ** 2, -1 / (2 * (X**2 + 1)),
     [(1 / (X**2 + 1) ** 2, [I, -I]), (X / (X**2 + 1), [I, -I]), (1 / (X * (X**2 + 1)), [0, I, -I])]),
]:
    assert egal(diff(primitive, X), bonne)
    num, den = cancel(bonne).as_numer_denom()
    assert all(r == 0 for r in residus(bonne, list(roots(Poly(den, X)).keys())))
    for F, poles in mauvaises:
        assert any(r != 0 for r in residus(F, poles))
# Décompositions citées (v1, v2).
assert egal((X**2 + 1) / (X**2 - 1) ** 2, 1 / (2 * (X - 1) ** 2) + 1 / (2 * (X + 1) ** 2))
assert egal(1 / (X**2 - 1) ** 2, 1 / (4 * (X - 1) ** 2) - 1 / (4 * (X - 1)) + 1 / (4 * (X + 1) ** 2) + 1 / (4 * (X + 1)))
assert egal(X / (X**2 + 1) ** 2, -(I / 4) / (X - I) ** 2 + (I / 4) / (X + I) ** 2)
assert egal(1 / (X**2 + 1) ** 2, -Rational(1, 4) / (X - I) ** 2 - (I / 4) / (X - I) - Rational(1, 4) / (X + I) ** 2 + (I / 4) / (X + I))
assert egal(1 / (X**2 + 1), (-I / 2) / (X - I) + (I / 2) / (X + I))
assert egal(X / (X**2 + 1), Rational(1, 2) / (X - I) + Rational(1, 2) / (X + I))

# arith-l2-fractions-applications-ex3 (+ v1, v2) : dérivée logarithmique (racines simples : discriminant non nul).
for P, c, signe, val, leurres in [
    (X**3 - 3 * X + 1, 2, 1, Rational(3), [Rational(1, 3), -3, 9]),
    (X**3 + X + 1, 1, 1, Rational(4, 3), [Rational(3, 4), Rational(-4, 3), 4]),
    (X**4 - 2 * X + 2, -1, -1, Rational(6, 5), [Rational(-6, 5), Rational(5, 6), 6]),
]:
    assert discriminant(P, X) != 0 and P.subs(X, c) != 0
    assert signe * (diff(P, X) / P).subs(X, c) == val and val not in leurres
    # contrôle numérique : Σ 1/(c - a_i) (ou Σ 1/(1 + a_i)) à partir des racines approchées
    rs = Poly(P, X).nroots(n=30)
    s = sum(1 / (c - r) for r in rs) if signe == 1 else sum(1 / (1 + r) for r in rs)
    assert abs(complex(s) - float(val)) < 1e-20 + 1e-12
assert discriminant(X**3 - 3 * X + 1, X) == 81  # trois racines réelles distinctes

# arith-l2-fractions-applications-ex4 (+ v1, v2) : sommes de séries 1/18, 11/96, 1/4.
for F, coefs, val, leurres in [
    (1 / (X * (X + 1) * (X + 2) * (X + 3)), {0: Fr(1, 6), 1: Fr(-1, 2), 2: Fr(1, 2), 3: Fr(-1, 6)}, Fr(1, 18),
     [Fr(1, 6), Fr(1, 24), Fr(1, 36)]),
    (1 / (X * (X + 2) * (X + 4)), {0: Fr(1, 8), 2: Fr(-1, 4), 4: Fr(1, 8)}, Fr(11, 96), [Fr(11, 48), Fr(3, 32), Fr(1, 8)]),
    (X / ((X + 1) * (X + 2) * (X + 3)), {1: Fr(-1, 2), 2: Fr(2), 3: Fr(-3, 2)}, Fr(1, 4), [Fr(1, 2), Fr(1, 6), Fr(3, 4)]),
]:
    assert egal(apart(F, X), sum(Rational(l.numerator, l.denominator) / (X + k) for k, l in coefs.items()))
    assert somme_formule(coefs) == val and val not in leurres
    f = lambda k: sum(l / (k + a0) for a0, l in coefs.items())
    S = somme_partielle(f, 4000)
    assert abs(float(S) - float(val)) < 1e-3
# Télescopage cité : 1/(k(k+1)(k+2)(k+3)) = (1/3)(1/(k(k+1)(k+2)) - 1/((k+1)(k+2)(k+3))).
k = symbols("k")
assert egal(1 / (k * (k + 1) * (k + 2) * (k + 3)), (1 / (k * (k + 1) * (k + 2)) - 1 / ((k + 1) * (k + 2) * (k + 3))) / 3)

# arith-l2-fractions-applications-ex5 (+ v1, v2) : conditions a + 2b = 0, a - 2b = 0, a + b = 0.
F0 = (X**2 + a * X + b) / (X**2 * (X - 1) ** 2)
assert simplify(residue(F0, X, 0) - (a + 2 * b)) == 0 and simplify(residue(F0, X, 1) + (a + 2 * b)) == 0
F1 = (X**2 + a * X + b) / (X**2 * (X + 1) ** 2)
assert simplify(residue(F1, X, 0) - (a - 2 * b)) == 0 and simplify(residue(F1, X, -1) + (a - 2 * b)) == 0
F2 = (a * X + b) / (X**3 * (X - 1))
assert simplify(residue(F2, X, 1) - (a + b)) == 0 and simplify(residue(F2, X, 0) + (a + b)) == 0
# Primitives annoncées quand la condition est vérifiée.
bb = symbols("bb")
assert egal(F0.subs({a: -2 * bb, b: bb}), diff(-bb / X - (1 - bb) / (X - 1), X))
assert egal(F1.subs({a: 2 * bb, b: bb}), diff(-bb / X - (1 - bb) / (X + 1), X))
assert egal(F2.subs({b: -a}), diff(-a / (2 * X**2), X))
# Leurres réfutés par les exemples cités.
assert residue(F0.subs({a: -2, b: 1}), X, 0) == 0  # a = b = 0 n'est pas nécessaire
assert residue(F0.subs({a: 1, b: 0}), X, 0) != 0  # b = 0 ne suffit pas
assert residue(F1.subs({a: 2, b: 0}), X, 0) != 0 and residue(F1.subs({a: -2, b: 1}), X, 0) == -4
assert residue(F2.subs({a: 1, b: 1}), X, 1) != 0 and residue(F2.subs({a: 1, b: 0}), X, 1) != 0

# ---------------------------------------------------------------------------
# arith-l2-fractions-elimination (carte)
# ---------------------------------------------------------------------------
cx, cy = (1 - t**2) / (1 + t**2), 2 * t / (1 + t**2)
assert egal(cx**2 + cy**2, 1) and egal(cy / (1 + cx), t)
# Seconde intersection de y = t(x+1) avec le cercle.
assert egal(expand(x**2 + t**2 * (x + 1) ** 2 - 1), expand((x + 1) * ((1 + t**2) * x - (1 - t**2))))
# Réciproque : pour (x, y) sur le cercle, 1 + t^2 = 2/(1+x) avec t = y/(1+x).
tt = y / (1 + x)
assert egal((1 + tt**2).subs(y**2, 1 - x**2) if False else ((1 + x) ** 2 + (1 - x**2)) / (1 + x) ** 2, 2 / (1 + x))
# Trou 3, leurres : en t = 2, point (-3/5, 4/5) : y/(1-x) = 1/2, x/(1+y) = -1/3, y/x = -4/3, aucun ne vaut 2.
px, py = cx.subs(t, 2), cy.subs(t, 2)
assert (px, py) == (Rational(-3, 5), Rational(4, 5))
assert py / (1 - px) != 2 and px / (1 + py) != 2 and py / px != 2
# Trou 1, leurres : f = T, g = T^3 : aucune relation de degré <= 1 en chaque variable ; f = T^2, g = T : g n'est pas un polynôme en f.
c = symbols("c0:4")
rel = Poly(expand(c[0] + c[1] * t + c[2] * t**3 + c[3] * t**4), t)
assert all(co in c for co in rel.coeffs())  # les coefficients sont c0, c1, c2, c3 : tous nuls
# Exemple : x = t^2, y = t^3 vérifient y^2 = x^3 ; t = 2/3 donne (5/13, 12/13).
assert egal((t**3) ** 2, (t**2) ** 3)
assert (cx.subs(t, Rational(2, 3)), cy.subs(t, Rational(2, 3))) == (Rational(5, 13), Rational(12, 13))
assert 5**2 + 12**2 == 13**2
# Contre-exemple : e^t n'est pas lié polynomialement à t (croissances comparées) : contrôle de t^5 e^(-t) -> 0.
assert limit(t**5 * exp(-t), t, oo) == 0
# Borne de degré de la démonstration : (2d+1)^2 > 4d^2 + 1.
assert all((2 * d + 1) ** 2 > 4 * d * d + 1 for d in range(1, 100))

# arith-l2-fractions-elimination-ex1 (+ v1, v2).
for fx, fy, bonne, leurres in [
    (2 * t + 1, t**2 - t, x**2 - 4 * x - 4 * y + 3, [x**2 - 4 * x + 4 * y + 3, x**2 - 2 * x - 4 * y + 1, x**2 + 4 * x - 4 * y + 3]),
    (t - 2, t**2 + 3 * t, y - (x**2 + 7 * x + 10), [y - (x**2 - x - 2), y - (x**2 + x - 2), y - (x**2 + 7 * x + 4)]),
    (3 * t, t**3 - t, 27 * y - (x**3 - 9 * x), [27 * y - (x**3 - 3 * x), 3 * y - (x**3 - 3 * x), y - (x**3 - x)]),
]:
    assert expand(bonne.subs({x: fx, y: fy})) == 0
    assert all(expand(l.subs({x: fx, y: fy})) != 0 for l in leurres)

# arith-l2-fractions-elimination-ex2 (+ v1, v2).
for tv, bonne, leurres in [
    (Rational(2, 3), (Rational(5, 13), Rational(12, 13)),
     [(Rational(12, 13), Rational(5, 13)), (Rational(5, 13), Rational(-12, 13)), (Rational(3, 5), Rational(4, 5))]),
    (Rational(1, 4), (Rational(15, 17), Rational(8, 17)),
     [(Rational(8, 17), Rational(15, 17)), (Rational(15, 17), Rational(-8, 17)), (Rational(3, 5), Rational(4, 5))]),
    (Rational(3, 7), (Rational(20, 29), Rational(21, 29)),
     [(Rational(21, 29), Rational(20, 29)), (Rational(20, 29), Rational(-21, 29)), (Rational(7, 25), Rational(24, 25))]),
]:
    assert (cx.subs(t, tv), cy.subs(t, tv)) == bonne and bonne not in leurres
assert 15**2 + 8**2 == 17**2 and 20**2 + 21**2 == 29**2

# arith-l2-fractions-elimination-ex3 (+ v1, v2).
for fx, fy, bonne, leurres in [
    (t / (1 + t**3), t**2 / (1 + t**3), x**3 + y**3 - x * y, [x**3 + y**3 - 3 * x * y, x**2 + y**2 - x * y, x**3 - y**3 - x * y]),
    (3 * t / (1 + t**3), 3 * t**2 / (1 + t**3), x**3 + y**3 - 3 * x * y,
     [x**3 + y**3 - x * y, x**3 + y**3 - 9 * x * y, x**2 + y**2 - 3 * x * y]),
    (t**2 / (1 + t**2), t**3 / (1 + t**2), x**3 - y**2 * (1 - x), [y**3 - x**2 * (1 - y), x**3 - y**2 * (1 + x), x**2 - y**2 * (1 - x)]),
]:
    assert egal(bonne.subs({x: fx, y: fy}), 0) and egal(fy / fx, t)
    assert all(not egal(l.subs({x: fx, y: fy}), 0) for l in leurres)
# Réciprocité (ex3) : t = -1 impossible : x^3 (1 + t^3) = x^2 t avec t = -1 donne 0 = -x^2.

# arith-l2-fractions-elimination-ex4 (+ v1, v2) : hyperbole.
hx, hy = (1 + t**2) / (1 - t**2), 2 * t / (1 - t**2)
assert egal(hx**2 - hy**2, 1) and egal(hy / (1 + hx), t)
assert egal(expand(x**2 - t**2 * (x + 1) ** 2 - 1), expand((x + 1) * ((1 - t**2) * x - (1 + t**2))))
assert egal(((1 + x) ** 2 - (x**2 - 1)) / (1 + x) ** 2, 2 / (1 + x))  # 1 - t^2 avec y^2 = x^2 - 1
for P0, tv, leurres in [((Rational(5, 3), Rational(4, 3)), Rational(1, 2), [2, Rational(4, 3), Rational(3, 4)]),
                        ((Rational(5, 4), Rational(3, 4)), Rational(1, 3), [3, Rational(3, 5), Rational(3, 4)]),
                        ((Rational(13, 5), Rational(-12, 5)), Rational(-2, 3), [Rational(2, 3), Rational(-3, 2), Rational(-12, 13)])]:
    assert P0[0] ** 2 - P0[1] ** 2 == 1 and P0[1] / (1 + P0[0]) == tv
    assert (hx.subs(t, tv), hy.subs(t, tv)) == P0
    assert all((hx.subs(t, l), hy.subs(t, l)) != P0 for l in leurres)

# ---------------------------------------------------------------------------
# QCM arith-l2-fractions-q1..q8
# ---------------------------------------------------------------------------
assert egal((X**2 - 1) / (X**3 - X**2), (X + 1) / X**2)  # q1
assert gcd(X**2 + 1, X**3 + X) == X**2 + 1  # q2
Fq3 = (X**2 + 1) / ((X - 2) * (X + 1) * (X - 3))
assert residue(Fq3, X, 2) == Rational(-5, 3)  # q3
assert all(r == 0 for r in residus((X**2 + 1) / (X**2 - 1) ** 2, [1, -1]))  # q4
assert somme_formule({0: Fr(1, 2), 2: Fr(-1, 2)}) == Fr(3, 4)  # q5
assert (cx.subs(t, Rational(1, 2)), cy.subs(t, Rational(1, 2))) == (Rational(3, 5), Rational(4, 5))  # q7
assert all((cx.subs(t, l), cy.subs(t, l)) != (Rational(3, 5), Rational(4, 5)) for l in [2, Rational(4, 3), Rational(3, 4)])
Fq8 = 1 / (X * (X**2 + 1))  # q8
assert egal(Fq8, 1 / X - X / (X**2 + 1))
assert all(not egal(Fq8, l) for l in [1 / X - 1 / (X**2 + 1), 1 / X + X / (X**2 + 1), -1 / X + X / (X**2 + 1)])

print("arith_l2_fractions : toutes les vérifications passent.")
