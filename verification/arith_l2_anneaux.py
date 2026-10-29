"""Vérifications calculatoires du chapitre « arith-l2-anneaux » (arithmétique, L2) : cartes, exercices et QCM de
content/*/arithmetique/l2-arith-l2-anneaux.json et content/qcm/arith-l2-anneaux.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/arith_l2_anneaux.py (doit finir sans erreur).
"""

import datetime as dt
from fractions import Fraction as F
from itertools import product
from math import comb, gcd

from sympy import (GF, Poly, Rational, divisors, expand, factorint, gcd as pgcd, invert, isprime, rem, sqrt, symbols,
                   totient)

X = symbols("X")
x = symbols("x")


def lcm(a, b):
    return a * b // gcd(a, b)


def somme2(n):
    """n est-il somme de deux carrés d'entiers (force brute) ?"""
    return any(int((n - a * a) ** 0.5 + 0.5) ** 2 == n - a * a for a in range(int(n ** 0.5) + 1))


def red(expr, P):
    """Reste de expr modulo P dans Q[X]."""
    return expand(rem(expand(expr), P, X))


# ---------------------------------------------------------------------------
# arith-l2-anneaux-ideaux (carte)
# ---------------------------------------------------------------------------
# Exemple : (1+x)(2-3x) = 5 - x dans R[X]/(X^2+1) (et (1+i)(2-3i) = 5 - i).
assert red((1 + X) * (2 - 3 * X), X**2 + 1) == 5 - X
# Contre-exemple / pièges 1-2 : Z sous-anneau de Q, multiplication mal définie modulo Z.
a, a2 = F(1, 2), F(3, 2)
assert (a - a2).denominator == 1 and (a * a - a2 * a).denominator != 1
# Piège 3 / contre-exemple : J = 2Z u 3Z absorbant, relation non transitive.
J = lambda n: n % 2 == 0 or n % 3 == 0
assert all(J(k * n) for k in range(-20, 21) for n in range(-30, 31) if J(n))
assert J(0 - 2) and J(2 - 5) and not J(0 - 5)
# Affirmation 2 : 2 et 3 dans 2Z u 3Z, pas 5.
assert J(2) and J(3) and not J(5)
# Trou 3, leurres : I = {0} dans Z (0 dans I, contient 0 non inversible, Q corps avec I = {0}) : quotient non nul.

# arith-l2-anneaux-ideaux-ex1 (+ v1, v2) : produits dans Q[X]/(P).
for P, A, B, bon, leurres in [
    (X**2 - 3, 1 + 2 * X, 4 - X, -2 + 7 * X, [10 + 7 * X, 4 + 7 * X, -2 + 9 * X]),
    (X**2 + X + 1, 2 + X, 3 - X, 7 + 2 * X, [5 + 0 * X, 6 + X, 5 + 2 * X]),
    (X**2 - 2 * X + 2, 1 + X, 2 + 3 * X, -4 + 11 * X, [8 + 11 * X, -4 - X, 2 + 5 * X]),
]:
    assert red(A * B, P) == expand(bon) and all(red(A * B, P) != expand(l) for l in leurres)

# arith-l2-anneaux-ideaux-ex2 (+ v1, v2) : nombre d'idéaux de Z/nZ = nombre de diviseurs (force brute sur les
# sous-groupes additifs stables par multiplication, c'est-à-dire les parties engendrées).
def ideaux_zn(n):
    vus = set()
    for g in range(n):
        vus.add(frozenset((g * k) % n for k in range(n)))
    # dans Z/nZ, tout idéal est engendré par un élément (image d'un idéal principal de Z)
    return len(vus)


for n, bon, leurres, divs in [(12, 6, [12, 4, 2], [1, 2, 3, 4, 6, 12]),
                              (36, 9, [36, 12, 2], [1, 2, 3, 4, 6, 9, 12, 18, 36]),
                              (60, 12, [60, 16, 2], [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60])]:
    assert divisors(n) == divs and len(divs) == bon == ideaux_zn(n) and bon not in leurres
    assert int(totient(n)) in leurres or n in (12,)  # leurres « phi(n) » : 4 (n = 12), 12, 16
assert totient(12) == 4 and totient(36) == 12 and totient(60) == 16

# arith-l2-anneaux-ideaux-ex4 (+ v1, v2) : nilpotents de Z/nZ = multiples du radical.
for n, rad, card, leurres in [(72, 6, 12, [6, 24, 36]), (200, 10, 20, [10, 80, 100]), (108, 6, 18, [6, 36, 54])]:
    nil = [a for a in range(n) if any(pow(a, k, n) == 0 for k in range(1, 10))]
    assert nil == list(range(0, n, rad)) and len(nil) == card == n // rad and card not in leurres
    assert rad**3 % n == 0 and int(totient(n)) == leurres[1]

# ---------------------------------------------------------------------------
# arith-l2-anneaux-isomorphisme (carte)
# ---------------------------------------------------------------------------
# Contre-exemples : h(a,b) = 2a - b additive, h(1,1) = 1, (1,2) dans ker, (1,0)(1,2) = (1,0) hors ker.
h = lambda a, b: 2 * a - b
assert h(1, 1) == 1 and h(1, 2) == 0 and h(1 * 1, 0 * 2) != 0
# n -> n^2 : 2 n'est pas un carré.
assert all(k * k != 2 for k in range(3))
# Exemple : P(i) = 0 <=> X^2+1 | P (reste aX+b avec ai+b = 0).
assert rem(X**5 + X**3, X**2 + 1, X) == 0
# Affirmation 1 : 6 * f(1) = 6 != 0 dans Z.

# arith-l2-anneaux-isomorphisme-ex1 (+ v1, v2) : P(i) et reste modulo X^2+1.
from sympy import I as Ii

for P, val, R, leurres in [
    (X**5 + 2 * X**3 - X + 4, 4 - 2 * Ii, -2 * X + 4, [4 + 2 * Ii, 4, 4 - 4 * Ii]),
    (X**6 - 3 * X**4 + X + 1, -3 + Ii, X - 3, [-1 + Ii, -5 + Ii, -3 - Ii]),
    (X**7 + X**3 + 5 * X**2 - 2, -7 - 2 * Ii, -2 * X - 7, [3 - 2 * Ii, -7 + 2 * Ii, -7]),
]:
    assert expand(P.subs(X, Ii)) == val and rem(P, X**2 + 1, X) == R
    assert all(expand(l - val) != 0 for l in leurres)
    # les restes des leurres sont les polynômes aX + b associés aux valeurs fausses : distincts du bon reste.

# arith-l2-anneaux-isomorphisme-ex2 (+ v1, v2) : morphismes Z/aZ -> Z/bZ (force brute sur f(1) = 1 imposé).
def nb_morph(a, b):
    """Nombre d'applications k mod a -> k mod b bien définies (f(1) = 1 force f(k) = k)."""
    return 1 if all((k - kp) % b == 0 for k in range(a * b) for kp in range(a * b) if (k - kp) % a == 0) else 0


assert (nb_morph(12, 4), nb_morph(4, 12)) == (1, 0)
assert (nb_morph(10, 5), nb_morph(5, 10)) == (1, 0)
assert (nb_morph(6, 4), nb_morph(4, 6)) == (0, 0)
# remarque : k -> 3k (Z/4Z -> Z/12Z) bien définie, additive, envoie 1 sur 3.
assert all((3 * k - 3 * kp) % 12 == 0 for k in range(40) for kp in range(40) if (k - kp) % 4 == 0)

# arith-l2-anneaux-isomorphisme-ex3 (+ v1, v2) : image de Z dans Z/m x Z/n de cardinal ppcm(m, n).
for m, n, l, leurres in [(4, 6, 12, [24, 2, 6]), (6, 9, 18, [54, 3, 9]), (10, 15, 30, [150, 5, 15])]:
    img = {(k % m, k % n) for k in range(m * n)}
    assert len(img) == l == lcm(m, n) and l not in leurres and leurres[0] == m * n

# arith-l2-anneaux-isomorphisme-ex4 (+ v1, v2) : idempotents de R[X]/(P).
e = symbols("e")
for P, (r1, r2), antec, idem, leurres in [
    (X**2 - 1, (1, -1), lambda u, v: Rational(1, 2) * (u + v) + Rational(1, 2) * (u - v) * X,
     [(1 + X) / 2, (1 - X) / 2], [[1 + X, 1 - X], [X, -X], [(1 + X) / 2, (X - 1) / 2]]),
    (X**2 - 3 * X + 2, (1, 2), lambda u, v: (2 * u - v) + (v - u) * X,
     [2 - X, X - 1], [[X - 2, 1 - X], [X, 1 - X], [X / 2, 1 - X / 2]]),
    (X**2 - 4, (2, -2), lambda u, v: Rational(1, 2) * (u + v) + Rational(1, 4) * (u - v) * X,
     [(2 + X) / 4, (2 - X) / 4], [[(2 + X) / 2, (2 - X) / 2], [X / 2, -X / 2], [2 + X, 2 - X]]),
]:
    u, v = symbols("u v")
    Q = antec(u, v)
    assert expand(Q.subs(X, r1) - u) == 0 and expand(Q.subs(X, r2) - v) == 0
    assert expand(P - (X - r1) * (X - r2)) == 0
    for E in idem:
        assert red(E * E - E, P) == 0 and expand(E) not in (0, 1)
    for paire in leurres:  # au moins un élément de chaque paire leurre n'est pas idempotent
        assert any(red(E * E - E, P) != 0 for E in paire)
    assert red((X - r1) * (X - r2), P) == 0  # non intègre

# ---------------------------------------------------------------------------
# arith-l2-anneaux-integres (carte)
# ---------------------------------------------------------------------------
# Contre-exemple : Z/6Z, (0,1) ~ (0,2) ~ (3,1) mais (0,1) !~ (3,1) ; 2*3 = 0 ; pas de nilpotent non nul.
rel = lambda p, q: (p[0] * q[1] - p[1] * q[0]) % 6 == 0
assert rel((0, 1), (0, 2)) and rel((0, 2), (3, 1)) and not rel((0, 1), (3, 1))
assert (2 * 3) % 6 == 0 and [a for a in range(6) if any(pow(a, k, 6) == 0 for k in range(1, 7))] == [0]
# Trou 1, leurre « ac = bd » : (1,2) non équivalent à lui-même dans Z.
assert 1 * 1 != 2 * 2
# Trou 1, leurre « a + d = b + c » : (1,1) ~ (2,2) mais les sommes (2,1) et (4,2) ne sont pas équivalentes.
assert 1 + 2 == 1 + 2 and not (2 + 2 == 1 + 4)
# Exemple : 1/(1+sqrt2) = -1 + sqrt2.
assert expand((1 + sqrt(2)) * (-1 + sqrt(2))) == 1
# Affirmation 5 : sqrt2 hors de Q(i) (sinon partie imaginaire nulle et sqrt2 rationnel) ; sqrt2 irrationnel.
assert not sqrt(2).is_rational

# arith-l2-anneaux-integres-ex1 (+ v1, v2) : diviseurs de zéro de Z/nZ.
for n, bon, leurres in [(12, 7, [8, 4, 6]), (15, 6, [7, 8, 2]), (20, 11, [12, 8, 9])]:
    dz = [a for a in range(1, n) if any((a * b) % n == 0 for b in range(1, n))]
    inv = [a for a in range(1, n) if gcd(a, n) == 1]
    assert len(dz) == bon == n - 1 - totient(n) and bon not in leurres and not set(dz) & set(inv)
    assert len(dz) + len(inv) == n - 1

# arith-l2-anneaux-integres-ex3 (+ v1, v2) : inverses dans Q(sqrt2).
for z, inv, leurres in [
    (1 + 2 * sqrt(2), (-1 + 2 * sqrt(2)) / 7, [(1 - 2 * sqrt(2)) / 7, (1 - 2 * sqrt(2)) / 9, 1 + sqrt(2) / 2]),
    (3 + sqrt(2), (3 - sqrt(2)) / 7, [(3 - sqrt(2)) / 11, (-3 + sqrt(2)) / 7, Rational(1, 3) + sqrt(2)]),
    (2 - 3 * sqrt(2), -(2 + 3 * sqrt(2)) / 14, [(2 + 3 * sqrt(2)) / 14, (2 + 3 * sqrt(2)) / 22, Rational(1, 2) - sqrt(2) / 3]),
]:
    assert expand(z * inv) == 1 and all(expand(z * l) != 1 for l in leurres)

# arith-l2-anneaux-integres-ex4 (+ v1, v2) : inversibles de D = Z[1/10] : forme irréductible en 2 et 5 seulement.
def inversible_D(q):
    return all(p in (2, 5) for p in factorint(abs(q.numerator))) and all(p in (2, 5) for p in factorint(q.denominator))


for nombres, bon in [([F(64, 100), F(12, 100), F(15, 10), F(7, 10)], F(64, 100)),
                     ([F(6, 10), F(25, 10), F(14, 10), F(35, 10)], F(25, 10)),
                     ([F(21, 100), F(45, 100), F(11, 10), F(32, 1000)], F(32, 1000))]:
    assert [q for q in nombres if inversible_D(q)] == [bon]
assert 1 / F(64, 100) == F(15625, 10000) and 1 / F(25, 10) == F(4, 10) and 1 / F(32, 1000) == F(3125, 100)

# ---------------------------------------------------------------------------
# arith-l2-anneaux-principaux (carte)
# ---------------------------------------------------------------------------
# Exemple : (X+1)(X^2-X+1) = (X^3-2) + 3 ; X^3 - 2 sans racine rationnelle.
assert expand((X + 1) * (X**2 - X + 1) - (X**3 - 2)) == 3
assert all((Rational(p, 1)) ** 3 != 2 for p in (-2, -1, 1, 2))
assert Poly(X**3 - 2, X, domain="QQ").is_irreducible
# Trou 2, leurre « sans racine » : (X^2+1)^2 réductible ; leurres « unitaire », « deg <= 2 » : X^2.
assert not Poly((X**2 + 1) ** 2, X, domain="QQ").is_irreducible
# Trou 3 : 9 impair, 6 sans facteur carré, 4 puissance de premier : non premiers (Z/nZ non intègre).
for n in (9, 6, 4):
    assert not isprime(n) and any((a * b) % n == 0 for a in range(1, n) for b in range(1, n))
# Affirmations : X^2 - 2 irréductible sur Q, X^2 - 1 réductible.
assert Poly(X**2 - 2, X, domain="QQ").is_irreducible and not Poly(X**2 - 1, X, domain="QQ").is_irreducible

# arith-l2-anneaux-principaux-ex1 (+ v1, v2) : PGCD de X^a - 1 et X^b - 1.
for A, B, D, leurres in [(X**4 - 1, X**6 - 1, X**2 - 1, [X**12 - 1, X - 1, X**2 + 1]),
                         (X**6 - 1, X**9 - 1, X**3 - 1, [X**18 - 1, X - 1, X**3 + 1]),
                         (X**8 - 1, X**12 - 1, X**4 - 1, [X**24 - 1, X**2 - 1, X**4 + 1])]:
    assert expand(pgcd(A, B) - D) == 0 and all(expand(l - D) != 0 for l in leurres)
assert expand((X**6 - 1) - X**2 * (X**4 - 1) - (X**2 - 1)) == 0
assert expand((X**9 - 1) - X**3 * (X**6 - 1) - (X**3 - 1)) == 0
assert expand((X**12 - 1) - X**4 * (X**8 - 1) - (X**4 - 1)) == 0

# arith-l2-anneaux-principaux-ex2 (+ v1, v2) : polynômes nuls en a et b.
for a_, b_, P, leurres in [(2, 3, X**2 - 5 * X + 6, [X**2 + 5 * X + 6, X**2 - 5 * X - 6, 2 * X - 5]),
                           (1, -1, X**2 - 1, [X**2 + 1, X - 1, X**2 - 2 * X + 1]),
                           (0, 4, X**2 - 4 * X, [X**2 + 4 * X, X - 4, X**2 - 4])]:
    assert expand(P - (X - a_) * (X - b_)) == 0
    for l in leurres:  # chaque leurre ne s'annule pas en a et b, ou n'est pas de degré 2
        assert l.subs(X, a_) != 0 or l.subs(X, b_) != 0

# arith-l2-anneaux-principaux-ex3 (+ v1, v2) : inverses dans Q[X]/(P).
for P, A, inv, leurres in [
    (X**3 - 2, X + 1, (X**2 - X + 1) / 3, [(X**2 + X + 1) / 3, X**2 - X + 1, (X**2 - X + 1) / 2]),
    (X**3 - 2, X - 1, X**2 + X + 1, [X**2 - X + 1, (X**2 + X + 1) / 3, -X**2 - X - 1]),
    (X**2 + X + 1, X + 2, (1 - X) / 3, [(1 + X) / 3, 1 - X, (X - 1) / 3]),
]:
    assert expand(invert(A, P, X) - inv) == 0 and red(A * inv, P) == 1
    assert all(red(A * l, P) != 1 for l in leurres)
    assert Poly(P, X, domain="QQ").is_irreducible
assert expand((X - 1) * (X**2 + X + 1) - 1 - (X**3 - 2)) == 0
assert expand((X + 2) * (X - 1) + 3 - (X**2 + X + 1)) == 0

# arith-l2-anneaux-principaux-ex5 (+ v1, v2) : nombre de diviseurs unitaires.
for P, dom, bon, leurres in [(X**3 - X, "QQ", 8, [3, 4, 6]), (X**4 - 1, "RR", 8, [16, 4, 5]), (X**3 - X**2, "QQ", 6, [4, 3, 8])]:
    fl = Poly(P, X, domain="QQ").factor_list()[1]  # sur Q ; sur R, X^2+1 reste irréductible (pas de racine réelle)
    nb = 1
    for _, mult in fl:
        nb *= mult + 1
    assert nb == bon and bon not in leurres
assert len(Poly(X**4 - 1, X, extension=Ii).factor_list()[1]) == 4  # sur C : 16 idéaux

# ---------------------------------------------------------------------------
# arith-l2-anneaux-gauss (carte)
# ---------------------------------------------------------------------------
N = lambda a, b: a * a + b * b
# Inversibles : norme 1.
assert [(a, b) for a in range(-3, 4) for b in range(-3, 4) if N(a, b) == 1] == [(-1, 0), (0, -1), (0, 1), (1, 0)]
# Trou 2, leurre « norme au plus 2 » : 1+i de norme 2, non inversible (1/(1+i) = (1-i)/2).
assert expand(1 / (1 + Ii) - (1 - Ii) / 2) == 0
# Deux carrés : premiers impairs < 2000.
for p in range(3, 2000, 2):
    if isprime(p):
        assert somme2(p) == (p % 4 == 1)
# Wilson : x = ((p-1)/2)! vérifie x^2 = -1 mod p pour p = 1 mod 4.
from math import factorial

for p in (5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97):
    assert (factorial((p - 1) // 2) ** 2 + 1) % p == 0
# Cas général : n <= 3000.
for n in range(1, 3001):
    cond = all(e % 2 == 0 for q, e in factorint(n).items() if q % 4 == 3)
    assert somme2(n) == cond
# Trou 1, leurres : 5 = 1 + 4 et 5 non congru à 1 mod 8 ; « -1 non carré mod p » faux pour 5 (2^2 = -1).
assert somme2(5) and 5 % 8 != 1 and (2 * 2 + 1) % 5 == 0
# Trou 3, leurres : 9 = 9 + 0 somme (3 divise 9) ; 2 = 1 + 1 ; 5 somme avec exposant impair de 5.
assert somme2(9) and somme2(2) and somme2(5)
# Pièges / contre-exemple : 21 = 1 mod 4 non somme ; 3, 7, 11 non sommes.
assert 21 % 4 == 1 and not somme2(21) and not any(somme2(k) for k in (3, 7, 11))
# Exemple : produits, division de 7+5i par 2+i.
assert expand((1 + 2 * Ii) * (2 + 3 * Ii)) == -4 + 7 * Ii and expand((1 + 2 * Ii) * (2 - 3 * Ii)) == 8 + Ii
assert N(4, 7) == N(8, 1) == 65
assert expand((7 + 5 * Ii) - (4 + Ii) * (2 + Ii)) == -Ii
# Cours : dans Z[i sqrt5], 2 irréductible (pas d'élément de norme 2 = a^2 + 5b^2).
assert all(a * a + 5 * b * b != 2 for a in range(3) for b in range(2))
# Affirmations : 3 irréductible (3 non somme), 5 = (1+2i)(1-2i).
assert not somme2(3) and expand((1 + 2 * Ii) * (1 - 2 * Ii)) == 5


def div_gauss(z, w):
    q0 = expand(z / w)
    re, im = q0.as_real_imag()
    q = round(float(re)) + round(float(im)) * Ii
    return q, expand(z - q * w)


# arith-l2-anneaux-gauss-ex1 (+ v1, v2).
for z, w, q, r, leurres in [
    (11 + 3 * Ii, 3 - 2 * Ii, 2 + 2 * Ii, 1 + Ii, [(2 + 2 * Ii, 1 - Ii), (2 - 2 * Ii, 1 + Ii), (3 + 2 * Ii, -2 + Ii)]),
    (8 - 3 * Ii, 1 + 3 * Ii, -3 * Ii, -1, [(-3 * Ii, 1), (3 * Ii, -1), (-1 - 3 * Ii, -1)]),
    (13 + Ii, 4 + 3 * Ii, 2 - Ii, 2 - Ii, [(2 - Ii, 2 + Ii), (2 + Ii, 2 - Ii), (3 - Ii, -2 + Ii)]),
]:
    assert div_gauss(z, w) == (q, r)
    nr = expand(r * r.conjugate())
    nw = expand(w * w.conjugate())
    assert 2 * nr <= nw
    for ql, rl in leurres:  # les leurres ne sont même pas des divisions : z != q w + r
        assert expand(ql * w + rl - z) != 0

# arith-l2-anneaux-gauss-ex2 (+ v1, v2) : décompositions.
for p, (a_, b_), faux_ in [(29, (2, 5), [(3 + 2 * Ii) * (3 - 2 * Ii), (1 + Ii) * (1 - Ii) * (5 + 2 * Ii)]),
                           (37, (1, 6), [(2 + 5 * Ii) * (2 - 5 * Ii), (1 + Ii) * (1 - Ii) * (6 + Ii)]),
                           (53, (2, 7), [(1 + 7 * Ii) * (1 - 7 * Ii), (4 + 5 * Ii) * (4 - 5 * Ii)])]:
    assert isprime(p) and p % 4 == 1 and N(a_, b_) == p
    assert all(expand(f) != p for f in faux_)
assert expand((1 + Ii) * (1 - Ii) * (5 + 2 * Ii)) == 10 + 4 * Ii and expand((1 + Ii) * (1 - Ii) * (6 + Ii)) == 12 + 2 * Ii

# arith-l2-anneaux-gauss-ex3 (+ v1, v2).
for nombres, bon, ecr in [([245, 147, 63, 99], 245, (7, 14)), ([135, 84, 325, 300], 325, (10, 15)),
                          ([126, 192, 231, 98], 98, (7, 7))]:
    assert [n for n in nombres if somme2(n)] == [bon] and N(*ecr) == bon
assert N(1, 18) == N(6, 17) == 325
assert factorint(245) == {5: 1, 7: 2} and factorint(147) == {3: 1, 7: 2} and factorint(63) == {3: 2, 7: 1}
assert factorint(99) == {3: 2, 11: 1} and factorint(325) == {5: 2, 13: 1} and factorint(135) == {3: 3, 5: 1}
assert factorint(84) == {2: 2, 3: 1, 7: 1} and factorint(300) == {2: 2, 3: 1, 5: 2} and factorint(98) == {2: 1, 7: 2}
assert factorint(126) == {2: 1, 3: 2, 7: 1} and factorint(192) == {2: 6, 3: 1} and factorint(231) == {3: 1, 7: 1, 11: 1}

# arith-l2-anneaux-gauss-ex4 (+ v1, v2) : Z[i]/(w) = Z/pZ, i = c, classe de z.
for w, p, c, idiff, z, k, leurres in [
    (2 + Ii, 5, -2, 2 + Ii, 3 + 2 * Ii, 4, [1, 0, 2]),
    (3 + 2 * Ii, 13, 5, (3 + 2 * Ii) * (-1 + Ii), 1 + 2 * Ii, 11, [4, 3, 7]),
    (1 + 4 * Ii, 17, 4, (1 + 4 * Ii) * Ii, 2 + 3 * Ii, 14, [7, 5, 3]),
]:
    assert expand(w * w.conjugate()) == p and expand(Ii - c - idiff) == 0 and (c * c + 1) % p == 0
    re, im = z.as_real_imag()
    assert (re + c * im) % p == k and k not in leurres
    # vérification directe : (z - k)/w est un entier de Gauss
    qq = expand((z - k) / w)
    assert all(t.is_integer for t in qq.as_real_imag())

# ---------------------------------------------------------------------------
# arith-l2-anneaux-caracteristique (carte)
# ---------------------------------------------------------------------------
def car_produit(*mods):
    k = 1
    while any(k % m for m in mods):
        k += 1
    return k


assert car_produit(4, 6) == 12  # exemple
assert car_produit(6) == 6 and car_produit(4) == 4 and car_produit(2, 2) == 2  # contre-exemples
# Exemple : (X+1)^9 = X^9 + 1 dans F_3[X].
assert Poly((X + 1) ** 9, X, modulus=3) == Poly(X**9 + 1, X, modulus=3)
# Frobenius : p | C(p,k) pour 0 < k < p.
for p in (2, 3, 5, 7, 11, 13):
    assert all(comb(p, k) % p == 0 for k in range(1, p))
# Trou 3, leurres faux pour x = y = 1 (ou x = 1, y = 0) dans F_p.
for p in (2, 3, 5, 7):
    assert (1 + 1) ** p % p != (1 + 1 + 1) % p and (1 + 0) ** p % p != (1 * 0) % p and 1 % p != 0

# arith-l2-anneaux-caracteristique-ex1 (+ v1, v2).
for m, n, l, leurres in [(6, 10, 30, [60, 2, 10]), (3, 9, 9, [27, 3, 12]), (8, 12, 24, [96, 4, 12])]:
    assert car_produit(m, n) == l == lcm(m, n) and l not in leurres

# arith-l2-anneaux-caracteristique-ex2 (+ v1, v2).
for p, A, B, leurres in [(5, (X + 1) ** 25, X**25 + 1, [X**5 + 1, X**25 + X**5 + X + 1, X**25 - 1]),
                         (3, (X + 2) ** 9, X**9 + 2, [X**9 + 1, X**3 + 2, X**9 + 2 * X**3 + 2]),
                         (2, (X**2 + X + 1) ** 4, X**8 + X**4 + 1, [X**8 + 1, X**6 + X**3 + 1, X**8 + X**2 + 1])]:
    assert Poly(A, X, modulus=p) == Poly(B, X, modulus=p)
    assert all(Poly(A, X, modulus=p) != Poly(l, X, modulus=p) for l in leurres)
assert all(comb(25, k) % 5 == 0 for k in range(1, 25)) and pow(2, 9, 3) == 2

# arith-l2-anneaux-caracteristique-ex3 (+ v1, v2) : cardinaux de corps finis = puissances de premiers.
est_ppremier = lambda n: len(factorint(n)) == 1
for nombres, bon, P, p in [([6, 10, 9, 12], 9, X**2 + 1, 3), ([15, 18, 20, 25], 25, X**2 + 2, 5),
                           ([8, 14, 21, 24], 8, X**3 + X + 1, 2)]:
    assert [n for n in nombres if est_ppremier(n)] == [bon]
    assert Poly(P, X, modulus=p).is_irreducible and p ** Poly(P, X).degree() == bon
    assert all(P.subs(X, t) % p != 0 for t in range(p))

# arith-l2-anneaux-caracteristique-ex4 (+ v1, v2) : racines carrées dans F_4 et F_8.
for P, y, rac, leurres in [(X**2 + X + 1, X, X + 1, [X, 1, 0]), (X**2 + X + 1, X + 1, X, [X + 1, 1, 0]),
                           (X**3 + X + 1, X, X**2 + X, [X**2, X + 1, X**2 + 1])]:
    car2 = lambda E: Poly(rem(Poly(E, X, modulus=2).as_expr() ** 2, P, X, modulus=2), X, modulus=2)
    assert car2(rac) == Poly(y, X, modulus=2)
    assert all(car2(l) != Poly(y, X, modulus=2) for l in leurres)
    q = 2 ** Poly(P, X).degree()
    # a^q = a pour tous les éléments
    for coefs in product(range(2), repeat=Poly(P, X).degree()):
        E = sum(c * X**j for j, c in enumerate(coefs))
        assert Poly(rem(Poly(E, X, modulus=2).as_expr() ** q, P, X, modulus=2), X, modulus=2) == Poly(E, X, modulus=2)

# ---------------------------------------------------------------------------
# arith-l2-anneaux-applications (carte)
# ---------------------------------------------------------------------------
# Exemple RSA : p = 5, q = 11, e = 3, d = 27, 7 -> 13 -> 7 ; calendrier 2027 -> 2030.
assert pow(3, -1, 40) == 27 and pow(7, 3, 55) == 13 and pow(13, 27, 55) == 7
assert dt.date(2027, 1, 1).weekday() == 4 and dt.date(2030, 1, 1).weekday() == 1  # vendredi, mardi
assert (dt.date(2030, 1, 1) - dt.date(2027, 1, 1)).days == 1096 == 156 * 7 + 4
# RSA pour tout x (petits cas).
for p, q in [(3, 5), (5, 11), (7, 13), (11, 17)]:
    n, m = p * q, (p - 1) * (q - 1)
    for e_ in range(1, m):
        if gcd(e_, m) == 1:
            d_ = pow(e_, -1, m)
            assert all(pow(t, e_ * d_, n) == t for t in range(n))
# Pièges et contre-exemples : p = q = 3 ; p = 9, q = 5 ; ed = 1 mod n ; affine 2x+1.
assert 9 % 4 == 1 and pow(3, 9, 9) == 0
assert 33 % 32 == 1 and pow(2, 33, 45) == 17
assert 16 % 15 == 1 and pow(2, 16, 15) == 1
assert (2 * 0 + 1) % 26 == (2 * 13 + 1) % 26
# Trou 1 : « premiers entre eux », « deux entiers impairs distincts » (9, 5) ; « deux premiers » (3, 3) : ci-dessus.
# Trou 2 : ed = 1 mod (p+q) avec p = 3, q = 7, e = 3, d = 7 : 3^21 != 3 mod 21 ; ed = 1 mod (p-1) avec p = 3, q = 5, e = 3, d = 1.
assert (3 * 7) % 10 == 1 and pow(3, 21, 21) != 3
assert (3 * 1) % 2 == 1 and pow(2, 3, 15) != 2
# Trou 3 : x^e = x^d faux (p = 3, q = 11, e = 3, d = 7, x = 2) ; x^e = x faux (n = 15, e = 3, x = 2) ; x^ed = 1 faux (x = 0).
assert (3 * 7) % 20 == 1 and pow(2, 3, 33) != pow(2, 7, 33)
assert (3 * 3) % 8 == 1 and pow(2, 3, 15) != 2
# Affirmations : 4x+3 non injective ; 1er mars ; 5^81 = 5 mod 55.
assert (4 * 13 + 3) % 26 == 3
assert pow(5, 81, 55) == 5

# arith-l2-anneaux-applications-ex1 (+ v1, v2) : jours de la semaine (0 = lundi).
for d1, d2, j1, j2, leurres in [
    (dt.date(2026, 7, 14), dt.date(2030, 7, 14), 1, 6, [5, 0, 1]),
    (dt.date(2026, 12, 25), dt.date(2035, 12, 25), 4, 1, [0, 2, 4]),
    (dt.date(2027, 3, 1), dt.date(2031, 3, 1), 0, 5, [4, 6, 0]),
]:
    assert d1.weekday() == j1 and d2.weekday() == j2 and j2 not in leurres
assert (dt.date(2030, 7, 14) - dt.date(2026, 7, 14)).days == 1461 == 208 * 7 + 5
assert (dt.date(2035, 12, 25) - dt.date(2026, 12, 25)).days == 3287 == 469 * 7 + 4
assert (dt.date(2031, 3, 1) - dt.date(2027, 3, 1)).days == 1461

# arith-l2-anneaux-applications-ex2 (+ v1, v2) : déchiffrement affine modulo 26.
for a_, b_, (ap, bp), leurres in [(7, 3, (15, 7), [(15, 23), (7, 23), (19, 7)]),
                                  (5, 8, (21, 14), [(21, 8), (5, 18), (21, 18)]),
                                  (9, 2, (3, 20), [(3, 2), (3, 24), (9, 24)])]:
    assert pow(a_, -1, 26) == ap
    assert all((ap * (a_ * t + b_) + bp) % 26 == t for t in range(26))
    for al, bl in leurres:
        assert any((al * (a_ * t + b_) + bl) % 26 != t for t in range(26))
assert 3 * 26 - 11 * 7 == 1 and 26 - 5 * 5 == 1 and 3 * 9 - 26 == 1

# arith-l2-anneaux-applications-ex3 (+ v1, v2) : clés ISBN.
for digits, s, key, leurres in [("284123905", 170, 5, [6, 0, 10]), ("207115864", None, 4, [7, 0, 10]),
                                ("101928374", None, 2, [9, 0, 10])]:
    tot = sum((k + 1) * int(c) for k, c in enumerate(digits))
    assert (s is None or tot == s) and tot % 11 == key and key not in leurres
    full = [int(c) for c in digits] + [key]
    assert sum((k + 1) * a for k, a in enumerate(full)) % 11 == 0
# détection de toute erreur simple et de tout échange (force brute sur un numéro)
full = [2, 8, 4, 1, 2, 3, 9, 0, 5, 5]
for k in range(10):
    for v in range(11 if k == 9 else 10):
        if v != full[k]:
            g = full[:k] + [v] + full[k + 1:]
            assert sum((j + 1) * a for j, a in enumerate(g)) % 11 != 0
for j in range(10):
    for k in range(j + 1, 10):
        if full[j] != full[k]:
            g = full[:]
            g[j], g[k] = g[k], g[j]
            assert sum((i + 1) * a for i, a in enumerate(g)) % 11 != 0

# arith-l2-anneaux-applications-ex4 (+ v1, v2) : RSA.
for p, q, e_, y, d_, x_, leurres in [(7, 11, 7, 5, 43, 26, [(43, 47), (17, 26), (13, 5)]),
                                     (3, 17, 5, 6, 13, 27, [(13, 6), (7, 27), (19, 12)]),
                                     (5, 13, 5, 7, 29, 37, [(29, 7), (19, 37), (13, 22)])]:
    n, m = p * q, (p - 1) * (q - 1)
    assert pow(e_, -1, m) == d_ and pow(y, d_, n) == x_ and pow(x_, e_, n) == y
    for dl, xl in leurres:
        assert (e_ * dl) % m != 1 or pow(xl, e_, n) != y
assert 2 * 60 - 17 * 7 == 1 and 13 * 5 - 2 * 32 == 1 and 2 * 48 - 19 * 5 == 1
assert pow(5, 43, 7) == 5 and pow(5, 43, 11) == 4 and pow(6, 13, 17) == 10 and pow(7, 29, 5) == 2 and pow(7, 29, 13) == 11
assert pow(4, 7, 11) == 5 and 16384 % 11 == 5

# arith-l2-anneaux-applications-ex5 (+ v1, v2) : factoriser n connaissant (p-1)(q-1).
for n, m, s, delta, p, q, e_, d_, leurres in [
    (323, 288, 36, 4, 17, 19, 5, 173, [(13, 23, 173), (17, 19, 115), (11, 29, 173)]),
    (391, 352, 40, 36, 17, 23, 5, 141, [(19, 21, 141), (17, 23, 211), (13, 29, 141)]),
    (437, 396, 42, 16, 19, 23, 5, 317, [(17, 25, 317), (19, 23, 79), (13, 29, 317)]),
]:
    assert p * q == n and (p - 1) * (q - 1) == m and n - m + 1 == s == p + q and s * s - 4 * n == delta
    assert isprime(p) and isprime(q) and pow(e_, -1, m) == d_
    for pl, ql, dl in leurres:
        assert pl * ql != n or (e_ * dl) % m != 1
assert 2 * 288 - 115 * 5 == 1 and 141 * 5 - 2 * 352 == 1 and 396 - 79 * 5 == 1

# ---------------------------------------------------------------------------
# QCM arith-l2-anneaux
# ---------------------------------------------------------------------------
# q2 : R[X]/(X^2-1) non intègre ; q3 : 3*3 = 0 mod 9 ; q4 : seul X^2 - 2 irréductible, 15 non premier.
assert (3 * 3) % 9 == 0
assert Poly(X**2 - 2, X, domain="QQ").is_irreducible and not Poly(X**2 - 4, X, domain="QQ").is_irreducible and not isprime(15)
# q5 : 41 = 16 + 25 ; 43, 47, 51 non sommes de deux carrés.
assert 41 == 4**2 + 5**2 and not any(somme2(k) for k in (43, 47, 51))
# q7 : n = 33, e = 3 -> d = 7 ; leurres 11, 3, 13 faux.
assert pow(3, -1, 20) == 7 and all((3 * dl) % 20 != 1 for dl in (11, 3, 13))

# ---------------------------------------------------------------------------
# Ajouts de la relecture indépendante
# ---------------------------------------------------------------------------
# arith-l2-anneaux-applications (cours) : cycle solaire (28 ans) et cycle lunaire (19 ans), premiers entre eux,
# combinés en 532 ans ; piège p = q : p^(ed) = 0 modulo p^2 dès que ed >= 2 (et ed = 1 mod (p-1)^2).
assert gcd(28, 19) == 1 and lcm(28, 19) == 532
for p in (3, 5, 7):
    for ed in range(2, 200):
        if ed % ((p - 1) ** 2) == 1:
            assert pow(p, ed, p * p) == 0 != p
# arith-l2-anneaux-applications (cours) : p + q = n - (p-1)(q-1) + 1 et p, q racines de X^2 - (p+q)X + n.
for p, q in [(5, 11), (17, 19), (101, 103)]:
    n, m = p * q, (p - 1) * (q - 1)
    assert sorted(int(r) for r in Poly(X**2 - (n - m + 1) * X + n, X).all_roots()) == [p, q]
# arith-l2-anneaux-ideaux (trou 1, leurre « partie stable par addition et multiplication ») : N dans Z,
# 1 = 2 modulo N ? (2 - 1 = 1 dans N) mais 1 - 2 = -1 hors de N : la relation n'est même pas symétrique.
assert (2 - 1) >= 0 and (1 - 2) < 0
# arith-l2-anneaux-integres (trou 2, leurres) : Frac(Q) = Q est un corps (« jamais un corps » faux) ; Frac(Z) = Q
# n'est pas isomorphe à Z (2 inversible) ni fini.
assert F(1, 2) * 2 == 1
# arith-l2-anneaux-isomorphisme-ex4 (+ v1, v2) : antécédents de (u, v) donnés dans la correction.
u, v = symbols("u v")
for P0, (r1, r2), ant in [(X**2 - 1, (1, -1), (u + v) / 2 + (u - v) / 2 * X),
                          (X**2 - 3 * X + 2, (1, 2), (2 * u - v) + (v - u) * X),
                          (X**2 - 4, (2, -2), (u + v) / 2 + (u - v) / 4 * X)]:
    assert expand(ant.subs(X, r1) - u) == 0 and expand(ant.subs(X, r2) - v) == 0
# arith-l2-anneaux-principaux-ex3 v2 : X^2 + X + 1 = (X + 2)(X - 1) + 3.
assert expand((X + 2) * (X - 1) + 3 - (X**2 + X + 1)) == 0
# arith-l2-anneaux-gauss-ex1 : z/w = (27 + 31i)/13, parties réelle et imaginaire non demi-entières (arrondi sans ambiguïté).
for z, w in [(11 + 3 * Ii, 3 - 2 * Ii), (8 - 3 * Ii, 1 + 3 * Ii), (13 + Ii, 4 + 3 * Ii)]:
    for t in expand(z / w).as_real_imag():
        assert (2 * t).q != 1 or (2 * t).p % 2 == 0

print("arith_l2_anneaux : toutes les vérifications passent.")
