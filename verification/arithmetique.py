"""Vérifications calculatoires du thème « arithmetique » (cartes et exercices terminale et l1).

Chaque assert cite l'identifiant de la carte ou de l'exercice concerné.
Lancer : python3 verification/arithmetique.py (doit finir sans erreur).
"""

from math import gcd, factorial

import sympy
from sympy import Poly, Quaternion, Rational, expand, interpolate, symbols, totient

X = symbols("X")


def lcm(a, b):
    return a * b // gcd(a, b)


def solutions_systeme(a, m, b, n, borne):
    """Solutions dans [0, borne[ du système x ≡ a [m], x ≡ b [n]."""
    return [x for x in range(borne) if (x - a) % m == 0 and (x - b) % n == 0]


def est_division_euclidienne(a, b, q, r):
    return a == b * q + r and 0 <= r < b


# ---------------------------------------------------------------------------
# arith-t-division-euclidienne
# ---------------------------------------------------------------------------

# Exemple / affirmation 1 : division de -17 par 5 : q = -4, r = 3 (et non q = -3, r = -2).
assert est_division_euclidienne(-17, 5, -4, 3)
assert divmod(-17, 5) == (-4, 3)
assert -17 == 5 * (-3) - 2 and not (0 <= -2 < 5)
# Exemple : 17 = 5 × 3 + 2.
assert est_division_euclidienne(17, 5, 3, 2)

# Affirmation 2 : 17 = 5 × 2 + 7 mais 7 n'est pas le reste (7 ≥ 5).
assert 17 == 5 * 2 + 7 and not (7 < 5)

# Contre-exemple et leurre {{3}} « 0 ≤ r ≤ b » : a = b = 5 donne deux couples (1, 0) et (0, 5).
assert 5 == 5 * 1 + 0 and 5 == 5 * 0 + 5 and 0 <= 5 <= 5
# Leurre « |r| < b » : a = 1, b = 2 donne (0, 1) et (1, -1).
assert 1 == 2 * 0 + 1 and 1 == 2 * 1 + (-1) and abs(-1) < 2
# Leurre « 0 < r < b » : a = b = 5, aucun reste r avec 0 < r < 5 (5 - r non multiple de 5).
assert all((5 - r) % 5 != 0 for r in range(1, 5))
# Contre-exemples b = 0 et b < 0 : aucun r avec 0 ≤ r < b.
assert list(range(0, 0)) == [] and list(range(0, -3)) == []

# Affirmation 3 (vraie) : si r ≠ 0, le reste de -a par b est b - r.
for b in range(1, 15):
    for a in range(-60, 61):
        r = a % b
        if r != 0:
            assert (-a) % b == b - r, (a, b)

# Démonstration : le quotient est la partie entière de a/b.
for b in range(1, 12):
    for a in range(-50, 51):
        assert a // b == sympy.floor(Rational(a, b))

# arith-t-division-euclidienne-ex1 : -2026 = 7 × (-290) + 4.
assert 2026 == 7 * 289 + 3
assert est_division_euclidienne(-2026, 7, -290, 4)
assert divmod(-2026, 7) == (-290, 4)
for q, r in [(-289, -3), (-289, 3), (-290, -4)]:  # leurres : tous faux
    assert not est_division_euclidienne(-2026, 7, q, r), (q, r)

# arith-t-division-euclidienne-ex2 : n² ≡ 0 ou 1 [4], une somme de deux carrés n'est jamais ≡ 3 [4].
assert {(n * n) % 4 for n in range(-100, 101)} == {0, 1}
assert all(((2 * k) ** 2 - 4 * k * k) == 0 and ((2 * k + 1) ** 2 - 4 * (k * k + k) - 1) == 0 for k in range(-50, 50))
assert all((x * x + y * y) % 4 != 3 for x in range(-60, 61) for y in range(-60, 61))
assert 2027 == 4 * 506 + 3
assert not any(x * x + y * y == 2027 for x in range(0, 46) for y in range(0, 46))

# arith-t-division-euclidienne-ex3 : a ≡ 7 [12].
restes_par_8 = set()
for k in range(-50, 51):
    a = 12 * k + 7
    assert (-a) % 12 == 5  # (a)
    assert -a == 12 * (-k - 1) + 5
    assert a % 4 == 3  # (b)
    assert (a * a) % 12 == 1  # (d)
    assert a * a == 12 * (12 * k * k + 14 * k + 4) + 1
    restes_par_8.add(a % 8)
assert restes_par_8 == {3, 7}  # (c) : non déterminé
assert 7 % 8 == 7 and 19 % 8 == 3 and 19 % 12 == 7
for leurre in [-7, 7, -5]:  # leurres de la réponse finale (bonne : 5)
    assert leurre != 5

# ---------------------------------------------------------------------------
# arith-t-bezout
# ---------------------------------------------------------------------------

# Exemple : algorithme d'Euclide pour 41 et 29, et 41 × (-12) + 29 × 17 = 1.
assert 41 == 29 * 1 + 12 and 29 == 12 * 2 + 5 and 12 == 5 * 2 + 2 and 5 == 2 * 2 + 1
assert 1 == 5 - 2 * 2
assert 1 == 5 * 5 - 2 * 12
assert 1 == 5 * 29 - 12 * 12
assert 1 == 17 * 29 - 12 * 41
assert 41 * (-12) + 29 * 17 == 1 and gcd(41, 29) == 1
# Exemple (identité de Bézout) : 12 ∧ 18 = 6 = 12 × (-1) + 18 × 1.
assert gcd(12, 18) == 6 == 12 * (-1) + 18 * 1

# Contre-exemple / affirmation 3 : 2 ∧ 3 = 1 mais 2u + 3v = 1 n'a pas de solution dans N².
assert gcd(2, 3) == 1
assert not any(2 * u + 3 * v == 1 for u in range(0, 50) for v in range(0, 50))
# Affirmation 2 : non-unicité, (-1, 1) et (2, -1).
assert 2 * (-1) + 3 * 1 == 1 and 2 * 2 + 3 * (-1) == 1
# Toutes les solutions de 2u + 3v = 1 : (-1 + 3k, 1 - 2k).
assert all(2 * (-1 + 3 * k) + 3 * (1 - 2 * k) == 1 for k in range(-20, 21))
# Affirmation 1 : 4 × 1 + 6 × 1 = 10 mais 4 ∧ 6 = 2 ; et 4 × (-1) + 6 × 1 = 2.
assert 4 * 1 + 6 * 1 == 10 and gcd(4, 6) == 2 and 4 * (-1) + 6 * 1 == 2
# Affirmation 4 (vraie) : a ∧ b divise toute combinaison au + bv.
for a in range(-12, 13):
    for b in range(-12, 13):
        if (a, b) != (0, 0):
            d = gcd(a, b)
            assert all((a * u + b * v) % d == 0 for u in range(-5, 6) for v in range(-5, 6))
# Démonstration : le plus petit élément > 0 de aZ + bZ est a ∧ b.
for a in range(-15, 16):
    for b in range(-15, 16):
        if (a, b) != (0, 0):
            combis = {a * u + b * v for u in range(-20, 21) for v in range(-20, 21)}
            assert min(c for c in combis if c > 0) == gcd(a, b)

# Leurres de trous : (u, v) ∈ Q² : 2 × 1/2 + 4 × 0 = 1 alors que 2 ∧ 4 = 2.
assert 2 * Rational(1, 2) + 4 * 0 == 1 and gcd(2, 4) == 2
# Leurre « au + bv premier » : 2 × 1 + 4 × 0 = 2 premier alors que 2 ∧ 4 = 2.
assert sympy.isprime(2 * 1 + 4 * 0) and gcd(2, 4) != 1
# Leurre « au + bv = a ∧ b » : toujours réalisable (identité de Bézout), même pour 2 et 4.
assert 2 * 1 + 4 * 0 == gcd(2, 4)
# Leurre « (u, v) ∈ N² » déjà vu avec 2 et 3 ; leurre « unique » vu avec (-1, 1) et (2, -1).

# arith-t-bezout-ex1 : 2(3n + 5) - 3(2n + 3) = 1.
n = symbols("n")
assert expand(2 * (3 * n + 5) - 3 * (2 * n + 3)) == 1
assert all(gcd(3 * k + 5, 2 * k + 3) == 1 for k in range(-300, 301))

# arith-t-bezout-ex2 : couple solution de 41u + 29v = 1 et leurres faux.
assert 41 * (-12) + 29 * 17 == 1
for u, v in [(12, -17), (17, -12), (-11, 16)]:
    assert 41 * u + 29 * v != 1, (u, v)
assert 41 * 12 + 29 * (-17) == -1  # erreur de signe typique

# arith-t-bezout-ex3 : (n + 4) ∧ (2n + 3) vaut 5 si n ≡ 1 [5], 1 sinon.
assert expand(2 * (n + 4) - (2 * n + 3)) == 5
for k in range(-500, 501):
    attendu = 5 if k % 5 == 1 else 1
    assert gcd(k + 4, 2 * k + 3) == attendu, k
assert 2026 % 5 == 1 and gcd(2030, 4055) == 5
assert 2027 % 5 == 2 and gcd(2031, 4057) == 1
assert gcd(406, 811) == 1 and 811 == 2 * 406 - 1

# ---------------------------------------------------------------------------
# arith-t-gauss
# ---------------------------------------------------------------------------

# Contre-exemple : 6 | 4 × 3 sans diviser 4 ni 3.
assert (4 * 3) % 6 == 0 and 4 % 6 != 0 and 3 % 6 != 0
assert gcd(6, 4) == 2 and gcd(6, 3) == 3
# Affirmation 1 (fausse) : 4 | 12 et 6 | 12 mais 24 ∤ 12.
assert 12 % 4 == 0 and 12 % 6 == 0 and 12 % 24 != 0
# Affirmation 2 (vraie) et théorème de Gauss vérifiés par force brute.
for a in range(1, 30):
    for b in range(-30, 31):
        for c in range(-30, 31):
            if gcd(a, b) == 1 and (b * c) % a == 0:
                assert c % a == 0, (a, b, c)
for a in range(1, 20):
    for b in range(1, 20):
        if gcd(a, b) == 1:
            for c in range(-80, 81):
                if c % a == 0 and c % b == 0:
                    assert c % (a * b) == 0
# Affirmation 4 (vraie, lemme d'Euclide).
for p in sympy.primerange(2, 30):
    for b in range(-30, 31):
        for c in range(-30, 31):
            if (b * c) % p == 0:
                assert b % p == 0 or c % p == 0

# Hypothèses pièges : « a premier » seul (a = 2, b = 2, c = 1) ; « b ∧ c = 1 » (6, 2, 3) ;
# « a ∧ c = 1 » (2, 2, 1) ; « a ne divise pas b » (6, 4, 3).
assert (2 * 1) % 2 == 0 and 1 % 2 != 0
assert (2 * 3) % 6 == 0 and gcd(2, 3) == 1 and 3 % 6 != 0
assert gcd(2, 1) == 1 and (2 * 1) % 2 == 0 and 1 % 2 != 0
assert 4 % 6 != 0 and (4 * 3) % 6 == 0 and 3 % 6 != 0
# Leurres du trou {{1}} : « a | b + c » (2, 1, 1), « a | b²c² » (4, 1, 2), « a | abc » (2, 1, 1).
assert (1 + 1) % 2 == 0 and gcd(2, 1) == 1 and 1 % 2 != 0
assert (1 * 1 * 2 * 2) % 4 == 0 and gcd(4, 1) == 1 and 2 % 4 != 0
assert (2 * 1 * 1) % 2 == 0 and 1 % 2 != 0
# Leurres du trou {{3}} : « a | b » (2, 1, 2), « ab | c » (2, 3, 2), « a | b et a | c » (2, 1, 2).
assert (1 * 2) % 2 == 0 and gcd(2, 1) == 1 and 1 % 2 != 0
assert (3 * 2) % 2 == 0 and gcd(2, 3) == 1 and 2 % 6 != 0

# Exemple : 7 | 12n ⇒ 7 | n.
assert gcd(7, 12) == 1
assert all(k % 7 == 0 for k in range(-200, 201) if (12 * k) % 7 == 0)

# arith-t-gauss-ex1 : 9 | 14n ⇒ 9 | n.
assert gcd(9, 14) == 1
assert all(k % 9 == 0 for k in range(-500, 501) if (14 * k) % 9 == 0)

# arith-t-gauss-ex2 : 6 | n(n + 1)(2n + 1).
for k in range(-300, 301):
    assert (k * (k + 1)) % 2 == 0
    assert (k * (k + 1) * (2 * k + 1)) % 3 == 0
    assert (k * (k + 1) * (2 * k + 1)) % 6 == 0
assert all((2 * k + 1) % 3 == 0 for k in range(-30, 30) if k % 3 == 1)

# arith-t-gauss-ex3 : solutions de 41x + 29y = 3 : (-7 + 29k, 10 - 41k).
assert 41 * (-7) + 29 * 10 == 3
assert 3 * (-12) == -36 and 3 * 17 == 51 and -36 + 29 == -7 and 51 - 41 == 10
boite = range(-400, 401)
sol = {(x, y) for x in boite for y in boite if 41 * x + 29 * y == 3}
param = {(-7 + 29 * k, 10 - 41 * k) for k in range(-40, 41)}
param_boite = {(x, y) for (x, y) in param if x in boite and y in boite}
assert sol == param_boite and len(sol) > 10
# Leurres : (-7 + 41k, 10 - 29k) et (-7 + 29k, 10 + 41k) contiennent des non-solutions ;
# (-7 + 87k, 10 - 123k) oublie la solution (22, -31).
assert 41 * (-7 + 41) + 29 * (10 - 29) != 3
assert 41 * (-7 + 29) + 29 * (10 + 41) != 3
assert (22, -31) in sol and (22 + 7) % 87 != 0

# arith-t-gauss-ex4 : 12 | n² ⇔ 6 | n.
assert {k for k in range(-600, 601) if (k * k) % 12 == 0} == {k for k in range(-600, 601) if k % 6 == 0}
assert (6 * 6) % 12 == 0 and 6 % 12 != 0  # leurre 12Z
assert (3 * 3) % 12 != 0  # leurre 3Z
assert (2 * 2) % 12 != 0  # leurre 2Z ∪ 3Z
assert all((36 * k * k) == 12 * (3 * k * k) for k in range(-10, 11))

# ---------------------------------------------------------------------------
# arith-t-fermat
# ---------------------------------------------------------------------------

# Théorème : a^p ≡ a [p] et, si p ∤ a, a^(p-1) ≡ 1 [p].
for p in sympy.primerange(2, 60):
    for a in range(-80, 81):
        assert (pow(a, p) - a) % p == 0
        if a % p != 0:
            assert pow(a, p - 1, p) == 1
    # Démonstration : p divise C(p, k) pour 1 ≤ k ≤ p - 1, et k C(p, k) = p C(p-1, k-1).
    for k in range(1, p):
        assert sympy.binomial(p, k) % p == 0
        assert k * sympy.binomial(p, k) == p * sympy.binomial(p - 1, k - 1)

# Exemple : 3^2026 ≡ 4 [7].
assert 2026 == 6 * 337 + 4 and pow(3, 6, 7) == 1 and 81 % 7 == 4
assert pow(3, 2026, 7) == 4

# Contre-exemple n = 4 : 3³ ≡ 3 [4] (et non 1), 2⁴ ≡ 0 [4] (et non 2).
assert gcd(3, 4) == 1 and pow(3, 3, 4) == 3
assert pow(2, 4, 4) == 0
# Contre-exemple a = p : p^(p-1) ≡ 0 [p].
assert all(pow(p, p - 1, p) == 0 for p in sympy.primerange(2, 30))
# Réciproque fausse : 341 = 11 × 31 et 2^340 ≡ 1 [341].
assert 341 == 11 * 31 and not sympy.isprime(341)
assert pow(2, 340, 341) == 1
assert 2 ** 10 == 1024 == 3 * 341 + 1 and 340 == 10 * 34
# Affirmation 4 (fausse) : 2³ ≡ 1 [7] et 6 ∤ 3.
assert pow(2, 3, 7) == 1 and 3 % 6 != 0

# Hypothèses pièges et leurres de trous.
assert pow(2, 8, 9) != 1 and pow(2, 9, 9) != 2  # p = 9, impair ≥ 3, non premier
assert pow(3, 340, 341) != 1  # « p tel que 2^(p-1) ≡ 1 [p] » : p = 341, a = 3
assert gcd(3, 341) == 1
assert pow(3, 2, 3) == 0  # « a quelconque », « a non nul », « a premier » : a = p = 3, a^(p-1) ≡ 0
assert pow(2, 3, 3) != 1  # « a^p ≡ 1 » : a = 2, p = 3 (8 ≡ 2)
assert pow(2, 2, 3) != 2  # « a^(p-1) ≡ a » : a = 2, p = 3 (4 ≡ 1)
assert pow(1, 2, 3) != (-1) % 3  # « a^(p-1) ≡ -1 » : a = 1, p = 3

# arith-t-fermat-ex1 : 5^2026 ≡ 5 [11].
assert pow(5, 10, 11) == 1 and 2026 == 10 * 202 + 6
assert pow(5, 6, 11) == 5 and (25 % 11) == 3 and (9 % 11) == 9 and 27 % 11 == 5
assert pow(5, 2026, 11) == 5
for leurre in [1, 3, 9]:
    assert leurre != pow(5, 2026, 11)

# arith-t-fermat-ex2 : inverse de 6 modulo 13 et équation 6x ≡ 5 [13].
assert pow(6, 11, 13) == 11 and (6 * 11) % 13 == 1
assert (11 * 5) % 13 == 3 and (6 * 3) % 13 == 5
assert [x for x in range(13) if (6 * x) % 13 == 5] == [3]
for leurre in [11, 5, 10]:
    assert (6 * leurre) % 13 != 5

# arith-t-fermat-ex3 : 42 | n^7 - n.
assert 42 == 2 * 3 * 7
for k in range(-300, 301):
    assert (k ** 7 - k) % 42 == 0
    assert (k ** 2 - k) % 2 == 0 and (k ** 3 - k) % 3 == 0
    assert (k ** 7 - k ** 3) % 3 == 0  # n^7 = n^3 · n^3 · n ≡ n · n · n = n^3 [3]

# arith-t-fermat-ex4 : 3^340 modulo 11, 31 et 341.
assert pow(3, 10, 11) == 1 and pow(3, 340, 11) == 1
assert 27 % 31 == 27 and (27 - 31) == -4 and 243 % 31 == 26 and (26 - 31) == -5
assert 340 == 30 * 11 + 10 and pow(3, 10, 31) == 25 and (-5) ** 2 == 25
assert pow(3, 340, 31) == 25
assert pow(3, 340, 341) != 1
for leurre in [1, 5, 9]:
    assert leurre != 25

# ---------------------------------------------------------------------------
# arith-l1-chinois
# ---------------------------------------------------------------------------

# Théorème vérifié par force brute : pour m ∧ n = 1, une seule solution dans [0, mn[.
for m in range(1, 16):
    for nn in range(1, 16):
        if gcd(m, nn) == 1:
            for a in range(m):
                for b in range(nn):
                    assert len(solutions_systeme(a, m, b, nn, m * nn)) == 1

# Démonstration : x0 = a n v + b m u avec mu + nv = 1.
for m in range(1, 13):
    for nn in range(1, 13):
        if gcd(m, nn) == 1:
            u, v, g = sympy.gcdex(m, nn)
            assert g == 1 and m * u + nn * v == 1
            for a in range(-5, 6):
                for b in range(-5, 6):
                    x0 = a * nn * v + b * m * u
                    assert (x0 - a) % m == 0 and (x0 - b) % nn == 0

# Exemple : x ≡ 2 [5], x ≡ 3 [7] ; 5 × 3 - 7 × 2 = 1, x0 = 17.
assert 5 * 3 + 7 * (-2) == 1
assert 2 * 7 * (-2) + 3 * 5 * 3 == 17
assert solutions_systeme(2, 5, 3, 7, 35) == [17]

# Contre-exemples sans coprimalité.
# Z/4Z : la classe de 1 est d'ordre additif 4 ; dans Z/2Z × Z/2Z, 2y = 0 pour tout y.
assert min(k for k in range(1, 5) if (k * 1) % 4 == 0) == 4
assert all(((2 * y1) % 2, (2 * y2) % 2) == (0, 0) for y1 in range(2) for y2 in range(2))
# x ≡ 0 [2], x ≡ 1 [4] : aucune solution.
assert solutions_systeme(0, 2, 1, 4, 8) == []
# x ≡ 1 [4], x ≡ 1 [6] : 1 et 13 solutions distinctes modulo 24.
assert solutions_systeme(1, 4, 1, 6, 24) == [1, 13]
# Affirmation 4 (fausse) : x ≡ 3 [6], x ≡ 5 [10] admet 15 pour solution.
assert 15 % 6 == 3 and 15 % 10 == 5
assert solutions_systeme(3, 6, 5, 10, 60) == [15, 45]
# Affirmation 3 (vraie) : (1, 1) engendre Z/mZ × Z/nZ quand m ∧ n = 1.
for m in range(1, 13):
    for nn in range(1, 13):
        ordre = min(k for k in range(1, m * nn + 1) if k % m == 0 and k % nn == 0)
        assert ordre == lcm(m, nn)
        if gcd(m, nn) == 1:
            assert ordre == m * nn
        else:
            assert ordre < m * nn
# Isomorphisme : l'application x mod mn ↦ (x mod m, x mod n) est bijective si m ∧ n = 1.
for m, nn in [(4, 9), (5, 7), (8, 15), (1, 6)]:
    images = {(x % m, x % nn) for x in range(m * nn)}
    assert len(images) == m * nn
# et ni injective ni surjective sinon (cas 2, 2 et 4, 6).
for m, nn in [(2, 2), (4, 6)]:
    images = {(x % m, x % nn) for x in range(m * nn)}
    assert len(images) < m * nn

# Leurres de trous et hypothèses pièges.
assert solutions_systeme(0, 2, 1, 2, 4) == []  # m = n = 2 premiers mais égaux
assert solutions_systeme(1, 2, 2, 4, 8) == [] and gcd(1, 2) == 1  # « a ∧ b = 1 » ne suffit pas
assert solutions_systeme(0, 3, 1, 9, 27) == []  # « m et n impairs »
assert solutions_systeme(0, 2, 1, 4, 8) == []  # « distincts », « l'un divise l'autre »
x0 = solutions_systeme(1, 2, 2, 3, 6)[0]
assert x0 == 5
assert ((x0 + 5) - 1) % 2 != 0  # leurre x0 + (m+n)Z : x0 + 5 = 10 n'est pas ≡ 1 [2]
assert ((x0 + 6) - 1) % 2 == 0 and ((x0 + 6) - 2) % 3 == 0  # leurre {x0} : x0 + mn est solution
assert ((x0 + 2) - 2) % 3 != 0  # leurre x0 + mZ

# arith-l1-chinois-ex1 : plus petit naturel avec x ≡ 1 [3], x ≡ 2 [4] : 10.
assert solutions_systeme(1, 3, 2, 4, 12) == [10]
for leurre in [2, 7, 14]:
    assert not ((leurre - 1) % 3 == 0 and (leurre - 2) % 4 == 0)

# arith-l1-chinois-ex2 : x ≡ 4 [9], x ≡ 6 [11] : 94 + 99Z.
assert 9 * 5 + 11 * (-4) == 1
assert 4 * 11 * (-4) + 6 * 9 * 5 == 94
assert solutions_systeme(4, 9, 6, 11, 99) == [94]
for leurre in [94 + 20, 22, 40]:
    assert not ((leurre - 4) % 9 == 0 and (leurre - 6) % 11 == 0)

# arith-l1-chinois-ex3 : φ(77) = 60.
assert sum(1 for x in range(77) if gcd(x, 77) == 1) == 60 == totient(77)
assert (7 - 1) * (11 - 1) == 60
for leurre in [76, 70, 66]:
    assert leurre != 60
for m in range(1, 20):
    for nn in range(1, 20):
        if gcd(m, nn) == 1:
            assert totient(m * nn) == totient(m) * totient(nn)

# arith-l1-chinois-ex4 : cas général, d = m ∧ n.
for m in range(1, 13):
    for nn in range(1, 13):
        d = gcd(m, nn)
        L = lcm(m, nn)
        for a in range(m):
            for b in range(nn):
                s = solutions_systeme(a, m, b, nn, L)
                if (a - b) % d == 0:
                    assert len(s) == 1
                    # construction x = a - m u k avec mu + nv = d et a - b = dk
                    u, v, g = sympy.gcdex(m, nn)
                    assert g == d
                    k = (a - b) // d
                    x = a - m * u * k
                    assert (x - a) % m == 0 and (x - b) % nn == 0
                else:
                    assert s == []
assert gcd(6, 10) == 2 and lcm(6, 10) == 30
assert solutions_systeme(3, 6, 5, 10, 30) == [15]
assert solutions_systeme(1, 6, 4, 10, 60) == []
assert (45 - 3) % 6 == 0 and (45 - 5) % 10 == 0 and (45 - 15) % 60 != 0  # leurre 15 + 60Z
assert (5 - 3) % 6 != 0  # leurre 5 + 30Z
# Calcul de la correction : 6 × 2 + 10 × (-1) = 2, a - b = -2 = 2 × (-1), x = 3 - 6 × 2 × (-1) = 15.
assert 6 * 2 + 10 * (-1) == 2 and 3 - 5 == 2 * (-1) and 3 - 6 * 2 * (-1) == 15

# ---------------------------------------------------------------------------
# arith-l1-racines-polynome
# ---------------------------------------------------------------------------

# Contre-exemple : X² - 1 a 4 racines dans Z/8Z.
assert [x for x in range(8) if (x * x - 1) % 8 == 0] == [1, 3, 5, 7]
# Contre-exemple : X² - X a 4 racines dans Z/6Z.
assert [x for x in range(6) if (x * x - x) % 6 == 0] == [0, 1, 3, 4]
# Contre-exemple : quaternions, (bi + cj + dk)² = -(b² + c² + d²).
for q in [Quaternion(0, 1, 0, 0), Quaternion(0, 0, 1, 0), Quaternion(0, 0, 0, 1),
          Quaternion(0, -1, 0, 0), Quaternion(0, Rational(3, 5), Rational(4, 5), 0),
          Quaternion(0, Rational(2, 3), Rational(2, 3), Rational(1, 3))]:
    assert q * q == Quaternion(-1, 0, 0, 0)
b_, c_, d_ = symbols("b c d", real=True)
carre = Quaternion(0, b_, c_, d_) * Quaternion(0, b_, c_, d_)
assert expand(carre.a + b_ ** 2 + c_ ** 2 + d_ ** 2) == 0 and carre.b == 0 and carre.c == 0 and carre.d == 0
# Exemple : X³ - 3X + 2 = (X - 1)²(X + 2).
assert expand((X - 1) ** 2 * (X + 2)) == X ** 3 - 3 * X + 2
# Affirmation 3 (fausse) : X² + 1 n'a pas de racine réelle.
assert sympy.real_roots(X ** 2 + 1) == []
# Affirmation « X^p - X » : non nul mais s'annule partout sur Z/pZ.
for p in sympy.primerange(2, 30):
    assert all((pow(x, p) - x) % p == 0 for x in range(p))
# Théorème vérifié sur Z/pZ : nombre de racines ≤ degré pour des polynômes aléatoires.
import random
random.seed(0)
for _ in range(2000):
    p = random.choice([2, 3, 5, 7, 11])
    deg = random.randint(0, 6)
    coeffs = [random.randrange(p) for _ in range(deg)] + [random.randrange(1, p)]
    nb = sum(1 for x in range(p) if sum(c * x ** i for i, c in enumerate(coeffs)) % p == 0)
    assert nb <= deg
# Leurres de trous : X (degré 1, une racine) contredit « strictement inférieure à n » ;
# X² + 1 sur R contredit « égale à n » et « supérieure ou égale à n ».
assert sympy.roots(X, X) == {0: 1}

# arith-l1-racines-polynome-ex1 : P de degré ≤ 4 avec P(k) = k² + 1 pour k = -2..2 : P = X² + 1.
pts = [(k, k * k + 1) for k in range(-2, 3)]
P = interpolate(pts, X)
assert expand(P - (X ** 2 + 1)) == 0
assert P.subs(X, 3) == 10
for leurre in [9, 28]:
    assert leurre != 10

# arith-l1-racines-polynome-ex2 : x² ≡ 1 [p] ⇔ x ≡ ±1 [p] ; 4 racines dans Z/8Z.
for p in sympy.primerange(3, 200):
    assert [x for x in range(p) if (x * x - 1) % p == 0] == [1, p - 1]
assert [x for x in range(2) if (x * x - 1) % 2 == 0] == [1]
assert len([x for x in range(8) if (x * x - 1) % 8 == 0]) == 4
assert (2 * 4) % 8 == 0  # 2 et 4 sont des diviseurs de zéro dans Z/8Z
for leurre in [2, 3, 8]:
    assert leurre != 4

# arith-l1-racines-polynome-ex3 : racines de X² - 1 dans Z/15Z.
assert [x for x in range(15) if (x * x - 1) % 15 == 0] == [1, 4, 11, 14]
for x, (r3, r5) in [(1, (1, 1)), (4, (1, 4)), (11, (2, 1)), (14, (2, 4))]:
    assert x % 3 == r3 and x % 5 == r5
assert (6 * 6 - 1) % 15 != 0
assert (3 * 5) % 15 == 0  # Z/15Z n'est pas intègre
assert ((4 - 1) * (4 + 1)) % 15 == 0 and (4 - 1) % 15 != 0 and (4 + 1) % 15 != 0

# arith-l1-racines-polynome-ex4 : X^(p-1) - 1 = ∏ (X - a) dans (Z/pZ)[X], puis Wilson.
for p in [2, 3, 5, 7, 11, 13]:
    gauche = Poly(X ** (p - 1) - 1, X, modulus=p)
    droite = Poly(1, X, modulus=p)
    for a in range(1, p):
        droite = droite * Poly(X - a, X, modulus=p)
    assert gauche == droite, p
for p in sympy.primerange(2, 80):
    assert (factorial(p - 1) + 1) % p == 0
assert factorial(12) % 13 == 12
for leurre in [1, 0, 6]:
    assert leurre != factorial(12) % 13
# Réciproque (remarque de la correction) : pour n composé > 4, (n-1)! ≡ 0 [n] ; pour n = 4, 3! ≡ 2.
assert factorial(3) % 4 == 2
assert all(factorial(k - 1) % k == 0 for k in range(5, 80) if not sympy.isprime(k))

# ---------------------------------------------------------------------------
# Compléments (valeurs citées dans les corrections et affirmations)
# ---------------------------------------------------------------------------

# arith-l1-chinois (affirmation 2) : 17 et 52 sont solutions de x ≡ 2 [5], x ≡ 3 [7].
assert all((x - 2) % 5 == 0 and (x - 3) % 7 == 0 for x in [17, 52])
# arith-l1-chinois-ex1 : restes modulo 4 de 1, 4, 7, 10.
assert [x % 4 for x in [1, 4, 7, 10]] == [1, 0, 3, 2]
assert 10 == 3 * 3 + 1 == 4 * 2 + 2
# arith-l1-chinois-ex2 : 114 = 9 × 12 + 6 et 94 = 9 × 10 + 4 = 11 × 8 + 6.
assert 114 == 9 * 12 + 6 and 94 == 9 * 10 + 4 == 11 * 8 + 6
# arith-l1-chinois-ex3 : contrôle par inclusion-exclusion.
assert len(range(0, 77, 7)) == 11 and len(range(0, 77, 11)) == 7 and 77 - 11 - 7 + 1 == 60
# arith-l1-racines-polynome-ex2 : carrés des classes paires modulo 8.
assert [(x * x) % 8 for x in [0, 2, 4, 6]] == [0, 4, 0, 4]
assert 9 == 8 + 1 and 25 == 3 * 8 + 1 and 49 == 6 * 8 + 1
# arith-l1-racines-polynome-ex3 : 4² = 16 = 15 + 1 et (4 - 1)(4 + 1) = 15.
assert 4 ** 2 == 15 + 1 and (4 - 1) * (4 + 1) == 15
# arith-t-fermat-ex4 : 243 = 7 × 31 + 26.
assert 3 ** 5 == 243 == 7 * 31 + 26
# arith-t-fermat-ex2 : 66 = 5 × 13 + 1, 55 ≡ 3 [13], 18 ≡ 5 [13].
assert 66 == 5 * 13 + 1 and 55 % 13 == 3 and 18 % 13 == 5
# arith-t-bezout-ex2 : -492 + 493 = 1.
assert 41 * (-12) == -492 and 29 * 17 == 493
# arith-t-gauss-ex3 : -287 + 290 = 3.
assert 41 * (-7) == -287 and 29 * 10 == 290
# arith-t-bezout-ex3 : n = 0 donne A = 4, B = 3, d = 1 ; 2026 = 5 × 405 + 1.
assert gcd(4, 3) == 1 and 2026 == 5 * 405 + 1 and 2030 == 5 * 406 and 4055 == 5 * 811


# ===========================================================================
# Enrichissement (lot enr-arithmetique) : variantes (variante k = variantes[k-1]) et exercices ajoutés.
# ===========================================================================
import json
import pathlib

from sympy import isprime, factorint, fibonacci, mod_inverse, Symbol

RACINE = pathlib.Path(__file__).resolve().parent.parent
EXOS = {}
for niv in ("terminale", "l1"):
    for e in json.loads((RACINE / f"content/exercices/arithmetique/{niv}.json").read_text())["exercices"]:
        EXOS[e["id"]] = e


def rf(ident, k=0):
    """reponse_finale de l'exercice (k = 0) ou de sa variante k."""
    e = EXOS[ident]
    return e["reponse_finale"] if k == 0 else e["variantes"][k - 1]["reponse_finale"]


def bonne(ident, k, attendu, leurres):
    """La proposition juste est `attendu`, les leurres sont exactement `leurres` (tous distincts)."""
    r = rf(ident, k)
    assert r["propositions"][r["bonne"]] == attendu, (ident, k, r["propositions"][r["bonne"]])
    autres = [x for i, x in enumerate(r["propositions"]) if i != r["bonne"]]
    assert autres == leurres and len(set(r["propositions"])) == 4, (ident, k)


def nb_variantes(ident):
    return len(EXOS[ident].get("variantes", []))


# Règle du lot : 3 à 6 exercices par carte, une difficulté 1-2 et une 4-5 ; 2 variantes pour chaque exercice
# à valeurs numériques (seuls arith-t-gauss-ex2 et arith-t-fermat-ex3, purs exercices de démonstration, n'en ont pas).
for carte in ["arith-t-division-euclidienne", "arith-t-bezout", "arith-t-gauss", "arith-t-fermat",
              "arith-l1-chinois", "arith-l1-racines-polynome"]:
    ds = [e["difficulte"] for e in EXOS.values() if e["carte"] == carte]
    assert 3 <= len(ds) <= 6 and min(ds) <= 2 and max(ds) >= 4, carte
for ident, e in EXOS.items():
    attendu = 0 if ident in ("arith-t-gauss-ex2", "arith-t-fermat-ex3") else 2
    assert nb_variantes(ident) == attendu, ident


# --- arith-t-division-euclidienne-ex1, variantes 1 et 2 : -1789 par 9, -3005 par 11.
assert 1789 == 9 * 198 + 7 and -1789 == 9 * (-198) - 7
assert est_division_euclidienne(-1789, 9, -199, 2) and divmod(-1789, 9) == (-199, 2)
for q, r in [(-198, -7), (-198, 7), (-199, -2)]:  # leurres v1
    assert not est_division_euclidienne(-1789, 9, q, r)
assert 3005 == 11 * 273 + 2 and -3005 == 11 * (-273) - 2
assert est_division_euclidienne(-3005, 11, -274, 9) and divmod(-3005, 11) == (-274, 9)
for q, r in [(-273, -2), (-273, 2), (-274, -9)]:  # leurres v2
    assert not est_division_euclidienne(-3005, 11, q, r)
bonne("arith-t-division-euclidienne-ex1", 1, "$q = -199$, $r = 2$", ["$q = -198$, $r = -7$", "$q = -198$, $r = 7$", "$q = -199$, $r = -2$"])
bonne("arith-t-division-euclidienne-ex1", 2, "$q = -274$, $r = 9$", ["$q = -273$, $r = -2$", "$q = -273$, $r = 2$", "$q = -274$, $r = -9$"])

# --- arith-t-division-euclidienne-ex2, v1 : 10^6 + 3 ≡ 3 [4] ; v2 : 2^10 - 1 = 1023 ≡ 3 [4] (donc pas sommes de deux carrés).
assert 10 ** 6 == 4 * 250000 and 10 ** 6 + 3 == 4 * 250000 + 3
assert 2 ** 10 - 1 == 1023 == 4 * 255 + 3
for N in (10 ** 6 + 3, 1023):
    assert not any(sympy.sqrt(N - x * x).is_integer for x in range(0, sympy.floor(sympy.sqrt(N)) + 1))

# --- arith-t-division-euclidienne-ex3, v1 : a ≡ 11 [18] ; v2 : a ≡ 8 [15].
k = symbols("k", integer=True)
for m, r, d, e, opp, rd, a2, carre in [(18, 11, 6, 4, 7, 5, 29, 13), (15, 8, 5, 6, 7, 3, 23, 4)]:
    for kk in range(-40, 41):
        a = m * kk + r
        assert (-a) % m == opp and a % d == rd and (a * a) % m == carre
    assert {(m * kk + r) % e for kk in range(-40, 41)} != {r % e}  # reste par e non déterminé
    assert r % e != a2 % e and a2 == m + r
assert expand(-(18 * k + 11) - (18 * (-k - 1) + 7)) == 0 and expand((18 * k + 11) - (6 * (3 * k + 1) + 5)) == 0
assert expand((18 * k + 11) ** 2 - (18 * (18 * k ** 2 + 22 * k + 6) + 13)) == 0
assert 11 == 4 * 2 + 3 and 29 == 4 * 7 + 1
assert expand(-(15 * k + 8) - (15 * (-k - 1) + 7)) == 0 and expand((15 * k + 8) - (5 * (3 * k + 1) + 3)) == 0
assert expand((15 * k + 8) ** 2 - (15 * (15 * k ** 2 + 16 * k + 4) + 4)) == 0
assert 8 == 6 * 1 + 2 and 23 == 6 * 3 + 5
bonne("arith-t-division-euclidienne-ex3", 1, "$7$", ["$-11$", "$11$", "$-7$"])
bonne("arith-t-division-euclidienne-ex3", 2, "$7$", ["$-8$", "$8$", "$-7$"])

# --- arith-t-division-euclidienne-ex4 (ajouté) : reste de 2^a - 1 par 2^b - 1 = 2^r - 1.
x, q = symbols("x q")
for qq in range(1, 8):
    assert expand(x ** qq - 1 - (x - 1) * sum(x ** i for i in range(qq))) == 0
for b in range(1, 13):
    for a in range(0, 80):
        assert (2 ** a - 1) % (2 ** b - 1) == 2 ** (a % b) - 1, (a, b)
        qq, r = divmod(a, b)
        m_ = sum(2 ** (b * i) for i in range(qq))
        assert 2 ** a - 1 == (2 ** b - 1) * 2 ** r * m_ + (2 ** r - 1)
# Exercice : 2^100 - 1 par 127 ; v1 : 2^75 - 1 par 255 ; v2 : 2^64 - 1 par 31.
for k_, (a, b, r, reste) in enumerate([(100, 7, 2, 3), (75, 8, 3, 7), (64, 5, 4, 15)]):
    assert a % b == r and (2 ** a - 1) % (2 ** b - 1) == reste == 2 ** r - 1
    assert pow(2, a, 2 ** b - 1) == reste + 1  # leurre « reste de 2^a »
    assert reste not in (reste + 1, r, 0)  # leurres faux
    bonne("arith-t-division-euclidienne-ex4", k_, f"${reste}$", [f"${reste + 1}$", f"${r}$", "$0$"])
assert 100 == 7 * 14 + 2 and 75 == 8 * 9 + 3 and 64 == 5 * 12 + 4
# Question 4 : l'affirmation sur 2^a + 1 est fausse (a, b) = (3, 2), (5, 3), (4, 3).
for a, b, val, quo, reste in [(3, 2, 9, 1, 4), (5, 3, 33, 3, 6), (4, 3, 17, 1, 8)]:
    r = a % b
    assert 2 ** a + 1 == val == (2 ** b + 1) * quo + reste and 0 <= reste < 2 ** b + 1 and reste != 2 ** r + 1
    assert 2 ** a + 1 == 2 ** r * (2 ** (b * (a // b)) - 1) + 2 ** r + 1
    assert (2 ** (b * (a // b)) - 1) % (2 ** b + 1) != 0

# --- arith-t-bezout-ex1, v1 : (-2)(5n + 7) + 5(2n + 3) = 1 ; v2 : 3(4n + 3) - 4(3n + 2) = 1.
n = symbols("n", integer=True)
assert expand(-2 * (5 * n + 7) + 5 * (2 * n + 3)) == 1
assert expand(3 * (4 * n + 3) - 4 * (3 * n + 2)) == 1
assert all(gcd(5 * t + 7, 2 * t + 3) == 1 and gcd(4 * t + 3, 3 * t + 2) == 1 for t in range(-200, 201))

# --- arith-t-bezout-ex2, v1 : 53u + 37v = 1 ; v2 : 47u + 30v = 1.
assert 53 == 37 * 1 + 16 and 37 == 16 * 2 + 5 and 16 == 5 * 3 + 1
assert 7 * 16 - 3 * 37 == 1 and 7 * 53 - 10 * 37 == 1 and 53 * 7 == 371 and 37 * 10 == 370
for u, v in [(-7, 10), (-10, 7), (6, -9)]:  # leurres v1
    assert 53 * u + 37 * v != 1
assert 47 == 30 + 17 and 30 == 17 + 13 and 17 == 13 + 4 and 13 == 4 * 3 + 1
assert 4 * 13 - 3 * 17 == 1 and 4 * 30 - 7 * 17 == 1 and 11 * 30 - 7 * 47 == 1 and 47 * 7 == 329 and 30 * 11 == 330
for u, v in [(7, -11), (11, -7), (-6, 10)]:  # leurres v2
    assert 47 * u + 30 * v != 1
bonne("arith-t-bezout-ex2", 1, "$(7, -10)$", ["$(-7, 10)$", "$(-10, 7)$", "$(6, -9)$"])
bonne("arith-t-bezout-ex2", 2, "$(-7, 11)$", ["$(7, -11)$", "$(11, -7)$", "$(-6, 10)$"])

# --- arith-t-bezout-ex3, v1 : (n + 3) ∧ (3n + 2) = 7 ssi n ≡ 4 [7] ; v2 : (n + 2) ∧ (4n + 5) = 3 ssi n ≡ 1 [3].
assert expand(3 * (n + 3) - (3 * n + 2)) == 7 and expand(4 * (n + 2) - (4 * n + 5)) == 3
for t in range(-300, 301):
    assert 3 * t + 2 != 0 and 4 * t + 5 != 0
    assert gcd(t + 3, 3 * t + 2) == (7 if t % 7 == 4 else 1)
    assert gcd(t + 2, 4 * t + 5) == (3 if t % 3 == 1 else 1)
assert gcd(3, 2) == 1 and gcd(2, 5) == 1  # n = 0
assert 2026 == 7 * 289 + 3 and 2029 == 7 * 289 + 6 and gcd(2029, 3 * 2026 + 2) == 1
assert 2026 == 3 * 675 + 1 and 2028 == 3 * 676 and 4 * 2026 + 5 == 8109 == 3 * 2703 and gcd(2028, 8109) == 3
bonne("arith-t-bezout-ex3", 1, "$1$", ["$7$", "$14$", "$49$"])
bonne("arith-t-bezout-ex3", 2, "$3$", ["$1$", "$9$", "$6$"])

# --- arith-t-bezout-ex4 (ajouté) : Cassini F(n+1)F(n-1) - F(n)^2 = (-1)^n et coefficients de Bézout.
F = [fibonacci(i) for i in range(40)]
assert F[0] == 0 and F[1] == 1 and all(F[i + 2] == F[i + 1] + F[i] for i in range(38))
for i in range(1, 39):
    assert F[i + 1] * F[i - 1] - F[i] ** 2 == (-1) ** i
    u, v = (-1) ** i * F[i - 1], (-1) ** (i + 1) * F[i]
    assert u * F[i + 1] + v * F[i] == 1 and gcd(F[i], F[i + 1]) == 1
    assert i + 2 >= 40 or (u * F[i + 2] + (v - u) * F[i] == 1 and gcd(F[i], F[i + 2]) == 1)
assert (F[9], F[10], F[11], F[12], F[8]) == (34, 55, 89, 144, 21)
for k_, (A, B, sol, leurres) in enumerate([
    (89, 55, (34, -55), [(-34, 55), (21, -34), (55, -89)]),
    (144, 89, (-55, 89), [(55, -89), (-34, 55), (89, -144)]),
    (55, 34, (-21, 34), [(21, -34), (-13, 21), (34, -55)]),
]):
    assert A * sol[0] + B * sol[1] == 1
    for u, v in leurres:
        assert A * u + B * v != 1, (A, B, u, v)
    bonne("arith-t-bezout-ex4", k_, f"$({sol[0]}, {sol[1]})$", [f"$({u}, {v})$" for u, v in leurres])
assert 89 * 34 == 3026 and 55 * 55 == 3025 and 89 * 21 == 1869 and 55 * 34 == 1870
assert 144 * 55 == 7920 and 89 * 89 == 7921 and 144 * 34 == 4896 and 89 * 55 == 4895
assert 55 * 21 == 1155 and 34 * 34 == 1156 and 55 * 13 == 715 and 34 * 21 == 714

# --- arith-t-gauss-ex1, v1 : 15 ∧ 22 = 1 ; v2 : 8 ∧ 35 = 1.
assert gcd(15, 22) == 1 and gcd(8, 35) == 1 and 22 == 2 * 11 and 35 == 5 * 7
assert all((22 * t) % 15 != 0 or t % 15 == 0 for t in range(-300, 301))
assert all((35 * t) % 8 != 0 or t % 8 == 0 for t in range(-300, 301))

# --- arith-t-gauss-ex3, v1 : 53x + 37y = 2 ; v2 : 47x + 30y = 4.
for a, b, c, (u, v), (x0, y0), oubli in [(53, 37, 2, (7, -10), (14, -20), (51, -73)), (47, 30, 4, (-7, 11), (2, -3), (32, -50))]:
    assert a * u + b * v == 1 and a * x0 + b * y0 == c
    sols = {(xx, yy) for xx in range(-400, 401) for yy in [(c - a * xx) // b] if a * xx + b * yy == c}
    attendu = {(x0 + b * t, y0 - a * t) for t in range(-30, 31)}
    assert {s for s in sols if -400 <= s[0] <= 400} == {s for s in attendu if -400 <= s[0] <= 400}
    assert a * oubli[0] + b * oubli[1] == c and (oubli[0] - x0) % (c * b) != 0  # solution manquée par le pas c(b, -a)
    assert a * (x0 + a) + b * (y0 - b) != c and a * (x0 + b) + b * (y0 + a) != c  # leurres 2 et 3 contiennent des non-solutions
assert 53 * 14 == 742 and 37 * 20 == 740 and (2 * 7, 2 * (-10)) == (14, -20)
assert (4 * (-7), 4 * 11) == (-28, 44) and (-28 + 30, 44 - 47) == (2, -3) and 47 * 2 == 94 and 30 * 3 == 90
bonne("arith-t-gauss-ex3", 1, r"$\{(14 + 37k,\ -20 - 53k) : k \in \mathbb{Z}\}$",
      [r"$\{(14 + 53k,\ -20 - 37k) : k \in \mathbb{Z}\}$", r"$\{(14 + 37k,\ -20 + 53k) : k \in \mathbb{Z}\}$",
       r"$\{(14 + 74k,\ -20 - 106k) : k \in \mathbb{Z}\}$"])
bonne("arith-t-gauss-ex3", 2, r"$\{(2 + 30k,\ -3 - 47k) : k \in \mathbb{Z}\}$",
      [r"$\{(2 + 47k,\ -3 - 30k) : k \in \mathbb{Z}\}$", r"$\{(2 + 30k,\ -3 + 47k) : k \in \mathbb{Z}\}$",
       r"$\{(2 + 120k,\ -3 - 188k) : k \in \mathbb{Z}\}$"])

# --- arith-t-gauss-ex4, v1 : 20 | n² ⇔ 10 | n ; v2 : 45 | n² ⇔ 15 | n ; leurres faux.
for N, m in [(20, 10), (45, 15)]:
    assert {t for t in range(-500, 501) if (t * t) % N == 0} == {t for t in range(-500, 501) if t % m == 0}
    assert (m * m) % N == 0 and m % N != 0  # pas N Z
assert 10 ** 2 == 100 == 5 * 20 and 25 % 20 != 0 and 4 % 20 != 0 and 15 ** 2 == 225 == 5 * 45 and 25 % 45 != 0 and 9 % 45 != 0
bonne("arith-t-gauss-ex4", 1, r"$10\mathbb{Z}$", [r"$20\mathbb{Z}$", r"$5\mathbb{Z}$", r"$2\mathbb{Z} \cup 5\mathbb{Z}$"])
bonne("arith-t-gauss-ex4", 2, r"$15\mathbb{Z}$", [r"$45\mathbb{Z}$", r"$5\mathbb{Z}$", r"$3\mathbb{Z} \cup 5\mathbb{Z}$"])

# --- arith-t-fermat-ex1, v1 : 3^2027 ≡ 9 [13] ; v2 : 7^2026 ≡ 2 [17].
assert isprime(13) and isprime(17) and 2027 == 12 * 168 + 11 and 2026 == 16 * 126 + 10
assert pow(3, 2027, 13) == 9 == pow(3, 11, 13) and 27 == 2 * 13 + 1
assert pow(7, 2026, 17) == 2 == pow(7, 10, 17) and 49 == 3 * 17 - 2 and pow(7, 4, 17) == 4 and pow(7, 8, 17) == 16
bonne("arith-t-fermat-ex1", 1, "$9$", ["$1$", "$3$", "$12$"])
bonne("arith-t-fermat-ex1", 2, "$2$", ["$1$", "$7$", "$15$"])

# --- arith-t-fermat-ex2, v1 : 5^{-1} ≡ 7 [17], 5x ≡ 3 ⇔ x ≡ 4 ; v2 : 4^{-1} ≡ 3 [11], 4x ≡ 9 ⇔ x ≡ 5.
assert mod_inverse(5, 17) == 7 == pow(5, 15, 17) and 35 == 2 * 17 + 1 and 21 % 17 == 4 and 20 % 17 == 3
assert [t for t in range(17) if (5 * t) % 17 == 3] == [4] and all((5 * t) % 17 != 3 for t in (7, 3, 15))
assert mod_inverse(4, 11) == 3 == pow(4, 9, 11) and 12 == 11 + 1 and 27 % 11 == 5 and 20 % 11 == 9
assert [t for t in range(11) if (4 * t) % 11 == 9] == [5] and all((4 * t) % 11 != 9 for t in (3, 9, 8))
bonne("arith-t-fermat-ex2", 1, "$4$", ["$7$", "$3$", "$15$"])
bonne("arith-t-fermat-ex2", 2, "$5$", ["$3$", "$9$", "$8$"])

# --- arith-t-fermat-ex4, v1 : 2047 = 23 × 89 ; v2 : 1387 = 19 × 73 (pseudo-premiers de Fermat en base 2).
assert factorint(2047) == {23: 1, 89: 1} and 2 ** 11 == 2048 == 2047 + 1 and 2046 == 11 * 186 and pow(2, 2046, 2047) == 1
assert 2046 == 88 * 23 + 22 and pow(3, 2046, 89) == 34 == pow(3, 22, 89)
assert 81 - 89 == -8 and 64 * 64 == 4096 == 46 * 89 + 2 and (2 * (-8) * 9) == -144 and -144 + 2 * 89 == 34
assert 2046 == 22 * 93 and pow(3, 2046, 23) == 1 and pow(3, 2046, 2047) != 1
assert factorint(1387) == {19: 1, 73: 1} and 2 ** 9 == 512 == 7 * 73 + 1 and 1386 == 9 * 154 == 18 * 77
assert pow(2, 1386, 73) == 1 and pow(2, 1386, 19) == 1 and pow(2, 1386, 1387) == 1
assert 1386 == 72 * 19 + 18 and pow(3, 1386, 73) == 72 == pow(3, 18, 73)
assert 81 % 73 == 8 and 64 - 73 == -9 and 81 % 73 == 8 and (8 * 9) % 73 == 72 and pow(3, 1386, 1387) != 1
bonne("arith-t-fermat-ex4", 1, "$34$", ["$1$", "$9$", "$55$"])
bonne("arith-t-fermat-ex4", 2, "$72$", ["$1$", "$9$", "$8$"])

# --- arith-l1-chinois-ex1, v1 : x ≡ 1 [3], x ≡ 2 [5] → 7 + 15Z ; v2 : x ≡ 3 [4], x ≡ 4 [5] → 19 + 20Z.
assert solutions_systeme(1, 3, 2, 5, 60) == [7, 22, 37, 52] and [t % 3 for t in (2, 7, 12)] == [2, 1, 0]
assert solutions_systeme(3, 4, 4, 5, 80) == [19, 39, 59, 79] and [t % 4 for t in (4, 9, 14, 19)] == [0, 1, 2, 3]
assert 7 == 3 * 2 + 1 == 5 * 1 + 2 and 19 == 4 * 4 + 3 == 5 * 3 + 4
bonne("arith-l1-chinois-ex1", 1, "$7$", ["$2$", "$12$", "$22$"])  # 22 est solution mais pas la plus petite
bonne("arith-l1-chinois-ex1", 2, "$19$", ["$3$", "$4$", "$39$"])

# --- arith-l1-chinois-ex2, v1 : x ≡ 3 [7], x ≡ 5 [8] → 45 + 56Z ; v2 : x ≡ 2 [5], x ≡ 9 [13] → 22 + 65Z.
assert 7 * (-1) + 8 * 1 == 1 and 3 * 8 * 1 + 5 * 7 * (-1) == -11 == 7 * (-2) + 3 == 8 * (-2) + 5 and -11 % 56 == 45
assert solutions_systeme(3, 7, 5, 8, 56) == [45] and 60 == 7 * 8 + 4
assert 5 * (-5) + 13 * 2 == 1 and 2 * 13 * 2 + 9 * 5 * (-5) == -173 == 5 * (-35) + 2 == 13 * (-14) + 9 and -173 % 65 == 22
assert solutions_systeme(2, 5, 9, 13, 65) == [22] and 40 == 5 * 8
for c, mod, (a, m, b, n) in [(45, 15, (3, 7, 5, 8)), (24, 56, (3, 7, 5, 8)), (35, 56, (3, 7, 5, 8)),
                             (22, 18, (2, 5, 9, 13)), (52, 65, (2, 5, 9, 13)), (9, 65, (2, 5, 9, 13))]:
    classe = {c + mod * t for t in range(-5, 6)}
    assert any((y - a) % m or (y - b) % n for y in classe)  # chaque leurre contient une non-solution
bonne("arith-l1-chinois-ex2", 1, r"$45 + 56\mathbb{Z}$", [r"$45 + 15\mathbb{Z}$", r"$24 + 56\mathbb{Z}$", r"$35 + 56\mathbb{Z}$"])
bonne("arith-l1-chinois-ex2", 2, r"$22 + 65\mathbb{Z}$", [r"$22 + 18\mathbb{Z}$", r"$52 + 65\mathbb{Z}$", r"$9 + 65\mathbb{Z}$"])

# --- arith-l1-chinois-ex3, v1 : φ(91) = 72 ; v2 : φ(85) = 64.
assert totient(91) == 72 == 6 * 12 and 91 - 13 - 7 + 1 == 72 and totient(85) == 64 == 4 * 16 and 85 - 17 - 5 + 1 == 64
assert len(range(0, 91, 7)) == 13 and len(range(0, 91, 13)) == 7 and len(range(0, 85, 5)) == 17 and len(range(0, 85, 17)) == 5
assert all(t not in (72, 64) for t in (90, 84, 78, 80, 68))
bonne("arith-l1-chinois-ex3", 1, "$72$", ["$90$", "$84$", "$78$"])
bonne("arith-l1-chinois-ex3", 2, "$64$", ["$84$", "$80$", "$68$"])

# --- arith-l1-chinois-ex4, v1 : x ≡ 5 [12], x ≡ 11 [18] → 29 + 36Z ; x ≡ 1 [12], x ≡ 4 [18] sans solution.
assert gcd(12, 18) == 6 and lcm(12, 18) == 36 and 12 * (-1) + 18 * 1 == 6 and 5 - 11 == 6 * (-1)
assert 5 - 12 * (-1) * (-1) == -7 == 12 * (-1) + 5 == 18 * (-1) + 11 and -7 % 36 == 29
assert solutions_systeme(5, 12, 11, 18, 216) == [29, 65, 101, 137, 173, 209] and solutions_systeme(1, 12, 4, 18, 216) == []
assert (1 - 4) % 6 != 0 and 65 % 12 == 5 and 65 % 18 == 11 and (65 - 29) % 216 != 0 and 11 % 12 != 5
# v2 : x ≡ 7 [15], x ≡ 2 [25] → 52 + 75Z ; x ≡ 1 [15], x ≡ 3 [25] sans solution.
assert gcd(15, 25) == 5 and lcm(15, 25) == 75 and 15 * 2 + 25 * (-1) == 5 and 7 - 2 == 5 * 1
assert 7 - 15 * 2 * 1 == -23 == 15 * (-2) + 7 == 25 * (-1) + 2 and -23 % 75 == 52
assert solutions_systeme(7, 15, 2, 25, 375) == [52, 127, 202, 277, 352] and solutions_systeme(1, 15, 3, 25, 375) == []
assert (1 - 3) % 5 != 0 and (127 - 52) % 375 != 0 and 2 % 15 != 7
bonne("arith-l1-chinois-ex4", 1, r"$29 + 36\mathbb{Z}$",
      [r"$29 + 216\mathbb{Z}$", "L'ensemble vide, car $12$ et $18$ ne sont pas premiers entre eux", r"$11 + 36\mathbb{Z}$"])
bonne("arith-l1-chinois-ex4", 2, r"$52 + 75\mathbb{Z}$",
      [r"$52 + 375\mathbb{Z}$", "L'ensemble vide, car $15$ et $25$ ne sont pas premiers entre eux", r"$2 + 75\mathbb{Z}$"])

# --- arith-l1-racines-polynome-ex1, v1 : P = X³ - X, P(3) = 24 ; v2 : P = X² - 2, P(0) = -2.
assert expand(interpolate([(t, t ** 3 - t) for t in (-1, 0, 1, 2)], X) - (X ** 3 - X)) == 0 and 3 ** 3 - 3 == 24
assert expand(interpolate([(t, t ** 2 - 2) for t in (1, 2, 3, 4)], X) - (X ** 2 - 2)) == 0
assert 24 not in (0, 6) and -2 not in (0, 2)
bonne("arith-l1-racines-polynome-ex1", 1, "$24$", ["$0$", "$6$", "On ne peut pas le déterminer avec ces données"])
bonne("arith-l1-racines-polynome-ex1", 2, "$-2$", ["$0$", "$2$", "On ne peut pas le déterminer avec ces données"])

# --- arith-l1-racines-polynome-ex2, v1 : 4 racines de X² - 1 dans Z/12Z ; v2 : 8 racines dans Z/24Z.
assert [t for t in range(12) if (t * t) % 12 == 1] == [1, 5, 7, 11] and 81 == 6 * 12 + 9 and (3 * 4) % 12 == 0
assert [t for t in range(24) if (t * t) % 24 == 1] == [1, 5, 7, 11, 13, 17, 19, 23] and (4 * 6) % 24 == 0
assert 25 == 2 * 12 + 1 and 49 == 4 * 12 + 1 and 121 == 10 * 12 + 1 and 25 == 24 + 1 and 49 == 2 * 24 + 1 and 121 == 5 * 24 + 1
bonne("arith-l1-racines-polynome-ex2", 1, "$4$", ["$2$", "$6$", "$3$"])
bonne("arith-l1-racines-polynome-ex2", 2, "$8$", ["$2$", "$4$", "$12$"])

# --- arith-l1-racines-polynome-ex3, v1 : Z/21Z → {1, 8, 13, 20} ; v2 : Z/35Z → {1, 6, 29, 34}.
assert [t for t in range(21) if (t * t) % 21 == 1] == [1, 8, 13, 20] and 64 == 3 * 21 + 1 and (7 * 9) % 21 == 0
assert (13 % 3, 13 % 7, 8 % 3, 8 % 7) == (1, 6, 2, 1) and (14 * 14) % 21 != 1
assert [t for t in range(35) if (t * t) % 35 == 1] == [1, 6, 29, 34] and 36 == 35 + 1 and (5 * 7) % 35 == 0
assert (6 % 5, 6 % 7, 29 % 5, 29 % 7) == (1, 6, 4, 1) and (11 * 11) % 35 != 1
ENS = lambda *cs: "$\\{" + ", ".join(rf"\overline{{{c}}}" for c in cs) + "\\}$"
bonne("arith-l1-racines-polynome-ex3", 1, ENS(1, 8, 13, 20), [ENS(1, 20), ENS(1, 8, 20), ENS(1, 8, 14, 20)])
bonne("arith-l1-racines-polynome-ex3", 2, ENS(1, 6, 29, 34), [ENS(1, 34), ENS(1, 6, 34), ENS(1, 6, 11, 34)])

# --- arith-l1-racines-polynome-ex4, v1 : 16! ≡ 16 [17] ; v2 : 17! ≡ 1 [19].
assert factorial(16) % 17 == 16 and factorial(18) % 19 == 18 and factorial(17) % 19 == 1
bonne("arith-l1-racines-polynome-ex4", 1, "$16$", ["$1$", "$0$", "$8$"])
bonne("arith-l1-racines-polynome-ex4", 2, "$1$", ["$18$", "$0$", "$17$"])

# Les exercices d'origine gardent leur bonne réponse (variante 0).
bonne("arith-t-division-euclidienne-ex1", 0, "$q = -290$, $r = 4$", ["$q = -289$, $r = -3$", "$q = -289$, $r = 3$", "$q = -290$, $r = -4$"])

# --- Phase 3 (relecture enr-arithmetique) : champs « cours » et introductions de chapitre.
# arith-t-bezout (cours) : 2 × 2 + 1 × 1 = 5 alors que 2 ∧ 1 = 1 ; si d > 0 divise a et b et d = au + bv, alors d = a ∧ b.
assert 2 * 2 + 1 * 1 == 5 and gcd(2, 1) == 1
for a in range(1, 30):
    for b in range(1, 30):
        g = gcd(a, b)
        assert all(not (a % d == 0 and b % d == 0 and any((d - a * u) % b == 0 for u in range(b))) or d == g for d in range(1, 30))
# arith-t-fermat (cours) : test de Fermat seulement si n ∤ a (5 premier mais 5^4 ≡ 0 [5]) ; 2^340 ≡ 1 [341], 341 composé.
assert isprime(5) and pow(5, 4, 5) == 0 and pow(2, 340, 341) == 1 and not isprime(341)
for nn in range(2, 400):
    for a in range(1, 12):
        if a % nn and pow(a, nn - 1, nn) != 1:
            assert not isprime(nn)
# arith-t-gauss (cours) : 6 | 4 × 3, 6 ∤ 4, 6 ∤ 3.
assert (4 * 3) % 6 == 0 and 4 % 6 and 3 % 6
# arith-l1-chinois (cours) : x ≡ 0 [2] et x ≡ 1 [4] sans solution.
assert solutions_systeme(0, 2, 1, 4, 100) == []
# arith-l1-racines-polynome (cours) : X² - 1 a 4 racines dans Z/8Z.
assert [t for t in range(8) if (t * t) % 8 == 1] == [1, 3, 5, 7]
# Introduction arith-l1-congruences : Euler exige a ∧ n = 1 (2^φ(4) ≡ 0 [4]) ; seuls ±1 sont leur propre inverse modulo p.
assert pow(2, int(totient(4)), 4) == 0
for nn in range(2, 60):
    for a in range(1, nn):
        if gcd(a, nn) == 1:
            assert pow(a, int(totient(nn)), nn) == 1
for pp in [3, 5, 7, 11, 13, 17, 19, 23]:
    assert [t for t in range(1, pp) if (t * t) % pp == 1] == [1, pp - 1]

# Titres et indices hérités par les variantes : tout nombre cité figure dans l'énoncé de la variante.
import re
for ident, e in EXOS.items():
    for v in e.get("variantes", []):
        herite = ([] if "titre" in v else [e["titre"]]) + ([] if "indices" in v else e["indices"])
        for t in herite:
            for nombre in re.findall(r"\d+", t):
                assert nombre in ("1", "2") or nombre in re.findall(r"\d+", v["enonce"]), (ident, t, nombre)

print("Toutes les vérifications du thème arithmetique passent.")
