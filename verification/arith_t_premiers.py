"""Vérifications calculatoires du chapitre « arith-t-premiers » (cartes et exercices du fichier
content/{cartes,exercices}/arithmetique/terminale-arith-t-premiers.json : infinité des nombres premiers,
décomposition en facteurs premiers, lemme d'Euclide).

Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise) concerné.
Lancer : python3 verification/arith_t_premiers.py (doit finir sans erreur).
"""

from itertools import product
from math import factorial, gcd, isqrt, log

from sympy import divisor_count, divisors, factorint, isprime, primepi, primerange


def plus_petit_diviseur(n):
    return min(d for d in range(2, n + 1) if n % d == 0)


# ---------------------------------------------------------------------------
# arith-t-premiers-infinitude (carte)
# ---------------------------------------------------------------------------

# Énoncé : le plus petit diviseur p ≥ 2 de n ≥ 2 est premier, et p ≤ √n si n n'est pas premier.
for m in range(2, 5000):
    p = plus_petit_diviseur(m)
    assert isprime(p)
    if not isprime(m):
        assert p * p <= m
# Cours : 25 nombres premiers inférieurs à 100, 143 entre 100 et 1000.
assert primepi(100) == 25 and primepi(1000) - primepi(100) == 143
# Cours / contre-exemple / affirmation 1 : 2·3·5·7·11·13 + 1 = 30031 = 59 × 509, non premier (59, 509 premiers).
assert 2 * 3 * 5 * 7 * 11 * 13 + 1 == 30031 == 59 * 509 and isprime(59) and isprime(509)
# Exemple : 211 premier ; 14² = 196 ≤ 211 < 225 = 15² ; restes cités ; 211 = 2·3·5·7 + 1.
assert isprime(211) and 14 ** 2 == 196 <= 211 < 225 == 15 ** 2 and 2 * 3 * 5 * 7 + 1 == 211
assert [divmod(211, p) for p in (2, 3, 5, 7, 11, 13)] == [(105, 1), (70, 1), (42, 1), (30, 1), (19, 2), (16, 3)]
# Contre-exemples : 1 n'a pas de diviseur premier ; n = 7 : plus petit diviseur 7 > √7.
assert [p for p in primerange(2, 100) if 1 % p == 0] == [] and plus_petit_diviseur(7) == 7 and 7 * 7 > 7
# Leurres {{1}} : n ≥ 1, n ∈ Z, n ∈ N contiennent 1 (ou -1), sans diviseur premier.
# Leurres {{2}} : n = 8 : p = 2 < √8 (donc « p ≥ √n » faux) ; n = 4 : p = 2 > 4/3 et 2 > ln 4.
assert plus_petit_diviseur(8) == 2 and 2 * 2 < 8
assert plus_petit_diviseur(4) == 2 and 2 > 4 / 3 and 2 > log(4)
# Leurres {{3}} : p_1 divise p_1⋯p_k, p_1 + ⋯ + p_k pour k = 1, et p_1⋯p_k - p_1.
for liste in ([2], [2, 3], [3, 5, 7], [2, 3, 5, 7, 11]):
    prod = 1
    for p in liste:
        prod *= p
    assert prod % liste[0] == 0 and (prod - liste[0]) % liste[0] == 0
    assert all((prod + 1) % p != 0 for p in liste)  # bonne réponse
assert sum([2]) % 2 == 0
# Pièges : « n impair » superflu (n = 6 a le diviseur premier 2) ; majoration avec n impair premier fausse (n = 7).
assert plus_petit_diviseur(6) == 2
# Affirmation 2 (vraie) : couverte par la boucle ci-dessus (contraposée).
# Affirmation 3 (fausse) : 40² + 40 + 41 = 1681 = 41².
assert 40 ** 2 + 40 + 41 == 1681 == 41 ** 2
# Affirmation 4 (fausse) : 2^11 - 1 = 2047 = 23 × 89, avec 11 premier.
assert 2 ** 11 - 1 == 2047 == 23 * 89 and isprime(11)

# ---------------------------------------------------------------------------
# arith-t-premiers-decomposition (carte)
# ---------------------------------------------------------------------------
# Énoncé : nombre de diviseurs = produit des (α_i + 1), diviseurs = produits à exposants ≤ α_i.
for m in range(2, 3000):
    f = factorint(m)
    nb = 1
    for a in f.values():
        nb *= a + 1
    assert nb == divisor_count(m)
    for d in divisors(m):
        fd = factorint(d)
        assert all(p in f and fd[p] <= f[p] for p in fd)
# Exemple : 360 = 2³ × 3² × 5, 24 diviseurs ; 84 = 2² × 3 × 7 ; 360 ∧ 84 = 12.
assert factorint(360) == {2: 3, 3: 2, 5: 1} and divisor_count(360) == 24 and (3 + 1) * (2 + 1) * (1 + 1) == 24
assert 360 == 2 * 180 == 2 ** 2 * 90 == 2 ** 3 * 45 and factorint(84) == {2: 2, 3: 1, 7: 1} and gcd(360, 84) == 12 == 2 ** 2 * 3
# Contre-exemples : 12 = 2² × 3 = 2² × 3 × 5^0 ; 12 = 3 × 4 avec 3 ∧ 4 = 1.
assert 12 == 2 ** 2 * 3 == 2 ** 2 * 3 * 5 ** 0 == 3 * 4 and gcd(3, 4) == 1
# Contre-exemple « entiers pairs » : les pairs indécomposables (non produits de deux pairs) sont les 2 mod 4 ;
# 2, 6, 10, 30 en sont, et 60 = 2 × 30 = 6 × 10.


def indecomposable_pair(m):
    return m % 2 == 0 and not any(m % a == 0 and a % 2 == 0 and (m // a) % 2 == 0 for a in range(2, m))


assert all(indecomposable_pair(m) == (m % 4 == 2) for m in range(2, 400, 2))
assert all(indecomposable_pair(m) for m in (2, 6, 10, 30)) and 60 == 2 * 30 == 6 * 10
# Leurres {{1}} : -2 (dans Z*), 0 (dans N et Z) n'ont pas de telle écriture (produit de premiers > 0, non nul).
# Leurres {{2}} : 12 = 12¹ = 3 × 4 (deux à deux premiers entre eux) ; 2 n'est pas un produit d'impairs ;
# 4 = 4¹ = 2² (entiers ≥ 2 quelconques).
assert 2 % 2 == 0 and 4 == 2 ** 2
# Leurres {{3}} : exposants naturels ou relatifs : 2 = 2¹ = 2¹ × 3⁰ ; exposants premiers : 2 = p^α avec α premier impossible.
assert not any(p ** a == 2 for p in (2, 3, 5) for a in (2, 3, 5, 7))
# Piège « n non premier » superflu : 7 = 7¹.
assert factorint(7) == {7: 1}
# Affirmation 1 (fausse) : 2³ × 3² × 5 = 360 a 24 diviseurs, et 1, 2, 3, 4, 5, 6, 8 en sont déjà sept.
assert divisor_count(2 ** 3 * 3 ** 2 * 5) == 24 != 6 and all(360 % d == 0 for d in (1, 2, 3, 4, 5, 6, 8))
# Affirmation 2 (vraie) : carré <=> exposants pairs.
for m in range(2, 5000):
    assert (isqrt(m) ** 2 == m) == all(a % 2 == 0 for a in factorint(m).values())
# Affirmation 3 (fausse) : PGCD de 2³·3·5² et 2·3²·7 = 6 ; le produit proposé est le PPCM et ne divise pas 2·3²·7.
x, y = 2 ** 3 * 3 * 5 ** 2, 2 * 3 ** 2 * 7
assert gcd(x, y) == 6 and x * y // gcd(x, y) == 2 ** 3 * 3 ** 2 * 5 ** 2 * 7 and y % (2 ** 3 * 3 ** 2 * 5 ** 2 * 7) != 0

# ---------------------------------------------------------------------------
# arith-t-premiers-euclide (carte)
# ---------------------------------------------------------------------------
# Énoncé et réciproque : p ≥ 2 vérifie (p | ab => p | a ou p | b) pour tous a, b si et seulement si p est premier.
for p in range(2, 150):
    propriete = all(not ((a * b) % p == 0) or a % p == 0 or b % p == 0 for a in range(1, p + 1) for b in range(1, p + 1))
    assert propriete == isprime(p)
# Plusieurs facteurs : p | a1 a2 a3 => p divise l'un d'eux (affirmation 4).
for p in (2, 3, 5, 7):
    for a1, a2, a3 in product(range(1, 15), repeat=3):
        if (a1 * a2 * a3) % p == 0:
            assert a1 % p == 0 or a2 % p == 0 or a3 % p == 0
# Exemple : a² = 2b² n'a pas de solution avec a, b ≥ 1 (√2 irrationnel ; contrôle pour b < 3000).
assert not any(isqrt(2 * b * b) ** 2 == 2 * b * b for b in range(1, 3000))
# Contre-exemples : 6 | 2 × 3 ; 4 | 2 × 2 ; 2 | 1 + 1 ; 3² ≡ 1 (mod 8), 3 ≢ ±1 (mod 8).
assert (2 * 3) % 6 == 0 and 2 % 6 and 3 % 6 and (2 * 2) % 4 == 0 and 2 % 4
assert (1 + 1) % 2 == 0 and 1 % 2 and pow(3, 2, 8) == 1 and 3 % 8 not in (1, 7)
# Leurres {{1}} : p = 4 (≥ 2, carré parfait ≥ 4), p = 9 (impair ≥ 3) : 4 | 2·2, 9 | 3·3.
assert (2 * 2) % 4 == 0 and 2 % 4 and (3 * 3) % 9 == 0 and 3 % 9
# Leurres {{2}} : p = 2, a = b = 1 : 2 | 1 + 1, 2 | 1² + 1², 2 | 1·1 + 1, mais 2 ∤ 1.
assert (1 + 1) % 2 == 0 and (1 + 1) % 2 == 0 and (1 * 1 + 1) % 2 == 0 and 1 % 2 == 1
# Leurres {{3}} : 2 | 2·1 mais 2 ∤ 1 (« et », « p | a » avec a = 1, b = 2) ; 4 ∤ 2·1 (« p² | ab »).
assert (2 * 1) % 2 == 0 and 1 % 2 and (2 * 1) % 4
# Affirmation 1 (vraie) : p | a² => p² | a² (contrôle). Affirmation 2 (fausse) : 4 | 2², 4 ∤ 2.
for p in (2, 3, 5, 7, 11):
    assert all((a * a) % (p * p) == 0 for a in range(1, 500) if (a * a) % p == 0)
assert (2 ** 2) % 4 == 0 and 2 % 4
# Affirmation 3 (fausse) : 2 | 2·2, 2 | 2 et 2 | 2.
assert (2 * 2) % 2 == 0 and 2 % 2 == 0
# Pièges : « a ∧ b = 1 » superflu (2 | 2·4 avec 2 ∧ 4 = 2 : conclusion vraie) ; « p | a + b » faux (2 | 1 + 1).
assert gcd(2, 4) == 2 and (2 * 4) % 2 == 0

# ---------------------------------------------------------------------------
# arith-t-premiers-infinitude-ex1 (et variantes) : essais de divisions
# ---------------------------------------------------------------------------
# Exercice : 391 = 17 × 23 ; 19² = 361 ≤ 391 < 400 ; restes cités.
assert plus_petit_diviseur(391) == 17 and 391 == 17 * 23 and 19 ** 2 == 361 <= 391 < 400
assert divmod(391, 7) == (55, 6) and divmod(391, 11) == (35, 6) and divmod(391, 13) == (30, 1) and 3 + 9 + 1 == 13
assert all(x != 17 for x in (391, 13, 23)) and not isprime(391)
# Variante 1 : 409 premier ; 20² = 400 ≤ 409 < 441 ; restes cités.
assert isprime(409) and 400 <= 409 < 441
assert [divmod(409, p) for p in (7, 11, 13, 17, 19)] == [(58, 3), (37, 2), (31, 6), (24, 1), (21, 10)] and 4 + 0 + 9 == 13
assert all(409 % x != 0 for x in (7, 17, 19))
# Variante 2 : 527 = 17 × 31 ; 22² = 484 ≤ 527 < 529 ; restes cités.
assert plus_petit_diviseur(527) == 17 and 527 == 17 * 31 and 484 <= 527 < 529
assert [divmod(527, p) for p in (7, 11, 13)] == [(75, 2), (47, 10), (40, 7)] and 5 + 2 + 7 == 14
assert not isprime(527) and all(x != 17 for x in (31, 13))

# ---------------------------------------------------------------------------
# arith-t-premiers-infinitude-ex2 (et variantes) : nombres d'Euclide non premiers
# ---------------------------------------------------------------------------
for (liste, N, facteurs, div_leurre) in [([2, 3, 5, 7, 11, 13], 30031, {59: 1, 509: 1}, 13),   # exercice
                                        ([3, 5, 7, 11], 1156, {2: 2, 17: 2}, 11),            # variante 1
                                        ([2, 7, 11], 155, {5: 1, 31: 1}, 7)]:                # variante 2
    prod = 1
    for p in liste:
        prod *= p
    assert prod + 1 == N and factorint(N) == facteurs and not isprime(N)
    assert all(N % p != 0 for p in liste) and all(q not in liste for q in facteurs)
    assert N % div_leurre != 0  # leurre « N divisible par … » faux
    assert plus_petit_diviseur(N) not in liste  # leurre « plus petit diviseur dans la liste » faux
assert 30030 + 1 == 30031 and 1155 + 1 == 1156 and 154 + 1 == 155 and 4 * 289 == 1156

# ---------------------------------------------------------------------------
# arith-t-premiers-infinitude-ex3 (et variantes) : diviseurs premiers de n! + 1
# ---------------------------------------------------------------------------
for m in range(1, 30):
    assert all(p > m for p in factorint(factorial(m) + 1))
for (m, N, ppd, leurres) in [(6, 721, 7, [103, 5]), (5, 121, 11, [7, 3]), (10, 3628801, 11, [13, 7])]:
    assert factorial(m) + 1 == N and plus_petit_diviseur(N) == ppd and not isprime(N)  # exercice, variantes 1, 2
    assert all(L != ppd for L in leurres)
assert 721 == 7 * 103 and isprime(103) and 10 ** 2 < 103 < 11 ** 2 and all(103 % p for p in (2, 3, 5, 7))
assert divmod(121, 7) == (17, 2) and 121 == 11 ** 2 and 3628801 == 11 * 329891 and factorial(10) == 3628800

# ---------------------------------------------------------------------------
# arith-t-premiers-infinitude-ex4 : nombres premiers congrus à 3 modulo 4
# ---------------------------------------------------------------------------
# Question 2 : tout m ≥ 2, m ≡ 3 (mod 4), a un diviseur premier ≡ 3 (mod 4).
assert all(any(p % 4 == 3 for p in factorint(m)) for m in range(3, 20000, 4))
# Question 3 : 4 × 3 - 1 = 11 ; 4·p1⋯pk - 1 ≡ 3 (mod 4) et premier avec chaque p_i.
assert 4 * 3 - 1 == 11 and all((4 * P - 1) % 4 == 3 for P in (3, 21, 231, 3 * 7 * 11 * 19))
# Question 4 : 21 = 3 × 7 et 9 = 3² sont ≡ 1 (mod 4) sans diviseur premier ≡ 1 (mod 4).
assert 21 % 4 == 1 and 9 % 4 == 1 and all(p % 4 == 3 for p in factorint(21)) and all(p % 4 == 3 for p in factorint(9))

# ---------------------------------------------------------------------------
# arith-t-premiers-decomposition-ex1 (et variantes)
# ---------------------------------------------------------------------------
assert factorint(1764) == {2: 2, 3: 2, 7: 2} and 1764 == 42 ** 2 == 2 ** 2 * 21 ** 2  # exercice
assert 1764 == 2 * 882 and 882 == 2 * 441 and 441 == 3 * 147 and 147 == 3 * 49
assert 2 * 3 ** 2 * 7 ** 2 != 1764 and 2 ** 2 * 3 * 7 ** 2 != 1764 and not isprime(21)
assert factorint(2250) == {2: 1, 3: 2, 5: 3} and 2 * 15 ** 2 * 5 == 2250  # variante 1
assert 2250 == 2 * 1125 and 1125 == 3 * 375 and 375 == 3 * 125
assert 2 * 3 * 5 ** 3 != 2250 and 2 ** 2 * 3 ** 2 * 5 ** 2 != 2250 and not isprime(15)
assert factorint(3136) == {2: 6, 7: 2} and 3136 == 56 ** 2 == 8 ** 2 * 7 ** 2  # variante 2
assert [3136 // 2 ** i for i in range(1, 7)] == [1568, 784, 392, 196, 98, 49]
assert 2 ** 5 * 7 ** 2 != 3136 and 2 ** 6 * 7 != 3136 and not isprime(8)

# ---------------------------------------------------------------------------
# arith-t-premiers-decomposition-ex2 (et variantes) : nombre de diviseurs
# ---------------------------------------------------------------------------
for (N, f, nb, leurres) in [(720, {2: 4, 3: 2, 5: 1}, 30, [8, 10, 15]),     # exercice
                            (1400, {2: 3, 5: 2, 7: 1}, 24, [6, 9, 12]),     # variante 1
                            (2016, {2: 5, 3: 2, 7: 1}, 36, [10, 11, 18])]:  # variante 2
    assert factorint(N) == f and divisor_count(N) == nb
    exps = list(f.values())
    prod_alpha = 1
    for a in exps:
        prod_alpha *= a
    assert leurres[0] == prod_alpha and leurres[1] == sum(a + 1 for a in exps)
    assert leurres[2] == nb // (exps[-1] + 1)  # oubli du dernier facteur
    assert nb not in leurres

# ---------------------------------------------------------------------------
# arith-t-premiers-decomposition-ex3 (et variantes) : plus petit entier ayant k diviseurs
# ---------------------------------------------------------------------------
for (kk, bon, candidats, leurres) in [(12, 60, [2048, 96, 72, 60], [72, 96, 2048]),           # exercice
                                     (10, 48, [512, 48], [512, 80, 162]),                     # variante 1
                                     (18, 180, [131072, 768, 288, 180], [288, 768, 300])]:    # variante 2
    assert min(m for m in range(1, 200000) if divisor_count(m) == kk) == bon
    assert all(divisor_count(c) == kk for c in candidats) and min(candidats) == bon
    assert all(L != bon and L > bon for L in leurres)
assert divisors(60) == [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60]
assert divisors(48) == [1, 2, 3, 4, 6, 8, 12, 16, 24, 48]
assert divisors(180) == [1, 2, 3, 4, 5, 6, 9, 10, 12, 15, 18, 20, 30, 36, 45, 60, 90, 180]
assert 2 ** 11 == 2048 and 2 ** 5 * 3 == 96 and 2 ** 3 * 3 ** 2 == 72 and 2 ** 2 * 3 * 5 == 60
assert 2 ** 9 == 512 and 2 ** 4 * 3 == 48 and 2 ** 17 == 131072 and 2 ** 8 * 3 == 768 and 2 ** 5 * 9 == 288 and 4 * 9 * 5 == 180
# Échange d'exposants : (p/q)^(β - α) < 1 pour p < q, α < β (exemples).
assert 2 ** 3 * 3 ** 1 < 2 ** 1 * 3 ** 3 and 2 ** 2 * 5 < 2 * 5 ** 2

# ---------------------------------------------------------------------------
# arith-t-premiers-decomposition-ex4 (et variantes) : zéros de n!
# ---------------------------------------------------------------------------


def v(m, p):
    s, q = 0, p
    while q <= m:
        s += m // q
        q *= p
    return s


def zeros(m):
    z, x = 0, factorial(m)
    while x % 10 == 0:
        x //= 10
        z += 1
    return z


for m in range(1, 300):
    for p in (2, 3, 5, 7):
        assert factorint(factorial(m)).get(p, 0) == v(m, p)
    assert zeros(m) == min(v(m, 2), v(m, 5))
for (m, t5, t2, z, leurres) in [(100, [20, 4], [50, 25, 12, 6, 3, 1], 24, [20, 10, 97]),                              # exercice
                                (200, [40, 8, 1], [100, 50, 25, 12, 6, 3, 1], 49, [40, 48, 20]),                       # variante 1
                                (2026, [405, 81, 16, 3], [1013, 506, 253, 126, 63, 31, 15, 7, 3, 1], 505, [405, 486, 202])]:  # variante 2
    assert [m // 5 ** j for j in range(1, len(t5) + 1)] == t5 and m // 5 ** (len(t5) + 1) == 0
    assert [m // 2 ** j for j in range(1, len(t2) + 1)] == t2 and m // 2 ** (len(t2) + 1) == 0
    assert sum(t5) == v(m, 5) == z and zeros(m) == z and z not in leurres
assert v(100, 2) == 97 and v(200, 2) == 197 and v(2026, 2) == 2018
assert 405 + 81 == 486  # leurre de la variante 2 (oubli de 125 et 625)

# ---------------------------------------------------------------------------
# arith-t-premiers-decomposition-ex5 : ab carré et a ∧ b = 1
# ---------------------------------------------------------------------------
for a, b in product(range(1, 200), repeat=2):
    if gcd(a, b) == 1 and isqrt(a * b) ** 2 == a * b:
        assert isqrt(a) ** 2 == a and isqrt(b) ** 2 == b
assert 2 * 8 == 16 == 4 ** 2 and gcd(2, 8) == 2  # question 3
assert not any(isqrt(m * (m + 1)) ** 2 == m * (m + 1) for m in range(1, 100000))  # question 4

# ---------------------------------------------------------------------------
# arith-t-premiers-euclide-ex1 et ex2 (et variantes)
# ---------------------------------------------------------------------------
for p in (13, 11, 17):  # ex1 : exercice, variante 1, variante 2
    assert isprime(p) and all((m * m) % (p * p) == 0 for m in range(-1000, 1000) if (m * m) % p == 0)
assert 6 ** 2 == 36 and 36 % 4 == 0 and 36 % 16 != 0  # contre-exemple cité
for p in (3, 5, 7):  # ex2 : √3, √5, √7 irrationnels (a² = p b² sans solution, contrôle b < 3000)
    assert isprime(p) and not any(isqrt(p * b * b) ** 2 == p * b * b for b in range(1, 3000))
assert 9 == 3 ** 2 and 25 == 5 ** 2 and 49 == 7 ** 2

# ---------------------------------------------------------------------------
# arith-t-premiers-euclide-ex3 (et variantes) : x² ≡ 1
# ---------------------------------------------------------------------------


def racines_de_1(m):
    return [x for x in range(m) if (x * x - 1) % m == 0]


for p in (13, 17, 11):  # question 2 : exercice, variante 1, variante 2
    assert racines_de_1(p) == [1, p - 1]
for p in primerange(2, 200):  # question 1
    assert racines_de_1(p) == sorted({1 % p, (p - 1) % p})
assert racines_de_1(8) == [1, 3, 5, 7] and [x * x % 8 for x in range(8)] == [0, 1, 4, 1, 0, 1, 4, 1]  # exercice
assert racines_de_1(24) == [1, 5, 7, 11, 13, 17, 19, 23]  # variante 1
assert racines_de_1(15) == [1, 4, 11, 14]  # variante 2
assert [x for x in range(15) if x % 5 in (1, 4)] == [1, 4, 6, 9, 11, 14]
for (m, nb, leurres) in [(8, 4, [2, 1, 3]), (24, 8, [2, 4, 6]), (15, 4, [2, 1, 3])]:
    assert len(racines_de_1(m)) == nb and nb not in leurres
# Contre-exemples cités : 8 | 2 × 4, 24 | 4 × 6, 15 | 3 × 5.
assert (2 * 4) % 8 == 0 and 2 % 8 and 4 % 8 and (4 * 6) % 24 == 0 and (3 * 5) % 15 == 0
# Variante 1 : x impair => 8 | x² - 1 ; 3 ∤ x => 3 | x² - 1.
assert all((x * x - 1) % 8 == 0 for x in range(1, 1000, 2)) and all((x * x - 1) % 3 == 0 for x in range(1, 1000) if x % 3)

# ---------------------------------------------------------------------------
# arith-t-premiers-euclide-ex4 (et variantes) : carrés modulo p
# ---------------------------------------------------------------------------
for p in primerange(3, 200):
    carres = {x * x % p for x in range(1, p)}
    assert len(carres) == (p - 1) // 2
    assert carres == {x * x % p for x in range(1, (p - 1) // 2 + 1)}
    assert len({x * x % p for x in range(1, (p - 1) // 2 + 1)}) == (p - 1) // 2
    for x, y in product(range(p), repeat=2):
        assert ((x * x - y * y) % p == 0) == ((x - y) % p == 0 or (x + y) % p == 0)
for (p, carres, non_carres, bon, leurres, m4, x4) in [
    (11, [1, 3, 4, 5, 9], [2, 6, 7, 8, 10], 5, [2, 6, 7], 15, 4),                                     # exercice
    (13, [1, 3, 4, 9, 10, 12], [2, 5, 6, 7, 8, 11], 10, [2, 5, 7], 21, 8),                           # variante 1
    (17, [1, 2, 4, 8, 9, 13, 15, 16], [3, 5, 6, 7, 10, 11, 12, 14], 13, [3, 5, 6], 35, 6),           # variante 2
]:
    C = sorted({x * x % p for x in range(1, p)})
    assert C == carres and sorted(set(range(1, p)) - set(C)) == non_carres
    assert bon in C and all(L not in C for L in leurres)
    # question 4 : x4² ≡ 1 = 1² (mod m4) sans x4 ≡ ±1 (mod m4)
    assert (x4 * x4) % m4 == 1 and x4 % m4 not in (1, m4 - 1) and not isprime(m4)
assert [x * x % 11 for x in range(1, 6)] == [1, 4, 9, 5, 3]
assert [x * x % 13 for x in range(1, 7)] == [1, 4, 9, 3, 12, 10]
assert [x * x % 17 for x in range(1, 9)] == [1, 4, 9, 16, 8, 2, 15, 13]

print("arith_t_premiers : toutes les vérifications passent.")
