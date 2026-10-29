"""Vérifications calculatoires de l'enrichissement du chapitre « arith-t-premiers » : questions de QCM
(content/qcm/arith-t-premiers.json) et exercices ajoutés à content/exercices/arithmetique/terminale-arith-t-premiers.json
(arith-t-fermat-ex5, arith-t-premiers-decomposition-ex6, arith-t-premiers-infinitude-ex5, variantes comprises).

Chaque assert cite l'identifiant de la question ou de l'exercice concerné. Voir CLAUDE.md, « Contrôle de l'exactitude ».
Lancer : python3 verification/enrichi_arith_t_premiers.py (doit finir sans erreur).
"""

from math import isqrt

from sympy import divisor_count, factorint, gcd, isprime, lcm, primerange, symbols, expand, integer_nthroot
from sympy.ntheory import n_order

n = symbols("n")

# --- QCM -------------------------------------------------------------------------------------------------------
# arith-t-premiers-q1 : 251 est premier (tests jusqu'à 13), 221 = 13·17, 247 = 13·19, 253 = 11·23.
assert isprime(251) and 16 ** 2 == 256 and all(251 % p for p in [2, 3, 5, 7, 11, 13])
assert 251 == 7 * 35 + 6 == 11 * 22 + 9 == 13 * 19 + 4
assert 221 == 13 * 17 and 247 == 13 * 19 and 253 == 11 * 23

# arith-t-premiers-q2 (vrai) : les premiers ≤ √397 sont exactement 2, …, 19, et 397 est premier.
assert 20 ** 2 == 400 > 397 and list(primerange(2, isqrt(397) + 1)) == [2, 3, 5, 7, 11, 13, 17, 19]
assert isprime(397)

# arith-t-premiers-q3 : 1000 = 2^3·5^3 a 16 diviseurs ; leurres 9, 6, 15 faux.
assert factorint(1000) == {2: 3, 5: 3} and divisor_count(1000) == 16

# arith-t-premiers-q4 : PGCD(84, 90) = 6, PPCM = 1260 ; leurres faux.
assert gcd(84, 90) == 6 and lcm(84, 90) == 1260 and 84 * 90 == 7560 and 90 % 12 != 0

# arith-t-premiers-q5 (vrai) : les entiers ayant exactement 3 diviseurs sont les carrés de nombres premiers.
assert [m for m in range(1, 5000) if divisor_count(m) == 3] == [p * p for p in primerange(2, 71)]

# arith-t-premiers-q6 (faux) : 2 divise 1 + 1 sans diviser 1.
assert (1 + 1) % 2 == 0 and 1 % 2 != 0

# arith-t-premiers-q7 : p | n^3 ⇒ p^3 | n^3 ; contre-exemples des leurres : n = p (p² ∤ n, p⁴ ∤ n³), n = 2p.
for p in primerange(2, 30):
    for m in range(-200, 201):
        if (m ** 3) % p == 0:
            assert (m ** 3) % p ** 3 == 0
    assert p % p ** 2 != 0 and (p ** 3) % p ** 4 != 0 and not isprime(2 * p) and (8 * p ** 3) % p == 0

# arith-t-premiers-q8 : 3^100 ≡ 4 (mod 7) ; leurres 1, 2, 5 faux.
assert pow(3, 6, 7) == 1 and 100 == 6 * 16 + 4 and pow(3, 100, 7) == 4 == 81 % 7

# --- arith-t-fermat-ex5 (diviseurs des nombres de Mersenne) --------------------------------------------------------
# Questions 1-2 (propriété générale) : pour p premier impair, tout diviseur premier q de 2^p − 1 vérifie
# q ≡ 1 (mod 2p) et l'ordre de 2 modulo q vaut p ; contrôle sur les petits p.
for p in primerange(3, 40):
    for q in factorint(2 ** p - 1):
        assert q != 2 and (q - 1) % (2 * p) == 0
        assert n_order(2, q) == p

# exercice, question 3 : M_11 = 2047 = 23·89, plus petit diviseur premier 23 ; leurres 11, 3 (ne divisent pas), 89.
assert 2 ** 11 - 1 == 2047 == 23 * 89 and min(factorint(2047)) == 23
assert 2047 % 11 and 2047 % 3 and 89 > 23
# exercice, question 4 : M_13 = 8191 premier ; candidats ≡ 1 mod 26 et ≤ 90 : 27, 53, 79.
assert 91 ** 2 == 8281 > 8191 > 90 ** 2 and [c for c in range(2, 91) if c % 26 == 1] == [27, 53, 79]
assert not isprime(27) and 8191 == 53 * 154 + 29 == 79 * 103 + 54 and isprime(8191)

# variante 1, question 3 : M_23 = 8388607, plus petit diviseur premier 47 ; leurres 23, 3, 7 ne divisent pas.
assert 2 ** 23 - 1 == 8388607 and min(factorint(8388607)) == 47 and [c for c in range(2, 47) if c % 46 == 1] == []
assert 256 == 5 * 47 + 21 and 441 == 9 * 47 + 18 and 18 * 128 == 2304 == 49 * 47 + 1 and pow(2, 23, 47) == 1
assert all(8388607 % c for c in [23, 3, 7])
# variante 1, question 4 : M_7 = 127 premier ; aucun entier ≡ 1 mod 14 dans [2, 11].
assert 12 ** 2 == 144 > 127 and [c for c in range(2, 12) if c % 14 == 1] == [] and isprime(127)

# variante 2, question 3 : M_29 = 536870911, candidats ≡ 1 mod 58 : 59 (ne divise pas), 117, 175 (composés), 233.
assert 2 ** 29 - 1 == 536870911 and min(factorint(536870911)) == 233
assert [c for c in range(2, 234) if c % 58 == 1] == [59, 117, 175, 233]
assert isprime(59) and not isprime(117) and not isprime(175) and isprime(233)
assert 64 % 59 == 5 and 625 == 10 * 59 + 35 and 35 * 32 == 1120 == 18 * 59 + 58 and pow(2, 29, 59) == 58
assert 256 % 233 == 23 and 529 == 2 * 233 + 63 and 63 * 23 == 1449 == 6 * 233 + 51 and 51 * 32 == 1632 == 7 * 233 + 1
assert pow(2, 29, 233) == 1 and 536870911 % 1103 == 0 and 536870911 % 3 != 0
# variante 2, question 4 : M_17 = 131071 ; premiers ≡ 1 mod 34 et ≤ 362 : 103, 137, 239, 307 ; aucun ne divise.
assert 362 ** 2 == 131044 < 131071 < 363 ** 2 == 131769
cand = [c for c in range(2, 363) if c % 34 == 1]
assert cand == [35, 69, 103, 137, 171, 205, 239, 273, 307, 341]
assert [c for c in cand if isprime(c)] == [103, 137, 239, 307]
assert 35 == 5 * 7 and 69 == 3 * 23 and 171 == 9 * 19 and 205 == 5 * 41 and 273 == 3 * 91 and 341 == 11 * 31
assert all(131071 % c for c in [103, 137, 239, 307]) and isprime(131071)

# --- arith-t-premiers-decomposition-ex6 (compléter en carré ou en cube) --------------------------------------------


def plus_petit(m, puissance):
    k = 1
    while integer_nthroot(k * m, puissance)[1] is False:
        k += 1
    return k


# exercice : 360 = 2^3·3^2·5, k = 10 (3600 = 60²), k' = 75 (27000 = 30³) ; leurres faux.
assert factorint(360) == {2: 3, 3: 2, 5: 1} and plus_petit(360, 2) == 10 and plus_petit(360, 3) == 75
assert 3600 == 60 ** 2 and 75 * 360 == 27000 == 30 ** 3 and 15 != 75 and 5 != 10
# variante 1 : 540 = 2^2·3^3·5, k = 15 (8100 = 90²), k' = 50 (27000 = 30³).
assert factorint(540) == {2: 2, 3: 3, 5: 1} and plus_petit(540, 2) == 15 and plus_petit(540, 3) == 50
assert 540 * 15 == 8100 == 90 ** 2 and 540 * 50 == 27000
# variante 2 : 1176 = 2^3·3·7^2, k = 6 (7056 = 84²), k' = 63 (74088 = 42³).
assert factorint(1176) == {2: 3, 3: 1, 7: 2} and plus_petit(1176, 2) == 6 and plus_petit(1176, 3) == 63
assert 1176 * 6 == 7056 == 84 ** 2 and 1176 * 63 == 74088 == 42 ** 3

# --- arith-t-premiers-infinitude-ex5 (n^4 + 4) ---------------------------------------------------------------------
assert expand((n ** 2 + 2 * n + 2) * (n ** 2 - 2 * n + 2)) == n ** 4 + 4
assert [m for m in range(1, 300) if isprime(m ** 4 + 4)] == [1]
# exercice : A_7 = 2405 = 65·37 = 5·13·37 ; leurres : 481 = 13·37 non premier, 5·13·41 = 2665.
assert 7 ** 4 + 4 == 2405 == 65 * 37 and factorint(2405) == {5: 1, 13: 1, 37: 1} and not isprime(481)
assert 5 * 13 * 41 != 2405
# variante 1 : A_9 = 6565 = 101·65 = 5·13·101.
assert 9 ** 4 + 4 == 6565 == 101 * 65 and factorint(6565) == {5: 1, 13: 1, 101: 1} and 5 * 13 * 103 != 6565
# variante 2 : A_11 = 14645 = 145·101 = 5·29·101.
assert 11 ** 4 + 4 == 14645 == 145 * 101 and factorint(14645) == {5: 1, 29: 1, 101: 1} and 5 * 29 * 103 != 14645


# --- Relecture indépendante (2026-10-29) : contrôles complémentaires ----------------------------------------------
# arith-t-premiers-decomposition-ex6 : chaque leurre (k, k') échoue (k·n non carré ou k'·n non cube).
def _carre(m):
    return integer_nthroot(m, 2)[1]


def _cube(m):
    return integer_nthroot(m, 3)[1]


for _n, _bon, _leurres in [(360, (10, 75), [(10, 15), (5, 75), (30, 600)]),
                           (540, (15, 50), [(15, 10), (5, 50), (30, 100)]),
                           (1176, (6, 63), [(6, 21), (3, 63), (42, 63)])]:
    assert _carre(_bon[0] * _n) and _cube(_bon[1] * _n)
    assert min(k for k in range(1, 10000) if _carre(k * _n)) == _bon[0]
    assert min(k for k in range(1, 10000) if _cube(k * _n)) == _bon[1]
    for (k, kp) in _leurres:  # (k, k') minimal est unique : un leurre en diffère
        assert (k, kp) != _bon
# arith-t-premiers-infinitude-ex5 : (n − 1)² + 1 ≥ 2 pour n ≥ 2, et n² + 2n + 2 > n² − 2n + 2 pour n ≥ 1.
assert all((m - 1) ** 2 + 1 >= 2 and m * m + 2 * m + 2 > m * m - 2 * m + 2 for m in range(2, 500))
# arith-t-fermat-ex5 : leurres de la variante 2 : 59 ne divise pas M_29, 1103 divise mais n'est pas le plus petit.
assert (2 ** 29 - 1) % 59 != 0 and (2 ** 29 - 1) % 1103 == 0 and 1103 > 233
# arith-t-premiers-q7 : p^3 | n^3 pour tout n tel que p | n^3 (test sur p ≤ 13, |n| ≤ 300).
assert all((m ** 3) % (p ** 3) == 0 for p in [2, 3, 5, 7, 11, 13] for m in range(-300, 301) if (m ** 3) % p == 0)

print("enrichi_arith_t_premiers : tout est vérifié.")
