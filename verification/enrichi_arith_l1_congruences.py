"""Vérifications calculatoires de l'enrichissement du chapitre « arith-l1-congruences » : questions de QCM
(content/qcm/arith-l1-congruences.json) et exercices ajoutés à
content/exercices/arithmetique/l1-arith-l1-congruences.json (arith-l1-congruences-euler-ex6, arith-l1-chinois-ex5,
arith-l1-congruences-inversibles-ex6, variantes comprises).

Chaque assert cite l'identifiant de la question ou de l'exercice concerné. Voir CLAUDE.md, « Contrôle de l'exactitude ».
Lancer : python3 verification/enrichi_arith_l1_congruences.py (doit finir sans erreur).
"""

from math import factorial, gcd

from sympy import isprime, totient

# --- QCM -------------------------------------------------------------------------------------------------------
# arith-l1-congruences-q1 : l'inverse de 7 modulo 24 est 7 ; leurres 17, 5, 3 faux.
assert (7 * 7) % 24 == 1 and 49 == 2 * 24 + 1 and pow(7, -1, 24) == 7
assert [(7 * x) % 24 for x in (17, 5, 3)] == [23, 11, 21] and 7 * 17 == 119 == 5 * 24 - 1

# arith-l1-congruences-q2 : 6x ≡ 4 (mod 10) a 2 solutions dans {0, …, 9} : 4 et 9 ; leurres 0, 1, 6 faux.
assert [x for x in range(10) if (6 * x - 4) % 10 == 0] == [4, 9]
assert all(((6 * x - 4) % 10 == 0) == ((3 * x - 2) % 5 == 0) == (x % 5 == 4) for x in range(-200, 200))
assert (3 * 2) % 5 == 1 and gcd(6, 10) == 2

# arith-l1-congruences-q3 (faux) : dans Z/12Z, les inversibles 1, 5, 7, 11 sont tous leur propre inverse.
assert [x for x in range(12) if gcd(x, 12) == 1] == [1, 5, 7, 11]
assert [x for x in range(12) if (x * x) % 12 == 1] == [1, 5, 7, 11] and 25 % 12 == 1 and 49 % 12 == 1

# arith-l1-congruences-q4 : φ(360) = 96 ; leurres 72, 120, 144 faux.
assert totient(360) == 96 == totient(8) * totient(9) * totient(5) == 4 * 6 * 4

# arith-l1-congruences-q5 (vrai) : a ∧ 15 = 1 ⇒ a^4 ≡ 1 (mod 15), alors que φ(15) = 8.
assert totient(15) == 8 and all(pow(a, 4, 15) == 1 for a in range(-300, 300) if gcd(a, 15) == 1)
assert all(pow(a, 2, 3) == 1 for a in (1, 2)) and all(pow(a, 4, 5) == 1 for a in (1, 2, 3, 4))

# arith-l1-congruences-q6 : 7^2026 ≡ 9 (mod 10) ; leurres 1, 3, 7 faux.
assert totient(10) == 4 and 7 ** 4 == 2401 and 2026 == 4 * 506 + 2 and pow(7, 2026, 10) == 9

# arith-l1-congruences-q7 (vrai) : 14! ≡ 8 (mod 17).
assert isprime(17) and factorial(16) % 17 == 16 and (16 * 15) % 17 == 2 and factorial(14) % 17 == 8
assert (2 * 9) % 17 == 1 and (2 * 8) % 17 == 16

# arith-l1-congruences-q8 : plus petit x ≥ 0 avec x ≡ 2 (5) et x ≡ 3 (7) : 17 ; 52 solution non minimale ; 12, 23 faux.
assert [x for x in range(0, 100) if x % 5 == 2 and x % 7 == 3] == [17, 52, 87]
assert 12 % 7 == 5 and 23 % 5 == 3

# --- arith-l1-congruences-euler-ex6 (RSA miniature) ----------------------------------------------------------------


def verifier_rsa(p, q, e, d, c, x, leurres):
    n, ph = p * q, (p - 1) * (q - 1)
    assert isprime(p) and isprime(q) and p != q and totient(n) == ph  # question 1
    assert gcd(e, ph) == 1 and (e * d) % ph == 1 and 1 <= d < ph  # question 2
    assert (e * d - 1) % (p - 1) == 0 and (e * d - 1) % (q - 1) == 0  # question 4
    assert all(pow(y, e * d, n) == y % n for y in range(-3 * n, 3 * n))  # questions 3-4 : pour tout entier
    assert pow(x, e, n) == c and pow(c, d, n) == x  # question 5
    for l in leurres:  # leurres : messages dont le chiffré n'est pas c
        assert 0 <= l < n and pow(l, e, n) != c


# exercice : p = 11, q = 17, n = 187, φ = 160, e = 7, d = 23, c = 175 → x = 10 ; leurres 175, 43, 98.
verifier_rsa(11, 17, 7, 23, 175, 10, [175, 43, 98])
assert 160 == 7 * 22 + 6 and 23 * 7 - 160 == 1 and 7 * 23 == 161
assert 175 == 15 * 11 + 10 and 175 == 10 * 17 + 5 and 25 % 17 == 8 and 64 % 17 == 13 and 104 % 17 == 2
assert pow(5, 16, 17) == 1 and pow(5, 23, 17) == 10 and pow(175, 23, 11) == 10
assert 10 ** 7 == 187 * 53475 + 175
# variante 1 : p = 5, q = 23, n = 115, φ = 88, e = 3, d = 59, c = 52 → x = 8 ; leurres 52, 23, 37.
verifier_rsa(5, 23, 3, 59, 52, 8, [52, 23, 37])
assert 88 == 3 * 29 + 1 and 3 * 59 == 177 == 2 * 88 + 1 and 176 == 4 * 44 == 22 * 8
assert 52 % 5 == 2 and pow(2, 59, 5) == 3 and 59 == 4 * 14 + 3 and 52 % 23 == 6 and 59 == 22 * 2 + 15
assert 36 % 23 == 13 and 169 % 23 == 8 and 64 % 23 == 18 and (18 * 8) % 23 == 6 and (6 * 13) % 23 == 9
assert (9 * 6) % 23 == 8 and pow(6, 15, 23) == 8 and 8 ** 3 == 512 == 4 * 115 + 52
# variante 2 : p = 7, q = 13, n = 91, φ = 72, e = 5, d = 29, c = 23 → x = 4 ; leurres 23, 10, 17.
verifier_rsa(7, 13, 5, 29, 23, 4, [23, 10, 17])
assert 72 == 5 * 14 + 2 and 29 * 5 - 2 * 72 == 1 and 144 == 2 * 72 == 6 * 24 == 12 * 12
assert 23 % 7 == 2 and 29 == 6 * 4 + 5 and pow(2, 29, 7) == 4 and 32 % 7 == 4
assert 23 % 13 == 10 and 29 == 12 * 2 + 5 and 243 == 18 * 13 + 9 and pow(-3, 29, 13) == 4 and 4 ** 5 == 1024 == 11 * 91 + 23

# --- arith-l1-chinois-ex5 (trois congruences) ----------------------------------------------------------------------


def sols(mods, restes, lo, hi):
    return [N for N in range(lo, hi + 1) if all(N % m == r for m, r in zip(mods, restes))]


# exercice : N ≡ 3 (4), 2 (5), 4 (7), 100 ≤ N ≤ 300 → 207 ; solutions 67 + 140Z ; leurres 67, 247, 187.
assert sols((4, 5), (3, 2), 0, 19) == [7] and sols((4, 5, 7), (3, 2, 4), 0, 139) == [67]
assert sols((4, 5, 7), (3, 2, 4), 100, 300) == [207]
assert 67 == 16 * 4 + 3 == 13 * 5 + 2 == 9 * 7 + 4 and 20 % 7 == 6
assert 247 % 7 != 4 and 187 % 7 != 4 and not 100 <= 67
# variante 1 : N ≡ 2 (3), 4 (5), 5 (8), 200 ≤ N ≤ 300 → 269 ; classe 29 + 120Z ; leurres 149, 254, 239.
assert sols((3, 5), (2, 4), 0, 14) == [14] and sols((3, 5, 8), (2, 4, 5), 0, 119) == [29]
assert sols((3, 5, 8), (2, 4, 5), 200, 300) == [269] and 149 < 200
assert 254 % 8 != 5 and 239 % 8 != 5 and 254 % 3 == 2 and 239 % 5 == 4
# variante 2 : N ≡ 4 (5), 1 (6), 2 (7), 200 ≤ N ≤ 400 → 289 ; classe 79 + 210Z ; leurres 79, 259, 349.
assert sols((5, 6), (4, 1), 0, 29) == [19] and sols((5, 6, 7), (4, 1, 2), 0, 209) == [79]
assert sols((5, 6, 7), (4, 1, 2), 200, 400) == [289]
assert 259 % 7 != 2 and 349 % 7 != 2 and 259 % 6 == 1 and 349 % 5 == 4

# --- arith-l1-congruences-inversibles-ex6 (congruence à coefficient non inversible) --------------------------------


def lin(a, b, n):
    return [x for x in range(n) if (a * x - b) % n == 0]


# exercice : 12x ≡ 18 (30) ⇔ 2x ≡ 3 (5) ⇔ x ≡ 4 (5) : 6 solutions.
assert gcd(12, 30) == 6 and lin(12, 18, 30) == [4, 9, 14, 19, 24, 29]
assert all(((12 * x - 18) % 30 == 0) == ((2 * x - 3) % 5 == 0) for x in range(-300, 300))
# variante 1 : 15x ≡ 21 (36) ⇔ 5x ≡ 7 (12) ⇔ x ≡ 11 (12) : 3 solutions.
assert gcd(15, 36) == 3 and lin(15, 21, 36) == [11, 23, 35] and (5 * 5) % 12 == 1 and 35 % 12 == 11
assert all(((15 * x - 21) % 36 == 0) == ((5 * x - 7) % 12 == 0) for x in range(-300, 300))
# variante 2 : 14x ≡ 20 (48) ⇔ 7x ≡ 10 (24) ⇔ x ≡ 22 (24) : 2 solutions.
assert gcd(14, 48) == 2 and lin(14, 20, 48) == [22, 46] and (7 * 7) % 24 == 1 and 70 % 24 == 22
assert all(((14 * x - 20) % 48 == 0) == ((7 * x - 10) % 24 == 0) for x in range(-300, 300))

print("enrichi_arith_l1_congruences : tout est vérifié.")
