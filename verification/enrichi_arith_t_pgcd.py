"""Vérifications calculatoires de l'enrichissement du chapitre « arith-t-pgcd » : questions de QCM
(content/qcm/arith-t-pgcd.json) et exercices ajoutés à content/exercices/arithmetique/terminale-arith-t-pgcd.json
(arith-t-bezout-ex5, arith-t-pgcd-euclide-ex6, arith-t-gauss-ex5, variantes comprises).

Chaque assert cite l'identifiant de la question ou de l'exercice concerné. Voir CLAUDE.md, « Contrôle de l'exactitude ».
Lancer : python3 verification/enrichi_arith_t_pgcd.py (doit finir sans erreur).
"""

from math import gcd


def euclide(a, b):
    etapes = []
    while b:
        q, r = divmod(a, b)
        etapes.append((a, b, q, r))
        a, b = b, r
    return a, etapes


# --- QCM -------------------------------------------------------------------------------------------------------
# arith-t-pgcd-q1 : 1071 ∧ 462 = 21 par les divisions 1071 = 462·2 + 147, 462 = 147·3 + 21, 147 = 21·7 ; leurres faux.
d, et = euclide(1071, 462)
assert d == 21 and et == [(1071, 462, 2, 147), (462, 147, 3, 21), (147, 21, 7, 0)]
assert 7 != d and 1071 % 42 != 0 and 462 % 63 == 21

# arith-t-pgcd-q2 : (2n+1) ∧ (n+5) prend exactement les valeurs 1, 3, 9 (diviseurs de 9), atteintes en n = 0, 1, 4.
valeurs = {gcd(2 * n + 1, n + 5) for n in range(0, 2000)}
assert valeurs == {1, 3, 9}
assert gcd(1, 5) == 1 and gcd(3, 6) == 3 and gcd(9, 9) == 9

# arith-t-pgcd-q3 (faux) : a = c = 1, b = 2 ; a∧b = c∧b = 1 mais (a+c)∧b = 2.
assert gcd(1, 2) == 1 and gcd(1 + 1, 2) == 2

# arith-t-pgcd-q4 (vrai) : si a∧b = 1 alors (a+b)∧ab = 1 (test exhaustif sur une plage, a, b non nuls).
for a in range(-40, 41):
    for b in range(-40, 41):
        if a and b and gcd(a, b) == 1:
            assert gcd(a + b, a * b) == 1

# arith-t-pgcd-q5 : seul (−2, 7) vérifie 17u + 5v = 1.
assert 17 * (-2) + 5 * 7 == 1
assert [17 * u + 5 * v for (u, v) in [(2, -7), (-3, 10), (1, -3)]] == [-1, -1, 2]

# arith-t-pgcd-q6 : solutions de 4x + 6y = 10 = {(1+3k, 1−2k)} ; les trois autres ensembles sont faux.
sol = {(x, y) for x in range(-60, 61) for y in range(-60, 61) if 4 * x + 6 * y == 10}
assert sol == {(1 + 3 * k, 1 - 2 * k) for k in range(-100, 101) if -60 <= 1 + 3 * k <= 60 and -60 <= 1 - 2 * k <= 60}
assert (4, -1) in sol and all(((x - 1) % 6 != 0) for (x, y) in [(4, -1)])  # oublié par {(1+6k, 1−4k)}
assert 4 * 4 + 6 * 3 != 10  # (4, 3) dans {(1+3k, 1+2k)}
assert 4 * (-1) + 6 * 4 == 20  # (−1, 4) dans {(1−2k, 1+3k)}

# arith-t-pgcd-q7 : 3 couples d'entiers naturels avec 3x + 5y = 47 : (14,1), (9,4), (4,7).
assert [(x, y) for x in range(48) for y in range(48) if 3 * x + 5 * y == 47] == [(4, 7), (9, 4), (14, 1)]

# arith-t-pgcd-q8 (faux) : 4 | 6·2 mais 4 ∤ 2 ; et 4 | 6n ⇔ 2 | n.
assert (6 * 2) % 4 == 0 and 2 % 4 != 0
assert all(((6 * n) % 4 == 0) == (n % 2 == 0) for n in range(-100, 100))

# --- arith-t-bezout-ex5 (chiffrement affine) ---------------------------------------------------------------------
ALPHA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def verifier_affine(a, b, u, v, ap, c, lettre_code, code_attendu, bonne, leurres):
    assert gcd(a, 26) == 1 and a * u + 26 * v == 1  # question 2 a)
    assert 0 <= ap <= 25 and (a * ap) % 26 == 1 and ap % 26 == u % 26  # question 2 b)
    assert c == (-ap * b) % 26  # question 3 : x ≡ a'y + c
    for x in range(26):  # question 3 : équivalence, testée sur toutes les classes
        for y in range(26):
            assert ((y - a * x - b) % 26 == 0) == ((x - ap * y - c) % 26 == 0)
    assert ALPHA[(a * ALPHA.index(lettre_code) + b) % 26] == code_attendu  # question 1
    assert ALPHA[(ap * 22 + c) % 26] == bonne  # question 4 : W = 22
    assert ALPHA[(a * ALPHA.index(bonne) + b) % 26] == "W"
    for l in leurres:  # leurres : lettres qui ne sont pas codées par W
        assert ALPHA[(a * ALPHA.index(l) + b) % 26] != "W"
    assert len({(a * x + b) % 26 for x in range(26)}) == 26  # question 5 : permutation


# exercice : clé (11, 8), (u, v) = (−7, 3), a' = 19, x ≡ 19y + 4, K → O, W ← G ; leurres C, Q, U.
verifier_affine(11, 8, -7, 3, 19, 4, "K", "O", "G", ["C", "Q", "U"])
assert euclide(26, 11)[1] == [(26, 11, 2, 4), (11, 4, 2, 3), (4, 3, 1, 1), (3, 1, 3, 0)]
assert 11 * 10 + 8 == 118 == 26 * 4 + 14 and 11 * 19 == 209 == 8 * 26 + 1 and -152 == -6 * 26 + 4
assert 209 % 26 == 1 and 44 % 26 == 18 and (-18) % 26 == 8 and 19 * 22 + 4 == 422 == 26 * 16 + 6 and 74 == 2 * 26 + 22
# Leurres : C = 19·22 mod 26 (oubli du + 4), Q = codage de W, U = 19·22 − 8.
assert (19 * 22) % 26 == 2 and (11 * 22 + 8) % 26 == 16 and (19 * 22 - 8) % 26 == 20
# variante 1 : clé (7, 3), (u, v) = (−11, 3), a' = 15, x ≡ 15y + 7, K → V, W ← Z ; leurres S, B, P.
verifier_affine(7, 3, -11, 3, 15, 7, "K", "V", "Z", ["S", "B", "P"])
assert euclide(26, 7)[1] == [(26, 7, 3, 5), (7, 5, 1, 2), (5, 2, 2, 1), (2, 1, 2, 0)]
assert 7 * 10 + 3 == 73 == 2 * 26 + 21 and 105 == 4 * 26 + 1 and -45 + 52 == 7 and -49 + 52 == 3 and 49 % 26 == 23
assert 15 * 22 + 7 == 337 == 12 * 26 + 25 and 7 * 25 + 3 == 178 == 6 * 26 + 22
assert (15 * 22) % 26 == 18 and (7 * 22 + 3) % 26 == 1 and (15 * 22 - 3) % 26 == 15
# variante 2 : clé (5, 12), (u, v) = (−5, 1), a' = 21, x ≡ 21y + 8, E → G, W ← C ; leurres U, S, I.
verifier_affine(5, 12, -5, 1, 21, 8, "E", "G", "C", ["U", "S", "I"])
assert euclide(26, 5)[1] == [(26, 5, 5, 1), (5, 1, 5, 0)]
assert 252 == 9 * 26 + 18 and -252 == -10 * 26 + 8 and 40 % 26 == 14 and -40 + 52 == 12
assert 21 * 22 + 8 == 470 == 18 * 26 + 2 and 5 * 2 + 12 == 22
assert (21 * 22) % 26 == 20 and (5 * 22 + 12) % 26 == 18 and (21 * 22 - 12) % 26 == 8

# --- arith-t-pgcd-euclide-ex6 (somme et PGCD donnés) ---------------------------------------------------------------


def couples(S, d):
    return [(a, S - a) for a in range(0, S + 1) if a <= S - a and gcd(a, S - a) == d]


# exercice : a + b = 168, a ∧ b = 12 → 3 couples ; leurres 7, 4, 6 faux.
assert couples(168, 12) == [(12, 156), (36, 132), (60, 108)]
assert [x for x in range(1, 8) if gcd(x, 14) == 1] == [1, 3, 5]
# variante 1 : a + b = 300, a ∧ b = 15 → 4 couples.
assert couples(300, 15) == [(15, 285), (45, 255), (105, 195), (135, 165)]
assert [x for x in range(1, 11) if gcd(x, 20) == 1] == [1, 3, 7, 9]
# variante 2 : a + b = 168, a ∧ b = 8 → 6 couples.
assert couples(168, 8) == [(8, 160), (16, 152), (32, 136), (40, 128), (64, 104), (80, 88)]
assert [x for x in range(1, 11) if gcd(x, 21) == 1] == [1, 2, 4, 5, 8, 10]
# Le lemme de la question 2 : a' ∧ (m − a') = a' ∧ m pour 0 ≤ a' ≤ m.
assert all(gcd(x, m - x) == gcd(x, m) for m in range(1, 60) for x in range(0, m + 1))

# --- arith-t-gauss-ex5 (ax = by avec contrainte) -------------------------------------------------------------------
# exercice : 8x = 15y ⇔ (x, y) = (15k, 8k) ; 100 ≤ x + y ≤ 200 → 4 couples.
sol = [(x, y) for x in range(-400, 401) for y in range(-400, 401) if 8 * x == 15 * y]
assert sol == [(15 * k, 8 * k) for k in range(-26, 27)]
assert [(x, y) for (x, y) in sol if 100 <= x + y <= 200] == [(75, 40), (90, 48), (105, 56), (120, 64)]
assert 23 * 4 == 92 and 23 * 9 == 207
# variante 1 : 7x = 12y, 100 ≤ x + y ≤ 300 → 10 couples (k de 6 à 15).
sol = [(x, y) for x in range(-400, 401) for y in range(-400, 401) if 7 * x == 12 * y]
assert sol == [(12 * k, 7 * k) for k in range(-33, 34)]
assert [(x, y) for (x, y) in sol if 100 <= x + y <= 300] == [(12 * k, 7 * k) for k in range(6, 16)]
assert 19 * 5 == 95 and 19 * 16 == 304
# variante 2 : 10x = 21y, 50 ≤ x − y ≤ 150 → 9 couples (k de 5 à 13).
sol = [(x, y) for x in range(-400, 401) for y in range(-400, 401) if 10 * x == 21 * y]
assert sol == [(21 * k, 10 * k) for k in range(-19, 20)]
assert [(x, y) for (x, y) in sol if 50 <= x - y <= 150] == [(21 * k, 10 * k) for k in range(5, 14)]
assert 11 * 4 == 44 and 11 * 14 == 154

print("enrichi_arith_t_pgcd : tout est vérifié.")
