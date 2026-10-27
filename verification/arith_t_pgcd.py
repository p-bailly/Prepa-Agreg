"""Vérifications calculatoires du chapitre « arith-t-pgcd » (cartes et exercices du fichier
content/{cartes,exercices}/arithmetique/terminale-arith-t-pgcd.json : algorithme d'Euclide, équations ax + by = c).

Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise) concerné.
Lancer : python3 verification/arith_t_pgcd.py (doit finir sans erreur).
"""

from itertools import product
from math import gcd

from sympy import expand, fibonacci, symbols

k = symbols("k")


def euclide(a, b):
    """Liste des divisions (a, b, q, r) de l'algorithme d'Euclide et PGCD."""
    etapes = []
    while b:
        q, r = divmod(a, b)
        etapes.append((a, b, q, r))
        a, b = b, r
    return a, etapes


def nb_divisions(a, b):
    return len(euclide(a, b)[1])


# ---------------------------------------------------------------------------
# arith-t-pgcd-euclide (carte)
# ---------------------------------------------------------------------------

# Lemme (quotient quelconque) : a = bq + r => a ∧ b = b ∧ r (exhaustif ; gcd de Python : gcd(x, 0) = |x|).
for a, b, q in product(range(-40, 41), range(1, 25), range(-5, 6)):
    r = a - b * q
    assert gcd(a, b) == gcd(b, r)
    # mêmes diviseurs communs
    assert {d for d in range(1, 50) if a % d == 0 and b % d == 0} == {d for d in range(1, 50) if b % d == 0 and r % d == 0}
# Algorithme : le dernier reste non nul est le PGCD, en au plus b divisions (a ∈ Z, b ≥ 1).
for a, b in product(range(-60, 61), range(1, 40)):
    g, et = euclide(a, b)
    assert g == gcd(a, b) and len(et) <= b + 1
    restes = [b] + [r for (_, _, _, r) in et]
    assert all(restes[i] > restes[i + 1] >= 0 for i in range(len(restes) - 1))
# Exemple : 2340 ∧ 1428 = 12, divisions citées, 195 ∧ 119 = 1.
g, et = euclide(2340, 1428)
assert g == 12 and et == [(2340, 1428, 1, 912), (1428, 912, 1, 516), (912, 516, 1, 396), (516, 396, 1, 120),
                          (396, 120, 3, 36), (120, 36, 3, 12), (36, 12, 3, 0)]
assert 2340 == 12 * 195 and 1428 == 12 * 119 and gcd(195, 119) == 1
# Cours, théorème de Lamé : nombre de divisions ≤ 5 × (nombre de chiffres de b) (vérifié pour b < 3000, a ≤ 3b + 50,
# et sur le pire cas de Fibonacci jusqu'à b < 10^12).
for b in range(1, 3000, 7):
    for a in range(b + 1, 3 * b + 50, 13):
        assert nb_divisions(a, b) <= 5 * len(str(b))
for m in range(3, 60):
    Fa, Fb = int(fibonacci(m + 1)), int(fibonacci(m))
    if Fb < 10 ** 12:
        assert nb_divisions(Fa, Fb) <= 5 * len(str(Fb))
# Cours / contre-exemple : (3 + 1) ∧ (3 - 1) = 2 ≠ 3 ∧ 1 = 1 (affirmation 3 aussi).
assert gcd(4, 2) == 2 and gcd(3, 1) == 1
# Contre-exemple « quotients quelconques » : 10 = 4 × 0 + 10, puis 4 = 10 × 0 + 4 : retour au couple (10, 4).
assert 10 == 4 * 0 + 10 and 4 == 10 * 0 + 4
# Contre-exemple « dernier reste » : 36 = 12 × 3 + 0 et 36 ∧ 12 = 12.
assert 36 == 12 * 3 + 0 and gcd(36, 12) == 12
# Contre-exemple / leurre {{1}} « q ∧ r » : 6 = 2 × 3 + 0, 3 ∧ 0 = 3 ≠ 6 ∧ 2 = 2.
assert divmod(6, 2) == (3, 0) and gcd(3, 0) == 3 and gcd(6, 2) == 2
# Leurre {{1}} « a ∧ b = r » : 12 = 5 × 2 + 2, mais 12 ∧ 5 = 1 ≠ 2.
assert divmod(12, 5) == (2, 2) and gcd(12, 5) == 1
# Leurre {{1}} « a ∧ b = a ∧ q » : 6 = 2 × 3, 6 ∧ 3 = 3 ≠ 2.
assert gcd(6, 3) == 3 != gcd(6, 2)
# Leurre {{3}} « le dernier quotient » : (12, 5) : divisions 12 = 5·2 + 2, 5 = 2·2 + 1, 2 = 1·2 + 0 : dernier quotient 2 ≠ 1.
assert euclide(12, 5)[1][-1][2] == 2 and gcd(12, 5) == 1
# Leurre {{3}} « le premier reste non nul » (r_1 = b) : (12, 5) : 5 ≠ 1 ; « le dernier reste » vaut 0.
assert euclide(12, 5)[1][-1][3] == 0
# Leurres {{2}} « strictement croissante », « constante », « périodique » : pour (12, 5), restes 5, 2, 1, 0.
assert [5] + [r for (_, _, _, r) in euclide(12, 5)[1]] == [5, 2, 1, 0]
# Affirmation 1 (vraie) : a ∧ b = b ∧ (a - kb) (couverte par le lemme ci-dessus). Affirmation 2 (fausse) : 30 = 12 + 18.
assert 30 == 12 * 1 + 18 and 18 >= 12 and gcd(12, 18) == 6 == gcd(30, 12)
# Affirmation 4 (vraie) : si 0 ≤ a < b, la première division échange les nombres.
assert all(euclide(a, b)[1][0] == (a, b, 0, a) for a, b in product(range(0, 20), range(1, 20)) if a < b)
# Pièges : « a ≥ b » et « a ∧ b = 1 » superflus (exemple (5, 12) et (2340, 1428)).
assert euclide(5, 12)[0] == 1 and euclide(2340, 1428)[0] == 12

# ---------------------------------------------------------------------------
# arith-t-pgcd-diophantiennes (carte)
# ---------------------------------------------------------------------------


def sols(a, b, c, borne=60):
    return {(x, y) for x in range(-borne, borne + 1) for y in range(-borne, borne + 1) if a * x + b * y == c}


# Énoncé : existence <=> d | c, et description des solutions (exhaustif sur un petit domaine).
for a, b in product([x for x in range(-9, 10) if x], repeat=2):
    d = gcd(a, b)
    a1, b1 = a // d, b // d
    assert gcd(a1, b1) == 1
    for c in range(-12, 13):
        S = sols(a, b, c, 40)
        if c % d:
            assert not S
        else:
            assert S
            x0, y0 = min(S)
            famille = {(x0 + t * b1, y0 - t * a1) for t in range(-200, 201)}
            assert S == {p for p in famille if abs(p[0]) <= 40 and abs(p[1]) <= 40}
# Exemple : 12x + 18y = 30 <=> 2x + 3y = 5, (1, 1), solutions (1 + 3k, 1 - 2k) ; (4, -1), (-2, 3) ; 12x + 18y = 20 : aucune.
assert gcd(12, 18) == 6 and 30 % 6 == 0 and 2 + 3 == 5
assert sols(12, 18, 30) == {(1 + 3 * t, 1 - 2 * t) for t in range(-30, 30)} & {(x, y) for x in range(-60, 61) for y in range(-60, 61)}
assert 12 * 4 + 18 * (-1) == 30 and 12 * (-2) + 18 * 3 == 30 and not sols(12, 18, 20)
# Contre-exemple « sans d | c » : 6x + 10y = 7 n'a pas de solution.
assert gcd(6, 10) == 2 and not sols(6, 10, 7)
# Contre-exemple « oubli de la division par d » (et leurre {{2}} (x0 + kb, y0 - ka), affirmation 3) :
# 2x + 4y = 6 : (3, 0) solution, pas de la forme (1 + 4k, 1 - 2k) ; bonne famille (1 + 2k, 1 - k).
assert 2 * 3 + 4 * 0 == 6 and all((1 + 4 * t, 1 - 2 * t) != (3, 0) for t in range(-50, 50))
assert sols(2, 4, 6) == {(1 + 2 * t, 1 - t) for t in range(-100, 100)} & {(x, y) for x in range(-60, 61) for y in range(-60, 61)}
# Leurre {{2}} (x0 + ka', y0 - kb') : (1 + 1, 1 - 2) = (2, -1) n'est pas solution de 2x + 4y = 6.
assert 2 * 2 + 4 * (-1) != 6
# Leurre {{2}} (x0 + kb', y0 + ka') : (1 + 2, 1 + 1) = (3, 2) n'est pas solution.
assert 2 * 3 + 4 * 2 != 6
# Contre-exemple « solutions positives » / affirmation 4 : 2x + 3y = 1 sans solution dans N².
assert not any(2 * x + 3 * y == 1 for x in range(0, 10) for y in range(0, 10)) and gcd(2, 3) == 1
# Leurre {{1}} « a ∧ b = 1 » : 2x + 4y = 6 a des solutions avec 2 ∧ 4 = 2.
assert sols(2, 4, 6) and gcd(2, 4) == 2
# Leurre {{1}} « c | d » : 2x + 4y = 1 : 1 | 2 mais aucune solution.
assert 2 % 1 == 0 and not sols(2, 4, 1)
# Leurre {{1}} « ab | c » : 2x + 3y = 1 a des solutions (dans Z²) alors que 6 ∤ 1.
assert sols(2, 3, 1) and 1 % 6 != 0
# Pièges : « c ≠ 0 » superflu (c = 0 : solution (0, 0)) ; « c multiple de ab » non nécessaire (2x + 3y = 1).
assert (0, 0) in sols(2, 4, 0)
# Affirmation 1 (vraie) : 6x + 15y = 4 sans solution. Affirmation 2 (fausse) : 2x + 3y = 5, (1, 1), 2 ∧ 3 = 1 ≠ 5.
assert not sols(6, 15, 4) and 2 * 1 + 3 * 1 == 5 and gcd(2, 3) == 1

# ---------------------------------------------------------------------------
# arith-t-pgcd-euclide-ex1 (et variantes)
# ---------------------------------------------------------------------------
g, et = euclide(4199, 1309)  # exercice
assert g == 17 and et == [(4199, 1309, 3, 272), (1309, 272, 4, 221), (272, 221, 1, 51), (221, 51, 4, 17), (51, 17, 3, 0)]
assert 4199 == 17 * 247 and 1309 == 17 * 77 and gcd(247, 77) == 1 and all(x != 17 for x in (51, 221, 1))
g, et = euclide(2457, 1386)  # variante 1
assert g == 63 and et == [(2457, 1386, 1, 1071), (1386, 1071, 1, 315), (1071, 315, 3, 126), (315, 126, 2, 63), (126, 63, 2, 0)]
assert 2457 == 63 * 39 and 1386 == 63 * 22 and gcd(39, 22) == 1 and all(x != 63 for x in (126, 315, 21))
g, et = euclide(1890, 1155)  # variante 2
assert g == 105 and et == [(1890, 1155, 1, 735), (1155, 735, 1, 420), (735, 420, 1, 315), (420, 315, 1, 105), (315, 105, 3, 0)]
assert 1890 == 105 * 18 and 1155 == 105 * 11 and gcd(18, 11) == 1 and all(x != 105 for x in (315, 35, 15))

# ---------------------------------------------------------------------------
# arith-t-pgcd-euclide-ex2 (et variantes) : pavage par des carrés
# ---------------------------------------------------------------------------
for (L, l, c, N, leurres, et_attendu) in [
    (1632, 1020, 204, 40, [(102, 160), (204, 13), (612, 4)],
     [(1632, 1020, 1, 612), (1020, 612, 1, 408), (612, 408, 1, 204), (408, 204, 2, 0)]),          # exercice
    (1078, 924, 154, 42, [(77, 168), (154, 13), (308, 10)],
     [(1078, 924, 1, 154), (924, 154, 6, 0)]),                                                      # variante 1
    (2405, 1443, 481, 15, [(37, 2535), (481, 8), (962, 4)],
     [(2405, 1443, 1, 962), (1443, 962, 1, 481), (962, 481, 2, 0)]),                                # variante 2
]:
    assert euclide(L, l) == (c, et_attendu)
    assert max(d for d in range(1, l + 1) if L % d == 0 and l % d == 0) == c
    assert (L // c) * (l // c) == N
    # Leurres faux : soit le côté n'est pas maximal (ou ne divise pas), soit le nombre de carreaux est faux.
    for (cl, Nl) in leurres:
        if cl == c:
            assert Nl != N
        elif L % cl == 0 and l % cl == 0:
            assert cl < c
        else:
            assert L % cl != 0 or l % cl != 0
assert (1632 // 204, 1020 // 204) == (8, 5) and (1078 // 154, 924 // 154) == (7, 6) and (2405 // 481, 1443 // 481) == (5, 3)

# ---------------------------------------------------------------------------
# arith-t-pgcd-euclide-ex3 (et variantes) : a = n² + k0, b = n + s
# ---------------------------------------------------------------------------
nn = symbols("n")
for (s, k0, p, res, n_test, A, B, qa, qb) in [(3, 2, 11, 8, 30, 902, 33, 82, 3),     # exercice
                                              (2, 3, 7, 5, 12, 147, 14, 21, 2),      # variante 1
                                              (4, 1, 17, 13, 30, 901, 34, 53, 2)]:   # variante 2
    assert expand((nn + s) * (nn - s) + p - (nn ** 2 + k0)) == 0
    for m in range(0, 400):
        g = gcd(m * m + k0, m + s)
        assert g == (p if m % p == res else 1)
    assert (-s) % p == res
    assert n_test % p == res and n_test ** 2 + k0 == A == p * qa and n_test + s == B == p * qb and gcd(qa, qb) == 1
    assert gcd(A, B) == p
    # Leurres : n ≡ s (mod p) et n ≡ 0 (mod p) ne donnent pas le PGCD p ; « aucun » est faux.
    assert gcd(s ** 2 + k0, s + s) == 1 and gcd(k0, s) == 1

# ---------------------------------------------------------------------------
# arith-t-pgcd-euclide-ex4 (et variantes) : (2^a - 1) ∧ (2^b - 1) = 2^(a ∧ b) - 1
# ---------------------------------------------------------------------------
for a, b in product(range(0, 40), range(1, 40)):
    q, r = divmod(a, b)
    Ma, Mb, Mbq, Mr = 2 ** a - 1, 2 ** b - 1, 2 ** (b * q) - 1, 2 ** r - 1
    assert Mbq % Mb == 0 and Ma == 2 ** r * Mbq + Mr and gcd(Ma, Mb) == gcd(Mb, Mr) == 2 ** gcd(a, b) - 1
for (a4, b4, g4, val, et4, a5, b5, et5, leurres) in [
    (91, 65, 13, 8191, [(91, 65, 1, 26), (65, 26, 2, 13), (26, 13, 2, 0)], 91, 64,
     [(91, 64, 1, 27), (64, 27, 2, 10), (27, 10, 2, 7), (10, 7, 1, 3), (7, 3, 2, 1), (3, 1, 3, 0)],
     [2 ** 13, 13, 2 ** 7 - 1]),                                                              # exercice
    (84, 60, 12, 4095, [(84, 60, 1, 24), (60, 24, 2, 12), (24, 12, 2, 0)], 84, 55,
     [(84, 55, 1, 29), (55, 29, 1, 26), (29, 26, 1, 3), (26, 3, 8, 2), (3, 2, 1, 1), (2, 1, 2, 0)],
     [2 ** 12, 12, 2 ** 24 - 1]),                                                             # variante 1
    (100, 75, 25, 33554431, [(100, 75, 1, 25), (75, 25, 3, 0)], 100, 63,
     [(100, 63, 1, 37), (63, 37, 1, 26), (37, 26, 1, 11), (26, 11, 2, 4), (11, 4, 2, 3), (4, 3, 1, 1), (3, 1, 3, 0)],
     [2 ** 25, 25, 2 ** 5 - 1]),                                                              # variante 2
]:
    assert euclide(a4, b4) == (g4, et4) and gcd(2 ** a4 - 1, 2 ** b4 - 1) == 2 ** g4 - 1 == val
    assert euclide(a5, b5) == (1, et5) and gcd(2 ** a5 - 1, 2 ** b5 - 1) == 1
    assert all(x != val for x in leurres)

# ---------------------------------------------------------------------------
# arith-t-pgcd-euclide-ex5 (et variantes) : nombre de divisions (Lamé)
# ---------------------------------------------------------------------------
F = [0] + [int(fibonacci(i)) for i in range(1, 30)]
assert F[1:18] == [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597]
# Questions 1 à 3 : q_k ≥ 1, q_N ≥ 2, b ≥ F_{N+1}, a ≥ F_{N+2} (exhaustif pour a > b ≥ 1, a ≤ 400).
for a in range(2, 401):
    for b in range(1, a):
        g, et = euclide(a, b)
        N = len(et)
        assert all(q >= 1 for (_, _, q, _) in et) and et[-1][2] >= 2
        assert b >= F[N + 1] and a >= F[N + 2]
# Question 4 : (F_{N+2}, F_{N+1}) demande exactement N divisions.
for N in range(1, 25):
    assert nb_divisions(F[N + 2], F[N + 1]) == N
# Question 5 : maximum pour 1 ≤ b ≤ B (exercice B = 100 : 10 ; variante 1 B = 50 : 8 ; variante 2 B = 1000 : 15).
for (B, maxi, couple, leurres) in [(100, 10, (144, 89), [11, 12, 9]), (50, 8, (55, 34), [9, 10, 7]),
                                    (1000, 15, (1597, 987), [16, 17, 14])]:
    assert nb_divisions(*couple) == maxi and couple[1] <= B
    # majoration : F_{N+1} ≤ B impose N ≤ maxi
    assert F[maxi + 1] <= B < F[maxi + 2]
    # contrôle par force brute (a parcourt un large intervalle ; la borne de Lamé vaut pour tout a)
    pas = 1 if B <= 100 else 37
    assert max(nb_divisions(a, b) for b in range(1, B + 1, pas if B > 100 else 1)
               for a in range(b + 1, 3 * B + 2 * b, 1 if B <= 100 else 29)) <= maxi
    assert all(x != maxi for x in leurres)

# ---------------------------------------------------------------------------
# arith-t-pgcd-diophantiennes-ex1 (et variantes) : quelle équation a des solutions ?
# ---------------------------------------------------------------------------
for (eqs, bonne, sol) in [([(6, 15, 20), (14, 21, 30), (10, 25, 35), (12, 18, 15)], 3, (1, 1)),   # exercice
                          ([(9, 12, 20), (8, 20, 26), (21, 28, 30), (15, 35, 45)], 4, (3, 0)),    # variante 1
                          ([(4, 10, 7), (18, 27, 36), (12, 20, 30), (15, 21, 10)], 2, (-1, 2))]:  # variante 2
    resolubles = [i + 1 for i, (a, b, c) in enumerate(eqs) if c % gcd(a, b) == 0]
    assert resolubles == [bonne]
    a, b, c = eqs[bonne - 1]
    assert a * sol[0] + b * sol[1] == c
    for i, (a, b, c) in enumerate(eqs):
        if i + 1 != bonne:
            assert not sols(a, b, c, 30)
# PGCD cités : 6 ∧ 15 = 3, 14 ∧ 21 = 7, 10 ∧ 25 = 5, 12 ∧ 18 = 6 ; 9 ∧ 12 = 3, 8 ∧ 20 = 4, 21 ∧ 28 = 7, 15 ∧ 35 = 5 ;
# 4 ∧ 10 = 2, 18 ∧ 27 = 9, 12 ∧ 20 = 4, 15 ∧ 21 = 3. Équations réduites : 2x + 5y = 7, 3x + 7y = 9, 2x + 3y = 4.
assert [gcd(*e) for e in [(6, 15), (14, 21), (10, 25), (12, 18), (9, 12), (8, 20), (21, 28), (15, 35), (4, 10), (18, 27), (12, 20), (15, 21)]] \
    == [3, 7, 5, 6, 3, 4, 7, 5, 2, 9, 4, 3]
assert 2 + 5 == 7 and 3 * 3 == 9 and 2 * (-1) + 3 * 2 == 4

# ---------------------------------------------------------------------------
# arith-t-pgcd-diophantiennes-ex2 (et variantes) : résolution complète
# ---------------------------------------------------------------------------


def famille(x0, kx, y0, ky, borne=60):
    return {(x0 + kx * t, y0 + ky * t) for t in range(-300, 301)} & {(x, y) for x in range(-borne, borne + 1) for y in range(-borne, borne + 1)}


for (a, b, c, bonne, leurres, temoin) in [
    (15, 6, 21, (1, 2, 1, -5), [(1, 6, 1, -15), (1, 2, 1, 5), (1, 5, 1, -2)], (3, -4)),            # exercice
    (20, 12, 28, (2, 3, -1, -5), [(2, 12, -1, -20), (2, 3, -1, 5), (2, 5, -1, -3)], (5, -6)),     # variante 1
    (14, -21, 35, (1, 3, -1, 2), [(1, 21, -1, 14), (1, 3, -1, -2), (1, 2, -1, 3)], (4, 1)),       # variante 2
]:
    S = sols(a, b, c)
    assert S == famille(*bonne)
    for L in leurres:
        assert famille(*L) != S
    assert temoin in S and temoin not in famille(*leurres[0])
    # la famille non réduite est formée de solutions
    assert famille(*leurres[0]) <= S
# Réductions et solutions particulières : 5x + 2y = 7 (1, 1) ; 5x + 3y = 7 (2, -1) ; 2x - 3y = 5 (1, -1).
assert 5 + 2 == 7 and 10 - 3 == 7 and 2 + 3 == 5
assert expand(5 * (1 + 2 * k) + 2 * (1 - 5 * k)) == 7 and expand(5 * (2 + 3 * k) + 3 * (-1 - 5 * k)) == 7
assert expand(2 * (1 + 3 * k) - 3 * (-1 + 2 * k)) == 5

# ---------------------------------------------------------------------------
# arith-t-pgcd-diophantiennes-ex3 (et variantes) : billetterie
# ---------------------------------------------------------------------------
for (p1, p2, S, u, v, kmin, kmax, compo, leurres) in [
    (9, 14, 300, -3, 2, 65, 66, [(10, 15), (24, 6)], [1, 3]),                 # exercice
    (7, 12, 250, -5, 3, 105, 107, [(10, 15), (22, 8), (34, 1)], [2, 4]),       # variante 1
    (11, 15, 500, -4, 3, 134, 136, [(10, 26), (25, 15), (40, 4)], [2, 4]),     # variante 2
]:
    assert p1 * u + p2 * v == 1 and gcd(p1, p2) == 1
    x0, y0 = u * S, v * S
    assert p1 * x0 + p2 * y0 == S
    assert expand(p1 * (x0 + p2 * k) + p2 * (y0 - p1 * k)) == S
    ks = [t for t in range(-10000, 10000) if x0 + p2 * t >= 0 and y0 - p1 * t >= 0]
    assert ks == list(range(kmin, kmax + 1))
    assert [(x0 + p2 * t, y0 - p1 * t) for t in ks] == compo
    assert [(x, y) for x in range(0, S // p1 + 1) for y in range(0, S // p2 + 1) if p1 * x + p2 * y == S] == compo
    assert all(len(compo) != L for L in leurres)  # « une infinité » est aussi faux (ensemble fini)
# Approximations citées : 900/14 ≈ 64,3 ; 600/9 ≈ 66,7 ; 1250/12 ≈ 104,2 ; 750/7 ≈ 107,1 ; 2000/15 ≈ 133,3 ; 1500/11 ≈ 136,4.
assert [round(x, 1) for x in (900 / 14, 600 / 9, 1250 / 12, 750 / 7, 2000 / 15, 1500 / 11)] == [64.3, 66.7, 104.2, 107.1, 133.3, 136.4]
# Remontées de l'algorithme d'Euclide citées.
assert euclide(14, 9)[1] == [(14, 9, 1, 5), (9, 5, 1, 4), (5, 4, 1, 1), (4, 1, 4, 0)] and 2 * 14 - 3 * 9 == 1
assert euclide(12, 7)[1][:3] == [(12, 7, 1, 5), (7, 5, 1, 2), (5, 2, 2, 1)] and 3 * 12 - 5 * 7 == 1
assert euclide(15, 11)[1][:3] == [(15, 11, 1, 4), (11, 4, 2, 3), (4, 3, 1, 1)] and 3 * 15 - 4 * 11 == 1

# ---------------------------------------------------------------------------
# arith-t-pgcd-diophantiennes-ex4 (et variantes) : points entiers d'une droite
# ---------------------------------------------------------------------------
for (a, b, c, d, X, x0, y0, nb, abscisses_sans_red, c4) in [
    (21, 15, 9, 3, 50, 4, 5, 10, [4, 19, 34, 49], 10),        # exercice
    (26, 10, 8, 2, 60, 3, 7, 12, [3, 13, 23, 33, 43, 53], 7),  # variante 1
    (18, 24, 30, 6, 45, 3, 1, 11, [3, 27], 20),                # variante 2
]:
    assert gcd(a, b) == d and c % d == 0
    a1, b1, c1 = a // d, b // d, c // d
    assert a1 * x0 - b1 * y0 == c1
    S = sols(a, -b, c, 200)
    assert S == famille(x0, b1, y0, a1, 200)
    pts = sorted(x for (x, y) in S if 0 <= x <= X)
    assert len(pts) == nb and pts == list(range(x0, x0 + b1 * nb, b1))
    assert sorted(x for x in range(0, X + 1) if (x - x0) % b == 0 and x >= x0) == abscisses_sans_red
    # leurres : nombre sans réduction, nb + 1, nb - 1
    assert all(L != nb for L in (len(abscisses_sans_red), nb + 1, nb - 1))
    # question 4 : aucune solution
    assert c4 % d != 0 and not sols(a, -b, c4, 100)

# ---------------------------------------------------------------------------
# arith-t-pgcd-diophantiennes-ex5 (et variantes) : montants payables
# ---------------------------------------------------------------------------
for (p, q, u, v, inv, F_, nb, liste, termes, leurres) in [
    (5, 8, -3, 2, 5, 27, 14, [1, 2, 3, 4, 6, 7, 9, 11, 12, 14, 17, 19, 22, 27], [0, 0, 1, 1, 2, 3, 3, 4], [27, 13, 20]),  # exercice
    (3, 7, -2, 1, 5, 11, 6, [1, 2, 4, 5, 8, 11], [0, 0, 0, 1, 1, 2, 2], [11, 10, 5]),                                      # variante 1
    (4, 9, -2, 1, 7, 23, 12, [1, 2, 3, 5, 6, 7, 10, 11, 14, 15, 19, 23], [0, 0, 0, 1, 1, 2, 2, 3, 3], [23, 13, 18]),      # variante 2
]:
    assert p * u + q * v == 1 and (p * inv) % q == 1

    def payable(m):
        return any((m - p * x) >= 0 and (m - p * x) % q == 0 for x in range(0, m // p + 1))

    def xn(m):
        return (u * m) % q

    for m in range(-50, 300):
        x = xn(m)
        assert 0 <= x <= q - 1 and (m - p * x) % q == 0
        assert [t for t in range(q) if (m - p * t) % q == 0] == [x]  # unicité
        if m >= 0:
            assert payable(m) == ((m - p * x) // q >= 0)
    non_payables = [m for m in range(0, 400) if not payable(m)]
    assert non_payables == liste and len(liste) == nb and max(liste) == F_ == p * q - p - q
    assert all(payable(m) for m in range(F_ + 1, 1000))
    assert [p * x // q for x in range(q)] == termes and sum(termes) == nb
    assert F_ % q * inv % q == q - 1 and (F_ - p * (q - 1)) // q == -1 and (F_ - p * (q - 1)) % q == 0
    assert (F_ + 1) - p * (q - 1) > -q
    assert all(L != nb for L in leurres)

print("arith_t_pgcd : toutes les vérifications passent.")
