"""Vérifications calculatoires du chapitre « arith-t-divisibilite » (cartes et exercices du fichier
content/{cartes,exercices}/arithmetique/terminale-arith-t-divisibilite.json).

Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise) concerné.
Lancer : python3 verification/arith_t_divisibilite.py (doit finir sans erreur).
"""

from fractions import Fraction
from itertools import product

from sympy import expand, symbols

n = symbols("n")


def divise(d, a):
    """d | a dans Z (0 ne divise que 0)."""
    return a == 0 if d == 0 else a % d == 0


def congru(a, b, m):
    return (a - b) % m == 0


# ---------------------------------------------------------------------------
# arith-t-divisibilite-combinaisons (carte)
# ---------------------------------------------------------------------------

# Énoncé : si d | a et d | b, alors d | au + bv (vérification exhaustive sur un petit domaine).
for d, a, b, u, v in product(range(-6, 7), range(-8, 9), range(-8, 9), range(-3, 4), range(-3, 4)):
    if divise(d, a) and divise(d, b):
        assert divise(d, a * u + b * v)

# Exemple : n - 1 | n + 5  <=>  n ∈ {-5, -2, -1, 0, 2, 3, 4, 7}.
assert {k for k in range(-200, 201) if divise(k - 1, k + 5)} == {-5, -2, -1, 0, 2, 3, 4, 7}
assert not divise(0, 6)  # cas n = 1 exclu

# Contre-exemple « coefficients rationnels » (et leurres {{2}} « Q », « R », « Q avec u + v entier ») :
# 2 | 2, 2 | 4, u = v = 1/2 donne 3, non divisible par 2 ; u + v = 1 est entier.
x = 2 * Fraction(1, 2) + 4 * Fraction(1, 2)
assert x == 3 and divise(2, 2) and divise(2, 4) and not divise(2, 3)
assert Fraction(1, 2) + Fraction(1, 2) == 1
# Contre-exemple « une seule des deux divisibilités » (piège, leurre {{1}} « a ou b ») : 2 | 2, 2 ∤ 3, 2 ∤ 5.
assert divise(2, 2) and not divise(2, 3) and not divise(2, 5)
# Contre-exemple « diviser le produit » (piège, leurre {{1}} « ab ») : 6 | 2 × 3 mais 6 ∤ 5.
assert divise(6, 2 * 3) and not divise(6, 2 + 3)
# Réciproque fausse / leurre {{1}} « a + b » / affirmation 1 : 2 | 1 + 1 mais 2 ∤ 1 (u = 1, v = 0 donne 1).
assert divise(2, 1 + 1) and not divise(2, 1)
# Leurre {{3}} « d² | au + bv » : d = 2, a = b = 2, u = 1, v = 0 : 4 ∤ 2.
assert divise(2, 2) and not divise(4, 2 * 1 + 2 * 0)
# Leurre {{3}} « au + bv | d » : u = v = 0, d = 2 : 0 ne divise pas 2.
assert not divise(0, 2)
# Leurre {{3}} « d | u + v » : a = b = 0, d = 2, u = 1, v = 0 : 2 ∤ 1.
assert divise(2, 0) and not divise(2, 1 + 0)
# Piège « d premier » : la propriété vaut pour d = 6 (non premier), par exemple 6 | 12, 6 | 18, 6 | 12·5 − 18·2.
assert divise(6, 12 * 5 - 18 * 2)
# Affirmation 2 (vraie) : d | a et d | a + b => d | b (exhaustif).
for d, a, b in product(range(1, 8), range(-10, 11), range(-10, 11)):
    if divise(d, a) and divise(d, a + b):
        assert divise(d, b)
# Affirmation 3 (fausse) : 6 | 2 × 3, 6 ∤ 2, 6 ∤ 3.
assert divise(6, 6) and not divise(6, 2) and not divise(6, 3)
# Affirmation 4 (fausse) : 2 | -2 et -2 | 2 mais 2 ≠ -2.
assert divise(2, -2) and divise(-2, 2) and 2 != -2

# ---------------------------------------------------------------------------
# arith-t-divisibilite-congruences (carte)
# ---------------------------------------------------------------------------

# Énoncé : n | a - b <=> même reste ; compatibilité somme, produit, puissances (exhaustif sur un petit domaine).
for m in range(1, 10):
    for a, b in product(range(-15, 16), repeat=2):
        assert congru(a, b, m) == (a % m == b % m)
    for a, b, c, d in product(range(-6, 7), repeat=4):
        if congru(a, b, m) and congru(c, d, m):
            assert congru(a + c, b + d, m) and congru(a * c, b * d, m)
            assert all(congru(a ** k, b ** k, m) for k in range(0, 6))

# Exemple : 7^2026 ≡ 9 (mod 10), via 49 ≡ -1, 7^4 ≡ 1 et 2026 = 4 × 506 + 2.
assert congru(49, -1, 10) and pow(7, 4, 10) == 1 and 2026 == 4 * 506 + 2 and pow(7, 2026, 10) == 9

# Contre-exemple « pas de simplification » / affirmation 1 : 2·3 ≡ 2·0 (mod 6), 2 ≢ 0, 3 ≢ 0 (mod 6).
assert congru(2 * 3, 2 * 0, 6) and not congru(2, 0, 6) and not congru(3, 0, 6)
# Contre-exemple « modules différents » : 1 ≡ 3 (mod 2), 0 ≡ 3 (mod 3), 1 et 6 non congrus mod 2, 3, 6.
assert congru(1, 3, 2) and congru(0, 3, 3) and all(not congru(1, 6, m) for m in (2, 3, 6))
# Contre-exemple et piège « modulo un diviseur m de n » : 0 ≡ 0 (mod 4), 0 ≡ 2 (mod 2), 0 ≢ 2 (mod 4).
assert congru(0, 0, 4) and congru(0, 2, 2) and not congru(0 + 0, 0 + 2, 4)
# Contre-exemple « exposants » / affirmation 2 : 3 ≡ 8 (mod 5), 2^3 ≡ 3, 2^8 ≡ 1 (mod 5).
assert congru(3, 8, 5) and pow(2, 3, 5) == 3 and pow(2, 8, 5) == 1
# Contre-exemple « diviseurs de zéro » / affirmation 4 : 2 × 3 ≡ 0 (mod 6).
assert congru(2 * 3, 0, 6) and not congru(2, 0, 6) and not congru(3, 0, 6)
# Affirmation 3 (vraie) : a ≡ b (mod n) et m | n => a ≡ b (mod m) (exhaustif).
for nn in range(1, 13):
    for m in [k for k in range(1, nn + 1) if nn % k == 0]:
        for a, b in product(range(-12, 13), repeat=2):
            if congru(a, b, nn):
                assert congru(a, b, m)


def contre_exemple_compat(relation, conclusion, borne=6, modules=range(1, 7)):
    """Cherche (m, a, b, c, d) avec a ~ b, c ~ d (pour « relation ») mais conclusion fausse."""
    for m in modules:
        for a, b, c, d in product(range(0, borne), repeat=4):
            if relation(a, b, m) and relation(c, d, m) and not conclusion(a, b, c, d, m):
                return (m, a, b, c, d)
    return None


# Leurre {{1}} « a - b divise n » : la compatibilité avec la somme devient fausse (n = 6 : 1 ~ 3, 0 ~ 3, 1 ≁ 6).
assert divise(1 - 3, 6) and divise(0 - 3, 6) and not divise(1 - 6, 6)
# Leurre {{1}} « n divise a + b » : la compatibilité avec le produit devient fausse (n = 3 : 1 ~ 2, 1 ~ 2, 1·1 ≁ 2·2).
assert divise(3, 1 + 2) and not divise(3, 1 * 1 + 2 * 2)
# Leurre {{1}} « n divise ab » : la compatibilité avec la somme devient fausse (n = 2 : 1 ~ 2, 2 ~ 1, 3 ≁ 3).
assert divise(2, 1 * 2) and divise(2, 2 * 1) and not divise(2, (1 + 2) * (2 + 1))
assert contre_exemple_compat(lambda a, b, m: divise(m, a * b), lambda a, b, c, d, m: divise(m, (a + c) * (b + d))) is not None
# Leurre {{2}} « ac ≡ bd (mod n²) » : n = 2, a = 0, b = 2, c = d = 1.
assert congru(0, 2, 2) and congru(1, 1, 2) and not congru(0 * 1, 2 * 1, 4)
# Leurre {{2}} « ac ≡ bd (mod 2n) » : n = 1, a = 0, b = 1, c = d = 1.
assert congru(0, 1, 1) and not congru(0 * 1, 1 * 1, 2)
# Leurre {{2}} « a^c ≡ b^d (mod n) » : n = 2, a = b = 2, c = 0, d = 2 : 1 ≢ 4 (mod 2).
assert congru(0, 2, 2) and not congru(2 ** 0, 2 ** 2, 2)
# Leurre {{3}} « a^k ≡ b^k (mod n^k) » : n = 2, a = 1, b = 3, k = 3 : 1 ≢ 27 (mod 8).
assert congru(1, 3, 2) and not congru(1 ** 3, 3 ** 3, 8)
# Leurre {{3}} « k^a ≡ k^b (mod n) » : n = 3, a = 0, b = 3, k = 2 : 1 ≢ 8 (mod 3).
assert congru(0, 3, 3) and not congru(2 ** 0, 2 ** 3, 3)
# Leurre {{3}} « a^k ≡ b (mod n) » : n = 5, a = b = 2, k = 2 : 4 ≢ 2 (mod 5).
assert not congru(2 ** 2, 2, 5)

# ---------------------------------------------------------------------------
# arith-t-divisibilite-combinaisons-ex1 (et variantes 1, 2)
# ---------------------------------------------------------------------------
for (m, a, b, c, k, comb) in [(7, 3, 2, 3, 5, 7 * (2 * n + 1)),      # exercice : 5(3n+2) - (n+3) = 14n + 7 = 7(2n+1)
                              (11, 4, 3, 9, 3, 11 * n),             # variante 1 : 3(4n+3) - (n+9) = 11n
                              (13, 5, 1, 8, 8, 13 * 3 * n)]:        # variante 2 : 8(5n+1) - (n+8) = 39n = 13 × 3n
    assert expand(k * (a * n + b) - (n + c) - comb) == 0
    assert all(divise(m, kk + c) for kk in range(-500, 500) if divise(m, a * kk + b))
# Coefficients intermédiaires de la correction : (3k - 1)n + (2k - 3) etc.
assert expand(5 * (3 * n + 2) - (n + 3)) == 14 * n + 7

# ---------------------------------------------------------------------------
# arith-t-divisibilite-combinaisons-ex2 (et variantes) : n + s | a n + B
# ---------------------------------------------------------------------------


def ensemble(cond, borne=400):
    return {k for k in range(-borne, borne + 1) if cond(k)}


for (s, a, B, c, bon, leurres) in [
    (3, 2, 13, 7, {-10, -4, -2, 4}, [{-2, 4}, {4}, {-10, -4, -3, -2, 4}]),                    # exercice
    (2, 3, 17, 11, {-13, -3, -1, 9}, [{-1, 9}, {9}, {-13, -3, -2, -1, 9}]),                    # variante 1
    (-4, 2, 1, 9, {-5, 1, 3, 5, 7, 13}, [{5, 7, 13}, {13}, {-5, 1, 3, 4, 5, 7, 13}]),          # variante 2
]:
    assert a * (n + s) + c == expand(a * n + B)  # aN + B = a(n + s) + c
    assert ensemble(lambda k: divise(k + s, a * k + B)) == bon
    assert all(L != bon for L in leurres)  # leurres faux

# ---------------------------------------------------------------------------
# arith-t-divisibilite-combinaisons-ex3 (et variantes) : n + s | n² + k0
# ---------------------------------------------------------------------------
for (s, k0, c, bon, leurres) in [
    (2, 5, 9, {-11, -5, -3, -1, 1, 7}, [{-1, 1, 7}, {1, 7}, {-11, -5, -3, -2, -1, 1, 7}]),                      # exercice
    (1, 7, 8, {-9, -5, -3, -2, 0, 1, 3, 7}, [{0, 1, 3, 7}, {1, 3, 7}, {-9, -5, -3, -2, -1, 0, 1, 3, 7}]),       # variante 1
    (-3, 1, 10, {-7, -2, 1, 2, 4, 5, 8, 13}, [{4, 5, 8, 13}, {5, 8, 13}, {-7, -2, 1, 2, 3, 4, 5, 8, 13}]),      # variante 2
]:
    assert expand((n + s) * (n - s) + c) == expand(n ** 2 + k0)
    assert ensemble(lambda k: divise(k + s, k * k + k0)) == bon
    assert all(L != bon for L in leurres)
# Vérifications numériques citées (exercice) : 9 | 54, 3 | 6, -9 | 126, -3 | 30.
assert 7 ** 2 + 5 == 54 and 1 + 5 == 6 and 121 + 5 == 126 and 25 + 5 == 30
# (variante 1) 8 | 56, 4 | 16, 2 | 8, -2 | 16, -4 | 32, -8 | 88.
assert [k * k + 7 for k in (7, 3, 1, -3, -5, -9)] == [56, 16, 8, 16, 32, 88]
# (variante 2) 10 | 170, 5 | 65, 2 | 26, -2 | 2, -5 | 5, -10 | 50.
assert [k * k + 1 for k in (13, 8, 5, 1, -2, -7)] == [170, 65, 26, 2, 5, 50]

# ---------------------------------------------------------------------------
# arith-t-divisibilite-combinaisons-ex4 (et variantes) : critères de divisibilité
# ---------------------------------------------------------------------------
A, Bb, N = symbols("a b N")
NN = 10 * A + Bb
# Exercice (7) : a - 2b = 21a - 2N, N = 10(a - 2b) + 21b ; faux critère a + 2b = 2N - 19a.
assert expand(21 * A - 2 * NN - (A - 2 * Bb)) == 0 and expand(10 * (A - 2 * Bb) + 21 * Bb - NN) == 0
assert expand(2 * NN - 19 * A - (A + 2 * Bb)) == 0 and 21 % 7 == 0 and 19 % 7 != 0
# Variante 1 (13) : a + 4b = 4N - 39a, N = 10(a + 4b) - 39b ; faux critère a - 4b = 41a - 4N.
assert expand(4 * NN - 39 * A - (A + 4 * Bb)) == 0 and expand(10 * (A + 4 * Bb) - 39 * Bb - NN) == 0
assert expand(41 * A - 4 * NN - (A - 4 * Bb)) == 0 and 39 % 13 == 0 and 41 % 13 != 0
# Variante 2 (17) : a - 5b = 51a - 5N, N = 10(a - 5b) + 51b ; faux critère a + 5b = 5N - 49a.
assert expand(51 * A - 5 * NN - (A - 5 * Bb)) == 0 and expand(10 * (A - 5 * Bb) + 51 * Bb - NN) == 0
assert expand(5 * NN - 49 * A - (A + 5 * Bb)) == 0 and 51 % 17 == 0 and 49 % 17 != 0


def etape(x, signe, coef):
    return x // 10 + signe * coef * (x % 10)


# Le critère est une équivalence (vérification exhaustive).
for x in range(10, 20000):
    assert (x % 7 == 0) == (etape(x, -1, 2) % 7 == 0)
    assert (x % 13 == 0) == (etape(x, +1, 4) % 13 == 0)
    assert (x % 17 == 0) == (etape(x, -1, 5) % 17 == 0)
# Exercice : 5488 -> 532 -> 49 ; 5478 -> 531 -> 51 ; 5486 -> 536 -> 41 ; 5498 -> 533 -> 47.
assert [etape(5488, -1, 2), etape(532, -1, 2)] == [532, 49] and 5488 == 7 * 784
assert [etape(5478, -1, 2), etape(531, -1, 2)] == [531, 51]
assert [etape(5486, -1, 2), etape(536, -1, 2)] == [536, 41]
assert [etape(5498, -1, 2), etape(533, -1, 2)] == [533, 47]
assert [x % 7 == 0 for x in (5478, 5486, 5488, 5498)] == [False, False, True, False]
assert 51 % 7 and 41 % 7 and 47 % 7
# Variante 1 : 5005 -> 520 -> 52 ; 4995 -> 519 -> 87 ; 5003 -> 512 -> 59 ; 5015 -> 521 -> 56.
assert [etape(5005, 1, 4), etape(520, 1, 4)] == [520, 52] and 5005 == 13 * 385 and 52 == 4 * 13
assert [etape(4995, 1, 4), etape(519, 1, 4)] == [519, 87]
assert [etape(5003, 1, 4), etape(512, 1, 4)] == [512, 59]
assert [etape(5015, 1, 4), etape(521, 1, 4)] == [521, 56]
assert [x % 13 == 0 for x in (4995, 5003, 5005, 5015)] == [False, False, True, False]
assert 87 % 13 and 59 % 13 and 56 % 13
# Variante 2 : 5457 -> 510 -> 51 ; 5447 -> 509 -> 5 ; 5455 -> 520 -> 52 ; 5467 -> 511 -> 46.
assert [etape(5457, -1, 5), etape(510, -1, 5)] == [510, 51] and 5457 == 17 * 321 and 51 == 3 * 17
assert [etape(5447, -1, 5), etape(509, -1, 5)] == [509, 5]
assert [etape(5455, -1, 5), etape(520, -1, 5)] == [520, 52]
assert [etape(5467, -1, 5), etape(511, -1, 5)] == [511, 46]
assert [x % 17 == 0 for x in (5447, 5455, 5457, 5467)] == [False, False, True, False]
assert 5 % 17 and 52 % 17 and 46 % 17
# Question 4 : faux critères testés sur 14 (1 + 8 = 9), 13 (1 - 12 = -11), 17 (1 + 35 = 36).
assert 14 % 7 == 0 and (1 + 2 * 4) == 9 and 9 % 7 != 0
assert 13 % 13 == 0 and (1 - 4 * 3) == -11 and -11 % 13 != 0
assert 17 % 17 == 0 and (1 + 5 * 7) == 36 and 36 % 17 != 0
# Exemple de l'énoncé : 2026 = 10 × 202 + 6.
assert 2026 == 10 * 202 + 6

# ---------------------------------------------------------------------------
# arith-t-divisibilite-congruences-ex1 (et variantes) : reste de 2a² - ab + 4
# ---------------------------------------------------------------------------
for (m, ra, rb, a2, ab, tot, res, leurres) in [(7, 5, 4, 4, 6, 6, 6, [2, 4, 5]),     # exercice
                                               (9, 7, 5, 4, 8, 4, 4, [0, 2, 7]),     # variante 1
                                               (11, 8, 6, 9, 4, 18, 7, [9, 4, 3])]:  # variante 2
    assert (ra * ra) % m == a2 and (ra * rb) % m == ab and 2 * a2 - ab + 4 == tot and tot % m == res
    # Indépendance du représentant : tous les a ≡ ra, b ≡ rb donnent le même reste.
    assert all((2 * x * x - x * y + 4) % m == res for x in range(ra - 5 * m, ra + 5 * m, m) for y in range(rb - 5 * m, rb + 5 * m, m))
    assert res not in leurres and len(set(leurres)) == 3

# ---------------------------------------------------------------------------
# arith-t-divisibilite-congruences-ex2 (et variantes) : grande puissance modulo un entier composé
# ---------------------------------------------------------------------------
# Exercice : 2^2026 mod 9 = 7 ; 2^3 ≡ -1, 2^6 ≡ 1, 2026 = 6 × 337 + 4, 16 ≡ 7, 2^8 ≡ 4.
assert congru(8, -1, 9) and pow(2, 6, 9) == 1 and 2026 == 6 * 337 + 4 and 16 % 9 == 7
assert pow(2, 2026, 9) == 7 and pow(2, 8, 9) == 4
assert all(x != 7 for x in (2, 4, 1))
# Variante 1 : 3^2027 mod 14 = 5 ; 27 ≡ -1, 2027 = 6 × 337 + 5, 243 = 14 × 17 + 5.
assert congru(27, -1, 14) and pow(3, 6, 14) == 1 and 2027 == 6 * 337 + 5 and 243 == 14 * 17 + 5
assert pow(3, 2027, 14) == 5 and all(x != 5 for x in (9, 3, 1))
# Variante 2 : 4^2026 mod 21 = 4 ; 64 = 3 × 21 + 1, 2026 = 3 × 675 + 1.
assert 64 == 3 * 21 + 1 and 2026 == 3 * 675 + 1 and pow(4, 2026, 21) == 4
assert all(x != 4 for x in (16, 1, 13))
# 9, 14, 21 ne sont pas premiers (Fermat ne s'applique pas).
assert 9 == 3 * 3 and 14 == 2 * 7 and 21 == 3 * 7

# ---------------------------------------------------------------------------
# arith-t-divisibilite-congruences-ex3 (et variantes) : divisibilité pour tout n
# ---------------------------------------------------------------------------
assert 64 == 5 * 11 + 9 and 9 % 5 == 4 and 25 % 17 == 8
for k in range(0, 60):
    assert (3 ** (2 * k + 2) + 2 ** (6 * k + 1)) % 11 == 0   # exercice
    assert (2 ** (2 * k + 1) + 3 ** (2 * k + 1)) % 5 == 0    # variante 1
    assert (3 * 5 ** (2 * k + 1) + 2 ** (3 * k + 1)) % 17 == 0  # variante 2
assert 9 + 2 == 11 and 2 + 3 == 5 and 15 + 2 == 17  # cas n = 0

# ---------------------------------------------------------------------------
# arith-t-divisibilite-congruences-ex4 (et variantes) : ax ≡ b (mod n) non simplifiable
# ---------------------------------------------------------------------------


def solutions(a, b, m):
    return {x for x in range(m) if (a * x - b) % m == 0}


# Exercice : 6x ≡ 4 (mod 10) <=> 3x ≡ 2 (mod 5) <=> x ≡ 4 (mod 5) : {4, 9} ; élève : {4} ; 6x ≡ 3 : aucune.
assert solutions(6, 4, 10) == {4, 9} == {x for x in range(10) if x % 5 == 4}
assert all(((6 * x - 4) % 10 == 0) == ((3 * x - 2) % 5 == 0) for x in range(-100, 100))
assert (3 * 2) % 5 == 1 and (2 * 2) % 5 == 4 and (3 * 4) % 5 == 2
assert solutions(3, 2, 10) == {4} and (3 * 7) % 10 == 1 and 14 % 10 == 4
assert solutions(6, 3, 10) == set()
assert 24 % 10 == 4 and 54 % 10 == 4
# Leurres : « x ≡ 4 (mod 10) » ({4}), « x ≡ 9 (mod 10) » ({9}), « x ≡ 2 (mod 5) » ({2, 7}) sont faux.
assert {4} != solutions(6, 4, 10) and {9} != solutions(6, 4, 10) and {2, 7} != solutions(6, 4, 10)
# Variante 1 : 9x ≡ 6 (mod 15) : {4, 9, 14} ; élève 3x ≡ 2 (mod 15) : aucune ; 9x ≡ 5 : aucune.
assert solutions(9, 6, 15) == {4, 9, 14} == {x for x in range(15) if x % 5 == 4}
assert all(((9 * x - 6) % 15 == 0) == ((3 * x - 2) % 5 == 0) for x in range(-100, 100))
assert solutions(3, 2, 15) == set() and solutions(9, 5, 15) == set()
assert 36 % 15 == 6 and 81 % 15 == 6 and 126 % 15 == 6
assert {4} != solutions(9, 6, 15) and set() != solutions(9, 6, 15) and {2, 7, 12} != solutions(9, 6, 15)
# Variante 2 : 8x ≡ 12 (mod 28) : {5, 12, 19, 26} ; élève 2x ≡ 3 (mod 28) : aucune ; 8x ≡ 6 : aucune.
assert solutions(8, 12, 28) == {5, 12, 19, 26} == {x for x in range(28) if x % 7 == 5}
assert all(((8 * x - 12) % 28 == 0) == ((2 * x - 3) % 7 == 0) for x in range(-100, 100))
assert (2 * 4) % 7 == 1 and (3 * 4) % 7 == 5 and (2 * 5) % 7 == 3
assert solutions(2, 3, 28) == set() and solutions(8, 6, 28) == set()
assert [8 * x % 28 for x in (5, 12, 19, 26)] == [12, 12, 12, 12]
assert {5} != solutions(8, 12, 28) and set() != solutions(8, 12, 28) and {3, 10, 17, 24} != solutions(8, 12, 28)

print("arith_t_divisibilite : toutes les vérifications passent.")
