"""Vérifications calculatoires du chapitre « arith-l1-congruences » (thème arithmétique, niveau l1).

Cartes content/cartes/arithmetique/l1-arith-l1-congruences.json (inversibles, Euler, Wilson) et exercices
content/exercices/arithmetique/l1-arith-l1-congruences.json (la carte arith-l1-chinois, dans l1.json, est vérifiée
ailleurs). Chaque assert cite l'identifiant de la carte ou de l'exercice vérifié (variante comprise : « v1 », « v2 »).
Les leurres des questions à choix sont vérifiés FAUX ; la bonne réponse du JSON est relue.
Lancer : python3 verification/arith_l1_congruences.py (doit finir sans erreur).
"""

import json
import pathlib
from math import factorial, gcd, prod

from sympy import factorint, isprime, totient

RACINE = pathlib.Path(__file__).resolve().parent.parent
CARTES = {c["id"]: c for c in json.loads((RACINE / "content/cartes/arithmetique/l1-arith-l1-congruences.json").read_text())["cartes"]}
EXOS = {e["id"]: e for e in json.loads((RACINE / "content/exercices/arithmetique/l1-arith-l1-congruences.json").read_text())["exercices"]}
INV, EUL, WIL = "arith-l1-congruences-inversibles", "arith-l1-congruences-euler", "arith-l1-congruences-wilson"


def phi(n):
    return int(totient(n))


def rf(eid, v=0):
    e = EXOS[eid]
    return e["reponse_finale"] if v == 0 else e["variantes"][v - 1]["reponse_finale"]


def bonne(eid, v=0):
    r = rf(eid, v)
    return r["propositions"][r["bonne"]]


def leurres(eid, v=0):
    r = rf(eid, v)
    return [p for i, p in enumerate(r["propositions"]) if i != r["bonne"]]


def entier(t):
    """'$6\\,667$' -> 6667."""
    return int(t.strip("$").replace("\\,", "").replace(" ", ""))


def inversibles(n):
    return [k for k in range(n) if gcd(k, n) == 1]


# ============================================================== Cartes
# arith-l1-congruences-inversibles : exemple 17 mod 60
assert 60 == 3 * 17 + 9 and 17 == 9 + 8 and 9 == 8 + 1  # arith-l1-congruences-inversibles (exemple, Euclide)
assert 2 * 60 - 7 * 17 == 1 and 17 * 53 == 901 == 15 * 60 + 1  # arith-l1-congruences-inversibles (exemple, inverse 53)
assert 53 * 5 == 265 and 265 % 60 == 25 and 17 * 25 == 425 == 7 * 60 + 5  # arith-l1-congruences-inversibles (exemple, 17x = 5)
# cours : 8 non inversible modulo 12, 8*3 = 8*6 mod 12
assert gcd(8, 12) == 4 and (8 * 3 - 8 * 6) % 12 == 0 and (3 - 6) % 12 != 0  # arith-l1-congruences-inversibles (cours)
# contre-exemples : 6 mod 15
assert gcd(6, 15) == 3 and all((6 * x) % 15 != 1 for x in range(15))  # arith-l1-congruences-inversibles (contre-ex. 1)
assert 6 * 7 == 42 == 2 * 15 + 12 and (2 - 7) % 15 != 0 and (6 * 5) % 15 == 0  # arith-l1-congruences-inversibles (contre-ex. 2)
assert gcd(3, 6) == 3 and all((3 * x) % 6 != 1 for x in range(6))  # arith-l1-congruences-inversibles (contre-ex. 3)
# pièges : chacun rend l'énoncé faux
assert 2 % 4 != 0 and all((2 * x) % 4 != 1 for x in range(4))  # piège « n ne divise pas a » : a = 2, n = 4
assert isprime(3) and all((3 * x) % 6 != 1 for x in range(6))  # piège « a premier » : a = 3, n = 6
assert isprime(5) and all((0 * x) % 5 != 1 for x in range(5))  # piège « n premier, a quelconque » : a = 0
assert all((3 * x) % 9 != 1 for x in range(9))  # piège « a et n impairs » : a = 3, n = 9
# trous : leurres faux
assert all((2 * x) % 4 != 1 for x in range(4))  # trou 1, « n ne divise pas a » (a = 2, n = 4)
assert isprime(5) and gcd(0, 5) != 1 and 1 * 1 % 4 == 1  # trou 1, « n premier » (0 mod 5 non inversible ; 1 inversible mod 4)
assert (1 * 1 + 3 * 0 == 1) and (1 * 0) % 3 != 1  # trou 2, « classe de v » : a = 1, n = 3, u = 1, v = 0
assert (1 * (-1)) % 3 != 1  # trou 2, « classe de -u » : a = 1, n = 3, u = 1
assert pow(3, 8 - 2, 8) == 1 and (3 * 1) % 8 != 1  # trou 2, « a^(n-2) » : 3^6 = 1 mod 8 n'est pas l'inverse de 3
assert len(inversibles(4)) == 2 != 4 - 1 and gcd((1 + 1) % 4, 4) != 1  # trou 3, leurres « n - 1 éléments », « sous-anneau »
# affirmations
assert (1 + 1) % 4 == 2 and gcd(2, 4) == 2  # arith-l1-congruences-inversibles (affirmation 4)
assert (6 * 2 - 6 * 7) % 15 == 0  # arith-l1-congruences-inversibles (affirmation 2)
for n in range(1, 40):  # arith-l1-congruences-inversibles (affirmation 3 et énoncé, vérifiés exhaustivement pour n < 40)
    for a in range(n):
        inv = [b for b in range(n) if (a * b) % n == 1 % n]
        assert (len(inv) >= 1) == (gcd(a, n) == 1)
        if gcd(a, n) == 1:
            assert len(inv) == 1
            for b in range(n):
                assert len([x for x in range(n) if (a * x - b) % n == 0]) == 1

# arith-l1-congruences-euler
assert phi(360) == 96 == 4 * 6 * 4 and 360 * 1 * 2 * 4 == 96 * 2 * 3 * 5  # arith-l1-congruences-euler (exemple)
assert phi(36) == 12 and 100 == 8 * 12 + 4 and pow(5, 100, 36) == 13 and 625 == 17 * 36 + 13  # arith-l1-congruences-euler (exemple)
assert all(pow(2, k, 12) != 1 for k in range(1, 30)) and pow(2, phi(12), 12) == 4  # arith-l1-congruences-euler (cours, contre-ex.)
assert phi(4) == 2 and phi(2) ** 2 == 1  # arith-l1-congruences-euler (contre-ex. multiplicativité, affirmation 2)
assert 6 - 1 == 5 and phi(6) == 2  # arith-l1-congruences-euler (contre-ex. p = 6)
assert pow(3, 2, 8) == 1 and phi(8) == 4 and 2 % 4 != 0  # arith-l1-congruences-euler (contre-ex., affirmation 1)
assert pow(2, phi(4), 4) == 0  # piège « n ne divise pas a » : a = 2, n = 4 ; et trou 3 « premier » (a = 2)
assert pow(0, phi(5), 5) == 0  # trou 3 « quelconque » : a = 0, n = 5
assert phi(4) == 2 != 4 - 1 and phi(4) != (2 - 1) ** 2 and phi(4) != 4 * (2 - 1)  # trou 1 : leurres faux pour p = 2, alpha = 2
assert phi(8) != phi(2) * phi(4)  # trou 2 : « distincts », « m divise n », « quelconques » (m = 2, n = 4)
assert phi(4) != 4 - 1  # piège « p entier >= 2 quelconque » : p = 4, alpha = 1
for n in range(1, 300):  # arith-l1-congruences-euler : définition, formules, théorème, affirmations 3 et 4
    assert phi(n) == len([k for k in range(1, n + 1) if gcd(k, n) == 1])
    f = factorint(n)
    assert phi(n) == prod(p ** a - p ** (a - 1) for p, a in f.items())
    for a in range(n):
        if gcd(a, n) == 1:
            assert pow(a, phi(n), n) == 1 % n
    if n >= 3:
        assert phi(n) % 2 == 0
    if n >= 2:
        assert (phi(n) == n - 1) == isprime(n)
for m in range(1, 40):  # arith-l1-congruences-euler : multiplicativité
    for n in range(1, 40):
        if gcd(m, n) == 1:
            assert phi(m * n) == phi(m) * phi(n)

# arith-l1-congruences-wilson
for p in range(2, 400):  # arith-l1-congruences-wilson : théorème, réciproque, cas composé
    r = factorial(p - 1) % p
    if isprime(p):
        assert r == p - 1
        assert factorial(p - 2) % p == 1 % p  # affirmation 2
    else:
        assert r == (2 if p == 4 else 0)
assert factorial(10) == 3628800 == 11 * 329891 - 1  # arith-l1-congruences-wilson (exemple)
assert all(x * y % 11 == 1 for x, y in [(2, 6), (3, 4), (5, 9), (7, 8)])  # arith-l1-congruences-wilson (exemple, paires)
assert factorial(5) == 120 == 20 * 6  # arith-l1-congruences-wilson (exemple n = 6, affirmation 1)
assert factorial(6) == 720 and (2 * 4) % 7 == 1 and (3 * 5) % 7 == 1  # arith-l1-congruences-wilson (cours, p = 7)
assert 1 * 3 * 5 * 7 == 105 == 13 * 8 + 1  # arith-l1-congruences-wilson (cours, affirmation 3)
assert factorial(8) == 40320 == 9 * 4480  # arith-l1-congruences-wilson (contre-ex. n = 9 ; piège « impair »)
assert pow(2, 340, 341) == 1 and 341 == 11 * 31 and factorial(340) % 341 == 0  # (contre-ex. 341 ; piège « 2^(p-1) »)
assert factorial(0) % 1 == (-1) % 1  # (contre-ex. n = 1 ; piège « n >= 1 »)
assert factorial(3) % 4 == 2  # (contre-ex. n = 4)
assert factorial(6) + 1 == 721 == 7 * 103 and 721 % 49 != 0 and (factorial(4) + 1) % 25 == 0  # affirmation 4
assert factorial(2) % 3 != 1 and factorial(2) % 3 != 0 and (3 - 2) % 3 != factorial(2) % 3  # trou 2 : leurres faux (p = 3)
assert factorial(5) % 6 not in ((-1) % 6, 1, 2)  # trou 3 : leurres faux (n = 6)

# ============================================================== Exercices
# ---- inversibles-ex1
for v, (n, cl, inv) in enumerate([(90, [21, 35, 49, 77], [49, 77]), (84, [25, 35, 55, 63], [25, 55]),
                                    (120, [49, 77, 85, 91], [49, 77, 91])]):
    assert [c for c in cl if gcd(c, n) == 1] == inv  # arith-l1-congruences-inversibles-ex1 (v{v})
    assert all(c % n != 0 for c in cl)  # leurre « les quatre » faux ; « aucune » faux car inv non vide
    assert all(f"\\overline{{{c}}}" in bonne(f"{INV}-ex1", v) for c in inv)
    assert not any(f"\\overline{{{c}}}" in bonne(f"{INV}-ex1", v) for c in cl if c not in inv)
    seul = leurres(f"{INV}-ex1", v)[1]
    assert sum(f"\\overline{{{c}}}" in seul for c in inv) < len(inv)  # leurre « ... seulement » incomplet
assert gcd(21, 3) == 3 and gcd(85, 5) == 5 and gcd(35, 7) == 7 and gcd(63, 84) == 21  # ex1 (corrections)

# ---- inversibles-ex2
for v, (a, n, b, u, x, lr) in enumerate([(37, 100, 12, 73, 76, [73, 44, 24]), (53, 90, 11, 17, 7, [17, 43, 83]),
                                         (43, 80, 6, 67, 2, [67, 18, 78])]):
    assert pow(a, -1, n) == u and (u * b) % n == x  # arith-l1-congruences-inversibles-ex2 (v{v})
    assert [y for y in range(n) if (a * y - b) % n == 0] == [x]
    assert f"x \\equiv {x} " in bonne(f"{INV}-ex2", v)
    for l, t in zip(lr, leurres(f"{INV}-ex2", v)):
        assert f"x \\equiv {l} " in t and (a * l - b) % n != 0  # leurre faux
assert 10 * 100 - 27 * 37 == 1 and 37 * 73 == 2701 and 73 * 12 == 876 and 37 * 76 == 2812  # ex2
assert [100 - 2 * 37, 37 - 26, 26 - 2 * 11, 11 - 2 * 4, 4 - 3] == [26, 11, 4, 3, 1]  # ex2 (Euclide)
assert 17 * 53 - 10 * 90 == 1 and 53 * 7 == 371 == 4 * 90 + 11 and 17 * 11 == 187  # ex2 v1
assert [90 - 53, 53 - 37, 37 - 2 * 16, 16 - 3 * 5] == [37, 16, 5, 1]  # ex2 v1 (Euclide)
assert 7 * 80 - 13 * 43 == 1 and 43 * 67 == 2881 == 36 * 80 + 1 and 67 * 6 == 402  # ex2 v2
assert [80 - 43, 43 - 37, 37 - 6 * 6] == [37, 6, 1]  # ex2 v2 (Euclide)

# ---- inversibles-ex3 (systèmes)
for v, (a, b, c, d, e, f, n, x0, y0, sous) in enumerate([(3, 5, 7, 2, 1, 4, 26, 6, 7, 13), (5, 2, 2, 3, 5, 4, 21, 14, 20, 7),
                                                         (4, 7, 6, 4, 1, 3, 15, 7, 9, 5)]):
    sol = [(x, y) for x in range(n) for y in range(n) if (a * x + b * y - e) % n == 0 and (c * x + d * y - f) % n == 0]
    assert sol == [(x0, y0)]  # arith-l1-congruences-inversibles-ex3 (v{v}) : solution unique modulo n
    assert gcd(a * d - b * c, n) == 1
    assert f"x \\equiv {x0}$ et $y \\equiv {y0} \\pmod{{{n}}}" in bonne(f"{INV}-ex3", v)
    l1, l2, l3 = leurres(f"{INV}-ex3", v)
    assert f"x \\equiv {y0}$ et $y \\equiv {x0} \\pmod{{{n}}}" in l1 and (y0, x0) not in sol
    xo, yo = (-x0) % n, (-y0) % n
    assert f"x \\equiv {xo}$ et $y \\equiv {yo} \\pmod{{{n}}}" in l2 and (xo, yo) not in sol
    assert f"\\pmod{{{sous}}}" in l3 and ((x0 + sous) % n, y0) not in sol  # « modulo sous » contient des non-solutions
assert 3 * 9 == 27 and 23 * 17 == 391 == 15 * 26 + 1 and 9 * 21 == 189 and 17 * 8 == 136  # ex3 (inverses, calculs)

# ---- inversibles-ex4 (Newton)
for v, (a, b0, b1, b2) in enumerate([(3, 7, 67, 6667), (7, 3, 43, 7143), (9, 9, 89, 8889)]):
    assert (a * b0) % 10 == 1 and (b0 * (2 - a * b0)) % 100 == b1 and (a * b1) % 100 == 1  # ...-ex4 (v{v})
    assert (b1 * (2 - a * b1)) % 10000 == b2 and (a * b2) % 10000 == 1 and pow(a, -1, 10000) == b2
    assert entier(bonne(f"{INV}-ex4", v)) == b2
    assert all((a * entier(t)) % 10000 != 1 for t in leurres(f"{INV}-ex4", v))  # leurres faux
assert 7 * (2 - 21) == -133 and 67 * (2 - 201) == -13333 and 3 * (2 - 21) == -57 and 43 * (2 - 301) == -12857  # ex4
assert 9 * (2 - 81) == -711 and 89 * (2 - 801) == -71111 and 3 * 6667 == 20001 and 7 * 7143 == 50001 and 9 * 8889 == 80001
for m in range(2, 30):  # ex4 Q1 : identité 1 - ab(2 - ab) = (1 - ab)^2, et conséquence modulo m^2
    for a in range(m):
        for b in range(m):
            if (a * b) % m == 1:
                assert (a * b * (2 - a * b)) % (m * m) == 1 % (m * m)

# ---- inversibles-ex5 (diviseurs de zéro)
for v, (n, nb) in enumerate([(12, 7), (20, 11), (15, 6)]):
    zd = [a for a in range(1, n) if any((a * b) % n == 0 for b in range(1, n))]
    assert len(zd) == nb == n - 1 - phi(n)  # arith-l1-congruences-inversibles-ex5 (v{v})
    assert all(gcd(a, n) > 1 for a in zd) and set(zd).isdisjoint(inversibles(n))
    assert entier(bonne(f"{INV}-ex5", v)) == nb and all(entier(t) != nb for t in leurres(f"{INV}-ex5", v))
assert [a for a in range(1, 12) if gcd(a, 12) > 1] == [2, 3, 4, 6, 8, 9, 10]  # ex5 (liste)
assert [a for a in range(1, 20) if gcd(a, 20) > 1] == [2, 4, 5, 6, 8, 10, 12, 14, 15, 16, 18]  # ex5 v1 (liste)
assert [a for a in range(1, 15) if gcd(a, 15) > 1] == [3, 5, 6, 9, 10, 12]  # ex5 v2 (liste)
for n in range(2, 60):  # ex5 Q1-Q3 vérifiées exhaustivement
    for a in range(1, n):
        d = gcd(a, n)
        if d > 1:
            assert (a * (n // d)) % n == 0 and 1 <= n // d < n
        simplifiable = all((b - c) % n == 0 for b in range(n) for c in range(n) if (a * b - a * c) % n == 0)
        assert simplifiable == (d == 1)

# ---- euler-ex1
for v, (n, val, faux) in enumerate([(756, 216, [468, 12, 755]), (1350, 360, [624, 8, 1349]), (1200, 320, [720, 8, 1199])]):
    assert phi(n) == val == entier(bonne(f"{EUL}-ex1", v))  # arith-l1-congruences-euler-ex1 (v{v})
    f = factorint(n)
    assert prod(p ** a - 1 for p, a in f.items()) == faux[0] and prod(p - 1 for p in f) == faux[1] and faux[2] == n - 1
    assert [entier(t) for t in leurres(f"{EUL}-ex1", v)] == faux and val not in faux
assert factorint(756) == {2: 2, 3: 3, 7: 1} and factorint(1350) == {2: 1, 3: 3, 5: 2} and factorint(1200) == {2: 4, 3: 1, 5: 2}

# ---- euler-ex2
for v, (a, N, res) in enumerate([(13, 2026, 9), (17, 2027, 73), (23, 2026, 89)]):
    assert gcd(a, 100) == 1 and phi(100) == 40 and pow(a, N, 100) == res == pow(a, N % 40, 100)  # euler-ex2 (v{v})
    assert int(bonne(f"{EUL}-ex2", v).strip("$")) == res
    assert all(int(t.strip("$")) != res for t in leurres(f"{EUL}-ex2", v))  # leurres faux
assert [pow(13, 2 ** i, 100) for i in range(1, 5)] == [69, 61, 21, 41] and 41 * 21 == 861 and 61 * 69 == 4209  # ex2
assert [pow(17, 2 ** i, 100) for i in range(1, 5)] == [89, 21, 41, 81] and 81 * 41 == 3321 and 21 * 89 == 1869 and 69 * 17 == 1173
assert [pow(23, 2 ** i, 100) for i in range(1, 5)] == [29, 41, 81, 61] and 61 * 81 == 4941 and 41 * 29 == 1189  # ex2 v2

# ---- euler-ex3 (tours)
for v, (a, e, res) in enumerate([(13, 13 ** 13, 53), (9, 9 ** 9, 89), (3, 3 ** 2026, 83)]):
    assert pow(a, e, 100) == res and pow(a, e % 40, 100) == res  # arith-l1-congruences-euler-ex3 (v{v})
    assert int(bonne(f"{EUL}-ex3", v).strip("$")) == res
    assert all(int(t.strip("$")) != res for t in leurres(f"{EUL}-ex3", v))
assert pow(13, 2, 40) == 9 and pow(13, 4, 40) == 1 and pow(13, 13, 40) == 13 and phi(40) == 16  # ex3
assert pow(13, 13, 100) == 53 and 21 * 61 == 1281 and 81 * 13 == 1053  # ex3
assert pow(9, 2, 40) == 1 and pow(9, 9, 40) == 9 and pow(9, 4, 100) == 61 and pow(9, 8, 100) == 21 and 21 * 9 == 189  # ex3 v1
assert pow(3, 4, 40) == 1 and 2026 == 4 * 506 + 2 and pow(3, 2026, 40) == 9 and 3 ** 9 == 19683  # ex3 v2

# ---- euler-ex4 (phi(n) = k)
for v, (k, ens) in enumerate([(12, {13, 21, 26, 28, 36, 42}), (8, {15, 16, 20, 24, 30}), (14, set())]):
    # borne : phi(n) >= sqrt(n/2), donc phi(n) = k entraîne n <= 2k^2 ; on cherche bien au-delà
    assert {n for n in range(1, 10 * k * k) if phi(n) == k} == ens  # arith-l1-congruences-euler-ex4 (v{v})
    pr = rf(f"{EUL}-ex4", v)["propositions"]
    for i, t in enumerate(pr):
        s = set() if "vide" in t else {int(x) for x in t.strip("$").strip("\\{").strip("\\}").replace("\\{", "").replace("\\}", "").split(",")}
        assert (s == ens) == (i == rf(f"{EUL}-ex4", v)["bonne"])  # seule la bonne proposition est l'ensemble exact
assert [p for p in range(2, 30) if isprime(p) and 12 % (p - 1) == 0] == [2, 3, 5, 7, 13]  # ex4 Q2
assert [p for p in range(2, 30) if isprime(p) and 8 % (p - 1) == 0] == [2, 3, 5]  # ex4 v1 Q2
assert [p for p in range(2, 30) if isprime(p) and 14 % (p - 1) == 0] == [2, 3]  # ex4 v2 Q2
assert phi(18) == 6 and phi(52) == 24 and phi(32) == 16  # ex4 (leurres et remarque)
for n in range(2, 2000):  # ex4 Q1 : borne utilisée ci-dessus, phi(n) >= sqrt(n/2)
    assert 2 * phi(n) ** 2 >= n


def R(k):
    return (10 ** k - 1) // 9


# ---- euler-ex5 (repunits)
for v, (n0, kmin, faux) in enumerate([(7, 6, [36, 7, 3]), (41, 5, [240, 40, 41]), (37, 3, [216, 36, 37])]):
    ks = [k for k in range(1, 300) if R(k) % n0 == 0]
    assert ks[0] == kmin == entier(bonne(f"{EUL}-ex5", v))  # arith-l1-congruences-euler-ex5 (v{v})
    assert [entier(t) for t in leurres(f"{EUL}-ex5", v)] == faux and all(f != kmin for f in faux)
    assert R(phi(9 * n0)) % n0 == 0 and phi(9 * n0) == faux[0]  # ex5 Q3 : k = phi(9 n0)
    if n0 % 3:
        assert R(phi(n0)) % n0 == 0  # ex5 Q2
for n in range(1, 400):  # ex5 Q1-Q3 vérifiées pour n < 400 premier avec 10
    if gcd(n, 10) == 1:
        assert (pow(10, phi(n), n) - 1) % n == 0 and R(phi(9 * n)) % n == 0
        if n % 3:
            assert R(phi(n)) % n == 0
assert [R(k) % 7 for k in range(1, 7)] == [1, 4, 6, 5, 2, 0] and R(6) == 111111 == 7 * 15873  # ex5
assert [R(k) % 41 for k in range(1, 6)] == [1, 11, 29, 4, 0] and R(5) == 11111 == 41 * 271 and 291 % 41 == 4  # ex5 v1
assert [R(k) % 37 for k in range(1, 4)] == [1, 11, 0] and R(3) == 111 == 3 * 37  # ex5 v2
assert phi(63) == 36 and phi(369) == 240 and phi(333) == 216  # ex5

# ---- wilson-ex1
for v, (p, inv2, res) in enumerate([(23, 12, 11), (19, 10, 9), (29, 15, 14)]):
    assert isprime(p) and factorial(p - 3) % p == res == (p - 1) // 2  # arith-l1-congruences-wilson-ex1 (v{v})
    assert (2 * inv2) % p == 1 and (inv2 * (p - 1)) % p == res and (2 * res) % p == p - 1
    assert entier(bonne(f"{WIL}-ex1", v)) == res and all(entier(t) != res for t in leurres(f"{WIL}-ex1", v))
assert 12 * 22 == 264 == 11 * 23 + 11 and 10 * 18 == 180 == 9 * 19 + 9 and 15 * 28 == 420 == 14 * 29 + 14  # ex1

# ---- wilson-ex2
for p in [2, 3, 5, 7, 11, 13, 71]:  # ex2 Q1 : (k-1)! (p-k)! = (-1)^k mod p
    for k in range(1, p + 1):
        assert (factorial(k - 1) * factorial(p - k) - (-1) ** k) % p == 0
for v, (p, k, res) in enumerate([(71, 8, 70), (61, 9, 1), (103, 7, 1)]):
    assert isprime(p) and factorial(k - 1) % p == p - 1  # arith-l1-congruences-wilson-ex2 (v{v})
    assert factorial(p - k) % p == res
    assert entier(bonne(f"{WIL}-ex2", v)) == res and all(entier(t) != res for t in leurres(f"{WIL}-ex2", v))
assert factorial(7) == 5040 == 71 * 71 - 1 and factorial(8) == 40320 == 661 * 61 - 1 and factorial(6) == 720 == 7 * 103 - 1

# ---- wilson-ex3
for p in range(3, 200):  # ex3 Q1-Q3 : (m!)^2 = (-1)^(m+1) ; -1 carré mod p ssi p = 1 mod 4
    if isprime(p):
        m = (p - 1) // 2
        assert (factorial(m) ** 2 - (-1) ** (m + 1)) % p == 0
        assert any((x * x + 1) % p == 0 for x in range(p)) == (p % 4 == 1)
for v, (p, res, autre) in enumerate([(13, 5, 8), (17, 13, 4), (29, 12, 17)]):
    m = (p - 1) // 2
    assert factorial(m) % p == res and (res * res + 1) % p == 0 and (autre * autre + 1) % p == 0  # wilson-ex3 (v{v})
    assert sorted(x for x in range(p) if (x * x + 1) % p == 0) == sorted([res, autre])
    assert entier(bonne(f"{WIL}-ex3", v)) == res and all(entier(t) != res for t in leurres(f"{WIL}-ex3", v))
assert 720 == 55 * 13 + 5 and 169 == 10 * 17 - 1 and 144 == 5 * 29 - 1  # ex3
assert [factorial(k) % 17 for k in range(4, 9)] == [7, 1, 6, 8, 13]  # ex3 v1 (produits successifs)
assert [factorial(k) % 29 for k in range(4, 15)] == [24, 4, 24, 23, 10, 3, 1, 11, 16, 5, 12]  # ex3 v2

# ---- wilson-ex4
for v, (n, E, res) in enumerate([(16, [1, 7, 9, 15], 1), (20, [1, 9, 11, 19], 1), (9, [1, 8], 8)]):
    G = inversibles(n)
    assert [x for x in G if (x * x) % n == 1] == E  # arith-l1-congruences-wilson-ex4 (v{v})
    assert prod(G) % n == res == prod(E) % n
    assert entier(bonne(f"{WIL}-ex4", v)) == res and all(entier(t) != res for t in leurres(f"{WIL}-ex4", v))
assert 945 == 59 * 16 + 1 and 1881 == 94 * 20 + 1 and 2240 == 248 * 9 + 8 and prod([1, 2, 4, 5, 7, 8]) == 2240  # ex4
for n in range(3, 80):  # ex4 Q1-Q2 : P = produit de E ; E = {1, -1} entraîne P = -1
    G = inversibles(n)
    E = [x for x in G if (x * x) % n == 1]
    assert prod(G) % n == prod(E) % n and 1 in E and n - 1 in E
    if len(E) == 2:
        assert prod(G) % n == n - 1

print("arith_l1_congruences : toutes les vérifications passent.")
