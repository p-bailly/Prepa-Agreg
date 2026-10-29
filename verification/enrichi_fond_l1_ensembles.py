"""Vérifications de l'enrichissement du chapitre « fond-l1-ensembles » (QCM et exercices ajoutés), par énumération.

Fichiers content/qcm/fond-l1-ensembles.json et fin de content/exercices/fondements/l1-fond-l1-ensembles.json
(fond-l1-ensembles-tiroirs-ex6, fond-l1-ensembles-recurrence-forte-ex5, fond-l1-ensembles-composition-ex6). Chaque
assert cite l'identifiant de la question ou de l'exercice (variantes : « v1 », « v2 ») ; les leurres sont vérifiés FAUX.
Lancer : python3 verification/enrichi_fond_l1_ensembles.py (doit finir sans erreur).
"""

import math
from fractions import Fraction as F
from itertools import combinations, product

# fond-l1-ensembles-q1 : (p, q) -> p + q est surjective (n = f(n, 0)) et pas injective (f(1, 0) = f(0, 1)).
f = lambda p, q: p + q
assert all(f(n, 0) == n for n in range(100))
assert f(1, 0) == f(0, 1)

# fond-l1-ensembles-q2 : f(n) = 2n, g(n) = n // 2 : g∘f = id, f∘g ≠ id (f(g(1)) = 0) ; ni f ni g bijective.
assert all((2 * n) // 2 == n for n in range(1000))
assert 2 * (1 // 2) == 0
assert 1 not in {2 * n for n in range(1000)} and 0 // 2 == 1 // 2

# fond-l1-ensembles-q3 : x R y ssi x - y ∈ Z ; [0, 1[ contient un représentant par classe (x - ⌊x⌋), deux
# éléments distincts de [0, 1[ ne sont pas équivalents ; 0 et 1 (dans [0, 1]) sont équivalents ; ]0, 1[ manque la
# classe de 0 ; Z est une seule classe.
for x in (F(-7, 3), F(5, 2), F(0), F(13, 7)):
    r = x - math.floor(x)
    assert 0 <= r < 1 and (x - r).denominator == 1
assert (F(1) - F(0)).denominator == 1
assert all(not (0 < k < 1) for k in range(-5, 6))

# fond-l1-ensembles-q4 (vrai) : u_0 = u_1 = 1, u_{n+2} = u_{n+1} + u_n : u_n ≤ 2^n.
u = [1, 1]
for n in range(300):
    u.append(u[-1] + u[-2])
assert all(u[n] <= 2**n for n in range(len(u)))

# fond-l1-ensembles-q5 (faux) : 1, 6, 11, 16, 21 : aucune somme de deux d'entre eux divisible par 5.
L = [1, 6, 11, 16, 21]
assert len(set(L)) == 5 and all((a + b) % 5 != 0 for a, b in combinations(L, 2))
# et parmi six entiers, deux ont une différence divisible par 5 (résidus) : vérifié sur toutes les classes.
assert all(len(set(t)) < 6 for t in product(range(5), repeat=6))

# fond-l1-ensembles-q6 : 3 couleurs (10 noires, 8 blanches, 6 grises) : 7 tirages assurent 3 de même couleur ;
# 6 ne suffisent pas (2 + 2 + 2). Leurres 6, 9, 4 : 6 et 4 insuffisants, 9 n'est pas le minimum.
stocks = (10, 8, 6)


def garanti(k):
    return all(max(t) >= 3 for t in product(*(range(s + 1) for s in stocks)) if sum(t) == k)


assert garanti(7) and not garanti(6) and not garanti(4) and garanti(9)
assert min(k for k in range(1, 25) if garanti(k)) == 7

# fond-l1-ensembles-q7 : rien de calculatoire (Cantor + réunion) ; contrôle fini : le nombre de parties finies de
# {0, ..., n-1} est 2^n (toutes), aucune contradiction à tester ici. On vérifie l'identification polynômes de degré
# ≤ k <-> Z^{k+1} sur un exemple : (a_0, ..., a_k) détermine le polynôme et réciproquement.
coeffs = (3, -1, 0, 2)
P = lambda x: sum(c * x**i for i, c in enumerate(coeffs))
assert [P(x) for x in range(4)] == [3, 4, 17, 54]

# fond-l1-ensembles-q8 (faux) : application constante de N dans N : image {0}, finie.
assert {0 * n for n in range(100)} == {0}


# fond-l1-ensembles-tiroirs-ex6 : approximation de Dirichlet. Questions 1-2 testées pour plusieurs x et N :
# il existe q ∈ {1..N} et p avec |qx - p| < 1/N.
def frac(t):
    return t - math.floor(t)


for x in (math.sqrt(2), math.pi, math.e, (1 + math.sqrt(5)) / 2, 0.123456789):
    for N in range(1, 60):
        rs = [frac(k * x) for k in range(N + 1)]
        cases = [min(int(r * N), N - 1) for r in rs]
        assert len(set(cases)) < N + 1  # principe des tiroirs
        ok = False
        for k in range(N + 1):
            for l in range(k):
                if cases[k] == cases[l]:
                    q, p = k - l, math.floor(k * x) - math.floor(l * x)
                    assert abs(q * x - p) < 1 / N and abs(x - p / q) < 1 / (q * N) <= 1 / q**2
                    ok = True
        assert ok


# Question 4 : ensemble A pour x rationnel (énumération exacte ; un rationnel est dans A ssi son écriture
# irréductible convient, car la condition est plus facile pour un dénominateur plus petit).
def A(x, qmax=100):
    s = set()
    for q in range(1, qmax):
        for p in range(-5 * q, 5 * q + 1):
            if abs(x - F(p, q)) < F(1, q * q):
                s.add(F(p, q))
    return s


assert A(F(1, 2)) == {F(0), F(1, 2), F(1)}  # exercice
assert A(F(1, 3)) == {F(0), F(1, 3), F(1, 2), F(1)}  # v1
assert A(F(2, 5)) == {F(0), F(1, 3), F(2, 5), F(1, 2), F(1)}  # v2
# Minoration utilisée : si p/q ≠ a/b, |a/b - p/q| ≥ 1/(bq).
for (a, b) in ((1, 2), (1, 3), (2, 5)):
    for q in range(1, 30):
        for p in range(-30, 31):
            if F(p, q) != F(a, b):
                assert abs(F(a, b) - F(p, q)) >= F(1, b * q)
# Pour x irrationnel (√2), A contient beaucoup de fractions (réduites de √2) : contrôle numérique.
x = math.sqrt(2)
fr = {F(p, q) for q in range(1, 500) for p in (math.floor(q * x), math.ceil(q * x)) if abs(x - p / q) < 1 / q**2}
assert len(fr) >= 5

# fond-l1-ensembles-recurrence-forte-ex5 : formules et u_10.
def suite(u0, u1, a, b, n):
    u = [u0, u1]
    while len(u) <= n:
        u.append(a * u[-1] + b * u[-2])
    return u


u = suite(2, 3, 3, -2, 30)
assert all(u[n] == 2**n + 1 for n in range(31)) and u[10] == 1025
assert all(v != 1025 for v in (1024, 2049, 1023))
u = suite(1, 4, 5, -6, 30)  # v1
assert all(u[n] == 2 * 3**n - 2**n for n in range(31)) and u[10] == 117074
assert all(v != 117074 for v in (118098, 58025, 352246))
u = suite(0, 1, 4, -4, 30)  # v2
assert all(u[n] == n * 2 ** (n - 1) if n >= 1 else u[n] == 0 for n in range(31)) and u[10] == 5120
assert all(v != 5120 for v in (10240, 1024, 4608))

# fond-l1-ensembles-composition-ex6 : φ(X) = (X ∩ A, X ∩ B), injective ssi A ∪ B = E, surjective ssi A ∩ B = ∅,
# énumération pour E de cardinal ≤ 4 et toutes les parties A, B.
for n in range(0, 5):
    E = frozenset(range(n))
    parties = [frozenset(c) for k in range(n + 1) for c in combinations(range(n), k)]
    for A_ in parties:
        for B_ in parties:
            images = [(X & A_, X & B_) for X in parties]
            injective = len(set(images)) == len(images)
            cible = {(Y, Z) for Y in parties if Y <= A_ for Z in parties if Z <= B_}
            surjective = set(images) == cible
            assert injective == (A_ | B_ == E)
            assert surjective == (not (A_ & B_))
            if injective and surjective:
                assert all((X & A_) | (X & B_) == X for X in parties)  # réciproque (Y, Z) -> Y ∪ Z
                assert all(((Y | Z) & A_, (Y | Z) & B_) == (Y, Z) for (Y, Z) in cible)
            # Leurre « A = B = E » : φ n'est pas bijective dès que E ≠ ∅.
            if n >= 1 and A_ == E and B_ == E:
                assert not surjective

print("enrichi_fond_l1_ensembles : tout est vérifié.")
