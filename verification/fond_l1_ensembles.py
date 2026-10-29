"""Vérifications calculatoires du chapitre « fond-l1-ensembles » (thème fondements, niveau l1).

Cartes content/cartes/fondements/l1-fond-l1-ensembles.json et exercices
content/exercices/fondements/l1-fond-l1-ensembles.json. Chaque assert cite l'identifiant de la carte ou de
l'exercice vérifié (variante comprise : « v1 », « v2 »). Les leurres des questions à choix sont vérifiés FAUX.
Lancer : python3 verification/fond_l1_ensembles.py (doit finir sans erreur).
"""

import itertools
import json
import pathlib
from math import ceil, comb, factorial

import sympy
from sympy import Rational, floor, simplify, symbols

RACINE = pathlib.Path(__file__).resolve().parent.parent
x, y = symbols("x y")


def est_injective(f, dom):
    images = [f(n) for n in dom]
    return len(set(images)) == len(images)


def surjections(n, m):
    """Nombre de surjections d'un ensemble à n éléments sur un ensemble à m éléments (force brute)."""
    return sum(1 for t in itertools.product(range(m), repeat=n) if len(set(t)) == m)


def partitions(elements):
    """Toutes les partitions d'une liste (récursif)."""
    if not elements:
        yield []
        return
    premier, reste = elements[0], elements[1:]
    for p in partitions(reste):
        for i in range(len(p)):
            yield p[:i] + [[premier] + p[i]] + p[i + 1:]
        yield [[premier]] + p


def nb_relations_equivalence_force_brute(n):
    """Compte les relations sur {0..n-1} réflexives, symétriques et transitives (parcours de toutes les relations)."""
    E = range(n)
    couples = [(a, b) for a in E for b in E]
    total = 0
    for masque in range(1 << len(couples)):
        R = {couples[i] for i in range(len(couples)) if masque >> i & 1}
        if all((a, a) in R for a in E) and all((b, a) in R for (a, b) in R) and all(
            (a, d) in R for (a, b) in R for (c, d) in R if b == c
        ):
            total += 1
    return total


# ---------------------------------------------------------------------------
# fond-l1-ensembles-composition (carte)
# ---------------------------------------------------------------------------
N = range(0, 50)
f_succ = lambda n: n + 1
g_pred = lambda n: 0 if n == 0 else n - 1
# Exemple / affirmation 3 : g∘f = id sur N, f injective non surjective (0 sans antécédent), g non injective, f∘g ≠ id.
assert all(g_pred(f_succ(n)) == n for n in N)
assert est_injective(f_succ, N) and all(f_succ(n) != 0 for n in N)
assert g_pred(0) == g_pred(1) == 0
assert f_succ(g_pred(0)) == 1
# Contre-exemple / affirmation 1 : E = G = {0}, F = {0,1}, f(0) = 0, g ≡ 0 : g∘f = id_E, f non surjective, g non injective.
F_ = {0, 1}
f0 = {0: 0}
g0 = {0: 0, 1: 0}
assert g0[f0[0]] == 0 and set(f0.values()) != F_ and g0[0] == g0[1]
# Contre-exemple : (g∘f)^{-1} = f^{-1}∘g^{-1} avec f(x) = x+1, g(x) = 2x ; g^{-1}∘f^{-1} diffère (en y = 0 : -1 contre -1/2).
gf = 2 * (x + 1)
inv_gf = sympy.solve(sympy.Eq(gf, y), x)[0]
assert simplify(inv_gf - (y / 2 - 1)) == 0
assert inv_gf.subs(y, 0) == -1 and ((y - 1) / 2).subs(y, 0) == Rational(-1, 2)

# fond-l1-ensembles-composition-ex1 : n ↦ 2n+3 injective non surjective (valeurs impaires) ; v1 : 5-n bijective
# (involution) ; v2 : ⌊n/2⌋ surjective non injective.
Z = range(-60, 61)
assert est_injective(lambda n: 2 * n + 3, Z) and all((2 * n + 3) % 2 == 1 for n in Z)
assert all(5 - (5 - n) == n for n in Z)
assert all((2 * m) // 2 == m for m in Z) and 0 // 2 == 1 // 2 == 0

# fond-l1-ensembles-composition-ex2 (et v1, v2) : réciproques des homographies, leurres faux.
cas_homographies = [
    # (f, pôle, valeur exclue, réciproque attendue, leurres)
    ((2 * x + 1) / (x - 3), 3, 2, (3 * y + 1) / (y - 2), [(y - 3) / (2 * y + 1), (3 * y - 1) / (y - 2), (2 * y + 1) / (y - 3)]),
    ((x + 4) / (x - 1), 1, 1, (y + 4) / (y - 1), [(y - 1) / (y + 4), (y - 4) / (y - 1), (y + 4) / (y + 1)]),
    ((3 * x - 2) / (x + 2), -2, 3, (2 * y + 2) / (3 - y), [(y + 2) / (3 * y - 2), (2 * y - 2) / (3 - y), (2 * y + 2) / (y - 3)]),
]
for f, pole, exclue, inv, leurres in cas_homographies:
    assert sympy.solve(sympy.Eq(f, exclue), x) == []  # la valeur exclue n'est pas atteinte
    assert sympy.solve(sympy.Eq(inv, pole), y) == []  # la réciproque ne prend pas la valeur du pôle
    assert simplify(f.subs(x, inv) - y) == 0 and simplify(inv.subs(y, f) - x) == 0
    for l in leurres:
        assert simplify(f.subs(x, l) - y) != 0
# Formes réduites citées dans les corrections.
assert simplify((2 * x + 1) / (x - 3) - (2 + 7 / (x - 3))) == 0 and simplify((3 * y + 1) / (y - 2) - (3 + 7 / (y - 2))) == 0
assert simplify((x + 4) / (x - 1) - (1 + 5 / (x - 1))) == 0
assert simplify((3 * x - 2) / (x + 2) - (3 - 8 / (x + 2))) == 0 and simplify((2 * y + 2) / (3 - y) - (-2 + 8 / (3 - y))) == 0

# fond-l1-ensembles-composition-ex3 : surjections 4 → 3 : 36 ; v1 : 5 → 3 : 150 ; v2 : 5 → 2 : 30 ; leurres faux.
assert surjections(4, 3) == 36 == comb(4, 2) * factorial(3) == 3**4 - 3 * 2**4 + 3
assert all(v != 36 for v in [81, 72, 33]) and 3**4 == 81 and 3**4 - 3 * 2**4 == 33 and 4 * 3 * 2 * 3 == 72
assert surjections(5, 3) == 150 == 3 * comb(5, 3) * 2 + 3 * 5 * comb(4, 2) == 3**5 - 3 * 2**5 + 3
assert all(v != 150 for v in [243, 540, 147]) and 3**5 - 3 * 2**5 == 147 and 5 * 4 * 3 * 3**2 == 540
assert surjections(5, 2) == 30 == 2**5 - 2
assert all(v != 30 for v in [32, 31, 160]) and 5 * 4 * 2**3 == 160

# ---------------------------------------------------------------------------
# fond-l1-ensembles-equivalence (carte)
# ---------------------------------------------------------------------------
# Contre-exemple / affirmation 1 : {(1,1),(1,2),(2,1),(2,2)} sur {1,2,3} symétrique, transitive, non réflexive.
Rel = {(1, 1), (1, 2), (2, 1), (2, 2)}
assert all((b, a) in Rel for (a, b) in Rel)
assert all((a, d) in Rel for (a, b) in Rel for (c, d) in Rel if b == c)
assert (3, 3) not in Rel and not any(a == 3 for (a, b) in Rel)
# Contre-exemple : |x - y| ≤ 1 non transitive (0~1, 1~2, 0≁2) ; [-1,1] ∩ [0,2] = [0,1] ≠ ∅ et [-1,1] ≠ [0,2].
assert abs(0 - 1) <= 1 and abs(1 - 2) <= 1 and not abs(0 - 2) <= 1
assert sympy.Interval(-1, 1).intersect(sympy.Interval(0, 2)) == sympy.Interval(0, 1)
# Affirmation 3 : réunion des congruences modulo 2 et modulo 3 non transitive (0~2, 2~5, 0≁5).
rel23 = lambda a, b: (a - b) % 2 == 0 or (a - b) % 3 == 0
assert rel23(0, 2) and rel23(2, 5) and not rel23(0, 5)

# fond-l1-ensembles-equivalence-ex1 : classe de -3 pour x² : {-3, 3} ; leurres faux.
assert set(sympy.solve(x**2 - 9, x)) == {-3, 3}
assert {-3} != {-3, 3} and 9**2 != 9 and 0**2 != 9  # {−3} incomplet, {9} faux, [−3,3] contient 0
# v1 : x² − 2x, classe de 3 : {-1, 3} ; leurres {3}, {1,3} (f(1) = −1 ≠ 3), {−3,3} (f(−3) = 15 ≠ 3).
fv1 = lambda t: t**2 - 2 * t
assert fv1(3) == 3 and set(sympy.solve(x**2 - 2 * x - 3, x)) == {-1, 3}
assert fv1(1) != 3 and fv1(-3) == 15
# v2 : partie entière, classe de 5/2 : [2, 3[ ; leurres : ]2,3] (contient 3, exclut 2), {2}, [2,3] (contient 3).
assert floor(Rational(5, 2)) == 2 and floor(2) == 2 and floor(3) == 3 and floor(Rational(29, 10)) == 2
assert floor(Rational(5, 2)) == 2 and Rational(5, 2) != 2  # {2} n'est pas la classe (5/2 y manque)

# fond-l1-ensembles-equivalence-ex2 (et v1, v2) : nombre de carrés modulo n : 7 → 4, 8 → 3, 11 → 6 ; leurres faux.
for n, attendu, leurres in [(7, 4, [7, 3, 2]), (8, 3, [8, 4, 2]), (11, 6, [11, 5, 3])]:
    carres = sorted({(t * t) % n for t in range(n)})
    assert len(carres) == attendu and attendu not in leurres
assert [(t * t) % 7 for t in range(7)] == [0, 1, 4, 2, 2, 4, 1]
assert [(t * t) % 8 for t in range(8)] == [0, 1, 4, 1, 0, 1, 4, 1]
assert [(t * t) % 11 for t in range(11)] == [0, 1, 4, 9, 5, 3, 3, 5, 9, 4, 1]

# fond-l1-ensembles-equivalence-ex3 (et v1, v2) : nombres de Bell 15, 5, 52 ; décomptes par type ; leurres faux.
assert nb_relations_equivalence_force_brute(3) == 5 and nb_relations_equivalence_force_brute(4) == 15
assert sum(1 for _ in partitions(list(range(5)))) == 52 == sympy.bell(5)
assert 1 + 4 + 3 + 6 + 1 == 15 and 1 + 3 + 1 == 5 and 1 + 5 + 10 + 10 + 15 + 10 + 1 == 52
assert 1 + 4 + 6 + 6 + 1 == 18 and 1 + 5 + 10 + 10 + 30 + 10 + 1 == 67  # leurres « paires comptées deux fois »
assert all(v != 15 for v in [16, 24, 18]) and all(v != 5 for v in [6, 8, 4]) and all(v != 52 for v in [67, 120, 32])
# Type 2+2 sur 4 éléments : 3 partitions ; type 2+2+1 sur 5 éléments : 15.
assert sum(1 for p in partitions(list(range(4))) if sorted(map(len, p)) == [2, 2]) == 3
assert sum(1 for p in partitions(list(range(5))) if sorted(map(len, p)) == [1, 2, 2]) == 15

# fond-l1-ensembles-equivalence-ex4 : {(1,2)} sur {1,2} est circulaire (vacuement), ni réflexive ni symétrique.
Rc = {(1, 2)}
assert all((c, a) in Rc for (a, b) in Rc for (b2, c) in Rc if b == b2)
assert (1, 1) not in Rc and (2, 1) not in Rc
# Et toute relation réflexive circulaire sur un ensemble à 3 éléments est une relation d'équivalence (force brute).
E3 = range(3)
cpl = [(a, b) for a in E3 for b in E3]
for masque in range(1 << 9):
    R = {cpl[i] for i in range(9) if masque >> i & 1}
    refl = all((a, a) in R for a in E3)
    circ = all((c, a) in R for (a, b) in R for (b2, c) in R if b == b2)
    sym = all((b, a) in R for (a, b) in R)
    trans = all((a, d) in R for (a, b) in R for (c, d) in R if b == c)
    if refl and circ:
        assert sym and trans
    if refl and sym and trans:
        assert circ

# fond-l1-ensembles-equivalence-ex5 : sur Z×Z, (1,0)~(0,0)~(0,1) mais (1,0)≁(0,1).
sim = lambda a, b: a[0] * b[1] == b[0] * a[1]
assert sim((1, 0), (0, 0)) and sim((0, 0), (0, 1)) and not sim((1, 0), (0, 1))
# Sur Z×(Z∖{0}) : transitivité vérifiée sur une grille.
grille = [(p, q) for p in range(-3, 4) for q in range(-3, 4) if q != 0]
assert all(sim(a, c) for a in grille for b in grille for c in grille if sim(a, b) and sim(b, c))

# ---------------------------------------------------------------------------
# fond-l1-ensembles-recurrence-forte (carte)
# ---------------------------------------------------------------------------
# Exemple : tout entier n ≥ 2 est produit de nombres premiers (contrôle sur 2..500).
assert all(sympy.prod([p**k for p, k in sympy.factorint(n).items()]) == n for n in range(2, 501))
# Contre-exemple : {1/(n+1)} sans plus petit élément (1/(n+2) < 1/(n+1)).
assert all(Rational(1, n + 2) < Rational(1, n + 1) for n in range(100))
# Contre-exemple / affirmation 3 : « n > n + 1 » n'est jamais vraie.
assert not any(n > n + 1 for n in range(100))

# fond-l1-ensembles-recurrence-forte-ex1 (et v1, v2) : u_{n+1} = u_0 + ... + u_n.
def u_somme(u0, nmax):
    u = [u0]
    for _ in range(nmax):
        u.append(sum(u))
    return u
for u0, N_, attendu, leurres in [(1, 10, 512, [1024, 256, 1023]), (3, 8, 384, [768, 192, 765]), (5, 7, 320, [640, 160, 635])]:
    u = u_somme(u0, N_)
    assert u[N_] == attendu == u0 * 2 ** (N_ - 1) and attendu not in leurres
    assert all(u[n] == u0 * 2 ** (n - 1) for n in range(1, N_ + 1))
assert u_somme(3, 4)[1:] == [3, 6, 12, 24] and u_somme(5, 4)[1:] == [5, 10, 20, 40] and u_somme(1, 4)[1:] == [1, 2, 4, 8]

# fond-l1-ensembles-recurrence-forte-ex2 (et v1, v2) : plus grand entier non représentable a x + b y.
def representable(n, a, b):
    return any((n - b * t) % a == 0 for t in range(n // b + 1))
for a, b, seuil, attendu, leurres in [(3, 5, 8, 7, [8, 4, 15]), (3, 7, 12, 11, [12, 8, 21]), (5, 7, 24, 23, [24, 18, 35])]:
    assert all(representable(n, a, b) for n in range(seuil, 300))
    assert not representable(attendu, a, b) and attendu == a * b - a - b and attendu == seuil - 1
    for l in leurres:  # chaque leurre est représentable, ou bien non représentable mais plus petit que la réponse
        assert representable(l, a, b) or l < attendu
assert not representable(4, 3, 5) and not representable(8, 3, 7) and not representable(18, 5, 7)
# Écritures explicites des corrections.
assert 8 == 3 + 5 and 9 == 3 * 3 and 10 == 5 * 2
assert 12 == 3 * 4 and 13 == 3 * 2 + 7 and 14 == 7 * 2
assert 24 == 5 * 2 + 7 * 2 and 25 == 5 * 5 and 26 == 5 + 7 * 3 and 27 == 5 * 4 + 7 and 28 == 7 * 4
assert [23 - 7 * t for t in range(4)] == [23, 16, 9, 2]

# fond-l1-ensembles-recurrence-forte-ex3 : 1 < √2 < 2 et (q√2 − q)√2 = 2q − q√2.
q_ = symbols("q", positive=True)
assert sympy.expand((q_ * sympy.sqrt(2) - q_) * sympy.sqrt(2)) == sympy.expand(2 * q_ - q_ * sympy.sqrt(2))
assert 1 < sympy.sqrt(2) < 2

# fond-l1-ensembles-recurrence-forte-ex4 (et v1, v2) : Zeckendorf de 100, 50, 200 ; leurres de bonne somme mais
# avec deux termes de Fibonacci d'indices consécutifs.
F = [0, 1]
while len(F) < 30:
    F.append(F[-1] + F[-2])
indice = {F[k]: k for k in range(2, 30)}  # F_k, k ≥ 2, strictement croissante : indice bien défini
assert all(F[k + 1] > F[k] for k in range(2, 29)) and all(F[k] >= k - 1 for k in range(2, 30))
def glouton(n):
    termes = []
    while n:
        k = max(k for k in range(2, 30) if F[k] <= n)
        termes.append(F[k])
        n -= F[k]
    return termes
def valide(termes):
    ks = sorted(indice[t] for t in termes)
    return len(set(ks)) == len(ks) and all(ks[i + 1] - ks[i] >= 2 for i in range(len(ks) - 1))
for n, attendu, leurres in [
    (100, [89, 8, 3], [[55, 34, 8, 3], [89, 8, 2, 1], [89, 5, 3, 2, 1]]),
    (50, [34, 13, 3], [[34, 8, 5, 3], [21, 13, 8, 5, 3], [34, 13, 2, 1]]),
    (200, [144, 55, 1], [[144, 34, 21, 1], [89, 55, 34, 21, 1], [144, 34, 13, 8, 1]]),
]:
    assert glouton(n) == attendu and valide(attendu) and sum(attendu) == n
    for l in leurres:
        assert sum(l) == n and not valide(l)
assert (F[11], F[12], F[6], F[7], F[4]) == (89, 144, 8, 13, 3)
assert (F[9], F[10], F[7], F[8]) == (34, 55, 13, 21) and (F[12], F[13], F[10], F[11], F[2]) == (144, 233, 55, 89, 1)
# Existence de la décomposition pour tout n de 1 à 2000 (glouton valide).
assert all(valide(glouton(n)) and sum(glouton(n)) == n for n in range(1, 2001))

# ---------------------------------------------------------------------------
# fond-l1-ensembles-tiroirs (carte)
# ---------------------------------------------------------------------------
# Exemple : ⌈25/12⌉ = 3 ; contre-exemple : 7 objets, 3 tiroirs, répartition 3,2,2, ⌈7/3⌉ = 3.
assert ceil(25 / 12) == 3 and ceil(Rational(7, 3)) == 3 and 3 + 2 + 2 == 7 and max(3, 2, 2) < 4
# Affirmations : ⌈10/3⌉ = 4 (vrai) ; 12 = 4+4+4 sans tiroir à 5 (faux).
assert ceil(Rational(10, 3)) == 4 and 3 * 3 < 10 and 4 + 4 + 4 == 12
# Trou 1 : leurres faux. ⌊n/m⌋+1 : n=4, m=2, répartition 2,2 (pas de tiroir à 3) ; n−m : n=10, m=2, 5,5 ;
# ⌈m/n⌉ : n=1, m=2 (un seul objet).
assert 4 // 2 + 1 == 3 and max(2, 2) < 3
assert 10 - 2 == 8 and max(5, 5) < 8
assert ceil(Rational(2, 1)) == 2 and 1 < 2
# Énoncé : vérification exhaustive pour n ≤ 6 et m ≤ 4 (toutes les applications).
for n in range(0, 7):
    for m in range(1, 5):
        for t in itertools.product(range(m), repeat=n):
            assert max([t.count(j) for j in range(m)]) >= ceil(Rational(n, m))

# fond-l1-ensembles-tiroirs-ex1 (et v1, v2) : minimum (k−1)m + 1 ; leurres faux.
for k, m, attendu, leurres in [(3, 7, 15, [21, 14, 22]), (4, 12, 37, [48, 36, 49]), (5, 7, 29, [35, 28, 36])]:
    assert attendu == (k - 1) * m + 1 and ceil(Rational(attendu, m)) >= k and ceil(Rational(attendu - 1, m)) < k
    assert attendu not in leurres

# fond-l1-ensembles-tiroirs-ex2 (et v1, v2) : plus petit k = m + 1 (restes modulo m) ; force brute pour les deux
# premiers cas, argument des restes pour le troisième ; leurres faux.
def plus_petit_k_diff(Nmax, m):
    for k in range(1, Nmax + 1):
        if all(any((a - b) % m == 0 for a, b in itertools.combinations(S, 2)) for S in itertools.combinations(range(1, Nmax + 1), k)):
            return k
assert plus_petit_k_diff(20, 7) == 8 and plus_petit_k_diff(30, 5) == 6
assert all((a - b) % 9 != 0 for a, b in itertools.combinations(range(1, 10), 2)) and len({r % 9 for r in range(1, 51)}) == 9
assert all(v != 8 for v in [7, 14, 3]) and all(v != 6 for v in [5, 10, 26]) and all(v != 10 for v in [9, 18, 6])

# fond-l1-ensembles-tiroirs-ex3 (et v1, v2) : plus petit k tel que toute k-partie de {1..N} contienne deux
# consécutifs = (taille maximale d'une partie sans consécutifs) + 1, calculée par programmation dynamique.
def max_sans_consecutifs(Nmax):
    pris, libre = 0, 0  # meilleure taille sur {1..i} selon que i est pris ou non
    for _ in range(Nmax):
        pris, libre = libre + 1, max(pris, libre)
    return max(pris, libre)
def plus_petit_k_cons_force_brute(Nmax):
    for k in range(1, Nmax + 1):
        if all(any(b - a == 1 for a, b in zip(S, S[1:])) for S in itertools.combinations(range(1, Nmax + 1), k)):
            return k
assert all(max_sans_consecutifs(n) + 1 == plus_petit_k_cons_force_brute(n) for n in range(2, 15))
for Nmax, attendu, leurres in [(30, 16, [15, 17, 30]), (40, 21, [20, 22, 40]), (25, 14, [13, 15, 12])]:
    assert max_sans_consecutifs(Nmax) + 1 == attendu and attendu not in leurres
    impairs = list(range(1, Nmax + 1, 2))
    assert len(impairs) == attendu - 1 and all(b - a == 2 for a, b in zip(impairs, impairs[1:]))

# fond-l1-ensembles-tiroirs-ex4 : pour n ≤ 6, toute (n+1)-partie de {1..2n} contient a | b (a ≠ b), et
# {n+1..2n} n'en contient pas.
for n in range(1, 7):
    for S in itertools.combinations(range(1, 2 * n + 1), n + 1):
        assert any(b % a == 0 for a, b in itertools.combinations(S, 2))
    assert not any(b % a == 0 for a, b in itertools.combinations(range(n + 1, 2 * n + 1), 2))

# fond-l1-ensembles-tiroirs-ex5 : Erdős–Szekeres. 3,2,1,6,5,4,9,8,7 n'a pas de sous-suite monotone de longueur 4 ;
# toute permutation de 5 (n = 2) ou de 10 éléments (n = 3, échantillon) en a une de longueur n + 1.
def plus_longue(seq, croissante=True):
    L_ = [1] * len(seq)
    for i in range(len(seq)):
        for j in range(i):
            if (seq[j] < seq[i]) == croissante and seq[j] != seq[i]:
                L_[i] = max(L_[i], L_[j] + 1)
    return max(L_)
s9 = [3, 2, 1, 6, 5, 4, 9, 8, 7]
assert len(set(s9)) == 9 and plus_longue(s9, True) == 3 and plus_longue(s9, False) == 3
assert all(max(plus_longue(p, True), plus_longue(p, False)) >= 3 for p in itertools.permutations(range(5)))
assert all(max(plus_longue(p, True), plus_longue(p, False)) >= 4 for p in itertools.islice(itertools.permutations(range(10)), 0, 200000, 7))
assert all(ceil(Rational(n * n + 1, n)) == n + 1 for n in range(1, 50))

# ---------------------------------------------------------------------------
# fond-l1-ensembles-denombrables (carte)
# ---------------------------------------------------------------------------
# (p, q) ↦ 2^p (2q + 1) − 1 est une bijection de N² sur N : contrôle sur {0..2^12 − 2}.
images = {2**p * (2 * q + 1) - 1 for p in range(13) for q in range(2**12)}
assert set(range(2**12 - 1)) <= images
assert len({(p, q) for p in range(13) for q in range(2**12)}) == len(images)  # injectivité sur ce domaine
# Bijection N → Z (démonstration, ex1) : n/2 si n pair, −(n+1)/2 si n impair.
phi = lambda n: n // 2 if n % 2 == 0 else -(n + 1) // 2
psi = lambda m: 2 * m if m >= 0 else -2 * m - 1
assert all(psi(phi(n)) == n for n in range(500)) and all(phi(psi(m)) == m for m in range(-250, 250))
assert [phi(n) for n in range(6)] == [0, -1, 1, -2, 2, -3]
# Affirmation 4 : diagonale t(k) = 1 − s_k(k) diffère de chaque s_k (contrôle sur des suites aléatoires tronquées).
import random
random.seed(1)
s = [[random.randint(0, 1) for _ in range(20)] for _ in range(20)]
t = [1 - s[k][k] for k in range(20)]
assert all(t != s[k] for k in range(20))
# Contre-exemples / affirmation 2 : n ↦ 1_{n} (indicatrice tronquée) est injective de N dans {0,1}^N.
indic = lambda n: tuple(1 if k == n else 0 for k in range(40))
assert len({indic(n) for n in range(40)}) == 40

# fond-l1-ensembles-denombrables-ex1 (et v1, v2) : φ^{-1}(−7) = 13, φ^{-1}(5) = 10, φ^{-1}(−20) = 39 ; leurres faux.
for m, attendu, leurres in [(-7, 13, [14, 12, 15]), (5, 10, [9, 11, 5]), (-20, 39, [40, 41, 38])]:
    assert phi(attendu) == m and psi(m) == attendu and all(phi(l) != m for l in leurres)

# fond-l1-ensembles-denombrables-ex2 (et v1, v2) : σ(p, q) = (p+q)(p+q+1)/2 + q.
sigma = lambda p, q: (p + q) * (p + q + 1) // 2 + q
assert len({sigma(p, q) for p in range(60) for q in range(60)}) == 3600 and set(range(1830)) <= {sigma(p, q) for p in range(60) for q in range(60)}
assert (sigma(0, 0), sigma(1, 0), sigma(0, 1), sigma(2, 0)) == (0, 1, 2, 3)
assert [s_ * (s_ + 1) // 2 for s_ in range(10)] == [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]
for n, attendu, leurres in [(25, (2, 4), [(4, 2), (1, 4), (3, 3)]), (40, (4, 4), [(3, 5), (5, 4), (4, 3)]), (17, (3, 2), [(2, 3), (2, 2), (4, 1)])]:
    assert sigma(*attendu) == n and all(sigma(*l) != n for l in leurres)

# fond-l1-ensembles-denombrables-ex3 (et v1, v2) : β(A) = Σ 2^a, bijection P_f(N) → N ; β^{-1}(100) = {2,5,6},
# β^{-1}(77) = {0,2,3,6}, β^{-1}(255) = {0,...,7} ; leurres faux.
beta = lambda A: sum(2**a for a in A)
assert len({beta(A) for r in range(11) for A in itertools.combinations(range(10), r)}) == 2**10
assert {beta(A) for r in range(11) for A in itertools.combinations(range(10), r)} == set(range(2**10))
for n, attendu, leurres in [
    (100, {2, 5, 6}, [{3, 6, 7}, {4, 32, 64}, {1, 4, 5}]),
    (77, {0, 2, 3, 6}, [{1, 3, 4, 7}, {1, 4, 8, 64}, {2, 3, 6}]),
    (255, set(range(8)), [set(range(1, 9)), {8}, set(range(7))]),
]:
    assert beta(attendu) == n and all(beta(l) != n for l in leurres)
assert bin(100) == "0b1100100" and bin(77) == "0b1001101" and bin(255) == "0b11111111"

# ---------------------------------------------------------------------------
# fond-l1-ensembles-cantor (carte)
# ---------------------------------------------------------------------------
# Exemple : E = {1,2}, f(1) = {1}, f(2) = {1,2} : D = ∅, qui n'est pas une valeur de f.
fE = {1: {1}, 2: {1, 2}}
D = {e for e in fE if e not in fE[e]}
assert D == set() and D not in fE.values()
# Contre-exemple : 0,0999... = 0,1 (somme géométrique 9 Σ_{n≥2} 10^{-n} = 1/10).
n_ = symbols("n", integer=True)
assert sympy.summation(9 * Rational(1, 10) ** n_, (n_, 2, sympy.oo)) == Rational(1, 10)
# Démonstration : le nombre diagonal (chiffres 4 ou 5) est dans [0, 1[ : au plus 0,555... = 5/9 < 1.
assert sympy.summation(5 * Rational(1, 10) ** n_, (n_, 1, sympy.oo)) == Rational(5, 9)
# Trous : f(x) = E pour tout x atteint E et {x | x ∈ f(x)} = E ; f ≡ ∅ atteint ∅ ; x ↦ {x} injective.
E2 = {1, 2}
fE2 = {1: E2, 2: E2}
assert {e for e in E2 if e in fE2[e]} == E2 and E2 in fE2.values()
assert len({frozenset({e}) for e in E2}) == len(E2)
# Exhaustif : pour |E| ≤ 3, aucune application E → P(E) n'est surjective, et D n'est jamais atteinte.
for taille in range(0, 4):
    Ens = list(range(taille))
    parties = [frozenset(c) for r in range(taille + 1) for c in itertools.combinations(Ens, r)]
    for valeurs in itertools.product(parties, repeat=taille):
        fx = dict(zip(Ens, valeurs))
        Dx = frozenset(e for e in Ens if e not in fx[e])
        assert Dx not in valeurs and set(valeurs) != set(parties)

# fond-l1-ensembles-cantor-ex1 (et v1, v2) : chiffres diagonaux ; leurres différents de la réponse.
def diagonale(xs):
    return "".join("4" if xs[i][i] == "5" else "5" for i in range(len(xs)))
for xs, attendu, leurres in [
    (["5172", "3591", "2245", "8135"], "4454", ["5545", "4455", "5454"]),
    (["2718", "1415", "6457", "3335"], "5544", ["4455", "5545", "2455"]),
    (["9050", "1528", "4053", "7777"], "5445", ["4554", "9557", "5545"]),
]:
    assert diagonale(xs) == attendu and attendu not in leurres
    assert all(attendu[i] != xs[i][i] for i in range(4))

# fond-l1-ensembles-cantor-ex2 (et v1, v2) : parties diagonales {1,3}, {2}, {1,3,5} ; non atteintes ; leurres faux.
for fx, attendu, leurres in [
    ({1: {2, 3}, 2: {1, 2, 4}, 3: set(), 4: {4}}, {1, 3}, [{2, 4}, {1, 2, 3, 4}, {3}]),
    ({1: {1, 3}, 2: {1}, 3: {1, 2, 3}}, {2}, [{1, 3}, set(), {1, 2}]),
    ({1: set(), 2: {2, 5}, 3: {1, 2}, 4: {1, 3, 4}, 5: {1, 2}}, {1, 3, 5}, [{2, 4}, {1, 3}, {1, 2, 3, 4, 5}]),
]:
    Dx = {e for e in fx if e not in fx[e]}
    assert Dx == attendu and Dx not in fx.values() and all(l != Dx for l in leurres)

# fond-l1-ensembles-cantor-ex3 : t ↦ a + (b − a)t et x ↦ (x − a)/(b − a) sont réciproques.
a_, b_, t_ = symbols("a b t")
assert simplify(((a_ + (b_ - a_) * t_) - a_) / (b_ - a_) - t_) == 0

# ---------------------------------------------------------------------------
# Cohérence avec les fichiers JSON : la proposition « bonne » est celle attendue ci-dessus.
# ---------------------------------------------------------------------------
BONNES = {
    "fond-l1-ensembles-composition-ex1": ["injective et non surjective", "bijective", "surjective et non injective"],
    "fond-l1-ensembles-composition-ex2": [r"$\frac{3y + 1}{y - 2}$", r"$\frac{y + 4}{y - 1}$", r"$\frac{2y + 2}{3 - y}$"],
    "fond-l1-ensembles-composition-ex3": ["$36$", "$150$", "$30$"],
    "fond-l1-ensembles-equivalence-ex1": [r"$\{-3, 3\}$", r"$\{-1, 3\}$", r"$\left[2, 3\right[$"],
    "fond-l1-ensembles-equivalence-ex2": ["$4$", "$3$", "$6$"],
    "fond-l1-ensembles-equivalence-ex3": ["$15$", "$5$", "$52$"],
    "fond-l1-ensembles-recurrence-forte-ex1": ["$512$", "$384$", "$320$"],
    "fond-l1-ensembles-recurrence-forte-ex2": ["$7$", "$11$", "$23$"],
    "fond-l1-ensembles-recurrence-forte-ex4": [r"$89 + 8 + 3$", r"$34 + 13 + 3$", r"$144 + 55 + 1$"],
    "fond-l1-ensembles-tiroirs-ex1": ["$15$", "$37$", "$29$"],
    "fond-l1-ensembles-tiroirs-ex2": ["$8$", "$6$", "$10$"],
    "fond-l1-ensembles-tiroirs-ex3": ["$16$", "$21$", "$14$"],
    "fond-l1-ensembles-denombrables-ex1": ["$13$", "$10$", "$39$"],
    "fond-l1-ensembles-denombrables-ex2": [r"$(2, 4)$", r"$(4, 4)$", r"$(3, 2)$"],
    "fond-l1-ensembles-denombrables-ex3": [r"$\{2, 5, 6\}$", r"$\{0, 2, 3, 6\}$", r"$\{0, 1, 2, 3, 4, 5, 6, 7\}$"],
    "fond-l1-ensembles-cantor-ex1": [r"$0{,}4454$", r"$0{,}5544$", r"$0{,}5445$"],
    "fond-l1-ensembles-cantor-ex2": [r"$\{1, 3\}$", r"$\{2\}$", r"$\{1, 3, 5\}$"],
}
exos = json.loads((RACINE / "content/exercices/fondements/l1-fond-l1-ensembles.json").read_text(encoding="utf-8"))["exercices"]
# Exercices ajoutés le 2026-10-29 (enrichissement) : vérifiés dans enrichi_fond_l1_ensembles.py, écartés ici.
AJOUTS = {"fond-l1-ensembles-tiroirs-ex6", "fond-l1-ensembles-recurrence-forte-ex5", "fond-l1-ensembles-composition-ex6"}
assert AJOUTS <= {e["id"] for e in exos}
vus = set()
for e in exos:
    if e["id"] in AJOUTS:
        continue
    rfs = [e.get("reponse_finale")] + [v.get("reponse_finale") for v in e.get("variantes", [])]
    if rfs[0] is None:
        assert all(r is None for r in rfs) and e["id"] not in BONNES
        continue
    assert e["id"] in BONNES, e["id"]
    bonnes = [r["propositions"][r["bonne"]] for r in rfs]
    assert bonnes == BONNES[e["id"]], (e["id"], bonnes)
    vus.add(e["id"])
assert vus == set(BONNES)

print("fond_l1_ensembles : toutes les vérifications passent.")
