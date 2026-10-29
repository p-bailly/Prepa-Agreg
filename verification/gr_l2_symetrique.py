"""Vérifications calculatoires du chapitre « gr-l2-symetrique » (groupes, L2) : cartes, exercices et QCM de
content/*/groupes/l2-gr-l2-symetrique.json et content/qcm/gr-l2-symetrique.json.

Convention du chapitre : sigma tau = sigma ∘ tau (on applique d'abord tau). Chaque assert cite l'identifiant de la carte,
de l'exercice (variante comprise : v1, v2) ou de la question. Calculs exacts par force brute sur des petits groupes
symétriques. Lancer : python3 verification/gr_l2_symetrique.py (doit finir sans erreur).
"""

from collections import deque
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial, gcd, lcm, prod

# ---------------------------------------------------------------------------
# Outils : une permutation de {1..n} est un tuple p avec p[i-1] = image de i.
# ---------------------------------------------------------------------------


def ident(n):
    return tuple(range(1, n + 1))


def cyc(n, *cs):
    """Produit de cycles à supports disjoints (ou un seul cycle), donné par listes."""
    p = list(range(1, n + 1))
    for c in cs:
        for a, b in zip(c, c[1:] + c[:1]):
            p[a - 1] = b
    return tuple(p)


def comp(*ps):
    """comp(s, t, ...) = s ∘ t ∘ ... (le dernier s'applique en premier)."""
    n = len(ps[0])
    r = ident(n)
    for p in reversed(ps):
        r = tuple(p[r[i] - 1] for i in range(n))
    return r


def prodc(n, *cs):
    """Produit (au sens de la composition) de cycles quelconques, chacun donné par une liste."""
    return comp(*[cyc(n, c) for c in cs])


def inverse(p):
    q = [0] * len(p)
    for i, x in enumerate(p, 1):
        q[x - 1] = i
    return tuple(q)


def puiss(p, k):
    r = ident(len(p))
    for _ in range(k % ordre(p) if k >= 0 else 0):
        r = comp(p, r)
    return r


def cycles(p):
    vu, res = set(), []
    for i in range(1, len(p) + 1):
        if i in vu:
            continue
        c, j = [i], p[i - 1]
        vu.add(i)
        while j != i:
            c.append(j)
            vu.add(j)
            j = p[j - 1]
        if len(c) > 1:
            res.append(tuple(c))
    return res


def norm(cs):
    """Forme normale d'une liste de cycles disjoints (chaque cycle commence par son plus petit élément)."""
    out = []
    for c in cs:
        k = c.index(min(c))
        out.append(tuple(c[k:] + c[:k]))
    return sorted(out)


def orbites(p):
    return len(cycles(p)) + sum(1 for i, x in enumerate(p, 1) if x == i)


def inv(p):
    n = len(p)
    return sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])


def eps(p):
    return (-1) ** inv(p)


def ordre(p):
    cs = cycles(p)
    return lcm(*[len(c) for c in cs]) if cs else 1


def type_(p):
    return tuple(sorted((len(c) for c in cycles(p)), reverse=True))


def engendre(gens):
    """Sous-groupe engendré (groupe fini : produits des générateurs)."""
    n = len(gens[0])
    vus = {ident(n)}
    file = deque([ident(n)])
    while file:
        x = file.popleft()
        for g in gens:
            y = comp(g, x)
            if y not in vus:
                vus.add(y)
                file.append(y)
    return vus


def distance(gens, n):
    """Longueur minimale d'écriture de chaque permutation comme produit de générateurs (involutions)."""
    d = {ident(n): 0}
    file = deque([ident(n)])
    while file:
        x = file.popleft()
        for g in gens:
            y = comp(x, g)
            if y not in d:
                d[y] = d[x] + 1
                file.append(y)
    return d


def S(n):
    return [tuple(p) for p in permutations(range(1, n + 1))]


def transp(n, a, b):
    return cyc(n, (a, b))


# ---------------------------------------------------------------------------
# gr-l2-symetrique-cycles (carte)
# ---------------------------------------------------------------------------
# Théorème (force brute, n <= 6) : produit des cycles des orbites = sigma ; supports disjoints => commutation.
for n in range(1, 7):
    for p in S(n):
        cs = cycles(p)
        assert comp(*([cyc(n, c) for c in cs] or [ident(n)])) == p
        assert sum(len(c) for c in cs) + (orbites(p) - len(cs)) == n
# Exemple : matrice de S_9 = (1 4 9 8)(2 7 5), 3 et 6 fixes.
sig_ex = (4, 7, 3, 9, 2, 6, 5, 1, 8)
assert sig_ex == cyc(9, (1, 4, 9, 8), (2, 7, 5)) and norm(cycles(sig_ex)) == [(1, 4, 9, 8), (2, 7, 5)]
assert comp(cyc(9, (1, 4, 9, 8)), cyc(9, (2, 7, 5))) == comp(cyc(9, (2, 7, 5)), cyc(9, (1, 4, 9, 8)))
# Contre-exemples : (1 2)(2 3) = (1 2 3), (2 3)(1 2) = (1 3 2), (1 2 3) = (1 3 2)(1 3 2).
assert prodc(3, (1, 2), (2, 3)) == cyc(3, (1, 2, 3)) and prodc(3, (2, 3), (1, 2)) == cyc(3, (1, 3, 2))
assert prodc(3, (1, 3, 2), (1, 3, 2)) == cyc(3, (1, 2, 3))
# Pièges des hypothèses : (1 3 2) et (1 3 2) commutent mais l'écriture (1 3 2)(1 3 2) n'est pas celle du théorème.
# Trou 1 : leurres « de longueur 2 », « de même longueur », « qui commutent » : (1 2 3) a deux écritures.
assert prodc(3, (1, 2), (2, 3)) == prodc(3, (1, 3), (1, 2)) == cyc(3, (1, 2, 3))  # deux écritures en transpositions
# Trou 2 : leurres faux pour sigma = (1 2 3) dans S_4 (orbite {4} sans cycle, réunion différente de {1..4},
# paire {1, 2} qui n'est pas un support).
p = cyc(4, (1, 2, 3))
assert orbites(p) == 2 and len(cycles(p)) == 1 and set(cycles(p)[0]) == {1, 2, 3} != {1, 2}
# Trou 3 : (1 2)(3 4 5) : longueurs différentes ; (1 2 3 4) d'ordre 4 non premier ; cycles d'une décomposition
# n'engendrent pas S_n (un seul cycle (1 2 3) engendre un groupe à 3 éléments).
assert ordre(cyc(4, (1, 2, 3, 4))) == 4 and len(engendre([cyc(3, (1, 2, 3))])) == 3
# Affirmations : (1 2 3) et son carré commutent ; (1 2)(3 4 5) sans point fixe n'est pas un cycle ;
# nombre de 3-cycles 2 C(n,3) = n(n-1)(n-2)/3 ; support d'un produit à supports disjoints.
c = cyc(3, (1, 2, 3))
assert comp(c, comp(c, c)) == comp(comp(c, c), c)
p = cyc(5, (1, 2), (3, 4, 5))
assert orbites(p) == 2 and all(p[i - 1] != i for i in range(1, 6)) and len(cycles(p)) == 2
for n in range(3, 8):
    assert sum(1 for q in S(n) if type_(q) == (3,)) == n * (n - 1) * (n - 2) // 3 == 2 * comb(n, 3)
for a, b in [((1, 2), (3, 4, 5)), ((1, 5, 2), (3, 4))]:
    q = comp(cyc(5, a), cyc(5, b))
    assert {i for i in range(1, 6) if q[i - 1] != i} == set(a) | set(b)

# gr-l2-symetrique-cycles-ex1 (+ v1, v2) : décompositions et leurres (inverses, cycles fusionnés, faux).
for images, bonne, leurres in [
    ((5, 8, 1, 7, 3, 6, 2, 4), [(1, 5, 3), (2, 8, 4, 7)],
     [[(1, 3, 5), (2, 7, 4, 8)], [(1, 5, 3), (2, 8, 4, 7, 6)], [(1, 5), (2, 8), (4, 7)]]),
    ((7, 9, 6, 1, 3, 8, 4, 5, 2), [(1, 7, 4), (2, 9), (3, 6, 8, 5)],
     [[(1, 4, 7), (2, 9), (3, 5, 8, 6)], [(1, 7, 4, 2, 9), (3, 6, 8, 5)], [(1, 7), (2, 9), (3, 6), (5, 8)]]),
    ((4, 6, 8, 1, 3, 7, 2, 5), [(1, 4), (2, 6, 7), (3, 8, 5)],
     [[(1, 4), (2, 7, 6), (3, 5, 8)], [(1, 4), (2, 6, 7, 3, 8, 5)], [(1, 4, 2, 6, 7), (3, 8, 5)]]),
]:
    n = len(images)
    assert norm(cycles(images)) == norm(bonne) and cyc(n, *bonne) == images
    for l in leurres:
        assert cyc(n, *l) != images
# Remarque de la correction d'ex1 : lire à l'envers donne l'inverse.
assert norm(cycles(inverse((5, 8, 1, 7, 3, 6, 2, 4)))) == [(1, 3, 5), (2, 7, 4, 8)]

# gr-l2-symetrique-cycles-ex2 (+ v1, v2) : produits de cycles non disjoints, et l'autre ordre.
for n, a, b, bonne, autre, leurres in [
    (5, (1, 2, 3, 4), (3, 4, 5), [(1, 2, 3), (4, 5)], [(1, 2, 4), (3, 5)], [[(1, 2, 3, 4, 5)], [(1, 3, 2), (4, 5)]]),
    (6, (1, 3, 5), (2, 3, 4, 6), [(1, 3, 4, 6, 2, 5)], [(1, 4, 6, 2, 3, 5)], [[(1, 5, 2, 6, 4, 3)], [(1, 3, 5), (2, 4, 6)]]),
    (5, (1, 4, 2), (2, 5, 3, 4), [(1, 4), (2, 5, 3)], [(1, 2), (3, 4, 5)], [[(1, 4), (2, 3, 5)], [(1, 4, 2, 5, 3)]]),
]:
    pi = prodc(n, a, b)
    assert pi == cyc(n, *bonne) and prodc(n, b, a) == cyc(n, *autre) != pi
    for l in leurres:
        assert cyc(n, *l) != pi


# gr-l2-symetrique-cycles-ex3 (+ v1, v2) : nombre de permutations d'un type donné, formule n!/prod k^m m!.
def nb_type(n, t):
    return sum(1 for q in S(n) if type_(q) == t)


def formule(n, t):
    t1 = list(t) + [1] * (n - sum(t))
    return factorial(n) // prod(k ** t1.count(k) * factorial(t1.count(k)) for k in set(t1))


for n, t, bonne, leurres in [(6, (3, 2), 120, [60, 240, 20]), (7, (3, 3), 280, [560, 140, 70]), (6, (4, 2), 90, [180, 45, 360])]:
    assert nb_type(n, t) == formule(n, t) == bonne and bonne not in leurres
assert comb(6, 3) * 2 * comb(3, 2) == 120 and comb(7, 3) * comb(4, 3) // 2 * 4 == 280 and comb(6, 4) * 6 == 90

# gr-l2-symetrique-cycles-ex4 : conjugaison d'un cycle ; commutant de (1 2 ... n) = puissances, n éléments.
for n in range(3, 7):
    c = cyc(n, tuple(range(1, n + 1)))
    comm = [t for t in S(n) if comp(t, c) == comp(c, t)]
    assert len(comm) == n and set(comm) == {puiss(c, k) for k in range(n)}
    assert n not in (factorial(n), factorial(n - 1), 2 * n)
for t in S(5):
    a = (1, 3, 4)
    assert comp(t, cyc(5, a), inverse(t)) == cyc(5, tuple(t[x - 1] for x in a))

# gr-l2-symetrique-cycles-ex5 : nombre de sigma avec l(1) = k égal à (n-1)! ; moyenne des orbites H_n.
import math

for n in range(2, 8):
    perms = S(n)
    for k in range(1, n + 1):
        assert sum(1 for q in perms if next(len(c) for c in cycles(q) + [(1,)] if 1 in c) == k) == factorial(n - 1)
    moy = F(sum(orbites(q) for q in perms), factorial(n))
    H = sum(F(1, k) for k in range(1, n + 1))
    assert moy == H and H != F(n, 2) and abs(float(H) - math.log(n)) > 0.1 and H * H != n  # H_n rationnel ; H_n^2 = n exclu exactement
assert sorted(orbites(q) for q in S(3)) == [1, 1, 2, 2, 2, 3]

# ---------------------------------------------------------------------------
# gr-l2-symetrique-signature (carte)
# ---------------------------------------------------------------------------
for n in range(2, 7):
    perms = S(n)
    for p in perms:
        e = eps(p)
        assert F(prod(p[j] - p[i] for i in range(n) for j in range(i + 1, n)),
                 prod(j - i for i in range(n) for j in range(i + 1, n))) == e  # formule produit
        assert e == (-1) ** (n - orbites(p)) == prod((-1) ** (len(c) - 1) for c in cycles(p))
    if n <= 5:
        for p in perms:
            for q in perms:
                assert eps(comp(p, q)) == eps(p) * eps(q)  # morphisme
    for a, b in combinations(range(1, n + 1), 2):
        assert inv(transp(n, a, b)) == 2 * (b - a - 1) + 1
# Décomposition d'un cycle : (a1 a2)(a2 a3)...(a_{k-1} a_k) = (a1 ... ak).
for k in range(2, 7):
    a = (2, 5, 1, 6, 3, 4)[:k]
    assert comp(*[transp(6, a[i], a[i + 1]) for i in range(k - 1)]) == cyc(6, a)
# Exemple : sigma = (1 4 9 8)(2 7 5), N = 19, p = 4, epsilon = -1.
assert inv(sig_ex) == 19 and orbites(sig_ex) == 4 and eps(sig_ex) == -1
# Trous : leurres (-1)^k, -1, (-1)^{k(k-1)/2} différents de (-1)^{k-1} pour k = 3 ; (-1)^p, (-1)^{n-p+1}, (-1)^n
# faux pour id de S_3 (p = 3) ; « constante », « additive », « eps(s^-1) = -eps(s) » faux (eps(id) = 1).
assert (-1) ** 3 != (-1) ** 2 and (-1) ** (3 * 2 // 2) != (-1) ** 2
assert eps(ident(3)) == 1 and (-1) ** 3 == -1 and (-1) ** (3 - 3 + 1) == -1
assert eps(transp(2, 1, 2)) == -1 and eps(ident(2)) != eps(ident(2)) + eps(ident(2))
# Pièges des hypothèses : (1 2 3) produit de m = 1 cycle (de longueur 3) ; (1 2)(3 4) d'ordre 2 et pair.
assert eps(cyc(3, (1, 2, 3))) == 1 and eps(cyc(4, (1, 2), (3, 4))) == 1 and ordre(cyc(4, (1, 2), (3, 4))) == 2
# Contre-exemple des descentes : (1 3) dans S_3, suite 3 2 1, 2 descentes, 3 inversions.
p = transp(3, 1, 3)
assert p == (3, 2, 1) and sum(1 for i in range(2) if p[i] > p[i + 1]) == 2 and inv(p) == 3
# Affirmations : eps(s^-1) = eps(s), conjugaison, produit de deux impaires pair.
for p in S(4):
    assert eps(inverse(p)) == eps(p)
    for t in S(4):
        assert eps(comp(t, p, inverse(t))) == eps(p)
assert eps(prodc(3, (1, 2), (2, 3))) == 1

# gr-l2-symetrique-signature-ex1 (+ v1, v2) : inversions, signature, décomposition.
for images, N, cs, leurres in [
    ((3, 1, 4, 6, 5, 2), 6, [(1, 3, 4, 6, 2)], [5, 7, 8]),
    ((2, 5, 1, 6, 4, 3), 7, [(1, 2, 5, 4, 6, 3)], [6, 8, 5]),
    ((5, 2, 6, 1, 4, 3), 9, [(1, 5, 4), (3, 6)], [8, 10, 7]),
]:
    assert inv(images) == N and N not in leurres and norm(cycles(images)) == norm(cs)
    assert eps(images) == prod((-1) ** (len(c) - 1) for c in cs)

# gr-l2-symetrique-signature-ex2 : retournement, N = n(n-1)/2, pair ssi n = 0, 1 mod 4 ; floor(n/2) transpositions.
for n in range(2, 30):
    w = tuple(n + 1 - i for i in range(1, n + 1))
    if n <= 9:
        assert inv(w) == n * (n - 1) // 2 and len(cycles(w)) == n // 2
    assert ((n * (n - 1) // 2) % 2 == 0) == (n % 4 in (0, 1)) == ((n // 2) % 2 == 0)
# Leurres « n pair », « n = 0 ou 3 mod 4 », « n = 1 ou 2 mod 4 » : chacun se trompe pour un n (n = 2, 3 ou 4).
pair_w = lambda n: (n * (n - 1) // 2) % 2 == 0
assert pair_w(2) != (2 % 2 == 0) and pair_w(3) != (3 % 4 in (0, 3)) and pair_w(4) != (4 % 4 in (1, 2))


# gr-l2-symetrique-signature-ex3 (+ v1, v2) : x -> ax mod p ; lemme de Zolotarev pour p <= 23.
def mult(a, pr):
    return tuple(a * x % pr for x in range(1, pr))


for a, pr, cs, e, leurres in [
    (2, 7, [(1, 2, 4), (3, 6, 5)], 1, [([(1, 4, 2), (3, 5, 6)], 1), ([(1, 2, 4, 3, 6, 5)], -1)]),
    (3, 7, [(1, 3, 2, 6, 4, 5)], -1, [([(1, 3, 2), (4, 5, 6)], 1), ([(1, 5, 4, 6, 2, 3)], 1)]),
    (3, 11, [(1, 3, 9, 5, 4), (2, 6, 7, 10, 8)], 1, [([(1, 3, 9, 5, 4, 2, 6, 7, 10, 8)], -1), ([(1, 4, 5, 9, 3), (2, 8, 10, 7, 6)], -1)]),
]:
    s = mult(a, pr)
    assert sorted(s) == list(range(1, pr)) and norm(cycles(s)) == norm(cs) and eps(s) == e
    for l, el in leurres:
        assert cyc(pr - 1, *l) != s or el != e
for pr in (3, 5, 7, 11, 13, 17, 19, 23):
    carres = {x * x % pr for x in range(1, pr)}
    for a in range(1, pr):
        assert eps(mult(a, pr)) == (1 if a in carres else -1)

# gr-l2-symetrique-signature-ex4 : sigma (a b) sigma^-1 = (sigma(a) sigma(b)) ; exactement 2 morphismes S_n -> {±1}.
for n in range(2, 6):
    perms = S(n)
    for s in perms:
        for a, b in combinations(range(1, n + 1), 2):
            assert comp(s, transp(n, a, b), inverse(s)) == transp(n, s[a - 1], s[b - 1])
    # Une fonction f à valeurs ±1 constante sur les transpositions et multiplicative : f = 1 ou eps.
    gen = [transp(n, k, k + 1) for k in range(1, n)]
    for v in (1, -1):
        f = {q: v ** d for q, d in distance(gen, n).items()}
        assert all(f[comp(p, q)] == f[p] * f[q] for p in perms for q in perms)
        assert all(f[q] == (1 if v == 1 else eps(q)) for q in perms)

# gr-l2-symetrique-signature-ex5 : taquin. Un coup u <-> v voisins change la couleur ; (14 15) est impaire.
assert eps(transp(16, 14, 15)) == -1
case = lambda c: ((c - 1) // 4 + 1, (c - 1) % 4 + 1)
for u in range(1, 17):
    for v in range(1, 17):
        (i1, j1), (i2, j2) = case(u), case(v)
        if abs(i1 - i2) + abs(j1 - j2) == 1:
            assert (i1 + j1) % 2 != (i2 + j2) % 2
pi = list(ident(16))
pi[15 - 1], pi[14 - 1] = 14, 15
assert tuple(pi) == transp(16, 14, 15)
# Le coup : pi' = pi ∘ (u v) échange les contenus des cases u et v.
pi0 = (3, 1, 16, 2) + tuple(range(5, 17))[:-1] + (4,)
assert sorted(pi0) == list(range(1, 17))
pi1 = comp(pi0, transp(16, 3, 4))
assert pi1[2] == pi0[3] and pi1[3] == pi0[2] and all(pi1[k] == pi0[k] for k in range(16) if k not in (2, 3))

# ---------------------------------------------------------------------------
# gr-l2-symetrique-transpositions (carte)
# ---------------------------------------------------------------------------
for n in range(2, 7):
    toutes = [transp(n, a, b) for a, b in combinations(range(1, n + 1), 2)]
    d = distance(toutes, n)
    assert len(d) == factorial(n) and all(d[p] == n - orbites(p) for p in d)  # n - p transpositions, minimal
    assert max(d.values()) == n - 1 and d[cyc(n, tuple(range(1, n + 1)))] == n - 1
    assert len(engendre([transp(n, 1, k) for k in range(2, n + 1)])) == factorial(n)
    assert len(engendre([transp(n, k, k + 1) for k in range(1, n)])) == factorial(n)
    if n >= 3:
        # Trou 1 : n - 2, floor(n/2), n faux (un n-cycle demande exactement n - 1).
        assert n - 1 not in (n - 2, n // 2, n) or n == 2
        # Trou 2 et trou 3 : leurres qui n'engendrent pas S_n.
        for fam in ([transp(n, 1, k) for k in range(2, n)], [transp(n, 2, k) for k in range(3, n + 1)],
                    [transp(n, 1, k) for k in range(2, n + 1) if k % 2 == 0],
                    [transp(n, k, k + 1) for k in range(1, n - 1)],
                    [transp(n, k, k + 1) for k in range(1, n) if k % 2 == 1],
                    [transp(n, k, k + 2) for k in range(1, n - 1)]):
            assert len(engendre(fam)) < factorial(n)
# Formules de la démonstration : (a b) = (1 a)(1 b)(1 a) ; (i j+1) = (j j+1)(i j)(j j+1).
for n in range(3, 7):
    for a, b in permutations(range(2, n + 1), 2):
        assert comp(transp(n, 1, a), transp(n, 1, b), transp(n, 1, a)) == transp(n, a, b)
    for i in range(1, n):
        for j in range(i + 1, n):
            assert comp(transp(n, j, j + 1), transp(n, i, j), transp(n, j, j + 1)) == transp(n, i, j + 1)
# Coupure et fusion d'orbites : le nombre d'orbites de sigma tau vaut p ± 1.
for p in S(5):
    for a, b in combinations(range(1, 6), 2):
        assert abs(orbites(comp(p, transp(5, a, b))) - orbites(p)) == 1
# Exemple : (1 3 5)(2 4) = (1 3)(3 5)(2 4), 3 = 5 - 2 transpositions ; (1 3) = (2 3)(1 2)(2 3).
assert prodc(5, (1, 3), (3, 5), (2, 4)) == cyc(5, (1, 3, 5), (2, 4)) and orbites(cyc(5, (1, 3, 5), (2, 4))) == 2
assert prodc(3, (2, 3), (1, 2), (2, 3)) == transp(3, 1, 3)
# Contre-exemples et affirmations : (1 2), (3 4) n'engendrent pas S_4 ; (1 2), (2 3) fixent 4 ; parité de (1 2).
assert len(engendre([transp(4, 1, 2), transp(4, 3, 4)])) == 4 and len(engendre([transp(4, 1, 2), transp(4, 2, 3)])) == 6
assert comp(*[transp(2, 1, 2)] * 5) == transp(2, 1, 2) and eps(transp(2, 1, 2)) == -1

# gr-l2-symetrique-transpositions-ex1 (+ v1, v2) : nombre minimal n - p et décompositions données.
T7 = lambda n: [transp(n, a, b) for a, b in combinations(range(1, n + 1), 2)]
for n, cs, bonne, leurres, ecr in [
    (7, [(1, 4, 2, 6), (3, 7)], 4, [5, 3, 6], [(1, 4), (4, 2), (2, 6), (3, 7)]),
    (8, [(1, 5), (2, 3, 6, 4, 7)], 5, [6, 4, 7], [(1, 5), (2, 3), (3, 6), (6, 4), (4, 7)]),
    (9, [(1, 2, 3), (4, 5, 6), (7, 8, 9)], 6, [3, 9, 7], [(1, 2), (2, 3), (4, 5), (5, 6), (7, 8), (8, 9)]),
]:
    s = cyc(n, *cs)
    assert n - orbites(s) == bonne and bonne not in leurres and prodc(n, *ecr) == s and len(ecr) == bonne
    assert n - (orbites(s) - 1) in leurres  # oublier le point fixe (ou compter un par cycle pour v2 : 3)
# Minimalité vérifiée par recherche exhaustive pour ex1 (S_7) : aucune écriture en 3 transpositions.
s = cyc(7, (1, 4, 2, 6), (3, 7))
assert all(comp(*w) != s for w in product(T7(7), repeat=3))

# gr-l2-symetrique-transpositions-ex2 (+ v1, v2) : écritures en transpositions voisines et leurres.
for n, cible, bonne, leurres in [
    (4, (1, 4), [(3, 4), (2, 3), (1, 2), (2, 3), (3, 4)], [[(1, 2), (2, 3), (3, 4)], [(1, 2), (3, 4), (1, 2)], [(2, 3), (1, 2), (3, 4), (2, 3)]]),
    (5, (2, 5), [(4, 5), (3, 4), (2, 3), (3, 4), (4, 5)], [[(2, 3), (3, 4), (4, 5)], [(2, 3), (4, 5), (2, 3)], [(3, 4), (2, 3), (4, 5), (3, 4)]]),
    (5, (1, 5), [(4, 5), (3, 4), (2, 3), (1, 2), (2, 3), (3, 4), (4, 5)],
     [[(1, 2), (2, 3), (3, 4), (4, 5)], [(4, 5), (3, 4), (1, 2), (3, 4), (4, 5)], [(1, 2), (2, 3), (4, 5), (2, 3), (1, 2)]]),
]:
    assert prodc(n, *bonne) == transp(n, *cible)
    for l in leurres:
        assert prodc(n, *l) != transp(n, *cible)
assert prodc(4, (1, 2), (2, 3), (3, 4)) == cyc(4, (1, 2, 3, 4)) and prodc(4, (1, 2), (3, 4), (1, 2)) == transp(4, 3, 4)
assert prodc(4, (2, 3), (1, 2), (3, 4), (2, 3)) == cyc(4, (1, 3), (2, 4))
assert prodc(5, (3, 4), (2, 3), (4, 5), (3, 4)) == cyc(5, (2, 4), (3, 5))
assert prodc(5, (4, 5), (3, 4), (1, 2), (3, 4), (4, 5)) == transp(5, 1, 2)
assert prodc(5, (1, 2), (2, 3), (4, 5), (2, 3), (1, 2)) == transp(5, 4, 5)

# gr-l2-symetrique-transpositions-ex3 (+ v1, v2) : nombre minimal de transpositions voisines = N(sigma).
for n in range(2, 7):
    d = distance([transp(n, k, k + 1) for k in range(1, n)], n)
    assert all(d[p] == inv(p) for p in d)
    for p in S(n):
        for k in range(1, n):
            assert inv(comp(p, transp(n, k, k + 1))) - inv(p) == (1 if p[k - 1] < p[k] else -1)
for images, N, cs, m, leurres in [
    ((4, 3, 1, 2), 5, [(1, 4, 2, 3)], 3, [3, 4, 6]),
    ((3, 5, 1, 4, 2), 6, [(1, 3), (2, 5)], 2, [2, 5, 7]),
    ((5, 4, 3, 2, 1), 10, [(1, 5), (2, 4)], 2, [2, 5, 8]),
]:
    n = len(images)
    assert inv(images) == N and norm(cycles(images)) == norm(cs) and n - orbites(images) == m and N not in leurres
assert prodc(4, (1, 4), (4, 2), (2, 3)) == (4, 3, 1, 2)

# gr-l2-symetrique-transpositions-ex4 : c^k (1 2) c^-k = (k+1 k+2) ; <(1 2), c> = S_n ; <(1 3), (1 2 3 4)> d'ordre 8.
for n in range(2, 8):
    c = cyc(n, tuple(range(1, n + 1)))
    for k in range(n - 1):
        assert comp(puiss(c, k), transp(n, 1, 2), puiss(c, n - k)) == transp(n, k + 1, k + 2)
    if n <= 6:
        assert len(engendre([transp(n, 1, 2), c])) == factorial(n)
c4, t = cyc(4, (1, 2, 3, 4)), transp(4, 1, 3)
assert comp(t, c4, inverse(t)) == cyc(4, (3, 2, 1, 4)) == inverse(c4) and len(engendre([t, c4])) == 8

# gr-l2-symetrique-transpositions-ex5 : T engendre S_n ssi le graphe est connexe ; minimum n - 1 (n = 3, 4, 5).


def connexe(n, aretes):
    vus, pile = {1}, [1]
    while pile:
        x = pile.pop()
        for a, b in aretes:
            for u, v in ((a, b), (b, a)):
                if u == x and v not in vus:
                    vus.add(v)
                    pile.append(v)
    return len(vus) == n


for n in (3, 4, 5):
    paires = list(combinations(range(1, n + 1), 2))
    mini = None
    for r in range(1, len(paires) + 1):
        for E in combinations(paires, r):
            if n == 5 and r > 5:
                break
            gen = len(engendre([transp(n, a, b) for a, b in E])) == factorial(n)
            assert gen == connexe(n, E)
            if gen and mini is None:
                mini = r
    assert mini == n - 1 and n - 1 not in (n, n * (n - 1) // 2, n // 2)
for a, b, c in permutations(range(1, 6), 3):
    assert comp(transp(5, b, c), transp(5, a, b), transp(5, b, c)) == transp(5, a, c)

# ---------------------------------------------------------------------------
# gr-l2-symetrique-alterne (carte)
# ---------------------------------------------------------------------------
for n in range(2, 7):
    A = {p for p in S(n) if eps(p) == 1}
    assert len(A) == factorial(n) // 2 and engendre(list(A)) == A
    t = transp(n, 1, 2)
    assert {comp(t, a) for a in A} == set(S(n)) - A
    if n <= 5:
        assert all({comp(s, a, inverse(s)) for a in A} == A for s in S(n))
    if n >= 3:
        trois = [p for p in S(n) if type_(p) == (3,)]
        assert engendre(trois) == A
        # Trou 2 : (n-1)!, 2^(n-1), n(n-1)/2 différents de n!/2 pour n = 4.
        if n == 4:
            assert 12 not in (factorial(3), 2 ** 3, 6)
        # Trou 3 : transpositions, 4-cycles, permutations d'ordre 2 n'engendrent pas A_n.
        for fam in ([p for p in S(n) if type_(p) == (2,)], [p for p in S(n) if type_(p) == (4,)],
                    [p for p in S(n) if ordre(p) == 2]):
            assert (engendre(fam) if fam else {ident(n)}) != A
    # Trou 1 : permutations d'ordre 2, permutations sans point fixe : différentes de A_n (id manque).
    assert ident(n) in A and ordre(ident(n)) != 2 and any(x == i for i, x in enumerate(ident(n), 1))
# Cycles pairs ssi longueur impaire.
for k in range(2, 8):
    assert (eps(cyc(8, tuple(range(1, k + 1)))) == 1) == (k % 2 == 1)
# Démonstration : (a b)(a c) = (a c b) ; (a b)(c d) = (a b c)(b c d) ; (x y)(y z) = (x y z).
for a, b, c, d in permutations(range(1, 6), 4):
    assert comp(transp(5, a, b), transp(5, a, c)) == cyc(5, (a, c, b))
    assert comp(transp(5, a, b), transp(5, c, d)) == comp(cyc(5, (a, b, c)), cyc(5, (b, c, d)))
    assert comp(transp(5, a, b), transp(5, b, c)) == cyc(5, (a, b, c))
# Exemple et contre-exemples : A_3, A_4, (1 2)(3 4) = (1 2 3)(2 3 4), (1 2)(3 4 5) impaire.
assert {p for p in S(3) if eps(p) == 1} == {ident(3), cyc(3, (1, 2, 3)), cyc(3, (1, 3, 2))}
A4 = [p for p in S(4) if eps(p) == 1]
assert sorted(type_(p) for p in A4).count((3,)) == 8 and sorted(type_(p) for p in A4).count((2, 2)) == 3
assert cyc(4, (1, 2), (3, 4)) == comp(cyc(4, (1, 2, 3)), cyc(4, (2, 3, 4)))
assert eps(cyc(5, (1, 2), (3, 4, 5))) == -1 and eps(transp(2, 1, 2)) == -1
# Affirmations : (1 2 3)(1 2 4) = (1 3)(2 4), (1 2 4)(1 2 3) = (1 4)(2 3).
assert prodc(4, (1, 2, 3), (1, 2, 4)) == cyc(4, (1, 3), (2, 4)) and prodc(4, (1, 2, 4), (1, 2, 3)) == cyc(4, (1, 4), (2, 3))

# gr-l2-symetrique-alterne-ex1 (+ v1, v2) : une seule permutation impaire dans chaque liste (la première).
for n, liste in [
    (6, [[(1, 2, 3), (4, 5)], [(1, 2, 3, 4), (5, 6)], [(1, 2, 3, 4, 5)], [(1, 2), (3, 4)]]),
    (7, [[(1, 2, 3, 4, 5, 6)], [(1, 2, 3), (4, 5, 6)], [(1, 2), (3, 4, 5, 6)], [(1, 2, 3, 4, 5, 6, 7)]]),
    (8, [[(1, 2), (3, 4), (5, 6)], [(1, 2, 3, 4), (5, 6, 7, 8)], [(1, 2, 3, 4, 5, 6, 7)], [(1, 2, 3), (4, 5, 6, 7, 8)]]),
]:
    assert [eps(cyc(n, *cs)) for cs in liste] == [-1, 1, 1, 1]

# gr-l2-symetrique-alterne-ex2 (+ v1, v2) : effectifs des types dans A_4, A_5, A_6.
for n, attendu in [(4, {(): 1, (3,): 8, (2, 2): 3}),
                   (5, {(): 1, (3,): 20, (2, 2): 15, (5,): 24}),
                   (6, {(): 1, (3,): 40, (2, 2): 45, (5,): 144, (3, 3): 40, (4, 2): 90})]:
    cpt = {}
    for p in S(n):
        if eps(p) == 1:
            cpt[type_(p)] = cpt.get(type_(p), 0) + 1
    assert cpt == attendu and sum(attendu.values()) == factorial(n) // 2
    for t, v in attendu.items():
        if t:
            assert formule(n, t) == v
# Types de S_5 et types impairs de S_6 cités dans les corrections.
assert [nb_type(5, t) for t in [(2,), (3,), (4,), (5,), (2, 2), (3, 2)]] == [10, 20, 30, 24, 15, 20]
assert [nb_type(6, t) for t in [(2,), (4,), (6,), (3, 2), (2, 2, 2)]] == [15, 90, 120, 120, 15]

# gr-l2-symetrique-alterne-ex3 (+ v1, v2) : décompositions en 3-cycles et leurres.
for n, cible, bonne, leurres in [
    (7, [(1, 2, 3, 4, 5, 6, 7)], [(1, 2, 3), (3, 4, 5), (5, 6, 7)], [[(5, 6, 7), (3, 4, 5), (1, 2, 3)], [(1, 2, 3), (4, 5, 6), (6, 7, 1)], [(1, 3, 5), (5, 7, 1)]]),
    (4, [(1, 2), (3, 4)], [(1, 2, 3), (2, 3, 4)], [[(2, 3, 4), (1, 2, 3)], [(1, 2, 3), (1, 3, 4)], [(1, 2, 4), (1, 3, 4)]]),
    (6, [(1, 2, 3, 4), (5, 6)], [(1, 2, 3), (3, 4, 5), (4, 5, 6)], [[(4, 5, 6), (3, 4, 5), (1, 2, 3)], [(1, 2, 3), (4, 5, 6)], [(1, 2, 3), (3, 4, 5), (4, 6, 5)]]),
]:
    assert prodc(n, *bonne) == cyc(n, *cible)
    for l in leurres:
        assert prodc(n, *l) != cyc(n, *cible)
assert prodc(7, (5, 6, 7), (3, 4, 5), (1, 2, 3)) == cyc(7, (1, 2, 4, 6, 7, 5, 3))
assert prodc(4, (2, 3, 4), (1, 2, 3)) == cyc(4, (1, 3), (2, 4))
assert prodc(6, (4, 5, 6), (3, 4, 5), (1, 2, 3)) == cyc(6, (1, 2, 5, 3), (4, 6))

# gr-l2-symetrique-alterne-ex4 : g_k^2 = (1 k 2), (1 a b) = g_b g_a^2, (2 a b) = g_b^2 g_a, (a b c) = (1 a b)(1 b c).
n = 7
g = lambda k: cyc(n, (1, 2, k))
for k in range(3, n + 1):
    assert comp(g(k), g(k)) == cyc(n, (1, k, 2))
for a, b in permutations(range(3, n + 1), 2):
    assert comp(g(b), g(a), g(a)) == cyc(n, (1, a, b)) and comp(g(b), g(b), g(a)) == cyc(n, (2, a, b))
for a, b, c in permutations(range(2, n + 1), 3):
    assert comp(cyc(n, (1, a, b)), cyc(n, (1, b, c))) == cyc(n, (a, b, c))
for m in range(3, 7):
    assert engendre([cyc(m, (1, 2, k)) for k in range(3, m + 1)]) == {p for p in S(m) if eps(p) == 1}

# gr-l2-symetrique-alterne-ex5 : carrés dans S_n engendrent A_n (n = 3..5) ; tout 3-cycle est un carré ;
# A_4 n'a pas de sous-groupe d'ordre 6 (recherche exhaustive).
for m in range(3, 6):
    carres = {comp(p, p) for p in S(m)}
    assert engendre(list(carres)) == {p for p in S(m) if eps(p) == 1}
for c in [p for p in S(4) if type_(p) == (3,)]:
    assert comp(comp(c, c), comp(c, c)) == c
autres = [p for p in A4 if p != ident(4)]
for E in combinations(autres, 5):
    H = set(E) | {ident(4)}
    assert not all(comp(x, y) in H for x in H for y in H)
# Sous-groupes de cardinal n!/2 de S_4 : un tel sous-groupe contient les carrés, donc A_4 (vérifié ci-dessus).

# ---------------------------------------------------------------------------
# gr-l2-symetrique-ordre (carte)
# ---------------------------------------------------------------------------
for n in range(1, 8):
    L = lcm(*range(1, n + 1))
    for p in S(n):
        o = next(k for k in range(1, 500) if puiss(p, k) == ident(n)) if p != ident(n) else 1
        assert o == ordre(p) and L % o == 0 and factorial(n) % o == 0
# Exemple : (1 6 3)(2 5 8 4)(7 9), ordre 12, sigma^100 = sigma^4 = (1 6 3).
s = cyc(9, (1, 6, 3), (2, 5, 8, 4), (7, 9))
assert ordre(s) == 12 != 24 and 100 % 12 == 4 and puiss(s, 100) == puiss(s, 4) == cyc(9, (1, 6, 3))
# Contre-exemples et pièges : (1 2)(2 3) d'ordre 3 ; (1 2)(3 4) d'ordre 2 ; c c^2 = id.
assert ordre(prodc(3, (1, 2), (2, 3))) == 3 and ordre(cyc(4, (1, 2), (3, 4))) == 2
c = cyc(3, (1, 2, 3))
assert comp(c, comp(c, c)) == ident(3)
# Trous : produit, somme, max des longueurs faux pour (1 2)(3 4) ou (1 2)(3 4 5) ; n, 2^n, n(n-1) faux.
assert ordre(cyc(4, (1, 2), (3, 4))) not in (4,) and ordre(cyc(5, (1, 2), (3, 4, 5))) != 3
assert 3 % 4 != 0 and 8 % 3 != 0 and 20 % 6 != 0  # S_4 : ordre 3 ∤ 4 ; S_3 : 3 ∤ 8 ; S_5 : 6 ∤ 20
# Affirmations : (1 2)(3 4 5) d'ordre 6 ; ordre premier => cycles de longueur p ; (1 2)^2 = id.
assert ordre(cyc(5, (1, 2), (3, 4, 5))) == 6
for p in S(6):
    if ordre(p) in (2, 3, 5):
        assert set(type_(p)) == {ordre(p)}
# Cours : plus grand ordre 4, 6, 6, 12 pour n = 4, 5, 6, 7.
assert [max(ordre(p) for p in S(n)) for n in (4, 5, 6, 7)] == [4, 6, 6, 12]

# gr-l2-symetrique-ordre-ex1 (+ v1, v2) : matrices, décompositions et ordres.
for images, cs, o, leurres in [
    ((6, 7, 9, 2, 8, 3, 4, 5, 1), [(1, 6, 3, 9), (2, 7, 4), (5, 8)], 12, [24, 9, 4]),
    ((4, 8, 6, 7, 3, 5, 2, 1), [(1, 4, 7, 2, 8), (3, 6, 5)], 15, [8, 5, 30]),
    ((3, 1, 5, 10, 7, 8, 9, 4, 2, 6), [(1, 3, 5, 7, 9, 2), (4, 10, 6, 8)], 12, [24, 10, 6]),
]:
    assert images == cyc(len(images), *cs) and ordre(images) == o and o not in leurres

# gr-l2-symetrique-ordre-ex2 (+ v1, v2) : grandes puissances.
for n, cs, m, bonne, leurres in [
    (7, [(1, 2, 3), (4, 5, 6, 7)], 2026, [(1, 2, 3), (4, 6), (5, 7)], [[(1, 3, 2), (4, 6), (5, 7)], [(1, 2, 3), (4, 7, 6, 5)], [(1, 2, 3), (4, 5, 6, 7)]]),
    (7, [(1, 2, 3, 4, 5), (6, 7)], 2027, [(1, 3, 5, 2, 4), (6, 7)], [[(1, 3, 5, 2, 4)], [(1, 2, 3, 4, 5), (6, 7)], [(1, 4, 2, 5, 3), (6, 7)]]),
    (10, [(1, 2, 3, 4), (5, 6, 7, 8, 9, 10)], 1000, [(5, 9, 7), (6, 10, 8)], [[(1, 3), (2, 4), (5, 9, 7), (6, 10, 8)], [(5, 7, 9), (6, 8, 10)], []]),
]:
    s = cyc(n, *cs)
    r = puiss(s, m)
    assert r == cyc(n, *bonne)
    for l in leurres:
        assert cyc(n, *l) != r
assert 2026 == 3 * 675 + 1 == 4 * 506 + 2 == 12 * 168 + 10 and 2027 == 5 * 405 + 2 and 1000 == 6 * 166 + 4


# gr-l2-symetrique-ordre-ex3 (+ v1, v2) : ordre maximal (fonction de Landau) par partitions.
def partitions(n, maxi=None):
    if maxi is None:
        maxi = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxi), 0, -1):
        for r in partitions(n - k, k):
            yield (k,) + r


def landau(n):
    return max(lcm(*p) for p in partitions(n))


ordres = lambda n: {lcm(*p) for p in partitions(n)}
assert landau(7) == 12 and 10 in ordres(7) and 21 not in ordres(7)
assert landau(9) == 20 and 18 not in ordres(9) and {15, 14} <= ordres(9)
assert landau(10) == 30 and 21 in ordres(10) and 24 not in ordres(10) and 20 in ordres(10)
assert ordre(cyc(7, (1, 2, 3, 4), (5, 6, 7))) == 12 and ordre(cyc(9, (1, 2, 3, 4), (5, 6, 7, 8, 9))) == 20
assert ordre(cyc(10, (1, 2), (3, 4, 5), (6, 7, 8, 9, 10))) == 30
assert max(ordre(p) for p in S(7)) == 12

# gr-l2-symetrique-ordre-ex4 (+ v1, v2) : éléments d'ordre 2 de S_5, S_6, S_4.
for n, bonne, leurres in [(5, 25, [10, 15, 30]), (6, 75, [60, 90, 45]), (4, 9, [6, 10, 12])]:
    assert sum(1 for p in S(n) if ordre(p) == 2) == bonne and bonne not in leurres
assert comb(5, 2) == 10 and comb(5, 2) * comb(3, 2) // 2 == 15 and comb(6, 2) * comb(4, 2) // 2 == 45
assert comb(6, 2) * comb(4, 2) * comb(2, 2) // 6 == 15


# gr-l2-symetrique-ordre-ex5 (+ v1, v2) : plus petit n tel que S_n contienne un élément d'ordre m.
def plus_petit(m):
    return next(n for n in range(1, 60) if m in ordres(n))


for m, bonne, leurres in [(12, 7, [12, 6, 5]), (30, 10, [30, 8, 11]), (18, 11, [9, 18, 8])]:
    assert plus_petit(m) == bonne and bonne not in leurres
assert lcm(3, 3, 2) == 6 and ordre(cyc(11, (1, 2), tuple(range(3, 12)))) == 18
assert all(a * b >= a + b for a, b in [(2, 3), (2, 5), (3, 5)]) and 2 * 3 * 5 >= 10 and 4 * 3 >= 7 and 9 * 2 >= 11

# ---------------------------------------------------------------------------
# QCM gr-l2-symetrique-q1 … q8
# ---------------------------------------------------------------------------
assert nb_type(5, (3,)) == 20 and 5 * 4 * 3 == 60  # q1
assert eps(prodc(3, (1, 2), (2, 3))) == 1  # q2
assert [eps(cyc(5, *cs)) for cs in ([(1, 2, 3, 4, 5)], [(1, 2, 3, 4)], [(1, 2), (3, 4, 5)], [(1, 2)])] == [1, -1, -1, -1]  # q3
assert len(engendre([transp(4, 1, 3), cyc(4, (1, 2, 3, 4))])) == 8 and len(engendre([transp(4, 1, 2), cyc(4, (1, 2, 3, 4))])) == 24  # q4
assert sum(1 for p in A4 if ordre(p) == 2) == 3  # q5
assert prodc(5, (1, 2, 3), (3, 4, 5)) == cyc(5, (1, 2, 3, 4, 5))  # q6
assert ordre(cyc(9, (1, 2, 3), (4, 5), (6, 7, 8, 9))) == 12  # q7
assert not any(ordre(p) == 6 for p in S(6) if eps(p) == 1) and any(ordre(p) == 6 for p in S(6))  # q8

# ---------------------------------------------------------------------------
# Relecture indépendante (2026-10-29) : faits des corrections et des cartes non encore contrôlés
# ---------------------------------------------------------------------------
# Produits « dans l'autre ordre » cités par les corrections :
# gr-l2-symetrique-cycles-ex2 (+ v1, v2), gr-l2-symetrique-alterne-ex3 (+ v1, v2).
for n, gauche, droite in [
    (5, [(3, 4, 5), (1, 2, 3, 4)], [(1, 2, 4), (3, 5)]),
    (6, [(2, 3, 4, 6), (1, 3, 5)], [(1, 4, 6, 2, 3, 5)]),
    (5, [(2, 5, 3, 4), (1, 4, 2)], [(1, 2), (3, 4, 5)]),
    (7, [(5, 6, 7), (3, 4, 5), (1, 2, 3)], [(1, 2, 4, 6, 7, 5, 3)]),
    (4, [(2, 3, 4), (1, 2, 3)], [(1, 3), (2, 4)]),
    (6, [(4, 5, 6), (3, 4, 5), (1, 2, 3)], [(1, 2, 5, 3), (4, 6)]),
]:
    assert prodc(n, *gauche) == cyc(n, *droite)
# gr-l2-symetrique-transpositions-ex4 (remarque) et q4 : t c t^-1 = c^-1 pour t = (1 3), c = (1 2 3 4).
t13, c4 = transp(4, 1, 3), cyc(4, (1, 2, 3, 4))
assert comp(t13, c4, inverse(t13)) == inverse(c4) == cyc(4, (3, 2, 1, 4))
# gr-l2-symetrique-signature (exemple) : la suite 4,7,3,9,2,6,5,1,8 a 19 inversions et 4 orbites.
ex_sig = (4, 7, 3, 9, 2, 6, 5, 1, 8)
assert ex_sig == cyc(9, (1, 4, 9, 8), (2, 7, 5)) and inv(ex_sig) == 19 and orbites(ex_sig) == 4
# gr-l2-symetrique-ordre (cours) : plus grand ordre 4, 6, 6, 12 pour n = 4, 5, 6, 7.
assert [landau(n) for n in (4, 5, 6, 7)] == [4, 6, 6, 12]
# gr-l2-symetrique-signature-ex5 (taquin) : un coup (case vide de u vers v voisine) change pi en pi o (u v),
# change la couleur (parité de i + j) de la case vide ; la position visée (14 15) est impaire.
def case(c):
    return (c - 1) // 4 + 1, (c - 1) % 4 + 1
for u in range(1, 17):
    for v in range(1, 17):
        (i, j), (k, l) = case(u), case(v)
        if abs(i - k) + abs(j - l) == 1:
            assert (i + j) % 2 != (k + l) % 2
            pi = tuple(range(1, 17))
            pi = comp(pi, transp(16, u, 16)) if u != 16 else pi  # une position quelconque où u est vide
            nouv = list(pi)
            nouv[u - 1], nouv[v - 1] = pi[v - 1], 16
            assert pi[u - 1] == 16 and tuple(nouv) == comp(pi, transp(16, u, v))
assert eps(transp(16, 14, 15)) == -1
# Trous : chaque leurre rend l'énoncé faux pour au moins un n (énoncés universels).
# gr-l2-symetrique-signature trou 3 : (-1)^p faux pour n = 3 (id : p = 3).
assert (-1) ** orbites(ident(3)) != eps(ident(3))
# gr-l2-symetrique-alterne trou 2 : (n-1)!, 2^(n-1), n(n-1)/2 faux pour n = 4 (|A_4| = 12).
assert 12 not in (factorial(3), 2 ** 3, 4 * 3 // 2)
# gr-l2-symetrique-ordre trou 3 : n, 2^n, n(n-1) ne sont pas multiples de l'ordre 6 de (1 2)(3 4 5) dans S_5.
assert all(x % 6 for x in (5, 2 ** 5, 5 * 4))

print("gr_l2_symetrique : toutes les vérifications passent.")
