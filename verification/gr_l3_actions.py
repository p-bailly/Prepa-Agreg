"""Vérifications calculatoires du chapitre « gr-l3-actions » (groupes, L3) : cartes, exercices et QCM de
content/*/groupes/l3-gr-l3-actions.json et content/qcm/gr-l3-actions.json.

Convention du chapitre : sigma tau = sigma ∘ tau (on applique d'abord tau). Chaque assert cite l'identifiant de la carte,
de l'exercice (variante comprise : v1, v2) ou de la question. Tout est calculé par force brute (énumération des groupes,
des orbites, des stabilisateurs, des classes de conjugaison, des coloriages), en entiers ou en fractions exactes, sauf
le groupe du dodécaèdre (matrices flottantes arrondies). Lancer : python3 verification/gr_l3_actions.py.
"""

from collections import deque
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, factorial, gcd, lcm

import numpy as np

# ---------------------------------------------------------------------------
# Outils : une permutation de {1..n} est un tuple p avec p[i-1] = image de i.
# ---------------------------------------------------------------------------


def ident(n):
    return tuple(range(1, n + 1))


def cyc(n, *cs):
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


def inverse(p):
    q = [0] * len(p)
    for i, x in enumerate(p, 1):
        q[x - 1] = i
    return tuple(q)


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
    out = []
    for c in cs:
        k = c.index(min(c))
        out.append(tuple(c[k:] + c[:k]))
    return sorted(out)


def type_(p):
    return tuple(sorted((len(c) for c in cycles(p)), reverse=True))


def ordre(p):
    cs = cycles(p)
    return lcm(*[len(c) for c in cs]) if cs else 1


def eps(p):
    n = len(p)
    return (-1) ** sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])


def S(n):
    return [tuple(p) for p in permutations(range(1, n + 1))]


def A(n):
    return [p for p in S(n) if eps(p) == 1]


def engendre(gens):
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


def conj(t, s):
    return comp(t, s, inverse(t))


def classes(G):
    G = list(G)
    reste, res = set(G), []
    while reste:
        x = next(iter(reste))
        cl = {conj(g, x) for g in G}
        res.append(cl)
        reste -= cl
    return res


def centre(G):
    G = list(G)
    return [z for z in G if all(comp(z, g) == comp(g, z) for g in G)]


def centralisateur(G, x):
    return [g for g in G if comp(g, x) == comp(x, g)]


def orbites_action(G, X, act):
    """Orbites (liste de frozensets) d'une action donnée par act(g, x)."""
    reste, res = set(X), []
    while reste:
        x = next(iter(reste))
        o = frozenset(act(g, x) for g in G)
        res.append(o)
        reste -= o
    return res


def stab(G, x, act):
    return [g for g in G if act(g, x) == x]


perm_act = lambda g, x: g[x - 1]  # noqa: E731


def nb_orbites_burnside(G, X, act):
    return F(sum(sum(1 for x in X if act(g, x) == x) for g in G), len(G))


S3, S4, S5, S6 = S(3), S(4), S(5), S(6)
A4, A5 = A(4), A(5)
D = engendre([cyc(4, [1, 2, 3, 4]), cyc(4, [1, 3])])  # groupe à 8 éléments (D_4 dans S_4)
assert len(D) == 8

# ===========================================================================
# Carte gr-l3-actions-orbites
# ===========================================================================
# énoncé : orbites = partition, stabilisateurs sous-groupes conjugués, |orbite| = |G|/|G_x| (S4 sur les paires, S3, D)
for G, n in [(S3, 3), (S4, 4), (list(D), 4), (A4, 4)]:
    X = range(1, n + 1)
    for x in X:
        Gx = stab(G, x, perm_act)
        o = {perm_act(g, x) for g in G}
        assert len(o) * len(Gx) == len(G)
        assert all(comp(a, inverse(b)) in Gx for a in Gx for b in Gx)  # sous-groupe
        for g in G:  # G_{g.x} = g G_x g^-1
            assert set(stab(G, perm_act(g, x), perm_act)) == {conj(g, h) for h in Gx}
        for g in G:  # bijection gG_x -> g.x : g.x = g'.x <=> g^-1 g' in G_x
            for h in G:
                assert (perm_act(g, x) == perm_act(h, x)) == (comp(inverse(g), h) in Gx)
# exemple : S4 sur les paires, orbite de {1,2} = 6 paires, stabilisateur à 4 éléments
paire_act = lambda g, P: frozenset(g[a - 1] for a in P)  # noqa: E731
P4 = [frozenset(c) for c in [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]]
assert {paire_act(g, frozenset({1, 2})) for g in S4} == set(P4)
st = stab(S4, frozenset({1, 2}), paire_act)
assert set(st) == {ident(4), cyc(4, [1, 2]), cyc(4, [3, 4]), cyc(4, [1, 2], [3, 4])}
# contre-exemple : G_3 = {id,(1 2)} non distingué, G_1 != G_3
G3 = stab(S3, 3, perm_act)
assert set(G3) == {ident(3), cyc(3, [1, 2])}
assert conj(cyc(3, [1, 3]), cyc(3, [1, 2])) == cyc(3, [2, 3]) and cyc(3, [2, 3]) not in G3
assert set(stab(S3, 1, perm_act)) == {ident(3), cyc(3, [2, 3])}
# piège 1 / contre-exemple : P = {id,(1 2),(2 3)} n'est pas un groupe ; « orbite » de 1 à 2 éléments, 2 fixant 1, 3/2 != 2
P = [ident(3), cyc(3, [1, 2]), cyc(3, [2, 3])]
assert {g[0] for g in P} == {1, 2} and sum(1 for g in P if g[0] == 1) == 2 and F(3, 2) != 2
assert comp(cyc(3, [1, 2]), cyc(3, [2, 3])) not in P
# piège 2 : g.x = a (constant) vérifie g.(h.x) = (gh).x, mais {g : g.b = b} est vide (pas un sous-groupe)
cst = lambda g, x: "a"  # noqa: E731
assert all(cst(g, cst(h, x)) == cst(comp(g, h), x) for g in S3 for h in S3 for x in "ab")
assert [g for g in S3 if cst(g, "b") == "b"] == []
# piège 3 : Z/3 sur {a, b}, 0 et 2 agissent par l'identité, 1 échange : e.x = x, mais {k : k.a = a} = {0, 2} n'est pas un sous-groupe
act3 = lambda k, x: x if k in (0, 2) else {"a": "b", "b": "a"}[x]  # noqa: E731
Sa = [k for k in range(3) if act3(k, "a") == "a"]
assert Sa == [0, 2] and (2 + 2) % 3 not in Sa
assert act3((1 + 2) % 3, "a") != act3(1, act3(2, "a"))  # la condition g.(h.x) = (gh).x est bien violée
# trous : leurres faux
G4 = stab(S4, 4, perm_act)
assert any(comp(a, b) != comp(b, a) for a in G4 for b in G4)  # « sous-groupe commutatif » faux (S4 sur {1..4})
assert len(G3) > 1  # « réduit à {e} » faux
g, h = cyc(3, [1, 3]), cyc(3, [1, 2])  # gG_3 = (gh)G_3, mais g^-1.3 != (gh)^-1.3 : g^-1 . x mal défini
assert comp(inverse(g), comp(g, h)) in G3 and perm_act(inverse(g), 3) != perm_act(inverse(comp(g, h)), 3)
# « G_x g -> g.x » mal défini : G_3 g = G_3 g' sans g.3 = g'.3
assert any(comp(g2, inverse(g1)) in G3 and perm_act(g1, 3) != perm_act(g2, 3) for g1 in S3 for g2 in S3)
o3 = len({perm_act(g, 3) for g in S3})
assert o3 == 3 and o3 != len(G3) and o3 != len(S3) - len(G3) and F(3, len(G3)) != o3  # leurres |G_x|, |G|-|G_x|, |X|/|G_x|
# affirmations : action fidèle sans stabilisateurs triviaux ; orbites de cardinaux différents
assert all(any(perm_act(g, x) != x for x in range(1, 4)) for g in S3 if g != ident(3))
assert sorted(len(o) for o in orbites_action([ident(3), cyc(3, [1, 2])], range(1, 4), perm_act)) == [1, 2]

# --- gr-l3-actions-orbites-ex1 (et v1, v2) : stabilisateurs dans <sigma>
for n, cs, x, attendu, orbs in [
    (6, [[1, 2, 3], [4, 5]], 4, 3, [[1, 2, 3], [4, 5], [6]]),
    (8, [[1, 2, 3, 4], [5, 6, 7]], 5, 4, [[1, 2, 3, 4], [5, 6, 7], [8]]),
    (8, [[1, 2, 3, 4, 5], [6, 7]], 1, 2, [[1, 2, 3, 4, 5], [6, 7], [8]]),
]:
    s = cyc(n, *cs)
    G = engendre([s])
    assert len(G) == ordre(s)
    assert sorted(sorted(o) for o in orbites_action(G, range(1, n + 1), perm_act)) == sorted(orbs)
    assert len(stab(G, x, perm_act)) == attendu
# descriptions explicites des stabilisateurs
s = cyc(6, [1, 2, 3], [4, 5])
assert set(stab(engendre([s]), 4, perm_act)) == {ident(6), cyc(6, [1, 3, 2]), cyc(6, [1, 2, 3])}
s = cyc(8, [1, 2, 3, 4], [5, 6, 7])
assert set(stab(engendre([s]), 5, perm_act)) == {ident(8), cyc(8, [1, 4, 3, 2]), cyc(8, [1, 3], [2, 4]), cyc(8, [1, 2, 3, 4])}
s = cyc(8, [1, 2, 3, 4, 5], [6, 7])
assert set(stab(engendre([s]), 1, perm_act)) == {ident(8), cyc(8, [6, 7])}

# --- gr-l3-actions-orbites-ex2 (et v1, v2) : S_n sur les paires, stabilisateur de {1,2} : 12, 48, 240
for n, attendu, leurres in [(5, 12, [6, 10, 60]), (6, 48, [24, 15, 360]), (7, 240, [120, 21, 2520])]:
    assert factorial(n) // comb(n, 2) == attendu == 2 * factorial(n - 2)
    assert leurres == [factorial(n - 2), comb(n, 2), factorial(n) // 2] and attendu not in leurres
for n in (5, 6):  # force brute
    Sn = S(n)
    orb = {paire_act(g, frozenset({1, 2})) for g in Sn}
    assert len(orb) == comb(n, 2)
    assert len(stab(Sn, frozenset({1, 2}), paire_act)) == 2 * factorial(n - 2)


# --- gr-l3-actions-orbites-ex3 (et v1, v2) : |GL_2(Z/pZ)| = 48, 480, 6 par l'action sur F_p^2 \ {0}
def gl2(p):
    return [(a, b, c, d) for a, b, c, d in product(range(p), repeat=4) if (a * d - b * c) % p]


for p, attendu, leurres in [(3, 48, [80, 64, 72]), (5, 480, [624, 576, 600]), (2, 6, [15, 9, 12])]:
    G = gl2(p)
    X = [(x, y) for x, y in product(range(p), repeat=2) if (x, y) != (0, 0)]
    act = lambda M, v, p=p: ((M[0] * v[0] + M[1] * v[1]) % p, (M[2] * v[0] + M[3] * v[1]) % p)  # noqa: E731
    assert {act(M, (1, 0)) for M in G} == set(X)  # transitive
    st = stab(G, (1, 0), act)
    assert len(st) == p * (p - 1) and all(M[0] == 1 and M[2] == 0 for M in st)
    assert len(G) == attendu == (p * p - 1) * (p * p - p) == len(X) * len(st)
    assert leurres == [p ** 4 - 1, (p * p - 1) ** 2, p * p * (p * p - 1)] and attendu not in leurres
# v2 : GL_2(F_2) opère fidèlement sur les 3 vecteurs non nuls
G = gl2(2)
X = [(1, 0), (0, 1), (1, 1)]
img = {tuple(((M[0] * v[0] + M[1] * v[1]) % 2, (M[2] * v[0] + M[3] * v[1]) % 2) for v in X) for M in G}
assert len(img) == 6


# --- gr-l3-actions-orbites-ex4 : lemme de Jordan (élément sans point fixe), et la majoration de la réunion des conjugués
def verifier_jordan(G, X, act):
    x = X[0]
    H = stab(G, x, act)
    conjugues = {frozenset(conj(g, h) for h in H) for g in G}
    indice = len(G) // len(H)
    assert len(conjugues) <= indice
    reunion = set().union(*conjugues)
    assert len(reunion) <= len(G) - indice + 1 < len(G)
    sans = [g for g in G if all(act(g, y) != y for y in X)]
    assert sans and set(sans) == set(G) - reunion


verifier_jordan(S4, [1, 2, 3, 4], perm_act)
verifier_jordan(A4, [1, 2, 3, 4], perm_act)
verifier_jordan(list(D), [1, 2, 3, 4], perm_act)
verifier_jordan(S4, P4, paire_act)


# --- gr-l3-actions-orbites-ex5 : |G| = 2m, m impair : H = noyau de eps∘rho a m éléments (S3, Z/6, D_5, Z/10, Z/2 x Z/3)
def regulier(elements, mult):
    idx = {g: i + 1 for i, g in enumerate(elements)}
    return {g: tuple(idx[mult(g, x)] for x in elements) for g in elements}


def groupe_diedral(n):
    """D_n comme permutations de {1..n} (sommet k+1 <-> omega^k) : r : k -> k+1, s : k -> -k."""
    r = tuple(((k + 1) % n) + 1 for k in range(n))
    s = tuple(((-k) % n) + 1 for k in range(n))
    return engendre([r, s]), r, s


for elements, mult, m in [
    (S3, comp, 3),
    (list(range(6)), lambda a, b: (a + b) % 6, 3),
    (list(groupe_diedral(5)[0]), comp, 5),
    (list(range(10)), lambda a, b: (a + b) % 10, 5),
]:
    rho = regulier(elements, mult)
    assert len(set(rho.values())) == 2 * m  # Cayley : injectif
    e = ident(len(elements[0])) if mult is comp else 0
    H = [g for g in elements if eps(rho[g]) == 1]
    assert len(H) == m
    involutions = [g for g in elements if g != e and mult(g, g) == e]
    assert involutions
    for t in involutions:
        assert type_(rho[t]) == (2,) * m and eps(rho[t]) == -1
assert set(H) == set(range(0, 10, 2))  # Z/10 : H = {pairs}
# remarque : A4 (12 = 2 x 6) n'a pas de sous-groupe à 6 éléments (un tel sous-groupe contiendrait les carrés des 3-cycles)
A4set = set(A4)
from itertools import combinations  # noqa: E402

trouve = False
for gens in combinations(A4, 2):
    if len(engendre(list(gens))) == 6:
        trouve = True
assert not trouve and all(len(engendre([g])) != 6 for g in A4)

# ===========================================================================
# Carte gr-l3-actions-formule-classes
# ===========================================================================
# exemple : rotations du carré sur les 16 coloriages des sommets : 2 points fixes, une orbite de 2, trois de 4
r4 = cyc(4, [1, 2, 3, 4])
R4 = list(engendre([r4]))
col_act = lambda g, c: tuple(c[inverse(g)[i] - 1] for i in range(len(c)))  # noqa: E731  (g.c)(A) = c(g^-1 A)
X16 = list(product((0, 1), repeat=4))
assert sorted(len(o) for o in orbites_action(R4, X16, col_act)) == [1, 1, 2, 4, 4, 4]
# contre-exemples et pièges : Z/6 sur Z/2 ⊔ Z/3 (5 points), sans point fixe
X5 = [("a", i) for i in range(2)] + [("b", i) for i in range(3)]
act6 = lambda k, x: (x[0], (x[1] + k) % (2 if x[0] == "a" else 3))  # noqa: E731
assert all(act6(0, x) == x for x in X5) and all(act6((k + l) % 6, x) == act6(k, act6(l, x)) for k in range(6) for l in range(6) for x in X5)
assert sorted(len(o) for o in orbites_action(range(6), X5, act6)) == [2, 3]
assert [x for x in X5 if all(act6(k, x) == x for k in range(6))] == []
assert 5 % 2 and 5 % 3 and 5 % 6  # 5 non congru à 0 modulo 2, 3 ou 6 (pièges « p divise |G| », « m quelconque », « pq »)
X2 = [("a", 0), ("a", 1)]
assert len(X2) % 6 != 0  # piège « m^n » : Z/6 (= 6^1) sur Z/2 sans point fixe, 2 non congru à 0 modulo 6
# leurres des trous : S3 sur {1,2,3} (somme des |G_x| = 2, r|G| = 6, somme des |Fix| = 6 != 3) ; Z/2 sur 2 points
assert sum(len(stab(S3, 1, perm_act)) for _ in [0]) == 2 and nb_orbites_burnside(S3, range(1, 4), perm_act) * 6 == 6
swap = lambda k, x: x if k == 0 else 1 - x  # noqa: E731
assert len(orbites_action(range(2), [0, 1], swap)) == 1 and 2 % 2 != 1  # nombre d'orbites 1, |X| = 2
assert [x for x in [0, 1] if all(swap(k, x) == x for k in range(2))] == []  # « p divise |X| », « |X| premier », « |X| >= p » : pas de point fixe


# --- gr-l3-actions-formule-classes-ex1 (et v1, v2) : nombre minimal de points fixes pour Z/N sur |X| points
def decompositions(total, tailles):
    """Toutes les listes (multiensembles) de tailles d'orbites de somme total."""
    if total == 0:
        yield []
        return
    for i, t in enumerate(tailles):
        if t <= total:
            for reste in decompositions(total - t, tailles[i:]):
                yield [t] + reste


for N, X, attendu, props in [(25, 12, 2, [2, 0, 1, 3]), (27, 22, 1, [1, 0, 2, 3]), (49, 17, 3, [3, 0, 1, 2])]:
    tailles = sorted(d for d in range(1, N + 1) if N % d == 0)
    mini = min(sum(1 for t in dec if t == 1) for dec in decompositions(X, tailles))
    assert mini == attendu == props[0] and len(set(props)) == 4
    # la construction proposée : Z/N agit sur Z/d par translation (bien défini), orbites de taille d
    for d in tailles:
        assert all(((k % N) - k) % d == 0 for k in range(-2 * N, 2 * N))

# --- gr-l3-actions-formule-classes-ex2 (et v1, v2) : mots circulaires, nombre d'orbites 51, 20, 24
for p, a, attendu, leurres in [(5, 3, 51, [48, 81, 243]), (7, 2, 20, [18, 64, 128]), (3, 4, 24, [20, 16, 64])]:
    X = list(product(range(a), repeat=p))
    dec = lambda k, w, p=p: tuple(w[(i + k) % p] for i in range(p))  # noqa: E731
    assert all(dec(k, dec(l, w)) == dec(k + l, w) for k in range(p) for l in range(p) for w in X[:50])
    orbs = orbites_action(range(p), X, dec)
    assert len(orbs) == attendu == a + (a ** p - a) // p
    assert sorted({len(o) for o in orbs}) == [1, p] and sum(1 for o in orbs if len(o) == 1) == a
    assert leurres == [(a ** p - a) // p, a ** (p - 1), a ** p] and attendu not in leurres
    assert (a ** p - a) % p == 0  # Fermat

# --- gr-l3-actions-formule-classes-ex3 (et v1, v2) : plus grand n non somme de diviseurs > 1 de N (Frobenius)
for N, u, v, attendu, props in [(15, 3, 5, 7, [7, 4, 8, 14]), (21, 3, 7, 11, [11, 8, 10, 20]), (35, 5, 7, 23, [23, 13, 24, 34])]:
    tailles = [d for d in range(2, N + 1) if N % d == 0]
    rep = lambda n: any(True for _ in decompositions(n, tailles))  # noqa: E731
    non = [n for n in range(1, 3 * N) if not rep(n)]
    assert max(non) == attendu == u * v - u - v == props[0]
    assert all(rep(n) == any((n - v * j) >= 0 and (n - v * j) % u == 0 for j in range(n // v + 1)) for n in range(1, 3 * N))
    assert all(x != attendu for x in props[1:])
    assert len(set(props)) == 4
assert [n for n in range(1, 45) if not any((n - 5 * j) % 3 == 0 for j in range(n // 5 + 1))] == [1, 2, 4, 7]
assert [n for n in range(1, 63) if not any((n - 7 * j) % 3 == 0 for j in range(n // 7 + 1))] == [1, 2, 4, 5, 8, 11]

# --- gr-l3-actions-formule-classes-ex4 : Burnside ; moyenne de f = 1, moyenne de f^2 = 2 pour n = 2..6
for n in range(2, 7):
    Sn = S(n)
    f = [sum(1 for i in range(n) if s[i] == i + 1) for s in Sn]
    assert F(sum(f), len(Sn)) == 1 and F(sum(x * x for x in f), len(Sn)) == 2
    if n <= 4:
        X = list(product(range(1, n + 1), repeat=2))
        act = lambda g, c: (g[c[0] - 1], g[c[1] - 1])  # noqa: E731
        assert len(orbites_action(Sn, X, act)) == 2 == nb_orbites_burnside(Sn, X, act)
# formule de Burnside contrôlée sur plusieurs actions
for G, X, act in [(S4, P4, paire_act), (R4, X16, col_act), (list(D), X16, col_act)]:
    assert len(orbites_action(G, X, act)) == nb_orbites_burnside(G, X, act)
assert F(3, 2) != 2 and 3 != 2  # leurres


# --- gr-l3-actions-formule-classes-ex5 : McKay / Cauchy ; |X| = |G|^(p-1), points fixes = solutions de x^p = e
def verifier_mckay(G, p):
    G = list(G)
    e = ident(len(G[0]))
    X = [t for t in product(G, repeat=p) if comp(*t) == e]
    assert len(X) == len(G) ** (p - 1)
    assert all(comp(*(t[1:] + t[:1])) == e for t in X)
    fixes = [t for t in X if len(set(t)) == 1]
    sol = [x for x in G if comp(*([x] * p)) == e]
    assert len(fixes) == len(sol) and len(sol) % p == 0 and (len(sol) - 1) % p == p - 1
    assert any(ordre(x) == p for x in G)


verifier_mckay(S3, 2)
verifier_mckay(S3, 3)
verifier_mckay(D, 2)
verifier_mckay(A4, 3)
verifier_mckay(S4, 3)
assert sum(1 for x in S4 if ordre(x) == 3) == 8 and 8 % 3 == 2

# ===========================================================================
# Carte gr-l3-actions-conjugaison
# ===========================================================================
# exemple : classes de S4 de cardinaux 1, 6, 3, 8, 6 ; centralisateur d'un 3-cycle = <c>
cl4 = classes(S4)
assert sorted(len(c) for c in cl4) == [1, 3, 6, 6, 8] and centre(S4) == [ident(4)]
c = cyc(4, [1, 2, 3])
assert set(centralisateur(S4, c)) == engendre([c])
# classes de S_n = types (n <= 6) ; formule de conjugaison d'un cycle
for n in range(2, 7):
    Sn = S(n)
    cl = classes(Sn)
    assert all(len({type_(x) for x in k}) == 1 for k in cl)
    assert len(cl) == len({type_(x) for x in Sn})
    for k in cl:
        x = next(iter(k))
        assert len(k) * len(centralisateur(Sn, x)) == len(Sn)
for t in S5:
    cy = [1, 3, 4, 5]
    assert conj(t, cyc(5, cy)) == cyc(5, [t[a - 1] for a in cy])
# équation aux classes et centre non trivial pour les groupes de cardinal p^n disponibles (D d'ordre 8, Z/9 via <(1..9)>)
for G in [list(D), list(engendre([cyc(9, [1, 2, 3, 4, 5, 6, 7, 8, 9])])), list(engendre([cyc(6, [1, 2, 3]), cyc(6, [4, 5, 6])]))]:
    assert len(centre(G)) > 1
assert centre(S3) == [ident(3)]  # pièges : 6 divisible par 2 et 3 (= 2·3, = 6^1), centre trivial
# contre-exemples : (1 2) et (1 2)(3 4) même ordre, non conjugués ; (1 2 3) et (1 3 2) conjugués dans S3 pas dans A3
assert ordre(cyc(4, [1, 2])) == ordre(cyc(4, [1, 2], [3, 4])) == 2 and type_(cyc(4, [1, 2])) != type_(cyc(4, [1, 2], [3, 4]))
assert conj(cyc(3, [2, 3]), cyc(3, [1, 2, 3])) == cyc(3, [1, 3, 2])
A3 = list(engendre([cyc(3, [1, 2, 3])]))
assert all(conj(g, cyc(3, [1, 2, 3])) == cyc(3, [1, 2, 3]) for g in A3)
# affirmation : centre de D = {id, (1 3)(2 4)}, D non commutatif
assert set(centre(D)) == {ident(4), cyc(4, [1, 3], [2, 4])}
# trous : leurres faux (signature, ordre, points fixes), cardinal de classe
assert eps(cyc(4, [1, 2, 3])) == eps(cyc(4, [1, 2], [3, 4])) and type_(cyc(4, [1, 2, 3])) != type_(cyc(4, [1, 2], [3, 4]))
fixes = lambda p: sum(1 for i, x in enumerate(p, 1) if x == i)  # noqa: E731
assert fixes(cyc(4, [1, 2, 3, 4])) == fixes(cyc(4, [1, 2], [3, 4])) == 0
k12 = next(k for k in classes(S3) if cyc(3, [1, 2]) in k)
assert len(k12) == 3 != ordre(cyc(3, [1, 2])) and len(k12) != len(centralisateur(S3, cyc(3, [1, 2])))
assert len(next(k for k in classes(S3) if cyc(3, [1, 2, 3]) in k)) == 2 != len(classes(S3))
# « de cardinal p » faux : Z/9 (via <(1..9)>) commutatif de cardinal 3^2, centre à 9 éléments
assert len(centre(list(engendre([cyc(9, list(range(1, 10)))])))) == 9

# --- gr-l3-actions-conjugaison-ex1 (et v1, v2)
for n, s, t, rep, leurres in [
    (5, [[1, 2, 3, 4]], [[1, 3, 5], [2, 4]], [(2, 3, 4, 5)], [[(1, 2, 5, 4)], [(2, 5, 4, 3)], [(1, 2, 3, 4)]]),
    (5, [[1, 2, 3], [4, 5]], [[1, 4], [2, 5, 3]], [(1, 3), (2, 4, 5)], [[(1, 2), (3, 5, 4)], [(1, 3), (2, 5, 4)], [(1, 2, 3), (4, 5)]]),
    (5, [[1, 5, 2, 4]], [[1, 2, 3, 4, 5]], [(1, 3, 5, 2)], [[(1, 3, 5, 4)], [(1, 2, 5, 3)], [(1, 5, 2, 4)]]),
]:
    sg, tg = cyc(n, *s), cyc(n, *t)
    assert norm(cycles(conj(tg, sg))) == norm(rep)
    assert norm(cycles(conj(inverse(tg), sg))) == norm(leurres[0])
    for l in leurres:
        assert norm(l) != norm(rep)

# --- gr-l3-actions-conjugaison-ex2 (et v1, v2) : centralisateurs 8, 6, 48
for G, x, ncl, attendu in [
    (S5, cyc(5, [1, 2], [3, 4]), 15, 8),
    (S5, cyc(5, [1, 2, 3]), 20, 6),
    (S6, cyc(6, [1, 2], [3, 4], [5, 6]), 15, 48),
]:
    assert sum(1 for y in G if type_(y) == type_(x)) == ncl
    assert len(centralisateur(G, x)) == attendu == len(G) // ncl
assert len(engendre([cyc(5, [1, 2, 3]), cyc(5, [4, 5])])) == 6
assert set(centralisateur(S5, cyc(5, [1, 2, 3]))) == engendre([cyc(5, [1, 2, 3]), cyc(5, [4, 5])])

# --- gr-l3-actions-conjugaison-ex3 : A4 a 4 classes (1, 3, 4, 4) ; v1 : D_6 a 6 classes ; v2 : A5 a 5 classes
clA4 = classes(A4)
assert sorted(len(k) for k in clA4) == [1, 3, 4, 4]
k123 = next(k for k in clA4 if cyc(4, [1, 2, 3]) in k)
assert k123 == {cyc(4, [1, 2, 3]), cyc(4, [1, 3, 4]), cyc(4, [1, 4, 2]), cyc(4, [2, 4, 3])}
assert cyc(4, [1, 3, 2]) not in k123  # QCM q4
assert {inverse(x) for x in k123} == next(k for k in clA4 if cyc(4, [1, 3, 2]) in k)
assert len({type_(x) for x in A4}) == 3  # leurre « 3 »
D6, r6, s6 = groupe_diedral(6)
D6 = list(D6)
assert len(D6) == 12
cl6 = classes(D6)
assert len(cl6) == 6 and sorted(len(k) for k in cl6) == [1, 1, 2, 2, 3, 3]
assert set(centre(D6)) == {ident(6), comp(r6, r6, r6)}
assert comp(s6, r6, s6) == inverse(r6)
clA5 = classes(A5)
assert sorted(len(k) for k in clA5) == [1, 12, 12, 15, 20]
assert len({type_(x) for x in A5}) == 4 and len(classes(S5)) == 7

# --- gr-l3-actions-conjugaison-ex4 : groupes de cardinal p^2 commutatifs (exemples), p^2 classes ; leurres distincts
for G in [list(engendre([cyc(4, [1, 2, 3, 4])])), list(engendre([cyc(4, [1, 2]), cyc(4, [3, 4])])),
          list(engendre([cyc(9, list(range(1, 10)))])), list(engendre([cyc(6, [1, 2, 3]), cyc(6, [4, 5, 6])]))]:
    assert len(classes(G)) == len(G) and len(centre(G)) == len(G)
for p in (2, 3, 5, 7):
    vals = [p * p, p, p * p - 1, p * p + p - 1]
    assert len(set(vals)) == 4

# --- gr-l3-actions-conjugaison-ex5 : P(G) = k(G)/|G| <= 5/8 si G non commutatif, égalité pour D
for G in [S3, list(D), A4, S4, D6, list(groupe_diedral(5)[0]), A5]:
    paires = sum(1 for g in G for h in G if comp(g, h) == comp(h, g))
    assert paires == sum(len(centralisateur(G, g)) for g in G) == len(classes(G)) * len(G)
    assert F(paires, len(G) ** 2) <= F(5, 8)
    assert len(G) // len(centre(G)) >= 4 and len(G) // len(centre(G)) not in (2, 3, 5, 7, 11, 13)
assert F(len(classes(D)), 8) == F(5, 8) and F(len(classes(S3)), 6) == F(1, 2)
sg, tg = cyc(4, [1, 2, 3, 4]), cyc(4, [1, 3])
assert conj(tg, sg) == inverse(sg)
puiss = [ident(4), sg, comp(sg, sg), comp(sg, sg, sg)]
assert set(D) == set(puiss) | {comp(p, tg) for p in puiss}
attendu = [{puiss[0]}, {puiss[2]}, {puiss[1], puiss[3]}, {tg, comp(puiss[2], tg)}, {comp(puiss[1], tg), comp(puiss[3], tg)}]
assert sorted(map(sorted, classes(D))) == sorted(map(sorted, attendu))

# ===========================================================================
# Carte gr-l3-actions-automorphismes
# ===========================================================================


def automorphismes(G, gens):
    """Automorphismes d'un groupe de permutations, déterminés par les images des générateurs (force brute)."""
    G = list(G)
    e = ident(len(G[0]))
    # écriture de chaque élément comme produit des générateurs
    mot = {e: []}
    file = deque([e])
    while file:
        x = file.popleft()
        for i, g in enumerate(gens):
            y = comp(x, g)
            if y not in mot:
                mot[y] = mot[x] + [i]
                file.append(y)
    res = []
    for imgs in product(G, repeat=len(gens)):
        if any(ordre(a) != ordre(b) for a, b in zip(gens, imgs)):
            continue
        f = {x: comp(*([e] + [imgs[i] for i in w])) for x, w in mot.items()}
        if len(set(f.values())) == len(G) and all(f[comp(x, y)] == comp(f[x], f[y]) for x in G for y in G):
            res.append(f)
    return res


i_ = lambda g: (lambda x: conj(g, x))  # noqa: E731

# énoncé et pièges sur S3
for g in S3:
    for h in S3:
        assert all(conj(g, conj(h, x)) == conj(comp(g, h), x) for x in S3)  # i_g ∘ i_h = i_gh
ig_inv = lambda g, x: comp(inverse(g), x, g)  # noqa: E731
assert any(any(ig_inv(comp(g, h), x) != ig_inv(g, ig_inv(h, x)) for x in S3) for g in S3 for h in S3)  # piège 1
assert any(comp(g, ident(3), g) != ident(3) for g in S3)  # piège 2 : x -> gxg n'envoie pas e sur e
p3 = lambda g, x: comp(g, inverse(x), inverse(g))  # noqa: E731
assert any(p3(ident(3), comp(x, y)) != comp(p3(ident(3), x), p3(ident(3), y)) for x in S3 for y in S3)  # piège 3
# noyau = centre ; |Int| = |G|/|Z|
for G in [S3, list(D), A4, S4]:
    Z = centre(G)
    imgs = {tuple(conj(g, x) for x in G) for g in G}
    assert len(imgs) == len(G) // len(Z)
    assert [g for g in G if all(conj(g, x) == x for x in G)] == Z
# Aut(S3) : 6 automorphismes, tous intérieurs ; phi ∘ i_g ∘ phi^-1 = i_phi(g) ; leurres i_g et i_{phi^-1(g)} faux
auts3 = automorphismes(S3, [cyc(3, [1, 2]), cyc(3, [1, 2, 3])])
assert len(auts3) == 6
for f in auts3:
    finv = {v: k for k, v in f.items()}
    for g in S3:
        assert all(f[conj(g, finv[x])] == conj(f[g], x) for x in S3)
f = {x: conj(cyc(3, [1, 2, 3]), x) for x in S3}
finv = {v: k for k, v in f.items()}
g = cyc(3, [1, 2])
assert f[g] != g and f[g] != finv[g]
assert any(conj(f[g], x) != conj(g, x) for x in S3) and any(conj(f[g], x) != conj(finv[g], x) for x in S3)
assert any(f[x] != conj(f[g], x) for x in S3)  # leurre « phi »
assert len([x for x in S3 if ordre(x) <= 2]) == 4 != len(centre(S3))  # leurre « éléments d'ordre au plus 2 »
# exemple V = (Z/2)^2 : 6 automorphismes ; contre-exemple x -> -x sur Z/3 ; x -> x^-1 non morphisme dans S3
V = [(a, b) for a in range(2) for b in range(2)]
add = lambda u, v: ((u[0] + v[0]) % 2, (u[1] + v[1]) % 2)  # noqa: E731
autV = [dict(zip(V, im)) for im in permutations(V) if dict(zip(V, im))[(0, 0)] == (0, 0)
        and all(dict(zip(V, im))[add(u, v)] == add(dict(zip(V, im))[u], dict(zip(V, im))[v]) for u in V for v in V)]
assert len(autV) == 6
assert all(((-(a + b)) % 3) == ((-a) % 3 + (-b) % 3) % 3 for a in range(3) for b in range(3))
a, b = cyc(3, [1, 2]), cyc(3, [1, 3])
assert inverse(comp(a, b)) == cyc(3, [1, 2, 3]) and comp(inverse(a), inverse(b)) == cyc(3, [1, 3, 2])
assert all((x + x) % 2 == 0 for x in range(2))  # x -> x + x nulle sur Z/2

# --- gr-l3-actions-automorphismes-ex1 (et v1, v2) : PAP^-1 (matrices entières, calcul exact)
for P, Am, rep, leurres in [
    ([[1, 1], [0, 1]], [[1, 0], [0, 2]], [[1, 1], [0, 2]], [[[1, -1], [0, 2]], [[1, 2], [0, 2]], [[2, 1], [0, 1]]]),
    ([[2, 1], [1, 1]], [[0, 1], [1, 0]], [[-1, 3], [0, 1]], [[[-1, 0], [3, 1]], [[1, 3], [0, -1]], [[-1, -3], [0, 1]]]),
    ([[1, 0], [2, 1]], [[3, 1], [0, 3]], [[1, 1], [-4, 5]], [[[5, 1], [-4, 1]], [[1, -1], [4, 5]], [[3, 1], [0, 3]]]),
]:
    Pn, An = np.array(P), np.array(Am)
    assert round(np.linalg.det(Pn)) == 1
    Pinv = np.round(np.linalg.inv(Pn)).astype(int)
    assert (Pn @ Pinv == np.eye(2, dtype=int)).all()
    assert (Pn @ An @ Pinv == np.array(rep)).all()
    assert (Pinv @ An @ Pn == np.array(leurres[0])).all()
    for L in leurres:  # même trace et même déterminant (l'indice de contrôle n'élimine aucun leurre), mais différents
        L = np.array(L)
        assert not (L == np.array(rep)).all()
        assert np.trace(L) == np.trace(An) and round(np.linalg.det(L)) == round(np.linalg.det(An))

# --- gr-l3-actions-automorphismes-ex2 (et v1, v2) : |Int| = 6, 4, 12
for G, attendu in [(S3, 6), (list(D), 4), (A4, 12)]:
    assert len({tuple(conj(g, x) for x in G) for g in G}) == attendu
assert set(centralisateur(S3, cyc(3, [1, 2]))) == {ident(3), cyc(3, [1, 2])}
assert set(centralisateur(list(D), sg)) == set(puiss)
assert conj(cyc(4, [1, 2, 3]), cyc(4, [1, 2], [3, 4])) == cyc(4, [2, 3], [1, 4])
assert conj(cyc(4, [1, 3, 2]), cyc(4, [1, 2], [3, 4])) == cyc(4, [1, 3], [2, 4])

# --- gr-l3-actions-automorphismes-ex3 (et v1, v2) : |Aut(Z/n)| = phi(n) = 4, 8, 6
for n, attendu, gens in [(12, 4, [1, 5, 7, 11]), (15, 8, [1, 2, 4, 7, 8, 11, 13, 14]), (9, 6, [1, 2, 4, 5, 7, 8])]:
    auts = [a for a in range(n) if len({(k * a) % n for k in range(n)}) == n]
    assert auts == gens and len(auts) == attendu == sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


# --- gr-l3-actions-automorphismes-ex4 : groupe non commutatif de cardinal p^3 : centre à p éléments (p = 2, 3)
def heisenberg(p):
    return [(a, b, c) for a in range(p) for b in range(p) for c in range(p)]


def hmul(x, y, p):
    # [[1,a,b],[0,1,c],[0,0,1]] * [[1,a',b'],[0,1,c'],[0,0,1]] = [[1,a+a',b+b'+ac'],[0,1,c+c'],[0,0,1]]
    return ((x[0] + y[0]) % p, (x[1] + y[1] + x[0] * y[2]) % p, (x[2] + y[2]) % p)


for p in (2, 3):
    H = heisenberg(p)
    assert any(hmul(x, y, p) != hmul(y, x, p) for x in H for y in H)
    Z = [z for z in H if all(hmul(z, g, p) == hmul(g, z, p) for g in H)]
    assert len(Z) == p and all(z[0] == 0 and z[2] == 0 for z in Z)
assert len(centre(D)) == 2

# --- gr-l3-actions-automorphismes-ex5 : Aut(S4) = Int(S4), 24 éléments
auts4 = automorphismes(S4, [cyc(4, [1, 2]), cyc(4, [1, 2, 3, 4])])
assert len(auts4) == 24 and centre(S4) == [ident(4)]
for f in auts4:
    assert any(all(f[x] == conj(t, x) for x in S4) for t in S4)
    tr = [x for x in S4 if type_(x) == (2,)]
    assert {f[x] for x in tr} == set(tr)
assert len(engendre([cyc(4, [1, 2]), cyc(4, [1, 3]), cyc(4, [2, 3])])) == 6  # cas « triangle » : sous-groupe à 6 éléments
assert len(engendre([cyc(4, [1, 2]), cyc(4, [1, 3]), cyc(4, [1, 4])])) == 24
assert comp(cyc(4, [1, 2]), cyc(4, [1, 3])) == cyc(4, [1, 3, 2])
assert sum(1 for x in S4 if type_(x) == (2,)) == 6 and sum(1 for x in S4 if type_(x) == (2, 2)) == 3

# ===========================================================================
# Carte gr-l3-actions-geometrie
# ===========================================================================
# énoncé : r_k, s_k (complexes) conservent S, sont 2n applications distinctes, forment un groupe ; action transitive,
# stabilisateur de 1 = {id, conj}, fidèle (n = 3..10)
for n in range(3, 11):
    w = np.exp(2j * np.pi / n)
    Sv = [w ** k for k in range(n)]
    idx = lambda z, Sv=Sv: min(range(len(Sv)), key=lambda j: abs(Sv[j] - z))  # noqa: E731
    maps = [(lambda z, k=k: w ** k * z) for k in range(n)] + [(lambda z, k=k: w ** k * np.conj(z)) for k in range(n)]
    perms = []
    for f in maps:
        imgs = [idx(f(z)) for z in Sv]
        assert all(abs(f(z) - Sv[i]) < 1e-9 for z, i in zip(Sv, imgs))
        perms.append(tuple(i + 1 for i in imgs))
    assert len(set(perms)) == 2 * n  # fidèle et 2n éléments distincts
    Gp = set(perms)
    assert all(comp(a, b) in Gp for a in Gp for b in Gp)
    assert {p[0] for p in Gp} == set(range(1, n + 1))  # transitive
    st = [p for p in Gp if p[0] == 1]
    assert len(st) == 2 and tuple(((-k) % n) + 1 for k in range(n)) in st
    assert Gp == set(groupe_diedral(n)[0])
    assert 2 * n != n and 2 * n != n * n and 2 * n != 4 * n  # leurres du trou 3
    # isométries « directes » seules : n éléments, stabilisateur trivial
    assert len(set(perms[:n])) == n and sum(1 for p in perms[:n] if p[0] == 1) == 1
    # z -> -z ne fixe pas le sommet 1
    assert abs(-Sv[0] - Sv[0]) > 1
assert len(groupe_diedral(3)[0]) == 6 == factorial(3)
# exemple : D_4 engendré par (1 2 3 4) et (2 4) ; la transposition (1 2) n'en fait pas partie
D4 = engendre([cyc(4, [1, 2, 3, 4]), cyc(4, [2, 4])])
assert D4 == set(groupe_diedral(4)[0]) and len(D4) == 8 and cyc(4, [1, 2]) not in D4
# contre-exemple S = {1, i, -1} : les permutations de S qui conservent les distances sont id et l'échange de 1 et -1
pts = [1, 1j, -1]
iso = [p for p in permutations(range(3)) if all(abs(abs(pts[p[a]] - pts[p[b]]) - abs(pts[a] - pts[b])) < 1e-12 for a in range(3) for b in range(3))]
assert len(iso) == 2 and abs(np.mean(pts) - 1j / 3) < 1e-12
# démonstration : |conj(u) - w|^2 - |u - w|^2 = 4 Im u Im w
for u, w0 in [(0.3 + 0.7j, np.exp(2j * np.pi / 5)), (-1.2 + 0.4j, np.exp(2j * np.pi / 3))]:
    assert abs(abs(np.conj(u) - w0) ** 2 - abs(u - w0) ** 2 - 4 * u.imag * w0.imag) < 1e-12

# --- gr-l3-actions-geometrie-ex1 (et v1, v2) : stabilisateurs dans les groupes de rotations (cube, tétraèdre, dodécaèdre)


def perm_signees():
    out = []
    for p in permutations(range(3)):
        for s in product((1, -1), repeat=3):
            M = np.zeros((3, 3), dtype=int)
            for j in range(3):
                M[p[j], j] = s[j]
            out.append(M)
    return out


W = perm_signees()
Gc = [M for M in W if round(np.linalg.det(M)) == 1]
assert len(W) == 48 and len(Gc) == 24
vec = lambda M, v: tuple(int(x) for x in M @ np.array(v))  # noqa: E731
sommets = [v for v in product((1, -1), repeat=3)]
aretes = [frozenset({a, b}) for a in sommets for b in sommets if a < b and sum(x != y for x, y in zip(a, b)) == 1]
assert len(aretes) == 12
act_ar = lambda M, e: frozenset(vec(M, v) for v in e)  # noqa: E731
assert {act_ar(M, aretes[0]) for M in Gc} == set(aretes)
assert len(stab(Gc, aretes[0], act_ar)) == 2
tetra = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
Gt = [M for M in Gc if {vec(M, v) for v in tetra} == set(tetra)]
assert len(Gt) == 12
faces_t = [frozenset(set(tetra) - {v}) for v in tetra]
act_f = lambda M, f: frozenset(vec(M, v) for v in f)  # noqa: E731
assert {act_f(M, faces_t[0]) for M in Gt} == set(faces_t) and len(stab(Gt, faces_t[0], act_f)) == 3
# dodécaèdre : ses 12 faces correspondent aux 12 sommets de l'icosaèdre dual (0, ±1, ±φ) et permutations circulaires
phi = (1 + 5 ** 0.5) / 2
ico = []
for a, b in product((1, -1), repeat=2):
    ico += [(0, a, b * phi), (a, b * phi, 0), (b * phi, 0, a)]
ico = np.array(ico)


def cle(v):
    return tuple(np.round(v, 6))


ens = {cle(v) for v in ico}
v0 = ico[0]
voisins = [v for v in ico if abs(np.linalg.norm(v - v0) - 2) < 1e-9]
v1 = voisins[0]
v2 = next(v for v in voisins if abs(np.linalg.norm(v - v1) - 2) < 1e-9)
B = np.array([v0, v1, v2]).T
Gd = []
for a in ico:
    for b in ico:
        for c in ico:
            M = np.array([a, b, c]).T @ np.linalg.inv(B)
            if np.allclose(M.T @ M, np.eye(3), atol=1e-9) and abs(np.linalg.det(M) - 1) < 1e-9 and {cle(M @ v) for v in ico} == ens:
                Gd.append(M)
assert len(Gd) == 60
assert {cle(M @ v0) for M in Gd} == ens  # transitive sur les 12 sommets de l'icosaèdre (faces du dodécaèdre)
assert sum(1 for M in Gd if np.allclose(M @ v0, v0)) == 5
for n, attendu, leurres in [(24, 2, [12, 4, 3]), (12, 3, [4, 12, 1]), (60, 5, [12, 60, 3])]:
    assert attendu not in leurres

# --- gr-l3-actions-geometrie-ex2 (et v1, v2) : orbites de D_n sur les paires de sommets : floor(n/2) = 3, 4, 5
for n, attendu, nb, tailles in [(6, 3, 15, [3, 6, 6]), (9, 4, 36, [9, 9, 9, 9]), (10, 5, 45, [5, 10, 10, 10, 10])]:
    Dn = list(groupe_diedral(n)[0])
    paires = [frozenset({a, b}) for a in range(1, n + 1) for b in range(a + 1, n + 1)]
    assert len(paires) == nb
    orbs = orbites_action(Dn, paires, paire_act)
    assert len(orbs) == attendu == n // 2 and sorted(len(o) for o in orbs) == tailles
    assert n - 1 not in (attendu,)  # leurre « distinguer d et n - d »


# --- gr-l3-actions-geometrie-ex3 (et v1, v2) : coloriages des sommets sous D_n : 13, 39, 30 ; rotations seules 14, 51, 36
def coloriages(G, nsom, k):
    X = list(product(range(k), repeat=nsom))
    return len(orbites_action(G, X, col_act))


for n, k, attendu, rot in [(6, 2, 13, 14), (5, 3, 39, 51), (8, 2, 30, 36)]:
    Dn, r, s = groupe_diedral(n)
    Rn = list(engendre([r]))
    assert coloriages(list(Dn), n, k) == attendu == nb_orbites_burnside(list(Dn), list(product(range(k), repeat=n)), col_act)
    assert coloriages(Rn, n, k) == rot
    # nombres de cycles annoncés dans la correction
    for g in Rn:
        kk = next(j for j in range(1, n) if comp(*([r] * j)) == g) if g != ident(n) else 0
        assert len(cycles(g)) + sum(1 for i in range(n) if g[i] == i + 1) == gcd(kk, n)
assert (64 + 2 + 4 + 8 + 4 + 2) == 84 and 3 * 16 + 3 * 8 == 72 and 156 // 12 == 13
assert 243 + 12 == 255 and 5 * 27 == 135 and 390 // 10 == 39
assert 256 + 8 + 8 + 16 == 288 and 4 * 32 + 4 * 16 == 192 and 480 // 16 == 30

# --- gr-l3-actions-geometrie-ex4 (et v1, v2) : faces du cube coloriées, à rotation près : 57, 10, 240 (isométries : 56, 10, 220)
faces = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
pf = lambda M: tuple(faces.index(vec(M, f)) + 1 for f in faces)  # noqa: E731
Gcf = [pf(M) for M in Gc]
Wf = [pf(M) for M in W]
types_rot = sorted(tuple(sorted(len(c) for c in cycles(g))) for g in Gcf)
assert types_rot.count(()) == 1 and types_rot.count((4,)) == 6 and types_rot.count((2, 2)) == 3
assert types_rot.count((3, 3)) == 8 and types_rot.count((2, 2, 2)) == 6
for k, attendu, iso_tot in [(3, 57, 56), (2, 10, 10), (4, 240, 220)]:
    assert coloriages(Gcf, 6, k) == attendu == (k ** 6 + 3 * k ** 4 + 12 * k ** 3 + 8 * k ** 2) // 24
    assert coloriages(Wf, 6, k) == iso_tot
assert (4 ** 6) % 24 != 0
assert [sum(1 for o in orbites_action(Gcf, list(product(range(2), repeat=6)), col_act) if sum(next(iter(o))) == j)
        for j in range(7)] == [1, 1, 2, 2, 2, 1, 1]

# --- gr-l3-actions-geometrie-ex5 : rotations du cube ≅ S4 par l'action sur les 4 grandes diagonales ; 9 éléments d'ordre 2
diag = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]
droite = lambda v: diag.index(v) if v in diag else diag.index(tuple(-x for x in v))  # noqa: E731
rho = {i: tuple(droite(vec(M, d)) + 1 for d in diag) for i, M in enumerate(Gc)}
assert all(sorted(p) == [1, 2, 3, 4] for p in rho.values())
assert len(set(rho.values())) == 24  # injectif, donc bijectif sur S4
for i, M in enumerate(Gc):
    for j, N in enumerate(Gc):
        k = next(t for t, Q in enumerate(Gc) if (Q == M @ N).all())
        assert rho[k] == comp(rho[i], rho[j])  # morphisme
assert round(np.linalg.det(np.array(diag[1:]))) == 4 and tuple(np.sum(np.array(diag[1:]), axis=0)) == diag[0]
ordre_mat = lambda M: next(k for k in range(1, 7) if (np.linalg.matrix_power(M, k) == np.eye(3, dtype=int)).all())  # noqa: E731
assert sum(1 for M in Gc if ordre_mat(M) == 2) == 9
M = np.diag([1, -1, -1])
assert rho[next(i for i, Q in enumerate(Gc) if (Q == M).all())] == cyc(4, [1, 2], [3, 4])

# ===========================================================================
# QCM gr-l3-actions-q1 … q8
# ===========================================================================
assert 12 // 4 == 3 and len(stab(A4, 4, perm_act)) == 3  # q1
assert all(len({g for g in G if all(g[x - 1] == x for x in range(1, len(g) + 1))}) == 1 for G in [S3, S4])  # q2 (fidèle)
assert [x for x in [4, 0, 2, 3] if x % 3 == 31 % 3] == [4] and 27 + 4 == 31  # q3
# q4 : vérifié plus haut (classe de (1 2 3) dans A4)
assert len(classes(S5)) == 7  # q5
assert (-1) % 3 != 1  # q6 : x -> -x n'est pas l'identité sur Z/3
assert len(centre(D)) == 2 and len({tuple(conj(g, x) for x in D) for g in D}) == 4  # q7
D4s = list(groupe_diedral(4)[0])
assert len(D4s) == 8 and len(stab(D4s, 1, perm_act)) == 2  # q8

print("gr_l3_actions : toutes les vérifications passent.")
