"""Vérifications calculatoires du chapitre « gr-l3-quotients » (groupes, L3) : cartes, exercices et QCM de
content/*/groupes/l3-gr-l3-quotients.json et content/qcm/gr-l3-quotients.json.

Convention du chapitre : sigma tau = sigma ∘ tau (on applique d'abord tau). Chaque assert cite l'identifiant de la carte,
de l'exercice (variante comprise : v1, v2) ou de la question. Calculs exacts : force brute dans de petits groupes
(sous-groupes, sous-groupes distingués, quotients, morphismes), fractions exactes pour Q et R/Z, racines de l'unité
représentées par leurs exposants k/n modulo 1. Lancer : python3 verification/gr_l3_quotients.py (doit finir sans erreur).
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm

import sympy as sp

# ---------------------------------------------------------------------------
# Permutations de {1..n} : tuple p avec p[i-1] = image de i.
# ---------------------------------------------------------------------------


def ident(n):
    return tuple(range(1, n + 1))


def cyc(n, *cs):
    """Produit de cycles à supports disjoints, chacun donné par une liste."""
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


def inv(p):
    q = [0] * len(p)
    for i, x in enumerate(p, 1):
        q[x - 1] = i
    return tuple(q)


def conj(s, x):
    return comp(s, x, inv(s))


def engendre(gens, mul, e):
    """Sous-groupe engendré dans un groupe fini : fermeture par produit (suffit, voir la carte engendre)."""
    H = {e} | set(gens)
    file = list(H)
    while file:
        x = file.pop()
        for g in list(gens):
            for y in (mul(x, g), mul(g, x)):
                if y not in H:
                    H.add(y)
                    file.append(y)
    return frozenset(H)


def tous_sous_groupes(G, mul, e):
    """Tous les sous-groupes d'un petit groupe fini : sous-groupes engendrés par au plus 2 éléments, puis fermeture
    par sous-groupe engendré par la réunion (suffisant pour les groupes d'ordre <= 24 testés ici)."""
    S = {engendre([a, b], mul, e) for a in G for b in G}
    S.add(frozenset([e]))
    change = True
    while change:
        change = False
        for A in list(S):
            for B in list(S):
                C = engendre(list(A | B), mul, e)
                if C not in S:
                    S.add(C)
                    change = True
    return S


def distingue(H, G, mul, inverse):
    return all(mul(mul(g, h), inverse(g)) in H for g in G for h in H)


def Sn(n):
    from itertools import permutations
    return [tuple(p) for p in permutations(range(1, n + 1))]


def signe(p):
    n = len(p)
    return (-1) ** sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])


def ordre(x, mul, e):
    k, y = 1, x
    while y != e:
        y = mul(y, x)
        k += 1
    return k


S3, S4 = Sn(3), Sn(4)
e3, e4 = ident(3), ident(4)
V4 = frozenset([e4, cyc(4, [1, 2], [3, 4]), cyc(4, [1, 3], [2, 4]), cyc(4, [1, 4], [2, 3])])
A4 = frozenset(p for p in S4 if signe(p) == 1)
A3 = frozenset(p for p in S3 if signe(p) == 1)
r, s = cyc(4, [1, 2, 3, 4]), cyc(4, [1, 3])
D = engendre([r, s], comp, e4)

# ===========================================================================
# Carte gr-l3-quotients-engendre
# ===========================================================================
# exemple : éléments de <r, s> et relations.
assert comp(r, r) == cyc(4, [1, 3], [2, 4]) and comp(r, r, r) == cyc(4, [1, 4, 3, 2])
assert comp(r, s) == cyc(4, [1, 4], [2, 3]) and comp(r, r, s) == cyc(4, [2, 4])
assert comp(r, r, r, s) == cyc(4, [1, 2], [3, 4]) and comp(s, r) == cyc(4, [1, 2], [3, 4])
assert comp(s, r) == comp(inv(r), s)
assert len(D) == 8 and D == {comp(*([r] * k + [s] * l)) if k + l else e4 for k in range(4) for l in range(2)}
# contre-exemple / affirmation 4 : a = (1 2), b = (1 3) dans S3.
a, b = cyc(3, [1, 2]), cyc(3, [1, 3])
assert comp(a, b) == cyc(3, [1, 3, 2]) and comp(b, a) == cyc(3, [1, 2, 3])
akbl = {comp(*([a] * k + [b] * l)) if k + l else e3 for k in range(2) for l in range(2)}
assert len(akbl) == 4 and comp(b, a) not in akbl and len(engendre([a, b], comp, e3)) == 6
# trou 3 (leurres « d'ordre 2 », « d'ordre fini », « engendrent un groupe fini ») : même contre-exemple.
assert ordre(a, comp, e3) == ordre(b, comp, e3) == 2
# contre-exemple : 2Z ∪ 3Z non stable (5), <2,3> = Z (gcd).
assert 5 % 2 and 5 % 3 and gcd(2, 3) == 1
# affirmation 2 / ex4 : s(x) = -x, t(x) = 1 - x ; t∘s = translation de 1 (vérifié symboliquement).
x, k_, l_ = sp.symbols("x k l")
s_f = lambda u: -u
t_f = lambda u: 1 - u
assert sp.simplify(t_f(s_f(x)) - (x + 1)) == 0 and s_f(s_f(x)) == x and sp.expand(t_f(t_f(x))) == x
tau = lambda k: (lambda u: u + k)
sig = lambda k: (lambda u: k - u)
# gr-l3-quotients-engendre-ex4 question 2 : tables de composition.
assert sp.expand(tau(k_)(tau(l_)(x)) - tau(k_ + l_)(x)) == 0
assert sp.expand(tau(k_)(sig(l_)(x)) - sig(k_ + l_)(x)) == 0
assert sp.expand(sig(k_)(tau(l_)(x)) - sig(k_ - l_)(x)) == 0
assert sp.expand(sig(k_)(sig(l_)(x)) - tau(k_ - l_)(x)) == 0
# ex4 question 4 : dans S_n, deux involutions a, b avec ab d'ordre m engendrent au plus 2m éléments (testé sur S4, S5).
for n in (4, 5):
    G = Sn(n)
    e = ident(n)
    invol = [p for p in G if p != e and comp(p, p) == e]
    for a_, b_ in combinations(invol, 2):
        m = ordre(comp(a_, b_), comp, e)
        assert len(engendre([a_, b_], comp, e)) <= 2 * m

# gr-l3-quotients-engendre-ex1 (et v1, v2) : sous-groupe engendré = PGCD Z ; relations de Bézout citées.
for (u, v, w), g, l, faux_ in [((84, 120, 70), 2, 840, (12, 14)), ((90, 150, 105), 15, 3150, (30, 5)),
                               ((60, 84, 35), 1, 420, (12, 7))]:
    assert gcd(gcd(u, v), w) == g and lcm(u, v, w) == l
    assert all(f != g for f in faux_) and l != g
assert 3 * 84 - 2 * 120 == 12 and 6 * 12 - 70 == 2 and 18 * 84 - 12 * 120 - 70 == 2
assert 2 * 90 - 150 == 30 and 105 - 3 * 30 == 15 and 105 - 6 * 90 + 3 * 150 == 15
assert 3 * 60 - 2 * 84 == 12 and 3 * 12 - 35 == 1 and 9 * 60 - 6 * 84 - 35 == 1
assert gcd(84, 120) == 12 and gcd(84, 70) == 14 and gcd(90, 150) == 30 and gcd(60, 84) == 12 and gcd(84, 35) == 7

# gr-l3-quotients-engendre-ex2 : |<(1234),(13)>| = 8 ; srs^-1 = r^-1. v1 : <(1234),(12)> = S4. v2 : V4.
assert conj(s, r) == inv(r) and len(D) == 8
r1, s1 = cyc(4, [1, 2, 3, 4]), cyc(4, [1, 2])
assert conj(r1, s1) == cyc(4, [2, 3]) and conj(comp(r1, r1), s1) == cyc(4, [3, 4])
assert len(engendre([r1, s1], comp, e4)) == 24
r2, s2 = cyc(4, [1, 2], [3, 4]), cyc(4, [1, 3], [2, 4])
assert comp(r2, s2) == comp(s2, r2) == cyc(4, [1, 4], [2, 3]) and engendre([r2, s2], comp, e4) == V4

# gr-l3-quotients-engendre-ex3 : c^k (1 2) c^-k = (k+1 k+2), <c, (1 2)> = S_n, c = (1 2)(2 3)...(n-1 n).
for n in range(2, 8):
    c = cyc(n, list(range(1, n + 1)))
    t12 = cyc(n, [1, 2])
    ck = ident(n)
    for k in range(n - 1):
        assert conj(ck, t12) == cyc(n, [k + 1, k + 2])
        if k:  # leurres (k k+1), (1 k+2), (1 2) faux pour k >= 1 (question du QCM : n >= 3, 1 <= k <= n - 2)
            assert conj(ck, t12) != cyc(n, [k, k + 1]) and conj(ck, t12) != cyc(n, [1, k + 2])
            assert conj(ck, t12) != t12
        else:  # pour k = 0, (1 k+2) et (1 2) coïncident avec la réponse : d'où la restriction k >= 1 du QCM
            assert conj(ck, t12) == cyc(n, [1, k + 2]) == t12
        ck = comp(c, ck)
    assert comp(*[cyc(n, [i, i + 1]) for i in range(1, n)]) == c
    if n <= 5:
        assert len(engendre([c, t12], comp, ident(n))) == len(Sn(n))
# question 3 : morphismes S_n -> C* : seules les deux applications 1 et signature (force brute sur S3, S4 :
# une telle application est déterminée par f(c), f(tau) ; on teste toutes les valeurs possibles dans U_12).
for n in (3, 4):
    G, e = Sn(n), ident(n)
    c, t12 = cyc(n, list(range(1, n + 1))), cyc(n, [1, 2])
    # écriture de chaque élément comme mot en c, t12 (BFS)
    mot = {e: ()}
    file = [e]
    while file:
        g = file.pop(0)
        for lettre, h in (("c", c), ("t", t12)):
            y = comp(g, h)
            if y not in mot:
                mot[y] = mot[g] + (lettre,)
                file.append(y)
    homs = []
    for fc, ft in product(range(12), repeat=2):  # valeurs e^{2iπ j/12}, en exposants mod 12
        f = {g: sum(fc if l == "c" else ft for l in mot[g]) % 12 for g in G}
        if all(f[comp(g, h)] == (f[g] + f[h]) % 12 for g in G for h in G):
            homs.append(f)
    assert len(homs) == 2 and all(set(f.values()) <= {0, 6} for f in homs)
    assert any(all(f[g] == (0 if signe(g) == 1 else 6) for g in G) for f in homs)

# gr-l3-quotients-engendre-ex5 (et v1, v2), affirmation 3, q8 : sous-groupes de Q.
for (x1, x2), q, d, leurres in [((F(1, 6), F(1, 10)), 30, F(1, 30), [F(1, 60), F(1, 15), F(1, 2)]),
                               ((F(1, 4), F(1, 6)), 12, F(1, 12), [F(1, 24), F(1, 2), F(1, 6)]),
                               ((F(2, 9), F(1, 6)), 18, F(1, 18), [F(1, 54), F(1, 9), F(1, 3)])]:
    a1, a2 = int(x1 * q), int(x2 * q)
    assert x1 * q == a1 and x2 * q == a2 and F(gcd(a1, a2), q) == d
    # LZ = dZ (L, d > 0) si et seulement si L = d : chaque leurre est un autre sous-groupe
    assert all(L != d for L in leurres)
assert 2 * F(1, 6) - 3 * F(1, 10) == F(1, 30) and F(1, 4) - F(1, 6) == F(1, 12) and F(2, 9) - F(1, 6) == F(1, 18)
assert (F(1, 6) / F(1, 15)).denominator != 1 and (F(1, 4) / F(1, 6)).denominator != 1 and (F(1, 6) / F(1, 9)).denominator != 1
assert (F(1, 6) / F(1, 2)).denominator != 1 and (F(1, 6) / F(1, 3)).denominator != 1
assert (F(1, 30) / F(1, 60)).denominator == 1 and (F(1, 60) / F(1, 30)).denominator != 1

# ===========================================================================
# Carte gr-l3-quotients-distingues
# ===========================================================================
inv_p = inv
ss3 = tous_sous_groupes(S3, comp, e3)
ss4 = tous_sous_groupes(S4, comp, e4)
ssD = tous_sous_groupes(list(D), comp, e4)
assert len(ss3) == 6 and len(ss4) == 30 and len(ssD) == 10
nd3 = [H for H in ss3 if distingue(H, S3, comp, inv_p)]
nd4 = [H for H in ss4 if distingue(H, S4, comp, inv_p)]
ndD = [H for H in ssD if distingue(H, D, comp, inv_p)]
# gr-l3-quotients-distingues-ex4 (question 3) : sous-groupes distingués de S3 = {id}, A3, S3.
assert sorted(map(len, nd3)) == [1, 3, 6] and A3 in nd3
# gr-l3-quotients-distingues-ex3 v1 : S4 a 4 sous-groupes distingués ({id}, V4, A4, S4) parmi 30.
assert sorted(map(len, nd4)) == [1, 4, 12, 24] and V4 in nd4 and A4 in nd4
# v1 : sommes 1 + parties de {6, 3, 8, 6} et diviseurs de 24.
sommes = {1 + sum(c) for k in range(5) for c in combinations([6, 3, 8, 6], k)}
assert sommes == {1, 4, 7, 9, 10, 12, 13, 15, 16, 18, 21, 24}
assert sorted(v for v in sommes if 24 % v == 0) == [1, 4, 12, 24]
from collections import Counter
cnt = Counter(1 + sum(c) for k in range(5) for c in combinations([6, 3, 8, 6], k))
assert cnt[4] == cnt[12] == cnt[24] == cnt[1] == 1
# types de cycles de S4 : 1, 6, 3, 8, 6 ; deux permutations de même type sont conjuguées.


def type_cycles(p):
    n, vu, t = len(p), set(), []
    for i in range(1, n + 1):
        if i not in vu:
            j, L = i, 0
            while j not in vu:
                vu.add(j)
                j = p[j - 1]
                L += 1
            t.append(L)
    return tuple(sorted(t))


assert sorted(Counter(map(type_cycles, S4)).values()) == [1, 3, 6, 6, 8]
for p in S4:
    for q in S4:
        if type_cycles(p) == type_cycles(q):
            assert any(conj(g, p) == q for g in S4)
# gr-l3-quotients-distingues-ex3 : D a 10 sous-groupes, dont 6 distingués ; liste admise dans l'énoncé.
assert len(ndD) == 6
r2_ = comp(r, r)
liste = {frozenset([e4]), D, frozenset(engendre([r], comp, e4)), frozenset([e4, r2_, s, comp(r2_, s)]),
         frozenset([e4, r2_, comp(r, s), comp(r, r, r, s)])}
liste |= {frozenset([e4, z]) for z in [r2_, s, comp(r, s), comp(r2_, s), comp(r, r, r, s)]}
assert liste == set(ssD)
assert frozenset([e4, r2_, comp(r, s), comp(r, r, r, s)]) == V4  # remarque : un des sous-groupes d'ordre 4 est V4
for k in range(4):  # r (r^k s) r^-1 = r^(k+2) s
    rk = comp(*([r] * k)) if k else e4
    assert conj(r, comp(rk, s)) == comp(comp(*([r] * ((k + 2) % 4))) if (k + 2) % 4 else e4, s)
assert comp(s, inv(r)) == comp(r, s)
# gr-l3-quotients-distingues-ex3 v2 : Q8 (matrices), 6 sous-groupes, tous distingués.
I2 = sp.eye(2)
A_ = sp.Matrix([[sp.I, 0], [0, -sp.I]])
B_ = sp.Matrix([[0, 1], [-1, 0]])
C_ = A_ * B_
assert A_ ** 2 == -I2 and B_ ** 2 == -I2 and C_ ** 2 == -I2 and B_ * A_ == -C_ and C_ == sp.Matrix([[0, sp.I], [sp.I, 0]])
Q8l = [sg * M for sg in (1, -1) for M in (I2, A_, B_, C_)]
cle = lambda M: tuple(M)
Q8 = [cle(M) for M in Q8l]
mat = {cle(M): M for M in Q8l}
mulQ = lambda u, v: cle(mat[u] * mat[v])
invQ = lambda u: cle(mat[u].inv())
assert all(mulQ(u, v) in Q8 for u in Q8 for v in Q8) and len(set(Q8)) == 8
ssQ = tous_sous_groupes(Q8, mulQ, cle(I2))
assert len(ssQ) == 6 and all(distingue(H, Q8, mulQ, invQ) for H in ssQ)
assert sum(1 for u in Q8 if u != cle(I2) and mulQ(u, u) == cle(I2)) == 1  # seul -I2 est d'ordre 2
assert mulQ(cle(A_), cle(B_)) != mulQ(cle(B_), cle(A_))  # Q8 non abélien

# gr-l3-quotients-distingues-ex1 : seul <(1 2 3)> est distingué dans S3 parmi les quatre.
cands = [engendre([cyc(3, [1, 2])], comp, e3), engendre([cyc(3, [1, 3])], comp, e3),
         engendre([cyc(3, [2, 3])], comp, e3), engendre([cyc(3, [1, 2, 3])], comp, e3)]
assert [distingue(H, S3, comp, inv) for H in cands] == [False, False, False, True]
assert conj(cyc(3, [1, 3]), cyc(3, [1, 2])) == cyc(3, [2, 3]) and conj(cyc(3, [1, 2]), cyc(3, [1, 3])) == cyc(3, [2, 3])
assert conj(cyc(3, [1, 2]), cyc(3, [2, 3])) == cyc(3, [1, 3])
# v1 : dans S4, seul V4 parmi V4, <(12)(34)>, <(1234)>, stabilisateur de 4.
stab4 = frozenset(p for p in S4 if p[3] == 4)
cands = [V4, engendre([cyc(4, [1, 2], [3, 4])], comp, e4), engendre([r], comp, e4), stab4]
assert [distingue(H, S4, comp, inv) for H in cands] == [True, False, False, False]
assert conj(cyc(4, [1, 3]), cyc(4, [1, 2], [3, 4])) == cyc(4, [1, 4], [2, 3])
assert conj(cyc(4, [1, 2]), r) == cyc(4, [1, 3, 4, 2]) and conj(cyc(4, [1, 4]), cyc(4, [1, 2])) == cyc(4, [2, 4])
# v2 : dans D, seul <r> parmi <r>, <(13)>, <(12)(34)>, <(14)(23)>.
cands = [engendre([r], comp, e4), engendre([s], comp, e4), engendre([cyc(4, [1, 2], [3, 4])], comp, e4),
         engendre([cyc(4, [1, 4], [2, 3])], comp, e4)]
assert [distingue(H, D, comp, inv) for H in cands] == [True, False, False, False]
assert conj(r, s) == cyc(4, [2, 4]) and conj(r, cyc(4, [1, 2], [3, 4])) == cyc(4, [1, 4], [2, 3])
assert conj(r, cyc(4, [1, 4], [2, 3])) == cyc(4, [1, 2], [3, 4])

# gr-l3-quotients-distingues-ex2 (et v1, v2) : conjugués de sous-groupes de S4.
for h, sig_, att, leurres in [
        (cyc(4, [1, 2, 3]), cyc(4, [3, 4]), cyc(4, [1, 2, 4]), [cyc(4, [1, 2, 3]), cyc(4, [1, 3, 4]), cyc(4, [2, 3, 4])]),
        (cyc(4, [1, 2, 3]), cyc(4, [1, 4]), cyc(4, [2, 3, 4]), [cyc(4, [1, 2, 3]), cyc(4, [1, 2, 4]), cyc(4, [1, 3, 4])]),
        (cyc(4, [1, 2], [3, 4]), cyc(4, [1, 2, 4]), cyc(4, [1, 3], [2, 4]),
         [cyc(4, [1, 2], [3, 4]), cyc(4, [1, 4], [2, 3]), cyc(4, [2, 4])])]:
    H = engendre([h], comp, e4)
    image = frozenset(conj(sig_, y) for y in H)
    assert image == engendre([att], comp, e4) and image != H
    assert all(image != engendre([L], comp, e4) for L in leurres)

# gr-l3-quotients-distingues-ex4 : centre de S_n trivial pour n = 3..6 ; argument tau = (j k).
for n in range(3, 7):
    G = Sn(n)
    assert [p for p in G if all(comp(p, q) == comp(q, p) for q in G)] == [ident(n)]
for p in S4:
    if p != e4:
        i = next(i for i in range(1, 5) if p[i - 1] != i)
        j = p[i - 1]
        k = next(k for k in range(1, 5) if k not in (i, j))
        t_ = cyc(4, [j, k])
        assert conj(t_, p)[i - 1] == k != j
# carte : exemple et contre-exemples ; trous et pièges.
K = frozenset([e4, cyc(4, [1, 2], [3, 4])])
assert distingue(K, V4, comp, inv) and distingue(V4, S4, comp, inv) and not distingue(K, S4, comp, inv)
assert conj(cyc(4, [1, 3]), cyc(4, [1, 2], [3, 4])) == cyc(4, [2, 3], [1, 4])
H12 = engendre([cyc(3, [1, 2])], comp, e3)
assert len(S3) // len(H12) == 3 and not distingue(H12, S3, comp, inv)  # piège « indice 3 », leurre « d'ordre 2 »
# leurre trou 1 « gH = Hg pour au moins un g hors de H » : <(12)> dans S4, g = (34).
H12_4 = engendre([cyc(4, [1, 2])], comp, e4)
g34 = cyc(4, [3, 4])
assert g34 not in H12_4 and {comp(g34, h) for h in H12_4} == {comp(h, g34) for h in H12_4}
assert not distingue(H12_4, S4, comp, inv)
# leurres trou 3 (« fini », « d'ordre pair », « engendré par deux éléments ») : S3.
assert len(engendre([cyc(3, [1, 2]), cyc(3, [1, 2, 3])], comp, e3)) == 6
# affirmations : centre distingué (S4 : centre trivial, D : {id, r^2}) ; (1 2)(1 2 3) = (2 3), (1 2 3)(1 2) = (1 3).
ZD = frozenset(z for z in D if all(comp(z, g) == comp(g, z) for g in D))
assert ZD == {e4, r2_} and distingue(ZD, D, comp, inv)
assert comp(cyc(3, [1, 2]), cyc(3, [1, 2, 3])) == cyc(3, [2, 3]) and comp(cyc(3, [1, 2, 3]), cyc(3, [1, 2])) == cyc(3, [1, 3])
# indice 2 => distingué (tous les sous-groupes d'indice 2 des groupes testés).
for G, ss in ((S3, ss3), (S4, ss4), (list(D), ssD)):
    assert all(distingue(H, G, comp, inv) for H in ss if 2 * len(H) == len(G))

# ===========================================================================
# Carte gr-l3-quotients-quotient
# ===========================================================================


def ordre_classe_zn(n, h, x):
    H = {(h * j) % n for j in range(n)}
    k = 1
    while (k * x) % n not in H:
        k += 1
    return k, len(H)


# exemple : Z/12, H = <4>, |G/H| = 4, classe de 1 d'ordre 4, classe de 8 = H.
assert ordre_classe_zn(12, 4, 1) == (4, 3) and ordre_classe_zn(12, 4, 8)[0] == 1
# gr-l3-quotients-quotient-ex1 (et v1, v2).
for n, h, x, o, oG, cardH, leurres in [(24, 8, 4, 2, 6, 3, [6, 8, 3]), (30, 10, 4, 5, 15, 3, [15, 10, 3]),
                                        (36, 9, 6, 3, 6, 4, [6, 9, 4])]:
    k, c = ordre_classe_zn(n, h, x)
    assert (k, c) == (o, cardH) and n // gcd(x, n) == oG and o not in leurres
    assert n // c == leurres[1]  # le leurre |G/H|
# contre-exemple : Z/4, H = {0, 2}, classe de 1 d'ordre 2 ; affirmation 3.
assert ordre_classe_zn(4, 2, 1)[0] == 2


def ordre_z2(a, b, x, y):
    k = 1
    while (k * x) % a or (k * y) % b:
        k += 1
    return k


# gr-l3-quotients-quotient-ex2 (et v1, v2).
assert ordre_z2(4, 6, 1, 1) == 12 and ordre_z2(6, 10, 1, 1) == 30 and ordre_z2(4, 9, 2, 3) == 6
assert max(ordre_z2(4, 6, x, y) for x in range(4) for y in range(6)) == 12  # G/H (24 éléments) non cyclique
assert max(ordre_z2(6, 10, x, y) for x in range(6) for y in range(10)) == 30  # 60 éléments, non cyclique
assert ordre_z2(4, 9, 1, 1) == 36  # v2 : G/H cyclique
# gr-l3-quotients-quotient-ex3 (et v1, v2), affirmation 4 : ordre de la classe de p/q dans R/Z.


def ordre_rz(x):
    k = 1
    while (k * x).denominator != 1:
        k += 1
    return k


assert ordre_rz(F(21, 12)) == 4 and ordre_rz(F(45, 18)) == 2 and ordre_rz(F(50, 35)) == 7
assert 4 * F(21, 12) == 7 and 2 * F(45, 18) == 5 and 7 * F(50, 35) == 10
assert all(ordre_rz(F(1, n)) == n for n in range(1, 30))
# contre-exemple : H = <(1 2)> non distingué dans S3, loi mal définie.
x_, xp, y_ = cyc(3, [1, 3]), comp(cyc(3, [1, 3]), cyc(3, [1, 2])), cyc(3, [1, 3])
assert {comp(x_, h) for h in H12} == {comp(xp, h) for h in H12}
assert comp(x_, y_) == e3 and comp(xp, y_) == cyc(3, [2, 3]) and comp(xp, y_) not in H12
# piège « loi (xH)(yH) = yxH » : la projection n'est pas un morphisme pour S3/{id} (S3 non abélien).
assert any(comp(p, q) != comp(q, p) for p in S3 for q in S3)
# gr-l3-quotients-quotient-ex5 : S4/V4.
assert len(S4) // len(V4) == 6 and distingue(V4, S4, comp, inv) and distingue(V4, list(A4), comp, inv)
assert comp(cyc(4, [1, 2], [3, 4]), cyc(4, [1, 3], [2, 4])) == cyc(4, [1, 4], [2, 3])
classe = lambda g: frozenset(comp(g, v) for v in V4)
assert classe(r) != classe(e4) and classe(comp(r, r)) == classe(e4)  # ordre 2
assert comp(r, cyc(4, [1, 2], [3, 4])) == cyc(4, [1, 3]) and cyc(4, [1, 3]) in classe(r)
S_ = [p for p in S4 if p[3] == 4]
assert len(S_) == 6 and len({classe(p) for p in S_}) == 6 and set(S_) & V4 == {e4}  # restriction de pi bijective
assert any(classe(comp(p, q)) != classe(comp(q, p)) for p in S4 for q in S4)  # quotient non abélien
assert len(A4) // len(V4) == 3 and classe(cyc(4, [1, 2, 3])) != classe(e4)
# gr-l3-quotients-quotient-ex4 : indices du centre 4 (D), 6 (S3), 24 (S4) ; groupes non abéliens.
assert len(D) // len(ZD) == 4 and any(comp(p, q) != comp(q, p) for p in D for q in D)
classesD = {frozenset(comp(g, z) for z in ZD) for g in D}
mulc = lambda X, Y: frozenset(comp(a_, b_) for a_ in X for b_ in Y)
assert all(mulc(X, X) == ZD for X in classesD)  # D/Z(D) : tous les éléments d'ordre <= 2, non cyclique

# ===========================================================================
# Carte gr-l3-quotients-isomorphisme
# ===========================================================================
# exemple : f(x) = 3x sur Z/12.
assert {x for x in range(12) if 3 * x % 12 == 0} == {0, 4, 8} and {3 * x % 12 for x in range(12)} == {0, 3, 6, 9}
# gr-l3-quotients-isomorphisme-ex1 (et v1, v2) : |ker| et |Im| ; leurres faux.
for n, a, ke, im, leurres in [(24, 18, 6, 4, [(4, 6), (3, 8), (2, 12)]), (30, 12, 6, 5, [(5, 6), (3, 10), (2, 15)]),
                              (36, 8, 4, 9, [(9, 4), (3, 12), (2, 18)])]:
    K_ = [x for x in range(n) if a * x % n == 0]
    I_ = {a * x % n for x in range(n)}
    assert (len(K_), len(I_)) == (ke, im) and ke * im == n and (ke, im) not in leurres
assert [x for x in range(24) if 18 * x % 24 == 0] == [0, 4, 8, 12, 16, 20]
assert sorted({12 * x % 30 for x in range(30)}) == [0, 6, 12, 18, 24]
# gr-l3-quotients-isomorphisme-ex2 (et v1, v2) : noyau de k -> e^{2iπ p k / q}.
for p, q, d, leurres in [(6, 15, 5, [15, 3, 6]), (10, 14, 7, [14, 2, 10]), (9, 12, 4, [12, 3, 9])]:
    noyau_min = next(k for k in range(1, 100) if (F(p, q) * k).denominator == 1)
    assert noyau_min == d == F(p, q).denominator and d not in leurres
# gr-l3-quotients-isomorphisme-ex3 (et v1, v2) : nombre de morphismes Z/m -> Z/n = pgcd (force brute).


def nb_hom(m, n):
    total = 0
    for img in product(range(n), repeat=m):
        if all(img[(x + y) % m] == (img[x] + img[y]) % n for x in range(m) for y in range(m)):
            total += 1
    return total


def nb_hom_rapide(m, n):  # morphismes déterminés par a = f(1) avec m a = 0 (question 1), recoupés ci-dessous
    return sum(1 for a in range(n) if m * a % n == 0)


for m, n in [(4, 6), (6, 4), (3, 6), (6, 9)]:
    assert nb_hom(m, n) == nb_hom_rapide(m, n) == gcd(m, n)
assert nb_hom_rapide(12, 18) == 6 and nb_hom_rapide(20, 15) == 5 and nb_hom_rapide(8, 12) == 4
assert [a for a in range(18) if 12 * a % 18 == 0] == [0, 3, 6, 9, 12, 15]
# gr-l3-quotients-isomorphisme-ex4 (et v1, v2), affirmation 3 : morphismes de S3 dans Z/6, Z/3, S3 (force brute).
idx = {p: i for i, p in enumerate(S3)}


def homs_S3(cible, mul):
    n = 0
    gens = [cyc(3, [1, 2]), cyc(3, [1, 3])]
    for ia, ib in product(range(len(cible)), repeat=2):
        f = {e3: None}
        # extension par mots : f(g h) = f(g) f(h)
        f = {}
        file = [(e3, None)]
        ok = True
        img = {gens[0]: cible[ia], gens[1]: cible[ib]}
        neutre = [c for c in cible if all(mul(c, d) == d for d in cible)][0]
        f[e3] = neutre
        file = [e3]
        while file and ok:
            g = file.pop()
            for h in gens:
                y = comp(g, h)
                v = mul(f[g], img[h])
                if y in f:
                    if f[y] != v:
                        ok = False
                else:
                    f[y] = v
                    file.append(y)
        if ok and all(f[comp(g, h)] == mul(f[g], f[h]) for g in S3 for h in S3):
            n += 1
    return n


assert homs_S3(list(range(6)), lambda u, v: (u + v) % 6) == 2
assert homs_S3(list(range(3)), lambda u, v: (u + v) % 3) == 1
assert homs_S3(S3, comp) == 10
# carte : contre-exemples (Z/4 et (Z/2)^2 ; Z/4 x Z/2 avec H = <(2,0)>, K = <(0,1)>).
G42 = [(a_, b_) for a_ in range(4) for b_ in range(2)]
add42 = lambda u, v: ((u[0] + v[0]) % 4, (u[1] + v[1]) % 2)
Hh, Kk = {(0, 0), (2, 0)}, {(0, 0), (0, 1)}
cl = lambda X, g: frozenset(add42(g, h) for h in X)
quotH = {cl(Hh, g) for g in G42}
quotK = {cl(Kk, g) for g in G42}
assert len(quotH) == len(quotK) == 4
assert all(cl(Hh, add42(g, g)) == frozenset(Hh) for g in G42)  # G/H d'exposant 2 : non cyclique
assert any(all(cl(Kk, (k * g[0] % 4, k * g[1] % 2)) != frozenset(Kk) for k in (1, 2, 3)) for g in G42)  # G/K cyclique
assert max(ordre(x, lambda u, v: (u + v) % 4, 0) for x in range(4)) == 4  # Z/4 cyclique, (Z/2)^2 non
# affirmation 2 : 4 ne divise pas 6.
assert 6 % 4
# piège « application avec f(e) = e' » : x -> x^2, l'image n'est pas stable (2 non carré).
assert all(k * k != 2 for k in range(3))

# ===========================================================================
# Carte gr-l3-quotients-racines : U_n représenté par les exposants k/n (fractions mod 1).
# ===========================================================================


def U(n):
    return {F(k, n) for k in range(n)}


def ordre_u(t):
    return t.denominator


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


# énoncé : générateurs de U_n <=> k ∧ n = 1 ; nombre phi(n) ; n = somme des phi(d) ; inclusion U_d ⊂ U_n <=> d | n.
for n in range(1, 40):
    gens = [k for k in range(n) if ordre_u(F(k, n)) == n]
    assert gens == [k for k in range(n) if gcd(k, n) == 1] and len(gens) == phi(n)
    assert sum(phi(d) for d in range(1, n + 1) if n % d == 0) == n
    for d in range(1, 40):
        assert (U(d) <= U(n)) == (n % d == 0)
# piège / leurre « k premier », « k divise n », « k non multiple de n ».
assert ordre_u(F(2, 6)) == 3 and gcd(1, 6) == 1 and ordre_u(F(2, 4)) == 2
# exemple : n = 12.
assert [k for k in range(1, 13) if gcd(k, 12) == 1] == [1, 5, 7, 11]
assert [phi(d) for d in (1, 2, 3, 4, 6, 12)] == [1, 1, 2, 2, 2, 4] and F(10, 12).denominator == 6
# affirmations : j * j^2 = 1 ; somme des racines primitives 6-ièmes = 1 ; U_m ∩ U_n = U_{m∧n} ; U_2 U_2 = U_2.
assert (F(1, 3) + F(2, 3)) % 1 == 0
assert sp.nsimplify(sp.simplify(sp.exp(sp.I * sp.pi / 3) + sp.exp(-sp.I * sp.pi / 3))) == 1
for m in range(1, 25):
    for n in range(1, 25):
        assert U(m) & U(n) == U(gcd(m, n))
        assert {(a_ + b_) % 1 for a_ in U(m) for b_ in U(n)} == U(lcm(m, n))
assert {(a_ + b_) % 1 for a_ in U(2) for b_ in U(2)} == U(2) != U(4)
# gr-l3-quotients-racines-ex1 (et v1, v2), q5.
assert phi(36) == 12 and phi(45) == 24 and phi(40) == 16 and phi(20) == 8
assert 36 * F(2, 3) == 24 and 36 * F(1, 2) == 18 and 45 * F(2, 3) == 30 and 45 * F(4, 5) == 36
assert 40 * F(1, 2) == 20 and 40 * F(4, 5) == 32 and 20 * F(1, 2) == 10 and phi(5) == 4
assert [k for k in range(1, 21) if gcd(k, 20) == 1] == [1, 3, 7, 9, 11, 13, 17, 19]
# gr-l3-quotients-racines-ex2 (et v1, v2).
for a_, b_, g, l, leurres in [(24, 36, 12, 72, [72, 6, 864]), (20, 30, 10, 60, [60, 5, 600]), (18, 45, 9, 90, [90, 3, 810])]:
    assert U(a_) & U(b_) == U(g) and lcm(a_, b_) == l and all(U(L) != U(g) for L in leurres)
assert 45 - 2 * 18 == 9 and 36 - 24 == 12 and 30 - 20 == 10
# gr-l3-quotients-racines-ex3 (et v1, v2) : sommes des racines primitives (calcul exact SymPy).


def S_prim(n):
    return sp.nsimplify(sp.simplify(sum(sp.exp(2 * sp.I * sp.pi * k / n) for k in range(1, n + 1) if gcd(k, n) == 1)))


assert S_prim(30) == -1 and S_prim(12) == 0 and S_prim(10) == 1
assert [S_prim(n) for n in (1, 2, 3, 4, 5, 6, 15)] == [1, -1, -1, 0, -1, 1, 1]
for n in range(2, 13):  # somme des S(d) pour d | n nulle
    assert sp.simplify(sum(S_prim(d) for d in range(1, n + 1) if n % d == 0)) == 0
assert phi(30) == 8 and phi(12) == 4 and phi(10) == 4
# gr-l3-quotients-racines-ex4 : U_m x U_n cyclique <=> m ∧ n = 1 ; noyau à m ∧ n éléments.
for m in range(1, 13):
    for n in range(1, 13):
        ordre_max = max(lcm(ordre_u(a_), ordre_u(b_)) for a_ in U(m) for b_ in U(n))
        assert (ordre_max == m * n) == (gcd(m, n) == 1)
        noyau = [(a_, b_) for a_ in U(m) for b_ in U(n) if (a_ + b_) % 1 == 0]
        assert len(noyau) == gcd(m, n) and gcd(m, n) * lcm(m, n) == m * n
# leurres : « m et n premiers » (m = n = 2), « m = 1 ou n = 1 » (2, 3), « toujours » (2, 2).
assert gcd(2, 2) != 1 and gcd(2, 3) == 1 and gcd(1, 4) == 1
# gr-l3-quotients-racines-ex5 (et v1, v2) : (Z/pZ)* cyclique, phi(p - 1) générateurs, générateurs cités.
for p, g, puiss, n_gen, leurres in [(13, 2, (4, 6), 4, [12, 6, 3]), (17, 3, (4, 8), 8, [16, 4, 2]),
                                    (19, 2, (6, 9), 6, [18, 9, 3])]:
    ordres = [next(k for k in range(1, p) if pow(x, k, p) == 1) for x in range(1, p)]
    assert ordres.count(p - 1) == phi(p - 1) == n_gen and n_gen not in leurres
    assert ordres[g - 1] == p - 1
    for d in range(1, p):  # psi(d) = phi(d) pour d | p - 1
        if (p - 1) % d == 0:
            assert ordres.count(d) == phi(d)
assert pow(2, 4, 13) == 3 and pow(2, 6, 13) == 12
assert pow(3, 4, 17) == 13 and pow(3, 8, 17) == 16 and 169 % 17 == 16
assert pow(2, 6, 19) == 7 and pow(2, 9, 19) == 18 and 512 % 19 == 18

# ===========================================================================
# QCM
# ===========================================================================
# q1 : <(1 2 3), (1 2)(3 4)> = A4 ; conjugué et produit cités.
assert engendre([cyc(4, [1, 2, 3]), cyc(4, [1, 2], [3, 4])], comp, e4) == A4
assert conj(cyc(4, [1, 2, 3]), cyc(4, [1, 2], [3, 4])) == cyc(4, [2, 3], [1, 4])
assert comp(cyc(4, [1, 2], [3, 4]), cyc(4, [1, 4], [2, 3])) == cyc(4, [1, 3], [2, 4])
# q2 : <(1 2 3 4)> non distingué dans S4, distingué dans D.
R = engendre([r], comp, e4)
assert not distingue(R, S4, comp, inv) and distingue(R, D, comp, inv)
# q3 : ordre de la classe de (1 2 3) dans S4/V4 = 3.
c123 = cyc(4, [1, 2, 3])
assert c123 not in V4 and comp(c123, c123) not in V4 and comp(c123, c123, c123) == e4
# q4 : |Im f| divise 12 et 8, donc 4 ; atteint (k -> 2k de Z/12 dans Z/8 a une image d'ordre 4).
assert gcd(12, 8) == 4 and len({2 * k % 8 for k in range(12)}) == 4 and all(24 % 12 == 0 for _ in [0])
assert all(((2 * (x + y)) - (2 * x + 2 * y)) % 8 == 0 for x in range(12) for y in range(12))
assert all((2 * x) % 8 == (2 * (x + 12)) % 8 for x in range(12))  # bien défini sur Z/12
# q6 : S3, A3 cycliques ; (Z/2)^2 non cyclique.
assert len(A3) == 3 and len(S3) // len(A3) == 2
# q7 : det surjectif sur R* (diag(t, 1, ..., 1)), -1 d'ordre 2 dans R*, pas dans R*_+ ni dans (R, +).
t_ = sp.symbols("t", nonzero=True)
assert sp.diag(t_, 1, 1).det() == t_
# q8 : 1/(2q) n'est pas dans (1/q) Z.
assert all((F(1, 2 * q) * q).denominator == 2 for q in range(1, 50))

print("gr_l3_quotients : toutes les vérifications passent.")
