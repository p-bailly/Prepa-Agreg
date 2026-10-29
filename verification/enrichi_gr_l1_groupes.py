"""Vérifications de l'enrichissement du chapitre « gr-l1-groupes » (QCM et exercices ajoutés), par énumération.

Fichiers content/qcm/gr-l1-groupes.json et fin de content/exercices/groupes/l1-gr-l1-groupes.json
(gr-l1-groupes-noyau-ex5, gr-l1-groupes-sous-groupes-ex5, gr-l1-groupes-lagrange-ex5). Chaque assert cite
l'identifiant de la question ou de l'exercice (variantes : « v1 », « v2 ») ; les leurres sont vérifiés FAUX.
Lancer : python3 verification/enrichi_gr_l1_groupes.py (doit finir sans erreur).
"""

import cmath
from itertools import permutations
from math import gcd

from sympy import totient


def ordre_additif(k, n):
    return n // gcd(k, n)


def sous_groupes_Zn(n):
    """Sous-groupes de Z/nZ par énumération des parties engendrées par une ou deux classes (suffit : cyclique)."""
    res = set()
    for a in range(n):
        res.add(frozenset((a * j) % n for j in range(n)))
    return res


# gr-l1-groupes-q1 : ordre de e^{3iπ/4} = 8 ; 4, 3, 6 ne conviennent pas.
z = cmath.exp(3j * cmath.pi / 4)
puissances = [k for k in range(1, 30) if abs(z**k - 1) < 1e-12]
assert puissances[0] == 8
assert abs(z**4 - (-1)) < 1e-12

# gr-l1-groupes-q2 : Z/12Z a 6 sous-groupes (leurres 4 = φ(12), 12, 5 faux). On vérifie par énumération complète
# des parties stables (tout sous-groupe d'un groupe fini est une partie non vide stable par +).
n = 12
stables = []
for masque in range(1, 2**n):
    P = {i for i in range(n) if masque >> i & 1}
    if all((x + y) % n in P for x in P for y in P):
        stables.append(frozenset(P))
assert len(stables) == 6 and set(stables) == sous_groupes_Zn(12)
assert totient(12) == 4

# gr-l1-groupes-q3 : éléments d'ordre 10 dans Z/30Z : 3, 9, 21, 27 (4 = φ(10)).
assert [k for k in range(30) if ordre_additif(k, 30) == 10] == [3, 9, 21, 27]
assert totient(10) == 4

# gr-l1-groupes-q4 : |H| divise 12, est divisible par 2, 2 < |H| < 12 : {4, 6} ; les deux sont atteintes dans Z/12Z.
possibles = [h for h in range(3, 12) if 12 % h == 0 and h % 2 == 0]
assert possibles == [4, 6]
K = frozenset({0, 6})
cardinaux = sorted({len(H) for H in sous_groupes_Zn(12) if K < H and len(H) < 12})
assert cardinaux == [4, 6]

# gr-l1-groupes-q5 (faux) : 2 est de module ≥ 1, son inverse 1/2 non.
assert abs(2) >= 1 and abs(1 / 2) < 1

# gr-l1-groupes-q6 (vrai) : x -> 5x est une bijection de Z/12Z, morphisme.
img = [(5 * x) % 12 for x in range(12)]
assert sorted(img) == list(range(12))
assert all((5 * (x + y)) % 12 == (img[x] + img[y]) % 12 for x in range(12) for y in range(12))

# gr-l1-groupes-q7 (faux) : aucun morphisme injectif de Z/6Z dans Z/9Z. Énumération : un morphisme est déterminé par
# a = f(1) avec 9 | 6a.
morph = [a for a in range(9) if (6 * a) % 9 == 0]
assert morph == [0, 3, 6]
assert all(len({(k * a) % 9 for k in range(6)}) < 6 for a in morph)
assert 9 % 6 != 0

# gr-l1-groupes-q8 : (12Z + 18Z) ∩ 10Z = 30Z ; leurres 2Z, 60Z, 6Z faux.
somme = {(12 * u + 18 * v) for u in range(-20, 21) for v in range(-20, 21)}
assert min(x for x in somme if x > 0) == 6
inter = [x for x in range(1, 200) if x % 6 == 0 and x % 10 == 0]
assert inter[0] == 30
assert 2 not in inter and 60 != inter[0] and 6 not in inter


# gr-l1-groupes-noyau-ex5 : morphismes de Z/mZ dans Z/nZ, énumérés directement comme applications
# f : Z/mZ -> Z/nZ (déterminées par f(1), testées sur toute la table d'addition).
def morphismes(m, n):
    res = []
    for a in range(n):
        f = [(k * a) % n for k in range(m)]
        # bien définie : f(k) = ka doit être compatible avec k + m ≡ k, c'est-à-dire n | ma
        if (m * a) % n != 0:
            continue
        assert all(f[(x + y) % m] == (f[x] + f[y]) % n for x in range(m) for y in range(m))
        res.append(a)
    # contrôle : tout morphisme vérifie f(k) = k f(1), donc il n'y en a pas d'autres que ceux-ci
    return res


def injectifs(m, n):
    return [a for a in morphismes(m, n) if len({(k * a) % n for k in range(m)}) == m]


for m in range(1, 13):
    for n in range(1, 13):
        assert len(morphismes(m, n)) == gcd(m, n)  # question 3
        if n % m != 0:
            assert injectifs(m, n) == []  # question 4
        d = gcd(m, n)
        assert morphismes(m, n) == [a for a in range(n) if a % (n // d) == 0]
# Exercice : m = 12, n = 18 : 6 morphismes (a multiple de 3), aucun injectif.
assert morphismes(12, 18) == [0, 3, 6, 9, 12, 15] and injectifs(12, 18) == []
assert (6 * 3) % 18 == 0  # f_3 envoie la classe de 6 sur 0
# v1 : m = 6, n = 18 : 6 morphismes, dont 2 injectifs (a = 3 et a = 15).
assert len(morphismes(6, 18)) == 6 and injectifs(6, 18) == [3, 15]
# v2 : m = 4, n = 10 : 2 morphismes (a = 0, 5), aucun injectif.
assert morphismes(4, 10) == [0, 5] and injectifs(4, 10) == []
assert (2 * 5) % 10 == 0


# gr-l1-groupes-sous-groupes-ex5 : H ∩ U = racines 2n-ièmes de l'unité, 2n éléments (n = 12 ; v1 : 6 ; v2 : 5).
def compte(n, N=2000):
    """Compte les e^{iθ}, θ = 2πk/N (N multiple de 2n), tels que z^n est réel ; la grille contient toutes les racines
    2n-ièmes de l'unité, et ce sont les seules solutions trouvées."""
    sols = []
    for k in range(N):
        z = cmath.exp(2j * cmath.pi * k / N)
        if abs((z**n).imag) < 1e-9:
            sols.append(k)
    return sols


for n, attendu in ((12, 24), (6, 12), (5, 10)):
    N = 2 * n * 7
    sols = compte(n, N)
    assert len(sols) == attendu
    assert all(abs(cmath.exp(2j * cmath.pi * k / N) ** (2 * n) - 1) < 1e-9 for k in sols)
    # leurres : n (seulement z^n = 1), 2, une infinité
    assert attendu != n and attendu != 2
# Stabilité : (z w^{-1})^n = z^n / w^n réel, testée sur des éléments de H.
z1, z2 = cmath.exp(1j * cmath.pi / 12) * 3, cmath.exp(5j * cmath.pi / 12) * 0.5
reel = lambda w: abs(w.imag) <= 1e-9 * abs(w)
assert reel(z1**12) and reel(z2**12) and reel((z1 / z2) ** 12)


# gr-l1-groupes-lagrange-ex5 : sous-groupe d'indice 2 : x^2 ∈ H pour tout x, et les éléments d'ordre impair sont
# dans H. Vérifié sur S3 ⊃ A3, S4 ⊃ A4, Z/12Z ⊃ <2> (notation additive), Z/6Z × Z/2Z ⊃ Z/6Z × {0}.
def compose(s, t):
    return tuple(s[t[i]] for i in range(len(t)))


def signature(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def ordre_perm(p):
    e = tuple(range(len(p)))
    q, k = p, 1
    while q != e:
        q, k = compose(p, q), k + 1
    return k


for n in (3, 4):
    G = list(permutations(range(n)))
    H = {p for p in G if signature(p) == 1}
    assert len(G) == 2 * len(H)
    for x in G:
        assert compose(x, x) in H
        if ordre_perm(x) % 2 == 1:
            assert x in H
        if x not in H:
            xH = {compose(x, h) for h in H}
            assert xH.isdisjoint(H) and xH | H == set(G)
H = {x for x in range(12) if x % 2 == 0}
for x in range(12):
    assert (2 * x) % 12 in H
    if ordre_additif(x, 12) % 2 == 1:
        assert x in H
G = [(a, b) for a in range(6) for b in range(2)]
H = {(a, 0) for a in range(6)}
for (a, b) in G:
    assert ((2 * a) % 6, 0) in H

# --- Relecture indépendante (2026-10-29) : recalculs à partir des seuls énoncés ---
from sympy import totient as _phi

# gr-l1-groupes-noyau-ex5 : pour 1 ≤ m, n ≤ 20, les morphismes Z/mZ -> Z/nZ (applications k -> ka avec n | ma,
# testées sur toute la table) sont au nombre de m ∧ n ; les injectifs sont φ(m) si m | n, aucun sinon.
for m_ in range(1, 21):
    for n_ in range(1, 21):
        # une application f est un morphisme ssi f((x + y) mod m) = f(x) + f(y) mod n ; un morphisme est
        # déterminé par a = f(1), et f(k) = ka : on énumère les applications distinctes ainsi obtenues.
        morph = {tuple((k * a) % n_ for k in range(m_)) for a in range(n_)}
        morph = [f_ for f_ in morph if all(f_[(x + y) % m_] == (f_[x] + f_[y]) % n_
                                          for x in range(m_) for y in range(m_))]
        assert len(morph) == gcd(m_, n_)
        inj = [f_ for f_ in morph if len(set(f_)) == m_]
        assert len(inj) == (_phi(m_) if n_ % m_ == 0 else 0)
# Valeurs de la question 5 : (12, 18) -> 6, 0 ; v1 (6, 18) -> 6, 2 ; v2 (4, 10) -> 2, 0.
assert (gcd(12, 18), gcd(6, 18), gcd(4, 10)) == (6, 6, 2) and _phi(6) == 2
# gr-l1-groupes-sous-groupes-ex5 : les z de module 1 avec z^n réel sont les racines 2n-ièmes de l'unité :
# 24 (n = 12), 12 (v1, n = 6), 10 (v2, n = 5) ; contrôle sur une grille fine d'arguments.
for n_, attendu in ((12, 24), (6, 12), (5, 10)):
    N_ = 2 * n_ * 997  # grille d'arguments 2πk/N_, qui contient les racines 2n-ièmes de l'unité
    pts = [k for k in range(N_) if abs((cmath.exp(2j * cmath.pi * k / N_) ** n_).imag) < 1e-9]
    assert len(pts) == attendu and all(k % 997 == 0 for k in pts)
    # leurres : n (racines n-ièmes seulement), 2 (seulement ±1), une infinité : faux
    assert attendu not in (n_, 2)

# gr-l1-groupes-lagrange-ex5 (remarque de la correction) : dans S3, {id, deux permutations circulaires} est un
# sous-groupe d'indice 2 qui contient tous les carrés et les éléments d'ordre impair.
S3 = list(permutations(range(3)))
comp = lambda s, t_: tuple(s[t_[i]] for i in range(3))
idt = (0, 1, 2)
H3 = {idt, (1, 2, 0), (2, 0, 1)}
assert all(comp(a, b) in H3 for a in H3 for b in H3) and len(S3) == 2 * len(H3)
assert all(comp(s, s) in H3 for s in S3)
def _pow(s, k):
    r = idt
    for _ in range(k):
        r = comp(r, s)
    return r
assert all(s in H3 for s in S3 if next(k for k in range(1, 7) if _pow(s, k) == idt) % 2 == 1)

print("enrichi_gr_l1_groupes : tout est vérifié.")
