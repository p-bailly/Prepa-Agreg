"""Vérifications calculatoires du chapitre « gr-l1-groupes » (thème groupes, niveau l1).

Cartes content/cartes/groupes/l1-gr-l1-groupes.json et exercices content/exercices/groupes/l1-gr-l1-groupes.json.
Chaque assert cite l'identifiant de la carte ou de l'exercice vérifié (variante comprise : « v1 », « v2 »). Pour
chaque question à choix, CHAQUE proposition est évaluée par un calcul indépendant (force brute dans Z/nZ, calcul
exact SymPy) : la bonne doit être la seule vraie, les leurres doivent être faux.
Lancer : python3 verification/gr_l1_groupes.py (doit finir sans erreur).
"""

import itertools
import json
import pathlib
import re
from fractions import Fraction
from math import gcd

import sympy
from sympy import I, Matrix, Rational, exp, eye, pi, sqrt, symbols

RACINE = pathlib.Path(__file__).resolve().parent.parent
CARTES = json.loads((RACINE / "content/cartes/groupes/l1-gr-l1-groupes.json").read_text(encoding="utf-8"))["cartes"]
EXOS = json.loads((RACINE / "content/exercices/groupes/l1-gr-l1-groupes.json").read_text(encoding="utf-8"))["exercices"]
EX = {e["id"]: e for e in EXOS}
VUES = set()


def lcm(a, b):
    return a * b // gcd(a, b)


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def ordre_add(x, n):
    """Ordre de la classe de x dans (Z/nZ, +), par force brute."""
    m = 1
    while (m * x) % n:
        m += 1
    return m


def ordre_mult(a, p):
    m, y = 1, a % p
    while y != 1:
        y = y * a % p
        m += 1
    return m


def question(id_, v):
    """(propositions, bonne) de l'exercice id_ (v = 0) ou de sa variante v."""
    e = EX[id_]
    rf = e["reponse_finale"] if v == 0 else e["variantes"][v - 1]["reponse_finale"]
    VUES.add((id_, v))
    return rf["propositions"], rf["bonne"]


def verifier(id_, v, verite):
    """verite : proposition (texte) -> bool, calculé indépendamment. Une seule vraie, à l'indice « bonne »."""
    props, bonne = question(id_, v)
    valeurs = [verite(p) for p in props]
    assert valeurs.count(True) == 1, (id_, v, list(zip(props, valeurs)))
    assert valeurs.index(True) == bonne, (id_, v, props[bonne])


def entier(p):
    m = re.fullmatch(r"\$(\d+)(\\mathbb\{Z\})?\$", p)
    assert m, p
    return int(m.group(1))


def sous_groupes_zn(n):
    """Tous les sous-groupes de Z/nZ (force brute sur les parties contenant 0)."""
    res = []
    for r in range(1, n + 1):
        for S in itertools.combinations(range(n), r):
            S = set(S)
            if 0 in S and all((x - y) % n in S for x in S for y in S):
                res.append(frozenset(S))
    return res


def est_sous_groupe(S, mul, inv, e):
    return bool(S) and all(mul(x, inv(y)) in S for x in S for y in S) and e in S


# ---------- petits groupes : permutations (S3, A4), produit Z/2 x Z/2
def compo(s, t):
    return tuple(s[t[i]] for i in range(len(t)))


def inv_perm(s):
    r = [0] * len(s)
    for i, x in enumerate(s):
        r[x] = i
    return tuple(r)


def signature(s):
    return (-1) ** sum(1 for i in range(len(s)) for j in range(i + 1, len(s)) if s[i] > s[j])


S3 = list(itertools.permutations(range(3)))
A4 = [s for s in itertools.permutations(range(4)) if signature(s) == 1]
assert len(A4) == 12

# =====================================================================  Carte gr-l1-groupes-sous-groupes
FEN = range(-60, 61)
# trou 2, leurre x^{-1}y^{-1} : H = {1} dans Z/3Z vérifie -x-y ∈ H, est non vide, mais n'est pas un sous-groupe
H = {1}
assert all((-x - y) % 3 in H for x in H for y in H) and 0 not in H
# trou 2, leurre x y x^{-1} : H = {a}, a ≠ e dans un groupe commutatif (ici {1} dans Z/3Z) : x + y - x = y ∈ H
assert all((x + y - x) % 3 in H for x in H for y in H)
# trou 2, leurre xy (et piège « stable par produit ») : N est stable par + dans Z mais 0 - 1 ∉ N
assert all(x + y >= 0 for x in range(30) for y in range(30)) and -1 < 0
# piège « couples distincts » : H = {2} dans Z, condition vide, 0 ∉ H
assert not [(x, y) for x in [2] for y in [2] if x != y] and 0 not in [2]
# contre-exemple : {2, 1/2} stable par inverse dans Q*, sans 1
assert {1 / Fraction(x) for x in [Fraction(2), Fraction(1, 2)]} == {Fraction(2), Fraction(1, 2)}
# contre-exemple / affirmation 1 : 2Z ∪ 3Z non stable (3 - 2 = 1)
U = {n for n in FEN if n % 2 == 0 or n % 3 == 0}
assert 2 in U and 3 in U and 1 not in U
# affirmation 4 (parties finies stables par produit) : vérifiée par force brute dans Z/6Z et S3
for n in range(1, 9):
    for r in range(1, n + 1):
        for S in itertools.combinations(range(n), r):
            S = set(S)
            if all((x + y) % n in S for x in S for y in S):
                assert 0 in S and all((-x) % n in S for x in S), (n, S)
for r in range(1, 7):
    for S in itertools.combinations(S3, r):
        S = set(S)
        if all(compo(x, y) in S for x in S for y in S):
            assert est_sous_groupe(S, compo, inv_perm, (0, 1, 2)), S
# exemple : U est un sous-groupe de C* (|z w^{-1}| = 1), contrôle sur des racines de l'unité
z, w = exp(2 * I * pi / 5), exp(2 * I * pi / 7)
assert sympy.simplify(abs(z / w) - 1) == 0

# gr-l1-groupes-sous-groupes-ex1 (et v1, v2) : on évalue chaque partie sur une fenêtre de Z
def est_sg_fenetre(pred):
    S = [n for n in FEN if pred(n)]
    if 0 not in S:
        return False
    # sous-groupe <=> égal à nZ pour n = plus petit élément > 0 (ou {0}) ; sinon un couple x - y sort de S
    pos = [n for n in S if n > 0]
    if not pos:
        return S == [0]
    n0 = min(pos)
    return all(pred(k) == (k % n0 == 0) for k in range(-200, 201))

PARTIES = {
    0: {"A": lambda n: n % 6 == 0, "B": lambda n: n % 6 == 1, "C": lambda n: n >= 0, "D": lambda n: n % 4 == 0 or n % 6 == 0},
    1: {"A": lambda n: n % 10 == 0, "B": lambda n: n % 3 == 0 and n % 5 == 0, "C": lambda n: abs(n) <= 100, "D": lambda n: n % 3 == 0 or n % 5 == 0},
    2: {"A": lambda n: n % 4 == 0, "B": lambda n: n % 4 == 2, "C": lambda n: n <= 0, "D": lambda n: n % 2 == 0 and n % 3 == 0},
}
# contre-exemples cités : C = {|n| ≤ 100} (v1) : 100 - (-1) = 101 ∉ C ; D (base) : 6 - 4 = 2 ∉ D ; D (v1) : 5 - 3 = 2
assert abs(101) > 100 and not PARTIES[0]["D"](2) and not PARTIES[1]["D"](2) and not PARTIES[2]["C"](1)
for v, parts in PARTIES.items():
    vrais = {k for k, pred in parts.items() if est_sg_fenetre(pred)}
    def verite(p, vrais=vrais):
        return set(re.findall(r"\$([A-D])\$", p)) == vrais
    verifier("gr-l1-groupes-sous-groupes-ex1", v, verite)
assert {k for k, p in PARTIES[0].items() if est_sg_fenetre(p)} == {"A"}

# =====================================================================  Carte gr-l1-groupes-noyau
# contre-exemple : g(n) = n^2 ne s'annule qu'en 0, n'est pas injective, n'est pas un morphisme
assert [n for n in FEN if n * n == 0] == [0] and 1 ** 2 == (-1) ** 2 and (1 + 1) ** 2 != 1 + 1
# n -> 2n : injectif (noyau {0}), non surjectif (1 n'est pas atteint)
assert [n for n in FEN if 2 * n == 0] == [0] and all(2 * n != 1 for n in FEN)
# exemple : noyau de z -> z^n a n éléments (n = 2..6)
X = symbols("X")
for n in range(2, 7):
    assert len(sympy.roots(X ** n - 1, X)) == n

# gr-l1-groupes-noyau-ex1 : noyau de t -> e^{i c t} : c t ∈ 2πZ ; on teste l'appartenance d'éléments témoins
def noyau_c(c):
    return lambda t: sympy.simplify(exp(I * c * t) - 1) == 0
NOY = {0: 1, 1: 3, 2: pi / 2}
ENS = {  # texte -> (générateur g de gZ) ; {0} et {0, 2π} traités à part
    r"$2\pi\mathbb{Z}$": 2 * pi, r"$\pi\mathbb{Z}$": pi, r"$\frac{2\pi}{3}\mathbb{Z}$": 2 * pi / 3, r"$6\pi\mathbb{Z}$": 6 * pi,
    r"$\frac{\pi}{3}\mathbb{Z}$": pi / 3, r"$4\mathbb{Z}$": 4, r"$2\mathbb{Z}$": 2, r"$\frac{1}{4}\mathbb{Z}$": Rational(1, 4), r"$4\pi\mathbb{Z}$": 4 * pi,
}
for v, c in NOY.items():
    k = noyau_c(c)
    vrai_g = 2 * pi / c  # générateur positif du noyau
    assert k(vrai_g) and not k(vrai_g / 2) and not k(vrai_g / 3)
    def verite(p, vrai_g=vrai_g):
        if p in (r"$\{0\}$", r"$\{0, 2\pi\}$"):
            return False  # le noyau contient 2π/c·k pour tout k : il est infini
        return sympy.simplify(ENS[p] - vrai_g) == 0
    verifier("gr-l1-groupes-noyau-ex1", v, verite)

# gr-l1-groupes-noyau-ex2 : ensembles de solutions de z^n = b, comparés à sympy.roots
def ens(liste):
    """Liste triée des points (partie réelle, partie imaginaire) arrondis, pour comparer des ensembles finis de complexes."""
    return sorted((round(complex(sympy.N(x)).real, 9), round(complex(sympy.N(x)).imag, 9)) for x in liste)


def proche(A, B):
    return len(A) == len(B) and all(abs(a[0] - b[0]) < 1e-7 and abs(a[1] - b[1]) < 1e-7 for a, b in zip(A, B))
EQ = {
    0: (6, 64, {
        r"$\{2e^{ik\pi/3} : 0 \leq k \leq 5\}$": [2 * exp(I * k * pi / 3) for k in range(6)],
        r"$\{-2, 2\}$": [-2, 2],
        r"$\{2e^{ik\pi/6} : 0 \leq k \leq 5\}$": [2 * exp(I * k * pi / 6) for k in range(6)],
        r"$\{2e^{ik\pi/6} : 0 \leq k \leq 11\}$": [2 * exp(I * k * pi / 6) for k in range(12)]}),
    1: (4, -81, {
        r"$\{3e^{i(2k+1)\pi/4} : 0 \leq k \leq 3\}$": [3 * exp(I * (2 * k + 1) * pi / 4) for k in range(4)],
        r"$\{3e^{ik\pi/2} : 0 \leq k \leq 3\}$": [3 * exp(I * k * pi / 2) for k in range(4)],
        r"$\{3i, -3i\}$": [3 * I, -3 * I],
        r"$\{3e^{i(2k+1)\pi/8} : 0 \leq k \leq 3\}$": [3 * exp(I * (2 * k + 1) * pi / 8) for k in range(4)]}),
    2: (3, 8 * I, {
        r"$\{-2i, \sqrt{3} + i, -\sqrt{3} + i\}$": [-2 * I, sqrt(3) + I, -sqrt(3) + I],
        r"$\{2i, \sqrt{3} - i, -\sqrt{3} - i\}$": [2 * I, sqrt(3) - I, -sqrt(3) - I],
        r"$\{-2i\}$": [-2 * I],
        r"$\{\sqrt{3} + i, -\sqrt{3} + i\}$": [sqrt(3) + I, -sqrt(3) + I]}),
}
for v, (n, b, props) in EQ.items():
    sol = ens(sympy.roots(X ** n - b, X).keys())
    assert len(sol) == n
    verifier("gr-l1-groupes-noyau-ex2", v, lambda p, props=props, sol=sol: proche(ens(props[p]), sol))
# valeurs citées : (2e^{iπ/6})^6 = -64 ; (3i)^4 = 81 ; (-2i)^3 = 8i ; -2i j = √3 + i ; 2i, √3 - i solutions de z^3 = -8i
assert sympy.simplify((2 * exp(I * pi / 6)) ** 6 + 64) == 0 and (3 * I) ** 4 == 81 and sympy.expand((-2 * I) ** 3) == 8 * I
j = Rational(-1, 2) + I * sqrt(3) / 2
assert sympy.expand(-2 * I * j) == sqrt(3) + I and sympy.expand(-2 * I * j ** 2) == -sqrt(3) + I
assert all(sympy.expand(u ** 3) == -8 * I for u in [2 * I, sqrt(3) - I, -sqrt(3) - I])

# gr-l1-groupes-noyau-ex3 : f(a, b) = p a + q b ; noyau et image par force brute
NI = {0: (4, 6), 1: (6, 15), 2: (10, -4)}
for v, (p, q) in NI.items():
    noyau = {(a, b) for a in range(-30, 31) for b in range(-30, 31) if p * a + q * b == 0}
    image = {p * a + q * b for a in range(-30, 31) for b in range(-30, 31)}
    g = gcd(p, q)
    assert all(m in image for m in range(-20, 21) if m % g == 0) and all(m % g == 0 for m in image)
    def verite(prop, noyau=noyau, g=g):
        mi = re.search(r"f\(\\mathbb\{Z\}\^2\) = (\d+)\\mathbb\{Z\}", prop)
        mk = re.search(r"\\ker f = \\\{\((-?\d*)k, (-?\d*)k\)", prop)
        ok_im = int(mi.group(1)) == g
        if mk is None:
            assert r"\ker f = \{(0, 0)\}" in prop
            ok_k = noyau == {(0, 0)}
        else:
            u, w = (int(t) for t in mk.groups())
            ok_k = {(u * k, w * k) for k in range(-60, 61)} & {(a, b) for a in range(-30, 31) for b in range(-30, 31)} == noyau
        return ok_im and ok_k
    verifier("gr-l1-groupes-noyau-ex3", v, verite)
# valeurs citées : f(-1, 1) = 2 ; f(-2, 1) = 3 (v1) ; f(1, 2) = 2 (v2) ; f(2, -5) = 40 (v2)
assert 4 * -1 + 6 == 2 and 6 * -2 + 15 == 3 and 10 - 8 == 2 and 10 * 2 - 4 * -5 == 40

# gr-l1-groupes-noyau-ex4 Q4 : 4 ne divise pas 6
assert 6 % 4 != 0

# =====================================================================  Carte gr-l1-groupes-lagrange
SG6 = sous_groupes_zn(6)
assert sorted(sorted(S) for S in SG6) == [[0], [0, 1, 2, 3, 4, 5], [0, 2, 4], [0, 3]]
assert {frozenset((x + h) % 6 for h in (0, 3)) for x in range(6)} == {frozenset({0, 3}), frozenset({1, 4}), frozenset({2, 5})}
# contre-exemples : {0, 1} et {1, 2} dans Z/3Z, de cardinal 2 qui ne divise pas 3 ; -1 = 2
assert 3 % 2 != 0 and (-1) % 3 == 2
# A4 n'a aucun sous-groupe d'ordre 6 (énumération des parties à 6 éléments contenant l'identité)
e4 = (0, 1, 2, 3)
autres = [s for s in A4 if s != e4]
assert not any(est_sous_groupe(set(S) | {e4}, compo, inv_perm, e4) for S in itertools.combinations(autres, 5))
# trou 3, leurres : Klein a 3 sous-groupes d'ordre 2 (indice 2) ; Z/4Z, H = {0} : indice 4, 2 générateurs
K4 = list(itertools.product(range(2), repeat=2))
add2 = lambda x, y: ((x[0] + y[0]) % 2, (x[1] + y[1]) % 2)
sg2 = [S for S in itertools.combinations(K4, 2) if est_sous_groupe(set(S), add2, lambda x: x, (0, 0))]
assert len(sg2) == 3 and 4 // 2 == 2 and phi(4) == 2
# trou 1 et affirmation : H = {e} dans Z/2Z : |H| = 1 alors que |G|/|H| = 2
assert 2 // 1 == 2

# gr-l1-groupes-lagrange-ex1 : un entier est possible ssi il divise n (réalisé par les multiples de k dans Z/nZ)
LG = {0: (36, 3, 12), 1: (40, 5, 8), 2: (45, 3, 15)}
for v, (n, k, o) in LG.items():
    assert len({k * j % n for j in range(n)}) == o  # multiples de k : sous-groupe (image de j -> kj) à o éléments
    verifier("gr-l1-groupes-lagrange-ex1", v, lambda p, n=n: n % entier(p) == 0)

# gr-l1-groupes-lagrange-ex4 : |H∩K| divise pgcd(|H|, |K|) et est ≥ |H||K|/|G| -> une seule valeur possible
LI = {0: (36, 12, 18, 3, 2), 1: (60, 20, 15, 3, 4), 2: (48, 16, 24, 3, 2)}
for v, (n, h, k, gh, gk) in LI.items():
    possibles = [x for x in range(1, min(h, k) + 1) if h % x == 0 and k % x == 0 and x * n >= h * k]
    assert len(possibles) == 1
    # la situation existe dans Z/nZ : multiples de gh (ordre h) et de gk (ordre k)
    Hs = {gh * j % n for j in range(n)}
    Ks = {gk * j % n for j in range(n)}
    assert len(Hs) == h and len(Ks) == k and len(Hs & Ks) == possibles[0]
    verifier("gr-l1-groupes-lagrange-ex4", v, lambda p, x=possibles[0]: entier(p) == x)
    # formule |HK| |H∩K| = |H| |K| vérifiée sur l'exemple
    assert len({(a + b) % n for a in Hs for b in Ks}) * len(Hs & Ks) == h * k
# formule du produit vérifiée dans S3 pour tous les couples de sous-groupes
SGS3 = [set(S) for r in range(1, 7) for S in itertools.combinations(S3, r) if est_sous_groupe(set(S), compo, inv_perm, (0, 1, 2))]
for Hh in SGS3:
    for Kk in SGS3:
        assert len({compo(a, b) for a in Hh for b in Kk}) * len(Hh & Kk) == len(Hh) * len(Kk)

# =====================================================================  Carte gr-l1-groupes-ordre
assert ordre_add(8, 12) == 3 and ordre_mult(2, 7) == 3 and ordre_mult(3, 7) == 6
assert [pow(2, k, 7) for k in (2, 3)] == [4, 1] and [pow(3, k, 7) for k in (2, 3)] == [2, 6]
A = Matrix([[0, -1], [1, 0]])
B = Matrix([[0, 1], [-1, -1]])
assert A ** 4 == eye(2) and A ** 2 != eye(2) and B ** 3 == eye(2) and B != eye(2)
assert A * B == Matrix([[1, 1], [0, 1]]) and A * B != B * A
assert all((A * B) ** k == Matrix([[1, k], [0, 1]]) for k in range(1, 30))
assert (-1) ** 4 == 1 and (-1) ** 2 == 1 and 2 % 4 != 0  # (-1)^2 = 1 alors que 4 ∤ 2
# méthode du cours : a^d = e et a^{d/p} ≠ e pour tout premier p | d  <=>  ordre d (testé dans Z/nZ, n ≤ 60)
for n in range(1, 61):
    for x in range(n):
        o = ordre_add(x, n)
        for d in [d for d in range(1, n + 1) if n % d == 0]:
            critere = (d * x) % n == 0 and all((d // p * x) % n != 0 for p in sympy.primefactors(d))
            assert critere == (o == d)

# gr-l1-groupes-ordre-ex1
for v, (x, n) in {0: (18, 30), 1: (12, 42), 2: (20, 36)}.items():
    o = ordre_add(x, n)
    assert o == n // gcd(x, n)
    verifier("gr-l1-groupes-ordre-ex1", v, lambda p, o=o: entier(p) == o)
# gr-l1-groupes-ordre-ex2 (+ puissances citées dans les corrections)
for v, (a, p_) in {0: (2, 13), 1: (3, 17), 2: (5, 11)}.items():
    o = ordre_mult(a, p_)
    verifier("gr-l1-groupes-ordre-ex2", v, lambda p, o=o: entier(p) == o)
assert [pow(2, k, 13) for k in (2, 3, 4, 6)] == [4, 8, 3, 12] and ordre_mult(2, 13) == 12
assert [pow(3, k, 17) for k in (2, 4, 8)] == [9, 13, 16] and 169 == 10 * 17 - 1 and ordre_mult(3, 17) == 16
assert pow(5, 2, 11) == 3 and pow(5, 5, 11) == 1 and 45 == 4 * 11 + 1 and ordre_mult(5, 11) == 5
# gr-l1-groupes-ordre-ex3 : ordre de a^k pour a d'ordre d (modèle : 1 dans Z/dZ)
for v, (d, k) in {0: (24, 18), 1: (30, 12), 2: (40, 12)}.items():
    o = ordre_add(k, d)
    assert o == d // gcd(d, k)
    verifier("gr-l1-groupes-ordre-ex3", v, lambda p, o=o: entier(p) == o)
assert 36 % 24 != 0 and 12 % 30 != 0 and 60 % 40 != 0 and 24 % 40 != 0
# formule d/(d ∧ k) testée pour d ≤ 40, 1 ≤ k ≤ 80
assert all(ordre_add(k, d) == d // gcd(d, k) for d in range(1, 41) for k in range(1, 81))
# gr-l1-groupes-ordre-ex4 : ordre d'un complexe e^{2iπ r} (r rationnel) = dénominateur de r mod 1
def ordre_racine(r):
    return Fraction(r).limit_denominator(10 ** 6).denominator
PR = {0: (Fraction(1, 4), Fraction(1, 5)), 1: (Fraction(1, 2), Fraction(1, 7)), 2: (Fraction(1, 3), Fraction(1, 8))}
for v, (r1, r2) in PR.items():
    m, n = ordre_racine(r1), ordre_racine(r2)
    assert gcd(m, n) == 1
    o = ordre_racine((r1 + r2) % 1)
    assert o == m * n
    verifier("gr-l1-groupes-ordre-ex4", v, lambda p, o=o: entier(p) == o)
assert (Fraction(1, 4) + Fraction(1, 5)) == Fraction(9, 20) and Fraction(1, 2) + Fraction(1, 7) == Fraction(9, 14) and Fraction(1, 3) + Fraction(1, 8) == Fraction(11, 24)
# les complexes de l'énoncé sont bien ces racines : i e^{2iπ/5}, -e^{2iπ/7}, e^{2iπ/3} e^{iπ/4}
assert sympy.simplify(I * exp(2 * I * pi / 5) - exp(2 * I * pi * Rational(9, 20))) == 0
assert sympy.simplify(-exp(2 * I * pi / 7) - exp(2 * I * pi * Rational(9, 14))) == 0
assert sympy.simplify(exp(2 * I * pi / 3) * exp(I * pi / 4) - exp(2 * I * pi * Rational(11, 24))) == 0
# gr-l1-groupes-ordre-ex5 : un groupe d'ordre pair a un élément d'ordre 2 (contrôle sur Z/nZ, n pair, et S3)
assert all(any(ordre_add(x, n) == 2 for x in range(n)) for n in range(2, 41, 2))
assert any(compo(s, s) == (0, 1, 2) and s != (0, 1, 2) for s in S3)

# =====================================================================  Carte gr-l1-groupes-sous-groupes-z
def ensemble_nZ(pred):
    """Renvoie n si {k : pred(k)} = nZ (sur une fenêtre), sinon None."""
    S = [k for k in range(-300, 301) if pred(k)]
    pos = [k for k in S if k > 0]
    n0 = min(pos) if pos else 0
    ok = all(pred(k) == (k % n0 == 0 if n0 else k == 0) for k in range(-300, 301))
    return n0 if ok else None
assert ensemble_nZ(lambda k: any(k == 12 * u + 18 * w for u in range(-30, 31) for w in range(-30, 31))) == 6
assert ensemble_nZ(lambda k: k % 12 == 0 and k % 18 == 0) == 36
assert ensemble_nZ(lambda k: (4 * k) % 6 == 0) == 3 and 4 % 6 and 8 % 6 and 12 % 6 == 0
assert ensemble_nZ(lambda k: k >= 0) is None  # N n'est pas un nZ
assert 0 not in {-1, 1}
# trou 3, leurres (a, b) = (2, 3) : pgcd 1, ppcm 6, |ab| = 6, |a| + |b| = 5
assert gcd(2, 3) == 1 and lcm(2, 3) == 6 and 2 + 3 == 5
# affirmation 2 : 4Z ⊂ 2Z mais 4 ∤ 2
assert all(k % 2 == 0 for k in range(-100, 101, 4)) and 2 % 4 != 0
# théorème vérifié sur les sous-groupes de Z/NZ « relevés » : tout sous-groupe de Z contenant N est nZ avec n | N
for N in range(1, 25):
    for S in sous_groupes_zn(N):
        n0 = min(x for x in S if x > 0) if len(S) > 1 else N
        assert N % n0 == 0 and S == frozenset(range(0, N, n0))

# gr-l1-groupes-sous-groupes-z-ex1
for v, (a, b) in {0: (84, 120), 1: (90, 126), 2: (70, 98)}.items():
    g = ensemble_nZ(lambda k, a=a, b=b: (k % gcd(a, b)) == 0)
    assert g == gcd(a, b) and any(a * u + b * w == g for u in range(-10, 11) for w in range(-10, 11))
    verifier("gr-l1-groupes-sous-groupes-z-ex1", v, lambda p, g=g: entier(p) == g)
assert 3 * 84 - 2 * 120 == 12 and 3 * 90 - 2 * 126 == 18 and 3 * 70 - 2 * 98 == 14
assert lcm(84, 120) == 840 and lcm(90, 126) == 630 and lcm(70, 98) == 490
assert 120 == 84 + 36 and 84 == 2 * 36 + 12 and 126 == 90 + 36 and 90 == 2 * 36 + 18 and 98 == 70 + 28 and 70 == 2 * 28 + 14
# gr-l1-groupes-sous-groupes-z-ex2
for v, (n, c) in {0: (35, 14), 1: (36, 15), 2: (84, 30)}.items():
    h = ensemble_nZ(lambda k, n=n, c=c: (c * k) % n == 0)
    assert h == n // gcd(n, c)
    verifier("gr-l1-groupes-sous-groupes-z-ex2", v, lambda p, h=h: entier(p) == h)
# gr-l1-groupes-sous-groupes-z-ex3 : sous-groupes nZ (n ≥ 0) contenant N
for v, N in {0: 60, 1: 48, 2: 100}.items():
    nb = sum(1 for n in range(0, N + 1) if (N == 0 if n == 0 else N % n == 0))
    assert nb == len(sympy.divisors(N))
    verifier("gr-l1-groupes-sous-groupes-z-ex3", v, lambda p, nb=nb: entier(p) == nb)
assert sympy.divisors(60) == [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60]
assert sympy.divisors(48) == [1, 2, 3, 4, 6, 8, 12, 16, 24, 48] and sympy.divisors(100) == [1, 2, 4, 5, 10, 20, 25, 50, 100]
# gr-l1-groupes-sous-groupes-z-ex4 : périodes ; T est « période de toutes ces fonctions » ssi g | T
for v, (p_, q_) in {0: (18, 30), 1: (28, 42), 2: (45, 75)}.items():
    g = gcd(p_, q_)
    ind = lambda k, g=g: 1 if k % g == 0 else 0
    periodes = [T for T in range(1, 300) if all(ind(k + T) == ind(k) for k in range(-100, 100))]
    assert periodes[0] == g and all(T % g == 0 for T in periodes) and p_ in periodes and q_ in periodes
    def verite(p, g=g):
        T = entier(p)
        return T % g == 0 and T == min(x for x in range(1, T + 1) if x % g == 0)
    verifier("gr-l1-groupes-sous-groupes-z-ex4", v, verite)
assert 2 * 18 - 30 == 6 and 42 - 28 == 14 and 2 * 45 - 75 == 15 and lcm(45, 75) == 225 and lcm(28, 42) == 84

# =====================================================================  Carte gr-l1-groupes-cycliques
assert [k for k in range(1, 13) if gcd(k, 12) == 1] == [1, 5, 7, 11]
SG12 = sous_groupes_zn(12)
assert len(SG12) == 6 and sorted(len(S) for S in SG12) == [1, 2, 3, 4, 6, 12]
assert all(frozenset(range(0, 12, 12 // d)) in SG12 for d in (1, 2, 3, 4, 6, 12))
assert [k for k in range(1, 7) if gcd(k, 6) == 1] == [1, 5]  # générateurs ω, ω^5 de U6
assert all(add2(x, x) == (0, 0) for x in K4)  # Klein : aucun élément d'ordre 4
assert all(pow(x, 2, 8) == 1 for x in (1, 3, 5, 7))  # (Z/8Z)* non cyclique
mult = [((k * 1) % 2, (k * 1) % 3) for k in range(1, 7)]
assert mult == [(1, 1), (0, 2), (1, 0), (0, 1), (1, 2), (0, 0)]  # (1, 1) d'ordre 6 dans Z/2 x Z/3
assert ordre_add(2, 4) == 2  # 2 n'engendre pas Z/4Z
# théorème vérifié par force brute pour n ≤ 30 : φ(n) générateurs, un unique sous-groupe par diviseur
for n in range(1, 31):
    assert sum(1 for x in range(n) if ordre_add(x, n) == n) == phi(n)
    if n <= 16:
        SG = sous_groupes_zn(n)
        assert sorted(len(S) for S in SG) == sympy.divisors(n)
    assert sum(phi(d) for d in sympy.divisors(n)) == n
    for d in sympy.divisors(n):
        assert sum(1 for x in range(n) if ordre_add(x, n) == d) == phi(d)
# trous : n = 4 : φ = 2 ≠ n - 1, n ; n = 3 : φ = 2 ≠ n/2 ; d = n = 3 : φ(3) = 2 ≠ 1 ; n = 2, d = 1 : n/d = 2 ≠ 1
assert phi(4) == 2 and phi(3) == 2 and Fraction(3, 2) != 2 and ordre_add(2, 4) != 4

# gr-l1-groupes-cycliques-ex1
for v, n in {0: 18, 1: 20, 2: 21}.items():
    nb = sum(1 for x in range(n) if ordre_add(x, n) == n)
    verifier("gr-l1-groupes-cycliques-ex1", v, lambda p, nb=nb: entier(p) == nb)
assert [k for k in range(1, 19) if gcd(k, 18) == 1] == [1, 5, 7, 11, 13, 17] and ordre_add(3, 18) == 6
assert [k for k in range(1, 21) if gcd(k, 20) == 1] == [1, 3, 7, 9, 11, 13, 17, 19] and ordre_add(5, 20) == 4
assert [k for k in range(1, 22) if gcd(k, 21) == 1] == [1, 2, 4, 5, 8, 10, 11, 13, 16, 17, 19, 20] and ordre_add(7, 21) == 3
# gr-l1-groupes-cycliques-ex2 : la proposition <x> est-elle d'ordre d ?
for v, (n, d) in {0: (24, 6), 1: (20, 4), 2: (36, 9)}.items():
    def verite(p, n=n, d=d):
        m = re.fullmatch(r"\$\\langle \\(?:bar (\d+)|overline\{(\d+)\}) \\rangle\$", p)
        x = int(m.group(1) or m.group(2))
        return ordre_add(x, n) == d
    verifier("gr-l1-groupes-cycliques-ex2", v, verite)
assert (ordre_add(6, 24), ordre_add(3, 24), ordre_add(8, 24)) == (4, 8, 3)
assert (ordre_add(4, 20), ordre_add(2, 20), ordre_add(10, 20)) == (5, 10, 2)
assert (ordre_add(9, 36), ordre_add(3, 36), ordre_add(6, 36)) == (4, 12, 6)
# gr-l1-groupes-cycliques-ex3 : nombre d'éléments d'ordre d dans Z/nZ (+ générateurs cités)
for v, (n, d, gen, liste) in {0: (36, 12, 3, [3, 15, 21, 33]), 1: (27, 9, 3, [3, 6, 12, 15, 21, 24]), 2: (60, 15, 4, [4, 8, 16, 28, 32, 44, 52, 56])}.items():
    elts = [x for x in range(n) if ordre_add(x, n) == d]
    assert elts == liste and len(elts) == phi(d) and set(elts) <= set(range(0, n, gen)) and len(set(range(0, n, gen))) == d
    verifier("gr-l1-groupes-cycliques-ex3", v, lambda p, nb=len(elts): entier(p) == nb)
# gr-l1-groupes-cycliques-ex4 : Z/m x Z/n cyclique ssi un élément d'ordre mn (force brute)
def cyclique_produit(m, n):
    return any(lcm(ordre_add(x, m), ordre_add(y, n)) == m * n for x in range(m) for y in range(n))
for v in (0, 1, 2):
    def verite(p):
        m_, n_ = (int(t) for t in re.findall(r"\\mathbb\{Z\}/(\d+)\\mathbb\{Z\}", p))
        assert cyclique_produit(m_, n_) == (gcd(m_, n_) == 1)
        return cyclique_produit(m_, n_)
    verifier("gr-l1-groupes-cycliques-ex4", v, verite)
# exemples cités : dans Z/4 x Z/6 (resp. Z/3 x Z/12, Z/9 x Z/12), ordre maximal 12 (resp. 12, 36)
assert max(lcm(ordre_add(x, 4), ordre_add(y, 6)) for x in range(4) for y in range(6)) == 12
assert max(lcm(ordre_add(x, 3), ordre_add(y, 12)) for x in range(3) for y in range(12)) == 12
assert max(lcm(ordre_add(x, 9), ordre_add(y, 12)) for x in range(9) for y in range(12)) == 36
# gr-l1-groupes-cycliques-ex5 Q4 : (Z/pZ)* cyclique pour p premier ≤ 100 ; puissances de 3 modulo 7
for p_ in sympy.primerange(2, 100):
    assert any(ordre_mult(a, p_) == p_ - 1 for a in range(1, p_))
assert [pow(3, k, 7) for k in range(1, 7)] == [3, 2, 6, 4, 5, 1]

# ---------- toutes les questions à choix du fichier ont été vérifiées
attendues = {(e["id"], 0) for e in EXOS if "reponse_finale" in e} | {
    (e["id"], i + 1) for e in EXOS for i, v in enumerate(e.get("variantes", [])) if "reponse_finale" in v}
assert attendues == VUES, attendues ^ VUES
# chaque carte a 3 à 6 exercices, et les identifiants suivent <carte>-exN
for c in CARTES:
    ids = [e["id"] for e in EXOS if e["carte"] == c["id"]]
    assert 3 <= len(ids) <= 6 and ids == [f"{c['id']}-ex{k}" for k in range(1, len(ids) + 1)], c["id"]
print(f"gr_l1_groupes : {len(VUES)} questions à choix vérifiées, toutes les assertions passent.")
