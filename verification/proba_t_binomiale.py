# Vérification calculatoire du chapitre « proba-t-binomiale » (cartes et exercices des fichiers
# content/cartes/probabilites/terminale-proba-t-binomiale.json et content/exercices/probabilites/terminale-proba-t-binomiale.json).
# Calculs exacts (fractions, entiers), souvent par énumération complète. Chaque assert cite l'identifiant concerné ;
# « v1 », « v2 » désignent les variantes. La fonction rep() vérifie, dans le fichier d'exercices, que la proposition
# marquée juste vaut la valeur calculée et que chaque leurre vaut une valeur différente.
import json
import pathlib
import random
import re
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial, perm

import sympy as sp

RACINE = pathlib.Path(__file__).resolve().parent.parent
EX = {e["id"]: e for e in json.loads((RACINE / "content/exercices/probabilites/terminale-proba-t-binomiale.json").read_text())["exercices"]}
CARTES = {c["id"]: c for c in json.loads((RACINE / "content/cartes/probabilites/terminale-proba-t-binomiale.json").read_text())["cartes"]}


def valeur(tex):
    """Valeur exacte d'une proposition simple : entier, décimal à virgule ou fraction \\frac{a}{b}."""
    s = tex.strip().strip("$").replace("{,}", ".").replace(" ", "")
    m = re.fullmatch(r"(-?)\\frac\{(\d+)\}\{(\d+)\}", s)
    if m:
        return (-1 if m.group(1) else 1) * F(int(m.group(2)), int(m.group(3)))
    if re.fullmatch(r"-?\d+(\.\d+)?", s):
        return F(s)
    raise ValueError(tex)


def rep(ident, k, attendu):
    """Exercice ident, variante k (0 = énoncé principal) : bonne proposition = attendu, leurres différents."""
    e = EX[ident] if k == 0 else EX[ident]["variantes"][k - 1]
    r = e["reponse_finale"]
    props = r["propositions"]
    assert len(set(props)) == 4, f"{ident} v{k} : propositions distinctes"
    assert valeur(props[r["bonne"]]) == attendu, f"{ident} v{k} : bonne réponse {props[r['bonne']]} != {attendu}"
    for i, x in enumerate(props):
        if i != r["bonne"]:
            assert valeur(x) != attendu, f"{ident} v{k} : le leurre {x} est juste"
    return [valeur(x) for i, x in enumerate(props) if i != r["bonne"]]


# ==========================================================================================
# Carte proba-t-binomiale-denombrement
# ==========================================================================================
for n in range(0, 7):
    E = list(range(n))
    for k in range(0, 6):
        # n^k k-uplets ; n!/(n-k)! k-uplets d'éléments distincts (0 si k > n) ; C(n,k) parties
        assert len(list(product(E, repeat=k))) == n ** k, "proba-t-binomiale-denombrement : n^k"
        nb_dist = len(list(permutations(E, k)))
        assert nb_dist == (factorial(n) // factorial(n - k) if k <= n else 0), "denombrement : k-uplets distincts"
        nb_parties = len(list(combinations(E, k)))
        if k <= n:
            assert nb_parties == factorial(n) // (factorial(k) * factorial(n - k)) == comb(n, k), "denombrement : parties"
            assert nb_dist == factorial(k) * nb_parties, "denombrement : a_k = k! c_k"
            # leurres du texte à trous : faux en général
        assert len(list(permutations(E))) == factorial(n), "denombrement : n! permutations"
# Leurres des trous : chacun est faux pour au moins un couple (n,k) avec k <= n
assert 2 ** 3 != 3 ** 2 and 3 * 2 != 3 ** 2, "denombrement trou 1 : k^n et nk différents de n^k (n=3,k=2)"
assert perm(3, 2) != 3 ** 2, "denombrement trou 1 : n!/(n-k)! != n^k (n=3,k=2)"
assert 3 ** 2 != perm(3, 2) and factorial(3) // factorial(1) != perm(3, 1) and comb(3, 2) != perm(3, 2), "denombrement trou 2"
assert perm(3, 2) != comb(3, 2) and F(2 ** 2, 2) != comb(2, 2) and factorial(3) // factorial(1) != comb(3, 1), "denombrement trou 3"
# Exemple : E = {a,b,c,d}, k = 2
assert (4 ** 2, perm(4, 2), factorial(4), comb(4, 2)) == (16, 12, 24, 6), "denombrement exemple"
assert 12 == 2 * 6
# Contre-exemples : 3 parties / 6 couples ; 8 triplets de {0,1}, aucun d'éléments distincts ; 6 multi-ensembles de taille 2
assert comb(3, 2) == 3 and perm(3, 2) == 6, "denombrement contre-exemple ordre"
assert len(list(product([0, 1], repeat=3))) == 8 and len(list(permutations([0, 1], 3))) == 0, "denombrement k > n"
multi = {tuple(sorted(t)) for t in product("abc", repeat=2)}
assert len(multi) == 6 and 6 not in (3 ** 2, comb(3, 2)), "denombrement contre-exemple combinaisons avec répétitions"
# Affirmations : 720 bureaux = 6 x 120 ; parties <-> {0,1}^n ; aucun triplet distinct de {a,b} ; 5^3 = 125
assert perm(10, 3) == 720 == 6 * comb(10, 3) and comb(10, 3) == 120, "denombrement affirmation 1"
for n in range(0, 8):
    parties = {frozenset(c) for k in range(n + 1) for c in combinations(range(n), k)}
    assert len(parties) == 2 ** n == len(list(product([0, 1], repeat=n))), "denombrement affirmation 2"
assert len(list(permutations("ab", 3))) == 0 and len(list(product(range(5), repeat=3))) == 125, "denombrement affirmations 3 et 4"

# ------------------------------------------------------------------------------------------ ex1 (codes)
for k, (taille, longueur, leurres) in enumerate([(10, 4, [10 ** 4, comb(10, 4), factorial(4)]),
                                                  (10, 5, [10 ** 5, comb(10, 5), factorial(5)]),
                                                  (26, 3, [26 ** 3, comb(26, 3), factorial(3)])]):
    juste = len(list(permutations(range(taille), longueur)))
    assert sorted(rep("proba-t-binomiale-denombrement-ex1", k, juste)) == sorted(leurres), f"denombrement-ex1 v{k} : leurres attendus"
    assert taille ** longueur in (10000, 100000, 17576)
assert perm(10, 4) == 5040 and perm(10, 5) == 30240 and perm(26, 3) == 15600, "denombrement-ex1 : valeurs de la correction"

# ------------------------------------------------------------------------------------------ ex2 (bureau / commission)
for k, (n, kk) in enumerate([(25, 3), (30, 4), (20, 3)]):
    bureaux, commissions = perm(n, kk), comb(n, kk)
    assert bureaux == factorial(kk) * commissions, f"denombrement-ex2 v{k} : bureaux = k! commissions"
    assert sorted(rep("proba-t-binomiale-denombrement-ex2", k, commissions)) == sorted([bureaux, n ** kk, n * kk]), f"denombrement-ex2 v{k}"
assert (perm(25, 3), comb(25, 3), perm(30, 4), comb(30, 4), perm(20, 3), comb(20, 3)) == (13800, 2300, 657720, 27405, 6840, 1140)

# ------------------------------------------------------------------------------------------ ex3 (anagrammes)
for k, (mot, leurres) in enumerate([("BALLON", [720, 180, 15]), ("RADAR", [120, 60, 10]), ("ANANAS", [720, 120, 360])]):
    juste = len(set(permutations(mot)))
    assert sorted(rep("proba-t-binomiale-denombrement-ex3", k, juste)) == sorted(leurres), f"denombrement-ex3 v{k}"
assert comb(6, 2) * factorial(4) == 360 == factorial(6) // 2, "denombrement-ex3 : 15 x 24 = 720 / 2"
assert comb(5, 2) * comb(3, 2) == 30 == factorial(5) // 4, "denombrement-ex3 v1"
assert comb(6, 3) * comb(3, 2) == 60 == factorial(6) // 12, "denombrement-ex3 v2"
# leurres : 6!/(2!2!) = 180, C(6,2) = 15 ; 5!/2! = 60 ; 6!/3! = 120, 6!/2! = 360
assert factorial(6) // 4 == 180 and factorial(5) // 2 == 60 and factorial(6) // 6 == 120 and factorial(6) // 2 == 360

# ------------------------------------------------------------------------------------------ ex4 (au moins un as / un cœur)
def mains(N, a, k, brute):
    total, sans = comb(N, k), comb(N - a, k)
    if brute:
        cartes = range(N)  # les cartes 0..a-1 sont les cartes « spéciales »
        nb_A = 0
        nb_total = 0
        multiplicite = 0
        for main in combinations(cartes, k):
            nb_total += 1
            j = sum(1 for c in main if c < a)
            nb_A += j >= 1
            multiplicite += j
        assert nb_total == total and nb_A == total - sans
        # le procédé de l'élève compte chaque main autant de fois qu'elle a de cartes spéciales
        assert multiplicite == a * comb(N - 1, k - 1)
    return F(total - sans, total), total, sans

for k, (N, a, kk, brute) in enumerate([(32, 4, 5, True), (52, 4, 5, False), (32, 8, 3, True)]):
    p, total, sans = mains(N, a, kk, brute)
    leurres = rep("proba-t-binomiale-denombrement-ex4", k, p)
    faux = [F(a * comb(N - 1, kk - 1), total), F(sans, total), 1 - (1 - F(a, N)) ** kk]
    assert sorted(leurres) == sorted(faux), f"denombrement-ex4 v{k} : leurres (double comptage, complémentaire, avec remise)"
    assert a * comb(N - 1, kk - 1) > total - sans, f"denombrement-ex4 v{k} : le décompte de l'élève est trop grand"
assert (comb(32, 5), comb(28, 5), comb(32, 5) - comb(28, 5), 4 * comb(31, 4)) == (201376, 98280, 103096, 125860)
assert F(103096, 201376) == F(1841, 3596) and abs(float(F(1841, 3596)) - 0.512) < 5e-4, "denombrement-ex4"
assert (comb(52, 5), comb(48, 5), comb(52, 5) - comb(48, 5), 4 * comb(51, 4)) == (2598960, 1712304, 886656, 999600)
assert F(886656, 2598960) == F(18472, 54145) and abs(float(F(18472, 54145)) - 0.341) < 5e-4, "denombrement-ex4 v1"
assert (comb(32, 3), comb(24, 3), comb(32, 3) - comb(24, 3), 8 * comb(31, 2)) == (4960, 2024, 2936, 3720)
assert F(2936, 4960) == F(367, 620) and abs(float(F(367, 620)) - 0.592) < 5e-4, "denombrement-ex4 v2"

# ==========================================================================================
# Carte proba-t-binomiale-pascal
# ==========================================================================================
for n in range(0, 25):
    assert comb(n, 0) == comb(n, n) == 1, "pascal : bords"
    if n >= 1:
        assert comb(n, 1) == n
    if n >= 2:
        assert comb(n, 2) == n * (n - 1) // 2
    for k in range(n + 1):
        assert comb(n, k) == comb(n, n - k), "pascal : symétrie"
    for k in range(n):
        assert comb(n, k) + comb(n, k + 1) == comb(n + 1, k + 1), "pascal : relation de Pascal"
        # même identité par le calcul symbolique : réduction au même dénominateur
    assert sum(comb(n, j) for j in range(n + 1)) == 2 ** n, "pascal : somme d'une ligne"
nn, kk = sp.symbols("n k", integer=True, nonnegative=True)
assert sp.simplify(sp.combsimp(sp.binomial(nn, kk) + sp.binomial(nn, kk + 1) - sp.binomial(nn + 1, kk + 1))) == 0, "pascal : calcul"
# Pascal combinatoire : parties à k+1 éléments de F = {0..n} contenant ou non a = n
for n in range(1, 8):
    for k in range(n):
        Fp = [frozenset(c) for c in combinations(range(n + 1), k + 1)]
        assert sum(1 for A in Fp if n in A) == comb(n, k) and sum(1 for A in Fp if n not in A) == comb(n, k + 1)
# Leurres des trous : faux pour au moins une valeur
assert comb(4, 1) != comb(4, 2) and comb(3, 1) != comb(4, 1) and factorial(4) // factorial(1) != comb(4, 1), "pascal trou 1"
assert comb(4, 1) + comb(4, 2) != comb(5, 1) and comb(4, 1) + comb(4, 2) != comb(4, 3) and comb(2, 0) + comb(2, 1) != comb(4, 1), "pascal trou 2"
assert factorial(3) != 2 ** 3 and 2 ** 3 - 1 != 2 ** 3 and 2 ** 4 != 2 ** 3, "pascal trou 3"
# Exemple : ligne 5
assert [comb(5, j) for j in range(6)] == [1, 5, 10, 10, 5, 1] and comb(4, 1) + comb(4, 2) == 10 == comb(5, 3) and 2 ** 5 == 32
# Contre-exemples et affirmations
assert comb(4, 2) == 6 and factorial(4) // factorial(2) == 12, "pascal contre-exemple n!/k!"
assert comb(4, 1) + comb(4, 2) == comb(5, 2) == 10 and comb(5, 1) == 5, "pascal contre-exemple mauvaise colonne"
assert comb(4, 2) == 6 and 2 * comb(2, 2) == 2, "pascal contre-exemple additivité"
for n in range(1, 40):
    assert comb(2 * n, n) % 2 == 0 and comb(2 * n, n) == 2 * comb(2 * n - 1, n - 1), "pascal affirmation 1 : C(2n,n) pair"
    if n >= 2:
        assert comb(2 * n, 2) == n * (2 * n - 1)
assert comb(4, 2) % 4 != 0, "pascal affirmation 2 : 4 ne divise pas C(4,2)"
for p in sp.primerange(2, 40):
    assert all(comb(p, k) % p == 0 for k in range(1, p)), "pascal affirmation 2 : vrai pour n premier"
assert comb(4, 2) == 6 != 2 * comb(2, 2), "pascal affirmation 3"

# ------------------------------------------------------------------------------------------ ex1 (triangle)
def lignes(n):
    L = [[1]]
    for _ in range(n):
        prec = L[-1]
        L.append([1] + [prec[i] + prec[i + 1] for i in range(len(prec) - 1)] + [1])
    return L

for k, (n, leurres) in enumerate([(6, [comb(6, 3), 6 * 5, perm(6, 4)]), (7, [comb(7, 3), 7 * 6, perm(7, 5)]), (8, [comb(8, 4), 8 * 7, perm(8, 6)])]):
    L = lignes(n)
    assert L[n] == [comb(n, j) for j in range(n + 1)] and L[n][2] == L[n][n - 2]
    assert sorted(rep("proba-t-binomiale-pascal-ex1", k, L[n][n - 2])) == sorted(leurres), f"pascal-ex1 v{k}"
    texte = (EX["proba-t-binomiale-pascal-ex1"] if k == 0 else EX["proba-t-binomiale-pascal-ex1"]["variantes"][k - 1])["correction"]
    for i, ligne in enumerate(L):
        assert f"ligne ${i}$ : " + ", ".join(f"${x}$" for x in ligne) in texte, f"pascal-ex1 v{k} : ligne {i} de la correction"
assert (perm(6, 4), perm(7, 5), perm(8, 6)) == (360, 2520, 20160)

# ------------------------------------------------------------------------------------------ ex2 (poignées de main)
x = sp.symbols("x")
for k, (m, n, leurres) in enumerate([(45, 10, [9, 11, 45]), (66, 12, [11, 13, 66]), (105, 15, [14, 16, 105])]):
    assert sp.solve(sp.Eq(x * (x - 1) / 2, m), x) == [-(n - 1), n] or set(sp.solve(sp.Eq(x * (x - 1) / 2, m), x)) == {-(n - 1), n}
    assert sp.expand((x - n) * (x + n - 1)) == x ** 2 - x - 2 * m, f"pascal-ex2 v{k} : factorisation"
    assert comb(n, 2) == m and len(list(combinations(range(n), 2))) == m
    assert sorted(rep("proba-t-binomiale-pascal-ex2", k, n)) == sorted(leurres), f"pascal-ex2 v{k}"
    assert all(comb(y, 2) != m for y in leurres), f"pascal-ex2 v{k} : leurres faux"

# ------------------------------------------------------------------------------------------ ex3 (chemins)
def nb_chemins(a, b, par=None):
    total = 0
    for pas in product([0, 1], repeat=a + b):  # 1 = pas vers le haut
        if sum(pas) != b:
            continue
        pos, passe = (0, 0), par is None
        for s in pas:
            pos = (pos[0] + (1 - s), pos[1] + s)
            passe = passe or pos == par
        total += passe
    return total

for k, (a, b, c, d, leurres) in enumerate([(5, 3, 2, 1, [56, 13, 10]), (4, 4, 2, 2, [70, 12, 6]), (6, 2, 3, 1, [28, 8, 4])]):
    total, via = nb_chemins(a, b), nb_chemins(a, b, (c, d))
    assert total == comb(a + b, b), f"pascal-ex3 v{k} : nombre total"
    assert via == comb(c + d, d) * comb(a - c + b - d, b - d), f"pascal-ex3 v{k} : par C"
    assert sorted(rep("proba-t-binomiale-pascal-ex3", k, via)) == sorted(leurres), f"pascal-ex3 v{k}"
    assert total in leurres and comb(c + d, d) + comb(a - c + b - d, b - d) in leurres, f"pascal-ex3 v{k} : leurres total et somme"
assert (comb(3, 1), comb(5, 2), comb(4, 2), comb(4, 1)) == (3, 10, 6, 4)
for xx in range(1, 6):
    for yy in range(1, 6):
        assert nb_chemins(xx, yy) == nb_chemins(xx - 1, yy) + nb_chemins(xx, yy - 1), "pascal-ex3 question 3"

# ------------------------------------------------------------------------------------------ ex4 (comité et président)
for n in range(1, 15):
    for k in range(1, n + 1):
        assert k * comb(n, k) == n * comb(n - 1, k - 1), "pascal-ex4 question 1"
    couples = sum(1 for k in range(n + 1) for C in combinations(range(n), k) for p in C) if n <= 10 else None
    s = sum(k * comb(n, k) for k in range(n + 1))
    assert s == n * 2 ** (n - 1), "pascal-ex4 question 3"
    if couples is not None:
        assert couples == s
    assert F(s, 2 ** n) == F(n, 2), "pascal-ex4 question 4 : E(|A|) = n/2"
# Leurres : n 2^n, (n-1) 2^n, n(n+1)/2 ne sont pas égaux à n 2^(n-1) pour tout n >= 1
e = EX["proba-t-binomiale-pascal-ex4"]["reponse_finale"]
assert e["propositions"][e["bonne"]] == "$n2^{n-1}$"
nsym = sp.symbols("n")
for leurre, expr in [("$n2^{n}$", nsym * 2 ** nsym), ("$(n-1)2^{n}$", (nsym - 1) * 2 ** nsym), (r"$\frac{n(n+1)}{2}$", nsym * (nsym + 1) / 2)]:
    assert leurre in e["propositions"] and any(expr.subs(nsym, m) != m * 2 ** (m - 1) for m in range(1, 6)), "pascal-ex4 : leurre faux"

# ==========================================================================================
# Carte proba-t-binomiale-produit
# ==========================================================================================
random.seed(1)
for _ in range(30):
    tailles = [random.randint(1, 4) for _ in range(random.randint(1, 3))]
    lois = []
    for t in tailles:
        poids = [random.randint(1, 9) for _ in range(t)]
        lois.append([F(w, sum(poids)) for w in poids])
    Omega = list(product(*[range(t) for t in tailles]))
    Pw = {w: sp.prod([lois[i][w[i]] for i in range(len(w))]) for w in Omega}
    assert sum(Pw.values()) == 1, "proba-t-binomiale-produit : la probabilité produit a une masse totale 1"
    A = [set(x for x in range(t) if random.random() < 0.5) for t in tailles]
    pave = sum(Pw[w] for w in Omega if all(w[i] in A[i] for i in range(len(w))))
    assert pave == sp.prod([sum(lois[i][x] for x in A[i]) for i in range(len(tailles))]), "produit : formule des pavés"
    for i in range(len(tailles)):
        assert sum(Pw[w] for w in Omega if w[i] in A[i]) == sum(lois[i][x] for x in A[i]), "produit : loi de la i-ème épreuve"
    # affirmation 4 : au moins un succès
    au_moins = sum(Pw[w] for w in Omega if any(w[i] in A[i] for i in range(len(w))))
    assert au_moins == 1 - sp.prod([1 - sum(lois[i][x] for x in A[i]) for i in range(len(tailles))]), "produit affirmation 4"
# affirmation 1 : produit de lois uniformes = loi uniforme
for tailles in [(2, 3), (6, 4, 2)]:
    Omega = list(product(*[range(t) for t in tailles]))
    assert all(sp.prod([F(1, t) for t in tailles]) == F(1, len(Omega)) for _ in Omega), "produit affirmation 1"
# Exemple : urne 2R 1N puis pièce P(pile) = 1/4
ex = {("R", "P"): F(2, 3) * F(1, 4), ("R", "F"): F(2, 3) * F(3, 4), ("N", "P"): F(1, 3) * F(1, 4), ("N", "F"): F(1, 3) * F(3, 4)}
assert ex == {("R", "P"): F(1, 6), ("R", "F"): F(1, 2), ("N", "P"): F(1, 12), ("N", "F"): F(1, 4)} and sum(ex.values()) == 1, "produit exemple"
# Contre-exemple / affirmation 2 : sans remise, 2 rouges et 2 noires
urne = ["R1", "R2", "N1", "N2"]
tir = list(permutations(urne, 2))
pR1 = F(sum(t[0][0] == "R" for t in tir), len(tir))
pR2 = F(sum(t[1][0] == "R" for t in tir), len(tir))
pRR = F(sum(t[0][0] == "R" and t[1][0] == "R" for t in tir), len(tir))
assert pR1 == pR2 == F(1, 2) and pRR == F(1, 6) != pR1 * pR2, "produit contre-exemple sans remise"
# Contre-exemple : uniforme abusive (dé + pièce 1/3) ; somme de deux dés
assert F(6, 12) == F(1, 2) != F(1, 3), "produit contre-exemple uniformité abusive"
deux_des = list(product(range(1, 7), repeat=2))
assert F(sum(a + b == 2 for a, b in deux_des), 36) == F(1, 36) and F(sum(a + b == 7 for a, b in deux_des), 36) == F(1, 6) != F(1, 11), "produit : sommes de deux dés"

# ------------------------------------------------------------------------------------------ ex1 (dé puis pièce, etc.)
def produit2(l1, l2):
    return {(x, y): p * q for x, p in l1.items() for y, q in l2.items()}

de6 = {i: F(1, 6) for i in range(1, 7)}
de4 = {i: F(1, 4) for i in range(1, 5)}
cas = [
    (de6, {"P": F(1, 3), "F": F(2, 3)}, (6, "P"), lambda x, y: x % 2 == 0 and y == "F", F(1, 18), F(1, 3), [F(1, 6), F(1, 4), F(1, 9)]),
    (de4, {"P": F(3, 5), "F": F(2, 5)}, (4, "P"), lambda x, y: x >= 3 and y == "P", F(3, 20), F(3, 10), [F(1, 5), F(1, 4), F(3, 20)]),
    ({"R": F(2, 5), "B": F(3, 5)}, de6, ("R", 6), lambda x, y: x == "B" and y >= 5, F(1, 15), F(1, 5), [F(2, 15), F(1, 6), F(1, 10)]),
]
for k, (l1, l2, issue, ev, p_issue, p_ev, leurres) in enumerate(cas):
    P = produit2(l1, l2)
    assert sum(P.values()) == 1 and P[issue] == p_issue, f"produit-ex1 v{k} : issue"
    juste = sum(p for w, p in P.items() if ev(*w))
    assert juste == p_ev and sorted(rep("proba-t-binomiale-produit-ex1", k, juste)) == sorted(leurres), f"produit-ex1 v{k}"
    # leurre « uniforme sur les issues » (v0 : 3/12 = 1/4 ; v1 : 2/8 = 1/4 ; v2 : 2/12 = 1/6)
    assert F(sum(1 for w in P if ev(*w)), len(P)) in leurres, f"produit-ex1 v{k} : leurre uniforme"
assert 6 * F(1, 18) + 6 * F(1, 9) == 1 and 4 * F(3, 20) + 4 * F(1, 10) == 1 and 6 * F(1, 15) + 6 * F(1, 10) == 1

# ------------------------------------------------------------------------------------------ ex2 (trois machines)
for k, (ps, termes, leurres) in enumerate([
    ((F("0.1"), F("0.2"), F("0.05")), (F("0.076"), F("0.171"), F("0.036")), None),
    ((F("0.2"), F("0.3"), F("0.5")), (F("0.07"), F("0.12"), F("0.28")), None),
    ((F("0.1"), F("0.1"), F("0.2")), (F("0.072"), F("0.072"), F("0.162")), None),
]):
    issues = list(product([1, 0], repeat=3))  # 1 = panne
    proba = {w: sp.prod([ps[i] if w[i] else 1 - ps[i] for i in range(3)]) for w in issues}
    une = sum(p for w, p in proba.items() if sum(w) == 1)
    assert [proba[w] for w in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]] == list(termes) and une == sum(termes), f"produit-ex2 v{k} : termes"
    au_moins = 1 - proba[(0, 0, 0)]
    faux = rep("proba-t-binomiale-produit-ex2", k, une)
    attendus = [au_moins, sum(ps), sp.prod(ps)] if k != 1 else [au_moins, proba[(0, 0, 0)], sp.prod(ps)]
    assert sorted(faux) == sorted(attendus), f"produit-ex2 v{k} : leurres (au moins une, somme ou aucune, produit)"
assert F("0.283") == F("0.076") + F("0.171") + F("0.036")

# ------------------------------------------------------------------------------------------ ex3 (deux dés différents)
for k, (fa, fb, juste, egal, leurres) in enumerate([(6, 4, F(7, 12), F(1, 6), None), (8, 6, F(9, 16), F(1, 8), None), (6, 6, F(5, 12), F(1, 6), None)]):
    issues = list(product(range(1, fa + 1), range(1, fb + 1)))
    q = F(1, len(issues))
    assert q == F(1, fa) * F(1, fb), f"produit-ex3 v{k} : probabilité d'une issue"
    sup = q * sum(a > b for a, b in issues)
    assert sup == juste and q * sum(a == b for a, b in issues) == egal, f"produit-ex3 v{k}"
    sup_eg, inf_eg = q * sum(a >= b for a, b in issues), q * sum(a <= b for a, b in issues)
    faux = rep("proba-t-binomiale-produit-ex3", k, sup)
    assert sorted(faux) == sorted({0: [sup_eg, inf_eg, F(1, 2)], 1: [sup_eg, inf_eg, F(1, 2)], 2: [F(1, 2), sup_eg, egal]}[k]), f"produit-ex3 v{k} : leurres"
assert sum(a > b for a, b in product(range(1, 7), range(1, 5))) == 14 and sum(a > b for a, b in product(range(1, 9), range(1, 7))) == 27

# ------------------------------------------------------------------------------------------ ex4 (match en deux parties gagnantes)
def match_reel(p):
    """Probabilité que A gagne, par simulation exacte des déroulements réels (arrêt à 2 victoires)."""
    total = F(0)
    def jouer(a, b, prob):
        nonlocal total
        if a == 2:
            total += prob
            return
        if b == 2:
            return
        jouer(a + 1, b, prob * p)
        jouer(a, b + 1, prob * (1 - p))
    jouer(0, 0, F(1))
    return total

for k, (p, leurres) in enumerate([(F("0.6"), [F("0.6"), F("0.504"), F("0.36")]), (F("0.7"), [F("0.7"), F("0.637"), F("0.49")]), (F(2, 3), [F(2, 3), F(16, 27), F(4, 9)])]):
    modele = sum(sp.prod([p if x else 1 - p for x in w]) for w in product([1, 0], repeat=3) if sum(w) >= 2)
    assert modele == match_reel(p) == p ** 3 + 3 * p ** 2 * (1 - p) == p ** 2 + 2 * p ** 2 * (1 - p), f"produit-ex4 v{k}"
    # équivalence « A gagne le match » <=> « A gagne au moins 2 des 3 parties » issue par issue
    for w in product([1, 0], repeat=3):
        a = b = 0
        for x in w:
            if a == 2 or b == 2:
                break
            a, b = a + x, b + (1 - x)
        assert (a == 2) == (sum(w) >= 2), "produit-ex4 question 1"
    faux = rep("proba-t-binomiale-produit-ex4", k, modele)
    assert sorted(faux) == sorted(leurres) and p ** 2 in leurres and p ** 2 + p ** 2 * (1 - p) in leurres, f"produit-ex4 v{k} : leurres"
assert (F("0.6") ** 3, 3 * F("0.6") ** 2 * F("0.4"), F("0.6") ** 2, 2 * F("0.6") ** 2 * F("0.4")) == (F("0.216"), F("0.432"), F("0.36"), F("0.288"))
assert (F("0.7") ** 3, 3 * F("0.7") ** 2 * F("0.3"), F("0.7") ** 2, 2 * F("0.7") ** 2 * F("0.3")) == (F("0.343"), F("0.441"), F("0.49"), F("0.294"))
assert (F(2, 3) ** 3, 3 * F(2, 3) ** 2 * F(1, 3), F(2, 3) ** 2, 2 * F(2, 3) ** 2 * F(1, 3)) == (F(8, 27), F(12, 27), F(4, 9), F(8, 27))

print("proba_t_binomiale : toutes les vérifications passent.")
