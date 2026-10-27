# Vérification calculatoire des cartes « proba-l1-conditionnement-composees » et « proba-l1-conditionnement-independance »
# et de leurs exercices (content/cartes/probabilites/l1-proba-l1-conditionnement.json,
# content/exercices/probabilites/l1-proba-l1-conditionnement.json). Calculs exacts en fractions, par énumération d'univers
# finis (tirages de boules étiquetées, permutations, univers produits), indépendamment des formules du cours.
# Chaque assert cite l'identifiant concerné ; « v0 » = exercice tel qu'écrit, « v1 », « v2 » = variantes.
import json
import pathlib
import re
from fractions import Fraction as F
from itertools import permutations, product

import sympy as sp

RACINE = pathlib.Path(__file__).resolve().parent.parent
EX = {e["id"]: e for e in json.loads((RACINE / "content/exercices/probabilites/l1-proba-l1-conditionnement.json").read_text())["exercices"]}
C = "proba-l1-conditionnement-composees"
I = "proba-l1-conditionnement-independance"


def valeur(tex):
    """Texte LaTeX d'un nombre ($\\frac{a}{b}$, $0{,}42$, $-3$) -> valeur exacte SymPy."""
    s = tex.strip().strip("$").replace("{,}", ".").replace(" ", "")
    s = re.sub(r"\\frac\{([^{}]*)\}\{([^{}]*)\}", r"((\1)/(\2))", s)
    s = re.sub(r"\\sqrt\{([^{}]*)\}", r"sqrt(\1)", s)
    s = re.sub(r"(\d+\.\d+)", r"Rational('\1')", s)
    return sp.nsimplify(sp.sympify(s, rational=True))


def version(ident, k):
    e = EX[ident]
    return e if k == 0 else e["variantes"][k - 1]


def verifie_rep(ident, k, attendu):
    """La bonne proposition vaut attendu, les trois autres en diffèrent, les quatre sont distinctes."""
    r = version(ident, k)["reponse_finale"]
    props = r["propositions"]
    assert len(set(props)) == 4, f"{ident} v{k} : propositions distinctes"
    vals = [valeur(x) for x in props]
    assert sp.simplify(vals[r["bonne"]] - sp.nsimplify(attendu)) == 0, f"{ident} v{k} : bonne réponse {props[r['bonne']]} != {attendu}"
    for i, v in enumerate(vals):
        if i != r["bonne"]:
            assert sp.simplify(v - sp.nsimplify(attendu)) != 0, f"{ident} v{k} : le leurre {props[i]} est juste"
    return [v for i, v in enumerate(vals) if i != r["bonne"]]


def verifie_texte(ident, k, attendu_texte):
    r = version(ident, k)["reponse_finale"]
    assert len(set(r["propositions"])) == 4, f"{ident} v{k} : propositions distinctes"
    assert r["propositions"][r["bonne"]] == attendu_texte, f"{ident} v{k} : bonne réponse"


def contient(ident, k, *morceaux):
    t = version(ident, k)["enonce"]
    for m in morceaux:
        assert m in t, f"{ident} v{k} : l'énoncé ne contient pas {m!r}"


def tirages(boules, n):
    """Tirages ordonnés sans remise de n boules étiquetées : liste équiprobable de n-uplets de couleurs."""
    return [tuple(boules[i] for i in t) for t in permutations(range(len(boules)), n)]


def proba(univers, evt):
    return F(sum(1 for w in univers if evt(w)), len(univers))


# ==========================================================================================
# Carte proba-l1-conditionnement-composees
# ==========================================================================================
U = tirages("RRRBB", 3)
pR1R2R3 = proba(U, lambda w: w == ("R", "R", "R"))
assert pR1R2R3 == F(1, 10) == F(3, 5) * F(2, 4) * F(1, 3), "composees exemple : 1/10"
for i in range(3):
    assert proba(U, lambda w: w[i] == "R") == F(3, 5), "composees contre-exemple : P(R_k) = 3/5"
assert F(3, 5) ** 3 == F(27, 125) != pR1R2R3, "composees contre-exemple : produit des P(R_k)"
pR2R3 = proba(U, lambda w: w[1] == "R" and w[2] == "R")
assert pR2R3 == F(3, 10) and pR2R3 / F(3, 5) == F(1, 2), "composees contre-exemple : P(R3 | R2) = 1/2"
assert pR1R2R3 / proba(U, lambda w: w[0] == "R" and w[1] == "R") == F(1, 3), "composees : P(R3 | R1 ∩ R2) = 1/3"
assert F(3, 5) * F(1, 2) * F(1, 2) == F(3, 20) != pR1R2R3, "composees : conditionner par l'étape précédente seulement (trou 2, affirmation 3)"
assert proba(U, lambda w: w[0] == "R" and w[1] == "R") == F(3, 10) != F(9, 25), "composees affirmation 1"
# Hypothèse P(A1 ∩ … ∩ A_{n-1}) > 0 nécessaire ; leurres du trou 1 (P(A1) > 0, P(A_n) > 0, produit > 0) insuffisants
De = range(1, 7)
A1, A2, A3 = {1}, {2}, {1, 2, 3}
assert len(A1) > 0 and len(A2) > 0 and len(A3) > 0 and len(A1 & A2) == 0 and len(A1 & A2 & A3) == 0, "composees contre-exemple : dé"
# Formule générale vérifiée sur des familles aléatoires d'événements d'un univers pondéré
import random
random.seed(65)
for _ in range(300):
    n_om = random.randint(2, 7)
    poids = [random.randint(0, 5) for _ in range(n_om)]
    if sum(poids) == 0:
        continue
    P = [F(p, sum(poids)) for p in poids]
    n = random.randint(2, 4)
    A = [{w for w in range(n_om) if random.random() < 0.6} for _ in range(n)]
    pr = lambda S: sum(P[w] for w in S)
    inter = [set(range(n_om))]
    for a in A:
        inter.append(inter[-1] & a)
    if pr(inter[n - 1]) > 0:
        prod = pr(A[0])
        for k in range(1, n):
            prod *= pr(inter[k] & A[k]) / pr(inter[k])
        assert prod == pr(inter[n]), "composees : formule des probabilités composées (aléatoire)"

# ------------------------------------------------------------------------------------------ ex1 (trois rouges)
for k, (r, b) in enumerate([(4, 6), (5, 7), (6, 4)]):
    contient(C + "-ex1", k, f"${r}$ boules rouges", f"${b}$ boules noires")
    U = tirages("R" * r + "N" * b, 3)
    juste = proba(U, lambda w: w == ("R", "R", "R"))
    faux = verifie_rep(C + "-ex1", k, juste)
    N = r + b
    assert sorted(faux) == sorted([sp.Rational(r, N) ** 3, sp.Rational(r * (r - 1), N * (N - 1)), sp.Rational(1, sp.binomial(N, 3))]), f"composees-ex1 v{k} : leurres"
assert [F(1, 30), F(1, 22), F(1, 6)] == [F(4 * 3 * 2, 720), F(5 * 4 * 3, 1320), F(6 * 5 * 4, 720)], "composees-ex1 : valeurs"

# ------------------------------------------------------------------------------------------ ex2 (clés)
for k, (n, rang) in enumerate([(6, 3), (8, 5), (10, 4)]):
    contient(C + "-ex2", k, f"${n}$ clés")
    assert f"numéro ${rang}$" in version(C + "-ex2", k)["reponse_finale"]["question"], f"composees-ex2 v{k} : rang"
    # suites ordonnées des « rang » premières clés essayées (équiprobables), la clé 0 étant la bonne
    ordres = list(permutations(range(n), rang))
    juste = F(sum(1 for o in ordres if 0 not in o[:-1] and o[-1] == 0), len(ordres))
    assert juste == F(1, n), f"composees-ex2 v{k} : loi uniforme du rang"
    faux = verifie_rep(C + "-ex2", k, juste)
    # leurre « avec remise » : vérifié par l'univers produit des essais avec remise
    avec = F(n - 1, n) ** (rang - 1) * F(1, n)
    assert sp.nsimplify(avec) in faux, f"composees-ex2 v{k} : leurre avec remise"
    del ordres

# ------------------------------------------------------------------------------------------ ex3 (hypergéométrique)
for k, (N, R, n, kk) in enumerate([(10, 4, 3, 2), (12, 5, 4, 2), (8, 3, 3, 1)]):
    contient(C + "-ex3", k, f"${N}$ boules", f"${R}$ rouges", f"${N - R}$ blanches", f"${n}$ boules sans remise")
    U = tirages("R" * R + "B" * (N - R), n)
    juste = proba(U, lambda w: w.count("R") == kk)
    assert juste == F(sp.binomial(R, kk) * sp.binomial(N - R, n - kk), sp.binomial(N, n)), f"composees-ex3 v{k} : hypergéométrique"
    seq = proba(U, lambda w: w == tuple("R" * kk + "B" * (n - kk)))
    # question 2 : toutes les suites à kk rouges ont la même probabilité
    suites = {w for w in U if w.count("R") == kk}
    assert all(proba(U, lambda x, s=s: x == s) == seq for s in suites), f"composees-ex3 v{k} : suites équiprobables"
    faux = verifie_rep(C + "-ex3", k, juste)
    binom = sp.binomial(n, kk) * sp.Rational(R, N) ** kk * sp.Rational(N - R, N) ** (n - kk)
    avec = list(product("R" * R + "B" * (N - R), repeat=n))
    assert proba(avec, lambda w: w.count("R") == kk) == binom, f"composees-ex3 v{k} : cas avec remise"
    assert sorted(faux) == sorted([sp.nsimplify(seq), binom, sp.Rational(sp.binomial(R, kk), sp.binomial(N, n))]), f"composees-ex3 v{k} : leurres"


# ------------------------------------------------------------------------------------------ ex4 (urne de Pólya)
def polya(r, b, c, n):
    """Loi exacte des n premiers tirages de l'urne de Pólya : dictionnaire suite -> probabilité."""
    loi = {(): F(1)}
    for _ in range(n):
        nv = {}
        for s, p in loi.items():
            rr, bb = r + c * s.count("R"), b + c * s.count("B")
            nv[s + ("R",)] = p * F(rr, rr + bb)
            nv[s + ("B",)] = p * F(bb, rr + bb)
        loi = nv
    return loi


for k, (r, b, c) in enumerate([(2, 1, 1), (3, 2, 1), (1, 1, 2)]):
    loi = polya(r, b, c, 3)
    pr = lambda evt: sum(p for s, p in loi.items() if evt(s))
    assert pr(lambda s: s[0] == "R") == pr(lambda s: s[1] == "R") == F(r, r + b), f"composees-ex4 v{k} : P(R1) = P(R2)"
    p12 = pr(lambda s: s[0] == "R" and s[1] == "R") / pr(lambda s: s[1] == "R")
    assert p12 == F(r + c, r + b + c), f"composees-ex4 v{k} : P(R1 | R2)"
    assert pr(lambda s: s == ("R", "R", "B")) == pr(lambda s: s == ("R", "B", "R")) == pr(lambda s: s == ("B", "R", "R")), f"composees-ex4 v{k} : échangeabilité"
    juste = pr(lambda s: s.count("B") == 1)
    faux = verifie_rep(C + "-ex4", k, juste)
    au_moins = pr(lambda s: s.count("B") >= 1)
    assert sorted(faux) == sorted([sp.nsimplify(pr(lambda s: s == ("R", "R", "B"))), 3 * sp.Rational(r, r + b) ** 2 * sp.Rational(b, r + b), sp.nsimplify(au_moins)]), f"composees-ex4 v{k} : leurres"
assert [F(3, 10), F(12, 35), F(3, 16)] == [3 * F(1, 10), 3 * F(4, 35), 3 * F(1, 16)], "composees-ex4 : valeurs"

# ==========================================================================================
# Carte proba-l1-conditionnement-independance
# ==========================================================================================
De = set(range(1, 7))
p6 = lambda S: F(len(S), 6)
A, B = {2, 4, 6}, {1, 2}
assert p6(A & B) == p6(A) * p6(B) == F(1, 6), "independance exemple : A, B indépendants"
Ab, Bb = De - A, De - B
assert Ab & Bb == {3, 5} and p6(Ab & Bb) == p6(Ab) * p6(Bb) == F(1, 3), "independance exemple : contraires"
assert p6(A & Bb) == p6(A) * p6(Bb) and p6(Ab & B) == p6(Ab) * p6(B), "independance : passage aux contraires"
# leurres du trou 2 sur cet exemple : P(A ∩ B̄) = 1/3
assert p6(A & Bb) == F(1, 3) and p6(A) - p6(B) != F(1, 3) and 1 - p6(A) * p6(B) != F(1, 3) and (1 - p6(A)) * p6(B) != F(1, 3), "independance trou 2"
# leurres du trou 1 : (P(A∩B)=0 n'équivaut pas) l'exemple est indépendant avec P(A∩B) ≠ 0 ; P(A∪B) = P(A)+P(B) faux ici ;
# P(A|B) = P(B|A) pour A = B de probabilité 1/2, non indépendants
assert p6(A & B) != 0 and p6(A | B) != p6(A) + p6(B), "independance trou 1"
assert p6(A & A) != p6(A) * p6(A), "independance trou 1 : A = A, P(A|A) = P(A|A) sans indépendance"
# incompatibles
assert p6({1} & {2}) == 0 != F(1, 36), "independance contre-exemple : incompatibles"
# dé truqué
Pt = {1: F(1, 10), 2: F(1, 10), 3: F(1, 10), 4: F(1, 10), 5: F(1, 10), 6: F(1, 2)}
assert sum(Pt.values()) == 1
pt = lambda S: sum(Pt[w] for w in S)
assert pt(A) == F(7, 10) and pt(B) == F(1, 5) and pt(A & B) == F(1, 10) != F(7, 50) == pt(A) * pt(B), "independance : dé truqué (contre-exemple, affirmation 4)"
# non-transitivité et affirmation 3 : A indépendant de lui-même ssi P(A) ∈ {0, 1}
assert F(1, 2) != F(1, 4), "independance : non transitive"
for x in [F(0), F(1, 3), F(1, 2), F(1)]:
    assert (x == x * x) == (x in (0, 1)), "independance affirmation 3"
# affirmation 2 : deux lancers d'une pièce
Om = list(product("PF", repeat=2))
pp = lambda S: F(len(S), 4)
Ae = {w for w in Om if w[0] == w[1]}
Be = {w for w in Om if w[0] == "P"}
Ce = {w for w in Om if w[1] == "P"}
assert pp(Ae & Be) == pp(Ae) * pp(Be) and pp(Ae & Ce) == pp(Ae) * pp(Ce), "independance affirmation 2 : A indépendant de B et de C"
assert pp(Ae & (Be | Ce)) == F(1, 4) != F(3, 8) == pp(Ae) * pp(Be | Ce), "independance affirmation 2 : pas de B ∪ C"
# affirmation 1 et passage aux contraires, sur des univers aléatoires
random.seed(651)
for _ in range(300):
    n_om = random.randint(1, 6)
    poids = [random.randint(0, 4) for _ in range(n_om)]
    if sum(poids) == 0:
        continue
    P = [F(p, sum(poids)) for p in poids]
    pr = lambda S: sum(P[w] for w in S)
    S1 = {w for w in range(n_om) if random.random() < 0.5}
    S2 = {w for w in range(n_om) if random.random() < 0.5}
    tout = set(range(n_om))
    if pr(S1 & S2) == pr(S1) * pr(S2):
        for X1 in (S1, tout - S1):
            for X2 in (S2, tout - S2):
                assert pr(X1 & X2) == pr(X1) * pr(X2), "independance : contraires (aléatoire)"
    if pr(S1) in (0, 1):
        assert pr(S1 & S2) == pr(S1) * pr(S2), "independance : événement de probabilité 0 ou 1"

# ------------------------------------------------------------------------------------------ ex1 (réunion)
for k, (pa, pb) in enumerate([(F(3, 10), F(1, 2)), (F(2, 5), F(1, 4)), (F(1, 3), F(2, 5))]):
    # univers produit de deux épreuves de Bernoulli indépendantes
    loi = {(x, y): (pa if x else 1 - pa) * (pb if y else 1 - pb) for x in (0, 1) for y in (0, 1)}
    juste = sum(p for (x, y), p in loi.items() if x or y)
    assert loi[(0, 0)] == 1 - juste, f"independance-ex1 v{k} : P(Ā ∩ B̄)"
    faux = verifie_rep(I + "-ex1", k, juste)
    assert sorted(faux) == sorted([sp.nsimplify(pa + pb), sp.nsimplify(pa * pb), sp.nsimplify(1 - juste)]), f"independance-ex1 v{k} : leurres"
    assert pa + pb <= 1, f"independance-ex1 v{k} : leurre P(A)+P(B) reste une probabilité"

# ------------------------------------------------------------------------------------------ ex2 (dé)
PROPS = EX[I + "-ex2"]["reponse_finale"]["propositions"]
for k, (f, Aset, Bset) in enumerate([(6, {2, 4, 6}, {1, 2, 3, 4}), (6, {2, 4, 6}, {1, 2, 3}), (8, {4, 8}, {1, 2, 3, 4, 5, 6})]):
    contient(I + "-ex2", k, f"${f}$ faces")
    pr = lambda S: F(len(S), f)
    indep = pr(Aset & Bset) == pr(Aset) * pr(Bset)
    assert Aset & Bset, f"independance-ex2 v{k} : A ∩ B non vide (les leurres « incompatibles » sont faux)"
    verifie_texte(I + "-ex2", k, PROPS[0] if indep else PROPS[1])
    comp = set(range(1, f + 1)) - Aset
    assert (pr(comp & Bset) == pr(comp) * pr(Bset)) == indep, f"independance-ex2 v{k} : question 2"
    assert version(I + "-ex2", k)["reponse_finale"]["propositions"] == PROPS
assert [F(1, 3) == F(1, 2) * F(2, 3), F(1, 6) != F(1, 4), F(1, 8) != F(3, 16)] == [True, True, True], "independance-ex2 : valeurs"

# ------------------------------------------------------------------------------------------ ex3 (parallèle)
for k, (p1, p2) in enumerate([(F(9, 10), F(4, 5)), (F(7, 10), F(3, 5)), (F(1, 2), F(2, 3))]):
    loi = {(x, y): (p1 if x else 1 - p1) * (p2 if y else 1 - p2) for x in (0, 1) for y in (0, 1)}
    marche = sum(p for (x, y), p in loi.items() if x or y)
    un = loi[(1, 0)] + loi[(0, 1)]
    juste = loi[(1, 1)] / marche
    faux = verifie_rep(I + "-ex3", k, juste)
    assert sorted(faux) == sorted([sp.nsimplify(loi[(1, 1)]), sp.nsimplify(marche), sp.nsimplify(un / marche)]), f"independance-ex3 v{k} : leurres"
assert [F(36, 49), F(21, 44), F(2, 5)] == [F(72, 98), F(42, 88), F(1, 3) / F(5, 6)], "independance-ex3 : valeurs"

# ------------------------------------------------------------------------------------------ ex4 (critère ad = bc)
a, b, c, d = sp.symbols("a b c d")
assert sp.expand((a - (a + b) * (a + c) - (a * d - b * c)).subs(d, 1 - a - b - c)) == 0, "independance-ex4 : identité ad - bc"
assert sp.expand((b * c - a * d) + (a * d - b * c)) == 0, "independance-ex4 : couple (A, B̄)"
PROPS = EX[I + "-ex4"]["reponse_finale"]["propositions"]
for k, (aa, bb, cc, dd) in enumerate([(12, 28, 18, 42), (10, 30, 20, 40), (25, 15, 25, 35)]):
    aa, bb, cc, dd = (F(x, 100) for x in (aa, bb, cc, dd))
    assert aa + bb + cc + dd == 1
    diff = aa - (aa + bb) * (aa + cc)
    assert diff == aa * dd - bb * cc, f"independance-ex4 v{k}"
    attendu = PROPS[0] if diff == 0 else (PROPS[1] if diff > 0 else PROPS[2])
    verifie_texte(I + "-ex4", k, attendu)
print("proba_l1_conditionnement : toutes les vérifications passent.")
