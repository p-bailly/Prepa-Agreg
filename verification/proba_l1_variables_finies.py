# Vérification calculatoire des cartes du chapitre « proba-l1-variables-finies » (transfert, variance, uniforme-bernoulli,
# couples, independance, covariance) et de leurs exercices (content/cartes/probabilites/l1-proba-l1-variables-finies.json,
# content/exercices/probabilites/l1-proba-l1-variables-finies.json). Calculs exacts (fractions, SymPy), par énumération
# d'univers finis, indépendamment des formules du cours. Chaque assert cite l'identifiant concerné ; « v0 » = exercice tel
# qu'écrit, « v1 », « v2 » = variantes.
import json
import pathlib
import random
import re
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial

import sympy as sp

RACINE = pathlib.Path(__file__).resolve().parent.parent
EX = {e["id"]: e for e in json.loads((RACINE / "content/exercices/probabilites/l1-proba-l1-variables-finies.json").read_text())["exercices"]}
P_ = "proba-l1-variables-finies-"
TR, VA, UB, CO, IN, CV = (P_ + s for s in ["transfert", "variance", "uniforme-bernoulli", "couples", "independance", "covariance"])


def valeur(tex):
    s = tex.strip().strip("$").replace("{,}", ".").replace(" ", "")
    s = re.sub(r"\\sqrt\{([^{}]*)\}", r"sqrt(\1)", s)
    s = re.sub(r"\\frac\{([^{}]*)\}\{([^{}]*)\}", r"((\1)/(\2))", s)
    s = re.sub(r"(\d)sqrt", r"\1*sqrt", s)
    s = re.sub(r"(\d+\.\d+)", r"Rational('\1')", s)
    return sp.nsimplify(sp.sympify(s, rational=True))


def version(ident, k):
    e = EX[ident]
    return e if k == 0 else e["variantes"][k - 1]


def verifie_rep(ident, k, attendu):
    r = version(ident, k)["reponse_finale"]
    props = r["propositions"]
    assert len(set(props)) == 4, f"{ident} v{k} : propositions distinctes"
    vals = [valeur(x) for x in props]
    att = sp.nsimplify(attendu)
    assert sp.simplify(vals[r["bonne"]] - att) == 0, f"{ident} v{k} : bonne réponse {props[r['bonne']]} != {attendu}"
    for i, v in enumerate(vals):
        if i != r["bonne"]:
            assert sp.simplify(v - att) != 0, f"{ident} v{k} : le leurre {props[i]} est juste"
    return sorted([v for i, v in enumerate(vals) if i != r["bonne"]], key=lambda x: float(x))


def memes(a, b):
    return sorted(sp.nsimplify(x) for x in a) == sorted(sp.nsimplify(x) for x in b)


def verifie_texte(ident, k, attendu_texte):
    r = version(ident, k)["reponse_finale"]
    assert len(set(r["propositions"])) == 4, f"{ident} v{k} : propositions distinctes"
    assert r["propositions"][r["bonne"]] == attendu_texte, f"{ident} v{k} : bonne réponse"


def contient(ident, k, *morceaux):
    t = version(ident, k)["enonce"]
    for m in morceaux:
        assert m in t, f"{ident} v{k} : l'énoncé ne contient pas {m!r}"


def E(loi, f=lambda x: x):
    """loi : dict valeur -> probabilité (ou issue -> probabilité avec f)."""
    return sum(f(x) * p for x, p in loi.items())


def V(loi, f=lambda x: x):
    m = E(loi, f)
    return E(loi, lambda x: (f(x) - m) ** 2)


def loi_de(univers, f):
    d = {}
    for w, p in univers.items():
        d[f(w)] = d.get(f(w), 0) + p
    return d


def unif(liste):
    return {w: F(1, len(liste)) for w in liste}


# ==========================================================================================
# Cartes
# ==========================================================================================
# transfert : exemple, contre-exemples, affirmations
Lx = {-1: F(1, 4), 0: F(1, 2), 2: F(1, 4)}
assert sum(Lx.values()) == 1 and E(Lx, lambda x: x * x) == F(5, 4) and E(Lx) == F(1, 4) and E(Lx) ** 2 == F(1, 16), "transfert exemple"
assert E(loi_de(Lx, lambda x: x * x)) == E(Lx, lambda x: x * x), "transfert : loi de X^2 et transfert"
piece = unif(["P", "F"])
Xp = {"P": 1, "F": 0}
Yp = {w: 1 - Xp[w] for w in piece}
assert loi_de(piece, Xp.get) == loi_de(piece, Yp.get) and all(Xp[w] != Yp[w] for w in piece), "transfert : même loi, variables différentes"
assert E(piece, lambda w: Xp[w] * Xp[w]) == F(1, 2) and E(piece, lambda w: Xp[w] * Yp[w]) == 0, "transfert : E(XX) != E(XY)"
assert loi_de(piece, lambda w: Xp[w] + Yp[w]) == {1: 1} and loi_de(piece, lambda w: 2 * Xp[w]) == {0: F(1, 2), 2: F(1, 2)}, "transfert affirmation 4"
assert E(unif([-1, 1]), lambda x: x * x) == 1 and E(unif([-1, 1])) ** 2 == 0, "transfert affirmation 3"
# leurres des trous du transfert, faux sur l'exemple avec f(x) = x^2 (valeur 5/4)
assert E(Lx) ** 2 != F(5, 4), "transfert trou 1 : f(E(X))"
assert sum(x * x * Lx.get(x * x, 0) for x in Lx) != F(5, 4), "transfert trou 1 : f(x)P(X=f(x))"
assert F(sum(x * x for x in Lx), len(Lx)) != F(5, 4), "transfert trou 1 : moyenne non pondérée"
# trou 2 : même espérance, même ensemble de valeurs, mais E(X^2) différents
L1, L2 = {-1: F(1, 2), 1: F(1, 2)}, {-1: F(1, 4), 0: F(1, 2), 1: F(1, 4)}
L3 = {-1: F(1, 3), 0: F(1, 3), 1: F(1, 3)}
assert E(L1) == E(L2) and E(L1, lambda x: x * x) != E(L2, lambda x: x * x), "transfert trou 2 : l'espérance ne suffit pas"
assert set(L2) == set(L3) and E(L2, lambda x: x * x) != E(L3, lambda x: x * x), "transfert trou 2 : les valeurs ne suffisent pas"

# variance
Lv = {0: F(1, 2), 1: F(1, 4), 4: F(1, 4)}
assert E(Lv) == F(5, 4) and E(Lv, lambda x: x * x) == F(17, 4) and V(Lv) == F(43, 16), "variance exemple"
assert V(unif([-1, 1])) == 1, "variance contre-exemple : E(X) = 0, V(X) = 1"
Om2 = {"a": F(1), "b": F(0)}
Xz = {"a": 0, "b": 1}
assert E(Om2, Xz.get) == 0 and V(Om2, Xz.get) == 0 and len(set(Xz.values())) == 2, "variance : variance nulle, X non constante (contre-exemple, affirmation 2, trou 2)"
U12 = unif([1, 2])
assert E(U12, lambda x: x * x) - E(U12) == 1 and V(U12) == F(1, 4), "variance contre-exemple : E(X^2) - E(X)"
# leurres du trou 1 sur U12 : V = 1/4
assert E(U12) ** 2 - E(U12, lambda x: x * x) != F(1, 4) and E(U12, lambda x: x * x) - E(U12) != F(1, 4) and E(U12, lambda x: x * x) - 2 * E(U12) ** 2 != F(1, 4), "variance trou 1"
# leurres du trou 2 : E(X) = 0 sans variance nulle (±1) ; X constante 1 : V = 0 mais E(X^2) = 1
assert V({1: F(1)}) == 0 and E({1: F(1)}, lambda x: x * x) == 1, "variance trou 2 : E(X^2) = 0"
assert V({10: F(1)}) == 0 <= 1 == V(unif([-1, 1])) and E({10: F(1)}, lambda x: x * x) == 100 > 1, "variance affirmation 4"
random.seed(652)
for _ in range(300):
    vals = random.sample(range(-6, 7), random.randint(1, 4))
    poids = [random.randint(1, 5) for _ in vals]
    L = {v: F(p, sum(poids)) for v, p in zip(vals, poids)}
    m = E(L)
    assert V(L) == E(L, lambda x: x * x) - m * m >= 0, "variance : König-Huygens"
    a = F(random.randint(-8, 8), random.randint(1, 3))
    assert E(L, lambda x: (x - a) ** 2) == V(L) + (m - a) ** 2, "variance : forme générale"
    assert (V(L) == 0) == (L.get(m, 0) == 1), "variance : cas de variance nulle"
    lo, hi = min(vals), max(vals)
    assert V(L) <= (hi - m) * (m - lo) <= F((hi - lo) ** 2, 4), "variance-ex4 : V <= (M-mu)(mu-m) <= (M-m)^2/4"
    assert (V(L) == (hi - m) * (m - lo)) == (len(vals) <= 2), "variance-ex4 : cas d'égalité"

# uniforme et Bernoulli
for n in range(1, 30):
    L = unif(list(range(1, n + 1)))
    assert E(L) == F(n + 1, 2) and V(L) == F(n * n - 1, 12) and E(L, lambda x: x * x) == F((n + 1) * (2 * n + 1), 6), "uniforme-bernoulli : formules"
    if n >= 1:
        assert F(n, 2) != E(L) and F(n - 1, 2) != E(L), "uniforme-bernoulli trou 1"
    if n >= 2:
        assert F(n * (n + 1), 2) != E(L), "uniforme-bernoulli trou 1 : n(n+1)/2 (faux dès n = 2)"
    assert F(n * n + 1, 12) != V(L) and F((n + 1) ** 2 - 1, 12) != V(L) and F((n + 1) * (2 * n + 1), 6) != V(L), "uniforme-bernoulli trou 2"
    for a in range(-5, 6):
        La = unif(list(range(a, a + n)))
        assert E(La) == F(2 * a + n - 1, 2) and V(La) == F(n * n - 1, 12), "uniforme-bernoulli : loi uniforme sur {a,…,b}"
assert V(unif([0, 1])) == F(1, 4) and F(1 * 1 - 1, 12) == 0, "uniforme-bernoulli contre-exemple : {0,1}"
assert V(unif([0, 2])) == 1, "uniforme-bernoulli contre-exemple : valeurs 0 et 2"
Oabc = {"a": F(1, 2), "b": F(1, 4), "c": F(1, 4)}
assert loi_de(Oabc, {"a": 0, "b": 1, "c": 1}.get) == {0: F(1, 2), 1: F(1, 2)}, "uniforme-bernoulli contre-exemple : X uniforme, Ω non équiprobable"
deux_des = unif(list(product(range(1, 7), repeat=2)))
LS = loi_de(deux_des, sum)
assert LS[2] == F(1, 36) and LS[7] == F(1, 6), "uniforme-bernoulli : somme de deux dés non uniforme (affirmation 1)"
for p in [F(0), F(1, 7), F(1, 2), F(3, 4), F(1)]:
    B = {1: p, 0: 1 - p}
    assert E(B) == p and V(B) == p * (1 - p) <= F(1, 4), "uniforme-bernoulli : Bernoulli"
    assert loi_de(B, lambda x: x * x) == loi_de(B, lambda x: x), "uniforme-bernoulli affirmation 2"
    if 0 < p < 1 and p != F(1, 2):
        assert V(B) != p and V(B) != p * p and sp.sqrt(p * (1 - p)) != V(B), "uniforme-bernoulli trou 3"
assert loi_de(unif([1, 2]), lambda x: 2 * x).get(3, 0) == 0, "uniforme-bernoulli affirmation 4"
assert F(1, 6) * F(5, 6) == F(5, 36), "uniforme-bernoulli exemple"
assert V(unif(list(range(1, 7)))) == F(35, 12), "uniforme-bernoulli exemple : dé"

# couples
urne3 = unif(list(permutations([1, 2, 3], 2)))
Lc = loi_de(urne3, lambda w: w)
assert all(Lc.get((i, j), 0) == (F(1, 6) if i != j else 0) for i in (1, 2, 3) for j in (1, 2, 3)), "couples exemple : loi conjointe"
assert loi_de(urne3, lambda w: w[0]) == unif([1, 2, 3]) == loi_de(urne3, lambda w: w[1]), "couples exemple : marginales uniformes"
assert E(urne3, lambda w: w[0] * w[1]) == F(11, 3), "couples exemple : E(XY) = 11/3"
assert Lc.get((1, 1), 0) == 0 != F(1, 9), "couples contre-exemple : pas le produit des marginales"
assert E(piece, lambda w: Xp[w] * Xp[w]) == F(1, 2) != 0 == E(piece, lambda w: Xp[w] * Yp[w]), "couples contre-exemple : mêmes marginales"
# affirmation 3 (fausse) : Y constante 0, X de Bernoulli 1/2 ; trou 1 : leurres faux sur l'urne
assert F(1, 2) != 1, "couples affirmation 3"
x0 = 1
assert sum(Lc.get((x0, y), 0) for y in (1, 2, 3)) == F(1, 3), "couples : marginale"
assert sum(Lc.get((x0, y), 0) / F(1, 3) for y in (1, 2, 3)) != F(1, 3), "couples trou 1 : somme des conditionnelles"
assert sp.prod([Lc.get((x0, y), 0) for y in (1, 2, 3)]) != F(1, 3) and max(Lc.get((x0, y), 0) for y in (1, 2, 3)) != F(1, 3), "couples trou 1 : produit, max"

# indépendance de variables
D6 = unif(list(range(1, 7)))
Xd = lambda w: w % 2 == 0
Yd = lambda w: (w - 1) // 2
assert all(sum(p for w, p in D6.items() if Xd(w) == x and Yd(w) == y) == F(1, 6) for x in (False, True) for y in (0, 1, 2)), "independance exemple"
U3 = unif([-1, 0, 1])
assert E(U3, lambda x: x * x * x) == 0 == E(U3) * E(U3, lambda x: x * x), "independance contre-exemple : E(XY) = E(X)E(Y)"
assert F(1, 3) * F(2, 3) == F(2, 9) != 0, "independance contre-exemple : P(X=0,Y=1) = 0"
pieces = unif(list(product([0, 1], repeat=2)))
Zx = lambda w: 1 if w[0] == w[1] else 0
fam = [lambda w: w[0], lambda w: w[1], Zx]


def indep(univ, f, g):
    lf, lg = loi_de(univ, f), loi_de(univ, g)
    lfg = loi_de(univ, lambda w: (f(w), g(w)))
    return all(lfg.get((a, b), 0) == lf[a] * lg[b] for a in lf for b in lg)


assert all(indep(pieces, fam[i], fam[j]) for i in range(3) for j in range(i + 1, 3)), "independance : deux à deux (contre-exemple, affirmation 2)"
assert sum(p for w, p in pieces.items() if w[0] == 1 and w[1] == 1 and Zx(w) == 0) == 0 != F(1, 8), "independance : pas mutuellement"
assert not indep(pieces, lambda w: w[0], lambda w: w[0] + w[1]), "independance affirmation 4 : X et X+Y"
random.seed(653)
for _ in range(100):
    # trois variables mutuellement indépendantes sur un univers produit : sous-familles et fonctions (affirmations 1 et 3)
    lois = [[random.randint(1, 4) for _ in range(random.randint(1, 3))] for _ in range(3)]
    univ = {}
    for w in product(*[range(len(l)) for l in lois]):
        univ[w] = sp.prod([F(lois[i][w[i]], sum(lois[i])) for i in range(3)])
    assert indep(univ, lambda w: w[0], lambda w: w[1]), "independance : sous-famille"
    assert indep(univ, lambda w: w[0] ** 2, lambda w: 3 - w[1]), "independance : fonctions"
assert indep(D6, lambda w: 5, lambda w: w), "independance : variable constante"

# covariance
assert E(urne3, lambda w: w[0] * w[1]) - E(urne3, lambda w: w[0]) * E(urne3, lambda w: w[1]) == F(-1, 3), "covariance exemple"
assert V(urne3, lambda w: w[0]) == F(2, 3) and V(urne3, lambda w: w[0] + w[1]) == F(2, 3), "covariance exemple : V(X+Y)"
assert loi_de(urne3, lambda w: w[0] + w[1]) == unif([3, 4, 5]), "covariance exemple : loi de X+Y"
quatre = unif([(1, 0), (-1, 0), (0, 1), (0, -1)])
cov = lambda u, f, g: E(u, lambda w: f(w) * g(w)) - E(u, f) * E(u, g)
fx, fy = (lambda w: w[0]), (lambda w: w[1])
assert cov(quatre, fx, fy) == 0 and loi_de(quatre, fx) == loi_de(quatre, fy) and E(quatre, fx) == 0, "covariance contre-exemple : covariance nulle"
assert not indep(quatre, fx, fy), "covariance contre-exemple : non indépendantes (trou 3)"
assert cov(U3, lambda x: x, lambda x: x * x) == 0 and not indep(U3, lambda x: x, lambda x: x * x), "covariance affirmation 1"
U11 = unif([-1, 1])
assert cov(U11, lambda x: x, lambda x: 2 * x) == 2 > 1 == V(U11), "covariance affirmation 4"
for _ in range(300):
    n_om = random.randint(1, 5)
    poids = [random.randint(1, 4) for _ in range(n_om)]
    u = {w: F(poids[w], sum(poids)) for w in range(n_om)}
    X = {w: random.randint(-4, 4) for w in u}
    Y = {w: random.randint(-4, 4) for w in u}
    Z = {w: random.randint(-4, 4) for w in u}
    a, b = random.randint(-3, 3), random.randint(-3, 3)
    cxy = cov(u, X.get, Y.get)
    assert cxy == E(u, lambda w: (X[w] - E(u, X.get)) * (Y[w] - E(u, Y.get))), "covariance : forme développée"
    assert cov(u, lambda w: a * X[w] + b * Z[w], Y.get) == a * cxy + b * cov(u, Z.get, Y.get), "covariance : bilinéarité"
    assert cov(u, lambda w: a * X[w] + b, Y.get) == a * cxy, "covariance affirmation 2"
    assert V(u, lambda w: X[w] + Y[w]) == V(u, X.get) + V(u, Y.get) + 2 * cxy, "covariance : V(X+Y)"
    assert V(u, lambda w: X[w] - Y[w]) == V(u, X.get) + V(u, Y.get) - 2 * cxy, "covariance affirmation 3"
    assert cxy ** 2 <= V(u, X.get) * V(u, Y.get), "covariance : Cauchy-Schwarz"
# leurres des trous de la covariance, faux sur l'exemple de l'urne (Cov = -1/3, V(X+Y) = 2/3)
mX = E(urne3, fx)
assert mX * mX - F(11, 3) != F(-1, 3) and 0 != F(-1, 3), "covariance trou 1"
# leurre E(XY) - E(X) - E(Y) (qui coïncide par hasard sur l'urne) : faux pour X = Y = 1 constantes (Cov = 0, leurre = -1)
assert cov({0: F(1)}, lambda w: 1, lambda w: 1) == 0 != 1 - 1 - 1, "covariance trou 1 : E(XY)-E(X)-E(Y)"
assert F(4, 3) != F(2, 3) and F(4, 3) - F(1, 3) != F(2, 3) and F(4, 3) + F(2, 3) != F(2, 3), "covariance trou 2"

# ==========================================================================================
# Exercices
# ==========================================================================================
# ------------------------------------------------------------------------------------------ transfert-ex1
for k, (vals, probs) in enumerate([([-2, 0, 1, 3], [F(1, 10), F(3, 10), F(4, 10), F(2, 10)]),
                                   ([-1, 0, 2, 4], [F(2, 10), F(4, 10), F(3, 10), F(1, 10)]),
                                   ([-3, -1, 1, 2], [F(1, 6), F(1, 3), F(1, 3), F(1, 6)])]):
    assert sum(probs) == 1
    contient(TR + "-ex1", k, " & ".join(map(str, vals)))
    L = dict(zip(vals, probs))
    faux = verifie_rep(TR + "-ex1", k, E(L, lambda x: x * x))
    assert memes(faux, [E(L) ** 2, E(L), F(sum(v * v for v in vals), len(vals))]), f"transfert-ex1 v{k} : leurres"

# ------------------------------------------------------------------------------------------ transfert-ex2
for k, (n, m) in enumerate([(6, 3), (6, 2), (4, 2)]):
    contient(TR + "-ex2", k, f"${n}$ faces", f"$Y=(X-{m})^2$")
    L = unif(list(range(1, n + 1)))
    EY = E(L, lambda x: (x - m) ** 2)
    assert E(loi_de(L, lambda x: (x - m) ** 2)) == EY, f"transfert-ex2 v{k}"
    assert EY - (E(L) - m) ** 2 == V(L), f"transfert-ex2 v{k} : l'écart est la variance de X"
    faux = verifie_rep(TR + "-ex2", k, EY)
    valsY = set((x - m) ** 2 for x in L)
    assert memes(faux, [(E(L) - m) ** 2, F(sum(valsY), len(valsY)), E(L, lambda x: x * x) - m * m]), f"transfert-ex2 v{k} : leurres"

# ------------------------------------------------------------------------------------------ transfert-ex3
for k, (f, nd) in enumerate([(6, 2), (4, 2), (6, 3)]):
    contient(TR + "-ex3", k, f"${f}$ faces", {2: "deux dés", 3: "trois dés"}[nd])
    U = unif(list(product(range(1, f + 1), repeat=nd)))
    EM, Emin = E(U, max), E(U, min)
    faux = verifie_rep(TR + "-ex3", k, EM)
    assert memes(faux, [Emin, F(f + 1, 2), F(f * nd, nd + 1)]), f"transfert-ex3 v{k} : leurres"
    LM = loi_de(U, max)
    assert all(LM[j] == F(j ** nd - (j - 1) ** nd, f ** nd) for j in range(1, f + 1)), f"transfert-ex3 v{k} : loi de M"
assert [F(161, 36), F(25, 8), F(119, 24)] == [F(161, 36), F(50, 16), F(1071, 216)], "transfert-ex3 : valeurs"

# ------------------------------------------------------------------------------------------ transfert-ex4
for k, (vals, probs, med) in enumerate([([0, 1, 2, 10], [F(3, 10), F(3, 10), F(2, 10), F(2, 10)], 1),
                                        ([0, 2, 3, 9], [F(2, 10), F(4, 10), F(1, 10), F(3, 10)], 2),
                                        ([-4, 0, 1, 2], [F(2, 10), F(1, 10), F(4, 10), F(3, 10)], 1)]):
    L = dict(zip(vals, probs))
    g = lambda a: E(L, lambda x: abs(x - a))
    grille = [F(i, 20) for i in range(20 * (min(vals) - 2), 20 * (max(vals) + 2) + 1)]
    gmin = min(g(a) for a in grille)
    assert g(F(med)) == gmin and all(g(a) > gmin for a in grille if a != med), f"transfert-ex4 v{k} : minimum unique en la médiane"
    faux = verifie_rep(TR + "-ex4", k, gmin)
    autre = vals[vals.index(med) + 1]
    assert memes(faux, [g(E(L)), E(L), g(F(autre))]), f"transfert-ex4 v{k} : leurres"
    # pentes annoncées : 2P(X < a) - 1 entre deux valeurs consécutives
    for x1, x2 in zip(vals, vals[1:]):
        a1, a2 = x1 + F(1, 3) * (x2 - x1), x1 + F(2, 3) * (x2 - x1)
        assert (g(a2) - g(a1)) / (a2 - a1) == 2 * sum(p for x, p in L.items() if x <= x1) - 1, f"transfert-ex4 v{k} : pentes"

# ------------------------------------------------------------------------------------------ variance-ex1
for k, (vals, probs) in enumerate([([0, 1, 2, 3], [F(1, 10), F(4, 10), F(3, 10), F(2, 10)]),
                                   ([-1, 0, 1, 2], [F(2, 10), F(3, 10), F(3, 10), F(2, 10)]),
                                   ([1, 2, 5], [F(1, 2), F(1, 4), F(1, 4)])]):
    assert sum(probs) == 1
    contient(VA + "-ex1", k, " & ".join(map(str, vals)))
    L = dict(zip(vals, probs))
    faux = verifie_rep(VA + "-ex1", k, V(L))
    assert memes(faux, [E(L, lambda x: x * x) - E(L), E(L, lambda x: x * x), E(L) ** 2]), f"variance-ex1 v{k} : leurres"

# ------------------------------------------------------------------------------------------ variance-ex2
for k, (gain, perte, p) in enumerate([(10, 2, F(1, 5)), (8, 1, F(1, 4)), (4, 2, F(1, 2))]):
    contient(VA + "-ex2", k, f"gagne ${gain}$ euros", f"perd ${perte}$ euros")
    L = {gain: p, -perte: 1 - p}
    var = V(L)
    assert var == p * (1 - p) * (gain + perte) ** 2, f"variance-ex2 v{k} : formule"
    faux = verifie_rep(VA + "-ex2", k, sp.sqrt(sp.nsimplify(var)))
    assert memes(faux, [var, p * (1 - p) * (gain + perte), sp.sqrt(sp.nsimplify(E(L, lambda x: x * x)))]), f"variance-ex2 v{k} : leurres"
a_, b_, p_ = sp.symbols("a b p")
assert sp.expand(p_ * a_ ** 2 + (1 - p_) * b_ ** 2 - (p_ * a_ + (1 - p_) * b_) ** 2 - p_ * (1 - p_) * (a_ - b_) ** 2) == 0, "variance-ex2 : identité"

# ------------------------------------------------------------------------------------------ variance-ex3
for k, (vals, probs, a) in enumerate([([0, 1, 4], [F(2, 10), F(5, 10), F(3, 10)], 3),
                                      ([-1, 0, 2], [F(1, 4), F(1, 4), F(1, 2)], 2),
                                      ([1, 2, 6], [F(1, 3), F(1, 3), F(1, 3)], 1)]):
    contient(VA + "-ex3", k, " & ".join(map(str, vals)), f"E\\left((X-{a})^2\\right)")
    L = dict(zip(vals, probs))
    juste = E(L, lambda x: (x - a) ** 2)
    faux = verifie_rep(VA + "-ex3", k, juste)
    assert memes(faux, [V(L), (E(L) - a) ** 2, V(L) - (E(L) - a) ** 2]), f"variance-ex3 v{k} : leurres"

# ------------------------------------------------------------------------------------------ uniforme-bernoulli-ex1
for k, n in enumerate([10, 20, 8]):
    contient(UB + "-ex1", k, f"\\{{1,\\dots,{n}\\}}")
    L = unif(list(range(1, n + 1)))
    faux = verifie_rep(UB + "-ex1", k, V(L))
    assert memes(faux, [F(n * n + 1, 12), E(L, lambda x: x * x), F((n + 1) ** 2 - 1, 12)]), f"uniforme-bernoulli-ex1 v{k} : leurres"

# ------------------------------------------------------------------------------------------ uniforme-bernoulli-ex2
for k, (a, b) in enumerate([(5, 14), (-3, 3), (2, 6)]):
    contient(UB + "-ex2", k, f"\\{{{a},\\dots,{b}\\}}")
    L = unif(list(range(a, b + 1)))
    assert E(L) == F(a + b, 2), f"uniforme-bernoulli-ex2 v{k} : E(X)"
    faux = verifie_rep(UB + "-ex2", k, V(L))
    troisieme = E(L, lambda x: x * x) if E(L) != 0 else F((b - a) ** 2, 12)
    assert memes(faux, [F((b - a) ** 2 - 1, 12), F(b * b - 1, 12), troisieme]), f"uniforme-bernoulli-ex2 v{k} : leurres"

# ------------------------------------------------------------------------------------------ uniforme-bernoulli-ex3
for k, (pa, pb, pab) in enumerate([(F(1, 2), F(2, 5), F(1, 10)), (F(3, 5), F(1, 2), F(2, 5)), (F(1, 2), F(1, 3), F(0))]):
    # univers à quatre issues : A∩B, A∩B̄, Ā∩B, Ā∩B̄
    u = {(1, 1): pab, (1, 0): pa - pab, (0, 1): pb - pab, (0, 0): 1 - pa - pb + pab}
    assert all(p >= 0 for p in u.values()) and sum(u.values()) == 1
    LX = loi_de(u, sum)
    juste = V(LX)
    faux = verifie_rep(UB + "-ex3", k, juste)
    sv = V(u, lambda w: w[0]) + V(u, lambda w: w[1])
    autres = [sv, E(LX)] + ([E(LX, lambda x: x * x)] if pab else [LX[0]])
    assert memes(faux, autres), f"uniforme-bernoulli-ex3 v{k} : leurres"
    assert juste - sv == 2 * (pab - pa * pb) != 0, f"uniforme-bernoulli-ex3 v{k} : écart non nul"

# ------------------------------------------------------------------------------------------ uniforme-bernoulli-ex4
assert loi_de(deux_des, sum)[2] == F(1, 36) != F(1, 11), "uniforme-bernoulli-ex4 : question 1"
for _ in range(200):
    pp = [random.randint(0, 5) for _ in range(6)]
    qq = [random.randint(0, 5) for _ in range(6)]
    if sum(pp) == 0 or sum(qq) == 0:
        continue
    pp = [F(x, sum(pp)) for x in pp]
    qq = [F(x, sum(qq)) for x in qq]
    LS = {}
    for i in range(6):
        for j in range(6):
            LS[i + j + 2] = LS.get(i + j + 2, 0) + pp[i] * qq[j]
    assert LS.get(7, 0) >= pp[0] * qq[5] + pp[5] * qq[0], "uniforme-bernoulli-ex4 : question 3"
    x, y = pp[0] * qq[5], pp[5] * qq[0]
    assert (x + y) ** 2 >= 4 * x * y, "uniforme-bernoulli-ex4 : x + y >= 2 sqrt(xy)"
    assert not all(LS.get(s, 0) == F(1, 11) for s in range(2, 13)), "uniforme-bernoulli-ex4 : jamais uniforme"
assert F(2, 11) > F(1, 11) and F(1, 11) ** 2 == F(1, 121), "uniforme-bernoulli-ex4 : conclusion"

# ------------------------------------------------------------------------------------------ couples-ex1
c = lambda x: F(x, 100)
for k, tab in enumerate([[[c(10), c(20)], [c(20), c(10)], [c(30), c(10)]],
                         [[c(5), c(15)], [c(25), c(15)], [c(10), c(30)]],
                         [[c(15), c(5)], [c(10), c(20)], [c(30), c(20)]]]):
    u = {(x, y): tab[x][y] for x in range(3) for y in range(2)}
    assert sum(u.values()) == 1
    juste = E(u, lambda w: w[0])
    faux = verifie_rep(CO + "-ex1", k, juste)
    assert memes(faux, [1, sum(x * tab[x][0] for x in range(3)), E(u, lambda w: w[1])]), f"couples-ex1 v{k} : leurres"

# ------------------------------------------------------------------------------------------ couples-ex2
for k, f in enumerate([4, 3, 6]):
    contient(CO + "-ex2", k, f"${f}$ faces")
    U = unif(list(product(range(1, f + 1), repeat=2)))
    LXY = loi_de(U, lambda w: (min(w), max(w)))
    assert all(LXY.get((i, j), 0) == (F(2, f * f) if i < j else (F(1, f * f) if i == j else 0)) for i in range(1, f + 1) for j in range(1, f + 1)), f"couples-ex2 v{k} : loi conjointe"
    assert all(loi_de(U, min)[i] == F(2 * f + 1 - 2 * i, f * f) for i in range(1, f + 1)), f"couples-ex2 v{k} : marginale de X"
    juste = E(U, lambda w: max(w) - min(w))
    faux = verifie_rep(CO + "-ex2", k, juste)
    assert memes(faux, [juste / 2, 0, 2 * juste]), f"couples-ex2 v{k} : leurres"

# ------------------------------------------------------------------------------------------ couples-ex3
for k, (R, N, B) in enumerate([(2, 3, 1), (3, 3, 2), (2, 2, 2)]):
    contient(CO + "-ex3", k, f"${R}$ boules rouges", f"${N}$ noires", f"${B}$ blanche")
    boules = "R" * R + "N" * N + "B" * B
    U = unif(list(combinations(range(len(boules)), 2)))
    Xc = lambda w: sum(boules[i] == "R" for i in w)
    Yc = lambda w: sum(boules[i] == "N" for i in w)
    EXY = E(U, lambda w: Xc(w) * Yc(w))
    faux = verifie_rep(CO + "-ex3", k, EXY)
    assert memes(faux, [E(U, Xc) * E(U, Yc), loi_de(U, Xc)[1] * loi_de(U, Yc)[1], 2 * EXY]), f"couples-ex3 v{k} : leurres"
    assert loi_de(U, lambda w: Xc(w) + Yc(w)) == loi_de(U, lambda w: 2 - sum(boules[i] == "B" for i in w)), f"couples-ex3 v{k} : X+Y"


# ------------------------------------------------------------------------------------------ couples-ex4 (bornes de Fréchet)
def intervalle(tex):
    s = tex.strip().strip("$")
    m = re.fullmatch(r"\\left\[(.*),(.*)\\right\]", s)
    if m:
        return ("I", valeur(m.group(1)), valeur(m.group(2)))
    m = re.fullmatch(r"\\left\\\{(.*)\\right\\\}", s)
    return ("S", valeur(m.group(1)))


for k, (p, q) in enumerate([(F(3, 5), F(7, 10)), (F(3, 4), F(2, 3)), (F(9, 10), F(1, 2))]):
    # ensemble des c possibles : les quatre cases de la loi conjointe positives (grille fine de valeurs de c)
    possibles = [F(i, 600) for i in range(601) if all(x >= 0 for x in (F(i, 600), p - F(i, 600), q - F(i, 600), 1 - p - q + F(i, 600)))]
    lo, hi = min(possibles), max(possibles)
    assert lo == max(0, p + q - 1) and hi == min(p, q), f"couples-ex4 v{k} : bornes"
    r = version(CO + "-ex4", k)["reponse_finale"]
    assert len(set(r["propositions"])) == 4
    ens = [intervalle(x) for x in r["propositions"]]
    assert ens[r["bonne"]] == ("I", lo, hi), f"couples-ex4 v{k} : bonne réponse"
    for i, e in enumerate(ens):
        if i != r["bonne"]:
            assert e != ("I", lo, hi), f"couples-ex4 v{k} : leurre juste"
    assert ("S", sp.nsimplify(p * q)) in ens and lo < p * q < hi, f"couples-ex4 v{k} : le cas indépendant n'est qu'un point"

# ------------------------------------------------------------------------------------------ independance-ex1
for k, (px, py) in enumerate([([F(2, 10), F(5, 10), F(3, 10)], [F(4, 10), F(6, 10)]),
                              ([F(1, 10), F(6, 10), F(3, 10)], [F(5, 10), F(5, 10)]),
                              ([F(3, 10), F(3, 10), F(4, 10)], [F(2, 10), F(8, 10)])]):
    u = {(x, y): px[x] * py[y] for x in range(3) for y in range(2)}
    juste = sum(p for w, p in u.items() if sum(w) == 2)
    faux = verifie_rep(IN + "-ex1", k, juste)
    assert memes(faux, [px[1] * py[1], px[2] * py[0], px[1] + py[1] - px[1] * py[1]]), f"independance-ex1 v{k} : leurres"

# ------------------------------------------------------------------------------------------ independance-ex2
PROPS = EX[IN + "-ex2"]["reponse_finale"]["propositions"]
for k, tab in enumerate([[[c(10), c(10)], [c(15), c(15)], [c(25), c(25)]],
                         [[c(10), c(10)], [c(20), c(10)], [c(20), c(30)]],
                         [[c(6), c(14)], [c(12), c(28)], [c(12), c(28)]]]):
    u = {(x, y): tab[x][y] for x in range(3) for y in range(2)}
    assert sum(u.values()) == 1
    ind = indep(u, lambda w: w[0], lambda w: w[1])
    verifie_texte(IN + "-ex2", k, PROPS[0] if ind else PROPS[1])
    assert version(IN + "-ex2", k)["reponse_finale"]["propositions"] == PROPS
# leurre 3 faux : une variable à trois valeurs peut être indépendante d'une variable à deux valeurs (version 0)
assert indep({(x, y): [[c(10), c(10)], [c(15), c(15)], [c(25), c(25)]][x][y] for x in range(3) for y in range(2)}, lambda w: w[0], lambda w: w[1])

# ------------------------------------------------------------------------------------------ independance-ex3
PROPS = EX[IN + "-ex3"]["reponse_finale"]["propositions"]
for k, (f, kk) in enumerate([(6, 3), (12, 6), (6, 2)]):
    contient(IN + "-ex3", k, f"${f}$ faces", f"D\\leq {kk}")
    u = unif(list(range(1, f + 1)))
    X, Y, Z = (lambda w: w % 2), (lambda w: w % 3), (lambda w: 1 if w <= kk else 0)
    assert indep(u, X, Y), f"independance-ex3 v{k} : X, Y"
    xz, yz = indep(u, X, Z), indep(u, Y, Z)
    attendu = PROPS[{(False, False): 0, (True, False): 1, (False, True): 2, (True, True): 3}[(xz, yz)]]
    verifie_texte(IN + "-ex3", k, attendu)
    lt = loi_de(u, lambda w: (X(w), Y(w), Z(w)))
    lx, ly, lz = loi_de(u, X), loi_de(u, Y), loi_de(u, Z)
    mut = all(lt.get((a, b, e), 0) == lx[a] * ly[b] * lz[e] for a in lx for b in ly for e in lz)
    assert mut == (xz and yz), f"independance-ex3 v{k} : question 4"
    correction = version(IN + "-ex3", k)["correction"]
    assert ("**4.** Oui" in correction) == mut, f"independance-ex3 v{k} : conclusion de la question 4"

# ------------------------------------------------------------------------------------------ independance-ex4
assert E(U3, lambda x: x ** 4) == F(2, 3) and E(U3, lambda x: x * x) * E(U3, lambda x: x * x) == F(4, 9), "independance-ex4 : k = 2, l = 1"
for _ in range(50):
    xs = random.sample(range(-6, 7), random.randint(2, 5))
    f = {x: random.randint(-5, 5) for x in xs}
    t = sp.symbols("t")
    Q = sum(f[xi] * sp.prod([(t - xj) / sp.Integer(xi - xj) for xj in xs if xj != xi]) for xi in xs)
    assert all(Q.subs(t, x) == f[x] for x in xs), "independance-ex4 : interpolation de Lagrange"

# ------------------------------------------------------------------------------------------ covariance-ex1
for k, (p00, p01, p10, p11) in enumerate([(c(30), c(20), c(10), c(40)), (c(10), c(40), c(30), c(20)), (c(20), c(10), c(30), c(40))]):
    u = {(0, 0): p00, (0, 1): p01, (1, 0): p10, (1, 1): p11}
    juste = cov(u, fx, fy)
    faux = verifie_rep(CV + "-ex1", k, juste)
    assert memes(faux, [p11, -juste, E(u, fx) * E(u, fy)]), f"covariance-ex1 v{k} : leurres"

# ------------------------------------------------------------------------------------------ covariance-ex2
for k, (vx, vy, cxy) in enumerate([(4, 9, -2), (1, 4, 1), (9, 1, 2)]):
    contient(CV + "-ex2", k, f"$V(X)={vx}$", f"$V(Y)={vy}$", f"\\mathrm{{Cov}}(X,Y)={cxy}$")
    # réalisation explicite : X = aU, Y = cU + dW, avec U, W indépendantes uniformes sur {-1, 1}
    a = sp.sqrt(vx)
    cc = sp.Rational(cxy) / a
    dd = sp.sqrt(vy - cc ** 2)
    u = unif(list(product([-1, 1], repeat=2)))
    Xr = lambda w: a * w[0]
    Yr = lambda w: cc * w[0] + dd * w[1]
    assert sp.simplify(V(u, Xr) - vx) == 0 and sp.simplify(V(u, Yr) - vy) == 0 and sp.simplify(cov(u, Xr, Yr) - cxy) == 0, f"covariance-ex2 v{k} : réalisation"
    juste = sp.simplify(V(u, lambda w: 2 * Xr(w) - 3 * Yr(w) + 1))
    faux = verifie_rep(CV + "-ex2", k, juste)
    assert memes(faux, [4 * vx + 9 * vy, 4 * vx + 9 * vy + 12 * cxy, juste + 1]), f"covariance-ex2 v{k} : leurres"
    assert sp.simplify(V(u, lambda w: Xr(w) + Yr(w)) - (vx + vy + 2 * cxy)) == 0, f"covariance-ex2 v{k} : V(X+Y)"

# ------------------------------------------------------------------------------------------ covariance-ex3
for k, (N, R, n) in enumerate([(8, 3, 4), (10, 5, 4), (12, 4, 3)]):
    contient(CV + "-ex3", k, f"${N}$ boules, dont ${R}$ rouges", f"${n}$ boules sans remise")
    boules = [1] * R + [0] * (N - R)
    u = unif(list(permutations(range(N), n)))
    S = lambda w: sum(boules[i] for i in w)
    juste = V(u, S)
    assert juste == n * F(R, N) * F(N - R, N) * F(N - n, N - 1), f"covariance-ex3 v{k} : variance hypergéométrique"
    cij = cov(u, lambda w: boules[w[0]], lambda w: boules[w[1]])
    assert cij == F(R * (R - 1), N * (N - 1)) - F(R, N) ** 2 < 0, f"covariance-ex3 v{k} : covariance"
    faux = verifie_rep(CV + "-ex3", k, juste)
    vb = n * F(R, N) * F(N - R, N)
    assert memes(faux, [vb, n * F(R, N), vb - n * (n - 1) * cij]), f"covariance-ex3 v{k} : leurres"
    del u

# ------------------------------------------------------------------------------------------ covariance-ex4
u = unif(list(product([-1, 1], repeat=2)))
Xr, Sr = (lambda w: w[0]), (lambda w: w[0] + w[1])
rho = sp.nsimplify(cov(u, Xr, Sr)) / sp.sqrt(sp.nsimplify(V(u, Xr) * V(u, Sr)))
faux = verifie_rep(CV + "-ex4", 0, rho)
assert sp.simplify(rho - sp.sqrt(2) / 2) == 0, "covariance-ex4 : rho(X, X+Y)"
# cas |rho| = 1 pour Y = aX + b
for a, b in [(2, 1), (-3, 5)]:
    u3 = unif([0, 1, 4])
    r = sp.nsimplify(cov(u3, lambda x: x, lambda x: a * x + b)) / sp.sqrt(sp.nsimplify(V(u3) * V(u3, lambda x: a * x + b)))
    assert r == (1 if a > 0 else -1), "covariance-ex4 : |rho| = 1"

# ------------------------------------------------------------------------------------------ covariance-ex5
for n in range(2, 7):
    u = unif(list(permutations(range(n))))
    Nfix = lambda s: sum(s[i] == i for i in range(n))
    assert E(u, Nfix) == 1 and V(u, Nfix) == 1, "covariance-ex5 : E(N) = V(N) = 1"
    assert cov(u, lambda s: int(s[0] == 0), lambda s: int(s[1] == 1)) == F(1, n * n * (n - 1)), "covariance-ex5 : covariance"
    for leurre in (F(n - 1, n), F(n - 2, n), F(1, n)):
        assert leurre != 1, "covariance-ex5 : leurres faux"
r = EX[CV + "-ex5"]["reponse_finale"]
assert r["propositions"][r["bonne"]] == "$1$" and len(set(r["propositions"])) == 4, "covariance-ex5 : bonne réponse"
print("proba_l1_variables_finies : toutes les vérifications passent.")
