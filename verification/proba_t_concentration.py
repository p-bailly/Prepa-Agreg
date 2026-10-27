# Vérification calculatoire des cartes « proba-t-concentration-linearite » et « proba-t-concentration-variance-somme »
# et de leurs exercices (fichiers content/cartes/probabilites/terminale-proba-t-concentration.json et
# content/exercices/probabilites/terminale-proba-t-concentration.json). Calculs exacts en fractions, par énumération
# d'univers finis. Chaque assert cite l'identifiant concerné ; « v1 », « v2 » désignent les variantes.
import json
import pathlib
import random
import re
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial

import sympy as sp

RACINE = pathlib.Path(__file__).resolve().parent.parent
EX = {e["id"]: e for e in json.loads((RACINE / "content/exercices/probabilites/terminale-proba-t-concentration.json").read_text())["exercices"]}
NPC = "On ne peut pas conclure sans savoir si $X$ et $Y$ sont indépendantes."


def valeur(tex):
    s = tex.strip().strip("$").replace("{,}", ".").replace(" ", "")
    m = re.fullmatch(r"(-?)\\frac\{(\d+)\}\{(\d+)\}", s)
    if m:
        return (-1 if m.group(1) else 1) * F(int(m.group(2)), int(m.group(3)))
    if re.fullmatch(r"-?\d+(\.\d+)?", s):
        return F(s)
    raise ValueError(tex)


def rep(ident, k, attendu, textes=()):
    """Exercice ident, variante k (0 = énoncé) : la bonne proposition vaut attendu ; les leurres numériques diffèrent ;
    les leurres textuels autorisés sont listés dans textes (ils doivent être faux, ce que vérifie l'appelant)."""
    e = EX[ident] if k == 0 else EX[ident]["variantes"][k - 1]
    r = e["reponse_finale"]
    props = r["propositions"]
    assert len(set(props)) == 4, f"{ident} v{k} : propositions distinctes"
    assert valeur(props[r["bonne"]]) == attendu, f"{ident} v{k} : bonne réponse"
    autres = []
    for i, x in enumerate(props):
        if i == r["bonne"]:
            continue
        if x in textes:
            continue
        autres.append(valeur(x))
        assert autres[-1] != attendu, f"{ident} v{k} : le leurre {x} est juste"
    return autres


def esperance(loi):
    return sum(F(v) * p for v, p in loi.items())


def variance(loi):
    m = esperance(loi)
    return sum((F(v) - m) ** 2 * p for v, p in loi.items())


def loi_de(issues, f):
    d = {}
    for w, p in issues.items():
        d[f(w)] = d.get(f(w), 0) + p
    return d


# ==========================================================================================
# Carte proba-t-concentration-linearite
# ==========================================================================================
random.seed(2)
for _ in range(40):
    N = random.randint(1, 6)
    poids = [random.randint(1, 9) for _ in range(N)]
    P = {w: F(poids[w], sum(poids)) for w in range(N)}
    X = {w: random.randint(-5, 5) for w in range(N)}
    Y = {w: random.randint(-5, 5) for w in range(N)}  # variables quelconques, en général non indépendantes
    a, b = F(random.randint(-4, 4)), F(random.randint(-4, 4))
    EXv = esperance(loi_de(P, lambda w: X[w]))
    EYv = esperance(loi_de(P, lambda w: Y[w]))
    assert sum(X[w] * P[w] for w in P) == EXv, "proba-t-concentration-linearite : espérance comme moyenne sur les issues"
    assert esperance(loi_de(P, lambda w: a * X[w] + b * Y[w])) == a * EXv + b * EYv, "linearite : E(aX+bY)"
    assert esperance(loi_de(P, lambda w: a * X[w] + b)) == a * EXv + b, "linearite : E(aX+b)"
    if all(X[w] <= Y[w] for w in P):
        assert EXv <= EYv, "linearite affirmation 3"
# Leurres des trous (faux en général) : X = Y = 1 constantes, a = 2, b = 3 donne E(aX+bY) = 5
assert 2 * 3 * 1 * 1 != 5 and 4 + 9 != 5 and 2 * 1 + 3 != 2 and 4 * 1 + 3 != 5 and 2 * 1 + 6 != 5, "linearite trous"
# leurre |a|E(X)+|b|E(Y) : a = -1, X = 1 constante : E(-X) = -1 != 1
assert -1 != abs(-1) * 1, "linearite trou 1 : |a|E(X)+|b|E(Y)"
# Exemple : lettres et enveloppes, E(N) = 1
for n in range(1, 7):
    perms = list(permutations(range(n)))
    assert F(sum(sum(s[i] == i for i in range(n)) for s in perms), len(perms)) == 1, "linearite exemple : E(N) = 1"
    assert F(factorial(n - 1), factorial(n)) == F(1, n)
# A_1 et A_2 non indépendants (n = 3) : P(A1 ∩ A2) = 1/6 != 1/9
perms = list(permutations(range(3)))
assert F(sum(s[0] == 0 and s[1] == 1 for s in perms), 6) == F(1, 6) != F(1, 9), "linearite exemple : A_k non indépendants"
# Contre-exemples et affirmations : Bernoulli(1/2) avec Y = X ; X uniforme sur {1,2}
B = {0: F(1, 2), 1: F(1, 2)}
assert esperance({v * v: p for v, p in B.items()}) == F(1, 2) != esperance(B) ** 2 == F(1, 4), "linearite : E(X^2) != E(X)^2, E(XY) != E(X)E(Y)"
assert esperance({1: F(1, 2), F(1, 2): F(1, 2)}) == F(3, 4) and 1 / esperance({1: F(1, 2), 2: F(1, 2)}) == F(2, 3), "linearite : E(1/X)"
assert variance({2 * v: p for v, p in B.items()}) == 4 * variance(B) != 2 * variance(B), "linearite : V(2X) = 4V(X)"
assert esperance({2 * v: p for v, p in B.items()}) == 1 == 2 * esperance(B), "linearite affirmation 1 : E(2X) = E(X) + E(X)"
assert F(1, 2) != F(1, 4), "linearite affirmation 1 : P(X=1,Y=1) = 1/2 != 1/4"

# ------------------------------------------------------------------------------------------ ex1
for k, (ex_, ey, a, b, c, juste, leurres) in enumerate([(3, -2, 2, -5, 4, 20, [0, 16]), (-1, 4, 3, 2, -1, 4, [5, -12]), (5, F(1, 2), 1, -6, 3, 5, [2, 11])]):
    assert a * ex_ + b * ey + c == juste, f"linearite-ex1 v{k}"
    # vérification sur un exemple où X et Y ne sont pas indépendantes (Y fonction affine de X) : même résultat
    Xv = {0: ex_ - 1, 1: ex_ + 1}
    Yv = {0: ey + 3, 1: ey - 3}
    P = {0: F(1, 2), 1: F(1, 2)}
    assert esperance(loi_de(P, lambda w: a * Xv[w] + b * Yv[w] + c)) == juste, f"linearite-ex1 v{k} : variables liées"
    assert sorted(rep("proba-t-concentration-linearite-ex1", k, juste, textes=(NPC,))) == sorted(leurres), f"linearite-ex1 v{k}"
    assert NPC in EX["proba-t-concentration-linearite-ex1"]["reponse_finale"]["propositions"]
# leurres : erreur de signe / oubli de la constante
assert 2 * 3 - 5 * 2 + 4 == 0 and 2 * 3 + 5 * 2 == 16 and -3 + 8 == 5 and -3 - 8 - 1 == -12 and 5 - 3 == 2 and 5 + 3 + 3 == 11

# ------------------------------------------------------------------------------------------ ex2
for k, (ps, leurres) in enumerate([([F("0.9"), F("0.6"), F("0.5")], None), ([F("0.8"), F("0.7"), F("0.3")], None),
                                   ([F("0.5"), F("0.25"), F("0.75"), F("0.9")], None)]):
    juste = sum(ps)
    # quelle que soit la dépendance : exemple indépendant et exemple « comonotone » (une seule uniforme U, R_i = {U < p_i})
    indep = {w: sp.prod([ps[i] if w[i] else 1 - ps[i] for i in range(len(ps))]) for w in product([1, 0], repeat=len(ps))}
    assert esperance(loi_de(indep, sum)) == juste, f"linearite-ex2 v{k} : cas indépendant"
    seuils = sorted(set([F(0)] + ps + [F(1)]))
    comono = {tuple(1 if g < p else 0 for p in ps): d - g for g, d in zip(seuils, seuils[1:])}
    assert esperance(loi_de(comono, sum)) == juste, f"linearite-ex2 v{k} : cas dépendant"
    faux = rep("proba-t-concentration-linearite-ex2", k, juste)
    assert sorted(faux) == sorted([sp.prod(ps), juste / len(ps), len(ps)]), f"linearite-ex2 v{k} : leurres"
assert sp.prod([F("0.9"), F("0.6"), F("0.5")]) == F("0.27") and sp.prod([F("0.8"), F("0.7"), F("0.3")]) == F("0.168")
assert sp.prod([F("0.5"), F("0.25"), F("0.75"), F("0.9")]) == F("0.084375")

# ------------------------------------------------------------------------------------------ ex3 (faces jamais obtenues)
for k, (f, n, pj, en, dist, approx, leurres) in enumerate([(6, 4, F(625, 1296), F(625, 216), F(671, 216), "2.89", [F(625, 1296), F(671, 216), 2]),
                                                           (6, 3, F(125, 216), F(125, 36), F(91, 36), "3.47", [F(125, 216), F(91, 36), 3]),
                                                           (4, 3, F(27, 64), F(27, 16), F(37, 16), "1.69", [F(27, 64), F(37, 16), 1])]):
    tirages = list(product(range(1, f + 1), repeat=n))
    q = F(1, len(tirages))
    assert q * sum(1 not in t for t in tirages) == pj == F(f - 1, f) ** n, f"linearite-ex3 v{k} : P(A_j)"
    en_brut = q * sum(f - len(set(t)) for t in tirages)
    assert en_brut == en == f * pj and f - en == dist, f"linearite-ex3 v{k} : E(N) par énumération"
    assert abs(float(en) - float(approx)) < 5e-3, f"linearite-ex3 v{k} : valeur approchée"
    assert sorted(rep("proba-t-concentration-linearite-ex3", k, en)) == sorted(leurres), f"linearite-ex3 v{k}"
    # les A_j ne sont pas indépendants : P(A_1 ∩ ... ∩ A_f) = 0 != produit
    assert all(len(set(t)) >= 1 for t in tirages) and pj ** f != 0

# ------------------------------------------------------------------------------------------ ex4 (anniversaires)
for k, (n, seuil, m_min, leurres) in enumerate([(30, 1, 28, [F(6, 73), F(174, 73), F(180, 73)]),
                                                (40, 2, 39, [F(8, 73), F(312, 73), F(320, 73)]),
                                                (23, F(1, 2), 20, [F(23, 365), F(506, 365), F(529, 365)])]):
    en = F(comb(n, 2), 365)
    assert sorted(rep("proba-t-concentration-linearite-ex4", k, en)) == sorted(leurres), f"linearite-ex4 v{k}"
    assert sorted(leurres) == sorted([F(n, 365), F(n * (n - 1), 365), F(n * n, 365)]), f"linearite-ex4 v{k} : origine des leurres"
    assert min(m for m in range(1, 200) if F(m * (m - 1), 730) >= seuil) == m_min, f"linearite-ex4 v{k} : question 4"
assert (comb(30, 2), comb(40, 2), comb(23, 2)) == (435, 780, 253) and F(435, 365) == F(87, 73) and F(780, 365) == F(156, 73)
assert (27 * 26, 28 * 27, 38 * 37, 39 * 38, 19 * 18, 20 * 19) == (702, 756, 1406, 1482, 342, 380)
# P(A_ij) = 1/365 ; linéarité vérifiée par énumération complète sur un petit calendrier (5 jours, 4 personnes)
J, n = 5, 4
cal = list(product(range(J), repeat=n))
assert F(sum(sum(t[i] == t[j] for i, j in combinations(range(n), 2)) for t in cal), len(cal)) == F(comb(n, 2), J), "linearite-ex4 : E(N) = C(n,2)/J"
# A_12 et A_13 réalisés implique A_23 réalisé
assert all(t[1] == t[2] for t in cal if t[0] == t[1] and t[0] == t[2]), "linearite-ex4 : non-indépendance"

# ==========================================================================================
# Carte proba-t-concentration-variance-somme
# ==========================================================================================
for _ in range(40):
    # variables indépendantes : construites sur un univers produit
    l1 = {random.randint(-4, 4): F(random.randint(1, 5)) for _ in range(3)}
    l2 = {random.randint(-4, 4): F(random.randint(1, 5)) for _ in range(3)}
    l1 = {v: p / sum(l1.values()) for v, p in l1.items()}
    l2 = {v: p / sum(l2.values()) for v, p in l2.items()}
    P = {(x, y): p * q for x, p in l1.items() for y, q in l2.items()}
    S = loi_de(P, lambda w: w[0] + w[1])
    assert variance(S) == variance(l1) + variance(l2), "proba-t-concentration-variance-somme : additivité"
    assert esperance(loi_de(P, lambda w: w[0] * w[1])) == esperance(l1) * esperance(l2), "variance-somme : E(XY) = E(X)E(Y)"
    a, b = F(random.randint(-4, 4)), F(random.randint(-4, 4))
    assert variance({a * v + b: p for v, p in l1.items()} if a != 0 else {b: 1}) == a ** 2 * variance(l1), "variance-somme : V(aX+b)"
    # formule générale V(X+Y) = V(X) + V(Y) + 2(E(XY) - E(X)E(Y)) pour des variables liées
    Q = {w: F(random.randint(1, 5)) for w in range(5)}
    Q = {w: p / sum(Q.values()) for w, p in Q.items()}
    Xv = {w: random.randint(-3, 3) for w in Q}
    Yv = {w: random.randint(-3, 3) for w in Q}
    EXY = esperance(loi_de(Q, lambda w: Xv[w] * Yv[w]))
    VX, VY = variance(loi_de(Q, lambda w: Xv[w])), variance(loi_de(Q, lambda w: Yv[w]))
    EXv, EYv = esperance(loi_de(Q, lambda w: Xv[w])), esperance(loi_de(Q, lambda w: Yv[w]))
    assert variance(loi_de(Q, lambda w: Xv[w] + Yv[w])) == VX + VY + 2 * (EXY - EXv * EYv), "variance-somme : terme croisé"
# n variables deux à deux indépendantes mais pas mutuellement : X1, X2 de Bernoulli(1/2) indépendantes, X3 = X1 xor X2
P3 = {(x1, x2): F(1, 4) for x1 in (0, 1) for x2 in (0, 1)}
lois = [loi_de(P3, f) for f in (lambda w: w[0], lambda w: w[1], lambda w: w[0] ^ w[1])]
assert variance(loi_de(P3, lambda w: w[0] + w[1] + (w[0] ^ w[1]))) == sum(variance(l) for l in lois) == F(3, 4), "variance-somme : deux à deux suffit"
assert F(sum(1 for w in P3 if w[0] == 1 and w[1] == 1 and (w[0] ^ w[1]) == 1), 4) == 0 != F(1, 8), "variance-somme : X1, X2, X3 pas mutuellement indépendantes"
# Leurres des trous : X = Y Rademacher (même loi, espérance nulle, même univers) : V(2X) = 4 != 2
R = {-1: F(1, 2), 1: F(1, 2)}
assert variance({2 * v: p for v, p in R.items()}) == 4 != 2 * variance(R), "variance-somme trou 1 : leurres faux"
Bi = {(x, y): F(1, 4) for x in (0, 1) for y in (0, 1)}
VS = variance(loi_de(Bi, lambda w: w[0] + w[1]))
assert VS == F(1, 2) and F(1, 16) != VS and 1 != VS and F(1, 2) + 2 * F(1, 4) != VS, "variance-somme trou 2 : V(X)V(Y), (σX+σY)^2, +2E(X)E(Y)"
assert variance({2 * v + 1: p for v, p in B.items()}) == 1 and 2 * F(1, 4) + 1 != 1 and 4 * F(1, 4) + 1 != 1 and 2 * F(1, 4) != 1, "variance-somme trou 3"
# Exemple : deux dés
de = {i: F(1, 6) for i in range(1, 7)}
assert esperance({v * v: p for v, p in de.items()}) == F(91, 6) and variance(de) == F(35, 12), "variance-somme exemple"
deux = loi_de({(a, b): F(1, 36) for a in range(1, 7) for b in range(1, 7)}, sum)
assert variance(deux) == F(35, 6) and variance({2 * v: p for v, p in de.items()}) == F(35, 3), "variance-somme exemple : S et 2X_1"
# Contre-exemples : Y = X, Y = -X ; réciproque fausse (X uniforme sur {-1,0,1}, Y = X^2)
assert variance({0: F(1)}) == 0, "variance-somme : X + (-X) = 0"
U = {-1: F(1, 3), 0: F(1, 3), 1: F(1, 3)}
VXu, VYu = variance(U), variance(loi_de(U, lambda x: x * x))
VXY = variance(loi_de(U, lambda x: x + x * x))
assert (VXu, VYu, VXY) == (F(2, 3), F(2, 9), F(8, 9)) and VXY == VXu + VYu, "variance-somme réciproque fausse"
assert loi_de(U, lambda x: x + x * x) == {0: F(2, 3), 2: F(1, 3)} and F(4, 3) - F(4, 9) == F(8, 9)
assert F(1, 3) * F(2, 3) == F(2, 9) != 0, "variance-somme : X et X^2 non indépendantes"
# Affirmations : V(X-Y), σ(X+Y) pour deux Bernoulli(1/2) indépendantes ; V(sum) = nV(X), V(nX) = n^2 V(X)
D = loi_de(Bi, lambda w: w[0] - w[1])
assert D == {0: F(1, 2), 1: F(1, 4), -1: F(1, 4)} and variance(D) == F(1, 2) and variance(B) - variance(B) == 0, "variance-somme affirmation 1"
assert sp.sqrt(sp.Rational(1, 2)) == sp.sqrt(2) / 2 != 1, "variance-somme affirmation 2"
for n in range(1, 5):
    Pn = {w: F(1, 2 ** n) for w in product((0, 1), repeat=n)}
    assert variance(loi_de(Pn, sum)) == n * variance(B) and variance({n * v: p for v, p in B.items()}) == n * n * variance(B), "variance-somme affirmation 4"

# ------------------------------------------------------------------------------------------ ex1
for k, (vx, vy, a, b, juste, vsum, leurres) in enumerate([(4, 9, 2, -3, 97, 13, [98, 35, -65]), (2, 5, 3, -2, 38, 7, [42, 16, -2]), (1, 3, 1, -4, 49, 4, [47, 13, -47])]):
    assert a ** 2 * vx + b ** 2 * vy == juste and vx + vy == vsum, f"variance-somme-ex1 v{k}"
    assert sorted(rep("proba-t-concentration-variance-somme-ex1", k, juste)) == sorted(leurres), f"variance-somme-ex1 v{k}"
    assert sorted(leurres) == sorted([juste + {0: 1, 1: 4, 2: -2}[k], abs(a) * vx + abs(b) * vy, a ** 2 * vx - b ** 2 * vy]), f"variance-somme-ex1 v{k} : origine des leurres"
# contrôle sur des lois concrètes : X de variance 4 (±2), Y de variance 9 (±3), indépendantes
X = {-2: F(1, 2), 2: F(1, 2)}
Y = {-3: F(1, 2), 3: F(1, 2)}
PXY = {(x, y): p * q for x, p in X.items() for y, q in Y.items()}
assert variance(loi_de(PXY, lambda w: 2 * w[0] - 3 * w[1] + 1)) == 97 and variance(loi_de(PXY, lambda w: w[0] - w[1])) == 13, "variance-somme-ex1"

# ------------------------------------------------------------------------------------------ ex2 (somme de deux dés)
def de_uniforme(f):
    return {i: F(1, f) for i in range(1, f + 1)}

for k, (fa, fb, es, vs, leurres) in enumerate([(6, 4, 6, F(25, 6), [F(35, 6), F(175, 48), F(25, 12)]),
                                               (6, 6, 7, F(35, 6), [F(35, 12), F(1225, 144), F(35, 3)]),
                                               (4, 4, 5, F(5, 2), [F(5, 4), F(25, 16), 5])]):
    A, Bd = de_uniforme(fa), de_uniforme(fb)
    S = loi_de({(a, b): p * q for a, p in A.items() for b, q in Bd.items()}, sum)
    assert esperance(S) == es and variance(S) == vs == variance(A) + variance(Bd), f"variance-somme-ex2 v{k}"
    assert sorted(rep("proba-t-concentration-variance-somme-ex2", k, vs)) == sorted(leurres), f"variance-somme-ex2 v{k}"
    assert variance(A) * variance(Bd) in leurres, f"variance-somme-ex2 v{k} : leurre produit"
assert variance(de_uniforme(4)) == F(5, 4) and esperance({v * v: p for v, p in de_uniforme(4).items()}) == F(15, 2)
assert F(91, 6) - F(49, 4) == F(35, 12) and F(15, 2) - F(25, 4) == F(5, 4)

# ------------------------------------------------------------------------------------------ ex3 (gain à un jeu)
for k, (n, p, gain, juste, eg, sigma, leurres) in enumerate([(20, F(1, 6), 5, 100, 0, sp.Integer(10), [F(50, 3), F(650, 9), F(25, 9)]),
                                                             (30, F(1, 6), 3, F(200, 3), -10, 10 * sp.sqrt(6) / 3, [F(50, 3), F(125, 3), F(25, 6)]),
                                                             (10, F(1, 2), 2, F(45, 2), 5, 3 * sp.sqrt(10) / 2, [F(15, 2), F(25, 2), F(5, 2)])]):
    X = {j: comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(n + 1)}
    G = {gain * j - (n - j): q for j, q in X.items()}
    assert esperance(G) == eg and variance(G) == juste, f"variance-somme-ex3 v{k}"
    assert sp.sqrt(sp.Rational(juste.numerator, juste.denominator) if isinstance(juste, F) else juste) == sigma, f"variance-somme-ex3 v{k} : écart-type"
    VX = variance(X)
    assert VX == n * p * (1 - p)
    attendus = [(gain + 1) * VX, (gain ** 2 + 1) * VX, VX]
    assert sorted(rep("proba-t-concentration-variance-somme-ex3", k, juste)) == sorted(attendus) == sorted(leurres), f"variance-somme-ex3 v{k} : leurres (aV, V(gX)+V(n-X), V(X))"

print("proba_t_concentration : toutes les vérifications passent.")
