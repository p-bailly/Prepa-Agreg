# Vérification calculatoire des questions de QCM (content/qcm/proba-l1-conditionnement.json) et des exercices ajoutés au
# chapitre proba-l1-conditionnement (proba-l1-conditionnement-composees-ex5, proba-l1-independance-mutuelle-ex5,
# proba-l1-probas-totales-bayes-ex5, dans content/exercices/probabilites/l1-proba-l1-conditionnement.json). Calculs
# exacts en fractions, par énumération d'univers finis. Chaque assert cite l'identifiant ; « v1 », « v2 » = variantes.
import json
import pathlib
import re
from fractions import Fraction as F
from itertools import permutations, product
from math import log

RACINE = pathlib.Path(__file__).resolve().parent.parent
QCM = {q["id"]: q for q in json.loads((RACINE / "content/qcm/proba-l1-conditionnement.json").read_text())["questions"]}
EX = {e["id"]: e for e in json.loads((RACINE / "content/exercices/probabilites/l1-proba-l1-conditionnement.json").read_text())["exercices"]}


def valeur(tex):
    s = tex.strip().strip("$").replace("{,}", ".").replace(" ", "")
    m = re.fullmatch(r"(-?)\\frac\{(\d+)\}\{(\d+)\}", s)
    if m:
        return (-1 if m.group(1) else 1) * F(int(m.group(2)), int(m.group(3)))
    if re.fullmatch(r"-?\d+(\.\d+)?", s):
        return F(s)
    raise ValueError(tex)


def choix(props, bonne, attendu):
    vals = [valeur(p) for p in props]
    assert vals[bonne] == attendu, (props, attendu)
    assert len(set(vals)) == 4, props
    assert all(v != attendu for i, v in enumerate(vals) if i != bonne), props


def qcm(ident, attendu):
    q = QCM[ident]
    choix(q["propositions"], q["bonne"], attendu)


def rep(ident, k, attendu):
    e = EX[ident]
    r = e["reponse_finale"] if k == 0 else e["variantes"][k - 1]["reponse_finale"]
    choix(r["propositions"], r["bonne"], attendu)


# q1 : prévalence 1 %, sensibilité 99 %, faux positifs 2 % : P(M | T) = 1/3.
pM, pTM, pTS = F(1, 100), F(99, 100), F(2, 100)
pT = pTM * pM + pTS * (1 - pM)
assert pT == F(297, 10000)
qcm("proba-l1-conditionnement-q1", pTM * pM / pT)

# q2 : 5 rouges, 3 noires, sans remise : première noire au 3e tirage (énumération des 8! ordres).
boules = "RRRRRNNN"
ordres = list(permutations(range(8)))
fav = sum(1 for o in ordres if boules[o[0]] == "R" and boules[o[1]] == "R" and boules[o[2]] == "N")
qcm("proba-l1-conditionnement-q2", F(fav, len(ordres)))
assert F(fav, len(ordres)) == F(5, 28) and F(5, 8) ** 2 * F(3, 8) == F(75, 512)

# q3 (vrai) : jeu de 32 cartes, « cœur » et « as » indépendants.
cartes = list(product(range(4), range(8)))  # (couleur, rang) ; couleur 0 = cœur, rang 0 = as
P = lambda ev: F(sum(1 for c in cartes if ev(c)), 32)
assert P(lambda c: c[0] == 0 and c[1] == 0) == P(lambda c: c[0] == 0) * P(lambda c: c[1] == 0)
assert QCM["proba-l1-conditionnement-q3"]["vrai"] is True

# q4 (faux) : P(B) = 0,5 P(A) + 0,3 (1 − P(A)) dépend de P(A) ; P(A) = 0,2 donne 0,34.
pB = lambda a: F(1, 2) * a + F(3, 10) * (1 - a)
assert pB(F(1, 5)) == F(34, 100) != F(4, 10) and pB(F(1, 2)) == F(4, 10)
assert QCM["proba-l1-conditionnement-q4"]["vrai"] is False

# q5 (faux) : trois événements mutuellement indépendants de probabilité 1/3 : P(réunion) = 19/27.
issues = list(product(range(3), repeat=3))  # A_i = {coordonnée i nulle}, équiprobabilité sur 27 issues
assert all(F(sum(1 for w in issues if all(w[i] == 0 for i in I)), 27) == F(1, 3) ** len(I) for I in [(0, 1), (0, 2), (1, 2), (0, 1, 2)])
assert F(sum(1 for w in issues if 0 in w), 27) == F(19, 27) != 1
assert QCM["proba-l1-conditionnement-q5"]["vrai"] is False

# q6 : A, B indépendants, P(A) = 0,4, P(A ∪ B) = 0,7 : P(B) = 0,5.
a = F(2, 5)
b = (F(7, 10) - a) / (1 - a)
assert a + b - a * b == F(7, 10)
qcm("proba-l1-conditionnement-q6", b)

# q7 : transfert d'une boule de U1 (3 R, 2 B) vers U2 (1 R, 4 B), puis tirage dans U2 (énumération).
U1 = "RRRBB"
U2 = "RBBBB"
tot = F(0)
for x in U1:
    u2 = U2 + x
    tot += F(1, 5) * F(u2.count("R"), 6)
qcm("proba-l1-conditionnement-q7", tot)
assert tot == F(4, 15)

# q8 : plus petit n tel que 1 − (5/6)^n > 0,9 : n = 13.
n = 1
while not (1 - F(5, 6) ** n > F(9, 10)):
    n += 1
assert n == 13 and 12 < log(10) / log(1.2) < 13
qcm("proba-l1-conditionnement-q8", F(n))

# proba-l1-conditionnement-composees-ex5 : énumération des tirages ordonnés de 3 pièces parmi N.
for v, (N, d, att) in enumerate([(20, 4, F(8, 57)), (24, 6, F(153, 1012)), (18, 3, F(35, 272))]):
    def_ = set(range(d))
    tirages = list(permutations(range(N), 3))
    fav = sum(1 for t in tirages if t[0] not in def_ and t[1] not in def_ and t[2] in def_)
    p = F(fav, len(tirages))
    assert p == att
    rep("proba-l1-conditionnement-composees-ex5", v, p)
    q = F(N - d, N)
    # leurres : avec remise, sans tenir compte de l'ordre (hypergéométrique), probabilité que la 3e soit défectueuse
    r = EX["proba-l1-conditionnement-composees-ex5"]["reponse_finale"] if v == 0 else EX["proba-l1-conditionnement-composees-ex5"]["variantes"][v - 1]["reponse_finale"]
    assert [valeur(x) for x in r["propositions"]][1:] == [q * q * (1 - q), 3 * att, F(d, N)]

# proba-l1-independance-mutuelle-ex5 : énumération des 8 états des composants (indépendance mutuelle).
for v, (p1, p2, p3, PS, cond) in enumerate([(F(9, 10), F(8, 10), F(7, 10), F(846, 1000), F(40, 47)),
                                            (F(95, 100), F(6, 10), F(5, 10), F(76, 100), F(3, 4)),
                                            (F(8, 10), F(9, 10), F(9, 10), F(792, 1000), F(10, 11))]):
    ps = [p1, p2, p3]
    proba = lambda w: (ps[0] if w[0] else 1 - ps[0]) * (ps[1] if w[1] else 1 - ps[1]) * (ps[2] if w[2] else 1 - ps[2])
    etats = list(product([0, 1], repeat=3))
    S = lambda w: w[0] and (w[1] or w[2])
    pS = sum(proba(w) for w in etats if S(w))
    pF2S = sum(proba(w) for w in etats if S(w) and w[1])
    assert pS == PS and pF2S / pS == cond
    rep("proba-l1-independance-mutuelle-ex5", v, pS)
    # leurres : trois composants en série, trois en parallèle, C1 oublié
    r = EX["proba-l1-independance-mutuelle-ex5"]["reponse_finale"] if v == 0 else EX["proba-l1-independance-mutuelle-ex5"]["variantes"][v - 1]["reponse_finale"]
    assert [valeur(x) for x in r["propositions"]][1:] == [p1 * p2 * p3, 1 - (1 - p1) * (1 - p2) * (1 - p3), 1 - (1 - p2) * (1 - p3)]

# proba-l1-probas-totales-bayes-ex5 : récurrence exacte, forme close, limite, p_2 et P(G_1 | G_2).
for v, (pg, pp, p1, a, b, ell, u1, p2, bayes) in enumerate([
        (F(6, 10), F(3, 10), F(1, 2), F(3, 10), F(3, 10), F(3, 7), F(1, 14), F(45, 100), F(2, 3)),
        (F(7, 10), F(2, 10), F(1, 2), F(5, 10), F(2, 10), F(2, 5), F(1, 10), F(45, 100), F(7, 9)),
        (F(8, 10), F(4, 10), F(1, 4), F(4, 10), F(4, 10), F(2, 3), F(-5, 12), F(5, 10), F(2, 5))]):
    # univers des deux premières parties
    loi12 = {(1, 1): p1 * pg, (1, 0): p1 * (1 - pg), (0, 1): (1 - p1) * pp, (0, 0): (1 - p1) * (1 - pp)}
    assert sum(loi12.values()) == 1
    assert loi12[(1, 1)] + loi12[(0, 1)] == p2
    assert loi12[(1, 1)] / p2 == bayes
    assert pg - pp == a and pp == b and ell == b / (1 - a) and u1 == p1 - ell
    p = p1
    for n in range(1, 40):
        assert 0 < p < 1 and p == ell + u1 * a ** (n - 1)
        suivant = pg * p + pp * (1 - p)  # probabilités totales
        assert suivant == a * p + b  # question 2
        p = suivant
    assert abs(float(p) - float(ell)) < 1e-12
    rep("proba-l1-probas-totales-bayes-ex5", v, ell)

print("enrichi_proba_l1_conditionnement : tout est vérifié")


# --- Relecture indépendante (2026-10-29) : résultats intermédiaires des corrections, recalculés à part. ---
from fractions import Fraction as _F
from itertools import permutations as _perm, product as _prod

# proba-l1-conditionnement-composees-ex5 (et v1, v2) : énumération des tirages ordonnés de 3 pièces distinctes :
# 8/57, 153/1012, 35/272 ; le leurre « avec remise » et le leurre « dans un ordre quelconque » (×3) sont faux.
for (N, D, att) in [(20, 4, _F(8, 57)), (24, 6, _F(153, 1012)), (18, 3, _F(35, 272))]:
    tir = list(_perm(range(N), 3))
    fav = sum(1 for t in tir if t[0] >= D and t[1] >= D and t[2] < D)  # pièces 0..D−1 défectueuses
    assert _F(fav, len(tir)) == att and att != _F(N - D, N) ** 2 * _F(D, N) and 3 * att < 1
# proba-l1-independance-mutuelle-ex5 (et v1, v2) : énumération des 8 états : P(S) et P(F_2 | S) = 40/47, 3/4, 10/11.
for (p1, p2, p3, pS, pc) in [(_F(9, 10), _F(8, 10), _F(7, 10), _F(846, 1000), _F(40, 47)),
                             (_F(95, 100), _F(6, 10), _F(5, 10), _F(76, 100), _F(3, 4)),
                             (_F(8, 10), _F(9, 10), _F(9, 10), _F(792, 1000), _F(10, 11))]:
    tot = s2 = _F(0)
    for e in _prod([0, 1], repeat=3):
        w = (p1 if e[0] else 1 - p1) * (p2 if e[1] else 1 - p2) * (p3 if e[2] else 1 - p3)
        if e[0] and (e[1] or e[2]):
            tot += w
            s2 += w if e[1] else 0
    assert tot == pS and s2 / tot == pc
# proba-l1-probas-totales-bayes-ex5 (et v1, v2) : p_2, forme close p_n = ℓ + (p_1 − ℓ) a^(n−1) vérifiée sur 30 termes
# par la récurrence exacte, limite ℓ = 3/7, 2/5, 2/3 et P(G_1 | G_2) = 2/3, 7/9, 2/5.
for (p1, g, q, p2_att, ell_att, bayes_att) in [(_F(1, 2), _F(6, 10), _F(3, 10), _F(45, 100), _F(3, 7), _F(2, 3)),
                                               (_F(1, 2), _F(7, 10), _F(2, 10), _F(45, 100), _F(2, 5), _F(7, 9)),
                                               (_F(1, 4), _F(8, 10), _F(4, 10), _F(1, 2), _F(2, 3), _F(2, 5))]:
    a, b = g - q, q
    ell = b / (1 - a)
    assert ell == ell_att
    p = p1
    for n in range(1, 31):
        assert p == ell + (p1 - ell) * a ** (n - 1) and 0 < p < 1
        if n == 1:
            p2 = g * p + q * (1 - p)
            assert p2 == p2_att and g * p1 / p2 == bayes_att
        p = g * p + q * (1 - p)
print("enrichi_proba_l1_conditionnement : relecture indépendante vérifiée.")
