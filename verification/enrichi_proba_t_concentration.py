# Vérification calculatoire des questions de QCM (content/qcm/proba-t-concentration.json) et des exercices ajoutés au
# chapitre proba-t-concentration (proba-t-concentration-linearite-ex5, proba-t-concentration-variance-somme-ex5,
# proba-t-lgn-ex5, dans content/exercices/probabilites/terminale-proba-t-concentration.json). Calculs exacts en
# fractions et par énumération. Chaque assert cite l'identifiant ; « v1 », « v2 » désignent les variantes.
import json
import pathlib
import re
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, ceil

RACINE = pathlib.Path(__file__).resolve().parent.parent
QCM = {q["id"]: q for q in json.loads((RACINE / "content/qcm/proba-t-concentration.json").read_text())["questions"]}
EX = {e["id"]: e for e in json.loads((RACINE / "content/exercices/probabilites/terminale-proba-t-concentration.json").read_text())["exercices"]}


def valeur(tex):
    s = tex.strip().strip("$").replace("{,}", ".").replace(" ", "")
    m = re.fullmatch(r"(-?)\\frac\{(\d+)\}\{(\d+)\}", s)
    if m:
        return (-1 if m.group(1) else 1) * F(int(m.group(2)), int(m.group(3)))
    if re.fullmatch(r"-?\d+(\.\d+)?", s):
        return F(s)
    raise ValueError(tex)


def choix(props, bonne, attendu):
    """La proposition « bonne » vaut attendu, les trois autres en diffèrent, toutes distinctes."""
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


def binom(n, p):
    return {k: comb(n, k) * p**k * (1 - p) ** (n - k) for k in range(n + 1)}


def esp(loi, f=lambda x: x):
    return sum(f(x) * q for x, q in loi.items())


def var(loi):
    m = esp(loi)
    return esp(loi, lambda x: (x - m) ** 2)


# q1 : X ~ B(20 ; 0,3), E(3X − 5) = 13.
L = binom(20, F(3, 10))
qcm("proba-t-concentration-q1", esp(L, lambda x: 3 * x - 5))
assert esp(L, lambda x: 3 * x - 5) == 13

# q2 : X ~ B(50 ; 0,2), V(−2X + 7) = 32.
L = binom(50, F(1, 5))
loiY = {}
for x, q in L.items():
    loiY[-2 * x + 7] = loiY.get(-2 * x + 7, 0) + q
assert var(L) == 8 and var(loiY) == 32
qcm("proba-t-concentration-q2", var(loiY))

# q3 (vrai) : exemple de X, Y indépendantes d'écarts-types 3 et 4 : σ(X+Y)² = 25.
loiS = {}
for x, y in product([-3, 3], [-4, 4]):
    loiS[x + y] = loiS.get(x + y, 0) + F(1, 4)
assert var({-3: F(1, 2), 3: F(1, 2)}) == 9 and var({-4: F(1, 2), 4: F(1, 2)}) == 16 and var(loiS) == 25
assert QCM["proba-t-concentration-q3"]["vrai"] is True

# q4 (faux) : X uniforme sur {8, 12} : E = 10, V = 4, P(8 < X < 12) = 0 < 1/2.
L = {8: F(1, 2), 12: F(1, 2)}
assert esp(L) == 10 and var(L) == 4 and sum(q for x, q in L.items() if 8 < x < 12) == 0
assert QCM["proba-t-concentration-q4"]["vrai"] is False

# q5 : V/(n δ²) = 9/(400 × 0,25) = 0,09 ; leurres 0,045 (δ), 0,0225 (sans δ²), 0,3 (σ/(√n δ)).
qcm("proba-t-concentration-q5", F(9) / (400 * F(1, 2) ** 2))
assert F(9) / (400 * F(1, 2)) == F(45, 1000) and F(9, 400) == F(225, 10000) and F(3) / (20 * F(1, 2)) == F(3, 10)

# q6 (faux) : doubler n divise le majorant par 2, pas par 4.
for V, n, d in [(1, 1, 1), (F(9), 400, F(1, 2))]:
    assert (V / (2 * n * d**2)) / (V / (n * d**2)) == F(1, 2)
assert QCM["proba-t-concentration-q6"]["vrai"] is False

# q7 : 5 cartes parmi 32 dont 4 as : espérance du nombre d'as = 5/8 (énumération complète).
total = 0
somme = 0
for main in combinations(range(32), 5):
    total += 1
    somme += sum(1 for c in main if c < 4)
assert total == comb(32, 5) and F(somme, total) == F(5, 8)
assert F(comb(31, 4), comb(32, 5)) == F(5, 32)
qcm("proba-t-concentration-q7", F(somme, total))

# q8 : σ = 5, 1 − 25/a² ≥ 0,96 ⇔ a ≥ 25.
a = 25
assert 1 - F(25, a * a) >= F(96, 100) and 1 - F(25, 24 * 24) < F(96, 100)
qcm("proba-t-concentration-q8", F(a))

# proba-t-concentration-linearite-ex5 : E(N_n) = (n − 1)/k vérifié par énumération pour de petits n.
for k in (2, 4, 6):
    for n in (2, 3, 4):
        issues = list(product(range(k), repeat=n))
        e = F(sum(sum(1 for i in range(n - 1) if w[i] == w[i + 1]) for w in issues), len(issues))
        assert e == F(n - 1, k)
for v, (n, k, seuil, mmin) in enumerate([(12, 6, 3, 19), (20, 2, 5, 11), (25, 4, 10, 41)]):
    rep("proba-t-concentration-linearite-ex5", v, F(n - 1, k))
    m = 2
    while F(m - 1, k) < seuil:
        m += 1
    assert m == mmin  # plus petit m tel que E(N_m) ≥ seuil (question 4)

# proba-t-concentration-variance-somme-ex5 : σ(T) = √n σ et minoration de Bienaymé-Tchebychev.
for v, (n, m, s, h, mino) in enumerate([(25, 200, 4, 50, F(84, 100)), (16, 500, 3, 40, F(91, 100)), (36, 250, 5, 60, F(75, 100))]):
    VT = n * s * s
    sT = round(VT**0.5)
    assert sT * sT == VT
    rep("proba-t-concentration-variance-somme-ex5", v, F(sT))
    assert 1 - F(VT, h * h) == mino
    # leurres : n σ (somme des écarts-types), V(T), √n σ²
    r = (EX["proba-t-concentration-variance-somme-ex5"]["reponse_finale"] if v == 0 else EX["proba-t-concentration-variance-somme-ex5"]["variantes"][v - 1]["reponse_finale"])
    assert [valeur(p) for p in r["propositions"]][1:] == [n * s, VT, round(n**0.5) * s * s]

# proba-t-lgn-ex5 : loi exacte de S = (b + m)N − m q, E(S) = 0, V(S), minoration de B.-T., plus petit n.
for v, (nq, k, b, a, delta, bt, nmin) in enumerate([(20, 4, 3, 10, F(1), F(2, 5), 1200), (30, 5, 4, 20, F(2), F(7, 10), 600), (10, 2, 1, 5, F(1, 2), F(3, 5), 800)]):
    LN = binom(nq, F(1, k))
    LS = {}
    for x, q in LN.items():
        s = b * x - (nq - x)
        LS[s] = LS.get(s, 0) + q
    assert esp(LS) == 0
    VS = var(LS)
    assert 1 - VS / a**2 == bt  # question 3
    n = ceil(VS / (delta**2 * F(5, 100)))
    assert VS / (n * delta**2) <= F(5, 100) < VS / ((n - 1) * delta**2) and n == nmin  # question 4
    rep("proba-t-lgn-ex5", v, F(n))
    # leurres : V(N)/0,05 (V(N) à la place de V(S), δ oublié) ; (b + 1)V(N)/0,05 (carré de b + 1 oublié) ; V(S)/δ² (risque oublié)
    r = EX["proba-t-lgn-ex5"]["reponse_finale"] if v == 0 else EX["proba-t-lgn-ex5"]["variantes"][v - 1]["reponse_finale"]
    VN = var(LN)
    assert [valeur(p) for p in r["propositions"]][1:] == [VN / F(5, 100), (b + 1) * VN / F(5, 100), VS / delta**2]

print("enrichi_proba_t_concentration : tout est vérifié")


# --- Relecture indépendante (2026-10-29) : résultats intermédiaires des corrections, recalculés à part. ---
# proba-t-concentration-linearite-ex5 (et v1, v2) : P(A_i) = 1/k par énumération des couples, plus petit m tel que
# (m − 1)/k ≥ seuil : 19, 11, 41.
for (k, seuil, m_attendu) in [(6, 3, 19), (2, 5, 11), (4, 10, 41)]:
    assert F(sum(1 for a, b in product(range(k), repeat=2) if a == b), k * k) == F(1, k)
    assert min(m for m in range(2, 1000) if F(m - 1, k) >= seuil) == m_attendu
# proba-t-concentration-variance-somme-ex5 (et v1, v2) : E(T), V(T) et minorant 1 − V(T)/a² : 0,84 ; 0,91 ; 0,75.
for (mu, s, n, a, minor) in [(200, 4, 25, 50, F(84, 100)), (500, 3, 16, 40, F(91, 100)), (250, 5, 36, 60, F(3, 4))]:
    assert 1 - F(n * s * s, a * a) == minor and n * mu - a > 0
# proba-t-lgn-ex5 (et v1, v2) : loi binomiale exacte de N, S = gain·N − perte·(q − N), E(S) = 0, V(S) = 60, 120, 10,
# minorant de P(|S| < a) : 0,4 ; 0,7 ; 0,6 ; plus petit n tel que V(S)/(n δ²) ≤ 0,05 : 1200, 600, 800.
for (q, p, g, pe, a, d, VS, minor, n_att) in [(20, F(1, 4), 3, 1, 10, 1, 60, F(2, 5), 1200),
                                              (30, F(1, 5), 4, 1, 20, 2, 120, F(7, 10), 600),
                                              (10, F(1, 2), 1, 1, 5, F(1, 2), 10, F(3, 5), 800)]:
    loi = {g * k - pe * (q - k): comb(q, k) * p**k * (1 - p) ** (q - k) for k in range(q + 1)}
    ES = sum(s * w for s, w in loi.items())
    assert ES == 0 and sum((s - ES) ** 2 * w for s, w in loi.items()) == VS
    assert 1 - F(VS) / a**2 == minor
    assert sum(w for s, w in loi.items() if abs(s) < a) >= minor  # la minoration est bien vérifiée par la vraie loi
    assert min(n for n in range(1, 10**5) if F(VS) / (n * F(d) ** 2) <= F(5, 100)) == n_att
# proba-t-concentration-q4 : la variable à deux valeurs 8 et 12 a bien E = 10, V = 4 et P(8 < X < 12) = 0.
assert F(8 + 12, 2) == 10 and F((8 - 10) ** 2 + (12 - 10) ** 2, 2) == 4
print("enrichi_proba_t_concentration : relecture indépendante vérifiée.")
