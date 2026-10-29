# Vérification calculatoire des questions de QCM (content/qcm/proba-l1-variables-finies.json) et des exercices ajoutés au
# chapitre proba-l1-variables-finies (proba-l1-variables-finies-transfert-ex5, -independance-ex5, -couples-ex5, dans
# content/exercices/probabilites/l1-proba-l1-variables-finies.json). Calculs exacts en fractions, par énumération.
# Chaque assert cite l'identifiant ; « v1 », « v2 » désignent les variantes.
import json
import pathlib
import re
from fractions import Fraction as F
from itertools import combinations, product

RACINE = pathlib.Path(__file__).resolve().parent.parent
QCM = {q["id"]: q for q in json.loads((RACINE / "content/qcm/proba-l1-variables-finies.json").read_text())["questions"]}
EX = {e["id"]: e for e in json.loads((RACINE / "content/exercices/probabilites/l1-proba-l1-variables-finies.json").read_text())["exercices"]}


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


def reponse(ident, k):
    e = EX[ident]
    return e["reponse_finale"] if k == 0 else e["variantes"][k - 1]["reponse_finale"]


def rep(ident, k, attendu):
    r = reponse(ident, k)
    choix(r["propositions"], r["bonne"], attendu)


def E(loi, f=lambda w: w):
    return sum(f(w) * p for w, p in loi.items())


# q1 : X uniforme sur {−2, …, 2} : E(X²) = 2.
L = {x: F(1, 5) for x in range(-2, 3)}
qcm("proba-l1-variables-finies-q1", E(L, lambda x: x * x))

# q2 : valeurs −1, 0, 2 de probabilités 1/4, 1/4, 1/2 : V(X) = 27/16 ; leurres E(X²), E(X²) + E(X)², E(X²) − E(X).
L = {-1: F(1, 4), 0: F(1, 4), 2: F(1, 2)}
m = E(L)
V = E(L, lambda x: (x - m) ** 2)
qcm("proba-l1-variables-finies-q2", V)
assert [valeur(p) for p in QCM["proba-l1-variables-finies-q2"]["propositions"]][1:] == [E(L, lambda x: x * x), E(L, lambda x: x * x) + m * m, E(L, lambda x: x * x) - m]

# q3 (vrai) : X premier dé, Y = 1 si la somme vaut 7 : indépendantes ; avec la somme 8, non.
des = {w: F(1, 36) for w in product(range(1, 7), repeat=2)}
def indep(f, g, loi):
    vf = {f(w) for w in loi}
    vg = {g(w) for w in loi}
    return all(sum(p for w, p in loi.items() if f(w) == x and g(w) == y) == sum(p for w, p in loi.items() if f(w) == x) * sum(p for w, p in loi.items() if g(w) == y) for x in vf for y in vg)
assert indep(lambda w: w[0], lambda w: int(w[0] + w[1] == 7), des)
assert not indep(lambda w: w[0], lambda w: int(w[0] + w[1] == 8), des)
assert QCM["proba-l1-variables-finies-q3"]["vrai"] is True

# q4 (vrai) : V(X + Y) = V(X) + V(Y) + 2 Cov ; exemple X, Y à valeurs ±1 avec P(X = Y) = 3/4.
L = {(1, 1): F(3, 8), (-1, -1): F(3, 8), (1, -1): F(1, 8), (-1, 1): F(1, 8)}
var_ = lambda f: E(L, lambda w: f(w) ** 2) - E(L, f) ** 2
cov = E(L, lambda w: w[0] * w[1]) - E(L, lambda w: w[0]) * E(L, lambda w: w[1])
assert var_(lambda w: w[0]) == 1 and var_(lambda w: w[1]) == 1 and var_(lambda w: w[0] + w[1]) == 3 and cov == F(1, 2)
assert (3 - 1 - 1) / F(2) == F(1, 2) and QCM["proba-l1-variables-finies-q4"]["vrai"] is True

# q5 (vrai) : marginales de Bernoulli(1/2) et c = P(X=1, Y=1) déterminent la loi conjointe, pour tout c de [0, 1/2].
for c in [F(0), F(1, 7), F(1, 4), F(1, 2)]:
    loi = {(1, 1): c, (1, 0): F(1, 2) - c, (0, 1): F(1, 2) - c, (0, 0): c}
    assert all(p >= 0 for p in loi.values()) and sum(loi.values()) == 1
    assert loi[(1, 1)] + loi[(1, 0)] == F(1, 2) and loi[(1, 1)] + loi[(0, 1)] == F(1, 2)
assert QCM["proba-l1-variables-finies-q5"]["vrai"] is True

# q6 (faux) : somme de deux dés uniformes indépendants non uniforme sur {2, …, 12}.
assert sum(p for w, p in des.items() if sum(w) == 2) == F(1, 36) and sum(p for w, p in des.items() if sum(w) == 7) == F(1, 6)
assert QCM["proba-l1-variables-finies-q6"]["vrai"] is False

# q7 : Cov(3X + 1, 5 − Y) = −3 Cov(X, Y) = −6 ; exemple X uniforme sur {0, 2}, Y = 2X (Cov(X, Y) = 2).
L = {0: F(1, 2), 2: F(1, 2)}
c = lambda f, g: E(L, lambda x: f(x) * g(x)) - E(L, f) * E(L, g)
assert c(lambda x: x, lambda x: 2 * x) == 2
qcm("proba-l1-variables-finies-q7", c(lambda x: 3 * x + 1, lambda x: 5 - 2 * x))

# q8 : plus petit de deux jetons tirés parmi 1..5 : E(X) = 2.
tir = list(combinations(range(1, 6), 2))
qcm("proba-l1-variables-finies-q8", F(sum(min(t) for t in tir), len(tir)))

# proba-l1-variables-finies-couples-ex5 : loi de D = |X − Y| et E(D) = (n² − 1)/(3n), vérifiées pour 2 ≤ n ≤ 15.
for n in range(2, 16):
    loi = {(x, y): F(1, n * n) for x in range(1, n + 1) for y in range(1, n + 1)}
    assert sum(p for (x, y), p in loi.items() if x < y) == F(n - 1, 2 * n)
    assert sum(p for (x, y), p in loi.items() if x == y) == F(1, n)
    for k in range(1, n):
        assert sum(p for (x, y), p in loi.items() if abs(x - y) == k) == F(2 * (n - k), n * n)
    assert E(loi, lambda w: abs(w[0] - w[1])) == F(n * n - 1, 3 * n)
for v, n in enumerate([6, 4, 10]):
    rep("proba-l1-variables-finies-couples-ex5", v, F(n * n - 1, 3 * n))
    # leurres : E(X − Y) = 0, (n − 1)/2, la moitié de la bonne valeur
    assert [valeur(p) for p in reponse("proba-l1-variables-finies-couples-ex5", v)["propositions"]][1:] == [0, F(n - 1, 2), F(n * n - 1, 6 * n)]

# proba-l1-variables-finies-transfert-ex5 : P(X = k) = k/S : E(X) et E(1/X) = m/S ≠ 1/E(X).
for v, mm in enumerate([4, 5, 6]):
    S = mm * (mm + 1) // 2
    loi = {k: F(k, S) for k in range(1, mm + 1)}
    assert sum(loi.values()) == 1
    EX_, Einv = E(loi), E(loi, lambda k: F(1, k))
    assert Einv != 1 / EX_
    rep("proba-l1-variables-finies-transfert-ex5", v, Einv)
    # leurres : 1/E(X), moyenne des 1/k, 1/S
    assert [valeur(p) for p in reponse("proba-l1-variables-finies-transfert-ex5", v)["propositions"]][1:] == [1 / EX_, sum(F(1, k) for k in range(1, mm + 1)) / mm, F(1, S)]

# proba-l1-variables-finies-independance-ex5 : énumération de l'urne pour 1 ≤ n ≤ 60 ; ensemble des n qui conviennent.
def ensemble(r, d):
    ok = set()
    for n in range(1, 61):
        boules = {("R", k): F(1, n + r) for k in range(1, r + 1)} | {("B", k): F(1, n + r) for k in range(1, n + 1)}
        if indep(lambda w: int(w[0] == "R"), lambda w: int(w[1] % d == 0), boules):
            ok.add(n)
    return ok
N = set(range(1, 61))
for v, (r, d, bon, faux) in enumerate([
        (4, 2, {n for n in N if n % 2 == 0}, [{n for n in N if n % 2}, {n for n in N if n % 4 == 0}, set()]),
        (6, 3, {n for n in N if n % 3 == 0}, [{n for n in N if n % 6 == 0}, {n for n in N if n % 3}, set()]),
        (5, 5, {n for n in N if n % 5 == 0}, [{n for n in N if n % 10 == 0}, {5}, set()])]):
    ens = ensemble(r, d)
    assert ens == bon and all(f != ens for f in faux)
    rf_ = reponse("proba-l1-variables-finies-independance-ex5", v)
    assert rf_["bonne"] == 0 and len(set(rf_["propositions"])) == 4

print("enrichi_proba_l1_variables_finies : tout est vérifié")
