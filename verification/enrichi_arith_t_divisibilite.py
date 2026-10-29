# Vérification calculatoire des questions de QCM (content/qcm/arith-t-divisibilite.json) et des exercices ajoutés au
# chapitre arith-t-divisibilite (arith-t-division-euclidienne-ex5, arith-t-divisibilite-combinaisons-ex5,
# arith-t-divisibilite-congruences-ex5, dans content/exercices/arithmetique/terminale-arith-t-divisibilite.json).
# Calculs en entiers et par force brute. Chaque assert cite l'identifiant ; « v1 », « v2 » désignent les variantes.
import json
import pathlib
import re

RACINE = pathlib.Path(__file__).resolve().parent.parent
QCM = {q["id"]: q for q in json.loads((RACINE / "content/qcm/arith-t-divisibilite.json").read_text())["questions"]}
EX = {e["id"]: e for e in json.loads((RACINE / "content/exercices/arithmetique/terminale-arith-t-divisibilite.json").read_text())["exercices"]}


def entier(tex):
    s = tex.strip().strip("$").replace(" ", "")
    assert re.fullmatch(r"-?\d+", s), tex
    return int(s)


def choix(props, bonne, attendu):
    vals = [entier(p) for p in props]
    assert vals[bonne] == attendu, (props, attendu)
    assert len(set(vals)) == 4, props
    assert all(v != attendu for i, v in enumerate(vals) if i != bonne), props


def reponse(ident, k):
    e = EX[ident]
    return e["reponse_finale"] if k == 0 else e["variantes"][k - 1]["reponse_finale"]


# q1 : chiffre des unités de 7^2026 = 9.
q = QCM["arith-t-divisibilite-q1"]
choix(q["propositions"], q["bonne"], pow(7, 2026, 10))

# q2 (vrai) : n impair ⇒ n² ≡ 1 (mod 8).
assert all((n * n) % 8 == 1 for n in range(-1001, 1002, 2)) and QCM["arith-t-divisibilite-q2"]["vrai"] is True

# q3 (faux) : 4^n + 1 divisible par 5 ⇔ n impair ; contre-exemple n = 2 (17).
assert (4**2 + 1) == 17 and 17 % 5 != 0
assert all(((4**n + 1) % 5 == 0) == (n % 2 == 1) for n in range(200))
assert QCM["arith-t-divisibilite-q3"]["vrai"] is False

# q4 (vrai) : 6 | n³ − n ; les restes 0, 1, …, 5 donnent 0, 0, 6, 24, 60, 120.
assert [r**3 - r for r in range(6)] == [0, 0, 6, 24, 60, 120]
assert all((n**3 - n) % 6 == 0 for n in range(-500, 501)) and QCM["arith-t-divisibilite-q4"]["vrai"] is True

# q5 (faux) : a = b = 2 : 4 | a + b et 4 | a − b, mais 4 ∤ a.
a = b = 2
assert (a + b) % 4 == 0 and (a - b) % 4 == 0 and a % 4 != 0 and QCM["arith-t-divisibilite-q5"]["vrai"] is False

# q6 : diviseurs naturels communs de 5n + 3 et 3n + 2 : {1} pour tout n (force brute), jamais 2, 3 ou 19.
for n in range(-300, 301):
    x, y = 5 * n + 3, 3 * n + 2
    communs = {d for d in range(1, max(abs(x), abs(y)) + 1) if x % d == 0 and y % d == 0}
    assert communs == {1}
assert 5 * (3 * 0 + 2) - 3 * (5 * 0 + 3) == 1
q = QCM["arith-t-divisibilite-q6"]
assert q["bonne"] == 0 and q["propositions"][0] == "$\\{1\\}$" and len(set(q["propositions"])) == 4

# q7 : 256 = 9b + r avec 0 ≤ r < b : b ∈ {26, 27, 28}.
bs = [b for b in range(1, 300) if 256 // b == 9]
assert bs == [26, 27, 28] and [256 % b for b in bs] == [22, 13, 4]
q = QCM["arith-t-divisibilite-q7"]
choix(q["propositions"], q["bonne"], len(bs))

# q8 : 7 | 2^n − 1 ⇔ 3 | n ; les autres caractérisations (multiples de 6, pairs, multiples de 7) sont fausses.
ens = {n for n in range(300) if (2**n - 1) % 7 == 0}
assert ens == {n for n in range(300) if n % 3 == 0}
for autre in ({n for n in range(300) if n % 6 == 0}, {n for n in range(300) if n % 2 == 0}, {n for n in range(300) if n % 7 == 0}):
    assert autre != ens
assert QCM["arith-t-divisibilite-q8"]["bonne"] == 0

# arith-t-division-euclidienne-ex5 : entiers a ≥ 1 dont quotient = reste dans la division par b (force brute).
for v, b in enumerate([23, 17, 31]):
    sols = [a for a in range(1, 5000) if a // b == a % b]
    assert sols == [(b + 1) * r for r in range(1, b)]
    choix(reponse("arith-t-division-euclidienne-ex5", v)["propositions"], reponse("arith-t-division-euclidienne-ex5", v)["bonne"], max(sols))
    # leurres : (b + 1) b, b², b(b − 1)
    assert [entier(p) for p in reponse("arith-t-division-euclidienne-ex5", v)["propositions"]][1:] == [(b + 1) * b, b * b, b * (b - 1)]

# arith-t-divisibilite-combinaisons-ex5 : couples d'entiers naturels tels que x² − y² = N (force brute).
for v, (N, att) in enumerate([(45, [(7, 2), (9, 6), (23, 22)]), (60, [(8, 2), (16, 14)]), (105, [(11, 4), (13, 8), (19, 16), (53, 52)])]):
    sols = [(x, y) for x in range(0, N + 2) for y in range(0, x + 1) if x * x - y * y == N]
    assert sols == att
    choix(reponse("arith-t-divisibilite-combinaisons-ex5", v)["propositions"], reponse("arith-t-divisibilite-combinaisons-ex5", v)["bonne"], len(sols))

# arith-t-divisibilite-congruences-ex5 : 7 | 3^n + 4^n ⇔ n impair, etc. ; restes de A_2026.
for v, (a, b, p, ta, tb, rest) in enumerate([(3, 4, 7, 6, 3, 1), (2, 3, 5, 4, 4, 3), (5, 8, 13, 4, 4, 11)]):
    assert pow(a, ta, p) == 1 and pow(b, tb, p) == 1 and (a + b) % p == 0  # questions 1, 2 et 5
    assert all(pow(a, k, p) != 1 for k in range(1, ta))  # ta est la période des puissances de a
    assert all(((pow(a, n, p) + pow(b, n, p)) % p == 0) == (n % 2 == 1) for n in range(400))  # question 3
    assert (pow(a, 2026, p) + pow(b, 2026, p)) % p == rest  # question 4
    choix(reponse("arith-t-divisibilite-congruences-ex5", v)["propositions"], reponse("arith-t-divisibilite-congruences-ex5", v)["bonne"], rest)
    # tableau de la question 3 : pour n pair, A_n ≢ 0 (2 a^n vaut 2, 4, 1 modulo 7, etc.)
    assert all((2 * pow(a, r, p)) % p != 0 for r in range(0, ta, 2))

print("enrichi_arith_t_divisibilite : tout est vérifié")


# --- Relecture indépendante (2026-10-29) : résultats intermédiaires des corrections, recalculés à part. ---
# arith-t-division-euclidienne-ex5 (et v1, v2) : quotient = reste dans la division par b ⇔ a = (b + 1) r, 1 ≤ r ≤ b − 1 ;
# 22, 16, 30 entiers, le plus grand 528, 288, 960.
for (b, nb, gd) in [(23, 22, 528), (17, 16, 288), (31, 30, 960)]:
    sols = [a for a in range(1, 5000) if a // b == a % b]
    assert sols == [(b + 1) * r for r in range(1, b)] and len(sols) == nb and max(sols) == gd
# arith-t-divisibilite-combinaisons-ex5 (et v1, v2) : solutions exactes de x² − y² = N.
for (N, att) in [(45, {(23, 22), (9, 6), (7, 2)}), (60, {(16, 14), (8, 2)}),
                 (105, {(53, 52), (19, 16), (13, 8), (11, 4)})]:
    assert {(x, y) for x in range(N + 1) for y in range(x + 1) if x * x - y * y == N} == att
# arith-t-divisibilite-congruences-ex5 (et v1, v2) : listes des restes, périodes, divisibilité ⇔ n impair, reste de A_2026.
for (u, v, m, per_u, per_v, restes_u, r2026) in [(3, 4, 7, 6, 3, [1, 3, 2, 6, 4, 5, 1], 1),
                                                 (2, 3, 5, 4, 4, [1, 2, 4, 3, 1], 3),
                                                 (5, 8, 13, 4, 4, [1, 5, 12, 8, 1], 11)]:
    assert [pow(u, k, m) for k in range(len(restes_u))] == restes_u
    assert pow(u, per_u, m) == 1 and pow(v, per_v, m) == 1 and (u + v) % m == 0
    assert all(((u**n + v**n) % m == 0) == (n % 2 == 1) for n in range(200))
    assert (pow(u, 2026, m) + pow(v, 2026, m)) % m == r2026
# arith-t-divisibilite-q7 : quotient 9 dans la division de 256 par b : b ∈ {26, 27, 28}, restes 22, 13, 4.
assert [(b, 256 % b) for b in range(1, 257) if 256 // b == 9] == [(26, 22), (27, 13), (28, 4)]
# arith-t-divisibilite-q8 : 7 | 2^n − 1 ⇔ 3 | n (n ≤ 300) ; 2^3 − 1 = 7 contredit « pairs », « multiples de 6 ou de 7 ».
assert all(((2**n - 1) % 7 == 0) == (n % 3 == 0) for n in range(301))
print("enrichi_arith_t_divisibilite : relecture indépendante vérifiée.")
