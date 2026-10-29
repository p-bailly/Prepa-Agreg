"""Vérifications de l'enrichissement du chapitre « proba-t-binomiale » (QCM et exercices ajoutés), calcul exact.

Fichiers content/qcm/proba-t-binomiale.json et fin de content/exercices/probabilites/terminale-proba-t-binomiale.json
(proba-t-binomiale-ex5, proba-t-binomiale-pascal-ex5, proba-t-binomiale-produit-ex5). Chaque assert cite
l'identifiant de la question ou de l'exercice (variantes : « v1 », « v2 ») ; les leurres sont vérifiés FAUX.
Lancer : python3 verification/enrichi_proba_t_binomiale.py (doit finir sans erreur).
"""

import math
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, perm


def loi_binomiale(n, p):
    return [comb(n, k) * p**k * (1 - p) ** (n - k) for k in range(n + 1)]


# proba-t-binomiale-q1 : P(X = 2) pour B(10, 0,3) = C(10,2) 0,3^2 0,7^8 ; les leurres donnent d'autres valeurs.
p = F(3, 10)
bonne = comb(10, 2) * p**2 * (1 - p) ** 8
assert loi_binomiale(10, p)[2] == bonne
assert all(v != bonne for v in (p**2 * (1 - p) ** 8, comb(10, 2) * p**8 * (1 - p) ** 2, comb(10, 2) * p**2))
assert round(float(bonne), 4) == 0.2335

# proba-t-binomiale-q2 : np = 12, np(1-p) = 3 => p = 3/4, n = 16 ; leurres faux.
def EV(n, p):
    return n * p, n * p * (1 - p)


assert EV(16, F(3, 4)) == (12, 3)
assert all(EV(n, q) != (12, 3) for n, q in ((48, F(1, 4)), (16, F(1, 4)), (4, F(3, 4))))

# proba-t-binomiale-q3 (vrai) : B(20, 1/2) symétrique : P(X ≤ 9) = P(X ≥ 11).
L = loi_binomiale(20, F(1, 2))
assert sum(L[:10]) == sum(L[11:])

# proba-t-binomiale-q4 : au moins un 6 en 4 lancers = 1 - (5/6)^4 = 671/1296, par énumération de {1..6}^4.
issues = list(product(range(1, 7), repeat=4))
assert F(sum(1 for t in issues if 6 in t), len(issues)) == 1 - F(5, 6) ** 4 == F(671, 1296)
assert all(v != F(671, 1296) for v in (F(4, 6), F(1, 6) ** 4, 1 - F(1, 6) ** 4))

# proba-t-binomiale-q5 : mains de 5 cartes parmi 32 avec exactement 2 cœurs (8 cœurs) = C(8,2) C(24,3) = 56672.
bonne = comb(8, 2) * comb(24, 3)
assert bonne == 56672
assert all(v != bonne for v in (comb(8, 2) * comb(32, 3), comb(8, 2) + comb(24, 3), 8 * 7 * comb(24, 3)))
# énumération sur un modèle réduit (12 cartes dont 4 cœurs, mains de 4, exactement 2 cœurs) pour la méthode
cartes = range(12)
coeurs = set(range(4))
assert sum(1 for m in combinations(cartes, 4) if len(coeurs & set(m)) == 2) == comb(4, 2) * comb(8, 2)

# proba-t-binomiale-q6 (vrai) : C(n+1, 2) - C(n, 2) = n.
assert all(comb(n + 1, 2) - comb(n, 2) == n for n in range(2, 200))

# proba-t-binomiale-q7 (faux) : parties à 3 éléments de {1..10} contenant 1 : C(9,2) = 36 ≠ C(9,3) = 84.
nb = sum(1 for s in combinations(range(1, 11), 3) if 1 in s)
assert nb == comb(9, 2) == 36 and comb(9, 3) == 84 and nb != comb(9, 3)

# proba-t-binomiale-q8 : C(21, 9) = C(20, 8) + C(20, 9) = 293930 ; leurres faux.
assert comb(20, 8) == 125970 and comb(20, 9) == 167960 and comb(21, 9) == 293930
assert all(v != 293930 for v in (41990, 335920, 251940))
assert comb(20, 9) - comb(20, 8) == 41990

# proba-t-binomiale-ex5 : contrôle de qualité. (p, n, seuil, attendu, P(X=0), P(X≥2), V)
cas = [
    (F(4, 100), 50, F(95, 100), 74, 0.1299, 0.2706, 0.5995, F(192, 100)),  # exercice
    (F(5, 100), 40, F(99, 100), 90, 0.1285, 0.2706, 0.6009, F(190, 100)),  # v1
    (F(2, 100), 100, F(90, 100), 114, 0.1326, 0.2707, 0.5967, F(196, 100)),  # v2
]
for p, n, seuil, attendu, p0, p1, p2, V in cas:
    L = loi_binomiale(n, p)
    assert round(float(L[0]), 4) == p0
    assert round(float(L[1]), 4) == p1
    assert round(float(1 - L[0] - L[1]), 4) == p2
    assert n * p == 2 and n * p * (1 - p) == V
    q = 1 - p
    # plus petit n tel que 1 - q^n ≥ seuil, en calcul exact
    m = 1
    while 1 - q**m < seuil:
        m += 1
    assert m == attendu
    assert math.ceil(math.log(1 - seuil) / math.log(q)) == attendu
    # leurres : attendu - 1 (insuffisant), attendu + 1 (pas le plus petit), seuil / p arrondi (raisonnement faux)
    assert 1 - q ** (attendu - 1) < seuil
    leurre_naif = math.ceil(seuil / p)
    assert leurre_naif != attendu and 1 - q**leurre_naif < seuil
# Valeurs citées dans la correction (exercice) : 1 - 0,96^73 ≈ 0,9491 et 1 - 0,96^74 ≈ 0,9512.
assert round(1 - 0.96**73, 4) == 0.9492 and round(1 - 0.96**74, 4) == 0.9512
assert round(1 - 0.95**89, 4) == 0.9896 and round(1 - 0.95**90, 4) == 0.9901  # v1
assert round(1 - 0.98**113, 4) == 0.8980 and round(1 - 0.98**114, 4) == 0.9001  # v2
assert round(math.log(0.05) / math.log(0.96), 2) == 73.39
assert round(math.log(0.01) / math.log(0.95), 2) == 89.78
assert round(math.log(0.1) / math.log(0.98), 2) == 113.97

# proba-t-binomiale-pascal-ex5 : C(n+1, 3) - C(n, 3) = C(n, 2) = valeur, n ≥ 3.
for valeur, attendu in ((28, 8), (45, 10), (66, 12)):
    sols = [n for n in range(3, 200) if comb(n + 1, 3) - comb(n, 3) == valeur]
    assert sols == [attendu]
    assert all(comb(n + 1, 3) - comb(n, 3) != valeur for n in (attendu - 1, attendu + 1, valeur))
assert comb(9, 3) - comb(8, 3) == 28 and comb(9, 3) == 84 and comb(8, 3) == 56
assert comb(11, 3) == 165 and comb(10, 3) == 120 and comb(13, 3) == 286 and comb(12, 3) == 220
assert 1 + 4 * 56 == 15**2 and 1 + 4 * 90 == 19**2 and 1 + 4 * 132 == 23**2


# proba-t-binomiale-produit-ex5 : coïncidences (m valeurs équiprobables, k personnes), énumération exacte.
def coincidence(m, k):
    issues = list(product(range(m), repeat=k))
    return F(sum(1 for t in issues if len(set(t)) < k), len(issues))


assert 12**5 == 248832 and perm(12, 5) == 95040 and F(95040, 248832) == F(55, 144)
assert coincidence(12, 5) == F(89, 144) and round(89 / 144, 3) == 0.618
assert coincidence(12, 4) == F(41, 96) and F(41, 96) < F(1, 2) and perm(12, 4) == 11880 and 12**4 == 20736
assert all(v != F(89, 144) for v in (F(55, 144), F(5, 6), F(41, 96)))
# v1 : jours de la semaine, 4 personnes ; question 4 avec 3 personnes.
assert 7**4 == 2401 and perm(7, 4) == 840 and F(840, 2401) == F(120, 343)
assert coincidence(7, 4) == F(223, 343) and round(223 / 343, 3) == 0.650
assert coincidence(7, 3) == F(19, 49) and F(19, 49) < F(1, 2) and perm(7, 3) == 210
assert all(v != F(223, 343) for v in (F(120, 343), F(6, 7), F(19, 49)))
# v2 : chiffres de 0 à 9, 4 personnes ; question 4 avec 5 personnes.
assert perm(10, 4) == 5040 and F(5040, 10000) == F(63, 125)
assert coincidence(10, 4) == F(62, 125) and F(62, 125) < F(1, 2)
assert coincidence(10, 5) == F(436, 625) and F(436, 625) > F(1, 2) and perm(10, 5) == 30240
assert all(v != F(62, 125) for v in (F(63, 125), F(3, 5), F(436, 625)))

print("enrichi_proba_t_binomiale : tout est vérifié.")
