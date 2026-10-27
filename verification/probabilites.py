# Vérification calculatoire du thème « probabilites » (cartes et exercices, niveaux terminale et l1).
# Toutes les probabilités sont calculées en fractions exactes, le plus souvent par énumération complète
# d'un univers fini équiprobable. Chaque assert cite l'identifiant de la carte ou de l'exercice concerné.
from fractions import Fraction as F
from itertools import product, permutations, combinations
from math import comb
import sympy as sp


def loi(valeurs_probas):
    """Regroupe une liste de couples (valeur, probabilité) en dictionnaire valeur -> probabilité."""
    d = {}
    for v, q in valeurs_probas:
        d[v] = d.get(v, 0) + q
    return d


def esperance(d):
    return sum(F(v) * q for v, q in d.items())


def variance(d):
    m = esperance(d)
    return sum((F(v) - m) ** 2 * q for v, q in d.items())


def binom_pmf(n, p, k):
    return comb(n, k) * p ** k * (1 - p) ** (n - k)


# ------------------------------------------------------------------------------------------
# proba-t-binomiale (carte)
# ------------------------------------------------------------------------------------------
# Formule générale : pour n petits et p rationnels, la loi du nombre de succès obtenue par énumération
# des issues de {0,1}^n (produit des probabilités, épreuves indépendantes) coïncide avec B(n,p),
# l'espérance vaut np et la variance np(1-p).
for n in range(1, 6):
    for p in [F(0), F(1, 3), F(1, 2), F(3, 4), F(1)]:
        d = loi([(sum(w), sp.prod([p if x else 1 - p for x in w])) for w in product([0, 1], repeat=n)])
        for k in range(n + 1):
            assert d.get(k, 0) == binom_pmf(n, p, k), "proba-t-binomiale : loi binomiale"
        assert esperance(d) == n * p, "proba-t-binomiale : espérance np"
        assert variance(d) == n * p * (1 - p), "proba-t-binomiale : variance np(1-p)"
        # Identités utilisées dans la démonstration : k C(n,k) = n C(n-1,k-1), k(k-1) C(n,k) = n(n-1) C(n-2,k-2)
        for k in range(1, n + 1):
            assert k * comb(n, k) == n * comb(n - 1, k - 1)
        for k in range(2, n + 1):
            assert k * (k - 1) * comb(n, k) == n * (n - 1) * comb(n - 2, k - 2)
        assert sum(F(k * (k - 1)) * binom_pmf(n, p, k) for k in range(n + 1)) == n * (n - 1) * p ** 2

# Contre-exemple « sans remise » : urne de 3 rouges et 2 blanches, 2 tirages successifs sans remise.
boules = ['R1', 'R2', 'R3', 'B1', 'B2']
tirages = list(permutations(boules, 2))  # 20 issues équiprobables
assert len(tirages) == 20
p_rouge_1 = F(sum(1 for t in tirages if t[0][0] == 'R'), len(tirages))
p_rouge_2 = F(sum(1 for t in tirages if t[1][0] == 'R'), len(tirages))
assert p_rouge_1 == F(3, 5) and p_rouge_2 == F(3, 5), "proba-t-binomiale : chaque tirage donne rouge avec proba 3/5"
p_deux_rouges = F(sum(1 for t in tirages if t[0][0] == 'R' and t[1][0] == 'R'), len(tirages))
assert p_deux_rouges == F(3, 10), "proba-t-binomiale : P(X=2)=3/10 sans remise"
assert binom_pmf(2, F(3, 5), 2) == F(9, 25) != F(3, 10), "proba-t-binomiale : B(2,3/5) donnerait 9/25"
assert p_deux_rouges != p_rouge_1 * p_rouge_2, "proba-t-binomiale : les deux tirages ne sont pas indépendants"

# Contre-exemple « paramètres différents » : épreuves indépendantes de paramètres 1/2 et 1/3.
d = loi([(a + b, (F(1, 2) if a else F(1, 2)) * (F(1, 3) if b else F(2, 3))) for a in [0, 1] for b in [0, 1]])
assert d == {0: F(1, 3), 1: F(1, 2), 2: F(1, 6)}, "proba-t-binomiale : loi pour paramètres 1/2 et 1/3"
q = sp.symbols('q', positive=True)
# Si X suivait B(2,q) : q^2 = 1/6 impose q = 1/sqrt(6), mais alors (1-q)^2 != 1/3.
sol = sp.solve(sp.Eq(q ** 2, sp.Rational(1, 6)), q)
assert sol == [1 / sp.sqrt(6)]
assert sp.simplify((1 - sol[0]) ** 2 - sp.Rational(1, 3)) != 0, "proba-t-binomiale : aucune loi B(2,q) ne convient"
assert esperance(d) == F(1, 2) + F(1, 3), "proba-t-binomiale : l'espérance reste la somme des paramètres"

# Contre-exemple « épreuves non indépendantes » : la même épreuve comptée deux fois, X = 2 X_1.
for p in [F(1, 4), F(1, 2), F(2, 3)]:
    d = loi([(0, 1 - p), (2, p)])
    assert d.get(1, 0) == 0 != binom_pmf(2, p, 1), "proba-t-binomiale : X=2X_1 n'est pas binomiale"

# Affirmation : X ~ B(1,1/2) et Y = X donnent X+Y = 2X, P(X+Y=1) = 0 alors que B(2,1/2) donne 1/2.
assert binom_pmf(2, F(1, 2), 1) == F(1, 2), "proba-t-binomiale : affirmation sur X+Y"
# Affirmation vraie : n - X suit B(n, 1-p).
for n in range(1, 7):
    for p in [F(1, 5), F(2, 3)]:
        for k in range(n + 1):
            assert binom_pmf(n, p, n - k) == binom_pmf(n, 1 - p, k), "proba-t-binomiale : n-X ~ B(n,1-p)"

# Exemple de la carte : QCM, X ~ B(10,1/4).
assert 10 * F(1, 4) == F(5, 2) and 10 * F(1, 4) * F(3, 4) == F(15, 8), "proba-t-binomiale : exemple E et V"
assert binom_pmf(10, F(1, 4), 0) == F(3, 4) ** 10 == F(59049, 1048576), "proba-t-binomiale : exemple P(X=0)"

# ------------------------------------------------------------------------------------------
# proba-t-binomiale-ex1 : X ~ B(4,1/3), P(X=2) = 8/27
# ------------------------------------------------------------------------------------------
assert binom_pmf(4, F(1, 3), 2) == F(8, 27)
for leurre in [F(4, 81), F(16, 81), F(4, 9)]:
    assert leurre != F(8, 27), "proba-t-binomiale-ex1 : leurre faux"
assert F(4, 81) == F(1, 3) ** 2 * F(2, 3) ** 2  # oubli du coefficient binomial
assert F(16, 81) == binom_pmf(4, F(1, 3), 0)    # confusion avec P(X=0)
assert 4 * F(1, 3) == F(4, 3) and 4 * F(1, 3) * F(2, 3) == F(8, 9), "proba-t-binomiale-ex1 : E et V"

# ------------------------------------------------------------------------------------------
# proba-t-binomiale-ex2 : plus petit n tel que 1 - (4/5)^n >= 9/10 : n = 11
# ------------------------------------------------------------------------------------------
n_min = next(n for n in range(1, 100) if 1 - F(4, 5) ** n >= F(9, 10))
assert n_min == 11, "proba-t-binomiale-ex2"
assert F(4, 5) ** 10 > F(1, 10) and F(4, 5) ** 11 <= F(1, 10)
assert 10 < sp.N(sp.log(10) / sp.log(sp.Rational(5, 4))) < 11  # ln 10 / ln(5/4) environ 10,32
for leurre in [10, 12, 5]:
    assert leurre != 11

# ------------------------------------------------------------------------------------------
# proba-t-binomiale-ex3 : urne 4 rouges, 2 blanches, 3 tirages successifs sans remise
# ------------------------------------------------------------------------------------------
boules = ['R1', 'R2', 'R3', 'R4', 'B1', 'B2']
tirages = list(permutations(boules, 3))
assert len(tirages) == 120
for i in range(3):
    assert F(sum(1 for t in tirages if t[i][0] == 'R'), 120) == F(2, 3), "proba-t-binomiale-ex3 : rouge au tirage i"
    assert sum(1 for t in tirages if t[i][0] == 'R') == 4 * 5 * 4
d = loi([(sum(1 for b in t if b[0] == 'R'), F(1, 120)) for t in tirages])
assert d == {1: F(1, 5), 2: F(3, 5), 3: F(1, 5)}, "proba-t-binomiale-ex3 : loi de X"
for k in range(4):
    perm = lambda m, r: sp.ff(m, r)  # nombre d'arrangements
    assert F(comb(3, k) * int(perm(4, k)) * int(perm(2, 3 - k)), 120) == d.get(k, 0)
    assert F(comb(4, k) * comb(2, 3 - k), comb(6, 3)) == d.get(k, 0)
assert esperance(d) == 2 == 3 * F(2, 3), "proba-t-binomiale-ex3 : même espérance que B(3,2/3)"
assert sum(F(v) ** 2 * q for v, q in d.items()) == F(22, 5)
assert variance(d) == F(2, 5), "proba-t-binomiale-ex3 : V(X) = 2/5"
assert 3 * F(2, 3) * F(1, 3) == F(2, 3) != F(2, 5), "proba-t-binomiale-ex3 : variance binomiale 2/3 différente"
assert binom_pmf(3, F(2, 3), 0) == F(1, 27) != 0
# Dépendance : P(R1 et R2) != P(R1) P(R2)
assert F(sum(1 for t in tirages if t[0][0] == 'R' and t[1][0] == 'R'), 120) == F(2, 5) != F(4, 9)
for leurre in [F(2, 3), F(2), F(4, 9)]:
    assert leurre != F(2, 5)

# ------------------------------------------------------------------------------------------
# proba-t-bienayme (carte)
# ------------------------------------------------------------------------------------------
def p_ecart(d, a, stricte=False):
    m = esperance(d)
    return sum(q for v, q in d.items() if (abs(F(v) - m) > a if stricte else abs(F(v) - m) >= a))

# Vérification de l'inégalité sur une famille de lois et de seuils a > 0.
lois_test = [
    {-1: F(1, 2), 1: F(1, 2)},
    {-1: F(1, 8), 0: F(3, 4), 1: F(1, 8)},
    {0: F(1, 3), 1: F(1, 6), 5: F(1, 2)},
    {k: F(1, 9) for k in range(1, 10)},
    {7: F(1)},
]
for d in lois_test:
    for a in [F(1, 10), F(1, 2), F(1), F(3, 2), F(4), F(7)]:
        assert p_ecart(d, a) <= variance(d) / a ** 2, "proba-t-bienayme : inégalité"
        assert p_ecart(d, a, stricte=True) <= variance(d) / a ** 2, "proba-t-bienayme : version stricte (affirmation vraie)"
        assert 1 - p_ecart(d, a) >= 1 - variance(d) / a ** 2, "proba-t-bienayme : passage au complémentaire"

# Cas d'égalité : X = -1, 0, 1 avec probabilités 1/8, 3/4, 1/8 et a = 1.
d_eg = {-1: F(1, 8), 0: F(3, 4), 1: F(1, 8)}
assert esperance(d_eg) == 0 and variance(d_eg) == F(1, 4)
assert p_ecart(d_eg, F(1)) == F(1, 4) == variance(d_eg) / 1, "proba-t-bienayme : égalité possible"

# Confusion avec Markov : X = -1 ou 1 avec probabilité 1/2 ; E(X) = 0 mais P(X >= 1) = 1/2.
d_pm = {-1: F(1, 2), 1: F(1, 2)}
assert esperance(d_pm) == 0
assert sum(q for v, q in d_pm.items() if v >= 1) == F(1, 2) > esperance(d_pm) / 1, "proba-t-bienayme : Markov sans positivité"

# Leurres du texte à trous : V/a, sigma/a^2 et E/a ne sont pas des majorants valables.
c = F(1, 2)
d_c = {-c: F(1, 2), c: F(1, 2)}
assert variance(d_c) == c ** 2 and p_ecart(d_c, c) == 1 > variance(d_c) / c, "proba-t-bienayme : leurre V/a"
c = 2
d_c = {-c: F(1, 2), c: F(1, 2)}
sigma = sp.sqrt(variance(d_c))
assert sigma == 2 and p_ecart(d_c, c) == 1 > sigma / c ** 2, "proba-t-bienayme : leurre sigma/a^2"
assert p_ecart(d_pm, 1) == 1 > esperance(d_pm) / 1, "proba-t-bienayme : leurre E/a"
# Leurres sur a : a = 0 n'a pas de sens ; a < 0 rend l'inégalité fausse (P = 1 > V/a^2 pour |a| grand).
assert sum(d_pm.values()) == 1 > variance(d_pm) / F(-3) ** 2, "proba-t-bienayme : a négatif"

# ------------------------------------------------------------------------------------------
# proba-t-bienayme-ex1 : E = 8, V = 2, a = 4 : borne 1/8
# ------------------------------------------------------------------------------------------
assert F(2) / 4 ** 2 == F(1, 8)
assert F(2) / 4 == F(1, 2) and 1 - F(1, 8) == F(7, 8)
assert sp.sqrt(2) / 16 != sp.Rational(1, 8)

# ------------------------------------------------------------------------------------------
# proba-t-bienayme-ex2 : X ~ B(200,1/4), P(35 < X < 65) >= 5/6
# ------------------------------------------------------------------------------------------
n, p = 200, F(1, 4)
assert n * p == 50 and n * p * (1 - p) == F(75, 2)
borne = 1 - F(75, 2) / 15 ** 2
assert borne == F(5, 6), "proba-t-bienayme-ex2"
exacte = sum(binom_pmf(n, p, k) for k in range(36, 65))
assert exacte >= borne, "proba-t-bienayme-ex2 : la minoration est bien vérifiée par la valeur exacte"
assert all(abs(k - 50) < 15 for k in range(36, 65)) and abs(35 - 50) == 15 and abs(65 - 50) == 15
assert 1 - F(50) / 225 == F(7, 9)            # leurre : V remplacée par np
assert 1 - F(75, 2) / 30 ** 2 == F(23, 24)   # leurre : a = 30 au lieu de 15
for leurre in [F(1, 6), F(7, 9), F(23, 24)]:
    assert leurre != F(5, 6)

# ------------------------------------------------------------------------------------------
# proba-t-bienayme-ex3 : X uniforme sur {1,...,9}
# ------------------------------------------------------------------------------------------
d_u = {k: F(1, 9) for k in range(1, 10)}
assert esperance(d_u) == 5 and variance(d_u) == F(20, 3) == F(9 ** 2 - 1, 12)
assert sum(F(k) ** 2 for k in range(1, 10)) / 9 == F(95, 3)
assert variance(d_u) / 16 == F(5, 12), "proba-t-bienayme-ex3 : borne 5/12"
assert p_ecart(d_u, 4) == F(2, 9), "proba-t-bienayme-ex3 : valeur exacte 2/9"
assert F(2, 9) < F(5, 12)
assert variance(d_u) / 4 == F(5, 3) and variance(d_u) / 64 == F(5, 48)
for leurre in [F(2, 9), F(5, 3), F(5, 48)]:
    assert leurre != F(5, 12)
# Quel que soit a dans {1,...,4}, comparaison exacte / borne (utilisée dans la correction : a = 2 donne 5/3 > 1).
assert p_ecart(d_u, 2) == F(6, 9) and variance(d_u) / 4 == F(5, 3)

# ------------------------------------------------------------------------------------------
# proba-t-bienayme-ex4 : cas d'égalité. Égalité si et seulement si P(X dans {m-a, m, m+a}) = 1.
# Vérification exhaustive sur des lois à 3 valeurs, puis sur la construction proposée.
# ------------------------------------------------------------------------------------------
valeurs = [F(-2), F(-1), F(0), F(1, 2), F(1), F(2), F(3)]
poids = [F(1, 6), F(1, 3), F(1, 2)]
for vs in combinations(valeurs, 3):
    for ps in set(permutations(poids)):
        d = dict(zip(vs, ps))
        m = esperance(d)
        for a in [F(1, 2), F(1), F(3, 2), F(2), F(5, 2)]:
            egalite = p_ecart(d, a) == variance(d) / a ** 2
            support = all(v in (m - a, m, m + a) for v in d)
            assert egalite == support, "proba-t-bienayme-ex4 : caractérisation de l'égalité"
# Construction : X = m - a, m, m + a avec probabilités v/(2a^2), 1 - v/a^2, v/(2a^2), pour 0 < v <= a^2.
for m in [F(0), F(3)]:
    for a in [F(1), F(2)]:
        for v in [a ** 2 / 4, a ** 2 / 2, a ** 2]:
            d = loi([(m - a, v / (2 * a ** 2)), (m, 1 - v / a ** 2), (m + a, v / (2 * a ** 2))])
            assert esperance(d) == m and variance(d) == v and p_ecart(d, a) == v / a ** 2, "proba-t-bienayme-ex4 : construction"
            assert all(q >= 0 for q in d.values()) and sum(d.values()) == 1

# ------------------------------------------------------------------------------------------
# proba-t-lgn (carte)
# ------------------------------------------------------------------------------------------
def loi_moyenne(d, n):
    """Loi de M_n pour un échantillon de taille n de la loi d (univers produit fini)."""
    return loi([(sum(F(x) for x in w) / n, sp.prod([d[x] for x in w])) for w in product(list(d), repeat=n)])

for d in [{0: F(1, 2), 1: F(1, 2)}, {0: F(2, 3), 1: F(1, 3)}, {-1: F(1, 4), 0: F(1, 4), 2: F(1, 2)}]:
    mu, V = esperance(d), variance(d)
    for n in range(1, 6):
        dm = loi_moyenne(d, n)
        assert esperance(dm) == mu and variance(dm) == V / n, "proba-t-lgn : E(M_n) = mu, V(M_n) = V/n"
        for delta in [F(1, 4), F(1, 2), F(1)]:
            assert p_ecart(dm, delta) <= V / (n * delta ** 2), "proba-t-lgn : inégalité de concentration"

# Leurres du texte à trous : X Bernoulli(1/2), n = 2, delta = 1/2 : P = 1/2 > V/(n^2 delta^2) = 1/4.
dB = {0: F(1, 2), 1: F(1, 2)}
dm2 = loi_moyenne(dB, 2)
assert dm2 == {0: F(1, 4), F(1, 2): F(1, 2), 1: F(1, 4)}
assert p_ecart(dm2, F(1, 2)) == F(1, 2) > variance(dB) / (4 * F(1, 4)), "proba-t-lgn : leurre V/(n^2 delta^2)"
# Leurres V/(n delta) et sigma/(n delta^2) : n = 1, X = -c ou c.
for c, test in [(F(1, 2), 'V'), (F(2), 'sigma')]:
    d_c = {-c: F(1, 2), c: F(1, 2)}
    if test == 'V':
        assert p_ecart(d_c, c) == 1 > variance(d_c) / c, "proba-t-lgn : leurre V/(n delta)"
    else:
        assert p_ecart(d_c, c) == 1 > sp.sqrt(variance(d_c)) / c ** 2, "proba-t-lgn : leurre sigma/(n delta^2)"
# L'inégalité sans le facteur n (V/delta^2) est vraie (plus faible) : ce n'est donc PAS un leurre utilisé.

# Contre-exemple : X_i = X_1 pour tout i, X_1 Bernoulli(1/2), delta = 1/4.
for n in range(1, 30):
    # M_n = X_1 vaut 0 ou 1 : |M_n - 1/2| = 1/2 >= 1/4 toujours.
    assert p_ecart(dB, F(1, 4)) == 1, "proba-t-lgn : X_i = X_1"
    borne = variance(dB) / (n * F(1, 4) ** 2)
    assert borne == F(4, n)
    if n >= 5:
        assert 1 > borne, "proba-t-lgn : l'inégalité tombe en défaut pour n >= 5"
assert all(F(4, n) >= 1 for n in range(1, 5))

# Affirmation fausse « M_n = mu pour n grand » : Bernoulli(1/2), n impair, P(M_n = 1/2) = 0.
for n in [1, 3, 5]:
    assert loi_moyenne(dB, n).get(F(1, 2), 0) == 0, "proba-t-lgn : M_n ne vaut jamais 1/2 pour n impair"

# Remarque de la carte : k événements mutuellement indépendants de probabilités dans ]0,1[ exigent 2^k atomes.
# (Argument : les 2^k intersections d'événements ou de complémentaires sont disjointes et de probabilité > 0.)
for k in range(1, 5):
    probs = [F(1, i + 2) for i in range(k)]
    atomes = [sp.prod([pr if s else 1 - pr for pr, s in zip(probs, signes)]) for signes in product([0, 1], repeat=k)]
    assert len(atomes) == 2 ** k and all(x > 0 for x in atomes) and sum(atomes) == 1

# ------------------------------------------------------------------------------------------
# proba-t-lgn-ex1 : V = 4, n = 100, delta = 1/2 : borne 4/25
# ------------------------------------------------------------------------------------------
assert F(4) / (100 * F(1, 2) ** 2) == F(4, 25)
assert F(4) / (100 * F(1, 2)) == F(2, 25)
assert F(4) / (100 ** 2 * F(1, 2) ** 2) == F(1, 625)
assert 1 - F(4, 25) == F(21, 25)
for leurre in [F(2, 25), F(1, 625), F(21, 25)]:
    assert leurre != F(4, 25)

# ------------------------------------------------------------------------------------------
# proba-t-lgn-ex2 : dé équilibré, fréquence des 6, delta = 1/10, risque 1/20 : n >= 2500/9, n = 278
# ------------------------------------------------------------------------------------------
V6 = F(1, 6) * F(5, 6)
assert V6 == F(5, 36)
nm = next(n for n in range(1, 10 ** 4) if V6 / (n * F(1, 10) ** 2) <= F(1, 20))
assert nm == 278 and F(2500, 9) > 277, "proba-t-lgn-ex2"
assert next(n for n in range(1, 10 ** 4) if V6 / (n * F(1, 10)) <= F(1, 20)) == 28       # leurre delta non élevé au carré
assert next(n for n in range(1, 10 ** 4) if F(1, 6) / (n * F(1, 10) ** 2) <= F(1, 20)) == 334  # leurre V = p
for leurre in [28, 334, 277]:
    assert leurre != 278

# ------------------------------------------------------------------------------------------
# proba-t-lgn-ex3 : p(1-p) <= 1/4 ; delta = 1/20, risque 1/10 : n >= 1000
# ------------------------------------------------------------------------------------------
pp = sp.symbols('p', real=True)
assert sp.expand(sp.Rational(1, 4) - pp * (1 - pp) - (pp - sp.Rational(1, 2)) ** 2) == 0, "proba-t-lgn-ex3 : identité"
nm = next(n for n in range(1, 10 ** 5) if F(1, 4) / (n * F(1, 20) ** 2) <= F(1, 10))
assert nm == 1000, "proba-t-lgn-ex3"
assert next(n for n in range(1, 10 ** 5) if F(1, 4) / (n * F(1, 20)) <= F(1, 10)) == 50
assert next(n for n in range(1, 10 ** 5) if F(1, 2) / (n * F(1, 20) ** 2) <= F(1, 10)) == 2000
assert next(n for n in range(1, 10 ** 5) if F(1) / (n * F(1, 20) ** 2) <= F(1, 10)) == 4000
# Pour n = 1000 la borne est atteinte pour p = 1/2 : 1/(4 n delta^2) = 1/10.
assert F(1, 4) / (1000 * F(1, 20) ** 2) == F(1, 10)

# ------------------------------------------------------------------------------------------
# proba-t-lgn-ex4 : X_1 Bernoulli(1/2), Omega = {0,1} équiprobable
# ------------------------------------------------------------------------------------------
Omega = [0, 1]
for n in range(1, 12):
    # (1) X_i = X_1 : M_n = X_1
    P1 = F(sum(1 for w in Omega if abs(F(w) - F(1, 2)) >= F(1, 4)), 2)
    assert P1 == 1, "proba-t-lgn-ex4 (1)"
for m in range(1, 6):
    # (2) n = 2m, X_{2k-1} = X_1, X_{2k} = 1 - X_1 : M_n = 1/2 identiquement
    for w in Omega:
        xs = [w if i % 2 == 0 else 1 - w for i in range(2 * m)]
        assert F(sum(xs), 2 * m) == F(1, 2), "proba-t-lgn-ex4 (2)"
# X_2 = 1 - X_1 a la loi de Bernoulli(1/2), mais P(X_1 = 1, X_2 = 1) = 0 != 1/4.
assert F(sum(1 for w in Omega if 1 - w == 1), 2) == F(1, 2)
assert F(sum(1 for w in Omega if w == 1 and 1 - w == 1), 2) == 0 != F(1, 4)

# ------------------------------------------------------------------------------------------
# proba-l1-probas-totales-bayes (carte)
# ------------------------------------------------------------------------------------------
# Test médical : P(M) = 1/100, P(T|M) = 99/100, P(T|non M) = 1/20.
PM, PT_M, PT_nM = F(1, 100), F(99, 100), F(1, 20)
PT = PT_M * PM + PT_nM * (1 - PM)
assert PT == F(594, 10000) == F(297, 5000), "proba-l1-probas-totales-bayes : P(T)"
PM_T = PT_M * PM / PT
assert PM_T == F(1, 6), "proba-l1-probas-totales-bayes : P(M|T) = 1/6"
assert PM_T != PT_M, "proba-l1-probas-totales-bayes : P(M|T) != P(T|M)"
# Leurres du 3e trou : ils donnent des valeurs fausses sur cet exemple.
assert PT_M != PM_T
assert PT_M * PT / PM != PM_T
assert PM / PT != PM_T
# Leurres du 2e trou : avec B = Omega et le système (M, non M).
assert PM ** 2 + (1 - PM) ** 2 != 1  # somme des P(A_i|B)P(A_i) avec B = Omega
assert 1 + 1 != 1                     # somme des P(B|A_i) avec B = Omega
assert PT_M * PM * PT_nM * (1 - PM) != PT  # produit

# Hypothèses nécessaires : dé équilibré, Omega = {1,...,6}.
def P(E):
    return F(len(E), 6)
B = {2, 5}
A1, A2 = {1, 2}, {3, 4}       # incompatibles, réunion différente de Omega
assert P(B & A1) + P(B & A2) == F(1, 6) != P(B) == F(1, 3), "proba-l1-probas-totales-bayes : réunion != Omega"
B = {3}
A1, A2 = {1, 2, 3, 4}, {3, 4, 5, 6}  # réunion Omega, non incompatibles
assert P(B & A1) + P(B & A2) == F(1, 3) != P(B) == F(1, 6), "proba-l1-probas-totales-bayes : non incompatibles"
# Un système complet de deux événements de probabilités non nulles n'est jamais formé d'événements indépendants.
A = {1, 2, 3}
Ac = {4, 5, 6}
assert P(A & Ac) == 0 != P(A) * P(Ac) == F(1, 4), "proba-l1-probas-totales-bayes : A et son contraire non indépendants"
# Exemple de la carte : pièce équilibrée ou pièce truquée (P(pile) = 3/4), choisie avec probabilité 1/2.
Ppile = F(1, 2) * F(1, 2) + F(1, 2) * F(3, 4)
assert Ppile == F(5, 8) and F(1, 2) * F(3, 4) / Ppile == F(3, 5), "proba-l1-probas-totales-bayes : exemple"
# Affirmation fausse : somme des P(B|A_i) = 1. Avec B = Omega et deux A_i, la somme vaut 2.
# Affirmation vraie (utilisée dans la justification) : P(A|B) P(B) = P(B|A) P(A) = P(A inter B).

# ------------------------------------------------------------------------------------------
# proba-l1-probas-totales-bayes-ex1 : urnes
# ------------------------------------------------------------------------------------------
PR = F(1, 3) * F(1, 2) + F(2, 3) * F(1, 4)
assert PR == F(1, 3), "proba-l1-probas-totales-bayes-ex1"
assert (F(1, 2) + F(1, 4)) / 2 == F(3, 8) and F(1, 2) + F(1, 4) == F(3, 4)
for leurre in [F(3, 8), F(3, 4), F(1, 6)]:
    assert leurre != F(1, 3)

# ------------------------------------------------------------------------------------------
# proba-l1-probas-totales-bayes-ex2 : trois machines
# ------------------------------------------------------------------------------------------
pm = [F(1, 2), F(3, 10), F(1, 5)]
pd = [F(1, 100), F(1, 50), F(1, 20)]
assert sum(pm) == 1
PD = sum(a * b for a, b in zip(pm, pd))
assert PD == F(21, 1000), "proba-l1-probas-totales-bayes-ex2 : P(D)"
assert pm[2] * pd[2] / PD == F(10, 21), "proba-l1-probas-totales-bayes-ex2 : P(M3|D)"
assert [a * b / PD for a, b in zip(pm, pd)] == [F(5, 21), F(6, 21), F(10, 21)]
for leurre in [F(1, 5), F(1, 20), F(1, 100)]:
    assert leurre != F(10, 21)

# ------------------------------------------------------------------------------------------
# proba-l1-probas-totales-bayes-ex3 : test de dépistage
# ------------------------------------------------------------------------------------------
PM, PT_M, PT_nM = F(1, 50), F(9, 10), F(1, 20)
PT = PT_M * PM + PT_nM * (1 - PM)
assert PT == F(67, 1000), "proba-l1-probas-totales-bayes-ex3 : P(T)"
assert PT_M * PM == F(18, 1000)
assert PT_M * PM / PT == F(18, 67), "proba-l1-probas-totales-bayes-ex3 : P(M|T)"
assert 1 - PT == F(933, 1000)
assert (1 - PT_nM) * (1 - PM) / (1 - PT) == F(931, 933), "proba-l1-probas-totales-bayes-ex3 : P(non M | non T)"
assert (1 - PT_M) * PM / (1 - PT) == F(2, 933), "proba-l1-probas-totales-bayes-ex3 : P(M | non T)"
assert F(931, 933) + F(2, 933) == 1
for leurre in [F(9, 10), F(9, 500), F(67, 1000)]:
    assert leurre != F(18, 67)
assert F(18, 67) < F(3, 10)  # environ 0,27

# ------------------------------------------------------------------------------------------
# proba-l1-probas-totales-bayes-ex4 : P(A|B) = P(B|A) et P(A inter B) > 0 impose P(A) = P(B)
# Vérification exhaustive sur un univers à 4 issues équiprobables, puis application à ex3.
# ------------------------------------------------------------------------------------------
Om = range(4)
parties = [set(s) for r in range(1, 5) for s in combinations(Om, r)]
Q = lambda E: F(len(E), 4)
for A in parties:
    for B in parties:
        pa, pb, pab = Q(A), Q(B), Q(A & B)
        # (a) P(A|B) > P(A) <=> P(B|A) > P(B)
        assert (pab / pb > pa) == (pab / pa > pb), "proba-l1-probas-totales-bayes-ex4 (a)"
        # (b)
        if pab > 0 and pab / pb == pab / pa:
            assert pa == pb, "proba-l1-probas-totales-bayes-ex4 (b)"
# (c) A, B incompatibles de probabilités différentes : P(A|B) = P(B|A) = 0 mais P(A) != P(B).
A, B = {0}, {1, 2}
assert Q(A & B) == 0 and Q(A) != Q(B)
# (d) Dans ex3, P(M inter T) > 0 et P(M) = 1/50 != 67/1000 = P(T).
assert F(1, 50) != F(67, 1000) and F(18, 1000) > 0
assert F(18, 67) != F(9, 10)

# ------------------------------------------------------------------------------------------
# proba-l1-independance-mutuelle (carte)
# ------------------------------------------------------------------------------------------
def mutuellement_indep(evts, Pr):
    for r in range(2, len(evts) + 1):
        for I in combinations(range(len(evts)), r):
            inter = set.intersection(*[evts[i] for i in I])
            if Pr(inter) != sp.prod([Pr(evts[i]) for i in I]):
                return False
    return True

# Deux lancers d'une pièce équilibrée.
Om2 = set(product('PF', repeat=2))
Pr2 = lambda E: F(len(E), 4)
A = {w for w in Om2 if w[0] == 'P'}
B = {w for w in Om2 if w[1] == 'P'}
C = {w for w in Om2 if w[0] == w[1]}
assert Pr2(A) == Pr2(B) == Pr2(C) == F(1, 2)
for X, Y in [(A, B), (A, C), (B, C)]:
    assert Pr2(X & Y) == F(1, 4) == Pr2(X) * Pr2(Y), "proba-l1-independance-mutuelle : indépendance deux à deux"
assert A & B & C == {('P', 'P')}
assert Pr2(A & B & C) == F(1, 4) != F(1, 8), "proba-l1-independance-mutuelle : pas d'indépendance mutuelle"
assert not mutuellement_indep([A, B, C], Pr2)
assert Pr2(A & (B & C)) == F(1, 4) != Pr2(A) * Pr2(B & C) == F(1, 8), "proba-l1-independance-mutuelle : A et B inter C"

# P(A inter B inter C) = P(A)P(B)P(C) seul ne suffit pas : Omega = {1,...,8}.
Pr8 = lambda E: F(len(E), 8)
A, B, C = {1, 2, 3, 4}, {1, 5, 6, 7}, {1, 5, 6, 8}
assert Pr8(A) == Pr8(B) == Pr8(C) == F(1, 2)
assert Pr8(A & B & C) == F(1, 8) == Pr8(A) * Pr8(B) * Pr8(C), "proba-l1-independance-mutuelle : produit triple"
assert Pr8(A & B) == F(1, 8) != F(1, 4), "proba-l1-independance-mutuelle : A, B non indépendants"
assert Pr8(B & C) == F(3, 8) and Pr8(A & C) == F(1, 8)
# Cas dégénéré cité : C = vide, A = B de probabilité 1/2.
A = B = {1, 2, 3, 4}
assert Pr8(A & B & set()) == 0 == Pr8(A) * Pr8(B) * 0 and Pr8(A & B) == F(1, 2) != F(1, 4)

# Exemple : trois lancers, A_i = pile au i-ème lancer : mutuellement indépendants.
Om3 = set(product('PF', repeat=3))
Pr3 = lambda E: F(len(E), 8)
Ai = [{w for w in Om3 if w[i] == 'P'} for i in range(3)]
assert mutuellement_indep(Ai, Pr3), "proba-l1-independance-mutuelle : exemple"
assert Pr3(Ai[0] & Ai[2]) == F(1, 4) and Pr3(Ai[0] & Ai[1] & Ai[2]) == F(1, 8)
# Démo : passer au complémentaire conserve l'indépendance mutuelle (vérifié sur l'exemple).
assert mutuellement_indep([Om3 - Ai[0], Ai[1], Ai[2]], Pr3)
assert mutuellement_indep([Om3 - Ai[0], Om3 - Ai[1], Om3 - Ai[2]], Pr3)
# Affirmation fausse : deux événements incompatibles de probabilités non nulles ne sont pas indépendants.
X, Y = {1}, {2}
assert Pr8(X & Y) == 0 != Pr8(X) * Pr8(Y)

# Nombre d'égalités dans la définition : 2^n - n - 1.
for n in range(2, 8):
    assert sum(comb(n, r) for r in range(2, n + 1)) == 2 ** n - n - 1

# ------------------------------------------------------------------------------------------
# proba-l1-independance-mutuelle-ex1 : n = 4 : 11 égalités
# ------------------------------------------------------------------------------------------
assert 2 ** 4 - 4 - 1 == 11 and comb(4, 2) == 6 and 2 ** 4 - 1 == 15
for leurre in [6, 15, 16]:
    assert leurre != 11

# ------------------------------------------------------------------------------------------
# proba-l1-independance-mutuelle-ex2 : tirs indépendants 1/2, 1/3, 1/4
# ------------------------------------------------------------------------------------------
pts = [F(1, 2), F(1, 3), F(1, 4)]
aucun = sp.prod([1 - x for x in pts])
assert aucun == F(1, 4) and 1 - aucun == F(3, 4), "proba-l1-independance-mutuelle-ex2"
# Vérification par énumération sur l'univers produit.
d = sum(sp.prod([x if s else 1 - x for x, s in zip(pts, w)]) for w in product([0, 1], repeat=3) if any(w))
assert d == F(3, 4)
assert sum(pts) == F(13, 12) and sp.prod(pts) == F(1, 24)
for leurre in [F(13, 12), F(1, 24), F(1, 4)]:
    assert leurre != F(3, 4)

# ------------------------------------------------------------------------------------------
# proba-l1-independance-mutuelle-ex3 : deux dés
# ------------------------------------------------------------------------------------------
Om36 = set(product(range(1, 7), repeat=2))
Pr36 = lambda E: F(len(E), 36)
A = {w for w in Om36 if w[0] % 2 == 0}
B = {w for w in Om36 if w[1] >= 4}
C = {w for w in Om36 if (w[0] + w[1]) % 2 == 1}
assert Pr36(A) == Pr36(B) == Pr36(C) == F(1, 2)
for X, Y in [(A, B), (A, C), (B, C)]:
    assert Pr36(X & Y) == F(1, 4), "proba-l1-independance-mutuelle-ex3 : deux à deux"
assert len(B & C) == 9 and len(A & C) == 9 and len(A & B) == 9
assert A & B & C == {(2, 5), (4, 5), (6, 5)}
assert Pr36(A & B & C) == F(1, 12) != F(1, 8), "proba-l1-independance-mutuelle-ex3 : P(A inter B inter C) = 1/12"
for leurre in [F(1, 8), F(1, 4), F(1, 6)]:
    assert leurre != F(1, 12)

# ------------------------------------------------------------------------------------------
# proba-l1-independance-mutuelle-ex4 : n lancers, A_i = pile au i-ème, B = nombre de piles pair.
# Toute sous-famille de n événements parmi A_1, ..., A_n, B est mutuellement indépendante,
# mais la famille des n + 1 ne l'est pas.
# ------------------------------------------------------------------------------------------
for n in range(2, 6):
    Om = set(product([0, 1], repeat=n))
    Pr = lambda E, N=len(Om): F(len(E), N)
    Ai = [{w for w in Om if w[i] == 1} for i in range(n)]
    B = {w for w in Om if sum(w) % 2 == 0}
    assert Pr(B) == F(1, 2)
    famille = Ai + [B]
    assert not mutuellement_indep(famille, Pr), "proba-l1-independance-mutuelle-ex4 : n + 1 événements"
    for retire in range(n + 1):
        sous = [E for j, E in enumerate(famille) if j != retire]
        assert mutuellement_indep(sous, Pr), "proba-l1-independance-mutuelle-ex4 : n événements"
    tous = set.intersection(*Ai)
    assert Pr(tous & B) == (F(1, 2 ** n) if n % 2 == 0 else 0)
    assert Pr(tous & B) != F(1, 2 ** (n + 1))

# ------------------------------------------------------------------------------------------
# Exemples sur un dé équilibré cités dans les textes (Omega = {1,...,6}, fonction P définie plus haut)
# ------------------------------------------------------------------------------------------
# proba-l1-probas-totales-bayes-ex4 (3) : A = {1}, B = {2,3} incompatibles, P(A) = 1/6 != 1/3 = P(B).
A, B = {1}, {2, 3}
assert P(A & B) == 0 and P(A) == F(1, 6) and P(B) == F(1, 3), "proba-l1-probas-totales-bayes-ex4 (3)"
# proba-l1-probas-totales-bayes-ex4 (4) : 1/50 = 20/1000.
assert F(1, 50) == F(20, 1000)
# proba-l1-probas-totales-bayes (affirmation) : A = {1,2,3} et son contraire.
assert P({1, 2, 3}) * P({4, 5, 6}) == F(1, 4)
# proba-l1-independance-mutuelle (affirmation) : {1} et {2} incompatibles, P(A)P(B) = 1/36.
assert P({1} & {2}) == 0 != P({1}) * P({2}) == F(1, 36), "proba-l1-independance-mutuelle : incompatibles"
# proba-t-bienayme (exemple) : E = 20, V = 4, a = 5.
assert F(4, 25) == F(4) / 5 ** 2 and 1 - F(4, 25) == F(21, 25), "proba-t-bienayme : exemple"
# proba-t-lgn (exemple) : V = 1/4, delta = 1/10 : borne 25/n.
for n in range(1, 50):
    assert F(1, 4) / (n * F(1, 10) ** 2) == F(25, n), "proba-t-lgn : exemple"
# proba-t-lgn (affirmation 1) : X_i = X_1 avec V(X_1) > 0 donne V(M_n) = V(X) != V(X)/n pour n >= 2.
assert variance(dB) == F(1, 4) and all(F(1, 4) != F(1, 4) / n for n in range(2, 10))

print("probabilites.py : toutes les vérifications passent.")
