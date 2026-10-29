"""Vérifications calculatoires du chapitre « proba-l2-discretes » (probabilités, L2) :
cartes, exercices et QCM de content/*/probabilites/l2-proba-l2-discretes.json et content/qcm/proba-l2-discretes.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variantes v1, v2) ou de la question.
Calculs exacts en fractions, séries par SymPy (sommes symboliques simples) ou par sommes partielles à 1e-12.
Lancer : python3 verification/proba_l2_discretes.py (doit finir sans erreur, en quelques secondes).
"""

import math
from fractions import Fraction as F
from itertools import product

from sympy import Rational, binomial, exp, log, oo, pi, simplify, summation, symbols

k = symbols("k", integer=True, positive=True)
C = math.comb


def somme(f, n=4000):
    """Somme partielle numérique sum_{k=1}^{n} f(k)."""
    return sum(f(j) for j in range(1, n + 1))


# ---------------------------------------------------------------------------
# proba-l2-discretes-espace (carte)
# ---------------------------------------------------------------------------
# Contre-exemple : P({k}) = 2^-(k+1) sur N, A_n = {k >= n} : P(A_n) = 2^-n.
assert summation(Rational(1, 2) ** (k + 1), (k, 0, oo)) == 1
assert all(sum(F(1, 2 ** (j + 1)) for j in range(n, 200)) - F(1, 2 ** n) < F(1, 2 ** 150) for n in range(10))
# Exemple : P(A_n) = 1 - 2^-n -> 1.
# Trou 2, leurre « somme » : A_n = Omega pour tout n donne une somme infinie ; leurre inf : A_n = {0..n}, P(A_0) = 1/2 != 1.
assert F(1, 2) != 1

# proba-l2-discretes-espace-ex1 (+ v1, v2) : P({k}) = c b^-k, c = b - 1, P(pair) = 1/(b+1).
for b, rep, leurres in [(3, F(1, 4), [F(1, 2), F(3, 4), F(1, 9)]), (4, F(1, 5), [F(1, 2), F(4, 5), F(1, 16)]),
                        (2, F(1, 3), [F(1, 2), F(2, 3), F(1, 4)])]:
    c = b - 1
    assert summation(c * Rational(1, b) ** k, (k, 1, oo)) == 1
    pair = summation(c * Rational(1, b) ** (2 * k), (k, 1, oo))
    assert pair == Rational(rep.numerator, rep.denominator) == Rational(c, b * b - 1)
    assert F(c, b) == {3: F(2, 3), 4: F(3, 4), 2: F(1, 2)}[b]  # P({1})
    assert rep not in leurres and len(set(leurres)) == 3

# proba-l2-discretes-espace-ex2 (+ v1, v2) : nombre d'atomes et cardinal de la tribu engendrée (force brute).


def tribu_engendree(omega, parties):
    T = {frozenset(), frozenset(omega)}
    T |= {frozenset(p) for p in parties}
    while True:
        nouveau = set(T)
        for a in T:
            nouveau.add(frozenset(omega) - a)
            for b2 in T:
                nouveau.add(a | b2)
        if nouveau == T:
            return T
        T = nouveau


for omega, parties, rep, leurres in [
    (range(1, 7), [{1, 2}, {2, 3}], 16, [4, 8, 64]),
    (range(1, 9), [{1, 2, 3}, {3, 4}, {4, 5, 6}], 32, [8, 16, 256]),
    (range(1, 7), [{1, 2, 3}, {4, 5, 6}], 4, [16, 8, 64]),
]:
    T = tribu_engendree(set(omega), parties)
    assert len(T) == rep and rep not in leurres
t0 = tribu_engendree(set(range(1, 7)), [{1, 2}, {2, 3}])
assert frozenset({4, 5, 6}) in t0 and frozenset({4}) not in t0

# proba-l2-discretes-espace-ex4 (+ v1, v2) : P(T multiple de 3) = p q^2/(1-q^3) ; résidus 1, 2, 0 proportionnels à 1, q, q^2.
for p, rep, trois, leurres in [
    (F(1, 2), F(1, 7), (F(4, 7), F(2, 7), F(1, 7)), [F(1, 3), F(2, 7), F(1, 8)]),
    (F(1, 3), F(4, 19), (F(9, 19), F(6, 19), F(4, 19)), [F(1, 3), F(6, 19), F(4, 27)]),
    (F(2, 3), F(1, 13), (F(9, 13), F(3, 13), F(1, 13)), [F(1, 3), F(3, 13), F(2, 27)]),
]:
    q = 1 - p
    assert p * q * q / (1 - q ** 3) == rep
    assert (p / (1 - q ** 3), p * q / (1 - q ** 3), p * q * q / (1 - q ** 3)) == trois and sum(trois) == 1
    assert leurres[2] == q * q * p  # P(T = 3)
    assert abs(sum(float(q) ** (3 * j - 1) * float(p) for j in range(1, 400)) - float(rep)) < 1e-12
    assert rep not in leurres

# ---------------------------------------------------------------------------
# proba-l2-discretes-repartition (carte et exercices)
# ---------------------------------------------------------------------------
# Exemple : P(X=k) = 2^-k, F(n) = 1 - 2^-n.
assert all(sum(F(1, 2 ** j) for j in range(1, n + 1)) == 1 - F(1, 2 ** n) for n in range(1, 30))
# Piège 1 (carte) : P({a}) = P({b}) = 1 donne une « loi » de somme 2.
# Affirmation 5 : X = 0, P(0 <= X <= 1) = 1 mais F(1) - F(0) = 0.

# proba-l2-discretes-repartition-ex1 (+ v1, v2)
p2 = lambda j: F(1, 2 ** j)
assert sum(p2(j) for j in (3, 4, 5)) == F(7, 32)
assert [sum(p2(j) for j in (2, 3, 4, 5)), sum(p2(j) for j in (3, 4)), 1 - F(1, 32)] == [F(15, 32), F(3, 16), F(31, 32)]
p3 = lambda j: F(2, 3 ** j)
assert sum(p3(j) for j in range(1, 30)) > 1 - F(1, 10 ** 12)
assert p3(2) + p3(3) == F(8, 27) and 1 - F(1, 27) - (1 - F(1, 3)) == F(8, 27)
assert [p3(1) + p3(2) + p3(3), p3(2), p3(3)] == [F(26, 27), F(2, 9), F(2, 27)]
pr = lambda j: F(1, (j + 1) * (j + 2))
Fr = lambda n: sum(pr(j) for j in range(0, n + 1))
assert all(Fr(n) == 1 - F(1, n + 2) for n in range(30))
assert Fr(4) - Fr(1) == F(1, 6) and pr(2) + pr(3) + pr(4) == F(1, 6)
assert [Fr(4) - Fr(0), Fr(3) - Fr(1), Fr(4)] == [F(1, 3), F(2, 15), F(5, 6)]

# proba-l2-discretes-repartition-ex2 (+ v1, v2) : sauts de F, P(X >= a), E(X).
for loi, a, rep, leurres, esp in [
    ({-1: F(2, 10), 0: F(5, 10), 2: F(3, 10)}, 0, F(8, 10), [F(3, 10), F(7, 10), F(5, 10)], F(4, 10)),
    ({0: F(1, 10), 1: F(3, 10), 3: F(5, 10), 5: F(1, 10)}, 3, F(6, 10), [F(1, 10), F(9, 10), F(5, 10)], F(23, 10)),
    ({-2: F(25, 100), 1: F(35, 100), 4: F(40, 100)}, 1, F(75, 100), [F(40, 100), F(60, 100), F(35, 100)], F(145, 100)),
]:
    assert sum(loi.values()) == 1
    assert sum(v for x, v in loi.items() if x >= a) == rep
    assert [sum(v for x, v in loi.items() if x > a), sum(v for x, v in loi.items() if x <= a), loi[a]] == leurres
    assert sum(x * v for x, v in loi.items()) == esp
# Valeurs de F dans les énoncés (cumuls).
assert [F(2, 10), F(7, 10)] == [F(2, 10), F(2, 10) + F(5, 10)]
assert [F(1, 10), F(4, 10), F(9, 10)] == [F(1, 10), F(1, 10) + F(3, 10), F(1, 10) + F(3, 10) + F(5, 10)]
assert [F(25, 100), F(60, 100)] == [F(25, 100), F(25, 100) + F(35, 100)]

# proba-l2-discretes-repartition-ex3 (+ v1, v2) : X de loi 1/(k(k+1)), P(X >= k) = 1/k, E(min(X, m)) = H_m.
px = lambda j: F(1, j * (j + 1))
assert all(1 - sum(px(j) for j in range(1, n)) == F(1, n) for n in range(1, 40))  # P(X >= n)
H = lambda m: sum(F(1, j) for j in range(1, m + 1))
for m, rep, leurres in [(5, F(137, 60), [F(77, 60), F(29, 20), F(25, 12)]),
                        (4, F(25, 12), [F(13, 12), F(77, 60), F(11, 6)]),
                        (6, F(49, 20), [F(29, 20), F(223, 140), F(137, 60)])]:
    loiY = {j: px(j) for j in range(1, m)}
    loiY[m] = F(1, m)
    assert sum(loiY.values()) == 1
    EY = sum(j * v for j, v in loiY.items())
    assert EY == rep == H(m)
    oubli = sum(j * px(j) for j in range(1, m))
    mauvais = oubli + m * px(m)
    assert leurres[0] == oubli and leurres[1] == mauvais and leurres[2] == H(m - 1)
    assert rep not in leurres

# ---------------------------------------------------------------------------
# proba-l2-discretes-lois-usuelles (carte)
# ---------------------------------------------------------------------------
N_, K_, n_ = symbols("N K n", positive=True)
pp = K_ / N_
Vh = K_ * (n_ / N_) * (1 - n_ / N_) + K_ * (K_ - 1) * (n_ * (n_ - 1) / (N_ * (N_ - 1)) - n_ ** 2 / N_ ** 2)
assert simplify(Vh - n_ * pp * (1 - pp) * (N_ - n_) / (N_ - 1)) == 0
# Hypergéométrique : loi, espérance, variance, contrôle exhaustif sur de petites urnes.
for N in range(2, 9):
    for K in range(0, N + 1):
        for n in range(1, N + 1):
            loi = {j: F(C(K, j) * C(N - K, n - j), C(N, n)) for j in range(0, n + 1) if j <= K and n - j <= N - K}
            assert sum(loi.values()) == 1
            E1 = sum(j * v for j, v in loi.items())
            V1 = sum(j * j * v for j, v in loi.items()) - E1 ** 2
            p = F(K, N)
            assert E1 == n * p and V1 == n * p * (1 - p) * F(N - n, N - 1)
            assert V1 <= n * p * (1 - p)  # affirmation 3
# Géométrique : E = 1/p, V = q/p^2, P(X > k) = q^k, absence de mémoire.
for p in [F(1, 2), F(1, 3), F(1, 6), F(1, 4), F(2, 3)]:
    q = 1 - p
    E1 = somme(lambda j: j * float(q) ** (j - 1) * float(p), 3000)
    E2 = somme(lambda j: j * j * float(q) ** (j - 1) * float(p), 3000)
    assert abs(E1 - float(1 / p)) < 1e-9 and abs(E2 - E1 ** 2 - float(q / p ** 2)) < 1e-8
    for kk in range(6):
        for mm in range(6):
            assert q ** (kk + mm) / q ** kk == q ** mm
# Trous : (1-p)^k p de somme 1-p, k p (1-p)^(k-1) de somme 1/p, e^-p p^k/k! (k >= 1) de somme 1 - e^-p : jamais 1 pour 0<p<1.
for p in [0.1, 0.3, 0.5, 0.7, 0.9]:
    assert abs(somme(lambda j: (1 - p) ** j * p, 2000) - (1 - p)) < 1e-9
    assert abs(somme(lambda j: j * p * (1 - p) ** (j - 1), 4000) - 1 / p) < 1e-9
    assert abs(somme(lambda j: math.exp(-p) * p ** j / math.factorial(j), 60) - (1 - math.exp(-p))) < 1e-12
    assert p != 1 / p and (1 - p) / p != 1 / p and 1 / p ** 2 != 1 / p
# Poisson : E = V = lambda ; leurres 2λ, λ + λ², λ(1-λ) différents de λ pour λ > 0.
for l in (0.5, 2, 4.5):
    pk = [math.exp(-l) * l ** j / math.factorial(j) for j in range(120)]
    E1 = sum(j * v for j, v in enumerate(pk))
    assert abs(sum(pk) - 1) < 1e-12 and abs(E1 - l) < 1e-10 and abs(sum(j * j * v for j, v in enumerate(pk)) - E1 ** 2 - l) < 1e-9
    assert len({l, 2 * l, l + l * l, l * (1 - l)}) == 4
# Exemple : C(28,5)/C(32,5) = 1755/3596, E = 5/8 ; e^-2 ≈ 0,135.
assert F(C(28, 5), C(32, 5)) == F(1755, 3596) and C(28, 5) == 98280 and C(32, 5) == 201376
assert abs(1755 / 3596 - 0.488) < 5e-4 and abs(math.exp(-2) - 0.135) < 5e-4
# Contre-exemples : avec remise P(X=2) = 1/4 ; X ∈ {0,2} équiprobable : E = V = 1.
assert F(1, 4) == F(1, 2) ** 2 and F(0 + 2, 2) == 1 and F(0 + 4, 2) - 1 == 1
# Affirmation 2 : λ = 3, P(X=1) = 3 e^-3 > e^-3.

# proba-l2-discretes-lois-usuelles-ex1 (+ v1, v2) : P(X >= 2) = 1 - (1 + λ) e^-λ.
for l, approx in [(2, 0.594), (3, 0.801), (1, 0.264)]:
    v = 1 - (1 + l) * math.exp(-l)
    assert abs(v - approx) < 5e-4
    assert abs(v - (1 - sum(math.exp(-l) * l ** j / math.factorial(j) for j in range(2)))) < 1e-15
    leurres = [(1 + l) * math.exp(-l), 1 - l * math.exp(-l), 1 - math.exp(-l)] if l != 1 else \
        [2 * math.exp(-1), 1 - math.exp(-1), math.exp(-1)]
    assert all(abs(v - x) > 1e-3 for x in leurres)

# proba-l2-discretes-lois-usuelles-ex2 (+ v1, v2) : hypergéométrique, leurres (binomiale, oubli, autre événement).
for N, K, n, j, rep, esp, leurres in [
    (10, 4, 3, 2, F(3, 10), F(6, 5), [F(36, 125), F(1, 20), F(1, 3)]),
    (12, 5, 4, 2, F(14, 33), F(5, 3), [F(1225, 3456), F(2, 99), F(19, 33)]),
    (8, 3, 3, 1, F(15, 28), F(9, 8), [F(225, 512), F(3, 56), F(5, 7)]),
]:
    h = lambda i: F(C(K, i) * C(N - K, n - i), C(N, n))
    assert h(j) == rep and F(n * K, N) == esp
    assert leurres[0] == C(n, j) * F(K, N) ** j * F(N - K, N) ** (n - j)
    assert leurres[1] == F(C(K, j), C(N, n))
    autre = sum(h(i) for i in range(j, n + 1)) if j == 2 else sum(h(i) for i in range(0, j + 1))
    assert leurres[2] == autre and rep not in leurres

# proba-l2-discretes-lois-usuelles-ex3 (+ v1, v2) : P(A gagne) = 1/(2-p) ; leurres 1/2, (1-p)/(2-p), 1/(1+p).
for p, rep, leurres in [(F(1, 6), F(6, 11), [F(1, 2), F(5, 11), F(6, 7)]), (F(1, 3), F(3, 5), [F(1, 2), F(2, 5), F(3, 4)]),
                        (F(1, 4), F(4, 7), [F(1, 2), F(3, 7), F(4, 5)])]:
    q = 1 - p
    assert p / (1 - q * q) == rep == 1 / (2 - p)
    assert leurres == [F(1, 2), 1 - rep, 1 / (1 + p)] and rep not in leurres
    assert abs(sum(float(q) ** (2 * i) * float(p) for i in range(500)) - float(rep)) < 1e-12

# proba-l2-discretes-lois-usuelles-ex4 (+ v1, v2) : amincissement de Poisson, X ~ P(λp), P(X=1) = λp e^-λp.
for l, p in [(5, F(2, 5)), (6, F(1, 2)), (10, F(2, 5))]:
    for kk in range(5):
        tot = sum(math.exp(-l) * l ** n / math.factorial(n) * C(n, kk) * float(p) ** kk * float(1 - p) ** (n - kk)
                  for n in range(kk, 150))
        lp = l * float(p)
        assert abs(tot - math.exp(-lp) * lp ** kk / math.factorial(kk)) < 1e-12
for lp, lam_, approx, leurres in [(2, 5, 0.271, [math.exp(-2), 2 * math.exp(-5), 2 * math.exp(-3)]),
                                  (3, 6, 0.149, [math.exp(-3), 3 * math.exp(-6), 6 * math.exp(-3)]),
                                  (4, 10, 0.073, [math.exp(-4), 4 * math.exp(-10), 4 * math.exp(-6)])]:
    v = lp * math.exp(-lp)
    assert abs(v - approx) < 5e-4 and all(abs(v - x) > 1e-3 for x in leurres)

# ---------------------------------------------------------------------------
# proba-l2-discretes-esperance (carte)
# ---------------------------------------------------------------------------
# Exemple : P(X=k) = 2^-k : E = 2, E(X(X-1)) = 4, E(X²) = 6, V = 2, E(1/X) = ln 2.
assert summation(k * Rational(1, 2) ** k, (k, 1, oo)) == 2
assert summation(k * (k - 1) * Rational(1, 2) ** k, (k, 1, oo)) == 4
assert summation(k ** 2 * Rational(1, 2) ** k, (k, 1, oo)) == 6
assert simplify(summation(Rational(1, 2) ** k / k, (k, 1, oo)) - log(2)) == 0
# Contre-exemples : 1/(k(k+1)) de somme 1 ; 4/(k(k+1)(k+2)) de somme 1, E = 2, E(X²) infinie.
assert summation(1 / (k * (k + 1)), (k, 1, oo)) == 1
assert summation(4 / (k * (k + 1) * (k + 2)), (k, 1, oo)) == 1
assert summation(4 / ((k + 1) * (k + 2)), (k, 1, oo)) == 2
assert summation(4 * k / ((k + 1) * (k + 2)), (k, 1, oo)) == oo
# Série semi-convergente : valeurs (-1)^k 2^k / k distinctes, sum (-1)^k/k = -ln 2, sum 2^-k = 1.
vals = [F((-1) ** j * 2 ** j, j) for j in range(1, 60)]
assert len(set(vals)) == len(vals)
assert abs(sum((-1) ** j / j for j in range(1, 200001)) + math.log(2)) < 1e-5
# Trous : φ(x) = x², X uniforme sur {0,1} : E(X²) = 1/2, φ(E X) = 1/4, Σ φ(x P) = 1/4, Σ x φ(P) = 1/4.
assert F(1, 2) != F(1, 4)

# proba-l2-discretes-esperance-ex1 (+ v1, v2)
assert summation(2 * 2 ** k * Rational(1, 3) ** k, (k, 1, oo)) == 4
assert summation(2 * k * Rational(1, 3) ** k, (k, 1, oo)) == Rational(3, 2)
assert summation(3 * 2 ** k * Rational(1, 4) ** k, (k, 1, oo)) == 3
assert summation(3 * k * Rational(1, 4) ** k, (k, 1, oo)) == Rational(4, 3)
assert summation((-1) ** k * Rational(1, 2) ** k, (k, 1, oo)) == Rational(-1, 3)
assert summation(Rational(1, 4) ** k, (k, 1, oo)) == Rational(1, 3)  # P(X pair), contrôle de v2
assert abs(2 ** 1.5 - 2.83) < 5e-3 and abs(2 ** (4 / 3) - 2.52) < 5e-3
assert 4 not in (2 * math.sqrt(2), 3, 6) and 3 not in (2 ** (4 / 3), 4, 1.5)

# proba-l2-discretes-esperance-ex2 (+ v1, v2)
assert summation(18 / (k * (k + 1) * (k + 2) * (k + 3)), (k, 1, oo)) == 1
assert summation(18 / ((k + 1) * (k + 2) * (k + 3)), (k, 1, oo)) == Rational(3, 2)
assert summation(18 * k / ((k + 1) * (k + 2) * (k + 3)), (k, 1, oo)) == Rational(9, 2)
assert simplify(k / ((k + 1) * (k + 2) * (k + 3)) - (-Rational(1, 2) / (k + 1) + 2 / (k + 2) - Rational(3, 2) / (k + 3))) == 0
assert -Rational(1, 2) * (Rational(1, 2) + Rational(1, 3)) + 2 * Rational(1, 3) == Rational(1, 4)
assert Rational(9, 2) - Rational(9, 4) == Rational(9, 4)
assert simplify(3 / (k * (k + 1) * (k + 2) * (k + 3)) - (1 / (k * (k + 1) * (k + 2)) - 1 / ((k + 1) * (k + 2) * (k + 3)))) == 0
assert simplify(2 / (k * (k + 1) * (k + 2)) - (1 / (k * (k + 1)) - 1 / ((k + 1) * (k + 2)))) == 0
assert summation(k / (k * (k + 1)), (k, 1, oo)) == oo  # v2 : pas d'espérance

# proba-l2-discretes-esperance-ex3 (+ v1, v2) : E(X) = sum P(X >= k).
assert summation(1 / k ** 2, (k, 1, oo)) == pi ** 2 / 6
assert simplify(1 / k ** 2 - 1 / (k + 1) ** 2 - (2 * k + 1) / (k ** 2 * (k + 1) ** 2)) == 0
assert summation(1 / k, (k, 1, oo)) == oo
assert simplify((k + 1) / 2 ** k - (k + 2) / 2 ** (k + 1) - k / 2 ** (k + 1)) == 0
assert summation(k / 2 ** (k + 1), (k, 1, oo)) == 1 and summation((k + 1) / 2 ** k, (k, 1, oo)) == 3
assert summation(k * k / 2 ** (k + 1), (k, 1, oo)) == 3  # calcul direct de E(X) pour v2
# Formule de queue contrôlée sur la loi 4/(k(k+1)(k+2)) : P(X >= k) = 2/(k(k+1)), somme 2 = E(X).
assert all(1 - sum(F(4, j * (j + 1) * (j + 2)) for j in range(1, n)) == F(2, n * (n + 1)) for n in range(1, 30))

# proba-l2-discretes-esperance-ex4 (+ v1, v2) : E(min(2^T, 2^m)) = m + 1.
for m in (20, 10, 30):
    e = sum(F(2 ** j, 2 ** j) for j in range(1, m + 1)) + 2 ** m * F(1, 2 ** m)
    assert e == m + 1 and m + 1 not in (m, 2 ** m, 2 * m + 1)
assert abs(2 ** 20 / 1e6 - 1) < 0.05 and abs(2 ** 10 / 1e3 - 1) < 0.05 and abs(2 ** 30 / 1e9 - 1) < 0.08

# ---------------------------------------------------------------------------
# proba-l2-discretes-independance (carte)
# ---------------------------------------------------------------------------
# Contre-exemple : X uniforme sur {-1,0,1}, Y = X² : Cov = 0, non indépendantes.
X3 = {-1: F(1, 3), 0: F(1, 3), 1: F(1, 3)}
assert sum(x ** 3 * v for x, v in X3.items()) - sum(x * v for x, v in X3.items()) * sum(x * x * v for x, v in X3.items()) == 0
assert F(1, 3) * F(2, 3) == F(2, 9) != 0
# Z = XY : deux à deux indépendantes, pas mutuellement.
issues = list(product([-1, 1], repeat=2))
for i, j in [(0, 1), (0, 2), (1, 2)]:
    for a, b in product([-1, 1], repeat=2):
        cpt = sum(1 for x, y in issues if (x, y, x * y)[i] == a and (x, y, x * y)[j] == b)
        assert F(cpt, 4) == F(1, 4)
assert sum(1 for x, y in issues if x == 1 and y == 1 and x * y == -1) == 0
# Exemple : V(X + Y) = 2 + 2 = 4 ; dé : V = 35/12, 10 lancers : 175/6.
assert F(1, 2) / F(1, 4) == 2 and 10 * F(35, 12) == F(175, 6)
assert sum(F(j * j, 6) for j in range(1, 7)) - F(7, 2) ** 2 == F(35, 12)
# Trous (leurres) : indépendantes Bernoulli 1/2 : V(X+Y) = 1/2, V V = 1/16, (σ+σ)² = 1.
assert F(1, 4) + F(1, 4) == F(1, 2) and F(1, 4) * F(1, 4) == F(1, 16) and (F(1, 2) + F(1, 2)) ** 2 == 1

# proba-l2-discretes-independance-ex1 (+ v1, v2)
for a, b, c0, ex_, vx, ey, vy, e_c, v_c, exy, leurres in [
    (2, -1, 1, 2, 1, 3, 4, 2, 8, 6, [6, 0, 9]),
    (3, 2, -5, 1, 2, -2, 1, -6, 22, -2, [8, 14, 47]),
    (1, -3, 0, 0, 5, 4, 2, -12, 23, 0, [11, -13, 7]),
]:
    assert a * ex_ + b * ey + c0 == e_c and a * a * vx + b * b * vy == v_c and ex_ * ey == exy
    assert v_c not in leurres
assert [2 * 1 + 4, 4 * 1 - 4, 8 + 1] == [6, 0, 9]
assert [3 * 2 + 2 * 1, 18 - 4, 22 + 25] == [8, 14, 47]
assert [5 + 3 * 2, 5 - 18, 5 + 2] == [11, -13, 7]

# proba-l2-discretes-independance-ex2 (+ v1, v2) : Cov(X, 1_{X pair}) = q/(1+q)^2.
for p, ey, exy, rep in [(F(1, 2), F(1, 3), F(8, 9), F(2, 9)), (F(1, 3), F(2, 5), F(36, 25), F(6, 25)),
                        (F(1, 4), F(3, 7), F(96, 49), F(12, 49))]:
    q = 1 - p
    assert p * q / (1 - q * q) == ey == q / (1 + q)
    assert 2 * p * q / (1 - q * q) ** 2 == exy
    assert exy - ey / p == rep == q / (1 + q) ** 2
    assert abs(sum(2 * j * float(q) ** (2 * j - 1) * float(p) for j in range(1, 2000)) - float(exy)) < 1e-10
    assert rep not in (0, exy, ey / p)

# proba-l2-discretes-independance-ex3 : X uniforme sur {-2,-1,1,2} : Cov(X, X²) = 0 et non indépendance.
X4 = {x: F(1, 4) for x in (-2, -1, 1, 2)}
assert sum(x ** 3 * v for x, v in X4.items()) == 0 and sum(x * v for x, v in X4.items()) == 0
assert F(1, 4) != F(1, 4) * F(1, 2)  # P(X=1, X²=1) = 1/4 et P(X=1) P(X²=1) = 1/8
assert all(abs(x) ** 3 <= 1 + x ** 4 and x * x <= 1 + x ** 4 for x in [i / 10 for i in range(-50, 51)])

# proba-l2-discretes-independance-ex4 (+ v1, v2) : collectionneur.
for N, esp, rep, leurres in [(3, F(11, 2), F(27, 4), [F(11, 2), F(5, 2), F(49, 4)]),
                             (4, F(25, 3), F(130, 9), [F(25, 3), F(13, 3), F(205, 9)]),
                             (5, F(137, 12), F(3625, 144), [F(137, 12), F(77, 12), F(5269, 144)])]:
    ps = [F(N - i + 1, N) for i in range(1, N + 1)]
    assert sum(1 / p for p in ps) == esp == N * H(N)
    assert sum((1 - p) / p ** 2 for p in ps) == rep == N * N * sum(F(1, j * j) for j in range(1, N + 1)) - N * H(N)
    assert leurres == [esp, sum((1 - p) / p for p in ps), N * N * sum(F(1, j * j) for j in range(1, N + 1))]
    assert rep not in leurres
assert F(1, 3) / F(2, 3) ** 2 == F(3, 4) and F(2, 3) / F(1, 3) ** 2 == 6
assert [F(1, 4) / F(3, 4) ** 2, F(1, 2) / F(1, 2) ** 2, F(3, 4) / F(1, 4) ** 2] == [F(4, 9), 2, 12]
assert [F(1, 5) / F(4, 5) ** 2, F(2, 5) / F(3, 5) ** 2, F(3, 5) / F(2, 5) ** 2, F(4, 5) / F(1, 5) ** 2] == [F(5, 16), F(10, 9), F(15, 4), 20]

# ---------------------------------------------------------------------------
# proba-l2-discretes-markov (carte)
# ---------------------------------------------------------------------------
# Exemple : X ~ P(4) : Markov 1/3, BT 1/16, vraie valeur ≈ 9,15e-4.
vraie = 1 - sum(math.exp(-4) * 4 ** j / math.factorial(j) for j in range(12))
assert abs(vraie - 9.15e-4) < 1e-6 and vraie < F(1, 16) < F(1, 3)
# Contre-exemples : X ∈ {-10, 10} ; égalité X ∈ {0, 2}, a = 2 ; X ∈ {-1, 1}, a = 1.
assert F(1, 2) > 0 and F(0 + 2, 2) / 2 == F(1, 2) and F(1 + 1, 2) - 0 == 1
# Trou 2 (leurres) : X = 2, a = 2 : 1 > 2/4 ; X = 1, a = 1/2 : 1 > 1/2 ; X = 1, a = 1 : 1 > V/a = 0.

# proba-l2-discretes-markov-ex1 (+ v1, v2)
for m, a, rep, leurres in [(20, 50, F(2, 5), [F(3, 5), F(4, 25), F(1, 125)]), (12, 48, F(1, 4), [F(3, 4), F(1, 16), F(1, 192)]),
                           (30, 40, F(3, 4), [F(1, 4), F(9, 16), F(3, 160)])]:
    assert F(m, a) == rep and leurres == [1 - rep, rep ** 2, F(m, a * a)]
    assert a * rep == m  # la variable a·1 avec probabilité m/a réalise l'égalité

# proba-l2-discretes-markov-ex2 (+ v1, v2)
for m, v, lo, hi, a, rep, leurres in [(50, 25, 40, 60, 10, F(3, 4), [F(1, 4), F(1, 2), F(15, 16)]),
                                      (100, 16, 90, 110, 10, F(21, 25), [F(4, 25), F(3, 5), F(24, 25)]),
                                      (10, 4, 7, 13, 3, F(5, 9), [F(4, 9), F(1, 3), F(8, 9)])]:
    assert m - lo == hi - m == a and 1 - F(v, a * a) == rep
    assert leurres == [F(v, a * a), 1 - F(math.isqrt(v), a), 1 - F(v, (2 * a) ** 2)]

# proba-l2-discretes-markov-ex3 (+ v1, v2) : géométrique, Markov, BT, valeur exacte.
for p, a, markov_, bt_, faux, approx in [(F(1, 4), 16, F(1, 4), F(1, 12), F(3, 64), 0.013), (F(1, 2), 10, F(1, 5), F(1, 32), F(1, 50), 0.002),
                                         (F(1, 3), 9, F(1, 3), F(1, 6), F(2, 27), 0.039)]:
    q = 1 - p
    E1, V1 = 1 / p, q / p ** 2
    assert F(E1, a) == markov_ and V1 / (a - E1) ** 2 == bt_ and V1 / a ** 2 == faux
    exact = q ** (a - 1)
    assert abs(float(exact) - approx) < 5e-4 and exact < bt_ < markov_ and bt_ != faux
assert F(1, 2) ** 9 == F(1, 512)

# proba-l2-discretes-markov-ex4 (+ v1, v2) : borne de Chernoff e^{a-λ} (λ/a)^a.
for l, a, d, approx, exact_approx, markov_, bt_ in [(1, 5, 4, 0.0175, 0.0037, F(1, 5), F(1, 16)),
                                                    (2, 6, 4, 0.075, 0.017, F(1, 3), F(1, 8)),
                                                    (2, 8, 6, 0.0062, 0.0011, F(1, 4), F(1, 18))]:
    borne = math.exp(a - l) * (l / a) ** a
    t0 = math.log(a / l)
    h = lambda t: l * (math.exp(t) - 1) - t * a
    assert abs(math.exp(h(t0)) - borne) < 1e-12 and all(h(t0) <= h(t0 + s) for s in (-0.3, -0.01, 0.01, 0.3))
    assert abs(borne - approx) < 5e-4 * max(1, approx * 100)
    exact = 1 - sum(math.exp(-l) * l ** j / math.factorial(j) for j in range(a))
    assert abs(exact - exact_approx) < 5e-4 and exact < borne < float(bt_) < float(markov_)
    assert F(l, a) == markov_ and F(l, d * d) == bt_ and a - l == d
    # E(e^{tX}) = e^{λ(e^t - 1)} (transfert), contrôlé en t = t0.
    assert abs(sum(math.exp(t0 * j) * math.exp(-l) * l ** j / math.factorial(j) for j in range(80)) - math.exp(l * (math.exp(t0) - 1))) < 1e-9
assert abs(math.exp(4) / 3125 - 0.01747) < 1e-5 and abs(math.exp(4) / 729 - 0.07489) < 1e-5 and abs(math.exp(6) / 65536 - 0.006156) < 1e-6
assert 4 ** 8 == 65536 and 3 ** 6 == 729 and 5 ** 5 == 3125
assert math.exp(5) / 3125 > math.exp(4) / 3125  # leurre e^5/3125 différent

# ---------------------------------------------------------------------------
# QCM proba-l2-discretes
# ---------------------------------------------------------------------------
# q2 : F(n) = n/(n+1) : P(X=3) = 1/12 ; leurres 3/4, 1/4, 1/6.
Fq = lambda n: F(n, n + 1)
assert Fq(3) - Fq(2) == F(1, 12) and [Fq(3), 1 - Fq(3), Fq(2) - Fq(1)] == [F(3, 4), F(1, 4), F(1, 6)]
# q3 : G(1/3) : P(X > 4) = 16/81.
assert F(2, 3) ** 4 == F(16, 81) and F(2, 3) ** 3 == F(8, 27) and 1 - F(16, 81) == F(65, 81)
# q4 : P(X pair) = (1 + e^{-2λ})/2 > 1/2.
for l in (0.1, 1, 3):
    assert abs(sum(math.exp(-l) * l ** (2 * j) / math.factorial(2 * j) for j in range(60)) - (1 + math.exp(-2 * l)) / 2) < 1e-12
# q6 : V(X - 2Y) = 2 + 4·2 = 10 ; leurres 6, -6, 18.
assert 2 + 4 * 2 == 10 and [2 + 2 * 2, 2 - 4 * 2, 2 + 4 * 4] == [6, -6, 18]
# q8 : Markov 3/12 = 1/4, atteint par X ∈ {0, 12}.
assert F(3, 12) == F(1, 4) and 12 * F(1, 4) == 3

print("proba_l2_discretes : toutes les vérifications passent.")
