"""Vérifications calculatoires du chapitre « proba-l2-generatrices » (probabilités, L2) :
cartes, exercices et QCM de content/*/probabilites/l2-proba-l2-generatrices.json et content/qcm/proba-l2-generatrices.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/proba_l2_generatrices.py (doit finir sans erreur).
"""

import itertools
import math
from fractions import Fraction as F

from sympy import (E, Poly, Rational, binomial, diff, exp, expand, factor, factorial, log, nsimplify, oo, pi, series,
                   simplify, solve, sqrt, symbols, zeta)

t = symbols("t")
lam, mu, p = symbols("lambda mu p", positive=True)
q = 1 - p


def coeffs(expr, n):
    """Coefficients 0..n du développement en série entière de expr en 0."""
    s = series(expr, t, 0, n + 1).removeO()
    return [simplify(s.coeff(t, k)) for k in range(n + 1)]


def moments(G):
    """(E, V) à partir de la fonction génératrice G (G'(1), G''(1) + G'(1) - G'(1)^2)."""
    g1 = simplify(diff(G, t).subs(t, 1))
    g2 = simplify(diff(G, t, 2).subs(t, 1))
    return g1, simplify(g2 + g1 - g1**2), g2


# ---------------------------------------------------------------------------
# proba-l2-generatrices-caracterisation (carte)
# ---------------------------------------------------------------------------
# Exemple : fonctions génératrices usuelles.
assert [simplify(c - exp(-lam) * lam**k / factorial(k)) for k, c in enumerate(coeffs(exp(lam * (t - 1)), 6))] == [0] * 7
pp = Rational(1, 3)
assert coeffs(pp * t / (1 - (1 - pp) * t), 6) == [0] + [pp * (1 - pp) ** (k - 1) for k in range(1, 7)]
assert expand((Rational(1, 2) + t / 2) ** 3 - (1 + t) ** 3 / 8) == 0  # exemple B(3, 1/2)
# Contre-exemple 1 (et QCM q2, piège « même espérance et variance ») : P(1) et {0,2} équiprobables.
GX, GY = exp(t - 1), (1 + t**2) / 2
assert moments(GX)[:2] == moments(GY)[:2] == (1, 1)
assert GX.subs(t, 0) != GY.subs(t, 0)
# Contre-exemple 2 et affirmation 2 (piège « G_X(t0) = G_Y(t0) ») : B(1/2) et P(Y=0)=2/3, P(Y=2)=1/3 en t = 1/2.
assert F(1, 2) + F(1, 2) * F(1, 2) == F(3, 4) == F(2, 3) + F(1, 3) * F(1, 4)
# Contre-exemple 3 : sum 6/(pi^2 n^2) = 1 ; rayon 1 (racine n-ième de n^-2 -> 1, la série diverge en t > 1).
assert abs(6 / math.pi**2 * sum(1 / k**2 for k in range(1, 200000)) - 1) < 1e-5
# Trou 1 : « égal à 1 » faux (Bernoulli : polynôme, rayon infini) ; « > 1 » et « infini » faux (6/(pi^2 n^2) : rayon 1).
assert all(abs((1 / k**2) ** (1 / k) - 1) < 0.01 for k in [10**4, 10**5])
# Trou 2 : X = 2 constante, G = t^2 : G''(0) = 2 != 1, n! G''(0) = 4 != 1 ; X = 1 constante, G(1)/0! = 1 != P(X=0) = 0.
G = t**2
assert diff(G, t, 2).subs(t, 0) == 2 and 2 * diff(G, t, 2).subs(t, 0) == 4 and diff(G, t, 2).subs(t, 0) / 2 == 1
assert (t).subs(t, 1) == 1  # X = 1 : G_X(1)/0! = 1 alors que P(X=0) = 0
# Trou 3 : X = 0 et Y = 1 (G(1) égaux), X = 1 et Y = 2 (G(0) = 0 et même rayon infini) : lois différentes.
assert (1 + 0 * t).subs(t, 1) == t.subs(t, 1) and t.subs(t, 0) == (t**2).subs(t, 0)
# Affirmation 4 : croissance sur [0, 1] (coefficients positifs) — rien à calculer.
# Affirmation 5 : P(X pair) = (1 + G(-1))/2, testé sur une loi finie.
loi = {0: F(1, 10), 1: F(2, 10), 2: F(3, 10), 5: F(4, 10)}
assert (1 + sum(v * (-1) ** k for k, v in loi.items())) / 2 == loi[0] + loi[2]

# proba-l2-generatrices-caracterisation-ex1 (+ v1, v2) : reconnaître la loi, P(X = 2).
c = coeffs((1 + 2 * t) ** 3 / 27, 3)
assert c == [binomial(3, k) * Rational(2, 3) ** k * Rational(1, 3) ** (3 - k) for k in range(4)] and c[2] == Rational(4, 9)
assert expand((1 + 2 * t) ** 3 / 27 - (Rational(2, 3) + t / 3) ** 3) != 0  # leurre B(3, 1/3)
assert expand((1 + 2 * t) ** 3 / 27 - (Rational(1, 3) + 2 * t / 3) ** 2) != 0  # leurre B(2, 2/3)
assert c != [Rational(1, 4)] * 4  # leurre uniforme
c = coeffs((3 + t) ** 4 / 256, 4)
assert c == [binomial(4, k) * Rational(1, 4) ** k * Rational(3, 4) ** (4 - k) for k in range(5)] and c[2] == Rational(27, 128)
assert c != [binomial(4, k) * Rational(3, 4) ** k * Rational(1, 4) ** (4 - k) for k in range(5)] and c != [Rational(1, 5)] * 5
assert c[4] != 0  # leurre B(3, 1/4) : degré 3
c = coeffs(exp(3 * t - 3), 4)
assert c == [exp(-3) * 3**k / factorial(k) for k in range(5)] and simplify(c[2] - Rational(9, 2) * exp(-3)) == 0
# Leurres P(1/3), G(1/3) sur N* (P(X=0)=0), B(3,1/2) (P(X=0)=1/8) : P(X=0) = e^-3 les distingue.
assert c[0] != exp(-Rational(1, 3)) and c[0] != 0 and c[0] != Rational(1, 8)

# proba-l2-generatrices-caracterisation-ex2 (+ v1, v2) : lire la loi sur le développement.
for G, k, rep, leurres in [(1 / (2 - t), 3, Rational(1, 16), [Rational(1, 8), Rational(1, 32), Rational(3, 16)]),
                           (2 / (3 - t), 2, Rational(2, 27), [Rational(2, 9), Rational(1, 9), Rational(2, 81)]),
                           (t / (2 - t), 4, Rational(1, 16), [Rational(1, 32), Rational(1, 8), Rational(1, 64)])]:
    c = coeffs(G, 8)
    assert all(x >= 0 for x in c) and G.subs(t, 1) == 1 and c[k] == rep and rep not in leurres
assert coeffs(1 / (2 - t), 5) == [Rational(1, 2 ** (n + 1)) for n in range(6)]

# proba-l2-generatrices-caracterisation-ex3 (+ v1, v2) : P(X pair) = (1 + G(-1))/2.
a = (1 + exp(-4)) / 2
assert abs(float(a) - sum(math.exp(-2) * 2 ** (2 * k) / math.factorial(2 * k) for k in range(40))) < 1e-14
for leurre in [(1 - exp(-4)) / 2, (1 + exp(-2)) / 2, Rational(1, 2)]:
    assert abs(float(a - leurre)) > 1e-3
b = F(1, 2) * (1 + F(1, 2) ** 10)
assert b == F(1025, 2048) == sum(F(math.comb(10, k)) * F(1, 4) ** k * F(3, 4) ** (10 - k) for k in range(0, 11, 2))
assert b not in [F(1023, 2048), F(1, 2), (1 + F(3, 4) ** 10) / 2]
g = F(1, 3) * (-1) / (1 - F(2, 3) * (-1))
assert g == F(-1, 5) and (1 + g) / 2 == F(2, 5)
assert abs(sum(1 / 3 * (2 / 3) ** (2 * k - 1) for k in range(1, 200)) - 0.4) < 1e-12
assert F(2, 5) not in [F(3, 5), F(1, 2), F(1, 3)]

# proba-l2-generatrices-caracterisation-ex4 (+ v1, v2) : P(X dans 3N) pour B(n, 1/2), n = 6, 9, 12.
j = complex(math.cos(2 * math.pi / 3), math.sin(2 * math.pi / 3))
assert all(abs(1 + j**n + j ** (2 * n) - (3 if n % 3 == 0 else 0)) < 1e-12 for n in range(30))
for n, rep, leurres in [(6, F(11, 32), [F(1, 3), F(21, 64), F(5, 16)]),
                        (9, F(85, 256), [F(1, 3), F(43, 256), F(257, 768)]),
                        (12, F(683, 2048), [F(1, 3), F(2047, 6144), F(231, 1024)])]:
    direct = F(sum(math.comb(n, k) for k in range(0, n + 1, 3)), 2**n)
    signe = 1 if n % 2 == 0 else -1  # (1+j)^n = e^{i n pi/3} = (-1)^(n/3)... : n = 6, 12 -> 1 ; n = 9 -> -1
    assert abs((1 + j) ** n - signe) < 1e-9 and abs((1 + j**2) ** n - signe) < 1e-9
    assert direct == rep == (1 + F(2 * signe, 2**n)) / 3 and rep not in leurres
assert F(5, 16) == F(math.comb(6, 3), 64) and F(21, 64) == F(1 + 20, 64)
assert F(43, 256) == F(1 + 84 + 1, 512) and F(257, 768) == (1 + F(2, 512)) / 3
assert F(2047, 6144) == (1 - F(2, 4096)) / 3 and F(231, 1024) == F(924, 4096)

# proba-l2-generatrices-caracterisation-ex5 : moments de X et Y (N = 3).
X = {0: F(1, 2), 2: F(1, 2)}
Y = {0: F(1, 3), 1: F(1, 2), 3: F(1, 6)}
mom = lambda L, k: sum(v * x**k for x, v in L.items())
assert mom(X, 1) == mom(Y, 1) == 1 and mom(X, 2) == mom(Y, 2) == 2 and mom(X, 3) == 4 and mom(Y, 3) == 5
assert sum(Y.values()) == 1
# Taylor exacte au point 1 et G^{(k)}(1) = moment factoriel, sur Y.
GYp = sum(v * t**x for x, v in Y.items())
assert expand(sum(diff(GYp, t, k).subs(t, 1) / factorial(k) * (t - 1) ** k for k in range(4)) - GYp) == 0
assert diff(GYp, t, 2).subs(t, 1) == sum(v * x * (x - 1) for x, v in Y.items())

# ---------------------------------------------------------------------------
# proba-l2-generatrices-moments (carte)
# ---------------------------------------------------------------------------
# Exemple : Poisson et géométrique.
assert moments(exp(lam * (t - 1)))[:2] == (lam, lam)
Eg, Vg, _ = moments(p * t / (1 - q * t))
assert simplify(Eg - 1 / p) == 0 and simplify(Vg - q / p**2) == 0
assert simplify(diff(p * t / (1 - q * t), t, 2) - 2 * p * q / (1 - q * t) ** 3) == 0
# Contre-exemple 2 : c/n^3, espérance c zeta(2) finie, sum n(n-1) c/n^3 diverge (série harmonique moins zeta(2)).
assert float(zeta(2) / zeta(3)) > 0
# Trous : B(1/2) et X = 2 réfutent les leurres.
GB = (1 + t) / 2
EB, VB, G2B = moments(GB)
assert (EB, VB, G2B) == (Rational(1, 2), Rational(1, 4), 0)
assert G2B - EB**2 != VB and G2B + EB != VB and G2B != VB  # trou 3 (et affirmation 3 : -1/4)
assert G2B - EB**2 == Rational(-1, 4)
G2 = t**2  # X = 2 : E = 2
assert G2.subs(t, 1) != 2 and diff(G2, t).subs(t, 0) != 2 and diff(G2, t).subs(t, 1) ** 2 != 2  # trou 2

# proba-l2-generatrices-moments-ex1 (+ v1, v2).
for G, (Ev, Vv, G2v), leurres in [
        (t * (1 + t) ** 2 / 4, (2, Rational(1, 2), Rational(5, 2)), [(2, Rational(5, 2)), (2, Rational(9, 2)), (1, Rational(1, 2))]),
        (exp(2 * (t**2 - 1)), (4, 8, 20), [(4, 20), (4, 24), (2, 2)]),
        (t**2 * (1 + t) / 2, (Rational(5, 2), Rational(1, 4), 4), [(Rational(5, 2), 4), (Rational(5, 2), Rational(13, 2)), (Rational(3, 2), Rational(1, 4))])]:
    assert moments(G) == (Ev, Vv, G2v) and (Ev, Vv) not in leurres
    assert leurres[1] == (Ev, G2v + Ev)
assert coeffs(exp(2 * (t**2 - 1)), 4) == [exp(-2), 0, 2 * exp(-2), 0, 2 * exp(-2)]  # X = 2Y, Y ~ P(2)

# proba-l2-generatrices-moments-ex2 : Y = nombre d'échecs avant le premier succès.
GY2 = p / (1 - q * t)
assert [simplify(x) for x in coeffs(GY2, 4)] == [simplify(p * q**k) for k in range(5)]
Ey, Vy, _ = moments(GY2)
assert simplify(Ey - q / p) == 0 and simplify(Vy - q / p**2) == 0
for (a, b) in [(1 / p, q / p**2), (q / p, q**2 / p**2), (q / p, 1 / p**2)]:
    assert simplify(a - Ey) != 0 or simplify(b - Vy) != 0

# proba-l2-generatrices-moments-ex3 : moments factoriels de la loi de Poisson, E(X^3).
x = symbols("x")
assert expand(x * (x - 1) * (x - 2) + 3 * x * (x - 1) + x - x**3) == 0
for k in range(1, 5):
    assert simplify(diff(exp(lam * (t - 1)), t, k).subs(t, 1) - lam**k) == 0
E3 = sum(n**3 * math.exp(-2.0) * 2.0**n / math.factorial(n) for n in range(80))
assert abs(E3 - (8 + 12 + 2)) < 1e-9  # lambda = 2
for f in [lambda l: l**3, lambda l: l**3 + l**2 + l, lambda l: l**3 + 3 * l**2 + 3 * l]:
    assert f(1) != 5  # vraie valeur 5 pour lambda = 1

# proba-l2-generatrices-moments-ex4 : p_n = 1/(n(n+1)), G(t) = 1 + (1-t) ln(1-t)/t, pas d'espérance.
for tv in [0.1, 0.5, 0.9]:
    s = sum(tv**n / (n * (n + 1)) for n in range(1, 4000))
    assert abs(s - (1 + (1 - tv) * math.log(1 - tv) / tv)) < 1e-12
assert simplify((1 - (1 + (1 - t) * log(1 - t) / t)) / (1 - t) + log(1 - t) / t) == 0

# proba-l2-generatrices-moments-ex5 (+ v1, v2) : deux piles consécutifs.
for pv, E_, V_, G2_, Gtxt, leurres in [
        (Rational(1, 2), 6, 22, 52, t**2 / (4 - 2 * t - t**2), [(6, 52), (6, 58), (4, 12)]),
        (Rational(1, 3), 12, 114, 246, t**2 / (9 - 6 * t - 2 * t**2), [(12, 246), (12, 258), (9, 72)]),
        (Rational(2, 3), Rational(15, 4), Rational(93, 16), Rational(129, 8), 4 * t**2 / (9 - 3 * t - 2 * t**2),
         [(Rational(15, 4), Rational(129, 8)), (Rational(15, 4), Rational(159, 8)), (Rational(9, 4), Rational(45, 16))])]:
    qv = 1 - pv
    G = pv**2 * t**2 / (1 - qv * t - pv * qv * t**2)
    assert simplify(G - Gtxt) == 0 and G.subs(t, 1) == 1
    assert moments(G) == (E_, V_, G2_) and (E_, V_) not in leurres
    assert (E_, G2_) in leurres and (E_, G2_ + E_) in leurres and leurres[2][0] == 1 / pv**2
    racines = solve(1 - qv * t - pv * qv * t**2, t)
    assert min(abs(float(r)) for r in racines) > 1
    # Récurrence P(T=n) = q P(T=n-1) + pq P(T=n-2) contre une énumération exacte des suites de lancers (n <= 12).
    PT = [F(0)] * 13
    for n in range(2, 13):
        tot = F(0)
        for s in itertools.product("PF", repeat=n):
            w = "".join(s)
            if w.endswith("PP") and "PP" not in w[:-1]:
                tot += F(pv.p, pv.q) ** w.count("P") * F(qv.p, qv.q) ** w.count("F")
        PT[n] = tot
        if n >= 3:
            assert PT[n] == F(qv.p, qv.q) * PT[n - 1] + F((pv * qv).p, (pv * qv).q) * PT[n - 2]
    assert [Rational(x.numerator, x.denominator) for x in PT] == coeffs(G, 12)
assert sorted([float(r) for r in solve(t**2 + 2 * t - 4, t)]) == sorted([-1 - math.sqrt(5), -1 + math.sqrt(5)])
assert sorted(solve(9 - 3 * t - 2 * t**2, t)) == [-3, Rational(3, 2)]
assert simplify(solve(9 - 6 * t - 2 * t**2, t)[0] ** 2 * 2 + 6 * solve(9 - 6 * t - 2 * t**2, t)[0] - 9) == 0
assert sorted([simplify(r) for r in solve(9 - 6 * t - 2 * t**2, t)], key=float) == sorted([(-3 - 3 * sqrt(3)) / 2, (-3 + 3 * sqrt(3)) / 2], key=float)

# ---------------------------------------------------------------------------
# proba-l2-generatrices-somme (carte)
# ---------------------------------------------------------------------------
# Exemple : Poisson et binomiales.
assert simplify(exp(lam * (t - 1)) * exp(mu * (t - 1)) - exp((lam + mu) * (t - 1))) == 0
# Contre-exemple 1 : X = Y ~ B(1/2).
assert expand((1 + t**2) / 2 - ((1 + t) / 2) ** 2) != 0
# Contre-exemple 2 (et piège « deux à deux indépendantes ») : X3 = X1 xor X2.
issues = [(a, b, a ^ b) for a in (0, 1) for b in (0, 1)]
for i, k in [(0, 1), (0, 2), (1, 2)]:
    for u in (0, 1):
        for v in (0, 1):
            assert F(sum(1 for w in issues if w[i] == u and w[k] == v), 4) == F(1, 4)
assert F(sum(1 for w in issues if sum(w) == 0), 4) == F(1, 4) and {sum(w) for w in issues} == {0, 2}
assert expand(((1 + t) / 2) ** 3).coeff(t, 0) == Rational(1, 8)
# Contre-exemple 3 et affirmation 2 (réciproque fausse).
loi2 = {(0, 0): F(1, 9), (1, 1): F(1, 9), (2, 2): F(1, 9), (1, 0): F(2, 9), (2, 1): F(2, 9), (0, 2): F(2, 9)}
assert sum(loi2.values()) == 1
for v in range(3):
    assert sum(w for (a, b), w in loi2.items() if a == v) == F(1, 3) == sum(w for (a, b), w in loi2.items() if b == v)
somme = [sum(w for (a, b), w in loi2.items() if a + b == s) for s in range(5)]
assert somme == [F(1, 9), F(2, 9), F(3, 9), F(2, 9), F(1, 9)]
assert expand(((1 + t + t**2) / 3) ** 2 - sum(Rational(c.numerator, c.denominator) * t**s for s, c in enumerate(somme))) == 0
assert loi2.get((0, 1), 0) == 0
# Piège « non corrélées » : (X, Y) uniforme sur {(0,1),(1,0),(1,2),(2,1)}.
pts = [(0, 1), (1, 0), (1, 2), (2, 1)]
EX = F(sum(a for a, b in pts), 4); EY = F(sum(b for a, b in pts), 4); EXY = F(sum(a * b for a, b in pts), 4)
assert EXY - EX * EY == 0
lx = sum(Rational(1, 4) * t**a for a, b in pts); ly = sum(Rational(1, 4) * t**b for a, b in pts)
assert expand(sum(Rational(1, 4) * t ** (a + b) for a, b in pts) - lx * ly) != 0
# Trou 2 : X1 = 0 et X2 = 1 (indépendantes) réfutent somme, puissance n-ième, moyenne.
assert t != 1 + t and t != 1**2 and t != (1 + t) / 2
# Affirmation 3 : X = Y = 1 : G_XY = t, G_X G_Y = t^2. Affirmation 4 : B(1,1/2) + B(1,1/4) non binomiale.
lo = [F(1, 2) * F(3, 4), F(1, 2) * F(3, 4) + F(1, 2) * F(1, 4), F(1, 2) * F(1, 4)]
assert lo == [F(3, 8), F(1, 2), F(1, 8)]
r = 1 / (2 * math.sqrt(2))
assert abs(r**2 - 1 / 8) < 1e-15 and abs((1 - r) ** 2 - 0.418) < 1e-3 and abs((1 - r) ** 2 - 3 / 8) > 0.04

# proba-l2-generatrices-somme-ex1 (+ v1, v2).
assert expand(((1 + t) / 2) ** 5).coeff(t, 2) == Rational(5, 16)
assert Rational(3, 8) == binomial(3, 2) / 8 and Rational(3, 32) == Rational(3, 8) * Rational(1, 4)
assert Rational(5, 16) not in [Rational(3, 8), Rational(1, 4), Rational(3, 32)]
c = coeffs(exp(t - 1) * exp(2 * (t - 1)), 2)
assert simplify(c[2] - Rational(9, 2) * exp(-3)) == 0
assert simplify(2 * exp(-3) - exp(-1) * 2 * exp(-2)) == 0  # leurre : P(X=0)P(Y=2)
P6 = expand((Rational(2, 3) + t / 3) ** 6)
assert P6.coeff(t, 1) == Rational(64, 243) and Rational(4, 9) == 2 * Rational(1, 3) * Rational(2, 3)
assert Rational(64, 243) not in [Rational(32, 243), Rational(4, 9), Rational(1, 3)]

# proba-l2-generatrices-somme-ex2 (+ v1, v2).
d6 = sum(t**k for k in range(1, 7)) / 6
assert expand(d6**2).coeff(t, 8) == Rational(5, 36) and expand(d6**2).coeff(t, 7) == Rational(1, 6)
assert expand(d6**2).coeff(t, 9) == Rational(1, 9) and expand(d6**2).coeff(t, 10) == Rational(1, 12)
assert expand(d6**3).coeff(t, 10) == Rational(1, 8) and expand(d6**3).coeff(t, 9) == Rational(25, 216)
assert expand(d6**3).coeff(t, 8) == Rational(7, 72)
d4 = sum(t**k for k in range(1, 5)) / 4
assert expand(d4**2).coeff(t, 5) == Rational(1, 4) and expand(d4**2).coeff(t, 4) == Rational(3, 16)
assert expand(d4**2).coeff(t, 3) == Rational(1, 8)
assert Rational(5, 16) not in [expand(d4**2).coeff(t, k) for k in range(9)]

# proba-l2-generatrices-somme-ex3 : loi de Pascal.
for r_ in range(1, 5):
    pv = Rational(2, 7); qv = 1 - pv
    c = coeffs((pv * t / (1 - qv * t)) ** r_, 12)
    assert all(c[n] == (binomial(n - 1, r_ - 1) * pv**r_ * qv ** (n - r_) if n >= r_ else 0) for n in range(13))
    Er, Vr, _ = moments((p * t / (1 - q * t)) ** r_)
    assert simplify(Er - r_ / p) == 0 and simplify(Vr - r_ * q / p**2) == 0
pv = Rational(2, 7); qv = 1 - pv
assert coeffs(pv * t / (1 - qv * t), 2)[2] != binomial(2, 1) * pv * qv  # leurre binom(n, r) (r = 1, n = 2)
assert binomial(0, 1) * pv == 0 != pv and qv != pv  # leurres binom(n-1, r) et p^{n-r} q^r (r = n = 1)
# Développement (1-u)^{-r}.
u = symbols("u")
assert all(series((1 - u) ** -3, u, 0, 8).removeO().coeff(u, k) == binomial(k + 2, 2) for k in range(8))

# proba-l2-generatrices-somme-ex4 : 1 + t + ... + t^10 sans racine réelle ; somme uniforme impossible.
P10 = Poly(sum(t**k for k in range(11)), t)
assert P10.count_roots() == 0  # nombre de racines réelles

# proba-l2-generatrices-somme-ex5 : somme aléatoire, Poisson aminci.
pv = symbols("pv", positive=True)
assert simplify(exp(lam * ((1 - pv + pv * t) - 1)) - exp(lam * pv * (t - 1))) == 0
# Contrôle numérique de G_S = G_N o G_X : N ~ P(1,5), X ~ B(0,3) ; P(S = k) par la formule de la question 1.
lv, pn = 1.5, 0.3
PS = [sum(math.exp(-lv) * lv**n / math.factorial(n) * math.comb(n, k) * pn**k * (1 - pn) ** (n - k) for n in range(k, 80))
      for k in range(6)]
assert all(abs(PS[k] - math.exp(-lv * pn) * (lv * pn) ** k / math.factorial(k)) < 1e-12 for k in range(6))

# ---------------------------------------------------------------------------
# proba-l2-generatrices-applications (carte)
# ---------------------------------------------------------------------------
# Exemple : x + 2y = n, décomposition en éléments simples et a_n = floor(n/2) + 1.
assert simplify(1 / ((1 - t) ** 2 * (1 + t)) - (1 / (2 * (1 - t) ** 2) + 1 / (4 * (1 - t)) + 1 / (4 * (1 + t)))) == 0
c = coeffs(1 / ((1 - t) * (1 - t**2)), 20)
assert all(c[n] == n // 2 + 1 == F(2 * n + 3 + (-1) ** n, 4) == sum(1 for y in range(n + 1) if 2 * y <= n) for n in range(21))
# Contre-exemples : 1/(1-t)^2 ; A = 1, B = 2t en 1/2 ; sum n! t^n de rayon nul (rapport (n+1) -> infini).
assert coeffs(1 / (1 - t) ** 2, 5) == [1, 2, 3, 4, 5, 6]
assert 2 * Rational(1, 2) == 1
# Trou 1 : leurres réfutés ; trou 2 : (1,0,...) et (0,1,0,...) ont la même somme en t = 1.
a_, b_ = [1, 1, 0], [1, 1, 0]
cn = lambda a, b, n: sum(a[k] * b[n - k] for k in range(n + 1))
assert cn(a_, b_, 1) == 2 != a_[1] * b_[1]  # a_n b_n (et affirmation 2 avec a = b = 1 : c_1 = 2)
assert cn([1, 0, 0], [0, 1, 0], 1) == 1 != sum([1, 0][k] * [0, 1][k] for k in range(2))  # sum a_k b_k
assert cn(a_, b_, 2) == 1 != a_[0] * b_[2] + a_[2] * b_[0]  # a_0 b_n + a_n b_0
# Affirmation 1 (Vandermonde) sur quelques valeurs.
assert all(sum(math.comb(m, k) * math.comb(pp_, n - k) for k in range(n + 1)) == math.comb(m + pp_, n)
           for m in range(6) for pp_ in range(6) for n in range(12))
# Affirmation 4 : (t - 1/2)(t - 1/3) = 1/6 - 5/6 t + t^2.
assert expand((t - Rational(1, 2)) * (t - Rational(1, 3)) - (Rational(1, 6) - Rational(5, 6) * t + t**2)) == 0

# proba-l2-generatrices-applications-ex1 (+ v1, v2).
c2 = coeffs(1 / ((1 - t) * (1 - t**2)), 15)
assert c2[10] == 6 and c2[15] == 8 and 6 not in [5, 11, 7] and 8 not in [7, 16, 9]
c3 = coeffs(1 / ((1 - t) * (1 - t**3)), 10)
assert c3[10] == 4 == sum(1 for y in range(4) if 3 * y <= 10) and 4 not in [3, 5, 11]

# proba-l2-generatrices-applications-ex2 (+ v1, v2).
assert sum(math.comb(5, k) ** 2 for k in range(6)) == 252 == math.comb(10, 5) and 252 not in [32, 1024, 126]
assert sum(math.comb(4, k) * math.comb(6, 3 - k) for k in range(4)) == 120 == 20 + 60 + 36 + 4
assert 120 not in [210, 84, 60] and math.comb(10, 4) == 210 and math.comb(9, 3) == 84
assert sum(math.comb(7, k) ** 2 for k in range(8)) == 3432 == math.comb(14, 7) and 3432 not in [128, 16384, 1716]

# proba-l2-generatrices-applications-ex3 (+ v1, v2) : carrelages.
phi, psi = (1 + sqrt(5)) / 2, (1 - sqrt(5)) / 2
assert simplify(phi + psi - 1) == 0 and simplify(phi * psi + 1) == 0
assert simplify((1 - phi * t) * (1 - psi * t) - (1 - t - t**2)) == 0
for (a1, a2, c0, c1, borne, formule, leurres) in [
        (1, 1, 1, 1, 2, lambda n: (phi ** (n + 1) - psi ** (n + 1)) / sqrt(5),
         [lambda n: (phi**n - psi**n) / sqrt(5), lambda n: phi ** (n + 1) - psi ** (n + 1), lambda n: (phi ** (n + 1) + psi ** (n + 1)) / sqrt(5)]),
        (1, 2, 1, 1, 2, lambda n: Rational(2 ** (n + 1) + (-1) ** n, 3),
         [lambda n: Rational(2**n + (-1) ** n, 3), lambda n: Rational(2 ** (n + 1) - (-1) ** n, 3), lambda n: 2**n]),
        (2, 3, 1, 2, 3, lambda n: Rational(3 ** (n + 1) + (-1) ** n, 4),
         [lambda n: Rational(3 ** (n + 1) - (-1) ** n, 4), lambda n: 3**n, lambda n: Rational(3**n + (-1) ** n, 2)])]:
    cs = [c0, c1]
    while len(cs) < 25:
        cs.append(a1 * cs[-1] + a2 * cs[-2])
    assert all(cs[n] <= borne**n for n in range(25))
    assert coeffs(1 / (1 - a1 * t - a2 * t**2), 24) == cs
    assert all(simplify(formule(n) - cs[n]) == 0 for n in range(25))
    assert all(any(simplify(f(n) - cs[n]) != 0 for n in range(4)) for f in leurres)
# Énumération directe des carrelages (petits n).
def carrelages(n, couleurs1, couleurs2):
    if n == 0:
        return 1
    return couleurs1 * carrelages(n - 1, couleurs1, couleurs2) + (couleurs2 * carrelages(n - 2, couleurs1, couleurs2) if n >= 2 else 0)
assert carrelages(3, 1, 1) == 3 and carrelages(2, 1, 2) == 3 and carrelages(1, 2, 3) == 2 and carrelages(2, 2, 3) == 7
assert sorted(float(r) for r in solve(1 - t - t**2, t)) == sorted([(-1 - math.sqrt(5)) / 2, (-1 + math.sqrt(5)) / 2])
assert simplify(1 / ((1 - 2 * t) * (1 + t)) - (Rational(2, 3) / (1 - 2 * t) + Rational(1, 3) / (1 + t))) == 0
assert simplify(1 / ((1 - 3 * t) * (1 + t)) - (Rational(3, 4) / (1 - 3 * t) + Rational(1, 4) / (1 + t))) == 0

# proba-l2-generatrices-applications-ex4 : dés de Sicherman.
assert expand(t * (1 + t) * (1 + t + t**2) * (1 - t + t**2) - sum(t**k for k in range(1, 7))) == 0
A = (t + 2 * t**2 + 2 * t**3 + t**4) / 6
B = (t + t**3 + t**4 + t**5 + t**6 + t**8) / 6
assert expand(A - t * (1 + t) * (1 + t + t**2) / 6) == 0
assert expand(B - t * (1 + t) * (1 + t + t**2) * (1 - t + t**2) ** 2 / 6) == 0
assert expand((1 + t) * (1 + t + t**2) * (1 - t + t**2) ** 2 - (1 + t**2 + t**3 + t**4 + t**5 + t**7)) == 0
assert expand(A * B - d6**2) == 0

# proba-l2-generatrices-applications-ex5 : nombres de Catalan.
cat = [1]
for n in range(15):
    cat.append(sum(cat[k] * cat[n - k] for k in range(n + 1)))
assert cat == [math.comb(2 * n, n) // (n + 1) for n in range(16)]
assert all(math.comb(2 * n, n) % (n + 1) == 0 for n in range(16))
s = series(sqrt(1 - 4 * t), t, 0, 12).removeO()
assert all(s.coeff(t, n) == -binomial(2 * n, n) / (2 * n - 1) for n in range(1, 12))
assert all(binomial(2 * m + 2, m + 1) / (2 * (2 * m + 1)) == binomial(2 * m, m) / (m + 1) for m in range(12))
assert all(Rational(math.comb(2 * k, k), 2 * k - 1) == Rational(2 * math.comb(2 * k - 2, k - 1), k) for k in range(1, 12))
f = (1 - sqrt(1 - 4 * t)) / (2 * t)
assert simplify(t * f**2 - f + 1) == 0
assert [simplify(x) for x in coeffs(f, 8)] == cat[:9]
for g in [lambda n: math.comb(2 * n, n), lambda n: F(math.comb(2 * n + 2, n + 1), n + 1), lambda n: F(math.comb(2 * n, n), 2 * n + 1)]:
    assert g(1) != cat[1]

# ---------------------------------------------------------------------------
# QCM proba-l2-generatrices-q1 à q8
# ---------------------------------------------------------------------------
assert expand((1 + t) ** 2 / 4 - (Rational(1, 2) + t / 2) ** 2) == 0  # q1
assert expand((1 + t) ** 2 / 4).coeff(t, 1) == Rational(1, 2) != expand((Rational(3, 4) + t / 4) ** 2).coeff(t, 1)
Eq3, Vq3, _ = moments(t / (3 - 2 * t))  # q3
assert (Eq3, Vq3) == (3, 6) and diff(t / (3 - 2 * t), t, 2).subs(t, 1) == 12
assert simplify(exp(2 * (t - 1)) * exp(3 * (t - 1))).subs(t, 0) == exp(-5)  # q5
assert simplify(diff(exp(5 * (t - 1)), t).subs(t, 0) - 5 * exp(-5)) == 0
assert sum(1 for x_ in range(6) for y_ in range(6) for z_ in range(6) if x_ + y_ + z_ == 5) == 21 == math.comb(7, 2)  # q7
assert math.comb(6, 2) == 15 and math.comb(5, 2) == 10
assert sorted((9 - 3 * y_, y_) for y_ in range(4)) == [(0, 3), (3, 2), (6, 1), (9, 0)]  # q8
assert coeffs(1 / ((1 - t) * (1 - t**3)), 9)[9] == 4
print("proba_l2_generatrices : toutes les vérifications passent.")
