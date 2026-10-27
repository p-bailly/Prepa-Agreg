"""Vérifications calculatoires du chapitre « int-l1-riemann » (thème intégration, L1).

Cartes content/cartes/integration/l1-int-l1-riemann.json et exercices
content/exercices/integration/l1-int-l1-riemann.json. Chaque assert cite l'identifiant de la carte ou de l'exercice
(« v1 », « v2 » = variantes 1 et 2). Pour chaque réponse proposée, la bonne est vérifiée égale à la valeur calculée et
chaque leurre est vérifié faux.
Lancer : python3 verification/int_l1_riemann.py (doit finir sans erreur).
"""

from sympy import (E, Rational, atan, cos, diff, exp, floor, integrate, limit, log, oo, pi, Piecewise, simplify,
                   sin, sqrt, symbols, summation, factor, apart, nsimplify, Abs, Integer, S, N)

x, t, u, s = symbols("x t u s", real=True)
n, k = symbols("n k", positive=True, integer=True)


def egal(a, b):
    return simplify(a - b) == 0


def different(a, b):
    """a et b diffèrent : contrôle exact, puis numérique (en n = 3 si l'expression dépend de n)."""
    d = simplify(S(a) - S(b))
    if d == 0:
        return False
    d = d.subs(n, 3)
    return abs(float(N(d, 30))) > 1e-9


def bonne_et_leurres(valeur, bonne, leurres, ou=""):
    assert egal(valeur, bonne), (ou, valeur, bonne)
    for l in leurres:
        assert different(valeur, l), (ou, valeur, l)


def num_integrale(f, a, b, m=4000):
    """Intégrale par la méthode de Simpson (contrôle numérique indépendant de SymPy)."""
    g = lambda z: float(f(z))
    h = (b - a) / m
    tot = g(a) + g(b)
    for i in range(1, m):
        tot += (4 if i % 2 else 2) * g(a + i * h)
    return tot * h / 3


# ===========================================================================
# Carte int-l1-riemann-cpm
# ===========================================================================
# Exemple : int_0^3 floor = 3 ; f = x^2 sur [0,1[, 5 sur [1,2] : 1/3 + 5 = 16/3.
assert sum(range(3)) == 3
assert egal(integrate(x**2, (x, 0, 1)) + 5, Rational(16, 3))
# Contre-exemple sin(1/x) : f(1/(pi/2+2k pi)) = 1 et f(1/(3pi/2+2k pi)) = -1.
for kk in range(1, 5):
    assert egal(sin(pi / 2 + 2 * kk * pi), 1) and egal(sin(3 * pi / 2 + 2 * kk * pi), -1)

# int-l1-riemann-cpm-ex1 : int_0^{7/2} floor = 0+1+2+3/2 = 9/2 ; leurres 6, 3, 49/8 (= int x).
def int_floor(a):
    a = Rational(a)
    return sum(Integer(j) * (min(a, j + 1) - j) for j in range(int(floor(a)) + 1) if j < a)


bonne_et_leurres(int_floor(Rational(7, 2)), Rational(9, 2), [6, 3, Rational(49, 8)], "cpm-ex1")
assert egal(integrate(x, (x, 0, Rational(7, 2))), Rational(49, 8))
# v1 : a = 5/2 -> 2 ; leurres 3, 1, 25/8.  v2 : a = 9/2 -> 8 ; leurres 10, 6, 81/8.
bonne_et_leurres(int_floor(Rational(5, 2)), 2, [3, 1, Rational(25, 8)], "cpm-ex1 v1")
bonne_et_leurres(int_floor(Rational(9, 2)), 8, [10, 6, Rational(81, 8)], "cpm-ex1 v2")
# contrôle numérique de la fonction int_floor
assert abs(num_integrale(lambda z: int(z), 0, 3.5, 7000) - 4.5) < 1e-2

# int-l1-riemann-cpm-ex2 : x^2 sur [-1,1[, 3-x sur [1,2] -> 13/6 ; f(1-) = 1 != f(1) = 2.
v = integrate(x**2, (x, -1, 1)) + integrate(3 - x, (x, 1, 2))
bonne_et_leurres(v, Rational(13, 6), [Rational(11, 6), Rational(7, 6), Rational(13, 3)], "cpm-ex2")
assert integrate(x**2, (x, -1, 1)) == Rational(2, 3) and integrate(3 - x, (x, 1, 2)) == Rational(3, 2)
# v1 : e^x sur [0,1[, 2x sur [1,3] -> e+7 ; f(1-) = e != 2.
v = integrate(exp(x), (x, 0, 1)) + integrate(2 * x, (x, 1, 3))
bonne_et_leurres(v, E + 7, [E + 8, E + 3, E + 9], "cpm-ex2 v1")
assert different(E, 2)
# v2 : 1/x sur [1,2[, x-1 sur [2,4] -> ln2 + 4 ; f(2-) = 1/2 != f(2) = 1.
v = integrate(1 / x, (x, 1, 2)) + integrate(x - 1, (x, 2, 4))
bonne_et_leurres(v, log(2) + 4, [log(2) + 2, log(2) + 6, 4 - log(2)], "cpm-ex2 v2")
assert integrate(x, (x, 2, 4)) == 6

# int-l1-riemann-cpm-ex3 : f = x^2, int psi_n - int phi_n = 1/n ; leurres 1/n^2, 2/n, 1/(2n).
kk = symbols("kk", integer=True, nonnegative=True)
psi = summation((kk / n) ** 2, (kk, 1, n)) / n
phi = summation((kk / n) ** 2, (kk, 0, n - 1)) / n
assert egal(psi, (n + 1) * (2 * n + 1) / (6 * n**2)) and egal(phi, (n - 1) * (2 * n - 1) / (6 * n**2))
bonne_et_leurres(psi - phi, 1 / n, [1 / n**2, 2 / n, 1 / (2 * n)], "cpm-ex3")  # contrôle symbolique en n
assert limit(psi, n, oo) == Rational(1, 3) and limit(phi, n, oo) == Rational(1, 3)
# v1 : f = x^3 sur [0,2], h = 2/n -> 16/n ; leurres 8/n, 16/n^2, 2/n.
psi = (2 / n) * summation((2 * kk / n) ** 3, (kk, 1, n))
phi = (2 / n) * summation((2 * kk / n) ** 3, (kk, 0, n - 1))
assert egal(psi, 4 * (n + 1) ** 2 / n**2) and egal(phi, 4 * (n - 1) ** 2 / n**2)
bonne_et_leurres(psi - phi, 16 / n, [8 / n, 16 / n**2, 2 / n], "cpm-ex3 v1")
assert limit(psi, n, oo) == 4 and integrate(x**3, (x, 0, 2)) == 4
# v2 : f = e^{-x} décroissante -> (1 - e^{-1})/n ; leurres (e-1)/n, 1/n, (1-e^{-1})/n^2.
psi = summation(exp(-kk / n), (kk, 0, n - 1)) / n
phi = summation(exp(-kk / n), (kk, 1, n)) / n
bonne_et_leurres(simplify(psi - phi), (1 - exp(-1)) / n, [(E - 1) / n, 1 / n, (1 - exp(-1)) / n**2], "cpm-ex3 v2")
assert egal(limit(n * (1 - exp(-1 / n)), n, oo), 1)
assert egal(integrate(exp(-x), (x, 0, 1)), 1 - exp(-1))

# int-l1-riemann-cpm-ex4 : int_0^eps 2 + eps(1-eps) <= 3 eps pour eps dans ]0,1[.
eps = symbols("eps", positive=True)
assert simplify(3 * eps - (2 * eps + eps * (1 - eps)) - eps**2) == 0  # différence = eps^2 >= 0
assert egal(sin(1 / (1 / (pi / 2 + 2 * pi))), 1) and egal(sin(1 / (1 / (3 * pi / 2 + 2 * pi))), -1)

# ===========================================================================
# Carte int-l1-riemann-sommes
# ===========================================================================
# Exemple : n/(n^2+k^2) -> pi/4 ; (1/n) sum_{k=0}^{n-1} e^{k/n} = (1/n)(e-1)/(e^{1/n}-1) -> e - 1.
assert egal(integrate(1 / (1 + x**2), (x, 0, 1)), pi / 4)
geo = summation(exp(kk / n), (kk, 0, n - 1))
assert egal(simplify(geo - (E - 1) / (exp(1 / n) - 1)), 0)
assert egal(limit((E - 1) / (n * (exp(1 / n) - 1)), n, oo), E - 1)
# Contre-exemple 1/x : (1/n) sum n/k = H_n -> +oo ; points irrationnels (k - sqrt2/2)/n dans ](k-1)/n, k/n[.
assert limit(summation(1 / kk, (kk, 1, n)), n, oo) == oo
assert 0 < float(sqrt(2) / 2) < 1
# Démonstration : R_n - R'_n = h (f(b) - f(a)) (vérifié sur un exemple, f = e^x, [0,1], n = 5).
Rn = sum(exp(Rational(j, 5)) for j in range(1, 6)) / 5
Rpn = sum(exp(Rational(j, 5)) for j in range(0, 5)) / 5
assert egal(Rn - Rpn, (E - 1) / 5)
# Affirmation 4 / ex1 : sum 1/(n+k) -> ln 2.

# int-l1-riemann-sommes-ex1 : -> ln 2 ; leurres 0, 1/2, ln 3.
bonne_et_leurres(integrate(1 / (1 + x), (x, 0, 1)), log(2), [0, Rational(1, 2), log(3)], "sommes-ex1")
assert abs(float(sum(1.0 / (10**5 + j) for j in range(1, 10**5 + 1))) - float(log(2))) < 1e-4
# v1 : sum 1/(2n+k) -> ln(3/2) ; leurres ln3, ln2, 1/3.
bonne_et_leurres(integrate(1 / (2 + x), (x, 0, 1)), log(Rational(3, 2)), [log(3), log(2), Rational(1, 3)], "sommes-ex1 v1")
# v2 : sum 1/(n+2k) -> (1/2) ln 3 ; leurres ln3, ln2, (1/2)ln2.
bonne_et_leurres(integrate(1 / (1 + 2 * x), (x, 0, 1)), log(3) / 2, [log(3), log(2), log(2) / 2], "sommes-ex1 v2")
M = 10**5
assert abs(sum(1.0 / (M + 2 * j) for j in range(1, M + 1)) - float(log(3) / 2)) < 1e-4

# int-l1-riemann-sommes-ex2 : sum_{k=1}^{2n} 1/(n+k) -> ln 3 ; leurres ln2, 2ln2, ln(3/2).
bonne_et_leurres(integrate(1 / (1 + x), (x, 0, 2)), log(3), [log(2), 2 * log(2), log(Rational(3, 2))], "sommes-ex2")
assert abs(sum(1.0 / (M + j) for j in range(1, 2 * M + 1)) - float(log(3))) < 1e-4
# v1 : sum_{k=1}^{3n} n/(n^2+k^2) -> arctan 3 ; leurres pi/4, 3pi/4, arctan(1/3).
bonne_et_leurres(integrate(1 / (1 + x**2), (x, 0, 3)), atan(3), [pi / 4, 3 * pi / 4, atan(Rational(1, 3))], "sommes-ex2 v1")
assert abs(sum(M / (M**2 + j**2) for j in range(1, 3 * M + 1)) - float(atan(3))) < 1e-4
# v2 : sum_{k=1}^{2n} k/n^2 = (2n+1)/n -> 2 ; leurres 1/2, 1, 4.
assert egal(summation(kk, (kk, 1, 2 * n)) / n**2, (2 * n + 1) / n)
bonne_et_leurres(limit((2 * n + 1) / n, n, oo), 2, [Rational(1, 2), 1, 4], "sommes-ex2 v2")

# int-l1-riemann-sommes-ex3 : prod (1+k/n)^{1/n} -> 4/e ; leurres 2ln2 - 1, 2, e/4.
I = integrate(log(1 + x), (x, 0, 1))
assert egal(I, 2 * log(2) - 1)
bonne_et_leurres(exp(I), 4 / E, [2 * log(2) - 1, 2, E / 4], "sommes-ex3")
import math
assert abs(math.exp(sum(math.log(1 + j / M) for j in range(1, M + 1)) / M) - 4 / math.e) < 1e-4
# v1 : prod (1+k^2/n^2)^{1/n} -> 2 e^{pi/2 - 2} ; leurres ln2 - 2 + pi/2, 2e^{pi/4-2}, 2e^{pi/2-1}.
I = integrate(log(1 + x**2), (x, 0, 1))
assert egal(I, log(2) - 2 + pi / 2)
bonne_et_leurres(exp(I), 2 * exp(pi / 2 - 2), [log(2) - 2 + pi / 2, 2 * exp(pi / 4 - 2), 2 * exp(pi / 2 - 1)], "sommes-ex3 v1")
# v2 : prod (2+k/n)^{1/n} -> 27/(4e) ; leurres 3ln3 - 2ln2 - 1, 27/4, 9/(4e).
I = integrate(log(2 + x), (x, 0, 1))
assert egal(I, 3 * log(3) - 1 - 2 * log(2))
bonne_et_leurres(exp(I), 27 / (4 * E), [3 * log(3) - 2 * log(2) - 1, Rational(27, 4), 9 / (4 * E)], "sommes-ex3 v2")

# int-l1-riemann-sommes-ex4 : n(R_n - int f) -> (f(1)-f(0))/2.
# Q1 : int_{x_{k-1}}^{x_k} (f(x_k) - f(t)) dt = int (t - x_{k-1}) f'(t) dt (contrôle sur f = e^t, [a, b]).
a_, b_ = symbols("a_ b_", real=True)
assert egal(integrate(exp(b_) - exp(t), (t, a_, b_)), integrate((t - a_) * exp(t), (t, a_, b_)))
assert egal(integrate(t - a_, (t, a_, a_ + 1 / n)), 1 / (2 * n**2))
# f = 1/(1+x) : (1/2 - 1)/2 = -1/4 ; leurres 1/4, -1/2, 0. Contrôle numérique pour n = 10^5.
f = lambda z: 1 / (1 + z)
bonne_et_leurres((f(Rational(1)) - f(Rational(0))) / 2, Rational(-1, 4), [Rational(1, 4), Rational(-1, 2), 0], "sommes-ex4")
nn = 10**5
assert abs(nn * (sum(1.0 / (nn + j) for j in range(1, nn + 1)) - math.log(2)) + 0.25) < 1e-3
# v1 : f = x^2 : R_n = 1/3 + 1/(2n) + 1/(6n^2) ; limite 1/2 ; leurres -1/2, 1/6, 1/3.
Rn = summation(kk**2, (kk, 1, n)) / n**3
assert egal(Rn, Rational(1, 3) + 1 / (2 * n) + 1 / (6 * n**2))
bonne_et_leurres(limit(n * (Rn - Rational(1, 3)), n, oo), Rational(1, 2), [Rational(-1, 2), Rational(1, 6), Rational(1, 3)], "sommes-ex4 v1")
# v2 : f = e^x : limite (e-1)/2 ; leurres -(e-1)/2, e/2, (e+1)/2.
Rn = summation(exp(kk / n), (kk, 1, n)) / n
bonne_et_leurres(limit(n * (Rn - (E - 1)), n, oo), (E - 1) / 2, [-(E - 1) / 2, E / 2, (E + 1) / 2], "sommes-ex4 v2")

# ===========================================================================
# Carte int-l1-riemann-changement-variable
# ===========================================================================
# Exemples : int_0^1 sqrt(1-x^2) = pi/4 ; int_{-1}^2 2t e^{t^2} = e^4 - e = int_1^4 e^x.
assert egal(integrate(cos(t) ** 2, (t, 0, pi / 2)), pi / 4)
assert egal(integrate(2 * t * exp(t**2), (t, -1, 2)), exp(4) - E) and egal(integrate(exp(x), (x, 1, 4)), exp(4) - E)
# Contre-exemple 1 : int_0^{sqrt2/2} dx/(1-x) = ln 2/(2 - sqrt 2) = ln(2 + sqrt 2) ; sin([0, 3pi/4]) contient 1.
assert egal(integrate(1 / (1 - x), (x, 0, sqrt(2) / 2)), log(2 + sqrt(2)))
assert egal(sin(3 * pi / 4), sqrt(2) / 2)
# cos t / (1 - sin t) ~ -2/(t - pi/2) : non borné.
assert egal(limit((t - pi / 2) * cos(t) / (1 - sin(t)), t, pi / 2), -2)
# Contre-exemple 2 : x = 1/u donne formellement -1/(1+u^2) ; l'intégrale vaut pi/2.
assert egal((1 / (1 + (1 / u) ** 2)) * diff(1 / u, u), -1 / (1 + u**2))
assert egal(integrate(1 / (1 + x**2), (x, -1, 1)), pi / 2)
# Contre-exemple 3 : phi = t^2 sin(1/t^2), phi'(1/sqrt(2k pi)) = -2 sqrt(2 k pi).
ph = t**2 * sin(1 / t**2)
dph = diff(ph, t)
assert egal(dph, 2 * t * sin(1 / t**2) - 2 / t * cos(1 / t**2))
for kk_ in (1, 2, 5):
    assert egal(dph.subs(t, 1 / sqrt(2 * kk_ * pi)), -2 * sqrt(2 * kk_ * pi))
# Affirmation 4 : int_0^pi f(sin t) cos t = 0 (contrôle sur f = exp).
assert egal(integrate(exp(sin(t)) * cos(t), (t, 0, pi)), 0)

# int-l1-riemann-changement-variable-ex1 : (ln x)^2/x -> 1/3 ; leurres 1/2, 1, e^3/3.
xp = symbols("xp", positive=True)
bonne_et_leurres(integrate(log(xp) ** 2 / xp, (xp, 1, E)), Rational(1, 3), [Rational(1, 2), 1, exp(3) / 3], "cv-ex1")
bonne_et_leurres(integrate(log(xp) ** 3 / xp, (xp, 1, E)), Rational(1, 4), [Rational(1, 3), 1, exp(4) / 4], "cv-ex1 v1")
bonne_et_leurres(integrate(log(xp) / xp, (xp, 1, exp(2))), 2, [1, 4, Rational(1, 2)], "cv-ex1 v2")

# int-l1-riemann-changement-variable-ex2 : sqrt(1-x^2) -> pi/4 ; leurres pi/2, 1, pi/8.
bonne_et_leurres(integrate(sqrt(1 - x**2), (x, 0, 1)), pi / 4, [pi / 2, 1, pi / 8], "cv-ex2")
v = integrate(cos(t) ** 2, (t, 0, pi / 6))
assert egal(v, integrate(sqrt(1 - x**2), (x, 0, Rational(1, 2))))
bonne_et_leurres(v, pi / 12 + sqrt(3) / 8, [pi / 12, pi / 12 + sqrt(3) / 4, pi / 6 + sqrt(3) / 8], "cv-ex2 v1")
bonne_et_leurres(integrate(sqrt(4 - x**2), (x, 0, 2)), pi, [pi / 2, 2 * pi, 4], "cv-ex2 v2")

# int-l1-riemann-changement-variable-ex3 : 1/(e^x+1) -> 1 - ln((1+e)/2).
up = symbols("up", positive=True)
assert egal(apart(1 / (up * (up + 1)), up), 1 / up - 1 / (up + 1))
v = integrate(1 / (up * (up + 1)), (up, 1, E))
assert abs(float(v) - num_integrale(lambda z: 1 / (math.exp(z) + 1), 0, 1)) < 1e-9
bonne_et_leurres(v, 1 - log((1 + E) / 2), [log((1 + E) / 2), 1 + log((1 + E) / 2), log(E / (1 + E))], "cv-ex3")
# v1 : 1/(e^x+2) -> (1/2)(1 - ln((e+2)/3)).
v = integrate(1 / (up * (up + 2)), (up, 1, E))
assert abs(float(v) - num_integrale(lambda z: 1 / (math.exp(z) + 2), 0, 1)) < 1e-9
bonne_et_leurres(v, (1 - log((E + 2) / 3)) / 2,
                 [1 - log((E + 2) / 3), log((E + 2) / 3) / 2, (1 + log((E + 2) / 3)) / 2], "cv-ex3 v1")
# v2 : e^{2x}/(e^x+1) -> e - 1 - ln((e+1)/2).
v = integrate(up / (up + 1), (up, 1, E))
assert abs(float(v) - num_integrale(lambda z: math.exp(2 * z) / (math.exp(z) + 1), 0, 1)) < 1e-9
bonne_et_leurres(v, E - 1 - log((E + 1) / 2), [E - 1 + log((E + 1) / 2), E - log((E + 1) / 2), log((E + 1) / 2)], "cv-ex3 v2")

# int-l1-riemann-changement-variable-ex4 : int_0^pi x sin x/(1+cos^2 x) = pi^2/4.
for den, val, leurres, ou in [
    (1 + cos(x) ** 2, pi**2 / 4, [pi**2 / 2, pi**2 / 8, pi / 2], "cv-ex4"),
    (3 + cos(x) ** 2, pi**2 / (6 * sqrt(3)), [pi**2 / (3 * sqrt(3)), pi**2 / 12, pi / (6 * sqrt(3))], "cv-ex4 v1"),
    (1 + sin(x) ** 2, pi / sqrt(2) * log(1 + sqrt(2)),
     [pi / (2 * sqrt(2)) * log(1 + sqrt(2)), pi * log(1 + sqrt(2)), pi / sqrt(2) * log(2)], "cv-ex4 v2"),
]:
    # intégrale réduite (pi/2) int_0^pi sin x/den, puis contrôle numérique de l'intégrale de départ.
    numv = num_integrale(lambda z: z * math.sin(z) / float(den.subs(x, z)), 0, math.pi)
    assert abs(numv - float(val)) < 1e-8, (ou, numv, float(val))
    for l in leurres:
        assert abs(numv - float(l)) > 1e-3, (ou, l)
assert egal(integrate(1 / (1 + u**2), (u, -1, 1)), pi / 2)
assert egal(integrate(1 / (3 + u**2), (u, -1, 1)), pi / (3 * sqrt(3)))
assert egal(simplify(1 / (2 - u**2) - (1 / (2 * sqrt(2))) * (1 / (sqrt(2) - u) + 1 / (sqrt(2) + u))), 0)
assert egal((sqrt(2) + 1) / (sqrt(2) - 1), (1 + sqrt(2)) ** 2)
# g(s) : sin/(1+cos^2) = s/(2-s^2), sin/(3+cos^2) = s/(4-s^2).
assert egal((sin(x) / (1 + cos(x) ** 2)) - sin(x) / (2 - sin(x) ** 2), 0)
assert egal((sin(x) / (3 + cos(x) ** 2)) - sin(x) / (4 - sin(x) ** 2), 0)

# int-l1-riemann-changement-variable-ex5 : I(a) = int_0^a ln(1+ax)/(1+x^2) = (1/2) arctan(a) ln(1+a^2).
ap = symbols("ap", positive=True)
phi_ = (ap - u) / (1 + ap * u)
assert egal(diff(phi_, u), -(1 + ap**2) / (1 + ap * u) ** 2)
assert egal(1 + ap * phi_, (1 + ap**2) / (1 + ap * u))
assert egal(diff(phi_, u) / (1 + phi_**2), -1 / (1 + u**2))
for a, val, leurres, ou in [
    (1, pi * log(2) / 8, [pi * log(2) / 4, pi * log(2) / 16, log(2) / 2], "cv-ex5"),
    (sqrt(3), pi / 3 * log(2), [pi / 6 * log(2), 2 * pi / 3 * log(2), pi / 4 * log(2)], "cv-ex5 v1"),  # leurre pi/4 ln2 (remplace pi/3 ln4 = 2pi/3 ln2, en double)
    (1 / sqrt(3), pi / 12 * log(Rational(4, 3)),
     [pi / 6 * log(Rational(4, 3)), pi / 24 * log(Rational(4, 3)), pi / 12 * log(Rational(3, 4))], "cv-ex5 v2"),
]:
    formule = atan(a) * log(1 + a**2) / 2
    bonne_et_leurres(formule, val, leurres, ou)
    # leurres deux à deux distincts (propositions de valeurs distinctes)
    for i_, l1 in enumerate(leurres):
        for l2 in leurres[i_ + 1:]:
            assert different(l1, l2), (ou, l1, l2)
    af = float(a)
    numv = num_integrale(lambda z: math.log(1 + af * z) / (1 + z * z), 0, af)
    assert abs(numv - float(val)) < 1e-9, (ou, numv)

# ===========================================================================
# Carte int-l1-riemann-nulle
# ===========================================================================
# Exemple : int_0^{pi/2} sin = 1 < pi^2/8.
assert integrate(sin(t), (t, 0, pi / 2)) == 1 and egal(integrate(t, (t, 0, pi / 2)), pi**2 / 8) and float(pi**2 / 8) > 1
# Contre-exemples : int_0^2 (1-t) = 0 ; int_0^{2pi} sin = 0 ; affirmation 1 : int_0^1 (2t-1) = 0.
assert integrate(1 - t, (t, 0, 2)) == 0 and integrate(sin(t), (t, 0, 2 * pi)) == 0 and integrate(2 * t - 1, (t, 0, 1)) == 0

# int-l1-riemann-nulle-ex1 : (a), (b) : 2x-1 ; (c) : 3x-2.
assert integrate(2 * x - 1, (x, 0, 1)) == 0 and integrate((2 * x - 1) ** 3, (x, 0, 1)) == 0
assert integrate(x * (3 * x - 2), (x, 0, 1)) == 0

# int-l1-riemann-nulle-ex2 : int f^2 = 2 int x f - 1/3 <=> int (f - x)^2 = 0.
def condition(fexpr, poids, cst):
    return simplify(integrate(fexpr**2, (x, 0, 1)) - 2 * integrate(poids * fexpr, (x, 0, 1)) + cst) == 0


assert integrate(x**2, (x, 0, 1)) == Rational(1, 3)
assert condition(x, x, Rational(1, 3))
for l in (2 * x, S(1), S(0)):  # leurres : 2x ; « toutes les affines » réfuté par f = 1 ; aucune réfuté par f = x
    assert not condition(l, x, Rational(1, 3))
# v1 : poids x^2, constante 1/5 -> f = x^2 ; leurres x, 2x^2.
assert integrate(x**4, (x, 0, 1)) == Rational(1, 5)
assert condition(x**2, x**2, Rational(1, 5))
for l in (x, 2 * x**2):
    assert not condition(l, x**2, Rational(1, 5))
# v2 : poids e^x, constante (e^2-1)/2 -> f = e^x ; leurres e^{2x}, 2e^x.
assert egal(integrate(exp(2 * x), (x, 0, 1)), (E**2 - 1) / 2)
assert condition(exp(x), exp(x), (E**2 - 1) / 2)
for l in (exp(2 * x), 2 * exp(x)):
    assert not condition(l, exp(x), (E**2 - 1) / 2)

# int-l1-riemann-nulle-ex3 : f^4 + c^2 f^2 = 2c f^3 <=> f^2 (f - c)^2 = 0, pour c = 1, 2, 3.
fs = symbols("fs")
for c, coef2, coef3 in [(1, 1, 2), (2, 4, 4), (3, 9, 6)]:
    assert egal(fs**4 + coef2 * fs**2 - coef3 * fs**3, fs**2 * (fs - c) ** 2)
    crit = lambda g: simplify(integrate(g**4, (x, 0, 1)) + coef2 * integrate(g**2, (x, 0, 1))
                              - coef3 * integrate(g**3, (x, 0, 1)))
    assert crit(S(0)) == 0 and crit(S(c)) == 0
    # leurres : « toutes les constantes » (c+1 ne convient pas), « valeurs dans [0,c] » (x ne convient pas),
    # « seule la constante c » (0 convient aussi).
    assert crit(S(c + 1)) != 0 and crit(x) != 0

# int-l1-riemann-nulle-ex4 : démonstration ; contrôle sur un exemple (n = 1) : f à moments d'ordre 0 et 1 nuls
# a bien deux zéros dans ]0,1[.
fex = x**2 - x + Rational(1, 6)  # int f = 1/3 - 1/2 + 1/6 = 0 ; int x f = 1/4 - 1/3 + 1/12 = 0
assert integrate(fex, (x, 0, 1)) == 0 and integrate(x * fex, (x, 0, 1)) == 0
from sympy import solve
racines = solve(fex, x)
assert len(racines) == 2 and all(0 < float(r) < 1 for r in racines)

# ===========================================================================
# Carte int-l1-riemann-fractions
# ===========================================================================
xg = symbols("xg", real=True)
pp, qq, al, be = symbols("pp qq al be", real=True)
D = 4 * qq - pp**2
# Formule de l'élément de seconde espèce (dérivée), et (2x+p)^2 + D = 4(x^2+px+q).
assert egal((2 * x + pp) ** 2 + D - 4 * (x**2 + pp * x + qq), 0)
Fprim = al / 2 * log(x**2 + pp * x + qq) + (be - al * pp / 2) * 2 / sqrt(D) * atan((2 * x + pp) / sqrt(D))
assert egal(diff(Fprim, x) - (al * x + be) / (x**2 + pp * x + qq), 0)
# Récurrence J_k : J_k = u/(2(k-1)(1+u^2)^{k-1}) + (2k-3)/(2(k-1)) J_{k-1}, contrôlée pour k = 2, 3, 4.
uu = symbols("uu", positive=True)
for kk_ in (2, 3, 4):
    Jk = integrate(1 / (1 + s**2) ** kk_, (s, 0, uu))
    Jk1 = integrate(1 / (1 + s**2) ** (kk_ - 1), (s, 0, uu))
    assert egal(Jk, uu / (2 * (kk_ - 1) * (1 + uu**2) ** (kk_ - 1)) + Rational(2 * kk_ - 3, 2 * (kk_ - 1)) * Jk1)
# Primitive de u/(1+u^2)^k.
for kk_ in (2, 3):
    assert egal(diff(-1 / (2 * (kk_ - 1) * (1 + u**2) ** (kk_ - 1)), u), u / (1 + u**2) ** kk_)
# Exemples : 1/(x^2-1) ; arctan(x+2) ; x^3/(x^2-1) = x + x/(x^2-1).
assert egal(diff(log((x - 1) / (x + 1)) / 2, x), 1 / (x**2 - 1))
assert egal(diff(atan(x + 2), x), 1 / (x**2 + 4 * x + 5)) and egal(integrate(1 / (x**2 + 4 * x + 5), (x, -2, -1)), pi / 4)
assert egal(apart(x**3 / (x**2 - 1), x), x + Rational(1, 2) / (x - 1) + Rational(1, 2) / (x + 1))
assert egal(Rational(1, 2) / (x - 1) + Rational(1, 2) / (x + 1), x / (x**2 - 1))
assert egal(diff(x**2 / 2 + log(x**2 - 1) / 2, x), x**3 / (x**2 - 1))
# Affirmations : ln(x^2+2x+2)' = 1 en 0 ; 1/(x^2+2x+2) = 1/2 en 0 ; primitive de 1/(1+x^2)^2.
assert diff(log(x**2 + 2 * x + 2), x).subs(x, 0) == 1 and (1 / (x**2 + 2 * x + 2)).subs(x, 0) == Rational(1, 2)
assert egal(diff((atan(x) + x / (1 + x**2)) / 2, x), 1 / (1 + x**2) ** 2)

# int-l1-riemann-fractions-ex1 : int_2^3 dx/(x^2-1) = (1/2) ln(3/2) ; leurres ln(3/2), (1/2)ln(2/3), (1/2)ln3.
bonne_et_leurres(integrate(1 / (x**2 - 1), (x, 2, 3)), log(Rational(3, 2)) / 2,
                 [log(Rational(3, 2)), log(Rational(2, 3)) / 2, log(3) / 2], "fr-ex1")
bonne_et_leurres(integrate(1 / (x**2 - 4), (x, 0, 1)), -log(3) / 4, [log(3) / 4, -log(3) / 2, -log(2) / 4], "fr-ex1 v1")
bonne_et_leurres(integrate(1 / (x * (x + 1)), (x, 1, 2)), log(Rational(4, 3)),
                 [log(Rational(3, 2)), log(2), log(Rational(3, 4))], "fr-ex1 v2")

# int-l1-riemann-fractions-ex2.
bonne_et_leurres(integrate(1 / (x**2 + x + 1), (x, 0, 1)), pi / (3 * sqrt(3)), [pi / 6, pi / (6 * sqrt(3)), pi / 3], "fr-ex2")
bonne_et_leurres(integrate(1 / (x**2 + 2 * x + 2), (x, 0, 1)), atan(2) - pi / 4,
                 [atan(2), log(Rational(5, 2)) / 2, atan(3) - pi / 4], "fr-ex2 v1")
bonne_et_leurres(integrate(1 / (x**2 + 2 * x + 5), (x, -1, 1)), pi / 8, [pi / 4, pi / 16, log(2) / 2], "fr-ex2 v2")

# int-l1-riemann-fractions-ex3.
assert egal(apart((x + 3) / ((x + 1) * (x**2 + 1)), x), 1 / (x + 1) + (2 - x) / (x**2 + 1))
bonne_et_leurres(integrate((x + 3) / ((x + 1) * (x**2 + 1)), (x, 0, 1)), log(2) / 2 + pi / 2,
                 [log(2) + pi / 2, log(2) / 2 + pi / 4, 3 * log(2) / 2 + pi / 2], "fr-ex3")
assert egal(apart(2 * x**2 / ((x + 1) * (x**2 + 1)), x), 1 / (x + 1) + (x - 1) / (x**2 + 1))
bonne_et_leurres(integrate(2 * x**2 / ((x + 1) * (x**2 + 1)), (x, 0, 1)), 3 * log(2) / 2 - pi / 4,
                 [3 * log(2) / 2 + pi / 4, log(2) / 2 - pi / 4, log(2) - pi / 4], "fr-ex3 v1")
assert egal(apart((3 * x + 1) / ((x + 1) * (x**2 + 1)), x), -1 / (x + 1) + (x + 2) / (x**2 + 1))
bonne_et_leurres(integrate((3 * x + 1) / ((x + 1) * (x**2 + 1)), (x, 0, 1)), pi / 2 - log(2) / 2,
                 [pi / 2 + log(2) / 2, pi / 4 - log(2) / 2, pi / 2 - 3 * log(2) / 2], "fr-ex3 v2")

# int-l1-riemann-fractions-ex4 : x^4 + 1 = (x^2 + sqrt2 x + 1)(x^2 - sqrt2 x + 1).
Qp, Qm = x**2 + sqrt(2) * x + 1, x**2 - sqrt(2) * x + 1
assert egal(Qp * Qm, x**4 + 1)
a4 = 1 / (2 * sqrt(2))
assert egal((a4 * x + Rational(1, 2)) / Qp + (-a4 * x + Rational(1, 2)) / Qm, 1 / (1 + x**4))
assert egal(a4 * x / Qm - a4 * x / Qp, x**2 / (1 + x**4))
assert egal(Rational(1, 2) / Qp + Rational(1, 2) / Qm, (1 + x**2) / (1 + x**4))
assert egal(diff(sqrt(2) * atan(sqrt(2) * x + 1), x), 1 / Qp) and egal(diff(sqrt(2) * atan(sqrt(2) * x - 1), x), 1 / Qm)
P0 = log(Qp / Qm) / (4 * sqrt(2)) + (atan(sqrt(2) * x + 1) + atan(sqrt(2) * x - 1)) / (2 * sqrt(2))
assert egal(diff(P0, x), 1 / (1 + x**4))
P1 = log(Qm / Qp) / (4 * sqrt(2)) + (atan(sqrt(2) * x - 1) + atan(sqrt(2) * x + 1)) / (2 * sqrt(2))
assert egal(diff(P1, x), x**2 / (1 + x**4))
P2 = (atan(sqrt(2) * x + 1) + atan(sqrt(2) * x - 1)) / sqrt(2)
assert egal(diff(P2, x), (1 + x**2) / (1 + x**4))
assert egal((2 + sqrt(2)) / (2 - sqrt(2)), 3 + 2 * sqrt(2)) and egal((1 + sqrt(2)) ** 2, 3 + 2 * sqrt(2))
assert abs(float(atan(sqrt(2) + 1) + atan(sqrt(2) - 1)) - math.pi / 2) < 1e-14
L = log(1 + sqrt(2))
for f4, prim, val, leurres, ou in [
    (lambda z: 1 / (1 + z**4), P0, (pi + 2 * L) / (4 * sqrt(2)), [pi / 4, (pi - 2 * L) / (4 * sqrt(2)), (pi + 2 * L) / (2 * sqrt(2))], "fr-ex4"),
    (lambda z: z**2 / (1 + z**4), P1, (pi - 2 * L) / (4 * sqrt(2)), [(pi + 2 * L) / (4 * sqrt(2)), pi / (4 * sqrt(2)), (pi - 2 * L) / (2 * sqrt(2))], "fr-ex4 v1"),
    (lambda z: (1 + z**2) / (1 + z**4), P2, pi / (2 * sqrt(2)), [pi / (4 * sqrt(2)), pi / 2, pi / sqrt(2)], "fr-ex4 v2"),
]:
    numv = num_integrale(f4, 0.0, 1.0)
    assert abs(numv - float(val)) < 1e-10, (ou, numv, float(val))
    assert abs(float(prim.subs(x, 1) - prim.subs(x, 0)) - float(val)) < 1e-12, ou
    for l in leurres:
        assert different(val, l), (ou, l)

# int-l1-riemann-fractions-ex5 : I_n(a) = int_0^a dx/(1+x^2)^n, I_{n+1} = a/(2n(1+a^2)^n) + (2n-1)/(2n) I_n.
for a, val, leurres, ou in [
    (S(1), Rational(1, 4) + 3 * pi / 32, [Rational(1, 4) + pi / 8, Rational(1, 8) + 3 * pi / 32, 3 * pi / 32], "fr-ex5"),
    (sqrt(3), 7 * sqrt(3) / 64 + pi / 8, [sqrt(3) / 8 + pi / 6, 7 * sqrt(3) / 64 + pi / 6, sqrt(3) / 64 + pi / 8], "fr-ex5 v1"),
    (1 / sqrt(3), 9 * sqrt(3) / 64 + pi / 16, [sqrt(3) / 8 + pi / 12, 9 * sqrt(3) / 64 + pi / 12, 3 * sqrt(3) / 64 + pi / 16], "fr-ex5 v2"),
]:
    I1, I2, I3 = (integrate(1 / (1 + x**2) ** m, (x, 0, a)) for m in (1, 2, 3))
    assert egal(I2, a / (2 * (1 + a**2)) + I1 / 2) and egal(I3, a / (4 * (1 + a**2) ** 2) + Rational(3, 4) * I2), ou
    bonne_et_leurres(I3, val, leurres, ou)
assert egal(integrate(1 / (1 + x**2) ** 2, (x, 0, 1)), Rational(1, 4) + pi / 8)
assert egal(integrate(1 / (1 + x**2) ** 2, (x, 0, sqrt(3))), sqrt(3) / 8 + pi / 6)
assert egal(integrate(1 / (1 + x**2) ** 2, (x, 0, 1 / sqrt(3))), sqrt(3) / 8 + pi / 12)
assert egal(9 / (64 * sqrt(3)), 3 * sqrt(3) / 64)

# ===========================================================================
# Carte int-l1-riemann-cauchy-schwarz
# ===========================================================================
# Exemple : (int sqrt t e^t)^2 <= (1/2)(e^2-1)/2 ; contrôle numérique.
assert egal(integrate(t, (t, 0, 1)) * integrate(exp(2 * t), (t, 0, 1)), (E**2 - 1) / 4)
assert num_integrale(lambda z: math.sqrt(z) * math.exp(z), 0, 1) ** 2 <= float((E**2 - 1) / 4)
# Contre-exemple : sur [0,2], f = 1 : 4 > 2.
assert integrate(1, (t, 0, 2)) ** 2 == 4 and integrate(1, (t, 0, 2)) == 2

# int-l1-riemann-cauchy-schwarz-ex1 : sqrt(5)/2 ; leurres 5/4, sqrt5/4, sqrt(5/2).
for poly, val, leurres, ou in [
    (1 + x**3, sqrt(5) / 2, [Rational(5, 4), sqrt(5) / 4, sqrt(Rational(5, 2))], "cs-ex1"),
    (1 + x**2, 2 / sqrt(3), [Rational(4, 3), 1 / sqrt(3), sqrt(Rational(8, 3))], "cs-ex1 v1"),
    (1 + x**5, sqrt(Rational(7, 6)), [Rational(7, 6), sqrt(Rational(7, 12)), sqrt(Rational(6, 5))], "cs-ex1 v2"),
]:
    bonne_et_leurres(sqrt(integrate(poly, (x, 0, 1))), val, leurres, ou)
    # la vraie intégrale est bien inférieure au majorant
    assert num_integrale(lambda z, poly=poly: math.sqrt(float(poly.subs(x, z))), 0, 1) <= float(val)

# int-l1-riemann-cauchy-schwarz-ex2 : (b-a)^2.
for a, b, val, leurres, ou in [(1, 3, 4, [2, 9, 1], "cs-ex2"), (0, 5, 25, [5, 10, 1], "cs-ex2 v1"),
                               (-1, 2, 9, [3, 4, 1], "cs-ex2 v2")]:
    bonne_et_leurres(Integer(b - a) ** 2, val, leurres, ou)
    # f = 1 réalise la valeur ; f = 1 + x^2 donne une valeur plus grande.
    assert integrate(1, (x, a, b)) * integrate(1, (x, a, b)) == val
    assert float(integrate(1 + x**2, (x, a, b)) * integrate(1 / (1 + x**2), (x, a, b))) > val

# int-l1-riemann-cauchy-schwarz-ex3 : max = sqrt(int g^2), atteint en g / ||g||.
for g, val, leurres, ou in [
    (x, 1 / sqrt(3), [Rational(1, 3), Rational(1, 2), sqrt(3)], "cs-ex3"),
    (x**2, 1 / sqrt(5), [Rational(1, 5), Rational(1, 3), sqrt(5)], "cs-ex3 v1"),
    (exp(x), sqrt((E**2 - 1) / 2), [(E**2 - 1) / 2, E - 1, sqrt(E**2 - 1)], "cs-ex3 v2"),
]:
    ng = sqrt(integrate(g**2, (x, 0, 1)))
    bonne_et_leurres(ng, val, leurres, ou)
    fopt = g / ng
    assert egal(integrate(fopt**2, (x, 0, 1)), 1) and egal(integrate(g * fopt, (x, 0, 1)), val)
assert egal(sqrt(3) / 3, 1 / sqrt(3)) and egal(sqrt(5) / 5, 1 / sqrt(5))

# int-l1-riemann-cauchy-schwarz-ex4 : min int_0^L f'^2 = (f(L)-f(0))^2 / L, atteint pour f affine.
for L_, f0, fL, val, leurres, ou in [(1, 0, 1, 1, [Rational(1, 2), 2], "cs-ex4"),
                                     (2, 1, 5, 8, [16, 4, 2], "cs-ex4 v1"),
                                     (3, 0, 2, Rational(4, 3), [Rational(4, 9), Rational(2, 3), 4], "cs-ex4 v2")]:
    bonne_et_leurres(Rational((fL - f0) ** 2, L_), val, leurres, ou)
    faff = f0 + Rational(fL - f0, L_) * x
    assert egal(integrate(diff(faff, x) ** 2, (x, 0, L_)), val)
    # une autre fonction admissible fait plus : f = faff + sin(pi x / L).
    autre = faff + sin(pi * x / L_)
    assert autre.subs(x, 0) == f0 and egal(autre.subs(x, L_), fL)
    assert float(integrate(diff(autre, x) ** 2, (x, 0, L_))) > float(val)
assert egal(1 + 2 * x - (1 + Rational(4, 2) * x), 0)

# int-l1-riemann-cauchy-schwarz-ex5 : C = L^2/8 ; rapport pour x(L-x) = L^2/10.
for L_, val, leurres, ou in [(1, Rational(1, 8), [Rational(1, 2), Rational(1, 4), Rational(1, 10)], "cs-ex5"),
                             (2, Rational(1, 2), [2, 1, Rational(2, 5)], "cs-ex5 v1"),
                             (3, Rational(9, 8), [Rational(9, 2), Rational(9, 4), Rational(9, 10)], "cs-ex5 v2")]:
    assert egal(integrate(x, (x, 0, Rational(L_, 2))), Rational(L_**2, 8))
    assert egal(integrate(L_ - x, (x, Rational(L_, 2), L_)), Rational(L_**2, 8))
    bonne_et_leurres(Rational(L_**2, 8), val, leurres, ou)
    g = x * (L_ - x)
    num_, den_ = integrate(g**2, (x, 0, L_)), integrate(diff(g, x) ** 2, (x, 0, L_))
    assert egal(num_ / den_, Rational(L_**2, 10)) and num_ / den_ <= val
assert integrate((x * (1 - x)) ** 2, (x, 0, 1)) == Rational(1, 30) and integrate((1 - 2 * x) ** 2, (x, 0, 1)) == Rational(1, 3)
assert integrate((x * (2 - x)) ** 2, (x, 0, 2)) == Rational(16, 15) and integrate((2 - 2 * x) ** 2, (x, 0, 2)) == Rational(8, 3)
assert integrate((x * (3 - x)) ** 2, (x, 0, 3)) == Rational(81, 10) and integrate((3 - 2 * x) ** 2, (x, 0, 3)) == 9

print("int_l1_riemann : toutes les vérifications passent.")
