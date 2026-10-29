"""Vérifications calculatoires du chapitre « int-l2-generalisees » (intégration, L2) :
cartes, exercices et QCM de content/*/integration/l2-int-l2-generalisees.json et content/qcm/int-l2-generalisees.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Calculs exacts avec SymPy quand ils sont rapides, sinon numériques à 30 chiffres (mpmath).
Lancer : python3 verification/int_l2_generalisees.py (doit finir sans erreur, en quelques secondes).
"""

from fractions import Fraction as F

import mpmath as mp
from sympy import (E, Rational, atan, diff, exp, integrate, limit, log, oo, pi, series, simplify, sin, cos, sqrt,
                   symbols)

mp.mp.dps = 30
t, x, u, s = symbols("t x u s", positive=True)
n = symbols("n", positive=True, integer=True)


def mpf(v):
    return mp.mpf(v.numerator) / v.denominator if isinstance(v, F) else mp.mpf(v)


def proche(a, b, tol=mp.mpf("1e-20")):
    return abs(mpf(a) - mpf(b)) < tol


def Q(f, a, b):
    return mp.quad(f, [a, b])


def Qinf(f, a=0):
    return mp.quad(f, [a, a + 1, a + 10, mp.inf])


def distincts(vals):
    """Les valeurs numériques des propositions sont deux à deux distinctes (leurres faux)."""
    vals = [mp.mpf(v.numerator) / v.denominator if isinstance(v, F) else mp.mpf(v) for v in vals]
    return all(abs(vals[i] - vals[j]) > mp.mpf("1e-10") for i in range(len(vals)) for j in range(i))


# ---------------------------------------------------------------------------
# Carte int-l2-generalisees-positive
# ---------------------------------------------------------------------------
# Exemple : Riemann sur [1,+oo[ et ]0,1], exponentielle.
for a in [Rational(3, 2), 2, 3]:
    assert simplify(integrate(t**(-a), (t, 1, x)) - (x**(1 - a) - 1) / (1 - a)) == 0
    assert integrate(t**(-a), (t, 1, oo)) == 1 / (Rational(a) - 1)
for a in [Rational(1, 2), Rational(2, 3)]:
    assert integrate(t**(-a), (t, 0, 1)) == 1 / (1 - a)
assert integrate(exp(-3 * t), (t, 0, oo)) == Rational(1, 3)
# Contre-exemple : int_0^{n pi} |cos| = 2n.
assert proche(mp.quad(lambda v: abs(mp.cos(v)), [0] + [mp.pi / 2 + j * mp.pi for j in range(5)] + [5 * mp.pi]), 10, mp.mpf("1e-15"))
# Affirmation 1 : pics de hauteur 1 sur [n - 2^-n, n + 2^-n], aire 2^-n, total 1 ; intervalles disjoints.
assert sum(F(1, 2**k) for k in range(1, 200)) == 1 - F(1, 2**199)
assert all(F(1, 2**k) + F(1, 2**(k + 1)) < 1 for k in range(1, 50))

# positive-ex1 (+ v1, v2) : 1/2, 1/3, 1/8 ; leurres distincts.
assert integrate(t**-3, (t, 1, oo)) == Rational(1, 2)
assert integrate(t**-4, (t, 1, oo)) == Rational(1, 3)
assert integrate(t**-3, (t, 2, oo)) == Rational(1, 8)
assert distincts([F(1, 2), F(1, 3), F(1, 4), 1]) and distincts([F(1, 3), F(1, 4), F(1, 2), 1]) and distincts([F(1, 8), F(1, 4), F(1, 16), F(1, 2)])
# positive-ex2 (+ v1, v2) : primitives et valeurs 1, 1/4, 4.
assert simplify(diff(t - t * log(t), t) + log(t)) == 0
assert limit(1 - x + x * log(x), x, 0) == 1
assert simplify(diff(-t**2 / 2 * log(t) + t**2 / 4, t) + t * log(t)) == 0
assert limit(Rational(1, 4) + x**2 / 2 * log(x) - x**2 / 4, x, 0) == Rational(1, 4)
G = -2 * sqrt(t) * log(t) + 4 * sqrt(t)
assert simplify(diff(G, t) + log(t) / sqrt(t)) == 0
assert limit(4 - G.subs(t, x), x, 0) == 4
assert proche(Q(lambda v: -mp.log(v) / mp.sqrt(v), 0, 1), 4, mp.mpf("1e-6"))
# positive-ex3 (+ v1, v2) : ln 2, ln(3)/2, 1 - ln 2 ; décompositions.
assert simplify(1 / ((1 + t) * (2 + t)) - (1 / (1 + t) - 1 / (2 + t))) == 0
assert simplify(1 / ((1 + t) * (3 + t)) - (1 / (1 + t) - 1 / (3 + t)) / 2) == 0
assert simplify(1 / (t**2 * (t + 1)) - (1 / t**2 - 1 / t + 1 / (t + 1))) == 0
assert proche(Qinf(lambda v: 1 / ((1 + v) * (2 + v))), mp.log(2))
assert proche(Qinf(lambda v: 1 / ((1 + v) * (3 + v))), mp.log(3) / 2)
assert proche(Qinf(lambda v: 1 / (v**2 * (v + 1)), 1), 1 - mp.log(2))
assert distincts([mp.log(2), 0.5, mp.log(3), 1]) and distincts([mp.log(3) / 2, mp.log(3), mp.log(2) / 2, mp.mpf(1) / 3])
assert distincts([1 - mp.log(2), mp.log(2), 1 + mp.log(2), 0.5])
# positive-ex4 (+ v1, v2) : sommes des aires 1, 4/9, 1/2.
assert sum(F(k, 4**k) for k in range(1, 300)) - F(4, 9) < F(1, 10**100)
assert F(1, 4) / (1 - F(1, 4))**2 == F(4, 9)
assert F(1, 3) / (1 - F(1, 3)) == F(1, 2)
assert distincts([F(4, 9), F(1, 3), F(4, 3)]) and distincts([F(1, 2), 1, F(1, 3)])

# ---------------------------------------------------------------------------
# Carte int-l2-generalisees-comparaison
# ---------------------------------------------------------------------------
# Affirmation 1 : primitive de 1/((t+2) ln(t+2)) = ln ln(t+2) ; t f(t) -> 0.
assert simplify(diff(log(log(t + 2)), t) - 1 / ((t + 2) * log(t + 2))) == 0
assert limit(t / ((t + 2) * log(t + 2)), t, oo) == 0
# Affirmation 4 : int e^-t = 1, int (e^-t + e^-2t) = 3/2.
assert integrate(exp(-t) + exp(-2 * t), (t, 0, oo)) == Rational(3, 2)
# Exemples : t^4 e^-t -> 0, sqrt(t) ln t -> 0, ln t / sqrt t -> 0.
assert limit(t**4 * exp(-t), t, oo) == 0 and limit(sqrt(t) * log(t), t, 0) == 0 and limit(log(t) / sqrt(t), t, oo) == 0

# comparaison-ex1 et variantes : seules t^2 e^-t, ln t (sur ]0,1]) et ln t / t^2 sont intégrables.
assert integrate(t**2 * exp(-t), (t, 1, oo)) == 5 * exp(-1)
assert integrate(log(t), (t, 0, 1)) == -1
assert integrate(log(t) / t**2, (t, 1, oo)) == 1
# les autres sont équivalentes à 1/t ou minorées par 1/t (exposant 1/2 pour 1/sqrt t) :
for g in [1 / (t + t**2), 1 / (t + sqrt(t)), t / (t**2 + 1)]:
    assert limit(g * t, t, 0 if g == 1 / (t + t**2) else oo) == 1
assert limit(t / sin(t), t, 0) == 1

# comparaison-ex2 (+ v1, v2) : exposants des équivalents aux bornes.
# en +oo, t^alpha/(1+t^3) ~ t^(alpha-3), t^alpha/(1+t^4) ~ t^(alpha-4), t^alpha/sqrt(1+t^2) ~ t^(alpha-1) :
assert limit(t**3 / (1 + t**3), t, oo) == 1 and limit(t**4 / (1 + t**4), t, oo) == 1 and limit(t / sqrt(1 + t**2), t, oo) == 1
# point intérieur : alpha = 1/2 (ex2), alpha = 3/2 (ex3), alpha = -1/2 (ex2 v2) donnent des intégrales finies.
assert integrate(sqrt(t) / (1 + t**3), (t, 0, oo)).is_finite
assert mp.isfinite(Qinf(lambda v: (1 - mp.exp(-v)) / v**1.5))
assert mp.isfinite(Qinf(lambda v: v**-0.5 / mp.sqrt(1 + v**2)))
# comparaison-ex3 (+ v1, v2) : équivalents en 0.
assert limit((1 - exp(-t)) / t, t, 0) == 1 and limit((1 - exp(-t**2)) / t**2, t, 0) == 1
# comparaison-ex4 : changement de variable u = ln t.
assert simplify(integrate(1 / (t * log(t)**2), (t, E, x)) - (1 - 1 / log(x))) == 0
assert simplify(diff(log(log(t)), t) - 1 / (t * log(t))) == 0

# ---------------------------------------------------------------------------
# Carte int-l2-generalisees-cauchy
# ---------------------------------------------------------------------------
assert integrate(1 / sqrt(1 - t), (t, 0, 1)) == 2          # exemple, cauchy-ex1
assert limit(log(1 + 1 / x), x, oo) == 0
assert simplify(integrate(1 / t, (t, x, 2 * x)) - log(2)) == 0     # contre-exemple
# Affirmation 1 : sin(t^2) = 1 en sqrt(pi/2 + 2k pi), 0 en sqrt(k pi).
assert all(abs(mp.sin(mp.sqrt(mp.pi / 2 + 2 * k * mp.pi)**2) - 1) < 1e-20 for k in range(5))
# cauchy-ex1 v1, v2 : 3/2 et 1.
assert integrate((1 - t)**Rational(-1, 3), (t, 0, 1)) == Rational(3, 2)
assert integrate(t / sqrt(1 - t**2), (t, 0, 1)) == 1
# cauchy-ex2 : pi/2 - arctan t = arctan(1/t) ~ 1/t (divergence) ; v1 : 1/(t^2+1) -> pi/4 ; v2 : 1/(e^t+1) -> ln 2.
assert limit((pi / 2 - atan(t)) * t, t, oo) == 1
assert simplify(1 - t**2 / (t**2 + 1) - 1 / (t**2 + 1)) == 0
assert integrate(1 / (t**2 + 1), (t, 1, oo)) == pi / 4
assert proche(Qinf(lambda v: 1 / (mp.exp(v) + 1)), mp.log(2))
# cauchy-ex3 (+ v1, v2) : int_{n pi}^{2n pi} |sin| = 2n ; minorations ; les leurres dépassent la valeur pour n grand.
for k in [1, 2, 5]:
    assert proche(mp.quad(lambda v: abs(mp.sin(v)), [j * mp.pi for j in range(k, 2 * k + 1)]), 2 * k, mp.mpf("1e-15"))
def I_arches(p, q, N, puis=1):
    return sum(Q(lambda v: abs(mp.sin(v)) / v**puis, j * mp.pi, (j + 1) * mp.pi) for j in range(p * N, q * N))
for N in [1, 3]:
    assert I_arches(1, 2, N) >= 1 / mp.pi and I_arches(1, 3, N) >= 4 / (3 * mp.pi)
    assert I_arches(1, 2, N, 0.5) >= mp.sqrt(2 * N / mp.pi)
    assert all(I_arches(1, 2, N, 0.5) < L for L in [2 * mp.sqrt(N / mp.pi), mp.sqrt(2 * N * mp.pi), 2 * N / mp.sqrt(mp.pi)])
# limites (2/pi) ln 2 et (2/pi) ln 3 : les leurres 2/pi, 4/pi, 1 et 4/pi, 8/(3pi), 1 sont trop grands pour n grand.
assert proche(I_arches(1, 2, 40), 2 / mp.pi * mp.log(2), mp.mpf("1e-2"))
assert 2 / mp.pi * mp.log(2) < 2 / mp.pi < 1 < 4 / mp.pi
assert I_arches(1, 3, 40) < 8 / (3 * mp.pi) and 2 / mp.pi * mp.log(3) < 8 / (3 * mp.pi)
# cauchy-ex4 : contre-exemples des leurres.
assert limit(t**2 / (1 + t)**Rational(3, 2), t, oo) == oo          # pas O(1/x^2)
assert integrate(exp(-t), (t, x, oo)) == exp(-x)                    # pas o(f(x))

# ---------------------------------------------------------------------------
# Carte int-l2-generalisees-absolue
# ---------------------------------------------------------------------------
# Trou 2, leurres : f = e^-t - 2 e^-2t d'intégrale nulle mais non nulle ; 2 e^{-t/2} d'intégrale 4 > sup = 2.
assert integrate(exp(-t) - 2 * exp(-2 * t), (t, 0, oo)) == 0
assert integrate(2 * exp(-t / 2), (t, 0, oo)) == 4
# Démonstration : identité d'intégration par parties.
assert simplify(diff(-cos(t) / t, t) - (sin(t) / t + cos(t) / t**2)) == 0
# Exemple : |sin| >= sin^2 = (1 - cos 2t)/2.
assert simplify(sin(t)**2 - (1 - cos(2 * t)) / 2) == 0
# absolue-ex3 : int (sin t / t)^2 = pi/2 ; v1 : sin^3 linéarisé, pi/4 ; v2 : sin t cos 2t = (sin 3t - sin t)/2, 0.
assert simplify(diff(-sin(t)**2 / t, t) - sin(t)**2 / t**2 + sin(2 * t) / t) == 0
assert integrate(sin(t)**2 / t**2, (t, 0, oo)) == pi / 2 and integrate(sin(t)**3 / t, (t, 0, oo)) == pi / 4
assert simplify(sin(t)**3 - (3 * sin(t) - sin(3 * t)) / 4) == 0
assert simplify(sin(t) * cos(2 * t) - (sin(3 * t) - sin(t)) / 2) == 0
assert Rational(3, 4) * pi / 2 - Rational(1, 4) * pi / 2 == pi / 4
assert distincts([mp.pi / 4, 3 * mp.pi / 8, mp.pi / 2, mp.pi / 8]) and distincts([0, mp.pi / 4, -mp.pi / 4, mp.pi / 2])
# absolue-ex4 : changement de variable et intégration par parties.
assert simplify(diff(-cos(u) / (2 * sqrt(u)), u) - (sin(u) / (2 * sqrt(u)) + cos(u) / (4 * u**Rational(3, 2)))) == 0
# absolue-ex5 : télescopage, développements de phi et phi'.
k = symbols("k", positive=True)
assert simplify((2 * sin(t / 2) * cos(k * t) - (sin((k + Rational(1, 2)) * t) - sin((k - Rational(1, 2)) * t))).rewrite(exp)) == 0
phi = 1 / t - 1 / (2 * sin(t / 2))
assert limit(phi / t, t, 0) == Rational(-1, 24)
assert limit(diff(phi, t), t, 0) == Rational(-1, 24)
th = mp.mpf("1.3")
assert proche(sum(mp.cos(j * th) for j in range(1, 8)) + mp.mpf(1) / 2, mp.sin(mp.mpf("7.5") * th) / (2 * mp.sin(th / 2)), mp.mpf("1e-20"))

# ---------------------------------------------------------------------------
# Carte int-l2-generalisees-relations
# ---------------------------------------------------------------------------
# Contre-exemple 1 : sommes alternées de ln((k+2)/(k+1)) sur 2m termes > 0.
assert all(sum((-1)**j * mp.log(mp.mpf(j + 2) / (j + 1)) for j in range(2 * m)) > 0 for m in range(1, 30))
# Trou 1, leurres « bornée », « continue », « de limite nulle » : int_0^oo cos t/(1+t) et cos t/(1+t)^2 non nuls.
L1 = mp.quadosc(lambda v: mp.cos(v) / (1 + v), [0, mp.inf], omega=1)
L2 = mp.quadosc(lambda v: mp.cos(v) / (1 + v)**2, [0, mp.inf], omega=1)
assert L1 > 0.3 and L2 > 0.3
# relations-ex1 (+ v1, v2) : équivalents des restes ; v2 forme exacte.
assert limit(t**3 / (t**3 + t + 1), t, oo) == 1 and integrate(t**-3, (t, x, oo)) == 1 / (2 * x**2)
assert integrate(t**-4, (t, x, oo)) == 1 / (3 * x**3)
R2 = log((1 + sqrt(1 + x**2)) / x)
assert simplify(diff(R2, x) + 1 / (x * sqrt(x**2 + 1))) == 0 and limit(R2 * x, x, oo) == 1
# relations-ex2 (+ v1, v2) : x^3/3, x^2/2, pi x^2/4.
assert limit(sqrt(1 + t**4) / t**2, t, oo) == 1 and limit(sqrt(t**2 + t) / t, t, oo) == 1
assert limit(t * atan(t) / t, t, oo) == pi / 2
assert proche(Q(lambda v: v * mp.atan(v), 1, 1000) / (mp.pi * 1000**2 / 4), 1, mp.mpf("1e-2"))
# relations-ex3 (+ v1, v2) : formes exactes.
assert simplify(integrate(atan(t), (t, 0, x)) - (x * atan(x) - log(1 + x**2) / 2)) == 0
assert limit((x * atan(x) - log(1 + x**2) / 2) / x, x, oo) == pi / 2
assert limit(integrate(t**2 / (t**2 + 1), (t, 0, x)) / x, x, oo) == 1
assert limit(integrate((2 * t + 1) / (t + 3), (t, 0, x)) / x, x, oo) == 2
# relations-ex4 (+ v1, v2) : intégrations par parties et maximum de (ln t)^2 / t.
assert simplify(diff(t / log(t), t) - (1 / log(t) - 1 / log(t)**2)) == 0
assert simplify(diff(t / log(t)**2, t) - (1 / log(t)**2 - 2 / log(t)**3)) == 0
assert simplify(diff(exp(t) / t, t) - (exp(t) / t - exp(t) / t**2)) == 0
assert (log(E**2)**2 / E**2).evalf() < 1
assert proche(Q(lambda v: 1 / mp.log(v), 2, 1e6) / (1e6 / mp.log(1e6)), 1.08, mp.mpf("0.02"))
# relations-ex5 (+ v1, v2) : intégrations par parties.
assert simplify(diff(-exp(-t**2) / (2 * t), t) - (exp(-t**2) + exp(-t**2) / (2 * t**2))) == 0
assert simplify(diff(-exp(-t) / t, t) - (exp(-t) / t + exp(-t) / t**2)) == 0
assert simplify(diff(-exp(-t**3) / (3 * t**2), t) - (exp(-t**3) + 2 * exp(-t**3) / (3 * t**3))) == 0
X = mp.mpf(6)
assert proche(Qinf(lambda v: mp.exp(-v**2), X) / (mp.exp(-X**2) / (2 * X)), 1, mp.mpf("0.02"))

# ---------------------------------------------------------------------------
# Carte int-l2-generalisees-rectangles
# ---------------------------------------------------------------------------
# Contre-exemples : f = -t donne w_n = 1/2 ; int_0^N |sin(pi t)| = 2N/pi ; pics à partir de n = 2 : aire 1/2.
assert simplify(-n - integrate(-t, (t, n, n + 1)) - Rational(1, 2)) == 0
assert proche(mp.quad(lambda v: abs(mp.sin(mp.pi * v)), [0, 1, 2, 3]), 6 / mp.pi, mp.mpf("1e-15"))
assert sum(F(1, 2**j) for j in range(2, 200)) < F(1, 2)
# Exemple : ln(N+1) <= H_N <= 1 + ln N.
H = lambda N: mp.fsum(1 / mp.mpf(j) for j in range(1, N + 1))
assert all(mp.log(N + 1) <= H(N) <= 1 + mp.log(N) for N in range(1, 200))
# rectangles-ex1 : leurres faux pour un n explicite.
assert H(1) > mp.log(2) and H(2) > mp.log(2) + mp.mpf(1) / 2 and mp.log(2) + mp.mpf(1) / 2 > H(1)
S = lambda N: sum(1 / mp.sqrt(j) for j in range(1, N + 1))
assert all(2 * mp.sqrt(N + 1) - 2 <= S(N) <= 2 * mp.sqrt(N) - 1 for N in range(1, 200))
assert S(1) < 2 and S(2) > mp.sqrt(2) and S(1) > 0
T = lambda N: sum(F(1, j * j) for j in range(1, N + 1))
assert all(F(N, N + 1) <= T(N) <= 2 - F(1, N) for N in range(1, 200))
assert T(2) > 1 and T(3) > 1 + F(1, 3)
# rectangles-ex2 (+ v1, v2) : primitives.
assert simplify(diff(-1 / log(t), t) - 1 / (t * log(t)**2)) == 0
assert simplify(diff(2 * sqrt(log(t)), t) - 1 / (t * sqrt(log(t)))) == 0
# rectangles-ex3 (+ v1, v2) : bornes des restes.
for N in [1, 5, 20]:
    R = mp.zeta(2, N + 1)
    assert mp.mpf(1) / (N + 1) <= R <= mp.mpf(1) / N
    R = mp.zeta(3, N + 1)
    assert 1 / (2 * mp.mpf(N + 1)**2) <= R <= 1 / (2 * mp.mpf(N)**2)
    R = mp.zeta(mp.mpf(3) / 2, N + 1)
    assert 2 / mp.sqrt(N + 1) <= R <= 2 / mp.sqrt(N)
# rectangles-ex4 (+ v1, v2) : gamma dans [0,1] ; limites ln 2, ln 3, ln 2.
assert 0 <= mp.euler <= 1
N = 10**4
assert proche(sum(mp.mpf(1) / j for j in range(N + 1, 2 * N + 1)), mp.log(2), mp.mpf("1e-4"))
assert proche(sum(mp.mpf(1) / j for j in range(N + 1, 3 * N + 1)), mp.log(3), mp.mpf("1e-4"))
assert proche(sum(mp.mpf((-1)**(j + 1)) / j for j in range(1, 2 * N + 1)), mp.log(2), mp.mpf("1e-4"))

# ---------------------------------------------------------------------------
# Carte int-l2-generalisees-dominee
# ---------------------------------------------------------------------------
# Exemple : int_0^n (1 - t/n)^n = n/(n+1) ; min(n, 1/sqrt t) : 2 - 1/n.
assert simplify(integrate((1 - t / 5)**5, (t, 0, 5)) - Rational(5, 6)) == 0
assert integrate(sqrt(t)**-1, (t, Rational(1, 25), 1)) + 5 * Rational(1, 25) == 2 - Rational(1, 5)
# Affirmation 2 : int_0^oo dt/(1+t^n) -> 1.
assert proche(Qinf(lambda v: 1 / (1 + v**60)), 1, mp.mpf("0.03"))
# dominee-ex1 (+ v1, v2) : n^2/(n^2+1) -> 1 ; 1/4 ; 4.
assert simplify(integrate(n * sin(t / n) * exp(-t), (t, 0, oo)) - n**2 / (n**2 + 1)) == 0
assert integrate(t * exp(-2 * t), (t, 0, oo)) == Rational(1, 4) and integrate(t * exp(-t / 2), (t, 0, oo)) == 4
# dominee-ex2 (+ v1, v2) : 1, 3/2, 2 ; forme exacte de v2.
assert integrate(exp(-2 * t / 3), (t, 0, oo)) == Rational(3, 2) and integrate(t**2 * exp(-t), (t, 0, oo)) == 2
assert simplify(integrate((1 - t / 7)**7 * t**2, (t, 0, 7)) - Rational(2 * 7**3, 8 * 9 * 10)) == 0
assert proche(Q(lambda v: (1 - v / 200)**200 * mp.exp(v / 3), 0, 200), 1.5, mp.mpf("0.02"))
# dominee-ex3 (+ v1, v2) : intégrales 1, n/(n+1), 2 ; minoration de sup_n f_n.
assert integrate(n**2 * t * exp(-n * t), (t, 0, oo)) == 1
assert integrate(n * t**n, (t, 0, 1)) == n / (n + 1)
assert integrate(n**3 * t**2 * exp(-n * t), (t, 0, oo)) == 2
for tt in [mp.mpf("0.37"), mp.mpf("0.05"), mp.mpf("0.9")]:
    m = int(mp.ceil(1 / tt))
    assert m**2 * tt * mp.exp(-m * tt) >= mp.exp(-2) / tt
    ss = 1 - (tt if tt >= 0.5 else mp.mpf("0.6"))
    m = int(mp.ceil(1 / ss))
    assert m * (1 - ss)**m >= mp.exp(-4) / ss
assert all(mp.log(1 - v) >= -2 * v for v in [mp.mpf(j) / 100 for j in range(51)])
# dominee-ex4 : f = cos : I_n = n^2/(n^2+1).
assert simplify(integrate(n * cos(t) * exp(-n * t), (t, 0, oo)) - n**2 / (n**2 + 1)) == 0
# QCM q7 : int n/(1+n^2 t^2) = pi/2.
assert integrate(n / (1 + n**2 * t**2), (t, 0, oo)) == pi / 2

# ---------------------------------------------------------------------------
# Carte int-l2-generalisees-terme-a-terme
# ---------------------------------------------------------------------------
# Contre-exemple 1 et terme-a-terme-ex3 : int n e^{-nt} = 1 ; int |u_n| = 2/(n+1) (n/(n+1))^n ~ 2/(e n).
assert integrate(n * exp(-n * t), (t, 0, oo)) == 1
for m in [1, 2, 7]:
    un = lambda v, m=m: (m + 1) * mp.exp(-(m + 1) * v) - m * mp.exp(-m * v)
    tn = mp.log(mp.mpf(m + 1) / m)
    assert proche(Q(lambda v: abs(un(v)), 0, tn) + Qinf(lambda v: abs(un(v)), tn), 2 / mp.mpf(m + 1) * (mp.mpf(m) / (m + 1))**m, mp.mpf("1e-15"))
assert limit(2 * n / (n + 1) * (n / (n + 1))**n, n, oo) == 2 / E
# Contre-exemple 2 : la série harmonique diverge (S = 1/(floor t + 1) non intégrable).
# Exemple et terme-a-terme-ex2 : int t e^{-nt} = 1/n^2 ; int t/(e^t-1) = pi^2/6 ; v1 2 zeta(3) ; v2 pi^2/12.
assert integrate(t * exp(-n * t), (t, 0, oo)) == 1 / n**2 and integrate(t**2 * exp(-n * t), (t, 0, oo)) == 2 / n**3
assert proche(Qinf(lambda v: v / mp.expm1(v)), mp.pi**2 / 6)
assert proche(Qinf(lambda v: v**2 / mp.expm1(v)), 2 * mp.zeta(3))
assert proche(Qinf(lambda v: v / (mp.exp(v) + 1)), mp.pi**2 / 12)
assert distincts([2 * mp.zeta(3), mp.zeta(3), 3 * mp.zeta(3), mp.pi**2 / 3])
assert distincts([mp.pi**2 / 12, mp.pi**2 / 6, mp.log(2), mp.pi**2 / 24]) and distincts([mp.pi**2 / 6, mp.pi**2 / 12, 1, mp.pi / 2])
# terme-a-terme-ex1 (+ v1, v2) : pi^2/6, pi^2/12, pi^2/8 ; int -t^a ln t = 1/(a+1)^2.
assert proche(Q(lambda v: -mp.log(1 - v) / v, 0, 1), mp.pi**2 / 6)
assert proche(Q(lambda v: mp.log(1 + v) / v, 0, 1), mp.pi**2 / 12)
assert proche(Q(lambda v: -mp.log(v) / (1 - v**2), 0, 1), mp.pi**2 / 8)
a = symbols("a", positive=True)
assert simplify(integrate(-t**a * log(t), (t, 0, 1)) - 1 / (a + 1)**2) == 0
assert distincts([mp.pi**2 / 8, mp.pi**2 / 6, mp.pi**2 / 12, mp.pi**2 / 4])
# terme-a-terme-ex4 (+ v1, v2) : int (t ln t)^n = (-1)^n n!/(n+1)^(n+1) ; rêve du bachelier.
for m in range(5):
    assert proche(Q(lambda v: (v * mp.log(v))**m, 0, 1), (-1)**m * mp.factorial(m) / mp.mpf(m + 1)**(m + 1))
assert proche(Q(lambda v: v**-v, 0, 1), mp.nsum(lambda j: j**-j, [1, mp.inf]))
assert proche(Q(lambda v: v**v, 0, 1), mp.nsum(lambda j: (-1)**(j + 1) * j**-j, [1, mp.inf]))
assert proche(Q(lambda v: v**(2 * v), 0, 1), mp.nsum(lambda j: (-2)**(j - 1) * j**-j, [1, mp.inf]))
# leurres : valeurs distinctes.
sk = lambda f: mp.nsum(f, [1, mp.inf])
assert distincts([sk(lambda j: j**-j), sk(lambda j: (-1)**(j + 1) * j**-j), mp.e, sk(lambda j: 1 / (j * mp.factorial(j)))])
assert distincts([sk(lambda j: (-1)**(j + 1) * j**-j), sk(lambda j: j**-j), 1 / mp.e, sk(lambda j: (-1)**(j + 1) / (j * mp.factorial(j)))])
assert distincts([sk(lambda j: (-2)**(j - 1) * j**-j), sk(lambda j: 2**(j - 1) * j**-j), sk(lambda j: (-1)**(j + 1) * j**-j), sk(lambda j: (-1)**(j + 1) * 2**j * j**-j)])

# ---------------------------------------------------------------------------
# QCM
# ---------------------------------------------------------------------------
assert integrate(t**Rational(-2, 3), (t, 0, 1)) == 3                       # q1
assert limit(sin(t) / t**Rational(3, 2) * sqrt(t), t, 0) == 1              # q2
assert limit(t**2 / (t**2 + sin(t)), t, oo) == 1                           # q5
assert integrate(t**Rational(-1, 3), (t, 1, x)) == Rational(3, 2) * (x**Rational(2, 3) - 1)   # q6
assert proche(mp.nsum(lambda j: (-1)**j / (j + 1), [0, mp.inf]), mp.log(2))  # q8

# ---------------------------------------------------------------------------
# Compléments du relecteur indépendant
# ---------------------------------------------------------------------------
# cauchy-ex3 (v0, v1) : n = 1 suffit à réfuter les leurres 2/pi et 8/(3pi) ; la borne 1/pi (resp. 4/(3pi)) tient
# pour n = 1..6, et la limite (2/pi) ln 2 (resp. (2/pi) ln 3) la dépasse.
assert I_arches(1, 2, 1) < 2 / mp.pi and I_arches(1, 3, 1) < 8 / (3 * mp.pi)
assert all(I_arches(1, 2, N) >= 1 / mp.pi and I_arches(1, 3, N) >= 4 / (3 * mp.pi) for N in range(1, 7))
assert 2 / mp.pi * mp.log(2) > 1 / mp.pi and 2 / mp.pi * mp.log(3) > 4 / (3 * mp.pi)
# cauchy-ex3 v2 : le minorant sqrt(2n/pi) dépasse 1 dès n = 2.
assert mp.sqrt(4 / mp.pi) > 1 > mp.sqrt(2 / mp.pi)
# positive-ex4 v1, v2 : intervalles disjoints (4^-n + 4^-(n+1) < 1, idem pour 3).
assert all(F(1, 4**k) + F(1, 4**(k + 1)) < 1 and F(1, 3**k) + F(1, 3**(k + 1)) < 1 for k in range(1, 50))
# dominee-ex2 : forme exacte n/(n+1) pour n symbolique.
assert simplify(n * integrate(u**n, (u, 0, 1)) - n / (n + 1)) == 0
# QCM q7 : avec 1/t <= n < 1/t + 1, n/(1 + n^2 t^2) >= 1/(5t) sur ]0,1].
for tt in [mp.mpf(j) / 97 for j in range(1, 98)]:
    m = int(mp.ceil(1 / tt))
    assert m / (1 + (m * tt)**2) >= 1 / (5 * tt)
# QCM q3 : pics de hauteur n sur [n - 4^-n, n + 4^-n], disjoints, aire n 4^-n, total 4/9 (déjà vérifié plus haut).
# terme-a-terme-ex4 (v0, v1) : valeurs approchées citées dans les corrections.
assert abs(mp.nsum(lambda j: j**-j, [1, mp.inf]) - mp.mpf("1.2913")) < 1e-4
assert abs(mp.nsum(lambda j: (-1)**(j + 1) * j**-j, [1, mp.inf]) - mp.mpf("0.7834")) < 1e-4
# relations-ex4 v1 : h(t) = (ln t)^2 / t a pour dérivée ln t (2 - ln t) / t^2.
assert simplify(diff(log(t)**2 / t, t) - log(t) * (2 - log(t)) / t**2) == 0
# rectangles-ex2 v0, v2 : encadrements des sommes partielles (n = 2..300).
f1 = lambda v: 1 / (v * mp.log(v))
f2 = lambda v: 1 / (v * mp.sqrt(mp.log(v)))
s1 = s2 = mp.mpf(0)
for N in range(2, 301):
    s1 += f1(N); s2 += f2(N)
    l2 = mp.log(2)
    eps = mp.mpf("1e-25")    # égalité à droite pour N = 2
    assert mp.log(mp.log(N + 1)) - mp.log(l2) <= s1 <= 1 / (2 * l2) + mp.log(mp.log(N)) - mp.log(l2) + eps
    assert 2 * mp.sqrt(mp.log(N + 1)) - 2 * mp.sqrt(l2) <= s2 <= 1 / (2 * mp.sqrt(l2)) + 2 * mp.sqrt(mp.log(N)) - 2 * mp.sqrt(l2) + eps

print("int_l2_generalisees : toutes les vérifications passent.")
