"""Vérifications calculatoires du chapitre « int-t-integrale » (thème intégration, Terminale).

Cartes content/cartes/integration/terminale-int-t-integrale.json et exercices
content/exercices/integration/terminale-int-t-integrale.json. Chaque assert cite l'identifiant de la carte ou de
l'exercice (variante comprise : « v1 », « v2 » = variantes 1 et 2). Pour chaque réponse proposée, la bonne est vérifiée
égale à la valeur calculée et chaque leurre est vérifié faux.
Lancer : python3 verification/int_t_integrale.py (doit finir sans erreur).
"""

import math

from sympy import (E, Abs, Rational, cbrt, cos, diff, exp, integrate, log, oo, pi, simplify, sin, solve, sqrt,
                   symbols, limit, Piecewise, nsimplify)

x, t, y = symbols("x t y", real=True)
h = symbols("h", real=True)


def egal(a, b):
    return simplify(a - b) == 0


def different(a, b):
    """a et b (expressions numériques) diffèrent : contrôle exact puis numérique."""
    return simplify(a - b) != 0 and abs(float(a) - float(b)) > 1e-9


def fonctions_differentes(f, g, var=x, points=(-1.3, -0.4, 0.3, 0.7, 1.9)):
    return any(abs(float((f - g).subs(var, p))) > 1e-9 for p in points)


def bonne_et_leurres(valeur, bonne, leurres):
    assert egal(valeur, bonne), (valeur, bonne)
    for l in leurres:
        assert different(valeur, l), (valeur, l)


# ===========================================================================
# Carte int-t-integrale-fondamental
# ===========================================================================
# Contre-exemple / affirmation 4 : f = indicatrice de [0, +oo[ sur [-1,1], F(x) = int_{-1}^x f = max(x, 0).
F_marche = Piecewise((0, x <= 0), (x, True))
for p in (-0.7, -0.2, 0.3, 0.9):
    num = integrate(Piecewise((0, t < 0), (1, True)), (t, -1, p))
    assert abs(float(num) - float(F_marche.subs(x, p))) < 1e-12
# taux d'accroissement en 0 : 0 à gauche, 1 à droite -> non dérivable.
hpos, hneg = symbols("hpos", positive=True), symbols("hneg", negative=True)
assert simplify(F_marche.subs(x, hneg) / hneg) == 0 and simplify(F_marche.subs(x, hpos) / hpos) == 1
# Contre-exemple / affirmation 3 : G = -1/x et H = -1/x + 1 (x > 0) sont deux primitives de 1/x^2 sur R*.
assert egal(diff(-1 / x, x), 1 / x**2) and egal(diff(-1 / x + 1, x), 1 / x**2)
# Contre-exemple / affirmation 1 : sin + 2 primitive de cos, à valeurs dans [1,3] (jamais nulle).
assert egal(diff(sin(x) + 2, x), cos(x))
# Affirmation 2 : d/dx int_0^{x^2} f = 2x f(x^2), contrôlé sur f = cos.
assert egal(diff(integrate(cos(t), (t, 0, x**2)), x), 2 * x * cos(x**2))
# Exemple : int_1^x dt/t = ln x.
xp = symbols("xp", positive=True)
assert egal(integrate(1 / t, (t, 1, xp)), log(xp))

# int-t-integrale-fondamental-ex1 : F(x)=x^3-2x^2+5x-2, primitive de 3x^2-4x+5 avec F(1)=2.
f1 = 3 * x**2 - 4 * x + 5
def est_solution(F, f, x0, v0):
    return egal(diff(F, x), f) and egal(F.subs(x, x0), v0)
assert est_solution(x**3 - 2 * x**2 + 5 * x - 2, f1, 1, 2)
for l in (x**3 - 2 * x**2 + 5 * x + 2, 6 * x - 4, x**3 - 4 * x**2 + 5 * x):
    assert not est_solution(l, f1, 1, 2)
# ... v1 : 6x^2+2x-1, F(0)=3.
f1b = 6 * x**2 + 2 * x - 1
assert est_solution(2 * x**3 + x**2 - x + 3, f1b, 0, 3)
for l in (2 * x**3 + x**2 - x, 12 * x + 3, 6 * x**3 + 2 * x**2 - x + 3):
    assert not est_solution(l, f1b, 0, 3)
# ... v2 : 4x^3-6x+2, F(-1)=1.
f1c = 4 * x**3 - 6 * x + 2
assert est_solution(x**4 - 3 * x**2 + 2 * x + 5, f1c, -1, 1)
for l in (x**4 - 3 * x**2 + 2 * x + 1, x**4 - 3 * x**2 + 2 * x - 3, x**4 - 6 * x**2 + 2 * x + 8):
    assert not est_solution(l, f1c, -1, 1)


# int-t-integrale-fondamental-ex2 (et v1, v2) : F' = e^{-x^2} P(x) ; plus grand intervalle où P <= 0.
def intervalle_negatif(P):
    racines = sorted(solve(P, x))
    assert len(racines) == 2
    a, b = racines
    milieu = (a + b) / 2
    assert P.subs(x, milieu) < 0 and P.subs(x, a - 1) > 0 and P.subs(x, b + 1) > 0
    return a, b


def decroissante_sur(P, a, b):
    """F décroissante sur [a,b] <=> P <= 0 sur [a,b] (échantillonnage + racines exactes)."""
    pas = [a + (b - a) * Rational(k, 50) for k in range(51)]
    return all(P.subs(x, p) <= 0 for p in pas)


for P, attendu, leurres in (
        (x**2 - 4, (-2, 2), [(-4, 4), (0, 2), (2, 10)]),            # ex2
        (x**2 - 9, (-3, 3), [(-9, 9), (0, 3), (3, 10)]),            # ex2 v1
        (x**2 - 6 * x + 5, (1, 5), [(-5, -1), (0, 5), (5, 12)])):   # ex2 v2
    assert intervalle_negatif(P) == attendu
    for a, b in leurres:
        # Leurre faux : soit F n'y est pas décroissante, soit l'intervalle n'est pas le plus grand ([0,2] ⊊ [-2,2]).
        assert (not decroissante_sur(P, a, b)) or (attendu[0] <= a and b <= attendu[1] and (a, b) != attendu)
# Le dernier leurre de chaque jeu ([2,+oo[, [3,+oo[, [5,+oo[) : F y est strictement croissante.
assert (x**2 - 4).subs(x, 3) > 0 and (x**2 - 9).subs(x, 4) > 0 and (x**2 - 6 * x + 5).subs(x, 6) > 0

# int-t-integrale-fondamental-ex3 (et v1, v2) : dérivée de int_{u(x)}^{v(x)} e^{-t^2} dt.
g = exp(-t**2)
def deriv_borne(bas, haut):
    return diff(haut, x) * g.subs(t, haut) - diff(bas, x) * g.subs(t, bas)
for bas, haut, bonne, leurres in (
        (x, 2 * x, 2 * exp(-4 * x**2) - exp(-x**2),
         [exp(-4 * x**2) - exp(-x**2), 2 * exp(-4 * x**2) + exp(-x**2), 2 * exp(-2 * x**2) - exp(-x**2)]),
        (x, 3 * x, 3 * exp(-9 * x**2) - exp(-x**2),
         [exp(-9 * x**2) - exp(-x**2), 3 * exp(-9 * x**2) + exp(-x**2), 3 * exp(-3 * x**2) - exp(-x**2)]),
        (-x, x**2, 2 * x * exp(-x**4) + exp(-x**2),
         [2 * x * exp(-x**4) - exp(-x**2), exp(-x**4) + exp(-x**2), 2 * x * exp(-x**2) + exp(-x**2)])):
    d = deriv_borne(bas, haut)
    assert egal(d, bonne)
    for l in leurres:
        assert fonctions_differentes(d, l)
# Contrôle numérique de ex3 : taux d'accroissement de G en x = 0,7.
import mpmath
Gnum = lambda a: mpmath.quad(lambda s: mpmath.e**(-s**2), [a, 2 * a])
assert abs((Gnum(0.7 + 1e-6) - Gnum(0.7 - 1e-6)) / 2e-6 - (2 * math.exp(-4 * 0.49) - math.exp(-0.49))) < 1e-6


# int-t-integrale-fondamental-ex4 (et v1, v2) : phi(x) = int_0^x (x-t) f(t) dt.
def phi(f):
    return simplify(integrate((x - t) * f, (t, 0, x)))
p = phi(exp(t))
assert egal(p, exp(x) - x - 1) and egal(diff(p, x, 2), exp(x))
assert egal(p.subs(x, 0), 0) and egal(diff(p, x).subs(x, 0), 0)
bonne_et_leurres(p.subs(x, 1), E - 2, [E - 1, 1, E])
assert egal(integrate((1 - t) * exp(t), (t, 0, 1)), E - 2)
p = phi(cos(t))                                   # v1
assert egal(p, 1 - cos(x))
bonne_et_leurres(p.subs(x, pi), 2, [0, 1, pi])
p = phi(t**2)                                     # v2
assert egal(p, x**4 / 12)
bonne_et_leurres(p.subs(x, 2), Rational(4, 3), [Rational(8, 3), Rational(16, 3), Rational(2, 3)])


# int-t-integrale-fondamental-ex5 (et v1, v2) : f(x) = alpha(x) + beta int_0^x f.
def verifie(f, alpha, beta):
    return egal(f - (alpha + beta * integrate(f.subs(x, t), (t, 0, x))), 0)
assert verifie(exp(x), 1, 1)
C_ = symbols("C")
assert not verifie(2 * exp(x), 1, 1)              # leurre « C e^x » (C = 2 ne convient pas)
assert not verifie(1 + x, 1, 1) and not verifie(exp(x) - 1, 1, 1)
assert verifie(2 * exp(3 * x), 2, 3)              # v1
for l in (exp(3 * x), 3 * exp(2 * x), 2 * exp(x / 3)):
    assert not verifie(l, 2, 3)
assert verifie(exp(x) - 1, x, 1)                  # v2
for l in (exp(x), x * exp(x), 2 * exp(x) - 1):    # 2e^x - 1 : un cas du leurre « C e^x - 1 »
    assert not verifie(l, x, 1)


# ===========================================================================
# Carte int-t-integrale-proprietes
# ===========================================================================
# Exemple : sin x <= x et 1 - cos x <= x^2/2 pour x >= 0 (contrôle numérique).
for p in [k / 10 for k in range(0, 60)]:
    assert math.sin(p) <= p + 1e-15 and 1 - math.cos(p) <= p * p / 2 + 1e-15
assert integrate(cos(t), (t, 0, x)) == sin(x) and egal(integrate(sin(t), (t, 0, x)), 1 - cos(x))
# Contre-exemples : int_1^0 1 = -1 ; int_0^1 sin(pi x) = 2/pi ; int x^2 = 1/3 vs 1/4.
assert integrate(1, (x, 1, 0)) == -1
assert egal(integrate(sin(pi * x), (x, 0, 1)), 2 / pi) and sin(0) == 0 and sin(pi) == 0
assert integrate(x**2, (x, 0, 1)) == Rational(1, 3) and integrate(x, (x, 0, 1))**2 == Rational(1, 4)
# Affirmation 1 : int_{-1}^2 x = 3/2 >= 0.
assert integrate(x, (x, -1, 2)) == Rational(3, 2)
# Affirmation 3 : f = -2, g = 1 sur [0,1].
assert abs(integrate(-2, (x, 0, 1))) == 2 and integrate(1, (x, 0, 1)) == 1
# Affirmation 2 (Chasles hors de l'intervalle) sur un exemple : f = e^x.
assert egal(integrate(exp(x), (x, -1, 1)), integrate(exp(x), (x, -1, 3)) + integrate(exp(x), (x, 3, 1)))

# int-t-integrale-proprietes-ex1 : 2*4 - 3*(-2) = 14.
If25 = 7 - 3
bonne_et_leurres(2 * If25 - 3 * (-2), 14, [2, 26, -2])
assert 2 * If25 + 3 * (-2) == 2 and 2 * (3 + 7) + 6 == 26 and 2 * (3 - 7) + 6 == -2   # origine des leurres
# ... v1 : int_4^6 f = 2 - 5 = -3 ; 3(-3) + 2*3 = -3.
bonne_et_leurres(3 * (2 - 5) + 2 * 3, -3, [27, 15, -15])
# ... v2 : int_{-1}^0 f = 4 - 1 = 3 ; 3 - 4*2 = -5.
bonne_et_leurres((4 - 1) - 4 * 2, -5, [-11, 11, 3])

# int-t-integrale-proprietes-ex2 (et v1, v2).
bonne_et_leurres(integrate(Abs(x), (x, -1, 2)), Rational(5, 2), [Rational(3, 2), 2, 3])
assert integrate(-x, (x, -1, 0)) == Rational(1, 2) and integrate(x, (x, 0, 2)) == 2
v = integrate(1 - x, (x, -2, 1)) + integrate(x - 1, (x, 1, 3))
assert integrate(1 - x, (x, -2, 1)) == Rational(9, 2) and integrate(x - 1, (x, 1, 3)) == 2
bonne_et_leurres(v, Rational(13, 2), [Rational(-5, 2), Rational(5, 2), Rational(9, 2)])
assert integrate(x - 1, (x, -2, 3)) == Rational(-5, 2)          # leurre : sans valeur absolue
v = integrate(1 - x**2, (x, 0, 1)) + integrate(x**2 - 1, (x, 1, 3))
assert integrate(x**2 - 1, (x, 1, 3)) == Rational(20, 3)
bonne_et_leurres(v, Rational(22, 3), [6, Rational(20, 3), 8])
assert integrate(x**2 - 1, (x, 0, 3)) == 6                         # leurre : sans valeur absolue

# int-t-integrale-proprietes-ex3 (et v1, v2) : u_n = int_0^1 x^n/(c+x) ; c u_n + u_{n+1} = 1/(n+1).
for c, bonne, leurres in (
        (1, log(2) - Rational(1, 2), [Rational(1, 2) - log(2), 1 - log(2), Rational(3, 2) - log(2)]),
        (2, 4 * log(Rational(3, 2)) - Rational(3, 2),
         [Rational(3, 2) - 4 * log(Rational(3, 2)), 1 - 2 * log(Rational(3, 2)), 2 * log(Rational(3, 2)) - Rational(1, 2)]),
        (3, 9 * log(Rational(4, 3)) - Rational(5, 2),
         [Rational(5, 2) - 9 * log(Rational(4, 3)), 1 - 3 * log(Rational(4, 3)), 3 * log(Rational(4, 3)) - Rational(1, 2)])):
    u = [integrate(x**n / (c + x), (x, 0, 1)) for n in range(6)]
    assert egal(u[0], log(Rational(c + 1, c)))
    for n in range(5):
        assert egal(c * u[n] + u[n + 1], Rational(1, n + 1))
        assert 0 <= float(u[n]) <= 1 / (c * (n + 1)) + 1e-15            # encadrement de la question 2
    assert egal(u[1], 1 - c * log(Rational(c + 1, c)))
    bonne_et_leurres(u[2], bonne, leurres)
    assert float(u[2]) > 0
# Limite : 1/(n+1) -> 0.
assert limit(1 / (h + 1), h, oo) == 0

# int-t-integrale-proprietes-ex4 (et v1, v2) : x - x^2/2 <= ln(1+x) <= x ; encadrements proposés.
assert egal(integrate(1 - t, (t, 0, x)), x - x**2 / 2) and egal(integrate(1 / (1 + t), (t, 0, xp)), log(1 + xp))
assert egal((1 - t) * (1 + t), 1 - t**2)
for p in [k / 20 for k in range(0, 100)]:
    assert p - p * p / 2 <= math.log(1 + p) <= p + 1e-15
for xv, bonne, leurres in (
        (Rational(2, 10), (Rational(18, 100), Rational(2, 10)),
         [(Rational(2, 10), Rational(22, 100)), (Rational(16, 100), Rational(18, 100)), (Rational(19, 100), Rational(21, 100))]),
        (Rational(1, 10), (Rational(95, 1000), Rational(1, 10)),
         [(Rational(1, 10), Rational(105, 1000)), (Rational(9, 100), Rational(95, 1000)), (Rational(975, 10000), Rational(1025, 10000))]),
        (Rational(3, 10), (Rational(255, 1000), Rational(3, 10)),
         [(Rational(3, 10), Rational(345, 1000)), (Rational(21, 100), Rational(255, 1000)), (Rational(285, 1000), Rational(3, 10))])):
    assert (xv - xv**2 / 2, xv) == bonne
    lv = math.log(1 + float(xv))
    assert float(bonne[0]) <= lv <= float(bonne[1])
    for a, b in leurres:
        assert not (float(a) <= lv <= float(b))        # chaque leurre est un encadrement faux


# ===========================================================================
# Carte int-t-integrale-ipp
# ===========================================================================
# Exemple : int_0^pi x cos x = -2.
assert integrate(x * cos(x), (x, 0, pi)) == -2
# Contre-exemple (signe) : u = v = x sur [0,1] : 1/2 contre 3/2.
assert integrate(x, (x, 0, 1)) == Rational(1, 2) and 1 + integrate(x, (x, 0, 1)) == Rational(3, 2)
# Affirmation 2 : int_0^1 x^2 e^x = e - 2.
assert egal(integrate(x**2 * exp(x), (x, 0, 1)), E - 2)

# int-t-integrale-ipp-ex1 (et v1, v2).
bonne_et_leurres(integrate(x * exp(x), (x, 0, 1)), 1, [E - 1, 2 * E - 1, E])
bonne_et_leurres(integrate(x * exp(2 * x), (x, 0, 1)), (E**2 + 1) / 4, [(E**2 - 1) / 4, E**2 / 2, (3 * E**2 - 1) / 4])
bonne_et_leurres(integrate(x * exp(-x), (x, 0, 1)), 1 - 2 / E, [1 - 1 / E, -1, 1 / E])

# int-t-integrale-ipp-ex2 (et v1, v2).
bonne_et_leurres(integrate(log(x), (x, 1, E)), 1, [E - 1, E, 0])
bonne_et_leurres(integrate(x * log(x), (x, 1, E)), (E**2 + 1) / 4, [(E**2 - 1) / 4, E**2 / 2, (E**2 + 1) / 2])
bonne_et_leurres(integrate(x**2 * log(x), (x, 1, E)), (2 * E**3 + 1) / 9, [(2 * E**3 - 1) / 9, E**3 / 3, (E**3 + 1) / 9])

# int-t-integrale-ipp-ex3 (et v1, v2).
bonne_et_leurres(integrate(x**2 * exp(x), (x, 0, 1)), E - 2, [E - 1, 2 - E, E + 2])
assert integrate(x * sin(x), (x, 0, pi / 2)) == 1
bonne_et_leurres(integrate(x**2 * cos(x), (x, 0, pi / 2)), pi**2 / 4 - 2, [pi**2 / 4, pi**2 / 4 + 2, 2 - pi**2 / 4])
bonne_et_leurres(integrate(x**2 * exp(-x), (x, 0, 1)), 2 - 5 / E, [2 - 3 / E, 5 / E - 2, 1 - 2 / E])

# int-t-integrale-ipp-ex4 (et v1, v2) : relations de récurrence et troisième terme.
I = [integrate(x**n * exp(x), (x, 0, 1)) for n in range(6)]
assert egal(I[0], E - 1)
for n in range(5):
    assert egal(I[n + 1], E - (n + 1) * I[n])
    assert 0 <= float(I[n]) <= float(E) / (n + 1)
assert egal(I[1], 1) and egal(I[2], E - 2)
bonne_et_leurres(I[3], 6 - 2 * E, [2 * E - 6, E - 3, 6 - 3 * E])
J = [integrate(x**n * exp(-x), (x, 0, 1)) for n in range(6)]
assert egal(J[0], 1 - 1 / E)
for n in range(5):
    assert egal(J[n + 1], (n + 1) * J[n] - 1 / E)
    assert 0 <= float(J[n]) <= 1 / (n + 1)
assert egal(J[1], 1 - 2 / E) and egal(J[2], 2 - 5 / E)
bonne_et_leurres(J[3], 6 - 16 / E, [16 / E - 6, 6 - 15 / E, 2 - 5 / E])
K = [integrate(x**n * exp(2 * x), (x, 0, 1)) for n in range(6)]
assert egal(K[0], (E**2 - 1) / 2)
for n in range(5):
    assert egal(K[n + 1], E**2 / 2 - Rational(n + 1, 2) * K[n])
    assert 0 <= float(K[n]) <= float(E**2) / (n + 1)
assert egal(K[1], (E**2 + 1) / 4) and egal(K[2], (E**2 - 1) / 4)
bonne_et_leurres(K[3], (E**2 + 3) / 8, [(E**2 - 3) / 8, (3 - E**2) / 4, (E**2 - 1) / 4])

# int-t-integrale-ipp-ex5 (et v1, v2) : intégrations par parties cycliques.
Iv = integrate(exp(x) * sin(x), (x, 0, pi))
Jv = integrate(exp(x) * cos(x), (x, 0, pi))
assert egal(Iv, -Jv) and egal(Jv, -exp(pi) - 1 + Iv)
bonne_et_leurres(Iv, (exp(pi) + 1) / 2, [(exp(pi) - 1) / 2, exp(pi) + 1, -(exp(pi) + 1) / 2])
Iv = integrate(exp(x) * cos(x), (x, 0, pi / 2))
Jv = integrate(exp(x) * sin(x), (x, 0, pi / 2))
assert egal(Iv, -1 + Jv) and egal(Jv, exp(pi / 2) - Iv)
bonne_et_leurres(Iv, (exp(pi / 2) - 1) / 2, [(exp(pi / 2) + 1) / 2, exp(pi / 2) - 1, (1 - exp(pi / 2)) / 2])
Iv = integrate(exp(2 * x) * sin(x), (x, 0, pi))
Jv = integrate(exp(2 * x) * cos(x), (x, 0, pi))
assert egal(Iv, exp(2 * pi) + 1 + 2 * Jv) and egal(Jv, -2 * Iv)
bonne_et_leurres(Iv, (exp(2 * pi) + 1) / 5, [(exp(2 * pi) + 1) / 2, (exp(2 * pi) - 1) / 5, (exp(2 * pi) + 1) / 3])


# ===========================================================================
# Carte int-t-integrale-moyenne
# ===========================================================================
# Exemple : moyenne de sin sur [0, pi] = 2/pi, dans [0,1], atteinte (2/pi < 1).
assert egal(integrate(sin(x), (x, 0, pi)) / pi, 2 / pi) and 0 < 2 / math.pi < 1
# Contre-exemples : indicatrice de [1/2,1] sur [0,1] : moyenne 1/2 (jamais prise) ; 4x(1-x) : 2/3 ; x^2 : 1/3.
assert integrate(Piecewise((0, x < Rational(1, 2)), (1, True)), (x, 0, 1)) == Rational(1, 2)
assert integrate(4 * x * (1 - x), (x, 0, 1)) == Rational(2, 3)
assert integrate(x**2, (x, 0, 1)) == Rational(1, 3) and (0 + 1) / 2 == 0.5
# Affirmations : f = x : moyenne de f^2 = 1/3, carré de la moyenne = 1/4 ; f = 1 sur [0,2] : 2.
assert integrate(1, (x, 0, 2)) == 2

# Cours de int-t-integrale-moyenne : mu = (f(a)+f(b))/2 pour f affine, faux pour x^2 sur [0,1] (1/3 contre 1/2) ;
# l'égalité peut aussi avoir lieu sans f affine (x^3 sur [-1,1]), d'où la formulation « en général ».
a_, b_, p_, q_ = symbols("a_ b_ p_ q_", real=True)
assert egal(integrate(p_ * x + q_, (x, a_, b_)) / (b_ - a_), ((p_ * a_ + q_) + (p_ * b_ + q_)) / 2)
assert integrate(x**2, (x, 0, 1)) != Rational(1, 2) and integrate(x**3, (x, -1, 1)) / 2 == 0

# int-t-integrale-moyenne-ex1 (et v1, v2).
moy = lambda f, a, b: integrate(f, (x, a, b)) / (b - a)
bonne_et_leurres(moy(x**2, 0, 3), 3, [9, Rational(9, 2), Rational(9, 4)])
assert egal(sqrt(3)**2, 3)                                      # atteinte en c = sqrt(3)
bonne_et_leurres(moy(x**3, 0, 2), 2, [4, 1, 8])
bonne_et_leurres(moy(3 * x**2 - 2 * x, 1, 3), 9, [18, 11, 8])
assert (1 + 21) / 2 == 11                                         # origine du leurre (f(1)+f(3))/2

# int-t-integrale-moyenne-ex2 (et v1, v2).
bonne_et_leurres(integrate(4 * exp(-t / 2), (t, 0, 4)) / 4, 2 * (1 - exp(-2)),
                 [8 * (1 - exp(-2)), 2 * (1 + exp(-2)), 4 * exp(-1)])
assert abs(float(2 * (1 - exp(-2))) - 1.73) < 0.005
assert egal(integrate(sin(pi * t / 12), (t, 0, 12)), 24 / pi)
bonne_et_leurres(integrate(20 + 5 * sin(pi * t / 12), (t, 0, 12)) / 12, 20 + 10 / pi, [20, 25, 20 + 5 / pi])
assert abs(float(20 + 10 / pi) - 23.2) < 0.05
bonne_et_leurres(integrate(10 * (1 - exp(-t)), (t, 0, 2)) / 2, 5 * (1 + exp(-2)),
                 [10 * (1 + exp(-2)), 5 * (1 - exp(-2)), 10 * (1 - exp(-1))])
assert abs(float(5 * (1 + exp(-2))) - 5.68) < 0.005


# int-t-integrale-moyenne-ex3 (et v1, v2) : encadrements de sommes ; chaque leurre est faux pour un n.
def verifie_encadrement(terme, bas, haut, leurres, nmax=300):
    s = 0.0
    contre = [False] * len(leurres)
    for n in range(1, nmax + 1):
        s += terme(n)
        assert bas(n) - 1e-12 <= s <= haut(n) + 1e-12, n
        for i, (lb, lh) in enumerate(leurres):
            if not (lb(n) - 1e-12 <= s <= lh(n) + 1e-12):
                contre[i] = True
    assert all(contre), contre


L = math.log
verifie_encadrement(lambda k: 1 / k, lambda n: L(n + 1), lambda n: 1 + L(n), [
    (lambda n: L(n), lambda n: L(n + 1)),
    (lambda n: 1 + L(n + 1), lambda n: 2 + L(n)),
    (lambda n: L(n + 1) - 1, lambda n: L(n))])
verifie_encadrement(lambda k: 1 / math.sqrt(k), lambda n: 2 * math.sqrt(n + 1) - 2, lambda n: 2 * math.sqrt(n) - 1, [
    (lambda n: 2 * math.sqrt(n) - 2, lambda n: 2 * math.sqrt(n + 1) - 2),
    (lambda n: 2 * math.sqrt(n + 1) - 1, lambda n: 2 * math.sqrt(n)),
    (lambda n: math.sqrt(n + 1) - 1, lambda n: math.sqrt(n))])
verifie_encadrement(lambda k: 1 / k**2, lambda n: 1 - 1 / (n + 1), lambda n: 2 - 1 / n, [
    (lambda n: 1 - 1 / n, lambda n: 1 - 1 / (n + 1)),
    (lambda n: 1 / (n + 1), lambda n: 1 / n),
    (lambda n: 2 - 1 / n, lambda n: 3 - 1 / n)])
# Intégrales utilisées : int_k^{k+1} dx/sqrt(x) = 2sqrt(k+1) - 2sqrt(k) ; int_1^n dx/x^2 = 1 - 1/n.
kk = symbols("k", positive=True)
assert egal(integrate(1 / sqrt(x), (x, kk, kk + 1)), 2 * sqrt(kk + 1) - 2 * sqrt(kk))
assert egal(integrate(1 / x**2, (x, 1, kk)), 1 - 1 / kk)


# int-t-integrale-moyenne-ex4 (et v1, v2) : limites de la valeur moyenne en 0+ et en +oo.
def moyenne_num(f, X):
    return float(mpmath.quad(f, [0, X])) / X


for f, l0, linf, borne in (
        (lambda s: mpmath.e**(-s**2), 1, 0, lambda X: (1 + math.exp(-1)) / X),
        (lambda s: 2 / (1 + s**2), 2, 0, lambda X: 4 / X),
        (lambda s: 3 + mpmath.e**(-s**2), 4, 3, lambda X: 3 + (1 + math.exp(-1)) / X)):
    assert abs(moyenne_num(f, 1e-4) - l0) < 1e-6
    assert abs(moyenne_num(f, 1e5) - linf) < 1e-4
    for X in (1, 2, 5, 20, 100):
        assert linf <= moyenne_num(f, X) <= borne(X) + 1e-12      # majoration de la question 2
    for X in (0.1, 0.5, 1, 3):
        assert float(f(X)) <= moyenne_num(f, X) <= l0 + 1e-12     # encadrement de la question 1
# Leurres (couples de limites) distincts du bon couple pour chaque jeu.
for bonne, leurres in (((1, 0), [(0, 0), (1, 1), (0, 1)]), ((2, 0), [(1, 0), (2, 2), (2, 1)]), ((4, 3), [(4, 0), (3, 3), (4, 4)])):
    assert all(l != bonne for l in leurres) and len(set(leurres)) == 3
# Majorations utilisées : 1 + int_1^x e^{-t} <= 1 + 1/e ; 2 + int_1^x 2/t^2 = 4 - 2/x.
assert egal(1 + integrate(exp(-t), (t, 1, xp)), 1 + exp(-1) - exp(-xp))
assert egal(2 + integrate(2 / t**2, (t, 1, xp)), 4 - 2 / xp)


# ===========================================================================
# Carte int-t-integrale-aire
# ===========================================================================
# Exemple : int_0^1 (x - x^3) = 1/4.
assert integrate(x - x**3, (x, 0, 1)) == Rational(1, 4)
# Contre-exemples / affirmations : int_{-1}^1 x = 0, aire 1 ; x et x^2 sur [0,2] : -2/3, aire 1 ; 4x(1-x) : -2/3.
assert integrate(x, (x, -1, 1)) == 0 and integrate(Abs(x), (x, -1, 1)) == 1
assert integrate(x - x**2, (x, 0, 2)) == Rational(-2, 3)
assert integrate(x - x**2, (x, 0, 1)) == Rational(1, 6) and integrate(x**2 - x, (x, 1, 2)) == Rational(5, 6)
assert abs(mpmath.quad(lambda s_: abs(s_ - s_**2), [0, 1, 2]) - 1) < 1e-12
assert integrate(-4 * x * (1 - x), (x, 0, 1)) == Rational(-2, 3)
# Leurre du trou 2 (|int g| - |int f|) : f = -1, g = 0 sur [0,1] : aire 1, formule -1.
assert abs(integrate(0, (x, 0, 1))) - abs(integrate(-1, (x, 0, 1))) == -1

# int-t-integrale-aire-ex1 (et v1, v2).
bonne_et_leurres(integrate(x - x**2, (x, 0, 1)), Rational(1, 6), [Rational(1, 2), Rational(1, 3), Rational(-1, 6)])
bonne_et_leurres(integrate(2 * x - x**2, (x, 0, 2)), Rational(4, 3), [4, Rational(8, 3), Rational(-4, 3)])
bonne_et_leurres(integrate(sqrt(x) - x**2, (x, 0, 1)), Rational(1, 3), [Rational(2, 3), 1, Rational(-1, 3)])

# int-t-integrale-aire-ex2 (et v1, v2).
bonne_et_leurres(integrate(4 - x**2, (x, 0, 2)) + integrate(x**2 - 4, (x, 2, 3)), Rational(23, 3), [3, Rational(16, 3), 9])
assert abs(mpmath.quad(lambda s_: abs(s_**2 - 4), [0, 2, 3]) - 23 / 3) < 1e-12
assert integrate(4 - x**2, (x, 0, 2)) == Rational(16, 3) and integrate(x**2 - 4, (x, 2, 3)) == Rational(7, 3)
assert integrate(x**2 - 4, (x, 0, 3)) == -3
bonne_et_leurres(integrate(x - x**3, (x, 0, 1)) + integrate(x**3 - x, (x, 1, 2)), Rational(5, 2),
                 [2, Rational(9, 4), Rational(1, 4)])
assert integrate(x**3 - x, (x, 1, 2)) == Rational(9, 4) and integrate(x**3 - x, (x, 0, 2)) == 2
bonne_et_leurres(integrate(sin(x), (x, 0, pi)) + integrate(-sin(x), (x, pi, 3 * pi / 2)), 3, [1, 2, -1])
assert integrate(sin(x), (x, 0, 3 * pi / 2)) == 1

# int-t-integrale-aire-ex3 (et v1, v2) : domaine borné entre deux courbes.
assert sorted(solve(x**2 - x - 2, x)) == [-1, 2]
bonne_et_leurres(integrate(x + 2 - x**2, (x, -1, 2)), Rational(9, 2), [Rational(-9, 2), Rational(10, 3), Rational(7, 6)])
assert integrate(x + 2 - x**2, (x, 0, 2)) == Rational(10, 3) and integrate(x + 2 - x**2, (x, -1, 0)) == Rational(7, 6)
assert sorted(solve(x**2 - 2 * x - 3, x)) == [-1, 3]
bonne_et_leurres(integrate(2 * x + 3 - x**2, (x, -1, 3)), Rational(32, 3), [Rational(-32, 3), 9, Rational(5, 3)])
assert integrate(2 * x + 3 - x**2, (x, 0, 3)) == 9 and integrate(2 * x + 3 - x**2, (x, -1, 0)) == Rational(5, 3)
assert sorted(solve((1 - x**2) - (x**2 - 1), x)) == [-1, 1]
bonne_et_leurres(integrate(2 - 2 * x**2, (x, -1, 1)), Rational(8, 3), [Rational(4, 3), Rational(-8, 3), 4])

# int-t-integrale-aire-ex4 (et v1, v2) : aire sous |ln| ; primitive x ln x - x.
Lp = x * log(x) - x
assert egal(diff(Lp, x), log(x))
assert egal(integrate(log(t), (t, 1, xp)), xp * log(xp) - xp + 1)
for a, b, bonne, leurres in (
        (1 / E, E, 2 - 2 / E, [2 / E, 1, 2]),
        (exp(-2), E, 2 - 3 / E**2, [3 / E**2, 1, 2 - 2 / E**2]),
        (Rational(1, 2), 2, Rational(3, 2) * log(2) - Rational(1, 2),
         [Rational(5, 2) * log(2) - Rational(3, 2), 2 * log(2) - 1, Rational(3, 2) * log(2) + Rational(1, 2)])):
    aire = -(Lp.subs(x, 1) - Lp.subs(x, a)) + (Lp.subs(x, b) - Lp.subs(x, 1))
    num = mpmath.quad(lambda s: abs(mpmath.log(s)), [float(a), 1, float(b)])
    assert abs(float(aire) - float(num)) < 1e-10
    bonne_et_leurres(aire, bonne, leurres)
# Intégrales algébriques (leurres) : int_{1/e}^e ln = 2/e ; int_{e^-2}^e ln = 3/e^2 ; int_{1/2}^2 ln = 5/2 ln 2 - 3/2.
assert egal(integrate(log(x), (x, 1 / E, E)), 2 / E)
assert egal(integrate(log(x), (x, exp(-2), E)), 3 / E**2)
assert egal(integrate(log(x), (x, Rational(1, 2), 2)), Rational(5, 2) * log(2) - Rational(3, 2))

# int-t-integrale-aire-ex5 (et v1, v2) : la droite y = m x partage en deux l'aire sous y = a x - x^2.
m = symbols("m", positive=True)
for a, bonne, leurres in (
        (2, 2 - cbrt(4), [1, 2 - sqrt(2), 2 - cbrt(2)]),
        (3, 3 - 3 / cbrt(2), [Rational(3, 2), 3 - 3 / sqrt(2), 3 - 3 / cbrt(4)]),
        (4, 4 - 2 * cbrt(4), [2, 4 - 2 * sqrt(2), 4 - 2 * cbrt(2)])):
    total = integrate(a * x - x**2, (x, 0, a))
    assert egal(total, Rational(a**3, 6))
    A = integrate((a - m) * x - x**2, (x, 0, a - m))
    assert egal(A, (a - m)**3 / 6)
    assert egal(A.subs(m, bonne), total / 2) and 0 < float(bonne) < a
    for l in leurres:
        assert 0 < float(l) < a and abs(float(A.subs(m, l) - total / 2)) > 1e-6
    # Contrôle numérique direct : aire de la partie de D située au-dessus de la droite.
    mb = float(bonne)
    haut = mpmath.quad(lambda s: max(a * s - s * s - mb * s, 0), [0, a - mb, a])
    assert abs(float(haut) - float(total) / 2) < 1e-9
assert abs(float(2 - cbrt(4)) - 0.413) < 5e-4 and abs(float(3 - 3 / cbrt(2)) - 0.619) < 5e-4
assert abs(float(4 - 2 * cbrt(4)) - 0.825) < 5e-4 and egal(cbrt(32), 2 * cbrt(4))

print("int_t_integrale : toutes les vérifications passent.")
