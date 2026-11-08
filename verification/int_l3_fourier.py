"""Vérifications calculatoires du chapitre « int-l3-fourier » (intégration, L3, séries de Fourier) :
cartes, exercices et QCM de content/*/integration/l3-int-l3-fourier.json et content/qcm/int-l3-fourier.json.

Conventions du chapitre : c_n(f) = (1/2π)∫_{-π}^{π} f(t)e^{-int} dt, a_n = (1/π)∫ f cos(nt), b_n = (1/π)∫ f sin(nt),
S_N(f) = Σ_{|n|≤N} c_n e_n, σ_N(f) = (1/N)Σ_{k<N} S_k(f).
Chaque assert cite l'identifiant de la carte, de l'exercice (variantes comprises : v1, v2) ou de la question.
Les intégrales exactes sont calculées en SymPy ; les sommes de séries obtenues par Parseval ou Dirichlet sont contrôlées
numériquement avec mpmath (30 chiffres, nsum, polylog), les leurres sont vérifiés différents de la bonne réponse.
Lancer : python3 verification/int_l3_fourier.py (doit finir sans erreur).
"""

from fractions import Fraction as Fr
import math

import mpmath as mp
from sympy import (Abs, Rational, chebyshevt, cos, cosh, coth, exp, expand, integrate, pi, simplify, sin, sinh, sqrt,
                   symbols, I, Integer, trigsimp, nsimplify)

mp.mp.dps = 30
t, x = symbols("t x", real=True)
n = symbols("n", integer=True, positive=True)
EPS = mp.mpf(10) ** -20
PI = mp.pi


def proche(a, b, tol=EPS):
    return abs(mp.mpc(a) - mp.mpc(b)) < tol


def distincts(vraie, leurres, tol=mp.mpf(10) ** -6):
    """La bonne valeur diffère numériquement de chaque leurre, et les leurres entre eux."""
    vals = [mp.mpc(vraie)] + [mp.mpc(l) for l in leurres]
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            assert abs(vals[i] - vals[j]) > tol, (vals[i], vals[j])


def c_num(f, k, pts=None):
    """c_k(f) numérique, f donnée sur [-π, π] (points de coupure facultatifs)."""
    pts = pts or [-PI, PI]
    return mp.quad(lambda s: f(s) * mp.expj(-k * s), pts) / (2 * PI)


def a_num(f, k, pts=None):
    pts = pts or [-PI, PI]
    return mp.quad(lambda s: f(s) * mp.cos(k * s), pts) / PI


def b_num(f, k, pts=None):
    pts = pts or [-PI, PI]
    return mp.quad(lambda s: f(s) * mp.sin(k * s), pts) / PI


def a_sym(expr, k, lo=-pi, hi=pi):
    return simplify(integrate(expr * cos(k * t), (t, lo, hi)) / pi)


def b_sym(expr, k, lo=-pi, hi=pi):
    return simplify(integrate(expr * sin(k * t), (t, lo, hi)) / pi)


creneau = lambda s: mp.sign(s)          # sur ]-π, π[ : 1 sur ]0,π[, -1 sur ]-π,0[
dent = lambda s: s                       # f(t) = t sur ]-π, π]

# ===========================================================================
# int-l3-fourier-coefficients (carte)
# ===========================================================================
# Orthonormalité : (1/2π)∫ e^{int} e^{-imt} = δ_{n,m}.
for p in range(-3, 4):
    for q in range(-3, 4):
        assert proche(mp.quad(lambda s: mp.expj((p - q) * s), [-PI, PI]) / (2 * PI), 1 if p == q else 0)
# Exemple (créneau) : b_n = 2(1-(-1)^n)/(πn), a_n = 0, c_n = -2i/(πn) pour n impair (aussi n < 0), 0 pour n pair ≠ 0.
for k in range(1, 9):
    assert proche(b_num(creneau, k, [-PI, 0, PI]), 2 * (1 - (-1) ** k) / (PI * k))
    assert proche(a_num(creneau, k, [-PI, 0, PI]), 0)
for k in list(range(-7, 0)) + list(range(1, 8)):
    assert proche(c_num(creneau, k, [-PI, 0, PI]), (-2j / (PI * k)) if k % 2 else 0)
# Contre-exemples et pièges : sin réelle avec b_1 = 1 ; f = i paire avec c_0 = i ≠ conj ; c_n avec 1/π donnerait c_0(1) = 2.
assert b_sym(sin(t), 1) == 1 and integrate(I, (t, -pi, pi)) / (2 * pi) == I
assert integrate(Integer(1), (t, -pi, pi)) / pi == 2
# Trous : pour P = e_1, c_1 = 1 ≠ γ_{-1} = 0 ; f = sin : c_1 = -i/2 = (a_1 - i b_1)/2, (a_1 + i b_1)/2 = i/2 ;
# f = cos : c_1 = 1/2 ≠ a_1 - i b_1 = 1 ; leurre (a+b)/2 pour sin : 1/2 ≠ -i/2 ; sin(t) impaire, réelle, continue, b_1 = 1.
assert simplify(integrate(sin(t) * exp(-I * t), (t, -pi, pi)) / (2 * pi) + I / 2) == 0
assert simplify(integrate(cos(t) * exp(-I * t), (t, -pi, pi)) / (2 * pi) - Rational(1, 2)) == 0
# Affirmation 3 : sin t = (e^{it} - e^{-it})/(2i), γ_1 = -i/2.
assert simplify((exp(I * t) - exp(-I * t)) / (2 * I) - sin(t)).rewrite(exp).simplify() == 0

# int-l3-fourier-coefficients-ex1 (+ v1, v2) : linéarisations.
assert simplify(cos(t) ** 3 - (Rational(3, 4) * cos(t) + Rational(1, 4) * cos(3 * t))) == 0
assert a_sym(cos(t) ** 3, 3) == Rational(1, 4) and a_sym(cos(t) ** 3, 1) == Rational(3, 4)
assert simplify(sin(t) ** 3 - (Rational(3, 4) * sin(t) - Rational(1, 4) * sin(3 * t))) == 0
assert b_sym(sin(t) ** 3, 3) == Rational(-1, 4) and b_sym(sin(t) ** 3, 1) == Rational(3, 4)
assert simplify(integrate(sin(t) ** 3 * exp(-3 * I * t), (t, -pi, pi)) / (2 * pi) - I / 8) == 0
assert a_sym(cos(t) ** 4, 0) == Rational(3, 4) and a_sym(cos(t) ** 4, 2) == Rational(1, 2) and a_sym(cos(t) ** 4, 4) == Rational(1, 8)
for vraie, leurres in [(Fr(1, 4), [Fr(3, 4), Fr(1, 8), Fr(1, 2)]), (Fr(-1, 4), [Fr(1, 4), Fr(-3, 4), Fr(3, 4)]),
                       (Fr(3, 4), [Fr(3, 8), Fr(3, 2), Fr(1, 8)])]:
    assert vraie not in leurres and len(set(leurres)) == 3

# int-l3-fourier-coefficients-ex2 (+ v1, v2).
for k in range(1, 9):
    assert b_sym(t, k) == 2 * (-1) ** (k + 1) / Integer(k) and a_sym(t, k) == 0
    # leurres 2/n, 2(-1)^n/n, (-1)^{n+1}/n : faux pour n = 1 ou n = 2.
    assert b_sym(t, k, 0, 2 * pi) == -2 / Integer(k) and a_sym(t, k, 0, 2 * pi) == 0      # v1
    assert simplify(a_sym(Abs(t), k) - 2 * ((-1) ** k - 1) / (pi * k ** 2)) == 0             # v2
assert a_sym(t, 0, 0, 2 * pi) == 2 * pi and a_sym(Abs(t), 0) == pi
# Leurres de ex2 (formules en n) : chacun diffère de la bonne formule pour n = 1 ou n = 2 ; idem pour v1 et v2.
def differe(vraie, leurre):
    return any(simplify(vraie.subs(n, k) - leurre.subs(n, k)) != 0 for k in (1, 2, 3))
for vraie, leurres in [(2 * (-1) ** (n + 1) / n, [2 / n, 2 * (-1) ** n / n, (-1) ** (n + 1) / n]),
                       (-2 / n, [2 / n, 2 * (-1) ** (n + 1) / n, -1 / n]),
                       (2 * ((-1) ** n - 1) / (pi * n**2), [2 * (1 - (-1) ** n) / (pi * n**2), ((-1) ** n - 1) / (pi * n**2),
                                                          2 * ((-1) ** n - 1) / (pi * n)])]:
    assert all(differe(vraie, l) for l in leurres)
for k in range(1, 6):   # c_n = i(-1)^n / n (main), i/n (v1)
    assert proche(c_num(dent, k), mp.mpc(0, (-1) ** k) / k) and proche(c_num(dent, -k), mp.mpc(0, (-1) ** k) / (-k))
    assert proche(mp.quad(lambda s: s * mp.expj(-k * s), [0, 2 * PI]) / (2 * PI), mp.mpc(0, 1) / k)
# Leurres de v2 faux pour n = 1 : 2(1-(-1)^n)/(πn²) = 4/π, ((-1)^n-1)/(πn²) = -2/π, 2((-1)^n-1)/(πn) = -4/π = vraie ?
# (pour n = 1, ce dernier coïncide ; il diffère pour n = 3 : -4/(3π) contre -4/(9π)).
assert 2 * ((-1) ** 3 - 1) / (PI * 3) != 2 * ((-1) ** 3 - 1) / (PI * 9)

# int-l3-fourier-coefficients-ex3 (+ v1, v2) : c_n(e^{at}) = (-1)^n sh(aπ)/(π(a - in)).
for a in [1, 2, -1]:
    for k in range(-4, 5):
        exact = (-1) ** k * mp.sinh(a * PI) / (PI * (a - 1j * k))
        assert proche(c_num(lambda s: mp.exp(a * s), k), exact)
        if k:   # leurres : sans (-1)^n (k impair), conjugué, facteur 2
            assert not proche(exact, mp.sinh(a * PI) / (PI * (a - 1j * k))) or k % 2 == 0
            assert not proche(exact, (-1) ** k * mp.sinh(a * PI) / (PI * (a + 1j * k)))
        assert not proche(exact, exact / 2)
    for k in range(0, 5):   # a_n = 2a(-1)^n sh(aπ)/(π(a² + n²))
        assert proche(a_num(lambda s: mp.exp(a * s), k), 2 * a * (-1) ** k * mp.sinh(a * PI) / (PI * (a * a + k * k)))
# v2 (a = -1) : c_n = (-1)^n sh(π)/(π(1+in)), et le leurre « -… » est faux.
assert proche(c_num(lambda s: mp.exp(-s), 2), mp.sinh(PI) / (PI * (1 + 2j)))

# int-l3-fourier-coefficients-ex4 (+ v1, v2) : Tchebychev T_4, T_5, T_6 et cos(nt) = T_n(cos t).
X = symbols("X")
T = {0: Integer(1), 1: X}
for k in range(1, 6):
    T[k + 1] = expand(2 * X * T[k] - T[k - 1])
assert T[4] == 8 * X**4 - 8 * X**2 + 1 and T[5] == 16 * X**5 - 20 * X**3 + 5 * X and T[6] == 32 * X**6 - 48 * X**4 + 18 * X**2 - 1
for k in range(7):
    assert expand(T[k] - chebyshevt(k, X)) == 0
    for s in [0.3, 1.1, 2.7]:
        assert abs(math.cos(k * s) - float(T[k].subs(X, math.cos(s)))) < 1e-12
for k in range(1, 7):
    assert T[k].as_poly(X).LC() == 2 ** (k - 1)
# Leurres : les polynômes proposés diffèrent de T_k (valeur en 1, en 1/2).
for vraie, leurres in [(T[5], [16 * X**5 + 20 * X**3 + 5 * X, 32 * X**5 - 40 * X**3 + 10 * X, 16 * X**5 - 20 * X**3 - 5 * X]),
                       (T[4], [8 * X**4 + 8 * X**2 + 1, 8 * X**4 - 8 * X**2 - 1, 16 * X**4 - 16 * X**2 + 2]),
                       (T[6], [32 * X**6 - 48 * X**4 + 18 * X**2 + 1, 64 * X**6 - 96 * X**4 + 36 * X**2 - 2, 32 * X**6 - 48 * X**4 - 18 * X**2 - 1])]:
    for l in leurres:
        assert expand(l - vraie) != 0

# ===========================================================================
# int-l3-fourier-approximation (carte)
# ===========================================================================
# Exemple : c_{±1}(dent) = ∓i, ‖f‖² = π²/3, distance² = π²/3 - 2.
assert proche(c_num(dent, 1), -1j) and proche(c_num(dent, -1), 1j)
assert integrate(t**2, (t, -pi, pi)) / (2 * pi) == pi**2 / 3
# Contre-exemple : f = |sin t|, sup|f - 2/π| = 2/π > 1/2 = sup|f - 1/2| (f prend toutes les valeurs de [0, 1]).
assert max(2 / PI, 1 - 2 / PI) == 2 / PI and 2 / PI > mp.mpf(1) / 2
# Trous : leurre L¹ — f = indicatrice de [0, π/2] (périodisée) : c_0 = 1/4, ∫|f - c| = π/2 + πc pour c ∈ [0,1],
# minimal en c = 0 ≠ 1/4.
L1 = lambda c: mp.quad(lambda s: abs((1 if 0 <= s <= PI / 2 else 0) - c), [-PI, 0, PI / 2, PI])
assert proche(L1(0), PI / 2, mp.mpf(10) ** -15) and L1(mp.mpf(1) / 4) > L1(0) + mp.mpf(1) / 2
# Leurres de la distance : f = e_{-1}, N = 1 : ‖f‖² - Σ_{n≥0}|c_n|² = 1 ≠ 0 ; f = 1 : ‖f‖² + |c_0|² = 2 ≠ 0 ;
# f = e_0 + e_1 : ‖f‖² - (Σ|c_n|)² = 2 - 4 = -2 ≠ 0. Leurre Σ c_n e_{-n} pour f = e_1 : distance² ‖e_1 - e_{-1}‖² = 2.
assert abs(mp.quad(lambda s: abs(mp.expj(s) - mp.expj(-s)) ** 2, [-PI, PI]) / (2 * PI) - 2) < EPS

# int-l3-fourier-approximation-ex1 (+ v1, v2) : S_1(dent) = 2 sin t, S_1(|t|) = π/2 - (4/π)cos t, S_2(créneau) = (4/π) sin t.
assert b_sym(t, 1) == 2 and a_sym(t, 0) == 0 and a_sym(t, 1) == 0
assert a_sym(Abs(t), 0) / 2 == pi / 2 and simplify(a_sym(Abs(t), 1) + 4 / pi) == 0 and b_sym(Abs(t), 1) == 0
assert proche(b_num(creneau, 1, [-PI, 0, PI]), 4 / PI) and proche(b_num(creneau, 2, [-PI, 0, PI]), 0)

# int-l3-fourier-approximation-ex2 (+ v1, v2) : distances².
d_dent = PI**2 / 3 - 2
d_cre = 1 - 8 / PI**2
d_abs = PI**2 / 12 - 8 / PI**2
assert proche(d_dent, mp.quad(lambda s: (s - 2 * mp.sin(s)) ** 2, [-PI, PI]) / (2 * PI))
assert proche(d_cre, mp.quad(lambda s: (mp.sign(s) - 4 / PI * mp.sin(s)) ** 2, [-PI, 0, PI]) / (2 * PI))
assert proche(d_abs, mp.quad(lambda s: (abs(s) - PI / 2 + 4 / PI * mp.cos(s)) ** 2, [-PI, 0, PI]) / (2 * PI))
distincts(d_dent, [PI**2 / 3 - 4, PI**2 / 3 - 1, PI**2 / 3])
distincts(d_cre, [1 - 16 / PI**2, 1 - 4 / PI**2, 8 / PI**2])
distincts(d_abs, [PI**2 / 3 - 8 / PI**2, PI**2 / 12 - 16 / PI**2, PI**2 / 12])
assert d_abs > 0 and d_cre > 0 and d_dent > 0

# int-l3-fourier-approximation-ex3 (+ v1, v2) : moindres carrés sur [0, π], minimisation exacte en SymPy.
A, B = symbols("A B", real=True)
from sympy import solve, diff
for g, u, v, a_att, b_att, m_att in [(Integer(1), sin(t), sin(3 * t), 4 / pi, 4 / (3 * pi), pi - 80 / (9 * pi)),
                                      (t, sin(t), sin(2 * t), Integer(2), Integer(-1), pi**3 / 3 - 5 * pi / 2),
                                      (t, Integer(1), cos(t), pi / 2, -4 / pi, pi**3 / 12 - 8 / pi)]:
    J = expand(integrate(expand((g - A * u - B * v) ** 2), (t, 0, pi)))
    sol = solve([diff(J, A), diff(J, B)], [A, B], dict=True)[0]
    assert simplify(sol[A] - a_att) == 0 and simplify(sol[B] - b_att) == 0
    assert simplify(J.subs(sol) - m_att) == 0
distincts(PI - 80 / (9 * PI), [PI - 8 / PI, PI - 160 / (9 * PI), 80 / (9 * PI)])
distincts(PI**3 / 3 - 5 * PI / 2, [PI**3 / 3 - 2 * PI, PI**3 / 3 - 5 * PI, PI**3 / 3 - 3 * PI / 2])
distincts(PI**3 / 12 - 8 / PI, [PI**3 / 12 - 4 / PI, PI**3 / 3 - 8 / PI, PI**3 / 12])

# int-l3-fourier-approximation-ex4 (+ v1, v2) : lim ∫_0^π g|sin(nt)| = (2/π)∫_0^π g ; contrôle numérique à n = 200.
# Développement admis de |sin u| : a_0 = 4/π, a_{2k} = -4/(π(4k²-1)) (vérifié dans la carte normale ci-dessous).
mp.mp.dps = 15
for g, lim, leurres in [(mp.exp, 2 / PI * (mp.e**PI - 1), [mp.e**PI - 1, (mp.e**PI - 1) / 2, 0]),
                        (lambda s: 1 / (1 + s), 2 / PI * mp.log(1 + PI), [mp.log(1 + PI), mp.log(1 + PI) / 2, 0]),
                        (lambda s: s * s, 2 * PI**2 / 3, [PI**3 / 3, PI**3 / 6, 0])]:
    nn = 200
    val = mp.fsum(mp.quad(lambda s: g(s) * abs(mp.sin(nn * s)), [k * PI / nn, (k + 1) * PI / nn]) for k in range(nn))
    assert abs(val - lim) < mp.mpf(10) ** -3, (val, lim)
    distincts(lim, leurres)
mp.mp.dps = 30

# ===========================================================================
# int-l3-fourier-parseval (carte)
# ===========================================================================
# Exemple : π²/3 = (1/2)Σ 4/n² ⇒ ζ(2) = π²/6.
assert proche(mp.zeta(2), PI**2 / 6)
# Trou 2, leurre |a_n + i b_n|² : pour f = e^{it}, a_1 = 1, b_1 = i ⇒ |a_1 + i b_1| = 0 alors que ‖f‖² = 1.
assert proche(a_num(lambda s: mp.expj(s), 1), 1) and proche(b_num(lambda s: mp.expj(s), 1), 1j)
# Affirmation 2 : f = 1 : a_0² = 4 ≠ 1.
assert a_sym(Integer(1), 0) ** 2 == 4
# Forme polarisée : f = dent, g = créneau : (1/2π)∫ t·signe(t) = π/2 = (1/2)Σ b_n(f) b_n(g) = (4/π)Σ_{n impair} 1/n².
assert proche(mp.quad(lambda s: s * mp.sign(s), [-PI, 0, PI]) / (2 * PI), PI / 2)
assert proche(4 / PI * mp.nsum(lambda k: 1 / (2 * k + 1) ** 2, [0, mp.inf]), PI / 2)

# int-l3-fourier-parseval-ex1 (+ v1, v2).
assert proche(mp.zeta(2), PI**2 / 6)
distincts(PI**2 / 6, [PI**2 / 8, PI**2 / 12, PI**2 / 3])
s2_imp = mp.nsum(lambda k: 1 / (2 * k + 1) ** 2, [0, mp.inf])
assert proche(1, mp.mpf(1) / 2 * 16 / PI**2 * s2_imp) and proche(s2_imp, PI**2 / 8)
distincts(PI**2 / 8, [PI**2 / 6, PI**2 / 4, PI / 4])
s4_imp = mp.nsum(lambda k: 1 / (2 * k + 1) ** 4, [0, mp.inf])
assert proche(PI**2 / 3, PI**2 / 4 + mp.mpf(1) / 2 * 16 / PI**2 * s4_imp) and proche(s4_imp, PI**4 / 96)
distincts(PI**4 / 96, [PI**4 / 90, PI**4 / 48, PI**2 / 8])

# int-l3-fourier-parseval-ex2 (+ v1, v2).
for k in range(1, 8):
    assert simplify(a_sym(t**2, k) - 4 * (-1) ** k / Integer(k) ** 2) == 0
    assert simplify(2 * integrate(t * (pi - t) * sin(k * t), (t, 0, pi)) / pi - 4 * (1 - (-1) ** k) / (pi * k**3)) == 0
    assert simplify(2 * integrate(cos(t / 2) * cos(k * t), (t, 0, pi)) / pi - 4 * (-1) ** (k + 1) / (pi * (4 * k**2 - 1))) == 0
assert a_sym(t**2, 0) == 2 * pi**2 / 3 and integrate(t**4, (t, -pi, pi)) / (2 * pi) == pi**4 / 5
assert simplify(2 * integrate(cos(t / 2), (t, 0, pi)) / pi - 4 / pi) == 0
assert integrate(t**2 * (pi - t) ** 2, (t, 0, pi)) / pi == pi**4 / 30
assert simplify(integrate(cos(t / 2) ** 2, (t, -pi, pi)) / (2 * pi) - Rational(1, 2)) == 0
assert proche(mp.zeta(4), PI**4 / 90) and proche(PI**4 / 5, PI**4 / 9 + 8 * mp.zeta(4))
distincts(PI**4 / 90, [PI**4 / 96, PI**4 / 45, PI**4 / 72])
s6_imp = mp.nsum(lambda k: 1 / (2 * k + 1) ** 6, [0, mp.inf])
assert proche(PI**4 / 30, mp.mpf(1) / 2 * 64 / PI**2 * s6_imp) and proche(s6_imp, PI**6 / 960)
distincts(PI**6 / 960, [PI**6 / 945, PI**6 / 480, PI**6 / 1920])
s_v2 = mp.nsum(lambda k: 1 / (4 * k * k - 1) ** 2, [1, mp.inf])
assert proche(mp.mpf(1) / 2, 4 / PI**2 + 8 / PI**2 * s_v2) and proche(s_v2, PI**2 / 16 - mp.mpf(1) / 2)
distincts(s_v2, [PI**2 / 8 - 1, PI**2 / 16 - mp.mpf(1) / 4, mp.mpf(1) / 2 - 4 / PI**2])

# int-l3-fourier-parseval-ex3 (+ v1, v2) : Σ_{n∈Z} 1/(a² + n²) = π coth(πa)/a par Parseval pour e^{at}.
for a in [1, 2, mp.mpf(1) / 2]:
    norme2 = mp.quad(lambda s: mp.exp(2 * a * s), [-PI, PI]) / (2 * PI)
    coef2 = lambda k: abs((-1) ** int(k) * mp.sinh(a * PI) / (PI * (a - 1j * k))) ** 2
    assert proche(norme2, mp.sinh(2 * a * PI) / (2 * a * PI))
    assert proche(mp.sinh(a * PI) ** 2 / PI**2 * (1 / a**2 + 2 * mp.nsum(lambda k: 1 / (a * a + k * k), [1, mp.inf])), norme2)
    assert proche(coef2(3), mp.sinh(a * PI) ** 2 / (PI**2 * (a * a + 9)))
S1 = mp.nsum(lambda k: 1 / (1 + k * k), [1, mp.inf])
assert proche(S1, (PI * mp.coth(PI) - 1) / 2)
distincts(S1, [PI * mp.coth(PI), PI * mp.coth(PI) / 2, (PI * mp.coth(PI) + 1) / 2])
S4 = mp.nsum(lambda k: 1 / (4 + k * k), [1, mp.inf])
assert proche(S4, PI / 4 * mp.coth(2 * PI) - mp.mpf(1) / 8)
distincts(S4, [PI / 2 * mp.coth(2 * PI), PI / 4 * mp.coth(2 * PI), PI / 2 * mp.coth(2 * PI) - mp.mpf(1) / 4])
Sq = mp.nsum(lambda k: 1 / (1 + 4 * k * k), [1, mp.inf])
assert proche(Sq, PI / 4 * mp.coth(PI / 2) - mp.mpf(1) / 2)
distincts(Sq, [PI / 2 * mp.coth(PI / 2), PI / 4 * mp.coth(PI / 2), PI / 2 * mp.coth(PI / 2) - 1])
assert proche(c_num(lambda s: mp.exp(s / 2), 3), 2 * (-1) ** 3 * mp.sinh(PI / 2) / (PI * (1 - 6j)))

# int-l3-fourier-parseval-ex4 (Wirtinger) : cas d'égalité a cos t + b sin t ; cos(2t) donne ∫f'² = 4∫f².
for f_, fp in [(cos(t), -sin(t)), (3 * cos(t) - 2 * sin(t), -3 * sin(t) - 2 * cos(t))]:
    assert simplify(integrate(f_**2, (t, -pi, pi)) - integrate(fp**2, (t, -pi, pi))) == 0
assert integrate(diff(cos(2 * t), t) ** 2, (t, -pi, pi)) == 4 * integrate(cos(2 * t) ** 2, (t, -pi, pi))
# Inégalité sur une fonction continue C¹ par morceaux de moyenne nulle : f = |t| - π/2.
assert integrate((Abs(t) - pi / 2) ** 2, (t, -pi, pi)) < integrate(Integer(1), (t, -pi, pi))

# ===========================================================================
# int-l3-fourier-dirichlet (carte)
# ===========================================================================
# Noyau : D_N(t) = sin((N+1/2)t)/sin(t/2), (1/2π)∫ D_N = 1.
for N in [1, 3, 6]:
    for s in [mp.mpf("0.3"), mp.mpf("1.7"), mp.mpf("-2.9")]:
        assert proche(mp.fsum(mp.cos(k * s) for k in range(-N, N + 1)), mp.sin((N + mp.mpf(1) / 2) * s) / mp.sin(s / 2))
# Exemple : Σ (-1)^k/(2k+1) = π/4 ; S_N(dent)(π) = 0.
assert proche(mp.nsum(lambda k: (-1) ** k / (2 * k + 1), [0, mp.inf]), PI / 4)
assert all(abs(mp.fsum(2 * (-1) ** (k + 1) / k * mp.sin(k * PI) for k in range(1, N))) < EPS for N in [5, 50])
# Gibbs (cours) : (2/π)Si(π) - 1 ≈ 0,179.
assert abs(2 / PI * mp.si(PI) - 1 - mp.mpf("0.17898")) < mp.mpf(10) ** -5

# int-l3-fourier-dirichlet-ex1 (+ v1, v2) : demi-sommes 1/2, ch π, π.
assert proche((mp.e**PI + mp.e**-PI) / 2, mp.cosh(PI))
distincts(mp.cosh(PI), [mp.e**PI, mp.e**-PI, mp.sinh(PI)])
distincts(PI, [0, 2 * PI, PI / 2])
# Contrôle direct de v1 : S_N(e^t)(π) = Σ c_n (-1)^n = (sh π/π) Σ 1/(1 - in) → ch π.
SN = lambda N: mp.sinh(PI) / PI * mp.fsum(1 / (1 - 1j * k) for k in range(-N, N + 1))
assert abs(SN(20000) - mp.cosh(PI)) < mp.mpf(10) ** -2

# int-l3-fourier-dirichlet-ex2 (+ v1, v2).
assert proche(mp.nsum(lambda k: (-1) ** k / (2 * k + 1), [0, mp.inf]), PI / 4)
distincts(PI / 4, [PI / 2, PI / 8, 1])
# v1 : ε_k = +,+,-,- ; somme par blocs de 4 (8m+1, 8m+3, 8m+5, 8m+7) = π√2/4 ; sin((2k+1)π/4) = ε_k √2/2.
blocs = mp.nsum(lambda m: 1 / (8 * m + 1) + 1 / (8 * m + 3) - 1 / (8 * m + 5) - 1 / (8 * m + 7), [0, mp.inf])
assert proche(blocs, PI * mp.sqrt(2) / 4)
for k in range(8):
    assert proche(mp.sin((2 * k + 1) * PI / 4), (1 if k % 4 in (0, 1) else -1) * mp.sqrt(2) / 2)
distincts(blocs, [PI / 4, PI * mp.sqrt(2) / 8, PI * mp.sqrt(2) / 2])
# v2 : η_n = 1 si n ≡ 1, -1 si n ≡ 2 (mod 3) ; (-1)^{n+1} sin(nπ/3) = η_n √3/2 ; somme = π√3/9.
for k in range(1, 13):
    eta = {1: 1, 2: -1, 0: 0}[k % 3]
    assert proche((-1) ** (k + 1) * mp.sin(k * PI / 3), eta * mp.sqrt(3) / 2)
blocs3 = mp.nsum(lambda m: 1 / (3 * m + 1) - 1 / (3 * m + 2), [0, mp.inf])
assert proche(blocs3, PI * mp.sqrt(3) / 9)
distincts(blocs3, [PI / 3, PI / 6, PI * mp.sqrt(3) / 6])

# int-l3-fourier-dirichlet-ex3 (+ v1, v2) : f = cos(αt), a_n = 2α(-1)^n sin(απ)/(π(α² - n²)).
al = symbols("alpha", positive=True)
for av in [Rational(1, 3), Rational(2, 3)]:
    for k in range(0, 5):
        assert simplify(2 * integrate(cos(av * t) * cos(k * t), (t, 0, pi)) / pi
                        - 2 * av * (-1) ** k * sin(av * pi) / (pi * (av**2 - k**2))) == 0
S_13 = mp.nsum(lambda k: (-1) ** (k + 1) / (9 * k * k - 1), [1, mp.inf])
assert proche(S_13, PI * mp.sqrt(3) / 9 - mp.mpf(1) / 2)
assert proche(1, 3 * mp.sqrt(3) / (2 * PI) + 3 * mp.sqrt(3) / PI * S_13)
distincts(S_13, [PI * mp.sqrt(3) / 9, mp.mpf(1) / 2 - PI * mp.sqrt(3) / 18, PI * mp.sqrt(3) / 18 - mp.mpf(1) / 2])
S_23 = mp.nsum(lambda k: (-1) ** (k + 1) / (9 * k * k - 4), [1, mp.inf])
assert proche(S_23, PI * mp.sqrt(3) / 18 - mp.mpf(1) / 8)
assert proche(1, 3 * mp.sqrt(3) / (4 * PI) + 6 * mp.sqrt(3) / PI * S_23)
distincts(S_23, [PI * mp.sqrt(3) / 18, PI * mp.sqrt(3) / 9 - mp.mpf(1) / 4, mp.mpf(1) / 8 - PI * mp.sqrt(3) / 36])
S_13pi = mp.nsum(lambda k: 1 / (9 * k * k - 1), [1, mp.inf])
assert proche(S_13pi, mp.mpf(1) / 2 - PI * mp.sqrt(3) / 18)
assert proche(mp.mpf(1) / 2, 3 * mp.sqrt(3) / (2 * PI) - 3 * mp.sqrt(3) / PI * S_13pi)
distincts(S_13pi, [PI * mp.sqrt(3) / 9 - mp.mpf(1) / 2, mp.mpf(1) / 2 - PI * mp.sqrt(3) / 9, PI * mp.sqrt(3) / 18])

# int-l3-fourier-dirichlet-ex4 (+ v1, v2) : Σ sin(na)/n = (π-a)/2 et Σ sin²(na)/n² = a(π-a)/2 (polylog).
for a in [mp.mpf(1), mp.mpf(2), mp.mpf(1) / 2]:
    s1 = mp.im(mp.polylog(1, mp.expj(a)))
    s2 = (mp.zeta(2) - mp.re(mp.polylog(2, mp.expj(2 * a)))) / 2
    assert proche(s1, (PI - a) / 2) and proche(s2, a * (PI - a) / 2)
    # c_n de l'indicatrice de [-a, a] : sin(na)/(πn), c_0 = a/π.
    assert proche(c_num(lambda s: 1 if abs(s) <= a else 0, 3, [-PI, -a, a, PI]), mp.sin(3 * a) / (3 * PI))
# Propositions (couples) : la bonne est la seule à valoir ((π-a)/2, a(π-a)/2).
for a, props in [(1, [((PI - 1) / 2, (PI - 1) / 2), ((PI - 1) / 2, (PI - 1) / 4), (PI / 2, (PI - 1) / 2), ((PI - 1) / 2, (PI**2 - 1) / 2)]),
                 (2, [((PI - 2) / 2, PI - 2), ((PI - 2) / 2, (PI - 2) / 2), (PI - 2, (PI - 2) / 2), ((PI - 2) / 2, 2 * PI - 4)]),
                 (mp.mpf(1) / 2, [((2 * PI - 1) / 4, (2 * PI - 1) / 8), ((2 * PI - 1) / 4, (2 * PI - 1) / 4),
                                  ((2 * PI - 1) / 8, (2 * PI - 1) / 4), ((PI - 1) / 4, (2 * PI - 1) / 8)])]:
    juste = ((PI - a) / 2, a * (PI - a) / 2)
    assert [proche(p[0], juste[0]) and proche(p[1], juste[1]) for p in props] == [True, False, False, False]

# int-l3-fourier-dirichlet-ex5 (+ v1, v2) : Gibbs, sommes de Riemann.
Si = mp.si(PI)
N = 3000
milieu = 2 / PI * PI / N * mp.fsum(mp.sinc((k + mp.mpf(1) / 2) * PI / N) for k in range(N))
assert abs(milieu - 2 / PI * Si) < mp.mpf(10) ** -6
cre_S2N = 4 / PI * mp.fsum(mp.sin((2 * k + 1) * PI / (2 * N)) / (2 * k + 1) for k in range(N))
assert proche(cre_S2N, milieu, mp.mpf(10) ** -20)
dent_S = 2 * mp.fsum((-1) ** (k + 1) * mp.sin(k * (PI - PI / N)) / k for k in range(1, N + 1))
assert abs(dent_S - 2 * Si) < 2 * PI / N   # extrémités droites : erreur O(1/N)
v2_S = mp.fsum(mp.sin(k * PI / N) / k for k in range(1, N + 1))
assert abs(v2_S - Si) < PI / N
assert abs(2 / PI * Si - mp.mpf("1.179")) < mp.mpf("0.0005") and abs(2 * Si - mp.mpf("3.704")) < mp.mpf("0.0005")
assert abs(Si - mp.mpf("1.852")) < mp.mpf("0.0005") and abs(Si / PI - mp.mpf("0.589")) < mp.mpf("0.0005")
assert abs(4 / PI * Si - mp.mpf("2.358")) < mp.mpf("0.0005") and 2 * Si > PI and Si > PI / 2
assert abs(2 * Si - PI - mp.mpf("0.562")) < mp.mpf("0.001") and abs(Si - PI / 2 - mp.mpf("0.281")) < mp.mpf("0.001")

# ===========================================================================
# int-l3-fourier-fejer (carte)
# ===========================================================================
for N in [1, 2, 5, 9]:
    for s in [mp.mpf("0.4"), mp.mpf("2.2"), mp.mpf("-1.3")]:
        Kc = mp.fsum((1 - mp.mpf(abs(k)) / N) * mp.cos(k * s) for k in range(-N + 1, N))
        Kd = mp.fsum(mp.sin((k + mp.mpf(1) / 2) * s) / mp.sin(s / 2) for k in range(N)) / N
        assert proche(Kc, mp.sin(N * s / 2) ** 2 / (N * mp.sin(s / 2) ** 2)) and proche(Kc, Kd) and Kc >= 0
    assert proche(mp.quad(lambda s: mp.sin(N * s / 2) ** 2 / (N * mp.sin(s / 2) ** 2), [-PI, 0, PI]) / (2 * PI), 1)
# Contre-exemple : D_1(π) = -1 ; exemple : S_5(créneau)(π/2) = 52/(15π) > 1.
assert 1 + 2 * mp.cos(PI) == -1 and 52 / (15 * PI) > 1
# Trou 3, leurre « sup|σ_N f - f| ≤ sup|f|/N » : f = cos 2t, N = 4 : σ_4 = cos(2t)/2, écart 1/2 > 1/4.
assert Fr(1, 2) > Fr(1, 4)
# Affirmation 1 (Fejér pour une fonction continue par morceaux) : créneau en 0, σ_N(0) = 0 = demi-somme.

# int-l3-fourier-fejer-ex1 (+ v1, v2).
def fejer_trig(coefs, N, s):
    """coefs : {n: (a_n, b_n)} pour n ≥ 1, et 0: a_0/2 ; renvoie σ_N en s."""
    tot = Fr(0)
    for k, val in coefs.items():
        if k == 0:
            tot += val
        elif k < N:
            tot += (1 - Fr(k, N)) * (val[0] * s[0][k] + val[1] * s[1][k])
    return tot
# valeurs exactes de cos(k s) et sin(k s) aux points utilisés
cs0 = ({k: 1 for k in range(6)}, {k: 0 for k in range(6)})
csp2 = ({1: 0, 2: -1, 3: 0}, {1: 1, 2: 0, 3: -1})
assert fejer_trig({0: Fr(1), 2: (4, 0), 3: (0, -2)}, 4, cs0) == 3
assert fejer_trig({1: (1, 0), 2: (1, 0)}, 3, cs0) == 1
assert fejer_trig({0: Fr(2), 1: (0, 1), 2: (3, 0)}, 2, csp2) == Fr(5, 2)
assert 3 not in (5, 2, 4) and 1 not in (2, Fr(5, 3), Fr(4, 3)) and Fr(5, 2) not in (0, 3, 1)
assert 1 + 4 - 0 == 5 and 2 + 1 - 3 == 0     # f(0) et f(π/2), leurres

# int-l3-fourier-fejer-ex2 (+ v1, v2) : (1/2π)∫K_N² = Σ(1-|n|/N)² = (2N²+1)/(3N).
for N, val in [(3, Fr(19, 9)), (2, Fr(3, 2)), (4, Fr(11, 4))]:
    assert sum((1 - Fr(abs(k), N)) ** 2 for k in range(-N + 1, N)) == val == Fr(2 * N * N + 1, 3 * N)
    assert sum(1 - Fr(abs(k), N) for k in range(-N + 1, N)) == N
    assert proche(mp.quad(lambda s: (mp.sin(N * s / 2) ** 2 / (N * mp.sin(s / 2) ** 2)) ** 2, [-PI, 0, PI]) / (2 * PI),
                  mp.mpf(val.numerator) / val.denominator)
assert Fr(19, 9) not in (3, 1, Fr(14, 9)) and Fr(3, 2) not in (2, 1, Fr(5, 4)) and Fr(11, 4) not in (4, 1, Fr(15, 8))

# int-l3-fourier-fejer-ex3 (+ v1, v2) : créneau, S_N et σ_N en un point.
bcre = lambda k: (4 / (PI * k)) if k % 2 else 0
S_cre = lambda N, s: mp.fsum(bcre(k) * mp.sin(k * s) for k in range(1, N + 1))
sig_cre = lambda N, s: mp.fsum((1 - mp.mpf(k) / N) * bcre(k) * mp.sin(k * s) for k in range(1, N))
assert proche(S_cre(5, PI / 2), 52 / (15 * PI)) and proche(sig_cre(6, PI / 2), 14 / (5 * PI))
assert proche(S_cre(3, PI / 4), 8 * mp.sqrt(2) / (3 * PI)) and proche(sig_cre(4, PI / 4), 5 * mp.sqrt(2) / (3 * PI))
assert proche(S_cre(5, PI / 6), 56 / (15 * PI)) and proche(sig_cre(6, PI / 6), 12 / (5 * PI))
for N, s in [(6, PI / 2), (4, PI / 4), (6, PI / 6)]:     # σ_N par la convolution avec K_N (|σ_N| ≤ 1)
    conv = mp.quad(lambda u: mp.sign(s - u) * mp.sin(N * u / 2) ** 2 / (N * mp.sin(u / 2) ** 2), [s - PI, 0, s, PI + s]) / (2 * PI)
    assert proche(conv, sig_cre(N, s), mp.mpf(10) ** -15) and abs(conv) <= 1
distincts(14 / (5 * PI), [52 / (15 * PI), 1, 4 / PI])
distincts(5 * mp.sqrt(2) / (3 * PI), [8 * mp.sqrt(2) / (3 * PI), 1, 2 * mp.sqrt(2) / PI])
distincts(12 / (5 * PI), [56 / (15 * PI), 1, 2 / PI])
assert S_cre(5, PI / 2) > 1 and S_cre(3, PI / 4) > 1 and S_cre(5, PI / 6) > 1

# int-l3-fourier-fejer-ex4 : L_N ≥ (4/π²) H_N, L_N → +∞ ; (1/2π)∫|K_N| = 1.
for N in [1, 5, 20]:
    L = mp.quad(lambda s: abs(mp.sin((N + mp.mpf(1) / 2) * s) / mp.sin(s / 2)), mp.linspace(0, PI, 2 * N + 2)) / PI
    borne = 2 / PI * mp.quad(lambda u: abs(mp.sin(u)) / u, mp.linspace(0, (N + mp.mpf(1) / 2) * PI, 2 * N + 2))
    assert L >= borne - EPS and borne >= 4 / PI**2 * mp.harmonic(N) - EPS
assert mp.harmonic(10**6) > 13

# int-l3-fourier-fejer-ex5 : preuve (Cesàro, Fejér) ; pas de calcul.

# ===========================================================================
# int-l3-fourier-normale (carte)
# ===========================================================================
# Exemple : |t|, c_n(f') (créneau) = -2i/(πn) = in·c_n(|t|) pour n impair.
for k in [-5, -3, -1, 1, 3, 5]:
    cf = c_num(lambda s: abs(s), k, [-PI, 0, PI])
    assert proche(cf, -2 / (PI * k * k)) and proche(1j * k * cf, -2j / (PI * k))
    assert proche(c_num(creneau, k, [-PI, 0, PI]), 1j * k * cf)
# Contre-exemple : dent, in·c_n = -(-1)^n alors que c_n(f') = 0.
for k in range(1, 5):
    assert proche(1j * k * c_num(dent, k), -(-1) ** k)
# Inégalité Σ_{n≠0}|c_n| ≤ (π/√3)‖f'‖_2 sur |t| et |sin t|, égalité pour t² (exercice ex4).
for cf, fp2 in [(lambda k: 2 / (PI * k * k) if k % 2 else 0, mp.mpf(1)),
                (lambda k: 2 / (PI * (k * k - 1)) if k % 2 == 0 else 0, mp.mpf(1) / 2)]:
    assert 2 * mp.nsum(lambda k: cf(int(k)), [1, mp.inf]) <= PI / mp.sqrt(3) * mp.sqrt(fp2)

# int-l3-fourier-normale-ex1 : |sin| continue et C¹ par morceaux ; autres discontinues (pas de calcul).

# int-l3-fourier-normale-ex2 (+ v1, v2) : coefficients de |sin t| et sommes.
for k in range(0, 9):
    att = 0 if k % 2 else -4 / (pi * (k * k - 1))
    assert simplify(2 * integrate(sin(t) * cos(k * t), (t, 0, pi)) / pi - att) == 0
f_absin = lambda s: 2 / PI - 4 / PI * mp.nsum(lambda k: mp.cos(2 * k * s) / (4 * k * k - 1), [1, mp.inf])
for s in [0, mp.mpf("0.7"), PI / 3, 2]:
    assert proche(f_absin(s), abs(mp.sin(s)), mp.mpf(10) ** -15)
Sa = mp.nsum(lambda k: (-1) ** (k + 1) / (4 * k * k - 1), [1, mp.inf])
assert proche(Sa, (PI - 2) / 4)
distincts(Sa, [(PI - 2) / 2, (PI + 2) / 4, mp.mpf(1) / 2])
Sb = mp.nsum(lambda m: (-1) ** (m + 1) / (16 * m * m - 1), [1, mp.inf])
assert proche(Sb, PI * mp.sqrt(2) / 8 - mp.mpf(1) / 2)
distincts(Sb, [PI * mp.sqrt(2) / 4 - mp.mpf(1) / 2, mp.mpf(1) / 2 - PI * mp.sqrt(2) / 8, PI * mp.sqrt(2) / 8])
Sc = mp.nsum(lambda k: mp.cos(2 * k * PI / 3) / (4 * k * k - 1), [1, mp.inf])
assert proche(Sc, mp.mpf(1) / 2 - PI * mp.sqrt(3) / 8, mp.mpf(10) ** -15)
distincts(Sc, [PI * mp.sqrt(3) / 8 - mp.mpf(1) / 2, mp.mpf(1) / 2 - PI * mp.sqrt(3) / 4, PI * mp.sqrt(3) / 8])
assert proche(mp.nsum(lambda k: 1 / (4 * k * k - 1), [1, mp.inf]), mp.mpf(1) / 2)

# int-l3-fourier-normale-ex3 (+ v1, v2) : vitesses de convergence uniforme.
for N in [10**3, 10**5]:
    reste_abs = 4 / PI * mp.zeta(2, N // 2 + (N % 2) + mp.mpf(1) / 2) / 4   # impairs > N (zêta de Hurwitz)
    assert abs(reste_abs * PI * N / 2 - 1) < mp.mpf(3) / N
    reste_t2 = 4 * mp.zeta(2, N + 1)
    assert abs(reste_t2 * N / 4 - 1) < mp.mpf(2) / N
for N in [0, 1, 4, 30]:   # |sin t| : (4/π)Σ_{k>N} 1/(4k²-1) = 2/(π(2N+1)) exactement
    assert proche(4 / PI * mp.nsum(lambda k: 1 / (4 * k * k - 1), [N + 1, mp.inf]), 2 / (PI * (2 * N + 1)))
for N in [3, 4]:          # le sup est atteint en 0 (|t|) et en π (t²) : écart calculé directement
    S_abs0 = PI / 2 - 4 / PI * mp.fsum(1 / mp.mpf(k * k) for k in range(1, N + 1, 2))
    assert proche(S_abs0, 4 / PI * mp.nsum(lambda j: 1 / (2 * j + 1) ** 2, [(N + 1) // 2, mp.inf]))
    S_t2pi = PI**2 / 3 + 4 * mp.fsum(mp.mpf((-1) ** k) / k**2 * (-1) ** k for k in range(1, N + 1))
    assert proche(PI**2 - S_t2pi, 4 * mp.nsum(lambda k: 1 / k**2, [N + 1, mp.inf]))
# Leurres : équivalents 4/(πN), 2/(πN²), 1/N ; 4/(π(2N+1)), 2/(π(2N+1)²), 1/(2N+1) ; 2/N, 4/N², 8/N (différents à N = 10).
distincts(2 / (PI * 10), [4 / (PI * 10), 2 / (PI * 100), mp.mpf(1) / 10])
distincts(2 / (PI * 21), [4 / (PI * 21), 2 / (PI * 441), mp.mpf(1) / 21])
distincts(mp.mpf(4) / 10, [mp.mpf(2) / 10, mp.mpf(4) / 100, mp.mpf(8) / 10])

# int-l3-fourier-normale-ex4 : constante optimale π/√3, atteinte pour t².
assert integrate(t**2, (t, -pi, pi)) / (2 * pi) == pi**2 / 3
assert simplify(sqrt(integrate((2 * t) ** 2, (t, -pi, pi)) / (2 * pi)) - 2 * pi / sqrt(3)) == 0
assert simplify(pi / sqrt(3) * 2 * pi / sqrt(3) - 2 * pi**2 / 3) == 0 and pi**2 - pi**2 / 3 == 2 * pi**2 / 3
assert proche(2 * mp.nsum(lambda k: 2 / k**2, [1, mp.inf]), 2 * PI**2 / 3)     # Σ_{n≠0}|c_n(t²)| = 2π²/3
distincts(PI / mp.sqrt(3), [PI / mp.sqrt(6), 1, 2 * PI / mp.sqrt(3)])

# ===========================================================================
# int-l3-fourier-weierstrass (carte)
# ===========================================================================
# Contre-exemple : |π - P(π)| + |-π - P(π)| ≥ 2π (inégalité triangulaire) ; QCM q7 : ∫ g cos t = 0 pour g ⊥ cos t.
for k in range(2, 6):
    assert integrate(cos(k * t) * cos(t), (t, -pi, pi)) == 0 and integrate(sin(k * t) * cos(t), (t, -pi, pi)) == 0
assert integrate(cos(t) * cos(t), (t, -pi, pi)) == pi

# int-l3-fourier-weierstrass-ex1 (+ v1, v2).
assert simplify(-I / 2 * exp(2 * I * t) + I / 2 * exp(-2 * I * t) - sin(2 * t)).rewrite(exp).simplify() == 0
assert cos(pi / 3) == Rational(1, 2) and sin(pi / 4) == sqrt(2) / 2 and 1 + cos(pi / 2) == 1
assert Rational(1, 2) not in (1, sqrt(3) / 2, 0) and sqrt(2) / 2 not in (-sqrt(2) / 2, 1, 0)

# int-l3-fourier-weierstrass-ex2 : cos((2n+1)π/2) = 0 ; sin(0) = 0.
assert all(cos((2 * k + 1) * pi / 2) == 0 for k in range(6))

# int-l3-fourier-weierstrass-ex3 (+ v1, v2) : polynômes P tels que S_N(f)(t) = P(cos t).
assert simplify(2 * integrate(Abs(cos(t)), (t, 0, pi)) / pi - 4 / pi) == 0
assert simplify(2 * integrate(Abs(cos(t)) * cos(t), (t, 0, pi)) / pi) == 0
assert simplify(4 * integrate(cos(t) * cos(2 * t), (t, 0, pi / 2)) / pi - 4 / (3 * pi)) == 0
assert proche(a_num(lambda s: abs(mp.cos(s)), 2, [-PI, -PI / 2, PI / 2, PI]), 4 / (3 * PI))
P1 = (2 + 8 * x**2) / (3 * pi)
assert simplify(P1.subs(x, cos(t)) - (2 / pi + 4 / (3 * pi) * cos(2 * t))) == 0
P2 = (10 - 8 * x**2) / (3 * pi)
assert simplify(P2.subs(x, cos(t)) - (2 / pi - 4 / (3 * pi) * cos(2 * t))) == 0
assert simplify(2 * integrate(sin(t) * cos(2 * t), (t, 0, pi)) / pi + 4 / (3 * pi)) == 0
P3 = pi / 2 - 8 * x / (3 * pi) - 16 * x**3 / (9 * pi)
assert simplify(P3.subs(x, cos(t)) - (pi / 2 - 4 / pi * cos(t) - 4 / (9 * pi) * cos(3 * t))) == 0
assert simplify(a_sym(Abs(t), 3) + 4 / (9 * pi)) == 0 and a_sym(Abs(t), 2) == 0
for vraie, leurres in [(P1, [2 / pi + 8 * x**2 / (3 * pi), (2 - 8 * x**2) / (3 * pi), (4 + 8 * x**2) / (3 * pi)]),
                       (P2, [2 / pi - 8 * x**2 / (3 * pi), (10 + 8 * x**2) / (3 * pi), (2 - 8 * x**2) / (3 * pi)]),
                       (P3, [pi / 2 - 4 * x / pi - 16 * x**3 / (9 * pi), pi - 8 * x / (3 * pi) - 16 * x**3 / (9 * pi),
                             pi / 2 + 8 * x / (3 * pi) + 16 * x**3 / (9 * pi)])]:
    for l in leurres:
        assert simplify(l - vraie) != 0

# int-l3-fourier-weierstrass-ex4 (+ v1, v2) : moyennes de Weyl le long des entiers.
Nw = 200000
for g, lim, leurres in [(lambda k: abs(math.sin(k)), 2 / math.pi, [0.5, math.pi / 4, 0]),
                        (lambda k: math.sin(k) ** 4, 3 / 8, [0.25, 0.5, 0]),
                        (lambda k: abs(math.cos(k)) ** 3, 4 / (3 * math.pi), [2 / math.pi, 2 / 3, 0])]:
    moy = math.fsum(g(k) for k in range(1, Nw + 1)) / Nw
    assert abs(moy - lim) < 1e-3, (moy, lim)
    assert all(abs(lim - l) > 0.05 for l in leurres)
assert proche(mp.quad(lambda s: abs(mp.cos(s)) ** 3, [-PI, -PI / 2, PI / 2, PI]) / (2 * PI), 4 / (3 * PI))
assert simplify(sin(t) ** 4 - (Rational(3, 8) - cos(2 * t) / 2 + cos(4 * t) / 8)) == 0

# ===========================================================================
# QCM int-l3-fourier-q1 … q8
# ===========================================================================
assert a_sym(sin(t) ** 2, 2) == Rational(-1, 2) and a_sym(sin(t) ** 2, 0) == 1                       # q1
assert simplify(integrate(sin(t) ** 2 * exp(-2 * I * t), (t, -pi, pi)) / (2 * pi) + Rational(1, 4)) == 0
assert 2 / PI > mp.mpf(1) / 2                                                                           # q2
assert 1 + 2 * Fr(1, 4) / (1 - Fr(1, 4)) == Fr(5, 3) and 1 + Fr(1, 4) / (1 - Fr(1, 4)) == Fr(4, 3)     # q3
assert proche(mp.quad(lambda s: abs(mp.fsum(mp.mpf(2) ** -abs(k) * mp.expj(k * s) for k in range(-60, 61))) ** 2,
                      [-PI, PI]) / (2 * PI), mp.mpf(5) / 3)
assert 1 + 2 * mp.nsum(lambda k: mp.mpf(2) ** -k, [1, mp.inf]) == 3
assert all(abs(mp.fsum(2 * (-1) ** (k + 1) / k * mp.sin(k * PI) for k in range(1, N))) < EPS for N in [3, 30])  # q4
assert (1 - Fr(2, 4)) == Fr(1, 2) and (Fr(0) + 0 + 1 + 1) / 4 == Fr(1, 2)                               # q5
# q6 : créneau, |c_n| = 2/(π|n|) (voir carte coefficients) ; q7 : orthogonalité ci-dessus ; q8 : modification en un point.

# ===========================================================================
# Relecture indépendante (2026-11-08) : chaque bonne réponse est recalculée par un autre chemin (intégrales numériques
# des coefficients, sommes directes de séries, polylogarithmes), et chaque leurre est vérifié FAUX (écart > 1e-6),
# avec mpmath à 30 chiffres.
# ===========================================================================
def juste_et_leurres(vraie, bonne, leurres, tol=mp.mpf(10) ** -15):
    """La bonne proposition vaut la valeur recalculée ; chaque leurre en diffère."""
    assert abs(mp.mpc(vraie) - mp.mpc(bonne)) < tol, (vraie, bonne)
    for l in leurres:
        assert abs(mp.mpc(vraie) - mp.mpc(l)) > mp.mpf(10) ** -6, (vraie, l)


SH, CH = mp.sinh, mp.cosh
COTH = lambda u: mp.cosh(u) / mp.sinh(u)
# coefficients-ex1, v1, v2 : a_3(cos³), b_3(sin³), a_0(cos⁴) par quadrature.
juste_et_leurres(a_num(lambda s: mp.cos(s) ** 3, 3), mp.mpf(1) / 4, [mp.mpf(3) / 4, mp.mpf(1) / 8, mp.mpf(1) / 2])
juste_et_leurres(b_num(lambda s: mp.sin(s) ** 3, 3), -mp.mpf(1) / 4, [mp.mpf(1) / 4, -mp.mpf(3) / 4, mp.mpf(3) / 4])
juste_et_leurres(a_num(lambda s: mp.cos(s) ** 4, 0), mp.mpf(3) / 4, [mp.mpf(3) / 8, mp.mpf(3) / 2, mp.mpf(1) / 8])


# coefficients-ex2 (+ v1, v2) : formules en n, comparées pour n = 1..4 ; chaque leurre faux pour au moins un n.
def formule_juste(calc, bonne, leurres, ns=range(1, 5)):
    vals = [calc(k) for k in ns]
    assert all(abs(v - bonne(k)) < mp.mpf(10) ** -15 for v, k in zip(vals, ns))
    for l in leurres:
        assert any(abs(v - l(k)) > mp.mpf(10) ** -6 for v, k in zip(vals, ns)), l


formule_juste(lambda k: b_num(dent, k), lambda k: 2 * (-1) ** (k + 1) / mp.mpf(k),
              [lambda k: 2 / mp.mpf(k), lambda k: 2 * (-1) ** k / mp.mpf(k), lambda k: (-1) ** (k + 1) / mp.mpf(k)])
formule_juste(lambda k: mp.quad(lambda s: s * mp.sin(k * s), [0, 2 * PI]) / PI, lambda k: -2 / mp.mpf(k),   # v1
              [lambda k: 2 / mp.mpf(k), lambda k: 2 * (-1) ** (k + 1) / mp.mpf(k), lambda k: -1 / mp.mpf(k)])
formule_juste(lambda k: a_num(abs, k, [-PI, 0, PI]), lambda k: 2 * ((-1) ** k - 1) / (PI * k ** 2),            # v2
              [lambda k: 2 * (1 - (-1) ** k) / (PI * k ** 2), lambda k: ((-1) ** k - 1) / (PI * k ** 2),
               lambda k: 2 * ((-1) ** k - 1) / (PI * k)])
# coefficients-ex3 (+ v1, v2) : c_n(e^{αt}), n = -2..2.
for al, bonne, leurres in [
        (1, lambda k: (-1) ** k * SH(PI) / (PI * (1 - 1j * k)),
         [lambda k: SH(PI) / (PI * (1 - 1j * k)), lambda k: (-1) ** k * SH(PI) / (PI * (1 + 1j * k)),
          lambda k: (-1) ** k * SH(PI) / (2 * PI * (1 - 1j * k))]),
        (2, lambda k: (-1) ** k * SH(2 * PI) / (PI * (2 - 1j * k)),
         [lambda k: SH(2 * PI) / (PI * (2 - 1j * k)), lambda k: (-1) ** k * SH(2 * PI) / (PI * (2 + 1j * k)),
          lambda k: (-1) ** k * SH(2 * PI) / (2 * PI * (2 - 1j * k))]),
        (-1, lambda k: (-1) ** k * SH(PI) / (PI * (1 + 1j * k)),
         [lambda k: (-1) ** k * SH(PI) / (PI * (1 - 1j * k)), lambda k: -(-1) ** k * SH(PI) / (PI * (1 + 1j * k)),
          lambda k: SH(PI) / (PI * (1 + 1j * k))])]:
    formule_juste(lambda k: c_num(lambda s: mp.exp(al * s), k), bonne, leurres, ns=range(-2, 3))
# coefficients-ex4 (+ v1, v2) : T_k(cos u) = cos(ku) en trois points ; chaque leurre faux en l'un d'eux.
for k, bonne, leurres in [
        (5, lambda X: 16 * X ** 5 - 20 * X ** 3 + 5 * X,
         [lambda X: 16 * X ** 5 + 20 * X ** 3 + 5 * X, lambda X: 32 * X ** 5 - 40 * X ** 3 + 10 * X,
          lambda X: 16 * X ** 5 - 20 * X ** 3 - 5 * X]),
        (4, lambda X: 8 * X ** 4 - 8 * X ** 2 + 1,
         [lambda X: 8 * X ** 4 + 8 * X ** 2 + 1, lambda X: 8 * X ** 4 - 8 * X ** 2 - 1, lambda X: 16 * X ** 4 - 16 * X ** 2 + 2]),
        (6, lambda X: 32 * X ** 6 - 48 * X ** 4 + 18 * X ** 2 - 1,
         [lambda X: 32 * X ** 6 - 48 * X ** 4 + 18 * X ** 2 + 1, lambda X: 64 * X ** 6 - 96 * X ** 4 + 36 * X ** 2 - 2,
          lambda X: 32 * X ** 6 - 48 * X ** 4 - 18 * X ** 2 - 1])]:
    formule_juste(lambda u: mp.cos(k * u), lambda u: bonne(mp.cos(u)), [lambda u, l=l: l(mp.cos(u)) for l in leurres],
                  ns=[mp.mpf(1) / 3, mp.mpf(1), mp.mpf(2)])


# approximation-ex1 (+ v1, v2) : chaque proposition P donne ∫|f - P|² ; la bonne donne strictement le plus petit.
def ecart2(f, P, pts):
    return mp.quad(lambda s: abs(f(s) - P(s)) ** 2, pts)


for f, pts, bonne, leurres in [
        (dent, [-PI, PI], lambda s: 2 * mp.sin(s), [mp.sin, lambda s: -2 * mp.sin(s), lambda s: 4 / PI * mp.sin(s)]),
        (abs, [-PI, 0, PI], lambda s: PI / 2 - 4 / PI * mp.cos(s),
         [lambda s: PI - 4 / PI * mp.cos(s), lambda s: PI / 2 + 4 / PI * mp.cos(s), lambda s: PI / 2 - 2 / PI * mp.cos(s)]),
        (creneau, [-PI, 0, PI], lambda s: 4 / PI * mp.sin(s),
         [lambda s: 4 / PI * mp.sin(s) + 2 / PI * mp.sin(2 * s), lambda s: 2 / PI * mp.sin(s), mp.sin])]:
    eb = ecart2(f, bonne, pts)
    assert all(ecart2(f, l, pts) > eb + mp.mpf(10) ** -6 for l in leurres)
# approximation-ex2 (+ v1, v2) : d² = ‖f‖² - Σ_{|n|≤1}|c_n|², coefficients par quadrature.
for f, pts, bonne, leurres in [
        (dent, [-PI, PI], PI ** 2 / 3 - 2, [PI ** 2 / 3 - 4, PI ** 2 / 3 - 1, PI ** 2 / 3]),
        (creneau, [-PI, 0, PI], 1 - 8 / PI ** 2, [1 - 16 / PI ** 2, 1 - 4 / PI ** 2, 8 / PI ** 2]),
        (abs, [-PI, 0, PI], PI ** 2 / 12 - 8 / PI ** 2, [PI ** 2 / 3 - 8 / PI ** 2, PI ** 2 / 12 - 16 / PI ** 2, PI ** 2 / 12])]:
    d2 = mp.quad(lambda s: abs(f(s)) ** 2, pts) / (2 * PI) - mp.fsum(abs(c_num(f, k, pts)) ** 2 for k in (-1, 0, 1))
    juste_et_leurres(d2, bonne, leurres)


# approximation-ex3 (+ v1, v2) : moindres carrés sur [0, π] par le système normal (matrice de Gram numérique).
def moindres_carres(g, u, v):
    G = mp.matrix([[mp.quad(lambda s: a(s) * b(s), [0, PI]) for b in (u, v)] for a in (u, v)])
    B = mp.matrix([mp.quad(lambda s: g(s) * a(s), [0, PI]) for a in (u, v)])
    sol = mp.lu_solve(G, B)
    return mp.quad(lambda s: (g(s) - sol[0] * u(s) - sol[1] * v(s)) ** 2, [0, PI]), sol


m, sol = moindres_carres(lambda s: 1, mp.sin, lambda s: mp.sin(3 * s))
juste_et_leurres(m, PI - 80 / (9 * PI), [PI - 8 / PI, PI - 160 / (9 * PI), 80 / (9 * PI)])
assert proche(sol[0], 4 / PI) and proche(sol[1], 4 / (3 * PI))
m, sol = moindres_carres(lambda s: s, mp.sin, lambda s: mp.sin(2 * s))                                     # v1
juste_et_leurres(m, PI ** 3 / 3 - 5 * PI / 2, [PI ** 3 / 3 - 2 * PI, PI ** 3 / 3 - 5 * PI, PI ** 3 / 3 - 3 * PI / 2])
assert proche(sol[0], 2) and proche(sol[1], -1)
m, sol = moindres_carres(lambda s: s, lambda s: 1, mp.cos)                                                # v2
juste_et_leurres(m, PI ** 3 / 12 - 8 / PI, [PI ** 3 / 12 - 4 / PI, PI ** 3 / 3 - 8 / PI, PI ** 3 / 12])
assert proche(sol[0], PI / 2) and proche(sol[1], -4 / PI)
# approximation-ex4 (+ v1, v2) : limite (2/π)∫_0^π g, et contrôle direct à n = 300 (intégrale découpée en arches).
for g, bonne, leurres in [(mp.exp, 2 / PI * (mp.e ** PI - 1), [mp.e ** PI - 1, (mp.e ** PI - 1) / 2, 0]),
                          (lambda s: 1 / (1 + s), 2 / PI * mp.log(1 + PI), [mp.log(1 + PI), mp.log(1 + PI) / 2, 0]),
                          (lambda s: s ** 2, 2 * PI ** 2 / 3, [PI ** 3 / 3, PI ** 3 / 6, 0])]:
    juste_et_leurres(2 / PI * mp.quad(g, [0, PI]), bonne, leurres)
    with mp.workdps(15):
        n0 = 300
        I_n = mp.quad(lambda s: g(s) * abs(mp.sin(n0 * s)), mp.linspace(0, PI, n0 + 1))
        assert abs(I_n - bonne) < mp.mpf(10) ** -3 * max(1, abs(bonne)), (I_n, bonne)

# parseval-ex1, ex2, ex3 (+ v1, v2) : sommes directes (nsum) comparées aux propositions.
Z = lambda F: mp.nsum(F, [1, mp.inf])
juste_et_leurres(Z(lambda k: 1 / k ** 2), PI ** 2 / 6, [PI ** 2 / 8, PI ** 2 / 12, PI ** 2 / 3])
juste_et_leurres(Z(lambda k: 1 / (2 * k - 1) ** 2), PI ** 2 / 8, [PI ** 2 / 6, PI ** 2 / 4, PI / 4])
juste_et_leurres(Z(lambda k: 1 / (2 * k - 1) ** 4), PI ** 4 / 96, [PI ** 4 / 90, PI ** 4 / 48, PI ** 2 / 8])
juste_et_leurres(Z(lambda k: 1 / k ** 4), PI ** 4 / 90, [PI ** 4 / 96, PI ** 4 / 45, PI ** 4 / 72])
juste_et_leurres(Z(lambda k: 1 / (2 * k - 1) ** 6), PI ** 6 / 960, [PI ** 6 / 945, PI ** 6 / 480, PI ** 6 / 1920])
juste_et_leurres(Z(lambda k: 1 / (4 * k ** 2 - 1) ** 2), PI ** 2 / 16 - mp.mpf(1) / 2,
                 [PI ** 2 / 8 - 1, PI ** 2 / 16 - mp.mpf(1) / 4, mp.mpf(1) / 2 - 4 / PI ** 2])
juste_et_leurres(Z(lambda k: 1 / (1 + k ** 2)), (PI * COTH(PI) - 1) / 2,
                 [PI * COTH(PI), PI * COTH(PI) / 2, (PI * COTH(PI) + 1) / 2])
juste_et_leurres(Z(lambda k: 1 / (4 + k ** 2)), PI / 4 * COTH(2 * PI) - mp.mpf(1) / 8,
                 [PI / 2 * COTH(2 * PI), PI / 4 * COTH(2 * PI), PI / 2 * COTH(2 * PI) - mp.mpf(1) / 4])
juste_et_leurres(Z(lambda k: 1 / (1 + 4 * k ** 2)), PI / 4 * COTH(PI / 2) - mp.mpf(1) / 2,
                 [PI / 2 * COTH(PI / 2), PI / 4 * COTH(PI / 2), PI / 2 * COTH(PI / 2) - 1])
# parseval-ex2 : coefficients a_n(t²) = 4(-1)^n/n² et b_n de v1 = 8/(πn³) (n impair) par quadrature, n = 1..4.
formule_juste(lambda k: a_num(lambda s: s ** 2, k), lambda k: 4 * (-1) ** k / mp.mpf(k) ** 2, [])
formule_juste(lambda k: b_num(lambda s: s * (PI - abs(s)), k, [-PI, 0, PI]),
              lambda k: 4 * (1 - (-1) ** k) / (PI * k ** 3), [])
formule_juste(lambda k: a_num(lambda s: mp.cos(s / 2), k), lambda k: 4 * (-1) ** (k + 1) / (PI * (4 * k ** 2 - 1)), [],
              ns=range(0, 4))
# parseval-ex4 : cos(2t) vérifie les hypothèses sans égalité ; a cos t + b sin t réalise l'égalité.
assert proche(mp.quad(lambda s: (-2 * mp.sin(2 * s)) ** 2, [-PI, PI]), 4 * mp.quad(lambda s: mp.cos(2 * s) ** 2, [-PI, PI]))
assert proche(mp.quad(lambda s: (2 * mp.cos(s) - 3 * mp.sin(s)) ** 2, [-PI, PI]),
              mp.quad(lambda s: (-2 * mp.sin(s) - 3 * mp.cos(s)) ** 2, [-PI, PI]))

# dirichlet-ex1 (+ v1, v2) : demi-somme des limites latérales ; v1 contrôlé par la série Σ c_n(-1)^n (convergente).
juste_et_leurres(mp.mpf(1) / 2, mp.mpf(1) / 2, [1, 0])
cn_exp = lambda k: (-1) ** k * SH(PI) / (PI * (1 - 1j * k))
S_pi = cn_exp(0) + Z(lambda k: (cn_exp(k) + cn_exp(-k)) * (-1) ** k)
juste_et_leurres(S_pi, CH(PI), [mp.e ** PI, mp.e ** -PI, SH(PI)])
juste_et_leurres((0 + 2 * PI) / 2, PI, [0, 2 * PI, PI / 2])
# dirichlet-ex2 (+ v1, v2) : sommes directes, regroupées par périodes pour accélérer.
juste_et_leurres(Z(lambda k: 1 / (4 * k - 3) - 1 / (4 * k - 1)), PI / 4, [PI / 2, PI / 8, 1])
juste_et_leurres(Z(lambda k: 1 / (8 * k - 7) + 1 / (8 * k - 5) - 1 / (8 * k - 3) - 1 / (8 * k - 1)), PI * mp.sqrt(2) / 4,
                 [PI / 4, PI * mp.sqrt(2) / 8, PI * mp.sqrt(2) / 2])
juste_et_leurres(Z(lambda k: 1 / (3 * k - 2) - 1 / (3 * k - 1)), PI * mp.sqrt(3) / 9,
                 [PI / 3, PI / 6, PI * mp.sqrt(3) / 6])
# dirichlet-ex3 (+ v1, v2) : sommes directes.
juste_et_leurres(Z(lambda k: (-1) ** (k + 1) / (9 * k ** 2 - 1)), PI * mp.sqrt(3) / 9 - mp.mpf(1) / 2,
                 [PI * mp.sqrt(3) / 9, mp.mpf(1) / 2 - PI * mp.sqrt(3) / 18, PI * mp.sqrt(3) / 18 - mp.mpf(1) / 2])
juste_et_leurres(Z(lambda k: (-1) ** (k + 1) / (9 * k ** 2 - 4)), PI * mp.sqrt(3) / 18 - mp.mpf(1) / 8,
                 [PI * mp.sqrt(3) / 18, PI * mp.sqrt(3) / 9 - mp.mpf(1) / 4, mp.mpf(1) / 8 - PI * mp.sqrt(3) / 36])
juste_et_leurres(Z(lambda k: 1 / (9 * k ** 2 - 1)), mp.mpf(1) / 2 - PI * mp.sqrt(3) / 18,
                 [PI * mp.sqrt(3) / 9 - mp.mpf(1) / 2, mp.mpf(1) / 2 - PI * mp.sqrt(3) / 9, PI * mp.sqrt(3) / 18])
# dirichlet-ex4 (+ v1, v2) : Σ sin(na)/n = Im(-log(1 - e^{ia})), Σ sin²(na)/n² = (ζ(2) - Re Li₂(e^{2ia}))/2.
for a, bonne, leurres in [
        (1, ((PI - 1) / 2, (PI - 1) / 2), [((PI - 1) / 2, (PI - 1) / 4), (PI / 2, (PI - 1) / 2), ((PI - 1) / 2, (PI ** 2 - 1) / 2)]),
        (2, ((PI - 2) / 2, PI - 2), [((PI - 2) / 2, (PI - 2) / 2), (PI - 2, (PI - 2) / 2), ((PI - 2) / 2, 2 * PI - 4)]),
        (mp.mpf(1) / 2, ((2 * PI - 1) / 4, (2 * PI - 1) / 8),
         [((2 * PI - 1) / 4, (2 * PI - 1) / 4), ((2 * PI - 1) / 8, (2 * PI - 1) / 4), ((PI - 1) / 4, (2 * PI - 1) / 8)])]:
    s1 = mp.im(-mp.log(1 - mp.expj(a)))
    s2 = (mp.zeta(2) - mp.re(mp.polylog(2, mp.expj(2 * a)))) / 2
    assert proche(s1, bonne[0]) and proche(s2, bonne[1])
    assert all(abs(s1 - l1) + abs(s2 - l2) > mp.mpf(10) ** -6 for l1, l2 in leurres)
# dirichlet-ex5 (+ v1, v2) : sommes partielles calculées directement pour N = 4000, proches de la limite (Si(π)) ;
# les leurres en sont éloignés de plus de 0,01.
SI = mp.si(PI)
with mp.workdps(20):
    N0 = 4000
    v0 = 4 / PI * mp.fsum(mp.sin((2 * k + 1) * PI / (2 * N0)) / (2 * k + 1) for k in range(N0))
    v1 = 2 * mp.fsum((-1) ** (k + 1) * mp.sin(k * (PI - PI / N0)) / k for k in range(1, N0 + 1))
    v2 = mp.fsum(mp.sin(k * PI / N0) / k for k in range(1, N0 + 1))
for v, bonne, leurres in [(v0, 2 / PI * SI, [1, SI / PI, 4 / PI * SI]), (v1, 2 * SI, [PI, SI, 2 / PI * SI]),
                          (v2, SI, [PI / 2, 2 / PI * SI, 2 * SI])]:
    assert abs(v - bonne) < 2 * mp.mpf(10) ** -3, (v, bonne)            # erreur en O(1/N) (sommes de Riemann)
    assert all(abs(bonne - l) > mp.mpf(10) ** -2 for l in leurres)
assert abs(2 / PI * SI - mp.mpf("1.17898")) < 1e-5 and abs(2 * SI - mp.mpf("3.70387")) < 1e-5

# fejer-ex1 (+ v1, v2) et fejer-ex3 (+ v1, v2), q5 : σ_N(f)(x) recalculé comme moyenne des S_k(f)(x),
# coefficients obtenus par quadrature.
def sigma(f, Nf, x0, pts=None):
    cs = {k: c_num(f, k, pts) for k in range(-Nf, Nf + 1)}
    S = lambda K: mp.fsum(cs[k] * mp.expj(k * x0) for k in range(-K, K + 1))
    return mp.re(mp.fsum(S(K) for K in range(Nf)) / Nf), mp.re(S(Nf - 1))


juste_et_leurres(sigma(lambda s: 1 + 4 * mp.cos(2 * s) - 2 * mp.sin(3 * s), 4, 0)[0], 3, [5, 2, 4])
juste_et_leurres(sigma(lambda s: mp.cos(s) + mp.cos(2 * s), 3, 0)[0], 1, [2, mp.mpf(5) / 3, mp.mpf(4) / 3])
juste_et_leurres(sigma(lambda s: 2 + mp.sin(s) + 3 * mp.cos(2 * s), 2, PI / 2)[0], mp.mpf(5) / 2, [0, 3, 1])
juste_et_leurres(sigma(lambda s: mp.cos(2 * s), 4, 0)[0], mp.mpf(1) / 2, [1, mp.mpf(3) / 4, mp.mpf(1) / 4])   # q5
for Nf, x0, bonne_sig, bonne_S, leurres in [
        (6, PI / 2, 14 / (5 * PI), 52 / (15 * PI), [52 / (15 * PI), 1, 4 / PI]),
        (4, PI / 4, 5 * mp.sqrt(2) / (3 * PI), 8 * mp.sqrt(2) / (3 * PI), [8 * mp.sqrt(2) / (3 * PI), 1, 2 * mp.sqrt(2) / PI]),
        (6, PI / 6, 12 / (5 * PI), 56 / (15 * PI), [56 / (15 * PI), 1, 2 / PI])]:
    sg, sp = sigma(creneau, Nf, x0, [-PI, 0, PI])
    juste_et_leurres(sg, bonne_sig, leurres)
    assert proche(sp, bonne_S) and sp > 1 > sg                          # la somme partielle dépasse 1, pas Fejér
# fejer-ex2 (+ v1, v2) : (1/2π)∫K_N² avec K_N = (1/N)(sin(Nt/2)/sin(t/2))².
for Nf, bonne, leurres in [(3, mp.mpf(19) / 9, [3, 1, mp.mpf(14) / 9]), (2, mp.mpf(3) / 2, [2, 1, mp.mpf(5) / 4]),
                           (4, mp.mpf(11) / 4, [4, 1, mp.mpf(15) / 8])]:
    KN = lambda s: (mp.sin(Nf * s / 2) / mp.sin(s / 2)) ** 2 / Nf
    juste_et_leurres(mp.quad(lambda s: KN(s) ** 2, [-PI, 0, PI]) / (2 * PI), bonne, leurres)
# fejer-ex4 : L_N croît (et dépasse les leurres finis 1, 4/π², 0) ; minoration (4/π²)Σ1/k vérifiée pour N = 5, 20.
with mp.workdps(15):
    LN = lambda Nf: mp.quad(lambda s: abs(mp.sin((Nf + mp.mpf(1) / 2) * s) / mp.sin(s / 2)),
                            mp.linspace(0, PI, 2 * Nf + 2)) / PI
    L5, L20 = LN(5), LN(20)
    assert L20 > L5 > 1 and L5 >= 4 / PI ** 2 * mp.fsum(mp.mpf(1) / k for k in range(1, 6))
    assert L20 >= 4 / PI ** 2 * mp.fsum(mp.mpf(1) / k for k in range(1, 21))

# normale-ex2 (+ v1, v2) : sommes directes.
juste_et_leurres(Z(lambda k: (-1) ** (k + 1) / (4 * k ** 2 - 1)), (PI - 2) / 4, [(PI - 2) / 2, (PI + 2) / 4, mp.mpf(1) / 2])
juste_et_leurres(Z(lambda k: (-1) ** (k + 1) / (16 * k ** 2 - 1)), PI * mp.sqrt(2) / 8 - mp.mpf(1) / 2,
                 [PI * mp.sqrt(2) / 4 - mp.mpf(1) / 2, mp.mpf(1) / 2 - PI * mp.sqrt(2) / 8, PI * mp.sqrt(2) / 8])
juste_et_leurres(Z(lambda k: mp.cos(2 * k * PI / 3) / (4 * k ** 2 - 1)), mp.mpf(1) / 2 - PI * mp.sqrt(3) / 8,
                 [PI * mp.sqrt(3) / 8 - mp.mpf(1) / 2, mp.mpf(1) / 2 - PI * mp.sqrt(3) / 4, PI * mp.sqrt(3) / 8])
# normale-ex3 (+ v2) : N·sup|f - S_N| tend vers 2/π (|t|) et 4 (t²), pas vers les constantes des leurres ;
# v1 : valeur exacte 2/(π(2N+1)), testée pour N = 1..5 contre la somme directe des |a_{2k}|.
for Nf in [10 ** 4, 10 ** 5]:
    r1 = Nf * 4 / PI * mp.nsum(lambda k: 1 / (2 * k + 1) ** 2, [(Nf + 1) // 2, mp.inf])
    r2 = Nf * 4 * mp.nsum(lambda k: 1 / k ** 2, [Nf + 1, mp.inf])
    assert abs(r1 - 2 / PI) < 10 / mp.mpf(Nf) and abs(r2 - 4) < 10 / mp.mpf(Nf)
    assert all(abs(r1 - l) > mp.mpf(10) ** -2 for l in [4 / PI, 1]) and all(abs(r2 - l) > 1 for l in [2, 8])
for Nf in range(1, 6):
    vraie = 4 / PI * mp.nsum(lambda k: 1 / (4 * k ** 2 - 1), [Nf + 1, mp.inf])
    juste_et_leurres(vraie, 2 / (PI * (2 * Nf + 1)), [4 / (PI * (2 * Nf + 1)), 2 / (PI * (2 * Nf + 1) ** 2), mp.mpf(1) / (2 * Nf + 1)])
# normale-ex4 : pour t², sup|f - c_0| / ‖f'‖₂ = π/√3 (donc π/√6 et 1 ne sont pas des constantes valables).
c0 = mp.quad(lambda s: s ** 2, [-PI, PI]) / (2 * PI)
supf = max(abs(PI ** 2 - c0), abs(0 - c0))
normd = mp.sqrt(mp.quad(lambda s: (2 * s) ** 2, [-PI, PI]) / (2 * PI))
assert proche(supf / normd, PI / mp.sqrt(3)) and PI / mp.sqrt(6) < supf / normd and 1 < supf / normd

# weierstrass-ex3 (+ v1, v2) : a_n par quadrature, P(x) recomposé et comparé aux propositions en x = 0, 1/2, 1.
for g, Nf, bonne, leurres in [
        (lambda s: abs(mp.cos(s)), 2, lambda X: (2 + 8 * X ** 2) / (3 * PI),
         [lambda X: 2 / PI + 8 * X ** 2 / (3 * PI), lambda X: (2 - 8 * X ** 2) / (3 * PI), lambda X: (4 + 8 * X ** 2) / (3 * PI)]),
        (lambda s: abs(mp.sin(s)), 2, lambda X: (10 - 8 * X ** 2) / (3 * PI),
         [lambda X: 2 / PI - 8 * X ** 2 / (3 * PI), lambda X: (10 + 8 * X ** 2) / (3 * PI), lambda X: (2 - 8 * X ** 2) / (3 * PI)]),
        (lambda s: mp.acos(mp.cos(s)), 3, lambda X: PI / 2 - 8 * X / (3 * PI) - 16 * X ** 3 / (9 * PI),
         [lambda X: PI / 2 - 4 * X / PI - 16 * X ** 3 / (9 * PI), lambda X: PI - 8 * X / (3 * PI) - 16 * X ** 3 / (9 * PI),
          lambda X: PI / 2 + 8 * X / (3 * PI) + 16 * X ** 3 / (9 * PI)])]:
    an = [a_num(g, k, [-PI, -PI / 2, 0, PI / 2, PI]) for k in range(Nf + 1)]
    SN = lambda u: an[0] / 2 + mp.fsum(an[k] * mp.cos(k * u) for k in range(1, Nf + 1))
    formule_juste(lambda X: SN(mp.acos(X)), bonne, leurres, ns=[mp.mpf(0), mp.mpf(1) / 2, mp.mpf(1)])
# weierstrass-ex4 (+ v1, v2) : moyennes M_N(g) pour N = 2·10⁵, et moyennes intégrales (1/2π)∫g.
with mp.workdps(15):
    Nw = 200000
    ks = range(1, Nw + 1)
    for g, bonne, leurres in [(lambda u: abs(math.sin(u)), 2 / PI, [mp.mpf(1) / 2, PI / 4, 0]),
                              (lambda u: math.sin(u) ** 4, mp.mpf(3) / 8, [mp.mpf(1) / 4, mp.mpf(1) / 2, 0]),
                              (lambda u: abs(math.cos(u)) ** 3, 4 / (3 * PI), [2 / PI, mp.mpf(2) / 3, 0])]:
        moy = math.fsum(g(k) for k in ks) / Nw
        assert abs(moy - bonne) < 1e-3 and all(abs(bonne - l) > 0.05 for l in leurres)
juste_et_leurres(mp.quad(lambda s: abs(mp.sin(s)), [-PI, 0, PI]) / (2 * PI), 2 / PI, [mp.mpf(1) / 2, PI / 4, 0])

# QCM : q1 (a_2(sin²) par quadrature), q3 (Σ 4^{-|n|}).
juste_et_leurres(a_num(lambda s: mp.sin(s) ** 2, 2), -mp.mpf(1) / 2, [mp.mpf(1) / 2, -mp.mpf(1) / 4, 0])
juste_et_leurres(1 + 2 * Z(lambda k: mp.mpf(4) ** -k), mp.mpf(5) / 3, [mp.mpf(4) / 3, 3, 9])

print("int_l3_fourier : toutes les vérifications passent.")
