"""Vérifications calculatoires du chapitre « num-l2-equations » (méthodes numériques, L2) :
cartes, exercices et QCM de content/*/methodes-numeriques/l2-num-l2-equations.json et content/qcm/num-l2-equations.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variantes v1, v2) ou de la question.
Lancer : python3 verification/num_l2_equations.py (doit finir sans erreur, en quelques secondes).
"""

import math
from fractions import Fraction as F

import mpmath
from sympy import Rational, diff, factor, limit, simplify, sqrt, symbols

mpmath.mp.dps = 50
x = symbols("x")


def dicho(f, a, b, n):
    """Segments [a_k, b_k], k = 1..n, de la dichotomie (on garde [a, c] si f(a) f(c) <= 0)."""
    a, b = F(a), F(b)
    out = []
    for _ in range(n):
        c = (a + b) / 2
        if f(a) * f(c) <= 0:
            b = c
        else:
            a = c
        out.append((a, b))
    return out


def newton(f, fp, x0, n, exact=True):
    xs = [F(x0) if exact else mpmath.mpf(x0)]
    for _ in range(n):
        xs.append(xs[-1] - f(xs[-1]) / fp(xs[-1]))
    return xs


# ---------------------------------------------------------------------------------------------------------------
# num-l2-equations-dichotomie (carte)
# ---------------------------------------------------------------------------------------------------------------
# Exemple : x^2 - 2 sur [1, 2], quatre étapes, puis n = 19 pour 10^-6.
assert dicho(lambda t: t * t - 2, 1, 2, 4) == [(1, F(3, 2)), (F(5, 4), F(3, 2)), (F(11, 8), F(3, 2)), (F(11, 8), F(23, 16))]
assert min(n for n in range(60) if F(1, 2 ** (n + 1)) <= F(1, 10**6)) == 19 and 2**20 == 1048576
# Contre-exemple sans continuité : f(0) = -1, f = 1 sur ]0,1] : a_n = 0, b_n = 2^-n.
saut = lambda t: -1 if t == 0 else 1
assert all(s == (0, F(1, 2 ** (k + 1))) for k, s in enumerate(dicho(saut, 0, 1, 30)))
# Contre-exemple sans changement de signe : (x - 1/3)^2, b_n = 1 et a_n -> 1.
segs = dicho(lambda t: (t - F(1, 3)) ** 2, 0, 1, 40)
assert all(b == 1 for _, b in segs) and 1 - segs[-1][0] < F(1, 10**12)
# Plusieurs zéros : x^3 - x converge vers 1 sur [-2,3] et vers -1 sur [-3,2] (aussi exercice dichotomie-ex3 et v1).
cube = lambda t: t**3 - t
assert abs(dicho(cube, -2, 3, 80)[-1][0] - 1) < F(1, 10**20) and dicho(cube, -2, 3, 1)[0] == (F(1, 2), 3)
assert abs(dicho(cube, -3, 2, 80)[-1][0] + 1) < F(1, 10**20) and dicho(cube, -3, 2, 1)[0] == (-3, F(-1, 2))
# Trou 4, leurres : f(x) = x sur [0,1] donne c_n = 2^-(n+1) et l = 0 (> 2^-(n+2), > 4^-n pour n = 2) ; sur [0,4], c_n = 4/2^(n+1) > 1/2^n.
for k, (a, b) in enumerate(dicho(lambda t: t, 0, 1, 10)):
    assert a == 0 and (a + b) / 2 == F(1, 2 ** (k + 2))  # c_{k+1}
assert F(1, 8) > F(1, 16)  # n = 2 : |c_2 - l| = 1/8 > 1/4^2
assert all((F(4, 2 ** (k + 1))) > F(1, 2**k) for k in range(20))
# Trou 3 : (b-a)/(2n) vaut (b-a)/2^n pour n = 1, 2 mais pas n = 3 ; (b-a)/n et (b-a)/n^2 diffèrent pour n = 3.
assert F(1, 6) != F(1, 8) and F(1, 3) != F(1, 8) and F(1, 9) != F(1, 8)
# Affirmations : x^2 - 1/4 sur [-1,1] ; |c_{n+1} - c_n| = (b-a)/2^(n+2) ; x - 1/2 : c_0 = 1/2, c_1 = 1/4.
assert (1 - F(1, 4)) ** 2 > 0
segs = [(F(-2), F(3))] + dicho(cube, -2, 3, 20)
assert all(abs((a2 + b2) / 2 - (a1 + b1) / 2) == F(5, 2 ** (k + 2)) for k, ((a1, b1), (a2, b2)) in enumerate(zip(segs, segs[1:])))
segs = dicho(lambda t: t - F(1, 2), 0, 1, 60)
assert segs[0] == (0, F(1, 2)) and (segs[0][0] + segs[0][1]) / 2 == F(1, 4) and abs(segs[-1][0] - F(1, 2)) < F(1, 10**15)

# num-l2-equations-dichotomie-ex1 (+ v1, v2)
assert dicho(lambda t: t**3 + t - 1, 0, 1, 3) == [(F(1, 2), 1), (F(1, 2), F(3, 4)), (F(5, 8), F(3, 4))]
assert F(1, 2) ** 3 + F(1, 2) - 1 == F(-3, 8) and F(3, 4) ** 3 + F(3, 4) - 1 == F(11, 64) and F(5, 8) ** 3 + F(5, 8) - 1 == F(-67, 512)
assert dicho(lambda t: t**3 - 2, 1, 2, 3) == [(1, F(3, 2)), (F(5, 4), F(3, 2)), (F(5, 4), F(11, 8))]
assert F(5, 4) ** 3 - 2 == F(-3, 64) and F(11, 8) ** 3 - 2 == F(307, 512)
ex = lambda t: mpmath.exp(t) - 3
assert dicho(ex, 0, 2, 3) == [(1, 2), (1, F(3, 2)), (1, F(5, 4))]
assert abs(mpmath.exp(1.5) - 3 - 1.48) < 0.01 and abs(mpmath.exp(1.25) - 3 - 0.49) < 0.01
for (a, b), r in [((F(5, 8), F(3, 4)), 0.6823278), ((F(5, 4), F(11, 8)), 2 ** (1 / 3)), ((1, F(5, 4)), math.log(3))]:
    assert a < r < b

# num-l2-equations-dichotomie-ex2 (+ v1, v2) : plus petit n avec (b-a)/2^n <= eps.
for L, eps, n, leurres in [(1, F(1, 10**6), 20, [19, 6, 21]), (3, F(1, 10**4), 15, [14, 4, 16]), (2, F(1, 10**8), 28, [27, 8, 29])]:
    assert min(k for k in range(100) if F(L, 2**k) <= eps) == n and n not in leurres
    assert F(L, 2 ** (n - 1)) / 2 < eps  # le milieu c_{n-1} approche déjà le zéro à eps près

# num-l2-equations-dichotomie-ex3 v2 : x^3 - 4x sur [-3,4] -> 2.
c4 = lambda t: t**3 - 4 * t
assert dicho(c4, -3, 4, 1)[0] == (F(1, 2), 4) and abs(dicho(c4, -3, 4, 80)[-1][0] - 2) < F(1, 10**20)

# num-l2-equations-dichotomie-ex4 (+ v1, v2) : un milieu tombe sur un zéro.
s = dicho(cube, -2, 2, 60)
assert s[0] == (-2, 0) and s[1] == (-2, -1) and all(b == -1 for _, b in s[1:]) and abs(s[-1][0] + 1) < F(1, 10**15)
s = dicho(c4, -4, 4, 60)
assert s[0] == (-4, 0) and s[1] == (-4, -2) and all(b == -2 for _, b in s[1:]) and abs(s[-1][0] + 2) < F(1, 10**15)
s = dicho(lambda t: mpmath.sin(t), -4, 4, 60)
assert s[:3] == [(-4, 0), (-4, -2), (-4, -3)] and abs(s[-1][0] + mpmath.pi) < 1e-15
assert mpmath.sin(-4) > 0 and mpmath.sin(-2) < 0 and mpmath.sin(-3) < 0

# ---------------------------------------------------------------------------------------------------------------
# num-l2-equations-point-fixe (carte)
# ---------------------------------------------------------------------------------------------------------------
# Exemple : cos sur [0,1] ; cos 1 ~ 0,54 ; sin 1 ~ 0,841 ; l ~ 0,739085 ; n = 81 pour (sin 1)^n <= 1e-6.
assert 0 < math.cos(1) < 1 and abs(math.sin(1) - 0.841) < 1e-3
ell_cos = mpmath.findroot(lambda t: mpmath.cos(t) - t, 0.7)
assert abs(ell_cos - mpmath.mpf("0.739085")) < 1e-6
assert min(n for n in range(200) if math.sin(1) ** n <= 1e-6) == 81
# Contre-exemples : x + 1/x sur [1, +oo[ (rapport 1 - 1/(xy) dans [0,1[) ; sqrt(1+x^2) > x ; -x/2 alterne.
for a, b in [(1, 2), (F(3, 2), 7), (10, 11), (1, F(101, 100))]:
    g = lambda t: t + 1 / F(t)
    assert abs(g(a) - g(b)) < abs(F(a) - b) and g(a) >= 2
assert all(math.sqrt(1 + t * t) > t for t in [-5, 0, 3, 1e6])
assert all(abs(t) / math.sqrt(1 + t * t) < 1 for t in [-5, 0, 3, 1e3])
assert F(5, 2) == F(1) / 2 + 2 and 2 == F(0) / 2 + 2 and 4 == 4 / 2 + 2  # x/2 + 2 : g([0,1]) = [2, 5/2], point fixe 4
# Trou 4 : avec g(x) = x/2 sur R, u_0 = 1 : |u_n| = 2^-n > 4^-n (n >= 1), > 2^-(n+1) ; n k |u_0| = 0 pour n = 0.
assert all(F(1, 2**n) > F(1, 4**n) for n in range(1, 30))

# num-l2-equations-point-fixe-ex1 (+ v1, v2) : sqrt(c + x) sur [0, B], k = 1/(2 sqrt c), l = B.
for c, B, other, k in [(2, 2, -1, sqrt(2) / 4), (6, 3, -2, sqrt(6) / 12), (12, 4, -3, sqrt(3) / 12)]:
    assert sqrt(c + B) == B and B * B - B - c == 0 and other**2 - other - c == 0
    assert simplify(1 / (2 * sqrt(c)) - k) == 0 and float(k) < 1 and sqrt(c) >= 0
    assert simplify(1 / (2 * sqrt(c + B)) - Rational(1, 2 * B)) == 0  # leurre g'(B) = 1/(2B)
    assert simplify(2 * k - 1 / sqrt(c)) == 0 and 2 * k != k  # leurre 2k

# num-l2-equations-point-fixe-ex2 (+ v1, v2) : plus petit n avec k^n/(1-k) |u_1 - u_0| <= 1e-6.
for g, k, u1, coef, n, leurres, l in [
    (lambda t: 1 / (2 + t), F(1, 4), F(1, 2), F(2, 3), 10, [9, 11, 20], math.sqrt(2) - 1),
    (lambda t: 1 / (3 + t), F(1, 9), F(1, 3), F(3, 8), 6, [5, 7, 12], (-3 + math.sqrt(13)) / 2),
    (lambda t: 2 / (4 + t), F(1, 8), F(1, 2), F(4, 7), 7, [6, 8, 14], math.sqrt(6) - 2),
]:
    assert g(F(0)) == u1 and u1 / (1 - k) == coef and 0 <= g(F(1)) < g(F(0)) <= 1
    assert min(m for m in range(100) if coef * k**m <= F(1, 10**6)) == n and n not in leurres
    assert abs(g(l) - l) < 1e-12
    # |g'| <= k sur [0,1] : maximum en 0
    assert abs((g(F(1, 10**9)) - g(F(0))) * 10**9) <= k + F(1, 10**6)

# num-l2-equations-point-fixe-ex3 (+ v1, v2) : (x^2 + ab)/(a+b) converge vers a ssi u_0 dans ]-b, b[.
for a, b in [(1, 3), (2, 4), (1, 5)]:
    g = lambda t: (t * t + a * b) / F(a + b)
    assert factor(x**2 - (a + b) * x + a * b) == (x - a) * (x - b)
    assert F(2 * a, a + b) < 1 < F(2 * b, a + b)
    for u0 in [F(-b) + F(1, 100), F(-a), F(0), F(a) + F(1, 3), F(b) - F(1, 100), F(b), F(-b), F(b) + F(1, 100)]:
        u = u0
        for _ in range(3000):
            u = g(u)
            if u > 10**6:
                break
        if abs(u0) < b:
            assert abs(u - a) < F(1, 10**6), (a, b, u0)
        elif abs(u0) == b:
            assert u == b
        else:
            assert u > 10**6

# num-l2-equations-point-fixe-ex4 : contre-exemple x + 1/x (déjà vérifié avec la carte).

# ---------------------------------------------------------------------------------------------------------------
# num-l2-equations-newton (carte)
# ---------------------------------------------------------------------------------------------------------------
heron = lambda a, n: newton(lambda t: t * t - a, lambda t: 2 * t, 2, n)
xs = heron(2, 4)
assert xs[:4] == [2, F(3, 2), F(17, 12), F(577, 408)]
errs = [mpmath.mpf(t.numerator) / t.denominator - mpmath.sqrt(2) for t in xs]
assert abs(errs[1] - 0.0858) < 1e-4 and abs(errs[2] - 0.00245) < 1e-5 and abs(errs[3] - 2.12e-6) < 1e-8
assert all(errs[k + 1] <= errs[k] ** 2 / 2 for k in range(4))  # constante M/(2m) = 1/2 sur [1,2]
# Contre-exemple « mauvais côté » : x^3 - 1 sur [1/10, 3/2], x_0 = 1/10 -> x_1 = 33,4.
x1 = newton(lambda t: t**3 - 1, lambda t: 3 * t * t, F(1, 10), 1)[1]
assert x1 == F(334, 10)
# Trou 2, leurre milieu : x_0 = 0,8 donne x_1 > x_0.
assert newton(lambda t: t**3 - 1, lambda t: 3 * t * t, F(8, 10), 1)[1] > F(8, 10)
# Contre-exemple arctan : x_0 = 2 -> -3,54, puis |x_n| croît en changeant de signe.
xa = newton(mpmath.atan, lambda t: 1 / (1 + t * t), 2, 6, exact=False)
assert abs(xa[1] + 3.5357) < 1e-3 and all(abs(xa[k + 1]) > abs(xa[k]) and xa[k + 1] * xa[k] < 0 for k in range(6))
assert all(abs(2 * t / (1 + t * t) ** 2) <= 0.65 for t in [i / 100 for i in range(-300, 301)])  # leurre |f''| <= 1
# Leurre f'' <= 0 : ln sur [1/2, 3], x_0 = 3 -> 3 - 3 ln 3 < 1/2.
assert 3 - 3 * math.log(3) < 0.5
# f décroissante : e^-x - 1/2, x_0 = 0 -> 1/2, et ln 2 > 1/2.
assert newton(lambda t: mpmath.exp(-t) - 0.5, lambda t: -mpmath.exp(-t), 0, 1, exact=False)[1] == 0.5 and math.log(2) > 0.5
# Cycle : x^3 - 2x + 2, x_0 = 0 -> 1 -> 0 ; zéro réel ~ -1,77 avec f' > 0.
cyc = newton(lambda t: t**3 - 2 * t + 2, lambda t: 3 * t * t - 2, 0, 4)
assert cyc == [0, 1, 0, 1, 0]
ell_c = mpmath.findroot(lambda t: t**3 - 2 * t + 2, -1.8)
assert abs(ell_c + 1.769) < 1e-3 and 3 * ell_c**2 - 2 > 0
# Trou 4, leurres : (x_n - l)^3 (x^2 - 2, n = 1) ; m/(2M) (x^3 - 2, n = 0) ; M/(4m) (x^2 - 2 sur [1,4 ; 1,5], x_0 = 1,5).
assert errs[2] > errs[1] ** 3 / 2
e0 = 2 - mpmath.cbrt(2)
assert F(3, 2) - mpmath.cbrt(2) > F(3, 24) * e0**2
e0 = mpmath.mpf(1.5) - mpmath.sqrt(2)
e1 = mpmath.mpf(17) / 12 - mpmath.sqrt(2)
assert e1 > 2 / (4 * 2.8) * e0**2 and 1.4**2 < 2 < 1.5**2
# Affirmation : x^2 - 2, x_0 = 1 -> 3/2.
assert newton(lambda t: t * t - 2, lambda t: 2 * t, 1, 1)[1] == F(3, 2)

# num-l2-equations-newton-ex1 (+ v1, v2) : Héron pour a = 2, 3, 5.
for a, x1, x2, x3, err in [(2, F(3, 2), F(17, 12), F(577, 408), 2.2e-6), (3, F(7, 4), F(97, 56), F(18817, 10864), 2.5e-9),
                           (5, F(9, 4), F(161, 72), F(51841, 23184), 4.2e-10)]:
    assert heron(a, 3)[1:] == [x1, x2, x3]
    assert 0 < mpmath.mpf(x3.numerator) / x3.denominator - mpmath.sqrt(a) < err
assert F(17, 12) not in (F(3, 2), F(577, 408), F(7, 5)) and F(97, 56) != F(5, 3) and F(161, 72) != F(11, 5)

# num-l2-equations-newton-ex2 : x_{n+1} = x_n (2 - a x_n), r_{n+1} = r_n^2.
a_, X_ = symbols("a X")
assert simplify(X_ - (1 / X_ - a_) / (-1 / X_**2) - X_ * (2 - a_ * X_)) == 0
assert simplify(1 - a_ * X_ * (2 - a_ * X_) - (1 - a_ * X_) ** 2) == 0
for x0, ok in [(F(1, 10), True), (F(19, 10), True), (F(2), False), (F(3), False)]:  # a = 1
    t = x0
    for _ in range(12):
        t = t * (2 - t)
        if t <= 0:
            break
    assert (t > 0 and abs(t - 1) < F(1, 10**6)) == ok

# num-l2-equations-newton-ex3 (+ v1, v2) : majorants de x_4 - l.
for f, fp, l, x1, adm, C, rep in [
    (lambda t: t**3 - 2, lambda t: 3 * t * t, mpmath.cbrt(2), F(3, 2), F(1, 4), F(2), F(1, 512)),
    (lambda t: t**3 - 5, lambda t: 3 * t * t, mpmath.cbrt(5), F(7, 4), F(1, 20), F(2), F(5, 10**9)),
    (lambda t: t**2 - 3, lambda t: 2 * t, mpmath.sqrt(3), F(7, 4), F(1, 50), F(1, 2), F(2, 10**16)),
]:
    xs = newton(f, fp, 2, 4)
    assert xs[1] == x1 and 0 <= x1 - l <= adm
    d = C * adm
    for _ in range(3):
        d = d * d
    assert d / C == rep
    assert all(0 <= (xs[k + 1] - l) <= C * (xs[k] - l) ** 2 for k in range(4))
    assert mpmath.mpf(xs[4].numerator) / xs[4].denominator - l <= rep
assert F(5, 4) ** 3 < 2 and F(17, 10) ** 3 < 5 and F(173, 100) ** 2 < 3
# Vraies erreurs citées : 8,1e-7 ; 1,4e-13 ; 1,7e-18.
for f, fp, l, v in [(lambda t: t**3 - 2, lambda t: 3 * t * t, mpmath.cbrt(2), 8.1e-7),
                    (lambda t: t**3 - 5, lambda t: 3 * t * t, mpmath.cbrt(5), 1.4e-13),
                    (lambda t: t**2 - 3, lambda t: 2 * t, mpmath.sqrt(3), 1.7e-18)]:
    x4 = newton(f, fp, 2, 4)[4]
    assert abs((mpmath.mpf(x4.numerator) / x4.denominator - l) / v - 1) < 0.05

# num-l2-equations-newton-ex4 (+ v1, v2) : cycles.
assert newton(lambda t: t**3 - 5 * t, lambda t: 3 * t * t - 5, 1, 4) == [1, -1, 1, -1, 1]
assert newton(lambda t: t**3 - 20 * t, lambda t: 3 * t * t - 20, 2, 4) == [2, -2, 2, -2, 2]
r = math.sqrt(2 / 3)
assert abs((r**3 - 2 * r + 2) - (2 - 4 * math.sqrt(6) / 9)) < 1e-12 and 2 - 4 * math.sqrt(6) / 9 > 0.9
assert (-2) ** 3 - 2 * (-2) + 2 == -2 and (-1) ** 3 - 2 * (-1) + 2 == 3
# h(x) = -f(-x) = x^3 - 2x - 2 : Newton de f issu de -2 converge vers l (croissante).
xs = newton(lambda t: t**3 - 2 * t + 2, lambda t: 3 * t * t - 2, -2, 8, exact=False)
assert all(xs[k + 1] >= xs[k] for k in range(7)) and abs(xs[-1] - ell_c) < 1e-20
assert 1 - 2 - 2 == -3 and 8 - 4 - 2 == 2
xs = newton(lambda t: t**3 - 5 * t, lambda t: 3 * t * t - 5, 3, 8, exact=False)
assert abs(xs[-1] - mpmath.sqrt(5)) < 1e-20 and 8 - 10 == -2 and 27 - 15 == 12
xs = newton(lambda t: t**3 - 20 * t, lambda t: 3 * t * t - 20, 5, 8, exact=False)
assert abs(xs[-1] - 2 * mpmath.sqrt(5)) < 1e-20 and 64 - 80 == -16 and 125 - 100 == 25

# num-l2-equations-newton-ex5 (+ v1, v2) : (x_n - √a)/(x_n + √a) = q^(2^n) ; plus petit n.
y = symbols("y", positive=True)
A = symbols("A", positive=True)
xn1 = (y + A / y) / 2
assert simplify((xn1 - sqrt(A)) / (xn1 + sqrt(A)) - ((y - sqrt(A)) / (y + sqrt(A))) ** 2) == 0
assert simplify((2 - sqrt(2)) / (2 + sqrt(2)) - (3 - 2 * sqrt(2))) == 0
assert simplify((2 - sqrt(3)) / (2 + sqrt(3)) - (7 - 4 * sqrt(3))) == 0
assert simplify((2 - sqrt(5)) / (2 + sqrt(5)) + (9 - 4 * sqrt(5))) == 0
for a, prec, n, lq, lc in [(2, 100, 8, -0.7656, 0.4515), (3, 100, 7, -1.1439, 0.5396), (5, 1000, 10, -1.2539, 0.6505)]:
    s_ = mpmath.sqrt(a)
    q = (2 - s_) / (2 + s_)
    assert abs(mpmath.log10(abs(q)) - lq) < 1e-4 and abs(mpmath.log10(2 * s_) - lc) < 1e-4
    mpmath.mp.dps = 2200
    s_ = mpmath.sqrt(a)
    t = mpmath.mpf(2)
    errs = []
    for k in range(n + 1):
        errs.append(t - s_)
        t = (t + a / t) / 2
    mpmath.mp.dps = 50
    assert errs[n] <= mpmath.mpf(10) ** (-prec) < errs[n - 1]
    q = (2 - mpmath.sqrt(a)) / (2 + mpmath.sqrt(a))
    assert q**2 < (F(1, 30) if a == 2 else F(1, 10))

# ---------------------------------------------------------------------------------------------------------------
# num-l2-equations-ordre (carte)
# ---------------------------------------------------------------------------------------------------------------
gh = (x + 2 / x) / 2
assert simplify(diff(gh, x).subs(x, sqrt(2))) == 0 and simplify(diff(gh, x, 2).subs(x, sqrt(2)) - 1 / sqrt(2)) == 0
assert abs(float(1 / (2 * sqrt(2))) - 0.354) < 1e-3
assert abs(mpmath.sin(ell_cos) - 0.674) < 1e-3 and abs(-mpmath.log10(mpmath.sin(ell_cos)) - 0.17) < 0.005
# x/2 + x^2 sur [0, 1/4] : stable, e_{n+1}/e_n^2 = 1/(2 e_n) + 1.
assert F(1, 4) / 2 + F(1, 16) == F(3, 16) <= F(1, 4)
u = F(1, 4)
for _ in range(10):
    v = u / 2 + u * u
    assert v / u**2 == 1 / (2 * u) + 1 and 0 < v < u
    u = v
# 1/n! : e_{n+1}/e_n^p = (n!)^(p-1)/(n+1) -> +oo (ici p = 1,1, n = 60).
assert math.factorial(60) ** 0.1 / 61 > 1e6
# sin x_n / x_n -> 1 ; x^3 : e_{n+1}/e_n^2 = e_n.
assert abs(mpmath.sin(mpmath.mpf("1e-5")) / mpmath.mpf("1e-5") - 1) < 1e-9
# Affirmation : 2^-2^n et q = 3/2 : exposant -2^(n+1) + (3/2) 2^n = -2^(n-1).
assert all(-(2 ** (n + 1)) + F(3, 2) * 2**n == -(2 ** (n - 1)) for n in range(1, 20))
# Newton : e_{n+1}/e_n^2 -> f''(l)/(2 f'(l)) (formule générale vérifiée sur x^2 - 2).
assert simplify(limit(((x - (x * x - 2) / (2 * x)) - sqrt(2)) / (x - sqrt(2)) ** 2, x, sqrt(2)) - sqrt(2) / 4) == 0

# num-l2-equations-ordre-ex1 (+ v1, v2)
for k in range(1, 6):
    assert F(1, 2 ** (2 ** (k + 1))) / F(1, 2 ** (2**k)) ** 2 == 1
    assert F(1, 3 ** (3 ** (k + 1))) / F(1, 3 ** (3**k)) ** 3 == 1
    assert F(2, 10 ** (2 ** (k + 1))) / F(2, 10 ** (2**k)) ** 2 == F(1, 2)

# num-l2-equations-ordre-ex2 (+ v1, v2) : limite (u_n + a)/den -> 2a/den.
for num, den, a, b, lim in [(2, 3, 1, 2, F(2, 3)), (4, 5, 1, 4, F(2, 5)), (6, 5, 2, 3, F(4, 5))]:
    g = lambda t: (t * t + num) / F(den)
    assert g(F(a)) == a and g(F(b)) == b and F(2 * a, den) == lim
    u = F(0)
    for _ in range(40):
        v = g(u)
        assert u < v < a
        u = v
    w = mpmath.mpf(0)
    for _ in range(200):
        w = (w * w + num) / den
    w1 = (w * w + num) / den
    assert abs((w1 - a) / (w - a) - lim) < 1e-6

# num-l2-equations-ordre-ex3 (+ v1, v2)
assert min(n for n in range(50) if 2**n >= 100) == 7 and min(n for n in range(500) if 2**n >= 10**99) == 329
assert min(n for n in range(50) if 2**n >= 1000) == 10 and min(n for n in range(5000) if 2**n >= 10**999) == 3319
assert min(n for n in range(50) if 3**n >= 100) == 5 and min(n for n in range(500) if 1 + n >= 100) == 99
assert abs(99 / math.log10(2) - 328.87) < 0.01 and abs(999 / math.log10(2) - 3318.6) < 0.05

# num-l2-equations-ordre-ex4 (+ v1, v2) : Newton en un zéro de multiplicité m ; méthode modifiée.
for f, r, m, gm1, h, lim, h2 in [
    ((x - 1) ** 2 * (x + 1), 1, 2, (x + 1) / (3 * x + 1), (x**2 + x + 2) / (3 * x + 1), Rational(1, 2), Rational(1, 2)),
    ((x - 1) ** 3 * (x + 2), 1, 3, (x + 2) / (4 * x + 5), (x**2 + 2 * x + 6) / (4 * x + 5), Rational(2, 3), Rational(2, 9)),
    ((x - 2) ** 4 * (x + 1), 2, 4, (x + 1) / (5 * x + 2), (x**2 + 6 * x + 8) / (5 * x + 2), Rational(3, 4), Rational(1, 6)),
]:
    fp = diff(f, x)
    g = x - f / fp
    assert simplify(g - r - (x - r) * (1 - gm1)) == 0
    assert limit((g - r) / (x - r), x, r) == lim == 1 - Rational(1, m)
    assert simplify(x - m * f / fp - h) == 0
    assert simplify(diff(h, x).subs(x, r)) == 0 and simplify(diff(h, x, 2).subs(x, r)) == h2 != 0

# ---------------------------------------------------------------------------------------------------------------
# num-l2-equations-erreur (carte)
# ---------------------------------------------------------------------------------------------------------------
ell_xc = ell_cos
assert abs((mpmath.mpf("0.74") - mpmath.cos(mpmath.mpf("0.74"))) - 0.00153) < 1e-5 and 0.00153 <= 0.0016
assert abs(mpmath.mpf("0.74") - ell_xc - 0.00091) < 1e-5
assert mpmath.mpf("0.739") - mpmath.cos(mpmath.mpf("0.739")) < 0 < mpmath.mpf("0.740") - mpmath.cos(mpmath.mpf("0.740"))
# Contre-exemples : arctan 10 ~ 1,47 < 10 ; /(1/2) ~ 2,94 ; f' >= 1/2 sur [-1,1] ; x/100 ; x(x-5).
assert abs(math.atan(10) - 1.47) < 0.01 and abs(2 * math.atan(10) - 2.94) < 0.01 and 1 / (1 + 1) == 0.5
assert (-1) * (-1 - 5) * (1 * (1 - 5)) < 0
# Leurres du trou 1 et 2 : x^2 - 2x sur [0,2] (f'' = 2) ; 2x avec m = 2 ; m |f| avec f = x/100.
assert 2**2 - 2 * 2 == 0 and F(2, 4) < 1 and F(1, 100) * F(1, 100) < 1
# Trou 3, leurres : g(x) = x/2, u_0 = 1 : |u_1| = 1/2 > k |u_1 - u_0| = 1/4 ; g = 0,9 x : |u_1| = 0,9 > 0,1 ;
# n = 2 : |u_2| = 1/4 > (k^2/(1-k)) |u_2 - u_1| = 1/8.
assert F(1, 2) > F(1, 4) and F(9, 10) > F(1, 10) and F(1, 4) > (F(1, 4) / F(1, 2)) * F(1, 4)

# num-l2-equations-erreur-ex1 (+ v1, v2)
for f, xt, m, rep, leurres, ell in [
    (lambda t: t**3 + t - 1, F(68, 100), 1, F(5568, 10**6), [F(1392, 10**6), F(2784, 10**6), F(11136, 10**6)], 0.6823278038),
    (lambda t: t**3 + 2 * t - 1, F(45, 100), 2, F(44375, 10**7), [F(8875, 10**6), F(1775, 10**6), F(221875, 10**8)], 0.4533976515),
    (lambda t: t**5 + t - 1, F(3, 4), 1, F(13, 1024), [F(13, 6144), F(13, 2048), F(13, 512)], 0.7548776662),
]:
    assert abs(f(xt)) / m == rep and rep not in leurres and abs(float(xt) - ell) <= rep
    assert abs(float(xt) - ell) > min(leurres)  # le majorant obtenu en divisant par f'(1) est faux
assert F(68, 100) ** 3 == F(314432, 10**6) and F(45, 100) ** 3 == F(91125, 10**6)
assert 3 * F(6, 10) ** 2 + 1 == F(208, 100)

# num-l2-equations-erreur-ex2 (+ v1, v2)
for f, pts, rep in [
    (lambda t: t**3 + t - 1, [F(68, 100), F(69, 100), F(682, 1000), F(683, 1000)], (F(682, 1000), F(683, 1000))),
    (lambda t: t**3 + 2 * t - 1, [F(45, 100), F(46, 100), F(453, 1000), F(454, 1000)], (F(453, 1000), F(454, 1000))),
    (lambda t: t**5 + t - 1, [F(75, 100), F(76, 100), F(754, 1000), F(755, 1000)], (F(754, 1000), F(755, 1000))),
]:
    assert f(pts[0]) < 0 < f(pts[1]) and f(pts[2]) < 0 < f(pts[3]) and f(F(0)) == -1 and f(F(1)) > 0
assert abs(float((F(69, 100)) ** 3 + F(69, 100) - 1) - 0.0185) < 1e-4
assert abs(float(F(682, 1000) ** 3 + F(682, 1000) - 1) + 0.00079) < 1e-5 and abs(float(F(683, 1000) ** 3 + F(683, 1000) - 1) - 0.00161) < 1e-5
assert abs(float(F(46, 100) ** 3 + 2 * F(46, 100) - 1) - 0.0173) < 1e-4
assert abs(float(F(453, 1000) ** 3 + 2 * F(453, 1000) - 1) + 0.00104) < 1e-5 and abs(float(F(454, 1000) ** 3 + 2 * F(454, 1000) - 1) - 0.00158) < 1e-5
assert abs(float(F(76, 100) ** 5 + F(76, 100) - 1) - 0.0136) < 1e-4
assert abs(float(F(754, 1000) ** 5 + F(754, 1000) - 1) + 0.0023) < 1e-4 and abs(float(F(755, 1000) ** 5 + F(755, 1000) - 1) - 0.00032) < 1e-5

# num-l2-equations-erreur-ex3 (+ v1, v2) : estimation a posteriori k/(1-k) |u_2 - u_1|.
for c3, k, g12, rep, leurres, apriori in [
    (3, F(1, 4), F(3, 8), F(1, 243), [F(1, 81), F(1, 324), F(4, 243)], F(1, 32)),
    (4, F(3, 16), F(9, 32), F(3, 3328), [F(1, 256), F(3, 4096), F(1, 208)], F(9, 512)),
    (5, F(3, 20), F(9, 40), F(3, 10625), [F(1, 625), F(3, 12500), F(4, 2125)], F(9, 800)),
]:
    g = lambda t: (t**3 + 1) / F(c3)
    u1, u2 = g(F(0)), g(g(F(0)))
    assert g(F(1, 2)) == g12 <= F(1, 2) and 3 * F(1, 4) / c3 == k  # max de |g'| = 3x^2/c3 en 1/2
    d = u2 - u1
    assert k / (1 - k) * d == rep and rep not in leurres
    assert leurres == [d, k * d, d / (1 - k)] and k**2 / 2 == apriori
    ell = mpmath.findroot(lambda t: t**3 - c3 * t + 1, 0.3)
    assert 0 < ell < 0.5 and abs(mpmath.mpf(u2.numerator) / u2.denominator - ell) <= rep
assert F(28, 81) == (F(1, 27) + 1) / 3 and F(65, 256) == (F(1, 64) + 1) / 4 and F(126, 625) == (F(1, 125) + 1) / 5

# num-l2-equations-erreur-ex4 (+ v1, v2) : résidu de x_3 pour Héron.
for a, x3, m, rep, leurres in [
    (2, F(577, 408), 2, F(1, 332928), [F(1, 166464), F(1, 665856), F(1, 408)]),
    (3, F(18817, 10864), 2, F(1, 236052992), [F(1, 118026496), F(1, 472105984), F(1, 10864)]),
    (5, F(51841, 23184), 4, F(1, 2149991424), [F(1, 537497856), F(1, 1074995712), F(1, 3224987136)]),
]:
    assert heron(a, 3)[3] == x3 and (x3 * x3 - a) / m == rep
    err = mpmath.mpf(x3.numerator) / x3.denominator - mpmath.sqrt(a)
    assert 0 < err <= rep
    assert err > (leurres[2] if a == 5 else leurres[1])  # diviser par f'(2) (ou f'(3)) donne un majorant faux
assert 577**2 == 332929 and 408**2 == 166464 and 18817**2 == 354079489 and 10864**2 == 118026496
assert 51841**2 == 2687489281 and 23184**2 == 537497856 and 2 <= 51841 / 23184 <= 3
assert mpmath.mpf(577) / 408 - mpmath.sqrt(2) < F(1, 166464 * 28 // 10) < F(215, 10**8)  # minorant 2,8 sur [1,4 ; 1,5]

# num-l2-equations-erreur-ex5 (+ v1, v2) : test |f(x_n)| <= 1e-6 pour x^p, x_n = ((p-1)/p)^n.
for p, q, n, err, resid, prev in [(3, F(2, 3), 12, 7.7e-3, 4.6e-7, 1.2e-2), (2, F(1, 2), 10, 9.8e-4, 9.5e-7, 2.0e-3),
                                   (5, F(4, 5), 13, 5.5e-2, 5.0e-7, 6.9e-2)]:
    assert newton(lambda t: t**p, lambda t: p * t ** (p - 1), 1, 3) == [1, q, q**2, q**3]
    assert min(k for k in range(100) if q ** (p * k) <= F(1, 10**6)) == n
    assert abs(float(q**n) / err - 1) < 0.02 and abs(float(q ** (p * n)) / resid - 1) < 0.02 and abs(float(q ** (n - 1)) / prev - 1) < 0.03
assert float(F(2, 3) ** 12) * 1e6 > 7000 and float(F(1, 2) ** 10) * 1e6 > 970 and float(F(4, 5) ** 13) * 1e6 > 50000
assert abs(2 / math.log10(1.5) - 11.36) < 0.01 and abs(3 / math.log10(2) - 9.97) < 0.01 and abs(1.2 / math.log10(1.25) - 12.38) < 0.01

# ---------------------------------------------------------------------------------------------------------------
# QCM num-l2-equations-q1 … q8
# ---------------------------------------------------------------------------------------------------------------
assert min(n for n in range(30) if 2**n >= 1000) == 10  # q1
assert F(1, 4) == F(1, 2) - F(1, 4)  # q2 : c_1 = 1/4 à distance 1/4 de l = 1/2
phi = (1 + math.sqrt(5)) / 2  # q3
for g in [lambda t: 1 + 1 / t, lambda t: t * t - 1, lambda t: 1 / (t - 1), lambda t: -t * t + 2 * t + 1]:
    assert abs(g(phi) - phi) < 1e-12
assert 1 + 1 / F(3, 2) == F(5, 3) and 1 + F(1, 2) == F(3, 2) and F(4, 9) < 1
assert F(3, 2) ** 2 - 1 == F(5, 4) and 1 / (F(2) - 1) == 1 and -4 + 4 + 1 == 1
assert newton(lambda t: t**3 - 8, lambda t: 3 * t * t, 2, 2) == [2, 2, 2]  # q4
assert (F(2) + F(8, 2)) / 2 != 2 and (2 * F(2) + F(8, 2)) / 3 != 2 and (F(2) + 2 * F(8, 4)) / 3 != 2
X_ = symbols("X", positive=True)
assert simplify(X_ - (X_**3 - a_) / (3 * X_**2) - (2 * X_ + a_ / X_**2) / 3) == 0
assert newton(lambda t: t * t - 2, lambda t: 2 * t, 1, 2)[1:] == [F(3, 2), F(17, 12)]  # q5
assert simplify(Rational(2, 1) / (2 * 2 * sqrt(2)) - sqrt(2) / 4) == 0  # q6
assert 0 < mpmath.sin(ell_cos) < 1  # q7
assert F(1, 10**4) / 4 == F(25, 10**6) and F(1, 10**4) / 16 == F(625, 10**8)  # q8

print("num_l2_equations : toutes les vérifications passent.")
