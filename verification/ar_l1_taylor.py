"""Vérifications calculatoires du chapitre « ar-l1-taylor » (formules de Taylor et développements limités) : cartes et
exercices de content/*/analyse-reelle/l1-ar-l1-taylor.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice concerné (v1, v2 = variantes 1 et 2).
Lancer : python3 verification/ar_l1_taylor.py (doit finir sans erreur).
"""

import math

from sympy import (E, N, Rational, asin, asinh, atan, binomial, cos, cosh, diff, exp, factorial, integrate, limit, log,
                   oo, pi, series, simplify, sin, sinh, solve, sqrt, symbols, tan)

x, h, t, a = symbols("x h t a", real=True)
n = symbols("n", positive=True)
R = Rational


def egal(u, v):
    return simplify(u - v) == 0


def dl(f, ordre, var=x, point=0):
    """Partie régulière du développement limité de f à l'ordre donné."""
    return series(f, var, point, ordre + 1).removeO()


def coeffs(P, ordre, var=x):
    return [P.expand().coeff(var, k) for k in range(ordre + 1)]


# ===========================================================================
# Carte ar-l1-taylor-young
# ===========================================================================
# Exemple : sqrt(1+x) = 1 + x/2 - x^2/8 + o(x^2), limite -1/8.
assert egal(dl(sqrt(1 + x), 2), 1 + x / 2 - x**2 / 8)
assert limit((sqrt(1 + x) - 1 - x / 2) / x**2, x, 0) == R(-1, 8)
# Contre-exemple : x|x| / x^2 vaut ±1 (pas de DL d'ordre 2) ; f' = 2|x|.
assert limit(x * abs(x) / x**2, x, 0, "+") == 1 and limit(x * abs(x) / x**2, x, 0, "-") == -1
# Contre-exemple et affirmation 1 (x^3 sin(1/x)) : f' et taux de f' le long des deux suites.
k = symbols("k", integer=True, positive=True)
g = x**3 * sin(1 / x)
gp = diff(g, x)
assert egal(gp, 3 * x**2 * sin(1 / x) - x * cos(1 / x))
taux = gp / x
assert simplify(taux.subs(x, 1 / (2 * k * pi))) == -1 and simplify(taux.subs(x, 1 / ((2 * k + 1) * pi))) == 1
assert limit(g / x**2, x, 0) == 0 and limit(x**2 * sin(1 / x), x, 0) == 0
# Affirmation 3 : e^{-1/x^2} -> dérivées nulles en 0 (premières dérivées : limites nulles).
fe = exp(-1 / x**2)
for j in range(4):
    assert limit(diff(fe, x, j), x, 0) == 0

# ar-l1-taylor-young-ex1 et variantes : DL à partir des dérivées ; leurres différents de la bonne réponse.
X1 = x - 1
assert egal(2 - X1 + R(4, 2) * X1**2, 2 - X1 + 2 * X1**2)                       # ex1
assert not egal(2 - X1 + 4 * X1**2, 2 - X1 + 2 * X1**2)                         # ex1 leurre (oubli du 1/2)
assert not egal(2 - x + 2 * x**2, 2 - X1 + 2 * X1**2)                           # ex1 leurre (DL en 0)
assert egal(-1 + 3 * x + R(-2, 2) * x**2, -1 + 3 * x - x**2)                    # ex1 v1
X2 = x + 2
assert egal(0 + X2 + 0 * X2**2 + R(12, 6) * X2**3, X2 + 2 * X2**3)              # ex1 v2 (12/3! = 2)

# ar-l1-taylor-young-ex2 et variantes : différences finies, exactes sur les polynômes de degré 2 et testées sur exp.
c0, c1, c2 = symbols("c0 c1 c2")
P2 = lambda u: c0 + c1 * (u - a) + c2 / 2 * (u - a)**2  # f(a)=c0, f'(a)=c1, f''(a)=c2
for expr, val in [(lambda f: f(a + 2 * h) - 2 * f(a + h) + f(a), 1),           # ex2 : f''(a)
                  (lambda f: f(a + 3 * h) - 3 * f(a + h) + 2 * f(a), 3),       # ex2 v1 : 3 f''(a)
                  (lambda f: f(a + 2 * h) + f(a - 2 * h) - 2 * f(a), 4)]:      # ex2 v2 : 4 f''(a)
    assert egal(expr(P2) / h**2, val * c2)
    assert limit(expr(lambda u: exp(u)).subs(a, 0) / h**2, h, 0) == val       # f = exp, f''(0) = 1
    for leurre in {1: [2, R(1, 2), 0], 3: [R(3, 2), R(9, 2), 6], 4: [2, 1, 8]}[val]:
        assert leurre != val

# ar-l1-taylor-young-ex3 : x^3, x^4, -x^4 ont f'(0) = f''(0) = 0.
for f in (x**3, x**4, -x**4):
    assert diff(f, x).subs(x, 0) == 0 and diff(f, x, 2).subs(x, 0) == 0

# ar-l1-taylor-young-ex4 : mêmes calculs que le contre-exemple ; |f'(x)| <= 3x^2 + |x|.
assert limit(gp, x, 0) == 0

# ar-l1-taylor-young-ex5 : theta_h -> 1/2 pour exp ; theta = 1/sqrt(3) pour x^3 ; leurres faux.
theta_exp = log((exp(h) - 1) / h) / h
assert limit(theta_exp, h, 0) == R(1, 2)
assert egal(dl(theta_exp, 1, h), R(1, 2) + h / 24)
th = symbols("theta", positive=True)
assert solve(h**3 - h * 3 * (th * h)**2, th) == [1 / sqrt(3)]
assert 1 / sqrt(3) != R(1, 2)
# Cas f(x) = x^2 + x^3 (f''(0) = 2 ≠ 0) : la limite de theta_h vaut bien 1/2 (contrôle numérique).
f5 = lambda u: u**2 + u**3
for hv in [R(1, 1000), R(-1, 1000)]:
    sols = [s for s in solve(f5(hv) - f5(0) - hv * diff(f5(x), x).subs(x, th * hv), th) if s.is_real and 0 < s < 1]
    assert sols and all(abs(N(s) - 0.5) < 0.01 for s in sols)

# ===========================================================================
# Carte ar-l1-taylor-lagrange
# ===========================================================================
# Exemple : |sin x - x| <= |x|^3/6 (contrôle sur une grille).
for xv in [R(j, 4) for j in range(-40, 41)]:
    assert abs(N(sin(xv) - xv)) <= N(abs(xv)**3 / 6) + 1e-15
# Contre-exemples : exp avec M = e sur [-1,1] et x = 3 ; sin(10x) avec |f'| <= 10.
assert N(exp(3) - 4) > N(E * 9 / 2) and abs(N(exp(3) - 4) - 16.0855) < 1e-3 and abs(N(E * 9 / 2) - 12.2323) < 1e-3
assert N(1 - sin(1)) > 0.05 and abs(N(1 - sin(1)) - 0.159) < 1e-3 and N(100 * R(1, 100) / 2) == 0.5
assert N(1 - sin(1)) <= 0.5
# Affirmation 1 (fausse) : x = 1, e - 2 > 1/2. Affirmation 4 : |x|^{n+1}/(n+1)! -> 0.
assert N(E - 2) > 0.5
assert limit(5**n / factorial(n), n, oo) == 0

# ar-l1-taylor-lagrange-ex1 et variantes : polynômes de Taylor, majorants, erreur réelle <= majorant.
assert egal(dl(cos(x), 3), 1 - x**2 / 2)                                       # ex1
borne = R(1, 10)**4 / 24
assert abs(N(cos(R(1, 10)) - R(995, 1000))) <= N(borne)
assert abs(N(cos(R(1, 10)) - R(995, 1000)) - 4.165e-6) < 1e-8 and abs(N(borne) - 4.1667e-6) < 1e-9
assert N(R(1, 10)**3 / 6) > 1.6e-4                                              # ex1 : majorant d'ordre 2
assert len({R(1, 10)**4 / 24, R(1, 10)**3 / 6, R(1, 10)**4 / 6, R(1, 10)**2 / 2}) == 4
assert egal(dl(sin(x), 4), x - x**3 / 6)                                        # ex1 v1
b1 = R(2, 10)**5 / 120
assert abs(N(sin(R(2, 10)) - R(2, 10) + R(2, 10)**3 / 6)) <= N(b1) and abs(N(b1) - 2.667e-6) < 1e-8
assert len({b1, R(2, 10)**4 / 24, R(2, 10)**5 / 24, R(2, 10)**3 / 6}) == 4
assert egal(dl(cos(x), 5), 1 - x**2 / 2 + x**4 / 24)                            # ex1 v2
b2 = R(1, 2)**6 / 720
assert abs(N(cos(R(1, 2)) - (1 - R(1, 8) + R(1, 16) / 24))) <= N(b2) and abs(N(b2) - 2.17e-5) < 1e-7
assert len({b2, R(1, 2)**5 / 120, R(1, 2)**6 / 120, R(1, 2)**4 / 24}) == 4
assert egal(diff(cos(x), x, 4), cos(x)) and egal(diff(sin(x), x, 5), cos(x)) and egal(diff(cos(x), x, 6), -cos(x))


def plus_petit_n(condition, depart=0):
    m = depart
    while not condition(m):
        m += 1
    return m


# ar-l1-taylor-lagrange-ex2 et variantes : plus petit n.
assert plus_petit_n(lambda m: R(3) / factorial(m + 1) <= R(1, 10**3)) == 6      # ex2 (leurres 5, 7, 4)
assert plus_petit_n(lambda m: R(3) / factorial(m + 1) <= R(1, 10**5)) == 8      # ex2 v1 (leurres 7, 9, 5)
assert plus_petit_n(lambda m: R(1) / factorial(m + 1) <= R(1, 10**4)) == 7      # ex2 v2 (leurres 6, 8, 4)
S6 = sum(R(1, factorial(j)) for j in range(7))
assert S6 == R(1957, 720) and abs(N(E - S6) - 2.26e-4) < 1e-6 and N(E - S6) <= N(R(3, 5040))
assert N(E) < 3
# Erreurs réelles <= majorants (ex2, v1, v2).
assert N(E - sum(R(1, factorial(j)) for j in range(9))) <= N(R(3) / factorial(9))
assert abs(N(exp(-1) - sum(R((-1)**j, factorial(j)) for j in range(8)))) <= N(R(1) / factorial(8))

# ar-l1-taylor-lagrange-ex3 et variantes : dérivées de ln(1+x), majorants, plus petit n.
fl = log(1 + t)
for kk in range(1, 7):
    assert egal(diff(fl, t, kk), (-1)**(kk - 1) * factorial(kk - 1) / (1 + t)**kk)
assert plus_petit_n(lambda m: R(1, m + 1) <= R(1, 100), 1) == 99                 # ex3 (leurres 100, 10, 9)
assert plus_petit_n(lambda m: R(1, (m + 1) * 2**(m + 1)) <= R(1, 1000), 1) == 7  # ex3 v1 (leurres 6, 8, 999)
assert plus_petit_n(lambda m: R(1, (m + 1) * 5**(m + 1)) <= R(1, 10**4), 1) == 4  # ex3 v2 (leurres 3, 5, 9999)
assert 7 * 128 == 896 and 8 * 256 == 2048 and 4 * 625 == 2500 and 5 * 3125 == 15625
for xv, borne_fn in [(1, lambda m: R(1, m + 1)), (R(1, 2), lambda m: R(1, (m + 1) * 2**(m + 1))),
                     (R(1, 5), lambda m: R(1, (m + 1) * 5**(m + 1)))]:
    for m in range(1, 8):
        somme = sum(R((-1)**(j - 1), j) * xv**j for j in range(1, m + 1))
        assert abs(N(log(1 + xv) - somme)) <= N(borne_fn(m))

# ar-l1-taylor-lagrange-ex4 : minimum de phi(h) = 2 M0/h + M2 h/2 ; leurres faux ; constante de Landau.
M0, M2, hh = symbols("M0 M2 hh", positive=True)
phi = 2 * M0 / hh + M2 * hh / 2
h0 = solve(diff(phi, hh), hh)
assert h0 == [2 * sqrt(M0 / M2)] and egal(phi.subs(hh, h0[0]), 2 * sqrt(M0 * M2))
for leurre in [sqrt(M0 * M2), sqrt(2 * M0 * M2), 2 * M0 + M2 / 2]:
    assert not egal(leurre, 2 * sqrt(M0 * M2))
psi = M0 / hh + M2 * hh / 2
assert egal(psi.subs(hh, solve(diff(psi, hh), hh)[0]), sqrt(2 * M0 * M2))
# Contrôle : f = sin, M0 = M2 = 1, |cos| <= 2 (et même <= sqrt 2).
assert 1 <= 2 * sqrt(1 * 1) and 1 <= sqrt(2)

# ===========================================================================
# Carte ar-l1-taylor-reste-integral
# ===========================================================================
# Formule vérifiée symboliquement sur quelques fonctions et ordres.
for f in (exp(x), sin(x), log(1 + x), x**5):
    for m in range(0, 4):
        T_m = sum(diff(f, x, j).subs(x, 0) / factorial(j) * x**j for j in range(m + 1))
        reste = integrate((x - t)**m / factorial(m) * diff(f, x, m + 1).subs(x, t), (t, 0, x))
        assert egal(T_m + reste, f) or abs(N((T_m + reste - f).subs(x, R(1, 3)))) < 1e-12
# Leurres du texte à trous (trou 2), testés avec f = x^3, n = 1, a = 0 (vraie valeur x^3 - 0 = x^3).
f3 = t**3
vrai = integrate((x - t) * diff(f3, t, 2), (t, 0, x))
assert egal(vrai, x**3)
assert not egal(integrate((x - t)**2 / 2 * diff(f3, t, 2), (t, 0, x)), x**3)
assert not egal(integrate((x - t) * diff(f3, t, 1), (t, 0, x)), x**3)
assert not egal(x**2 / 2 * diff(f3, t, 2).subs(t, 0), x**3)
# Pour n = 0 aussi (f = x) : les deux leurres intégraux sont faux.
assert not egal(integrate((x - t) * 1, (t, 0, x)), x) and not egal(integrate(t, (t, 0, x)), x)
# Contre-exemple : x^2 sin(1/x^2), f'(x_k) = -2 sqrt(2 k pi) en x_k = 1/sqrt(2 k pi).
fq = x**2 * sin(1 / x**2)
fqp = diff(fq, x)
assert egal(fqp, 2 * x * sin(1 / x**2) - 2 / x * cos(1 / x**2))
assert egal(fqp.subs(x, 1 / sqrt(2 * k * pi)), -2 * sqrt(2 * k * pi))
# Affirmation 1 (fausse) : exp, n = 0, x = -1.
assert N(exp(-1)) < 1
# Affirmation 3 : intégrale de |x-t|^n/n! sur [a,x] = |x-a|^{n+1}/(n+1)!.
assert egal(integrate((x - t)**3 / factorial(3), (t, 0, x)), x**4 / factorial(4))

# ar-l1-taylor-reste-integral-ex1 et variantes : intégrales.
assert egal(integrate((1 - t)**2 * exp(t), (t, 0, 1)), 2 * E - 5)               # ex1
assert egal(2 * (E - R(5, 2)), 2 * E - 5)
for leurre in [E - R(5, 2), 2 * E - 4, 5 - 2 * E]:
    assert not egal(leurre, 2 * E - 5)
assert egal(integrate((1 - t)**3 * exp(t), (t, 0, 1)), 6 * E - 16)              # ex1 v1
for leurre in [E - R(8, 3), 6 * E - 15, 3 * E - 8]:
    assert not egal(leurre, 6 * E - 16)
assert egal(integrate((pi / 2 - t)**2 * cos(t), (t, 0, pi / 2)), pi - 2)        # ex1 v2
for leurre in [pi / 2 - 1, 2 - pi, pi + 2]:
    assert not egal(leurre, pi - 2)
assert egal(diff(sin(x), x, 3), -cos(x))

# ar-l1-taylor-reste-integral-ex2 et variantes : dérivées, signe du reste, inégalités et leurres (faux en x = 1).
fl = log(1 + x)
assert egal(diff(fl, x, 4), -6 / (1 + x)**4) and egal(diff(fl, x, 3), 2 / (1 + x)**3)
assert egal(dl(fl, 3), x - x**2 / 2 + x**3 / 3)
xs = [R(j, 10) for j in range(1, 60)]
T3 = x - x**2 / 2 + x**3 / 3
T2 = x - x**2 / 2
assert all(N(log(1 + v) - T3.subs(x, v)) < 0 for v in xs)                       # ex2 : ln(1+x) <= T3
assert N(log(2)) < N(T3.subs(x, 1)) and N(log(2)) > N(T2.subs(x, 1)) and N(log(2)) < 1   # ex2 leurres
assert all(N(log(1 + v) - T2.subs(x, v)) > 0 for v in xs)                       # ex2 v1 : ln(1+x) >= T2
assert N(log(2)) > 0.5 and N(log(2)) < 1 and N(log(2)) < N(R(5, 6))             # ex2 v1 leurres
E2 = 1 - x + x**2 / 2
assert egal(dl(exp(-x), 2), E2) and egal(diff(exp(-x), x, 3), -exp(-x))
assert all(N(exp(-v) - E2.subs(x, v)) < 0 for v in xs)                          # ex2 v2 : e^{-x} <= T2
assert N(exp(-1)) < 0.5 and N(exp(-1)) > 0 and N(exp(-1)) > N(R(1, 3))          # ex2 v2 leurres
assert all(N(exp(-v) - (E2 - x**3 / 6).subs(x, v)) > 0 for v in xs)             # ex2 v2 : e^{-x} >= T3

# ar-l1-taylor-reste-integral-ex3 : signe de e^x - T_n(x) selon la parité de n et le signe de x.
for m in range(0, 7):
    Tm = sum(x**j / factorial(j) for j in range(m + 1))
    for v in [R(-3), R(-1), R(-1, 10), R(1, 10), R(2)]:
        d = N(exp(v) - Tm.subs(x, v))
        if v > 0:
            assert d > 0
        else:
            assert (d > 0) == (m % 2 == 1) and d != 0
assert abs(N(exp(-3)) - 0.0498) < 1e-3 and 1 - 3 + R(9, 2) == R(5, 2)

# ar-l1-taylor-reste-integral-ex4 : intégrale de (1-t)^n/n! = 1/(n+1)!, encadrement de e - S_n, n! S_n entier.
for m in range(0, 8):
    assert egal(integrate((1 - t)**m / factorial(m), (t, 0, 1)), R(1, factorial(m + 1)))
    Sm = sum(R(1, factorial(j)) for j in range(m + 1))
    assert 0 < N(E - Sm) < N(E / factorial(m + 1)) < N(3 / factorial(m + 1))
    assert (factorial(m) * Sm).is_integer
    if m >= 2:
        assert R(3, m + 1) <= 1

# ===========================================================================
# Carte ar-l1-taylor-dl-usuels
# ===========================================================================
for m in range(0, 5):
    assert egal(dl(exp(x), m), sum(x**j / factorial(j) for j in range(m + 1)))
    assert egal(dl(cos(x), 2 * m + 1), sum((-1)**j * x**(2 * j) / factorial(2 * j) for j in range(m + 1)))
    assert egal(dl(sin(x), 2 * m + 2), sum((-1)**j * x**(2 * j + 1) / factorial(2 * j + 1) for j in range(m + 1)))
    assert egal(dl(cosh(x), 2 * m + 1), sum(x**(2 * j) / factorial(2 * j) for j in range(m + 1)))
    assert egal(dl(sinh(x), 2 * m + 2), sum(x**(2 * j + 1) / factorial(2 * j + 1) for j in range(m + 1)))
    assert egal(dl(1 / (1 - x), m), sum(x**j for j in range(m + 1)))
    assert egal(dl(atan(x), 2 * m + 2), sum((-1)**j * x**(2 * j + 1) / (2 * j + 1) for j in range(m + 1)))
    if m >= 1:
        assert egal(dl(log(1 + x), m), sum((-1)**(j - 1) * x**j / j for j in range(1, m + 1)))
    for al in [R(1, 2), R(-1, 2), R(-2), R(1, 3), pi]:
        assert egal(dl((1 + x)**al, m), sum(binomial(al, j) * x**j for j in range(m + 1)))
# Reste de 1/(1-x) et de arctan : identités exactes.
assert egal(1 / (1 - x) - sum(x**j for j in range(4)), x**4 / (1 - x))
assert egal(1 / (1 + t**2) - sum((-1)**j * t**(2 * j) for j in range(3)), (-1)**3 * t**6 / (1 + t**2))
# Trous : leurres faux (ln, (1+x)^alpha pour alpha générique, sin).
assert not egal(x + x**2 / 2 + x**3 / 3, dl(log(1 + x), 3)) and not egal(x - x**2 / 2 + x**3 / 6, dl(log(1 + x), 3))
assert not egal(1 + x - x**2 / 2 + x**3 / 3, dl(log(1 + x), 3))
al = R(3, 7)
vrai2 = dl((1 + x)**al, 2)
for L in [1 + al * x + al * (al + 1) / 2 * x**2, 1 + al * x + al**2 / 2 * x**2, 1 + al * x + al * (al - 1) * x**2]:
    assert not egal(L, vrai2)
assert egal(dl(sin(x), 4), x - x**3 / 6)
for L in [x - x**3 / 3, x + x**3 / 6, x - x**3 / 6 + x**4 / 24]:
    assert not egal(L, dl(sin(x), 4))
# Exemples : 1/(1+x), 1/sqrt(1+x), limite (cos x - 1 + x^2/2)/x^4 = 1/24.
assert egal(dl(1 / (1 + x), 3), 1 - x + x**2 - x**3)
assert egal(dl(1 / sqrt(1 + x), 2), 1 - x / 2 + R(3, 8) * x**2) and binomial(R(-1, 2), 2) == R(3, 8)
assert limit((cos(x) - 1 + x**2 / 2) / x**4, x, 0) == R(1, 24)
# Contre-exemple : ln(1+x) - x + x^2/2 ~ x^3/3 ; sqrt(x)/x -> +oo.
assert limit((log(1 + x) - x + x**2 / 2) / x**3, x, 0) == R(1, 3)
assert limit(sqrt(x) / x, x, 0, "+") == oo
# Affirmations : ln(1.1) > 0.095 ; (1+x)^2 ≠ 1 + 2x + 2x^2 ; coefficients x^3 de cos et x^4 de arctan nuls.
assert N(log(R(11, 10))) > 0.095 and abs(N(log(R(11, 10))) - 0.09531) < 1e-5
assert not egal((1 + x)**2, 1 + 2 * x + 2 * x**2)
assert dl(cos(x), 3).coeff(x, 3) == 0 and dl(atan(x), 4).coeff(x, 4) == 0

# ar-l1-taylor-dl-usuels-ex1 et variantes : limites ; leurres distincts.
assert limit((exp(x) - 1 - x) / x**2, x, 0) == R(1, 2)                          # ex1 (leurres 1, 0, 1/6)
assert limit((x - sin(x)) / x**3, x, 0) == R(1, 6)                              # ex1 v1 (leurres -1/6, 0, 1/3)
assert limit((log(1 + x) - x) / x**2, x, 0) == R(-1, 2)                         # ex1 v2 (leurres 1/2, -1, 0)

# ar-l1-taylor-dl-usuels-ex2 et variantes.
assert egal(dl(sqrt(1 + 2 * x), 3), 1 + x - x**2 / 2 + x**3 / 2)                 # ex2
for L in [1 + x - x**2 / 8 + x**3 / 16, 1 + x + x**2 / 2 + x**3 / 2, 1 + 2 * x - x**2 / 2 + x**3 / 2]:
    assert not egal(L, dl(sqrt(1 + 2 * x), 3))
assert [binomial(R(1, 2), j) for j in range(4)] == [1, R(1, 2), R(-1, 8), R(1, 16)]
assert egal(dl((1 - x)**-2, 3), 1 + 2 * x + 3 * x**2 + 4 * x**3)                # ex2 v1
for L in [1 - 2 * x + 3 * x**2 - 4 * x**3, 1 + 2 * x + 4 * x**2 + 8 * x**3, 1 + 2 * x + x**2]:
    assert not egal(L, dl((1 - x)**-2, 3))
assert all(binomial(-2, j) == (-1)**j * (j + 1) for j in range(6))
assert egal(dl(log(1 - 3 * x), 3), -3 * x - R(9, 2) * x**2 - 9 * x**3)          # ex2 v2
for L in [-3 * x + R(9, 2) * x**2 - 9 * x**3, -3 * x - R(3, 2) * x**2 - x**3, 3 * x - R(9, 2) * x**2 + 9 * x**3]:
    assert not egal(L, dl(log(1 - 3 * x), 3))

# ar-l1-taylor-dl-usuels-ex3 et variantes : DL en un point a ≠ 0.
v = dl(log(2 + h), 3, h)
assert egal(v, log(2) + h / 2 - h**2 / 8 + h**3 / 24)                           # ex3
for L in [log(2) + h / 2 - h**2 / 4 + h**3 / 8, h / 2 - h**2 / 8 + h**3 / 24, log(2) + h - h**2 / 2 + h**3 / 3]:
    assert not egal(L, v)
assert diff(log(x), x).subs(x, 2) == R(1, 2) and diff(log(x), x, 2).subs(x, 2) / 2 == R(-1, 8)
assert diff(log(x), x, 3).subs(x, 2) / 6 == R(1, 24)
v = dl(sqrt(4 + h), 3, h)
assert egal(v, 2 + h / 4 - h**2 / 64 + h**3 / 512)                               # ex3 v1
for L in [2 + h / 8 - h**2 / 128 + h**3 / 1024, 2 + h / 4 + h**2 / 64 + h**3 / 512, 2 + h / 4 - h**2 / 16 + h**3 / 64]:
    assert not egal(L, v)
v = dl(cos(pi / 3 + h), 3, h)
assert egal(v, R(1, 2) - sqrt(3) / 2 * h - h**2 / 4 + sqrt(3) / 12 * h**3)       # ex3 v2
for L in [R(1, 2) + sqrt(3) / 2 * h - h**2 / 4 - sqrt(3) / 12 * h**3,
          R(1, 2) - sqrt(3) / 2 * h - h**2 / 2 + sqrt(3) / 6 * h**3,
          sqrt(3) / 2 - h / 2 - sqrt(3) / 4 * h**2 + h**3 / 12]:
    assert not egal(L, v)

# ar-l1-taylor-dl-usuels-ex4 et variantes : binom(-1/2, k) = (-1)^k C(2k,k)/4^k ; coefficients.
for kk in range(1, 12):
    assert binomial(R(-1, 2), kk) == R((-1)**kk * binomial(2 * kk, kk), 4**kk)
    assert factorial(2 * kk) / (2**kk * factorial(kk)) == math.prod(range(1, 2 * kk, 2))
P = dl(1 / sqrt(1 - x), 5)
assert coeffs(P, 5) == [1, R(1, 2), R(3, 8), R(5, 16), R(35, 128), R(63, 256)]
assert P.coeff(x, 4) == R(35, 128) and R(35, 128) not in (R(35, 64), R(-35, 128), R(5, 16))   # ex4
assert P.coeff(x, 5) == R(63, 256) and R(63, 256) not in (R(63, 128), R(-63, 256), R(35, 128))  # ex4 v1
Q = dl(1 / sqrt(1 + 4 * x), 3)
assert coeffs(Q, 3) == [1, -2, 6, -20] and -20 not in (20, R(-5, 16), -80)      # ex4 v2

# ===========================================================================
# Carte ar-l1-taylor-dl-operations
# ===========================================================================
# Exemples : 1/cos x = 1 + x^2/2 + o(x^3), tan x = x + x^3/3 + o(x^3), arctan à l'ordre 6, 1/(1+x^2) exact.
assert egal(dl(1 / cos(x), 3), 1 + x**2 / 2) and egal(dl(tan(x), 3), x + x**3 / 3)
assert egal(dl(atan(x), 6), x - x**3 / 3 + x**5 / 5)
assert egal(1 / (1 + x**2) - (1 - x**2 + x**4), -x**6 / (1 + x**2))
u = symbols("u")
assert egal(1 / (1 - u) - (1 + u), u**2 / (1 - u))
# Contre-exemple « composée sans g -> 0 » : f(y) = |y - 1| vaut 1 - y près de 0 ; f(1 + x) = |x|.
assert all(abs(R(j, 10) - 1) == 1 - R(j, 10) for j in range(-5, 6))
assert limit(abs(x) / x, x, 0, "+") == 1 and limit(abs(x) / x, x, 0, "-") == -1
# Affirmation 1 (fausse) : e^x - (1+x) = x^2/2 + o(x^2).
assert limit((exp(x) - 1 - x) / x**2, x, 0) == R(1, 2)
# Trou 1 : leurres faux (exemple P = 1 + x, Q = 1 + x, n = 2 : [PQ]_2 = 1 + 2x + x^2).
Pp, Qq = 1 + x, 1 + x
assert egal((Pp * Qq).expand(), 1 + 2 * x + x**2)
for L in [1 + 2 * x, Pp + Qq, Pp * 1 + 1 * Qq]:
    assert not egal(L, 1 + 2 * x + x**2)
# Trou 2 : Q(0) = 1 donne un calcul faux (même exemple que |y - 1|).

# ar-l1-taylor-dl-operations-ex1 et variantes : produits.
assert egal(dl(exp(x) * cos(x), 4), 1 + x - x**3 / 3 - x**4 / 6)                # ex1
for L in [1 + x + x**3 / 3 - x**4 / 6, 1 + x - x**3 / 3 + x**4 / 24, 1 + x + x**2 - x**3 / 3 - x**4 / 6]:
    assert not egal(L, dl(exp(x) * cos(x), 4))
assert R(1, 24) - R(1, 4) + R(1, 24) == R(-1, 6)
assert egal(dl(exp(x) * sin(x), 4), x + x**2 + x**3 / 3)                        # ex1 v1
for L in [x + x**2 + x**3 / 3 + x**4 / 6, x + x**2 + 2 * x**3 / 3, x + x**2 - x**3 / 3]:
    assert not egal(L, dl(exp(x) * sin(x), 4))
assert egal(dl(sin(x) * log(1 + x), 4), x**2 - x**3 / 2 + x**4 / 6)             # ex1 v2
for L in [x**2 - x**3 / 2 + x**4 / 3, x**2 + x**3 / 2 + x**4 / 6, x**2 - x**3 / 2 - x**4 / 6]:
    assert not egal(L, dl(sin(x) * log(1 + x), 4))

# ar-l1-taylor-dl-operations-ex2 et variantes : quotients.
assert egal(dl(tan(x), 5), x + x**3 / 3 + 2 * x**5 / 15) and dl(tan(x), 6).coeff(x, 6) == 0     # ex2
for L in [x + x**3 / 3 + x**5 / 15, x - x**3 / 3 + 2 * x**5 / 15, x + x**3 / 6 + x**5 / 120]:
    assert not egal(L, dl(tan(x), 5))
assert egal(dl(1 / cos(x), 5), 1 + x**2 / 2 + R(5, 24) * x**4) and R(25 - 10 + 1, 120) == R(2, 15)
assert egal(dl(1 / cos(x), 4), 1 + x**2 / 2 + R(5, 24) * x**4)                  # ex2 v1
for L in [1 + x**2 / 2 - x**4 / 24, 1 + x**2 / 2 + x**4 / 4, 1 - x**2 / 2 + R(5, 24) * x**4]:
    assert not egal(L, dl(1 / cos(x), 4))
assert egal(dl(x / sin(x), 4), 1 + x**2 / 6 + R(7, 360) * x**4)                 # ex2 v2
for L in [1 + x**2 / 6 - x**4 / 120, 1 + x**2 / 6 + x**4 / 36, 1 - x**2 / 6 + R(7, 360) * x**4]:
    assert not egal(L, dl(x / sin(x), 4))
assert R(1, 36) - R(1, 120) == R(7, 360)

# ar-l1-taylor-dl-operations-ex3 et variantes : composées.
assert egal(dl(exp(sin(x)), 4), 1 + x + x**2 / 2 - x**4 / 8)                    # ex3
for L in [dl(exp(x), 4), 1 + x + x**2 / 2 - x**3 / 6, 1 + x + x**2 / 2 + x**4 / 8]:
    assert not egal(L, dl(exp(sin(x)), 4))
assert egal(dl((x - x**3 / 6)**2, 4), x**2 - x**4 / 3)
assert egal(dl(log(cos(x)), 4), -x**2 / 2 - x**4 / 12)                          # ex3 v1
for L in [-x**2 / 2 + x**4 / 24, -x**2 / 2 - x**4 / 8, x**2 / 2 - x**4 / 12]:
    assert not egal(L, dl(log(cos(x)), 4))
assert egal(dl(exp(cos(x)), 4), E - E / 2 * x**2 + E / 6 * x**4)                # ex3 v2
for L in [1 - x**2 / 2 + x**4 / 6, E - E / 2 * x**2 + E / 24 * x**4, E - E / 2 * x**2 + E / 8 * x**4]:
    assert not egal(L, dl(exp(cos(x)), 4))

# ar-l1-taylor-dl-operations-ex4 et variantes : primitives.
assert egal(dl(1 / sqrt(1 - x**2), 4), 1 + x**2 / 2 + R(3, 8) * x**4)
assert egal(dl(asin(x), 5), x + x**3 / 6 + R(3, 40) * x**5) and egal(diff(asin(x), x), 1 / sqrt(1 - x**2))   # ex4
for L in [x + x**3 / 2 + R(3, 8) * x**5, x - x**3 / 6 + R(3, 40) * x**5, x + x**3 / 6 + x**5 / 120]:
    assert not egal(L, dl(asin(x), 5))
F = integrate(exp(-t**2), (t, 0, x))
assert egal(dl(F, 5), x - x**3 / 3 + x**5 / 10)                                 # ex4 v1
for L in [x - x**3 / 3 + x**5 / 2, x - x**3 + x**5 / 2, x + x**3 / 3 + x**5 / 10]:
    assert not egal(L, dl(F, 5))
assert egal(diff(asinh(x), x), 1 / sqrt(1 + x**2)) and asinh(0) == 0
assert egal(dl(asinh(x), 5), x - x**3 / 6 + R(3, 40) * x**5)                     # ex4 v2
for L in [x + x**3 / 6 + R(3, 40) * x**5, x - x**3 / 2 + R(3, 8) * x**5, x - x**3 / 6 + x**5 / 120]:
    assert not egal(L, dl(asinh(x), 5))

# ar-l1-taylor-dl-operations-ex5 et variantes : suites.
y = symbols("y", positive=True)  # y = 1/n
assert egal(dl(log(1 + y) / y, 1, y), 1 - y / 2)                                # ex5 Q1
assert limit(n * (E - (1 + 1 / n)**n), n, oo) == E / 2                          # ex5 (leurres -e/2, 1/2, 0)
assert egal(dl(log(1 + 2 * y) / y, 1, y), 2 - 2 * y)                            # ex5 v1 Q1
assert limit(n * ((1 + 2 / n)**n - E**2), n, oo) == -2 * E**2                   # ex5 v1 (leurres 2e^2, -e^2, -2)
assert egal(dl(log(1 - y) / y, 1, y), -1 - y / 2)                               # ex5 v2 Q1
assert limit(n * ((1 - 1 / n)**n - exp(-1)), n, oo) == -exp(-1) / 2             # ex5 v2 (leurres 1/(2e), -1/2, -e/2)
for vrai_, leurres in [(E / 2, [-E / 2, R(1, 2), 0]), (-2 * E**2, [2 * E**2, -E**2, -2]),
                       (-exp(-1) / 2, [exp(-1) / 2, R(-1, 2), -E / 2])]:
    assert all(not egal(L, vrai_) for L in leurres)

# ===========================================================================
# Carte ar-l1-taylor-dl-position
# ===========================================================================
# Exemples : ln(1+x) sous sa tangente, sin traverse.
assert egal(dl(log(1 + x), 2), x - x**2 / 2) and egal(dl(sin(x), 3), x - x**3 / 6)
# Cours : ch x + cos x = 2 + x^4/12 + o(x^4), dérivée seconde nulle en 0.
assert egal(dl(cosh(x) + cos(x), 4), 2 + x**4 / 12) and diff(cosh(x) + cos(x), x, 2).subs(x, 0) == 0
# Contre-exemples : x^3, x^4, -x^4 ont un DL d'ordre 2 nul ; x^4'' (0) = 0.
for f in (x**3, x**4, -x**4):
    assert dl(f, 2) == 0
assert diff(x**4, x, 2).subs(x, 0) == 0

# ar-l1-taylor-dl-position-ex1 et variantes.
assert egal(dl(x / (1 + x**2), 3), x - x**3)                                    # ex1 : traverse, signe de -x^3
assert egal(x / (1 + x**2) - x, -x**3 / (1 + x**2))
assert egal(dl(cos(x) + x, 2), 1 + x - x**2 / 2)                                # ex1 v1 : en dessous
assert egal(dl(log(1 + x) + x**2, 2), x + x**2 / 2)                             # ex1 v2 : au-dessus

# ar-l1-taylor-dl-position-ex2 et variantes : f'(0) = f''(0) = 0 et DL.
for f, ordre, partie in [(2 * cos(x) + x**2, 4, 2 + x**4 / 12),                 # ex2 : minimum local strict
                         (exp(x) - x - x**2 / 2, 3, 1 + x**3 / 6),              # ex2 v1 : pas d'extremum
                         (log(1 + x**2) - x**2, 4, -x**4 / 2)]:                 # ex2 v2 : maximum local strict
    assert diff(f, x).subs(x, 0) == 0 and diff(f, x, 2).subs(x, 0) == 0
    assert egal(dl(f, ordre), partie)
assert egal(diff(log(1 + x**2), x, 2), (2 - 2 * x**2) / (1 + x**2)**2)

# ar-l1-taylor-dl-position-ex3 et variantes : valeur de a pour une inflexion.
for f, aval, c3 in [(log(1 + x) - x + a * x**2, R(1, 2), R(1, 3)),              # ex3 (leurres -1/2, 0)
                    (exp(x) - 1 - x + a * x**2, R(-1, 2), R(1, 6)),             # ex3 v1 (leurres 1/2, 0)
                    (sin(x) - x + a * x**2, 0, R(-1, 6))]:                      # ex3 v2 (leurres 1/6, -1/6)
    P3 = dl(f, 3)
    assert P3.coeff(x, 0) == 0 and P3.coeff(x, 1) == 0
    assert solve(P3.coeff(x, 2), a) == [aval] and P3.coeff(x, 3) == c3

# ar-l1-taylor-dl-position-ex4 et variantes : prolongement, DL, f'(x) -> f'(0), position.
fA = (exp(x) - 1) / x
assert limit(fA, x, 0) == 1 and egal(dl(fA, 2), 1 + x / 2 + x**2 / 6)          # ex4 : au-dessus de y = 1 + x/2
assert limit(diff(fA, x), x, 0) == R(1, 2)
assert egal(dl(x * exp(x) - exp(x) + 1, 3), x**2 / 2 + x**3 / 3)
fB = sin(x) / x
assert limit(fB, x, 0) == 1 and egal(dl(fB, 2), 1 - x**2 / 6)                  # ex4 v1 : en dessous de y = 1
assert limit(diff(fB, x), x, 0) == 0 and egal(dl(x * cos(x) - sin(x), 3), -x**3 / 3)
fC = (1 - cos(x)) / x
assert limit(fC, x, 0) == 0 and egal(dl(fC, 3), x / 2 - x**3 / 24)             # ex4 v2 : traverse y = x/2
assert dl(fC, 2).coeff(x, 2) == 0
assert limit(diff(fC, x), x, 0) == R(1, 2) and egal(dl(x * sin(x) - 1 + cos(x), 4), x**2 / 2 - x**4 / 8)

print("ar_l1_taylor : toutes les vérifications passent.")
