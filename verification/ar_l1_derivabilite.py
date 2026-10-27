"""Vérifications calculatoires du chapitre « ar-l1-derivabilite » (limite de la dérivée, dérivée d'une réciproque,
formule de Leibniz, dérivabilité et continuité) : cartes et exercices de content/*/analyse-reelle/l1-ar-l1-derivabilite.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice concerné (v1, v2 = variantes 1 et 2).
Lancer : python3 verification/ar_l1_derivabilite.py (doit finir sans erreur).
"""

from sympy import (Abs, E, Poly, Rational, S, acos, asin, asinh, atan, binomial, cos, cosh, diff, exp, factorial,
                   limit, log, oo, pi, simplify, sin, sinh, solve, sqrt, symbols)

x, y = symbols("x y", real=True)
h = symbols("h", positive=True)
n = symbols("n", integer=True, positive=True)


def egal(u, v):
    return simplify(u - v) == 0


# ===========================================================================
# Carte ar-l1-derivabilite-limite-derivee
# ===========================================================================
# Exemple : (x|x|)' = 2|x| sur R* ; arcsin' → +∞ en 1^-.
assert egal(diff(h * h, h), 2 * h) and egal(diff(-h * h, h), -2 * h)
assert limit(1 / sqrt(1 - x**2), x, 1, "-") == oo
# Contre-exemple : f'(x) = 2x sin(1/x) - cos(1/x) ; valeurs -1 en 1/(2nπ), 1 en 1/((2n+1)π).
fp = diff(x**2 * sin(1 / x), x)
assert egal(fp, 2 * x * sin(1 / x) - cos(1 / x))
assert simplify(fp.subs(x, 1 / (2 * n * pi))) == -1 and simplify(fp.subs(x, 1 / ((2 * n + 1) * pi))) == 1
assert limit(x * sin(1 / x), x, 0) == 0

# ar-l1-derivabilite-limite-derivee-ex1 et variantes : continuité en 0 et limite de f'.
for f, lim_fp in [(x**2 * log(x), 0),                 # ex1
                  (x * log(x), -oo),                  # ex1 v1
                  (sin(x) + x**2 * log(x), 1)]:       # ex1 v2
    assert limit(f, x, 0, "+") == 0
    assert limit(diff(f, x), x, 0, "+") == lim_fp
assert limit(x * log(x) / x, x, 0, "+") == -oo        # ex1 v1 : taux = ln x

# ar-l1-derivabilite-limite-derivee-ex2 et variantes : (a, b) = (2, 1), (0, 1), (1/2, 1) ; leurres faux.
for g, a_ok, b_ok, leurres in [(exp(2 * x), 2, 1, [(1, 2), (2, 0), (1, 1)]),                       # ex2
                               (cos(x), 0, 1, [(1, 1), (-1, 1), (0, 0)]),                          # ex2 v1
                               (sqrt(1 + x), Rational(1, 2), 1, [(1, 1), (2, 1), (Rational(1, 2), 0)])]:  # ex2 v2
    assert g.subs(x, 0) == b_ok and egal(diff(g, x).subs(x, 0), a_ok)
    assert all((a, b) != (a_ok, b_ok) and (b != g.subs(x, 0) or a != diff(g, x).subs(x, 0)) for a, b in leurres)
assert limit(diff(sqrt(1 + x), x), x, -1, "+") == oo   # ex2 v2 : pas dérivable en -1

# ar-l1-derivabilite-limite-derivee-ex3 : f', f'' et leurs limites en 0.
f = exp(-1 / x**2)
assert limit(f, x, 0) == 0
assert egal(diff(f, x), 2 / x**3 * exp(-1 / x**2)) and limit(diff(f, x), x, 0) == 0
assert egal(diff(f, x, 2), (4 / x**6 - 6 / x**4) * exp(-1 / x**2)) and limit(diff(f, x, 2), x, 0) == 0

# ar-l1-derivabilite-limite-derivee-ex4 et variantes : formules de la dérivée et valeur aux points de contrôle.
a = symbols("a", positive=True)
for k, seuil in [(1, 2), (2, 3), (3, 4)]:            # ex4 : sin(1/x), v1 : sin(1/x^2), v2 : sin(1/x^3)
    g = x**a * sin(1 / x**k)
    assert egal(diff(g, x), a * x**(a - 1) * sin(1 / x**k) - k * x**(a - k - 1) * cos(1 / x**k))
    xn = (2 * n * pi)**Rational(-1, k)
    # g'(x_n) = -k x_n^(a-k-1) : vaut -k pour a = k+1 (seuil), tend vers -∞ pour a < seuil.
    assert simplify(diff(g, x).subs(x, xn).subs(a, seuil)) == -k
    assert limit((-k * xn**(Rational(seuil * 2 - 1, 2) - k - 1)), n, oo) == -oo   # a = seuil - 1/2
    assert limit(xn, n, oo) == 0
# Taux en 0 pour a ≤ 1 : en y_n = (2nπ+π/2)^(-1/k), x^(a-1) sin(1/x^k) vaut y_n^(a-1) (1 si a = 1).
for k in (1, 2, 3):
    yn = (2 * n * pi + pi / 2)**Rational(-1, k)
    assert simplify(sin(1 / yn**k)) == 1
    assert limit(yn**Rational(-1, 2), n, oo) == oo       # a = 1/2 : y_n^(a-1) → +∞

# ===========================================================================
# Carte ar-l1-derivabilite-reciproque
# ===========================================================================
# Exemple : arctan' = 1/(1+y^2) ; ln' = 1/y.
assert egal(diff(atan(y), y), 1 / (1 + y**2)) and egal(diff(log(h), h), 1 / h)
# Contre-exemples : taux de ∛ en 0 → +∞ ; f(x) = x^3+1 : f(0)=1, f'(1)=3, f'(0)=0.
assert limit(h**Rational(1, 3) / h, h, 0, "+") == oo
assert (x**3 + 1).subs(x, 0) == 1 and diff(x**3 + 1, x).subs(x, 1) == 3 and diff(x**3 + 1, x).subs(x, 0) == 0
# Affirmation 2 : f = x^3+x, f(1) = 2, 1/f'(1) = 1/4 ≠ 1/f'(2) = 1/13.
f = x**3 + x
assert f.subs(x, 1) == 2 and 1 / diff(f, x).subs(x, 1) == Rational(1, 4) and 1 / diff(f, x).subs(x, 2) == Rational(1, 13)

# ar-l1-derivabilite-reciproque-ex1 et variantes : antécédent, 1/f'(a), et l'erreur 1/f'(b).
for f, b, a0, bonne, leurres in [
        (x**5 + x**3 + x, 3, 1, Rational(1, 9), [9, Rational(1, 433), Rational(1, 3)]),                   # ex1
        (2 * x + sin(x), 2 * pi, pi, 1, [Rational(1, 3), 3, Rational(1, 2)]),                              # ex1 v1
        (x + log(x), E + 1, E, E / (E + 1), [(E + 1) / E, (E + 1) / (E + 2), 1 / E])]:                     # ex1 v2
    assert egal(f.subs(x, a0), b)
    assert egal(1 / diff(f, x).subs(x, a0), bonne)
    assert all(not egal(l, bonne) for l in leurres)
    assert egal(1 / diff(f, x).subs(x, b), leurres[1] if f != 2 * x + sin(x) else leurres[0])  # erreur « au mauvais point »

# ar-l1-derivabilite-reciproque-ex2 et variantes : dérivées de arcsin, arccos, argsh.
assert egal(1 / cos(asin(y)), 1 / sqrt(1 - y**2))                                                         # ex2
assert egal(diff(asin(y), y), 1 / sqrt(1 - y**2))
assert egal(diff(acos(y), y), -1 / sqrt(1 - y**2))                                                        # ex2 v1
assert egal(diff(asinh(y), y), 1 / sqrt(1 + y**2)) and egal(cosh(x)**2 - sinh(x)**2, 1)                   # ex2 v2
assert not egal(diff(asinh(y), y), 1 / cosh(y)) and not egal(diff(asin(y), y), 1 / cos(y))

# ar-l1-derivabilite-reciproque-ex3 et variantes : g' = 1/f'(g), g'' = -f''(g)/f'(g)^3.
for f, b, a0, gp, gpp in [(x + exp(x), 1, 0, Rational(1, 2), Rational(-1, 8)),                            # ex3
                          (x**3 + 3 * x, 4, 1, Rational(1, 6), Rational(-1, 36)),                         # ex3 v1
                          (x + log(x), 1, 1, Rational(1, 2), Rational(1, 8))]:                            # ex3 v2
    assert egal(f.subs(x, a0), b)
    f1, f2 = diff(f, x).subs(x, a0), diff(f, x, 2).subs(x, a0)
    assert egal(1 / f1, gp) and egal(-f2 / f1**3, gpp)
# Contrôle de la formule g'' sur une réciproque explicite : f = exp, g = ln : g'' = -1/y^2 = -f''(g)/f'(g)^3.
assert egal(diff(log(h), h, 2), -exp(log(h)) / exp(log(h))**3)

# ar-l1-derivabilite-reciproque-ex4 et variantes : W' = 1/((1+W) e^W).
f = x * exp(x)
assert egal(diff(f, x), (1 + x) * exp(x)) and f.subs(x, -1) == -1 / E
for b, w, bonne, leurres in [(E, 1, 1 / (2 * E), [2 * E, 1 / E, 1 / ((1 + E) * exp(E))]),                     # ex4
                             (2 * E**2, 2, 1 / (3 * E**2), [3 * E**2, 1 / (2 * E**2), 1 / ((1 + 2 * E**2) * exp(2 * E**2))]),  # ex4 v1
                             (-1 / (2 * sqrt(E)), Rational(-1, 2), 2 * sqrt(E), [1 / (2 * sqrt(E)), -2 * sqrt(E), sqrt(E)])]:  # ex4 v2
    assert egal(f.subs(x, w), b)
    assert egal(1 / ((1 + w) * exp(w)), bonne)
    assert all(not egal(l, bonne) for l in leurres)
    if w != 0:
        assert egal(w / (b * (1 + w)), bonne)     # forme W/(y(1+W))
assert egal(1 / diff(f, x).subs(x, E), 1 / ((1 + E) * exp(E)))

# ===========================================================================
# Carte ar-l1-derivabilite-leibniz
# ===========================================================================
# Exemple : (x sin x)^(n) = x sin(x + nπ/2) + n sin(x + (n-1)π/2), contrôlé pour n = 0..8.
for k in range(9):
    assert egal(diff(x * sin(x), x, k), x * sin(x + k * pi / 2) + k * sin(x + (k - 1) * pi / 2))
# Contre-exemple : (x^2)'' = 2 ≠ 0 ; g = x|x| : g' = 2|x| (vérifié de chaque côté).
assert diff(x**2, x, 2) == 2 and egal(diff(h**2, h), 2 * h)

# ar-l1-derivabilite-leibniz-ex1 et variantes : formule fermée = dérivée n-ième pour n = 0..8 ; leurres faux.
for f, bonne, leurres in [
        (x * exp(2 * x), lambda k: 2**(k - 1) * (2 * x + k) * exp(2 * x),
         [lambda k: 2**k * (x + k) * exp(2 * x), lambda k: (2 * x + k) * exp(2 * x), lambda k: 2**k * (2 * x + k) * exp(2 * x)]),   # ex1
        (x * exp(-x), lambda k: (-1)**k * (x - k) * exp(-x),
         [lambda k: (-1)**k * (x + k) * exp(-x), lambda k: (x - k) * exp(-x), lambda k: (-1)**(k + 1) * (x - k) * exp(-x)]),       # ex1 v1
        (x**2 * exp(3 * x), lambda k: S(3)**(k - 2) * (9 * x**2 + 6 * k * x + k * (k - 1)) * exp(3 * x),
         [lambda k: S(3)**(k - 2) * (9 * x**2 + 6 * k * x + 2 * k * (k - 1)) * exp(3 * x),
          lambda k: S(3)**k * (x**2 + 2 * k * x + k * (k - 1)) * exp(3 * x),
          lambda k: S(3)**(k - 1) * (3 * x**2 + 2 * k * x) * exp(3 * x)])]:                                                         # ex1 v2
    for k in range(9):
        assert egal(diff(f, x, k), bonne(k))
    for l in leurres:
        assert any(not egal(diff(f, x, k), l(k)) for k in range(9))

# ar-l1-derivabilite-leibniz-ex2 et variantes.
for f, k, bonne, leurres in [(x**2 * sin(x), 9, -72, [72, -36, 0]),          # ex2
                             (x**2 * cos(x), 10, 90, [-90, 45, 0]),          # ex2 v1
                             (x**3 * exp(2 * x), 5, 240, [40, 480, 120])]:   # ex2 v2
    assert diff(f, x, k).subs(x, 0) == bonne and bonne not in leurres
assert binomial(9, 2) * 2 == 72 and sin(7 * pi / 2) == -1

# ar-l1-derivabilite-leibniz-ex3 et variantes : relations de récurrence et valeurs en 0.
f = atan(x)                                                                                    # ex3
for k in range(1, 7):
    assert egal((1 + x**2) * diff(f, x, k + 1) + 2 * k * x * diff(f, x, k) + k * (k - 1) * diff(f, x, k - 1), 0)
assert diff(f, x, 7).subs(x, 0) == -720 and diff(f, x, 3).subs(x, 0) == -2 and diff(f, x, 5).subs(x, 0) == 24
assert -720 not in (720, -5040, Rational(-1, 7)) and factorial(6) == 720
f = exp(x**2)                                                                                  # ex3 v1
for k in range(1, 7):
    assert egal(diff(f, x, k + 1), 2 * x * diff(f, x, k) + 2 * k * diff(f, x, k - 1))
assert diff(f, x, 6).subs(x, 0) == 120 and diff(f, x, 4).subs(x, 0) == 12 and diff(f, x, 2).subs(x, 0) == 2
f = asin(x)                                                                                    # ex3 v2
assert egal((1 - x**2) * diff(f, x, 2), x * diff(f, x))
for k in range(0, 6):
    assert egal((1 - x**2) * diff(f, x, k + 2) - (2 * k + 1) * x * diff(f, x, k + 1) - k**2 * diff(f, x, k), 0)
assert diff(f, x, 7).subs(x, 0) == 225 and diff(f, x, 5).subs(x, 0) == 9 and diff(f, x, 3).subs(x, 0) == 1

# ar-l1-derivabilite-leibniz-ex4 et variantes : valeurs aux bornes et coefficient dominant, n = 1..7 ; leurres faux.
for k in range(1, 8):
    P = diff((x**2 - 1)**k, x, k)                                                              # ex4
    assert Poly(P, x).degree() == k and Poly(P, x).LC() == factorial(2 * k) / factorial(k)
    assert P.subs(x, 1) == 2**k * factorial(k) and P.subs(x, -1) == (-1)**k * 2**k * factorial(k)
    Q = diff((x**2 - 4)**k, x, k)                                                              # ex4 v1
    assert Q.subs(x, 2) == 4**k * factorial(k) and Q.subs(x, -2) == (-1)**k * 4**k * factorial(k)
    R = diff(x**k * (1 - x)**k, x, k)                                                          # ex4 v2
    assert Poly(R, x).LC() == (-1)**k * factorial(2 * k) / factorial(k)
    assert R.subs(x, 0) == factorial(k) and R.subs(x, 1) == (-1)**k * factorial(k)
# Leurres : faux pour au moins une valeur de n (n = 2).
k = 2
assert 2**k * factorial(k) not in (factorial(k), 2**(k - 1) * factorial(k), factorial(k)**2)
assert 4**k * factorial(k) not in (2**k * factorial(k), factorial(k), 4**(k - 1) * factorial(k))
assert (-1)**k * factorial(k) not in (0, (-1)**k) and (-1)**3 * factorial(3) != factorial(3)

# ===========================================================================
# Carte ar-l1-derivabilite-continuite
# ===========================================================================
# Exemple : taux de x^2 sin(1/x) en 0 → 0 ; contre-exemples : √x (taux → +∞), fonction signe (taux 1/|x|).
assert limit(x * sin(1 / x), x, 0) == 0 and limit(sqrt(h) / h, h, 0, "+") == oo
# Affirmation 1 : |sin(πx)| en un entier k : taux → ±π (calcul en k = 0 de chaque côté).
assert limit(sin(pi * h) / h, h, 0, "+") == pi and limit(sin(pi * h) / (-h), h, 0, "+") == -pi

# ar-l1-derivabilite-continuite-ex1 et variantes : racines simples (g(a) = 0, g'(a) ≠ 0), taux à droite = |g'(a)|.
for g, racines in [(x**2 - 1, {-1, 1}),                 # ex1
                   (x**2 - 4 * x + 3, {1, 3}),          # ex1 v1
                   (x**3 - x, {-1, 0, 1})]:             # ex1 v2
    assert set(solve(g, x)) == racines
    assert all(diff(g, x).subs(x, r) != 0 for r in racines)
assert all(diff(x**2 - 1, x).subs(x, r) in (-2, 2) for r in (-1, 1))
assert (x**2 - 4 * x + 3).subs(x, 2) == -1             # ex1 v1 : en 2, f = -(x-1)(x-3) localement

# ar-l1-derivabilite-continuite-ex2 et variantes : taux |x|^(β-1) avec β = α, 2α, α/2 ; seuil β = 1.
for beta_de_alpha, seuil in [(lambda al: al, 1), (lambda al: 2 * al, Rational(1, 2)), (lambda al: al / 2, 2)]:
    assert beta_de_alpha(seuil) == 1
    assert limit(h**(beta_de_alpha(seuil + 1) - 1), h, 0, "+") == 0              # α > seuil : dérivable
    assert limit(h**(beta_de_alpha(seuil / 2) - 1), h, 0, "+") == oo             # α < seuil : taux → +∞

# ar-l1-derivabilite-continuite-ex3 et variantes : valeurs de F et G aux racines.
F, G = Abs(x + 2) * (x - 2), Abs(x - 2) * (x - 2)                                      # ex3
assert F.subs(x, 2) == 0 and G.subs(x, -2) == -16
F, G = Abs(x - 1) * (x + 1), Abs(x + 1) * (x + 1)                                      # ex3 v1
assert F.subs(x, -1) == 0 and G.subs(x, 1) == 4
F, G = Abs(x + 3) * (x**2 - 9), Abs(x - 3) * (x**2 - 9)                                # ex3 v2
assert F.subs(x, 3) == 0 and G.subs(x, -3) == 0
# Contrôle direct (ex3) : taux de φ en -2 à droite et à gauche : -16 et 16 ; en 2 : 0.
phi = Abs(x**2 - 4) * (x - 2)
assert limit(phi / (x + 2), x, -2, "+") == -16 and limit(phi / (x + 2), x, -2, "-") == 16
assert limit(phi / (x - 2), x, 2, "+") == 0 and limit(phi / (x - 2), x, 2, "-") == 0

# ar-l1-derivabilite-continuite-ex4 Q4 : P(x) = (x-1)^2 (x-2)^2 ne s'annule qu'en 1 et 2 ; P/(x-a) → 0 en a.
P = (x - 1)**2 * (x - 2)**2
assert set(solve(P, x)) == {1, 2}
assert limit(P / (x - 1), x, 1) == 0 and limit(P / (x - 2), x, 2) == 0

print("ar_l1_derivabilite : toutes les vérifications passent.")
