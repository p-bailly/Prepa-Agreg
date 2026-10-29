"""Vérifications SymPy des questions de QCM et des exercices ajoutés au chapitre ar-l1-taylor.

Voir CLAUDE.md, « Contrôle de l'exactitude ». Chaque assert cite l'identifiant vérifié.
"""
from sympy import (E, Rational, cos, exp, integrate, limit, log, nsimplify, oo, pi, series, simplify, sin,
                   sqrt, symbols, diff, Abs, N, summation, Product, factorial, Symbol)

x, t, u = symbols('x t u', real=True)
n, k = symbols('n k', positive=True, integer=True)

# ar-l1-taylor-q1 : lim (cos x - 1 + x²/2)/x⁴ = 1/24 ; leurres -1/24, 1/12, 0 faux.
L = limit((cos(x) - 1 + x**2/2)/x**4, x, 0)
assert L == Rational(1, 24)
assert all(L != v for v in (Rational(-1, 24), Rational(1, 12), 0))

# ar-l1-taylor-q2 : e^x cos x = 1 + x + 0·x² + o(x²) ; leurres (coefficient de x² 1, -1, 1/2) faux.
s = series(exp(x)*cos(x), x, 0, 3).removeO()
assert simplify(s - (1 + x)) == 0
assert s.coeff(x, 2) == 0

# ar-l1-taylor-q3 : coefficient de x² = f''(0)/2 = 3, donc f''(0) = 6 (et non 3, 3/2 ni -6 = f'''(0)).
f = 2 - x + 3*x**2 - x**3
assert diff(f, x, 2).subs(x, 0) == 6
assert diff(f, x, 3).subs(x, 0) == -6

# ar-l1-taylor-q4 (vrai) : 0 <= e^x - 1 - x <= e x²/2 sur [0,1]. Taylor-Lagrange, |exp''| <= e sur [0,1].
# Contrôle numérique sur une grille fine, et minoration par la convexité (e^x >= 1 + x).
for i in range(0, 1001):
    xv = Rational(i, 1000)
    g = exp(xv) - 1 - xv
    assert N(g) >= 0 and N(E*xv**2/2 - g) >= 0

# ar-l1-taylor-q5 (faux) : ln(1+x) - x = -x²/2 + o(x²) < 0 près de 0 ; exemple x = 0,1.
assert series(log(1 + x) - x, x, 0, 3).removeO() == -x**2/2
assert N(log(Rational(11, 10))) < Rational(1, 10)

# ar-l1-taylor-q6 : cos x + x²/2 = 1 + x⁴/24 + o(x⁴) : minimum local strict en 0.
s = series(cos(x) + x**2/2, x, 0, 6).removeO()
assert simplify(s - (1 + x**4/24)) == 0

# ar-l1-taylor-q7 (vrai) : cos x - (1 - x²/2) = ∫_0^x (x-t)(1 - cos t) dt (reste intégral à l'ordre 1).
R = integrate((x - t)*(1 - cos(t)), (t, 0, x))
assert simplify(R - (cos(x) - 1 + x**2/2)) == 0
for i in range(-500, 501):
    xv = Rational(i, 50)
    assert N(cos(xv) - 1 + xv**2/2) >= 0

# ar-l1-taylor-q8 : plus petite constante C telle que |√(1+x) - 1 - x/2| <= C x² pour x >= 0 : C = 1/8.
g = sqrt(1 + x) - 1 - x/2
assert diff(sqrt(1 + x), x, 2) == -1/(4*(1 + x)**Rational(3, 2))
assert limit(Abs(g)/x**2, x, 0, '+') == Rational(1, 8)  # aucune constante < 1/8 ne convient
for i in range(1, 2001):
    xv = Rational(i, 100)
    assert N(Abs(g.subs(x, xv)) - xv**2/8) <= 0      # 1/8 convient
# leurres 1/4, 1/2 (pas les plus petites), 1/16 (fausse près de 0)
assert N(Abs(g.subs(x, Rational(1, 100))) - Rational(1, 100)**2/16) > 0

# ar-l1-taylor-reste-integral-ex5, question 1 : ln(1+x) = x - ∫_0^x (x-t)/(1+t)² dt, reste dans [0, x²/2].
R = integrate((x - t)/(1 + t)**2, (t, 0, x))
assert simplify(diff(x - R - log(1 + x), x)) == 0 and (x - R - log(1 + x)).subs(x, 0) == 0
assert abs(N((x - R - log(1 + x)).subs(x, Rational(7, 3)))) < 1e-12
assert integrate(x - t, (t, 0, x)) == x**2/2
# question 2 : bornes de ln P_n.
S1 = summation(k, (k, 1, n)); S2 = summation(k**2, (k, 1, n)); S4 = summation(k**4, (k, 1, n))
haut = S1/n**2
bas = haut - S2/(2*n**4)
assert limit(haut, n, oo) == Rational(1, 2) and limit(bas, n, oo) == Rational(1, 2)
# valeurs numériques : P_n -> √e
Pn = lambda m, a: N(Product(1 + a(k, m), (k, 1, m)).doit())
assert abs(Pn(2000, lambda kk, m: kk/m**2) - N(sqrt(E))) < 1e-3
# variante 1 : ∏(1 + 2k/n²) -> e
haut = 2*S1/n**2; bas = haut - 4*S2/(2*n**4)
assert limit(haut, n, oo) == 1 and limit(bas, n, oo) == 1
assert abs(Pn(2000, lambda kk, m: 2*kk/m**2) - N(E)) < 2e-3
# variante 2 : ∏(1 + k²/n³) -> e^(1/3) ; Σ k⁴/(2n⁶) <= 1/(2n) -> 0
haut = S2/n**3; bas = haut - S4/(2*n**6)
assert limit(haut, n, oo) == Rational(1, 3) and limit(bas, n, oo) == Rational(1, 3)
assert abs(Pn(2000, lambda kk, m: kk**2/m**3) - N(exp(Rational(1, 3)))) < 2e-3
# leurres distincts de la bonne réponse
assert len({N(sqrt(E)), N(E), 1.0, N(exp(Rational(1, 3)))}) == 4
assert len({N(E), N(sqrt(E)), N(E**2), 1.0}) == 4
assert len({N(exp(Rational(1, 3))), N(sqrt(E)), N(E), 1.0}) == 4
# question 3 : Q_n = ∏ (n+k+1)/(n+k) = (2n+1)/(n+1) -> 2, et 0 <= S_n - ln Q_n <= Σ 1/(2(n+k)²) <= 1/(2n).
for m in (1, 2, 5, 17):
    Q = 1
    for kk in range(1, m + 1):
        Q *= Rational(m + kk + 1, m + kk)
    assert Q == Rational(2*m + 1, m + 1)
    reste = sum(Rational(1, 2*(m + kk)**2) for kk in range(1, m + 1))
    assert reste <= Rational(1, 2*m)
assert limit((2*n + 1)/(n + 1), n, oo) == 2

# ar-l1-taylor-dl-usuels-ex5 : lim (x - sin x)/(x(1 - cos x)) = 1/3 ; leurres 1/6, -1/3, 1/2 faux.
assert limit((x - sin(x))/(x*(1 - cos(x))), x, 0) == Rational(1, 3)
# variante 1 : lim (ln(1+x) - x)/(1 - cos x) = -1 ; leurres 1, -1/2, 0.
assert limit((log(1 + x) - x)/(1 - cos(x)), x, 0) == -1
# variante 2 : lim (e^x - 1 - sin x)/(1 - cos x) = 1 ; leurres 1/2, 0, 2.
assert limit((exp(x) - 1 - sin(x))/(1 - cos(x)), x, 0) == 1

# ar-l1-taylor-dl-usuels-ex6 : x e^{1/x} = x + 1 + 1/(2x) + o(1/x) en +oo : asymptote y = x+1, courbe au-dessus.
h = symbols('h', positive=True)
F = (1/h)*exp(h)  # x = 1/h
s = series(F, h, 0, 2).removeO()
assert simplify(s - (1/h + 1 + h/2)) == 0
assert limit(x*exp(1/x) - x - 1, x, oo) == 0
assert limit(x*(x*exp(1/x) - x - 1), x, oo) == Rational(1, 2)
# variante 1 : (x+1)e^{1/x} = x + 2 + 3/(2x) + o(1/x) : asymptote y = x+2, au-dessus.
assert limit((x + 1)*exp(1/x) - x - 2, x, oo) == 0
assert limit(x*((x + 1)*exp(1/x) - x - 2), x, oo) == Rational(3, 2)
# variante 2 : √(x² + 2x) = x + 1 - 1/(2x) + o(1/x) : asymptote y = x+1, au-dessous.
assert limit(sqrt(x**2 + 2*x) - x - 1, x, oo) == 0
assert limit(x*(sqrt(x**2 + 2*x) - x - 1), x, oo) == Rational(-1, 2)

print("enrichi_ar_l1_taylor : tout est vérifié.")
