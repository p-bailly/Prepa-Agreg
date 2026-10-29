"""Vérifications SymPy des questions de QCM et des exercices ajoutés au chapitre int-l1-riemann.

Voir CLAUDE.md, « Contrôle de l'exactitude ». Chaque assert cite l'identifiant vérifié.
"""
from sympy import (E, Rational, cos, exp, floor, integrate, log, N, pi, simplify, sin, sqrt, symbols, diff,
                   Piecewise, factorial, apart, Integral)

import math

x, t = symbols('x t', real=True)

def distincts(*v):
    vals = [N(a, 30) for a in v]
    return all(abs(vals[i] - vals[j]) > 1e-12 for i in range(len(vals)) for j in range(i))

def riemann(f, m):
    # somme de Riemann (points de droite) calculée en flottants
    return sum(f(k/m) for k in range(1, m + 1))/m

# int-l1-riemann-q1 : Σ k/(n²+k²) = (1/n) Σ (k/n)/(1+(k/n)²) -> ∫_0^1 x/(1+x²) = ln2/2.
I = integrate(x/(1 + x**2), (x, 0, 1))
assert simplify(I - log(2)/2) == 0 and distincts(I, log(2), pi/4, Rational(1, 2))
m = 4000
assert abs(sum(k/(m**2 + k**2) for k in range(1, m + 1)) - float(N(I))) < 1e-3

# int-l1-riemann-q2 : Σ (k/n²) e^{k/n} -> ∫_0^1 x e^x = 1 ; leurres e-1, e, 1/2.
I = integrate(x*exp(x), (x, 0, 1))
assert I == 1 and distincts(I, E - 1, E, Rational(1, 2))
assert abs(riemann(lambda u: u*math.exp(u), 2000) - 1) < 2e-3

# int-l1-riemann-q3 : ∫_1^4 dx/(x+√x) = 2 ln(3/2) (x = t², t ∈ [1,2]) ; leurres ln(3/2), 2 ln 3, ln(5/2).
J = integrate(2*t/(t**2 + t), (t, 1, 2))
assert simplify(J - 2*log(Rational(3, 2))) == 0
assert abs(N(Integral(1/(x + sqrt(x)), (x, 1, 4))) - N(J)) < 1e-12
assert distincts(J, log(Rational(3, 2)), 2*log(3), log(Rational(5, 2)))

# int-l1-riemann-q4 (vrai) : ∫_0^{π/2} sin/(sin+cos) = π/4.
assert abs(N(Integral(sin(t)/(sin(t) + cos(t)), (t, 0, pi/2))) - N(pi/4)) < 1e-12
assert simplify(sin(pi/2 - t) - cos(t)) == 0

# int-l1-riemann-q5 (faux) : f(t) = t - 1/2 continue, non nulle, ∫_0^1 f³ = 0.
assert integrate((t - Rational(1, 2))**3, (t, 0, 1)) == 0

# int-l1-riemann-q6 (vrai) : (∫ t f)² <= ∫t² ∫f² = (1/3)∫f² (Cauchy-Schwarz) ; égalité pour f(t) = t.
assert integrate(t**2, (t, 0, 1)) == Rational(1, 3)
for f in (t, 1 + 0*t, exp(t), t**2 - 3*t + 1, sin(5*t)):
    g = integrate(t*f, (t, 0, 1))**2; d = integrate(f**2, (t, 0, 1))/3
    assert N(d - g) >= -1e-14

# int-l1-riemann-q7 : ∫_0^1 dx/(x²+3x+2) = ln(4/3) ; leurres ln(3/2), ln 6, ln(2/3).
assert apart(1/(x**2 + 3*x + 2)) == 1/(x + 1) - 1/(x + 2)
I = integrate(1/(x**2 + 3*x + 2), (x, 0, 1))
assert simplify(I - log(Rational(4, 3))) == 0 and distincts(I, log(Rational(3, 2)), log(6), log(Rational(2, 3)))

# int-l1-riemann-q8 (faux) : sin(1/x) n'a pas de limite en 0⁺ (valeurs 1 et -1 sur deux suites tendant vers 0).
for k in range(0, 5):
    assert simplify(sin(1/(1/(pi/2 + 2*k*pi)))) == 1
    assert simplify(sin(1/(1/(3*pi/2 + 2*k*pi)))) == -1

# int-l1-riemann-cpm-ex5 : ∫_0^2 ⌊x²⌋ = 5 - √2 - √3.
I = 1*(sqrt(2) - 1) + 2*(sqrt(3) - sqrt(2)) + 3*(2 - sqrt(3))
assert simplify(I - (5 - sqrt(2) - sqrt(3))) == 0
# contrôle par une somme de Riemann fine (la fonction est en escalier)
mm = 200000
approx = sum(int((2*k/mm)**2) for k in range(mm))*2/mm
assert abs(approx - float(N(I))) < 1e-4
assert distincts(I, 6 - sqrt(2) - sqrt(3), Rational(8, 3), 3 - sqrt(2))
# variante 1 : ∫_0^2 ⌊e^x⌋ = 14 - ln(7!) = 14 - ln 5040.
I1 = sum(k*(log(k + 1) - log(k)) for k in range(1, 7)) + 7*(2 - log(7))
assert simplify(I1 - (14 - log(5040))) == 0 and factorial(7) == 5040
approx = sum(math.floor(math.exp(2*k/20000)) for k in range(20000))*2/20000
assert abs(approx - float(N(I1))) < 2e-3
assert distincts(I1, 14 - log(720), exp(2) - 1, 13 - log(5040))
# variante 2 : ∫_0^1 ⌊3x²⌋ = 2 - (√6 + √3)/3.
I2 = 1*(sqrt(Rational(2, 3)) - 1/sqrt(3)) + 2*(1 - sqrt(Rational(2, 3)))
assert simplify(I2 - (2 - (sqrt(6) + sqrt(3))/3)) == 0
approx = sum(int(3*(k/mm)**2) for k in range(mm))/mm
assert abs(approx - float(N(I2))) < 1e-4
assert distincts(I2, 1, 3 - (sqrt(6) + sqrt(3))/3, 2 - sqrt(6)/3)

# int-l1-riemann-sommes-ex5 : Jensen intégral (preuve) ; contrôle numérique sur des exemples.
for f, phi in ((t, exp), (sin(3*t), lambda u: u**2), (log(1 + t), exp)):
    gauche = phi(integrate(f, (t, 0, 1))); droite = integrate(phi(f), (t, 0, 1))
    assert N(droite - gauche) >= -1e-14
# application : ∫ ln g <= ln ∫ g, exemple g = 1 + t : 2 ln 2 - 1 <= ln(3/2).
assert simplify(integrate(log(1 + t), (t, 0, 1)) - (2*log(2) - 1)) == 0
assert N(log(Rational(3, 2)) - (2*log(2) - 1)) > 0

# int-l1-riemann-nulle-ex5, question 3 : f(t) = 2t - 1 ne convient pas (|∫f| = 0 ≠ ∫|f| = 1/2).
assert integrate(2*t - 1, (t, 0, 1)) == 0
assert integrate(1 - 2*t, (t, 0, Rational(1, 2))) + integrate(2*t - 1, (t, Rational(1, 2), 1)) == Rational(1, 2)
# Relecture indépendante (résolution à l'aveugle) :
# int-l1-riemann-sommes-ex5, question 3, exemple f(t) = t : (∫t)² = 1/4 <= ∫t² = 1/3.
assert integrate(t, (t, 0, 1))**2 == Rational(1, 4) and integrate(t**2, (t, 0, 1)) == Rational(1, 3)
# int-l1-riemann-q1 : la limite de Σ n/(n²+k²) (citée dans l'explication) est ∫_0^1 dx/(1+x²) = π/4.
assert simplify(integrate(1/(1 + x**2), (x, 0, 1)) - pi/4) == 0
print("enrichi_int_l1_riemann : tout est vérifié.")
