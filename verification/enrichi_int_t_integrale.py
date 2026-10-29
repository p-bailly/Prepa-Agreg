"""Vérifications SymPy des questions de QCM et des exercices ajoutés au chapitre int-t-integrale.

Voir CLAUDE.md, « Contrôle de l'exactitude ». Chaque assert cite l'identifiant vérifié.
"""
from sympy import (E, Rational, cos, exp, integrate, log, N, pi, simplify, sin, sqrt, symbols, diff, oo)

x, t = symbols('x t', real=True)

def distincts(*v):
    vals = [N(a, 30) for a in v]
    return all(abs(vals[i] - vals[j]) > 1e-12 for i in range(len(vals)) for j in range(i))

# int-t-integrale-q1 : ∫_0^π x sin x dx = π ; leurres 0, -π, 2.
I = integrate(x*sin(x), (x, 0, pi))
assert I == pi and distincts(pi, 0, -pi, 2)

# int-t-integrale-q2 : F(x) = ∫_0^x e^{-t²} dt : F' = e^{-x²} > 0 ; F(-1) < 0 < F(1) (donc F non paire).
F1 = N(integrate(exp(-t**2), (t, 0, 1))); Fm1 = N(integrate(exp(-t**2), (t, 0, -1)))
assert F1 > 0 and Fm1 < 0

# int-t-integrale-q3 (vrai) : x^{n+1}/(1+x²) <= x^n/(1+x²) sur [0,1] ; contrôle n = 0..6.
prev = None
for m in range(0, 7):
    v = N(integrate(x**m/(1 + x**2), (x, 0, 1)))
    if prev is not None:
        assert v <= prev
    prev = v

# int-t-integrale-q4 (faux) : valeur moyenne de √x sur [0,4] = 4/3 ≠ 1.
mu = integrate(sqrt(x), (x, 0, 4))/4
assert mu == Rational(4, 3)

# int-t-integrale-q5 (vrai) : ∫_0^{ln 2} e^{2x} dx = 3/2.
assert simplify(integrate(exp(2*x), (x, 0, log(2))) - Rational(3, 2)) == 0

# int-t-integrale-q6 : aire entre y = x² et y = 2 - x² : 8/3 ; leurres 4/3, 4, 2.
A = integrate((2 - x**2) - x**2, (x, -1, 1))
assert A == Rational(8, 3) and distincts(A, Rational(4, 3), 4, 2)

# int-t-integrale-q7 : valeur moyenne de 1/x sur [1,e] = 1/(e-1) ; leurres (e+1)/(2e), 1, e-1.
mu = simplify(integrate(1/x, (x, 1, E))/(E - 1))
assert simplify(mu - 1/(E - 1)) == 0 and distincts(mu, (E + 1)/(2*E), 1, E - 1)

# int-t-integrale-q8 : ∫_0^1 (2x+1) e^x dx = e + 1 ; leurres 3e-1, e-1, 2e.
I = integrate((2*x + 1)*exp(x), (x, 0, 1))
assert simplify(I - (E + 1)) == 0 and distincts(I, 3*E - 1, E - 1, 2*E)

# int-t-integrale-ipp-ex6 : intégrales de Wallis.
def W(m, f=cos):
    return integrate(f(x)**m, (x, 0, pi/2))
assert W(0) == pi/2 and W(1) == 1 and W(0, sin) == pi/2 and W(1, sin) == 1
for m in range(0, 6):
    for f in (cos, sin):
        assert simplify((m + 2)*W(m + 2, f) - (m + 1)*W(m, f)) == 0      # relation de récurrence
        assert simplify((m + 1)*W(m + 1, f)*W(m, f) - pi/2) == 0          # produit constant
        assert N(W(m + 1, f)) <= N(W(m, f))                               # décroissance
# IPP : W_{n+2} = (n+1)∫ sin² cos^n = (n+1)(W_n - W_{n+2}) ; contrôle du terme tout intégré et de l'identité
n = 3
assert simplify(integrate(sin(x)**2*cos(x)**n, (x, 0, pi/2)) - (W(n) - W(n + 2))) == 0
assert W(2) == pi/4 and W(3) == Rational(2, 3) and W(4) == 3*pi/16
assert W(5, sin) == Rational(8, 15) and W(3, sin) == Rational(2, 3)
assert W(6) == 5*pi/32
assert distincts(3*pi/16, 3*pi/8, pi/4, Rational(8, 15))
assert distincts(Rational(8, 15), Rational(4, 15), 3*pi/16, Rational(2, 3))
assert distincts(5*pi/32, 3*pi/16, 5*pi/16, Rational(16, 35))

# int-t-integrale-moyenne-ex5 : moyenne de 2x/(1+x²) sur [0,a] = ln(1+a²)/a ; 0 <= f <= 1 ; ln(1+a²) <= a.
f = 2*x/(1 + x**2)
assert simplify(1 - f - (x - 1)**2/(1 + x**2)) == 0
for a, mu_attendu, leurres in ((2, log(5)/2, (log(5), Rational(2, 5), 2*log(5))),
                               (3, log(10)/3, (log(10), Rational(3, 10), 3*log(10))),
                               (1, log(2), (log(2)/2, Rational(1, 2), 2*log(2)))):
    mu = integrate(f, (x, 0, a))/a
    assert simplify(mu - mu_attendu) == 0 and distincts(mu, *leurres)
    assert (f.subs(x, 0) + f.subs(x, a))/2 == leurres[1]  # leurre « moyenne des valeurs aux bornes »
    assert N(log(1 + a**2)) <= a

# int-t-integrale-proprietes-ex5 : encadrements.
for g, bas, haut, faux in (
        (exp(x**2), Rational(4, 3), E - 1, [(Rational(5, 3), E), (1, Rational(4, 3)), (E - 1, E)]),
        (exp(x**3), Rational(5, 4), E - 1, [(1, Rational(5, 4)), (E - 1, E), (Rational(3, 2), E)]),
        (exp(x**2/2), Rational(7, 6), 2*(sqrt(E) - 1), [(1, Rational(7, 6)), (2*(sqrt(E) - 1), sqrt(E)), (Rational(4, 3), 2)])):
    I = N(integrate(g, (x, 0, 1)), 30).as_real_imag()[0]  # partie imaginaire numérique parasite (≈ 1e-36)
    assert N(bas) <= I <= N(haut)
    for (u, v) in faux:
        assert not (N(u) <= I <= N(v))   # chaque leurre est un encadrement faux
assert integrate(1 + x**2, (x, 0, 1)) == Rational(4, 3)
assert integrate(1 + x**3, (x, 0, 1)) == Rational(5, 4)
assert integrate(1 + x**2/2, (x, 0, 1)) == Rational(7, 6)
assert simplify(integrate(exp(x/2), (x, 0, 1)) - 2*(sqrt(E) - 1)) == 0
assert simplify(integrate(exp(x), (x, 0, 1)) - (E - 1)) == 0
# Relecture indépendante (résolution à l'aveugle) :
# int-t-integrale-q4 : la moyenne des valeurs aux bornes de √x sur [0,4] vaut 1 (le piège), la valeur moyenne 4/3.
assert (sqrt(0) + sqrt(4))/2 == 1 and integrate(sqrt(x), (x, 0, 4))/4 == Rational(4, 3)
# int-t-integrale-moyenne-ex5 variante 2 (correction allégée) : μ = ∫_0^1 2x/(1+x²) dx = ln 2 <= 1.
assert simplify(integrate(2*x/(1 + x**2), (x, 0, 1)) - log(2)) == 0 and N(log(2)) <= 1
print("enrichi_int_t_integrale : tout est vérifié.")
