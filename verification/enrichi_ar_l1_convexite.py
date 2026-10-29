"""Vérifications SymPy des questions de QCM et des exercices ajoutés au chapitre ar-l1-convexite.

Voir CLAUDE.md, « Contrôle de l'exactitude ». Chaque assert cite l'identifiant vérifié.
"""
from sympy import (Rational, sqrt, symbols, diff, simplify, solve, exp, log, sin, cos, pi, integrate, N,
                   factor, Max, Abs, discriminant, Interval, solveset, S, oo, E)

x, a, b, c = symbols('x a b c', real=True)

# ar-l1-convexite-q1 : f = 1/(1+x²), f'' = (6x²-2)/(1+x²)³ ; sur [0,+oo[, f'' >= 0 exactement sur [1/√3,+oo[.
f = 1/(1 + x**2)
f2 = simplify(diff(f, x, 2))
assert simplify(f2 - (6*x**2 - 2)/(1 + x**2)**3) == 0
assert solveset(6*x**2 - 2 >= 0, x, Interval(0, oo)) == Interval(1/sqrt(3), oo)
# leurres : [1,+oo[ est contenu (pas le plus grand), [0,+oo[ et [1/3,+oo[ contiennent des points où f'' < 0.
assert f2.subs(x, 0) < 0 and f2.subs(x, Rational(1, 2)) < 0 and Rational(1, 2) > Rational(1, 3)

# ar-l1-convexite-q2 : (x e^x)'' = (x+2)e^x change de signe en -2.
assert simplify(diff(x*exp(x), x, 2) - (x + 2)*exp(x)) == 0

# ar-l1-convexite-q3 (vrai) : ((a+b+c)/3)² <= (a²+b²+c²)/3, l'écart vaut ((a-b)²+(b-c)²+(c-a)²)/9.
assert simplify((a**2 + b**2 + c**2)/3 - ((a + b + c)/3)**2 - ((a - b)**2 + (b - c)**2 + (c - a)**2)/9) == 0

# ar-l1-convexite-q4 (faux) : x = 1, y = 3 : ln 2 > (ln 1 + ln 3)/2.
assert N(log(2) - log(3)/2) > 0

# ar-l1-convexite-q5 : (a+b+c)(1/a+1/b+1/c) >= 9, égalité pour a = b = c ; leurres 3, 6, 1 non atteints (minimum 9).
expr = (a + b + c)*(1/a + 1/b + 1/c)
assert expr.subs({a: 1, b: 1, c: 1}) == 9
import random
random.seed(1)
for _ in range(2000):
    va, vb, vc = (Rational(random.randint(1, 1000), random.randint(1, 1000)) for _ in range(3))
    assert expr.subs({a: va, b: vb, c: vc}) >= 9

# ar-l1-convexite-q6 : f(0)=0, f(1)=1 convexe => f(3) >= 3 ; f(3) = 4 atteint par max(x, (3x-1)/2).
g = Max(x, (3*x - 1)/2)
assert g.subs(x, 0) == 0 and g.subs(x, 1) == 1 and g.subs(x, 3) == 4
for v in (2, Rational(5, 2), 0):
    assert v < 3  # leurres impossibles

# ar-l1-convexite-q7 (vrai) : f(0)=f(1)=0, x>1 : pente(0,1) = 0 <= pente(1,x) donc f(x) >= 0 (exemple x(x-1)).
assert all((t*(t - 1)) >= 0 for t in range(1, 50))

# ar-l1-convexite-q8 (faux) : ln(1+x²)'' = 2(1-x²)/(1+x²)², négative en 2 ; corde : ln 5 > (ln 2 + ln 10)/2.
h2 = simplify(diff(log(1 + x**2), x, 2))
assert simplify(h2 - 2*(1 - x**2)/(1 + x**2)**2) == 0 and h2.subs(x, 2) < 0
assert N(log(5) - (log(2) + log(10))/2) > 0

# ar-l1-convexite-derivables-ex5 : sin'' = -sin <= 0 sur [0,π/2] ; corde y = 2x/π ; intégrales.
n = symbols('n', positive=True)
assert integrate(exp(-n*x), (x, 0, pi/2)).equals((1 - exp(-n*pi/2))/n)
assert simplify(integrate(exp(-2*n*x/pi), (x, 0, pi/2)) - pi/(2*n)*(1 - exp(-n))) == 0
for i in range(0, 201):
    xv = pi/2*Rational(i, 200)
    assert N(sin(xv) - 2*xv/pi) >= -1e-15 and N(xv - sin(xv)) >= -1e-15
# encadrement numérique de I_n pour quelques n
for m in (1, 3, 10):
    I = N(integrate(exp(-m*sin(x)), (x, 0, pi/2)))
    assert N((1 - exp(-m*pi/2))/m) <= I <= N(pi/(2*m))
# leurres de la réponse finale : faux
assert N(sin(1) - 2/pi) > 0                    # sin x <= 2x/π faux en x = 1 (sin 1 > 2/π)
assert N(1 - pi/2/pi - sin(pi/2)) < 0          # sin(π/2) = 1 > 1 - 1/2
assert N(sin(Rational(1, 2)) - 2*Rational(1, 2)/pi) > 0  # sin x <= 2x/π faux en x = 1/2

# ar-l1-convexite-pentes-ex5 : meilleur encadrement de f(1), minorations de f(5) et f(-1).
def verifie(f0, f2_, f3, bas1, haut1, min5, minm1, g_haut, g_bas):
    s02 = Rational(f2_ - f0, 2); s23 = f3 - f2_
    assert s02 <= s23
    assert Rational(f0 + f2_, 2) == haut1 and f2_ - s23 == bas1
    assert f3 + 2*s23 == min5 and f0 - s02 == minm1
    for g, v1 in ((g_haut, haut1), (g_bas, bas1)):
        assert (g.subs(x, 0), g.subs(x, 1), g.subs(x, 2), g.subs(x, 3)) == (f0, v1, f2_, f3)
    assert g_haut.subs(x, 5) == min5 and g_haut.subs(x, -1) == minm1
# exercice : f(0)=2, f(2)=0, f(3)=1 : -1 <= f(1) <= 1, f(5) >= 3, f(-1) >= 3.
verifie(2, 0, 1, -1, 1, 3, 3, Max(2 - x, x - 2), Max(2 - 3*x, x - 2))
# variante 1 : f(0)=0, f(2)=4, f(3)=7 : 1 <= f(1) <= 2, f(5) >= 13, f(-1) >= -2.
verifie(0, 4, 7, 1, 2, 13, -2, Max(2*x, 3*x - 2), Max(x, 3*x - 2))
# variante 2 (remplacée à la relecture : l'ancienne, f(0)=4, f(2)=2, f(3)=3, n'était que l'exercice décalé de 2) :
# f(0)=1, f(2)=-1, f(3)=2 : -4 <= f(1) <= 0, f(5) >= 8, f(-1) >= 2.
verifie(1, -1, 2, -4, 0, 8, 2, Max(1 - x, 3*x - 7), Max(1 - 5*x, 3*x - 7))
# leurres de la variante 2 : -1 <= f(1) et -2 <= f(1) sont contredits par la borne -4 atteinte ;
# f(1) <= 1 est vrai mais pas le meilleur (0 est atteint).
assert Max(1 - 5*x, 3*x - 7).subs(x, 1) == -4 and Max(1 - x, 3*x - 7).subs(x, 1) == 0
# (les fonctions max de deux fonctions affines sont convexes)

# ar-l1-convexite-derivables-ex6 : f_a = x⁴ + a x³ + c x², f'' = 12x² + 6ax + 2c >= 0 sur R ssi Δ <= 0.
for cc, borne in ((6, 4), (24, 8), (Rational(3, 2), 2)):
    fa = x**4 + a*x**3 + cc*x**2
    d2 = diff(fa, x, 2)
    D = discriminant(d2, x)
    assert solveset(D <= 0, a, S.Reals) == Interval(-borne, borne)
    # aux bornes, f'' est un carré parfait (convexe), au-delà f'' prend une valeur < 0
    assert factor(d2.subs(a, borne)) == 12*(x + Rational(borne, 4))**2 or simplify(d2.subs(a, borne) - 12*(x + Rational(borne, 4))**2) == 0
    xm = -Rational(6*(borne + 1), 24)
    assert d2.subs({a: borne + 1, x: xm}) < 0
# ar-l1-convexite-derivables-ex5, réponse finale (remplacée à la relecture : l'ancienne demandait l'encadrement de
# sin x, déjà écrit dans l'énoncé) : 1 - e^{-π/2} <= n I_n <= π/2 pour tout n >= 1.
# Leurres faux : n I_1 ≈ 0,873 < 1 ; n I_2 ≈ 1,075 > 1 ; n I_n <= π/2 contredit « π/2 <= n I_n » (n I_1 < π/2).
from mpmath import mp, quad, exp as mexp, sin as msin, pi as mpi
mp.dps = 30
def nI(m):
    return m*quad(lambda s: mexp(-m*msin(s)), [0, mpi/2])
for m in range(1, 40):
    assert 1 - mexp(-mpi/2) <= nI(m) <= mpi/2
assert nI(1) < 1 and nI(2) > 1 and nI(1) < mpi/2
print("enrichi_ar_l1_convexite : tout est vérifié.")
