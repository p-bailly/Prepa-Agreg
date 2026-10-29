"""Vérifications calculatoires du chapitre « ar-l2-suites-fonctions » (analyse réelle, L2) :
cartes, exercices et QCM de content/*/analyse-reelle/l2-ar-l2-suites-fonctions.json et content/qcm/.

Chaque assert cite l'identifiant (sans le préfixe « ar-l2-suites-fonctions- ») de la carte, de l'exercice
(variante comprise : v1, v2) ou de la question. Lancer : python3 verification/ar_l2_suites_fonctions.py
"""

import math
from fractions import Fraction as F

from sympy import (Rational, binomial, diff, exp, expand, integrate, limit, log, oo, pi, simplify, sin, sqrt,
                   symbols, zeta)

n = symbols("n", positive=True, integer=True)
x = symbols("x", nonnegative=True)
t = symbols("t", real=True)
E1 = math.e


def grille(a, b, m=4001):
    return [a + (b - a) * i / (m - 1) for i in range(m)]


def sup_num(f, a, b, m=20001):
    return max(abs(f(v)) for v in grille(a, b, m))


def distinctes(formules, N=200):
    """Formules en n (leurres et bonne réponse) : deux à deux différentes pour chaque n de 1 à N."""
    for k in range(1, N + 1):
        vals = [f(k) for f in formules]
        for i in range(len(vals)):
            for j in range(i + 1, len(vals)):
                assert abs(vals[i] - vals[j]) > 1e-12 * max(1, abs(vals[i])), (k, i, j)


# ---------------------------------------------------------------------------
# Carte « uniforme »
# ---------------------------------------------------------------------------
# Exemple : x/(1+nx) < 1/n sur [0,+inf[, de limite 1/n en +inf (sup non atteint).
assert limit(x / (1 + n * x), x, oo) == 1 / n
assert simplify(1 / n - x / (1 + n * x) - 1 / (n * (1 + n * x))) == 0
# Contre-exemple : n x e^{-nx} maximal en 1/n, valeur e^{-1}.
g = n * x * exp(-n * x)
assert simplify(diff(g, x).subs(x, 1 / n)) == 0 and simplify(g.subs(x, 1 / n) - exp(-1)) == 0
# Pièges et leurres : H_n constantes (écarts consécutifs 1/(n+1) -> 0, divergence) ; f_n = -f (normes égales) ;
# f_n = 1 - 1/(n+1) (normes croissantes, convergence uniforme vers 1) ; x sur [0,1[ bornée sans maximum.
H = [sum(F(1, k) for k in range(1, m + 1)) for m in range(1, 60)]
assert all(H[i + 1] - H[i] == F(1, i + 2) for i in range(58)) and H[-1] > 4
assert all(1 - F(1, k + 1) < 1 - F(1, k + 2) for k in range(50))
# Affirmation 3 : (x + 1/n) x - x^2 = x/n, non borné sur R.
assert expand((x + 1 / n) * x - x**2) == x / n

# uniforme-ex1 (+ v1, v2) : sup = 1/n (sin(nx)/n, atteint en pi/(2n)), 3/n (x/n sur [0,3]), 1/n (1/(n+x^2)).
assert abs(math.sin(5 * math.pi / 10) / 5 - 1 / 5) < 1e-15
for k in (1, 3, 7):
    assert abs(sup_num(lambda v: math.sin(k * v) / k, 0, 2 * math.pi) - 1 / k) < 1e-6
    assert abs(sup_num(lambda v: v / k, 0, 3) - 3 / k) < 1e-12
    assert abs(sup_num(lambda v: 1 / (k + v * v), -5, 5) - 1 / k) < 1e-12
distinctes([lambda k: 1 / k, lambda k: 2 / k, lambda k: 1 / (2 * k), lambda k: math.pi / (2 * k)])
distinctes([lambda k: 3 / k, lambda k: 1 / k, lambda k: 9 / k, lambda k: 3 / (2 * k)])
distinctes([lambda k: 1 / k, lambda k: 1 / (k + 1), lambda k: 2 / k])

# uniforme-ex2 (+ v1, v2) : x e^{-nx} -> 1/(ne) ; x/(1+nx^2) -> 1/(2 sqrt n) ; x^2 e^{-nx} -> 4/(e^2 n^2).
for f, x0, val in [(x * exp(-n * x), 1 / n, 1 / (n * exp(1))),
                   (x / (1 + n * x**2), 1 / sqrt(n), 1 / (2 * sqrt(n))),
                   (x**2 * exp(-n * x), 2 / n, 4 / (exp(2) * n**2))]:
    assert simplify(diff(f, x).subs(x, x0)) == 0 and simplify(f.subs(x, x0) - val) == 0
for k in (1, 2, 5):
    assert abs(sup_num(lambda v: v * math.exp(-k * v), 0, 20) - 1 / (k * E1)) < 1e-7
    assert abs(sup_num(lambda v: v / (1 + k * v * v), 0, 20) - 1 / (2 * math.sqrt(k))) < 1e-7
    assert abs(sup_num(lambda v: v * v * math.exp(-k * v), 0, 20) - 4 / (E1**2 * k * k)) < 1e-7
distinctes([lambda k: 1 / (k * E1), lambda k: 1 / k, lambda k: E1 / k, lambda k: 1 / (2 * k * E1)])
distinctes([lambda k: 1 / (2 * k**0.5), lambda k: 1 / k**0.5, lambda k: 1 / (4 * k**0.5), lambda k: 1 / (2 * k) ** 0.5])
distinctes([lambda k: 4 / (E1**2 * k * k), lambda k: 4 / k**2, lambda k: 2 / (E1 * k), lambda k: 1 / (E1**2 * k * k)])

# uniforme-ex3 (+ v1, v2) : n x^n (1-x) -> (n/(n+1))^{n+1} -> 1/e ; v1 n x^n (1-x)^2 -> 0 ; v2 n^2 x^n (1-x) -> +inf.
X = symbols("X", positive=True)
assert simplify(diff(X**n * (1 - X), X) - X**(n - 1) * (n - (n + 1) * X)) == 0
assert simplify(diff(X**n * (1 - X)**2, X) - X**(n - 1) * (1 - X) * (n - (n + 2) * X)) == 0
assert limit((n / (n + 1))**(n + 1), n, oo) == exp(-1)
assert limit(4 * n / (n + 2)**2 * (n / (n + 2))**n, n, oo) == 0
assert limit(n * (n / (n + 1))**(n + 1), n, oo) == oo
for k in (3, 10):
    assert abs(sup_num(lambda v: k * v**k * (1 - v), 0, 1) - (k / (k + 1))**(k + 1)) < 1e-7
    assert abs(sup_num(lambda v: k * v**k * (1 - v)**2, 0, 1) - 4 * k / (k + 2)**2 * (k / (k + 2))**k) < 1e-7
    assert 4 * k / (k + 2)**2 <= 4 / k

# uniforme-ex4 : |sin u - u| <= |u|^3/6 ; f_n(n pi) - n pi = -n pi.
assert all(abs(math.sin(u) - u) <= abs(u)**3 / 6 + 1e-15 for u in grille(-20, 20, 40001))
assert abs(5 * math.sin(math.pi) - 5 * math.pi + 5 * math.pi) < 1e-12

# uniforme-ex5 : u - u^2/2 <= ln(1+u) <= u ; 1 - e^{-t} <= t ; e^n - 2^n -> +inf.
assert all(u - u * u / 2 <= math.log1p(u) <= u for u in grille(0, 50, 20001))
assert all(1 - math.exp(-s) <= s for s in grille(0, 30, 3001))
for a in (0.5, 2.0):
    for k in (1, 4, 30):
        assert sup_num(lambda v: math.exp(v) - (1 + v / k)**k, 0, a) <= a * a * math.exp(a) / (2 * k) + 1e-12
assert math.exp(10) - 2**10 > 20000

# ---------------------------------------------------------------------------
# Carte « continuite »
# ---------------------------------------------------------------------------
# Exemple et contre-exemples : 0 <= sqrt(x^2+1/n^2) - |x| <= 1/n, égalité en 0.
for k in (1, 4, 9):
    assert abs(sup_num(lambda v: math.sqrt(v * v + 1 / k**2) - abs(v), -3, 3) - 1 / k) < 1e-12
# Trou 3, leurre « lipschitzienne » (carte integration, trou 2) : sqrt(x+1/n) - sqrt(x) <= 1/sqrt(n).
for k in (1, 4, 16):
    assert sup_num(lambda v: math.sqrt(v + 1 / k) - math.sqrt(v), 0, 1) <= 1 / math.sqrt(k) + 1e-12

# continuite-ex1 (+ v1, v2) : limites simples.
assert abs(2.0**40 / (1 + 2.0**40) - 1) < 1e-11 and 0.5**40 / (1 + 0.5**40) < 1e-11
assert abs(math.atan(1e6) - math.pi / 2) < 1e-5 and abs(1e6 * 3 / (1 + 1e6 * 3) - 1) < 1e-6

# continuite-ex2 (+ v1, v2) : sup de sqrt(x^2+c) - |x| = sqrt(c) : 1/sqrt(n), 2/n, 3/sqrt(n).
for k in (1, 2, 9):
    for c, val in [(1 / k, 1 / math.sqrt(k)), (4 / k**2, 2 / k), (9 / k, 3 / math.sqrt(k))]:
        assert abs(sup_num(lambda v: math.sqrt(v * v + c) - abs(v), -4, 4) - val) < 1e-12
distinctes([lambda k: k**-0.5, lambda k: 0.5 * k**-0.5, lambda k: 2 * k**-0.5])
distinctes([lambda k: 2 / k, lambda k: 1 / k, lambda k: 4 / k, lambda k: math.sqrt(2) / k])
distinctes([lambda k: 3 * k**-0.5, lambda k: math.sqrt(3) * k**-0.5, lambda k: 9 * k**-0.5, lambda k: k**-0.5])

# continuite-ex3 : f_n(1/n) = 1, dérivée, f_n(a) <= 2/(na).
fn = 2 * n * x / (1 + n**2 * x**2)
assert simplify(fn.subs(x, 1 / n)) == 1
assert simplify(diff(fn, x) - 2 * n * (1 - n**2 * x**2) / (1 + n**2 * x**2)**2) == 0

# continuite-ex4 (+ v1, v2) : (1-1/n)^n -> 1/e, (1-2/n)^n -> e^{-2}, (1-1/n^2)^n -> 1.
assert limit((1 - 1 / n)**n, n, oo) == exp(-1)
assert limit((1 - 2 / n)**n, n, oo) == exp(-2)
assert limit((1 - 1 / n**2)**n, n, oo) == 1

# ---------------------------------------------------------------------------
# Carte « integration »
# ---------------------------------------------------------------------------
# Contre-exemples : int_0^1 n^2 x e^{-nx} = 1 - (n+1) e^{-n} -> 1 ; int_0^inf e^{-x/n}/n = 1 ; int_0^1 x^n = 1/(n+1).
I2 = integrate(n**2 * x * exp(-n * x), (x, 0, 1))
assert simplify(I2 - (1 - (n + 1) * exp(-n))) == 0 and limit(I2, n, oo) == 1
assert simplify(integrate(exp(-x / n) / n, (x, 0, oo))) == 1
assert simplify(integrate(x**n, (x, 0, 1)) - 1 / (n + 1)) == 0
# Piège « segments de ]a,b] » : sur [c,1] (c>0), sup n^2 x e^{-nx} = n^2 c e^{-nc} -> 0 dès que 1/n < c.
assert limit(n**2 * Rational(1, 10) * exp(-n / 10), n, oo) == 0
# Exemple : int_0^1 arctan = pi/4 - ln2/2.
from sympy import atan  # noqa: E402
assert simplify(integrate(atan(x), (x, 0, 1)) - (pi / 4 - log(2) / 2)) == 0
# Trou 3 : f_n = f + 1/n sur [0,2] (écart 2/n > 1/n et > 1/(2n)) ; f_n = -f, f = 1 sur [0,1] (écart 2 > 0).
assert F(2, 7) > F(1, 7) > F(1, 14)

# integration-ex1 (+ v1, v2).
assert simplify((n * x + 1) / (n + x) - x - (1 - x**2) / (n + x)) == 0
assert integrate(x, (x, 0, 1)) == Rational(1, 2) and integrate(x, (x, 0, 2)) == 2
assert integrate(sin(x + pi / 2), (x, 0, pi / 2)) == 1
assert all(abs(k * math.log1p(v / k) - v) <= v * v / (2 * k) + 1e-13 for k in (1, 3, 10) for v in grille(0, 2, 201))
for k in (1, 5):
    assert sup_num(lambda v: (math.exp(v / k) - 1) * math.cos(v), 0, math.pi / 2) <= math.exp(math.pi / (2 * k)) - 1
assert abs(math.log(2) - 0.5) > 0.1

# integration-ex2 (+ v1, v2).
J1 = integrate(n**2 * x * (1 - x)**n, (x, 0, 1))
assert simplify(J1 - n**2 / ((n + 1) * (n + 2))) == 0 and limit(J1, n, oo) == 1
J2 = integrate(n**3 * x * exp(-n * x), (x, 0, 1))
assert simplify(J2 - n * (1 - (n + 1) * exp(-n))) == 0 and limit(J2, n, oo) == oo
assert simplify((n**2 * x * exp(-n * x)).subs(x, 1 / n) - n / exp(1)) == 0

# integration-ex3 (+ v1, v2) : sum 1/(k 2^k) = ln 2, sum 1/(k 3^k) = ln(3/2), sum (-1)^k/((2k+1)2^{2k+1}) = arctan(1/2).
s1 = sum(1 / (k * 2**k) for k in range(1, 80))
s2 = sum(1 / (k * 3**k) for k in range(1, 80))
s3 = sum((-1)**k / ((2 * k + 1) * 2**(2 * k + 1)) for k in range(0, 60))
assert abs(s1 - math.log(2)) < 1e-15 and abs(s2 - math.log(1.5)) < 1e-15 and abs(s3 - math.atan(0.5)) < 1e-15
for bonne, leurres in [(math.log(2), [1, math.log(1.5), 0.5]), (math.log(1.5), [math.log(3), 0.5, math.log(2)]),
                       (math.atan(0.5), [math.log(1.25), math.pi / 8, 0.5 * math.log(1.5)])]:
    assert all(abs(bonne - v) > 1e-3 for v in leurres)
assert simplify(1 / (1 + t**2) - sum((-t**2)**k for k in range(6)) - (-1)**6 * t**12 / (1 + t**2)) == 0

# integration-ex4 : f_{n+1} - f_n = x^{n+1}/(n+1)! ; sum 1/k! converge (vers e).
assert abs(sum(1 / math.factorial(k) for k in range(25)) - math.e) < 1e-15

# ---------------------------------------------------------------------------
# Carte « derivation »
# ---------------------------------------------------------------------------
# Contre-exemples : f_n'(x) = x/sqrt(x^2+1/n^2) ; (sin(nx)/n)' = cos(nx) ; (x|x|)' = 2|x| (sur x>0 : x^2 -> 2x).
assert simplify(diff(sqrt(t**2 + 1 / n**2), t) - t / sqrt(t**2 + 1 / n**2)) == 0
from sympy import cos  # noqa: E402
assert simplify(diff(sin(n * t) / n, t) - cos(n * t)) == 0
assert diff(t**2, t) == 2 * t
# Trou 2, leurre « segments évitant c » : sur [d,1], 1 - x/sqrt(x^2+1/n^2) <= 1 - d/sqrt(d^2+1/n^2) -> 0.
assert limit(1 - Rational(1, 10) / sqrt(Rational(1, 100) + 1 / n**2), n, oo) == 0

# derivation-ex1 (+ v1, v2) : f_n'(0) = 1 ; sqrt(n) ; -sin(n pi/2) = -1, 0, 1, 0, ...
assert [round(-math.sin(k * math.pi / 2)) for k in range(1, 5)] == [-1, 0, 1, 0]

# derivation-ex2 (+ v1, v2) : f'(0) = zeta(2) = pi^2/6 ; v1 f'(0) = 0 ; v2 f'(pi) = sum (-1)^k/k^2 = -pi^2/12.
assert zeta(2) == pi**2 / 6
alt = sum((-1)**k / k**2 for k in range(1, 200001))
assert abs(alt + math.pi**2 / 12) < 1e-9
assert abs(math.pi**2 / 6 - float(zeta(3))) > 0.1 and abs(math.pi**2 / 12 - float(zeta(3))) > 0.1

# derivation-ex3 : f_n(x) = x/(1+n^2x^2), max 1/(2n) en 1/n, f_n'(0) = 1.
h = x / (1 + n**2 * x**2)
assert simplify(diff(h, x) - (1 - n**2 * x**2) / (1 + n**2 * x**2)**2) == 0
assert simplify(h.subs(x, 1 / n) - 1 / (2 * n)) == 0 and diff(h, x).subs(x, 0) == 1

# derivation-ex4 (+ v1, v2) : f'(0) = pi^2/12, f'(1) = 1 - pi^2/12, (v2) f'(0) = pi^2/48.
d0 = sum((-1)**(k + 1) / k**2 for k in range(1, 200001))
d1 = sum((-1)**(k + 1) / (k + 1)**2 for k in range(1, 200001))
d2 = sum((-1)**(k + 1) / (2 * k)**2 for k in range(1, 200001))
assert abs(d0 - math.pi**2 / 12) < 1e-9 and abs(d1 - (1 - math.pi**2 / 12)) < 1e-9 and abs(d2 - math.pi**2 / 48) < 1e-9
for bonne, leurres in [(math.pi**2 / 12, [-math.pi**2 / 12, math.pi**2 / 6, math.log(2)]),
                       (1 - math.pi**2 / 12, [math.pi**2 / 12 - 1, math.pi**2 / 12, math.pi**2 / 6 - 1]),
                       (math.pi**2 / 48, [math.pi**2 / 12, -math.pi**2 / 48, math.pi**2 / 24])]:
    assert all(abs(bonne - v) > 1e-2 for v in leurres)

# ---------------------------------------------------------------------------
# Carte « normale »
# ---------------------------------------------------------------------------
# Exemple : sum_{k>n} 1/k^2 <= 1/n.
for m in (1, 5, 20):
    assert sum(F(1, k * k) for k in range(m + 1, 3000)) < F(1, m)
# Trou 3 : u_k = 2^{-k} (reste 2^{-n} > 2^{-n-1}) ; u_k = 10.2^{-k} (reste en n=0 : 10 > 1) ;
# carrés : sum_{k>=1} 4^{-k} = 1/3 < 1 = reste en n=0 pour u_k = 2^{-k}.
assert F(1, 3) < 1

# normale-ex1 (+ v1, v2) : 1/n^2 (oui), 1/(n 2^n) (oui), 1/n (non).
distinctes([lambda k: 1 / k**2, lambda k: 1 / (k * 2**k)])

# normale-ex2 (+ v1, v2) : x/(n(1+nx^2)) -> 1/(2 n sqrt n) ; x/(1+n^2x^2) -> 1/(2n) ; x^2 e^{-nx^2} -> 1/(en).
for f, x0, val in [(x / (n * (1 + n * x**2)), 1 / sqrt(n), 1 / (2 * n * sqrt(n))),
                   (x / (1 + n**2 * x**2), 1 / n, 1 / (2 * n)),
                   (x**2 * exp(-n * x**2), 1 / sqrt(n), 1 / (n * exp(1)))]:
    assert simplify(diff(f, x).subs(x, x0)) == 0 and simplify(f.subs(x, x0) - val) == 0
for k in (1, 3):
    assert abs(sup_num(lambda v: v * v * math.exp(-k * v * v), 0, 10) - 1 / (E1 * k)) < 1e-7
distinctes([lambda k: 1 / (2 * k**1.5), lambda k: 1 / k**0.5])
distinctes([lambda k: 1 / (2 * k), lambda k: 1 / k, lambda k: 1 / (4 * k * k)])
distinctes([lambda k: 1 / (E1 * k), lambda k: 1 / (E1**2 * k * k), lambda k: 2 / (E1 * k)])

# normale-ex3 : sup_{x>1} n^{-x} = 1/n.
assert all(abs(sup_num(lambda v: k**-v, 1, 3) - 1 / k) < 1e-15 for k in (1, 2, 7))

# normale-ex4 : |R_n(x)| <= 1/(n+1+x^2) (séries alternées), vérifié numériquement.
for xx in (0.0, 0.7, 3.0):
    total = sum((-1)**k / (k + xx * xx) for k in range(1, 400001))
    for m in (1, 2, 10):
        reste = total - sum((-1)**k / (k + xx * xx) for k in range(1, m + 1))
        assert abs(reste) <= 1 / (m + 1 + xx * xx) + 1e-5

# ---------------------------------------------------------------------------
# Carte « cauchy-uniforme »
# ---------------------------------------------------------------------------
# Contre-exemple : x^n = 1/2 donne x^n - x^{2n} = 1/4 ; troncatures décimales de sqrt(2).
assert F(1, 2) - F(1, 4) == F(1, 4)
tr = [F(math.isqrt(2 * 10**(2 * p)), 10**p) for p in range(12)]
assert all(abs(tr[p] - tr[q]) < F(1, 10**min(p, q)) for p in range(12) for q in range(12) if p != q)
# Affirmation 4, cauchy-uniforme-ex2 : sum_{k=N+1}^{2N} 1/k >= 1/2.
assert all(sum(F(1, k) for k in range(N + 1, 2 * N + 1)) >= F(1, 2) for N in range(1, 60))
# Trous 2 et 3 : |H_p - H_q| <= |p-q|/(min(p,q)+1) <= eps |p-q| pour p, q >= N >= 1/eps.
assert all(abs(H[p] - H[q]) <= F(abs(p - q), min(p, q) + 2) for p in range(40) for q in range(40))

# cauchy-uniforme-ex1 (+ v1, v2) : sup |x^n| = 1 ; sup n^2x^2e^{-nx} = 4/e^2 ; sup nxe^{-nx} = 1/e.
assert abs(sup_num(lambda s: s * s * math.exp(-s), 0, 20) - 4 / E1**2) < 1e-7
assert abs(sup_num(lambda s: s * math.exp(-s), 0, 20) - 1 / E1) < 1e-7

# cauchy-uniforme-ex3 : formule de S_n, majoration 1/sin(a/2), identité d'Abel, minoration sqrt(2)/4.
for xx in (0.3, 1.0, 2.5, 4.0):
    for m in (1, 5, 13):
        S = sum(math.sin(k * xx) for k in range(1, m + 1))
        assert abs(S - math.sin(m * xx / 2) * math.sin((m + 1) * xx / 2) / math.sin(xx / 2)) < 1e-12
a = 0.4
assert all(abs(sum(math.sin(k * v) for k in range(1, m + 1))) <= 1 / math.sin(a / 2) + 1e-12
           for v in grille(a, 2 * math.pi - a, 400) for m in range(1, 60))
Sx = lambda m, v: sum(math.sin(k * v) for k in range(1, m + 1))  # noqa: E731
for v in (0.5, 2.0):
    for p, q in [(1, 2), (3, 9), (5, 30)]:
        gauche = sum(math.sin(k * v) / k for k in range(p + 1, q + 1))
        droite = Sx(q, v) / q - Sx(p, v) / (p + 1) + sum(Sx(k, v) * (1 / k - 1 / (k + 1)) for k in range(p + 1, q))
        assert abs(gauche - droite) < 1e-12
for N in range(1, 80):
    v = math.pi / (4 * N)
    assert sum(math.sin(k * v) / k for k in range(N + 1, 2 * N + 1)) >= math.sqrt(2) / 4

# ---------------------------------------------------------------------------
# Carte « approximation »
# ---------------------------------------------------------------------------
# Exemple : corde de x^2 entre a et b : (a+b)x - ab - x^2 = (x-a)(b-x) <= (b-a)^2/4.
aa, bb = symbols("a b", real=True)
assert expand((aa + bb) * t - aa * bb - t**2 - (t - aa) * (bb - t)) == 0
# Bernstein : moments binomiaux (ordre 0, 1, 2) ; B_n(t^2) = t^2 + t(1-t)/n (approximation-ex2).
for m in (1, 2, 5, 8):
    s0 = sum(binomial(m, k) * t**k * (1 - t)**(m - k) for k in range(m + 1))
    s1 = sum(k * binomial(m, k) * t**k * (1 - t)**(m - k) for k in range(m + 1))
    s2 = sum(k**2 * binomial(m, k) * t**k * (1 - t)**(m - k) for k in range(m + 1))
    assert expand(s0) == 1 and expand(s1 - m * t) == 0 and expand(s2 - m * t * (1 - t) - m**2 * t**2) == 0
    B = lambda f: expand(sum(binomial(m, k) * f(Rational(k, m)) * t**k * (1 - t)**(m - k) for k in range(m + 1)))  # noqa: E731
    assert expand(B(lambda u: u**2) - t**2 - t * (1 - t) / m) == 0
    assert expand(B(lambda u: 2 * u**2 - u) - (2 * t**2 - t) - 2 * t * (1 - t) / m) == 0
    assert expand(B(lambda u: 3 * u**2) - 3 * t**2 - 3 * t * (1 - t) / m) == 0
distinctes([lambda k: 1 / (4 * k), lambda k: 1 / k, lambda k: 1 / (2 * k), lambda k: 1 / (8 * k)])
distinctes([lambda k: 1 / (2 * k), lambda k: 1 / (4 * k), lambda k: 1 / k, lambda k: 3 / (4 * k)])
distinctes([lambda k: 3 / (4 * k), lambda k: 3 / k, lambda k: 1 / (4 * k), lambda k: 3 / (2 * k)])
# Contre-exemples : sin(1/x) prend 1 et -1 en des points arbitrairement proches (pas uniformément continue).
for k in (10, 100):
    p1, p2 = 1 / (2 * k * math.pi + math.pi / 2), 1 / (2 * k * math.pi + 3 * math.pi / 2)
    assert abs(math.sin(1 / p1) - 1) < 1e-9 and abs(math.sin(1 / p2) + 1) < 1e-9 and abs(p1 - p2) < 1e-3

# approximation-ex1 (+ v1, v2) : sup de f - phi_n = (2n-1)/n^2 (x^2), 1/sqrt(n) (sqrt), (3n^2-3n+1)/n^3 (x^3).
for k in (1, 2, 5, 10):
    for f, val in [(lambda v: v * v, (2 * k - 1) / k**2), (math.sqrt, 1 / math.sqrt(k)),
                   (lambda v: v**3, (3 * k * k - 3 * k + 1) / k**3)]:
        ecart = max(f((j + 1) / k) - f(j / k) for j in range(k))  # borne supérieure sur chaque morceau
        assert abs(ecart - val) < 1e-12
distinctes([lambda k: (2 * k - 1) / k**2, lambda k: (2 * k + 1) / k**2, lambda k: 2 / k, lambda k: (k - 1) / k**2])
distinctes([lambda k: k**-0.5, lambda k: 0.5 * k**-0.5, lambda k: (math.sqrt(2) - 1) * k**-0.5, lambda k: 2 * k**-0.5])
distinctes([lambda k: (3 * k * k - 3 * k + 1) / k**3, lambda k: (3 * k * k + 3 * k + 1) / k**3, lambda k: 3 / k,
            lambda k: (3 * k * k - 3 * k) / k**3])

# ---------------------------------------------------------------------------
# QCM
# ---------------------------------------------------------------------------
# q1 : x/(n+x) -> 1 en +inf (sup non atteint), vaut 1/2 en x = n.
assert limit(x / (n + x), x, oo) == 1 and simplify((x / (n + x)).subs(x, n)) == Rational(1, 2)
# q3 : écart <= pi^2/n, limite int_0^pi sin = 2.
assert integrate(sin(x), (x, 0, pi)) == 2
for k in (1, 10):
    assert sup_num(lambda v: math.sin(v) - math.sin(v) / (1 + v * v / k), 0, math.pi) <= math.pi**2 / k
# q4 : sup x^{n+1}/(n+1) = 1/(n+1) sur [0,1] ; dérivée x^n.
assert simplify(diff(x**(n + 1) / (n + 1), x) - x**n) == 0
# q8 : max de x^n(1-x) = n^n/(n+1)^{n+1} <= 1/(n+1).
assert all(F(k**k, (k + 1)**(k + 1)) <= F(1, k + 1) for k in range(1, 40))

print("ar_l2_suites_fonctions : toutes les vérifications passent.")
