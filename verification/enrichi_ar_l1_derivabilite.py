"""Vérifications calculatoires de l'enrichissement du chapitre « ar-l1-derivabilite » (analyse réelle, L1) :
questions de content/qcm/ar-l1-derivabilite.json et exercices ajoutés à
content/exercices/analyse-reelle/l1-ar-l1-derivabilite.json.

Chaque assert cite l'identifiant de la question ou de l'exercice (variante comprise : v1, v2).
Lancer : python3 verification/enrichi_ar_l1_derivabilite.py (doit finir sans erreur).
"""

from sympy import (Abs, Poly, Rational, atan, cos, diff, exp, expand, limit, log, oo, pi, real_roots, simplify, sin,
                   sqrt, symbols)

x = symbols("x", real=True)

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-q1 : f(x) = x^3 e^{2x}, f^(5)(0) = C(5,3) 3! 2^2 = 240 ; leurres 480, 60, 120.
# ---------------------------------------------------------------------------
v = diff(x**3 * exp(2 * x), x, 5).subs(x, 0)
assert v == 240 and all(v != l for l in (480, 60, 120))

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-q2 : f(x) = x + sin x ; f'(pi) = 0 donc f^-1 n'est pas dérivable en f(pi) = pi.
# ---------------------------------------------------------------------------
f = x + sin(x)
assert f.subs(x, pi) == pi and diff(f, x).subs(x, pi) == 0
assert simplify(diff(f, x) - (1 + cos(x))) == 0  # f' >= 0, nulle seulement en pi + 2k pi : f strictement croissante
# taux d'accroissement de f^-1 en pi : (f^-1(pi + h) - pi)/h = t/(f(pi + t) - pi) avec t -> 0 ; f(pi+t) - pi = t - sin t
t = symbols("t", positive=True)
assert limit(t / (t - sin(t)), t, 0) == oo  # tangente verticale : leurres 1/2, 1, 0 faux

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-q3 (vrai) : g(x) = (e^x - 1)/x, g(0) = 1 : g'(x) -> 1/2, g de classe C^1.
# ---------------------------------------------------------------------------
g = (exp(x) - 1) / x
assert limit(g, x, 0) == 1
gp = simplify(diff(g, x))
assert simplify(gp - (x * exp(x) - exp(x) + 1) / x**2) == 0
assert limit(gp, x, 0) == Rational(1, 2) and limit((g - 1) / x, x, 0) == Rational(1, 2)

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-q5 (faux) : arctan dérivable, de dérivée > 0, mais limite pi/2 en +oo.
# ---------------------------------------------------------------------------
assert simplify(diff(atan(x), x) - 1 / (1 + x**2)) == 0 and limit(atan(x), x, oo) == pi / 2

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-q6 (faux) : |x| et -|x| non dérivables en 0, somme nulle dérivable.
# ---------------------------------------------------------------------------
assert limit(Abs(x) / x, x, 0, "+") == 1 and limit(Abs(x) / x, x, 0, "-") == -1

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-q7 : TAF pour sqrt sur [1, 4] : 1 = 3/(2 sqrt c), c = 9/4 ; leurres 5/2, 2, 3/2.
# ---------------------------------------------------------------------------
c = Rational(9, 4)
assert sqrt(4) - sqrt(1) == 3 * diff(sqrt(x), x).subs(x, c)
for l in (Rational(5, 2), 2, Rational(3, 2)):
    assert 3 * diff(sqrt(x), x).subs(x, l) != 1

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-q8 : x|x| : dérivée 2|x|, continue, non dérivable en 0.
# ---------------------------------------------------------------------------
xp = symbols("xp", positive=True)
assert diff(xp**2, xp) == 2 * xp and diff(-(xp**2), xp) == -2 * xp  # x|x| = x^2 (x > 0), -x^2 (x < 0)
assert limit((x * Abs(x)) / x, x, 0) == 0  # dérivable en 0, f'(0) = 0
assert limit(2 * Abs(x) / x, x, 0, "+") == 2 and limit(2 * Abs(x) / x, x, 0, "-") == -2  # f' non dérivable en 0

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-taf-ex5 et variantes : règle de L'Hôpital.
# ---------------------------------------------------------------------------
# exercice : (x - sin x)/x^3 -> 1/6 ; f'/g' = (1 - cos x)/(3x^2), puis sin x/(6x) -> 1/6
assert limit((x - sin(x)) / x**3, x, 0) == Rational(1, 6)
assert limit(sin(x) / (6 * x), x, 0) == Rational(1, 6)
assert all(l != Rational(1, 6) for l in (Rational(1, 3), 0, Rational(-1, 6)))
# v1 : (e^x - 1 - x)/x^2 -> 1/2 (des deux côtés) ; f'/g' = (e^x - 1)/(2x)
assert limit((exp(x) - 1 - x) / x**2, x, 0, "+") == Rational(1, 2) == limit((exp(x) - 1 - x) / x**2, x, 0, "-")
assert limit((exp(x) - 1) / (2 * x), x, 0) == Rational(1, 2)
assert all(l != Rational(1, 2) for l in (1, 0, Rational(1, 6)))
# v2 : (x - arctan x)/x^3 -> 1/3 ; f'/g' = 1/(3(1+x^2))
assert limit((x - atan(x)) / x**3, x, 0) == Rational(1, 3)
assert simplify(diff(x - atan(x), x) / diff(x**3, x) - 1 / (3 * (1 + x**2))) == 0
assert all(l != Rational(1, 3) for l in (Rational(1, 6), Rational(-1, 3), 1))
# question 1 : h = (f(b)-f(a)) g - (g(b)-g(a)) f vérifie h(a) = h(b)
F_, G_ = symbols("F G", cls=__import__("sympy").Function)
a, b = symbols("a b", real=True)
h = (F_(b) - F_(a)) * G_(x) - (G_(b) - G_(a)) * F_(x)
assert simplify(h.subs(x, a) - h.subs(x, b)) == 0

# ---------------------------------------------------------------------------
# ar-l1-derivabilite-leibniz-ex5 et variantes : dérivée n-ième d'un produit polynôme de degré 2 x exponentielle.
# ---------------------------------------------------------------------------
def verifie(expr, formule, leurres, N=9):
    for k in range(N):
        assert simplify(diff(expr, x, k) - formule(k)) == 0
    for L in leurres:
        assert any(simplify(diff(expr, x, k) - L(k)) != 0 for k in range(4))


# exercice : (x^2 + 1) e^{-x} -> (-1)^n e^{-x} (x^2 - 2nx + n^2 - n + 1)
verifie((x**2 + 1) * exp(-x),
        lambda k: (-1) ** k * exp(-x) * (x**2 - 2 * k * x + k**2 - k + 1),
        [lambda k: (-1) ** k * exp(-x) * (x**2 + 2 * k * x + k**2 - k + 1),
         lambda k: exp(-x) * (x**2 - 2 * k * x + k**2 - k + 1),
         lambda k: (-1) ** k * exp(-x) * (x**2 - 2 * k * x + k**2 + k + 1)])
# v1 : (x^2 - x) e^{2x} -> 2^{n-2} e^{2x} (4x^2 + 4(n-1)x + n^2 - 3n)
verifie((x**2 - x) * exp(2 * x),
        lambda k: Rational(2) ** (k - 2) * exp(2 * x) * (4 * x**2 + 4 * (k - 1) * x + k**2 - 3 * k),
        [lambda k: 2**k * exp(2 * x) * (x**2 + (k - 1) * x + k**2 - 3 * k),
         lambda k: Rational(2) ** (k - 2) * exp(2 * x) * (4 * x**2 + 4 * (k - 1) * x + k**2 - k),
         lambda k: Rational(2) ** (k - 2) * exp(2 * x) * (4 * x**2 - 4 * x + k**2 - 3 * k)])
# v2 : (x^2 + 2x) e^{-x} -> (-1)^n e^{-x} (x^2 + 2(1-n)x + n^2 - 3n)
verifie((x**2 + 2 * x) * exp(-x),
        lambda k: (-1) ** k * exp(-x) * (x**2 + 2 * (1 - k) * x + k**2 - 3 * k),
        [lambda k: exp(-x) * (x**2 + 2 * (1 - k) * x + k**2 - 3 * k),
         lambda k: (-1) ** k * exp(-x) * (x**2 + 2 * (1 + k) * x + k**2 + k),
         lambda k: (-1) ** k * exp(-x) * (x**2 + 2 * (1 - k) * x + k**2 - k)])

# ---------------------------------------------------------------------------
# ar-l1-rolle-ex5 : Q_n = ((x^2-1)^n)^(n) a n racines réelles distinctes dans ]-1, 1[ ;
# les dérivées d'ordre k < n de (x^2-1)^n s'annulent en -1 et 1.
# ---------------------------------------------------------------------------
for m in range(1, 8):
    Q = (x**2 - 1) ** m
    for k in range(m):
        dk = diff(Q, x, k)
        assert dk.subs(x, 1) == 0 and dk.subs(x, -1) == 0
    Pn = Poly(expand(diff(Q, x, m)), x)
    assert Pn.degree() == m
    rr = real_roots(Pn)
    assert len(rr) == m and len(set(rr)) == m and all(-1 < r.evalf() < 1 for r in rr)
# question 1 : un polynôme scindé à racines simples ; exemple (x-1)(x-2)(x-4)(x-7) : P' a 3 racines dans les intervalles
P = Poly((x - 1) * (x - 2) * (x - 4) * (x - 7), x)
rp = [r.evalf() for r in real_roots(P.diff(x))]
assert len(rp) == 3 and 1 < rp[0] < 2 < rp[1] < 4 < rp[2] < 7

# ===========================================================================
# Relecture indépendante (2026-10-29) : réponses retrouvées à l'aveugle, recalculées autrement.
# ===========================================================================
from sympy import Symbol as _Sym

_x = _Sym("_x", real=True)

# ar-l1-derivabilite-leibniz-ex5 : formules vérifiées par dérivation directe pour n = 0..12 ; chaque leurre
# est faux pour au moins un n.
_cas = [
    ((_x**2 + 1) * exp(-_x),
     lambda k: (-1) ** k * exp(-_x) * (_x**2 - 2 * k * _x + k**2 - k + 1),
     [lambda k: (-1) ** k * exp(-_x) * (_x**2 + 2 * k * _x + k**2 - k + 1),
      lambda k: exp(-_x) * (_x**2 - 2 * k * _x + k**2 - k + 1),
      lambda k: (-1) ** k * exp(-_x) * (_x**2 - 2 * k * _x + k**2 + k + 1)]),
    ((_x**2 - _x) * exp(2 * _x),
     lambda k: Rational(2) ** (k - 2) * exp(2 * _x) * (4 * _x**2 + 4 * (k - 1) * _x + k**2 - 3 * k),
     [lambda k: Rational(2) ** k * exp(2 * _x) * (_x**2 + (k - 1) * _x + k**2 - 3 * k),
      lambda k: Rational(2) ** (k - 2) * exp(2 * _x) * (4 * _x**2 + 4 * (k - 1) * _x + k**2 - k),
      lambda k: Rational(2) ** (k - 2) * exp(2 * _x) * (4 * _x**2 - 4 * _x + k**2 - 3 * k)]),
    ((_x**2 + 2 * _x) * exp(-_x),
     lambda k: (-1) ** k * exp(-_x) * (_x**2 + 2 * (1 - k) * _x + k**2 - 3 * k),
     [lambda k: exp(-_x) * (_x**2 + 2 * (1 - k) * _x + k**2 - 3 * k),
      lambda k: (-1) ** k * exp(-_x) * (_x**2 + 2 * (1 + k) * _x + k**2 + k),
      lambda k: (-1) ** k * exp(-_x) * (_x**2 + 2 * (1 - k) * _x + k**2 - k)]),
]
for _f, _bonne, _leurres in _cas:
    _d = [diff(_f, _x, k) for k in range(13)]
    assert all(simplify(_d[k] - _bonne(k)) == 0 for k in range(13))
    for _l in _leurres:
        assert any(simplify(_d[k] - _l(k)) != 0 for k in range(13))

# ar-l1-taf-ex5 : limites 1/6, 1/2, 1/3 ; quotients des dérivées utilisés dans la correction.
assert limit((_x - sin(_x)) / _x**3, _x, 0) == Rational(1, 6)
assert limit((1 - cos(_x)) / (3 * _x**2), _x, 0) == Rational(1, 6)
assert limit((exp(_x) - 1 - _x) / _x**2, _x, 0) == Rational(1, 2)
assert limit((_x - atan(_x)) / _x**3, _x, 0) == Rational(1, 3)
assert simplify(diff(_x - atan(_x), _x) / diff(_x**3, _x) - 1 / (3 * (1 + _x**2))) == 0
# question 1 : h(a) = h(b) pour h = (f(b)-f(a)) g - (g(b)-g(a)) f (symbolique)
_fa, _fb, _ga, _gb = symbols("_fa _fb _ga _gb")
assert expand(((_fb - _fa) * _ga - (_gb - _ga) * _fa) - ((_fb - _fa) * _gb - (_gb - _ga) * _fb)) == 0

# ar-l1-rolle-ex5 : pour n = 1..8, Q^(k)(±1) = 0 pour k < n et Q^(n) a n racines réelles distinctes dans ]-1,1[.
for _n in range(1, 9):
    _Q = Poly((_x**2 - 1) ** _n, _x)
    assert all(_Q.diff((_x, k)).eval(1) == 0 and _Q.diff((_x, k)).eval(-1) == 0 for k in range(_n))
    _r = real_roots(_Q.diff((_x, _n)))
    assert len(_r) == _n and len(set(_r)) == _n and all(-1 < r.evalf() < 1 for r in _r)

# ar-l1-derivabilite-q1 : f^(5)(0) = 240 pour x³ e^{2x}.
assert diff(_x**3 * exp(2 * _x), _x, 5).subs(_x, 0) == 240
# ar-l1-derivabilite-q7 : c = 9/4 pour sqrt sur [1,4].
assert 3 / (2 * sqrt(Rational(9, 4))) == 1

print("enrichi_ar_l1_derivabilite : toutes les vérifications passent.")
