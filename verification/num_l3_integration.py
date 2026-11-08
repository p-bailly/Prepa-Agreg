"""Vérifications calculatoires du chapitre « num-l3-integration » (méthodes numériques, L3) : cartes, exercices et
QCM de content/*/methodes-numeriques/l3-num-l3-integration.json et content/qcm/num-l3-integration.json.

Conventions du chapitre : h = (b-a)/n, x_k = a + kh ; R_n rectangles à gauche, M_n point milieu, T_n trapèzes,
S_n = (T_n + 2 M_n)/3 Simpson à n PANNEAUX (2n+1 points). Chaque assert cite l'identifiant de la carte, de l'exercice
(variante comprise : v1, v2) ou de la question. Calculs en fractions exactes quand c'est possible, sinon en SymPy ou
avec mpmath à 40 chiffres. Lancer : python3 verification/num_l3_integration.py.
"""

from fractions import Fraction as F
from math import ceil, comb

import mpmath as mp
import sympy as sp

mp.mp.dps = 40
CH = "num-l3-integration"

# ---------------------------------------------------------------------------
# Méthodes composées, génériques (fractions, mpmath ou SymPy selon les arguments).
# ---------------------------------------------------------------------------


def R(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + k * h) for k in range(n))


def Rd(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + k * h) for k in range(1, n + 1))


def M(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + (2 * k + 1) * h / 2) for k in range(n))


def T(f, a, b, n):
    h = (b - a) / n
    return h * (f(a) / 2 + sum(f(a + k * h) for k in range(1, n)) + f(b) / 2)


def S(f, a, b, n):
    return (T(f, a, b, n) + 2 * M(f, a, b, n)) / 3


def romberg(f, a, b, kmax):
    Rt = [[T(f, a, b, 2 ** k)] for k in range(kmax + 1)]
    for k in range(1, kmax + 1):
        for j in range(1, k + 1):
            Rt[k].append((4 ** j * Rt[k][j - 1] - Rt[k - 1][j - 1]) / (4 ** j - 1))
    return Rt


def plus_petit_n(seuil, p):
    """Plus petit entier n >= 1 tel que n**p >= seuil (seuil exact : Fraction ou SymPy)."""
    seuil = sp.nsimplify(seuil) if isinstance(seuil, F) else sp.sympify(seuil)
    n = max(1, int(sp.N(seuil ** sp.Rational(1, p), 50)) - 2)
    while bool(sp.Integer(n) ** p < seuil):
        n += 1
    return n


t = sp.symbols("t", real=True)
x, s, u, r, h = sp.symbols("x s u r h", positive=True)
e = sp.E
fl = lambda q: mp.mpf(q.numerator) / q.denominator  # noqa: E731

# ===========================================================================
# Carte num-l3-integration-rectangles
# ===========================================================================
c, al, be, hh = sp.symbols("c alpha beta hh", real=True)
# ∫ (t-c) sur un segment centré en c est nul ; ∫ (t-c)^2 = h^3/12 ; h^3/12 * M2/2 = M2 h^3/24.
assert sp.expand(sp.integrate(t - c, (t, c - hh / 2, c + hh / 2))) == 0
assert sp.simplify(sp.integrate((t - c) ** 2, (t, c - hh / 2, c + hh / 2)) - hh ** 3 / 12) == 0
# Rectangles : ∫ (t - x_k) = h^2/2.
assert sp.simplify(sp.integrate(t - al, (t, al, al + hh)) - hh ** 2 / 2) == 0
# Optimalité : f(x)=x, I - R_n = (b-a)^2/(2n) ; f(x)=x^2, I - M_n = (b-a)^3/(12 n^2) = majoration avec M2 = 2.
for (a_, b_) in [(F(0), F(1)), (F(-1), F(2)), (F(1, 3), F(5, 2))]:
    for n in range(1, 7):
        assert (b_ ** 2 - a_ ** 2) / 2 - R(lambda y: y, a_, b_, n) == (b_ - a_) ** 2 / (2 * n)
        assert (b_ ** 3 - a_ ** 3) / 3 - M(lambda y: y * y, a_, b_, n) == (b_ - a_) ** 3 / (12 * n * n)
        assert (b_ - a_) ** 3 * 2 / (24 * n * n) == (b_ - a_) ** 3 / (12 * n * n)
        # rectangles à droite : même majoration (égalité en valeur absolue pour f(x) = x)
        assert abs((b_ ** 2 - a_ ** 2) / 2 - Rd(lambda y: y, a_, b_, n)) == (b_ - a_) ** 2 / (2 * n)
# Exactitude : point milieu exact pour les fonctions affines.
for n in range(1, 5):
    assert M(lambda y: 3 * y - F(7, 2), F(-2), F(3), n) == sp.Rational(1) * F(3, 2) * (9 - 4) - F(7, 2) * 5
# exemple : e - 1, n = 4.
I_e = mp.e - 1
R4, M4, T4 = R(mp.exp, mp.mpf(0), mp.mpf(1), 4), M(mp.exp, mp.mpf(0), mp.mpf(1), 4), T(mp.exp, mp.mpf(0), mp.mpf(1), 4)
assert abs(R4 - mp.mpf("1.512437")) < 1e-6 and abs(I_e - R4 - mp.mpf("0.206")) < 1e-3 and abs(mp.e / 8 - mp.mpf("0.340")) < 1e-3
assert abs(M4 - mp.mpf("1.713815")) < 1e-6 and abs(I_e - M4 - mp.mpf("4.5e-3")) < 5e-5 and abs(mp.e / 384 - mp.mpf("7.1e-3")) < 5e-5
assert round(float((I_e - R4) / (I_e - M4))) == 46
# pièges et contre-exemples : x^3 sur [0,1], f''(0)=0 mais I - M_1 = 1/8 ; x^2 sur [0,1/4] : 1/768 > 1/3072.
assert F(1, 4) - M(lambda y: y ** 3, F(0), F(1), 1) == F(1, 8)
assert F(1, 3) * F(1, 4) ** 3 - M(lambda y: y * y, F(0), F(1, 4), 1) == F(1, 768)
assert F(1, 4) ** 3 * F(1, 2) / 24 == F(1, 3072) and F(1, 768) > F(1, 3072)
assert M(lambda y: y * y, F(0), F(1), 1) == F(1, 4)  # pas exact pour x^2
for n in range(1, 6):  # f(x)=x : R_n = (n-1)/(2n), I - R_n = 1/(2n) > 1/(4n) (affirmation fausse)
    assert R(lambda y: y, F(0), F(1), n) == F(n - 1, 2 * n)
# trous : leurres réfutés. Trou 1 : h Σ f(a+kh) (f(x)=x), (h/2) Σ f(milieux) (f=1), Σ_{k=1}^{n} (f(x)=x, n=1).
assert R(lambda y: y, F(0), F(1), 1) != F(1, 2)
assert F(1, 2) * 1 * 1 != 1
assert 1 * (F(0) + F(3, 2)) != F(1, 2)
# Trou 2 : /48 (x^2, [0,1]) ; (b-a)^2 (x^2 sur [0,2], n=1) ; n^3 (x^2, [0,1], n=2).
assert F(1, 12) > F(2, 48)
assert F(8, 3) - M(lambda y: y * y, F(0), F(2), 1) == F(2, 3) and F(2, 3) > F(4 * 2, 24)
assert F(1, 3) - M(lambda y: y * y, F(0), F(1), 2) == F(1, 48) and F(1, 48) > F(2, 24 * 8)
# Trou 3 : pas exact pour x^2 ni x^3 ; exact pour x (« constantes seulement » faux).
assert M(lambda y: y ** 3, F(0), F(1), 1) != F(1, 4) and M(lambda y: y, F(0), F(1), 1) == F(1, 2)
# affirmation convexe : M_n <= I pour e^x.
assert M4 <= I_e

# --- exercices de la carte rectangles ---
# ex1 : x^2 sur [0,1] ; v1 : x^3 sur [0,2] ; v2 : 1/x sur [1,3] (n = 2).
for (f, a_, b_, Rv, Mv, faux, I_) in [
    (lambda y: y * y, F(0), F(1), F(1, 8), F(5, 16), [(F(5, 8), F(5, 16)), (F(1, 8), F(1, 3)), (F(1, 8), F(3, 8))], F(1, 3)),
    (lambda y: y ** 3, F(0), F(2), F(1), F(7, 2), [(F(9), F(7, 2)), (F(1), F(5)), (F(1), F(4))], F(4)),
    (lambda y: 1 / y, F(1), F(3), F(3, 2), F(16, 15), [(F(5, 6), F(16, 15)), (F(3, 2), F(7, 6)), (F(3, 2), None)], None),
]:
    assert R(f, a_, b_, 2) == Rv and M(f, a_, b_, 2) == Mv  # ex1 (et v1, v2)
    for (rr, mm) in faux:
        assert (rr, mm) != (Rv, Mv)
assert Rd(lambda y: y * y, F(0), F(1), 2) == F(5, 8) and T(lambda y: y * y, F(0), F(1), 2) == F(3, 8)  # leurres ex1
assert F(1, 3) - F(1, 8) == F(5, 24) and F(1, 3) - F(5, 16) == F(1, 48) and F(2, 96) == F(1, 48)  # ex1 correction
assert 4 - 1 == 3 and 4 - F(7, 2) == F(1, 2) and F(8 * 12, 96) == 1 and F(4 * 12, 4) == 12  # ex1 v1
assert Rd(lambda y: y ** 3, F(0), F(2), 2) == 9 and T(lambda y: y ** 3, F(0), F(2), 2) == 5  # ex1 v1 leurres
assert Rd(lambda y: 1 / y, F(1), F(3), 2) == F(5, 6) and T(lambda y: 1 / y, F(1), F(3), 2) == F(7, 6)  # ex1 v2 leurres
assert abs(mp.log(3) - F(16, 15)) > 0.03 and abs(mp.log(3) - mp.mpf(3) / 2 + mp.mpf("0.401")) < 1e-3  # ex1 v2
assert abs(mp.log(3) - mp.mpf(16) / 15 - mp.mpf("0.032")) < 1e-3 and F(8 * 2, 96) == F(1, 6)

# ex2 : point milieu, plus petit n.
def n_milieu(L, M2, eps):
    return plus_petit_n(L ** 3 * M2 / (24 * eps), 2)


assert n_milieu(1, sp.E, sp.Rational(1, 1000)) == 11  # ex2
assert plus_petit_n(sp.E * 1000 / 12, 2) == 16 and plus_petit_n(sp.Rational(1000, 24), 2) == 7  # leurres ex2
assert sp.ceiling(sp.E * 1000 / 24) == 114 and abs(float(sp.E * 1000 / 24) - 113.3) < 0.05
assert n_milieu(sp.pi, 1, sp.Rational(1, 1000)) == 36  # ex2 v1
assert plus_petit_n(sp.pi ** 3 * 1000 / 12, 2) == 51 and plus_petit_n(sp.pi * 1000 / 24, 2) == 12
assert sp.ceiling(sp.pi ** 3 * 1000 / 24) == 1292 and abs(float(sp.pi ** 3 * 1000 / 24) - 1291.9) < 0.05
assert n_milieu(1, 1, sp.Rational(1, 10 ** 4)) == 21  # ex2 v2 (M2 = 1 en t = 1)
assert plus_petit_n(sp.Rational(10 ** 4, 12), 2) == 29 and plus_petit_n(sp.Rational(10 ** 4, 24 * 4), 2) == 11
assert sp.ceiling(sp.Rational(10 ** 4, 24)) == 417
assert sp.diff(sp.ln(t), t, 2) == -1 / t ** 2

# ex3 : formule à un point g(θ) sur [0,1] : exacte pour t ssi θ = 1/2 ; pour t^2 ssi θ = 1/√3.
th = sp.symbols("theta", real=True)
assert sp.solve(sp.Eq(th, sp.Rational(1, 2)), th) == [sp.Rational(1, 2)]
assert sp.solve(sp.Eq(th ** 2, sp.Rational(1, 3)), th) == [-1 / sp.sqrt(3), 1 / sp.sqrt(3)]
assert sp.Rational(1, 2) ** 2 != sp.Rational(1, 3) and 0 != sp.Rational(1, 2)  # leurres θ = 0, θ = 1/√3, tout θ

# ex4 : n (I - R_n) → (b-a)/2 (f(b) - f(a)).
nn = sp.symbols("n", positive=True, integer=True)
Rn_exp = (e - 1) / (nn * (sp.exp(1 / nn) - 1))
assert sp.simplify(sp.summation(sp.exp(sp.Symbol("k") / nn), (sp.Symbol("k"), 0, nn - 1)) / nn - Rn_exp) == 0
lim = sp.limit(nn * ((e - 1) - Rn_exp), nn, sp.oo)
assert sp.simplify(lim - (e - 1) / 2) == 0  # ex4
for faux in [e - 1, -(e - 1) / 2, (e - 1) / 12]:
    assert sp.simplify(lim - faux) != 0
k = sp.Symbol("k", integer=True)
Rn2 = sp.summation(k ** 2, (k, 0, nn - 1)) / nn ** 3  # ex4 v1
assert sp.simplify(Rn2 - (nn - 1) * (2 * nn - 1) / (6 * nn ** 2)) == 0
assert sp.simplify(sp.Rational(1, 3) - Rn2 - (3 * nn - 1) / (6 * nn ** 2)) == 0
assert sp.limit(nn * (sp.Rational(1, 3) - Rn2), nn, sp.oo) == sp.Rational(1, 2)
Rn3 = (2 / nn) ** 4 * sp.summation(k ** 3, (k, 0, nn - 1))  # ex4 v2
assert sp.simplify(Rn3 - 4 * (nn - 1) ** 2 / nn ** 2) == 0
assert sp.simplify(4 - Rn3 - 4 * (2 * nn - 1) / nn ** 2) == 0
assert sp.limit(nn * (4 - Rn3), nn, sp.oo) == 8 and sp.Rational(2 - 0, 2) * (8 - 0) == 8
for n in range(1, 8):  # contrôle direct en fractions
    assert F(1, 3) - R(lambda y: y * y, F(0), F(1), n) == F(3 * n - 1, 6 * n * n)
    assert 4 - R(lambda y: y ** 3, F(0), F(2), n) == F(4 * (2 * n - 1), n * n)

# ex5 : n^2 (I - M_n) → (b-a)^2 (f'(b) - f'(a))/24 ; bornes de Taylor.
assert sp.integrate(s ** 2 / 2, (s, -hh / 2, hh / 2)) == hh ** 3 / 24
assert sp.simplify(2 * sp.integrate(s ** 3 / 6, (s, 0, hh / 2)) - hh ** 4 / 192) == 0
for (f, a_, b_, attendu, faux) in [
    (mp.exp, 0, 1, (mp.e - 1) / 24, [(mp.e - 1) / 12, -(mp.e - 1) / 24, (mp.e - 1) / 2]),
    (mp.log, 1, 2, -mp.mpf(1) / 48, [mp.mpf(1) / 48, -mp.mpf(1) / 24, -mp.mpf(1) / 96]),
    (mp.sin, 0, mp.pi, -mp.pi ** 2 / 12, [mp.pi ** 2 / 12, -mp.pi / 12, -mp.pi ** 2 / 24]),
]:
    a_, b_ = mp.mpf(a_), mp.mpf(b_)
    I_ = mp.quad(f, [a_, b_])
    n = 4000
    val = n ** 2 * (I_ - M(f, a_, b_, n))
    assert abs(val - attendu) < 1e-6  # ex5 (v1, v2)
    for fx in faux:
        assert abs(fx - attendu) > 1e-3
# ex5 : (e-1)/384 ≈ 4,47e-3 et erreur effective ≈ 4,47e-3.
assert abs((mp.e - 1) / 384 - mp.mpf("4.47e-3")) < 1e-5 and abs(I_e - M4 - mp.mpf("4.47e-3")) < 1e-5

# ===========================================================================
# Carte num-l3-integration-trapezes
# ===========================================================================
ff = sp.Function("f")
phi = (t - al) * (t - be) / 2
assert sp.diff(phi, t).expand() == (t - (al + be) / 2).expand()
assert sp.simplify(sp.integrate(phi, (t, al, be)) + (be - al) ** 3 / 12) == 0
# Identité ∫ φ f'' = ∫ f - (β-α)(f(α)+f(β))/2, testée sur une base de polynômes et sur exp.
for g in [t ** j for j in range(7)] + [sp.exp(t), sp.sin(3 * t)]:
    gauche = sp.integrate(phi * sp.diff(g, t, 2), (t, al, be))
    droite = sp.integrate(g, (t, al, be)) - (be - al) * (g.subs(t, al) + g.subs(t, be)) / 2
    assert sp.simplify(gauche - droite) == 0
# T_{2n} = (T_n + M_n)/2 et formule exacte pour x^2 : I - T_n = -(b-a) h^2 / 6.
for n in range(1, 7):
    for f in [lambda y: y * y, lambda y: 1 / (1 + y), lambda y: y ** 5 - 3 * y]:
        assert T(f, F(0), F(1), 2 * n) == (T(f, F(0), F(1), n) + M(f, F(0), F(1), n)) / 2
    assert F(1, 3) - T(lambda y: y * y, F(0), F(1), n) == -F(1, 6 * n * n)
    assert F(1, 3) - T(lambda y: y * y, F(0), F(1), n) == -2 * (F(1, 3) - M(lambda y: y * y, F(0), F(1), n))  # q4, T-ex4 Q4
# exemple : T_4 ≈ 1,727222, I - T_4 ≈ -8,9e-3, e/192 ≈ 1,42e-2, T_8 ≈ 1,720519.
assert abs(T4 - mp.mpf("1.727222")) < 1e-6 and abs(I_e - T4 + mp.mpf("8.9e-3")) < 5e-5
assert abs(mp.e / 192 - mp.mpf("1.42e-2")) < 5e-5 and abs((T4 + M4) / 2 - mp.mpf("1.720519")) < 1e-6
assert abs(T(mp.exp, mp.mpf(0), mp.mpf(1), 8) - (T4 + M4) / 2) < 1e-30
# pièges / contre-exemples : sin(πx), T_1 = 0, I = 2/π ; x^2 sur [0,1/4] : 1/384 > 1/1536 ; x^3 : T_1 = 1/2 ≠ 1/4.
assert T(lambda y: sp.sin(sp.pi * y), sp.Integer(0), sp.Integer(1), 1) == 0
assert sp.integrate(sp.sin(sp.pi * t), (t, 0, 1)) == 2 / sp.pi
assert abs(F(1, 4) ** 3 / 3 - T(lambda y: y * y, F(0), F(1, 4), 1)) == F(1, 384)
assert F(1, 4) ** 3 * F(1, 2) / 12 == F(1, 1536)
assert T(lambda y: y ** 3, F(0), F(1), 1) == F(1, 2) and T(lambda y: y * y, F(0), F(1), 1) == F(1, 2)
assert abs(F(1, 3) - F(1, 2)) == F(1, 6) and abs(F(1, 3) - F(1, 4)) == F(1, 12)
# trous : leurres réfutés (f = 1 ; x^2 sur [0,1] et [0,2]).
assert 2 * 1 != 1 and F(1, 2) * 1 != 1 and 0 != 1  # trou 1 (n = 1, f = 1)
assert F(1, 6) > F(2, 24) and F(1, 24) > F(2, 12 * 8)  # trou 2 : /24, n^3 (n = 2)
assert T(lambda y: y * y, F(0), F(2), 1) - F(8, 3) == F(4, 3) and F(4, 3) > F(4 * 2, 12)  # trou 2 : (b-a)^2
assert T(lambda y: y * y, F(0), F(1), 1) > F(1, 3) > M(lambda y: y * y, F(0), F(1), 1)  # trou 3

# --- exercices de la carte trapèzes ---
assert T(lambda y: y * y, F(0), F(1), 2) == F(3, 8) and F(3, 8) - F(1, 3) == F(1, 24)  # ex1
assert F(2, 12 * 4) == F(1, 24)
assert T(lambda y: y ** 3, F(0), F(2), 2) == 5 and 5 - 4 == 1 and F(8 * 12, 48) == 2  # ex1 v1
assert T(lambda y: 1 / y, F(1), F(3), 2) == F(7, 6) and F(8 * 2, 48) == F(1, 3)  # ex1 v2
assert abs(mp.mpf(7) / 6 - mp.log(3) - mp.mpf("0.068")) < 1e-3
for v in [F(16, 15), F(5, 6)]:
    assert v != F(7, 6)
# ex2 : trapèzes, plus petit n.
f2 = sp.exp(-t ** 2)
assert sp.simplify(sp.diff(f2, t, 2) - (4 * t ** 2 - 2) * sp.exp(-t ** 2)) == 0
assert sp.simplify(sp.diff(f2, t, 3) - 4 * t * (3 - 2 * t ** 2) * sp.exp(-t ** 2)) == 0
assert sp.diff(f2, t, 2).subs(t, 0) == -2 and sp.simplify(sp.diff(f2, t, 2).subs(t, 1) - 2 / e) == 0
assert plus_petit_n(sp.Rational(2 * 10 ** 4, 12), 2) == 41  # ex2
assert plus_petit_n(sp.Rational(2 * 10 ** 4, 24), 2) == 29 and plus_petit_n(2 / e * 10 ** 4 / 12, 2) == 25
assert sp.ceiling(sp.Rational(2 * 10 ** 4, 12)) == 1667
assert plus_petit_n(sp.Rational(10 ** 4, 12), 2) == 29  # ex2 v1
assert plus_petit_n(sp.Rational(10 ** 4, 24), 2) == 21 and plus_petit_n(sp.Rational(10 ** 4, 48), 2) == 15
assert sp.ceiling(sp.Rational(10 ** 4, 12)) == 834
assert plus_petit_n(sp.pi ** 3 * 10 ** 4 / 96, 2) == 57  # ex2 v2
assert plus_petit_n(sp.pi ** 3 * 10 ** 4 / 192, 2) == 41 and plus_petit_n(sp.pi / 2 * 10 ** 4 / 12, 2) == 37
assert sp.ceiling(sp.pi ** 3 * 10 ** 4 / 96) == 3230 and abs(float(sp.pi ** 3 * 10 ** 4 / 96) - 3229.8) < 0.05
# ex3 : T_4 via T_1, M_1, M_2.
for (f, a_, b_, T1, M1, T2, M2, T4v, S2v, I_) in [
    (lambda y: 1 / (1 + y), F(0), F(1), F(3, 4), F(2, 3), F(17, 24), F(24, 35), F(1171, 1680), F(1747, 2520), mp.log(2)),
    (lambda y: 1 / (1 + y * y), F(0), F(1), F(3, 4), F(4, 5), F(31, 40), F(336, 425), F(5323, 6800), F(8011, 10200), mp.pi / 4),
    (lambda y: y ** 4, F(0), F(1), F(1, 2), F(1, 16), F(9, 32), F(41, 256), F(113, 512), F(77, 384), mp.mpf(1) / 5),
]:
    assert T(f, a_, b_, 1) == T1 and M(f, a_, b_, 1) == M1 and T(f, a_, b_, 2) == T2 and M(f, a_, b_, 2) == M2
    assert T(f, a_, b_, 4) == T4v == (T2 + M2) / 2 and S(f, a_, b_, 2) == S2v  # ex3 (v1, v2)
    assert len({T4v, S2v, T2, M2}) == 4
    ratio = (mp.mpf(T2.numerator) / T2.denominator - I_) / (mp.mpf(T4v.numerator) / T4v.denominator - I_)
    assert 3.8 < ratio < 4.2  # « divisée par environ 4 »
assert abs(mp.mpf(17) / 24 - mp.log(2) - mp.mpf("1.52e-2")) < 5e-5 and abs(mp.mpf(1171) / 1680 - mp.log(2) - mp.mpf("3.88e-3")) < 5e-6
assert abs(mp.pi / 4 - mp.mpf(31) / 40 - mp.mpf("1.04e-2")) < 5e-5 and abs(mp.pi / 4 - mp.mpf(5323) / 6800 - mp.mpf("2.60e-3")) < 5e-6
assert F(9, 32) - F(1, 5) == F(13, 160) and F(113, 512) - F(1, 5) == F(53, 2560)
# ex4 : limites -2, sommes ; contrôle numérique sur exp.
for n in [200, 400]:
    q = (I_e - T(mp.exp, mp.mpf(0), mp.mpf(1), n)) / (I_e - M(mp.exp, mp.mpf(0), mp.mpf(1), n))
    assert abs(q + 2) < 1e-4
assert F(-1, 12) / F(1, 24) == -2 and -F(1, 12) + 2 * F(1, 24) == 0
# ex5 : trapèzes et fonctions périodiques.
for n in range(1, 9):
    for kk in range(-20, 21):
        val = 2 * mp.pi / n * sum(mp.expj(kk * 2 * mp.pi * j / n) for j in range(n))
        attendu = 2 * mp.pi if kk % n == 0 else 0
        assert abs(val - attendu) < 1e-30
rr = 2 - mp.sqrt(3)
cas_per = [
    (lambda y: 1 / (2 - mp.cos(y)), 1 / mp.sqrt(3), rr, 2 * mp.pi / mp.sqrt(3),
     lambda n: 4 * mp.pi / mp.sqrt(3) * rr ** n / (1 - rr ** n),
     [lambda n: 4 * mp.pi / mp.sqrt(3) * rr ** n, lambda n: 2 * mp.pi / mp.sqrt(3) * rr ** n / (1 - rr ** n),
      lambda n: 4 * mp.pi / mp.sqrt(3) * rr ** (2 * n) / (1 - rr ** (2 * n))]),
    (lambda y: 1 / (5 - 4 * mp.cos(y)), mp.mpf(1) / 3, mp.mpf(1) / 2, 2 * mp.pi / 3,
     lambda n: 4 * mp.pi / (3 * (2 ** n - 1)),
     [lambda n: 4 * mp.pi / (3 * 2 ** n), lambda n: 2 * mp.pi / (3 * (2 ** n - 1)), lambda n: 4 * mp.pi / (3 * (4 ** n - 1))]),
    (lambda y: 1 / (5 + 3 * mp.cos(y)), mp.mpf(1) / 4, -mp.mpf(1) / 3, mp.pi / 2,
     lambda n: (-1) ** n * mp.pi / (3 ** n - (-1) ** n),
     [lambda n: mp.pi / (3 ** n + 1), lambda n: (-1) ** n * mp.pi / 3 ** n,
      lambda n: (-1) ** n * mp.pi / (2 * (3 ** n - (-1) ** n))]),
]
for (f, cst, q_, I_, erreur, faux) in cas_per:  # ex5 (v1, v2)
    # développement admis : f(x) = cst (1 + 2 Σ q^k cos kx)
    for y in [mp.mpf("0.3"), mp.mpf("1.7"), mp.mpf("4.1")]:
        serie = cst * (1 + 2 * mp.nsum(lambda kk: q_ ** kk * mp.cos(kk * y), [1, mp.inf]))
        assert abs(serie - f(y)) < 1e-25
    assert abs(mp.quad(f, [0, mp.pi, 2 * mp.pi]) - I_) < 1e-25 and abs(2 * mp.pi * cst - I_) < 1e-30
    for n in range(1, 13):
        Tn = 2 * mp.pi / n * sum(f(2 * mp.pi * j / n) for j in range(n))
        assert abs(T(f, mp.mpf(0), 2 * mp.pi, n) - Tn) < 1e-30
        assert abs(Tn - I_ - erreur(n)) < 1e-25
        for g in faux:  # leurres : jamais égaux à la bonne réponse (et nettement différents pour n = 1)
            assert abs(g(n) - erreur(n)) > mp.mpf(10) ** -35 * abs(erreur(n))
            if n == 1:
                assert abs(g(n) - erreur(n)) > 1e-2 * abs(erreur(n))
assert abs(4 * mp.pi / mp.sqrt(3) * rr ** 10 / (1 - rr ** 10) - mp.mpf("1.4e-5")) < 5e-7  # ex5 Q3
assert abs(1 / rr - mp.mpf("3.73")) < 0.01 and abs(rr - mp.mpf("0.268")) < 1e-3
assert abs(4 * mp.pi / (3 * (2 ** 10 - 1)) - mp.mpf("4.1e-3")) < 5e-5 and abs(4 * mp.pi / (3 * (2 ** 30 - 1)) - mp.mpf("3.9e-9")) < 5e-11
assert abs(mp.pi / (3 ** 10 - 1) - mp.mpf("5.3e-5")) < 5e-7
assert abs(2 * mp.pi * (1 - 1 / mp.sqrt(3)) - 4 * mp.pi / mp.sqrt(3) * rr / (1 - rr)) < 1e-30  # n = 1
# affirmation (carte erreur) : T_n(cos) = 0 sur [0, 2π] pour n >= 2.
for n in range(2, 9):
    assert abs(T(mp.cos, mp.mpf(0), 2 * mp.pi, n)) < 1e-30

# ===========================================================================
# Carte num-l3-integration-simpson
# ===========================================================================
mm = sp.symbols("m", real=True)
U = sp.symbols("U", positive=True)
Fu = sp.integrate(ff(t), (t, mm - U, mm + U)) - U / 3 * (ff(mm - U) + 4 * ff(mm) + ff(mm + U))
F1 = sp.diff(Fu, U)
F2 = sp.diff(Fu, U, 2)
F3 = sp.diff(Fu, U, 3)
fp = lambda y, kk: sp.diff(ff(t), t, kk).subs(t, y)  # noqa: E731
attF1 = sp.Rational(2, 3) * (ff(mm + U) + ff(mm - U)) - sp.Rational(4, 3) * ff(mm) - U / 3 * (fp(mm + U, 1) - fp(mm - U, 1))
attF2 = (fp(mm + U, 1) - fp(mm - U, 1)) / 3 - U / 3 * (fp(mm + U, 2) + fp(mm - U, 2))
attF3 = -U / 3 * (fp(mm + U, 3) - fp(mm - U, 3))
for (calc, att) in [(F1, attF1), (F2, attF2), (F3, attF3)]:
    for g in [sp.exp(2 * t), sp.sin(t) * t, t ** 7 - t ** 3]:
        # substitution explicite de f par g
        rep = {ff: sp.Lambda(t, g)}
        assert sp.simplify((calc.subs(rep) - att.subs(rep)).doit()) == 0
assert sp.integrate((r - u) ** 2 * u ** 2, (u, 0, r)) == r ** 5 / 30
assert sp.simplify(((h / 2) ** 5) / 90 - h ** 5 / 2880) == 0
assert sp.integrate((r - u) ** 2 / 2 * (-sp.Rational(2, 3) * u ** 2), (u, 0, r)) == -r ** 5 / 90
# F(r) exact pour f = x^4 : erreur -(2r)^5/2880 * 24.
F4 = sp.integrate(t ** 4, (t, mm - U, mm + U)) - U / 3 * ((mm - U) ** 4 + 4 * mm ** 4 + (mm + U) ** 4)
assert sp.simplify(F4.subs(U, r) + (2 * r) ** 5 / 2880 * 24) == 0
# Exactitude : degré 3, pas 4 (sur [-r, r] avec la formule r/3 (g(-r) + 4 g(0) + g(r))).
for j in range(5):  # (les degrés impairs sont toujours exacts par symétrie ; le premier échec est j = 4)
    exact = sp.integrate(s ** j, (s, -r, r))
    formule = r / 3 * ((-r) ** j + 4 * (0 if j else 1) + r ** j)
    assert (sp.simplify(exact - formule) == 0) == (j <= 3)
# S_n = (T_n + 2 M_n)/3 = (4 T_{2n} - T_n)/3 ; borne et formule exacte pour x^4.
for n in range(1, 6):
    for f in [lambda y: 1 / (1 + y), lambda y: y ** 4, lambda y: y ** 6 - y]:
        h_ = F(1, n)
        direct = h_ / 6 * sum(f(k * h_) + 4 * f(k * h_ + h_ / 2) + f((k + 1) * h_) for k in range(n))
        assert direct == S(f, F(0), F(1), n) == (4 * T(f, F(0), F(1), 2 * n) - T(f, F(0), F(1), n)) / 3
    assert F(1, 5) - S(lambda y: y ** 4, F(0), F(1), n) == -F(24, 2880 * n ** 4)
for n in range(1, 4):  # polynômes de degré 3 : exacte
    assert S(lambda y: 2 * y ** 3 - y + 5, F(-1), F(2), n) == F(2 * 16, 4) - F(2 * 1, 4) - (F(4, 2) - F(1, 2)) + 15
# exemples
S2e = S(mp.exp, mp.mpf(0), mp.mpf(1), 2)
assert abs(S2e - mp.mpf("1.718319")) < 1e-6 and abs(I_e - S2e + mp.mpf("3.7e-5")) < 5e-7
assert abs(mp.e / (2880 * 16) - mp.mpf("5.9e-5")) < 5e-7
assert S(lambda y: 1 / (1 + y), F(0), F(1), 2) == F(1747, 2520)
assert abs(mp.mpf(1747) / 2520 - mp.mpf("0.693254")) < 1e-6
# contre-exemples : x^4, poids 1/4-1/2-1/4 = T_{2n} (pas exact pour x^2), max|f''| (1/120 > 1/240), x^6 (23/672 > 1/32).
assert S(lambda y: y ** 4, F(0), F(1), 1) == F(5, 24)
assert (T(lambda y: y * y, F(0), F(1), 1) + M(lambda y: y * y, F(0), F(1), 1)) / 2 == F(3, 8) != F(1, 3)
assert F(12, 2880) == F(1, 240) and F(1, 120) > F(1, 240)
assert F(1, 7) - S(lambda y: y ** 6, F(0), F(1), 1) == -F(23, 672) and F(23, 672) > F(1, 32)
assert sp.diff(t ** 6, t, 4).subs(t, sp.Rational(1, 2)) == 90 and F(90, 2880) == F(1, 32)
g5 = lambda y: 5 * y ** 4 - 9 * y ** 2  # noqa: E731
assert T(g5, F(0), F(1), 1) == -2 == F(1) - F(3) and S(g5, F(0), F(1), 1) == F(-47, 24)
# trous : leurres réfutés.
T1q, M1q = T(lambda y: y * y, F(0), F(1), 1), M(lambda y: y * y, F(0), F(1), 1)
for comb_ in [(2 * T1q + M1q) / 3, (T1q + M1q) / 2, (4 * T1q - M1q) / 3]:
    assert comb_ != F(1, 3)  # trou 1 (x^2, n = 1)
assert S(lambda y: y ** 3, F(0), F(1), 1) == F(1, 4) and S(lambda y: y ** 5, F(0), F(1), 1) != F(1, 6)  # trou 2
assert abs(F(1, 5) - S(lambda y: y ** 4, F(0), F(1), 2)) == F(1, 1920) and F(1, 1920) > F(24, 2880 * 32)  # trou 3 : n^5
assert abs(F(32, 5) - S(lambda y: y ** 4, F(0), F(2), 1)) == F(4, 15) and F(4, 15) > F(16 * 24, 2880)  # trou 3 : (b-a)^4
assert F(24, 5760) == F(1, 240) < F(1, 120)  # trou 3 : 5760

# --- exercices de la carte Simpson ---
for (f, a_, b_, I_, S1v, err, faux) in [
    (lambda y: y ** 4, F(0), F(1), F(1, 5), F(5, 24), F(-1, 120), [F(1, 120), F(0), F(-1, 2880)]),
    (lambda y: y ** 4, F(0), F(2), F(32, 5), F(20, 3), F(-4, 15), [F(4, 15), F(0), F(-1, 120)]),
    (lambda y: y ** 5, F(0), F(1), F(1, 6), F(3, 16), F(-1, 48), [F(1, 48), F(0), F(-1, 120)]),
]:
    assert S(f, a_, b_, 1) == S1v and I_ - S1v == err  # ex1 (v1, v2)
    assert err not in faux and len(set(faux)) == 3
assert -F(2 * 16 * 24, 2880) == F(-4, 15) and -F(120, 2880) * F(1, 2) == F(-1, 48)  # ξ = 1/2 pour t^5
# ex2 : Simpson, plus petit n.
assert plus_petit_n(sp.E * 10 ** 8 / 2880, 4) == 18 and 17 ** 4 == 83521 and 18 ** 4 == 104976  # ex2
assert abs(float(sp.E * 10 ** 8 / 2880) - 94384) < 1
assert plus_petit_n(sp.Integer(10 ** 8) / 2880, 4) == 14 and plus_petit_n(sp.E * 10 ** 8 / 2880, 2) == 308
assert plus_petit_n(sp.E * 10 ** 8 / 12, 2) == 4760
assert plus_petit_n(sp.Rational(24 * 10 ** 8, 2880), 4) == 31 and 30 ** 4 == 810000 and 31 ** 4 == 923521  # ex2 v1
assert sp.diff(1 / t, t, 4) == 24 / t ** 5
assert plus_petit_n(sp.Rational(3, 4) * 10 ** 8 / 2880, 4) == 13 and plus_petit_n(sp.Rational(24 * 10 ** 8, 2880), 2) == 913
assert plus_petit_n(sp.Rational(2 * 10 ** 8, 12), 2) == 4083
assert plus_petit_n(sp.pi ** 5 * 10 ** 8 / 2880, 4) == 58 and 57 ** 4 == 10556001 and 58 ** 4 == 11316496  # ex2 v2
assert abs(float(sp.pi ** 5 * 10 ** 8 / 2880) - 10625684) < 1
assert plus_petit_n(sp.pi * 10 ** 8 / 2880, 4) == 19 and plus_petit_n(sp.pi ** 5 * 10 ** 8 / 2880, 2) == 3260
assert plus_petit_n(sp.pi ** 3 * 10 ** 8 / 12, 2) == 16075
# ex3 : poids (1/3, 4/3, 1/3), degré 3.
l0, l1, l2 = sp.symbols("l0 l1 l2")
sol = sp.solve([l0 + l1 + l2 - 2, -l0 + l2, l0 + l2 - sp.Rational(2, 3)], [l0, l1, l2])
assert sol == {l0: sp.Rational(1, 3), l1: sp.Rational(4, 3), l2: sp.Rational(1, 3)}
assert sp.Rational(1, 3) * (-1) ** 3 + sp.Rational(1, 3) == 0 and sp.integrate(s ** 4, (s, -1, 1)) == sp.Rational(2, 5)
assert sp.Rational(2, 3) != sp.Rational(2, 5)
assert sp.integrate(s ** 2, (s, -1, 1)) != sp.Rational(1, 2) + sp.Rational(1, 2)  # leurre (1/2, 1, 1/2) : degré 1
# ex4 : S_2 = (T_2 + 2 M_2)/3.
for (f, a_, b_, T2, M2, S2v, faux, I_) in [
    (lambda y: 4 / (1 + y * y), F(0), F(1), F(31, 10), F(1344, 425), F(8011, 2550), [F(3979, 1275), F(5323, 1700), F(8117, 2550)], mp.pi),
    (lambda y: 1 / (1 + y), F(0), F(1), F(17, 24), F(24, 35), F(1747, 2520), [F(883, 1260), F(1171, 1680), F(1709, 2520)], mp.log(2)),
    (lambda y: 1 / (y * y), F(1), F(3), F(29, 36), F(136, 225), F(1813, 2700), [F(997, 1350), F(141, 200), F(1709 * 0 + 1451, 2700)], mp.mpf(2) / 3),
]:
    assert T(f, a_, b_, 2) == T2 and M(f, a_, b_, 2) == M2 and S(f, a_, b_, 2) == S2v  # ex4 (v1, v2)
    assert faux == [(2 * T2 + M2) / 3, (T2 + M2) / 2, (4 * M2 - T2) / 3] and S2v not in faux
    eT, eM = I_ - mp.mpf(T2.numerator) / T2.denominator, I_ - mp.mpf(M2.numerator) / M2.denominator
    assert eT * eM < 0
assert [4 / (1 + F(j, 4) ** 2) for j in range(5)] == [4, F(64, 17), F(16, 5), F(64, 25), 2]
assert abs(mp.pi - F(31, 10) - mp.mpf("4.16e-2")) < 5e-5 and abs(mp.pi - mp.mpf(1344) / 425 + mp.mpf("2.08e-2")) < 5e-5
assert abs(mp.pi - mp.mpf(8011) / 2550 - mp.mpf("2.4e-5")) < 5e-7 and abs(mp.mpf(8011) / 2550 - mp.mpf("3.1415686")) < 1e-7
assert abs((mp.pi - mp.mpf(31) / 10) / (mp.pi - mp.mpf(1344) / 425) + 2) < 0.01
assert abs(mp.log(2) - mp.mpf(17) / 24 + mp.mpf("1.52e-2")) < 5e-5 and abs(mp.log(2) - mp.mpf(24) / 35 - mp.mpf("7.4e-3")) < 5e-5
assert abs(mp.log(2) - mp.mpf(1747) / 2520 + mp.mpf("1.07e-4")) < 5e-7
assert abs(mp.mpf(2) / 3 - mp.mpf(29) / 36 + mp.mpf("0.139")) < 5e-4 and abs(mp.mpf(2) / 3 - mp.mpf(136) / 225 - mp.mpf("0.062")) < 5e-4
assert abs(mp.mpf(2) / 3 - mp.mpf(1813) / 2700 + mp.mpf("4.8e-3")) < 5e-5
assert F(29, 36) == F(1, 2) + F(1, 4) + F(1, 18) and F(136, 225) == F(4, 9) + F(4, 25)
# ex5 : n^4 (I - S_n) → -(b-a)^4 (f'''(b) - f'''(a))/2880.
for (f, a_, b_, attendu, faux) in [
    (mp.exp, 0, 1, -(mp.e - 1) / 2880, [(mp.e - 1) / 2880, -(mp.e - 1) / 180, -mp.e / 2880]),
    (mp.log, 1, 2, mp.mpf(7) / 11520, [-mp.mpf(7) / 11520, mp.mpf(7) / 720, mp.mpf(7) / 5760]),
    (mp.sin, 0, mp.pi, -mp.pi ** 4 / 1440, [mp.pi ** 4 / 1440, -mp.pi ** 4 / 90, -mp.pi ** 4 / 2880]),
]:
    a_, b_ = mp.mpf(a_), mp.mpf(b_)
    I_ = mp.quad(f, [a_, b_])
    n = 400
    assert abs(n ** 4 * (I_ - S(f, a_, b_, n)) - attendu) < 1e-6 * max(1, abs(attendu)) + 1e-8  # ex5 (v1, v2)
    for fx in faux:
        assert abs(fx - attendu) > 1e-5
assert sp.diff(sp.ln(t), t, 3) == 2 / t ** 3 and sp.Rational(2, 8) - 2 == -sp.Rational(7, 4)
assert -(-sp.Rational(7, 4)) / 2880 == sp.Rational(7, 11520) and sp.diff(sp.ln(t), t, 4) == -6 / t ** 4
assert abs((mp.e - 1) / 2880 - mp.mpf("5.97e-4")) < 5e-7 and abs((mp.e - 1) / (2880 * 16) - mp.mpf("3.73e-5")) < 5e-8
assert abs(I_e - S2e + mp.mpf("3.70e-5")) < 5e-8 and abs(mp.mpf(7) / 11520 - mp.mpf("6.08e-4")) < 5e-7
assert abs(mp.pi ** 4 / 1440 - mp.mpf("0.0676")) < 5e-5

# ===========================================================================
# Carte num-l3-integration-erreur
# ===========================================================================
dd = sp.symbols("d", integer=True, nonnegative=True)
assert sp.simplify(sp.integrate(s ** (dd + 1) / sp.factorial(dd + 1), (s, 0, h)) - h ** (dd + 2) / sp.factorial(dd + 2)) == 0
# Majoration générale : vérifiée numériquement pour le point milieu (d=1, W=1), Simpson (d=3, W=1) et Radau (d=2, W=1).
for n in range(1, 9):
    hn = mp.mpf(1) / n
    assert abs(I_e - M(mp.exp, mp.mpf(0), mp.mpf(1), n)) <= (1 + mp.mpf(1) / 3) / 2 * mp.e * hn ** 2
    assert abs(I_e - S(mp.exp, mp.mpf(0), mp.mpf(1), n)) <= (1 + mp.mpf(1) / 5) / 24 * mp.e * hn ** 4
# degrés : rectangles 0, point milieu 1, trapèzes 1, Simpson 3.
def degre(formule):
    j = 0
    while formule(lambda y: y ** j) == F(1, j + 1):
        j += 1
    return j - 1


assert degre(lambda g: g(F(0))) == 0 and degre(lambda g: g(F(1, 2))) == 1
assert degre(lambda g: (g(F(0)) + g(F(1))) / 2) == 1 and degre(lambda g: (g(F(0)) + 4 * g(F(1, 2)) + g(F(1))) / 6) == 3
# encadrement convexe et exemple : 24/35 <= ln 2 <= 17/24, (T_2 - M_2)/2 = 19/1680, |ln 2 - T_4| ≈ 0,0039.
assert fl(F(24, 35)) < mp.log(2) < fl(F(17, 24)) and (F(17, 24) - F(24, 35)) / 2 == F(19, 1680)
assert abs(mp.mpf(1171) / 1680 - mp.log(2) - mp.mpf("0.0039")) < 5e-5 and abs(mp.mpf(19) / 1680 - mp.mpf("0.0113")) < 5e-5
# contre-exemples : √t sur [0,1] : T_1 = 1/2 < 2/3 < M_1 = √2/2 ; premier panneau de Simpson.
assert sp.Rational(1, 2) < sp.Rational(2, 3) < sp.sqrt(2) / 2
assert sp.simplify(sp.Rational(2, 3) - (4 * sp.sqrt(sp.Rational(1, 2)) + 1) / 6 - (3 - 2 * sp.sqrt(2)) / 6) == 0
# Gauss à deux points : degré 3 (exemple) ; borne 2r+1 (affirmation).
gs = 1 / sp.sqrt(3)
for j in range(5):  # premier échec en j = 4 (les degrés impairs sont exacts par symétrie)
    assert (sp.simplify(gs ** j + (-gs) ** j - sp.integrate(s ** j, (s, -1, 1))) == 0) == (j <= 3)
# trous : O(h^{d+2}) faux pour le point milieu (erreur ~ (e-1)h^2/24) ; O(h^{2d+2}) faux pour les rectangles (d=0) ;
# R_n, M_n (x^2) ; T_n, S_n (x^4) ; R_n, T_n ((1-x)^2).
assert abs(1000 ** 2 * (I_e - M(mp.exp, mp.mpf(0), mp.mpf(1), 1000)) - (mp.e - 1) / 24) < 1e-6
assert abs(1000 * (I_e - R(mp.exp, mp.mpf(0), mp.mpf(1), 1000)) - (mp.e - 1) / 2) < 1e-3
assert not (R(lambda y: y * y, F(0), F(1), 1) <= F(1, 3) <= M(lambda y: y * y, F(0), F(1), 1))
assert not (S(lambda y: y ** 4, F(0), F(1), 1) <= F(1, 5) <= T(lambda y: y ** 4, F(0), F(1), 1))
assert R(lambda y: (1 - y) ** 2, F(0), F(1), 1) == 1 and T(lambda y: (1 - y) ** 2, F(0), F(1), 1) == F(1, 2)
assert not (min(F(1), F(1, 2)) <= F(1, 3) <= max(F(1), F(1, 2)))
# pièges : d-1 (rectangles avec d = 1 : majoration 0 mais erreur 1/(2n)) ; |f''(a)| (x^3, M_1).
assert F(1, 2) - R(lambda y: y, F(0), F(1), 1) == F(1, 2)

# --- exercices de la carte erreur ---
# ex1 : encadrements.
assert M(lambda y: 1 / (1 + y), F(0), F(1), 2) == F(24, 35) and T(lambda y: 1 / (1 + y), F(0), F(1), 2) == F(17, 24)  # ex1
assert not (mp.log(2) <= fl(F(24, 35))) and not (mp.log(2) >= fl(F(17, 24)))
assert M(lambda y: 1 / y, F(1), F(3), 2) == F(16, 15) and T(lambda y: 1 / y, F(1), F(3), 2) == F(7, 6)  # ex1 v1
assert fl(F(16, 15)) < mp.log(3) < fl(F(7, 6)) and (F(7, 6) - F(16, 15)) / 2 == F(1, 20) and mp.log(3) > 1
assert abs(fl(T(lambda y: 1 / y, F(1), F(3), 4)) - mp.log(3)) < fl(F(1, 20)) and abs(mp.mpf(67) / 60 - mp.log(3) - mp.mpf("0.018")) < 5e-4
assert M(lambda y: 1 / y, F(1), F(4), 3) == F(142, 105) and T(lambda y: 1 / y, F(1), F(4), 3) == F(35, 24)  # ex1 v2
assert fl(F(142, 105)) < mp.log(4) < fl(F(35, 24)) and (F(35, 24) - F(142, 105)) / 2 == F(89, 1680)
assert fl(F(4, 3)) < mp.log(4)
T6 = T(lambda y: 1 / y, F(1), F(4), 6)
assert T6 == (F(35, 24) + F(142, 105)) / 2 and abs(fl(T6) - mp.log(4)) <= fl(F(89, 1680))
assert abs(mp.mpf(T6.numerator) / T6.denominator - mp.log(4) - mp.mpf("0.019")) < 5e-4
assert abs(mp.mpf(89) / 1680 - mp.mpf("0.053")) < 5e-4
# ex2 : degrés d'exactitude.
assert degre(lambda g: g(F(0)) / 4 + 3 * g(F(2, 3)) / 4) == 2 and F(3, 4) * F(8, 27) == F(2, 9)  # ex2
assert (1 + F(1, 4)) / 6 == F(5, 24)
uu = sp.sqrt(3) / 6
noeuds = [sp.Rational(1, 2) - uu, sp.Rational(1, 2) + uu]
moments = [sp.simplify(sum(sp.Rational(1, 2) * y ** j for y in noeuds)) for j in range(6)]
assert [mm_ == sp.Rational(1, j + 1) for j, mm_ in enumerate(moments)] == [True] * 4 + [False] * 2  # ex2 v1 : d = 3
assert moments[4] == sp.Rational(7, 36) and uu ** 2 == sp.Rational(1, 12) and (1 + F(1, 5)) / 24 == F(1, 20)
assert degre(lambda g: (g(F(1, 4)) + g(F(3, 4))) / 2) == 1 and (1 + F(1, 3)) / 2 == F(2, 3)  # ex2 v2
for n in range(1, 5):  # v2 = point milieu à 2n sous-intervalles ; constante optimale (b-a)^3 M2 / (24 (2n)^2) = (b-a) M2 h^2 / 96
    f = lambda y: 1 / (1 + y)  # noqa: E731
    h_ = F(1, n)
    assert sum(h_ * (f(k * h_ + h_ / 4) + f(k * h_ + 3 * h_ / 4)) / 2 for k in range(n)) == M(f, F(0), F(1), 2 * n)
# ex3 : ordre observé ; les valeurs de l'énoncé sont M_n, R_n, S_n arrondies à 12 décimales.
vals = {"M": ["1.713815279771", "1.717163664996", "1.718002192053"],
        "R": ["1.512436676000", "1.613125977886", "1.665144821441"],
        "S": ["1.718284154700", "1.718281974052", "1.718281837562"]}
for cle, meth, p, rapport_txt in [("M", M, 2, "3.99"), ("R", R, 1, "1.94"), ("S", S, 4, "15.98")]:
    exacts = [meth(mp.exp, mp.mpf(0), mp.mpf(1), n) for n in (4, 8, 16)]
    for v, ex_ in zip(vals[cle], exacts):
        assert abs(mp.mpf(v) - ex_) <= mp.mpf("5e-13")  # ex3 (v1 : R, v2 : S)
    a4, a8, a16 = (mp.mpf(v) for v in vals[cle])
    rapport = (a4 - a8) / (a8 - a16)
    assert abs(rapport - mp.mpf(rapport_txt)) < 0.005
    puissances = {q: abs(rapport - 2 ** q) for q in (1, 2, 3, 4)}
    assert min(puissances, key=puissances.get) == p
assert abs(mp.mpf(vals["M"][0]) - mp.mpf(vals["M"][1]) + mp.mpf("3.348e-3")) < 5e-7
assert abs(mp.mpf(vals["M"][1]) - mp.mpf(vals["M"][2]) + mp.mpf("8.385e-4")) < 5e-8
assert abs(mp.mpf(vals["R"][0]) - mp.mpf(vals["R"][1]) + mp.mpf("0.1007")) < 5e-5
assert abs(mp.mpf(vals["R"][1]) - mp.mpf(vals["R"][2]) + mp.mpf("0.0520")) < 5e-5
assert abs(mp.mpf(vals["S"][0]) - mp.mpf(vals["S"][1]) - mp.mpf("2.1806e-6")) < 5e-11
assert abs(mp.mpf(vals["S"][1]) - mp.mpf(vals["S"][2]) - mp.mpf("1.3649e-7")) < 5e-12
# ex4 : √t ; Simpson, point milieu (v1), trapèzes (v2).
hs = sp.symbols("hs", positive=True)
sq = sp.sqrt
I0 = sp.integrate(sq(t), (t, 0, hs))
assert sp.simplify(I0 - sp.Rational(2, 3) * hs ** sp.Rational(3, 2)) == 0
cS = sp.simplify((I0 - hs / 6 * (0 + 4 * sq(hs / 2) + sq(hs))) / hs ** sp.Rational(3, 2))
cM = sp.simplify((hs * sq(hs / 2) - I0) / hs ** sp.Rational(3, 2))
cT = sp.simplify((I0 - hs / 2 * (0 + sq(hs))) / hs ** sp.Rational(3, 2))
assert sp.simplify(cS - (3 - 2 * sp.sqrt(2)) / 6) == 0 and abs(float(cS) - 0.0286) < 1e-4  # ex4
assert sp.simplify(cM - (sp.sqrt(2) / 2 - sp.Rational(2, 3))) == 0 and abs(float(cM) - 0.0404) < 1e-4  # ex4 v1
assert cT == sp.Rational(1, 6)  # ex4 v2
for c_, faux in [(cS, [(2 * sp.sqrt(2) - 3) / 6, (3 - 2 * sp.sqrt(2)) / 3, (1 + 2 * sp.sqrt(2)) / 6]),
                 (cM, [sp.Rational(2, 3) - sp.sqrt(2) / 2, sp.sqrt(2) / 2, sp.Rational(1, 24)]),
                 (cT, [-sp.Rational(1, 6), sp.Rational(1, 12), sp.Rational(2, 3)])]:
    for fx in faux:
        assert sp.simplify(c_ - fx) != 0
assert sp.diff(sq(t), t, 4) == -sp.Rational(15, 16) * t ** sp.Rational(-7, 2)
assert sp.diff(sq(t), t, 2) == -sp.Rational(1, 4) * t ** sp.Rational(-3, 2)
assert sp.Rational(1, 24) * sp.Rational(1, 4) == sp.Rational(1, 96) and sp.Rational(1, 12) * sp.Rational(1, 4) == sp.Rational(1, 48)
for n in [1, 2, 5, 10, 40, 160]:  # encadrements numériques c h^{3/2} <= erreur <= C h^{3/2}
    hn = mp.mpf(1) / n
    sqrtf = lambda y: mp.sqrt(y)  # noqa: E731
    eS = mp.mpf(2) / 3 - S(sqrtf, mp.mpf(0), mp.mpf(1), n)
    eM = M(sqrtf, mp.mpf(0), mp.mpf(1), n) - mp.mpf(2) / 3
    eT = mp.mpf(2) / 3 - T(sqrtf, mp.mpf(0), mp.mpf(1), n)
    cSm, cMm = (3 - 2 * mp.sqrt(2)) / 6, mp.sqrt(2) / 2 - mp.mpf(2) / 3
    tol = mp.mpf(10) ** -30  # égalité pour n = 1 (un seul panneau)
    assert cSm * hn ** 1.5 - tol <= eS <= (cSm + mp.mpf(15) / (16 * 2880) * mp.zeta(3.5)) * hn ** 1.5
    assert cMm * hn ** 1.5 - tol <= eM <= (cMm + mp.zeta(1.5) / 96) * hn ** 1.5
    assert mp.mpf(1) / 6 * hn ** 1.5 - tol <= eT <= (mp.mpf(1) / 6 + mp.zeta(1.5) / 48) * hn ** 1.5
# ex5 : Gauss à deux points.
w_, t_ = sp.symbols("w t_", positive=True)
assert sp.solve([2 * w_ - 2, 2 * w_ * t_ ** 2 - sp.Rational(2, 3)], [w_, t_]) == [(1, 1 / sp.sqrt(3))]
assert 2 * (1 / sp.sqrt(3)) ** 4 == sp.Rational(2, 9) and sp.integrate(s ** 4, (s, -1, 1)) == sp.Rational(2, 5)
assert 2 * sp.Rational(1, 2) ** 2 != sp.Rational(2, 3) and 2 * (1 / sp.sqrt(2)) ** 2 != sp.Rational(2, 3)  # leurres t = 1/2, 1/√2
s0, s1 = sp.symbols("s0 s1", real=True)
Pq = (s - s0) ** 2 * (s - s1) ** 2
assert sp.expand(Pq.subs(s, s0)) == 0 and sp.expand(Pq.subs(s, s1)) == 0
assert (1 + F(1, 5)) / 24 == F(1, 20)

# ===========================================================================
# Carte num-l3-integration-romberg
# ===========================================================================
# Richardson : (2^p A(h/2) - A(h))/(2^p - 1) élimine c h^p.
A_, cc, p_ = sp.symbols("A c p", positive=True)
assert sp.simplify((2 ** p_ * (A_ + cc * (h / 2) ** p_) - (A_ + cc * h ** p_)) / (2 ** p_ - 1) - A_) == 0
# Leurres du trou 1 et pièges : A(h) = A + h (p = 1) ; A(h) = A + h^2 (pas h/3).
Ah = lambda hv: A_ + hv  # noqa: E731
assert sp.simplify((2 * Ah(h) - Ah(h / 2)) / 1 - A_) == sp.Rational(3, 2) * h
assert sp.simplify((4 * Ah(h / 2) - Ah(h)) / 3 - A_) == h / 3
assert sp.simplify((A_ + h ** 2 + A_ + (h / 2) ** 2) / 2 - A_) == sp.Rational(5, 8) * h ** 2
assert sp.simplify((4 * (A_ + (h / 3) ** 2) - (A_ + h ** 2)) / 3 - A_) == -sp.Rational(5, 27) * h ** 2
# Développement des trapèzes : constantes 1/60, 1/12, 1/10, 1/320, 1/288, 19/2880.
assert sp.integrate(s ** 4 / 24, (s, -r, r)) == r ** 5 / 60 and sp.integrate(s ** 2 / 2, (s, -r, r)) == r ** 3 / 3
assert r * 2 * r ** 4 / 24 == r ** 5 / 12 and sp.Rational(1, 60) + sp.Rational(1, 12) == sp.Rational(1, 10)
assert sp.Rational(1, 3) - 1 == -sp.Rational(2, 3) and sp.Rational(2, 3) * (h / 2) ** 3 == h ** 3 / 12
assert sp.Rational(1, 10) / 32 == sp.Rational(1, 320) and sp.Rational(1, 12) / 24 == sp.Rational(1, 288)
assert sp.Rational(1, 288) + sp.Rational(1, 320) == sp.Rational(19, 2880)
for (f, df, d4max, a_, b_) in [(mp.exp, mp.exp, mp.e, 0, 1), (mp.sin, mp.cos, 1, 0, 3),
                               (lambda y: 1 / (1 + y), lambda y: -1 / (1 + y) ** 2, 24, 0, 1)]:
    a_, b_ = mp.mpf(a_), mp.mpf(b_)
    I_ = mp.quad(f, [a_, b_])
    for n in range(1, 12):
        hn = (b_ - a_) / n
        reste = T(f, a_, b_, n) - I_ - hn ** 2 / 12 * (df(b_) - df(a_))
        assert abs(reste) <= mp.mpf(19) / 2880 * (b_ - a_) * d4max * hn ** 4
# Euler-Maclaurin (admis) : coefficients 1/12 et -1/720 (nombres de Bernoulli), vérifiés sur exp.
assert sp.bernoulli(2) / sp.factorial(2) == sp.Rational(1, 12) and sp.bernoulli(4) / sp.factorial(4) == -sp.Rational(1, 720)
for n in [10, 20, 40]:
    hn = mp.mpf(1) / n
    em = (mp.e - 1) * (hn ** 2 / 12 - hn ** 4 / 720)
    assert abs(T(mp.exp, mp.mpf(0), mp.mpf(1), n) - I_e - em) < 3 * (mp.e - 1) * hn ** 6 / 30240
# Romberg : R_{k,j} = I + O(h^{2j+2}) sur exp ; la colonne 1 est Simpson ; la colonne 2 est Boole.
Rt = romberg(mp.exp, mp.mpf(0), mp.mpf(1), 6)
for j in range(0, 3):
    rapports = [(Rt[k][j] - I_e) / (Rt[k + 1][j] - I_e) for k in range(j + 2, 6)]
    assert abs(rapports[-1] - 4 ** (j + 1)) < 0.05 * 4 ** (j + 1)
for k in range(1, 5):
    assert abs(Rt[k][1] - S(mp.exp, mp.mpf(0), mp.mpf(1), 2 ** (k - 1))) < 1e-30
# exemple : ln 2.
Rl = romberg(lambda y: 1 / (1 + y), F(0), F(1), 2)
assert [Rl[0][0], Rl[1][0], Rl[2][0]] == [F(3, 4), F(17, 24), F(1171, 1680)]
assert Rl[1][1] == F(25, 36) and Rl[2][1] == F(1747, 2520) and Rl[2][2] == F(4367, 6300)
assert abs(mp.mpf(4367) / 6300 - mp.log(2) - mp.mpf("2.7e-5")) < 5e-7 and abs(mp.mpf(1171) / 1680 - mp.log(2) - mp.mpf("3.9e-3")) < 5e-5
# cours : pour 4/(1+t^2), R_{2,2} moins précis que R_{2,1}.
Rp = romberg(lambda y: 4 / (1 + y * y), F(0), F(1), 2)
assert abs(fl(Rp[2][2]) - mp.pi) > abs(fl(Rp[2][1]) - mp.pi)
# contre-exemples : rectangles f(x) = x ; cos sur [0, 2π] ; √t ; affirmations.
assert R(lambda y: y, F(0), F(1), 1) == 0 and R(lambda y: y, F(0), F(1), 2) == F(1, 4)
assert (4 * F(1, 4) - 0) / 3 == F(1, 3) and 2 * F(1, 4) - 0 == F(1, 2)
assert abs(T(mp.cos, mp.mpf(0), 2 * mp.pi, 1) - 2 * mp.pi) < 1e-30 and abs(T(mp.cos, mp.mpf(0), 2 * mp.pi, 2)) < 1e-30
assert abs((4 * T(mp.cos, mp.mpf(0), 2 * mp.pi, 2) - T(mp.cos, mp.mpf(0), 2 * mp.pi, 1)) / 3 + 2 * mp.pi / 3) < 1e-30
sqrtf = lambda y: mp.sqrt(y)  # noqa: E731
eT_ = [mp.mpf(2) / 3 - T(sqrtf, mp.mpf(0), mp.mpf(1), 2 ** k) for k in range(6, 13)]
eR_ = [mp.mpf(2) / 3 - (4 * T(sqrtf, mp.mpf(0), mp.mpf(1), 2 ** (k + 1)) - T(sqrtf, mp.mpf(0), mp.mpf(1), 2 ** k)) / 3 for k in range(6, 12)]
for liste in (eT_, eR_):
    for k in range(len(liste) - 1):
        assert abs(liste[k] / liste[k + 1] - 2 ** 1.5) < 0.05
for n in range(1, 7):  # affirmation : 2 R_{2n} - R_n = M_n ; R_{2n} = (R_n + M_n)/2
    f = lambda y: 1 / (2 + y * y)  # noqa: E731
    assert 2 * R(f, F(0), F(1), 2 * n) - R(f, F(0), F(1), n) == M(f, F(0), F(1), n)
# trous : c_1 et identité de Simpson (x^2 sur [0,1]) ; leurres.
for n in range(1, 6):
    hn = F(1, n)
    assert T(lambda y: y * y, F(0), F(1), n) - F(1, 3) == hn ** 2 * F(2 - 0, 12)
assert F(2, 12) != -F(2, 12) and F(2, 12) != F(2, 24) and F(2, 12) != 0
T1x, T2x = T(lambda y: y * y, F(0), F(1), 1), T(lambda y: y * y, F(0), F(1), 2)
assert (4 * T2x - T1x) / 3 == F(1, 3) == S(lambda y: y * y, F(0), F(1), 1)
assert M(lambda y: y * y, F(0), F(1), 1) != F(1, 3) and T(lambda y: y * y, F(0), F(1), 4) == F(11, 32) != F(1, 3)
assert R(lambda y: y * y, F(0), F(1), 2) != F(1, 3)
assert (T1x + M(lambda y: y * y, F(0), F(1), 1)) / 2 == F(3, 8)  # q8

# --- exercices de la carte Romberg ---
for (f, a_, b_, T1, T2, val, faux) in [
    (lambda y: 1 / (1 + y), F(0), F(1), F(3, 4), F(17, 24), F(25, 36), [F(2, 3), F(127, 180), F(35, 48)]),
    (lambda y: y ** 4, F(0), F(1), F(1, 2), F(9, 32), F(5, 24), [F(1, 16), F(4, 15), F(25, 64)]),
    (lambda y: 4 / (1 + y * y), F(0), F(1), F(3), F(31, 10), F(47, 15), [F(16, 5), F(233, 75), F(61, 20)]),
]:
    assert T(f, a_, b_, 1) == T1 and T(f, a_, b_, 2) == T2  # ro-ex1 (v1, v2)
    assert (4 * T2 - T1) / 3 == val == S(f, a_, b_, 1)
    assert faux == [2 * T2 - T1, (16 * T2 - T1) / 15, (T1 + T2) / 2] and val not in faux
assert abs(mp.mpf(17) / 24 - mp.log(2) - mp.mpf("1.5e-2")) < 5e-4 and abs(mp.mpf(25) / 36 - mp.log(2) - mp.mpf("1.3e-3")) < 5e-5
assert F(9, 32) - F(1, 5) == F(13, 160) and F(5, 24) - F(1, 5) == F(1, 120)
assert abs(mp.pi - mp.mpf(31) / 10 - mp.mpf("4.2e-2")) < 5e-4 and abs(mp.pi - mp.mpf(47) / 15 - mp.mpf("8.3e-3")) < 5e-5
# ro-ex2 : R_n = T_n - (h/2)(f(b) - f(a)) ; (4 R_{2n} - R_n)/3 pour f(x) = x, n = 1.
for n in range(1, 6):
    f = lambda y: y ** 3 - y  # noqa: E731
    assert R(f, F(-1), F(2), n) == T(f, F(-1), F(2), n) - F(3, n) / 2 * (f(F(2)) - f(F(-1)))
    assert 2 * R(f, F(-1), F(2), 2 * n) - R(f, F(-1), F(2), n) == M(f, F(-1), F(2), n)
    assert 2 * R(f, F(-1), F(2), 2 * n) - R(f, F(-1), F(2), n) not in (T(f, F(-1), F(2), n), S(f, F(-1), F(2), n), T(f, F(-1), F(2), 2 * n))
assert (R(lambda y: y, F(0), F(1), 1) + 2 * M(lambda y: y, F(0), F(1), 1)) / 3 == F(1, 3)
# ro-ex3 : tableau de Romberg.
for (f, a_, b_, Ts, R11, R21, R22, faux, I_) in [
    (lambda y: 1 / (1 + y), F(0), F(1), [F(3, 4), F(17, 24), F(1171, 1680)], F(25, 36), F(1747, 2520), F(4367, 6300),
     [F(97, 140), F(1747, 2520), F(2917, 4200)], mp.log(2)),
    (lambda y: 1 / y, F(1), F(3), [F(4, 3), F(7, 6), F(67, 60)], F(10, 9), F(11, 10), F(742, 675),
     [F(148, 135), F(11, 10), F(1501, 1350)], mp.log(3)),
    (lambda y: 1 / (y * y), F(1), F(3), [F(10, 9), F(29, 36), F(141, 200)], F(19, 27), F(1813, 2700), F(251, 375),
     [F(446, 675), F(1813, 2700), F(9529, 13500)], mp.mpf(2) / 3),
]:
    Rt_ = romberg(f, a_, b_, 2)
    assert [Rt_[0][0], Rt_[1][0], Rt_[2][0]] == Ts  # ro-ex3 (v1, v2)
    assert Rt_[1][1] == R11 and Rt_[2][1] == R21 and Rt_[2][2] == R22
    assert faux == [(4 * R21 - R11) / 3, R21, (16 * R11 - R21) / 15] and R22 not in faux
    eTs = [abs(mp.mpf(v.numerator) / v.denominator - I_) for v in (Ts[2], R21, R22)]
    assert eTs[0] > eTs[1] > eTs[2]
assert F(1171, 420) - F(17, 24) == F(2342 - 595, 840) and 16 * F(1747, 2520) == F(3494, 315)
assert F(3494, 315) - F(25, 36) == F(13976 - 875, 1260)
assert F(67, 15) - F(7, 6) == F(134 - 35, 30) and 16 * F(11, 10) == F(88, 5) and F(88, 5) - F(10, 9) == F(792 - 50, 45)
assert F(141, 50) - F(29, 36) == F(2538 - 725, 900) and 16 * F(1813, 2700) == F(7252, 675)
assert F(7252, 675) - F(19, 27) == F(7252 - 475, 675)
assert abs(mp.mpf(1747) / 2520 - mp.log(2) - mp.mpf("1.07e-4")) < 5e-7
assert abs(mp.mpf(67) / 60 - mp.log(3) - mp.mpf("1.8e-2")) < 5e-4 and abs(mp.mpf(11) / 10 - mp.log(3) - mp.mpf("1.4e-3")) < 5e-5
assert abs(mp.mpf(742) / 675 - mp.log(3) - mp.mpf("6.5e-4")) < 5e-6 and abs(mp.mpf(742) / 675 - mp.mpf("1.099259")) < 1e-6
assert abs(F(141, 200) - F(2, 3) - F(38, 1000)) < F(5, 10000) and abs(F(1813, 2700) - F(2, 3) - F(48, 10000)) < F(5, 10000)
assert abs(F(251, 375) - F(2, 3) - F(27, 10000)) < F(5, 10000)
# ro-ex4 : λ = 1/3, puis 1/15 ; illustration numérique sur exp.
assert (-F(1, 4)) / (F(1, 4) - 1) == F(1, 3) and (-F(1, 16)) / (F(1, 16) - 1) == F(1, 15)
T8 = T(mp.exp, mp.mpf(0), mp.mpf(1), 8)
assert abs((T8 - T4) / 3 + mp.mpf("2.2344e-3")) < 5e-8 and abs(I_e - T8 + mp.mpf("2.2368e-3")) < 5e-8
assert F(1, 3) / F(1, 15) == 5
for n in [50, 100]:
    hn = mp.mpf(1) / n
    assert abs((I_e - T(mp.exp, mp.mpf(0), mp.mpf(1), 2 * n)) / (T(mp.exp, mp.mpf(0), mp.mpf(1), 2 * n) - T(mp.exp, mp.mpf(0), mp.mpf(1), n)) - mp.mpf(1) / 3) < 1e-3
# Cas f'(b) = f'(a) : f(t) = t^4 - 2t^3 + t^2 ... on prend f = cos(2πt) + t^3(1-t)^3 sur [0,1] : f'(0) = f'(1) = 0,
# f'''(1) - f'''(0) ≠ 0 ; le rapport tend vers 1/15.
fq = lambda y: y ** 3 * (1 - y) ** 3  # noqa: E731
assert sp.diff(t ** 3 * (1 - t) ** 3, t).subs(t, 0) == 0 and sp.diff(t ** 3 * (1 - t) ** 3, t).subs(t, 1) == 0
assert sp.diff(t ** 3 * (1 - t) ** 3, t, 3).subs(t, 1) - sp.diff(t ** 3 * (1 - t) ** 3, t, 3).subs(t, 0) != 0
Iq = F(1, 140)
assert sp.integrate(t ** 3 * (1 - t) ** 3, (t, 0, 1)) == sp.Rational(1, 140)
n = 64
rq = (Iq - T(fq, F(0), F(1), 2 * n)) / (T(fq, F(0), F(1), 2 * n) - T(fq, F(0), F(1), n))
assert abs(rq - F(1, 15)) < F(1, 1000)
# ro-ex5 : formule de Boole.
w = [F(16, 15 * 12) * a - F(1, 15 * 6) * b for a, b in zip([1, 4, 2, 4, 1], [1, 0, 4, 0, 1])]
assert w == [F(7, 90), F(32, 90), F(12, 90), F(32, 90), F(7, 90)] and sum(w) == 1
deg_boole = degre(lambda g: sum(wj * g(F(j, 4)) for j, wj in enumerate(w)))
assert deg_boole == 5  # ro-ex5 (et pas 7)
for f in [lambda y: 1 / (1 + y), lambda y: y ** 7 - y ** 2]:
    assert romberg(f, F(0), F(1), 2)[2][2] == sum(wj * f(F(j, 4)) for j, wj in enumerate(w))
assert degre(lambda g: sum(F(c_, 12) * g(F(j, 4)) for j, c_ in enumerate([1, 4, 2, 4, 1]))) == 3
assert degre(lambda g: sum(F(c_, 8) * g(F(j, 4)) for j, c_ in enumerate([1, 2, 2, 2, 1]))) == 1
vv = [F(-1, 2), F(-1, 4), F(0), F(1, 4), F(1, 2)]
mom = lambda j: sum(wj * v ** j for wj, v in zip(w, vv))  # noqa: E731
assert mom(2) == F(1, 12) and mom(4) == F(1, 80) and mom(6) == F(1, 384) and 2 * F(1, 2) ** 7 / 7 == F(1, 448)
assert F(2, 90) * (F(7, 4) + 2) == F(1, 12) and F(2, 90) * (F(7, 16) + F(32, 256)) == F(1, 80)
assert F(2, 90) * (F(7, 64) + F(32, 4096)) == F(1, 384)
# Composite : R_{k,2} = Boole composée sur 2^{k-2} panneaux (k = 3, [0, 2]).
fb = lambda y: 1 / (1 + y * y)  # noqa: E731
Rb = romberg(fb, F(0), F(2), 3)
H = F(2, 2)
boole_comp = sum(H * sum(wj * fb(i * H + F(j, 4) * H) for j, wj in enumerate(w)) for i in range(2))
assert Rb[3][2] == boole_comp

# ===========================================================================
# QCM
# ===========================================================================
assert M(lambda y: y * y, F(0), F(1), 1) == F(1, 4) and degre(lambda g: g(F(1, 2))) == 1  # q1
for n in range(1, 6):  # q2
    assert F(1, 2) - R(lambda y: y, F(0), F(1), n) == F(1, 2 * n) > F(1, 4 * n)
assert M(lambda y: y * y, F(0), F(1), 1) <= F(1, 3) <= T(lambda y: y * y, F(0), F(1), 1)  # q3
assert not (F(1, 2) <= F(1, 3)) and not (F(1, 3) <= F(1, 4))
assert F(1, 5) - S(lambda y: y ** 4, F(0), F(1), 1) == F(-1, 120)  # q5, q6
# q7 : erreur d'une formule de degré 5 (Boole) pour t^6 : e ≠ 0, erreur totale (b-a) e h^6.
e6 = F(1, 7) - sum(wj * F(j, 4) ** 6 for j, wj in enumerate(w))
assert e6 != 0
for n in range(1, 5):
    hn = F(1, n)
    tot = F(1, 7) - sum(hn * sum(wj * (i * hn + F(j, 4) * hn) ** 6 for j, wj in enumerate(w)) for i in range(n))
    assert tot == e6 * hn ** 6
# q8 : T_4 pour x^2 = 11/32.
assert T(lambda y: y * y, F(0), F(1), 4) == F(11, 32)

print("num_l3_integration : toutes les vérifications passent.")
