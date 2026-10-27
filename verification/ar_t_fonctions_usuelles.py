"""Vérifications calculatoires du chapitre « ar-t-fonctions-usuelles » (analyse réelle, terminale) : cartes et
exercices de content/*/analyse-reelle/terminale-ar-t-fonctions-usuelles.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise : v1, v2). Les leurres des
questions à choix sont vérifiés faux.
Lancer : python3 verification/ar_t_fonctions_usuelles.py (doit finir sans erreur).
"""

import math

from sympy import (E, Abs, Rational, cos, diff, exp, expand, factor, limit, log, nsolve, oo, pi, simplify, sin,
                   solve, solveset, sqrt, symbols, S, Interval, Union, FiniteSet, EmptySet, tan)

x = symbols("x", real=True)
y = symbols("y", real=True)
h = symbols("h", real=True)
t = symbols("t", positive=True)
n = symbols("n", positive=True, integer=True)
P = "ar-t-fonctions-usuelles"


def egal(a, b):
    return simplify(a - b) == 0


def diff_(a, b):
    return simplify(a - b) != 0


# ===========================================================================
# Carte ar-t-fonctions-usuelles-exp
# ===========================================================================
K = f"{P}-exp"
# Démo : phi(x) = e^x e^-x a une dérivée nulle ; psi(x) = e^(x+b) e^-x aussi.
b = symbols("b", real=True)
assert diff(exp(x) * exp(-x), x) == 0  # K demo
assert simplify(diff(exp(x + b) * exp(-x), x)) == 0  # K demo
# Réciproque : h(x) = g(x) e^(-kx) avec g' = k g a une dérivée nulle (vérifié pour g = C e^(kx)).
k, C = symbols("k C", real=True)
assert simplify(diff(C * exp(k * x) * exp(-k * x), x)) == 0  # K demo
# Trou 1 : leurres faux. e^a + e^b en a = b = 0 ; e^(ab) en a = b = 1 ; e^a + e^b - 1 en a = b = 1.
assert exp(0) != 2 * exp(0)  # K trou 1, leurre e^a+e^b
assert exp(2) != exp(1)  # K trou 1, leurre e^(ab)
assert diff_(exp(2), 2 * E - 1)  # K trou 1, leurre e^a+e^b-1
# Trou 2 : la fonction nulle vérifie g(0) = 0, g(1) = 0 et la relation ; elle n'est pas e^(kx) (qui vaut 1 en 0).
assert exp(k * 0) == 1  # K trou 2 et contre-exemple 1
# Trou 3, leurre additif : g(x) = x (et contre-exemple g(x) = 2x) : additive, g(0) = 0 != 1.
assert egal((x + y), x + y) and egal(2 * (x + y), 2 * x + 2 * y)  # K trou 3 / contre-exemple 3
# Trou 3, leurre multiplicatif : g(x) = x^2 vérifie g(xy) = g(x)g(y), g(0) = 0 != 1.
assert expand((x * y) ** 2 - x**2 * y**2) == 0  # K trou 3, leurre g(xy)
# Trou 3 / piège / contre-exemple 2 : g(x) = exp(x + sin(2 pi x)) vérifie la relation sur les entiers, pas partout.
g = lambda u: exp(u + sin(2 * pi * u))
assert all(egal(g(S(m1) + m2), g(S(m1)) * g(S(m2))) for m1 in range(-4, 5) for m2 in range(-4, 5))  # K contre 2
assert egal(g(S(1)), E) and egal(g(Rational(1, 4)), exp(Rational(5, 4)))  # K contre 2
assert diff_(exp(Rational(5, 4)), exp(Rational(1, 4)))  # K contre 2 : g(1/4) != e^(1/4)
assert diff_(g(Rational(1, 4)) * g(Rational(1, 4)), g(Rational(1, 2)))  # K piège : la relation échoue hors des entiers
# Exemple : e^(2x) e^(1-x) / e^(x+1) = 1 ; e^(2x-1) = e^(x+3) <=> x = 4 ; g'(0) = 2 => g = e^(2x).
assert simplify(exp(2 * x) * exp(1 - x) / exp(x + 1)) == 1  # K exemple
assert solve(2 * x - 1 - (x + 3), x) == [4]  # K exemple
# Affirmations : (e^1)^2 = e^2 != e^(1^2) = e ; e^(0+0) = 1 != 2.
assert diff_(exp(1) ** 2, exp(1))  # K affirmation 2
assert exp(0) != 2  # K affirmation 1

# ---------------------------------------------------------------------------
# K-ex1 (+ variantes) : simplifications
# ---------------------------------------------------------------------------
A = exp(3 * x + 1) * exp(-x) / exp(x - 2) ** 2
assert egal(A, exp(5))  # K-ex1
for faux in [exp(3), exp(-3), exp(4 * x - 3)]:
    assert diff_(A, faux)  # K-ex1 leurres
# Erreurs annoncées : (e^(x-2))^2 = e^(2x-2) donne e^3 ; signe sur -4 donne e^-3 ; addition donne e^(4x-3).
assert egal((2 * x + 1) - (2 * x - 2), 3) and egal((2 * x + 1) - 2 * x - 4, -3) and egal((2 * x + 1) + (2 * x - 4), 4 * x - 3)
A1 = exp(2 * x - 3) * exp(x + 1) ** 3 / exp(5 * x)
assert egal(A1, 1)  # K-ex1 v1
for faux in [exp(-2), S(0), exp(10 * x)]:
    assert diff_(A1, faux)  # K-ex1 v1 leurres
assert egal(2 * x - 3 + 3 * x + 1 - 5 * x, -2) and egal(2 * x - 3 + 3 * x + 3 + 5 * x, 10 * x)  # K-ex1 v1 erreurs
A2 = exp(x + 2) ** 2 * exp(-3 * x) / exp(1 - x)
assert egal(A2, exp(3))  # K-ex1 v2
for faux in [exp(3 - x), exp(3 - 2 * x), E]:
    assert diff_(A2, faux)  # K-ex1 v2 leurres
assert egal(x + 4 - 3 * x - (1 - x), 3 - x) and egal(2 * x + 4 - 3 * x - 1 - x, 3 - 2 * x) and egal(2 * x + 2 - 3 * x - (1 - x), 1)

# ---------------------------------------------------------------------------
# K-ex2 (+ variantes) : équations
# ---------------------------------------------------------------------------
assert set(solve((x**2 + 1) - (3 * x - 1), x)) == {1, 2}  # K-ex2
assert egal(exp(x**2 + 1) - exp(x) ** 3 * exp(-1), exp(x**2 + 1) - exp(3 * x - 1))  # K-ex2 réécriture
for s_ in [1, 2]:
    assert egal(exp(S(s_) ** 2 + 1), exp(S(s_)) ** 3 * exp(-1))  # K-ex2 vérification
for faux in [{-2, -1}, {0, 3}, set()]:
    assert faux != {1, 2}  # K-ex2 leurres
assert set(solve(x**2 - 3 * x, x)) == {0, 3} and set(solve(x**2 + 3 * x + 2, x)) == {-2, -1}  # K-ex2 origine des leurres
assert set(solve(x**2 - 4 - (-3 * x), x)) == {-4, 1}  # K-ex2 v1
for s_ in [-4, 1]:
    assert egal(exp(S(s_) ** 2) * exp(-4), 1 / exp(S(s_)) ** 3)  # K-ex2 v1 vérification
assert set(solve(x**2 + 3 * x, x)) == {-3, 0} and set(solve(x**2 - 3 * x - 4, x)) == {-1, 4}  # K-ex2 v1 leurres
assert solve(2 * x + 3 + x**2, x) == [] and (2**2 - 4 * 3) == -8  # K-ex2 v2 : aucune solution réelle, Delta = -8
assert egal(x**2 + 2 * x + 3, (x + 1) ** 2 + 2)  # K-ex2 v2
assert set(solve(x**2 + 2 * x - 3, x)) == {-3, 1} and set(solve(x**2 - 2 * x - 3, x)) == {-1, 3} \
    and set(solve(x**2 + 2 * x, x)) == {-2, 0}  # K-ex2 v2 leurres (non vides donc faux)

# ---------------------------------------------------------------------------
# K-ex3 (+ variantes) : inéquations
# ---------------------------------------------------------------------------
X = symbols("X", positive=True)
assert expand((X - 1) * (X - E)) == expand(X**2 - (1 + E) * X + E)  # K-ex3


def sol_ineq(expr, rel):
    return solveset(rel(expr, 0), x, S.Reals)


from sympy import Le, Gt
s3 = sol_ineq(exp(2 * x) - (1 + E) * exp(x) + E, Le)
assert s3 == Interval(0, 1)  # K-ex3
for faux in [Interval(1, E), Union(Interval(-oo, 0), Interval(1, oo)), Interval(0, E)]:
    assert faux != s3  # K-ex3 leurres
assert expand((X - exp(-1)) * (X - exp(2))) == expand(X**2 - (exp(-1) + exp(2)) * X + E)  # K-ex3 v1
ee = math.e
assert all(((math.exp(2 * v) - (ee**-1 + ee**2) * math.exp(v) + ee) <= 0) == (-1 <= v <= 2)
           for v in [k_ / 100 for k_ in range(-400, 401) if k_ not in (-100, 200)])  # K-ex3 v1 : solution [-1, 2] (grille)
assert simplify(exp(-2) - (exp(-1) + exp(2)) * exp(-1) + E) == 0 and simplify(exp(4) - (exp(-1) + exp(2)) * exp(2) + E) == 0  # bornes
for faux in [(exp(-1), exp(2)), (1, 2)]:
    assert faux != (-1, 2)  # K-ex3 v1 leurres (bornes différentes)
assert ((exp(2 * Rational(-3, 2)) - (exp(-1) + exp(2)) * exp(Rational(-3, 2)) + E) > 0)  # K-ex3 v1 : x = -3/2 hors solution (leurre complémentaire faux)
assert ((exp(2 * Rational(1, 2)) - (exp(-1) + exp(2)) * exp(Rational(1, 2)) + E) <= 0)  # K-ex3 v1 : x = 1/2 solution, hors [1,2]
assert expand((X - 1) * (X - exp(3))) == expand(X**2 - (1 + exp(3)) * X + exp(3))  # K-ex3 v2
assert all(((math.exp(2 * v) - (1 + ee**3) * math.exp(v) + ee**3) > 0) == (v < 0 or v > 3)
           for v in [k_ / 100 for k_ in range(-400, 701) if k_ not in (0, 300)])  # K-ex3 v2 : ]-oo,0[ u ]3,+oo[ (grille)
assert simplify(exp(6) - (1 + exp(3)) * exp(3) + exp(3)) == 0  # K-ex3 v2 : 3 est racine (exclu)
assert not ((exp(0) ** 2 - (1 + exp(3)) * exp(0) + exp(3)) > 0)  # K-ex3 v2 : 0 exclu (leurre avec bornes fermées faux)
assert ((exp(2 * Rational(1, 2)) - (1 + exp(3)) * exp(Rational(1, 2)) + exp(3)) < 0)  # K-ex3 v2 : 1/2 non solution (leurres ]0,3[, ]-oo,1[)

# ---------------------------------------------------------------------------
# K-ex4 (+ variantes) : équation fonctionnelle
# ---------------------------------------------------------------------------
for (a_, val, kk, faux1, faux2) in [(1, exp(3), 3, exp(x + 2), exp(3) * x),
                                    (2, E, Rational(1, 2), exp(x - 1), E / 2 * x),
                                    (-1, exp(2), -2, exp(x + 3), -exp(2) * x)]:
    f_ = exp(kk * x)
    assert egal(f_.subs(x, a_), val)  # K-ex4 (v1, v2) : f(a) = valeur imposée
    assert egal(f_.subs(x, x + y), f_ * f_.subs(x, y))  # K-ex4 : relation fonctionnelle
    assert solve(exp(kk * a_) - val, kk) is not None
    for faux in [faux1, faux2]:
        assert egal(faux.subs(x, a_), val)  # K-ex4 leurres : vérifient f(a) = valeur...
        assert diff_(faux.subs(x, 2 * a_), faux.subs(x, a_) ** 2) or diff_(faux.subs(x, 0), faux.subs(x, 0) ** 2)  # ... mais pas la relation
    assert val != 0  # K-ex4 : la fonction nulle ne vérifie pas f(a) = valeur (leurre « ou f nulle »)
# Valeurs citées dans les pièges : e^(x+2) en 0 : e^2 != e^4 ; e^3 x : f(2) = 2e^3 != e^6, etc.
assert diff_(exp(2), exp(4)) and diff_(2 * exp(3), exp(6))  # K-ex4
assert diff_(exp(-1), exp(-2)) and diff_(2 * E, exp(2))  # K-ex4 v1
assert diff_(exp(3), exp(6)) and diff_(2 * exp(2), exp(4))  # K-ex4 v2
assert egal((E / 2 * x).subs(x, 4), 2 * E) and egal((-exp(2) * x).subs(x, -2), 2 * exp(2))

# ---------------------------------------------------------------------------
# K-ex5 : encadrement de e
# ---------------------------------------------------------------------------
for m in range(2, 60):
    lo, hi = (1 + Rational(1, m)) ** m, (1 + Rational(1, m - 1)) ** m
    assert lo <= E <= hi  # K-ex5 : bonne réponse
assert not ((1 + Rational(1, 1)) ** 2 <= E <= (1 + Rational(1, 2)) ** 2)  # K-ex5 leurre 2 (n = 2 : 4 <= e faux)
assert not (E <= (1 + Rational(1, 2)) ** 2)  # K-ex5 leurre 4 (n = 2 : e <= 9/4 faux)
assert not (E <= (1 + Rational(1, 2)) ** 2) and (1 - Rational(1, 2)) ** 2 <= E  # K-ex5 leurre 3 (borne supérieure fausse)
assert abs(1.1**10 - 2.594) < 1e-3 and abs((10 / 9) ** 10 - 2.868) < 1e-3  # K-ex5 valeurs pour n = 10
assert 2.59 <= 1.1**10 and (10 / 9) ** 10 <= 2.87  # K-ex5 arrondis dans le bon sens

# ===========================================================================
# Carte ar-t-fonctions-usuelles-ln
# ===========================================================================
K = f"{P}-ln"
assert egal(diff(log(x), x), 1 / x)  # K énoncé
assert egal(log(12) - log(3), 2 * log(2))  # K exemple
assert solve(2 * x - 1 - (x + 2), x) == [3] and 3 > Rational(1, 2)  # K exemple
assert egal(diff(log(x**2 + 1), x), 2 * x / (x**2 + 1))  # K exemple
# Trous : leurres faux.
assert diff_(log(E) * log(E), log(E * E))  # K trou 3 : ln a ln b (a = b = e) : 1 != 2
assert diff_(log(2), 0)  # K trou 3 : ln(a+b) (a = b = 1) ; contre-exemple 3 ; affirmation 3
assert diff_(2 * log(1), log(2))  # K trou 3 : a ln b (a = 2, b = 1) : 0 != ln 2
assert diff_(log(1), E) and diff_(log(E), 1 / E - 0) and diff_(log(10), 1)  # K trou 2 : e^x = y (x = 1 : y = e != ln 1) ; e^y = 1/x (x = e) ; 10^y = x (x = 10 : y = 1 != ln 10)
assert egal(exp(-log(E)), 1 / E)  # K trou 2 : avec y = ln x, e^y = x et non 1/x (sauf x = 1)
# Contre-exemples : a = b = -1 ; x = -e.
assert (-1) * (-1) == 1 and log(1) == 0  # K contre 1 / affirmation 1
assert egal(log(E**2), 2) and egal(log(Abs(-E)) * 2, 2)  # K contre 2 / affirmation 2
assert all(egal(log(S(v) ** 2), 2 * log(Abs(S(v)))) for v in [-5, -2, -1, 3, 7])  # K affirmation 2 (exemples)
assert all(egal(log(exp(S(m))), m) for m in range(0, 20))  # K affirmation 4

# K-ex1 (+ v1, v2)
A = log(16) - 3 * log(2) + log(Rational(1, 8))
assert egal(A, -2 * log(2))  # K-ex1
for faux in [4 * log(2), -10 * log(2), -log(2)]:
    assert diff_(A, faux)  # K-ex1 leurres
assert egal(log(16) - 3 * log(2) + log(8), 4 * log(2))  # K-ex1 erreur de signe annoncée
B = log(27) + 2 * log(Rational(1, 3)) - log(9)
assert egal(B, -log(3))  # K-ex1 v1
for faux in [3 * log(3), -7 * log(3), log(3)]:
    assert diff_(B, faux)  # K-ex1 v1 leurres
assert egal(log(27) + 2 * log(3) - log(9), 3 * log(3))  # K-ex1 v1 erreur annoncée
Cc = 2 * log(25) - log(125) + log(Rational(1, 5))
assert egal(Cc, 0)  # K-ex1 v2
for faux in [2 * log(5), -8 * log(5), log(5)]:
    assert diff_(Cc, faux)  # K-ex1 v2 leurres
assert egal(2 * log(25) - log(125) + log(5), 2 * log(5))  # K-ex1 v2 erreur annoncée

# K-ex2 (+ v1, v2) : ln(a) + ln(b) = ln(c)
for (a_, b_, c_, dom, bonne, autre, somme) in [(x + 3, x - 1, 5, 1, 2, -4, Rational(3, 2)),
                                                (x, x + 4, 21, 0, 3, -7, Rational(17, 2)),
                                                (x - 2, x + 1, 10, 2, 4, -3, Rational(11, 2))]:
    assert set(solve(a_ * b_ - c_, x)) == {bonne, autre}  # K-ex2 (v1, v2) : racines du trinôme
    assert bonne > dom and not autre > dom  # K-ex2 : seule la bonne est dans le domaine
    assert egal(log(a_.subs(x, bonne)) + log(b_.subs(x, bonne)), log(c_))  # K-ex2 : vérification
    assert solve(a_ + b_ - c_, x) == [somme] and not somme == bonne  # K-ex2 leurre ln(a+b)
    assert min(a_.subs(x, autre), b_.subs(x, autre)) <= 0  # K-ex2 leurres contenant autre : un ln non défini
assert expand((x - 2) * (x + 4)) == x**2 + 2 * x - 8 and expand((x - 3) * (x + 7)) == x**2 + 4 * x - 21 \
    and expand((x - 4) * (x + 3)) == x**2 - x - 12  # K-ex2 factorisations
assert expand((x + 3) * (x - 1) - 5) == x**2 + 2 * x - 8 and expand(x * (x + 4) - 21) == x**2 + 4 * x - 21 \
    and expand((x - 2) * (x + 1) - 10) == x**2 - x - 12  # K-ex2 trinômes

# K-ex3 (+ v1, v2) : minimum
for (f_, xm, vm, faux) in [(x - 2 * log(x), 2, 2 - 2 * log(2), [(1, 1), (2, 2 - log(2)), (Rational(1, 2), Rational(1, 2) + 2 * log(2))]),
                           (x - 3 * log(x), 3, 3 - 3 * log(3), [(1, 1), (3, 3 - log(3)), (Rational(1, 3), Rational(1, 3) + 3 * log(3))]),
                           (2 * x - log(x), Rational(1, 2), 1 + log(2), [(1, 2), (Rational(1, 2), 1 - log(2)), (2, 4 - log(2))])]:
    xp = symbols("xp", positive=True)
    fp = f_.subs(x, xp)
    assert solve(diff(fp, xp), xp) == [xm]  # K-ex3 : unique point critique
    assert diff(fp, xp, 2).subs(xp, xm) > 0  # K-ex3 : minimum
    assert egal(fp.subs(xp, xm), vm)  # K-ex3 : valeur du minimum
    for (xf, vf) in faux:
        assert xf != xm or diff_(vf, vm)  # K-ex3 leurres faux
        if xf != xm:
            assert egal(fp.subs(xp, xf), vf) or diff_(vf, vm)
    assert fp.subs(xp, 1) > fp.subs(xp, xm)  # K-ex3 : f(1) n'est pas le minimum
assert abs(float(2 - 2 * log(2)) - 0.614) < 1e-3 and abs(float(3 - 3 * log(3)) + 0.296) < 1e-3 \
    and abs(float(1 + log(2)) - 1.693) < 1e-3  # K-ex3 valeurs approchées
assert egal(1 - log(Rational(1, 2)), 1 + log(2))  # K-ex3 v2

# K-ex4 (+ v1, v2) : seuils
def seuil(q, s, inf):
    m = 0
    while not ((q**m < s) if inf else (q**m > s)):
        m += 1
    return m


from fractions import Fraction as Fr
assert seuil(Fr(9, 10), Fr(1, 100), True) == 44  # K-ex4
assert abs(math.log(0.01) / math.log(0.9) - 43.7) < 0.05  # K-ex4 quotient
assert abs(0.9**43 - 0.0108) < 1e-4 and abs(0.9**44 - 0.0097) < 1e-4  # K-ex4 vérification
assert seuil(Fr(9, 10), Fr(1, 10), True) == 22  # K-ex4 leurre 22 (confusion 0,01 / 0,1)
assert all(v != 44 for v in [43, 1, 22])  # K-ex4 leurres
assert seuil(Fr(105, 100), 2, False) == 15  # K-ex4 v1
assert abs(math.log(2) / math.log(1.05) - 14.2) < 0.05  # K-ex4 v1 quotient
assert abs(1.05**14 - 1.980) < 1e-3 and abs(1.05**15 - 2.079) < 1e-3  # K-ex4 v1 vérification
assert Fr(105, 100) ** 20 > 2 and 1 + Fr(5, 100) * 21 > 2 and not 1 + Fr(5, 100) * 20 > 2  # K-ex4 v1 leurres 20, 21 (non minimaux)
assert all(v != 15 for v in [14, 20, 21])  # K-ex4 v1 leurres
assert seuil(Fr(8, 10), Fr(1, 1000), True) == 31  # K-ex4 v2
assert abs(math.log(0.001) / math.log(0.8) - 30.96) < 0.005  # K-ex4 v2 quotient
assert abs(0.8**30 - 0.00124) < 1e-5 and abs(0.8**31 - 0.00099) < 1e-5  # K-ex4 v2 vérification
assert seuil(Fr(8, 10), Fr(1, 100), True) == 21  # K-ex4 v2 leurre 21
assert all(v != 31 for v in [30, 1, 21])  # K-ex4 v2 leurres

# K-ex5 : somme harmonique
for m in range(1, 300):
    H = sum(Fr(1, j) for j in range(1, m + 1))
    assert math.log(m + 1) <= float(H) <= 1 + math.log(m) + 1e-12  # K-ex5 bonne réponse
H1, H2 = 1, 1.5
assert not (math.log(1) <= H1 <= math.log(2))  # K-ex5 leurre 2 (n = 1)
assert not (1 + math.log(2) <= H2)  # K-ex5 leurre 3 (n = 2)
assert not (H1 <= 0.5 + math.log(1))  # K-ex5 leurre 4 (n = 1)
kk = symbols("kk", positive=True, integer=True)
assert egal(log(1 + 1 / kk), log(kk + 1) - log(kk))  # K-ex5 télescopage
assert egal(log(kk / (kk + 1)), -log(1 + 1 / kk))  # K-ex5 question 2
assert egal(kk / (kk + 1) - 1, -1 / (kk + 1))  # K-ex5 question 2

# ===========================================================================
# Carte ar-t-fonctions-usuelles-croissances
# ===========================================================================
K = f"{P}-croissances"
for m in range(1, 7):
    assert limit(exp(x) / x**m, x, oo) == oo  # K énoncé
    assert limit(x**m * exp(x), x, -oo) == 0  # K énoncé
    assert limit(log(x) / x**m, x, oo) == 0  # K énoncé
    assert limit(x**m * log(x), x, 0, "+") == 0  # K énoncé
assert limit(log(x), x, oo) == oo and limit(log(x), x, 0, "+") == -oo  # K contre 1 (n = 0) ; trou 1 leurre n >= 0
assert limit(log(x) / x, x, 0, "+") == -oo  # K contre 2 (n = -1) ; affirmation 3 ; trous leurres Z, n <= -1
assert limit(x * exp(x), x, oo) == oo  # K trou 1, n = -1 : e^x / x^(-1) -> +oo (mais la 2e limite échoue)
assert E**2 < 8 and abs(float(E**2) - 7.39) < 0.005  # K contre 3 ; affirmation 1
assert limit(exp(x) / x**100, x, oo) == oo  # K affirmation 2
assert limit(log(x) ** 2 / x, x, oo) == 0  # K affirmation 4
assert egal((2 * log(t) / t) ** 2, log(t**2) ** 2 / t**2)  # K affirmation 4 (X = sqrt x)
assert limit(x**3 - exp(x), x, oo) == -oo and limit(log(x) / sqrt(x), x, oo) == 0  # K exemple
assert limit(x * log(x), x, 0, "+") == 0  # K exemple
# Démo : phi = e^x - 1 - x a pour minimum 0 en 0 ; g = e^x - x^2/2, g' = e^x - x > 0.
assert solve(diff(exp(x) - 1 - x, x), x) == [0] and (exp(x) - 1 - x).subs(x, 0) == 0  # K démo
assert egal(diff(exp(x) - x**2 / 2, x), exp(x) - x)  # K démo
# Trous 2 et 3 : leurres faux (la limite vaut +oo, resp. 0).
assert limit(exp(x) / x**2, x, oo) not in (0, 1, 2) and limit(x**2 * log(x), x, 0, "+") not in (-oo, oo, 1)

# K-ex1 (+ v1, v2)
assert limit(exp(x) - x**3, x, oo) == oo  # K-ex1
assert limit(x**4 - 2 * exp(x), x, oo) == -oo  # K-ex1 v1
assert limit(3 * exp(x) - x**5 + x, x, oo) == oo  # K-ex1 v2
# Leurres : les autres propositions (-oo / +oo, 0, 3, « n'existe pas ») sont différentes de la limite, qui existe.

# K-ex2 (+ v1, v2)
assert limit((x**2 + exp(x)) / (exp(x) + x), x, oo) == 1  # K-ex2
assert limit((3 * exp(x) + x**3) / (exp(x) - x**2), x, oo) == 3  # K-ex2 v1
assert limit((x**4 - 2 * exp(x)) / (5 * exp(x) + x), x, oo) == Rational(-2, 5)  # K-ex2 v2
# v1 : e^x - x^2 > 0 sur [0, +oo[ (admis dans l'énoncé) : vrai, car (e^x - x^2)' = e^x - 2x > 0 et valeur 1 en 0.
assert solve(diff(exp(x) - 2 * x, x), x) == [log(2)] and (exp(x) - 2 * x).subs(x, log(2)) > 0  # K-ex2 v1
assert all(math.exp(v / 10) - (v / 10) ** 2 > 0 for v in range(0, 300))  # K-ex2 v1 (grille)
assert all(math.exp(v / 10) + v / 10 > 0 for v in range(0, 300))  # K-ex2 : dénominateur non nul sur [0, +oo[

# K-ex3 (+ v1, v2) : maximum
for (f_, xm, vm, faux) in [(x**2 * exp(-x), 2, 4 * exp(-2), [exp(-1), 2 * exp(-2)]),
                           (x**3 * exp(-x), 3, 27 * exp(-3), [exp(-1), 9 * exp(-3)]),
                           (x * exp(-2 * x), Rational(1, 2), 1 / (2 * E), [exp(-2), exp(-1)])]:
    crit = [c for c in solve(diff(f_, x), x) if c > 0]
    assert crit == [xm]  # K-ex3 (v1, v2) : point critique dans ]0, +oo[
    assert egal(f_.subs(x, xm), vm)  # K-ex3 : valeur maximale
    assert limit(f_, x, oo) == 0 and f_.subs(x, 0) == 0  # K-ex3 : limite et valeur en 0
    assert diff(f_, x).subs(x, xm / 2) > 0 and diff(f_, x).subs(x, 2 * xm) < 0  # K-ex3 : croissante puis décroissante
    for fv in faux:
        assert diff_(fv, vm)  # K-ex3 leurres
assert egal(diff(x**2 * exp(-x), x), x * (2 - x) * exp(-x))  # K-ex3
assert egal(diff(x**3 * exp(-x), x), x**2 * (3 - x) * exp(-x))  # K-ex3 v1
assert egal(diff(x * exp(-2 * x), x), (1 - 2 * x) * exp(-2 * x))  # K-ex3 v2
assert egal((x * exp(-x)).subs(x, 1) * 1, exp(-1)) and egal((x**2 * exp(-x)).subs(x, 1), exp(-1))  # K-ex3 : f(1) = e^-1
assert abs(float(4 * exp(-2)) - 0.541) < 1e-3 and abs(float(exp(-1)) - 0.368) < 1e-3  # K-ex3 valeurs approchées
assert abs(float(27 * exp(-3)) - 1.344) < 1e-3 and abs(float(1 / (2 * E)) - 0.184) < 1e-3  # K-ex3 v1, v2
assert egal((x * exp(-2 * x)).subs(x, 1), exp(-2))  # K-ex3 v2 : f(1) = e^-2

# K-ex4 (+ v1, v2) : nombre de solutions de e^x = x^m
for (m, nb, approx) in [(3, 2, [1.857, 4.536]), (2, 1, [-0.703]), (4, 3, [-0.815, 1.430, 8.613])]:
    racines = set()
    for x0 in [-3, -1, -0.5, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 9, 10, 12]:
        try:
            r = float(nsolve(exp(x) - x**m, x, x0))
            racines.add(round(r, 6))
        except Exception:
            pass
    assert len(racines) == nb  # K-ex4 (v1 : m = 2, v2 : m = 4) : nombre de solutions trouvées numériquement
    assert all(any(abs(r - a) < 2e-3 for r in racines) for a in approx)  # K-ex4 : valeurs approchées citées
    # Changements de signe de e^x - x^m sur une grille fine de [-20, 40] : exactement nb.
    pts = [i / 1000 for i in range(-20000, 40001)]
    signes = [math.exp(p) - p**m for p in pts]
    assert sum(1 for a_, b_ in zip(signes, signes[1:]) if a_ * b_ < 0) == nb
    gm = m - m * math.log(m)  # minimum de g(x) = x - m ln x, atteint en m
    assert (gm < 0) == (m > math.e)  # K-ex4 : signe du minimum selon ln m > 1
assert 3 * (1 - math.log(3)) < 0 and 2 * (1 - math.log(2)) > 0 and 4 * (1 - math.log(4)) < 0  # K-ex4, v1, v2
assert egal(diff(x - 3 * log(x), x), (x - 3) / x)  # K-ex4

# ===========================================================================
# Carte ar-t-fonctions-usuelles-limites-zero
# ===========================================================================
K = f"{P}-limites-zero"
assert limit((exp(x) - 1) / x, x, 0) == 1 and limit(log(1 + x) / x, x, 0) == 1  # K énoncé
assert limit(sin(x) / x, x, 0) == 1 and limit((cos(x) - 1) / x, x, 0) == 0  # K énoncé
assert egal((cos(h) - 1) / h, -(sin(h) / h) * (sin(h) / (1 + cos(h))))  # K démo cosinus
# Démo sinus : sin x < x < tan x et cos x < sin x / x < 1 sur ]0, pi/2[ (grille).
assert all(math.sin(v) < v < math.tan(v) and math.cos(v) < math.sin(v) / v < 1
           for v in [k_ / 1000 * math.pi / 2 for k_ in range(1, 1000)])  # K démo
assert limit(Abs(h) / h, h, 0, "+") == 1 and limit(Abs(h) / h, h, 0, "-") == -1  # K contre 1 ; affirmation 4
assert limit(1 / h, h, 0, "+") == oo  # K contre 2 (taux 1/h)
assert limit((1 - cos(x)) / x**2, x, 0) == Rational(1, 2)  # K contre 3 ; affirmation 2 ; cours
assert egal(1 - cos(x), sin(x) ** 2 / (1 + cos(x)))  # K contre 3
assert limit((exp(2 * x) - 1) / x, x, 0) == 2  # K affirmation 1
assert limit(sin(x) / x, x, oo) == 0  # K affirmation 3
assert limit(sin(4 * x) / x, x, 0) == 4 and limit((exp(x) - 1) / sin(x), x, 0) == 1  # K exemple
assert limit(log(x) / (x - 1), x, 1) == 1  # K exemple
# Trous : leurres faux. f'(0) au lieu de f'(a) : f = x^2, a = 1 ; f(a) au lieu de f'(a) : f = x^2 + 1, a = 0.
assert limit(((1 + h) ** 2 - 1) / h, h, 0) == 2 != 0  # K trou 2 (f'(a) = 2, f'(0) = 0, f(a) = 1)
assert limit((exp(h) - 1) / h, h, 0) not in (0, E, oo)  # K trou 3

# K-ex1 (+ v1, v2)
assert limit((exp(3 * x) - 1) / x, x, 0) == 3  # K-ex1
assert limit((exp(-2 * x) - 1) / x, x, 0) == -2  # K-ex1 v1
assert limit((exp(5 * x) - 1) / (2 * x), x, 0) == Rational(5, 2)  # K-ex1 v2
assert all(v != 3 for v in [1, Rational(1, 3), 0]) and all(v != -2 for v in [2, 1, Rational(-1, 2)]) \
    and all(v != Rational(5, 2) for v in [5, Rational(2, 5), Rational(1, 2)])  # K-ex1 leurres
# K-ex2 (+ v1, v2)
assert limit(sin(4 * x) / log(1 + x), x, 0) == 4  # K-ex2
assert limit(sin(3 * x) / (exp(2 * x) - 1), x, 0) == Rational(3, 2)  # K-ex2 v1
assert limit(log(1 + 5 * x) / sin(2 * x), x, 0) == Rational(5, 2)  # K-ex2 v2
assert all(v != 4 for v in [1, Rational(1, 4), 0]) and all(v != Rational(3, 2) for v in [Rational(2, 3), 1, 6]) \
    and all(v != Rational(5, 2) for v in [Rational(2, 5), 10, 1])  # K-ex2 leurres
# K-ex3 (+ v1, v2)
assert limit((1 - cos(x)) / x**2, x, 0) == Rational(1, 2)  # K-ex3
assert limit((1 - cos(3 * x)) / x**2, x, 0) == Rational(9, 2)  # K-ex3 v1
assert limit((1 - cos(2 * x)) / (x * sin(x)), x, 0) == 2  # K-ex3 v2
assert egal(1 - cos(2 * x), 2 * sin(x) ** 2)  # K-ex3 v2 (variante de calcul)
assert egal((1 - cos(2 * x)) / (x * sin(x)),
            4 * (sin(2 * x) / (2 * x)) ** 2 * (x / sin(x)) / (1 + cos(2 * x)))  # K-ex3 v2 (factorisation)
assert all(v != Rational(1, 2) for v in [0, 1, 2]) and all(v != Rational(9, 2) for v in [Rational(3, 2), 9, 0]) \
    and all(v != 2 for v in [1, 4, Rational(1, 2)])  # K-ex3 leurres
# K-ex4 (+ v1, v2) : (1 + a/n)^n
assert limit(x * log(1 + 2 / x), x, oo) == 2 and limit((1 + Rational(2) / n) ** n, n, oo) == exp(2)  # K-ex4
assert limit(x * log(1 - 3 / x), x, oo) == -3 and limit((1 - Rational(3) / n) ** n, n, oo) == exp(-3)  # K-ex4 v1
assert limit(x * log(1 + 1 / (2 * x)), x, oo) == Rational(1, 2) \
    and egal(limit((1 + 1 / (2 * n)) ** n, n, oo), sqrt(E))  # K-ex4 v2
assert 1 - Rational(3, 4) > 0 and not (1 - Rational(3, 3) > 0)  # K-ex4 v1 : n >= 4 nécessaire pour 1 - 3/n > 0
for bonne, faux in [(exp(2), [1, oo, E]), (exp(-3), [1, 0, exp(3)]), (sqrt(E), [1, exp(2), E])]:
    assert all(diff_(bonne, f_) if f_ != oo else True for f_ in faux)  # K-ex4 leurres
assert abs(float(exp(2)) - 7.39) < 0.005 and abs(float(exp(-3)) - 0.0498) < 5e-5 \
    and abs(float(sqrt(E)) - 1.649) < 5e-4  # K-ex4 valeurs approchées

# ===========================================================================
# Carte ar-t-fonctions-usuelles-trigo
# ===========================================================================
K = f"{P}-trigo"
assert diff(sin(x), x) == cos(x) and diff(cos(x), x) == -sin(x)  # K énoncé
assert egal(cos(x) ** 2 + sin(x) ** 2, 1)  # K énoncé
assert egal(cos(3 * (x + 2 * pi / 3)), cos(3 * x))  # K exemple ; affirmation 2
assert egal(diff(cos(3 * x), x), -3 * sin(3 * x))  # K exemple
assert diff(cos(x), x).subs(x, pi / 2) == -1 and sin(pi / 2) == 1  # K affirmation 1 ; trous 1-2
assert limit(sin(Abs(h)) / h, h, 0, "+") == 1 and limit(sin(Abs(h)) / h, h, 0, "-") == -1  # K contre 1
kq = symbols("kq", positive=True, integer=True)
assert egal(diff(sin(x**2), x).subs(x, sqrt(2 * kq * pi)), 2 * sqrt(2 * kq * pi))  # K contre 2 ; affirmation 3
assert (cos(pi / 2) + sin(pi / 2)) == 1 and (cos(-pi / 2) + sin(-pi / 2)) == -1 and cos(0) + sin(0) == 1  # K contre 3
assert cos(pi) != cos(0)  # K trou 3 : cos n'est pas pi-périodique
assert egal(diff(sin(-x), x), -cos(-x))  # K affirmation 4

# K-ex1 (+ v1, v2)
f_ = cos(2 * x - pi / 3)
assert egal(diff(f_, x), -2 * sin(2 * x - pi / 3)) and egal(diff(f_, x).subs(x, 0), sqrt(3))  # K-ex1
assert all(diff_(sqrt(3), v) for v in [-sqrt(3), sqrt(3) / 2, -1])  # K-ex1 leurres
assert egal(-2 * cos(-pi / 3), -1)  # K-ex1 erreur annoncée
f_ = sin(3 * x + pi / 4)
assert egal(diff(f_, x).subs(x, 0), 3 * sqrt(2) / 2)  # K-ex1 v1
assert all(diff_(3 * sqrt(2) / 2, v) for v in [sqrt(2) / 2, -3 * sqrt(2) / 2, 3])  # K-ex1 v1 leurres
f_ = cos(x / 2 + pi / 3)
assert egal(diff(f_, x).subs(x, 0), -sqrt(3) / 4)  # K-ex1 v2
assert all(diff_(-sqrt(3) / 4, v) for v in [-sqrt(3) / 2, sqrt(3) / 4, Rational(-1, 4)])  # K-ex1 v2 leurres
assert egal(-cos(pi / 3) / 2, Rational(-1, 4))  # K-ex1 v2 erreur annoncée

# K-ex2 (+ v1, v2) : parité et pi-périodicité (bonne : 0, 2, 1)
for (f_, impaire, pi_per) in [(cos(3 * x) * sin(x), True, True), (sin(2 * x) * cos(x), True, False),
                              (sin(3 * x) * sin(x), False, True)]:
    assert egal(f_.subs(x, -x), -f_) == impaire and egal(f_.subs(x, -x), f_) == (not impaire)  # K-ex2 parité
    assert egal(f_.subs(x, x + pi), f_) == pi_per  # K-ex2 période pi
    assert egal(f_.subs(x, x + pi), -f_) == (not pi_per)  # K-ex2 f(x + pi) = -f(x) sinon
    assert egal(f_.subs(x, x + 2 * pi), f_)  # K-ex2 2pi-périodique
assert egal((sin(2 * x) * cos(x)).subs(x, pi / 4), sqrt(2) / 2)  # K-ex2 v1 : f(pi/4) != 0

# K-ex3 (+ v1, v2) : inéquations sur [-pi, pi], vérifiées sur une grille fine (bornes exactes vérifiées à part).
def appartient(v, intervalles):
    return any((v > a if oa else v >= a) and (v < b if ob else v <= b) for (a, b, oa, ob) in intervalles)


grille = [-math.pi + k_ * 2 * math.pi / 7200 for k_ in range(7201)] + [-math.pi, math.pi, 0.0]
p = math.pi
cas = [  # (inégalité, bonne réponse, leurres)  intervalles : (a, b, a exclu, b exclu)
    (lambda c: 2 * c * c - c - 1 <= 1e-12,
     [(-2 * p / 3, 2 * p / 3, False, False)],
     [[(-p / 3, p / 3, False, False)], [(-p, -2 * p / 3, False, False), (2 * p / 3, p, False, False)],
      [(0, 2 * p / 3, False, False)]]),
    (lambda c: 2 * c * c + c - 1 < -1e-12,
     [(-p, -p / 3, True, True), (p / 3, p, True, True)],
     [[(-p, -p / 3, False, True), (p / 3, p, True, False)], [(-p / 3, p / 3, True, True)],
      [(-2 * p / 3, 2 * p / 3, True, True)]]),
    (lambda c: 4 * c * c - 3 >= -1e-12,
     [(-p, -5 * p / 6, False, False), (-p / 6, p / 6, False, False), (5 * p / 6, p, False, False)],
     [[(-p / 6, p / 6, False, False)], [(-5 * p / 6, -p / 6, False, False), (p / 6, 5 * p / 6, False, False)],
      [(-p / 3, p / 3, False, False)]]),
]
for (ineg, bonne, leurres) in cas:
    pts_test = [v for v in grille if min(abs(v - a) for iv in bonne for a in iv[:2]) > 1e-9]
    assert all(ineg(math.cos(v)) == appartient(v, bonne) for v in pts_test)  # K-ex3 (v1, v2) : bonne réponse
    for L in leurres:
        assert any(ineg(math.cos(v)) != appartient(v, L) for v in grille)  # K-ex3 leurres faux
# Bornes exactes : cos(2pi/3) = -1/2, cos(pi/3) = 1/2, cos(pi/6) = sqrt3/2, cos(5pi/6) = -sqrt3/2, cos(pi) = -1.
assert cos(2 * pi / 3) == Rational(-1, 2) and cos(pi / 3) == Rational(1, 2)  # K-ex3, v1
assert cos(pi / 6) == sqrt(3) / 2 and cos(5 * pi / 6) == -sqrt(3) / 2 and cos(pi) == -1  # K-ex3 v1, v2
c_ = symbols("c_")
assert expand((c_ - 1) * (2 * c_ + 1)) == 2 * c_**2 - c_ - 1  # K-ex3
assert expand((c_ + 1) * (2 * c_ - 1)) == 2 * c_**2 + c_ - 1  # K-ex3 v1
assert expand((2 * c_ - sqrt(3)) * (2 * c_ + sqrt(3))) == 4 * c_**2 - 3  # K-ex3 v2
assert (2 * c_**2 + c_ - 1).subs(c_, -1) == 0  # K-ex3 v1 : x = ±pi exclus (inégalité stricte)
assert (4 * cos(pi) ** 2 - 3) >= 0 and (2 * cos(2 * pi / 3) ** 2 - cos(2 * pi / 3) - 1) == 0  # bornes incluses

# K-ex4 : encadrements de sin
assert egal(diff(x - sin(x), x), 1 - cos(x)) and egal(diff(cos(x) - 1 + x**2 / 2, x), x - sin(x))  # K-ex4
assert egal(diff(sin(x) - x + x**3 / 6, x), cos(x) - 1 + x**2 / 2)  # K-ex4
assert all(v - v**3 / 6 <= math.sin(v) <= v for v in [k_ / 100 for k_ in range(0, 3000)])  # K-ex4 bonne réponse
assert not (1 <= math.sin(1))  # K-ex4 leurre 2 (x = 1)
assert math.sin(1) > 1 - 1 / 6 and abs(math.sin(1) - 0.841) < 1e-3 and abs(1 - 1 / 6 - 0.833) < 1e-3  # K-ex4 leurre 3 (x = 1)
assert not (1 - 0 <= math.sin(0))  # K-ex4 leurre 4 (x = 0)
assert limit(1 - x**2 / 6, x, 0) == 1  # K-ex4 question 4

# K-ex5 (+ v1, v2) : triangle d'aire maximale
th = symbols("th", real=True)
for (r, maxi, approx, droit, faux) in [(1, 3 * sqrt(3) / 4, 1.299, 1, [1, 3 * sqrt(3) / 2, sqrt(3) / 2]),
                                       (2, 3 * sqrt(3), 5.196, 4, [4, 6 * sqrt(3), 3 * sqrt(3) / 2]),
                                       (3, 27 * sqrt(3) / 4, 11.69, 9, [9, 27 * sqrt(3) / 2, 9 * sqrt(3) / 4])]:
    A_, B_, C_ = (0, r), (-r * sin(th), -r * cos(th)), (r * sin(th), -r * cos(th))
    assert egal(B_[0] ** 2 + B_[1] ** 2, r**2) and egal(C_[0] ** 2 + C_[1] ** 2, r**2)  # K-ex5 : points sur le cercle
    aire = Rational(1, 2) * ((B_[0] - A_[0]) * (C_[1] - A_[1]) - (C_[0] - A_[0]) * (B_[1] - A_[1]))
    f_ = r**2 * sin(th) * (1 + cos(th))
    assert egal(aire, f_)  # K-ex5 : aire (orientation directe, f > 0 sur ]0, pi[)
    assert egal(diff(f_, th), r**2 * (2 * cos(th) ** 2 + cos(th) - 1))  # K-ex5 question 2
    assert egal(r**2 * (2 * cos(th) - 1) * (cos(th) + 1), r**2 * (2 * cos(th) ** 2 + cos(th) - 1))  # K-ex5 factorisation
    assert egal(f_.subs(th, pi / 3), maxi) and abs(float(maxi) - approx) < 5e-3  # K-ex5 : maximum
    assert egal(f_.subs(th, pi / 2), droit)  # K-ex5 : aire du triangle rectangle
    assert all(float(f_.subs(th, k_ * math.pi / 1000)) <= float(maxi) + 1e-12 for k_ in range(1, 1000))  # K-ex5 : max global
    for fv in faux:
        assert diff_(fv, maxi)  # K-ex5 leurres
    # Triangle équilatéral pour theta = pi/3 : AB^2 = BC^2 = 3 r^2.
    AB2 = ((B_[0] - A_[0]) ** 2 + (B_[1] - A_[1]) ** 2).subs(th, pi / 3)
    BC2 = ((C_[0] - B_[0]) ** 2 + (C_[1] - B_[1]) ** 2).subs(th, pi / 3)
    assert egal(AB2, 3 * r**2) and egal(BC2, 3 * r**2)  # K-ex5 : équilatéral
    assert egal(((B_[0] - A_[0]) ** 2 + (B_[1] - A_[1]) ** 2), r**2 * (2 + 2 * cos(th)))  # K-ex5 : AB^2 = r^2 (2 + 2 cos)
    # Triangle rectangle en A pour theta = pi/2 : produit scalaire AB.AC nul.
    assert egal(((B_[0] - A_[0]) * (C_[0] - A_[0]) + (B_[1] - A_[1]) * (C_[1] - A_[1])).subs(th, pi / 2), 0)

print("ar_t_fonctions_usuelles : toutes les vérifications passent.")
