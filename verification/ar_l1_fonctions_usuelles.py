"""Vérifications calculatoires du chapitre « ar-l1-fonctions-usuelles » (analyse réelle, L1) :
cartes, exercices et QCM de content/*/analyse-reelle/l1-ar-l1-fonctions-usuelles.json et content/qcm/.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/ar_l1_fonctions_usuelles.py (doit finir sans erreur).
"""

from fractions import Fraction as F
from math import gcd

import mpmath as mp
import sympy as sp
from sympy import E, Rational, asinh, cosh, diff, exp, expand, limit, log, oo, pi, simplify, sinh, sqrt, symbols, tanh

mp.mp.dps = 30
x, y, t = symbols("x y t", real=True)
xp = symbols("xp", positive=True)


def proche(a, b, eps=mp.mpf("1e-25")):
    return abs(mp.mpf(a) - mp.mpf(b)) < eps


# ===========================================================================
# ar-l1-fonctions-usuelles-puissances (carte)
# ===========================================================================
a_, b_ = symbols("a b", real=True)
# Dérivée de x^a = exp(a ln x) : a x^(a-1).
assert simplify(diff(exp(a_ * log(xp)), xp) - a_ * xp ** (a_ - 1)) == 0
# Prolongement en 0 : taux x^(a-1) -> 0 (a = 3/2), 1 (a = 1), +oo (a = 1/2) ; trous 2 (leurres a > 1, a >= 1/2, entier).
assert limit(xp ** Rational(1, 2), xp, 0, "+") == 0
assert limit(xp ** Rational(1, 2) / xp, xp, 0, "+") == oo        # a = 1/2 : pas dérivable (leurre a >= 1/2 faux)
assert limit(xp ** 1 / xp, xp, 0, "+") == 1                      # a = 1 : dérivable (leurre a > 1 faux)
assert limit(xp ** Rational(3, 2) / xp, xp, 0, "+") == 0         # a = 3/2 non entier, dérivable (leurre « entier » faux)
assert limit(xp ** pi / xp, xp, 0, "+") == 0                     # q2 : x^pi dérivable en 0
# Croissances comparées (énoncé) et leurres du trou 1, sur un exemple alpha = 1/2, beta = 3.
assert limit(log(xp) ** 3 / xp ** Rational(1, 2), xp, oo) == 0
assert limit(xp ** Rational(1, 2) / log(xp) ** 3, xp, oo) == oo
assert limit(xp ** Rational(1, 2) / exp(3 * xp), xp, oo) == 0
assert limit(exp(3 * xp) / xp ** Rational(1, 2), xp, oo) == oo
assert limit(xp ** Rational(1, 2) * (-log(xp)) ** 3, xp, 0, "+") == 0
# Pièges : alpha = 0 (ln x -> +oo), beta = -1 (x e^x -> +oo), alpha + beta > 0 avec alpha = -1, beta = 2 (x (ln x)^2 -> +oo).
assert limit(log(xp), xp, oo) == oo and limit(xp * exp(xp), xp, oo) == oo and limit(xp * log(xp) ** 2, xp, oo) == oo
# Trou 3 : (x^a)^b = x^(ab) ; leurres x^(a+b) (a=b=1, x=2 : 2 != 4), x^(a^b) (a=2, b=3, x=2 : 64 != 256), x^a x^b.
assert (2 ** 1) ** 1 != 2 ** 2 and (2 ** 2) ** 3 == 64 and 2 ** (2 ** 3) == 256
# Contre-exemple : ((-1)^2)^(1/2) = 1 alors que (-1)^1 = -1.
assert ((-1) ** 2) ** F(1, 2) == 1
# Exemples : x^(1/x) -> 1, x^x -> 1 en 0+, (ln x)^100 / sqrt(x) > 10^100 en x = 10^10 (environ 10^131).
assert limit(xp ** (1 / xp), xp, oo) == 1 and limit(xp ** xp, xp, 0, "+") == 1
q = mp.log(mp.mpf(10) ** 10) ** 100 / mp.sqrt(mp.mpf(10) ** 10)
assert q > mp.mpf(10) ** 100 and 131 < mp.log10(q) < 132
# Affirmation 2 (fausse) : (1/4)^sqrt2 ≈ 0,141 > (1/4)^(3/2) = 0,125.
v = mp.mpf(1) / 4
assert proche(v ** mp.sqrt(2), mp.mpf("0.1408"), mp.mpf("1e-4")) and v ** mp.sqrt(2) > v ** 1.5 == mp.mpf("0.125")
# Affirmation 4 (fausse) : (1/4)^(1/4) = 1/sqrt2 ≈ 0,707 > e^(-1/e) ≈ 0,692, avec 1/4 < 1/e.
assert v < 1 / mp.e and proche(v ** v, 1 / mp.sqrt(2)) and v ** v > mp.exp(-1 / mp.e)
assert proche(mp.exp(-1 / mp.e), mp.mpf("0.6922"), mp.mpf("1e-4"))

# ar-l1-fonctions-usuelles-puissances-ex1 (+ v1, v2) : exposants, puis leurres distincts de la bonne réponse.
r = F(2, 1) * F(3, 4) + F(3, 1) * F(1, 3) - F(4, 1) * F(-1, 2)          # 4^(3/4) 8^(1/3) / 16^(-1/2)
assert r == F(9, 2) and set([F(1, 2), F(15, 4), F(7, 2)]).isdisjoint({r})
assert F(3, 2) + 1 - 2 == F(1, 2)                                        # leurre : multiplier au lieu de diviser
r = 2 * F(1, 4) + 3 * F(2, 3) - 4 * F(-1, 4)                             # v1 : 9^(1/4) 27^(2/3) / 81^(-1/4)
assert r == F(7, 2) and r not in (F(3, 2), F(13, 4), F(5, 2))
r = 2 * F(3, 4) + 3 * F(-1, 3) - F(-1, 2)                                # v2 : 25^(3/4) 125^(-1/3) / 5^(-1/2)
assert r == 1 and r not in (0, F(1, 4), 3)
assert proche(mp.mpf(4) ** 0.75 * mp.mpf(8) ** (mp.mpf(1) / 3) / mp.mpf(16) ** -0.5, 16 * mp.sqrt(2))
assert proche(mp.mpf(25) ** 0.75 * mp.mpf(125) ** (-mp.mpf(1) / 3) / mp.mpf(5) ** -0.5, 5)

# ar-l1-fonctions-usuelles-puissances-ex2 (+ v1, v2) : quotients successifs -> 0 ; les leurres sont d'autres ordres.
for ordre in [
    [log(xp) ** 10, xp ** Rational(1, 10), xp ** 10, exp(sqrt(xp))],
    [log(xp) ** 3, xp ** Rational(1, 3), xp ** 3, exp(log(xp) ** 2)],
    [log(xp) ** 100, xp ** 100, exp(log(xp) ** 2), 2 ** xp],
]:
    for f1, f2 in zip(ordre, ordre[1:]):
        # On passe au logarithme du quotient, qui doit tendre vers -oo.
        assert limit(log(f1) - log(f2), xp, oo) == -oo, (f1, f2)
# x^(ln x) = e^((ln x)^2) (v2).
assert simplify(log(xp ** log(xp)) - log(xp) ** 2) == 0
# Remarque de la correction : (ln 10^10)^10 ≈ 4·10^13.
assert 4e13 < float(mp.log(mp.mpf(10) ** 10) ** 10) < 4.3e13

# ar-l1-fonctions-usuelles-puissances-ex3 (+ v1, v2).
for f, xm, ext, l0, taux in [
    (xp ** xp, 1 / E, exp(-1 / E), 1, -oo),
    (xp ** (2 * xp), 1 / E, exp(-2 / E), 1, -oo),
    (xp ** (1 / xp), E, exp(1 / E), 0, 0),
]:
    d = diff(f, xp)
    assert simplify(d.subs(xp, xm)) == 0 and simplify(f.subs(xp, xm) - ext) == 0
    assert limit(f, xp, 0, "+") == l0
    assert limit((f - l0) / xp, xp, 0, "+") == taux
assert limit(diff(xp ** xp, xp) / xp ** xp - (log(xp) + 1), xp, 2) == 0
# Signes des dérivées : x^x et x^(2x) décroissantes puis croissantes (min), x^(1/x) croissante puis décroissante (max).
assert diff(xp ** xp, xp).subs(xp, Rational(1, 10)) < 0 < diff(xp ** xp, xp).subs(xp, 1)
assert diff(xp ** (1 / xp), xp).subs(xp, 1) > 0 > diff(xp ** (1 / xp), xp).subs(xp, 4)
# Leurres : 1/e, e^(1/e), e^(-e) (ex3) ; e^(-1/e), e^(2/e), e^-2 (v1) ; e^(-1/e), e, 1/e (v2) différents des bonnes valeurs.
assert len({float(exp(-1 / E)), float(1 / E), float(exp(1 / E)), float(exp(-E))}) == 4
assert len({float(exp(-2 / E)), float(exp(-1 / E)), float(exp(2 / E)), float(exp(-2))}) == 4
assert len({float(exp(1 / E)), float(exp(-1 / E)), float(E), float(1 / E)}) == 4
assert abs(float(exp(-2 / E)) - 0.479) < 1e-3 and abs(float(exp(1 / E)) - 1.445) < 1e-3

# ar-l1-fonctions-usuelles-puissances-ex4 (+ v1, v2).
g = log(xp) / xp
assert simplify(diff(g, xp) - (1 - log(xp)) / xp ** 2) == 0 and g.subs(xp, E) == 1 / E
assert limit(g, xp, 0, "+") == -oo and limit(g, xp, oo) == 0
assert 4 ** 3 == 64 < 81 == 3 ** 4                                       # « tout a >= e » faux (a = 4, x = 3)
assert 2 ** 3 < 3 ** 2                                                   # « tout a > 1 » faux (a = 2, x = 3)
assert mp.sqrt(3) < mp.e and 3 < mp.e ** 2
s2, s3 = mp.sqrt(2), mp.sqrt(3)
assert s2 ** s3 < s3 ** s2 and proche(s2 ** s3, mp.mpf("1.8226"), mp.mpf("1e-3")) and proche(s3 ** s2, mp.mpf("2.1746"), mp.mpf("1e-3"))
# v1 : aucun a ; exemple a = e, x = 1 ; « a dans ]0,1] » faux (a = 1, x = 1/2 : (1/2)^1 < 1^(1/2)) ; 3^pi > pi^3.
assert 1 < mp.e and F(1, 2) < 1
assert mp.mpf(3) ** mp.pi > mp.pi ** 3 and proche(mp.mpf(3) ** mp.pi, mp.mpf("31.544"), mp.mpf("1e-2"))
assert proche(mp.pi ** 3, mp.mpf("31.006"), mp.mpf("1e-2"))
# v2 : a >= 1 ; g croissante sur ]0,1] de maximum g(1) = 0 ; a = 1/2, x = 1 : 1/2 < 1 ; « tout a >= e » exclut a = 1 qui convient.
assert g.subs(xp, 1) == 0 and all(1 ** F(1, k) == 1 for k in range(1, 5))
assert all(mp.mpf(1) ** xx >= xx ** 1 for xx in [mp.mpf(k) / 10 for k in range(1, 11)])   # a = 1 convient
assert (mp.mpf(1) / 2) ** (mp.mpf(1) / 3) > (mp.mpf(1) / 3) ** (mp.mpf(1) / 2)
assert proche((mp.mpf(1) / 3) ** 0.5, mp.mpf("0.5774"), mp.mpf("1e-4")) and proche((mp.mpf(1) / 2) ** (mp.mpf(1) / 3), mp.mpf("0.7937"), mp.mpf("1e-4"))
# Contrôle numérique global : e^x >= x^e sur une grille, et pour a = 3 l'inégalité a^x >= x^a échoue près de e.
assert all(mp.e ** xx >= xx ** mp.e for xx in [mp.mpf(k) / 7 for k in range(1, 100)])
assert mp.mpf(3) ** mp.e < mp.e ** 3

# ===========================================================================
# ar-l1-fonctions-usuelles-circulaires (carte)
# ===========================================================================
assert simplify(diff(sp.asin(x), x) - 1 / sqrt(1 - x ** 2)) == 0
assert simplify(diff(sp.acos(x), x) + 1 / sqrt(1 - x ** 2)) == 0
assert simplify(diff(sp.atan(x), x) - 1 / (1 + x ** 2)) == 0
assert simplify(diff(sp.atan(xp) + sp.atan(1 / xp), xp)) == 0
# Exemples.
assert sp.asin(Rational(1, 2)) == pi / 6 and sp.acos(Rational(-1, 2)) == 2 * pi / 3 and sp.atan(sqrt(3)) == pi / 3
assert sp.asin(sp.sin(5 * pi / 6)) == pi / 6
assert all(proche(mp.tan(mp.asin(u)), u / mp.sqrt(1 - u ** 2)) for u in [mp.mpf(k) / 10 for k in range(-9, 10)])
# Contre-exemples et affirmations : arcsin(sin pi) = 0 ; 2 arctan 2 ≈ 2,214 et arctan(-4/3) ≈ -0,927 ; arccos(cos 3pi/2) = pi/2.
assert sp.asin(sp.sin(pi)) == 0 and sp.acos(sp.cos(3 * pi / 2)) == pi / 2
assert proche(2 * mp.atan(2), mp.mpf("2.2143"), mp.mpf("1e-4")) and proche(mp.atan(-mp.mpf(4) / 3), mp.mpf("-0.9273"), mp.mpf("1e-4"))
assert all(proche(mp.acos(mp.cos(th)), abs(th)) for th in [mp.mpf(k) / 5 for k in range(-15, 16)])  # arccos(cos θ) = |θ| sur [-pi, pi]
# Piège 2 : theta = 3pi/4 dans [0, pi] mais arcsin(sin 3pi/4) = pi/4.
assert sp.asin(sp.sin(3 * pi / 4)) == pi / 4
# Taux de arcsin en 1 : (arcsin x - pi/2)/(x - 1) -> +oo.
assert limit((sp.asin(x) - pi / 2) / (x - 1), x, 1, "-") == oo
# Trous : leurres de la dérivée (arccos', argsh', sqrt(1-x^2)) et de arctan' (1/(1-x^2), 1+tan^2, -1/(1+x^2)) faux en x = 1/2.
u = Rational(1, 2)
assert len({1 / sqrt(1 - u ** 2), -1 / sqrt(1 - u ** 2), 1 / sqrt(1 + u ** 2), sqrt(1 - u ** 2)}) == 4
assert len({1 / (1 + u ** 2), 1 / (1 - u ** 2), 1 + sp.tan(u) ** 2, -1 / (1 + u ** 2)}) == 4

# ar-l1-fonctions-usuelles-circulaires-ex1 (+ v1, v2) et leurres (même sinus / cosinus / tangente, hors intervalle).
P = mp.pi
assert proche(mp.asin(mp.sin(9 * P / 7)), -2 * P / 7) and proche(mp.sin(-5 * P / 7), mp.sin(9 * P / 7))
assert proche(mp.acos(mp.cos(13 * P / 8)), 3 * P / 8) and proche(mp.cos(-3 * P / 8), mp.cos(13 * P / 8))
assert proche(mp.atan(mp.tan(11 * P / 9)), 2 * P / 9) and proche(mp.tan(-7 * P / 9), mp.tan(11 * P / 9))
assert not proche(mp.cos(5 * P / 8), mp.cos(13 * P / 8))

# ar-l1-fonctions-usuelles-circulaires-ex2 (+ v1, v2) et valeurs des leurres.
assert sp.simplify(sp.expand_trig(sp.cos(2 * sp.atan(Rational(1, 2))))) == Rational(3, 5)
assert sp.simplify(sp.expand_trig(sp.sin(2 * sp.atan(Rational(1, 2))))) == Rational(4, 5)
assert sp.simplify(sp.expand_trig(sp.sin(2 * sp.asin(Rational(1, 3)))) - 4 * sqrt(2) / 9) == 0
assert sp.simplify(sp.expand_trig(sp.cos(2 * sp.asin(Rational(1, 3))))) == Rational(7, 9)
assert sp.simplify(sp.expand_trig(sp.tan(2 * sp.atan(Rational(1, 3))))) == Rational(3, 4)
assert proche(mp.sin(2 * mp.atan(mp.mpf(1) / 3)), mp.mpf(3) / 5) and proche(mp.cos(2 * mp.atan(mp.mpf(1) / 3)), mp.mpf(4) / 5)

# ar-l1-fonctions-usuelles-circulaires-ex3 (+ v1, v2).
assert proche(mp.atan(mp.mpf(1) / 2) + mp.atan(mp.mpf(1) / 3), P / 4)
assert proche(mp.atan(2) + mp.atan(3), 3 * P / 4) and not proche(mp.atan(2) + mp.atan(3), 2 * P / 3)
assert F(1, 2) * F(1, 5) < 1 and (F(1, 2) + F(1, 5)) / (1 - F(1, 10)) == F(7, 9) and (F(7, 9) + F(1, 8)) / (1 - F(7, 72)) == 1
assert proche(mp.atan(mp.mpf(1) / 2) + mp.atan(mp.mpf(1) / 5) + mp.atan(mp.mpf(1) / 8), P / 4)
assert proche(mp.atan(2) + mp.atan(5) + mp.atan(8), 5 * P / 4)
assert (2 * F(1, 3)) / (1 - F(1, 9)) == F(3, 4) and (F(3, 4) + F(1, 7)) / (1 - F(3, 28)) == 1
assert proche(2 * mp.atan(mp.mpf(1) / 3) + mp.atan(mp.mpf(1) / 7), P / 4)
assert proche(2 * mp.atan(3) + mp.atan(7), 5 * P / 4)
# Leurre « formule appliquée sans précaution » : (2+3)/(1-6) = -1 ; (2+5)/(1-10) = -7/9 puis avec 8 : tangente 1.
assert F(2 + 3, 1 - 6) == -1 and (F(-7, 9) + 8) / (1 - F(-7, 9) * 8) == 1
# Formule de cos(alpha+beta) = cos alpha cos beta (1 - ab) : contrôle numérique.
for A, B in [(0.3, 2.5), (-4, 0.2), (1.5, 0.6)]:
    al, be = mp.atan(A), mp.atan(B)
    assert proche(mp.cos(al + be), mp.cos(al) * mp.cos(be) * (1 - mp.mpf(A) * B))

# ar-l1-fonctions-usuelles-circulaires-ex4 : f(x) = arcsin(2x/(1+x^2)) par morceaux, et leurres faux en sqrt(3).
fx = lambda z: mp.asin(2 * z / (1 + z * z))
for z in [mp.mpf(k) / 4 for k in range(-40, 41)]:
    if abs(z) <= 1:
        assert proche(fx(z), 2 * mp.atan(z))
    elif z > 1:
        assert proche(fx(z), P - 2 * mp.atan(z))
    else:
        assert proche(fx(z), -P - 2 * mp.atan(z))
z = mp.sqrt(3)
assert proche(fx(z), P / 3)
assert not any(proche(fx(z), w) for w in [2 * mp.atan(z), 2 * mp.atan(z) - P, -P - 2 * mp.atan(z)])
assert simplify(1 - (2 * x / (1 + x ** 2)) ** 2 - (1 - x ** 2) ** 2 / (1 + x ** 2) ** 2) == 0
assert simplify(diff(2 * x / (1 + x ** 2), x) - 2 * (1 - x ** 2) / (1 + x ** 2) ** 2) == 0

# ar-l1-fonctions-usuelles-circulaires-ex5 (+ v1, v2) : identités et sommes télescopiques.
k = symbols("k", positive=True, integer=True)
assert simplify(((2 * k + 1) - (2 * k - 1)) / (1 + (2 * k + 1) * (2 * k - 1)) - 1 / (2 * k ** 2)) == 0
assert simplify(1 / (1 + k * (k + 1)) - 1 / (k ** 2 + k + 1)) == 0
assert simplify(2 / (1 + (k - 1) * (k + 1)) - 2 / k ** 2) == 0
for n in [1, 2, 5, 30]:
    S = mp.fsum(mp.atan(mp.mpf(1) / (2 * j * j)) for j in range(1, n + 1))
    assert proche(S, mp.atan(2 * n + 1) - P / 4)
    S = mp.fsum(mp.atan(mp.mpf(1) / (j * j + j + 1)) for j in range(0, n + 1))
    assert proche(S, mp.atan(n + 1))
    S = mp.fsum(mp.atan(mp.mpf(2) / (j * j)) for j in range(1, n + 1))
    assert proche(S, mp.atan(n + 1) + mp.atan(n) - P / 4)
# Limites : pi/4, pi/2, 3pi/4 (vérifiées à n = 10^8 à 10^-7 près).
n = 10 ** 8
assert abs(mp.atan(2 * n + 1) - P / 4 - P / 4) < 1e-7 and abs(mp.atan(n + 1) - P / 2) < 1e-7
assert abs(mp.atan(n + 1) + mp.atan(n) - P / 4 - 3 * P / 4) < 1e-7

# QCM q3 : arctan x + arctan(1/x) = -pi/2 pour x < 0 ; q4 : arcsin(sin pi) = 0.
assert all(proche(mp.atan(z) + mp.atan(1 / z), -P / 2) for z in [-mp.mpf(k) / 3 for k in range(1, 20)])

# ===========================================================================
# ar-l1-fonctions-usuelles-hyperboliques (carte)
# ===========================================================================
assert simplify((cosh(x) ** 2 - sinh(x) ** 2).rewrite(exp)) == 1
assert simplify(diff(tanh(x), x) - (1 - tanh(x) ** 2)) == 0
assert simplify((cosh(a_ + b_) - cosh(a_) * cosh(b_) - sinh(a_) * sinh(b_)).rewrite(exp)) == 0
assert simplify((sinh(a_ + b_) - sinh(a_) * cosh(b_) - cosh(a_) * sinh(b_)).rewrite(exp)) == 0
assert simplify(asinh(x).rewrite(log) - log(x + sqrt(x ** 2 + 1))) == 0
assert simplify(diff(log(x + sqrt(x ** 2 + 1)), x) - 1 / sqrt(x ** 2 + 1)) == 0
X1 = symbols("X1", positive=True)
assert simplify(diff(log(X1 + 1 + sqrt((X1 + 1) ** 2 - 1)), X1) - 1 / sqrt((X1 + 1) ** 2 - 1)) == 0   # argch' sur ]1,+oo[
assert all(proche(mp.acosh(mp.cosh(s)), abs(s)) for s in [mp.mpf(k) / 3 for k in range(-12, 13)])        # argch(ch t) = |t| (q6)
assert all(proche(mp.atanh(s), mp.log((1 + s) / (1 - s)) / 2) for s in [mp.mpf(k) / 10 for k in range(-9, 10)])
assert simplify(diff(log((1 + x) / (1 - x)) / 2, x) - 1 / (1 - x ** 2)) == 0
# argch non dérivable en 1 : argch(x)/(x-1) -> +oo.
assert limit(sp.acosh(x) / (x - 1), x, 1, "+") == oo
# Exemple : ch(ln 2) = 5/4, sh(ln 2) = 3/4, th(ln 2) = 3/5, argsh(3/4) = argth(3/5) = ln 2.
assert simplify(cosh(log(2)).rewrite(exp)) == Rational(5, 4) and simplify(sinh(log(2)).rewrite(exp)) == Rational(3, 4)
assert simplify(tanh(log(2)).rewrite(exp)) == Rational(3, 5)
assert proche(mp.asinh(mp.mpf(3) / 4), mp.log(2)) and proche(mp.atanh(mp.mpf(3) / 5), mp.log(2))
# Trous : leurres -1, 0, ch(2x) de ch^2 - sh^2 ; leurres de argsh faux en x = 1/2 (ln(x+sqrt(x^2-1)) non défini) ;
# leurres de th' (1 + th^2, th^2 - 1, 1/sh^2) faux en x = 1.
assert proche(mp.cosh(1) ** 2 + mp.sinh(1) ** 2, mp.cosh(2)) and mp.cosh(2) > 3.7
h = mp.mpf(1) / 2
assert not proche(mp.log(mp.sqrt(h * h + 1) - h), mp.asinh(h)) and not proche(mp.log((1 + h) / (1 - h)) / 2, mp.asinh(h))
assert h * h - 1 < 0
thp = 1 - mp.tanh(1) ** 2
assert all(not proche(thp, w) for w in [1 + mp.tanh(1) ** 2, mp.tanh(1) ** 2 - 1, 1 / mp.sinh(1) ** 2])
# Affirmations : ch(2x) = 2 ch^2 x - 1 ; sh(-1) ≈ -1,175 < -1.
assert simplify((cosh(2 * x) - 2 * cosh(x) ** 2 + 1).rewrite(exp)) == 0
assert proche(mp.sinh(-1), mp.mpf("-1.1752"), mp.mpf("1e-4")) and mp.sinh(-1) < -1
assert all(mp.sinh(s) >= s for s in [mp.mpf(k) / 4 for k in range(0, 40)])

# ar-l1-fonctions-usuelles-hyperboliques-ex1 (+ v1, v2) : solutions et leurres.
Xs = symbols("Xs", positive=True)
assert sp.solve(sp.Eq((Xs - 1 / Xs) / 2, Rational(3, 4)), Xs) == [2]
assert sorted(sp.solve(sp.Eq((Xs + 1 / Xs) / 2, Rational(5, 4)), Xs)) == [Rational(1, 2), 2]
assert sp.solve(sp.Eq((Xs ** 2 - 1) / (Xs ** 2 + 1), Rational(3, 5)), Xs) == [2]          # e^x = 2 pour th x = 3/5
assert not proche(mp.sinh(mp.log(1.5)), 0.75) and not proche(mp.sinh(mp.log(7) / 2), 0.75) and not proche(mp.sinh(-mp.log(2)), 0.75)
assert not proche(mp.cosh(mp.log(1.25)), 1.25)
assert not proche(mp.tanh(mp.log(4)), 0.6) and not proche(mp.tanh(mp.log(0.6) / 2), 0.6)

# ar-l1-fonctions-usuelles-hyperboliques-ex2 (+ v1, v2).
for A, B, C, sols in [(3, -1, 3, [1, 2]), (5, 3, 5, [Rational(1, 4), 1]), (2, 1, 2, [Rational(1, 3), 1])]:
    eq = sp.Eq(A * (Xs + 1 / Xs) / 2 + B * (Xs - 1 / Xs) / 2, C)
    assert sorted(sp.solve(eq, Xs)) == sols
    f_ = lambda s: A * mp.cosh(s) + B * mp.sinh(s) - C
    assert all(proche(f_(mp.log(mp.mpf(sp.Rational(r).p) / sp.Rational(r).q)), 0) for r in sols)
assert not proche(3 * mp.cosh(mp.log(3)) - mp.sinh(mp.log(3)), 3) and not proche(3 * mp.cosh(-mp.log(2)) - mp.sinh(-mp.log(2)), 3)
assert not proche(5 * mp.cosh(2 * mp.log(2)) + 3 * mp.sinh(2 * mp.log(2)), 5) and not proche(5 * mp.cosh(-mp.log(2)) + 3 * mp.sinh(-mp.log(2)), 5)
assert not proche(2 * mp.cosh(mp.log(3)) + mp.sinh(mp.log(3)), 2) and not proche(2 * mp.cosh(-mp.log(2)) + mp.sinh(-mp.log(2)), 2)
assert proche(mp.cosh(-2 * mp.log(2)), mp.mpf(17) / 8) and proche(mp.sinh(-2 * mp.log(2)), -mp.mpf(15) / 8)
assert proche(mp.cosh(-mp.log(3)), mp.mpf(5) / 3) and proche(mp.sinh(-mp.log(3)), -mp.mpf(4) / 3)

# ar-l1-fonctions-usuelles-hyperboliques-ex3 (+ v1, v2).
for tt, sh2, ch2, th2 in [(F(1, 3), F(3, 4), F(5, 4), F(3, 5)), (F(1, 2), F(4, 3), F(5, 3), F(4, 5))]:
    yy = mp.atanh(mp.mpf(tt.numerator) / tt.denominator)
    assert proche(mp.sinh(2 * yy), mp.mpf(sh2.numerator) / sh2.denominator)
    assert proche(mp.cosh(2 * yy), mp.mpf(ch2.numerator) / ch2.denominator)
    assert proche(mp.tanh(2 * yy), mp.mpf(th2.numerator) / th2.denominator)
    assert sh2 == 2 * tt / (1 - tt ** 2)
yy = mp.atanh(mp.mpf(1) / 3)
assert proche(mp.cosh(yy), 3 / (2 * mp.sqrt(2))) and proche(mp.sinh(yy), 1 / (2 * mp.sqrt(2)))
yy = mp.asinh(mp.mpf(1) / 2)
assert proche(mp.cosh(yy), mp.sqrt(5) / 2) and proche(mp.cosh(2 * yy), 1.5)
assert not any(proche(mp.cosh(2 * yy), w) for w in [1.25, mp.sqrt(5) / 2, mp.sqrt(5)])

# ar-l1-fonctions-usuelles-hyperboliques-ex4 (+ v1, v2).
for z in [mp.mpf(k) / 4 for k in range(-20, 21) if k != 0]:
    assert proche(mp.atanh((z * z - 1) / (z * z + 1)), mp.log(abs(z)))
    assert proche(mp.asinh((z * z - 1) / (2 * z)), mp.sign(z) * mp.log(abs(z)))
    if z > 0:
        assert proche(mp.atanh((z - 1) / (z + 1)), mp.log(z) / 2)
    else:
        assert abs((z - 1) / (z + 1)) >= 1 if z != -1 else True       # hors du domaine de argth pour x <= 0
z = mp.mpf(2)
assert not proche(mp.atanh((z * z - 1) / (z * z + 1)), 2 * mp.log(z)) and not proche(mp.asinh((z * z - 1) / (2 * z)), mp.log((z * z + 1) / (2 * z)))

# ar-l1-fonctions-usuelles-hyperboliques-ex5 : θ(x) = 2 arctan(e^x) - pi/2.
th_ = lambda s: 2 * mp.atan(mp.e ** s) - P / 2
for s in [mp.mpf(k) / 3 for k in range(-15, 16)]:
    T_ = th_(s)
    assert -P / 2 < T_ < P / 2
    assert proche(mp.cos(T_), 1 / mp.cosh(s)) and proche(mp.sin(T_), mp.tanh(s)) and proche(mp.tan(T_), mp.sinh(s))
    assert proche(T_, mp.atan(mp.sinh(s))) and proche(mp.atanh(mp.sin(T_)), s) and proche(mp.asinh(mp.tan(T_)), s)
    assert proche(mp.diff(th_, s), 1 / mp.cosh(s))
# Leurres : argsh(sin θ) faux, argch(1/cos θ) = |x| faux pour x < 0, argth(tan θ) non défini quand sh x > 1.
s = mp.mpf(-1)
assert not proche(mp.asinh(mp.sin(th_(s))), s) and not proche(mp.acosh(1 / mp.cos(th_(s))), s) and mp.tan(th_(mp.mpf(2))) > 1

# QCM q5 : ch(ln 3) = 5/3, sh(ln 3) = 4/3.
assert simplify(cosh(log(3)).rewrite(exp)) == Rational(5, 3) and simplify(sinh(log(3)).rewrite(exp)) == Rational(4, 3)

# ===========================================================================
# ar-l1-fonctions-usuelles-equations-fonctionnelles (carte)
# ===========================================================================
# Contre-exemple : sin(2 pi x) nul sur les entiers, f(1/4) = 1 != 0.
assert sp.sin(2 * pi * Rational(1, 4)) == 1 and all(sp.sin(2 * pi * m) == 0 for m in range(-5, 6))
# Contre-exemple : f(x) = x (x >= 0), 2x (x < 0) additive sur [0, +oo[, f(-1) = -2.
fd = lambda s: s if s >= 0 else 2 * s
assert all(fd(p + q_) == fd(p) + fd(q_) for p in range(0, 6) for q_ in range(0, 6)) and fd(-1) == -2 != -fd(1)
# Trous : leurres du trou 1 réfutés par x^2, e^x, |x| ; trou 3 : e^(2x) réfute « f(x) = e^x », e^x réfute f(1) x et x^f(1).
assert expand((x + y) ** 2 - x ** 2 * y ** 2) != 0 and (x ** 2).subs(x, 2) != 2          # x^2 multiplicative, pas linéaire
assert all(abs(p + q_) <= abs(p) + abs(q_) for p in range(-3, 4) for q_ in range(-3, 4)) and abs(-1) != 1 * (-1)
assert simplify(exp(2 * (x + y)) - exp(2 * x) * exp(2 * y)) == 0 and exp(2 * x) != exp(x)
# Exemples : f(2) = 6 -> 3x ; f(1) = 2 -> 2^x ; f(4) = 2 -> sqrt(x).
assert 4 ** Rational(1, 2) == 2
# Affirmation 3 : 2 ln vérifie f(xy) = f(x) + f(y).
assert simplify(2 * log(xp * X1) - 2 * log(xp) - 2 * log(X1)) == 0

# ar-l1-fonctions-usuelles-equations-fonctionnelles-ex1 (+ v1, v2).
assert 4 * F(-5, 2) == -10 and F(12, 3) == 4 and -10 not in (10, -30, F(-24, 5))
assert proche(mp.mpf(3) ** -0.5, 1 / mp.sqrt(3)) and not any(proche(1 / mp.sqrt(3), w) for w in [-2.25, mp.mpf(1) / 3, mp.sqrt(3)])
c = 6 / mp.log(8)
assert proche(c, 2 / mp.log(2)) and proche(c * mp.log(mp.mpf(1) / 4), -4)

# ar-l1-fonctions-usuelles-equations-fonctionnelles-ex2 (+ v1, v2) : la solution vérifie la relation, pas les leurres.
for lam, f1, sol, leurres in [
    (2, 3, x ** 2 + 2 * x, [3 * x, 3 * x ** 2, 2 * x ** 2 + x]),
    (6, 1, 3 * x ** 2 - 2 * x, [x, x ** 2, 2 * x ** 2 - x]),
    (-4, 0, 2 * x - 2 * x ** 2, [0 * x, 2 * x ** 2 - 2 * x, x - x ** 2]),
]:
    rel = lambda f: expand(f.subs(x, x + y) - f - f.subs(x, y) - lam * x * y)
    assert rel(sol) == 0 and sol.subs(x, 1) == f1
    assert all(rel(l) != 0 for l in leurres)

# ar-l1-fonctions-usuelles-equations-fonctionnelles-ex3 : affines solutions de Jensen, x^2 ne l'est pas.
jen = lambda f: simplify(f.subs(x, (x + y) / 2) - (f + f.subs(x, y)) / 2)
assert jen(a_ * x + b_) == 0 and jen(x ** 2) != 0

# ar-l1-fonctions-usuelles-equations-fonctionnelles-ex4 : e^(kx) - 1 et -1 solutions ; e^(kx) + 1 et 1 non.
kk = symbols("kk", real=True)
rel4 = lambda f: simplify(f.subs(x, x + y) - f - f.subs(x, y) - f * f.subs(x, y))
assert rel4(exp(kk * x) - 1) == 0 and rel4(sp.Integer(-1) + 0 * x) == 0
assert rel4(exp(x) + 1) != 0 and rel4(sp.Integer(1) + 0 * x) != 0 and rel4(2 * x) != 0

# QCM q7 : c ln dérivable de dérivée c/x.
assert simplify(diff(kk * log(xp), xp) - kk / xp) == 0

# ===========================================================================
# ar-l1-fonctions-usuelles-decimal (carte)
# ===========================================================================


def chiffres(p, q_, n):
    """Les n premiers chiffres après la virgule du rationnel p/q (p >= 0), par la division posée."""
    r_, out = p % q_, []
    for _ in range(n):
        r_ *= 10
        out.append(r_ // q_)
        r_ %= q_
    return out


def periode_pure(p, q_):
    """Plus petite période T si p/q (0 < p < q) est périodique dès le premier chiffre, sinon None
    (exact : les chiffres sont périodiques dès le début si et seulement si la suite des restes revient au premier)."""
    N0, T = periode_tardive(p, q_)
    return T if N0 == 0 else None


def periode_tardive(p, q_):
    """(rang de début, période) minimaux du développement de p/q, par les restes successifs."""
    vus, r_, k_ = {}, p % q_, 0
    while r_ not in vus:
        vus[r_] = k_
        r_ = (r_ * 10) % q_
        k_ += 1
    return vus[r_], k_ - vus[r_]


# Exemple : 1/7 période 6 (<= 7), 1/6 périodique à partir du deuxième chiffre, 0,1212... = 4/33, sqrt2 entre 1,414 et 1,415.
assert chiffres(1, 7, 12) == [1, 4, 2, 8, 5, 7] * 2 and periode_pure(1, 7) == 6 and 10 ** 6 - 1 == 7 * 142857
assert chiffres(1, 6, 5) == [1, 6, 6, 6, 6] and periode_pure(1, 6) is None and periode_tardive(1, 6) == (1, 1)
assert F(12, 99) == F(4, 33) and chiffres(4, 33, 6) == [1, 2, 1, 2, 1, 2]
assert F(1414, 1000) ** 2 < 2 < F(1415, 1000) ** 2
# Période au plus q pour p/q (énoncé), contrôlé pour q <= 200.
for q_ in range(1, 201):
    for p in range(0, q_):
        N0, T = periode_tardive(p, q_)
        assert T <= q_
# 0,4999... = 0,5 : 0,4 + somme 9·10^-k (k = 2..n) = 0,5 - 10^-n.
assert all(F(4, 10) + sum(F(9, 10 ** j) for j in range(2, n + 1)) == F(1, 2) - F(1, 10 ** n) for n in range(2, 12))
assert all(sum(F(9, 10 ** j) for j in range(1, n + 1)) == 1 - F(1, 10 ** n) for n in range(1, 12))   # 0,999... = 1
# Algorithme des chiffres : a_k = floor(10^k x) - 10 floor(10^(k-1) x) est dans {0..9}, et d_n = somme ; sur des rationnels.
for xr in [F(22, 7), F(-5, 3), F(1, 2), F(355, 113), F(-1, 2)]:
    fl = lambda z: z.numerator // z.denominator
    d = [F(fl(10 ** m * xr), 10 ** m) for m in range(12)]
    a = [fl(xr)] + [fl(10 ** m * xr) - 10 * fl(10 ** (m - 1) * xr) for m in range(1, 12)]
    assert all(0 <= a[m] <= 9 for m in range(1, 12))
    assert all(d[m] == a[0] + sum(F(a[j], 10 ** j) for j in range(1, m + 1)) for m in range(12))
    assert all(d[m] <= xr < d[m] + F(1, 10 ** m) and d[m] <= d[m + 1] for m in range(11))
# Trou 1 : leurres réfutés par x = 0,5 (n = 0) et x = 0,55 (n = 0).
assert not (F(1, 2) < 0 + F(1, 10)) and not (F(1, 2) < F(0 + 1, 10)) and not (F(55, 100) < F(5, 10))
# Contre-exemple / affirmation : 0,1010010001... non périodique (vérifié sur ses 2000 premiers chiffres).
ch_ = []
m = 0
while len(ch_) < 2000:
    ch_ += [1] + [0] * m
    m += 1
ch_ = ch_[:2000]
assert not any(all(ch_[i + T] == ch_[i] for i in range(N, 2000 - T)) for T in range(1, 60) for N in range(0, 800))
# QCM q8 : 0,123123... = 41/333 ; 0,1999... = 0,2 ; 22/7 rationnel.
assert F(123, 999) == F(41, 333) and chiffres(41, 333, 9) == [1, 2, 3] * 3
assert sum(F(9, 10 ** j) for j in range(2, 30)) + F(1, 10) == F(2, 10) - F(1, 10 ** 29)

# ar-l1-fonctions-usuelles-decimal-ex1 (+ v1, v2) et leurres (autres développements).
assert F(36, 99) == F(4, 11) and chiffres(4, 11, 6) == [3, 6] * 3
assert chiffres(9, 25, 4) == [3, 6, 0, 0] and chiffres(2, 55, 6) != chiffres(4, 11, 6) and chiffres(3, 8, 3) == [3, 7, 5]
assert F(243, 990) == F(27, 110) and chiffres(27, 110, 7) == [2, 4, 5, 4, 5, 4, 5]
assert chiffres(245, 999, 6) == [2, 4, 5, 2, 4, 5] and chiffres(49, 198, 5) == [2, 4, 7, 4, 7] and chiffres(9, 37, 6) == [2, 4, 3, 2, 4, 3]
assert F(1285713, 999999) == F(9, 7) and 1285713 == 9 * 142857 and chiffres(9, 7, 6) == [2, 8, 5, 7, 1, 4] and 9 // 7 == 1
assert all(F(z) != F(9, 7) for z in (F(8, 7), F(10, 7), F(7, 5)))

# ar-l1-fonctions-usuelles-decimal-ex2 (+ v1, v2).
assert chiffres(1, 7, 2026)[2025] == 8 and 2026 == 6 * 337 + 4
assert chiffres(1, 13, 6) == [0, 7, 6, 9, 2, 3] and chiffres(1, 13, 100)[99] == 9 and 100 == 6 * 16 + 4
assert chiffres(3, 11, 2027)[2026] == 2

# ar-l1-fonctions-usuelles-decimal-ex3 : écriture de 10^m dans la constante de Champernowne (chiffres vérifiés jusqu'à 10^4).
champ = "".join(str(i) for i in range(1, 20001))
assert champ.startswith("123456789101112131415")
for mm in range(1, 5):
    pos = champ.find("1" + "0" * mm + str(10 ** mm + 1))
    assert pos >= 0 and champ[pos + 1:pos + 1 + mm] == "0" * mm

# ar-l1-fonctions-usuelles-decimal-ex4 (+ v1, v2) : périodique dès le premier chiffre <=> dénominateur premier avec 10.
for q_ in range(2, 300):
    for p in range(1, q_):
        if gcd(p, q_) == 1:
            assert (periode_pure(p, q_) is not None) == (gcd(q_, 10) == 1), (p, q_)
for bonne, autres in [((1, 21), [(1, 6), (1, 12), (1, 15)]),
                      ((5, 37), [(7, 30), (3, 44), (11, 50)]),
                      ((13, 99), [(4, 45), (2, 55), (9, 40)])]:
    assert periode_pure(*bonne) is not None and all(periode_pure(*z) is None for z in autres)
    assert all(gcd(*z) == 1 for z in [bonne] + autres)
assert chiffres(1, 21, 6) == [0, 4, 7, 6, 1, 9] and chiffres(5, 37, 3) == [1, 3, 5] and chiffres(3, 44, 5) == [0, 6, 8, 1, 8]

# ar-l1-fonctions-usuelles-decimal-ex5 (+ v1, v2) : plus petite période = ordre de 10 modulo q.
def ordre10(q_):
    return next(T for T in range(1, q_ + 1) if (10 ** T - 1) % q_ == 0)


for q_ in range(2, 400):
    if gcd(q_, 10) == 1:
        assert periode_pure(1, q_) == ordre10(q_)
assert ordre10(41) == 5 and 99999 == 41 * 2439 and 9999 == 9 * 11 * 101 and 999 == 27 * 37
assert ordre10(37) == 3 and 99 == 2 * 37 + 25
assert ordre10(101) == 4 and 9999 == 99 * 101 and 999 == 9 * 101 + 90
assert chiffres(1, 41, 10) == [0, 2, 4, 3, 9] * 2 and chiffres(1, 101, 8) == [0, 0, 9, 9] * 2
# Leurres : 40, 10 sont des périodes non minimales ; 8 n'en est pas une (41) ; 36, 6, 9 (37) ; 100, 8 (101) : 2 n'en est pas une.
assert (10 ** 40 - 1) % 41 == 0 and (10 ** 10 - 1) % 41 == 0 and (10 ** 8 - 1) % 41 != 0
assert all((10 ** T - 1) % 37 == 0 for T in (36, 6, 9)) and (10 ** 2 - 1) % 101 != 0 and (10 ** 100 - 1) % 101 == 0

# ===========================================================================
# Compléments de la relecture indépendante
# ===========================================================================
# QCM q1 : x^sqrt(x) = e^(sqrt(x) ln x) -> 1 en 0+ (leurres 0, +oo, e faux).
assert limit(sqrt(xp) * log(xp), xp, 0, "+") == 0 and limit(xp ** sqrt(xp), xp, 0, "+") == 1
# ar-l1-fonctions-usuelles-hyperboliques-ex4 : (x^2-1)/(x^2+1) vaut -1 en x = 0 (exclu) et reste dans ]-1,1[ ailleurs.
assert ((x ** 2 - 1) / (x ** 2 + 1)).subs(x, 0) == -1
assert simplify(1 + (x ** 2 - 1) / (x ** 2 + 1) - 2 * x ** 2 / (x ** 2 + 1)) == 0
assert simplify(1 - (x ** 2 - 1) / (x ** 2 + 1) - 2 / (x ** 2 + 1)) == 0
# ar-l1-fonctions-usuelles-hyperboliques-ex4 v2 : u + sqrt(u^2 + 1) = x pour x > 0, avec u = (x^2 - 1)/(2x).
u_ = (xp ** 2 - 1) / (2 * xp)
assert simplify(u_ ** 2 + 1 - ((xp ** 2 + 1) / (2 * xp)) ** 2) == 0 and simplify(u_ + (xp ** 2 + 1) / (2 * xp) - xp) == 0
# ar-l1-fonctions-usuelles-circulaires-ex5 v2 : télescopage à deux rangs, S_n = arctan(n+1) + arctan n - pi/4 (n <= 30).
for n_ in range(1, 31):
    assert proche(mp.fsum(mp.atan(mp.mpf(2) / k ** 2) for k in range(1, n_ + 1)), mp.atan(n_ + 1) + mp.atan(n_) - P / 4)

print("ar_l1_fonctions_usuelles : tout est vérifié.")
