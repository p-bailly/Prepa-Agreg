"""Vérifications calculatoires du chapitre « proba-l3-densites » (probabilités, L3, variables à densité) :
cartes, exercices et QCM de content/*/probabilites/l3-proba-l3-densites.json et content/qcm/proba-l3-densites.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Pour les réponses à choix, `choix(...)` relit le fichier JSON : il vérifie que les propositions affichées sont exactement
celles qu'on lui donne, que la bonne vaut la valeur calculée ici indépendamment, et que chaque leurre en diffère.
Lancer : python3 verification/proba_l3_densites.py (doit finir sans erreur).
"""

import json
import math
import pathlib
from fractions import Fraction as Fr

import mpmath as mp
from sympy import (E, Rational as R, Symbol, atan, cos, diff, erf, exp, factorial, gamma, integrate, limit, log, oo,
                   pi, simplify, sin, sqrt, symbols, tan)

mp.mp.dps = 30
RACINE = pathlib.Path(__file__).resolve().parent.parent / "content"
CH = "proba-l3-densites"
EXOS = {e["id"]: e for e in json.load(open(RACINE / "exercices/probabilites/l3-proba-l3-densites.json"))["exercices"]}
CARTES = {c["id"]: c for c in json.load(open(RACINE / "cartes/probabilites/l3-proba-l3-densites.json"))["cartes"]}
QCM = {q["id"]: q for q in json.load(open(RACINE / f"qcm/{CH}.json"))["questions"]}

x, y, z, t, u, s = symbols("x y z t u s", real=True)
lam, a, b, sig, m = symbols("lambda a b sigma m", positive=True)


def num(v):
    return complex(v.evalf(40)) if hasattr(v, "evalf") else complex(float(v))


def egal(v, w):
    return abs(num(v) - num(w)) < 1e-12


def choix(ident, v, vrai, table):
    """ident : identifiant court de l'exercice ; v : 0 (exercice) ou 1, 2 (variantes) ; table : {latex: valeur}."""
    e = EXOS[f"{CH}-{ident}"]
    rep = e["reponse_finale"] if v == 0 else e["variantes"][v - 1]["reponse_finale"]
    props = rep["propositions"]
    assert set(props) == set(table) and len(props) == 4, (ident, v, props)
    bonne = props[rep["bonne"]]
    assert egal(table[bonne], vrai), (ident, v, bonne)
    for p in props:
        if p != bonne:
            assert not egal(table[p], vrai), (ident, v, p)


def choix_qcm(ident, vrai, table):
    qq = QCM[f"{CH}-{ident}"]
    assert set(qq["propositions"]) == set(table)
    bonne = qq["propositions"][qq["bonne"]]
    assert egal(table[bonne], vrai), ident
    for p in qq["propositions"]:
        if p != bonne:
            assert not egal(table[p], vrai), (ident, p)


def Phi(v):
    return (1 + erf(v / sqrt(2))) / 2


def I(f, *bornes):
    return simplify(integrate(f, *bornes))


# ===========================================================================================
# Carte proba-l3-densites-densite
# ===========================================================================================
# Exemple : F = x^2 sur [0,1] -> densité 2x ; P(1/4 < X <= 1/2) = 3/16 ; f(1) = 2 ; dérivées à gauche/droite en 1 : 2 et 0.
assert I(2 * x, (x, 0, 1)) == 1 and R(1, 4) - R(1, 16) == R(3, 16)
assert diff(x**2, x).subs(x, 1) == 2
# Contre-exemple « mixte » : P(max(U,1/2) = 1/2) = P(U <= 1/2) = 1/2.
# Affirmation 1 : 2*1_[0,1/2] d'intégrale 1. Affirmation 3 : 1/(2 sqrt x) d'intégrale 1 sur ]0,1], non bornée.
assert 2 * R(1, 2) == 1 and I(1 / (2 * sqrt(x)), (x, 0, 1)) == 1 and limit(1 / (2 * sqrt(x)), x, 0, "+") == oo
# Affirmation 5 : F_{aX}(x) = F_X(x/a) a pour dérivée f(x/a)/a (testé avec f = densité de E(1), a = 3).
F1 = 1 - exp(-x)
assert simplify(diff(F1.subs(x, x / 3), x) - exp(-x / 3) / 3) == 0
# Pièges et trou 2 (loi uniforme, x0 = 0) : f continue à droite en 0, f(0) = 1 > 0, F continue en 0,
# mais dérivée à gauche 0 et à droite 1 : F non dérivable en 0.
assert limit((0 - 0) / t, t, 0, "-") == 0 and limit((t - 0) / t, t, 0, "+") == 1
# Trou 1 (U([0,1]), a = 0, b = 1/2) : int f = 1/2 ; leurres f(b)-f(a) = 0, int F = 1/8, int t f = 1/8.
assert I(1, (t, 0, R(1, 2))) == R(1, 2)
assert 1 - 1 == 0 != R(1, 2) and I(t, (t, 0, R(1, 2))) == R(1, 8) != R(1, 2)

# --- densite-ex1 et variantes : constante et probabilité.
c1 = 1 / I(x * (2 - x), (x, 0, 2))
assert c1 == R(3, 4)
v = c1 * I(x * (2 - x), (x, 0, R(1, 2)))
choix("densite-ex1", 0, v, {r"$\frac5{32}$": R(5, 32), r"$\frac5{24}$": I(x * (2 - x), (x, 0, R(1, 2))),
                            r"$\frac{27}{32}$": 1 - v, r"$\frac14$": R(1, 4)})
c2 = 1 / I(x**2, (x, 0, 3))
assert c2 == R(1, 9)
v = c2 * I(x**2, (x, 0, 2))
choix("densite-ex1", 1, v, {r"$\frac8{27}$": R(8, 27), r"$\frac83$": I(x**2, (x, 0, 2)), r"$\frac{19}{27}$": 1 - v,
                            r"$\frac23$": R(2, 3)})
c3 = 1 / I(4 - x**2, (x, -2, 2))
assert c3 == R(3, 32)
v = c3 * I(4 - x**2, (x, -2, 1))
choix("densite-ex1", 2, v, {r"$\frac{27}{32}$": R(27, 32), r"$\frac5{32}$": 1 - v, r"$\frac34$": R(3, 4),
                            r"$\frac{11}{32}$": c3 * I(4 - x**2, (x, 0, 1))})

# --- densite-ex2 et variantes : densité par dérivation (continuité en 0), médiane.
for F, dens, med, ver in [(1 - 1 / (1 + x)**2, 2 / (1 + x)**3, sqrt(2) - 1, 0),
                          (1 - 1 / (1 + x)**3, 3 / (1 + x)**4, 2**R(1, 3) - 1, 1),
                          (x**2 / (1 + x**2), 2 * x / (1 + x**2)**2, 1, 2)]:
    assert F.subs(x, 0) == 0 and limit(F, x, oo) == 1 and simplify(diff(F, x) - dens) == 0
    assert simplify(F.subs(x, med) - R(1, 2)) == 0
choix("densite-ex2", 0, sqrt(2) - 1, {r"$\sqrt2-1$": sqrt(2) - 1, r"$1$": 1, r"$\frac12$": R(1, 2), r"$\sqrt2$": sqrt(2)})
choix("densite-ex2", 1, 2**R(1, 3) - 1, {r"$\sqrt[3]2-1$": 2**R(1, 3) - 1, r"$\sqrt2-1$": sqrt(2) - 1,
                                         r"$\frac12$": R(1, 2), r"$\sqrt[3]2$": 2**R(1, 3)})
choix("densite-ex2", 2, 1, {r"$1$": 1, r"$\sqrt2-1$": sqrt(2) - 1, r"$\frac12$": R(1, 2), r"$\sqrt2$": sqrt(2)})


# --- densite-ex3 et variantes : Y = X^2, X uniforme sur [p, q] ; P(Y <= y) = longueur([p,q] ∩ [-√y, √y]) / (q - p).
def FY(p, q, yy):
    lo, hi = max(p, -sqrt(yy)), min(q, sqrt(yy))
    return (hi - lo) / (q - p)


for p, q, rac in [(-1, 2, 1), (-1, 3, 1), (-2, 3, 4)]:
    # formules par morceaux de la correction et continuité au raccord y = rac, densités qui s'intègrent à 1
    hi = q**2
    assert simplify(FY(p, q, R(rac, 2)) - 2 * sqrt(R(rac, 2)) / (q - p)) == 0
    assert simplify(FY(p, q, rac + R(1, 2)) - (sqrt(rac + R(1, 2)) - p) / (q - p)) == 0
    assert 2 * sqrt(rac) / (q - p) == (sqrt(rac) - p) / (q - p)
    assert I(diff(2 * sqrt(y) / (q - p), y), (y, 0, rac)) + I(diff((sqrt(y) - p) / (q - p), y), (y, rac, hi)) == 1
choix("densite-ex3", 0, FY(-1, 2, 2), {r"$\frac{1+\sqrt2}3$": (1 + sqrt(2)) / 3, r"$\frac{2\sqrt2}3$": 2 * sqrt(2) / 3,
                                       r"$\frac{\sqrt2}3$": sqrt(2) / 3, r"$\frac23$": R(2, 3)})
choix("densite-ex3", 1, FY(-1, 3, 4), {r"$\frac34$": R(3, 4), r"$1$": 1, r"$\frac12$": R(1, 2), r"$\frac14$": R(1, 4)})
choix("densite-ex3", 2, FY(-2, 3, 6), {r"$\frac{2+\sqrt6}5$": (2 + sqrt(6)) / 5, r"$\frac{2\sqrt6}5$": 2 * sqrt(6) / 5,
                                       r"$\frac{\sqrt6}5$": sqrt(6) / 5, r"$\frac45$": R(4, 5)})
assert FY(-1, 2, 1) == R(2, 3) and FY(-2, 3, 4) == R(4, 5)  # valeurs F_Y(1), F_Y(4) citées

# --- densite-ex4 : g(u) = tan(pi u - pi/2) vérifie F(g(u)) = u ; les leurres ne donnent pas la loi de Cauchy.
Fc = R(1, 2) + atan(x) / pi
for uu in [R(1, 10), R(1, 3), R(1, 2), R(4, 5)]:
    assert egal(Fc.subs(x, tan(pi * uu - pi / 2)), uu)
# tan(pi u / 2) >= 0 sur ]0,1[ (P(g(U) < 0) = 0 au lieu de 1/2) ; 1/(pi(1+u^2)) dans [1/(2pi), 1/pi] ;
# arctan(u)/pi + 1/2 dans [1/2, 3/4] : ces variables sont bornées ou positives, pas de Cauchy.
assert all(math.tan(math.pi * (k / 200) / 2) >= 0 for k in range(1, 200))

# ===========================================================================================
# Carte proba-l3-densites-lois-usuelles
# ===========================================================================================
assert I(lam * exp(-lam * x), (x, 0, oo)) == 1 and I(lam * exp(-lam * x), (x, t, oo)).subs(t, 2) == exp(-2 * lam)
assert I(1 / (pi * (1 + x**2)), (x, -oo, oo)) == 1
assert simplify(I(1 / (sig * sqrt(2 * pi)) * exp(-(x - m)**2 / (2 * sig**2)), (x, -oo, oo))) == 1
for aa in [R(1, 2), 1, R(5, 2), 3]:
    assert simplify(I(lam**aa / gamma(aa) * x**(aa - 1) * exp(-lam * x), (x, 0, oo))) == 1
# Exemple : e^{-1/2} ≈ 0,61 ; Cauchy P(X > 1) = 1/4 ; gamma(2,1) : P(X > t) = (1+t) e^{-t}.
assert round(math.exp(-0.5), 2) == 0.61 and R(1, 2) - atan(1) / pi == R(1, 4)
assert simplify(I(x * exp(-x), (x, t, oo)) - (1 + t) * exp(-t)) == 0
# Absence de mémoire de E(lambda).
assert simplify(exp(-lam * (s + t)) / exp(-lam * s) - exp(-lam * t)) == 0
# Pièges : loi géométrique (p = 1/3) sans mémoire aux entiers, et pour P(X > s+1 | X > s) avec s réel.
pg = Fr(1, 3)
G = lambda r: (1 - pg) ** math.floor(r)
assert all(G(k + n) / G(k) == G(n) for k in range(6) for n in range(6))
assert all(G(sv + 1) / G(sv) == G(1) for sv in [0, 0.3, 1.5, 2.99, 7.2])
# Piège « inégalité » : gamma(2,1), (1+s+t)/(1+s) e^{-t} <= (1+t) e^{-t}, car (1+s)(1+t) - (1+s+t) = s t >= 0,
# et ce n'est pas une loi exponentielle (P(X > 1) = 2/e, P(X > 2) = 3/e^2 != (2/e)^2).
assert simplify((1 + s) * (1 + t) - (1 + s + t) - s * t) == 0 and 3 * math.exp(-2) != (2 * math.exp(-1))**2
# Trou 1 : leurres. λe^{-λ|x|} d'intégrale 2 ; λ^2 x e^{-λx} (gamma(2,λ)) vaut 0 en 0 ; λe^{λx} sur ]-∞,0] porté par ℝ-.
assert I(lam * exp(-lam * abs(x)), (x, -oo, oo)) == 2
assert simplify(I(lam**2 * x * exp(-lam * x), (x, 0, oo))) == 1  # densité, mais de gamma(2,λ) : P(X>t) = (1+λt)e^{-λt}
assert simplify(I(lam**2 * x * exp(-lam * x), (x, t, oo)) - (1 + lam * t) * exp(-lam * t)) == 0
# Trou 2 : leurres faux (λ = 1) : 1 - e^{-t} vaut 0 en 0 ; 1/(1+t) != e^{-1} en t = 1 ; e^{-t^2} != e^{-t} en t = 2.
assert 1 - math.exp(0) == 0 and 1 / 2 != math.exp(-1) and math.exp(-4) != math.exp(-2)
# Affirmations : max de N(0, 1/100) ; ⌊X⌋+1 géométrique ; X^2 pour U([0,1]) ; 1/X Cauchy ; cX.
assert abs(10 / math.sqrt(2 * math.pi) - 3.99) < 0.005
k = Symbol("k", positive=True, integer=True)
assert simplify(exp(-lam * (k - 1)) - exp(-lam * k) - exp(-lam)**(k - 1) * (1 - exp(-lam))) == 0
for yv in [R(1, 3), 2, 7]:  # y > 0 : 1 - arctan(1/y)/pi = 1/2 + arctan(y)/pi
    assert egal(1 - atan(1 / yv) / pi, R(1, 2) + atan(yv) / pi)
for yv in [-R(1, 3), -2, -7]:  # y < 0 : -arctan(1/y)/pi = 1/2 + arctan(y)/pi
    assert egal(-atan(1 / yv) / pi, R(1, 2) + atan(yv) / pi)
assert math.exp(-0.5) != math.exp(-2)

# --- lois-usuelles-ex1 et variantes (absence de mémoire).
choix("lois-usuelles-ex1", 0, exp(-R(8, 2)) / exp(-R(5, 2)),
      {r"$e^{-3/2}$": exp(-R(3, 2)), r"$e^{-4}$": exp(-4), r"$e^{-5/2}$": exp(-R(5, 2)), r"$1-e^{-3/2}$": 1 - exp(-R(3, 2))})
choix("lois-usuelles-ex1", 1, exp(-R(10, 3)) / exp(-R(4, 3)),
      {r"$e^{-2}$": exp(-2), r"$e^{-10/3}$": exp(-R(10, 3)), r"$e^{-4/3}$": exp(-R(4, 3)), r"$1-e^{-2}$": 1 - exp(-2)})
choix("lois-usuelles-ex1", 2, exp(-R(6, 4)) / exp(-R(2, 4)),
      {r"$e^{-1}$": exp(-1), r"$e^{-3/2}$": exp(-R(3, 2)), r"$e^{-1/2}$": exp(-R(1, 2)), r"$1-e^{-1}$": 1 - exp(-1)})
assert round(math.exp(-1.5), 2) == 0.22 and round(math.exp(-2), 2) == 0.14 and round(math.exp(-1), 2) == 0.37

# --- lois-usuelles-ex2 et variantes (Cauchy).
Fcau = lambda v: R(1, 2) + atan(v) / pi
choix("lois-usuelles-ex2", 0, 1 - Fcau(sqrt(3)), {r"$\frac16$": R(1, 6), r"$\frac13$": R(1, 3), r"$\frac14$": R(1, 4),
                                                  r"$\frac56$": R(5, 6)})
choix("lois-usuelles-ex2", 1, Fcau(1 / sqrt(3)) - Fcau(-1 / sqrt(3)),
      {r"$\frac13$": R(1, 3), r"$\frac16$": R(1, 6), r"$\frac23$": R(2, 3), r"$\frac12$": R(1, 2)})
choix("lois-usuelles-ex2", 2, Fcau(sqrt(3)) - Fcau(-1),
      {r"$\frac7{12}$": R(7, 12), r"$\frac1{12}$": R(1, 12), r"$\frac5{12}$": R(5, 12), r"$\frac34$": R(3, 4)})
assert simplify(I(1 / (pi * (1 + t**2)), (t, -oo, x)) - Fcau(x)) == 0

# --- lois-usuelles-ex3 et variantes (lois Gamma).
for aa, ll, cc in [(3, 2, 4), (2, 3, 9), (4, 1, R(1, 6))]:
    assert ll**aa / gamma(aa) == cc and I(cc * x**(aa - 1) * exp(-ll * x), (x, 0, oo)) == 1
queue = lambda aa, ll, cc, tt: I(cc * x**(aa - 1) * exp(-ll * x), (x, tt, oo))
assert simplify(I(4 * x**2 * exp(-2 * x), (x, t, oo)) - (2 * t**2 + 2 * t + 1) * exp(-2 * t)) == 0
assert simplify(I(9 * x * exp(-3 * x), (x, t, oo)) - (3 * t + 1) * exp(-3 * t)) == 0
assert simplify(I(x**3 * exp(-x) / 6, (x, t, oo)) - (t**3 / 6 + t**2 / 2 + t + 1) * exp(-t)) == 0
choix("lois-usuelles-ex3", 0, queue(3, 2, 4, 1), {r"$5e^{-2}$": 5 * exp(-2), r"$e^{-2}$": exp(-2), r"$3e^{-2}$": 3 * exp(-2),
                                                  r"$1-5e^{-2}$": 1 - 5 * exp(-2)})
choix("lois-usuelles-ex3", 1, queue(2, 3, 9, 1), {r"$4e^{-3}$": 4 * exp(-3), r"$e^{-3}$": exp(-3), r"$3e^{-3}$": 3 * exp(-3),
                                                  r"$1-4e^{-3}$": 1 - 4 * exp(-3)})
choix("lois-usuelles-ex3", 2, queue(4, 1, R(1, 6), 2),
      {r"$\frac{19}3e^{-2}$": R(19, 3) * exp(-2), r"$5e^{-2}$": 5 * exp(-2), r"$\frac43e^{-2}$": R(4, 3) * exp(-2),
       r"$1-\frac{19}3e^{-2}$": 1 - R(19, 3) * exp(-2)})
assert round(5 * math.exp(-2), 2) == 0.68 and round(4 * math.exp(-3), 2) == 0.20 and round(19 / 3 * math.exp(-2), 2) == 0.86

# --- lois-usuelles-ex4 et variantes : N = ⌊X⌋, R = X - N.
q_ = exp(-lam)
n_ = Symbol("n", nonnegative=True, integer=True)
r_ = Symbol("r", positive=True)
assert simplify(I(lam * exp(-lam * x), (x, n_, n_ + 1)) - q_**n_ * (1 - q_)) == 0
assert simplify(I(lam * exp(-lam * x), (x, n_, n_ + r_)) - q_**n_ * (1 - exp(-lam * r_))) == 0
assert simplify(I(lam * exp(-lam * r_) / (1 - exp(-lam)), (r_, 0, 1))) == 1  # densité de R
PN = lambda ll, nn: exp(-ll * nn) * (1 - exp(-ll))
choix("lois-usuelles-ex4", 0, PN(log(2), 2), {r"$\frac18$": R(1, 8), r"$\frac14$": R(1, 4), r"$\frac34$": R(3, 4),
                                              r"$\frac1{16}$": R(1, 16)})
choix("lois-usuelles-ex4", 1, PN(log(3), 1), {r"$\frac29$": R(2, 9), r"$\frac13$": R(1, 3), r"$\frac23$": R(2, 3),
                                              r"$\frac19$": R(1, 9)})
choix("lois-usuelles-ex4", 2, PN(log(R(3, 2)), 2), {r"$\frac4{27}$": R(4, 27), r"$\frac49$": R(4, 9), r"$\frac59$": R(5, 9),
                                                    r"$\frac8{27}$": R(8, 27)})

# --- lois-usuelles-ex5 : phare. Q1 : dérivée de 1/2 + arctan(x/a)/pi. Q3 : tan(2θ) et angle ramené uniforme.
assert simplify(diff(R(1, 2) + atan(x / a) / pi, x) - a / (pi * (a**2 + x**2))) == 0
th = Symbol("theta", real=True)
assert all(abs(2 * math.tan(w) / (1 - math.tan(w)**2) - math.tan(2 * w)) < 1e-9 for w in [-1.2, -0.5, 0.1, 0.3, 1.1])
for tv in [-1.4, -0.7, -0.1, 0, 0.2, 0.9, 1.5]:
    L1 = tv / 2 + math.pi / 4
    L2 = math.pi / 4 if tv >= 0 else tv / 2 + math.pi / 4
    L3 = tv / 2 if tv >= 0 else 0
    # mesure exacte de {θ ∈ ]-π/2, π/2[ : θ' <= t} par discrétisation fine
    N = 200000
    pas = math.pi / N
    cpt = 0
    for i in range(N):
        tt = -math.pi / 2 + (i + 0.5) * pas
        tp = 2 * tt if abs(tt) <= math.pi / 4 else (2 * tt - math.pi if tt > 0 else 2 * tt + math.pi)
        cpt += tp <= tv
    assert abs(cpt * pas - (L1 + L2 + L3)) < 1e-4 and abs(L1 + L2 + L3 - (tv + math.pi / 2)) < 1e-12

# ===========================================================================================
# Carte proba-l3-densites-esperance
# ===========================================================================================
# Lois usuelles : moments.
aa_, bb_ = symbols("a_ b_", real=True)
EU = I(x / (bb_ - aa_), (x, aa_, bb_))
VU = simplify(I(x**2 / (bb_ - aa_), (x, aa_, bb_)) - EU**2)
assert simplify(EU - (aa_ + bb_) / 2) == 0 and simplify(VU - (bb_ - aa_)**2 / 12) == 0
EG = I(x * lam**a / gamma(a) * x**(a - 1) * exp(-lam * x), (x, 0, oo))
EG2 = I(x**2 * lam**a / gamma(a) * x**(a - 1) * exp(-lam * x), (x, 0, oo))
assert simplify(EG - a / lam) == 0 and simplify(EG2 - EG**2 - a / lam**2) == 0
phi = exp(-x**2 / 2) / sqrt(2 * pi)
assert I(x * phi, (x, -oo, oo)) == 0 and I(x**2 * phi, (x, -oo, oo)) == 1
assert I(x**2 * exp(-x**2 / 2), (x, -oo, oo)) == sqrt(2 * pi)
# Cauchy : ∫_0^A x/(π(1+x^2)) = ln(1+A^2)/(2π) -> +∞.
A = Symbol("A", positive=True)
assert simplify(I(x / (pi * (1 + x**2)), (x, 0, A)) - log(1 + A**2) / (2 * pi)) == 0
# Contre-exemples : 2/x^3 (E = 2, pas de moment d'ordre 2) ; x sin x contre 1/x^2 (semi-convergence) ; e^X pour E(1).
assert I(2 / x**2, (x, 1, oo)) == 2 and I(2 / x, (x, 1, oo)) == oo
assert abs(float(mp.quadosc(lambda v: mp.sin(v) / v, [1, mp.inf], omega=1)) - float(mp.pi / 2 - mp.si(1))) < 1e-15
assert I(sin(x)**2 / x, (x, 1, oo)) == oo  # |sin x|/x >= sin^2 x / x : pas de convergence absolue
assert I(exp(x) * exp(-x), (x, 0, oo)) == oo
# Exemple : E(e^X) = e - 1 pour U([0,1]) ; ∫ dx/x diverge ; E(e^{tX}) = λ/(λ-t) pour t < λ.
assert I(exp(x), (x, 0, 1)) == E - 1 and I(1 / x, (x, 0, 1)) == oo
assert simplify(I(lam * exp((2 - lam) * x), (x, 0, oo)).subs(lam, 5) - R(5, 3)) == 0
# Trou 2 (φ = x^2, U([0,1])) : ∫φf = 1/3 ; leurre φ(E X) = 1/4 ; (φ = id) ∫φ(f) = 1 != E X = 1/2 ; (φ = 1) ∫φ diverge.
assert I(x**2, (x, 0, 1)) == R(1, 3) != R(1, 4)
# Trou 1 / pièges : E(1) et φ = exp (∫|x| f converge, φ(X) sans espérance) ; φ = 1 (∫|φ| diverge, φ(X) = 1 intégrable) ;
# f = 1/x^2 sur [1,∞[, φ = x : φ f -> 0 mais ∫ 1/x diverge.
assert I(x * exp(-x), (x, 0, oo)) == 1 and I(1 / x, (x, 1, oo)) == oo
# Affirmations : E(X^2) = 2/λ^2 pour E(λ) ; formule de la queue (testée avec gamma(2,1) : ∫(1+t)e^{-t} = 2 = E(X)).
assert I(x**2 * lam * exp(-lam * x), (x, 0, oo)) == 2 / lam**2
assert I((1 + t) * exp(-t), (t, 0, oo)) == 2 == I(x * x * exp(-x), (x, 0, oo))

# --- esperance-ex1 et variantes.
for dens, E1, E2, V, ver, leur in [(3 * x**2, R(3, 4), R(3, 5), R(3, 80), 0, (R(3, 5), R(3, 4), R(9, 16))),
                                   (2 * x, R(2, 3), R(1, 2), R(1, 18), 1, (R(1, 2), R(2, 3), R(4, 9))),
                                   (4 * x**3, R(4, 5), R(2, 3), R(2, 75), 2, (R(2, 3), R(4, 5), R(16, 25)))]:
    assert I(dens, (x, 0, 1)) == 1 and I(x * dens, (x, 0, 1)) == E1 and I(x**2 * dens, (x, 0, 1)) == E2
    assert E2 - E1**2 == V
choix("esperance-ex1", 0, R(3, 80), {r"$\frac3{80}$": R(3, 80), r"$\frac35$": R(3, 5), r"$\frac34$": R(3, 4),
                                     r"$\frac9{16}$": R(9, 16)})
choix("esperance-ex1", 1, R(1, 18), {r"$\frac1{18}$": R(1, 18), r"$\frac12$": R(1, 2), r"$\frac23$": R(2, 3),
                                     r"$\frac49$": R(4, 9)})
choix("esperance-ex1", 2, R(2, 75), {r"$\frac2{75}$": R(2, 75), r"$\frac23$": R(2, 3), r"$\frac45$": R(4, 5),
                                     r"$\frac{16}{25}$": R(16, 25)})

# --- esperance-ex2 et variantes : E(e^{tX}) = λ/(λ - t) pour t < λ, divergence pour t >= λ.
mgf = lambda ll, tt: I(ll * exp((tt - ll) * x), (x, 0, oo))
assert mgf(2, 2) == oo and mgf(3, 3) == oo and mgf(2, R(5, 2)) == oo
choix("esperance-ex2", 0, mgf(2, 1), {r"$2$": 2, r"$e^{1/2}$": exp(R(1, 2)), r"$\frac12$": R(1, 2), r"$\frac23$": R(2, 3)})
choix("esperance-ex2", 1, mgf(3, 1), {r"$\frac32$": R(3, 2), r"$e^{1/3}$": exp(R(1, 3)), r"$\frac34$": R(3, 4), r"$3$": 3})
choix("esperance-ex2", 2, mgf(3, 2), {r"$3$": 3, r"$e^{2/3}$": exp(R(2, 3)), r"$\frac35$": R(3, 5), r"$\frac32$": R(3, 2)})

# --- esperance-ex3 et variantes : loi de Pareto.
rr = Symbol("r", positive=True)
assert I(a / x**(a + 1), (x, 1, oo)) == 1
assert I(x**3 * 3 / x**4, (x, 1, oo)) == oo  # a = 3, r = 3 : X^3 sans espérance
PV = lambda aa: I(x**2 * aa / x**(aa + 1), (x, 1, oo)) - I(x * aa / x**(aa + 1), (x, 1, oo))**2
NON = 99  # valeur sentinelle pour « X n'admet pas de variance »
choix("esperance-ex3", 0, PV(3), {r"$\frac34$": R(3, 4), r"$3$": 3, r"$\frac94$": R(9, 4), r"$X$ n'admet pas de variance": NON})
choix("esperance-ex3", 1, PV(4), {r"$\frac29$": R(2, 9), r"$2$": 2, r"$\frac{16}9$": R(16, 9), r"$X$ n'admet pas de variance": NON})
choix("esperance-ex3", 2, PV(R(5, 2)), {r"$\frac{20}9$": R(20, 9), r"$5$": 5, r"$\frac{25}9$": R(25, 9),
                                        r"$X$ n'admet pas de variance": NON})

# --- esperance-ex4 et variantes : E(X) = ∫ P(X > t) dt ; les survies données sont continues en 0 (valeur 1).
for Gs, Ev, ver in [(1 / (1 + t**2), pi / 2, 0), (1 / (1 + t)**3, R(1, 2), 1), (exp(-t**2), sqrt(pi) / 2, 2)]:
    assert Gs.subs(t, 0) == 1 and limit(Gs, t, oo) == 0
    fd = -diff(Gs, t)
    assert I(t * fd, (t, 0, oo)) == Ev == I(Gs, (t, 0, oo))
choix("esperance-ex4", 0, pi / 2, {r"$\frac\pi2$": pi / 2, r"$\pi$": pi, r"$1$": 1, r"$\frac\pi4$": pi / 4})
choix("esperance-ex4", 1, R(1, 2), {r"$\frac12$": R(1, 2), r"$\frac13$": R(1, 3), r"$1$": 1, r"$\frac32$": R(3, 2)})
choix("esperance-ex4", 2, sqrt(pi) / 2, {r"$\frac{\sqrt\pi}2$": sqrt(pi) / 2, r"$\sqrt\pi$": sqrt(pi), r"$\frac12$": R(1, 2),
                                         r"$\frac\pi2$": pi / 2})

# --- esperance-ex5 et variantes : φ(a) = E|X - a| minimal en la médiane.
aS = Symbol("aS", positive=True)
phiE = I((aS - x) * exp(-x), (x, 0, aS)) + I((x - aS) * exp(-x), (x, aS, oo))
assert simplify(phiE - (aS - 1 + 2 * exp(-aS))) == 0
assert simplify(phiE.subs(aS, log(2)) - log(2)) == 0 and simplify(phiE.subs(aS, 1) - 2 / E) == 0
assert round(math.log(2), 3) == 0.693 and round(2 / math.e, 3) == 0.736
choix("esperance-ex5", 0, log(2), {r"$\ln2$": log(2), r"$\frac2e$": 2 / E, r"$1$": 1, r"$\frac12$": R(1, 2)})
phiU = lambda av: I(abs(x - av), (x, 0, 1))
assert phiU(R(1, 2)) == R(1, 4)
choix("esperance-ex5", 1, phiU(R(1, 2)), {r"$\frac14$": R(1, 4), r"$\frac12$": R(1, 2), r"$\frac13$": R(1, 3),
                                          r"$\frac18$": R(1, 8)})
phiT = lambda av: I((x - av) * 2 * x, (x, av, 1)) + I((av - x) * 2 * x, (x, 0, av))
assert simplify(phiT(aS) - (R(2, 3) * aS**3 - aS + R(2, 3))) == 0
assert simplify(phiT(1 / sqrt(2)) - (2 - sqrt(2)) / 3) == 0 and phiT(R(2, 3)) == R(16, 81)
assert (2 - math.sqrt(2)) / 3 < 16 / 81
choix("esperance-ex5", 2, phiT(1 / sqrt(2)), {r"$\frac{2-\sqrt2}3$": (2 - sqrt(2)) / 3, r"$\frac{16}{81}$": R(16, 81),
                                              r"$\frac{\sqrt2}2$": sqrt(2) / 2, r"$\frac13$": R(1, 3)})

# ===========================================================================================
# Carte proba-l3-densites-vecteurs
# ===========================================================================================
# Exemple : densité e^{-y} sur 0 < x < y.
assert I(I(exp(-y), (y, x, oo)), (x, 0, oo)) == 1
assert I(exp(-y), (y, x, oo)) == exp(-x) and I(exp(-y), (x, 0, y)) == y * exp(-y)
assert I(I(exp(-y), (y, 2 * x, oo)), (x, 0, oo)) == R(1, 2)
# Disque unité : marginale, covariance nulle (E(XY) = 0 en polaires), 9/16 + 9/16 > 1.
assert I(1 / pi, (y, -sqrt(1 - x**2), sqrt(1 - x**2))) == 2 * sqrt(1 - x**2) / pi
rho, th2 = symbols("rho th2", nonnegative=True)
assert I(I(rho * cos(th2) * rho * sin(th2) / pi * rho, (rho, 0, 1)), (th2, 0, 2 * pi)) == 0
assert R(9, 16) + R(9, 16) > 1
# Densité 1 + (2x-1)(2y-1) : marginales uniformes, P(X <= 1/2, Y <= 1/2) = 5/16.
fF = 1 + (2 * x - 1) * (2 * y - 1)
assert I(fF, (y, 0, 1)) == 1 and I(I(fF, (x, 0, R(1, 2))), (y, 0, R(1, 2))) == R(5, 16)
# Contre-exemple Y = X : P(X <= 1/2, Y <= 1/2) = 1/2, produit 1/4.
# Trou 1, leurres : X, Y indépendantes de lois E(1) et E(2) : ∫ f(y,x) dy = f_Y(x) = 2e^{-2x} != e^{-x} ;
# X, Y indépendantes de loi E(1) : f(x,x) = e^{-2x} != e^{-x} ; X, Y uniformes : ∫_{-∞}^x f(x,y) dy = x != 1 en x = 1/2.
assert I(exp(-y) * 2 * exp(-2 * x), (y, 0, oo)) == 2 * exp(-2 * x)
# Trou 2, leurres : f_X(x) + f_Y(y) et f_X(x) f_Y(x) ont une intégrale double infinie (exemple : uniformes sur [0,1]) ;
# f_X(x) f_Y(x - y) pour E(1) indépendantes est nulle près de (1, 2), où le produit vaut e^{-3} > 0.
assert (1 - 2) < 0  # f_Y(x - y) = 0 en (1, 2) puisque x - y = -1 < 0

# --- vecteurs-ex1 et variantes : constante et marginale.
for dens, c, fx, ver in [(x + y, 1, x + R(1, 2), 0), (x * y, 4, 2 * x, 1), (x**2 + y, R(6, 5), (6 * x**2 + 3) / 5, 2)]:
    assert 1 / I(I(dens, (x, 0, 1)), (y, 0, 1)) == c
    assert simplify(I(c * dens, (y, 0, 1)) - fx) == 0
for fx_l in [x, 2 * x, x + 1]:  # leurres de l'exercice : pas la marginale x + 1/2
    assert simplify(fx_l - (x + R(1, 2))) != 0
assert I(x**2 + R(1, 2), (x, 0, 1)) == R(5, 6) and I(6 * x**2 / 5, (x, 0, 1)) == R(2, 5) and I((6 * x**2 + 6) / 5, (x, 0, 1)) == R(8, 5)
assert I(I(x + y, (x, 0, R(1, 2))), (y, 0, R(1, 2))) == R(1, 8) and (R(1, 8) + R(1, 4))**2 == R(9, 64)
# Les propositions sont des fonctions : on les compare par leur valeur en x = 1/3.
choix("vecteurs-ex1", 0, R(1, 3) + R(1, 2), {r"$x\mapsto x+\frac12$": R(1, 3) + R(1, 2), r"$x\mapsto x$": R(1, 3),
                                             r"$x\mapsto2x$": R(2, 3), r"$x\mapsto x+1$": R(4, 3)})
choix("vecteurs-ex1", 1, R(2, 3), {r"$x\mapsto2x$": R(2, 3), r"$x\mapsto4x$": R(4, 3), r"$x\mapsto x$": R(1, 3),
                                   r"$x\mapsto x+\frac12$": R(5, 6)})
X3 = R(1, 3)
choix("vecteurs-ex1", 2, (6 * X3**2 + 3) / 5, {r"$x\mapsto\frac{6x^2+3}5$": (6 * X3**2 + 3) / 5,
                                               r"$x\mapsto x^2+\frac12$": X3**2 + R(1, 2),
                                               r"$x\mapsto\frac{6x^2}5$": 6 * X3**2 / 5,
                                               r"$x\mapsto\frac{6x^2+6}5$": (6 * X3**2 + 6) / 5})

# --- vecteurs-ex2 et variantes : P(X < Y) = λ/(λ+μ).
PXY = lambda l1, l2: I(I(l1 * l2 * exp(-l1 * x - l2 * y), (y, x, oo)), (x, 0, oo))
choix("vecteurs-ex2", 0, PXY(2, 3), {r"$\frac25$": R(2, 5), r"$\frac35$": R(3, 5), r"$\frac12$": R(1, 2), r"$\frac23$": R(2, 3)})
choix("vecteurs-ex2", 1, PXY(1, 4), {r"$\frac15$": R(1, 5), r"$\frac45$": R(4, 5), r"$\frac12$": R(1, 2), r"$\frac14$": R(1, 4)})
choix("vecteurs-ex2", 2, PXY(3, 1), {r"$\frac34$": R(3, 4), r"$\frac14$": R(1, 4), r"$\frac12$": R(1, 2), r"$\frac13$": R(1, 3)})

# --- vecteurs-ex3 et variantes : disque uniforme.
fXd = 2 * sqrt(1 - x**2) / pi
assert I(fXd, (x, -1, 1)) == 1
assert simplify(I(sqrt(1 - x**2), (x, R(1, 2), 1)) - (pi / 6 - sqrt(3) / 8)) == 0
choix("vecteurs-ex3", 0, I(fXd, (x, R(1, 2), 1)),
      {r"$\frac13-\frac{\sqrt3}{4\pi}$": R(1, 3) - sqrt(3) / (4 * pi), r"$\frac14$": R(1, 4), r"$\frac13$": R(1, 3),
       r"$\frac16-\frac{\sqrt3}{8\pi}$": R(1, 6) - sqrt(3) / (8 * pi)})
choix("vecteurs-ex3", 1, I(fXd, (x, -R(1, 2), R(1, 2))),
      {r"$\frac13+\frac{\sqrt3}{2\pi}$": R(1, 3) + sqrt(3) / (2 * pi), r"$\frac12$": R(1, 2),
       r"$\frac23-\frac{\sqrt3}{4\pi}$": R(2, 3) - sqrt(3) / (4 * pi), r"$\frac{\sqrt3}{2\pi}$": sqrt(3) / (2 * pi)})
choix("vecteurs-ex3", 2, I(fXd, (x, sqrt(2) / 2, 1)),
      {r"$\frac14-\frac1{2\pi}$": R(1, 4) - 1 / (2 * pi), r"$\frac{2-\sqrt2}4$": (2 - sqrt(2)) / 4, r"$\frac14$": R(1, 4),
       r"$\frac12-\frac1\pi$": R(1, 2) - 1 / pi})
assert round(float(R(1, 3) - sqrt(3) / (4 * pi)), 3) == 0.196 and round(float(R(1, 3) + sqrt(3) / (2 * pi)), 3) == 0.609
assert round(float(R(1, 4) - 1 / (2 * pi)), 3) == 0.091

# --- vecteurs-ex4 et variantes : min et max de n uniformes, E(M - m) = (n-1)/(n+1).
for nn, rep in [(2, R(1, 3)), (3, R(1, 2)), (4, R(3, 5))]:
    EM = I(x * nn * x**(nn - 1), (x, 0, 1))
    Em = I(x * nn * (1 - x)**(nn - 1), (x, 0, 1))
    assert EM == R(nn, nn + 1) and Em == R(1, nn + 1) and EM - Em == rep
choix("vecteurs-ex4", 0, R(1, 3), {r"$\frac13$": R(1, 3), r"$\frac23$": R(2, 3), r"$\frac12$": R(1, 2), r"$\frac16$": R(1, 6)})
choix("vecteurs-ex4", 1, R(1, 2), {r"$\frac12$": R(1, 2), r"$\frac34$": R(3, 4), r"$\frac14$": R(1, 4), r"$\frac13$": R(1, 3)})
choix("vecteurs-ex4", 2, R(3, 5), {r"$\frac35$": R(3, 5), r"$\frac45$": R(4, 5), r"$\frac15$": R(1, 5), r"$\frac12$": R(1, 2)})

# --- vecteurs-ex5 et variantes : aiguille de Buffon, p = 2ℓ/(πd) (d = 1).
buffon = lambda l: I(I(4 / pi, (x, 0, l / 2 * sin(th))), (th, 0, pi / 2))
choix("vecteurs-ex5", 0, buffon(1), {r"$\frac2\pi$": 2 / pi, r"$\frac1\pi$": 1 / pi, r"$\frac12$": R(1, 2), r"$\frac\pi4$": pi / 4})
choix("vecteurs-ex5", 1, buffon(R(1, 2)), {r"$\frac1\pi$": 1 / pi, r"$\frac2\pi$": 2 / pi, r"$\frac14$": R(1, 4),
                                           r"$\frac\pi8$": pi / 8})
choix("vecteurs-ex5", 2, buffon(R(3, 4)), {r"$\frac3{2\pi}$": 3 / (2 * pi), r"$\frac3{4\pi}$": 3 / (4 * pi),
                                           r"$\frac38$": R(3, 8), r"$\frac{3\pi}{16}$": 3 * pi / 16})
assert round(2 / math.pi, 3) == 0.637 and round(1 / math.pi, 3) == 0.318 and round(1.5 / math.pi, 3) == 0.477

# ===========================================================================================
# Carte proba-l3-densites-somme
# ===========================================================================================
# Exemple : U([0,1]) * U([0,1]) triangulaire ; E(λ) * E(λ) = λ^2 z e^{-λz}.
zp = Symbol("zp", positive=True)
assert I(lam * exp(-lam * x) * lam * exp(-lam * (zp - x)), (x, 0, zp)) == lam**2 * zp * exp(-lam * zp)
assert I(z, (z, 0, R(1, 2))) == R(1, 8)
assert simplify(I(z**2 * lam**2 * z * exp(-lam * z), (z, 0, oo)) - I(z * lam**2 * z * exp(-lam * z), (z, 0, oo))**2 - 2 / lam**2) == 0
# Lois Gamma : convolution proportionnelle à z^{a+b-1}e^{-z}, et B(a,b) = Γ(a)Γ(b)/Γ(a+b) (valeurs testées).
for aa, bb in [(R(1, 2), R(3, 2)), (2, 3), (R(1, 2), R(1, 2)), (R(3, 2), R(5, 2))]:
    Bab = integrate(u**(aa - 1) * (1 - u)**(bb - 1), (u, 0, 1))
    assert simplify(Bab - gamma(aa) * gamma(bb) / gamma(aa + bb)) == 0
    hz = integrate(x**(aa - 1) * exp(-x) / gamma(aa) * (zp - x)**(bb - 1) * exp(-(zp - x)) / gamma(bb), (x, 0, zp))
    assert simplify(hz - zp**(aa + bb - 1) * exp(-zp) / gamma(aa + bb)) == 0
# Contre-exemples : Y = X uniforme (densité 1/2 de 2X contre triangle h(1) = 1) ; disque : x + y <= √2 ;
# f_X * f_X > 0 sur ]-2,2[ (valeur en 1,9) ; E(1) + E(2) : 2(e^{-z} - e^{-2z}) contre 2z e^{-√2 z}.
hdisque = mp.quad(lambda v: 2 / mp.pi * mp.sqrt(1 - v**2) * 2 / mp.pi * mp.sqrt(max(0, 1 - (1.9 - v)**2)), [0.9, 1])
assert hdisque > 0
assert simplify(I(exp(-x) * 2 * exp(-2 * (zp - x)), (x, 0, zp)) - 2 * (exp(-zp) - exp(-2 * zp))) == 0
assert abs(2 * (math.exp(-1) - math.exp(-2)) - 2 * math.exp(-math.sqrt(2))) > 0.01
assert limit(2 * (exp(-z) - exp(-2 * z)) / (2 * z), z, 0) == 1
# Trou 1, leurres : f + g d'intégrale 2 ; f g pour U([0,1]) : 1_[0,1] (vaut 0 en 3/2, h(3/2) = 1/2) ;
# ∫ f(x) g(z+x) dx pour E(1) : e^{-z}/2 en z > 0, au lieu de z e^{-z} (en z = 1 : e^{-1}/2 != e^{-1}).
assert I(exp(-x) * exp(-(1 + x)), (x, 0, oo)) == exp(-1) / 2
# Trou 2, leurres (X, Y indépendantes de loi E(1)) : V(X+Y) = 2 ; 2 + 2·1·1 = 4 ; (1+1)^2 = 4 ; 1·1 = 1.
assert 1 + 1 == 2 != 4 and 2 != 1
# Affirmations : V(X - Y) = 1/6 pour U([0,1]) ; Cauchy (X+Y)/2 ; E(X+Y) = 2/λ contre 1/(2λ).
assert 2 * R(1, 12) == R(1, 6)
for zv in [0, 0.5, 1, 3, -2.2]:
    hc = mp.quad(lambda v: 1 / (mp.pi**2 * (1 + v**2) * (1 + (zv - v)**2)), [-mp.inf] + sorted({0, zv}) + [mp.inf])
    assert abs(hc - 2 / (mp.pi * (4 + zv**2))) < 1e-12
assert simplify(2 * 2 / (pi * (4 + (2 * u)**2)) - 1 / (pi * (1 + u**2))) == 0

# --- somme-ex1 et variantes.
choix("somme-ex1", 0, 4 * R(1, 3) + 1, {r"$\frac73$": R(7, 3), r"$\frac13$": R(1, 3), r"$\frac53$": R(5, 3), r"$\frac43$": R(4, 3)})
choix("somme-ex1", 1, R(9, 12) + 4 * R(1, 4), {r"$\frac74$": R(7, 4), r"$-\frac14$": -R(1, 4), r"$\frac54$": R(5, 4), r"$1$": 1})
choix("somme-ex1", 2, 9 * R(4, 12) + 4, {r"$7$": 7, r"$5$": 5, r"$\frac{13}3$": R(13, 3), r"$-1$": -1})
assert I(x**2 / 2, (x, 0, 2)) - I(x / 2, (x, 0, 2))**2 == R(1, 3)

# --- somme-ex2 et variantes : covariance sur un triangle.
def cov_tri(dens, haut):
    EX = I(I(x * dens, (x, 0, y)), (y, 0, haut))
    EY = I(I(y * dens, (x, 0, y)), (y, 0, haut))
    EXY = I(I(x * y * dens, (x, 0, y)), (y, 0, haut))
    assert I(I(dens, (x, 0, y)), (y, 0, haut)) == 1
    return EX, EY, EXY, EXY - EX * EY


assert cov_tri(2, 1) == (R(1, 3), R(2, 3), R(1, 4), R(1, 36))
assert cov_tri(8 * x * y, 1) == (R(8, 15), R(4, 5), R(4, 9), R(4, 225))
assert cov_tri(R(1, 2), 2) == (R(2, 3), R(4, 3), 1, R(1, 9))
choix("somme-ex2", 0, R(1, 36), {r"$\frac1{36}$": R(1, 36), r"$\frac14$": R(1, 4), r"$-\frac1{36}$": -R(1, 36), r"$\frac29$": R(2, 9)})
choix("somme-ex2", 1, R(4, 225), {r"$\frac4{225}$": R(4, 225), r"$\frac49$": R(4, 9), r"$\frac{32}{75}$": R(32, 75),
                                  r"$-\frac4{225}$": -R(4, 225)})
choix("somme-ex2", 2, R(1, 9), {r"$\frac19$": R(1, 9), r"$1$": 1, r"$\frac89$": R(8, 9), r"$-\frac19$": -R(1, 9)})

# --- somme-ex3 et variantes : loi triangulaire.
def htri(lo, hi):
    """P(lo <= X + Y <= hi) pour la densité triangulaire (morceaux vides ignorés)."""
    tot = 0
    for g, p, q in [(z, 0, 1), (2 - z, 1, 2)]:
        if max(lo, p) < min(hi, q):
            tot += I(g, (z, max(lo, p), min(hi, q)))
    return tot


assert htri(0, 2) == 1
choix("somme-ex3", 0, htri(0, R(3, 2)), {r"$\frac78$": R(7, 8), r"$\frac34$": R(3, 4), r"$\frac18$": R(1, 8), r"$\frac12$": R(1, 2)})
choix("somme-ex3", 1, htri(0, R(1, 2)), {r"$\frac18$": R(1, 8), r"$\frac14$": R(1, 4), r"$\frac78$": R(7, 8), r"$\frac12$": R(1, 2)})
choix("somme-ex3", 2, htri(R(1, 2), R(3, 2)), {r"$\frac34$": R(3, 4), r"$\frac12$": R(1, 2), r"$\frac14$": R(1, 4),
                                               r"$\frac78$": R(7, 8)})

# --- somme-ex4 et variantes : E(λ) + E(μ), λ != μ.
def P_somme(l1, l2):
    h = I(l1 * exp(-l1 * x) * l2 * exp(-l2 * (zp - x)), (x, 0, zp))
    assert simplify(h - l1 * l2 / (l2 - l1) * (exp(-l1 * zp) - exp(-l2 * zp))) == 0
    return I(h, (zp, 0, 1))


v = P_somme(1, 2)
assert simplify(v - (1 - exp(-1))**2) == 0
choix("somme-ex4", 0, v, {r"$\left(1-e^{-1}\right)^2$": (1 - exp(-1))**2, r"$1-e^{-3}$": 1 - exp(-3),
                          r"$\left(1-e^{-1}\right)\left(1-e^{-2}\right)$": (1 - exp(-1)) * (1 - exp(-2)),
                          r"$1-2e^{-1}$": 1 - 2 * exp(-1)})
choix("somme-ex4", 1, P_somme(1, 3), {r"$1-\frac32e^{-1}+\frac12e^{-3}$": 1 - R(3, 2) * exp(-1) + exp(-3) / 2,
                                      r"$1-e^{-4}$": 1 - exp(-4),
                                      r"$\left(1-e^{-1}\right)\left(1-e^{-3}\right)$": (1 - exp(-1)) * (1 - exp(-3)),
                                      r"$1-\frac32e^{-1}$": 1 - R(3, 2) * exp(-1)})
choix("somme-ex4", 2, P_somme(2, 3), {r"$1-3e^{-2}+2e^{-3}$": 1 - 3 * exp(-2) + 2 * exp(-3), r"$1-e^{-5}$": 1 - exp(-5),
                                      r"$\left(1-e^{-2}\right)\left(1-e^{-3}\right)$": (1 - exp(-2)) * (1 - exp(-3)),
                                      r"$1-3e^{-2}$": 1 - 3 * exp(-2)})
assert [round(float(P_somme(*p)), 2) for p in [(1, 2), (1, 3), (2, 3)]] == [0.40, 0.47, 0.69]

# --- somme-ex5 et variantes : fonction Bêta.
choix("somme-ex5", 0, I(sqrt(x * (1 - x)), (x, 0, 1)), {r"$\frac\pi8$": pi / 8, r"$\frac\pi4$": pi / 4, r"$\frac16$": R(1, 6),
                                                       r"$\frac\pi{16}$": pi / 16})
assert I(x * (1 - x), (x, 0, 1)) == R(1, 6) and gamma(R(3, 2))**2 / gamma(3) == pi / 8
choix("somme-ex5", 1, integrate(sqrt(x) * (1 - x)**R(3, 2), (x, 0, 1)),
      {r"$\frac\pi{16}$": pi / 16, r"$\frac\pi8$": pi / 8, r"$\frac{3\pi}{16}$": 3 * pi / 16, r"$\frac\pi{32}$": pi / 32})
choix("somme-ex5", 2, integrate(1 / sqrt(x * (1 - x)), (x, 0, 1)),
      {r"$\pi$": pi, r"$\frac\pi2$": pi / 2, r"$2$": 2, r"$\sqrt\pi$": sqrt(pi)})

# --- somme-ex6 : identité de décomposition, h(0) = 1/(2π), P(|(X+Y)/2| > 1) = 1/2.
zz = Symbol("zz", real=True, nonzero=True)
gauche = 1 / ((1 + x**2) * (1 + (zz - x)**2))
droite = ((1 + 2 * x / zz) / (1 + x**2) + (1 + 2 * (zz - x) / zz) / (1 + (zz - x)**2)) / (zz**2 + 4)
assert simplify(gauche - droite) == 0
assert I(1 / (pi**2 * (1 + x**2)**2), (x, -oo, oo)) == 1 / (2 * pi) == R(2, 4) / pi
assert 1 - 2 * atan(1) / pi == R(1, 2)
choix("somme-ex6", 0, R(1, 2), {r"$\frac12$": R(1, 2), r"$\frac14$": R(1, 4), r"$\frac2\pi\arctan\frac12$": 2 * atan(R(1, 2)) / pi,
                                r"$\frac1{2\pi}$": 1 / (2 * pi)})
# Résultat admis de la question 3 (densités a/(π(a^2+x^2)) et b/(π(b^2+x^2))), contrôlé numériquement pour (a, b) = (2, 1).
for zv in [0, 1.3, -4]:
    hc = mp.quad(lambda v: 2 / (mp.pi * (4 + v**2)) / (mp.pi * (1 + (zv - v)**2)), [-mp.inf] + sorted({0, zv}) + [mp.inf])
    assert abs(hc - 3 / (mp.pi * (9 + zv**2))) < 1e-12

# ===========================================================================================
# Carte proba-l3-densites-normales
# ===========================================================================================
s1, s2 = symbols("s1 s2", positive=True)
S2 = s1**2 + s2**2
assert simplify(x**2 / s1**2 + (z - x)**2 / s2**2 - (S2 / (s1**2 * s2**2) * (x - s1**2 * z / S2)**2 + z**2 / S2)) == 0
for (v1, v2) in [(1, 2), (R(1, 2), 3)]:  # convolution de N(0, v1^2) et N(0, v2^2) = N(0, v1^2 + v2^2)
    hz = integrate(exp(-x**2 / (2 * v1**2)) / (v1 * sqrt(2 * pi)) * exp(-(z - x)**2 / (2 * v2**2)) / (v2 * sqrt(2 * pi)), (x, -oo, oo))
    ss = sqrt(v1**2 + v2**2)
    assert simplify(hz - exp(-z**2 / (2 * ss**2)) / (ss * sqrt(2 * pi))) == 0
# Transformation affine avec a < 0 : densité de aX + b pour X ~ N(m, σ^2) (a = -2, b = 1, m = 3, σ = 2).
Fx = lambda v: Phi((v - 3) / 2)
Fa = lambda v: 1 - Fx((v - 1) / -2)
for v in [-9, -5, 0.3]:
    assert egal(Fa(v), Phi((v - (-2 * 3 + 1)) / 4))
# Contre-exemple « couple à densité » : marginale N(0,1) ; h(0) = 0 ; h(z) -> 0 en 0.
fxy = exp(-(x**2 + y**2) / 2) / pi
assert simplify(I(fxy.subs(x, R(7, 10)), (y, 0, oo)) - exp(-R(49, 200)) / sqrt(2 * pi)) == 0
hpos = lambda zv: mp.quad(lambda v: mp.exp(-(v**2 + (zv - v)**2) / 2) / mp.pi, [0, zv])
assert hpos(1e-8) < 1e-8 and hpos(1) > 0.1
# Exemple : Φ(2) ≈ 0,977 ; affirmation : 2Φ(1) - 1 ≈ 0,683 ; contre-exemple εX : P(εX <= y) = Φ(y).
assert round(float(Phi(2)), 3) == 0.977 and round(float(2 * Phi(1) - 1), 3) == 0.683
for yv in [-1.2, 0.4, 2]:
    assert egal(Phi(yv) / 2 + (1 - Phi(-yv)) / 2, Phi(yv))
# Trous : (X - m)/σ^2, X - m, σ(X - m) ont pour variance 1/σ^2, σ^2, σ^4 (≠ 1 pour σ = 2) ;
# (σ1 + σ2)^2, σ1 σ2, (σ1^2 + σ2^2)/2 diffèrent de σ1^2 + σ2^2 pour σ1, σ2 > 0.
assert simplify((s1 + s2)**2 - S2 - 2 * s1 * s2) == 0 and simplify(S2 - s1 * s2 - ((s1 - s2)**2 + s1 * s2)) == 0
assert [R(1, 4), 4, 16] == [R(1, 2**2), 2**2, 2**4]

# --- normales-ex1 et variantes : Φ.
choix("normales-ex1", 0, Phi(R(5 - 3, 2)), {r"$\Phi(1)$": Phi(1), r"$\Phi\left(\frac12\right)$": Phi(R(1, 2)),
                                           r"$\Phi(2)$": Phi(2), r"$1-\Phi(1)$": 1 - Phi(1)})
choix("normales-ex1", 1, 1 - Phi(R(2 + 1, 3)), {r"$1-\Phi(1)$": 1 - Phi(1), r"$\Phi(1)$": Phi(1),
                                               r"$1-\Phi\left(\frac13\right)$": 1 - Phi(R(1, 3)), r"$1-\Phi(3)$": 1 - Phi(3)})
choix("normales-ex1", 2, Phi((9 - 10) / R(1, 2)), {r"$1-\Phi(2)$": 1 - Phi(2), r"$\Phi(2)$": Phi(2),
                                                  r"$1-\Phi\left(\frac12\right)$": 1 - Phi(R(1, 2)), r"$1-\Phi(4)$": 1 - Phi(4)})
assert [round(float(w), 3) for w in [Phi(1), 1 - Phi(1), 1 - Phi(2)]] == [0.841, 0.159, 0.023]

# --- normales-ex2 et variantes : E|X| = σ sqrt(2/π).
Eabs = lambda sg: I(abs(x) * exp(-x**2 / (2 * sg**2)) / (sg * sqrt(2 * pi)), (x, -oo, oo))
choix("normales-ex2", 0, Eabs(2), {r"$2\sqrt{\frac2\pi}$": 2 * sqrt(2 / pi), r"$0$": 0, r"$\sqrt{\frac2\pi}$": sqrt(2 / pi),
                                   r"$4\sqrt{\frac2\pi}$": 4 * sqrt(2 / pi)})
choix("normales-ex2", 1, Eabs(3), {r"$3\sqrt{\frac2\pi}$": 3 * sqrt(2 / pi), r"$0$": 0, r"$9\sqrt{\frac2\pi}$": 9 * sqrt(2 / pi),
                                   r"$\frac13\sqrt{\frac2\pi}$": sqrt(2 / pi) / 3})
choix("normales-ex2", 2, Eabs(R(1, 2)), {r"$\frac1{\sqrt{2\pi}}$": 1 / sqrt(2 * pi), r"$0$": 0,
                                         r"$\frac14\sqrt{\frac2\pi}$": sqrt(2 / pi) / 4, r"$\sqrt{\frac2\pi}$": sqrt(2 / pi)})
assert [round(float(Eabs(w)), 2) for w in [2, 3, R(1, 2)]] == [1.60, 2.39, 0.40]

# --- normales-ex3 et variantes : P(X > Y) avec X - Y ~ N(m1 - m2, v1 + v2).
PXgY = lambda m1, v1, m2, v2: 1 - Phi((0 - (m1 - m2)) / sqrt(v1 + v2))
choix("normales-ex3", 0, PXgY(50, 1, 51, 4),
      {r"$1-\Phi\left(\frac1{\sqrt5}\right)$": 1 - Phi(1 / sqrt(5)), r"$1-\Phi\left(\frac1{\sqrt3}\right)$": 1 - Phi(1 / sqrt(3)),
       r"$1-\Phi\left(\frac13\right)$": 1 - Phi(R(1, 3)), r"$\Phi\left(\frac1{\sqrt5}\right)$": Phi(1 / sqrt(5))})
choix("normales-ex3", 1, PXgY(5, 9, 3, 16),
      {r"$\Phi\left(\frac25\right)$": Phi(R(2, 5)), r"$1-\Phi\left(\frac25\right)$": 1 - Phi(R(2, 5)),
       r"$\Phi\left(\frac27\right)$": Phi(R(2, 7)), r"$\Phi\left(\frac2{\sqrt7}\right)$": Phi(2 / sqrt(7))})
choix("normales-ex3", 2, PXgY(0, 2, 1, 2),
      {r"$1-\Phi\left(\frac12\right)$": 1 - Phi(R(1, 2)), r"$\Phi\left(\frac12\right)$": Phi(R(1, 2)),
       r"$1-\Phi\left(\frac14\right)$": 1 - Phi(R(1, 4)), r"$1-\Phi\left(\frac1{2\sqrt2}\right)$": 1 - Phi(1 / (2 * sqrt(2)))})
assert [round(float(w), 2) for w in [PXgY(50, 1, 51, 4), PXgY(5, 9, 3, 16), PXgY(0, 2, 1, 2)]] == [0.33, 0.66, 0.31]

# --- normales-ex4 : X^2 ~ gamma(1/2, 1/2), E(X^4) = 3, X^2 + Y^2 ~ E(1/2).
yp = Symbol("yp", positive=True)
dX2 = diff(2 * Phi(sqrt(yp)) - 1, yp)
assert simplify(dX2 - exp(-yp / 2) / sqrt(2 * pi * yp)) == 0
assert simplify(R(1, 2)**R(1, 2) / gamma(R(1, 2)) * yp**(-R(1, 2)) * exp(-yp / 2) - dX2) == 0
assert I(x**4 * phi, (x, -oo, oo)) == 3 == R(1, 2) * R(3, 2) / R(1, 4)
# Leurres : e^{-y/2}/√(2π) (intégrale √2·... ≠ 1), 2 e^{-y/2}/√(2πy) (intégrale 2), e^{-y^2/2}/√(2π) (intégrale 1/2).
assert I(exp(-yp / 2) / sqrt(2 * pi), (yp, 0, oo)) != 1 and I(2 * exp(-yp / 2) / sqrt(2 * pi * yp), (yp, 0, oo)) == 2
assert I(exp(-yp**2 / 2) / sqrt(2 * pi), (yp, 0, oo)) == R(1, 2)
# Les propositions sont des fonctions : on les compare par leur valeur en y = 4.
choix("normales-ex4", 0, dX2.subs(yp, 4), {r"$y\mapsto\frac{e^{-y/2}}{\sqrt{2\pi y}}$": exp(-2) / sqrt(8 * pi),
                                           r"$y\mapsto\frac{e^{-y/2}}{\sqrt{2\pi}}$": exp(-2) / sqrt(2 * pi),
                                           r"$y\mapsto\frac{2e^{-y/2}}{\sqrt{2\pi y}}$": 2 * exp(-2) / sqrt(8 * pi),
                                           r"$y\mapsto\frac{e^{-y^2/2}}{\sqrt{2\pi}}$": exp(-8) / sqrt(2 * pi)})

# --- normales-ex5 et variantes : P(X^2 + Y^2 <= r^2) = 1 - e^{-r^2/2}.
rp = Symbol("rp", positive=True)
assert simplify(I(I(exp(-rho**2 / 2) * rho / (2 * pi), (rho, 0, rp)), (th2, 0, 2 * pi)) - (1 - exp(-rp**2 / 2))) == 0
choix("normales-ex5", 0, exp(-R(4, 2)), {r"$e^{-2}$": exp(-2), r"$1-e^{-2}$": 1 - exp(-2), r"$e^{-4}$": exp(-4), r"$e^{-1}$": exp(-1)})
choix("normales-ex5", 1, exp(-R(2, 2)), {r"$e^{-1}$": exp(-1), r"$1-e^{-1}$": 1 - exp(-1), r"$e^{-2}$": exp(-2),
                                         r"$e^{-1/2}$": exp(-R(1, 2))})
choix("normales-ex5", 2, 1 - exp(-R(6, 2)), {r"$1-e^{-3}$": 1 - exp(-3), r"$e^{-3}$": exp(-3), r"$1-e^{-6}$": 1 - exp(-6),
                                             r"$1-e^{-3/2}$": 1 - exp(-R(3, 2))})
assert [round(math.exp(-2), 3), round(math.exp(-1), 3), round(1 - math.exp(-3), 3)] == [0.135, 0.368, 0.950]

# ===========================================================================================
# QCM
# ===========================================================================================
# q1 : seule ½e^{-|x|} est une densité ; ∫ 1/(1+x^2) = π ; x 1_[-1,√3] d'intégrale 1 mais négative.
assert I(exp(-abs(x)) / 2, (x, -oo, oo)) == 1 and I(1 / (1 + x**2), (x, -oo, oo)) == pi
assert I(x, (x, -1, sqrt(3))) == 1 and I(exp(-x), (x, -oo, 0)) == oo
assert QCM[f"{CH}-q1"]["propositions"][QCM[f"{CH}-q1"]["bonne"]] == r"$x\mapsto\frac12e^{-|x|}$"
assert QCM[f"{CH}-q2"]["vrai"] is False and QCM[f"{CH}-q4"]["vrai"] is False and QCM[f"{CH}-q6"]["vrai"] is True
choix_qcm("q3", log(2), {r"$\ln2$": log(2), r"$\frac12$": R(1, 2), r"$2$": 2, r"$\frac1{\ln2}$": 1 / log(2)})
assert exp(-log(2)) == R(1, 2)
choix_qcm("q5", I(exp(x), (x, 0, 1)), {r"$e-1$": E - 1, r"$\sqrt e$": sqrt(E), r"$e$": E, r"$\frac e2$": E / 2})
# q7 : comparées en z = 3.
h3 = I(exp(-x) * exp(-(zp - x)), (x, 0, zp)).subs(zp, 3)
choix_qcm("q7", h3, {r"$ze^{-z}$": 3 * exp(-3), r"$e^{-z}$": exp(-3), r"$2e^{-2z}$": 2 * exp(-6), r"$\frac{z^2}2e^{-z}$": R(9, 2) * exp(-3)})
assert simplify(1 - exp(-2 * zp) - I(2 * exp(-2 * z), (z, 0, zp))) == 0  # densité de min(X,Y) : P(min > z) = e^{-2z}
# q8 : X + 2Y ~ N(1 - 2, 4 + 4·5) ; propositions codées par (moyenne, variance).
q8 = QCM[f"{CH}-q8"]
code = {r"$\mathcal{N}(-1,24)$": (-1, 24), r"$\mathcal{N}(-1,14)$": (-1, 14), r"$\mathcal{N}(-1,9)$": (-1, 9),
        r"$\mathcal{N}(1,24)$": (1, 24)}
assert set(q8["propositions"]) == set(code) and code[q8["propositions"][q8["bonne"]]] == (1 - 2, 4 + 4 * 5)

print("proba_l3_densites : toutes les vérifications passent.")
