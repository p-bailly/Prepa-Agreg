"""Vérifications calculatoires du chapitre « int-l3-multiples » (intégration, L3) :
cartes, exercices et QCM de content/*/integration/l3-int-l3-multiples.json et content/qcm/int-l3-multiples.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Les identités sont vérifiées en SymPy ; les intégrales généralisées délicates numériquement (mpmath, 20 chiffres).
Pour chaque réponse à choix, la bonne valeur est recalculée et chaque leurre est vérifié différent.
Lancer : python3 verification/int_l3_multiples.py (doit finir sans erreur).
"""

import mpmath as mp
from sympy import (E, Abs, Matrix, Rational, atan, cos, cosh, diff, exp, gamma, integrate, limit, log, oo, pi, simplify,
                   sin, sinh, sqrt, symbols)

mp.mp.dps = 20
x, y, z, t, s, u, v, r, th, ph, rho = symbols("x y z t s u v r theta phi rho", real=True)
xp, yp, tp, Rp, ap, dp = symbols("xp yp tp R a d", positive=True)
R_ = Rational


def egal(a, b):
    return simplify(a - b) == 0


def choix(vraie, bonne, leurres, nom):
    """La valeur calculée est la bonne proposition, et chaque leurre en diffère."""
    assert egal(vraie, bonne), (nom, vraie, bonne)
    for l in leurres:
        assert abs(float((vraie - l).evalf(30))) > 1e-9, (nom, l)


def norme(vec):
    return sqrt(sum(c**2 for c in vec))


def longueur(gam, a, b):
    d = [diff(c, t) for c in gam]
    return integrate(simplify(norme(d)), (t, a, b))


def travail(V, gam, a, b):
    vals = {k: c for k, c in zip((x, y, z), gam)}
    d = [diff(c, t) for c in gam]
    return integrate(simplify(sum(Vi.subs(vals) * di for Vi, di in zip(V, d))), (t, a, b))


# ===========================================================================
# int-l3-multiples-curvilignes (carte)
# ===========================================================================
# Exemple : arche de cycloïde, longueur 8 (|γ'|² = 4 sin²(t/2), sin(t/2) ≥ 0 sur [0, 2π]).
assert egal((1 - cos(t)) ** 2 + sin(t) ** 2 - 4 * sin(t / 2) ** 2, 0)
assert integrate(2 * sin(t / 2), (t, 0, 2 * pi)) == 8
# Exemple : V = (y, x) = ∇(xy), le long de (t, 2t²) : 6t², intégrale 2 = 1·2 − 0.
assert travail((y, x), (t, 2 * t**2), 0, 1) == 2
# Contre-exemple : V = (−y, x)/(x²+y²) vérifie ∂P/∂y = ∂Q/∂x mais a une circulation 2π sur le cercle unité.
P, Q = -y / (x**2 + y**2), x / (x**2 + y**2)
assert egal(diff(P, y), diff(Q, x)) and egal(diff(P, y), (y**2 - x**2) / (x**2 + y**2) ** 2)
assert travail((P, Q), (cos(t), sin(t)), 0, 2 * pi) == 2 * pi
# Contre-exemple : φ(s) = (1 − cos 3πs)/2 surjective sur [0,1], ∫|φ'| = 3 (le segment parcouru trois fois).
phi_ = (1 - cos(3 * pi * s)) / 2
assert phi_.subs(s, 0) == 0 and phi_.subs(s, R_(1, 3)) == 1 and phi_.subs(s, R_(2, 3)) == 0 and phi_.subs(s, 1) == 1
assert abs(mp.quad(lambda w: abs(3 * mp.pi / 2 * mp.sin(3 * mp.pi * w)), [0, mp.mpf(1) / 3, mp.mpf(2) / 3, 1]) - 3) < 1e-15
# Contre-exemple : V = (y, x) sur (t, t) : travail 1, et −1 en sens inverse (φ(s) = 1 − s).
assert travail((y, x), (t, t), 0, 1) == 1 and travail((y, x), (1 - t, 1 - t), 0, 1) == -1
# Trou 1 : leurres faux. γ(t) = (t, 0) : ∫‖γ‖ = 1/2 ≠ 1 ; demi-cercle : ‖∫γ'‖ = 2 ≠ π ; cercle de rayon 2 : ∫‖γ'‖² = 8π ≠ 4π.
assert integrate(t, (t, 0, 1)) == R_(1, 2) and longueur((t, 0 * t), 0, 1) == 1
assert longueur((cos(t), sin(t)), 0, pi) == pi
assert integrate(4, (t, 0, 2 * pi)) == 8 * pi and longueur((2 * cos(t), 2 * sin(t)), 0, 2 * pi) == 4 * pi
# Trou 2 : un φ décroissant (bijectif, d'image [a, b]) change le signe du travail : voir ci-dessus (1 et −1).
# Trou 3 : f(x, y) = x, segment de (0,0) à (1,0) : f(B) − f(A) = 1 ; V = −∇f donne −1 ; V = (0, 1) (même norme
# que ∇f, et orthogonal à ∇f) donne 0.
assert travail((-1 + 0 * x, 0 * x), (t, 0 * t), 0, 1) == -1 and travail((0 * x, 1 + 0 * x), (t, 0 * t), 0, 1) == 0
# Affirmation 1 : le cercle parcouru deux fois a pour longueur 4π.
assert longueur((cos(t), sin(t)), 0, 4 * pi) == 4 * pi

# --- Exercices -------------------------------------------------------------
# int-l3-multiples-curvilignes-ex1 : (t², t³) sur [0,1] : (13√13 − 8)/27 ; leurres 13√13/27, (13√13 − 8)/18, √2.
L1 = integrate(tp * sqrt(4 + 9 * tp**2), (tp, 0, 1))
choix(L1, (13 * sqrt(13) - 8) / 27, [13 * sqrt(13) / 27, (13 * sqrt(13) - 8) / 18, sqrt(2)], "curvilignes-ex1")
assert egal(sqrt(4 * tp**2 + 9 * tp**4) - tp * sqrt(4 + 9 * tp**2), 0)
assert abs(float(L1) - 1.4397) < 1e-3 and float(L1) > float(sqrt(2))
# curvilignes-ex1 v1 : astroïde (cos³t, sin³t) sur [0, π/2] : 3/2 ; leurres 3/4, √2, 3.
assert egal((diff(cos(t) ** 3, t)) ** 2 + (diff(sin(t) ** 3, t)) ** 2 - 9 * cos(t) ** 2 * sin(t) ** 2, 0)
choix(integrate(3 * cos(t) * sin(t), (t, 0, pi / 2)), R_(3, 2), [R_(3, 4), sqrt(2), 3], "curvilignes-ex1 v1")
# curvilignes-ex1 v2 : hélice (cos t, sin t, t) sur [0, 2π] : 2π√2 ; leurres 2π, 4π, √(4π²+1).
choix(longueur((cos(t), sin(t), t), 0, 2 * pi), 2 * pi * sqrt(2), [2 * pi, 4 * pi, sqrt(4 * pi**2 + 1)], "curvilignes-ex1 v2")

# int-l3-multiples-curvilignes-ex2 : V = (−y, x) ; parabole 1/3, segment 0 ; leurres 0, 2/3, −1/3.
Vr = (-y, x)
choix(travail(Vr, (t, t**2), 0, 1), R_(1, 3), [0, R_(2, 3), R_(-1, 3)], "curvilignes-ex2")
assert travail(Vr, (t, t), 0, 1) == 0 and diff(-y, y) == -1 and diff(x, x) == 1
# curvilignes-ex2 v1 : cubique 1/2 ; leurres 0, 1/4, 3/4.
choix(travail(Vr, (t, t**3), 0, 1), R_(1, 2), [0, R_(1, 4), R_(3, 4)], "curvilignes-ex2 v1")
# curvilignes-ex2 v2 : quart de cercle π/2, segment (1−t, t) : 1 ; leurres 0, π, 1.
choix(travail(Vr, (cos(t), sin(t)), 0, pi / 2), pi / 2, [0, pi, 1], "curvilignes-ex2 v2")
assert travail(Vr, (1 - t, t), 0, 1) == 1

# int-l3-multiples-curvilignes-ex3 : V = (2xy + z, x², x) = ∇(x²y + xz) ; de (0,0,0) à (1,2,3) : 5 ; leurres 2, 3, 7.
f3 = x**2 * y + x * z
assert [diff(f3, w) for w in (x, y, z)] == [2 * x * y + z, x**2, x]
gam3 = (t, 2 * t**3, 3 * sin(pi * t / 2))
assert [c.subs(t, 0) for c in gam3] == [0, 0, 0] and [c.subs(t, 1) for c in gam3] == [1, 2, 3]
choix(f3.subs({x: 1, y: 2, z: 3}) - f3.subs({x: 0, y: 0, z: 0}), 5, [2, 3, 7], "curvilignes-ex3")
# Contrôle direct par l'intégrale curviligne (numérique).
_g = [c for c in gam3]
_d = [diff(c, t) for c in gam3]
_V = (2 * x * y + z, x**2, x)
_int = sum(Vi.subs({x: _g[0], y: _g[1], z: _g[2]}) * di for Vi, di in zip(_V, _d))
assert abs(mp.quad(lambda w: mp.mpf(str(_int.subs(t, w).evalf(25))), [0, 1]) - 5) < 1e-12
# curvilignes-ex3 v1 : V = (yz, xz, xy) = ∇(xyz), de (1,1,1) à (2,3,4) : 23 ; leurres 24, 9, 22.
assert [diff(x * y * z, w) for w in (x, y, z)] == [y * z, x * z, x * y]
choix(R_(2 * 3 * 4 - 1), 23, [24, 9, 22], "curvilignes-ex3 v1")
assert [c.subs(t, 1) for c in (1 + t, 1 + 2 * t**2, 1 + 3 * t)] == [2, 3, 4]
# curvilignes-ex3 v2 : V = (e^y, xe^y + 2y) = ∇(xe^y + y²), de (0,0) à (2, ln 3) : 6 + ln²3 ; leurres 6, 2 + ln²3, 6 + 2 ln 3.
f3b = x * exp(y) + y**2
assert egal(diff(f3b, x) - exp(y), 0) and egal(diff(f3b, y) - (x * exp(y) + 2 * y), 0)
choix(f3b.subs({x: 2, y: log(3)}) - f3b.subs({x: 0, y: 0}), 6 + log(3) ** 2, [6, 2 + log(3) ** 2, 6 + 2 * log(3)], "curvilignes-ex3 v2")

# int-l3-multiples-curvilignes-ex4 : V = (−y, x)/(x²+y²) ; circulation 2π sur tout cercle de centre 0 ; ∇ arctan(y/x) = V.
assert travail((P, Q), (Rp * cos(t), Rp * sin(t)), 0, 2 * pi) == 2 * pi
assert egal(diff(atan(y / x), x) - P, 0) and egal(diff(atan(y / x), y) - Q, 0)
# Question 4 : cercle de centre (3,0) et de rayon 1 (contenu dans x ≥ 2) : 0 ; leurres 2π, π, −2π.
choix(travail((P, Q), (3 + cos(t), sin(t)), 0, 2 * pi), 0, [2 * pi, pi, -2 * pi], "curvilignes-ex4")
# curvilignes-ex4 v1 : cercle de rayon 2, deux tours en sens indirect : −4π ; leurres 4π, 0, −8π.
assert egal(sum(Vi.subs({x: 2 * cos(t), y: -2 * sin(t)}) * di for Vi, di in zip((P, Q), (diff(2 * cos(t), t), diff(-2 * sin(t), t)))), -1)
choix(R_(-1) * 4 * pi, -4 * pi, [4 * pi, 0, -8 * pi], "curvilignes-ex4 v1")
# curvilignes-ex4 v2 : arc (2cos t, 2 sin t), t ∈ [−π/3, π/3] : 2π/3 = arctan √3 − arctan(−√3) ; leurres 4π/3, π/3, 0.
choix(travail((P, Q), (2 * cos(t), 2 * sin(t)), -pi / 3, pi / 3), 2 * pi / 3, [4 * pi / 3, pi / 3, 0], "curvilignes-ex4 v2")
assert egal(atan(sqrt(3)) - atan(-sqrt(3)), 2 * pi / 3)

# QCM q1 : longueur de y = ch x sur [0,1] : sh 1 ; leurres ch 1, ch 1 − 1, √2 ; corde ≈ 1,14.
choix(integrate(cosh(x), (x, 0, 1)), sinh(1), [cosh(1), cosh(1) - 1, sqrt(2)], "q1")
assert egal(1 + sinh(x) ** 2 - cosh(x) ** 2, 0)
assert abs(float(sqrt(1 + (cosh(1) - 1) ** 2)) - 1.138) < 1e-3
# QCM q8 : (2xy, x²) = ∇(x²y).
assert [diff(x**2 * y, w) for w in (x, y)] == [2 * x * y, x**2]

# ===========================================================================
# int-l3-multiples-fubini-rectangle (carte)
# ===========================================================================
# Exemple : ∬ xy² sur [0,1]×[0,2] = 4/3 ; ∬ x e^{xy} sur [0,1]² = e − 2 (dans les deux ordres).
assert integrate(x * y**2, (x, 0, 1), (y, 0, 2)) == R_(4, 3)
assert egal(integrate(x * exp(x * y), (y, 0, 1)), exp(x) - 1)
assert egal(integrate(exp(x) - 1, (x, 0, 1)), E - 2)
assert abs(mp.quad(lambda a, b: a * mp.e ** (a * b), [0, 1], [0, 1]) - (mp.e - 2)) < 1e-15
# Cours : ∫_0^1 (x⁴ − x)/ln x dx = ln(5/2).
assert abs(mp.quad(lambda w: (w**4 - w) / mp.log(w), [0, 1]) - mp.log(mp.mpf(5) / 2)) < 1e-15
# Contre-exemple 1 (et q2, ex4) : f = (x² − y²)/(x² + y²)², ∫_0^1 f dy = 1/(1+x²), A = π/4, B = −π/4.
f4 = (x**2 - y**2) / (x**2 + y**2) ** 2
assert egal(diff(y / (x**2 + y**2), y), f4)
assert egal(integrate(f4.subs(x, xp), (y, 0, 1)), 1 / (1 + xp**2))
assert integrate(1 / (1 + x**2), (x, 0, 1)) == pi / 4
assert egal(f4.subs({x: y, y: x}, simultaneous=True), -f4)
# Contre-exemple 2 (et q3) : (2 − xy)xy e^{−xy} = ∂_y(xy² e^{−xy}) = ∂_x(x²y e^{−xy}) ; A = 1, B = 0.
g2 = (2 - x * y) * x * y * exp(-x * y)
assert egal(diff(x * y**2 * exp(-x * y), y), g2) and egal(diff(x**2 * y * exp(-x * y), x), g2)
assert egal(integrate(g2, (y, 0, 1)), x * exp(-x)) and integrate(x * exp(-x), (x, 0, oo)) == 1
assert limit(xp**2 * yp * exp(-xp * yp), xp, oo) == 0 and g2.subs(y, 0) == 0
# Trou 1, leurre « séparément continue » : f = xy/(x²+y²)² (f(0,0) = 0), ∫_0^1 f dy = 1/(2x) − x/(2(1+x²)), non intégrable en 0.
assert egal(integrate(xp * y / (xp**2 + y**2) ** 2, (y, 0, 1)), 1 / (2 * xp) - xp / (2 * (1 + xp**2)))
# Relecture : même leurre avec des intégrales itérées finies et différentes. f = xy(x² − y²)/(x² + y²)³, f(0,0) = 0 :
# f(x, 0) = f(0, y) = 0, donc f est continue en chaque variable séparément sur [0,1]², mais les itérées valent 1/8 et −1/8.
fsep = x * y * (x**2 - y**2) / (x**2 + y**2) ** 3
assert fsep.subs(y, 0) == 0 and fsep.subs(x, 0) == 0
assert egal(integrate(fsep.subs(x, xp), (y, 0, 1)), xp / (2 * (1 + xp**2) ** 2))
assert integrate(x / (2 * (1 + x**2) ** 2), (x, 0, 1)) == R_(1, 8)
assert egal(fsep.subs({x: y, y: x}, simultaneous=True), -fsep)  # l'autre ordre donne −1/8
# Trou 2, leurres : f = x sur [0,1]×[2,3] (∬ = 1/2, leurre 1 = 5/2) ; f = y sur [0,1]² (∬ = 1/2, leurres 2 et 3 = 0).
assert integrate(x, (x, 0, 1), (y, 2, 3)) == R_(1, 2) and integrate(integrate(x, (x, 2, 3)), (y, 0, 1)) == R_(5, 2)
assert integrate(y, (x, 0, 1), (y, 0, 1)) == R_(1, 2)
# Affirmations : ∬ xy sur [0,1]×[0,2] = 1 ; ∬ (x − y) sur [0,1]² = 0.
assert integrate(x * y, (x, 0, 1), (y, 0, 2)) == 1 and integrate(x - y, (x, 0, 1), (y, 0, 1)) == 0

# --- Exercices -------------------------------------------------------------
# int-l3-multiples-fubini-rectangle-ex1 : ∬ (x + y²) sur [0,1]×[0,2] = 11/3 ; leurres 19/6, 7/3, 14/3.
choix(integrate(x + y**2, (y, 0, 2), (x, 0, 1)), R_(11, 3), [R_(19, 6), R_(7, 3), R_(14, 3)], "fubini-rectangle-ex1")
assert integrate(x + y**2, (x, 0, 1), (y, 0, 2)) == R_(11, 3)
# v1 : ∬ xy² sur [0,2]×[0,1] = 2/3 ; leurres 4/3, 1/6, 1/3.
choix(integrate(x * y**2, (x, 0, 2), (y, 0, 1)), R_(2, 3), [R_(4, 3), R_(1, 6), R_(1, 3)], "fubini-rectangle-ex1 v1")
# v2 : ∬ sin(x + y) sur [0, π/2]² = 2 ; leurres 1, 0, π/2.
choix(integrate(sin(x + y), (y, 0, pi / 2), (x, 0, pi / 2)), 2, [1, 0, pi / 2], "fubini-rectangle-ex1 v2")

# int-l3-multiples-fubini-rectangle-ex2 : ∬ x/(1+xy) sur [0,2]×[0,1] = 3 ln 3 − 2 ; leurres 3 ln 3, 2 ln 3 − 1, ln 3.
assert egal(integrate(xp / (1 + xp * y), (y, 0, 1)), log(1 + xp))
choix(integrate(log(1 + x), (x, 0, 2)), 3 * log(3) - 2, [3 * log(3), 2 * log(3) - 1, log(3)], "fubini-rectangle-ex2")
assert egal(1 / y - 1 / (y * (1 + x * y)) - x / (1 + x * y), 0)
# v1 : ∬ y sin(xy) sur [0,π]×[0,1] = 1 ; leurres 2, 2/π, 1 − 2/π.
assert egal(integrate(yp * sin(x * yp), (x, 0, pi)), 1 - cos(pi * yp))
choix(integrate(1 - cos(pi * y), (y, 0, 1)), 1, [2, 2 / pi, 1 - 2 / pi], "fubini-rectangle-ex2 v1")
# v2 : ∬ y/(1+xy) sur [0,1]² = 2 ln 2 − 1 ; leurres ln 2, 2 ln 2, 1 − ln 2.
assert egal(integrate(yp / (1 + x * yp), (x, 0, 1)), log(1 + yp))
choix(integrate(log(1 + y), (y, 0, 1)), 2 * log(2) - 1, [log(2), 2 * log(2), 1 - log(2)], "fubini-rectangle-ex2 v2")


# int-l3-multiples-fubini-rectangle-ex3 (et v1, v2) : ∫_0^1 (x^n − x^a)/ln x = ∫_a^b dy/(y+1) = ln((b+1)/(a+1)).
def frullani_log(n, a, b, attendu, leurres, nom):
    assert n == b
    Ln = symbols("Ln", negative=True)  # Ln = ln x pour 0 < x < 1
    assert egal(integrate(exp(y * Ln), (y, a, b)) - (exp(b * Ln) - exp(a * Ln)) / Ln, 0)
    assert limit((x**n - x**a) / log(x), x, 1) == b - a and limit((x**n - x**a) / log(x), x, 0, "+") == 0
    val = integrate(1 / (y + 1), (y, a, b))
    choix(val, attendu, leurres, nom)
    assert abs(mp.quad(lambda w: (w**n - w**a) / mp.log(w), [0, 1]) - mp.mpf(str(val.evalf(25)))) < 1e-15


frullani_log(2, 1, 2, log(R_(3, 2)), [log(2), R_(1, 2), log(3)], "fubini-rectangle-ex3")
frullani_log(3, 1, 3, log(2), [log(R_(3, 2)), log(3), R_(2, 3)], "fubini-rectangle-ex3 v1")
frullani_log(4, 2, 4, log(R_(5, 3)), [log(2), log(R_(5, 4)), R_(2, 5)], "fubini-rectangle-ex3 v2")

# int-l3-multiples-fubini-rectangle-ex4 : A = π/4 (leurres −π/4, 0, π/2) ; ∫|f| dy = 1/x − 1/(1+x²).
choix(pi / 4, pi / 4, [-pi / 4, 0, pi / 2], "fubini-rectangle-ex4")
assert egal(integrate(f4.subs(x, xp), (y, 0, xp)) - integrate(f4.subs(x, xp), (y, xp, 1)), 1 / xp - 1 / (1 + xp**2))
# v1 : f = (x − y)/(x + y)³ = ∂_y(y/(x+y)²) ; ∫_0^1 f dy = 1/(1+x)² ; A = 1/2 (leurres −1/2, 0, 1) ; ∫|f| dy = 1/(2x) − 1/(1+x)².
f4b = (x - y) / (x + y) ** 3
assert egal(diff(y / (x + y) ** 2, y), f4b) and egal(f4b.subs({x: y, y: x}, simultaneous=True), -f4b)
assert egal(integrate(f4b.subs(x, xp), (y, 0, 1)), 1 / (1 + xp) ** 2)
choix(integrate(1 / (1 + x) ** 2, (x, 0, 1)), R_(1, 2), [R_(-1, 2), 0, 1], "fubini-rectangle-ex4 v1")
assert egal(integrate(f4b.subs(x, xp), (y, 0, xp)) - integrate(f4b.subs(x, xp), (y, xp, 1)), 1 / (2 * xp) - 1 / (1 + xp) ** 2)
assert egal(f4b.subs(y, 0), 1 / x**2)
# v2 : f = (y² − x²)/(x² + y²)² = ∂_y(−y/(x²+y²)) ; A = −π/4 (leurres π/4, 0, −π/2).
assert egal(diff(-y / (x**2 + y**2), y), -f4)
choix(-integrate(1 / (1 + x**2), (x, 0, 1)), -pi / 4, [pi / 4, 0, -pi / 2], "fubini-rectangle-ex4 v2")

# int-l3-multiples-fubini-rectangle-ex5 : ∬_R ∂x∂y f = f(b,d) − f(a,d) − f(b,c) + f(a,c), sur un exemple polynomial.
F5 = x**3 * y**2 + sin(x * y)
a5, b5, c5, d5 = 1, 2, 1, 3
I5 = integrate(diff(diff(F5, y), x), (x, a5, b5), (y, c5, d5))
coins = F5.subs({x: b5, y: d5}) - F5.subs({x: a5, y: d5}) - F5.subs({x: b5, y: c5}) + F5.subs({x: a5, y: c5})
assert egal(I5, coins)
assert not egal(coins, F5.subs({x: b5, y: d5}) - F5.subs({x: a5, y: c5}))
assert not egal(coins, F5.subs({x: b5, y: d5}) + F5.subs({x: a5, y: d5}) - F5.subs({x: b5, y: c5}) - F5.subs({x: a5, y: c5}))
assert not egal(coins, diff(F5, y).subs({x: b5, y: d5}) - diff(F5, y).subs({x: a5, y: c5}))

# ===========================================================================
# int-l3-multiples-fubini-tonelli (carte)
# ===========================================================================
# Exemple (Frullani) : ∫_0^∞ (e^{−at} − e^{−bt})/t dt = ln(b/a), ici a = 1, b = 3 ; e^{−at} − e^{−bt} = t ∫_a^b e^{−xt} dx.
assert egal(tp * integrate(exp(-x * tp), (x, ap, 3 * ap)) - (exp(-ap * tp) - exp(-3 * ap * tp)), 0)
assert abs(mp.quad(lambda w: (mp.e**-w - mp.e ** (-3 * w)) / w, [0, 1, mp.inf]) - mp.log(3)) < 1e-15
# Contre-exemple : f = (x − y)/(x + y)³ sur [1, +∞[² : ∫_1^∞ f dy = −1/(x+1)², itérées −1/2 et 1/2 ; |f| ≤ 1/4.
assert egal(integrate(f4b.subs(x, xp + 1), (y, 1, oo)), -1 / (xp + 2) ** 2)
assert integrate(-1 / (x + 1) ** 2, (x, 1, oo)) == R_(-1, 2)
# ∫_1^∞ |f(x,y)| dy ~ c/x avec c = ∫_0^∞ |1 − w|/(1 + w)³ dw > 0.
c_ = mp.quad(lambda w: abs(1 - w) / (1 + w) ** 3, [0, 1, mp.inf])
X = mp.mpf(10) ** 4
I_abs = mp.quad(lambda w: abs(X - w) / (X + w) ** 3, [1, X, mp.inf])
assert c_ > 0.1 and abs(X * I_abs / c_ - 1) < 1e-3
# Relecture, second contre-exemple (2 − xy)xy e^{−xy} : pour y > 0, ∫_0^∞ |f(x,y)| dx = c/y avec c = ∫_0^∞ |2 − w| w e^{−w} dw > 0
# (w = xy), et 0 pour y = 0 : intégrale partielle de |f| non continue en 0, intégrale itérée de |f| infinie.
c2 = mp.quad(lambda w: abs(2 - w) * w * mp.e ** (-w), [0, 2, mp.inf])
for yy in (mp.mpf(1) / 3, mp.mpf(2), mp.mpf(1) / 50):
    J_abs = mp.quad(lambda w: abs((2 - w * yy) * w * yy * mp.e ** (-w * yy)), [0, 2 / yy, mp.inf])
    assert abs(J_abs - c2 / yy) < 1e-12
assert abs(c2 - 8 * mp.e ** -2) < 1e-15  # c = 8 e^{−2} (∫_0^∞ (2 − w) w e^{−w} dw = 0)
# Trou 1, leurres : ce contre-exemple est borné, à valeurs réelles, de classe C^∞.
# Trou 2, leurres : 1/(x+y)² sur [1,∞[² : ∫_1^∞ dy = 1/(x+1), intégrale itérée infinie (pas « finies ») ;
# 1/√x sur ]0,1]×[0,1] : non bornée, itérées égales à 2 ; e^{−x−y} sur [0,∞[² : égales à 1 sur des intervalles non bornés.
assert egal(integrate(1 / (xp + 1 + y) ** 2, (y, 1, oo)), 1 / (xp + 2)) and integrate(1 / (x + 1), (x, 1, oo)) == oo
assert integrate(1 / sqrt(x), (x, 0, 1)) == 2
# Affirmation 3 : ∫∫ e^{−x−y} = 1.
assert integrate(exp(-x - y), (y, 0, oo), (x, 0, oo)) == 1

# --- Exercices -------------------------------------------------------------
# int-l3-multiples-fubini-tonelli-ex1 : ∫_1^∞ ∫_1^∞ (x+y)^{-3} = 1/4 ; leurres 1/2, 1/8, +∞ (évalué comme 10^6).
assert egal(integrate((xp + 1 + y) ** -3, (y, 1, oo)), 1 / (2 * (xp + 2) ** 2))
choix(integrate(1 / (2 * (x + 1) ** 2), (x, 1, oo)), R_(1, 4), [R_(1, 2), R_(1, 8), 10**6], "fubini-tonelli-ex1")
# v1 : (1+x+y)^{-3} sur [0,∞[² : 1/2 ; leurres 1/4, 1, +∞.
assert egal(integrate((1 + xp + y) ** -3, (y, 0, oo)), 1 / (2 * (1 + xp) ** 2))
choix(integrate(1 / (2 * (1 + x) ** 2), (x, 0, oo)), R_(1, 2), [R_(1, 4), 1, 10**6], "fubini-tonelli-ex1 v1")
# v2 : (1+x+y)^{-4} sur [0,∞[² : 1/6 ; leurres 1/3, 1/12, 1/9.
assert egal(integrate((1 + xp + y) ** -4, (y, 0, oo)), 1 / (3 * (1 + xp) ** 3))
choix(integrate(1 / (3 * (1 + x) ** 3), (x, 0, oo)), R_(1, 6), [R_(1, 3), R_(1, 12), R_(1, 9)], "fubini-tonelli-ex1 v2")

# int-l3-multiples-fubini-tonelli-ex2 : ∫_0^∞ (arctan 2t − arctan t)/t = (π/2) ln 2 ; leurres π ln 2, ln 2, (π/4) ln 2.
assert egal(integrate(tp / (1 + x**2 * tp**2), (x, 1, 2)) - (atan(2 * tp) - atan(tp)), 0)
assert egal(integrate(1 / (1 + xp**2 * t**2), (t, 0, oo)), pi / (2 * xp))
choix(integrate(pi / (2 * x), (x, 1, 2)), pi / 2 * log(2), [pi * log(2), log(2), pi / 4 * log(2)], "fubini-tonelli-ex2")
assert abs(mp.quad(lambda w: (mp.atan(2 * w) - mp.atan(w)) / w, [0, 1, mp.inf]) - mp.pi / 2 * mp.log(2)) < 1e-12
# v1 : (arctan 3t − arctan t)/t : (π/2) ln 3 ; leurres π ln 3, ln 3, (π/2) ln 2 ; g(0) = 2.
choix(integrate(pi / (2 * x), (x, 1, 3)), pi / 2 * log(3), [pi * log(3), log(3), pi / 2 * log(2)], "fubini-tonelli-ex2 v1")
assert abs(mp.quad(lambda w: (mp.atan(3 * w) - mp.atan(w)) / w, [0, 1, mp.inf]) - mp.pi / 2 * mp.log(3)) < 1e-12
assert integrate(1, (x, 1, 3)) == 2
# v2 : (e^{−t²} − e^{−2t²})/t : (ln 2)/2 ; ∫_1^2 t² e^{−xt²} dx = e^{−t²} − e^{−2t²} ; ∫_0^∞ t e^{−xt²} dt = 1/(2x).
assert egal(integrate(tp**2 * exp(-x * tp**2), (x, 1, 2)) - (exp(-tp**2) - exp(-2 * tp**2)), 0)
assert egal(integrate(t * exp(-xp * t**2), (t, 0, oo)), 1 / (2 * xp))
choix(integrate(1 / (2 * x), (x, 1, 2)), log(2) / 2, [log(2), sqrt(pi) / 2 * log(2), log(2) / 4], "fubini-tonelli-ex2 v2")
assert abs(mp.quad(lambda w: (mp.e ** (-w * w) - mp.e ** (-2 * w * w)) / w, [0, 1, mp.inf]) - mp.log(2) / 2) < 1e-15

# int-l3-multiples-fubini-tonelli-ex3 : f = x e^{−x²(1+y²)} ; ∫ dx = 1/(2(1+y²)) ; ∫ dy = G e^{−x²} ; G = √π/2.
assert egal(integrate(x * exp(-x**2 * (1 + yp**2)), (x, 0, oo)), 1 / (2 * (1 + yp**2)))
assert egal(integrate(xp * exp(-xp**2 * (1 + y**2)), (y, 0, oo)), exp(-xp**2) * sqrt(pi) / 2)
assert integrate(1 / (2 * (1 + y**2)), (y, 0, oo)) == pi / 4
choix(sqrt(pi / 4), sqrt(pi) / 2, [sqrt(pi), pi / 4, sqrt(pi) / 4], "fubini-tonelli-ex3")

# int-l3-multiples-fubini-tonelli-ex4 : ∫_0^A e^{−xt} sin t dt = (1 − e^{−xA}(x sin A + cos A))/(1 + x²) ; valeur π/2.
hA = (1 - exp(-xp * Rp) * (xp * sin(Rp) + cos(Rp))) / (1 + xp**2)
assert egal(integrate(exp(-xp * t) * sin(t), (t, 0, Rp)) - hA, 0)
assert egal(integrate(exp(-x * tp), (x, 0, oo)), 1 / tp) and egal(integrate(exp(-x * tp) * Abs(sin(tp)), (x, 0, oo)), Abs(sin(tp)) / tp)
# Domination : (1 + x) e^{−x} ≤ 1 pour x ≥ 0 (maximum 1 en 0, dérivée −x e^{−x}).
assert egal(diff((1 + x) * exp(-x), x), -x * exp(-x))
choix(integrate(1 / (1 + x**2), (x, 0, oo)), pi / 2, [pi, pi / 4, 1], "fubini-tonelli-ex4")
mp.mp.dps = 20
assert abs(mp.quadosc(lambda w: mp.sin(w) / w if w else mp.mpf(1), [0, mp.inf], omega=1) - mp.pi / 2) < 1e-12

# QCM q3 : voir le contre-exemple 2 de fubini-rectangle (A = 1, B = 0).

# ===========================================================================
# int-l3-multiples-domaines (carte)
# ===========================================================================
# Exemple : ∬_T xy = 1/8 dans les deux ordres.
assert integrate(x * y, (y, 0, x), (x, 0, 1)) == R_(1, 8) and integrate(x * y, (x, y, 1), (y, 0, 1)) == R_(1, 8)
# Contre-exemples : f = x, 1/3 contre 1/6 ; φ1 = x, φ2 = 0 : ∫(φ2 − φ1) = −1/2.
assert integrate(x, (y, 0, x), (x, 0, 1)) == R_(1, 3) and integrate(x, (x, 0, y), (y, 0, 1)) == R_(1, 6)
assert integrate(-x, (x, 0, 1)) == R_(-1, 2)
# Trou 1, leurre « φ1(a) ≤ φ2(a), φ1(b) ≤ φ2(b) » : φ1 = 0, φ2 = (2x−1)² − 1/2 ≥ 0 en 0 et 1, < 0 en 1/2 ;
# aire vraie ∫ max(φ2, 0) ≠ ∫ φ2.
p2 = (2 * x - 1) ** 2 - R_(1, 2)
assert p2.subs(x, 0) > 0 and p2.subs(x, 1) > 0 and p2.subs(x, R_(1, 2)) < 0
assert abs(mp.quad(lambda w: max((2 * w - 1) ** 2 - 0.5, 0), [0, 0.5 - mp.sqrt(2) / 4, 0.5 + mp.sqrt(2) / 4, 1])
           - float(integrate(p2, (x, 0, 1)))) > 0.1
# Trou 2, leurres : triangle (φ1 = 0, φ2 = x), f = 1 : 1 ≠ 1/2 ; f = y sur [0,1]² : 0 ≠ 1/2 ; f = 1, φ2 = 2 sur [0,1] : 4 ≠ 2.
assert integrate(1, (x, 0, 1), (y, 0, 1)) == 1 and integrate(1, (y, 0, x), (x, 0, 1)) == R_(1, 2)
assert integrate(2 * 2, (x, 0, 1)) == 4
# Affirmations : f = x : 1/3 contre 1/6 (ci-dessus) ; ∫_0^1 ∫_x^1 f dy dx = ∫_0^1 ∫_0^y f dx dy (exemple f = x²y).
assert integrate(x**2 * y, (y, x, 1), (x, 0, 1)) == integrate(x**2 * y, (x, 0, y), (y, 0, 1))
# QCM q4 : ∫_0^1∫_0^x f dy dx = ∫_0^1∫_y^1 f dx dy (exemple f = x y² + e^x).
assert egal(integrate(x * y**2 + exp(x), (y, 0, x), (x, 0, 1)), integrate(x * y**2 + exp(x), (x, y, 1), (y, 0, 1)))

# --- Exercices -------------------------------------------------------------
# int-l3-multiples-domaines-ex1 : ∬_T (x + y) = 1/2 ; leurres 1, 1/3, 1/4.
choix(integrate(x + y, (y, 0, x), (x, 0, 1)), R_(1, 2), [1, R_(1, 3), R_(1, 4)], "domaines-ex1")
# v1 : ∬ x sur {x² ≤ y ≤ x} = 1/12 ; leurres 1/6 (aire), 1/20, 1/15 (∬ y).
choix(integrate(x, (y, x**2, x), (x, 0, 1)), R_(1, 12), [R_(1, 6), R_(1, 20), R_(1, 15)], "domaines-ex1 v1")
assert integrate(1, (y, x**2, x), (x, 0, 1)) == R_(1, 6) and integrate(y, (y, x**2, x), (x, 0, 1)) == R_(1, 15)
# v2 : ∬ y sur le triangle x, y ≥ 0, x + y ≤ 1 : 1/6 ; leurres 1/3, 1/2, 1/24.
choix(integrate(y, (y, 0, 1 - x), (x, 0, 1)), R_(1, 6), [R_(1, 3), R_(1, 2), R_(1, 24)], "domaines-ex1 v2")

# int-l3-multiples-domaines-ex2 : ∫_0^1 ∫_x^1 e^{y²} = (e − 1)/2 ; leurres e − 1, e/2, (e − 2)/2.
choix(integrate(y * exp(y**2), (y, 0, 1)), (E - 1) / 2, [E - 1, E / 2, (E - 2) / 2], "domaines-ex2")
assert abs(mp.quad(lambda a: mp.quad(lambda b: mp.e ** (b * b), [a, 1]), [0, 1]) - (mp.e - 1) / 2) < 1e-15
# v1 : ∫_0^1 ∫_y^1 sin(x²) dx dy = (1 − cos 1)/2 ; leurres 1 − cos 1, (sin 1)/2, (1 + cos 1)/2.
choix(integrate(x * sin(x**2), (x, 0, 1)), (1 - cos(1)) / 2, [1 - cos(1), sin(1) / 2, (1 + cos(1)) / 2], "domaines-ex2 v1")
assert abs(mp.quad(lambda b: mp.quad(lambda a: mp.sin(a * a), [b, 1]), [0, 1]) - (1 - mp.cos(1)) / 2) < 1e-15
# v2 : ∫_0^1 ∫_{√x}^1 e^{y³} dy dx = (e − 1)/3 ; leurres (e − 1)/2, e − 1, e/3.
choix(integrate(y**2 * exp(y**3), (y, 0, 1)), (E - 1) / 3, [(E - 1) / 2, E - 1, E / 3], "domaines-ex2 v2")
assert abs(mp.quad(lambda a: mp.quad(lambda b: mp.e ** (b**3), [mp.sqrt(a), 1]), [0, 1]) - (mp.e - 1) / 3) < 1e-12

# int-l3-multiples-domaines-ex3 : parabole y = x², droite y = x + 2 : aire 9/2, ∬ x = 9/4 ; leurres 9/2, 1/2, 27/4.
assert integrate(x + 2 - x**2, (x, -1, 2)) == R_(9, 2)
choix(integrate(x * (x + 2 - x**2), (x, -1, 2)), R_(9, 4), [R_(9, 2), R_(1, 2), R_(27, 4)], "domaines-ex3")
# v1 : droite y = 2x : aire 4/3, ∬ y = 32/15 ; leurres 4/3, 16/15, 64/15.
assert integrate(2 * x - x**2, (x, 0, 2)) == R_(4, 3)
choix(integrate(y, (y, x**2, 2 * x), (x, 0, 2)), R_(32, 15), [R_(4, 3), R_(16, 15), R_(64, 15)], "domaines-ex3 v1")
# v2 : droite y = 1 : aire 4/3, ∬ y = 4/5 ; leurres 4/3, 2/5, 3/5 (ordonnée du centre de gravité).
assert integrate(1 - x**2, (x, -1, 1)) == R_(4, 3)
choix(integrate(y, (y, x**2, 1), (x, -1, 1)), R_(4, 5), [R_(4, 3), R_(2, 5), R_(3, 5)], "domaines-ex3 v2")

# int-l3-multiples-domaines-ex4 : ∫_0^x ∫_0^t f(s) ds dt = ∫_0^x (x − s) f(s) ds (vérifié pour f = cos, e^s, 1/(1+s)).
for f_ in (cos(s), exp(s), 1 / (1 + s)):
    gauche = integrate(integrate(f_, (s, 0, t)), (t, 0, xp))
    droite = integrate((xp - s) * f_, (s, 0, xp))
    assert egal(gauche - droite, 0)
# Application : f = 1/(1+s) : 2 ln 2 − 1 ; leurres ln 2, 1 − ln 2, 2 ln 2.
choix(integrate((1 - s) / (1 + s), (s, 0, 1)), 2 * log(2) - 1, [log(2), 1 - log(2), 2 * log(2)], "domaines-ex4")
assert egal(integrate(log(1 + t), (t, 0, 1)), 2 * log(2) - 1)
# v1 : f = e^s : e − 2 ; leurres e − 1, 1, e − 3.
choix(integrate((1 - s) * exp(s), (s, 0, 1)), E - 2, [E - 1, 1, E - 3], "domaines-ex4 v1")
# v2 : f = 1/(1+s)² : 1 − ln 2 ; leurres ln 2, 2 ln 2 − 1, 1/2.
choix(integrate((1 - s) / (1 + s) ** 2, (s, 0, 1)), 1 - log(2), [log(2), 2 * log(2) - 1, R_(1, 2)], "domaines-ex4 v2")
assert egal(integrate(1 - 1 / (1 + t), (t, 0, 1)), 1 - log(2))

# ===========================================================================
# int-l3-multiples-changement-variables (carte)
# ===========================================================================
Jpol = Matrix([r * cos(th), r * sin(th)]).jacobian([r, th])
assert simplify(Jpol.det()) == r
# Exemple : ∬_{disque R} e^{−(x²+y²)} = π(1 − e^{−R²}).
assert egal(2 * pi * integrate(exp(-r**2) * r, (r, 0, Rp)), pi * (1 - exp(-Rp**2)))
# Contre-exemples : sans r, 2π ≠ π ; Δ = [1,2]×[0,4π] : 6π, couronne d'aire 3π ; Φ(u,v) = (v,u) : jacobien −1.
assert integrate(1, (r, 0, 1), (th, 0, 2 * pi)) == 2 * pi and integrate(r, (r, 0, 1), (th, 0, 2 * pi)) == pi
assert integrate(r, (r, 1, 2), (th, 0, 4 * pi)) == 6 * pi and pi * (4 - 1) == 3 * pi
assert Matrix([v, u]).jacobian([u, v]).det() == -1
# Trou 2, leurres : Φ(u,v) = (2u, v) (det 2) sur [0,1]² : aire 2, alors que 1/|det| et det² donnent 1/2 et 4.
assert Matrix([2 * u, v]).jacobian([u, v]).det() == 2
# Affirmations : (2u, 3v) a pour jacobien 6 ; disque de centre (1,0) : (r cos θ − 1)² + (r sin θ)² ≤ 1 ⟺ r² ≤ 2r cos θ.
assert Matrix([2 * u, 3 * v]).jacobian([u, v]).det() == 6
assert egal((r * cos(th) - 1) ** 2 + (r * sin(th)) ** 2 - 1, r**2 - 2 * r * cos(th))
# Cours : le disque de centre (1,0) de rayon 1, en polaires 0 ≤ r ≤ 2cos θ, θ ∈ [−π/2, π/2], a pour aire π.
assert integrate(integrate(r, (r, 0, 2 * cos(th))), (th, -pi / 2, pi / 2)) == pi

# --- Exercices -------------------------------------------------------------
# int-l3-multiples-changement-variables-ex1 : ∬_{r≤2} (x²+y²) = 8π ; leurres 4π, 16π/3 (sans r), 16π.
choix(2 * pi * integrate(r**3, (r, 0, 2)), 8 * pi, [4 * pi, 16 * pi / 3, 16 * pi], "changement-variables-ex1")
assert 2 * pi * integrate(r**2, (r, 0, 2)) == 16 * pi / 3
# v1 : ∬_{r≤1} √(x²+y²) = 2π/3 ; leurres π (sans r), π/3, 4π/3.
choix(2 * pi * integrate(r**2, (r, 0, 1)), 2 * pi / 3, [pi, pi / 3, 4 * pi / 3], "changement-variables-ex1 v1")
# v2 : ∬_{1≤r≤2} 1/(x²+y²) = 2π ln 2 ; leurres π ln 2, π (sans r), 4π ln 2.
choix(2 * pi * integrate(1 / r, (r, 1, 2)), 2 * pi * log(2), [pi * log(2), pi, 4 * pi * log(2)], "changement-variables-ex1 v2")
assert 2 * pi * integrate(1 / r**2, (r, 1, 2)) == pi

# int-l3-multiples-changement-variables-ex2 : ∬ xy sur le quart de disque = 1/8 ; leurres 1/4, 1/6, 1/16.
choix(integrate(r**3 * cos(th) * sin(th), (r, 0, 1), (th, 0, pi / 2)), R_(1, 8), [R_(1, 4), R_(1, 6), R_(1, 16)], "changement-variables-ex2")
assert integrate(x * y, (y, 0, sqrt(1 - x**2)), (x, 0, 1)) == R_(1, 8)
# v1 : ∬ x sur le demi-disque x ≥ 0 = 2/3 ; leurres 1/3, π/4, 1.
choix(integrate(r**2 * cos(th), (r, 0, 1), (th, -pi / 2, pi / 2)), R_(2, 3), [R_(1, 3), pi / 4, 1], "changement-variables-ex2 v1")
# v2 : ∬ x sur le quart de couronne = 7/3 ; leurres 3/2 (sans r), 8/3, 7π/6.
choix(integrate(r**2 * cos(th), (r, 1, 2), (th, 0, pi / 2)), R_(7, 3), [R_(3, 2), R_(8, 3), 7 * pi / 6], "changement-variables-ex2 v2")

# int-l3-multiples-changement-variables-ex3 : Φ(u,v) = ((u+v)/2, (u−v)/2), det −1/2.
Phi = Matrix([(u + v) / 2, (u - v) / 2])
assert Phi.jacobian([u, v]).det() == R_(-1, 2)
assert egal(Phi[0] + Phi[1] - u, 0) and egal(Phi[0] - Phi[1] - v, 0)
# ∬_D (x+y)² sur |x+y| ≤ 1, |x−y| ≤ 1 : 2/3 ; leurres 4/3, 1/3, 8/3.
choix(R_(1, 2) * integrate(u**2, (u, -1, 1), (v, -1, 1)), R_(2, 3), [R_(4, 3), R_(1, 3), R_(8, 3)], "changement-variables-ex3")
# Contrôle direct : D est le carré de sommets (±1, 0), (0, ±1) ; ∬ (x+y)² en cartésiennes.
assert abs(mp.quad(lambda a: mp.quad(lambda b: (a + b) ** 2, [-(1 - abs(a)), 1 - abs(a)]), [-1, 0, 1]) - mp.mpf(2) / 3) < 1e-15
# v1 : ∬ e^{x+y} sur 0 ≤ x+y ≤ 1, 0 ≤ x−y ≤ 2 : e − 1 ; leurres 2(e−1), (e−1)/2, 4(e−1).
choix(R_(1, 2) * integrate(exp(u), (u, 0, 1), (v, 0, 2)), E - 1, [2 * (E - 1), (E - 1) / 2, 4 * (E - 1)], "changement-variables-ex3 v1")
# v2 : ∬ (x+y)(x−y)² sur 0 ≤ x+y ≤ 2, 0 ≤ x−y ≤ 1 : 1/3 ; leurres 2/3, 1/6, 4/3.
choix(R_(1, 2) * integrate(u * v**2, (u, 0, 2), (v, 0, 1)), R_(1, 3), [R_(2, 3), R_(1, 6), R_(4, 3)], "changement-variables-ex3 v2")

# int-l3-multiples-changement-variables-ex4 : Γ(a) = 2∫ u^{2a−1} e^{−u²} du ; Γ(a)Γ(b) = Γ(a+b) B(a,b) ; application 1/24.
for a_ in (R_(1, 2), 1, R_(3, 2), 2, 3):
    assert egal(2 * integrate(u ** (2 * a_ - 1) * exp(-u**2), (u, 0, oo)), gamma(a_))
for a_, b_ in ((2, 3), (R_(1, 2), R_(1, 2)), (1, R_(3, 2))):
    B_ = 2 * integrate(cos(th) ** (2 * a_ - 1) * sin(th) ** (2 * b_ - 1), (th, 0, pi / 2))
    assert egal(gamma(a_) * gamma(b_), gamma(a_ + b_) * B_)
choix(integrate(cos(th) ** 3 * sin(th) ** 5, (th, 0, pi / 2)), R_(1, 24), [R_(1, 12), R_(1, 48), R_(1, 60)], "changement-variables-ex4")
assert integrate((1 - s**2) * s**5, (s, 0, 1)) == R_(1, 24) and gamma(2) * gamma(3) / gamma(5) == R_(1, 12)
assert egal((4 * (r * cos(th)) ** 3 * (r * sin(th)) ** 5 * exp(-r**2)) * r,
            4 * r ** (2 * 2 + 2 * 3 - 1) * exp(-r**2) * cos(th) ** 3 * sin(th) ** 5)

# QCM q5 : ∬_{r≤1} e^{x²+y²} = π(e − 1) ; leurres 2π(e−1), πe, π(e−1)/2.
choix(2 * pi * integrate(exp(r**2) * r, (r, 0, 1)), pi * (E - 1), [2 * pi * (E - 1), pi * E, pi * (E - 1) / 2], "q5")

# ===========================================================================
# int-l3-multiples-aires-volumes (carte)
# ===========================================================================
# Exemples : cardioïde 3π/2 ; volume sous le paraboloïde 1 − r² : π/2 ; boule 4πR³/3 par révolution.
assert integrate((1 + cos(th)) ** 2, (th, 0, 2 * pi)) / 2 == 3 * pi / 2
assert 2 * pi * integrate((1 - r**2) * r, (r, 0, 1)) == pi / 2
assert egal(pi * integrate(Rp**2 - x**2, (x, -Rp, Rp)), R_(4, 3) * pi * Rp**3)
# Contre-exemples (et q6) : ρ = 2cos θ : π sur [−π/2, π/2], 2π sur [0, 2π] ; signe : ∬ x sur [−1,1]×[0,1] = 0, volume 1.
assert integrate(4 * cos(th) ** 2, (th, -pi / 2, pi / 2)) / 2 == pi
assert integrate(4 * cos(th) ** 2, (th, 0, 2 * pi)) / 2 == 2 * pi
assert integrate(x, (x, -1, 1), (y, 0, 1)) == 0 and integrate(Abs(x), (x, -1, 1), (y, 0, 1)) == 1
# Trous, leurres (ρ = 1 sur [0, 2π] : aire π ; g = 1 sur [0, 1] : volume π) : 2π, 2π, 2π² ; 2π, 1, π/3.
assert integrate(1, (th, 0, 2 * pi)) / 2 == pi
for l in (2 * pi, 2 * pi, 2 * pi**2):
    assert l != pi
for l in (2 * pi, 1, pi / 3):
    assert l != pi
# Pièges : ρ = 1 sur [0, 4π] donne 2π ; f = −1 donne −1 ; g = −1 donne π alors que le domaine est vide.
assert integrate(1, (th, 0, 4 * pi)) / 2 == 2 * pi

# --- Exercices -------------------------------------------------------------
# int-l3-multiples-aires-volumes-ex1 : aire entre x² et √x : 1/3 ; leurres 2/3, 1/6, 1.
choix(integrate(sqrt(x) - x**2, (x, 0, 1)), R_(1, 3), [R_(2, 3), R_(1, 6), 1], "aires-volumes-ex1")
# v1 : entre x³ et x : 1/4 ; leurres 1/2, 3/4, 1/12.
choix(integrate(x - x**3, (x, 0, 1)), R_(1, 4), [R_(1, 2), R_(3, 4), R_(1, 12)], "aires-volumes-ex1 v1")
# v2 : sous sin sur [0, π] : 2 ; leurres 1, π, π/2.
choix(integrate(sin(x), (x, 0, pi)), 2, [1, pi, pi / 2], "aires-volumes-ex1 v2")

# int-l3-multiples-aires-volumes-ex2 : pétale de ρ = cos 2θ : π/8 ; leurres π/4, 1 (= ∫ρ), π/16.
choix(integrate(cos(2 * th) ** 2, (th, -pi / 4, pi / 4)) / 2, pi / 8, [pi / 4, 1, pi / 16], "aires-volumes-ex2")
assert integrate(cos(2 * th), (th, -pi / 4, pi / 4)) == 1
# v1 : spirale ρ = θ sur [0, 2π] : 4π³/3 ; leurres 8π³/3, 2π² (= ∫ρ), 2π³/3.
choix(integrate(th**2, (th, 0, 2 * pi)) / 2, 4 * pi**3 / 3, [8 * pi**3 / 3, 2 * pi**2, 2 * pi**3 / 3], "aires-volumes-ex2 v1")
# v2 : ρ = 2 + cos θ : 9π/2 ; leurres 9π, 4π (= ∫ρ), 5π/2.
choix(integrate((2 + cos(th)) ** 2, (th, 0, 2 * pi)) / 2, 9 * pi / 2, [9 * pi, 4 * pi, 5 * pi / 2], "aires-volumes-ex2 v2")
assert integrate(2 + cos(th), (th, 0, 2 * pi)) == 4 * pi

# int-l3-multiples-aires-volumes-ex3 : entre z = r² et z = 2 − r² : π ; leurres 2π, π/2, 4π/3.
choix(2 * pi * integrate((2 - 2 * r**2) * r, (r, 0, 1)), pi, [2 * pi, pi / 2, 4 * pi / 3], "aires-volumes-ex3")
# v1 : cône √(x²+y²) ≤ z ≤ 1 : π/3 ; leurres π/2, 2π/3, π.
choix(2 * pi * integrate((1 - r) * r, (r, 0, 1)), pi / 3, [pi / 2, 2 * pi / 3, pi], "aires-volumes-ex3 v1")
# v2 : 0 ≤ z ≤ 4 − r² : 8π ; leurres 16π, 4π, 32π/3 (sans r).
choix(2 * pi * integrate((4 - r**2) * r, (r, 0, 2)), 8 * pi, [16 * pi, 4 * pi, 32 * pi / 3], "aires-volumes-ex3 v2")
assert 2 * pi * integrate(4 - r**2, (r, 0, 2)) == 32 * pi / 3

# int-l3-multiples-aires-volumes-ex4 : deux cylindres (Steinmetz) : 16/3 ; leurres 4π/3, 8/3, 4π.
choix(integrate(4 * (1 - x**2), (x, -1, 1)), R_(16, 3), [R_(4, 3) * pi, R_(8, 3), 4 * pi], "aires-volumes-ex4")
assert float(R_(4, 3) * pi) < 16 / 3 and egal((R_(4, 3) * pi) / R_(16, 3), pi / 4)
# Contrôle Monte-Carlo déterministe (grille) du volume de Steinmetz.
_n = 60
_pts = [(-1 + (2 * i + 1) / _n) for i in range(_n)]
_cnt = sum(1 for a in _pts for b in _pts for c in _pts if a * a + b * b <= 1 and a * a + c * c <= 1)
assert abs(_cnt * (2 / _n) ** 3 - 16 / 3) < 0.05
# v1 : boule de rayon 2 percée d'un trou de rayon 1 : 4√3 π ; leurres 32π/3, 32π/3 − 2√3π, 8π.
V1 = 2 * pi * integrate(2 * sqrt(4 - r**2) * r, (r, 1, 2))
choix(V1, 4 * sqrt(3) * pi, [R_(32, 3) * pi, R_(32, 3) * pi - 2 * sqrt(3) * pi, 8 * pi], "aires-volumes-ex4 v1")
assert egal(R_(4, 3) * pi * sqrt(3) ** 3, 4 * sqrt(3) * pi)
# v2 : partie de la boule de rayon 2 dans le cylindre de rayon 1 : 4π(8 − 3√3)/3 ; leurres 4π, 32π/3, 2√3π.
V2 = 2 * pi * integrate(2 * sqrt(4 - r**2) * r, (r, 0, 1))
choix(V2, 4 * pi * (8 - 3 * sqrt(3)) / 3, [4 * pi, R_(32, 3) * pi, 2 * sqrt(3) * pi], "aires-volumes-ex4 v2")
assert egal(V1 + V2, R_(32, 3) * pi) and float(2 * sqrt(3) * pi) < float(V2) < float(4 * pi)
assert abs(float(V2) - 11.74) < 0.01

# int-l3-multiples-aires-volumes-ex5 : tore : 2π²Ra² ; leurres π²Ra², 4π²Ra², 2πRa².
sz = sqrt(ap**2 - z**2)
assert egal(pi * ((Rp + sz) ** 2 - (Rp - sz) ** 2), 4 * pi * Rp * sz)
Vt = integrate(4 * pi * Rp * sqrt(ap**2 - z**2), (z, -ap, ap))
assert egal(Vt, 2 * pi**2 * Rp * ap**2)
_n35 = {Rp: 3, ap: 1}
choix(Vt.subs(_n35), (2 * pi**2 * Rp * ap**2).subs(_n35),
      [(pi**2 * Rp * ap**2).subs(_n35), (4 * pi**2 * Rp * ap**2).subs(_n35), (2 * pi * Rp * ap**2).subs(_n35)], "aires-volumes-ex5")
# Guldin : 2π ∬_D x dx dz sur le disque de centre (R, 0) (en polaires centrées) = 2π · R · πa².
assert egal(2 * pi * integrate((Rp + r * cos(th)) * r, (r, 0, ap), (th, 0, 2 * pi)), 2 * pi**2 * Rp * ap**2)

# ===========================================================================
# int-l3-multiples-triples (carte)
# ===========================================================================
Jcyl = Matrix([r * cos(th), r * sin(th), z]).jacobian([r, th, z])
assert simplify(Jcyl.det()) == r
Jsph = Matrix([rho * sin(ph) * cos(th), rho * sin(ph) * sin(th), rho * cos(ph)]).jacobian([rho, ph, th])
assert simplify(Jsph.det() - rho**2 * sin(ph)) == 0
# Développement selon la dernière ligne (démonstration) : cos φ · ρ² sin φ cos φ + ρ sin φ · ρ sin² φ.
c31 = Jsph[2, 0] * Jsph.minor(2, 0)
c32 = -Jsph[2, 1] * Jsph.minor(2, 1)
assert egal(c31, cos(ph) * rho**2 * sin(ph) * cos(ph)) and egal(c32, rho * sin(ph) * rho * sin(ph) ** 2)
# Exemples : boule 4πR³/3 ; cône ∭ z = π/4.
assert egal(integrate(rho**2 * sin(ph), (rho, 0, Rp), (ph, 0, pi), (th, 0, 2 * pi)), R_(4, 3) * pi * Rp**3)
assert integrate(z * pi * z**2, (z, 0, 1)) == pi / 4
# Contre-exemples et pièges : latitude avec sin φ : 0 ; avec cos φ : 4π/3 ; φ ∈ [0, 2π] : 0 ; sans jacobien : 2π².
assert integrate(sin(ph), (ph, -pi / 2, pi / 2)) == 0 and integrate(sin(ph), (ph, 0, 2 * pi)) == 0
assert 2 * pi * R_(1, 3) * integrate(cos(ph), (ph, -pi / 2, pi / 2)) == R_(4, 3) * pi
assert integrate(1, (rho, 0, 1), (ph, 0, pi), (th, 0, 2 * pi)) == 2 * pi**2
# Trou 2, leurres (boule unité, volume 4π/3) : ρ² cos φ → 0 ; ρ sin φ → 2π ; ρ² → 2π²/3.
assert integrate(rho**2 * cos(ph), (rho, 0, 1), (ph, 0, pi), (th, 0, 2 * pi)) == 0
assert integrate(rho * sin(ph), (rho, 0, 1), (ph, 0, pi), (th, 0, 2 * pi)) == 2 * pi
assert integrate(rho**2, (rho, 0, 1), (ph, 0, pi), (th, 0, 2 * pi)) == 2 * pi**2 / 3
# Trou 1, leurres (cylindre r ≤ 1, 0 ≤ z ≤ 1, volume π) : 1 → 2π ; r² → 2π/3 ; rz → π/2.
assert integrate(r, (r, 0, 1), (th, 0, 2 * pi), (z, 0, 1)) == pi
assert integrate(1, (r, 0, 1), (th, 0, 2 * pi), (z, 0, 1)) == 2 * pi
assert integrate(r**2, (r, 0, 1), (th, 0, 2 * pi), (z, 0, 1)) == 2 * pi / 3
assert integrate(r * z, (r, 0, 1), (th, 0, 2 * pi), (z, 0, 1)) == pi / 2
# Piège « par piles » : h1 = 1, h2 = 0 sur [0,1]² : formule −1 pour un ensemble vide.
assert integrate(0 - 1, (x, 0, 1), (y, 0, 1)) == -1
# Affirmations : ∭_B z = 0 ; tétraèdre 1/6.
assert integrate(rho * cos(ph) * rho**2 * sin(ph), (rho, 0, 1), (ph, 0, pi), (th, 0, 2 * pi)) == 0
assert integrate(1 - x - y, (y, 0, 1 - x), (x, 0, 1)) == R_(1, 6)

# --- Exercices -------------------------------------------------------------
# int-l3-multiples-triples-ex1 : ∭_{[0,1]³} (x + yz) = 3/4 ; leurres 1/2, 1, 5/4.
choix(integrate(x + y * z, (x, 0, 1), (y, 0, 1), (z, 0, 1)), R_(3, 4), [R_(1, 2), 1, R_(5, 4)], "triples-ex1")
# v1 : ∭ xyz sur [0,1]×[0,2]×[0,3] = 9/2 ; leurres 9/4, 9, 3/2.
choix(integrate(x * y * z, (x, 0, 1), (y, 0, 2), (z, 0, 3)), R_(9, 2), [R_(9, 4), 9, R_(3, 2)], "triples-ex1 v1")
# v2 : ∭ x²yz sur [0,1]³ = 1/12 ; leurres 1/8, 1/24, 1/6.
choix(integrate(x**2 * y * z, (x, 0, 1), (y, 0, 1), (z, 0, 1)), R_(1, 12), [R_(1, 8), R_(1, 24), R_(1, 6)], "triples-ex1 v2")

# int-l3-multiples-triples-ex2 : ∭_T x = 1/24 ; leurres 1/6, 1/12, 1/8.
choix(integrate(x, (z, 0, 1 - x - y), (y, 0, 1 - x), (x, 0, 1)), R_(1, 24), [R_(1, 6), R_(1, 12), R_(1, 8)], "triples-ex2")
# v1 : volume de {x/2 + y + z ≤ 1} = 1/3 ; leurres 1/6, 2/3, 1.
choix(integrate(1, (z, 0, 1 - x / 2 - y), (y, 0, 1 - x / 2), (x, 0, 2)), R_(1, 3), [R_(1, 6), R_(2, 3), 1], "triples-ex2 v1")
# v2 : ∭_T xyz = 1/720 ; leurres 1/120, 1/48, 1/360 ; ∫_0^c y (c − y)²/2 dy = c⁴/24.
choix(integrate(x * y * z, (z, 0, 1 - x - y), (y, 0, 1 - x), (x, 0, 1)), R_(1, 720), [R_(1, 120), R_(1, 48), R_(1, 360)], "triples-ex2 v2")
assert egal(integrate(y * (ap - y) ** 2 / 2, (y, 0, ap)), ap**4 / 24)

# int-l3-multiples-triples-ex3 : ∭ (x²+y²) sur le cylindre r ≤ 1, 0 ≤ z ≤ 2 : π ; leurres 4π/3 (sans r), 2π, π/2.
choix(integrate(r**3, (r, 0, 1), (th, 0, 2 * pi), (z, 0, 2)), pi, [R_(4, 3) * pi, 2 * pi, pi / 2], "triples-ex3")
assert integrate(r**2, (r, 0, 1), (th, 0, 2 * pi), (z, 0, 2)) == R_(4, 3) * pi
# v1 : ∭ z sur x²+y² ≤ z ≤ 1 : π/3 ; leurres π/2, π/6, 2π/3.
choix(integrate(z * pi * z, (z, 0, 1)), pi / 3, [pi / 2, pi / 6, 2 * pi / 3], "triples-ex3 v1")
assert integrate(z * r, (z, r**2, 1), (r, 0, 1), (th, 0, 2 * pi)) == pi / 3
# v2 : ∭ z√(x²+y²) sur r ≤ 2, 0 ≤ z ≤ 1 : 8π/3 ; leurres 16π/3, 2π (sans r), 4π/3.
choix(integrate(z * r * r, (r, 0, 2), (th, 0, 2 * pi), (z, 0, 1)), R_(8, 3) * pi, [R_(16, 3) * pi, 2 * pi, R_(4, 3) * pi], "triples-ex3 v2")
assert integrate(z * r, (r, 0, 2), (th, 0, 2 * pi), (z, 0, 1)) == 2 * pi


# int-l3-multiples-triples-ex4 (et v1, v2) : boule de rayon R ; ∭ (x²+y²+z²) = 4πR⁵/5 ; chaque ∭ x_i² = 4πR⁵/15.
def boule(Rb, expr):
    sub = {x: rho * sin(ph) * cos(th), y: rho * sin(ph) * sin(th), z: rho * cos(ph)}
    return integrate(simplify(expr.subs(sub, simultaneous=True)) * rho**2 * sin(ph), (rho, 0, Rb), (ph, 0, pi), (th, 0, 2 * pi))


assert boule(1, x**2 + y**2 + z**2) == R_(4, 5) * pi
assert boule(1, x**2) == boule(1, y**2) == boule(1, z**2) == R_(4, 15) * pi
choix(boule(1, x**2 + y**2), R_(8, 15) * pi, [R_(4, 5) * pi, R_(4, 15) * pi, R_(2, 5) * pi], "triples-ex4")
assert boule(2, x**2 + y**2 + z**2) == R_(128, 5) * pi
choix(boule(2, x**2 + y**2), R_(256, 15) * pi, [R_(128, 5) * pi, R_(128, 15) * pi, R_(64, 15) * pi], "triples-ex4 v1")
choix(boule(1, x**2 + 2 * y**2 + 3 * z**2), R_(8, 5) * pi, [R_(12, 5) * pi, R_(4, 5) * pi, R_(24, 5) * pi], "triples-ex4 v2")

# int-l3-multiples-triples-ex5 : ‖M − P‖² = ρ² + d² − 2ρd cos φ ; ∫_0^π sin φ/√(…) dφ = 2/d ; U(d) = 4π/(3d) ; U(2) = 2π/3.
Msph = (rho * sin(ph) * cos(th), rho * sin(ph) * sin(th), rho * cos(ph))
assert egal(Msph[0] ** 2 + Msph[1] ** 2 + (Msph[2] - dp) ** 2, rho**2 + dp**2 - 2 * rho * dp * cos(ph))
for rr, dd in ((R_(1, 3), 2), (R_(1, 2), 3), (1, R_(3, 2)), (R_(9, 10), R_(11, 10))):
    val = mp.quad(lambda p: mp.sin(p) / mp.sqrt(rr**2 + dd**2 - 2 * rr * dd * mp.cos(p)), [0, mp.pi])
    assert abs(val - mp.mpf(2) / dd) < 1e-15
U2 = 2 * pi * integrate(rho**2 * R_(2, 2), (rho, 0, 1))  # 2/d avec d = 2
choix(U2, 2 * pi / 3, [R_(4, 3) * pi, pi / 3, pi / 6], "triples-ex5")
# Contrôle numérique direct de U(2) (intégrale triple en sphériques, axe Oz).
U2num = mp.quad(lambda p_, q_: 2 * mp.pi * p_**2 * mp.sin(q_) / mp.sqrt(p_**2 + 4 - 4 * p_ * mp.cos(q_)), [0, 1], [0, mp.pi])
assert abs(U2num - 2 * mp.pi / 3) < 1e-15

# QCM q7 : cornet φ ∈ [0, π/4] : π(2 − √2)/3 ; leurres 2π/3, π√2/3, π/3.
choix(integrate(rho**2 * sin(ph), (rho, 0, 1), (ph, 0, pi / 4), (th, 0, 2 * pi)), pi * (2 - sqrt(2)) / 3,
      [2 * pi / 3, pi * sqrt(2) / 3, pi / 3], "q7")

print("int_l3_multiples : toutes les vérifications passent.")
