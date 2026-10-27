# Vérification calculatoire du chapitre « geo-l1-complexes » (géométrie, L1 : complexes et géométrie plane).
# Fichiers : content/cartes/geometrie/l1-geo-l1-complexes.json, content/exercices/geometrie/l1-geo-l1-complexes.json.
# Chaque assert cite l'identifiant de la carte ou de l'exercice (et la variante : v1, v2) concerné.
# Les réponses justes sont vérifiées en calcul exact (SymPy) ; chaque leurre est vérifié FAUX.
import cmath
import math
import random

import sympy as sp

I, pi, sqrt = sp.I, sp.pi, sp.sqrt
z, Z, t = sp.symbols("z Z t")
th = sp.symbols("theta", real=True)
random.seed(1)


def egal(a, b):
    """Égalité exacte de deux complexes (expressions SymPy)."""
    d = sp.simplify(sp.expand_complex(sp.sympify(a) - sp.sympify(b)))
    return d == 0


def racines_carrees(w):
    return [sp.nsimplify(r) for r in sp.solve(sp.Eq(z**2, w), z)]


def meme_ensemble(A, B):
    A, B = list(A), list(B)
    return all(any(egal(a, b) for b in B) for a in A) and all(any(egal(a, b) for a in A) for b in B)


def solutions(P):
    return [r for r in sp.roots(sp.Poly(sp.expand(P), z)).keys()]


def identite_trigo(f, g, n=12):
    """f(θ) = g(θ) pour tout θ réel : réécriture en exponentielles (exacte), sinon contrôle numérique en n points."""
    d = sp.expand((f - g).rewrite(sp.exp))
    if sp.simplify(d) == 0:
        return True
    h = sp.lambdify(th, f - g, "cmath")
    return all(abs(h(random.uniform(-7, 7))) < 1e-10 for _ in range(n))


def jamais_identique(f, g):
    """f et g diffèrent en au moins un point (un leurre de linéarisation est faux)."""
    return any(abs(complex(sp.N((f - g).subs(th, v)))) > 1e-9 for v in [0, 0.3, 1, pi / 2, 2.1])


# ======================================================================================
# Carte geo-l1-complexes-racines-carrees
# ======================================================================================
C = "geo-l1-complexes-racines-carrees"
# Énoncé : formules algébriques x^2 = (|w|+u)/2, y^2 = (|w|-u)/2, 2xy = v, testées sur des w aléatoires.
for _ in range(20):
    u, v = random.randint(-9, 9), random.randint(-9, 9)
    if u == v == 0:
        continue
    w = u + I * v
    for r in racines_carrees(w):
        x, y = sp.re(r), sp.im(r)
        assert egal(x**2, (sp.Abs(w) + u) / 2) and egal(y**2, (sp.Abs(w) - u) / 2) and egal(2 * x * y, v)  # carte, calcul algébrique
    assert len(racines_carrees(w)) == 2  # carte : deux racines carrées pour w ≠ 0
    r1, r2 = racines_carrees(w)
    assert egal(r1 + r2, 0)  # carte : elles sont opposées
assert racines_carrees(0) == [0]  # carte : 0 n'a que la racine carrée 0
# Forme exponentielle : ±sqrt(r) e^{iθ/2}.
for r_, t_ in [(4, pi / 3), (2, -3 * pi / 4), (9, 5 * pi / 6)]:
    assert egal((sqrt(r_) * sp.exp(I * t_ / 2)) ** 2, r_ * sp.exp(I * t_))  # carte, forme exponentielle
# Équation à coefficients complexes : formule (-b ± δ)/(2a), test aléatoire.
for _ in range(10):
    a_, b_, c_ = [random.randint(-5, 5) + I * random.randint(-5, 5) for _ in range(3)]
    if a_ == 0:
        continue
    D = b_**2 - 4 * a_ * c_
    d_ = racines_carrees(D)[0]
    z1, z2 = (-b_ + d_) / (2 * a_), (-b_ - d_) / (2 * a_)
    assert egal(a_ * z1**2 + b_ * z1 + c_, 0) and egal(a_ * z2**2 + b_ * z2 + c_, 0)  # carte, équation
    assert egal(sp.expand(a_ * (z - z1) * (z - z2)) - (a_ * z**2 + b_ * z + c_), 0)  # carte, factorisation
# Cours : racines carrées de 3+4i = ±(2+i).
assert meme_ensemble(racines_carrees(3 + 4 * I), [2 + I, -2 - I])  # carte (cours)
# Exemple : -7+24i, |w| = 25, racines ±(3+4i) ; équation z^2-3z+3-i = 0, Δ = -3+4i = (1+2i)^2, solutions 2+i, 1-i.
assert sp.Abs(-7 + 24 * I) == 25 and meme_ensemble(racines_carrees(-7 + 24 * I), [3 + 4 * I, -3 - 4 * I])  # carte (exemple)
assert egal(9 - 4 * (3 - I), -3 + 4 * I) and egal((1 + 2 * I) ** 2, -3 + 4 * I)  # carte (exemple, Δ et δ)
assert meme_ensemble(solutions(z**2 - 3 * z + 3 - I), [2 + I, 1 - I])  # carte (exemple, contre-exemple, affirmation 2)
assert not egal(sp.conjugate(2 + I), 1 - I)  # carte (contre-exemple) : solutions non conjuguées
# Trous : leurres faux.
assert len(racines_carrees(4)) == 2  # trou 1 : ni une seule, ni quatre, ni une infinité (w = 4)
assert not egal(sp.conjugate(1 + I), -1 - I) and meme_ensemble(racines_carrees(2 * I), [1 + I, -1 - I])  # trou 2, « conjuguées » (w = 2i)
assert not egal((1 + I) * (-1 - I), 1)  # trou 2, « inverses » (w = 2i)
assert sp.Abs(racines_carrees(4)[0]) != 4  # trou 2, « de même module » (w = 4)
# Trou 3 : sur z^2 + 1 = 0 (a=1, b=0, c=1, Δ=-4, δ=2i), les leurres ne donnent pas les solutions ±i.
P0 = z**2 + 1
assert not egal(P0.subs(z, (0 + 2 * I) / 1), 0)  # trou 3, leurre (-b±δ)/a
assert not egal((z**2 + 2 * z).subs(z, (2 + 0) / 2), 0)  # trou 3, leurre (b±δ)/(2a) sur z^2+2z (b=2, δ=2) : donne 2 et 0 au lieu de 0 et -2
assert not egal(P0.subs(z, 2 / 2), 0)  # trou 3, leurre sqrt(|Δ|) : donne ±1
# Contre-exemple : racine carrée multiplicative impossible (s(1) = 1 mais s(-1)^2 = -1).
assert egal((I) ** 2, -1) and egal((-I) ** 2, -1)  # contre-exemple : s(-1) ∈ {i, -i}, et s(-1)^2 = -1 ≠ 1
# Affirmation 1 : z^2-(1+i)z+i = (z-1)(z-i), Δ = -2i non réel, 1 solution réelle.
assert egal(sp.expand((z - 1) * (z - I)) - (z**2 - (1 + I) * z + I), 0) and egal((1 + I) ** 2 - 4 * I, -2 * I)  # affirmation 1
# Affirmation 4 : racines carrées de i.
assert meme_ensemble(racines_carrees(I), [sqrt(2) / 2 * (1 + I), -sqrt(2) / 2 * (1 + I)])  # affirmation 4

# --- Exercice ex1 et variantes : racines carrées.
for nom, w, bonne, leurres in [
    ("ex1", 5 - 12 * I, [3 - 2 * I, -3 + 2 * I], [[3 + 2 * I, -3 - 2 * I], [2 - 3 * I, -2 + 3 * I], [3 - 2 * I, 3 + 2 * I]]),
    ("ex1 v1", 8 + 6 * I, [3 + I, -3 - I], [[3 - I, -3 + I], [1 + 3 * I, -1 - 3 * I], [3 + I, 3 - I]]),
    ("ex1 v2", -15 - 8 * I, [1 - 4 * I, -1 + 4 * I], [[1 + 4 * I, -1 - 4 * I], [4 - I, -4 + I], [1 - 4 * I, 1 + 4 * I]]),
]:
    assert meme_ensemble(racines_carrees(w), bonne), nom  # racines-carrees-ex1 (et variantes) : réponse
    for L in leurres:
        assert not meme_ensemble(racines_carrees(w), L), nom  # leurre faux
assert sp.Abs(5 - 12 * I) == 13 and sp.Abs(8 + 6 * I) == 10 and sp.Abs(-15 - 8 * I) == 17  # ex1 : modules cités

# --- Exercice ex2 et variantes : équations du second degré.
for nom, P, D, d, bonne, leurres in [
    ("ex2", z**2 - (5 + I) * z + 8 + 4 * I, -8 - 6 * I, 1 - 3 * I, [3 - I, 2 + 2 * I],
     [[3 + I, 2 - 2 * I], [-3 + I, -2 - 2 * I], [6 - 2 * I, 4 + 4 * I]]),
    ("ex2 v1", z**2 - (4 + 2 * I) * z + 3 + 6 * I, -8 * I, 2 - 2 * I, [3, 1 + 2 * I],
     [[3, 1 - 2 * I], [-3, -1 - 2 * I], [6, 2 + 4 * I]]),
    ("ex2 v2", z**2 - (1 + 2 * I) * z - 3 + I, 9, 3, [2 + I, -1 + I],
     [[2 - I, -1 - I], [-2 - I, 1 - I], [4 + 2 * I, -2 + 2 * I]]),
]:
    Pp = sp.Poly(P, z)
    a_, b_, c_ = Pp.all_coeffs()
    assert egal(b_**2 - 4 * a_ * c_, D) and egal(d**2, D), nom  # racines-carrees-ex2 : Δ et δ cités
    assert meme_ensemble(solutions(P), bonne), nom  # racines-carrees-ex2 : réponse
    for L in leurres:
        assert not meme_ensemble(solutions(P), L), nom  # leurre faux
assert egal((5 + I) ** 2, 24 + 10 * I) and egal((4 + 2 * I) ** 2, 12 + 16 * I)  # ex2, ex2 v1 : calculs intermédiaires

# --- Exercice ex3 et variantes : cos(π/8), cos(3π/8), cos(π/12).
for nom, w, r_, ang, xx, yy, bonne, leurres in [
    ("ex3", 1 + I, sqrt(2), pi / 4, (sqrt(2) + 1) / 2, (sqrt(2) - 1) / 2, sqrt(2 + sqrt(2)) / 2,
     [sqrt(2 - sqrt(2)) / 2, sqrt(2 + sqrt(2)) / 4, sqrt(1 + sqrt(2)) / 2]),
    ("ex3 v1", -1 + I, sqrt(2), 3 * pi / 4, (sqrt(2) - 1) / 2, (sqrt(2) + 1) / 2, sqrt(2 - sqrt(2)) / 2,
     [sqrt(2 + sqrt(2)) / 2, sqrt(2 - sqrt(2)) / 4, -sqrt(2 - sqrt(2)) / 2]),
    ("ex3 v2", sqrt(3) + I, 2, pi / 6, (2 + sqrt(3)) / 2, (2 - sqrt(3)) / 2, sqrt(2 + sqrt(3)) / 2,
     [sqrt(2 - sqrt(3)) / 2, sqrt(2 + sqrt(3)) / 4, sqrt(3) / 2]),
]:
    assert egal(sp.Abs(w), r_) and egal(r_ * sp.exp(I * ang), w), nom  # racines-carrees-ex3 : forme exponentielle
    z0 = sqrt(xx) + I * sqrt(yy)
    assert egal(sp.expand(z0**2), w), nom  # ex3 : z0 est une racine carrée (x, y > 0, xy > 0)
    assert egal(sqrt(r_) * sp.exp(I * ang / 2), z0), nom  # ex3 : identification des deux écritures
    assert egal(sp.cos(ang / 2), bonne), nom  # ex3 : réponse
    for L in leurres:
        assert not egal(sp.cos(ang / 2), L), nom  # leurre faux
    assert egal(sqrt(xx) * sqrt(yy), sp.Rational(1, 2)), nom  # ex3 : contrôle 2xy = 1

# --- Exercice ex4 et variantes : équations bicarrées.
for nom, c2, c0, Zs, sols, bonne, leurres in [
    ("ex4", 6, 25, [-3 + 4 * I, -3 - 4 * I], [1 + 2 * I, 1 - 2 * I, -1 + 2 * I, -1 - 2 * I],
     (z**2 - 2 * z + 5) * (z**2 + 2 * z + 5),
     [(z**2 - 4 * z + 5) * (z**2 + 4 * z + 5), (z**2 - 2 * z - 5) * (z**2 + 2 * z - 5), (z**2 - 2 * z + 5) ** 2]),
    ("ex4 v1", -6, 25, [3 + 4 * I, 3 - 4 * I], [2 + I, 2 - I, -2 + I, -2 - I],
     (z**2 - 4 * z + 5) * (z**2 + 4 * z + 5),
     [(z**2 - 2 * z + 5) * (z**2 + 2 * z + 5), (z**2 - 4 * z - 5) * (z**2 + 4 * z - 5), (z**2 - 4 * z + 5) ** 2]),
    ("ex4 v2", 10, 169, [-5 + 12 * I, -5 - 12 * I], [2 + 3 * I, 2 - 3 * I, -2 + 3 * I, -2 - 3 * I],
     (z**2 - 4 * z + 13) * (z**2 + 4 * z + 13),
     [(z**2 - 6 * z + 13) * (z**2 + 6 * z + 13), (z**2 - 4 * z - 13) * (z**2 + 4 * z - 13), (z**2 - 4 * z + 13) ** 2]),
]:
    P = z**4 + c2 * z**2 + c0
    assert meme_ensemble([r for r in sp.roots(sp.Poly(Z**2 + c2 * Z + c0, Z)).keys()], Zs), nom  # ex4 : question 1
    assert meme_ensemble(solutions(P), sols), nom  # ex4 : question 2
    assert sp.expand(bonne - P) == 0, nom  # ex4 : réponse (factorisation)
    for L in leurres:
        assert sp.expand(L - P) != 0, nom  # leurre faux
assert 36 - 100 == -64 and 100 - 676 == -576 and 24**2 == 576  # ex4 : discriminants cités

# --- Exercice ex5 : s(w) = sqrt|w| (w+|w|)/|w+|w||.
def s_princ(w):
    return cmath.sqrt(abs(w)) * (w + abs(w)) / abs(w + abs(w))


for _ in range(200):
    w = complex(random.uniform(-5, 5), random.uniform(-5, 5))
    sw = s_princ(w)
    assert abs(sw**2 - w) < 1e-9 and sw.real > 0  # racines-carrees-ex5, question 3 (contrôle numérique)
    assert abs(sw - cmath.sqrt(w)) < 1e-9  # ex5 : s est la racine principale
# Question 2 : w + |w| = 2r cos(θ/2) e^{iθ/2}.
r, tt = sp.symbols("r t", positive=True)
assert egal(sp.expand_complex(r * sp.exp(I * tt) + r - 2 * r * sp.cos(tt / 2) * sp.exp(I * tt / 2)), 0)  # ex5, question 2
# Question 4 : w = w' = e^{2iπ/3} : s(w)s(w') = -s(ww').
w = cmath.exp(2j * math.pi / 3)
assert abs(s_princ(w) * s_princ(w) + s_princ(w * w)) < 1e-12  # ex5 : bonne réponse
assert abs(s_princ(w) - cmath.exp(1j * math.pi / 3)) < 1e-12 and abs(s_princ(w * w) - cmath.exp(-1j * math.pi / 3)) < 1e-12  # ex5 : valeurs citées
# Leurre « Oui » : faux (contre-exemple ci-dessus). Leurre w = w' = i : ww' = -1 n'est pas dans D.
assert abs(1j * 1j - (-1)) < 1e-15  # ex5, leurre 4 : i·i = -1 ∈ R_-, hors du domaine


# ======================================================================================
# Carte geo-l1-complexes-triangulaire
# ======================================================================================
C = "geo-l1-complexes-triangulaire"
# Énoncé : inégalité, cas d'égalité, inégalité renversée (contrôles numériques).
for _ in range(500):
    a_ = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    b_ = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    assert abs(a_ + b_) <= abs(a_) + abs(b_) + 1e-12  # carte, inégalité triangulaire
    assert abs(abs(a_) - abs(b_)) <= abs(a_ - b_) + 1e-12  # carte, inégalité renversée
    lam = random.uniform(0, 3)
    assert abs(abs(lam * b_ + b_) - (abs(lam * b_) + abs(b_))) < 1e-9  # carte, cas d'égalité (λ ≥ 0)
# Identité |z+w|^2 = |z|^2 + |w|^2 + 2 Re(z w̄).
x1, y1, x2, y2 = sp.symbols("x1 y1 x2 y2", real=True)
z1s, z2s = x1 + I * y1, x2 + I * y2
assert sp.expand(sp.Abs(z1s + z2s) ** 2 - (sp.Abs(z1s) ** 2 + sp.Abs(z2s) ** 2 + 2 * sp.re(z1s * sp.conjugate(z2s)))) == 0  # carte (démonstration)
# Contre-exemples et affirmations.
assert abs(1 + (-1)) == 0 < 2  # carte (contre-exemple λ < 0 ; affirmation 1 ; trou 2 leurres 1 et 2)
assert egal(sp.Abs(1 + I), sqrt(2)) and sqrt(2) < 2  # carte (contre-exemple même module ; trou 2 leurre 3)
assert abs(1 + 1 - 1) == 1 < 3  # carte (contre-exemple n termes)
assert 2 > math.sqrt(2)  # affirmation 3 : |Re|+|Im| = 2 > sqrt(2) pour z = 1+i
# Exemple : |z^2+3z+1| ≤ 5 sur le disque unité, égalité seulement en z = 1 (contrôle sur une grille).
M = max(abs(complex(rr * math.cos(a), rr * math.sin(a)) ** 2 + 3 * complex(rr * math.cos(a), rr * math.sin(a)) + 1)
        for rr in [k / 40 for k in range(41)] for a in [2 * math.pi * k / 720 for k in range(720)])
assert abs(M - 5) < 1e-12  # carte (exemple) : maximum 5 atteint (en z = 1)
assert all(abs(complex(math.cos(a), math.sin(a)) ** 2 + 3 * complex(math.cos(a), math.sin(a)) + 1) < 5 - 1e-6
           for a in [2 * math.pi * k / 720 for k in range(1, 720)])  # carte (exemple) : seulement en z = 1
# Trous : leurres faux (z = w = 1, z = 1 et w = -1, z = 0 et w = 1).
assert not abs(1 + 1) <= abs(abs(1) - abs(1))  # trou 1, leurre |z+w| ≤ ||z|-|w|| (z = w = 1)
assert not abs(1 + 1) ** 2 <= 1 + 1  # trou 1, leurre |z+w|^2 ≤ |z|^2+|w|^2 (z = w = 1)
assert not abs(1 + 1j) >= 1 + 1  # trou 1, leurre ≥ (z = 1, w = i)
assert not abs(1 - (-1)) <= abs(abs(1) - abs(-1))  # trou 3, leurre |z-w| ≤ ||z|-|w|| (z = 1, w = -1)
assert not abs(0 - 1) <= 0 - 1  # trou 3, leurre |z-w| ≤ |z|-|w| (z = 0, w = 1)
assert not abs(1 - 1) >= 1 + 1  # trou 3, leurre |z-w| ≥ |z|+|w| (z = w = 1)

# --- ex1 et variantes : min et max de |z + w| pour |z| = r.
for nom, rr, w, mini, maxi, leurres in [
    ("ex1", 3, 4 - 3 * I, 2, 8, [(0, 8), (3, 8), (2, sqrt(34))]),
    ("ex1 v1", 5, 12 + 5 * I, 8, 18, [(0, 18), (5, 18), (8, sqrt(194))]),
    ("ex1 v2", 7, -3 + 4 * I, 2, 12, [(0, 12), (7, 12), (2, sqrt(74))]),
]:
    W = sp.Abs(w)
    assert abs(rr - W) == mini and rr + W == maxi, nom  # triangulaire-ex1 : bornes
    zmax, zmin = sp.Rational(rr) / W * w, -sp.Rational(rr) / W * w
    assert egal(sp.Abs(zmax), rr) and egal(sp.Abs(zmax + w), maxi) and egal(sp.Abs(zmin), rr) and egal(sp.Abs(zmin + w), mini), nom  # atteintes
    vals = [abs(complex(rr * math.cos(a), rr * math.sin(a)) + complex(w)) for a in [2 * math.pi * k / 3600 for k in range(3600)]]
    assert abs(min(vals) - mini) < 1e-3 and abs(max(vals) - maxi) < 1e-3, nom  # contrôle numérique
    for lo, hi in leurres:
        assert (lo, hi) != (mini, maxi) and not (egal(lo, mini) and egal(hi, maxi)), nom  # leurre faux
assert egal(sp.Rational(3, 5) * (4 - 3 * I), (12 - 9 * I) / 5) and egal(-sp.Rational(7, 5) * (-3 + 4 * I), (21 - 28 * I) / 5)  # ex1, ex1 v2 : points cités

# --- ex2 et variantes : |z-a| + |z-b| = |a-b| ⇔ segment [AB].
for nom, a_, b_, d in [("ex2", 1, 3 * I, sqrt(10)), ("ex2 v1", -2, 2 * I, 2 * sqrt(2)), ("ex2 v2", 1 + I, -3 - 2 * I, 5)]:
    assert egal(sp.Abs(b_ - a_), d), nom  # triangulaire-ex2 : AB cité
    A, B = complex(a_), complex(b_)
    D = float(d)
    for k in range(11):
        P_ = A + k / 10 * (B - A)
        assert abs(abs(P_ - A) + abs(P_ - B) - D) < 1e-12, nom  # segment ⊂ E
    # Hors du segment : sur la droite (AB) prolongée, sur le cercle de diamètre [AB] (hors A, B) : somme > AB.
    assert abs(A + 1.5 * (B - A) - A) + abs(A + 1.5 * (B - A) - B) > D + 1e-9, nom  # leurre « droite (AB) »
    Mi = (A + B) / 2
    Q = Mi + (B - A) / 2 * 1j
    assert abs(Q - A) + abs(Q - B) > D + 1e-9, nom  # leurre « cercle de diamètre [AB] »
    assert abs(Mi - A) + abs(Mi - B) - D < 1e-12, nom  # leurre « seuls A et B » : le milieu convient aussi
    for _ in range(300):
        P_ = complex(random.uniform(-6, 6), random.uniform(-6, 6))
        if abs(abs(P_ - A) + abs(P_ - B) - D) < 1e-12:
            tpar = ((P_ - A) / (B - A))
            assert abs(tpar.imag) < 1e-9 and -1e-9 <= tpar.real <= 1 + 1e-9  # E ⊂ segment (rarement testé)

# --- ex3 et variantes : |P(z)| ≥ m sur le disque unité, égalité seulement en z = 1.
for nom, P, m, cube in [("ex3", z**3 + 2 * z - 6, 3, False), ("ex3 v1", z**4 + 3 * z - 8, 4, True), ("ex3 v2", z**5 + 4 * z**2 - 10, 5, True)]:
    f = sp.lambdify(z, P, "math")
    f = sp.lambdify(z, P)
    pts = [complex(rr * math.cos(a), rr * math.sin(a)) for rr in [k / 50 for k in range(51)] for a in [2 * math.pi * k / 900 for k in range(900)]]
    vals = [abs(f(q)) for q in pts]
    assert min(vals) >= m - 1e-12, nom  # triangulaire-ex3 : minoration
    assert abs(abs(f(1)) - m) < 1e-12, nom  # égalité en z = 1
    assert all(abs(f(q)) > m + 1e-9 for q in pts if abs(q - 1) > 1e-3), nom  # seulement en z = 1 (bonne réponse)
    assert abs(abs(f(-1)) - m) > 1e-6, nom  # leurres « z = -1 » faux
    assert abs(abs(f(1j)) - m) > 1e-6, nom  # leurre « tous les z de module 1 » faux
    j = cmath.exp(2j * math.pi / 3)
    assert abs(abs(f(j)) - m) > 1e-6 and abs(abs(f(j.conjugate())) - m) > 1e-6, nom  # leurre « racines cubiques de l'unité » faux
assert sp.lambdify(z, z**3 + 2 * z)(-1) == -3 and sp.lambdify(z, z**3 + 2 * z)(1) == 3  # ex3 : u = ±3 cité
assert abs(abs(4 * cmath.exp(2j * math.pi / 3) - 8) - 4 * math.sqrt(7)) < 1e-12  # ex3 v1 : |P(j)| = 4 sqrt 7 (contrôle)

# --- ex4 et variantes : triplets (a, b, c).
for nom, ra, rb, rc, S, sol, leurres in [
    ("ex4", 1, 2, 3, 6 * I, (I, 2 * I, 3 * I), [(3 * I, 2 * I, I)]),
    ("ex4 v1", 2, 3, 5, 6 + 8 * I, ((6 + 8 * I) / 5, (9 + 12 * I) / 5, 3 + 4 * I), [(2, 3, 5)]),
    ("ex4 v2", 1, 1, 2, -4, (-1, -1, -2), [(1, 1, 2)]),
]:
    assert sp.Abs(S) == ra + rb + rc, nom  # triangulaire-ex4 : |S| = somme des bornes (d'où les égalités)
    u = S / sp.Abs(S)
    assert all(egal(x_, r_ * u) for x_, r_ in zip(sol, (ra, rb, rc))), nom  # réponse : z_k = |z_k| u
    assert egal(sum(sol), S) and all(egal(sp.Abs(x_), r_) for x_, r_ in zip(sol, (ra, rb, rc))), nom  # le triplet convient
    for L in leurres:
        ok = egal(sum(L), S) and all(sp.Abs(x_) <= r_ for x_, r_ in zip(L, (ra, rb, rc)))
        assert not ok, nom  # leurre : ce triplet ne convient pas (donc l'ensemble proposé est faux)
    # « l'ensemble vide » et « un ensemble infini » sont faux : il existe une solution et elle est unique (question 1).

# ======================================================================================
# Carte geo-l1-complexes-similitudes
# ======================================================================================
C = "geo-l1-complexes-similitudes"
a, b = sp.symbols("a b")
# Énoncé : point fixe b/(1-a), f(z) - ω = a(z - ω), composée, réciproque, deux points.
w0 = b / (1 - a)
assert sp.simplify(a * w0 + b - w0) == 0  # carte : ω est fixe
assert sp.simplify((a * z + b - w0) - a * (z - w0)) == 0  # carte : f(z) - ω = a(z - ω)
a2, b2 = sp.symbols("a2 b2")
assert sp.expand(a2 * (a * z + b) + b2 - (a2 * a * z + a2 * b + b2)) == 0  # carte : composée
assert sp.simplify(a * (z / a - b / a) + b - z) == 0  # carte : réciproque
al, be, al2, be2 = sp.symbols("alpha beta alpha2 beta2")
sol = sp.solve([a * al + b - al2, a * be + b - be2], [a, b], dict=True)[0]
assert sp.simplify(sol[a] - (be2 - al2) / (be - al)) == 0 and sp.simplify(sol[b] - (al2 - sol[a] * al)) == 0  # carte : deux points
# Trous : leurres faux sur a = 2, b = 1 (point fixe -1) et a = i (angle π/2, rapport 1).
fA = lambda q: 2 * q + 1
assert fA(1 / (1 - 2)) == 1 / (1 - 2)  # trou 1 : b/(1-a) = -1 est fixe
for cand in [1 / (2 - 1), -1 / 2, 1 / (1 + 2)]:
    assert fA(cand) != cand  # trou 1 : leurres non fixes
assert egal(sp.arg(I), pi / 2) and not egal(-sp.arg(I), pi / 2) and not egal(sp.arg(1 - I), pi / 2) and not egal(sp.arg(1), pi / 2)  # trou 2 (a = i, b = 1)
assert sp.Abs(2 * I) == 2 and sp.Rational(1, 2) != 2 and sp.Abs(1) != 2  # trou 3 (a = 2i, b = 1) : 1/|a| et |b| faux
# Exemple : f(z) = (1+i)z + 2 - i, centre 1+2i, rapport sqrt 2, angle π/4.
assert egal((2 - I) / (1 - (1 + I)), 1 + 2 * I) and egal((1 + I) * (1 + 2 * I) + 2 - I, 1 + 2 * I)  # carte (exemple)
assert egal(sp.Abs(1 + I), sqrt(2)) and egal(sp.arg(1 + I), pi / 4)  # carte (exemple)
# Contre-exemples : z ↦ z̄ renverse l'angle (quotients i et -i) ; z ↦ z+1 sans point fixe.
assert egal((I - 0) / (1 - 0), I) and egal((sp.conjugate(I) - 0) / (1 - 0), -I)  # contre-exemple, affirmation 3
assert sp.solve(sp.Eq(z + 1, z), z) == []  # contre-exemple, affirmation 1
# Affirmation 4 : réciproque de rapport 1/k et d'angle -θ.
k_, t_ = sp.symbols("k t", positive=True)
assert sp.simplify(1 / (k_ * sp.exp(I * t_)) - sp.exp(-I * t_) / k_) == 0  # affirmation 4

# --- ex1 et variantes.
def elements(a_, b_):
    return sp.Abs(a_), sp.arg(a_), sp.simplify(b_ / (1 - a_))


for nom, a_, b_, bonne, leurres in [
    ("ex1", 1 + I * sqrt(3), sqrt(3) - I * sqrt(3), (2, pi / 3, 1 + I), [(2, pi / 6, 1 + I), (4, pi / 3, 1 + I), (2, pi / 3, -1 - I)]),
    ("ex1 v1", -3 * I, 5 + 5 * I, (3, -pi / 2, 2 - I), [(3, pi / 2, 2 - I), (3, -pi / 2, -1 + 2 * I), (-3, pi / 2, 2 - I)]),
    ("ex1 v2", 1 - I, -3, (sqrt(2), -pi / 4, 3 * I), [(sqrt(2), pi / 4, 3 * I), (2, -pi / 4, 3 * I), (sqrt(2), -pi / 4, -3 * I)]),
]:
    k_, t_, w_ = elements(a_, b_)
    assert egal(k_, bonne[0]) and egal(t_, bonne[1]) and egal(w_, bonne[2]), nom  # similitudes-ex1 : réponse
    assert egal(a_ * w_ + b_, w_), nom  # ex1 : vérification du point fixe
    for L in leurres:
        assert not (egal(k_, L[0]) and egal(t_, L[1]) and egal(w_, L[2])), nom  # leurre faux
assert egal((5 + 5 * I) / (1 - 3 * I), -1 + 2 * I)  # ex1 v1 : le leurre -1+2i provient de b/(1+a) au lieu de b/(1-a)

# --- ex2 et variantes.
for nom, A, B, A2, B2, bonne, centre, leurres in [
    ("ex2", 1, I, 2 + I, I, (1 + I, 1), I, [(1 - I, 1 + 2 * I), (-(1 + I), 3 + 2 * I), (1 + I, -1)]),
    ("ex2 v1", 0, 2, 1 + I, 1 + 3 * I, (I, 1 + I), I, [(-I, 1 + I), (I, 1 + 3 * I), (2 * I, 1 + I)]),
    ("ex2 v2", 1 + I, 3 + I, -1, -1 + 4 * I, (2 * I, 1 - 2 * I), 1, [(-2 * I, -3 + 2 * I), (I, -I), (2 * I, -1)]),
]:
    a_, b_ = bonne
    assert egal(a_ * A + b_, A2) and egal(a_ * B + b_, B2), nom  # similitudes-ex2 : réponse
    assert egal(b_ / (1 - a_), centre), nom  # ex2 : centre
    for la, lb in leurres:
        assert not (egal(la * A + lb, A2) and egal(la * B + lb, B2)), nom  # leurre faux
assert egal(2 / (1 - I), 1 + I) and egal((1 + I) / (1 - I), I)  # ex2, ex2 v1 : calculs cités

# --- ex3 et variantes : suites z_{n+1} = a z_n + b.
for nom, a_, b_, N, w_, k_, t_, bonne, leurres in [
    ("ex3", (1 + I) / 2, 1, 8, 1 + I, sqrt(2) / 2, pi / 4, sp.Rational(15, 16) * (1 + I),
     [sp.Rational(1, 16) * (1 + I), sp.Rational(17, 16) * (1 + I), sp.Rational(15, 16) * (1 - I)]),
    ("ex3 v1", (1 + I * sqrt(3)) / 4, (3 - I * sqrt(3)) / 2, 6, 2, sp.Rational(1, 2), pi / 3, sp.Rational(63, 32),
     [sp.Rational(1, 32), sp.Rational(65, 32), sp.Rational(63, 64)]),
    ("ex3 v2", I / 3, 1 + 3 * I, 4, 3 * I, sp.Rational(1, 3), pi / 2, sp.Rational(80, 27) * I,
     [I / 27, sp.Rational(82, 27) * I, sp.Rational(80, 81) * I]),
]:
    assert egal(b_ / (1 - a_), w_) and egal(sp.Abs(a_), k_) and egal(sp.arg(a_), t_), nom  # similitudes-ex3 : question 1
    zn = sp.Integer(0)
    for n in range(N):
        zn = sp.expand(a_ * zn + b_)
        assert egal(zn - w_, a_ ** (n + 1) * (0 - w_)), nom  # ex3 : question 2
    assert egal(zn, bonne), nom  # ex3 : réponse
    for L in leurres:
        assert not egal(zn, L), nom  # leurre faux
    assert k_ < 1, nom  # ex3 : |z_n - ω| → 0
assert egal(((1 + I) / 2) ** 8, sp.Rational(1, 16)) and egal(((1 + I * sqrt(3)) / 4) ** 6, sp.Rational(1, 64)) and egal((I / 3) ** 4, sp.Rational(1, 81))  # ex3 : a^N cités

# --- ex4 et variantes : similitude qui permute A, B, C.
al, be, a = sp.symbols("alpha beta a")
for aa in [sp.exp(2 * I * pi / 3), sp.exp(-2 * I * pi / 3)]:
    assert egal(1 + aa + aa**2, 0)  # similitudes-ex4, question 2
    f = lambda q, aa=aa: aa * (q - al) + be
    ga = f(be)
    assert sp.simplify(sp.expand_complex(f(ga) - al)) == 0  # ex4, question 3 : f(γ) = α
    assert egal(sp.Abs(aa), 1)  # ex4 : |a| = 1, triangle équilatéral
assert sp.factor(a**3 - 1) == (a - 1) * (a**2 + a + 1)  # ex4, question 2
for nom, A, B, bonne, leurres in [
    ("ex4", 0, 2, [1 + I * sqrt(3), 1 - I * sqrt(3)], [[1 + I * sqrt(3)], [-1 + I * sqrt(3), -1 - I * sqrt(3)], []]),
    ("ex4 v1", 0, 2 * I, [-sqrt(3) + I, sqrt(3) + I], [[-sqrt(3) + I], [-sqrt(3) - I, sqrt(3) - I], []]),
    ("ex4 v2", -1, 1, [I * sqrt(3), -I * sqrt(3)], [[I * sqrt(3)], [sqrt(3), -sqrt(3)], []]),
]:
    gam = [B + aa * (B - A) for aa in [sp.exp(2 * I * pi / 3), sp.exp(-2 * I * pi / 3)]]
    assert meme_ensemble(gam, bonne), nom  # similitudes-ex4 : réponse
    for G in bonne:
        assert egal(sp.Abs(G - A), sp.Abs(B - A)) and egal(sp.Abs(G - B), sp.Abs(B - A)), nom  # triangle équilatéral
    for L in leurres:
        assert not meme_ensemble(gam, L), nom  # leurre faux (dont « aucune »)

# ======================================================================================
# Carte geo-l1-complexes-euler
# ======================================================================================
C = "geo-l1-complexes-euler"
E = sp.exp(I * th)
assert egal((E + 1 / E) / 2, sp.cos(th)) and egal((E - 1 / E) / (2 * I), sp.sin(th))  # carte : formules d'Euler
A_, B_ = sp.symbols("A B", real=True)
assert egal(sp.exp(I * A_) + sp.exp(I * B_), 2 * sp.cos((A_ - B_) / 2) * sp.exp(I * (A_ + B_) / 2))  # carte : angle moitié (+)
assert egal(sp.exp(I * A_) - sp.exp(I * B_), 2 * I * sp.sin((A_ - B_) / 2) * sp.exp(I * (A_ + B_) / 2))  # carte : angle moitié (-)
# Linéarisation : pour p + q ≤ 6, cos^p sin^q est combinaison des cos(kθ), sin(kθ), k ≤ p+q, cos seuls si q pair, sin seuls si q impair.
for p_ in range(5):
    for q_ in range(5):
        g = sp.cos(th) ** p_ * sp.sin(th) ** q_
        n = p_ + q_
        # Coefficients de Fourier exacts : g = Σ α_k cos kθ + β_k sin kθ, k ≤ n.
        # Développement en exponentielles : c_m = coefficient de e^{imθ} (u = e^{iθ}).
        uu = sp.symbols("u")
        dev = sp.expand(((uu + 1 / uu) / 2) ** p_ * ((uu - 1 / uu) / (2 * I)) ** q_)
        c = {m: dev.coeff(uu, m) for m in range(-n - 1, n + 2)}
        assert c[n + 1] == 0 and c[-n - 1] == 0  # carte : pas d'exposant au-delà de p+q
        alpha = [c[0]] + [sp.simplify(c[k] + c[-k]) for k in range(1, n + 1)]
        beta = [sp.simplify(I * (c[k] - c[-k])) for k in range(1, n + 1)]
        rec = sum(alpha[k] * sp.cos(k * th) for k in range(n + 1)) + sum(beta[k - 1] * sp.sin(k * th) for k in range(1, n + 1))
        assert identite_trigo(g, rec)  # carte : linéarisation avec 0 ≤ k ≤ p+q
        if q_ % 2 == 0:
            assert all(bb == 0 for bb in beta)  # carte : q pair → cosinus seuls
        else:
            assert all(aa == 0 for aa in alpha)  # carte : q impair → sinus seuls
# Exemple : cos^3 = (cos 3θ + 3 cos θ)/4 et ∫_0^{π/2} cos^3 = 2/3 ; 1 + e^{iθ} = 2 cos(θ/2) e^{iθ/2}.
assert identite_trigo(sp.cos(th) ** 3, (sp.cos(3 * th) + 3 * sp.cos(th)) / 4)  # carte (exemple)
assert sp.integrate(sp.cos(th) ** 3, (th, 0, pi / 2)) == sp.Rational(2, 3)  # carte (exemple)
assert egal(1 + E, 2 * sp.cos(th / 2) * sp.exp(I * th / 2))  # carte (exemple)
# Contre-exemples.
assert egal(sp.exp(4 * I * pi / 3) + 1, sp.exp(-I * pi / 3)) and egal(2 * sp.cos(2 * pi / 3) * sp.exp(2 * I * pi / 3), sp.exp(-I * pi / 3))  # contre-exemple (module)
assert egal(sp.Abs(sp.exp(4 * I * pi / 3) + 1), 1) and sp.cos(2 * pi / 3) < 0  # contre-exemple (module)
assert identite_trigo(sp.cos(th) * sp.sin(th), sp.sin(2 * th) / 2)  # contre-exemple (indices jusqu'à p+q)
al_, be_, ga_ = sp.symbols("al be ga")
eqs = [sp.Eq(sp.sin(2 * v), al_ + be_ * sp.cos(v) + ga_ * sp.sin(v)) for v in [0, pi, pi / 2, pi / 4]]
assert sp.solve(eqs, [al_, be_, ga_], dict=True) == []  # contre-exemple : sin 2θ n'est pas combinaison de 1, cos, sin
# Trous : leurres faux (θ = π/2 pour cos, θ = π/2 pour sin).
for L in [(E - 1 / E) / 2, E + 1 / E, (E + 1 / E) / (2 * I)]:
    assert jamais_identique(L, sp.cos(th))  # trou 1 : leurre faux
for L in [(E - 1 / E) / 2, (1 / E - E) / (2 * I), (E + 1 / E) / (2 * I)]:
    assert jamais_identique(L, sp.sin(th))  # trou 2 : leurre faux
# Trou 3 : leurres max(p,q) (p = q = 1, voir contre-exemple), pq (p = 1, q = 0 : cos θ n'est pas constant), p (p = 0, q = 1).
assert jamais_identique(sp.cos(th), sp.cos(0 * th) * sp.cos(th).subs(th, 0))  # trou 3, leurre « k ≤ pq »
assert jamais_identique(sp.sin(th), sp.sin(th).subs(th, 0))  # trou 3, leurre « k ≤ p »
# Affirmations.
assert egal(sp.Abs(sp.exp(2 * I * pi) + 1), 2) and 2 * sp.cos(pi) == -2  # affirmation 1
assert identite_trigo(sp.cos(th) ** 4, (sp.cos(4 * th) + 4 * sp.cos(2 * th) + 3) / 8)  # affirmation 2
assert sp.sin(0) ** 2 == 0 and (1 + sp.cos(0)) / 2 == 1  # affirmation 3
assert sp.sin(pi / 2) ** 3 == 1 and (3 * sp.sin(pi / 2) + sp.sin(3 * pi / 2)) / 4 == sp.Rational(1, 2)  # affirmation 4

# --- ex1 et variantes : linéarisations.
for nom, f, bonne, leurres in [
    ("ex1", sp.cos(th) ** 2 * sp.sin(th), (sp.sin(th) + sp.sin(3 * th)) / 4,
     [(3 * sp.sin(th) - sp.sin(3 * th)) / 4, (sp.sin(th) - sp.sin(3 * th)) / 4, (sp.sin(th) + sp.sin(3 * th)) / 2]),
    ("ex1 v1", sp.sin(th) ** 2 * sp.cos(th), (sp.cos(th) - sp.cos(3 * th)) / 4,
     [(sp.cos(th) + sp.cos(3 * th)) / 4, (3 * sp.cos(th) + sp.cos(3 * th)) / 4, (sp.cos(th) - sp.cos(3 * th)) / 2]),
    ("ex1 v2", sp.sin(th) ** 4, (sp.cos(4 * th) - 4 * sp.cos(2 * th) + 3) / 8,
     [(sp.cos(4 * th) + 4 * sp.cos(2 * th) + 3) / 8, (sp.cos(4 * th) - 4 * sp.cos(2 * th) + 3) / 16, (3 - 4 * sp.cos(2 * th) - sp.cos(4 * th)) / 8]),
]:
    assert identite_trigo(f, bonne), nom  # euler-ex1 : réponse
    for L in leurres:
        assert jamais_identique(f, L), nom  # leurre faux

# --- ex2 et variantes : équations cos px + cos qx = 0 sur [0, π].
for nom, f, bonne, leurres in [
    ("ex2", sp.cos(th) + sp.cos(3 * th), [pi / 4, pi / 2, 3 * pi / 4], [[pi / 4, 3 * pi / 4], [pi / 2], [0, pi / 4, pi / 2, 3 * pi / 4, pi]]),
    ("ex2 v1", sp.cos(th) + sp.cos(5 * th), [pi / 6, pi / 4, pi / 2, 3 * pi / 4, 5 * pi / 6],
     [[pi / 6, pi / 2, 5 * pi / 6], [pi / 4, 3 * pi / 4], [pi / 6, pi / 4, 3 * pi / 4, 5 * pi / 6]]),
    ("ex2 v2", sp.cos(2 * th) + sp.cos(4 * th), [pi / 6, pi / 2, 5 * pi / 6], [[pi / 6, 5 * pi / 6], [pi / 2], [pi / 3, pi / 2, 2 * pi / 3]]),
]:
    S = sp.solveset(f, th, sp.Interval(0, pi))
    assert S == sp.FiniteSet(*bonne), nom  # euler-ex2 : réponse
    for L in leurres:
        assert S != sp.FiniteSet(*L), nom  # leurre faux
assert identite_trigo(sp.cos(th) + sp.cos(3 * th), 2 * sp.cos(th) * sp.cos(2 * th))  # ex2 : factorisation citée
assert identite_trigo(sp.cos(th) + sp.cos(5 * th), 2 * sp.cos(2 * th) * sp.cos(3 * th))  # ex2 v1
assert identite_trigo(sp.cos(2 * th) + sp.cos(4 * th), 2 * sp.cos(th) * sp.cos(3 * th))  # ex2 v2

# --- ex3 et variantes : intégrales.
for nom, f, a0, b0, bonne, leurres in [
    ("ex3", sp.sin(th) ** 2 * sp.cos(th) ** 2, 0, pi / 2, pi / 16, [pi / 8, pi / 32, sp.Rational(1, 8)]),
    ("ex3 v1", sp.cos(th) ** 4, 0, pi, 3 * pi / 8, [3 * pi / 16, pi / 2, 3 * pi / 4]),
    ("ex3 v2", sp.sin(th) ** 3, 0, pi / 3, sp.Rational(5, 24), [sp.Rational(5, 12), sp.Rational(7, 24), sp.Rational(1, 8)]),
]:
    val = sp.integrate(f, (th, a0, b0))
    assert egal(val, bonne), nom  # euler-ex3 : réponse
    for L in leurres:
        assert not egal(val, L), nom  # leurre faux
assert identite_trigo(sp.sin(th) ** 2 * sp.cos(th) ** 2, (1 - sp.cos(4 * th)) / 8)  # ex3 : linéarisation citée
assert identite_trigo(sp.sin(th) ** 3, (3 * sp.sin(th) - sp.sin(3 * th)) / 4)  # ex3 v2 : linéarisation citée

# --- ex4 et variantes : sommes trigonométriques.
n_ = sp.symbols("n", integer=True, nonnegative=True)
for nom, terme, debut, formule, leurres in [
    ("ex4", lambda k: sp.cos(k * th), 0, lambda n: sp.cos(n * th / 2) * sp.sin((n + 1) * th / 2) / sp.sin(th / 2),
     [lambda n: sp.sin(n * th / 2) * sp.sin((n + 1) * th / 2) / sp.sin(th / 2),
      lambda n: sp.cos((n + 1) * th / 2) * sp.sin(n * th / 2) / sp.sin(th / 2),
      lambda n: sp.cos(n * th / 2) * sp.sin((n + 1) * th / 2) / sp.cos(th / 2)]),
    ("ex4 v1", lambda k: sp.sin(k * th), 0, lambda n: sp.sin(n * th / 2) * sp.sin((n + 1) * th / 2) / sp.sin(th / 2),
     [lambda n: sp.cos(n * th / 2) * sp.sin((n + 1) * th / 2) / sp.sin(th / 2),
      lambda n: sp.sin((n + 1) * th / 2) * sp.sin((n + 2) * th / 2) / sp.sin(th / 2),
      lambda n: sp.sin(n * th / 2) * sp.sin((n + 1) * th / 2) / sp.cos(th / 2)]),
    ("ex4 v2", lambda k: sp.cos((2 * k - 1) * th), 1, lambda n: sp.sin(2 * n * th) / (2 * sp.sin(th)),
     [lambda n: sp.sin(2 * n * th) / sp.sin(th), lambda n: sp.sin(n * th) / (2 * sp.sin(th)),
      lambda n: sp.cos(n * th) * sp.sin(n * th) / sp.sin(2 * th)]),
]:
    for n in range(1, 7):
        somme = sum(terme(k) for k in range(debut, n + 1))
        assert identite_trigo(somme, formule(n)), nom  # euler-ex4 : réponse
    for L in leurres:
        assert any(jamais_identique(sum(terme(k) for k in range(debut, n + 1)), L(n)) for n in range(1, 4)), nom  # leurre faux
assert egal(1 - E, -2 * I * sp.sin(th / 2) * sp.exp(I * th / 2))  # ex4 : factorisation citée dans l'indice
assert identite_trigo(sp.cos(th / 2) * sp.sin(th) / sp.sin(th / 2), 1 + sp.cos(th))  # ex4 : contrôle n = 1

# --- ex5 et variantes : ∫_0^π cos^{2m} = π C(2m,m)/4^m.
m_ = sp.symbols("m", integer=True, positive=True)
for m in range(1, 7):
    g = sum(sp.binomial(2 * m, k) * 2 * sp.cos(2 * (m - k) * th) for k in range(m)) + sp.binomial(2 * m, m)
    assert identite_trigo(sp.cos(th) ** (2 * m), g / 4**m)  # euler-ex5 : question 1
    assert sp.integrate(sp.cos(th) ** (2 * m), (th, 0, pi)) == pi * sp.binomial(2 * m, m) / 4**m  # euler-ex5 : question 2
for nom, pp, bonne, leurres in [
    ("ex5", 6, 5 * pi / 16, [5 * pi / 32, 3 * pi / 8, 15 * pi / 16]),
    ("ex5 v1", 8, 35 * pi / 128, [35 * pi / 256, 5 * pi / 16, 35 * pi / 64]),
    ("ex5 v2", 10, 63 * pi / 256, [63 * pi / 512, 35 * pi / 128, 63 * pi / 128]),
]:
    val = sp.integrate(sp.cos(th) ** pp, (th, 0, pi))
    assert val == bonne, nom  # euler-ex5 : réponse
    for L in leurres:
        assert val != L, nom  # leurre faux
assert sp.binomial(6, 3) == 20 and sp.binomial(8, 4) == 70 and sp.binomial(10, 5) == 252  # ex5 : coefficients cités

print("geo_l1_complexes : toutes les vérifications passent.")
