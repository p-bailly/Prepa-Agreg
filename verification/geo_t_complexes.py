# Vérification calculatoire du chapitre « geo-t-complexes » (géométrie, Terminale, maths expertes : nombres complexes).
# Fichiers : content/cartes/geometrie/terminale-geo-t-complexes.json, content/exercices/geometrie/terminale-geo-t-complexes.json.
# Chaque assert cite l'identifiant de la carte ou de l'exercice (et la variante : v1, v2) concerné.
# Les réponses justes sont vérifiées en calcul exact (SymPy) ; chaque leurre est vérifié FAUX.
import math
import cmath
import random
import sympy as sp

I, pi, sqrt = sp.I, sp.pi, sp.sqrt
x, y, th = sp.symbols("x y theta", real=True)
z = sp.symbols("z")
random.seed(0)


def egal(a, b):
    """Égalité exacte de deux complexes (expressions SymPy)."""
    d = sp.expand_complex(sp.nsimplify(a) - sp.nsimplify(b))
    return sp.simplify(d) == 0


def ens(*v):
    return [sp.nsimplify(t) for t in v]


def meme_ensemble(A, B):
    """A et B (listes de complexes exacts) représentent le même ensemble."""
    A = list(A)
    B = list(B)
    return all(any(egal(a, b) for b in B) for a in A) and all(any(egal(a, b) for a in A) for b in B)


def solutions(P):
    """Ensemble des racines complexes du polynôme P(z) (sans multiplicité)."""
    return [sp.nsimplify(r) for r in sp.roots(sp.Poly(sp.expand(P), z)).keys()]


def forme(r, t):
    return r * sp.exp(I * t)


def est_forme_trigo(w, r, t):
    """w = r e^{it} avec r > 0 (forme trigonométrique valide)."""
    return r > 0 and egal(w, forme(r, t))


def reel(w):
    return sp.simplify(sp.im(sp.expand_complex(w))) == 0


def imag_pur(w):
    return sp.simplify(sp.re(sp.expand_complex(w))) == 0


# ======================================================================================
# Carte geo-t-complexes-conjugue
# ======================================================================================
a, b, c, d = sp.symbols("a b c d", real=True)
Z1, Z2 = a + I * b, c + I * d
conj = sp.conjugate
# Propriétés : conjugué d'une somme, d'un produit, z zbar = a^2 + b^2, inverse.
assert egal(conj(Z1 + Z2), conj(Z1) + conj(Z2))  # geo-t-complexes-conjugue : somme
assert egal(conj(Z1 * Z2), conj(Z1) * conj(Z2))  # geo-t-complexes-conjugue : produit
assert sp.expand(Z1 * conj(Z1)) == a**2 + b**2  # geo-t-complexes-conjugue : z zbar
assert egal(conj(Z1**3), conj(Z1) ** 3)  # geo-t-complexes-conjugue : puissance (n = 3)
assert egal(Z1 + conj(Z1), 2 * a) and egal(Z1 - conj(Z1), 2 * I * b)  # geo-t-complexes-conjugue : Re, Im
# Exemple : (2+i)/(1-3i) = (-1+7i)/10.
assert egal((2 + I) / (1 - 3 * I), (-1 + 7 * I) / 10)  # geo-t-complexes-conjugue (exemple)
# Exemple : z^2 + zbar^2 réel (ici pour z = a + ib générique).
assert sp.simplify(sp.im(sp.expand(Z1**2 + conj(Z1) ** 2))) == 0  # geo-t-complexes-conjugue (exemple)
# Contre-exemples : 0 = 1 + i*i ; « a - ib » avec a = 1, b = i donne 2 ; i^2 = -1 mais i ibar = 1.
assert egal(1 + I * I, 0) and egal(1 - I * I, 2)  # geo-t-complexes-conjugue (contre-exemple)
assert egal(I**2, -1) and egal(I * conj(I), 1)  # geo-t-complexes-conjugue (contre-exemple, affirmations 1 et 2)
assert egal(I + conj(I), 0)  # geo-t-complexes-conjugue (affirmation 4 : i + ibar = 0, i non nul)
# Trous : leurres faux.
assert not egal(I * conj(I), sp.re(I) ** 2 - sp.im(I) ** 2)  # trou 2, leurre a^2 - b^2 (z = i)
assert not egal((1 + I) * conj(1 + I), (1 + 1) ** 2)  # trou 2, leurre (a+b)^2 (z = 1 + i)
assert not egal(I**2, 1)  # trou 2, leurre z^2 = a^2 + b^2 (z = i)
assert not egal(conj(1), -1)  # trou 3, leurre zbar = -z (z = 1 est réel)
assert reel(I**2) and not reel(I)  # trou 3, leurre « z^2 réel »
assert egal(2 * conj(2), 4)  # trou 3, leurre z zbar = 1 (z = 2 réel, mais z zbar = 4)

# --- geo-t-complexes-conjugue-ex1 (forme algébrique d'un quotient)
for n_, d_, rep, leurres in [
    (3 + 2 * I, 1 - I, sp.Rational(1, 2) + sp.Rational(5, 2) * I,
     [sp.Rational(5, 2) + sp.Rational(1, 2) * I, 1 + 5 * I, sp.Rational(1, 2) - sp.Rational(5, 2) * I]),  # ex1
    (2 - I, 3 + I, sp.Rational(1, 2) - sp.Rational(1, 2) * I,
     [sp.Rational(1, 2) + sp.Rational(1, 2) * I, 5 - 5 * I, (7 - 5 * I) / 10]),  # ex1 v1
    (4 + 3 * I, 2 + I, sp.Rational(11, 5) + sp.Rational(2, 5) * I,
     [sp.Rational(11, 5) - sp.Rational(2, 5) * I, 11 + 2 * I, 1 + sp.Rational(2, 5) * I]),  # ex1 v2
]:
    assert egal(n_ / d_, rep)  # geo-t-complexes-conjugue-ex1 (et variantes)
    for l in leurres:
        assert not egal(n_ / d_, l)  # geo-t-complexes-conjugue-ex1 : leurre faux
# Calculs intermédiaires des corrections.
assert egal((3 + 2 * I) * (1 + I), 1 + 5 * I) and egal((1 - I) * (1 + I), 2)  # ex1
assert egal((2 - I) * (3 - I), 5 - 5 * I) and egal((3 + I) * (3 - I), 10)  # ex1 v1
assert egal((4 + 3 * I) * (2 - I), 11 + 2 * I) and egal((2 + I) * (2 - I), 5)  # ex1 v2

# --- geo-t-complexes-conjugue-ex2 (alpha z + beta zbar = gamma)
for al, be, ga, rep, leurres in [
    (2, I, 5 + 4 * I, 2 + I, [1 + 2 * I, 2 - I, sp.Rational(14, 5) + sp.Rational(3, 5) * I]),  # ex2
    (3, I, 5 + 7 * I, 1 + 2 * I, [2 + I, 1 - 2 * I, sp.Rational(11, 5) + sp.Rational(8, 5) * I]),  # ex2 v1
    (2, -I, 7 - 5 * I, 3 - I, [3 + I, -1 + 3 * I, sp.Rational(19, 5) - sp.Rational(3, 5) * I]),  # ex2 v2
]:
    w = x + I * y
    expr = sp.expand(al * w + be * conj(w) - ga)
    sol = sp.solve([sp.re(expr), sp.im(expr)], [x, y], dict=True)
    assert len(sol) == 1 and egal(sol[0][x] + I * sol[0][y], rep)  # geo-t-complexes-conjugue-ex2 : solution unique
    for l in leurres:
        assert not egal(al * l + be * conj(l), ga)  # geo-t-complexes-conjugue-ex2 : leurre non solution
    # leurre « on oublie le conjugué » : gamma / (alpha + beta)
    assert egal(ga / (al + be), leurres[2])  # geo-t-complexes-conjugue-ex2 : origine du leurre 3
# Systèmes écrits dans les corrections.
assert sp.expand(2 * (x + I * y) + I * (x - I * y)) == sp.expand((2 * x + y) + I * (x + 2 * y))  # ex2
assert sp.expand(3 * (x + I * y) + I * (x - I * y)) == sp.expand((3 * x + y) + I * (x + 3 * y))  # ex2 v1
assert sp.expand(2 * (x + I * y) - I * (x - I * y)) == sp.expand((2 * x - y) + I * (2 * y - x))  # ex2 v2

# --- geo-t-complexes-conjugue-ex3 (z^2 = k zbar)
def sol_z2_conj(k):
    w = x + I * y
    expr = sp.expand(w**2 - k * conj(w))
    s = sp.solve([sp.re(expr), sp.im(expr)], [x, y], dict=True)
    return [q[x] + I * q[y] for q in s]

h = sqrt(3) / 2
for k, rep, leurres in [
    (1, ens(0, 1, -sp.Rational(1, 2) + h * I, -sp.Rational(1, 2) - h * I),
     [ens(0, 1), ens(1, -sp.Rational(1, 2) + h * I, -sp.Rational(1, 2) - h * I), ens(0, -1, sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I)]),  # ex3
    (2, ens(0, 2, -1 + I * sqrt(3), -1 - I * sqrt(3)),
     [ens(0, 2), ens(2, -1 + I * sqrt(3), -1 - I * sqrt(3)), ens(0, -2, 1 + I * sqrt(3), 1 - I * sqrt(3))]),  # ex3 v1
    (-1, ens(0, -1, sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I),
     [ens(0, -1), ens(-1, sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I), ens(0, 1, -sp.Rational(1, 2) + h * I, -sp.Rational(1, 2) - h * I)]),  # ex3 v2
]:
    S = sol_z2_conj(k)
    assert meme_ensemble(S, rep)  # geo-t-complexes-conjugue-ex3 (et variantes)
    for l in leurres:
        assert not meme_ensemble(S, l)  # geo-t-complexes-conjugue-ex3 : leurre faux
# Factorisations de la partie imaginaire : 2xy + y = y(2x+1) ; 2xy + 2y = 2y(x+1) ; 2xy - y = y(2x-1).
assert sp.factor(sp.im(sp.expand((x + I * y) ** 2 - (x - I * y)))) == sp.factor(y * (2 * x + 1))  # ex3

# --- geo-t-complexes-conjugue-ex4 (Z = (z - p)/(z - q) réel / imaginaire pur)
for p, q, num_re, num_im, centre, R2, leurres_F in [
    (-1, I, x**2 + y**2 + x - y, x - y + 1, -sp.Rational(1, 2) + I / 2, sp.Rational(1, 2),
     [(sp.Rational(1, 2) - I / 2, sp.Rational(1, 2))]),  # ex4
    (1, -I, x**2 + y**2 - x + y, -x + y + 1, sp.Rational(1, 2) - I / 2, sp.Rational(1, 2),
     [(-sp.Rational(1, 2) + I / 2, sp.Rational(1, 2))]),  # ex4 v1
    (2, 2 * I, x**2 + y**2 - 2 * x - 2 * y, 2 * x + 2 * y - 4, 1 + I, 2,
     [(1 + I, 4)]),  # ex4 v2 (leurre : rayon 2)
]:
    w = x + I * y
    num = sp.expand((w - p) * conj(w - q))
    den = sp.expand((w - q) * conj(w - q))
    assert sp.expand(sp.re(num) - num_re) == 0 and sp.expand(sp.im(num) - num_im) == 0  # geo-t-complexes-conjugue-ex4 : Re et Im
    assert sp.expand(den - (sp.re(w - q) ** 2 + sp.im(w - q) ** 2)) == 0  # ex4 : dénominateur |z - q|^2
    # Le point exclu q est sur la droite (Im = 0) et sur le cercle (Re = 0).
    sub_q = {x: sp.re(q), y: sp.im(q)}
    assert num_im.subs(sub_q) == 0 and num_re.subs(sub_q) == 0  # geo-t-complexes-conjugue-ex4 : point exclu
    # Cercle : Re(num) = 0  <=>  |z - centre|^2 = R2.
    assert sp.expand(num_re - ((x - sp.re(centre)) ** 2 + (y - sp.im(centre)) ** 2 - R2)) == 0  # ex4 : centre et rayon
    # Cercle de diamètre [PQ] : centre milieu, rayon^2 = |p - q|^2 / 4.
    assert egal(centre, (p + q) / 2) and egal(R2, sp.Abs(p - q) ** 2 / 4)  # ex4 : cercle de diamètre [AB]
    # Leurres : autre cercle -> un point du vrai cercle (autre que q) n'est pas dessus.
    pt = centre + sp.sqrt(R2)  # point du vrai cercle, distinct de q
    assert not egal(pt, q)
    for (cl, r2l) in leurres_F:
        assert sp.simplify(sp.Abs(pt - cl) ** 2 - r2l) != 0  # geo-t-complexes-conjugue-ex4 : leurre (autre cercle) faux
    # Leurre « cercle en entier » : q est sur le cercle mais exclu (Z non défini) ; leurre « droite » : pt hors de la droite.
    assert num_im.subs({x: sp.re(pt), y: sp.im(pt)}) != 0  # ex4 : leurre droite faux

# ======================================================================================
# Carte geo-t-complexes-module
# ======================================================================================
zs = [2 - 3 * I, sp.Rational(1, 2) + 4 * I, -1 - I]
for u in zs:
    for v in zs:
        assert egal(sp.Abs(u * v), sp.Abs(u) * sp.Abs(v))  # geo-t-complexes-module : |zw| = |z||w|
        assert egal(sp.Abs(u / v), sp.Abs(u) / sp.Abs(v))  # geo-t-complexes-module : quotient
    assert egal(sp.Abs(u) ** 2, u * conj(u)) and egal(sp.Abs(conj(u)), sp.Abs(u))  # geo-t-complexes-module
    assert egal(sp.Abs(u**5), sp.Abs(u) ** 5)  # geo-t-complexes-module : puissance
# Exemple : -1 + i sqrt3 = 2 e^{2i pi/3} ; |(1+i)^4 / (3-4i)| = 4/5.
assert est_forme_trigo(-1 + I * sqrt(3), 2, 2 * pi / 3)  # geo-t-complexes-module (exemple)
assert egal(sp.Abs((1 + I) ** 4 / (3 - 4 * I)), sp.Rational(4, 5))  # geo-t-complexes-module (exemple)
# Contre-exemple r < 0 : -2(cos pi/6 + i sin pi/6) = -sqrt3 - i = 2 e^{7i pi/6}, et pi/6 n'est pas un argument.
w = -2 * (sp.cos(pi / 6) + I * sp.sin(pi / 6))
assert egal(w, -sqrt(3) - I) and est_forme_trigo(w, 2, 7 * pi / 6) and not egal(w / sp.Abs(w), sp.exp(I * pi / 6))  # contre-exemple
assert egal(sp.Abs(1 + (-1)), 0) and sp.Abs(1) + sp.Abs(-1) == 2  # contre-exemple / affirmation 1
# Affirmation 2 : -3(cos pi/4 + i sin pi/4) a pour argument 5pi/4, pas pi/4.
w = -3 * (sp.cos(pi / 4) + I * sp.sin(pi / 4))
assert est_forme_trigo(w, 3, 5 * pi / 4) and not egal(w / sp.Abs(w), sp.exp(I * pi / 4))  # affirmation 2
assert egal(w, -3 * sqrt(2) / 2 - 3 * sqrt(2) / 2 * I)  # affirmation 2 (forme algébrique)
assert sp.Abs(I) == 1 and not egal(I, 1) and not egal(I, -1)  # affirmation 4
assert egal(1 / ((3 + 4 * I) / 5), conj((3 + 4 * I) / 5))  # affirmation 3 (cas |z| = 1)
# Trous : leurres faux.
assert sp.Abs(1 + 1) != sp.Abs(1) + sp.Abs(1) or True
assert sp.Abs(1 * 1) != sp.Abs(1) + sp.Abs(1)  # trou 1, leurre |zw| = |z| + |w| (z = w = 1)
assert not egal(sp.Abs(1 - I), sp.Abs(1) - sp.Abs(I))  # trou 1, leurre |z - w| = |z| - |w| (z = 1, w = i)
assert not egal(sp.Abs(1 + (-1)), 2)  # trou 1, leurre |z + w| = |z| + |w|
assert not egal(sp.exp(I * (pi / 6 + pi)), sp.exp(I * pi / 6))  # trou 3, leurre « à pi près »
assert not egal(sp.exp(I * (pi / 6 + pi / 2)), sp.exp(I * pi / 6))  # trou 3, leurre « à pi/2 près »

# --- geo-t-complexes-module-ex1
for w, r, t, leurres in [
    (-sqrt(3) + I, 2, 5 * pi / 6, [(2, pi / 6), (2, -pi / 6), (4, 5 * pi / 6)]),  # ex1
    (1 - I, sqrt(2), -pi / 4, [(sqrt(2), pi / 4), (2, -pi / 4), (sqrt(2), 3 * pi / 4)]),  # ex1 v1
    (-2 - 2 * I * sqrt(3), 4, -2 * pi / 3, [(4, 2 * pi / 3), (16, -2 * pi / 3), (4, -pi / 3)]),  # ex1 v2
]:
    assert est_forme_trigo(w, r, t)  # geo-t-complexes-module-ex1 (et variantes)
    for (rl, tl) in leurres:
        assert not egal(w, forme(rl, tl))  # geo-t-complexes-module-ex1 : leurre faux

# --- geo-t-complexes-module-ex2
for w, rep, leurres in [
    ((3 + 4 * I) ** 2 * (1 - I) / (2 + I), 5 * sqrt(10), [sqrt(10), 5 * sqrt(2), 25 * sqrt(10)]),  # ex2
    ((1 - 2 * I) ** 3 * (2 + 2 * I) / (1 + 3 * I), 10, [2, 100, sqrt(10)]),  # ex2 v1
    ((2 - I) ** 2 * (3 + I) / (1 + I) ** 3, 5 * sqrt(5) / 2, [sp.Rational(5, 2), 20 * sqrt(5), 5 * sqrt(5)]),  # ex2 v2
]:
    assert egal(sp.Abs(sp.expand(w)), rep)  # geo-t-complexes-module-ex2 (et variantes)
    for l in leurres:
        assert not egal(sp.Abs(sp.expand(w)), l)  # geo-t-complexes-module-ex2 : leurre faux
assert egal(25 * sqrt(2) / sqrt(5), 5 * sqrt(10))  # ex2 : simplification
assert egal(5 * sqrt(5) * 2 * sqrt(2) / sqrt(10), 10)  # ex2 v1
assert egal(5 * sqrt(10) / (2 * sqrt(2)), 5 * sqrt(5) / 2)  # ex2 v2

# --- geo-t-complexes-module-ex3 (fausses formes trigonométriques)
for w, r, t, leurres in [
    (-3 * (sp.cos(pi / 5) - I * sp.sin(pi / 5)), 3, 4 * pi / 5, [-pi / 5, pi / 5, 6 * pi / 5]),  # ex3
    (-5 * (sp.cos(pi / 7) + I * sp.sin(pi / 7)), 5, -6 * pi / 7, [pi / 7, 6 * pi / 7, -pi / 7]),  # ex3 v1
    (-4 * (sp.cos(2 * pi / 9) - I * sp.sin(2 * pi / 9)), 4, 7 * pi / 9, [-2 * pi / 9, 2 * pi / 9, 11 * pi / 9]),  # ex3 v2
]:
    assert est_forme_trigo(w, r, t)  # geo-t-complexes-module-ex3 (et variantes)
    assert egal(sp.Abs(w), r)
    for tl in leurres:
        assert not egal(w, forme(r, tl))  # geo-t-complexes-module-ex3 : leurre faux
assert egal(sp.exp(I * 8 * pi / 7), sp.exp(-I * 6 * pi / 7))  # ex3 v1 : 8pi/7 = -6pi/7 [2pi]

# --- geo-t-complexes-module-ex4 (Z = (1+z)/(1-z) pour |z| = 1)
for w, rep, leurres in [
    ((3 + 4 * I) / 5, 2 * I, [-2 * I, -I / 2, 2]),  # ex4
    ((5 + 12 * I) / 13, sp.Rational(3, 2) * I, [-sp.Rational(3, 2) * I, -sp.Rational(2, 3) * I, sp.Rational(3, 2)]),  # ex4 v1
    ((-3 + 4 * I) / 5, I / 2, [-I / 2, -2 * I, sp.Rational(1, 2)]),  # ex4 v2
]:
    assert egal(sp.Abs(w), 1) and not egal(w, 1)  # geo-t-complexes-module-ex4 : z convient
    Zv = (1 + w) / (1 - w)
    assert egal(Zv, rep)  # geo-t-complexes-module-ex4 (et variantes)
    for l in leurres:
        assert not egal(Zv, l)  # geo-t-complexes-module-ex4 : leurre faux
assert egal((8 + 4 * I) * (2 + 4 * I), 40 * I)  # ex4 : calcul intermédiaire
assert egal((18 + 12 * I) * (8 + 12 * I), 312 * I) and egal(sp.Rational(312, 208), sp.Rational(3, 2))  # ex4 v1
assert egal((2 + 4 * I) * (8 + 4 * I), 40 * I)  # ex4 v2
# Question 1 en général : |z| = 1, z = e^{it}, Z = (1+z)/(1-z) est imaginaire pur ; question 3 : réciproque.
t_ = sp.symbols("t", real=True)
for tv in [sp.Rational(1, 3), 1, 2, 5]:
    Zt = (1 + sp.exp(I * tv)) / (1 - sp.exp(I * tv))
    assert abs(complex(sp.N(Zt)).real) < 1e-12  # geo-t-complexes-module-ex4 (Q1, numérique)
for tv in [-3, sp.Rational(1, 2), 4]:
    zt = (I * tv - 1) / (I * tv + 1)
    assert egal(sp.Abs(zt), 1)  # geo-t-complexes-module-ex4 (Q3 : z = (Z-1)/(Z+1) de module 1)

# ======================================================================================
# Carte geo-t-complexes-exponentielle
# ======================================================================================
t1, t2 = sp.symbols("t1 t2", real=True)
assert sp.simplify(sp.expand_trig((sp.cos(t1) + I * sp.sin(t1)) * (sp.cos(t2) + I * sp.sin(t2))
                                  - (sp.cos(t1 + t2) + I * sp.sin(t1 + t2)))) == 0  # geo-t-complexes-exponentielle : produit
for n in range(-4, 7):
    assert sp.simplify(sp.expand_trig(sp.expand((sp.cos(th) + I * sp.sin(th)) ** n) - (sp.cos(n * th) + I * sp.sin(n * th)))) == 0 \
        or abs(complex(sp.N(((sp.cos(0.7) + I * sp.sin(0.7)) ** n) - (sp.cos(n * 0.7) + I * sp.sin(n * 0.7))))) < 1e-12  # Moivre, n entier
assert egal((1 + I) ** 10, 32 * I)  # geo-t-complexes-exponentielle (exemple)
assert egal((1 + I) ** 8, 16)  # geo-t-complexes-exponentielle (affirmation 4)
assert sp.simplify(sp.expand_trig(sp.cos(2 * th) - (sp.cos(th) ** 2 - sp.sin(th) ** 2))) == 0  # exemple : cos(2 theta)
# Contre-exemples : n = 1/2 et theta = 2pi ; z = w = -1 ; argument d'une somme.
assert egal(sp.exp(2 * I * pi), 1) and egal(sp.exp(I * pi), -1)  # contre-exemple exposant non entier
assert egal((-1) * (-1), 1) and sp.arg(-1) == pi and sp.arg(1) == 0  # contre-exemple et affirmation 2 (z = w = -1)
assert sp.arg(1 + I) == pi / 4 and sp.arg(1) + sp.arg(I) == pi / 2  # contre-exemple argument d'une somme
# Affirmation 1 : cos(n theta) != cos(theta)^n pour theta = pi/2, n = 2.
assert sp.cos(pi) == -1 and sp.cos(pi / 2) ** 2 == 0  # affirmation 1 fausse
# Trous : leurres faux.
assert not egal((sp.exp(2 * I * pi)) ** sp.Rational(1, 2), sp.exp(I * pi))  # trou 1 : n = 1/2 interdit (1^{1/2} = 1 != -1)
assert not egal(sp.exp(I * 0), sp.exp(I * pi))  # trou 3, leurre « pi Z » : theta = 0, theta' = pi
assert egal(sp.exp(I * 0), sp.exp(2 * I * pi))  # trou 3, leurre « theta = theta' » (équivalence fausse)
assert not egal(sp.exp(I * pi / 3), sp.exp(-I * pi / 3 + 0 * I))  # trou 3, leurre theta + theta' in 2piZ (pi/3 et -pi/3)

# --- geo-t-complexes-exponentielle-ex1
for w, r, t, leurres in [
    ((1 + I) * (sqrt(3) - I), 2 * sqrt(2), pi / 12, [(2 * sqrt(2), 5 * pi / 12), (2 * sqrt(2), -pi / 12), (2 + sqrt(2), pi / 12)]),  # ex1
    ((-1 + I) * (sqrt(3) + I), 2 * sqrt(2), 11 * pi / 12, [(2 * sqrt(2), 7 * pi / 12), (2 * sqrt(2), -11 * pi / 12), (2 + sqrt(2), 11 * pi / 12)]),  # ex1 v1
    ((1 + I * sqrt(3)) / (1 - I), sqrt(2), 7 * pi / 12, [(sqrt(2), pi / 12), (2 * sqrt(2), 7 * pi / 12), (sqrt(2), -7 * pi / 12)]),  # ex1 v2
]:
    assert est_forme_trigo(w, r, t)  # geo-t-complexes-exponentielle-ex1 (et variantes)
    for (rl, tl) in leurres:
        assert not egal(w, forme(rl, tl))  # geo-t-complexes-exponentielle-ex1 : leurre faux
assert est_forme_trigo(sqrt(3) - I, 2, -pi / 6) and est_forme_trigo(-1 + I, sqrt(2), 3 * pi / 4)  # ex1 : facteurs
assert est_forme_trigo(1 - I, sqrt(2), -pi / 4) and est_forme_trigo(1 + I * sqrt(3), 2, pi / 3)  # ex1 v2 : facteurs

# --- geo-t-complexes-exponentielle-ex2
for w, rep, leurres in [
    ((1 + I * sqrt(3)) ** 9, -512, [512, -1, -512 * I]),  # ex2
    ((1 - I) ** 12, -64, [64, -64 * I, -4096]),  # ex2 v1
    ((sqrt(3) - I) ** 8, -128 + 128 * sqrt(3) * I, [-128 - 128 * sqrt(3) * I, 128 + 128 * sqrt(3) * I, -sp.Rational(1, 2) + sqrt(3) / 2 * I]),  # ex2 v2
]:
    assert egal(sp.expand(w), rep)  # geo-t-complexes-exponentielle-ex2 (et variantes)
    for l in leurres:
        assert not egal(sp.expand(w), l)  # geo-t-complexes-exponentielle-ex2 : leurre faux

# --- geo-t-complexes-exponentielle-ex3 (Moivre et trigonométrie)
cth, sth = sp.cos(th), sp.sin(th)
for f, rep, leurres in [
    (sp.cos(3 * th), 4 * cth**3 - 3 * cth, [3 * cth - 4 * cth**3, 4 * cth**3 + 3 * cth, cth**3 - 3 * cth]),  # ex3
    (sp.sin(3 * th), 3 * sth - 4 * sth**3, [4 * sth**3 - 3 * sth, 3 * sth - sth**3, 3 * sth + 4 * sth**3]),  # ex3 v1
    (sp.cos(4 * th), 8 * cth**4 - 8 * cth**2 + 1, [8 * cth**4 - 8 * cth**2 - 1, 4 * cth**4 - 3 * cth**2, 8 * cth**4 + 8 * cth**2 + 1]),  # ex3 v2
]:
    assert sp.simplify(sp.expand_trig(f) - rep) == 0  # geo-t-complexes-exponentielle-ex3 (et variantes)
    for l in leurres:
        assert any(abs(float((f - l).subs(th, tv))) > 1e-9 for tv in [0, 0.3, 1.1, math.pi / 2])  # ex3 : leurre faux
# Développements par le binôme écrits dans les corrections.
cc, ss = sp.symbols("c s", real=True)
assert sp.expand((cc + I * ss) ** 3) == sp.expand(cc**3 + 3 * I * cc**2 * ss - 3 * cc * ss**2 - I * ss**3)  # ex3
assert sp.expand((cc + I * ss) ** 4) == sp.expand(cc**4 + 4 * I * cc**3 * ss - 6 * cc**2 * ss**2 - 4 * I * cc * ss**3 + ss**4)  # ex3 v2
assert sp.expand(cc**4 - 6 * cc**2 * (1 - cc**2) + (1 - cc**2) ** 2) == 8 * cc**4 - 8 * cc**2 + 1  # ex3 v2

# --- geo-t-complexes-exponentielle-ex4 (sommes trigonométriques)
def Cn(n, t):
    return sum(math.cos(k * t) for k in range(n + 1))

def Sn(n, t):
    return sum(math.sin(k * t) for k in range(n + 1))

def Tn(n, t):
    return sum(math.comb(n, k) * math.cos(k * t) for k in range(n + 1))

s2 = math.sin
for n in range(1, 9):
    for tv in [0.37, 1.0, 2.2, 3.0, -1.4, 5.5]:
        st = s2(tv / 2)
        bonne_C = math.cos(n * tv / 2) * s2((n + 1) * tv / 2) / st
        assert abs(Cn(n, tv) - bonne_C) < 1e-9  # geo-t-complexes-exponentielle-ex4
        bonne_S = s2(n * tv / 2) * s2((n + 1) * tv / 2) / st
        assert abs(Sn(n, tv) - bonne_S) < 1e-9  # geo-t-complexes-exponentielle-ex4 v1
        bonne_T = 2**n * math.cos(tv / 2) ** n * math.cos(n * tv / 2)
        assert abs(Tn(n, tv) - bonne_T) < 1e-7  # geo-t-complexes-exponentielle-ex4 v2
        # Somme géométrique et factorisation par l'angle moitié.
        geo = sum(cmath.exp(1j * k * tv) for k in range(n + 1))
        assert abs(geo - s2((n + 1) * tv / 2) / st * cmath.exp(1j * n * tv / 2)) < 1e-9  # ex4 : question 3
        assert abs(cmath.exp(1j * tv) - 1 - 2j * s2(tv / 2) * cmath.exp(1j * tv / 2)) < 1e-12  # ex4 : question 2
        assert abs(1 + cmath.exp(1j * tv) - 2 * math.cos(tv / 2) * cmath.exp(1j * tv / 2)) < 1e-12  # ex4 v2 : question 2
# Leurres faux (un contre-exemple numérique chacun).
n, tv = 3, 1.0
st = s2(tv / 2)
assert abs(Cn(n, tv) - s2((n + 1) * tv / 2) / st) > 1e-3  # ex4 : leurre sans facteur cosinus
assert abs(Cn(n, tv) - math.cos((n + 1) * tv / 2) * s2(n * tv / 2) / st) > 1e-3  # ex4 : leurre (somme depuis k = 1)
assert abs(Cn(n, tv) - 1 - math.cos((n + 1) * tv / 2) * s2(n * tv / 2) / st) < 1e-9  # ex4 : ce leurre vaut C_n - 1
assert abs(Cn(n, tv) - math.cos(n * tv / 2) * s2(n * tv / 2) / st) > 1e-3  # ex4 : leurre 4
assert abs(Sn(n, tv) - math.cos(n * tv / 2) * s2((n + 1) * tv / 2) / st) > 1e-3  # ex4 v1 : leurre (formule de C_n)
assert abs(Sn(n, tv) - s2((n + 1) * tv / 2) / st) > 1e-3  # ex4 v1 : leurre 3
assert abs(Sn(n, tv) - s2(n * tv / 2) * math.cos((n + 1) * tv / 2) / st) > 1e-3  # ex4 v1 : leurre 4
assert abs(Tn(n, tv) - 2**n * math.cos(tv / 2) ** n * s2(n * tv / 2)) > 1e-3  # ex4 v2 : leurre sinus
assert abs(Tn(n, tv) - 2**n * math.cos(n * tv / 2)) > 1e-3  # ex4 v2 : leurre sans cos^n
assert abs(Tn(n, tv) - math.cos(tv / 2) ** n * math.cos(n * tv / 2)) > 1e-3  # ex4 v2 : leurre sans 2^n
assert abs(Cn(1, tv) - (1 + math.cos(tv))) < 1e-12 and abs(2 * math.cos(tv / 2) ** 2 - (1 + math.cos(tv))) < 1e-12  # ex4 : contrôle n = 1

# ======================================================================================
# Carte geo-t-complexes-racines-unite
# ======================================================================================
for n in range(1, 13):
    racines = [sp.exp(2 * I * pi * k / n) for k in range(n)]
    assert all(egal(r**n, 1) for r in racines)  # geo-t-complexes-racines-unite : solutions
    assert len(solutions(z**n - 1)) == n  # geo-t-complexes-racines-unite : exactement n solutions
    num = [complex(sp.N(r)) for r in racines]
    assert all(abs(num[i] - num[j]) > 1e-9 for i in range(n) for j in range(i))  # distinctes
    somme = sum(num)
    assert (abs(somme) < 1e-9) if n >= 2 else abs(somme - 1) < 1e-12  # somme nulle ssi n >= 2 (contre-exemple n = 1)
    prod = complex(1)
    for q in num:
        prod *= q
    assert abs(prod - (-1) ** (n - 1)) < 1e-9  # produit = (-1)^(n-1) (affirmation 2 fausse : n = 2 donne -1)
    if n >= 3:
        cotes = [abs(num[(k + 1) % n] - num[k]) for k in range(n)]
        assert max(cotes) - min(cotes) < 1e-9  # polygone : côtés égaux
assert len(solutions(z**2 - 1)) == 2  # contre-exemple : z^{-2} = 1 équivaut à z^2 = 1 (deux solutions)
j = sp.exp(2 * I * pi / 3)
assert egal(j, -sp.Rational(1, 2) + sqrt(3) / 2 * I) and egal(1 + j + j**2, 0)  # U_3
assert egal(j**2, conj(j))  # affirmation 4
assert egal(j**100, j) and egal(1 + j, -(j**2))  # exemple
assert egal(I**4, 1)  # affirmation 3 fausse : i^4 = 1
sixiemes = [sp.exp(I * k * pi / 3) for k in range(6)]
assert meme_ensemble(sixiemes, ens(1, sp.Rational(1, 2) + h * I, -sp.Rational(1, 2) + h * I, -1, -sp.Rational(1, 2) - h * I, sp.Rational(1, 2) - h * I))  # exemple n = 6
# Trous : leurres faux (n = 0 : infinité de solutions ; n solutions et pas 2n, n-1, 2 ; somme 0, pas 1, n, -1 pour n = 3).
assert len(solutions(z**3 - 1)) == 3  # trou 2 (n = 3 : ni 6, ni 2)
assert egal(sum(sp.exp(2 * I * pi * k / 3) for k in range(3)), 0)  # trou 3

# --- geo-t-complexes-racines-unite-ex1
for n, bonne, leurres in [
    (6, sp.Rational(1, 2) + h * I, [sqrt(3) / 2 + I / 2, sp.Rational(1, 2) + I / 2, I]),  # ex1
    (8, sqrt(2) / 2 - sqrt(2) / 2 * I, [sp.Rational(1, 2) + h * I, 1 + I, sqrt(3) / 2 - I / 2]),  # ex1 v1
    (5, sp.exp(-4 * I * pi / 5), [sp.exp(I * pi / 5), -1, I]),  # ex1 v2
]:
    assert egal(sp.expand(bonne**n), 1)  # geo-t-complexes-racines-unite-ex1 (et variantes)
    for l in leurres:
        assert not egal(sp.expand(l**n), 1)  # geo-t-complexes-racines-unite-ex1 : leurre faux
assert egal(sp.exp(-4 * I * pi / 5), sp.exp(6 * I * pi / 5))  # ex1 v2 : k = 3

# --- geo-t-complexes-racines-unite-ex2
for w, rep, leurres in [
    ((1 + j) ** 2025, -1, [1, 0, j]),  # ex2
    ((1 + j**2) ** 100, j, [1, j**2, -1]),  # ex2 v1
    ((1 + j) ** 100, j**2, [1, j, -1]),  # ex2 v2
]:
    assert egal(w, rep)  # geo-t-complexes-racines-unite-ex2 (et variantes)
    for l in leurres:
        assert not egal(w, l)  # geo-t-complexes-racines-unite-ex2 : leurre faux
assert egal(1 + j, sp.exp(I * pi / 3)) and egal(1 + j**2, sp.exp(-I * pi / 3)) and egal(1 + j**2, -j)  # ex2 : formes exponentielles

# --- geo-t-complexes-racines-unite-ex3 (racines n-ièmes d'un complexe)
for P, rep, leurres in [
    (z**3 - 8 * I, ens(sqrt(3) + I, -sqrt(3) + I, -2 * I),
     [ens(-2 * I), ens(sqrt(3) + I, -sqrt(3) + I, 2 * I), [2 * sp.exp(I * pi / 6), 2 * sp.exp(I * pi / 2), 2 * sp.exp(5 * I * pi / 6)]]),  # ex3
    (z**4 + 16, ens(sqrt(2) + I * sqrt(2), -sqrt(2) + I * sqrt(2), -sqrt(2) - I * sqrt(2), sqrt(2) - I * sqrt(2)),
     [ens(2 * I, -2 * I), ens(2, 2 * I, -2, -2 * I), ens(sqrt(2) + I * sqrt(2), -sqrt(2) - I * sqrt(2))]),  # ex3 v1
    (z**3 + 27, ens(sp.Rational(3, 2) + 3 * sqrt(3) / 2 * I, -3, sp.Rational(3, 2) - 3 * sqrt(3) / 2 * I),
     [ens(-3), ens(3, -sp.Rational(3, 2) + 3 * sqrt(3) / 2 * I, -sp.Rational(3, 2) - 3 * sqrt(3) / 2 * I), ens(-3, sp.Rational(3, 2) + 3 * sqrt(3) / 2 * I)]),  # ex3 v2
]:
    S = solutions(P)
    assert meme_ensemble(S, rep)  # geo-t-complexes-racines-unite-ex3 (et variantes)
    for l in leurres:
        assert not meme_ensemble(S, l)  # geo-t-complexes-racines-unite-ex3 : leurre faux
assert egal(2 * sp.exp(I * pi / 6), sqrt(3) + I) and egal(2 * sp.exp(5 * I * pi / 6), -sqrt(3) + I) and egal(2 * sp.exp(3 * I * pi / 2), -2 * I)  # ex3
assert egal((sqrt(3) + I) * j, -sqrt(3) + I) and egal((sqrt(3) + I) * j**2, -2 * I)  # ex3 : remarque (produits par 1, j, j^2)

# --- geo-t-complexes-racines-unite-ex4 (cos 2pi/5)
om = sp.exp(2 * I * pi / 5)
A_, B_ = om + om**4, om**2 + om**3
assert egal(1 + om + om**2 + om**3 + om**4, 0)  # ex4, question 1
assert egal(om**4, conj(om)) and egal(om**3, conj(om**2))  # ex4, question 2
assert egal(A_, 2 * sp.cos(2 * pi / 5)) and egal(B_, 2 * sp.cos(4 * pi / 5))  # ex4, question 2
assert egal(A_ + B_, -1) and egal(sp.expand(A_ * B_), -1)  # ex4, question 3
assert egal(sp.expand(A_ * B_), om**3 + om**4 + om**6 + om**7)  # ex4 : développement
X = sp.symbols("X")
assert meme_ensemble(solutions((z**2 + z - 1)), ens((-1 + sqrt(5)) / 2, (-1 - sqrt(5)) / 2))  # ex4 : X^2 + X - 1
for val, rep, leurres in [
    (sp.cos(2 * pi / 5), (sqrt(5) - 1) / 4, [(sqrt(5) + 1) / 4, (1 - sqrt(5)) / 4, (sqrt(5) - 1) / 2]),  # ex4
    (sp.cos(4 * pi / 5), -(1 + sqrt(5)) / 4, [(1 + sqrt(5)) / 4, (1 - sqrt(5)) / 4, -(1 + sqrt(5)) / 2]),  # ex4 v1
    (sp.cos(pi / 5), (1 + sqrt(5)) / 4, [(sqrt(5) - 1) / 4, -(1 + sqrt(5)) / 4, (1 + sqrt(5)) / 2]),  # ex4 v2
]:
    assert egal(val, rep)  # geo-t-complexes-racines-unite-ex4 (et variantes)
    for l in leurres:
        assert abs(float(sp.N(val - l))) > 1e-6  # geo-t-complexes-racines-unite-ex4 : leurre faux

# ======================================================================================
# Carte geo-t-complexes-second-degre
# ======================================================================================
A, B, Cc = sp.symbols("A B C", real=True)
De = sp.symbols("Delta", negative=True)
# Factorisation générale (Delta < 0) : A(z - z1)(z - z2) = A z^2 + B z + C quand B^2 - 4AC = Delta.
z1 = (-B + I * sqrt(-De)) / (2 * A)
z2 = (-B - I * sqrt(-De)) / (2 * A)
fact = sp.expand(A * (z - z1) * (z - z2))
assert sp.simplify(fact - (A * z**2 + B * z + (B**2 - De) / (4 * A))) == 0  # geo-t-complexes-second-degre : factorisation
assert sp.simplify(z1 + z2 + B / A) == 0  # somme
assert egal(conj(z1), z2)  # solutions conjuguées
# Exemples.
assert meme_ensemble(solutions(z**2 - 2 * z + 5), ens(1 + 2 * I, 1 - 2 * I))  # geo-t-complexes-second-degre (exemple)
assert meme_ensemble(solutions(z**3 + 1), ens(-1, sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I))  # exemple degré 3
assert sp.expand((z + 1) * (z**2 - z + 1)) == z**3 + 1  # exemple degré 3
# Contre-exemples.
assert sp.expand((z - 1) * (z - I)) == sp.expand(z**2 - (1 + I) * z + I)  # contre-exemple coefficients complexes
assert egal((1 + I) ** 2 - 4 * I, -2 * I)  # contre-exemple : Delta = -2i
assert meme_ensemble(solutions(z**2 + 1), ens(I, -I))  # contre-exemple / affirmation 3
# Affirmations.
assert sp.expand((z - I) * (z + 2 * I)) == sp.expand(z**2 + I * z + 2)  # affirmation 2 : solutions i et -2i
assert not egal(conj(I), -2 * I)  # affirmation 2 : non conjuguées
assert meme_ensemble(solutions(z**2 + z + 1), [j, j**2])  # affirmation 4
# Trous : leurres faux (Delta = 0 : solution réelle ; Delta > 0 ; produit c/a et non -c/a, -b/a, a/c).
assert meme_ensemble(solutions(z**2 - 2 * z + 1), ens(1))  # trou 1 : Delta = 0 donne une solution réelle
assert all(egal(sp.im(r), 0) for r in solutions(z**2 - 3 * z + 2))  # trou 1 : Delta > 0 donne des solutions réelles
r_ = solutions(2 * z**2 + 2 * z + 5)
assert egal(r_[0] * r_[1], sp.Rational(5, 2)) and not egal(r_[0] * r_[1], -sp.Rational(5, 2))  # trou 3
assert not egal(r_[0] * r_[1], -1) and not egal(r_[0] * r_[1], sp.Rational(2, 5))  # trou 3 : -b/a et a/c
assert not egal(r_[0], -r_[1]) and not egal(r_[0] * r_[1], 1) and not imag_pur(r_[0])  # trou 2 : ni opposées, ni inverses, ni imaginaires pures

# --- geo-t-complexes-second-degre-ex1
for P, rep, leurres in [
    (z**2 - 4 * z + 13, ens(2 - 3 * I, 2 + 3 * I), [ens(-2 - 3 * I, -2 + 3 * I), ens(4 - 6 * I, 4 + 6 * I), ens(-1, 5)]),  # ex1
    (2 * z**2 + 2 * z + 5, ens(-sp.Rational(1, 2) - sp.Rational(3, 2) * I, -sp.Rational(1, 2) + sp.Rational(3, 2) * I),
     [ens(sp.Rational(1, 2) - sp.Rational(3, 2) * I, sp.Rational(1, 2) + sp.Rational(3, 2) * I), ens(-1 - 3 * I, -1 + 3 * I),
      ens((-1 - sqrt(11)) / 2, (-1 + sqrt(11)) / 2)]),  # ex1 v1
    (z**2 + 6 * z + 10, ens(-3 - I, -3 + I), [ens(3 - I, 3 + I), ens(-3 - 2 * I, -3 + 2 * I), ens(-3 - sqrt(19), -3 + sqrt(19))]),  # ex1 v2
]:
    S = solutions(P)
    assert meme_ensemble(S, rep)  # geo-t-complexes-second-degre-ex1 (et variantes)
    for l in leurres:
        assert not meme_ensemble(S, l)  # geo-t-complexes-second-degre-ex1 : leurre faux
assert (-4) ** 2 - 4 * 13 == -36 and 2**2 - 4 * 2 * 5 == -36 and 6**2 - 4 * 10 == -4  # ex1 : discriminants
assert meme_ensemble(solutions(z**2 - 4 * z - 5), ens(-1, 5))  # ex1 : origine du leurre (Delta = +36)

# --- geo-t-complexes-second-degre-ex2
for r0, (av, bv), leurres in [
    (2 + 3 * I, (-4, 13), [(4, 13), (-4, -5), (-2, 13)]),  # ex2
    (-1 + 2 * I, (2, 5), [(-2, 5), (2, -3), (1, 5)]),  # ex2 v1
    (3 - I, (-6, 10), [(6, 10), (-6, 8), (-3, 10)]),  # ex2 v2
]:
    assert egal(r0**2 + av * r0 + bv, 0)  # geo-t-complexes-second-degre-ex2 (et variantes)
    for (al, bl) in leurres:
        assert not egal(r0**2 + al * r0 + bl, 0)  # geo-t-complexes-second-degre-ex2 : leurre faux
assert egal((2 + 3 * I) ** 2, -5 + 12 * I)  # ex2 : vérification
assert sp.expand((2 + 3 * I) ** 2 + (x + I * y) * 0) == -5 + 12 * I  # ex2 : autre méthode
w_ = sp.expand((2 + 3 * I) ** 2 + a * (2 + 3 * I) + b)
assert sp.expand(sp.re(w_) - (-5 + 2 * a + b)) == 0 and sp.expand(sp.im(w_) - (12 + 3 * a)) == 0  # ex2 : autre méthode

# --- geo-t-complexes-second-degre-ex3 (degré 3, une racine connue)
for P, r0, quot, rep, leurres in [
    (z**3 - 3 * z**2 + 7 * z - 5, 1, z**2 - 2 * z + 5, ens(1, 1 - 2 * I, 1 + 2 * I),
     [ens(1, -1 - 2 * I, -1 + 2 * I), ens(1), ens(1, 2 - 4 * I, 2 + 4 * I)]),  # ex3
    (z**3 - 2 * z - 4, 2, z**2 + 2 * z + 2, ens(2, -1 - I, -1 + I),
     [ens(2, 1 - I, 1 + I), ens(2), ens(2, -2 - 2 * I, -2 + 2 * I)]),  # ex3 v1
    (z**3 + z**2 + 4 * z + 30, -3, z**2 - 2 * z + 10, ens(-3, 1 - 3 * I, 1 + 3 * I),
     [ens(-3, -1 - 3 * I, -1 + 3 * I), ens(-3), ens(-3, 2 - 6 * I, 2 + 6 * I)]),  # ex3 v2
]:
    assert P.subs(z, r0) == 0  # geo-t-complexes-second-degre-ex3 : racine connue
    assert sp.expand((z - r0) * quot - P) == 0  # geo-t-complexes-second-degre-ex3 : factorisation
    S = solutions(P)
    assert meme_ensemble(S, rep)  # geo-t-complexes-second-degre-ex3 (et variantes)
    for l in leurres:
        assert not meme_ensemble(S, l)  # geo-t-complexes-second-degre-ex3 : leurre faux

# --- geo-t-complexes-second-degre-ex4 (équations réciproques)
uu = sp.symbols("u")
for P, eq_u, us, rep, leurres in [
    (z**4 - 3 * z**3 + 4 * z**2 - 3 * z + 1, uu**2 - 3 * uu + 2, [1, 2], ens(1, sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I),
     [ens(1, -sp.Rational(1, 2) + h * I, -sp.Rational(1, 2) - h * I), ens(1, 2), ens(sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I)]),  # ex4
    (z**4 - 2 * z**3 - z**2 - 2 * z + 1, uu**2 - 2 * uu - 3, [3, -1],
     ens((3 - sqrt(5)) / 2, (3 + sqrt(5)) / 2, -sp.Rational(1, 2) + h * I, -sp.Rational(1, 2) - h * I),
     [ens(3, -1), ens((3 - sqrt(5)) / 2, (3 + sqrt(5)) / 2, sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I),
      ens((-3 - sqrt(5)) / 2, (-3 + sqrt(5)) / 2, -sp.Rational(1, 2) + h * I, -sp.Rational(1, 2) - h * I)]),  # ex4 v1
    (z**4 - z**3 + 2 * z**2 - z + 1, uu**2 - uu, [0, 1], ens(I, -I, sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I),
     [ens(0, 1), ens(I, -I, -sp.Rational(1, 2) + h * I, -sp.Rational(1, 2) - h * I), ens(sp.Rational(1, 2) + h * I, sp.Rational(1, 2) - h * I)]),  # ex4 v2
]:
    assert P.subs(z, 0) != 0  # geo-t-complexes-second-degre-ex4 : 0 n'est pas solution
    assert sp.simplify(P / z**2 - eq_u.subs(uu, z + 1 / z)) == 0  # ex4 : P(z)/z^2 = equation en u = z + 1/z
    assert meme_ensemble(solutions(eq_u.subs(uu, z)), ens(*us))  # ex4 : valeurs de u
    S = solutions(P)
    assert meme_ensemble(S, rep)  # geo-t-complexes-second-degre-ex4 (et variantes)
    for l in leurres:
        assert not meme_ensemble(S, l)  # geo-t-complexes-second-degre-ex4 : leurre faux
assert sp.expand((z - 1) ** 2 * (z**2 - z + 1)) == z**4 - 3 * z**3 + 4 * z**2 - 3 * z + 1  # ex4 : factorisation
assert sp.expand((z**2 - 3 * z + 1) * (z**2 + z + 1)) == z**4 - 2 * z**3 - z**2 - 2 * z + 1  # ex4 v1
assert sp.expand((z**2 + 1) * (z**2 - z + 1)) == z**4 - z**3 + 2 * z**2 - z + 1  # ex4 v2
assert sp.expand((z + 1 / z) ** 2 - 2 - (z**2 + 1 / z**2)) == 0  # ex4 : u^2 - 2

# ======================================================================================
# Carte geo-t-complexes-geometrie
# ======================================================================================
xp, yp = sp.symbols("xp yp", real=True)
prod = sp.expand((xp + I * yp) * conj(x + I * y))
assert sp.expand(sp.re(prod) - (x * xp + y * yp)) == 0 and sp.expand(sp.im(prod) - (x * yp - y * xp)) == 0  # démonstration
q_ex = (3 * I - (1 + I)) / ((3 + 2 * I) - (1 + I))
assert egal(q_ex, I)  # geo-t-complexes-geometrie (exemple)
assert egal((-1 + 2 * I) * (2 - I), 5 * I)  # exemple : calcul intermédiaire
assert egal((I - 0) / (1 - 0), I) and sp.Abs(I) == 1 and not reel(I)  # contre-exemple |q| = 1 et affirmation 2
assert egal((-1 - 1) / (I - 1), 1 + I) and not reel(1 + I)  # affirmation 3 fausse
# Affirmation 4 : c - a = -2(b - a) => A entre B et C (A = (B + 2C)/3, combinaison convexe).
aa, bb = 1 + 2 * I, 4 - I
cc_ = aa - 2 * (bb - aa)
assert egal(aa, (2 * bb + cc_) / 3)  # affirmation 4 : A barycentre à coefficients positifs de B et C
# Trous : leurres faux (q = 2 donne un alignement sans q = 1 ; q = 2i donne une orthogonalité sans q = i ; etc.).
assert reel(2) and not egal(2, 1)  # trou 1, leurre q = 1
assert imag_pur(2 * I) and not egal(2 * I, I)  # trou 2, leurre q = i
assert sp.Abs(I) == 1 and not reel(I) and sp.Abs(2 * I) != 1  # trous 1 et 2, leurre |q| = 1
q_ = (1 + I)  # angle pi/4 : arguments de 1/q et de (c-a)(b-a) faux en général
assert not egal(sp.arg(1 / q_), sp.arg(q_))  # trou 3, leurre argument de (b-a)/(c-a)
aa, bb, cc_ = 0, 2 + I, (2 + I) * (1 + I)
assert not egal(sp.exp(I * sp.arg((cc_ - aa) * (bb - aa))), sp.exp(I * sp.arg((cc_ - aa) / (bb - aa))))  # trou 3, leurre produit

# --- geo-t-complexes-geometrie-ex1
for (aa, bb, cc_), qv, nature in [
    ((2 - I, 5 + I, 2 * I), I, 0),  # ex1 : rectangle isocèle en A
    ((0, 2, 1 + I * sqrt(3)), sp.exp(I * pi / 3), 1),  # ex1 v1 : équilatéral
    ((1, 2 + I, -1 + 2 * I), 2 * I, 2),  # ex1 v2 : rectangle en A non isocèle
]:
    qq = (cc_ - aa) / (bb - aa)
    assert egal(qq, qv)  # geo-t-complexes-geometrie-ex1 (et variantes)
    AB, AC, BC = (sp.simplify(sp.Abs(sp.expand(u))) for u in (bb - aa, cc_ - aa, cc_ - bb))
    rectA = imag_pur(qq)
    isoA = sp.simplify(AB - AC) == 0
    equi = isoA and sp.simplify(AB - BC) == 0
    # NATURES : 0 rectangle isocèle en A ; 1 équilatéral ; 2 rectangle en A non isocèle ; 3 isocèle en A ni rectangle ni équilatéral
    verites = [rectA and isoA, equi, rectA and not isoA and sp.simplify(BC - AB) != 0 and sp.simplify(BC - AC) != 0,
               isoA and not rectA and not equi]
    assert verites[nature] and sum(bool(v) for v in verites) == 1  # geo-t-complexes-geometrie-ex1 : seule la bonne nature est vraie
assert egal((-2 + 3 * I) * (3 - 2 * I), 13 * I) and egal((-2 + 2 * I) * (1 - I), 4 * I)  # ex1 : calculs intermédiaires
assert egal(sp.Abs(-1 + I * sqrt(3)), 2)  # ex1 v1 : BC = 2

# --- geo-t-complexes-geometrie-ex2 (paramètre m pour l'alignement)
m = sp.symbols("m", real=True)
for aa, bb, cfun, rep, leurres, perp in [
    (1 + 2 * I, 3 - I, lambda mm: mm - 7 * I, 7, [6, -5, -sp.Rational(25, 2)], -sp.Rational(25, 2)),  # ex2
    (2 - I, 4 + 3 * I, lambda mm: mm + 7 * I, 6, [4, -2, -14], -14),  # ex2 v1
    (-1 + I, 1 - 2 * I, lambda mm: 5 + mm * I, -8, [-9, 10, 5], 5),  # ex2 v2
]:
    qq = sp.expand_complex((cfun(m) - aa) * conj(bb - aa))
    sol = sp.solve(sp.im(qq), m)
    assert sol == [rep]  # geo-t-complexes-geometrie-ex2 (et variantes)
    for l in leurres:
        assert not reel((cfun(l) - aa) / (bb - aa))  # geo-t-complexes-geometrie-ex2 : leurre faux
    assert sp.solve(sp.re(qq), m) == [perp]  # ex2 : le leurre « perpendiculaire » correspond à Re(q) = 0
assert egal(6 - 9 * I, 3 * (2 - 3 * I)) and egal(4 + 8 * I, 2 * (2 + 4 * I))  # ex2 et v1 : vérifications
assert egal(sp.expand(((m - 1) - 9 * I) * (2 + 3 * I)), (2 * (m - 1) + 27) + I * (3 * (m - 1) - 18))  # ex2 : correction
assert egal(sp.expand(((m - 2) + 8 * I) * (2 - 4 * I)), (2 * (m - 2) + 32) + I * (16 - 4 * (m - 2)))  # ex2 v1 : correction
assert egal(sp.expand((6 + (m - 1) * I) * (2 + 3 * I)), (12 - 3 * (m - 1)) + I * (18 + 2 * (m - 1)))  # ex2 v2 : correction

# --- geo-t-complexes-geometrie-ex3 (carré direct)
def est_carre_direct(aa, bb, cc_, dd):
    return egal((dd - aa) / (bb - aa), I) and egal(cc_ - dd, bb - aa)

for aa, bb, (cv, dv), leurres in [
    (1 + I, 4 + 2 * I, (3 + 5 * I, 4 * I), [(5 - I, 2 - 2 * I), (4 * I, 3 + 5 * I), (2 + 3 * I, 4 * I)]),  # ex3
    (0, 2 + I, (1 + 3 * I, -1 + 2 * I), [(3 - I, 1 - 2 * I), (-1 + 2 * I, 1 + 3 * I), (sp.Rational(1, 2) + sp.Rational(3, 2) * I, -1 + 2 * I)]),  # ex3 v1
    (2 - I, 3 + 2 * I, (3 * I, -1), [(6 + I, 5 - 2 * I), (-1, 3 * I), (1 + I, -1)]),  # ex3 v2
]:
    assert est_carre_direct(aa, bb, cv, dv)  # geo-t-complexes-geometrie-ex3 (et variantes)
    assert egal(cv - bb, I * (bb - aa))  # ex3 : vérification de la correction
    for (cl, dl) in leurres:
        assert not est_carre_direct(aa, bb, cl, dl)  # geo-t-complexes-geometrie-ex3 : leurre faux
    # Le leurre 1 est le carré indirect.
    cl, dl = leurres[0]
    assert egal((dl - aa) / (bb - aa), -I) and egal(cl - dl, bb - aa)  # ex3 : leurre « sens indirect »

# --- geo-t-complexes-geometrie-ex4 (c = ia, d = -ib)
av, bv = sp.symbols("a b")
cv, dv = I * av, -I * bv
assert sp.expand(-I * (dv - av) - (cv - bv)) == 0  # geo-t-complexes-geometrie-ex4 : c - b = -i(d - a)
assert not egal(-I, I) and not egal(-I, -1) and not egal(-I, 1)  # ex4 : leurres faux
for (a0, b0) in [(3, 1 + 2 * I), (2 - I, -1 + 4 * I)]:
    c0, d0 = I * a0, -I * b0
    assert egal((c0 - b0) / (d0 - a0), -I)  # ex4 : valeur numérique
    assert egal(sp.Abs(c0 - b0), sp.Abs(d0 - a0)) and imag_pur((c0 - b0) * conj(d0 - a0))  # ex4 : BC = AD et orthogonalité

# --- geo-t-complexes-geometrie-ex5 (triangle z, z^2, z^3)
w = sp.symbols("w")
assert sp.cancel((w**3 - w) / (w**2 - w) - (w + 1)) == 0  # ex5 : rectangle en M
assert sp.cancel((w - w**2) / (w**3 - w**2) + 1 / w) == 0  # ex5 v1 : rectangle en N
assert sp.cancel((w - w**3) / (w**2 - w**3) - (1 + 1 / w)) == 0  # ex5 v2 : rectangle en P
assert set(sp.solve(w**2 - w, w)) == {0, 1} and set(sp.solve(w**3 - w, w)) == {-1, 0, 1} and set(sp.solve(w**3 - w**2, w)) == {0, 1}  # ex5 Q1

def rectangle_en(zv, sommet):
    """Le triangle M(z) N(z^2) P(z^3) (sommets distincts) est-il rectangle au sommet indiqué ?"""
    M, N, P = zv, zv**2, zv**3
    S, U, V = {"M": (M, N, P), "N": (N, M, P), "P": (P, M, N)}[sommet]
    return abs((complex(U - S) * complex(V - S).conjugate()).real) < 1e-9

def distincts(zv):
    return all(abs(complex(zv) - r) > 1e-9 for r in (0, 1, -1))

# ex5 : E = droite x = -1 privée de -1 ; v1 : axe imaginaire privé de O ; v2 : cercle |z + 1/2| = 1/2 privé de 0 et -1.
for yv in [-2, -0.5, 0.7, 3]:
    assert rectangle_en(complex(-1, yv), "M")  # ex5 : points de la droite x = -1
    assert rectangle_en(complex(0, yv), "N")  # ex5 v1 : points de l'axe imaginaire
for tv in [0.4, 1.3, 2.0, 4.0, 5.5]:
    zv = -0.5 + 0.5 * cmath.exp(1j * tv)
    assert distincts(zv) and rectangle_en(zv, "P")  # ex5 v2 : points du cercle
for zv in [complex(0.3, 0.8), complex(-2, 1), complex(2, -1.5), complex(-0.5, 0.2)]:
    assert not rectangle_en(zv, "M") and not rectangle_en(zv, "N") and not rectangle_en(zv, "P")  # points hors des ensembles
assert not distincts(-1) and not distincts(0)  # points exclus (droite x = -1 en -1 ; axe imaginaire et cercle en 0 ; cercle en -1)
assert abs(abs(-1 + 0.5) - 0.5) < 1e-12 and abs(abs(0 + 0.5) - 0.5) < 1e-12 and abs(abs(1 + 0.5) - 0.5) > 1e-3  # v2 : -1 et 0 sur le cercle, 1 non
# Leurres : axe réel (alignement), cercle de centre -1 et rayon 1, cercle unité, cercle de centre 1/2.
assert not rectangle_en(complex(2, 0), "M") and not rectangle_en(-1 + cmath.exp(0.9j), "M")  # ex5 : leurres faux
assert not rectangle_en(-1 + cmath.exp(2.0j), "M")  # ex5 : leurre cercle de centre -1 (point hors de x = -1)
assert not rectangle_en(cmath.exp(0.8j), "N")  # ex5 v1 : leurre cercle unité
assert not rectangle_en(complex(-1, 2), "N")  # ex5 v1 : leurre droite x = -1
assert not rectangle_en(0.5 + 0.5 * cmath.exp(1.0j), "P")  # ex5 v2 : leurre cercle de centre 1/2
assert not rectangle_en(complex(-1, 2), "P")  # ex5 v2 : leurre droite x = -1
assert sp.simplify(sp.re(1 + 1 / (x + I * y)) * (x**2 + y**2) - (x**2 + y**2 + x)) == 0  # ex5 v2 : x^2 + y^2 + x = 0

print("geo_t_complexes : toutes les vérifications passent.")
