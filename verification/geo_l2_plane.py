"""Vérifications calculatoires du chapitre « geo-l2-plane » (géométrie, L2 : triangle et cercle) :
cartes, exercices et QCM de content/*/geometrie/l2-geo-l2-plane.json et content/qcm/geo-l2-plane.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/geo_l2_plane.py (doit finir sans erreur).
"""

import cmath
import math
import random
from fractions import Fraction as F

from sympy import I, Rational, conjugate, cos, expand, im, pi, re, simplify, sin, sqrt, symbols

random.seed(2027)
EPS = 1e-9


# ---------------------------------------------------------------------------
# Outils (nombres complexes en flottants pour les angles, fractions pour l'exact)
# ---------------------------------------------------------------------------
def ang(u, v):
    """Mesure dans ]-pi, pi] de l'angle orienté (u, v), u et v complexes non nuls."""
    return cmath.phase(v / u)


def mod_pi(x):
    return (x + math.pi / 2) % math.pi - math.pi / 2


def egal_mod_pi(x, y):
    return abs(mod_pi(x - y)) < 1e-9


def geo(u, v):
    """Angle géométrique entre u et v, dans [0, pi]."""
    return abs(ang(u, v))


def cercle3(a, b, c):
    """Centre et rayon du cercle passant par a, b, c (complexes non alignés)."""
    ax, ay, bx, by, cx, cy = a.real, a.imag, b.real, b.imag, c.real, c.imag
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    o = complex(ux, uy)
    return o, abs(a - o)


def birapport(a, b, c, d):
    return (c - a) / (c - b) * (d - b) / (d - a)


def orthocentre(a, b, c):
    o, _ = cercle3(a, b, c)
    return a + b + c - 2 * o


def dot(u, v):
    return (u.conjugate() * v).real


# ===========================================================================
# Carte geo-l2-plane-angles-inscrits
# ===========================================================================
# Théorème (angle au centre, angle inscrit) testé sur des points aléatoires d'un cercle.
for _ in range(200):
    o = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    r = random.uniform(0.5, 4)
    al, be, mu, nu = (random.uniform(0, 2 * math.pi) for _ in range(4))
    a, b, m, n = (o + r * cmath.exp(1j * t) for t in (al, be, mu, nu))
    centre = ang(a - o, b - o)
    assert abs(((2 * ang(a - m, b - m) - centre + math.pi) % (2 * math.pi)) - math.pi) < 1e-7
    assert egal_mod_pi(ang(a - m, b - m), ang(a - n, b - n))
    # Angle tangentiel (cours) : (T, AB) = (MA, MB) mod pi, T dirigée par i(a - o).
    assert egal_mod_pi(ang(1j * (a - o), b - a), ang(a - m, b - m))
    # Forme géométrique (cours) : M du même côté de (AB) que O => AMB = AOB / 2, sinon pi - AOB / 2.
    cote = lambda p: ((b - a).conjugate() * (p - a)).imag
    g = geo(a - m, b - m)
    if cote(m) * cote(o) > 0:
        assert abs(g - geo(a - o, b - o) / 2) < 1e-7
    elif cote(m) * cote(o) < 0:
        assert abs(g - (math.pi - geo(a - o, b - o) / 2)) < 1e-7
# Démonstration : (1 - iy)/(-1 - iy) = -exp(-2i arctan y), de mesure pi - 2 arctan y.
for y in (-5, -1, -0.3, 0, 0.7, 2, 10):
    assert abs((1 - 1j * y) / (-1 - 1j * y) + cmath.exp(-2j * math.atan(y))) < 1e-12
# Exemple : A(1,0), B(0,1), M(-1,0) : AOB = pi/2, AMB = pi/4.
A, B, M = 1, 1j, -1
assert abs(geo(A, B) - math.pi / 2) < EPS and abs(geo(A - M, B - M) - math.pi / 4) < EPS
# Trous : configuration A = 1, B = i, M = -1 du cercle unité ; chaque leurre du trou 1 est faux.
th = ang(A - M, B - M)
assert abs(ang(A, B) - 2 * th) < EPS
for leurre in (th, math.pi - th, th + math.pi / 2):
    assert abs(((leurre - ang(A, B) + math.pi) % (2 * math.pi)) - math.pi) > 0.1
# Trou 3 et pièges : cocycliques avec (CA,CB) = pi/4 : la condition (DB,DA) et le double sont faux.
a, b, c, d = 1, 1j, -1, -1j
t1, t2 = ang(a - c, b - c), ang(a - d, b - d)
assert egal_mod_pi(t1, t2) and not egal_mod_pi(t1, -t2) and not egal_mod_pi(t1, 2 * t2)
# Angles géométriques : C, D sur des arcs opposés : angles supplémentaires, non égaux (trou 3, leurre 1).
a, b, c, d = 1, 1j, -1, cmath.exp(1j * math.pi / 4)
assert abs(geo(a - c, b - c) - math.pi / 4) < EPS and abs(geo(a - d, b - d) - 3 * math.pi / 4) < EPS
# Contre-exemple et pièges : A(-1,0), B(1,0), C(0,2), D(0,-2).
a, b, c, d = -1, 1, 2j, -2j
assert abs(geo(a - c, b - c) - geo(a - d, b - d)) < EPS                       # angles géométriques égaux
assert abs(abs(c - a) * abs(c - b) - abs(d - a) * abs(d - b)) < EPS          # produits égaux
assert egal_mod_pi(ang(a - c, b - c), ang(b - d, a - d))                     # (CA,CB) = (DB,DA)
assert abs(ang(a - c, b - c) - math.atan(4 / 3)) < EPS and abs(ang(a - d, b - d) + math.atan(4 / 3)) < EPS
o, r = cercle3(a, b, c)
assert abs(o - 0.75j) < EPS and abs(r - 1.25) < EPS and abs(abs(d - o) - 2.75) < EPS
assert abs(birapport(a, b, c, d).imag) > 0.1                                  # ni cocycliques ni alignés
# Affirmation 2 : M(-1,0), N(sqrt2/2, sqrt2/2) : pi/4 et 3pi/4.
n = cmath.exp(1j * math.pi / 4)
assert abs(geo(1 - (-1), 1j - (-1)) - math.pi / 4) < EPS and abs(geo(1 - n, 1j - n) - 3 * math.pi / 4) < EPS
# Affirmation 4 : parallélogramme A(0,0), B(2,0), C(3,1), D(1,1) : angles en A et C valant pi/4, non inscrit.
a, b, c, d = 0, 2, 3 + 1j, 1 + 1j
assert abs(geo(b - a, d - a) - math.pi / 4) < EPS and abs(geo(b - c, d - c) - math.pi / 4) < EPS
assert abs(birapport(a, b, c, d).imag) > 0.1

# geo-l2-plane-angles-inscrits-ex1 (+ v1, v2) : AMB = t/2, ANB = 180 - t/2 ; leurres distincts et faux.
for t in (110, 80, 150):
    o = 0
    a, b = cmath.exp(1j * math.radians(90 - t / 2)), cmath.exp(1j * math.radians(90 + t / 2))  # AOB = t, O sous (AB)
    m, n = cmath.exp(-1j * math.pi / 2), cmath.exp(1j * math.pi / 2)                          # M côté O, N de l'autre
    assert abs(math.degrees(geo(a, b)) - t) < 1e-7
    assert abs(math.degrees(geo(a - m, b - m)) - t / 2) < 1e-7
    assert abs(math.degrees(geo(a - n, b - n)) - (180 - t / 2)) < 1e-7
    vrai = (t / 2, 180 - t / 2)
    for l in ((t / 2, t / 2), (t, 180 - t), (2 * t, 360 - 2 * t)):
        assert l != vrai

# geo-l2-plane-angles-inscrits-ex2 (+ v1, v2) : D(t, yD) sur le cercle de diamètre [BC].
for Bx, Cy, yD, sols, leurres in [(6, 8, 7, {-1, 7}, [{7}, {1, 5}, {-2, 8}]),
                                  (8, 6, -1, {1, 7}, [{7}, {-1, 9}, {3, 5}]),
                                  (10, 24, 7, {-7, 17}, [{17}, {-8, 18}, {0, 10}])]:
    cx, cy, r2 = F(Bx, 2), F(Cy, 2), F(Bx * Bx + Cy * Cy, 4)
    sol = {t for t in range(-60, 60) if (t - cx) ** 2 + (yD - cy) ** 2 == r2}
    assert sol == sols
    assert all(l != sols for l in leurres)
    # Le cercle de diamètre [BC] passe bien par A, B, C.
    assert cx**2 + cy**2 == r2 and (Bx - cx) ** 2 + cy**2 == r2

# geo-l2-plane-angles-inscrits-ex3 : le symétrique de H par rapport à (BC) est sur le cercle circonscrit.
def sym_droite(p, u, v):
    w = (v - u) / abs(v - u)
    return u + w * w * (p - u).conjugate()


for _ in range(200):
    a, b, c = (complex(random.uniform(-5, 5), random.uniform(-5, 5)) for _ in range(3))
    if abs(((b - a).conjugate() * (c - a)).imag) < 0.5:
        continue
    o, r = cercle3(a, b, c)
    h1 = sym_droite(orthocentre(a, b, c), b, c)
    assert abs(abs(h1 - o) - r) < 1e-6

# geo-l2-plane-angles-inscrits-ex4 (Miquel) : les trois cercles sont concourants.
for _ in range(200):
    a, b, c = (complex(random.uniform(-5, 5), random.uniform(-5, 5)) for _ in range(3))
    if abs(((b - a).conjugate() * (c - a)).imag) < 1:
        continue
    p = b + random.uniform(-0.5, 1.5) * (c - b)
    q = c + random.uniform(-0.5, 1.5) * (a - c)
    rr = a + random.uniform(-0.5, 1.5) * (b - a)
    try:
        o1, r1 = cercle3(a, rr, q)
        o2, r2 = cercle3(b, p, rr)
        o3, r3 = cercle3(c, q, p)
    except ZeroDivisionError:
        continue
    # Second point commun M des deux premiers cercles : symétrique de R par rapport à (o1 o2).
    m = sym_droite(rr, o1, o2)
    if min(abs(m - p), abs(m - q), abs(m - rr)) < 1e-3:
        continue
    assert abs(abs(m - o3) - r3) < 1e-6 * max(1, r3)

# ===========================================================================
# Carte geo-l2-plane-remarquables
# ===========================================================================
# Relation d'Euler et hauteurs, sur des triangles aléatoires.
for _ in range(200):
    a, b, c = (complex(random.uniform(-5, 5), random.uniform(-5, 5)) for _ in range(3))
    if abs(((b - a).conjugate() * (c - a)).imag) < 0.5:
        continue
    o, R = cercle3(a, b, c)
    h = a + b + c - 2 * o
    g = (a + b + c) / 3
    assert abs(dot(h - a, c - b)) < 1e-6 and abs(dot(h - b, a - c)) < 1e-6
    assert abs((h - o) - 3 * (g - o)) < 1e-9
    # Centre du cercle inscrit : barycentre (a, b, c), équidistant des trois côtés.
    la, lb, lc = abs(c - b), abs(a - c), abs(b - a)
    i = (la * a + lb * b + lc * c) / (la + lb + lc)
    dist = lambda p, u, v: abs(((v - u).conjugate() * (p - u)).imag) / abs(v - u)
    assert abs(dist(i, a, b) - dist(i, b, c)) < 1e-9 and abs(dist(i, a, b) - dist(i, c, a)) < 1e-9
    # Cercle d'Euler (ex4) : milieux des côtés, pieds des hauteurs, milieux de [AH] à la distance R/2 du milieu de [OH].
    om = (o + h) / 2
    pied = lambda p, u, v: u + (v - u) * dot(v - u, p - u) / abs(v - u) ** 2
    for pt_ in ((b + c) / 2, (a + c) / 2, (a + b) / 2, (a + h) / 2, (b + h) / 2, (c + h) / 2,
                pied(a, b, c), pied(b, c, a), pied(c, a, b)):
        assert abs(abs(pt_ - om) - R / 2) < 1e-6
    # ex3 : symétrique de H par rapport au milieu de [BC] = point diamétralement opposé à A ; AH = 2 OA'.
    assert abs((b + c - h) - (2 * o - a)) < 1e-9 and abs((h - a) - 2 * ((b + c) / 2 - o)) < 1e-9
# Exemple : A(0,0), B(6,0), C(2,4) : O(3,1), H(2,2), G(8/3,4/3).
a, b, c = 0, 6, 2 + 4j
o, R = cercle3(a, b, c)
assert abs(o - (3 + 1j)) < EPS and abs(orthocentre(a, b, c) - (2 + 2j)) < EPS
# Contre-exemple 1 : triangle A(0,0), B(4,0), C(0,3) : I(1,1), I + IA + IB + IC = (2,1), orthocentre A.
a, b, c = 0, 4, 3j
i = (5 * a + 3 * b + 4 * c) / 12
assert abs(i - (1 + 1j)) < EPS and abs(i + (a - i) + (b - i) + (c - i) - (2 + 1j)) < EPS
assert abs(orthocentre(a, b, c) - a) < EPS
# Trou 1 et piège « O = G » : H := G + (GA+GB+GC) = G n'est pas l'orthocentre (triangle non équilatéral).
g = (a + b + c) / 3
assert abs(g - orthocentre(a, b, c)) > 0.1
# Leurres du trou 1 : dans ce triangle OH = 3 OG, et OH != 2 OG, OG != 3 OH, GH != 3 GO.
o, _ = cercle3(a, b, c)
h = orthocentre(a, b, c)
assert abs((h - o) - 2 * (g - o)) > 0.1 and abs((g - o) - 3 * (h - o)) > 0.1 and abs((h - g) - 3 * (o - g)) > 0.1
# Trou 3 : dans ce triangle, O n'est ni G, ni H, ni I.
assert min(abs(o - g), abs(o - h), abs(o - i)) > 0.1
# Contre-exemple 2 / affirmation 4 / QCM q2 : I hors de la droite d'Euler x + y = 4 pour A(0,0), B(6,0), C(2,4).
s = sqrt(2) * 4 + sqrt(5) * 2 + 6
somme = (12 * sqrt(5) + 36) / s
assert abs(float(somme) - 3.8957) < 1e-3 and abs(float(somme) - 4) > 0.05
la, lb, lc = 4 * math.sqrt(2), 2 * math.sqrt(5), 6
i = (la * 0 + lb * 6 + lc * (2 + 4j)) / (la + lb + lc)
assert abs(i.real + i.imag - float(somme)) < 1e-9
assert abs((2 + 2j).real + (2 + 2j).imag - 4) < EPS and abs(3 + 1 - 4) < EPS     # O et H sur x + y = 4
# Contre-exemple 3 / affirmation 2 : A(0,0), B(4,0), C(1,1) : angle obtus en C, H = (1,3).
a, b, c = 0, 4, 1 + 1j
assert dot(a - c, b - c) < 0 and abs(orthocentre(a, b, c) - (1 + 3j)) < EPS

# geo-l2-plane-remarquables-ex1 (+ v1, v2) : O, H, G et leurres.
for Bx, C, O_, H_, G_, faux in [(8, 2 + 6j, 4 + 2j, 2 + 2j, F(10, 3) + 2j, 6 + 4j),
                               (10, 4 + 6j, 5 + 1j, 4 + 4j, F(14, 3) + 2j, 9 + 5j),
                               (4, 1 + 3j, 2 + 1j, 1 + 1j, F(5, 3) + 1j, 3 + 2j)]:
    o, _ = cercle3(0, Bx, C)
    assert abs(o - O_) < EPS and abs(orthocentre(0, Bx, C) - H_) < EPS
    assert abs((0 + Bx + C) / 3 - complex(G_)) < EPS
    assert abs((0 + Bx + C - o) - faux) < EPS                       # leurre : H = A + B + C - O
    assert len({H_, O_, complex(G_), faux}) == 4

# geo-l2-plane-remarquables-ex2 (+ v1, v2) : centre du cercle inscrit et leurres (poids permutés, O).
for b_, c_, I_, l1, l2, O_ in [(3, 4, (1, 1), (F(4, 3), F(5, 4)), (F(5, 3), F(3, 4)), (2, F(3, 2))),
                               (5, 12, (2, 2), (F(24, 5), F(13, 6)), (F(26, 5), F(5, 6)), (6, F(5, 2))),
                               (8, 15, (3, 3), (F(45, 8), F(17, 5)), (F(51, 8), F(8, 5)), (F(15, 2), 4))]:
    a_ = math.isqrt(b_ * b_ + c_ * c_)
    assert a_ * a_ == b_ * b_ + c_ * c_
    s_ = a_ + b_ + c_
    A_, B_, C_ = (0, 0), (c_, 0), (0, b_)
    bar = lambda wa, wb, wc: (F(wa * A_[0] + wb * B_[0] + wc * C_[0], s_), F(wa * A_[1] + wb * B_[1] + wc * C_[1], s_))
    assert bar(a_, b_, c_) == I_ and bar(b_, c_, a_) == l1 and bar(c_, a_, b_) == l2
    assert (F(c_, 2), F(b_, 2)) == O_
    r_ = I_[0]
    assert F(abs(b_ * r_ + c_ * r_ - b_ * c_), a_) == r_ and F(b_ * c_, 2) / F(s_, 2) == r_
    assert len({I_, l1, l2, O_}) == 4

# ===========================================================================
# Carte geo-l2-plane-relations-metriques
# ===========================================================================
for _ in range(200):
    a, b, c = (complex(random.uniform(-5, 5), random.uniform(-5, 5)) for _ in range(3))
    if abs(((b - a).conjugate() * (c - a)).imag) < 0.5:
        continue
    la, lb, lc = abs(c - b), abs(a - c), abs(b - a)
    A_ = geo(b - a, c - a)
    S = abs(((b - a).conjugate() * (c - a)).imag) / 2
    _, R = cercle3(a, b, c)
    p = (la + lb + lc) / 2
    assert abs(la**2 - (lb**2 + lc**2 - 2 * lb * lc * math.cos(A_))) < 1e-8
    assert abs(S - lb * lc * math.sin(A_) / 2) < 1e-8
    assert abs(la / math.sin(A_) - 2 * R) < 1e-7 and abs(2 * R - la * lb * lc / (2 * S)) < 1e-7
    assert abs(S - math.sqrt(p * (p - la) * (p - lb) * (p - lc))) < 1e-7
    # Pièges : R du cercle inscrit (r = S/p), angle orienté (signe), côté AB à la place de BC.
    assert abs(la / math.sin(A_) - 2 * S / p) > 1e-6
    # Weitzenböck (ex4).
    assert la**2 + lb**2 + lc**2 >= 4 * math.sqrt(3) * S - 1e-9
# Héron : identité de factorisation (démonstration).
x, y, z = symbols("a b c", positive=True)
assert expand((2 * y * z) ** 2 - (y**2 + z**2 - x**2) ** 2 - (x - y + z) * (x + y - z) * (y + z - x) * (x + y + z)) == 0
pp = (x + y + z) / 2
assert expand((x - y + z) * (x + y - z) * (y + z - x) * (x + y + z) - 16 * pp * (pp - x) * (pp - y) * (pp - z)) == 0
# Trous : triangle équilatéral de côté 1 : a^2 = 1 ; leurres 3, 3/2, 2 - sqrt3 ; a/sinA = 2/sqrt3 = 2R.
A_ = pi / 3
assert 1 + 1 - 2 * cos(A_) == 1
assert [simplify(e) for e in (2 + 2 * cos(A_), 2 - cos(A_), 2 - 2 * sin(A_))] == [3, Rational(3, 2), 2 - sqrt(3)]
assert simplify(1 / sin(A_) - 2 / sqrt(3)) == 0
# Trou 3, leurre (1/2) ab sin A : triangle 3-4-5 rectangle en A (contre-exemple « angle non compris »).
assert Rational(1, 2) * 4 * 3 * 1 == 6 and Rational(1, 2) * 5 * 4 * 1 == 10
# Exemple : b = 5, c = 8, A = pi/3.
a2 = 25 + 64 - 2 * 5 * 8 * cos(pi / 3)
assert a2 == 49 and simplify(Rational(1, 2) * 40 * sin(pi / 3) - 10 * sqrt(3)) == 0
assert simplify(7 / (2 * sin(pi / 3)) - 7 * sqrt(3) / 3) == 0 and sqrt(10 * 3 * 5 * 2) == 10 * sqrt(3)
# Contre-exemple « cas ambigu » : A = pi/6, b = 2, a = sqrt2 : c = sqrt3 +/- 1, sin B = sqrt2/2.
cc = symbols("cc")
for cv in (sqrt(3) + 1, sqrt(3) - 1):
    assert simplify(2 - (4 + cv**2 - 2 * 2 * cv * cos(pi / 6))) == 0
assert simplify(2 * sin(pi / 6) / sqrt(2) - sqrt(2) / 2) == 0
assert pi / 6 + 3 * pi / 4 < pi
# Affirmation 3 : a = 2, b = 4, c = 3.
assert 2**2 < 4**2 + 3**2 and 4**2 > 2**2 + 3**2
# Weitzenböck (ex4) : cos A + sqrt3 sin A = 2 cos(A - pi/3).
t = symbols("t")
assert simplify(expand((cos(t) + sqrt(3) * sin(t)) - 2 * cos(t - pi / 3), trig=True)) == 0

# geo-l2-plane-relations-metriques-ex1 (+ v1, v2) : troisième côté ; leurres a^2 = s + 2bc cos, s, s - bc cos.
for b_, c_, Adeg, a_, leurres in [(3, 8, 60, 7, [97, 73, 61]), (5, 3, 120, 7, [19, 34, Rational(83, 2)]),
                                  (7, 15, 60, 13, [379, 274, Rational(443, 2)])]:
    co = cos(pi * Adeg / 180)
    assert b_**2 + c_**2 - 2 * b_ * c_ * co == a_**2
    assert leurres == [b_**2 + c_**2 + 2 * b_ * c_ * co, b_**2 + c_**2, b_**2 + c_**2 - b_ * c_ * co]
    assert all(l != a_**2 for l in leurres)

# geo-l2-plane-relations-metriques-ex2 (+ v1, v2) : loi des sinus.
for Ad, Bd, a_, b_, R_, binv in [(45, 60, 4, 2 * sqrt(6), 2 * sqrt(2), 4 * sqrt(6) / 3),
                                 (30, 45, 3, 3 * sqrt(2), 3, 3 * sqrt(2) / 2),
                                 (60, 45, 6, 2 * sqrt(6), 2 * sqrt(3), 3 * sqrt(6))]:
    sA, sB = sin(pi * Ad / 180), sin(pi * Bd / 180)
    assert simplify(a_ / sA - 2 * R_) == 0 and simplify(a_ * sB / sA - b_) == 0
    assert simplify(a_ * sA / sB - binv) == 0                       # leurre : rapport inversé
    assert simplify(binv - b_) != 0 and R_ != 2 * R_ and R_ != R_ / 2

# geo-l2-plane-relations-metriques-ex3 (+ v1, v2) : Héron, R = abc/(4S), r = S/p ; leurres abc/(2S), S/2.
for a_, b_, c_, S_, R_, r_ in [(13, 14, 15, 84, Rational(65, 8), 4), (5, 5, 6, 12, Rational(25, 8), Rational(3, 2)),
                               (7, 8, 9, 12 * sqrt(5), 21 * sqrt(5) / 10, sqrt(5))]:
    p = Rational(a_ + b_ + c_, 2)
    assert simplify(sqrt(p * (p - a_) * (p - b_) * (p - c_)) - S_) == 0
    assert simplify(Rational(a_ * b_ * c_, 4) / S_ - R_) == 0 and simplify(S_ / p - r_) == 0
    assert simplify(Rational(a_ * b_ * c_, 2) / S_ - 2 * R_) == 0   # leurre R = abc/(2S) = 2R
    assert simplify(Rational(a_ * b_ * c_, 4) / (S_ / 2) - 2 * R_) == 0
    assert simplify(r_ - R_) != 0

# ===========================================================================
# Carte geo-l2-plane-complexes-cercles
# ===========================================================================
zz, w = symbols("z w")
xr, yr = symbols("x y", real=True)
Z = xr + I * yr
# Équation d'un cercle : |z - w|^2 = z zb - wb z - w zb + |w|^2.
om = symbols("om")
wa, wb = symbols("wa wb", real=True)
W = wa + I * wb
assert expand((Z - W) * conjugate(Z - W) - (Z * conjugate(Z) - conjugate(W) * Z - W * conjugate(Z) + W * conjugate(W))) == 0
# Réciproque : z zb + bb z + b zb + g = |z + b|^2 - (|b|^2 - g).
ba, bb_ = symbols("ba bb", real=True)
g = symbols("g", real=True)
Bt = ba + I * bb_
assert expand(Z * conjugate(Z) + conjugate(Bt) * Z + Bt * conjugate(Z) + g - ((Z + Bt) * conjugate(Z + Bt) - (Bt * conjugate(Bt) - g))) == 0
# Trous 1-2 : leurres faux pour beta = i (centre) et pour gamma != 0, |beta|^2 - gamma not in {0,1} (rayon).
# (chaque leurre diffère de -beta pour un beta convenable : beta = 1 pour beta et conj(beta), beta = i pour -conj(beta))
assert 1 != -1 and complex(1).conjugate() != -1 and -(1j).conjugate() != -1j
assert abs(math.sqrt(5 - 1) - (5 - 1)) > 0.1 and abs(math.sqrt(5 - 1) - math.sqrt(5 + 1)) > 0.1
# Démonstration de la cocyclicité : (b-d)/(a-d) * (a-c)/(b-c) = birapport.
a, b, c, d = (complex(random.uniform(-3, 3), random.uniform(-3, 3)) for _ in range(4))
assert abs((b - d) / (a - d) * (a - c) / (b - c) - birapport(a, b, c, d)) < 1e-9
# Critère : birapport réel pour 4 points cocycliques aléatoires, non réel sinon.
for _ in range(100):
    o = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    r = random.uniform(0.5, 3)
    pts = [o + r * cmath.exp(1j * random.uniform(0, 2 * math.pi)) for _ in range(4)]
    assert abs(birapport(*pts).imag) < 1e-7 * max(1, abs(birapport(*pts)))
# Trou 3, leurre « de module 1 » : 1, i, -1, -i cocycliques, birapport 2.
assert abs(birapport(1, 1j, -1, -1j) - 2) < EPS
# Pièges 2, 3 et trou 3 leurre 3 : a=-1, b=1, c=2i, d=-2i : (c-a)/(c-b)*(d-a)/(d-b) = 1 réel, modules égaux, non cocycliques.
a, b, c, d = -1, 1, 2j, -2j
assert abs((c - a) / (c - b) * (d - a) / (d - b) - 1) < EPS
assert abs(abs((c - a) / (c - b)) - abs((d - a) / (d - b))) < EPS and abs(birapport(a, b, c, d).imag) > 0.1
# Piège 1 : un birapport imaginaire pur n'est pas réel (exemple : 1, -1, i, d avec birapport imaginaire pur).
a, b, c = 1, -1, 1j
k = (c - a) / (c - b)                    # = i
assert abs(k - 1j) < EPS
# d = 0 : (d+1)/(d-1) = -1, birapport = -i, imaginaire pur ; 0 n'est pas sur le cercle unité passant par 1, -1, i
d = 0
assert abs(birapport(a, b, c, d) + 1j) < EPS and abs(abs(d) - 1) > 0.5        # 0 n'est pas sur le cercle unité
# Exemple : z zb - (2-i) z - (2+i) zb + 1 = 0 : centre 2+i, rayon 2 ; birapport de 1, 3, 2+i, 2-i = -1.
for tt in (0.3, 1.7, 4.0):
    zt = 2 + 1j + 2 * cmath.exp(1j * tt)
    assert abs(zt * zt.conjugate() - (2 - 1j) * zt - (2 + 1j) * zt.conjugate() + 1) < 1e-9
assert abs(birapport(1, 3, 2 + 1j, 2 - 1j) + 1) < EPS
assert abs((1 + 1j) / (-1 + 1j) + 1j) < EPS and abs((-1 - 1j) / (1 - 1j) + 1j) < EPS
# Contre-exemples : z zb + z + 2 zb = 0 -> x^2+y^2+3x et -y ; birapport de 0,1,2,3 = 4/3.
e = expand(Z * conjugate(Z) + Z + 2 * conjugate(Z))
assert simplify(re(e) - (xr**2 + yr**2 + 3 * xr)) == 0 and simplify(im(e) + yr) == 0
assert birapport(F(0), F(1), F(2), F(3)) == F(4, 3)
# Affirmation 1 : |z-1| = 2|z+1| : cercle de centre -5/3 et de rayon 4/3.
for tt in (0.1, 2.0, 3.5):
    zt = -5 / 3 + 4 / 3 * cmath.exp(1j * tt)
    assert abs(abs(zt - 1) - 2 * abs(zt + 1)) < 1e-9
assert abs(math.sqrt(25 / 9 - 1) - 4 / 3) < EPS
# Affirmation 2 : 1, i, -1, 1+i non cocycliques.
assert abs(birapport(1, 1j, -1, 1 + 1j).imag) > 0.1

# geo-l2-plane-complexes-cercles-ex1 (+ v1, v2) : centre -beta, rayon sqrt(|beta|^2 - gamma).
for coef_z, coef_zb, gam, centre, r in [(2 - 4j, 2 + 4j, 11, -2 - 4j, 3), (-(3 + 1j), -(3 - 1j), 6, 3 - 1j, 2),
                                        (1 + 1j, 1 - 1j, -7, -1 + 1j, 3)]:
    bet = coef_zb
    assert abs(coef_z - bet.conjugate()) < EPS and abs(-bet - centre) < EPS
    assert abs(abs(bet) ** 2 - gam - r * r) < 1e-9
    for tt in (0.2, 1.9, 4.4):
        zt = centre + r * cmath.exp(1j * tt)
        assert abs(zt * zt.conjugate() + coef_z * zt + coef_zb * zt.conjugate() + gam) < 1e-9
    # Leurres : centre beta, rayon r^2, centre -conj(beta) : distincts de la bonne réponse.
    assert abs(bet - centre) > 0.1 and r * r != r and abs(-bet.conjugate() - centre) > 0.1

# geo-l2-plane-complexes-cercles-ex2 (+ v1, v2) : |z - a| = 2 |z - b|.
for a_, b_, centre, r, leurres in [(4, 1, 0, 2, [(1, 2), (0, 4)]), (2j, -1j, -2j, 2, [(2j, 2), (-1j, 2)]),
                                   (0, 3, 4, 2, [(-4, 2), (3, 2)])]:
    for tt in (0.4, 2.2, 5.0):
        zt = centre + r * cmath.exp(1j * tt)
        assert abs(abs(zt - a_) - 2 * abs(zt - b_)) < 1e-9
    for (cl, rl) in leurres:
        zt = cl + rl                                   # un point du cercle leurre qui ne convient pas
        assert abs(abs(zt - a_) - 2 * abs(zt - b_)) > 0.1 or abs(abs(cl + 1j * rl - a_) - 2 * abs(cl + 1j * rl - b_)) > 0.1
    # Leurre médiatrice : le milieu de [a, b] ne convient pas.
    mil = (a_ + b_) / 2
    assert abs(abs(mil - a_) - 2 * abs(mil - b_)) > 0.1

# geo-l2-plane-complexes-cercles-ex3 (+ v1, v2) : D_t = u + t i cocyclique avec u, -u, v i.
for u, v, t_ok, k, leurres in [(1, 2, F(3, 2), F(3, 4), [F(3, 4), F(-3, 2)]), (1, 3, F(8, 3), F(4, 3), [F(4, 3), F(-8, 3)]),
                               (2, 4, F(3), F(3, 2), [F(3, 2), F(-3)])]:
    a, b, c = complex(u), complex(-u), complex(0, v)
    o, r = cercle3(a, b, c)
    assert abs(o - complex(0, float(k))) < EPS
    assert abs(birapport(a, b, c, complex(u, float(t_ok))).imag) < 1e-9
    assert abs(abs(complex(u, float(t_ok)) - o) - r) < 1e-9
    for l in leurres:
        assert abs(birapport(a, b, c, complex(u, float(l))).imag) > 0.01
    # t = 0 donne D_t = A (points non distincts) : la proposition {0, t} est fausse.
    assert complex(u, 0) == a
# Formules de la correction : (-1+2i)/(1+2i) = (3+4i)/5 ; (-1+3i)/(1+3i) = (4+3i)/5.
assert abs((-1 + 2j) / (1 + 2j) - (3 + 4j) / 5) < EPS and abs((-1 + 3j) / (1 + 3j) - (4 + 3j) / 5) < EPS

# geo-l2-plane-complexes-cercles-ex4 (Ptolémée) : identité, inégalité, égalité pour un quadrilatère convexe inscrit direct.
A_, B_, C_, D_ = symbols("a b c d")
assert expand((C_ - A_) * (D_ - B_) - (B_ - A_) * (D_ - C_) - (D_ - A_) * (C_ - B_)) == 0
for _ in range(200):
    a, b, c, d = (complex(random.uniform(-3, 3), random.uniform(-3, 3)) for _ in range(4))
    assert abs(c - a) * abs(d - b) <= abs(b - a) * abs(d - c) + abs(d - a) * abs(c - b) + 1e-9
    ts = sorted(random.uniform(0, 2 * math.pi) for _ in range(4))
    o = complex(random.uniform(-2, 2), random.uniform(-2, 2))
    r = random.uniform(0.5, 3)
    a, b, c, d = (o + r * cmath.exp(1j * t) for t in ts)       # sens direct
    q = (b - a) * (d - c) / ((d - a) * (c - b))
    assert abs(q.imag) < 1e-7 * max(1, abs(q)) and q.real > 0
    assert abs(abs(c - a) * abs(d - b) - (abs(b - a) * abs(d - c) + abs(d - a) * abs(c - b))) < 1e-8
    assert abs(ang(b - a, d - a) - geo(b - a, d - a)) < 1e-9 and abs(ang(d - c, b - c) - geo(d - c, b - c)) < 1e-9

# ===========================================================================
# Carte geo-l2-plane-puissance
# ===========================================================================
for _ in range(200):
    o = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    R = random.uniform(0.5, 3)
    m = complex(random.uniform(-6, 6), random.uniform(-6, 6))
    u = cmath.exp(1j * random.uniform(0, math.pi))
    # Intersections de la droite m + s u avec le cercle : s^2 + 2 s Re(conj(u)(m-o)) + |m-o|^2 - R^2 = 0.
    pb = dot(u, m - o)
    disc = pb * pb - (abs(m - o) ** 2 - R * R)
    if disc < 0:
        continue
    s1, s2 = -pb + math.sqrt(disc), -pb - math.sqrt(disc)
    assert abs(dot(s1 * u, s2 * u) - (abs(m - o) ** 2 - R * R)) < 1e-8
    # Affirmation 1 : diamètre quelconque.
    e = cmath.exp(1j * random.uniform(0, 2 * math.pi))
    assert abs(dot(o + R * e - m, o - R * e - m) - (abs(m - o) ** 2 - R * R)) < 1e-8
# Pièges : M = O, A = B (points quelconques) ; centre remplacé par un point intérieur P ; B symétrique de A par rapport à M.
R = 2
assert dot(R, R) == R * R != -R * R
P = 0.5                     # point intérieur, corde passant par P : produit = OP^2 - R^2
assert abs((P**2 - R * R) - (0 - R * R)) > 0.1
mm, aa = 5, 2                # M extérieur, A sur le cercle, B = 2M - A
assert abs(-abs(aa - mm) ** 2 - (mm**2 - R * R)) > 0.1
# Trous : leurres faux (exemple M(13,0), cercle de centre O et de rayon 5 : p = 144 ; M = O pour MA x MB).
assert 13**2 - 25 == 144 and 13 - 5 != 144 and 25 - 169 != 144 and 169 + 25 != 144
assert 5 * 5 != -25 and (8 + 18) != 144
# Exemple : M(13,0), A(5,0), B(-5,0).
assert (5 - 13) * (-5 - 13) == 144 and math.isqrt(144) == 12
# Exemple : axe radical de x^2+y^2-1 et x^2+y^2-6x : 6x - 1 = 0, points (1/6, +/- sqrt35/6).
assert F(1, 36) + F(35, 36) == 1 and F(1, 36) + F(35, 36) - 6 * F(1, 6) == 0
# Contre-exemples : concentriques (différence 3) ; x^2+y^2=1 et (x-4)^2+y^2=9 : axe x = 1 (tangents en (1,0)).
X = symbols("X")
assert expand((X**2 - 1) - ((X - 4) ** 2 - 9)) == 8 * X - 8 and (1 - 4) ** 2 == 9
# Affirmation 3 : x^2+y^2=1 et (x-5)^2+y^2=1 : axe x = 5/2.
assert expand((X**2 - 1) - ((X - 5) ** 2 - 1)) == 10 * X - 25

# Axe radical (carte) : droite perpendiculaire à (OO'), point K + t OO'.
for _ in range(100):
    o1 = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    o2 = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    r1, r2 = random.uniform(0.5, 3), random.uniform(0.5, 3)
    kk = (o1 + o2) / 2
    t = (r1**2 - r2**2) / (2 * abs(o2 - o1) ** 2)
    p0 = kk + t * (o2 - o1)
    for s in (-2, 0.3, 5):
        mpt = p0 + s * 1j * (o2 - o1)
        assert abs((abs(mpt - o1) ** 2 - r1**2) - (abs(mpt - o2) ** 2 - r2**2)) < 1e-8

# geo-l2-plane-puissance-ex1 (+ v1, v2) : puissance et tangente ; leurres (p, p), (MO - R), (MO^2 + R^2).
for O_, R_, M_, p_, mt in [((1, 2), 3, (5, -1), 16, 4), ((-2, 1), 5, (4, 9), 75, 5 * sqrt(3)), ((1, 1), 2, (7, 9), 96, 4 * sqrt(6))]:
    mo2 = (M_[0] - O_[0]) ** 2 + (M_[1] - O_[1]) ** 2
    assert mo2 - R_**2 == p_ and simplify(mt - sqrt(p_)) == 0
    mo = sqrt(mo2)
    assert mo.is_integer and (mo - R_, mo2 + R_**2) in [(2, 34), (5, 125), (8, 104)]
    assert p_ != mt and mo - R_ != p_

# geo-l2-plane-puissance-ex2 (+ v1, v2) : MD = MA.MB / MC ; leurres MB.MC/MA, MA + MB - MC, MA.MC/MB.
for MA, MB, MC, MD, leurres in [(3, 12, 4, 9, [16, 11, 1]), (2, 6, 3, 4, [9, 5, 1]), (4, 10, 5, 8, [F(25, 2), 9, 2])]:
    assert F(MA * MB, MC) == MD
    assert leurres == [F(MB * MC, MA), MA + MB - MC, F(MA * MC, MB)] and MD not in leurres

# geo-l2-plane-puissance-ex3 (+ v1, v2) : centre radical.
def pcercle(cf, x_, y_):
    a_, b_, c_ = cf
    return x_ * x_ + y_ * y_ + 2 * a_ * x_ + 2 * b_ * y_ + c_


for cercles, K, leurres in [
    ([(-2, -1, -4), (1, -3, 6), (0, 0, -8)], (F(-1, 7), F(16, 7)), [(F(1, 3), F(4, 3)), (F(1, 7), F(16, 7)), (F(-1, 7), F(-16, 7))]),
    ([(-1, -2, -4), (2, 1, -4), (-3, 1, 1)], (F(1, 2), F(-1, 2)), [(F(2, 3), 0), (F(-1, 2), F(1, 2)), (F(1, 2), F(1, 2))]),
    ([(-1, -1, -7), (2, 1, 1), (3, -1, 9)], (-2, 1), [(F(-4, 3), F(1, 3)), (2, -1), (-2, -1)]),
]:
    ps = [pcercle(cf, *K) for cf in cercles]
    assert ps[0] == ps[1] == ps[2]
    rayons2 = [a_ * a_ + b_ * b_ - c_ for a_, b_, c_ in cercles]
    assert all(r2 > 0 for r2 in rayons2)
    (a1, b1, _), (a2, b2, _), (a3, b3, _) = cercles
    assert (a2 - a1) * (b3 - b1) - (a3 - a1) * (b2 - b1) != 0           # centres non alignés
    for L in leurres:
        q = [pcercle(cf, *L) for cf in cercles]
        assert not (q[0] == q[1] == q[2])
    # Leurre « isobarycentre des centres » (premier leurre).
    assert leurres[0] == (F(-(a1 + a2 + a3), 3), F(-(b1 + b2 + b3), 3))
# Rayons cités : 3, 2, 2 sqrt2 ; 3, 3, 3 ; 3, 2, 1.
assert [a_ * a_ + b_ * b_ - c_ for a_, b_, c_ in [(-2, -1, -4), (1, -3, 6), (0, 0, -8)]] == [9, 4, 8]
assert [a_ * a_ + b_ * b_ - c_ for a_, b_, c_ in [(-1, -2, -4), (2, 1, -4), (-3, 1, 1)]] == [9, 9, 9]
assert [a_ * a_ + b_ * b_ - c_ for a_, b_, c_ in [(-1, -1, -7), (2, 1, 1), (3, -1, 9)]] == [9, 4, 1]
# Équations des axes citées dans la correction principale.
assert pcercle((-2, -1, -4), X, 0) - pcercle((0, 0, -8), X, 0) == -4 * X + 4

# geo-l2-plane-puissance-ex4 : HA.HA1 = HB.HB1 = HC.HC1.
for _ in range(100):
    a, b, c = (complex(random.uniform(-5, 5), random.uniform(-5, 5)) for _ in range(3))
    if abs(((b - a).conjugate() * (c - a)).imag) < 0.5:
        continue
    h = orthocentre(a, b, c)
    pied = lambda p, u, v: u + (v - u) * dot(v - u, p - u) / abs(v - u) ** 2
    a1, b1, c1 = pied(a, b, c), pied(b, c, a), pied(c, a, b)
    v1, v2, v3 = dot(a - h, a1 - h), dot(b - h, b1 - h), dot(c - h, c1 - h)
    assert abs(v1 - v2) < 1e-6 and abs(v2 - v3) < 1e-6

# ===========================================================================
# Carte geo-l2-plane-orthogonaux
# ===========================================================================
def inter(o1, r1, o2, r2):
    d = abs(o2 - o1)
    a_ = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    h_ = math.sqrt(max(r1 * r1 - a_ * a_, 0))
    u = (o2 - o1) / d
    return o1 + a_ * u + 1j * h_ * u, o1 + a_ * u - 1j * h_ * u


for _ in range(100):
    o1 = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    r1, r2 = random.uniform(0.5, 3), random.uniform(0.5, 3)
    o2 = o1 + math.sqrt(r1 * r1 + r2 * r2) * cmath.exp(1j * random.uniform(0, 2 * math.pi))
    for A_ in inter(o1, r1, o2, r2):
        assert abs(abs(A_ - o1) - r1) < 1e-8 and abs(abs(A_ - o2) - r2) < 1e-8
        assert abs(dot(A_ - o1, A_ - o2)) < 1e-8
    # Droite passant par O1 coupant C2 en M, M' : OM . OM' = r1^2.
    u = cmath.exp(1j * random.uniform(0, math.pi))
    pb = dot(u, o1 - o2)
    disc = pb * pb - (abs(o1 - o2) ** 2 - r2 * r2)
    if disc >= 0:
        s1, s2 = -pb + math.sqrt(disc), -pb - math.sqrt(disc)
        assert abs(s1 * s2 - r1 * r1) < 1e-8
# Pièges : OO' = R + R' (tangents) ; R=5, R'=3, d=4 (angle droit en O', pas en A) ; d^2 = 2R'^2 avec R=1, R'=2.
for (r1, r2, d) in [(5, 3, 4), (1, 2, math.sqrt(8))]:
    A_, _ = inter(0, r1, d, r2)
    assert abs(abs(A_) - r1) < 1e-9 and abs(abs(A_ - d) - r2) < 1e-9 and abs(dot(A_, A_ - d)) > 0.1
assert 3 + 1 == 4 and (3 + 1) ** 2 != 9 + 1
# Trous : leurres faux sur l'exemple R = 1, R' = sqrt3, OO' = 2.
assert 4 == 1 + 3 and 2 != 1 + math.sqrt(3) and 4 != 1 - 3 and 4 != math.sqrt(3)
assert 4 - 3 == 1 and 4 - 3 != 3 and 4 - 3 != 0 and 4 - 3 != -1
# Exemple : A(1/2, sqrt3/2), AO . AO' = 0.
A_ = complex(0.5, math.sqrt(3) / 2)
assert abs(abs(A_) - 1) < EPS and abs(abs(A_ - 2) - math.sqrt(3)) < EPS and abs(dot(-A_, 2 - A_)) < EPS
# Contre-exemples : x^2+y^2=1 et (x-1)^2+y^2=1 : rayons à pi/3 ; non-transitivité.
A_ = complex(0.5, math.sqrt(3) / 2)
assert abs(geo(0 - A_, 1 - A_) - math.pi / 3) < EPS
assert 2**2 == 1 + 3 and 4**2 != 3 + 3
# Les trous 3 : à un point commun, O, A, O' ne sont pas alignés, A n'est pas le milieu, (OA) non perpendiculaire à (OO').
A_ = complex(0.5, math.sqrt(3) / 2)
assert abs((A_.conjugate() * 2).imag) > 0.1 and abs(A_ - 1) > 0.1 and abs(dot(A_, 2)) > 0.1

# geo-l2-plane-orthogonaux-ex1 (+ v1, v2) : R'^2 = OO'^2 - R^2 ; leurres OO' - R, sqrt(OO'^2 + R^2), R'^2.
for O_, R_, O2, Rp, leurres in [((0, 0), 6, (6, 8), 8, [4, 2 * sqrt(34), 64]), ((1, -1), 2, (4, 3), sqrt(21), [3, sqrt(29), 21]),
                                ((0, 0), 5, (-9, 12), 10 * sqrt(2), [10, 5 * sqrt(10), 200])]:
    d2 = (O2[0] - O_[0]) ** 2 + (O2[1] - O_[1]) ** 2
    assert simplify(Rp**2 - (d2 - R_**2)) == 0
    assert simplify(leurres[0] - (sqrt(d2) - R_)) == 0 and simplify(leurres[1] - sqrt(d2 + R_**2)) == 0
    assert simplify(leurres[2] - Rp**2) == 0 and all(simplify(l - Rp) != 0 for l in leurres)

# geo-l2-plane-orthogonaux-ex2 (+ v1, v2) : 2aa' + 2bb' = c + c' ; leurres (oubli du 2, signe, R^2 - R'^2).
for (a1, b1, c1), (a2, b2), m, leurres in [((-2, 1, -4), (1, -3), -6, [14, -1, 26]), ((-1, 0, -3), (2, -3), -1, [7, 1, 27]),
                                           ((3, -1, 6), (-1, 2), -16, [4, -11, 26])]:
    assert 2 * a1 * a2 + 2 * b1 * b2 == c1 + m
    d2 = (a1 - a2) ** 2 + (b1 - b2) ** 2
    R2, Rp2 = a1 * a1 + b1 * b1 - c1, a2 * a2 + b2 * b2 - m
    assert R2 > 0 and Rp2 > 0 and d2 == R2 + Rp2
    sgn = -(2 * a1 * a2 + 2 * b1 * b2) - c1           # 2aa' + 2bb' = -(c + c')
    oubli = a1 * a2 + b1 * b2 - c1                     # aa' + bb' = c + c'
    diff = d2 - R2 + (a2 * a2 + b2 * b2)               # d^2 = R^2 - R'^2
    assert sorted(leurres) == sorted([sgn, oubli, diff]) and m not in leurres

# geo-l2-plane-orthogonaux-ex3 (+ v1, v2) : centre sur l'axe radical x = x0, rho^2 = p(Omega).
for R1sq, c2, R2sq, om, rho2, leurres2 in [(1, 4, 1, (2, 1), 4, [5, 16, 6]), (4, 6, 28, (1, 3), 6, [10, 36, 14]),
                                           (4, 8, 20, (3, 2), 9, [13, 81, 17])]:
    x0 = F(c2 * c2 + R1sq - R2sq, 2 * c2)
    assert x0 == om[0]
    p1 = om[0] ** 2 + om[1] ** 2 - R1sq
    p2 = (om[0] - c2) ** 2 + om[1] ** 2 - R2sq
    assert p1 == p2 == rho2
    # Leurres (carrés) : O Omega^2, rho^4 (rho^2 pris pour rho), O Omega^2 + R1^2.
    oo = om[0] ** 2 + om[1] ** 2
    assert leurres2 == [oo, rho2 * rho2, oo + R1sq] and rho2 not in leurres2

# geo-l2-plane-orthogonaux-ex4 : cercles passant par M et son inverse M'.
for _ in range(100):
    o = complex(random.uniform(-2, 2), random.uniform(-2, 2))
    R = random.uniform(0.5, 2)
    m = o + random.uniform(0.2, 4) * cmath.exp(1j * random.uniform(0, 2 * math.pi))
    if abs(abs(m - o) - R) < 0.05:
        continue
    mp = o + (m - o) * R * R / abs(m - o) ** 2
    n = complex(random.uniform(-4, 4), random.uniform(-4, 4))
    if abs(((m - o).conjugate() * (n - o)).imag) < 0.3:
        continue
    om_, rho = cercle3(m, mp, n)
    assert abs(abs(om_ - o) ** 2 - (R * R + rho * rho)) < 1e-6 * max(1, abs(om_ - o) ** 2)

# ===========================================================================
# QCM geo-l2-plane
# ===========================================================================
# q1 : AOB = 2 x 35.
assert 2 * 35 == 70
# q3 : cosinus du plus grand angle du triangle 4, 5, 6 ; leurres = cosinus des autres angles.
assert F(16 + 25 - 36, 40) == F(1, 8) and F(25 + 36 - 16, 60) == F(3, 4) and F(16 + 36 - 25, 48) == F(9, 16)
# q4 : triangle 2, 3, 4 : cos = -1/4.
assert F(4 + 9 - 16, 12) == F(-1, 4)
# q5 : (z - i)/(z + 1) imaginaire pur sur le cercle de diamètre [i, -1] (privé de -1).
ctr, rr = (1j - 1) / 2, abs(1j + 1) / 2
for tt in (0.3, 1.0, 2.5, 4.0, 5.5):
    zt = ctr + rr * cmath.exp(1j * tt)
    if abs(zt + 1) > 1e-6:
        assert abs(((zt - 1j) / (zt + 1)).real) < 1e-9
zt = complex(-0.5, 0.5)            # milieu de [AB], sur la droite (AB) : quotient réel, non imaginaire pur
assert abs(((zt - 1j) / (zt + 1)).imag) < 1e-9 and ((zt - 1j) / (zt + 1)).real < 0
# q6 : OM = 3, R = 5 : MA x MB = 16.
assert 25 - 9 == 16
# q7 : x^2+y^2=4 et (x-3)^2+(y-4)^2=21.
assert 3**2 + 4**2 == 4 + 21
# q8 : unité et centre (2,0), rayon sqrt3 : O extérieur.
assert 2 > math.sqrt(3) and 4 == 1 + 3

# ===========================================================================
# Compléments du relecteur
# ===========================================================================
# geo-l2-plane-orthogonaux (cours) : pour x^2+y^2+2ax+2by+c=0 et x^2+y^2+2a'x+2b'y+c'=0,
# OO'^2 - R^2 - R'^2 = c + c' - 2aa' - 2bb' (identité symbolique).
a_, b_, c_, a2_, b2_, c2_ = symbols("a_ b_ c_ a2_ b2_ c2_", real=True)
assert expand((a_ - a2_)**2 + (b_ - b2_)**2 - (a_**2 + b_**2 - c_) - (a2_**2 + b2_**2 - c2_)
              - (c_ + c2_ - 2*a_*a2_ - 2*b_*b2_)) == 0
# geo-l2-plane-remarquables (démonstration, bissectrices) : le pied D de la bissectrice issue de A,
# barycentre de (B, b) et (C, c), vérifie BD/DC = c/b et voit [AB], [AC] sous des angles égaux.
for _ in range(50):
    pa, pb, pc = (complex(random.uniform(-3, 3), random.uniform(-3, 3)) for _ in range(3))
    if abs(((pb - pa).conjugate() * (pc - pa)).imag) < 0.5:
        continue
    la, lb, lc = abs(pc - pb), abs(pc - pa), abs(pb - pa)
    pd = (lb * pb + lc * pc) / (lb + lc)
    assert abs(abs(pd - pb) / abs(pc - pd) - lc / lb) < 1e-9
    assert abs(abs(cmath.phase((pd - pa) / (pb - pa))) - abs(cmath.phase((pc - pa) / (pd - pa)))) < 1e-9
# geo-l2-plane-complexes-cercles-ex4 (question 3) : pour ABCD convexe inscrit direct, q = (b-a)(d-c)/((d-a)(c-b)) > 0.
for _ in range(50):
    ts = sorted(random.uniform(0, 2 * math.pi) for _ in range(4))
    if min(ts[i + 1] - ts[i] for i in range(3)) < 0.1 or 2 * math.pi - ts[3] + ts[0] < 0.1:
        continue
    pa, pb, pc, pd = (cmath.exp(1j * u) for u in ts)
    q = (pb - pa) * (pd - pc) / ((pd - pa) * (pc - pb))
    assert abs(q.imag) < 1e-9 and q.real > 0

print("geo_l2_plane : toutes les vérifications passent")
