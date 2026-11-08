"""Vérifications calculatoires du chapitre « geo-l3-euclidienne » (géométrie, L3) : cartes, exercices et QCM de
content/*/geometrie/l3-geo-l3-euclidienne.json et content/qcm/geo-l3-euclidienne.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question. Les
applications affines sont des couples (M, c) : X -> M X + c, en calcul exact (SymPy) ; les isométries du plan sont aussi
écrites en complexes. Les groupes d'isométries des figures finies sont obtenus par force brute : on énumère les
permutations des sommets qui conservent les distances, et l'on construit pour chacune l'application linéaire (dans un
repère centré à l'isobarycentre) qui la réalise, dont on contrôle qu'elle est orthogonale et qu'elle conserve la figure.
Lancer : python3 verification/geo_l3_euclidienne.py.
"""

from itertools import permutations

from sympy import (I, Matrix, expand_complex, Rational, cos, eye, conjugate, expand, nsimplify, pi, simplify, sin, sqrt, symbols,
                   zeros, re, im, exp)

R = Rational
x, y, z, t, a_, b_ = symbols("x y z t a b", real=True)
th = symbols("theta", real=True)


# ---------------------------------------------------------------------------
# Outils : applications affines (M, c), décomposition canonique
# ---------------------------------------------------------------------------
def aff(exprs, vars_):
    """Application affine donnée par ses expressions : renvoie (M, c)."""
    n = len(vars_)
    M = Matrix(n, n, lambda i, j: exprs[i].coeff(vars_[j]))
    c = Matrix([exprs[i].subs({v: 0 for v in vars_}) for i in range(n)])
    assert all(expand(exprs[i] - (M.row(i) * Matrix(vars_))[0] - c[i]) == 0 for i in range(n))
    return M, c


def appl(f, X):
    M, c = f
    return M * Matrix(X) + c


def compo(f, g):
    """f ∘ g."""
    return (f[0] * g[0], f[0] * g[1] + f[1])


def trans(u):
    u = Matrix(u)
    return (eye(len(u)), u)


def egal(f, g):
    return simplify(f[0] - g[0]) == zeros(*f[0].shape) and simplify(f[1] - g[1]) == zeros(*f[1].shape)


def orthogonale(M):
    return simplify(M.T * M - eye(M.shape[0])) == zeros(*M.shape)


def proj_ortho(v, base):
    """Projeté orthogonal de v sur Vect(base) (base quelconque d'un sous-espace)."""
    if not base:
        return zeros(len(v), 1)
    B = Matrix.hstack(*base)
    return simplify(B * (B.T * B).inv() * B.T * Matrix(v))


def canonique(f):
    """Décomposition canonique de l'isométrie f = (M, c) : renvoie (u, g, base de ker(M - I), points fixes de g)."""
    M, c = f
    assert orthogonale(M)
    n = M.shape[0]
    noyau = (M - eye(n)).nullspace()
    u = proj_ortho(c, noyau)
    g = compo(trans(-u), f)
    assert simplify(M * u - u) == zeros(n, 1)  # u fixé par la partie linéaire
    return u, g, noyau


def points_fixes(f):
    """Solutions de M X + c = X : None si aucune, sinon (point particulier, base de la direction)."""
    M, c = f
    n = M.shape[0]
    A = M - eye(n)
    from sympy import linsolve
    X = symbols(f"p0:{n}")
    sol = linsolve((A, -c), *X)
    if not sol:
        return None
    (s,) = list(sol)
    s = Matrix(s)
    libres = sorted(s.free_symbols, key=str)
    part = s.subs({v: 0 for v in libres})
    return part, (A).nullspace()


def sur_droite_plan(P, a, b, cst):
    """P = (x, y) vérifie a x + b y = cst."""
    return simplify(a * P[0] + b * P[1] - cst) == 0


def reflexion_droite(a, b, cst):
    """Réflexion du plan d'axe a x + b y = cst."""
    nn = a * a + b * b
    exprs = [x - 2 * a * (a * x + b * y - cst) / nn, y - 2 * b * (a * x + b * y - cst) / nn]
    return aff(exprs, [x, y])


def reflexion_plan(n, cst):
    """Réflexion de l'espace par rapport au plan <n, X> = cst."""
    n = Matrix(n)
    X = Matrix([x, y, z])
    e = X - 2 * (n.dot(X) - cst) / n.dot(n) * n
    return aff(list(e), [x, y, z])


def demi_tour(P, b):
    """Demi-tour d'axe la droite passant par P dirigée par b."""
    P, b = Matrix(P), Matrix(b)
    b = b / sqrt(b.dot(b))
    X = Matrix([x, y, z])
    e = P + 2 * (X - P).dot(b) * b - (X - P)
    return aff([simplify(v) for v in e], [x, y, z])


def rotation_axe(a, theta):
    """Matrice de la rotation d'axe orienté par a, d'angle theta (formule de Rodrigues)."""
    a = Matrix(a) / sqrt(Matrix(a).dot(Matrix(a)))
    A = Matrix([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return simplify(cos(theta) * eye(3) + (1 - cos(theta)) * a * a.T + sin(theta) * A)


def angle_signe(M, a):
    """Pour M dans SO(3) ou O(3)^-, d'axe a : renvoie (cos θ, signe de sin θ) avec l'axe orienté par a."""
    a = Matrix(a)
    d = M.det()
    cth = (M.trace() - d) / 2  # tr = d + 2 cos θ
    # vecteur orthogonal à a
    for v in (Matrix([1, -1, 0]), Matrix([0, 1, -1]), Matrix([1, 0, 0]), Matrix([0, 1, 0])):
        w = v - a.dot(v) / a.dot(a) * a
        if w != zeros(3, 1):
            break
    s = Matrix.hstack(a, w, M * w).det()
    return simplify(cth), (0 if s == 0 else (1 if s > 0 else -1))


def cz(e):
    return simplify(expand_complex(e))


def nul(e):
    return simplify(expand_complex(e)) == 0


Z = symbols("Z")  # variable complexe


def compl(fun):
    """Isométrie / similitude complexe -> fonction sympy."""
    return fun


def conj(w):
    return conjugate(w)


# ===========================================================================
# Carte decomposition
# ===========================================================================
# geo-l3-euclidienne-decomposition (exemple) : f(x,y) = (x+2, -y+1), u = (2,0), g réflexion d'axe y = 1/2
f = aff([x + 2, -y + 1], [x, y])
u, g, ker = canonique(f)
assert u == Matrix([2, 0])
assert egal(g, reflexion_droite(0, 1, R(1, 2)))
s0 = reflexion_droite(0, 1, 0)
assert egal(compo(trans([2, 1]), s0), f)  # écriture exacte mais non canonique
assert not egal(compo(trans([2, 1]), s0), compo(s0, trans([2, 1])))  # pas de commutation
assert points_fixes(f) is None

# geo-l3-euclidienne-decomposition (contre-exemple) : cisaillement (x+y, y+1) : ker = im, pas de décomposition
f = aff([x + y, y + 1], [x, y])
A = f[0] - eye(2)
assert A.nullspace() == [Matrix([1, 0])] and A.columnspace() == [Matrix([1, 0])]
assert points_fixes(compo(trans([-a_, 0]), f)) is None  # pour tout a
# sans la commutation : la réflexion y=0 vaut t_(0,2) ∘ (réflexion d'axe y = -1)
assert egal(compo(trans([0, 2]), reflexion_droite(0, 1, -1)), reflexion_droite(0, 1, 0))
# affirmation : réflexion d'axe x=1 après celle d'axe x=0 = translation (2,0)
assert egal(compo(reflexion_droite(1, 0, 1), reflexion_droite(1, 0, 0)), trans([2, 0]))
# affirmation : rotation d'angle pi/2 de centre O : u = 0 mais M f(M) != 0 hors de O
rot = aff([-y, x], [x, y])
assert canonique(rot)[0] == zeros(2, 1) and appl(rot, [1, 0]) != Matrix([1, 0])
# trou 1 : pour f = id, ker = plan, im = {0} : ni égaux, ni de même dimension ; jamais tous deux nuls (somme = n)
assert len((eye(2) - eye(2)).nullspace()) == 2 and (eye(2) - eye(2)).rank() == 0
# trou 2 (« exactement un point fixe ») : pour une réflexion, ker(M - I) est une droite, donc g fixe une droite
assert len((reflexion_droite(0, 1, 0)[0] - eye(2)).nullspace()) == 1
# trou 3 (« quelconque », « orthogonal à ker », « dans Im ») : pour une rotation, t_{-v} ∘ f a un point fixe pour tout v
assert points_fixes(compo(trans([-a_, -b_]), rot)) is not None
# trou 4 : rotation (point fixe, M != I), réflexion (points fixes, 1 valeur propre), identité (u = 0)
assert points_fixes(rot) is not None and rot[0] != eye(2)
# piège 3 : pour la symétrie glissée (x+2, -y+1), aucun u = (0, b) orthogonal à ker ne donne t_{-u} ∘ f avec point fixe
assert points_fixes(compo(trans([0, -b_]), aff([x + 2, -y + 1], [x, y]))) is None
# lemme : ker(φ - id) ⊥ im(φ - id) pour une rotation d'axe (1,1,1) d'angle 2π/3
Mc = Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
for k_ in (Mc - eye(3)).nullspace():
    for i_ in (Mc - eye(3)).columnspace():
        assert k_.dot(i_) == 0


def verif_deco_plan(fexprs, bonne, leurres):
    """bonne / leurres : (u, (a, b, cst) axe de la réflexion g) ; bonne doit être canonique, les leurres non."""
    f = aff(fexprs, [x, y])
    u, g, _ = canonique(f)
    ub, (aa, bb, cc) = bonne
    assert u == Matrix(ub) and egal(g, reflexion_droite(aa, bb, cc))
    for ul, (aa, bb, cc) in leurres:
        cand = compo(trans(ul), reflexion_droite(aa, bb, cc))
        canon = egal(cand, f) and f[0] * Matrix(ul) == Matrix(ul)
        assert not canon


# geo-l3-euclidienne-decomposition-ex1 (et v1, v2)
verif_deco_plan([x + 3, -y + 2], ((3, 0), (0, 1, 1)), [((3, 2), (0, 1, 0)), ((3, 1), (0, 1, R(1, 2))), ((3, 0), (0, 1, 2))])
assert egal(compo(trans([3, 2]), reflexion_droite(0, 1, 0)), aff([x + 3, -y + 2], [x, y]))
assert egal(compo(trans([3, 1]), reflexion_droite(0, 1, R(1, 2))), aff([x + 3, -y + 2], [x, y]))
assert egal(compo(trans([3, 0]), reflexion_droite(0, 1, 2)), aff([x + 3, 4 - y], [x, y]))
verif_deco_plan([x - 4, -y + 6], ((-4, 0), (0, 1, 3)), [((-4, 6), (0, 1, 0)), ((-4, 2), (0, 1, 2)), ((-4, 0), (0, 1, 6))])
assert egal(compo(trans([-4, 2]), reflexion_droite(0, 1, 2)), aff([x - 4, -y + 6], [x, y]))
assert egal(compo(trans([-4, 0]), reflexion_droite(0, 1, 6)), aff([x - 4, 12 - y], [x, y]))
verif_deco_plan([-x + 2, y + 5], ((0, 5), (1, 0, 1)), [((2, 5), (1, 0, 0)), ((1, 5), (1, 0, R(1, 2))), ((0, 5), (1, 0, 2))])
assert egal(compo(trans([1, 5]), reflexion_droite(1, 0, R(1, 2))), aff([-x + 2, y + 5], [x, y]))
assert egal(compo(trans([0, 5]), reflexion_droite(1, 0, 2)), aff([4 - x, y + 5], [x, y]))


def verif_deco_espace(fexprs, ub, point, direction, gzero):
    f = aff(fexprs, [x, y, z])
    u, g, ker = canonique(f)
    assert u == Matrix(ub)
    pf = points_fixes(g)
    assert pf is not None and len(pf[1]) == 1  # une droite de points fixes
    assert appl(g, point) == Matrix(point) and pf[1][0].cross(Matrix(direction)) == zeros(3, 1)
    assert f[0] * Matrix(f[1]) != Matrix(f[1])  # le vecteur Of(O) n'est pas fixé par la partie linéaire
    assert appl(g, gzero) == zeros(3, 1)  # point de la « droite leurre » envoyé sur O
    assert points_fixes(f) is None
    # partie linéaire : demi-tour (rotation d'angle π)
    assert f[0].det() == 1 and f[0].trace() == -1


# geo-l3-euclidienne-decomposition-ex2 (et v1, v2) : vissages d'angle π
verif_deco_espace([-x + 2, -y + 4, z + 3], (0, 0, 3), (1, 2, 0), (0, 0, 1), (2, 4, 0))
verif_deco_espace([-x + 4, -y - 2, z - 1], (0, 0, -1), (2, -1, 0), (0, 0, 1), (4, -2, 0))
verif_deco_espace([x + 5, -y + 2, -z + 6], (5, 0, 0), (0, 1, 3), (1, 0, 0), (0, 2, 6))


def verif_glisse(fexprs, ub, axe, ff, leurres):
    """axe = (a, b, cst) ; leurres = liste (u, axe) ; ff = vecteur de f∘f."""
    f = aff(fexprs, [x, y])
    assert f[0].det() == -1 and points_fixes(f) is None
    u, g, _ = canonique(f)
    assert u == Matrix(ub) and egal(g, reflexion_droite(*axe))
    assert egal(compo(f, f), trans(ff)) and Matrix(ff) == 2 * u
    for ul, ax in leurres:
        assert not (Matrix(ul) == u and egal(reflexion_droite(*ax), g))


# geo-l3-euclidienne-decomposition-ex3 (et v1, v2)
verif_glisse([y + 1, x + 3], (2, 2), (-1, 1, 1), (4, 4), [((1, 3), (-1, 1, 0)), ((2, 2), (-1, 1, 2)), ((1, 3), (-1, 1, 1))])
assert Matrix([1, 3]) == 2 * Matrix([1, 1]) - Matrix([1, -1])
verif_glisse([y + 3, x - 1], (1, 1), (-1, 1, -2), (2, 2), [((3, -1), (-1, 1, 0)), ((1, 1), (-1, 1, 2)), ((3, -1), (-1, 1, -2))])
assert Matrix([3, -1]) == Matrix([1, 1]) + 2 * Matrix([1, -1])
verif_glisse([-y + 2, -x + 4], (-1, 1), (1, 1, 3), (-2, 2), [((2, 4), (1, 1, 0)), ((-1, 1), (1, 1, 2)), ((2, 4), (1, 1, 3))])
assert Matrix([2, 4]) == -Matrix([1, -1]) + 3 * Matrix([1, 1])

# geo-l3-euclidienne-decomposition-ex4 : ||Mf(M)||² = ||u||² + ||Mg(M)||² et f^k = t_{ku} ∘ g^k (exemple : vissage)
f = aff([-y + 1, x, z + 2], [x, y, z])
u, g, _ = canonique(f)
X = Matrix([x, y, z])
d1 = appl(f, X) - X
d2 = appl(g, X) - X
assert expand(d1.dot(d1) - u.dot(u) - d2.dot(d2)) == 0
fk, gk = f, g
for k in range(2, 6):
    fk, gk = compo(fk, f), compo(gk, g)
    assert egal(fk, compo(trans(k * u), gk)) and points_fixes(fk) is None

# geo-l3-euclidienne-decomposition-ex5 : h ∘ t_v ∘ h^{-1} = t_{h(v)} (h = réflexion d'axe x + 2y = 3)
h = reflexion_droite(1, 2, 3)
hinv = (h[0].inv(), -h[0].inv() * h[1])
assert egal(compo(compo(h, trans([a_, b_])), hinv), trans(h[0] * Matrix([a_, b_])))
# Q4 : la réflexion vectorielle s_{v-w} envoie v sur w quand ||v|| = ||w|| (v = (3,4), w = (5,0))
v_, w_ = Matrix([3, 4]), Matrix([5, 0])
n_ = v_ - w_
assert v_ - 2 * v_.dot(n_) / n_.dot(n_) * n_ == w_

# ===========================================================================
# Carte plan (complexes)
# ===========================================================================
a, b, w = symbols("a b w")
# geo-l3-euclidienne-plan : f∘f = translation de a b̄ + b quand |a| = 1 (a = e^{iθ})
aa = exp(I * th)
fz = lambda u: aa * conjugate(u) + b  # noqa: E731
assert nul(expand(fz(fz(w)) - w - (aa * conjugate(b) + b)))
# matrice de z -> a z̄ : (c+is)(x-iy) = cx+sy + i(sx-cy)
c_, s_ = symbols("c s", real=True)
assert expand((c_ + I * s_) * (x - I * y)) == expand(c_ * x + s_ * y + I * (s_ * x - c_ * y))


def nature_anti(av, bv):
    """Pour f(z) = a z̄ + b, |a| = 1 : renvoie (vecteur u, fonction test d'appartenance à l'axe)."""
    assert nul(abs(av) - 1)
    uu = cz((av * conjugate(bv) + bv) / 2)
    return uu


def axe_plan(av, bv, uu, P):
    """P (complexe) est-il sur l'axe : f(P) - P = u ?"""
    return cz(av * conjugate(P) + bv - P - uu) == 0


# geo-l3-euclidienne-plan (exemple) : i z̄ + 3 + i : glissée de vecteur 2+2i, axe y = x - 1 ; iz + 2 : centre 1+i
uu = nature_anti(I, 3 + I)
assert uu == 2 + 2 * I and axe_plan(I, 3 + I, uu, 1) and axe_plan(I, 3 + I, uu, 2 + I) and not axe_plan(I, 3 + I, uu, 0)
assert cz(2 / (1 - I)) == 1 + I
# exemple du cours : -z̄ + 2 est la réflexion d'axe x = 1
assert nature_anti(-1, 2) == 0 and cz(-conjugate(1 + 5 * I) + 2) == 1 + 5 * I
# contre-exemples : (2x, y/2) conserve les aires, n'est pas une isométrie ; (-y, x, z+1) sans point fixe
M2 = Matrix([[2, 0], [0, R(1, 2)]])
assert M2.det() == 1 and not orthogonale(M2)
assert points_fixes(aff([-y, x, z + 1], [x, y, z])) is None
# affirmation : f1 = z̄ + 1, f2 = -z̄ + i, f2∘f1 = -z - 1 + i, centre (-1+i)/2
assert nature_anti(1, 1) == 1 and nature_anti(-1, I) == I
assert cz(-conjugate(conjugate(w) + 1) + I) == cz(-w - 1 + I) and cz((-1 + I) / (1 - (-1))) == (-1 + I) / 2
# trou 3 : -z̄ + 2 réflexion avec b != 0 ; z̄ réflexion avec |b| != 1 ; -z̄ réflexion avec a != 1
assert nature_anti(-1, 2) == 0 and nature_anti(1, 0) == 0 and nature_anti(-1, 0) == 0


def verif_anti_reflexion(av, bv, sur_axe, hors_axe):
    assert nature_anti(av, bv) == 0
    for P in sur_axe:
        assert cz(av * conjugate(P) + bv) == cz(P)
    for P in hors_axe:
        assert cz(av * conjugate(P) + bv) != cz(P)


# geo-l3-euclidienne-plan-ex1 : -z̄ + 4 réflexion d'axe x = 2 ; leurres : x = 4 non fixe (f(4) = 0) ; -z + 4 != f
verif_anti_reflexion(-1, 4, [2, 2 + 3 * I], [4, 0])
assert cz(-conjugate(4) + 4) == 0 and cz(-conjugate(I) + 4) != cz(-I + 4)
# v1 : z̄ + 2i réflexion d'axe y = 1 ; f(2i) = 0
verif_anti_reflexion(1, 2 * I, [I, 3 + I], [2 * I, 0])
assert cz(conjugate(2 * I) + 2 * I) == 0
# v2 : i z̄ + 1 - i réflexion d'axe y = x - 1 ; f(0) = 1 - i
verif_anti_reflexion(I, 1 - I, [1, 3 + 2 * I], [0])
assert cz(I * conjugate(I * conjugate(w) + 1 - I) + 1 - I) == cz(w)

# geo-l3-euclidienne-plan-ex2 : rotations composées
r1 = lambda u: I * u  # noqa: E731
r2 = lambda u: 2 + I * (u - 2)  # noqa: E731
assert cz(r2(r1(w))) == cz(-w + 2 - 2 * I) and cz((2 - 2 * I) / 2) == 1 - I
assert cz(-(1 + I) + 2 - 2 * I) != 1 + I and cz(-1 + 2 - 2 * I) != 1  # leurres non fixes
# v1 : angle -π/2 : translation 2 + 2i
assert cz(2 - I * (I * w - 2)) == cz(w + 2 + 2 * I)
# v2 : angles π/3, centres 0 et 1 : rotation 2π/3, centre 1/2 - (√3/6) i
e3 = exp(I * pi / 3)
comp = cz(1 + e3 * (e3 * w - 1))
assert nul(comp - (exp(2 * I * pi / 3) * w + 1 - e3))
om = R(1, 2) - sqrt(3) / 6 * I
assert nul(exp(2 * I * pi / 3) * om + 1 - e3 - om)
assert nul((0 + 1 + (R(1, 2) - sqrt(3) / 2 * I)) / 3 - om)
assert not nul(exp(2 * I * pi / 3) * conjugate(om) + 1 - e3 - conjugate(om))
assert not nul(exp(I * pi / 3) * om + 1 - e3 - om)  # angle π/3 : faux

# geo-l3-euclidienne-plan-ex3 : f = r ∘ s
assert cz(I + I * (w - I)) == cz(I * w + 1 + I)
uu = nature_anti(I, 1 + I)
assert uu == 1 + I and axe_plan(I, 1 + I, uu, 0) and axe_plan(I, 1 + I, uu, 1 + I) and not axe_plan(I, 1 + I, uu, I)
# v1 : r de centre 1+i d'angle π : f = -z̄ + 2 + 2i, glissée de vecteur 2i, axe x = 1
assert cz((1 + I) - (w - (1 + I))) == cz(-w + 2 + 2 * I)
uu = nature_anti(-1, 2 + 2 * I)
assert uu == 2 * I and axe_plan(-1, 2 + 2 * I, uu, 1) and axe_plan(-1, 2 + 2 * I, uu, 1 + 7 * I) and not axe_plan(-1, 2 + 2 * I, uu, 2)
# v2 : r de centre 2 d'angle π/2 : f = i z̄ + 2 - 2i, réflexion d'axe y = x - 2
assert cz(2 + I * (w - 2)) == cz(I * w + 2 - 2 * I)
verif_anti_reflexion(I, 2 - 2 * I, [2, -2 * I, 5 + 3 * I], [0, 2 * I])

# geo-l3-euclidienne-plan-ex4 : s_D(z) = p + e^{2iβ} conj(z - p) ; produits de réflexions ; triangle
p1, p2, t1, s1 = symbols("p1 p2 t1 s1")
be, be1, be2 = symbols("beta beta1 beta2", real=True)
sD = lambda u, p, bb: p + exp(2 * I * bb) * conjugate(u - p)  # noqa: E731
tt = symbols("tt", real=True)
assert nul(sD(p1 + tt * exp(I * be), p1, be) - (p1 + tt * exp(I * be)))
om_ = symbols("omega")
assert nul(sD(sD(w, om_, be1), om_, be2) - (om_ + exp(2 * I * (be2 - be1)) * (w - om_)))
tr, sr = symbols("tr sr", real=True)
q2 = p1 + (tr + I * sr) * exp(I * be)
assert nul(sD(sD(w, p1, be), q2, be) - (w + 2 * I * sr * exp(I * be)))
# Q4 : r_A ∘ r_B ∘ r_C = id pour des triangles directs (angles calculés numériquement)
import cmath  # noqa: E402

for (A_, B_, C_) in [(0, 4, 1 + 3j), (0, 1, 1j), (-1 + 0.5j, 3 - 1j, 2 + 2j)]:
    assert ((B_ - A_).conjugate() * (C_ - A_)).imag > 0  # sens direct
    ang = lambda P, Q, S: abs(cmath.phase((S - P) / (Q - P)))  # noqa: E731
    al, bet, ga = ang(A_, B_, C_), ang(B_, C_, A_), ang(C_, A_, B_)
    rot_ = lambda c0, th0: (lambda u: c0 + cmath.exp(1j * th0) * (u - c0))  # noqa: E731
    rA, rB, rC = rot_(A_, 2 * al), rot_(B_, 2 * bet), rot_(C_, 2 * ga)
    for P in (0.3 + 0.7j, -2 + 1j, 5j):
        assert abs(rA(rB(rC(P))) - P) < 1e-9
    # contrôle de la correction (relecture) : 2(Â + B̂ + Ĉ) = 2π, partie linéaire identité
    assert abs(2 * (al + bet + ga) - 2 * cmath.pi) < 1e-12

# geo-l3-euclidienne-plan-ex5 : antidéplacement envoyant [AB] sur [A'B'] ; axe = (IJ), vecteur = projeté de AA'


def verif_segments(A_, B_, A2, B2, av, bv, uu, axe_pts, leurres):
    assert abs(B2 - A2) == abs(B_ - A_)
    a_c = cz((B2 - A2) / conjugate(B_ - A_))
    b_c = cz(A2 - a_c * conjugate(A_))
    assert a_c == av and b_c == bv
    assert nature_anti(av, bv) == uu and uu != 0
    I_, J_ = cz((A_ + A2) / 2), cz((B_ + B2) / 2)
    for P in (I_, J_) + tuple(axe_pts):
        assert axe_plan(av, bv, uu, P)
    # projeté de AA' sur la direction u
    d = A2 - A_
    proj = cz(re(d * conjugate(uu)) / (uu * conjugate(uu)) * uu)
    assert proj == uu
    for (ul, Pl) in leurres:
        assert not (cz(ul - uu) == 0 and axe_plan(av, bv, uu, Pl))


verif_segments(0, 2, 1 + I, 1 + 3 * I, I, 1 + I, 1 + I, [0, 5 + 5 * I],
               [(2 + 2 * I, 0), (1 + I, I), (0, 0)])
verif_segments(0, 1, 2, 2 + I, I, 2, 1 + I, [1, 3 + 2 * I], [(2, 1), (1 + I, 0), (2 + 2 * I, 1)])
verif_segments(1, -1, I, 3 * I, -I, 2 * I, -1 + I, [1, 2 - I], [(-2 + 2 * I, 1), (-1 + I, 0), (0, 1)])
assert cz(nature_anti(-I, 2 * I) * 2) == -2 + 2 * I

# ===========================================================================
# Carte espace
# ===========================================================================
# geo-l3-euclidienne-espace (démo) : det(φ - id) = 0 pour φ dans SO(3) ; trace = 1 + 2 cos θ
Rt = rotation_axe([1, 2, 2], th)
assert simplify(Rt.T * Rt - eye(3)) == zeros(3, 3) and simplify(Rt.det()) == 1
assert simplify((Rt - eye(3)).det()) == 0 and simplify(Rt.trace() - 1 - 2 * cos(th)) == 0
assert simplify(Rt * Matrix([1, 2, 2]) - Matrix([1, 2, 2])) == zeros(3, 1)
# trace d'une antirotation diag(-1,1,1) R_θ : -1 + 2 cos θ
Rx = Matrix([[1, 0, 0], [0, cos(th), -sin(th)], [0, sin(th), cos(th)]])
assert simplify((Matrix.diag(-1, 1, 1) * Rx).trace() - (-1 + 2 * cos(th))) == 0
assert simplify(Matrix.diag(-1, 1, 1) * Rx - Rx * Matrix.diag(-1, 1, 1)) == zeros(3, 3)
# énoncé de la carte (relecture) : dans la base (a, e2, e3), φ(a) = -a pour φ = diag(-1,1,1) R_θ
assert Matrix.diag(-1, 1, 1) * Rx * Matrix([1, 0, 0]) == Matrix([-1, 0, 0])
assert simplify(-Rx.subs(th, t) - Matrix.diag(-1, 1, 1) * Rx.subs(th, t + pi)) == zeros(3, 3)
# trou 4 : les leurres 2cosθ, -1+2cosθ, 2+2cosθ diffèrent de 1+2cosθ d'une constante non nulle
for leurre in (2 * cos(th), -1 + 2 * cos(th), 2 + 2 * cos(th)):
    dlt = simplify(1 + 2 * cos(th) - leurre)
    assert dlt.is_number and dlt != 0

# geo-l3-euclidienne-espace (exemple) : vissage (-y+1, x, z+2), axe x = y = 1/2, angle π/2, vecteur (0,0,2)
f = aff([-y + 1, x, z + 2], [x, y, z])
u, g, _ = canonique(f)
assert u == Matrix([0, 0, 2]) and appl(g, [R(1, 2), R(1, 2), 7]) == Matrix([R(1, 2), R(1, 2), 7])
assert angle_signe(f[0], [0, 0, 1]) == (0, 1)
# contre-exemples et affirmations
anti = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, -1]])
assert anti.det() == -1 and (anti - eye(3)).nullspace() == [] and anti != -eye(3)
r90, rm90 = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]), Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
assert r90.trace() == rm90.trace() == 1 and r90 != rm90
assert points_fixes(aff([x + y, y, z + 1], [x, y, z])) is None
assert points_fixes(aff([x + 1, y, -z], [x, y, z])) is None and aff([x + 1, y, -z], [x, y, z])[0].det() == -1
cc = symbols("cc", real=True)
assert egal(compo(reflexion_plan([0, 0, 1], cc), reflexion_plan([0, 0, 1], 0)), trans([0, 0, 2 * cc]))

# geo-l3-euclidienne-espace-ex1 (et v1, v2)
M = Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
assert orthogonale(M) and M.det() == 1 and (M - eye(3)).nullspace()[0] == Matrix([1, 1, 1])
assert angle_signe(M, [1, 1, 1]) == (R(-1, 2), 1)  # 2π/3
assert Matrix.hstack(Matrix([1, 1, 1]), Matrix([1, -1, 0]), M * Matrix([1, -1, 0])).det() == 3
assert simplify(rotation_axe([1, 1, 1], 2 * pi / 3) - M) == zeros(3, 3)
for leurre in (-2 * pi / 3, pi / 3):
    assert simplify(rotation_axe([1, 1, 1], leurre) - M) != zeros(3, 3)
assert M * Matrix([1, -1, 0]) != Matrix([1, -1, 0])
M = M.T  # v1
assert angle_signe(M, [1, 1, 1]) == (R(-1, 2), -1)
assert simplify(rotation_axe([1, 1, 1], -2 * pi / 3) - M) == zeros(3, 3)
assert Matrix.hstack(Matrix([1, 1, 1]), Matrix([1, -1, 0]), M * Matrix([1, -1, 0])).det() == -3
for leurre in (2 * pi / 3, -pi / 3):
    assert simplify(rotation_axe([1, 1, 1], leurre) - M) != zeros(3, 3)
M = Matrix([[0, 0, -1], [0, -1, 0], [-1, 0, 0]])  # v2
assert orthogonale(M) and M.det() == 1 and M.trace() == -1 and M == M.T
assert (M - eye(3)).nullspace()[0].cross(Matrix([1, 0, -1])) == zeros(3, 1)
assert M * Matrix([1, 0, 1]) == Matrix([-1, 0, -1]) and M * Matrix([0, 1, 0]) != Matrix([0, 1, 0])
assert simplify(rotation_axe([1, 0, -1], pi / 2) - M) != zeros(3, 3)

# geo-l3-euclidienne-espace-ex2 (et v1, v2) : antirotations
for (Mm, ax, cth, sg, plan_n) in [
        (Matrix([[0, 0, -1], [-1, 0, 0], [0, -1, 0]]), [1, 1, 1], R(1, 2), -1, [1, 1, 1]),
        (Matrix([[0, -1, 0], [1, 0, 0], [0, 0, -1]]), [0, 0, 1], 0, 1, [0, 0, 1]),
        (Matrix([[0, 0, 1], [-1, 0, 0], [0, 1, 0]]), [1, 1, -1], R(1, 2), 1, [1, 1, -1])]:
    assert orthogonale(Mm) and Mm.det() == -1
    assert Mm * Matrix(ax) == -Matrix(ax) and (Mm - eye(3)).nullspace() == []
    assert angle_signe(Mm, ax) == (cth, sg)
    theta_ok = {(R(1, 2), -1): -pi / 3, (0, 1): pi / 2, (R(1, 2), 1): pi / 3}[(cth, sg)]
    S = reflexion_plan(plan_n, 0)[0]
    assert simplify(rotation_axe(ax, theta_ok) * S - Mm) == zeros(3, 3)
    assert simplify(S * rotation_axe(ax, theta_ok) - Mm) == zeros(3, 3)  # commutation
    for leurre in (-theta_ok, 2 * pi / 3, -2 * pi / 3):
        assert simplify(rotation_axe(ax, leurre) * S - Mm) != zeros(3, 3)
    assert simplify(rotation_axe(ax, theta_ok) - Mm) != zeros(3, 3) and Mm != -eye(3)
# signes cités : bases et v2
assert Matrix.hstack(Matrix([1, 1, 1]), Matrix([1, -1, 0]), Matrix([0, -1, 1])).det() == -3
assert Matrix.hstack(Matrix([1, 1, -1]), Matrix([1, -1, 0]), Matrix([0, -1, -1])).det() == 3

# geo-l3-euclidienne-espace-ex3 (et v1, v2) : vissages


def verif_vissage(fexprs, ub, point, direction, theta_v, leurre_point, gO=None):
    f = aff(fexprs, [x, y, z])
    assert orthogonale(f[0]) and f[0].det() == 1
    u, g, ker = canonique(f)
    assert u == Matrix(ub) and len(ker) == 1 and ker[0].cross(Matrix(direction)) == zeros(3, 1)
    assert simplify(rotation_axe(direction, theta_v) - f[0]) == zeros(3, 3)
    assert appl(g, point) == Matrix(point) and points_fixes(f) is None
    assert f[1] != u  # Of(O) n'est pas le vecteur du vissage
    assert appl(g, leurre_point) != Matrix(leurre_point)  # l'axe leurre n'est pas l'axe
    if gO is not None:
        assert appl(g, [0, 0, 0]) == Matrix(gO)


verif_vissage([z + 1, x, y + 2], (1, 1, 1), (0, -1, 0), (1, 1, 1), 2 * pi / 3, (0, 0, 0), (0, -1, 1))
verif_vissage([y + 2, z - 1, x + 2], (1, 1, 1), (1, 0, 2), (1, 1, 1), -2 * pi / 3, (0, 0, 0), (1, -2, 1))
verif_vissage([-y + 3, x + 1, z - 2], (0, 0, -2), (1, 2, 0), (0, 0, 1), pi / 2, (3, 1, 0))
assert appl(canonique(aff([-y + 3, x + 1, z - 2], [x, y, z]))[1], [3, 1, 0]) == Matrix([2, 4, 0])

# geo-l3-euclidienne-espace-ex4 : R_θ^3 = I pour θ = ±2π/3
for th0 in (2 * pi / 3, -2 * pi / 3):
    assert simplify(Rx.subs(th, th0) ** 3 - eye(3)) == zeros(3, 3)
assert simplify(Rx.subs(th, pi) - Matrix.diag(1, -1, -1)) == zeros(3, 3)

# geo-l3-euclidienne-espace-ex5 : formule de Rodrigues (matrice de r dans la base (a, e2, a∧e2)) et applications
a1, a2_, a3 = symbols("a1 a2 a3", real=True)
av = Matrix([a1, a2_, a3])
Aw = Matrix([[0, -a3, a2_], [a3, 0, -a1], [-a2_, a1, 0]])
Mrod = cos(th) * eye(3) + (1 - cos(th)) * av * av.T + sin(th) * Aw
assert simplify(Mrod.trace() - (1 + 2 * cos(th)) - (1 - cos(th)) * (a1**2 + a2_**2 + a3**2 - 1)) == 0
assert simplify(Mrod - Mrod.T - 2 * sin(th) * Aw) == zeros(3, 3)
e2v = Matrix([1, -1, 0]) / sqrt(2)
a0 = Matrix([1, 1, 1]) / sqrt(3)
e3v = a0.cross(e2v)
P = Matrix.hstack(a0, e2v, e3v)
assert simplify(P.det()) == 1 and simplify(P.T * P - eye(3)) == zeros(3, 3)
assert simplify(P * Rx * P.T - rotation_axe([1, 1, 1], th)) == zeros(3, 3)
assert simplify(a0.cross(e3v) + e2v) == zeros(3, 1)
Mb = Matrix([[2, -1, 2], [2, 2, -1], [-1, 2, 2]]) / 3
assert orthogonale(Mb) and Mb.det() == 1 and Mb.trace() == 2
assert Mb - Mb.T == Matrix([[0, -1, 1], [1, 0, -1], [-1, 1, 0]])
assert simplify(rotation_axe([1, 1, 1], pi / 3) - Mb) == zeros(3, 3)
for leurre, axl in ((-pi / 3, [1, 1, 1]), (2 * pi / 3, [1, 1, 1]), (pi / 3, [1, -1, 1])):
    assert simplify(rotation_axe(axl, leurre) - Mb) != zeros(3, 3)
Mv1 = Mb.T  # v1
assert Mv1 == Matrix([[2, 2, -1], [-1, 2, 2], [2, -1, 2]]) / 3
assert Mv1 - Mv1.T == Matrix([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
assert simplify(rotation_axe([1, 1, 1], -pi / 3) - Mv1) == zeros(3, 3)
for leurre, axl in ((pi / 3, [1, 1, 1]), (-2 * pi / 3, [1, 1, 1]), (-pi / 3, [1, -1, 1])):
    assert simplify(rotation_axe(axl, leurre) - Mv1) != zeros(3, 3)
Mv2 = Matrix([[-1, 2, 2], [2, -1, 2], [2, 2, -1]]) / 3  # v2
assert orthogonale(Mv2) and Mv2.det() == 1 and Mv2.trace() == -1 and Mv2 == Mv2.T
assert simplify(rotation_axe([1, 1, 1], pi) - Mv2) == zeros(3, 3)
assert Mv2 + eye(3) == R(2, 3) * Matrix.ones(3, 3)
for leurre, axl in ((pi / 3, [1, 1, 1]), (pi, [1, -1, 0]), (2 * pi / 3, [1, 1, 1])):
    assert simplify(rotation_axe(axl, leurre) - Mv2) != zeros(3, 3)
assert (Mb.col(0).dot(Mb.col(1)), Mv2.col(0).dot(Mv2.col(1))) == (0, 0)

# ===========================================================================
# Carte generation
# ===========================================================================
# geo-l3-euclidienne-generation (démo) : ρ_c ∘ ρ_b = rotation d'axe a, d'angle θ (a = k, b = i, c = e^{iθ/2})
rho = lambda bb: 2 * bb * bb.T / bb.dot(bb) - eye(3)  # noqa: E731
sref = lambda bb: eye(3) - 2 * bb * bb.T / bb.dot(bb)  # noqa: E731
bvec = Matrix([1, 0, 0])
cvec = Matrix([cos(th / 2), sin(th / 2), 0])
Rz = Matrix([[cos(th), -sin(th), 0], [sin(th), cos(th), 0], [0, 0, 1]])
assert simplify(rho(cvec) * rho(bvec) - Rz) == zeros(3, 3)
assert simplify(sref(cvec) * sref(bvec) - Rz) == zeros(3, 3) and simplify(rho(bvec) + sref(bvec)) == zeros(3, 3)
# affine : ρ_{Δ2} ∘ ρ_{Δ1} = t_u ∘ r avec P1 = O sur l'axe (Oz), P2 = u/2
uz = symbols("uz", real=True)
X = Matrix([x, y, z])
D1 = lambda M_: rho(bvec) * M_  # noqa: E731
P2 = Matrix([0, 0, uz / 2])
D2 = lambda M_: P2 + rho(cvec) * (M_ - P2)  # noqa: E731
assert simplify(D2(D1(X)) - (Rz * X + Matrix([0, 0, uz]))) == zeros(3, 1)
# exemple : ρ_c ∘ ρ_i = (−y, x, z) avec c = (1,1,0)/√2 ; aussi s_c ∘ s_i (plans x = 0 puis x + y = 0)
c1 = Matrix([1, 1, 0])
assert rho(c1) * Matrix([x, y, z]) == Matrix([y, x, -z]) and rho(bvec) * X == Matrix([x, -y, -z])
assert rho(c1) * rho(bvec) == Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) == sref(c1) * sref(bvec)
assert (Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) - eye(3)).rank() == 2
# contre-exemple : transvection non orthogonale ; affirmation : translations par deux demi-tours
assert not orthogonale(Matrix([[1, 1], [0, 1]])) and Matrix([[1, 1], [0, 1]]).det() == 1
assert egal(compo(demi_tour([0, 0, cc], [1, 0, 0]), demi_tour([0, 0, 0], [1, 0, 0])), trans([0, 0, 2 * cc]))
# trou 1 : identité (0 réflexion, pas n ; dim ker = n) ; rotation plane d'angle π/2 : pas de valeur propre -1, rang 2
assert (Matrix([[0, -1], [1, 0]]) - eye(2)).rank() == 2 and (Matrix([[0, -1], [1, 0]]) + eye(2)).det() != 0
# trou 2 : une symétrie glissée du plan (antidéplacement sans point fixe) demande 3 > 2 réflexions ; un vissage 4 > 3
assert aff([x + 1, -y], [x, y])[0].det() == -1 and points_fixes(aff([x + 1, -y], [x, y])) is None

# geo-l3-euclidienne-generation-ex1 (et v1, v2) : rang de φ - id ; décompositions
D = Matrix.diag
assert (D(1, -1, -1) - eye(3)).rank() == 2 and D(1, -1, 1) * D(1, 1, -1) == D(1, -1, -1)
assert (-eye(3) - eye(3)).rank() == 3 and D(-1, 1, 1) * D(1, -1, 1) * D(1, 1, -1) == -eye(3)
Pm = Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
assert (Pm - eye(3)).rank() == 1 and Pm == sref(Matrix([1, -1, 0]))

# geo-l3-euclidienne-generation-ex2 (et v1, v2) : s2 ∘ s1 = translation de 2 H1H2 ; s1 ∘ s2 opposée
for (d1, d2, vec, h1, h2, leurres) in [
        ((1, 1, 1), (1, 1, 3), (2, 2), (R(1, 2), R(1, 2)), (R(3, 2), R(3, 2)), [(1, 1), (-2, -2), (4, 4)]),
        ((1, -2, 0), (1, -2, 5), (2, -4), (0, 0), (1, -2), [(1, -2), (-2, 4), (4, -8)]),
        ((2, -1, -1), (2, -1, 4), (4, -2), (0, 1), (2, 0), [(2, -1), (-4, 2), (8, -4)])]:
    S1, S2 = reflexion_droite(*d1), reflexion_droite(*d2)
    assert egal(compo(S2, S1), trans(vec)) and egal(compo(S1, S2), trans([-v for v in vec]))
    assert sur_droite_plan(h1, *d1) and sur_droite_plan(h2, *d2)
    assert Matrix(vec) == 2 * (Matrix(h2) - Matrix(h1))
    for lv in leurres:
        assert not egal(compo(S2, S1), trans(lv))

# geo-l3-euclidienne-generation-ex3 (et v1, v2) : composées de deux demi-tours
r_1 = demi_tour([0, 0, 0], [1, 0, 0])
assert egal(r_1, aff([x, -y, -z], [x, y, z]))
for (P_, b0, attendu, rho2, theta_v, vec, leurres) in [
        ((0, 0, 1), (1, 1, 0), [-y, x, z + 2], [y, x, 2 - z], pi / 2, (0, 0, 2),
         [(pi / 4, (0, 0, 1)), (pi / 2, (0, 0, 1)), (pi / 2, (0, 0, 0))]),
        ((0, 0, 2), (0, 1, 0), [-x, -y, z + 4], [-x, y, 4 - z], pi, (0, 0, 4),
         [(pi / 2, (0, 0, 2)), (pi, (0, 0, 2)), (0, (0, 0, 4))]),
        ((0, 0, -1), (-1, 1, 0), [y, -x, z - 2], [-y, -x, -2 - z], -pi / 2, (0, 0, -2),
         [(3 * pi / 4, (0, 0, -1)), (pi / 2, (0, 0, -2)), (-pi / 2, (0, 0, -1))])]:
    r_2 = demi_tour(P_, b0)
    assert egal(r_2, aff(rho2, [x, y, z]))
    h = compo(r_2, r_1)
    assert egal(h, aff(attendu, [x, y, z]))
    u, g, ker = canonique(h)
    assert u == Matrix(vec) and simplify(rotation_axe([0, 0, 1], theta_v) - h[0]) == zeros(3, 3)
    assert appl(g, [0, 0, 5]) == Matrix([0, 0, 5])  # g fixe l'axe (Oz)
    for (thl, vl) in leurres:
        assert not egal(h, compo(trans(vl), (rotation_axe([0, 0, 1], thl), zeros(3, 1))))

# geo-l3-euclidienne-generation-ex4 : (-x, -y, z+2) = s'_z ∘ s_z ∘ s_x ∘ s_y
sx, sy, szz, sz1 = (reflexion_plan([1, 0, 0], 0), reflexion_plan([0, 1, 0], 0), reflexion_plan([0, 0, 1], 0),
                    reflexion_plan([0, 0, 1], 1))
assert egal(compo(compo(sz1, szz), compo(sx, sy)), aff([-x, -y, z + 2], [x, y, z]))
assert egal(compo(sz1, szz), trans([0, 0, 2]))

# geo-l3-euclidienne-generation-ex5 : s_a ∘ s_c retournement pour a ⊥ c (dimension 4) ; s_a ∘ s_b en deux retournements
def sref_n(v):
    v = Matrix(v)
    return eye(len(v)) - 2 * v * v.T / v.dot(v)


av4, cv4, bv4 = Matrix([1, 1, 0, 0]), Matrix([1, -1, 2, 0]), Matrix([0, 1, 1, 1])
assert av4.dot(cv4) == 0
Rt4 = sref_n(av4) * sref_n(cv4)
assert Rt4 * av4 == -av4 and Rt4 * cv4 == -cv4 and (Rt4 - eye(4)).rank() == 2 and Rt4.det() == 1
cperp = Matrix.hstack(av4, bv4).T.nullspace()[0]
assert cperp.dot(av4) == 0 and cperp.dot(bv4) == 0
assert sref_n(av4) * sref_n(bv4) == (sref_n(av4) * sref_n(cperp)) * (sref_n(cperp) * sref_n(bv4))
# Q5 : u ρ_D u^{-1} = ρ_{u(D)} ; -s_{a-a'} envoie a sur -a'
uu3 = rotation_axe([1, 2, 2], pi / 5)
dD = Matrix([1, 0, 1])
assert simplify(uu3 * rho(dD) * uu3.T - rho(uu3 * dD)) == zeros(3, 3)
a_1, a_2 = Matrix([1, 0, 0]), Matrix([0, 3, 4]) / 5
U = -sref(a_1 - a_2)
assert simplify(U * a_1 + a_2) == zeros(3, 1) and simplify(U.det()) == 1

# ===========================================================================
# Carte similitudes
# ===========================================================================
# geo-l3-euclidienne-similitudes : centre d'une similitude indirecte et forme réduite (symbolique, k = |a| != 1)
kk = symbols("k", positive=True)
aind = kk * exp(I * th)
om_ind = (aind * conjugate(b) + b) / (1 - kk**2)
assert nul(aind * conjugate(om_ind) + b - om_ind)
f_ = lambda u: aind * conjugate(u) + b  # noqa: E731
h_ = lambda u: om_ind + kk * (u - om_ind)  # noqa: E731
s_ = lambda u: om_ind + exp(I * th) * conjugate(u - om_ind)  # noqa: E731
assert nul(h_(s_(w)) - f_(w)) and nul(s_(h_(w)) - f_(w))
assert nul(s_(om_ind + tt * exp(I * th / 2)) - (om_ind + tt * exp(I * th / 2)))
assert simplify(f_(f_(w)) - (om_ind + kk**2 * (w - om_ind))) == 0  # f∘f homothétie de rapport k² (q6, ex4)
# inverse de a z̄ + b : (z̄ - b̄)/ā
assert nul(f_((conjugate(w) - conjugate(b)) / conjugate(aind)) - w)


def centre_ind(av, bv):
    return cz((av * conjugate(bv) + bv) / (1 - av * conjugate(av)))


def fixe(av, bv, P):
    return cz(av * conjugate(P) + bv - P) == 0


# exemple de la carte : 2i z̄ + 1, centre -(1+2i)/3, axe y = x - 1/3 ; (1+2i)/5 a pour image (9+2i)/5
assert centre_ind(2 * I, 1) == cz(-(1 + 2 * I) / 3) and fixe(2 * I, 1, -(1 + 2 * I) / 3)
assert cz(2 * I * conjugate((1 + 2 * I) / 5) + 1) == cz((9 + 2 * I) / 5)
assert simplify(-R(2, 3) - (-R(1, 3) - R(1, 3))) == 0  # Ω(-1/3, -2/3) sur y = x - 1/3
# affirmation : z ↦ 2z multiplie l'aire du carré unité par 4
assert Matrix([[2, 0], [0, 2]]).det() == 4

# geo-l3-euclidienne-similitudes-ex1 (et v1, v2)
for (av, bv, om0, leurres, naive) in [
        (2 * I, 1, -(1 + 2 * I) / 3, [(1 + 2 * I) / 5, (1 + 2 * I) / 3, -(1 - 2 * I) / 3], (1 + 2 * I) / 5),
        (3, 2 * I, I / 2, [-I, -I / 2, I / 4], -I),
        (-2, 3 * I, -3 * I, [I, 3 * I, -I], I)]:
    assert centre_ind(av, bv) == cz(om0) and fixe(av, bv, om0)
    assert cz(bv / (1 - av) - naive) == 0
    for L in leurres:
        assert not fixe(av, bv, L)
assert cz(3 * conjugate(-I) + 2 * I) == 5 * I and cz(-2 * conjugate(I) + 3 * I) == 5 * I

# geo-l3-euclidienne-similitudes-ex2 (et v1, v2) : similitude indirecte par deux points, centre, axe


def droite_c(P, dvec):
    """Droite passant par P (complexe) dirigée par dvec : test d'appartenance."""
    return lambda Q: cz(im((Q - P) * conjugate(dvec))) == 0


for (A_, B_, A2, B2, av, bv, om0, axe_eq, perp_eq, directe, leurres) in [
        (0, 1, I, 2 + I, 2, I, I / 3, (0, 1, R(1, 3)), (1, 0, 0), (2, I, -I),
         [(I / 3, (1, 0, 0)), (-I, (0, 1, -1))]),
        (0, I, 1, 4, 3 * I, 1, -(1 + 3 * I) / 8, (-1, 1, R(-1, 4)), (1, 1, R(-1, 2)), (-3 * I, 1, (1 - 3 * I) / 10),
         [(-(1 + 3 * I) / 8, (1, 1, R(-1, 2))), ((1 - 3 * I) / 10, (-1, 1, R(-2, 5)))]),
        (1, -1, 0, 4 * I, -2 * I, 2 * I, (4 - 2 * I) / 3, (1, 1, R(2, 3)), (-1, 1, -2), (-2 * I, 2 * I, (4 + 2 * I) / 5),
         [((4 - 2 * I) / 3, (-1, 1, -2)), ((4 + 2 * I) / 3, (1, 1, 2))])]:
    assert cz(av * conjugate(A_) + bv) == cz(A2) and cz(av * conjugate(B_) + bv) == cz(B2)
    assert centre_ind(av, bv) == cz(om0)
    # axe : passe par Ω, dirigé par e^{iθ/2} : on teste que les points de l'axe sont envoyés sur l'axe avec rapport k
    pa, pb, pc = axe_eq
    Om = (re(om0), im(om0))
    assert simplify(pa * Om[0] + pb * Om[1] - pc) == 0
    dvec = -pb + pa * I  # vecteur directeur de la droite pa x + pb y = pc
    img = cz(av * conjugate(om0 + dvec) + bv - om0)
    assert cz(img - abs(av) * dvec) == 0  # f agit sur l'axe comme l'homothétie de rapport k
    qa, qb, qc = perp_eq
    assert simplify(qa * Om[0] + qb * Om[1] - qc) == 0
    dperp = -qb + qa * I
    assert cz(av * conjugate(om0 + dperp) + bv - om0 + abs(av) * dperp) == 0  # rapport -k sur la perpendiculaire
    ad, bd, omd = directe
    assert cz(ad * A_ + bd - A2) == 0 and cz(ad * B_ + bd - B2) == 0 and cz(bd / (1 - ad) - omd) == 0
    for (oml, ax) in leurres:
        la, lb, lc = ax
        bon_centre = fixe(av, bv, oml)
        bon_axe = bon_centre and cz(av * conjugate(oml + (-lb + la * I)) + bv - oml - abs(av) * (-lb + la * I)) == 0
        assert not bon_axe

# geo-l3-euclidienne-similitudes-ex3 (et v1, v2) : f = s ∘ h
for (hc, hk, sfun, av, bv, om0, dvec, leurres) in [
        (1, 3, lambda u: -conjugate(u), -3, 2, R(1, 2), I, [(1, I), (R(1, 2), 1), (-1, I)]),
        (I, 2, lambda u: I * conjugate(u), 2 * I, -1, (1 + 2 * I) / 3, 1 + I, [(I, 1 + I), ((1 + 2 * I) / 3, 1 - I), (-(1 + 2 * I) / 5, 1 + I)]),
        (2, -2, lambda u: conjugate(u), -2, 6, 2, I, [(2, 1), (-2, I), (6, I)])]:
    comp_ = cz(sfun(hc + hk * (w - hc)))
    assert cz(comp_ - (av * conjugate(w) + bv)) == 0
    assert centre_ind(av, bv) == cz(om0) and fixe(av, bv, om0)
    assert cz(av * conjugate(om0 + dvec) + bv - om0 - abs(av) * dvec) == 0  # axe dirigé par dvec
    for (oml, dl) in leurres:
        ok = fixe(av, bv, oml) and cz(av * conjugate(oml + dl) + bv - oml - abs(av) * dl) == 0
        assert not ok
assert centre_ind(-3, 2) == R(1, 2) and fixe(-3, 2, 1) is False

# geo-l3-euclidienne-similitudes-ex4 (et v1, v2) : z_{2n} par itération
nn = symbols("n", integer=True, nonnegative=True)
for (av, bv, om0, formule, leurres) in [
        (-3 * I, 2, (-1 + 3 * I) / 4, lambda n: (1 - 3 * I) * (9**n - 1) / 4,
         [lambda n: (1 - 3 * I) * (3**n - 1) / 4, lambda n: (1 + 3 * I) * (9**n - 1) / 4, lambda n: (1 - 3 * I) * (9**n - 1) / 8]),
        (2, 3 * I, I, lambda n: -I * (4**n - 1),
         [lambda n: -I * (2**n - 1), lambda n: I * (4**n - 1), lambda n: -3 * I * (4**n - 1)]),
        (I / 2, 1, (4 + 2 * I) / 3, lambda n: (4 + 2 * I) * (1 - R(1, 4)**n) / 3,
         [lambda n: (4 + 2 * I) * (1 - R(1, 2)**n) / 3, lambda n: (4 - 2 * I) * (1 - R(1, 4)**n) / 3,
          lambda n: (2 + I) * (1 - R(1, 4)**n) / 2])]:
    assert centre_ind(av, bv) == cz(om0)
    k_ = abs(av)
    zc = 0
    suite_ = [zc]
    for _ in range(14):
        zc = cz(av * conjugate(zc) + bv)
        suite_.append(zc)
    for n in range(7):
        assert cz(suite_[2 * n] - formule(n)) == 0
        assert cz(abs(suite_[2 * n + 1] - om0) - k_ * abs(suite_[2 * n] - om0)) == 0
    for L in leurres:
        assert cz(suite_[2] - L(1)) != 0
assert cz((-3 * I) * conjugate((-1 + 3 * I) / 4) + 2 - (-1 + 3 * I) / 4) == 0

# geo-l3-euclidienne-similitudes-ex5 (et v1, v2) : barycentres et droites invariantes
for (A_, B_, A2, B2, av, bv, k_, om0, g1, g2, leurres) in [
        (0, 2, 2, 2 + 4 * I, 2 * I, 2, 2, -(2 + 4 * I) / 3, (R(2, 3), 2 + R(4, 3) * I), (-2, 2 - 4 * I),
         [R(2, 3), -2, (2 + 4 * I) / 5]),
        (1, 2 + I, 3 + 3 * I, 6, 3, 3 * I, 3, 3 * I / 4, ((6 + 3 * I) / 4, (12 + 3 * I) / 4), (-3 * I / 2, 3 * I / 2),
         [(6 + 3 * I) / 4, -3 * I / 2, 3 * I / 2]),
        (I, -I, 0, 4, -2 * I, 2, 2, (-2 + 4 * I) / 3, (2 * I / 3, (4 - 2 * I) / 3), (2 * I, -4 - 2 * I),
         [2 * I / 3, 2 * I, (2 - 4 * I) / 5])]:
    assert cz(av * conjugate(A_) + bv - A2) == 0 and cz(av * conjugate(B_) + bv - B2) == 0
    assert cz(abs(B2 - A2) / abs(B_ - A_)) == k_ and abs(av) == k_
    assert centre_ind(av, bv) == cz(om0)
    G1 = [cz((k_ * M_ + cz(av * conjugate(M_) + bv)) / (k_ + 1)) for M_ in (A_, B_)]
    G2 = [cz((k_ * M_ - cz(av * conjugate(M_) + bv)) / (k_ - 1)) for M_ in (A_, B_)]
    assert [cz(v) for v in G1] == [cz(v) for v in g1] and [cz(v) for v in G2] == [cz(v) for v in g2]
    d1v, d2v = G1[1] - G1[0], G2[1] - G2[0]
    assert d1v != 0 and d2v != 0 and cz(re(d1v * conjugate(d2v))) == 0  # Δ ⊥ Δ'
    assert droite_c(G1[0], d1v)(om0) and droite_c(G2[0], d2v)(om0)  # Ω = Δ ∩ Δ'
    for L in leurres:
        assert not fixe(av, bv, L)
# formule générale de la question 2 (k réel, p sur Δ, q sur Δ')
pp, qq = symbols("p q")
assert simplify((kk * (pp + qq) + kk * pp - kk * qq) / (kk + 1) - 2 * kk * pp / (kk + 1)) == 0
assert simplify((kk * (pp + qq) - kk * pp + kk * qq) / (kk - 1) - 2 * kk * qq / (kk - 1)) == 0

# ===========================================================================
# Carte figures : force brute
# ===========================================================================


def groupe(points):
    """Isométries conservant l'ensemble fini 'points' (qui engendre l'espace) : liste de (perm, matrice linéaire en
    coordonnées centrées à l'isobarycentre). Chaque permutation qui conserve les distances est réalisée."""
    pts = [Matrix(p) for p in points]
    n = len(pts)
    G = sum(pts, zeros(len(pts[0]), 1)) / n
    cen = [p - G for p in pts]
    d2 = [[(pts[i] - pts[j]).dot(pts[i] - pts[j]) for j in range(n)] for i in range(n)]
    dim = len(pts[0])
    # base : dim vecteurs indépendants parmi les cen
    base = []
    for i in range(n):
        cand = base + [i]
        if Matrix.hstack(*[cen[j] for j in cand]).rank() == len(cand):
            base = cand
        if len(base) == dim:
            break
    assert len(base) == dim
    Bm = Matrix.hstack(*[cen[j] for j in base])
    res = []
    for perm in permutations(range(n)):
        if any(d2[i][j] != d2[perm[i]][perm[j]] for i in range(n) for j in range(i + 1, n)):
            continue
        M_ = Matrix.hstack(*[cen[perm[j]] for j in base]) * Bm.inv()
        assert simplify(M_.T * M_ - eye(dim)) == zeros(dim, dim)  # orthogonale
        assert all(simplify(M_ * cen[i] - cen[perm[i]]) == zeros(dim, 1) for i in range(n))  # réalise perm
        res.append((perm, simplify(M_)))
    return res


def est_reflexion(M_):
    d = M_.shape[0]
    return M_ * M_ == eye(d) and M_.det() == -1 and (M_ - eye(d)).rank() == 1


def signature(p):
    s, vu = 1, set()
    for i in range(len(p)):
        if i in vu:
            continue
        j, lg = i, 0
        while j not in vu:
            vu.add(j)
            j = p[j]
            lg += 1
        s *= (-1) ** (lg - 1)
    return s


# geo-l3-euclidienne-figures (exemple) : rectangle (±2, ±1) : 4 isométries, toutes d'ordre ≤ 2 (non cyclique)
Gr = groupe([(2, 1), (-2, 1), (-2, -1), (2, -1)])
assert len(Gr) == 4 and all(M_ * M_ == eye(2) for _, M_ in Gr)
# carré : 8 (groupe diédral)
assert len(groupe([(1, 1), (-1, 1), (-1, -1), (1, -1)])) == 8
# contre-exemple : triangle équilatéral de l'espace (plan z = 0) : 12 isométries (6 du plan × {id, s_{z=0}})
tri = [Matrix([1, 0]), Matrix([-R(1, 2), sqrt(3) / 2]), Matrix([-R(1, 2), -sqrt(3) / 2])]
G2d = groupe([tuple(p) for p in tri])
assert len(G2d) == 6
iso3 = set()
for _, M_ in G2d:
    for eps in (1, -1):
        M3 = Matrix.diag(M_, eps)
        assert all(any(simplify(M3 * Matrix([p[0], p[1], 0]) - Matrix([q[0], q[1], 0])) == zeros(3, 1) for q in tri) for p in tri)
        iso3.add(tuple(simplify(M3)))
assert len(iso3) == 12 and 6 % 12 != 0
# affirmation : losange non carré (±3,0), (0,±1) : 4 éléments, aucun d'ordre 4
Gl = groupe([(3, 0), (0, 1), (-3, 0), (0, -1)])
assert len(Gl) == 4 and all(M_ * M_ == eye(2) for _, M_ in Gl)

# geo-l3-euclidienne-figures : tétraèdre régulier : 24 permutations réalisées, det = signature (12 déplacements)
TET = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
Gt = groupe(TET)
assert len(Gt) == 24 and all(M_.det() == signature(p) for p, M_ in Gt)
assert sum(1 for p, M_ in Gt if M_.det() == 1) == 12
assert sum(1 for p, M_ in Gt if est_reflexion(M_)) == 6  # figures-ex4 v1 : 6 réflexions = transpositions
assert all(sum(1 for i in range(4) if p[i] != i) == 2 for p, M_ in Gt if est_reflexion(M_))
plans_t = {(0, 1, 1), (1, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, -1), (0, 1, -1)}
for p, M_ in Gt:
    if est_reflexion(M_):
        nv = (M_ - eye(3)).columnspace()[0]
        nv = nv / max(abs(v) for v in nv)
        nv = nv if nv[[i for i in range(3) if nv[i] != 0][0]] > 0 else -nv
        assert tuple(nv) in plans_t

# geo-l3-euclidienne-figures-ex1 (et v1, v2) : 4, 4, 2
assert len(groupe([(3, 1), (-3, 1), (-3, -1), (3, -1)])) == 4
assert len(groupe([(3, 0), (0, 1), (-3, 0), (0, -1)])) == 4
Gi = groupe([(0, 3), (-1, 0), (1, 0)])
assert len(Gi) == 2 and Gi[1][1] == Matrix.diag(-1, 1)

# geo-l3-euclidienne-figures-ex2 (et v1, v2) : prisme triangulaire 12 (6 déplacements), pyramide carrée 8, prisme carré 16
s3 = sqrt(3)
# triangle équilatéral inscrit dans le cercle de rayon 2/√3 (côté 2), hauteur 3
rr = 2 / s3
prisme = [(rr, 0, 0), (-rr / 2, 1, 0), (-rr / 2, -1, 0), (rr, 0, 3), (-rr / 2, 1, 3), (-rr / 2, -1, 3)]
assert simplify((Matrix(prisme[0]) - Matrix(prisme[1])).norm() - 2) == 0
Gp = groupe(prisme)
assert len(Gp) == 12 and sum(1 for _, M_ in Gp if M_.det() == 1) == 6
stab = [p for p, _ in Gp if p[0] == 0]
assert len(stab) == 2 and len({p[0] for p, _ in Gp}) == 6
pyr = [(1, 1, 0), (-1, 1, 0), (-1, -1, 0), (1, -1, 0), (0, 0, 1)]
Gy = groupe(pyr)
assert len(Gy) == 8 and all(p[4] == 4 for p, _ in Gy)
# figures-ex2 v1 (relecture) : 4 réflexions, de plans x = 0, y = 0, x = y, x = -y, et 4 rotations d'axe (Oz)
refl_y = [M_ for _, M_ in Gy if est_reflexion(M_)]
assert len(refl_y) == 4 and all(M_ * Matrix([0, 0, 1]) == Matrix([0, 0, 1]) for _, M_ in Gy)
assert {tuple((M_ - eye(3)).columnspace()[0].normalized()) for M_ in refl_y} <= {
    tuple(Matrix(v).normalized()) for v in ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [1, 1, 0], [-1, -1, 0], [1, -1, 0], [-1, 1, 0])}
# figures-ex2 (relecture) : prisme triangulaire, 6 déplacements = id, 2 rotations d'axe vertical, 3 demi-tours horizontaux
dep_p = [M_ for _, M_ in Gp if M_.det() == 1]
assert sum(1 for M_ in dep_p if M_ == eye(3)) == 1
assert sum(1 for M_ in dep_p if M_ != eye(3) and M_ * Matrix([0, 0, 1]) == Matrix([0, 0, 1])) == 2
assert sum(1 for M_ in dep_p if M_.trace() == -1 and M_ * Matrix([0, 0, 1]) == Matrix([0, 0, -1])) == 3
pc = [(sx_, sy_, sz_) for sx_ in (1, -1) for sy_ in (1, -1) for sz_ in (0, 3)]
Gc = groupe(pc)
assert len(Gc) == 16 and len([p for p, _ in Gc if p[0] == 0]) == 2
assert sum(1 for _, M_ in Gc if est_reflexion(M_)) == 5  # figures-ex4 v2
dist_A = sorted((Matrix(pc[0]) - Matrix(q)).dot(Matrix(pc[0]) - Matrix(q)) for q in pc[1:])
assert dist_A == [4, 4, 8, 9, 13, 13, 17]

# geo-l3-euclidienne-figures-ex3 (et v1, v2) : isométries réalisant (A B C D), (A B)(C D), (A B)
Tm = [Matrix(p) for p in TET]
assert (Tm[0] + Tm[1]) / 2 == Matrix([1, 0, 0]) and (Tm[0] + Tm[2]) / 2 == Matrix([0, 1, 0]) and (Tm[0] + Tm[3]) / 2 == Matrix([0, 0, 1])
for perm, M_att in [((1, 2, 3, 0), Matrix([[0, 0, 1], [0, -1, 0], [-1, 0, 0]])),
                    ((1, 0, 3, 2), Matrix.diag(1, -1, -1)),
                    ((1, 0, 2, 3), Matrix([[1, 0, 0], [0, 0, -1], [0, -1, 0]]))]:
    (M_,) = [M_ for p, M_ in Gt if p == perm]
    assert M_ == M_att
M4 = Matrix([[0, 0, 1], [0, -1, 0], [-1, 0, 0]])
assert M4.det() == -1 and M4 * Matrix([0, 1, 0]) == Matrix([0, -1, 0]) and angle_signe(M4, [0, 1, 0]) == (0, 1)
assert simplify(rotation_axe([0, 1, 0], pi / 2) * Matrix.diag(1, -1, 1) - M4) == zeros(3, 3)
assert M4 * M4 == Matrix.diag(-1, 1, -1)
for Ml in (rotation_axe([0, 1, 0], pi / 2), rotation_axe([1, 1, 1], 2 * pi / 3), -eye(3)):
    assert simplify(Ml - M4) != zeros(3, 3)
assert est_reflexion(Matrix([[1, 0, 0], [0, 0, -1], [0, -1, 0]])) and Matrix([[1, 0, 0], [0, 0, -1], [0, -1, 0]]) * Matrix([0, 1, 1]) == -Matrix([0, 1, 1])
assert Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]]) * Tm[0] == Tm[0]  # réflexion y = z : fixe A
assert rho(Matrix([0, 1, 1])) * Tm[0] == Matrix([-1, 1, 1])  # demi-tour d'axe R(0,1,1)
assert Matrix.diag(-1, 1, 1) * Tm[0] == Matrix([-1, 1, 1]) and -Tm[0] == Matrix([-1, -1, -1])

# geo-l3-euclidienne-figures-ex4 : cube : 48 isométries (matrices de permutation signées), 9 réflexions
CUBE = [(sx_, sy_, sz_) for sx_ in (1, -1) for sy_ in (1, -1) for sz_ in (1, -1)]
Gcu = groupe(CUBE)
assert len(Gcu) == 48
assert all(all(sum(1 for j in range(3) if M_[i, j] != 0) == 1 for i in range(3)) for _, M_ in Gcu)
refl = [M_ for _, M_ in Gcu if est_reflexion(M_)]
assert len(refl) == 9 and all(M_ == M_.T and M_.trace() == 1 for M_ in refl)
assert all(est_reflexion(M_) == (M_ == M_.T and M_.trace() == 1) for _, M_ in Gcu)
assert sum(1 for _, M_ in Gcu if M_.det() == 1) == 24  # q8 : 48 isométries, 24 déplacements

# geo-l3-euclidienne-figures-ex5 : centre de Is(C) = {±id} ; Is(T) ⊂ Is(C) d'indice 2 ; Is(C) = Is(T) ∪ (-Is(T))
mats_c = [M_ for _, M_ in Gcu]
centre_c = [M_ for M_ in mats_c if all(M_ * N_ == N_ * M_ for N_ in mats_c)]
assert sorted(map(tuple, centre_c)) == sorted([tuple(eye(3)), tuple(-eye(3))])
mats_t = [M_ for _, M_ in Gt]
assert all(any(M_ == N_ for N_ in mats_c) for M_ in mats_t)
assert all(any(M_ == N_ for N_ in mats_t) or any(-M_ == N_ for N_ in mats_t) for M_ in mats_c)
assert not any(M_ == -eye(3) for M_ in mats_t)
assert Matrix.hstack(Matrix([0, -2, -2]), Matrix([-2, 0, -2]), Matrix([-2, -2, 0])).det() == -16
assert all(p[0] * p[1] * p[2] == 1 for p in TET)

# ===========================================================================
# QCM
# ===========================================================================
# q1 : (-x+4, y+1) : glissée d'axe x = 2, vecteur (0,1)
f = aff([-x + 4, y + 1], [x, y])
u, g, _ = canonique(f)
assert u == Matrix([0, 1]) and egal(g, reflexion_droite(1, 0, 2)) and points_fixes(f) is None
assert not egal(f, reflexion_droite(1, 0, 2)) and not egal(compo(trans([0, 1]), reflexion_droite(1, 0, 4)), f)
# q2 : z̄ + 1 antidéplacement sans point fixe (vecteur 1)
assert nature_anti(1, 1) == 1
# q3 : 1 + 2 cos θ = 2 : θ = π/3 dans [0, π] ; les leurres ont une autre trace
for th0, ok in ((pi / 3, True), (2 * pi / 3, False), (pi / 6, False), (pi / 4, False)):
    assert (simplify(1 + 2 * cos(th0) - 2) == 0) == ok
# q6 : f∘f homothétie de rapport k² = 9 (vérifié symboliquement plus haut) ; exemple a = 3i
assert simplify((lambda u: 3 * I * conjugate(3 * I * conjugate(u) + b) + b)(w) - (9 * w + 3 * I * conjugate(b) + b)) == 0
# q7 : triangle isocèle non équilatéral : 2 isométries (figures-ex1 v2)
assert len(groupe([(0, 5), (-2, 0), (2, 0)])) == 2

print("geo_l3_euclidienne : toutes les vérifications passent.")
