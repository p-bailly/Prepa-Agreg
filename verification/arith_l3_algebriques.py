"""Vérifications calculatoires du chapitre « arith-l3-algebriques » (arithmétique, L3) : cartes, exercices et QCM de
content/*/arithmetique/l3-arith-l3-algebriques.json et content/qcm/arith-l3-algebriques.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/arith_l3_algebriques.py (doit finir sans erreur).
"""

from fractions import Fraction as F
from itertools import product
from math import comb, factorial

from sympy import (I, Poly, Rational, binomial, cbrt, diff, expand, factor_list, gcdex, minimal_polynomial, rem,
                   roots, sqrt, symbols, nroots, real_roots, simplify, exp, pi, conjugate)

X = symbols("X")
Y = symbols("Y")
n_ = symbols("n", integer=True, positive=True)

PREFIXE = "arith-l3-algebriques"


def minpoly(alpha):
    """Polynôme minimal unitaire sur Q, sous forme d'expression."""
    return Poly(minimal_polynomial(alpha, X), X, domain="QQ").monic().as_expr()


def irreductible_q(P):
    """P irréductible sur Q (un seul facteur, de multiplicité 1, de même degré)."""
    _, facteurs = factor_list(P, X)
    return len(facteurs) == 1 and facteurs[0][1] == 1 and Poly(facteurs[0][0], X).degree() == Poly(P, X).degree()


def fonction(coeffs_ou_P, p):
    """Table des valeurs sur Z/pZ d'un polynôme de Z[X] (expression SymPy)."""
    P = Poly(coeffs_ou_P, X)
    return tuple(int(P.eval(a)) % p for a in range(p))


# ---------------------------------------------------------------------------
# arith-l3-algebriques-algebriques (carte)
# ---------------------------------------------------------------------------
# Exemple : polynôme minimal de 2^(1/3) = X^3 - 2 ; inverse de 1 + alpha = (1 - alpha + alpha^2)/3.
assert expand(minpoly(cbrt(2)) - (X**3 - 2)) == 0
assert all(Poly(X**3 - 2, X).eval(r) != 0 for r in (1, -1, 2, -2))
assert rem(expand((1 + X) * (1 - X + X**2) - 3), X**3 - 2, X) == 0
# Piège 1 / contre-exemple : 1/2 annule 2X - 1 mais aucun polynôme unitaire de Z[X] (racine rationnelle entière) ;
# 3 non inversible dans Z[1/2] : 1/3 n'est pas de la forme a/2^k.
assert all(F(1, 3) * 2**k != int(F(1, 3) * 2**k) for k in range(30))
# Piège 3 : (X - 1)^2 unitaire, annule 1, non irréductible.
assert not irreductible_q((X - 1)**2)
# Contre-exemple : X^4 - 4 = (X^2 - 2)(X^2 + 2), polynôme minimal de sqrt(2) = X^2 - 2 ; i : X^2 + 1 sur Q.
assert expand((X**2 - 2) * (X**2 + 2)) == X**4 - 4 and expand(minpoly(sqrt(2)) - (X**2 - 2)) == 0
assert expand(minpoly(I) - (X**2 + 1)) == 0
# Trou 1 (relecture) : les anciens leurres « scindé sur K », « de degré 1 » laissaient l'énoncé vrai pour alpha dans K
# (pi_alpha = X - alpha), « de degré pair » pour sqrt(2) ; nouveaux leurres « réductible », « inversible », « divisible
# par le carré d'un polynôme non constant », faux pour tout alpha : pi_alpha est irréductible, non constant, et premier
# avec sa dérivée (racines simples). Vérifié sur des exemples, dont le degré 1 (alpha = 3/2).
for alpha in (Rational(3, 2), sqrt(2), cbrt(2), sqrt(2) + sqrt(3), 1 + cbrt(2), I * sqrt(2)):
    m = Poly(minpoly(alpha), X)
    assert irreductible_q(m.as_expr()) and m.degree() >= 1 and m.gcd(m.diff(X)).degree() == 0
# Trou 2 : la constante 1 divise pi_alpha sans annuler alpha, X^2 n'annule pas sqrt(2) ;
# trou 3 : d - 1, d + 1, 2^d sont différents de d pour tout d >= 1.
assert irreductible_q(X**2 - 2) and Poly(minpoly(cbrt(2)), X).degree() == 3  # X^2 - 2 irréductible, donc non scindé sur Q
# Piège 2 (relecture : remplace « P dans C[X] », dont la réfutation demandait un nombre transcendant) et contre-exemple
# « hors d'un anneau intègre » : N = [[0,1],[0,0]] vérifie N^2 = 0, N != 0, polynôme minimal X^2 (réductible),
# N non inversible.
from sympy import Matrix as _Mat
N_ = _Mat([[0, 1], [0, 0]])
assert N_ != _Mat.zeros(2, 2) and N_**2 == _Mat.zeros(2, 2) and N_.det() == 0
assert not irreductible_q(X**2)
assert all(2**d != d and d - 1 != d and d + 1 != d for d in range(1, 50))
# Affirmation 2 : P(X)P(-X) pair, exemple P = X^3 - X - 1.
P = X**3 - X - 1
PP = expand(P * P.subs(X, -X))
assert all(c == 0 for k, c in enumerate(Poly(PP, X).all_coeffs()[::-1]) if k % 2 == 1)
# Affirmation 3 : (X^2 - 2)(X^2 - 3) sans racine rationnelle, annule sqrt(2), non minimal.
P = expand((X**2 - 2) * (X**2 - 3))
assert P.subs(X, sqrt(2)) == 0 and all(P.subs(X, r) != 0 for r in (1, -1, 2, -2, 3, -3, 6, -6))
assert not irreductible_q(P)
# Affirmation 4 : racines simples du polynôme minimal (exemples).
for a in (sqrt(2) + sqrt(3), cbrt(2) + 1, sqrt(2) + cbrt(3)):
    m = Poly(minpoly(a), X)
    assert m.gcd(m.diff(X)).degree() == 0

# arith-l3-algebriques-algebriques-ex1 (+ v1, v2) : polynômes minimaux de nombres quadratiques ; leurres.
for alpha, bon, leurres in [
    (1 + sqrt(3), X**2 - 2 * X - 2, [X**2 - 2 * X + 2, X**2 + 2 * X - 2, X**2 - 2 * X - 4]),
    (2 - sqrt(5), X**2 - 4 * X - 1, [X**2 + 4 * X - 1, X**2 - 4 * X + 1, X**2 - 4 * X - 9]),
    ((1 + sqrt(7)) / 2, X**2 - X - Rational(3, 2), [4 * X**2 - 4 * X - 6, X**2 - X - 3, X**2 + X - Rational(3, 2)]),
]:
    assert expand(minpoly(alpha) - (bon)) == 0
    for L in leurres:
        assert Poly(L, X) != Poly(bon, X)  # 4X^2 - 4X - 6 annule (1 + sqrt 7)/2 mais n'est pas unitaire
    assert [simplify(L.subs(X, alpha)) == 0 for L in leurres] == ([True, False, False] if alpha == (1 + sqrt(7)) / 2
                                                                  else [False, False, False])

# arith-l3-algebriques-algebriques-ex2 (+ v1, v2) : polynômes minimaux de degré 3, racines rationnelles testées.
for alpha, bon, leurres, valeurs in [
    (1 + cbrt(2), X**3 - 3 * X**2 + 3 * X - 3,
     [X**3 - 3 * X**2 + 3 * X - 1, X**3 + 3 * X**2 + 3 * X - 1, X**3 - 3 * X**2 + 3 * X + 1],
     {1: -2, -1: -10, 3: 6, -3: -66}),
    (cbrt(3) - 1, X**3 + 3 * X**2 + 3 * X - 2,
     [X**3 - 3 * X**2 + 3 * X - 4, X**3 + 3 * X**2 + 3 * X + 4, X**3 + 3 * X**2 + 3 * X - 3],
     {1: 5, -1: -3, 2: 24, -2: -4}),
    (2 + cbrt(5), X**3 - 6 * X**2 + 12 * X - 13,
     [X**3 - 6 * X**2 + 12 * X - 3, X**3 + 6 * X**2 + 12 * X + 3, X**3 - 5],
     {1: -6, -1: -32, 13: 1326, -13: -3380}),
]:
    assert expand(minpoly(alpha) - (bon)) == 0 and irreductible_q(bon)
    assert all(bon.subs(X, r) == v for r, v in valeurs.items())
    assert all(simplify(L.subs(X, alpha)) != 0 for L in leurres)

# arith-l3-algebriques-algebriques-ex3 (+ v1, v2) : inverses dans Q[alpha] ; irréductibilité et leurres (coordonnées
# distinctes dans la base (1, alpha, alpha^2)).
for P, A, bon, leurres, valeurs in [
    (X**3 - X - 1, X**2 + 1, Rational(4, 5) + X / 5 - Rational(2, 5) * X**2,
     [Rational(4, 5) - X / 5 + Rational(2, 5) * X**2, 1 - X**2, Rational(4, 5) + X / 5 + Rational(2, 5) * X**2],
     {1: -1, -1: -1}),
    (X**3 + X + 1, X + 2, Rational(5, 9) - Rational(2, 9) * X + X**2 / 9,
     [Rational(5, 9) + Rational(2, 9) * X + X**2 / 9, -Rational(5, 9) + Rational(2, 9) * X - X**2 / 9,
      Rational(1, 5) - Rational(2, 5) * X + X**2], {1: 3, -1: -1}),
    (X**3 - 2 * X - 2, X**2 + 1, Rational(9, 13) + Rational(2, 13) * X - Rational(3, 13) * X**2,
     [Rational(9, 13) - Rational(2, 13) * X + Rational(3, 13) * X**2, Rational(13, 9) + Rational(2, 9) * X - X**2 / 3,
      Rational(9, 13) + Rational(2, 13) * X + Rational(3, 13) * X**2], {1: -3, -1: -1, 2: 2, -2: -6}),
]:
    assert irreductible_q(P) and all(P.subs(X, r) == v for r, v in valeurs.items())
    s, t, g = gcdex(A, P, X)
    assert g == 1 and expand(s - bon) == 0
    assert rem(expand(A * bon - 1), P, X) == 0
    assert all(rem(expand(A * L - 1), P, X) != 0 for L in leurres)
# Divisions écrites dans les corrections.
assert expand(X * (X**2 + 1) - (2 * X + 1) - (X**3 - X - 1)) == 0
assert expand((2 * X + 1) * (X / 2 - Rational(1, 4)) + Rational(5, 4) - (X**2 + 1)) == 0
assert expand((X + 2) * (X**2 - 2 * X + 5) - 9 - (X**3 + X + 1)) == 0
assert expand((3 * X + 2) * (X / 3 - Rational(2, 9)) + Rational(13, 9) - (X**2 + 1)) == 0
assert expand((X**2 + 1) * (1 - X**2 / 2 + X / 4) + (X**3 - X - 1) * (X / 2 - Rational(1, 4)) - Rational(5, 4)) == 0
assert expand((X**2 + 1) * (1 - X**2 / 3 + 2 * X / 9) + (X**3 - 2 * X - 2) * (X / 3 - Rational(2, 9))
              - Rational(13, 9)) == 0

# arith-l3-algebriques-algebriques-ex4 (+ v1, v2) : sqrt(p) + sqrt(q), polynôme de degré 4 irréductible.
for p, q, bon, leurres, rac_cand, b_cand in [
    (2, 3, X**4 - 10 * X**2 + 1, [X**4 + 10 * X**2 + 1, X**4 - 10 * X**2 - 1], (1,), (1, -1)),
    (3, 5, X**4 - 16 * X**2 + 4, [X**4 + 16 * X**2 + 4, X**4 - 16 * X**2 - 4], (1, 2, 4), (2, -2)),
    (2, 5, X**4 - 14 * X**2 + 9, [X**4 + 14 * X**2 + 9, X**4 - 14 * X**2 - 9], (1, 3, 9), (3, -3)),
]:
    alpha = sqrt(p) + sqrt(q)
    assert expand(minpoly(alpha) - (bon)) == 0 and irreductible_q(bon)
    assert all(bon.subs(X, s * r) != 0 for r in rac_cand for s in (1, -1))
    assert all(simplify(L.subs(X, alpha)) != 0 for L in leurres)
    # leurre « X^2 - (p + q) - 2 sqrt(pq) » : annule alpha, coefficients non rationnels.
    assert simplify(alpha**2 - (p + q) - 2 * sqrt(p * q)) == 0
    # Question 3 : cas a = 0 (b + d = -(2(p+q)), bd = (p-q)^2) : discriminant non carré ;
    # cas d = b : a^2 = 2(p+q) + 2b avec b^2 = (p-q)^2, jamais un carré.
    s_ = 2 * (p + q)
    pr = (p - q)**2
    disc = s_**2 - 4 * pr
    assert int(disc**0.5)**2 != disc
    assert all(int((s_ + 2 * b)**0.5)**2 != s_ + 2 * b for b in b_cand)
# Valeurs citées : alpha^4 - 10 alpha^2 - 1 = -2, etc.
assert simplify((sqrt(3) + sqrt(5))**4 - 16 * (sqrt(3) + sqrt(5))**2 - 4) == -8
assert simplify((sqrt(2) + sqrt(5))**4 - 14 * (sqrt(2) + sqrt(5))**2 - 9) == -18
assert [Poly(X**4 - 14 * X**2 + 9, X).eval(r) for r in (1, 3, 9)] == [-4, -36, 5436]
assert [Poly(X**4 - 16 * X**2 + 4, X).eval(r) for r in (1, 2, 4)] == [-11, -44, 4]

# arith-l3-algebriques-algebriques-ex5 : nombre de Liouville. Minoration k! >= (n+1)! + (k-n-1) pour k >= n+1,
# majoration du reste, et sommes partielles exactes pour n petit.
for n in range(1, 8):
    assert all(factorial(k) >= factorial(n + 1) + (k - n - 1) for k in range(n + 1, n + 12))
for n in range(1, 4):
    qn = 10**factorial(n)
    pn = sum(10**(factorial(n) - factorial(k)) for k in range(1, n + 1))
    reste = sum(F(1, 10**factorial(k)) for k in range(n + 1, n + 4))  # termes suivants, négligeables au-delà
    assert 0 < reste <= F(2, qn**(n + 1)) and F(10, 9) <= 2
    assert qn**(n + 1) == 10**factorial(n + 1)

# ---------------------------------------------------------------------------
# arith-l3-algebriques-denombrabilite (carte)
# ---------------------------------------------------------------------------
def polys_hauteur(H):
    res = []
    for n in range(H):
        S = H - n
        for c in product(range(-S, S + 1), repeat=n + 1):
            if c[-1] != 0 and sum(abs(x) for x in c) == S:
                res.append(c)
    return res


nb = {H: len(polys_hauteur(H)) for H in range(1, 6)}
# Exemple : hauteurs 1 et 2 ; racines des hauteurs <= 3 puis <= 4.
assert nb[1] == 2 and nb[2] == 4
rac = {}
for Hmax in (3, 4):
    R_ = set()
    for H in range(1, Hmax + 1):
        for c in polys_hauteur(H):
            P = sum(ci * X**i for i, ci in enumerate(c))
            R_ |= set(roots(P, X).keys())
    rac[Hmax] = R_
assert rac[3] == {0, 1, -1}
assert rac[4] == {0, 1, -1, Rational(1, 2), Rational(-1, 2), 2, -2, I, -I}
# Piège 1 / contre-exemple : tout complexe z est racine de X^2 - 2 Re(z) X + |z|^2.
for z in (2 + 3 * I, -1 + I / 2, Rational(5)):
    assert simplify((X**2 - 2 * z.as_real_imag()[0] * X + abs(z)**2).subs(X, z)) == 0
# Démonstration : bijection ]a,b[ -> R par x -> tan(pi (x-a)/(b-a) - pi/2) (bornes -pi/2, pi/2).
# Trous (relecture) : l'ancien leurre « fini » du trou 2 laissait l'énoncé vrai pour A = {0} (aucun polynôme non nul,
# ensemble vide de racines) ; « infinie », « bornée », « fermée » du trou 1 le laissaient vrai pour A = Z, {0, 1}, {0}.
# Nouveaux leurres : A non dénombrable, ouverte ou d'intérieur non vide (donc non dénombrable) ; alors, pour a1 dans A
# non nul, les racines -a0/a1 des polynômes a1 X + a0 (a0 dans A) forment une partie en bijection avec A : l'énoncé est
# faux pour toute telle A. Trou 2 : un ouvert non vide de C n'est pas dénombrable, l'ensemble des racines l'est toujours.
assert all(roots(a1 * X + a0, X) == {-Rational(a0) / a1: 1} for a1 in (1, 2, -3) for a0 in (0, 5, Rational(7, 2)))
assert [ci for H in range(1, 5) for ci in polys_hauteur(H) if all(x in (0,) for x in ci)] == []  # A = {0} : aucun
# Affirmation 4 : si P(t+q) = 0 alors P(X+q) annule t (identité de substitution).
assert expand((X**2 - 2).subs(X, X + Rational(1, 3))) == X**2 + 2 * X / 3 - Rational(17, 9)

# arith-l3-algebriques-denombrabilite-ex1 (+ v1, v2) : hauteur exactement 4 (24), exactement 3 (10), au plus 4 (40).
assert nb[4] == 24 and nb[3] == 10 and nb[1] + nb[2] + nb[3] + nb[4] == 40
assert all(x != 24 for x in (8, 20, 40)) and all(x != 10 for x in (6, 16, 8)) and all(x != 40 for x in (24, 38, 16))
# Détail de la correction (hauteur 4) : 2 + 10 + 10 + 2.
par_degre = {}
for c in polys_hauteur(4):
    par_degre[len(c) - 1] = par_degre.get(len(c) - 1, 0) + 1
assert par_degre == {0: 2, 1: 10, 2: 10, 3: 2}
# arith-l3-algebriques-denombrabilite-ex2 : P(3X+1) de coefficient dominant 3^n * a_n ; X^n P(1/X) de terme constant a_n.
P = 2 * X**3 - X + 5
assert Poly(P.subs(X, 3 * X + 1), X).LC() == 27 * 2
assert Poly(expand(X**3 * P.subs(X, 1 / X)), X).eval(0) == 2
# arith-l3-algebriques-denombrabilite-ex3 et ex4 : démonstrations (dénombrabilité) ; ex4, vecteur propre de C^T
# pour la matrice compagnon (vérifié pour un polynôme de degré 4).
from sympy import Matrix, zeros
c = [3, -1, 0, 2]  # P = X^4 + 2X^3 - X + 3
d = 4
C = zeros(d, d)
for i in range(d - 1):
    C[i + 1, i] = 1
for i in range(d):
    C[i, d - 1] = -c[i]
P = X**4 + sum(c[i] * X**i for i in range(d))
assert expand(C.charpoly(X).as_expr() - P) == 0
for lam in nroots(P, n=30):
    v = Matrix([lam**j for j in range(d)])
    assert max(abs(complex(z.evalf(30))) for z in (C.T * v - lam * v)) < 1e-20

# ---------------------------------------------------------------------------
# arith-l3-algebriques-fonctions-zp (carte)
# ---------------------------------------------------------------------------
for p in (2, 3, 5, 7, 11, 13):
    # X^p - X = prod (X - a) dans F_p[X].
    prod_ = Poly(1, X, modulus=p)
    for a in range(p):
        prod_ = prod_ * Poly(X - a, X, modulus=p)
    assert prod_ == Poly(X**p - X, X, modulus=p)
    # X^{p-1} - 1 = prod_{a != 0} (X - a) (affirmation 5) ; X^p - X + 1 sans racine (affirmation 3).
    prod2 = Poly(1, X, modulus=p)
    for a in range(1, p):
        prod2 = prod2 * Poly(X - a, X, modulus=p)
    assert prod2 == Poly(X**(p - 1) - 1, X, modulus=p)
    assert all((a**p - a + 1) % p == 1 % p for a in range(p))
    # E_a = 1 - (X - a)^{p-1} vaut 1 en a, 0 ailleurs (démonstration).
    for a in range(p):
        assert fonction(1 - (X - a)**(p - 1), p) == tuple(1 if x == a else 0 for x in range(p))
# Nombre de fonctions polynômes = p^p (force brute sur les polynômes de degré < p + 2 pour p = 2, 3).
for p in (2, 3):
    fs = {tuple(sum(c[k] * a**k for k in range(len(c))) % p for a in range(p)) for c in product(range(p), repeat=p + 2)}
    assert len(fs) == p**p
# Noyau : X^p et X même fonction ; deux polynômes de même fonction diffèrent d'un multiple de X^p - X (p = 3, degré <= 4).
p = 3
polys = list(product(range(p), repeat=5))
for c1 in polys[::7]:
    for c2 in polys[::11]:
        P1 = sum(c1[k] * X**k for k in range(5))
        P2 = sum(c2[k] * X**k for k in range(5))
        meme = fonction(P1, p) == fonction(P2, p)
        assert meme == (Poly(P1 - P2, X, modulus=p).rem(Poly(X**3 - X, X, modulus=p)).is_zero)
# Pièges et contre-exemples sur Z/4Z : X^4 - X vaut 2 en 2 ; 64 fonctions polynômes ; l'indicatrice de 0 n'en est pas une.
assert (2**4 - 2) % 4 == 2
fs4 = {tuple(sum(c[k] * a**k for k in range(6)) % 4 for a in range(4)) for c in product(range(4), repeat=6)}
assert len(fs4) == 64 and (1, 0, 0, 0) not in fs4
# Piège 2 : X^{p-1} - 1 ne s'annule pas en 0 ; piège 3 : X^p et X, de degré <= p, même fonction.
assert all(fonction(X**p_ - X, p_) == fonction(0 * X, p_) for p_ in (2, 3, 5, 7))
# Trou 1 : n = 9 (impair) : X^9 - X ne s'annule pas en 3 modulo 9 ; et toutes les fonctions ne sont pas polynomiales.
assert (3**9 - 3) % 9 != 0
fs9 = {tuple(sum(c[k] * a**k for k in range(4)) % 9 for a in range(9)) for c in product(range(9), repeat=4)}
assert (1, 0, 0, 0, 0, 0, 0, 0, 0) not in fs9  # f(3) - f(0) doit être multiple de 3 pour une fonction polynôme
assert all((P_3 - P_0) % 3 == 0 for (P_0, _, _, P_3, *_) in fs9)
# Trou 1 (relecture) : les anciens leurres « entier >= 2 », « impair >= 3 », « puissance de 2 » laissaient l'énoncé
# vrai pour n = 3, 3, 2 (premiers) ; nouveaux leurres, tous composés : pour n composé et d diviseur propre de n,
# d divise P(d) - P(0) pour tout P de Z[X], alors que l'indicatrice de 0 vérifie f(d) - f(0) = -1 : elle n'est pas
# polynomiale, et l'énoncé est faux pour tout n composé (4, 8, 16, 6, 10, 15... testés ci-dessous).
import random
random.seed(1)
for n in (4, 6, 8, 9, 10, 14, 15, 16, 21, 25, 32, 35):
    d = next(x for x in range(2, n) if n % x == 0)
    for _ in range(200):
        c_ = [random.randint(-50, 50) for _ in range(random.randint(1, 8))]
        Pd, P0 = sum(ci * d**i for i, ci in enumerate(c_)), c_[0]
        assert (Pd - P0) % d == 0
    assert (0 - 1) % d != 0
# Trou 3 : p!, p^{p-1} différents de p^p pour p premier.
assert all(factorial(p_) != p_**p_ and p_**(p_ - 1) != p_**p_ for p_ in (2, 3, 5, 7))
# Exemple : X^5 + 2X^4 + X et 2X^2 + 2X, même fonction sur F_3 ; indicatrice de 0 = 1 - X^2.
assert fonction(X**5 + 2 * X**4 + X, 3) == fonction(2 * X**2 + 2 * X, 3) == (0, 1, 0)
assert fonction(1 - X**2, 3) == (1, 0, 0)

# arith-l3-algebriques-fonctions-zp-ex1 (+ v1, v2) : réduction modulo X^p - X ; leurres de fonctions différentes.
for p, P, bon, leurres, table in [
    (5, X**13 + 2 * X**6, 2 * X**2 + X, [X**3 + 2 * X, X**2 + 2 * X, X + 2], (0, 3, 0, 1, 1)),
    (3, X**10 + X**5 + 1, X**2 + X + 1, [X**2 + 2 * X + 1, X + 1, 2 * X**2 + 1], (1, 0, 1)),
    (7, X**20 + 3 * X**8, 4 * X**2, [X**6 + 3 * X, 3 * X**2 + 1, 4 * X], (0, 4, 2, 1, 1, 2, 4)),
]:
    assert Poly(P, X, modulus=p).rem(Poly(X**p - X, X, modulus=p)) == Poly(bon, X, modulus=p)
    assert fonction(P, p) == fonction(bon, p) == table
    assert all(fonction(L, p) != table for L in leurres)
assert fonction(X**3 + 2 * X, 5)[2] == 2 and fonction(4 * X, 7)[2] == 1

# arith-l3-algebriques-fonctions-zp-ex2 (+ v1, v2) : interpolation ; leurres.
for p, table, bon, leurres in [
    (3, (2, 1, 1), 2 * X**2 + 2, [X**2 + 2, 2 * X + 2, 2 * X**2 + X + 2]),
    (3, (2, 0, 0), X**2 + 2, [2 * X**2 + 2, 2 * X + 2, X**2 + X + 2]),
    (5, (0, 0, 1, 0, 0), 4 * X**4 + 3 * X**3 + X**2 + 2 * X,
     [4 * X**4 + 2 * X**3 + X**2 + 3 * X, 4 * X**4 + 1, X**4 + 3 * X**3 + X**2 + 2 * X]),
]:
    interp = sum(table[a] * (1 - (X - a)**(p - 1)) for a in range(p))
    assert Poly(interp, X, modulus=p) == Poly(bon, X, modulus=p)
    assert fonction(bon, p) == table and all(fonction(L, p) != table for L in leurres)
assert expand(1 - (X - 2)**4) == -X**4 + 8 * X**3 - 24 * X**2 + 32 * X - 15

# arith-l3-algebriques-fonctions-zp-ex3 (+ v1, v2) : sommes de puissances (théorème vérifié par force brute).
for p in (2, 3, 5, 7, 11, 13):
    for k in range(1, 40):
        S = sum(pow(x, k, p) for x in range(p)) % p
        assert S == ((p - 1) % p if k % (p - 1) == 0 else 0)
for p, termes, bon, leurres in [
    (11, lambda x: x**3 + 2 * x**10 + x**20, 8, (3, 0, 9)),
    (7, lambda x: x**5 + 3 * x**6 + x**12 + 2, 3, (4, 5, 0)),
    (13, lambda x: x**12 + x**24 + 5 * x**7, 11, (2, 12, 0)),
]:
    S = sum(termes(x) for x in range(p)) % p
    assert S == bon and bon not in leurres

# arith-l3-algebriques-fonctions-zp-ex4 (+ v1, v2) : coefficient de X^{p-1} = -somme des valeurs.
for p, table, bon, leurres in [
    (7, (2, 5, 1, 1, 0, 3, 6), 3, (4, 0, 1)),
    (5, (4, 4, 1, 0, 3), 3, (2, 0, 4)),
    (11, (1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0), 8, (3, 0, 10)),
]:
    interp = Poly(sum(table[a] * (1 - (X - a)**(p - 1)) for a in range(p)), X, modulus=p)
    coef = interp.coeff_monomial(X**(p - 1)) % p
    assert coef == bon == (-sum(table)) % p and bon not in leurres
# Question 2 : bijections de Z/5Z et Z/7Z -> degré <= p - 2 (échantillon) ; p = 2, x -> x + 1 de degré 1.
from itertools import permutations
for p in (3, 5):
    for perm in permutations(range(p)):
        assert Poly(sum(perm[a] * (1 - (X - a)**(p - 1)) for a in range(p)), X, modulus=p).degree() <= p - 2
assert Poly(sum(((a + 1) % 2) * (1 - (X - a)) for a in range(2)), X, modulus=2) == Poly(X + 1, X, modulus=2)

# arith-l3-algebriques-fonctions-zp-ex5 : Z/4Z, P(x+2) = P(x) + 2P'(x) mod 4, 64 fonctions, construction de P.
for c in product(range(4), repeat=4):
    P = sum(c[k] * X**k for k in range(4))
    for x in range(-5, 6):
        assert (P.subs(X, x + 2) - P.subs(X, x) - 2 * diff(P, X).subs(X, x)) % 4 == 0
cond = {f for f in product(range(4), repeat=4) if (f[2] - f[0]) % 2 == 0 and (f[3] - f[1]) % 2 == 0}
cubiques = {tuple(sum(c[k] * a**k for k in range(4)) % 4 for a in range(4)) for c in product(range(4), repeat=4)}
assert cond == fs4 == cubiques and len(cond) == 64 and all(x != 64 for x in (256, 16, 128))

# ---------------------------------------------------------------------------
# arith-l3-algebriques-alembert-gauss (carte)
# ---------------------------------------------------------------------------
# Exemple : P = X^2 + 1, omega = i, P(ti) = 1 - t^2 ; coercivité |z|^2 - 1 > 1 pour |z| > sqrt(2).
t = symbols("t", positive=True)
assert expand((X**2 + 1).subs(X, t * I)) == 1 - t**2
# Démonstration : formule P(z0 + t w) = P(z0)(1 - t^k) + t^{k+1} w^{k+1} Q(t w) et majoration, sur un exemple.
P = X**3 - 2 * X + 4 + I
z0 = 1 + I
Pz = expand(P.subs(X, z0 + Y))
coeffs = Poly(Pz, Y).all_coeffs()[::-1]
k = next(j for j in range(1, len(coeffs)) if coeffs[j] != 0)
w = (-coeffs[0] / coeffs[k])**Rational(1, k)
Q = expand((Pz - coeffs[0] - coeffs[k] * Y**k) / Y**(k + 1))
M = sum(abs(complex(cc)) for cc in Poly(expand(w**(k + 1) * Q.subs(Y, w * Y)), Y).all_coeffs())
for tt in (0.5 * min(1, abs(complex(coeffs[0])) / (M + 1)), 1e-3, 1e-2):
    assert abs(complex(P.subs(X, z0 + tt * w))) < abs(complex(P.subs(X, z0)))
# Contre-exemples : |e^z| = e^{Re z} ; 1 + |z|^2 >= 1. Affirmation 2 : P = X, P'(0) = 1.
# Trous (relecture) : anciens leurres du trou 3, « P'(z1) = 0 », « z1 = 0 », « P(z1) = P(0) », vrais pour P = X^2
# (z1 = 0) ; trou 1, « non nul », « unitaire », « à coefficients réels », vrais pour P = X. Nouveaux leurres : trou 1,
# P constant (|P| ne tend pas vers +infini) ; trou 3, « P(z1) != 0 », « |P(z1)| > |P(0)| », « |P(z1)| = 1 », jamais vrais
# puisque z1 est un point de minimum (|P(z1)| <= |P(0)|) et une racine.
for P_ in (X, X**2, X**3 - 6 * X + 3, X**2 + 2 * X + 5, (X - 3)**2 * (X + I)):
    z1 = min(nroots(P_, n=30), key=lambda r_: abs(complex(P_.subs(X, r_))))
    assert abs(complex(P_.subs(X, z1))) < 1e-20 <= abs(complex(P_.subs(X, 0))) + 1e-20
assert expand(X**2).subs(X, 0) == 0 and diff(X**2, X).subs(X, 0) == 0  # P'(z1) = 0 pour P = X^2 : ancien leurre vrai

# arith-l3-algebriques-alembert-gauss-ex1 (+ v1, v2) : directions de descente ; leurres (|P| augmente).
for P, z0, bon, leurres in [
    (X**4 + 4, 0, 1 + I, [sqrt(2), I, sqrt(2) * I]),
    (X**3 + 8, 0, 1 + I * sqrt(3), [2, 2 * I, -1 + I * sqrt(3)]),
    (X**2 + 2 * X + 5, -1, 2 * I, [2, 1 + I, -2]),
]:
    v0 = abs(complex(P.subs(X, z0)))
    for tt in [j / 20 for j in range(1, 21)]:
        assert abs(complex(P.subs(X, z0 + tt * bon))) < v0
        assert all(abs(complex(P.subs(X, z0 + tt * L))) > v0 for L in leurres)
assert expand((1 + I)**4) == -4 and expand((1 + I * sqrt(3))**3) == -8 and expand((-1 + I * sqrt(3))**3) == 8
assert expand((X**2 + 2 * X + 5).subs(X, -1 + Y)) == Y**2 + 4

# arith-l3-algebriques-alembert-gauss-ex2 (+ v1, v2) : plus petit entier R tel que m(r) > |P(0)| pour r >= R.
for m, P0, bon, leurres in [
    (lambda r: r**3 - 6 * r - 3, 3, 3, (2, 4, 6)),
    (lambda r: r**4 - 3 * r - 1, 1, 2, (1, 3, 4)),
    (lambda r: r**3 - 10 * r - 4, 4, 4, (3, 5, 10)),
]:
    ok = lambda R: all(m(R + j / 100) > P0 for j in range(0, 2000))
    assert ok(bon) and not any(m(r) > P0 for r in [bon - 1]) and not ok(bon - 1)
    assert all(L != bon for L in leurres)

# arith-l3-algebriques-alembert-gauss-ex3 (+ v1, v2) : minimum de |P| sur le disque unité (grille fine).
import cmath
for P, bon, leurres, rac in [
    (lambda z: z**2 + 3, 2, (3, 4, 1), [cmath.sqrt(-3)]),
    (lambda z: z**3 - 4, 3, (4, 5, 1), [4**(1 / 3)]),
    (lambda z: 2 * z + 5, 3, (5, 7, 2), [-2.5]),
]:
    m = min(abs(P(r / 60 * cmath.exp(2j * cmath.pi * s / 720))) for r in range(61) for s in range(720))
    assert abs(m - bon) < 1e-6 and all(abs(L - bon) > 0.5 for L in leurres) and all(abs(z) > 1 for z in rac)

# arith-l3-algebriques-alembert-gauss-ex4 (+ v1, v2) : sommes de deux carrés ; leurres.
for P, Qc, bon, leurres in [
    (X**4 + 2 * X**3 + 3 * X**2 + 2 * X + 2, (X - I) * (X + 1 - I), (X**2 + X - 1)**2 + (2 * X + 1)**2,
     [(X**2 + X - 1)**2 + (2 * X - 1)**2, (X**2 + X + 1)**2 + (2 * X + 1)**2, (X**2 - X + 1)**2 + 1]),
    (X**4 + 5 * X**2 + 4, (X - I) * (X - 2 * I), (X**2 - 2)**2 + (3 * X)**2,
     [(X**2 + 2)**2 + (3 * X)**2, (X**2 - 2)**2 + X**2, (X**2 + 1)**2 + 3]),
    (X**4 - 2 * X**3 + 2 * X**2 - 2 * X + 1, (X - 1) * (X - I), (X**2 - X)**2 + (X - 1)**2,
     [(X**2 + X)**2 + (X - 1)**2, (X**2 - X)**2 + (X + 1)**2, (X**2 - 1)**2 + (X - 1)**2]),
]:
    Qe = expand(Qc)
    A, B = [expand(z) for z in Qe.as_real_imag()] if False else (None, None)
    Qp = Poly(Qe, X)
    A = sum(complex(cc).real * X**i for i, cc in enumerate(Qp.all_coeffs()[::-1]))
    B = sum(complex(cc).imag * X**i for i, cc in enumerate(Qp.all_coeffs()[::-1]))
    assert expand(Qe * expand(conjugate(Qe)).subs(conjugate(X), X) - P) == 0
    assert expand(A**2 + B**2 - P).simplify() == 0 or max(abs(cc) for cc in Poly(expand(A**2 + B**2 - P), X).all_coeffs()) < 1e-12
    assert expand(bon - P) == 0 and all(expand(L - P) != 0 for L in leurres)
# Autres écritures citées : (X^2 + X + 1)^2 + 1 et (X^2 + 2)^2 + X^2.
assert expand((X**2 + X + 1)**2 + 1) == X**4 + 2 * X**3 + 3 * X**2 + 2 * X + 2
assert expand((X**2 + 2)**2 + X**2) == X**4 + 5 * X**2 + 4

# ---------------------------------------------------------------------------
# arith-l3-algebriques-taylor (carte)
# ---------------------------------------------------------------------------
def taylor(P, a):
    d = Poly(P, X).degree()
    return [diff(P, X, k).subs(X, a) / factorial(k) for k in range(d + 1)]


assert taylor(X**3 - 2 * X + 1, 1) == [0, 1, 3, 1]
assert expand((X - 1) + 3 * (X - 1)**2 + (X - 1)**3) == X**3 - 2 * X + 1
# Identité générale (degré <= 6, coefficients symboliques) et piège « somme tronquée », « factorielles oubliées ».
cs = symbols("c0:7")
a = symbols("a")
P = sum(cs[k] * X**k for k in range(7))
assert expand(sum(diff(P, X, k).subs(X, a) / factorial(k) * (X - a)**k for k in range(7)) - P) == 0
assert sum(diff(X**4, X, k).subs(X, 0) / factorial(k) * X**k for k in range(4)) == 0
assert expand(sum(diff(X**2, X, k).subs(X, 0) * X**k for k in range(3))) == 2 * X**2
# Contre-exemple caractéristique 2 : dérivée formelle de X^2 = 2X = 0 dans F_2[X].
assert Poly(diff(X**2, X), X, modulus=2).is_zero
# Trous : leurres P^{(k)}(0)/k! et P^{(k)}(a)/n! faux (P = X, a = 1 ; P = X^2, n = 2).
assert taylor(X, 1)[0] == 1 and Poly(X, X).eval(0) == 0
assert taylor(X**2, 3)[0] == 9 and F(9, 2) != 9
# Trou 2 (relecture) : n >= 2 dans le texte à trous (pour n <= 1, P^{(k)}(a) et P^{(k)}(a)/n! coïncidaient avec la
# réponse) ; leurre P^{(k)}(0)/k! (vrai pour a = 0) remplacé par P^{(k+1)}(a)/(k+1)!, faux pour P = 1 (c_0 = 1, et 0).
assert taylor(Rational(1) + 0 * X, 5) == [1] and diff(1 + 0 * X, X) == 0
assert all(taylor(X**2, a_)[2] == 1 and diff(X**2, X, 2) == 2 and F(1, factorial(n)) != 1
           for a_ in (0, 1, -2) for n in range(2, 8))
# Affirmation 1 : P(X + a) à coefficients entiers (exemple).
assert all(cc.is_integer for cc in taylor(7 * X**5 - 3 * X**2 + X - 4, -3))

# arith-l3-algebriques-taylor-ex1 (+ v1, v2) : coefficients de Taylor ; leurres différents de P.
for P, a, bon, leurres in [
    (X**3 + 2 * X**2 - X + 3, 1, [5, 6, 5, 1], [[5, 6, 10, 6], "base", [5, 6, 5, 1, "plus"]]),
    (X**3 - 3 * X**2 + 4 * X - 2, 2, [2, 4, 3, 1], [[2, 4, 6, 6], "base", [2, 4, 3, 1, "plus"]]),
    (X**4 + X**3 - X + 2, 1, [3, 6, 9, 5, 1], [[3, 6, 18, 30, 24], "base", [3, 6, 9, 5, 1, "plus"]]),
]:
    assert taylor(P, a) == bon
    faux1 = sum(c_ * (X - a)**k for k, c_ in enumerate(leurres[0]))
    faux2 = sum(c_ * (X - a)**k for k, c_ in enumerate(Poly(P, X).all_coeffs()[::-1]))
    faux3 = sum(c_ * (X + a)**k for k, c_ in enumerate(bon))
    assert all(expand(f - P) != 0 for f in (faux1, faux2, faux3))
# Leurre de v1 « -2 + 4(X-2) - 3(X-2)^2 + (X-2)^3 » : coefficients de P (signes) dans la base en X - 2.
assert expand(-2 + 4 * (X - 2) - 3 * (X - 2)**2 + (X - 2)**3 - (X**3 - 3 * X**2 + 4 * X - 2)) != 0
assert expand(2 - (X - 1) + (X - 1)**3 + (X - 1)**4 - (X**4 + X**3 - X + 2)) != 0

# arith-l3-algebriques-taylor-ex2 (+ v1, v2) : restes de X^n par (X-1)^3, (X+1)^3, (X-2)^2, n de 2 à 12 ; leurres.
for n in range(2, 13):
    r1 = 1 + n * (X - 1) + F(n * (n - 1), 2) * (X - 1)**2
    assert expand(rem(X**n, (X - 1)**3, X) - r1) == 0
    for L in (1 + n * (X - 1) + n * (n - 1) * (X - 1)**2, 1 + (X - 1) + (X - 1)**2,
              n * (X - 1) + F(n * (n - 1), 2) * (X - 1)**2):
        assert expand(L - r1) != 0
    r2 = (-1)**n * (1 - n * (X + 1) + F(n * (n - 1), 2) * (X + 1)**2)
    assert expand(rem(X**n, (X + 1)**3, X) - r2) == 0
    for L in (-r2, (-1)**n * (1 + n * (X + 1) + F(n * (n - 1), 2) * (X + 1)**2),
              (-1)**n * (1 - n * (X + 1) + n * (n - 1) * (X + 1)**2)):
        assert expand(L - r2) != 0
    r3 = 2**n + n * 2**(n - 1) * (X - 2)
    assert expand(rem(X**n, (X - 2)**2, X) - r3) == 0
    for L in (2**n + n * 2**n * (X - 2), 2**n + n * (X - 2), 1 + n * 2**(n - 1) * (X - 2)):
        assert expand(L - r3) != 0

# arith-l3-algebriques-taylor-ex3 (+ v1, v2) : relèvement de Hensel ; leurres.
for p, P, a0, bon, leurres in [
    (5, lambda x: x**2 - 6, 1, 16, (6, 11, 21)),
    (7, lambda x: x**2 - 2, 3, 10, (3, 17, 24)),
    (5, lambda x: x**3 - 3, 2, 12, (7, 17, 22)),
]:
    sol = [x for x in range(p * p) if P(x) % (p * p) == 0 and x % p == a0]
    assert sol == [bon] and all(P(L) % (p * p) != 0 for L in leurres)
# Congruence P(a + m) = P(a) + m P'(a) mod m^2 (force brute).
for c in product(range(-3, 4), repeat=4):
    Pp = lambda x: sum(c[k] * x**k for k in range(4))
    dP = lambda x: sum(k * c[k] * x**(k - 1) for k in range(1, 4))
    for a_ in range(-3, 4):
        for m in range(1, 6):
            assert (Pp(a_ + m) - Pp(a_) - m * dP(a_)) % (m * m) == 0

# arith-l3-algebriques-taylor-ex4 (+ v1, v2) : puissances p-ièmes dans F_p[X] ; leurres.
for p, P, bon, leurres in [
    (5, X**10 + 2 * X**5 + 1, (X + 1)**10, [(X**2 + 1)**5, (X + 2)**10, (X + 1)**5]),
    (3, X**6 + X**3 + 1, (X + 2)**6, [(X + 1)**6, (X**2 + 1)**3, (X + 2)**3]),
    (7, X**14 + 3 * X**7 + 2, (X + 1)**7 * (X + 2)**7, [(X + 1)**14, (X + 3)**7 * (X + 2)**7, (X + 5)**7 * (X + 6)**7]),
]:
    assert Poly(bon, X, modulus=p) == Poly(P, X, modulus=p)
    assert all(Poly(L, X, modulus=p) != Poly(P, X, modulus=p) for L in leurres)
    assert Poly(diff(P, X), X, modulus=p).is_zero
# Question 2 : Q^p = Q(X^p) dans F_p[X] (exemples).
for p in (2, 3, 5, 7):
    Qx = 3 * X**3 + X + 2
    assert Poly(Qx**p, X, modulus=p) == Poly(Qx.subs(X, X**p), X, modulus=p)

# arith-l3-algebriques-taylor-ex5 (+ v1, v2) : plus petit entier a, toutes les racines réelles < a.
for P, bon, leurres in [
    (X**3 - 6 * X**2 + 5 * X + 3, 5, (4, 6, 3)),
    (X**3 - 7 * X**2 + 10 * X - 1, 6, (5, 7, 4)),
    (X**4 - 4 * X**3 + 2 * X**2 - X - 1, 4, (3, 5, 2)),
]:
    d = Poly(P, X).degree()
    ok = [a_ for a_ in range(0, 20) if P.subs(X, a_) > 0 and all(diff(P, X, k).subs(X, a_) >= 0 for k in range(1, d + 1))]
    assert ok[0] == bon and bon not in leurres
    assert max(real_roots(P)) < bon and P.subs(X, bon - 1) < 0
# Valeurs citées dans les corrections.
assert [(X**3 - 6 * X**2 + 5 * X + 3).subs(X, v) for v in (4, 5)] == [-9, 3]
assert [(3 * X**2 - 12 * X + 5).subs(X, v) for v in (2, 3, 5)] == [-7, -4, 20]
assert [(X**3 - 7 * X**2 + 10 * X - 1).subs(X, v) for v in (4, 5, 6)] == [-9, -1, 23]
assert [(4 * X**3 - 12 * X**2 + 4 * X - 1).subs(X, v) for v in (2, 3, 4)] == [-9, 11, 79]
assert [(12 * X**2 - 24 * X + 4).subs(X, v) for v in (1, 2, 4)] == [-8, 4, 100]

# ---------------------------------------------------------------------------
# QCM
# ---------------------------------------------------------------------------
# arith-l3-algebriques-q1 : polynôme minimal de i sqrt(2) = X^2 + 2 ; X^2 - 2 vaut -4 en i sqrt(2).
assert expand(minpoly(I * sqrt(2)) - (X**2 + 2)) == 0 and expand((I * sqrt(2))**2 - 2) == -4
# arith-l3-algebriques-q5 : 27 fonctions polynômes sur Z/3Z (force brute).
assert len({tuple(sum(c[k] * a**k for k in range(5)) % 3 for a in range(3)) for c in product(range(3), repeat=5)}) == 27
# arith-l3-algebriques-q6 : X^8 et X^2 même fonction sur Z/7Z.
assert fonction(X**8, 7) == fonction(X**2, 7)
# arith-l3-algebriques-q8 : (X - 2)^3 non nul, P(2) = P'(2) = P''(2) = 0.
assert taylor((X - 2)**3, 2)[:3] == [0, 0, 0]

print("arith_l3_algebriques : toutes les vérifications passent.")
