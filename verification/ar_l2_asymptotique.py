"""Vérifications calculatoires du chapitre « ar-l2-asymptotique » (analyse réelle, L2) :
cartes, exercices et QCM de content/*/analyse-reelle/l2-ar-l2-asymptotique.json et content/qcm/ar-l2-asymptotique.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variantes v1, v2) ou de la question.
Lancer : python3 verification/ar_l2_asymptotique.py (doit finir sans erreur, en quelques secondes).
"""

import math
from fractions import Fraction as F

import mpmath as mp
from sympy import (E, Rational, cos, exp, expand, factor, limit, log, oo, pi, series, simplify, sin, sqrt, symbols,
                   tan, Poly, solve, diff, atan)

mp.mp.dps = 40
n = symbols("n", positive=True)
x, y, a, b, c, d = symbols("x y a b c d")


def it(f, u0, N):
    u = [u0]
    for _ in range(N):
        u.append(f(u[-1]))
    return u


# ---------------------------------------------------------------------------
# ar-l2-asymptotique-comparaison
# ---------------------------------------------------------------------------
# Exemples de la carte.
assert limit((3 * n**2 + n * log(n)) / (2 * n**3 - 1) / (3 / (2 * n)), n, oo) == 1
assert limit((sqrt(n + 1) - sqrt(n)) * 2 * sqrt(n), n, oo) == 1
assert limit(n * log(1 + 1 / n), n, oo) == 1
# Contre-exemples / pièges : 1/n et 2/n (différence -> 0, rapport 1/2) ; n et 2n ; somme n+1 ~ n+2, -n.
assert limit((n + 1) / (n + 2), n, oo) == 1 and F(1, 2) != 1
assert limit(exp(n + 1) / exp(n), n, oo) == E
# Trou 3 : les leurres (u^2 ~ v^2, |u| ~ |v|, u v' ~ v u') sont vrais en général (produits d'équivalents) ;
# test numérique sur u = n+1, v = n+2, u' = -n + 1/n, v' = -n.
uu, vv, up, vp = n + 1, n + 2, -n + 1 / n, -n
assert limit(uu**2 / vv**2, n, oo) == 1 and limit(uu * vp / (vv * up), n, oo) == 1
assert limit((uu + vp) / (vv + vp), n, oo) == Rational(1, 2)  # somme : 1 / 2, pas équivalentes
# Affirmation 3 : ln u ~ ln v si u ~ v, v -> +oo (exemple u = 2n + 5, v = 2n).
assert limit(log(2 * n + 5) / log(2 * n), n, oo) == 1

# ar-l2-asymptotique-comparaison-ex1 (+ v1, v2) et leurres.
for expr, eq, leurres in [
    ((3 * n**2 - n * log(n) + 1) / (2 * n**3 + sqrt(n)), 3 / (2 * n), [3 / (2 * n**2), 2 / (3 * n), Rational(3, 2)]),
    ((5 * n**3 + n**2 * log(n)) / (4 * n**5 - n), 5 / (4 * n**2), [5 / (4 * n), 4 / (5 * n**2), 5 / (4 * n**3)]),
    ((2 * n**2 + n * sqrt(n)) / (7 * n**4 - 3 * n), 2 / (7 * n**2), [2 / (7 * n), 7 / (2 * n**2), 3 / (7 * n**2)]),
]:
    assert limit(expr / eq, n, oo) == 1
    assert all(limit(expr / l, n, oo) != 1 for l in leurres)

# ar-l2-asymptotique-comparaison-ex2 (+ v1, v2).
for expr, eq, leurres in [
    (sqrt(n + 2) - sqrt(n), 1 / sqrt(n), [2 / sqrt(n), 1 / (2 * sqrt(n)), 1 / n]),
    ((n + 3) ** Rational(1, 3) - n ** Rational(1, 3), n ** Rational(-2, 3),
     [3 * n ** Rational(-2, 3), n ** Rational(-2, 3) / 3, n ** Rational(-1, 3)]),
    (sqrt(4 * n + 1) - 2 * sqrt(n), 1 / (4 * sqrt(n)), [1 / (2 * sqrt(n)), 1 / sqrt(n), 1 / (4 * n)]),
]:
    for N in (10**8, 10**12):  # rapports évalués numériquement à 40 chiffres
        r = mp.mpf(str((expr / eq).subs(n, N).evalf(40)))
        assert abs(r - 1) < 1e-3
        assert all(abs(mp.mpf(str((expr / l).subs(n, N).evalf(40))) - 1) > 0.2 for l in leurres)

# ar-l2-asymptotique-comparaison-ex3 (+ v1, v2) : n^2 ln(1 + a/n) - a n -> -a^2/2.
for k, lim, leurres in [(1, -Rational(1, 2), [0, Rational(1, 2), -oo]), (-2, -2, [0, 2, -4]),
                        (3, -Rational(9, 2), [0, -Rational(3, 2), Rational(9, 2)])]:
    s = series(n**2 * log(1 + k / n) - k * n, n, oo, 2).removeO()
    assert limit(s, n, oo) == lim and lim == -Rational(k * k, 2)
    assert lim not in leurres  # leurres : exposants 0 (limite 1), etc.
# Valeurs numériques directes (u_n -> e^{-1/2} etc.).
for k, lim in [(1, -0.5), (-2, -2.0), (3, -4.5)]:
    N = mp.mpf(10) ** 6
    w = N**2 * mp.log(1 + k / N) - k * N
    assert abs(w - lim) < 1e-4

# ar-l2-asymptotique-comparaison-ex4 : contre-exemple l = 1 et e^{n+1/n} ~ e^n.
assert limit(log(1 + 2 / n) / log(1 + 1 / n), n, oo) == 2
assert limit(exp(n + 1 / n) / exp(n), n, oo) == 1

# QCM q1 : ln(1+1/n) - 1/n ~ -1/(2n^2) ; leurres faux.
s = (log(1 + 1 / n) - 1 / n)
assert limit(s / (-1 / (2 * n**2)), n, oo) == 1
assert all(limit(s / l, n, oo) != 1 for l in [1 / (2 * n**2), -1 / n**2, -1 / (3 * n**3)])

# ---------------------------------------------------------------------------
# ar-l2-asymptotique-recurrente
# ---------------------------------------------------------------------------
# Exemple : u_{n+1} = sqrt(2 + u_n), u_0 = 0 -> 2.
u = it(lambda t: mp.sqrt(2 + t), mp.mpf(0), 80)
assert all(u[i] < u[i + 1] for i in range(40)) and abs(u[-1] - 2) < 1e-30
# Contre-exemple / pièges : f(x) = 1 - x, u_0 = 0 : 0, 1, 0, 1 ; (u_{3n}) et (u_{3n+1}) non monotones.
u = it(lambda t: 1 - t, 0, 12)
assert u[:4] == [0, 1, 0, 1] and [u[0], u[3], u[6]] == [0, 1, 0] and [u[1], u[4], u[7]] == [1, 0, 1]
# Trou 2 : f(x) = x/2 sur [-1, 1], u_0 = 1 : u_0, u_1, f(0) >= 0 mais suite décroissante.
u = it(lambda t: t / 2, F(1), 5)
assert u[0] >= 0 and u[1] >= 0 and all(u[i + 1] < u[i] for i in range(5))
# Contre-exemple de continuité : f(0) = 1, f(x) = x/2 sinon ; u_n = 2^-n -> 0, f(0) = 1.
assert it(lambda t: 1 if t == 0 else t / 2, F(1), 10)[-1] == F(1, 1024)

# ar-l2-asymptotique-recurrente-ex1 (+ v1, v2) : sqrt(m + u) -> racine positive de l^2 - l - m.
for m, l, leurres in [(6, 3, [-2, mp.sqrt(6)]), (12, 4, [-3, mp.sqrt(12)]), (20, 5, [-4, mp.sqrt(20)])]:
    assert l * l - l - m == 0 and leurres[0] ** 2 - leurres[0] - m == 0
    u = it(lambda t: mp.sqrt(m + t), mp.mpf(0), 100)
    assert abs(u[-1] - l) < 1e-30 and all(u[i] <= u[i + 1] <= l for i in range(100))
    assert mp.sqrt(m) <= l  # f(0) = sqrt(m) dans [0, l]

# ar-l2-asymptotique-recurrente-ex2 : u - u^2 -> 0, décroissante.
u = it(lambda t: t - t * t, F(1, 2), 10)
assert all(0 < u[i + 1] < u[i] < 1 for i in range(10))

# ar-l2-asymptotique-recurrente-ex3 (+ v1, v2) : u^2 + c, points fixes a < b ; discussion.
for cc, pa, pb in [(F(3, 16), F(1, 4), F(3, 4)), (F(2, 9), F(1, 3), F(2, 3)), (F(6, 25), F(2, 5), F(3, 5))]:
    assert expand((x - pa) * (x - pb)) == expand(x**2 - x + cc) and 0 < cc < pa
    f = lambda t: t * t + cc
    for u0 in [F(0), pa / 2, pa, (pa + pb) / 2, pb - F(1, 1000)]:
        v = mp.mpf(u0.numerator) / u0.denominator
        for _ in range(3000):
            v = v * v + mp.mpf(cc.numerator) / cc.denominator
        assert abs(v - mp.mpf(pa.numerator) / pa.denominator) < 1e-20
    assert f(pb) == pb
    t = mp.mpf(pb.numerator) / pb.denominator + mp.mpf(1) / 100
    for _ in range(60):
        t = t * t + mp.mpf(cc.numerator) / cc.denominator
        if t > 10**6:
            break
    assert t > 10**6  # u_0 > b : divergence vers +oo
    # Leurres : seuil a (ex : u_0 dans ]a, b] -> b faux), seuil 1/2 faux, limite b sur [0, b] faux (cf. ci-dessus).

# ar-l2-asymptotique-recurrente-ex4 (+ v1, v2) : f(x) = 1 - x^2/k.
for k, pf, cyc, lim, leurres in [
    (2, x**2 + 2 * x - 2, x**2 - 2 * x + 2, sqrt(3) - 1, [(sqrt(5) - 1) / 2, Rational(1, 2)]),
    (3, x**2 + 3 * x - 3, x**2 - 3 * x + 6, (sqrt(21) - 3) / 2, [(3 - sqrt(5)) / 2, Rational(2, 3)]),
    (4, x**2 + 4 * x - 4, x**2 - 4 * x + 12, 2 * sqrt(2) - 2, [sqrt(2) - 1, Rational(3, 4)]),
]:
    f = 1 - x**2 / k
    ff = f.subs(x, f)
    assert expand(ff - x + pf * cyc / k**3) == 0  # f(f(x)) - x = -(1/k^3) pf cyc
    assert Poly(cyc, x).discriminant() < 0
    assert simplify(pf.subs(x, lim)) == 0 and 0 < lim.evalf() < 1
    assert all(abs(float(l) - float(lim)) > 0.05 for l in leurres)
    u = it(lambda t: 1 - t * t / k, mp.mpf(0), 200)
    assert abs(u[-1] - mp.mpf(str(lim.evalf(40)))) < 1e-20
    assert [F(1 - F(0) ** 2 / k)] == [1]
# Premiers termes cités : u_3 = 7/8, 23/27, 55/64.
for k, u3 in [(2, F(7, 8)), (3, F(23, 27)), (4, F(55, 64))]:
    assert it(lambda t: 1 - t * t / k, F(0), 3)[3] == u3

# QCM q3 : f(x) = 1/(1+x), u_0 = 0 ; f o f = (1+x)/(2+x).
assert simplify((1 / (1 + 1 / (1 + x))) - (1 + x) / (2 + x)) == 0
assert it(lambda t: 1 / (1 + t), F(0), 3) == [0, 1, F(1, 2), F(2, 3)]
# QCM q4 : (u^2 + 1)/2 - u = (u - 1)^2 / 2.
assert expand((x**2 + 1) / 2 - x - (x - 1) ** 2 / 2) == 0

# ---------------------------------------------------------------------------
# ar-l2-asymptotique-homographique
# ---------------------------------------------------------------------------
# Identité clé h(x) - g = (a - c g)(x - g)/(cx + d) quand g est point fixe (b = c g^2 + (d - a) g).
g = symbols("g")
bb = c * g**2 + (d - a) * g
assert simplify((a * x + bb) / (c * x + d) - g - (a - c * g) * (x - g) / (c * x + d)) == 0
# Le polynôme des points fixes vaut (ad - bc)/c en -d/c.
assert simplify((c * x**2 + (d - a) * x - b).subs(x, -d / c) - (a * d - b * c) / c) == 0
# Point fixe double : 1/(h(x) - alpha) = 2c/(a+d) + 1/(x - alpha), alpha = (a-d)/(2c), b = -(a-d)^2/(4c).
al = (a - d) / (2 * c)
bd = -(a - d) ** 2 / (4 * c)
assert simplify(1 / ((a * x + bd) / (c * x + d) - al) - 2 * c / (a + d) - 1 / (x - al)) == 0
# k = (a - c alpha)/(a - c beta) = (c beta + d)/(c alpha + d) (alpha + beta = (a - d)/c).
be = symbols("be")
alp = (a - d) / c - be
assert simplify((a - c * alp) - (c * be + d)) == 0


def hom(A, B, Cc, D):
    return lambda t: (A * t + B) / (Cc * t + D)


def est_geom(v):
    return all(v[i] != 0 for i in range(len(v) - 1)) and len({v[i + 1] / v[i] for i in range(len(v) - 1)}) == 1


def est_arith(v):
    return len({v[i + 1] - v[i] for i in range(len(v) - 1)}) == 1


# Exemple : (4x - 2)/(x + 1), u_0 = 0 ; v_n = 2 (2/3)^n ; u_1 = -2, u_2 = 10, u_3 = 38/11.
h = hom(4, -2, 1, 1)
u = it(h, F(0), 6)
assert u[1:4] == [-2, 10, F(38, 11)]
assert all((u[i] - 2) / (u[i] - 1) == 2 * F(2, 3) ** i for i in range(7))
# Pièges et leurres (avec u_0 = 3 pour la même h ; points fixes 2 et 1) :
u = it(h, F(3), 6)
assert not est_geom([(t + 1) / (t + 2) for t in u])  # racines de x^2 + (a - d) x - b
assert not est_geom([(t - 2) / t for t in u])  # beta = 0 quelconque
assert not est_geom([(t - 4) / (t + 1) for t in u])  # a/c et -d/c
assert not est_geom([t - 2 for t in u]) and not est_geom([(t - 2) * (t - 1) for t in u])
assert not est_geom([1 / (t - 2) for t in u])
vv = [(t - 2) / (t - 1) for t in u]
assert est_geom(vv) and vv[1] / vv[0] == F(2, 3) and not est_arith(vv)  # raison 2/3, ni 3/2, ni 2
# Trou 3 : h(x) = 1/(2 - x), u_0 = 0 : u_n = n/(n+1) ; 1/(u-1) arithmétique, les leurres non.
u = it(hom(0, 1, -1, 2), F(0), 8)
assert u == [F(k, k + 1) for k in range(9)]  # QCM q6
assert est_arith([1 / (t - 1) for t in u]) and (1 / (u[1] - 1)) - (1 / (u[0] - 1)) == -1
assert not est_arith([t - 1 for t in u]) and not est_arith([(t - 1) / (t + 1) for t in u])
lnv = [mp.log(abs(t - 1)) for t in u]
assert abs((lnv[2] - lnv[1]) - (lnv[1] - lnv[0])) > 0.1
# Contre-exemples : 1/(1-x) depuis 0 tombe sur le pôle ; -1/x et 1/x périodiques.
assert hom(0, 1, -1, 1)(F(0)) == 1
assert it(hom(0, -1, 1, 0), F(1), 4) == [1, -1, 1, -1, 1]
assert it(hom(0, 1, 1, 0), F(2), 4) == [2, F(1, 2), 2, F(1, 2), 2]
assert F(0 - 1, 0 + 1) == -1  # k = (a - c alpha)/(a - c beta) pour 1/x avec (1, -1)

# ar-l2-asymptotique-homographique-ex1 (+ v1, v2).
for (A, B, Cc, D), al_, be_, k, formule in [
    ((5, -3, 1, 1), 3, 1, F(1, 2), lambda m: F(3 * 2**m - 3, 2**m - 3)),
    ((3, 2, 1, 2), 2, -1, F(1, 4), lambda m: F(2 * 4**m - 2, 4**m + 2)),
    ((2, 3, 1, 4), 1, -3, F(1, 5), lambda m: F(3 * 5**m - 3, 3 * 5**m + 1)),
]:
    h = hom(A, B, Cc, D)
    assert h(F(al_)) == al_ and h(F(be_)) == be_ and A * D - B * Cc != 0
    assert F(A - Cc * al_, A - Cc * be_) == k
    u = it(h, F(0), 12)
    assert all(t != -F(D, Cc) for t in u)
    assert u == [formule(m) for m in range(13)]
    v = [(t - al_) / (t - be_) for t in u]
    assert all(v[i + 1] == k * v[i] for i in range(12))
assert it(hom(5, -3, 1, 1), F(0), 2)[1:] == [-3, 9]

# ar-l2-asymptotique-homographique-ex2 (+ v1, v2) et leurres.
for (A, B, Cc, D), u0, al_, formule, leurres in [
    ((3, -1, 1, 1), 2, 1, lambda m: F(m + 4, m + 2),
     [lambda m: F(m + 6, m + 3), lambda m: 1 + F(1, 2**m), lambda m: F(m + 2, m + 1)]),
    ((5, -4, 1, 1), 3, 2, lambda m: F(2 * m + 9, m + 3),
     [lambda m: F(2 * m + 12, m + 4), lambda m: 2 + F(1, 3**m), lambda m: F(3 * m + 9, m + 3)]),
    ((1, -4, 1, -3), 1, 2, lambda m: F(2 * m + 1, m + 1),
     [lambda m: 2 - F(1, 2**m), lambda m: F(2 * m + 2, m + 2), lambda m: F(3 * m + 2, 2 * m + 2)]),
]:
    h = hom(A, B, Cc, D)
    assert solve(Cc * x**2 + (D - A) * x - B, x) == [al_]
    u = it(h, F(u0), 15)
    assert u == [formule(m) for m in range(16)]
    w = [1 / (t - al_) for t in u]
    assert est_arith(w) and w[1] - w[0] == F(2 * Cc, A + D)
    for l in leurres:
        assert any(l(m) != u[m] for m in range(16))

# ar-l2-asymptotique-homographique-ex3 (+ v1, v2) : 1/(m + u) -> alpha = racine positive, k = -alpha^2.
for m, u0, termes, lim in [(1, F(1), [F(1, 2), F(2, 3), F(3, 5), F(5, 8)], (sqrt(5) - 1) / 2),
                           (3, F(0), [F(1, 3), F(3, 10), F(10, 33)], (sqrt(13) - 3) / 2),
                           (2, F(0), [F(1, 2), F(2, 5), F(5, 12)], sqrt(2) - 1)]:
    u = it(lambda t: 1 / (m + t), u0, 60)
    assert u[1:len(termes) + 1] == termes
    rac = solve(x**2 + m * x - 1, x)
    assert any(simplify(r - lim) == 0 for r in rac)
    beta_ = [r for r in rac if simplify(r - lim) != 0][0]
    k = simplify(lim / beta_)
    assert simplify(k + lim**2) == 0 and abs(float(k)) < 1
    assert abs(float(u[-1]) - float(lim)) < 1e-12
assert simplify(-((sqrt(5) - 1) / 2) ** 2 + (3 - sqrt(5)) / 2) == 0
# Fibonacci : u_n = F_{n+1}/F_{n+2}.
fib = [0, 1]
while len(fib) < 30:
    fib.append(fib[-1] + fib[-2])
assert it(lambda t: 1 / (1 + t), F(1), 20) == [F(fib[k + 1], fib[k + 2]) for k in range(21)]
# Pell : 1 + u_n = 3/2, 7/5, 17/12.
assert [1 + t for t in it(lambda t: 1 / (2 + t), F(0), 3)[1:]] == [F(3, 2), F(7, 5), F(17, 12)]

# ar-l2-asymptotique-homographique-ex4 : (u + m)/(u + 1).
mm = symbols("m", positive=True)
gm = symbols("gm")
assert simplify((x + gm**2) / (x + 1) - gm - (1 - gm) * (x - gm) / (x + 1)) == 0
u = it(lambda t: (t + 2) / (t + 1), F(1), 3)
assert u == [1, F(3, 2), F(7, 5), F(17, 12)]
k = (1 - sqrt(2)) / (1 + sqrt(2))
assert simplify(k + (3 - 2 * sqrt(2))) == 0
v3 = (3 - 2 * sqrt(2)) ** 4
assert simplify(v3 * (Rational(17, 12) + sqrt(2)) - (Rational(17, 12) - sqrt(2))) == 0  # égalité exacte
assert 8.6e-4 < float(v3) < 8.7e-4 and float(Rational(17, 12) - sqrt(2)) < 2.5e-3
assert abs(float((2 + sqrt(2)) * abs(k) ** 3) - 1.7e-2) < 1e-3
for m in [F(1, 4), F(2), F(9)]:  # convergence vers sqrt(m), majoration (max(1, m) + sqrt m) |k|^n
    sm = math.sqrt(m)
    kk = abs((1 - sm) / (1 + sm))
    for u0 in [F(0), F(5)]:
        u = it(lambda t: (t + m) / (t + 1), u0, 30)
        assert all(abs(float(u[j]) - sm) <= (max(1, float(m)) + sm) * kk**j + 1e-15 for j in range(1, 31))

# QCM q5 : (2u+1)/(u+2), points fixes ±1, k = 1/3, limite 1.
u = it(hom(2, 1, 1, 2), F(0), 40)
assert abs(float(u[-1]) - 1) < 1e-15 and F(2 - 1, 2 + 1) == F(1, 3)

# ---------------------------------------------------------------------------
# ar-l2-asymptotique-vitesse
# ---------------------------------------------------------------------------
# Exemple Héron : u_1..u_3, f'(sqrt 2) = 0, identité exacte, erreurs.
f = (x + 2 / x) / 2
assert simplify(diff(f, x).subs(x, sqrt(2))) == 0
assert simplify(f - sqrt(2) - (x - sqrt(2)) ** 2 / (2 * x)) == 0
u = it(lambda t: (t + 2 / t) / 2, F(1), 3)
assert u == [1, F(3, 2), F(17, 12), F(577, 408)]
err = [float(t) - math.sqrt(2) for t in u[1:]]
assert 0.08 < err[0] < 0.09 and 2.4e-3 < err[1] < 2.5e-3 and 2.0e-6 < err[2] < 2.2e-6
# Cosinus : l ≈ 0,739, -sin l ≈ -0,674.
l = mp.findroot(lambda t: mp.cos(t) - t, 0.7)
assert abs(l - 0.739) < 1e-3 and abs(-mp.sin(l) + 0.674) < 1e-3
# Contre-exemples : sin, u_n ~ sqrt(3/n) (vérifié numériquement) ; x + x^3 ; -2x ; x^2.
u = mp.mpf(1)
for _ in range(10**5):
    u = mp.sin(u)
assert abs(u / mp.sqrt(3 / mp.mpf(10**5)) - 1) < 1e-3
assert it(lambda t: t + t**3, F(1, 10), 5)[-1] > F(1, 10)
assert it(lambda t: -2 * t, 1, 4)[-1] == 16 and it(lambda t: t * t, 2, 4)[-1] == 2**16
# Trou 2 : f(x) = -x/2 : rapport -1/2 (et non |f'| = 1/2) ; f = x/2 : rapport 1/2 (ni 0, ni 2).
assert F(-1, 2) != abs(F(-1, 2))

# ar-l2-asymptotique-vitesse-ex1 (+ v1, v2) : Héron et leurres, erreurs < 1e-4.
for A, u0, u1, u2, leurres, e1 in [(3, F(2), F(7, 4), F(97, 56), [F(26, 15), F(19, 11)], 0.018),
                                   (5, F(2), F(9, 4), F(161, 72), [F(38, 17), F(47, 21)], 0.014),
                                   (7, F(3), F(8, 3), F(127, 48), [F(45, 17), F(37, 14)], 0.021)]:
    u = it(lambda t: (t + A / t) / 2, u0, 2)
    assert u[1:] == [u1, u2] and u1 not in leurres and u2 not in leurres
    r = math.sqrt(A)
    assert 0 < float(u1) - r < e1 and 0 <= float(u2) - r < 1e-4
    assert (e1**2) / (2 * float(u1)) < 1e-4
assert abs((97 / 56 - math.sqrt(3)) - 9.2e-5) < 1e-6 and abs((161 / 72 - math.sqrt(5)) - 4.3e-5) < 1e-6
assert abs((127 / 48 - math.sqrt(7)) - 8.2e-5) < 5e-7  # ar-l2-asymptotique-vitesse-ex1 v2 (8,2·10^-5)
assert math.sqrt(3) > 1.732 and math.sqrt(5) > 2.236 and math.sqrt(7) > 2.6457

# ar-l2-asymptotique-vitesse-ex2 (+ v1, v2) : (u^2 + a)/(a+1), rapport (u+1)/(a+1) -> 2/(a+1).
for A, lim, leurres in [(2, F(2, 3), [F(4, 3), 0, F(1, 3)]), (3, F(1, 2), [F(3, 2), 0, F(1, 4)]),
                        (4, F(2, 5), [F(8, 5), 0, F(1, 5)])]:
    C = A + 1
    assert expand((x**2 + A) / C - x - (x - 1) * (x - A) / C) == 0
    assert diff((x**2 + A) / C, x).subs(x, 1) == lim and diff((x**2 + A) / C, x).subs(x, A) == leurres[0]
    u = it(lambda t: (t * t + A) / C, F(0), 12)
    assert all(t < 1 for t in u)
    assert all((u[i + 1] - 1) / (u[i] - 1) == (u[i] + 1) / C for i in range(12))
    assert lim not in leurres
# f(0) = 2/3, 3/4, 4/5.
assert [F(A, A + 1) for A in (2, 3, 4)] == [F(2, 3), F(3, 4), F(4, 5)]

# ar-l2-asymptotique-vitesse-ex3 (+ v1, v2) : plus petit n avec (sin 1)^n <= 10^-p ; leurres (cos 1, etc.).
for p, nmin, borne, leurres in [(3, 41, 40.02, [40, 12, 1000]), (2, 27, 26.68, [26, 8, 100]),
                                (4, 54, 53.36, [53, 15, 10000])]:
    q = -p * mp.log(10) / mp.log(mp.sin(1))
    assert abs(q - borne) < 0.01 and int(mp.ceil(q)) == nmin
    assert mp.sin(1) ** nmin <= mp.mpf(10) ** -p < mp.sin(1) ** (nmin - 1)
    assert int(mp.ceil(-p * mp.log(10) / mp.log(mp.cos(1)))) == leurres[1]
    assert nmin not in leurres
# Vérification de la majoration |u_n - l| <= (sin 1)^n.
u = it(lambda t: mp.cos(t), mp.mpf(0), 60)
assert all(abs(u[j] - l) <= mp.sin(1) ** j for j in range(61))
assert all(u[j] != l for j in range(61)) and 0 < mp.cos(1)

# ar-l2-asymptotique-vitesse-ex4 : Newton, f' = g g''/g'^2, f''(l) = g''(l)/g'(l).
G = symbols("G", cls=__import__("sympy").Function)
fN = x - G(x) / diff(G(x), x)
f1 = simplify(diff(fN, x) - G(x) * diff(G(x), x, 2) / diff(G(x), x) ** 2)
assert f1 == 0
f2 = diff(fN, x, 2)
s0, s1, s2, s3 = symbols("s0 s1 s2 s3")
f2s = f2.subs(diff(G(x), x, 3), s3).subs(diff(G(x), x, 2), s2).subs(diff(G(x), x), s1).subs(G(x), s0)
assert simplify(f2s - ((s1 * s2 + s0 * s3) / s1**2 - 2 * s0 * s2**2 / s1**3)) == 0
assert simplify(f2s.subs(s0, 0) - s2 / s1) == 0  # en un point où g = 0 : f''(l) = g''(l)/g'(l)
# Exemple g = x^2 - 2 : limite 1/(2 sqrt 2) ; les leurres diffèrent.
gg = x**2 - 2
lim = (diff(gg, x, 2) / (2 * diff(gg, x))).subs(x, sqrt(2))
assert simplify(lim - 1 / (2 * sqrt(2))) == 0
autres = [(diff(gg, x, 2) / diff(gg, x)).subs(x, sqrt(2)), -lim, (diff(gg, x) / (2 * diff(gg, x, 2))).subs(x, sqrt(2))]
assert all(simplify(o - lim) != 0 for o in autres)
u = it(lambda t: (t + 2 / t) / 2, mp.mpf(1), 4)
assert abs((u[3] - mp.sqrt(2)) / (u[2] - mp.sqrt(2)) ** 2 - 1 / (2 * mp.sqrt(2))) < 1e-2

# ar-l2-asymptotique-vitesse-ex5 : 1/sin^2 - 1/x^2 -> 1/3 ; u_n ~ sqrt(3/n) (cf. contre-exemple ci-dessus).
assert limit(1 / sin(x) ** 2 - 1 / x**2, x, 0) == Rational(1, 3)
# Ordre de grandeur de la remarque : sqrt(3/n) = 10^-3 pour n = 3·10^6.
assert abs(math.sqrt(3 / 3e6) - 1e-3) < 1e-12

# QCM q7 : f'(l) = -sin l pour f = cos.
assert diff(cos(x), x) == -sin(x)

# ---------------------------------------------------------------------------
# ar-l2-asymptotique-implicite
# ---------------------------------------------------------------------------
# Contre-exemples de la carte.
assert sorted(solve(x**3 - x / 16, x)) == [-Rational(1, 4), 0, Rational(1, 4)]
assert (x**3 - x / 16).subs(x, 1) == Rational(15, 16)
# Exemple tan x = x : x_n = n pi + arctan x_n ; développement à O(1/n^3).


def xt(N):
    t = N * mp.pi + mp.pi / 2
    for _ in range(100):
        t = N * mp.pi + mp.atan(t)
    return t


for N in [10, 100, 1000]:
    r = xt(N)
    assert abs(mp.tan(r) - r) < 1e-25 and N * mp.pi - mp.pi / 2 < r < N * mp.pi + mp.pi / 2
    d1 = r - (N * mp.pi + mp.pi / 2 - 1 / (N * mp.pi))
    assert abs(d1 * N**2 - 1 / (2 * mp.pi)) < 0.02  # ar-l2-asymptotique-implicite-ex4 : + 1/(2 pi n^2)
    d2 = d1 - 1 / (2 * mp.pi * N**2)
    assert abs(d2 * N**3) < 0.2  # reste O(1/n^3)
    for coef in [-1 / (2 * mp.pi), 1 / mp.pi, 1 / (2 * mp.pi**2)]:  # leurres de ex4
        assert abs(d1 - coef / N**2) > 0.5 * abs(d2) + 1e-3 / N**2
r10 = xt(10)
assert abs(r10 - (10 * mp.pi + mp.pi / 2 - 1 / (10 * mp.pi) + 1 / (200 * mp.pi))) < 1.2e-4
assert abs(r10 - (10 * mp.pi + mp.pi / 2 - 1 / (10 * mp.pi))) > 1.4e-3
# x + ln x = n : e^{x_n} = e^n / x_n ~ e^n / n.

# ar-l2-asymptotique-implicite-ex1 (+ v1, v2) : x^p + a n x - nb = 0.
for p, A, nb, eq in [(5, 1, 1, 1), (3, 1, 2, 2), (7, 2, 1, 0.5)]:
    for N in [10, 1000, 10**6]:
        P = lambda t: t**p + A * N * t - nb
        borne = F(nb, A * N)
        assert P(F(0)) < 0 < P(borne) and P(borne) == borne**p
        r = mp.findroot(lambda t: t**p + A * N * t - nb, mp.mpf(eq) / N)
        assert 0 < r < mp.mpf(borne.numerator) / borne.denominator
    rN = mp.findroot(lambda t: t**p + A * 10**6 * t - nb, mp.mpf(eq) / 10**6)
    assert abs(rN * 10**6 / mp.mpf(eq) - 1) < 1e-6
# Terme suivant de ex1 : x_n = 1/n - 1/n^6 + o(1/n^6).
N = mp.mpf(100)
r = mp.findroot(lambda t: t**5 + N * t - 1, 1 / N)
assert abs((r - 1 / N) * N**6 + 1) < 1e-6
# QCM q8 : x^2 + n x - 1 = 0, x_n = 1/n - 1/n^3 + o(1/n^3).
xs = 2 / (n + sqrt(n**2 + 4))
assert simplify(xs**2 + n * xs - 1) == 0
assert limit((xs - 1 / n) * n**3, n, oo) == -1

# ar-l2-asymptotique-implicite-ex2 (+ v1, v2).
for coef_x, coef_ln, dev, leurres, nmin in [
    (1, 1, lambda N: N - mp.log(N), [lambda N: N - mp.log(N) + 1, lambda N: N + mp.log(N), lambda N: N - mp.log(N) / 2], 1),
    (1, 2, lambda N: N - 2 * mp.log(N), [lambda N: N - mp.log(N), lambda N: N + 2 * mp.log(N),
                                        lambda N: N - 2 * mp.log(N) + 2], 1),
    (2, 1, lambda N: N / 2 - mp.log(N) / 2 + mp.log(2) / 2,
     [lambda N: N / 2 - mp.log(N) / 2, lambda N: N / 2 - mp.log(N), lambda N: N / 2 - mp.log(N) / 2 - mp.log(2) / 2], 2),
]:
    for N in [nmin, 10, 1000]:
        r = mp.findroot(lambda t: coef_x * t + coef_ln * mp.log(t) - N, max(mp.mpf(N) / coef_x, 1))
        assert 1 <= r + 1e-30 and r <= mp.mpf(N) / coef_x + 1e-30
    N = mp.mpf(10) ** 8
    r = mp.findroot(lambda t: coef_x * t + coef_ln * mp.log(t) - N, N / coef_x)
    assert abs(r - dev(N)) < 1e-5
    assert all(abs(r - l(N)) > 0.3 for l in leurres)
# « Pour aller plus loin » : x_n = n - ln n + ln n / n + o(ln n / n).
N = mp.mpf(10) ** 6
r = mp.findroot(lambda t: t + mp.log(t) - N, N)
assert abs((r - N + mp.log(N)) / (mp.log(N) / N) - 1) < 1e-4

# ar-l2-asymptotique-implicite-ex3 : x^n + x - 1, (x_n) croissante, 1 - x_n ~ ln n / n (lentement).


def xn(N):
    lo, hi = mp.mpf(0), mp.mpf(1)
    for _ in range(200):
        m = (lo + hi) / 2
        if m**N + m - 1 < 0:
            lo = m
        else:
            hi = m
    return lo


vals = [xn(N) for N in range(1, 40)]
assert all(vals[i] < vals[i + 1] for i in range(38)) and abs(vals[0] - 0.5) < 1e-30
for N in [10**4, 10**8, 10**12]:
    yN = 1 - xn(N)
    ratio = yN / (mp.log(N) / N)
    approx = (mp.log(N) - mp.log(mp.log(N))) / mp.log(N)
    assert abs(ratio - approx) < 0.05  # terme suivant -ln ln n / n
assert abs((1 - xn(10**8)) / (mp.log(10**8) / 10**8) - 0.85) < 0.01  # remarque « environ 0,85 »
# Leurres : 1/n (rapport -> +oo), ln n/(2n) (rapport -> 2), e^-n : faux.

# ---------------------------------------------------------------------------
# Ajouts de la relecture indépendante (2026-10-29).
# ---------------------------------------------------------------------------
# lecons/ar-l2-asymptotique (introduction) : un point répulsif peut être atteint, la suite est alors stationnaire :
# f(x) = 4x(1 - x), f'(0) = 4 > 1, u_0 = 1 donne u_1 = 0 puis 0, 0, ...
assert it(lambda t: 4 * t * (1 - t), F(1), 4) == [1, 0, 0, 0, 0] and diff(4 * x * (1 - x), x).subs(x, 0) == 4
# ar-l2-asymptotique-homographique (cours) : points fixes non réels => |k| = 1 et k != 1.
for A, B, Cc, D in [(0, -1, 1, 0), (1, -2, 1, 1), (2, -5, 1, -1), (3, -7, 2, 1)]:
    r = solve(Cc * x**2 + (D - A) * x - B, x)
    assert all(not t.is_real for t in r) and A * D - B * Cc != 0
    kk = (A - Cc * r[0]) / (A - Cc * r[1])
    assert abs(complex(kk.evalf())) - 1 < 1e-12 and abs(abs(complex(kk.evalf())) - 1) < 1e-12
    assert abs(complex(kk.evalf()) - 1) > 1e-6
# ar-l2-asymptotique-vitesse-ex3 : « environ 13 itérations par décimale » : ln 10 / (-ln sin 1) ≈ 13,3.
assert 13 < -mp.log(10) / mp.log(mp.sin(1)) < 13.5
# ar-l2-asymptotique-comparaison-ex4 (leurre « (v_n) bornée ») : suffisant, pas nécessaire.
assert limit((n + 1 / n) - n, n, oo) == 0 and limit((n + 1 / n) / n, n, oo) == 1  # e^{n+1/n} ~ e^n, v non bornée
assert limit((1 + 1 / n) * (1 + 1 / n**2) - (1 + 1 / n), n, oo) == 0  # v bornée, u ~ v => u - v -> 0
# ar-l2-asymptotique-homographique-ex4 (indice 2) : h(x) - g = (1 - g)(x - g)/(x + 1) pour g^2 = m.
assert simplify((x + gm**2) / (x + 1) - gm - (1 - gm) * (x - gm) / (x + 1)) == 0

print("ar_l2_asymptotique : toutes les vérifications passent.")
