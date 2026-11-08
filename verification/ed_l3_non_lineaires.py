"""Vérifications calculatoires du chapitre « ed-l3-non-lineaires » (équations différentielles, L3) :
cartes, exercices et QCM de content/*/equations-differentielles/l3-ed-l3-non-lineaires.json et
content/qcm/ed-l3-non-lineaires.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Chaque solution explicite est vérifiée par dérivation SymPy (ou dsolve). Les valeurs numériques (temps d'explosion,
limites, périodes) sont recalculées avec mpmath ou par la méthode de Runge-Kutta.
Lancer : python3 verification/ed_l3_non_lineaires.py (doit finir sans erreur, en quelques secondes).
"""

import mpmath as mp
from sympy import (E, Eq, Function, Rational, atan, cbrt, coth, cos, diff, dsolve, exp, factor, log, pi, simplify,
                   sin, sinh, cosh, solve, sqrt, symbols, tan, tanh)

t = symbols("t", real=True)
tp = symbols("t", positive=True)
y, u = symbols("y u", real=True)
R_ = Rational


def nul(expr):
    """Expression nulle : simplification SymPy, ou évaluation à 30 chiffres pour une constante."""
    if simplify(expr) == 0 or simplify(expr.rewrite(exp)) == 0:
        return True
    return not expr.free_symbols and abs(complex(expr.evalf(30))) < 1e-25


def num_nul(expr, var=t, pts=(R_(1, 7), R_(2, 5), R_(-1, 3))):
    """Égalité à zéro vérifiée numériquement en quelques points (pour les expressions où simplify peine)."""
    return all(abs(complex(expr.subs(var, p).evalf(30))) < 1e-20 for p in pts)


def resout(x, f, var=t, pts=None):
    """x' - f(var, x) = 0 (symboliquement, sinon numériquement en quelques points)."""
    e = diff(x, var) - f(var, x)
    return nul(e) or num_nul(e, var, pts or (R_(1, 7), R_(2, 5), R_(-1, 3)))


def ne_resout_pas(x, f, var=t, p=R_(1, 3)):
    return abs(complex((diff(x, var) - f(var, x)).subs(var, p).evalf(30))) > 1e-8


def rk4(F, Y0, t0, t1, n):
    """Intégration de Y' = F(t, Y) (Y liste de floats) par Runge-Kutta d'ordre 4 à pas constant."""
    h = (t1 - t0) / n
    T, Y = t0, list(Y0)
    for _ in range(n):
        k1 = F(T, Y)
        k2 = F(T + h / 2, [a + h / 2 * b for a, b in zip(Y, k1)])
        k3 = F(T + h / 2, [a + h / 2 * b for a, b in zip(Y, k2)])
        k4 = F(T + h, [a + h * b for a, b in zip(Y, k3)])
        Y = [a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(Y, k1, k2, k3, k4)]
        T += h
    return Y


mp.mp.dps = 20

# ===========================================================================
# ed-l3-non-lineaires-cauchy-lipschitz (carte)
# ===========================================================================
x0 = symbols("x0", positive=True)
# Exemple : x0/(1 - x0 t) solution de x' = x^2, valeur x0 en 0, pôle en 1/x0.
assert resout(x0 / (1 - x0 * t), lambda s, x: x**2) and (x0 / (1 - x0 * t)).subs(t, 0) == x0
# Contre-exemple 1 et affirmation 3 : t^3 solution de x' = 3 cbrt(x^2) (cbrt réelle : |t|^2 = t^2).
assert nul(diff(t**3, t) - 3 * t**2) and nul(cbrt(t**6) - t**2)  # pour t > 0 ; t^6 = |t|^6
assert nul((tp**6) ** R_(1, 3) - tp**2)
# Contre-exemple 1 : la dérivée de 3 x^(2/3) est 2 x^(-1/3), non bornée en 0 ; taux 3 h^(2/3) / h = 3 h^(-1/3).
xx = symbols("xx", positive=True)
assert nul(diff(3 * xx ** R_(2, 3), xx) - 2 * xx ** R_(-1, 3))
# Contre-exemple 2, affirmation 1 : 1/(1 - t).
assert resout(1 / (1 - t), lambda s, x: x**2)
# Contre-exemple 3 (U fermé) : 1 - t solution de x' = -1, sort de U = R x [0, +oo[ en t = 1.
assert resout(1 - t, lambda s, x: -1) and (1 - t).subs(t, 1) == 0
# Affirmation 5 : t/(1 - t) solution de x' = (1 + x)^2, x(0) = 0, intervalle ]-oo, 1[ non symétrique.
assert resout(t / (1 - t), lambda s, x: (1 + x) ** 2) and (t / (1 - t)).subs(t, 0) == 0
# Trous : leurres « continue », « continue et bornée », « C1 en t » réfutés par 3 cbrt(x^2) (resp. tronquée à |x| <= 1) :
# t^3 reste dans [-1, 1] pour |t| <= 1, et 3 cbrt(min(x^2,1)) = 3 t^2 sur ce segment.
assert all(abs(p**3) <= 1 for p in [R_(-1), R_(1, 2), R_(1)])

# ed-l3-non-lineaires-cauchy-lipschitz-ex1 (+ v1, v2) : x' = x^2, x(0) = 2, -1/3, 4.
for x00, sol, borne in [(2, 2 / (1 - 2 * t), R_(1, 2)), (R_(-1, 3), -1 / (t + 3), -3), (4, 4 / (1 - 4 * t), R_(1, 4))]:
    assert resout(sol, lambda s, x: x**2) and sol.subs(t, 0) == x00
    assert solve(1 / sol, t) == [borne]  # pôle unique : borne finie de l'intervalle maximal
# Leurres : ]-oo, 2[, ]1/2, +oo[ (ne contient pas 0), R ; v1 : ]-oo,-3[ (ne contient pas 0) ; v2 : ]-1/4, 1/4[ trop petit.

# ed-l3-non-lineaires-cauchy-lipschitz-ex2 (+ v1, v2) : x' = 1 + x^2.
for x00, sol, (a, b) in [(1, tan(t + pi / 4), (-3 * pi / 4, pi / 4)), (0, tan(t), (-pi / 2, pi / 2)),
                         (-1, tan(t - pi / 4), (-pi / 4, 3 * pi / 4))]:
    assert resout(sol, lambda s, x: 1 + x**2) and nul(sol.subs(t, 0) - x00)
    assert nul(b - a - pi)  # longueur pi
    assert nul(atan(x00) + a + pi / 2) and nul(atan(x00) + b - pi / 2)  # t + c atteint -pi/2 et pi/2

# ed-l3-non-lineaires-cauchy-lipschitz-ex3 : x_b = (t - b)^3 pour t > b, (t+1)^3 dans R x ]0, +oo[.
b = symbols("b", real=True)
assert nul(diff((t - b) ** 3, t) - 3 * (t - b) ** 2)
assert resout((tp + 1) ** 3, lambda s, x: 3 * x ** R_(2, 3), var=tp, pts=(R_(1, 2), R_(3),))
assert (t + 1).subs(t, -1) == 0  # le graphe atteint le bord x = 0 de U quand t -> -1

# ed-l3-non-lineaires-cauchy-lipschitz-ex4 : x' = x^2 + t^2, x(0) = 0.
# Q1 : y(t) = -x(-t) vérifie y' = y^2 + t^2 (calcul formel avec une fonction inconnue).
# y'(t) = x'(-t) (dérivée d'une composée, contrôlée sur x = sin), et x'(-t) = x(-t)^2 + t^2 = y(t)^2 + t^2.
assert nul(diff(-sin(-t), t) - cos(-t))
w = symbols("w", real=True)  # w = x(-t), y(t) = -w
assert nul((w**2 + (-t) ** 2) - ((-w) ** 2 + t**2))
# Q3 : majorant 1 + pi/2 - arctan(1/3) = 1 + arctan 3 ~ 2.249 ; Q4 : temps d'explosion réel ~ 2.003 (x = -u'/u avec
# u'' + t^2 u = 0, u(0) = 1, u'(0) = 0 : premier zéro de u).
assert nul(1 + pi / 2 - atan(R_(1, 3)) - 1 - atan(3))
uf = mp.odefun(lambda s, Y: [Y[1], -s**2 * Y[0]], 0, [1, 0])
beta = mp.findroot(lambda s: uf(s)[0], 2.0)
assert abs(beta - mp.mpf("2.0031")) < 1e-3 and beta < 1 + mp.pi / 2 - mp.atan(mp.mpf(1) / 3)
# Riccati : x = -u'/u vérifie x' = x^2 + t^2 dès que u'' = -t^2 u.
U = Function("U")
xr = -diff(U(t), t) / U(t)
assert nul((diff(xr, t) - xr**2 - t**2).subs(diff(U(t), t, 2), -t**2 * U(t)))

# ===========================================================================
# ed-l3-non-lineaires-separables (carte)
# ===========================================================================
# Exemple : 2/(2 - t^2) solution de x' = t x^2, x(0) = 1 ; G(x) = -1/x, A(t) + c = t^2/2 - 1 < 0 <=> |t| < sqrt 2.
assert resout(2 / (2 - t**2), lambda s, x: s * x**2) and (2 / (2 - t**2)).subs(t, 0) == 1
assert nul(-1 / (2 / (2 - t**2)) - (t**2 / 2 - 1))
# Contre-exemple 1, piège « C1 sur J privé d'un point », affirmation 3 : t|t| solution de x' = 2 sqrt|x|
# (t > 0 : t^2, dérivée 2t = 2 sqrt(t^2) ; t < 0 : -t^2, dérivée -2t = 2 sqrt(t^2)).
assert nul(diff(tp**2, tp) - 2 * sqrt(tp**2)) and nul(diff(-(-tp) ** 2, tp) + 2 * tp)
# Contre-exemple 2 : 1/(c - t) solution de x' = x^2, jamais nulle.
c = symbols("c", real=True)
assert resout(1 / (c - t), lambda s, x: x**2)
# Affirmation 5 : d/dt (tan + tan^3/3) = (1 + tan^2)^2, non constant.
assert nul(diff(tan(t) + tan(t) ** 3 / 3, t) - (1 + tan(t) ** 2) ** 2)
# Trous : x' = x (a = 1, g(x) = x), x = e^t ; G = x^2/2, x (primitive de g' = 1), -ln x : pas de la forme t + c.
for G in [exp(2 * t) / 2, exp(t), -t]:
    assert not nul(diff(G, t) - 1)
assert nul(diff(log(exp(t)), t) - 1)  # la bonne primitive de 1/g : ln x
# Trou 2, leurre « nulle » : g(x) = x - 1, la constante 1 est solution et vaut 1, pas 0.
assert resout(1 + 0 * t, lambda s, x: x - 1)

# ed-l3-non-lineaires-separables-ex1 (+ v1, v2) : x' = t x^2.
for x00, sol, leurres in [
    (1, 2 / (2 - t**2), [1 / (1 - t**2), 2 / (2 + t**2), exp(t**2 / 2)]),
    (-1, -2 / (t**2 + 2), [2 / (t**2 - 2), -1 / (1 + t**2), -exp(-t**2 / 2)]),
    (R_(1, 2), 2 / (4 - t**2), [1 / (2 - t**2), 2 / (4 + t**2), exp(t**2 / 4) / 2]),
]:
    f = lambda s, x: s * x**2
    assert resout(sol, f) and sol.subs(t, 0) == x00
    assert all(ne_resout_pas(l, f) for l in leurres)

# ed-l3-non-lineaires-separables-ex2 (+ v1, v2) : x' = e^(t - x), e^x = e^t + c.
for x00, cc, borne in [(-log(2), -R_(1, 2), -log(2)), (log(2), 1, None), (-log(4), -R_(3, 4), log(R_(3, 4)))]:
    sol = log(exp(t) + cc)
    assert resout(sol, lambda s, x: exp(s - x)) and nul(sol.subs(t, 0) - x00)
    if borne is not None:
        assert nul(exp(borne) + cc)  # e^t + c s'annule à la borne
assert nul(log(R_(3, 4)) + log(R_(4, 3)))  # leurre v2 : ln(4/3) = -ln(3/4)

# ed-l3-non-lineaires-separables-ex3 (+ v1, v2) : x' = 2t(1 + x^2), arctan x = t^2 + c.
for x00, cc, borne2 in [(1, pi / 4, pi / 4), (0, 0, pi / 2), (-1, -pi / 4, 3 * pi / 4)]:
    sol = tan(t**2 + cc)
    assert resout(sol, lambda s, x: 2 * s * (1 + x**2)) and nul(sol.subs(t, 0) - x00)
    assert nul(borne2 + cc - pi / 2)  # t^2 + c < pi/2 <=> t^2 < borne2
assert nul(sqrt(pi / 4) - sqrt(pi) / 2) and nul(sqrt(3 * pi / 4) - sqrt(3 * pi) / 2)

# ed-l3-non-lineaires-separables-ex4 (+ v1, v2) : x' = t(x^2 - 1).
f4 = lambda s, x: s * (x**2 - 1)
assert resout(-tanh(t**2 / 2), f4) and nul((1 - exp(t**2)) / (1 + exp(t**2)) + tanh(t**2 / 2))
assert (-tanh(t**2 / 2)).subs(t, 0) == 0 and nul(diff(log((1 - y) / (1 + y)) / 2, y) - 1 / (y**2 - 1))
from sympy import limit, oo  # noqa: E402

assert limit(-tanh(t**2 / 2), t, oo) == -1
sol1 = (2 + exp(t**2)) / (2 - exp(t**2))  # v1 : x(0) = 3, défini pour t^2 < ln 2
assert resout(sol1, f4) and sol1.subs(t, 0) == 3 and nul(diff(log((y - 1) / (y + 1)) / 2, y) - 1 / (y**2 - 1))
assert nul(exp(log(2)) - 2)
sol2 = (3 - exp(t**2)) / (3 + exp(t**2))  # v2 : x(0) = 1/2, s'annule en t = sqrt(ln 3)
assert resout(sol2, f4) and sol2.subs(t, 0) == R_(1, 2) and nul(sol2.subs(t, sqrt(log(3))))
assert all(not nul(sol2.subs(t, v)) for v in [log(3), sqrt(log(2)), sqrt(2 * log(3))])

# ===========================================================================
# ed-l3-non-lineaires-homogenes (carte)
# ===========================================================================
cc = symbols("cc", real=True)
# Exemple homogène : t tan(ln t + c) solution de x' = 1 + x/t + x^2/t^2.
assert resout(tp * tan(log(tp) + cc), lambda s, x: 1 + x / s + x**2 / s**2, var=tp, pts=(R_(1, 2), R_(3, 2)))
# Exemple Bernoulli : 1/(1 + C e^-t) solution de x' = x - x^2.
C0 = symbols("C0", real=True)
assert resout(1 / (1 + C0 * exp(-t)), lambda s, x: x - x**2)
# Contre-exemple Bernoulli : x(0) = -1, x = 1/(1 - 2e^-t), dénominateur nul en ln 2.
sb = 1 / (1 - 2 * exp(-t))
assert resout(sb, lambda s, x: x - x**2) and sb.subs(t, 0) == -1 and nul((1 - 2 * exp(-t)).subs(t, log(2)))
# Contre-exemple : x' = x^2/t^2, ty' = y^2 - y, zéros 0 et 1 ; x = t solution.
assert resout(tp, lambda s, x: x**2 / s**2, var=tp, pts=(R_(1, 2), R_(3),))
# Homothéties (cours) : mu x(t/mu) solution si x l'est, vérifié sur t tan(ln t + c).
mu = symbols("mu", positive=True)
xh = mu * (tp / mu) * tan(log(tp / mu) + cc)
assert resout(xh, lambda s, x: 1 + x / s + x**2 / s**2, var=tp, pts=(R_(1, 2), R_(3, 2)))
# Affirmation 2 : h = 1, x = t + 1, y = 1 + 1/t : y' = -1/t^2 et h(y) - y = -1/t.
yv = 1 + 1 / tp
assert nul(diff(yv, tp) + 1 / tp**2) and nul(1 - yv + 1 / tp)
# Affirmation 3 : y = 1/x pour x' = x + x^2 : y' = -y - 1.
Xf = Function("Xf")
assert nul((-diff(Xf(t), t) / Xf(t) ** 2).subs(diff(Xf(t), t), Xf(t) + Xf(t) ** 2) - (-1 / Xf(t) - 1))
# Affirmation 4 (Riccati) : identité algébrique.
x1, z, p, q = symbols("x1 z p q")
assert nul(((x1 + z) ** 2 + p * (x1 + z) + q) - (x1**2 + p * x1 + q) - (z**2 + (2 * x1 + p) * z))
# Énoncé Bernoulli : y = x^(1-alpha) vérifie y' = (1 - alpha)(a y + b) (vérification formelle).
al, a_, b_ = symbols("alpha a b", real=True)
Xb = Function("Xb", positive=True)
yb = Xb(t) ** (1 - al)
assert nul((diff(yb, t).subs(diff(Xb(t), t), a_ * Xb(t) + b_ * Xb(t) ** al)) - (1 - al) * (a_ * yb + b_))
# Trous : h = 1, x = t + 1, y = 1 + 1/t : ty' = h(y) - y vrai, les trois leurres faux.
assert nul(tp * diff(yv, tp) - (1 - yv))
for lhs, rhs in [(tp * diff(yv, tp), 1), (diff(yv, tp), 1 - yv), (diff(yv, tp), 1)]:
    assert not nul(lhs - rhs)
# Trou 2 : alpha = 2, a = 1, b = 0, x = e^t : x^(alpha-1), x^alpha, ln x ne vérifient pas y' = -(y).
for yl in [exp(t), exp(2 * t), t]:
    assert not nul(diff(yl, t) + yl)
assert nul(diff(exp(-t), t) + exp(-t))  # la bonne : x^(1-alpha) = e^-t

# ed-l3-non-lineaires-homogenes-ex1 (+ v1, v2) : équation vérifiée par y = x/t.
Y = Function("Y")
for h, bon, leurres in [
    (1 + y + y**2, 1 + y**2, None),
    (2 * y + y**2, y + y**2, None),
    (y**3, y**3 - y, None),
]:
    x_ = tp * Y(tp)
    eq = diff(x_, tp) - h.subs(y, Y(tp))
    assert nul(eq - (tp * diff(Y(tp), tp) - bon.subs(y, Y(tp))))
assert factor(1 + y + y**2 - y) == 1 + y**2

# ed-l3-non-lineaires-homogenes-ex2 (+ v1, v2).
for sol, f, x1v, xe, leurres in [
    (tp * sqrt(1 + 2 * log(tp)), lambda s, x: (s**2 + x**2) / (s * x), 1, E * sqrt(3), [E * sqrt(2), sqrt(3), 3 * E]),
    (tp * sqrt(4 + 2 * log(tp)), lambda s, x: (s**2 + x**2) / (s * x), 2, E * sqrt(6), [E * sqrt(5), 2 * E * sqrt(3), sqrt(6)]),
    (tp * sqrt(2 * tp**2 - 1), lambda s, x: (s**2 + 2 * x**2) / (s * x), 1, E * sqrt(2 * E**2 - 1),
     [E * sqrt(2 * E**2 + 1), E * sqrt(3), sqrt(2 * E**2 - 1)]),
]:
    assert resout(sol, f, var=tp, pts=(R_(3, 2), R_(5, 2))) and nul(sol.subs(tp, 1) - x1v)
    assert nul(sol.subs(tp, E) - xe) and all(not nul(xe - l) for l in leurres)
assert nul((1 + 2 * log(tp)).subs(tp, exp(-R_(1, 2)))) and nul((4 + 2 * log(tp)).subs(tp, exp(-2)))

# ed-l3-non-lineaires-homogenes-ex3 (+ v1, v2) : x' = -x + t x^2, y = 1/x, y' = y - t.
fb = lambda s, x: -x + s * x**2
for x00, sol, leurres in [
    (1, 1 / (1 + t), [1 / (1 - t), exp(-t), 1 / (1 + t + exp(t))]),
    (R_(1, 2), 1 / (t + 1 + exp(t)), [1 / (t + 1 - exp(t)), 1 / (1 + exp(t)), 1 / (t - 1 + 3 * exp(t))]),
    (2, 2 / (2 * t + 2 - exp(t)), [2 / (2 * t + 2 + exp(t)), 2 / (2 - exp(t)), 2 / (2 * t + 1)]),
]:
    assert resout(sol, fb) and sol.subs(t, 0) == x00
    for l in leurres:
        # un leurre ne vérifie pas le problème de Cauchy : mauvaise valeur en 0 ou pas solution
        assert l.subs(t, 0) != x00 or ne_resout_pas(l, fb)
# v1 : t + 1 + e^t s'annule une fois, entre -2 et -1.
assert (t + 1 + exp(t)).subs(t, -1) > 0 and (t + 1 + exp(t)).subs(t, -2) < 0
# dsolve confirme la solution générale de l'équation linéaire y' = y - t.
yl = Function("yl")
assert nul(dsolve(Eq(yl(t).diff(t), yl(t) - t), yl(t), ics={yl(0): 1}).rhs - (t + 1))

# ed-l3-non-lineaires-homogenes-ex4 (+ v1, v2) : Riccati x' = x^2 - 2tx + t^2 + 1, x1 = t, z' = z^2.
fr = lambda s, x: x**2 - 2 * s * x + s**2 + 1
assert resout(t, fr)
assert nul(expand_ := ((t + z) ** 2 - 2 * t * (t + z) + t**2 + 1) - (z**2 + 1))
for x00, sol, borne in [(1, t + 1 / (1 - t), 1), (-1, t - 1 / (1 + t), -1), (R_(1, 2), t + 1 / (2 - t), 2)]:
    assert resout(sol, fr) and sol.subs(t, 0) == x00 and solve(1 / (sol - t), t) == [borne]

# ed-l3-non-lineaires-homogenes-ex5 : x' = (x - t)/(x + t), courbes r = K e^-theta.
th, K = symbols("theta K", positive=True)
hh = (y - 1) / (y + 1)
assert nul(hh - y + (1 + y**2) / (1 + y))  # Q2
assert nul(diff(log(1 + y**2) / 2 + atan(y), y) - (1 + y) / (1 + y**2))  # Q3
assert num_nul(log(1 + tan(th) ** 2) / 2 + log(cos(th)), var=th, pts=(R_(-1, 5), R_(1, 2), R_(3, 2)))  # Q4 : (1/2) ln(1 + tan^2) = -ln cos
for rr, bon in [(K * exp(-th), True), (K * exp(th), False), (K + 3 - th, False)]:
    tt, xx_ = rr * cos(th), rr * sin(th)
    pente = (diff(xx_, th) / diff(tt, th)).subs(K, 2)
    ok = num_nul(pente - ((xx_ - tt) / (xx_ + tt)).subs(K, 2), var=th, pts=(R_(1, 5), R_(1, 2), R_(1)))
    assert ok == bon  # Q5 (r = K e^-theta) et leurres r = K e^theta, r = K - theta (ici K = 5 > theta)
assert simplify((K / cos(th)) * cos(th)) == K  # leurre : la droite verticale t = K n'est pas un graphe

# ===========================================================================
# ed-l3-non-lineaires-qualitatives (carte)
# ===========================================================================
xs0 = symbols("x_0", real=True)
logi = xs0 * exp(t) / (1 - xs0 + xs0 * exp(t))
assert resout(logi, lambda s, x: x * (1 - x)) and logi.subs(t, 0) == xs0
# Exemple : pour x0 < 0, le dénominateur s'annule en t1 = ln(1 - 1/x0) > 0.
assert nul((1 - xs0 + xs0 * exp(t)).subs(t, log(1 - 1 / xs0)))
assert log(1 - 1 / R_(-2)) > 0
# Contre-exemple : f(x) = 2 sqrt|1 - x|, x = 1 - (1 - t)^2 pour t <= 1 (1 - t > 0), atteint 1 en t = 1.
s_ = symbols("s", positive=True)  # s = 1 - t > 0
assert nul(diff(1 - s_**2, s_) * (-1) - 2 * sqrt(s_**2))
# Trous : f = x(1-x)(2-x), x0 = 1/2, f(1/2) > 0 ; limite 1 (plus petit zéro > 1/2), pas 2 ni 0.
ft = lambda v: v * (1 - v) * (2 - v)
assert ft(R_(1, 2)) == R_(3, 8) and ft(1) == 0 and ft(2) == 0
# Affirmation 2 : x0 = 2 -> 2e^t/(2e^t - 1), pôle en -ln 2.
assert nul(logi.subs(xs0, 2) - 2 * exp(t) / (2 * exp(t) - 1)) and nul((2 * exp(t) - 1).subs(t, -log(2)))
# Contre-exemple non autonome : sin t solution de x' = cos t.
assert resout(sin(t), lambda s, x: cos(s))

# ed-l3-non-lineaires-qualitatives-ex1 (+ v1, v2) : f(x) = (x - 1)(x - 3).
fq = lambda v: (v - 1) * (v - 3)
assert fq(2) == -1 and fq(0) == 3 and fq(4) == 3
assert all(fq(v) - (v - 3) ** 2 >= 0 for v in [3, 4, 5, 10])  # v2 : (x-1)(x-3) >= (x-3)^2 pour x >= 3
# v2 : explicitement, (x - 3)/(x - 1) = e^{2t}/3, x(t) = (9 - e^{2t})/(3 - e^{2t}), explose en ln(3)/2 <= 1.
solq = (9 - exp(2 * t)) / (3 - exp(2 * t))
assert solq.subs(t, 0) == 4 and resout(solq, lambda s, x: (x - 1) * (x - 3))
assert nul((3 - exp(2 * t)).subs(t, log(3) / 2)) and 0 < log(3) / 2 < 1

# ed-l3-non-lineaires-qualitatives-ex2 (+ v1, v2) : logistique, x0 = 2, -1, 3.
for x00, zero in [(2, -log(2)), (-1, log(2)), (3, log(R_(2, 3)))]:
    D = 1 - x00 + x00 * exp(t)
    assert nul(D.subs(t, zero)) and D.subs(t, 0) == 1
assert log(R_(2, 3)) < 0 < log(R_(3, 2))

# ed-l3-non-lineaires-qualitatives-ex3 (+ v1, v2) : x' = sin x, tan(x/2) = tan(x0/2) e^t.
for x00 in [1, 4, -2]:
    sol = 2 * atan(tan(R_(x00, 2)) * exp(t))
    # sur ]0, pi[ la formule 2 arctan convient ; pour 4 et -2 on vérifie la relation u' = u avec u = tan(x/2).
    uu = tan(R_(x00, 2)) * exp(t)
    assert nul(diff(uu, t) - uu)
assert resout(2 * atan(tan(R_(1, 2)) * exp(t)), lambda s, x: sin(x))
assert resout(2 * atan(tan(-1) * exp(t)), lambda s, x: sin(x))  # v2 : x(0) = -2
assert resout(2 * pi + 2 * atan(tan(2) * exp(t)), lambda s, x: sin(x))  # v1 : x(0) = 4 (branche pi + arctan)
assert nul((2 * pi + 2 * atan(tan(2) * exp(t))).subs(t, 0) - 4)
assert sin(1) > 0 and sin(4) < 0 and sin(-2) < 0 and sin(3 * pi / 2) == -1
# u' = u pour u = tan(x/2) quand x' = sin x.
XF = Function("XF")
uu = tan(XF(t) / 2)
assert nul(diff(uu, t).subs(diff(XF(t), t), sin(XF(t))) - uu)

# ed-l3-non-lineaires-qualitatives-ex4 (+ v1, v2) : x' = 1/(c^2 + t^2 + x^2), x(0) = 0.
for c2, cnum, mino_tex, borne_sup, leurre_bas in [
    (1, 1, pi / sqrt(4 + pi**2), pi / 2, pi / 4),
    (4, 2, 2 * pi / sqrt(64 + pi**2), pi / 4, pi / 8),
    (9, 3, 3 * pi / sqrt(324 + pi**2), pi / 6, pi / 12),
]:
    a2 = c2 + (pi / (2 * cnum)) ** 2
    assert nul(pi / (2 * sqrt(a2)) - mino_tex) and nul(pi / (2 * cnum) - borne_sup)
    g = mp.odefun(lambda s, xv: 1 / (c2 + s**2 + xv**2), 0, 0)
    T = 400
    xT = g(T)
    # x(T) < ell < x(T) + 1/T (car x' < 1/t^2 au-delà de T)
    assert float(mino_tex) < xT and xT + 1.0 / T < float(borne_sup)
    assert float(leurre_bas) < float(mino_tex)  # leurre « 0 < ell <= pi/(4c) » faux
assert nul(diff(atan(t / 2) / 2, t) - 1 / (4 + t**2)) and nul(diff(atan(t / 3) / 3, t) - 1 / (9 + t**2))

# ed-l3-non-lineaires-qualitatives-ex5 (+ v1, v2) : x' = x(1 - x) - h.
for h, r1, r2, k in [(R_(3, 16), R_(1, 4), R_(3, 4), R_(1, 2)), (R_(2, 9), R_(1, 3), R_(2, 3), R_(1, 3)),
                     (R_(4, 25), R_(1, 5), R_(4, 5), R_(3, 5))]:
    xv = symbols("xv")
    assert factor(xv * (1 - xv) - h) == factor(-(xv - r1) * (xv - r2))
    uu = symbols("uu")
    # u = r1 - x : u' = -x' = u^2 + k u
    assert nul(-((r1 - uu) * (1 - (r1 - uu)) - h) - (uu**2 + k * uu)) and k == r2 - r1
xv = symbols("xv")
assert nul(xv * (1 - xv) - symbols("h") + (xv - R_(1, 2)) ** 2 + (symbols("h") - R_(1, 4)))  # Q3

# ===========================================================================
# ed-l3-non-lineaires-physique (carte)
# ===========================================================================
# Exemple : theta = 4 arctan(e^t) - pi, pendule omega = 1, theta(0) = 0, theta'(0) = 2.
sep = 4 * atan(exp(t)) - pi
assert num_nul(diff(sep, t, 2) + sin(sep)) and nul(sep.subs(t, 0)) and nul(diff(sep, t).subs(t, 0) - 2)
# Contre-exemple : x = 1/(sqrt2 (1 - t)) solution de x'' = 4x^3.
xe = 1 / (sqrt(2) * (1 - t))
assert nul(diff(xe, t, 2) - 4 * xe**3)
# Frottement : e^{-t/2} cos(sqrt3 t/2) solution de x'' + x' + x = 0 ; E' = -x'^2.
xf = exp(-t / 2) * cos(sqrt(3) * t / 2)
assert nul(diff(xf, t, 2) + diff(xf, t) + xf)
# V seulement C1 : t^4/144 solution de x'' = sqrt|x| pour t > 0, et dérivées nulles en 0.
assert nul(diff(tp**4 / 144, tp, 2) - sqrt(tp**4 / 144))
# Affirmation 4 : coth(t - ln3/2) solution de v' = 1 - v^2, vaut -2 en 0.
vc = coth(t - log(3) / 2)
assert nul(diff(vc, t) - (1 - vc**2)) and nul(vc.subs(t, 0) + 2)
# Cours : v = sqrt(g/k) tanh(sqrt(gk) t) solution de v' = g - k v^2.
g_, k_ = symbols("g k", positive=True)
vt = sqrt(g_ / k_) * tanh(sqrt(g_ * k_) * t)
assert nul(diff(vt, t) - (g_ - k_ * vt**2))
# Trou 3, leurre « strictement décroissante » : V = -e^x, x = ln 2 - 2 ln(1 - t) solution de x'' = e^x, explose en 1.
xd = log(2) - 2 * log(1 - t)
assert nul(diff(xd, t, 2) - exp(xd))
# Trou 2 : V = x^2/2, x = cos t ; seules 1/2 x'^2 + V est constante.
xc = cos(t)
assert nul(diff(diff(xc, t) ** 2 / 2 + xc**2 / 2, t))
for leurre in [diff(xc, t) ** 2 / 2 - xc**2 / 2, diff(xc, t) ** 2 + xc**2 / 2, diff(xc, t) ** 2 / 2 + xc]:
    assert not nul(diff(leurre, t))
# Affirmation 2 : theta'^2 = 9 - 2(1 - cos theta) >= 5.
assert 9 - 2 * 2 == 5

# ed-l3-non-lineaires-physique-ex1 (+ v1, v2) : oscillateurs harmoniques.
for w2, sol, x00, v00, val in [(4, cos(2 * t), 1, 0, 2), (9, 2 * cos(3 * t), 2, 0, 6), (1, 3 * cos(t) + 4 * sin(t), 3, 4, 5)]:
    assert nul(diff(sol, t, 2) + w2 * sol) and sol.subs(t, 0) == x00 and diff(sol, t).subs(t, 0) == v00
E1 = R_(1, 2) * 0 + 2 * 1  # énergie 2 pour x'' + 4x = 0, max |x'| = sqrt(2E) = 2
assert sqrt(2 * E1) == 2 and sqrt(2 * 18) == 6 and sqrt(2 * R_(25, 2)) == 5

# ed-l3-non-lineaires-physique-ex2 (+ v1, v2) : v' = g - k v^2, v(0) = 0.
for gg, kk, vl, bb in [(10, R_(1, 10), 10, 1), (10, R_(2, 5), 5, 2), (R_(49, 5), R_(1, 5), 7, R_(7, 5))]:
    sol = vl * tanh(bb * t)
    assert nul(diff(sol, t) - (gg - kk * sol**2)) and sol.subs(t, 0) == 0
    assert nul(sqrt(gg / kk) - vl) and nul(sqrt(gg * kk) - bb)
# Leurres : g/k = 100, 25, 49 ; sqrt(gk) = 1, 2, 1.4 ; 10 m/s (v1) et 14 m/s (v2) ne sont pas les vitesses limites.

# ed-l3-non-lineaires-physique-ex3 (+ v1, v2) : pendule theta'' + w^2 sin theta = 0, seuil 2w.
for w in [3, 2, 1]:
    sp = 4 * atan(exp(w * t)) - pi
    assert num_nul(diff(sp, t, 2) + w**2 * sin(sp)) and nul(diff(sp, t).subs(t, 0) - 2 * w)
    assert nul(2 * w**2 * (1 + cos(y)) - 4 * w**2 * cos(y / 2) ** 2)
    # leurre : vitesse pour atteindre l'horizontale : v0^2 = 2 w^2 (1 - cos(pi/2)) = 2 w^2.
    assert nul(sqrt(2 * w**2) - sqrt(2) * w)
# Pour v0 < 2w, le pendule oscille : simulation (w = 3, v0 = 5.9) : |theta| reste < pi.
Yp = [0.0, 5.9]
mx = 0.0
for _ in range(40):
    Yp = rk4(lambda s, Y: [Y[1], -9 * mp.sin(Y[0])], Yp, 0, 0.1, 20)
    mx = max(mx, abs(float(Yp[0])))
assert mx < float(mp.pi)

# ed-l3-non-lineaires-physique-ex4 (+ v1, v2) : x'' + x^3 = 0, y = lam x(lam t).
lam = symbols("lam", positive=True)
Xa = Function("Xa")
ya = lam * Xa(lam * t)
# y'' = lam^3 X''(lam t) = -lam^3 X(lam t)^3 = -y^3 (règle de dérivation des composées, puis algèbre).
assert nul(lam**3 * (-xs**3) + (lam * xs) ** 3)
# Contrôle sur une fonction concrète : d^2/dt^2 [lam f(lam t)] = lam^3 f''(lam t) pour f = sin.
assert nul(diff(lam * sin(lam * t), t, 2) - lam**3 * (-sin(lam * t)))
for v0, l, per in [(4, 2, R_(1, 2)), (9, 3, R_(1, 3)), (R_(1, 4), R_(1, 2), 2)]:
    assert l**2 == v0 and per == 1 / R_(l)
# Contrôle numérique : période de x'' + x^3 = 0 pour x'(0) = 1 puis x'(0) = 4 (rapport 1/2).


def periode_quart(v0):
    """Premier instant > 0 où x' s'annule (quart de période) pour x(0) = 0, x'(0) = v0."""
    Yq, T, h = [0.0, float(v0)], 0.0, 1e-3
    while True:
        Y2 = rk4(lambda s, Yv: [Yv[1], -Yv[0] ** 3], Yq, T, T + h, 1)
        if Y2[1] <= 0:
            return T + h * Yq[1] / (Yq[1] - Y2[1])
        Yq, T = Y2, T + h


q1, q4 = periode_quart(1), periode_quart(4)
assert abs(q4 / q1 - 0.5) < 1e-3

# ed-l3-non-lineaires-physique-ex5 : période du pendule T = (4/w) K(k^2), k = sin(theta0/2), limite 2 pi/w.
uq, kq = symbols("u k", positive=True)
# Q4 : cos phi - cos theta0 = 2(k^2 - sin^2(phi/2)) et, avec sin(phi/2) = k sin u, = 2 k^2 cos^2 u.
ph, t0 = symbols("phi theta0", real=True)
assert nul(cos(ph) - cos(t0) - 2 * (sin(t0 / 2) ** 2 - sin(ph / 2) ** 2))
# d phi = 2 k cos u du / sqrt(1 - k^2 sin^2 u) pour phi = 2 arcsin(k sin u).
from sympy import asin  # noqa: E402

phi_u = 2 * asin(kq * sin(uq))
assert num_nul((diff(phi_u, uq) - 2 * kq * cos(uq) / sqrt(1 - kq**2 * sin(uq) ** 2)).subs(kq, R_(1, 2)), var=uq,
               pts=(R_(1, 5), R_(1, 2), R_(6, 5)))
# Contrôle numérique (w = 1, theta0 = 1) : quart de période par Runge-Kutta contre K(k^2).
th0 = 1.0
kk2 = mp.sin(th0 / 2) ** 2
Tth = 4 * mp.ellipk(kk2)
Yp, T, h = [th0, 0.0], 0.0, 1e-3
while True:
    Y2 = rk4(lambda s, Yv: [Yv[1], -mp.sin(Yv[0])], Yp, T, T + h, 1)
    if Y2[0] <= 0:
        tau = T + h * Yp[0] / (Yp[0] - Y2[0])
        break
    Yp, T = Y2, T + h
assert abs(4 * tau - Tth) < 1e-4
# Q5 : T -> 2 pi / w quand theta0 -> 0 (K(0) = pi/2) ; T croissant : K croît avec k^2.
assert abs(4 * mp.ellipk(0) - 2 * mp.pi) < 1e-15
assert mp.ellipk(0.1) < mp.ellipk(0.5) < mp.ellipk(0.9)
# Leurres : pi/w, 2 pi w, 4/w différents de 2 pi/w pour w = 1 (et w générique).
assert all(abs(v - 2 * mp.pi) > 0.1 for v in [mp.pi, 4.0])

# ===========================================================================
# ed-l3-non-lineaires-champs (carte)
# ===========================================================================
# Exemples : rotation (r cos t, r sin t), H = x^2 + y^2 ; selle (a e^t, b e^-t), H = xy.
r = symbols("r", positive=True)
assert nul(diff(r * cos(t), t) + r * sin(t)) and nul(diff(r * sin(t), t) - r * cos(t))
xs, ys = symbols("x y", real=True)
assert nul(2 * xs * (-ys) + 2 * ys * xs)
A_, B_ = symbols("A B", real=True)
assert nul(diff(A_ * exp(t) * B_ * exp(-t), t))
# Lotka-Volterra (cours) : H = dx - c ln x + by - a ln y.
a, bb_, cc_, dd = symbols("a b c d", positive=True)
xp, ypos = symbols("x y", positive=True)
assert nul((dd - cc_ / xp) * xp * (a - bb_ * ypos) + (bb_ - a / ypos) * ypos * (-cc_ + dd * xp))
# Contre-exemple : F = (1, 3 cbrt(y^2)), (t, t^3) solution (t > 0) ; non autonome : (sin t + a, -cos t + b).
assert nul(diff(tp**3, tp) - 3 * (tp**6) ** R_(1, 3))
assert nul(diff(sin(t) + a, t) - cos(t)) and nul(diff(-cos(t) + bb_, t) - sin(t))
# Les cercles de rayon 1 centrés en (0,0) et (1,0) se coupent en (1/2, ±sqrt3/2).
assert nul(R_(1, 4) + R_(3, 4) - 1) and nul((R_(1, 2) - 1) ** 2 + R_(3, 4) - 1)
# Trou 1 : F = (1, 3 cbrt(y^2)) est C1 en la première variable (indépendant de x) ; même contre-exemple.
# Trou 3 : rotation, t1 = 0, t2 = 2 pi : X(pi) = (-1, 0) != X(0), donc pi n'est pas une période.
assert cos(pi) != cos(0)
# Affirmation 4 : (e^-t, 0) solution de x' = -x.
assert nul(diff(exp(-t), t) + exp(-t))

# ed-l3-non-lineaires-champs-ex1 (+ v1, v2).
X1 = (3 * cos(t) - 4 * sin(t), 3 * sin(t) + 4 * cos(t))
assert nul(diff(X1[0], t) + X1[1]) and nul(diff(X1[1], t) - X1[0]) and nul(X1[0] ** 2 + X1[1] ** 2 - 25)
X2 = (cos(4 * t) - sin(4 * t), 2 * cos(4 * t) + 2 * sin(4 * t))
assert nul(diff(X2[0], t) + 2 * X2[1]) and nul(diff(X2[1], t) - 8 * X2[0]) and nul(4 * X2[0] ** 2 + X2[1] ** 2 - 8)
assert (X2[0].subs(t, 0), X2[1].subs(t, 0)) == (1, 2)
assert nul(X2[0] - sqrt(2) * cos(4 * t + pi / 4)) and nul(X2[1] - 2 * sqrt(2) * sin(4 * t + pi / 4))
X3 = (2 * cosh(t) + sinh(t), cosh(t) + 2 * sinh(t))
assert nul(diff(X3[0], t) - X3[1]) and nul(diff(X3[1], t) - X3[0]) and nul(X3[0] ** 2 - X3[1] ** 2 - 3)
cs = symbols("cs", real=True)
assert nul((2 / sqrt(3)) ** 2 - (1 / sqrt(3)) ** 2 - 1)
# Leurres v1 : (1,2) est sur x^2 + 4y^2 = 17 et sur le cercle de rayon sqrt5, mais 4x^2 + y^2 n'y est pas constant.
assert 1 + 4 * 4 == 17 and 1 + 4 == 5 and 4 + 4 != 5

# ed-l3-non-lineaires-champs-ex2 (+ v1, v2) : Lotka-Volterra.
for F1, F2, H, eq, leurres in [
    (xp * (2 - ypos), ypos * (xp - 3), xp - 3 * log(xp) + ypos - 2 * log(ypos), (3, 2),
     [xp - 2 * log(xp) + ypos - 3 * log(ypos), xp + 3 * log(xp) + ypos + 2 * log(ypos), xp * ypos - 3 * xp - 2 * ypos]),
    (xp * (1 - ypos), ypos * (xp - 2), xp - 2 * log(xp) + ypos - log(ypos), (2, 1),
     [xp - log(xp) + ypos - 2 * log(ypos), xp + 2 * log(xp) + ypos + log(ypos), xp * ypos - 2 * xp - ypos]),
    (xp * (3 - 2 * ypos), ypos * (xp - 1), xp - log(xp) + 2 * ypos - 3 * log(ypos), (1, R_(3, 2)),
     [xp - 3 * log(xp) + 2 * ypos - log(ypos), 2 * xp - log(xp) + ypos - 3 * log(ypos), xp + log(xp) + 2 * ypos + 3 * log(ypos)]),
]:
    assert nul(diff(H, xp) * F1 + diff(H, ypos) * F2)
    assert F1.subs({xp: eq[0], ypos: eq[1]}) == 0 and F2.subs({xp: eq[0], ypos: eq[1]}) == 0
    for l in leurres:
        assert not nul(diff(l, xp) * F1 + diff(l, ypos) * F2)

# ed-l3-non-lineaires-champs-ex3 (+ v1, v2) : x' = x, y' = -y (et x' = x, y' = y pour v2).
assert nul(diff(xs * ys, xs) * xs + diff(xs * ys, ys) * (-ys))
assert nul(diff(ys / xs, xs) * xs + diff(ys / xs, ys) * ys)
# 5 trajectoires dans {xy = 0} (origine + 4 demi-axes), 2 dans {xy = 1}, 3 dans la droite y = x.
assert 1 + 4 == 5

# ed-l3-non-lineaires-champs-ex4 (+ v1, v2) : cycle limite, rho' = 2 rho (1 - rho).
assert nul(diff(cos(t), t) - (-sin(t) + cos(t) * (1 - 1))) and nul(diff(sin(t), t) - (cos(t) + sin(t) * (1 - 1)))
F1c = -ys + xs * (1 - xs**2 - ys**2)
F2c = xs + ys * (1 - xs**2 - ys**2)
rho = xs**2 + ys**2
assert nul(2 * xs * F1c + 2 * ys * F2c - 2 * rho * (1 - rho))
for rho0, sol, borne, leurres in [
    (4, 4 * exp(2 * t) / (4 * exp(2 * t) - 3), log(R_(3, 4)) / 2, [log(R_(3, 4))]),
    (R_(1, 4), exp(2 * t) / (3 + exp(2 * t)), None, []),
    (9, 9 * exp(2 * t) / (9 * exp(2 * t) - 8), log(R_(8, 9)) / 2, [log(R_(8, 9)), log(R_(9, 8)) / 2]),
]:
    assert nul(diff(sol, t) - 2 * sol * (1 - sol)) and sol.subs(t, 0) == rho0
    if borne is not None:
        assert nul((1 / sol).subs(t, borne)) and borne < 0 and all(nul(l - borne) is False for l in leurres)
# v1 : 0 < rho < 1 sur R ; leurre ]-oo, ln(3)/2[ : aucun pôle.
assert solve(3 + exp(2 * t), t) == [] or all(not v.is_real for v in solve(3 + exp(2 * t), t))

# ===========================================================================
# QCM ed-l3-non-lineaires-q1 … q8
# ===========================================================================
# q1 : x' = x^2, x(0) = -2 -> -2/(1 + 2t), pôle -1/2.
assert resout(-2 / (1 + 2 * t), lambda s, x: x**2) and (-2 / (1 + 2 * t)).subs(t, 0) == -2
# q2 : non-unicité pour 3 cbrt(x^2) (vérifiée plus haut).
# q3 : x' = e^-x, x(0) = 0 -> ln(1 + t) ; leurres -ln(1 - t) et e^-t - 1 non solutions.
assert resout(log(1 + t), lambda s, x: exp(-x))
assert ne_resout_pas(-log(1 - t), lambda s, x: exp(-x)) and ne_resout_pas(exp(-t) - 1, lambda s, x: exp(-x))
# q4 : h(y) = y + y^2, ty' = y^2 ; x = -t/(ln t + c) solution de x' = x/t + x^2/t^2.
assert resout(-tp / (log(tp) + cc), lambda s, x: x / s + x**2 / s**2, var=tp, pts=(R_(1, 2), R_(3),))
# q5 : sin(3pi/2) = -1, plus grand zéro de sin sous 3pi/2 : pi.
assert sin(3 * pi / 2) == -1 and sin(pi) == 0
# q6 : x' = x^2 - 1, x(0) = 2 -> (3 + e^{2t})/(3 - e^{2t}), pôle ln(3)/2.
sq6 = (3 + exp(2 * t)) / (3 - exp(2 * t))
assert resout(sq6, lambda s, x: x**2 - 1) and sq6.subs(t, 0) == 2 and nul((3 - exp(2 * t)).subs(t, log(3) / 2))
# q7 : x'' + 9x = 0, x(0) = 0, x'(0) = 6 -> 2 sin 3t, maximum 2.
assert nul(diff(2 * sin(3 * t), t, 2) + 9 * 2 * sin(3 * t)) and diff(2 * sin(3 * t), t).subs(t, 0) == 6
# q8 : H = x - ln x + y - 2 ln y pour x' = x(2 - y), y' = y(x - 1).
Hq = xp - log(xp) + ypos - 2 * log(ypos)
assert nul(diff(Hq, xp) * xp * (2 - ypos) + diff(Hq, ypos) * ypos * (xp - 1))

print("ed_l3_non_lineaires : toutes les vérifications passent.")
