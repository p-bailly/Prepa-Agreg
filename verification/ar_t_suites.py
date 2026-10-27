"""Vérifications calculatoires du chapitre « ar-t-suites » (analyse réelle, terminale) :
cartes et exercices de content/*/analyse-reelle/terminale-ar-t-suites.json.

Chaque assert cite l'identifiant de la carte ou de l'exercice (variante comprise : v1, v2).
Lancer : python3 verification/ar_t_suites.py (doit finir sans erreur).
"""

import math
from fractions import Fraction as F

from sympy import Rational, factorial, isprime, limit, log, oo, simplify, sqrt, symbols

n = symbols("n", positive=True)
x = symbols("x", real=True)


def suite(u0, f, N):
    """Premiers termes u_0..u_N de u_{n+1} = f(u_n), en fractions exactes."""
    u = [F(u0)]
    for _ in range(N):
        u.append(f(u[-1]))
    return u


def suite2(u0, u1, f, N):
    """Premiers termes de u_{n+2} = f(u_{n+1}, u_n)."""
    u = [F(u0), F(u1)]
    while len(u) <= N:
        u.append(f(u[-1], u[-2]))
    return u


# ---------------------------------------------------------------------------
# ar-t-suites-recurrence (carte)
# ---------------------------------------------------------------------------
# Exemple : 3 divise 4^n - 1 ; hérédité 4^{n+1} - 1 = 4(4^n - 1) + 3.
assert all((4**k - 1) % 3 == 0 for k in range(200))
# Contre-exemple / affirmation 1 : 4^n + 1 est héréditaire (4^{n+1}+1 = 4(4^n+1) - 3) mais jamais divisible par 3.
assert all(4 ** (k + 1) + 1 == 4 * (4**k + 1) - 3 for k in range(50))
assert all((4**k + 1) % 3 == 2 for k in range(200))
# Contre-exemple / piège 1 : n^2 - n + 41 premier pour n = 0..40, pas pour n = 41 (= 41^2).
assert all(isprime(k * k - k + 41) for k in range(41))
assert 41 * 41 - 41 + 41 == 41**2 and not isprime(41**2)
# Contre-exemple / affirmation 3 / piège 3 : P(n) = « n = 0 ou n >= 2 » : P(0), hérédité pour n >= 1, P(1) fausse.
P = lambda k: k == 0 or k >= 2
assert P(0) and not P(1) and all((not P(k)) or P(k + 1) for k in range(1, 100))
# Piège 2 (« il existe n ») et piège 4 (P(n+1) => P(n)) : P(n) = « n = 0 » : P(0) vraie, P(1) fausse,
# P(1) => P(2) vraie (prémisse fausse) ; P(n+1) => P(n) vraie pour tout n (prémisse toujours fausse).
Q = lambda k: k == 0
assert Q(0) and not Q(1) and ((not Q(1)) or Q(2)) and all((not Q(k + 1)) or Q(k) for k in range(100))
# Trous : leurres « pour une infinité d'entiers » : P(n) = « n = 0 », P(n) => P(n+1) vrai pour tout n >= 1 (infinité).
assert all((not Q(k)) or Q(k + 1) for k in range(1, 100)) and not Q(1)
# Leurre « P(n0+1) » / « P(n0) ou P(n0+1) » / « P(n1), n1 > n0 » : R(n) = « n >= 1 » (n0 = 0) est héréditaire,
# vraie en 1 mais fausse en 0.
R = lambda k: k >= 1
assert R(1) and not R(0) and all((not R(k)) or R(k + 1) for k in range(100))
# Affirmation 2 (récurrence double) : illustration sur Fibonacci F_n <= 2^n.
fib = suite2(0, 1, lambda a, b: a + b, 60)
assert all(fib[k] <= 2**k for k in range(61))

# ---------------------------------------------------------------------------
# ar-t-suites-recurrence-ex1 : 4 divise 5^n - 1 ; 5^{n+1} - 1 = 5(5^n - 1) + 4.
# ---------------------------------------------------------------------------
assert all((5**k - 1) % 4 == 0 for k in range(200))
assert all(5 ** (k + 1) - 1 == 5 * (5**k - 1) + 4 for k in range(50))

# ---------------------------------------------------------------------------
# ar-t-suites-recurrence-ex2 (+ v1, v2) : conjecture de u_n ; leurres faux sur un terme calculé.
# ---------------------------------------------------------------------------
def verifier_formule(u, bonne, leurres, termes_cites):
    assert all(u[k] == bonne(k) for k in range(len(u)))
    for le in leurres:
        assert any(u[k] != le(k) for k in range(4))
    for k, val in termes_cites.items():
        assert u[k] == val


# ex2 : u0 = 2, u_{n+1} = 2u_n - 1 : u_n = 2^n + 1 ; leurres n+2, 2^{n+1}, 3^n + 1.
u = suite(2, lambda a: 2 * a - 1, 30)
verifier_formule(u, lambda k: 2**k + 1, [lambda k: k + 2, lambda k: 2 ** (k + 1), lambda k: 3**k + 1],
                 {1: 3, 2: 5, 3: 9})
assert 2 + 2 == 4 and 2**2 == 4 and 3 + 1 == 4  # valeurs « 4 » citées pour n+2 (n=2), 2^{n+1}, 3^n+1 (n=1)
# ex2 v1 : u0 = 3, u_{n+1} = 3u_n - 4 : u_n = 3^n + 2 ; leurres 2n+3, 3^{n+1}, 2*3^n + 1.
u = suite(3, lambda a: 3 * a - 4, 30)
verifier_formule(u, lambda k: 3**k + 2, [lambda k: 2 * k + 3, lambda k: 3 ** (k + 1), lambda k: 2 * 3**k + 1],
                 {1: 5, 2: 11, 3: 29})
assert 2 * 2 + 3 == 7 and 3**2 == 9 and 2 * 3 + 1 == 7
# ex2 v2 : u0 = 0, u_{n+1} = 2u_n + 3 : u_n = 3*2^n - 3 ; leurres 3n, 3*2^{n+1} - 6, 2^{n+1} + n - 2.
u = suite(0, lambda a: 2 * a + 3, 30)
verifier_formule(u, lambda k: 3 * 2**k - 3, [lambda k: 3 * k, lambda k: 3 * 2 ** (k + 1) - 6,
                                            lambda k: 2 ** (k + 1) + k - 2], {1: 3, 2: 9, 3: 21})
assert 3 * 2 == 6 and 3 * 4 - 6 == 6 and 8 + 2 - 2 == 8

# ---------------------------------------------------------------------------
# ar-t-suites-recurrence-ex3 (+ v1, v2) : plus petit n0 tel que l'inégalité vaut pour tout n >= n0.
# ---------------------------------------------------------------------------
def plus_petit_n0(ok, N=400):
    """Plus petit n0 tel que ok(n) pour tout n0 <= n <= N (N grand, l'hérédité étant prouvée à part)."""
    n0 = N
    while n0 > 0 and ok(n0 - 1):
        n0 -= 1
    return n0


# ex3 : 2^n > n^2 -> n0 = 5 ; leurres 4, 0, 3 faux.
n0 = plus_petit_n0(lambda k: 2**k > k**2)
assert n0 == 5 and all(2**k > k**2 for k in range(5, 400))
assert not (2**4 > 4**2) and not (2**3 > 3**2) and (2**0 > 0) and (2**1 > 1)
# Hérédité : 2n^2 - (n+1)^2 = (n-1)^2 - 2 >= 14 > 0 pour n >= 5.
assert simplify(2 * x**2 - (x + 1) ** 2 - ((x - 1) ** 2 - 2)) == 0 and (5 - 1) ** 2 - 2 == 14
# ex3 v1 : 3^n > n^3 -> n0 = 4 ; leurres 3 (3^3 = 27 non > 27), 0, 5.
n0 = plus_petit_n0(lambda k: 3**k > k**3)
assert n0 == 4 and not (3**3 > 3**3) and 3**2 > 2**3
assert F(5, 4) ** 3 == F(125, 64) and F(125, 64) < 3  # (1 + 1/n)^3 <= (5/4)^3 < 3 pour n >= 4
# ex3 v2 : 2^n > n^3 -> n0 = 10 ; 2^9 = 512 < 729 ; leurres 9, 0, 11.
n0 = plus_petit_n0(lambda k: 2**k > k**3)
assert n0 == 10 and 2**9 == 512 and 9**3 == 729 and 2**10 == 1024
assert all(not (2**k > k**3) for k in range(2, 10))
assert F(11, 10) ** 3 == F(1331, 1000) and F(1331, 1000) < 2

# ---------------------------------------------------------------------------
# ar-t-suites-recurrence-ex4 (+ v1, v2) : récurrences d'ordre 2 ; leurres qui collent aux premiers termes.
# ---------------------------------------------------------------------------
# ex4 : u0=1, u1=3, u_{n+2} = 3u_{n+1} - 2u_n -> 2^{n+1} - 1 ; u2,u3,u4 = 7,15,31.
u = suite2(1, 3, lambda a, b: 3 * a - 2 * b, 40)
assert u[2:5] == [7, 15, 31] and all(u[k] == 2 ** (k + 1) - 1 for k in range(41))
assert [3**2, 2 * 2 + 1, 2**2 + 2 + 1] == [9, 5, 7] and [3**3, 2 * 3 + 1, 3**2 + 3 + 1] == [27, 7, 13]
for le in (lambda k: 3**k, lambda k: 2 * k + 1, lambda k: k * k + k + 1):
    assert le(0) == u[0] and le(1) == u[1] and any(le(k) != u[k] for k in range(5))
# Hérédité : 3(2^{n+2}-1) - 2(2^{n+1}-1) = 2^{n+3} - 1.
assert simplify(3 * (2 ** (n + 2) - 1) - 2 * (2 ** (n + 1) - 1) - (2 ** (n + 3) - 1)) == 0
# ex4 v1 : u0=2, u1=5, u_{n+2} = 5u_{n+1} - 6u_n -> 2^n + 3^n ; u2,u3,u4 = 13,35,97.
u = suite2(2, 5, lambda a, b: 5 * a - 6 * b, 40)
assert u[2:5] == [13, 35, 97] and all(u[k] == 2**k + 3**k for k in range(41))
assert [3 * 2 + 2, 4**2 + 1, F(5 * 4 + 2 + 4, 2)] == [8, 17, 13] and F(5 * 9 + 3 + 4, 2) == 26
for le in (lambda k: 3 * k + 2, lambda k: 4**k + 1, lambda k: F(5 * k * k + k + 4, 2)):
    assert le(0) == u[0] and le(1) == u[1] and any(le(k) != u[k] for k in range(5))
assert simplify(5 * (2 ** (n + 1) + 3 ** (n + 1)) - 6 * (2**n + 3**n) - (2 ** (n + 2) + 3 ** (n + 2))) == 0
# ex4 v2 : u0=0, u1=1, u_{n+2} = 4u_{n+1} - 3u_n -> (3^n - 1)/2 ; u2,u3,u4 = 4,13,40.
u = suite2(0, 1, lambda a, b: 4 * a - 3 * b, 40)
assert u[2:5] == [4, 13, 40] and all(u[k] == F(3**k - 1, 2) for k in range(41))
assert [2**2, 2**2 - 1, F(2 * 5, 2)] == [4, 3, 5] and 3**2 == 9
for le in (lambda k: k * k, lambda k: 2**k - 1, lambda k: F(k * (3 * k - 1), 2)):
    assert le(0) == u[0] and le(1) == u[1] and any(le(k) != u[k] for k in range(5))
assert simplify((4 * (3 ** (n + 1) - 1) - 3 * (3**n - 1)) / 2 - (3 ** (n + 2) - 1) / 2) == 0

# ---------------------------------------------------------------------------
# ar-t-suites-recurrence-ex5 : P(n) « n réels quelconques sont égaux » : P(1) vraie, P(n) fausse pour n >= 2,
# donc P(n) => P(n+1) fausse seulement pour n = 1.
# ---------------------------------------------------------------------------
Pn = lambda k: k <= 1  # valeur de vérité de P(k) pour k >= 1
assert [k for k in range(1, 100) if Pn(k) and not Pn(k + 1)] == [1]

# ---------------------------------------------------------------------------
# ar-t-suites-geometrique (carte)
# ---------------------------------------------------------------------------
# Démonstration : (1+a)^n >= 1 + n a (Bernoulli), contrôle numérique.
assert all((1 + a) ** k >= 1 + k * a for a in (F(0), F(1, 10), F(1), F(3)) for k in range(40))
# Contre-exemple / affirmation 1 / piège 2 : q_n = 2^{-1/n} dans ]0,1[ et q_n^n = 1/2.
for k in range(1, 50):
    q = 2 ** (-1 / k)
    assert 0 < q < 1 and math.isclose(q**k, 0.5, rel_tol=1e-12)
# Contre-exemple / affirmation 2 / pièges 1 et 4 : (-2)^n : pairs >= 1, impairs <= -2.
assert all((-2) ** (2 * k) >= 1 and (-2) ** (2 * k + 1) <= -2 for k in range(30))
# Démonstration cas q <= -1 : q^{2n} >= 1 et q^{2n+1} <= -1 (contrôle sur quelques q).
for q in (F(-1), F(-3, 2), F(-2), F(-7)):
    assert all(q ** (2 * k) >= 1 and q ** (2 * k + 1) <= -1 for k in range(20))
# Exemple / affirmation 3 : somme géométrique -> 1/(1-q) ; 1 + 1/2 + ... -> 2.
q = symbols("q")
assert simplify(sum(q**k for k in range(6)) - (1 - q**6) / (1 - q)) == 0
assert limit((1 - Rational(1, 2) ** (n + 1)) / (1 - Rational(1, 2)), n, oo) == 2
# Trous : leurres « converge vers q » (q = 1/2), « 1/(1-q) », « |q^n| -> +oo » faux pour q = -1.
assert limit(Rational(1, 2) ** n, n, oo) == 0 and Rational(1, 2) != 0 and 1 / (1 - Rational(1, 2)) != 0
assert all(abs((-1) ** k) == 1 for k in range(10))

# ---------------------------------------------------------------------------
# ar-t-suites-geometrique-ex1 (+ v1, v2)
# ---------------------------------------------------------------------------
# |u_n - L| = c |q|^n -> 0 (base négative : on passe par la valeur absolue).
assert limit(5 * Rational(3, 4) ** n, n, oo) == 0 and 2 - 5 == -3  # ex1 : 2 ; leurre -3 (= u_0)
assert all(abs(2 - 5 * F(-3, 4) ** k - 2) == 5 * F(3, 4) ** k for k in range(30))
assert limit(4 * Rational(2, 5) ** n, n, oo) == 0 and 1 + 4 == 5  # v1 : 1 ; leurre 5 (= u_0)
# v2 : q = -6/5 <= -1 : pas de limite ; termes pairs de u_n = 3 - q^n -> -oo, impairs -> +oo.
qq = F(-6, 5)
assert qq <= -1 and 3 - qq**40 < -1000 and 3 - qq**41 > 1000

# ---------------------------------------------------------------------------
# ar-t-suites-geometrique-ex2 (+ v1, v2) : seuils.
# ---------------------------------------------------------------------------
def seuil(q, eps):
    k = 0
    while not (q**k < eps):
        k += 1
    return k


assert seuil(F(8, 10), F(1, 1000)) == 31  # ex2 : 31 (leurres 30, 32, 29 faux)
assert math.isclose(-3 * math.log(10) / math.log(0.8), 30.957, abs_tol=1e-3)
assert math.isclose(0.8**30, 1.238e-3, rel_tol=1e-3) and math.isclose(0.8**31, 9.904e-4, rel_tol=1e-3)
assert seuil(F(9, 10), F(1, 100)) == 44  # v1 : 44 (leurres 43, 45, 42)
assert math.isclose(-2 * math.log(10) / math.log(0.9), 43.709, abs_tol=1e-3)
assert math.isclose(0.9**43, 1.0775e-2, rel_tol=1e-3) and math.isclose(0.9**44, 9.698e-3, rel_tol=1e-3)
assert seuil(F(1, 2), F(1, 10000)) == 14  # v2 : 14 (leurres 13, 15, 12)
assert math.isclose(4 * math.log(10) / math.log(2), 13.288, abs_tol=1e-3)
assert math.isclose(0.5**13, 1.2207e-4, rel_tol=1e-3) and math.isclose(0.5**14, 6.1035e-5, rel_tol=1e-3)

# ---------------------------------------------------------------------------
# ar-t-suites-geometrique-ex3 (+ v1, v2) : arithmético-géométriques.
# ---------------------------------------------------------------------------
def arith_geo(u0, a, b, c, N=60):
    u = suite(u0, lambda t: a * t + b, N)
    v = [t - c for t in u]
    assert all(v[k + 1] == a * v[k] for k in range(N))  # (v_n) géométrique de raison a
    return u


u = arith_geo(1, F(6, 10), 2, 5)  # ex3 : c = 5, u_n = 5 - 4*0,6^n -> 5
assert all(u[k] == 5 - 4 * F(6, 10) ** k for k in range(61))
assert limit(5 - 4 * Rational(3, 5) ** n, n, oo) == 5 and F(2) / F(16, 10) == F(5, 4)  # leurre 5/4
u = arith_geo(0, F(1, 2), 3, 6)  # v1 : c = 6 -> 6 ; leurre 2 = 3/1,5
assert all(u[k] == 6 - 6 * F(1, 2) ** k for k in range(61)) and F(3) / F(3, 2) == 2
u = arith_geo(1, F(3, 2), -1, 2)  # v2 : c = 2, u_n = 2 - 1,5^n -> -oo ; leurre -2/5 = -1/2,5
assert all(u[k] == 2 - F(3, 2) ** k for k in range(61)) and u[60] < -10**9 and F(-1) / F(5, 2) == F(-2, 5)
assert limit(2 - Rational(3, 2) ** n, n, oo) == -oo

# ---------------------------------------------------------------------------
# ar-t-suites-geometrique-ex4 (+ v1, v2) : raison dépendant de x.
# ---------------------------------------------------------------------------
def nature(qv):
    """'0', '1' (convergentes) ou 'div' selon la valeur de q."""
    if -1 < qv < 1:
        return "0"
    if qv == 1:
        return "1"
    return "div"


q1 = lambda t: F(2 * t, t * t + 1)
assert simplify(x**2 + 1 - 2 * abs(x) - (abs(x) - 1) ** 2) == 0
assert nature(q1(1)) == "1" and nature(q1(-1)) == "div"
pts = [F(k, 7) for k in range(-70, 71)]
assert all(nature(q1(t)) != "div" for t in pts if t != -1)  # bonne réponse R \ {-1}
assert nature(q1(1)) != "div"  # leurre R \ {-1, 1} et ]-1,1[ : x = 1 (et x = 2) convergent aussi
assert nature(q1(F(2))) == "0"
q2 = lambda t: F(t * t - 1, t * t + 1)
assert nature(q2(0)) == "div" and all(nature(q2(t)) == "0" for t in pts if t != 0)  # v1 : R \ {0}
q3 = lambda t: F(1 - t, 1 + t)
assert nature(q3(0)) == "1" and all(nature(q3(t)) == "0" for t in pts if t > 0)  # v2 : [0, +oo[
assert all(nature(q3(t)) == "div" for t in pts if t < 0 and t != -1)  # x < 0 diverge (]-1,0[ : q > 1)
assert all(q3(t) > 1 for t in pts if -1 < t < 0) and all(q3(t) < -1 for t in pts if t < -1)

# ---------------------------------------------------------------------------
# ar-t-suites-geometrique-ex5 : (1+a)^n >= 1 + na + n(n-1)/2 a^2 ; n q^n -> 0.
# ---------------------------------------------------------------------------
a = symbols("a", nonnegative=True)
k_ = symbols("k", nonnegative=True)
dev = ((1 + k_ * a + k_ * (k_ - 1) / 2 * a**2) * (1 + a)).expand()
cible = (1 + (k_ + 1) * a + (k_ + 1) * k_ / 2 * a**2).expand()
assert simplify(dev - cible - k_ * (k_ - 1) / 2 * a**3) == 0  # hérédité : reste = n(n-1)/2 a^3 >= 0
assert all((1 + b) ** k >= 1 + k * b + F(k * (k - 1), 2) * b * b for b in (F(1, 3), F(2)) for k in range(30))
assert limit(n * Rational(9, 10) ** n, n, oo) == 0 and limit(2 / ((n - 1) * Rational(1, 4)), n, oo) == 0

# ---------------------------------------------------------------------------
# ar-t-suites-operations (carte)
# ---------------------------------------------------------------------------
assert limit((3 * n**2 - n) / (2 * n**2 + 5), n, oo) == Rational(3, 2)  # exemple
assert limit(sqrt(n + 1) - sqrt(n), n, oo) == 0
assert limit(n**2 - n, n, oo) == oo and limit(n - n**2, n, oo) == -oo  # contre-exemples oo - oo
assert limit((1 / n) * n**2, n, oo) == oo and limit((1 / n**2) * n, n, oo) == 0  # 0 x oo
assert limit(n**2 / n, n, oo) == oo  # affirmation 1 : n^2 / n -> +oo (pas 1)
# Signe constant : 1 / ((-1)^n/(n+1)) = (-1)^n (n+1) : valeurs paires > 0 grandes, impaires < 0.
assert all(F(1) / (F((-1) ** k, k + 1)) == (-1) ** k * (k + 1) for k in range(50))
# Démonstration du quotient : u >= l/2 et 0 < v < l/(2A) donnent u/v > A (contrôle numérique).
for l_, A in ((F(1), F(10)), (F(3), F(1000))):
    u_, v_ = l_ / 2, l_ / (2 * A) * F(99, 100)
    assert u_ / v_ > A

# ---------------------------------------------------------------------------
# ar-t-suites-operations-ex1..ex4 (+ v1, v2) : limites et leurres.
# ---------------------------------------------------------------------------
def lim_et_leurres(expr, bonne, leurres):
    L = limit(expr, n, oo)
    assert simplify(L - bonne) == 0 if bonne not in (oo, -oo) else L == bonne
    for le in leurres:
        assert L != le


lim_et_leurres((2 * n**2 - 3 * n + 1) / (5 * n**2 + 4), Rational(2, 5), [Rational(1, 4), oo, 0])  # ex1
lim_et_leurres((4 * n**3 + n) / (1 - 2 * n**3), -2, [2, -oo, 0])  # ex1 v1
lim_et_leurres((3 * n**2 + 1) / (n**3 - 2), 0, [3, oo, Rational(-1, 2)])  # ex1 v2
lim_et_leurres(sqrt(n**2 + 4 * n) - n, 2, [0, 4, oo])  # ex2
lim_et_leurres(sqrt(n**2 + 6 * n + 1) - n, 3, [0, 6, 1])  # ex2 v1
lim_et_leurres(sqrt(4 * n**2 + n) - 2 * n, Rational(1, 4), [0, Rational(1, 2), 1])  # ex2 v2
lim_et_leurres((2 * n + 1) / (n - sqrt(n**2 + 1)), -oo, [oo, 0, 2])  # ex3
assert simplify((n - sqrt(n**2 + 1)) * (n + sqrt(n**2 + 1)) + 1) == 0  # ex3 : conjuguée
lim_et_leurres((n - 3) / (sqrt(n + 1) - sqrt(n)), oo, [-oo, 0, 1])  # ex3 v1
lim_et_leurres((1 - 2 * n) / (sqrt(n**2 + 3) - n), -oo, [oo, 0, -2])  # ex3 v2
assert simplify((sqrt(n**2 + 3) - n) * (sqrt(n**2 + 3) + n) - 3) == 0
# ex4 (+ v1, v2) : u_n = sqrt(...) - a n ; convergence seulement pour a = a0, vers L ; leurres.
for rad, a0, L in ((n**2 + 1, 1, 0), (4 * n**2 + n, 2, Rational(1, 4)), (n**2 + 6 * n, 1, 3)):
    assert limit(sqrt(rad) - a0 * n, n, oo) == L
    for av in (a0 - 1, a0 - Rational(1, 2), -1, 0):  # a < a0 : +oo
        assert limit(sqrt(rad) - av * n, n, oo) == oo
    for av in (a0 + Rational(1, 2), a0 + 1, 4):  # a > a0 : -oo
        if av != a0:
            assert limit(sqrt(rad) - av * n, n, oo) == -oo
assert limit(sqrt(n**2 + 1) - 1 * n, n, oo) != 1  # ex4 leurre « vers 1 »
assert limit(sqrt(4 * n**2 + n) - 2 * n, n, oo) != 0  # ex4 v1 leurre « vers 0 »
assert limit(sqrt(4 * n**2 + n) - 4 * n, n, oo) == -oo  # ex4 v1 leurre « a = 4 »
assert limit(sqrt(n**2 + 6 * n) - n, n, oo) not in (6, 0)  # ex4 v2 leurres « vers 6 », « vers 0 »
assert limit(sqrt(n**2 + 6 * n) + n, n, oo) == oo  # ex4 v2 leurre « a = -1 »

# ---------------------------------------------------------------------------
# ar-t-suites-comparaison (carte)
# ---------------------------------------------------------------------------
assert all(2 * k + (-1) ** k >= 2 * k - 1 for k in range(100))  # exemple
assert all(F(k * k, k + 1) >= F(k, 2) for k in range(1, 200))  # exemple n^2/(n+1) >= n/2
assert limit(n**2 - 10 * n, n, oo) == oo  # affirmation 2
# Contre-exemples : v_n = n (pair) / 0 (impair) : u_n = n <= v_n pour n pair ; v_n = 0 pour n impair.
v = [k if k % 2 == 0 else 0 for k in range(100)]
assert all(k <= v[k] for k in range(0, 100, 2)) and all(v[k] == 0 for k in range(1, 100, 2)) and max(v) == 98
# Trou 1, leurre |u_n| <= |v_n| : u_n = n, v_n = -n (|u| <= |v|) et v_n -> -oo.
assert all(abs(k) <= abs(-k) for k in range(50)) and limit(-n, n, oo) == -oo

# ---------------------------------------------------------------------------
# ar-t-suites-comparaison-ex1..ex3 (+ v1, v2)
# ---------------------------------------------------------------------------
assert limit(n**2 - 5, n, oo) == oo and all(k * k - 5 * math.sin(k) >= k * k - 5 for k in range(200))  # ex1
assert limit(-2 * n + 3, n, oo) == -oo and all(-2 * k + 3 * math.cos(k * k) <= -2 * k + 3 for k in range(200))  # v1
assert limit(sqrt(n) - 1, n, oo) == oo  # ex1 v2
# ex2 : S_n >= sqrt(n) ; v1 : S_n >= n/2 ; v2 : T_n >= sqrt(n/2).
for N in range(1, 300):
    assert sum(1 / math.sqrt(k) for k in range(1, N + 1)) >= math.sqrt(N) - 1e-12
    assert sum(F(k, k + 1) for k in range(1, min(N, 60) + 1)) >= F(min(N, 60), 2)
    assert sum(1 / math.sqrt(k) for k in range(N, 2 * N + 1)) >= math.sqrt(N / 2) - 1e-12
# ex3 (+ v1, v2) : u_n = n!/b^n ; n0 = plus petit n tel que (n+1)/b >= 2 ; leurres faux ; u_{n0}.
for b, n0, leurres, u_n0 in ((2, 3, (1, 2, 4), F(3, 4)), (3, 5, (2, 6, 3), F(40, 81)), (4, 7, (3, 8, 4), F(315, 1024))):
    un = lambda k: F(math.factorial(k), b**k)
    ok = lambda k: un(k + 1) >= 2 * un(k)
    assert all(un(k + 1) / un(k) == F(k + 1, b) for k in range(1, 40))
    assert ok(n0) and not ok(n0 - 1) and all(ok(k) for k in range(n0, 200)) and n0 == 2 * b - 1
    for le in leurres:
        assert le != n0 and (not all(ok(k) for k in range(le, 200)) or le > n0)
    assert un(n0) == u_n0 and all(un(k) >= u_n0 * 2 ** (k - n0) for k in range(n0, 200))

# ---------------------------------------------------------------------------
# ar-t-suites-comparaison-ex4 : u_{n+1} = u_n + 1/u_n, u_0 = 1 : u_n^2 >= 2n + 1.
# ---------------------------------------------------------------------------
u = suite(1, lambda t: t + 1 / t, 12)  # fractions exactes : la taille double à chaque pas
assert all(u[k] >= 1 and u[k] ** 2 >= 2 * k + 1 for k in range(13))
assert all(u[k + 1] ** 2 == u[k] ** 2 + 2 + 1 / u[k] ** 2 for k in range(12))

print("ar_t_suites : toutes les vérifications passent.")
