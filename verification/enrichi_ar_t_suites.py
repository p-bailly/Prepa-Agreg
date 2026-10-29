"""Vérifications calculatoires des questions de QCM (content/qcm/ar-t-suites.json) et des exercices ajoutés au chapitre
ar-t-suites (content/exercices/analyse-reelle/terminale-ar-t-suites.json). Voir CLAUDE.md, « Contrôle de l'exactitude »."""
from sympy import (E, Rational as R, cos, factorial, limit, oo, simplify, sin, sqrt, symbols, summation, nsimplify, S)

n, k = symbols('n k', positive=True, integer=True)

# ar-t-suites-q1 : u_n = (4^n - 3^(n+1)) / 2^(2n+1) = 1/2 - (3/2)(3/4)^n, de limite 1/2 ; leurres 0, 2, +oo faux.
u = (4**n - 3**(n + 1)) / 2**(2 * n + 1)
for m in range(8):
    assert simplify(u.subs(n, m) - (R(1, 2) - R(3, 2) * R(3, 4)**m)) == 0
assert limit(u, n, oo) == R(1, 2)
assert all(R(1, 2) != x for x in (0, 2)) and limit(u, n, oo) != oo

# ar-t-suites-q2 : u_n = 2(-1,5)^n ; u_{2n} = 2·2,25^n -> +oo, u_{2n+1} = -3·2,25^n -> -oo : pas de limite (affirmation fausse).
q = R(-3, 2)
for m in range(6):
    assert 2 * q**(2 * m) == 2 * R(9, 4)**m and 2 * q**(2 * m + 1) == -3 * R(9, 4)**m
assert 2 * q**10 > 0 and 2 * q**11 < 0

# ar-t-suites-q3 : (sqrt(4n^2+1) - n)/(n+5) -> 1 ; leurres 3, 2, 0 faux.
u = (sqrt(4 * n**2 + 1) - n) / (n + 5)
assert limit(u, n, oo) == 1

# ar-t-suites-q4 : contre-exemples u_n = n, v_n = 1/sqrt(n) (produit -> +oo) et v_n = 3/n (produit = 3).
assert limit(1 / sqrt(n), n, oo) == 0 and limit(n * (1 / sqrt(n)), n, oo) == oo
assert simplify(n * (3 / n)) == 3

# ar-t-suites-q5 : n/2 - cos n >= n/2 - 1 -> +oo. Contre-exemple de monotonie : u_n = n + 1 + (-1)^n vérifie u_n >= n/2 - cos n.
for m in range(1, 200):
    um = m + 1 + (-1)**m
    assert um >= m / 2 - cos(m).evalf()
for j in range(1, 50):
    assert (2 * j + 1 + 1 + 1) > (2 * j + 1 + 1 - 1)  # u_{2j} = 2j+2 > u_{2j+1} = 2j+1

# ar-t-suites-q6 : u_0 = 1, u_{n+1} = u_n + 1/2^n ; u_n = 3 - 1/2^(n-1) pour n >= 1, < 3, limite 3 (affirmation vraie).
um = S(1)
for m in range(0, 30):
    if m >= 1:
        assert um == 3 - R(1, 2**(m - 1))
    assert um < 3
    um = um + R(1, 2**m)
assert limit(3 - 1 / 2**(n - 1), n, oo) == 3

# ar-t-suites-q7 : (2n ± sqrt(n))/(n+3) -> 2 (gendarmes) ; leurres 0, 2/3 faux.
assert limit((2 * n - sqrt(n)) / (n + 3), n, oo) == 2 and limit((2 * n + sqrt(n)) / (n + 3), n, oo) == 2

# ar-t-suites-q8 : u_0 = 1, u_{n+1} = u_n + 2n + 3 donne u_n = (n+1)^2 ; leurres 3n+1, 4^n, n^2+3n+1 faux.
termes = [1]
for m in range(10):
    termes.append(termes[-1] + 2 * m + 3)
assert all(termes[m] == (m + 1)**2 for m in range(11))
assert simplify((n + 1)**2 + 2 * n + 3 - (n + 2)**2) == 0
assert 3 * 2 + 1 == 7 != termes[2] and 4**2 == 16 != termes[2] and 1 + 3 + 1 == 5 != termes[1]

# ar-t-limite-monotone-ex4 : u_n = somme des 1/k!, v_n = u_n + 1/(n n!).
v_diff = 1 / factorial(n + 1) + 1 / ((n + 1) * factorial(n + 1)) - 1 / (n * factorial(n))
assert simplify(v_diff + 1 / (n * (n + 1) * factorial(n + 1))) == 0
assert simplify(n * (n + 1) + n - (n + 1)**2 + 1) == 0
uq = lambda m: sum(R(1, factorial(j)) for j in range(m + 1))
vq = lambda m: uq(m) + R(1, m * factorial(m))
assert uq(1) == 2 and vq(1) == 3
assert uq(4) == R(65, 24) and vq(4) == R(87, 32) == R(261, 96)
assert uq(3) == R(8, 3) and vq(3) == R(49, 18)
# leurres faux : (8/3, 49/18) sont u_3, v_3 ; v_4 != 11/4 ; v_4 != 65/24.
assert vq(4) != R(11, 4) and vq(4) != R(65, 24) and uq(4) != R(8, 3)
assert R(2708, 1000) < R(65, 24) < E < R(87, 32) == R(271875, 100000)
for m in range(1, 12):
    assert uq(m) < E.evalf(50) < vq(m)
    assert factorial(m) * (vq(m) - uq(m)) == R(1, m)
    assert (factorial(m) * uq(m)).is_integer

# ar-t-gendarmes-ex4 : S_n = somme n/(n^2+k), k = 1..n ; bornes n^2/(n^2+n), n^2/(n^2+1) -> 1.
assert limit(n**2 / (n**2 + n), n, oo) == 1 and limit(n**2 / (n**2 + 1), n, oo) == 1
for m in (1, 2, 5, 20):
    s = sum(R(m, m**2 + j) for j in range(1, m + 1))
    assert R(m**2, m**2 + m) <= s <= R(m**2, m**2 + 1)
# variante 1 : somme sur k = 1..2n ; bornes 2n^2/(n^2+2n), 2n^2/(n^2+1) -> 2.
assert limit(2 * n**2 / (n**2 + 2 * n), n, oo) == 2 and limit(2 * n**2 / (n**2 + 1), n, oo) == 2
for m in (1, 2, 5, 20):
    s = sum(R(m, m**2 + j) for j in range(1, 2 * m + 1))
    assert R(2 * m**2, m**2 + 2 * m) <= s <= R(2 * m**2, m**2 + 1)
# variante 2 : somme des 1/sqrt(n^2+k) ; bornes n/sqrt(n^2+n), n/sqrt(n^2+1) -> 1.
assert limit(n / sqrt(n**2 + n), n, oo) == 1 and limit(n / sqrt(n**2 + 1), n, oo) == 1
for m in (1, 2, 5, 20):
    s = sum(1 / sqrt(m**2 + j) for j in range(1, m + 1)).evalf()
    assert (m / sqrt(m**2 + m)).evalf() <= s <= (m / sqrt(m**2 + 1)).evalf()
assert abs(sum(R(400, 400**2 + j) for j in range(1, 401)) - 1) < R(1, 100)

# ar-t-suites-comparaison-ex5 : u_{n+1} = u_n + 3 + 2cos n >= u_n + 1, u_n >= n + 1 ; variantes : >= 2 + n/2 et <= -n/2.
a, b, c = 1.0, 2.0, 0.0
for m in range(300):
    assert a >= m + 1 - 1e-9 and b >= 2 + m / 2 - 1e-9 and c <= -m / 2 + 1e-9
    assert 3 + 2 * float(cos(m)) >= 1 and 1 - float(cos(m)) / 2 >= 0.5 and -1 + float(sin(m)) / 2 <= -0.5
    a += 3 + 2 * float(cos(m)); b += 1 - float(cos(m)) / 2; c += -1 + float(sin(m)) / 2
# Relecture indépendante (2026-10-29) : contrôles complémentaires.
from sympy import Rational as _R, factorial as _fact, limit as _lim, oo as _oo, symbols as _sym
_n = _sym('n', positive=True, integer=True)
# ar-t-limite-monotone-ex4 : v_n - u_n = 1/(n n!) -> 0, donc même limite.
assert _lim(1 / (_n * _fact(_n)), _n, _oo) == 0
# ar-t-suites-q6 : itération directe, u_n = 3 - 1/2^(n-1) pour n >= 1.
_u = _R(1)
for _m in range(12):
    _u += _R(1, 2**_m)
    assert _u == 3 - _R(1, 2**_m)
# ar-t-suites-q8 : itération directe, u_n = (n+1)^2 jusqu'à n = 30.
_u = 1
for _m in range(30):
    assert _u == (_m + 1)**2
    _u += 2 * _m + 3
print('enrichi_ar_t_suites : tout est vérifié')
