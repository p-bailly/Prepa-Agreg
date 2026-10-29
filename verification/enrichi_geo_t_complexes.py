"""Vérifications calculatoires de l'enrichissement du chapitre « geo-t-complexes » (Terminale, maths expertes).

Questions de content/qcm/geo-t-complexes.json et exercices ajoutés (exponentielle-ex5, module-ex5, second-degre-ex5)
de content/exercices/geometrie/terminale-geo-t-complexes.json.
Chaque assert cite l'identifiant de la question ou de l'exercice (v1, v2 = variantes).
Lancer : python3 verification/enrichi_geo_t_complexes.py (doit finir sans erreur).
"""

from sympy import (I, Abs, Rational, arg, conjugate, cos, exp, expand, expand_complex, im, log, nsimplify, pi, re,
                   simplify, sin, solve, sqrt, symbols, N, ceiling, floor)

x, y = symbols("x y", real=True)
z = symbols("z")
th = symbols("theta", real=True)


def nul(e):
    return simplify(expand_complex(e)) == 0


# --- QCM ---------------------------------------------------------------------
# geo-t-complexes-q1 : z = 2 − i, z² + conj(z) = 5 − 3i ; leurres faux.
w = 2 - I
assert expand(w**2 + conjugate(w)) == 5 - 3 * I
assert expand(w**2 + w) == 5 - 5 * I
assert all(v != 5 - 3 * I for v in (5 - 5 * I, 7 - 3 * I, 3 - 3 * I))

# geo-t-complexes-q2 (vrai) : z² + conj(z)² = 2 Re(z²) est réel.
zz = x + I * y
assert im(expand(zz**2 + conjugate(zz) ** 2)) == 0

# geo-t-complexes-q3 : |z − 1 + 2i| = 3 : cercle de centre 1 − 2i, rayon 3.
assert expand(Abs(zz - 1 + 2 * I) ** 2 - ((x - 1) ** 2 + (y + 2) ** 2)) == 0
for (centre, rayon) in ((-1 + 2 * I, 3), (1 - 2 * I, 9)):
    pt = centre + rayon  # point du leurre
    assert Abs(pt - 1 + 2 * I) != 3
# la « droite x − 2y = 3 » : le point (3, 0) y est, mais |3 − 1 + 2i| = √8 ≠ 3
assert Abs(3 - 1 + 2 * I) == 2 * sqrt(2)

# geo-t-complexes-q4 : arg (1 − i)^5 ≡ −5π/4 ≡ 3π/4 [2π] ; leurres faux.
assert nul((1 - I) ** 5 - 4 * sqrt(2) * exp(3 * I * pi / 4))
assert arg(expand((1 - I) ** 5)) == 3 * pi / 4
for a in (-pi / 4, pi / 4, -3 * pi / 4):
    assert simplify((a - 3 * pi / 4) / (2 * pi)).is_integer is False

# geo-t-complexes-q5 (faux) : e^{iθ} + e^{−iθ} = 2 cos θ ; en θ = 0 : 2 ≠ 0 = 2i sin 0.
assert nul(exp(I * th) + exp(-I * th) - 2 * cos(th))
assert (exp(0) + exp(0)) != 2 * I * sin(0)

# geo-t-complexes-q6 : racines sixièmes de 1 de partie réelle > 0 : k = 0, 1, 5, soit 3.
assert sum(1 for k in range(6) if cos(2 * k * pi / 6) > 0) == 3

# geo-t-complexes-q7 : z² + 2z + 4 = 0 : z = −1 ± i√3, de module 2.
sols = solve(z**2 + 2 * z + 4, z)
assert set(sols) == {-1 + I * sqrt(3), -1 - I * sqrt(3)} and all(Abs(r) == 2 for r in sols)

# geo-t-complexes-q8 (vrai) : (c − a)/(b − a) avec a = 1 + i, b = 3 + 2i, c = −3 − i vaut −2, réel.
a, b, c = 1 + I, 3 + 2 * I, -3 - I
assert simplify((c - a) / (b - a)) == -2

# --- geo-t-complexes-exponentielle-ex5 : suite géométrique de points ----------------------------------------------
def suite(q, module, theta, z0, seuil, n_min, periode, rapport, leurres):
    assert nul(q - module * exp(I * theta))
    assert simplify(Abs(q) - module) == 0
    # |z_n| = z0 · module^n < seuil : plus petit n
    absz = lambda k: z0 * module**k
    assert absz(n_min) < seuil and absz(n_min - 1) >= seuil
    assert all(not (absz(L) < seuil and absz(L - 1) >= seuil) for L in leurres)
    # z_n réel ⟺ n θ ≡ 0 [π] ⟺ n multiple de periode (vérifié pour n < 4·periode)
    for k in range(4 * periode):
        reel = simplify(im(expand_complex((module * exp(I * theta)) ** k))) == 0
        assert reel == (k % periode == 0)
    # (z_n − z_{n+1})/(−z_{n+1}) = (1 − q)/(−q) = 1 − 1/q
    assert nul((1 - q) / (-q) - rapport)
    assert simplify(re(expand_complex(rapport))) == 0  # imaginaire pur : angle droit en A_{n+1}


# main : q = (1+i)/2 = (√2/2)e^{iπ/4}, z0 = 4, seuil 10^{−2} ⇒ n = 18 ; réel ⟺ 4 | n ; rapport i.
suite((1 + I) / 2, sqrt(2) / 2, pi / 4, 4, Rational(1, 100), 18, 4, I, [17, 9, 14])
assert N(2 * log(400) / log(2)) > 17 and N(2 * log(400) / log(2)) < 18
# v1 : q = (1 + i√3)/4 = (1/2)e^{iπ/3}, z0 = 2, seuil 10^{−3} ⇒ n = 11 ; réel ⟺ 3 | n ; rapport i√3.
suite((1 + I * sqrt(3)) / 4, Rational(1, 2), pi / 3, 2, Rational(1, 1000), 11, 3, I * sqrt(3), [10, 12, 7])
# v2 : q = (3 + i√3)/4 = (√3/2)e^{iπ/6}, z0 = 1, seuil 10^{−1} ⇒ n = 17 ; réel ⟺ 6 | n ; rapport i√3/3.
suite((3 + I * sqrt(3)) / 4, sqrt(3) / 2, pi / 6, 1, Rational(1, 10), 17, 6, I * sqrt(3) / 3, [16, 8, 24])
assert Rational(3, 4) ** 8 == Rational(6561, 65536) and Rational(6561, 65536) > Rational(1, 10)

# --- geo-t-complexes-module-ex5 : |z − a| = k|z − b| ------------------------------------------------------------------
def apollonius(a, k, b, cx, cy, r, leurres_cercles):
    zz = x + I * y
    e = expand(Abs(zz - a) ** 2 - k**2 * Abs(zz - b) ** 2)
    cercle = (x - cx) ** 2 + (y - cy) ** 2 - r**2
    # e est un multiple constant (non nul) de l'équation du cercle
    q = simplify(e / cercle)
    assert q.is_number and q != 0
    # chaque leurre « cercle » contient un point qui ne vérifie pas l'égalité
    for (lx, ly, lr) in leurres_cercles:
        pt = lx + lr + I * ly
        assert simplify(Abs(pt - a) - k * Abs(pt - b)) != 0


# main : |z − 4| = 2|z − 1| ⟺ x² + y² = 4.
apollonius(4, 2, 1, 0, 0, 2, [(0, 0, 4), (2, 0, 2)])
assert Abs(Rational(5, 2) - 4) != 2 * Abs(Rational(5, 2) - 1)  # le milieu (médiatrice) ne convient pas
# v1 : |z + 3| = 2|z| ⟺ (x − 1)² + y² = 4.
apollonius(-3, 2, 0, 1, 0, 2, [(-1, 0, 2), (0, 0, sqrt(3))])
assert Abs(Rational(-3, 2) + 3) != 2 * Abs(Rational(-3, 2))
# v2 : |z − 9| = 3|z − 1| ⟺ x² + y² = 9.
apollonius(9, 3, 1, 0, 0, 3, [(0, 0, 9), (3, 0, 3)])
assert Abs(5 - 9) != 3 * Abs(5 - 1)

# --- geo-t-complexes-second-degre-ex5 ---------------------------------------------------------------------------------
# Δ = 4cos²θ − 4 = −4 sin²θ ; solutions e^{±iθ} ; puis z⁴ − 2z² cos θ + 1 = 0 : ±e^{±iθ/2}.
assert simplify((2 * cos(th)) ** 2 - 4 + 4 * sin(th) ** 2) == 0
for r in (exp(I * th), exp(-I * th)):
    assert nul(r**2 - 2 * r * cos(th) + 1)
P4 = lambda u: u**4 - 2 * u**2 * cos(th) + 1
for r in (exp(I * th / 2), -exp(I * th / 2), exp(-I * th / 2), -exp(-I * th / 2)):
    assert nul(P4(r))
# Leurre : e^{iθ} n'est pas solution de l'équation de degré 4 en général (θ = π/2 : 1 + 2·0... P4(i) = 2 ≠ 0).
assert P4(exp(I * pi / 2)).subs(th, pi / 2) == 2
assert simplify(P4(exp(2 * I * pi / 3)).subs(th, pi / 3)) != 0  # e^{2iθ} pour θ = π/3
# Les quatre solutions sont distinctes pour θ = π/3 (exemple numérique).
vals = [complex(N(v.subs(th, pi / 3))) for v in (exp(I * th / 2), -exp(I * th / 2), exp(-I * th / 2), -exp(-I * th / 2))]
assert all(abs(vals[i] - vals[j]) > 1e-9 for i in range(4) for j in range(i))

print("enrichi_geo_t_complexes : tout est vérifié.")
