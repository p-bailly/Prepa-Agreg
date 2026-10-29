# Vérification calculatoire de l'enrichissement du chapitre « geo-l1-complexes » (thème geometrie, niveau l1) :
# questions de content/qcm/geo-l1-complexes.json et exercices ajoutés à
# content/exercices/geometrie/l1-geo-l1-complexes.json (similitudes-ex5, euler-ex6, racines-carrees-ex6).
# Calcul exact avec SymPy. Chaque assert cite l'identifiant de la question ou de l'exercice. Voir CLAUDE.md,
# « Contrôle de l'exactitude ».
import random
import sympy as sp
from sympy import I, sqrt, pi, exp, cos, sin, binomial, simplify, expand, Abs, symbols

z, w = symbols("z w")
x, th = symbols("x theta", real=True)
a, b, c = symbols("a b c")


def nul(e):
    return simplify(sp.expand_complex(sp.nsimplify(e))) == 0


# geo-l1-complexes-q1 : racines carrées de -5+12i = ±(2+3i) ; carrés des leurres
assert expand((2 + 3 * I) ** 2) == -5 + 12 * I
for l in (3 + 2 * I, 2 - 3 * I, 3 - 2 * I):
    assert expand(l ** 2) != -5 + 12 * I
assert expand((3 + 2 * I) ** 2) == 5 + 12 * I and expand((2 - 3 * I) ** 2) == -5 - 12 * I
assert expand((3 - 2 * I) ** 2) == 5 - 12 * I

# geo-l1-complexes-q2 : z^2-(3+i)z+2+i = (z-1)(z-2-i) ; Δ = 2i = (1+i)^2 ; leurres non solutions
P = z ** 2 - (3 + I) * z + 2 + I
assert expand(P - (z - 1) * (z - 2 - I)) == 0
assert expand((3 + I) ** 2 - 4 * (2 + I)) == 2 * I and expand((1 + I) ** 2) == 2 * I
for paire in ((-1, -2 - I), (1 + I, 2), (I, 3)):
    assert not all(expand(P.subs(z, r)) == 0 for r in paire)
assert expand((1 + I) * 2) == 2 + 2 * I and expand(I * 3) == 3 * I

# geo-l1-complexes-q3 (faux) : racines carrées de 3-4i = ±(2-i), parties réelle et imaginaire de signes contraires
assert expand((2 - I) ** 2) == 3 - 4 * I
u_, v_ = symbols("u v", real=True)
sols = sp.solve([u_ ** 2 - v_ ** 2 - 3, 2 * u_ * v_ + 4], [u_, v_], dict=True)
sols_reelles = [s for s in sols if s[u_].is_real and s[v_].is_real]
assert len(sols_reelles) == 2 and all(s[u_] * s[v_] < 0 for s in sols_reelles)

# geo-l1-complexes-q4 (vrai) : identité du parallélogramme (contrôle sur des complexes au hasard, exacts)
random.seed(1)
for _ in range(20):
    zz = random.randint(-9, 9) + random.randint(-9, 9) * I
    ww = random.randint(-9, 9) + random.randint(-9, 9) * I
    assert simplify(Abs(zz + ww) ** 2 + Abs(zz - ww) ** 2 - 2 * (Abs(zz) ** 2 + Abs(ww) ** 2)) == 0
zr, zi, wr, wi = symbols("zr zi wr wi", real=True)
Z, W = zr + I * zi, wr + I * wi
assert expand(sp.expand_complex(Abs(Z + W) ** 2 + Abs(Z - W) ** 2 - 2 * (Abs(Z) ** 2 + Abs(W) ** 2))) == 0

# geo-l1-complexes-q5 : |z+w| = |z|+|w| pour w = 3+6i seulement
z0 = 1 + 2 * I
egal = lambda ww: simplify(Abs(z0 + ww) - Abs(z0) - Abs(ww)) == 0
assert egal(3 + 6 * I)
assert not egal(-2 - 4 * I) and not egal(2 + I) and not egal(-4 + 2 * I)
assert Abs(-1 - 2 * I) == sqrt(5) and Abs(3 + 3 * I) == 3 * sqrt(2) and Abs(-3 + 4 * I) == 5

# geo-l1-complexes-q6 : h∘r(z) = 2iz - 1, r∘h(z) = 2iz - i
r = lambda t: I * t
h = lambda t: 2 * (t - 1) + 1
assert expand(h(r(z)) - (2 * I * z - 1)) == 0
assert expand(r(h(z)) - (2 * I * z - I)) == 0
for l in (2 * I * z - I, 2 * I * z + 1, 2 * I * z - 2):
    assert expand(h(r(z)) - l) != 0

# geo-l1-complexes-q7 (vrai) : r_B∘r_A(z) = z + (b-a)(1-conj(ρ)), translation non triviale si a ≠ b
rho = exp(I * pi / 3)
rA = lambda t: rho * (t - a) + a
rB = lambda t: sp.conjugate(rho) * (t - b) + b
assert simplify(sp.expand_complex(rB(rA(z)) - (z + (b - a) * (1 - sp.conjugate(rho))))) == 0
assert simplify(1 - sp.conjugate(rho)) != 0

# geo-l1-complexes-q8 : sin θ cos^3 θ = (sin 4θ + 2 sin 2θ)/8 ; leurres faux (valeurs en π/4)
f = sin(th) * cos(th) ** 3
assert simplify(sp.expand_trig(f - (sin(4 * th) + 2 * sin(2 * th)) / 8)) == 0
for l in ((sin(4 * th) - 2 * sin(2 * th)) / 8, (sin(4 * th) + 2 * sin(2 * th)) / 4, (sin(4 * th) + 3 * sin(2 * th)) / 8):
    assert simplify((f - l).subs(th, pi / 4)) != 0
assert f.subs(th, pi / 4) == sp.Rational(1, 4)

# geo-l1-complexes-similitudes-ex5
j = exp(2 * I * pi / 3)
# Q2 : 1 + j + j^2 = 0, -j^2 = e^{iπ/3}
assert nul(1 + j + j ** 2) and nul(-j ** 2 - exp(I * pi / 3))
# Q2 : c - a = -j^2 (b - a)  <=>  a + j b + j^2 c = 0 (la seconde vaut j^2 fois la première, au signe près)
assert nul(expand(j ** 2 * ((c - a) + j ** 2 * (b - a)) - (a + j * b + j ** 2 * c)).subs(j ** 3, 1)) or \
    nul(sp.expand_complex(j ** 2 * ((c - a) + j ** 2 * (b - a)) - (a + j * b + j ** 2 * c)))
# Q3 : (a + jb + j^2 c)(a + j^2 b + j c) = a^2+b^2+c^2-ab-bc-ca
assert nul(sp.expand_complex(expand((a + j * b + j ** 2 * c) * (a + j ** 2 * b + j * c)) - (a ** 2 + b ** 2 + c ** 2 - a * b - b * c - c * a)))
# Q4 : |u| = 1 et |u-1| = 1  <=>  u = e^{±iπ/3} ; triangle (0, 2, 1+i√3) équilatéral et vérifie la relation
phi = symbols("phi", real=True)
assert simplify(Abs(exp(I * phi) - 1) ** 2 - 4 * sin(phi / 2) ** 2) == 0 or \
    simplify(sp.expand_complex((exp(I * phi) - 1) * sp.conjugate(exp(I * phi) - 1)) - 4 * sin(phi / 2) ** 2) == 0
assert set(sp.solveset(2 * sin(phi / 2) - 1, phi, sp.Interval(-pi, pi))) | set(sp.solveset(2 * sin(phi / 2) + 1, phi, sp.Interval(-pi, pi))) == {pi / 3, -pi / 3}
A, B, C = 0, 2, 1 + I * sqrt(3)
assert nul(A ** 2 + B ** 2 + C ** 2 - A * B - B * C - C * A) and nul(A + j * B + j ** 2 * C)
# contre-épreuve : un triangle non équilatéral ne vérifie pas la relation
A, B, C = 0, 2, I
assert not nul(A ** 2 + B ** 2 + C ** 2 - A * B - B * C - C * A)

# geo-l1-complexes-euler-ex6 : C_n + i S_n = (1+e^{ix})^n, C_n = 2^n cos^n(x/2) cos(nx/2) ; leurres faux
assert nul(sp.expand_complex(1 + exp(I * x) - 2 * cos(x / 2) * exp(I * x / 2)))
for n in range(1, 8):
    Cn = sum(binomial(n, k) * cos(k * x) for k in range(n + 1))
    Sn = sum(binomial(n, k) * sin(k * x) for k in range(n + 1))
    bonne = 2 ** n * cos(x / 2) ** n * cos(n * x / 2)
    for xv in (sp.Rational(1, 3), sp.Rational(7, 5), 2):
        assert abs(sp.N((Cn - bonne).subs(x, xv), 30)) < 1e-25
        assert abs(sp.N((Sn - 2 ** n * cos(x / 2) ** n * sin(n * x / 2)).subs(x, xv), 30)) < 1e-25
    leurres = (2 ** n * cos(x / 2) ** n * sin(n * x / 2), 2 ** n * cos(x) ** n * cos(n * x), 2 ** n * cos(n * x / 2))
    for l in leurres:
        assert abs(sp.N((Cn - l).subs(x, sp.Rational(1, 3)), 30)) > 1e-6
    # Q3 : somme pour x = 2π/3 égale à cos(nπ/3)
    assert simplify(Cn.subs(x, 2 * pi / 3) - cos(n * pi / 3)) == 0

# geo-l1-complexes-racines-carrees-ex6 (et variantes v1, v2) : solutions de z^4 = w
cas = [
    (-119 - 120 * I, 5 - 12 * I, 3 - 2 * I, 2 + 3 * I, [(3 + 2 * I, 2 - 3 * I), (5 - 12 * I, 12 + 5 * I), (3 - 2 * I, 2 - 3 * I)]),
    (-7 - 24 * I, 3 - 4 * I, 2 - I, 1 + 2 * I, [(2 + I, 1 - 2 * I), (3 - 4 * I, 4 + 3 * I), (2 - I, 1 - 2 * I)]),
    (161 + 240 * I, 15 + 8 * I, 4 + I, 1 - 4 * I, [(4 - I, 1 + 4 * I), (15 + 8 * I, 8 - 15 * I), (4 + I, 1 + 4 * I)]),
]
for ww, d, r1, r2, leurres in cas:
    assert expand(d ** 2) == ww                      # racine carrée de w
    assert expand(r1 ** 2) == d and expand(r2 ** 2) == -d   # z^2 = δ et z^2 = -δ
    # les quatre nombres ±r1, ±r2 sont distincts et z^4 - w = (z-r1)(z+r1)(z-r2)(z+r2) : ce sont toutes les solutions
    assert len({r1, -r1, r2, -r2}) == 4
    assert expand((z - r1) * (z + r1) * (z - r2) * (z + r2) - (z ** 4 - ww)) == 0
    assert expand(I * r1) in (r2, -r2)               # Q3 : on passe de l'une à l'autre par multiplication par i
    for p, q in leurres:                             # chaque proposition fausse contient un non-solution
        assert expand(p ** 4) != ww or expand(q ** 4) != ww
    assert Abs(ww) == {(-119 - 120 * I): 169, (-7 - 24 * I): 25, (161 + 240 * I): 289}[ww]
    assert Abs(d) ** 2 == Abs(ww)

print("enrichi_geo_l1_complexes : tout est vérifié.")
