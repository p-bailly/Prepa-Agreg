"""Vérifications calculatoires du chapitre « arith-l1-polynomes » (thème arithmétique, niveau l1).

Cartes content/cartes/arithmetique/l1-arith-l1-polynomes.json et exercices
content/exercices/arithmetique/l1-arith-l1-polynomes.json (la carte arith-l1-racines-polynome, dans l1.json, est
vérifiée ailleurs). Chaque assert cite l'identifiant de la carte ou de l'exercice vérifié (variante comprise :
« v1 », « v2 »). Les leurres des questions à choix sont vérifiés FAUX ; la bonne réponse du JSON est relue.
Lancer : python3 verification/arith_l1_polynomes.py (doit finir sans erreur).
"""

import json
import pathlib
from itertools import permutations

import sympy
from sympy import (I, Poly, Rational, binomial, cos, diff, div, exp, expand, factor, factorial, gcd, nsimplify, pi, rem,
                   simplify, sqrt, symbols)

RACINE = pathlib.Path(__file__).resolve().parent.parent
X, n = symbols("X n")

CARTES = {c["id"]: c for c in json.loads((RACINE / "content/cartes/arithmetique/l1-arith-l1-polynomes.json").read_text())["cartes"]}
EXOS = {e["id"]: e for e in json.loads((RACINE / "content/exercices/arithmetique/l1-arith-l1-polynomes.json").read_text())["exercices"]}


def reponse(eid, v=0):
    """Texte de la bonne proposition de l'exercice eid (v = 0 : énoncé, v = 1, 2 : variantes)."""
    e = EXOS[eid]
    r = e["reponse_finale"] if v == 0 else e["variantes"][v - 1]["reponse_finale"]
    return r["propositions"][r["bonne"]]


def propositions(eid, v=0):
    e = EXOS[eid]
    r = e["reponse_finale"] if v == 0 else e["variantes"][v - 1]["reponse_finale"]
    return r["propositions"]


def egal(a, b):
    d = expand(a - b)
    if d.has(I) or d.has(exp) or d.has(cos):
        d = sympy.expand_complex(d)
    return simplify(d) == 0


def multiplicite(P, a):
    k = 0
    while egal(diff(P, X, k).subs(X, a), 0):
        k += 1
    return k


# Exercices ajoutés le 2026-10-29 (enrichissement) : vérifiés dans enrichi_arith_l1_polynomes.py, écartés ici.
AJOUTS = {"arith-l1-polynomes-division-ex5", "arith-l1-polynomes-viete-ex5", "arith-l1-polynomes-bezout-ex5"}
assert AJOUTS <= set(EXOS)
EXOS = {k: v for k, v in EXOS.items() if k not in AJOUTS}
# Structure : 6 cartes, 4 exercices par carte, 2 variantes par exercice sauf la démonstration pure.
assert len(CARTES) == 6 and len(EXOS) == 24
for eid, e in EXOS.items():
    nb = len(e.get("variantes", []))
    assert nb == (0 if eid == "arith-l1-polynomes-multiplicite-ex3" else 2), eid

# ===========================================================================
# arith-l1-polynomes-division (carte)
# ===========================================================================
# exemple : X^3 + 2X + 1 = (X^2 + 1) X + (X + 1) ; reste de X^5 + 1 par X - 2 = 33
assert div(X**3 + 2 * X + 1, X**2 + 1) == (X, X + 1)
assert rem(X**5 + 1, X - 2) == 33
# contre-exemple / affirmation 1 : X par 2X impossible dans Z[X] (X - c = 2X Q => 1 = 2 q0), possible dans Q[X]
assert div(X, 2 * X) == (Rational(1, 2), 0)
# contre-exemple B = 0 (corrigé après relecture) : deg R < deg 0 = -oo n'est vérifié par aucun R, même R = 0
assert sympy.degree(0, X) == -sympy.oo and not (sympy.degree(0, X) < sympy.degree(0, X))
# contre-exemple : deg R <= deg B casse l'unicité : X = X*1 + 0 = X*0 + X
assert egal(X * 1 + 0, X) and egal(X * 0 + X, X)
# affirmation 3 (fausse) : X^2 = (X+1)(X-2) + (X+2) mais le vrai reste est 1
assert egal((X + 1) * (X - 2) + (X + 2), X**2) and rem(X**2, X + 1) == 1
# trou 2, leurre « deg R < deg A » : A = 1, B = X, pas de R de degré < 0 non nul -> 1 = X Q impossible
assert rem(1, X) == 1


def verifie_division(eid, v, A, B, Q, R, leurres):
    assert div(A, B) == (expand(Q), expand(R)), (eid, v)
    for (q, r) in leurres:  # chaque leurre : BQ + R != A ou deg R >= deg B
        ok = egal(B * q + r, A) and Poly(r, X).degree() < Poly(B, X).degree()
        assert not ok, (eid, v, q, r)


# arith-l1-polynomes-division-ex1 (et v1, v2)
verifie_division("ex1", 0, X**4 + 2 * X**3 - X + 3, X**2 + X + 1, X**2 + X - 2, 5,
                 [(X**2 + X - 2, -2 * X + 5), (X**2 + X, -2 * X**2 - 2 * X + 3), (X**2 + 3 * X + 2, -6 * X + 1)])
assert reponse("arith-l1-polynomes-division-ex1") == "$Q = X^2 + X - 2$ et $R = 5$"
verifie_division("ex1", 1, X**5 - 3 * X**3 + 2 * X - 1, X**2 - 2 * X + 2, X**3 + 2 * X**2 - X - 6, -8 * X + 11,
                 [(X**3 + 2 * X**2 - X - 6, 8 * X - 11), (X**3 - 2 * X**2 - X + 6, -8 * X + 11),
                  (X**3 + 2 * X**2 - X, -6 * X**2 + 4 * X - 1)])
assert reponse("arith-l1-polynomes-division-ex1", 1) == "$Q = X^3 + 2X^2 - X - 6$ et $R = -8X + 11$"
verifie_division("ex1", 2, 2 * X**4 - X**3 + X**2 + 4, X**2 + 1, 2 * X**2 - X - 1, X + 5,
                 [(2 * X**2 - X - 1, -X + 5), (2 * X**2 - X + 3, X + 1), (2 * X**2 - X, -X**2 + X + 4)])
assert reponse("arith-l1-polynomes-division-ex1", 2) == "$Q = 2X^2 - X - 1$ et $R = X + 5$"
# étapes intermédiaires de la correction ex1
assert expand(X**4 + 2 * X**3 - X + 3 - X**2 * (X**2 + X + 1)) == X**3 - X**2 - X + 3
assert expand(X**3 - X**2 - X + 3 - X * (X**2 + X + 1)) == -2 * X**2 - 2 * X + 3
# étapes v1
B1 = X**2 - 2 * X + 2
assert expand(X**5 - 3 * X**3 + 2 * X - 1 - X**3 * B1) == 2 * X**4 - 5 * X**3 + 2 * X - 1
assert expand(2 * X**4 - 5 * X**3 + 2 * X - 1 - 2 * X**2 * B1) == -X**3 - 4 * X**2 + 2 * X - 1
assert expand(-X**3 - 4 * X**2 + 2 * X - 1 + X * B1) == -6 * X**2 + 4 * X - 1
assert expand(-6 * X**2 + 4 * X - 1 + 6 * B1) == -8 * X + 11


def verifie_reste(eid, v, A, B, R, leurres):
    assert egal(rem(A, B), R), (eid, v)
    for r in leurres:
        assert not egal(rem(A, B), r), (eid, v, r)


# arith-l1-polynomes-division-ex2 (et v1, v2) : leurres faux
verifie_reste("division-ex2", 0, X**10, (X - 1) * (X - 2), 1023 * X - 1022, [1024 * X - 1023, 1023 * X + 1, -1023 * X + 1024])
verifie_reste("division-ex2", 1, X**9, (X + 1) * (X - 3), 4921 * X + 4920, [4920 * X + 4921, 4921 * X - 4920, -4921 * X - 4922])
verifie_reste("division-ex2", 2, X**7, (X - 2) * (X + 1), 43 * X + 42, [42 * X + 43, 43 * X - 42, -43 * X + 214])
assert 3**9 == 19683 and 2**10 == 1024 and 19684 // 4 == 4921
assert reponse("arith-l1-polynomes-division-ex2") == "$1023X - 1022$"
assert reponse("arith-l1-polynomes-division-ex2", 1) == "$4921X + 4920$"
assert reponse("arith-l1-polynomes-division-ex2", 2) == "$43X + 42$"

# arith-l1-polynomes-division-ex3 (et v1, v2) : réduction modulo X^2 + X + 1
B3 = X**2 + X + 1
assert 2027 == 3 * 675 + 2 and 1000 == 3 * 333 + 1 and 500 == 3 * 166 + 2 and 1001 == 3 * 333 + 2 and 100 == 3 * 33 + 1
verifie_reste("division-ex3", 0, X**2027 + X**1000, B3, -1, [X**2 + X, 1, -2 * X - 1])
verifie_reste("division-ex3", 1, X**500 + 2 * X**300, B3, -X + 1, [X**2 + 2, -X - 1, X + 1])
verifie_reste("division-ex3", 2, X**1001 - X**100, B3, -2 * X - 1, [X**2 - X, -1, 2 * X + 1])
assert reponse("arith-l1-polynomes-division-ex3") == "$-1$"
assert reponse("arith-l1-polynomes-division-ex3", 1) == "$-X + 1$"
assert reponse("arith-l1-polynomes-division-ex3", 2) == "$-2X - 1$"

# arith-l1-polynomes-division-ex4 (et v1, v2) : reste de X^a - 1 par X^b - 1 = X^(a mod b) - 1
for a in range(0, 40):
    for b in range(1, 12):
        assert egal(rem(X**a - 1, X**b - 1), X**(a % b) - 1), (a, b)
verifie_reste("division-ex4", 0, X**100 - 1, X**12 - 1, X**4 - 1, [X**8 - 1, X**4, 0])
verifie_reste("division-ex4", 1, X**75 - 1, X**9 - 1, X**3 - 1, [X**8 - 1, X**3, X**6 - 1])
verifie_reste("division-ex4", 2, X**58 - 1, X**7 - 1, X**2 - 1, [X**8 - 1, X**2, X**3 - 1])
assert reponse("arith-l1-polynomes-division-ex4") == "$X^4 - 1$"

# ===========================================================================
# arith-l1-polynomes-multiplicite (carte)
# ===========================================================================
# exemple : X^4 - 2X^2 + 1, racine 1 double ; P''(1) = 8
P = X**4 - 2 * X**2 + 1
assert multiplicite(P, 1) == 2 and diff(P, X, 2).subs(X, 1) == 8 and egal(P, (X - 1)**2 * (X + 1)**2)
# contre-exemple : dans (Z/2Z)[X], X^2 a toutes ses dérivées nulles (coefficients 2 et 2 = 0 mod 2) ;
# dans (Z/pZ)[X], dérivées de X^p : p X^(p-1), p(p-1) X^(p-2)... toutes à coefficients multiples de p
for p in [2, 3, 5, 7]:
    for k in range(1, p + 1):
        assert (factorial(p) // factorial(p - k)) % p == 0, (p, k)
# contre-exemple / affirmation 1 : X^2 + 1, P'(0) = 0 mais P(0) = 1
assert diff(X**2 + 1, X).subs(X, 0) == 0 and (X**2 + 1).subs(X, 0) == 1
# trous, leurre « P^(m-1)(α) = 0 » : X^2 - 1, α = 0, m = 2 : P'(0) = 0 mais X^2 ne divise pas X^2 - 1
assert diff(X**2 - 1, X).subs(X, 0) == 0 and rem(X**2 - 1, X**2) != 0
# leurre « P(α) = 0 et P^(m)(α) ≠ 0 » : X + X^2, α = 0, m = 2 : P(0) = 0, P''(0) = 2, X^2 ne divise pas
assert (X + X**2).subs(X, 0) == 0 and diff(X + X**2, X, 2) == 2 and rem(X + X**2, X**2) != 0
# leurre « jusqu'à P^(m) » : P = X, m = 1 : X | X mais P'(0) = 1
assert diff(X, X) == 1
# affirmation 2 (vraie) : multiplicité m dans P => m - 1 dans P' ; test sur des exemples
for m in range(1, 6):
    Pm = (X - 2)**m * (X**2 + 1)
    assert multiplicite(Pm, 2) == m and multiplicite(diff(Pm, X), 2) == m - 1

# arith-l1-polynomes-multiplicite-ex1 (et v1, v2)
P0 = X**5 - 3 * X**4 + 4 * X**3 - 4 * X**2 + 3 * X - 1
assert multiplicite(P0, 1) == 3 and diff(P0, X, 3).subs(X, 1) == 12 and egal(P0, (X - 1)**3 * (X**2 + 1))
assert egal(diff(P0, X), 5 * X**4 - 12 * X**3 + 12 * X**2 - 8 * X + 3)
assert egal(diff(P0, X, 2), 20 * X**3 - 36 * X**2 + 24 * X - 8) and egal(diff(P0, X, 3), 60 * X**2 - 72 * X + 24)
P1 = X**5 - 5 * X + 4
assert multiplicite(P1, 1) == 2 and diff(P1, X, 2).subs(X, 1) == 20 and egal(P1, (X - 1)**2 * (X**3 + 2 * X**2 + 3 * X + 4))
P2 = X**5 + 7 * X**4 + 16 * X**3 + 8 * X**2 - 16 * X - 16
assert multiplicite(P2, -2) == 4 and diff(P2, X, 4).subs(X, -2) == -72 and egal(P2, (X + 2)**4 * (X - 1))
assert egal(diff(P2, X), 5 * X**4 + 28 * X**3 + 48 * X**2 + 16 * X - 16)
assert egal(diff(P2, X, 2), 20 * X**3 + 84 * X**2 + 96 * X + 16) and egal(diff(P2, X, 3), 60 * X**2 + 168 * X + 96)
assert reponse("arith-l1-polynomes-multiplicite-ex1") == "$3$"
assert reponse("arith-l1-polynomes-multiplicite-ex1", 1) == "$2$"
assert reponse("arith-l1-polynomes-multiplicite-ex1", 2) == "$4$"
# les leurres sont d'autres entiers : distincts de la multiplicité (vérifié par la distinction des propositions)
for v, m in [(0, 3), (1, 2), (2, 4)]:
    assert [p for p in propositions("arith-l1-polynomes-multiplicite-ex1", v) if p == f"${m}$"] == [f"${m}$"]


def racine_double(P, a):
    return egal(P.subs(X, a), 0) and egal(diff(P, X).subs(X, a), 0)


# arith-l1-polynomes-multiplicite-ex2 (et v1, v2) : (a, b) solution unique, leurres faux
a_, b_ = symbols("a b")
for v, Pab, r, bonne, leurres in [
    (0, a_ * X**5 + b_ * X**4 + 1, 1, (4, -5), [(-5, 4), (5, -6), (-4, 3)]),
    (1, X**4 + a_ * X + b_, -1, (4, 3), [(-4, -5), (4, -3), (-4, 3)]),
    (2, X**3 + a_ * X + b_, 2, (-12, 16), [(12, -32), (-12, -16), (-4, 0)]),
]:
    sol = sympy.solve([Pab.subs(X, r), diff(Pab, X).subs(X, r)], [a_, b_], dict=True)
    assert sol == [{a_: bonne[0], b_: bonne[1]}], (v, sol)
    for (aa, bb) in leurres:
        assert not racine_double(Pab.subs({a_: aa, b_: bb}), r), (v, aa, bb)
assert egal(4 * X**5 - 5 * X**4 + 1, (X - 1)**2 * (4 * X**3 + 3 * X**2 + 2 * X + 1))
assert egal(X**4 + 4 * X + 3, (X + 1)**2 * (X**2 - 2 * X + 3))
assert egal(X**3 - 12 * X + 16, (X - 2)**2 * (X + 4))
assert reponse("arith-l1-polynomes-multiplicite-ex2") == "$(4, -5)$"
assert reponse("arith-l1-polynomes-multiplicite-ex2", 1) == "$(4, 3)$"
assert reponse("arith-l1-polynomes-multiplicite-ex2", 2) == "$(-12, 16)$"

# arith-l1-polynomes-multiplicite-ex3 : P_n - P_n' = X^n / n!, racines simples (contrôle numérique n <= 12)
for k in range(1, 13):
    Pk = sum(X**j / factorial(j) for j in range(k + 1))
    assert egal(Pk - diff(Pk, X), X**k / factorial(k))
    assert gcd(Pk, diff(Pk, X)) == 1, k

# arith-l1-polynomes-multiplicite-ex4 (et v1, v2) : familles à paramètre n
for k in range(2, 12):
    Pk = X**(2 * k) - k * X**(k + 1) + k * X**(k - 1) - 1
    assert multiplicite(Pk, 1) == 3 and diff(Pk, X, 3).subs(X, 1) == 2 * k * (k - 1) * (k + 1)
    Qk = sympy.quo(Pk, (X - 1)**3)
    assert rem(Pk, (X - 1)**3) == 0 and Qk.subs(X, 1) == Rational(k * (k**2 - 1), 3)
    # leurres : 2n(n^2 - 1), n(n^2 - 1)/6, n(n+1)(n+2)/3 sont faux pour tout n >= 2
    for faux in [2 * k * (k**2 - 1), Rational(k * (k**2 - 1), 6), Rational(k * (k + 1) * (k + 2), 3)]:
        assert Qk.subs(X, 1) != faux, k
assert egal(X**4 - 2 * X**3 + 2 * X - 1, (X - 1)**3 * (X + 1))  # cas n = 2 cité
# calcul littéral de P_n'''(1) cité dans la correction
assert expand(2 * n * (2 * n - 1) * (2 * n - 2) - n**2 * (n + 1) * (n - 1) + n * (n - 1) * (n - 2) * (n - 3)
              - 2 * n * (n - 1) * (n + 1)) == 0
assert expand(4 * (2 * n - 1) - n * (n + 1) + (n - 2) * (n - 3) - (2 * n + 2)) == 0
for k in range(1, 12):
    Ak = k * X**(k + 1) - (k + 1) * X**k + 1  # v1
    assert multiplicite(Ak, 1) == 2 and diff(Ak, X, 2).subs(X, 1) == k * (k + 1)
    Bk = sympy.quo(Ak, (X - 1)**2)
    assert rem(Ak, (X - 1)**2) == 0 and Bk.subs(X, 1) == Rational(k * (k + 1), 2)
    for faux in [k * (k + 1), Rational(k * (k - 1), 2), Rational((k + 1) * (k + 2), 2)]:
        assert Bk.subs(X, 1) != faux, k
    Rk = X**(2 * k + 1) - (2 * k + 1) * X**(k + 1) + (2 * k + 1) * X**k - 1  # v2
    assert multiplicite(Rk, 1) == 3 and diff(Rk, X, 3).subs(X, 1) == (2 * k + 1) * k * (k + 1)
    Sk = sympy.quo(Rk, (X - 1)**3)
    assert rem(Rk, (X - 1)**3) == 0 and Sk.subs(X, 1) == Rational(k * (k + 1) * (2 * k + 1), 6)
    for faux in [k * (k + 1) * (2 * k + 1), Rational(k * (k + 1) * (2 * k + 1), 3), Rational(k * (k + 1) * (2 * k - 1), 6)]:
        assert Sk.subs(X, 1) != faux, k
assert expand(2 * (2 * n - 1) - (n + 1) * (n - 1) + (n - 1) * (n - 2) - (n + 1)) == 0

# ===========================================================================
# arith-l1-polynomes-irreductibles (carte)
# ===========================================================================
assert egal(X**3 + 1, (X + 1) * (X**2 - X + 1)) and (-1)**2 - 4 == -3  # exemple
assert egal((X - exp(I * pi / 3)) * (X - exp(-I * pi / 3)), X**2 - X + 1)
# contre-exemples : X^4 + 1 factorisé dans R[X], sans racine réelle ; (X^2+1)^2 réductible
assert egal((X**2 + sqrt(2) * X + 1) * (X**2 - sqrt(2) * X + 1), X**4 + 1)
assert sympy.real_roots(Poly(X**4 + 1, X)) == [] and sympy.real_roots(Poly((X**2 + 1)**2, X)) == []
# X^2 - 2 et X^3 - 2 irréductibles sur Q
assert Poly(X**2 - 2, X, domain="QQ").is_irreducible and Poly(X**3 - 2, X, domain="QQ").is_irreducible
# trous : discriminant nul (X^2) réductible ; « degré pair » : X^4 + 1 réductible ; degré 3 réductible sur R
assert sympy.real_roots(Poly(X**3 + X + 1, X)) != []
# affirmation 4 : z et conj(z) même multiplicité, exemple (X^2 + 1)^2 (X - 3)
Pz = expand((X**2 + 1)**2 * (X - 3))
assert multiplicite(Pz, I) == 2 and multiplicite(Pz, -I) == 2

# arith-l1-polynomes-irreductibles-ex1 (et v1, v2) : la bonne factorisation est égale au polynôme et ses facteurs
# sont irréductibles sur R ; chaque leurre est faux (produit différent, facteur réductible ou coefficients non réels)
def irr_reel(f):
    p = Poly(f, X)
    if not all(c.is_real for c in p.all_coeffs()):
        return False
    if p.degree() == 1:
        return True
    if p.degree() == 2:
        a2, a1, a0 = p.all_coeffs()
        return all(c.is_real for c in (a2, a1, a0)) and simplify(a1**2 - 4 * a2 * a0) < 0
    return False


def bonne_factorisation(P, facteurs):
    prod = 1
    for f in facteurs:
        prod *= f
    return egal(prod, P) and all(irr_reel(f) for f in facteurs)


s3 = sqrt(3)
for v, P, bonne, leurres in [
    (0, X**4 - 1, [X - 1, X + 1, X**2 + 1],
     [[X - 1, X + 1, X - I, X + I], [X**2 - 1, X**2 + 1], [X - 1, X - 1, X**2 + 1]]),
    (1, X**3 - 8, [X - 2, X**2 + 2 * X + 4],
     [[X - 2, X**2 - 2 * X + 4], [X - 2, X - 2, X - 2], [X - 2, X + 1 - I * s3, X + 1 + I * s3]]),
    (2, X**4 + 3 * X**2 - 4, [X - 1, X + 1, X**2 + 4],
     [[X**2 - 1, X**2 + 4], [X - 1, X + 1, X - 2 * I, X + 2 * I], [X - 2, X + 2, X**2 + 1]]),
]:
    assert bonne_factorisation(P, bonne), v
    for l in leurres:
        assert not bonne_factorisation(P, l), (v, l)
assert reponse("arith-l1-polynomes-irreductibles-ex1") == "$(X - 1)(X + 1)(X^2 + 1)$"
assert reponse("arith-l1-polynomes-irreductibles-ex1", 1) == "$(X - 2)(X^2 + 2X + 4)$"
assert reponse("arith-l1-polynomes-irreductibles-ex1", 2) == "$(X - 1)(X + 1)(X^2 + 4)$"

# arith-l1-polynomes-irreductibles-ex2 (et v1, v2) : bicarrés ; leurres « irréductible », carré, mauvais signes
for v, P, bonne, leurres in [
    (0, X**4 + 4, [X**2 - 2 * X + 2, X**2 + 2 * X + 2], [[X**4 + 4], [X**2 + 2, X**2 + 2], [X**2 - 2 * X - 2, X**2 + 2 * X - 2]]),
    (1, X**4 + 64, [X**2 - 4 * X + 8, X**2 + 4 * X + 8], [[X**4 + 64], [X**2 + 8, X**2 + 8], [X**2 - 4 * X - 8, X**2 + 4 * X - 8]]),
    (2, X**4 + X**2 + 1, [X**2 - X + 1, X**2 + X + 1], [[X**4 + X**2 + 1], [X**2 + 1, X**2 + 1], [X**2 - X - 1, X**2 + X - 1]]),
]:
    assert bonne_factorisation(P, bonne), v
    assert sympy.real_roots(Poly(P, X)) == []  # sans racine réelle, pourtant réductible
    for l in leurres:
        assert not bonne_factorisation(P, l), (v, l)

# arith-l1-polynomes-irreductibles-ex3 (et v1, v2) : nombre de facteurs irréductibles réels de X^n - 1
for v, k, nb, leurres in [(0, 6, 4, [2, 6, 3]), (1, 5, 3, [5, 2, 1]), (2, 8, 5, [8, 4, 3])]:
    f = sympy.factor_list(X**k - 1, extension=[sqrt(2), sqrt(5)])[1]  # facteurs sur Q(√2, √5) ⊂ R
    reels = [g for g, e in f]
    # chaque facteur trouvé est irréductible sur R (degré 1, ou degré 2 sans racine réelle)
    assert all(irr_reel(g) for g in reels), (v, reels)
    assert len(reels) == nb and nb not in leurres, v
    assert reponse("arith-l1-polynomes-irreductibles-ex3", v) == f"${nb}$"
assert egal(X**6 - 1, (X - 1) * (X + 1) * (X**2 - X + 1) * (X**2 + X + 1))
assert egal(2 * cos(2 * pi / 5), (sqrt(5) - 1) / 2) and egal(2 * cos(4 * pi / 5), -(1 + sqrt(5)) / 2)
assert egal(X**5 - 1, (X - 1) * (X**2 - (sqrt(5) - 1) / 2 * X + 1) * (X**2 + (1 + sqrt(5)) / 2 * X + 1))
assert egal(X**8 - 1, (X - 1) * (X + 1) * (X**2 + 1) * (X**2 - sqrt(2) * X + 1) * (X**2 + sqrt(2) * X + 1))

# arith-l1-polynomes-irreductibles-ex4 (et v1, v2) : l'énoncé (corrigé après relecture) impose les racines de Im > 0 ;
# d'autres choix de racines donneraient d'autres écritures justes, qui ne figurent pas parmi les propositions
assert egal((X**2)**2 + 1, X**4 + 1) and egal((X**2 + 2)**2 + X**2, X**4 + 5 * X**2 + 4) and egal((X**2)**2 + 2**2, X**4 + 4)
for v in range(3):
    assert all("(X^2)^2" not in p and "(X^2 + 2)^2 + X^2" not in p for p in propositions("arith-l1-polynomes-irreductibles-ex4", v))
# arith-l1-polynomes-irreductibles-ex4 (et v1, v2) : Q = produit des (X - w), w racines de Im > 0, puis A^2 + B^2
for v, P, w, A, B, leurres in [
    (0, X**4 + 1, [exp(I * pi / 4), exp(3 * I * pi / 4)], X**2 - 1, sqrt(2) * X,
     [(X**2 + 1, sqrt(2) * X), (X**2 - 1, 2 * X), (X**2 - 1, sqrt(2))]),
    (1, X**4 + 5 * X**2 + 4, [I, 2 * I], X**2 - 2, 3 * X,
     [(X**2 + 2, 3 * X), (X**2 - 2, 2 * X), (X**2 - 3, 2 * X)]),
    (2, X**4 + 4, [1 + I, -1 + I], X**2 - 2, 2 * X,
     [(X**2 + 2, 2 * X), (X**2 - 2, sqrt(2) * X), (X**2 - 1, 2 * X)]),
]:
    for z in w:
        assert egal(P.subs(X, z), 0) and sympy.im(z) > 0
    Q = expand((X - w[0]) * (X - w[1]))
    A1 = sum(sympy.re(c) * X**(3 - i - 1) for i, c in enumerate(Poly(Q, X).all_coeffs()))
    B1 = sum(sympy.im(c) * X**(3 - i - 1) for i, c in enumerate(Poly(Q, X).all_coeffs()))
    assert egal(A1, A) and egal(B1**2, B**2), v  # B au signe près
    assert egal(A**2 + B**2, P), v
    for (a, b) in leurres:
        assert not egal(a**2 + b**2, P), (v, a, b)

# ===========================================================================
# arith-l1-polynomes-viete (carte)
# ===========================================================================
assert egal(X**3 - 6 * X**2 + 11 * X - 6, (X - 1) * (X - 2) * (X - 3))  # exemple
assert egal(X**3 - 3 * X + 2, (X - 1)**2 * (X + 2)) and 1 + 1 - 2 == 0 and 1 - 2 == -1  # contre-exemple multiplicités
assert Rational(-4, 2) == -2  # affirmation 2 (fausse) : produit des racines de 2X^3 - X + 4 = -2, pas 2
assert sympy.real_roots(Poly(X**3 - 1, X)) == [1]  # affirmation 3 (fausse)
# trou 1 (corrigé après relecture) : leurres « non nul », « unitaire », « coefficients tous non nuls » réfutés par
# X^2 + 1 ∈ R[X], non nul, unitaire, à coefficients non nuls, sans aucune racine réelle (pas n = 2 racines dans R)
assert sympy.real_roots(Poly(X**2 + 1, X)) == []


def sigmas(P):
    c = Poly(P, X).all_coeffs()
    an = c[0]
    return [(-1)**k * c[k] / an for k in range(1, len(c))]


def racines(P):
    return [r for r, m in sympy.roots(Poly(P, X), multiple=False).items() for _ in range(m)]


# arith-l1-polynomes-viete-ex1 (et v1, v2) : somme des carrés
for v, P, val, leurres in [
    (0, X**3 - 2 * X**2 + 3 * X - 5, -2, [4, 10, -8]),
    (1, X**3 + X**2 - 4 * X + 1, 9, [1, -7, 3]),
    (2, 2 * X**3 - 6 * X**2 + X + 7, 8, [9, 10, 34]),
]:
    s1, s2, s3_ = sigmas(P)
    assert s1**2 - 2 * s2 == val, v
    # contrôle numérique sur les racines
    rs = sympy.Poly(P, X).nroots(n=30)
    assert abs(sum(r**2 for r in rs) - val) < 1e-20, v
    assert val not in leurres
    assert reponse("arith-l1-polynomes-viete-ex1", v) == f"${val}$"
# les racines de X^3 - 2X^2 + 3X - 5 ne sont pas toutes réelles
assert len(sympy.real_roots(Poly(X**3 - 2 * X**2 + 3 * X - 5, X))) == 1

# arith-l1-polynomes-viete-ex2 (et v1, v2) : somme des inverses = σ_{n-1}/σ_n
for v, P, val, leurres in [
    (0, X**3 - 3 * X**2 + 5 * X + 2, Rational(-5, 2), [Rational(5, 2), Rational(-2, 5), Rational(-3, 2)]),
    (1, X**3 + 2 * X**2 - X - 4, Rational(-1, 4), [Rational(1, 4), -4, Rational(-1, 2)]),
    (2, X**4 - X**3 + 2 * X - 3, Rational(2, 3), [Rational(-2, 3), Rational(3, 2), Rational(-1, 3)]),
]:
    assert P.subs(X, 0) != 0
    s = sigmas(P)
    assert s[-2] / s[-1] == val, v
    rs = sympy.Poly(P, X).nroots(n=30)
    assert abs(sum(1 / r for r in rs) - val) < 1e-20, v
    assert val not in leurres
assert reponse("arith-l1-polynomes-viete-ex2") == "$-\\frac{5}{2}$"
assert reponse("arith-l1-polynomes-viete-ex2", 1) == "$-\\frac{1}{4}$"
assert reponse("arith-l1-polynomes-viete-ex2", 2) == "$\\frac{2}{3}$"

# arith-l1-polynomes-viete-ex3 (et v1, v2) : λ tels qu'une racine soit le double d'une autre
lam = symbols("lam")


def double_d_une_autre(P):
    rs = racines(P)
    return any(i != j and egal(rs[i], 2 * rs[j]) for i in range(3) for j in range(3))


for v, c, sol, leurres in [
    (0, 7, {-6, 6}, [{6}, {-2, 2}, {Rational(-7, 3), Rational(7, 3)}]),
    (1, 28, {-48, 48}, [{48}, {-24, 24}, {-12, 12}]),
    (2, 63, {-162, 162}, [{162}, {-54, 54}, {-81, 81}]),
]:
    # condition nécessaire : a^2 = c/7, λ = 6a^3
    a = symbols("a")
    trouve = {6 * r**3 for r in sympy.solve(7 * a**2 - c, a)}
    assert trouve == sol, v
    for l in sol:
        assert double_d_une_autre(X**3 - c * X + l), (v, l)
    # chaque leurre est un ensemble faux : il contient une valeur qui ne convient pas, ou il en manque une
    for L in leurres:
        assert L != sol
        assert any(not double_d_une_autre(X**3 - c * X + l) for l in L) or not sol <= L, (v, L)
assert egal(X**3 - 7 * X + 6, (X - 1) * (X - 2) * (X + 3)) and egal(X**3 - 7 * X - 6, (X + 1) * (X + 2) * (X - 3))
assert egal(X**3 - 28 * X + 48, (X - 2) * (X - 4) * (X + 6)) and egal(X**3 - 28 * X - 48, (X + 2) * (X + 4) * (X - 6))
assert egal(X**3 - 63 * X + 162, (X - 3) * (X - 6) * (X + 9)) and egal(X**3 - 63 * X - 162, (X + 3) * (X + 6) * (X - 9))
assert reponse("arith-l1-polynomes-viete-ex3") == "$\\{-6, 6\\}$"

# arith-l1-polynomes-viete-ex4 (et v1, v2) : systèmes de Newton
for v, (p1, p2, p3), triplet, T, faux in [
    (0, (2, 6, 8), (-1, 1, 2), X**3 - 2 * X**2 - X + 2, [(1, -2, 3), (0, 1, 1)]),
    (1, (0, 26, -36), (1, 3, -4), X**3 - 13 * X + 12, [(-1, -3, 4)]),
    (2, (4, 14, 34), (-1, 2, 3), X**3 - 4 * X**2 + X + 6, [(1, -2, 5)]),
]:
    s1 = p1
    s2 = Rational(s1 * p1 - p2, 2)
    s3_ = Rational(p3 - s1 * p2 + s2 * p1, 3)
    assert egal(X**3 - s1 * X**2 + s2 * X - s3_, T), v
    assert sorted(racines(T)) == sorted(triplet), v
    # résolution directe du système : exactement les 6 permutations
    x, y, z = symbols("x y z")
    sols = sympy.solve([x + y + z - p1, x**2 + y**2 + z**2 - p2, x**3 + y**3 + z**3 - p3], [x, y, z], dict=False)
    assert set(map(tuple, sols)) == set(permutations(triplet)) and len(set(permutations(triplet))) == 6, v
    for f in faux:  # leurres : permutations d'un autre triplet, qui n'est pas solution
        assert (sum(f), sum(t**2 for t in f), sum(t**3 for t in f)) != (p1, p2, p3), (v, f)
# identités de Newton utilisées (question 1)
x, y, z = symbols("x y z")
S1, S2, S3 = x + y + z, x * y + x * z + y * z, x * y * z
P1, P2, P3 = x + y + z, x**2 + y**2 + z**2, x**3 + y**3 + z**3
assert expand(P2 - (S1 * P1 - 2 * S2)) == 0 and expand(P3 - (S1 * P2 - S2 * P1 + 3 * S3)) == 0

# ===========================================================================
# arith-l1-polynomes-bezout (carte)
# ===========================================================================
assert gcd(X**3 - 1, X**2 - 1) == X - 1 and egal((X**3 - 1) - X * (X**2 - 1), X - 1)  # exemple
assert gcd(X**2 - 1, 2 * X - 2) == X - 1  # contre-exemple normalisation
# trou 2, leurres : constantes U, V pour A = X^2, B = X + 1 ; U, V unitaires pour A = B = 1 ; même degré pour A = 1, B = X
u, w = symbols("u w")
assert sympy.solve(Poly(u * X**2 + w * (X + 1) - 1, X).all_coeffs(), [u, w]) == []
# affirmation 3 (fausse) : X*1 + 1*(-1) = X - 1 mais X ∧ 1 = 1
assert gcd(X, 1) == 1
# affirmation 2 (fausse) : X^2 + 1 et X^2 + 1
assert gcd(X**2 + 1, X**2 + 1) == X**2 + 1

# arith-l1-polynomes-bezout-ex1 (et v1, v2)
for v, A, B, D, leurres in [
    (0, X**4 - 1, X**3 - 1, X - 1, [1, X + 1, X**2 + X + 1]),
    (1, X**3 - X**2 - X + 1, X**3 + X**2 - X - 1, X**2 - 1, [X - 1, -2 * X**2 + 2, (X - 1)**2 * (X + 1)**2]),
    (2, X**3 + 2 * X**2 + 2 * X + 1, X**3 - 1, X**2 + X + 1, [2 * X**2 + 2 * X + 2, X + 1, 1]),
]:
    assert gcd(A, B) == D, v
    for l in leurres:
        assert not egal(gcd(A, B), l), (v, l)
assert rem(X**4 - 1, X**3 - 1) == X - 1
assert expand((X**3 - X**2 - X + 1) - (X**3 + X**2 - X - 1)) == -2 * X**2 + 2
assert expand((-2 * X**2 + 2) * (-X / 2 - Rational(1, 2)) - (X**3 + X**2 - X - 1)) == 0
assert expand((2 * X**2 + 2 * X + 2) * (X - 1) / 2 - (X**3 - 1)) == 0

# arith-l1-polynomes-bezout-ex2 (et v1, v2) : couple (U, V) unique avec deg U < deg B, deg V < deg A
for v, A, B, U, V, leurres in [
    (0, X**2 + 1, X - 1, Rational(1, 2), -(X + 1) / 2,
     [(Rational(1, 2), (X + 1) / 2), (1, -(X + 1)), (2, -2 * (X + 1))]),
    (1, X**2 + X + 1, X + 2, Rational(1, 3), -(X - 1) / 3,
     [(Rational(1, 3), (X - 1) / 3), (1, -(X - 1)), (Rational(1, 3), -(X + 1) / 3)]),
    (2, X**3 + 1, X**2 + 1, (X + 1) / 2, -(X**2 + X - 1) / 2,
     [((X + 1) / 2, (X**2 + X - 1) / 2), (X + 1, -(X**2 + X - 1)), ((X - 1) / 2, -(X**2 - X - 1) / 2)]),
]:
    su, sv, g = sympy.gcdex(A, B)
    assert g == 1 and egal(su, U) and egal(sv, V), v
    assert egal(A * U + B * V, 1), v
    for (uu, vv) in leurres:
        assert not egal(A * uu + B * vv, 1), (v, uu, vv)
assert egal((X + 2) * (X - 1) + 3, X**2 + X + 1) and egal((1 - X) * (-X - 1) + 2, X**2 + 1)
assert egal((X**3 + 1) * (X + 1) - (X**2 + 1) * (X**2 + X - 1), 2)

# arith-l1-polynomes-bezout-ex3 (et v1, v2) : (X^a - 1) ∧ (X^b - 1) = X^(a∧b) - 1
from math import gcd as pgcd_entier

for a in range(1, 25):
    for b in range(1, 25):
        assert gcd(X**a - 1, X**b - 1) == X**pgcd_entier(a, b) - 1, (a, b)
for v, a, b, bonne, leurres in [
    (0, 12, 18, X**6 - 1, [X**36 - 1, X - 1, X**3 - 1]),
    (1, 15, 25, X**5 - 1, [X**75 - 1, X**10 - 1, X - 1]),
    (2, 14, 9, X - 1, [X**126 - 1, X**5 - 1, X**2 - 1]),
]:
    assert gcd(X**a - 1, X**b - 1) == bonne, v
    for l in leurres:
        assert gcd(X**a - 1, X**b - 1) != l, (v, l)
assert [14 % 9, 9 % 5, 5 % 4, 4 % 1] == [5, 4, 1, 0]  # restes successifs cités en v2

# arith-l1-polynomes-bezout-ex4 (et v1, v2) : n tels que le PGCD ne soit pas 1 (contrôle n = 1..48)
def ensemble(D):
    return {k for k in range(1, 49) if gcd(D, X**k + 1) != 1}


assert ensemble(X**2 + X + 1) == set()                                  # énoncé : aucun n
assert ensemble(X**2 + 1) == {k for k in range(1, 49) if k % 4 == 2}     # v1 : n ≡ 2 mod 4
assert ensemble(X**2 - X + 1) == {k for k in range(1, 49) if k % 6 == 3}  # v2 : n ≡ 3 mod 6
# leurres de l'énoncé, v1, v2 : ensembles faux (non vides alors que la réponse est vide, ou différents)
N = range(1, 49)
for D, faux in [
    (X**2 + X + 1, [{k for k in N if k % 2 == 1}, {k for k in N if k % 3 == 0}, {k for k in N if k % 3 == 1}]),
    (X**2 + 1, [{k for k in N if k % 2 == 0}, {k for k in N if k % 4 == 0}, set()]),
    (X**2 - X + 1, [{k for k in N if k % 2 == 1}, {k for k in N if k % 3 == 0}, set()]),
]:
    for F in faux:
        assert ensemble(D) != F
assert egal(X**6 + 1, (X**2 + 1) * (X**4 - X**2 + 1)) and egal(X**3 + 1, (X + 1) * (X**2 - X + 1))
assert egal(exp(I * pi / 3) + exp(-I * pi / 3), 1) and egal(exp(I * pi / 3) * exp(-I * pi / 3), 1)

# ===========================================================================
# arith-l1-polynomes-lagrange (carte)
# ===========================================================================
def lagrange(xs, ys):
    tot = 0
    for i, xi in enumerate(xs):
        L = 1
        for j, xj in enumerate(xs):
            if j != i:
                L *= (X - xj) / sympy.Integer(xi - xj)
        tot += ys[i] * L
    return expand(tot)


def L_i(xs, i):
    L = 1
    for j, xj in enumerate(xs):
        if j != i:
            L *= (X - xj) / sympy.Integer(xs[i] - xj)
    return expand(L)


assert lagrange([-1, 0, 1], [0, 1, 0]) == 1 - X**2  # exemple
assert L_i([-1, 0, 1], 0) == expand(X * (X - 1) / 2) and L_i([-1, 0, 1], 2) == expand(X * (X + 1) / 2)
assert expand(sum(L_i([-1, 0, 1], i) for i in range(3))) == 1
assert lagrange([0, 2], [0, 1]) == X / 2  # contre-exemple Z : X/2 ∉ Z[X]
assert egal(X * (X - 1), 0) is False and [(X * (X - 1)).subs(X, t) for t in (0, 1)] == [0, 0]  # degré trop grand
# affirmation 2 (fausse) : P = (X^2 + X + 2)/2 vérifie P(0)=1, P(1)=2, P(2)=4 et P(3)=7
Paff = (X**2 + X + 2) / 2
assert [Paff.subs(X, t) for t in range(4)] == [1, 2, 4, 7]
# affirmation 3 (fausse) : X + 1 et X + 1 + X(X - 1)
assert [(X + 1 + X * (X - 1)).subs(X, t) for t in (0, 1)] == [1, 2]
# trous : « de degré exactement n » faux (valeurs constantes -> P constant)
assert lagrange([0, 1, 2], [5, 5, 5]) == 5

# arith-l1-polynomes-lagrange-ex1 (et v1, v2)
for v, xs, ys, P, leurres in [
    (0, [0, 1, 2], [1, 3, 7], X**2 + X + 1, [2 * X + 1, X**2 + 2 * X + 1, 2 * X**2 + 1]),
    (1, [-1, 0, 1], [2, 1, 4], 2 * X**2 + X + 1, [2 * X**2 - X + 1, X**2 + X + 2, 3 * X**2 + 1]),
    (2, [1, 2, 3], [0, 1, 4], X**2 - 2 * X + 1, [X**2 - 1, X**2 - 3 * X + 2, 2 * X**2 - 5 * X + 3]),
]:
    assert lagrange(xs, ys) == expand(P), v
    for l in leurres:
        assert [l.subs(X, t) for t in xs] != ys, (v, l)

# arith-l1-polynomes-lagrange-ex2 (et v1, v2) : L_i(4) pour les abscisses 0, 1, 2, 3
xs = [0, 1, 2, 3]
assert [L_i(xs, i).subs(X, 4) for i in range(4)] == [-1, 4, -6, 4]
for v, ys, val, leurres in [(0, [1, 0, -1, 4], 21, [13, -11, 9]), (1, [2, 1, 3, 0], -16, [16, -4, 8]),
                            (2, [1, 1, 2, 5], 11, [10, 8, 13])]:
    assert lagrange(xs, ys).subs(X, 4) == val and val not in leurres, v
    assert reponse("arith-l1-polynomes-lagrange-ex2", v) == f"${val}$"

# arith-l1-polynomes-lagrange-ex3 (et v1, v2) : Σ x_i^3 L_i = X^3 - N
for v, xs, val, leurres in [
    (0, [1, 2, 3], 6 * X**2 - 11 * X + 6, [X**3, -6 * X**2 + 11 * X - 6, 6 * X**2 + 11 * X + 6]),
    (1, [0, 1, -1], X, [X**3, -X, 0]),
    (2, [1, 2, 4], 7 * X**2 - 14 * X + 8, [X**3, -7 * X**2 + 14 * X - 8, 7 * X**2 - 14 * X - 8]),
]:
    S = lagrange(xs, [t**3 for t in xs])
    Nn = expand((X - xs[0]) * (X - xs[1]) * (X - xs[2]))
    assert S == expand(val) and egal(S, X**3 - Nn), v
    for l in leurres:
        assert not egal(S, l), (v, l)
    for k in range(3):  # question 1
        assert lagrange(xs, [t**k for t in xs]) == X**k

# arith-l1-polynomes-lagrange-ex4 (et v1, v2) : sommes alternées
def somme(nn, k):
    return sum((-1)**i * binomial(nn, i) * sympy.Integer(i)**k if i > 0 or k > 0 else (-1)**i * binomial(nn, i)
               for i in range(nn + 1))


for nn in range(1, 9):
    for k in range(0, nn + 1):
        attendu = (-1)**nn * factorial(nn) if k == nn else 0
        assert somme(nn, k) == attendu, (nn, k)
        # question 1 avec x_i = i : Σ x_i^k / w_i
        w = [sympy.prod([i - j for j in range(nn + 1) if j != i]) for i in range(nn + 1)]
        assert all(w[i] == (-1)**(nn - i) * factorial(i) * factorial(nn - i) for i in range(nn + 1))
for v, nn, k, val, leurres in [(0, 5, 5, -120, [120, 0, -1]), (1, 4, 4, 24, [-24, 0, 4]), (2, 6, 4, 0, [24, 720, -24])]:
    assert somme(nn, k) == val and val not in leurres, v
    assert reponse("arith-l1-polynomes-lagrange-ex4", v) == f"${val}$"
assert 0 - 4 + 6 * 16 - 4 * 81 + 256 == 24  # contrôle direct de v1

print("arith_l1_polynomes : toutes les vérifications passent.")
