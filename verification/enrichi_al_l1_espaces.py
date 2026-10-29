# Vérification calculatoire de l'enrichissement du chapitre « al-l1-espaces » (thème algebre-lineaire, niveau l1) :
# questions de content/qcm/al-l1-espaces.json et exercices ajoutés à
# content/exercices/algebre-lineaire/l1-al-l1-espaces.json (familles-ex5, dimension-ex6, inclusion-ex5, variantes
# « v1 », « v2 » comprises). Calcul exact avec SymPy. Chaque assert cite l'identifiant de la question ou de
# l'exercice. Voir CLAUDE.md, « Contrôle de l'exactitude ».
import sympy as sp
from sympy import Matrix as M, symbols, Abs, Poly, expand, zeros

X, a = symbols("X a")

# al-l1-espaces-q1 : f(1) = 2 f(0) stable par combinaison ; contre-exemples des autres
f, g = (lambda x: x), (lambda x: 1 - x)
assert f(0) * f(1) == 0 and g(0) * g(1) == 0 and (f(0) + g(0)) * (f(1) + g(1)) == 1
assert not (0 == 0 + 1)                      # la fonction nulle ne vérifie pas f(0) = f(1) + 1
l1_, l2_, u0, v0 = symbols("l1 l2 u0 v0")    # f(0) = u0, f(1) = 2 u0 ; g(0) = v0, g(1) = 2 v0
assert expand(l1_ * 2 * u0 + l2_ * 2 * v0 - 2 * (l1_ * u0 + l2_ * v0)) == 0

# al-l1-espaces-q2 (faux) : X^3 + X et -X^3 de degré 3, somme X de degré 1
assert Poly(X ** 3 + X, X).degree() == 3 and Poly(-X ** 3, X).degree() == 3 and Poly(X, X).degree() == 1

# al-l1-espaces-q3 : coordonnées (3, 3, 1) de X^2 + X + 1 dans (1, X-1, (X-1)^2) ; leurres faux
P = X ** 2 + X + 1
base = [1, X - 1, (X - 1) ** 2]
coord = lambda cs: expand(sum(c * b for c, b in zip(cs, base)))
assert coord((3, 3, 1)) == P
for l in ((1, 1, 1), (1, 3, 3), (3, 1, 1)):
    assert coord(l) != expand(P)
assert P.subs(X, 1) == 3 and sp.diff(P, X).subs(X, 1) == 3 and sp.diff(P, X, 2).subs(X, 1) / 2 == 1

# al-l1-espaces-q4 (faux) : (u-v) + (v-w) + (w-u) = 0
e = [M([1, 0, 0]), M([0, 1, 0]), M([0, 0, 1])]
assert (e[0] - e[1]) + (e[1] - e[2]) + (e[2] - e[0]) == zeros(3, 1)
assert M.hstack(e[0] - e[1], e[1] - e[2], e[2] - e[0]).rank() == 2

# al-l1-espaces-q5 : dimension de S_3(R) = 6 (rang des 9 coefficients des matrices symétriques génériques)
s = symbols("s0:6")
S = M([[s[0], s[3], s[4]], [s[3], s[1], s[5]], [s[4], s[5], s[2]]])
jac = M([[sp.diff(c, v) for v in s] for c in S])
assert jac.rank() == 6
assert 6 not in (3, 9, 5)

# al-l1-espaces-q6 (vrai) : (1, X+1, (X+1)^2, (X+1)^3) base de R_3[X]
coef = lambda p: [Poly(p, X).coeff_monomial(X ** k) for k in range(4)]
assert M([coef((X + 1) ** k) for k in range(4)]).rank() == 4

# al-l1-espaces-q7 : dim(F ∩ G) ∈ {1, 2, 3} (Grassmann) et chaque valeur atteinte
E = sp.eye(6)
F = [E.col(k) for k in range(4)]
def dim_inter(Fv, Gv):
    return len(Fv) + len(Gv) - M.hstack(*(Fv + Gv)).rank()
assert dim_inter(F, [E.col(3), E.col(4), E.col(5)]) == 1
assert dim_inter(F, [E.col(2), E.col(3), E.col(4)]) == 2
assert dim_inter(F, [E.col(0), E.col(1), E.col(2)]) == 3
assert 4 + 3 - 6 == 1                          # minoration ; 0 impossible, donc les leurres sont faux

# al-l1-espaces-q8 : seule (1, 1, 1) n'est pas dans le plan x - y + 2z = 0
plan = lambda v: v[0] - v[1] + 2 * v[2]
assert plan((1, 1, 1)) == 2
assert all(plan(v) == 0 for v in ((1, 1, 0), (2, 0, -1), (0, 2, 1)))

# al-l1-espaces-familles-ex5 : (|x - a_k|)_k et 1 libres ; contrôle par la matrice des valeurs en n+1 points
x = symbols("x", real=True)
for n in range(1, 6):
    fs = [Abs(x - k) for k in range(1, n + 1)] + [sp.Integer(1)]
    pts = [sp.Rational(k, 2) for k in range(0, 2 * n + 3, 2)][: n + 1]  # 0, 1, ..., n
    V = M([[fk.subs(x, p) for fk in fs] for p in pts])
    assert V.rank() == n + 1
    # taux d'accroissement en a_k : |t|/t vaut 1 à droite et -1 à gauche
    t = symbols("t")
    assert sp.limit(Abs(t) / t, t, 0, "+") == 1 and sp.limit(Abs(t) / t, t, 0, "-") == -1

# al-l1-espaces-dimension-ex6 (exercice, v1, v2) : relations, dimension et base extraite
cas = [
    ([(1, 0, 1, 2), (0, 1, 1, 1), (1, 1, 2, 3), (2, -1, 1, 3)], 2, [0, 1], [(2, (1, 1, 0, 0)), (3, (2, -1, 0, 0))]),
    ([(1, 1, 0, 0), (0, 1, 1, 0), (1, 0, -1, 0), (0, 0, 1, 1)], 3, [0, 1, 3], [(2, (1, -1, 0, 0))]),
    ([(1, 2, 0, -1), (2, 4, 0, -2), (0, 1, 1, 1), (1, 3, 1, 0)], 2, [0, 2], [(1, (2, 0, 0, 0)), (3, (1, 0, 1, 0))]),
]
for u, dimF, base_idx, relations in cas:
    U = [M(v) for v in u]
    assert M.hstack(*U).rank() == dimF
    assert M.hstack(*[U[k] for k in base_idx]).rank() == len(base_idx) == dimF
    for k, cs in relations:
        assert U[k] == sum((c * U[i] for i, c in enumerate(cs)), zeros(4, 1))
# leurres de dimension : toute autre valeur que 2, 3, 2 est fausse (contrôlé par le rang ci-dessus)

# al-l1-espaces-inclusion-ex5 (exercice, v1, v2) : H hyperplan, base de H, liberté, valeur de a, dim(F_a ∩ H)
cas = [
    ((1, 2, -1, 1), [(1, 0, 1, 0), (0, 1, 2, 0)], lambda a: (1, 1, a, -1), 2,
     [(1, 0, 1, 0), (0, 1, 2, 0), (0, 0, 1, 1)], [-2]),
    ((1, -1, 2, -1), [(1, 1, 0, 0), (0, 2, 1, 0)], lambda a: (1, 0, a, 3), 1,
     [(1, 0, 0, 1), (0, 1, 0, -1), (0, 0, 1, 2)], [-1]),
    ((2, 1, 1, -3), [(1, -2, 0, 0), (0, 1, -1, 0)], lambda a: (2, 0, a, 1), -1,
     [(1, -2, 0, 0), (0, -1, 1, 0), (0, 3, 0, 1)], [1]),
]
for n_h, (u1, u2), u3, aval, base_h, a_leurres in cas:
    eqH = lambda v: sum(c * x_ for c, x_ in zip(n_h, v))
    assert eqH(u1) == 0 and eqH(u2) == 0
    assert all(eqH(v) == 0 for v in base_h) and M.hstack(*[M(v) for v in base_h]).rank() == 3   # dim H = 3
    assert sp.solve(eqH(u3(a)), a) == [aval]
    Fa = M.hstack(M(u1), M(u2), M(u3(a)))
    for av in [aval, aval + 1, aval - 3, 0, 7]:
        assert Fa.subs(a, av).rank() == 3                                                        # liberté
    # a = aval : F_a ⊂ H, même dimension, donc F_a = H (rang de la réunion = 3)
    assert M.hstack(Fa.subs(a, aval), *[M(v) for v in base_h]).rank() == 3
    for av in a_leurres + [aval + 1, 0]:
        if av == aval:
            continue
        R = M.hstack(Fa.subs(a, av), *[M(v) for v in base_h]).rank()
        assert R == 4                                    # F_a + H = R^4, donc F_a ≠ H
        assert 3 + 3 - R == 2                            # dim(F_a ∩ H) = 2
    # les propositions « pour tout a » et « aucune valeur » sont fausses
    assert eqH(u3(aval + 1)) != 0 and eqH(u3(aval)) == 0

# --- Relecture indépendante (2026-10-29) : réponses retrouvées à l'aveugle, recalculées ici.
from sympy import Matrix as _M, symbols as _sym, expand as _expand
# al-l1-espaces-q3 : (X-1)² + 3(X-1) + 3 = X² + X + 1
_X = _sym('X')
assert _expand((_X - 1) ** 2 + 3 * (_X - 1) + 3) == _X ** 2 + _X + 1
# al-l1-espaces-q8 : (1,1,1) hors du plan x - y + 2z = 0, les trois autres vecteurs dedans
assert [v[0] - v[1] + 2 * v[2] for v in [(1, 1, 1), (1, 1, 0), (2, 0, -1), (0, 2, 1)]] == [2, 0, 0, 0]
# al-l1-espaces-dimension-ex6 (exercice, v1, v2) : rang de la famille et de la base extraite
for _u, _base, _r in [([(1, 0, 1, 2), (0, 1, 1, 1), (1, 1, 2, 3), (2, -1, 1, 3)], [0, 1], 2),
                      ([(1, 1, 0, 0), (0, 1, 1, 0), (1, 0, -1, 0), (0, 0, 1, 1)], [0, 1, 3], 3),
                      ([(1, 2, 0, -1), (2, 4, 0, -2), (0, 1, 1, 1), (1, 3, 1, 0)], [0, 2], 2)]:
    assert _M(_u).rank() == _r and _M([_u[i] for i in _base]).rank() == _r
# al-l1-espaces-inclusion-ex5 (exercice, v1, v2) : F_a = H exactement pour a = 2, 1, -1 ; sinon dim(F_a ∩ H) = 2
_a = _sym('a')
for _h, _u1, _u2, _u3, _bon in [((1, 2, -1, 1), (1, 0, 1, 0), (0, 1, 2, 0), (1, 1, _a, -1), 2),
                                ((1, -1, 2, -1), (1, 1, 0, 0), (0, 2, 1, 0), (1, 0, _a, 3), 1),
                                ((2, 1, 1, -3), (1, -2, 0, 0), (0, 1, -1, 0), (2, 0, _a, 1), -1)]:
    _H = _M([_h])
    for _val in [_bon, _bon + 3]:
        _F = _M([_u1, _u2, [c.subs(_a, _val) if hasattr(c, 'subs') else c for c in _u3]])
        assert _F.rank() == 3
        _dimsomme = _M.vstack(_F, _M(_H.nullspace()[0].T), _M(_H.nullspace()[1].T), _M(_H.nullspace()[2].T)).rank()
        assert (_dimsomme == 3) == (_val == _bon) and 3 + 3 - _dimsomme == (3 if _val == _bon else 2)

print("enrichi_al_l1_espaces : tout est vérifié.")
