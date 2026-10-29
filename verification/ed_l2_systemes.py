"""Vérifications calculatoires du chapitre « ed-l2-systemes » (équations différentielles, L2) :
cartes, exercices et QCM de content/*/equations-differentielles/l2-ed-l2-systemes.json et content/qcm/ed-l2-systemes.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question.
Lancer : python3 verification/ed_l2_systemes.py (doit finir sans erreur, en quelques secondes).
"""

from sympy import (E, Matrix, Piecewise, Rational, cos, cosh, diag, diff, eye, exp, factor, integrate, log, pi,
                   simplify, sin, sinh, sqrt, symbols, zeros)

t, s = symbols("t s", real=True)
R_ = Rational


def nul(expr):
    return simplify(expr) == 0


def mat_nul(M):
    return all(nul(x) for x in M)


def est_solution(X, A, B=None):
    """X' - A X - B = 0 (vecteurs colonnes SymPy)."""
    B = B if B is not None else zeros(X.shape[0], 1)
    return mat_nul(X.diff(t) - A * X - B)


def num_egal(f, g, pts=(R_(1, 3), R_(7, 5), R_(5, 2))):
    return all(abs(float((f - g).subs(t, p))) < 1e-12 for p in pts)


# ---------------------------------------------------------------------------
# ed-l2-systemes-cauchy (carte)
# ---------------------------------------------------------------------------
A_rot = Matrix([[0, 1], [-1, 0]])
# Exemple : (cos t, -sin t) solution de x' = y, y' = -x, valeur (1, 0) en 0.
assert est_solution(Matrix([cos(t), -sin(t)]), A_rot)
# Affirmation 5 : (cos + sin, cos - sin) est aussi solution, de première composante 1 en 0, mais X(0) = (1, 1).
X2 = Matrix([cos(t) + sin(t), cos(t) - sin(t)])
assert est_solution(X2, A_rot) and X2.subs(t, 0) == Matrix([1, 1])
# Contre-exemple : t x' = 2x, a t^2 (t >= 0), b t^2 (t < 0) : solution sur chaque demi-droite, dérivée nulle en 0.
assert nul(t * diff(t**2, t) - 2 * t**2)
# Contre-exemple : x' = x^2, x(0) = 1 : x = 1/(1 - t).
assert nul(diff(1 / (1 - t), t) - (1 / (1 - t)) ** 2)
# Affirmation 2 : x = t|t|/2 a pour dérivée |t| (calcul sur chaque demi-droite) ; x'' vaut 1 à droite, -1 à gauche.
assert nul(diff(t**2 / 2, t) - t) and nul(diff(-t**2 / 2, t) + t)
# Piège 2 : a(t) = 2/(t-1) (t != 1), a(1) = 0 sur [0, 2] : x = (t-1)^2 et x = d (t-1)^2 conviennent, x(1) = x'(1) = 0.
assert nul(diff((t - 1) ** 2, t) - 2 / (t - 1) * (t - 1) ** 2)
# Exemple 2 : t x' = x + t^2 y, t y' = -y se normalise en A(t) = [[1/t, t], [0, -1/t]].
x_, y_ = symbols("x y")
assert nul((x_ + t**2 * y_) / t - (x_ / t + t * y_))

# ed-l2-systemes-cauchy-ex1 (+ v1, v2) : domaines des coefficients normalisés.
# Coefficients 1/(t-2), t/(t-2), ln(t+1), 1/(t+3) : définis pour t > -1, t != 2 ; t0 = 1 -> ]-1, 2[, t0 = 5 -> ]2, +oo[.
assert -1 < 1 < 2 and 5 > 2 and -3 < -1  # ]-3,2[ déborde du domaine du logarithme
# v2 : 2/(t+1), -1/(t+1), 1/(4-t), e^t/(4-t) : t != -1, 4 ; t0 = 0 -> ]-1, 4[.
assert -1 < 0 < 4

# ed-l2-systemes-cauchy-ex2 (+ v1, v2) : x' = x/t + y, y' = 2y/t ; (x(1), y(1)) -> x(2).
A_c = Matrix([[1 / t, 1], [0, 2 / t]])
for (x0, y0), xs, x2, leurres in [
    ((1, 1), t**3 / 2 + t / 2, 5, [4, R_(9, 2), 6]),
    ((0, 2), t**3 - t, 6, [8, 7, 4]),
    ((2, -1), -t**3 / 2 + R_(5, 2) * t, 1, [-4, 0, 3]),
]:
    X = Matrix([xs, y0 * t**2])
    assert est_solution(X, A_c) and X.subs(t, 1) == Matrix([x0, y0])
    assert xs.subs(t, 2) == x2 and x2 not in leurres and len(set(leurres)) == 3

# ed-l2-systemes-cauchy-ex3 : solutions explicites de x' = sin(t) x, cos(t) x, (1 + cos t) x.
assert nul(diff(exp(-cos(t)), t) - sin(t) * exp(-cos(t)))
assert nul(diff(exp(sin(t)), t) - cos(t) * exp(sin(t)))
assert nul(diff(exp(t + sin(t)), t) - (1 + cos(t)) * exp(t + sin(t)))
assert nul(exp(2 * pi + sin(2 * pi)) - exp(2 * pi) * exp(sin(0)))  # x(2 pi) = e^{2 pi} x(0)

# ed-l2-systemes-cauchy-ex4 (+ v1, v2) : t x' = k x avec k = 2, 1, -1 ; dimensions 2, 1, 0.
for k, f in [(2, t**2), (1, t), (-1, 1 / t)]:
    assert nul(t * diff(f, t) - k * f)
# k = 2 : dérivée en 0 de a t^2 / b t^2 nulle des deux côtés ; k = 1 : dérivées a et b (raccord a = b) ;
# k = -1 : a/t non borné en 0 sauf a = 0.
assert diff(t**2, t).subs(t, 0) == 0 and diff(t, t) == 1

# ---------------------------------------------------------------------------
# ed-l2-systemes-dimension (carte)
# ---------------------------------------------------------------------------
A_e = Matrix([[0, 1], [-2 / t**2, 2 / t]])
X1, X2 = Matrix([t, 1]), Matrix([t**2, 2 * t])
assert est_solution(X1, A_e) and est_solution(X2, A_e)
W = Matrix.hstack(X1, X2).det()
assert nul(W - t**2)
assert nul(W.subs(t, 1) * exp(integrate(A_e.trace().subs(t, s), (s, 1, t))) - t**2)  # Liouville
# Contre-exemples : W((1,t),(t,t^2)) = 0 ; W((1,0),(0,t)) = t.
assert Matrix([[1, t], [t, t**2]]).det() == 0
assert Matrix([[1, 0], [0, t]]).det() == t
# Trou 3 : n = 1, A = 1, X1 = e^t, W(t) = e^t, t0 = 0 ; les trois leurres donnent e^{-t}, 1 + t, 1 en t = 1.
Wv = exp(t)
for leurre in [exp(-t), 1 + t, exp(1 - 1) * 1]:
    assert not nul((Wv - leurre).subs(t, 1))
# Affirmation 3 : A = [[0,1],[0,0]] de trace nulle, X = (t, 1) solution non bornée.
N2 = Matrix([[0, 1], [0, 0]])
assert est_solution(Matrix([t, 1]), N2) and N2.trace() == 0

# ed-l2-systemes-dimension-ex1 (+ v1, v2) : A(t) = [[0,1],[-k/t^2, k/t]] (k = 3, 4) ou [[0,1],[-1/t^2, 1/t]].
for A, Y2, Wbon, leurres in [
    (Matrix([[0, 1], [-3 / t**2, 3 / t]]), Matrix([t**3, 3 * t**2]), 2 * t**3, [4 * t**3, -2 * t**3, 3 * t**2]),
    (Matrix([[0, 1], [-1 / t**2, 1 / t]]), Matrix([t * log(t), log(t) + 1]), t, [t * log(t), -t, 2 * t]),
    (Matrix([[0, 1], [-4 / t**2, 4 / t]]), Matrix([t**4, 4 * t**3]), 3 * t**4, [5 * t**4, -3 * t**4, 4 * t**3]),
]:
    assert est_solution(Matrix([t, 1]), A) and est_solution(Y2, A)
    Wc = Matrix.hstack(Matrix([t, 1]), Y2).det()
    assert nul(Wc - Wbon)
    assert nul(Wc.subs(t, 1) * exp(integrate(A.trace().subs(t, s), (s, 1, t))) - Wbon)
    assert all(not nul(Wbon - l) for l in leurres)

# ed-l2-systemes-dimension-ex2 (+ v1, v2) : Liouville.
for tr, t0, w0, t1, bon, leurres in [
    (2 * t, 0, 2, 1, 2 * E, [2 * E**2, 2, E]),
    (2 * cos(t), 0, 3, pi / 2, 3 * E**2, [3 * E, 3, 3 * exp(pi)]),
    (2 / (1 + t), 0, -1, 1, -4, [-2, 4, -E**2]),
]:
    val = w0 * exp(integrate(tr.subs(t, s), (s, t0, t1)))
    assert nul(val - bon) and all(not nul(bon - l) for l in leurres)
assert Matrix([[t, exp(t)], [1, t]]).trace() == 2 * t
assert Matrix([[cos(t), t**2], [-1, cos(t)]]).trace() == 2 * cos(t)

# ed-l2-systemes-dimension-ex3 (+ v1, v2) : wronskiens nuls, familles libres ; A = R' R^{-1}.
assert Matrix([[1, t], [t**2, t**3]]).det() == 0
assert nul(Matrix([[1, exp(t)], [exp(t), exp(2 * t)]]).det())
for R, Abon, leurres in [
    (diag(1, t), Matrix([[0, 0], [0, 1 / t]]), [Matrix([[0, 0], [0, t]]), Matrix([[0, 0], [1 / t, 0]]), Matrix([[1 / t, 0], [0, 0]])]),
    (diag(1, t**2), Matrix([[0, 0], [0, 2 / t]]), [Matrix([[0, 0], [0, 2 * t]]), Matrix([[0, 0], [0, 1 / t**2]]), Matrix([[2 / t, 0], [0, 0]])]),
    (diag(t, 1), Matrix([[1 / t, 0], [0, 0]]), [Matrix([[0, 0], [0, 1 / t]]), Matrix([[t, 0], [0, 0]]), Matrix([[1 / t, 0], [0, 1]])]),
]:
    assert mat_nul(R.diff(t) * R.inv() - Abon)
    for L in leurres:  # un leurre ne rend pas les deux colonnes solutions
        assert not (est_solution(R[:, 0], L) and est_solution(R[:, 1], L))

# ed-l2-systemes-dimension-ex4 : trace de [[1/(1+t), sin t], [-sin t, 0]] et ∫ = ln(1 + t).
assert nul(integrate(1 / (1 + s), (s, 0, t)) - log(1 + t))

# ---------------------------------------------------------------------------
# ed-l2-systemes-variation-constantes (carte)
# ---------------------------------------------------------------------------
A_n = Matrix([[1, 1], [0, 1]])
Rf = exp(t) * Matrix([[1, t], [0, 1]])
assert mat_nul(Rf.diff(t) - A_n * Rf) and nul(Rf.det() - exp(2 * t))
B_ex = Matrix([0, exp(t)])
Cp = (Rf.inv() * B_ex).applyfunc(simplify)
assert mat_nul(Cp - Matrix([-t, 1]))
Xp = Rf * Matrix([-t**2 / 2, t])
assert mat_nul(Xp - exp(t) * Matrix([t**2 / 2, t])) and est_solution(Xp, A_n, B_ex)
# Trous et contre-exemple : n = 1, x' = x + 1, R(t) = 2 e^t, t0 = 0, X0 = 1 ; solution 2e^t - 1.
Rs = 2 * exp(t)
bonne = Rs * (R_(1, 2) + integrate(exp(-s) / 2, (s, 0, t)))
assert nul(bonne - (2 * exp(t) - 1)) and bonne.subs(t, 0) == 1 and nul(diff(bonne, t) - bonne - 1)
leurres = [Rs * (1 + integrate(exp(-s) / 2, (s, 0, t))),
           Rs * R_(1, 2) * 1 + integrate(exp(-s) / 2, (s, 0, t)),
           Rs * (R_(1, 2) + integrate(2 * exp(s), (s, 0, t)))]
for L in leurres:
    assert not (L.subs(t, 0) == 1 and nul(diff(L, t) - L - 1))
# Trou 1 (R = 2e^t, x' = x + 1) : C' = B, R'C = B, C' = R B ne donnent pas de solution.
for C in [t, exp(-t) / 2, integrate(2 * exp(s), (s, 0, t))]:
    assert not nul(diff(Rs * C, t) - Rs * C - 1)
# Affirmation 3 : C = -e^{-t} et 1 - e^{-t} avec R = e^t donnent -1 et e^t - 1.
for C in [-exp(-t), 1 - exp(-t)]:
    assert nul(diff(exp(t) * C, t) - exp(t) * C - 1)

# ed-l2-systemes-variation-constantes-ex1 (+ v1, v2) : A = diag(1, 2), X(0) = 0.
A_d = diag(1, 2)
for B, Xb, leurres in [
    (Matrix([exp(t), 0]), Matrix([t * exp(t), 0]),
     [Matrix([exp(t) - 1, 0]), Matrix([t * exp(2 * t), 0]), Matrix([exp(2 * t) - exp(t), 0])]),
    (Matrix([0, exp(2 * t)]), Matrix([0, t * exp(2 * t)]),
     [Matrix([0, exp(2 * t) - 1]), Matrix([0, t * exp(t)]), Matrix([0, (exp(2 * t) - 1) / 2])]),
    (Matrix([exp(2 * t), 0]), Matrix([exp(2 * t) - exp(t), 0]),
     [Matrix([t * exp(2 * t), 0]), Matrix([exp(2 * t), 0]), Matrix([t * exp(t), 0])]),
]:
    assert est_solution(Xb, A_d, B) and Xb.subs(t, 0) == zeros(2, 1)
    for L in leurres:
        assert not (est_solution(L, A_d, B) and L.subs(t, 0) == zeros(2, 1))

# ed-l2-systemes-variation-constantes-ex2 (+ v1, v2) : R = [[t, t^2], [1, 2t]], X(1) = 0.
Rt = Matrix([[t, t**2], [1, 2 * t]])
assert mat_nul(Rt.diff(t) - A_e * Rt) and nul(Rt.det() - t**2)
for b, Cp_att, xb, leurres in [
    (t, [-t, 1], t * (t - 1) ** 2 / 2, [t**3 / 2, t * (t**2 - 1) / 2, t**2 * (t - 1) / 2]),
    (t**2, [-t**2, t], t * (t - 1) ** 2 * (t + 2) / 6, [t**4 / 6, t * (t - 1) ** 2 / 2, t**2 * (t - 1) * (t + 2) / 6]),
    (1, [-1, 1 / t], t**2 * log(t) - t**2 + t, [t**2 * log(t), t**2 * log(t) + t**2 - t, t * log(t)]),
]:
    B = Matrix([0, b])
    assert mat_nul((Rt.inv() * B) - Matrix(Cp_att))
    X = Matrix([xb, diff(xb, t)])
    assert est_solution(X, A_e, B) and mat_nul(X.subs(t, 1))
    for L in leurres:  # même problème de Cauchy : unicité, donc un leurre en diffère
        XL = Matrix([L, diff(L, t)])
        assert not (est_solution(XL, A_e, B) and mat_nul(XL.subs(t, 1)))

# ed-l2-systemes-variation-constantes-ex3 (+ v1, v2) : X' = J X + B, X(0) = 0.
J = Matrix([[0, -1], [1, 0]])
Rr = Matrix([[cos(t), -sin(t)], [sin(t), cos(t)]])
assert mat_nul(Rr.diff(t) - J * Rr) and nul(Rr.det() - 1)
for B, Xb, leurres in [
    (Matrix([cos(t), sin(t)]), Matrix([t * cos(t), t * sin(t)]),
     [Matrix([sin(t), 1 - cos(t)]), Matrix([t * sin(t), -t * cos(t)]), Matrix([cos(t) - 1, sin(t)])]),
    (Matrix([-sin(t), cos(t)]), Matrix([-t * sin(t), t * cos(t)]),
     [Matrix([t * cos(t), t * sin(t)]), Matrix([cos(t) - 1, sin(t)]), Matrix([-sin(t), cos(t) - 1])]),
    (Matrix([1, 0]), Matrix([sin(t), 1 - cos(t)]),
     [Matrix([t, 0]), Matrix([1 - cos(t), sin(t)]), Matrix([t * cos(t), t * sin(t)])]),
]:
    C = (Rr.T * B).applyfunc(lambda f: integrate(simplify(f).subs(t, s), (s, 0, t)))
    assert mat_nul(Rr * C - Xb)
    assert est_solution(Xb, J, B) and mat_nul(Xb.subs(t, 0))
    for L in leurres:
        assert not (est_solution(L, J, B) and mat_nul(L.subs(t, 0)))

# ed-l2-systemes-variation-constantes-ex4 (+ v1, v2) : t X' = M X + (t, 0), X(1) = 0.
for M, (l1, l2), Xp_att, xb, yb, leurres in [
    (Matrix([[1, 2], [2, 1]]), (3, -1), Matrix([0, -t / 2]), (t**4 - 1) / (4 * t), (t**2 - 1) ** 2 / (4 * t),
     [(t**2 - 1) ** 2 / (4 * t), (t**3 - t) / 4, (t**4 - 1) / (2 * t)]),
    (Matrix([[1, 1], [1, 1]]), (2, 0), Matrix([0, -t]), (t**2 - 1) / 2, (t - 1) ** 2 / 2,
     [(t - 1) ** 2 / 2, t**2 - 1, (t**2 - t) / 2]),
    (Matrix([[3, 1], [1, 3]]), (4, 2), Matrix([-2 * t / 3, t / 3]), (t**4 + 3 * t**2 - 4 * t) / 6,
     t**4 / 6 - t**2 / 2 + t / 3, [(t**4 - 3 * t**2 + 2 * t) / 6, (t**4 + 3 * t**2 - 4 * t) / 3, (t**4 - t**2) / 6]),
]:
    V1, V2 = Matrix([1, 1]), Matrix([1, -1])
    assert M * V1 == l1 * V1 and M * V2 == l2 * V2
    Rm = Matrix.hstack(t**l1 * V1, t**l2 * V2)
    Bn = Matrix([1, 0])
    assert est_solution(t**l1 * V1, M / t) and est_solution(t**l2 * V2, M / t)
    assert est_solution(Xp_att, M / t, Bn)
    X = Matrix([xb, yb])
    assert est_solution(X, M / t, Bn) and mat_nul(X.subs(t, 1))
    assert all(not nul(xb - L) for L in leurres)
    # Variation des constantes : R C' = B.
    Cpr = (Rm.inv() * Bn).applyfunc(simplify)
    Cint = Cpr.applyfunc(lambda f: integrate(f, t))
    assert mat_nul((Rm * Cint) - Xp_att)

# ---------------------------------------------------------------------------
# ed-l2-systemes-exponentielle (carte)
# ---------------------------------------------------------------------------
assert mat_nul((t * N2).exp() - Matrix([[1, t], [0, 1]]))
assert mat_nul((t * J).exp() - Rr)
A_, B_ = Matrix([[0, 1], [0, 0]]), Matrix([[0, 0], [1, 0]])
eAeB, eBeA = A_.exp() * B_.exp(), B_.exp() * A_.exp()
assert eAeB == Matrix([[2, 1], [1, 1]]) and eBeA == Matrix([[1, 1], [1, 2]])
assert (A_ + B_) ** 2 == eye(2)
assert mat_nul((A_ + B_).exp() - Matrix([[cosh(1), sinh(1)], [sinh(1), cosh(1)]]))
assert not mat_nul((A_ + B_).exp() - eAeB)
# Piège / trou 3 « diagonalisables » : A = diag(1, 0), B = [[0,1],[0,1]] (valeurs propres 0, 1 distinctes).
Ad, Bd = diag(1, 0), Matrix([[0, 1], [0, 1]])
assert Ad * Bd != Bd * Ad and not mat_nul((Ad + Bd).exp() - Ad.exp() * Bd.exp())
# Trou 3 « A inversible » : A = diag(1, 2), B = N.
Ai = diag(1, 2)
assert Ai * N2 != N2 * Ai and not mat_nul((Ai + N2).exp() - Ai.exp() * N2.exp())
# Trou 1 : pour A = I, sum (-1)^k A^k / k! = e^{-1} I != e I ; sum A^k et sum A^k / k divergent (termes 1, 1/k).
assert abs(float(sum(R_((-1) ** k, 1) / __import__("math").factorial(k) for k in range(30))) - float(exp(-1))) < 1e-12
assert not nul(exp(-1) - E)
# Trou 2 (A = I, t0 = 1, X0 = 1, n = 1) : e^{t-1} ; leurres e^t, e^{1-t}, 1 + (t - 1).
for L in [exp(t), exp(1 - t), 1 + (t - 1)]:
    assert not nul(L - exp(t - 1))
# Contre-exemple « coefficients variables » : A(t) = [[1, t], [0, 0]], X(0) = (0, 1).
At = Matrix([[1, t], [0, 0]])
Xvrai = Matrix([exp(t) - 1 - t, 1])
assert est_solution(Xvrai, At) and Xvrai.subs(t, 0) == Matrix([0, 1])
Mt = Matrix([[t, t**2 / 2], [0, 0]])
assert mat_nul(Mt**2 - t * Mt)  # d'où exp(Mt) = I + (e^t - 1)/t Mt
Xfaux = (eye(2) + (exp(t) - 1) / t * Mt) * Matrix([0, 1])
assert mat_nul(Xfaux - Matrix([t * (exp(t) - 1) / 2, 1])) and not est_solution(Xfaux, At)
assert abs(float(E - 2 - (E - 1) / 2)) > 0.1
# Affirmation 3 : e^{2 pi J} = I.
assert mat_nul((2 * pi * J).exp() - eye(2))
# Affirmation 4 : det e^A = e^{tr A} sur un exemple.
Aex = Matrix([[1, 2], [3, -1]])
assert abs(float(Aex.exp().det() - exp(Aex.trace()))) < 1e-10

# ed-l2-systemes-exponentielle-ex1 (+ v1, v2) : e^{tN} e_3.
for N, Xb, leurres in [
    (Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]]), Matrix([t**2 / 2, t, 1]),
     [Matrix([t**2, t, 1]), Matrix([0, t, 1]), Matrix([exp(t**2 / 2), exp(t), 1])]),
    (Matrix([[0, 2, 0], [0, 0, 3], [0, 0, 0]]), Matrix([3 * t**2, 3 * t, 1]),
     [Matrix([6 * t**2, 3 * t, 1]), Matrix([2 * t, 3 * t, 1]), Matrix([0, 3 * t, 1])]),
    (Matrix([[0, 1, 1], [0, 0, 1], [0, 0, 0]]), Matrix([t + t**2 / 2, t, 1]),
     [Matrix([t**2 / 2, t, 1]), Matrix([t + t**2, t, 1]), Matrix([t, t, 1])]),
]:
    assert N**3 == zeros(3, 3)
    E3 = eye(3) + t * N + t**2 / 2 * N**2
    assert mat_nul(E3[:, 2] - Xb) and est_solution(Xb, N)
    assert all(not mat_nul(L - Xb) for L in leurres)

# ed-l2-systemes-exponentielle-ex2 (+ v1, v2) : A = l I + N, N^2 = 0.
for A, X0, Xb, leurres in [
    (Matrix([[2, 1], [0, 2]]), Matrix([1, 1]), exp(2 * t) * Matrix([1 + t, 1]),
     [exp(2 * t) * Matrix([1, 1]), Matrix([exp(2 * t) + exp(t), exp(2 * t)]), exp(2 * t) * Matrix([1 + 2 * t, 1])]),
    (Matrix([[-1, 3], [0, -1]]), Matrix([0, 1]), exp(-t) * Matrix([3 * t, 1]),
     [exp(-t) * Matrix([t, 1]), exp(-t) * Matrix([0, 1]), Matrix([exp(3 * t), exp(-t)])]),
    (Matrix([[1, -2], [0, 1]]), Matrix([2, -1]), exp(t) * Matrix([2 + 2 * t, -1]),
     [exp(t) * Matrix([2 - 2 * t, -1]), exp(t) * Matrix([2, -1]), exp(t) * Matrix([2 + t, -1])]),
]:
    assert mat_nul((t * A).exp() * X0 - Xb) and est_solution(Xb, A) and Xb.subs(t, 0) == X0
    for L in leurres:
        assert not (est_solution(L, A) and L.subs(t, 0) == X0)

# ed-l2-systemes-exponentielle-ex3 (+ v1, v2) : A^2 scalaire.
for A, Eb, leurres in [
    (Matrix([[0, -2], [2, 0]]), Matrix([[cos(2 * t), -sin(2 * t)], [sin(2 * t), cos(2 * t)]]),
     [Matrix([[cos(2 * t), sin(2 * t)], [-sin(2 * t), cos(2 * t)]]), Rr,
      Matrix([[cosh(2 * t), -sinh(2 * t)], [sinh(2 * t), cosh(2 * t)]])]),
    (Matrix([[0, 1], [1, 0]]), Matrix([[cosh(t), sinh(t)], [sinh(t), cosh(t)]]),
     [Matrix([[cos(t), sin(t)], [sin(t), cos(t)]]), Matrix([[1, exp(t)], [exp(t), 1]]),
      Matrix([[cosh(t), -sinh(t)], [-sinh(t), cosh(t)]])]),
    (Matrix([[0, -3], [3, 0]]), Matrix([[cos(3 * t), -sin(3 * t)], [sin(3 * t), cos(3 * t)]]),
     [Matrix([[cos(3 * t), sin(3 * t)], [-sin(3 * t), cos(3 * t)]]), Rr,
      Matrix([[cosh(3 * t), -sinh(3 * t)], [sinh(3 * t), cosh(3 * t)]])]),
]:
    assert (A**2 - A[1, 0] * A[0, 1] * eye(2)) == zeros(2, 2)
    assert mat_nul(Eb.diff(t) - A * Eb) and Eb.subs(t, 0) == eye(2)  # caractérise e^{tA}
    for L in leurres:
        assert not (mat_nul(L.diff(t) - A * L) and L.subs(t, 0) == eye(2))
assert mat_nul((pi * Matrix([[0, -2], [2, 0]])).exp() - eye(2))  # e^{(2pi/2) A} = I

# ed-l2-systemes-exponentielle-ex4 : A = [[-1,1],[0,-1]], ||e^{tA}||_oo = (1+t)e^{-t} <= (2/sqrt(e)) e^{-t/2}.
g = (1 + t) * exp(-t / 2)
assert nul(diff(g, t) - (1 - t) / 2 * exp(-t / 2)) and nul(g.subs(t, 1) - 2 / sqrt(E))
assert float(2 / sqrt(E)) > 1
assert mat_nul((t * Matrix([[-1, 1], [0, -1]])).exp() - exp(-t) * Matrix([[1, t], [0, 1]]))

# ed-l2-systemes-exponentielle-ex5 (+ v1, v2) : e^A e^B, e^B e^A, e^{A+B}.
for A, B, eAB, eBA, S2, leurres in [
    (Matrix([[0, 1], [0, 0]]), Matrix([[0, 0], [1, 0]]), Matrix([[2, 1], [1, 1]]), Matrix([[1, 1], [1, 2]]), eye(2),
     [Matrix([[1, 1], [1, 2]]), Matrix([[cosh(1), sinh(1)], [sinh(1), cosh(1)]]), Matrix([[1, 1], [1, 1]])]),
    (Matrix([[0, 2], [0, 0]]), Matrix([[0, 0], [1, 0]]), Matrix([[3, 2], [1, 1]]), Matrix([[1, 2], [1, 3]]), 2 * eye(2),
     [Matrix([[1, 2], [1, 3]]), Matrix([[2, 2], [1, 1]]), Matrix([[3, 1], [2, 1]])]),
    (Matrix([[0, 1], [0, 0]]), Matrix([[0, 0], [-1, 0]]), Matrix([[0, 1], [-1, 1]]), Matrix([[1, 1], [-1, 0]]), -eye(2),
     [Matrix([[1, 1], [-1, 0]]), Matrix([[cos(1), sin(1)], [-sin(1), cos(1)]]), Matrix([[1, 1], [-1, 1]])]),
]:
    assert A.exp() * B.exp() == eAB and B.exp() * A.exp() == eBA and (A + B) ** 2 == S2
    eS = (A + B).exp()
    assert not mat_nul(eS - eAB) and not mat_nul(eS - eBA)
    assert all(not mat_nul(L - eAB) for L in leurres)
eS1 = (Matrix([[0, 2], [1, 0]])).exp()
assert mat_nul(eS1 - (cosh(sqrt(2)) * eye(2) + sinh(sqrt(2)) / sqrt(2) * Matrix([[0, 2], [1, 0]])))
assert mat_nul(Matrix([[0, 1], [-1, 0]]).exp() - Matrix([[cos(1), sin(1)], [-sin(1), cos(1)]]))
# Question 3 : dérivées secondes en 0.
P_, Q_ = symbols("P Q", commutative=False)
assert (((P_ + Q_) ** 2).expand() - (P_**2 + 2 * P_ * Q_ + Q_**2) - (Q_ * P_ - P_ * Q_)).expand() == 0

# ---------------------------------------------------------------------------
# ed-l2-systemes-reduction (carte)
# ---------------------------------------------------------------------------
A12 = Matrix([[1, 2], [2, 1]])
assert A12 * Matrix([1, 1]) == 3 * Matrix([1, 1]) and A12 * Matrix([1, -1]) == -Matrix([1, -1])
Xex = (exp(3 * t) * Matrix([1, 1]) + exp(-t) * Matrix([1, -1])) / 2
assert est_solution(Xex, A12) and Xex.subs(t, 0) == Matrix([1, 0])
# Exemple 2 (cascade) : A = [[2,1],[0,2]], x = (a + bt) e^{2t}, y = b e^{2t}.
a_, b_ = symbols("a b")
assert est_solution(Matrix([(a_ + b_ * t) * exp(2 * t), b_ * exp(2 * t)]), Matrix([[2, 1], [0, 2]]))
# Trous : A = [[4,-2],[1,1]] = P D P^{-1}, P = [[1,2],[1,1]], D = diag(2,3) ; leurres faux.
A4 = Matrix([[4, -2], [1, 1]])
P4, D4 = Matrix([[1, 2], [1, 1]]), diag(2, 3)
assert P4 * D4 * P4.inv() == A4
eD = diag(exp(2 * t), exp(3 * t))
assert mat_nul((t * A4).exp() - P4 * eD * P4.inv())
for L in [P4.inv() * eD * P4, eD, P4 * eD]:
    assert not mat_nul((t * A4).exp() - L)
# Trou 2 : t e^{2t} V1, e^t * 2 V1, e^{2t} P^{-1} V1 = e^{2t} e_1 ne sont pas solutions.
V1 = P4[:, 0]
for L in [t * exp(2 * t) * V1, exp(t) * 2 * V1, exp(2 * t) * P4.inv() * V1]:
    assert not est_solution(L, A4)
# Contre-exemples : (t e^t, e^t) solution pour [[1,1],[0,1]] ; J et e^{it}(1, -i).
assert est_solution(Matrix([t * exp(t), exp(t)]), A_n)
assert est_solution(Matrix([cos(t), sin(t)]), J) and est_solution(Matrix([sin(t), -cos(t)]), J)
from sympy import I as Ic  # noqa: E402

assert J * Matrix([1, -Ic]) == Ic * Matrix([1, -Ic])
# Piège 1 : P = [[1,1],[0,0]] (colonnes propres de [[1,1],[0,1]]) non inversible.
assert Matrix([[1, 1], [0, 0]]).det() == 0
# Piège 2 : D = N triangulaire, P = I : e^{tD} = I + tN != I.
assert not mat_nul((t * N2).exp() - eye(2))

# ed-l2-systemes-reduction-ex1 (+ v1, v2) : cascades.
for A, X0, xb, leurres in [
    (Matrix([[1, 1], [0, 2]]), Matrix([0, 1]), exp(2 * t) - exp(t), [exp(2 * t), t * exp(t), exp(t) - exp(2 * t)]),
    (Matrix([[-1, 2], [0, 1]]), Matrix([3, 1]), exp(t) + 2 * exp(-t), [3 * exp(-t), exp(t) + 3 * exp(-t), 2 * exp(t) + exp(-t)]),
    (Matrix([[2, 3], [0, -1]]), Matrix([0, 1]), exp(2 * t) - exp(-t), [exp(-t) - exp(2 * t), 3 * t * exp(-t), exp(2 * t) - 3 * exp(-t)]),
]:
    X = (t * A).exp() * X0
    assert nul(X[0] - xb) and all(not nul(xb - L) for L in leurres)

# ed-l2-systemes-reduction-ex2 (+ v1, v2) : diagonalisation.
for A, X0, Xb, leurres in [
    (A4, Matrix([1, 0]), Matrix([2 * exp(3 * t) - exp(2 * t), exp(3 * t) - exp(2 * t)]),
     [exp(3 * t), 2 * exp(2 * t) - exp(3 * t), exp(3 * t) - exp(2 * t)]),
    (Matrix([[3, 1], [2, 2]]), Matrix([0, 3]), Matrix([exp(4 * t) - exp(t), exp(4 * t) + 2 * exp(t)]),
     [exp(t) - exp(4 * t), 3 * exp(4 * t), exp(4 * t) + 2 * exp(t)]),
    (Matrix([[0, 1], [6, 1]]), Matrix([2, 1]), Matrix([exp(3 * t) + exp(-2 * t), 3 * exp(3 * t) - 2 * exp(-2 * t)]),
     [2 * exp(3 * t), exp(3 * t) - exp(-2 * t), exp(2 * t) + exp(-3 * t)]),
]:
    assert est_solution(Xb, A) and Xb.subs(t, 0) == X0
    assert all(not nul(Xb[0] - L) for L in leurres)
assert Matrix([[3, 1], [2, 2]]) * Matrix([1, -2]) == Matrix([1, -2])
assert Matrix([[0, 1], [6, 1]]) * Matrix([1, 3]) == 3 * Matrix([1, 3])
assert Matrix([[0, 1], [6, 1]]) * Matrix([1, -2]) == -2 * Matrix([1, -2])

# ed-l2-systemes-reduction-ex3 (+ v1, v2) : valeur propre double, N^2 = 0.
for A, l, X0, Xb, leurres in [
    (Matrix([[3, 1], [-1, 1]]), 2, Matrix([1, 0]), exp(2 * t) * Matrix([1 + t, -t]),
     [exp(2 * t) * Matrix([1, 0]), exp(2 * t) * Matrix([1 - t, t]), Matrix([exp(3 * t), exp(3 * t) - exp(t)])]),
    (Matrix([[-2, 1], [-1, 0]]), -1, Matrix([0, 1]), exp(-t) * Matrix([t, 1 + t]),
     [exp(-t) * Matrix([0, 1]), exp(-t) * Matrix([-t, 1 - t]), exp(-t) * Matrix([t, 1])]),
    (Matrix([[5, -4], [1, 1]]), 3, Matrix([1, 1]), exp(3 * t) * Matrix([1 - 2 * t, 1 - t]),
     [exp(3 * t) * Matrix([1, 1]), exp(3 * t) * Matrix([1 + 2 * t, 1 + t]), exp(3 * t) * Matrix([1 - 2 * t, 1 - 2 * t])]),
]:
    lam = symbols("lam")
    assert factor((lam * eye(2) - A).det()) == (lam - l) ** 2
    N = A - l * eye(2)
    assert N**2 == zeros(2, 2) and N != zeros(2, 2)
    assert mat_nul((t * A).exp() - exp(l * t) * (eye(2) + t * N))
    assert est_solution(Xb, A) and Xb.subs(t, 0) == X0
    for L in leurres:
        assert not (est_solution(L, A) and L.subs(t, 0) == X0)

# ed-l2-systemes-reduction-ex4 (+ v1, v2) : matrices J - I et J (taille 3).
Jm = Matrix(3, 3, lambda i, j: 1)
for A, X0, Xb, leurres in [
    (Jm - eye(3), Matrix([1, 0, 0]), Matrix([exp(2 * t) + 2 * exp(-t), exp(2 * t) - exp(-t), exp(2 * t) - exp(-t)]) / 3,
     [(exp(2 * t) - exp(-t)) / 3, (2 * exp(2 * t) + exp(-t)) / 3, exp(-t)]),
    (Jm, Matrix([1, 0, 0]), Matrix([exp(3 * t) + 2, exp(3 * t) - 1, exp(3 * t) - 1]) / 3,
     [(exp(3 * t) - 1) / 3, exp(3 * t), (2 * exp(3 * t) + 1) / 3]),
    (Jm - eye(3), Matrix([0, 0, 3]), Matrix([exp(2 * t) - exp(-t), exp(2 * t) - exp(-t), exp(2 * t) + 2 * exp(-t)]),
     [exp(2 * t) + 2 * exp(-t), 3 * exp(2 * t), exp(-t) - exp(2 * t)]),
]:
    assert est_solution(Xb, A) and Xb.subs(t, 0) == X0
    assert all(not nul(Xb[0] - L) for L in leurres)
U = Matrix([1, 1, 1])
assert (Jm - eye(3)) * U == 2 * U and (Jm - eye(3)) * Matrix([1, -1, 0]) == -Matrix([1, -1, 0])

# ed-l2-systemes-reduction-ex5 (+ v1, v2) : conditions de stabilité.
lam, a = symbols("lam a", real=True)
# Principal : racines a ± i.
assert nul((lam * eye(2) - Matrix([[a, 1], [-1, a]])).det() - ((lam - a) ** 2 + 1))
assert est_solution(Matrix([cos(t), -sin(t)]), Matrix([[0, 1], [-1, 0]]))
# v1 : (lam + 1)^2 - a ; a = 1 : solution constante (1, 1).
assert nul((lam * eye(2) - Matrix([[-1, a], [1, -1]])).det() - ((lam + 1) ** 2 - a))
assert Matrix([[-1, 1], [1, -1]]) * Matrix([1, 1]) == zeros(2, 1)


def stable(M):
    return all(complex(v).real < 0 for v in M.eigenvals())


for av in [R_(-3), R_(-1, 2), 0, R_(1, 2), R_(9, 10)]:
    assert stable(Matrix([[-1, av], [1, -1]]))
for av in [1, R_(11, 10), 4]:
    assert not stable(Matrix([[-1, av], [1, -1]]))
# v2 : lam^2 + a lam + 2, stable ssi a > 0 (essais de part et d'autre de 2 sqrt 2).
assert nul((lam * eye(2) - Matrix([[0, 1], [-2, -a]])).det() - (lam**2 + a * lam + 2))
for av in [R_(1, 10), 1, 2, 3, 10]:
    assert stable(Matrix([[0, 1], [-2, -av]]))
for av in [0, R_(-1, 10), -3]:
    assert not stable(Matrix([[0, 1], [-2, -av]]))
for av in [R_(-1, 10), R_(-1, 2), -2]:
    assert stable(Matrix([[av, 1], [-1, av]]))
for av in [0, R_(1, 10)]:
    assert not stable(Matrix([[av, 1], [-1, av]]))
assert est_solution(Matrix([cos(sqrt(2) * t), -sqrt(2) * sin(sqrt(2) * t)]), Matrix([[0, 1], [-2, 0]]))

# Affirmations de la carte : (e^{3t} + e^{-t})/2 ; A = N non bornée.
assert nul(Xex[0] - (exp(3 * t) + exp(-t)) / 2)

# ---------------------------------------------------------------------------
# ed-l2-systemes-ordre2 (carte)
# ---------------------------------------------------------------------------
xp = cos(t) * log(cos(t)) + t * sin(t)
assert nul(diff(xp, t, 2) + xp - 1 / cos(t))
assert nul(-sin(t) / cos(t) - (-sin(t) * (1 / cos(t)) / 1)) and nul(cos(t) * (1 / cos(t)) - 1)  # λ1', λ2' (w = 1)
# Contre-exemple : t^2 x'' - 2t x' + 2x = 0 a pour solutions t et t^2.
for f in [t, t**2]:
    assert nul(t**2 * diff(f, t, 2) - 2 * t * diff(f, t) + 2 * f)
# Trou 2 et affirmation 3 : x'' + x' = 0, solutions 1 et e^{-t}, w = -e^{-t}, w' = -a w avec a = 1.
w = 1 * diff(exp(-t), t) - 0 * exp(-t)
assert nul(w + exp(-t)) and nul(diff(w, t) + w)
for L in [w, 0 * w, 0]:  # w' = a w, w' = -b w (b = 0), w' = 0 : faux
    assert not nul(diff(w, t) - L)
# Abaissement : x = x1 z, identité x'' + a x' + b x = x1 z'' + (2 x1' + a x1) z' quand x1 est solution.
af, bf, x1f, zf = (symbols("a"), symbols("b"), symbols("x1", cls=__import__("sympy").Function),
                   symbols("z", cls=__import__("sympy").Function))
aa, bb = __import__("sympy").Function("a")(t), __import__("sympy").Function("b")(t)
X1f, Zf = x1f(t), zf(t)
lhs = diff(X1f * Zf, t, 2) + aa * diff(X1f * Zf, t) + bb * X1f * Zf
rhs = X1f * diff(Zf, t, 2) + (2 * diff(X1f, t) + aa * X1f) * diff(Zf, t) + Zf * (diff(X1f, t, 2) + aa * diff(X1f, t) + bb * X1f)
assert nul(lhs - rhs)

# ed-l2-systemes-ordre2-ex1 (+ v1, v2) : matrices des systèmes associés.
for (P2, P1, P0), Abon, leurres in [
    ((1 + t**2, t, -1), Matrix([[0, 1], [1 / (1 + t**2), -t / (1 + t**2)]]),
     [Matrix([[0, 1], [-1 / (1 + t**2), t / (1 + t**2)]]), Matrix([[0, 1], [1, -t]]), Matrix([[1, 0], [1 / (1 + t**2), -t / (1 + t**2)]])]),
    ((t, -1, t**3), Matrix([[0, 1], [-t**2, 1 / t]]),
     [Matrix([[0, 1], [t**2, -1 / t]]), Matrix([[0, 1], [-t**3, 1]]), Matrix([[0, 1], [-t**2, -1 / t]])]),
    ((1, sin(t), exp(t)), Matrix([[0, 1], [-exp(t), -sin(t)]]),
     [Matrix([[0, 1], [exp(t), sin(t)]]), Matrix([[0, 1], [-sin(t), -exp(t)]]), Matrix([[1, 0], [-exp(t), -sin(t)]])]),
]:
    assert mat_nul(Abon - Matrix([[0, 1], [-P0 / P2, -P1 / P2]]))
    assert all(not mat_nul(L - Abon) for L in leurres)

# ed-l2-systemes-ordre2-ex2 (+ v1, v2) : formule d'Abel.
for a_c, t0, w0, t1, bon, leurres in [
    (2 * t / (1 + t**2), 0, 2, 1, 1, [4, 2, 2 * exp(-1)]),
    (-1 / t, 1, 3, 2, 6, [R_(3, 2), 3, 3 * exp(R_(3, 2))]),
    (R_(2), 0, -1, 1, -exp(-2), [-exp(2), -1, exp(-2)]),
]:
    val = w0 * exp(-integrate(a_c.subs(t, s), (s, t0, t1)))
    assert nul(val - bon) and all(not nul(bon - L) for L in leurres)
# Leurres « non normalisé » : 2 e^{-1} (coefficient 2t), 3 e^{3/2} (coefficient -t).
assert nul(2 * exp(-integrate(2 * s, (s, 0, 1))) - 2 * exp(-1))
assert nul(3 * exp(integrate(s, (s, 1, 2))) - 3 * exp(R_(3, 2)))

# ed-l2-systemes-ordre2-ex3 (+ v1, v2) : seconde solution par abaissement.
for (P2, P1, P0), x1, x2, faux in [
    ((t, -(t + 1), 1), exp(t), t + 1, [t, t * exp(t), exp(-t)]),
    ((2 * t + 1, -4 * (t + 1), 4), exp(2 * t), t + 1, [t, t * exp(2 * t), 2 * t + 1]),
    ((t**2, -t * (t + 2), t + 2), t, t * exp(t), [exp(t), t**2, t * log(t)]),
]:
    op = lambda f: P2 * diff(f, t, 2) + P1 * diff(f, t) + P0 * f  # noqa: E731
    assert nul(op(x1)) and nul(op(x2))
    for f in faux:
        assert not nul(op(f))
# Équation en u = z' : u' = (1/t - 1) u, u = t e^{-t} ; ∫ t e^{-t} = -(t+1) e^{-t}.
assert nul(diff(t * exp(-t), t) - (1 / t - 1) * t * exp(-t))
assert nul(diff(-(t + 1) * exp(-t), t) - t * exp(-t))
assert nul(diff(-(t + 1) * exp(-2 * t), t) - (2 * t + 1) * exp(-2 * t))
assert nul(diff((2 * t + 1) * exp(-2 * t), t) - (-4 * t / (2 * t + 1)) * (2 * t + 1) * exp(-2 * t))
# Wronskien de (e^t, t + 1) : -t e^t, et Abel w' = ((t+1)/t) w.
wq = exp(t) * 1 - exp(t) * (t + 1)
assert nul(wq + t * exp(t)) and nul(diff(wq, t) - (t + 1) / t * wq)

# ed-l2-systemes-ordre2-ex4 (+ v1, v2) : variation des constantes ; propositions.
for (P2, P1, P0), c, base, xp_att, xbon, leurres in [
    ((t, -(t + 1), 1), t**2 * exp(t), (exp(t), t + 1), (t**2 / 2 - 1) * exp(t), t**2 / 2 * exp(t),
     [t**2 * exp(t), (t**2 / 2 - t) * exp(t), t**3 / 3 * exp(t)]),
    ((t, -(t + 1), 1), t**2, (exp(t), t + 1), -t**2 - 2 * t - 2, -t**2, [t**2, -t**2 - 2 * t, -t**3 / 3]),
    ((t**2, -t * (t + 2), t + 2), t**3, (t, t * exp(t)), -t**2 - t, -t**2, [t**2, -t**2 - 1, t**2 * exp(t)]),
]:
    op = lambda f: P2 * diff(f, t, 2) + P1 * diff(f, t) + P0 * f  # noqa: E731
    x1, x2 = base
    cn = c / P2
    wv = simplify(x1 * diff(x2, t) - diff(x1, t) * x2)
    l1 = integrate(simplify(-x2 * cn / wv), t)
    l2 = integrate(simplify(x1 * cn / wv), t)
    assert nul(l1 * x1 + l2 * x2 - xp_att)
    assert nul(op(xp_att) - c) and nul(op(xbon) - c)
    for L in leurres:
        assert not nul(op(L) - c)
assert nul(simplify(exp(t) * diff(t + 1, t) - exp(t) * (t + 1)) + t * exp(t))
assert nul(t * diff(t * exp(t), t) - 1 * t * exp(t) - t**2 * exp(t))  # w = t^2 e^t (v2)

# ed-l2-systemes-ordre2-ex5 : exemple sin / cos.
assert sin(pi / 2) == 1 and cos(pi / 2) == 0

# ---------------------------------------------------------------------------
# QCM ed-l2-systemes
# ---------------------------------------------------------------------------
# q1 : points exclus -2, -1, 2 -> ]-1, 2[ autour de 0.
assert sorted([-2, 2, -1]) == [-2, -1, 2] and -1 < 0 < 2
# q2 : W(1) = e W(0).
assert nul(exp(integrate(1, (s, 0, 1))) - E)
# q3 : e^{tA} pour A = [[0,1],[-1,0]].
assert mat_nul((t * A_rot).exp() - Matrix([[cos(t), sin(t)], [-sin(t), cos(t)]]))
for L in [Rr, Matrix([[cosh(t), sinh(t)], [-sinh(t), cosh(t)]]), Matrix([[1, t], [-t, 1]])]:
    assert not mat_nul((t * A_rot).exp() - L)
# q5 : (t e^t, e^t) seule solution pour [[1,1],[0,1]].
assert est_solution(Matrix([t * exp(t), exp(t)]), A_n)
for L in [Matrix([exp(t), t * exp(t)]), Matrix([t * exp(t), 0]), Matrix([exp(t), exp(t)])]:
    assert not est_solution(L, A_n)
# q6 : valeurs propres -1, -2 ; solution générale.
assert est_solution(Matrix([a_ * exp(-t) - 5 * b_ * exp(-2 * t), b_ * exp(-2 * t)]), Matrix([[-1, 5], [0, -2]]))
# q7 : w = 1/t^2 ; vérifié sur sin t / t et cos t / t (w = -1/t^2, même loi).
assert nul(exp(-integrate(2 / s, (s, 1, t))) - 1 / t**2)
f1, f2 = sin(t) / t, cos(t) / t
for f in (f1, f2):
    assert nul(t * diff(f, t, 2) + 2 * diff(f, t) + t * f)
assert nul(f1 * diff(f2, t) - diff(f1, t) * f2 + 1 / t**2)
# q8 : formule sans R(t0)^{-1} : vaut 2 en 0 pour R = 2e^t, x' = x + 1, x(0) = 1.
assert (2 * exp(t) * (1 + integrate(exp(-s) / 2, (s, 0, t)))).subs(t, 0) == 2

# ---------------------------------------------------------------------------
# Relecture indépendante (2026-10-29) : faits re-dérivés qui n'étaient pas encore contrôlés.
# ---------------------------------------------------------------------------
# ed-l2-systemes-cauchy, trou 3 : le leurre « au moins une solution non nulle » est faux (B = 0, X0 = 0 :
# la solution nulle est la seule, car toute solution s'écrit e^{(t - t0)A} X0 à coefficients constants).
assert mat_nul((Matrix([[0, 1], [-1, 0]]) * t).exp() * zeros(2, 1))
# ed-l2-systemes-dimension, piège 2 : solutions t et t + 1 de x' = 1 (n = 1) ; « wronskien » t, nul en 0 seulement.
assert nul(diff(t, t) - 1) and t.subs(t, 0) == 0 and t.subs(t, 1) != 0
# ed-l2-systemes-ordre2, trou 1 : le leurre « x'(t0) = α et x''(t0) = β » est faux (x'' = 0, β = 1 : x'' est nulle).
assert diff(a_ + b_ * t, t, 2) == 0  # solutions a + b t : x''(0) = 1 est impossible
# ed-l2-systemes-ordre2, contre-exemple : raccord de λ1 t + μ1 t^2 (t > 0) et λ2 t + μ2 t^2 (t < 0) en 0 ;
# x'(0) impose λ1 = λ2, x''(0) impose μ1 = μ2.
l1_, l2_, m1_, m2_ = symbols("l1 l2 m1 m2")
assert diff(l1_ * t + m1_ * t**2, t).subs(t, 0) == l1_ and diff(l2_ * t + m2_ * t**2, t, 2) == 2 * m2_
# ed-l2-systemes-reduction, contre-exemple : parties réelle et imaginaire de e^{it}(1, -i) pour J, wronskien -1.
Re_, Im_ = Matrix([cos(t), sin(t)]), Matrix([sin(t), -cos(t)])
assert est_solution(Re_, J) and est_solution(Im_, J) and nul(Matrix.hstack(Re_, Im_).det() + 1)
# ed-l2-systemes-reduction-ex5 (+ v1) : solutions qui ne tendent pas vers 0 aux valeurs limites a = 0 et a = 1.
assert est_solution(Matrix([cos(t), -sin(t)]), Matrix([[0, 1], [-1, 0]]))
assert Matrix([[-1, 1], [1, -1]]) * Matrix([1, 1]) == zeros(2, 1)
# q7 : leurres t^2, 1/t, e^{1 - t^2} ne vérifient pas w' = -(2/t) w.
for L in [t**2, 1 / t, exp(1 - t**2)]:
    assert not nul(diff(L, t) + 2 / t * L)
# ed-l2-systemes-cauchy-ex4 : base de l'espace des solutions sur R pour k = 2 : max(t,0)^2, min(t,0)^2,
# dérivables en 0 (dérivées à gauche et à droite nulles) et solutions de t x' = 2x de chaque côté.
for f in [t**2, 0 * t]:
    assert nul(t * diff(f, t) - 2 * f) and diff(f, t).subs(t, 0) == 0

print("ed_l2_systemes : toutes les vérifications passent.")
