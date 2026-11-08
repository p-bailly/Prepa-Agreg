"""Vérifications calculatoires du chapitre « proba-l3-limites » (probabilités, L3, théorèmes limites) : cartes,
exercices et QCM de content/*/probabilites/l3-proba-l3-limites.json et content/qcm/proba-l3-limites.json.

Chaque assert cite l'identifiant de la carte, de l'exercice (variante comprise : v1, v2) ou de la question. Les calculs
sont faits en fractions exactes, avec SymPy (sommes, produits, séries), par énumération (petits univers finis) ou avec
mpmath à 30 chiffres (lois binomiales et de Poisson exactes, fonction de répartition Φ de la loi normale).
Lancer : python3 verification/proba_l3_limites.py.
"""

from fractions import Fraction as F
from itertools import product
from math import ceil, comb, floor, log

import mpmath as mp
import sympy as sp

mp.mp.dps = 30
Phi = mp.ncdf
k, n, N, x = sp.symbols("k n N x", positive=True, integer=False)
ki, ni, Ni = sp.symbols("k_i n_i N_i", positive=True, integer=True)


def proche(a, b, tol=mp.mpf("1e-12")):
    return abs(mp.mpf(a) - mp.mpf(b)) < tol


def binom_cdf(nn, p, m):
    """P(X <= m) pour X de loi B(nn, p), exacte (mpmath)."""
    p = mp.mpf(p)
    return mp.fsum(mp.binomial(nn, j) * p**j * (1 - p) ** (nn - j) for j in range(0, m + 1))


def binom_pmf(nn, p, j):
    p = mp.mpf(p)
    return mp.binomial(nn, j) * p**j * (1 - p) ** (nn - j)


# =====================================================================================================
# Carte proba-l3-limites-convergences
# =====================================================================================================
# trou 1 (leurres « uniformément », « en tout point », « en moyenne ») : lancers d'une pièce, D_n = « n premiers
# lancers faces », X_n = 2^n 1_{D_n}. P(D_n) = 2^-n -> 0 (convergence p.s. vers 0, car hors de « jamais pile », de
# probabilité nulle, X_n finit par valoir 0), E(X_n) = 1, et sur l'issue « que des faces », X_n = 2^n -> +inf.
for m in range(1, 40):
    assert F(2**m) * F(1, 2**m) == 1  # E(X_n) = 1 : pas de convergence en moyenne
    assert 2**m >= 2  # sup des X_n et valeur sur l'issue « que des faces » : ne tend pas vers 0

# trou 2 (leurres « en tout réel t », « uniformément », « où F_X est discontinue ») : X_n = 1/n, X = 0.
def F_const(c, t):
    return 1 if t >= c else 0


for m in range(1, 50):
    assert F_const(F(1, m), 0) == 0 and F_const(0, 0) == 1  # F_{X_n}(0) = 0 ne tend pas vers F_X(0) = 1
    for t in [F(1, 7), F(1, 3), 1, 5]:
        if m >= 1 / t:
            assert F_const(F(1, m), t) == F_const(0, t)  # convergence aux points de continuité t > 0
    assert F_const(F(1, m), F(-1, 5)) == F_const(0, F(-1, 5)) == 0

# trou 3 (leurres « bornée », « deux valeurs », « admet une espérance ») et contre-exemple « en loi sans être en
# probabilité » : X = ±1, X_n = -X. Même loi, |X_n - X| = 2.
loi_X = {1: F(1, 2), -1: F(1, 2)}
loi_Xn = {}
for v, p in loi_X.items():
    loi_Xn[-v] = loi_Xn.get(-v, 0) + p
assert loi_Xn == loi_X
assert all(abs(-v - v) == 2 for v in loi_X)

# exemple : X_n uniforme sur {1/n, ..., n/n} : F_{X_n}(t) = floor(nt)/n sur [0,1] et |F_{X_n}(t) - t| <= 1/n.
for m in [1, 2, 3, 7, 10, 25]:
    for t in [F(0), F(1, 10), F(1, 3), F(1, 2), F(5, 7), F(1)]:
        Fn = F(sum(1 for j in range(1, m + 1) if F(j, m) <= t), m)
        assert Fn == F(floor(m * t), m)
        assert abs(Fn - t) <= F(1, m)

# contre-exemple « en probabilité sans être p.s. » : prod_{k=n}^N (1 - 1/k) = (n-1)/N.
assert sp.simplify(sp.product(1 - 1 / sp.Symbol("j"), (sp.Symbol("j"), ni, Ni)) - (ni - 1) / Ni) == 0
for a in range(2, 8):
    for b in range(a, 30):
        pr = F(1)
        for j in range(a, b + 1):
            pr *= 1 - F(1, j)
        assert pr == F(a - 1, b)

# contre-exemple « un seul ε » : X_n = 1/2, X = 0 : P(|X_n - X| >= 1) = 0 mais F_{X_n}(1/4) = 0 != 1 = F_X(1/4).
assert F_const(F(1, 2), F(1, 4)) == 0 and F_const(0, F(1, 4)) == 1

# démonstration (limite constante) : {|y - c| >= eps} ⊂ {y <= c - eps} ∪ {y > c + eps/2}, testé sur une grille.
for c in [F(0), F(3, 2)]:
    for eps in [F(1, 10), F(1, 2), F(2)]:
        for j in range(-80, 81):
            y = c + F(j, 20)
            if abs(y - c) >= eps:
                assert y <= c - eps or y > c + eps / 2

# affirmation 3 : X_n = n avec probabilité 1/n : E(X_n) = 1.
assert all(F(m) * F(1, m) == 1 for m in range(1, 100))

# =====================================================================================================
# Carte proba-l3-limites-loi-faible
# =====================================================================================================
def loi_moyenne(loi, nn):
    """Loi exacte de la moyenne de nn variables indépendantes de loi « loi » (dict valeur -> probabilité)."""
    res = {F(0): F(1)}
    for _ in range(nn):
        nouv = {}
        for s, ps in res.items():
            for v, pv in loi.items():
                nouv[s + v] = nouv.get(s + v, 0) + ps * pv
        res = nouv
    return {s / nn: p for s, p in res.items()}


def proba_ecart(loi_m, mu, eps):
    return sum(p for v, p in loi_m.items() if abs(v - mu) >= eps)


# trou 2, leurre σ²/(nε) : X = ±1/2, n = 1, ε = 1/2 : P = 1 > (1/4)/(1/2).
L = {F(1, 2): F(1, 2), F(-1, 2): F(1, 2)}
assert proba_ecart(loi_moyenne(L, 1), 0, F(1, 2)) == 1 > F(1, 4) / F(1, 2)
# leurre σ/(nε²) : X = ±2, n = 1, ε = 2 : P = 1 > 2/4.
L = {F(2): F(1, 2), F(-2): F(1, 2)}
assert proba_ecart(loi_moyenne(L, 1), 0, 2) == 1 > F(2, 4)
# leurre σ²/(n²ε²) : X = ±1, n = 2, ε = 1 : P(|M_2| >= 1) = 1/2 > 1/4.
L = {F(1): F(1, 2), F(-1): F(1, 2)}
assert proba_ecart(loi_moyenne(L, 2), 0, 1) == F(1, 2) > F(1, 4)
# la vraie borne σ²/(nε²) est respectée sur ces exemples (et d'autres).
for (L, eps) in [({F(1, 2): F(1, 2), F(-1, 2): F(1, 2)}, F(1, 2)), ({F(2): F(1, 2), F(-2): F(1, 2)}, F(2)),
                 ({F(1): F(1, 2), F(-1): F(1, 2)}, F(1)), ({F(0): F(2, 3), F(3): F(1, 3)}, F(1, 2))]:
    mu = sum(v * p for v, p in L.items())
    s2 = sum((v - mu) ** 2 * p for v, p in L.items())
    for m in range(1, 7):
        assert proba_ecart(loi_moyenne(L, m), mu, eps) <= s2 / (m * eps**2)

# trou 1, leurre « toutes égales à X_1 » : X = ±1, n = 4, ε = 1 : M_n = X et P = 1 > σ²/(nε²) = 1/4.
assert 1 > F(1, 4)
# leurres « covariances > 0 » et « covariances = σ²/2 » : X_k = Z + Y_k, Z, Y_k = ±1 indépendantes.
loiXk = {}
for z, y in product([-1, 1], repeat=2):
    loiXk[z + y] = loiXk.get(z + y, 0) + F(1, 4)
assert loiXk == {-2: F(1, 4), 0: F(1, 2), 2: F(1, 4)}  # même loi pour tous les k
var_Xk = sum(v**2 * p for v, p in loiXk.items())
assert var_Xk == 2
cov = sum(F(1, 8) * (z + y1) * (z + y2) for z, y1, y2 in product([-1, 1], repeat=3))
assert cov == 1 == var_Xk / 2  # Cov(X_j, X_k) = σ²/2 > 0
# n = 100, ε = 1/2 : P(|M_n| >= 1/2) = P(|Z + Ybar| >= 1/2) calculée exactement, contre la borne 2/(100/4) = 2/25.
nn = 100
pY = {}
for j in range(nn + 1):  # j valeurs +1 parmi les Y_k : somme 2j - nn
    pY[F(2 * j - nn, nn)] = F(comb(nn, j), 2**nn)
P_ecart = sum(F(1, 2) * pv for z in [-1, 1] for v, pv in pY.items() if abs(z + v) >= F(1, 2))
assert P_ecart > F(99, 100) and P_ecart > F(2, 25)

# exemple : p(1-p)/(n ε²) <= 1/(4 n ε²) = 100/n pour ε = 1/20, <= 1/20 dès n >= 2000.
assert F(1, 4) / (F(1, 20) ** 2) == 100 and F(100, 2000) == F(1, 20)

# contre-exemple « même espérance, lois différentes » : |±3 ± 9 ± ... ± 3^n| >= (3^n + 3)/2 >= n.
for m in range(1, 11):
    mini = min(abs(sum(s * 3**j for s, j in zip(signes, range(1, m + 1)))) for signes in product([-1, 1], repeat=m))
    assert mini == (3**m + 3) // 2 and F(3**m + 3, 2 * m) >= 1

# affirmation 3 : Bernoulli 1/2, n impair : P(M_n = 1/2) = 0 ; n pair : P(M_n = 1/2) -> 0.
assert all(F(1, 2) not in loi_moyenne({F(0): F(1, 2), F(1): F(1, 2)}, m) for m in [1, 3, 5, 7])
assert F(comb(200, 100), 2**200) < F(1, 10)

# =====================================================================================================
# Carte proba-l3-limites-borel-cantelli
# =====================================================================================================
# démonstration : 1 - x <= e^{-x}.
assert all(1 - mp.mpf(j) / 100 <= mp.e ** (-mp.mpf(j) / 100) for j in range(0, 101))
# trou 1, leurres : P(A_n) = 1/(n+1) décroissante, -> 0, < 1, mais série divergente (second lemme : probabilité 1).
assert all(F(1, m + 1) > F(1, m + 2) and F(1, m + 1) < 1 for m in range(1, 1000))
assert sp.Sum(1 / (ki + 1), (ki, 1, sp.oo)).is_convergent() is sp.false
# trou 2, leurres : P(A_1) = 1, P(A_n) = 2^-n (n >= 2) : > 0, décroissante, somme >= 1, mais série convergente.
seq = [F(1)] + [F(1, 2**m) for m in range(2, 60)]
assert all(a > 0 for a in seq) and all(a > b for a, b in zip(seq, seq[1:])) and sum(seq) >= 1
assert sp.Sum(2 ** (-ki), (ki, 2, sp.oo)).doit() == sp.Rational(1, 2)
# exemple : P(A_n) = 2^-n (série convergente), P(B_n) = 1/4 (série divergente).
assert sp.Sum(2 ** (-ki), (ki, 1, sp.oo)).doit() == 1
# contre-exemple 3 : P(A_n) = 2^-n : série convergente (ci-dessus).
# trou 3, leurre « tous réalisés à partir d'un certain rang » : P(A_n) = 1/2 indépendants : P(A_n ∩ ... ∩ A_N) = 2^-(N-n+1) -> 0.
assert all(F(1, 2) ** (b - a + 1) <= F(1, 2**(b - a + 1)) for a in range(1, 5) for b in range(a, 30))

# =====================================================================================================
# Carte proba-l3-limites-loi-forte
# =====================================================================================================
def E_S4_brut(loi, m):
    """E(S_m^4) par énumération exacte pour m variables indépendantes de loi « loi »."""
    tot = F(0)
    for vals in product(list(loi.items()), repeat=m):
        p = F(1)
        s = 0
        for v, pv in vals:
            p *= pv
            s += v
        tot += p * s**4
    return tot


# démonstration : E(S_n^4) = n m4 + 3 n (n-1) σ⁴ pour des variables centrées (et 6 quadruplets par paire {i,j}).
assert len([q for q in product([0, 1], repeat=4) if sum(q) == 2]) == 6
for loi in [{1: F(1, 2), -1: F(1, 2)}, {-1: F(1, 3), 0: F(1, 3), 1: F(1, 3)}, {2: F(1, 2), -2: F(1, 2)},
            {-2: F(1, 3), 1: F(2, 3)}]:
    assert sum(v * p for v, p in loi.items()) == 0
    m2 = sum(v**2 * p for v, p in loi.items())
    m4 = sum(v**4 * p for v, p in loi.items())
    assert m2**2 <= m4
    for m in range(1, 6):
        assert E_S4_brut(loi, m) == m * m4 + 3 * m * (m - 1) * m2**2 <= m * m * (m4 + 3 * m2**2)
# |x|^r <= 1 + x^4 pour 0 <= r <= 4.
assert all(abs(mp.mpf(j) / 10) ** (mp.mpf(r) / 4) <= 1 + (mp.mpf(j) / 10) ** 4 for j in range(-60, 61) for r in range(0, 17))

# contre-exemple « sans espérance » : P(X = ±k) = 1/(2k(k+1)) : loi, P(|X| >= n) = 1/n, E|X| infinie.
assert sp.Sum(1 / (ki * (ki + 1)), (ki, 1, sp.oo)).doit() == 1
assert sp.simplify(sp.Sum(1 / (ki * (ki + 1)), (ki, ni, sp.oo)).doit() - 1 / ni) == 0
assert sp.Sum(1 / ni, (ni, 1, sp.oo)).is_convergent() is sp.false
# contre-exemple « même espérance, lois différentes » : |M_n| >= (3^n+3)/(2n) -> +inf (vérifié plus haut).
# affirmation 2 : P(X >= k) = 1/k² : E(X) = π²/6, E(X²) = Σ (2k-1)/k² = +inf ; P(X = k) >= 0.
assert sp.Sum(1 / ki**2, (ki, 1, sp.oo)).doit() == sp.pi**2 / 6
assert sp.Sum((2 * ki - 1) / ki**2, (ki, 1, sp.oo)).is_convergent() is sp.false
assert all(F(1, j * j) - F(1, (j + 1) ** 2) > 0 for j in range(1, 100))
# E(Y²) = Σ (2k-1) P(Y >= k) pour Y à valeurs dans N* : contrôle sur une loi finie.
loiY = {1: F(1, 2), 2: F(1, 3), 5: F(1, 6)}
assert sum(v * v * p for v, p in loiY.items()) == sum((2 * j - 1) * sum(p for v, p in loiY.items() if v >= j) for j in range(1, 6))
# exemple : X de loi P(λ) : E(X²) = λ + λ².
lam = sp.Symbol("lam", positive=True)
assert sp.simplify(sp.Sum(ki**2 * sp.exp(-lam) * lam**ki / sp.factorial(ki), (ki, 0, sp.oo)).doit() - (lam + lam**2)) == 0
# trou 2, leurres : pièce équilibrée : P(M_n = 1/2) -> 0 (n = 1000 : < 0,03) ; sup |M_n - 1/2| = 1/2 (que des piles).
assert binom_pmf(1000, mp.mpf(1) / 2, 500) < mp.mpf("0.03")

# =====================================================================================================
# Carte proba-l3-limites-tcl
# =====================================================================================================
assert abs(Phi(mp.mpf("1.96")) - mp.mpf("0.975")) < mp.mpf("0.0001")  # Φ(1,96) ≈ 0,975
for t in [0, 0.5, 1, 2, 3.7]:
    assert proche(Phi(-t), 1 - Phi(t))  # Φ(-t) = 1 - Φ(t)
assert abs(2 * Phi(1) - 1 - mp.mpf("0.683")) < mp.mpf("0.0005")  # exemple
assert abs(1.96 * 30 - 59) < 0.5  # exemple : 1,96 × 30 ≈ 59
# exemple : valeurs exactes pour S de loi B(3600, 1/2).
P1 = binom_cdf(3600, mp.mpf(1) / 2, 1830) - binom_cdf(3600, mp.mpf(1) / 2, 1769)
assert abs(P1 - mp.mpf("0.683")) < mp.mpf("0.01")
P2 = binom_cdf(3600, mp.mpf(1) / 2, 1859) - binom_cdf(3600, mp.mpf(1) / 2, 1740)
assert abs(P2 - mp.mpf("0.95")) < mp.mpf("0.01")
# contre-exemple 1 : X_k = ±k avec probabilité 1/(2k²) : E = 0, V = 1, Σ 1/k² < ∞ ; Φ(-1) > 0.
for j in range(1, 40):
    assert 2 * F(1, 2 * j * j) <= 1 and F(j * j) * 2 * F(1, 2 * j * j) == 1
assert Phi(-1) > mp.mpf("0.15")
# contre-exemple 2 : X_k = ε_k Y, Y ∈ {0, √2} : même loi (0 : 1/2, ±√2 : 1/4), V = 1, covariance nulle.
r2 = sp.sqrt(2)
loiX = {}
for e, y in product([-1, 1], [0, r2]):
    loiX[e * y] = loiX.get(e * y, 0) + sp.Rational(1, 4)
assert loiX == {0: sp.Rational(1, 2), r2: sp.Rational(1, 4), -r2: sp.Rational(1, 4)}
assert sum(v**2 * p for v, p in loiX.items()) == 1
assert sum(sp.Rational(1, 8) * e1 * e2 * y**2 for e1, e2, y in product([-1, 1], [-1, 1], [0, r2])) == 0
assert Phi(0.5) - Phi(-0.5) < mp.mpf("0.3830") and Phi(0.5) - Phi(-0.5) > mp.mpf("0.3829")  # ≈ 0,383 < 1/2
# contre-exemple 3 : Φ(1) ≈ 0,841 != 1/2.
assert abs(Phi(1) - mp.mpf("0.841")) < mp.mpf("0.0005")
# trou 1, leurre σ²√n (X = ±2) : P(-1 <= Z_n/2 <= 1) -> Φ(2) - Φ(-2) != Φ(1) - Φ(-1).
assert abs((Phi(2) - Phi(-2)) - (Phi(1) - Phi(-1))) > mp.mpf("0.2")
# trou 2, leurres : a = -1, b = 1 : Φ(b - a) = Φ(2) ; a = 0, b = 1 : 2Φ(1) - 1 ; intégrale de e^{-x²}.
assert abs(Phi(2) - (Phi(1) - Phi(-1))) > mp.mpf("0.2")
assert abs((2 * Phi(1) - 1) - (Phi(1) - Phi(0))) > mp.mpf("0.2")
I = mp.quad(lambda u: mp.e ** (-u * u), [-1, 1]) / mp.sqrt(2 * mp.pi)
assert abs(I - mp.mpf("0.596")) < mp.mpf("0.001") and abs(I - (Phi(1) - Phi(-1))) > mp.mpf("0.05")
# affirmation 1 : Φ(0) = 1/2.
assert proche(Phi(0), mp.mpf(1) / 2)

# =====================================================================================================
# Carte proba-l3-limites-binomiale
# =====================================================================================================
# démonstration : -x - x² <= ln(1 - x) <= -x sur [0, 1/2].
for j in range(0, 501):
    u = mp.mpf(j) / 1000
    assert -u - u * u <= mp.log(1 - u) <= -u
# convergence des probabilités ponctuelles (contrôle numérique, n = 10^6) : P(B(n, λ/n) = k) ≈ e^{-λ} λ^k / k!.
for lam_ in [1, 2, mp.mpf("3.5")]:
    for kk in range(0, 6):
        q = mp.e ** (-lam_) * mp.mpf(lam_) ** kk / mp.factorial(kk)
        assert abs(binom_pmf(10**6, mp.mpf(lam_) / 10**6, kk) - q) < mp.mpf("1e-5")
# exemple : 0,998^500 ≈ 0,3675 ; e^{-1} ≈ 0,3679 ; B(100, 1/2) : P(45 <= X <= 55) ≈ 0,729 = 2Φ(1,1) - 1 ≈ 0,729.
assert abs(mp.mpf("0.998") ** 500 - mp.mpf("0.3675")) < mp.mpf("0.00005")
assert abs(mp.e**-1 - mp.mpf("0.3679")) < mp.mpf("0.00005")
Pex = binom_cdf(100, mp.mpf(1) / 2, 55) - binom_cdf(100, mp.mpf(1) / 2, 44)
assert abs(Pex - mp.mpf("0.729")) < mp.mpf("0.0005") and abs(2 * Phi(1.1) - 1 - mp.mpf("0.729")) < mp.mpf("0.0005")
# contre-exemples : p = 1/2 : 2^-n -> 0 ; p_n = 1/n² : (1 - 1/n²)^n -> 1 ; p fixé : n^k (1-p)^(n-k) -> 0.
assert (1 - mp.mpf(1) / 10**12) ** (10**6) > mp.mpf("0.999999")
assert mp.mpf(10**4) ** 3 * mp.mpf("0.9") ** (10**4 - 3) < mp.mpf("1e-400")
# trou 1, leurres : p_n = 1/2 : P(X_n = 0) -> 0 != e^{-1/2} ; n² p_n -> 1 : -> 1 != e^{-1} ; √n p_n -> 1 : -> 0.
assert binom_pmf(10**4, mp.mpf(1) / 10**8, 0) > mp.mpf("0.9998")
assert binom_pmf(10**6, mp.mpf(1) / 10**3, 0) < mp.mpf("1e-400")
# trou 2, leurres : k = 0 : λ^0/0! = 1 != e^{-λ} ; k = 2 : e^{-λ}λ² != e^{-λ}λ²/2 ; k = 1, λ = 2 : e^{-1}·2 != e^{-2}·2.
assert mp.e**-2 * 2 != mp.e**-1 * 2 and 1 != mp.e**-2
# trou 3, leurres : n = 10^4, p = 1/2 : P(|X - np|/(np(1-p)) >= 1/10) = P(|X - 5000| >= 250) minuscule (limite 0) ;
# P(|X - np|/√(p(1-p)) <= 1) = P(X = 5000) ≈ 0,008 (limite 0, et non Φ(1) - Φ(-1)).
Pa = 1 - (binom_cdf(10**4, mp.mpf(1) / 2, 5249) - binom_cdf(10**4, mp.mpf(1) / 2, 4750))
assert Pa < mp.mpf("1e-6")
assert binom_pmf(10**4, mp.mpf(1) / 2, 5000) < mp.mpf("0.01")
# affirmation 3 : p_n = 1/√n : P(X_n = k) -> 0 (n = 10^6).
for kk in range(0, 6):
    assert binom_pmf(10**6, mp.mpf(1) / 1000, kk) < mp.mpf("1e-400")
# affirmation 1 : (1 - 3/n)^n -> e^{-3}.
assert abs((1 - mp.mpf(3) / 10**8) ** (10**8) - mp.e**-3) < mp.mpf("1e-8")

# =====================================================================================================
# Exercices de la carte convergences
# =====================================================================================================
# convergences-ex1 (et v1, v2) : E(X_n) = v/n pour v = n², n, √n ; P(X_n != 0) = 1/n.
nn_ = sp.Symbol("n", positive=True)
assert sp.simplify(nn_**2 / nn_) == nn_  # ex1 : E = n -> +inf
assert sp.simplify(nn_ / nn_) == 1  # v1 : E = 1
assert sp.simplify(sp.sqrt(nn_) / nn_ - 1 / sp.sqrt(nn_)) == 0 and sp.limit(1 / sp.sqrt(nn_), nn_, sp.oo) == 0  # v2

# convergences-ex2 : X Bernoulli 1/2, X_n = 1 - X : |X_n - X| = 1 ; P(|X_n - X| >= 1/2) = 1 ; P(X = Y) = 0.
assert all(abs((1 - v) - v) == 1 for v in [0, 1])
assert sum(F(1, 2) for v in [0, 1] if v == 1 - v) == 0
# v1 : X = ±1, X_n = -X : |X_n - X| = 2 ; P = 1 ; P(X = -X) = 0.
assert all(abs(-v - v) == 2 for v in [-1, 1])
# v2 : X uniforme sur {1,2,3}, X_n = 4 - X : même loi ; P(|X_n - X| >= 1) = 2/3 ; P(X = Y) = 1/3.
assert sorted(4 - v for v in [1, 2, 3]) == [1, 2, 3]
assert sum(F(1, 3) for v in [1, 2, 3] if abs(4 - 2 * v) >= 1) == F(2, 3)
assert sum(F(1, 3) for v in [1, 2, 3] if v == 4 - v) == F(1, 3)

# convergences-ex3 (et v1, v2) : P(T_n > m) = (1 - p)^m ; limite de (1 - λ/n)^{floor(n a)} = e^{-λ a}.
p_ = sp.Symbol("p", positive=True)
assert sp.simplify(sp.Sum((1 - p_) ** (ki - 1) * p_, (ki, Ni + 1, sp.oo)).doit().subs(p_, sp.Rational(1, 3)).subs(Ni, 5)
                   - sp.Rational(2, 3) ** 5) == 0
for lam_, a, attendu, leurres in [(1, mp.mpf(2), mp.e**-2, [mp.e ** (-0.5), 1 - mp.e**-2, 0]),
                                  (2, mp.mpf(3), mp.e**-6, [mp.e ** (-1.5), 1 - mp.e**-6, mp.e ** (-mp.mpf(2) / 3)]),
                                  (3, mp.mpf(1) / 3, mp.e**-1, [mp.e**-9, 1 - mp.e**-1, mp.e ** (-mp.mpf(1) / 9)])]:
    m = 10**7
    val = (1 - mp.mpf(lam_) / m) ** int(mp.floor(m * a))
    assert abs(val - attendu) < mp.mpf("1e-5")
    assert all(abs(attendu - l) > mp.mpf("0.01") for l in leurres)

# convergences-ex4 : produits télescopiques et probabilités de la question 3.
j_ = sp.Symbol("j", positive=True, integer=True)
assert sp.simplify(sp.product(1 - 1 / j_**2, (j_, ni, Ni)) - (ni - 1) * (Ni + 1) / (ni * Ni)) == 0
assert sp.limit((2 - 1) * (N + 1) / (2 * N), N, sp.oo) == sp.Rational(1, 2)  # ex4 : P(G_2) = 1/2
assert sp.limit((3 - 1) * (N + 1) / (3 * N), N, sp.oo) == sp.Rational(2, 3)  # v1 : P(G_3) = 2/3
assert sp.simplify(sp.product(1 - 4 / j_**2, (j_, ni, Ni))
                   - (ni - 2) * (ni - 1) * (Ni + 1) * (Ni + 2) / ((Ni - 1) * Ni * ni * (ni + 1))) == 0
assert sp.limit(1 * 2 * (N + 1) * (N + 2) / ((N - 1) * N * 3 * 4), N, sp.oo) == sp.Rational(1, 6)  # v2
assert F(4, 9) <= 1  # v2 : p_3 = 4/9 est une probabilité
assert float(6 / sp.pi**2) != 0.5  # ex4 : leurre 6/π² ≈ 0,608 != 1/2
for nn in range(3, 12):
    pr = F(1)
    for j in range(nn, 400):
        pr *= 1 - F(4, j * j)
    assert pr == F((nn - 2) * (nn - 1) * 400 * 401, 398 * 399 * nn * (nn + 1))  # v2 : formule exacte à N = 399

# =====================================================================================================
# Exercices de la carte loi-faible
# =====================================================================================================
# loi-faible-ex1 (et v1, v2) : majorant σ²/(nε²) et leurres (ε au lieu de ε², n², variance 1/4).
s2_de = F(91, 6) - F(49, 4)
assert s2_de == F(35, 12)
for (s2, nn, eps, bonne, leurres) in [(s2_de, 1000, F(1, 10), F(7, 24), [F(7, 240), F(7, 24000), F(1, 40)]),
                                      (s2_de, 3500, F(1, 10), F(1, 12), [F(1, 120), F(1, 42000), F(1, 140)]),
                                      (F(3, 10) * F(7, 10), 2100, F(1, 20), F(1, 25), [F(1, 500), F(1, 52500), F(1, 21)])]:
    assert s2 / (nn * eps**2) == bonne
    assert [s2 / (nn * eps), s2 / (nn * nn * eps**2), F(1, 4) / (nn * eps**2)] == leurres
    assert bonne not in leurres

# loi-faible-ex2 (et v1, v2) : X_k = ε_k Y, loi, covariance nulle, dépendance, plus petit n.
for (y1, y2, eps, alpha, n_att, leurres) in [(1, 3, F(1, 10), F(1, 10), 5000, [500, 4000, 50]),
                                              (1, 2, F(1, 10), F(1, 10), 2500, [250, 2250, 25]),
                                              (0, 2, F(1, 10), F(1, 20), 4000, [400, 2000, 10])]:
    loi = {}
    for e, y in product([-1, 1], [y1, y2]):
        loi[e * y] = loi.get(e * y, 0) + F(1, 4)
    EY2 = F(y1 * y1 + y2 * y2, 2)
    EY = F(y1 + y2, 2)
    assert sum(v * p for v, p in loi.items()) == 0 and sum(v * v * p for v, p in loi.items()) == EY2
    assert sum(F(1, 8) * e1 * e2 * y * y for e1, e2, y in product([-1, 1], [-1, 1], [y1, y2])) == 0  # Cov = 0
    # dépendance : P(X_1 = a, X_2 = b) = 0 avec |a| != |b|, alors que le produit est > 0.
    a, b = (1, 3) if y1 == 1 and y2 == 3 else ((1, 2) if y1 == 1 else (0, 2))
    joint = sum(F(1, 8) for e1, e2, y in product([-1, 1], [-1, 1], [y1, y2]) if e1 * y == a and e2 * y == b)
    assert joint == 0 and loi[a] * loi[b] > 0
    plus_petit = lambda var, e2, al: next(m for m in range(1, 10**7) if var / (m * e2) <= al)
    assert plus_petit(EY2, eps**2, alpha) == n_att
    assert [plus_petit(EY2, eps, alpha), plus_petit(EY**2, eps**2, alpha), EY2 / eps**2 * alpha] == leurres

# loi-faible-ex3 (et v1, v2) : moyenne des paramètres et sa limite.
def pbar(pk, m):
    return sum(pk(j) for j in range(1, m + 1)) / m


for (pk, ell, borne, leurres) in [(lambda j: F(1, 3) if j % 2 else F(2, 3), F(1, 2), lambda m: F(1, 6 * m), [F(1, 3), F(2, 3)]),
                                  (lambda j: F(9, 10) if j % 4 == 0 else F(1, 10), F(3, 10), lambda m: F(3, 5 * m), [F(1, 2), F(1, 10), F(9, 10)]),
                                  (lambda j: F(1, 2) if j % 3 == 0 else F(1, 4), F(1, 3), lambda m: F(1, 6 * m), [F(3, 8), F(1, 4), F(1, 2)])]:
    for m in range(1, 300):
        assert abs(pbar(pk, m) - ell) <= borne(m)
    assert ell not in leurres
    assert all(pk(j) * (1 - pk(j)) <= F(1, 4) for j in range(1, 50))
# formules exactes des corrections : n = 2m+1 (ex3) ; n = 4m + r (v1) ; n = 3m + r (v2).
assert all(pbar(lambda j: F(1, 3) if j % 2 else F(2, 3), 2 * m + 1) == F(1, 2) - F(1, 6 * (2 * m + 1)) for m in range(0, 50))
assert all(pbar(lambda j: F(9, 10) if j % 4 == 0 else F(1, 10), 4 * m + r) - F(3, 10) == -F(r, 5 * (4 * m + r))
           for m in range(0, 30) for r in range(0, 4) if 4 * m + r > 0)
assert all(pbar(lambda j: F(1, 2) if j % 3 == 0 else F(1, 4), 3 * m + r) - F(1, 3) == -F(r, 12 * (3 * m + r))
           for m in range(0, 30) for r in range(0, 3) if 3 * m + r > 0)

# loi-faible-ex4 (et v1, v2) : e_n et sa comparaison à ln n ; V(R_n) <= e_n.
for m in range(1, 60):
    H = sum(F(1, j) for j in range(1, m + 1))
    H1 = H + F(1, m + 1)
    assert sum(F(2, j + 1) for j in range(1, m + 1)) == 2 * (H1 - 1)  # v1
    assert sum(F(1, 2 * j) for j in range(1, m + 1)) == H / 2  # v2
    for pk in [lambda j: F(1, j), lambda j: F(2, j + 1), lambda j: F(1, 2 * j)]:
        assert sum(pk(j) * (1 - pk(j)) for j in range(1, m + 1)) <= sum(pk(j) for j in range(1, m + 1))
gam = mp.euler
for (en, lim) in [(lambda m: mp.harmonic(m), 1), (lambda m: 2 * (mp.harmonic(m + 1) - 1), 2), (lambda m: mp.harmonic(m) / 2, mp.mpf(1) / 2)]:
    assert abs(en(10**30) / mp.log(10**30) - lim) < mp.mpf("0.02")
assert len({F(1), F(0), F(1, 2)}) == 3 and float(mp.e) not in (1, 0, 0.5)  # ex4 : leurres distincts de 1

# loi-faible-ex5 : E(T_n) = n H_n ; V(T_n) = Σ (1-p)/p² <= n² Σ 1/k² <= π²n²/6.
for m in range(1, 40):
    ps = [F(m - i + 1, m) for i in range(1, m + 1)]
    assert sum(1 / p for p in ps) == m * sum(F(1, j) for j in range(1, m + 1))
    V = sum((1 - p) / p**2 for p in ps)
    assert V <= m * m * sum(F(1, j * j) for j in range(1, m + 1)) <= m * m * mp.pi**2 / 6

# =====================================================================================================
# Exercices de la carte borel-cantelli
# =====================================================================================================
# borel-cantelli-ex1 (et v1, v2) : nature des séries, et p_n <= 1.
assert sp.Sum(1 / ki**2, (ki, 1, sp.oo)).is_convergent() is sp.true
for t in [1 / ki, 1 / sp.sqrt(ki), 1 / sp.log(ki + 2)]:
    assert sp.Sum(t, (ki, 1, sp.oo)).is_convergent() is sp.false
assert sp.Sum(1 / ki**sp.Rational(3, 2), (ki, 1, sp.oo)).is_convergent() is sp.true  # v1
for t in [1 / (ki + 1), 1 / (2 * ki)]:
    assert sp.Sum(t, (ki, 1, sp.oo)).is_convergent() is sp.false
assert sp.limit(sp.log(sp.log(x)) - sp.log(sp.log(2)), x, sp.oo) == sp.oo  # v1 : Σ 1/((n+1)ln(n+1)) diverge
assert sp.diff(sp.log(sp.log(x)), x).equals(1 / (x * sp.log(x)))
assert sp.Sum(ki / 2**ki, (ki, 1, sp.oo)).doit() == 2  # v2 : convergente
for t in [1 / ki ** sp.Rational(1, 3), ki / (ki**2 + 1)]:
    assert sp.Sum(t, (ki, 1, sp.oo)).is_convergent() is sp.false
for m in range(1, 5000):
    assert 1 / log(m + 2) <= 1 and 1 / ((m + 1) * log(m + 1)) <= 1 and m / 2**m <= 0.5 and m / (m * m + 1) <= 0.5
    assert 1 / (2 + mp.sin(m)) >= mp.mpf(1) / 3 and 1 / (2 + mp.sin(m)) <= 1
    assert 1 / log(m + 2) >= 1 / (m + 2)

# borel-cantelli-ex2 (et v1, v2) : sommes et sommes partielles.
mm = sp.Symbol("m", positive=True, integer=True)
assert sp.simplify(sp.Sum(1 / (ki * (ki + 1)), (ki, 1, mm)).doit() - (1 - 1 / (mm + 1))) == 0
assert sp.simplify(sp.Sum(sp.Rational(1, 4) ** ki, (ki, 1, mm)).doit() - sp.Rational(1, 3) * (1 - sp.Rational(1, 4) ** mm)) == 0
assert sp.simplify(sp.Sum(2 / (ki * (ki + 1) * (ki + 2)), (ki, 1, mm)).doit() - (sp.Rational(1, 2) - 1 / ((mm + 1) * (mm + 2)))) == 0
assert sp.simplify(2 / (ki * (ki + 1) * (ki + 2)) - (1 / (ki * (ki + 1)) - 1 / ((ki + 1) * (ki + 2)))) == 0
for bonne, leurres in [(1, [sp.Rational(1, 2), 2, sp.pi**2 / 6]), (sp.Rational(1, 3), [sp.Rational(1, 4), sp.Rational(4, 3), 1]),
                       (sp.Rational(1, 2), [sp.Rational(1, 4), 1, sp.Rational(1, 6)])]:
    assert all(sp.simplify(bonne - l) != 0 for l in leurres)


# borel-cantelli-ex3 (et v1, v2) : ⌈log_b k⌉ calculé en entiers ; comparaisons à 1/(2k), 1/k², 1/(3k).
def plafond_log(b, kk, facteur=1):
    """Plus petit entier m tel que b^m >= kk^facteur."""
    m = 0
    while b**m < kk**facteur:
        m += 1
    return m


for kk in range(2, 3000):
    l1 = plafond_log(2, kk)
    assert l1 >= 1 and F(1, 2**l1) > F(1, 2 * kk)  # ex3
    l2 = plafond_log(2, kk, 2)  # ⌈2 log_2 k⌉ : 2^m >= k²
    assert F(1, 2**l2) <= F(1, kk * kk)  # v1
    l3 = plafond_log(3, kk)
    assert l3 >= 1 and F(1, 3**l3) > F(1, 3 * kk)  # v2
assert sp.Sum(sp.Rational(1, 2) ** ki, (ki, 1, sp.oo)).doit() == 1 and sp.Sum(sp.Rational(1, 3) ** ki, (ki, 1, sp.oo)).doit() == sp.Rational(1, 2)

# borel-cantelli-ex4 : ⌊k ln k⌋ strictement croissante (k >= 2), minorée par 1 ; Σ 1/k² convergente ; les autres divergent.
vals = [floor(j * log(j)) for j in range(2, 200000)]
assert vals[0] == 1 and all(a < b for a, b in zip(vals, vals[1:]))
assert all((j + 1) * log(j + 1) - j * log(j) >= log(j + 1) > 1 for j in range(2, 1000))
for t in [1 / (2 * ki), 1 / (ki + 1000)]:
    assert sp.Sum(t, (ki, 1, sp.oo)).is_convergent() is sp.false

# borel-cantelli-ex5 (et v1, v2) : encadrement de q^{floor(x)}, identité q^{c ln n} = n^{-c ln(1/q)}, seuil 1/ln(1/q).
for q, inv, bonne, leurres in [(F(1, 2), 2, 1 / mp.log(2), [mp.log(2), 1, 2]),
                               (F(1, 3), 3, 1 / mp.log(3), [mp.log(3), 1, 3]),
                               (F(2, 3), F(3, 2), 1 / mp.log(1.5), [mp.log(1.5), 1 / mp.log(3), 1.5])]:
    for j in range(0, 200):
        xx = mp.mpf(j) / 7
        val = mp.mpf(q) ** int(mp.floor(xx))
        assert mp.mpf(q) ** xx <= val + mp.mpf("1e-25") and val <= mp.mpf(inv) * mp.mpf(q) ** xx + mp.mpf("1e-25")
    for c in [0.3, 1, 2.5]:
        for nn in [2, 10, 1000]:
            assert proche(mp.mpf(q) ** (c * mp.log(nn)), mp.mpf(nn) ** (-c * mp.log(mp.mpf(inv))))
    assert all(abs(bonne - l) > mp.mpf("0.05") for l in leurres)
    assert proche(1 / mp.log(1 / mp.mpf(q)), bonne)

# =====================================================================================================
# Exercices de la carte loi-forte
# =====================================================================================================
# loi-forte-ex1 (et v1, v2) : espérances par transfert.
de = range(1, 7)
assert sum(F(d * d, 6) for d in de) == F(91, 6) and F(91, 6) not in [F(49, 4), F(35, 12), F(7, 2)]
assert sum(F((d - 3) ** 2, 6) for d in de) == F(19, 6) and F(19, 6) not in [F(1, 4), F(35, 12), F(3, 2)]
assert sum(F(1, 6 * d) for d in de) == F(49, 120) and F(49, 120) not in [F(2, 7), F(49, 20), F(1, 6)]
assert sum(F(abs(d - 3), 6) for d in de) == F(3, 2)  # v1 : le leurre 3/2 est E|D - 3|

# loi-forte-ex2 (et v1, v2) : E(C_n) = (E X)^n, limite p.s. de C_n^{1/n} = exp(E ln X) = √(ab).
for a, b, EX, lim, compare, leurres in [(F(1, 2), F(3, 2), 1, sp.sqrt(3) / 2, -1, [1, sp.Rational(3, 4), 0]),
                                        (F(1, 2), F(5, 2), F(3, 2), sp.sqrt(5) / 2, 1, [sp.Rational(3, 2), sp.Rational(5, 4), 1]),
                                        (F(1, 4), F(2), F(9, 8), sp.sqrt(2) / 2, -1, [sp.Rational(9, 8), sp.Rational(1, 2), 1])]:
    assert (a + b) / 2 == EX
    Eln = (sp.log(sp.Rational(a.numerator, a.denominator)) + sp.log(sp.Rational(b.numerator, b.denominator))) / 2
    assert sp.simplify(sp.exp(Eln) - lim) == 0
    assert sp.sign(lim - 1) == compare  # < 1 : C_n -> 0 ; > 1 : C_n -> +inf
    assert all(sp.simplify(lim - l) != 0 for l in leurres)

# loi-forte-ex3 (et v1, v2) : lois, P(|X| >= n), espérance infinie.
assert sp.Sum(2 / ((ki + 1) * (ki + 2)), (ki, 1, sp.oo)).doit() == 1
assert sp.simplify(sp.Sum(2 / ((ki + 1) * (ki + 2)), (ki, ni, sp.oo)).doit() - 2 / (ni + 1)) == 0
assert sp.Sum(3 / ((ki + 2) * (ki + 3)), (ki, 1, sp.oo)).doit() == 1  # v2
assert sp.simplify(sp.Sum(3 / ((ki + 2) * (ki + 3)), (ki, ni, sp.oo)).doit() - 3 / (ni + 2)) == 0
for K in range(1, 60):  # v1 : somme télescopique
    for nn in range(1, K + 1):
        s = mp.fsum(1 / mp.sqrt(j) - 1 / mp.sqrt(j + 1) for j in range(nn, K + 1))
        assert proche(s, 1 / mp.sqrt(nn) - 1 / mp.sqrt(K + 1))
for t in [2 / (ni + 1), 1 / sp.sqrt(ni), 3 / (ni + 2)]:
    assert sp.Sum(t, (ni, 1, sp.oo)).is_convergent() is sp.false  # E|X| = Σ P(|X| >= n) = +inf
# leurres : formules fausses (elles diffèrent de la bonne réponse en n = 2, ou pour tout n).
nn = 2
assert F(2, nn + 1) not in [F(1, nn + 1), F(2, (nn + 1) * (nn + 2)), F(1, nn)]  # ex3
assert F(3, nn + 2) not in [F(3, 2 * (nn + 2)), F(3, (nn + 2) * (nn + 3)), F(1, nn)]  # v2
assert all(not proche(1 / mp.sqrt(m), v) for m in [2, 3] for v in [mp.mpf(1) / m, 1 / mp.sqrt(m) - 1 / mp.sqrt(m + 1), 1 / (2 * mp.sqrt(m))])  # v1
# E(Y) = Σ_{n>=1} P(Y >= n) pour Y à valeurs dans N : contrôle sur une loi finie.
assert sum(v * p for v, p in loiY.items()) == sum(sum(p for v, p in loiY.items() if v >= j) for j in range(1, 6))

# loi-forte-ex4 (et v1, v2) : E(S_n^4) par énumération et leurres.
for loi, formule, leurres in [({1: F(1, 2), -1: F(1, 2)}, lambda m: 3 * m * m - 2 * m, [lambda m: 3 * m * m, lambda m: m * m, lambda m: 3 * m * m + m]),
                              ({-1: F(1, 3), 0: F(1, 3), 1: F(1, 3)}, lambda m: F(4 * m * m - 2 * m, 3), [lambda m: F(4 * m * m, 3), lambda m: F(2 * m, 3), lambda m: F(4 * m * m + 2 * m, 3)]),
                              ({2: F(1, 2), -2: F(1, 2)}, lambda m: 48 * m * m - 32 * m, [lambda m: 48 * m * m, lambda m: 16 * m * m, lambda m: 48 * m * m + 16 * m])]:
    for m in range(1, 6):
        assert E_S4_brut(loi, m) == formule(m)
    assert all(l(2) != formule(2) for l in leurres)
assert F(4 * 4 - 2 * 2, 3) == F(2, 3) * 2 + 3 * 2 * 1 * F(4, 9)  # v1 : expression de la correction en n = 2

# =====================================================================================================
# Exercices de la carte tcl
# =====================================================================================================
# tcl-ex1 (et v1, v2) : z = (s - nμ)/(σ√n) et leurres (sans σ, sans √n, n au lieu de √n).
for (nn, mu, sig, s, z, leurres) in [(100, 2, 2, 220, F(1), [F(1, 10), F(10), F(2)]),
                                     (36, 5, 3, 189, F(1, 2), [F(3, 2), F(3), F(1, 12)]),
                                     (64, 1, 4, 112, F(3, 2), [F(6), F(12), F(3, 16)])]:
    r = int(nn**0.5)
    assert r * r == nn
    assert F(s - nn * mu, sig * r) == z
    assert sorted([F(s - nn * mu, sig * nn), F(s - nn * mu, r), F(s - nn * mu, sig)]) == sorted(leurres)
    assert z not in leurres
assert abs(Phi(1) - mp.mpf("0.841")) < 0.0005 and abs(Phi(0.5) - mp.mpf("0.691")) < 0.0005 and abs(Phi(1.5) - mp.mpf("0.933")) < 0.0005

# tcl-ex2 (et v1, v2) : h = 1,96 σ/√n et leurres (variance, n, double).
for (nn, sig, h, leurres) in [(10000, F(1, 2), F("0.0098"), [F("0.0049"), F("0.000098"), F("0.0196")]),
                              (2500, F(1, 2), F("0.0196"), [F("0.0098"), F("0.000392"), F("0.0392")]),
                              (400, F(3, 10), F("0.0294"), [F("0.00882"), F("0.00147"), F("0.0588")])]:
    r = int(nn**0.5)
    assert F("1.96") * sig / r == h
    assert [F("1.96") * sig**2 / r, F("1.96") * sig / nn, 2 * F("1.96") * sig / r] == leurres
assert abs(2 * Phi(1.96) - 1 - mp.mpf("0.95")) < 0.0005
assert F(1, 10) * F(9, 10) == F(3, 10) ** 2  # v2 : σ = 3/10

# tcl-ex3 (et v1, v2) : tailles par Bienaymé-Tchebychev et par le TCL ; leurres.
def plus_petit_entier(borne):
    return int(mp.ceil(borne))


for (s2, eps, nt, nl, leurres) in [(F(5, 36), F(1, 100), 27778, 5336, [1389, 2723]),
                                  (F(1, 4), F(1, 50), 12500, 2401, [625, 1225]),
                                  (F(4, 25), F(1, 100), 32000, 6147, [1600, 3136])]:
    assert plus_petit_entier(s2 / (F(1, 20) * eps**2)) == nt
    assert plus_petit_entier(F("3.8416") * s2 / eps**2) == nl
    assert [plus_petit_entier(s2 / eps**2), plus_petit_entier(F("1.96") * s2 / eps**2)] == leurres
    assert len({nl, nt, *leurres}) == 4
assert F("1.96") ** 2 == F("3.8416")
assert abs(3.8416 * 5 * 10000 / 36 - 5335.6) < 0.05
assert abs(20 * 5 * 10000 / 36 - 27777.8) < 0.05 and F("3.8416") * F(4, 25) * 10000 == F("6146.56")
assert F("3.8416") * F(1, 4) * 2500 == 2401

# tcl-ex4 (et v1, v2) : E = 0, V = 1, sommes des P(X_k != 0) et leurres.
for j in range(1, 40):
    assert F(j * j) * 2 * F(1, 2 * j * j) == 1 and F(j**4) * 2 * F(1, 2 * j**4) == 1 and F(4**j) * 2 * F(1, 2 * 4**j) == 1
    assert 2 * F(1, 2 * j * j) <= 1
assert sp.Sum(1 / ki**2, (ki, 1, sp.oo)).doit() == sp.pi**2 / 6
assert sp.Sum(1 / ki**4, (ki, 1, sp.oo)).doit() == sp.pi**4 / 90
assert sp.Sum(sp.Rational(1, 4) ** ki, (ki, 1, sp.oo)).doit() == sp.Rational(1, 3)
assert sp.pi**2 / 6 != 1 and sp.pi**4 / 90 != sp.pi**2 / 6 and sp.Rational(1, 3) not in [sp.Rational(4, 3), sp.Rational(1, 4)]
assert abs(Phi(-1) - mp.mpf("0.159")) < 0.0005

# tcl-ex5 (et v1, v2) : u_n = P(Y_n <= m), Y_n de loi P(n) ; valeurs exactes pour n = 10^6 (convergence en 1/√n).
def poisson_cdf(lam_, m):
    return mp.gammainc(m + 1, lam_, mp.inf, regularized=True)


nn = 10**6
assert abs(poisson_cdf(nn, nn) - mp.mpf(1) / 2) < mp.mpf("0.001")  # ex5 : -> 1/2
assert abs(poisson_cdf(nn, nn + 1000) - Phi(1)) < mp.mpf("0.001")  # v1 : -> Φ(1)
assert abs(poisson_cdf(nn, nn - 2000) - (1 - Phi(2))) < mp.mpf("0.001")  # v2 : -> 1 - Φ(2)
assert proche(poisson_cdf(5, 3), mp.fsum(mp.e**-5 * mp.mpf(5) ** j / mp.factorial(j) for j in range(4)))
assert abs(Phi(1) - mp.mpf(1) / 2 - Phi(1)) > 0.4 and abs((1 - Phi(2)) - Phi(2)) > 0.9

# =====================================================================================================
# Exercices de la carte binomiale
# =====================================================================================================
# binomiale-ex1 (et v1, v2) : approximations de Poisson et valeurs exactes.
assert abs(mp.e**-2 - mp.mpf("0.1353")) < 0.00005 and abs(mp.mpf("0.998") ** 1000 - mp.mpf("0.1351")) < 0.00005
ap1 = 3 * mp.e**-2
ex1 = mp.mpf("0.996") ** 500 + 500 * mp.mpf("0.004") * mp.mpf("0.996") ** 499
assert abs(ap1 - mp.mpf("0.4060")) < 0.00005 and abs(ex1 - mp.mpf("0.4057")) < 0.00005
ap2 = mp.mpf(9) / 2 * mp.e**-3
ex2 = binom_pmf(2000, mp.mpf("0.0015"), 2)
assert abs(ap2 - mp.mpf("0.2240")) < 0.00005 and abs(ex2 - mp.mpf("0.2241")) < 0.00005
assert len({float(v) for v in [mp.e**-2, mp.e ** -0.002, 2 * mp.e**-2, 1 - mp.e**-2]}) == 4
assert len({float(v) for v in [3 * mp.e**-2, 2 * mp.e**-2, mp.e**-2, 1 - 2 * mp.e**-2]}) == 4
assert len({float(v) for v in [4.5 * mp.e**-3, 9 * mp.e**-3, 3 * mp.e**-3, 1.5 * mp.e**-3]}) == 4

# binomiale-ex2 (et v1, v2) : correction de continuité, valeurs exactes et leurres.
assert F(57.5 - 50) / 5 == F(3, 2) and [F(57 - 50, 5), F(56.5 - 50) / 5, F(57.5 - 50) / 25] == [F(7, 5), F(13, 10), F(3, 10)]
Pex = binom_cdf(100, mp.mpf(1) / 2, 57)
assert abs(Pex - mp.mpf("0.933")) < 0.0005 and abs(Phi(1.5) - mp.mpf("0.933")) < 0.0005 and abs(Phi(1.4) - mp.mpf("0.919")) < 0.0005
assert F(209.5 - 200) / 10 == F(19, 20) and [F(210 - 200, 10), F(210.5 - 200) / 10, F(209.5 - 200) / 100] == [1, F(21, 20), F(19, 200)]
Pex = 1 - binom_cdf(400, mp.mpf(1) / 2, 209)
assert abs(Pex - mp.mpf("0.171")) < 0.0005 and abs(1 - Phi(0.95) - mp.mpf("0.171")) < 0.0005 and abs(1 - Phi(1) - mp.mpf("0.159")) < 0.0005
assert F(94.5 - 90) / 9 == F(1, 2) and [F(94 - 90, 9), F(93.5 - 90) / 9, F(94.5 - 90) / 81] == [F(4, 9), F(7, 18), F(1, 18)]
Pex = binom_cdf(900, mp.mpf(1) / 10, 94)
assert abs(Pex - mp.mpf("0.695")) < 0.0005 and abs(Phi(0.5) - mp.mpf("0.691")) < 0.0005 and abs(Phi(mp.mpf(4) / 9) - mp.mpf("0.672")) < 0.0005
assert (900 * F(1, 10) * F(9, 10)) == 81 and (400 * F(1, 4)) == 100 and (100 * F(1, 4)) == 25

# binomiale-ex3 (et v1, v2) : (1 - λ/n)^n <= e^{-λ} et n(e^{-λ} - (1 - λ/n)^n) -> λ² e^{-λ}/2.
nsym = sp.Symbol("n", positive=True)
for lam_, lim, leurres in [(2, 2 * sp.exp(-2), [4 * sp.exp(-2), sp.exp(-2), 0]),
                           (3, sp.Rational(9, 2) * sp.exp(-3), [9 * sp.exp(-3), 3 * sp.exp(-3), 0]),
                           (1, 1 / (2 * sp.E), [1 / sp.E, sp.Rational(1, 2), 0])]:
    assert sp.simplify(sp.limit(nsym * (sp.exp(-lam_) - (1 - sp.Integer(lam_) / nsym) ** nsym), nsym, sp.oo) - lim) == 0
    assert sp.series(nsym * sp.log(1 - lam_ / nsym), nsym, sp.oo, 2).removeO().equals(-lam_ - sp.Rational(lam_**2, 2) / nsym)
    for m in range(lam_ + 1, 200):
        assert (1 - mp.mpf(lam_) / m) ** m <= mp.e ** (-lam_)
    assert all(sp.simplify(lim - l) != 0 for l in leurres)

# binomiale-ex4 (et v1, v2) : plus grand N (approximation normale, z = 1,645, correction de continuité),
# contrôle exact, et leurres (z = 1,96 ; arrondi supérieur de la racine ; N = places/p).
def plus_grand_N(s, p, z, corr=0.5):
    return max(NN for NN in range(s, 3 * s) if (s + corr - NN * p) / (NN * p * (1 - p)) ** 0.5 >= z)


def P_exacte_depasse(NN, p, s):
    return 1 - binom_cdf(NN, p, s)


for (s, p, nmax, racine, r2, leurres, ex_a, ex_b) in [(100, 0.8, 116, 10.805, 116.7, [117, 115, 125], "0.032", "0.051"),
                                                      (300, 0.95, 309, 17.598, 309.7, [310, 308, 315], "0.027", "0.051"),
                                                      (400, 0.9, 433, 20.823, 433.6, [434, 431, 444], "0.038", "0.052")]:
    assert plus_grand_N(s, p, 1.645) == nmax
    a, b, c = p, 1.645 * (p * (1 - p)) ** 0.5, -(s + 0.5)
    rac = (-b + (b * b - 4 * a * c) ** 0.5) / (2 * a)
    assert abs(rac - racine) < 0.0005 and abs(rac * rac - r2) < 0.05 and floor(rac * rac) == nmax
    assert [nmax + 1, plus_grand_N(s, p, 1.96), floor(s / p)] == leurres
    pa, pb = P_exacte_depasse(nmax, mp.mpf(str(p)), s), P_exacte_depasse(nmax + 1, mp.mpf(str(p)), s)
    assert pa <= 0.05 < pb and abs(pa - mp.mpf(ex_a)) < 0.0005 and abs(pb - mp.mpf(ex_b)) < 0.0005
assert abs(1.645 * 0.4 - 0.658) < 1e-12 and abs(1.645 * 0.3 - 0.4935) < 1e-12 and abs(1.645 * 0.0475**0.5 - 0.3585) < 0.00005
assert abs(Phi(1.645) - mp.mpf("0.95")) < 0.0005

# binomiale-ex5 : |x| = 2x⁺ - x, et contrôle numérique de Σ|a_n(k) - q(k)| -> 0 pour p_n = 2/n.
for xx in [F(-3, 2), F(0), F(5, 7)]:
    assert abs(xx) == 2 * max(xx, 0) - xx
tv = []
for nn in [50, 500, 5000]:
    a_ = [binom_pmf(nn, mp.mpf(2) / nn, j) for j in range(0, 80)]
    q_ = [mp.e**-2 * mp.mpf(2) ** j / mp.factorial(j) for j in range(0, 80)]
    tv.append(mp.fsum(abs(u - v) for u, v in zip(a_, q_)))
    assert proche(mp.fsum(abs(u - v) for u, v in zip(a_, q_)), 2 * mp.fsum(max(v - u, 0) for u, v in zip(a_, q_)), mp.mpf("1e-20"))
assert tv[0] > tv[1] > tv[2] and tv[2] < mp.mpf("0.001")

# =====================================================================================================
# QCM
# =====================================================================================================
# q1 : X_n = 1/n (voir la carte des convergences).
# q3 : (4/100)/(1/4) = 4/25 ; leurres : ε non élevé au carré (2/25), n² (1/625), sans ε² (1/25).
assert F(4, 100) / F(1, 4) == F(4, 25)
assert [F(4, 100) / F(1, 2), F(4, 100 * 100) / F(1, 4), F(4, 100)] == [F(2, 25), F(1, 625), F(1, 25)]
# q5 : Σ 1/(n+1) diverge (vérifié plus haut).
# q7 : (√n)/(√n/2) = 2 ; contrôle exact pour n = 10^4 : P(S_n <= 5100) ≈ Φ(2).
assert abs(binom_cdf(10**4, mp.mpf(1) / 2, 5100) - Phi(2)) < mp.mpf("0.005")
assert len({float(v) for v in [Phi(2), Phi(1), Phi(0.5), 0.5]}) == 4
# q8 : n·(3/n)(1 - 3/n)^(n-1) = 3(1 - 3/n)^(n-1) -> 3e^{-3}.
assert abs(3 * (1 - mp.mpf(3) / 10**8) ** (10**8 - 1) - 3 * mp.e**-3) < mp.mpf("1e-6")
assert len({float(v) for v in [3 * mp.e**-3, mp.e**-3, 3 * mp.e**-1, 1 - mp.e**-3]}) == 4

print("proba_l3_limites : toutes les vérifications passent.")
