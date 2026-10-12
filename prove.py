"""STRONGER tier: proof of Barker's recurrence for all n, from the transfer structure.

Route
-----
1. A maximum placement of 2n nonattacking kings on a 4 x 2n board has exactly
   one king in each of the 2n disjoint 2x2 blocks  (pigeonhole: a block holds at
   most one king, and there are 2n blocks for 2n kings).
2. The non-attacking condition is then exactly the four local constraints
   encoded by trans()  (case analysis on the two blocks' offsets).
3. Hence with M the 12x12 transfer matrix and Burnside over the rectangle
   group {id, h, v, rot180},
        a(n) = ( T(n) + Fh(n) + Fv(n) + Fr(n) ) / 4   (n >= 1)
   where T, Fh, Fr are 1^T M^(k) w for explicit vectors and Fv(n) = n+1.
4. Each of those is a linear recurrence sequence, so a(n) is one.  Computing
   the rational generating function exactly (no floating point, no fitting)
   proves the recurrence for ALL n.

Everything below is exact integer/rational arithmetic.
"""
import sympy as sp
import kings as K

x, y = sp.symbols("x y")
N = K.K
M = sp.Matrix(K.build_M())
ONE = [1] * N
J_h = [1 if K.trans(s, K.sigma_h(s)) else 0 for s, _ in sorted(K.IDX.items(),
                                                              key=lambda t: t[1])]
J_re = [1 if K.trans(s, K.sigma_rot(s)) else 0
        for s, _ in sorted(K.IDX.items(), key=lambda t: t[1])]
J_r = [0] * N
for s, i in K.IDX.items():
    J_r[i] = sum(1 for m, j in K.IDX.items()
                 if K.sigma_rot(m) == m and K.trans(s, m) and K.trans(m, K.sigma_rot(s)))
C1 = sum(1 for s in K.IDX if K.sigma_rot(s) == s)


def rational_from_series(vec_u, vec_w, var):
    """exact p/q with sum_{k>=0} (u^T M^k w) var^k = p(var)/q(var), deg p < deg q"""
    q = sp.expand((sp.eye(N) - var * M).det())
    dq = sp.degree(q, var)
    u = sp.Matrix(1, N, vec_u)
    w = sp.Matrix(N, 1, vec_w)
    coeffs = []
    P = sp.eye(N)
    for k in range(dq):
        coeffs.append(sp.expand((u * P * w)[0, 0]))
        P = P * M
    # p = q * (infinite series).  The product is a polynomial because c_k obeys
    # the recurrence for k >= deg q; take its first deg q coefficients.
    qc = sp.Poly(q, var).all_coeffs()[::-1]        # qc[i] = coefficient of var^i
    dq_q = len(qc)
    p_coeffs = []
    for j in range(dq):
        s = 0
        for i in range(min(j, dq_q - 1) + 1):
            if j - i < len(coeffs):
                s += qc[i] * coeffs[j - i]
        p_coeffs.append(sp.expand(s))
    p = sum(p_coeffs[k] * var**k for k in range(dq))
    return sp.expand(p), q


# ---- the four pieces, each as an exact rational function of x ----
q_x = sp.expand((sp.eye(N) - x * M).det())
p_T, q_T = rational_from_series(ONE, ONE, x)
P_T = sp.expand(x * p_T / q_T)                     # sum_{n>=1} T(n) x^n : T(n)=1^T M^{n-1} 1

q_y = sp.expand((sp.eye(N) - y * M).det())
p_h, _ = rational_from_series(ONE, J_h, y)
p_re, _ = rational_from_series(ONE, J_re, y)
p_r, _ = rational_from_series(ONE, J_r, y)

P_Fh = sp.expand(x**2 * p_h.subs(y, x**2) / q_y.subs(y, x**2))
P_Fre = sp.expand(x**2 * p_re.subs(y, x**2) / q_y.subs(y, x**2))
P_Fro = sp.expand(C1 * x + x**3 * p_r.subs(y, x**2) / q_y.subs(y, x**2))
P_Fv = sp.expand(1 / (1 - x)**2 - 1)               # sum_{n>=1} (n+1) x^n

A = sp.together(1 + (P_T + P_Fh + P_Fv + P_Fre + P_Fro) / 4)
A = sp.cancel(sp.together(A))
num, den = sp.fraction(A)
num = sp.expand(num)
den = sp.expand(den)
print("computed generating function:")
print("  numerator  :", num)
print("  denominator:", den)
print("  factor(num):", sp.factor(num))
print("  factor(den):", sp.factor(den))

# ---- Barker's conjectured generating function ----
B_NUM = sp.expand((1 - 2*x) * (1 - 6*x + 17*x**2 - 18*x**3 - 2*x**4 + 7*x**5
                               + 6*x**6 - 3*x**7))
B_DEN = sp.expand((1 - x)**2 * (1 - 3*x)**2 * (1 - 3*x + x**2)
                  * (1 - x - x**2) * (1 - 3*x**2))
print("\nBarker numerator  :", B_NUM)
print("Barker denominator:", B_DEN)

diff = sp.cancel(sp.together(A - sp.together(B_NUM / B_DEN)))
print("\nA(x) - Barker(x) =", diff)
print("\n*** GENERATING FUNCTIONS ARE IDENTICAL:", sp.simplify(diff) == 0, "***")

# ---- and a direct series check ----
ser = sp.series(A, x, 0, 30).removeO()
got = [sp.expand(ser).coeff(x, n) for n in range(30)]
print("\nfirst 30 coefficients of the computed G.f.:")
print(got[:14], "...")
print("match b-file:", all(got[n] == K.B_FILE[n] for n in range(22)))
