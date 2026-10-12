"""A321614 -- maximum nonattacking kings on a 4 x 2n board, counted up to
the symmetry group of the rectangle.

Structure used
--------------
A 4 x 2n board splits into 2n disjoint 2x2 blocks (2 row-pairs x n column-pairs).
Any two cells of one block are within Chebyshev distance 1, so a block holds at
most one king.  Since 2n kings is the maximum and there are exactly 2n blocks,
EVERY maximum placement has exactly one king per block.  So a placement is a
choice, for each block (i,j), of an offset (r,c) in {0,1}^2.

Local constraints (all verified by the brute force below):
  within a column  : not(r0=1 and r1=0)
  same row-pair, adjacent columns : c_{i,j} <= c_{i,j+1}
  diagonal (i,j)-(1-i,j+1)        : not(r=1 and r'=0 and c=1 and c'=0)

Counting is then a 16-state transfer matrix over columns.
"""
import itertools, json, os, sys
from math import gcd

ROOT = os.path.dirname(os.path.abspath(__file__))
STATES = [(r0, c0, r1, c1) for r0 in (0, 1) for c0 in (0, 1)
          for r1 in (0, 1) for c1 in (0, 1)]


def allowed(s):
    r0, c0, r1, c1 = s
    return not (r0 == 1 and r1 == 0)


def trans(s, u):
    """valid transition from column-state s to column-state u (u is the next column)"""
    r0, c0, r1, c1 = s
    R0, C0, R1, C1 = u
    if c0 > C0 or c1 > C1:
        return False
    # (0,j) vs (1,j+1): row gap 2+R1-r0 == 1  iff r0=1,R1=0
    #                   col gap 2+C1-c0 == 1  iff c0=1,C1=0
    if r0 == 1 and R1 == 0 and c0 == 1 and C1 == 0:
        return False
    # (1,j) vs (0,j+1): row gap 2+r1-R0 == 1  iff r1=0,R0=1
    #                   col gap 2+c1-C0 == 1  iff c1=1,C0=0
    if r1 == 0 and R0 == 1 and c1 == 1 and C0 == 0:
        return False
    return True


# ---------------------------------------------------------------- brute force
def brute_placements(n):
    """enumerate every maximum placement as a frozenset of (row,col)"""
    blocks = [(i, j) for i in range(2) for j in range(n)]
    out = []
    for choice in itertools.product([(r, c) for r in (0, 1) for c in (0, 1)],
                                    repeat=2 * n):
        cells = set()
        for (i, j), (r, c) in zip(blocks, choice):
            cells.add((2 * i + r, 2 * j + c))
        ok = True
        for a, b in itertools.combinations(cells, 2):
            if max(abs(a[0] - b[0]), abs(a[1] - b[1])) <= 1:
                ok = False
                break
        if ok:
            out.append(frozenset(cells))
    return out


def reflect_h(cells, n):
    return frozenset((r, 2 * n - 1 - c) for r, c in cells)


def reflect_v(cells):
    return frozenset((3 - r, c) for r, c in cells)


def rot180(cells, n):
    return frozenset((3 - r, 2 * n - 1 - c) for r, c in cells)


def brute_a(n):
    """Burnside over the rectangle symmetry group, computed by brute force"""
    if n == 0:
        return 1
    P = brute_placements(n)
    # Burnside counts placements FIXED by the symmetry, not merely mapped to
    # another valid placement (the latter is always true and gives |P|).
    fh = sum(1 for p in P if reflect_h(p, n) == p)
    fv = sum(1 for p in P if reflect_v(p) == p)
    fr = sum(1 for p in P if rot180(p, n) == p)
    # symmetry group of an m x n rectangle with m != n is Z2 x Z2 (order 4)
    return (len(P) + fh + fv + fr) // 4


# ------------------------------------------------------------ transfer matrix
ALLOWED = [s for s in STATES if allowed(s)]
IDX = {s: i for i, s in enumerate(ALLOWED)}          # compact 0..K-1 indices
K = len(ALLOWED)


def mat():
    M = [[0] * K for _ in range(K)]
    for s, i in IDX.items():
        for u, j in IDX.items():
            if trans(s, u):
                M[i][j] += 1
    return M


def matmul(A, B, mod=None):
    n = len(A)
    C = [[0] * n for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        Ci = C[i]
        for k in range(n):
            a = Ai[k]
            if a:
                Bk = B[k]
                for j in range(n):
                    Ci[j] += a * Bk[j]
        if mod:
            Ci[:] = [x % mod for x in Ci]
    return C


def matvec(A, v, mod=None):
    n = len(A)
    out = [0] * n
    for i in range(n):
        s = 0
        Ai = A[i]
        for j in range(n):
            if Ai[j]:
                s += Ai[j] * v[j]
        out[i] = s % mod if mod else s
    return out


def flip_h_state(s):
    r0, c0, r1, c1 = s
    return (r0, 1 - c0, r1, 1 - c1)


def build_M():
    return mat()


def total_placements(n, mod=None):
    """T(n) = all maximum placements, no symmetry quotient"""
    if n == 0:
        return 1
    M = build_M()
    v = [1] * K
    for _ in range(n - 1):
        v = matvec(M, v, mod)
    return sum(v) % mod if mod else sum(v)


def sigma_h(s):
    """state of the left-right mirrored column"""
    r0, c0, r1, c1 = s
    return (r0, 1 - c0, r1, 1 - c1)


def sigma_rot(s):
    """state of the 180-degree rotated column (row-pair and column index both flip)"""
    r0, c0, r1, c1 = s
    return (1 - r1, 1 - c1, 1 - r0, 1 - c0)


def fix_colreverse(n, sigma, mid_ok, mod=None):
    """placements invariant under a symmetry that reverses the column index.

    Columns j and n-1-j carry states related by sigma.  Walking from the middle
    outward: the last free column must connect to its image, and each earlier
    column is prepended.  matvec PREPENDS a state, so we start from the weight
    vector on the junction column rather than from a vector of ones.
    """
    if n == 1:
        return sum(1 for s in IDX if mid_ok(s)) % mod if mod else \
            sum(1 for s in IDX if mid_ok(s))
    if n % 2 == 1:
        h = (n - 1) // 2
        J = [0] * K
        for s, i in IDX.items():
            t = 0
            for m, j in IDX.items():
                # the middle column must connect on BOTH sides to its neighbour
                if mid_ok(m) and trans(s, m) and trans(m, sigma(s)):
                    t += 1
            J[i] = t
    else:
        h = n // 2
        J = [0] * K
        for s, i in IDX.items():
            if trans(s, sigma(s)):
                J[i] = 1
    M = build_M()
    u = J[:]
    for _ in range(h - 1):
        u = matvec(M, u, mod)
    tot = sum(u)
    return tot % mod if mod else tot


def fix_h(n, mod=None):
    """left-right reflection: middle column would need c = 1-c, so n odd gives 0"""
    return fix_colreverse(n, sigma_h, lambda s: False, mod)


def fix_v(n, mod=None):
    """top-bottom reflection: every column must be sigma_v-fixed -> n+1 of them"""
    return (n + 1) % mod if mod else n + 1


def fix_rot(n, mod=None):
    """180-degree rotation"""
    return fix_colreverse(n, sigma_rot, lambda s: sigma_rot(s) == s, mod)


def a_fast(n, mod=None):
    if n == 0:
        return 1
    tot = total_placements(n, mod) + fix_h(n, mod) + fix_v(n, mod) + fix_rot(n, mod)
    if mod:
        return (tot * pow(4, -1, mod)) % mod
    return tot // 4


B_FILE = [1, 4, 23, 106, 473, 1939, 7618, 28703, 105112, 375597, 1316944,
          4544124, 15474559, 52108212, 173799309, 574908646, 1888125243,
          6162032375, 19998659760, 64584817367, 207655073310, 665017743665]

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "check"

    if mode == "brute":
        print("brute force vs transfer matrix vs published b-file")
        print(f"{'n':>3} {'brute':>12} {'fast':>12} {'b-file':>12}  match")
        for n in range(0, 6):
            b = brute_a(n) if n <= 5 else None
            f = a_fast(n)
            ref = B_FILE[n] if n < len(B_FILE) else None
            ok = (b in (None, f)) and (ref in (None, f))
            print(f"{n:>3} {str(b):>12} {f:>12} {str(ref):>12}  {ok}")

    elif mode == "check":
        print("regression gate: transfer matrix vs the 22 published b-file terms")
        bad = 0
        for n, ref in enumerate(B_FILE):
            f = a_fast(n)
            flag = "OK" if f == ref else "MISMATCH"
            if f != ref:
                bad += 1
            print(f"  n={n:>2}  got {f:>14}  b-file {ref:>14}  {flag}")
        print(f"\n{bad} mismatches out of {len(B_FILE)}")
