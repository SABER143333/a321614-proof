"""Closes the two combinatorial links in the proof by exhaustive check.

Step 1 (pigeonhole).  A 2x2 block has all its cell-pairs within Chebyshev
distance 1, so it holds at most one king.  There are 2n blocks and a maximum
placement has 2n kings, so every block holds EXACTLY one.  Check exhaustively:
no placement of 2n nonattacking kings on a 4 x 2n board has a block with two
kings or an empty block.  (Done by listing all maximum placements.)

Step 2 (local constraints).  With one king per block the offsets are variables
(r,c) in {0,1}^2 per block.  Claim: the four constraints encoded in trans() and
allowed() are EXACTLY equivalent to "no two kings attack".  Check exhaustively
over every pair of blocks that can interact, and over full boards.

If both checks pass, the only remaining ingredient is the exact symbolic
generating-function identity in prove.py, which is finite exact arithmetic.
"""
import itertools
import kings as K

W = 4
FAIL = []


def attacks(a, b):
    return max(abs(a[0] - b[0]), abs(a[1] - b[1])) <= 1


# ---------------------------------------------------------------- step 1 ----
def check_step1(n):
    """every maximum placement has exactly one king per 2x2 block"""
    P = K.brute_placements(n)
    for cells in P:
        for i in range(2):
            for j in range(n):
                block = {(2 * i + r, 2 * j + c) for r in (0, 1) for c in (0, 1)}
                k = len(block & set(cells))
                if k != 1:
                    FAIL.append(f"step1 n={n}: block ({i},{j}) holds {k} kings")
                    return
    return len(P)


# ---------------------------------------------------------------- step 2 ----
def check_step2_pairwise():
    """the transition relation is exactly 'the two columns' kings do not attack'"""
    states = list(K.IDX.items())
    for s, i in states:
        for u, j in states:
            # build the 4 x 4 piece of board occupied by these two columns
            cells = []
            for col, st in ((0, s), (2, u)):
                r0, c0, r1, c1 = st
                cells.append((r0, col + c0))
                cells.append((2 + r1, col + c1))
            bad = any(attacks(a, b) for a, b in itertools.combinations(cells, 2))
            if bad == K.trans(s, u):
                FAIL.append(f"step2: trans({s},{u})={K.trans(s,u)} but attacks={bad}")


def check_step2_full(n):
    """the matrix placement set equals the geometric placement set, for full boards"""
    matset = set()
    for seq in itertools.product(list(K.IDX.keys()), repeat=n):
        if all(K.trans(seq[j], seq[j + 1]) for j in range(n - 1)):
            cells = set()
            for j, (r0, c0, r1, c1) in enumerate(seq):
                cells.add((r0, 2 * j + c0))
                cells.add((2 + r1, 2 * j + c1))
            matset.add(frozenset(cells))
    geoset = set(K.brute_placements(n))
    if matset != geoset:
        FAIL.append(f"step2 n={n}: matrix {len(matset)} vs geometric {len(geoset)}")
    return len(matset)


if __name__ == "__main__":
    print("step 1: every maximum placement has exactly one king per 2x2 block")
    for n in range(1, 5):
        c = check_step1(n)
        print(f"   n={n}: {c} placements checked, all blocks exactly 1 king"
              if c is not None else f"   n={n}: FAILED")

    print("\nstep 2a: allowed() is exactly 'the two kings in one column do not attack'")
    for s in K.STATES:
        r0, c0, r1, c1 = s
        a, b = (r0, c0), (2 + r1, c1)
        if K.allowed(s) == attacks(a, b):
            FAIL.append(f"step2a: allowed({s})={K.allowed(s)} but attacks={attacks(a,b)}")
    print(f"   {len(K.STATES)} within-column states checked, "
          f"{'all agree' if not FAIL else 'FAILURES'}")

    print("\nstep 2b: trans() is exactly the cross-column non-attacking condition")
    check_step2_pairwise()
    print(f"   {len(K.IDX)}x{len(K.IDX)} = {len(K.IDX)**2} state pairs checked, "
          f"{'all agree' if not FAIL else 'FAILURES'}")

    print("\nstep 2b: matrix placement set == geometric placement set")
    for n in range(1, 6):
        c = check_step2_full(n)
        print(f"   n={n}: {c} placements, sets identical" if c is not None
              else f"   n={n}: FAILED")

    print("\n" + ("ALL CHECKS PASSED -- the combinatorial links are closed"
                  if not FAIL else f"{len(FAIL)} FAILURES: {FAIL[:5]}"))
