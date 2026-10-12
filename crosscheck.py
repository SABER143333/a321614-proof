"""Structurally different second implementation, as SciNet's acceptance
criteria demand.

kings.py exploits the fact that a maximum placement puts exactly one king in
each 2x2 block and builds a 12-state transfer matrix on block offsets.

This file knows nothing about blocks.  It is the generic grid-independent-set
column-profile DP: sweep the 4 x 2n board column by column, the state is the
4-bit mask of occupied cells in the current column, and the count is tracked as
a polynomial in the number of kings.  We read off the coefficient of x^(2n).

If the two implementations agree, the shared assumption they could both be
wrong about is essentially nil.
"""
import itertools, sys

W = 4                      # board height


def masks_ok():
    return [m for m in range(1 << W) if not (m & (m >> 1))]


MASKS = masks_ok()


def compatible(m1, m2):
    """m1 = previous column, m2 = current column; no king attacks across"""
    if m1 & m2:
        return False
    if (m1 << 1) & m2:
        return False
    if (m1 >> 1) & m2:
        return False
    return True


def T_cell(n, cap=None):
    """number of ways to place exactly 2n nonattacking kings on 4 x 2n"""
    if n == 0:
        return 1
    cols = 2 * n
    top = cap if cap is not None else cols
    cur = {m: [0] * (top + 1) for m in MASKS}
    for m in MASKS:
        cur[m][bin(m).count("1")] = 1
    for _ in range(cols - 1):
        nxt = {m: [0] * (top + 1) for m in MASKS}
        for m2 in MASKS:
            k2 = bin(m2).count("1")
            t = nxt[m2]
            for m1 in MASKS:
                if compatible(m1, m2):
                    v = cur[m1]
                    for k in range(top + 1 - k2):
                        if v[k]:
                            t[k + k2] += v[k]
        cur = nxt
    return sum(cur[m][cols] for m in MASKS)


def a_cell(n):
    """Burnside over the rectangle group, with every fixed-count obtained by
    direct enumeration over the profile DP state space."""
    if n == 0:
        return 1
    cols = 2 * n
    total = T_cell(n)
    fh = fixed_under(n, lambda r, c: (r, cols - 1 - c)) if n % 2 == 0 else 0
    fv = fixed_under(n, lambda r, c: (W - 1 - r, c))
    fr = fixed_under(n, lambda r, c: (W - 1 - r, cols - 1 - c))
    return (total + fh + fv + fr) // 4


def fixed_under(n, g):
    """count maximum placements invariant under the cell map g, by direct
    construction: build the orbit-constrained DP over columns."""
    cols = 2 * n
    # a placement is a set of cells; g acts on cells.  Enumerate by DP over
    # columns, but only over placements that are g-invariant.  For the maps
    # used here an invariant placement is determined by its left half, so we
    # enumerate the left half by profile DP and test the constraint directly.
    halfcols = cols // 2
    count = 0
    # iterate over left-half column masks
    for left in itertools.product(MASKS, repeat=halfcols):
        cells = set()
        for j, m in enumerate(left):
            for r in range(W):
                if m >> r & 1:
                    cells.add((r, j))
        # mirror the left half to get the full placement
        full = set(cells)
        for (r, c) in cells:
            full.add(g(r, c))
        if any(g(a) not in full for a in full):
            continue
        if len(full) != 2 * n:
            continue
        if any(max(abs(a[0] - b[0]), abs(a[1] - b[1])) <= 1
               for a, b in itertools.combinations(full, 2)):
            continue
        count += 1
    return count


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    import kings as K
    print(f"{'n':>3} {'blocks T(n)':>16} {'cell DP T(n)':>16}  match")
    print("-" * 48)
    bad = 0
    for n in range(1, N + 1):
        a = K.total_placements(n)
        b = T_cell(n)
        ok = a == b
        bad += (not ok)
        print(f"{n:>3} {a:>16} {b:>16}  {ok}")
    print(f"\n{bad} disagreements between the two independent implementations")
