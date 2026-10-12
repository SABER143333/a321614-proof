"""Close the Burnside step of the A321614 proof COMPONENT BY COMPONENT.

The existing evidence (22 published terms; brute force n<=6; independent
column-profile DP n<=26) only checks the SUM  T + Fh + Fv + Fr  divided by 4.
That leaves the classic failure mode wide open: two components can be wrong in
opposite directions and cancel, or all four can be wrong together while their
average still matches at every computed point.

Burnside's lemma is a statement about each FIXED-POINT COUNT separately, so
each of T(n), Fh(n), Fv(n), Fr(n) must be checked individually against a
definition-level enumeration of placements.  That is what this file does.

If every component matches for n = 1..5, then the only remaining content of
the Burnside step is the transfer-matrix claim that M^{n-1} counts column
walks -- which verify_model.py already closed exhaustively.
"""
import sys
import kings as K

print("component-by-component Burnside check")
print(f"{'n':>3} | {'component':<4} {'brute':>8} {'formula':>8}  ok")
print("-" * 40)

allok = True
for n in range(1, 6):
    P = K.brute_placements(n)

    # definition-level fixed-point counts, straight from the board
    t = len(P)
    fh = sum(1 for p in P if K.reflect_h(p, n) == p)
    fv = sum(1 for p in P if K.reflect_v(p) == p)
    fr = sum(1 for p in P if K.rot180(p, n) == p)

    # the formula's values
    ft = K.total_placements(n)
    fh_ = K.fix_h(n)
    fv_ = K.fix_v(n)
    fr_ = K.fix_rot(n)

    for name, b, f in (("T", t, ft), ("Fh", fh, fh_),
                       ("Fv", fv, fv_), ("Fr", fr, fr_)):
        ok = (b == f)
        allok &= ok
        print(f"{n:>3} | {name:<4} {b:>8} {f:>8}  {'OK' if ok else 'MISMATCH'}")

    # and the Burnside average
    avg = (t + fh + fv + fr) // 4
    favg = K.a_fast(n)
    ok = (avg == favg) and (avg == K.B_FILE[n])
    allok &= ok
    print(f"{n:>3} | {'a(n)':<4} {avg:>8} {favg:>8}  "
          f"{'OK' if ok else 'MISMATCH'}   (b-file {K.B_FILE[n]})")
    print()

print("EVERY COMPONENT MATCHES" if allok else "*** COMPONENT MISMATCH ***")
sys.exit(0 if allok else 1)
