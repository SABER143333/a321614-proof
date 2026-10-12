"""Dump the explicit objects of the A321614 proof, for typesetting into the note.

Nothing is computed here that the note does not state in full; this is purely
so the 12x12 matrix, the state list and the three weight vectors can be printed
exactly rather than transcribed by hand.
"""
import kings as K

ALLOWED = [s for s in K.STATES if K.allowed(s)]
IDX = {s: i for i, s in enumerate(ALLOWED)}
print("allowed column states (12), in index order:")
for i, s in enumerate(ALLOWED):
    print(f"  s{i:<2} = {s}")

print("\nM[i][j] = 1 iff trans(s_i, s_j);  M =")
M = K.build_M()
for i, row in enumerate(M):
    print("  [" + ", ".join(str(x) for x in row) + "],")


def vec(name, pred):
    v = [1 if pred(s) else 0 for s in ALLOWED]
    print(f"{name} = {v}")
    return v


print()
vec("J_h  (trans(s, sigma_h(s)))", lambda s: K.trans(s, K.sigma_h(s)))
vec("J_re (trans(s, sigma_rot(s)))", lambda s: K.trans(s, K.sigma_rot(s)))
vec("sigma_rot-fixed indicator",
    lambda s: K.sigma_rot(s) == s)
print("C1 = #sigma_rot-fixed allowed states =",
      sum(1 for s in ALLOWED if K.sigma_rot(s) == s))
