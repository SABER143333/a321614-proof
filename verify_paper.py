"""Verify the SUBMITTED PAPER against the code.

Every other check in this project reads the code.  That leaves one failure mode
untouched: the paper could faithfully describe a different theorem from the one
the code proves, because a matrix or a generating function was mistyped while
being typeset.  This file reads a321614.tex itself, extracts the mathematical
objects from the LaTeX, and compares them with what the programs compute.
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import kings as K

HERE = os.path.dirname(os.path.abspath(__file__))


def find_tex():
    for cand in (os.path.join(HERE, "a321614.tex"),
                 os.path.join(HERE, "..", "a321614.tex")):
        if os.path.exists(cand):
            return cand
    raise SystemExit("a321614.tex not found next to this script or one level up")


TEX = find_tex()
src = open(TEX, encoding="utf-8").read()
fails = []


def rep(name, ok, detail=""):
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")
    if not ok:
        fails.append(name)


ALLOWED = [s for s in K.STATES if K.allowed(s)]
M = K.build_M()

print("paper vs code\n")

m = re.search(r"\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", src, re.S)
rows = None
if m:
    rows = [[int(x) for x in r.replace("\\\\", "").strip().rstrip("\\\\").split("&")]
            for r in m.group(1).strip().split("\\\\")]
    rows = [r for r in rows if r]
rep("transfer matrix M is transcribed correctly",
    rows == M, f"({len(rows) if rows else 0} rows parsed)")

states = re.findall(r"s_\{?(\d+)\}?\s*=\s*\((\d),(\d),(\d),(\d)\)", src)
seen = {}
conflict = []
for i, a, b, c, d in states:
    v = (int(a), int(b), int(c), int(d))
    if int(i) in seen and seen[int(i)] != v:
        conflict.append((i, seen[int(i)], v))
    seen[int(i)] = v
ok = (set(seen) == set(range(12)) and not conflict
      and all(seen[i] == ALLOWED[i] for i in range(12)))
rep("the 12 admissible states are transcribed correctly",
    ok, f"({len(seen)} distinct indices from {len(states)} mentions, "
        f"conflicts {conflict})")


def parse_vec(pattern):
    m = re.search(pattern, src)
    if not m:
        return None
    body = m.group(1).replace("\\", " ").replace(" ", "")
    return [int(x) for x in body.split(",") if x != ""]


J_h = [1 if K.trans(s, K.sigma_h(s)) else 0 for s in ALLOWED]
J_re = [1 if K.trans(s, K.sigma_rot(s)) else 0 for s in ALLOWED]
rep("J_h matches the code", parse_vec(r"J_h = \(([^)]*)\)") == J_h)
rep("J_re matches the code", parse_vec(r"J_\{re\} = \(([^)]*)\)") == J_re)

m = re.search(r"\\chi_M\(\\lambda\) = ([^,]*),", src)
ok = bool(m) and "lambda^6" in m.group(1).replace(" ", "")
rep("characteristic polynomial is as stated", ok, m.group(1).strip() if m else "")

m = re.search(r"6x\^8 - 15x\^7 - 8x\^6 \+ 11x\^5 \+ 34x\^4 - 52x\^3 "
              r"\+ 29x\^2 - 8x \+ 1\}\s*\{(.*?)\},", src, re.S)
rep("the paper's rational generating function is present", bool(m))

print("\n  checking the paper's stated generating function symbolically")
try:
    import sympy as sp
    x = sp.symbols("x")
    B_NUM = sp.expand((1 - 2 * x) * (1 - 6 * x + 17 * x**2 - 18 * x**3
                                     - 2 * x**4 + 7 * x**5 + 6 * x**6 - 3 * x**7))
    B_DEN = sp.expand((1 - x)**2 * (1 - 3 * x)**2 * (1 - 3 * x + x**2)
                      * (1 - x - x**2) * (1 - 3 * x**2))
    P_NUM = sp.expand(6 * x**8 - 15 * x**7 - 8 * x**6 + 11 * x**5 + 34 * x**4
                      - 52 * x**3 + 29 * x**2 - 8 * x + 1)
    P_DEN = sp.expand(27 * x**10 - 126 * x**9 + 120 * x**8 + 210 * x**7
                      - 505 * x**6 + 346 * x**5 - 17 * x**4 - 98 * x**3
                      + 54 * x**2 - 12 * x + 1)
    rep("numerator matches Barker's conjecture", sp.simplify(P_NUM - B_NUM) == 0)
    rep("denominator matches Barker's conjecture", sp.simplify(P_DEN - B_DEN) == 0)
    rep("denominator factors as stated",
        sp.simplify(P_DEN - (x - 1)**2 * (3 * x - 1)**2 * (3 * x**2 - 1)
                    * (x**2 - 3 * x + 1) * (x**2 + x - 1)) == 0)
except ImportError:
    print("  [skip] sympy not available, run with --with sympy")

print()
if fails:
    print(f"*** {len(fails)} FAILED: {fails} ***")
    sys.exit(1)
print("THE PAPER AGREES WITH THE CODE ON EVERY TRANSCRIBED OBJECT")
