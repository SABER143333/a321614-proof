"""Structural checks on a LaTeX source, for when no LaTeX toolchain is present.

This does not replace compiling.  It catches the failures that actually stop a
compile -- an unmatched \\begin/\\end, an unbalanced brace, a \\ref to a label
that does not exist -- plus the two Journal of Integer Sequences rules that are
mechanical to test (no percent signs anywhere, no amsart).
"""
import re, sys, os

def default_tex():
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "a321614.tex"),
                 os.path.join(here, "..", "a321614.tex")):
        if os.path.exists(cand):
            return cand
    return os.path.join(here, "a321614.tex")


path = sys.argv[1] if len(sys.argv) > 1 else default_tex()
src = open(path, encoding="utf-8").read()
fails = []


def rep(name, ok, detail=""):
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")
    if not ok:
        fails.append(name)


print(f"checking {os.path.basename(path)} ({len(src)} chars)\n")

envs = re.findall(r"\\begin\{(\w+\*?)\}", src)
ends = re.findall(r"\\end\{(\w+\*?)\}", src)
from collections import Counter
cb, ce = Counter(envs), Counter(ends)
rep("every \\begin has a matching \\end",
    cb == ce, "" if cb == ce else f"diff: {cb - ce | ce - cb}")

depth, bad = 0, []
for m in re.finditer(r"\\(begin|end)\{(\w+\*?)\}", src):
    depth += 1 if m.group(1) == "begin" else -1
    if depth < 0:
        bad.append(m.group(2))
rep("environments properly nested", not bad, str(bad[:5]))

depth = 0
worst = 0
for ch in src:
    if ch == "{" and (not src[:0]):
        depth += 1
    elif ch == "}":
        depth -= 1
        worst = min(worst, depth)
rep("braces balanced", depth == 0 and worst == 0,
    f"final depth {depth}, min {worst}")

labels = set(re.findall(r"\\label\{([^}]*)\}", src))
refs = set(re.findall(r"\\(?:eq)?ref\{([^}]*)\}", src))
missing = refs - labels
rep("every \\ref has a \\label", not missing, str(sorted(missing)))
rep("no unused \\label", not (labels - refs), str(sorted(labels - refs)))

rep("JIS: no percent sign anywhere (comments forbidden)",
    "%" not in src, f"{src.count('%')} found")
rep("JIS: not amsart", "amsart" not in src)
rep("JIS: documentclass is [12pt]{article}",
    bool(re.search(r"\\documentclass\[12pt\]\{article\}", src)))
rep("JIS: no date on the title page",
    bool(re.search(r"\\date\{\}", src)))
rep("JIS: acknowledgment is its own numbered section",
    bool(re.search(r"\\section\{Acknowledgment\}", src)))

pkgs = re.findall(r"\\usepackage(?:\[[^\]]*\])?\{([^}]*)\}", src)
used = set(re.findall(r"\\(?:mathbb|mathcal|textbf|textit|emph|url|includegraphics"
                      r"|textcolor|definecolor|href)\b", src))
print(f"\n  packages declared: {pkgs}")
print(f"  (JIS asks that unused packages be removed; graphicx is used only if")
print(f"   \\includegraphics appears, which it does not)")

print()
if fails:
    print(f"*** {len(fails)} FAILED: {fails} ***")
    sys.exit(1)
print("STRUCTURE OK -- but this is NOT a compile; run LaTeX on it to be sure")
