#!/usr/bin/env python
"""Build every Abstract that review/POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md section 7.4
can produce, from the Abstract of C1_LUAD_ADDENDA.md Entry 12, and count each one under the
four conventions of Entries 11 and 12.  Reads no data.  Output: ABSTRACT_VARIANTS.md beside
this file (deterministic; the same input gives the same bytes).

Conventions (Entry 12, "Deviation"): tokens = whitespace-separated tokens after removing TeX
markup; ranges = the same with number ranges split; typeset = tokens with every "P = x"
counted as three words; typeset+ranges = both.
"""
import hashlib
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "QUALITY_REVIEW_2026-11-20", "abstract_v12", "ABSTRACT_FINAL.tex")
SRC_SHA = "6a1bfede2f9a56dd811cf46b345ecd2989f7ffb5697b1578e574fa8fec298ddc"
BS = chr(92)

NAME = {"C": "composition scores", "R": "expression components"}
ADJ = {"u": "uterine", "l": "lung"}
NOUN = {"u": "uterus", "l": "lung"}

THESIS_OLD = "Here, we show that histology carries a set-level signal of which proteins depart from their own transcript."
THESIS_CASE2 = "We test whether histology signals which proteins depart from their own transcript."
# ordered cuts; (label, old, new, cases in which the cut exists)
CUTS = [
    ("c1", "At these cohort sizes, morphology adds", "Morphology adds", (1, 2)),
    ("c2", "Here, we show that histology carries", "Histology carries", (1,)),
    ("c3", "Elsewhere, four protein-array", "Four protein-array", (1, 2)),
    ("c5", "2{,}310 genes)", "2{,}310)", (1, 2)),
    ("c6", "(138 patients)", "(138)", (1, 2)),
]
ANCHOR = "Elsewhere, four protein-array"        # the qualifier sentence goes immediately before it


def counts(tex):
    s = tex.replace(BS + "noindent", " ")
    n_p = len(re.findall(r"\$P=[0-9.]+\$", s))
    s = re.sub(r"\$([^$]*)\$", r"\1", s)
    plain = re.sub(r"[{}" + BS + BS + r"]", "", s)
    tok = len(plain.split())
    rng = len(plain.replace("--", " ").replace("/", " ").split())
    return tok, rng, tok + 2 * n_p, rng + 2 * n_p


def qualifier(state):
    """state: dict cohort -> 'C' | 'R' (not shown, with the covariate set that set the flag)
    or 'NC:C' | 'NC:R' (that reading was closed before freezing and the test was not run)."""
    shown = {c: s for c, s in state.items() if not s.startswith("NC")}
    notc = {c: s[3:] for c, s in state.items() if s.startswith("NC")}
    if len(shown) == 2:
        if shown["u"] == shown["l"]:
            return f"Post hoc, adjusted for {NAME[shown['u']]}, neither uncorrected set did."
        return (f"Post hoc, adjusted (uterus, {NAME[shown['u']]}; lung, {NAME[shown['l']]}), "
                f"neither uncorrected set did.")
    if len(notc) == 2:
        if notc["u"] == notc["l"]:
            return f"Post hoc, the tests of the uncorrected sets adjusted for {NAME[notc['u']]} were not run."
        return (f"Post hoc, the tests of the uncorrected sets after adjustment (uterus, {NAME[notc['u']]}; "
                f"lung, {NAME[notc['l']]}) were not run.")
    if len(shown) == 1 and not notc:
        (c, s), = shown.items()
        return f"Post hoc, adjusted for {NAME[s]}, the uncorrected {ADJ[c]} set did not."
    if len(notc) == 1 and not shown:
        (c, s), = notc.items()
        return f"Post hoc, the test of the uncorrected {ADJ[c]} set adjusted for {NAME[s]} was not run."
    (c, s), = shown.items()
    (d, t), = notc.items()
    if s == t:
        return (f"Post hoc, adjusted for {NAME[s]}, the uncorrected {ADJ[c]} set did not; that test of the "
                f"uncorrected {ADJ[d]} set was not run.")
    return (f"Post hoc, adjusted for {NAME[s]}, the uncorrected {ADJ[c]} set did not; the test of the "
            f"uncorrected {ADJ[d]} set adjusted for {NAME[t]} was not run.")


def build(base, state):
    case = len(state)
    text = base
    assert text.count(ANCHOR) == 1
    text = text.replace(ANCHOR, qualifier(state) + " " + ANCHOR)
    if case == 2:
        assert text.count(THESIS_OLD) == 1
        text = text.replace(THESIS_OLD, THESIS_CASE2)
    target = counts(base)
    used = []
    for label, old, new, cases in CUTS:
        if all(a <= b for a, b in zip(counts(text), target)):
            break
        if case not in cases:
            continue
        assert text.count(old) == 1, (label, old)
        text = text.replace(old, new)
        used.append(label)
    return text, used


def main():
    raw = open(SRC, "rb").read()
    assert hashlib.sha256(raw).hexdigest() == SRC_SHA
    base = " ".join(raw.decode("utf-8").replace(BS + "noindent", "").split())
    assert base.count(THESIS_OLD) == 1
    b = counts(base)
    kinds = ("", "C", "R", "NC:C", "NC:R")
    states = []
    for ku in kinds:
        for kl in kinds:
            st = {}
            if ku:
                st["u"] = ku
            if kl:
                st["l"] = kl
            if st:
                states.append((f"{len(st)}-u[{ku or '-'}]-l[{kl or '-'}]", st))
    L = ["# Abstracts fixed by POSTHOC_C1_COVARIATES_PRESPEC_2026-11-21.md section 7.4", "",
         "Generated by `make_abstract_variants.py` from the Entry 12 Abstract (TeX source sha256",
         f"`{SRC_SHA}`). Nothing here depends on a result. Counts are tokens / ranges split /",
         "typeset / typeset with ranges split, as in Entries 11 and 12.", "",
         f"## Case 0 (no cohort flagged): the Entry 12 Abstract, unchanged ({b[0]} / {b[1]} / {b[2]} / {b[3]})", "",
         "```tex", base, "```", ""]
    for name, st in states:
        text, used = build(base, st)
        c = counts(text)
        over = max(x - y for x, y in zip(c, b))
        desc = ", ".join(f"{NOUN[k]}: type {v}" for k, v in st.items())
        L += [f"## Variant {name} ({desc})", "",
              f"Counts {c[0]} / {c[1]} / {c[2]} / {c[3]}; cuts applied: {', '.join(used) if used else 'none'}"
              + ("; thesis sentence replaced" if len(st) == 2 else "")
              + (f"; **longer than the Entry 12 Abstract by {over}**" if over > 0 else "") + ".", "",
              "Qualifier: " + qualifier(st), "", "```tex", text, "```", ""]
    out = os.path.join(HERE, "ABSTRACT_VARIANTS.md")
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L))
    print("written", out)
    for name, st in states:
        text, used = build(base, st)
        print(name, counts(text), used, "| Q words:", len(qualifier(st).split()))
    print("base", b)


if __name__ == "__main__":
    main()
