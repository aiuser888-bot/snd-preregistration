r"""Rewrite two formula patterns that github.com (MathJax) typesets badly (V2.0.12).  The mathematical meaning does not change.

  1. \bar X              ->  \overline{X}        MathJax draws \bar as a tiny dash floating at the top-left of the letter; \overline draws a full bar.
  2. accent-atom ^ sup   ->  (accent-atom)^sup   e.g. \widetilde f_{b,i}^{\,2}  ->  (\widetilde f_{b,i})^{2}
                                                 MathJax lifts the superscript above the accent; with parentheses it sits at the normal height.
     accents: \widehat \widetilde \hat \tilde \overline.  The subscript, if any, stays inside the parentheses.  A thin space \, that
     only separated the superscript from the accent, written as ^{\,2}, is dropped (the parenthesis does that job).

Only text inside math spans ($`...`$ and ```math blocks) is touched; everything else is copied byte for byte.

  python scripts/render_fix.py selftest        # built-in checks (exit 0 = pass)
  python scripts/render_fix.py apply FILE      # print FILE with the rewrites applied
  python scripts/render_fix.py diff FILE       # list every formula that would change:  line: old  =>  new
check_math.py uses fix_tex(): a formula that fix_tex() would still change counts as a rendering hazard.
"""
import re, sys

ACCENTS = {"\\widehat", "\\widetilde", "\\hat", "\\tilde", "\\overline"}
TOK = re.compile(r"\\[A-Za-z]+|\\.|\s+|.", re.S)
INLINE = re.compile(r"\$`(.+?)`\$")
DBL_CODE = re.compile(r"``(.+?)``")


def _tokens(s):
    return TOK.findall(s)


def _skip_ws(t, i):
    while i < len(t) and t[i].isspace():
        i += 1
    return i


def _arg_end(t, i):
    """t[i] starts a TeX argument (whitespace already skipped); return the index just after it."""
    if i >= len(t):
        raise ValueError("missing argument")
    if t[i] == "{":
        d = 0
        for j in range(i, len(t)):
            if t[j] == "{":
                d += 1
            elif t[j] == "}":
                d -= 1
                if d == 0:
                    return j + 1
        raise ValueError("unbalanced braces")
    return i + 1


def _fix_bar_once(s):
    t = _tokens(s)
    for i, tok in enumerate(t):
        if tok == "\\bar":
            j = _skip_ws(t, i + 1)
            k = _arg_end(t, j)
            arg = t[j:k]
            inner = arg[1:-1] if arg[0] == "{" else arg
            return "".join(t[:i]) + "\\overline{" + "".join(inner) + "}" + "".join(t[k:]), True
    return s, False


def _fix_accent_once(s):
    t = _tokens(s)
    for i, tok in enumerate(t):
        if tok not in ACCENTS:
            continue
        j = _skip_ws(t, i + 1)
        k = _arg_end(t, j)                         # end of the accent's base
        sub = sup = None                          # (start, end) token ranges of the script arguments
        m = k
        for _ in range(2):
            n = _skip_ws(t, m)
            if n < len(t) and t[n] == "_" and sub is None:
                a = _skip_ws(t, n + 1); b = _arg_end(t, a); sub = (n, b); m = b
            elif n < len(t) and t[n] == "^" and sup is None:
                a = _skip_ws(t, n + 1); b = _arg_end(t, a); sup = (n, b); m = b
            else:
                break
        if sup is None:
            continue
        atom = "".join(t[i:k]) + ("".join(t[sub[0]:sub[1]]) if sub else "")
        supx = "".join(t[sup[0]:sup[1]])
        supx = re.sub(r"^\^\s*\{\s*\\,\s*", "^{", supx)          # ^{\,2} -> ^{2}
        return "".join(t[:i]) + "(" + atom + ")" + supx + "".join(t[m:]), True
    return s, False


def fix_tex(tex):
    changed = True
    while changed:
        tex, changed = _fix_bar_once(tex)
    changed = True
    while changed:
        tex, changed = _fix_accent_once(tex)
    return tex


def fix_markdown(text):
    """-> (new_text, [(line_number, old_tex, new_tex)]).  Fenced code blocks and double-backtick literal examples are not touched."""
    assert "\r" not in text, "CRLF file: normalise first"
    out, changes = [], []
    fence, buf, start = None, [], 0
    for n, line in enumerate(text.split("\n"), 1):
        s = line.strip()
        if fence is None and s.startswith("```"):
            fence = "math" if s[3:].strip() == "math" else "code"
            buf, start = [], n
            out.append(line)
            continue
        if fence is not None:
            if s.startswith("```") and s.strip("`") == "":
                if fence == "math":
                    old = "\n".join(buf)
                    new = fix_tex(old)
                    if new != old:
                        changes.append((start, old, new))
                    out.extend(new.split("\n"))
                else:
                    out.extend(buf)
                out.append(line)
                fence = None
            else:
                buf.append(line)
            continue
        keep = []
        def hold(m):
            keep.append(m.group(0))
            return "\x00%d\x00" % (len(keep) - 1)
        masked = DBL_CODE.sub(hold, line)

        def inl(m):
            new = fix_tex(m.group(1))
            if new != m.group(1):
                changes.append((n, m.group(1), new))
            return "$`" + new + "`$"
        masked = INLINE.sub(inl, masked)
        out.append(re.sub(r"\x00(\d+)\x00", lambda m: keep[int(m.group(1))], masked))
    if fence is not None:
        raise ValueError("unclosed fenced block starting at line %d" % start)
    return "\n".join(out), changes


def selftest():
    cases = [
        (r"\bar f_b=L^{-1}\sum_if_{b,i}", r"\overline{f}_b=L^{-1}\sum_if_{b,i}"),
        (r"\bar{\widetilde f}_b", r"\overline{\widetilde f}_b"),
        (r"(\mathrm{Re}\,\bar wU)^2", r"(\mathrm{Re}\,\overline{w}U)^2"),
        (r"\bar q_{c(i)}", r"\overline{q}_{c(i)}"),
        (r"\widetilde f_{b,i}^{\,2}", r"(\widetilde f_{b,i})^{2}"),
        (r"\sigma_f^2=\int\widehat p_0\widehat f^2", r"\sigma_f^2=\int\widehat p_0(\widehat f)^2"),
        (r"\widehat\sigma_{\rm rep}^2+\widehat\sigma_{\rm orig}^2", r"(\widehat\sigma_{\rm rep})^2+(\widehat\sigma_{\rm orig})^2"),
        (r"\bar w^2U^2", r"(\overline{w})^2U^2"),
        (r"\widehat f^2_j", r"(\widehat f_j)^2"),
        (r"\sum_j\widetilde f_j^2\bigr]", r"\sum_j(\widetilde f_j)^2\bigr]"),
        (r"\widehat x^2\sim", r"(\widehat x)^2\sim"),
        # unchanged
        (r"\widehat g_j=K^{-1}\sum_b(f_{b,j}-\overline{f}_b)", None),
        (r"\widehat{\sigma^2}", None),
        (r"(\widehat f)^2", None),
        (r"\overline{AB}", None),
        (r"x\\bar y", None),                       # \\ then the letters "bar", not the macro
        (r"\barwedge", None),
        (r"\widehat p_0\widehat f", None),
        (r"\widehat\sigma_f", None),
        (r"\int\widehat p_0\widehat f=0", None),
    ]
    for old, new in cases:
        got = fix_tex(old)
        exp = old if new is None else new
        assert got == exp, (old, got, exp)
        assert fix_tex(got) == got, ("not idempotent", old, got)
    md = "a ``$`\\bar x`$`` b $`\\bar x^2`$\n```python\n$`\\bar y`$\n```\n```math\n\\bar z\n```\nend $`\\hat u`$"
    got, ch = fix_markdown(md)
    assert got == "a ``$`\\bar x`$`` b $`(\\overline{x})^2`$\n```python\n$`\\bar y`$\n```\n```math\n\\overline{z}\n```\nend $`\\hat u`$", got
    assert len(ch) == 2, ch
    print("selftest OK: %d formula cases (rewrites, unchanged, idempotence) + markdown masking" % len(cases))


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "selftest":
        selftest()
    elif len(sys.argv) == 3 and sys.argv[1] in ("apply", "diff"):
        txt = open(sys.argv[2], encoding="utf8", newline="").read()
        new, ch = fix_markdown(txt)
        if sys.argv[1] == "apply":
            sys.stdout.write(new)
        else:
            for ln, o, n in ch:
                print("%d: %s  =>  %s" % (ln, o.replace("\n", " "), n.replace("\n", " ")))
            print("%d formula(s) would change" % len(ch))
    else:
        raise SystemExit(__doc__)
