r"""Convert Markdown math between the two delimiter styles used in this repository, in both directions (V2.0.11).

  legacy (V2.0.01-V2.0.10 as first uploaded):  \( X \)  inline,   \[ X \]  display on its own line(s)
  github (V2.0.11 and later):                  $`X`$     inline,   ```math fenced block (X on the lines in between)

Only the delimiters change; the LaTeX between them is copied byte for byte.  The two conversions are exact inverses, so the original
upload can always be reconstructed from the current files:

  python scripts/convert_math.py to-legacy FILE      # print FILE converted to the legacy style
  python scripts/convert_math.py to-github FILE      # print FILE converted to the GitHub style
  python scripts/convert_math.py selftest            # round-trip checks on built-in samples (exit 0 = pass)

Example (verify that V2.0.11 equals the uploaded V2.0.10 apart from the title line):
  git show 25498f6:protocol/SND_protocol_V2.0.10.md | tail -n +2 > /tmp/old_body
  python scripts/convert_math.py to-legacy ../protocol/SND_protocol_V2.0.11.md | tail -n +2 | cmp - /tmp/old_body && echo identical
(`print` adds no trailing newline here, so the bytes match.)
"""
import re, sys

INL_OLD = re.compile(r"\\\((.+?)\\\)", re.S)
DSP_OLD = re.compile(r"(?m)^\\\[(.+?)\\\]$", re.S)
INL_NEW = re.compile(r"\$`(.+?)`\$")
DSP_NEW = re.compile(r"(?m)^```math\n(.+?)\n```$", re.S)


def to_github(s):
    def inl(m):
        x = m.group(1)
        assert x == x.strip() and "`" not in x and "$" not in x and "\n" not in x, "unsafe inline formula: %r" % x
        return "$`" + x + "`$"

    def dsp(m):
        x = m.group(1)
        assert "`" not in x and "$" not in x, "unsafe display formula: %r" % x
        return "```math\n" + x + "\n```"
    return INL_OLD.sub(inl, DSP_OLD.sub(dsp, s))      # display first, so a display block is not caught as inline


def to_legacy(s):
    s = INL_NEW.sub(lambda m: "\\(" + m.group(1) + "\\)", s)
    return DSP_NEW.sub(lambda m: "\\[" + m.group(1) + "\\]", s)


def selftest():
    samples = [
        r"全體：\(U=\sum_bU_b\)，\(V=\sum_bV_b\)。" + "\n",
        "x\n\\[\\beta=\\omega\\sqrt{\\tfrac{2(L-1)}{L}}\\,\\gamma .\\]\ny\n",
        r"\(\{a\}\,b\\c\)" + "\n",
        "no math here, price 5 and a (paren)\n",
    ]
    for s in samples:
        g = to_github(s)
        assert "\\(" not in g and "\\[" not in g, g
        assert to_legacy(g) == s, (s, g)
    for bad in (r"\( x \)", "\\(a`b\\)", r"\(a$b\)"):
        try:
            to_github(bad)
        except AssertionError:
            continue
        raise SystemExit("selftest FAIL: unsafe input accepted: %r" % bad)
    print("selftest OK: %d round-trip samples, 3 unsafe inputs rejected" % len(samples))


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "selftest":
        selftest()
    elif len(sys.argv) == 3 and sys.argv[1] in ("to-github", "to-legacy"):
        txt = open(sys.argv[2], encoding="utf8", newline="").read()
        sys.stdout.write((to_github if sys.argv[1] == "to-github" else to_legacy)(txt))
    else:
        raise SystemExit(__doc__)
