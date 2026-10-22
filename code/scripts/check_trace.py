"""Traceability between code and protocol: every protocol reference in the code exists in the protocol, every public function of snd/ cites one,
the constants shared by code and protocol agree, and docs/TRACEABILITY.md is the table generated from the docstrings.

Implements: protocol 3.1, 3.2, 4, 5.1, 7(D) (the constants compared below).  Checks that the code's references to sections, lemmas and sub-items are real, and that the
constants that appear in both (M_perm, level, tolerances, key size, table grid, check points, Philox variant) are equal.
Verified by: tests/test_snd.py::test_trace_labels_parse

Docstring convention (snd/*.py, scripts/*.py).  A line "Implements: ..." in the docstring of a module, class or function lists what it implements:
  protocol 3.1 | protocol 7(D) | protocol 12.5.1 step 3 | Lemma 5(c) | Lemma 2 (iv)       (sub-items: letters in 7, lemma items a, b, d', i, ii, ...)
and an optional line "Verified by: tests/test_snd.py::name, name, ..." names the tests.  Section numbers are the numbers of the headings of the
protocol; lemma numbers and items are those printed in the protocol text.

The same run checks that the version label (snd.__version__) and the number of tests agree in every file that states them (CITATION.cff,
the Zenodo files, the protocol file name and title, VERSIONING.md, README.md, CHANGELOG.md and its copy, the evolution documents, the audit file).

Usage:  python scripts/check_trace.py [--protocol PATH] [--write]      check; --write regenerates docs/TRACEABILITY.md
        python scripts/check_trace.py --check                         also fail if docs/TRACEABILITY.md is not the generated table
"""
import ast, glob, importlib, inspect, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.normpath(os.path.join(HERE, ".."))
ROOT = os.path.normpath(os.path.join(CODE, ".."))
sys.path.insert(0, CODE)
PRIME = "\u2032"

SEC_RE = re.compile(r"\bprotocol\s+(\d+(?:\.\d+){0,2})(?:\(([A-Za-z][\u2032']?)\))?(?:\s+steps?\s+(\d))?")
LEM_RE = re.compile(r"\bLemma\s+(\d)((?:\s*\(([a-z]+[\u2032']?)\))*)")
ITEM_RE = re.compile(r"\(([a-z]+[\u2032']?)\)")


def norm(x):
    return x.replace("'", PRIME)


def read_protocol(path):
    lines = open(path, encoding="utf8").read().split("\n")
    secs, lemma_items, steps_1251, seven = set(), {}, 0, set()
    lemma_start = re.compile(r"^(?:#{2,3} \d+(?:\.\d+)* |\*\*(?:\d+\.\d+\.\d+ )?)引理 (\d)")
    heading = re.compile(r"^(?:#{2,3} (\d+(?:\.\d+)*)\b|\*\*(12\.5\.\d)\b)")
    cur_lemma, cur_sec = None, None
    for ln in lines:
        h = heading.match(ln)
        if h:
            cur_sec = h.group(1) or h.group(2)
            secs.add(cur_sec)
            if cur_lemma is not None:
                cur_lemma = None
        lm = lemma_start.match(ln)
        if lm:
            cur_lemma = int(lm.group(1))
            lemma_items.setdefault(cur_lemma, set())
        if cur_lemma is not None:
            it = re.match(r"^\(([a-z]+[\u2032']?)\)", ln)
            if it:
                lemma_items[cur_lemma].add(norm(it.group(1)))
        if cur_sec == "7":
            sv = re.match(r"^\*\*\(([A-F][\u2032']?)\)", ln)
            if sv:
                seven.add(norm(sv.group(1)))
        if cur_sec == "12.5.1":
            if re.match(r"^\d\. \*\*", ln):
                steps_1251 += 1
    return dict(secs=secs, lemma_items=lemma_items, steps_1251=steps_1251, seven=seven)


def docstrings(path, rel):
    tree = ast.parse(open(path, encoding="utf8").read())
    out = []                                           # (kind, qualified name, lineno, docstring)
    out.append(("module", rel, 1, ast.get_docstring(tree) or ""))
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            out.append(("class" if isinstance(node, ast.ClassDef) else "function", node.name, node.lineno, ast.get_docstring(node) or ""))
    return out


def tag(doc, key):
    m = re.search(r"^%s:\s*(.*?)(?=\n\s*\n|\n[A-Z][A-Za-z ]*:|\Z)" % key, doc, re.S | re.M)
    return " ".join(m.group(1).split()) if m else None


def parse_refs(text):
    """-> list of ('sec', number, item or None, step or None) / ('lem', n, [items])"""
    refs = []
    for m in SEC_RE.finditer(text):
        refs.append(("sec", m.group(1), norm(m.group(2)) if m.group(2) else None, int(m.group(3)) if m.group(3) else None))
    for m in LEM_RE.finditer(text):
        refs.append(("lem", int(m.group(1)), [norm(i) for i in ITEM_RE.findall(m.group(2) or "")], None))
    return refs


def validate(refs, P):
    bad = []
    for r in refs:
        if r[0] == "sec":
            _, sec, item, step = r
            if sec not in P["secs"]:
                bad.append("section %s does not exist" % sec)
            elif sec == "7" and item and item not in P["seven"]:
                bad.append("7(%s) does not exist" % item)
            elif sec == "12.5.1" and step and not (1 <= step <= P["steps_1251"]):
                bad.append("12.5.1 step %d does not exist (has %d)" % (step, P["steps_1251"]))
        else:
            _, n, items, _ = r
            if n not in P["lemma_items"]:
                bad.append("Lemma %d does not exist" % n)
            else:
                for i in items:
                    if i not in P["lemma_items"][n]:
                        bad.append("Lemma %d(%s) does not exist" % (n, i))
    return bad


def collect():
    rows = []
    files = sorted(glob.glob(os.path.join(CODE, "snd", "*.py"))) + sorted(glob.glob(os.path.join(CODE, "scripts", "*.py")))
    for f in files:
        rel = os.path.relpath(f, CODE).replace(os.sep, "/")
        for kind, name, line, doc in docstrings(f, rel):
            rows.append(dict(file=rel, kind=kind, name=name, line=line, impl=tag(doc, "Implements"), tests=tag(doc, "Verified by")))
    return rows


def constants_checks(ptxt):
    """Constants that appear both in the protocol text and in the code; returns list of (ok, description)."""
    import numpy as np
    from snd import design, fcrice, numerics, randtest
    out = []
    add = lambda ok, d: out.append((bool(ok), d))
    sig = inspect.signature(randtest.randomization_test).parameters
    add(re.search(r"M_\{\\rm perm\}=1999", ptxt) and sig["M"].default == 1999, "M_perm = 1999 (protocol 4) = randomization_test default")
    add(re.search(r"\\widehat p\\le0\.05", ptxt) and randtest.ALPHA == 0.05, "reject iff p_hat <= 0.05 (protocol 4) = randtest.ALPHA")
    add(re.search(r"10\^\{-9\}L", ptxt) and inspect.signature(design.design_check).parameters["tol"].default == 1e-9, "design check tolerance 1e-9 L (protocol 3.2)")
    add("≥512 位元" in ptxt and design.MIN_RAW_BYTES * 8 == 512, "raw QRNG output >= 512 bits (protocol 3.1, 12.5.1) = MIN_RAW_BYTES")
    add("前 128 位元" in ptxt and len(design.seed_to_key(bytes(64))) * 64 == 128, "key = first 128 bits of SHA-256 (protocol 3.1) = 2 x uint64")
    add(re.search(r"\\Delta U\|/\\sqrt V\\le10\^\{-6\}", ptxt) and randtest.TOL_U == 1e-6, "|dU|/sqrt(V) <= 1e-6 (protocol 7(D)) = randtest.TOL_U")
    add(re.search(r"\\Delta V\|/V\\le10\^\{-9\}", ptxt) and numerics.sum_check.__doc__ and "dV_rel <= 1e-9" in numerics.sum_check.__doc__, "|dV|/V <= 1e-9 (protocol 7(D)) documented in numerics.sum_check")
    xs = np.array([1.0, 20.0, 25.0]); ucb = np.array([2.0, 3.0, 4.0])
    add("x+1.96" in ptxt and fcrice.XMAX_VALID == 18.0 and abs(float(fcrice.u_of_x(25.0, xs, ucb)) - 26.96) < 1e-12 and "\\widehat x\\le18" in ptxt,
        "u(x) = x + 1.96 beyond x_hat > 18 (protocol 5.1 step 4) = fcrice.u_of_x")
    p = inspect.signature(fcrice.build_table).parameters
    xs0, _ = fcrice.load_table(os.path.join(ROOT, "results", "fc_rice_table.npz"))
    add("x=0.004,0.008,\\dots,21.996" in ptxt.replace("`", "").replace("$", "") and "共 5499 點" in ptxt and "共 4001 點" in ptxt and
        (p["dx"].default, p["xmax"].default, p["dn"].default, p["numax"].default) == (0.004, 22.0, 0.005, 20.0) and
        (len(xs0), round(float(xs0[0]), 6), round(float(xs0[-1]), 6)) == (5499, 0.004, 21.996) and
        len(np.arange(0, p["numax"].default + p["dn"].default / 2, p["dn"].default)) == 4001, "table grid (protocol 5.1 step 4) = build_table defaults = stored table")
    cps = {float(a): float(b) for a, b in re.findall(r"u\((\d+(?:\.\d+)?)\)=(\d+\.\d+)", ptxt.split("檢核點")[1].split("\n")[0])}
    add(cps == fcrice.TABLE_CHECKPOINTS, "table check points (protocol 5.1 step 5) = fcrice.TABLE_CHECKPOINTS")
    add("Philox4x64-10" in ptxt and "4x32" not in ptxt and "Threefry" not in ptxt, "protocol names exactly Philox4x64-10 (the code's generator)")
    add(re.search(r"\(m\s*,\s*b\s*,\s*0\s*,\s*0\)", ptxt.replace("`", "").replace("$", "")) is not None, "counter layout (m, b, 0, 0) stated in the protocol (3.1)")
    out += second_tier_checks(ptxt)
    out += number_checks(ptxt, ROOT)
    return out


def sci(x, digits=1):
    """Mantissa and exponent as the protocol writes them: 7.35e-16 -> '7.4\\times10^{-16}'."""
    m, e = ("%.*e" % (digits, x)).split("e")
    return "%s\\times10^{%d}" % (m, int(e))


def second_tier_checks(ptxt):
    """Constants of the second tier (protocol 4, 7(D), 10) against highprec.py, decide.py, randtest.py and numerics.py.  Returns [(ok, description)]."""
    import inspect, math
    from decimal import Decimal, getcontext
    from snd import decide, highprec, numerics, randtest
    out = []
    add = lambda ok, d: out.append((bool(ok), d))
    flat = ptxt.replace("`", "").replace("$", "")
    add(re.search(r"\\delta=2\\times10\^\{-6\}\\sqrt V", ptxt) and "2.0 * tol_u * np.sqrt(V)" in inspect.getsource(randtest.randomization_test) and randtest.TOL_U == 1e-6,
        "near-tie width delta = 2e-6 sqrt(V) (protocol 4) = randomization_test")
    add("P=128" in flat and "P=1024" in flat and "128 \u4f4d\u5143" in flat and "1024 \u4f4d\u5143" in flat and highprec.P_START == 128 and highprec.P_CAP == 1024,
        "precision schedule 128 ... 1024 bits (protocol 4, 7(D)) = highprec.P_START, P_CAP")
    add(numerics.LIBM_ABS_TOL == 2.0 ** -50 and "2^{-50}" in ptxt and len(numerics._libm_sample()) == 3088 and "3088 \u500b\u89d2\u5ea6" in ptxt,
        "library qualification: absolute error 2^-50 and a sample of 3088 angles (protocol 7(D)) = numerics.LIBM_ABS_TOL, _libm_sample")
    getcontext().prec = 60
    two_ln20 = 2 * Decimal(20).ln()
    lo, hi = (Decimal(q.numerator) / Decimal(q.denominator) for q in (decide.CHI2_LO, decide.CHI2_HI))
    add(lo <= two_ln20 <= hi and hi - lo < Decimal("1e-23") and abs(Decimal(decide.CHI2_FLOAT) - two_ln20) < Decimal("1e-15") and "2\\ln20=5.99146" in ptxt,
        "chi-square quantile 2 ln 20 = 5.99146... (protocol 10 (ii)) = decide.CHI2_LO / CHI2_HI / CHI2_FLOAT")
    add(decide.R_MIN == 1e-4 and decide.MARGIN == 1e-12 and "10^{-12}" in ptxt and "UCB}_{95}<10^{-4}" in ptxt,
        "r_min = 1e-4 and the relative margin 1e-12 (protocol 10) = decide.R_MIN, MARGIN")
    add("d = 2.0 * tol_u" in inspect.getsource(decide.exclusion) and "max(1e-4, 100 * tol_u)" in inspect.getsource(decide.reproduction_ii) and
        "2\\times10^{-6}" in ptxt.split("\u9580\u6abb\u6bd4\u8f03\u7684\u6578\u503c\u8a8d\u8b49")[1] and "10^{-4}" in ptxt.split("\u9580\u6abb\u6bd4\u8f03\u7684\u6578\u503c\u8a8d\u8b49")[1],
        "windows that trigger the certified comparison (protocol 10: 2e-6 in x_hat, 1e-4 in the chi-square statistic) = decide.exclusion, reproduction_ii")
    xq = math.sqrt(2 * math.log(20))
    f = lambda x: x * math.exp(-x * x / 2)
    reg = 2.0 / 2000 * (1 - math.exp(-2000 * f(xq) * 2e-6))
    reg_mean = 1999 * 2 * 2e-6 * math.sqrt(math.pi) / 4
    add(sci(reg) in ptxt and sci(reg_mean) in ptxt and ("\\approx%d" % round(2000 * f(xq))) in ptxt,
        "trigger probability %s, mean number of near ties %s, lambda %d (protocol 7(D)) = the formulas evaluated at M_perm = 1999, delta_x = 2e-6"
        % (sci(reg), sci(reg_mean), round(2000 * f(xq))))
    from snd import design
    tb = design.theta_rounding_bound()
    eff = math.ceil(tb * 6 * 1.2e10 / math.sqrt(1.5e9) * 1e10) / 1e10                    # effect on U in units of sqrt(V) at sum|f~| = 6 N, V = 1.5e9, rounded up
    add(sci(tb, 2) in ptxt and ("\\lesssim%s\\sqrt V" % sci(eff, 1)) in ptxt,
        "binary64 angle error bound %s rad and its effect %s sqrt(V) at the registered scale (protocol 3.1 item 7) = design.theta_rounding_bound" % (sci(tb, 2), sci(eff, 1)))
    missing = []
    for m in sorted(set(re.findall(r"`(numerics|highprec|decide|randtest|design|stats)\.(\w+)`", ptxt))):
        mod = __import__("snd." + m[0], fromlist=[m[1]])
        if not hasattr(mod, m[1]):
            missing.append("%s.%s" % m)
    add(not missing, "code identifiers named in the protocol exist%s" % ((": missing " + ", ".join(missing)) if missing else ""))
    return out


def number_checks(ptxt, root):
    """Measured and bounded numbers quoted in protocol 7(D) against results/numerics_check.json.  Returns [(ok, description)]."""
    import json
    out = []
    add = lambda ok, d: out.append((bool(ok), d))
    path = os.path.join(root, "results", "numerics_check.json")
    if not os.path.exists(path):
        return [(False, "results/numerics_check.json exists")]
    d = json.load(open(path))
    rows, wc, tm = d["rows"], d["worst_case_bound"], d["timing"]
    q = lambda key: max(r[key] for r in rows)
    add(sci(q("dU_rel_sum")) in ptxt and sci(q("dV_rel_sum")) in ptxt and sci(q("dU_rel_vs_reference")) in ptxt and max(r["N"] for r in rows) == 10 ** 6,
        "reduced-scale summation differences %s (U), %s (V), %s (U vs reference) (protocol 7(D)) = results/numerics_check.json" % (sci(q("dU_rel_sum")), sci(q("dV_rel_sum")), sci(q("dU_rel_vs_reference"))))
    add(sci(wc["ulp4"]) in ptxt and sci(wc["tree_sum"]) in ptxt and wc["N"] == 1.2e10 and wc["V"] == 1.5e9 and wc["tree_depth"] == 35 and wc["formula_total"] < 1e-6,
        "registered-scale bounds %s (cos/sin) and %s (tree, depth 35) (protocol 7(D)) = results/numerics_check.json" % (sci(wc["ulp4"]), sci(wc["tree_sum"])))
    us, days = tm["reference_U_us_per_term"], tm["extrapolated_core_days_per_evaluation_N_1p2e10"]
    m1, m2 = re.search(r"\u6bcf\u9805\u7d04 (\d+)\u2013(\d+) \u5fae\u79d2", ptxt), re.search(r"\u7d04 (\d\.\d)\u2013(\d\.\d) \u6838\u00b7\u65e5", ptxt.split("\u7b97\u529b\u3002")[1])
    ok1 = m1 and int(m1.group(1)) / 1.5 <= us <= int(m1.group(2)) * 1.5
    ok2 = m2 and float(m2.group(1)) / 1.5 <= days <= float(m2.group(2)) * 1.5
    add(ok1 and ok2 and abs(days - us * 1e-6 * 1.2e10 / 86400) < 1e-9 * days,
        "cost %s us per term, %s core-days per evaluation (protocol 7(D)) agree within a factor 1.5 with the last measurement (%.1f us, %.2f core-days) in results/numerics_check.json (timing depends on the machine and its load)"
        % (m1 and "%s-%s" % m1.groups(), m2 and "%s-%s" % m2.groups(), us, days))
    return out


def internal_ref_checks(ptxt, P):
    """Cross-references inside the protocol resolve: sections, 7(X), lemma items, '12.5.1 step n', '3.1 item n', '12.3 property n', '12.5.3(g)n'."""
    lines = ptxt.split("\n")
    def block(start_pat, end_pat):
        i = next(k for k, l in enumerate(lines) if re.search(start_pat, l))
        j = next(k for k in range(i + 1, len(lines)) if re.search(end_pat, lines[k]))
        return lines[i:j]
    n31 = sum(1 for l in block(r"^\*\*規範性定義", r"^\*\*位元級可重現的範圍") if re.match(r"^\d+\. \*\*", l))
    n123 = sum(1 for l in block(r"^\*\*性質。\*\*", r"^### 12\.4") if re.match(r"^\d+\. ", l))
    ng = sum(1 for l in block(r"^\(g\)〔\*\*限制", r"^\*\*12\.5\.4") if re.match(r"^\d+\. ", l))
    bad = []
    body = "\n".join(l for l in lines if not l.startswith("```"))
    for m in re.finditer(r"§(\d+(?:\.\d+){0,2})", body):
        if m.group(1) not in P["secs"]:
            bad.append("§" + m.group(1))
    for m in re.finditer(r"§?\b7\(([A-Z][\u2032']?)\)", body):
        if norm(m.group(1)) not in P["seven"]:
            bad.append("7(%s)" % m.group(1))
    for m in re.finditer(r"引理\s*(\d)((?:\([a-z]+[\u2032']?\))*)", body):
        n = int(m.group(1))
        if n not in P["lemma_items"]:
            bad.append("引理 %d" % n)
        else:
            for it in ITEM_RE.findall(m.group(2)):
                if norm(it) not in P["lemma_items"][n]:
                    bad.append("引理 %d(%s)" % (n, it))
    for m in re.finditer(r"12\.5\.1 第 (\d) 步", body):
        if not 1 <= int(m.group(1)) <= P["steps_1251"]:
            bad.append(m.group(0))
    for m in re.finditer(r"§3\.1 第 (\d+) 點", body):
        if not 1 <= int(m.group(1)) <= n31:
            bad.append(m.group(0))
    for m in re.finditer(r"12\.3 性質 (\d)", body):
        if not 1 <= int(m.group(1)) <= n123:
            bad.append(m.group(0))
    for m in re.finditer(r"12\.5\.3\(g\)(\d)", body):
        if not 1 <= int(m.group(1)) <= ng:
            bad.append(m.group(0))
    out = [(False, "dangling reference in the protocol: " + b) for b in sorted(set(bad))]
    out.append((not bad, "protocol-internal references resolve (3.1 has %d numbered items, 12.3 has %d properties, 12.5.3(g) has %d points)" % (n31, n123, ng)))
    return out


FORBIDDEN = [(r"建議", "advice wording"), (r"沙盒", "environment-specific wording"), (r"已更正|更正後|已修正|舊式|舊換算|不再引用", "correction history (belongs in CHANGELOG.md)"),
             (r"新增（|新增\)|（V\d|\bV\d\.\d\.\d+", "version history inside the text (only the pointer in the first lines is allowed)"),
             (r"TODO|FIXME|草稿|草案|待補|待定|暫定", "draft marker"), (r"\bAI\b|ChatGPT|Claude|Gemini|語言模型", "reference to an assistant")]


def prose_checks(ptxt):
    """The text contains no draft markers, advice wording, correction history or assistant references (the first 3 lines hold the version pointer)."""
    out = []
    for pat, why in FORBIDDEN:
        for i, ln in enumerate(ptxt.split("\n"), 1):
            if i <= 3:
                continue
            m = re.search(pat, ln)
            if m:
                out.append((False, "prose: line %d %s: '%s' ..." % (i, why, ln[max(0, m.start() - 15):m.end() + 25])))
    out.append((not out, "no draft / advice / correction-history wording in the protocol text"))
    return out


def citation_checks(ptxt):
    """Files and tests that the protocol cites by name must exist.  Returns list of (ok, description)."""
    out = []
    for m in sorted(set(re.findall(r"`((?:legacy|code|results|docs|protocol|metadata)/[\w./-]+\.\w+)`", ptxt))):
        out.append((os.path.exists(os.path.join(ROOT, m)), "cited file exists: " + m))
    ttxt = open(os.path.join(CODE, "tests", "test_snd.py"), encoding="utf8").read()
    tests = set(re.findall(r"^def (test_\w+)", ttxt, re.M))
    for t in sorted(set(re.findall(r"`?(test_\w+\*?)`?", ptxt.replace("test_snd.py", "")))):
        ok = any(x.startswith(t[:-1]) for x in tests) if t.endswith("*") else t in tests
        out.append((ok, "cited test exists: " + t))
    for cpath in sorted(glob.glob(os.path.join(CODE, "crosscheck", "*.c"))):          # C sources carry the same Implements / Verified by header
        head = open(cpath, encoding="utf8").read()[:2000]
        named = re.findall(r"tests/test_snd\.py::(test_\w+)", head)
        out.append((bool(re.search(r"Implements: protocol", head)) and bool(named) and all(t in tests for t in named),
                    "C source %s: header names the protocol sections and existing tests" % os.path.basename(cpath)))
    bare = sorted(set(re.findall(r"(?<![/\w])`(r\d+[a-z]?\.py)`", ptxt)))
    out += [(False, "evidence script cited without path (use legacy/...): " + b) for b in bare]
    return out


def _private_words():
    """Words that must not appear anywhere in the repository (the author's personal identifiers).  They are NOT stored in the repository -- a list of them
    would itself publish them -- but read from a local file, one word per line: the path in the environment variable SND_IDENTITY_WORDS, else ~/.snd_private_words."""
    path = os.environ.get("SND_IDENTITY_WORDS") or os.path.join(os.path.expanduser("~"), ".snd_private_words")
    if not os.path.isfile(path):
        return None
    return [w.strip() for w in open(path, encoding="utf8").read().splitlines() if w.strip()]


def identity_checks(root, words=None):
    """No text file of the repository contains one of the author's personal identifiers (case-insensitive).  `words`: the list to look for; by default the local
    file of _private_words(); without any list the check reports 'SKIPPED' (it is the author's own hygiene check and cannot be reproduced by others).
    Returns [(ok, description)]."""
    words = _private_words() if words is None else words
    if not words:
        return [(True, "SKIPPED: personal-identifier check (no local word list; set SND_IDENTITY_WORDS to a file with one word per line)")]
    pat = re.compile("|".join(re.escape(w) for w in words), re.I)
    hits = []
    for dirpath, dirs, names in os.walk(root):
        dirs[:] = [d for d in dirs if d not in (".git", "__pycache__", ".pytest_cache")]
        for n in names:
            path = os.path.join(dirpath, n)
            if n.endswith((".npz", ".pyc", ".zip")) or os.path.getsize(path) > 5 * 1024 * 1024:
                continue
            if pat.search(open(path, "rb").read().decode("utf8", "ignore")):
                hits.append(os.path.relpath(path, root))
    return [(not hits, "no personal identifier in any repository text file" + ("" if not hits else ": " + ", ".join(sorted(hits))))]


def version_checks(root, v, ntests):
    """Version label and test count agree in every file that states them.  Returns [(ok, description)].
    Implements: protocol 3.1 (the protocol version string is part of the seed material, so the label is part of the definition)."""
    import json
    out = []

    def rd(rel):
        path = os.path.join(root, *rel.split("/"))
        return open(path, encoding="utf8").read() if os.path.exists(path) else ""

    def add(ok, d):
        out.append((bool(ok), d))

    cit = rd("CITATION.cff")
    add(re.search(r'^version:\s*"%s"' % re.escape(v), cit, re.M), "CITATION.cff: version %s" % v)
    add(("Protocol V%s " % v) in cit and ("(V%s)" % v) in cit, "CITATION.cff: title and abstract name V%s" % v)
    zen = []
    for rel in (".zenodo.json", "metadata/zenodo.json"):
        try:
            z = json.loads(rd(rel))
        except ValueError:
            z = {}
        zen.append(rd(rel))
        add(z.get("version") == v and ("V%s" % v) in z.get("title", "") and ("(V%s)" % v) in z.get("description", ""), "%s: version, title and description name %s" % (rel, v))
    add(zen[0] == zen[1] and zen[0] != "", ".zenodo.json and metadata/zenodo.json are identical")
    pdir = os.path.join(root, "protocol")
    cur = sorted(f for f in os.listdir(pdir) if f.endswith(".md")) if os.path.isdir(pdir) else []
    add(cur == ["SND_protocol_V%s.md" % v], "protocol/ holds exactly SND_protocol_V%s.md (found %s)" % (v, cur))
    add(rd("protocol/SND_protocol_V%s.md" % v).startswith("# Project SND \u2014 \u9810\u8a3b\u518a\u898f\u7bc4 V%s" % v), "protocol title names V%s" % v)
    add(not os.path.exists(os.path.join(pdir, "history", "SND_protocol_V%s.md" % v)), "protocol/history/ does not contain the current version")
    m = re.match(r"2\.0\.(\d+)$", v)
    n = int(m.group(1)) if m else None
    ver = rd("VERSIONING.md")
    add(("\u76ee\u524d\u7248\u672c\u70ba **V%s**" % v) in ver and ("`protocol/SND_protocol_V%s.md`" % v) in ver, "VERSIONING.md: current version and file are V%s" % v)
    if n:
        add(("**V2.0.01 \u2013 V%s** | %d |" % (v, n)) in ver and ("\u5408\u8a08 %d \u500b\u6709\u7248\u865f\u7684\u7248\u672c\uff0816 + 8 + %d\uff09" % (24 + n, n)) in ver,
            "VERSIONING.md: period six has %d versions, total %d" % (n, 24 + n))
    rm = rd("README.md")
    add(rm.startswith("# Project SND \u2014 Pre-registration Protocol V%s," % v) and ("SND_protocol_V%s.md" % v) in rm and ("SND_audit_V1.0.08_to_V%s.md" % v) in rm,
        "README.md: title, protocol file and audit file name V%s" % v)
    add(("| %d tests " % ntests) in rm and ("# %d PASS" % ntests) in rm, "README.md: %d tests" % ntests)
    cl, cl2 = rd("CHANGELOG.md"), rd("docs/CHANGELOG.md")
    add(cl == cl2 and cl != "", "CHANGELOG.md and docs/CHANGELOG.md are identical")
    heads = re.findall(r"^### (V[\d.]+)\uff08[^\n]*\u76ee\u524d\u7248\u672c\uff09", cl, re.M)
    add(heads == ["V" + v], "CHANGELOG.md: exactly one entry is marked as the current version, and it is V%s (found %s)" % (v, heads))
    if n:
        add(("V2.0.01 \u2013 V%s\uff08\u5171 %d \u7248\uff09" % (v, n)) in cl, "CHANGELOG.md: period six header has %d versions" % n)
    entry = cl.split("### V%s\uff08" % v)[-1].split("\n## ")[0]
    add(("%d \u9805\u6e2c\u8a66" % ntests) in entry, "CHANGELOG.md entry V%s: %d tests" % (v, ntests))
    aud = rd("docs/audits/SND_audit_V1.0.08_to_V%s.md" % v)
    add(aud.startswith("# SND \u7a3d\u6838 V%s" % v) and ("%d \u9805\u6e2c\u8a66" % ntests) in aud and ("%d PASS" % ntests) in aud, "audit file for V%s exists and states %d tests" % (v, ntests))
    add(("**\u76ee\u524d\u7248\u672c\uff1aV%s\u3002**" % v) in rd("docs/SND_evolution_for_everyone.md"), "docs/SND_evolution_for_everyone.md: current version V%s" % v)
    summ = rd("docs/SND_evolution_summary.md")
    add(("\u76ee\u524d\u7248\u672c V%s" % v) in summ and ("SND_audit_V1.0.08_to_V%s.md" % v) in summ, "docs/SND_evolution_summary.md: current version V%s" % v)
    osf = rd("OSF_REGISTRATION_GUIDE.md")
    add(("SND_release_V%s.zip" % v) in osf and ("version `%s`" % v) in osf and ("(protocol V%s)" % v) in osf, "OSF_REGISTRATION_GUIDE.md names V%s" % v)
    return out


def table(rows, P):
    by = {}
    for r in rows:
        if not r["impl"]:
            continue
        for ref in parse_refs(r["impl"]):
            label = ("%s%s" % (ref[1], "(%s)" % ref[2] if ref[2] else "")) if ref[0] == "sec" else "Lemma %d%s" % (ref[1], "".join("(%s)" % i for i in ref[2]))
            by.setdefault(label, []).append(r)
    def key(lbl):
        m = re.match(r"(Lemma )?(\d+(?:\.\d+)*)(.*)", lbl)
        return (1 if m.group(1) else 0, [int(x) for x in m.group(2).split(".")], m.group(3))
    lines = ["# 規範與程式碼對照表（SND_protocol 與 code/）", "",
             "本表由 `code/scripts/check_trace.py --write` 自 docstring 產生，請勿手改；`--check` 會比對本檔。程式碼的註解以 `Implements:` 標出所實作的規範章節（§ 為規範標題編號；`7(D)` 為 §7 的 (D) 項；`Lemma 5(c)` 為引理 5 的 (c) 項），`Verified by:` 標出測試。規範本身沒有公式編號。",
             "", "## 規範項目 → 程式碼", "", "| 規範項目 | 程式碼 | 測試 |", "|---|---|---|"]
    for lbl in sorted(by, key=key):
        for r in by[lbl]:
            sym = "%s" % r["file"] if r["kind"] == "module" else "%s::%s" % (r["file"], r["name"])
            lines.append("| %s | `%s` | %s |" % (lbl, sym, (r["tests"] or "").replace("tests/test_snd.py::", "").replace("|", "/") or "—"))
    cited = {l for l in by}
    uncited = sorted((s for s in P["secs"] if not any(c == s or c.startswith(s + "(") for c in cited)), key=lambda s: [int(x) for x in s.split(".")])
    lines += ["", "## 沒有程式碼對應的規範章節", "",
              "下列章節沒有被任何 `Implements:` 引用（敘述性、程序性或硬體章節，或其內容由上列相鄰項目涵蓋）：" + "、".join("§" + s for s in uncited) + "。", ""]
    return "\n".join(lines)


def main(argv):
    import snd
    ppath = None
    if "--protocol" in argv:
        ppath = argv[argv.index("--protocol") + 1]
    else:
        ppath = os.path.join(ROOT, "protocol", "SND_protocol_V%s.md" % snd.__version__)
    ptxt = open(ppath, encoding="utf8").read()
    P = read_protocol(ppath)
    rows = collect()
    fails = []
    # 1. every reference exists
    nref = 0
    for r in rows:
        if r["impl"]:
            refs = parse_refs(r["impl"])
            nref += len(refs)
            if not refs:
                fails.append("%s::%s: 'Implements:' has no protocol reference" % (r["file"], r["name"]))
            for b in validate(refs, P):
                fails.append("%s::%s: %s" % (r["file"], r["name"], b))
    # 2. every public function / class in snd/ cites the protocol
    for r in rows:
        if r["file"].startswith("snd/") and r["kind"] != "module" and not r["name"].startswith("_") and not r["impl"]:
            fails.append("%s::%s has no 'Implements:' line" % (r["file"], r["name"]))
        if r["file"].startswith("snd/") and r["kind"] == "module" and not r["impl"] and not r["file"].endswith("__init__.py"):
            fails.append("%s has no 'Implements:' line" % r["file"])
    # 3. named tests exist
    ttxt = open(os.path.join(CODE, "tests", "test_snd.py"), encoding="utf8").read()
    tests = set(re.findall(r"^def (test_\w+)", ttxt, re.M))
    for r in rows:
        for t in re.findall(r"(?:::|,)\s*(test_\w+\*?)", r["tests"] or ""):
            if t.endswith("*"):
                if not any(x.startswith(t[:-1]) for x in tests):
                    fails.append("%s::%s: no test matches %s" % (r["file"], r["name"], t))
            elif t not in tests:
                fails.append("%s::%s: test %s does not exist" % (r["file"], r["name"], t))
    # 4. constants
    cons = constants_checks(ptxt)
    for ok, d in cons:
        print("%s  %s" % ("OK  " if ok else "FAIL", d))
        if not ok:
            fails.append("constant mismatch: " + d)
    for ok, d in internal_ref_checks(ptxt, P):
        print("%s  %s" % ("OK  " if ok else "FAIL", d))
        if not ok:
            fails.append(d)
    for ok, d in prose_checks(ptxt):
        print("%s  %s" % ("OK  " if ok else "FAIL", d))
        if not ok:
            fails.append(d)
    for ok, d in citation_checks(ptxt):
        if not ok:
            print("FAIL  " + d)
            fails.append(d)
    print("citations checked in the protocol: %d" % len(citation_checks(ptxt)))
    ntests = len(tests)
    for ok, d in identity_checks(ROOT):
        print("%s  %s" % (("SKIP" if d.startswith("SKIPPED") else "OK  ") if ok else "FAIL", d))
        if not ok:
            fails.append(d)
    for ok, d in version_checks(ROOT, snd.__version__, ntests):
        print("%s  %s" % ("OK  " if ok else "FAIL", d))
        if not ok:
            fails.append("version label: " + d)
    # 5. generated table
    gen = table(rows, P) + "\n"
    tpath = os.path.join(ROOT, "docs", "TRACEABILITY.md")
    if "--write" in argv:
        open(tpath, "w", encoding="utf8", newline="\n").write(gen)
        print("wrote", os.path.relpath(tpath, ROOT))
    elif "--check" in argv:
        if not os.path.exists(tpath) or open(tpath, encoding="utf8").read() != gen:
            fails.append("docs/TRACEABILITY.md is not the generated table (run: python scripts/check_trace.py --write)")
    print("protocol: %s; sections %d, lemmas %s, 7-items %s, 12.5.1 steps %d; code references checked: %d in %d docstrings"
          % (os.path.basename(ppath), len(P["secs"]), sorted(P["lemma_items"]), sorted(P["seven"]), P["steps_1251"], nref, sum(1 for r in rows if r["impl"])))
    for f in fails:
        print("FAIL", f)
    print("trace check:", "OK" if not fails else "%d problem(s)" % len(fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
