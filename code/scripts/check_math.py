"""Formula hygiene check for every Markdown file in the repository (protocol V2.0.11 and later).

Formulas are written in GitHub's math syntax: inline  $`...`$  and fenced display blocks  ```math ... ```.
The older  \\( ... \\)  /  \\[ ... \\]  delimiters do NOT render on GitHub (Markdown treats  \\(  as an escaped parenthesis).

Checks (exit 0 = every check that ran passed):
  1. no raw  \\(  \\)  \\[  \\]  delimiters outside fenced code and inline code;
  2. every math span is non-empty, has balanced { } (escaped \\{ \\} ignored), matching \\left/\\right and \\begin/\\end;
  3. no stray '$' outside math spans and code;
  4. every formula parses in KaTeX -- only if node and the 'katex' package are available
     (set KATEX_DIR to the package directory, or make 'katex' resolvable from node); otherwise reported as SKIPPED.

Usage:  python scripts/check_math.py [ROOT]      (default ROOT: repository root)
Note: KaTeX is a close but not identical parser to the MathJax that github.com uses; a KaTeX pass is strong evidence, not proof of GitHub rendering.
"""
import json, os, re, subprocess, sys, tempfile

SKIP_DIRS = {".git", "node_modules", "__pycache__"}
INLINE = re.compile(r"\$`(.+?)`\$")
DBL_CODE = re.compile(r"``(.+?)``")
SGL_CODE = re.compile(r"`[^`]*`")


def md_files(root):
    for d, ds, fs in os.walk(root):
        ds[:] = [x for x in ds if x not in SKIP_DIRS]
        for f in sorted(fs):
            if f.endswith(".md"):
                yield os.path.join(d, f)


def extract(path):
    """-> (formulas [(line, tex, display)], problems [(line, msg)])"""
    forms, probs = [], []
    fence = None            # None | 'math' | 'code'
    buf, start = [], 0
    for n, line in enumerate(open(path, encoding="utf8").read().split("\n"), 1):
        s = line.strip()
        if fence is None and s.startswith("```"):
            fence = "math" if s[3:].strip() == "math" else "code"
            buf, start = [], n
            continue
        if fence is not None:
            if s.startswith("```") and s.strip("`") == "":
                if fence == "math":
                    forms.append((start, "\n".join(buf).strip(), True))
                fence = None
            else:
                buf.append(line)
            continue
        t = DBL_CODE.sub("", line)                 # literal examples written as ``...``
        for m in INLINE.finditer(t):
            forms.append((n, m.group(1), False))
        t = INLINE.sub("", t)
        t = SGL_CODE.sub("", t)
        for d in ("\\(", "\\)", "\\[", "\\]"):
            if d in t:
                probs.append((n, "raw delimiter %s outside code (does not render on GitHub)" % d))
        if "$" in t:
            probs.append((n, "stray '$' outside math spans and code"))
    if fence is not None:
        probs.append((start, "unclosed fenced block"))
    return forms, probs


def structure_problems(tex):
    out, depth, i = [], 0, 0
    while i < len(tex):
        c = tex[i]
        if c == "\\":
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth < 0:
                out.append("unbalanced braces (extra '}')")
                depth = 0
        i += 1
    if depth:
        out.append("unbalanced braces (%d unclosed '{')" % depth)
    if len(re.findall(r"\\left(?![a-zA-Z])", tex)) != len(re.findall(r"\\right(?![a-zA-Z])", tex)):
        out.append("\\left/\\right count mismatch")
    stack = []
    for m in re.finditer(r"\\(begin|end)\{([^}]*)\}", tex):
        if m.group(1) == "begin":
            stack.append(m.group(2))
        elif not stack or stack.pop() != m.group(2):
            out.append("\\begin/\\end mismatch")
            break
    if stack:
        out.append("\\begin without \\end")
    if not tex.strip():
        out.append("empty formula")
    return out


KATEX_JS = r"""
const path=process.env.KATEX_DIR||'katex';let katex;try{katex=require(path);}catch(e){console.log('NOKATEX');process.exit(3);}
const items=JSON.parse(require('fs').readFileSync(process.argv[2],'utf8'));const seen=new Map();const bad=[];
for(const [file,line,tex,disp] of items){const k=(disp?'D:':'I:')+tex;if(!seen.has(k)){let err=null;
 try{katex.renderToString(tex,{displayMode:disp,throwOnError:true,strict:'ignore',trust:false});}catch(e){err=String(e.message).split('\n')[0];}seen.set(k,err);}
 if(seen.get(k))bad.push([file,line,tex.slice(0,80),seen.get(k)]);}
console.log(JSON.stringify({total:items.length,unique:seen.size,bad}));process.exit(bad.length?1:0);
"""


def katex_check(items):
    with tempfile.TemporaryDirectory() as td:
        jf, sf = os.path.join(td, "f.json"), os.path.join(td, "k.js")
        json.dump(items, open(jf, "w", encoding="utf8"), ensure_ascii=False)
        open(sf, "w").write(KATEX_JS)
        try:
            r = subprocess.run(["node", sf, jf], capture_output=True, text=True, timeout=600, env=dict(os.environ))
        except FileNotFoundError:
            return None, "node not found"
        if r.returncode == 3 or "NOKATEX" in r.stdout:
            return None, "katex package not found (set KATEX_DIR)"
        try:
            return json.loads(r.stdout.strip().splitlines()[-1]), None
        except Exception:
            # katex and node ARE present but the runner broke: that is a failure, not a skip
            return {"total": 0, "unique": 0, "bad": [["<runner>", 0, "", "KaTeX runner failed: " + (r.stdout + r.stderr)[:200]]]}, None


def main(root):
    root = os.path.abspath(root)
    bad = 0
    allf, nfiles, nform = [], 0, 0
    for p in md_files(root):
        rel = os.path.relpath(p, root)
        forms, probs = extract(p)
        for (ln, tex, disp) in forms:
            for msg in structure_problems(tex):
                probs.append((ln, msg + ": " + tex[:70]))
        if forms:
            nfiles += 1
            nform += len(forms)
        for ln, msg in probs:
            bad += 1
            print("FAIL %s:%d %s" % (rel, ln, msg))
        allf += [[rel, ln, tex, disp] for ln, tex, disp in forms]
    print("%-5s hygiene: %d files with formulas, %d formulas, %d problems" % ("OK" if not bad else "FAIL", nfiles, nform, bad))
    res, why = katex_check(allf)
    if res is None:
        print("SKIPPED KaTeX parse: %s" % why)
    else:
        for f, ln, tex, err in res["bad"]:
            print("FAIL %s:%d KaTeX: %s | %s" % (f, ln, err, tex))
        bad += len(res["bad"])
        print("%-5s KaTeX parse: %d formulas (%d distinct), %d failures" % ("OK" if not res["bad"] else "FAIL", res["total"], res["unique"], len(res["bad"])))
    return 1 if bad else 0


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "..")))
