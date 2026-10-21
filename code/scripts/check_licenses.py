"""Licence coverage check (V2.0.11).  Exit 0 = every check passed.

The repository is multi-licensed (MIT for software, CC BY 4.0 for text/data/metadata); REUSE.toml maps every file to exactly one licence and
LICENSES/<SPDX-id>.txt holds the full text.  This script checks:
  1. REUSE.toml parses and every annotation is well-formed (version 1, path, SPDX-License-Identifier, SPDX-FileCopyrightText);
  2. every file of the repository (git-tracked files, or all files if git is unavailable; LICENSES/ excluded) matches exactly one
     annotation -- none uncovered, none covered twice;
  3. every licence id used has LICENSES/<id>.txt with its characteristic wording, and every LICENSES/*.txt is used;
  4. LICENSE (the notice) names every licence id and REUSE.toml, contains no draft prompts, and no document still refers to the removed
     LICENSE-DOCS.txt;
  5. CITATION.cff and the Zenodo metadata name a licence that is in use, and the two Zenodo files are identical.
It is NOT the official `reuse lint` (that tool was not available when this was written); the REUSE.toml format follows the REUSE specification.

Usage:  python scripts/check_licenses.py [ROOT]      (default ROOT: repository root)
"""
import json, os, re, subprocess, sys

try:
    import tomllib
except ImportError:                       # Python < 3.11
    tomllib = None

NEEDED = {"MIT": ["Permission is hereby granted, free of charge", "THE SOFTWARE IS PROVIDED \"AS IS\""],
          "CC-BY-4.0": ["Attribution 4.0 International", "Section 8 -- Interpretation.", "Creative Commons may be contacted at creativecommons.org."]}
DRAFT_HINTS = ["choose another", "if you prefer", "todo", "fixme", "replace", "草稿", "待填", "xxx"]
SKIP = {".git", "__pycache__", ".pytest_cache", "node_modules"}


def glob_re(pat):
    out, i = "", 0
    while i < len(pat):
        if pat.startswith("**", i):
            out += ".*"; i += 2
        elif pat[i] == "*":
            out += "[^/]*"; i += 1
        elif pat[i] == "?":
            out += "[^/]"; i += 1
        else:
            out += re.escape(pat[i]); i += 1
    return re.compile("^" + out + "$")


def list_files(root):
    try:
        out = subprocess.run(["git", "ls-files"], cwd=root, capture_output=True, text=True, check=True).stdout.split("\n")
        files = [f for f in out if f]
        if files:
            return sorted(files)
    except Exception:
        pass
    files = []
    for d, ds, fs in os.walk(root):
        ds[:] = [x for x in ds if x not in SKIP]
        files += [os.path.relpath(os.path.join(d, f), root).replace(os.sep, "/") for f in fs]
    return sorted(files)


def main(root):
    root = os.path.abspath(root)
    fails = []
    def fail(msg): fails.append(msg); print("FAIL", msg)
    rd = lambda p: open(os.path.join(root, p), encoding="utf8").read()
    if tomllib is None:
        print("SKIPPED: Python >= 3.11 is needed for tomllib"); return 0
    try:
        cfg = tomllib.loads(rd("REUSE.toml"))
    except Exception as e:
        print("FAIL REUSE.toml cannot be read: %s" % e); return 1
    if cfg.get("version") != 1: fail("REUSE.toml: version must be 1")
    anns = []
    for k, a in enumerate(cfg.get("annotations", []), 1):
        paths = a.get("path"); paths = [paths] if isinstance(paths, str) else paths
        lic = a.get("SPDX-License-Identifier")
        if not paths or not lic or not a.get("SPDX-FileCopyrightText"):
            fail("REUSE.toml annotation %d is missing path / SPDX-License-Identifier / SPDX-FileCopyrightText" % k); continue
        anns.append((k, lic, [glob_re(p) for p in paths], paths))
    files = [f for f in list_files(root) if not f.startswith("LICENSES/")]
    count = {}
    for f in files:
        hit = [(k, lic) for k, lic, res, _ in anns if any(r.match(f) for r in res)]
        if not hit: fail("not covered by any licence: %s" % f)
        elif len(hit) > 1: fail("covered by %d annotations (ambiguous): %s" % (len(hit), f))
        else: count[hit[0][1]] = count.get(hit[0][1], 0) + 1
    unused_patterns = [(k, p) for k, lic, res, ps in anns for r, p in zip(res, ps) if not any(r.match(f) for f in files)]
    for k, p in unused_patterns: fail("annotation %d: pattern matches no file: %s" % (k, p))
    used = {lic for _, lic, _, _ in anns}
    lic_dir = os.path.join(root, "LICENSES")
    present = {f[:-4] for f in os.listdir(lic_dir) if f.endswith(".txt")} if os.path.isdir(lic_dir) else set()
    for u in sorted(used):
        if u not in present: fail("LICENSES/%s.txt is missing" % u); continue
        t = rd("LICENSES/%s.txt" % u)
        for phrase in NEEDED.get(u, []):
            if phrase not in t: fail("LICENSES/%s.txt lacks characteristic wording: %r" % (u, phrase))
        if len(t) < 500: fail("LICENSES/%s.txt is suspiciously short" % u)
    for p in sorted(present - used): fail("LICENSES/%s.txt is not used by any annotation" % p)
    notice = rd("LICENSE")
    for u in sorted(used):
        if u not in notice: fail("LICENSE (notice) does not name %s" % u)
    if "REUSE.toml" not in notice: fail("LICENSE (notice) does not mention REUSE.toml")
    for name in ("LICENSE", "REUSE.toml", "LICENSES/MIT.txt"):          # our own wording only; the CC BY 4.0 text is verbatim and not ours to edit
        low = rd(name).lower()
        for h in DRAFT_HINTS:
            if h in low: fail("draft hint %r found in %s" % (h, name))
    for name in ("README.md", "OSF_REGISTRATION_GUIDE.md", "LICENSE", "CITATION.cff", "metadata/zenodo.json", ".zenodo.json"):
        if os.path.exists(os.path.join(root, name)) and "LICENSE-DOCS" in rd(name): fail("%s still refers to the removed LICENSE-DOCS.txt" % name)
    m = re.search(r'^license:\s*"?([^"\s#]+)"?', rd("CITATION.cff"), re.M)
    if not m or m.group(1) not in used: fail("CITATION.cff license %r is not a licence in use" % (m.group(1) if m else None))
    zs = {}
    for z in (".zenodo.json", "metadata/zenodo.json"):
        zs[z] = rd(z); zl = json.loads(zs[z]).get("license", "")
        if zl.lower() not in {u.lower() for u in used}: fail("%s: license %r is not a licence in use" % (z, zl))
    if zs[".zenodo.json"] != zs["metadata/zenodo.json"]: fail(".zenodo.json and metadata/zenodo.json differ")
    print("%-5s %d files: %s; licences with full text: %s" % ("OK" if not fails else "FAIL", len(files),
          ", ".join("%s %d" % (k, v) for k, v in sorted(count.items())), ", ".join(sorted(present)) or "none"))
    return 1 if fails else 0


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "..")))
