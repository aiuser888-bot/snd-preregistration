"""Write or verify SHA256SUMS.txt for the whole release tree.

    python scripts/make_checksums.py --write        # rewrite SHA256SUMS.txt from the files now on disk
    python scripts/make_checksums.py --check        # exit 0 only if the list is exactly the files on disk and every hash matches
    python scripts/make_checksums.py --check --root DIR   # DIR: a checked-out or exported tree (for example `git archive HEAD | tar -x -C DIR`)

The list covers every file below the root except: SHA256SUMS.txt itself, `.git/`, `.github/`, `.gitignore`, and generated caches
(`__pycache__/`, `*.pyc`, `.pytest_cache/`).  Format: "<sha256><two spaces>./<path>", sorted by path (byte order), as `sha256sum -c` reads it.

Order of work at a release: finish every edit, run --write, then export the tree that will be tagged and run --check on the export.
Writing the list before the last edit leaves a stale hash (this happened once: see CHANGELOG, V2.0.13, erratum for v2.0.12).
"""
import argparse
import hashlib
import os
import sys

NAME = "SHA256SUMS.txt"
SKIP_DIRS = {".git", ".github", "__pycache__", ".pytest_cache"}
SKIP_FILES = {NAME, ".gitignore"}


def tree_files(root):
    """Relative POSIX paths of all files that the list must cover, sorted in byte order."""
    out = []
    for d, dirs, files in os.walk(root):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
        for f in files:
            rel = os.path.relpath(os.path.join(d, f), root).replace(os.sep, "/")
            if rel in SKIP_FILES or f.endswith(".pyc"):
                continue
            out.append(rel)
    return sorted(out, key=lambda s: s.encode("utf-8"))


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def render(root):
    return "".join("%s  ./%s\n" % (sha256_of(os.path.join(root, rel)), rel) for rel in tree_files(root))


def check(root):
    """Return a list of problem strings (empty = OK)."""
    path = os.path.join(root, NAME)
    if not os.path.exists(path):
        return ["%s is missing" % NAME]
    listed, order = {}, []
    problems = []
    for n, line in enumerate(open(path, encoding="utf-8").read().splitlines(), 1):
        parts = line.split("  ", 1)
        if len(parts) != 2 or len(parts[0]) != 64 or not parts[1].startswith("./"):
            problems.append("line %d is not '<sha256>  ./<path>'" % n)
            continue
        if parts[1][2:] in listed:
            problems.append("%s listed twice" % parts[1])
        listed[parts[1][2:]] = parts[0]
        order.append(parts[1][2:])
    files = tree_files(root)
    problems += ["not listed: %s" % f for f in files if f not in listed]
    problems += ["listed but missing: %s" % f for f in listed if f not in set(files)]
    problems += ["hash differs: %s" % f for f in files if f in listed and sha256_of(os.path.join(root, f)) != listed[f]]
    if not problems and order != files:
        problems.append("list is not sorted in byte order")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    a = ap.parse_args(argv)
    root = os.path.abspath(a.root)
    if a.write:
        text = render(root)
        with open(os.path.join(root, NAME), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        print("wrote %s: %d files" % (NAME, text.count("\n")))
        return 0
    problems = check(root)
    for p in problems:
        print("PROBLEM", p)
    n = len(tree_files(root))
    print("checksums: %s (%d files)" % ("OK" if not problems else "%d problem(s)" % len(problems), n))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
