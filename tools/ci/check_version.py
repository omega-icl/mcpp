#!/usr/bin/env python3
"""check_version.py -- the version is written in three places; they must agree.

    pyproject.toml            version = "X.Y.Z"        what pip and the wheels report
    src/pymcpp/pymcpp.cpp     m.attr("__version__")    what `pymcpp.__version__` reports
    doc/MC++.dox              PROJECT_NUMBER           what the Doxygen pages are titled with

They drifted apart twice in 5.0.x -- MC++.dox stayed at 5.0.3 through 5.0.4, and the tag was the only
place that knew (2026-10-26).  Run from the repository root, or with the root as the first argument;
exits non-zero on a mismatch, naming each file's value.  With --expect X.Y.Z, also requires that value
(the release workflow passes the tag's).
"""
import os
import re
import sys


def read(path, pattern):
    with open(path, encoding="utf-8") as f:
        m = re.search(pattern, f.read(), re.M)
    return m.group(1) if m else None


def main(argv):
    root = "."
    expect = None
    args = list(argv[1:])
    while args:
        a = args.pop(0)
        if a == "--expect":
            expect = args.pop(0)
        else:
            root = a
    sources = [
        ("pyproject.toml", os.path.join(root, "pyproject.toml"), r'^version\s*=\s*"([^"]+)"'),
        ("pymcpp.__version__", os.path.join(root, "src", "pymcpp", "pymcpp.cpp"), r'm\.attr\(\s*"__version__"\s*\)\s*=\s*"([^"]+)"'),
        ("doc/MC++.dox PROJECT_NUMBER", os.path.join(root, "doc", "MC++.dox"), r'^PROJECT_NUMBER\s*=\s*(\S+)'),
    ]
    found = []
    for name, path, pat in sources:
        v = read(path, pat) if os.path.exists(path) else None
        found.append((name, v))
        print("  %-30s %s" % (name, v if v else "(not found)"))
    values = {v for _, v in found}
    ok = len(values) == 1 and None not in values
    if ok and expect is not None and expect not in values:
        print("check_version.py: FAIL -- the files say %s, the tag says %s" % (values.pop(), expect))
        return 1
    if not ok:
        print("check_version.py: FAIL -- the version is not the same everywhere")
        return 1
    print("check_version.py: the version is %s everywhere" % values.pop())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
