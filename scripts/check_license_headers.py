#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check that Maxent source files carry the ALPS project and SPDX notices.

Maxent is MIT-licensed (see LICENSE.txt) and follows the ALPS header
convention: every C/C++ source keeps its copyright attribution line and ends
its header with

    ALPS Project: https://alps.comp-phys.org/
    SPDX-License-Identifier: MIT

It also reports leftover phrases from the superseded GPL and "All rights
reserved" notices.  Third-party files keep their own notices and are skipped.

    python3 scripts/check_license_headers.py      # exit 1 if a file fails
"""

import re
import subprocess
import sys
from pathlib import Path

SOURCE_PATTERNS = ["*.cpp", "*.hpp", "*.h", "*.cc", "*.in"]
THIRD_PARTY = re.compile(r"^(test/gtest[^/]*|cmake/Find[^/]*\.cmake)$")
HEADER_LINES = 15
REQUIRED = [
    "ALPS Project: https://alps.comp-phys.org/",
    "SPDX-License-Identifier: MIT",
]
LEGACY = re.compile(
    r"GNU General Public|All rights reserved\. Use is subject|ACKNOWLEDGE\.TXT|LICENSE\.TXT"
)


def main():
    root = Path(__file__).resolve().parent.parent
    files = subprocess.run(
        ["git", "ls-files", *SOURCE_PATTERNS],
        cwd=root, capture_output=True, text=True, check=True,
    ).stdout.split()
    failures = []
    for rel in files:
        if THIRD_PARTY.match(rel):
            continue
        text = (root / rel).read_text(errors="replace")
        head = "\n".join(text.splitlines()[:HEADER_LINES])
        missing = [r for r in REQUIRED if r not in head]
        if missing:
            failures.append(f"{rel}: missing {', '.join(repr(m) for m in missing)}")
        if LEGACY.search(text):
            failures.append(f"{rel}: contains a superseded license notice")
    for f in failures:
        print(f)
    print(f"checked {len(files)} files, {len(failures)} problem(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
