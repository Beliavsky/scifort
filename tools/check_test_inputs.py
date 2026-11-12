#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Keep explicit console reads out of automatically discovered FPM tests."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
CONSOLE_READ = re.compile(r"\bread\s*\(\s*(?:\*|5\b)", re.IGNORECASE)


def main() -> int:
    failures = []
    for path in sorted((ROOT / "test").rglob("*.f90")):
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if CONSOLE_READ.search(line.split("!", 1)[0]):
                failures.append(f"{path.relative_to(ROOT)}:{number}: console input in FPM test")
    for failure in failures:
        print(failure)
    if failures:
        print("Move input-driven comparison programs to validation/reference_drivers/.")
    else:
        print("FPM test console-input check passed")
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
