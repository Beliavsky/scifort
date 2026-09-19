#!/usr/bin/env python3
"""Fail when project text files contain non-ASCII characters."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SUFFIXES = {".f90", ".h", ".md", ".toml", ".yml", ".yaml", ".txt"}
SKIP_PARTS = {"build", ".git"}


def main() -> int:
    failures: list[tuple[Path, int, str]] = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in SUFFIXES:
            continue
        if any(part in SKIP_PARTS for part in path.parts):
            continue
        text = path.read_text(encoding="utf-8")
        for line_number, line in enumerate(text.splitlines(), start=1):
            try:
                line.encode("ascii")
            except UnicodeEncodeError:
                failures.append((path.relative_to(ROOT), line_number, line))

    for path, line_number, line in failures:
        print(f"{path}:{line_number}: non-ASCII text: {line!r}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
