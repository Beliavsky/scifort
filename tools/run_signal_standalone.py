#!/usr/bin/env python3
"""Compile and run signal tests without fpm, using gfortran.

Usage: python tools/run_signal_standalone.py checked|optimized BUILD_DIR [--match SUBSTRING]
"""
from pathlib import Path
import re
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
mode = sys.argv[1] if len(sys.argv) > 1 else "checked"
if mode not in ("checked", "optimized"):
    raise SystemExit("expected checked or optimized")
build = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else root.parent / f"build_signal_{mode}"
build.mkdir(parents=True, exist_ok=True)
match = None
if len(sys.argv) > 3:
    if len(sys.argv) != 5 or sys.argv[3] != "--match" or not sys.argv[4]:
        raise SystemExit("optional arguments: --match TEST_NAME_SUBSTRING")
    match = sys.argv[4].lower()
flags = ["-std=f2018", "-Wimplicit-interface", "-Werror=implicit-interface"]
if mode == "checked":
    flags += ["-O0", "-g", "-fcheck=all", "-fbacktrace"]
else:
    flags += ["-O3"]
modules = {}
for source in (root / "src").rglob("*.f90"):
    for name in re.findall(r"^\s*module\s+(?!procedure|function|subroutine)([a-z]\w*)",
                           source.read_text(), re.I | re.M):
        modules[name.lower()] = source
seen, current, ordered = set(), set(), []

def visit(module):
    if module in seen:
        return
    if module in current:
        raise RuntimeError(f"circular module dependency: {module}")
    if module not in modules:
        raise RuntimeError(f"unknown module: {module}")
    current.add(module)
    source = modules[module]
    for dep in re.findall(r"^\s*use\s+(scifort_\w+)",source.read_text(),re.I|re.M):
        visit(dep.lower())
    current.remove(module)
    seen.add(module)
    if source not in ordered:
        ordered.append(source)

visit("scifort_signal")
tests = sorted((root / "test").glob("test_signal_*.f90"))
if match is not None:
    tests = [source for source in tests if match in source.stem.lower()]
    if not tests:
        raise SystemExit(f"no signal tests matched {match!r}")
for source in tests:
    for dep in re.findall(r"^\s*use\s+(scifort_\w+)",source.read_text(),re.I|re.M):
        visit(dep.lower())
print(f"{mode}: compiling {len(ordered)} sources and {len(tests)} signal tests",flush=True)
objects = []
for i, source in enumerate(ordered):
    obj = build / (source.stem + ".o")
    command = ["gfortran",*flags,"-I",str(build),"-J",str(build),"-c",str(source),"-o",str(obj)]
    result = subprocess.run(command,capture_output=True,text=True,cwd=build)
    if result.returncode:
        print(result.stderr,file=sys.stderr)
        raise SystemExit(f"compilation failed: {source.name}")
    objects.append(str(obj))
    if (i+1)%15 == 0:
        print(f"compiled {i+1}/{len(ordered)}",flush=True)
print("all modules compiled",flush=True)
passed = 0
for source in tests:
    exe = build / source.stem
    cmd = ["gfortran",*flags,"-I",str(build),str(source),*objects,"-o",str(exe)]
    link = subprocess.run(cmd,capture_output=True,text=True,cwd=build)
    if link.returncode:
        print(link.stderr,file=sys.stderr)
        raise SystemExit(f"test link failed: {source.name}")
    run = subprocess.run([str(exe)],capture_output=True,text=True,cwd=build,timeout=90)
    if run.returncode:
        print(run.stdout,run.stderr,file=sys.stderr)
        raise SystemExit(f"test failed: {source.name}")
    passed += 1
    if passed%20 == 0:
        print(f"passed {passed}/{len(tests)}",flush=True)
print(f"PASS: {passed}/{len(tests)} signal tests ({mode})",flush=True)
