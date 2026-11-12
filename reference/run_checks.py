"""Compile new signal references with gfortran and compare with installed SciPy.

From the SciFort root:
    python reference/run_checks.py
Requires NumPy and SciPy for the reference comparisons, not for SciFort.
"""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
BUILD = ROOT / 'build' / 'signal_reference'
BUILD.mkdir(parents=True, exist_ok=True)
flags = ['-std=f2018', '-ffree-line-length-none', '-O3', '-fcheck=bounds',
         '-J' + str(BUILD), '-I' + str(BUILD)]
objects = {}

for name in ('kinds', 'signal_spectral', 'signal_envelope', 'signal_deconvolve'):
    obj = BUILD / (name + '.o')
    subprocess.run(['gfortran', *flags, '-c', str(ROOT / 'src' / (name + '.f90')),
                    '-o', str(obj)], check=True)
    objects[name] = obj

programs = {
    'ref_envelope': ('signal_envelope_reference_driver', ('kinds', 'signal_spectral', 'signal_envelope')),
    'ref_envelope_2d': ('signal_envelope_2d_driver', ('kinds', 'signal_spectral', 'signal_envelope')),
    'ref_envelope_3d': ('signal_envelope_3d_driver', ('kinds', 'signal_spectral', 'signal_envelope')),
    'ref_deconvolve': ('signal_deconvolve_reference_driver', ('kinds', 'signal_deconvolve')),
}
for exe, (source, deps) in programs.items():
    output = BUILD / (exe + ('.exe' if os.name == 'nt' else ''))
    subprocess.run(['gfortran', *flags, str(ROOT / 'reference' / (source + '.f90')),
                    *(str(objects[dep]) for dep in deps), '-o', str(output)], check=True)

for name in ('check_deconvolve', 'check_envelope', 'check_envelope_2d', 'check_envelope_3d'):
    subprocess.run([sys.executable, str(ROOT / 'reference' / (name + '.py'))], check=True)
print('PASS: 1,900 SciPy comparisons')
