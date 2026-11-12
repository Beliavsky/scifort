#!/usr/bin/env python3
"""Compile/test arbitrary-center exponential windows against SciPy."""
from pathlib import Path
import subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parents[1]
MODS=('kinds','signal_windows_extra','signal_dpss','signal_spectral')
def run(build,optimized):
    build.mkdir()
    flags=['-std=f2018','-fimplicit-none','-ffree-line-length-none','-I.','-J.']
    flags+=['-O3'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    objs=[]
    for mod in MODS:
        subprocess.run(['gfortran',*flags,'-c',str(ROOT/'src'/(mod+'.f90')),
                        '-o',mod+'.o'],cwd=build,check=True)
        objs.append(mod+'.o')
    for name, folder in [('test_signal_exponential_center','test'),
                         ('exponential_center_driver','tools/signal_reference')]:
        exe=build/name
        subprocess.run(['gfortran',*flags,str(ROOT/folder/(name+'.f90')),*objs,
                        '-o',str(exe)],cwd=build,check=True)
        if folder=='test': subprocess.run([str(exe)],cwd=build,check=True)
        else: subprocess.run([sys.executable,str(ROOT/'tools/signal_reference/verify_exponential_center.py'),
                             str(exe)],cwd=build,check=True)
    print('PASS: arbitrary-center exponential window', 'optimized' if optimized else 'checked')
if __name__=='__main__':
    with tempfile.TemporaryDirectory(prefix='scifort_exponential_center_') as t:
        run(Path(t)/'checked',False)
        if '--checked-only' not in sys.argv: run(Path(t)/'optimized',True)
