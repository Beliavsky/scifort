#!/usr/bin/env python3
"""Build multirate-only module dependencies for rapid regression testing."""
from pathlib import Path
import sys,re,subprocess
root=Path(__file__).resolve().parents[1]
build=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else root.parent/'multirate_checked'
mode=sys.argv[2] if len(sys.argv)>2 else 'checked'
build.mkdir(parents=True,exist_ok=True)
flags=['-std=f2018','-Wimplicit-interface','-Werror=implicit-interface','-ffree-line-length-132']
if mode=='checked': flags+=['-O0','-g','-fcheck=all','-fbacktrace']
else: flags+=['-O3']
modules={}
for f in (root/'src').rglob('*.f90'):
    for module in re.findall(r'^\s*module\s+(?!procedure|function|subroutine)([a-z]\w*)', f.read_text(), re.I|re.M):
        modules[module.lower()]=f
ordered=[]; seen=set()
def visit(m):
    if m in seen: return
    seen.add(m)
    f=modules[m]
    for dep in re.findall(r'^\s*use\s+(scifort_\w+)', f.read_text(), re.I|re.M):
        visit(dep.lower())
    if f not in ordered: ordered.append(f)
for m in ('scifort_signal_resample','scifort_signal_fourier_resample'): visit(m)
objs=[]
for i,f in enumerate(ordered):
    out=build/(f.stem+'.o')
    p=subprocess.run(['gfortran',*flags,'-I',str(build),'-J',str(build),'-c',str(f),'-o',str(out)],capture_output=True,text=True,cwd=build)
    if p.returncode:
        print(p.stdout+p.stderr,file=sys.stderr)
        raise SystemExit('ERROR '+str(f))
    objs.append(str(out))
print('Compiled',len(objs),'modules,',mode,flush=True)
for test in ['test_signal_multirate_complex_taps.f90']:
    f=root/'test'/test
    if not f.exists(): continue
    # Lightweight dependency build: test the numerical module directly; the
    # standalone full-suite runner separately checks the umbrella export.
    local=build/('local_'+test)
    source=f.read_text().replace("    use scifort_signal, only : dp, resample_poly, upfirdn, &\n        resample_status_success, resample_status_invalid_input",
        "    use scifort_kinds, only : dp\n    use scifort_signal_resample, only : resample_poly, upfirdn, &\n        resample_status_success, resample_status_invalid_input")
    local.write_text(source)
    f=local
    exe=build/test.removesuffix('.f90')
    p=subprocess.run(['gfortran',*flags,'-I',str(build),str(f),*objs,'-o',str(exe)],capture_output=True,text=True,cwd=build)
    if p.returncode:
        print(p.stderr,file=sys.stderr);raise SystemExit(1)
    subprocess.run([str(exe)],check=True,cwd=build)
