"""Compare streamed multitaper cross-spectra against independent SciPy/NumPy.

Run with: python test/compare_signal_multitaper_stream.py [gfortran]
The reference does not invoke the Fortran batch spectrogram routine.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from compare_signal_multitaper_tf_cross import reference

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / f'{name}.f90' for name in
           ('kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral', 'signal_multitaper')]
DRIVER = ROOT / 'reference' / 'signal_multitaper_stream_reference_driver.f90'


def main() -> None:
    rng = np.random.default_rng(20261020)
    cases = []
    lines = ['80']
    for i in range(80):
        n = int(rng.integers(25, 83))
        seg = int(rng.integers(10, min(n, 35) + 1))
        ov = int(rng.integers(0, seg))
        k = int(rng.integers(2, min(seg, 7) + 1))
        nfft = seg + int(rng.integers(0, 10))
        nw = float(rng.uniform(1.1, min(3.5, seg/2-0.15)))
        rate = float(rng.uniform(1, 101))
        kind = i % 2
        center = i % 3 != 0
        mode = ('equal', 'eigen', 'adaptive')[i % 3]
        x = rng.normal(size=n)
        y = 0.6*x + rng.normal(size=n)
        if kind:
            x = x + 1j*rng.normal(size=n)
            y = y + 1j*rng.normal(size=n)
        if i % 17 == 0:
            x = x*0
        if i % 19 == 0:
            y = y*0
        cases.append((x,y,nw,seg,ov,k,nfft,rate,center,mode))
        lines.append(f'{n} {seg} {ov} {k} {nfft} {kind} {int(center)} '
                     f'{nw:.17g} {rate:.17g} {mode}')
        for a,b in zip(x,y):
            lines.append(f'{a.real:.17g} {a.imag:.17g} {b.real:.17g} {b.imag:.17g}')

    fc = sys.argv[1] if len(sys.argv)>1 else os.environ.get('FC','gfortran')
    with tempfile.TemporaryDirectory(prefix='scifort-stream-reference-') as tmp:
        build = Path(tmp)
        objs = []
        for src in SOURCES:
            obj = build / (src.stem+'.o')
            subprocess.run([fc,'-std=f2018','-O2','-I',tmp,'-J',tmp,
                            '-c',str(src),'-o',str(obj)],cwd=build,check=True)
            objs.append(obj)
        exe = build/'compare_stream'
        subprocess.run([fc,'-std=f2018','-O2','-I',tmp,str(DRIVER),
                        *map(str,objs),'-o',str(exe)],cwd=build,check=True)
        result = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',
                                text=True,capture_output=True,cwd=build,check=True)
    rows = result.stdout.splitlines()
    cursor = 0
    worst = 0.0
    for i,(x,y,nw,seg,ov,k,nfft,rate,center,mode) in enumerate(cases):
        status,nfreq,nframe,nchannel = map(int,rows[cursor].split()); cursor+=1
        assert status == 0 and nchannel == 2, (i,status,nchannel)
        freq = np.array([float(v) for v in rows[cursor:cursor+nfreq]])
        cursor += nfreq
        times = np.array([float(v) for v in rows[cursor:cursor+nframe]])
        cursor += nframe
        values = np.array([[float(v) for v in r.split()]
                           for r in rows[cursor:cursor+nfreq*nframe]])
        cursor += nfreq*nframe
        values = values.reshape(nframe,nfreq,4).transpose(1,0,2)
        ref_freq,ref_times,ref_values = reference(x,y,nw,seg,ov,k,nfft,rate,center,mode)
        np.testing.assert_allclose(freq,ref_freq,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(times,ref_times,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(values,ref_values[...,:4],rtol=5e-8,atol=1e-9,
                                   err_msg=f'streamed spectrum {i}')
        worst = max(worst,float(np.max(abs(values-ref_values[...,:4]))))
    assert cursor == len(rows), (cursor,len(rows))
    print(f'PASS: {len(cases)} streamed spectra against independent SciPy/NumPy reference')
    print(f'Maximum absolute spectral discrepancy: {worst:.3g}')


if __name__ == '__main__':
    main()
