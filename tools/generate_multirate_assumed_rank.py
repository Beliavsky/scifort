#!/usr/bin/env python3
"""Regenerate portable SELECT RANK pack/unpack cases through Fortran rank 15.

The resampling algorithm is not duplicated; these cases only convert between
assumed-rank arguments and a temporary Fortran-order lane buffer. Idempotent.
"""
from pathlib import Path

file = Path(__file__).resolve().parents[1] / 'src' / 'signal_multirate_nd_into.f90'
s = file.read_text()
for section in ('FLATTEN REAL', 'FLATTEN COMPLEX', 'SCATTER REAL', 'SCATTER COMPLEX'):
    begin = '        ! BEGIN GENERATED ' + section
    end = '        ! END GENERATED ' + section
    statements=[]
    for n in range(1,16):
        if section=='FLATTEN REAL':
            body='v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)'
        elif section=='FLATTEN COMPLEX':
            body='v = reshape(x,[size(x)])'
        elif section=='SCATTER REAL':
            body='y = reshape(real(v,dp),shape(y))'
        else:
            body='y = reshape(v,shape(y))'
        statements += [f'        rank ({n})', f'            {body}']
    start = s.index(begin) + len(begin)
    stop = s.index(end,start)
    s = s[:start] + '\n' + '\n'.join(statements) + '\n        ' + s[stop:]
file.write_text(s)
print('Generated portable assumed-rank real/complex pack/unpack cases for ranks 1-15.')
