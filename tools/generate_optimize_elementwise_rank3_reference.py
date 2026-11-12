"""Generate SciPy 1.17 rank-three independent elementwise references."""
from pathlib import Path

import numpy as np
from scipy.optimize import elementwise


def format_array(name: str, value: np.ndarray) -> str:
    """Write a rank-3 Fortran parameter array in column-major order."""
    data = np.asarray(value)
    assert data.ndim == 3
    dims = ','.join(map(str, data.shape))
    is_int = data.dtype.kind in 'iu'
    fmt = (lambda n: str(int(n))) if is_int else (lambda n: f'{float(n):.16e}_dp')
    flat = [fmt(n) for n in data.flatten(order='F')]
    lines = [', '.join(flat[i:i+2]) for i in range(0, len(flat), 2)]
    dtype = 'integer' if is_int else 'real(dp)'
    return f'    {dtype}, parameter :: {name}({dims}) = &\n' + \
        '        reshape([ &\n            ' + ', &\n            '.join(lines) + f'], [{dims}])\n'


def generate() -> None:
    """Compute starting grids and scalar-compatible batched SciPy results."""
    f = lambda x: x*x - 4.0
    g = lambda x: (x-1.5)**2
    left = np.stack([
        [[-4., 0., 3.], [-4., 0., 3.]],
        [[-3., -0.5, 3.2], [-4., 0.2, 3.1]]
    ], axis=2)
    right = np.stack([
        [[-3., 1., 3.5], [-3., 1., 3.5]],
        [[-2.5, 1., 3.7], [-3., 1., 3.8]]
    ], axis=2)
    lower = np.stack([
        [[-10., 0., 1.], [-10., 0., 1.]],
        [[-10., -1., 1.], [-10., 0., 1.]]
    ], axis=2)
    upper = np.stack([
        [[-1., 10., 4.], [-1., 10., 4.]],
        [[-1., 10., 4.], [-1., 10., 4.]]
    ], axis=2)
    root_pairs = np.empty((2,2,3,2))
    root_pairs[:,:,:,0] = np.array([
        [[-4., 1., 3.], [4., -4., 1.]],
        [[-1., 4., 4.], [1., -1., 4.]]
    ])
    root_pairs[:,:,:,1] = np.array([
        [[1., -4., 3.], [-4., 3., -4.]],
        [[4., -1., 4.], [-1., 4., -1.]]
    ])
    middle = np.stack([
        [[-3., 0.5, 1.5], [-3., 0.5, 1.5]],
        [[3.5, -0.2, 2.0], [4., 0.2, 1.0]],
    ], axis=2)
    minimum_triples = np.empty((3,2,3,2))
    minimum_triples[:,:,:,0] = np.array([
        [[0., -3., 3.], [3., 0., -3.]],
        [[1., 0., 4.], [1., 1., 0.]],
        [[3., 4., 5.], [0., 3., 4.]],
    ])
    minimum_triples[:,:,:,1] = np.array([
        [[0., -2., 3.], [4., 0., -2.]],
        [[1., 0., 4.], [1., 1., 0.]],
        [[3., 4., 5.], [0., 3., 4.]],
    ])
    br = elementwise.bracket_root(f, left, xr0=right, xmin=lower,
                                  xmax=upper, factor=2.0, maxiter=15)
    rr = elementwise.find_root(f, (root_pairs[0], root_pairs[1]),
                               tolerances={'xatol':1.e-12})
    bm = elementwise.bracket_minimum(g, middle, xl0=middle-0.3,
                                     xr0=middle+0.3, xmin=-10., xmax=10.)
    rm = elementwise.find_minimum(g, tuple(minimum_triples),
                                  tolerances={'xatol':1.e-9})
    parts = [
        '! Generated from SciPy 1.17.0 scipy.optimize.elementwise; do not edit.\n',
        '! SPDX-License-Identifier: MIT\n',
        'module optimize_elementwise_rank3_reference\n',
        '    use scifort_kinds, only : dp\n',
        '    implicit none\n',
    ]
    for label,data in [('root_left',left),('root_right',right),('root_lower',lower),
                       ('root_upper',upper),('minimum_middle',middle)]:
        parts.append(format_array(label+'_ref', data))
    for label,data in [('root_pairs',root_pairs),('minimum_triples',minimum_triples)]:
        # A rank-4 input bracket is stored as a 1-D parameter list and reshaped
        # in the Fortran regression program.
        shape = ','.join(map(str, data.shape))
        values = [f'{float(n):.16e}_dp' for n in data.flatten(order='F')]
        lines = [', '.join(values[i:i+2]) for i in range(0,len(values),2)]
        parts.append('    real(dp), parameter :: '+label+'_ref('+shape+') = &\n'
                     '        reshape([ &\n            ' + ', &\n            '.join(lines)
                     + '], [' + shape + '])\n')
    for prefix,res,data in [
        ('br_root',br,{'xl':br.bracket[0],'xr':br.bracket[1]}),
        ('find_root',rr,{'x':np.where(rr.success,rr.x,0.)}),
        ('br_min',bm,{'xm':bm.bracket[1]}),
        ('find_min',rm,{'x':np.where(rm.success,rm.x,0.)}),
    ]:
        for key in ('status','nit','nfev'):
            parts.append(format_array(prefix+'_'+key+'_ref',getattr(res,key)))
        for key,v in data.items():
            parts.append(format_array(prefix+'_'+key+'_ref',v))
    parts.append('end module optimize_elementwise_rank3_reference\n')
    output = Path(__file__).resolve().parents[1]/'test'/'optimize_elementwise_rank3_reference.f90'
    output.write_text(''.join(parts),encoding='ascii')
    print(f'Wrote {output}; rank-three grid {left.shape}; scipy {np.__version__}')

if __name__ == '__main__':
    generate()
