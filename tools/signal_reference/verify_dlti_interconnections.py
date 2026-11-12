#!/usr/bin/env python3
"""Compare SciFort MIMO feedback, Kalman matrices, and Gramians with SciPy.

Usage: python verify_dlti_interconnections.py /path/to/dlti_interconnections_driver
"""
from __future__ import annotations

import sys
import subprocess

import numpy as np
from scipy import linalg


def make_case(rng, i):
    np_, nk = (i % 5, (i * 3) % 4) if i % 17 == 0 else (1 + i % 5, 1 + (i * 3) % 4)
    nu, ny = 1 + (i * 7) % 3, 1 + (i * 11) % 3
    sign = -1 if i % 2 == 0 else 1
    dt = [0.05, 0.125, 1.0, 3.0][i % 4]

    def complex_normal(shape, factor=1):
        return factor * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)

    def transition(n):
        x = complex_normal((n, n))
        if n:
            radius = max(abs(np.linalg.eigvals(x)))
            x *= 0.69 / max(1., radius)
        return x

    ap, ak = transition(np_), transition(nk)
    bp, cp, dp = complex_normal((np_, nu), 0.22), complex_normal((ny, np_), 0.25), complex_normal((ny, nu), 0.16)
    bk, ck, dk = complex_normal((nk, ny), 0.24), complex_normal((nu, nk), 0.21), complex_normal((nu, ny), 0.15)
    return (np_, nk, nu, ny, sign, dt), (ap, bp, cp, dp, ak, bk, ck, dk)


def feedback_states(p, k, sign):
    a, b, c, d = p
    ak, bk, ck, dk = k
    np_, nk = a.shape[0], ak.shape[0]
    nu, ny = b.shape[1], c.shape[0]
    # Independently compute closed-loop matrices by eliminating the algebraic loop.
    e = linalg.solve(np.eye(nu) - sign * dk @ d, np.eye(nu))
    xp = sign * e @ dk @ c
    xc = sign * e @ ck
    aa = np.block([[a + b @ xp, b @ xc],
                   [bk @ (c + d @ xp), ak + bk @ d @ xc]])
    bb = np.concatenate((b @ e, bk @ d @ e), axis=0)
    cc = np.concatenate((c + d @ xp, d @ xc), axis=1)
    dd = d @ e
    return aa, bb, cc, dd


def gramian_reference(a, b, c, d):
    wc = linalg.solve_discrete_lyapunov(a, b @ b.conj().T)
    wo = linalg.solve_discrete_lyapunov(a.conj().T, c.conj().T @ c)
    energy = np.trace(c @ wc @ c.conj().T + d @ d.conj().T).real
    return wc, wo, np.sqrt(max(0., energy))


def read_matrix(lines, pos, shape):
    count = int(np.prod(shape))
    vals = []
    for ln in lines[pos:pos + count]:
        x, y = map(float, ln.split())
        vals.append(x + 1j * y)
    return np.asarray(vals, dtype=complex).reshape(shape, order='F'), pos + count


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    rng = np.random.default_rng(20261025)
    cases = [make_case(rng, i) for i in range(240)]
    tokens = [str(len(cases))]
    for (np_, nk, nu, ny, sign, dt), arrays in cases:
        tokens.append(f'{np_} {nk} {nu} {ny} {sign} {dt:.17g}')
        for mat in arrays:
            for value in mat.ravel(order='F'):
                tokens.append(f'({value.real:.17g},{value.imag:.17g})')
    output = subprocess.run([sys.argv[1]], input='\n'.join(tokens)+'\n', text=True,
                            capture_output=True, check=True)
    lines = output.stdout.splitlines()
    pos = 0
    values = 0
    largest = 0.0
    for i, ((np_, nk, nu, ny, sign, dt), arrays) in enumerate(cases):
        header = [int(x) for x in lines[pos].split()]
        pos += 1
        if header != [i + 1, 0, 0, 0, 0]:
            raise AssertionError((i, 'Fortran statuses', header))
        aa, bb, cc, dd = feedback_states(arrays[:4], arrays[4:], sign)
        n = np_ + nk
        shapes = [(n,n), (n,nu), (ny,n), (ny,nu), (ny+nu,nu+ny), (2,3),
                  (n,n*nu), (n*ny,n), (n,n), (n,n)]
        actual = []
        for shape in shapes:
            matrix, pos = read_matrix(lines, pos, shape)
            actual.append(matrix)
        expected_append_d = linalg.block_diag(arrays[3], arrays[7])
        expected_sel_d = expected_append_d[np.ix_([ny,0], [nu,0,nu])]
        ctr = np.concatenate([np.linalg.matrix_power(aa,j) @ bb for j in range(n)], axis=1) if n else np.zeros((0,0))
        obs = np.concatenate([cc @ np.linalg.matrix_power(aa,j) for j in range(n)], axis=0) if n else np.zeros((0,0))
        if n:
            wc, wo, h2 = gramian_reference(aa,bb,cc,dd)
        else:
            wc, wo = np.zeros((0,0)),np.zeros((0,0))
            h2 = np.linalg.norm(dd, 'fro')
        expected = (aa,bb,cc,dd,expected_append_d,expected_sel_d,ctr,obs,wc,wo)
        for name, x, y in zip(('A','B','C','D','append','select','ctr','obs','wc','wo'),actual,expected):
            err = 0. if x.size == 0 else np.max(abs(x-y))
            bound = 1.e-10*(1.+ (0. if y.size == 0 else np.max(abs(y))))
            if err > bound:
                raise AssertionError((i,name,err,bound,x,y))
            largest = max(largest, err)
            values += x.size
        energy, nit = lines[pos].split()
        pos += 1
        energy = float(energy)
        if abs(energy - h2) > 1.e-10*(1.+h2):
            raise AssertionError((i, 'h2', energy, h2))
        if not 0 <= int(nit) <= 60:
            raise AssertionError((i,'iterations',nit))
        values += 1
        # Independent frequency-domain feedback identity over five complex z.
        a_f, b_f, c_f, d_f = actual[:4]
        for w in np.linspace(0.07,2.8,5):
            z = np.exp(1j*w)
            def eval_ss(abcd):
                a_,b_,c_,d_ = abcd
                return c_ @ linalg.solve(z*np.eye(len(a_))-a_,b_) + d_ if len(a_) else d_
            gp = eval_ss(arrays[:4]); gk=eval_ss(arrays[4:])
            got = eval_ss((a_f,b_f,c_f,d_f))
            want = linalg.solve(np.eye(ny)-sign*gp@gk,gp)
            err=np.max(abs(got-want))
            if err>1.e-9:
                raise AssertionError((i,'feedback frequency',err))
            largest=max(largest,err)
            values += got.size
    if pos != len(lines):
        raise AssertionError(('unconsumed output',len(lines)-pos))
    print(f'PASS: {len(cases)} MIMO cases, {values} values, max absolute error {largest:.3e}')

if __name__ == '__main__':
    main()
