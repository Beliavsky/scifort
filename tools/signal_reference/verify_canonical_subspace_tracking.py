#!/usr/bin/env python3
"""Independent complex-SVD, principal-angle, and assignment reference tests."""
from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
import subprocess

import numpy as np
from scipy.optimize import linear_sum_assignment


@dataclass
class Scenario:
    freq: np.ndarray
    values: np.ndarray
    a: np.ndarray
    b: np.ndarray
    valid: np.ndarray
    gap: float
    minimum: float


def unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))[None, :]


def near_identity(rng: np.random.Generator, n: int, amplitude: float) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    hermitian = (z + z.conj().T) / 2
    eig, basis = np.linalg.eigh(hermitian)
    return (basis * np.exp(1j * amplitude * eig)[None, :]) @ basis.conj().T


def blocks(values: np.ndarray, gap: float) -> list[tuple[int, int]]:
    output = []
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values):
            v, w = values[j - 1:j + 1]
            if not (np.isfinite(v) and np.isfinite(w) and min(v, w) >= 0):
                break
            if abs(v - w) > gap * max(v, w, 0.05):
                break
            j += 1
        output.append((i, j))
        i = j
    return output


def stacked(a: np.ndarray, b: np.ndarray) -> np.ndarray | None:
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return None
    sa = np.sqrt(np.sum(np.abs(a) ** 2) / a.shape[1])
    sb = np.sqrt(np.sum(np.abs(b) ** 2) / b.shape[1])
    if min(sa, sb) <= np.finfo(float).tiny:
        return None
    return np.concatenate((a / (np.sqrt(2) * sa), b / (np.sqrt(2) * sb)))


def comparison(pa, pb, ca, cb):
    xp, xc = stacked(pa, pb), stacked(ca, cb)
    if xp is None or xc is None:
        return None
    qp, rp = np.linalg.qr(xp, mode='reduced')
    qc, rc = np.linalg.qr(xc, mode='reduced')
    if min(np.abs(np.diag(rp)).min(), np.abs(np.diag(rc)).min()) < 1.e-11:
        return None
    cos = np.clip(np.linalg.svd(qp.conj().T @ qc, compute_uv=False), 0, 1)
    u, singular, vh = np.linalg.svd(xp.conj().T @ xc)
    if singular.min() <= 1.e-12:
        return None
    rot = vh.conj().T @ u.conj().T
    return cos, rot


def reference(s: Scenario):
    nf, nm = s.values.shape
    a, b = s.a.copy(), s.b.copy()
    ids = np.zeros((nf, nm), dtype=int)
    first = np.zeros_like(ids)
    size = np.zeros_like(ids)
    matched = np.zeros_like(ids, dtype=bool)
    new = np.zeros_like(ids, dtype=bool)
    cosines = np.full((nf, nm), np.nan)
    similarity = np.full((nf, nm), np.nan)
    prev = []
    next_id = 0
    for f in range(nf):
        if not s.valid[f]:
            prev = []
            continue
        curr = blocks(s.values[f], s.gap)
        curr_good = {}
        for j, (lo, hi) in enumerate(curr):
            first[f, lo:hi] = lo + 1
            size[f, lo:hi] = hi - lo
            curr_good[j] = (np.all(np.isfinite(s.values[f, lo:hi]))
                            and s.values[f, lo:hi].min() > 256 * np.finfo(float).eps
                            and stacked(a[f, :, lo:hi], b[f, :, lo:hi]) is not None)
        links = {}
        if prev and np.isfinite(s.freq[f]) and np.isfinite(s.freq[f-1]) \
                and s.freq[f] > s.freq[f-1]:
            nmax = max(len(prev), len(curr))
            rewards = np.zeros((nmax, nmax))
            for i, (plo, phi, good) in enumerate(prev):
                if not good:
                    continue
                for j, (clo, chi) in enumerate(curr):
                    if not curr_good[j] or phi - plo != chi - clo:
                        continue
                    candidate = comparison(a[f-1, :, plo:phi],b[f-1, :, plo:phi],
                                           a[f, :, clo:chi],b[f, :, clo:chi])
                    if candidate is None:
                        continue
                    vmin = candidate[0].min()
                    if vmin > 64*np.finfo(float).eps and vmin >= s.minimum:
                        rewards[i,j] = 1+vmin/(nmax+1)
            ii,jj = linear_sum_assignment(-rewards)
            for i,j in zip(ii,jj):
                if i < len(prev) and j < len(curr) and rewards[i,j] > 0:
                    links[j] = i
        for j,(lo,hi) in enumerate(curr):
            if j in links:
                plo, phi, _ = prev[links[j]]
                val,rot = comparison(a[f-1,:,plo:phi],b[f-1,:,plo:phi],
                                     a[f,:,lo:hi],b[f,:,lo:hi])
                ids[f,lo:hi] = ids[f-1,plo]
                aa = a[f,:,lo:hi].copy()
                bb = b[f,:,lo:hi].copy()
                a[f,:,lo:hi] = aa @ rot
                b[f,:,lo:hi] = bb @ rot
                cosines[f,lo:hi] = val
                similarity[f,lo:hi] = val[-1]
                matched[f,lo:hi] = True
            else:
                next_id += 1
                ids[f,lo:hi] = next_id
                new[f,lo:hi] = True
        prev = [(lo,hi,curr_good[j]) for j,(lo,hi) in enumerate(curr)]
    return ids, first, size, matched, new, cosines, similarity, a, b


def cases(seed=20261030, ncases=100):
    rng = np.random.default_rng(seed)
    output = []
    for case in range(ncases):
        nm = int(rng.integers(2,6))
        na, nb = nm + int(rng.integers(0,3)), nm + int(rng.integers(0,3))
        nf = int(rng.integers(4,10))
        gap = float(rng.choice([0, 0.001, 0.01, 0.05]))
        minimum = float(rng.choice([0.2, 0.5, 0.7, 0.9]))
        freq = np.arange(nf, dtype=float)
        valid = np.ones(nf, dtype=bool)
        if case % 9 == 0:
            freq[nf//2] = -3.0
        if case % 11 == 0:
            valid[nf//2] = False
        values = np.zeros((nf,nm))
        a = np.zeros((nf,na,nm),dtype=complex)
        b = np.zeros((nf,nb,nm),dtype=complex)
        aa = unitary(rng,na)[:,:nm]
        bb = unitary(rng,nb)[:,:nm]
        for f in range(nf):
            # Vary clusters independently. Keep descending ranks, but sometimes
            # exchange physically continuous canonical vector labels.
            numbers = 0.92 - 0.14*np.arange(nm,dtype=float)
            if case % 17 == 0 and f in (nf//2,nf//2+1):
                numbers[-1] = 0.0
            else:
                k = 0
                while k < nm-1:
                    if rng.random() < 0.48:
                        numbers[k+1] = numbers[k]-float(rng.uniform(0,0.2))*gap
                        k += 1
                    k += 1
            values[f] = numbers
            if f:
                aa = near_identity(rng,na,0.06) @ aa
                bb = near_identity(rng,nb,0.06) @ bb
                if nm > 1 and case % 3 == 0:
                    index = int(rng.integers(0,nm-1))
                    aa[:,[index,index+1]] = aa[:,[index+1,index]]
                    bb[:,[index,index+1]] = bb[:,[index+1,index]]
            aa_new,bb_new = aa.copy(),bb.copy()
            for lo,hi in blocks(values[f],gap):
                r = unitary(rng,hi-lo)
                aa_new[:,lo:hi] = aa[:,lo:hi] @ r
                bb_new[:,lo:hi] = bb[:,lo:hi] @ r
            a[f],b[f] = aa_new,bb_new
        output.append(Scenario(freq,values,a,b,valid,gap,minimum))
    return output


def protocol(scenarios):
    chunks = [str(len(scenarios))]
    for s in scenarios:
        nf,nm = s.values.shape
        na,nb = s.a.shape[1],s.b.shape[1]
        chunks.append(f'{nf} {nm} {na} {nb} {s.gap:.17g} {s.minimum:.17g}')
        chunks.append(' '.join(f'{x:.17g}' for x in s.freq))
        chunks.append(' '.join(str(int(x)) for x in s.valid))
        for f in range(nf):
            chunks.append(' '.join(f'{x:.17g}' for x in s.values[f]))
            for k in range(nm):
                for matrix in (s.a,s.b):
                    for v in matrix[f,:,k]:
                        chunks.append(f'({v.real:.17g},{v.imag:.17g})')
    return '\n'.join(chunks)+'\n'


def run(binary: str, scenarios):
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
    proc = subprocess.run([binary],input=protocol(scenarios),text=True,check=True,
                          capture_output=True,env=env)
    lines = iter(proc.stdout.splitlines())
    maxerr = 0.0
    observed = 0
    for case,s in enumerate(scenarios):
        nf,nm = map(int,next(lines).split())
        assert (nf,nm) == s.values.shape
        truth = reference(s)
        na,nb = s.a.shape[1],s.b.shape[1]
        for f in range(nf):
            for k in range(nm):
                fields = next(lines).split()
                ii = np.asarray(fields[:7],dtype=int)
                assert (ii[0],ii[1]) == (f+1,k+1), (case, f, k, ii)
                expected_ints = [truth[0][f,k],truth[1][f,k],truth[2][f,k],
                                 int(truth[3][f,k]),int(truth[4][f,k])]
                if ii[2:].tolist() != expected_ints:
                    raise AssertionError(('block assignment',case,f,k,ii[2:].tolist(),expected_ints))
                realparts = np.asarray(fields[7:],dtype=float)
                numbers = [truth[5][f,k],truth[6][f,k]]
                for got,want in zip(realparts[:2],numbers):
                    if np.isnan(want):
                        assert np.isnan(got), (case,f,k,got,want)
                    else:
                        err = abs(got-want)
                        maxerr = max(maxerr,err)
                        assert err < 1.0e-10, ('principal cosine',case,f,k,err)
                expected_a = truth[7][f,:,k]
                expected_b = truth[8][f,:,k]
                result_complex = realparts[2::2] + 1j*realparts[3::2]
                for got,want in ((result_complex[:na],expected_a),(result_complex[na:],expected_b)):
                    err = np.max(np.abs(got-want))
                    maxerr = max(maxerr,err)
                    assert err < 1.0e-9, ('unitary Procrustes',case,f,k,err)
                observed += 1
    try:
        next(lines)
        raise AssertionError('unexpected driver output')
    except StopIteration:
        pass
    print(f'PASS {len(scenarios)} complex-SVD/assignment scenarios, '
          f'{observed} modes, max absolute error {maxerr:.3g}',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    args=parser.parse_args()
    run(args.driver,cases())
