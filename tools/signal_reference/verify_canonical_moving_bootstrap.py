#!/usr/bin/env python3
"""Independent NumPy/SciPy moving-block bootstrap of conditional canonical coherence.

The reference independently generates Park-Miller draws, constructs the full
resampled cross-spectral means, applies NumPy linear solves/SVD and uses NumPy
percentile/std computations. The Fortran RNG state and Fortran solver are not
queried or imported by this reference.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import subprocess

import numpy as np

from verify_canonical_jackknife import Case, canonical, integrated, regular


@dataclass
class BootstrapCase:
    case: Case
    block_length: int
    n_resamples: int
    seed: int
    circular: bool


def cases() -> list[BootstrapCase]:
    rng = np.random.default_rng(20261030)
    samples: list[BootstrapCase] = []
    for j in range(56):
        nc = int(rng.integers(3, 7))
        na = int(rng.integers(1, nc))
        nb = int(rng.integers(1, nc-na+1))
        cond = np.arange(na+nb, nc)
        nf = int(rng.integers(2, 6))
        nt = int(rng.integers(7, 16))
        bl = min(nt-1, (1, 2, 3, 4, 6, 8)[j % 6])
        nr = (33, 49, 70)[j % 3]
        ff = np.cumsum(rng.uniform(.14, .36, nf))
        h = rng.normal(size=(nc, nc)) + 1j*rng.normal(size=(nc, nc))
        shared = h @ h.conj().T + np.eye(nc)
        x = np.empty((nf, nt, nc, nc), dtype=complex)
        # Include neighboring-window dependence, while all CSDs remain HPD.
        for f in range(nf):
            previous = np.zeros((nc,nc), dtype=complex)
            for t in range(nt):
                z = rng.normal(size=(nc,nc)) + 1j*rng.normal(size=(nc,nc))
                previous = .65*previous+.35*z
                x[f,t] = ((.86 * shared+.14*previous@previous.conj().T) *
                          (1 + .06*f + .015*t))
        bands = None
        if j % 4 == 0:
            bands = np.array([[ff[0]+.08*(ff[1]-ff[0]),
                               ff[-2]+.77*(ff[-1]-ff[-2])]])
        case = Case(ff,x,np.arange(na),np.arange(na,na+nb),cond,bands,
                    float(rng.uniform(.02,.09)),float(rng.uniform(.04,.31)),
                    1e-6)
        samples.append(BootstrapCase(case, bl, nr, 23+1029*j, bool(j%2)))
    return samples


def resampling_indices(nt: int, length: int, nrep: int, seed: int,
                       circular: bool) -> list[np.ndarray]:
    """Independent Python implementation of documented Park-Miller contract."""
    state = seed
    starts = nt if circular else nt-length+1
    output = []
    for _ in range(nrep):
        indices: list[int] = []
        while len(indices) < nt:
            state = (state*16807) % 2147483647
            start = (state-1) % starts
            for offset in range(min(length,nt-len(indices))):
                indices.append((start+offset)%nt if circular else start+offset)
        output.append(np.asarray(indices,dtype=int))
    return output


def reference(item: BootstrapCase):
    case = item.case
    x = integrated(case)
    nf, nw = x.shape[:2]
    nm = min(len(case.a),len(case.b))
    indices = resampling_indices(nw,item.block_length,item.n_resamples,
                                 item.seed,item.circular)
    estimate = np.empty((nf,nm))
    modes = np.empty((nf,item.n_resamples,nm))
    groups = np.empty((nf,item.n_resamples))
    coherences = np.empty(nf)
    for f in range(nf):
        estimate[f],coherences[f] = canonical(x[f].mean(axis=0),case.a,case.b,
                                               case.c,case.ridge,case.shrink)
        for r,idx in enumerate(indices):
            modes[f,r],groups[f,r] = canonical(x[f,idx].mean(axis=0),case.a,
                                                case.b,case.c,case.ridge,case.shrink)
    group_stats = np.stack((coherences,groups.std(axis=1,ddof=1),
                            groups.mean(axis=1)-coherences,
                            np.percentile(groups,2.5,axis=1,method='linear'),
                            np.percentile(groups,97.5,axis=1,method='linear')),axis=1)
    mode_stats = np.zeros((nf,nm,5))
    validity = np.zeros((nf,nm),dtype=bool)
    for f in range(nf):
        for k in range(nm):
            vals = modes[f,:,k]
            validity[f,k] = regular(estimate[f],k,case.tol) and all(
                regular(modes[f,r],k,case.tol) for r in range(item.n_resamples))
            if validity[f,k]:
                mode_stats[f,k] = (estimate[f,k],vals.std(ddof=1),
                    vals.mean()-estimate[f,k],np.percentile(vals,2.5,method='linear'),
                    np.percentile(vals,97.5,method='linear'))
    return group_stats,mode_stats,validity,estimate


def payload(items: list[BootstrapCase]) -> str:
    lines = [str(len(items))]
    for item in items:
        s = item.case
        nf,nt,nc = s.samples.shape[:3]
        nband = 0 if s.bands is None else len(s.bands)
        lines.append(f'{nf} {nt} {nc} {len(s.a)} {len(s.b)} {len(s.c)} '
                     f'{nband} {item.block_length} {item.n_resamples} {item.seed} '
                     f'{int(item.circular)} {s.ridge:.17g} {s.shrink:.17g} {s.tol:.17g}')
        lines.extend(' '.join(map(str,arr+1)) for arr in (s.a,s.b))
        if len(s.c):
            lines.append(' '.join(map(str,s.c+1)))
        if nband:
            lines.append(' '.join(map(str,s.bands.ravel())))
        lines.append(' '.join(f'{v:.17g}' for v in s.frequencies))
        lines.append(' '.join(map(str,np.arange(nt,dtype=float)+1)))
        for f in range(nf):
            for t in range(nt):
                for row in s.samples[f,t]:
                    lines.append(' '.join(f'({v.real:.17e},{v.imag:.17e})' for v in row))
    return '\n'.join(lines)+'\n'


def main(driver: str) -> None:
    items = cases()
    run = subprocess.run([driver],input=payload(items),capture_output=True,
                         text=True,check=True)
    lines = iter(run.stdout.splitlines())
    maximum_error = 0.0
    valid_modes = 0
    for i,item in enumerate(items):
        nf,nm,nw,bl,nr,ng,wrap = map(int,next(lines).split())
        ref_group,ref_modes,ref_valid,roots = reference(item)
        assert (nf,nm,nw,bl,nr,ng,wrap) == (
            len(ref_group),len(item.case.a) if len(item.case.a) < len(item.case.b)
            else len(item.case.b),item.case.samples.shape[1],item.block_length,
            item.n_resamples,(nw-1)//bl+1,int(item.circular)), i
        for f in range(nf):
            row = next(lines).split()
            assert int(row[0]) == f+1 and int(row[-1]) == 1,(i,f,'group validity')
            actual = np.asarray([float(x) for x in row[1:-1]])
            maximum_error = max(maximum_error,float(np.max(np.abs(actual-ref_group[f]))))
            for k in range(nm):
                row = next(lines).split()
                assert (int(row[0]),int(row[1])) == (f+1,k+1)
                assert bool(int(row[-1])) == ref_valid[f,k],(i,f,k,'mode validity')
                if ref_valid[f,k]:
                    actual = np.asarray([float(x) for x in row[2:-1]])
                    maximum_error = max(maximum_error,float(np.max(
                        np.abs(actual-ref_modes[f,k]))))
                    valid_modes += 1
                else:
                    assert abs(float(row[2])-roots[f,k]) < 1e-10
                    assert all(np.isnan(float(x)) for x in row[3:-1])
    assert maximum_error < 1e-10,maximum_error
    print(f'PASS {len(items)} independent NumPy/SciPy moving-block scenarios, '
          f'{valid_modes} isolated modes, maximum absolute error {maximum_error:.3g}',
          flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    args = parser.parse_args()
    main(args.driver)
