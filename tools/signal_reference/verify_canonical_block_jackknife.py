#!/usr/bin/env python3
"""Independent NumPy/SciPy comparison of equal-group delete-block canonical jackknife."""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass

import numpy as np
from verify_canonical_jackknife import Case, canonical, integrated, regular, stdin_payload


@dataclass
class BlockCase:
    case: Case
    length: int


def cases() -> list[BlockCase]:
    rng = np.random.default_rng(20261031)
    samples = []
    for j in range(80):
        nc = int(rng.integers(3, 7))
        na = int(rng.integers(1, nc))
        nb = int(rng.integers(1, nc - na + 1))
        c = np.arange(na + nb, nc)
        nf = int(rng.integers(2, 6))
        length = [1, 2, 3, 4][j % 4]
        nblocks = int(rng.integers(4, 8))
        nt = length * nblocks
        ff = np.cumsum(rng.uniform(0.15, 0.40, nf))
        shared = rng.normal(size=(nc, nc)) + 1j*rng.normal(size=(nc, nc))
        shared = shared @ shared.conj().T + 1.0*np.eye(nc)
        # CSDs with adjacent window dependence but still strictly PD.
        fluctuations = rng.normal(size=(nf, nblocks, nc, nc)) + 1j*rng.normal(
            size=(nf, nblocks, nc, nc))
        x = np.empty((nf, nt, nc, nc), dtype=complex)
        for f in range(nf):
            for t in range(nt):
                g = t // length
                local = rng.normal(size=(nc, nc)) + 1j*rng.normal(size=(nc, nc))
                h = .45 * fluctuations[f, g] + .18 * local
                x[f, t] = (.85 * shared + .15 * h @ h.conj().T) * (1 + .08 * f)
        bands = None
        if j % 3 == 0:
            bands = np.array([[ff[0] + .07*(ff[1]-ff[0]),
                               ff[-2] + .83*(ff[-1]-ff[-2])]])
        samples.append(BlockCase(Case(ff,x,np.arange(na),np.arange(na,na+nb),c,
                                      bands,float(rng.uniform(.005,.08)),
                                      float(rng.uniform(.02,.35))),length))
    return samples


def reference(bc: BlockCase):
    c = bc.case
    x = integrated(c)
    nf, nw = x.shape[:2]
    nm = min(len(c.a),len(c.b))
    ng = nw // bc.length
    estimates = np.zeros((nf,nm))
    se = estimates.copy()
    se_z = estimates.copy()
    lows = estimates.copy()
    highs = estimates.copy()
    bias = estimates.copy()
    valid = np.ones((nf,nm), dtype=bool)
    group = np.zeros(nf)
    se_group = np.zeros(nf)
    for f in range(nf):
        estimates[f], group[f] = canonical(x[f].mean(axis=0),c.a,c.b,c.c,c.ridge,c.shrink)
        loo_roots = np.zeros((ng,nm))
        loo_groups = np.zeros(ng)
        for k in range(ng):
            left = k*bc.length
            right = (k+1)*bc.length
            mean = (x[f].sum(axis=0) - x[f,left:right].sum(axis=0))/(nw-bc.length)
            loo_roots[k], loo_groups[k] = canonical(mean,c.a,c.b,c.c,c.ridge,c.shrink)
        se_group[f] = np.sqrt((ng-1)/ng*np.sum((loo_groups-loo_groups.mean())**2))
        for k in range(nm):
            valid[f,k] = regular(estimates[f],k,c.tol) and all(
                regular(loo_roots[j],k,c.tol) for j in range(ng))
            if not valid[f,k]:
                continue
            roots = loo_roots[:,k]
            se[f,k] = np.sqrt((ng-1)/ng*np.sum((roots-roots.mean())**2))
            bias[f,k] = (ng-1)*(roots.mean()-estimates[f,k])
            zz = np.arctanh(np.sqrt(np.clip(roots,0,1-1e-12)))
            se_z[f,k] = np.sqrt((ng-1)/ng*np.sum((zz-zz.mean())**2))
            center = np.arctanh(np.sqrt(min(1-1e-12,estimates[f,k])))
            lows[f,k] = np.tanh(max(0.0, center-1.959963984540054*se_z[f,k]))**2
            highs[f,k] = np.tanh(center+1.959963984540054*se_z[f,k])**2
    return estimates,se,se_z,lows,highs,bias,valid,group,se_group


def payload(items: list[BlockCase]) -> str:
    lines = [str(len(items))]
    for item in items:
        s = item.case
        nf, nt, nc = s.samples.shape[:3]
        nbands = 0 if s.bands is None else len(s.bands)
        lines.append(f'{nf} {nt} {nc} {len(s.a)} {len(s.b)} {len(s.c)} '
                     f'{nbands} {item.length} {s.ridge:.17g} {s.shrink:.17g} {s.tol:.17g}')
        for arr in (s.a+1, s.b+1):
            lines.append(' '.join(map(str,arr)))
        if len(s.c):
            lines.append(' '.join(map(str,s.c+1)))
        if nbands:
            lines.append(' '.join(map(str,s.bands.ravel())))
        lines.append(' '.join(f'{v:.17g}' for v in s.frequencies))
        lines.append(' '.join(map(str, np.arange(nt,dtype=float)+1)))
        for f in range(nf):
            for t in range(nt):
                for row in s.samples[f,t]:
                    lines.append(' '.join(f'({v.real:.17e},{v.imag:.17e})' for v in row))
    return '\n'.join(lines)+'\n'


def main(binary: str) -> None:
    items = cases()
    run = subprocess.run([binary], input=payload(items),capture_output=True,text=True,check=True)
    lines = iter(run.stdout.splitlines())
    count = 0
    max_error = 0.0
    for i,item in enumerate(items):
        c = item.case
        nf,nm,nw,ng,bl = map(int,next(lines).split())
        ref = reference(item)
        assert (nf,nm,nw,ng,bl) == (ref[0].shape[0],ref[0].shape[1],
                                     c.samples.shape[1],c.samples.shape[1]//item.length,
                                     item.length),i
        for f in range(nf):
            row = next(lines).split()
            assert len(row)==4 and int(row[0])==f+1 and int(row[3])==1
            max_error = max(max_error,abs(float(row[1])-ref[7][f]),
                            abs(float(row[2])-ref[8][f]))
            for k in range(nm):
                row = next(lines).split()
                assert int(row[0])==f+1 and int(row[1])==k+1
                valid = int(row[-1])!=0
                assert valid==ref[6][f,k],(i,f,k,valid,ref[6][f,k])
                if valid:
                    actual = np.asarray([float(v) for v in row[2:-1]])
                    expected = np.asarray([x[f,k] for x in ref[:6]])
                    max_error = max(max_error,float(np.max(np.abs(actual-expected))))
                    count += 1
    assert max_error < 1e-9,max_error
    print(f'PASS {len(items)} independent NumPy/SciPy block scenarios, '
          f'{count} canonical modes, max absolute error {max_error:.3g}',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('driver')
    args = parser.parse_args()
    main(args.driver)
