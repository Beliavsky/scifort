#!/usr/bin/env python3
"""Numerically validate the real optional FITPACK periodic/parametric adapters.

Comparisons with SciPy are classified by whether the adaptive knot sets agree.
When knots differ, check FITPACK's weighted residual and separately reconstruct
SciFort's TCK using SciPy BSpline. Thus different valid solutions are not
incorrectly counted as implementation failures.
"""
import argparse
import subprocess
import numpy as np
from scipy.interpolate import BSpline, splrep, splprep, splev


def run_case(exe, mode, u, data, weight, k, s, q):
    n, ndim = data.shape
    lines = [f'{mode} {n} {ndim} {k} {len(q)} {s:.17g}']
    lines.extend(' '.join(f'{a:.17g}' for a in v) for v in [u, weight])
    lines.extend(' '.join(f'{a:.17g}' for a in data[:, j]) for j in range(ndim))
    lines.append(' '.join(f'{a:.17g}' for a in q))
    proc = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                          capture_output=True, text=True, check=True)
    out = proc.stdout.splitlines()
    if not out or int(out[0]) != 0:
        raise AssertionError(f'FITPACK failure (mode={mode}, k={k}, s={s}): {out[:3]}')
    residual = float(out[1])
    nt, nc, returned_k = map(int, out[2].split())
    assert returned_k == k
    knots = np.fromstring(out[3], sep=' ')
    coefs = np.stack([np.fromstring(line, sep=' ') for line in out[4:4+ndim]], axis=1)
    params = np.fromstring(out[4+ndim], sep=' ')
    val = np.stack([np.fromstring(line, sep=' ') for line in out[5+ndim:5+2*ndim]], axis=1)
    assert len(knots) == nt and coefs.shape == (nc, ndim)
    assert np.array_equal(params, u) or np.allclose(params, u, atol=1e-13)
    deriv = np.stack([np.fromstring(line, sep=' ') for line in out[5+2*ndim:]], axis=1)
    assert val.shape == (len(q), ndim) and deriv.shape == val.shape
    return residual, knots, coefs, val, deriv


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--executable', required=True)
    ap.add_argument('--seed', type=int, default=20261101)
    ap.add_argument('--cases', type=int, default=120)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    same_knots, different_knots = 0, 0
    max_tck_err, max_deriv_err, max_resid_err, max_scipy_err = 0., 0., 0., 0.
    max_scipy_in_domain_err = 0.0
    max_scipy_at_s0_err = 0.0
    counts = [0,0,0]
    for case in range(args.cases):
        mode = case % 3
        k = [1,2,3,4,5][case % 5]
        n = 24 + (case*7) % 48
        ndim = 1 if mode == 0 else (2 + case % 3)
        u = np.r_[0, np.cumsum(rng.uniform(.4,1.3,size=n-1))]
        u /= u[-1]
        w = np.exp(rng.uniform(-.3,.3,size=n))
        theta = 2*np.pi*u
        data = np.empty((n,ndim))
        for j in range(ndim):
            if mode == 1:
                data[:,j] = np.sin(theta*(j+1)*.37) + .25*np.cos(2*theta+.7*j)
            else:
                data[:,j] = np.sin(theta+(j+1)*.5) + .2*np.cos(3*theta+.3*j)
            data[:,j] += .025*rng.standard_normal(n)
        if mode == 2:
            # Exercise SciPy's convention that the closing observation is ignored.
            data[-1] += .13*rng.normal(size=ndim)
        s = [0.0, .005, .05, .25, 1.0][(case//3)%5]
        q = np.linspace(u[0]-.08,u[-1]+.08,79)
        r,t,c,v,dv = run_case(args.executable, mode, u, data, w, k, s, q)
        rebuilt = BSpline(t,c,k,axis=0)(q)
        tck_err = float(np.max(np.abs(rebuilt-v)))
        max_tck_err = max(max_tck_err,tck_err)
        assert tck_err < 5e-10, (case,mode,tck_err)
        derivative_err = float(np.max(np.abs(BSpline(t,c,k,axis=0).derivative(1)(q)-dv)))
        max_deriv_err = max(max_deriv_err, derivative_err)
        assert derivative_err < 5e-8, (case,mode,derivative_err)
        fitted = BSpline(t,c,k,axis=0)(u)
        if mode == 0 or mode == 2:
            # Periodic FITPACK ignores last ordinate; first m-1 are weighted.
            rr = np.sum(((fitted[:-1]-data[:-1])*w[:-1,None])**2)
        else:
            rr = np.sum(((fitted-data)*w[:,None])**2)
        resid_err = abs(rr-r)
        max_resid_err = max(max_resid_err,resid_err)
        assert resid_err < 5e-8*(1+abs(r)), (case,mode,rr,r)
        assert r <= s*(1.0015)+1e-8, (case,mode,s,r)
        if s == 0:
            assert r < 1e-8, (case,mode,r)
        if mode == 0:
            fit_tck = splrep(u,data[:,0],w=w,k=k,s=s,per=1)
            t_ref = fit_tck[0]
            v_ref = np.asarray(splev(q,fit_tck))[:,None]
        else:
            fit_tck, u_ref = splprep(data.T, u=u, w=w, k=k, s=s,per=int(mode==2))
            t_ref = fit_tck[0]
            v_ref = np.stack(splev(q,fit_tck),axis=1)
        if len(t_ref)==len(t) and np.allclose(t_ref,t,atol=2e-10,rtol=2e-10):
            same_knots += 1
            e = float(np.max(np.abs(v_ref-v)))
            max_scipy_err = max(max_scipy_err,e)
            inside = (q >= u[0]) & (q <= u[-1])
            max_scipy_in_domain_err = max(max_scipy_in_domain_err,
                                          float(np.max(np.abs(v_ref[inside]-v[inside]))))
            if s == 0.0:
                max_scipy_at_s0_err = max(max_scipy_at_s0_err, e)
                assert e < 1e-8, (case,mode,k,s,e)
        else:
            different_knots += 1
        counts[mode]+=1
    print('PASS:',args.cases,'comparisons (scalar periodic, open parametric, closed parametric)',counts)
    print('SciPy agreeing knots:',same_knots,'different knots:',different_knots)
    print('Maximum BSpline TCK discrepancy:',max_tck_err)
    print('Maximum independently computed derivative discrepancy:',max_deriv_err)
    print('Maximum independently computed weighted residual discrepancy:',max_resid_err)
    print('Maximum SciPy discrepancy at s=0 with agreeing knots:',max_scipy_at_s0_err)
    print('Maximum SciPy discrepancy within fit domain (s may be >0):',max_scipy_in_domain_err)
    print('Maximum SciPy discrepancy including extrapolation (s may be >0):',max_scipy_err)


if __name__ == '__main__':
    main()
