#!/usr/bin/env python3
"""Compare SciFort mixed-integer/constrained DE against SciPy 1.17+.

Comparison is based on objective quality and feasibility, never identical RNG
trajectories. The expected mathematical optima are checked independently.
"""
from __future__ import annotations
import argparse
import math
import pathlib
import subprocess
import tempfile

import numpy as np
from scipy.optimize import differential_evolution, LinearConstraint, NonlinearConstraint

ROOT = pathlib.Path(__file__).resolve().parents[1]
STRATEGIES = ('best1bin', 'rand1bin', 'best2bin', 'currenttobest1bin')
SEEDS = (11, 29, 71, 107)


def problem_data(name):
    if name == 'integer':
        target = np.array([2., -0.75, 3.])
        def fun(x):
            return float(np.sum((x-target)**2))
        return fun, [(-3.6, 4.4), (-3., 4.), (1.2, 5.7)], dict(integrality=[True, False, True]), 0.
    if name == 'linear':
        def fun(x):
            return float((x[0]-.1)**2 + (x[1]-.2)**2)
        return fun, [(-1., 2.), (-1., 2.)], dict(constraints=[
            LinearConstraint([[1., 1.]], lb=[1.5], ub=[np.inf])]), 0.72
    if name == 'circle':
        def fun(x):
            return float((x[0]-1.5)**2+(x[1]-1.3)**2)
        def radius(x):
            return float(np.dot(x, x))
        return fun, [(-1., 2.), (-1., 2.)], dict(constraints=[
            NonlinearConstraint(radius, -np.inf, 1.)]), (math.hypot(1.5,1.3)-1.)**2
    if name == 'combined':
        def fun(x):
            return float((x[0]-1.)**2+(x[1]-.7)**2+x[2]**2)
        def radius(x):
            return float(np.dot(x[:2],x[:2]))
        return fun, [(-2., 2.)]*3, dict(integrality=[True, False, False], constraints=[
            LinearConstraint([[1., 1., 0.]], .5, np.inf),
            NonlinearConstraint(radius, -np.inf, .81)]), 1.
    if name == 'infeasible':
        def fun(x):
            return float((x[0]-.1)**2+(x[1]-.2)**2)
        return fun, [(-1., 1.), (-1., 1.)], dict(constraints=[
            LinearConstraint([[1., 1.]], 10., np.inf)]), None
    raise ValueError(name)


def check_feasibility(problem, x, cv, feasible):
    if problem == 'integer':
        assert abs(x[0]-round(x[0])) <= 1e-12
        assert abs(x[2]-round(x[2])) <= 1e-12
        assert x[0] >= -3 and x[0] <= 4
        assert x[2] >= 2 and x[2] <= 5
    elif problem == 'linear':
        assert x[0]+x[1] >= 1.5-2e-9, (problem,x)
    elif problem == 'circle':
        assert np.dot(x[:2],x[:2]) <= 1.+2e-9, (problem,x)
    elif problem == 'combined':
        assert abs(x[0]-round(x[0])) <= 1e-12
        assert x[0]+x[1] >= .5-2e-9
        assert np.dot(x[:2],x[:2]) <= .81+2e-9
    else:
        assert not feasible
        assert cv >= 8.-1e-8
        return
    assert feasible and cv <= 2e-9, (problem,x,cv,feasible)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fc', default='gfortran')
    parser.add_argument('--mode', choices=('checked','optimized'), default='checked')
    args = parser.parse_args()
    opts = (['-O0','-fcheck=all','-fbacktrace'] if args.mode=='checked'
            else ['-O3','-march=native'])
    nchecked = 0
    worst_diff = 0.
    worst_error = 0.
    with tempfile.TemporaryDirectory(prefix='scifort_de_extended_') as temp:
        exe = pathlib.Path(temp)/'driver'
        subprocess.run([args.fc, '-std=f2018','-Wall','-Wextra','-Werror',
                        '-J',temp,'-I',temp,*opts,'-o',str(exe),
                        str(ROOT/'src/kinds.f90'),
                        str(ROOT/'src/optimize/differential_evolution.f90'),
                        str(ROOT/'tools/reference_de_constraints_driver.f90')],check=True)
        for problem in ('integer','linear','circle','combined','infeasible'):
            fun,bounds,opts_kw,optimal = problem_data(problem)
            for strategy in STRATEGIES:
                for seed in SEEDS:
                    proc=subprocess.run([str(exe),problem,str(seed),strategy],
                                        check=True,capture_output=True,text=True)
                    result=next(line for line in proc.stdout.splitlines() if line.startswith('RESULT'))
                    fields=result.split()
                    status=int(fields[1]); f=float(fields[2]); cv=float(fields[3]);
                    feasible=fields[4]=='T'
                    x=np.array([float(z) for z in fields[5:8]])[:len(bounds)]
                    nfev=int(fields[8]); ncev=int(fields[9]); nit=int(fields[10])
                    check_feasibility(problem,x,cv,feasible)
                    assert status != -1 and nfev >= 0 and nit <= 180
                    assert ncev > 0 or problem == 'integer'
                    if optimal is not None:
                        recalc=fun(x)
                        assert math.isclose(recalc,f,abs_tol=2e-12,rel_tol=2e-12), (problem,recalc,f)
                        quality=f-optimal
                        assert quality > -1e-7, (problem,f,optimal)
                        assert quality < (.06 if problem=='combined' else .035), (problem,f,optimal)
                        worst_error=max(worst_error,quality)
                    else:
                        assert status==4 and nfev==0
                    with np.errstate(over='ignore',invalid='ignore'):
                        ref=differential_evolution(fun,bounds,strategy=strategy,maxiter=(60 if optimal is None else 180),
                             popsize=12,tol=0,atol=0,updating='deferred',polish=False,
                             rng=np.random.default_rng(seed),**opts_kw)
                    if optimal is not None:
                        assert float(ref.fun)-optimal < (.06 if problem=='combined' else .035), (problem,ref.fun)
                        worst_diff=max(worst_diff,abs(f-float(ref.fun)))
                    else:
                        assert not ref.success and ref.maxcv>=8.-1e-7
                    nchecked+=1
    print('EXTENDED SciPy DE COMPARISONS:',nchecked)
    print('WORST ABSOLUTE OBJECTIVE DIFFERENCE:',format(worst_diff,'.8g'))
    print('WORST DISTANCE ABOVE KNOWN OPTIMUM:',format(worst_error,'.8g'))

if __name__=='__main__':
    main()
