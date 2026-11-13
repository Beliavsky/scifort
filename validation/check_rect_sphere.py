#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare native rectangular spherical interpolation and pole options with SciPy."""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.interpolate import RectSphereBivariateSpline
from scipy.interpolate import BSpline


def native(exe,u,v,r,continuity,poles,exact,flat,qt,qp,s=0.0,metadata=False):
    lines = [f"{len(u)} {len(v)} {len(qt)}"]
    for row in (u,v,*r):
        lines.append(" ".join(f"{x:.17g}" for x in row))
    lines.append(" ".join(str(int(x)) for x in (*continuity,*exact,*flat)))
    lines.append(" ".join("NaN" if x is None else f"{x:.17g}" for x in poles))
    lines += [f"{t:.17g} {p:.17g}" for t,p in zip(qt,qp)]
    lines.append(f"{s:.17g}")
    result = subprocess.run([str(exe.resolve())],input="\n".join(lines)+"\n",
                            text=True,capture_output=True,check=True)
    out = result.stdout.splitlines()
    assert int(out[0].split()[0]) == 0, result.stdout
    knots = [np.fromstring(line,sep=" ") for line in out[1:3]]
    coef = np.fromstring(out[3],sep=" ")
    values = np.array([np.fromstring(line,sep=" ") for line in out[4:]])
    assert values.shape == (len(qt),6)
    if metadata:
        return knots,coef,values,np.fromstring(out[0],sep=" ")
    return knots,coef,values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable",type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261114)
    configs = [((False,False),(None,None),(False,False),(False,False)),
        ((True,True),(None,None),(False,False),(False,False)),
        ((True,False),(None,None),(False,False),(False,False)),
        ((False,True),(None,None),(False,False),(False,False)),
        ((True,True),(None,None),(False,False),(True,True)),
        ((True,False),(1.2,None),(True,False),(True,False)),
        ((False,True),(None,-0.8),(False,True),(False,True)),
        ((True,True),(1.2,-0.8),(False,False),(False,False))]
    max_error = max_coef = max_node = 0.0
    count = comparisons = 0
    for mu,mv,origin in ((4,4,0.0),(7,9,-1.2),(8,12,0.65)):
        u = np.linspace(0.15,np.pi-0.2,mu)
        v = origin+np.arange(mv)*2*np.pi/mv
        if mu == 8:
            u = np.linspace(0.2,np.pi-0.25,mu)**1.05
            u *= (np.pi-0.15)/u[-1]
            v = origin+np.r_[0,np.sort(rng.uniform(0.2,2*np.pi-0.2,mv-1))]
        nodes_t,nodes_p = np.meshgrid(u,v,indexing="ij")
        qt = np.r_[nodes_t.ravel(),0,np.pi,rng.uniform(0,np.pi,25)]
        qp = np.r_[nodes_p.ravel(),origin,origin+2*np.pi,rng.uniform(origin,origin+2*np.pi,25)]
        for kind in ("constant","harmonic"):
            r = np.full((mu,mv),2.5) if kind == "constant" else (
                np.cos(nodes_t)+0.7*np.sin(nodes_t)*np.cos(nodes_p)
                +0.2*np.sin(nodes_t)**2*np.sin(2*nodes_p))
            for continuity,poles,exact,flat in configs:
                if kind == "constant":
                    poles = tuple(None if x is None else 2.5 for x in poles)
                ref = RectSphereBivariateSpline(u,v,r,pole_continuity=continuity,
                    pole_values=poles,pole_exact=exact,pole_flat=flat)
                knots,coef,values = native(args.executable,u,v,r,continuity,poles,exact,flat,qt,qp)
                for a,b in zip(knots,ref.get_knots()):
                    assert np.allclose(a,b,rtol=0,atol=3e-15)
                ce = np.max(abs(coef-ref.get_coeffs()))/(1+np.max(abs(ref.get_coeffs())))
                expected = np.column_stack([ref.ev(qt,qp,dtheta=a,dphi=b)
                    for a,b in ((0,0),(1,0),(0,1),(2,0),(0,2),(1,1))])
                error = np.max(abs(values-expected)/(1+abs(expected)))
                node_error = np.max(abs(values[:mu*mv,0]-r.ravel()))
                assert ce < 3e-9, (mu,mv,kind,continuity,poles,flat,"coef",ce)
                assert error < 5e-8, (mu,mv,kind,continuity,poles,flat,"value",error)
                assert node_error < 3e-11*(1+np.max(abs(r))), node_error
                max_error = max(max_error,error)
                max_coef = max(max_coef,ce)
                max_node = max(max_node,node_error)
                count += 1
                comparisons += values.size
    print(f"PASS: {count} grids, {comparisons} value/derivative references; SciPy {scipy.__version__}")
    print(f"Maximum normalized errors: value {max_error:.6g}, coefficients {max_coef:.6g}; "
          f"maximum absolute node error {max_node:.6g}")
    # Modern FITPACK includes supplied pole values in the optimization step
    # range. SciPy 1.15.3 omits a low south pole for constant grid data, so it
    # leaves the unknown north pole at its row-mean seed. Verify the modern
    # result independently by constrained cubic least squares instead.
    u = np.linspace(0.15,np.pi-0.2,4)
    v = np.arange(4)*2*np.pi/4
    r = np.full((4,4),2.5)
    _,_,values = native(args.executable,u,v,r,(False,True),(None,-0.8),
        (False,True),(False,True),np.array([0,np.pi]),np.array([0,0]))
    basis = BSpline.design_matrix(u,np.r_[[0]*4,[np.pi]*4],3).toarray()
    expected = np.linalg.lstsq(basis[:,:2],2.5+0.8*basis[:,2:].sum(axis=1),rcond=None)[0][0]
    assert abs(values[0,0]-expected) < 3e-12
    assert abs(values[1,0]+0.8) < 3e-13
    print(f"PASS: pole-estimation difference independently verified; north={expected:.16g}")


if __name__ == "__main__":
    main()
