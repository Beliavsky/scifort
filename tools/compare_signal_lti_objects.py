#!/usr/bin/env python3
"""Compare SciFort balanced reduction and continuous LTI with SciPy/NumPy.

Pass the compiled tools/signal_lti_objects_reference_driver executable as argv[1].
No additional Python dependencies beyond NumPy and SciPy are required.
"""
import sys
import subprocess
import numpy as np
from scipy import linalg

rng = np.random.default_rng(20261025)


def complex_matrix(n, m, scale=1.0):
    return (rng.normal(size=(n, m)) + 1j*rng.normal(size=(n, m)))*scale


def pack(matrix):
    return [f"{v.real:.18e} {v.imag:.18e}" for v in matrix.flat]


def random_system(mode, index):
    n = 1 + index % 5
    nu = 1 + (index//2) % 3
    ny = 1 + (index//3) % 3
    q, _ = np.linalg.qr(complex_matrix(n,n))
    if mode == 1:
        poles = rng.uniform(-0.65,0.65,n) + 1j*rng.uniform(-0.3,0.3,n)
    else:
        poles = -rng.uniform(0.4,2.6,n) + 1j*rng.uniform(-1.0,1.0,n)
    a = q @ np.diag(poles) @ q.conj().T
    b = complex_matrix(n,nu,scale=0.6)
    c = complex_matrix(ny,n,scale=0.5)
    d = complex_matrix(ny,nu,scale=0.1)
    order = index%(n+1)
    dt = 0.04+index%6*0.021
    return (mode,n,nu,ny,order,dt,a,b,c,d)


class Reader:
    def __init__(self, payload):
        self.tokens = iter(payload.split())
    def integer(self):
        return int(next(self.tokens))
    def real(self):
        return float(next(self.tokens))
    def matrix(self, n, m):
        return np.array([complex(self.real(),self.real()) for _ in range(n*m)],
                        dtype=np.complex128).reshape(n,m)


def transfer(a,b,c,d,z):
    return c @ linalg.solve(z*np.eye(len(a))-a,b) + d


def independent_hsv(a,b,c):
    wc = linalg.solve_discrete_lyapunov(a,b@b.conj().T)
    wo = linalg.solve_discrete_lyapunov(a.conj().T,c.conj().T@c)
    ec,uc = linalg.eigh(wc)
    eo,uo = linalg.eigh(wo)
    lc = uc @ np.diag(np.sqrt(np.maximum(0.0,ec)))
    lo = uo @ np.diag(np.sqrt(np.maximum(0.0,eo)))
    return linalg.svdvals(lo.conj().T@lc)


def check_discrete(case, reader):
    mode,n,nu,ny,order,dt,a,b,c,d = case
    status = reader.integer()
    assert status == 0, (mode,n,nu,ny,order,status)
    hsv = np.array([reader.real() for _ in range(n)])
    bound = reader.real()
    ar = reader.matrix(order,order)
    br = reader.matrix(order,nu)
    cr = reader.matrix(ny,order)
    dr = reader.matrix(ny,nu)
    expected = independent_hsv(a,b,c)
    np.testing.assert_allclose(hsv,expected,rtol=2e-7,atol=3e-10)
    np.testing.assert_allclose(bound,2*sum(expected[order:]),rtol=2e-7,atol=3e-10)
    np.testing.assert_allclose(dr,d,atol=2e-12,rtol=0)
    largest_error = 0.0
    for w in np.linspace(0,np.pi,21):
        h = transfer(a,b,c,d,np.exp(1j*w))
        hr = transfer(ar,br,cr,dr,np.exp(1j*w))
        largest_error = max(largest_error,abs(h-hr).max())
    assert largest_error <= bound+5e-7, (n,order,bound,largest_error)
    if order == n:
        assert largest_error < 1e-7, (n,order,largest_error)
    return max(abs(hsv-expected)),largest_error


def check_analog(case,reader):
    mode,n,nu,ny,order,dt,a,b,c,d = case
    assert reader.integer() == 0
    assert reader.integer() == 0
    w = [0.0,0.3,0.9,2.4,5.0]
    differences = []
    for f in w:
        received = reader.matrix(ny,nu)
        expected = transfer(a,b,c,d,1j*f)
        differences.append(abs(received-expected).max())
    assert reader.integer() == 0
    actual_y = reader.matrix(10,ny)
    u = np.array([[.5*complex(np.cos(.3*(i+1)+.2*(j+1)),
                  np.sin(.2*(i+1)-.12*(j+1))) for j in range(nu)]
                  for i in range(10)])
    aug = np.zeros((n+nu,n+nu),dtype=complex)
    aug[:n,:n] = a
    aug[:n,n:] = b
    expm = linalg.expm(aug*dt)
    ad = expm[:n,:n]
    bd = expm[:n,n:]
    x = np.zeros(n,dtype=complex)
    reference_y = np.empty((10,ny),dtype=complex)
    for i in range(10):
        reference_y[i,:] = c@x+d@u[i,:]
        if i < 9:
            x = ad@x+bd@u[i,:]
    differences.append(abs(reference_y-actual_y).max())
    assert reader.integer() == 0
    found_ad = reader.matrix(n,n)
    found_bd = reader.matrix(n,nu)
    found_cd = reader.matrix(ny,n)
    found_dd = reader.matrix(ny,nu)
    differences.append(abs(found_ad-ad).max())
    differences.append(abs(found_bd-bd).max())
    differences.append(abs(found_cd-c).max())
    differences.append(abs(found_dd-d).max())
    assert max(differences) < 3e-9, (n,nu,ny,differences)
    return max(differences)


def main():
    cases = [random_system(mode,i) for mode in (1,2) for i in range(100)]
    data = [str(len(cases))]
    for mode,n,nu,ny,order,dt,a,b,c,d in cases:
        data.append(f"{mode} {n} {nu} {ny} {order} {dt:.18e}")
        data.extend(pack(a))
        data.extend(pack(b))
        data.extend(pack(c))
        data.extend(pack(d))
    finished = subprocess.run([sys.argv[1]],input='\n'.join(data)+'\n',
                              capture_output=True,text=True,check=True)
    reader = Reader(finished.stdout)
    discretes = []
    analogs = []
    for case in cases:
        if case[0] == 1:
            discretes.append(check_discrete(case,reader))
        else:
            analogs.append(check_analog(case,reader))
    assert next(reader.tokens,None) is None
    print(f"PASS: 100 balanced reductions, max Hankel error "
          f"{max(d[0] for d in discretes):.3e}; "
          f"100 continuous LTI cases, max analog error {max(analogs):.3e}")


if __name__ == '__main__':
    main()
