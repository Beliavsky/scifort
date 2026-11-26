# SPDX-License-Identifier: MIT
"""Explicit MRG32k3a stream with caller-owned, snapshot-capable serialized state."""
import ctypes as C
import numpy as np
from ._runtime import StateBuffer


class RandomState(StateBuffer):
    """Native SciFort stream; reproducible within the documented state/algorithm contract."""
    def __init__(self, owner, seed=12345):
        if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)):
            raise TypeError('seed must be a signed int64 integer')
        if not -(2**63) <= int(seed) <= 2**63-1:
            raise ValueError('seed must fit signed int64')
        self._owner = owner
        super().__init__(owner._library, np.zeros(6,dtype=np.int64))
        self._seed = self._library.scifort_rng_seed_i64
        self._seed.argtypes = [C.c_int64,C.c_void_p,C.POINTER(C.c_int)]
        self._seed.restype = None
        self._uniform = self._library.scifort_rng_uniform_f64
        self._uniform.argtypes = [C.c_int,C.c_void_p,C.c_void_p,C.POINTER(C.c_int)]
        self._uniform.restype = None
        status = C.c_int(-1)
        self._seed(int(seed),self._state.ctypes.data_as(C.c_void_p),C.byref(status))
        if status.value:
            raise RuntimeError(f'RNG seed failed with status {status.value}')

    def set_state(self, values):
        """Restore a validated six-component snapshot; invalid input preserves the live stream."""
        candidate = StateBuffer(self._library,values)._state
        if candidate.shape != (6,):
            raise ValueError('state must have six components')
        scratch = np.empty(1,dtype=np.float64)
        status = C.c_int(-1)
        self._uniform(0,candidate.ctypes.data_as(C.c_void_p),scratch.ctypes.data_as(C.c_void_p),C.byref(status))
        if status.value:
            raise ValueError('invalid MRG32k3a state components')
        self._state[:] = candidate

    def uniform(self, n, *, out=None):
        """Draw n values in (0,1); an empty call validates state and does not advance it."""
        if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)):
            raise TypeError('n must be an integer')
        n = int(n)
        if not 0 <= n <= np.iinfo(np.intc).max:
            raise ValueError('n must be a nonnegative C-int extent')
        if out is None:
            out = np.empty(n,dtype=np.float64)
        if not isinstance(out,np.ndarray) or out.dtype != np.dtype(np.float64):
            raise TypeError('out must be a native float64 NumPy array')
        if (out.shape != (n,) or not out.flags.c_contiguous or not out.flags.aligned
                or not out.flags.writeable or np.may_share_memory(out,self._state)):
            raise ValueError('out must be aligned, writable, contiguous, disjoint and have shape (n,)')
        scratch = out if n else np.empty(1,dtype=np.float64)
        status = C.c_int(-1)
        self._uniform(n,self._state.ctypes.data_as(C.c_void_p),scratch.ctypes.data_as(C.c_void_p),C.byref(status))
        if status.value:
            raise ValueError(f'RNG call failed with status {status.value}')
        return out
