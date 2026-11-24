# SPDX-License-Identifier: MIT
"""Reusable runtime for generated scalar/elemental binary64 C interfaces."""
from __future__ import annotations

import ctypes as C
import os
from pathlib import Path

import numpy as np


class Distribution:
    """Own typed ctypes functions; input/output arrays remain alive during calls."""
    _specs = {}

    def __init__(self, library):
        self._library = library
        self._functions = {}
        for name, spec in self._specs.items():
            scalar = getattr(library, spec['scalar'])
            scalar.argtypes = [C.c_double] * len(spec['arguments'])
            scalar.restype = C.c_double
            bulk = getattr(library, spec['bulk'])
            bulk.argtypes = [C.c_size_t, *[C.c_void_p if a == spec['array'] else C.c_double
                                         for a in spec['arguments']],
                             C.c_void_p, C.POINTER(C.c_int)]
            bulk.restype = None
            self._functions[name] = (scalar, bulk)
            if spec.get('elementwise'):
                elementwise = getattr(library, spec['elementwise'])
                elementwise.argtypes = [C.c_size_t, *[C.c_void_p]*len(spec['arguments']),
                                        C.c_void_p, C.POINTER(C.c_int)]
                elementwise.restype = None
                self._functions[name] += (elementwise,)

    @staticmethod
    def _real_array(value, label):
        array = np.asarray(value)
        if array.dtype.kind not in 'biuf':
            raise TypeError(f'{label} must contain real numbers')
        return array

    def _evaluate(self, name, value, parameters, out):
        spec = self._specs[name]
        parameter_arrays = [self._real_array(p, label)
                            for label, p in zip(spec['parameters'], parameters)]
        if any(a.ndim for a in parameter_arrays):
            if not spec.get('elementwise'):
                raise TypeError('parameters must be scalar; broadcasting is not supported for this family')
            return self._broadcast(name, value, parameter_arrays, out)
        params = [float(a) for a in parameter_arrays]
        original = self._real_array(value, spec['array'])
        scalar, bulk = self._functions[name][:2]
        if original.ndim == 0:
            if out is not None:
                raise ValueError('out is supported only for array inputs')
            values = dict(zip(spec['parameters'], params), **{spec['array']: float(original)})
            return scalar(*(values[a] for a in spec['arguments']))
        # Require alignment as well as contiguity; arbitrary ndarray views can be unaligned.
        inputs = np.require(original, dtype=np.float64, requirements=['C', 'A'])
        if out is None:
            output = np.empty(inputs.shape, dtype=np.float64)
        else:
            if not isinstance(out, np.ndarray) or out.dtype != np.dtype(np.float64):
                raise TypeError('out must be a native-endian float64 ndarray')
            if (out.shape != inputs.shape or not out.flags.c_contiguous
                    or not out.flags.aligned or not out.flags.writeable):
                raise ValueError('out must match input shape and be contiguous, aligned and writable')
            if any(np.may_share_memory(a, out) for a in [original, *parameter_arrays]):
                raise ValueError('input and out must not overlap')
            output = out
        # Empty ndarray pointers need not be non-NULL; provide valid disjoint scratch buffers.
        xbuffer = inputs if inputs.size else np.empty(1, dtype=np.float64)
        ybuffer = output if output.size else np.empty(1, dtype=np.float64)
        values = dict(zip(spec['parameters'], params),
                      **{spec['array']: xbuffer.ctypes.data_as(C.c_void_p)})
        status = C.c_int(-1)
        bulk(inputs.size, *(values[a] for a in spec['arguments']),
             ybuffer.ctypes.data_as(C.c_void_p), C.byref(status))
        if status.value == 1:
            raise ValueError('C ABI rejected the length or uniform parameters')
        if status.value:
            raise RuntimeError(f'C ABI returned unexpected status {status.value}')
        return output

    def _broadcast(self, name, value, parameters, out):
        """Broadcast in Python; C receives only equally sized contiguous arrays."""
        spec = self._specs[name]
        originals = dict(zip(spec['parameters'], parameters),
                         **{spec['array']: self._real_array(value, spec['array'])})
        ordered = [originals[a] for a in spec['arguments']]
        arrays = np.broadcast_arrays(*ordered)
        shape = arrays[0].shape
        if out is None:
            output = np.empty(shape, dtype=np.float64)
        else:
            if not isinstance(out, np.ndarray) or out.dtype != np.dtype(np.float64):
                raise TypeError('out must be a native-endian float64 ndarray')
            if (out.shape != shape or not out.flags.c_contiguous
                    or not out.flags.aligned or not out.flags.writeable):
                raise ValueError('out must match broadcast shape and be contiguous, aligned and writable')
            if any(np.may_share_memory(a, out) for a in ordered):
                raise ValueError('input and out must not overlap')
            output = out
        buffers = [np.require(a, dtype=np.float64, requirements=['C', 'A'])
                   if a.size else np.empty(1, dtype=np.float64) for a in arrays]
        ybuffer = output if output.size else np.empty(1, dtype=np.float64)
        status = C.c_int(-1)
        self._functions[name][2](output.size,
                                *(a.ctypes.data_as(C.c_void_p) for a in buffers),
                                ybuffer.ctypes.data_as(C.c_void_p), C.byref(status))
        if status.value == 1:
            raise ValueError('C ABI rejected the length')
        if status.value:
            raise RuntimeError(f'C ABI returned unexpected status {status.value}')
        return output


class Library:
    """Load a caller-selected shared library; never search the working directory."""
    def __init__(self, path, *, dll_directories=()):
        library_path = Path(path).expanduser().resolve(strict=True)
        # add_dll_directory handles must remain live through library loading on Windows.
        handles = []
        try:
            if os.name == 'nt':
                handles = [os.add_dll_directory(str(Path(p).resolve(strict=True)))
                           for p in dll_directories]
            self._library = C.CDLL(str(library_path))
        finally:
            for handle in handles:
                handle.close()
        self.path = library_path
