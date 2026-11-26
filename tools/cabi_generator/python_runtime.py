# SPDX-License-Identifier: MIT
"""Reusable runtime for generated scalar/elemental binary64 C interfaces."""
from __future__ import annotations

import ctypes as C
from collections.abc import Mapping
from typing import NamedTuple
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


class ArrayResult(NamedTuple):
    """Native or synthetic call status and named outputs; numerical NaNs are retained."""
    status: int
    outputs: dict


class ArrayGroup:
    """Typed, descriptor-free C array calls driven by reviewed shape metadata."""
    _specs = {}

    @staticmethod
    def _integer_array(value, name):
        array = np.asarray(value)
        if array.dtype.kind not in 'iu':
            raise TypeError(f'{name} must contain integers, without floats or booleans')
        return array

    def __init__(self, library):
        self._library = library
        self._functions = {}
        for name, spec in self._specs.items():
            function = getattr(library, spec['symbol'])
            types = [C.c_int]*len(spec['dimensions'])
            for arg in spec['arguments']:
                if arg == spec['status']:
                    types.append(C.POINTER(C.c_int))
                    continue
                info = spec['metadata'][arg]
                if info['rank']:
                    types.append(C.c_void_p)
                elif info['category'] == 'logical':
                    types.append(C.c_int)
                else:
                    scalar = C.c_int if info['category'] == 'integer' else C.c_double
                    types.append(scalar if info['intent'] == 'in' else C.POINTER(scalar))
            function.argtypes, function.restype = types, None
            self._functions[name] = function

    def _invoke(self, name, inputs, out):
        spec = self._specs[name]
        originals, values = {}, {}
        for arg in spec['inputs']:
            info, value = spec['metadata'][arg], inputs[arg]
            if value is None and arg in spec.get('array_defaults', {}):
                continue
            if info['category'] == 'logical':
                if not isinstance(value, (bool, np.bool_)):
                    raise TypeError(f'{arg} must be a boolean')
                values[arg] = spec['logical_inputs'][arg]['true' if value else 'false']
                continue
            array = (self._integer_array(value, arg) if info['category'] == 'integer'
                     else Distribution._real_array(value, arg))
            if array.ndim != info['rank']:
                raise ValueError(f'{arg} must have rank {info["rank"]}')
            originals[arg] = array
            if not info['rank']:
                values[arg] = int(array) if info['category'] == 'integer' else float(array)
        dimensions = {}
        for dim in spec['dimensions']:
            definition = spec['definitions'][dim]
            array = originals.get(definition['argument'])
            axis = definition['axis']
            length = (spec['dimension_defaults'][dim] if array is None else
                      array.shape[axis-1] if axis else array.size)
            if not spec['minimums'][dim] <= length <= np.iinfo(np.intc).max:
                raise ValueError(f'{dim} is outside the supported C integer extent range')
            dimensions[dim] = length
        for dim, rule in spec.get('length_choices', {}).items():
            if dimensions[dim] not in [dimensions[rule['dimension']] + offset for offset in rule['offsets']]:
                raise ValueError(f'{dim} does not satisfy the required dimension choices')
        for arg, array in originals.items():
            if spec['metadata'][arg]['category'] != 'integer' or not array.size:
                continue
            minimum, maximum = int(np.min(array)), int(np.max(array))
            limits = spec.get('integer_limits', {}).get(arg, dict(minimum=-2147483648, maximum=2147483647))
            if minimum < limits['minimum'] or maximum > limits['maximum']:
                raise ValueError(f'{arg} is outside the supported C integer range or configured limits')
        for arg, ceiling in spec.get('integer_sums', {}).items():
            # Rank one and C-int extents/elements bound the widened sum below 2**63.
            total = np.sum(np.abs(originals[arg].astype(np.int64)), dtype=np.int64)
            if total > ceiling:
                raise ValueError(f'{arg} exceeds the configured absolute integer sum limit')
        shapes = {a: tuple(dimensions[d] for d in dims) for a, dims in spec['shapes'].items()}
        for arg, fill in spec.get('array_defaults', {}).items():
            if inputs[arg] is None:
                originals[arg] = (np.eye(shapes[arg][0], dtype=np.float64, order='F') if fill == 'identity'
                                  else np.zeros(shapes[arg], dtype=np.float64, order='F'))
        input_arrays = [a for arg, a in originals.items() if spec['metadata'][arg]['rank']]
        for index, array in enumerate(input_arrays):
            if any(np.may_share_memory(array, other) for other in input_arrays[:index]):
                raise ValueError('input buffers must not overlap')
        for arg, original in originals.items():
            if arg in shapes and original.shape != shapes[arg]:
                raise ValueError(f'{arg} must have shape {shapes[arg]}')
        result = spec.get('result_argument')
        if result and spec['metadata'][result]['rank'] and out is not None:
            if not isinstance(out, np.ndarray):
                raise TypeError('array-function out must be a NumPy ndarray')
            out = {result: out}
        if out is None:
            out = {}
        if not isinstance(out, Mapping):
            raise TypeError('out must map array output names to NumPy buffers')
        array_outputs = {a for a in spec['outputs'] if spec['metadata'][a]['rank']}
        if set(out) - array_outputs:
            raise ValueError('out contains an unknown or scalar output name')
        outputs, scalars, buffers = {}, {}, {}
        for arg in spec['outputs']:
            info = spec['metadata'][arg]
            if not info['rank']:
                cls = C.c_int if info['category'] == 'integer' else C.c_double
                scalars[arg] = cls(0 if cls is C.c_int else float('nan'))
                values[arg] = C.byref(scalars[arg])
                continue
            if arg in out:
                array = out[arg]
                if not isinstance(array, np.ndarray) or array.dtype != np.dtype(np.float64):
                    raise TypeError(f'out[{arg!r}] must be a native-endian float64 ndarray')
                if (array.shape != shapes[arg] or not array.flags.f_contiguous
                        or not array.flags.aligned or not array.flags.writeable):
                    raise ValueError(f'out[{arg!r}] must match shape and be Fortran-contiguous, aligned and writable')
                if any(np.may_share_memory(array, other) for other in originals.values()):
                    raise ValueError('output must not overlap any input')
                if any(np.may_share_memory(array, other) for other in outputs.values()):
                    raise ValueError('output buffers must not overlap each other')
            else:
                array = np.full(shapes[arg], np.nan, dtype=np.float64, order='F')
            outputs[arg] = array
        # Prepare C memory only after all Python validation. Keep every buffer alive.
        for arg, original in originals.items():
            if spec['metadata'][arg]['rank']:
                dtype = np.intc if spec['metadata'][arg]['category'] == 'integer' else np.float64
                buffers[arg] = (np.require(original, dtype=dtype, requirements=['F', 'A'])
                                if original.size else np.empty(1, dtype=dtype))
                values[arg] = buffers[arg].ctypes.data_as(C.c_void_p)
        for arg, array in outputs.items():
            buffers[arg] = array if array.size else np.empty(1, dtype=np.float64)
            values[arg] = buffers[arg].ctypes.data_as(C.c_void_p)
        status = C.c_int(-1)
        values[spec['status']] = C.byref(status)
        self._functions[name](*(dimensions[d] for d in spec['dimensions']),
                              *(values[a] for a in spec['arguments']))
        for arg, scalar in scalars.items():
            outputs[arg] = scalar.value
        if 'result_argument' in spec:
            if status.value:
                raise RuntimeError(f'C ABI returned unexpected status {status.value}')
            return outputs[spec['result_argument']]
        return ArrayResult(status.value, outputs)


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
