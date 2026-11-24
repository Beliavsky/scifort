# SPDX-License-Identifier: MIT
"""Independent integer-output and explicitly mapped logical-input adapters."""
import ctypes as C
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from test_generator import array_fixture_config, run
import generate

SOURCE = """module example_arrays
use iso_c_binding, only: c_double
implicit none
private
public :: summarize
contains
pure subroutine summarize(a,total,count,status,include_negative)
implicit none
real(c_double), intent(in) :: a(:,:) !! shape (nrows,ncols); input matrix
real(c_double), intent(out) :: total !! selected sum
integer, intent(out) :: count !! selected count
integer, intent(out) :: status !! zero on success
logical, intent(in), optional :: include_negative !! count negative entries, default true
logical :: include
include=.true.
if (present(include_negative)) include=include_negative
if (include) then
 total=sum(a)
 count=size(a)
else
 total=sum(a,mask=a>=0)
 count=count_nonnegative(a)
end if
status=0
end subroutine
pure function count_nonnegative(a) result(n)
implicit none
real(c_double), intent(in) :: a(:,:) !! input matrix
integer :: n
n=count(a>=0)
end function
end module
"""


class AdapterTypeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'array_ops.f90').write_text(SOURCE)
        self.config = array_fixture_config()
        entry = self.config['subroutines'].pop('row_sums')
        entry.update(symbol='example_summary',shapes=dict(a=['nrows','ncols']),
                     logical_inputs=dict(include_negative=dict(false=0,true=1)))
        self.config['subroutines']['summarize'] = entry
        self.config['argument_docs'].update(total='selected sum',count='selected count',
                                            include_negative='0 omit negatives; 1 include negatives')

    def test_explicit_mapping_validation_and_coverage(self):
        outputs = generate.generate(self.config,self.root)
        record = json.loads(outputs[self.root/'array_coverage.json'])['procedures'][0]
        self.assertEqual(record['optional_promoted_to_required'],['include_negative'])
        self.assertEqual(record['native_shape_documentation_review'],[])
        header = outputs[self.root/'array_api.h']
        self.assertIn('int * arg_count',header)
        self.assertIn('int arg_include_negative',header)
        policy = self.config['subroutines']['summarize']
        for mapping in ({},dict(status=dict(false=0,true=1)),dict(total=dict(false=0,true=1)),
                        dict(include_negative=dict(false=False,true=1)),
                        dict(include_negative=dict(false=1,true=1)),
                        dict(include_negative=dict(false=0,true=2147483648))):
            policy['logical_inputs'] = mapping
            with self.assertRaises(ValueError):
                generate.generate(self.config,self.root)

    def test_checked_and_optimized_conversion(self):
        fc = shutil.which('gfortran')
        if not fc:
            self.skipTest('gfortran required')
        for path,text in generate.generate(self.config,self.root).items():
            path.write_text(text)
        dll_dir = os.add_dll_directory(str(Path(fc).parent)) if os.name == 'nt' else None
        try:
            for profile,flags in (('checked',['-O0','-fcheck=all']),('optimized',['-O3'])):
                folder = self.root/profile
                folder.mkdir()
                path = folder/('example.dll' if os.name=='nt' else 'example.so')
                run([fc,*flags,'-std=f2018','-fimplicit-none','-Werror=implicit-interface',
                     '-shared','-fPIC','-J',folder,self.root/'array_ops.f90',self.root/'array_api.f90',
                     '-o',path],folder)
                library = C.CDLL(str(path))
                try:
                    func = library.example_summary
                    func.argtypes = [C.c_int,C.c_int,C.c_void_p,C.POINTER(C.c_double),
                                     C.POINTER(C.c_int),C.POINTER(C.c_int),C.c_int]
                    func.restype = None
                    a = (C.c_double*4)(1,-2,3,-4)
                    total,count,status = C.c_double(99),C.c_int(99),C.c_int(-1)
                    for flag,expected,n in ((0,4,2),(1,-2,4)):
                        func(2,2,a,C.byref(total),C.byref(count),C.byref(status),flag)
                        self.assertEqual((total.value,count.value,status.value),(expected,n,0))
                    for rows,flag in ((0,1),(2,2),(2,-1)):
                        total.value=count.value=99
                        func(rows,2,a,C.byref(total),C.byref(count),C.byref(status),flag)
                        self.assertEqual((total.value,count.value,status.value),(99,99,1))
                finally:
                    import _ctypes
                    (_ctypes.FreeLibrary if os.name=='nt' else _ctypes.dlclose)(library._handle)
        finally:
            if dll_dir:
                dll_dir.close()
