# SPDX-License-Identifier: MIT
"""Explicit Python-only index translation never narrows or mutates unchecked input."""
import shutil
import unittest
import test_python_arrays as fixtures
from test_generator import run
import generate

SOURCE='''module example_arrays
use iso_c_binding, only: c_double
implicit none
private
public :: row_sums
contains
subroutine row_sums(a, indices, y, status)
    !! @dimension nrows = size(a,1)
    !! @dimension ncols = size(a,2)
    !! @require nrows >= 1
    !! @require ncols >= 1
    implicit none
    real(c_double), intent(in) :: a(:,:) !! shape (nrows,ncols); matrix
    integer, intent(in) :: indices(:) !! shape (nrows); one-based column indices
    real(c_double), intent(out) :: y(:) !! shape (nrows); selected entries
    integer, intent(out) :: status !! zero on success
    integer :: i
    status=1
    if (any(indices<1) .or. any(indices>size(a,2))) return
    do i=1,size(a,1)
        y(i)=a(i,indices(i))
    end do
    status=0
end subroutine
end module
'''


class IndexPolicyTests(fixtures.PythonArrayGeneratorTests):
    test_generic_group_signatures_and_metadata=None
    test_missing_contract_bad_defaults_and_dimension_sources=None
    test_group_collisions_and_invalid_shapes_of_policy=None
    test_pipeline_build_and_coverage_include_array_groups=None

    def setUp(self):
        super().setUp()
        (self.root/'array_ops.f90').write_text(SOURCE,encoding='ascii')
        self.array['argument_docs']['indices']='one-based column indices'
        self.array['subroutines']['row_sums']['shapes']['indices']=['nrows']
        self.group=self.policy['array_groups']['matrices']
        self.group['defaults']={}
        self.group['zero_based_indices']={'row_sums':{'indices':'ncols'}}

    def test_metadata_and_rejected_index_policies(self):
        text=self.generate()[self.root/'python/example/_generated.py']
        self.assertIn("'zero_based_indices': {'indices': 'ncols'}",text)
        for rule in (None,[],{'unknown':{}},{'row_sums':{}},
                     {'row_sums':{'a':'ncols'}},{'row_sums':{'indices':'unknown'}},
                     {'row_sums':{'indices':True}},{'row_sums':{'y':'ncols'}}):
            self.group['zero_based_indices']=rule
            with self.assertRaises(ValueError):self.generate()
        self.group['zero_based_indices']={'row_sums':{'indices':'ncols'}}
        self.array['subroutines']['row_sums']['integer_limits']={'indices':dict(minimum=1,maximum=10)}
        with self.assertRaises(ValueError):self.generate()

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('python'),'compiler/Python required')
    def test_compiled_conversion_bounds_layout_and_aliasing(self):
        outputs=self.generate();outputs.update(generate.generate(self.array,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        script=self.root/'check.py'
        script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Matrices
api=Matrices(C.CDLL(sys.argv[2]))
a=np.array([[10.,20.,30.],[40.,50.,60.]],order='F')
indices=np.array([2,0],dtype=np.int64);indices.flags.writeable=False
for item in (indices,indices.astype('>i8'),np.array([2,9,0])[::2],indices.astype(np.uint64)):
    out=np.empty(2)
    result=api.row_sums(a,item,out={'y':out})
    assert result.status==0 and result.outputs['y'] is out
    np.testing.assert_array_equal(out,[30.,40.])
    np.testing.assert_array_equal(item,[2,0])
for bad in ([-1,0],[3,0],np.array([2**64-1,0],dtype=np.uint64),[2**40,0],
            [0.,1.],[True,False],[[0,1]],np.array([],dtype=int)):
    out=np.array([99.,98.])
    try:api.row_sums(a,bad,out={'y':out})
    except (ValueError,TypeError):pass
    else:raise AssertionError(bad)
    np.testing.assert_array_equal(out,[99.,98.])
spec=api._specs['row_sums']
assert spec['zero_based_indices']=={'indices':'ncols'}
''',encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            library=folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-shared','-fPIC',self.root/'array_ops.f90',
                 self.root/'array_api.f90','-o',library],folder)
            run(['python',script,self.root/'python',library],folder)
