# SPDX-License-Identifier: MIT
"""Rank-one/two default-integer outputs use C int and native NumPy intc."""
import shutil
import unittest
import test_optional_arrays as fixtures
from test_generator import run
import generate

SOURCE='''module example_optional
use iso_c_binding, only: c_double
implicit none
private
public :: counts
contains
subroutine counts(x,vector,matrix,status)
    !! @dimension n = size(x)
    !! @require n >= 0
    implicit none
    real(c_double),intent(in) :: x(:) !! shape (n); input
    integer,intent(out) :: vector(:) !! shape (n); integer vector
    integer,intent(out) :: matrix(:,:) !! shape (n,n); integer matrix
    integer,intent(out) :: status !! zero on success
    integer :: i
    vector=int(x)
    do i=1,size(x)
        matrix(:,i)=vector+i
    end do
    status=0
end subroutine
end module
'''

class IntegerOutputTests(fixtures.OptionalArrayTests):
    test_metadata_presence_forwarding_and_python_signature=None
    test_rejects_unsupported_native_and_flag_policies=None
    test_compiled_c_and_python_all_presence_combinations=None

    def setUp(self):
        super().setUp()
        self.source.write_text(SOURCE,encoding='ascii')
        self.config['argument_docs']=dict(n='length',x='values',vector='counts',matrix='columns',status='status')
        self.config['subroutines']=dict(counts=dict(symbol='example_counts',
            dimensions=[dict(name='n',minimum=0)],shapes=dict(x=['n'],vector=['n'],matrix=['n','n']),
            status_argument='status',invalid_status=2,description='Integer outputs.',shape_contract=True))

    def test_integer_declarations_and_inout_rejection(self):
        outputs=generate.generate(self.config,self.root)
        self.assertIn('int * arg_vector',outputs[self.root/'adapter.h'])
        self.assertIn('integer(c_int), intent(inout) :: matrix(n,n)',outputs[self.root/'adapter.f90'])
        self.source.write_text(SOURCE.replace('intent(out) :: vector','intent(inout) :: vector'),encoding='ascii')
        with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc') and shutil.which('python'),'compilers required')
    def test_compiled_integer_outputs(self):
        outputs=self.python_outputs();outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        client=self.root/'client.c'
        client.write_text('''#include "adapter.h"
int main(void) {
 double x[2]={3,7}; int v[2]={0},m[4]={0},s;
 example_counts(2,x,v,m,&s);
 if(s || v[0]!=3 || v[1]!=7 || m[0]!=4 || m[3]!=9) return 1;
 example_counts(-1,x,v,m,&s);
 return s!=2 || v[0]!=3;
}
''',encoding='ascii')
        script=self.root/'check.py'
        script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
api=Optional(C.CDLL(sys.argv[2]))
r=api.counts([3.,7.]);assert r.status==0
assert r.outputs['vector'].dtype==np.dtype(np.intc)
assert r.outputs['matrix'].flags.f_contiguous
np.testing.assert_array_equal(r.outputs['matrix'],[[4,5],[8,9]])
v=np.empty(2,dtype=np.intc);m=np.empty((2,2),dtype=np.intc,order='F')
r=api.counts([3.,7.],out={'vector':v,'matrix':m})
assert r.outputs['vector'] is v and r.outputs['matrix'] is m
for bad in (np.empty(2),np.empty(2,dtype=np.int64),np.empty(4,dtype=np.intc)[::2],np.empty(2,dtype='>i4')):
 try:api.counts([3.,7.],out={'vector':bad})
 except (ValueError,TypeError):pass
 else:raise AssertionError(bad.dtype)
v.flags.writeable=False
try:api.counts([3.,7.],out={'vector':v})
except ValueError:pass
else:raise AssertionError('readonly')
r=api.counts([]);assert r.status==0 and r.outputs['matrix'].shape==(0,0)
''',encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            lib=folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-std=f2018','-shared','-fPIC',self.source,self.root/'adapter.f90','-o',lib],folder)
            exe=folder/'client.exe'
            run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',self.root,client,lib,'-o',exe],folder)
            run([exe],folder);run(['python',script,self.root/'python',lib],folder)
        if shutil.which('g++'):
            run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',self.root,'-fsyntax-only',client],self.root)
