# SPDX-License-Identifier: MIT
"""Rank-three real/integer output and third-axis dimension contracts."""
import shutil
import unittest
import test_optional_arrays as fixtures
from test_generator import run
import generate

SOURCE='''module example_optional
use iso_c_binding, only: c_double
implicit none
private
public :: stack
contains
subroutine stack(x,cube,status)
    !! @dimension m = size(x,1)
    !! @dimension n = size(x,2)
    !! @dimension b = size(cube,3)
    !! @require m >= 0
    !! @require n >= 0
    !! @require b >= 0
    implicit none
    real(c_double),intent(in) :: x(:,:) !! shape (m,n); base matrix
    real(c_double),intent(out) :: cube(:,:,:) !! shape (m,n,b); repeated matrices
    integer,intent(out) :: status !! zero on success
    integer :: k
    do k=1,size(cube,3)
        cube(:,:,k)=x+100*k
    end do
    status=0
end subroutine
end module
'''

class RankThreeOutputTests(fixtures.OptionalArrayTests):
    test_metadata_presence_forwarding_and_python_signature=None
    test_rejects_unsupported_native_and_flag_policies=None
    test_compiled_c_and_python_all_presence_combinations=None

    def setUp(self):
        super().setUp();self.source.write_text(SOURCE,encoding='ascii')
        self.config['argument_docs']=dict(m='rows',n='columns',b='batch size',x='matrix',cube='stack',status='status')
        self.config['subroutines']=dict(stack=dict(symbol='example_stack',
            dimensions=[dict(name=d,minimum=0) for d in ('m','n','b')],
            shapes=dict(x=['m','n'],cube=['m','n','b']),status_argument='status',invalid_status=1,
            description='Build a matrix stack.',shape_contract=True))
        self.policy['array_groups']['optional']['output_dimensions']={'stack':{'b':1}}

    def test_metadata_and_rank_restrictions(self):
        outputs=generate.generate(self.config,self.root)
        self.assertIn('cube(m,n,b)',outputs[self.root/'adapter.f90'])
        self.assertIn('def stack(self, x, *, b=1, out=None):',self.python_outputs()[self.root/'python/example/_generated.py'])
        for source in (SOURCE.replace('intent(out) :: cube','intent(in) :: cube'),
                       SOURCE.replace('intent(out) :: cube','intent(inout) :: cube'),
                       SOURCE.replace('real(c_double),intent(out) :: cube','integer,intent(in) :: cube'),
                       SOURCE.replace('cube(:,:,:)', 'cube(:,:,:,:)'),
                       SOURCE.replace('size(cube,3)', 'size(cube,4)')):
            self.source.write_text(source,encoding='ascii')
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc') and shutil.which('python'),'compilers required')
    def test_compiled_layout_empty_and_buffers(self):
        outputs=self.python_outputs();outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        client=self.root/'client.c'
        client.write_text('''#include "adapter.h"
int main(void) {
 double x[6]={1,2,3,4,5,6},out[12]; int status;
 example_stack(2,3,2,x,out,&status);
 if(status || out[0]!=101 || out[5]!=106 || out[6]!=201 || out[11]!=206) return 1;
 out[0]=77; example_stack(2,3,-1,x,out,&status);
 return status!=1 || out[0]!=77;
}
''',encoding='ascii')
        script=self.root/'check.py'
        script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
api=Optional(C.CDLL(sys.argv[2]));x=np.arange(6.).reshape(2,3)
out=np.empty((2,3,4),order='F');r=api.stack(x,b=4,out={'cube':out})
assert r.status==0 and r.outputs['cube'] is out and out.flags.f_contiguous
for k in range(4):np.testing.assert_array_equal(out[:,:,k],x+100*(k+1))
for bad in (np.empty((2,3,4)),np.empty((2,3,4),dtype=np.intc,order='F'),np.empty((3,2,4),order='F')):
 try:api.stack(x,b=4,out={'cube':bad})
 except (ValueError,TypeError):pass
 else:raise AssertionError('invalid output')
assert api.stack(x,b=0).outputs['cube'].shape==(2,3,0)
assert api.stack(np.empty((0,3)),b=2).outputs['cube'].shape==(0,3,2)
for bad in (-1,True,2**31,1.5):
 try:api.stack(x,b=bad)
 except (ValueError,TypeError):pass
 else:raise AssertionError('invalid batch size')
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

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc'), 'GNU compilers required')
    def test_integer_cube_layout_and_python_buffers(self):
        self.source.write_text(SOURCE.replace('real(c_double),intent', 'integer,intent'), encoding='ascii')
        outputs = self.python_outputs()
        outputs.update(generate.generate(self.config, self.root))
        for path, text in outputs.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding='utf-8')
        client = self.root/'client.c'
        client.write_text('''#include "adapter.h"
int main(void) {
 int x[6]={1,2,3,4,5,6},out[12],status;
 example_stack(2,3,2,x,out,&status);
 if(status || out[0]!=101 || out[5]!=106 || out[6]!=201 || out[11]!=206) return 1;
 out[0]=77; example_stack(2,3,-1,x,out,&status);
 return status!=1 || out[0]!=77;
}
''', encoding='ascii')
        script = self.root/'check.py'
        script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
api=Optional(C.CDLL(sys.argv[2]));x=np.arange(6,dtype=np.intc).reshape(2,3)
out=np.empty((2,3,4),dtype=np.intc,order='F');r=api.stack(x,b=4,out={'cube':out})
assert r.status==0 and r.outputs['cube'] is out
for k in range(4): np.testing.assert_array_equal(out[:,:,k],x+100*(k+1))
assert api.stack(x,b=0).outputs['cube'].shape==(2,3,0)
for bad in (np.empty((2,3,4),order='F'),np.empty((2,3,4),dtype=np.intc)):
 try: api.stack(x,b=4,out={'cube':bad})
 except (ValueError,TypeError): pass
 else: raise AssertionError('invalid integer cube buffer')
''', encoding='ascii')
        for index, flags in enumerate((['-O0','-g','-fcheck=all','-ftrapv'], ['-O3'])):
            folder = self.root/str(index); folder.mkdir()
            lib = folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-std=f2018','-shared','-fPIC',self.source,self.root/'adapter.f90','-o',lib],folder)
            exe = folder/'client.exe'
            run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',self.root,client,lib,'-o',exe],folder)
            run([exe],folder)
            run(['python',script,self.root/'python',lib],folder)
