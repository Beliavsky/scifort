# SPDX-License-Identifier: MIT
"""Portable interleaved binary64 complex arrays and success-only output publication."""
import json
import shutil
import unittest
import test_optional_arrays as fixtures
from test_generator import run
import generate

SOURCE = '''module example_optional
use iso_c_binding, only: c_double
implicit none
private
public :: stack
contains
subroutine stack(x,fail,cube,status)
    !! @dimension m = size(x,1)
    !! @dimension n = size(x,2)
    !! @dimension b = size(cube,3)
    !! @require m >= 0
    !! @require n >= 0
    !! @require b >= 0
    complex(c_double),intent(in) :: x(:,:) !! shape (m,n); base matrix
    integer,intent(in) :: fail !! failure selector
    complex(c_double),intent(out) :: cube(:,:,:) !! shape (m,n,b); repeated matrices
    integer,intent(out) :: status !! zero on success
    integer :: k
    status=2
    if(fail<0) return
    do k=1,size(cube,3)
        cube(:,:,k)=conjg(x)+cmplx(100*k,-k,kind=c_double)
        if(fail>0) return
    end do
    status=0
end subroutine
end module
'''


class ComplexArrayTests(fixtures.OptionalArrayTests):
    test_metadata_presence_forwarding_and_python_signature = None
    test_rejects_unsupported_native_and_flag_policies = None
    test_compiled_c_and_python_all_presence_combinations = None

    def setUp(self):
        super().setUp()
        self.source.write_text(SOURCE, encoding='ascii')
        self.config['outputs']['shape_documentation'] = 'shapes.md'
        self.config['argument_docs'] = dict(m='rows', n='columns', b='batch count',
            x='base matrix pairs', fail='failure selector', cube='result pairs', status='status')
        self.config['subroutines'] = dict(stack=dict(symbol='example_stack',
            dimensions=[dict(name=d, minimum=0) for d in ('m', 'n', 'b')],
            shapes=dict(x=['m','n'], cube=['m','n','b']), status_argument='status', invalid_status=1,
            description='Build a complex matrix stack.', shape_contract=True))
        self.policy['array_groups']['optional']['output_dimensions'] = dict(stack=dict(b=1))

    def test_metadata_transport_and_commit_guard(self):
        outputs = generate.generate(self.config, self.root)
        self.assertEqual(outputs, generate.generate(self.config, self.root))
        source = outputs[self.root/'adapter.f90']
        self.assertIn('x(2,m,n)', source)
        self.assertIn('cube(2,m,n,b)', source)
        self.assertIn('complex(cg_kind), allocatable :: cg_complex_cube(:,:,:)', source)
        self.assertIn('cg_complex_x=cmplx(x(1,:,:), x(2,:,:), kind=cg_kind)', source)
        self.assertIn('if (status == 0_c_int) then', source)
        self.assertLess(source.index('if (m <'), source.index('allocate(cg_complex_'))
        self.assertIn('const double * arg_x', outputs[self.root/'adapter.h'])
        self.assertIn('`(2,m,n,b)`', outputs[self.root/'shapes.md'])
        self.assertIn('errors preserve it', outputs[self.root/'shapes.md'])
        report = json.loads(outputs[self.root/'coverage.json'])['procedures'][0]
        self.assertEqual(report['complex_arrays'], dict(x='interleaved_f64', cube='interleaved_f64'))
        self.assertIn("'category': 'complex'", self.python_outputs()[self.root/'python/example/_generated.py'])

    def test_rejects_complex_scalars_inout_and_unsupported_ranks(self):
        for source in (SOURCE.replace('x(:,:)', 'x'),
                       SOURCE.replace('intent(out) :: cube', 'intent(inout) :: cube'),
                       SOURCE.replace('x(:,:)', 'x(:,:,:)'),
                       SOURCE.replace('cube(:,:,:)', 'cube(:,:,:,:)'),
                       SOURCE.replace('complex(c_double)', 'complex(kind=4)')):
            self.source.write_text(source, encoding='ascii')
            with self.assertRaises(ValueError): generate.generate(self.config, self.root)

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc'), 'GNU compilers required')
    def test_checked_optimized_c_cpp_and_python(self):
        outputs = self.python_outputs()
        outputs.update(generate.generate(self.config, self.root))
        for path, text in outputs.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding='utf-8')
        client = self.root/'client.c'
        client.write_text('''#include "adapter.h"
int main(void) {
 double x[12]={1,11,2,12,3,13,4,14,5,15,6,16},out[24]; int status;
 example_stack(2,3,2,x,0,out,&status);
 if(status || out[0]!=101 || out[1]!=-12 || out[10]!=106 || out[11]!=-17) return 1;
 if(out[12]!=201 || out[13]!=-13 || out[22]!=206 || out[23]!=-18) return 2;
 for(int i=0;i<24;++i) out[i]=77;
 example_stack(2,3,2,x,1,out,&status); if(status!=2) return 3;
 for(int i=0;i<24;++i) if(out[i]!=77) return 4;
 example_stack(2,3,2,x,-1,out,&status); if(status!=2 || out[0]!=77) return 5;
 example_stack(2,3,-1,x,0,out,&status);
 return status!=1 || out[0]!=77;
}
''', encoding='ascii')
        script = self.root/'check.py'
        script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
api=Optional(C.CDLL(sys.argv[2]))
x=(np.arange(6.).reshape(2,3)+1j*np.arange(11.,17.).reshape(2,3))[:,::-1]
out=np.empty((2,3,4),dtype=np.complex128,order='F')
r=api.stack(x,0,b=4,out={'cube':out})
assert r.status==0 and r.outputs['cube'] is out
for k in range(4): np.testing.assert_array_equal(out[:,:,k],x.conj()+100*(k+1)-1j*(k+1))
for flag in (-1,1):
 out[:]=77+88j
 r=api.stack(x,flag,b=4,out={'cube':out})
 assert r.status==2
 np.testing.assert_array_equal(out,77+88j)
 r=api.stack(x,flag,b=4)
 assert r.status==2 and np.all(np.isnan(r.outputs['cube'].real)) and np.all(np.isnan(r.outputs['cube'].imag))
for bad in (np.empty((2,3,4),order='F'),np.empty((2,3,4),dtype=np.complex64,order='F'),
            np.empty((2,3,4),dtype=np.complex128),np.empty((2,3,4),dtype='>c16',order='F')):
 try: api.stack(x,0,b=4,out={'cube':bad})
 except (TypeError,ValueError): pass
 else: raise AssertionError('invalid complex buffer')
assert api.stack(x,0,b=0).outputs['cube'].shape==(2,3,0)
r=api.stack(x.real,0,b=1)
np.testing.assert_array_equal(r.outputs['cube'][:,:,0],x.real+100-1j)
''', encoding='ascii')
        for index, flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'], ['-O3'])):
            folder = self.root/str(index); folder.mkdir()
            lib = folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-std=f2018','-shared','-fPIC',self.source,self.root/'adapter.f90','-o',lib],folder)
            for compiler, standard in (('gcc','c11'), ('g++','c++17')):
                if not shutil.which(compiler): continue
                exe = folder/(compiler+'.exe')
                run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',self.root,client,lib,'-o',exe],folder)
                run([exe],folder)
            run(['python',script,self.root/'python',lib],folder)
