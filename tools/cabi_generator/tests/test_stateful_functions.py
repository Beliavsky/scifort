# SPDX-License-Identifier: MIT
"""Scalar-result functions combine explicit state and native optional arrays."""
import copy
import json
import shutil
import unittest
import test_optional_arrays as fixtures
from test_state_codecs import SOURCE as COUNTER_SOURCE
from test_generator import run
import generate

SOURCE=COUNTER_SOURCE[:COUNTER_SOURCE.index('module example_samples')]+'''module example_samples
use iso_c_binding, only: c_double
use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
use example_codec, only: counter_state
implicit none
private
public :: evaluate
contains
function evaluate(x,shift,weight,state) result(y)
    !! @dimension n = size(x)
    !! @require n >= 0
    implicit none
    real(c_double),intent(in) :: x(:) !! shape (n); input
    real(c_double),intent(in),optional :: shift(:) !! shape (n); optional shift
    real(c_double),intent(in),optional :: weight !! optional weight
    type(counter_state),intent(inout),optional :: state !! caller counter
    real(c_double) :: y
    y=sum(x)
    if(present(shift)) y=y+100+sum(shift)
    if(present(weight)) y=y*weight
    if(present(state)) then
        state%value=state%value+1
        y=y+real(state%value,c_double)
    else
        y=-999
    end if
    if(any(x<0)) y=ieee_value(y,ieee_quiet_nan)
end function
end module
'''

class StatefulFunctionTests(fixtures.OptionalArrayTests):
    test_metadata_presence_forwarding_and_python_signature=None
    test_rejects_unsupported_native_and_flag_policies=None
    test_compiled_c_and_python_all_presence_combinations=None

    def setUp(self):
        super().setUp();self.source.write_text(SOURCE,encoding='ascii')
        self.config.update(module='example_samples',wrapper_module='example_function_api',
            argument_docs=dict(n='length',x='input',shift='shift',weight='weight',state='counter',
                               value='result',status='status',has_shift='presence'))
        self.config.pop('subroutines')
        self.config['array_functions']=dict(evaluate=dict(symbol='example_evaluate',
            dimensions=[dict(name='n',minimum=0)],shapes=dict(x=['n'],shift=['n']),
            result_argument='value',status_argument='status',invalid_status=1,optional_policy='required',
            description='Evaluate with caller state.',shape_contract=True,optional_arrays={'shift':'has_shift'},
            state_codecs={'state':dict(module='example_codec',type_module='example_codec',type='counter_state',
                load='decode',store='encode',length=2,invalid_status=3,promote_optional=True)}))
        self.policy['array_groups']['optional']['defaults']={'evaluate':{'weight':2.}}

    def test_metadata_and_promotion_validation(self):
        outputs=generate.generate(self.config,self.root)
        row=json.loads(outputs[self.root/'coverage.json'])['procedures'][0]
        self.assertEqual(row['optional_promoted_to_required'],['weight','state'])
        text=self.python_outputs()[self.root/'python/example/_generated.py']
        self.assertIn('def evaluate(self, state, x, shift=None, weight=2.0):',text)
        codec=self.config['array_functions']['evaluate']['state_codecs']['state']
        for bad in (False,1,'required',None):
            codec['promote_optional']=bad
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        codec.pop('promote_optional')
        with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        codec['promote_optional']=True
        self.source.write_text(SOURCE.replace('intent(inout),optional :: state','intent(inout) :: state'),encoding='ascii')
        with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    def test_long_formal_argument_wrapping(self):
        # Force a long interoperable formal list without changing native mathematics.
        names=['long_control_parameter_'+str(i) for i in range(5)]
        source=SOURCE.replace('x,shift,weight,state)',','.join(['x','shift','weight','state',*names])+')')
        source=source.replace('    real(c_double) :: y',''.join('    real(c_double),intent(in) :: '+n+' !! control\n' for n in names)+'    real(c_double) :: y')
        self.source.write_text(source,encoding='ascii')
        self.config['argument_docs'].update({n:'control' for n in names})
        outputs=generate.generate(self.config,self.root)
        self.assertTrue(all(len(line)<=132 for line in outputs[self.root/'adapter.f90'].splitlines()))

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc') and shutil.which('python'),'compilers required')
    def test_compiled_function_state_presence_and_error_export(self):
        outputs=self.python_outputs();outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        client=self.root/'client.c';client.write_text('''#include "adapter.h"
int main(void) {
 int64_t state[2]={2,3}; double x[2]={1,2},shift[2]={0},value=77; int status;
 example_evaluate(2,x,0,shift,2,state,&value,&status);
 if(status || value!=9 || state[0]!=3) return 1;
 state[0]=-1;value=77;
 example_evaluate(2,x,0,shift,2,state,&value,&status);
 return status!=3 || value!=77 || state[0]!=-1;
}
''',encoding='ascii')
        script=self.root/'check.py';script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
from example._runtime import StateBuffer
lib=C.CDLL(sys.argv[2]);api=Optional(lib)
state=StateBuffer(lib,[2,3]);assert api.evaluate(state,[1.,2.])==9
assert state.state[0]==3
state=StateBuffer(lib,[2,3]);assert api.evaluate(state,[1.,2.],shift=[0.,0.])==209
state=StateBuffer(lib,[2,3]);assert np.isnan(api.evaluate(state,[-1.,2.]))
assert state.state[0]==3
state=StateBuffer(lib,[-1,1])
try:api.evaluate(state,[1.,2.])
except RuntimeError:pass
else:raise AssertionError('invalid state')
assert state.state[0]==-1
''',encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            lib=folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-std=f2018','-shared','-fPIC',self.source,self.root/'adapter.f90','-o',lib],folder)
            exe=folder/'client.exe';run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',self.root,client,lib,'-o',exe],folder)
            run([exe],folder);run(['python',script,self.root/'python',lib],folder)
        if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',self.root,'-fsyntax-only',client],self.root)
