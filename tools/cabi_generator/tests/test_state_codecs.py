# SPDX-License-Identifier: MIT
"""Independent two-int64 state codec and explicit output-dimension fixtures."""
import copy
import json
import shutil
import unittest
import test_optional_arrays as fixtures
from test_generator import run
import generate

SOURCE='''module example_codec
use iso_c_binding, only: c_int64_t
implicit none
type :: counter_state
    integer(c_int64_t) :: value=0
end type
contains
subroutine decode(values,state,status)
    implicit none
    integer(c_int64_t),intent(in) :: values(2)
    type(counter_state),intent(out) :: state
    integer,intent(out) :: status
    status=1
    if(any(values<0)) return
    state%value=values(1)
    status=0
end subroutine
subroutine encode(state,values)
    implicit none
    type(counter_state),intent(in) :: state
    integer(c_int64_t),intent(out) :: values(2)
    values=[state%value,state%value+1]
end subroutine
end module
module example_samples
use iso_c_binding,only:c_double
use example_codec,only:counter_state
implicit none
private
public :: draw
contains
subroutine draw(state,samples,scale,status)
    !! @dimension k = size(scale)
    !! @dimension n = size(samples,2)
    !! @require k >= 0
    !! @require n >= 0
    implicit none
    type(counter_state),intent(inout) :: state !! explicit counter
    real(c_double),intent(out) :: samples(:,:) !! shape (k,n); output columns
    real(c_double),intent(in) :: scale(:) !! shape (k); offsets
    integer,intent(out),optional :: status !! zero on success
    integer :: i,j
    if(present(status)) status=0
    if(size(samples,2)==0) return
    if(any(scale<0)) then
        if(present(status)) status=1
        return
    end if
    do j=1,size(samples,2)
        do i=1,size(scale)
            state%value=state%value+1
            samples(i,j)=real(state%value,c_double)+scale(i)
        end do
    end do
    if(any(scale==99)) then
        if(present(status)) status=2
    end if
end subroutine
end module
'''


class StateCodecTests(fixtures.OptionalArrayTests):
    test_metadata_presence_forwarding_and_python_signature=None
    test_rejects_unsupported_native_and_flag_policies=None
    test_compiled_c_and_python_all_presence_combinations=None

    def setUp(self):
        super().setUp()
        self.source.write_text(SOURCE,encoding='ascii')
        self.config.update(module='example_samples',wrapper_module='example_sample_api',
            argument_docs=dict(k='rows',n='columns',state='two components',samples='columns',scale='offsets',status='status'),
            subroutines=dict(draw=dict(symbol='example_draw',
                dimensions=[dict(name='k',minimum=0),dict(name='n',minimum=0)],
                shapes=dict(samples=['k','n'],scale=['k']),status_argument='status',invalid_status=4,
                optional_status=True,description='Draw stateful samples.',shape_contract=True,
                state_codecs=dict(state=dict(module='example_codec',type_module='example_codec',
                    type='counter_state',load='decode',store='encode',length=2,invalid_status=3)))))
        group=self.policy['array_groups']['optional']
        group['output_dimensions']={'draw':{'n':1}}

    def test_metadata_and_invalid_codec_policies(self):
        outputs=generate.generate(self.config,self.root)
        row=json.loads(outputs[self.root/'coverage.json'])['procedures'][0]
        self.assertEqual(row['state_codecs']['state']['length'],2)
        self.assertTrue(row['optional_status_promoted_to_required'])
        self.assertIn('int64_t * arg_state',outputs[self.root/'adapter.h'])
        text=self.python_outputs()[self.root/'python/example/_generated.py']
        self.assertIn('def draw(self, state, scale, *, n=1, out=None):',text)
        baseline=copy.deepcopy(self.config)
        for update in ({'length':0},{'length':True},{'invalid_status':0},{'type':'other'},
                       {'module':'bad-name'},{'extra':1}):
            self.config=copy.deepcopy(baseline)
            self.config['subroutines']['draw']['state_codecs']['state'].update(update)
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        self.config=copy.deepcopy(baseline)
        for value in (None,{},[],{'status':baseline['subroutines']['draw']['state_codecs']['state']}):
            self.config['subroutines']['draw']['state_codecs']=value
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        self.config=copy.deepcopy(baseline)
        self.config['subroutines']['draw']['optional_status']=False
        with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    def test_invalid_output_dimensions_and_state_defaults(self):
        group=self.policy['array_groups']['optional']
        for values in (None,[],{}, {'unknown':{}},{'draw':{}},{'draw':{'k':1}},
                       {'draw':{'unknown':1}},{'draw':{'n':-1}},{'draw':{'n':True}},
                       {'draw':{'n':2147483648}}):
            group['output_dimensions']=values
            with self.assertRaises(ValueError):self.python_outputs()
        group['output_dimensions']={'draw':{'n':1}}
        group['defaults']={'draw':{'state':1}}
        with self.assertRaises(ValueError):self.python_outputs()

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc') and shutil.which('python'),
                         'compilers and Python required')
    def test_compiled_transport_errors_and_empty_output(self):
        outputs=self.python_outputs();outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        client=self.root/'client.c'
        client.write_text('''#include "adapter.h"
int main(void) {
 int64_t state[2]={1099511627776LL,1}; double out[2]={77,88},scale[1]={1}; int status;
 example_draw(1,2,state,out,scale,&status);
 if(status || state[0]!=1099511627778LL || out[0]!=1099511627778.) return 1;
 state[0]=-1;out[0]=77;
 example_draw(1,2,state,out,scale,&status);
 return status!=3 || out[0]!=77 || state[0]!=-1;
}
''',encoding='ascii')
        script=self.root/'check.py'
        script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
from example._runtime import StateBuffer
lib=C.CDLL(sys.argv[2]);api=Optional(lib)
state=StateBuffer(lib,[2**40,1]);before=state.state
r=api.draw(state,[1.,2.],n=2)
assert r.status==0 and r.outputs['samples'].shape==(2,2)
assert state.state[0]==2**40+4
assert r.outputs['samples'][0,0]==2**40+2
saved=state.state
assert api.draw(state,[1.],n=0).outputs['samples'].shape==(1,0)
np.testing.assert_array_equal(state.state,saved)
out=np.full((1,1),77.,order='F')
invalid=StateBuffer(lib,[-1,1])
assert api.draw(invalid,[1.],out={'samples':out}).status==3
assert out[0,0]==77 and invalid.state[0]==-1
assert api.draw(state,[-1.]).status==1
np.testing.assert_array_equal(state.state,saved)
assert api.draw(state,[99.]).status==2
assert state.state[0]==saved[0]+1
for bad in (-1,2**31,True,1.5):
 try:api.draw(state,[1.],n=bad)
 except (ValueError,TypeError):pass
 else:raise AssertionError(bad)
for bad in ([1,2],StateBuffer(lib,[1]),StateBuffer(C.CDLL(sys.argv[2]),[1,2])):
 try:api.draw(bad,[1.])
 except (ValueError,TypeError):pass
 else:raise AssertionError(bad)
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

    def test_promoted_optional_scalar_validation(self):
        self.source.write_text(SOURCE.replace('scale,status)', 'scale,factor,status)').replace(
            '    integer,intent(out),optional :: status',
            '    real(c_double),intent(in),optional :: factor !! scale factor\n    integer,intent(out),optional :: status'),encoding='ascii')
        policy=self.config['subroutines']['draw']
        policy['required_optional_scalars']=['factor']
        self.config['argument_docs']['factor']='scale factor'
        outputs=generate.generate(self.config,self.root)
        row=json.loads(outputs[self.root/'coverage.json'])['procedures'][0]
        self.assertEqual(row['optional_promoted_to_required'],['factor'])
        for bad in ([],None,{},['unknown'],['factor','factor'],['scale'],['status'],['state'],['factor',1]):
            policy['required_optional_scalars']=bad
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('python'),'compiler and Python required')
    def test_compiled_state_presence_and_promoted_scalar(self):
        source=SOURCE.replace('scale,status)', 'scale,factor,shift,status)').replace(
            '    integer,intent(out),optional :: status',
            '    real(c_double),intent(in),optional :: factor !! scale factor\n'
            '    real(c_double),intent(in),optional :: shift(:) !! shape (k); offsets\n'
            '    integer,intent(out),optional :: status')
        source=source.replace('    if(any(scale==99)) then',
            '    if(present(factor)) samples=samples*factor\n'
            '    if(present(shift)) then\n'
            '        do j=1,size(samples,2)\n'
            '            samples(:,j)=samples(:,j)+shift\n'
            '        end do\n'
            '    end if\n'
            '    if(any(scale==99)) then')
        self.source.write_text(source,encoding='ascii')
        policy=self.config['subroutines']['draw'];policy['required_optional_scalars']=['factor']
        policy['optional_arrays']={'shift':'has_shift'};policy['shapes']['shift']=['k']
        self.config['argument_docs'].update(factor='scale factor',shift='offset',has_shift='presence')
        self.policy['array_groups']['optional']['defaults']={'draw':{'factor':1.}}
        outputs=self.python_outputs();outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        script=self.root/'check.py'
        script.write_text("""import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
from example._runtime import StateBuffer
lib=C.CDLL(sys.argv[2]);api=Optional(lib)
for shift in (None,[10.,20.]):
 state=StateBuffer(lib,[1,2]);r=api.draw(state,[1.,2.],factor=2.,shift=shift,n=2)
 assert r.status==0
 expected=np.array([[6.,10.],[10.,14.]])
 if shift is not None:expected+=np.array(shift)[:,None]
 np.testing.assert_array_equal(r.outputs['samples'],expected)
 assert state.state[0]==5
""",encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            lib=folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-std=f2018','-shared','-fPIC',self.source,self.root/'adapter.f90','-o',lib],folder)
            run(['python',script,self.root/'python',lib],folder)
