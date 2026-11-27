# SPDX-License-Identifier: MIT
"""Independent heap object transport, ownership, validation and close/call locking."""
import copy
import json
import shutil
import unittest
import test_optional_arrays as fixtures
from test_generator import run
import generate

SOURCE = """module example_optional
use iso_c_binding, only: c_double,c_ptr,c_loc,c_f_pointer,c_associated,c_null_ptr
implicit none
private
public :: make,read,release,combine
type :: counter
 real(c_double) :: value
end type
contains
subroutine make(x,handle,status)
 !! @dimension n = size(x)
 !! @require n >= 1
 real(c_double),intent(in) :: x(:) !! shape (n); values to sum
 type(c_ptr),intent(out) :: handle !! owning result
 integer,intent(out) :: status !! call status
 type(counter),pointer :: object
 handle=c_null_ptr
 status=2
 if(sum(x)<0) return
 allocate(object)
 object%value=sum(x)
 handle=c_loc(object)
 status=0
end subroutine
subroutine read(handle,x,y,status)
 !! @dimension n = size(x)
 !! @require n >= 0
 type(c_ptr),intent(in) :: handle !! live counter
 real(c_double),intent(in) :: x(:) !! shape (n); queries
 real(c_double),intent(out) :: y(:) !! shape (n); shifted queries
 integer,intent(out) :: status !! call status
 type(counter),pointer :: object
 y=-1
 status=5
 if(.not.c_associated(handle)) return
 call c_f_pointer(handle,object)
 y=x+object%value
 status=0
end subroutine
subroutine combine(handle,other,child,status)
 type(c_ptr),intent(in) :: handle !! first live counter
 type(c_ptr),intent(in) :: other !! second live counter, possibly identical
 type(c_ptr),intent(out) :: child !! independently owning result
 integer,intent(out) :: status !! call status
 type(counter),pointer :: first,second,result
 child=c_null_ptr
 status=5
 if(.not.c_associated(handle).or..not.c_associated(other)) return
 call c_f_pointer(handle,first)
 call c_f_pointer(other,second)
 allocate(result)
 result%value=first%value+second%value
 child=c_loc(result)
 status=0
end subroutine
subroutine release(handle,status)
 type(c_ptr),intent(inout) :: handle !! owning slot, nulled on release
 integer,intent(out) :: status !! call status
 type(counter),pointer :: object
 if(c_associated(handle)) then
  call c_f_pointer(handle,object)
  deallocate(object)
  handle=c_null_ptr
 end if
 status=0
end subroutine
end module
"""

class OpaqueHandleTests(fixtures.OptionalArrayTests):
    test_metadata_presence_forwarding_and_python_signature = None
    test_rejects_unsupported_native_and_flag_policies = None
    test_compiled_c_and_python_all_presence_combinations = None

    def setUp(self):
        super().setUp()
        self.source.write_text(SOURCE,encoding='ascii')
        self.config['argument_docs']=dict(n='count',x='input values',y='shifted values',
                                         handle='owning or borrowed object',status='call status')
        self.config['subroutines']={}
        for name in ('make','read','release'):
            policy=dict(symbol='example_'+name,dimensions=[],shapes={},status_argument='status',
                invalid_status=1,description='Counter '+name+'.',shape_contract=True,
                opaque_handles=dict(handle=dict(type='counter',release_symbol='example_release')))
            if name=='release':policy['scalar_only']=True
            else:
                policy.update(dimensions=[dict(name='n',minimum=1 if name=='make' else 0)],shapes=dict(x=['n']))
                if name=='read':policy['shapes']['y']=['n']
            self.config['subroutines'][name]=policy
        self.config['argument_docs'].update(other='second counter',child='owning independent result')
        tag=dict(type='counter',release_symbol='example_release')
        self.config['subroutines']['combine']=dict(symbol='example_combine',dimensions=[],shapes={},
            scalar_only=True,status_argument='status',invalid_status=1,
            description='Combine two counters into an independent object.',shape_contract=True,
            opaque_handles={name:copy.deepcopy(tag) for name in ('handle','other','child')})

    def test_metadata_and_policy_rejections(self):
        outputs=generate.generate(self.config,self.root)
        self.assertEqual(outputs,generate.generate(self.config,self.root))
        self.assertIn('void ** arg_handle',outputs[self.root/'adapter.h'])
        self.assertIn('type(c_ptr), value, intent(in) :: handle',outputs[self.root/'adapter.f90'])
        coverage=json.loads(outputs[self.root/'coverage.json'])
        self.assertTrue(all(p['opaque_handles']['handle']['type']=='counter' for p in coverage['procedures']))
        baseline=copy.deepcopy(self.config)
        for mapping in ({},None,{'wrong':dict(type='counter',release_symbol='example_release')},
                        {'handle':dict(type='counter',release_symbol='missing')},
                        {'handle':dict(type='counter',release_symbol='example_read')},
                        {'handle':dict(type='bad-tag',release_symbol='example_release')}):
            self.config=copy.deepcopy(baseline)
            self.config['subroutines']['make']['opaque_handles']=mapping
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        self.config=copy.deepcopy(baseline)
        self.config['subroutines']['release']['opaque_handles']['handle']['type']='different'
        with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        self.config=baseline
        for source in (SOURCE.replace('type(c_ptr),intent(in) :: handle','type(c_ptr),intent(in) :: handle(:)'),
                       SOURCE.replace('type(c_ptr),intent(in) :: handle','type(c_ptr),intent(in),optional :: handle'),
                       SOURCE.replace('type(c_ptr),intent(inout) :: handle','type(c_ptr),intent(in) :: handle')):
            self.source.write_text(source,encoding='ascii')
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc'),'GNU compilers required')
    def test_checked_optimized_c_cpp_python_ownership(self):
        outputs=self.python_outputs();outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        client=self.root/'client.c'
        client.write_text("""#include "adapter.h"
#include <stddef.h>
int main(void) {
 void *handle=NULL,*child=NULL; double x[2]={1,2},y[2]; int status;
 example_make(2,x,&handle,&status); if(status || !handle) return 1;
 example_read(2,handle,x,y,&status); if(status || y[0]!=4 || y[1]!=5) return 2;
 example_combine(handle,handle,&child,&status); if(status || !child) return 7;
 example_release(&handle,&status); if(status || handle) return 3;
 example_read(2,child,x,y,&status); if(status || y[0]!=7 || y[1]!=8) return 8;
 example_release(&child,&status); if(status || child) return 9;
 example_release(&handle,&status); if(status || handle) return 4;
 example_read(2,NULL,x,y,&status); if(status!=5 || y[0]!=-1) return 5;
 x[0]=-100; example_make(2,x,&handle,&status);
 return status!=2 || handle!=NULL;
}
""",encoding='ascii')
        script=self.root/'check.py'
        script.write_text("""import copy,ctypes as C,gc,sys,weakref,pickle
from concurrent.futures import ThreadPoolExecutor
from threading import Event
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
lib=C.CDLL(sys.argv[2]); api=Optional(lib)
other=Optional(C.CDLL(sys.argv[2]))
h=api.make([1.,2.]).outputs['handle']
np.testing.assert_array_equal(api.read(h,[10.,20.]).outputs['y'],[13.,23.])
child=api.combine(h,h).outputs['child']
assert api.read(h,[]).outputs['y'].shape==(0,)
for operation in (lambda:copy.copy(h),lambda:copy.deepcopy(h),lambda:pickle.dumps(h),
                  lambda:api.read(123,[1.]),lambda:other.read(h,[1.])):
 try:operation()
 except (TypeError,ValueError):pass
 else:raise AssertionError('unsafe ownership accepted')
saved=h._policy; h._policy=dict(saved,type='other')
try:api.read(h,[1.])
except ValueError:pass
else:raise AssertionError('wrong object type accepted')
h._policy=saved
# A concurrent close must wait until the native call has finished using the pointer.
entered,proceed,closing,closed=Event(),Event(),Event(),Event()
original=api._functions['read']
def blocked(*args):
 entered.set(); assert proceed.wait(5); return original(*args)
def close():
 closing.set();h.close();closed.set()
api._functions['read']=blocked
with ThreadPoolExecutor(max_workers=2) as pool:
 read=pool.submit(api.read,h,[10.]); assert entered.wait(5)
 release=pool.submit(close); assert closing.wait(5)
 assert not closed.wait(.05)
 proceed.set();np.testing.assert_array_equal(read.result().outputs['y'],[13.])
 release.result()
assert h.closed;h.close()
api._functions['read']=original
np.testing.assert_array_equal(api.read(child,[0.]).outputs['y'],[6.]);child.close()
try:api.read(h,[1.])
except ValueError:pass
else:raise AssertionError('closed handle accepted')
api._functions['read']=original
h=api.make([1.]).outputs['handle'];pointer=h._pointer;reference=weakref.ref(h)
del h;gc.collect();assert reference() is None and pointer.value is None
failed=api.make([-1.]);assert failed.status==2 and failed.outputs['handle'] is None
first=api.make([1.]).outputs['handle'];second=api.make([2.]).outputs['handle']
def combined(reverse):
 child=api.combine(second if reverse else first,first if reverse else second).outputs['child']
 try:return api.read(child,[0.]).outputs['y'][0]
 finally:child.close()
with ThreadPoolExecutor(max_workers=4) as pool:
 futures=[pool.submit(combined,i%2) for i in range(24)]
 assert all(f.result(timeout=10)==3 for f in futures)
first.close();second.close()
for _ in range(100):
 h=api.make([1.]).outputs['handle'];assert api.release(h).status==0 and h.closed
""",encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            lib=folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-std=f2018','-fimplicit-none','-Werror=implicit-interface',
                 '-shared','-fPIC',self.source,self.root/'adapter.f90','-o',lib],folder)
            for compiler,standard in (('gcc','c11'),('g++','c++17')):
                if not shutil.which(compiler):continue
                exe=folder/(compiler+'.exe')
                run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',self.root,client,lib,'-o',exe],folder)
                run([exe],folder)
            run(['python',script,self.root/'python',lib],folder)
