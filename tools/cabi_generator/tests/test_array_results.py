# SPDX-License-Identifier: MIT
"""Independent result-array and finite dimension-choice generation tests."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from test_generator import run
import generate

SOURCE = """module example_results
use iso_c_binding, only: c_double
implicit none
private
public :: doubled, table, score
contains
function doubled(v) result(y)
    !! @dimension n = size(v)
    !! @require n >= 1
    implicit none
    real(c_double), intent(in) :: v(:) !! shape (n); vector
    real(c_double) :: y(size(v)) !! shape (n); doubled vector
    y=2*v
end function
function table(v) result(y)
    !! @dimension n = size(v)
    !! @require n >= 1
    implicit none
    real(c_double), intent(in) :: v(:) !! shape (n); vector
    real(c_double) :: y(size(v),size(v)) !! shape (n,n); result matrix
    integer :: i,j
    do j=1,size(v)
        do i=1,size(v)
            y(i,j)=v(i)+10*v(j)
        end do
    end do
end function
function score(v,b) result(y)
    !! @dimension n = size(v)
    !! @dimension m = size(b)
    !! @require n >= 1
    !! @require m >= 1
    !! @require m in (n,n-1)
    implicit none
    real(c_double), intent(in) :: v(:) !! shape (n); vector
    real(c_double), intent(in) :: b(:) !! shape (m); full or shortened vector
    real(c_double) :: y
    y=sum(v)+sum(b)
end function
end module
"""


class ArrayResultGeneratorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root/'results.f90'
        self.source.write_text(SOURCE,encoding='ascii')
        self.config = dict(schema_version=2,source='results.f90',module='example_results',
            kind_module='iso_c_binding',kind_name='c_double',wrapper_module='example_results_c',
            header_guard='EXAMPLE_RESULTS_H',functions={},excluded={},
            argument_docs=dict(n='vector length',m='second length',v='input vector',b='second vector',
                               value='result',status='call status'),
            outputs=dict(fortran='adapter.f90',header='adapter.h',coverage='coverage.json',shape_documentation='shapes.md'),
            array_functions={})
        for name in ('doubled','table','score'):
            entry = dict(symbol='example_'+name,dimensions=[dict(name='n',minimum=1)],shapes=dict(v=['n']),
                status_argument='status',invalid_status=1,result_argument='value',optional_policy='required',
                description='Compute '+name+'.',shape_contract=True)
            if name=='score':
                entry['dimensions'].append(dict(name='m',minimum=1))
                entry['shapes']['b']=['m']
                entry['length_choices']=dict(m=dict(dimension='n',offsets=[0,-1]))
            else:
                entry['shapes']['value']=['n','n'] if name=='table' else ['n']
            self.config['array_functions'][name]=entry

    def test_result_shapes_contracts_and_determinism(self):
        artifacts = generate.generate(self.config,self.root)
        self.assertEqual(artifacts,generate.generate(self.config,self.root))
        report = json.loads(artifacts[self.root/'coverage.json'])
        record = next(r for r in report['procedures'] if r['procedure']=='table')
        self.assertEqual(record['native_result_extents'],['size(v)','size(v)'])
        self.assertIn('| `value` | `(size(v),size(v))` | `(n,n)` |',artifacts[self.root/'shapes.md'])
        self.assertIn('m in (n,n-1)',artifacts[self.root/'shapes.md'])
        self.assertIn('m /= n - 1_c_int',artifacts[self.root/'adapter.f90'])

    def test_unsupported_or_mismatched_result_extents(self):
        for source in (SOURCE.replace('y(size(v),size(v))','y(size(v)+1,size(v))'),
                       SOURCE.replace('y(size(v),size(v))','y(0:size(v)-1,size(v))'),
                       SOURCE.replace('y(size(v),size(v))','y(size(v),size(v),size(v))'),
                       SOURCE.replace('y(size(v),size(v))','y(2,size(v))'),
                       SOURCE.replace('y(size(v),size(v))','y(size(v,2),size(v))'),
                       SOURCE.replace('!! shape (n,n); result matrix','!! shape (n); wrong')):
            self.source.write_text(source,encoding='ascii')
            with self.assertRaises(ValueError):
                generate.generate(self.config,self.root)

    def test_invalid_choice_policies_and_source_drift(self):
        original = copy.deepcopy(self.config)
        for rule in (dict(dimension='n',offsets=[1]),dict(dimension='n',offsets=[0,0]),
                     dict(dimension='m',offsets=[0]),dict(dimension='unknown',offsets=[0]),
                     dict(dimension='n',offsets=[-2147483648]),dict(dimension='n',offsets=[True])):
            self.config = copy.deepcopy(original)
            self.config['array_functions']['score']['length_choices']['m']=rule
            with self.assertRaises(ValueError):
                generate.generate(self.config,self.root)
        self.config=original
        for source in (SOURCE.replace('m in (n,n-1)','m in (n,n-2)'),
                       SOURCE.replace('m in (n,n-1)','m in (n,n+1)'),
                       SOURCE.replace('m in (n,n-1)','m in (n,__import__(os))'),
                       SOURCE.replace('    !! @require m in (n,n-1)\n','')):
            self.source.write_text(source,encoding='ascii')
            with self.assertRaises(ValueError):
                generate.generate(self.config,self.root)

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc'),'GNU compilers required')
    def test_checked_optimized_c_buffers_and_choice_guards(self):
        for path,text in generate.generate(self.config,self.root).items():
            path.write_text(text,encoding='ascii')
        client = self.root/'client.c'
        client.write_text("""#include "adapter.h"
int main(void) {
 double v[2]={2,3}, b[3]={1,4,5}, result[4]={99,99,99,99}; int status;
 example_doubled(2,v,result,&status);
 if (status || result[0]!=4 || result[1]!=6) return 1;
 example_table(2,v,result,&status);
 if (status || result[0]!=22 || result[1]!=23 || result[2]!=32 || result[3]!=33) return 2;
 result[0]=99; example_table(0,v,result,&status);
 if (status!=1 || result[0]!=99) return 3;
 example_score(2,1,v,b,result,&status);
 if (status || result[0]!=6) return 4;
 example_score(2,2,v,b,result,&status);
 if (status || result[0]!=10) return 5;
 result[0]=99; example_score(2,3,v,b,result,&status);
 if (status!=1 || result[0]!=99) return 6;
 example_score(1,0,v,b,result,&status);
 return status!=1 || result[0]!=99;
}
""",encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            run(['gfortran',*flags,'-std=f2018','-Werror=implicit-interface','-c',self.source,self.root/'adapter.f90'],folder)
            run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',self.root,'-c',client],folder)
            exe=folder/'client.exe'
            run(['gfortran','results.o','adapter.o','client.o','-o',exe],folder)
            run([exe],folder)
