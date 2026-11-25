# SPDX-License-Identifier: MIT
"""Independent integer scalar/array input, bounds and widened sum guards."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from test_generator import run
import generate

SOURCE="""module example_counts
use iso_c_binding, only: c_double
implicit none
private
public :: score, shifted
contains
function score(n,counts,weights) result(y)
    !! @dimension k = size(counts)
    !! @require k >= 1
    implicit none
    integer, intent(in) :: n !! scalar count
    integer, intent(in) :: counts(:) !! shape (k); integer counts
    real(c_double), intent(in) :: weights(:) !! shape (k); weights
    real(c_double) :: y
    y=real(n,c_double)+sum(real(counts,c_double)*weights)
end function
function shifted(n,counts) result(y)
    !! @dimension k = size(counts)
    !! @require k >= 1
    implicit none
    integer, intent(in) :: n !! scalar count
    integer, intent(in) :: counts(:) !! shape (k); integer counts
    real(c_double) :: y(size(counts)) !! shape (k); shifted values
    y=real(counts,c_double)+real(n,c_double)
end function
end module
"""


class IntegerInputTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.source=self.root/'counts.f90'
        self.source.write_text(SOURCE,encoding='ascii')
        self.config=dict(schema_version=2,source='counts.f90',module='example_counts',
            kind_module='iso_c_binding',kind_name='c_double',wrapper_module='example_counts_c',
            header_guard='EXAMPLE_COUNTS_H',functions={},excluded={},
            outputs=dict(fortran='adapter.f90',header='adapter.h',coverage='coverage.json'),
            argument_docs=dict(k='length',n='count',counts='integer counts',weights='weights',value='result',status='status'),
            array_functions={})
        for name in ('score','shifted'):
            shapes=dict(counts=['k'])
            if name=='score':shapes['weights']=['k']
            else:shapes['value']=['k']
            self.config['array_functions'][name]=dict(symbol='example_'+name,
                dimensions=[dict(name='k',minimum=1)],shapes=shapes,status_argument='status',
                invalid_status=1,result_argument='value',optional_policy='required',
                description='Evaluate '+name+'.',shape_contract=True,
                integer_limits=dict(n=dict(minimum=-8,maximum=8)),integer_sums=dict(counts=10))

    def test_types_metadata_and_guard_order(self):
        outputs=generate.generate(self.config,self.root)
        self.assertEqual(outputs,generate.generate(self.config,self.root))
        header=outputs[self.root/'adapter.h'];source=outputs[self.root/'adapter.f90']
        self.assertIn('int arg_n, const int * arg_counts',header)
        self.assertIn('integer(c_int), value, intent(in) :: n',source)
        self.assertIn('sum(abs(int(counts,c_int64_t)))',source)
        self.assertLess(source.index('if (k <'),source.index('sum(abs('))
        self.assertLess(source.index('sum(abs('),source.index('value=cg_native_0'))
        row=json.loads(outputs[self.root/'coverage.json'])['procedures'][0]
        self.assertEqual(row['integer_sums'],dict(counts=10))

    def test_rejects_invalid_integer_policies(self):
        original=copy.deepcopy(self.config)
        for rules in (dict(integer_limits=dict(weights=dict(minimum=0,maximum=3))),
                      dict(integer_limits=dict(n=dict(minimum=8,maximum=2))),
                      dict(integer_limits=dict(n=dict(minimum=0,maximum=2147483648))),
                      dict(integer_limits=dict(n=dict(minimum=False,maximum=2))),
                      dict(integer_sums=dict(n=10)),dict(integer_sums=dict(counts=-1)),
                      dict(integer_sums=dict(counts=True))):
            self.config=copy.deepcopy(original)
            self.config['array_functions']['score'].update(rules)
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    def test_integer_kinds_and_outputs_stay_restricted(self):
        for source in (SOURCE.replace('integer, intent(in) :: counts','integer(kind=8), intent(in) :: counts'),
                       SOURCE.replace('integer, intent(in) :: counts','integer, intent(inout) :: counts'),
                       SOURCE.replace('integer, intent(in) :: counts','integer, intent(out) :: counts')):
            self.source.write_text(source,encoding='ascii')
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc'),'GNU compilers required')
    def test_checked_optimized_c_integer_calls(self):
        for path,text in generate.generate(self.config,self.root).items():path.write_text(text,encoding='ascii')
        client=self.root/'client.c'
        client.write_text("""#include "adapter.h"
#include <limits.h>
int main(void) {
 int counts[2]={2,-3},status; double weights[2]={2,4},value[2]={99,99};
 example_score(2,5,counts,weights,value,&status);
 if (status || value[0]!=-3) return 1;
 example_shifted(2,5,counts,value,&status);
 if (status || value[0]!=7 || value[1]!=2) return 2;
 value[0]=99; example_score(2,9,counts,weights,value,&status);
 if (status!=1 || value[0]!=99) return 3;
 counts[0]=INT_MIN; example_score(2,5,counts,weights,value,&status);
 if (status!=1 || value[0]!=99) return 4;
 counts[0]=INT_MAX;counts[1]=INT_MAX;
 example_shifted(2,5,counts,value,&status);
 return status!=1 || value[0]!=99;
}
""",encoding='ascii')
        for i,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace','-ftrapv'],['-O3'])):
            folder=self.root/str(i);folder.mkdir()
            run(['gfortran',*flags,'-std=f2018','-Werror=implicit-interface','-c',self.source,self.root/'adapter.f90'],folder)
            run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',self.root,'-c',client],folder)
            exe=folder/'client.exe';run(['gfortran','counts.o','adapter.o','client.o','-o',exe],folder)
            run([exe],folder)
