# SPDX-License-Identifier: MIT
"""Explicit opt-in scalar adapters retain paired outputs and numerical NaNs."""
import copy
import json
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
public :: pair
contains
pure subroutine pair(gain, first, second)
 real(c_double), intent(in), optional :: gain !! multiplier, default 1
 real(c_double), intent(out) :: first !! multiplier plus one
 real(c_double), intent(out) :: second !! multiplier squared
 first=2
 second=1
 if(present(gain)) then
  first=gain+1
  second=gain*gain
 end if
end subroutine
end module
"""

class ScalarOutputTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        (self.root/'array_ops.f90').write_text(SOURCE,encoding='ascii')
        self.config=array_fixture_config()
        self.config['subroutines']={'pair':dict(symbol='example_pair',dimensions=[],shapes={},
            status_argument='status',invalid_status=1,description='Return paired values.',
            scalar_only=True,synthetic_status=True,shape_contract=True,
            required_optional_scalars=['gain'])}
        self.config['argument_docs'].update(first='multiplier plus one',second='multiplier squared')

    def test_contract_and_rejection_of_ambiguous_dimensions(self):
        outputs=generate.generate(self.config,self.root)
        record=json.loads(outputs[self.root/'array_coverage.json'])['procedures'][0]
        self.assertEqual(record['interface'],'scalar_subroutine')
        self.assertEqual(record['dimensions'],{})
        header=outputs[self.root/'array_api.h']
        self.assertNotIn('Invalid extents',header)
        self.assertIn('double * arg_first',header)
        self.assertEqual(outputs,generate.generate(self.config,self.root))
        original=copy.deepcopy(self.config)
        for update in ({'scalar_only':False},{'scalar_only':1},{'scalar_only':None},
                       {'dimensions':[dict(name='n',minimum=0)]},{'shapes':{'first':['n']}}):
            self.config=copy.deepcopy(original)
            self.config['subroutines']['pair'].update(update)
            with self.assertRaises(ValueError): generate.generate(self.config,self.root)
        self.config=original
        del self.config['subroutines']['pair']['scalar_only']
        with self.assertRaises(ValueError): generate.generate(self.config,self.root)

    def test_checked_and_optimized_c_clients(self):
        if not shutil.which('gfortran') or not shutil.which('gcc'): self.skipTest('compilers required')
        for path,text in generate.generate(self.config,self.root).items():
            path.write_text(text,encoding='ascii')
        client=self.root/'client.c'
        client.write_text("""#include "array_api.h"
#include <math.h>
int main(void) {
 double first=99,second=98; int status=-1;
 example_pair(3,&first,&second,&status);
 if(status || first!=4 || second!=9) return 1;
 example_pair(NAN,&first,&second,&status);
 return status || !isnan(first) || !isnan(second);
}
""",encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'],['-O3'])):
            folder=self.root/str(index); folder.mkdir()
            run(['gfortran',*flags,'-std=f2018','-fimplicit-none','-Werror=implicit-interface',
                 '-c',self.root/'array_ops.f90',self.root/'array_api.f90'],folder)
            run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',self.root,'-c',client],folder)
            exe=folder/'client.exe'
            run(['gfortran','array_ops.o','array_api.o','client.o','-o',exe],folder)
            run([exe],folder)
