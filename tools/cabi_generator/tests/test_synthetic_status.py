# SPDX-License-Identifier: MIT
"""Status-free native subroutines keep numerical outputs and gain call status."""
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
use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
implicit none
private
public :: row_sums
contains
pure subroutine row_sums(a, sums)
    !! @dimension nrows = size(a,1)
    !! @dimension ncols = size(a,2)
    !! @require nrows >= 1
    !! @require ncols >= 1
    implicit none
    real(c_double), intent(in) :: a(:,:) !! shape (nrows,ncols); input matrix
    real(c_double), intent(out) :: sums(:) !! shape (nrows); row sums
    sums=sum(a,dim=2)
    if (a(1,1)<0) sums=ieee_value(0.0_c_double,ieee_quiet_nan)
end subroutine
end module
"""


class SyntheticStatusTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'array_ops.f90').write_text(SOURCE, encoding='ascii')
        self.config = array_fixture_config()
        policy = self.config['subroutines']['row_sums']
        policy.update(synthetic_status=True, shape_contract=True)
        policy['shapes']['sums'] = policy['shapes'].pop('y')
        self.config['argument_docs']['sums'] = 'row sums'
        for dim in policy['dimensions']:
            dim['minimum'] = 1

    def test_metadata_and_status_not_forwarded(self):
        outputs = generate.generate(self.config, self.root)
        record = json.loads(outputs[self.root/'array_coverage.json'])['procedures'][0]
        self.assertTrue(record['synthetic_status'])
        self.assertEqual(record['arguments'], ['a', 'sums', 'status'])
        source = outputs[self.root/'array_api.f90']
        self.assertIn('status=0_c_int', source)
        self.assertNotIn('status=status', source)
        self.assertIn('native numerical NaN/infinity is retained', outputs[self.root/'array_api.h'])
        self.assertEqual(outputs, generate.generate(self.config, self.root))

    def test_rejects_ambiguous_or_missing_status_policies(self):
        original = copy.deepcopy(self.config)
        for update in ({'synthetic_status': 1}, {'synthetic_status': None},
                       {'synthetic_status': 'true'}, {'status_argument': 'sums'},
                       {'status_argument': 'a'}, {'synthetic_status': False}):
            self.config = copy.deepcopy(original)
            self.config['subroutines']['row_sums'].update(update)
            with self.assertRaises(ValueError):
                generate.generate(self.config, self.root)
        self.config = copy.deepcopy(original)
        del self.config['subroutines']['row_sums']['synthetic_status']
        with self.assertRaises(ValueError):
            generate.generate(self.config, self.root)

    def test_checked_and_optimized_c_clients(self):
        if not shutil.which('gfortran') or not shutil.which('gcc'):
            self.skipTest('gfortran and gcc required')
        for path, text in generate.generate(self.config, self.root).items():
            path.write_text(text, encoding='ascii')
        client = self.root/'client.c'
        client.write_text('''#include "array_api.h"
#include <math.h>
int main(void) {
 double a[4]={1,2,3,4}, sums[2]={99,98}; int status=-1;
 example_row_sums(2,2,a,sums,&status);
 if(status || sums[0]!=4 || sums[1]!=6) return 1;
 sums[0]=99; sums[1]=98;
 example_row_sums(0,2,a,sums,&status);
 if(status!=1 || sums[0]!=99 || sums[1]!=98) return 2;
 example_row_sums(2,-1,a,sums,&status);
 if(status!=1 || sums[0]!=99 || sums[1]!=98) return 3;
 a[0]=-1; example_row_sums(2,2,a,sums,&status);
 return status || !isnan(sums[0]) || !isnan(sums[1]);
}
''', encoding='ascii')
        for index, flags in enumerate((['-O0', '-g', '-fcheck=all', '-fbacktrace'], ['-O3'])):
            folder = self.root/str(index)
            folder.mkdir()
            run(['gfortran', *flags, '-std=f2018', '-fimplicit-none', '-Werror=implicit-interface',
                 '-c', self.root/'array_ops.f90', self.root/'array_api.f90'], folder)
            run(['gcc', '-std=c11', '-Wall', '-Wextra', '-Werror', '-I', self.root,
                 '-c', client], folder)
            exe = folder/'client.exe'
            run(['gfortran', 'array_ops.o', 'array_api.o', 'client.o', '-o', exe], folder)
            run([exe], folder)
        if shutil.which('g++'):
            run(['g++', '-std=c++17', '-Wall', '-Wextra', '-Werror', '-I', self.root,
                 '-fsyntax-only', client], self.root)
