# Studentized-range censored finite-df extension through 100

Date: 2026-11-18.

This extension permits real 1<=df<=100 for censored studentized-range
likelihoods, NNLFs and fits, retaining 2<=k<=10 and existing affine limits.
The df>=100000 asymptotic switch is unchanged. At df=100 the code continues
to integrate the finite-df chi-scale mixture. The gap 100<df<100000 remains
outside the supported domain. See the [API contract](../docs/censored_fitting.md#studentized-range-numerical-domain).

The existing positive GL16 log-scale integration has maximum panel width 0.25.
Independent checks determine whether that resolution remains adequate as the
chi-scale density becomes more concentrated. No new algorithm or coefficients
are assumed necessary merely to raise the documented limit.

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Python 3.13.3, SciPy 1.15.3, mpmath 1.3.0.

```cmd
C:\Programs\Python313\python.exe validation/check_censored_final.py
fpm test
```

The checker prints all compile and link commands. Both -O0 and -O3 builds use
-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace. They compile
empirical_distribution, censored, fit and stats and link the comparison driver
and native tests using the existing local FPM module/archive cache. Cache paths
are printed in the commands and may require adjustment for other compilers.

New references include df=31,50,75,100 for the exact folded-t identity and
(k,df)=(2.5,31),(5,50),(10,100) for ordinary/adaptive range checks. Tiny CDFs
use independent lower-power coefficients; very remote tails use multiprecision
incomplete-beta identities. The native tests also exercise monotonicity,
complementarity, adjacent intervals, rejection of df=101 and fixed/free fits.

For the mixed sample exact [.5,1,2], left [.8], right [3], interval [1.2,1.4],
an independent folded-t scale optimization at k=2, df=100 gives scale
1.3538868383774807 and NNLF 9.32101963355031. Native fit checks use parameter
tolerance 3e-6 and scaled objective tolerance 3e-10.

The existing integration resolution is retained. A tight positive adaptive
reference at q=10, k=10, df=100 gives log density -17.252824406403576;
the native value was -17.25282440640359. SciPy's default value was
-17.252823669981833. The tighter adaptive oracle, not a relaxed tolerance,
is used for the newly added small-density comparisons at q=10.

Both -O0 and -O3 builds passed all 277 native checks and 395 reference
comparisons. The maximum comparison error over all reference types was
3.04e-9; each met its case-specific criterion: 3e-10 scaled log error for
independent multiprecision/adaptive references and 3e-8 for ordinary SciPy
comparisons. Tiny complementary multiprecision logs use relative error.
Invalid inputs and endpoint semantics use exact checks.

Commands run:

```cmd
C:\Programs\Python313\python.exe validation/check_censored_final.py
C:\Programs\Python313\python.exe validation/check_censored_final.py --reuse
```

The first run identified the default SciPy small-density reference issue
above; the corrected tighter oracle passed in both profiles on the second
run. The native finite-df, tail, interval and fit checks passed throughout.
The final full `fpm test` run passed with exit code 0 in 4382.913 seconds
(73 minutes 3 seconds), including all 277 native checks and existing fit and
distribution regressions. Existing array-temporary runtime warnings were emitted
by unrelated interpolation, optimization and spatial tests. The full suite used
the normal FPM configuration; focused -O0/-O3 builds supplied runtime checks
and floating-point traps. No numerical code changed after the successful run;
only documentation, including dummy-argument comments, was finalized.

Final audit passed: ASCII and standard line lengths in changed Fortran,
Python syntax, changed relative documentation links, consistent 109-family
inventory, and `git diff --check`. No generated fixtures, compiled artifacts,
new dependencies, or license-notice changes are included.
