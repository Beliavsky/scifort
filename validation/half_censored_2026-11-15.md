# Half-Cauchy and half-logistic censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
References: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical checks

All 288 new native checks pass at both `-O0` and `-O3`, with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` in both builds.
They cover both families, mixed censoring, negative likelihood, fixed/free
location and scale, fit objective and iteration-limit forwarding, unchanged
array fit/likelihood calls and scalar-kernel agreement, and positive-side
identities with the Cauchy/logistic parent censored likelihoods.

Support tests include finite density at the endpoint, zero CDF and unit SF
at/below support, impossible observations, clipping finite/infinite intervals,
equal endpoints as exact observations, and whole-support intervals. Numerical
stress tests include the smallest interval width, scaled-width underflow,
adjacent far-tail intervals, overflowing raw displacement/width, tiny boundary
CDFs, near-one SF/CDF logs, almost-full-support intervals, and unrepresentable
exponential log densities and interval masses. Invalid parameters/data, empty
data, invalid/no-finite-objective fit statuses and flat one-sided fits are
also checked. The half-logistic endpoint density is evaluated directly.

Near-one probabilities are tested by ratios of their tiny nonzero logarithms,
so rounded zero or a small positive log probability fails. These checks use
`3e-13` relative tolerance, including almost-full-support intervals with both
lower and upper excluded mass. Ordinary scalar agreement uses `8e-15`
absolute tolerance; mixed likelihoods use `4e-14`; log-scale extremes allow
up to `8e-13`. Fixed-location fitted scales allow `4e-4` for optimizer
variation. No blanket tolerance is used for all numerical regimes.

## Independent references

The common mixed sample has exact values `[0.5,1.2,2.5]`, left limit `[0.3]`,
right limit `[3.5]`, and intervals `[1.1,1.4]`, `[4.0,4.4]`. At location 0
and scale 1.3, mpmath gives log-likelihoods `-14.499835370067348` for half-Cauchy
and `-14.320770635057839` for half-logistic. SciPy's fixed-location censored
fits return scales `1.4831054687500016` and `1.575110174798021`, respectively.
Native tests also check free location/scale fits converge and improve the
specified fixed candidate, without claiming a global solution.

The reference harness passes all 1,403 comparisons at checked `-O0` and
optimized `-O3`. Of these, 230 are new: 100 mixed samples per family plus
17 half-Cauchy and 13 half-logistic interval stress cases. Random samples
include positive intervals and intervals crossing support. The seed is
20261112. Half-family references use 700-digit arithmetic and independently
evaluate arctangent/tanh CDF differences, not the production parent interval
identities or wide-interval branch.

The maximum scaled error for the new families is `3.35e-15`. Across all
families it remains `1.08e-12`, in an existing gamma stress case. The harness
scales by `max(1,abs(reference))`, with threshold `2e-11`; the native relative
near-one checks above cover the tiny logs this absolute scale does not resolve.

Commands (CMD; matching local FPM build hashes are required):

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The harness explicitly compiles the changed censored module and driver before
linking unchanged archive dependencies. Its configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` and
`-O3 -fcheck=bounds`.

## Native builds and guide example

Create `build\half-native` before compiling:

```cmd
mkdir build\half-native
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\half-native -Ibuild\half-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_half_families.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\half-native\test.exe
build\half-native\test.exe
```

The optimized command uses `-O3` and `build\half-optimized`, retaining the
same runtime checks and traps. The new guide example is wrapped in a program,
compiled with the checked configuration against the changed modules, and run.

## Scope and limits

The change adds family-specific callbacks and fit overloads. Generic
censoring accumulation, the parent kernels, all existing scalar/array APIs,
and the bounded optimizer are unchanged. Fits are local box-constrained
searches; success does not guarantee global optimality, identifiability or
existence of an unconstrained finite MLE. Support-endpoint data can favor the
lower scale bound, and all-censored samples can have flat optima or boundary
solutions. Sums must fit in binary64. Half-Cauchy callbacks preserve finite
log density and SF beyond the representable standardized range; half-logistic
returns negative infinity for unrepresentable log magnitudes. No weights,
new runtime dependencies or C ABI additions are introduced.

## Full-suite verification and final audit

The full `fpm test` run passed with exit code 0 in 2272.504 seconds
(about 37.9 minutes), including all 288 new half-family checks.
The run used the compiler above and FPM's default debug profile. The
equivalent CMD command for the captured invocation is:

```cmd
fpm test > build\half-full-test.log 2>&1
```

The actual runner used Python `subprocess.run(['fpm','test'])` with stdout
and stderr redirected to that log, preserving the native exit code. Existing
array-temporary warnings were emitted; there were no test failures. This
successful full run covers the final numerical code and tests.

`python tools\check_test_inputs.py` passed. Changed Fortran ASCII/132-column
limits, Python syntax, newly added local Markdown links, the updated inventory
counts (13 supported families and 96 remaining), and `git diff --check` passed.
The guide example compiled and ran. No build artifacts, imported third-party
code, new dependencies, or unrelated edits are included. Subsequent edits only
record these validation results.
