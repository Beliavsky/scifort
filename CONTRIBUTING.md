# Contributing to SciFort

Read `AGENTS.md` before contributing.

## Development workflow

1. Create a focused branch.
2. Add or update tests for the intended behavior.
3. Implement the smallest coherent change.
4. Update documentation and provenance.
5. Run `fpm test` with at least one supported compiler.
6. Open a pull request describing numerical methods, sources, licenses,
   compatibility effects, and test coverage.

## Definition of done

### Unattended tests and reference comparisons

Every program discovered under `test/` by `fpm test` must run without keyboard
input, required command-line arguments, or files prepared by an external script.
Standalone tests should generate their inputs or use committed fixtures.
Optional arguments used by a test to launch its own error-path subprocesses
are permitted when the no-argument parent runs unattended.

Put programs driven by Python comparison scripts in
`validation/reference_drivers/`. These programs may read standard input,
accept file paths, or consume files produced by the comparison script. Compile
and run them through their documented comparison workflow, outside `fpm test`.
Keep the comparison scripts' source paths and any source generators in sync.

Run `python tools/check_test_inputs.py` to check for explicit console reads in
automatically discovered tests. CI runs this check before building.

### Numerical contributions

A numerical contribution is not complete until:

- the mathematical parameterization is documented;
- boundary and invalid-input behavior is documented;
- central and difficult numerical regimes are tested;
- source provenance and license obligations are recorded;
- the native Fortran API follows project conventions;
- any C ABI addition has a declaration in `include/scifort.h`;
- all tests pass under the compiler used by the contributor.

## Commit and pull-request guidance

Keep commits reviewable. Separate mechanical cleanup from numerical changes.
Do not combine a new algorithm, an API redesign, broad formatting, and build
system changes in one commit unless they are inseparable.

Pull requests that adapt third-party code must identify:

- upstream project and file;
- exact release or commit;
- original language;
- copyright holder;
- original license;
- nature of the adaptation;
- verification method.

## Numerical references

SciPy may be used as one reference implementation, but it is not an infallible
oracle. Difficult cases should also use mathematical identities,
high-precision reference values, published test tables, or independent
implementations.
