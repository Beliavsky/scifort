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
