# SciFort repository rules for coding agents

These instructions apply to Codex, Claude Code, and every other automated or
human contributor. Read them before changing the repository.

## Project objective

Build a trustworthy modern Fortran statistics library that:

1. has an idiomatic native Fortran API;
2. builds with the Fortran Package Manager;
3. exposes a small, stable C ABI for Python, R, MATLAB, and Octave;
4. uses permissively licensed code with complete provenance; and
5. treats numerical behavior as a documented public contract.

The immediate scope is probability distributions and their required special
functions. Do not expand into unrelated SciPy areas without an approved issue
or roadmap change.

## Priority order

When goals conflict, use this order:

1. mathematical and numerical correctness;
2. licensing and provenance correctness;
3. API and ABI stability;
4. portability and reproducibility;
5. clarity and maintainability;
6. performance;
7. breadth of features.

Never trade away correctness or provenance merely to add more functions.

## Licensing and provenance

- Original SciFort code uses `SPDX-License-Identifier: MIT`.
- Preserve every applicable third-party copyright and license notice.
- Record imported, translated, or substantially adapted code in
  `CODE_PROVENANCE.md` before merging it.
- Record full third-party license text in `THIRD_PARTY_LICENSES.md` or an
  appropriate file under `LICENSES/`.
- Do not import GPL, AGPL, or LGPL implementation code into the permissive core.
- Do not assume code is reusable merely because it is online, old, on Netlib,
  or described as public domain by an unrelated source.
- Do not remove attribution from BSD, MIT, Apache, Boost, or other permissive
  components.
- A clean implementation from a paper or mathematical specification must cite
  the paper or specification in provenance documentation.
- Never invent a source URL, commit hash, author, test result, or license.
- When provenance is uncertain, stop the import and leave a documented issue.

## AI-generated implementation policy

- Do not ask a model to reproduce a known library routine from memory and then
  label the result as an original implementation.
- When a task names an upstream routine, inspect the exact upstream source and
  license before adapting it, or implement independently from a cited public
  mathematical specification.
- Treat close structural or textual correspondence to third-party code as an
  adaptation requiring provenance and retained notices, even when an AI tool
  produced the intermediate draft.
- Review generated constants, coefficients, stopping rules, and edge cases
  against authoritative references. Never accept plausible-looking numerical
  coefficients without verification.
- Do not introduce generated vendored code, downloaded archives, or binary
  artifacts unless the task explicitly requires them and their origin and
  license are recorded.
- Never change the repository license or remove a notice merely to make an
  imported implementation appear MIT licensed.

## Fortran style

- Use lowercase Fortran.
- Use ASCII characters only in source code, comments, and identifiers.
- Use `implicit none` in every program, module, and procedure scope where it is
  permitted.
- Define the working real kind once as `dp = kind(1.0d0)` in
  `scifort_kinds`; do not use `kind=8`.
- Use explicit `public` lists and make modules `private` by default.
- Prefix public modules with `scifort_`.
- Prefer one distribution family per module.
- Use `result(...)` on functions.
- Use `newunit=` for file I/O and `i0` for unrestricted integer output.
- Put one statement on each line.
- Keep lines at or below 100 columns when practical and never exceed the
  standard free-form limit without a compelling reason.
- Avoid preprocessor use in numerical source unless portability requires it.
- Do not rely on compiler-specific name mangling, module files, or descriptors
  in the C ABI.
- Do not add global mutable state.
- Random-number APIs must use explicit state objects; never hide an RNG in a
  module variable.
- Do not use `stop` or `error stop` inside library code. Return NaN, status
  information, or a documented error object as appropriate.

## Numerical API rules

- Scalar distribution functions should be `pure elemental` when practical.
- Add separate bulk procedures when they reduce foreign-function overhead or
  permit better vectorized implementations.
- Use explicit and consistent parameterizations. Document `loc`, `scale`,
  shape parameters, support, and endpoint behavior.
- Invalid scalar parameters return a quiet NaN unless an API explicitly
  documents another behavior.
- Procedural and C ABI routines use integer status codes for call-level errors.
- Implement `sf`, `logcdf`, and `logsf` directly or through stable identities.
  Do not generally compute them as `1-cdf`, `log(cdf)`, or `log(1-cdf)` in
  cancellation-prone regions.
- Use log-domain formulas for densities and likelihoods where overflow or
  underflow is plausible.
- Handle NaN, positive infinity, negative infinity, support endpoints, and
  probability endpoints deliberately.
- Do not silently change behavior by compiler or optimization level.
- Binary64 is the initial public precision. Do not add generic multi-kind APIs
  until their ABI, testing, and maintenance costs are addressed.

## API and ABI rules

- The native Fortran API may use modules, generics, optional arguments, and
  derived types.
- The C ABI may use only interoperable scalar types, explicit lengths,
  assumed-size arrays, caller-owned memory, opaque `c_ptr` handles, and stable
  symbol names.
- Do not expose allocatable arrays, assumed-shape descriptors, Fortran
  `logical`, compiler-dependent character data, or noninteroperable derived
  types through `bind(c)`.
- Existing public Fortran names and C symbols are compatibility commitments.
  Change them only with a documented deprecation or major-version decision.
- Keep language-specific conventions such as NumPy broadcasting, R `NA`
  handling, and MATLAB dimensions in their wrapper layers.

## Testing requirements

Every numerical change must add or update tests covering the relevant subset
of:

- ordinary reference values;
- support and probability endpoints;
- invalid parameters;
- NaN and infinity behavior;
- extreme tails;
- monotonicity and nonnegativity;
- CDF/PPF and SF/ISF round trips;
- scalar and array agreement;
- identities between related distributions;
- native Fortran and C ABI agreement.

Use absolute, relative, and tail-appropriate error measures. Do not use one
blanket tolerance for every regime.

Before reporting completion, run at least:

```text
fpm test
```

When available, also run a debug configuration with bounds and runtime checks.
Report the exact compiler and commands used. If a command could not be run,
say so plainly.

## Change workflow

1. Read the relevant implementation, tests, documentation, and provenance.
2. Make the smallest coherent change that satisfies the task.
3. Update tests before or with implementation changes.
4. Update API documentation, changelog, and provenance when applicable.
5. Run the full test suite.
6. Inspect the diff for accidental generated files, build artifacts, copied
   notices, API changes, and non-ASCII text.
7. Summarize what changed, what was tested, and any remaining numerical or
   portability limitations.

Do not rewrite unrelated files, perform broad formatting churn, or add a new
dependency without a clear technical and licensing justification.
