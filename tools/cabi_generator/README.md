# C ABI generator prototype

This is an original MIT-licensed, project-independent generator under development
inside SciFort. All project names, source paths, symbol names, kind names, and
argument policies come from a JSON manifest. Its independent arithmetic fixture
does not depend on SciFort. The intended extraction boundary is this directory;
SciFort configuration and integration artifacts remain in the parent project.

Version 0.3 supports private-by-default modules with explicit public exports,
directly defined functions with `result(...)`, explicit scalar `real(kind)`
argument/result declarations, and `intent(in)` inputs. Elemental functions can
also receive a bulk wrapper with one varying input and uniform scalar parameters.
The native real kind and its module are configured; conversion to/from
`c_double` is explicit. SciFort's initial integration uses binary64 throughout.
Rank-one/rank-two arrays and rank-three real output arrays can receive explicit-shape
`bind(c)` adapters, with a native integer status output or an explicitly
configured synthetic call status. Their native
implementations are unchanged; no C implementation or pointer conversion is needed.
Scalar-result functions with array inputs use the separate `array_functions`
pattern described below, including explicit promotion of optional real arrays.

The generator uses [fparser 0.2.5](https://pypi.org/project/fparser/0.2.5/)
to parse Fortran syntax. It does not parse declarations using regular expressions
or copy numerical implementations. Python 3.10+ and the pinned parser are needed
only to regenerate/check artifacts. Ordinary FPM builds use checked-in wrappers
and do not need Python or fparser.

## Regeneration

For normal development, use the unified reviewed-interface workflow:

```text
python tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json
python tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
```

This refreshes/checks C wrappers, Python artifacts, dependency build inputs,
aggregate header and combined coverage together. The individual commands below
remain useful for debugging but do not update the complete build plan.


From the SciFort repository root, in CMD:

```text
python -m pip install -r tools/cabi_generator/requirements.txt
python tools/cabi_generator/generate.py bindings/cabi_config/manifests.json --batch
python tools/cabi_generator/generate.py bindings/cabi_config/manifests.json --batch --check
```

The configuration path is relative to the current directory. `--root PATH`
sets the base for all source and output paths in the manifest. The command
produces a Fortran module, a C/C++ header, and a JSON coverage record. Output
contains no timestamps or absolute paths; the source SHA-256 uses UTF-8 with
line endings normalized to LF, so Windows and Unix checkouts agree. `--check`
exits nonzero if any artifact
is missing or differs, and never writes files. Invalid manifests and unsupported
selected interfaces are rejected before outputs are written.

`--batch` reads a nonempty JSON list of manifest paths relative to the list file.
It validates every manifest and rejects duplicate outputs/symbols or outputs that
would overwrite any source/configuration before writing. `--check` works with
single manifests and batches and never changes artifacts. File-system failures
during writing are not transactional; correct the failure and regenerate.

## Configuration and contract

See the [normal manifest](../../bindings/cabi_config/normal.json) for a complete
schema-version-1 example, and the [linear-algebra manifest](../../bindings/cabi_config/linalg.json)
for schema version 2 with explicit-shape adapters. The top-level fields are:

| Field | Meaning |
| --- | --- |
| `source`, `module` | Fortran source and native module to inspect |
| `kind_module`, `kind_name` | Module exporting the real kind used by inputs/results |
| `wrapper_module`, `header_guard` | Generated Fortran module and C include guard |
| `outputs` | Paths for `fortran`, `header`, and `coverage` artifacts |
| `argument_docs` | One ASCII description per native input argument |
| `functions` | Public functions to wrap and their C interface policies |
| `subroutines` | Optional mapping of assumed-shape subroutines to explicit-shape adapters |
| `excluded` | Every remaining explicit public export, with a reason |

A function policy requires `symbol` or `vector_symbol` (or both).
`vector_symbol` and `array_argument` request a bulk entry point for an elemental function. `finite` and `positive`
list uniform parameters checked before any output is changed. These constraints
are explicit decisions, not guesses based on argument names. Native optional
arguments require `optional_policy: "required"`: callers must supply them at
the C boundary, and no defaults are inferred from procedure implementations.
Scalar wrappers delegate numerical validation to the native procedure.

A function policy may provide `symbol` for a generated scalar wrapper,
`vector_symbol` for a generated bulk wrapper, or both. At least one is required.
For bulk-only generation, `existing_scalar_symbol` optionally records a scalar
entry point maintained elsewhere. It is mutually exclusive with `symbol`, is
never emitted or redeclared, and does not contribute to `generated_symbol_count`.
The coverage record separates this existing symbol from generated `symbols`.
This is an explicit compatibility declaration, not automatic discovery or
signature verification of an external ABI; integration tests must verify it.
Existing symbols cannot collide with symbols generated in the same batch.
See the [exponential manifest](../../bindings/cabi_config/exponential.json).

For elemental bulk functions, `array_style` selects `explicit` or `pointer`.
Schema version 2 defaults to `explicit`; version 1 retains the original `pointer`
default. Existing normal manifests explicitly select `pointer` to preserve their
published null-buffer checks. Both styles have the same C pointer-and-length
prototype, but different buffer preconditions. Explicit wrappers declare ordinary
Fortran arrays such as `real(c_double) :: x(n)` and require valid non-NULL buffers,
even for zero-length calls. Pointer wrappers validate null buffers and allow NULL
for zero-length calls.

Elemental bulk calls use `size_t`, uniform doubles, and a caller-owned integer
status. Status 0 means success; 1 means an invalid uniform parameter, invalid
pointer where checks are provided, or length exceeding the compiler's default
Fortran integer maximum. Output remains unchanged on call-level errors. Empty
calls still validate uniform parameters. Invalid individual input values retain
the native function's NaN/infinity/endpoint behavior and do not become call-level
status failures. The status pointer must always be valid. Arrays must be suitably
aligned, hold the declared number of doubles, and must not overlap.

### Assumed-shape subroutine adapters

Each `subroutines` entry supplies `symbol`, `dimensions` (an ordered list of
objects with `name` and `minimum` fields), `shapes` (array names mapped to lists of dimension names),
`status_argument`, `invalid_status`, and a one-line `description`. Thus the
linear-solver manifest maps native `l(:,:)`, `b(:)`, and `x(:)` to adapter
arguments `l(n,n)`, `b(n)`, and `x(n)` with a by-value `integer(c_int)` order.
Dimension names must differ from native argument names. Shared dimensions express
required shape relationships; bounds and output sizes are explicit policy.
Dimension-list order defines the leading C arguments and is recorded explicitly
in coverage, so reordering JSON object keys does not change the ABI.

Adapters accept explicit real scalars and assumed-shape arrays declared on the
entity, for example `real(dp), intent(in) :: a(:,:)`, plus default-integer scalar
`intent(out)` results (C `int *`). A native status argument must be
default-integer `intent(out)`. Default-integer scalar/rank-one/two array inputs
are also supported. Allocatable, pointer and callback arguments remain unsupported.

An optional `logical_inputs` object explicitly maps native scalar `logical,
intent(in)` arguments to distinct integer `false` and `true` codes. For example,
`"logical_inputs": {"allow_singular": {"false": 0, "true": 1}}` emits a required,
by-value C `int` and converts it to a local native logical. Native optional logical
inputs are supported only with this explicit mapping; their defaults are never
inferred or applied at the C boundary. Other optional argument types remain
unsupported for subroutines; `array_functions` has an explicit required-input policy.
Unrecognized flag values return the configured invalid status
before touching any output apart from status. All scalar output pointers must be
valid, non-NULL and disjoint from input and other output storage. No Fortran
logical representation is exposed through C.
Compile-time checks require the native real kind to equal `c_double` and the
native default integer kind to equal `c_int`. Arrays are passed directly to the
native procedure, without conversion buffers or wrapper allocation.

Required shapes belong in the documentation of every constrained assumed-shape
procedure, public or private, independently of C wrapper selection. Define extent
names, minimum sizes, allowed empty extents and capacity inequalities, without
inventing restrictions for unconstrained arrays. Before generating an adapter,
review these native requirements against the manifest. For example:

```fortran
    !! Solve L*x = b; n = size(b), n >= 1.
    real(dp), intent(in) :: l(:,:) !! shape (n,n); lower-triangular matrix
    real(dp), intent(in) :: b(:) !! shape (n); right-hand side
    real(dp), intent(out) :: x(:) !! shape (n); solution
```

Define dimension names and their minimum sizes in the procedure documentation,
then record matching dimensions/shapes in the manifest. The generator compares
trailing FORD comments against the reviewed manifest, accepting `shape (n,n)` or
`(n,n)` at the beginning of a comment (whitespace and case are ignored). It reads
only declarations identified by the Fortran AST, including continued declarations.

Missing, unsupported or mismatched tuples appear in the coverage record's
`native_shape_documentation_review` list. Generation and `--check`, including the
unified pipeline, print a `documentation review` advisory to stderr. An empty list
means the documented tuples match the configuration. Advisories do not fail a
build or rewrite native source; resolve them through human review. A successful
`--check` verifies artifact freshness and can still print documentation advisories.

The generator audits only procedures selected in its adapter manifests; it does
not audit the documentation of the entire native library. The contributor rule
also applies to procedures outside those manifests.

This is a deliberately narrow documentation check, not semantic proof: it does
not infer shapes from executable code, validate prose dimension definitions or
minimum sizes, or understand alternative mathematical notation. More complicated
relationships and any intentional dimension renaming need manual review.

All arrays use contiguous column-major storage. Callers must provide valid
non-NULL buffers, including for permitted empty arrays, and adequate storage for
every declared extent; these adapters cannot check buffer addresses or capacities.
Output and mutable arguments use `intent(inout)` at the adapter boundary so an
invalid-extent return preserves caller data. After validation, native status
codes and output-on-failure semantics are preserved. Numerical failures can
therefore leave native partial results. No leading dimensions, arbitrary strides,
or pointer-style subroutine adapters are implemented yet.

Every explicit public export must be generated or explicitly excluded. Public
generics, reexports, derived types, callbacks and unsupported function forms
cannot silently become scalar wrappers. Reviewed explicit-array function policies
are handled separately. This is not yet
a whole-project scanner or a full Fortran semantic resolver. Default-public
modules are rejected because complete reexport discovery needs more analysis.
Python scalar/elemental and selected explicit-array binding generation are
described below. Other interface patterns remain unsupported.

## Checks

```text
python -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/cabi_config/tests -v
```

The first command tests only the reusable generator and independent arithmetic
fixture. The second is SciFort-specific integration and remains outside this
directory when the generator is extracted into its own repository.

The tests require fparser. Compiled tests additionally require GNU Fortran and
GCC (override with `FC` and `CC`); they are explicitly skipped if unavailable.
If `g++`/`CXX` is available, the umbrella header is also checked as C++17.
Compiled tests run both `-O0 -g -fcheck=all -fbacktrace` and `-O3`, with
`-std=f2018 -fimplicit-none -Werror=implicit-interface`. A C11 client is compiled
with `-Wall -Wextra -Werror`, linked through GNU Fortran, and executed.

Coverage includes an independent arithmetic module, source discovery, unsupported
interfaces, explicit optional policies, duplicate symbols, stale generated files,
and errors that must not overwrite existing outputs. SciFort integration compares
all ten new normal-distribution symbols against direct native calls, including
endpoints, tails, NaN/infinity, invalid scales, scalar/bulk agreement, and
independent standard-library references. It also checks null/empty/oversized
bulk calls and output preservation.

The family integration test additionally compiles unchanged production scalar
wrappers for the seven benchmark families, plus direct native-call oracles.
It checks 53 scalar/bulk pairs (the five generated normal methods and all eight
methods of the six added families), endpoints, tails, invalid uniform parameters,
empty/oversized calls, output preservation, C clients and C++ declarations.
The three original normal bulk methods remain unchanged.

Windows validation of version 0.3 on 2026-11-23 passed 21 focused tests
(18 reusable tests and three SciFort integration tests), using
Python 3.13.3, fparser 0.2.5, and GNU Fortran/GCC/G++ 17.0.0 20260510
(experimental). The commands above and the generator's `--check` command were
run with `build/cabi-generator-venv/Scripts/python.exe`. Other operating systems
and compilers have not yet been validated for this generator.

These focused checks cover the generator and configured ABI integrations; the
numerical implementations, existing C symbols, and global build configuration
are unchanged. Run the full project suite before release or merge as required
by the repository rules.

## Next interface patterns

1. Additional scalar/elemental distribution families and bulk parameters.
2. Additional explicit-rank array procedures, strides, and output sizing policies.
3. Allocatable results using size-query/fill interfaces.
4. Derived types using opaque handles with defined lifetimes.
5. Callbacks, strings, and generic resolution with explicit policies.

Keep these rules reusable; project-specific mathematical constraints and ABI
compatibility decisions belong in manifests and integration tests. fparser is
an external build tool, not vendored code or a linked numerical dependency;
see [provenance](../../CODE_PROVENANCE.md) and its retained
[license text](../../LICENSES/fparser-BSD-3-Clause.txt).

## Generated Python layer

`generate_python.py` reads a separate Python policy referencing the same C
manifests. The policy maps native procedure prefixes to Python method names,
provides explicit finite numeric defaults and identifies any existing excluded
C functions. Native scalar signatures still come from fparser; parameter defaults
are never inferred. Each selected function needs scalar and bulk C symbols.

```text
python tools/cabi_generator/generate_python.py bindings/cabi_config/python.json
python tools/cabi_generator/generate_python.py bindings/cabi_config/python.json --check
```

It emits `_generated.py`, a copy of the reusable `python_runtime.py`, and
supplementary Fortran adapters for the focused standalone development build.
The runtime depends on NumPy and Python ctypes, not fparser or SciFort-specific
module names. Generated code exposes a generic `API` owner; SciFort's package
exports it as `Scifort`. All numerical functions remain native. Both scalar and
array paths, alignment, byte order, output aliasing, statuses and empty inputs
are tested in the [Python integration](../../bindings/python/README.md).

The supplements must not be linked alongside a full C ABI object defining the
same symbols. They let a small distribution-only library build without compiling
unrelated modules. The regular SciFort source and C ABI remain unchanged. Python
policy/runtime and supplement outputs have their own read-only drift check.

## API inventory and draft manifests

`discover.py` inventories explicit public exports in private-by-default modules,
including dummy types, ranks, shape declarations, optional flags and function
results. Unsupported signatures, unresolved generics/types/constants/reexports,
and default-public modules are reported explicitly. Scalar/elemental support is
checked by the same signature validator used for generation. Array subroutines
need reviewed dimension/status policies and are not automatically drafted as
scalar functions; configured array adapters are counted separately.

```text
python tools/cabi_generator/discover.py src/stats/gamma.f90 src/stats/beta.f90 src/stats/student_t.f90 --abi-source src/c_api --manifest-batch bindings/cabi_config/manifests.json --python-policy bindings/cabi_config/python.json --output build/cabi-discovery/inventory.json
```

Inputs may be files or directories (recursive `.f90` discovery); `--abi-source`
can be repeated. JSON records full argument details and reasons. Sibling Markdown
shows separate counts for native exports, supported scalar/elemental patterns,
observed C function/procedural bindings, configured scalar/bulk wrappers and
Python methods. Counts apply only to inspected native files. Parser/module-access
failures make the report incomplete and produce a nonzero exit status. Large
data-initialization modules can take a long time in fparser; prefer a bounded set
of numerical API files when reviewing a family. `--max-source-bytes N` explicitly
skips larger native files, records their names/reasons and returns a nonzero
incomplete report; the default has no size limit. It does not limit ABI-source
parsing. `--progress` prints each native filename before parsing. The demonstrated
complete report covers the three added families, not the whole repository.

Observed C symbols come from `bind(c)` module procedures, linked to native
exports through explicit `use ..., only:` imports, including renames. Both
positional and keyword function references are recognized. These are syntactic
links: no interprocedural tracing, implicit imports, external procedures or ABI
compatibility proof is claimed. Configured coverage is read from explicit C
manifests and Python policy; it does not assert that generated files, installed
packages or loaded libraries exist. Other integer/array/callback patterns remain
explicitly outside scalar/elemental draft generation.

To create drafts in a new directory:

```text
python tools/cabi_generator/discover.py src/stats/gamma.f90 src/stats/beta.f90 src/stats/student_t.f90 --kind-module scifort_kinds --symbol-prefix scifort_ --draft-directory build/cabi-discovery/new-drafts --output build/cabi-discovery/draft-inventory.json
```

`--kind` defaults to `dp`; the kind module and symbol prefix must be supplied for
drafting. Each draft includes `draft: true` and review instructions, which make
it invalid input to the wrapper generator. Argument meanings, varying inputs,
optional handling, finite/positive constraints and Python defaults remain
unresolved. Symbol/output names are proposals, not API commitments. Existing
C symbols must be reviewed before retaining or replacing proposed scalar symbols.

Review and replace every placeholder, choose exclusions and output paths, then
remove the draft/review metadata. The gamma/beta/t manifests show the result
of that process. Existing drafts are never overwritten; inputs and referenced
policies are protected from report writes. Normal filesystem failures during
report writing are not transactional. The tool generates no numerical code.

## Unified reviewed interface workflow

`pipeline.py` reads [pipeline.json](../../bindings/cabi_config/pipeline.json),
which names the C batch, Python policy, source search directories, build-plan
output, coverage output and aggregate-header output/guard. Project names and
mathematical policies remain in configuration. The reusable tool does not contain
a family list or a special case for the chi-square distribution.

The command validates and computes all artifacts before writing, and leaves
unchanged outputs untouched. `--check` is
read-only and reports missing/stale C, Python, header, build-plan or coverage
artifacts. Invalid dependencies, source/config/template collisions or unresolved
review policies fail before any files change. Filesystem write failures are not
transactional; correct the failure and regenerate. Existing handwritten source
and aggregate-header targets cannot be overwritten by regeneration.

The build planner indexes MODULE statements with fparser's logical-line reader,
then parses only the required source closure. It follows module and procedure-local
USE statements, chooses providers by module declaration rather than filename,
and sorts dependencies before dependents. Multiple modules in one file must
already appear in a compilable order. Missing or ambiguous providers, dependency
cycles and forward references within a file produce explicit diagnostics.
Intrinsic USE statements do not need source providers. For implicit module nature,
a configured source provider takes precedence; otherwise the standard intrinsic
module set is recognized. NON_INTRINSIC imports always require a source provider.

The planner currently supports free-form Fortran module sources accepted by the
configured fparser grammar. Submodules, preprocessing, INCLUDE inputs, separately
compiled external modules, non-module procedure dependencies and non-Fortran link
inputs need explicit future support; they are not inferred from filenames.
Required submodule/preprocessor/INCLUDE sources are rejected. Other standalone
external procedure calls are outside module-dependency analysis and can still
produce compiler/linker errors. The source index does not fully parse unrelated
numerical bodies or large data initializers.

The checked-in [build_sources.json](../../bindings/python/build_sources.json)
records ordered files, dependency edges, module providers, Python group names
and normalized UTF-8/LF source hashes. The development builder reads that plan
without fparser and rejects changed source hashes with a regeneration diagnostic.
Consequently even a dependency-body edit requires refreshing the reviewed
pipeline before using this optional development build. Ordinary FPM builds are
unchanged. This guard prevents silently compiling against stale build metadata.

[interface_coverage.md](../../bindings/cabi_config/interface_coverage.md) and its
JSON sibling report reviewed exports, scalar C symbols, bulk C symbols, array
adapters and Python methods. They do not claim whole-project coverage or SciPy
parity. The generated aggregate header keeps `scifort.h` declarations synchronized
with the C batch without a handwritten list of family includes.

To extend the workflow: review a native module's C manifest, list it in the C
batch, and add its Python mapping/defaults when desired. Run the unified command
and focused tests. Source dependencies, header aggregation, build ordering and
coverage then update automatically. The chi-square addition exercises this path.

## Per-element parameters

A function policy may additionally specify `elementwise_symbol`. This emits a new
explicit-shape C entry point with `size_t n`, one `const double *` per native
argument in declaration order, `double *output`, and `int *status`. The native
function must be elemental. All input arrays have length n; input arrays may
share storage, but output must not overlap any input. Every buffer must be valid,
aligned and non-NULL even for n=0. The status pointer is always required.

Status 1 rejects lengths exceeding the default Fortran integer maximum and leaves
output unchanged. Status 0 includes invalid per-element parameters, whose results
follow the native scalar function. The existing `finite` and `positive` policies
apply only to uniform-parameter bulk symbols. No broadcast rules live in C.

The generated Python runtime uses this optional symbol when any parameter has
nonzero rank, broadcasts all inputs in NumPy, validates output shape and overlap
against every original input, and prepares contiguous aligned binary64 buffers.
All-scalar and uniform-parameter calls retain their existing paths and contracts.
All eleven reviewed distribution configurations provide this pattern for all eight methods.
Regenerate with the unified pipeline so the development build hashes and combined
coverage also include the new symbols.

## Native shape-documentation review validation (2026-11-24)

The following CMD command passed all 43 focused generator/workflow tests:

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
```

New tests cover matching/missing/mismatched tuples, rank/order differences,
continued declarations, bare tuples, case/whitespace, coverage records, CLI
advisories during generation and `--check`, and unchanged native source content.
The existing independent ABI fixtures also compiled and ran with GNU Fortran,
GCC and G++ 17.0.0 20260510 (experimental) on Windows under runtime checks and
`-O3`, using the flags documented above. The native solver and fixture changes
are comments only; no numerical suite was rerun. Other platforms remain unverified.

## Structured shape-contract prototype

The four [matrix-helper adapters](../../bindings/cabi_config/matrix_helpers.json)
pilot a SciFort-specific annotation convention. It is not FORD or Sphinx syntax.
Readable trailing `!! shape (...)` comments remain the exact array-shape contract.
Initial procedure documentation may add:

```fortran
    !! @dimension n = size(l,1)
    !! @dimension k = size(b,2)
    !! @require n >= 1
    !! @require k >= 0
    !! @when k == 0 : Right-hand sides and triangular diagonals are not evaluated.
```

Opt in with `"shape_contract": true` on a subroutine policy. Every configured
dimension then needs exactly one `@dimension` definition and one `@require`
minimum. Definitions accept `size(array,axis)` or `size(vector)`. The referenced
array axis must carry the same dimension in the manifest. Minimums and trailing
shape tuples must match that manifest exactly. Missing, conflicting, repeated or
unsupported annotations are errors before any output is written. Existing
non-opted-in adapters retain the advisory shape-comment review.

`@when dimension == integer : description` is a structured conditional behavior
note. Its dimension/value are validated, and its text is rendered in the generated
documentation. It does not introduce a conditional runtime assertion or prove the
behavior described; focused native/C tests check that behavior. This prototype
does not yet support arithmetic extent expressions, capacity inequalities,
optional-array conditions or general logical expressions. Unsupported constructs
are rejected rather than evaluated or guessed. Source annotations never replace
explicitly reviewed ABI policies, and executable source is not analyzed for them.

An optional `outputs.shape_documentation` path ending in `.md` generates a
[dimension and argument table](../../docs/matrix_helper_shapes.md), alongside
Fortran/header/coverage outputs. The source-derived contract also appears in
coverage JSON. The unified pipeline generates/checks the Markdown file, aggregates
the new header and updates reviewed-interface coverage. The Markdown can be
included in documentation builds that support Markdown tables; no automatic
FORD/Sphinx integration is introduced in this prototype.

## Python array groups

Python policy schema 1 optionally accepts `array_groups` alongside distribution
`groups`. Each array group supplies a list of reviewed C `manifests` and an explicit
`defaults` dictionary keyed by native procedure name. The SciFort `linalg` group
combines `linalg.json` and `matrix_helpers.json`, with
`"defaults": {"psd_decompose": {"allow_singular": true}}`.

The reusable `generate_arrays.py` validates native declarations, C policies and
opted-in structured contracts, then emits typed method metadata/signatures.
Dimension definitions must reference input arrays so dimensions can be determined
before allocating outputs. Exact array shapes and minimum sizes come from those
checked contracts. All selected subroutines in the supplied manifests are exposed;
names/defaults and group collisions are validated. Native scalar real inputs,
scalar logical inputs with explicit C mappings, and real/integer scalar and
rank-one/two array outputs are supported, plus rank-three real subroutine outputs. Native `intent(inout)` arguments need
a reviewed state codec; optional real input arrays need a presence policy.
Methods return `ArrayResult(status, outputs)`;
numerical statuses and partial outputs are retained.

The runtime validates shapes, C integer extent limits, real numeric dtypes,
column-major storage and buffer overlap. Inputs can be converted; caller-provided
outputs must already satisfy the exact native-endian float64 or `np.intc`
Fortran layout declared by the native output category.
Input/input, input/output and output/output overlap are conservatively rejected
against original arrays. It supplies valid disjoint buffers for zero-size arrays
and keeps all inputs/outputs and the shared-library owner alive through each call.

The pipeline includes these primary C adapter roots in the optional development
library's dependency closure and records array groups separately from distribution
families. Array adapters use their primary C implementations, with no duplicate
compatibility adapters. Combined coverage and discovery include the eight Python
matrix methods. Regenerate/check with the unified workflow so the source plan,
Python runtime copy and all contract tables stay current.

## Functions with array inputs and scalar/array results

An optional `array_functions` manifest mapping uses the explicit-array policy
fields (`symbol`, ordered `dimensions`, `shapes`, `status_argument`,
`invalid_status`, `description`, optional `logical_inputs` and `shape_contract`)
plus two required fields: `result_argument` and `optional_policy: "required"`.
The result/status names are synthetic C output arguments, distinct from all
native dummy and dimension names. Selected names must be disjoint from
`functions`, `subroutines` and `excluded`.

This pattern accepts directly defined non-bind(c), non-impure functions with an
explicit scalar or rank-one/two `real(kind)` result declaration and `result(...)`. Native inputs
must be `intent(in)` real/default-integer scalars or rank-one/rank-two assumed-shape arrays,
or scalar logicals with explicit integer mappings. Reviewed state codecs also
allow scalar derived-type intent(inout) arguments. Optional real/logical inputs
are promoted to required C inputs unless optional_arrays explicitly preserves
array absence through flags. Arrays still require declared shapes; no callbacks,
C descriptors or native algorithm changes are introduced.

The adapter is a `void` C subroutine with ordered dimensions, native inputs,
result buffer pointer and status pointer. It checks extents/flags before modifying
the result, calls the native function and sets status 0. A valid call can return
NaN/infinity: status denotes call validity, not numerical parameter validity.
All pointer/alignment/column-major/nonoverlap contracts of explicit-array adapters
apply. See the reviewed [multivariate normal policy](../../bindings/cabi_config/multivariate_normal.json).

Python `array_groups` may specify a `prefix` stripped from native method names;
default policies remain keyed by full native names. Scalar-result methods return
Python float and expose no `out` argument. Explicit defaults for optional real
arrays may be `{"fill": "zeros"}` or `{"fill": "identity"}` (square matrices only).
They generate `None` defaults and materialize native binary64 arrays after input
rank/dimension validation. Dimension definitions normally use required inputs.
A defaulted dimension source requires explicit `dimension_defaults` as described
below. Defaults are configured explicitly, not inferred from source comments,
procedure names or executable branches. Existing subroutine methods
retain `ArrayResult(status, outputs)` and their caller-buffer behavior.

## Explicit array-valued results and dimension choices

An `array_functions` policy may include the synthetic `result_argument` in
`shapes`, with rank one or two. The native result must be explicitly declared
with default lower bounds and extents exactly of the form `size(input)` (rank-one
input) or `size(input,axis)` (axis 1/2). The AST-identified extent must map to the
same reviewed dimension as the input axis. Allocatable/deferred/assumed-shape
results, literal or general expression bounds and ranks above two remain
unsupported. Native result declarations need matching trailing shape comments.
Generated shape tables show their explicit extents and identify the adapter's
result buffer. No mathematical output size is inferred from a function name.

The C adapter assigns the native function result to caller-owned explicit-shape
storage. It preserves the output on call-validation errors. The compiler/native
function may create temporary storage for an array result; caller-owned C output
does not promise an allocation-free numerical call. Python returns the ndarray,
with optional `out=ndarray` identity preservation, exact shape/native-float64/F-layout
checks and input/output nonoverlap checks. Matrix subroutines retain their existing
`out={name: buffer}` and `ArrayResult` conventions.

An optional `length_choices` policy maps a declared dimension to a finite list
of another dimension plus nonpositive integer offsets, for example:

```json
"length_choices": {"m": {"dimension": "k", "offsets": [0, -1]}}
```

Its matching source directive is `!! @require m in (k,k-1)`. Both dimensions must
be declared and distinct; offsets must be unique nonpositive C-int values.
Positive offsets and arbitrary expressions are rejected. This restricted grammar
uses no eval, source execution or inference. C guards run after extent minima,
using safe subtractions; invalid choices preserve output and return the configured
invalid status. Python checks the same choices before allocating outputs or
calling C. Structured contracts require exact source/policy agreement and render
the choices in Markdown. `@when` prose remains descriptive, not a checked assertion.
See the [reviewed Dirichlet manifest](../../bindings/cabi_config/dirichlet.json).

## Integer inputs and explicit overflow policy

Default-integer native `intent(in)` scalars map to by-value C int; rank-one/two
integer inputs map to const C-int buffers with reviewed explicit shapes. Integer
kinds other than the native default and integer inout arguments remain unsupported.
Default-integer scalar and rank-one/two array outputs map to caller-owned C-int
buffers with explicit shapes. Python allocates `np.intc` arrays initialized to zero
or accepts exact native-endian, aligned, writable, Fortran-contiguous `np.intc`
outputs. Empty outputs use correctly typed nonnull scratch buffers. Existing
scalar integer outputs are unchanged.
Direct association requires native default integer kind to match `c_int`.
Real-valued function results may depend on integer input sizes.

Python accepts integer dtype inputs, including signed/unsigned and different
endianness, checks their range, then converts arrays to aligned Fortran-contiguous
`np.intc` buffers. It rejects floats, booleans, strings and object arrays; integer
scalars also require integer dtype. Dimension checks precede element scans and
output allocation. Output arrays remain binary64. Original-input nonoverlap
checks apply even if conversion would create separate storage.

Optional `integer_limits` maps integer input names to inclusive signed-C-int
minimum/maximum values. Optional `integer_sums` maps rank-one integer input names
to a nonnegative C-int absolute-sum ceiling, for example:

```json
"integer_limits": {"n": {"minimum": -2147483648, "maximum": 2147483646}},
"integer_sums": {"m": 2147483646, "x": 2147483646}
```

These are explicit reviewed safety policies, not inferred from executable code.
C adapters check them before the native call, returning configured invalid status
and preserving results. Comparisons and absolute sums widen to `c_int64_t` before
arithmetic, including before abs(INT_MIN). Rank-one C-int extents/elements bound
the widened sum below 2**63. Python checks bounds before narrowing and uses the
same widened sum. Signed 32-bit limits are the initial supported integer contract.
Guards, bounds and sum limits appear in C headers, coverage and shape tables.

The multinomial/hypergeometric policies cap n and hypergeometric absolute sums at
2147483646 to protect existing default-integer sums and n+1/count+1 expressions.
This deliberately bounds foreign calls without changing the native numerical
API. Numerical parameter errors inside these bounds retain native NaN/zero
semantics and valid-call status 0. See the reviewed
[multinomial](../../bindings/cabi_config/multinomial.json) and
[hypergeometric](../../bindings/cabi_config/multivariate_hypergeom.json) policies.

## Optional-array dimension defaults

An array group may opt into fallback extents when a dimension-defining optional
array is omitted. Configure `dimension_defaults` by full native procedure name,
for example:

```json
"dimension_defaults": {
  "matrix_normal_entropy": {"m": 1, "n": 1}
}
```

The procedure must also configure identity/zero fills for its optional arrays.
The fallback keys must cover exactly the dimensions whose declared sources are
defaulted arrays. Values must be integers (not booleans), satisfy the declared
minimum, and fit a signed C int. Unknown procedures/dimensions, missing fallbacks
and fallbacks for required input sources are rejected. Existing policies without
defaulted dimension sources need no change.

Python uses the provided source array's extent when present, otherwise the
configured fallback. Each dimension resolves independently, before default
arrays are allocated. Another argument's size never overrides the declared
source or its explicit fallback. Normal shape, extent and alias checks still
apply. C callers always supply explicit dimensions and arrays; fallback policies
exist only in the Python wrapper layer.

## Reviewed public generic selections

A manifest may supply `generic_specifics`, mapping a public named generic to
one locally defined module procedure declared in its interface:

```json
"generic_specifics": {
  "poisson_binom_pmf": "poisson_binom_pmf_real"
}
```

The AST supplies the selected procedure's argument/result types and source shape
annotations. The emitted Fortran imports and calls the public generic, allowing
the compiler to resolve the reviewed specific from those argument types; private
specific procedures are never imported by name. Per-manifest coverage records
`native_specific`, and shape tables identify the selection. Public export
coverage remains keyed by the generic name, not its private members.

Selections must name generated public generics and actual, local module-procedure
members. Unknown/private generic names, nonmember or unresolved imported targets,
and selections for excluded procedures are rejected. Unconfigured generics remain
unsupported and require explicit exclusion. This adds no automatic overload
choice, descriptor ABI or general Fortran semantic resolver. The existing
scalar/array signature restrictions still apply to the selected procedure.

## Subroutines without native status

A reviewed `subroutines` entry may set `"synthetic_status": true`. Its
`status_argument` then names a new, distinct output appended after all native
arguments at the C boundary. The generator never forwards this argument to the
native subroutine. The opt-in must be a JSON boolean; omission or false retains
the existing native-status requirement. Array functions already synthesize a
status and reject this subroutine-only option.

Dimension, flag and integer guards return `invalid_status` and preserve all
other outputs. After a valid native call, the adapter sets status 0 even when
numerical inputs produce NaN or infinity. It does not infer numerical success
or intercept native failures. This distinction is recorded in coverage and
headers; the native Fortran signature remains unchanged.

Python returns `ArrayResult(status, outputs)`, with named scalar/array outputs
and optional caller-owned array buffers through `out={"name": array}`. The same
shape, type and alias checks apply as for native-status subroutines. See the
[Poisson-binomial manifest](../../bindings/cabi_config/poisson_binom.json) and
[derivative contract](../../bindings/python/README.md#poisson-binomial-score).

## Ordered optional dimension sources

A Python array group may configure `dimension_sources` by full procedure name,
then dimension. The ordered candidates have the same `argument`/`axis` format
as structured dimension definitions (`axis: null` for a vector's length):

```json
"dimension_sources": {
  "multivariate_normal_entropy": {
    "n": [{"argument": "mean", "axis": null}, {"argument": "cov", "axis": 1}]
  }
}
```

The first candidate must equal the declared source. Every candidate must be a
distinct configured optional/defaulted array axis with the same required shape
dimension. The dimension must also have an explicit `dimension_defaults` value.
Missing, malformed, repeated, mismatched or unused policies are rejected.
Python uses the first supplied candidate, then the fallback if none is supplied.
It validates every provided shape against that dimension; conflicting arrays
raise rather than silently changing the selected dimension. Rank and C-int
extent checks precede allocation. C still receives every input explicitly.
This policy must be reviewed against native inference order; annotations alone
do not prove native optional-argument behavior.

## Python zero-based index inputs

A Python array group may configure `zero_based_indices`, mapping full procedure
names to `{ "integer_input": "bound_dimension" }`. For example:

```json
"zero_based_indices": {
  "multivariate_normal_marginal": {"dimensions": "n"}
}
```

Only required rank-one integer input arrays qualify. The bound must be a declared
dimension. Unknown/unused policies, defaults and combinations with native integer
range/sum guards are rejected. Python checks every index in `[0,n)` before
narrowing, copies it into widened signed storage, adds one, and passes C-int
indices to C. Caller arrays remain unchanged. Invalid ranges raise `ValueError`;
floats and booleans raise `TypeError`. Empty arrays still obey normal minimum
shape rules. The C ABI, native annotations and Fortran procedures remain
one-based; uniqueness and other native constraints remain native validations.

## Presence-controlled optional arrays

A `subroutines` policy may explicitly map native optional array inputs to new
C presence flags with `optional_arrays`:

```json
"optional_arrays": {
  "fixed_mean": "has_fixed_mean",
  "fixed_cov": "has_fixed_cov"
}
```

Only real rank-one/two `optional, intent(in)` arrays are supported, in selected
subroutines and array-input functions. Optional scalar/output arrays, unknown
arguments, flag-name collisions and
non-identifier flags are rejected. The mapping must be nonempty when supplied.
Each by-value C-int flag precedes its corresponding explicit-shape array in the
C signature: 0 means absent, 1 means present. Other values return the configured
invalid status before native calls or optional-input copies, preserving outputs.
Dimension guards run first. Every C buffer remains required, nonnull, fully
sized, aligned and disjoint, including absent-array placeholders. Placeholder
values are ignored, not inspected for numerical validity.

After validating all flags, the adapter copies present optional inputs into
local allocatable arrays. It leaves absent ones unallocated and forwards those
to nonallocatable optional native dummies, preserving Fortran 2008+ absence
semantics. Descriptors remain inside Fortran; no descriptors or allocatables
cross the C boundary. This avoids generating exponentially many call branches,
at the cost of copies for present optional inputs. No numerical routine is
rewritten, and no synthetic numerical success is inferred from native status.

Coverage records `optional_arrays` separately from optional arguments promoted
to required. Shape tables retain each array's shape and name its presence flag.
Python hides these flags and exposes `None` defaults, passing None as absent
and validating supplied arrays normally. Optional arrays cannot define Python
dimensions or also configure fill/default policies. The Python wrapper allocates
full-sized uninitialized placeholders when absent; it never substitutes zeros or
identity as a native default. These placeholders and native present-input copies
add memory overhead. Native partial outputs and statuses remain unchanged.

## Reviewed serialized state codecs

A `subroutines` or `array_functions` policy may opt into `state_codecs` for scalar
`type(name), intent(inout)` arguments. Each entry declares `module`,
`type_module`, `type`, `load`, `store`, fixed positive `length`, and positive
`invalid_status`. For example the Dirichlet policy names the native `rng_state`,
SciFort's RNG bridge load/store procedures, six int64 components, and status 3.
No automatic derived-type layout inference or arbitrary structure exposure is
performed. Names and declaration type/intent are validated, while helper
signatures/behavior require review and compiled tests.

The C boundary is a caller-owned `int64_t[length]`; descriptors stay in Fortran.
After extent/flag guards, the adapter decodes state into a local native object.
A nonzero decoder status returns the configured invalid-state status with state
and outputs unchanged. Valid state is required even for empty calls. After the
native call, state is always serialized back, including on numerical failure;
this preserves native stream advancement. Independent streams require separate
buffers, with no shared mutation during a call. The codec helpers must preserve
state on rejected input and export only initialized native states.

Python accepts a library-bound `StateBuffer` (SciFort exposes `RandomState`),
checks its int64 storage, extent, ownership and buffer disjointness, and passes
that live state without copying. Snapshots are copies. States cannot have Python
defaults. A required native integer status may opt into `optional_status: true`
when its declaration is optional; C always supplies it. Other optional status
forms remain unsupported. A codec can explicitly set `promote_optional: true`
for a native optional state argument. This flag must be true when supplied and
must match an optional native declaration; omit it for required native state.
C/Python still require the state buffer and the native function sees it present.
No absent-state branch or hidden RNG is introduced.

## Output-defined dimensions

An array group's `output_dimensions` maps native procedure names to dimension
names and nonnegative C-int defaults, for example
`{"dirichlet_rvs_array": {"n_samples": 1}}`. The dimension must be declared
in the native structured contract using an output array axis; all array shapes
must still agree with the C policy. Python exposes it as a keyword-only integer
parameter, validates it before allocation and uses it for caller-owned output
shape checks. Boolean, negative, oversized, unknown and reserved dimensions are
rejected. Required input-defined dimensions continue to be inferred normally.
No dimension is inferred from an unvalidated output buffer or hidden default.

## Required optional scalar inputs

A subroutine policy can explicitly list `required_optional_scalars`, for example
`["df"]` for multivariate Student t sampling. Each name must identify a native
optional real or default-integer scalar `intent(in)`. C always supplies the
by-value argument; the native procedure always sees it as present. This policy
is distinct from presence-controlled optional arrays. It does not infer a
native default, generate a scalar presence flag, or allow None. Python defaults
must be reviewed and configured separately (Student t uses `df=1.0`).

Coverage lists these arguments in `optional_promoted_to_required`. Unknown or
repeated names, nonoptional dummies, arrays, logicals, derived types, status
outputs, other intents and function policies are rejected. Existing scalar
function promotion and other interfaces retain their behavior. See the
[multivariate sampler contract](../../bindings/python/README.md#multivariate-vector-batch-sampling).

## Rank-three real outputs

Reviewed subroutines may expose `real(kind), intent(out) :: array(:,:,:)` using
three declared extent names and caller-owned C double storage. Native structured
annotations accept `size(array,3)` for an output-defined batch extent. Python
allocates native-endian float64 Fortran-contiguous tensors, or accepts exact
shape/dtype/alignment/writeability through the existing out mapping. The first
axis varies fastest in C storage; a matrix batch `(p,p,n_samples)` stores each
matrix contiguously. Empty extents use valid nonnull scratch buffers.

Rank-three inputs, integer tensors, real intent(inout) tensors, higher ranks and
rank-three function results remain unsupported. Optional real inputs remain
rank one/two. Shape annotations still validate the referenced axis against the
actual rank and configured extent. See the independent rank-three fixture and
[Wishart sampling contract](../../bindings/python/README.md#wishart-and-inverse-wishart-batch-sampling).

## Stateful functions and optional arrays

Array-input functions can combine serialized state codecs with optional-array
presence flags. They return their native scalar/array result through caller-owned
C storage and append call-level status, as other array functions do. State is
decoded before the function call and exported afterward, including when the
native result is NaN. Status 0 means a valid call, not numerical success or
convergence. Extent/flag errors and codec errors preserve result/state.

Python puts state parameters first so required state can precede native optional
inputs; C keeps native argument order with array flags and result/status appended.
Scalar functions return float, with no out mapping. A nonzero C call status raises
RuntimeError; Python storage/ownership validation raises before C. State objects
have no default and must belong to the same loaded library instance. Coverage
and shape tables record native optional state promoted to required C state.
The native optional_policy field remains required: C array pointers are always
required, even when a configured flag makes the corresponding native array absent.

Long Fortran adapter formal lists wrap below the free-form limit. Other existing
short signatures are unchanged. The independent stateful-function fixture tests
optional presence, scalar results, state export on NaN, guards, promotion review
and compiled C/Python behavior. See the
[multivariate CDF contract](../../bindings/python/README.md#multivariate-box-cdf-and-normal-log-cdf).


## Module scope in generated Fortran

Generated adapters use one module-level `implicit none`; contained wrappers
inherit that rule. Intrinsic kinds and shared RNG codec imports appear before
`contains`, alongside native-procedure imports. Codec imports are grouped by
module and reused for identical bindings. Distinct codec bindings receive
separate private aliases, including when different procedures use the same
state argument name; aliases also avoid public wrapper names.

This changes source organization without changing C symbols, prototypes,
argument order, numerical formulas, or state transport. Generator regression
fixtures compile shared and distinct codecs with checked and optimized builds.

Validation commands from the repository root (CMD):

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

## Comments on explicit-shape arguments

Generated declarations omit a leading `shape (...);` or `(...);` comment when
it matches the explicit extents already shown in the declaration. Descriptions,
value constraints, and native assumed-shape annotations are retained. RNG buffer
comments describe the caller-owned components without repeating the declared
component count.

## Rank-three integer outputs

Explicit-array adapters now accept default-integer `intent(out)` arrays with
rank three, alongside real rank-three outputs. C uses `int *`; Python allocates
or validates Fortran-contiguous `np.intc` output cubes. Dimension and state
policies are unchanged. Rank-three input/inout arrays and higher ranks remain
unsupported. The independent rank-three fixture compiles real and integer C and
Python clients; the random_table manifest exercises this pattern in SciFort.

## Complex explicit-array transport

Complex arguments of the configured native kind are supported for rank-one/two
intent(in/out) arrays and rank-three intent(out) arrays. Native complex scalars,
intent(inout), higher ranks and complex function results are rejected. Existing
real/integer interfaces keep their transport and error behavior.

The C ABI uses `const double *` inputs and `double *` outputs with interleaved
(real, imaginary) pairs. Generated real(c_double) bind(c) dummies prepend an axis
of length 2 to the logical shape. Local complex(cg_kind) allocatable arrays are
allocated only after extent, input-policy and state validation. Inputs are
converted with cmplx; output real/imaginary components are exported only when the
native status is zero. State codecs export state even after native errors.
This avoids copying unspecified native outputs on failure, at the cost of
workspace and copies. Caller complex buffers remain unchanged on every error.

Coverage records label complex arguments `interleaved_f64`; generated shape
documentation shows physical pair-axis shapes. The Python runtime allocates and
validates native-endian, aligned, Fortran-contiguous complex128 buffers. Complex
inputs accept numeric real/complex arrays; all input/output overlap remains
forbidden. No language complex object layout or Fortran descriptors cross the
C ABI. The independent complex fixture compiles and runs C11, C++17 and Python
clients with bounds-checked and optimized Fortran builds, including deliberately
partial/undefined native output on error.

## Reviewed scalar-only adapters with paired outputs

A `subroutines` or `array_functions` entry may set `"scalar_only": true` with
`"dimensions": []` and `"shapes": {}`. This explicitly selects the existing
status/output-pointer adapter convention for procedures without native arrays.
Without that opt-in, an empty dimension policy remains an error. Scalar-only
policies reject native arrays or declared dimensions. Serialized RNG state is a
scalar native object with its existing separately reviewed state codec.

Scalar subroutines may return multiple real/integer outputs. Status-free natives
still require `synthetic_status`; optional real/integer inputs still require
`required_optional_scalars`. Function policies retain `result_argument` and
`optional_policy`. Numerical NaN/infinity is retained on call-level status 0.
Coverage records distinguish `scalar_subroutine` and `scalar_function`; the
aggregate report's historical `array_adapters` category includes these reviewed
output/status adapters as well as adapters with array arguments.

The normal-inverse-gamma manifest demonstrates scalar joint densities, paired
moments, paired single sampling and matching batch output vectors. Python uses
one group with scalar density calls and `ArrayResult` for paired outputs.
See [the native shape contracts](../../docs/normal_inverse_gamma_shapes.md).
