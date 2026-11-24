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
Rank-one/rank-two assumed-shape array subroutines can now receive explicit-shape
`bind(c)` adapters, with a required native integer status output. Their native
implementations are unchanged; no C implementation or pointer conversion is needed.

The generator uses [fparser 0.2.5](https://pypi.org/project/fparser/0.2.5/)
to parse Fortran syntax. It does not parse declarations using regular expressions
or copy numerical implementations. Python 3.10+ and the pinned parser are needed
only to regenerate/check artifacts. Ordinary FPM builds use checked-in wrappers
and do not need Python or fparser.

## Regeneration

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
entity, for example `real(dp), intent(in) :: a(:,:)`. Optional, allocatable,
pointer, callback, and additional integer arguments are not supported in this
pattern. The required native status argument is default-integer `intent(out)`.
Compile-time checks require the native real kind to equal `c_double` and the
native default integer kind to equal `c_int`. Arrays are passed directly to the
native procedure, without conversion buffers or wrapper allocation.

All arrays use contiguous column-major storage. Callers must provide valid
non-NULL buffers, including for permitted empty arrays, and adequate storage for
every declared extent; these adapters cannot check buffer addresses or capacities.
Output and mutable arguments use `intent(inout)` at the adapter boundary so an
invalid-extent return preserves caller data. After validation, native status
codes and output-on-failure semantics are preserved. Numerical failures can
therefore leave native partial results. No leading dimensions, arbitrary strides,
or pointer-style subroutine adapters are implemented yet.

Every explicit public export must be generated or explicitly excluded. Public
generics, reexports, derived types, callbacks, array-valued functions, and other
unsupported interfaces cannot silently become scalar wrappers. This is not yet
a whole-project scanner or a full Fortran semantic resolver. Default-public
modules are rejected because complete reexport discovery needs more analysis.
Python scalar/elemental binding generation is described below; other interface
patterns are not yet exposed through Python.

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
