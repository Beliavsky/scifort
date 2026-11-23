# C ABI generator prototype

This is an original MIT-licensed, project-independent generator under development
inside SciFort. All project names, source paths, symbol names, kind names, and
argument policies come from a JSON manifest. Its independent arithmetic fixture
does not depend on SciFort. The intended extraction boundary is this directory;
SciFort configuration and integration artifacts remain in the parent project.

Version 0.1 supports private-by-default modules with explicit public exports,
directly defined functions with `result(...)`, explicit scalar `real(kind)`
argument/result declarations, and `intent(in)` inputs. Elemental functions can
also receive a bulk wrapper with one varying input and uniform scalar parameters.
The native real kind and its module are configured; conversion to/from
`c_double` is explicit. SciFort's initial integration uses binary64 throughout.

The generator uses [fparser 0.2.5](https://pypi.org/project/fparser/0.2.5/)
to parse Fortran syntax. It does not parse declarations using regular expressions
or copy numerical implementations. Python 3.10+ and the pinned parser are needed
only to regenerate/check artifacts. Ordinary FPM builds use checked-in wrappers
and do not need Python or fparser.

## Regeneration

From the SciFort repository root, in CMD:

```text
python -m pip install -r tools/cabi_generator/requirements.txt
python tools/cabi_generator/generate.py bindings/cabi_config/normal.json
python tools/cabi_generator/generate.py bindings/cabi_config/normal.json --check
```

The configuration path is relative to the current directory. `--root PATH`
sets the base for all source and output paths in the manifest. The command
produces a Fortran module, a C/C++ header, and a JSON coverage record. Output
contains no timestamps or absolute paths; the source SHA-256 uses UTF-8 with
line endings normalized to LF, so Windows and Unix checkouts agree. `--check`
exits nonzero if any artifact
is missing or differs, and never writes files. Invalid manifests and unsupported
selected interfaces are rejected before outputs are written.

## Configuration and contract

See the [normal manifest](../../bindings/cabi_config/normal.json) for a complete
schema-version-1 example. The top-level fields are:

| Field | Meaning |
| --- | --- |
| `source`, `module` | Fortran source and native module to inspect |
| `kind_module`, `kind_name` | Module exporting the real kind used by inputs/results |
| `wrapper_module`, `header_guard` | Generated Fortran module and C include guard |
| `outputs` | Paths for `fortran`, `header`, and `coverage` artifacts |
| `argument_docs` | One ASCII description per native input argument |
| `functions` | Public functions to wrap and their C interface policies |
| `excluded` | Every remaining explicit public export, with a reason |

Each function policy requires `symbol`. `vector_symbol` and `array_argument`
request a bulk entry point for an elemental function. `finite` and `positive`
list uniform parameters checked before any output is changed. These constraints
are explicit decisions, not guesses based on argument names. Native optional
arguments require `optional_policy: "required"`: callers must supply them at
the C boundary, and no defaults are inferred from procedure implementations.
Scalar wrappers delegate numerical validation to the native procedure.

Bulk calls use `size_t`, raw input/output pointers, uniform doubles, and a
caller-owned integer status. Status 0 means success; 1 means an invalid uniform
parameter, null nonempty buffer, or length exceeding the compiler's default
Fortran integer maximum. Output remains unchanged on those errors. Zero-length
calls allow null buffers, but still validate uniform parameters. The status
pointer must always be valid. Nonempty arrays must be suitably aligned, contain
at least `n` doubles, and must not overlap. Invalid individual input values retain
the native function's NaN/infinity/endpoint behavior and do not become call-level
status failures.

Every explicit public export must be generated or explicitly excluded. Public
generics, reexports, derived types, arrays, callbacks, subroutines, and other
unsupported interfaces cannot silently become scalar wrappers. This is not yet
a whole-project scanner or a full Fortran semantic resolver. Default-public
modules are rejected because complete reexport discovery needs more analysis.
There is no Python binding generation in this first version.

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

Initial Windows validation on 2026-11-23 passed all 12 focused tests (11 reusable
tests and one SciFort integration test with multiple methods/input cases), using
Python 3.13.3, fparser 0.2.5, and GNU Fortran/GCC/G++ 17.0.0 20260510
(experimental). The commands above and the generator's `--check` command were
run with `build/cabi-generator-venv/Scripts/python.exe`. Other operating systems
and compilers have not yet been validated for this generator.

These focused checks cover this additive normal-family interface change; the
numerical implementations, existing C symbols, and global build configuration
are unchanged. Run the full project suite before release or merge as required
by the repository rules.

## Next interface patterns

1. Additional scalar/elemental distribution families and bulk parameters.
2. Explicit-rank array procedures and caller-owned output sizing.
3. Allocatable results using size-query/fill interfaces.
4. Derived types using opaque handles with defined lifetimes.
5. Callbacks, strings, and generic resolution with explicit policies.

Keep these rules reusable; project-specific mathematical constraints and ABI
compatibility decisions belong in manifests and integration tests. fparser is
an external build tool, not vendored code or a linked numerical dependency;
see [provenance](../../CODE_PROVENANCE.md) and its retained
[license text](../../LICENSES/fparser-BSD-3-Clause.txt).
