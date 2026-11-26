#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate explicit C wrappers for a deliberately small Fortran API subset."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

from fparser.common.readfortran import FortranFileReader
from fparser.two import Fortran2003 as F
from fparser.two.parser import ParserFactory
from fparser.two.utils import FortranSyntaxError, walk

import explicit_arrays

VERSION = "0.3.0"


def identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r"[a-z][a-z0-9_]{0,62}", value):
        raise ValueError(f"expected a lowercase Fortran identifier, got {value!r}")
    return value


def public_api(path, module_name, *, tree=None, specifics=None):
    """Read an AST; never infer a procedure signature from its implementation."""
    if tree is None:
        tree = ParserFactory().create(std="f2008")(FortranFileReader(str(path)))
    modules = [m for m in walk(tree, F.Module)
               if str(m.content[0].items[1]).lower() == module_name]
    if len(modules) != 1:
        raise ValueError(f"expected exactly one module {module_name} in {path}")
    module = modules[0]
    spec = next((n for n in module.content if isinstance(n, F.Specification_Part)), None)
    body = next((n for n in module.content if isinstance(n, F.Module_Subprogram_Part)), None)
    procedures = {}
    if body:
        for node in body.content:
            if isinstance(node, (F.Function_Subprogram, F.Subroutine_Subprogram)):
                procedures[str(node.content[0].items[1]).lower()] = node
    default = "public"
    access = {}
    if spec:
        for node in spec.content:
            if isinstance(node, F.Access_Stmt):
                mode, names = node.items
                if names is None:
                    default = str(mode).lower()
                else:
                    access.update((str(n).lower(), str(mode).lower()) for n in names.items)
            elif isinstance(node, F.Type_Declaration_Stmt) and node.items[1]:
                attrs = {str(a).lower() for a in node.items[1].items}
                for mode in ("public", "private"):
                    if mode in attrs:
                        access.update((str(e.items[0]).lower(), mode)
                                      for e in node.items[2].items)
            elif isinstance(node, F.Derived_Type_Def):
                attributes, name, _ = node.content[0].items
                attrs = {str(a).lower() for a in attributes.items} if attributes else set()
                for mode in ("public", "private"):
                    if mode in attrs:
                        access[str(name).lower()] = mode
    if default != "private":
        raise ValueError("first release requires a private-by-default module with explicit public exports")
    names = {name for name in procedures if access.get(name, default) == "public"}
    names.update(name for name, mode in access.items() if mode == "public")
    api = {name: procedures.get(name) for name in sorted(names)}
    specifics = {} if specifics is None else specifics
    if not isinstance(specifics, dict):
        raise ValueError('generic_specifics must map public generic names to local procedures')
    interfaces = {}
    for block in (spec.content if spec else ()):
        if not isinstance(block, F.Interface_Block):
            continue
        generic = block.content[0].items[0]
        if not isinstance(generic, F.Name):
            continue
        members = set()
        for statement in block.content[1:-1]:
            if not isinstance(statement, F.Procedure_Stmt) or str(statement.items[1]).lower() != 'module':
                continue
            members.update(str(n).lower() for n in statement.items[0].items)
        interfaces[str(generic).lower()] = members
    for generic, specific in specifics.items():
        identifier(generic)
        identifier(specific)
        if (generic not in api or generic not in interfaces or specific not in interfaces[generic]
                or specific not in procedures):
            raise ValueError('generic selection requires a public named interface and a declared local module procedure')
        api[generic] = procedures[specific]
    return api


def signature(node, kind):
    if not isinstance(node, F.Function_Subprogram):
        raise ValueError("only directly defined functions are supported")
    statement = node.content[0]
    prefix, name, args, suffix = statement.items
    if suffix is None or suffix.items[0] is None or suffix.items[1] is not None:
        raise ValueError("a function needs explicit result(...) and must not already bind(c)")
    if any(str(p).lower() == "impure" for p in (prefix.items if prefix else ())):
        raise ValueError("impure functions are not supported")
    args = [identifier(str(a).lower()) for a in (args.items if args else ())]
    result = str(suffix.items[0]).lower()
    spec = next((n for n in node.content if isinstance(n, F.Specification_Part)), None)
    declarations = {}
    if spec:
        for decl in spec.content:
            if isinstance(decl, F.Type_Declaration_Stmt):
                dtype, attributes, entities = decl.items
                attrs = {str(a).lower() for a in attributes.items} if attributes else set()
                for entity in entities.items:
                    declarations[str(entity.items[0]).lower()] = (dtype, attrs, entity)
    for arg in [*args, result]:
        if arg not in declarations:
            raise ValueError(f"{arg}: explicit local type declaration required")
        dtype, attrs, entity = declarations[arg]
        if str(dtype).lower().replace(" ", "") not in (f"real({kind})", f"real(kind={kind})"):
            raise ValueError(f"{arg}: only real({kind}) is supported")
        allowed = {"intent(in)", "optional"} if arg != result else set()
        if attrs - allowed or any(item is not None for item in entity.items[1:]):
            raise ValueError(f"{arg}: only scalar intent(in) inputs and scalar results are supported")
        if arg != result and "intent(in)" not in attrs:
            raise ValueError(f"{arg}: intent(in) is required")
    optional = [arg for arg in args if "optional" in declarations[arg][1]]
    elemental = prefix is not None and any(str(p).lower() == "elemental" for p in prefix.items)
    return args, optional, elemental


def generate(config, root):
    """Return deterministic artifacts; validation completes before any files are written."""
    required = {"schema_version", "source", "module", "kind_module", "kind_name",
                "wrapper_module", "header_guard", "functions", "excluded", "outputs",
                "argument_docs"}
    if (not isinstance(config, dict) or set(config) - {"subroutines", "array_functions", "generic_specifics"} != required
            or type(config["schema_version"]) is not int or config["schema_version"] not in (1, 2)):
        raise ValueError(f"schema versions 1/2 require these fields (plus optional subroutines): {sorted(required)}")
    default_array_style = "pointer" if config["schema_version"] == 1 else "explicit"
    for key in ("functions", "excluded", "outputs", "argument_docs"):
        if not isinstance(config[key], dict):
            raise ValueError(f"{key} must be an object")
    if not isinstance(config["source"], str) or not all(
            isinstance(value, str) for value in config["outputs"].values()):
        raise ValueError("source and output paths must be strings")
    for key in ("module", "kind_module", "kind_name", "wrapper_module"):
        identifier(config[key])
    if config["wrapper_module"] in (config["module"], config["kind_module"]):
        raise ValueError("wrapper_module must differ from the native and kind modules")
    guard = config["header_guard"]
    if not isinstance(guard, str) or not re.fullmatch(r"[A-Z][A-Z0-9_]*", guard):
        raise ValueError("header_guard must be an uppercase C identifier")
    if set(config["outputs"]) - {"shape_documentation"} != {"fortran", "header", "coverage"}:
        raise ValueError("outputs must specify fortran, header, and coverage")
    source = (root / config["source"]).resolve()
    outputs = {key: (root / value).resolve() for key, value in config["outputs"].items()}
    if 'shape_documentation' in outputs and outputs['shape_documentation'].suffix.lower() != '.md':
        raise ValueError('shape_documentation must name a Markdown (.md) output')
    if len(set(outputs.values())) != len(outputs) or source in outputs.values():
        raise ValueError("output paths must be distinct and must not overwrite the source")
    for key in ('subroutines', 'array_functions'):
        if not isinstance(config.get(key, {}), dict):
            raise ValueError(f'{key} must be an object')
    if not isinstance(config.get("generic_specifics", {}), dict):
        raise ValueError("generic_specifics must be an object")
    api = public_api(source, config["module"], specifics=config.get("generic_specifics"))
    selected, excluded = config["functions"], config["excluded"]
    subroutines = config.get("subroutines", {})
    array_functions = config.get('array_functions', {})
    if (set(selected) & (set(subroutines) | set(array_functions))
            or set(subroutines) & set(array_functions)):
        raise ValueError("a procedure cannot be both a function and subroutine policy")
    covered = set(selected) | set(subroutines) | set(array_functions)
    if set(config.get("generic_specifics", {})) - covered:
        raise ValueError("generic selections must name generated procedures")
    if covered & set(excluded):
        raise ValueError("a procedure cannot be both generated and excluded")
    if covered | set(excluded) != set(api):
        raise ValueError(f"public API coverage mismatch: unaccounted={sorted(set(api)-covered-set(excluded))}, "
                         f"unknown={sorted((covered|set(excluded))-set(api))}")
    if not all(isinstance(reason, str) and reason.strip() for reason in excluded.values()):
        raise ValueError("every exclusion needs a nonempty reason")
    entries, symbols = [], set()
    for name, policy in sorted(selected.items()):
        identifier(name)
        allowed = {"symbol", "vector_symbol", "array_argument", "optional_policy", "finite", "positive", "array_style", "existing_scalar_symbol", "elementwise_symbol"}
        if not isinstance(policy, dict) or set(policy) - allowed or not (policy.get("symbol") or policy.get("vector_symbol") or policy.get("elementwise_symbol")):
            raise ValueError(f"{name}: invalid or missing policy fields")
        for key in ("symbol", "vector_symbol", "existing_scalar_symbol", "elementwise_symbol"):
            if key in policy:
                identifier(policy[key])
        existing = policy.get("existing_scalar_symbol")
        if existing is not None:
            identifier(existing)
            if policy.get("symbol"):
                raise ValueError(f"{name}: existing_scalar_symbol and symbol are mutually exclusive")
        args, optional, elemental = signature(api[name], config["kind_name"])
        if optional and policy.get("optional_policy") != "required":
            raise ValueError(f"{name}: explicitly set optional_policy=required to promote optional inputs")
        vector = policy.get("vector_symbol")
        if policy.get("elementwise_symbol") and not elemental:
            raise ValueError(f"{name}: elementwise wrappers require an elemental function")
        if policy.get("array_style", default_array_style) not in ("explicit", "pointer"):
            raise ValueError(f"{name}: array_style must be explicit or pointer")
        if "array_style" in policy and not vector:
            raise ValueError(f"{name}: array_style needs a vector_symbol")
        if vector and (not elemental or policy.get("array_argument") not in args):
            raise ValueError(f"{name}: bulk wrappers require an elemental function and an array_argument")
        if not vector and ("array_argument" in policy or policy.get("finite") or policy.get("positive")):
            raise ValueError(f"{name}: bulk policy supplied without a vector_symbol")
        for check in ("finite", "positive"):
            values = policy.get(check, [])
            if not isinstance(values, list) or not all(a in args and a != policy.get("array_argument") for a in values):
                raise ValueError(f"{name}: {check} must list uniform parameter names")
        for arg in args:
            doc = config["argument_docs"].get(arg)
            if not isinstance(doc, str) or not doc or "\n" in doc or not doc.isascii():
                raise ValueError(f"{name}: argument_docs needs one ASCII line for {arg}")
        for symbol in (policy.get("symbol"), vector, policy.get("elementwise_symbol")):
            if symbol:
                identifier(symbol)
                if symbol in symbols:
                    raise ValueError(f"duplicate C symbol: {symbol}")
                symbols.add(symbol)
        entries.append((name, args, optional, policy))
    array_entries = explicit_arrays.prepare(config, api, identifier, source.read_text(encoding="utf-8").splitlines())
    for entry in array_entries:
        symbol = entry["policy"]["symbol"]
        if symbol in symbols:
            raise ValueError(f"duplicate C symbol: {symbol}")
        symbols.add(symbol)
    existing_symbols = {p["existing_scalar_symbol"] for _, _, _, p in entries
                        if p.get("existing_scalar_symbol")}
    if existing_symbols & symbols:
        raise ValueError("an existing scalar symbol cannot also be generated")
    module = config["wrapper_module"]
    lines = ["! SPDX-License-Identifier: MIT", "! Generated by cabi_generator; do not edit.",
             f"module {module}",
             "    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_ptr, &",
             "        c_associated, c_f_pointer",
             "    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite",
             f"    use {config['kind_module']}, only : cg_kind => {config['kind_name']}"]
    if not any(policy.get("vector_symbol") and policy.get("array_style", default_array_style) == "pointer"
               for _, _, _, policy in entries):
        lines[3:5] = ["    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t"]
    for index, (name, _, _, _) in enumerate(entries):
        lines.append(f"    use {config['module']}, only : cg_native_{index} => {name}")
    for index, entry in enumerate(array_entries, start=len(entries)):
        lines.append(f"    use {config['module']}, only : cg_native_{index} => {entry['name']}")
    lines += ["    implicit none", "    private"]
    for symbol in sorted(symbols):
        lines.append(f"    public :: {symbol}")
    lines += ["contains", ""]
    header = ["/* SPDX-License-Identifier: MIT */", "/* Generated by cabi_generator; do not edit. */",
              f"#ifndef {guard}", f"#define {guard}", "#include <stddef.h>",
              "#ifdef __cplusplus", 'extern "C" {', "#endif", "",
              "/* All native optional arguments are required at the C boundary.",
              " * Elemental bulk: uniform scalar parameters, one input/output array.",
              " * Status 0: success; 1: invalid uniform parameter, null nonempty buffer,",
              " * or length exceeding the default Fortran integer maximum.",
              " * On call-level errors output is unchanged. Empty calls still validate parameters.",
              " * Buffer rules are stated per entry point. Nonempty buffers must not overlap.",
              " * Invalid individual values follow the native function's behavior.",
              " * The status pointer must always be valid. */", ""]
    coverage = []
    for index, (name, args, optional, policy) in enumerate(entries):
        aliases = {arg: f"cg_arg_{i}" for i, arg in enumerate(args)}
        symbol = policy.get("symbol")
        arglist = list(aliases.values())
        if symbol:
            lines += [f"    function {symbol}( &", "        " + ", ".join(arglist) + ") &",
                      f'        result(cg_result) bind(c, name="{symbol}")',
                      "        !! Forward one scalar evaluation to the native Fortran function.", "        implicit none"]
            for arg in args:
                lines.append(f"        real(c_double), value, intent(in) :: {aliases[arg]} "
                             f"!! {config['argument_docs'][arg]}; required")
            lines += ["        real(c_double) :: cg_result", ""]
            lines += call_lines(index, args, aliases, "cg_result")
            lines += [f"    end function {symbol}", ""]
            header.append(f"double {symbol}(" + (", ".join(f"double arg_{a}" for a in args) or "void") + ");")
        vector = policy.get("vector_symbol")
        if vector:
            array = policy["array_argument"]
            pointer = policy.get("array_style", default_array_style) == "pointer"
            vecargs = ["cg_count", *arglist, "cg_output", "cg_status"]
            lines += [f"    subroutine {vector}( &", "        " + ", ".join(vecargs) + ") &",
                      f'        bind(c, name="{vector}")',
                      "        !! Evaluate independent inputs using uniform scalar parameters.",
                      "        implicit none",
                      "        integer(c_size_t), value, intent(in) :: cg_count !! array length, >= 0"]
            for arg in args:
                if arg == array and pointer:
                    declaration, suffix = "type(c_ptr), value, intent(in)", ""
                elif arg == array:
                    declaration, suffix = "real(c_double), intent(in)", "(cg_count)"
                else:
                    declaration, suffix = "real(c_double), value, intent(in)", ""
                lines.append(f"        {declaration} :: {aliases[arg]}{suffix} "
                             f"!! {config['argument_docs'][arg]}; required")
            if pointer:
                lines.append("        type(c_ptr), value, intent(in) :: cg_output !! caller-owned output buffer")
                lines.append("        real(c_double), pointer :: cg_x(:), cg_y(:)")
            else:
                lines.append("        real(c_double), intent(inout) :: cg_output(cg_count) !! output buffer")
            lines += ["        integer(c_int), intent(out) :: cg_status !! 0 success; 1 invalid call",
                      "        integer :: cg_i", "", "        cg_status=1_c_int",
                      "        if (cg_count < 0_c_size_t .or. cg_count > int(huge(0),c_size_t)) return"]
            for arg in policy.get("finite", []):
                lines.append(f"        if (.not. ieee_is_finite({aliases[arg]})) return")
            for arg in policy.get("positive", []):
                lines.append(f"        if (.not. ({aliases[arg]} > 0.0_c_double)) return")
            if pointer:
                lines += ["        if (cg_count == 0_c_size_t) then", "            cg_status=0_c_int",
                          "            return", "        end if",
                          f"        if (.not. c_associated({aliases[array]})) return",
                          "        if (.not. c_associated(cg_output)) return",
                          f"        call c_f_pointer({aliases[array]},cg_x,[int(cg_count)])",
                          "        call c_f_pointer(cg_output,cg_y,[int(cg_count)])"]
            lines += ["        do cg_i=1,int(cg_count)"]
            element = "cg_x(cg_i)" if pointer else f"{aliases[array]}(cg_i)"
            output_element = "cg_y(cg_i)" if pointer else "cg_output(cg_i)"
            actuals = dict(aliases, **{array: element})
            lines += call_lines(index, args, actuals, output_element, indent="            ")
            lines += ["        end do", "        cg_status=0_c_int", f"    end subroutine {vector}", ""]
            cargs = ["size_t n", *[("const double *" if a == array else "double ") + f"arg_{a}" for a in args],
                     "double *output", "int *status"]
            header.append("/* Pointer bulk: NULL buffers are allowed only for n=0. */" if pointer else
                          "/* Explicit-shape bulk: valid non-NULL buffers are required even for n=0. */")
            header.append(f"void {vector}(" + ", ".join(cargs) + ");")
        elementwise = policy.get("elementwise_symbol")
        if elementwise:
            lines += [f"    subroutine {elementwise}( &",
                      "        " + ", ".join(["cg_count", *arglist, "cg_output", "cg_status"]) + ") &",
                      f'        bind(c, name="{elementwise}")',
                      "        !! Evaluate corresponding elements of equally sized input arrays.",
                      "        implicit none",
                      "        integer(c_size_t), value, intent(in) :: cg_count !! shared length, >= 0"]
            for arg in args:
                lines.append(f"        real(c_double), intent(in) :: {aliases[arg]}(cg_count) "
                             f"!! {config['argument_docs'][arg]}; per element")
            lines += ["        real(c_double), intent(inout) :: cg_output(cg_count) !! output buffer",
                      "        integer(c_int), intent(out) :: cg_status !! 0 success; 1 invalid length",
                      "        integer :: cg_i", "", "        cg_status=1_c_int",
                      "        if (cg_count < 0_c_size_t .or. cg_count > int(huge(0),c_size_t)) return",
                      "        do cg_i=1,int(cg_count)"]
            lines += call_lines(index, args, {a: f"{aliases[a]}(cg_i)" for a in args},
                                "cg_output(cg_i)", indent="            ")
            lines += ["        end do", "        cg_status=0_c_int",
                      f"    end subroutine {elementwise}", ""]
            header += ["/* Elementwise: all inputs are arrays of n doubles, in native argument order.",
                       " * Valid aligned non-NULL buffers required even for n=0; output must not",
                       " * overlap any input. Inputs may overlap each other. No uniform validation:",
                       " * invalid parameter elements retain native results (usually NaN), status 0.",
                       " * Status 1: rejected length, output unchanged. */",
                       f"void {elementwise}(" + ", ".join(["size_t n",
                           *[f"const double *arg_{a}" for a in args],
                           "double *output", "int *status"]) + ");"]
        header.append("")
        coverage.append({"procedure": name, "state": "generated", "arguments": args,
                         "optional_promoted_to_required": optional,
                         "symbols": ([symbol] if symbol else []) + ([vector] if vector else [])
                                    + ([elementwise] if elementwise else [])})
        if policy.get("existing_scalar_symbol"):
            coverage[-1]["existing_scalar_symbol"] = policy["existing_scalar_symbol"]
    for index, entry in enumerate(array_entries, start=len(entries)):
        body, prototypes, record = explicit_arrays.emit(entry, index, config)
        lines += body
        header += prototypes
        coverage.append(record)
    for record in coverage:
        if record["procedure"] in config.get("generic_specifics", {}):
            record["native_specific"] = config["generic_specifics"][record["procedure"]]
    for name, reason in sorted(excluded.items()):
        coverage.append({"procedure": name, "state": "excluded", "reason": reason})
    lines.append(f"end module {module}")
    if any(len(line) > 132 for line in lines):
        raise ValueError("generated Fortran exceeds 132 columns; shorten configured names/docs")
    header += ["#ifdef __cplusplus", "}", "#endif", f"#endif /* {guard} */"]
    report = {"generator_version": VERSION, "source": config["source"],
              "source_sha256": hashlib.sha256(source.read_text(encoding="utf-8").encode()).hexdigest(),
              "public_api_count": len(api), "generated_symbol_count": len(symbols),
              "procedures": sorted(coverage, key=lambda row: row["procedure"])}
    artifacts = {outputs["fortran"]: "\n".join(lines) + "\n",
            outputs["header"]: "\n".join(header) + "\n",
            outputs["coverage"]: json.dumps(report, indent=2, sort_keys=True) + "\n"}
    if 'shape_documentation' in outputs:
        if not any(entry.get('shape_contract') for entry in array_entries):
            raise ValueError('shape_documentation needs at least one opted-in shape contract')
        artifacts[outputs['shape_documentation']] = explicit_arrays.shape_contracts.markdown(config, array_entries)
    return artifacts


def call_lines(index, args, actuals, output, indent="        "):
    lines = [f"{indent}{output}=real(cg_native_{index}( &"]
    for i, arg in enumerate(args):
        ending = ", &" if i < len(args) - 1 else "),c_double)"
        lines.append(f"{indent}    {arg}=real({actuals[arg]},cg_kind){ending}")
    if not args:
        lines = [f"{indent}{output}=real(cg_native_{index}(),c_double)"]
    return lines


def collect(config_path, root, batch=False):
    """Validate a whole batch, including collisions, before writing any output."""
    value = json.loads(config_path.read_text(encoding="utf-8"))
    if batch:
        if not isinstance(value, list) or not value or not all(isinstance(p, str) for p in value):
            raise ValueError("batch must be a nonempty JSON list of manifest paths")
        paths = [(config_path.parent / p).resolve() for p in value]
    else:
        paths = [config_path.resolve()]
    if len(set(paths)) != len(paths):
        raise ValueError("duplicate manifest in batch")
    artifacts, sources, symbols, existing_symbols = {}, set(), set(), set()
    protected = {config_path.resolve(), *paths}
    for path in paths:
        config = json.loads(path.read_text(encoding="utf-8"))
        generated = generate(config, root)
        sources.add((root / config["source"]).resolve())
        if artifacts.keys() & generated.keys():
            raise ValueError("duplicate output path in batch")
        report = json.loads(generated[(root / config["outputs"]["coverage"]).resolve()])
        current = {s for row in report["procedures"] for s in row.get("symbols", [])}
        existing = {row["existing_scalar_symbol"] for row in report["procedures"]
                    if "existing_scalar_symbol" in row}
        if symbols & current:
            raise ValueError("duplicate generated C symbol in batch")
        symbols.update(current)
        existing_symbols.update(existing)
        artifacts.update(generated)
    if artifacts.keys() & (sources | protected):
        raise ValueError("output paths must not overwrite sources or configurations")
    if symbols & existing_symbols:
        raise ValueError("an existing scalar symbol cannot also be generated in batch")
    return artifacts


def print_documentation_reviews(artifacts):
    """Surface advisory native shape-documentation issues without rewriting source."""
    for path, content in artifacts.items():
        if path.suffix != '.json':
            continue
        report = json.loads(content)
        if not isinstance(report, dict):
            continue
        for row in report.get('procedures', []):
            for issue in row.get('native_shape_documentation_review', []):
                shape = '(' + ','.join(issue['required_shape']) + ')'
                print(f"cabi_generator documentation review: {report['source']}:"
                      f"{row['procedure']}:{issue['argument']}: expected shape {shape}; "
                      f"{issue['reason']}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--batch", action="store_true",
                        help="config is a JSON list of manifest paths relative to that list")
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="base for configured source/output paths")
    parser.add_argument("--check", action="store_true", help="check generated artifacts without writing")
    args = parser.parse_args()
    try:
        artifacts = collect(args.config, args.root.resolve(), args.batch)
        print_documentation_reviews(artifacts)
        if args.check:
            stale = [str(path) for path, value in artifacts.items()
                     if not path.exists() or path.read_text(encoding="utf-8") != value]
            if stale:
                raise ValueError("generated files are missing or stale: " + ", ".join(stale))
        else:
            for path, value in artifacts.items():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(value, encoding="utf-8", newline="\n")
        print(f"{'Checked' if args.check else 'Generated'} {len(artifacts)} artifacts")
        return 0
    except (ValueError, OSError, FortranSyntaxError) as error:
        print(f"cabi_generator: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
