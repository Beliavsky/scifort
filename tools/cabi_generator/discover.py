#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Inventory public Fortran APIs and C bindings; emit deliberately non-executable drafts."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

from fparser.common.readfortran import FortranFileReader
from fparser.two import Fortran2003 as F
from fparser.two.parser import ParserFactory
from fparser.two.utils import walk
import generate
import generate_python


def declarations(node):
    """Describe local declarations without interpreting mathematical constraints."""
    spec = next((n for n in node.content if isinstance(n, F.Specification_Part)), None)
    result = {}
    if spec:
        for decl in spec.content:
            if not isinstance(decl, F.Type_Declaration_Stmt):
                continue
            dtype, attrs, entities = decl.items
            attributes = [str(a).lower() for a in attrs.items] if attrs else []
            dimension = next((a.items[1] for a in (attrs.items if attrs else ())
                              if isinstance(a, F.Dimension_Attr_Spec)), None)
            for entity in entities.items:
                shape = entity.items[1] or dimension
                rank = (1 + (len(shape.items[0].items) if shape.items[0] else 0)
                        if isinstance(shape, F.Assumed_Size_Spec) else len(shape.items) if shape else 0)
                result[str(entity.items[0]).lower()] = dict(
                    type=str(dtype).lower(), rank=rank,
                    shape=str(shape).lower() if shape else None, attributes=attributes,
                    optional='optional' in attributes)
    return result


def arguments(node):
    names = node.content[0].items[2]
    declared = declarations(node)
    return [dict(name=str(n).lower(), **declared.get(str(n).lower(),
                 dict(type=None, rank=None, shape=None, attributes=[], optional=None)))
            for n in (names.items if names else ())]


def imports(scope):
    """Resolve explicit ONLY imports, including local => original renames."""
    spec = next((n for n in scope.content if isinstance(n, F.Specification_Part)), None)
    result = {}
    if spec:
        for use in spec.content:
            if not isinstance(use, F.Use_Stmt) or str(use.items[3]).strip().lower() != ', only:':
                continue
            module = str(use.items[2]).lower()
            for item in use.items[4].items:
                if isinstance(item, F.Name):
                    local = original = str(item).lower()
                elif isinstance(item, F.Rename) and item.items[0] is None:
                    local, original = (str(n).lower() for n in item.items[1:])
                else:
                    continue
                result[local] = (module, original)
    return result


def binding(node):
    statement = node.content[0]
    if isinstance(node, F.Function_Subprogram):
        suffix = statement.items[3]
        spec = suffix.items[1] if suffix else None
    else:
        spec = statement.items[3]
    if spec is None:
        return None
    # Language_Binding_Spec stores its optional NAME literal as its only item.
    literal = spec.items[0]
    return str(literal.items[0])[1:-1] if literal else str(statement.items[1]).lower()


def abi_symbols(paths, root):
    rows = []
    parser = ParserFactory().create(std='f2008')
    for path in paths:
        tree = parser(FortranFileReader(str(path)))
        for module in walk(tree, F.Module):
            common = imports(module)
            body = next((n for n in module.content if isinstance(n, F.Module_Subprogram_Part)), None)
            for node in (body.content if body else ()):
                if not isinstance(node, (F.Function_Subprogram, F.Subroutine_Subprogram)):
                    continue
                symbol = binding(node)
                if symbol is None:
                    continue
                used = dict(common, **imports(node))
                calls = set()
                for cls in (F.Part_Ref, F.Function_Reference, F.Structure_Constructor, F.Call_Stmt):
                    for call in walk(node, cls):
                        local = str(call.items[0]).lower()
                        if local in used:
                            calls.add(used[local])
                declared_args = arguments(node)
                scalar = isinstance(node, F.Function_Subprogram) and all(
                    a['rank'] == 0 and 'value' in a['attributes'] and a['type'] != 'type(c_ptr)'
                    for a in declared_args)
                rows.append(dict(symbol=symbol, source=relative(path, root),
                                 procedure=str(node.content[0].items[1]).lower(),
                                 interface='scalar' if scalar else 'procedural',
                                 arguments=declared_args,
                                 imported_calls=[dict(module=m, procedure=p) for m, p in sorted(calls)]))
    return sorted(rows, key=lambda r: (r['symbol'], r['source']))


def relative(path, root):
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def files(paths):
    result = set()
    for value in paths:
        path = Path(value)
        if path.is_dir():
            result.update(p.resolve() for p in path.rglob('*.f90'))
        elif path.is_file():
            result.add(path.resolve())
        else:
            raise ValueError(f'source path does not exist: {value}')
    return sorted(result)


def configured(root, batch=None, python=None):
    native = {}
    if batch:
        generate.collect(batch, root, True)  # Validate policies and collision rules first.
        for name in json.loads(batch.read_text(encoding='utf-8')):
            config = json.loads((batch.parent / name).read_text(encoding='utf-8'))
            for name, policy in config['functions'].items():
                row = native.setdefault((config['module'], name), {})
                row['scalar'] = policy.get('symbol') or policy.get('existing_scalar_symbol')
                row['bulk'] = policy.get('vector_symbol')
                if policy.get('elementwise_symbol'):
                    row['elementwise'] = policy['elementwise_symbol']
            for name, policy in config.get('subroutines', {}).items():
                native.setdefault((config['module'], name), {})['array_subroutine'] = policy['symbol']
    if python:
        generate_python.artifacts(python, root)  # Validate signatures and explicit defaults without writing.
        policy = json.loads(python.read_text(encoding='utf-8'))
        for group, options in policy['groups'].items():
            config = json.loads((root / options['manifest']).read_text(encoding='utf-8'))
            selected = dict(config['functions'], **options['existing_functions'])
            for name, entry in selected.items():
                row = native.setdefault((config['module'], name), {})
                row.update(scalar=entry.get('symbol') or entry.get('existing_scalar_symbol'),
                           bulk=entry.get('vector_symbol'),
                           python=f'{group}.{name[len(options["prefix"]):]}')
        for group, settings in policy.get('array_groups', {}).items():
            for manifest in settings['manifests']:
                config = json.loads((root/manifest).read_text(encoding='utf-8'))
                for name in config.get('subroutines', {}):
                    native.setdefault((config['module'], name), {})['python'] = f'{group}.{name}'
    return native


def inventory(paths, root, kind='dp', abi_paths=(), policies=None, max_source_bytes=None, progress=False):
    modules, errors = [], []
    symbols = abi_symbols(abi_paths, root)
    policies = policies or {}
    parser = ParserFactory().create(std='f2008')
    for index, path in enumerate(paths, 1):
        if progress:
            print(f'Inspecting {index}/{len(paths)}: {relative(path, root)}', file=sys.stderr, flush=True)
        if max_source_bytes is not None and path.stat().st_size > max_source_bytes:
            errors.append(dict(source=relative(path, root), reason='exceeds explicit max-source-bytes limit'))
            continue
        try:
            tree = parser(FortranFileReader(str(path)))
            names = [str(m.content[0].items[1]).lower() for m in walk(tree, F.Module)]
        except (generate.FortranSyntaxError, OSError) as error:
            errors.append(dict(source=relative(path, root), reason=str(error)))
            continue
        for module in names:
            record = dict(source=relative(path, root), module=module, procedures=[])
            try:
                api = generate.public_api(path, module, tree=tree)
            except ValueError as error:
                record['error'] = str(error)
                modules.append(record)
                continue
            for name, node in api.items():
                row = dict(name=name, kind='unresolved_public_export', arguments=[],
                           scalar_supported=False, bulk_supported=False)
                if isinstance(node, (F.Function_Subprogram, F.Subroutine_Subprogram)):
                    row.update(kind='function' if isinstance(node, F.Function_Subprogram) else 'subroutine',
                               arguments=arguments(node))
                    try:
                        args, optional, elemental = generate.signature(node, kind)
                        row.update(scalar_supported=True, bulk_supported=elemental, elemental=elemental,
                                   result=dict(name=str(node.content[0].items[3].items[0]).lower(),
                                               **declarations(node)[str(node.content[0].items[3].items[0]).lower()]),
                                   optional_arguments=optional)
                        if not elemental:
                            row['bulk_reason'] = 'bulk pattern requires an elemental scalar function'
                    except ValueError as error:
                        row['unsupported_reason'] = str(error)
                else:
                    row['unsupported_reason'] = 'generic, type, constant or reexport needs explicit resolution'
                if isinstance(node, F.Subroutine_Subprogram):
                    row['unsupported_reason'] = ('scalar/elemental pattern does not apply; array adapters require '
                                                 'reviewed dimensions, shapes and status policy')
                key = (module, name)
                observed = [s for s in symbols if dict(module=module, procedure=name) in s['imported_calls']]
                row['observed_c_symbols'] = [dict(symbol=s['symbol'], interface=s['interface'], source=s['source'])
                                              for s in observed]
                row['configured'] = policies.get(key, {})
                record['procedures'].append(row)
            modules.append(record)
    rows = [p for m in modules for p in m['procedures']]
    return dict(schema_version=1, kind_name=kind, modules=modules, errors=errors, c_symbols=symbols,
                incomplete=bool(errors or any('error' in m for m in modules)),
                coverage=dict(public_exports=len(rows), scalar_supported=sum(p['scalar_supported'] for p in rows),
                              bulk_supported=sum(p['bulk_supported'] for p in rows),
                              observed_c_scalar=sum(any(s['interface']=='scalar' for s in p['observed_c_symbols']) for p in rows),
                              observed_c_procedural=sum(any(s['interface']=='procedural' for s in p['observed_c_symbols']) for p in rows),
                              configured_scalar=sum(bool(p['configured'].get('scalar')) for p in rows),
                              configured_bulk=sum(bool(p['configured'].get('bulk')) for p in rows),
                              python_methods=sum(bool(p['configured'].get('python')) for p in rows),
                              configured_array_subroutines=sum(bool(p['configured'].get('array_subroutine')) for p in rows)))


def draft(module, kind, kind_module, symbol_prefix):
    """No mathematical policies or Python defaults can be executable until reviewed."""
    generate.identifier(module['module'])
    if 'error' in module:
        raise ValueError(f"cannot draft unresolved module {module['module']}")
    functions, excluded, docs = {}, {}, {}
    for row in module['procedures']:
        name = row['name']
        if not row['scalar_supported']:
            excluded[name] = row['unsupported_reason']
            continue
        scalar = symbol_prefix + name + '_f64'
        generate.identifier(scalar)
        entry = dict(symbol=scalar)
        if row['optional_arguments']:
            entry['optional_policy'] = 'REVIEW_REQUIRED'
        if row['bulk_supported']:
            entry.update(vector_symbol=symbol_prefix+name+'_vec_f64', array_argument=None,
                         finite='REVIEW_REQUIRED', positive='REVIEW_REQUIRED')
        functions[name] = entry
        docs.update((a['name'], 'REVIEW_REQUIRED: meaning, valid range and units') for a in row['arguments'])
    name = module['module']
    return dict(schema_version=2, draft=True,
                review_required=['Resolve existing C symbols from the inventory; do not redeclare them.',
                                 'Choose the varying input, uniform constraints and optional-argument policy.',
                                 'Supply argument documentation and explicit Python defaults separately.',
                                 'Review exclusions and all output/module/symbol names.',
                                 'Remove draft and review_required only after every policy is resolved.'],
                source=module['source'], module=name, kind_module=kind_module, kind_name=kind,
                wrapper_module='cg_'+name, header_guard=('CG_'+name+'_H').upper(),
                outputs=dict(fortran=f'REVIEW_REQUIRED/{name}.f90', header=f'REVIEW_REQUIRED/{name}.h',
                             coverage=f'REVIEW_REQUIRED/{name}.json'),
                argument_docs=docs, functions=functions, excluded=excluded)


def markdown(report):
    lines = ['# Fortran / C / Python API inventory', '',
             'Counts apply only to the inspected sources. Unresolved modules make the report incomplete.',
             'Observed C calls are syntactic links through explicit ONLY imports, not ABI or mathematical validation.',
             'Configured coverage records policy, not proof that binaries or Python packages are installed.', '',
             f"Incomplete: {report['incomplete']}", '', '| Coverage | Count |', '| --- | ---: |']
    lines += [f'| {k} | {v} |' for k, v in report['coverage'].items()]
    lines += ['', '| Module | Public export | Scalar pattern | Bulk pattern | Configured scalar | Configured bulk | Array adapter | Python |',
              '| --- | --- | --- | --- | --- | --- | --- | --- |']
    for module in report['modules']:
        for row in module['procedures']:
            configured = row['configured']
            lines.append(f"| {module['module']} | {row['name']} | {row['scalar_supported']} | "
                         f"{row['bulk_supported']} | {configured.get('scalar') or ''} | "
                         f"{configured.get('bulk') or ''} | {configured.get('array_subroutine') or ''} | {configured.get('python') or ''} |")
        if 'error' in module:
            lines += ['', f"{module['module']}: {module['error']}"]
    for error in report['errors']:
        lines += ['', f"{error['source']}: {error['reason']}"]
    lines += ['', 'Argument types, ranks, optional flags, exclusion reasons and observed C symbols are in the JSON.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sources', nargs='+', type=Path, help='native source files or directories (recursive *.f90)')
    parser.add_argument('--abi-source', action='append', type=Path, default=[])
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--kind', default='dp')
    parser.add_argument('--kind-module', default=None)
    parser.add_argument('--symbol-prefix', default=None)
    parser.add_argument('--manifest-batch', type=Path)
    parser.add_argument('--python-policy', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--draft-directory', type=Path)
    parser.add_argument('--progress', action='store_true', help='print each native file before parsing')
    parser.add_argument('--max-source-bytes', type=int, help='skip larger native files and mark inventory incomplete')
    args = parser.parse_args()
    try:
        if args.max_source_bytes is not None and args.max_source_bytes <= 0:
            raise ValueError('max-source-bytes must be positive')
        root = args.root.resolve()
        sources, abi = files(args.sources), files(args.abi_source)
        policies = configured(root, args.manifest_batch, args.python_policy)
        report = inventory(sources, root, args.kind, abi, policies, args.max_source_bytes, args.progress)
        outputs = {args.output.resolve(): json.dumps(report, indent=2, sort_keys=True)+'\n',
                   args.output.with_suffix('.md').resolve(): markdown(report)}
        if len(outputs) != 2:
            raise ValueError('output must have a non-Markdown suffix (normally .json)')
        if args.draft_directory:
            if not args.kind_module or not args.symbol_prefix:
                raise ValueError('drafts require explicit --kind-module and --symbol-prefix')
            generate.identifier(args.kind_module)
            for module in report['modules']:
                if 'error' not in module:
                    path = (args.draft_directory / (module['module']+'.json')).resolve()
                    if path in outputs:
                        raise ValueError(f'duplicate draft/output path: {path}')
                    if path.exists():
                        raise ValueError(f'refusing to overwrite draft: {path}')
                    outputs[path] = json.dumps(draft(module, args.kind, args.kind_module, args.symbol_prefix), indent=2)+'\n'
        protected = set(sources+abi)
        protected.update(p.resolve() for p in (args.manifest_batch, args.python_policy) if p)
        if args.manifest_batch:
            protected.update((args.manifest_batch.parent / p).resolve()
                             for p in json.loads(args.manifest_batch.read_text(encoding='utf-8')))
        if args.python_policy:
            python = json.loads(args.python_policy.read_text(encoding='utf-8'))
            protected.update((root / p['manifest']).resolve() for p in python['groups'].values())
        if protected & outputs.keys():
            raise ValueError('inventory outputs must not overwrite sources or policies')
        for path, value in outputs.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(value, encoding='utf-8', newline='\n')
        print(json.dumps(report['coverage'], sort_keys=True))
        return 1 if report['incomplete'] else 0
    except (ValueError, OSError, generate.FortranSyntaxError) as error:
        print(f'cabi_generator inventory: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
