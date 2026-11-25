#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Regenerate/check reviewed C/Python interfaces, dependency build plan and coverage together."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

from fparser.common.readfortran import FortranStringReader
from fparser.two import Fortran2003 as F
from fparser.two import Fortran2008 as F08
from fparser.two.parser import ParserFactory
from fparser.two.utils import NoMatchError, walk
import generate
import generate_python

INTRINSIC_MODULES = {'iso_c_binding', 'iso_fortran_env', 'ieee_arithmetic',
                     'ieee_exceptions', 'ieee_features'}


def relative(path, root):
    return Path(os.path.relpath(path, root)).as_posix()


def read(path, overlay):
    return overlay[path] if path in overlay else path.read_text(encoding='utf-8')


def providers(source_directories, root, overlay):
    """Index MODULE statements with the logical-line reader, without parsing numerical bodies."""
    paths = set()
    for directory in source_directories:
        directory = (root / directory).resolve()
        if not directory.is_dir():
            raise ValueError(f'source directory does not exist: {directory}')
        paths.update(p.resolve() for p in directory.rglob('*') if p.is_file() and p.suffix.lower()=='.f90')
    paths.update(p for p in overlay if p.suffix.lower()=='.f90')
    index = {}
    for path in sorted(paths):
        reader = FortranStringReader(read(path, overlay))
        for item in reader:
            line = getattr(item, 'line', '').strip()
            if not line.lower().startswith('module'):
                continue
            try:
                statement = F.Module_Stmt(line)
            except NoMatchError:
                continue  # MODULE PROCEDURE/FUNCTION is not a module definition.
            name = str(statement.items[1]).lower()
            index.setdefault(name, []).append(path)
    return index, paths


def build_plan(roots, source_directories, root, overlay):
    """Resolve a deterministic file-level module closure; never guess a provider from its filename."""
    # The fparser statement classes must be initialized before indexing statements.
    parser = ParserFactory().create(std='f2008')
    index, indexed = providers(source_directories, root, overlay)
    ordered, edges, definitions, visiting, visited, hashes = [], {}, {}, [], set(), {}

    def visit(path):
        path = path.resolve()
        if path in visiting:
            chain = visiting[visiting.index(path):] + [path]
            raise ValueError('module dependency cycle: '+' -> '.join(relative(p, root) for p in chain))
        if path in visited:
            return
        visiting.append(path)
        text = read(path, overlay)
        # INCLUDE and preprocessing need separately tracked inputs; do not silently omit them.
        if any(line.lstrip().startswith('#') for line in text.splitlines()):
            raise ValueError(f'{relative(path, root)}: preprocessing is not supported by this build planner')
        reader = FortranStringReader(text)
        for item in reader:
            line = getattr(item, 'line', '').strip()
            if line.lower().startswith('include'):
                try:
                    F.Include_Stmt(line)
                except NoMatchError:
                    continue
                raise ValueError(f'{relative(path, root)}: INCLUDE inputs need explicit build-planner support')
        try:
            tree = parser(FortranStringReader(text))
        except generate.FortranSyntaxError as error:
            raise ValueError(f'cannot inspect dependency source {relative(path, root)}: {error}') from error
        if walk(tree, F08.Submodule):
            raise ValueError(f'{relative(path, root)}: submodule compilation is not supported')
        modules = walk(tree, F.Module)
        if not modules:
            raise ValueError(f'{relative(path, root)}: a build source must define a module')
        positions = {str(m.content[0].items[1]).lower(): i for i, m in enumerate(modules)}
        dependencies = set()
        for position, module in enumerate(modules):
            name = str(module.content[0].items[1]).lower()
            if len(index.get(name, [])) != 1:
                candidates = ', '.join(relative(p, root) for p in index.get(name, []))
                raise ValueError(f'ambiguous provider for module {name}: {candidates}')
            definitions[name] = relative(path, root)
            for use in walk(module, F.Use_Stmt):
                nature = str(use.items[0]).lower() if use.items[0] else None
                dependency = str(use.items[2]).lower()
                if nature == 'intrinsic':
                    continue
                candidates = index.get(dependency, [])
                if not candidates:
                    if nature != 'non_intrinsic' and dependency in INTRINSIC_MODULES:
                        continue
                    raise ValueError(f'{relative(path, root)}: missing provider for module {dependency}')
                if len(candidates) != 1:
                    raise ValueError(f'{relative(path, root)}: ambiguous provider for module {dependency}: '+
                                     ', '.join(relative(p, root) for p in candidates))
                provider = candidates[0]
                if provider == path:
                    if dependency not in positions or positions[dependency] >= position:
                        raise ValueError(f'{relative(path, root)}: module {dependency} must precede {name} in the same file')
                else:
                    dependencies.add(provider)
        for dependency in sorted(dependencies):
            visit(dependency)
        edges[relative(path, root)] = [relative(p, root) for p in sorted(dependencies)]
        hashes[relative(path, root)] = hashlib.sha256(text.encode('utf-8')).hexdigest()
        visiting.pop()
        visited.add(path)
        ordered.append(relative(path, root))

    for path in sorted(set(roots)):
        visit(path)
    return dict(schema_version=1, sources=ordered, dependencies=edges,
                modules=dict(sorted(definitions.items())), source_sha256=hashes), indexed


def coverage(configs, python, root):
    rows = {}
    for config in configs:
        module = config['module']
        for name, reason in config['excluded'].items():
            rows[(module, name)] = dict(module=module, procedure=name, c_manifest_exclusion=reason,
                                       scalar_generated=False, bulk_generated=False)
        for name, entry in config['functions'].items():
            row = rows.setdefault((module, name), dict(module=module, procedure=name))
            row.update(scalar=entry.get('symbol') or entry.get('existing_scalar_symbol'),
                       bulk=entry.get('vector_symbol'), scalar_generated=bool(entry.get('symbol')),
                       bulk_generated=bool(entry.get('vector_symbol')))
            if entry.get('elementwise_symbol'):
                row['elementwise'] = entry['elementwise_symbol']
        for name, entry in dict(config.get('subroutines', {}), **config.get('array_functions', {})).items():
            rows[(module, name)] = dict(module=module, procedure=name, array_adapter=entry['symbol'])
    for group, options in sorted(python['groups'].items()):
        config = json.loads((root/options['manifest']).read_text(encoding='utf-8'))
        for name, entry in dict(config['functions'], **options['existing_functions']).items():
            key = (config['module'], name)
            if key not in rows:
                raise ValueError(f'Python manifest must be included in the C batch: {options["manifest"]}')
            row = rows[key]
            row.update(scalar=entry.get('symbol') or entry.get('existing_scalar_symbol'), bulk=entry['vector_symbol'],
                       python=f'{group}.{name[len(options["prefix"]):]}')
    for group, settings in python.get('array_groups', {}).items():
        for manifest in settings['manifests']:
            native = json.loads((root/manifest).read_text(encoding='utf-8'))
            for name in dict(native.get('subroutines', {}), **native.get('array_functions', {})):
                key = (native['module'], name)
                if key not in rows:
                    raise ValueError(f'Python array manifest must be included in the C batch: {manifest}')
                rows[key]['python'] = f"{group}.{name[len(settings.get('prefix', '')):]}"
    procedures = [rows[k] for k in sorted(rows)]
    counts = dict(exports=len(rows), scalar=sum(bool(r.get('scalar')) for r in procedures),
                  bulk=sum(bool(r.get('bulk')) for r in procedures),
                  array_adapters=sum(bool(r.get('array_adapter')) for r in procedures),
                  elementwise=sum(bool(r.get('elementwise')) for r in procedures),
                  python=sum(bool(r.get('python')) for r in procedures))
    return dict(schema_version=1, scope='Reviewed C batch and Python policy only; not a whole-project inventory',
                counts=counts, procedures=procedures)


def coverage_markdown(report):
    lines = ['# Reviewed interface coverage', '', report['scope']+'.', '',
             '| Interface | Count |', '| --- | ---: |']
    lines += [f'| {key} | {value} |' for key, value in report['counts'].items()]
    lines += ['', '| Module | Export | Scalar C | Bulk C | Array adapter | Elementwise C | Python |',
              '| --- | --- | --- | --- | --- | --- | --- |']
    for row in report['procedures']:
        lines.append('| '+' | '.join(str(row.get(k) or '') for k in
                                    ('module','procedure','scalar','bulk','array_adapter','elementwise','python'))+' |')
    lines += ['', 'C-manifest exclusion reasons and generated/existing symbol distinctions are in the JSON.', '']
    return '\n'.join(lines)


def artifacts(config_path, root, *, check=False):
    config_path = config_path.resolve()
    config = json.loads(config_path.read_text(encoding='utf-8'))
    fields = {'schema_version','c_manifests','python_policy','source_directories','build_manifest',
              'coverage_report','aggregate_header','aggregate_header_guard'}
    if not isinstance(config, dict) or set(config) != fields or type(config['schema_version']) is not int or config['schema_version'] != 1:
        raise ValueError('pipeline requires schema_version=1 and explicit C/Python/build/coverage/header settings')
    if (not isinstance(config['source_directories'], list) or not config['source_directories']
            or not all(isinstance(p, str) for p in config['source_directories'])):
        raise ValueError('source_directories must be a nonempty list of directories')
    if not all(isinstance(config[k], str) for k in fields - {'schema_version','source_directories'}):
        raise ValueError('pipeline paths and header guard must be strings')
    paths = {k:(root/config[k]).resolve() for k in
             ('c_manifests','python_policy','build_manifest','coverage_report','aggregate_header')}
    results = generate.collect(paths['c_manifests'], root, True)
    batch = json.loads(paths['c_manifests'].read_text(encoding='utf-8'))
    configs_paths = [(paths['c_manifests'].parent/p).resolve() for p in batch]
    configs = [json.loads(p.read_text(encoding='utf-8')) for p in configs_paths]
    python = json.loads(paths['python_policy'].read_text(encoding='utf-8'))
    python_outputs = generate_python.artifacts(paths['python_policy'], root)
    if results.keys() & python_outputs.keys():
        raise ValueError('C and Python output paths collide')
    results.update(python_outputs)
    # A Python group must refer to a reviewed manifest in the same batch.
    for options in python['groups'].values():
        if (root/options['manifest']).resolve() not in configs_paths:
            raise ValueError('each Python manifest must appear in the C batch')
    roots = [p for p in python_outputs if p.suffix=='.f90']
    for options in python['groups'].values():
        native = json.loads((root/options['manifest']).read_text(encoding='utf-8'))
        roots.append((root/native['outputs']['fortran']).resolve())
    for settings in python.get('array_groups', {}).values():
        for manifest in settings['manifests']:
            path = (root/manifest).resolve()
            if path not in configs_paths:
                raise ValueError('each Python array manifest must appear in the C batch')
            native = json.loads(path.read_text(encoding='utf-8'))
            roots.append((root/native['outputs']['fortran']).resolve())
    plan, indexed = build_plan(roots, config['source_directories'], root, results)
    plan.update(families=sorted(python['groups']), pipeline_config=relative(config_path, root))
    if python.get('array_groups'):
        plan['array_groups'] = sorted(python['array_groups'])
    report = coverage(configs, python, root)
    headers = sorted((root/c['outputs']['header']).resolve() for c in configs)
    guard = config['aggregate_header_guard']
    if not isinstance(guard, str) or not guard.isascii() or not guard.replace('_','A').isalnum() or not guard[0].isalpha():
        raise ValueError('aggregate_header_guard must be a C identifier')
    includes = [relative(p, paths['aggregate_header'].parent) for p in headers]
    if any(any(c in name for c in ('"','\n','\r')) for name in includes):
        raise ValueError('header paths cannot contain quotes or newlines')
    header = ['/* SPDX-License-Identifier: MIT */','/* Generated by cabi_generator; do not edit. */',
              f'#ifndef {guard}',f'#define {guard}',*[f'#include "{p}"' for p in includes],f'#endif /* {guard} */','']
    additional = [(paths['build_manifest'],json.dumps(plan, indent=2, sort_keys=True)+'\n'),
                  (paths['coverage_report'],json.dumps(report, indent=2, sort_keys=True)+'\n'),
                  (paths['coverage_report'].with_suffix('.md'),coverage_markdown(report)),
                  (paths['aggregate_header'],'\n'.join(header))]
    protected = {config_path,paths['c_manifests'],paths['python_policy'],*configs_paths,
                 Path(__file__).resolve(),Path(generate.__file__).resolve(),
                 Path(generate_python.__file__).resolve(),Path(generate_python.generate_arrays.__file__).resolve(),Path(generate_python.__file__).with_name('python_runtime.py').resolve()}
    protected.update((root/c['source']).resolve() for c in configs)
    # Original Fortran files discovered in the source roots must never become outputs.
    protected.update(path for path in indexed if path not in results)
    for path in indexed:
        if not check and path in results and path.exists() and 'Generated by cabi_generator; do not edit.' not in path.read_text(encoding='utf-8'):
            raise ValueError(f'refusing to overwrite original source: {path}')
    for path, text in additional:
        if (not check and path.suffix.lower()=='.h' and path.exists()
                and 'Generated by cabi_generator; do not edit.' not in path.read_text(encoding='utf-8')):
            raise ValueError(f'refusing to overwrite manual header: {path}')
        if path in results:
            raise ValueError(f'pipeline output path collision: {path}')
        results[path] = text
    if results.keys() & protected:
        raise ValueError('pipeline outputs must not overwrite sources, templates or policies')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path)
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    try:
        outputs = artifacts(args.config, args.root.resolve(), check=args.check)
        generate.print_documentation_reviews(outputs)
        if args.check:
            stale = [str(p) for p, text in outputs.items() if not p.exists() or p.read_text(encoding='utf-8') != text]
            if stale:
                raise ValueError('missing or stale pipeline artifacts: '+', '.join(stale))
        else:
            for path, text in outputs.items():
                if path.exists() and path.read_text(encoding='utf-8') == text:
                    continue
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding='utf-8', newline='\n')
        print(f'{"Checked" if args.check else "Generated"} {len(outputs)} pipeline artifacts')
        return 0
    except (ValueError,OSError,generate.FortranSyntaxError) as error:
        print(f'cabi_generator pipeline: {error}', file=sys.stderr)
        return 1


if __name__=='__main__':
    sys.exit(main())
