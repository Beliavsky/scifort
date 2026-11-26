# SPDX-License-Identifier: MIT
"""Generate Python array-procedure metadata from reviewed native/C contracts."""
import json
import math
from pprint import pformat

import generate


def prepare(policy, root, py_name):
    options = policy.get('array_groups', {})
    if not isinstance(options, dict):
        raise ValueError('array_groups must be an object')
    lines, groups, protected = [], [], set()
    for group, settings in sorted(options.items()):
        py_name(group)
        if group in policy['groups'] or group in {'api', 'library', 'distribution', 'path', 'arraygroup'}:
            raise ValueError(f'{group}: reserved or duplicate array group')
        if (not isinstance(settings, dict) or set(settings) -
                {'prefix', 'dimension_defaults', 'dimension_sources', 'zero_based_indices', 'output_dimensions'} !=
                {'manifests', 'defaults'}):
            raise ValueError('array groups need manifests and defaults')
        prefix = settings.get('prefix', '')
        if not isinstance(prefix, str):
            raise ValueError('array group prefix must be a string')
        paths, defaults = settings['manifests'], settings['defaults']
        if (not isinstance(paths, list) or not paths or not all(isinstance(p, str) for p in paths)
                or len(set(paths)) != len(paths) or not isinstance(defaults, dict)):
            raise ValueError('array manifests must be a unique nonempty list; defaults must be an object')
        dimension_defaults = settings.get('dimension_defaults', {})
        if not isinstance(dimension_defaults, dict):
            raise ValueError('dimension_defaults must map native procedures to dimension values')
        sources = settings.get('dimension_sources', {})
        indices = settings.get('zero_based_indices', {})
        output_dimensions = settings.get('output_dimensions', {})
        if not isinstance(output_dimensions, dict):
            raise ValueError('output_dimensions must map native procedures to dimension defaults')
        if not isinstance(sources, dict) or not isinstance(indices, dict):
            raise ValueError('dimension_sources and zero_based_indices must map native procedures')
        used_sources, used_indices = set(), set()
        used_output_dimensions = set()
        specs, methods, used_defaults, used_dimension_defaults = {}, [], set(), set()
        for path in paths:
            manifest = (root/path).resolve()
            config = json.loads(manifest.read_text(encoding='utf-8'))
            generate.generate(config, root)
            source = (root/config['source']).resolve()
            protected.update({manifest, source, *((root/p).resolve() for p in config['outputs'].values())})
            api = generate.public_api(source, config['module'], specifics=config.get('generic_specifics'))
            entries = generate.explicit_arrays.prepare(config, api, generate.identifier,
                                                       source.read_text(encoding='utf-8').splitlines())
            if not entries:
                raise ValueError('array manifests need selected subroutines')
            for entry in entries:
                native_name = entry['name']
                if not native_name.startswith(prefix):
                    raise ValueError('array method must start with its configured prefix')
                name = py_name(native_name[len(prefix):])
                if name in specs or name in {'self', 'out'} or hasattr(object, name):
                    raise ValueError(f'duplicate or reserved array method: {name}')
                contract = entry.get('shape_contract')
                if not contract:
                    raise ValueError(f'{name}: Python arrays require a reviewed structured shape contract')
                status = entry['policy']['status_argument']
                metadata = entry['metadata']
                optional_arrays = entry.get('optional_arrays', {})
                flags = set(optional_arrays.values())
                inputs = [a for a in entry['args'] if a != status and a not in flags
                          and (metadata[a]['intent'] == 'in' or metadata[a]['category'] == 'state')]
                outputs = [a for a in entry['args'] if a != status and metadata[a]['intent'] == 'out']
                if set(inputs + outputs + [status]) | flags != set(entry['args']):
                    raise ValueError('Python array adapters do not yet support native intent(inout)')
                for arg in entry['args']:
                    py_name(arg)
                    if arg in {'self', 'out'}:
                        raise ValueError(f'{arg}: reserved Python argument')
                output_sizes = output_dimensions.get(native_name, {})
                if native_name in output_dimensions:
                    used_output_dimensions.add(native_name)
                if not isinstance(output_sizes, dict) or (native_name in output_dimensions and not output_sizes):
                    raise ValueError('output dimensions need a nonempty dimension/default object')
                for dim, value in output_sizes.items():
                    py_name(dim)
                    if dim in {'self', 'out'}:
                        raise ValueError('output dimension name conflicts with a Python parameter')
                    definition = contract['definitions'].get(dim, {})
                    if (definition.get('argument') not in outputs or dim in entry['args']
                            or type(value) is not int or not entry['dimensions'][dim] <= value <= 2147483647):
                        raise ValueError('output dimension must be defined by an output array and have a valid default')
                for dim, definition in contract['definitions'].items():
                    if dim in output_sizes:
                        continue
                    if definition['argument'] not in inputs or definition['argument'] in optional_arrays:
                        raise ValueError(f'{name}: dimension {dim} must be determined by an input array')
                values = defaults.get(native_name, {})
                if native_name in defaults:
                    used_defaults.add(native_name)
                if not isinstance(values, dict) or set(values) - set(inputs):
                    raise ValueError(f'{name}: defaults must name scalar inputs')
                for arg, value in values.items():
                    if metadata[arg]['category'] == 'state':
                        raise ValueError('state arguments cannot have defaults')
                    if arg in optional_arrays:
                        raise ValueError('presence-controlled arrays use None, not configured fills/defaults')
                    info = metadata[arg]
                    if info['rank']:
                        fill = value.get('fill') if isinstance(value, dict) else None
                        shape = entry['policy']['shapes'][arg]
                        if (not info.get('optional') or not isinstance(value, dict) or set(value) != {'fill'}
                                or not (fill == 'zeros' or (fill == 'identity' and len(shape) == 2
                                                           and shape[0] == shape[1]))):
                            raise ValueError(f'{name}: array defaults need an optional array and zeros or square identity fill')
                        continue
                    valid = (type(value) is bool if info['category'] == 'logical' else
                             type(value) is int and -2147483648 <= value <= 2147483647
                             if info['category'] == 'integer' else
                             type(value) in (int, float) and math.isfinite(value))
                    if not valid:
                        raise ValueError(f'{name}: defaults must match scalar input types')
                parameters, seen_default = [], False
                for arg in inputs:
                    if arg in optional_arrays:
                        parameters.append(f'{arg}=None')
                        seen_default = True
                    elif arg in values:
                        default = None if metadata[arg]['rank'] else values[arg]
                        parameters.append(f'{arg}={default!r}')
                        seen_default = True
                    else:
                        if seen_default:
                            raise ValueError('required array inputs must precede defaults')
                        parameters.append(arg)
                spec = dict(symbol=entry['policy']['symbol'], arguments=entry['args'], status=status,
                            metadata=metadata, inputs=inputs, outputs=outputs,
                            dimensions=list(entry['dimensions']), shapes=entry['policy']['shapes'],
                            definitions=contract['definitions'], minimums=contract['minimums'],
                            logical_inputs=entry['policy'].get('logical_inputs', {}))
                if optional_arrays:
                    spec['optional_arrays'] = optional_arrays
                if entry.get('state_codecs'):
                    spec['state_codecs'] = entry['state_codecs']
                if output_sizes:
                    spec['output_dimensions'] = output_sizes
                array_defaults = {a: v['fill'] for a, v in values.items() if metadata[a]['rank']}
                if array_defaults:
                    spec['array_defaults'] = array_defaults
                fallback = dimension_defaults.get(native_name, {})
                if native_name in dimension_defaults:
                    used_dimension_defaults.add(native_name)
                expected_fallback = {dim for dim, definition in contract['definitions'].items()
                                     if definition['argument'] in array_defaults}
                if (not isinstance(fallback, dict) or set(fallback) != expected_fallback
                        or any(type(value) is not int or not contract['minimums'][dim] <= value <= 2147483647
                               for dim, value in fallback.items())):
                    raise ValueError('dimension defaults must cover exactly defaulted dimension sources with valid C extents')
                if fallback:
                    spec['dimension_defaults'] = fallback
                if native_name in sources:
                    used_sources.add(native_name)
                    rules = sources[native_name]
                    if not isinstance(rules, dict) or not rules or set(rules) - set(fallback):
                        raise ValueError('dimension sources need dimensions with explicit fallbacks')
                    for dim, candidates in rules.items():
                        if (not isinstance(candidates, list) or len(candidates) < 2
                                or candidates[0] != contract['definitions'][dim]):
                            raise ValueError('ordered dimension sources must begin with the declared source')
                        seen = []
                        for candidate in candidates:
                            if not isinstance(candidate, dict) or set(candidate) != {'argument', 'axis'}:
                                raise ValueError('dimension sources need argument/axis objects')
                            arg, axis = candidate['argument'], candidate['axis']
                            info = metadata.get(arg, {}) if isinstance(arg, str) else {}
                            rank = info.get('rank', 0)
                            if (not isinstance(arg, str) or arg not in array_defaults
                                    or not rank or candidate in seen
                                    or not (axis is None and rank == 1 or
                                            type(axis) is int and 1 <= axis <= rank)
                                    or entry['policy']['shapes'][arg][(axis or 1)-1] != dim):
                                raise ValueError('dimension source must be a distinct matching optional array axis')
                            seen.append(candidate)
                    spec['dimension_sources'] = rules
                if native_name in indices:
                    used_indices.add(native_name)
                    rules = indices[native_name]
                    if not isinstance(rules, dict) or not rules:
                        raise ValueError('zero_based_indices needs input-array/dimension mappings')
                    for arg, dim in rules.items():
                        info = metadata.get(arg, {})
                        if (arg not in inputs or info.get('category') != 'integer'
                                or info.get('rank') != 1 or not isinstance(dim, str)
                                or dim not in entry['dimensions'] or arg in values):
                            raise ValueError('zero-based indices need required rank-one integer inputs and a bound dimension')
                        if any(arg in entry['policy'].get(key, {}) for key in ('integer_limits', 'integer_sums')):
                            raise ValueError('zero-based indices cannot combine with native integer guards')
                    spec['zero_based_indices'] = rules
                if 'result_argument' in entry:
                    spec['result_argument'] = entry['result_argument']
                if entry.get('synthetic_status'):
                    spec['synthetic_status'] = True
                if entry.get('length_choices'):
                    spec['length_choices'] = entry['length_choices']
                for key in ('integer_limits', 'integer_sums'):
                    if entry['policy'].get(key):
                        spec[key] = entry['policy'][key]
                specs[name] = spec
                actuals = ', '.join(f'{a!r}: {a}' for a in inputs)
                if output_sizes:
                    actuals += ', ' + ', '.join(f'{d!r}: {d}' for d in output_sizes)
                function_result = 'result_argument' in entry
                scalar_result = function_result and not metadata[entry['result_argument']]['rank']
                keywords = [f'{d}={v}' for d, v in output_sizes.items()]
                if not scalar_result:
                    keywords.append('out=None')
                signature = ', '.join(['self', *parameters] + (['*', *keywords] if keywords else []))
                doc = ('Return the native scalar value.' if scalar_result else
                       'Return the native array result; optionally use caller-owned out.' if function_result else
                       'Return ArrayResult(call status, outputs); native NaN/infinity is retained.' if entry.get('synthetic_status') else
                       'Return ArrayResult(status, outputs); retain native failure outputs.')
                output = 'None' if scalar_result else 'out'
                methods += [f'    def {name}({signature}):',
                            f'        """{doc}"""',
                            f'        return self._invoke({name!r}, {{{actuals}}}, {output})', '']
        if set(dimension_defaults) != used_dimension_defaults:
            raise ValueError('unused dimension defaults')
        if set(sources) != used_sources or set(indices) != used_indices:
            raise ValueError('unused dimension sources or zero-based index policies')
        if set(output_dimensions) != used_output_dimensions:
            raise ValueError('unused output dimension policies')
        if set(defaults) != used_defaults:
            raise ValueError('unused array method defaults')
        cls = group.capitalize()
        lines += [f'class {cls}(ArrayGroup):',
                  '    _specs = ' + pformat(specs, width=96).replace('\n', '\n    '), '', *methods]
        groups.append((group, cls))
    return lines, groups, protected
