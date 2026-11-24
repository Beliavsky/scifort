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
        if not isinstance(settings, dict) or set(settings) != {'manifests', 'defaults'}:
            raise ValueError('array groups need manifests and defaults')
        paths, defaults = settings['manifests'], settings['defaults']
        if (not isinstance(paths, list) or not paths or not all(isinstance(p, str) for p in paths)
                or len(set(paths)) != len(paths) or not isinstance(defaults, dict)):
            raise ValueError('array manifests must be a unique nonempty list; defaults must be an object')
        specs, methods, used_defaults = {}, [], set()
        for path in paths:
            manifest = (root/path).resolve()
            config = json.loads(manifest.read_text(encoding='utf-8'))
            generate.generate(config, root)
            source = (root/config['source']).resolve()
            protected.update({manifest, source, *((root/p).resolve() for p in config['outputs'].values())})
            api = generate.public_api(source, config['module'])
            entries = generate.explicit_arrays.prepare(config, api, generate.identifier,
                                                       source.read_text(encoding='utf-8').splitlines())
            if not entries:
                raise ValueError('array manifests need selected subroutines')
            for entry in entries:
                name = py_name(entry['name'])
                if name in specs or name in {'self', 'out'} or hasattr(object, name):
                    raise ValueError(f'duplicate or reserved array method: {name}')
                contract = entry.get('shape_contract')
                if not contract:
                    raise ValueError(f'{name}: Python arrays require a reviewed structured shape contract')
                status = entry['policy']['status_argument']
                metadata = entry['metadata']
                inputs = [a for a in entry['args'] if a != status and metadata[a]['intent'] == 'in']
                outputs = [a for a in entry['args'] if a != status and metadata[a]['intent'] == 'out']
                if set(inputs + outputs + [status]) != set(entry['args']):
                    raise ValueError('Python array adapters do not yet support native intent(inout)')
                for arg in entry['args']:
                    py_name(arg)
                    if arg in {'self', 'out'}:
                        raise ValueError(f'{arg}: reserved Python argument')
                for dim, definition in contract['definitions'].items():
                    if definition['argument'] not in inputs:
                        raise ValueError(f'{name}: dimension {dim} must be determined by an input array')
                values = defaults.get(name, {})
                if name in defaults:
                    used_defaults.add(name)
                if not isinstance(values, dict) or set(values) - set(inputs):
                    raise ValueError(f'{name}: defaults must name scalar inputs')
                for arg, value in values.items():
                    info = metadata[arg]
                    valid = type(value) is bool if info['category'] == 'logical' else (
                        type(value) in (int, float) and math.isfinite(value))
                    if info['rank'] or not valid:
                        raise ValueError(f'{name}: defaults must match scalar input types')
                parameters, seen_default = [], False
                for arg in inputs:
                    if arg in values:
                        parameters.append(f'{arg}={values[arg]!r}')
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
                specs[name] = spec
                actuals = ', '.join(f'{a!r}: {a}' for a in inputs)
                signature = ', '.join(['self', *parameters, '*', 'out=None'])
                methods += [f'    def {name}({signature}):',
                            '        """Return ArrayResult(status, outputs); retain native failure outputs."""',
                            f'        return self._invoke({name!r}, {{{actuals}}}, out)', '']
        if set(defaults) != used_defaults:
            raise ValueError('unused array method defaults')
        cls = group.capitalize()
        lines += [f'class {cls}(ArrayGroup):',
                  '    _specs = ' + pformat(specs, width=96).replace('\n', '\n    '), '', *methods]
        groups.append((group, cls))
    return lines, groups, protected
