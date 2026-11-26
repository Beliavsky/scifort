# SPDX-License-Identifier: MIT
"""Small, non-executing shape-contract prototype for explicit-array adapters."""
import re


def length_choices(policy, dimensions):
    """Validate finite choices of another dimension minus a nonnegative offset."""
    choices = policy.get('length_choices', {})
    if not isinstance(choices, dict):
        raise ValueError('length_choices must be an object')
    for dim, rule in choices.items():
        if (dim not in dimensions or not isinstance(rule, dict)
                or set(rule) != {'dimension', 'offsets'} or rule['dimension'] not in dimensions
                or rule['dimension'] == dim):
            raise ValueError('length choices must relate two distinct declared dimensions')
        offsets = rule['offsets']
        if (not isinstance(offsets, list) or not offsets
                or any(type(v) is not int or not -2147483647 <= v <= 0 for v in offsets)
                or len(set(offsets)) != len(offsets)):
            raise ValueError('length choice offsets must be distinct nonpositive C int values')
    return choices


def choice_text(rule):
    return ','.join(rule['dimension'] + (str(v) if v else '') for v in rule['offsets'])


def prepare(node, source_lines, entry):
    """Validate source annotations against an explicitly reviewed wrapper policy."""
    name, policy = entry['name'], entry['policy']
    definitions, minimums, conditions, choices = {}, {}, [], {}
    # Only the procedure's initial documentation block is eligible for directives.
    for line in source_lines[node.content[0].item.span[1]:]:
        text = line.strip()
        if text and not text.startswith('!'):
            break
        if not text.startswith('!! @'):
            continue
        text = text[3:].strip()
        match = re.fullmatch(r'@dimension\s+(\w+)\s*=\s*size\(\s*(\w+)\s*(?:,\s*([12]))?\s*\)',
                             text, re.I)
        if match:
            dim, arg = match[1].lower(), match[2].lower()
            axis = int(match[3]) if match[3] else None
            if dim in definitions or dim not in entry['dimensions']:
                raise ValueError(f'{name}: duplicate or unknown dimension {dim}')
            rank = entry['metadata'].get(arg, {}).get('rank', 0)
            if not rank or (axis is None and rank != 1) or (axis and axis > rank):
                raise ValueError(f'{name}: dimension definition must reference a valid array axis')
            if policy['shapes'][arg][(axis or 1)-1] != dim:
                raise ValueError(f'{name}: dimension definition disagrees with array shape for {arg}')
            definitions[dim] = dict(argument=arg, axis=axis)
            continue
        match = re.fullmatch(r'@require\s+(\w+)\s*>=\s*(0|[1-9][0-9]*)', text, re.I)
        if match:
            dim, minimum = match[1].lower(), int(match[2])
            if dim in minimums or dim not in entry['dimensions']:
                raise ValueError(f'{name}: duplicate or unknown minimum dimension {dim}')
            if minimum != entry['dimensions'][dim]:
                raise ValueError(f'{name}: minimum for {dim} disagrees with wrapper policy')
            minimums[dim] = minimum
            continue
        match = re.fullmatch(r'@require\s+(\w+)\s+in\s*\(([^)]+)\)', text, re.I)
        if match:
            dim = match[1].lower()
            terms = [re.fullmatch(r'([a-z]\w*)(?:\s*-\s*([1-9][0-9]*))?', t.strip(), re.I)
                     for t in match[2].split(',')]
            if dim in choices or not all(terms) or len({t[1].lower() for t in terms}) != 1:
                raise ValueError(f'{name}: unsupported or repeated dimension choices')
            rule = dict(dimension=terms[0][1].lower(), offsets=[-int(t[2] or 0) for t in terms])
            length_choices(dict(length_choices={dim:rule}), entry['dimensions'])
            choices[dim] = rule
            continue
        match = re.fullmatch(r'@when\s+(\w+)\s*==\s*(0|[1-9][0-9]*)\s*:\s*(.+)', text, re.I)
        if match:
            dim, value, description = match[1].lower(), int(match[2]), match[3].strip()
            if dim not in entry['dimensions'] or value < entry['dimensions'][dim]:
                raise ValueError(f'{name}: conditional note must reference an allowed dimension value')
            if not description.isascii() or any(c in description for c in '|`<>'):
                raise ValueError(f'{name}: conditional note must be plain ASCII text')
            record = dict(dimension=dim, value=value, description=description)
            if record in conditions:
                raise ValueError(f'{name}: repeated conditional note')
            conditions.append(record)
            continue
        raise ValueError(f'{name}: unsupported shape-contract annotation: {text}')
    expected = set(entry['dimensions'])
    if set(definitions) != expected or set(minimums) != expected:
        raise ValueError(f'{name}: shape contract needs one definition and minimum for every dimension')
    if entry['documentation_review']:
        raise ValueError(f'{name}: shape contract requires matching trailing shape comments')
    if choices != entry.get('length_choices', {}):
        raise ValueError(f'{name}: source length choices disagree with wrapper policy')
    contract = dict(schema_version=1, definitions=definitions, minimums=minimums, conditions=conditions)
    if choices:
        contract['length_choices'] = choices
    return contract


def markdown(config, entries):
    lines = ['# Native shape contracts', '', '<!-- Generated by cabi_generator; do not edit. -->', '',
             f"Source: `{config['source']}`; module: `{config['module']}`.", '',
             'These declared contracts are checked against reviewed C adapter policies.',
             'They do not prove that the numerical implementation obeys the contract.', '']
    for entry in entries:
        specific = config.get('generic_specifics', {}).get(entry['name'])
        if specific:
            lines += [f"Public generic `{entry['name']}`: reviewed specific `{specific}`.", '']
        contract = entry.get('shape_contract')
        if not contract:
            continue
        lines += [f"## {entry['name']}", '', entry['policy']['description'], '',
                  '| Dimension | Definition | Minimum |', '| --- | --- | ---: |']
        for dim in entry['dimensions']:
            definition = contract['definitions'][dim]
            axis = f",{definition['axis']}" if definition['axis'] else ''
            lines.append(f"| `{dim}` | `size({definition['argument']}{axis})` | {contract['minimums'][dim]} |")
        lines += ['', '| Argument | Declared shape | Required shape | Meaning |',
                  '| --- | --- | --- | --- |']
        for arg in entry['args']:
            if arg not in entry['policy']['shapes']:
                continue
            rank = entry['metadata'][arg]['rank']
            required = ','.join(entry['policy']['shapes'][arg])
            declared = ','.join(entry['result_extents'] if arg == entry.get('result_argument')
                                and 'result_extents' in entry else [':']*rank)
            meaning = config['argument_docs'][arg].replace('|', r'\|')
            lines.append(f"| `{arg}` | `({declared})` | `({required})` | {meaning} |")
        lines.append('')
        if entry.get('synthetic_status'):
            lines += ['The adapter appends a call-level status; the native subroutine has no status argument.',
                      'Status 0 retains native numerical NaN/infinity; invalid extents preserve outputs.', '']
        if entry.get('result_extents'):
            lines += [f"`{entry['result_argument']}` is the adapter output for the native array-valued result.", '']
        for dim, rule in contract.get('length_choices', {}).items():
            lines += [f"Required length choice: `{dim} in ({choice_text(rule)})`.", '']
        for arg, bounds in entry['policy'].get('integer_limits', {}).items():
            lines += [f"Adapter integer bounds: `{bounds['minimum']} <= {arg} <= {bounds['maximum']}`.", '']
        for arg, ceiling in entry['policy'].get('integer_sums', {}).items():
            lines += [f"Adapter widened sum bound: `sum(abs({arg})) <= {ceiling}`.", '']
        for note in contract['conditions']:
            lines += [f"When `{note['dimension']} == {note['value']}`: {note['description']}", '']
    return '\n'.join(lines)
