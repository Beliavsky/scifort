# SPDX-License-Identifier: MIT
"""Validated explicit-shape adapters for native assumed-shape subroutines."""
import re

import shape_contracts

from fparser.two import Fortran2003 as F


def prepare(config, api, identifier, source_lines):
    entries = []
    policies = dict(config.get("subroutines", {}), **config.get("array_functions", {}))
    for name, policy in sorted(policies.items()):
        identifier(name)
        fields = {"symbol", "dimensions", "shapes", "status_argument", "invalid_status", "description"}
        is_function = name in config.get('array_functions', {})
        if is_function:
            fields |= {'result_argument', 'optional_policy'}
        optional_fields = {'logical_inputs', 'shape_contract', 'length_choices',
                           'integer_limits', 'integer_sums', 'synthetic_status', 'optional_arrays',
                           'state_codecs', 'optional_status', 'required_optional_scalars'}
        if not isinstance(policy, dict) or set(policy) - optional_fields != fields:
            raise ValueError(f"{name}: subroutine policy requires {sorted(fields)}")
        synthetic_status = policy.get('synthetic_status', False)
        if type(synthetic_status) is not bool or (is_function and 'synthetic_status' in policy):
            raise ValueError(f'{name}: synthetic_status must be boolean and subroutine-only')
        optional_status = policy.get('optional_status', False)
        if ('optional_status' in policy and (optional_status is not True or is_function or synthetic_status)):
            raise ValueError(f'{name}: optional_status must be true for a native subroutine status')
        codecs = policy.get('state_codecs', {})
        if (not isinstance(codecs, dict) or ('state_codecs' in policy and not codecs)):
            raise ValueError(f'{name}: state_codecs must be a nonempty policy')
        for arg, codec in codecs.items():
            identifier(arg)
            if (not isinstance(codec, dict) or set(codec) - {'promote_optional'} !=
                    {'module', 'type_module', 'type', 'load', 'store', 'length', 'invalid_status'}):
                raise ValueError(f'{name}: state codec needs module/type/load/store/length/status')
            if 'promote_optional' in codec and codec['promote_optional'] is not True:
                raise ValueError(f'{name}: promote_optional must be true when specified')
            for field in ('module', 'type_module', 'type', 'load', 'store'):
                identifier(codec[field])
            if (type(codec['length']) is not int or not 1 <= codec['length'] <= 2147483647
                    or type(codec['invalid_status']) is not int
                    or not 1 <= codec['invalid_status'] <= 2147483647):
                raise ValueError(f'{name}: state length and invalid status must be positive C ints')
        identifier(policy["symbol"])
        node = api[name]
        if is_function:
            if not isinstance(node, F.Function_Subprogram):
                raise ValueError(f'{name}: expected a directly defined function')
            prefix, _, _, suffix = node.content[0].items
            if (suffix is None or suffix.items[0] is None or suffix.items[1] is not None
                    or any(str(p).lower() == 'impure' for p in (prefix.items if prefix else ()))):
                raise ValueError(f'{name}: needs result(...), no bind(c) or impure prefix')
            native_result = identifier(str(suffix.items[0]).lower())
            if policy['optional_policy'] != 'required':
                raise ValueError(f'{name}: optional_policy must be required')
        elif not isinstance(node, F.Subroutine_Subprogram) or node.content[0].items[3] is not None:
            raise ValueError(f"{name}: expected a directly defined non-bind(c) subroutine")
        args = [identifier(str(a).lower()) for a in node.content[0].items[2].items]
        native_args = args[:]
        if set(codecs) - set(native_args):
            raise ValueError(f'{name}: state codecs must name native arguments')
        optional_arrays = policy.get('optional_arrays', {})
        if (not isinstance(optional_arrays, dict)
                or ('optional_arrays' in policy and not optional_arrays)
                or set(optional_arrays) - set(native_args)):
            raise ValueError(f'{name}: optional_arrays must map native array inputs to presence flags')
        promoted = policy.get('required_optional_scalars', [])
        if (not isinstance(promoted, list) or ('required_optional_scalars' in policy and not promoted)
                or any(not isinstance(a, str) or a not in native_args for a in promoted)
                or len(set(promoted)) != len(promoted) or is_function and promoted):
            raise ValueError(f'{name}: required_optional_scalars needs distinct native subroutine inputs')
        flags = []
        for arg, flag in optional_arrays.items():
            identifier(flag)
            if flag in args or flag in flags or flag.startswith('cg_'):
                raise ValueError(f'{name}: optional presence flags must be distinct new argument names')
            flags.append(flag)
        args = [item for arg in args for item in
                ([optional_arrays[arg], arg] if arg in optional_arrays else [arg])]
        if synthetic_status:
            status = identifier(policy['status_argument'])
            if status in args:
                raise ValueError(f'{name}: synthetic status must not name a native argument')
            args.append(status)
        if is_function:
            result = identifier(policy['result_argument'])
            status = identifier(policy['status_argument'])
            if result == status or result in args or status in args:
                raise ValueError(f'{name}: synthetic result/status must be distinct from native arguments')
            args += [result, status]
        dimensions = {}
        if not isinstance(policy["dimensions"], list) or not policy["dimensions"]:
            raise ValueError(f"{name}: dimensions must be an ordered list of name/minimum objects")
        for item in policy["dimensions"]:
            if not isinstance(item, dict) or set(item) != {"name", "minimum"}:
                raise ValueError(f"{name}: each dimension needs name and minimum")
            dim = identifier(item["name"])
            if dim in dimensions:
                raise ValueError(f"{name}: repeated dimension name")
            dimensions[dim] = item["minimum"]
        for dim, minimum in dimensions.items():
            identifier(dim)
            if dim in args or dim.startswith("cg_"):
                raise ValueError(f"{name}: dimension name conflicts with an argument or reserved cg_ prefix")
            if type(minimum) is not int or minimum < 0 or minimum > 2147483647:
                raise ValueError(f"{name}: minimum extents must be nonnegative C int values")
        if any(arg.startswith("cg_") for arg in args):
            raise ValueError(f"{name}: cg_ argument names are reserved")
        invalid = policy["invalid_status"]
        if type(invalid) is not int or not 1 <= invalid <= 2147483647:
            raise ValueError(f"{name}: invalid_status must be a positive C int value")
        description = policy["description"]
        if not isinstance(description, str) or not description or not description.isascii() or "\n" in description:
            raise ValueError(f"{name}: description must be one nonempty ASCII line")
        spec = next(n for n in node.content if isinstance(n, F.Specification_Part))
        declared = {}
        declarations = {}
        for decl in spec.content:
            if isinstance(decl, F.Type_Declaration_Stmt):
                dtype, attributes, entities = decl.items
                attrs = {str(a).lower() for a in attributes.items} if attributes else set()
                for entity in entities.items:
                    declared[str(entity.items[0]).lower()] = (dtype, attrs, entity)
                    declarations[str(entity.items[0]).lower()] = decl
        logical_inputs = policy.get("logical_inputs", {})
        if not isinstance(logical_inputs, dict) or set(logical_inputs) - set(native_args):
            raise ValueError(f"{name}: logical_inputs must name native arguments")
        for arg, mapping in logical_inputs.items():
            if (not isinstance(mapping, dict) or set(mapping) != {"false", "true"}
                    or any(type(v) is not int or not -2147483647 <= v <= 2147483647
                           for v in mapping.values()) or mapping['false'] == mapping['true']):
                raise ValueError(f"{name}: logical mapping needs distinct C int false/true values")
        metadata = {}
        for flag in flags:
            metadata[flag] = dict(intent='in', rank=0, category='integer')
        status = policy["status_argument"]
        if status in codecs:
            raise ValueError(f'{name}: a status output cannot be a state codec')
        if status in logical_inputs:
            raise ValueError(f"{name}: status cannot be a logical mapping")
        if status not in args:
            raise ValueError(f"{name}: status_argument must name a native argument")
        if is_function:
            kind = config['kind_name']
            if native_result not in declared:
                raise ValueError(f'{name}: result needs an explicit declaration')
            dtype, attrs, entity = declared[native_result]
            result_shape = entity.items[1]
            if (str(dtype).lower().replace(' ', '') not in (f'real({kind})', f'real(kind={kind})')
                    or attrs or any(entity.items[2:])
                    or (result_shape and (not isinstance(result_shape, F.Explicit_Shape_Spec_List)
                                          or len(result_shape.items) > 2))):
                raise ValueError(f'{name}: only scalar or explicit rank-one/two real results are supported')
            metadata[result] = dict(intent='out', rank=len(result_shape.items) if result_shape else 0,
                                    category='real')
            declarations[result] = declarations[native_result]
        for arg in args:
            if arg in flags:
                continue
            if (is_function and arg in (result, status)) or (synthetic_status and arg == status):
                continue
            if arg not in declared:
                raise ValueError(f"{name}: {arg} needs an explicit declaration")
            dtype, attrs, entity = declared[arg]
            if arg == status:
                expected_attrs = {'intent(out)', 'optional'} if optional_status else {'intent(out)'}
                if str(dtype).lower() != "integer" or attrs != expected_attrs or any(entity.items[1:]):
                    raise ValueError(f"{name}: status must be a default-integer scalar intent(out)")
                continue
            if arg in codecs:
                codec = codecs[arg]
                codec_attrs = ({'intent(inout)', 'optional'} if codec.get('promote_optional')
                               else {'intent(inout)'})
                if (str(dtype).lower().replace(' ', '') != f"type({codec['type']})"
                        or attrs != codec_attrs or any(entity.items[1:])
                        or arg in optional_arrays or arg in logical_inputs):
                    raise ValueError(f'{name}: state codec requires a matching scalar derived-type intent(inout) policy')
                metadata[arg] = dict(intent='inout', rank=0, category='state')
                if codec.get('promote_optional'):
                    metadata[arg]['optional'] = True
                continue
            intent = next((value for value in ("in", "out", "inout") if f"intent({value})" in attrs), None)
            optional_policy = (arg in logical_inputs or is_function or arg in optional_arrays
                               or arg in promoted)
            base_attrs = attrs - {"optional"} if optional_policy else attrs
            if intent is None or base_attrs != {f"intent({intent})"}:
                raise ValueError(f"{name}: {arg} requires intent and no unconfigured optional attributes")
            if is_function and intent != 'in':
                raise ValueError(f'{name}: array functions require intent(in) inputs')
            if str(dtype).lower() == "logical":
                if (arg not in logical_inputs or intent != "in" or any(entity.items[1:])):
                    raise ValueError(f"{name}: logical inputs need an explicit scalar mapping")
                metadata[arg] = dict(intent=intent, rank=0, category='logical',
                                     optional='optional' in attrs)
                continue
            if arg in logical_inputs:
                raise ValueError(f"{name}: logical_inputs must map logical arguments")
            category = 'integer' if str(dtype).lower() == 'integer' else 'real'
            kind = config["kind_name"]
            if category != 'integer' and str(dtype).lower().replace(" ", "") not in (f"real({kind})", f"real(kind={kind})"):
                raise ValueError(f"{name}: {arg} must be real({kind}) or default integer")
            shape = entity.items[1]
            if any(entity.items[2:]) or (shape and (not isinstance(shape, F.Assumed_Shape_Spec_List)
                                                  or any(str(item) != ":" for item in shape.items))):
                raise ValueError(f"{name}: {arg} must be scalar or assumed-shape with ':' bounds")
            rank = len(shape.items) if shape else 0
            if rank > 2 and not (rank == 3 and category in ('real', 'integer') and intent == 'out'):
                raise ValueError(f"{name}: rank-three arrays require real/integer intent(out); higher ranks are unsupported")
            if category == 'integer' and intent == 'inout':
                raise ValueError(f'{name}: integer intent(inout) arguments are unsupported')
            if arg in promoted and (rank or intent != 'in' or 'optional' not in attrs):
                raise ValueError(f'{name}: promoted optional inputs must be real/integer scalars')
            if arg in optional_arrays and (category != 'real' or not rank or intent != 'in'
                                          or 'optional' not in attrs):
                raise ValueError(f'{name}: optional arrays must be real rank-one/two optional intent(in)')
            metadata[arg] = dict(intent=intent, rank=rank, category=category)
            if (is_function or arg in optional_arrays or arg in promoted) and 'optional' in attrs:
                metadata[arg]['optional'] = True
        for arg in promoted:
            info = metadata.get(arg, {})
            if (info.get('category') not in ('real', 'integer') or info.get('rank') != 0
                    or info.get('intent') != 'in' or not info.get('optional')):
                raise ValueError(f'{name}: promoted optional inputs must be real/integer scalars')
        for arg in optional_arrays:
            info = metadata.get(arg, {})
            if (info.get('category') != 'real' or not info.get('rank')
                    or info.get('intent') != 'in' or not info.get('optional')):
                raise ValueError(f'{name}: optional arrays must be real rank-one/two optional intent(in)')
        shapes = policy["shapes"]
        array_names = {a for a, value in metadata.items() if value["rank"]}
        if not isinstance(shapes, dict) or set(shapes) != array_names or not array_names:
            raise ValueError(f"{name}: shapes must account for every native array")
        used = set()
        for arg, shape in shapes.items():
            if not isinstance(shape, list) or len(shape) != metadata[arg]["rank"]:
                raise ValueError(f"{name}: shape rank mismatch for {arg}")
            if not all(isinstance(dim, str) and dim in dimensions for dim in shape):
                raise ValueError(f"{name}: shapes must reference declared dimension names")
            used.update(shape)
        if used != set(dimensions):
            raise ValueError(f"{name}: unused dimension names")
        choices = shape_contracts.length_choices(policy, dimensions)
        for key in ('integer_limits', 'integer_sums'):
            rules = policy.get(key, {})
            if not isinstance(rules, dict):
                raise ValueError(f'{name}: {key} must be an object')
            for arg, rule in rules.items():
                info = metadata.get(arg, {})
                if info.get('category') != 'integer' or info.get('intent') != 'in':
                    raise ValueError(f'{name}: {key} must name integer inputs')
                if key == 'integer_limits':
                    if (not isinstance(rule, dict) or set(rule) != {'minimum','maximum'}
                            or any(type(v) is not int or not -2147483648 <= v <= 2147483647 for v in rule.values())
                            or rule['minimum'] > rule['maximum']):
                        raise ValueError(f'{name}: integer limits need ordered C int minimum/maximum')
                elif (info['rank'] != 1 or type(rule) is not int or not 0 <= rule <= 2147483647):
                    raise ValueError(f'{name}: integer sums need rank-one inputs and nonnegative C int absolute-sum limits')
        result_extents = []
        if is_function and metadata[result]['rank']:
            for dim, extent in zip(shapes[result], result_shape.items):
                lower, upper = extent.items
                text = str(upper).lower().replace(' ', '')
                match = re.fullmatch(r'size\(([a-z]\w*)(?:,([12]))?\)', text)
                if lower is not None or not match:
                    raise ValueError(f'{name}: result extents must be size(input[,axis]) with default lower bounds')
                source, axis = match[1], int(match[2]) if match[2] else None
                info = metadata.get(source, {})
                if (info.get('intent') != 'in' or not info.get('rank')
                        or (axis is None and info['rank'] != 1) or (axis and axis > info['rank'])
                        or shapes[source][(axis or 1)-1] != dim):
                    raise ValueError(f'{name}: native result extent disagrees with reviewed result shape')
                result_extents.append(text)
        for arg in [*dimensions, *args]:
            doc = config["argument_docs"].get(arg)
            if not isinstance(doc, str) or not doc or not doc.isascii() or "\n" in doc:
                raise ValueError(f"{name}: argument_docs needs one ASCII line for {arg}")
        review = []
        for arg, expected in shapes.items():
            declaration = declarations[arg]
            start, end = declaration.item.span
            # Read only comments on this AST-identified declaration, never executable code.
            comments = [line.split('!!', 1)[1].strip() for line in source_lines[start-1:end]
                        if '!!' in line]
            match = re.match(r"(?:shape\s*)?\(\s*([a-z]\w*(?:\s*,\s*[a-z]\w*)*)\s*\)(?:\s|[;:]|$)",
                             ' '.join(comments), re.I)
            documented = [v.strip().lower() for v in match[1].split(',')] if match else None
            if documented != expected:
                review.append(dict(argument=arg, required_shape=expected, documented_shape=documented,
                                   reason='missing or unsupported trailing shape comment' if documented is None
                                   else 'documented shape differs from wrapper policy'))
        entry = dict(name=name, policy=policy, args=args, metadata=metadata,
                     dimensions=dimensions, documentation_review=review)
        if codecs:
            entry['state_codecs'] = codecs
        if synthetic_status:
            entry.update(native_args=native_args, synthetic_status=True)
        if optional_arrays:
            entry.update(native_args=native_args, optional_arrays=optional_arrays)
        if is_function:
            entry.update(native_args=native_args, result_argument=result)
            if result_extents:
                entry['result_extents'] = result_extents
        if choices:
            entry['length_choices'] = choices
        if 'shape_contract' in policy:
            if type(policy['shape_contract']) is not bool or not policy['shape_contract']:
                raise ValueError(f'{name}: shape_contract must be true when specified')
            entry['shape_contract'] = shape_contracts.prepare(node, source_lines, entry)
        entries.append(entry)
    return entries


def module_imports(entries, reserved=()):
    """Share codec imports through host association, separating distinct bindings."""
    bindings, aliases, modules = {}, set(reserved), {}
    for entry in entries:
        entry['state_import_aliases'] = {}
        for arg, codec in entry.get('state_codecs', {}).items():
            names = {}
            for role, module, remote in (('type', codec['type_module'], codec['type']),
                                         ('load', codec['module'], codec['load']),
                                         ('store', codec['module'], codec['store'])):
                key = (role, module, remote)
                if key not in bindings:
                    base = f'cg_{role}_{arg}'
                    alias, number = base, 1
                    while alias in aliases:
                        number += 1
                        suffix = f'_{number}'
                        alias = base[:63-len(suffix)] + suffix
                    aliases.add(alias)
                    bindings[key] = alias
                    modules.setdefault(module, []).append(f'{alias} => {remote}')
                names[role] = bindings[key]
            entry['state_import_aliases'][arg] = names
    lines = []
    for module, names in modules.items():
        line = f'    use {module}, only : '
        for i, name in enumerate(names):
            term = name + (', ' if i < len(names)-1 else '')
            if len(line) + len(term) > 100 and not line.endswith('only : '):
                lines.append(line.rstrip() + ' &')
                line = '        '
            line += term
        lines.append(line.rstrip())
    return lines


def emit(entry, index, config):
    name, policy, args, metadata = (entry[key] for key in ("name", "policy", "args", "metadata"))
    symbol, status = policy["symbol"], policy["status_argument"]
    dimensions, shapes = entry["dimensions"], policy["shapes"]
    formal = "        " + ", ".join([*dimensions, *args]) + ") &"
    formal_lines = [formal]
    if len(formal) > 132:
        formal_lines, line = [], "        "
        names = [*dimensions, *args]
        for i, arg in enumerate(names):
            term = arg + (", " if i < len(names)-1 else ") &")
            if len(line) + len(term) > 120 and line.strip():
                formal_lines.append(line.rstrip() + " &")
                line = "        "
            line += term
        formal_lines.append(line)
    lines = [f"    subroutine {symbol}( &", *formal_lines,
             f'        bind(c, name="{symbol}")', f"        !! {policy['description']}"]
    state_aliases = entry.get('state_import_aliases', {})
    cargs = []
    for dim in dimensions:
        lines.append(f"        integer(c_int), value, intent(in) :: {dim} !! {config['argument_docs'][dim]}")
        cargs.append(f"int arg_{dim}")
    for arg in args:
        if arg == status:
            declaration = "integer(c_int), intent(out)"
            cargs.append(f"int *arg_{arg}")
        else:
            info = metadata[arg]
            # Preserve output on adapter validation errors; native calls retain native intent.
            intent = "in" if info["intent"] == "in" else "inout"
            value = ", value" if info["rank"] == 0 and intent == "in" else ""
            if info['category'] == 'state':
                declaration = 'integer(c_int64_t), intent(inout)'
                ctype = 'int64_t *'
            elif info['category'] == 'logical':
                declaration = "integer(c_int), value, intent(in)"
                ctype = "int"
            elif info['category'] == 'integer':
                declaration = f"integer(c_int){value}, intent({intent})"
                ctype = 'int' if value else ('const int *' if intent == 'in' else 'int *')
            else:
                declaration = f"real(c_double){value}, intent({intent})"
                ctype = "double" if value else ("const double *" if intent == "in" else "double *")
            cargs.append(f"{ctype} arg_{arg}")
        shape = "(" + ",".join(shapes[arg]) + ")" if arg in shapes else ""
        if arg in entry.get('state_codecs', {}):
            shape = f"({entry['state_codecs'][arg]['length']})"
        doc = config['argument_docs'][arg]
        # Explicit extents already document the shape; retain the argument's meaning.
        match = re.match(r"(?:shape\s*)?(\([^()]+\))\s*;\s*(.+)", doc, re.I) if shape else None
        if match and re.sub(r"\s+", "", match[1]).lower() == shape.lower():
            doc = match[2]
        lines.append(f"        {declaration} :: {arg}{shape} !! {doc}")
    for arg in policy.get('logical_inputs', {}):
        lines.append(f"        logical :: cg_logical_{arg}")
    if entry.get('state_codecs'):
        lines.append('        integer :: cg_state_status')
        for arg in entry['state_codecs']:
            lines.append(f"        type({state_aliases[arg]['type']}) :: cg_state_{arg}")
    for arg in entry.get('optional_arrays', {}):
        shape = ','.join([':']*metadata[arg]['rank'])
        lines.append(f'        real(c_double), allocatable :: cg_optional_{arg}({shape})')
    lines += ["        ! Direct array association requires matching native and C kinds.",
              "        integer, parameter :: cg_real_kind_check=1/merge(1,0,cg_kind == c_double)",
              "        integer, parameter :: cg_int_kind_check=1/merge(1,0,kind(0) == c_int)", "",
              "        ! Reject invalid extents before calling the native procedure."]
    for dim, minimum in dimensions.items():
        lines += [f"        if ({dim} < {minimum}_c_int) then",
                  f"            {status}={policy['invalid_status']}_c_int", "            return", "        end if"]
    for dim, rule in entry.get('length_choices', {}).items():
        terms = [f"{dim} /= {rule['dimension']}" + (f" - {-offset}_c_int" if offset else '')
                 for offset in rule['offsets']]
        lines += ['        if ( &']
        lines += ['            ' + term + (' .and. &' if i < len(terms)-1 else ') then')
                  for i, term in enumerate(terms)]
        lines += [f"            {status}={policy['invalid_status']}_c_int", '            return', '        end if']
    for arg, mapping in policy.get('logical_inputs', {}).items():
        lines += [f"        if ({arg} /= {mapping['false']}_c_int .and. &",
                  f"                {arg} /= {mapping['true']}_c_int) then",
                  f"            {status}={policy['invalid_status']}_c_int", "            return", "        end if",
                  f"        cg_logical_{arg}=({arg} == {mapping['true']}_c_int)"]
    for arg, rule in policy.get('integer_limits', {}).items():
        array = metadata[arg]['rank'] > 0
        checks = [f"{arg} < {rule['minimum']}_c_int64_t", f"{arg} > {rule['maximum']}_c_int64_t"]
        if array:
            checks = [f'any({check})' for check in checks]
        lines += [f'        if ({checks[0]} .or. &', f'                {checks[1]}) then',
                  f"            {status}={policy['invalid_status']}_c_int", '            return', '        end if']
    for arg, ceiling in policy.get('integer_sums', {}).items():
        lines += [f'        if (sum(abs(int({arg},c_int64_t))) > {ceiling}_c_int64_t) then',
                  f"            {status}={policy['invalid_status']}_c_int", '            return', '        end if']
    for flag in entry.get('optional_arrays', {}).values():
        lines += [f'        if ({flag} /= 0_c_int .and. {flag} /= 1_c_int) then',
                  f"            {status}={policy['invalid_status']}_c_int", '            return', '        end if']
    for arg, codec in entry.get('state_codecs', {}).items():
        lines += [f"        call {state_aliases[arg]['load']}({arg},cg_state_{arg},cg_state_status)",
                  '        if (cg_state_status /= 0) then',
                  f"            {status}={codec['invalid_status']}_c_int", '            return', '        end if']
    if entry.get('optional_arrays'):
        lines += ['        ! An unallocated allocatable actual makes an optional array dummy absent.']
    for arg, flag in entry.get('optional_arrays', {}).items():
        lines += [f'        if ({flag} == 1_c_int) then',
                  f'            allocate(cg_optional_{arg}, source={arg})', '        end if']
    native_args = entry.get('native_args', args)
    call = (f"{entry['result_argument']}=cg_native_{index}" if 'result_argument' in entry
            else f'call cg_native_{index}')
    lines += ["        ! The compiler supplies assumed-shape descriptors inside Fortran.",
              f"        {call}( &"]
    for i, arg in enumerate(native_args):
        actual = f"cg_logical_{arg}" if arg in policy.get('logical_inputs', {}) else arg
        if arg in entry.get('optional_arrays', {}):
            actual = f'cg_optional_{arg}'
        if arg in entry.get('state_codecs', {}):
            actual = f'cg_state_{arg}'
        lines.append(f"            {arg}={actual}" + (", &" if i < len(native_args)-1 else ")"))
    for arg in entry.get('state_codecs', {}):
        lines.append(f"        call {state_aliases[arg]['store']}(cg_state_{arg},{arg})")
    if 'result_argument' in entry or entry.get('synthetic_status'):
        lines.append(f'        {status}=0_c_int')
    lines += [f"    end subroutine {symbol}", ""]
    header = [f"/* {policy['description']}",
              " * Explicit-shape arrays: contiguous column-major storage; no overlap.",
              " * All array pointers must be valid and non-NULL, including empty arrays.",
              " * All scalar output pointers must also be valid and non-NULL.",
              f" * Invalid extents return status {policy['invalid_status']} without changing output.",
              " * Other status/output behavior is inherited from the native procedure. */",
              f"void {symbol}(" + ", ".join(cargs) + ");", ""]
    if policy.get('logical_inputs'):
        header[-3:-3] = [f" * {arg}: false={mapping['false']}, true={mapping['true']}; required, no default."
                         for arg, mapping in policy['logical_inputs'].items()]
        header[-3:-3] = [" * Other flag values return invalid status with all outputs unchanged."]
    record = dict(procedure=name, state="generated", symbols=[symbol], interface="explicit_shape_subroutine",
                  dimensions=dimensions, dimension_order=list(dimensions), shapes=shapes, status_argument=status,
                  optional_promoted_to_required=[a for a in args if metadata.get(a, {}).get('optional')
                                                and a not in entry.get('optional_arrays', {})],
                  arguments=args,
                  native_shape_documentation_review=entry["documentation_review"])
    if entry.get('optional_arrays'):
        record['optional_arrays'] = entry['optional_arrays']
        header[-3:-3] = [' * Optional-array flags: 0 absent, 1 present; other values preserve outputs and return invalid status.',
                        ' * Buffers remain required and fully sized when absent, but are not read.']
    if entry.get('state_codecs'):
        record['state_codecs'] = entry['state_codecs']
        for arg, codec in entry['state_codecs'].items():
            header[-3:-3] = [f" * {arg}: caller-owned int64_t[{codec['length']}]; invalid state returns {codec['invalid_status']} unchanged.",
                            ' * State is validated even for empty calls and exported after native numerical errors.']
    if policy.get('optional_status'):
        record['optional_status_promoted_to_required'] = True
    if policy.get('logical_inputs'):
        record['logical_inputs'] = policy['logical_inputs']
    if entry.get('shape_contract'):
        record['shape_contract'] = entry['shape_contract']
    if entry.get('length_choices'):
        record['length_choices'] = entry['length_choices']
        header[-3:-3] = [f" * Required: {dim} in ({shape_contracts.choice_text(rule)}); invalid status preserves output."
                         for dim, rule in entry['length_choices'].items()]
    for key in ('integer_limits', 'integer_sums'):
        if policy.get(key):
            record[key] = policy[key]
    if policy.get('integer_limits') or policy.get('integer_sums'):
        header[-3:-3] = [f" * Required integer range: {r['minimum']} <= {a} <= {r['maximum']}."
                         for a, r in policy.get('integer_limits', {}).items()]
        header[-3:-3] = [f" * Required widened bound: sum(abs({a})) <= {ceiling}."
                         for a, ceiling in policy.get('integer_sums', {}).items()]
        header[-3:-3] = [f" * Rejected integer guards return status {policy['invalid_status']} with output unchanged."]
    if 'result_argument' in entry:
        record.update(interface='explicit_shape_function', result_argument=entry['result_argument'])
        if entry.get('result_extents'):
            record['native_result_extents'] = entry['result_extents']
    if entry.get('synthetic_status'):
        record['synthetic_status'] = True
    if 'result_argument' in entry or entry.get('synthetic_status'):
        header[-3:-3] = [' * Status 0 means the call was valid; native numerical NaN/infinity is retained.']
        header = [line.replace(' * Other status/output behavior is inherited from the native procedure. */',
                               ' * Numerical parameter errors use native values, not call-level statuses. */')
                  for line in header]
    return lines, header, record
