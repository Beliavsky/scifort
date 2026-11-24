# SPDX-License-Identifier: MIT
"""Validated explicit-shape adapters for native assumed-shape subroutines."""
import re

import shape_contracts

from fparser.two import Fortran2003 as F


def prepare(config, api, identifier, source_lines):
    entries = []
    for name, policy in sorted(config.get("subroutines", {}).items()):
        identifier(name)
        fields = {"symbol", "dimensions", "shapes", "status_argument", "invalid_status", "description"}
        if not isinstance(policy, dict) or set(policy) - {"logical_inputs", "shape_contract"} != fields:
            raise ValueError(f"{name}: subroutine policy requires {sorted(fields)}")
        identifier(policy["symbol"])
        node = api[name]
        if not isinstance(node, F.Subroutine_Subprogram) or node.content[0].items[3] is not None:
            raise ValueError(f"{name}: expected a directly defined non-bind(c) subroutine")
        args = [identifier(str(a).lower()) for a in node.content[0].items[2].items]
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
        if not isinstance(logical_inputs, dict) or set(logical_inputs) - set(args):
            raise ValueError(f"{name}: logical_inputs must name native arguments")
        for arg, mapping in logical_inputs.items():
            if (not isinstance(mapping, dict) or set(mapping) != {"false", "true"}
                    or any(type(v) is not int or not -2147483647 <= v <= 2147483647
                           for v in mapping.values()) or mapping['false'] == mapping['true']):
                raise ValueError(f"{name}: logical mapping needs distinct C int false/true values")
        metadata = {}
        status = policy["status_argument"]
        if status in logical_inputs:
            raise ValueError(f"{name}: status cannot be a logical mapping")
        if status not in args:
            raise ValueError(f"{name}: status_argument must name a native argument")
        for arg in args:
            if arg not in declared:
                raise ValueError(f"{name}: {arg} needs an explicit declaration")
            dtype, attrs, entity = declared[arg]
            if arg == status:
                if str(dtype).lower() != "integer" or attrs != {"intent(out)"} or any(entity.items[1:]):
                    raise ValueError(f"{name}: status must be a default-integer scalar intent(out)")
                continue
            intent = next((value for value in ("in", "out", "inout") if f"intent({value})" in attrs), None)
            base_attrs = attrs - {"optional"} if arg in logical_inputs else attrs
            if intent is None or base_attrs != {f"intent({intent})"}:
                raise ValueError(f"{name}: {arg} requires intent and no unconfigured optional attributes")
            if str(dtype).lower() == "logical":
                if (arg not in logical_inputs or intent != "in" or any(entity.items[1:])):
                    raise ValueError(f"{name}: logical inputs need an explicit scalar mapping")
                metadata[arg] = dict(intent=intent, rank=0, category='logical',
                                     optional='optional' in attrs)
                continue
            if arg in logical_inputs:
                raise ValueError(f"{name}: logical_inputs must map logical arguments")
            if str(dtype).lower() == "integer":
                if intent != "out" or any(entity.items[1:]):
                    raise ValueError(f"{name}: only default-integer scalar outputs are supported")
                metadata[arg] = dict(intent=intent, rank=0, category='integer')
                continue
            kind = config["kind_name"]
            if str(dtype).lower().replace(" ", "") not in (f"real({kind})", f"real(kind={kind})"):
                raise ValueError(f"{name}: {arg} must be real({kind})")
            shape = entity.items[1]
            if any(entity.items[2:]) or (shape and (not isinstance(shape, F.Assumed_Shape_Spec_List)
                                                  or any(str(item) != ":" for item in shape.items))):
                raise ValueError(f"{name}: {arg} must be scalar or assumed-shape with ':' bounds")
            rank = len(shape.items) if shape else 0
            if rank > 2:
                raise ValueError(f"{name}: only rank-one and rank-two arrays are supported")
            metadata[arg] = dict(intent=intent, rank=rank, category='real')
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
        if 'shape_contract' in policy:
            if type(policy['shape_contract']) is not bool or not policy['shape_contract']:
                raise ValueError(f'{name}: shape_contract must be true when specified')
            entry['shape_contract'] = shape_contracts.prepare(node, source_lines, entry)
        entries.append(entry)
    return entries


def emit(entry, index, config):
    name, policy, args, metadata = (entry[key] for key in ("name", "policy", "args", "metadata"))
    symbol, status = policy["symbol"], policy["status_argument"]
    dimensions, shapes = entry["dimensions"], policy["shapes"]
    lines = [f"    subroutine {symbol}( &", "        " + ", ".join([*dimensions, *args]) + ") &",
             f'        bind(c, name="{symbol}")', f"        !! {policy['description']}",
             "        implicit none"]
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
            if info['category'] == 'logical':
                declaration = "integer(c_int), value, intent(in)"
                ctype = "int"
            elif info['category'] == 'integer':
                declaration = "integer(c_int), intent(inout)"
                ctype = "int *"
            else:
                declaration = f"real(c_double){value}, intent({intent})"
                ctype = "double" if value else ("const double *" if intent == "in" else "double *")
            cargs.append(f"{ctype} arg_{arg}")
        shape = "(" + ",".join(shapes[arg]) + ")" if arg in shapes else ""
        lines.append(f"        {declaration} :: {arg}{shape} !! {config['argument_docs'][arg]}")
    for arg in policy.get('logical_inputs', {}):
        lines.append(f"        logical :: cg_logical_{arg}")
    lines += ["        ! Direct array association requires matching native and C kinds.",
              "        integer, parameter :: cg_real_kind_check=1/merge(1,0,cg_kind == c_double)",
              "        integer, parameter :: cg_int_kind_check=1/merge(1,0,kind(0) == c_int)", "",
              "        ! Reject invalid extents before calling the native procedure."]
    for dim, minimum in dimensions.items():
        lines += [f"        if ({dim} < {minimum}_c_int) then",
                  f"            {status}={policy['invalid_status']}_c_int", "            return", "        end if"]
    for arg, mapping in policy.get('logical_inputs', {}).items():
        lines += [f"        if ({arg} /= {mapping['false']}_c_int .and. &",
                  f"                {arg} /= {mapping['true']}_c_int) then",
                  f"            {status}={policy['invalid_status']}_c_int", "            return", "        end if",
                  f"        cg_logical_{arg}=({arg} == {mapping['true']}_c_int)"]
    lines += ["        ! The compiler supplies assumed-shape descriptors inside Fortran.",
              f"        call cg_native_{index}( &"]
    for i, arg in enumerate(args):
        actual = f"cg_logical_{arg}" if arg in policy.get('logical_inputs', {}) else arg
        lines.append(f"            {arg}={actual}" + (", &" if i < len(args)-1 else ")"))
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
                  optional_promoted_to_required=[a for a in args if metadata.get(a, {}).get('optional')],
                  arguments=args,
                  native_shape_documentation_review=entry["documentation_review"])
    if policy.get('logical_inputs'):
        record['logical_inputs'] = policy['logical_inputs']
    if entry.get('shape_contract'):
        record['shape_contract'] = entry['shape_contract']
    return lines, header, record
