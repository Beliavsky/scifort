#!/usr/bin/env python3
"""Generate rank-2/3 axis-aware wrappers for existing 1-D multirate kernels.

The generated wrappers preserve the original 1-D interfaces and status semantics.
Positive axes are one-based (Fortran), and negative axes count from the end.
Run from any directory; repeated invocations are idempotent.
"""
from pathlib import Path

src = Path(__file__).resolve().parents[1] / 'src'
begin = '    ! BEGIN GENERATED MULTIRATE ND WRAPPERS'
end = '    ! END GENERATED MULTIRATE ND WRAPPERS'


def generate(name, kind, rank, xt, htype=None):
    dims = ','.join(':' for _ in range(rank))
    qname = f'{name}_nd_{kind}{rank}'
    if name == 'resample':
        result_t = xt
        args = 'x, num, window, t, output_times, ' + ('domain, ' if kind == 'c' else '') + 'axis, status'
        declarations = [
            f'        {xt}, intent(in) :: x({dims}) !! input samples',
            '        integer, intent(in) :: num !! output length along the transform axis',
            '        character(*), intent(in), optional :: window !! periodic spectral window',
            '        real(dp), intent(in), optional :: t(:) !! original sample times',
            '        real(dp), allocatable, intent(out), optional :: output_times(:) !! resampled times',
        ]
        if kind == 'c':
            declarations.append('        character(*), intent(in), optional :: domain !! time or freq')
        firstargs = 'num, window, t, output_times, ' + ('domain, ' if kind == 'c' else '')
        otherargs = 'num, window=window, ' + ('domain=domain, ' if kind == 'c' else '')
        default_axis = 1
        status_constant = 'fourier_resample_status_invalid_input'
        success_constant = 'fourier_resample_status_success'
    elif name == 'upfirdn':
        result_t = 'real(dp)' if kind == 'rr' else 'complex(dp)'
        args = 'h, x, up, down, mode, cval, axis, status'
        declarations = [
            f'        {htype}, intent(in) :: h(:) !! FIR coefficients',
            f'        {xt}, intent(in) :: x({dims}) !! input samples',
            '        integer, intent(in), optional :: up !! integer upsampling factor',
            '        integer, intent(in), optional :: down !! integer downsampling factor',
            '        character(*), intent(in), optional :: mode !! extension mode',
            f'        {result_t}, intent(in), optional :: cval !! constant padding value',
        ]
        firstargs = 'up, down, mode, cval, '
        otherargs = firstargs
        default_axis = rank
        status_constant = 'resample_status_invalid_input'
        success_constant = 'resample_status_success'
    elif name in ('resample_poly', 'resample_poly_ctaps'):
        is_complex_taps = name == 'resample_poly_ctaps'
        result_t = 'complex(dp)' if is_complex_taps else xt
        args = 'x, up, down, taps, window, beta, padtype, cval, axis, status'
        declarations = [
            f'        {xt}, intent(in) :: x({dims}) !! input samples',
            '        integer, intent(in) :: up !! integer upsampling factor',
            '        integer, intent(in) :: down !! integer downsampling factor',
            ('        complex(dp), intent(in) :: taps(:) !! required complex FIR coefficients'
             if is_complex_taps else '        real(dp), intent(in), optional :: taps(:) !! custom real FIR coefficients'),
            '        character(*), intent(in), optional :: window !! FIR window',
            '        real(dp), intent(in), optional :: beta !! Kaiser parameter',
            '        character(*), intent(in), optional :: padtype !! boundary extension',
            f'        {result_t}, intent(in), optional :: cval !! constant padding value',
        ]
        firstargs = 'up, down, taps, window, beta, padtype, cval, '
        otherargs = firstargs
        default_axis = 1
        status_constant = 'resample_status_invalid_input'
        success_constant = 'resample_status_success'
    elif name == 'decimate':
        result_t = xt
        args = 'x, q, n, ftype, zero_phase, axis, status, taps, sos, b, a'
        declarations = [
            f'        {xt}, intent(in) :: x({dims}) !! input samples',
            '        integer, intent(in) :: q !! decimation factor',
            '        integer, intent(in), optional :: n !! filter order',
            '        character(*), intent(in), optional :: ftype !! fir or iir',
            '        logical, intent(in), optional :: zero_phase !! noncausal if true',
        ]
        firstargs = 'q, n, ftype, zero_phase, '
        otherargs = firstargs
        default_axis = rank
        status_constant = 'resample_status_invalid_input'
        success_constant = 'resample_status_success'
    else:
        raise AssertionError(name)

    declarations += [
        '        integer, intent(in), optional :: axis !! one-based or negative axis',
        '        integer, intent(out), optional :: status !! zero on success',
    ]
    axis_var = 'ax' if name == 'decimate' else 'a'
    if name == 'decimate':
        dtype = 'real(dp)' if kind == 'r' else 'complex(dp)'
        declarations += [
            '        real(dp), intent(in), optional :: taps(:) !! supplied real FIR coefficients',
            '        real(dp), intent(in), optional :: sos(:,:) !! supplied real IIR SOS coefficients',
            f'        {dtype}, intent(in), optional :: b(:) !! transfer-function numerator',
            f'        {dtype}, intent(in), optional :: a(:) !! transfer-function denominator',
        ]
    declarations += [
        f'        {result_t}, allocatable :: y({dims}), line(:)',
        f'        integer :: {axis_var}, code, sizes({rank})',
        '        integer :: i, j, k' if rank == 3 else '        integer :: i, j',
    ]
    srcname = (f'{name}_{kind}' if name == 'upfirdn' else
               f'resample_poly_{"real" if kind=="r" else "complex"}_taps'
               if name == 'resample_poly_ctaps' else
               f'{name}_{"real" if kind=="r" else "complex"}')
    extra = ', taps=taps, sos=sos, b=b, a=a' if name == 'decimate' else ''
    def first_call(s):
        if name == 'upfirdn':
            return f'{srcname}(h, {s}, {firstargs}status=code)'
        return f'{srcname}({s}, {firstargs}status=code{extra})'
    def other_call(s):
        if name == 'upfirdn':
            return f'{srcname}(h, {s}, {otherargs}status=code)'
        return f'{srcname}({s}, {otherargs}status=code{extra})'
    # For the first slice, all other axes are set to 1.
    sfirst = {2:['x(:,1)', 'x(1,:)'], 3:['x(:,1,1)', 'x(1,:,1)', 'x(1,1,:)']}[rank]
    call_cases = []
    for a in range(1,rank+1):
        fixed = [p for p in range(1,rank+1) if p!=a]
        sym = {fixed[0]:'i'} if len(fixed)==1 else {fixed[0]:'i',fixed[1]:'j'}
        input_slice = ','.join(':' if p==a else sym[p] for p in range(1,rank+1))
        input_slice = f'x({input_slice})'
        output_slice = f'y({",".join(":" if p==a else sym[p] for p in range(1,rank+1))})'
        body = [f'                line = {other_call(input_slice)}',
                f'                if (code /= {success_constant} .or. size(line) /= sizes({axis_var})) then',
                '                    deallocate(y)',
                '                    allocate(y(' + ','.join('0' for _ in range(rank)) + '))',
                '                    return',
                '                end if',
                f'                {output_slice} = line']
        for p in reversed(fixed):
            body = [f'            do {sym[p]} = 1, size(x,{p})'] + ['    '+line for line in body] + ['            end do']
        call_cases.extend([f'        case ({a})'] + body)
    lines = [f'    function {qname}({args}) result(y)',*declarations,
             f'        allocate(y({",".join("0" for _ in range(rank))}))',
             f'        if (present(status)) status = {status_constant}',
             '        if (any(shape(x) == 0)) return',
             f'        {axis_var} = {default_axis}',
             f'        if (present(axis)) {axis_var} = axis',
             f'        if ({axis_var} < 0) {axis_var} = {axis_var} + {rank} + 1',
             f'        if ({axis_var} < 1 .or. {axis_var} > {rank}) return',
             f'        select case ({axis_var})']
    for idx, slic in enumerate(sfirst,1):
        lines += [f'        case ({idx})', f'            line = {first_call(slic)}']
    lines += ['        end select',
              f'        if (code /= {success_constant} .or. size(line) < 1) return',
              '        sizes = shape(x)',
              f'        sizes({axis_var}) = size(line)',
              '        deallocate(y)',
              '        allocate(y(' + ','.join('sizes('+str(p)+')' for p in range(1,rank+1)) + '))',
              f'        select case ({axis_var})',*call_cases,
              '        end select',
              f'        if (present(status)) status = {success_constant}',
              f'    end function {qname}', '']
    return '\n'.join(lines)


def patch_file(filename, routines):
    """Update decimate wrappers without erasing later API extensions.

    Other existing wrappers have acquired parameterized windows and callable
    Fourier-window interfaces since the original generator was written;
    retain their complete declarations and implementations. A missing routine
    can still be generated from this template.
    """
    path = src/filename
    text = path.read_text()
    if begin not in text or end not in text:
        raise ValueError(f'missing wrapper markers in {path}')
    import re
    for generic, procs in routines.items():
        for proc in procs:
            name, rest = proc.split('_nd_')
            kind, rank = rest[:-1], int(rest[-1])
            xtype = 'real(dp)' if kind in ('rr','cr','r') else 'complex(dp)'
            htype = ('real(dp)' if kind[0]=='r' else 'complex(dp)') if name == 'upfirdn' else None
            match = re.search(r'^    function ' + re.escape(proc) +
                              r'\(.*?^    end function ' + re.escape(proc) + r'$',
                              text, flags=re.MULTILINE | re.DOTALL)
            if match is not None:
                if name == 'decimate':
                    text = text[:match.start()] + generate(name,kind,rank,xtype,htype).rstrip() + text[match.end():]
            else:
                # First creation of a named generated wrapper.
                code = generate(name,kind,rank,xtype,htype)
                text = text.replace(end,code + '\n' + end,1)
                marker = f'    interface {generic}\n'
                if marker not in text:
                    raise ValueError(marker)
                text = text.replace(marker,marker + f'        module procedure {proc}\n',1)
    path.write_text(text)


patch_file('signal_fourier_resample.f90', {
    'resample': [f'resample_nd_{k}{r}' for r in (2,3) for k in 'rc']
})
patch_file('signal_resample.f90', {
    'upfirdn': [f'upfirdn_nd_{k}{r}' for r in (2,3) for k in ('rr','rc','cr','cc')],
    'resample_poly': ([f'resample_poly_nd_{k}{r}' for r in (2,3) for k in 'rc'] +
                      [f'resample_poly_ctaps_nd_{k}{r}' for r in (2,3) for k in 'rc']),
    'decimate': [f'decimate_nd_{k}{r}' for r in (2,3) for k in 'rc'],
})
print('Generated 24 rank-2/3 multidimensional resampling wrappers.')
