#!/usr/bin/env python3
"""Generate only the repeated SELECT RANK dispatch for spectral array I/O.

The transform and channel traversal stay handwritten in signal_spectral_nd.f90.
"""
from pathlib import Path
import re
path = Path(__file__).resolve().parents[1] / 'src' / 'signal_spectral_nd.f90'
source = path.read_text()
for kind in ('real','complex'):
    for action in ('flatten','scatter'):
        start = f'        ! BEGIN GENERATED {action.upper()} {kind.upper()}'
        end = f'        ! END GENERATED {action.upper()} {kind.upper()}'
        pieces = []
        for rank in range(1,16):
            if action == 'flatten':
                expr = 'cmplx(reshape(x,[size(x)]),0.0_dp,dp)' if kind == 'real' else 'reshape(x,[size(x)])'
                statement = f'            flat = {expr}'
            else:
                expr = 'reshape(real(flat,dp),shape(y))' if kind == 'real' else 'reshape(flat,shape(y))'
                statement = f'            y = {expr}'
            pieces += [f'        rank ({rank})', statement]
        marker = start + '\n' + '\n'.join(pieces) + '\n' + end
        source, changes = re.subn(re.escape(start) + r'.*?' + re.escape(end), lambda _:marker, source, flags=re.S)
        if changes != 1: raise RuntimeError(f'bad generator marker {action} {kind}: {changes}')
path.write_text(source)

# Public signatures are mechanically repeated for real/complex combinations.
wrappers = []
for opname, method in [('periodogram',1),('welch',2),('csd',3),('coherence',4)]:
    pairs = [('real',None),('complex',None)] if method <= 2 else [
        ('real','real'),('complex','complex'),('real','complex'),('complex','real')]
    for xtype,ytype in pairs:
        key = {'real':'r','complex':'c'}
        sub = opname + '_into_' + (xtype if ytype is None else key[xtype]+key[ytype])
        first = ('real(dp)' if xtype == 'real' else 'complex(dp)')
        second = ('real(dp)' if ytype == 'real' else 'complex(dp)')
        output = 'complex(dp)' if method == 3 else 'real(dp)'
        vars = ['x'] + (['y'] if method >= 3 else [])
        kwargs = ['fs','window']
        if method >= 2: kwargs += ['nperseg','noverlap']
        kwargs += ['nfft','detrend']
        if method != 4: kwargs += ['scaling']
        if method in (2,3): kwargs += ['average']
        # Append new options after axis/status: all old positional calls stay valid.
        window_options = ['window_parameter', 'window_parameter2',
                          'window_coefficients', 'window_values']
        if method <= 3: window_options += ['return_onesided']
        signature = vars + ['result','frequencies'] + kwargs + ['axis','status']
        n = f'    subroutine {sub}(' + ','.join(signature) + ', &\n'
        n += '        ' + ','.join(window_options) + ')\n'
        n += f'        {first}, intent(in) :: x(..) !! first input, rank 1-15\n'
        if method >= 3:
            n += f'        {second}, intent(in) :: y(..) !! second input, rank 1-15\n'
        n += f'        {output}, intent(inout) :: result(..) !! caller-allocated spectrum\n'
        n += '        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies\n'
        for field in kwargs:
            t = 'real(dp)' if field=='fs' else ('integer' if field in ('nperseg','noverlap','nfft') else 'character(*)')
            n+= f'        {t}, intent(in), optional :: {field} !! SciPy spectral option\n'
        n+='        integer, intent(in), optional :: axis !! one-based or negative; default last\n'
        n+='        integer, intent(out), optional :: status !! zero on success\n'
        n+='        real(dp), intent(in), optional :: window_parameter !! named window shape parameter\n'
        n+='        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter\n'
        n+='        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients\n'
        n+='        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples\n'
        if method <= 3:
            n+='        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies\n' 
        n+='        type(spectral_job) :: job\n'
        n+='        complex(dp), allocatable :: u(:), v(:)\n'
        n+='        integer :: st\n'
        n+=f'        job%method = {method}\n'
        n+=f"        job%first_real = {'.true.' if xtype=='real' else '.false.'}\n"
        n+=f"        job%second_real = {'.true.' if ytype in (None,'real') else '.false.'}\n"
        for field in kwargs: n+=f'        if (present({field})) job%{field} = {field}\n'
        for field in window_options:
            n += f'        if (present({field})) job%{field} = {field}\n'
        n+=f'        call flatten_{xtype}(x,u)\n'
        if method>=3: n+=f'        call flatten_{ytype}(y,v)\n'
        else: n+='        allocate(v(0))\n'
        dims2 = 'shape(y)' if method>=3 else 'shape(x)'
        to = 'complex' if method==3 else 'real'
        n+=f'        call apply_to_{to}(shape(x),{dims2},u,v,result,frequencies,job,axis,st)\n'
        n+='        if (present(status)) status = st\n'
        n+=f'    end subroutine {sub}\n'
        wrappers.append(n)
start = '    ! BEGIN GENERATED PUBLIC WRAPPERS'
end = '    ! END GENERATED PUBLIC WRAPPERS'
source, changes = re.subn(re.escape(start)+r'.*?'+re.escape(end),
    lambda _: start+'\n'+'\n'.join(wrappers)+'\n'+end, source, flags=re.S)
if changes != 1: raise RuntimeError('public wrapper generation markers missing')
path.write_text(source)
