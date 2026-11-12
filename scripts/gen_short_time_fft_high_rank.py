#!/usr/bin/env python3
"""Generate rank-6..14 ShortTimeFFT overloads and result types.

The SciPy API replaces a time axis with a frequency axis, and appends a time
slice axis. Fortran's explicit-result-rank interfaces therefore require a
small rank-specific adapter, while all transforms use one shared 1-D kernel.
With maximum Fortran array rank 15, native inputs can have ranks at most 14.
Run from the repository root: python scripts/gen_short_time_fft_high_rank.py
This script refuses to regenerate an already extended module.
"""
from pathlib import Path
import re

p = Path(__file__).resolve().parents[1] / 'src' / 'signal_short_time_fft.f90'
s = p.read_text()
ranks = tuple(range(6, 15))
if 'type, public :: short_time_fft_result_7d' in s:
    raise SystemExit('ShortTimeFFT rank 6..14 support already generated')

def subs(n): return ','.join(f'idx({i})' for i in range(1,n+1))
def osubs(n): return ','.join(f'oidx({i})' for i in range(1,n+1))
def slc(n): return ','.join(':' for _ in range(n))
def shape(n,var='dout'): return ','.join(f'{var}({i})' for i in range(1,n+1))

newtypes = ''
for r in ranks:
    out = r+1
    newtypes += f'''    ! Rank-{r} signals yield rank-{out} STFT coefficients.
    type, public :: short_time_fft_result_{out}d
        complex(dp), allocatable :: values({slc(out)}) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_{out}d

    type, public :: short_time_spectrogram_result_{out}d
        real(dp), allocatable :: power({slc(out)}) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_{out}d

'''
s=s.replace('    type, public :: short_time_fft\n',newtypes+'    type, public :: short_time_fft\n',1)

names = ('stft', 'stft_detrend', 'istft', 'spectrogram', 'cross_spectrogram')
for method in names:
    methods = []
    for r in ranks:
        n = r+1 if method == 'istft' else r
        methods.append(f'        procedure, private :: {method}_nd{n}_real, {method}_nd{n}_complex\n')
    pattern = re.compile(r'(?m)^        generic :: ' + method + r' =>[^\n]*(?:\n            [^\n]*)*')
    old = pattern.search(s)
    if old is None:
        raise RuntimeError('missing type-bound generic '+method)
    parts = [f'{method}_real',f'{method}_complex']
    for r in range(2, 15):
        n = r+1 if method == 'istft' else r
        parts.extend([f'{method}_nd{n}_real',f'{method}_nd{n}_complex'])
    generic = '        generic :: '+method+' => '+', '.join(parts[:2])+', &\n'
    for i in range(2,len(parts),2):
        generic += '            '+', '.join(parts[i:i+2])
        if i+2 < len(parts): generic += ', &'
        generic += '\n'
    replacement = ''.join(methods)+generic.rstrip('\n')
    s = s[:old.start()] + replacement + s[old.end():]

code = ''
for r in ranks:
    out=r+1
    axes='['+','.join(map(str,range(1,r+1)))+']'
    axes2='['+','.join(map(str,range(1,out+1)))+']'
    for kind in ('real','complex'):
        xt='real(dp)' if kind=='real' else 'complex(dp)'
        arg='cmplx(x,0.0_dp,dp)' if kind=='real' else 'x'
        is_complex = '.false.' if kind=='real' else '.true.'
        for detrending in (False,True):
            method='stft_detrend' if detrending else 'stft'
            name=f'{method}_nd{r}_{kind}'
            params='self, x, result, '+('detrend, ' if detrending else '')+'p0, p1, k_offset, padding, axis'
            code+=f'''    subroutine {name}({params})
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        {xt}, intent(in) :: x({slc(r)}) !! input time series
        type(short_time_fft_result_{out}d), intent(out) :: result !! spectrum
'''
            if detrending: code+='        character(*), intent(in) :: detrend !! constant or linear segment detrending\n'
            code+='''        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
'''
            code+=f'''        call stft_nd{r}_kernel(self, {arg}, result, axis, p0, p1, &
            k_offset, padding, {is_complex}'''+(', detrend)' if detrending else ')')+f'''
    end subroutine {name}

'''
    code+=f'''    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd{r}_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x({slc(r)}) !! input signals
        type(short_time_fft_result_{out}d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims({r}), dout({out}), idx({r}), oidx({out})
        integer :: batch_axes({r-1}), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,{r},{r})
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack({axes},{axes} /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,{r-1}
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x({subs(r)})
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:{r}) = dims
                dout(a) = size(frame%values,1)
                dout({out}) = size(frame%values,2)
                allocate(result%values({shape(out)}))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:{r}) = idx
                    oidx(a) = q
                    oidx({out}) = t
                    result%values({osubs(out)}) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd{r}_kernel

'''
    for kind in ('real','complex'):
        xt='real(dp)' if kind=='real' else 'complex(dp)'
        code+=f'''    subroutine istft_nd{out}_{kind}(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s({slc(out)}) !! STFT coefficients
        {xt}, allocatable, intent(out) :: x({slc(r)}) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z({slc(r)})

        allocate(x({','.join(['0']*r)}))
        status = short_fft_invalid_input
'''
        if kind=='real':
            code+='        if (self%fft_mode /= \'onesided\' .and. self%fft_mode /= \'onesided2X\') return\n'
        else:
            code+='        if (self%fft_mode == \'onesided\' .or. self%fft_mode == \'onesided2X\') return\n'
        code+=f'''        call inverse_nd{out}_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = '''+ ('real(z,dp)' if kind=='real' else 'z')+f'''
    end subroutine istft_nd{out}_{kind}

'''
    code+=f'''    subroutine inverse_nd{out}_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s({slc(out)}) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x({slc(r)}) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims({out}), dout({r}), batch_axes({r-1}), out_axes({r})
        integer :: idx({out}), oidx({r}), b, ib, j, q, p, k, nb

        allocate(x({','.join(['0']*r)}))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,{out},{r})
        t = fft_resolve_axis(t_axis,{out},{out})
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack({axes2},{axes2} /= f .and. {axes2} /= t)
        out_axes = pack({axes2},{axes2} /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,{r-1}
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s({subs(out)})
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,{axes2} /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x({shape(r)}))
            end if
            do j=1,{r-1}
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x({osubs(r)}) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd{out}_kernel

'''
    for kind in ('real','complex'):
        xt='real(dp)' if kind=='real' else 'complex(dp)'
        code+=f'''    subroutine spectrogram_nd{r}_{kind}(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        {xt}, intent(in) :: x({slc(r)}) !! signal array
        type(short_time_spectrogram_result_{out}d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_{out}d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd{r}_{kind}

'''
    for kind in ('real','complex'):
        xt='real(dp)' if kind=='real' else 'complex(dp)'
        code+=f'''    subroutine cross_spectrogram_nd{r}_{kind}(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        {xt}, intent(in) :: x({slc(r)}) !! first array of signals
        {xt}, intent(in) :: y({slc(r)}) !! second array of signals
        complex(dp), allocatable, intent(out) :: values({slc(out)}) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_{out}d) :: sx, sy

        allocate(values({','.join(['0']*out)}))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd{r}_{kind}

'''
s=s.replace('end module scifort_signal_short_time_fft',code+'end module scifort_signal_short_time_fft',1)
# Keep standard free-form line lengths without compiler-specific options.
wrapped = []
for line in s.splitlines():
    while len(line) > 132:
        pos = line.rfind(',', 0, 109)
        if pos < 20:
            raise RuntimeError('Cannot wrap generated Fortran line: '+line)
        wrapped.append(line[:pos+1]+' &')
        line = '            '+line[pos+1:].lstrip()
    wrapped.append(line)
s = '\n'.join(wrapped)+'\n'
p.write_text(s)
print('Generated ShortTimeFFT ranks 6..14; total lines:',len(s.splitlines()))
