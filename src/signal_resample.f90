! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! One-dimensional polyphase/Fourier resampling and FIR/IIR decimation.
! Independent implementation of the standard multirate DSP formulas.
module scifort_signal_resample
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window, spectral_named_window, spectral_status_success
    use scifort_signal_design, only : firwin, cheby1, butter_result, design_status_success
    use scifort_signal_sos, only : sosfilt, sos_status_success
    use scifort_signal_lfilter, only : lfilter, lfilter_status_success
    use scifort_signal_filtfilt, only : filtfilt, sosfiltfilt, filtfilt_status_success
    implicit none
    private

    integer, parameter, public :: resample_status_success = 0
    integer, parameter, public :: resample_status_invalid_input = 1
    real(dp), parameter :: pi = acos(-1.0_dp)

    interface upfirdn
        module procedure upfirdn_nd_rr2
        module procedure upfirdn_nd_rc2
        module procedure upfirdn_nd_cr2
        module procedure upfirdn_nd_cc2
        module procedure upfirdn_nd_rr3
        module procedure upfirdn_nd_rc3
        module procedure upfirdn_nd_cr3
        module procedure upfirdn_nd_cc3
        module procedure upfirdn_rr, upfirdn_rc, upfirdn_cr, upfirdn_cc
    end interface upfirdn
    interface resample_poly
        module procedure resample_poly_nd_r2
        module procedure resample_poly_nd_c2
        module procedure resample_poly_nd_r3
        module procedure resample_poly_nd_c3
        module procedure resample_poly_ctaps_nd_r2
        module procedure resample_poly_ctaps_nd_c2
        module procedure resample_poly_ctaps_nd_r3
        module procedure resample_poly_ctaps_nd_c3
        module procedure resample_poly_complex_taps, resample_poly_real_taps
        module procedure resample_poly_real, resample_poly_complex
    end interface resample_poly
    interface decimate
        module procedure decimate_nd_r2
        module procedure decimate_nd_c2
        module procedure decimate_nd_r3
        module procedure decimate_nd_c3
        module procedure decimate_real, decimate_complex
        module procedure decimate_real_complex_tf
        module procedure decimate_nd_rc2, decimate_nd_rc3
        module procedure decimate_complex_sos, decimate_real_complex_sos
        module procedure decimate_nd_complex_sos2, decimate_nd_complex_sos3
        module procedure decimate_nd_real_complex_sos2, decimate_nd_real_complex_sos3
    end interface decimate
    public :: upfirdn, resample_poly, decimate

contains

    pure integer function gcd_positive(a, b) result(g)
        integer, intent(in) :: a !! positive numerator
        integer, intent(in) :: b !! positive denominator
        integer :: p, q, r
        p = a
        q = b
        do while (q /= 0)
            r = mod(p, q)
            p = q
            q = r
        end do
        g = p
    end function gcd_positive

    ! Values outside the signal are extended before insertion of upsampled zeros.
    pure function extended_value(x, k, mode, cval) result(value)
        complex(dp), intent(in) :: x(:) !! original signal
        integer, intent(in) :: k !! zero-based sample index, possibly out of bounds
        character(*), intent(in) :: mode !! padding convention
        complex(dp), intent(in) :: cval !! constant extrapolation value
        complex(dp) :: value
        integer :: n, j, p

        n = size(x)
        if (k >= 0 .and. k < n) then
            value = x(k+1)
            return
        end if
        select case (mode)
        case ('edge')
            value = x(max(1, min(n, k+1)))
        case ('wrap')
            value = x(modulo(k, n)+1)
        case ('reflect')
            if (n == 1) then
                value = x(1)
            else
                p = 2*n-2
                j = modulo(k, p)
                if (j >= n) j = p-j
                value = x(j+1)
            end if
        case ('symmetric')
            p = 2*n
            j = modulo(k, p)
            if (j >= n) j = p-j-1
            value = x(j+1)
        case ('antisymmetric')
            ! Odd extension about boundaries halfway between adjacent samples.
            p = 2*n
            j = modulo(k, p)
            if (j < n) then
                value = x(j+1)
            else
                value = -x(p-j)
            end if
        case ('antireflect')
            ! Odd reflection around endpoint *values*; each full period
            ! adds 2*(last-first). For n==1 use a constant extension.
            if (n == 1) then
                value = x(1)
            else
                p = 2*n-2
                j = modulo(k,p)
                if (j < n) then
                    value = x(j+1)
                else
                    value = 2.0_dp*x(n)-x(p-j+1)
                end if
                value = value + 2.0_dp*real((k-j)/p,dp)*(x(n)-x(1))
            end if
        case ('line')
            if (n == 1) then
                value = x(1)
            else
                if (k < 0) then
                    value = x(1) + real(k,dp)*(x(n)-x(1))/real(n-1,dp)
                else
                    value = x(n) + real(k-(n-1),dp)*(x(n)-x(1))/real(n-1,dp)
                end if
            end if
        case ('smooth')
            if (n == 1) then
                value = x(1)
            else
                if (k < 0) then
                    value = x(1) + real(k,dp)*(x(2)-x(1))
                else
                    value = x(n) + real(k-(n-1),dp)*(x(n)-x(n-1))
                end if
            end if
        case default
            value = cval
        end select
    end function extended_value

    pure logical function mode_valid(mode) result(ok)
        character(*), intent(in) :: mode !! padding selection
        select case (mode)
        case ('constant', 'edge', 'wrap', 'reflect', 'symmetric', 'antisymmetric', &
              'antireflect', 'line', 'smooth')
            ok = .true.
        case default
            ok = .false.
        end select
    end function mode_valid

    function upfirdn_cc(h, x, up, down, mode, cval, status) result(y)
        complex(dp), intent(in) :: h(:) !! FIR coefficients in increasing time order
        complex(dp), intent(in) :: x(:) !! input signal
        integer, intent(in), optional :: up !! integer upsampling factor, default 1
        integer, intent(in), optional :: down !! integer downsampling factor, default 1
        character(*), intent(in), optional :: mode !! boundary mode (see mode_valid)
        complex(dp), intent(in), optional :: cval !! constant extension value, default zero
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        character(:), allocatable :: boundary
        complex(dp) :: outside
        integer :: nu, nd, nx, nh, ny, m, j, k, t

        if (present(status)) status = resample_status_invalid_input
        allocate(y(0))
        nu = 1
        nd = 1
        if (present(up)) nu = up
        if (present(down)) nd = down
        nx = size(x)
        nh = size(h)
        boundary = 'constant'
        if (present(mode)) boundary = trim(mode)
        outside = cmplx(0.0_dp,0.0_dp,dp)
        if (present(cval)) outside = cval
        if (.not. mode_valid(boundary) .or. nu < 1 .or. nd < 1 .or. nx < 1 .or. nh < 1) return
        if (.not. all(ieee_is_finite(real(h,dp))) .or. .not. all(ieee_is_finite(aimag(h)))) return
        if (.not. all(ieee_is_finite(real(x,dp))) .or. .not. all(ieee_is_finite(aimag(x)))) return
        if (.not. ieee_is_finite(real(outside,dp)) .or. .not. ieee_is_finite(aimag(outside))) return
        if ((nx-1) > (huge(nx)-nh)/nu) return
        ny = ((nx-1)*nu+nh-1)/nd+1
        deallocate(y)
        allocate(y(ny))
        y = cmplx(0.0_dp,0.0_dp,dp)
        do m = 0, ny-1
            t = m*nd
            ! Only filter taps congruent with t modulo the upsampling factor contribute.
            k = mod(t,nu)
            do while (k < nh)
                j = (t-k)/nu
                y(m+1) = y(m+1) + h(k+1)*extended_value(x,j,boundary,outside)
                k = k+nu
            end do
        end do
        if (present(status)) status = resample_status_success
    end function upfirdn_cc

    function upfirdn_rr(h, x, up, down, mode, cval, status) result(y)
        real(dp), intent(in) :: h(:) !! real FIR coefficients
        real(dp), intent(in) :: x(:) !! real input signal
        integer, intent(in), optional :: up !! upsampling factor
        integer, intent(in), optional :: down !! downsampling factor
        character(*), intent(in), optional :: mode !! boundary mode
        real(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: z(:)
        complex(dp) :: outer
        outer = cmplx(0.0_dp,0.0_dp,dp)
        if (present(cval)) outer = cmplx(cval,0.0_dp,dp)
        z = upfirdn_cc(cmplx(h,0.0_dp,dp),cmplx(x,0.0_dp,dp), &
            up,down,mode,outer,status)
        y = real(z,dp)
    end function upfirdn_rr

    function upfirdn_rc(h, x, up, down, mode, cval, status) result(y)
        real(dp), intent(in) :: h(:) !! real FIR coefficients
        complex(dp), intent(in) :: x(:) !! complex input signal
        integer, intent(in), optional :: up !! upsampling factor
        integer, intent(in), optional :: down !! downsampling factor
        character(*), intent(in), optional :: mode !! boundary mode
        complex(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        y = upfirdn_cc(cmplx(h,0.0_dp,dp),x,up,down,mode,cval,status)
    end function upfirdn_rc

    function upfirdn_cr(h, x, up, down, mode, cval, status) result(y)
        complex(dp), intent(in) :: h(:) !! complex FIR coefficients
        real(dp), intent(in) :: x(:) !! real input signal
        integer, intent(in), optional :: up !! upsampling factor
        integer, intent(in), optional :: down !! downsampling factor
        character(*), intent(in), optional :: mode !! boundary mode
        complex(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        y = upfirdn_cc(h,cmplx(x,0.0_dp,dp),up,down,mode,cval,status)
    end function upfirdn_cr

    ! Modified Bessel function I0 for the Kaiser lowpass prototype.
    pure real(dp) function i0_series(t) result(value)
        real(dp), intent(in) :: t !! nonnegative argument
        real(dp) :: term
        integer :: j
        value = 1.0_dp
        term = 1.0_dp
        do j = 1, 100
            term = term*(t/(2.0_dp*real(j,dp)))**2
            value = value+term
            if (term <= 8.0_dp*epsilon(value)*value) exit
        end do
    end function i0_series

    function default_taps(rate, window, beta, status, &
            window_parameter, window_parameter2, window_coefficients) result(taps)
        integer, intent(in) :: rate !! largest reduced up/down factor
        character(*), intent(in) :: window !! Kaiser or existing symmetric window
        real(dp), intent(in) :: beta !! Kaiser shape
        integer, intent(out) :: status !! design result status
        real(dp), intent(in), optional :: window_parameter !! first named window parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! cosine-series weights
        real(dp), allocatable :: taps(:), win(:)
        real(dp) :: w, u, r, norm
        integer :: half, i, st
        status = resample_status_invalid_input
        allocate(taps(0))
        if (rate < 1 .or. rate > (huge(rate)-1)/20) return
        if (.not. ieee_is_finite(beta) .or. beta < 0.0_dp .or. beta > 100.0_dp) return
        half = 10*rate
        if (present(window_parameter) .or. present(window_parameter2) .or. &
                present(window_coefficients)) then
            ! FIR prototype with a caller-parameterized symmetric window.
            ! Preserve the previous (unparameterized) Kaiser implementation.
            win = spectral_named_window(window, 2*half+1, st, &
                parameter=window_parameter, parameter2=window_parameter2, &
                coefficients=window_coefficients, fftbins=.false.)
            if (st /= spectral_status_success) return
            deallocate(taps)
            allocate(taps(2*half+1))
            do i = 1, size(taps)
                u = real(i-1-half,dp)/real(rate,dp)
                if (i == half+1) then
                    taps(i) = win(i)/real(rate,dp)
                else
                    taps(i) = win(i)*sin(pi*u)/(pi*real(i-1-half,dp))
                end if
            end do
            norm = sum(taps)
            if (abs(norm) <= tiny(norm) .or. .not. ieee_is_finite(norm)) then
                deallocate(taps)
                allocate(taps(0))
                return
            end if
            taps = taps/norm
        else if (trim(window) /= 'kaiser') then
            taps = firwin(2*half+1,[1.0_dp/real(rate,dp)],window=window,status=st)
            if (st /= design_status_success) then
                if (allocated(taps)) deallocate(taps)
                allocate(taps(0))
                return
            end if
        else
            deallocate(taps)
            allocate(taps(2*half+1))
            norm = i0_series(beta)
            do i=1,size(taps)
                r = real(i-1-half,dp)/real(half,dp)
                w = i0_series(beta*sqrt(max(0.0_dp,1.0_dp-r*r)))/norm
                u = real(i-1-half,dp)/real(rate,dp)
                if (i==half+1) then
                    taps(i)=w/real(rate,dp)
                else
                    taps(i)=w*sin(pi*u)/(pi*real(i-1-half,dp))
                end if
            end do
            taps = taps / sum(taps)
        end if
        status = resample_status_success
    end function default_taps

    pure logical function complex_less(a, b) result(less)
        complex(dp), intent(in) :: a !! first value
        complex(dp), intent(in) :: b !! second value
        less = real(a,dp) < real(b,dp)
        if (real(a,dp) == real(b,dp)) less = aimag(a) < aimag(b)
    end function complex_less

    recursive subroutine sort_complex(v, left, right)
        complex(dp), intent(inout) :: v(:) !! complex values, sorted lexicographically
        integer, intent(in) :: left !! first index
        integer, intent(in) :: right !! last index
        integer :: i, j
        complex(dp) :: pivot, temp
        if (left >= right) return
        i = left
        j = right
        pivot = v(left + (right-left)/2)
        do
            do while (complex_less(v(i),pivot))
                i = i+1
            end do
            do while (complex_less(pivot,v(j)))
                j = j-1
            end do
            if (i > j) exit
            temp=v(i)
            v(i)=v(j)
            v(j)=temp
            i=i+1
            j=j-1
            if (i > j) exit
        end do
        if (left < j) call sort_complex(v,left,j)
        if (i < right) call sort_complex(v,i,right)
    end subroutine sort_complex

    function padding_background(x, mode) result(background)
        complex(dp), intent(in) :: x(:) !! signal samples, finite and nonempty
        character(*), intent(in) :: mode !! mean, median, minimum, or maximum
        complex(dp) :: background
        complex(dp), allocatable :: work(:)
        integer :: i, n
        n = size(x)
        select case (mode)
        case ('mean')
            background = sum(x)/real(n,dp)
        case ('median')
            work = x
            call sort_complex(work,1,n)
            if (mod(n,2)==1) then
                background = work((n+1)/2)
            else
                background = 0.5_dp*work(n/2) + 0.5_dp*work(n/2+1)
            end if
        case ('minimum')
            background = x(1)
            do i = 2, n
                if (complex_less(x(i),background)) background = x(i)
            end do
        case ('maximum')
            background = x(1)
            do i = 2, n
                if (complex_less(background,x(i))) background = x(i)
            end do
        case default
            background = cmplx(0.0_dp,0.0_dp,dp)
        end select
    end function padding_background

    function resample_poly_complex(x, up, down, taps, window, beta, padtype, cval, status, ctaps, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:) !! input samples
        integer, intent(in) :: up !! integer upsampling ratio numerator
        integer, intent(in) :: down !! integer downsampling ratio denominator
        real(dp), intent(in), optional :: taps(:) !! optional user-designed real FIR kernel
        complex(dp), intent(in), optional :: ctaps(:) !! optional user-designed complex FIR kernel
        character(*), intent(in), optional :: window !! 'kaiser' (default) or window name
        real(dp), intent(in), optional :: beta !! Kaiser parameter, default 5
        character(*), intent(in), optional :: padtype !! constant (default), mean, edge, line, etc.
        complex(dp), intent(in), optional :: cval !! constant boundary value
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:), z(:), signal(:)
        real(dp), allocatable :: real_taps(:)
        complex(dp), allocatable :: h(:), padded(:)
        complex(dp) :: outer, background
        character(:), allocatable :: mode, wname
        real(dp) :: kb
        integer :: g, nu, nd, n, nout, half, nprepad, npostpad, first, st

        allocate(y(0))
        if (present(status)) status = resample_status_invalid_input
        if (up<1 .or. down<1 .or. size(x)<1) return
        ! Do not silently ignore new named-window options for explicit taps.
        if (present(taps) .or. present(ctaps)) then
            if (present(window_parameter) .or. present(window_parameter2) .or. &
                    present(window_coefficients)) return
        end if
        if (present(beta)) then
            ! beta is the legacy Kaiser shape, not an additional named-window
            ! parameter. Reject combinations rather than ignoring one value.
            if (present(window_parameter) .or. present(window_parameter2) .or. &
                    present(window_coefficients)) return
        end if
        if (.not. all(ieee_is_finite(real(x,dp))) .or. .not. all(ieee_is_finite(aimag(x)))) return
        g = gcd_positive(up,down)
        nu=up/g
        nd=down/g
        if (nu==1 .and. nd==1) then
            y = x
            if (present(status)) status = resample_status_success
            return
        end if
        if (size(x) > (huge(n)-nd)/nu) return
        n=size(x)
        nout=(n*nu+nd-1)/nd
        mode='constant'
        if (present(padtype)) mode=trim(padtype)
        if (.not. mode_valid(mode) .and. mode/='mean' .and. mode/='median' .and. &
            mode/='minimum' .and. mode/='maximum') return
        if (present(cval) .and. mode/='constant') return
        outer=cmplx(0.0_dp,0.0_dp,dp)
        if (present(cval)) outer=cval
        background=cmplx(0.0_dp,0.0_dp,dp)
        if (mode=='mean' .or. mode=='median' .or. mode=='minimum' .or. mode=='maximum') then
            background=padding_background(x,mode)
            mode='constant'
        end if
        wname='kaiser'
        if (present(window)) wname=trim(window)
        kb=5.0_dp
        if (present(beta)) kb=beta
        if (present(taps) .and. present(ctaps)) return
        if (present(ctaps)) then
            if (size(ctaps)<1) return
            if (.not. all(ieee_is_finite(real(ctaps,dp))) .or. &
                .not. all(ieee_is_finite(aimag(ctaps)))) return
            h=ctaps
        else if (present(taps)) then
            if (size(taps)<1) return
            if (.not. all(ieee_is_finite(taps))) return
            h=cmplx(taps,0.0_dp,dp)
        else
            real_taps=default_taps(max(nu,nd),wname,kb,st, &
                window_parameter=window_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients)
            if (st/=resample_status_success) return
            h=cmplx(real_taps,0.0_dp,dp)
        end if
        half=(size(h)-1)/2
        ! Match SciPy's alignment of the first output with the first input.
        nprepad=nd-mod(half,nd)
        first=(half+nprepad)/nd
        npostpad=0
        do while (((n-1)*nu+size(h)+nprepad+npostpad-1)/nd+1 < first+nout)
            npostpad=npostpad+1
        end do
        if (size(h)+nprepad > huge(n)-npostpad) return
        allocate(padded(nprepad+size(h)+npostpad))
        padded=cmplx(0.0_dp,0.0_dp,dp)
        padded(nprepad+1:nprepad+size(h))=real(nu,dp)*h
        signal=x-background
        z=upfirdn(padded,signal,up=nu,down=nd,mode=mode,cval=outer, status=st)
        if (st /= resample_status_success) return
        if (size(z)<first+nout) return
        y=z(first+1:first+nout)+background
        if (present(status)) status=resample_status_success
    end function resample_poly_complex

    function resample_poly_real(x, up, down, taps, window, beta, padtype, cval, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:) !! input signal
        integer, intent(in) :: up !! upsampling factor
        integer, intent(in) :: down !! downsampling factor
        real(dp), intent(in), optional :: taps(:) !! optional FIR taps
        character(*), intent(in), optional :: window !! FIR window name
        real(dp), intent(in), optional :: beta !! Kaiser shape
        character(*), intent(in), optional :: padtype !! boundary extension mode
        real(dp), intent(in), optional :: cval !! constant boundary value
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: z(:)
        complex(dp) :: outer
        outer=cmplx(0.0_dp,0.0_dp,dp)
        if (present(cval)) outer=cmplx(cval,0.0_dp,dp)
        if (present(cval)) then
            z=resample_poly_complex(cmplx(x,0.0_dp,dp),up,down,taps,window,beta, &
                padtype,outer,status, &
                window_parameter=window_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients)
        else
            z=resample_poly_complex(cmplx(x,0.0_dp,dp),up,down,taps,window,beta, &
                padtype,status=status, &
                window_parameter=window_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients)
        end if
        y=real(z,dp)
    end function resample_poly_real

    ! Complex custom-tap overloads. A complex FIR applied to real input
    ! necessarily produces a complex output, as in SciPy resample_poly.
    function resample_poly_complex_taps(x, up, down, taps, window, beta, padtype, cval, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:) !! complex input samples
        integer, intent(in) :: up !! upsampling factor
        integer, intent(in) :: down !! downsampling factor
        complex(dp), intent(in) :: taps(:) !! complex FIR coefficients
        character(*), intent(in), optional :: window !! used only for designed filters
        real(dp), intent(in), optional :: beta !! used only for designed filters
        character(*), intent(in), optional :: padtype !! boundary mode
        complex(dp), intent(in), optional :: cval !! constant boundary value
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:)
        y = resample_poly_complex(x,up,down,window=window,beta=beta, &
            padtype=padtype,cval=cval,status=status,ctaps=taps, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients)
    end function resample_poly_complex_taps

    function resample_poly_real_taps(x, up, down, taps, window, beta, padtype, cval, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:) !! real input samples
        integer, intent(in) :: up !! upsampling factor
        integer, intent(in) :: down !! downsampling factor
        complex(dp), intent(in) :: taps(:) !! complex FIR coefficients
        character(*), intent(in), optional :: window !! used only for designed filters
        real(dp), intent(in), optional :: beta !! used only for designed filters
        character(*), intent(in), optional :: padtype !! boundary mode
        complex(dp), intent(in), optional :: cval !! constant boundary value
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:)
        y = resample_poly_complex(cmplx(x,0.0_dp,dp),up,down,window=window, &
            beta=beta,padtype=padtype,cval=cval,status=status,ctaps=taps, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients)
    end function resample_poly_real_taps

    function decimate_real(x, q, n, ftype, zero_phase, status, taps, sos, b, a) result(y)
        real(dp), intent(in) :: x(:) !! finite input samples
        integer, intent(in) :: q !! positive downsampling factor
        integer, intent(in), optional :: n !! designed-filter order
        character(*), intent(in), optional :: ftype !! 'fir' or 'iir' (default)
        logical, intent(in), optional :: zero_phase !! compensate FIR delay or filter IIR both ways
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: taps(:) !! explicit FIR coefficients (alternative to ftype/n)
        real(dp), intent(in), optional :: sos(:,:) !! explicit IIR sections (nsections,6), a0=1
        real(dp), intent(in), optional :: b(:) !! numerator coefficients for a custom transfer function
        real(dp), intent(in), optional :: a(:) !! denominator coefficients, starting with nonzero a(1)
        real(dp), allocatable :: y(:), h(:), temp(:), sections(:,:)
        type(butter_result) :: filt
        character(:), allocatable :: method
        integer :: order, st, nout
        logical :: centered

        allocate(y(0))
        if (present(status)) status = resample_status_invalid_input
        if (q < 1 .or. size(x) < 1 .or. q > (huge(q)-1)/20) return
        if (.not. all(ieee_is_finite(x))) return
        ! Explicit taps or SOS replace the designed-filter options, not augment them.
        if (present(b) .neqv. present(a)) return
        if (present(taps) .and. present(sos)) return
        if (present(b)) then
            if (present(taps) .or. present(sos)) return
        end if
        if (present(taps) .or. present(sos) .or. present(b)) then
            if (present(n) .or. present(ftype)) return
        end if
        method = 'iir'
        if (present(ftype)) method = trim(ftype)
        if (present(taps)) then
            if (size(taps) < 1) return
            if (.not. all(ieee_is_finite(taps))) return
            method = 'custom_fir'
        end if
        if (present(sos)) then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            if (.not. all(ieee_is_finite(sos))) return
            if (any(sos(:,4) /= 1.0_dp)) return
            method = 'custom_sos'
        end if
        if (present(b)) then
            if (size(b) < 1 .or. size(a) < 1) return
            if (abs(a(1)) < tiny(1.0_dp)) return
            if (.not. all(ieee_is_finite(real(b,dp)))) return
            if (.not. all(ieee_is_finite(real(a,dp)))) return
            method = 'custom_tf'
            if (size(a) == 1) method = 'custom_fir'
        end if
        if (method /= 'fir' .and. method /= 'iir' .and. &
                method /= 'custom_fir' .and. method /= 'custom_sos' .and. &
                method /= 'custom_tf') return
        centered = .true.
        if (present(zero_phase)) centered = zero_phase
        st = resample_status_invalid_input
        if (method == 'fir' .or. method == 'custom_fir') then
            if (present(b)) then
                h = b / a(1)
            else if (present(taps)) then
                h = taps
            else
                order = 20*q
                if (present(n)) order = n
                if (order < 1) return
                h = firwin(order+1,[1.0_dp/real(q,dp)],window='hamming',status=st)
                if (st /= design_status_success) return
            end if
            if (centered) then
                y = resample_poly(x,1,q,taps=h,status=st)
            else
                temp = upfirdn(h,x,down=q,status=st)
                nout = (size(x)-1)/q+1
                if (st == resample_status_success .and. size(temp) >= nout) y = temp(:nout)
            end if
        else if (method == 'custom_tf') then
            if (centered) then
                temp = filtfilt(b,a,x,status=st)
                if (st /= filtfilt_status_success) return
            else
                temp = lfilter(b,a,x,status=st)
                if (st /= lfilter_status_success) return
            end if
            y = temp(1:size(temp):q)
            st = resample_status_success
        else
            if (present(sos)) then
                sections = sos
            else
                order = 8
                if (present(n)) order = n
                if (order < 1 .or. order > 20) return
                filt = cheby1(order,0.05_dp,0.8_dp/real(q,dp))
                if (filt%status /= design_status_success) return
                sections = filt%sos
            end if
            if (centered) then
                temp = sosfiltfilt(sections,x,status=st)
                if (st /= filtfilt_status_success) return
            else
                temp = sosfilt(sections,x,status=st)
                if (st /= sos_status_success) return
            end if
            y = temp(1:size(temp):q)
            st = resample_status_success
        end if
        if (st /= resample_status_success .or. size(y) /= (size(x)-1)/q+1) then
            if (allocated(y)) deallocate(y)
            allocate(y(0))
            return
        end if
        if (present(status)) status = resample_status_success
    end function decimate_real

    function decimate_complex(x, q, n, ftype, zero_phase, status, taps, sos, b, a) result(y)
        complex(dp), intent(in) :: x(:) !! finite complex signal
        integer, intent(in) :: q !! positive downsampling factor
        integer, intent(in), optional :: n !! designed-filter order
        character(*), intent(in), optional :: ftype !! 'fir' or 'iir' (default)
        logical, intent(in), optional :: zero_phase !! compensate group delay or filter both ways
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: taps(:) !! explicit real FIR coefficients
        real(dp), intent(in), optional :: sos(:,:) !! explicit real IIR sections (nsections,6)
        complex(dp), intent(in), optional :: b(:) !! numerator coefficients for a custom transfer function
        complex(dp), intent(in), optional :: a(:) !! denominator coefficients, starting with nonzero a(1)
        complex(dp), allocatable :: y(:), temp(:)
        complex(dp), allocatable :: h(:)
        real(dp), allocatable :: sections(:,:)
        type(butter_result) :: filt
        character(:), allocatable :: method
        integer :: order, st, nout
        logical :: centered

        allocate(y(0))
        if (present(status)) status = resample_status_invalid_input
        if (q < 1 .or. size(x) < 1 .or. q > (huge(q)-1)/20) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        if (present(b) .neqv. present(a)) return
        if (present(taps) .and. present(sos)) return
        if (present(b)) then
            if (present(taps) .or. present(sos)) return
        end if
        if (present(taps) .or. present(sos) .or. present(b)) then
            if (present(n) .or. present(ftype)) return
        end if
        method = 'iir'
        if (present(ftype)) method = trim(ftype)
        if (present(taps)) then
            if (size(taps) < 1) return
            if (.not. all(ieee_is_finite(taps))) return
            method = 'custom_fir'
        end if
        if (present(sos)) then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            if (.not. all(ieee_is_finite(sos))) return
            if (any(sos(:,4) /= 1.0_dp)) return
            method = 'custom_sos'
        end if
        if (present(b)) then
            if (size(b) < 1 .or. size(a) < 1) return
            if (abs(a(1)) < tiny(1.0_dp)) return
            if (.not. all(ieee_is_finite(real(b,dp)))) return
            if (.not. all(ieee_is_finite(real(a,dp)))) return
            if (.not. all(ieee_is_finite(aimag(b)))) return
            if (.not. all(ieee_is_finite(aimag(a)))) return
            method = 'custom_tf'
            if (size(a) == 1) method = 'custom_fir'
        end if
        if (method /= 'fir' .and. method /= 'iir' .and. &
                method /= 'custom_fir' .and. method /= 'custom_sos' .and. &
                method /= 'custom_tf') return
        centered = .true.
        if (present(zero_phase)) centered = zero_phase
        st = resample_status_invalid_input
        if (method == 'fir' .or. method == 'custom_fir') then
            if (present(b)) then
                h = b / a(1)
            else if (present(taps)) then
                h = taps
            else
                order = 20*q
                if (present(n)) order = n
                if (order < 1) return
                h = firwin(order+1,[1.0_dp/real(q,dp)],window='hamming',status=st)
                if (st /= design_status_success) return
            end if
            if (centered) then
                y = resample_poly(x,1,q,taps=h,status=st)
            else
                temp = upfirdn(h,x,down=q,status=st)
                nout = (size(x)-1)/q+1
                if (st == resample_status_success .and. size(temp) >= nout) y = temp(:nout)
            end if
        else if (method == 'custom_tf') then
            if (centered) then
                temp = filtfilt(b,a,x,status=st)
                if (st /= filtfilt_status_success) return
            else
                temp = lfilter(b,a,x,status=st)
                if (st /= lfilter_status_success) return
            end if
            y = temp(1:size(temp):q)
            st = resample_status_success
        else
            if (present(sos)) then
                sections = sos
            else
                order = 8
                if (present(n)) order = n
                if (order < 1 .or. order > 20) return
                filt = cheby1(order,0.05_dp,0.8_dp/real(q,dp))
                if (filt%status /= design_status_success) return
                sections = filt%sos
            end if
            if (centered) then
                temp = sosfiltfilt(cmplx(sections,0.0_dp,dp),x,status=st)
                if (st /= filtfilt_status_success) return
            else
                temp = sosfilt(cmplx(sections,0.0_dp,dp),x,status=st)
                if (st /= sos_status_success) return
            end if
            y = temp(1:size(temp):q)
            st = resample_status_success
        end if
        if (st /= resample_status_success .or. size(y) /= (size(x)-1)/q+1) then
            if (allocated(y)) deallocate(y)
            allocate(y(0))
            return
        end if
        if (present(status)) status = resample_status_success
    end function decimate_complex


    ! A complex SOS has a different coefficient type than the legacy real-SOS
    ! keyword, so dedicated generic specifics preserve all earlier call forms.
    function decimate_complex_sos(x, q, sos, zero_phase, status) result(y)
        complex(dp), intent(in) :: x(:) !! finite signal to downsample
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: sos(:,:) !! complex SOS, (nsections,6), a0=1
        logical, intent(in), optional :: zero_phase !! forward-backward IIR filtering (default true)
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: filtered(:)
        integer :: st
        logical :: centered

        allocate(y(0))
        if (present(status)) status = resample_status_invalid_input
        if (q < 1 .or. size(x) < 1) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        centered = .true.
        if (present(zero_phase)) centered = zero_phase
        if (centered) then
            filtered = sosfiltfilt(sos, x, status=st)
            if (st /= filtfilt_status_success) return
        else
            filtered = sosfilt(sos, x, status=st)
            if (st /= sos_status_success) return
        end if
        y = filtered(1:size(filtered):q)
        if (present(status)) status = resample_status_success
    end function decimate_complex_sos

    function decimate_real_complex_sos(x, q, sos, zero_phase, status) result(y)
        real(dp), intent(in) :: x(:) !! finite real signal
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: sos(:,:) !! complex SOS, (nsections,6), a0=1
        logical, intent(in), optional :: zero_phase !! forward-backward filtering (default true)
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)

        y = decimate_complex_sos(cmplx(x,0.0_dp,dp), q, sos, zero_phase, status)
    end function decimate_real_complex_sos

    function decimate_nd_complex_sos2(x, q, sos, zero_phase, axis, status) result(y)
        complex(dp), intent(in) :: x(:,:) !! rank-two signal
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: sos(:,:) !! complex SOS, (nsections,6), a0=1
        logical, intent(in), optional :: zero_phase !! forward-backward filtering (default true)
        integer, intent(in), optional :: axis !! one-based or negative axis (default last)
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: filtered(:,:)
        integer :: along, st
        logical :: centered

        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along < 0) along = along + 3
        if (along < 1 .or. along > 2 .or. q < 1) return
        if (any(shape(x) < 1)) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        centered = .true.
        if (present(zero_phase)) centered = zero_phase
        if (centered) then
            filtered = sosfiltfilt(sos, x, axis=along, status=st)
            if (st /= filtfilt_status_success) return
        else
            filtered = sosfilt(sos, x, axis=along, status=st)
            if (st /= sos_status_success) return
        end if
        if (along == 1) then
            y = filtered(1:size(filtered,1):q,:)
        else
            y = filtered(:,1:size(filtered,2):q)
        end if
        if (present(status)) status = resample_status_success
    end function decimate_nd_complex_sos2

    function decimate_nd_complex_sos3(x, q, sos, zero_phase, axis, status) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! rank-three signal
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: sos(:,:) !! complex SOS, (nsections,6), a0=1
        logical, intent(in), optional :: zero_phase !! forward-backward filtering (default true)
        integer, intent(in), optional :: axis !! one-based or negative axis (default last)
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: filtered(:,:,:)
        integer :: along, st
        logical :: centered

        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 0) along = along + 4
        if (along < 1 .or. along > 3 .or. q < 1) return
        if (any(shape(x) < 1)) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        centered = .true.
        if (present(zero_phase)) centered = zero_phase
        if (centered) then
            filtered = sosfiltfilt(sos, x, axis=along, status=st)
            if (st /= filtfilt_status_success) return
        else
            filtered = sosfilt(sos, x, axis=along, status=st)
            if (st /= sos_status_success) return
        end if
        select case(along)
        case (1)
            y = filtered(1:size(filtered,1):q,:,:)
        case (2)
            y = filtered(:,1:size(filtered,2):q,:)
        case (3)
            y = filtered(:,:,1:size(filtered,3):q)
        end select
        if (present(status)) status = resample_status_success
    end function decimate_nd_complex_sos3

    function decimate_nd_real_complex_sos2(x, q, sos, zero_phase, axis, status) result(y)
        real(dp), intent(in) :: x(:,:) !! rank-two real signal
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: sos(:,:) !! complex SOS, (nsections,6), a0=1
        logical, intent(in), optional :: zero_phase !! forward-backward filtering (default true)
        integer, intent(in), optional :: axis !! one-based or negative axis (default last)
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)

        y = decimate_nd_complex_sos2(cmplx(x,0.0_dp,dp),q,sos,zero_phase,axis,status)
    end function decimate_nd_real_complex_sos2

    function decimate_nd_real_complex_sos3(x, q, sos, zero_phase, axis, status) result(y)
        real(dp), intent(in) :: x(:,:,:) !! rank-three real signal
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: sos(:,:) !! complex SOS, (nsections,6), a0=1
        logical, intent(in), optional :: zero_phase !! forward-backward filtering (default true)
        integer, intent(in), optional :: axis !! one-based or negative axis (default last)
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:)

        y = decimate_nd_complex_sos3(cmplx(x,0.0_dp,dp),q,sos,zero_phase,axis,status)
    end function decimate_nd_real_complex_sos3

    function decimate_real_complex_tf(x, q, b, a, zero_phase, status) result(y)
        real(dp), intent(in) :: x(:) !! real input signal
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: b(:) !! complex transfer-function numerator, ascending lag
        complex(dp), intent(in) :: a(:) !! complex denominator, a(1) nonzero
        logical, intent(in), optional :: zero_phase !! zero-phase filtering if true
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        y = decimate_complex(cmplx(x,0.0_dp,dp),q,zero_phase=zero_phase, &
            status=status,b=b,a=a)
    end function decimate_real_complex_tf

    function decimate_nd_rc2(x, q, b, a, zero_phase, axis, status) result(y)
        real(dp), intent(in) :: x(:,:) !! real rank-2 input
        integer, intent(in) :: q !! downsampling factor
        complex(dp), intent(in) :: b(:) !! complex numerator coefficients
        complex(dp), intent(in) :: a(:) !! complex denominator coefficients
        logical, intent(in), optional :: zero_phase !! zero-phase filtering if true
        integer, intent(in), optional :: axis !! selected axis (one-based or negative)
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        y = decimate_nd_c2(cmplx(x,0.0_dp,dp),q,zero_phase=zero_phase, &
            axis=axis,status=status,b=b,a=a)
    end function decimate_nd_rc2

    function decimate_nd_rc3(x, q, b, a, zero_phase, axis, status) result(y)
        real(dp), intent(in) :: x(:,:,:) !! real rank-3 input
        integer, intent(in) :: q !! downsampling factor
        complex(dp), intent(in) :: b(:) !! complex numerator coefficients
        complex(dp), intent(in) :: a(:) !! complex denominator coefficients
        logical, intent(in), optional :: zero_phase !! zero-phase filtering if true
        integer, intent(in), optional :: axis !! selected axis (one-based or negative)
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:)
        y = decimate_nd_c3(cmplx(x,0.0_dp,dp),q,zero_phase=zero_phase, &
            axis=axis,status=status,b=b,a=a)
    end function decimate_nd_rc3

    ! BEGIN GENERATED MULTIRATE ND WRAPPERS
    function upfirdn_nd_rr2(h, x, up, down, mode, cval, axis, status) result(y)
        real(dp), intent(in) :: h(:) !! FIR coefficients
        real(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in), optional :: up !! integer upsampling factor
        integer, intent(in), optional :: down !! integer downsampling factor
        character(*), intent(in), optional :: mode !! extension mode
        real(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 2
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        select case (a)
        case (1)
            line = upfirdn_rr(h, x(:,1), up, down, mode, cval, status=code)
        case (2)
            line = upfirdn_rr(h, x(1,:), up, down, mode, cval, status=code)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = upfirdn_rr(h, x(:,i), up, down, mode, cval, status=code)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = upfirdn_rr(h, x(i,:), up, down, mode, cval, status=code)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function upfirdn_nd_rr2

    function upfirdn_nd_rc2(h, x, up, down, mode, cval, axis, status) result(y)
        real(dp), intent(in) :: h(:) !! FIR coefficients
        complex(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in), optional :: up !! integer upsampling factor
        integer, intent(in), optional :: down !! integer downsampling factor
        character(*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 2
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        select case (a)
        case (1)
            line = upfirdn_rc(h, x(:,1), up, down, mode, cval, status=code)
        case (2)
            line = upfirdn_rc(h, x(1,:), up, down, mode, cval, status=code)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = upfirdn_rc(h, x(:,i), up, down, mode, cval, status=code)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = upfirdn_rc(h, x(i,:), up, down, mode, cval, status=code)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function upfirdn_nd_rc2

    function upfirdn_nd_cr2(h, x, up, down, mode, cval, axis, status) result(y)
        complex(dp), intent(in) :: h(:) !! FIR coefficients
        real(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in), optional :: up !! integer upsampling factor
        integer, intent(in), optional :: down !! integer downsampling factor
        character(*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 2
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        select case (a)
        case (1)
            line = upfirdn_cr(h, x(:,1), up, down, mode, cval, status=code)
        case (2)
            line = upfirdn_cr(h, x(1,:), up, down, mode, cval, status=code)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = upfirdn_cr(h, x(:,i), up, down, mode, cval, status=code)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = upfirdn_cr(h, x(i,:), up, down, mode, cval, status=code)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function upfirdn_nd_cr2

    function upfirdn_nd_cc2(h, x, up, down, mode, cval, axis, status) result(y)
        complex(dp), intent(in) :: h(:) !! FIR coefficients
        complex(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in), optional :: up !! integer upsampling factor
        integer, intent(in), optional :: down !! integer downsampling factor
        character(*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 2
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        select case (a)
        case (1)
            line = upfirdn_cc(h, x(:,1), up, down, mode, cval, status=code)
        case (2)
            line = upfirdn_cc(h, x(1,:), up, down, mode, cval, status=code)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = upfirdn_cc(h, x(:,i), up, down, mode, cval, status=code)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = upfirdn_cc(h, x(i,:), up, down, mode, cval, status=code)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function upfirdn_nd_cc2

    function upfirdn_nd_rr3(h, x, up, down, mode, cval, axis, status) result(y)
        real(dp), intent(in) :: h(:) !! FIR coefficients
        real(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in), optional :: up !! integer upsampling factor
        integer, intent(in), optional :: down !! integer downsampling factor
        character(*), intent(in), optional :: mode !! extension mode
        real(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 3
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        select case (a)
        case (1)
            line = upfirdn_rr(h, x(:,1,1), up, down, mode, cval, status=code)
        case (2)
            line = upfirdn_rr(h, x(1,:,1), up, down, mode, cval, status=code)
        case (3)
            line = upfirdn_rr(h, x(1,1,:), up, down, mode, cval, status=code)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = upfirdn_rr(h, x(:,i,j), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = upfirdn_rr(h, x(i,:,j), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = upfirdn_rr(h, x(i,j,:), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function upfirdn_nd_rr3

    function upfirdn_nd_rc3(h, x, up, down, mode, cval, axis, status) result(y)
        real(dp), intent(in) :: h(:) !! FIR coefficients
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in), optional :: up !! integer upsampling factor
        integer, intent(in), optional :: down !! integer downsampling factor
        character(*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 3
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        select case (a)
        case (1)
            line = upfirdn_rc(h, x(:,1,1), up, down, mode, cval, status=code)
        case (2)
            line = upfirdn_rc(h, x(1,:,1), up, down, mode, cval, status=code)
        case (3)
            line = upfirdn_rc(h, x(1,1,:), up, down, mode, cval, status=code)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = upfirdn_rc(h, x(:,i,j), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = upfirdn_rc(h, x(i,:,j), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = upfirdn_rc(h, x(i,j,:), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function upfirdn_nd_rc3

    function upfirdn_nd_cr3(h, x, up, down, mode, cval, axis, status) result(y)
        complex(dp), intent(in) :: h(:) !! FIR coefficients
        real(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in), optional :: up !! integer upsampling factor
        integer, intent(in), optional :: down !! integer downsampling factor
        character(*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 3
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        select case (a)
        case (1)
            line = upfirdn_cr(h, x(:,1,1), up, down, mode, cval, status=code)
        case (2)
            line = upfirdn_cr(h, x(1,:,1), up, down, mode, cval, status=code)
        case (3)
            line = upfirdn_cr(h, x(1,1,:), up, down, mode, cval, status=code)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = upfirdn_cr(h, x(:,i,j), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = upfirdn_cr(h, x(i,:,j), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = upfirdn_cr(h, x(i,j,:), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function upfirdn_nd_cr3

    function upfirdn_nd_cc3(h, x, up, down, mode, cval, axis, status) result(y)
        complex(dp), intent(in) :: h(:) !! FIR coefficients
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in), optional :: up !! integer upsampling factor
        integer, intent(in), optional :: down !! integer downsampling factor
        character(*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 3
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        select case (a)
        case (1)
            line = upfirdn_cc(h, x(:,1,1), up, down, mode, cval, status=code)
        case (2)
            line = upfirdn_cc(h, x(1,:,1), up, down, mode, cval, status=code)
        case (3)
            line = upfirdn_cc(h, x(1,1,:), up, down, mode, cval, status=code)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = upfirdn_cc(h, x(:,i,j), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = upfirdn_cc(h, x(i,:,j), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = upfirdn_cc(h, x(i,j,:), up, down, mode, cval, status=code)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function upfirdn_nd_cc3

    function resample_poly_nd_r2(x, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in) :: up !! integer upsampling factor
        integer, intent(in) :: down !! integer downsampling factor
        real(dp), intent(in), optional :: taps(:) !! custom real FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser parameter
        character(*), intent(in), optional :: padtype !! boundary extension
        real(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        real(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        select case (a)
        case (1)
            line = resample_poly_real(x(:,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (2)
            line = resample_poly_real(x(1,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = resample_poly_real(x(:,i), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = resample_poly_real(x(i,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function resample_poly_nd_r2

    function resample_poly_nd_c2(x, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in) :: up !! integer upsampling factor
        integer, intent(in) :: down !! integer downsampling factor
        real(dp), intent(in), optional :: taps(:) !! custom real FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser parameter
        character(*), intent(in), optional :: padtype !! boundary extension
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        select case (a)
        case (1)
            line = resample_poly_complex(x(:,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (2)
            line = resample_poly_complex(x(1,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = resample_poly_complex(x(:,i), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = resample_poly_complex(x(i,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function resample_poly_nd_c2

    function resample_poly_nd_r3(x, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in) :: up !! integer upsampling factor
        integer, intent(in) :: down !! integer downsampling factor
        real(dp), intent(in), optional :: taps(:) !! custom real FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser parameter
        character(*), intent(in), optional :: padtype !! boundary extension
        real(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        real(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        select case (a)
        case (1)
            line = resample_poly_real(x(:,1,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (2)
            line = resample_poly_real(x(1,:,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (3)
            line = resample_poly_real(x(1,1,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = resample_poly_real(x(:,i,j), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = resample_poly_real(x(i,:,j), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = resample_poly_real(x(i,j,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function resample_poly_nd_r3

    function resample_poly_nd_c3(x, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in) :: up !! integer upsampling factor
        integer, intent(in) :: down !! integer downsampling factor
        real(dp), intent(in), optional :: taps(:) !! custom real FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser parameter
        character(*), intent(in), optional :: padtype !! boundary extension
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        select case (a)
        case (1)
            line = resample_poly_complex(x(:,1,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (2)
            line = resample_poly_complex(x(1,:,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (3)
            line = resample_poly_complex(x(1,1,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = resample_poly_complex(x(:,i,j), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = resample_poly_complex(x(i,:,j), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = resample_poly_complex(x(i,j,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function resample_poly_nd_c3

    function resample_poly_ctaps_nd_r2(x, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in) :: up !! integer upsampling factor
        integer, intent(in) :: down !! integer downsampling factor
        complex(dp), intent(in) :: taps(:) !! required complex FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser parameter
        character(*), intent(in), optional :: padtype !! boundary extension
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        select case (a)
        case (1)
            line = resample_poly_real_taps(x(:,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (2)
            line = resample_poly_real_taps(x(1,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = resample_poly_real_taps(x(:,i), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = resample_poly_real_taps(x(i,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function resample_poly_ctaps_nd_r2

    function resample_poly_ctaps_nd_c2(x, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in) :: up !! integer upsampling factor
        integer, intent(in) :: down !! integer downsampling factor
        complex(dp), intent(in) :: taps(:) !! required complex FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser parameter
        character(*), intent(in), optional :: padtype !! boundary extension
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        select case (a)
        case (1)
            line = resample_poly_complex_taps(x(:,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (2)
            line = resample_poly_complex_taps(x(1,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = resample_poly_complex_taps(x(:,i), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = resample_poly_complex_taps(x(i,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function resample_poly_ctaps_nd_c2

    function resample_poly_ctaps_nd_r3(x, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in) :: up !! integer upsampling factor
        integer, intent(in) :: down !! integer downsampling factor
        complex(dp), intent(in) :: taps(:) !! required complex FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser parameter
        character(*), intent(in), optional :: padtype !! boundary extension
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        select case (a)
        case (1)
            line = resample_poly_real_taps(x(:,1,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (2)
            line = resample_poly_real_taps(x(1,:,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (3)
            line = resample_poly_real_taps(x(1,1,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = resample_poly_real_taps(x(:,i,j), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = resample_poly_real_taps(x(i,:,j), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = resample_poly_real_taps(x(i,j,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function resample_poly_ctaps_nd_r3

    function resample_poly_ctaps_nd_c3(x, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in) :: up !! integer upsampling factor
        integer, intent(in) :: down !! integer downsampling factor
        complex(dp), intent(in) :: taps(:) !! required complex FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser parameter
        character(*), intent(in), optional :: padtype !! boundary extension
        complex(dp), intent(in), optional :: cval !! constant padding value
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        complex(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        select case (a)
        case (1)
            line = resample_poly_complex_taps(x(:,1,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (2)
            line = resample_poly_complex_taps(x(1,:,1), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        case (3)
            line = resample_poly_complex_taps(x(1,1,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = resample_poly_complex_taps(x(:,i,j), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = resample_poly_complex_taps(x(i,:,j), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = resample_poly_complex_taps(x(i,j,:), up, down, taps, window, beta, padtype, cval, status=code, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                        if (code /= resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function resample_poly_ctaps_nd_c3

    function decimate_nd_r2(x, q, n, ftype, zero_phase, axis, status, taps, sos, b, a) result(y)
        real(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in) :: q !! decimation factor
        integer, intent(in), optional :: n !! filter order
        character(*), intent(in), optional :: ftype !! fir or iir
        logical, intent(in), optional :: zero_phase !! noncausal if true
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: taps(:) !! supplied real FIR coefficients
        real(dp), intent(in), optional :: sos(:,:) !! supplied real IIR SOS coefficients
        real(dp), intent(in), optional :: b(:) !! transfer-function numerator
        real(dp), intent(in), optional :: a(:) !! transfer-function denominator
        real(dp), allocatable :: y(:,:), line(:)
        integer :: ax, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        ax = 2
        if (present(axis)) ax = axis
        if (ax < 0) ax = ax + 2 + 1
        if (ax < 1 .or. ax > 2) return
        select case (ax)
        case (1)
            line = decimate_real(x(:,1), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        case (2)
            line = decimate_real(x(1,:), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(ax) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (ax)
        case (1)
            do i = 1, size(x,2)
                    line = decimate_real(x(:,i), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                    if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = decimate_real(x(i,:), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                    if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function decimate_nd_r2

    function decimate_nd_c2(x, q, n, ftype, zero_phase, axis, status, taps, sos, b, a) result(y)
        complex(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in) :: q !! decimation factor
        integer, intent(in), optional :: n !! filter order
        character(*), intent(in), optional :: ftype !! fir or iir
        logical, intent(in), optional :: zero_phase !! noncausal if true
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: taps(:) !! supplied real FIR coefficients
        real(dp), intent(in), optional :: sos(:,:) !! supplied real IIR SOS coefficients
        complex(dp), intent(in), optional :: b(:) !! transfer-function numerator
        complex(dp), intent(in), optional :: a(:) !! transfer-function denominator
        complex(dp), allocatable :: y(:,:), line(:)
        integer :: ax, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        ax = 2
        if (present(axis)) ax = axis
        if (ax < 0) ax = ax + 2 + 1
        if (ax < 1 .or. ax > 2) return
        select case (ax)
        case (1)
            line = decimate_complex(x(:,1), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        case (2)
            line = decimate_complex(x(1,:), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(ax) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (ax)
        case (1)
            do i = 1, size(x,2)
                    line = decimate_complex(x(:,i), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                    if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = decimate_complex(x(i,:), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                    if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = resample_status_success
    end function decimate_nd_c2

    function decimate_nd_r3(x, q, n, ftype, zero_phase, axis, status, taps, sos, b, a) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in) :: q !! decimation factor
        integer, intent(in), optional :: n !! filter order
        character(*), intent(in), optional :: ftype !! fir or iir
        logical, intent(in), optional :: zero_phase !! noncausal if true
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: taps(:) !! supplied real FIR coefficients
        real(dp), intent(in), optional :: sos(:,:) !! supplied real IIR SOS coefficients
        real(dp), intent(in), optional :: b(:) !! transfer-function numerator
        real(dp), intent(in), optional :: a(:) !! transfer-function denominator
        real(dp), allocatable :: y(:,:,:), line(:)
        integer :: ax, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        ax = 3
        if (present(axis)) ax = axis
        if (ax < 0) ax = ax + 3 + 1
        if (ax < 1 .or. ax > 3) return
        select case (ax)
        case (1)
            line = decimate_real(x(:,1,1), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        case (2)
            line = decimate_real(x(1,:,1), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        case (3)
            line = decimate_real(x(1,1,:), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(ax) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (ax)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = decimate_real(x(:,i,j), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                        if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = decimate_real(x(i,:,j), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                        if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = decimate_real(x(i,j,:), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                        if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function decimate_nd_r3

    function decimate_nd_c3(x, q, n, ftype, zero_phase, axis, status, taps, sos, b, a) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in) :: q !! decimation factor
        integer, intent(in), optional :: n !! filter order
        character(*), intent(in), optional :: ftype !! fir or iir
        logical, intent(in), optional :: zero_phase !! noncausal if true
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: taps(:) !! supplied real FIR coefficients
        real(dp), intent(in), optional :: sos(:,:) !! supplied real IIR SOS coefficients
        complex(dp), intent(in), optional :: b(:) !! transfer-function numerator
        complex(dp), intent(in), optional :: a(:) !! transfer-function denominator
        complex(dp), allocatable :: y(:,:,:), line(:)
        integer :: ax, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = resample_status_invalid_input
        if (any(shape(x) == 0)) return
        ax = 3
        if (present(axis)) ax = axis
        if (ax < 0) ax = ax + 3 + 1
        if (ax < 1 .or. ax > 3) return
        select case (ax)
        case (1)
            line = decimate_complex(x(:,1,1), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        case (2)
            line = decimate_complex(x(1,:,1), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        case (3)
            line = decimate_complex(x(1,1,:), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
        end select
        if (code /= resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(ax) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (ax)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = decimate_complex(x(:,i,j), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                        if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = decimate_complex(x(i,:,j), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                        if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = decimate_complex(x(i,j,:), q, n, ftype, zero_phase, status=code, taps=taps, sos=sos, b=b, a=a)
                        if (code /= resample_status_success .or. size(line) /= sizes(ax)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = resample_status_success
    end function decimate_nd_c3
    ! END GENERATED MULTIRATE ND WRAPPERS
end module scifort_signal_resample
