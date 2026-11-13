! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
!
! Adapted from SciPy v1.15.3 scipy/fft/_fftlog_backend.py (BSD-3-Clause).
! Full upstream license: LICENSES/scipy-fftlog-BSD-3-Clause.txt.
! Hamilton's FFTLog equations specify the Mellin coefficients and bias factors.
! SciFort adds explicit status/range checks, logarithmic scaling, native FFT and
! log-gamma kernels, and selectable-axis batched routines. See CODE_PROVENANCE.md.
module scifort_fftlog
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi, scifort_log_two
    use scifort_gamma_family, only : loggamma, gammasgn
    use scifort_signal_spectral, only : fft_plan
    implicit none
    private
    public :: fht, ifht, fhtoffset, fht_nd, ifht_nd
    integer, parameter, public :: fft_status_singular = 3
    integer, parameter, public :: fft_status_range_error = 4
    real(dp), parameter :: log_limit = log(huge(1.0_dp))
    real(dp), parameter :: parameter_limit = sqrt(huge(1.0_dp))/64.0_dp
contains

    subroutine fht_nd(a,b,dln,mu,offset,bias,axis,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: a(..) !! finite rank 1-15 logarithmic samples
        real(dp), contiguous, target, intent(out) :: b(..) !! same shape as a; must not overlap a
        real(dp), intent(in) :: dln !! finite nonzero natural-log grid spacing
        real(dp), intent(in) :: mu !! finite real Bessel order
        real(dp), optional, intent(in) :: offset !! log center product, default 0
        real(dp), optional, intent(in) :: bias !! power-law bias exponent, default 0
        integer, optional, intent(in) :: axis !! one-based transform axis, default last
        integer, optional, intent(out) :: status !! 0 success, 1 input, 2 setup, 3 singular, 4 range
        call hankel_nd(a,b,dln,mu,offset,bias,axis,.false.,status)
    end subroutine fht_nd

    subroutine ifht_nd(a,b,dln,mu,offset,bias,axis,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: a(..) !! finite rank 1-15 logarithmic samples
        real(dp), contiguous, target, intent(out) :: b(..) !! same shape as a; must not overlap a
        real(dp), intent(in) :: dln !! finite nonzero natural-log grid spacing
        real(dp), intent(in) :: mu !! finite real Bessel order
        real(dp), optional, intent(in) :: offset !! log center product, default 0
        real(dp), optional, intent(in) :: bias !! power-law bias exponent, default 0
        integer, optional, intent(in) :: axis !! one-based transform axis, default last
        integer, optional, intent(out) :: status !! 0 success, 1 input, 2 setup, 3 singular, 4 range
        call hankel_nd(a,b,dln,mu,offset,bias,axis,.true.,status)
    end subroutine ifht_nd

    subroutine hankel_nd(a,b,dln,mu,offset,bias,axis,inverse,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: a(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: b(..) !! expected result shape
        real(dp), intent(in) :: dln !! logarithmic spacing
        real(dp), intent(in) :: mu !! Bessel order
        real(dp), optional, intent(in) :: offset !! log center product
        real(dp), optional, intent(in) :: bias !! power-law bias
        integer, optional, intent(in) :: axis !! one-based axis, default last
        logical, intent(in) :: inverse !! inverse transform when true
        integer, optional, intent(out) :: status !! completion status
        real(dp), pointer :: input(:), output(:)
        real(dp), allocatable :: buffer(:), line(:), transformed(:)
        complex(dp), allocatable :: weights(:)
        type(fft_plan) :: plan
        real(dp) :: q, shift, coefficient_scale
        integer :: dims(15), r, dim, n, stride, blocks, block, off, base, j, ierr, code
        code = 1
        if (present(status)) status = code
        r = rank(b)
        if (r < 1 .or. r > 15) return
        if (size(b,kind=int64) > int(huge(0),int64)) return
        if (size(b) > 0) then
            call c_f_pointer(c_loc(b),output,[size(b)])
            output = ieee_value(0.0_dp,ieee_quiet_nan)
        end if
        if (rank(a) /= r) return
        if (any(shape(a) /= shape(b))) return
        dims(:r) = shape(a)
        dim = r
        if (present(axis)) dim = axis
        if (dim < 1 .or. dim > r) return
        n = dims(dim)
        call parameters(n,dln,mu,offset,bias,q,shift,code)
        if (code /= 0) then
            if (present(status)) status = code
            return
        end if
        if (size(b) == 0) then
            if (present(status)) status = 0
            return
        end if
        call c_f_pointer(c_loc(a),input,[size(a)])
        if (.not. all(ieee_is_finite(input))) return
        call coefficients(n,dln,mu,q,shift,inverse,weights,coefficient_scale,code)
        if (code == 0) then
            call plan%initialize(n,ierr)
            if (ierr /= 0) code = 2
        end if
        if (code == 0) then
            allocate(buffer(size(input)),line(n),transformed(n),stat=ierr)
            if (ierr /= 0) code = 2
        end if
        if (code == 0) then
            stride = product(dims(:dim-1))
            blocks = product(dims(dim+1:r))
            outer: do block = 0, blocks-1
                do off = 1, stride
                    base = block*stride*n+off
                    do j = 1, n
                        line(j) = input(base+(j-1)*stride)
                    end do
                    call transform_line(line,transformed,dln,q,shift,inverse,weights, &
                        coefficient_scale,plan,code)
                    if (code /= 0) exit outer
                    do j = 1, n
                        buffer(base+(j-1)*stride) = transformed(j)
                    end do
                end do
            end do outer
            if (code == 0) output = buffer
        end if
        if (present(status)) status = code
    end subroutine hankel_nd

    function fht(a,dln,mu,offset,bias,status) result(y)
        implicit none
        real(dp), intent(in) :: a(:) !! finite samples on a uniform logarithmic grid
        real(dp), intent(in) :: dln !! finite nonzero natural-log grid spacing
        real(dp), intent(in) :: mu !! finite real Bessel order
        real(dp), optional, intent(in) :: offset !! log of grid center product, default 0
        real(dp), optional, intent(in) :: bias !! power-law bias exponent, default 0
        integer, optional, intent(out) :: status !! 0 success, 1 input, 2 setup, 3 singular, 4 range
        real(dp), allocatable :: y(:)
        y = transform(a,dln,mu,offset,bias,.false.,status)
    end function fht

    function ifht(a,dln,mu,offset,bias,status) result(y)
        implicit none
        real(dp), intent(in) :: a(:) !! finite Hankel coefficients on the logarithmic output grid
        real(dp), intent(in) :: dln !! finite nonzero natural-log grid spacing
        real(dp), intent(in) :: mu !! finite real Bessel order of the forward transform
        real(dp), optional, intent(in) :: offset !! forward grid center log product, default 0
        real(dp), optional, intent(in) :: bias !! forward power-law bias exponent, default 0
        integer, optional, intent(out) :: status !! 0 success, 1 input, 2 setup, 3 singular, 4 range
        real(dp), allocatable :: y(:)
        y = transform(a,dln,mu,offset,bias,.true.,status)
    end function ifht

    function transform(a,dln,mu,offset,bias,inverse,status) result(y)
        implicit none
        real(dp), intent(in) :: a(:) !! finite logarithmically sampled data
        real(dp), intent(in) :: dln !! logarithmic spacing
        real(dp), intent(in) :: mu !! Bessel order
        real(dp), optional, intent(in) :: offset !! log grid center product
        real(dp), optional, intent(in) :: bias !! bias exponent
        logical, intent(in) :: inverse !! inverse transform when true
        integer, optional, intent(out) :: status !! completion code
        real(dp), allocatable :: y(:), work(:)
        complex(dp), allocatable :: weights(:)
        type(fft_plan) :: plan
        real(dp) :: q, shift, coefficient_scale
        integer :: code, ierr
        allocate(y(0))
        call parameters(size(a),dln,mu,offset,bias,q,shift,code)
        if (.not. all(ieee_is_finite(a))) code = 1
        if (code == 0) then
            call coefficients(size(a),dln,mu,q,shift,inverse,weights,coefficient_scale,code)
        end if
        if (code == 0) then
            call plan%initialize(size(a),ierr)
            if (ierr /= 0) code = 2
        end if
        if (code == 0) then
            allocate(work(size(a)),stat=ierr)
            if (ierr /= 0) code = 2
        end if
        if (code == 0) then
            call transform_line(a,work,dln,q,shift,inverse,weights,coefficient_scale,plan,code)
            if (code == 0) call move_alloc(work,y)
        end if
        if (present(status)) status = code
    end function transform

    subroutine parameters(n,dln,mu,offset,bias,q,shift,code)
        implicit none
        integer, intent(in) :: n !! positive transform length
        real(dp), intent(in) :: dln !! finite nonzero logarithmic spacing
        real(dp), intent(in) :: mu !! finite order
        real(dp), optional, intent(in) :: offset !! log center product
        real(dp), optional, intent(in) :: bias !! bias exponent
        real(dp), intent(out) :: q !! resolved bias
        real(dp), intent(out) :: shift !! resolved offset
        integer, intent(out) :: code !! completion code
        real(dp) :: span
        code = 1
        q = 0.0_dp
        shift = 0.0_dp
        if (present(bias)) q = bias
        if (present(offset)) shift = offset
        if (n < 1) return
        if (.not. all(ieee_is_finite([dln,mu,q,shift]))) return
        if (dln == 0.0_dp) return
        code = fft_status_range_error
        if (max(abs(dln),abs(mu),abs(q),abs(shift)) > parameter_limit) return
        ! Keep recurrence-index conversions inside the existing log-gamma
        ! kernel within default-integer range, including reflection.
        if (max(abs(mu),abs(q)) > 0.25_dp*real(huge(0),dp)) return
        if (abs(dln) < 1.0_dp/parameter_limit) return
        span = 0.5_dp*real(n-1,dp)*abs(dln)+abs(shift)
        if (abs(q) > 1.0_dp) then
            if (span > log_limit/abs(q)) return
        else if (span > log_limit) then
            if (abs(q) > log_limit/span) return
        end if
        code = 0
    end subroutine parameters

    subroutine coefficients(n,dln,mu,q,shift,inverse,weights,scale_log,code)
        implicit none
        integer, intent(in) :: n !! positive FFT length
        real(dp), intent(in) :: dln !! validated grid spacing
        real(dp), intent(in) :: mu !! validated Bessel order
        real(dp), intent(in) :: q !! validated bias exponent
        real(dp), intent(in) :: shift !! validated logarithmic offset
        logical, intent(in) :: inverse !! reciprocal conjugate multiplier when true
        complex(dp), allocatable, intent(out) :: weights(:) !! normalized Hermitian multipliers
        real(dp), intent(out) :: scale_log !! common logarithmic multiplier magnitude
        integer, intent(out) :: code !! completion code
        real(dp), allocatable :: logs(:), phases(:)
        complex(dp) :: plus, minus
        real(dp) :: xp, xm, y, projection
        integer :: m, half, ierr, first
        logical :: pole_plus, pole_minus, zero_dc
        code = 2
        half = n/2
        allocate(weights(n),logs(half+1),phases(half+1),stat=ierr)
        if (ierr /= 0) return
        xp = 0.5_dp*mu+0.5_dp+0.5_dp*q
        xm = 0.5_dp*mu+0.5_dp-0.5_dp*q
        pole_plus = xp <= 0.0_dp .and. xp == anint(xp)
        pole_minus = xm <= 0.0_dp .and. xm == anint(xm)
        code = fft_status_singular
        zero_dc = pole_plus .neqv. pole_minus
        if (zero_dc) then
            if (pole_plus .and. .not. inverse) return
            if (pole_minus .and. inverse) return
            logs(1) = 0.0_dp
            phases(1) = 0.0_dp
        else if (pole_plus) then
            ! Simultaneous gamma poles have a removable ratio: residues give
            ! (-1)**(xp-xm) * factorial(-xm)/factorial(-xp).
            logs(1) = q*scifort_log_two+log_gamma(1.0_dp-xm)-log_gamma(1.0_dp-xp)
            phases(1) = 0.0_dp
            if (modulo(xp-xm,2.0_dp) /= 0.0_dp) phases(1) = scifort_pi
        else
            plus = loggamma(cmplx(xp,0.0_dp,dp))
            minus = loggamma(cmplx(xm,0.0_dp,dp))
            logs(1) = q*scifort_log_two+real(plus,dp)-real(minus,dp)
            phases(1) = 0.0_dp
            if (gammasgn(xp)*gammasgn(xm) < 0.0_dp) phases(1) = scifort_pi
        end if
        do m = 1, half
            y = (scifort_pi*real(m,dp)/real(n,dp))/dln
            plus = loggamma(cmplx(xp,y,dp))
            minus = loggamma(cmplx(xm,y,dp))
            logs(m+1) = q*scifort_log_two+real(plus,dp)-real(minus,dp)
            phases(m+1) = 2.0_dp*y*(scifort_log_two-shift)+aimag(plus)+aimag(minus)
        end do
        code = fft_status_range_error
        if (.not. all(ieee_is_finite(logs)) .or. .not. all(ieee_is_finite(phases))) return
        if (mod(n,2) == 0) then
            projection = cos(phases(half+1))
            code = fft_status_singular
            if (abs(projection) <= 32.0_dp*epsilon(1.0_dp)) return
            logs(half+1) = logs(half+1)+log(abs(projection))
            phases(half+1) = 0.0_dp
            if (projection < 0.0_dp) phases(half+1) = scifort_pi
        end if
        if (inverse) logs = -logs
        first = 1
        weights = (0.0_dp,0.0_dp)
        if (zero_dc) first = 2
        if (first > size(logs)) then
            scale_log = 0.0_dp
            code = 0
            return
        end if
        scale_log = maxval(logs(first:))
        code = fft_status_range_error
        if (scale_log-minval(logs(first:)) > log_limit) return
        do m = first-1, half
            weights(m+1) = exp(logs(m+1)-scale_log)* &
                cmplx(cos(phases(m+1)),sin(phases(m+1)),dp)
        end do
        weights(1) = cmplx(real(weights(1),dp),0.0_dp,dp)
        if (mod(n,2) == 0) weights(half+1) = cmplx(real(weights(half+1),dp),0.0_dp,dp)
        do m = 1, (n-1)/2
            weights(n-m+1) = conjg(weights(m+1))
        end do
        code = 0
    end subroutine coefficients

    subroutine transform_line(a,b,dln,q,shift,inverse,weights,coefficient_scale,plan,code)
        implicit none
        real(dp), intent(in) :: a(:) !! finite input line
        real(dp), intent(out) :: b(:) !! output line, defined on success
        real(dp), intent(in) :: dln !! validated grid spacing
        real(dp), intent(in) :: q !! validated bias
        real(dp), intent(in) :: shift !! validated offset
        logical, intent(in) :: inverse !! inverse transform when true
        complex(dp), intent(in) :: weights(:) !! scaled multiplier spectrum
        real(dp), intent(in) :: coefficient_scale !! common log multiplier scale
        type(fft_plan), intent(inout) :: plan !! local reusable Fourier plan
        integer, intent(out) :: code !! completion code
        complex(dp), allocatable :: values(:)
        real(dp), allocatable :: magnitudes(:)
        real(dp) :: center, coordinate, input_scale, logn, value, exponent_value, direction
        integer :: n, j, ierr
        n = size(a)
        code = 2
        allocate(values(n),magnitudes(n),stat=ierr)
        if (ierr /= 0) return
        if (all(a == 0.0_dp)) then
            b = 0.0_dp
            code = 0
            return
        end if
        center = 0.5_dp*real(n-1,dp)
        logn = log(real(n,dp))
        direction = -1.0_dp
        if (inverse) direction = 1.0_dp
        do j = 1, n
            magnitudes(j) = -huge(1.0_dp)
            if (a(j) == 0.0_dp) cycle
            coordinate = (real(j-1,dp)-center)*dln
            exponent_value = direction*q*coordinate
            if (inverse) exponent_value = exponent_value+q*shift
            magnitudes(j) = log(abs(a(j)))+exponent_value
        end do
        input_scale = maxval(magnitudes)
        do j = 1, n
            values(j) = (0.0_dp,0.0_dp)
            if (a(j) /= 0.0_dp) values(j) = cmplx( &
                sign(exp(magnitudes(j)-input_scale)/real(n,dp),a(j)),0.0_dp,dp)
        end do
        call plan%execute(values,ierr)
        if (ierr /= 0) return
        values = values*weights
        call plan%inverse(values,ierr)
        if (ierr /= 0) return
        code = fft_status_range_error
        do j = 1, n
            value = real(values(n-j+1),dp)
            b(j) = 0.0_dp
            if (value == 0.0_dp) cycle
            coordinate = (real(j-1,dp)-center)*dln
            exponent_value = direction*q*coordinate
            if (.not. inverse) exponent_value = exponent_value-q*shift
            exponent_value = log(abs(value))+input_scale+coefficient_scale+logn+exponent_value
            if (exponent_value > log_limit-4.0_dp*epsilon(1.0_dp)*log_limit) return
            b(j) = sign(exp(exponent_value),value)
        end do
        code = 0
    end subroutine transform_line

    function fhtoffset(dln,mu,initial,bias,status) result(offset)
        implicit none
        real(dp), intent(in) :: dln !! finite nonzero logarithmic spacing
        real(dp), intent(in) :: mu !! finite real Bessel order
        real(dp), optional, intent(in) :: initial !! requested nearby offset, default 0
        real(dp), optional, intent(in) :: bias !! bias exponent, default 0
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 4 numerical range
        real(dp) :: offset, q, shift, xp, xm, y, arg, nearest
        complex(dp) :: plus, minus
        integer :: code
        offset = ieee_value(0.0_dp,ieee_quiet_nan)
        call parameters(1,dln,mu,bias=bias,q=q,shift=shift,code=code)
        if (present(initial)) shift = initial
        if (.not. ieee_is_finite(shift)) code = 1
        if (ieee_is_finite(shift)) then
            if (abs(shift) > parameter_limit) code = fft_status_range_error
        end if
        if (code == 0) then
            xp = 0.5_dp*mu+0.5_dp+0.5_dp*q
            xm = 0.5_dp*mu+0.5_dp-0.5_dp*q
            y = (0.5_dp*scifort_pi)/dln
            plus = loggamma(cmplx(xp,y,dp))
            minus = loggamma(cmplx(xm,y,dp))
            arg = (scifort_log_two-shift)/dln+(aimag(plus)+aimag(minus))/scifort_pi
            code = fft_status_range_error
            if (ieee_is_finite(arg)) then
                ! Match nearest-even ties without converting a potentially
                ! large phase index to a default integer.
                nearest = anint(arg)
                if (abs(arg-nearest) == 0.5_dp .and. modulo(nearest,2.0_dp) /= 0.0_dp) then
                    nearest = nearest-sign(1.0_dp,arg)
                end if
                offset = shift+(arg-nearest)*dln
                code = 0
            end if
        end if
        if (present(status)) status = code
    end function fhtoffset
end module scifort_fftlog
