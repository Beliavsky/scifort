! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Savitzky-Golay polynomial least-squares filtering along array axes.
module scifort_signal_savgol
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: savgol_status_success = 0
    integer, parameter, public :: savgol_status_invalid_input = 1
    integer, parameter, public :: savgol_status_singular = 2
    interface savgol_filter
        module procedure savgol_filter_real_1d
        module procedure savgol_filter_real_2d, savgol_filter_real_3d, savgol_filter_real_4d
        module procedure savgol_filter_complex_1d
        module procedure savgol_filter_complex_2d, savgol_filter_complex_3d, savgol_filter_complex_4d
    end interface savgol_filter

    ! Immutable weights shared by all lines along one filtered axis.  The
    ! center weights and each interpolated-edge fit are computed only once.
    type :: savgol_kernel
        real(dp), allocatable :: center(:), left(:,:), right(:,:)
        character(len=:), allocatable :: boundary
        integer :: n = 0, window = 0, half = 0
    end type savgol_kernel

    interface savgol_filter_nd
        module procedure savgol_filter_nd_real, savgol_filter_nd_complex
    end interface savgol_filter_nd

    interface savgol_irregular
        module procedure savgol_irregular_real, savgol_irregular_complex
    end interface savgol_irregular

    public :: savgol_coeffs, savgol_filter, savgol_filter_nd, savgol_irregular

contains

    function polynomial_weights(window_length, polyorder, deriv, delta, position, ok) result(w)
        integer, intent(in) :: window_length !! number of points in the fit
        integer, intent(in) :: polyorder !! polynomial degree, less than window length
        integer, intent(in) :: deriv !! required derivative order, nonnegative
        real(dp), intent(in) :: delta !! sample spacing, nonzero
        real(dp), intent(in) :: position !! zero-based evaluation position in the window
        logical, intent(out) :: ok !! true when the polynomial fit is numerically nonsingular
        real(dp), allocatable :: w(:)

        real(dp), allocatable :: q(:, :), r(:, :), v(:), rhs(:), z(:)
        real(dp) :: scale, projection, norm_v, nanval
        integer :: i, j, k, pass, p

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(w(window_length))
        w = nanval
        ok = .false.
        p = polyorder + 1
        allocate(q(window_length, p), r(p, p), v(window_length), rhs(p), z(p))
        q = 0.0_dp
        r = 0.0_dp
        scale = max(1.0_dp, abs(position), abs(real(window_length - 1, dp) - position))

        ! Form an orthonormal basis for the scaled Vandermonde design.
        ! Two passes of modified Gram-Schmidt avoid squaring its condition number.
        do j = 1, p
            do i = 1, window_length
                v(i) = ((real(i - 1, dp) - position) / scale)**(j - 1)
            end do
            do pass = 1, 2
                do k = 1, j - 1
                    projection = dot_product(q(:, k), v)
                    r(k, j) = r(k, j) + projection
                    v = v - projection * q(:, k)
                end do
            end do
            norm_v = sqrt(dot_product(v, v))
            if (.not. ieee_is_finite(norm_v)) return
            if (norm_v <= 64.0_dp * epsilon(norm_v) * sqrt(real(window_length, dp))) return
            r(j, j) = norm_v
            q(:, j) = v / norm_v
        end do

        rhs = 0.0_dp
        if (deriv <= polyorder) then
            rhs(deriv + 1) = 1.0_dp
            do i = 1, deriv
                rhs(deriv + 1) = rhs(deriv + 1) * real(i, dp) / (delta * scale)
            end do
        end if
        if (.not. all(ieee_is_finite(rhs))) return

        ! If A = Q R, the minimum-norm weights solve A^T w = rhs:
        ! R^T z = rhs, w = Q z.  All matrix dimensions are window x (degree+1).
        do j = 1, p
            z(j) = rhs(j)
            do k = 1, j - 1
                z(j) = z(j) - r(k, j) * z(k)
            end do
            z(j) = z(j) / r(j, j)
        end do
        if (.not. all(ieee_is_finite(z))) return
        w = matmul(q, z)
        if (.not. all(ieee_is_finite(w))) then
            w = nanval
            return
        end if
        ok = .true.
    end function polynomial_weights

    function savgol_coeffs(window_length, polyorder, deriv, delta, pos, use, status) result(coeffs)
        integer, intent(in) :: window_length !! number of samples, at least one
        integer, intent(in) :: polyorder !! fitting polynomial degree, [0, window_length)
        integer, intent(in), optional :: deriv !! derivative order, default zero
        real(dp), intent(in), optional :: delta !! sampling interval, default 1
        integer, intent(in), optional :: pos !! one-based evaluation position, default center
        character(len=*), intent(in), optional :: use !! 'conv' (default) or 'dot'
        integer, intent(out), optional :: status !! zero on success; nonzero on error
        real(dp), allocatable :: coeffs(:)

        integer :: d, ierr
        real(dp) :: spacing, position, nanval
        logical :: ok
        character(len=:), allocatable :: ordering

        if (present(status)) status = savgol_status_invalid_input
        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(coeffs(max(0, window_length)))
        coeffs = nanval
        if (window_length < 1) return
        if (polyorder < 0 .or. polyorder >= window_length) return
        d = 0
        if (present(deriv)) d = deriv
        if (d < 0) return
        spacing = 1.0_dp
        if (present(delta)) spacing = delta
        if (.not. ieee_is_finite(spacing)) return
        if (abs(spacing) < tiny(spacing)) return
        ordering = 'conv'
        if (present(use)) ordering = use
        if (ordering /= 'conv' .and. ordering /= 'dot') return
        position = 0.5_dp * real(window_length - 1, dp)
        if (present(pos)) then
            if (pos < 1 .or. pos > window_length) return
            position = real(pos - 1, dp)
        end if
        deallocate(coeffs)
        coeffs = polynomial_weights(window_length, polyorder, d, spacing, position, ok)
        if (.not. ok) then
            ierr = savgol_status_singular
            coeffs = nanval
        else
            ierr = savgol_status_success
            if (ordering == 'conv') coeffs = coeffs(window_length:1:-1)
        end if
        if (present(status)) status = ierr
    end function savgol_coeffs

    function make_savgol_kernel(n, window_length, polyorder, deriv, delta, mode, status) result(kernel)
        integer, intent(in) :: n !! number of samples on the filtered axis
        integer, intent(in) :: window_length !! length of the polynomial window
        integer, intent(in) :: polyorder !! degree of the fitted polynomial
        integer, intent(in), optional :: deriv !! requested derivative order
        real(dp), intent(in), optional :: delta !! spacing between samples
        character(len=*), intent(in), optional :: mode !! signal extension mode
        integer, intent(out) :: status !! zero on success
        type(savgol_kernel) :: kernel
        integer :: d, i, code, start
        real(dp), allocatable :: weights(:)
        real(dp) :: spacing
        logical :: ok

        status = savgol_status_invalid_input
        kernel%n = n
        kernel%window = window_length
        kernel%half = max(0,window_length / 2)
        kernel%boundary = 'interp'
        if (present(mode)) kernel%boundary = mode
        select case (kernel%boundary)
        case ('interp','mirror','nearest','wrap','constant')
        case default
            return
        end select
        if (kernel%boundary == 'interp' .and. window_length > n) return
        d = 0
        if (present(deriv)) d = deriv
        spacing = 1.0_dp
        if (present(delta)) spacing = delta
        weights = savgol_coeffs(window_length, polyorder, deriv=d, delta=spacing, status=code)
        if (code /= savgol_status_success) then
            status = code
            return
        end if
        call move_alloc(weights,kernel%center)
        allocate(kernel%left(0,0),kernel%right(0,0))
        if (kernel%boundary == 'interp' .and. kernel%half > 0) then
            deallocate(kernel%left,kernel%right)
            allocate(kernel%left(window_length,kernel%half))
            allocate(kernel%right(window_length,kernel%half))
            start = n - window_length + 1
            do i = 1, kernel%half
                kernel%left(:,i) = polynomial_weights(window_length, polyorder, d, spacing, &
                    real(i-1,dp), ok)
                if (.not. ok) then
                    status = savgol_status_singular
                    return
                end if
                kernel%right(:,i) = polynomial_weights(window_length, polyorder, d, spacing, &
                    real(n-kernel%half+i-start,dp), ok)
                if (.not. ok) then
                    status = savgol_status_singular
                    return
                end if
            end do
        end if
        status = savgol_status_success
    end function make_savgol_kernel

    subroutine apply_savgol_real(x, kernel, exterior, y)
        real(dp), intent(in) :: x(:) !! input channel, length kernel%n
        type(savgol_kernel), intent(in) :: kernel !! precomputed polynomial weights
        real(dp), intent(in) :: exterior !! constant extension value
        real(dp), intent(out) :: y(:) !! filtered output with the same length as x
        real(dp) :: value
        integer :: i, j, idx, n, nw, half

        n = size(x)
        nw = kernel%window
        half = kernel%half
        do i = 1, n
            if (kernel%boundary == 'interp') then
                if (i <= half) then
                    y(i) = dot_product(kernel%left(:,i),x(1:nw))
                    cycle
                else if (i > n-half) then
                    y(i) = dot_product(kernel%right(:,i-(n-half)),x(n-nw+1:n))
                    cycle
                end if
            end if
            value = 0.0_dp
            do j = 1, nw
                idx = i + half + 1 - j
                if (idx >= 1 .and. idx <= n) then
                    value = value + kernel%center(j)*x(idx)
                else
                    select case (kernel%boundary)
                    case ('nearest')
                        value = value + kernel%center(j)*x(max(1,min(n,idx)))
                    case ('wrap')
                        value = value + kernel%center(j)*x(modulo(idx-1,n)+1)
                    case ('mirror')
                        if (n == 1) then
                            idx = 1
                        else
                            idx = modulo(idx-1,2*(n-1))
                            if (idx >= n) idx = 2*(n-1)-idx
                            idx = idx+1
                        end if
                        value = value + kernel%center(j)*x(idx)
                    case ('constant')
                        value = value + kernel%center(j)*exterior
                    end select
                end if
            end do
            y(i) = value
        end do
    end subroutine apply_savgol_real

    subroutine apply_savgol_complex(x, kernel, exterior, y)
        complex(dp), intent(in) :: x(:) !! complex input channel
        type(savgol_kernel), intent(in) :: kernel !! precomputed polynomial weights
        complex(dp), intent(in) :: exterior !! complex constant extension
        complex(dp), intent(out) :: y(:) !! filtered output with the same length as x
        complex(dp) :: value
        integer :: i, j, idx, n, nw, half

        n = size(x)
        nw = kernel%window
        half = kernel%half
        do i = 1, n
            if (kernel%boundary == 'interp') then
                if (i <= half) then
                    y(i) = sum(kernel%left(:,i)*x(1:nw))
                    cycle
                else if (i > n-half) then
                    y(i) = sum(kernel%right(:,i-(n-half))*x(n-nw+1:n))
                    cycle
                end if
            end if
            value = cmplx(0.0_dp,0.0_dp,dp)
            do j = 1, nw
                idx = i + half + 1 - j
                if (idx >= 1 .and. idx <= n) then
                    value = value + kernel%center(j)*x(idx)
                else
                    select case (kernel%boundary)
                    case ('nearest')
                        value = value + kernel%center(j)*x(max(1,min(n,idx)))
                    case ('wrap')
                        value = value + kernel%center(j)*x(modulo(idx-1,n)+1)
                    case ('mirror')
                        if (n == 1) then
                            idx = 1
                        else
                            idx = modulo(idx-1,2*(n-1))
                            if (idx >= n) idx = 2*(n-1)-idx
                            idx = idx+1
                        end if
                        value = value + kernel%center(j)*x(idx)
                    case ('constant')
                        value = value + kernel%center(j)*exterior
                    end select
                end if
            end do
            y(i) = value
        end do
    end subroutine apply_savgol_complex


    ! A contiguous assumed-rank dummy is the standard Fortran 2018 bridge
    ! between arbitrary array rank and a single column-major linear buffer.
    ! Noncontiguous actual arguments are handled by the compiler's copy-in/
    ! copy-out mechanism. The C interoperability intrinsics only remap the
    ! already-contiguous storage; no C code or external library is used.
    ! dp corresponds to the interoperable C double kind on supported targets.
    subroutine savgol_filter_nd_real(x, y, window_length, polyorder, deriv, delta, mode, cval, status, axis)
        real(dp), intent(in), target, contiguous :: x(..) !! input array, ranks 1 through 15
        real(dp), intent(out), target, contiguous :: y(..) !! output, same shape as x
        integer, intent(in) :: window_length !! number of samples in a fitting window
        integer, intent(in) :: polyorder !! fitting polynomial degree
        integer, intent(in), optional :: deriv !! derivative order, default zero
        real(dp), intent(in), optional :: delta !! sample spacing, default one
        character(len=*), intent(in), optional :: mode !! interp, mirror, nearest, wrap or constant
        real(dp), intent(in), optional :: cval !! constant extension, default zero
        integer, intent(out), optional :: status !! success or invalid/singular status
        integer, intent(in), optional :: axis !! one-based axis, default last
        type(savgol_kernel) :: kernel
        real(dp), pointer :: input(:) => null(), output(:) => null()
        real(dp) :: exterior, nanval
        integer :: code, ax, dimensions(rank(x))
        integer :: n
        integer(int64) :: stride, block, outer, i, j, first, last

        if (present(status)) status = savgol_status_invalid_input
        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(y, kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y), output, [size(y, kind=int64)])
            output = nanval
        end if
        if (rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        dimensions = shape(x)
        n = dimensions(ax)
        kernel = make_savgol_kernel(n, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, status=code)
        if (code /= savgol_status_success) then
            if (present(status)) status = code
            return
        end if
        if (size(x, kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(x), input, [size(x, kind=int64)])
            exterior = 0.0_dp
            if (present(cval)) exterior = cval
            stride = product(int(dimensions(:ax-1), int64))
            block = stride*n
            outer = product(int(dimensions(ax+1:), int64))
            do j = 0_int64, outer-1_int64
                do i = 0_int64, stride-1_int64
                    first = j*block+i+1
                    last = first+(n-1)*stride
                    call apply_savgol_real(input(first:last:stride), kernel, exterior, &
                        output(first:last:stride))
                end do
            end do
        end if
        if (present(status)) status = savgol_status_success
    end subroutine savgol_filter_nd_real

    subroutine savgol_filter_nd_complex(x, y, window_length, polyorder, deriv, delta, mode, cval, status, axis)
        complex(dp), intent(in), target, contiguous :: x(..) !! complex input, ranks 1 through 15
        complex(dp), intent(out), target, contiguous :: y(..) !! output, same shape as x
        integer, intent(in) :: window_length !! number of samples in a fitting window
        integer, intent(in) :: polyorder !! fitting polynomial degree
        integer, intent(in), optional :: deriv !! derivative order, default zero
        real(dp), intent(in), optional :: delta !! sample spacing, default one
        character(len=*), intent(in), optional :: mode !! interp, mirror, nearest, wrap or constant
        complex(dp), intent(in), optional :: cval !! constant extension, default zero
        integer, intent(out), optional :: status !! success or invalid/singular status
        integer, intent(in), optional :: axis !! one-based axis, default last
        type(savgol_kernel) :: kernel
        complex(dp), pointer :: input(:) => null(), output(:) => null()
        complex(dp) :: exterior
        real(dp) :: nanval
        integer :: code, ax, dimensions(rank(x))
        integer :: n
        integer(int64) :: stride, block, outer, i, j, first, last

        if (present(status)) status = savgol_status_invalid_input
        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(y, kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y), output, [size(y, kind=int64)])
            output = cmplx(nanval, nanval, dp)
        end if
        if (rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        dimensions = shape(x)
        n = dimensions(ax)
        kernel = make_savgol_kernel(n, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, status=code)
        if (code /= savgol_status_success) then
            if (present(status)) status = code
            return
        end if
        if (size(x, kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(x), input, [size(x, kind=int64)])
            exterior = cmplx(0.0_dp, 0.0_dp, dp)
            if (present(cval)) exterior = cval
            stride = product(int(dimensions(:ax-1), int64))
            block = stride*n
            outer = product(int(dimensions(ax+1:), int64))
            do j = 0_int64, outer-1_int64
                do i = 0_int64, stride-1_int64
                    first = j*block+i+1
                    last = first+(n-1)*stride
                    call apply_savgol_complex(input(first:last:stride), kernel, exterior, &
                        output(first:last:stride))
                end do
            end do
        end if
        if (present(status)) status = savgol_status_success
    end subroutine savgol_filter_nd_complex

    function savgol_filter_real_1d(x, window_length, polyorder, deriv, delta, mode, cval, status) result(y)
        real(dp), intent(in) :: x(:) !! rank-1 input array
        integer, intent(in) :: window_length !! number of samples in the fitting window
        integer, intent(in) :: polyorder !! polynomial degree
        integer, intent(in), optional :: deriv !! derivative order
        real(dp), intent(in), optional :: delta !! sample spacing
        character(len=*), intent(in), optional :: mode !! extension mode
        real(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)

        allocate(y(size(x,1)))
        call savgol_filter_nd(x, y, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, cval=cval, status=status)
    end function savgol_filter_real_1d

    function savgol_filter_real_2d(x, window_length, polyorder, deriv, delta, mode, cval, status, axis) result(y)
        real(dp), intent(in) :: x(:,:) !! rank-2 input array
        integer, intent(in) :: window_length !! number of samples in the fitting window
        integer, intent(in) :: polyorder !! polynomial degree
        integer, intent(in), optional :: deriv !! derivative order
        real(dp), intent(in), optional :: delta !! sample spacing
        character(len=*), intent(in), optional :: mode !! extension mode
        real(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! one-based axis, default 2
        real(dp), allocatable :: y(:,:)

        allocate(y(size(x,1),size(x,2)))
        call savgol_filter_nd(x, y, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, cval=cval, status=status, axis=axis)
    end function savgol_filter_real_2d

    function savgol_filter_real_3d(x, window_length, polyorder, deriv, delta, mode, cval, status, axis) result(y)
        real(dp), intent(in) :: x(:,:,:) !! rank-3 input array
        integer, intent(in) :: window_length !! number of samples in the fitting window
        integer, intent(in) :: polyorder !! polynomial degree
        integer, intent(in), optional :: deriv !! derivative order
        real(dp), intent(in), optional :: delta !! sample spacing
        character(len=*), intent(in), optional :: mode !! extension mode
        real(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! one-based axis, default 3
        real(dp), allocatable :: y(:,:,:)

        allocate(y(size(x,1),size(x,2),size(x,3)))
        call savgol_filter_nd(x, y, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, cval=cval, status=status, axis=axis)
    end function savgol_filter_real_3d

    function savgol_filter_real_4d(x, window_length, polyorder, deriv, delta, mode, cval, status, axis) result(y)
        real(dp), intent(in) :: x(:,:,:,:) !! rank-4 input array
        integer, intent(in) :: window_length !! number of samples in the fitting window
        integer, intent(in) :: polyorder !! polynomial degree
        integer, intent(in), optional :: deriv !! derivative order
        real(dp), intent(in), optional :: delta !! sample spacing
        character(len=*), intent(in), optional :: mode !! extension mode
        real(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! one-based axis, default 4
        real(dp), allocatable :: y(:,:,:,:)

        allocate(y(size(x,1),size(x,2),size(x,3),size(x,4)))
        call savgol_filter_nd(x, y, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, cval=cval, status=status, axis=axis)
    end function savgol_filter_real_4d

    function savgol_filter_complex_1d(x, window_length, polyorder, deriv, delta, mode, cval, status) result(y)
        complex(dp), intent(in) :: x(:) !! rank-1 input array
        integer, intent(in) :: window_length !! number of samples in the fitting window
        integer, intent(in) :: polyorder !! polynomial degree
        integer, intent(in), optional :: deriv !! derivative order
        real(dp), intent(in), optional :: delta !! sample spacing
        character(len=*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)

        allocate(y(size(x,1)))
        call savgol_filter_nd(x, y, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, cval=cval, status=status)
    end function savgol_filter_complex_1d

    function savgol_filter_complex_2d(x, window_length, polyorder, deriv, delta, mode, cval, status, axis) result(y)
        complex(dp), intent(in) :: x(:,:) !! rank-2 input array
        integer, intent(in) :: window_length !! number of samples in the fitting window
        integer, intent(in) :: polyorder !! polynomial degree
        integer, intent(in), optional :: deriv !! derivative order
        real(dp), intent(in), optional :: delta !! sample spacing
        character(len=*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! one-based axis, default 2
        complex(dp), allocatable :: y(:,:)

        allocate(y(size(x,1),size(x,2)))
        call savgol_filter_nd(x, y, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, cval=cval, status=status, axis=axis)
    end function savgol_filter_complex_2d

    function savgol_filter_complex_3d(x, window_length, polyorder, deriv, delta, mode, cval, status, axis) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! rank-3 input array
        integer, intent(in) :: window_length !! number of samples in the fitting window
        integer, intent(in) :: polyorder !! polynomial degree
        integer, intent(in), optional :: deriv !! derivative order
        real(dp), intent(in), optional :: delta !! sample spacing
        character(len=*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! one-based axis, default 3
        complex(dp), allocatable :: y(:,:,:)

        allocate(y(size(x,1),size(x,2),size(x,3)))
        call savgol_filter_nd(x, y, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, cval=cval, status=status, axis=axis)
    end function savgol_filter_complex_3d

    function savgol_filter_complex_4d(x, window_length, polyorder, deriv, delta, mode, cval, status, axis) result(y)
        complex(dp), intent(in) :: x(:,:,:,:) !! rank-4 input array
        integer, intent(in) :: window_length !! number of samples in the fitting window
        integer, intent(in) :: polyorder !! polynomial degree
        integer, intent(in), optional :: deriv !! derivative order
        real(dp), intent(in), optional :: delta !! sample spacing
        character(len=*), intent(in), optional :: mode !! extension mode
        complex(dp), intent(in), optional :: cval !! constant extension value
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! one-based axis, default 4
        complex(dp), allocatable :: y(:,:,:,:)

        allocate(y(size(x,1),size(x,2),size(x,3),size(x,4)))
        call savgol_filter_nd(x, y, window_length, polyorder, &
            deriv=deriv, delta=delta, mode=mode, cval=cval, status=status, axis=axis)
    end function savgol_filter_complex_4d

    ! Local-polynomial smoothing/derivatives for irregular observation times.
    ! This is a SciFort extension: scipy.signal.savgol_filter assumes regular spacing.
    ! Each output uses the nearest consecutive window_length observations and
    ! a weighted, centered Vandermonde QR; no normal equations are formed.
    pure integer function irregular_window_start(x, at, count) result(first)
        real(dp), intent(in) :: x(:) !! strictly increasing observation times
        real(dp), intent(in) :: at !! evaluation time
        integer, intent(in) :: count !! number of nearest observations
        integer :: lo, hi, mid, left, right, j, n

        n = size(x)
        lo = 0
        hi = n + 1
        do while (hi - lo > 1)
            mid = lo + (hi - lo) / 2
            if (x(mid) <= at) then
                lo = mid
            else
                hi = mid
            end if
        end do
        left = lo
        right = hi
        do j = 1, count
            if (left < 1) then
                right = right + 1
            else if (right > n) then
                left = left - 1
            else if (at - x(left) <= x(right) - at) then
                left = left - 1
            else
                right = right + 1
            end if
        end do
        first = left + 1
    end function irregular_window_start

    function irregular_savgol_weights(x, weights, at, polyorder, deriv, ok) result(coeff)
        real(dp), intent(in) :: x(:) !! observation times for one local window
        real(dp), intent(in) :: weights(:) !! nonnegative precision weights for this window
        real(dp), intent(in) :: at !! requested evaluation location
        integer, intent(in) :: polyorder !! fitted polynomial degree
        integer, intent(in) :: deriv !! derivative order, nonnegative
        logical, intent(out) :: ok !! true when the weighted design has full rank
        real(dp), allocatable :: coeff(:)
        real(dp), allocatable :: q(:,:), r(:,:), v(:), root_w(:), rhs(:), z(:)
        real(dp) :: scale, largest, projection, norm_v, test_scale
        integer :: n, p, i, j, k, pass

        ok = .false.
        n = size(x)
        p = polyorder + 1
        allocate(coeff(n))
        coeff = 0.0_dp
        largest = maxval(weights)
        if (largest <= 0.0_dp) return
        scale = maxval(abs(x - at))
        if (scale <= 0.0_dp .or. .not. ieee_is_finite(scale)) return
        allocate(root_w(n), q(n,p), r(p,p), v(n), rhs(p), z(p))
        root_w = sqrt(weights / largest)
        q = 0.0_dp
        r = 0.0_dp
        do j = 1, p
            do i = 1, n
                v(i) = root_w(i) * ((x(i) - at) / scale)**(j - 1)
            end do
            test_scale = sqrt(dot_product(v, v))
            do pass = 1, 2
                do k = 1, j - 1
                    projection = dot_product(q(:,k), v)
                    r(k,j) = r(k,j) + projection
                    v = v - projection * q(:,k)
                end do
            end do
            norm_v = sqrt(dot_product(v, v))
            if (.not. ieee_is_finite(norm_v)) return
            if (norm_v <= 128.0_dp * epsilon(norm_v) * test_scale) return
            r(j,j) = norm_v
            q(:,j) = v / norm_v
        end do
        rhs = 0.0_dp
        if (deriv <= polyorder) then
            rhs(deriv+1) = 1.0_dp
            do i = 1, deriv
                rhs(deriv+1) = rhs(deriv+1) * real(i,dp) / scale
            end do
        end if
        if (.not. all(ieee_is_finite(rhs))) return
        ! A=QR, A^T (sqrt(W)*Q*z) = rhs => R^T*z=rhs.
        do j = 1, p
            z(j) = rhs(j)
            do k = 1, j - 1
                z(j) = z(j) - r(k,j) * z(k)
            end do
            z(j) = z(j) / r(j,j)
        end do
        coeff = root_w * matmul(q,z)
        if (.not. all(ieee_is_finite(coeff))) return
        ok = .true.
    end function irregular_savgol_weights

    function savgol_irregular_real(x, y, window_length, polyorder, deriv, weights, query, status) result(out)
        real(dp), intent(in) :: x(:) !! strictly increasing, possibly nonuniform observation times
        real(dp), intent(in) :: y(:) !! observed real signal, same length as x
        integer, intent(in) :: window_length !! number of nearest neighbors in each polynomial fit
        integer, intent(in) :: polyorder !! polynomial degree in [0,window_length)
        integer, intent(in), optional :: deriv !! derivative with respect to time, default zero
        real(dp), intent(in), optional :: weights(:) !! nonnegative relative precision weights
        real(dp), intent(in), optional :: query(:) !! output times, default x(:)
        integer, intent(out), optional :: status !! success, invalid input, or singular local fit
        real(dp), allocatable :: out(:)
        real(dp), allocatable :: weight_values(:), locations(:), coeff(:)
        real(dp) :: nanval
        integer :: n, nout, d, j, first
        logical :: ok

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        n = size(x)
        nout = n
        if (present(query)) nout = size(query)
        allocate(out(nout))
        out = nanval
        if (present(status)) status = savgol_status_invalid_input
        d = 0
        if (present(deriv)) d = deriv
        if (window_length < 2 .or. window_length > n) return
        if (polyorder < 0 .or. polyorder >= window_length .or. d < 0) return
        if (n /= size(y)) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) return
        if (any(x(2:) <= x(:n-1))) return
        allocate(weight_values(n), locations(nout))
        weight_values = 1.0_dp
        if (present(weights)) then
            if (size(weights) /= n) return
            if (.not. all(ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            weight_values = weights
        end if
        if (present(query)) then
            locations = query
        else
            locations = x
        end if
        if (.not. all(ieee_is_finite(locations))) return
        do j = 1, nout
            first = irregular_window_start(x, locations(j), window_length)
            coeff = irregular_savgol_weights(x(first:first+window_length-1), &
                weight_values(first:first+window_length-1), locations(j), polyorder, d, ok)
            if (.not. ok) then
                out = nanval
                if (present(status)) status = savgol_status_singular
                return
            end if
            out(j) = dot_product(coeff, y(first:first+window_length-1))
            if (.not. ieee_is_finite(out(j))) then
                out = nanval
                if (present(status)) status = savgol_status_singular
                return
            end if
        end do
        if (present(status)) status = savgol_status_success
    end function savgol_irregular_real

    function savgol_irregular_complex(x, y, window_length, polyorder, deriv, weights, query, status) result(out)
        real(dp), intent(in) :: x(:) !! strictly increasing, possibly nonuniform observation times
        complex(dp), intent(in) :: y(:) !! complex observations, same length as x
        integer, intent(in) :: window_length !! number of nearest neighbors in each polynomial fit
        integer, intent(in) :: polyorder !! polynomial degree in [0,window_length)
        integer, intent(in), optional :: deriv !! derivative with respect to time, default zero
        real(dp), intent(in), optional :: weights(:) !! nonnegative relative precision weights
        real(dp), intent(in), optional :: query(:) !! output times, default x(:)
        integer, intent(out), optional :: status !! success, invalid input, or singular local fit
        complex(dp), allocatable :: out(:)
        real(dp), allocatable :: weight_values(:), locations(:), coeff(:)
        real(dp) :: nanval
        integer :: n, nout, d, j, first
        logical :: ok

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        n = size(x)
        nout = n
        if (present(query)) nout = size(query)
        allocate(out(nout))
        out = cmplx(nanval,nanval,dp)
        if (present(status)) status = savgol_status_invalid_input
        d = 0
        if (present(deriv)) d = deriv
        if (window_length < 2 .or. window_length > n) return
        if (polyorder < 0 .or. polyorder >= window_length .or. d < 0) return
        if (n /= size(y)) return
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(real(y,dp))) .or. .not. all(ieee_is_finite(aimag(y)))) return
        if (any(x(2:) <= x(:n-1))) return
        allocate(weight_values(n), locations(nout))
        weight_values = 1.0_dp
        if (present(weights)) then
            if (size(weights) /= n) return
            if (.not. all(ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            weight_values = weights
        end if
        if (present(query)) then
            locations = query
        else
            locations = x
        end if
        if (.not. all(ieee_is_finite(locations))) return
        do j = 1, nout
            first = irregular_window_start(x, locations(j), window_length)
            coeff = irregular_savgol_weights(x(first:first+window_length-1), &
                weight_values(first:first+window_length-1), locations(j), polyorder, d, ok)
            if (.not. ok) then
                out = cmplx(nanval,nanval,dp)
                if (present(status)) status = savgol_status_singular
                return
            end if
            out(j) = sum(coeff * y(first:first+window_length-1))
            if (.not. ieee_is_finite(real(out(j),dp)) .or. .not. ieee_is_finite(aimag(out(j)))) then
                out = cmplx(nanval,nanval,dp)
                if (present(status)) status = savgol_status_singular
                return
            end if
        end do
        if (present(status)) status = savgol_status_success
    end function savgol_irregular_complex

end module scifort_signal_savgol
