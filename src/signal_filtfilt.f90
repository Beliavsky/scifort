! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! 1D zero-phase forward/backward filtering with the SciPy padding convention.
! Supports both reflection padding and Gustafsson initial conditions.
module scifort_signal_filtfilt
    use scifort_kinds, only : dp
    use scifort_signal_filtfilt_gust, only : filtfilt_gust, sosfiltfilt_gust, &
        gust_status_success, gust_status_invalid_input, gust_status_numerical_failure
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use scifort_signal_lfilter, only : lfilter, lfilter_zi, lfilter_status_success, &
        lfilter_status_no_steady_state
    use scifort_signal_sos, only : sosfilt, sosfilt_zi, sos_status_success, &
        sos_status_no_steady_state
    implicit none
    private
    integer, parameter, public :: filtfilt_status_success = 0
    integer, parameter, public :: filtfilt_status_invalid_input = 1
    integer, parameter, public :: filtfilt_status_no_steady_state = 2
    integer, parameter, public :: filtfilt_status_numerical_failure = 3
    interface filtfilt
        module procedure filtfilt_real
        module procedure filtfilt_complex
        module procedure filtfilt_matrix_real, filtfilt_matrix_complex
        module procedure filtfilt_rank3_real, filtfilt_rank3_complex
    end interface
    interface sosfiltfilt
        module procedure sosfiltfilt_real
        module procedure sosfiltfilt_complex
        module procedure sosfiltfilt_matrix_real, sosfiltfilt_matrix_complex
        module procedure sosfiltfilt_rank3_real, sosfiltfilt_rank3_complex
    end interface
    public :: filtfilt, sosfiltfilt
    public :: valid_sos_gust_real, valid_sos_gust_complex
contains
    pure integer function pad_width(n, ntaps, padtype, padlen) result(edge)
        integer, intent(in) :: n !! number of samples
        integer, intent(in) :: ntaps !! nominal number of filter taps
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, or none
        integer, intent(in), optional :: padlen !! number of samples padded at each edge
        character(len=:), allocatable :: mode
        mode = 'odd'
        if (present(padtype)) mode = trim(padtype)
        edge = -1
        select case (mode)
        case ('odd', 'even', 'constant', 'none')
        case default
            return
        end select
        edge = 3 * ntaps
        if (present(padlen)) edge = padlen
        if (mode == 'none') edge = 0
        if (edge < 0 .or. n <= edge) edge = -1
    end function


    pure logical function valid_sos_gust_real(sos,n,irlen) result(ok)
        real(dp), intent(in) :: sos(:,:) !! SOS matrix
        integer, intent(in) :: n !! selected dimension length
        integer, optional, intent(in) :: irlen !! optional positive impulse response length
        ok = .false.
        if (size(sos,1) < 1 .or. size(sos,2) /= 6 .or. n < 1) return
        if (present(irlen)) then
            if (irlen < 1) return
        end if
        if (any(sos(:,4) /= 1.0_dp)) return
        if (.not. all(ieee_is_finite(sos))) return
        ok = .true.
    end function valid_sos_gust_real

    pure logical function valid_sos_gust_complex(sos,n,irlen) result(ok)
        complex(dp), intent(in) :: sos(:,:) !! SOS matrix
        integer, intent(in) :: n !! selected dimension length
        integer, optional, intent(in) :: irlen !! optional positive impulse response length
        ok = .false.
        if (size(sos,1) < 1 .or. size(sos,2) /= 6 .or. n < 1) return
        if (present(irlen)) then
            if (irlen < 1) return
        end if
        if (any(sos(:,4) /= (1.0_dp,0.0_dp))) return
        if (.not. all(ieee_is_finite(real(sos,dp)))) return
        if (.not. all(ieee_is_finite(aimag(sos)))) return
        ok = .true.
    end function valid_sos_gust_complex

    function extend_real(x, edge, padtype) result(out)
        real(dp), intent(in) :: x(:) !! signal to extend
        integer, intent(in) :: edge !! samples added at each end
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, or none
        real(dp), allocatable :: out(:)
        character(len=:), allocatable :: mode
        integer :: n, i

        n = size(x)
        mode = 'odd'
        if (present(padtype)) mode = trim(padtype)
        allocate(out(n + 2*edge))
        out(edge+1:edge+n) = x
        do i = 1, edge
            select case (mode)
            case ('even')
                out(i) = x(edge-i+2)
                out(edge+n+i) = x(n-i)
            case ('constant')
                out(i) = x(1)
                out(edge+n+i) = x(n)
            case default
                out(i) = 2*x(1) - x(edge-i+2)
                out(edge+n+i) = 2*x(n) - x(n-i)
            end select
        end do
    end function extend_real

    function filtfilt_real(b, a, x, padtype, padlen, status, method, irlen) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients, ascending lag
        real(dp), intent(in) :: a(:) !! denominator coefficients, ascending lag
        real(dp), intent(in) :: x(:) !! input samples
        character(len=*), intent(in), optional :: padtype !! odd (default), even, constant, none
        integer, intent(in), optional :: padlen !! padding length, default 3*max(size(a),size(b))
        character(len=*), intent(in), optional :: method !! 'pad' (default) or 'gust'
        integer, intent(in), optional :: irlen !! impulse-response length for Gustafsson method
        integer, intent(out), optional :: status !! zero on success
        character(:), allocatable :: mode
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: zi(:), ext(:), first(:), second(:)
        real(dp) :: nanval
        integer :: edge, ls, n

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        n = size(x)
        allocate(y(n))
        y = nanval
        if (present(status)) status = filtfilt_status_invalid_input
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            y = filtfilt_gust(b, a, x, irlen=irlen, status=ls)
            if (present(status)) then
                status = filtfilt_status_invalid_input
                if (ls == gust_status_success) status = filtfilt_status_success
                if (ls == gust_status_numerical_failure) status = filtfilt_status_numerical_failure
            end if
            return
        end if
        if (mode /= 'pad') return
        if (size(a) < 1 .or. size(b) < 1) return
        edge = pad_width(n, max(size(a),size(b)), padtype, padlen)
        if (edge < 0) return
        zi = lfilter_zi(b, a, status=ls)
        if (ls /= lfilter_status_success) then
            if (present(status)) then
                if (ls == lfilter_status_no_steady_state) status = filtfilt_status_no_steady_state
            end if
            return
        end if
        ext = extend_real(x, edge, padtype)
        first = lfilter(b, a, ext, zi=zi*ext(1), status=ls)
        if (ls /= lfilter_status_success) return
        second = lfilter(b, a, first(size(first):1:-1), zi=zi*first(size(first)), status=ls)
        if (ls /= lfilter_status_success) return
        y = second(size(second)-edge:edge+1:-1)
        if (present(status)) status = filtfilt_status_success
    end function filtfilt_real

    function sosfiltfilt_real(sos, x, padtype, padlen, status, method, irlen) result(y)
        real(dp), intent(in) :: sos(:,:) !! SOS rows (nsections,6), with a0=1
        real(dp), intent(in) :: x(:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd (default), even, constant, none
        integer, intent(in), optional :: padlen !! optional padding length
        integer, intent(out), optional :: status !! zero on success
        character(len=*), intent(in), optional :: method !! pad (default) or gust
        integer, intent(in), optional :: irlen !! positive Gustafsson truncation length
        character(:), allocatable :: mode
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: zi(:,:), ext(:), first(:), second(:)
        real(dp) :: nanval
        integer :: edge, ntaps, n, ls

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        n = size(x)
        allocate(y(n))
        y = nanval
        if (present(status)) status = filtfilt_status_invalid_input
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            y = sosfiltfilt_gust(sos, x, irlen=irlen, status=ls)
            if (present(status)) then
                status = filtfilt_status_invalid_input
                if (ls == gust_status_success) status = filtfilt_status_success
                if (ls == gust_status_numerical_failure) status = filtfilt_status_numerical_failure
            end if
            return
        end if
        if (mode /= 'pad') return
        if (size(sos,2) /= 6 .or. size(sos,1) < 1) return
        ntaps = 2*size(sos,1) + 1 - min(count(sos(:,3) == 0.0_dp), &
            count(sos(:,6) == 0.0_dp))
        edge = pad_width(n, ntaps, padtype, padlen)
        if (edge < 0) return
        zi = sosfilt_zi(sos, status=ls)
        if (ls /= sos_status_success) then
            if (present(status)) then
                if (ls == sos_status_no_steady_state) status = filtfilt_status_no_steady_state
            end if
            return
        end if
        ext = extend_real(x, edge, padtype)
        first = sosfilt(sos, ext, zi=zi*ext(1), status=ls)
        if (ls /= sos_status_success) return
        second = sosfilt(sos, first(size(first):1:-1), zi=zi*first(size(first)), status=ls)
        if (ls /= sos_status_success) return
        y = second(size(second)-edge:edge+1:-1)
        if (present(status)) status = filtfilt_status_success
    end function sosfiltfilt_real

    function extend_complex(x, edge, padtype) result(out)
        complex(dp), intent(in) :: x(:) !! signal to extend
        integer, intent(in) :: edge !! samples added at each end
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, or none
        complex(dp), allocatable :: out(:)
        character(len=:), allocatable :: mode
        integer :: n, i

        n = size(x)
        mode = 'odd'
        if (present(padtype)) mode = trim(padtype)
        allocate(out(n + 2*edge))
        out(edge+1:edge+n) = x
        do i = 1, edge
            select case (mode)
            case ('even')
                out(i) = x(edge-i+2)
                out(edge+n+i) = x(n-i)
            case ('constant')
                out(i) = x(1)
                out(edge+n+i) = x(n)
            case default
                out(i) = 2*x(1) - x(edge-i+2)
                out(edge+n+i) = 2*x(n) - x(n-i)
            end select
        end do
    end function extend_complex

    function filtfilt_complex(b, a, x, padtype, padlen, status, method, irlen) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients, ascending lag
        complex(dp), intent(in) :: a(:) !! denominator coefficients, ascending lag
        complex(dp), intent(in) :: x(:) !! input samples
        character(len=*), intent(in), optional :: padtype !! odd (default), even, constant, none
        integer, intent(in), optional :: padlen !! padding length, default 3*max(size(a),size(b))
        character(len=*), intent(in), optional :: method !! 'pad' (default) or 'gust'
        integer, intent(in), optional :: irlen !! impulse-response length for Gustafsson method
        integer, intent(out), optional :: status !! zero on success
        character(:), allocatable :: mode
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: zi(:), ext(:), first(:), second(:)
        real(dp) :: nanval
        integer :: edge, ls, n

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        n = size(x)
        allocate(y(n))
        y = cmplx(nanval,nanval,dp)
        if (present(status)) status = filtfilt_status_invalid_input
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            y = filtfilt_gust(b, a, x, irlen=irlen, status=ls)
            if (present(status)) then
                status = filtfilt_status_invalid_input
                if (ls == gust_status_success) status = filtfilt_status_success
                if (ls == gust_status_numerical_failure) status = filtfilt_status_numerical_failure
            end if
            return
        end if
        if (mode /= 'pad') return
        if (size(a) < 1 .or. size(b) < 1) return
        edge = pad_width(n, max(size(a),size(b)), padtype, padlen)
        if (edge < 0) return
        zi = lfilter_zi(b, a, status=ls)
        if (ls /= lfilter_status_success) then
            if (present(status)) then
                if (ls == lfilter_status_no_steady_state) status = filtfilt_status_no_steady_state
            end if
            return
        end if
        ext = extend_complex(x, edge, padtype)
        first = lfilter(b, a, ext, zi=zi*ext(1), status=ls)
        if (ls /= lfilter_status_success) return
        second = lfilter(b, a, first(size(first):1:-1), zi=zi*first(size(first)), status=ls)
        if (ls /= lfilter_status_success) return
        y = second(size(second)-edge:edge+1:-1)
        if (present(status)) status = filtfilt_status_success
    end function filtfilt_complex

    function sosfiltfilt_complex(sos, x, padtype, padlen, status, method, irlen) result(y)
        complex(dp), intent(in) :: sos(:,:) !! SOS rows (nsections,6), with a0=1
        complex(dp), intent(in) :: x(:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd (default), even, constant, none
        integer, intent(in), optional :: padlen !! optional padding length
        integer, intent(out), optional :: status !! zero on success
        character(len=*), intent(in), optional :: method !! pad (default) or gust
        integer, intent(in), optional :: irlen !! positive Gustafsson truncation length
        character(:), allocatable :: mode
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: zi(:,:), ext(:), first(:), second(:)
        real(dp) :: nanval
        integer :: edge, ntaps, n, ls

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        n = size(x)
        allocate(y(n))
        y = cmplx(nanval,nanval,dp)
        if (present(status)) status = filtfilt_status_invalid_input
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            y = sosfiltfilt_gust(sos, x, irlen=irlen, status=ls)
            if (present(status)) then
                status = filtfilt_status_invalid_input
                if (ls == gust_status_success) status = filtfilt_status_success
                if (ls == gust_status_numerical_failure) status = filtfilt_status_numerical_failure
            end if
            return
        end if
        if (mode /= 'pad') return
        if (size(sos,2) /= 6 .or. size(sos,1) < 1) return
        ntaps = 2*size(sos,1) + 1 - min(count(sos(:,3) == (0.0_dp,0.0_dp)), &
            count(sos(:,6) == (0.0_dp,0.0_dp)))
        edge = pad_width(n, ntaps, padtype, padlen)
        if (edge < 0) return
        zi = sosfilt_zi(sos, status=ls)
        if (ls /= sos_status_success) then
            if (present(status)) then
                if (ls == sos_status_no_steady_state) status = filtfilt_status_no_steady_state
            end if
            return
        end if
        ext = extend_complex(x, edge, padtype)
        first = sosfilt(sos, ext, zi=zi*ext(1), status=ls)
        if (ls /= sos_status_success) return
        second = sosfilt(sos, first(size(first):1:-1), zi=zi*first(size(first)), status=ls)
        if (ls /= sos_status_success) return
        y = second(size(second)-edge:edge+1:-1)
        if (present(status)) status = filtfilt_status_success
    end function sosfiltfilt_complex


    function filtfilt_matrix_real(b, a, x, padtype, padlen, status, axis, method, irlen) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:,:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, none
        integer, intent(in), optional :: padlen !! number of samples padded
        character(len=*), intent(in), optional :: method !! 'pad' (default) or 'gust'
        integer, intent(in), optional :: irlen !! impulse-response length for Gustafsson method
        integer, intent(out), optional :: status !! success or error
        integer, intent(in), optional :: axis !! one-based axis, default 2
        character(:), allocatable :: mode
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        real(dp), allocatable :: steady(:)
        integer :: edge, ntaps, steady_status

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2)))
        y = nanval
        if (present(status)) status = filtfilt_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along < 1 .or. along > 2) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode /= 'pad' .and. mode /= 'gust') return
        if (size(b) < 1 .or. size(a) < 1) return
        if (mode == 'gust') then
            if (size(x,along) < 1) return
            if (abs(a(1)) <= tiny(1.0_dp)) return
            if (present(irlen)) then
                if (irlen < 1) return
            end if
            if (.not. all(ieee_is_finite(b))) return
            if (.not. all(ieee_is_finite(a))) return
        end if
        if (mode == 'pad') then
            ntaps = max(size(b),size(a))
            edge = pad_width(size(x,along), ntaps, padtype, padlen)
            if (edge < 0) return
            steady = lfilter_zi(b,a,status=steady_status)
            if (steady_status /= lfilter_status_success) then
                if (present(status)) then
                    if (steady_status == lfilter_status_no_steady_state) status = filtfilt_status_no_steady_state
                end if
                return
            end if
        end if
        select case (along)
        case (1)
            do i = 1, size(x,2)
                line = filtfilt_real(b, a, x(:,i), padtype=padtype, &
                    padlen=padlen, status=code, method=method, irlen=irlen)
                if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                line = filtfilt_real(b, a, x(i,:), padtype=padtype, &
                    padlen=padlen, status=code, method=method, irlen=irlen)
                if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                y(i,:) = line
            end do
        end select
        if (present(status)) status = filtfilt_status_success
    end function filtfilt_matrix_real

    function filtfilt_matrix_complex(b, a, x, padtype, padlen, status, axis, method, irlen) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:,:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, none
        integer, intent(in), optional :: padlen !! number of samples padded
        character(len=*), intent(in), optional :: method !! 'pad' (default) or 'gust'
        integer, intent(in), optional :: irlen !! impulse-response length for Gustafsson method
        integer, intent(out), optional :: status !! success or error
        integer, intent(in), optional :: axis !! one-based axis, default 2
        character(:), allocatable :: mode
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        complex(dp), allocatable :: steady(:)
        integer :: edge, ntaps, steady_status

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2)))
        y = cmplx(nanval,nanval,dp)
        if (present(status)) status = filtfilt_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along < 1 .or. along > 2) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode /= 'pad' .and. mode /= 'gust') return
        if (size(b) < 1 .or. size(a) < 1) return
        if (mode == 'gust') then
            if (size(x,along) < 1) return
            if (abs(a(1)) <= tiny(1.0_dp)) return
            if (present(irlen)) then
                if (irlen < 1) return
            end if
            if (.not. all(ieee_is_finite(real(b,dp))) .or. &
                .not. all(ieee_is_finite(aimag(b)))) return
            if (.not. all(ieee_is_finite(real(a,dp))) .or. &
                .not. all(ieee_is_finite(aimag(a)))) return
        end if
        if (mode == 'pad') then
            ntaps = max(size(b),size(a))
            edge = pad_width(size(x,along), ntaps, padtype, padlen)
            if (edge < 0) return
            steady = lfilter_zi(b,a,status=steady_status)
            if (steady_status /= lfilter_status_success) then
                if (present(status)) then
                    if (steady_status == lfilter_status_no_steady_state) status = filtfilt_status_no_steady_state
                end if
                return
            end if
        end if
        select case (along)
        case (1)
            do i = 1, size(x,2)
                line = filtfilt_complex(b, a, x(:,i), padtype=padtype, &
                    padlen=padlen, status=code, method=method, irlen=irlen)
                if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                line = filtfilt_complex(b, a, x(i,:), padtype=padtype, &
                    padlen=padlen, status=code, method=method, irlen=irlen)
                if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                y(i,:) = line
            end do
        end select
        if (present(status)) status = filtfilt_status_success
    end function filtfilt_matrix_complex

    function filtfilt_rank3_real(b, a, x, padtype, padlen, status, axis, method, irlen) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:,:,:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, none
        integer, intent(in), optional :: padlen !! number of samples padded
        character(len=*), intent(in), optional :: method !! 'pad' (default) or 'gust'
        integer, intent(in), optional :: irlen !! impulse-response length for Gustafsson method
        integer, intent(out), optional :: status !! success or error
        integer, intent(in), optional :: axis !! one-based axis, default 3
        character(:), allocatable :: mode
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        real(dp), allocatable :: steady(:)
        integer :: edge, ntaps, steady_status

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2),size(x,3)))
        y = nanval
        if (present(status)) status = filtfilt_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 1 .or. along > 3) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode /= 'pad' .and. mode /= 'gust') return
        if (size(b) < 1 .or. size(a) < 1) return
        if (mode == 'gust') then
            if (size(x,along) < 1) return
            if (abs(a(1)) <= tiny(1.0_dp)) return
            if (present(irlen)) then
                if (irlen < 1) return
            end if
            if (.not. all(ieee_is_finite(b))) return
            if (.not. all(ieee_is_finite(a))) return
        end if
        if (mode == 'pad') then
            ntaps = max(size(b),size(a))
            edge = pad_width(size(x,along), ntaps, padtype, padlen)
            if (edge < 0) return
            steady = lfilter_zi(b,a,status=steady_status)
            if (steady_status /= lfilter_status_success) then
                if (present(status)) then
                    if (steady_status == lfilter_status_no_steady_state) status = filtfilt_status_no_steady_state
                end if
                return
            end if
        end if
        select case (along)
        case (1)
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    line = filtfilt_real(b, a, x(:,i,j), padtype=padtype, &
                        padlen=padlen, status=code, method=method, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                    y(:,i,j) = line
                end do
            end do
        case (2)
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    line = filtfilt_real(b, a, x(i,:,j), padtype=padtype, &
                        padlen=padlen, status=code, method=method, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                    y(i,:,j) = line
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    line = filtfilt_real(b, a, x(i,j,:), padtype=padtype, &
                        padlen=padlen, status=code, method=method, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                    y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = filtfilt_status_success
    end function filtfilt_rank3_real

    function filtfilt_rank3_complex(b, a, x, padtype, padlen, status, axis, method, irlen) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:,:,:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, none
        integer, intent(in), optional :: padlen !! number of samples padded
        character(len=*), intent(in), optional :: method !! 'pad' (default) or 'gust'
        integer, intent(in), optional :: irlen !! impulse-response length for Gustafsson method
        integer, intent(out), optional :: status !! success or error
        integer, intent(in), optional :: axis !! one-based axis, default 3
        character(:), allocatable :: mode
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        complex(dp), allocatable :: steady(:)
        integer :: edge, ntaps, steady_status

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2),size(x,3)))
        y = cmplx(nanval,nanval,dp)
        if (present(status)) status = filtfilt_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 1 .or. along > 3) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode /= 'pad' .and. mode /= 'gust') return
        if (size(b) < 1 .or. size(a) < 1) return
        if (mode == 'gust') then
            if (size(x,along) < 1) return
            if (abs(a(1)) <= tiny(1.0_dp)) return
            if (present(irlen)) then
                if (irlen < 1) return
            end if
            if (.not. all(ieee_is_finite(real(b,dp))) .or. &
                .not. all(ieee_is_finite(aimag(b)))) return
            if (.not. all(ieee_is_finite(real(a,dp))) .or. &
                .not. all(ieee_is_finite(aimag(a)))) return
        end if
        if (mode == 'pad') then
            ntaps = max(size(b),size(a))
            edge = pad_width(size(x,along), ntaps, padtype, padlen)
            if (edge < 0) return
            steady = lfilter_zi(b,a,status=steady_status)
            if (steady_status /= lfilter_status_success) then
                if (present(status)) then
                    if (steady_status == lfilter_status_no_steady_state) status = filtfilt_status_no_steady_state
                end if
                return
            end if
        end if
        select case (along)
        case (1)
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    line = filtfilt_complex(b, a, x(:,i,j), padtype=padtype, &
                        padlen=padlen, status=code, method=method, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                    y(:,i,j) = line
                end do
            end do
        case (2)
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    line = filtfilt_complex(b, a, x(i,:,j), padtype=padtype, &
                        padlen=padlen, status=code, method=method, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                    y(i,:,j) = line
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    line = filtfilt_complex(b, a, x(i,j,:), padtype=padtype, &
                        padlen=padlen, status=code, method=method, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                    if (present(status)) status = code
                    y = nanval
                    return
                end if
                    y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = filtfilt_status_success
    end function filtfilt_rank3_complex

    function sosfiltfilt_matrix_real(sos, x, padtype, padlen, status, axis, method, irlen) result(y)
        real(dp), intent(in) :: sos(:,:) !! SOS coefficients
        real(dp), intent(in) :: x(:,:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, none
        integer, intent(in), optional :: padlen !! number of samples padded
        integer, intent(out), optional :: status !! success or error
        integer, intent(in), optional :: axis !! one-based axis, default 2
        character(len=*), intent(in), optional :: method !! pad (default) or gust
        integer, intent(in), optional :: irlen !! positive truncation length for gust
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        real(dp), allocatable :: steady(:)
        integer :: edge, ntaps, steady_status
        character(:), allocatable :: mode

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2)))
        y = nanval
        if (present(status)) status = filtfilt_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along < 1 .or. along > 2) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            if (.not. valid_sos_gust_real(sos,size(x,along),irlen)) return
        else if (mode == 'pad') then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            ntaps = 2*size(sos,1) + 1 - min(count(sos(:,3) == 0.0_dp), &
                count(sos(:,6) == 0.0_dp))
            edge = pad_width(size(x,along), ntaps, padtype, padlen)
            if (edge < 0) return
            ! Validate steady-state initialization even if there are no channels.
            block
                real(dp), allocatable :: state(:,:)
                state = sosfilt_zi(sos, status=steady_status)
            end block
            if (steady_status /= sos_status_success) then
                if (present(status)) then
                    if (steady_status == sos_status_no_steady_state) status = filtfilt_status_no_steady_state
                end if
                return
            end if
        else
            return
        end if
        select case (along)
        case (1)
            do i = 1, size(x,2)
                line = sosfiltfilt_real(sos, x(:,i), padtype=padtype, &
                    padlen=padlen, status=code, method=mode, irlen=irlen)
                if (code /= filtfilt_status_success) then
                    y = nanval
                    if (present(status)) status = code
                    return
                end if
                y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                line = sosfiltfilt_real(sos, x(i,:), padtype=padtype, &
                    padlen=padlen, status=code, method=mode, irlen=irlen)
                if (code /= filtfilt_status_success) then
                    y = nanval
                    if (present(status)) status = code
                    return
                end if
                y(i,:) = line
            end do
        end select
        if (present(status)) status = filtfilt_status_success
    end function sosfiltfilt_matrix_real

    function sosfiltfilt_matrix_complex(sos, x, padtype, padlen, status, axis, method, irlen) result(y)
        complex(dp), intent(in) :: sos(:,:) !! SOS coefficients
        complex(dp), intent(in) :: x(:,:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, none
        integer, intent(in), optional :: padlen !! number of samples padded
        integer, intent(out), optional :: status !! success or error
        integer, intent(in), optional :: axis !! one-based axis, default 2
        character(len=*), intent(in), optional :: method !! pad (default) or gust
        integer, intent(in), optional :: irlen !! positive truncation length for gust
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        complex(dp), allocatable :: steady(:)
        integer :: edge, ntaps, steady_status
        character(:), allocatable :: mode

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2)))
        y = cmplx(nanval,nanval,dp)
        if (present(status)) status = filtfilt_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along < 1 .or. along > 2) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            if (.not. valid_sos_gust_complex(sos,size(x,along),irlen)) return
        else if (mode == 'pad') then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            ntaps = 2*size(sos,1) + 1 - min(count(sos(:,3) == (0.0_dp,0.0_dp)), &
                count(sos(:,6) == (0.0_dp,0.0_dp)))
            edge = pad_width(size(x,along), ntaps, padtype, padlen)
            if (edge < 0) return
            ! Validate steady-state initialization even if there are no channels.
            block
                complex(dp), allocatable :: state(:,:)
                state = sosfilt_zi(sos, status=steady_status)
            end block
            if (steady_status /= sos_status_success) then
                if (present(status)) then
                    if (steady_status == sos_status_no_steady_state) status = filtfilt_status_no_steady_state
                end if
                return
            end if
        else
            return
        end if
        select case (along)
        case (1)
            do i = 1, size(x,2)
                line = sosfiltfilt_complex(sos, x(:,i), padtype=padtype, &
                    padlen=padlen, status=code, method=mode, irlen=irlen)
                if (code /= filtfilt_status_success) then
                    y = cmplx(nanval,nanval,dp)
                    if (present(status)) status = code
                    return
                end if
                y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                line = sosfiltfilt_complex(sos, x(i,:), padtype=padtype, &
                    padlen=padlen, status=code, method=mode, irlen=irlen)
                if (code /= filtfilt_status_success) then
                    y = cmplx(nanval,nanval,dp)
                    if (present(status)) status = code
                    return
                end if
                y(i,:) = line
            end do
        end select
        if (present(status)) status = filtfilt_status_success
    end function sosfiltfilt_matrix_complex

    function sosfiltfilt_rank3_real(sos, x, padtype, padlen, status, axis, method, irlen) result(y)
        real(dp), intent(in) :: sos(:,:) !! SOS coefficients
        real(dp), intent(in) :: x(:,:,:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, none
        integer, intent(in), optional :: padlen !! number of samples padded
        integer, intent(out), optional :: status !! success or error
        integer, intent(in), optional :: axis !! one-based axis, default 3
        character(len=*), intent(in), optional :: method !! pad (default) or gust
        integer, intent(in), optional :: irlen !! positive truncation length for gust
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        real(dp), allocatable :: steady(:)
        integer :: edge, ntaps, steady_status
        character(:), allocatable :: mode

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2),size(x,3)))
        y = nanval
        if (present(status)) status = filtfilt_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 1 .or. along > 3) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            if (.not. valid_sos_gust_real(sos,size(x,along),irlen)) return
        else if (mode == 'pad') then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            ntaps = 2*size(sos,1) + 1 - min(count(sos(:,3) == 0.0_dp), &
                count(sos(:,6) == 0.0_dp))
            edge = pad_width(size(x,along), ntaps, padtype, padlen)
            if (edge < 0) return
            ! Validate steady-state initialization even if there are no channels.
            block
                real(dp), allocatable :: state(:,:)
                state = sosfilt_zi(sos, status=steady_status)
            end block
            if (steady_status /= sos_status_success) then
                if (present(status)) then
                    if (steady_status == sos_status_no_steady_state) status = filtfilt_status_no_steady_state
                end if
                return
            end if
        else
            return
        end if
        select case (along)
        case (1)
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    line = sosfiltfilt_real(sos, x(:,i,j), padtype=padtype, &
                        padlen=padlen, status=code, method=mode, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                        y = nanval
                        if (present(status)) status = code
                        return
                    end if
                    y(:,i,j) = line
                end do
            end do
        case (2)
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    line = sosfiltfilt_real(sos, x(i,:,j), padtype=padtype, &
                        padlen=padlen, status=code, method=mode, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                        y = nanval
                        if (present(status)) status = code
                        return
                    end if
                    y(i,:,j) = line
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    line = sosfiltfilt_real(sos, x(i,j,:), padtype=padtype, &
                        padlen=padlen, status=code, method=mode, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                        y = nanval
                        if (present(status)) status = code
                        return
                    end if
                    y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = filtfilt_status_success
    end function sosfiltfilt_rank3_real

    function sosfiltfilt_rank3_complex(sos, x, padtype, padlen, status, axis, method, irlen) result(y)
        complex(dp), intent(in) :: sos(:,:) !! SOS coefficients
        complex(dp), intent(in) :: x(:,:,:) !! signal to filter
        character(len=*), intent(in), optional :: padtype !! odd, even, constant, none
        integer, intent(in), optional :: padlen !! number of samples padded
        integer, intent(out), optional :: status !! success or error
        integer, intent(in), optional :: axis !! one-based axis, default 3
        character(len=*), intent(in), optional :: method !! pad (default) or gust
        integer, intent(in), optional :: irlen !! positive truncation length for gust
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        complex(dp), allocatable :: steady(:)
        integer :: edge, ntaps, steady_status
        character(:), allocatable :: mode

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2),size(x,3)))
        y = cmplx(nanval,nanval,dp)
        if (present(status)) status = filtfilt_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 1 .or. along > 3) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            if (.not. valid_sos_gust_complex(sos,size(x,along),irlen)) return
        else if (mode == 'pad') then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            ntaps = 2*size(sos,1) + 1 - min(count(sos(:,3) == (0.0_dp,0.0_dp)), &
                count(sos(:,6) == (0.0_dp,0.0_dp)))
            edge = pad_width(size(x,along), ntaps, padtype, padlen)
            if (edge < 0) return
            ! Validate steady-state initialization even if there are no channels.
            block
                complex(dp), allocatable :: state(:,:)
                state = sosfilt_zi(sos, status=steady_status)
            end block
            if (steady_status /= sos_status_success) then
                if (present(status)) then
                    if (steady_status == sos_status_no_steady_state) status = filtfilt_status_no_steady_state
                end if
                return
            end if
        else
            return
        end if
        select case (along)
        case (1)
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    line = sosfiltfilt_complex(sos, x(:,i,j), padtype=padtype, &
                        padlen=padlen, status=code, method=mode, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                        y = cmplx(nanval,nanval,dp)
                        if (present(status)) status = code
                        return
                    end if
                    y(:,i,j) = line
                end do
            end do
        case (2)
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    line = sosfiltfilt_complex(sos, x(i,:,j), padtype=padtype, &
                        padlen=padlen, status=code, method=mode, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                        y = cmplx(nanval,nanval,dp)
                        if (present(status)) status = code
                        return
                    end if
                    y(i,:,j) = line
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    line = sosfiltfilt_complex(sos, x(i,j,:), padtype=padtype, &
                        padlen=padlen, status=code, method=mode, irlen=irlen)
                    if (code /= filtfilt_status_success) then
                        y = cmplx(nanval,nanval,dp)
                        if (present(status)) status = code
                        return
                    end if
                    y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = filtfilt_status_success
    end function sosfiltfilt_rank3_complex
end module scifort_signal_filtfilt
