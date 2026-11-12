! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rank-independent forward/backward digital filtering with SciPy padding and Gustafsson states.
! The 1D implementations remain the single source of numerical semantics.
module scifort_signal_filtfilt_nd
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_filtfilt, only : filtfilt, sosfiltfilt, &
        valid_sos_gust_real, valid_sos_gust_complex, &
        filtfilt_status_success, filtfilt_status_invalid_input, &
        filtfilt_status_no_steady_state
    use scifort_signal_lfilter, only : lfilter_zi, lfilter_status_success, &
        lfilter_status_no_steady_state
    use scifort_signal_sos, only : sosfilt_zi, sos_status_success, &
        sos_status_no_steady_state
    implicit none
    private
    interface filtfilt_nd
        module procedure filtfilt_nd_real, filtfilt_nd_complex
        module procedure filtfilt_nd_complex_real_input
        module procedure filtfilt_nd_mixed_rrc, filtfilt_nd_mixed_rcr
        module procedure filtfilt_nd_mixed_rcc, filtfilt_nd_mixed_crr
        module procedure filtfilt_nd_mixed_crc
    end interface filtfilt_nd
    interface sosfiltfilt_nd
        module procedure sosfiltfilt_nd_real, sosfiltfilt_nd_complex
        module procedure sosfiltfilt_nd_real_sos_complex_input
        module procedure sosfiltfilt_nd_complex_sos_real_input
    end interface sosfiltfilt_nd
    public :: filtfilt_nd, sosfiltfilt_nd
contains
    pure logical function nd_valid_padding(n, taps, padtype, padlen) result(ok)
        integer, intent(in) :: n !! extent of the filtered axis
        integer, intent(in) :: taps !! default padding unit (3 times this)
        character(len=*), optional, intent(in) :: padtype !! odd, even, constant or none
        integer, optional, intent(in) :: padlen !! explicit number of samples on each side
        integer :: edge
        character(len=:), allocatable :: mode
        mode = 'odd'
        if (present(padtype)) mode = trim(padtype)
        ok = .false.
        select case (mode)
        case ('odd', 'even', 'constant', 'none')
        case default
            return
        end select
        if (taps < 1) return
        if (taps > huge(taps)/3) return
        edge = 3*taps
        if (present(padlen)) edge = padlen
        if (mode == 'none') edge = 0
        ok = edge >= 0 .and. n > edge
    end function nd_valid_padding

    ! Gustafsson does not need a steady-state solution or any signal padding.
    ! These checks also apply when the array has zero outer channels, and
    ! prevent such cases from silently accepting invalid filter parameters.
    pure logical function nd_valid_gust_real(b, a, n, irlen) result(ok)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        integer, intent(in) :: n !! extent of filtered axis
        integer, optional, intent(in) :: irlen !! positive truncation length
        ok = .false.
        if (n < 1 .or. size(b) < 1 .or. size(a) < 1) return
        if (abs(a(1)) <= tiny(1.0_dp)) return
        if (present(irlen)) then
            if (irlen < 1) return
        end if
        if (.not. all(ieee_is_finite(b))) return
        if (.not. all(ieee_is_finite(a))) return
        ok = .true.
    end function nd_valid_gust_real

    pure logical function nd_valid_gust_complex(b, a, n, irlen) result(ok)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        integer, intent(in) :: n !! extent of filtered axis
        integer, optional, intent(in) :: irlen !! positive truncation length
        ok = .false.
        if (n < 1 .or. size(b) < 1 .or. size(a) < 1) return
        if (abs(a(1)) <= tiny(1.0_dp)) return
        if (present(irlen)) then
            if (irlen < 1) return
        end if
        if (.not. all(ieee_is_finite(real(b,dp)))) return
        if (.not. all(ieee_is_finite(aimag(b)))) return
        if (.not. all(ieee_is_finite(real(a,dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        ok = .true.
    end function nd_valid_gust_complex

    ! With a contiguous assumed-rank actual, this reproduces the same
    ! column-major striding used for arbitrary-rank lfilter_nd/sosfilt_nd.
    ! A line occupies addresses outer*block + inner + 1 + sample*stride.
    subroutine filtfilt_nd_real(b, a, x, y, padtype, padlen, status, axis, method, irlen)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), contiguous, target, intent(in) :: x(..) !! rank 1 through 15 signal
        real(dp), contiguous, target, intent(out) :: y(..) !! same shape as x
        character(len=*), optional, intent(in) :: padtype !! odd, even, constant, none
        integer, optional, intent(in) :: padlen !! padding on each end
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: axis !! one-based axis, default last
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response length for gust
        real(dp), pointer :: xin(:) => null()
        real(dp), pointer :: yout(:) => null()
        real(dp), allocatable :: line(:)
        real(dp), allocatable :: filtered(:)
        real(dp), allocatable :: steady(:)
        real(dp) :: nanval
        character(len=:), allocatable :: mode
        integer :: dims(rank(x)), ax, taps, n, code, k
        integer(int64) :: stride, block, outer, i, j, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = filtfilt_status_invalid_input
        if (rank(y) == 0) then
            select rank (y)
            rank (0)
                y = nanval
            end select
            return
        end if
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = nanval
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode /= 'pad' .and. mode /= 'gust') return
        if (size(b) < 1 .or. size(a) < 1) return
        if (mode == 'gust') then
            if (.not. nd_valid_gust_real(b,a,size(x,ax),irlen)) return
        else
            taps = max(size(a),size(b))
            if (.not. nd_valid_padding(size(x,ax),taps,padtype,padlen)) return
            steady = lfilter_zi(b,a,status=code)
            if (code /= lfilter_status_success) then
                if (present(status)) then
                    if (code == lfilter_status_no_steady_state) &
                        status = filtfilt_status_no_steady_state
                end if
                return
            end if
        end if
        dims = shape(x)
        n = dims(ax)
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(n,int64)
        outer = product(int(dims(ax+1:),int64))
        if (size(x,kind=int64) > 0_int64) &
            call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(line(n))
        do j = 0_int64,outer-1_int64
            do i = 0_int64,stride-1_int64
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    line(k+1) = xin(pos)
                end do
                filtered = filtfilt(b, a,line, &
                    padtype=padtype,padlen=padlen,status=code, &
                    method=method,irlen=irlen)
                if (code /= filtfilt_status_success) then
                    yout = nanval
                    if (present(status)) status = code
                    return
                end if
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    yout(pos) = filtered(k+1)
                end do
            end do
        end do
        if (present(status)) status = filtfilt_status_success
    end subroutine filtfilt_nd_real

    subroutine filtfilt_nd_complex(b, a, x, y, padtype, padlen, status, axis, method, irlen)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), contiguous, target, intent(in) :: x(..) !! rank 1 through 15 signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! same shape as x
        character(len=*), optional, intent(in) :: padtype !! odd, even, constant, none
        integer, optional, intent(in) :: padlen !! padding on each end
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: axis !! one-based axis, default last
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response length for gust
        complex(dp), pointer :: xin(:) => null()
        complex(dp), pointer :: yout(:) => null()
        complex(dp), allocatable :: line(:)
        complex(dp), allocatable :: filtered(:)
        complex(dp), allocatable :: steady(:)
        real(dp) :: nanval
        character(len=:), allocatable :: mode
        integer :: dims(rank(x)), ax, taps, n, code, k
        integer(int64) :: stride, block, outer, i, j, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = filtfilt_status_invalid_input
        if (rank(y) == 0) then
            select rank (y)
            rank (0)
                y = cmplx(nanval,nanval,dp)
            end select
            return
        end if
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = cmplx(nanval,nanval,dp)
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode /= 'pad' .and. mode /= 'gust') return
        if (size(b) < 1 .or. size(a) < 1) return
        if (mode == 'gust') then
            if (.not. nd_valid_gust_complex(b,a,size(x,ax),irlen)) return
        else
            taps = max(size(a),size(b))
            if (.not. nd_valid_padding(size(x,ax),taps,padtype,padlen)) return
            steady = lfilter_zi(b,a,status=code)
            if (code /= lfilter_status_success) then
                if (present(status)) then
                    if (code == lfilter_status_no_steady_state) &
                        status = filtfilt_status_no_steady_state
                end if
                return
            end if
        end if
        dims = shape(x)
        n = dims(ax)
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(n,int64)
        outer = product(int(dims(ax+1:),int64))
        if (size(x,kind=int64) > 0_int64) &
            call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(line(n))
        do j = 0_int64,outer-1_int64
            do i = 0_int64,stride-1_int64
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    line(k+1) = xin(pos)
                end do
                filtered = filtfilt(b, a,line, &
                    padtype=padtype,padlen=padlen,status=code, &
                    method=method,irlen=irlen)
                if (code /= filtfilt_status_success) then
                    yout = cmplx(nanval,nanval,dp)
                    if (present(status)) status = code
                    return
                end if
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    yout(pos) = filtered(k+1)
                end do
            end do
        end do
        if (present(status)) status = filtfilt_status_success
    end subroutine filtfilt_nd_complex

    subroutine filtfilt_nd_complex_real_input(b, a, x, y, padtype, padlen, status, axis, method, irlen)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), contiguous, target, intent(in) :: x(..) !! rank 1 through 15 signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! same shape as x
        character(len=*), optional, intent(in) :: padtype !! odd, even, constant, none
        integer, optional, intent(in) :: padlen !! padding on each end
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: axis !! one-based axis, default last
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response length for gust
        real(dp), pointer :: xin(:) => null()
        complex(dp), pointer :: yout(:) => null()
        real(dp), allocatable :: line(:)
        complex(dp), allocatable :: filtered(:)
        complex(dp), allocatable :: steady(:)
        real(dp) :: nanval
        character(len=:), allocatable :: mode
        integer :: dims(rank(x)), ax, taps, n, code, k
        integer(int64) :: stride, block, outer, i, j, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = filtfilt_status_invalid_input
        if (rank(y) == 0) then
            select rank (y)
            rank (0)
                y = cmplx(nanval,nanval,dp)
            end select
            return
        end if
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = cmplx(nanval,nanval,dp)
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode /= 'pad' .and. mode /= 'gust') return
        if (size(b) < 1 .or. size(a) < 1) return
        if (mode == 'gust') then
            if (.not. nd_valid_gust_complex(b,a,size(x,ax),irlen)) return
        else
            taps = max(size(a),size(b))
            if (.not. nd_valid_padding(size(x,ax),taps,padtype,padlen)) return
            steady = lfilter_zi(b,a,status=code)
            if (code /= lfilter_status_success) then
                if (present(status)) then
                    if (code == lfilter_status_no_steady_state) &
                        status = filtfilt_status_no_steady_state
                end if
                return
            end if
        end if
        dims = shape(x)
        n = dims(ax)
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(n,int64)
        outer = product(int(dims(ax+1:),int64))
        if (size(x,kind=int64) > 0_int64) &
            call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(line(n))
        do j = 0_int64,outer-1_int64
            do i = 0_int64,stride-1_int64
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    line(k+1) = xin(pos)
                end do
                filtered = filtfilt(b, a,cmplx(line,0.0_dp,dp), &
                    padtype=padtype,padlen=padlen,status=code, &
                    method=method,irlen=irlen)
                if (code /= filtfilt_status_success) then
                    yout = cmplx(nanval,nanval,dp)
                    if (present(status)) status = code
                    return
                end if
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    yout(pos) = filtered(k+1)
                end do
            end do
        end do
        if (present(status)) status = filtfilt_status_success
    end subroutine filtfilt_nd_complex_real_input

    subroutine sosfiltfilt_nd_real(sos, x, y, padtype, padlen, status, axis, method, irlen)
        real(dp), intent(in) :: sos(:,:) !! sections, with a0=1
        real(dp), contiguous, target, intent(in) :: x(..) !! rank 1 through 15 signal
        real(dp), contiguous, target, intent(out) :: y(..) !! same shape as x
        character(len=*), optional, intent(in) :: padtype !! odd, even, constant, none
        integer, optional, intent(in) :: padlen !! padding on each end
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: axis !! one-based axis, default last
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response cutoff
        real(dp), pointer :: xin(:) => null()
        real(dp), pointer :: yout(:) => null()
        real(dp), allocatable :: line(:)
        real(dp), allocatable :: filtered(:)
        real(dp), allocatable :: steady(:,:)
        real(dp) :: nanval
        integer :: dims(rank(x)), ax, taps, n, code, k
        character(:), allocatable :: mode
        integer(int64) :: stride, block, outer, i, j, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = filtfilt_status_invalid_input
        if (rank(y) == 0) then
            select rank (y)
            rank (0)
                y = nanval
            end select
            return
        end if
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = nanval
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            if (.not. valid_sos_gust_real(sos,size(x,ax),irlen)) return
        else if (mode == 'pad') then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            taps = 2*size(sos,1)+1-min(count(sos(:,3) == 0.0_dp), &
                count(sos(:,6) == 0.0_dp))
            if (.not. nd_valid_padding(size(x,ax),taps,padtype,padlen)) return
            steady = sosfilt_zi(sos,status=code)
            if (code /= sos_status_success) then
                if (present(status)) then
                    if (code == sos_status_no_steady_state) &
                        status = filtfilt_status_no_steady_state
                end if
                return
            end if
        else
            return
        end if
        dims = shape(x)
        n = dims(ax)
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(n,int64)
        outer = product(int(dims(ax+1:),int64))
        if (size(x,kind=int64) > 0_int64) &
            call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(line(n))
        do j = 0_int64,outer-1_int64
            do i = 0_int64,stride-1_int64
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    line(k+1) = xin(pos)
                end do
                filtered = sosfiltfilt(sos,line, &
                    padtype=padtype,padlen=padlen,status=code, &
                    method=mode,irlen=irlen)
                if (code /= filtfilt_status_success) then
                    yout = nanval
                    if (present(status)) status = code
                    return
                end if
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    yout(pos) = filtered(k+1)
                end do
            end do
        end do
        if (present(status)) status = filtfilt_status_success
    end subroutine sosfiltfilt_nd_real

    subroutine sosfiltfilt_nd_complex(sos, x, y, padtype, padlen, status, axis, method, irlen)
        complex(dp), intent(in) :: sos(:,:) !! sections, with a0=1
        complex(dp), contiguous, target, intent(in) :: x(..) !! rank 1 through 15 signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! same shape as x
        character(len=*), optional, intent(in) :: padtype !! odd, even, constant, none
        integer, optional, intent(in) :: padlen !! padding on each end
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: axis !! one-based axis, default last
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response cutoff
        complex(dp), pointer :: xin(:) => null()
        complex(dp), pointer :: yout(:) => null()
        complex(dp), allocatable :: line(:)
        complex(dp), allocatable :: filtered(:)
        complex(dp), allocatable :: steady(:,:)
        real(dp) :: nanval
        integer :: dims(rank(x)), ax, taps, n, code, k
        character(:), allocatable :: mode
        integer(int64) :: stride, block, outer, i, j, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = filtfilt_status_invalid_input
        if (rank(y) == 0) then
            select rank (y)
            rank (0)
                y = cmplx(nanval,nanval,dp)
            end select
            return
        end if
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = cmplx(nanval,nanval,dp)
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            if (.not. valid_sos_gust_complex(sos,size(x,ax),irlen)) return
        else if (mode == 'pad') then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            taps = 2*size(sos,1)+1-min(count(sos(:,3) == (0.0_dp,0.0_dp)), &
                count(sos(:,6) == (0.0_dp,0.0_dp)))
            if (.not. nd_valid_padding(size(x,ax),taps,padtype,padlen)) return
            steady = sosfilt_zi(sos,status=code)
            if (code /= sos_status_success) then
                if (present(status)) then
                    if (code == sos_status_no_steady_state) &
                        status = filtfilt_status_no_steady_state
                end if
                return
            end if
        else
            return
        end if
        dims = shape(x)
        n = dims(ax)
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(n,int64)
        outer = product(int(dims(ax+1:),int64))
        if (size(x,kind=int64) > 0_int64) &
            call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(line(n))
        do j = 0_int64,outer-1_int64
            do i = 0_int64,stride-1_int64
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    line(k+1) = xin(pos)
                end do
                filtered = sosfiltfilt(sos,line, &
                    padtype=padtype,padlen=padlen,status=code, &
                    method=mode,irlen=irlen)
                if (code /= filtfilt_status_success) then
                    yout = cmplx(nanval,nanval,dp)
                    if (present(status)) status = code
                    return
                end if
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    yout(pos) = filtered(k+1)
                end do
            end do
        end do
        if (present(status)) status = filtfilt_status_success
    end subroutine sosfiltfilt_nd_complex

    subroutine filtfilt_nd_mixed_rrc(b, a, x, y, padtype, padlen, status, axis, method, irlen)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), contiguous, target, intent(in) :: x(..) !! signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output
        character(len=*), optional, intent(in) :: padtype !! padding mode
        integer, optional, intent(in) :: padlen !! padding length
        integer, optional, intent(out) :: status !! result status
        integer, optional, intent(in) :: axis !! one-based axis
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response length for gust
        call filtfilt_nd_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),x,y,padtype,padlen,status,axis,method,irlen)
    end subroutine filtfilt_nd_mixed_rrc

    subroutine filtfilt_nd_mixed_rcr(b, a, x, y, padtype, padlen, status, axis, method, irlen)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), contiguous, target, intent(in) :: x(..) !! signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output
        character(len=*), optional, intent(in) :: padtype !! padding mode
        integer, optional, intent(in) :: padlen !! padding length
        integer, optional, intent(out) :: status !! result status
        integer, optional, intent(in) :: axis !! one-based axis
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response length for gust
        call filtfilt_nd_complex_real_input(cmplx(b,0.0_dp,dp),a,x,y,padtype,padlen,status,axis,method,irlen)
    end subroutine filtfilt_nd_mixed_rcr

    subroutine filtfilt_nd_mixed_rcc(b, a, x, y, padtype, padlen, status, axis, method, irlen)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), contiguous, target, intent(in) :: x(..) !! signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output
        character(len=*), optional, intent(in) :: padtype !! padding mode
        integer, optional, intent(in) :: padlen !! padding length
        integer, optional, intent(out) :: status !! result status
        integer, optional, intent(in) :: axis !! one-based axis
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response length for gust
        call filtfilt_nd_complex(cmplx(b,0.0_dp,dp),a,x,y,padtype,padlen,status,axis,method,irlen)
    end subroutine filtfilt_nd_mixed_rcc

    subroutine filtfilt_nd_mixed_crr(b, a, x, y, padtype, padlen, status, axis, method, irlen)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), contiguous, target, intent(in) :: x(..) !! signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output
        character(len=*), optional, intent(in) :: padtype !! padding mode
        integer, optional, intent(in) :: padlen !! padding length
        integer, optional, intent(out) :: status !! result status
        integer, optional, intent(in) :: axis !! one-based axis
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response length for gust
        call filtfilt_nd_complex_real_input(b,cmplx(a,0.0_dp,dp),x,y,padtype,padlen,status,axis,method,irlen)
    end subroutine filtfilt_nd_mixed_crr

    subroutine filtfilt_nd_mixed_crc(b, a, x, y, padtype, padlen, status, axis, method, irlen)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), contiguous, target, intent(in) :: x(..) !! signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output
        character(len=*), optional, intent(in) :: padtype !! padding mode
        integer, optional, intent(in) :: padlen !! padding length
        integer, optional, intent(out) :: status !! result status
        integer, optional, intent(in) :: axis !! one-based axis
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response length for gust
        call filtfilt_nd_complex(b,cmplx(a,0.0_dp,dp),x,y,padtype,padlen,status,axis,method,irlen)
    end subroutine filtfilt_nd_mixed_crc

    subroutine sosfiltfilt_nd_real_sos_complex_input(sos,x,y,padtype,padlen,status,axis,method,irlen)
        real(dp), intent(in) :: sos(:,:) !! real SOS coefficients
        complex(dp), contiguous, target, intent(in) :: x(..) !! complex samples
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output
        character(len=*), optional, intent(in) :: padtype !! padding mode
        integer, optional, intent(in) :: padlen !! padding length
        integer, optional, intent(out) :: status !! result status
        integer, optional, intent(in) :: axis !! one-based axis
        character(len=*), optional, intent(in) :: method !! pad or gust
        integer, optional, intent(in) :: irlen !! positive impulse response cutoff
        call sosfiltfilt_nd_complex(cmplx(sos,0.0_dp,dp),x,y, &
            padtype,padlen,status,axis,method,irlen)
    end subroutine sosfiltfilt_nd_real_sos_complex_input

    subroutine sosfiltfilt_nd_complex_sos_real_input(sos, x, y, padtype, padlen, status, axis, method, irlen)
        complex(dp), intent(in) :: sos(:,:) !! sections, with a0=1
        real(dp), contiguous, target, intent(in) :: x(..) !! rank 1 through 15 signal
        complex(dp), contiguous, target, intent(out) :: y(..) !! same shape as x
        character(len=*), optional, intent(in) :: padtype !! odd, even, constant, none
        integer, optional, intent(in) :: padlen !! padding on each end
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: axis !! one-based axis, default last
        character(len=*), optional, intent(in) :: method !! 'pad' (default) or 'gust'
        integer, optional, intent(in) :: irlen !! positive impulse response cutoff
        real(dp), pointer :: xin(:) => null()
        complex(dp), pointer :: yout(:) => null()
        real(dp), allocatable :: line(:)
        complex(dp), allocatable :: filtered(:)
        complex(dp), allocatable :: steady(:,:)
        real(dp) :: nanval
        integer :: dims(rank(x)), ax, taps, n, code, k
        character(:), allocatable :: mode
        integer(int64) :: stride, block, outer, i, j, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = filtfilt_status_invalid_input
        if (rank(y) == 0) then
            select rank (y)
            rank (0)
                y = cmplx(nanval,nanval,dp)
            end select
            return
        end if
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = cmplx(nanval,nanval,dp)
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        mode = 'pad'
        if (present(method)) mode = trim(method)
        if (mode == 'gust') then
            if (.not. valid_sos_gust_complex(sos,size(x,ax),irlen)) return
        else if (mode == 'pad') then
            if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
            taps = 2*size(sos,1)+1-min(count(sos(:,3) == (0.0_dp,0.0_dp)), &
                count(sos(:,6) == (0.0_dp,0.0_dp)))
            if (.not. nd_valid_padding(size(x,ax),taps,padtype,padlen)) return
            steady = sosfilt_zi(sos,status=code)
            if (code /= sos_status_success) then
                if (present(status)) then
                    if (code == sos_status_no_steady_state) &
                        status = filtfilt_status_no_steady_state
                end if
                return
            end if
        else
            return
        end if
        dims = shape(x)
        n = dims(ax)
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(n,int64)
        outer = product(int(dims(ax+1:),int64))
        if (size(x,kind=int64) > 0_int64) &
            call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(line(n))
        do j = 0_int64,outer-1_int64
            do i = 0_int64,stride-1_int64
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    line(k+1) = xin(pos)
                end do
                filtered = sosfiltfilt(sos,cmplx(line,0.0_dp,dp), &
                    padtype=padtype,padlen=padlen,status=code, &
                    method=mode,irlen=irlen)
                if (code /= filtfilt_status_success) then
                    yout = cmplx(nanval,nanval,dp)
                    if (present(status)) status = code
                    return
                end if
                do k = 0,n-1
                    pos = j*block+i+1_int64+int(k,int64)*stride
                    yout(pos) = filtered(k+1)
                end do
            end do
        end do
        if (present(status)) status = filtfilt_status_success
    end subroutine sosfiltfilt_nd_complex_sos_real_input

end module scifort_signal_filtfilt_nd
