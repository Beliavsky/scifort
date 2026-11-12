! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Fast linear convolution using reusable FFT plans. The caller-output fftconvolve_nd
! handles ranks 1-15 without assuming a statically known output rank.
module scifort_signal_fftconvolve
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft, fft_plan
    implicit none
    private
    integer, parameter :: convolution_status_success = 0
    integer, parameter :: convolution_status_invalid_input = 1

    interface fftconvolve
        module procedure fftconvolve_real, fftconvolve_complex
        module procedure fftconvolve_real_complex, fftconvolve_complex_real
        module procedure fftconvolve_r2_real, fftconvolve_r2_complex
        module procedure fftconvolve_r2_real_complex, fftconvolve_r2_complex_real
        module procedure fftconvolve_r3_real, fftconvolve_r3_complex
        module procedure fftconvolve_r3_real_complex, fftconvolve_r3_complex_real
    end interface
    interface fftconvolve_nd
        module procedure fftconvolve_nd_rr, fftconvolve_nd_cc
        module procedure fftconvolve_nd_rc, fftconvolve_nd_cr
    end interface fftconvolve_nd
    interface oaconvolve
        module procedure oaconvolve_real, oaconvolve_complex
        module procedure oaconvolve_real_complex, oaconvolve_complex_real
    end interface
    interface oaconvolve_nd
        module procedure oaconvolve_nd_rr, oaconvolve_nd_cc
        module procedure oaconvolve_nd_rc, oaconvolve_nd_cr
    end interface oaconvolve_nd
    public :: fftconvolve, fftconvolve_nd, oaconvolve, oaconvolve_nd
    ! Internal packed-array kernel shared with rank-independent direct/correlation
    ! dispatch. Do not expose this low-level routine through scifort_signal.
    public :: fftconvolve_flat

contains

    pure subroutine conv_slice(n,m,mode,offset,count,ok)
        integer, intent(in) :: n !! first signal length
        integer, intent(in) :: m !! second signal length
        character(*), intent(in) :: mode !! full, same, valid
        integer, intent(out) :: offset !! one-based slice offset
        integer, intent(out) :: count !! returned elements
        logical, intent(out) :: ok !! valid input and mode

        offset = 1
        count = 0
        ok = n>0 .and. m>0
        if (.not. ok) return
        select case (mode)
        case ('full')
            count = n+m-1
        case ('same')
            offset = 1+(m-1)/2
            count = n
        case ('valid')
            offset = min(n,m)
            count = abs(n-m)+1
        case default
            ok = .false.
        end select
    end subroutine conv_slice

    subroutine backward_fft(x)
        complex(dp), intent(inout) :: x(:) !! transform samples in place
        x = conjg(x)
        call forward_fft(x)
        x = conjg(x)/real(size(x),dp)
    end subroutine backward_fft

    pure function power_two_at_least(n) result(nfft)
        integer, intent(in) :: n !! minimum transform length
        integer :: nfft
        nfft = 1
        do while (nfft < n)
            if (nfft > huge(nfft)/2) then
                nfft = 0
                return
            end if
            nfft = 2*nfft
        end do
    end function power_two_at_least

    ! Find the smallest supported mixed-radix transform length >= n.
    ! Enumerating 3, 5, 7 factors and rounding each candidate up by powers
    ! of two avoids scanning large prime gaps. Small direct-DFT lengths are
    ! intentionally left unchanged. Return zero if no supported length fits
    ! in a default integer (the caller already treats zero as invalid).
    pure function fast_fft_length(n) result(nfft)
        integer, intent(in) :: n !! smallest acceptable transform length
        integer :: nfft, quotient, factor
        integer(int64) :: best, p3, p5, p7, base, candidate

        nfft = 0
        if (n < 1) return
        if (n <= 16) then
            nfft = n
            return
        end if
        quotient = n
        do factor = 2,7
            if (factor /= 2 .and. factor /= 3 .and. factor /= 5 .and. factor /= 7) cycle
            do while (mod(quotient,factor) == 0)
                quotient = quotient/factor
            end do
        end do
        if (quotient == 1) then
            nfft = n
            return
        end if

        nfft = power_two_at_least(n)
        best = int(huge(nfft),int64)+1_int64
        if (nfft > 0) best = int(nfft,int64)
        nfft = 0
        p7 = 1_int64
        do
            p5 = 1_int64
            do
                p3 = 1_int64
                do
                    base = p7*p5*p3
                    candidate = base
                    do while (candidate < int(n,int64))
                        candidate = 2_int64*candidate
                    end do
                    best = min(best,candidate)
                    if (best == int(n,int64)) then
                        nfft = n
                        return
                    end if
                    if (base > best/3_int64) exit
                    p3 = 3_int64*p3
                end do
                if (p7*p5 > best/5_int64) exit
                p5 = 5_int64*p5
            end do
            if (p7 > best/7_int64) exit
            p7 = 7_int64*p7
        end do
        if (best <= int(huge(nfft),int64)) nfft = int(best)
    end function fast_fft_length

    ! Relative cost per complex sample for this implementation's planned FFT.
    ! Radix-2 is considerably faster per point than the recursive mixed-radix
    ! kernel. The fixed term approximates workspace initialization, packing,
    ! and cropping per element. This deliberately favors predictable speed
    ! over choosing the smallest possible padded tensor.
    pure function convolution_padding_cost(padded,selected) result(cost)
        integer, intent(in) :: padded(:) !! candidate padded dimensions
        logical, intent(in) :: selected(:) !! axes requiring FFTs
        real(dp) :: cost, per_point
        integer :: j

        cost = 0.0_dp
        per_point = 2.5_dp
        do j=1,size(padded)
            if (.not. selected(j) .or. padded(j)==1) cycle
            if (iand(padded(j),padded(j)-1)==0) then
                per_point = per_point+1.0_dp
            else
                per_point = per_point+3.8_dp
            end if
        end do
        cost = product(real(padded,dp))*per_point
    end function convolution_padding_cost

    pure function oa_chunk_length(long_len,short_len) result(chunk)
        integer, intent(in) :: long_len !! length of candidate split operand
        integer, intent(in) :: short_len !! opposing operand's length
        integer :: chunk, nf

        chunk = 0
        if (short_len < 1 .or. long_len < 1) return
        if (int(long_len,int64) < 4_int64*int(short_len,int64)) return
        if (short_len > huge(nf)/4) return
        nf = power_two_at_least(max(64,4*short_len))
        if (nf == 0) return
        if (long_len > nf-short_len+1) chunk = nf-short_len+1
    end function oa_chunk_length

    function complex_full_fft(a,b,overlap_add) result(y)
        complex(dp), intent(in) :: a(:) !! first signal
        complex(dp), intent(in) :: b(:) !! second signal
        logical, intent(in) :: overlap_add !! true uses blocked overlap-add
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: fx(:), fy(:)
        type(fft_plan) :: plan
        integer :: fft_status
        integer :: n,m,nf,first,last,chunk,k
        complex(dp), allocatable :: long(:),short(:)

        n = size(a)
        m = size(b)
        allocate(y(n+m-1))
        y = (0.0_dp,0.0_dp)
        if (.not. overlap_add .or. max(n,m) <= 4*min(n,m)) then
            nf = power_two_at_least(n+m-1)
            if (nf == 0) return
            call plan%initialize(nf,fft_status)
            if (fft_status /= 0) return
            allocate(fx(nf),fy(nf))
            fx = (0.0_dp,0.0_dp)
            fy = (0.0_dp,0.0_dp)
            fx(1:n) = a
            fy(1:m) = b
            call plan%execute(fx,fft_status)
            call plan%execute(fy,fft_status)
            fx = fx*fy
            call plan%inverse(fx,fft_status)
            y = fx(1:n+m-1)
            return
        end if
        if (n >= m) then
            long = a
            short = b
        else
            long = b
            short = a
        end if
        ! Aim for four kernel lengths per FFT: a bounded, reusable kernel transform.
        if (size(short) > huge(nf)/4) return
        nf = power_two_at_least(max(64,4*size(short)))
        if (nf == 0) return
        chunk = nf-size(short)+1
        call plan%initialize(nf,fft_status)
        if (fft_status /= 0) return
        allocate(fx(nf),fy(nf))
        fy = (0.0_dp,0.0_dp)
        fy(1:size(short)) = short
        call plan%execute(fy,fft_status)
        do first = 1,size(long),chunk
            last = min(size(long),first+chunk-1)
            fx = (0.0_dp,0.0_dp)
            fx(1:last-first+1) = long(first:last)
            call plan%execute(fx,fft_status)
            fx = fx*fy
            call plan%inverse(fx,fft_status)
            do k = 1,last-first+size(short)
                y(first+k-1) = y(first+k-1)+fx(k)
            end do
        end do
    end function complex_full_fft

    function fftconvolve_real(in1,in2,mode,status,axes) result(y)
        real(dp), intent(in) :: in1(:) !! first real array
        real(dp), intent(in) :: in2(:) !! second real array
        character(*), intent(in), optional :: mode !! full (default), same, valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based selected axis; only [1] is valid
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: full(:)
        character(:), allocatable :: selected
        integer :: offset,count
        logical :: ok

        if (present(status)) status = convolution_status_invalid_input
        allocate(y(0))
        if (present(axes)) then
            if (size(axes) /= 1) return
            if (axes(1) /= 1) return
        end if
        selected = 'full'
        if (present(mode)) selected = trim(mode)
        call conv_slice(size(in1),size(in2),selected,offset,count,ok)
        if (.not. ok) return
        if (.not. all(ieee_is_finite(in1)) .or. .not. all(ieee_is_finite(in2))) return
        full = complex_full_fft(cmplx(in1,0.0_dp,dp),cmplx(in2,0.0_dp,dp),.false.)
        y = real(full(offset:offset+count-1),dp)
        if (present(status)) status = convolution_status_success
    end function fftconvolve_real

    function fftconvolve_complex(in1,in2,mode,status,axes) result(y)
        complex(dp), intent(in) :: in1(:) !! first complex signal
        complex(dp), intent(in) :: in2(:) !! second complex signal
        character(*), intent(in), optional :: mode !! full, same, valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based selected axis; only [1] is valid
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: full(:)
        character(:), allocatable :: selected
        integer :: offset,count
        logical :: ok

        if (present(status)) status = convolution_status_invalid_input
        allocate(y(0))
        if (present(axes)) then
            if (size(axes) /= 1) return
            if (axes(1) /= 1) return
        end if
        selected = 'full'
        if (present(mode)) selected = trim(mode)
        call conv_slice(size(in1),size(in2),selected,offset,count,ok)
        if (.not. ok) return
        if (.not. all(ieee_is_finite(real(in1,dp))) .or. &
            .not. all(ieee_is_finite(aimag(in1))) .or. &
            .not. all(ieee_is_finite(real(in2,dp))) .or. &
            .not. all(ieee_is_finite(aimag(in2)))) return
        full = complex_full_fft(in1,in2,.false.)
        y = full(offset:offset+count-1)
        if (present(status)) status = convolution_status_success
    end function fftconvolve_complex

    function fftconvolve_real_complex(in1, in2, mode, status, axes) result(y)
        real(dp), intent(in) :: in1(:) !! real signal
        complex(dp), intent(in) :: in2(:) !! complex signal
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! selected one-based axis; only [1] is valid
        complex(dp), allocatable :: y(:)

        y = fftconvolve_complex(cmplx(in1, 0.0_dp, dp), in2, mode, status, axes)
    end function fftconvolve_real_complex

    function fftconvolve_complex_real(in1, in2, mode, status, axes) result(y)
        complex(dp), intent(in) :: in1(:) !! complex signal
        real(dp), intent(in) :: in2(:) !! real signal
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! selected one-based axis; only [1] is valid
        complex(dp), allocatable :: y(:)

        y = fftconvolve_complex(in1, cmplx(in2, 0.0_dp, dp), mode, status, axes)
    end function fftconvolve_complex_real

    function oaconvolve_real(in1,in2,mode,status) result(y)
        real(dp), intent(in) :: in1(:) !! first real signal
        real(dp), intent(in) :: in2(:) !! second real signal
        character(*), intent(in), optional :: mode !! full (default), same, valid
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: full(:)
        character(:), allocatable :: selected
        integer :: offset,count
        logical :: ok

        if (present(status)) status = convolution_status_invalid_input
        allocate(y(0))
        selected = 'full'
        if (present(mode)) selected = trim(mode)
        call conv_slice(size(in1),size(in2),selected,offset,count,ok)
        if (.not. ok) return
        if (.not. all(ieee_is_finite(in1)) .or. .not. all(ieee_is_finite(in2))) return
        full = complex_full_fft(cmplx(in1,0.0_dp,dp),cmplx(in2,0.0_dp,dp),.true.)
        y = real(full(offset:offset+count-1),dp)
        if (present(status)) status = convolution_status_success
    end function oaconvolve_real

    function oaconvolve_complex(in1,in2,mode,status) result(y)
        complex(dp), intent(in) :: in1(:) !! first complex signal
        complex(dp), intent(in) :: in2(:) !! second complex signal
        character(*), intent(in), optional :: mode !! full, same, valid
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: full(:)
        character(:), allocatable :: selected
        integer :: offset,count
        logical :: ok

        if (present(status)) status = convolution_status_invalid_input
        allocate(y(0))
        selected = 'full'
        if (present(mode)) selected = trim(mode)
        call conv_slice(size(in1),size(in2),selected,offset,count,ok)
        if (.not. ok) return
        if (.not. all(ieee_is_finite(real(in1,dp))) .or. &
            .not. all(ieee_is_finite(aimag(in1))) .or. &
            .not. all(ieee_is_finite(real(in2,dp))) .or. &
            .not. all(ieee_is_finite(aimag(in2)))) return
        full = complex_full_fft(in1,in2,.true.)
        y = full(offset:offset+count-1)
        if (present(status)) status = convolution_status_success
    end function oaconvolve_complex


    function oaconvolve_real_complex(in1, in2, mode, status) result(y)
        real(dp), intent(in) :: in1(:) !! real signal
        complex(dp), intent(in) :: in2(:) !! complex signal
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)

        y = oaconvolve_complex(cmplx(in1, 0.0_dp, dp), in2, mode, status)
    end function oaconvolve_real_complex

    function oaconvolve_complex_real(in1, in2, mode, status) result(y)
        complex(dp), intent(in) :: in1(:) !! complex signal
        real(dp), intent(in) :: in2(:) !! real signal
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)

        y = oaconvolve_complex(in1, cmplx(in2, 0.0_dp, dp), mode, status)
    end function oaconvolve_complex_real

    ! Multidimensional FFTs are separable. Non-convolved axes use singleton
    ! broadcasting; the optional axes use one-based Fortran indices.
    subroutine fft_axis_set(rank, axes, selected, ok)
        integer, intent(in) :: rank !! rank of both input arrays
        integer, intent(in), optional :: axes(:) !! selected one-based axes, default all
        logical, intent(out) :: selected(rank) !! true on a convolution axis
        logical, intent(out) :: ok !! selection validity
        integer :: k, idx

        selected = .true.
        ok = .true.
        if (.not. present(axes)) return
        selected = .false.
        ok = size(axes) > 0
        do k = 1, size(axes)
            idx = axes(k)
            if (idx < 1 .or. idx > rank) then
                ok = .false.
                return
            end if
            if (selected(idx)) then
                ok = .false.
                return
            end if
            selected(idx) = .true.
        end do
    end subroutine fft_axis_set

    subroutine fft_nd_layout(n, m, selected, mode, padded, first, count, ok)
        integer, intent(in) :: n(:) !! shape of the first input
        integer, intent(in) :: m(:) !! shape of the second input
        logical, intent(in) :: selected(:) !! axes on which convolution is performed
        character(*), intent(in) :: mode !! full, same or valid
        integer, intent(out) :: padded(:) !! FFT working-array shape
        integer, intent(out) :: first(:) !! one-based crop origin
        integer, intent(out) :: count(:) !! output shape
        logical, intent(out) :: ok !! valid shape and mode
        integer :: j, extent, compact(size(n)), trial(size(n))
        real(dp) :: best_cost, cost

        padded = 0
        count = 0
        first = 1
        ok = all(n > 0) .and. all(m > 0)
        if (.not. ok) return
        select case (mode)
        case ('full','same','valid')
        case default
            ok = .false.
            return
        end select
        if (mode == 'valid') then
            ok = all(pack(n >= m, selected)) .or. all(pack(n <= m, selected))
            if (.not. ok) return
        end if
        do j = 1, size(n)
            if (selected(j)) then
                if (n(j) > huge(extent) - m(j)) then
                    ok = .false.
                    return
                end if
                extent = n(j) + m(j) - 1
                padded(j) = power_two_at_least(extent)
                select case (mode)
                case ('full')
                    count(j) = extent
                case ('same')
                    first(j) = 1 + (m(j) - 1)/2
                    count(j) = n(j)
                case ('valid')
                    first(j) = min(n(j),m(j))
                    count(j) = abs(n(j)-m(j)) + 1
                end select
            else
                if (n(j) /= m(j) .and. n(j) /= 1 .and. m(j) /= 1) then
                    ok = .false.
                    return
                end if
                padded(j) = max(n(j),m(j))
                count(j) = padded(j)
                if (mode == 'same') then
                    count(j) = n(j)
                    first(j) = 1 + (padded(j) - n(j))/2
                end if
            end if
            if (padded(j) == 0) then
                ok = .false.
                return
            end if
        end do

        ! Evaluate mixed-radix padding only when its smaller tensor outweighs
        ! the slower per-point transform. Use a conservative 12% cost margin;
        ! on small or moderately padded problems radix-2 usually wins.
        compact = padded
        do j=1,size(n)
            if (.not. selected(j)) cycle
            compact(j) = fast_fft_length(n(j)+m(j)-1)
            if (compact(j) == 0) then
                ok = .false.
                return
            end if
        end do
        best_cost = convolution_padding_cost(padded,selected)
        cost = convolution_padding_cost(compact,selected)
        if (cost < 0.88_dp*best_cost) then
            padded = compact
            best_cost = cost
        end if
        ! Consider single-axis reductions too, which avoid paying the
        ! recursive mixed-radix penalty on every transformed dimension.
        do j=1,size(n)
            if (.not. selected(j)) cycle
            if (padded(j) == compact(j)) cycle
            trial = padded
            trial(j) = compact(j)
            cost = convolution_padding_cost(trial,selected)
            if (cost < 0.88_dp*best_cost) then
                padded = trial
                best_cost = cost
            end if
        end do
    end subroutine fft_nd_layout

    subroutine fft_nd2(x, selected, backward, plans)
        complex(dp), intent(inout) :: x(:,:) !! transform data in place
        logical, intent(in) :: selected(2) !! transform axes
        logical, intent(in) :: backward !! inverse FFT if true
        type(fft_plan), intent(inout) :: plans(2) !! one plan per axis
        integer :: fft_status
        integer :: i

        if (selected(1)) then
            do i = 1, size(x,2)
                if (backward) then
                    call plans(1)%inverse(x(:,i),fft_status)
                else
                    call plans(1)%execute(x(:,i),fft_status)
                end if
            end do
        end if
        if (selected(2)) then
            do i = 1, size(x,1)
                if (backward) then
                    call plans(2)%inverse(x(i,:),fft_status)
                else
                    call plans(2)%execute(x(i,:),fft_status)
                end if
            end do
        end if
    end subroutine fft_nd2

    subroutine fft_nd3(x, selected, backward, plans)
        complex(dp), intent(inout) :: x(:,:,:) !! transform data in place
        logical, intent(in) :: selected(3) !! transform axes
        logical, intent(in) :: backward !! inverse FFT if true
        type(fft_plan), intent(inout) :: plans(3) !! one plan per axis
        integer :: fft_status
        integer :: i, j

        if (selected(1)) then
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    if (backward) then
                        call plans(1)%inverse(x(:,i,j),fft_status)
                    else
                        call plans(1)%execute(x(:,i,j),fft_status)
                    end if
                end do
            end do
        end if
        if (selected(2)) then
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    if (backward) then
                        call plans(2)%inverse(x(i,:,j),fft_status)
                    else
                        call plans(2)%execute(x(i,:,j),fft_status)
                    end if
                end do
            end do
        end if
        if (selected(3)) then
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    if (backward) then
                        call plans(3)%inverse(x(i,j,:),fft_status)
                    else
                        call plans(3)%execute(x(i,j,:),fft_status)
                    end if
                end do
            end do
        end if
    end subroutine fft_nd3

    function fftconvolve_core2(a,b,mode,axes,status) result(y)
        complex(dp), intent(in) :: a(:,:) !! first complex matrix
        complex(dp), intent(in) :: b(:,:) !! second complex matrix
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(in), optional :: axes(:) !! one-based convolution axes
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: fx(:,:), fy(:,:)
        integer :: n(2),m(2),p(2),first(2),cnt(2),i,j,stat
        logical :: selected(2),ok
        type(fft_plan) :: plans(2)
        integer :: ax, fft_status
        character(:), allocatable :: selected_mode

        if (present(status)) status = convolution_status_invalid_input
        allocate(y(0,0))
        if (.not. all(ieee_is_finite(real(a,dp))) .or. .not. all(ieee_is_finite(aimag(a))) .or. &
            .not. all(ieee_is_finite(real(b,dp))) .or. .not. all(ieee_is_finite(aimag(b)))) return
        call fft_axis_set(2,axes,selected,ok)
        if (.not. ok) return
        selected_mode = 'full'
        if (present(mode)) selected_mode = trim(mode)
        n = shape(a)
        m = shape(b)
        call fft_nd_layout(n,m,selected,selected_mode,p,first,cnt,ok)
        if (.not. ok) return
        do ax = 1,2
            if (.not. selected(ax)) cycle
            call plans(ax)%initialize(p(ax),fft_status)
            if (fft_status /= 0) return
        end do
        allocate(fx(p(1),p(2)),fy(p(1),p(2)),stat=stat)
        if (stat /= 0) return
        fx = (0.0_dp,0.0_dp)
        fy = (0.0_dp,0.0_dp)
        do j = 1, p(2)
            do i = 1, p(1)
                if (any(([i,j] > n) .and. selected)) cycle
                fx(i,j) = a(min(i,n(1)),min(j,n(2)))
            end do
        end do
        do j = 1, p(2)
            do i = 1, p(1)
                if (any(([i,j] > m) .and. selected)) cycle
                fy(i,j) = b(min(i,m(1)),min(j,m(2)))
            end do
        end do
        call fft_nd2(fx,selected,.false.,plans)
        call fft_nd2(fy,selected,.false.,plans)
        fx = fx*fy
        call fft_nd2(fx,selected,.true.,plans)
        y = fx(first(1):first(1)+cnt(1)-1,first(2):first(2)+cnt(2)-1)
        if (present(status)) status = convolution_status_success
    end function fftconvolve_core2

    function fftconvolve_core3(a,b,mode,axes,status) result(y)
        complex(dp), intent(in) :: a(:,:,:) !! first complex rank-three array
        complex(dp), intent(in) :: b(:,:,:) !! second complex rank-three array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(in), optional :: axes(:) !! one-based convolution axes
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: fx(:,:,:), fy(:,:,:)
        integer :: n(3),m(3),p(3),first(3),cnt(3),i,j,k,stat
        logical :: selected(3),ok
        type(fft_plan) :: plans(3)
        integer :: ax, fft_status
        character(:), allocatable :: selected_mode

        if (present(status)) status = convolution_status_invalid_input
        allocate(y(0,0,0))
        if (.not. all(ieee_is_finite(real(a,dp))) .or. .not. all(ieee_is_finite(aimag(a))) .or. &
            .not. all(ieee_is_finite(real(b,dp))) .or. .not. all(ieee_is_finite(aimag(b)))) return
        call fft_axis_set(3,axes,selected,ok)
        if (.not. ok) return
        selected_mode = 'full'
        if (present(mode)) selected_mode = trim(mode)
        n = shape(a)
        m = shape(b)
        call fft_nd_layout(n,m,selected,selected_mode,p,first,cnt,ok)
        if (.not. ok) return
        do ax = 1,3
            if (.not. selected(ax)) cycle
            call plans(ax)%initialize(p(ax),fft_status)
            if (fft_status /= 0) return
        end do
        allocate(fx(p(1),p(2),p(3)),fy(p(1),p(2),p(3)),stat=stat)
        if (stat /= 0) return
        fx = (0.0_dp,0.0_dp)
        fy = (0.0_dp,0.0_dp)
        do k = 1, p(3)
            do j = 1, p(2)
                do i = 1, p(1)
                    if (any(([i,j,k] > n) .and. selected)) cycle
                    fx(i,j,k) = a(min(i,n(1)),min(j,n(2)),min(k,n(3)))
                end do
            end do
        end do
        do k = 1, p(3)
            do j = 1, p(2)
                do i = 1, p(1)
                    if (any(([i,j,k] > m) .and. selected)) cycle
                    fy(i,j,k) = b(min(i,m(1)),min(j,m(2)),min(k,m(3)))
                end do
            end do
        end do
        call fft_nd3(fx,selected,.false.,plans)
        call fft_nd3(fy,selected,.false.,plans)
        fx = fx*fy
        call fft_nd3(fx,selected,.true.,plans)
        y = fx(first(1):first(1)+cnt(1)-1,first(2):first(2)+cnt(2)-1, &
            first(3):first(3)+cnt(3)-1)
        if (present(status)) status = convolution_status_success
    end function fftconvolve_core3

    function fftconvolve_r2_real(a,b,mode,status,axes) result(y)
        real(dp), intent(in) :: a(:,:) !! first input array
        real(dp), intent(in) :: b(:,:) !! second input array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based convolution axes, default all
        real(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: work(:,:)

        work = fftconvolve_core2(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),mode,axes,status)
        y = real(work,dp)
    end function fftconvolve_r2_real

    function fftconvolve_r2_complex(a,b,mode,status,axes) result(y)
        complex(dp), intent(in) :: a(:,:) !! first input array
        complex(dp), intent(in) :: b(:,:) !! second input array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based convolution axes, default all
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: work(:,:)

        work = fftconvolve_core2(a,b,mode,axes,status)
        y = work
    end function fftconvolve_r2_complex

    function fftconvolve_r2_real_complex(a,b,mode,status,axes) result(y)
        real(dp), intent(in) :: a(:,:) !! first input array
        complex(dp), intent(in) :: b(:,:) !! second input array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based convolution axes, default all
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: work(:,:)

        work = fftconvolve_core2(cmplx(a,0.0_dp,dp),b,mode,axes,status)
        y = work
    end function fftconvolve_r2_real_complex

    function fftconvolve_r2_complex_real(a,b,mode,status,axes) result(y)
        complex(dp), intent(in) :: a(:,:) !! first input array
        real(dp), intent(in) :: b(:,:) !! second input array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based convolution axes, default all
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: work(:,:)

        work = fftconvolve_core2(a,cmplx(b,0.0_dp,dp),mode,axes,status)
        y = work
    end function fftconvolve_r2_complex_real

    function fftconvolve_r3_real(a,b,mode,status,axes) result(y)
        real(dp), intent(in) :: a(:,:,:) !! first input array
        real(dp), intent(in) :: b(:,:,:) !! second input array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based convolution axes, default all
        real(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: work(:,:,:)

        work = fftconvolve_core3(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),mode,axes,status)
        y = real(work,dp)
    end function fftconvolve_r3_real

    function fftconvolve_r3_complex(a,b,mode,status,axes) result(y)
        complex(dp), intent(in) :: a(:,:,:) !! first input array
        complex(dp), intent(in) :: b(:,:,:) !! second input array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based convolution axes, default all
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: work(:,:,:)

        work = fftconvolve_core3(a,b,mode,axes,status)
        y = work
    end function fftconvolve_r3_complex

    function fftconvolve_r3_real_complex(a,b,mode,status,axes) result(y)
        real(dp), intent(in) :: a(:,:,:) !! first input array
        complex(dp), intent(in) :: b(:,:,:) !! second input array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based convolution axes, default all
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: work(:,:,:)

        work = fftconvolve_core3(cmplx(a,0.0_dp,dp),b,mode,axes,status)
        y = work
    end function fftconvolve_r3_real_complex

    function fftconvolve_r3_complex_real(a,b,mode,status,axes) result(y)
        complex(dp), intent(in) :: a(:,:,:) !! first input array
        real(dp), intent(in) :: b(:,:,:) !! second input array
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axes(:) !! one-based convolution axes, default all
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: work(:,:,:)

        work = fftconvolve_core3(a,cmplx(b,0.0_dp,dp),mode,axes,status)
        y = work
    end function fftconvolve_r3_complex_real



    ! The ND entry point uses caller-owned output because the rank of an
    ! allocatable function result must be known at compile time in Fortran.
    ! We store column-major arrays as rank-one contiguous views, and perform
    ! separable FFTs along the selected dimensions. This also permits rank 15.
    pure function nd_checked_product(extents) result(total)
        integer, intent(in) :: extents(:) !! positive dimensions to multiply
        integer(int64) :: total
        integer :: j

        total = 1_int64
        do j = 1,size(extents)
            if (extents(j) < 1) then
                total = 0_int64
                return
            end if
            if (total > huge(total)/int(extents(j),int64)) then
                total = 0_int64
                return
            end if
            total = total*int(extents(j),int64)
        end do
    end function nd_checked_product

    subroutine fft_flat_axis(x, dims, axis, inverse, line, plan)
        complex(dp), intent(inout) :: x(:) !! packed column-major transform data
        integer, intent(in) :: dims(:) !! padded dimensions of the data
        integer, intent(in) :: axis !! one-based transform axis
        logical, intent(in) :: inverse !! true for inverse FFT
        complex(dp), intent(inout) :: line(:) !! scratch, at least dims(axis)
        type(fft_plan), intent(inout) :: plan !! initialized axis plan
        integer :: fft_status
        integer(int64) :: inner, outer, stride, block, base, j, k
        integer :: n

        n = dims(axis)
        stride = nd_checked_product(dims(:axis-1))
        outer = nd_checked_product(dims(axis+1:))
        block = stride*int(n,int64)
        do j = 0_int64,outer-1_int64
            do inner = 0_int64,stride-1_int64
                base = j*block+inner+1_int64
                do k = 0_int64,int(n,int64)-1_int64
                    line(k+1) = x(base+k*stride)
                end do
                if (inverse) then
                    call plan%inverse(line(:n),fft_status)
                else
                    call plan%execute(line(:n),fft_status)
                end if
                do k = 0_int64,int(n,int64)-1_int64
                    x(base+k*stride) = line(k+1)
                end do
            end do
        end do
    end subroutine fft_flat_axis

    subroutine fftconvolve_flat(a, b, n, m, out_dims, result, mode, axes, status)
        complex(dp), intent(in) :: a(:) !! packed first input, column-major
        complex(dp), intent(in) :: b(:) !! packed second input, column-major
        integer, intent(in) :: n(:) !! shape of first input
        integer, intent(in) :: m(:) !! shape of second input
        integer, intent(in) :: out_dims(:) !! shape of caller-provided output
        complex(dp), allocatable, intent(out) :: result(:) !! packed output on success
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, intent(out) :: status !! zero for success; otherwise invalid input
        integer :: r, ax, istat, p(size(n)), first(size(n)), cnt(size(n))
        integer :: j
        integer(int64) :: np, ny, k, q, qout, ka, kb, sa, sb, pos, pstride
        integer(int64) :: coord
        complex(dp), allocatable :: fx(:), fy(:), line(:)
        type(fft_plan) :: plans(size(n))
        integer :: fft_status
        logical :: selected(size(n)), ok, usea, useb
        character(:), allocatable :: selected_mode

        status = convolution_status_invalid_input
        r = size(n)
        if (r < 1 .or. r > 15) return
        if (size(m) /= r .or. size(out_dims) /= r) return
        if (nd_checked_product(n) /= size(a,kind=int64)) return
        if (nd_checked_product(m) /= size(b,kind=int64)) return
        if (.not. all(ieee_is_finite(real(a,dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        if (.not. all(ieee_is_finite(real(b,dp)))) return
        if (.not. all(ieee_is_finite(aimag(b)))) return
        call fft_axis_set(r,axes,selected,ok)
        if (.not. ok) return
        selected_mode = 'full'
        if (present(mode)) selected_mode = trim(mode)
        call fft_nd_layout(n,m,selected,selected_mode,p,first,cnt,ok)
        if (.not. ok) return
        if (any(out_dims /= cnt)) return
        np = nd_checked_product(p)
        ny = nd_checked_product(cnt)
        if (np == 0_int64 .or. ny == 0_int64) return
        do ax = 1,r
            if (.not. selected(ax)) cycle
            call plans(ax)%initialize(p(ax),fft_status)
            if (fft_status /= 0) return
        end do
        allocate(fx(np),fy(np),line(maxval(p,mask=selected)),stat=istat)
        if (istat /= 0) return
        fx = (0.0_dp,0.0_dp)
        fy = (0.0_dp,0.0_dp)
        ! Zero-pad on transformed axes and broadcast singleton dimensions on
        ! untransformed axes. Every address is computed in Fortran order.
        do k = 0_int64,np-1_int64
            q = k
            ka = 1_int64
            kb = 1_int64
            sa = 1_int64
            sb = 1_int64
            usea = .true.
            useb = .true.
            do j = 1,r
                coord = mod(q,int(p(j),int64))+1_int64
                q = q/int(p(j),int64)
                if (selected(j)) then
                    if (coord > int(n(j),int64)) usea = .false.
                    if (coord > int(m(j),int64)) useb = .false.
                end if
                ka = ka + (min(coord,int(n(j),int64))-1_int64)*sa
                kb = kb + (min(coord,int(m(j),int64))-1_int64)*sb
                sa = sa*int(n(j),int64)
                sb = sb*int(m(j),int64)
            end do
            if (usea) fx(k+1) = a(ka)
            if (useb) fy(k+1) = b(kb)
        end do
        do ax = 1,r
            if (.not. selected(ax)) cycle
            call fft_flat_axis(fx,p,ax,.false.,line,plans(ax))
            call fft_flat_axis(fy,p,ax,.false.,line,plans(ax))
        end do
        fx = fx*fy
        do ax = 1,r
            if (.not. selected(ax)) cycle
            call fft_flat_axis(fx,p,ax,.true.,line,plans(ax))
        end do
        allocate(result(ny),stat=istat)
        if (istat /= 0) return
        do k = 0_int64,ny-1_int64
            qout = k
            pos = 1_int64
            pstride = 1_int64
            do j = 1,r
                coord = mod(qout,int(cnt(j),int64))
                qout = qout/int(cnt(j),int64)
                pos = pos+(coord+int(first(j)-1,int64))*pstride
                pstride = pstride*int(p(j),int64)
            end do
            result(k+1) = fx(pos)
        end do
        status = convolution_status_success
    end subroutine fftconvolve_flat

    ! Populate one padded FFT workspace directly from a packed operand.
    ! starts(j) is one-based in the source; tile_dims(j) is the number
    ! of available samples from that start on each transformed axis.
    ! Untransformed singleton axes are broadcast, not zero-padded.
    subroutine oa_pad_flat(src,dims,padded,selected,starts,tile_dims,work)
        complex(dp), intent(in) :: src(:) !! source in Fortran storage order
        integer, intent(in) :: dims(:) !! complete source shape
        integer, intent(in) :: padded(:) !! shape of FFT workspace
        logical, intent(in) :: selected(:) !! transformed dimensions
        integer, intent(in) :: starts(:) !! first source element on each axis
        integer, intent(in) :: tile_dims(:) !! active tile shape
        complex(dp), intent(out) :: work(:) !! padded, broadcast workspace
        integer(int64) :: k, q, coord, source_coord, pos, stride
        integer :: j
        logical :: inside

        work = (0.0_dp,0.0_dp)
        do k = 0_int64,size(work,kind=int64)-1_int64
            q = k
            pos = 1_int64
            stride = 1_int64
            inside = .true.
            do j = 1,size(dims)
                coord = mod(q,int(padded(j),int64))+1_int64
                q = q/int(padded(j),int64)
                if (selected(j)) then
                    if (coord > int(tile_dims(j),int64)) then
                        inside = .false.
                        exit
                    end if
                    source_coord = coord+int(starts(j)-1,int64)
                else
                    source_coord = min(coord,int(dims(j),int64))
                end if
                pos = pos+(source_coord-1_int64)*stride
                stride = stride*int(dims(j),int64)
            end do
            if (inside) work(k+1_int64) = src(pos)
        end do
    end subroutine oa_pad_flat

    ! Multidimensional overlap-add using a fixed, reusable kernel spectrum.
    ! For unequal input sizes, partition multiple transformed dimensions of
    ! one operand.  A single FFT of the other operand is then reused across
    ! every tile, and contributions are cropped as they are accumulated.
    ! Splitting dimensions on opposing operands is deliberately not combined:
    ! doing so would invalidate reuse of the unchanged kernel spectrum.
    subroutine oaconvolve_flat(a,b,n,m,out_dims,result,mode,axes,status)
        complex(dp), intent(in) :: a(:) !! packed first operand, column-major
        complex(dp), intent(in) :: b(:) !! packed second operand, column-major
        integer, intent(in) :: n(:) !! first operand extents
        integer, intent(in) :: m(:) !! second operand extents
        integer, intent(in) :: out_dims(:) !! required output extents
        complex(dp), allocatable, intent(out) :: result(:) !! packed convolution
        character(*), intent(in), optional :: mode !! full, same or valid
        integer, intent(in), optional :: axes(:) !! one-based convolved axes
        integer, intent(out) :: status !! success or invalid input
        integer :: r, ax, j, code
        integer :: full_p(size(n)), first(size(n)), cnt(size(n))
        integer :: long_dims(size(n)), kernel_dims(size(n))
        integer :: tile_dims(size(n)), actual_dims(size(n))
        integer :: padded(size(n)), tile_first(size(n)), tile_count(size(n))
        integer :: starts(size(n)), chunks(size(n)), last_valid(size(n))
        integer :: ones(size(n))
        integer(int64) :: ny, nwork, k, q, pos, stride, coord, offset
        real(dp) :: score_a, score_b
        logical :: selected(size(n)), split(size(n)), ok, split_first, inside
        character(:), allocatable :: selected_mode
        complex(dp), allocatable :: fx(:), fy(:), line(:)
        type(fft_plan) :: plans(size(n))
        integer :: fft_status

        status = convolution_status_invalid_input
        r = size(n)
        if (r < 1 .or. r > 15) return
        if (size(m) /= r .or. size(out_dims) /= r) return
        if (nd_checked_product(n) /= size(a,kind=int64)) return
        if (nd_checked_product(m) /= size(b,kind=int64)) return
        if (.not. all(ieee_is_finite(real(a,dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        if (.not. all(ieee_is_finite(real(b,dp)))) return
        if (.not. all(ieee_is_finite(aimag(b)))) return
        call fft_axis_set(r,axes,selected,ok)
        if (.not. ok) return
        selected_mode = 'full'
        if (present(mode)) selected_mode = trim(mode)
        call fft_nd_layout(n,m,selected,selected_mode,full_p,first,cnt,ok)
        if (.not. ok) return
        if (any(out_dims /= cnt)) return
        ny = nd_checked_product(cnt)
        if (ny == 0_int64) return

        ! Select the side with the greater combined asymmetric length ratios.
        ! Only sufficiently long axes benefit from the minimum 64-point FFT.
        score_a = 0.0_dp
        score_b = 0.0_dp
        do ax = 1,r
            if (.not. selected(ax)) cycle
            if (oa_chunk_length(n(ax),m(ax)) > 0) then
                score_a = score_a+log(real(n(ax),dp)/real(m(ax),dp))
            end if
            if (oa_chunk_length(m(ax),n(ax)) > 0) then
                score_b = score_b+log(real(m(ax),dp)/real(n(ax),dp))
            end if
        end do
        if (score_a+score_b <= 0.0_dp) then
            call fftconvolve_flat(a,b,n,m,out_dims,result,mode,axes,status)
            return
        end if
        split_first = score_a >= score_b
        if (split_first) then
            long_dims = n
            kernel_dims = m
        else
            long_dims = m
            kernel_dims = n
        end if
        chunks = long_dims
        split = .false.
        tile_dims = long_dims
        do ax = 1,r
            if (.not. selected(ax)) cycle
            chunks(ax) = oa_chunk_length(long_dims(ax),kernel_dims(ax))
            if (chunks(ax) == 0) then
                chunks(ax) = long_dims(ax)
            else
                tile_dims(ax) = chunks(ax)
                split(ax) = .true.
            end if
        end do
        if (.not. any(split)) then
            call fftconvolve_flat(a,b,n,m,out_dims,result,mode,axes,status)
            return
        end if

        ! The FFT layout remains constant for all tiles. The last (shorter)
        ! tiles are zero-padded to the same size rather than retransforming
        ! the kernel for every different block length.
        call fft_nd_layout(tile_dims,kernel_dims,selected,'full',padded,tile_first,tile_count,ok)
        if (.not. ok) return
        nwork = nd_checked_product(padded)
        if (nwork == 0_int64) return
        do ax = 1,r
            if (.not. selected(ax)) cycle
            call plans(ax)%initialize(padded(ax),fft_status)
            if (fft_status /= 0) return
        end do
        allocate(fx(nwork),fy(nwork),line(maxval(padded,mask=selected)), &
            result(ny),stat=code)
        if (code /= 0) return
        result = (0.0_dp,0.0_dp)
        starts = 1
        ones = 1
        if (split_first) then
            call oa_pad_flat(b,m,padded,selected,ones,m,fy)
        else
            call oa_pad_flat(a,n,padded,selected,ones,n,fy)
        end if
        do ax = 1,r
            if (selected(ax)) call fft_flat_axis(fy,padded,ax,.false.,line,plans(ax))
        end do

        do
            actual_dims = long_dims
            where (split) actual_dims = min(chunks,long_dims-starts+1)
            if (split_first) then
                call oa_pad_flat(a,n,padded,selected,starts,actual_dims,fx)
            else
                call oa_pad_flat(b,m,padded,selected,starts,actual_dims,fx)
            end if
            do ax = 1,r
                if (selected(ax)) call fft_flat_axis(fx,padded,ax,.false.,line,plans(ax))
            end do
            fx = fx*fy
            do ax = 1,r
                if (selected(ax)) call fft_flat_axis(fx,padded,ax,.true.,line,plans(ax))
            end do

            last_valid = padded
            where (selected) last_valid = actual_dims+kernel_dims-1
            do k = 0_int64,nwork-1_int64
                q = k
                pos = 1_int64
                stride = 1_int64
                inside = .true.
                do j = 1,r
                    coord = mod(q,int(padded(j),int64))+1_int64
                    q = q/int(padded(j),int64)
                    if (coord > int(last_valid(j),int64)) then
                        inside = .false.
                        exit
                    end if
                    if (split(j)) coord = coord+int(starts(j)-1,int64)
                    offset = coord-int(first(j),int64)
                    if (offset < 0_int64 .or. offset >= int(cnt(j),int64)) then
                        inside = .false.
                        exit
                    end if
                    pos = pos+offset*stride
                    stride = stride*int(cnt(j),int64)
                end do
                if (inside) result(pos) = result(pos)+fx(k+1_int64)
            end do
            ! Advance a mixed-radix tile counter over every split dimension.
            do j = 1,r
                if (.not. split(j)) cycle
                starts(j) = starts(j)+chunks(j)
                if (starts(j) <= long_dims(j)) exit
                starts(j) = 1
            end do
            if (j > r) exit
        end do
        status = convolution_status_success
    end subroutine oaconvolve_flat

    subroutine fftconvolve_nd_rr(in1,in2,y,mode,axes,status)
        real(dp), contiguous, target, intent(in) :: in1(..) !! first rank 1-15 array
        real(dp), contiguous, target, intent(in) :: in2(..) !! second rank 1-15 array
        real(dp), contiguous, target, intent(out) :: y(..) !! requested output shape
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, optional, intent(out) :: status !! zero for success
        real(dp), pointer :: a(:) => null()
        real(dp), pointer :: b(:) => null()
        real(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = convolution_status_invalid_input
        if (rank(y) > 0) then
            if (size(y,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
                out = nanval
            end if
        end if
        if (rank(in1) < 1 .or. rank(in2) /= rank(in1)) return
        if (rank(y) /= rank(in1) .or. rank(y) > 15) return
        if (size(in1,kind=int64) == 0_int64) return
        if (size(in2,kind=int64) == 0_int64) return
        if (size(y,kind=int64) == 0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call fftconvolve_flat(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            shape(in1),shape(in2),shape(y),tmp,mode,axes,code)
        if (code /= convolution_status_success) return
        out = real(tmp,dp)
        if (present(status)) status = convolution_status_success
    end subroutine fftconvolve_nd_rr

    subroutine fftconvolve_nd_cc(in1,in2,y,mode,axes,status)
        complex(dp), contiguous, target, intent(in) :: in1(..) !! first rank 1-15 array
        complex(dp), contiguous, target, intent(in) :: in2(..) !! second rank 1-15 array
        complex(dp), contiguous, target, intent(out) :: y(..) !! requested output shape
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, optional, intent(out) :: status !! zero for success
        complex(dp), pointer :: a(:) => null()
        complex(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = convolution_status_invalid_input
        if (rank(y) > 0) then
            if (size(y,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
                out = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(in1) < 1 .or. rank(in2) /= rank(in1)) return
        if (rank(y) /= rank(in1) .or. rank(y) > 15) return
        if (size(in1,kind=int64) == 0_int64) return
        if (size(in2,kind=int64) == 0_int64) return
        if (size(y,kind=int64) == 0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call fftconvolve_flat(a,b, &
            shape(in1),shape(in2),shape(y),tmp,mode,axes,code)
        if (code /= convolution_status_success) return
        out = tmp
        if (present(status)) status = convolution_status_success
    end subroutine fftconvolve_nd_cc

    subroutine fftconvolve_nd_rc(in1,in2,y,mode,axes,status)
        real(dp), contiguous, target, intent(in) :: in1(..) !! first rank 1-15 array
        complex(dp), contiguous, target, intent(in) :: in2(..) !! second rank 1-15 array
        complex(dp), contiguous, target, intent(out) :: y(..) !! requested output shape
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, optional, intent(out) :: status !! zero for success
        real(dp), pointer :: a(:) => null()
        complex(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = convolution_status_invalid_input
        if (rank(y) > 0) then
            if (size(y,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
                out = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(in1) < 1 .or. rank(in2) /= rank(in1)) return
        if (rank(y) /= rank(in1) .or. rank(y) > 15) return
        if (size(in1,kind=int64) == 0_int64) return
        if (size(in2,kind=int64) == 0_int64) return
        if (size(y,kind=int64) == 0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call fftconvolve_flat(cmplx(a,0.0_dp,dp),b, &
            shape(in1),shape(in2),shape(y),tmp,mode,axes,code)
        if (code /= convolution_status_success) return
        out = tmp
        if (present(status)) status = convolution_status_success
    end subroutine fftconvolve_nd_rc

    subroutine fftconvolve_nd_cr(in1,in2,y,mode,axes,status)
        complex(dp), contiguous, target, intent(in) :: in1(..) !! first rank 1-15 array
        real(dp), contiguous, target, intent(in) :: in2(..) !! second rank 1-15 array
        complex(dp), contiguous, target, intent(out) :: y(..) !! requested output shape
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, optional, intent(out) :: status !! zero for success
        complex(dp), pointer :: a(:) => null()
        real(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = convolution_status_invalid_input
        if (rank(y) > 0) then
            if (size(y,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
                out = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(in1) < 1 .or. rank(in2) /= rank(in1)) return
        if (rank(y) /= rank(in1) .or. rank(y) > 15) return
        if (size(in1,kind=int64) == 0_int64) return
        if (size(in2,kind=int64) == 0_int64) return
        if (size(y,kind=int64) == 0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call fftconvolve_flat(a,cmplx(b,0.0_dp,dp), &
            shape(in1),shape(in2),shape(y),tmp,mode,axes,code)
        if (code /= convolution_status_success) return
        out = tmp
        if (present(status)) status = convolution_status_success
    end subroutine fftconvolve_nd_cr

    subroutine oaconvolve_nd_rr(in1,in2,y,mode,axes,status)
        real(dp), contiguous, target, intent(in) :: in1(..) !! first rank 1-15 array
        real(dp), contiguous, target, intent(in) :: in2(..) !! second rank 1-15 array
        real(dp), contiguous, target, intent(out) :: y(..) !! requested output shape
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, optional, intent(out) :: status !! zero for success
        real(dp), pointer :: a(:) => null()
        real(dp), pointer :: b(:) => null()
        real(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = convolution_status_invalid_input
        if (rank(y) > 0) then
            if (size(y,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
                out = nanval
            end if
        end if
        if (rank(in1) < 1 .or. rank(in2) /= rank(in1)) return
        if (rank(y) /= rank(in1) .or. rank(y) > 15) return
        if (size(in1,kind=int64) == 0_int64) return
        if (size(in2,kind=int64) == 0_int64) return
        if (size(y,kind=int64) == 0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call oaconvolve_flat(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            shape(in1),shape(in2),shape(y),tmp,mode,axes,code)
        if (code /= convolution_status_success) return
        out = real(tmp,dp)
        if (present(status)) status = convolution_status_success
    end subroutine oaconvolve_nd_rr

    subroutine oaconvolve_nd_cc(in1,in2,y,mode,axes,status)
        complex(dp), contiguous, target, intent(in) :: in1(..) !! first rank 1-15 array
        complex(dp), contiguous, target, intent(in) :: in2(..) !! second rank 1-15 array
        complex(dp), contiguous, target, intent(out) :: y(..) !! requested output shape
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, optional, intent(out) :: status !! zero for success
        complex(dp), pointer :: a(:) => null()
        complex(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = convolution_status_invalid_input
        if (rank(y) > 0) then
            if (size(y,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
                out = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(in1) < 1 .or. rank(in2) /= rank(in1)) return
        if (rank(y) /= rank(in1) .or. rank(y) > 15) return
        if (size(in1,kind=int64) == 0_int64) return
        if (size(in2,kind=int64) == 0_int64) return
        if (size(y,kind=int64) == 0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call oaconvolve_flat(a,b, &
            shape(in1),shape(in2),shape(y),tmp,mode,axes,code)
        if (code /= convolution_status_success) return
        out = tmp
        if (present(status)) status = convolution_status_success
    end subroutine oaconvolve_nd_cc

    subroutine oaconvolve_nd_rc(in1,in2,y,mode,axes,status)
        real(dp), contiguous, target, intent(in) :: in1(..) !! first rank 1-15 array
        complex(dp), contiguous, target, intent(in) :: in2(..) !! second rank 1-15 array
        complex(dp), contiguous, target, intent(out) :: y(..) !! requested output shape
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, optional, intent(out) :: status !! zero for success
        real(dp), pointer :: a(:) => null()
        complex(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = convolution_status_invalid_input
        if (rank(y) > 0) then
            if (size(y,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
                out = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(in1) < 1 .or. rank(in2) /= rank(in1)) return
        if (rank(y) /= rank(in1) .or. rank(y) > 15) return
        if (size(in1,kind=int64) == 0_int64) return
        if (size(in2,kind=int64) == 0_int64) return
        if (size(y,kind=int64) == 0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call oaconvolve_flat(cmplx(a,0.0_dp,dp),b, &
            shape(in1),shape(in2),shape(y),tmp,mode,axes,code)
        if (code /= convolution_status_success) return
        out = tmp
        if (present(status)) status = convolution_status_success
    end subroutine oaconvolve_nd_rc

    subroutine oaconvolve_nd_cr(in1,in2,y,mode,axes,status)
        complex(dp), contiguous, target, intent(in) :: in1(..) !! first rank 1-15 array
        real(dp), contiguous, target, intent(in) :: in2(..) !! second rank 1-15 array
        complex(dp), contiguous, target, intent(out) :: y(..) !! requested output shape
        character(*), optional, intent(in) :: mode !! full (default), same, valid
        integer, optional, intent(in) :: axes(:) !! one-based axes, default all
        integer, optional, intent(out) :: status !! zero for success
        complex(dp), pointer :: a(:) => null()
        real(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = convolution_status_invalid_input
        if (rank(y) > 0) then
            if (size(y,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
                out = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(in1) < 1 .or. rank(in2) /= rank(in1)) return
        if (rank(y) /= rank(in1) .or. rank(y) > 15) return
        if (size(in1,kind=int64) == 0_int64) return
        if (size(in2,kind=int64) == 0_int64) return
        if (size(y,kind=int64) == 0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call oaconvolve_flat(a,cmplx(b,0.0_dp,dp), &
            shape(in1),shape(in2),shape(y),tmp,mode,axes,code)
        if (code /= convolution_status_success) return
        out = tmp
        if (present(status)) status = convolution_status_success
    end subroutine oaconvolve_nd_cr

end module scifort_signal_fftconvolve
