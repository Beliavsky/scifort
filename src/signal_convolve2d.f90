! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Direct two-dimensional convolution and cross-correlation.
! Indices and output centering follow scipy.signal.convolve2d/correlate2d.
module scifort_signal_convolve2d
    use scifort_kinds, only : dp
    use scifort_signal_convolve, only : convolution_status_success, convolution_status_invalid_input
    implicit none
    private

    interface convolve2d
        module procedure convolve2d_real
        module procedure convolve2d_complex
        module procedure convolve2d_real_complex
        module procedure convolve2d_complex_real
    end interface

    interface correlate2d
        module procedure correlate2d_real
        module procedure correlate2d_complex
        module procedure correlate2d_real_complex
        module procedure correlate2d_complex_real
    end interface

    public :: convolve2d, correlate2d

contains

    pure subroutine layout(n, m, mode, correlation, first, count, ok)
        integer, intent(in) :: n(2) !! dimensions of the first array
        integer, intent(in) :: m(2) !! dimensions of the second array
        character(len=*), intent(in) :: mode !! full, same or valid
        logical, intent(in) :: correlation !! true for two-dimensional correlation
        integer, intent(out) :: first(2) !! one-based starting index in full output
        integer, intent(out) :: count(2) !! number of output elements per axis
        logical, intent(out) :: ok !! input and mode validity

        ok = all(n > 0) .and. all(m > 0)
        first = 1
        count = 0
        if (.not. ok) return
        select case (mode)
        case ('full')
            count = n + m - 1
        case ('same')
            if (correlation) then
                ! Unlike 1D correlate, scipy.signal.correlate2d centers even
                ! kernels with an offset of floor(m/2), not floor((m-1)/2).
                first = 1 + m / 2
            else
                first = 1 + (m - 1) / 2
            end if
            count = n
        case ('valid')
            ! The larger operand must dominate along every axis.
            ok = all(n >= m) .or. all(m >= n)
            if (.not. ok) return
            first = min(n, m)
            count = abs(n - m) + 1
        case default
            ok = .false.
        end select
    end subroutine layout

    pure integer function extended_index(k, n, boundary) result(mapped)
        integer, intent(in) :: k !! one-based possibly exterior coordinate
        integer, intent(in) :: n !! positive axis length
        character(len=*), intent(in) :: boundary !! wrap or symm
        integer :: j

        select case (boundary)
        case ('wrap')
            mapped = modulo(k - 1, n) + 1
        case ('symm')
            j = modulo(k - 1, 2*n)
            if (j < n) then
                mapped = j + 1
            else
                mapped = 2*n - j
            end if
        case default
            mapped = k
        end select
    end function extended_index

    function convolve2d_real(in1, in2, mode, boundary, fillvalue, status) result(y)
        real(dp), intent(in) :: in1(:,:) !! first input matrix
        real(dp), intent(in) :: in2(:,:) !! second input matrix
        character(len=*), intent(in), optional :: mode !! full (default), same or valid
        character(len=*), intent(in), optional :: boundary !! fill (default), wrap or symm
        real(dp), intent(in), optional :: fillvalue !! outside value when boundary=fill
        integer, intent(out), optional :: status !! zero if successful
        real(dp), allocatable :: y(:,:)

        integer :: n(2), m(2), first(2), count(2), i, j, ki, kj, ii, jj, li, lj
        real(dp) :: fill, sample, total
        logical :: ok, within, use_fill
        character(len=:), allocatable :: smode, sboundary

        if (present(status)) status = convolution_status_invalid_input
        n = shape(in1)
        m = shape(in2)
        smode = 'full'
        if (present(mode)) smode = mode
        sboundary = 'fill'
        if (present(boundary)) sboundary = boundary
        call layout(n, m, smode, .false., first, count, ok)
        ok = ok .and. (sboundary == 'fill' .or. sboundary == 'wrap' .or. sboundary == 'symm')
        if (.not. ok) then
            allocate(y(0,0))
            return
        end if
        fill = 0.0_dp
        if (present(fillvalue)) fill = fillvalue
        use_fill = sboundary == 'fill' .or. smode == 'valid'
        allocate(y(count(1),count(2)))
        y = 0.0_dp
        do j = 1, count(2)
            lj = first(2) + j - 1
            do i = 1, count(1)
                li = first(1) + i - 1
                total = 0.0_dp
                do kj = 1, m(2)
                    jj = lj - kj + 1
                    do ki = 1, m(1)
                        ii = li - ki + 1
                        within = ii >= 1 .and. ii <= n(1) .and. jj >= 1 .and. jj <= n(2)
                        if (within) then
                            sample = in1(ii,jj)
                        else if (use_fill) then
                            if (smode == 'valid') then
                                sample = 0.0_dp
                            else
                                sample = fill
                            end if
                        else
                            sample = in1(extended_index(ii,n(1),sboundary), &
                                         extended_index(jj,n(2),sboundary))
                        end if
                        total = total + sample * in2(ki,kj)
                    end do
                end do
                y(i,j) = total
            end do
        end do
        if (present(status)) status = convolution_status_success
    end function convolve2d_real

    function convolve2d_complex(in1, in2, mode, boundary, fillvalue, status) result(y)
        complex(dp), intent(in) :: in1(:,:) !! first input matrix
        complex(dp), intent(in) :: in2(:,:) !! second input matrix
        character(len=*), intent(in), optional :: mode !! full (default), same or valid
        character(len=*), intent(in), optional :: boundary !! fill (default), wrap or symm
        complex(dp), intent(in), optional :: fillvalue !! outside value when boundary=fill
        integer, intent(out), optional :: status !! zero if successful
        complex(dp), allocatable :: y(:,:)

        integer :: n(2), m(2), first(2), count(2), i, j, ki, kj, ii, jj, li, lj
        complex(dp) :: fill, sample, total
        logical :: ok, within, use_fill
        character(len=:), allocatable :: smode, sboundary

        if (present(status)) status = convolution_status_invalid_input
        n = shape(in1)
        m = shape(in2)
        smode = 'full'
        if (present(mode)) smode = mode
        sboundary = 'fill'
        if (present(boundary)) sboundary = boundary
        call layout(n, m, smode, .false., first, count, ok)
        ok = ok .and. (sboundary == 'fill' .or. sboundary == 'wrap' .or. sboundary == 'symm')
        if (.not. ok) then
            allocate(y(0,0))
            return
        end if
        fill = (0.0_dp,0.0_dp)
        if (present(fillvalue)) fill = fillvalue
        use_fill = sboundary == 'fill' .or. smode == 'valid'
        allocate(y(count(1),count(2)))
        y = (0.0_dp,0.0_dp)
        do j = 1, count(2)
            lj = first(2) + j - 1
            do i = 1, count(1)
                li = first(1) + i - 1
                total = (0.0_dp,0.0_dp)
                do kj = 1, m(2)
                    jj = lj - kj + 1
                    do ki = 1, m(1)
                        ii = li - ki + 1
                        within = ii >= 1 .and. ii <= n(1) .and. jj >= 1 .and. jj <= n(2)
                        if (within) then
                            sample = in1(ii,jj)
                        else if (use_fill) then
                            if (smode == 'valid') then
                                sample = (0.0_dp,0.0_dp)
                            else
                                sample = fill
                            end if
                        else
                            sample = in1(extended_index(ii,n(1),sboundary), &
                                         extended_index(jj,n(2),sboundary))
                        end if
                        total = total + sample * in2(ki,kj)
                    end do
                end do
                y(i,j) = total
            end do
        end do
        if (present(status)) status = convolution_status_success
    end function convolve2d_complex

    function correlate2d_real(in1, in2, mode, boundary, fillvalue, status) result(y)
        real(dp), intent(in) :: in1(:,:) !! first real input matrix
        real(dp), intent(in) :: in2(:,:) !! second real input matrix
        character(len=*), intent(in), optional :: mode !! full, same, valid
        character(len=*), intent(in), optional :: boundary !! fill, wrap, symm
        real(dp), intent(in), optional :: fillvalue !! boundary fill value
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: full(:,:)
        integer :: n(2), m(2), first(2), count(2)
        logical :: ok
        character(len=:), allocatable :: selected

        selected = 'full'
        if (present(mode)) selected = mode
        n = shape(in1)
        m = shape(in2)
        call layout(n, m, selected, .true., first, count, ok)
        if (present(status)) status = convolution_status_invalid_input
        allocate(y(0,0))
        if (.not. ok) return
        if (selected == 'valid') then
            ! No boundary extension may contribute when either operand is
            ! the larger one; in particular, ignore a nonzero fill value.
            full = convolve2d_real(in1, in2(m(1):1:-1,m(2):1:-1), &
                mode='full', boundary='fill', fillvalue=0.0_dp, status=status)
        else
            full = convolve2d_real(in1, in2(m(1):1:-1,m(2):1:-1), &
                mode='full', boundary=boundary, fillvalue=fillvalue, status=status)
        end if
        if (present(status)) then
            if (status /= convolution_status_success) return
        else
            if (size(full) == 0) return
        end if
        y = full(first(1):first(1)+count(1)-1, first(2):first(2)+count(2)-1)
    end function correlate2d_real

    function correlate2d_complex(in1, in2, mode, boundary, fillvalue, status) result(y)
        complex(dp), intent(in) :: in1(:,:) !! first complex input matrix
        complex(dp), intent(in) :: in2(:,:) !! second complex input matrix
        character(len=*), intent(in), optional :: mode !! full, same, valid
        character(len=*), intent(in), optional :: boundary !! fill, wrap, symm
        complex(dp), intent(in), optional :: fillvalue !! boundary fill value
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: full(:,:)
        integer :: n(2), m(2), first(2), count(2)
        logical :: ok
        character(len=:), allocatable :: selected

        selected = 'full'
        if (present(mode)) selected = mode
        n = shape(in1)
        m = shape(in2)
        call layout(n, m, selected, .true., first, count, ok)
        if (present(status)) status = convolution_status_invalid_input
        allocate(y(0,0))
        if (.not. ok) return
        if (selected == 'valid') then
            full = convolve2d_complex(in1, conjg(in2(m(1):1:-1,m(2):1:-1)), &
                mode='full', boundary='fill', fillvalue=(0.0_dp,0.0_dp), status=status)
        else
            full = convolve2d_complex(in1, conjg(in2(m(1):1:-1,m(2):1:-1)), &
                mode='full', boundary=boundary, fillvalue=fillvalue, status=status)
        end if
        if (present(status)) then
            if (status /= convolution_status_success) return
        else
            if (size(full) == 0) return
        end if
        y = full(first(1):first(1)+count(1)-1, first(2):first(2)+count(2)-1)
        ! SciPy swaps operands for valid mode if the second matrix dominates.
        ! For complex inputs it returns the conjugate of the unswapped
        ! correlation (unlike full and same); match that convention.
        if (selected == 'valid' .and. .not. all(n >= m)) y = conjg(y)
    end function correlate2d_complex

    function convolve2d_real_complex(in1, in2, mode, boundary, fillvalue, status) result(y)
        real(dp), intent(in) :: in1(:,:) !! first real input matrix
        complex(dp), intent(in) :: in2(:,:) !! second complex input matrix
        character(len=*), intent(in), optional :: mode !! full, same or valid
        character(len=*), intent(in), optional :: boundary !! fill, wrap or symm
        complex(dp), intent(in), optional :: fillvalue !! optional complex boundary value
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        y = convolve2d_complex(cmplx(in1,kind=dp), in2, mode, boundary, fillvalue, status)
    end function convolve2d_real_complex

    function convolve2d_complex_real(in1, in2, mode, boundary, fillvalue, status) result(y)
        complex(dp), intent(in) :: in1(:,:) !! first complex input matrix
        real(dp), intent(in) :: in2(:,:) !! second real input matrix
        character(len=*), intent(in), optional :: mode !! full, same or valid
        character(len=*), intent(in), optional :: boundary !! fill, wrap or symm
        complex(dp), intent(in), optional :: fillvalue !! optional complex boundary value
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        y = convolve2d_complex(in1, cmplx(in2,kind=dp), mode, boundary, fillvalue, status)
    end function convolve2d_complex_real

    function correlate2d_real_complex(in1, in2, mode, boundary, fillvalue, status) result(y)
        real(dp), intent(in) :: in1(:,:) !! first real input matrix
        complex(dp), intent(in) :: in2(:,:) !! second complex input matrix
        character(len=*), intent(in), optional :: mode !! full, same or valid
        character(len=*), intent(in), optional :: boundary !! fill, wrap or symm
        complex(dp), intent(in), optional :: fillvalue !! optional complex boundary value
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        y = correlate2d_complex(cmplx(in1,kind=dp), in2, mode, boundary, fillvalue, status)
    end function correlate2d_real_complex

    function correlate2d_complex_real(in1, in2, mode, boundary, fillvalue, status) result(y)
        complex(dp), intent(in) :: in1(:,:) !! first complex input matrix
        real(dp), intent(in) :: in2(:,:) !! second real input matrix
        character(len=*), intent(in), optional :: mode !! full, same or valid
        character(len=*), intent(in), optional :: boundary !! fill, wrap or symm
        complex(dp), intent(in), optional :: fillvalue !! optional complex boundary value
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        y = correlate2d_complex(in1, cmplx(in2,kind=dp), mode, boundary, fillvalue, status)
    end function correlate2d_complex_real

end module scifort_signal_convolve2d
