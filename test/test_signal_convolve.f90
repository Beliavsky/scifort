! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy 1.17.0 reference cases for 1D convolution and correlation.
program test_signal_convolve
    use scifort_signal
    implicit none

    integer :: ierr, n, m, i
    integer, allocatable :: lags(:)
    real(dp), allocatable :: y(:)
    real(dp) :: a(3), b(2)
    complex(dp), allocatable :: z(:)
    complex(dp) :: ca(2), cb(2)

    a = [1.0_dp, 2.0_dp, 3.0_dp]
    b = [4.0_dp, 5.0_dp]
    y = convolve(a, b)
    call check_real(y, [4.0_dp, 13.0_dp, 22.0_dp, 15.0_dp])
    y = convolve(a, b, mode='same')
    call check_real(y, [4.0_dp, 13.0_dp, 22.0_dp])
    y = convolve(a, b, mode='valid')
    call check_real(y, [13.0_dp, 22.0_dp])
    y = correlate(a, b)
    call check_real(y, [5.0_dp, 14.0_dp, 23.0_dp, 12.0_dp])
    y = correlate(a, b, mode='same')
    call check_real(y, [5.0_dp, 14.0_dp, 23.0_dp])
    y = correlate(a, b, mode='valid')
    call check_real(y, [14.0_dp, 23.0_dp])

    y = convolve(b, a, mode='same')
    call check_real(y, [13.0_dp, 22.0_dp])
    y = correlate(b, a, mode='valid')
    call check_real(y, [23.0_dp, 14.0_dp])

    lags = correlation_lags(3, 2, mode='same', status=ierr)
    if (ierr /= convolution_status_success) error stop 'lags status'
    call check_integer(lags, [0, 1, 2])
    lags = correlation_lags(4, 2, mode='same')
    call check_integer(lags, [-1, 0, 1, 2])
    lags = correlation_lags(3, 2, mode='full')
    call check_integer(lags, [-1, 0, 1, 2])
    lags = correlation_lags(2, 3, mode='valid')
    call check_integer(lags, [-1, 0])

    ca = [cmplx(1.0_dp, 2.0_dp, dp), cmplx(3.0_dp, -1.0_dp, dp)]
    cb = [cmplx(-2.0_dp, 1.0_dp, dp), cmplx(0.0_dp, 4.0_dp, dp)]
    z = convolve(ca, cb)
    call check_complex(z, [ca(1)*cb(1), ca(1)*cb(2)+ca(2)*cb(1), ca(2)*cb(2)])
    z = correlate(ca, cb)
    call check_complex(z, [ca(1)*conjg(cb(2)), &
        ca(1)*conjg(cb(1))+ca(2)*conjg(cb(2)), ca(2)*conjg(cb(1))])

    y = convolve(a, [real(dp) :: ], status=ierr)
    if (ierr /= convolution_status_invalid_input .or. size(y) /= 0) &
        error stop 'empty convolve status'
    y = correlate(a, b, mode='invalid', status=ierr)
    if (ierr /= convolution_status_invalid_input) error stop 'bad mode status'
    y = convolve(a, b, method='fft', status=ierr)
    if (ierr /= convolution_status_success) error stop 'fft status'
    call check_real(y, [4.0_dp,13.0_dp,22.0_dp,15.0_dp])
    lags = correlation_lags(0, 3, status=ierr)
    if (ierr /= convolution_status_invalid_input) error stop 'bad lag lengths'

    do n = 1, 20
        do m = 1, 20
            block
                real(dp), allocatable :: x(:), q(:), direct(:)
                allocate(x(n), q(m))
                do i = 1, n
                    x(i) = real(i, dp) / real(n, dp)
                end do
                do i = 1, m
                    q(i) = real(i, dp) / real(m, dp)
                end do
                direct = convolve(x, q)
                y = convolve(q, x)
                call check_real(y, direct)
                y = correlate(x, q)
                direct = convolve(x, q(m:1:-1))
                call check_real(y, direct)
                if (size(convolve(x, q, 'valid')) /= abs(n-m)+1) &
                    error stop 'valid output size'
            end block
        end do
    end do
    print '(a)', 'test_signal_convolve passed'

contains
    subroutine check_real(actual, expected)
        real(dp), intent(in) :: actual(:), expected(:)
        if (size(actual) /= size(expected)) error stop 'real size mismatch'
        if (any(abs(actual - expected) > 1.0e-12_dp)) then
            print *, actual, expected
            error stop 'real convolution mismatch'
        end if
    end subroutine

    subroutine check_complex(actual, expected)
        complex(dp), intent(in) :: actual(:), expected(:)
        if (size(actual) /= size(expected)) error stop 'complex size mismatch'
        if (any(abs(actual - expected) > 1.0e-12_dp)) error stop 'complex mismatch'
    end subroutine

    subroutine check_integer(actual, expected)
        integer, intent(in) :: actual(:), expected(:)
        if (size(actual) /= size(expected)) error stop 'lag size mismatch'
        if (any(actual /= expected)) error stop 'lag mismatch'
    end subroutine
end program test_signal_convolve
