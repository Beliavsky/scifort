! SPDX-License-Identifier: MIT
! Regression for the non-power-of-two shared Fourier transform.
program test_signal_fft_bluestein
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft
    implicit none
    integer, parameter :: lengths(*) = [1,2,3,5,7,11,15,16,17,19,31,33,35,45, &
        63,64,65,97,101,127,129,257,509,1009]
    complex(dp), allocatable :: x(:), y(:), reference(:), work(:)
    real(dp) :: pi, bound
    integer :: i, j, n

    pi = acos(-1.0_dp)
    allocate(x(0))
    call forward_fft(x)
    if (size(x) /= 0) error stop 'empty FFT'

    do i = 1, size(lengths)
        n = lengths(i)
        x = [(cmplx(sin(0.13_dp*real(j,dp)) + cos(0.071_dp*real(j*j,dp)), &
             cos(0.17_dp*real(j,dp)) - sin(0.041_dp*real(j*j,dp)),dp),j=1,n)]
        y = x
        reference = direct_transform(x)
        call forward_fft(y)
        bound = 2.0e-11_dp*max(1.0_dp,maxval(abs(reference)))
        if (maxval(abs(y-reference)) > bound) then
            print *, 'length, difference, bound:',n,maxval(abs(y-reference)),bound
            error stop 'FFT forward mismatch'
        end if
        work = conjg(y)
        call forward_fft(work)
        work = conjg(work)/real(n,dp)
        if (maxval(abs(work-x)) > 1.0e-11_dp*max(1.0_dp,maxval(abs(x)))) then
            error stop 'FFT inverse via conjugation'
        end if
        if (abs(sum(abs(x)**2)-sum(abs(y)**2)/real(n,dp)) > &
             1.0e-9_dp*max(1.0_dp,sum(abs(x)**2))) error stop 'FFT Parseval'
        ! Ensure strided input/output works as for any Fortran assumed-shape array.
        deallocate(work)
        allocate(work(2*n))
        work = (-500.0_dp, 3.0_dp)
        work(1:2*n:2) = x
        call forward_fft(work(1:2*n:2))
        if (maxval(abs(work(1:2*n:2)-y)) > bound) error stop 'strided FFT'
        if (any(work(2:2*n:2) /= cmplx(-500.0_dp,3.0_dp,dp))) error stop 'FFT modified gap'
        deallocate(work)
    end do
    print '(a,i0,a)', 'test_signal_fft_bluestein passed (',size(lengths),' lengths)'
contains
    function direct_transform(input) result(output)
        complex(dp), intent(in) :: input(:) !! source samples
        complex(dp) :: output(size(input))
        complex(dp) :: factor
        integer :: j, k, n
        n = size(input)
        do k = 0, n-1
            output(k+1) = (0.0_dp,0.0_dp)
            do j = 0, n-1
                factor = cmplx(cos(2.0_dp*pi*real(j,dp)*real(k,dp)/real(n,dp)), &
                    -sin(2.0_dp*pi*real(j,dp)*real(k,dp)/real(n,dp)),dp)
                output(k+1) = output(k+1) + input(j+1)*factor
            end do
        end do
    end function direct_transform
end program test_signal_fft_bluestein
