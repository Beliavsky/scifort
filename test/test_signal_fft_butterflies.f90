! SPDX-License-Identifier: MIT
! Regression for the pair-sum radix-3, radix-5, and radix-7 butterflies.
program test_signal_fft_butterflies
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, forward_fft, spectral_status_success
    implicit none
    integer, parameter :: lengths(*) = [21,25,27,35,49,63,75,105,125,147,175, &
        243,343,375,625,729,875,1225,1701,2401,3125,5040]
    type(fft_plan) :: plan
    complex(dp), allocatable :: x(:), y(:), z(:), strided(:)
    complex(dp) :: total, root
    real(dp) :: pi, tol, theta
    integer :: i, j, k, n, m, status, cases, mode, bins

    pi = acos(-1.0_dp)
    cases = 0
    do i = 1, size(lengths)
        n = lengths(i)
        allocate(x(n), y(n), z(n), strided(2*n))
        call plan%initialize(n,status)
        if (status /= spectral_status_success) error stop 'butterfly plan init'
        do mode = 1, 3
            select case (mode)
            case (1)
                do j = 1, n
                    x(j) = cmplx(sin(0.023_dp*real(j*j,dp))+cos(0.31_dp*real(j,dp)), &
                        cos(0.017_dp*real(j*j,dp))-sin(0.37_dp*real(j,dp)),dp)
                end do
            case (2)
                x = (0.0_dp,0.0_dp)
                x(1+mod(n/3,n)) = (1.25_dp,-0.75_dp)
            case (3)
                do j = 1, n
                    theta = 2.0_dp*pi*real(11*(j-1),dp)/real(n,dp)
                    x(j) = cmplx(cos(theta),sin(theta),dp)
                end do
            end select
            y = x
            z = x
            call plan%execute(y,status)
            if (status /= spectral_status_success) error stop 'butterfly execute'
            call forward_fft(z)
            tol = 4.0e-11_dp*max(1.0_dp,maxval(abs(y)))
            if (maxval(abs(y-z)) > tol) error stop 'planned/unplanned butterflies'
            if (mode == 1) then
                ! Compare multiple bins against a direct DFT, independently
                ! of the new Cooley-Tukey combination formula.
                bins = n
                if (n > 175) bins = 24
                do k = 0, bins-1
                    total = (0.0_dp,0.0_dp)
                    do j = 0, n-1
                        theta = -2.0_dp*pi*real(j,dp)*real(k,dp)/real(n,dp)
                        root = cmplx(cos(theta),sin(theta),dp)
                        total = total + x(j+1)*root
                    end do
                    if (abs(y(k+1)-total) > 6.0e-11_dp*real(n,dp)) &
                        error stop 'independent DFT butterfly'
                end do
            end if
            ! Noncontiguous array sections must not overwrite intervening data.
            strided = (-999.0_dp,999.0_dp)
            strided(2:2*n:2) = x
            call plan%execute(strided(2:2*n:2),status)
            if (status /= spectral_status_success) error stop 'strided butterfly'
            if (maxval(abs(strided(2:2*n:2)-y)) > tol) error stop 'strided outputs butterfly'
            if (any(strided(1:2*n:2) /= cmplx(-999.0_dp,999.0_dp,dp))) &
                error stop 'strided gaps butterfly'
            call plan%inverse(y,status)
            if (status /= spectral_status_success) error stop 'inverse butterfly'
            if (maxval(abs(y-x)) > 6.0e-11_dp*max(1.0_dp,maxval(abs(x)))) &
                error stop 'roundtrip butterfly'
            cases = cases + 1
        end do
        deallocate(x, y, z, strided)
    end do
    call plan%clear()
    if (plan%length() /= 0) error stop 'clear butterfly'
    print '(a,i0,a)', 'test_signal_fft_butterflies passed (', cases, ' cases)'
end program test_signal_fft_butterflies
