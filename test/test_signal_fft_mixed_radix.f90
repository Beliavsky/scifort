! SPDX-License-Identifier: MIT
! Regression for recursive mixed-radix transforms and reusable plans.
program test_signal_fft_mixed_radix
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, forward_fft, spectral_status_success, &
        spectral_status_invalid_input
    implicit none
    integer, parameter :: lengths(*) = [18,21,25,27,28,35,42,49,63,75,81,98,105,120, &
        125,147,175,210,225,256,343,375,420,480,625,729,960,1000,1024,1225, &
        1728,2048,2401,3072,4096,5040,8192,16807]
    type(fft_plan) :: plan
    complex(dp), allocatable :: x(:), y(:), z(:), strided(:), baseline(:)
    real(dp) :: pi, tol
    integer :: n, i, j, status, k, cases

    pi = acos(-1.0_dp)
    cases = 0
    do i = 1, size(lengths)
        n = lengths(i)
        x = [(cmplx(sin(0.031_dp*real(j,dp))+0.6_dp*cos(0.022_dp*real(j,dp)), &
              cos(0.071_dp*real(j,dp))-0.2_dp*sin(0.063_dp*real(j,dp)), dp),j=1,n)]
        tol = 4.0e-12_dp*max(1.0_dp,real(n,dp))
        call plan%initialize(n,status)
        if (status /= spectral_status_success) error stop 'init mixed'
        if (plan%length() /= n) error stop 'plan length mixed'
        do k = 1, 3
            baseline = x*cmplx(1.0_dp+0.1_dp*real(k,dp), 0.11_dp*real(k,dp),dp)
            y = baseline
            z = baseline
            call plan%execute(y,status)
            if (status /= spectral_status_success) error stop 'execute mixed'
            call forward_fft(z)
            if (maxval(abs(y-z)) > tol) error stop 'planned vs direct mixed'
            call plan%inverse(y,status)
            if (status /= spectral_status_success) error stop 'inverse status mixed'
            if (maxval(abs(y-baseline)) > tol) error stop 'inverse data mixed'
            cases = cases + 1
        end do
        ! Verify forward DFT at a few bins against an independent direct sum.
        do k = 0, min(n-1,7)
            baseline = x
            call plan%execute(baseline,status)
            z = cmplx(0.0_dp,0.0_dp,dp)
            do j = 0, n-1
                z(1) = z(1) + x(j+1)*cmplx(cos(-2*pi*real(j,dp)*real(k,dp)/real(n,dp)), &
                    sin(-2*pi*real(j,dp)*real(k,dp)/real(n,dp)),dp)
            end do
            if (abs(z(1)-baseline(k+1)) > tol) error stop 'direct bin mixed'
        end do
        ! A strided view must preserve the untouched gaps.
        allocate(strided(2*n))
        strided = (-777.0_dp,4.0_dp)
        strided(1:2*n:2) = x
        y = x
        call plan%execute(y,status)
        call plan%execute(strided(1:2*n:2),status)
        if (status /= spectral_status_success) error stop 'strided mixed'
        if (maxval(abs(strided(1:2*n:2)-y)) > tol) error stop 'strided value mixed'
        if (any(strided(2:2*n:2) /= cmplx(-777.0_dp,4.0_dp,dp))) error stop 'strided gaps mixed'
        deallocate(strided)
        ! Invalid plan reinitialization and mismatched input are transactional.
        call plan%initialize(-1,status)
        if (status /= spectral_status_invalid_input .or. plan%length() /= n) error stop 'bad init mixed'
        baseline = [(cmplx(real(j,dp),-real(j,dp),dp),j=1,3)]
        if (n /= 3) then
            z = baseline
            call plan%execute(z,status)
            if (status /= spectral_status_invalid_input .or. any(z /= baseline)) &
                error stop 'invalid execute mixed'
        end if
    end do
    ! Reinitialize a mixed-radix plan for a Bluestein length and back again.
    call plan%initialize(1009,status)
    if (status /= spectral_status_success) error stop 'switch bluestein'
    call plan%initialize(210,status)
    if (status /= spectral_status_success) error stop 'switch mixed'
    call plan%clear()
    if (plan%length() /= 0) error stop 'clear mixed'
    print '(a,i0,a,i0,a)', 'test_signal_fft_mixed_radix passed (', cases, &
        ' repeated transforms, ', size(lengths),' lengths)'
end program test_signal_fft_mixed_radix
