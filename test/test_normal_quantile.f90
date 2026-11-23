! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! References: mpmath 1.3.0, 400 digits, sqrt(2)*erfinv(2*p-1), binary64 p.
program test_normal_quantile
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, &
        ieee_quiet_nan, ieee_positive_inf, ieee_invalid, ieee_get_halting_mode, &
        ieee_set_halting_mode, ieee_set_flag
    use scifort_kinds, only : dp
    use scifort_normal, only : normal_ppf, normal_isf, normal_cdf, normal_sf, normal_logcdf
    implicit none
    real(dp), parameter :: probabilities(32) = [ &
        4.94065645841246544e-324_dp, &
        9.99988867182683005e-321_dp, &
        9.99999999999996945e-311_dp, &
        2.22507385850720138e-308_dp, &
        1.00000000000000003e-300_dp, &
        9.99999999999999982e-201_dp, &
        1.00000000000000005e-101_dp, &
        1.00000000000000002e-100_dp, &
        1.00000000000000001e-50_dp, &
        1.26641655490941740e-14_dp, &
        1.26641655490941755e-14_dp, &
        1.26641655490941771e-14_dp, &
        1.00000000000000004e-10_dp, &
        9.99999999999999955e-07_dp, &
        1.00000000000000002e-03_dp, &
        1.00000000000000002e-02_dp, &
        1.00000000000000006e-01_dp, &
        1.35335283236612675e-01_dp, &
        1.35335283236612702e-01_dp, &
        1.35335283236612730e-01_dp, &
        2.50000000000000000e-01_dp, &
        4.99999999999999944e-01_dp, &
        5.00000000000000000e-01_dp, &
        5.00000000000000111e-01_dp, &
        7.50000000000000000e-01_dp, &
        8.64664716763387187e-01_dp, &
        8.64664716763387298e-01_dp, &
        8.64664716763387409e-01_dp, &
        9.00000000000000022e-01_dp, &
        9.74999999999999978e-01_dp, &
        9.99998999999999971e-01_dp, &
        9.99999999999999889e-01_dp ]
    real(dp), parameter :: reference(32) = [ &
        -3.84674056171443439e+01_dp, &
        -3.82691253430326483e+01_dp, &
        -3.76630603319495236e+01_dp, &
        -3.75193793471445005e+01_dp, &
        -3.70470962993612005e+01_dp, &
        -3.02055941795796414e+01_dp, &
        -2.13811818695259106e+01_dp, &
        -2.12734535609653257e+01_dp, &
        -1.49333375347884889e+01_dp, &
        -7.62019982525522543e+00_dp, &
        -7.62019982525522543e+00_dp, &
        -7.62019982525522543e+00_dp, &
        -6.36134090240405659e+00_dp, &
        -4.75342430882289868e+00_dp, &
        -3.09023230616781364e+00_dp, &
        -2.32634787404084120e+00_dp, &
        -1.28155156554460037e+00_dp, &
        -1.10151962849875029e+00_dp, &
        -1.10151962849875029e+00_dp, &
        -1.10151962849875007e+00_dp, &
        -6.74489750196081705e-01_dp, &
        -1.39145821233588356e-16_dp, &
        0.00000000000000000e+00_dp, &
        2.78291642467176713e-16_dp, &
        6.74489750196081705e-01_dp, &
        1.10151962849874963e+00_dp, &
        1.10151962849875029e+00_dp, &
        1.10151962849875074e+00_dp, &
        1.28155156554460059e+00_dp, &
        1.95996398454005383e+00_dp, &
        4.75342430881708733e+00_dp, &
        8.20953615160138739e+00_dp ]
    real(dp) :: x, previous, p, nan, inf, result, z(10001), grid(10001)
    logical :: halt_invalid
    integer :: i, checks
    checks=0
    do i=1,size(probabilities)
        p=probabilities(i)
        x=normal_ppf(p)
        call check(ieee_is_finite(x),'finite quantile')
        if (reference(i) == 0.0_dp) then
            call check(x == 0.0_dp,'exact median')
        else
            call check(abs(x-reference(i)) <= 8*epsilon(x)*abs(reference(i)), &
                'independent high-precision quantile')
        end if
        call check(normal_isf(p) == -x,'ISF reflection uses the same small probability')
        call check(normal_ppf(p,2.0_dp,3.0_dp) == 2.0_dp+3.0_dp*x,'location and scale')
        if (p < 1e-100_dp) then
            call check(abs(normal_logcdf(x)-log(p)) < 5e-13_dp*abs(log(p)), &
                'extreme log-CDF round trip')
        else if (p <= 0.5_dp) then
            call check(abs(normal_cdf(x)/p-1.0_dp) < 1e-12_dp,'CDF round trip')
            call check(abs(normal_sf(normal_isf(p))/p-1.0_dp) < 1e-12_dp,'SF round trip')
        end if
    end do
    grid=[(real(i,dp)/10002.0_dp,i=1,10001)]
    z=normal_ppf(grid)
    previous=-huge(1.0_dp)
    do i=1,size(grid)
        call check(z(i) > previous,'dense-grid monotonicity')
        call check(z(i) == normal_ppf(grid(i)),'scalar and elemental array agreement')
        previous=z(i)
    end do
    call check(normal_ppf(0.0_dp) < -huge(1.0_dp),'PPF zero endpoint')
    call check(normal_ppf(1.0_dp) > huge(1.0_dp),'PPF one endpoint')
    call check(normal_isf(0.0_dp) > huge(1.0_dp),'ISF zero endpoint')
    call check(normal_isf(1.0_dp) < -huge(1.0_dp),'ISF one endpoint')
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    ! Invalid inputs return NaN, including when ordinary arithmetic traps are enabled.
    call ieee_get_halting_mode(ieee_invalid,halt_invalid)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call check(ieee_is_nan(normal_ppf(nan)),'NaN probability')
    call check(ieee_is_nan(normal_isf(nan)),'NaN ISF probability')
    call check(ieee_is_nan(normal_ppf(inf)),'infinite probability')
    call check(ieee_is_nan(normal_ppf(-inf)),'negative infinite probability')
    call check(ieee_is_nan(normal_ppf(-epsilon(1.0_dp))),'negative probability')
    call check(ieee_is_nan(normal_ppf(1.0_dp+epsilon(1.0_dp))),'probability above one')
    call check(ieee_is_nan(normal_ppf(0.5_dp,scale=0.0_dp)),'zero scale')
    call check(ieee_is_nan(normal_isf(0.5_dp,scale=-1.0_dp)),'negative scale')
    call check(ieee_is_nan(normal_ppf(0.5_dp,loc=nan)),'NaN location')
    call check(ieee_is_nan(normal_ppf(0.5_dp,scale=inf)),'infinite scale')
    result=normal_ppf(0.5_dp,scale=nan)
    call check(ieee_is_nan(result),'NaN scale')
    call ieee_set_flag(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_invalid,halt_invalid)
    print '(a,i0,a)','Normal quantile: ',checks,' checks passed.'
contains
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! whether the assertion passed
        character(*), intent(in) :: label !! assertion context
        if (.not. ok) then
            print '(a,a)','FAILED: ',label
            error stop 1
        end if
        checks=checks+1
    end subroutine check
end program test_normal_quantile
