! SPDX-License-Identifier: MIT
! Native regression tests for classical analog lowpass prototypes.
program test_signal_analog_prototypes
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_analog_prototypes, only : analog_prototype_result, &
        buttap, cheb1ap, cheb2ap, besselap, prototype_success, &
        prototype_invalid_input, prototype_too_expensive
    implicit none
    type(analog_prototype_result) :: prototype
    complex(dp) :: h, p
    real(dp) :: rp, rs, freq, gain_dc
    integer :: n, i, j

    do n = 0, 32
        prototype = buttap(n)
        call validate(prototype, n, 0)
        if (abs(prototype%gain-1.0_dp) > 1.0e-13_dp) error stop 'Butterworth gain'
        do i = 1, size(prototype%poles)
            if (abs(abs(prototype%poles(i))-1.0_dp) > 2.0e-14_dp) error stop 'Butterworth circle'
        end do
    end do

    do n = 0, 16
        do j = 1, 4
            select case (j)
            case (1)
                rp = 0.01_dp
                rs = 0.5_dp
            case (2)
                rp = 0.5_dp
                rs = 15.0_dp
            case (3)
                rp = 3.0_dp
                rs = 40.0_dp
            case default
                rp = 40.0_dp
                rs = 80.0_dp
            end select
            prototype = cheb1ap(n,rp)
            call validate(prototype,n,0)
            gain_dc = 1.0_dp
            if (mod(n,2) == 0) gain_dc = 10.0_dp**(-rp/20.0_dp)
            if (n == 0) then
                h = cmplx(prototype%gain,0.0_dp,dp)
            else
                h = response(prototype,cmplx(0.0_dp,0.0_dp,dp))
            end if
            if (abs(real(h,dp)-gain_dc) > 2.0e-10_dp) error stop 'Chebyshev I DC gain'
            prototype = cheb2ap(n,rs)
            call validate(prototype,n,n-mod(n,2))
            h = response(prototype,cmplx(0.0_dp,0.0_dp,dp))
            if (abs(h-1.0_dp) > 1.0e-10_dp) error stop 'Chebyshev II DC gain'
            do i = 1, size(prototype%zeros)
                if (abs(real(prototype%zeros(i),dp)) > 1.0e-13_dp) error stop 'Chebyshev II zeros'
            end do
        end do
    end do

    do n = 0, 20
        do j = 1, 3
            select case(j)
            case (1)
                prototype = besselap(n,'phase')
            case (2)
                prototype = besselap(n,'delay')
            case default
                prototype = besselap(n,'mag')
            end select
            call validate(prototype,n,0)
            if (n == 0) cycle
            h = response(prototype,cmplx(0.0_dp,0.0_dp,dp))
            if (abs(h-1.0_dp) > 5.0e-9_dp) error stop 'Bessel DC gain'
            if (j == 1) then
                if (abs(prototype%gain-1.0_dp) > 1.0e-12_dp) error stop 'Bessel phase normalization'
            else if (j == 2) then
                freq = -sum(real(1.0_dp/prototype%poles,dp))
                if (abs(freq-1.0_dp) > 1.0e-8_dp) error stop 'Bessel DC group delay'
            else
                h = response(prototype,cmplx(0.0_dp,1.0_dp,dp))
                if (abs(abs(h)-1.0_dp/sqrt(2.0_dp)) > 5.0e-10_dp) error stop 'Bessel half power'
            end if
            ! All conjugate pairs must yield real coefficients.
            do i = 1, n
                p = conjg(prototype%poles(i))
                if (minval(abs(prototype%poles-p)) > 5.0e-11_dp) error stop 'Bessel conjugate symmetry'
            end do
        end do
    end do

    prototype = buttap(-1)
    if (prototype%status /= prototype_invalid_input) error stop 'negative Butterworth order'
    prototype = buttap(65)
    if (prototype%status /= prototype_invalid_input) error stop 'oversized Butterworth order'
    prototype = cheb1ap(2,0.0_dp)
    if (prototype%status /= prototype_invalid_input) error stop 'Chebyshev I zero ripple'
    prototype = cheb1ap(2,ieee_value(0.0_dp,ieee_quiet_nan))
    if (prototype%status /= prototype_invalid_input) error stop 'Chebyshev I NaN ripple'
    prototype = cheb2ap(3,-4.0_dp)
    if (prototype%status /= prototype_invalid_input) error stop 'negative stopband attenuation'
    prototype = cheb2ap(-3,20.0_dp)
    if (prototype%status /= prototype_invalid_input) error stop 'negative Chebyshev II order'
    prototype = besselap(-2)
    if (prototype%status /= prototype_invalid_input) error stop 'negative Bessel order'
    prototype = besselap(21)
    if (prototype%status /= prototype_too_expensive) error stop 'Bessel order limit'
    prototype = besselap(3,'unknown')
    if (prototype%status /= prototype_invalid_input) error stop 'Bessel unknown normalization'

    print '(a)', 'test_signal_analog_prototypes passed'
contains
    subroutine validate(out, expected_poles, expected_zeros)
        type(analog_prototype_result), intent(in) :: out !! prototype being tested
        integer, intent(in) :: expected_poles !! number of poles
        integer, intent(in) :: expected_zeros !! number of finite zeros
        if (out%status /= prototype_success) error stop 'prototype status'
        if (.not. allocated(out%zeros) .or. .not. allocated(out%poles)) error stop 'allocatable roots'
        if (size(out%zeros) /= expected_zeros) error stop 'prototype zero count'
        if (size(out%poles) /= expected_poles) error stop 'prototype pole count'
        if (any(real(out%poles,dp) >= 0.0_dp)) error stop 'unstable analog prototype'
        if (out%gain <= 0.0_dp) error stop 'nonpositive prototype gain'
    end subroutine validate

    pure function response(out, s) result(value)
        type(analog_prototype_result), intent(in) :: out !! lowpass prototype
        complex(dp), intent(in) :: s !! complex continuous-time frequency
        complex(dp) :: value
        integer :: idx
        value = cmplx(out%gain,0.0_dp,dp)
        do idx = 1, size(out%zeros)
            value = value*(s-out%zeros(idx))
        end do
        do idx = 1, size(out%poles)
            value = value/(s-out%poles(idx))
        end do
    end function response
end program test_signal_analog_prototypes
