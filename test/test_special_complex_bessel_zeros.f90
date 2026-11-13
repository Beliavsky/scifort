! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_special_complex_bessel_zeros
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_special, only : y0_zeros, y1_zeros, y1p_zeros, yv, yvp, &
        bessel_zeros_status_success, bessel_zeros_status_invalid_argument
    implicit none
    complex(dp), allocatable :: roots(:), values(:), prefix(:), pv(:)
    complex(dp) :: f, expected
    integer :: family, k, status, checks
    real(dp) :: amplitude
    checks = 0
    do family = 0, 2
        select case (family)
        case (0)
            call y0_zeros(16,roots,values,status)
            expected = cmplx(-2.40301663199377604042_dp,0.53988231304266412662_dp,dp)
            call y0_zeros(1,prefix,pv)
        case (1)
            call y1_zeros(16,roots,values,status)
            expected = cmplx(-0.50274327295272190335_dp,0.78624371447749664787_dp,dp)
            call y1_zeros(1,prefix,pv)
        case default
            call y1p_zeros(16,roots,values,status)
            expected = cmplx(0.57678512860684465965_dp,0.90398479216641068853_dp,dp)
            call y1p_zeros(1,prefix,pv)
        end select
        call require(status == bessel_zeros_status_success,'complex status')
        call require(size(roots) == 16 .and. size(values) == 16,'complex sizes')
        call require(abs(roots(1)-expected) < 3e-14_dp,'multiprecision first root')
        call require(size(prefix) == 1 .and. size(pv) == 1,'optional status')
        call require(prefix(1) == roots(1) .and. pv(1) == values(1),'prefix agreement')
        call require(all(aimag(roots) > 0.0_dp),'upper half-plane')
        call require(all(real(roots(2:),dp) < real(roots(:15),dp)),'ordered unique roots')
        if (family == 2) then
            call require(real(roots(1),dp) > 0.0_dp,'exceptional positive real root')
            call require(all(real(roots(2:),dp) < 0.0_dp),'negative real string')
        else
            call require(all(real(roots,dp) < 0.0_dp),'negative real roots')
        end if
        do k = 1, size(roots)
            call require(finite(roots(k)) .and. finite(values(k)),'finite outputs')
            select case (family)
            case (0)
                f = yv(0.0_dp,roots(k))
                expected = -yv(1.0_dp,roots(k))
            case (1)
                f = yv(1.0_dp,roots(k))
                expected = yvp(1.0_dp,roots(k))
            case default
                f = yvp(1.0_dp,roots(k))
                expected = yv(1.0_dp,roots(k))
            end select
            amplitude = max(abs(values(k)),1.0_dp/sqrt(abs(roots(k))))
            call require(abs(f) <= 2e-13_dp*max(1.0_dp,abs(roots(k)))*amplitude,'residual')
            call require(abs(values(k)-expected) <= 3e-12_dp*amplitude,'companion value')
            if (family == 0) then
                expected = yvp(0.0_dp,roots(k))
                call require(abs(values(k)-expected) <= 3e-12_dp*amplitude,'Y0 derivative sign')
            end if
        end do
    end do
    do k = -1, 0
        call y0_zeros(k,roots,values,status)
        call invalid_result()
        call y1_zeros(k,roots,values,status)
        call invalid_result()
        call y1p_zeros(k,roots,values,status)
        call invalid_result()
    end do
    call y0_zeros(1001,roots,values,status)
    call invalid_result()
    call y1_zeros(huge(0),roots,values,status)
    call invalid_result()
    call y1p_zeros(1001,roots,values,status)
    call invalid_result()
    print '(a,i0,a)', 'Complex Bessel zeros: ',checks,' checks passed.'

contains
    pure logical function finite(z) result(ok)
        implicit none
        complex(dp), intent(in) :: z !! Value to classify.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite

    subroutine invalid_result()
        implicit none
        call require(status == bessel_zeros_status_invalid_argument,'invalid count status')
        call require(size(roots) == 0 .and. size(values) == 0,'empty invalid outputs')
    end subroutine invalid_result

    subroutine require(ok, label)
        implicit none
        logical, intent(in) :: ok !! Assertion that must hold.
        character(len=*), intent(in) :: label !! Failure description.
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',label
            error stop 1
        end if
        checks = checks+1
    end subroutine require
end program test_special_complex_bessel_zeros
