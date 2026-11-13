! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_special_complex_zeta
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, &
        ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_special, only : zeta, loggamma
    implicit none
    complex(dp) :: points(5), values(5), s, w, factor
    real(dp) :: nan_value, inf, x, y
    integer :: i, checks
    checks = 0
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    inf = ieee_value(0.0_dp,ieee_positive_inf)

    call close(zeta(cmplx(2.0_dp,3.0_dp,dp)), &
        cmplx(0.79802198514627572062_dp,-0.11374430805293850022_dp,dp),3e-14_dp,'right')
    call close(zeta(cmplx(-2.0_dp,1.0_dp,dp)), &
        cmplx(0.02919559324346427876_dp,-0.01778710697084208175_dp,dp),3e-14_dp,'left')
    call close(zeta(cmplx(0.5_dp,14.0_dp,dp)), &
        cmplx(0.02224114260999358925_dp,-0.10325812326645005790_dp,dp),3e-14_dp,'critical')
    call close(zeta(cmplx(0.5_dp,14.134725141734695_dp,dp)), &
        cmplx(0.0_dp,0.0_dp,dp),3e-14_dp,'first nontrivial zero')
    call close(zeta(cmplx(-20.0_dp,10.0_dp,dp)), &
        cmplx(-9609819.6252705349031_dp,24850620.200734894776_dp,dp),2e-13_dp,'reflection')
    call close(zeta(cmplx(0.5_dp,1000.0_dp,dp)), &
        cmplx(0.35633436719439605507_dp,0.93199783123299366512_dp,dp),3e-12_dp,'high')

    points = [cmplx(-5.0_dp,2.0_dp,dp),cmplx(-0.25_dp,3.0_dp,dp), &
        cmplx(0.5_dp,100.0_dp,dp),cmplx(2.0_dp,-10.0_dp,dp),cmplx(40.0_dp,1.0_dp,dp)]
    values = zeta(points)
    do i = 1, size(points)
        call close(values(i),zeta(points(i)),0.0_dp,'elemental scalar agreement')
        call close(zeta(conjg(points(i))),conjg(values(i)),2e-13_dp,'conjugate symmetry')
    end do
    do i = -20, 20
        x = real(i,dp)/4.0_dp
        if (x == 1.0_dp) cycle
        call close(zeta(cmplx(x,0.0_dp,dp)),cmplx(zeta(x),0.0_dp,dp),0.0_dp,'real axis')
    end do
    call close(zeta(cmplx(0.0_dp,0.0_dp,dp)),cmplx(-0.5_dp,0.0_dp,dp),0.0_dp,'origin')
    do i = 1, 10
        call close(zeta(cmplx(-real(2*i,dp),0.0_dp,dp)), &
            cmplx(0.0_dp,0.0_dp,dp),0.0_dp,'trivial zero')
    end do
    w = zeta(cmplx(1.0_dp,0.0_dp,dp))
    call require(.not. ieee_is_finite(real(w,dp)) .and. .not. ieee_is_nan(real(w,dp)), &
        'real-axis pole infinity')
    call require(aimag(w) == 0.0_dp,'pole imaginary zero')
    do i = 1, 3
        y = 10.0_dp**(-100*i)
        s = cmplx(1.0_dp,y,dp)
        call close((s-1.0_dp)*zeta(s),cmplx(1.0_dp,0.0_dp,dp),3e-14_dp,'pole residue')
    end do
    call close(zeta(cmplx(0.0_dp,1e-100_dp,dp)), &
        cmplx(-0.5_dp,-9.1893853320467274e-101_dp,dp),3e-14_dp,'near origin')
    w = zeta(cmplx(-2.0_dp,1e-100_dp,dp))
    call require(abs(real(w,dp)) < 1e-114_dp,'near trivial zero real part')
    call require(abs(aimag(w)/(-3.0448457058393271e-102_dp)-1.0_dp) < 3e-13_dp, &
        'near trivial zero derivative')
    s = cmplx(-3.2_dp,4.5_dp,dp)
    factor = 2.0_dp*exp((s-1.0_dp)*log(2.0_dp*acos(-1.0_dp))) &
        *sin(0.5_dp*acos(-1.0_dp)*s)*exp(loggamma(1.0_dp-s))
    call close(zeta(s),factor*zeta(1.0_dp-s),3e-13_dp,'functional equation')
    call nan_case(cmplx(nan_value,1.0_dp,dp),'NaN real')
    call nan_case(cmplx(1.0_dp,nan_value,dp),'NaN imaginary')
    call nan_case(cmplx(inf,1.0_dp,dp),'infinite real with nonzero imaginary')
    call nan_case(cmplx(1.0_dp,inf,dp),'infinite imaginary')
    call nan_case(cmplx(10001.0_dp,1.0_dp,dp),'real range limit')
    call nan_case(cmplx(0.5_dp,100001.0_dp,dp),'imaginary range limit')
    call nan_case(cmplx(-1000.0_dp,1.0_dp,dp),'unrepresentable reflected value')
    call nan_case(cmplx(1.0_dp,1e-310_dp,dp),'pole overflow guard')
    call close(zeta(cmplx(inf,0.0_dp,dp)),cmplx(1.0_dp,0.0_dp,dp),0.0_dp,'axis infinity')
    call require(zeta(2.0_dp,1.0_dp) == zeta(2.0_dp),'real Hurwitz API unchanged')
    print '(a,i0,a)', 'Complex zeta: ',checks,' checks passed.'

contains
    subroutine close(got, ref, tol, label)
        implicit none
        complex(dp), intent(in) :: got !! Computed result.
        complex(dp), intent(in) :: ref !! Independent reference or identity result.
        real(dp), intent(in) :: tol !! Bound on abs(error)/(1+abs(reference)).
        character(len=*), intent(in) :: label !! Assertion description.
        call require(ieee_is_finite(real(got,dp)) .and. ieee_is_finite(aimag(got)),label)
        call require(abs(got-ref) <= tol*(1.0_dp+abs(ref)),label)
    end subroutine close

    subroutine nan_case(s, label)
        implicit none
        complex(dp), intent(in) :: s !! Input expected to yield NaNs.
        character(len=*), intent(in) :: label !! Assertion description.
        complex(dp) :: result
        result = zeta(s)
        call require(ieee_is_nan(real(result,dp)) .and. ieee_is_nan(aimag(result)),label)
    end subroutine nan_case

    subroutine require(ok, label)
        implicit none
        logical, intent(in) :: ok !! Assertion that must hold.
        character(len=*), intent(in) :: label !! Assertion description.
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',label
            error stop 1
        end if
        checks = checks+1
    end subroutine require
end program test_special_complex_zeta
