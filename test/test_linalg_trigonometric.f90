! SPDX-License-Identifier: MIT
program test_linalg_trigonometric
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_trigonometric, only : sinm, cosm, tanm, sinhm, coshm, tanhm
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: threshold = 2.0e-12_dp
    real(dp) :: a(2, 2), x(2, 2), s(2, 2), c(2, 2), t(2, 2), eye(2, 2)
    complex(dp) :: ac(2, 2), xc(2, 2), sc(2, 2), cc(2, 2), tc(2, 2)
    real(dp) :: tiny_a(2, 2)
    integer :: status

    eye = 0.0_dp
    eye(1, 1) = 1.0_dp
    eye(2, 2) = 1.0_dp
    a = 0.0_dp
    call sinm(a, s, status)
    call verify(status == linalg_status_success, 'zero sinm status')
    call verify(maxval(abs(s)) == 0.0_dp, 'zero sinm value')
    call cosm(a, c, status)
    call verify(status == linalg_status_success .and. all(c == eye), 'zero cosm')
    call sinhm(a, s, status)
    call verify(status == linalg_status_success .and. all(s == 0.0_dp), 'zero sinhm')
    call coshm(a, c, status)
    call verify(status == linalg_status_success .and. all(c == eye), 'zero coshm')
    call tanm(a, t, status)
    call verify(status == linalg_status_success .and. all(t == 0.0_dp), 'zero tanm')
    call tanhm(a, t, status)
    call verify(status == linalg_status_success .and. all(t == 0.0_dp), 'zero tanhm')

    a(1, 2) = 1.5_dp
    call sinm(a, s, status)
    call verify(status == linalg_status_success .and. maxval(abs(s-a)) < threshold, 'nilpotent sinm')
    call cosm(a, c, status)
    call verify(status == linalg_status_success .and. maxval(abs(c-eye)) < threshold, 'nilpotent cosm')
    call tanm(a, t, status)
    call verify(status == linalg_status_success .and. maxval(abs(t-a)) < threshold, 'nilpotent tanm')
    call sinhm(a, s, status)
    call verify(status == linalg_status_success .and. maxval(abs(s-a)) < threshold, 'nilpotent sinhm')
    call coshm(a, c, status)
    call verify(status == linalg_status_success .and. maxval(abs(c-eye)) < threshold, 'nilpotent coshm')
    call tanhm(a, t, status)
    call verify(status == linalg_status_success .and. maxval(abs(t-a)) < threshold, 'nilpotent tanhm')

    a = 0.0_dp
    a(1, 1) = 0.3_dp
    a(2, 2) = -0.8_dp
    a(1, 2) = 0.4_dp
    call sinm(a, s, status)
    call verify(status == linalg_status_success, 'real sinm')
    call cosm(a, c, status)
    call verify(status == linalg_status_success, 'real cosm')
    call verify(maxval(abs(matmul(s, s)+matmul(c, c)-eye)) < threshold, 'sin^2+cos^2')
    call tanm(a, t, status)
    call verify(status == linalg_status_success, 'real tanm')
    call verify(maxval(abs(matmul(t, c)-s)) < threshold, 'tan*cos=sin')
    call sinhm(a, s, status)
    call verify(status == linalg_status_success, 'real sinhm')
    call coshm(a, c, status)
    call verify(status == linalg_status_success, 'real coshm')
    call verify(maxval(abs(matmul(c, c)-matmul(s, s)-eye)) < threshold, 'cosh^2-sinh^2')
    call tanhm(a, t, status)
    call verify(status == linalg_status_success, 'real tanhm')
    call verify(maxval(abs(matmul(t, c)-s)) < threshold, 'tanh*cosh=sinh')

    ac(1, 1) = cmplx(0.2_dp, 0.4_dp, dp)
    ac(1, 2) = cmplx(-0.7_dp, 0.3_dp, dp)
    ac(2, 1) = cmplx(0.6_dp, -0.1_dp, dp)
    ac(2, 2) = cmplx(-0.3_dp, 0.1_dp, dp)
    call sinm(ac, sc, status)
    call verify(status == linalg_status_success, 'complex sinm')
    call cosm(ac, cc, status)
    call verify(status == linalg_status_success, 'complex cosm')
    call verify(maxval(abs(matmul(sc, sc)+matmul(cc, cc)-cmplx(eye,kind=dp))) < &
        threshold, 'complex trigonometric identity')
    call tanm(ac, tc, status)
    call verify(status == linalg_status_success, 'complex tanm')
    call verify(maxval(abs(matmul(tc, cc)-sc)) < threshold, 'complex tangent identity')
    call sinhm(cmplx(0.0_dp, 1.0_dp, dp)*ac, xc, status)
    call verify(status == linalg_status_success .and. &
        maxval(abs(xc-cmplx(0.0_dp, 1.0_dp, dp)*sc)) < threshold, 'sinh(iA)=i*sin(A)')
    call coshm(cmplx(0.0_dp, 1.0_dp, dp)*ac, xc, status)
    call verify(status == linalg_status_success .and. maxval(abs(xc-cc)) < threshold, &
        'cosh(iA)=cos(A)')

    tiny_a = 0.0_dp
    tiny_a(1, 1) = 1.0e-120_dp
    tiny_a(2, 1) = -3.0e-120_dp
    call sinm(tiny_a, s, status)
    call verify(status == linalg_status_success .and. all(s == tiny_a), 'tiny sinm')
    call sinhm(tiny_a, s, status)
    call verify(status == linalg_status_success .and. all(s == tiny_a), 'tiny sinhm')

    a = 0.0_dp
    a(1, 1) = pi/2.0_dp
    call tanm(a, t, status)
    call verify(status == linalg_status_singular, 'tan pole')
    ac = cmplx(0.0_dp, 0.0_dp, dp)
    ac(1, 1) = cmplx(0.0_dp, pi/2.0_dp, dp)
    call tanhm(ac, tc, status)
    call verify(status == linalg_status_singular, 'complex tanh pole')

    a = 0.0_dp
    call sinm(a(:, 1:1), x, status)
    call verify(status == linalg_status_invalid_shape, 'nonsquare input')
    call sinhm(a, x(:, 1:1), status)
    call verify(status == linalg_status_invalid_shape, 'wrong output shape')
    a(1, 1) = ieee_value(1.0_dp, ieee_quiet_nan)
    call cosm(a, x, status)
    call verify(status == linalg_status_no_convergence, 'reject NaN input')
    a(1, 1) = ieee_value(1.0_dp, ieee_positive_inf)
    call sinhm(a, x, status)
    call verify(status == linalg_status_no_convergence, 'reject infinite input')
    a = 0.0_dp
    a(1, 1) = 1.0e100_dp
    call coshm(a, x, status)
    call verify(status == linalg_status_no_convergence, 'reject excessive scaling')
    block
        real(dp) :: large_a(41, 41), large_x(41, 41)
        large_a = 0.0_dp
        call sinm(large_a, large_x, status)
        call verify(status == linalg_status_invalid_shape, 'order larger than 40')
    end block

    print '(a)', 'linalg_trigonometric_native PASS'
contains
    subroutine verify(ok, message)
        logical, intent(in) :: ok !! expression expected to be true
        character(*), intent(in) :: message !! regression-case description
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine verify
end program test_linalg_trigonometric
