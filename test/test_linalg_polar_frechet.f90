! SPDX-License-Identifier: MIT
program test_linalg_polar_frechet
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_polar, only : polar
    use scifort_linalg_expm_frechet, only : expm_frechet
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(3,2), ur(3,2), pr(2,2), pl(3,3), err, orth
    real(dp) :: ar(2,2), er(2,2), xr(2,2), lr(2,2)
    complex(dp) :: ac(2,3), uc(2,3), pc(3,3), tc(2,2), &
        ec(2,2), xc(2,2), lc(2,2), zero(2,2), xp(2,2), xm(2,2), work(2,2)
    integer :: status, tcount, scount
    real(dp), parameter :: h = 1.0e-5_dp

    a(:,1) = [1.0_dp, 2.0_dp, -0.5_dp]
    a(:,2) = [0.75_dp, 3.0_dp, 1.0_dp]
    call polar(a, ur, pr, status, residual=err, orthogonality=orth)
    if (status /= linalg_status_success .or. err > 1.0e-11_dp .or. orth > 1.0e-11_dp) error stop 'tall polar'
    call polar(a, ur, pl, status, side='left', residual=err)
    if (status /= linalg_status_success .or. err > 1.0e-11_dp) error stop 'left polar'
    call polar(a, ur, pr, status, side='broken')
    if (status /= linalg_status_invalid_shape) error stop 'invalid side'
    call polar(a, ur, pr, status, max_sweeps=1)
    if (status /= linalg_status_no_convergence) error stop 'polar sweep limit'
    a(:,2) = 2.0_dp*a(:,1)
    call polar(a, ur, pr, status, orthogonality=orth, residual=err)
    if (status /= linalg_status_success .or. err > 1.0e-10_dp .or. orth > 1.0e-10_dp) error stop 'rank deficient'
    a = 0.0_dp
    call polar(a, ur, pr, status, residual=err, orthogonality=orth)
    if (status /= linalg_status_success .or. err > 1.0e-14_dp .or. orth > 1.0e-14_dp) error stop 'zero polar'

    ac(1,:) = [cmplx(1.0_dp, 1.0_dp,dp), cmplx(0.0_dp,2.0_dp,dp), cmplx(2.0_dp,-1.0_dp,dp)]
    ac(2,:) = [cmplx(0.3_dp,-1.0_dp,dp), cmplx(1.2_dp,0.5_dp,dp), cmplx(1.0_dp,1.0_dp,dp)]
    call polar(ac, uc, pc, status, residual=err, orthogonality=orth)
    if (status /= linalg_status_success .or. err > 1.0e-10_dp .or. orth > 1.0e-10_dp) error stop 'wide complex polar'

    ar = 0.0_dp
    er = 0.0_dp
    ar(1,1) = 1.0_dp
    ar(2,2) = -0.5_dp
    er(1,2) = 1.0_dp
    call expm_frechet(ar, er, xr, lr, status, terms=tcount, scalings=scount)
    if (status /= linalg_status_success .or. tcount <= 0 .or. scount < 0) error stop 'real Frechet'
    if (abs(xr(1,1)-exp(1.0_dp)) > 1.0e-11_dp .or. &
            abs(xr(2,2)-exp(-0.5_dp)) > 1.0e-11_dp) error stop 'exponential values'
    if (abs(lr(1,2)-(exp(1.0_dp)-exp(-0.5_dp))/1.5_dp) > 1.0e-11_dp) error stop 'divided difference'
    if (maxval(abs(lr(2,:))) > 1.0e-14_dp) error stop 'Frechet structure'

    tc(1,:) = [cmplx(0.3_dp,0.7_dp,dp), cmplx(0.2_dp,-0.1_dp,dp)]
    tc(2,:) = [cmplx(0.6_dp,-0.3_dp,dp), cmplx(-0.8_dp,0.3_dp,dp)]
    ec(1,:) = [cmplx(0.0_dp,0.2_dp,dp), cmplx(-1.0_dp,0.3_dp,dp)]
    ec(2,:) = [cmplx(0.5_dp,0.4_dp,dp), cmplx(0.1_dp,-0.6_dp,dp)]
    call expm_frechet(tc, ec, xc, lc, status)
    if (status /= linalg_status_success) error stop 'complex Frechet'
    zero = cmplx(0.0_dp,0.0_dp,dp)
    call expm_frechet(tc+h*ec, zero, xp, work, status)
    if (status /= linalg_status_success) error stop 'forward exponential'
    call expm_frechet(tc-h*ec, zero, xm, work, status)
    if (status /= linalg_status_success) error stop 'backward exponential'
    if (maxval(abs((xp-xm)/(2.0_dp*h)-lc)) > 1.0e-9_dp) error stop 'Frechet finite difference'
    call expm_frechet(tc, zero, xc, lc, status)
    if (status /= linalg_status_success .or. maxval(abs(lc)) > 1.0e-14_dp) error stop 'zero direction'
    call expm_frechet(tc, ec, xc, lc, status, max_terms=1)
    if (status /= linalg_status_no_convergence) error stop 'Frechet limit'
    call expm_frechet(tc, ec, xc, lc, status, tol=-1.0_dp)
    if (status /= linalg_status_invalid_shape) error stop 'Frechet invalid tolerance'
    ar(1,1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call expm_frechet(ar, er, xr, lr, status)
    if (status /= linalg_status_invalid_shape) error stop 'Frechet nonfinite input'
    ar(1,1) = 0.0_dp
    a(1,1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call polar(a, ur, pr, status)
    if (status /= linalg_status_invalid_shape) error stop 'polar nonfinite input'
    print '(a)', 'PASS test_linalg_polar_frechet'
end program test_linalg_polar_frechet
