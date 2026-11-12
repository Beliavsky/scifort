! SPDX-License-Identifier: MIT
! Cubic/quadratic spline coefficients and mirror-boundary evaluation.
program test_signal_spline
    use scifort_signal, only : dp, cspline1d, qspline1d, cspline1d_eval, qspline1d_eval, &
        spline_status_success, spline_status_invalid_input
    implicit none
    real(dp), parameter :: x(6) = [1.0_dp,0.25_dp,2.0_dp,-1.0_dp,0.75_dp,3.0_dp]
    real(dp), parameter :: coords(7) = [-2.0_dp,-0.5_dp,0.0_dp,0.75_dp,2.5_dp,5.0_dp,7.0_dp]
    real(dp), parameter :: refc(6) = [ &
        1.390846497464884_dp,-0.9520468370355863_dp,3.9173408506774603_dp, &
        -2.717316565674251_dp,0.9519254120195397_dp,3.409614917596093_dp]
    real(dp), parameter :: refq(6) = [ &
        1.1965543244256016_dp,-0.37576058306747806_dp,3.0580091739792694_dp, &
        -1.972294460808157_dp,0.775757590869685_dp,3.3177489155900424_dp]
    real(dp), parameter :: valc(7) = [ &
        2.0_dp,0.320845407042004_dp,1.000364275048139_dp,0.13468872123865755_dp, &
        0.5750091068762035_dp,3.0_dp,-1.0_dp]
    real(dp), parameter :: valq(7) = [ &
        2.0_dp,0.4103968706790618_dp,1.0000149609889666_dp,0.17375828957266146_dp, &
        0.5428573565855562_dp,3.0_dp,-1.0_dp]
    real(dp), allocatable :: rc(:), rq(:), value(:)
    complex(dp), allocatable :: cc(:), qc(:), cv(:)
    complex(dp) :: sample(6)
    real(dp) :: v
    integer :: status

    rc = cspline1d(x, status=status)
    if (status /= spline_status_success) error stop 'real cubic status'
    if (maxval(abs(rc-refc)) > 3e-13_dp) error stop 'real cubic coefficients'
    rq = qspline1d(x, status=status)
    if (status /= spline_status_success) error stop 'real quadratic status'
    if (maxval(abs(rq-refq)) > 3e-13_dp) error stop 'real quadratic coefficients'

    value = cspline1d_eval(rc, coords, status=status)
    if (status /= spline_status_success) error stop 'cubic eval status'
    if (maxval(abs(value-valc)) > 2e-13_dp) error stop 'cubic eval samples'
    value = qspline1d_eval(rq, coords, status=status)
    if (status /= spline_status_success) error stop 'quadratic eval status'
    if (maxval(abs(value-valq)) > 2e-13_dp) error stop 'quadratic eval samples'

    v = cspline1d_eval(rc, -2.0_dp, status=status)
    if (status /= spline_status_success .or. abs(v-2.0_dp) > 3e-13_dp) &
        error stop 'scalar cubic eval'
    v = qspline1d_eval(rq, 5.0_dp, status=status)
    if (status /= spline_status_success .or. abs(v-3.0_dp) > 3e-13_dp) &
        error stop 'scalar quadratic eval'

    ! An origin and knot spacing affect coordinates, not spline weights.
    value = cspline1d_eval(rc, 1.25_dp + 2.0_dp * coords, dx=2.0_dp, x0=1.25_dp)
    if (maxval(abs(value-valc)) > 2e-13_dp) error stop 'knot spacing'

    sample = cmplx(x,2.0_dp*x,dp)
    cc = cspline1d(sample, status=status)
    qc = qspline1d(sample, status=status)
    if (maxval(abs(cc-cmplx(refc,2.0_dp*refc,dp))) > 5e-13_dp) &
        error stop 'complex cubic coefficients'
    if (maxval(abs(qc-cmplx(refq,2.0_dp*refq,dp))) > 5e-13_dp) &
        error stop 'complex quadratic coefficients'
    cv = qspline1d_eval(qc, coords, status=status)
    if (status /= spline_status_success) error stop 'complex quadratic eval status'
    if (maxval(abs(cv-cmplx(valq,2.0_dp*valq,dp))) > 5e-13_dp) &
        error stop 'complex quadratic eval'
    cv = cspline1d_eval(cc, coords, status=status)
    if (maxval(abs(cv-cmplx(valc,2.0_dp*valc,dp))) > 5e-13_dp) &
        error stop 'complex cubic eval'

    rc = cspline1d(x, lamb=0.2_dp, status=status)
    if (status /= spline_status_success) error stop 'smoothed cubic status'
    if (maxval(abs(rc-[0.9976871901547613_dp,1.0578290903976972_dp, &
                     1.140487936275733_dp,-0.2738979384314078_dp, &
                     0.8402073858316423_dp,2.7707535216354655_dp])) > 2e-12_dp) &
        error stop 'cubic smooth shape'

    rc = cspline1d(x, lamb=0.02_dp, status=status)
    if (status /= spline_status_invalid_input .or. size(rc) /= 0) &
        error stop 'unstable cubic smoothing'

    rq = qspline1d(x, lamb=0.2_dp, status=status)
    if (status /= spline_status_invalid_input .or. size(rq) /= 0) &
        error stop 'unsupported quadratic smoothing'
    rc = cspline1d([real(dp)::], status=status)
    if (status /= spline_status_invalid_input .or. size(rc) /= 0) &
        error stop 'empty input'
    value = cspline1d_eval(refc, [0.0_dp,2.0_dp], dx=0.0_dp, status=status)
    if (status /= spline_status_invalid_input .or. size(value) /= 0) &
        error stop 'zero sample spacing'
    value = qspline1d_eval(refq, [real(dp)::], status=status)
    if (status /= spline_status_success .or. size(value) /= 0) &
        error stop 'empty query'

    print '(a)', 'test_signal_spline passed'
end program test_signal_spline
