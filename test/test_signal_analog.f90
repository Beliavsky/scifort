! SPDX-License-Identifier: MIT
! Deterministic SciPy-reference tests for analog and digital ZPK responses.
program test_signal_analog
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_design, only : freqz, freqz_zpk, frequency_response, &
        design_status_success, design_status_invalid_input, design_status_singular
    use scifort_signal_analog, only : freqs, freqs_zpk
    implicit none

    real(dp), parameter :: tol = 8.0e-13_dp
    real(dp), parameter :: w(6) = [0.0_dp,0.2_dp,0.5_dp,1.0_dp,3.0_dp,5.0_dp]
    complex(dp), parameter :: p(2) = [cmplx(-1.0_dp,2.0_dp,dp), &
        cmplx(-1.0_dp,-2.0_dp,dp)]
    complex(dp), parameter :: analog_values(6) = [ &
        cmplx(1.8_dp,0.0_dp,dp), &
        cmplx(1.8175239079865593_dp,0.33729645903334193_dp,dp), &
        cmplx(1.9177718832891248_dp,0.8594164456233423_dp,dp), &
        cmplx(2.4_dp,1.8_dp,dp), &
        cmplx(5.5384615384615383_dp,-0.69230769230769229_dp,dp), &
        cmplx(3.84_dp,-1.08_dp,dp) ]
    complex(dp), parameter :: tf_values(6) = [ &
        cmplx(0.0_dp,0.0_dp,dp), &
        cmplx(0.028560548362528568_dp,0.093297791317593304_dp,dp), &
        cmplx(0.1411764705882353_dp,0.16470588235294117_dp,dp), &
        cmplx(0.3_dp,0.1_dp,dp), &
        cmplx(0.2076923076923077_dp,-0.16153846153846155_dp,dp), &
        cmplx(0.09946949602122017_dp,-0.15251989389920426_dp,dp) ]
    type(frequency_response) :: r, s
    real(dp) :: pi

    pi = acos(-1.0_dp)

    r = freqs_zpk([-3.0_dp,-1.0_dp],p,3.0_dp,w=w)
    if (r%status /= design_status_success) error stop 'complex-pole analog ZPK status'
    call equal_complex(r%response,analog_values,tol)
    call equal_real(r%frequencies,w,tol)

    r = freqs([1.0_dp,0.0_dp],[1.0_dp,3.0_dp,2.0_dp],w=w)
    if (r%status /= design_status_success) error stop 'analog TF status'
    call equal_complex(r%response,tf_values,tol)

    ! Log-spaced grid follows the pole and zero locations, not a fixed [0,pi] range.
    r = freqs([1.0_dp,0.0_dp],[1.0_dp,8.0_dp,25.0_dp],worn=9)
    if (r%status /= design_status_success) error stop 'analog automatic grid status'
    call equal_real(r%frequencies, &
        [0.01_dp,0.03162277660168379_dp,0.1_dp,0.31622776601683794_dp, &
        1.0_dp,3.1622776601683795_dp,10.0_dp,31.622776601683793_dp,100.0_dp],tol)
    r = freqs([1.0_dp],w=[-1.0_dp,0.0_dp,1.0_dp])
    if (r%status /= design_status_success) error stop 'analog default denominator'
    call equal_complex(r%response,[cmplx(1.0_dp,0.0_dp,dp), &
        cmplx(1.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.0_dp,dp)],tol)
    r = freqs([1.0_dp],[1.0_dp],worn=9)
    if (r%status /= design_status_success) error stop 'constant analog TF status'
    if (abs(r%frequencies(1)-100.0_dp) > tol) error stop 'constant analog start'
    if (abs(r%frequencies(9)-10000.0_dp) > tol) error stop 'constant analog end'
    r = freqs_zpk([real(dp) ::],[real(dp) ::],1.0_dp,worn=9)
    if (r%status /= design_status_success) error stop 'empty ZPK status'
    if (abs(r%frequencies(1)-100.0_dp) > tol) error stop 'empty ZPK start'

    ! Repeated poles have slow root convergence but still yield a usable grid.
    r = freqs([1.0_dp],[1.0_dp,4.0_dp,6.0_dp,4.0_dp,1.0_dp],worn=7)
    if (r%status /= design_status_success) error stop 'repeated pole grid status'
    if (size(r%response) /= 7) error stop 'repeated pole grid count'

    r = freqs_zpk([cmplx(-2.0_dp,0.0_dp,dp)],p,1.0_dp,worn=7)
    if (r%status /= design_status_success) error stop 'ZPK auto grid status'
    if (abs(r%frequencies(1)-0.1_dp) > tol) error stop 'ZPK auto start'
    if (abs(r%frequencies(7)-10.0_dp) > tol) error stop 'ZPK auto end'

    ! Coefficients of (s+2)/(s^2+2*s+5) agree with its roots.
    r = freqs([1.0_dp,2.0_dp],[1.0_dp,2.0_dp,5.0_dp],w=w)
    s = freqs_zpk([cmplx(-2.0_dp,0.0_dp,dp)],p,1.0_dp,w=w)
    call equal_complex(r%response,s%response,tol)

    ! Input gain may be complex, and real roots may be passed without conversion.
    r = freqs_zpk([-2.0_dp],[-1.0_dp],cmplx(2.0_dp,3.0_dp,dp),w=w)
    s = freqs([cmplx(2.0_dp,3.0_dp,dp),cmplx(4.0_dp,6.0_dp,dp)], &
        [cmplx(1.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.0_dp,dp)],w=w)
    call equal_complex(r%response,s%response,tol)

    ! Digital ZPK uses z=exp(i*2*pi*w/fs), with matching units for w and fs.
    r = freqz_zpk([1.0_dp,-1.0_dp],[0.5_dp],0.25_dp, &
        w=[0.0_dp,0.3_dp,0.9_dp,1.7_dp,3.2_dp])
    call equal_complex(r%response,[ &
        cmplx(0.0_dp,0.0_dp,dp), &
        cmplx(0.07409484829493179_dp,0.26192518349918714_dp,dp), &
        cmplx(0.24411631962791774_dp,0.4295620352396589_dp,dp), &
        cmplx(0.1783013059772498_dp,0.3827661822366977_dp,dp), &
        cmplx(0.0003789027864127623_dp,-0.019461737309172226_dp,dp) ],tol)
    if (r%status /= design_status_success) error stop 'digital ZPK status'
    r = freqz_zpk([-1.0_dp],[0.5_dp],0.25_dp,w=w)
    s = freqz([0.25_dp,0.25_dp],[1.0_dp,-0.5_dp],w=w)
    call equal_complex(r%response,s%response,tol)
    r = freqz_zpk([real(dp) ::],[real(dp) ::],1.0_dp,worn=5,fs=400.0_dp)
    if (r%status /= design_status_success) error stop 'digital empty ZPK status'
    call equal_real(r%frequencies,[0.0_dp,40.0_dp,80.0_dp,120.0_dp,160.0_dp],tol)
    r = freqz_zpk([real(dp) ::],[real(dp) ::],1.0_dp,worn=5,whole=.true.,fs=400.0_dp)
    call equal_real(r%frequencies,[0.0_dp,80.0_dp,160.0_dp,240.0_dp,320.0_dp],tol)
    r = freqz_zpk([0.0_dp],[-0.5_dp],cmplx(1.0_dp,2.0_dp,dp),w=[0.0_dp,pi/2])
    if (r%status /= design_status_success) error stop 'digital complex gain status'
    call equal_complex(r%response,[cmplx(1.0_dp,2.0_dp,dp)/cmplx(1.5_dp,0.0_dp,dp), &
        cmplx(1.0_dp,2.0_dp,dp)*cmplx(0.0_dp,1.0_dp,dp)/cmplx(0.5_dp,1.0_dp,dp)],tol)

    r = freqs([1.0_dp], [1.0_dp,0.0_dp],w=[0.0_dp,1.0_dp])
    if (r%status /= design_status_singular) error stop 'analog zero pole singularity'
    if (.not. ieee_is_nan(real(r%response(1),dp))) error stop 'analog singular NaN'
    r = freqz_zpk([real(dp) ::],[1.0_dp],1.0_dp,w=[0.0_dp,1.0_dp])
    if (r%status /= design_status_singular) error stop 'digital unit pole singularity'
    if (.not. ieee_is_nan(real(r%response(1),dp))) error stop 'digital singular NaN'

    r = freqs([1.0_dp],[0.0_dp],w=w)
    if (r%status /= design_status_invalid_input) error stop 'invalid analog denominator'
    r = freqs([1.0_dp],[1.0_dp],worn=0)
    if (r%status /= design_status_invalid_input) error stop 'invalid analog count'
    r = freqs_zpk([real(dp) ::],[real(dp) ::],1.0_dp,w=[ &
        ieee_value(0.0_dp,ieee_quiet_nan)])
    if (r%status /= design_status_invalid_input) error stop 'invalid analog frequencies'
    r = freqz_zpk([real(dp) ::],[real(dp) ::],1.0_dp,fs=-1.0_dp)
    if (r%status /= design_status_invalid_input) error stop 'invalid digital rate'
    r = freqs([1.0_dp],[1.0_dp],w=[real(dp) ::])
    if (r%status /= design_status_success) error stop 'empty analog frequency success'
    if (size(r%response) /= 0) error stop 'empty analog frequency response'

    print '(a)', 'test_signal_analog passed'
contains
    subroutine equal_complex(actual, expected, eps)
        complex(dp), intent(in) :: actual(:) !! observed complex values
        complex(dp), intent(in) :: expected(:) !! SciPy reference complex values
        real(dp), intent(in) :: eps !! maximum absolute discrepancy
        if (size(actual) /= size(expected)) error stop 'complex shape mismatch'
        if (any(abs(actual-expected) > eps)) then
            print *, 'actual:',actual
            print *, 'expected:',expected
            error stop 'complex response mismatch'
        end if
    end subroutine equal_complex
    subroutine equal_real(actual, expected, eps)
        real(dp), intent(in) :: actual(:) !! observed real values
        real(dp), intent(in) :: expected(:) !! SciPy reference real values
        real(dp), intent(in) :: eps !! maximum absolute discrepancy
        if (size(actual) /= size(expected)) error stop 'frequency shape mismatch'
        if (any(abs(actual-expected) > eps)) then
            print *, 'actual:',actual
            print *, 'expected:',expected
            error stop 'frequency mismatch'
        end if
    end subroutine equal_real
end program test_signal_analog
