! SPDX-License-Identifier: MIT
! Self-contained tests of the SciPy-compatible signal waveforms.
program test_signal_waveforms
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_waveforms, only : chirp, chirp_complex, sweep_poly, &
        gausspulse, gausspulse_cutoff, waveform_status_success, waveform_status_invalid_input
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: t(5) = [-0.25_dp, 0.0_dp, 0.3_dp, 1.1_dp, 1.75_dp]
    real(dp), allocatable :: y(:), q(:), env(:)
    complex(dp), allocatable :: z(:)
    real(dp) :: cutoff, nan
    integer :: status

    nan = ieee_value(0.0_dp, ieee_quiet_nan)

    y = chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, phi=17.0_dp, status=status)
    if (status /= waveform_status_success) error stop 'chirp linear status'
    call check_real(y, [-0.1856666153855782_dp, 0.9563047559630354_dp, &
        0.391245174417658_dp, -0.74770538162227_dp, 0.9690157314068709_dp])
    call check_real(y, chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='lin', phi=17.0_dp))
    call check_real(y, chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='li', phi=17.0_dp))

    y = chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='quadratic', phi=17.0_dp)
    call check_real(y, [0.3580420466734027_dp, 0.9563047559630354_dp, &
        0.8860301524662683_dp, -0.03602138639618806_dp, 0.7447775998868769_dp])
    call check_real(y, chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='q', phi=17.0_dp))
    y = chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='quadratic', phi=17.0_dp, vertex_zero=.false.)
    call check_real(y, [-0.674044576966651_dp, 0.9563047559630354_dp, &
        -0.280910605518205_dp, -0.9966091837231887_dp, 0.3344648657410446_dp])
    call check_real(y, chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='quad', phi=17.0_dp, vertex_zero=.false.))

    y = chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='logarithmic', phi=17.0_dp)
    call check_real(y, [0.05359251353878756_dp, 0.9563047559630354_dp, &
        0.6556892347574105_dp, 0.5019831913954672_dp, 0.9165791216701943_dp])
    call check_real(y, chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='lo', phi=17.0_dp))
    y = chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='hyperbolic', phi=17.0_dp)
    call check_real(y, [0.162041625176535_dp, 0.9563047559630354_dp, &
        0.7709848387491898_dp, -0.9763804980474584_dp, 0.9814047077801192_dp])
    call check_real(y, chirp(t, 3.4_dp, 2.5_dp, 10.0_dp, method='hyp', phi=17.0_dp))
    call check_real(chirp(t, -4.0_dp, 2.0_dp, -4.0_dp, method='log'), cos(-8.0_dp*pi*t))
    call check_real(chirp(t, -4.0_dp, 2.0_dp, -4.0_dp, method='hyp'), cos(-8.0_dp*pi*t))

    z = chirp_complex(t, 3.4_dp, 2.5_dp, 10.0_dp, method='hyperbolic', phi=17.0_dp, status=status)
    if (status /= waveform_status_success) error stop 'chirp complex status'
    call check_real(real(z,dp), y)
    call check_real(abs(z), spread(1.0_dp,1,size(t)))
    call check_complex(z(2), cmplx(cos(17.0_dp*pi/180.0_dp),sin(17.0_dp*pi/180.0_dp),dp))
    z = chirp_complex(t, 3.4_dp, 2.5_dp, 10.0_dp, method='quadratic', vertex_zero=.false.)
    call check_real(real(z,dp), chirp(t,3.4_dp,2.5_dp,10.0_dp,method='quad',vertex_zero=.false.))

    y = sweep_poly(t, [0.125_dp,-0.3_dp,1.7_dp,2.3_dp], phi=-15.0_dp, status=status)
    if (status /= waveform_status_success) error stop 'sweep_poly status'
    call check_real(y, [-0.9254178337724231_dp, 0.9659258262890683_dp, &
        -0.1726318367440497_dp, -0.9034484921337013_dp, -0.5548472634484219_dp])
    y = sweep_poly(t, [2.5_dp])
    call check_real(y, cos(2.0_dp*pi*2.5_dp*t))
    y = sweep_poly(t, [real(dp) ::], phi=60.0_dp)
    call check_real(y, spread(0.5_dp,1,size(t)))

    y = gausspulse(t, fc=3.4_dp, bw=0.63_dp, bwr=-10.0_dp, &
        quadrature=q, envelope=env, status=status)
    if (status /= waveform_status_success) error stop 'gausspulse status'
    call check_real(y, [0.3179166756634496_dp, 1.0_dp, 0.4094662796703614_dp, &
        -4.271751516937617e-07_dp, 7.940408840717027e-14_dp])
    call check_real(q, [0.4375747645994702_dp, 0.0_dp, 0.05172762060539703_dp, &
        -6.789754604781591e-06_dp, -2.579995228490889e-14_dp])
    call check_real(env, [0.5408721542834136_dp, 1.0_dp, 0.412720705708573_dp, &
        6.803179124745827e-06_dp, 8.349039941021535e-14_dp])
    call check_real(y*y + q*q, env*env)
    cutoff = gausspulse_cutoff(fc=3.4_dp, bw=0.63_dp, bwr=-10.0_dp, tpr=-40.0_dp, status=status)
    if (status /= waveform_status_success) error stop 'gausspulse cutoff status'
    call check_scalar(cutoff, 0.6843469644065618_dp)
    y = gausspulse([0.0_dp,0.1_dp],fc=0.0_dp,quadrature=q,envelope=env)
    call check_real(y,[1.0_dp,1.0_dp])
    call check_real(q,[0.0_dp,0.0_dp])
    call check_real(env,[1.0_dp,1.0_dp])

    y = chirp([real(dp) ::],1.0_dp,1.0_dp,2.0_dp,status=status)
    if (size(y) /= 0 .or. status /= waveform_status_success) error stop 'empty chirp'
    y = sweep_poly([real(dp) ::], [1.0_dp],status=status)
    if (size(y) /= 0 .or. status /= waveform_status_success) error stop 'empty sweep'
    y = gausspulse([real(dp) ::],quadrature=q,envelope=env,status=status)
    if (size(y) /= 0 .or. size(q) /= 0 .or. size(env) /= 0) error stop 'empty pulse size'
    if (status /= waveform_status_success) error stop 'empty pulse status'

    y = chirp(t,1.0_dp,0.0_dp,2.0_dp,status=status)
    call reject_real(y,status)
    y = chirp(t,1.0_dp,1.0_dp,2.0_dp,method='bad',status=status)
    call reject_real(y,status)
    y = chirp(t,1.0_dp,1.0_dp,-2.0_dp,method='logarithmic',status=status)
    call reject_real(y,status)
    y = chirp(t,0.0_dp,1.0_dp,2.0_dp,method='hyperbolic',status=status)
    call reject_real(y,status)
    y = chirp([2.0_dp],2.0_dp,1.0_dp,4.0_dp,method='hyp',status=status)
    call reject_real(y,status)
    y = chirp([nan],1.0_dp,1.0_dp,2.0_dp,status=status)
    call reject_real(y,status)
    z = chirp_complex(t,1.0_dp,1.0_dp,2.0_dp,phi=nan,status=status)
    if (status /= waveform_status_invalid_input) error stop 'chirp complex nan status'
    if (.not. all(ieee_is_nan(real(z)) .and. ieee_is_nan(aimag(z)))) error stop 'chirp complex nan'
    y = sweep_poly(t,[1.0_dp,nan],status=status)
    call reject_real(y,status)
    y = sweep_poly(t,[1.0_dp],phi=nan,status=status)
    call reject_real(y,status)
    y = gausspulse(t,fc=-1.0_dp,status=status)
    call reject_real(y,status)
    y = gausspulse(t,bw=0.0_dp,status=status)
    call reject_real(y,status)
    y = gausspulse(t,bwr=0.0_dp,status=status)
    call reject_real(y,status)
    y = gausspulse(t,fc=nan,quadrature=q,envelope=env,status=status)
    call reject_real(y,status)
    if (.not. all(ieee_is_nan(q))) error stop 'pulse bad quadrature'
    if (.not. all(ieee_is_nan(env))) error stop 'pulse bad envelope'
    cutoff = gausspulse_cutoff(fc=0.0_dp,status=status)
    if (status /= waveform_status_invalid_input .or. .not. ieee_is_nan(cutoff)) error stop 'cutoff zero frequency'
    cutoff = gausspulse_cutoff(tpr=0.0_dp,status=status)
    if (status /= waveform_status_invalid_input .or. .not. ieee_is_nan(cutoff)) error stop 'cutoff nonnegative'

    print '(a)', 'test_signal_waveforms passed'

contains

    subroutine check_real(actual, expected)
        real(dp), intent(in) :: actual(:) !! observed array
        real(dp), intent(in) :: expected(:) !! reference array
        if (size(actual) /= size(expected)) error stop 'waveform real size mismatch'
        if (any(.not. abs(actual-expected) < 8.0e-13_dp)) then
            print '(a,*(1x,es22.14))', 'obtained:', actual
            print '(a,*(1x,es22.14))', 'expected:', expected
            error stop 'waveform real mismatch'
        end if
    end subroutine check_real

    subroutine check_scalar(actual, expected)
        real(dp), intent(in) :: actual !! observed scalar
        real(dp), intent(in) :: expected !! reference scalar
        if (.not. abs(actual-expected) < 8.0e-13_dp) error stop 'waveform scalar mismatch'
    end subroutine check_scalar

    subroutine check_complex(actual, expected)
        complex(dp), intent(in) :: actual !! observed complex number
        complex(dp), intent(in) :: expected !! reference complex number
        if (.not. abs(actual-expected) < 8.0e-13_dp) error stop 'waveform complex mismatch'
    end subroutine check_complex

    subroutine reject_real(actual, actual_status)
        real(dp), intent(in) :: actual(:) !! invalid array
        integer, intent(in) :: actual_status !! expected nonzero status
        if (actual_status /= waveform_status_invalid_input) error stop 'waveform invalid status'
        if (.not. all(ieee_is_nan(actual))) error stop 'waveform invalid value'
    end subroutine reject_real
end program test_signal_waveforms
