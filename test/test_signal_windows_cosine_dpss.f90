! SPDX-License-Identifier: MIT
! Regression for coefficient-vector general cosine and get_window DPSS.
program test_signal_windows_cosine_dpss
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window, general_cosine, &
        spectral_status_success, spectral_status_invalid_input
    use scifort_signal_dpss, only : dpss, dpss_result, dpss_status_success
    implicit none
    real(dp), allocatable :: w(:), v(:), coeff(:)
    type(dpss_result) :: taper
    integer :: code, n

    coeff = [0.42_dp, 0.5_dp, 0.08_dp]
    w = general_cosine(9, coeff, status=code)
    if (code /= spectral_status_success) error stop 'cosine direct status'
    call close(w, [0.0_dp, 0.06644660940672624_dp, 0.34_dp, &
                   0.7735533905932737_dp, 1.0_dp, 0.7735533905932737_dp, &
                   0.34_dp, 0.06644660940672624_dp, 0.0_dp], 1.0e-12_dp)
    v = get_window('general_cosine', 9, coefficients=coeff, fftbins=.false., status=code)
    if (code /= spectral_status_success) error stop 'cosine get_window status'
    call close(w,v,1.0e-14_dp)
    v = get_window('general_cosine', 10, coefficients=coeff, status=code)
    w = general_cosine(10, coeff, sym=.false.)
    if (code /= spectral_status_success) error stop 'cosine periodic status'
    call close(w,v,1.0e-14_dp)
    w = general_cosine(8,[2.0_dp])
    if (any(w /= 2.0_dp)) error stop 'cosine constant term'
    w = general_cosine(4,[real(dp)::])
    if (any(w /= 0.0_dp)) error stop 'cosine empty coefficients'
    w = general_cosine(0,coeff,status=code)
    if (size(w) /= 0 .or. code /= spectral_status_success) error stop 'cosine zero-length'
    w = general_cosine(1,[real(dp)::],status=code)
    if (size(w) /= 1 .or. w(1) /= 1.0_dp) error stop 'cosine unit-length'
    w = general_cosine(-1,coeff,status=code)
    call invalid(w,code)
    w = general_cosine(4,[ieee_value(0.0_dp,ieee_quiet_nan)],status=code)
    call invalid(w,code)
    w = get_window('general_cosine',10,status=code)
    call invalid(w,code)
    w = get_window('general_cosine',10,coefficients=coeff,parameter=1.0_dp,status=code)
    call invalid(w,code)
    w = get_window('hann',10,coefficients=coeff,status=code)
    call invalid(w,code)

    do n = 2, 17
        taper = dpss(n,0.7_dp,kmax=1,sym=.false.,norm='approximate')
        if (taper%status /= dpss_status_success) error stop 'direct DPSS status'
        w = get_window('dpss',n,parameter=0.7_dp,status=code)
        if (code /= spectral_status_success) error stop 'get_window DPSS status'
        call close(w,taper%tapers(:,1),1.0e-13_dp)
    end do
    w = get_window('dpss', 1, parameter=2.0_dp,status=code)
    if (code /= spectral_status_success .or. size(w) /= 1) error stop 'single DPSS'
    if (w(1) /= 1.0_dp) error stop 'single DPSS value'
    w = get_window('slepian',12,fftbins=.false.,parameter=2.5_dp,status=code)
    taper = dpss(12,2.5_dp,sym=.true.,norm='approximate')
    call close(w,taper%tapers(:,1),1.0e-12_dp)
    w = get_window('dpss', 0, parameter=2.0_dp,status=code)
    call invalid(w,code)
    w = get_window('dpss', 12, status=code)
    call invalid(w,code)
    w = get_window('dpss', 12, parameter=6.0_dp,status=code)
    call invalid(w,code)
    w = get_window('dpss', 12, parameter=2.0_dp,parameter2=0.5_dp,status=code)
    call invalid(w,code)
    w = get_window('dpss', 12, parameter=2.0_dp,coefficients=coeff,status=code)
    call invalid(w,code)

    print '(a)', 'test_signal_windows_cosine_dpss passed'
contains
    subroutine close(actual,expected,tol)
        real(dp), intent(in) :: actual(:) !! computed samples
        real(dp), intent(in) :: expected(:) !! reference samples
        real(dp), intent(in) :: tol !! absolute tolerance
        if (size(actual) /= size(expected)) error stop 'window length mismatch'
        if (any(abs(actual-expected) > tol)) error stop 'window values mismatch'
    end subroutine close

    subroutine invalid(actual,stat)
        real(dp), intent(in) :: actual(:) !! invalid result
        integer, intent(in) :: stat !! failure status
        if (size(actual) /= 0 .or. stat /= spectral_status_invalid_input) error stop 'expected invalid'
    end subroutine invalid
end program test_signal_windows_cosine_dpss
