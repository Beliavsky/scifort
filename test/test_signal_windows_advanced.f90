! SPDX-License-Identifier: MIT
! Regression and edge checks for SciPy-style Chebyshev, Taylor, and KBD windows.
program test_signal_windows_advanced
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window, chebwin, taylor, &
        kaiser_bessel_derived, spectral_status_success, spectral_status_invalid_input
    implicit none
    real(dp), allocatable :: w(:), v(:)
    integer :: st, i

    w = chebwin(9, 85.0_dp, status=st)
    if (st /= spectral_status_success) error stop 'chebwin direct status'
    call close(w, [ &
        2.72604933080068217e-02_dp,1.62892271132077909e-01_dp,4.67062894701336273e-01_dp, &
        8.31319293339228604e-01_dp,1.00000000000000000e+00_dp,8.31319293339228604e-01_dp, &
        4.67062894701336273e-01_dp,1.62892271132077909e-01_dp,2.72604933080068217e-02_dp])
    w = get_window('chebwin', 8, parameter=75.0_dp, status=st)
    if (st /= spectral_status_success) error stop 'chebwin get_window status'
    call close(w, [ &
        3.35456752513922721e-02_dp,1.81762495262913409e-01_dp,4.89529484194940301e-01_dp, &
        8.40948866178676480e-01_dp,1.00000000000000000e+00_dp,8.40948866178676480e-01_dp, &
        4.89529484194940301e-01_dp,1.81762495262913409e-01_dp])

    w = taylor(11, nbar=5, sll=48.0_dp, status=st)
    if (st /= spectral_status_success) error stop 'taylor direct status'
    call close(w, [ &
        8.27813225407344061e-02_dp,2.26434284361667515e-01_dp,4.56362623158694469e-01_dp, &
        7.13805128271927658e-01_dp,9.21236880857823071e-01_dp,1.00000000000000000e+00_dp, &
        9.21236880857823071e-01_dp,7.13805128271927658e-01_dp,4.56362623158694469e-01_dp, &
        2.26434284361667459e-01_dp,8.27813225407344061e-02_dp])
    w = taylor(12, nbar=4, sll=40.0_dp, norm=.false., sym=.false.)
    call close(w, [ &
        2.37204017511008503e-01_dp,4.32019557331678405e-01_dp,7.51773451194250208e-01_dp, &
        1.10844459307354026e+00_dp,1.43214583080759872e+00_dp,1.66381550107736143e+00_dp, &
        1.74919409800912495e+00_dp,1.66381550107736143e+00_dp,1.43214583080759872e+00_dp, &
        1.10844459307354026e+00_dp,7.51773451194250208e-01_dp,4.32019557331678405e-01_dp])
    w = get_window('taylor', 12, parameter=4.0_dp, parameter2=40.0_dp)
    v = taylor(12, nbar=4, sll=40.0_dp, sym=.false.)
    call close(w, v)
    w = get_window('taylor', 11)
    v = taylor(11, sym=.false.)
    call close(w, v)

    w = kaiser_bessel_derived(10, 6.0_dp, status=st)
    if (st /= spectral_status_success) error stop 'kbd direct status'
    call close(w, [ &
        7.71537974818444355e-02_dp,3.76346686479013537e-01_dp,7.07106781186547462e-01_dp, &
        9.26478910486507989e-01_dp,9.97019203192260783e-01_dp,9.97019203192260783e-01_dp, &
        9.26478910486507989e-01_dp,7.07106781186547462e-01_dp,3.76346686479013537e-01_dp, &
        7.71537974818444355e-02_dp])
    do i = 1, 5
        if (abs(w(i)**2 + w(i+5)**2 - 1.0_dp) > 1.0e-13_dp) error stop 'kbd PB condition'
    end do
    v = get_window('kbd', 10, fftbins=.false., parameter=6.0_dp, status=st)
    if (st /= spectral_status_success) error stop 'get_window kbd status'
    call close(w,v)

    w = chebwin(0, 60.0_dp, status=st)
    if (size(w) /= 0 .or. st /= spectral_status_success) error stop 'empty chebwin'
    w = taylor(1, status=st)
    if (size(w) /= 1 .or. w(1) /= 1.0_dp) error stop 'unit taylor'
    w = kaiser_bessel_derived(2, 0.0_dp)
    call close(w,[sqrt(0.5_dp),sqrt(0.5_dp)])
    w = chebwin(-1, 60.0_dp, status=st)
    call invalid(w,st)
    w = get_window('chebwin', 7, status=st)
    call invalid(w,st)
    w = get_window('taylor', 8, parameter=3.5_dp, status=st)
    call invalid(w,st)
    w = taylor(10, nbar=0, status=st)
    call invalid(w,st)
    w = taylor(10, sll=0.0_dp, status=st)
    call invalid(w,st)
    w = kaiser_bessel_derived(9, 4.0_dp, status=st)
    call invalid(w,st)
    w = get_window('kbd', 10, parameter=4.0_dp, status=st)
    call invalid(w,st)
    w = get_window('kbd', 10, fftbins=.false., status=st)
    call invalid(w,st)
    w = chebwin(8, ieee_value(0.0_dp,ieee_quiet_nan), status=st)
    call invalid(w,st)

    print '(a)', 'test_signal_windows_advanced passed'
contains
    subroutine close(actual, expected)
        real(dp), intent(in) :: actual(:) !! actual window values
        real(dp), intent(in) :: expected(:) !! trusted SciPy 1.17 reference
        if (size(actual) /= size(expected)) error stop 'size mismatch'
        if (any(abs(actual-expected) > 5.0e-11_dp)) error stop 'numerical window mismatch'
    end subroutine close

    subroutine invalid(actual, status)
        real(dp), intent(in) :: actual(:) !! zero-length invalid result
        integer, intent(in) :: status !! invalid-input code
        if (status /= spectral_status_invalid_input .or. size(actual) /= 0) error stop 'invalid input'
    end subroutine invalid
end program test_signal_windows_advanced
