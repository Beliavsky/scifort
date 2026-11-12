! SPDX-License-Identifier: MIT
! Complex continuous-wavelet transform and Morlet regressions.
program test_signal_cwt_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_cwt_peaks, only : cwt_complex_wavelet, cwt_morlet2, &
        morlet2_wavelet, cwt_peak_success, cwt_peak_invalid_input
    implicit none

    real(dp) :: x(24), widths(4), nanval
    complex(dp) :: z(24), expected
    complex(dp), allocatable :: y(:,:), yr(:,:), yc(:,:), w(:)
    integer :: i, j, k, m, off, ix, status

    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    widths = [0.5_dp, 1.0_dp, 2.5_dp, 4.7_dp]
    do i = 1, size(x)
        x(i) = sin(0.13_dp*i)+0.2_dp*cos(0.37_dp*i)
        z(i) = cmplx(x(i),cos(0.11_dp*i)-0.5_dp*sin(0.21_dp*i),dp)
    end do

    y = cwt_morlet2(x,widths,status=status)
    call require(status == cwt_peak_success, 'real input status')
    call require(all(shape(y)==[size(widths),size(x)]), 'real input shape')
    do i = 1, size(widths)
        w = morlet2_wavelet(min(10.0_dp*widths(i),real(size(x),dp)),widths(i))
        call require(size(w) == ceiling(min(10.0_dp*widths(i),real(size(x),dp))), &
            'wavelet shape')
        m = size(w)
        off = (m-1)/2
        do j = 1, size(x)
            expected = (0.0_dp,0.0_dp)
            do k = 1, m
                ix = j+off+1-k
                if (ix < 1 .or. ix > size(x)) cycle
                expected = expected + cmplx(x(ix),0.0_dp,dp)*conjg(w(m+1-k))
            end do
            call require(abs(y(i,j)-expected) < 1.0e-13_dp, 'Morlet direct convolution')
        end do
    end do
    yr = cwt_complex_wavelet(x,widths,morlet2_wavelet,status)
    call require(status == cwt_peak_success, 'wavelet callback status')
    call require(maxval(abs(y-yr)) < 1.0e-13_dp,'Morlet callback and direct mode differ')

    yc = cwt_morlet2(z(1:24:2),widths,omega0=2.3_dp,status=status)
    call require(status == cwt_peak_success, 'complex stride status')
    call require(all(shape(yc)==[4,12]), 'complex stride shape')
    yr = cwt_morlet2(z(1:24:2),widths,omega0=2.3_dp)
    call require(maxval(abs(yc-yr)) < 1.0e-13_dp,'complex stride input')
    yc = cwt_complex_wavelet(z,widths,asymmetric_wavelet,status)
    call require(status == cwt_peak_success, 'custom complex wavelet status')
    call require(all(shape(yc)==[4,24]), 'custom complex wavelet shape')

    ! At a single scale an asymmetric complex wavelet must be conjugated
    ! before the convolution, not just time-reversed.
    w = asymmetric_wavelet(5.0_dp,0.5_dp)
    m = size(w)
    expected = (0.0_dp,0.0_dp)
    do k = 1, m
        ix = 12+(m-1)/2+1-k
        expected = expected + z(ix)*conjg(w(m+1-k))
    end do
    call require(abs(yc(1,12)-expected) < 1.0e-13_dp, 'complex conjugation')

    yc = cwt_morlet2(z,[-1.0_dp],status=status)
    call require(status == cwt_peak_invalid_input .and. size(yc)==0, 'negative scale rejection')
    yc = cwt_morlet2(z,widths,omega0=nanval,status=status)
    call require(status == cwt_peak_invalid_input .and. size(yc)==0, 'nonfinite frequency rejection')
    z(3) = cmplx(nanval,0.0_dp,dp)
    yc = cwt_complex_wavelet(z,widths,asymmetric_wavelet,status)
    call require(status == cwt_peak_invalid_input .and. size(yc)==0, 'nonfinite signal rejection')
    yc = cwt_morlet2(x, [real(dp)::],status=status)
    call require(status == cwt_peak_invalid_input .and. size(yc)==0, 'empty scales rejection')
    yc = cwt_complex_wavelet(x,widths,wrong_length_wavelet,status)
    call require(status == cwt_peak_invalid_input .and. size(yc)==0, 'short callback rejection')
    yc = cwt_complex_wavelet(x,widths,nan_wavelet,status)
    call require(status == cwt_peak_invalid_input .and. size(yc)==0, 'NaN callback rejection')
    w = morlet2_wavelet(5.0_dp,-2.0_dp)
    call require(size(w)==0,'negative Morlet width')

    print '(a)', 'PASS complex CWT native tests'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! test assertion
        character(*), intent(in) :: message !! diagnostic on assertion failure
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require

    function asymmetric_wavelet(points,width) result(wavelet)
        real(dp), intent(in) :: points !! wavelet length before ceiling
        real(dp), intent(in) :: width !! wavelet scale
        complex(dp), allocatable :: wavelet(:)
        real(dp) :: r
        integer :: t
        allocate(wavelet(ceiling(points)))
        do t = 1, size(wavelet)
            r = (real(t-1,dp)-0.5_dp*(points-1.0_dp))/width
            wavelet(t) = exp(-0.5_dp*r*r)*cmplx(1.0_dp+0.1_dp*r,0.7_dp*sin(r)+0.2_dp,dp)
        end do
    end function asymmetric_wavelet

    function wrong_length_wavelet(points,width) result(wavelet)
        real(dp), intent(in) :: points !! wavelet length before ceiling
        real(dp), intent(in) :: width !! wavelet scale
        complex(dp), allocatable :: wavelet(:)
        allocate(wavelet(max(0,ceiling(points)-1)))
        wavelet = cmplx(width,0.0_dp,dp)
    end function wrong_length_wavelet

    function nan_wavelet(points,width) result(wavelet)
        real(dp), intent(in) :: points !! wavelet length before ceiling
        real(dp), intent(in) :: width !! wavelet scale
        complex(dp), allocatable :: wavelet(:)
        wavelet = asymmetric_wavelet(points,width)
        wavelet(1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    end function nan_wavelet
end program test_signal_cwt_complex
