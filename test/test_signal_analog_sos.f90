! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_analog_sos
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal, only : analog_zpk_result, analog_sos_result, &
        analog_output_status_success, analog_output_status_invalid_input, &
        iirfilter_zpk, iirfilter_sos, zpk2sos_analog, sosfreqs, freqs_zpk, &
        frequency_response, design_status_success, design_status_invalid_input, &
        design_status_singular
    implicit none

    character(7), parameter :: families(4) = ['butter ', 'cheby1 ', 'cheby2 ', 'ellip  ']
    character(8), parameter :: shapes(4) = ['lowpass ', 'highpass', 'bandpass', 'bandstop']
    real(dp), parameter :: w(10) = [0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp, 1.5_dp, &
        2.0_dp, 2.8_dp, 4.5_dp, 7.0_dp, 12.0_dp]
    real(dp) :: nan
    complex(dp), allocatable :: z(:), p(:)
    type(analog_zpk_result) :: design
    type(analog_sos_result) :: sections
    type(frequency_response) :: h_sos, h_zpk
    integer :: i, j, n

    do i = 1, size(families)
        do j = 1, size(shapes)
            do n = 1, 8
                if (j < 3) then
                    design = iirfilter_zpk(n, 3.0_dp, ftype=trim(families(i)), &
                        btype=trim(shapes(j)), rp=1.0_dp, rs=45.0_dp)
                    sections = iirfilter_sos(n, 3.0_dp, ftype=trim(families(i)), &
                        btype=trim(shapes(j)), rp=1.0_dp, rs=45.0_dp)
                else
                    design = iirfilter_zpk(n, [2.0_dp,5.0_dp], ftype=trim(families(i)), &
                        btype=trim(shapes(j)), rp=1.0_dp, rs=45.0_dp)
                    sections = iirfilter_sos(n, [2.0_dp,5.0_dp], ftype=trim(families(i)), &
                        btype=trim(shapes(j)), rp=1.0_dp, rs=45.0_dp)
                end if
                if (design%status /= analog_output_status_success) error stop 'analog ZPK design failed'
                if (sections%status /= analog_output_status_success) error stop 'analog SOS design failed'
                if (size(design%poles) /= n*merge(1,2,j<3)) error stop 'wrong pole count'
                if (size(design%zeros) > size(design%poles)) error stop 'improper filter'
                if (size(sections%sos,2) /= 6) error stop 'SOS columns'
                if (size(sections%sos,1) /= (size(design%poles)+1)/2) error stop 'SOS rows'
                h_sos = sosfreqs(sections%sos,w)
                h_zpk = freqs_zpk(design%zeros,design%poles,design%gain,w=w)
                if (h_sos%status /= design_status_success) error stop 'SOS frequency response'
                if (h_zpk%status /= design_status_success) error stop 'ZPK frequency response'
                if (any(abs(h_sos%response-h_zpk%response) > &
                    2.0e-9_dp*max(1.0_dp,abs(h_zpk%response)))) error stop 'SOS/ZPK mismatch'
                sections = zpk2sos_analog(design%zeros,design%poles,design%gain)
                if (sections%status /= analog_output_status_success) error stop 'roundtrip failed'
            end do
        end do
    end do

    ! Bessel/Thomson designs need a ZPK path at higher orders because
    ! expanded transfer-function polynomials become ill-conditioned.
    do j = 1, size(shapes)
        do n = 1, 20
            if (j < 3) then
                design = iirfilter_zpk(n,3.0_dp,ftype='bessel',btype=trim(shapes(j)))
                sections = iirfilter_sos(n,3.0_dp,ftype='thomson',btype=trim(shapes(j)))
            else
                design = iirfilter_zpk(n,[2.0_dp,5.0_dp],ftype='bessel',btype=trim(shapes(j)))
                sections = iirfilter_sos(n,[2.0_dp,5.0_dp],ftype='thomson',btype=trim(shapes(j)))
            end if
            if (design%status /= analog_output_status_success) error stop 'Bessel analog ZPK'
            if (sections%status /= analog_output_status_success) error stop 'Bessel analog SOS'
            if (size(design%poles) /= n*merge(1,2,j<3)) error stop 'Bessel analog roots'
            h_sos = sosfreqs(sections%sos,w)
            h_zpk = freqs_zpk(design%zeros,design%poles,design%gain,w=w)
            if (h_sos%status /= design_status_success) error stop 'Bessel SOS response'
            if (h_zpk%status /= design_status_success) error stop 'Bessel ZPK response'
            if (any(abs(h_sos%response-h_zpk%response) > &
                1.0e-7_dp*max(1.0_dp,abs(h_zpk%response)))) error stop 'Bessel SOS/ZPK mismatch'
        end do
    end do
    design = iirfilter_zpk(5,3.0_dp,ftype='bessel',rp=1.0_dp)
    if (design%status /= analog_output_status_invalid_input) error stop 'Bessel extra ripple'
    design = iirfilter_zpk(21,3.0_dp,ftype='bessel')
    if (design%status /= analog_output_status_invalid_input) error stop 'Bessel order limit'

    ! A very high-frequency, high-order bandpass has finite ZPK/SOS data
    ! even though its expanded degree-40 BA polynomials overflow real(dp).
    design = iirfilter_zpk(20, [1.0e10_dp,1.0e11_dp], btype='bandpass')
    if (design%status /= analog_output_status_success) error stop 'high-order roots'
    sections = iirfilter_sos(20, [1.0e10_dp,1.0e11_dp], btype='bandpass')
    if (sections%status /= analog_output_status_success) error stop 'high-order sections'
    if (.not. all(ieee_is_finite(sections%sos))) error stop 'nonfinite high-order SOS'

    ! A polynomial with real roots must support more than one real root.
    z = [cmplx(-0.7_dp,0.0_dp,dp)]
    p = [cmplx(-1.0_dp,0.0_dp,dp), cmplx(-2.0_dp,0.0_dp,dp), cmplx(-3.0_dp,0.0_dp,dp)]
    sections = zpk2sos_analog(z,p,4.0_dp)
    if (sections%status /= analog_output_status_success) error stop 'real roots'
    h_sos = sosfreqs(sections%sos,w)
    h_zpk = freqs_zpk(z,p,4.0_dp,w=w)
    if (any(abs(h_sos%response-h_zpk%response) > 1.0e-12_dp)) error stop 'real roots response'

    ! Degenerate constant transfer function: no poles and no zeros.
    sections = zpk2sos_analog([complex(dp)::], [complex(dp)::], 2.5_dp)
    if (sections%status /= analog_output_status_success) error stop 'constant transfer function'
    h_sos = sosfreqs(sections%sos,w)
    if (any(abs(h_sos%response-cmplx(2.5_dp,0.0_dp,dp)) > 1.0e-12_dp)) &
        error stop 'constant response'

    sections = zpk2sos_analog([cmplx(1.0_dp,0.0_dp,dp)], [complex(dp)::], 1.0_dp)
    if (sections%status /= analog_output_status_invalid_input) error stop 'improper filter accepted'
    sections = zpk2sos_analog([cmplx(0.0_dp,1.0_dp,dp)], &
        [cmplx(-1.0_dp,0.0_dp,dp)], 1.0_dp)
    if (sections%status /= analog_output_status_invalid_input) error stop 'unpaired complex zero'
    design = iirfilter_zpk(3, 3.0_dp, ftype='unsupported')
    if (design%status /= analog_output_status_invalid_input) error stop 'unknown family'
    design = iirfilter_zpk(3, 3.0_dp, ftype='cheby1')
    if (design%status /= analog_output_status_invalid_input) error stop 'missing ripple'
    design = iirfilter_zpk(3, [3.0_dp,2.0_dp], btype='bandpass')
    if (design%status /= analog_output_status_invalid_input) error stop 'reversed cutoff'
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    design = iirfilter_zpk(3, nan)
    if (design%status /= analog_output_status_invalid_input) error stop 'NaN cutoff'
    sections = zpk2sos_analog([complex(dp)::], [cmplx(-1.0_dp,0.0_dp,dp)], nan)
    if (sections%status /= analog_output_status_invalid_input) error stop 'NaN gain'
    h_sos = sosfreqs(reshape([0.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp], [1,6]), [0.0_dp])
    if (h_sos%status /= design_status_singular) error stop 'singular response status'
    h_sos = sosfreqs(reshape([0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp], [1,6]),w)
    if (h_sos%status /= design_status_invalid_input) error stop 'zero denominator'

    print '(a)', 'PASS test_signal_analog_sos'
end program test_signal_analog_sos
