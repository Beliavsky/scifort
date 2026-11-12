! SPDX-License-Identifier: MIT
! Transfer-function coefficients from cascaded real second-order sections.
program test_signal_sos2tf
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, sos2tf, cheby2, butter, butter_result, &
        design_status_success, design_status_invalid_input, sosfreqz, freqz, &
        frequency_response
    implicit none

    real(dp), parameter :: two_sections(2,6) = reshape([ &
        0.5_dp, 2.0_dp, &
       -1.0_dp, 1.0_dp, &
        0.25_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, &
       -0.5_dp, 0.25_dp, &
        0.2_dp, 0.0_dp], [2,6])
    real(dp) :: not_finite(1,6)
    type(butter_result) :: result, designed
    type(frequency_response) :: tf_response, sos_response
    integer :: order, style
    character(8), parameter :: modes(4) = ['lowpass ', 'highpass', 'bandpass', 'bandstop']

    result = sos2tf(two_sections)
    if (result%status /= design_status_success) error stop 'two section status'
    call check(result%b, [1.0_dp, -1.5_dp, -0.5_dp, 0.25_dp, 0.0_dp], 1.0e-14_dp)
    call check(result%a, [1.0_dp, -0.25_dp, 0.075_dp, 0.05_dp, 0.0_dp], 1.0e-14_dp)
    if (any(result%sos /= two_sections)) error stop 'sos copy differs'

    ! Like SciPy, preserve the original section coefficients, including a0
    ! values other than unity, and retain trailing zeros for first-order rows.
    result = sos2tf(reshape([1.0_dp,2.0_dp,3.0_dp,2.0_dp,5.0_dp,6.0_dp], [1,6]))
    if (result%status /= design_status_success) error stop 'unnormalized section'
    call check(result%b, [1.0_dp,2.0_dp,3.0_dp], 0.0_dp)
    call check(result%a, [2.0_dp,5.0_dp,6.0_dp], 0.0_dp)

    result = sos2tf(reshape([real(dp)::], [0,6]))
    if (result%status /= design_status_success) error stop 'empty cascade status'
    call check(result%b, [1.0_dp], 0.0_dp)
    call check(result%a, [1.0_dp], 0.0_dp)

    result = sos2tf(reshape([real(dp)::], [0,5]))
    if (result%status /= design_status_invalid_input) error stop 'bad column count'
    if (allocated(result%b)) error stop 'bad shape allocated numerator'
    not_finite = 1.0_dp
    not_finite(1,3) = ieee_value(0.0_dp, ieee_quiet_nan)
    result = sos2tf(not_finite)
    if (result%status /= design_status_invalid_input) error stop 'nonfinite input'

    do style = 1, 4
        do order = 1, 6
            if (style <= 2) then
                designed = cheby2(order, 35.0_dp, 0.28_dp, btype=trim(modes(style)))
            else
                designed = butter(order, [0.18_dp, 0.64_dp], btype=trim(modes(style)))
            end if
            if (designed%status /= design_status_success) error stop 'filter design'
            result = sos2tf(designed%sos)
            if (result%status /= design_status_success) error stop 'design cascade conversion'
            ! An odd-order SOS contains a padded first-order section, so
            ! sos2tf retains a trailing zero that the design result omits.
            if (size(result%b) < size(designed%b)) error stop 'polynomial length'
            if (size(result%a) < size(designed%a)) error stop 'denominator length'
            call check(result%b(:size(designed%b)), designed%b, 2.0e-9_dp)
            call check(result%a(:size(designed%a)), designed%a, 2.0e-9_dp)
            if (size(result%b) > size(designed%b)) then
                if (maxval(abs(result%b(size(designed%b)+1:))) > 2.0e-9_dp) &
                    error stop 'numerator trailing coefficients'
            end if
            if (size(result%a) > size(designed%a)) then
                if (maxval(abs(result%a(size(designed%a)+1:))) > 2.0e-9_dp) &
                    error stop 'denominator trailing coefficients'
            end if
            tf_response = freqz(result%b, result%a, w=[0.02_dp,0.18_dp,0.29_dp,0.42_dp])
            sos_response = sosfreqz(designed%sos, w=[0.02_dp,0.18_dp,0.29_dp,0.42_dp])
            if (maxval(abs(tf_response%response-sos_response%response)) > 2.0e-7_dp) &
                error stop 'frequency response mismatch'
        end do
    end do

    print '(a)', 'test_signal_sos2tf passed'
contains
    subroutine check(actual, expected, tolerance)
        real(dp), intent(in) :: actual(:) !! computed coefficients
        real(dp), intent(in) :: expected(:) !! SciPy reference or filter design
        real(dp), intent(in) :: tolerance !! maximum absolute difference allowed
        if (size(actual) /= size(expected)) error stop 'coefficient shape mismatch'
        if (maxval(abs(actual - expected)) > tolerance) &
            error stop 'coefficient value mismatch'
    end subroutine check
end program test_signal_sos2tf
