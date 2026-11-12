! SPDX-License-Identifier: MIT
! Regression checks for Bessel/Thomson analog and bilinear digital IIR filters.
program test_signal_bessel_design
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : bessel, iirfilter, butter_result, design_status_success, &
        design_status_invalid_input
    implicit none

    type(butter_result) :: filter, selected
    complex(dp) :: z
    real(dp) :: lowpass_gain, highpass_gain
    integer :: n, i
    character(8), parameter :: modes(4) = [character(8) :: &
        'lowpass', 'highpass', 'bandpass', 'bandstop']
    real(dp), parameter :: single(1) = [0.23_dp]
    real(dp), parameter :: band(2) = [0.18_dp,0.37_dp]

    do n = 1, 20
        do i = 1, 4
            if (i <= 2) then
                filter = bessel(n,single,btype=trim(modes(i)),fs=2.0_dp)
            else
                filter = bessel(n,band,btype=trim(modes(i)),fs=2.0_dp)
            end if
            if (filter%status /= design_status_success) then
                print *, 'Bessel design failed', n, trim(modes(i))
                error stop 'Bessel filter design'
            end if
            if (size(filter%a) /= size(filter%b)) error stop 'digital polynomial sizes'
            if (size(filter%a) /= 1+n*merge(2,1,i>2)) error stop 'digital degree'
            if (size(filter%sos,2) /= 6) error stop 'digital SOS size'
            if (size(filter%sos,1) /= (size(filter%a))/2) then
                ! Odd orders round upward to form a complete last SOS row.
                if (size(filter%sos,1) /= (size(filter%a)+1)/2) error stop 'SOS count'
            end if
            if (any(abs(filter%sos(:,4)-1.0_dp) > 1.0e-10_dp)) &
                error stop 'monic digital SOS denominator'
            if (i == 1 .or. i == 4) then
                z = cmplx(1.0_dp,0.0_dp,dp)
                lowpass_gain = abs(sos_response(filter%sos,z))
                if (abs(lowpass_gain-1.0_dp) > 5.0e-8_dp) error stop 'DC gain'
            else if (i == 2) then
                z = cmplx(-1.0_dp,0.0_dp,dp)
                highpass_gain = abs(sos_response(filter%sos,z))
                if (abs(highpass_gain-1.0_dp) > 5.0e-8_dp) error stop 'Nyquist gain'
            end if
        end do
    end do

    do n = 1, 20
        do i = 1, 4
            if (i <= 2) then
                filter = bessel(n,[3.5_dp],btype=trim(modes(i)),analog=.true.)
            else
                filter = bessel(n,[2.0_dp,7.0_dp],btype=trim(modes(i)),analog=.true.)
            end if
            if (filter%status /= design_status_success) then
                print *, 'Bessel analog design failed', n, trim(modes(i))
                error stop 'Bessel analog filter design'
            end if
            if (allocated(filter%sos)) error stop 'analog SOS must remain unallocated'
            if (size(filter%a) /= 1+n*merge(2,1,i>2)) error stop 'analog degree'
        end do
    end do

    ! All three SciPy normalization conventions, including high-order designs.
    do i = 1, 3
        select case (i)
        case (1)
            filter = bessel(17,0.25_dp,norm='phase')
        case (2)
            filter = bessel(17,0.25_dp,norm='delay')
        case (3)
            filter = bessel(17,0.25_dp,norm='mag')
        end select
        if (filter%status /= design_status_success) error stop 'Bessel normalization'
    end do

    filter = bessel(6,0.25_dp,btype='highpass')
    selected = iirfilter(6,[0.25_dp],btype='highpass',ftype='bessel')
    if (selected%status /= design_status_success) error stop 'iirfilter Bessel'
    if (maxval(abs(selected%a-filter%a)) > 5.0e-13_dp) error stop 'iirfilter a'
    if (maxval(abs(selected%b-filter%b)) > 5.0e-13_dp) error stop 'iirfilter b'
    filter = bessel(6,0.25_dp,btype='highpass',norm='delay',analog=.true.)
    if (filter%status /= design_status_success) error stop 'analog delay highpass'
    filter = bessel(21,0.25_dp)
    if (filter%status /= design_status_invalid_input) error stop 'order limit'
    filter = bessel(-1,0.25_dp)
    if (filter%status /= design_status_invalid_input) error stop 'negative order'
    filter = bessel(5,0.25_dp,norm='unknown')
    if (filter%status /= design_status_invalid_input) error stop 'unknown norm'
    filter = bessel(5,0.25_dp,analog=.true.,fs=20.0_dp)
    if (filter%status /= design_status_invalid_input) error stop 'analog with fs'
    filter = bessel(5,[0.2_dp,0.1_dp],btype='bandpass')
    if (filter%status /= design_status_invalid_input) error stop 'reversed band'
    filter = bessel(5,0.0_dp)
    if (filter%status /= design_status_invalid_input) error stop 'zero frequency'
    filter = bessel(5,1.0_dp)
    if (filter%status /= design_status_invalid_input) error stop 'Nyquist frequency'
    filter = bessel(5,ieee_value(0.0_dp,ieee_quiet_nan))
    if (filter%status /= design_status_invalid_input) error stop 'NaN frequency'
    filter = bessel(5,[0.15_dp,0.28_dp],btype='lowpass')
    if (filter%status /= design_status_invalid_input) error stop 'bad cutoff rank'
    selected = iirfilter(5,[0.25_dp],ftype='bessel',rp=0.5_dp)
    if (selected%status /= design_status_invalid_input) error stop 'unexpected ripple'
    print '(a)', 'test_signal_bessel_design passed'

contains

    pure function sos_response(sos, z) result(response)
        real(dp), intent(in) :: sos(:,:) !! digital sections, six coefficients per row
        complex(dp), intent(in) :: z !! point on the unit circle
        complex(dp) :: response
        complex(dp) :: zi
        integer :: j

        zi = 1.0_dp/z
        response = cmplx(1.0_dp,0.0_dp,dp)
        do j = 1, size(sos,1)
            response = response * (sos(j,1)+zi*(sos(j,2)+zi*sos(j,3))) / &
                (sos(j,4)+zi*(sos(j,5)+zi*sos(j,6)))
        end do
    end function sos_response
end program test_signal_bessel_design
