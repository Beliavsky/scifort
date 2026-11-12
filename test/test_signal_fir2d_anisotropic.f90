! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_fir2d_anisotropic
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_signal_fir2d, only : firwin_2d
    use scifort_signal_design, only : firwin, &
        design_status_success, design_status_invalid_input
    implicit none
    real(dp), allocatable :: taps(:,:), row(:), col(:), expected(:,:), baseline(:,:)
    real(dp) :: bad
    integer :: st, i, j

    row = firwin(9,[0.17_dp],window='hann',status=st)
    if (st /= design_status_success) error stop 'first-axis FIR'
    col = firwin(11,[0.45_dp],window='hamming',fs=3.0_dp,status=st)
    if (st /= design_status_success) error stop 'second-axis FIR'
    taps = firwin_2d([9,11],'hann',0.17_dp,window2='hamming', &
        fc2=[0.45_dp],fs2=3.0_dp,status=st)
    if (st /= design_status_success) error stop 'anisotropic scalar status'
    if (any(shape(taps) /= [9,11])) error stop 'anisotropic scalar shape'
    allocate(expected(9,11))
    do j = 1,11
        expected(:,j) = row*col(j)
    end do
    if (maxval(abs(taps-expected)) > 3.0e-14_dp) error stop 'anisotropic scalar values'

    row = firwin(9,[0.12_dp,0.30_dp],window='blackman',pass_zero='bandpass',status=st)
    col = firwin(11,[0.37_dp],window='kaiser',window_parameter=5.0_dp, &
        pass_zero='highpass',scale=.false.,status=st)
    taps = firwin_2d([9,11],'blackman',[0.12_dp,0.30_dp], &
        window2='kaiser',window_parameter2=5.0_dp, &
        fc2=[0.37_dp],pass_zero='bandpass',pass_zero2='highpass', &
        scale2=.false.,status=st)
    if (st /= design_status_success) error stop 'different passbands status'
    do j = 1,11
        if (maxval(abs(taps(:,j)-row*col(j))) > 3.0e-14_dp) &
            error stop 'different passbands values'
    end do

    baseline = firwin_2d([7,9],'hamming',[0.14_dp,0.31_dp],status=st)
    taps = firwin_2d([7,9],'hamming',[0.14_dp,0.31_dp], &
        fc2=[0.14_dp,0.31_dp],fs2=2.0_dp,status=st)
    if (st /= design_status_success) error stop 'legacy default status'
    if (maxval(abs(taps-baseline)) > 2.0e-15_dp) error stop 'legacy default changed'

    taps = firwin_2d([9,13],'hamming',0.21_dp,circular=.true., &
        ellipse_axes=[1.4_dp,0.6_dp],ellipse_angle=0.5_dp,status=st)
    if (st /= design_status_success) error stop 'elliptical radial status'
    if (any(shape(taps) /= [13,9])) error stop 'elliptical radial shape'
    if (any(taps /= taps)) error stop 'elliptical radial nan'
    if (maxval(abs(taps-taps(13:1:-1,9:1:-1))) > 2.0e-13_dp) &
        error stop 'elliptical radial reversal symmetry'
    ! Angle only rotates the coordinates; a unit-circle profile is invariant.
    baseline = firwin_2d([11,11],'hamming',0.17_dp,circular=.true.,status=st)
    taps = firwin_2d([11,11],'hamming',0.17_dp,circular=.true., &
        ellipse_axes=[1.0_dp,1.0_dp],ellipse_angle=0.53_dp,status=st)
    if (maxval(abs(taps-baseline)) > 4.0e-14_dp) error stop 'unit-circle rotation'
    taps = firwin_2d([9,7],'hamming',0.21_dp,circular=.true., &
        ellipse_axes=[tiny(1.0_dp),tiny(1.0_dp)],status=st)
    if (st /= design_status_success) error stop 'tiny but positive radial scales'
    if (any(taps /= taps)) error stop 'tiny radial scales nan'

    taps = firwin_2d([9,11],'hamming',0.20_dp,fc2=[-0.1_dp],status=st)
    call must_reject(taps,st,'negative secondary cutoff')
    taps = firwin_2d([9,11],'hamming',0.20_dp,fc2=[0.4_dp,0.2_dp],status=st)
    call must_reject(taps,st,'unsorted secondary cutoffs')
    taps = firwin_2d([9,11],'hamming',0.20_dp,fc2=[real(dp)::],status=st)
    call must_reject(taps,st,'empty secondary cutoffs')
    taps = firwin_2d([9,11],'hamming',0.20_dp,fs2=0.0_dp,status=st)
    call must_reject(taps,st,'invalid secondary sampling rate')
    taps = firwin_2d([9,11],'hamming',0.20_dp,pass_zero2='none',status=st)
    call must_reject(taps,st,'invalid secondary passband')
    taps = firwin_2d([9,11],'hamming',0.20_dp,circular=.true.,fc2=[0.3_dp],status=st)
    call must_reject(taps,st,'axis-specific cutoffs in radial design')
    taps = firwin_2d([9,11],'hamming',0.20_dp,circular=.true.,scale2=.false.,status=st)
    call must_reject(taps,st,'axis-specific scaling in radial design')
    taps = firwin_2d([9,11],'hamming',0.20_dp,ellipse_axes=[1.0_dp,0.5_dp],status=st)
    call must_reject(taps,st,'radial axes in separable design')
    taps = firwin_2d([9,11],'hamming',0.20_dp,circular=.true., &
        ellipse_axes=[1.0_dp,0.0_dp],status=st)
    call must_reject(taps,st,'nonpositive ellipse axis')
    bad = ieee_value(0.0_dp,ieee_quiet_nan)
    taps = firwin_2d([9,11],'hamming',0.20_dp,circular=.true., &
        ellipse_axes=[1.0_dp,bad],status=st)
    call must_reject(taps,st,'NaN ellipse axis')
    bad = ieee_value(0.0_dp,ieee_positive_inf)
    taps = firwin_2d([9,11],'hamming',0.20_dp,circular=.true., &
        ellipse_angle=bad,status=st)
    call must_reject(taps,st,'nonfinite ellipse rotation')
    print '(a)', 'PASS: anisotropic and elliptical 2D FIR design'
contains
    subroutine must_reject(filter,code,message)
        real(dp), intent(in) :: filter(:,:) !! rejected FIR result, expected empty
        integer, intent(in) :: code !! expected invalid-input status
        character(*), intent(in) :: message !! case description for failure output
        if (code /= design_status_invalid_input .or. size(filter) /= 0) then
            print '(a)', message
            error stop 'invalid FIR2D input must yield an empty result'
        end if
    end subroutine must_reject
end program test_signal_fir2d_anisotropic
