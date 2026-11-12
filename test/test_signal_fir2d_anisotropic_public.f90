! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_fir2d_anisotropic_public
    use scifort_signal, only : dp, firwin_2d, firwin, &
        design_status_success, design_status_invalid_input
    implicit none
    real(dp), allocatable :: h(:,:), a(:), b(:)
    integer :: st, j

    h = firwin_2d([7,13],'hamming',0.19_dp,window2='hann', &
        fc2=[0.16_dp,0.48_dp],pass_zero2='bandpass',status=st)
    if (st /= design_status_success) error stop 'public anisotropic status'
    a = firwin(7,[0.19_dp],window='hamming',status=st)
    b = firwin(13,[0.16_dp,0.48_dp],window='hann',pass_zero='bandpass',status=st)
    if (any(shape(h) /= [7,13])) error stop 'public anisotropic shape'
    do j = 1,13
        if (maxval(abs(h(:,j)-a*b(j))) > 4.0e-14_dp) &
            error stop 'public anisotropic values'
    end do
    h = firwin_2d([11,15],'hamming',0.21_dp,circular=.true., &
        ellipse_axes=[0.65_dp,1.4_dp],ellipse_angle=0.4_dp,status=st)
    if (st /= design_status_success) error stop 'public elliptical status'
    if (any(shape(h) /= [15,11])) error stop 'public elliptical shape'
    h = firwin_2d([9,11],'hamming',0.21_dp,circular=.true., &
        fs2=2.0_dp,status=st)
    if (st /= design_status_invalid_input .or. size(h) /= 0) &
        error stop 'public invalid combinations'
    print '(a)', 'PASS: anisotropic FIR2D via public scifort_signal'
end program test_signal_fir2d_anisotropic_public
