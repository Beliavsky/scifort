! SPDX-License-Identifier: MIT
! Assumed-rank IIR filtering, state propagation, and mixed types through rank 15.
program test_signal_lfilter_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal, only : dp, lfilter, lfilter_nd, &
        lfilter_status_success, lfilter_status_invalid_input
    implicit none

    real(dp), parameter :: b(3) = [0.18_dp, 0.27_dp, 0.12_dp]
    real(dp), parameter :: a(2) = [2.0_dp, -0.8_dp]
    complex(dp) :: cb(3), ca(2)
    real(dp) :: x5(7,4,3,2,5), y5(7,4,3,2,5), zi5(7,2,3,2,5), zf5(7,2,3,2,5)
    real(dp) :: split5(7,4,3,2,5), split_state(7,2,3,2,5), split_end(7,2,3,2,5)
    real(dp) :: bad5(7,4,3,2,4), badzi(7,1,3,2,5), short_state(7,1,3,2,5)
    real(dp) :: x15(5,2,3,2,4,1,1,1,1,1,1,1,1,1,6)
    real(dp) :: y15(5,2,3,2,4,1,1,1,1,1,1,1,1,1,6)
    real(dp) :: zi15(5,2,3,2,4,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: zf15(5,2,3,2,4,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: empty(0,3), empty_out(0,3), empty_state(2,3), empty_final(2,3)
    real(dp) :: one(8), one_y(8), no_state(0), no_final_state(0), one_b(1)
    complex(dp) :: cx5(7,4,3,2,5), cy5(7,4,3,2,5), czi5(7,2,3,2,5), czf5(7,2,3,2,5)
    complex(dp) :: cx15(5,2,3,2,4,1,1,1,1,1,1,1,1,1,6)
    complex(dp) :: cy15(5,2,3,2,4,1,1,1,1,1,1,1,1,1,6)
    complex(dp) :: czi15(5,2,3,2,4,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: czf15(5,2,3,2,4,1,1,1,1,1,1,1,1,1,2)
    real(dp), allocatable :: reference(:), final_state(:)
    complex(dp), allocatable :: creference(:), cfinal_state(:)
    integer :: i, code

    x5 = reshape([(sin(0.017_dp*i)+0.003_dp*i, i=1,size(x5))], shape(x5))
    zi5 = reshape([(0.02_dp*cos(0.19_dp*i), i=1,size(zi5))],shape(zi5))
    cx5 = cmplx(x5,0.4_dp-0.25_dp*x5,dp)
    czi5 = cmplx(zi5,-0.8_dp*zi5,dp)
    cb = cmplx(b, 0.15_dp*b, dp)
    ca = cmplx(a, -0.1_dp*a, dp)

    call lfilter_nd(b,a,x5,y5,zi=zi5,zf=zf5,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'rank5 real status'
    reference = lfilter(b,a,x5(3,:,2,1,4),zi=zi5(3,:,2,1,4),zf=final_state)
    if (maxval(abs(reference-y5(3,:,2,1,4))) > 1.e-13_dp) error stop 'rank5 real result'
    if (maxval(abs(final_state-zf5(3,:,2,1,4))) > 1.e-13_dp) error stop 'rank5 real state'
    call lfilter_nd(b,a,x5(:,1:2,:,:,:),split5(:,1:2,:,:,:), &
        zi=zi5,zf=split_state,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'first streaming block'
    call lfilter_nd(b,a,x5(:,3:4,:,:,:),split5(:,3:4,:,:,:), &
        zi=split_state,zf=split_end,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'second streaming block'
    if (maxval(abs(split5-y5)) > 1.e-13_dp) error stop 'streaming output'
    if (maxval(abs(split_end-zf5)) > 1.e-13_dp) error stop 'streaming final state'

    call lfilter_nd(cb,ca,cx5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'rank5 complex status'
    creference = lfilter(cb,ca,cx5(3,:,2,1,4),zi=czi5(3,:,2,1,4),zf=cfinal_state)
    if (maxval(abs(creference-cy5(3,:,2,1,4))) > 1.e-13_dp) error stop 'rank5 complex result'
    if (maxval(abs(cfinal_state-czf5(3,:,2,1,4))) > 1.e-13_dp) error stop 'rank5 complex state'

    ! Exercise all six mixed coefficient/input combinations at rank 5.
    call lfilter_nd(b,a,cx5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'rrc status'
    creference = lfilter(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp), &
        cx5(3,:,2,1,4),zi=czi5(3,:,2,1,4),zf=cfinal_state)
    if (maxval(abs(creference-cy5(3,:,2,1,4))) > 1.e-13_dp) error stop 'rrc result'
    call lfilter_nd(b,ca,x5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'rcr status'
    creference = lfilter(cmplx(b,0.0_dp,dp),ca,cmplx(x5(3,:,2,1,4),0.0_dp,dp), &
        zi=czi5(3,:,2,1,4),zf=cfinal_state)
    if (maxval(abs(creference-cy5(3,:,2,1,4))) > 1.e-13_dp) error stop 'rcr result'
    call lfilter_nd(b,ca,cx5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'rcc status'
    call lfilter_nd(cb,a,x5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'crr status'
    call lfilter_nd(cb,a,cx5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'crc status'
    call lfilter_nd(cb,ca,x5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'ccr status'
    creference = lfilter(cb,ca,cmplx(x5(3,:,2,1,4),0.0_dp,dp), &
        zi=czi5(3,:,2,1,4),zf=cfinal_state)
    if (maxval(abs(creference-cy5(3,:,2,1,4))) > 1.e-13_dp) error stop 'ccr result'
    if (maxval(abs(cfinal_state-czf5(3,:,2,1,4))) > 1.e-13_dp) error stop 'ccr state'

    x15 = reshape([(0.4_dp*sin(0.061_dp*i)+0.01_dp*i, i=1,size(x15))],shape(x15))
    zi15 = 0.04_dp
    cx15 = cmplx(x15, 0.22_dp-0.31_dp*x15, dp)
    czi15 = cmplx(zi15,-0.02_dp,dp)
    call lfilter_nd(b,a,x15,y15,zi=zi15,zf=zf15,axis=15,status=code)
        if (code /= lfilter_status_success) error stop 'rank15 real status'
        reference = lfilter(b,a,x15(1,1,2,1,3,1,1,1,1,1,1,1,1,1,:), &
            zi=zi15(1,1,2,1,3,1,1,1,1,1,1,1,1,1,:),zf=final_state)
        if (maxval(abs(reference-y15(1,1,2,1,3,1,1,1,1,1,1,1,1,1,:))) > 1.e-13_dp) &
            error stop 'rank15 real result'
        if (maxval(abs(final_state-zf15(1,1,2,1,3,1,1,1,1,1,1,1,1,1,:))) > 1.e-13_dp) &
            error stop 'rank15 real final state'
        call lfilter_nd(cb,ca,cx15,cy15,zi=czi15,zf=czf15,axis=15,status=code)
        if (code /= lfilter_status_success) error stop 'rank15 complex status'
        creference = lfilter(cb,ca,cx15(1,1,2,1,3,1,1,1,1,1,1,1,1,1,:), &
            zi=czi15(1,1,2,1,3,1,1,1,1,1,1,1,1,1,:),zf=cfinal_state)
        if (maxval(abs(creference-cy15(1,1,2,1,3,1,1,1,1,1,1,1,1,1,:))) > 1.e-13_dp) &
            error stop 'rank15 complex result'
        call lfilter_nd(cb,ca,x15,cy15,zi=czi15,zf=czf15,axis=15,status=code)
        if (code /= lfilter_status_success) error stop 'rank15 mixed status'
    call lfilter_nd(b,a,x15,y15,axis=1,status=code)
    if (code /= lfilter_status_success) error stop 'rank15 first axis status'
    reference = lfilter(b,a,x15(:,1,2,1,3,1,1,1,1,1,1,1,1,1,4))
    if (maxval(abs(reference-y15(:,1,2,1,3,1,1,1,1,1,1,1,1,1,4))) > 1.e-13_dp) &
        error stop 'rank15 first axis result'
    call lfilter_nd(b,a,x15,y15,axis=3,status=code)
    if (code /= lfilter_status_success) error stop 'rank15 middle axis status'
    reference = lfilter(b,a,x15(1,1,:,1,3,1,1,1,1,1,1,1,1,1,4))
    if (maxval(abs(reference-y15(1,1,:,1,3,1,1,1,1,1,1,1,1,1,4))) > 1.e-13_dp) &
        error stop 'rank15 middle axis result'

    ! Confirm noncontiguous input, output, and state sections copy correctly.
    y5 = -999.0_dp
    zf5 = -888.0_dp
    call lfilter_nd(b,a,x5(1:7:2,:,:,:,:),y5(1:7:2,:,:,:,:), &
        zi=zi5(1:7:2,:,:,:,:),zf=zf5(1:7:2,:,:,:,:),axis=2,status=code)
    if (code /= lfilter_status_success) error stop 'strided real status'
    reference = lfilter(b,a,x5(3,:,2,1,4),zi=zi5(3,:,2,1,4),zf=final_state)
    if (maxval(abs(reference-y5(3,:,2,1,4))) > 1.e-13_dp) error stop 'strided real result'
    if (any(y5(2,:,:,:,:) /= -999.0_dp)) error stop 'strided real output outside section'
    if (any(zf5(2,:,:,:,:) /= -888.0_dp)) error stop 'strided real state outside section'

    ! No filter state and zero-length sample axes must be well-defined.
    one = [(0.1_dp*i,i=1,size(one))]
    one_b = [2.0_dp]
    call lfilter_nd(one_b,[1.0_dp],one,one_y,zi=no_state,zf=no_final_state,status=code)
    if (code /= lfilter_status_success) error stop 'zero-state status'
    if (maxval(abs(one_y-2.0_dp*one)) > 1.e-14_dp) error stop 'zero-state output'
    empty_state = 0.6_dp
    call lfilter_nd(b,a,empty,empty_out,zi=empty_state,zf=empty_final,axis=1,status=code)
    if (code /= lfilter_status_success) error stop 'empty axis status'
    if (any(empty_final /= empty_state)) error stop 'empty axis final state'

    call lfilter_nd(b,a,x5,bad5,axis=2,status=code)
    if (code /= lfilter_status_invalid_input) error stop 'invalid output shape status'
    if (.not. all(ieee_is_nan(bad5))) error stop 'invalid output shape NaNs'
    call lfilter_nd(b,a,x5,y5,zi=badzi,axis=2,status=code)
    if (code /= lfilter_status_invalid_input) error stop 'invalid initial shape'
    call lfilter_nd(b,a,x5,y5,zf=short_state,axis=2,status=code)
    if (code /= lfilter_status_invalid_input) error stop 'invalid final shape'
    call lfilter_nd(b,a,x15,y15,axis=16,status=code)
    if (code /= lfilter_status_invalid_input) error stop 'invalid axis'
    call lfilter_nd(b,[0.0_dp],x5,y5,status=code)
    if (code /= lfilter_status_invalid_input) error stop 'invalid denominator'

    print '(a)', 'test_signal_lfilter_nd passed'
end program test_signal_lfilter_nd
