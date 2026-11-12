! SPDX-License-Identifier: MIT
! Rank-independent padded filtfilt / sosfiltfilt regression tests.
program test_signal_filtfilt_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_filtfilt, only : filtfilt, sosfiltfilt, &
        filtfilt_status_success, filtfilt_status_invalid_input, &
        filtfilt_status_no_steady_state
    use scifort_signal_filtfilt_nd, only : filtfilt_nd, sosfiltfilt_nd
    implicit none
    real(dp) :: b(3), a(3), sos(2,6), bad_sos(2,6)
    complex(dp) :: cb(3), ca(3), csos(2,6)
    real(dp) :: x1(30), y1(30), x5(4,3,24,2,2), y5(4,3,24,2,2)
    real(dp) :: yother(4,3,24,2,2), wrong(4,3,24,2,1)
    real(dp) :: x15(24,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: y15(24,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: xe(0,24), ye(0,24), flat_in(24,2), flat_out(24,2)
    complex(dp) :: c1(30), cy1(30), cx5(4,3,24,2,2), cy5(4,3,24,2,2)
    complex(dp) :: cx15(24,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: cy15(24,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: expected(30)
    real(dp) :: integral_sos(1,6), scalar_input, scalar_output
    integer :: i, code

    b = [0.20_dp,0.12_dp,0.08_dp]
    a = [1.0_dp,-0.3_dp,0.05_dp]
    sos(1,:) = [0.30_dp,0.20_dp,0.10_dp,1.0_dp,-0.27_dp,0.04_dp]
    sos(2,:) = [0.90_dp,-0.12_dp,0.08_dp,1.0_dp,-0.18_dp,0.09_dp]
    cb = cmplx(b,0.0_dp,dp)
    ca = cmplx(a,0.0_dp,dp)
    cb(2) = cb(2)+(0.0_dp,0.06_dp)
    ca(2) = ca(2)-(0.0_dp,0.03_dp)
    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = csos(1,2)+(0.0_dp,0.04_dp)
    csos(2,5) = csos(2,5)-(0.0_dp,0.04_dp)

    x1 = [(sin(0.21_dp*i)+0.02_dp*i,i=1,30)]
    c1 = cmplx(x1,0.2_dp*x1,dp)
    call filtfilt_nd(b,a,x1,y1,status=code)
    if (code /= filtfilt_status_success) error stop 'rank1 filtfilt status'
    if (any(abs(y1-filtfilt(b,a,x1)) > 1.e-13_dp)) error stop 'rank1 filtfilt values'
    call sosfiltfilt_nd(sos,x1,y1,status=code)
    if (code /= filtfilt_status_success) error stop 'rank1 SOS status'
    if (any(abs(y1-sosfiltfilt(sos,x1)) > 1.e-13_dp)) error stop 'rank1 SOS values'
    call filtfilt_nd(cb,ca,c1,cy1,status=code,padtype='even',padlen=3)
    if (code /= filtfilt_status_success) error stop 'rank1 complex status'
    if (any(abs(cy1-filtfilt(cb,ca,c1,padtype='even',padlen=3)) > 1.e-13_dp)) &
        error stop 'rank1 complex values'

    x5 = reshape([(0.2_dp*cos(0.11_dp*i)+0.003_dp*i,i=1,size(x5))],shape(x5))
    cx5 = cmplx(x5,0.1_dp-0.2_dp*x5,dp)
    call filtfilt_nd(b,a,x5,y5,axis=3,padtype='odd',status=code)
    if (code /= filtfilt_status_success) error stop 'rank5 real status'
    if (any(abs(y5(2,1,:,2,1)-filtfilt(b,a,x5(2,1,:,2,1))) > 1.e-13_dp)) &
        error stop 'rank5 real values'
    call filtfilt_nd(b,ca,x5,cy5,axis=3,padtype='constant',padlen=2,status=code)
    if (code /= filtfilt_status_success) error stop 'rank5 rcr status'
    expected = filtfilt(cmplx(b,0.0_dp,dp),ca,cmplx(x5(2,1,:,2,1),0.0_dp,dp), &
        padtype='constant',padlen=2)
    if (any(abs(cy5(2,1,:,2,1)-expected(1:24)) > 1.e-13_dp)) &
        error stop 'rank5 rcr values'
    call filtfilt_nd(cb,a,cx5,cy5,axis=3,padtype='even',padlen=2,status=code)
    if (code /= filtfilt_status_success) error stop 'rank5 crc status'
    expected = filtfilt(cb,cmplx(a,0.0_dp,dp),cx5(2,1,:,2,1),padtype='even',padlen=2)
    if (any(abs(cy5(2,1,:,2,1)-expected(1:24)) > 1.e-13_dp)) &
        error stop 'rank5 crc values'
    call sosfiltfilt_nd(csos,x5,cy5,axis=3,padtype='constant',padlen=2,status=code)
    if (code /= filtfilt_status_success) error stop 'rank5 SOS cr status'
    expected = sosfiltfilt(csos,cmplx(x5(2,1,:,2,1),0.0_dp,dp), &
        padtype='constant',padlen=2)
    if (any(abs(cy5(2,1,:,2,1)-expected(1:24)) > 1.e-13_dp)) &
        error stop 'rank5 SOS cr values'
    call sosfiltfilt_nd(sos,cx5,cy5,axis=3,padtype='none',status=code)
    if (code /= filtfilt_status_success) error stop 'rank5 SOS rc status'
    expected = sosfiltfilt(cmplx(sos,0.0_dp,dp),cx5(2,1,:,2,1),padtype='none')
    if (any(abs(cy5(2,1,:,2,1)-expected(1:24)) > 1.e-13_dp)) &
        error stop 'rank5 SOS rc values'

    call filtfilt_nd(b,a,x5,yother,axis=1,padtype='none',status=code)
    if (code /= filtfilt_status_success) error stop 'axis1 filtfilt status'
    if (any(abs(yother(:,2,3,1,2)-filtfilt(b,a,x5(:,2,3,1,2), &
        padtype='none')) > 1.e-13_dp)) error stop 'axis1 values'
    call sosfiltfilt_nd(sos,x5,yother,axis=2,padtype='none',status=code)
    if (code /= filtfilt_status_success) error stop 'axis2 SOS status'
    if (any(abs(yother(2,:,3,1,2)-sosfiltfilt(sos,x5(2,:,3,1,2), &
        padtype='none')) > 1.e-13_dp)) error stop 'axis2 SOS values'

    ! Noncontiguous actual sections must copy out only into the selected cells.
    y5 = -777.0_dp
    call filtfilt_nd(b,a,x5(1:4:2,:,:,:,:),y5(1:4:2,:,:,:,:), &
        axis=3,padlen=2,status=code)
    if (code /= filtfilt_status_success) error stop 'noncontiguous status'
    if (any(abs(y5(2,:,:,:,:)+777.0_dp) > 0.0_dp)) error stop 'noncontiguous clobber'
    if (any(abs(y5(3,2,:,1,2)-filtfilt(b,a,x5(3,2,:,1,2),padlen=2)) > 1.e-13_dp)) &
        error stop 'noncontiguous values'

    ! Rank 15 works even when the filtering axis is not last.
    x15 = reshape([(0.4_dp*cos(0.06_dp*i),i=1,size(x15))],shape(x15))
    cx15 = cmplx(x15,x15*0.17_dp,dp)
    flat_in = reshape(x15,[24,2])
    call filtfilt_nd(b,a,x15,y15,axis=1,status=code)
    if (code /= filtfilt_status_success) error stop 'rank15 status'
    flat_out = reshape(y15,[24,2])
    if (any(abs(flat_out(:,2)-filtfilt(b,a,flat_in(:,2))) > 1.e-13_dp)) &
        error stop 'rank15 values'
    call sosfiltfilt_nd(csos,cx15,cy15,axis=1,status=code)
    if (code /= filtfilt_status_success) error stop 'rank15 complex SOS status'
    block
        complex(dp) :: flat_cin(24,2), flat_cout(24,2)
        flat_cin = reshape(cx15,[24,2])
        flat_cout = reshape(cy15,[24,2])
        expected(1:24) = sosfiltfilt(csos,flat_cin(:,1))
        if (any(abs(flat_cout(:,1)-expected(1:24)) > 1.e-13_dp)) &
            error stop 'rank15 complex SOS values'
    end block
    call sosfiltfilt_nd(sos,x15,y15,axis=15,padtype='none',status=code)
    if (code /= filtfilt_status_success) error stop 'rank15 axis15 status'
    flat_out = reshape(y15,[24,2])
    if (any(abs(flat_out(1,:)-sosfiltfilt(sos,flat_in(1,:), &
        padtype='none')) > 1.e-13_dp)) error stop 'rank15 axis15 values'

    ! Empty outer lanes are valid if the filtered extent/padding are valid.
    call filtfilt_nd(b,a,xe,ye,axis=2,padtype='none',status=code)
    if (code /= filtfilt_status_success) error stop 'empty channels status'
    call sosfiltfilt_nd(sos,xe,ye,axis=2,padtype='none',status=code)
    if (code /= filtfilt_status_success) error stop 'empty SOS status'
    call filtfilt_nd(b,a,x5,wrong,axis=3,status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'mismatched shape status'
    if (.not. all(ieee_is_nan(wrong))) error stop 'mismatched shape nan'
    call filtfilt_nd(b,a,x5,y5,axis=6,status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'invalid axis'
    if (.not. all(ieee_is_nan(y5))) error stop 'invalid axis nan'
    call filtfilt_nd(b,a,x5,y5,axis=1,status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'short axis default pad'
    call sosfiltfilt_nd(sos,x5,y5,axis=1,padtype='even',padlen=4,status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'short SOS pad'
    call sosfiltfilt_nd(sos,x5,y5,axis=3,padtype='invalid',status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'invalid padding mode'
    call filtfilt_nd(b,a,x5,y5,axis=3,padlen=-1,status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'negative padding length'
    bad_sos = sos
    bad_sos(1,4) = 2.0_dp
    call sosfiltfilt_nd(bad_sos,x5,y5,axis=3,padlen=1,status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'invalid SOS normalization'

    scalar_input = 1.0_dp
    call filtfilt_nd(b,a,scalar_input,scalar_output,status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'scalar input invalid'
    if (.not. ieee_is_nan(scalar_output)) error stop 'scalar output nan'
    call sosfiltfilt_nd(sos,scalar_input,scalar_output,status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'scalar SOS invalid'
    if (.not. ieee_is_nan(scalar_output)) error stop 'scalar SOS nan'
    call filtfilt_nd([1.0_dp],[1.0_dp,-1.0_dp],x1,y1, &
        padtype='none',status=code)
    if (code /= filtfilt_status_no_steady_state) error stop 'IIR pole at one'
    integral_sos(1,:) = [1.0_dp,0.0_dp,0.0_dp,1.0_dp,-1.0_dp,0.0_dp]
    call sosfiltfilt_nd(integral_sos,x1,y1,padtype='none',status=code)
    if (code /= filtfilt_status_no_steady_state) error stop 'SOS pole at one'
    print '(a)', 'test_signal_filtfilt_nd passed'
end program test_signal_filtfilt_nd
