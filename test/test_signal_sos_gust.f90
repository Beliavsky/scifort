! SPDX-License-Identifier: MIT
! SOS-native Gustafsson initialization: scalar, mixed types and ranks 1-15.
program test_signal_sos_gust
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : sosfiltfilt, sosfiltfilt_nd, sosfiltfilt_gust, &
        filtfilt, filtfilt_status_success, filtfilt_status_invalid_input, &
        filtfilt_status_numerical_failure
    implicit none
    real(dp) :: sos(2,6), one(1,6), high(6,6), x(45), x2(5,22), x3(2,8,3)
    real(dp) :: x5(2,1,11,1,2), y5(2,1,11,1,2), xr(2), yr(2)
    real(dp) :: x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,13)
    real(dp) :: y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,13)
    real(dp) :: xe(0,12), ye(0,12), nanval
    real(dp), allocatable :: y(:), ref(:), ypad(:), ymat(:,:), ycube(:,:,:), xlong(:)
    complex(dp) :: csos(2,6), cx(45), cx5(2,1,11,1,2), cy5(2,1,11,1,2)
    complex(dp), allocatable :: cy(:), cref(:)
    integer :: i, status
    real(dp), parameter :: tol = 2.e-10_dp

    one(1,:) = [0.1_dp,0.2_dp,0.1_dp,1.0_dp,-0.35_dp,0.08_dp]
    sos(1,:) = [0.32_dp,0.12_dp,0.08_dp,1.0_dp,-0.4_dp,0.13_dp]
    sos(2,:) = [0.7_dp,0.2_dp,0.1_dp,1.0_dp,-0.25_dp,0.09_dp]
    ! Fixed coefficients from scipy.signal.butter(12,0.16,output='sos').
    high(1,:) = [1.16468612940542510e-08_dp, &
        2.32937225881085021e-08_dp, &
        1.16468612940542510e-08_dp, &
        1.00000000000000000e+00_dp, &
        -1.18609580556573513e+00_dp, &
        3.53516791069497571e-01_dp]
    high(2,:) = [1.00000000000000000e+00_dp, &
        2.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        -1.21281209262021861e+00_dp, &
        3.84004162286553707e-01_dp]
    high(3,:) = [1.00000000000000000e+00_dp, &
        2.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        -1.26798743738845299e+00_dp, &
        4.46967672693062035e-01_dp]
    high(4,:) = [1.00000000000000000e+00_dp, &
        2.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        -1.35517658332736501e+00_dp, &
        5.46463828461893786e-01_dp]
    high(5,:) = [1.00000000000000000e+00_dp, &
        2.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        -1.47979889439721690e+00_dp, &
        6.88676953053861785e-01_dp]
    high(6,:) = [1.00000000000000000e+00_dp, &
        2.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        1.00000000000000000e+00_dp, &
        -1.64892643738303146e+00_dp, &
        8.81677356722299965e-01_dp]
    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = csos(1,2) + (0.0_dp,0.02_dp)
    csos(2,5) = csos(2,5) - (0.0_dp,0.04_dp)
    x = [(sin(0.18_dp*i)+0.3_dp*cos(0.072_dp*i),i=1,size(x))]
    cx = cmplx(x,0.15_dp*x+0.1_dp,dp)
    nanval = ieee_value(0.0_dp,ieee_quiet_nan)

    y = sosfiltfilt(one,x,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'sos rank1 gust status'
    ref = filtfilt(one(1,1:3),one(1,4:6),x,method='gust')
    if (maxval(abs(y-ref)) > tol) error stop 'SOS Gust differs from first-order cascade'
    ref = sosfiltfilt_gust(one,x,status=status)
    if (status /= filtfilt_status_success .or. maxval(abs(y-ref)) > tol) &
        error stop 'direct SOS Gust differs'
    y = sosfiltfilt(sos,x,method='gust',irlen=8,padtype='ignored',padlen=-8,status=status)
    if (status /= filtfilt_status_success .or. .not. all(abs(y) < 1.0e3_dp)) &
        error stop 'irlen filtering status'
    ref = sosfiltfilt_gust(sos,x,irlen=8)
    if (maxval(abs(y-ref)) > tol) error stop 'direct irlen comparison'
    ypad = sosfiltfilt(sos,x,padlen=7,status=status)
    if (status /= filtfilt_status_success) error stop 'legacy pad status'
    ref = sosfiltfilt(sos,x,padlen=7,method='pad')
    if (maxval(abs(ypad-ref)) > tol) error stop 'legacy pad changed'

    cy = sosfiltfilt(csos,cx,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'complex SOS Gust status'
    cref = sosfiltfilt_gust(csos,cx)
    if (maxval(abs(cy-cref)) > tol) error stop 'complex SOS direct differs'
    cy = sosfiltfilt_gust(sos,cx,irlen=4,status=status)
    if (status /= filtfilt_status_success) error stop 'mixed real SOS complex input status'
    cref = sosfiltfilt_gust(cmplx(sos,0.0_dp,dp),cx,irlen=4)
    if (maxval(abs(cy-cref)) > tol) error stop 'mixed real SOS complex data'
    cy = sosfiltfilt_gust(csos,x,irlen=4,status=status)
    if (status /= filtfilt_status_success) error stop 'mixed complex SOS real input status'
    cref = sosfiltfilt_gust(csos,cmplx(x,0.0_dp,dp),irlen=4)
    if (maxval(abs(cy-cref)) > tol) error stop 'mixed complex SOS real data'

    x2 = reshape([(sin(0.053_dp*i),i=1,size(x2))],shape(x2))
    x3 = reshape([(cos(0.037_dp*i),i=1,size(x3))],shape(x3))
    ymat = sosfiltfilt(sos,x2,axis=1,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'rank2 status'
    ref = sosfiltfilt_gust(sos,x2(:,3))
    if (any(abs(ymat(:,3)-ref)>tol)) error stop 'rank2 values'
    ycube = sosfiltfilt(sos,x3,axis=2,method='gust',irlen=3,status=status)
    if (status /= filtfilt_status_success) error stop 'rank3 status'
    ref = sosfiltfilt_gust(sos,x3(1,:,2),irlen=3)
    if (any(abs(ycube(1,:,2)-ref)>tol)) error stop 'rank3 values'

    x5 = reshape([(sin(0.061_dp*i)+0.001_dp*i,i=1,size(x5))],shape(x5))
    cx5 = cmplx(x5,0.2_dp*x5,dp)
    call sosfiltfilt_nd(sos,x5,y5,axis=3,method='gust',irlen=4, &
        padtype='unknown',padlen=-5,status=status)
    if (status /= filtfilt_status_success) error stop 'rank5 SOS gust status'
    ref = sosfiltfilt_gust(sos,x5(2,1,:,1,2),irlen=4)
    if (any(abs(y5(2,1,:,1,2)-ref)>tol)) error stop 'rank5 real values'
    call sosfiltfilt_nd(sos,cx5,cy5,axis=3,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'rank5 mixed real SOS status'
    cref = sosfiltfilt_gust(sos,cx5(2,1,:,1,2))
    if (any(abs(cy5(2,1,:,1,2)-cref)>tol)) error stop 'rank5 mixed values'
    call sosfiltfilt_nd(csos,cx5,cy5,axis=3,method='gust',irlen=3,status=status)
    if (status /= filtfilt_status_success) error stop 'rank5 all complex status'
    cref = sosfiltfilt_gust(csos,cx5(2,1,:,1,2),irlen=3)
    if (any(abs(cy5(2,1,:,1,2)-cref)>tol)) error stop 'rank5 complex values'
    call sosfiltfilt_nd(csos,x5,cy5,axis=3,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'rank5 complex SOS real input status'
    cref = sosfiltfilt_gust(csos,x5(1,1,:,1,1))
    if (any(abs(cy5(1,1,:,1,1)-cref)>tol)) error stop 'rank5 complex SOS values'

    x15 = reshape([(sin(0.063_dp*i),i=1,size(x15))],shape(x15))
    call sosfiltfilt_nd(sos,x15,y15,axis=15,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'rank15 gust status'
    ref = sosfiltfilt_gust(sos,x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,:))
    if (any(abs(y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-ref)>tol)) &
        error stop 'rank15 output'
    call sosfiltfilt_nd(sos,x5(:,:,1:11:2,:,:),y5(:,:,1:11:2,:,:), &
        axis=3,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'strided rank5'
    ref = sosfiltfilt_gust(sos,x5(2,1,1:11:2,1,2))
    if (any(abs(y5(2,1,1:11:2,1,2)-ref)>tol)) error stop 'strided values' 

    xr = [0.75_dp,-0.25_dp]
    yr = 0.0_dp
    call sosfiltfilt_nd(sos,xr,yr,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'short signal status'
    ref = sosfiltfilt_gust(sos,xr)
    if (maxval(abs(yr-ref))>tol) error stop 'short signal reference'
    call sosfiltfilt_nd(sos,xe,ye,axis=2,method='gust',status=status)
    if (status /= filtfilt_status_success) error stop 'zero channels status'
    call sosfiltfilt_nd(sos,xe,ye,axis=2,method='gust',irlen=0,status=status)
    if (status /= filtfilt_status_invalid_input) error stop 'empty invalid irlen'
    y = sosfiltfilt(sos,x,method='unsupported',status=status)
    if (status /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(y))) &
        error stop 'invalid method'
    y = sosfiltfilt(sos,x,method='gust',irlen=0,status=status)
    if (status /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(y))) &
        error stop 'invalid irlen'
    call sosfiltfilt_nd(sos,x5,y5,axis=3,method='unsupported',status=status)
    if (status /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(y5))) &
        error stop 'rank5 invalid method'
    ! For near-singular, untruncated high-order initial-state systems,
    ! return numerical_failure rather than an inaccurate estimate.
    y = sosfiltfilt(high,x(1:20),method='gust',status=status)
    if (status /= filtfilt_status_numerical_failure .or. &
        .not. all(ieee_is_nan(y))) error stop 'ill-conditioned SOS failure status'
    xlong = [(sin(0.08_dp*i)+cos(0.13_dp*i),i=1,240)]
    y = sosfiltfilt(high,xlong,method='gust',irlen=32,status=status)
    if (status /= filtfilt_status_success .or. &
        .not. all(abs(y)<1.0e5_dp)) error stop 'long high-order SOS Gust'

    x5(2,1,7,1,2) = nanval
    call sosfiltfilt_nd(sos,x5,y5,axis=3,method='gust',status=status)
    if (status /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(y5))) &
        error stop 'later channel failed to invalidate output'
    print '(a)', 'test_signal_sos_gust passed'
end program test_signal_sos_gust
