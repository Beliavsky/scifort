! SPDX-License-Identifier: MIT
! Gustafsson inverse-state filtering for rank-independent caller-output arrays.
program test_signal_filtfilt_nd_gust
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_filtfilt, only : filtfilt, filtfilt_status_success, &
        filtfilt_status_invalid_input
    use scifort_signal_filtfilt_nd, only : filtfilt_nd
    implicit none
    real(dp) :: b(3), a(3), xr(4,2,7,3,2), yr(4,2,7,3,2)
    real(dp) :: short(2), short_out(2), scalar_in, scalar_out
    real(dp) :: x15(5,1,2,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: y15(5,1,2,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: xe(0,7), ye(0,7), bad_shape(4,2,7,3,1)
    real(dp) :: flat_in(5,2,3), flat_out(5,2,3), row_in(10,3), row_out(10,3)
    complex(dp) :: cb(3), ca(3), xc(4,2,7,3,2), yc(4,2,7,3,2)
    complex(dp) :: complex_out(2), reference(7), cr(2)
    real(dp) :: ref(7), nanval, tol
    integer :: i, k, code

    b = [0.31_dp,0.12_dp,0.03_dp]
    a = [1.0_dp,-0.23_dp,0.07_dp]
    cb = cmplx(b,0.0_dp,dp)
    ca = cmplx(a,0.0_dp,dp)
    cb(2) = cb(2)+(0.0_dp,0.035_dp)
    ca(2) = ca(2)-(0.0_dp,0.02_dp)
    tol = 2.e-10_dp
    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    xr = reshape([(sin(0.071_dp*i)+0.003_dp*i,i=1,size(xr))],shape(xr))
    xc = cmplx(xr,0.07_dp+0.11_dp*xr,dp)

    ! Gustafsson needs no padding, even when signals are shorter than tap count.
    short = [0.7_dp,-0.4_dp]
    call filtfilt_nd(b,a,short,short_out,method='gust',status=code)
    if (code /= filtfilt_status_success) error stop 'short rank1 gust status'
    if (any(abs(short_out-filtfilt(b,a,short,method='gust')) > tol)) &
        error stop 'short rank1 gust values'
    call filtfilt_nd(cb,ca,short,complex_out,method='gust',irlen=1,status=code)
    if (code /= filtfilt_status_success) error stop 'short complex coeff status'
    cr = filtfilt(cb,ca,cmplx(short,0.0_dp,dp),method='gust',irlen=1)
    if (any(abs(complex_out-cr) > tol)) error stop 'short complex coeff values'

    ! Every lane along an interior axis agrees with 1D Gustafsson filtering.
    call filtfilt_nd(b,a,xr,yr,axis=3,method='gust',irlen=3, &
        padtype='invalid',padlen=-25,status=code)
    if (code /= filtfilt_status_success) error stop 'rank5 real gust status'
    do k=1,2
        ref = filtfilt(b,a,xr(3,2,:,2,k),method='gust',irlen=3)
        if (any(abs(yr(3,2,:,2,k)-ref)>tol)) error stop 'rank5 real gust values'
    end do
    ! A one-based first axis with less than the usual pad width is valid.
    call filtfilt_nd(b,a,xr,yr,axis=1,method='gust',status=code)
    if (code /= filtfilt_status_success) error stop 'rank5 first axis status'
    if (any(abs(yr(:,2,4,2,1)-filtfilt(b,a,xr(:,2,4,2,1), &
        method='gust'))>tol)) error stop 'rank5 first axis values'
    call filtfilt_nd(cb,ca,xc,yc,axis=3,method='gust',status=code)
    if (code /= filtfilt_status_success) error stop 'rank5 complex gust status'
    reference = filtfilt(cb,ca,xc(3,2,:,2,1),method='gust')
    if (any(abs(yc(3,2,:,2,1)-reference)>tol)) error stop 'rank5 complex values'

    ! Test all five remaining mixed-coefficient promotions plus RR complex data.
    call filtfilt_nd(b,a,xc,yc,axis=3,method='gust',status=code)
    if (code /= 0) error stop 'mixed rrc status'
    reference = filtfilt(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp), &
        xc(1,2,:,3,2),method='gust')
    if (any(abs(yc(1,2,:,3,2)-reference)>tol)) error stop 'mixed rrc values'
    call filtfilt_nd(b,ca,xr,yc,axis=3,method='gust',status=code)
    if (code /= 0) error stop 'mixed rcr status'
    reference = filtfilt(cmplx(b,0.0_dp,dp),ca,cmplx(xr(1,2,:,3,2),0.0_dp,dp), &
        method='gust')
    if (any(abs(yc(1,2,:,3,2)-reference)>tol)) error stop 'mixed rcr values'
    call filtfilt_nd(b,ca,xc,yc,axis=3,method='gust',status=code)
    if (code /= 0) error stop 'mixed rcc status'
    reference = filtfilt(cmplx(b,0.0_dp,dp),ca,xc(1,2,:,3,2),method='gust')
    if (any(abs(yc(1,2,:,3,2)-reference)>tol)) error stop 'mixed rcc values'
    call filtfilt_nd(cb,a,xr,yc,axis=3,method='gust',status=code)
    if (code /= 0) error stop 'mixed crr status'
    reference = filtfilt(cb,cmplx(a,0.0_dp,dp),cmplx(xr(1,2,:,3,2),0.0_dp,dp), &
        method='gust')
    if (any(abs(yc(1,2,:,3,2)-reference)>tol)) error stop 'mixed crr values'
    call filtfilt_nd(cb,a,xc,yc,axis=3,method='gust',status=code)
    if (code /= 0) error stop 'mixed crc status'
    reference = filtfilt(cb,cmplx(a,0.0_dp,dp),xc(1,2,:,3,2),method='gust')
    if (any(abs(yc(1,2,:,3,2)-reference)>tol)) error stop 'mixed crc values'
    call filtfilt_nd(cb,ca,xr,yc,axis=3,method='gust',status=code)
    if (code /= 0) error stop 'mixed ccr status'
    reference = filtfilt(cb,ca,cmplx(xr(1,2,:,3,2),0.0_dp,dp),method='gust')
    if (any(abs(yc(1,2,:,3,2)-reference)>tol)) error stop 'mixed ccr values'

    ! Sections with discontiguous strides must write only selected cells.
    yr = -987.0_dp
    call filtfilt_nd(b,a,xr(1:4:2,:,:,:,:),yr(1:4:2,:,:,:,:), &
        axis=3,method='gust',irlen=2,status=code)
    if (code /= 0) error stop 'strided section status'
    if (any(yr(2,:,:,:,:) /= -987.0_dp) .or. &
        any(yr(4,:,:,:,:) /= -987.0_dp)) error stop 'strided section overwrite'
    ref = filtfilt(b,a,xr(3,1,:,2,2),method='gust',irlen=2)
    if (any(abs(yr(3,1,:,2,2)-ref)>tol)) error stop 'strided section values'

    ! Rank 15 includes both the first and the last axis, without truncation.
    x15 = reshape([(cos(0.077_dp*i)+0.001_dp*i,i=1,size(x15))],shape(x15))
    flat_in = reshape(x15,[5,2,3])
    call filtfilt_nd(b,a,x15,y15,axis=1,method='gust',status=code)
    if (code /= 0) error stop 'rank15 axis1 status'
    flat_out = reshape(y15,[5,2,3])
    if (any(abs(flat_out(:,2,3)-filtfilt(b,a,flat_in(:,2,3),method='gust'))>tol)) &
        error stop 'rank15 axis1 values'
    call filtfilt_nd(b,a,x15,y15,axis=15,method='gust',irlen=1,status=code)
    if (code /= 0) error stop 'rank15 axis15 status'
    row_in = reshape(x15,[10,3])
    row_out = reshape(y15,[10,3])
    if (any(abs(row_out(4,:)-filtfilt(b,a,row_in(4,:),method='gust',irlen=1))>tol)) &
        error stop 'rank15 axis15 values'

    ! Empty outer lanes are valid, but invalid coefficients are not ignored.
    call filtfilt_nd(b,a,xe,ye,axis=2,method='gust',status=code)
    if (code /= 0) error stop 'empty channels status'
    call filtfilt_nd(b,[0.0_dp],xe,ye,axis=2,method='gust',status=code)
    if (code /= filtfilt_status_invalid_input) error stop 'empty channels invalid coeff'

    call filtfilt_nd(b,a,xr,yr,axis=3,method='missing',status=code)
    if (code /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(yr))) &
        error stop 'invalid method response'
    call filtfilt_nd(b,a,xr,yr,axis=3,method='gust',irlen=0,status=code)
    if (code /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(yr))) &
        error stop 'invalid irlen response'
    call filtfilt_nd(b,a,xr,yr,axis=6,method='gust',status=code)
    if (code /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(yr))) &
        error stop 'invalid axis response'
    call filtfilt_nd(b,a,xr,bad_shape,axis=3,method='gust',status=code)
    if (code /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(bad_shape))) &
        error stop 'wrong output shape response'
    xr(4,2,3,3,2) = nanval
    call filtfilt_nd(b,a,xr,yr,axis=3,method='gust',status=code)
    if (code /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(yr))) &
        error stop 'late-channel invalid input rollback'
    xr(4,2,3,3,2) = 0.0_dp
    ca(2) = cmplx(real(ca(2),dp),nanval,dp)
    call filtfilt_nd(cb,ca,xc,yc,axis=3,method='gust',status=code)
    if (code /= filtfilt_status_invalid_input .or. .not. all(ieee_is_nan(real(yc,dp)))) &
        error stop 'nonfinite complex filter response'
    scalar_in = 1.0_dp
    call filtfilt_nd(b,a,scalar_in,scalar_out,method='gust',status=code)
    if (code /= filtfilt_status_invalid_input .or. .not. ieee_is_nan(scalar_out)) &
        error stop 'scalar response'

    ! Padding method is unchanged and irlen is inert when method='pad'.
    xr = reshape([(sin(0.041_dp*i),i=1,size(xr))],shape(xr))
    call filtfilt_nd(b,a,xr,yr,axis=3,method='pad',padtype='none', &
        irlen=-1,status=code)
    if (code /= 0) error stop 'explicit pad status'
    ref = filtfilt(b,a,xr(2,1,:,2,1),padtype='none')
    if (any(abs(yr(2,1,:,2,1)-ref)>tol)) error stop 'explicit pad values'
    call filtfilt_nd(b,a,xr,yr,axis=3,method='gust',status=code)
    if (code /= 0) error stop 'final gust status'
    print '(a)', 'test_signal_filtfilt_nd_gust passed'
end program test_signal_filtfilt_nd_gust
