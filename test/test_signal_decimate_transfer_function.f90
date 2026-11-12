! SPDX-License-Identifier: MIT
program test_signal_decimate_transfer_function
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : decimate, upfirdn, resample_poly, &
        resample_status_success, resample_status_invalid_input
    use scifort_signal_lfilter, only : lfilter
    use scifort_signal_filtfilt, only : filtfilt
    use scifort_signal_multirate_nd_into, only : decimate_into
    implicit none
    real(dp), parameter :: tol = 5.0e-12_dp
    real(dp) :: x(91), h(5), ra(3), rb(3), xr(2,91), yr(2,31)
    real(dp) :: high_in(1,1,1,1,1,1,1,1,1,1,1,1,1,1,91)
    real(dp) :: high_out(1,1,1,1,1,1,1,1,1,1,1,1,1,1,31)
    real(dp) :: scratch_in(4,91), scratch_out(4,31), bad(2,30)
    complex(dp) :: cx(91), cb(3), ca(3), ch(5), cz(2,91,2), cy(2,31,2)
    complex(dp) :: c_high_out(1,1,1,1,1,1,1,1,1,1,1,1,1,1,31)
    complex(dp) :: mixed(2,31), invalid_output(2,31)
    real(dp), allocatable :: ry(:), reference(:)
    complex(dp), allocatable :: zy(:), expected(:), matrix(:,:), cube(:,:,:)
    integer :: i, j, st

    x = [(sin(0.11_dp*real(i,dp)) + 0.3_dp*cos(0.17_dp*real(i,dp)), i=1,91)]
    cx = cmplx(x,0.3_dp*x+0.08_dp,dp)
    h = [0.04_dp,0.2_dp,0.52_dp,0.2_dp,0.04_dp]
    ra = [1.0_dp,-0.56_dp,0.18_dp]
    rb = [0.22_dp,0.09_dp,0.06_dp]
    cb = cmplx(rb, [0.06_dp,0.04_dp,-0.02_dp], dp)
    ca = cmplx(ra, [0.0_dp,0.12_dp,-0.04_dp], dp)
    ch = cmplx(h, [0.03_dp,-0.02_dp,0.0_dp,0.02_dp,-0.03_dp],dp)

    ry = decimate(x,3,b=h,a=[2.0_dp],status=st)
    reference = resample_poly(x,1,3,taps=h/2.0_dp)
    if (st /= resample_status_success .or. size(ry) /= 31) error stop 'real normalized FIR status'
    if (maxval(abs(ry-reference)) > tol) error stop 'real normalized FIR'
    ry = decimate(x,3,b=h,a=[2.0_dp],zero_phase=.false.,status=st)
    reference = upfirdn(h/2.0_dp,x,down=3)
    if (st /= resample_status_success) error stop 'causal real FIR status'
    if (maxval(abs(ry-reference(:31))) > tol) error stop 'causal normalized FIR'

    ry = decimate(x,3,b=rb,a=ra,zero_phase=.false.,status=st)
    reference = lfilter(rb,ra,x)
    if (st /= resample_status_success) error stop 'real IIR status'
    if (maxval(abs(ry-reference(1:91:3))) > tol) error stop 'real IIR'
    ry = decimate(x,3,b=rb,a=ra,status=st)
    reference = filtfilt(rb,ra,x)
    if (st /= resample_status_success) error stop 'zero phase real IIR status'
    if (maxval(abs(ry-reference(1:91:3))) > tol) error stop 'zero phase real IIR'

    zy = decimate(cx,3,b=cb,a=ca,zero_phase=.false.,status=st)
    expected = lfilter(cb,ca,cx)
    if (st /= resample_status_success .or. size(zy) /= 31) error stop 'complex IIR status'
    if (maxval(abs(zy-expected(1:91:3))) > tol) error stop 'complex IIR'
    zy = decimate(cx,3,b=cb,a=ca,status=st)
    expected = filtfilt(cb,ca,cx)
    if (st /= resample_status_success) error stop 'complex zero phase IIR status'
    if (maxval(abs(zy-expected(1:91:3))) > tol) error stop 'complex zero phase IIR'
    zy = decimate(x,3,b=cb,a=ca,zero_phase=.false.,status=st)
    expected = lfilter(cb,ca,cmplx(x,0.0_dp,dp))
    if (st /= resample_status_success) error stop 'real to complex status'
    if (maxval(abs(zy-expected(1:91:3))) > tol) error stop 'real to complex IIR'
    zy = decimate(cx,3,b=ch,a=[cmplx(2.0_dp,0.3_dp,dp)],status=st)
    expected = resample_poly(cx,1,3,taps=ch/cmplx(2.0_dp,0.3_dp,dp))
    if (st /= resample_status_success) error stop 'complex FIR status'
    if (maxval(abs(zy-expected)) > tol) error stop 'complex FIR'

    xr(1,:) = x
    xr(2,:) = -2.0_dp*x
    matrix = decimate(xr,3,b=cb,a=ca,zero_phase=.false.,axis=-1,status=st)
    zy = decimate(x,3,b=cb,a=ca,zero_phase=.false.)
    if (st /= resample_status_success .or. any(shape(matrix)/=[2,31])) error stop 'rank2 mixed shape'
    if (maxval(abs(matrix(2,:)+2.0_dp*zy)) > tol) error stop 'rank2 mixed values'
    cz = (0.0_dp,0.0_dp)
    do j = 1,2
        cz(1,:,j) = real(j,dp)*cx
        cz(2,:,j) = -real(j,dp)*cx
    end do
    cube = decimate(cz,3,b=cb,a=ca,zero_phase=.false.,axis=2,status=st)
    expected = decimate(cx,3,b=cb,a=ca,zero_phase=.false.)
    if (st /= resample_status_success .or. any(shape(cube)/=[2,31,2])) error stop 'rank3 complex shape'
    if (maxval(abs(cube(2,:,2)+2.0_dp*expected)) > tol) error stop 'rank3 complex values'
    call decimate_into(xr,mixed,3,b=cb,a=ca,zero_phase=.false.,axis=2,status=st)
    if (st /= resample_status_success) error stop 'rank2 mixed into status'
    if (maxval(abs(mixed(2,:)+2.0_dp*zy)) > tol) error stop 'rank2 mixed into'

    scratch_in = 0.0_dp
    scratch_in(1,:) = x
    scratch_in(3,:) = -2.0_dp*x
    scratch_out = -19.0_dp
    call decimate_into(scratch_in(1:4:2,:), scratch_out(1:4:2,:), 3, &
        b=rb,a=ra,zero_phase=.false.,axis=-1,status=st)
    ry = decimate(x,3,b=rb,a=ra,zero_phase=.false.)
    if (st /= resample_status_success) error stop 'strided into status'
    if (maxval(abs(scratch_out(3,:)+2.0_dp*ry)) > tol) error stop 'strided into'
    if (any(scratch_out(2,:)/=-19.0_dp)) error stop 'strided sentinel'
    high_in(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:) = x
    high_out = -8.0_dp
    call decimate_into(high_in,high_out,3,b=rb,a=ra,axis=-1,status=st)
    ry = decimate(x,3,b=rb,a=ra)
    if (st /= resample_status_success) error stop 'rank15 real status'
    if (maxval(abs(high_out(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-ry)) > tol) error stop 'rank15 real'
    c_high_out = cmplx(-8.0_dp,2.0_dp,dp)
    call decimate_into(high_in,c_high_out,3,b=cb,a=ca,axis=-1,status=st)
    zy = decimate(x,3,b=cb,a=ca)
    if (st /= resample_status_success) error stop 'rank15 complex status'
    if (maxval(abs(c_high_out(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-zy)) > tol) error stop 'rank15 complex'

    invalid_output = cmplx(-88.0_dp,4.0_dp,dp)
    call decimate_into(xr,invalid_output,3,b=cb,a=[(0.0_dp,0.0_dp)],status=st)
    if (st /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'bad a preservation'
    call decimate_into(cmplx(xr,0.0_dp,dp),invalid_output,3,b=cb,a=ca,taps=h,status=st)
    if (st /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'not valid combination'
    call decimate_into(cmplx(xr,0.0_dp,dp),invalid_output,3,b=cb,status=st)
    if (st /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'missing denominator'
    ry = decimate(x,3,b=rb,status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'missing real a'
    ry = decimate(x,3,b=[real(dp)::],a=ra,status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'empty b'
    ry = decimate(x,3,b=rb,a=[real(dp)::],status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'empty a'
    ry = decimate(x,3,b=rb,a=ra,ftype='iir',status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'ambiguous iir'
    ry = decimate(x,3,b=rb,a=ra,n=2,status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'ambiguous n'
    ry = decimate(x,3,b=rb,a=ra,taps=h,status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'ambiguous FIR'
    ry = decimate(x,3,b=rb,a=ra,sos=reshape([1.0_dp,0.0_dp,0.0_dp, &
        1.0_dp,0.0_dp,0.0_dp],[1,6]),status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'ambiguous SOS'
    ry = decimate(x,3,b=rb,a=[0.0_dp,0.2_dp],status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'a0 zero'
    ry = decimate(x,3,b=[ieee_value(1.0_dp,ieee_quiet_nan)],a=[1.0_dp],status=st)
    if (st /= resample_status_invalid_input .or. size(ry) /= 0) error stop 'nan numerator'
    bad = -11.0_dp
    call decimate_into(xr,bad,3,b=rb,a=ra,status=st)
    if (st /= resample_status_invalid_input .or. any(bad/=-11.0_dp)) error stop 'output shape preservation'
    print '(a)', 'PASS: transfer-function decimation native regressions'
end program test_signal_decimate_transfer_function
