! SPDX-License-Identifier: MIT
! Native regression for caller-supplied decimate FIR taps and real IIR SOS.
program test_signal_decimate_custom_filters
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : decimate, resample_poly, upfirdn, &
        resample_status_success, resample_status_invalid_input
    use scifort_signal_sos, only : sosfilt
    use scifort_signal_multirate_nd_into, only : decimate_into
    implicit none

    real(dp), parameter :: tol = 2.0e-12_dp
    real(dp) :: x(61), a(2,61), b(2,1,61,2,1), output(2,1,21,2,1), original(2,1,21,2,1)
    real(dp) :: invalid_shape(2,1,22,2,1), h(5), sections(2,6), wrong(1,5)
    real(dp) :: strided_data(4,61), strided_output(4,21)
    real(dp) :: x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,61)
    real(dp) :: y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,21)
    complex(dp) :: z(61), c(2,61,2), outc(2,21,2)
    real(dp), allocatable :: r(:), expected(:), r2(:,:), t(:)
    complex(dp), allocatable :: w(:), w3(:,:,:)
    integer :: i, j, st

    x = [(sin(0.19_dp*real(i,dp)) + 0.15_dp*cos(0.39_dp*real(i,dp)), i=1,61)]
    z = cmplx(x,0.35_dp*x,dp)
    h = [0.07_dp, 0.18_dp, 0.50_dp, 0.18_dp, 0.07_dp]
    sections(1,:) = [0.11_dp, 0.22_dp, 0.11_dp, 1.0_dp, -0.68_dp, 0.21_dp]
    sections(2,:) = [0.18_dp, 0.36_dp, 0.18_dp, 1.0_dp, -0.55_dp, 0.27_dp]

    r = decimate(x,3,taps=h,status=st)
    expected = resample_poly(x,1,3,taps=h)
    if (st /= resample_status_success .or. size(r) /= 21) error stop 'explicit FIR length'
    if (maxval(abs(r-expected)) > tol) error stop 'zero-phase FIR value'
    r = decimate(x,3,taps=h,zero_phase=.false.,status=st)
    t = upfirdn(h,x,down=3)
    if (st /= resample_status_success) error stop 'causal FIR status'
    if (maxval(abs(r-t(:21))) > tol) error stop 'causal FIR value'
    r = decimate(x,3,sos=sections,zero_phase=.false.,status=st)
    t = sosfilt(sections,x)
    if (st /= resample_status_success) error stop 'causal SOS status'
    if (maxval(abs(r-t(1:61:3))) > tol) error stop 'causal SOS value'
    w = decimate(z,3,taps=h,zero_phase=.false.,status=st)
    if (st /= resample_status_success) error stop 'complex FIR status'
    if (maxval(abs(w-cmplx(r*0.0_dp+decimate(x,3,taps=h,zero_phase=.false.), &
        0.35_dp*decimate(x,3,taps=h,zero_phase=.false.), dp))) > tol) error stop 'complex FIR linearity'

    a(1,:) = x
    a(2,:) = 2.0_dp*x
    r2 = decimate(a,3,taps=h,axis=-1,status=st)
    if (st /= resample_status_success .or. any(shape(r2)/=[2,21])) error stop 'rank2 custom FIR shape'
    if (maxval(abs(r2(2,:)-2.0_dp*expected)) > tol) error stop 'rank2 custom FIR value'
    strided_data = 0.0_dp
    strided_data(1,:) = x
    strided_data(3,:) = 2.0_dp*x
    strided_output = -11.0_dp
    call decimate_into(strided_data(1:4:2,:),strided_output(1:4:2,:), &
        3,taps=h,axis=2,status=st)
    if (st /= resample_status_success) error stop 'strided custom FIR status'
    if (maxval(abs(strided_output(3,:)-2.0_dp*expected)) > tol) error stop 'strided output values'
    if (any(strided_output(2,:)/=-11.0_dp)) error stop 'strided output neighbor preservation'
    do j=1,2
        c(1,:,j) = z*real(j,dp)
        c(2,:,j) = -z*real(j,dp)
    end do
    w3 = decimate(c,3,sos=sections,zero_phase=.false.,axis=2,status=st)
    if (st /= resample_status_success .or. any(shape(w3)/=[2,21,2])) error stop 'rank3 SOS shape'
    w = decimate(z,3,sos=sections,zero_phase=.false.)
    if (maxval(abs(w3(2,:,2)+2.0_dp*w)) > tol) error stop 'rank3 SOS values'

    do j=1,2
        b(1,1,:,j,1) = x*real(j,dp)
        b(2,1,:,j,1) = -x*real(j,dp)
    end do
    output = -17.0_dp
    call decimate_into(b,output,3,sos=sections,zero_phase=.false.,axis=3,status=st)
    if (st /= resample_status_success) error stop 'rank5 SOS status'
    t = decimate(x,3,sos=sections,zero_phase=.false.)
    if (maxval(abs(output(2,1,:,2,1)+2.0_dp*t)) > tol) error stop 'rank5 SOS value'
    original=output
    x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:) = x
    y15 = -17.0_dp
    call decimate_into(x15,y15,3,sos=sections,zero_phase=.false.,axis=-1,status=st)
    if (st /= resample_status_success) error stop 'rank15 custom SOS status'
    if (maxval(abs(y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-t)) > tol) error stop 'rank15 custom SOS'
    y15 = -17.0_dp
    call decimate_into(x15,y15,3,taps=h,axis=14,status=st)
    if (st /= resample_status_invalid_input .or. any(y15/=-17.0_dp)) error stop 'rank15 invalid axis'
    call decimate_into(b,output,3,taps=h,sos=sections,axis=3,status=st)
    if (st /= resample_status_invalid_input .or. any(output/=original)) error stop 'conflicting filters'
    call decimate_into(b,output,3,taps=h,ftype='fir',axis=3,status=st)
    if (st /= resample_status_invalid_input .or. any(output/=original)) error stop 'conflicting design'
    invalid_shape = 31.0_dp
    call decimate_into(b,invalid_shape,3,taps=h,axis=3,status=st)
    if (st /= resample_status_invalid_input .or. any(invalid_shape/=31.0_dp)) error stop 'shape preservation'
    wrong=1.0_dp
    r=decimate(x,3,sos=wrong,status=st)
    if (st /= resample_status_invalid_input .or. size(r)/=0) error stop 'bad SOS shape'
    sections(1,4)=2.0_dp
    r=decimate(x,3,sos=sections,status=st)
    if (st /= resample_status_invalid_input .or. size(r)/=0) error stop 'bad SOS a0'
    r=decimate(x,3,taps=[real(dp)::],status=st)
    if (st /= resample_status_invalid_input .or. size(r)/=0) error stop 'empty FIR'
    r=decimate(x,3,taps=h,n=3,status=st)
    if (st /= resample_status_invalid_input .or. size(r)/=0) error stop 'ambiguous order'
    print '(a)', 'PASS: explicit FIR and SOS decimation native regressions'
end program test_signal_decimate_custom_filters
