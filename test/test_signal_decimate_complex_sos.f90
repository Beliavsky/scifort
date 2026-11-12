! SPDX-License-Identifier: MIT
program test_signal_decimate_complex_sos
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : decimate, resample_status_success, resample_status_invalid_input
    use scifort_signal_sos, only : sosfilt
    use scifort_signal_filtfilt, only : sosfiltfilt
    use scifort_signal_multirate_nd_into, only : decimate_into
    implicit none

    real(dp), parameter :: tol = 5.0e-12_dp
    real(dp) :: x(91), matrix(2,91), tensor(2,91,2), noncontig(4,91)
    real(dp) :: high_input(1,1,1,1,1,1,1,1,1,1,1,1,1,1,91)
    complex(dp) :: cx(91), sos(2,6), bad_sos(2,6)
    complex(dp) :: output(2,31), scratch(4,31), invalid_output(2,31)
    complex(dp) :: high_output(1,1,1,1,1,1,1,1,1,1,1,1,1,1,31)
    complex(dp) :: high_complex_input(1,1,1,1,1,1,1,1,1,1,1,1,1,1,91)
    complex(dp), allocatable :: y(:), ref(:), y2(:,:), y3(:,:,:)
    integer :: i, j, phase, status

    x = [(sin(0.18_dp*real(i,dp)) + 0.2_dp*cos(0.11_dp*real(i,dp)),i=1,91)]
    cx = cmplx(x,0.27_dp*x+0.08_dp,dp)
    sos(1,:) = [cmplx(0.28_dp,0.07_dp,dp),cmplx(0.22_dp,-0.04_dp,dp), &
                cmplx(0.07_dp,0.02_dp,dp),(1.0_dp,0.0_dp), &
                cmplx(-0.44_dp,0.07_dp,dp),cmplx(0.12_dp,-0.02_dp,dp)]
    sos(2,:) = [cmplx(0.39_dp,-0.03_dp,dp),cmplx(0.16_dp,0.05_dp,dp), &
                cmplx(0.02_dp,-0.01_dp,dp),(1.0_dp,0.0_dp), &
                cmplx(-0.37_dp,-0.06_dp,dp),cmplx(0.06_dp,0.02_dp,dp)]

    do phase=0,1
        y = decimate(cx,3,sos=sos,zero_phase=phase==1,status=status)
        if (status /= resample_status_success .or. size(y) /= 31) error stop 'complex vector status'
        if (phase == 1) then
            ref = sosfiltfilt(sos,cx)
        else
            ref = sosfilt(sos,cx)
        end if
        if (maxval(abs(y-ref(1:91:3))) > tol) error stop 'complex vector values'
        y = decimate(x,3,sos=sos,zero_phase=phase==1,status=status)
        if (status /= resample_status_success .or. size(y) /= 31) error stop 'mixed vector status'
        if (phase == 1) then
            ref = sosfiltfilt(sos,cmplx(x,0.0_dp,dp))
        else
            ref = sosfilt(sos,cmplx(x,0.0_dp,dp))
        end if
        if (maxval(abs(y-ref(1:91:3))) > tol) error stop 'mixed vector values'
    end do

    matrix(1,:) = x
    matrix(2,:) = -2.0_dp*x
    y = decimate(x,3,sos=sos,zero_phase=.false.)
    y2 = decimate(matrix,3,sos=sos,zero_phase=.false.,axis=-1,status=status)
    if (status /= resample_status_success .or. any(shape(y2) /= [2,31])) error stop 'mixed matrix shape'
    if (maxval(abs(y2(1,:)-y)) > tol) error stop 'mixed matrix first lane'
    if (maxval(abs(y2(2,:)+2.0_dp*y)) > tol) error stop 'mixed matrix second lane'
    y2 = decimate(cmplx(matrix,0.14_dp,dp),3,sos=sos, &
                  zero_phase=.false.,axis=2,status=status)
    y = decimate(cmplx(matrix(1,:),0.14_dp,dp),3,sos=sos,zero_phase=.false.)
    if (status /= resample_status_success .or. any(shape(y2) /= [2,31])) error stop 'complex matrix shape'
    if (maxval(abs(y2(1,:)-y)) > tol) error stop 'complex matrix values'

    do j=1,2
        tensor(1,:,j) = real(j,dp)*x
        tensor(2,:,j) = -real(j,dp)*x
    end do
    y3 = decimate(tensor,3,sos=sos,zero_phase=.true.,axis=2,status=status)
    y = decimate(x,3,sos=sos)
    if (status /= resample_status_success .or. any(shape(y3) /= [2,31,2])) error stop 'mixed tensor shape'
    if (maxval(abs(y3(2,:,2)+2.0_dp*y)) > tol) error stop 'mixed tensor values'
    y3 = decimate(cmplx(tensor,0.17_dp,dp),3,sos=sos,axis=2,status=status)
    y = decimate(cmplx(tensor(1,:,1),0.17_dp,dp),3,sos=sos)
    if (status /= resample_status_success .or. any(shape(y3) /= [2,31,2])) error stop 'complex tensor shape'
    if (maxval(abs(y3(1,:,1)-y)) > tol) error stop 'complex tensor values'

    ! Confirm that the first axis is handled, not only the last axis.
    y2 = decimate(transpose(matrix),3,sos=sos,zero_phase=.false.,axis=1,status=status)
    if (status /= resample_status_success .or. any(shape(y2) /= [31,2])) error stop 'axis-one shape'
    y = decimate(x,3,sos=sos,zero_phase=.false.)
    if (maxval(abs(y2(:,2)+2.0_dp*y)) > tol) error stop 'axis-one values'

    call decimate_into(matrix,output,3,sos=sos,axis=-1,status=status)
    y = decimate(x,3,sos=sos)
    if (status /= resample_status_success) error stop 'mixed assumed-rank status'
    if (maxval(abs(output(2,:)+2.0_dp*y)) > tol) error stop 'mixed assumed-rank values'
    call decimate_into(cmplx(matrix,0.15_dp,dp),output,3,sos=sos,axis=2,status=status)
    y = decimate(cmplx(matrix(1,:),0.15_dp,dp),3,sos=sos)
    if (status /= resample_status_success) error stop 'complex assumed-rank status'
    if (maxval(abs(output(1,:)-y)) > tol) error stop 'complex assumed-rank values'

    noncontig = 0.0_dp
    noncontig(1,:) = x
    noncontig(3,:) = -2.0_dp*x
    scratch = cmplx(-9.0_dp,2.0_dp,dp)
    call decimate_into(noncontig(1:4:2,:),scratch(1:4:2,:),3,sos=sos, &
                       zero_phase=.false.,axis=2,status=status)
    y = decimate(x,3,sos=sos,zero_phase=.false.)
    if (status /= resample_status_success) error stop 'strided assumed-rank status'
    if (maxval(abs(scratch(3,:)+2.0_dp*y)) > tol) error stop 'strided assumed-rank values'
    if (any(scratch(2,:) /= cmplx(-9.0_dp,2.0_dp,dp))) error stop 'strided output sentinel'

    high_input(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:) = x
    high_output = cmplx(-12.0_dp,3.0_dp,dp)
    call decimate_into(high_input,high_output,3,sos=sos,status=status)
    y = decimate(x,3,sos=sos)
    if (status /= resample_status_success) error stop 'rank-15 status'
    if (maxval(abs(high_output(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-y)) > tol) error stop 'rank-15 values'
    high_complex_input(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:) = cx
    high_output = cmplx(-12.0_dp,3.0_dp,dp)
    call decimate_into(high_complex_input,high_output,3,sos=sos,status=status)
    y = decimate(cx,3,sos=sos)
    if (status /= resample_status_success) error stop 'rank-15 complex status'
    if (maxval(abs(high_output(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-y)) > tol) error stop 'rank-15 complex values'

    invalid_output = cmplx(-88.0_dp,4.0_dp,dp)
    bad_sos = sos
    bad_sos(1,4) = cmplx(2.0_dp,0.0_dp,dp)
    call decimate_into(matrix,invalid_output,3,sos=bad_sos,status=status)
    if (status /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'invalid a0 preservation'
    bad_sos = sos
    bad_sos(2,2) = cmplx(ieee_value(1.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call decimate_into(matrix,invalid_output,3,sos=bad_sos,status=status)
    if (status /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'invalid nan preservation'
    call decimate_into(matrix,invalid_output,3,sos=sos(:,1:5),status=status)
    if (status /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'invalid shape preservation'
    call decimate_into(matrix,invalid_output,0,sos=sos,status=status)
    if (status /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'invalid q preservation'
    call decimate_into(matrix,invalid_output,3,sos=sos,axis=4,status=status)
    if (status /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'invalid axis preservation'
    call decimate_into(matrix,invalid_output(:,1:30),3,sos=sos,status=status)
    if (status /= resample_status_invalid_input .or. &
        any(invalid_output /= cmplx(-88.0_dp,4.0_dp,dp))) error stop 'invalid output size preservation'
    y = decimate(x,3,sos=bad_sos,status=status)
    if (status /= resample_status_invalid_input .or. size(y) /= 0) error stop 'invalid vector SOS'
    y2 = decimate(matrix,3,sos=sos(:,1:5),status=status)
    if (status /= resample_status_invalid_input .or. size(y2) /= 0) error stop 'invalid matrix SOS'
    y = decimate(x(:5),3,sos=sos,status=status)
    if (status /= resample_status_invalid_input .or. size(y) /= 0) error stop 'short vector SOS'

    print '(a)', 'PASS: complex SOS decimation native regressions'
end program test_signal_decimate_complex_sos
