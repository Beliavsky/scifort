! SPDX-License-Identifier: MIT
! Native regression for the parameterized resample_poly FIR window interface.
program test_signal_resample_poly_param_windows
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : resample_poly, resample_status_success, resample_status_invalid_input
    use scifort_signal_multirate_nd_into, only : resample_poly_into
    implicit none

    integer :: st, j
    real(dp), parameter :: tol = 5.0e-12_dp
    real(dp) :: x(19), a(2,19,3), b(2,2,19,3,2), y(2,2,29,3,2), saved(2,2,29,3,2)
    real(dp) :: incorrect(2,2,30,3,2)
    real(dp), allocatable :: r(:), t(:), r3(:,:,:)
    complex(dp) :: z(19), zin(2,19)
    complex(dp), allocatable :: q(:), cq(:,:), reference(:)

    x = [(sin(real(j,dp)*0.21_dp)+0.1_dp*cos(real(j,dp)*0.47_dp), j=1,19)]
    z = cmplx(x, 0.4_dp*x, dp)
    r = resample_poly(x,3,2,window='gaussian',window_parameter=4.0_dp,status=st)
    if (st /= resample_status_success .or. size(r) /= 29) error stop 'Gaussian scalar'
    t = resample_poly(x,3,2,window='kaiser',window_parameter=5.0_dp,status=st)
    if (st /= resample_status_success) error stop 'Kaiser parameter status'
    r = resample_poly(x,3,2,window='kaiser',beta=5.0_dp,status=st)
    if (maxval(abs(t-r)) > tol) error stop 'Kaiser beta and shape disagree'
    q = resample_poly(z,3,2,window='general_gaussian', &
        window_parameter=1.3_dp,window_parameter2=4.4_dp,status=st)
    if (st /= resample_status_success .or. size(q) /= 29) error stop 'general Gaussian complex'
    if (maxval(abs(aimag(q)-0.4_dp*real(q,dp))) > tol) error stop 'complex linearity'
    do j = 1, 3
        a(:, :, j) = spread(x,1,2)*real(j,dp)
    end do
    r3 = resample_poly(a,3,2,window='tukey',window_parameter=0.65_dp,axis=2,status=st)
    if (st /= resample_status_success .or. any(shape(r3)/=[2,29,3])) error stop 'rank3 shape'
    t = resample_poly(x,3,2,window='tukey',window_parameter=0.65_dp)
    if (maxval(abs(r3(1,:,3)-3.0_dp*t)) > tol) error stop 'rank3 axis values'
    zin(1,:) = z
    zin(2,:) = 2.0_dp*z
    cq = resample_poly(zin,3,2,window='general_cosine', &
        window_coefficients=[0.42_dp,0.50_dp,0.08_dp],axis=-1,status=st)
    if (st /= resample_status_success .or. any(shape(cq)/=[2,29])) error stop 'rank2 complex cosine'
    reference = resample_poly(z,3,2,window='general_cosine', &
        window_coefficients=[0.42_dp,0.50_dp,0.08_dp])
    if (maxval(abs(cq(2,:)-2.0_dp*reference)) > tol) error stop 'rank2 complex values'
    do j = 1, 2
        b(:,j,:,:,1) = a*real(j,dp)
        b(:,j,:,:,2) = -b(:,j,:,:,1)
    end do
    y = -17.0_dp
    call resample_poly_into(b,y,3,2,window='gaussian',window_parameter=4.0_dp,axis=3,status=st)
    if (st /= resample_status_success) error stop 'rank5 status'
    r = resample_poly(x,3,2,window='gaussian',window_parameter=4.0_dp)
    if (maxval(abs(y(2,2,:,3,2)+6.0_dp*r)) > tol) error stop 'rank5 values'
    ! Validation failures must leave caller-allocated rank-N outputs untouched.
    saved = y
    call resample_poly_into(b,y,3,2,window='gaussian', &
        window_parameter2=2.0_dp,axis=3,status=st)
    if (st /= resample_status_invalid_input .or. any(y/=saved)) error stop 'rank5 invalid shape parameter'
    incorrect = 63.0_dp
    call resample_poly_into(b,incorrect,3,2,window='gaussian', &
        window_parameter=3.0_dp,axis=3,status=st)
    if (st /= resample_status_invalid_input .or. any(incorrect/=63.0_dp)) error stop 'invalid output shape'
    r = resample_poly(x,3,2,window='general_gaussian',window_parameter=1.5_dp,status=st)
    if (st /= resample_status_invalid_input .or. size(r)/=0) error stop 'missing second shape parameter'
    r = resample_poly(x,3,2,window='kaiser',beta=5.0_dp,window_parameter=4.0_dp,status=st)
    if (st /= resample_status_invalid_input .or. size(r)/=0) error stop 'conflicting Kaiser beta'
    r = resample_poly(x,3,2,taps=[0.2_dp,0.6_dp,0.2_dp], &
        window='gaussian',window_parameter=3.0_dp,status=st)
    if (st /= resample_status_invalid_input .or. size(r)/=0) error stop 'ambiguous supplied taps'
    print '(a)', 'PASS: parameterized polyphase resampling native regressions'
end program test_signal_resample_poly_param_windows
