! SPDX-License-Identifier: MIT
! User-specified Fourier-bin window weights for periodic resampling.
program test_signal_resample_frequency_window
    use scifort_kinds, only : dp
    use scifort_signal, only : resample, resample_into, &
        fourier_resample_status_success, fourier_resample_status_invalid_input
    implicit none
    real(dp), parameter :: x(7) = [1.0_dp,2.0_dp,-1.0_dp,4.0_dp,0.0_dp,1.0_dp,3.0_dp]
    real(dp), parameter :: w(7) = [0.6_dp,1.1_dp,0.1_dp,1.3_dp,0.8_dp,0.4_dp,1.7_dp]
    real(dp), parameter :: expected(11) = [ &
        -0.08580939769522658_dp,2.356169399146915_dp,1.134013794531163_dp, &
        -1.584826087107771_dp,0.332524147959336_dp,2.819507086115153_dp, &
        0.5135909507745023_dp,-1.15405759438125_dp,1.818795753212278_dp, &
        3.058311042102508_dp,0.2203523339138203_dp]
    complex(dp), parameter :: z(6) = [ &
        cmplx(1.0_dp,1.0_dp,dp),cmplx(2.0_dp,-0.7_dp,dp), &
        cmplx(3.0_dp,2.0_dp,dp),cmplx(-1.0_dp,0.2_dp,dp), &
        cmplx(2.0_dp,-0.5_dp,dp),cmplx(0.2_dp,3.0_dp,dp)]
    complex(dp), parameter :: cw(6) = [ &
        cmplx(1.0_dp,0.0_dp,dp),cmplx(0.4_dp,-0.2_dp,dp), &
        cmplx(0.8_dp,0.1_dp,dp),cmplx(1.1_dp,0.3_dp,dp), &
        cmplx(-0.7_dp,0.2_dp,dp),cmplx(1.3_dp,-0.4_dp,dp)]
    complex(dp), parameter :: expected_complex(4) = [ &
        cmplx(0.3478272291180741_dp,1.167040531072257_dp,dp), &
        cmplx(3.322010153131195_dp,0.4417800386861434_dp,dp), &
        cmplx(-0.7055999916174193_dp,-0.05113877580391557_dp,dp), &
        cmplx(1.835762609368151_dp,1.775651539378848_dp,dp)]
    real(dp) :: r2(7,3), r3(2,7,2), r4(2,3,7,2), out4(2,3,11,2)
    real(dp) :: r15(1,1,1,1,1,1,1,7,1,1,1,1,1,1,1)
    real(dp) :: out15(1,1,1,1,1,1,1,11,1,1,1,1,1,1,1)
    real(dp) :: buffer(2,3,22,2), old(2,3,22,2)
    complex(dp) :: c2(6,2), o2(4,2), o_single(4,1)
    real(dp), allocatable :: y(:), t(:), y2(:,:), y3(:,:,:)
    complex(dp), allocatable :: yc(:), yc2(:,:)
    integer :: i, j, k, st
    real(dp), parameter :: tol = 5.0e-12_dp

    y = resample(x,11,frequency_window=w,status=st)
    if (st /= fourier_resample_status_success) error stop 'real window status'
    call check(y,expected,tol)
    yc = resample(z,4,complex_window=cw,status=st)
    if (st /= fourier_resample_status_success) error stop 'complex window status'
    if (maxval(abs(yc-expected_complex)) > tol) error stop 'complex window values'
    ! Real complex-input weights also use the two-sided (unfolded) spectrum.
    yc = resample(z,4,frequency_window=real(cw,dp),status=st)
    if (st /= fourier_resample_status_success) error stop 'complex/real window'
    yc2 = reshape([yc,yc],[4,2])
    c2(:,1) = z
    c2(:,2) = z
    call resample_into(c2,o2,4,frequency_window=real(cw,dp),status=st)
    if (st /= fourier_resample_status_success) error stop 'complex into real weights'
    if (maxval(abs(o2-yc2)) > tol) error stop 'complex into real values'

    do j=1,3
        r2(:,j) = x*real(j,dp)
    end do
    y2 = resample(r2,11,frequency_window=w,axis=1,status=st)
    if (st /= fourier_resample_status_success) error stop 'rank2 array window'
    do j=1,3
        call check(y2(:,j),real(j,dp)*expected,tol)
    end do
    do i=1,2
        do k=1,2
            r3(i,:,k) = x*real(i+k,dp)
        end do
    end do
    y3 = resample(r3,11,frequency_window=w,axis=2,status=st)
    if (st /= fourier_resample_status_success) error stop 'rank3 array window'
    do i=1,2
        do k=1,2
            call check(y3(i,:,k),real(i+k,dp)*expected,tol)
        end do
    end do
    do i=1,2
        do j=1,3
            do k=1,2
                r4(i,j,:,k) = x*real(i+j+k,dp)
            end do
        end do
    end do
    call resample_into(r4,out4,11,axis=3,frequency_window=w,status=st)
    if (st /= fourier_resample_status_success) error stop 'rank4 real window'
    do i=1,2
        do j=1,3
            do k=1,2
                call check(out4(i,j,:,k),real(i+j+k,dp)*expected,tol)
            end do
        end do
    end do
    r15(1,1,1,1,1,1,1,:,1,1,1,1,1,1,1) = x
    call resample_into(r15,out15,11,axis=8,frequency_window=w,status=st)
    if (st /= fourier_resample_status_success) error stop 'rank15 real window'
    call check(out15(1,1,1,1,1,1,1,:,1,1,1,1,1,1,1),expected,tol)

    ! A noncontiguous output view is writable, even with noncontiguous input.
    buffer = -17.0_dp
    call resample_into(r4(:, :, :, :),buffer(:,:,2:22:2,:),11,axis=3, &
        frequency_window=w,status=st)
    if (st /= fourier_resample_status_success) error stop 'strided output'
    if (maxval(abs(buffer(:,:,2:22:2,:)-out4)) > tol) error stop 'strided results'
    if (any(buffer(:,:,1:21:2,:) /= -17.0_dp)) error stop 'strided touched gaps'

    ! Invalid weights and conflicting forms leave the caller output unchanged.
    buffer = -13.0_dp
    old = buffer
    call resample_into(r4,buffer(:,:,2:22:2,:),11,axis=3, &
        frequency_window=w(:6),status=st)
    if (st /= fourier_resample_status_invalid_input) error stop 'bad weights status'
    if (any(buffer/=old)) error stop 'invalid resample mutated output'
    call resample_into(r4,out4,11,axis=3,frequency_window=w,window='hann',status=st)
    if (st /= fourier_resample_status_invalid_input) error stop 'conflicting window status'
    y = resample(x,11,frequency_window=w(:6),status=st)
    if (st /= fourier_resample_status_invalid_input .or. size(y) /= 0) &
        error stop 'invalid window length'
    yc = resample(z,4,frequency_window=real(cw,dp),complex_window=cw,status=st)
    if (st /= fourier_resample_status_invalid_input .or. size(yc) /= 0) &
        error stop 'conflicting bin weights'
    y = resample(x,11,window='boxcar',status=st)
    if (st /= fourier_resample_status_success) error stop 'named window backward compatibility'
    call check(y,resample(x,11),tol)
    yc = resample(z,4,complex_window=cw,status=st)
    call resample_into(reshape(z,[6,1]),o_single,4,axis=1,complex_window=cw,status=st)
    if (st /= fourier_resample_status_success) error stop 'complex into complex weight'
    if (maxval(abs(o_single(:,1)-yc)) > tol) error stop 'complex into complex values'
    print '(a)', 'test_signal_resample_frequency_window passed'
contains
    subroutine check(a,b,eps)
        real(dp), intent(in) :: a(:) !! calculated array
        real(dp), intent(in) :: b(:) !! reference array
        real(dp), intent(in) :: eps !! tolerance
        if (size(a) /= size(b)) error stop 'resampling output shape'
        if (size(a) > 0) then
            if (maxval(abs(a-b)) > eps) then
                print *, 'calculated: ',a
                print *, 'expected: ',b
                error stop 'resampling numerical mismatch'
            end if
        end if
    end subroutine check
end program test_signal_resample_frequency_window
