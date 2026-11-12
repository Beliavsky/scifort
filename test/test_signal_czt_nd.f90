! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_czt_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : czt, zoom_fft, czt_plan, zoom_fft_plan, &
        czt_into, czt_plan_into, zoom_fft_into, zoom_fft_plan_into, &
        czt_status_success, czt_status_invalid_input, czt_status_nonfinite_result
    implicit none
    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: x(2,7,3), xt(2,3,7), xr15(1,1,7,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: xe(0,7)
    complex(dp) :: z(2,7,3), bad(2,7,3), yr(2,11,3), yr2(2,11,3)
    complex(dp) :: out(2,11,3), buf(2,22,3), saved(2,22,3)
    complex(dp) :: out15(1,1,11,1,1,1,1,1,1,1,1,1,1,1,1), ye(0,9)
    complex(dp) :: yt(2,3,7), w, a
    complex(dp), allocatable :: expect(:)
    type(czt_plan) :: plan
    type(zoom_fft_plan) :: zoom
    integer :: i, j, k, st

    do k=1,3
        do j=1,7
            do i=1,2
                x(i,j,k) = cos(0.2_dp*real(i*j+k,dp)) + 0.15_dp*real(i*j-k,dp)
            end do
        end do
    end do
    z = cmplx(x,0.3_dp*sin(2.0_dp*x),dp)
    w = cmplx(0.94_dp,-0.05_dp,dp)
    a = cmplx(1.03_dp,0.04_dp,dp)
    call czt_into(x,yr,m=11,axis=2,status=st)
    if (st/=czt_status_success) error stop 'real rank-3 transform'
    do k=1,3
        do i=1,2
            expect = czt(x(i,:,k),m=11)
            if (maxval(abs(yr(i,:,k)-expect))>tol) error stop 'real rank-3 values'
        end do
    end do
    call czt_into(z,out,m=11,w=w,a=a,axis=-2,status=st)
    if (st/=czt_status_success) error stop 'complex spiral transform'
    do k=1,3
        do i=1,2
            expect = czt(z(i,:,k),m=11,w=w,a=a)
            if (maxval(abs(out(i,:,k)-expect))>tol) error stop 'spiral values'
        end do
    end do
    call plan%initialize(7,m=11,w=w,a=a,status=st)
    if (st/=0) error stop 'plan initialization'
    call czt_plan_into(plan,z,yr2,axis=2,status=st)
    if (st/=0 .or. maxval(abs(yr2-out))>tol) error stop 'planned complex transform'
    call czt_plan_into(plan,x,yr2,axis=2,status=st)
    if (st/=0) error stop 'planned real transform'
    do k=1,3
        do i=1,2
            expect = czt(x(i,:,k),m=11,w=w,a=a)
            if (maxval(abs(yr2(i,:,k)-expect))>tol) error stop 'planned real values'
        end do
    end do
    ! Rank-15 axis handling is independent of singleton channel dimensions.
    xr15(1,1,:,1,1,1,1,1,1,1,1,1,1,1,1) = x(1,:,1)
    call czt_into(xr15,out15,m=11,axis=3,status=st)
    if (st/=0) error stop 'rank-15 CZT'
    expect = czt(x(1,:,1),m=11)
    if (maxval(abs(out15(1,1,:,1,1,1,1,1,1,1,1,1,1,1,1)-expect))>tol) error stop 'rank-15 values'
    ! Negative axes and the SciPy default (last axis).
    do i=1,2
        do k=1,3
            xt(i,k,:) = x(i,:,k)
        end do
    end do
    call czt_into(xt,yt,status=st)
    if (st/=0) error stop 'default last-axis CZT'
    do i=1,2
        do k=1,3
            expect = czt(xt(i,k,:))
            if (maxval(abs(yt(i,k,:)-expect))>tol) error stop 'default last-axis values'
        end do
    end do
    buf = cmplx(-37.0_dp,-21.0_dp,dp)
    call czt_into(z,buf(:,2:22:2,:),m=11,w=w,a=a,axis=2,status=st)
    if (st/=0) error stop 'noncontiguous CZT output'
    if (maxval(abs(buf(:,2:22:2,:)-out))>tol) error stop 'noncontiguous output values'
    if (any(buf(:,1:21:2,:)/=cmplx(-37.0_dp,-21.0_dp,dp))) error stop 'noncontiguous gaps'
    saved=buf
    call czt_into(z,buf(:,2:22:2,:),m=11,axis=4,status=st)
    if (st/=czt_status_invalid_input .or. any(buf/=saved)) error stop 'invalid axis preservation'
    call czt_into(z,buf(:,2:22:2,:),m=10,axis=2,status=st)
    if (st/=czt_status_invalid_input .or. any(buf/=saved)) error stop 'wrong m preservation'
    call czt_into(z,buf(:,2:22:2,:),m=11,w=(0.0_dp,0.0_dp),axis=2,status=st)
    if (st/=czt_status_invalid_input .or. any(buf/=saved)) error stop 'invalid spiral preservation'
    bad=z
    bad(2,3,2)=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call czt_plan_into(plan,bad,buf(:,2:22:2,:),axis=2,status=st)
    if (st/=czt_status_invalid_input .or. any(buf/=saved)) error stop 'invalid later channel preservation'
    call czt_plan_into(plan,z,buf(:,2:22:2,:),axis=2,status=st)
    if (st/=0 .or. maxval(abs(buf(:,2:22:2,:)-out))>tol) error stop 'valid plan after failure'
    ! Empty channel dimensions are legal and preserve shape.
    call czt_into(xe,ye,m=9,axis=2,status=st)
    if (st/=0 .or. size(ye)/=0) error stop 'empty channels'

    call zoom_fft_into(z,out,[-1.25_dp,2.0_dp],m=11,fs=12.0_dp,endpoint=.true.,axis=2,status=st)
    if (st/=0) error stop 'zoom complex rank-3'
    do k=1,3
        do i=1,2
            expect=zoom_fft(z(i,:,k),[-1.25_dp,2.0_dp],m=11,fs=12.0_dp,endpoint=.true.)
            if (maxval(abs(out(i,:,k)-expect))>tol) error stop 'zoom complex values'
        end do
    end do
    call zoom%initialize(7,[-1.25_dp,2.0_dp],m=11,fs=12.0_dp,endpoint=.true.,status=st)
    if (st/=0) error stop 'zoom plan initialization'
    call zoom_fft_plan_into(zoom,z,yr2,axis=2,status=st)
    if (st/=0 .or. maxval(abs(out-yr2))>tol) error stop 'zoom planned complex'
    call zoom_fft_plan_into(zoom,x,yr2,axis=2,status=st)
    if (st/=0) error stop 'zoom planned real'
    do k=1,3
        do i=1,2
            expect=zoom_fft(x(i,:,k),[-1.25_dp,2.0_dp],m=11,fs=12.0_dp,endpoint=.true.)
            if (maxval(abs(yr2(i,:,k)-expect))>tol) error stop 'zoom planned real values'
        end do
    end do
    call zoom_fft_into(x,yr2,2.5_dp,m=11,fs=8.0_dp,axis=2,status=st)
    if (st/=0) error stop 'zoom scalar'
    do k=1,3
        do i=1,2
            expect=zoom_fft(x(i,:,k),2.5_dp,m=11,fs=8.0_dp)
            if (maxval(abs(yr2(i,:,k)-expect))>tol) error stop 'zoom scalar values'
        end do
    end do
    call zoom_fft_into(z,yr2,2.5_dp,m=11,fs=8.0_dp,axis=2,status=st)
    if (st/=0) error stop 'zoom scalar complex'
    call zoom_fft_into(xr15,out15,[0.0_dp,2.5_dp],m=11,fs=8.0_dp,axis=3,status=st)
    if (st/=0) error stop 'zoom rank-15'
    expect=zoom_fft(x(1,:,1),[0.0_dp,2.5_dp],m=11,fs=8.0_dp)
    if (maxval(abs(out15(1,1,:,1,1,1,1,1,1,1,1,1,1,1,1)-expect))>tol) error stop 'zoom rank-15 values'
    buf=cmplx(-10.0_dp,-7.0_dp,dp)
    saved=buf
    call zoom_fft_into(x,buf(:,2:22:2,:),[-1.0_dp,1.0_dp],m=11,fs=0.0_dp,axis=2,status=st)
    if (st/=czt_status_invalid_input .or. any(buf/=saved)) error stop 'zoom invalid fs preservation'
    call zoom_fft_into(x,buf(:,2:22:2,:),[-1.0_dp,1.0_dp],m=11,fs=8.0_dp,axis=2,status=st)
    if (st/=0) error stop 'zoom noncontiguous output'
    if (any(buf(:,1:21:2,:)/=cmplx(-10.0_dp,-7.0_dp,dp))) error stop 'zoom noncontiguous gaps'
    print '(a)', 'PASS: rank-independent CZT and zoom FFT (rank 1-15, real/complex, transactional)'
end program test_signal_czt_nd
