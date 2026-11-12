! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_czt_plans
    use scifort_signal, only : dp, czt, zoom_fft, czt_plan, zoom_fft_plan, &
        czt_status_success, czt_status_invalid_input
    implicit none
    type(czt_plan) :: plan, duplicate
    type(zoom_fft_plan) :: zoom
    complex(dp), allocatable :: y(:), expected(:), points(:), matrix(:,:)
    real(dp) :: x(9), xr(9,3), xt(3,9), f0(2)
    complex(dp) :: z(9), zmat(9,3), ztran(3,9), w, a
    integer :: i, j, stat, mode

    x = [(cos(0.73_dp*real(i,dp))+real(i,dp)*0.2_dp,i=1,9)]
    z = cmplx(x,[(sin(0.23_dp*real(i,dp)),i=1,9)],dp)
    do j=1,3
        xr(:,j) = x*real(j,dp)
        zmat(:,j) = z*real(j,dp)
    end do
    xt = transpose(xr)
    ztran = transpose(zmat)

    ! Unit-contour FFT plan, default DFT and custom rotated contours.
    do mode=1,4
        select case(mode)
        case(1)
            w = exp(cmplx(0.0_dp,-2.0_dp*acos(-1.0_dp)/9.0_dp,dp))
            a = cmplx(1.0_dp,0.0_dp,dp)
            call plan%initialize(9,status=stat)
            expected = czt(z)
        case(2)
            w = exp(cmplx(0.0_dp,0.18_dp,dp))
            a = exp(cmplx(0.0_dp,-0.32_dp,dp))
            call plan%initialize(9,m=13,w=w,a=a,status=stat)
            expected = czt(z,m=13,w=w,a=a)
        case(3)
            w = cmplx(0.977_dp,-0.13_dp,dp)
            a = cmplx(1.02_dp,0.07_dp,dp)
            call plan%initialize(9,m=11,w=w,a=a,status=stat)
            expected = czt(z,m=11,w=w,a=a)
        case default
            call plan%initialize(1,m=12,status=stat)
            y = plan%transform([z(1)],status=stat)
            if (stat/=0 .or. any(y/=z(1))) error stop 'single sample plan'
            cycle
        end select
        if (stat /= czt_status_success .or. .not. plan%is_valid()) error stop 'plan initialization'
        if (plan%input_length()/=9 .or. plan%output_length()/=size(expected)) error stop 'plan lengths'
        y = plan%transform(z,status=stat)
        if (stat/=0 .or. maxval(abs(y-expected)) > 2.e-12_dp) error stop 'planned CZT'
        y = plan%transform(x,status=stat)
        if (stat/=0 .or. maxval(abs(y-czt(x,m=plan%output_length(), &
            w=merge(w,exp(cmplx(0.0_dp,-2.0_dp*acos(-1.0_dp)/9.0_dp,dp)), mode/=1), &
            a=merge(a,cmplx(1.0_dp,0.0_dp,dp),mode/=1)))) > 2.e-12_dp) error stop 'real CZT'
        duplicate = plan
        y = duplicate%transform(z,status=stat)
        if (stat/=0 .or. maxval(abs(y-expected)) > 2.e-12_dp) error stop 'copied plan'
        matrix = plan%transform(zmat,axis=1,status=stat)
        if (stat/=0 .or. size(matrix,1)/=size(expected) .or. size(matrix,2)/=3) error stop 'axis 1 shape'
        do j=1,3
            if (maxval(abs(matrix(:,j)-real(j,dp)*expected)) > 2.e-12_dp) error stop 'axis 1'
        end do
        matrix = plan%transform(ztran,axis=2,status=stat)
        if (stat/=0 .or. size(matrix,1)/=3 .or. size(matrix,2)/=size(expected)) error stop 'axis 2 shape'
        do j=1,3
            if (maxval(abs(matrix(j,:)-real(j,dp)*expected)) > 2.e-12_dp) error stop 'axis 2'
        end do
        matrix = plan%transform(xr,axis=1,status=stat)
        if (stat/=0) error stop 'real rank-2 axis 1'
        matrix = plan%transform(xt,axis=2,status=stat)
        if (stat/=0) error stop 'real rank-2 axis 2'
        points = plan%points(status=stat)
        if (stat/=0 .or. size(points)/=plan%output_length()) error stop 'CZT points'
    end do

    call plan%initialize(9,m=0,status=stat)
    if (stat /= czt_status_invalid_input .or. plan%is_valid()) error stop 'invalid reinitialize'
    y = plan%transform(z,status=stat)
    if (stat /= czt_status_invalid_input .or. size(y)/=0) error stop 'uninitialized transform'
    call plan%initialize(9,status=stat)
    y = plan%transform(z(1:5),status=stat)
    if (stat /= czt_status_invalid_input .or. size(y)/=0) error stop 'wrong size'
    matrix = plan%transform(zmat,axis=3,status=stat)
    if (stat /= czt_status_invalid_input .or. size(matrix)/=0) error stop 'invalid axis'

    f0 = [-0.75_dp, 2.25_dp]
    do mode=1,2
        call zoom%initialize(9,f0,m=15,fs=12.0_dp,endpoint=mode==2,status=stat)
        if (stat/=0 .or. .not. zoom%is_valid()) error stop 'zoom initialization'
        y = zoom%transform(z,status=stat)
        expected = zoom_fft(z,f0,m=15,fs=12.0_dp,endpoint=mode==2)
        if (stat/=0 .or. maxval(abs(y-expected))>2.e-12_dp) error stop 'zoom plan'
        matrix = zoom%transform(ztran,axis=2,status=stat)
        if (stat/=0 .or. maxval(abs(matrix(2,:)-2.0_dp*expected))>2.e-12_dp) error stop 'zoom batch'
        if (mode==2) then
            block
                real(dp), allocatable :: fr(:)
                fr = zoom%frequencies()
                if (abs(fr(1)-f0(1))>1.e-14_dp .or. &
                    abs(fr(15)-f0(2))>1.e-14_dp) error stop 'frequency coordinates'
            end block
        end if
    end do
    call zoom%initialize(9,2.0_dp,m=5,fs=10.0_dp,status=stat)
    if (stat/=0) error stop 'zoom scalar initialization'
    y = zoom%transform(x,status=stat)
    expected = zoom_fft(x,2.0_dp,m=5,fs=10.0_dp)
    if (stat/=0 .or. maxval(abs(y-expected))>2.e-12_dp) error stop 'zoom scalar'
    call zoom%initialize(9,f0,fs=0.0_dp,status=stat)
    if (stat/=czt_status_invalid_input .or. zoom%is_valid()) error stop 'invalid zoom fs'
    print '(a)', 'PASS: reusable CZT and ZoomFFT plans'
end program test_signal_czt_plans
