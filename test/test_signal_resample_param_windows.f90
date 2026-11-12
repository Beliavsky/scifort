! SPDX-License-Identifier: MIT
! End-to-end native regression for shaped frequency-domain resampling tapers.
program test_signal_resample_param_windows
    use scifort_kinds, only : dp
    use scifort_signal_fourier_resample, only : resample
    use scifort_signal_multirate_nd_into, only : resample_into
    use scifort_signal_spectral, only : get_window
    implicit none
    integer, parameter :: n=9, m=13
    real(dp), parameter :: eps=2.0e-12_dp
    real(dp) :: x(n), x2(n,3), x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,n)
    real(dp) :: y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,m)
    real(dp) :: ybad(m), yout(m,2), ybuffer(m,3)
    real(dp), allocatable :: p(:), r(:), w(:), raw(:), y2(:,:), r2(:,:)
    complex(dp) :: z(n), z3(2,n,2), zout(m), zf(n)
    complex(dp), allocatable :: q(:), zref(:), zres3(:,:,:)
    integer :: i,j,k,st

    x = [(cos(.41_dp*i)+.14_dp*sin(.32_dp*i),i=1,n)]
    z = cmplx(x, [(sin(.21_dp*i), i=1,n)],dp)
    raw=get_window('gaussian',n,parameter=3.2_dp,status=st)
    call require(st==0 .and. size(raw)==n,'gaussian creation')
    w = [(raw(modulo(k+(n+1)/2,n)+1),k=0,n-1)]
    p=resample(x,m,window='gaussian',window_parameter=3.2_dp,status=st)
    call require(st==0 .and. size(p)==m,'real gaussian resample')
    r=resample(x,m,frequency_window=w,status=st)
    call require(st==0 .and. maxval(abs(p-r))<eps,'real window equivalence')
    q=resample(z,m,window='gaussian',window_parameter=3.2_dp,status=st)
    zref=resample(z,m,frequency_window=w,status=st)
    call require(st==0 .and. maxval(abs(q-zref))<eps,'complex window equivalence')

    do j=1,3
        x2(:,j) = x + real(j-1,dp)
    end do
    y2=resample(x2,m,window='gaussian',window_parameter=3.2_dp,axis=1,status=st)
    call require(st==0 .and. size(y2,1)==m .and. size(y2,2)==3,'rank2 result shape')
    do j=1,3
        p=resample(x2(:,j),m,window='gaussian',window_parameter=3.2_dp)
        call require(maxval(abs(p-y2(:,j)))<eps,'rank2 lane')
    end do
    ybuffer = -999.0_dp
    call resample_into(x2,ybuffer,m,window='gaussian',window_parameter=3.2_dp,axis=1,status=st)
    call require(st==0 .and. maxval(abs(ybuffer-y2))<eps,'rank2 caller output')
    yout = -999.0_dp
    call resample_into(x2(1:n:2,1:3:2),yout,m,window='tukey',window_parameter=.36_dp, &
        axis=1,status=st)
    call require(st==0,'strided caller output status')
    do j=1,2
        p=resample(x2(1:n:2,2*j-1),m,window='tukey',window_parameter=.36_dp)
        call require(maxval(abs(yout(:,j)-p))<eps,'strided caller output lane')
    end do
    do j=1,2
        do k=1,2
            z3(j,:,k)=z*cmplx(real(j,dp),real(k-1,dp),dp)
        end do
    end do
    zres3=resample(z3,m,window='general_gaussian',window_parameter=1.9_dp, &
        window_parameter2=3.7_dp,axis=2,status=st)
    call require(st==0 .and. size(zres3,2)==m,'complex rank3 resample')
    do j=1,2
        do k=1,2
            q=resample(z3(j,:,k),m,window='general_gaussian', &
                window_parameter=1.9_dp,window_parameter2=3.7_dp)
            call require(maxval(abs(q-zres3(j,:,k)))<eps,'complex rank3 lane')
        end do
    end do
    x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)=x
    y15=-999.0_dp
    call resample_into(x15,y15,m,axis=15,window='taylor',window_parameter=3.0_dp, &
        window_parameter2=34.0_dp,status=st)
    call require(st==0,'rank15 caller output status')
    p=resample(x,m,window='taylor',window_parameter=3.0_dp,window_parameter2=34.0_dp)
    call require(maxval(abs(y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-p))<eps,'rank15 values')

    ! Frequency-domain input with a shaped window (use independent direct DFT).
    zf=(0.0_dp,0.0_dp)
    do i=1,n
        do j=1,n
            zf(i)=zf(i)+z(j)*exp(cmplx(0.0_dp,-2.0_dp*acos(-1.0_dp)* &
                real((i-1)*(j-1),dp)/real(n,dp),dp))
        end do
    end do
    q=resample(zf,m,window='kaiser',window_parameter=4.2_dp,domain='freq',status=st)
    zref=resample(z,m,window='kaiser',window_parameter=4.2_dp,status=st)
    call require(st==0 .and. maxval(abs(q-zref))<eps,'frequency-domain resample')

    ybad=-789.0_dp
    call resample_into(x,ybad,m,window='boxcar',window_parameter=2.0_dp,status=st)
    call require(st/=0 .and. all(ybad==-789.0_dp),'reject unsupported shape and preserve y')
    call resample_into(x,ybad,m,window_parameter=2.0_dp,status=st)
    call require(st/=0 .and. all(ybad==-789.0_dp),'reject shape without family')
    call resample_into(x,ybad,m,window='gaussian',window_parameter=3.0_dp, &
        frequency_window=w,status=st)
    call require(st/=0 .and. all(ybad==-789.0_dp),'reject conflicting window forms')
    call resample_into(x,ybad,m,window='tukey',window_parameter=3.0_dp, &
        window_parameter2=1.0_dp,status=st)
    call require(st/=0 .and. all(ybad==-789.0_dp),'reject extra family parameter')
    call resample_into(x,ybad,m,window='general_cosine', &
        window_coefficients=[.6_dp,.3_dp,.1_dp], status=st)
    call require(st==0,'general cosine caller-output path')
    p=resample(x,m,window='general_cosine', &
        window_coefficients=[.6_dp,.3_dp,.1_dp])
    call require(maxval(abs(ybad-p))<eps,'general cosine caller output')

    print '(a)', 'PASS: shaped FFT resampling, rank-2/3/15, strided, invalid input'
contains
    subroutine require(condition, description)
        logical, intent(in) :: condition !! true when test passes
        character(*), intent(in) :: description !! diagnostic label
        if (.not. condition) then
            print '(a)', 'FAILED: '//description
            error stop 1
        end if
    end subroutine require
end program test_signal_resample_param_windows
