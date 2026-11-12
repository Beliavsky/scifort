! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multirate_assumed_rank
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, resample, upfirdn, resample_poly, decimate, &
        resample_into, resample_poly_into, upfirdn_into, decimate_into, &
        resample_status_success, resample_status_invalid_input
    implicit none
    real(dp) :: x4(3,7,2,2), y4(3,11,2,2), bad4(3,10,2,2), v4(3,7,2,5)
    real(dp) :: scratch(6,11,2,2), h(3)
    real(dp) :: x6(2,2,2,8,1,2), y6(2,2,2,4,1,2)
    real(dp) :: x8(1,1,2,7,1,1,1,2), y8(1,1,2,11,1,1,1,2)
    real(dp) :: x14(1,1,1,1,1,1,1,1,1,1,1,1,2,7)
    real(dp) :: y14(1,1,1,1,1,1,1,1,1,1,1,1,2,11)
    real(dp) :: x15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,7)
    real(dp) :: y15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,15)
    complex(dp) :: z4(3,7,2,2), z4out(3,11,2,2), hc(3)
    complex(dp) :: z5(2,1,7,2,2), z5out(2,1,11,2,2)
    complex(dp) :: v5(2,1,11,2,2)
    real(dp), allocatable :: times(:), line(:)
    complex(dp), allocatable :: cline(:)
    integer :: i,j,k,l,st
    real(dp), parameter :: tol=3.0e-11_dp
    h=[0.2_dp,0.6_dp,0.2_dp]
    hc=cmplx(h,[0.03_dp,-0.1_dp,0.2_dp],dp)
    do l=1,2
        do k=1,2
            do j=1,7
                do i=1,3
                    x4(i,j,k,l)=sin(0.17_dp*real(i+3*j+k*7+l,dp))
                    z4(i,j,k,l)=cmplx(x4(i,j,k,l),cos(0.1_dp*real(i+j+k+l,dp)),dp)
                end do
            end do
        end do
    end do
    y4=-100.0_dp
    call resample_into(x4,y4,11,axis=2,t=[0.0_dp,0.25_dp],output_times=times,status=st)
    if (st/=resample_status_success .or. size(times)/=11) error stop 'resample rank4 status'
    if (abs(times(11)-real(70,dp)/44.0_dp)>tol) error stop 'resample rank4 times'
    do l=1,2
        do k=1,2
            do i=1,3
                line=resample(x4(i,:,k,l),11)
                if (maxval(abs(line-y4(i,:,k,l)))>tol) error stop 'resample rank4 values'
            end do
        end do
    end do
    call resample_into(z4,z4out,11,window='hann',axis=-3,status=st)
    if (st/=resample_status_success) error stop 'complex resample rank4 status'
    cline=resample(z4(2,:,1,2),11,window='hann')
    if (maxval(abs(cline-z4out(2,:,1,2)))>tol) error stop 'complex resample rank4 values'

    ! Shape and invalid lane failures must not modify the caller's output.
    bad4=-200.0_dp
    call resample_into(x4,bad4,11,axis=2,status=st)
    if (st/=resample_status_invalid_input .or. any(bad4/=-200.0_dp)) error stop 'wrong output extent'
    y4=-300.0_dp
    call resample_into(x4,y4,11,axis=-5,status=st)
    if (st/=resample_status_invalid_input .or. any(y4/=-300.0_dp)) error stop 'invalid axis'
    z4(3,7,2,2)=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    z4out=cmplx(-301.0_dp,0.0_dp,dp)
    call resample_into(z4,z4out,11,axis=2,status=st)
    if (st/=resample_status_invalid_input .or. any(abs(z4out+301.0_dp)>0.0_dp)) &
        error stop 'atomic failure behavior'
    z4(3,7,2,2)=cmplx(x4(3,7,2,2),0.0_dp,dp)

    ! Noncontiguous input and output sections.
    scratch=-999.0_dp
    call resample_into(x4(1:3:2,:,:,:),scratch(1:6:3,:,:,:),11,axis=2,status=st)
    if (st/=resample_status_success) error stop 'strided resample'
    line=resample(x4(3,:,2,1),11)
    if (maxval(abs(scratch(4,:,2,1)-line))>tol) error stop 'strided values'
    if (any(scratch(2,:,:,:)/=-999.0_dp)) error stop 'strided overwrite'

    ! Complex taps and nonstandard axis with antireflect padding.
    do l=1,2
        do k=1,2
            do j=1,7
                do i=1,2
                    z5(i,1,j,k,l)=cmplx(x4(i,j,k,l),0.2_dp*real(i+j,dp),dp)
                end do
            end do
        end do
    end do
    call resample_poly_into(z5,z5out,3,2,taps=hc,axis=3,padtype='antireflect',status=st)
    if (st/=resample_status_success) error stop 'complex poly rank5 status'
    cline=resample_poly(z5(1,1,:,2,1),3,2,taps=hc,padtype='antireflect')
    if (maxval(abs(z5out(1,1,:,2,1)-cline))>tol) error stop 'complex poly rank5 values'
    call resample_poly_into(real(z5,dp),v5,3,2,taps=hc,axis=3,padtype='median',status=st)
    if (st/=resample_status_success) error stop 'real input complex taps'
    cline=resample_poly(real(z5(2,1,:,1,2),dp),3,2,taps=hc,padtype='median')
    if (maxval(abs(v5(2,1,:,1,2)-cline))>tol) error stop 'real complex poly values'

    ! All four upfirdn real/complex combinations.
    call upfirdn_into(h,x4,v4,up=2,down=1,axis=4,status=st)
    if (st/=resample_status_success) error stop 'upfirdn rank4 status'
    line=upfirdn(h,x4(2,3,1,:),up=2,down=1)
    if (maxval(abs(v4(2,3,1,:)-line))>tol) error stop 'upfirdn rank4'
    block
        complex(dp) :: cu4(3,7,2,5)
        call upfirdn_into(hc,x4,cu4,up=2,axis=4,status=st)
        if (st/=resample_status_success) error stop 'complex taps real upfirdn'
        cline=upfirdn(hc,x4(1,2,1,:),up=2)
        if (maxval(abs(cu4(1,2,1,:)-cline))>tol) error stop 'complex taps upfirdn values'
        call upfirdn_into(h,z4,cu4,up=2,axis=4,status=st)
        if (st/=resample_status_success) error stop 'real taps complex upfirdn'
        call upfirdn_into(hc,z4,cu4,up=2,axis=4,status=st)
        if (st/=resample_status_success) error stop 'complex taps complex upfirdn'
    end block

    do l=1,2
        do k=1,2
            do j=1,8
                do i=1,2
                    x6(i,k,1,j,1,l)=sin(0.12_dp*real(i+k+j+l,dp))
                    x6(i,k,2,j,1,l)=x6(i,k,1,j,1,l)+0.3_dp
                end do
            end do
        end do
    end do
    call decimate_into(x6,y6,2,n=10,ftype='fir',axis=4,status=st)
    if (st/=resample_status_success) error stop 'decimate rank6'
    line=decimate(x6(2,1,2,:,1,1),2,n=10,ftype='fir')
    if (maxval(abs(y6(2,1,2,:,1,1)-line))>tol) error stop 'decimate rank6 values'
    do l=1,2
        do k=1,2
            do j=1,7
                x8(1,1,k,j,1,1,1,l)=sin(0.17_dp*real(j+k+l,dp))
            end do
        end do
    end do
    call resample_poly_into(x8,y8,3,2,taps=h,axis=4,status=st)
    if (st/=resample_status_success) error stop 'poly rank8'
    line=resample_poly(x8(1,1,2,:,1,1,1,1),3,2,taps=h)
    if (maxval(abs(y8(1,1,2,:,1,1,1,1)-line))>tol) error stop 'poly rank8 values'
    do k=1,2
        do j=1,7
            x14(1,1,1,1,1,1,1,1,1,1,1,1,k,j)=sin(0.07_dp*real(k+j,dp))
            x15(1,1,1,1,1,1,1,1,1,1,1,1,1,k,j)=x14(1,1,1,1,1,1,1,1,1,1,1,1,k,j)
        end do
    end do
    call resample_into(x14,y14,11,axis=-1,status=st)
    if (st/=resample_status_success) error stop 'resample rank14'
    line=resample(x14(1,1,1,1,1,1,1,1,1,1,1,1,2,:),11)
    if (maxval(abs(y14(1,1,1,1,1,1,1,1,1,1,1,1,2,:)-line))>tol) error stop 'rank14 values'
    call upfirdn_into(h,x15,y15,up=2,down=1,status=st)
    if (st/=resample_status_success) error stop 'upfirdn rank15'
    line=upfirdn(h,x15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,:),up=2)
    if (maxval(abs(y15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,:)-line))>tol) error stop 'rank15 values'
    ! Default Kaiser taps, complex constant padding, and default IIR decimate.
    block
        real(dp) :: default_poly(3,11,2,2), xlong(2,48,1,1), ylong(2,24,1,1)
        real(dp) :: wrongrank(3,11,2)
        complex(dp) :: complex_poly(3,11,2,2)
        integer :: n
        call resample_poly_into(x4,default_poly,3,2,axis=2,status=st)
        if (st/=resample_status_success) error stop 'Kaiser default status'
        line=resample_poly(x4(1,:,2,2),3,2)
        if (maxval(abs(default_poly(1,:,2,2)-line))>tol) error stop 'Kaiser default values'
        call resample_poly_into(z4,complex_poly,3,2,taps=h, &
            padtype='constant',cval=cmplx(0.25_dp,-0.15_dp,dp),axis=2,status=st)
        if (st/=resample_status_success) error stop 'complex constant padding status'
        cline=resample_poly(z4(2,:,2,1),3,2,taps=h, &
            padtype='constant',cval=cmplx(0.25_dp,-0.15_dp,dp))
        if (maxval(abs(complex_poly(2,:,2,1)-cline))>tol) error stop 'complex padding values'
        do n=1,48
            xlong(1,n,1,1)=sin(0.11_dp*real(n,dp))
            xlong(2,n,1,1)=cos(0.07_dp*real(n,dp))
        end do
        call decimate_into(xlong,ylong,2,axis=2,status=st)
        if (st/=resample_status_success) error stop 'default IIR decimation'
        line=decimate(xlong(2,:,1,1),2)
        if (maxval(abs(ylong(2,:,1,1)-line))>tol) error stop 'default IIR values'
        wrongrank=-777.0_dp
        call resample_into(x4,wrongrank,11,axis=2,status=st)
        if (st/=resample_status_invalid_input .or. any(wrongrank/=-777.0_dp)) &
            error stop 'rank mismatch no mutation'
    end block
    print '(a)', 'test_signal_multirate_assumed_rank passed'
end program test_signal_multirate_assumed_rank
