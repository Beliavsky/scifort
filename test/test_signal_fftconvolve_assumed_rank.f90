! SPDX-License-Identifier: MIT
! Regression for caller-output assumed-rank FFT convolution.
program test_signal_fftconvolve_assumed_rank
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : fftconvolve, fftconvolve_nd
    use scifort_signal_convolve, only : convolve, convolution_status_success, &
        convolution_status_invalid_input
    implicit none
    real(dp) :: a1(3), b1(2), y1(4), ysame(3), yvalid(2), bad(3)
    real(dp) :: a2(3,2), b2(2,2), y2(4,3), r2(4,3)
    real(dp) :: a4(3,2,2,2), b4(2,1,1,1), y4(4,2,2,2)
    real(dp) :: y4same(3,2,2,2), y4valid(2,2,2,2)
    real(dp) :: a5(3,2,3,1,2), b5(2,1,2,1,1), y5(4,2,4,1,2)
    real(dp) :: a15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: b15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: y15(4,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: aempty(0,3), yempty(1,3), incompatible(2,3), unbroadcastable(2,3,1,1)
    real(dp) :: along(127,2), blong(35,1), ylong(161,2)
    real(dp) :: k15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: y15axis(3,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: y5valid(2,2,2,1,2)
    complex(dp) :: ac3(3,1,2), bc3(2,2,2), yc3(4,2,3), refc3(4,2,3)
    complex(dp) :: cy1(4)
    real(dp) :: reference(4), expected
    integer :: code, i, j, k, l, u, v

    a1 = [1.0_dp,2.0_dp,-1.0_dp]
    b1 = [0.5_dp,-0.25_dp]
    call fftconvolve_nd(a1,b1,y1,status=code)
    if (code /= convolution_status_success) error stop 'rank1 real status'
    if (maxval(abs(y1-convolve(a1,b1))) > 1.e-12_dp) error stop 'rank1 values'
    call fftconvolve_nd(a1,b1,ysame,mode='same',status=code)
    if (code /= convolution_status_success) error stop 'rank1 same status'
    if (maxval(abs(ysame-fftconvolve(a1,b1,'same'))) > 1.e-12_dp) &
        error stop 'rank1 same values'
    call fftconvolve_nd(a1,b1,yvalid,mode='valid',status=code)
    if (code /= convolution_status_success) error stop 'rank1 valid status'
    if (maxval(abs(yvalid-fftconvolve(a1,b1,'valid'))) > 1.e-12_dp) &
        error stop 'rank1 valid values'

    do j = 1,2
        do i = 1,3
            a2(i,j) = sin(real(i+j,dp))
        end do
    end do
    b2 = reshape([0.5_dp,0.25_dp,-0.2_dp,0.1_dp],[2,2])
    call fftconvolve_nd(a2,b2,y2,status=code)
    if (code /= convolution_status_success) error stop 'rank2 status'
    r2 = fftconvolve(a2,b2)
    if (maxval(abs(y2-r2)) > 1.e-12_dp) error stop 'rank2 values'

    ! Longer FFTs, padded to a length greater than 128, retain accuracy.
    along = reshape([(sin(0.17_dp*i),i=1,size(along))],shape(along))
    blong(:,1) = [(cos(0.03_dp*i),i=1,size(blong))]
    call fftconvolve_nd(along,blong,ylong,axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'long FFT status'
    if (maxval(abs(ylong(:,2)-convolve(along(:,2),blong(:,1)))) &
        > 2.e-12_dp) error stop 'long FFT values'

    ac3 = reshape([(cmplx(0.1_dp*i,0.03_dp*i,dp),i=1,size(ac3))],shape(ac3))
    bc3 = reshape([(cmplx(0.02_dp*i,-0.013_dp*i,dp),i=1,size(bc3))],shape(bc3))
    call fftconvolve_nd(ac3,bc3,yc3,axes=[1,3],status=code)
    if (code /= convolution_status_success) error stop 'rank3 complex status'
    refc3 = fftconvolve(ac3,bc3,axes=[1,3])
    if (maxval(abs(yc3-refc3)) > 1.e-12_dp) error stop 'rank3 complex values'

    a4 = reshape([(0.08_dp*i + 0.01_dp*sin(real(i,dp)),i=1,size(a4))],shape(a4))
    b4(:,1,1,1) = [0.4_dp,-0.15_dp]
    call fftconvolve_nd(a4,b4,y4,axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'rank4 status'
    do l = 1,2
        do k = 1,2
            do j = 1,2
                reference = convolve(a4(:,j,k,l),b4(:,1,1,1))
                if (maxval(abs(y4(:,j,k,l)-reference)) > 1.e-12_dp) &
                    error stop 'rank4 axis1 values'
            end do
        end do
    end do
    call fftconvolve_nd(a4,b4,y4same,mode='same',axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'rank4 same status'
    if (maxval(abs(y4same-y4(1:3,:,:,:))) > 1.e-12_dp) error stop 'rank4 same values'
    call fftconvolve_nd(a4,b4,y4valid,mode='valid',axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'rank4 valid status'
    if (maxval(abs(y4valid-y4(2:3,:,:,:))) > 1.e-12_dp) error stop 'rank4 valid values'

    a5 = reshape([(0.02_dp*i-0.25_dp*cos(0.13_dp*i),i=1,size(a5))],shape(a5))
    b5 = reshape([(0.09_dp*i-0.17_dp,i=1,size(b5))],shape(b5))
    call fftconvolve_nd(a5,b5,y5,axes=[1,3],status=code)
    if (code /= convolution_status_success) error stop 'rank5 status'
    do l = 1,2
        do k = 1,4
            do j = 1,2
                do i = 1,4
                    expected = 0.0_dp
                    do v = 1,3
                        do u = 1,3
                            if (i-u+1 < 1 .or. i-u+1 > 2) cycle
                            if (k-v+1 < 1 .or. k-v+1 > 2) cycle
                            expected = expected+a5(u,j,v,1,l)*b5(i-u+1,1,k-v+1,1,1)
                        end do
                    end do
                    if (abs(expected-y5(i,j,k,1,l)) > 1.e-12_dp) &
                        error stop 'rank5 direct convolution'
                end do
            end do
        end do
    end do

    call fftconvolve_nd(a5,b5,y5valid,mode='valid',axes=[1,3],status=code)
    if (code /= convolution_status_success) error stop 'rank5 valid status'
    if (maxval(abs(y5valid-y5(2:3,:,2:3,:,:))) > 1.e-12_dp) &
        error stop 'rank5 valid values'

    a15 = reshape([(sin(0.3_dp*i),i=1,size(a15))],shape(a15))
    b15 = reshape([0.75_dp,-0.5_dp],shape(b15))
    call fftconvolve_nd(a15,b15,y15,axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'rank15 status'
    do l = 1,2
        reference = convolve(a15(:,1,1,1,1,1,1,1,1,1,1,1,1,1,l), &
            b15(:,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
        if (maxval(abs(y15(:,1,1,1,1,1,1,1,1,1,1,1,1,1,l)-reference)) &
            > 1.e-12_dp) error stop 'rank15 values'
    end do

    ! The final axis is also valid at the maximum Fortran 2018 rank.
    k15 = reshape([0.25_dp,0.75_dp],shape(k15))
    call fftconvolve_nd(a15,k15,y15axis,axes=[15],status=code)
    if (code /= convolution_status_success) error stop 'rank15 axis15 status'
    reference(1:3) = 0.0_dp
    block
        real(dp) :: axis_ref(3)
        axis_ref = convolve(a15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,:), &
            k15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:))
        if (maxval(abs(y15axis(2,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-axis_ref)) &
            > 1.e-12_dp) error stop 'rank15 axis15 values'
    end block

    ! Real/complex combinations (including rank-1 combinations without a
    ! legacy mixed-type function) must use complex outputs.
    call fftconvolve_nd(a1,cmplx(b1,0.1_dp*b1,dp),cy1,status=code)
    if (code /= convolution_status_success) error stop 'real/complex status'
    if (maxval(abs(cy1 - &
        [cmplx(y1(1),0.1_dp*y1(1),dp),cmplx(y1(2),0.1_dp*y1(2),dp), &
         cmplx(y1(3),0.1_dp*y1(3),dp),cmplx(y1(4),0.1_dp*y1(4),dp)])) &
        > 1.e-12_dp) error stop 'real/complex values'
    call fftconvolve_nd(cmplx(a1,0.3_dp*a1,dp),b1,cy1,status=code)
    if (code /= convolution_status_success) error stop 'complex/real status'
    if (maxval(abs(cy1-cmplx(y1,0.3_dp*y1,dp))) > 1.e-12_dp) &
        error stop 'complex/real values'
    call fftconvolve_nd(cmplx(a1,0.3_dp*a1,dp),cmplx(b1,0.1_dp*b1,dp), &
        cy1,status=code)
    if (code /= convolution_status_success) error stop 'complex/complex status'
    if (maxval(abs(cy1-cmplx(0.97_dp*y1,0.4_dp*y1,dp))) > 1.e-12_dp) &
        error stop 'complex/complex values'

    ! The CONTIGUOUS dummies perform Fortran-defined copy-in/out for sections.
    y4 = -777.0_dp
    call fftconvolve_nd(a4(:,1:2:2,:,:),b4,y4(:,1:2:2,:,:),axes=[1],status=code)
    if (code /= convolution_status_success) error stop 'strided output status'
    if (any(y4(:,2,:,:) /= -777.0_dp)) error stop 'strided output clobber'
    reference = convolve(a4(:,1,2,2),b4(:,1,1,1))
    if (maxval(abs(y4(:,1,2,2)-reference)) > 1.e-12_dp) &
        error stop 'strided output values'

    incompatible = 1.0_dp
    call fftconvolve_nd(a2,incompatible,y2,mode='valid',status=code)
    if (code /= convolution_status_invalid_input) error stop 'mixed valid ordering'
    unbroadcastable = 1.0_dp
    call fftconvolve_nd(a4,unbroadcastable,y4,axes=[1],status=code)
    if (code /= convolution_status_invalid_input) error stop 'invalid broadcast'
    call fftconvolve_nd(a1,b1,bad,status=code)
    if (code /= convolution_status_invalid_input) error stop 'wrong output status'
    if (.not. all(ieee_is_nan(bad))) error stop 'wrong output NaN'
    call fftconvolve_nd(a1,b1,y1,axes=[1,1],status=code)
    if (code /= convolution_status_invalid_input) error stop 'duplicate axes'
    call fftconvolve_nd(a1,b1,y1,axes=[2],status=code)
    if (code /= convolution_status_invalid_input) error stop 'invalid axis'
    call fftconvolve_nd(a1,b1,y1,axes=[integer::],status=code)
    if (code /= convolution_status_invalid_input) error stop 'empty axes'
    call fftconvolve_nd(a1,b1,y1,mode='nonsense',status=code)
    if (code /= convolution_status_invalid_input) error stop 'invalid mode'
    call fftconvolve_nd(a1,b2,y1,status=code)
    if (code /= convolution_status_invalid_input) error stop 'invalid rank'
    call fftconvolve_nd(aempty,b2,yempty,status=code)
    if (code /= convolution_status_invalid_input) error stop 'empty input'
    block
        real(dp) :: nan_a(3)
        nan_a = a1
        nan_a(2) = ieee_value(0.0_dp,ieee_quiet_nan)
        call fftconvolve_nd(nan_a,b1,y1,status=code)
        if (code /= convolution_status_invalid_input) error stop 'nonfinite input'
    end block
    print '(a)', 'test_signal_fftconvolve_assumed_rank passed'
end program test_signal_fftconvolve_assumed_rank
