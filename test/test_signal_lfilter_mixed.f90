! SPDX-License-Identifier: MIT
! Test real/complex input and coefficient promotion for IIR filtering.
program test_signal_lfilter_mixed
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfilter, &
        lfilter_status_success, lfilter_status_invalid_input
    implicit none
    real(dp) :: b(3), a(2), x1(16), x2(11,8), x3(5,8,4)
    complex(dp) :: cb(3), ca(2), cx1(16), cx2(11,8), cx3(5,8,4)
    complex(dp) :: zi2(11,2), zi3(5,2,4)
    complex(dp), allocatable :: y1(:), y2(:,:), y3(:,:,:), zf2(:,:), zf3(:,:,:)
    integer :: i, j, k, code

    b = [0.5_dp, -0.25_dp, 0.1_dp]
    a = [1.0_dp, -0.6_dp]
    cb = cmplx(b,0.05_dp*b,dp)
    ca = cmplx(a,-0.08_dp*a,dp)
    do i = 1, 16
        x1(i) = sin(0.2_dp*i)
    end do
    do j = 1, 8
        do i = 1, 11
            x2(i,j) = sin(0.17_dp*(i+j))
        end do
    end do
    do k = 1, 4
        do j = 1, 8
            do i = 1, 5
                x3(i,j,k) = cos(0.09_dp*(i+j+k))
            end do
        end do
    end do
    cx1 = cmplx(x1,-0.5_dp*x1,dp)
    cx2 = cmplx(x2,0.2_dp*x2,dp)
    cx3 = cmplx(x3,0.1_dp*x3,dp)

    y1 = lfilter(b,a,cx1,status=code)
    if (code /= lfilter_status_success) error stop 'rank1 real coeffs/complex signal status'
    call check1(y1,lfilter(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),cx1))
    call check1(lfilter(b,ca,x1),lfilter(cmplx(b,0.0_dp,dp),ca,cmplx(x1,0.0_dp,dp)))
    call check1(lfilter(cb,a,x1),lfilter(cb,cmplx(a,0.0_dp,dp),cmplx(x1,0.0_dp,dp)))
    call check1(lfilter(cb,ca,x1),lfilter(cb,ca,cmplx(x1,0.0_dp,dp)))
    call check1(lfilter(cb,a,cx1),lfilter(cb,cmplx(a,0.0_dp,dp),cx1))
    call check1(lfilter(b,ca,cx1),lfilter(cmplx(b,0.0_dp,dp),ca,cx1))

    zi2 = cmplx(0.1_dp,-0.02_dp,dp)
    zi3 = cmplx(0.2_dp,0.01_dp,dp)
    y2 = lfilter(b,a,cx2,axis=2,zi=zi2,zf=zf2,status=code)
    if (code /= lfilter_status_success) error stop 'rank2 mixed status'
    call check2(y2,lfilter(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),cx2,axis=2,zi=zi2))
    call check2(lfilter(b,ca,x2,axis=1), &
        lfilter(cmplx(b,0.0_dp,dp),ca,cmplx(x2,0.0_dp,dp),axis=1))
    call check2(lfilter(cb,a,x2,axis=2), &
        lfilter(cb,cmplx(a,0.0_dp,dp),cmplx(x2,0.0_dp,dp),axis=2))
    call check2(lfilter(cb,ca,x2,axis=2),lfilter(cb,ca,cmplx(x2,0.0_dp,dp),axis=2))
    call check2(lfilter(b,ca,cx2,axis=1),lfilter(cmplx(b,0.0_dp,dp),ca,cx2,axis=1))
    call check2(lfilter(cb,a,cx2,axis=1),lfilter(cb,cmplx(a,0.0_dp,dp),cx2,axis=1))
    if (any(shape(zf2) /= [11,2])) error stop 'rank2 state output shape'

    y3 = lfilter(cb,a,x3,axis=2,zi=zi3,zf=zf3,status=code)
    if (code /= lfilter_status_success) error stop 'rank3 mixed status'
    call check3(y3,lfilter(cb,cmplx(a,0.0_dp,dp),cmplx(x3,0.0_dp,dp),axis=2,zi=zi3))
    call check3(lfilter(b,a,cx3,axis=1), &
        lfilter(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),cx3,axis=1))
    call check3(lfilter(b,ca,x3,axis=3), &
        lfilter(cmplx(b,0.0_dp,dp),ca,cmplx(x3,0.0_dp,dp),axis=3))
    call check3(lfilter(cb,ca,x3),lfilter(cb,ca,cmplx(x3,0.0_dp,dp)))
    call check3(lfilter(b,ca,cx3,axis=2),lfilter(cmplx(b,0.0_dp,dp),ca,cx3,axis=2))
    call check3(lfilter(cb,a,cx3,axis=2),lfilter(cb,cmplx(a,0.0_dp,dp),cx3,axis=2))
    if (any(shape(zf3) /= [5,2,4])) error stop 'rank3 state output shape'
    y3 = lfilter(b,ca,cx3,axis=0,status=code)
    if (code /= lfilter_status_invalid_input) error stop 'rank3 bad axis status'
    print '(a)', 'test_signal_lfilter_mixed passed'
contains
    subroutine check1(got,expected)
        complex(dp), intent(in) :: got(:) !! computed samples
        complex(dp), intent(in) :: expected(:) !! reference samples
        if (size(got) /= size(expected)) error stop 'rank1 size mismatch'
        if (any(abs(got-expected)>2.0e-13_dp)) error stop 'rank1 data mismatch'
    end subroutine check1

    subroutine check2(got,expected)
        complex(dp), intent(in) :: got(:,:) !! computed samples
        complex(dp), intent(in) :: expected(:,:) !! reference samples
        if (any(shape(got) /= shape(expected))) error stop 'rank2 shape mismatch'
        if (any(abs(got-expected)>2.0e-13_dp)) error stop 'rank2 data mismatch'
    end subroutine check2

    subroutine check3(got,expected)
        complex(dp), intent(in) :: got(:,:,:) !! computed samples
        complex(dp), intent(in) :: expected(:,:,:) !! reference samples
        if (any(shape(got) /= shape(expected))) error stop 'rank3 shape mismatch'
        if (any(abs(got-expected)>2.0e-13_dp)) error stop 'rank3 data mismatch'
    end subroutine check3
end program test_signal_lfilter_mixed
