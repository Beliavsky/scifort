! SPDX-License-Identifier: MIT
program test_ndimage_vectorized_mixed
    use scifort_kinds, only : dp
    use scifort_ndimage, only : vectorized_filter_nd, vectorized_filter_multi_nd, &
        ndimage_status_invalid_input
    implicit none
    real(dp) :: x(4,3), r(4,3), r2(4,3,2), vr(3,2), &
        vr2(3,2,2), noncontig(8,6), flat15(4,2)
    complex(dp) :: z(4,3), c(4,3), c2(4,3,2), vc(3,2), vc2(3,2,2)
    real(dp) :: x15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: i,j,code
    real(dp), parameter :: tol=1.0e-12_dp

    do j=1,3
        do i=1,4
            x(i,j)=0.4_dp*i+0.2_dp*j+0.1_dp*i*j
        end do
    end do
    z=cmplx(x,0.5_dp*x+0.1_dp,dp)
    call vectorized_filter_nd(x,rc_feature,c,status=code,sizes=[2,2], &
        mode='nearest',batch_size=3)
    call check(code==0,'real to complex scalar status')
    call vectorized_filter_multi_nd(x,rc_features,c2,status=code,sizes=[2,2], &
        mode='nearest',batch_size=2)
    call check(code==0,'real to complex multi status')
    call check(maxval(abs(c-c2(:,:,1)))<tol,'real to complex first channel')
    call check(abs(c2(4,3,2)-cmplx(sum(x(3:4,2:3)**2), &
        sum(x(3:4,2:3)),dp)/4.0_dp)<tol, &
        'second channel nearest boundary')

    call vectorized_filter_nd(z,cr_feature,r,status=code,sizes=[2,2], &
        mode='constant',cval=(-2.0_dp,0.5_dp),batch_size=2)
    call check(code==0,'complex to real scalar status')
    call vectorized_filter_multi_nd(z,cr_features,r2,status=code,sizes=[2,2], &
        mode='constant',cval=(-2.0_dp,0.5_dp),batch_size=4)
    call check(code==0,'complex to real multi status')
    call check(maxval(abs(r-r2(:,:,1)))<tol,'complex to real first channel')
    call check(all(r>=0.0_dp),'positive squared magnitudes')

    call vectorized_filter_nd(x,rc_feature,vc,status=code,sizes=[2,2],mode='valid')
    call check(code==0,'real to complex valid mode')
    call check(abs(vc(1,1)-cmplx(sum(x(1:2,1:2))/4.0_dp, &
        sum(x(1:2,1:2)**2)/4.0_dp,dp))<tol,'real to complex valid reference')
    call vectorized_filter_multi_nd(x,rc_features,vc2,status=code, &
        sizes=[2,2],mode='valid',batch_size=1)
    call check(code==0,'real to complex multi valid mode')
    call check(maxval(abs(vc-vc2(:,:,1)))<tol,'valid multi real to complex')
    call vectorized_filter_nd(z,cr_feature,vr,status=code,sizes=[2,2],mode='valid')
    call check(code==0,'complex to real valid mode')
    call check(abs(vr(1,1)-sum(abs(z(1:2,1:2))**2)/4.0_dp)<tol, &
        'complex to real valid reference')
    call vectorized_filter_multi_nd(z,cr_features,vr2,status=code, &
        sizes=[2,2],mode='valid')
    call check(code==0,'complex to real multi valid mode')
    call check(maxval(abs(vr-vr2(:,:,1)))<tol,'valid multi complex to real')

    noncontig=-77.0_dp
    call vectorized_filter_nd(z,cr_feature,noncontig(1:8:2,1:6:2), &
        status=code,sizes=[2,2],mode='constant',cval=(-2.0_dp,0.5_dp))
    call check(code==0,'complex to real strided destination')
    call check(maxval(abs(noncontig(1:8:2,1:6:2)-r))<tol,'strided results')
    call check(all(noncontig(2:8:2,:)==-77.0_dp),'untouched strided rows')

    x15=1.0_dp
    x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=3.0_dp
    x15(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=5.0_dp
    x15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=7.0_dp
    call vectorized_filter_multi_nd(cmplx(x15,0.0_dp,dp),cr_features,flat15, &
        status=code,sizes=[2],axes=[2],mode='valid')
    call check(code==0,'rank15 complex to real flattened output')
    call check(abs(flat15(1,1)-13.0_dp)<tol .and. abs(flat15(2,1)-29.0_dp)<tol, &
        'rank15 complex to real values')

    call vectorized_filter_nd(x,rc_feature,c,status=code,sizes=[2,2],batch_size=0)
    call check(code==ndimage_status_invalid_input,'invalid mixed batch')
    call vectorized_filter_multi_nd(z,cr_features,r2,status=code, &
        sizes=[2,2],axes=[1,1])
    call check(code==ndimage_status_invalid_input,'invalid mixed axes')
    print '(a)', 'PASS: mixed-type scalar and multi-result vectorized filters'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! assertion to check
        character(len=*), intent(in) :: message !! assertion description
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
    subroutine rc_feature(w,v)
        real(dp), intent(in) :: w(:,:) !! real neighborhoods
        complex(dp), intent(out) :: v(:) !! means and mean squares
        v=cmplx(sum(w,dim=1),sum(w*w,dim=1),dp)/real(size(w,1),dp)
    end subroutine rc_feature
    subroutine rc_features(w,v)
        real(dp), intent(in) :: w(:,:) !! real neighborhoods
        complex(dp), intent(out) :: v(:,:) !! channels by neighborhood
        v(1,:)=cmplx(sum(w,dim=1),sum(w*w,dim=1),dp)/real(size(w,1),dp)
        v(2,:)=cmplx(sum(w*w,dim=1),sum(w,dim=1),dp)/real(size(w,1),dp)
    end subroutine rc_features
    subroutine cr_feature(w,v)
        complex(dp), intent(in) :: w(:,:) !! complex neighborhoods
        real(dp), intent(out) :: v(:) !! mean squared magnitude
        v=sum(abs(w)**2,dim=1)/real(size(w,1),dp)
    end subroutine cr_feature
    subroutine cr_features(w,v)
        complex(dp), intent(in) :: w(:,:) !! complex neighborhoods
        real(dp), intent(out) :: v(:,:) !! channels by neighborhood
        v(1,:)=sum(abs(w)**2,dim=1)/real(size(w,1),dp)
        v(2,:)=sum(real(w,dp),dim=1)/real(size(w,1),dp)
    end subroutine cr_features
end program test_ndimage_vectorized_mixed
