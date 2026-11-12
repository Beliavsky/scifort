! SPDX-License-Identifier: MIT
! Multi-channel batched neighborhood filtering, valid mode, rank 15 and guards.
program test_ndimage_vectorized_multi
    use scifort_kinds, only : dp
    use scifort_ndimage, only : vectorized_filter_multi_nd, vectorized_filter_nd, &
        ndimage_status_too_expensive
    implicit none
    real(dp) :: x(4,3),y(4,3,2),flat(12,2),strided(8,6,2),baseline(4,3)
    real(dp) :: valid_y(3,2,2),invalid_y(4,2,2),empty_y(4,3,0)
    real(dp) :: x15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1),flat15(4,2)
    complex(dp) :: z(4,3),zy(4,3,2)
    integer :: i,j,code,calls,seen
    character(len=16) :: mixed_modes(2)
    real(dp), parameter :: tol=1.0e-12_dp

    do j=1,3
        do i=1,4
            x(i,j)=real(i*i+3*j,dp)/7.0_dp
        end do
    end do
    z=cmplx(x,0.4_dp*x+0.2_dp,dp)
    calls=0
    seen=0
    call vectorized_filter_multi_nd(x,real_features,y,status=code,sizes=[2,3], &
        mode='nearest',batch_size=5)
    call check(code==0,'multichannel real status')
    call check(calls==3 .and. seen==12,'batch grouping is preserved')
    call vectorized_filter_nd(x,scalar_mean,baseline,status=code,sizes=[2,3],mode='nearest')
    call check(code==0,'scalar mean reference')
    call check(maxval(abs(y(:,:,1)-baseline))<tol,'first real feature matches mean')
    call check(all(y(:,:,2)>=y(:,:,1)**2-tol),'mean square is at least squared mean')

    call vectorized_filter_multi_nd(x,real_features,flat,status=code,sizes=[2,3], &
        mode='nearest',batch_size=1)
    call check(code==0,'flat output status')
    call check(maxval(abs(flat-reshape(y,[12,2])))<tol,'flat and shaped output equal')
    strided=-999.0_dp
    call vectorized_filter_multi_nd(x,real_features,strided(1:8:2,2:6:2,:),status=code, &
        sizes=[2,3],mode='nearest',batch_size=2)
    call check(code==0,'strided output status')
    call check(maxval(abs(strided(1:8:2,2:6:2,:)-y))<tol,'noncontiguous output values')
    call check(all(strided(:,1,:)==-999.0_dp),'untouched output first column')
    call check(all(strided(2:8:2,:,:)==-999.0_dp),'untouched output even rows')

    call vectorized_filter_multi_nd(x,real_features,valid_y,status=code,sizes=[2,2],mode='valid')
    call check(code==0,'valid mode shape')
    do j=1,2
        do i=1,3
            call check(abs(valid_y(i,j,1)-sum(x(i:i+1,j:j+1))/4.0_dp)<tol,'valid means')
            call check(abs(valid_y(i,j,2)-sum(x(i:i+1,j:j+1)**2)/4.0_dp)<tol,'valid squares')
        end do
    end do
    mixed_modes=[character(len=16) :: 'wrap','constant']
    call vectorized_filter_multi_nd(x,real_features,y,status=code,sizes=[2,2], &
        modes=mixed_modes,cval=-2.0_dp)
    call check(code==0,'mixed-mode multichannel filtering')
    call check(abs(y(1,1,1)-(-4.0_dp+x(4,1)+x(1,1))/4.0_dp)<tol, &
        'mixed-mode real mean corner')
    call check(abs(y(1,1,2)-(8.0_dp+x(4,1)**2+x(1,1)**2)/4.0_dp)<tol, &
        'mixed-mode square corner')
    call vectorized_filter_multi_nd(x,real_features,y,status=code, &
        sizes=[integer::],axes=[integer::])
    call check(code==0,'empty-axis multichannel identity')
    call check(maxval(abs(y(:,:,1)-x))<tol .and. &
        maxval(abs(y(:,:,2)-x*x))<tol,'empty-axis channels')
    call vectorized_filter_multi_nd(x,real_features,invalid_y,status=code,sizes=[2,2],mode='valid')
    call check(code/=0,'invalid channel-last output shape')
    call vectorized_filter_multi_nd(x,real_features,y,status=code,sizes=[2,3],batch_size=0)
    call check(code/=0,'invalid batch size')
    call vectorized_filter_multi_nd(x,real_features,y,status=code,sizes=[2,3],axes=[1,1])
    call check(code/=0,'repeated axes')
    call vectorized_filter_multi_nd(x,real_features,valid_y,status=code,sizes=[2,2], &
        mode='valid',origin=[1,0])
    call check(code/=0,'shifted origin with valid mode rejected')
    call vectorized_filter_multi_nd(x,real_features,empty_y,status=code,sizes=[2,2])
    call check(code/=0,'zero channels rejected')
    call vectorized_filter_multi_nd(x,real_features,y,status=code,sizes=[60000,1])
    call check(code==ndimage_status_too_expensive,'oversized stencil guarded')

    call vectorized_filter_multi_nd(z,complex_features,zy,status=code,sizes=[2,2], &
        mode='constant',cval=cmplx(-1.0_dp,0.5_dp,dp),batch_size=3)
    call check(code==0,'complex multichannel callback')
    call check(abs(zy(1,1,1)-sum([cmplx(-1.0_dp,0.5_dp,dp), &
        cmplx(-1.0_dp,0.5_dp,dp),cmplx(-1.0_dp,0.5_dp,dp),z(1,1)])/4.0_dp)<tol, &
        'complex border mean')

    x15=1.0_dp
    x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=3.0_dp
    x15(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=5.0_dp
    x15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=7.0_dp
    call vectorized_filter_multi_nd(x15,real_features,flat15,status=code, &
        sizes=[2],axes=[2],mode='valid',batch_size=1)
    call check(code==0,'rank15 flat output')
    call check(abs(flat15(1,1)-3.0_dp)<tol .and. abs(flat15(2,1)-5.0_dp)<tol, &
        'rank15 first channel')
    call check(abs(flat15(1,2)-13.0_dp)<tol .and. abs(flat15(2,2)-29.0_dp)<tol, &
        'rank15 second channel')
    print '(a)', 'PASS: vectorized multichannel callbacks, valid, complex, rank15'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! condition
        character(len=*), intent(in) :: message !! diagnostic
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
    subroutine real_features(w,r)
        real(dp), intent(in) :: w(:,:) !! selected windows
        real(dp), intent(out) :: r(:,:) !! mean and mean square, channel by window
        calls=calls+1
        seen=seen+size(w,2)
        r(1,:)=sum(w,dim=1)/real(size(w,1),dp)
        r(2,:)=sum(w*w,dim=1)/real(size(w,1),dp)
    end subroutine real_features
    subroutine scalar_mean(w,r)
        real(dp), intent(in) :: w(:,:) !! windows
        real(dp), intent(out) :: r(:) !! means
        r=sum(w,dim=1)/real(size(w,1),dp)
    end subroutine scalar_mean
    subroutine complex_features(w,r)
        complex(dp), intent(in) :: w(:,:) !! complex windows
        complex(dp), intent(out) :: r(:,:) !! complex mean and mean square
        r(1,:)=sum(w,dim=1)/real(size(w,1),dp)
        r(2,:)=sum(w*w,dim=1)/real(size(w,1),dp)
    end subroutine complex_features
end program test_ndimage_vectorized_multi
