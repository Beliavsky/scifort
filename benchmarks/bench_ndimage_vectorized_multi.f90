! SPDX-License-Identifier: MIT
! Compare two separate scalar batches against one multichannel batch.
program bench_ndimage_vectorized_multi
    use scifort_kinds, only : dp
    use scifort_ndimage, only : vectorized_filter_nd, vectorized_filter_multi_nd
    implicit none
    integer, parameter :: nx=160,ny=128,repeats=4
    real(dp) :: x(nx,ny),a(nx,ny),b(nx,ny),pair(nx,ny,2)
    real(dp) :: t0,t1,t2,err
    integer :: i,j,k,code

    do j=1,ny
        do i=1,nx
            x(i,j)=sin(real(i,dp)*0.17_dp)+cos(real(j,dp)*0.13_dp)
        end do
    end do
    call cpu_time(t0)
    do k=1,repeats
        call vectorized_filter_nd(x,mean_callback,a,status=code,sizes=[7,7],mode='nearest')
        if (code/=0) error stop 'first scalar batch failed'
        call vectorized_filter_nd(x,square_callback,b,status=code,sizes=[7,7],mode='nearest')
        if (code/=0) error stop 'second scalar batch failed'
    end do
    call cpu_time(t1)
    do k=1,repeats
        call vectorized_filter_multi_nd(x,both_callback,pair,status=code,sizes=[7,7],mode='nearest')
        if (code/=0) error stop 'multichannel batch failed'
    end do
    call cpu_time(t2)
    err=max(maxval(abs(a-pair(:,:,1))),maxval(abs(b-pair(:,:,2))))
    if (err>1.0e-12_dp) error stop 'batch outputs differ'
    print '(a,i0,a,i0,a,i0)', 'Image ',nx,' x ',ny,', repeats ',repeats
    print '(a,f9.4)', 'Separate callbacks CPU seconds: ',t1-t0
    print '(a,f9.4)', 'Multi-result callback CPU seconds: ',t2-t1
    if (t2>t1) print '(a,f8.3)', 'Separate / multichannel CPU ratio: ',(t1-t0)/(t2-t1)
    print '(a,es12.3)', 'Maximum absolute output difference: ',err
contains
    subroutine mean_callback(w,r)
        real(dp), intent(in) :: w(:,:) !! windows
        real(dp), intent(out) :: r(:) !! mean
        r=sum(w,dim=1)/real(size(w,1),dp)
    end subroutine mean_callback
    subroutine square_callback(w,r)
        real(dp), intent(in) :: w(:,:) !! windows
        real(dp), intent(out) :: r(:) !! mean square
        r=sum(w*w,dim=1)/real(size(w,1),dp)
    end subroutine square_callback
    subroutine both_callback(w,r)
        real(dp), intent(in) :: w(:,:) !! windows
        real(dp), intent(out) :: r(:,:) !! channel by batch
        r(1,:)=sum(w,dim=1)/real(size(w,1),dp)
        r(2,:)=sum(w*w,dim=1)/real(size(w,1),dp)
    end subroutine both_callback
end program bench_ndimage_vectorized_multi
