! SPDX-License-Identifier: MIT
! Selected-axis order statistics and linear-time uniform windows.
program test_ndimage_order_axes_running
    use scifort_kinds, only : dp
    use scifort_ndimage, only : median_filter_nd, rank_filter_nd, percentile_filter_nd, &
        uniform_filter1d_nd, uniform_filter_nd, correlate1d_nd, &
        ndimage_status_success, ndimage_status_invalid_input, ndimage_status_too_expensive
    implicit none
    real(dp) :: x(4,5,3),y(4,5,3),tmp(4,5,3),ref(4,5,3)
    complex(dp) :: z(4,5,3),w(4,5,3),cz(4,5,3)
    real(dp) :: x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,5)
    real(dp) :: y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,5)
    real(dp) :: wide(4,3),big(8,10,3)
    logical :: fp(2,3)
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'reflect','mirror','nearest','wrap','grid-wrap','constant','grid-constant','grid-mirror']
    integer :: i,j,k,m,n,status,code
    real(dp), allocatable :: weights(:)

    do k=1,3
        do j=1,5
            do i=1,4
                x(i,j,k)=sin(.17_dp*i+.31_dp*j-.19_dp*k)+.05_dp*(i*i-2*j+k)
            end do
        end do
    end do
    z=cmplx(x,.4_dp*x+.2_dp,dp)
    fp=reshape([.true.,.false.,.true.,.true.,.false.,.true.],shape(fp))

    call median_filter_nd(x,3,y,status=status,axes=[integer ::])
    if (status/=0 .or. any(y/=x)) error stop 'empty selected median axes'
    call rank_filter_nd(x,-1,3,y,status=status,axes=[integer ::])
    if (status/=0 .or. any(y/=x)) error stop 'empty selected rank axes'
    call percentile_filter_nd(x,33.0_dp,3,y,status=status,axes=[integer ::])
    if (status/=0 .or. any(y/=x)) error stop 'empty selected percentile axes'

    call median_filter_nd(x,[3,4],y,status=status,axes=[1,1])
    if (status/=ndimage_status_invalid_input) error stop 'duplicate median axes'
    call rank_filter_nd(x,2,[2,4],y,status=status,axes=[1,4])
    if (status/=ndimage_status_invalid_input) error stop 'out-of-bounds rank axes'
    call percentile_filter_nd(x,55.0_dp,[3],y,status=status,axes=[1,3])
    if (status/=ndimage_status_invalid_input) error stop 'percentile size mismatch'
    call median_filter_nd(x,[2,3],y,status=status,axes=[3,1],origin=[0,1])
    if (status/=ndimage_status_invalid_input) error stop 'invalid canonical origin'
    call rank_filter_nd(x,4,[2,2],y,status=status,axes=[3,1])
    if (status/=ndimage_status_invalid_input) error stop 'rank too large'
    call percentile_filter_nd(x,-1.0_dp,3,y,status=status,axes=[2])
    if (status/=ndimage_status_invalid_input) error stop 'invalid percentile'
    call median_filter_nd(x,[2,3],y,status=status,axes=[3,1], &
        modes=[character(len=16) :: 'nearest'])
    if (status/=ndimage_status_invalid_input) error stop 'invalid modes count'
    call median_filter_nd(x,[2,3],y,status=status,axes=[1,3],mode='reflect', &
        modes=[character(len=16) :: 'nearest','mirror'])
    if (status/=ndimage_status_invalid_input) error stop 'conflicting modes'
    call rank_filter_nd(x,0,y=y,status=status,axes=[3,1],footprint=fp,origin=[0,1])
    if (status/=ndimage_status_invalid_input) error stop 'invalid footprint origin'
    call rank_filter_nd(x,0,y=y,status=status,axes=[3,1],footprint=fp,origin=[-1,1])
    if (status/=ndimage_status_invalid_input) error stop 'reversed footprint origin invalid'

    do m=1,size(modes)
        do n=2,9
            call uniform_filter1d_nd(x,n,2,y,status=status,mode=trim(modes(m)))
            allocate(weights(n))
            weights=1.0_dp/real(n,dp)
            call correlate1d_nd(x,weights,2,ref,status=code,mode=trim(modes(m)))
            deallocate(weights)
            if (status/=ndimage_status_success .or. code/=0 .or. maxval(abs(y-ref))>2.e-13_dp) &
                error stop 'uniform running sum differs from direct convolution'
        end do
    end do
    call uniform_filter_nd(x,[3,4],y,status=status,axes=[3,1], &
        modes=[character(len=16) :: 'constant','wrap'],cval=-1.25_dp)
    call uniform_filter1d_nd(x,3,3,tmp,status=code,mode='constant',cval=-1.25_dp)
    call uniform_filter1d_nd(tmp,4,1,ref,status=code,mode='wrap',cval=-1.25_dp)
    if (status/=0 .or. code/=0 .or. maxval(abs(y-ref))>1.e-13_dp) &
        error stop 'uniform selected axis order'
    call uniform_filter_nd(z,3,w,status=status,axes=[2],mode='nearest')
    call uniform_filter1d_nd(z,3,2,cz,status=code,mode='nearest')
    if (status/=0 .or. code/=0 .or. maxval(abs(w-cz))>2.e-13_dp) &
        error stop 'complex moving average'
    call uniform_filter1d_nd(x,50000,2,y,status=status,mode='constant',cval=3.5_dp)
    do k=1,3
        do i=1,4
            if (maxval(abs(y(i,:,k) - &
                (sum(x(i,:,k))+real(50000-5,dp)*3.5_dp)/50000.0_dp))>2.e-10_dp) &
                error stop 'wide constant moving mean'
        end do
    end do
    if (status/=0) error stop 'wide uniform status'
    call uniform_filter1d_nd(x,50001,2,y,status=status)
    if (status/=ndimage_status_too_expensive) error stop 'oversize uniform guard'
    call uniform_filter1d_nd(x,3,4,y,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid uniform axis'
    call uniform_filter1d_nd(x,4,2,y,status=status,origin=2)
    if (status/=ndimage_status_invalid_input) error stop 'invalid uniform origin'
    call uniform_filter_nd(x,[3,4],y,status=status,axes=[1,1])
    if (status/=ndimage_status_invalid_input) error stop 'invalid uniform axes duplicate'

    x15=reshape([3.0_dp,1.0_dp,5.0_dp,2.0_dp,4.0_dp],shape(x15))
    call rank_filter_nd(x15,-1,3,y15,status=status,axes=[15],mode='nearest')
    if (status/=0 .or. any(reshape(y15,[5])/=[3.0_dp,5.0_dp,5.0_dp,5.0_dp,4.0_dp])) &
        error stop 'rank-15 selected maximum'
    call uniform_filter_nd(x15,3,y15,status=status,axes=[15])
    if (status/=0 .or. any(abs(y15)>10.0_dp)) error stop 'rank-15 running uniform'
    big=0.0_dp
    call uniform_filter_nd(x,3,big(1:8:2,1:10:2,:),status=status,axes=[2])
    call uniform_filter1d_nd(x,3,2,ref,status=code)
    if (status/=0 .or. code/=0 .or. maxval(abs(big(1:8:2,1:10:2,:)-ref))>1.e-14_dp) &
        error stop 'noncontiguous uniform output'
    print '(a)', 'PASS: order-statistic selected axes and running uniform sums'
end program test_ndimage_order_axes_running
