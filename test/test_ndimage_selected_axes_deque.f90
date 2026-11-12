! SPDX-License-Identifier: MIT
! Native tests for O(n+width) extrema and selected-axis Gaussian/uniform/extrema.
program test_ndimage_selected_axes_deque
    use scifort_kinds, only : dp
    use scifort_ndimage, only : gaussian_filter_nd, gaussian_filter1d_nd, &
        minimum_filter_nd, maximum_filter_nd, minimum_filter1d_nd, &
        maximum_filter1d_nd, uniform_filter_nd, ndimage_status_success, &
        ndimage_status_invalid_input, ndimage_status_too_expensive
    implicit none
    real(dp) :: x(5,7,3),y(5,7,3),tmp(5,7,3),ref(5,7,3), &
        extended(10,14,3),x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,5), &
        y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,5)
    complex(dp) :: z(5,7,3),w(5,7,3)
    integer :: i,j,k,status,code

    do k=1,3
        do j=1,7
            do i=1,5
                x(i,j,k)=sin(.12_dp*i+.4_dp*j+.3_dp*k) + .2_dp*i - .17_dp*k
            end do
        end do
    end do
    z=cmplx(x,0.5_dp*x,dp)
    call maximum_filter_nd(x,[3,4],y,status=status,axes=[3,1], &
        modes=[character(len=16) :: 'nearest','mirror'])
    call maximum_filter1d_nd(x,3,3,tmp,status=code,mode='nearest')
    call maximum_filter1d_nd(tmp,4,1,ref,status=code,mode='mirror')
    if (status/=ndimage_status_success .or. code/=0 .or. any(y/=ref)) &
        error stop 'selected maximum must equal 1D passes'
    call minimum_filter_nd(x,4,y,status=status,mode='nearest',axes=[2])
    call minimum_filter1d_nd(x,4,2,ref,status=code,mode='nearest')
    if (status/=0 .or. code/=0 .or. any(y/=ref)) &
        error stop 'selected minimum equals axis filter'
    call maximum_filter_nd(x,[7,1,8],y,status=status,mode='wrap')
    call maximum_filter1d_nd(x,7,1,tmp,status=code,mode='wrap')
    call maximum_filter1d_nd(tmp,8,3,ref,status=code,mode='wrap')
    if (status/=0 .or. code/=0 .or. any(y/=ref)) &
        error stop 'rectangular max equals 1D composition'
    call gaussian_filter_nd(x,[.7_dp,1.1_dp],y,order=[1,0],status=status, &
        cval=.4_dp,modes=[character(len=16) :: 'constant','mirror'],axes=[3,1])
    call gaussian_filter1d_nd(x,.7_dp,3,tmp,order=1,status=code,mode='constant',cval=.4_dp)
    call gaussian_filter1d_nd(tmp,1.1_dp,1,ref,order=0,status=code,mode='mirror',cval=.4_dp)
    if (status/=0 .or. code/=0 .or. maxval(abs(y-ref))>3.e-13_dp) &
        error stop 'selected Gaussian derivative composition'
    call gaussian_filter_nd(z,.6_dp,w,status=status,axes=[integer ::])
    if (status/=0 .or. any(w/=z)) error stop 'zero selected axes identity'
    call uniform_filter_nd(x,3,y,status=status,axes=[integer ::])
    if (status/=0 .or. any(y/=x)) error stop 'empty uniform axes identity'
    call minimum_filter_nd(x,3,y,status=status,axes=[integer ::])
    if (status/=0 .or. any(y/=x)) error stop 'empty extrema axes identity'
    call minimum_filter1d_nd(x,50000,2,y,status=status,mode='constant',cval=-2.0_dp)
    if (status/=0 .or. any(y/=(-2.0_dp))) error stop 'very wide sliding minimum'
    call maximum_filter1d_nd(x,50000,2,y,status=status,mode='constant',cval=2.0_dp)
    if (status/=0 .or. any(y/=2.0_dp)) error stop 'very wide sliding maximum'
    call maximum_filter1d_nd(x,50001,2,y,status=status)
    if (status/=ndimage_status_too_expensive) error stop 'resource-limit regression'

    call gaussian_filter_nd(x,[.5_dp,.8_dp],y,status=status,axes=[1,1])
    if (status/=ndimage_status_invalid_input) error stop 'duplicate Gaussian axes'
    call uniform_filter_nd(x,[3],y,status=status,axes=[1,4])
    if (status/=ndimage_status_invalid_input) error stop 'invalid uniform axes'
    call maximum_filter_nd(x,[3,4],y,status=status,axes=[3,1,2])
    if (status/=ndimage_status_invalid_input) error stop 'sizes mismatch'
    call minimum_filter_nd(x,3,y,status=status,axes=[3],footprint=[.true.,.true.,.true.])
    if (status/=ndimage_status_invalid_input) error stop 'axes plus footprint must reject'
    call uniform_filter_nd(x,[3],y,status=status,axes=[2], &
        modes=[character(len=16) :: 'nearest','reflect'])
    if (status/=ndimage_status_invalid_input) error stop 'wrong number of selected modes'
    call gaussian_filter_nd(x,.7_dp,y,status=status,mode='wrap', &
        modes=[character(len=16) :: 'reflect'],axes=[1])
    if (status/=ndimage_status_invalid_input) error stop 'conflicting modes'
    call minimum_filter_nd(x,3,extended(1:10:2,1:14:2,:),status=status,axes=[2])
    call minimum_filter_nd(x,3,ref,status=code,axes=[2])
    if (status/=0 .or. code/=0 .or. any(extended(1:10:2,1:14:2,:)/=ref)) &
        error stop 'selected axes strided output'
    x15=0.0_dp
    x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,4)=8.0_dp
    call maximum_filter_nd(x15,3,y15,status=status,axes=[15])
    if (status/=0 .or. any(y15<0.0_dp) .or. maxval(y15)/=8.0_dp) &
        error stop 'rank-15 selected maximum'
    call gaussian_filter_nd(x15,.7_dp,y15,status=status,axes=[15])
    if (status/=0 .or. any(y15<0.0_dp)) error stop 'rank-15 selected Gaussian'
    print '(a)', 'PASS: selected axes and deque extrema'
end program test_ndimage_selected_axes_deque
