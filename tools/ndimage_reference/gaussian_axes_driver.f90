! SPDX-License-Identifier: MIT
! Gaussian selected-axis reference output, including validation and rank-15.
program gaussian_axes_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : gaussian_filter_nd, &
        ndimage_status_invalid_input,ndimage_status_success
    implicit none
    integer :: i,j,k,unit,status
    real(dp) :: x(4,5,3),y(4,5,3),x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: z(4,5,3),out(4,5,3)
    do k=1,3
        do j=1,5
            do i=1,4
                x(i,j,k)=sin(0.31_dp*real(i+2*j+4*k,dp))+0.3_dp*real(i-j,dp)
            end do
        end do
    end do
    z=cmplx(x,0.73_dp*x-0.2_dp,dp)
    open(newunit=unit,file='gaussian_axes_output.bin',status='replace', &
        action='write',access='stream',form='unformatted')
    call gaussian_filter_nd(x,[0.7_dp,1.2_dp],y,order=[1,0],radius=[3,2], &
        modes=[character(len=16) :: 'constant','mirror'],cval=-1.3_dp,axes=[3,1],status=status)
    if (status/=ndimage_status_success) error stop 'Gaussian selected vector'
    write(unit)y
    call gaussian_filter_nd(x,1.1_dp,y,order=2,mode='nearest',axes=[2],status=status)
    if (status/=ndimage_status_success) error stop 'Gaussian selected scalar'
    write(unit)y
    call gaussian_filter_nd(x,[1.0_dp,0.4_dp],y, &
        modes=[character(len=16) :: 'wrap','reflect'],axes=[3,1],status=status)
    if (status/=ndimage_status_success) error stop 'Gaussian mode sequence'
    write(unit)y
    call gaussian_filter_nd(x,2.0_dp,y,axes=[integer ::],status=status)
    if (status/=ndimage_status_success) error stop 'Gaussian no axes'
    write(unit)y
    call gaussian_filter_nd(z,[0.8_dp,1.1_dp],out,order=[1,2], &
        modes=[character(len=16) :: 'constant','nearest'], &
        cval=cmplx(0.4_dp,-0.8_dp,dp),axes=[3,1],status=status)
    if (status/=ndimage_status_success) error stop 'complex Gaussian selected'
    write(unit)out
    call gaussian_filter_nd(z,1.2_dp,out,mode='mirror',axes=[1,3],status=status)
    if (status/=ndimage_status_success) error stop 'complex Gaussian scalar'
    write(unit)out
    close(unit)

    call gaussian_filter_nd(x,[0.8_dp,1.2_dp],y,axes=[3,3],status=status)
    if (status/=ndimage_status_invalid_input) error stop 'duplicate axes accepted'
    call gaussian_filter_nd(x,[0.8_dp],y,axes=[3,1],status=status)
    if (status/=ndimage_status_invalid_input) error stop 'bad sigma length accepted'
    call gaussian_filter_nd(x,[0.8_dp,1.2_dp],y,order=[1],axes=[3,1],status=status)
    if (status/=ndimage_status_invalid_input) error stop 'bad order length accepted'
    call gaussian_filter_nd(x,[0.8_dp,1.2_dp],y,radius=[3],axes=[3,1],status=status)
    if (status/=ndimage_status_invalid_input) error stop 'bad radius length accepted'
    call gaussian_filter_nd(x,[0.8_dp,1.2_dp],y,mode='nearest', &
        modes=[character(len=16) :: 'reflect','mirror'],axes=[3,1],status=status)
    if (status/=ndimage_status_invalid_input) error stop 'conflicting modes accepted'
    x15=0.0_dp
    x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2)=2.0_dp
    call gaussian_filter_nd(x15,0.7_dp,y15,axes=[15],status=status)
    if (status/=ndimage_status_success) error stop 'Gaussian rank 15'
    if (.not.all(abs(y15)<2.01_dp)) error stop 'Gaussian rank 15 bounds'
    print '(a)', 'PASS: Gaussian axes reference driver'
end program gaussian_axes_driver
