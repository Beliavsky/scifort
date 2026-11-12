! SPDX-License-Identifier: MIT
! Reference array sequences for axis-selective filters and deque extrema.
program selected_axes_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : minimum_filter_nd, maximum_filter_nd, &
        minimum_filter1d_nd, maximum_filter1d_nd, uniform_filter_nd, &
        ndimage_status_success, ndimage_status_invalid_input
    implicit none
    integer :: i,j,k,status,unit
    real(dp) :: x(4,5,3),y(4,5,3)
    complex(dp) :: z(4,5,3),zout(4,5,3)
    character(len=16) :: modes(2)
    do k=1,3
        do j=1,5
            do i=1,4
                x(i,j,k)=sin(0.31_dp*real(i+2*j+4*k,dp))+0.3_dp*real(i-j,dp)
            end do
        end do
    end do
    modes=[character(len=16) :: 'wrap','mirror']
    open(newunit=unit,file='selected_axes_output.bin',status='replace', &
        access='stream',form='unformatted',action='write')
    call minimum_filter_nd(x,[2,4],y,status,modes=modes,axes=[3,1])
    if (status/=ndimage_status_success) error stop 'minimum axis sequence'
    write(unit)y
    call maximum_filter_nd(x,[2,4],y,status,modes=modes,axes=[3,1])
    if (status/=ndimage_status_success) error stop 'maximum axis sequence'
    write(unit)y
    call uniform_filter_nd(x,[2,4],y,status,modes=modes,axes=[3,1])
    if (status/=ndimage_status_success) error stop 'uniform axis sequence'
    write(unit)y
    call minimum_filter_nd(x,4,y,status,mode='constant',cval=-3.5_dp,axes=[2])
    if (status/=ndimage_status_success) error stop 'minimum scalar axis'
    write(unit)y
    call maximum_filter_nd(x,3,y,status,mode='nearest',axes=[1,3])
    if (status/=ndimage_status_success) error stop 'maximum scalar axes'
    write(unit)y
    call uniform_filter_nd(x,3,y,status,mode='reflect',axes=[1,3])
    if (status/=ndimage_status_success) error stop 'uniform scalar axes'
    write(unit)y
    call minimum_filter1d_nd(x,5000,axis=2,y=y,status=status,mode='grid-wrap',origin=1)
    if (status/=ndimage_status_success) error stop 'large minimum 1d'
    write(unit)y
    call maximum_filter1d_nd(x,5000,axis=2,y=y,status=status,mode='constant', &
        cval=2.5_dp,origin=-2500)
    if (status/=ndimage_status_success) error stop 'large maximum 1d'
    write(unit)y
    call maximum_filter_nd(x,[8,1,7],y,status,mode='nearest')
    if (status/=ndimage_status_success) error stop 'large rectangle'
    write(unit)y
    call minimum_filter_nd(x,[1,1,1],y,status,mode='constant',cval=-77.0_dp)
    if (status/=ndimage_status_success) error stop 'unit rectangle'
    write(unit)y
    call minimum_filter_nd(x,[4,3],y,status,axes=[3,1], &
        modes=[character(len=16) :: 'constant','nearest'],cval=-3.2_dp)
    if (status/=ndimage_status_success) error stop 'minimum mixed constant'
    write(unit)y
    call maximum_filter_nd(x,[4,3],y,status,axes=[3,1], &
        modes=[character(len=16) :: 'constant','mirror'],cval=3.2_dp)
    if (status/=ndimage_status_success) error stop 'maximum mixed constant'
    write(unit)y
    call uniform_filter_nd(x,[4,3],y,status,axes=[3,1], &
        modes=[character(len=16) :: 'constant','nearest'],cval=-3.2_dp)
    if (status/=ndimage_status_success) error stop 'uniform mixed constant'
    write(unit)y
    call uniform_filter_nd(x,4,y,status,axes=[integer ::])
    if (status/=ndimage_status_success) error stop 'uniform empty axes'
    write(unit)y
    call maximum_filter_nd(x,[3,3],y,status,axes=[1,1])
    if (status/=ndimage_status_invalid_input) error stop 'duplicate axes accepted'
    call minimum_filter_nd(x,[3],y,status,axes=[1,3])
    if (status/=ndimage_status_invalid_input) error stop 'mismatched sizes accepted'
    call uniform_filter_nd(x,[2],y,status,axes=[0])
    if (status/=ndimage_status_invalid_input) error stop 'zero axis accepted'
    call minimum_filter_nd(x,[2],y,status,axes=[2],origin=[0,0,0])
    if (status/=ndimage_status_invalid_input) error stop 'footprint origin with axes accepted'
    close(unit)
    z=cmplx(x,0.3_dp*x+0.4_dp,dp)
    open(newunit=unit,file='selected_axes_complex_output.bin',status='replace', &
        access='stream',form='unformatted',action='write')
    call uniform_filter_nd(z,[2,3],zout,status=status,axes=[3,1], &
        modes=[character(len=16) :: 'constant','wrap'], &
        cval=cmplx(1.25_dp,-0.5_dp,dp))
    if (status/=ndimage_status_success) error stop 'complex uniform selected axes'
    write(unit)zout
    call uniform_filter_nd(z,3,zout,status=status,axes=[integer ::])
    if (status/=ndimage_status_success .or. any(zout/=z)) &
        error stop 'complex uniform empty selected axes'
    write(unit)zout
    close(unit)
    print '(a)', 'PASS: selected axes driver'
end program selected_axes_driver
