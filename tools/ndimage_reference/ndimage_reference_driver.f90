! SPDX-License-Identifier: MIT
! Reproducible binary results for SciPy ndimage comparisons.
program ndimage_reference_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : gaussian_filter_nd, gaussian_filter1d_nd, &
        gaussian_laplace_nd, gaussian_gradient_magnitude_nd, &
        uniform_filter_nd, minimum_filter_nd, maximum_filter_nd, &
        sobel_nd, laplace_nd
    implicit none
    character(len=16), parameter :: modes(8) = [character(len=16) :: &
        'reflect','mirror','nearest','wrap','constant','grid-constant', &
        'grid-mirror','grid-wrap']
    real(dp) :: x(4,5,3),y(4,5,3),v(4,5,3)
    complex(dp) :: z(4,5,3),w(4,5,3)
    integer :: i,j,k,m,s,unit

    do k=1,3
        do j=1,5
            do i=1,4
                x(i,j,k)=sin(0.31_dp*i+0.42_dp*j)+cos(0.23_dp*i-0.19_dp*j)+0.11_dp*k
            end do
        end do
    end do
    z=cmplx(x,0.37_dp*x+0.21_dp,dp)
    open(newunit=unit,file='ndimage_reference_results.bin',access='stream', &
         form='unformatted',status='replace',action='write')
    do m=1,8
        call gaussian_filter_nd(x,[0.6_dp,1.3_dp,0.9_dp],y, &
            status=s,mode=trim(modes(m)),cval=-1.75_dp)
        if (s/=0) error stop 'gaussian real'
        write(unit) y
        call gaussian_filter_nd(z,[0.6_dp,1.3_dp,0.9_dp],w, &
            status=s,mode=trim(modes(m)),cval=cmplx(-1.75_dp,0.6_dp,dp))
        if (s/=0) error stop 'gaussian complex'
        write(unit) w
        call gaussian_filter_nd(x,[1.1_dp,0.8_dp,1.5_dp],y, &
            order=[1,0,2],status=s,mode=trim(modes(m)),cval=-1.75_dp, &
            truncate=3.0_dp,radius=[3,2,4])
        if (s/=0) error stop 'gaussian derivative'
        write(unit) y
        call gaussian_filter1d_nd(x,1.2_dp,2,y,order=3,status=s, &
            mode=trim(modes(m)),cval=-1.75_dp,radius=4)
        if (s/=0) error stop 'gaussian 1d'
        write(unit) y
        call gaussian_gradient_magnitude_nd(x,[1.0_dp,0.8_dp,1.2_dp],y, &
            status=s,mode=trim(modes(m)),cval=-1.75_dp)
        if (s/=0) error stop 'gradient magnitude'
        write(unit) y
        call gaussian_laplace_nd(x,[1.0_dp,0.8_dp,1.2_dp],y, &
            status=s,mode=trim(modes(m)),cval=-1.75_dp)
        if (s/=0) error stop 'gaussian laplace'
        write(unit) y
        call uniform_filter_nd(x,[2,3,4],y,status=s, &
            mode=trim(modes(m)),cval=-1.75_dp)
        if (s/=0) error stop 'uniform real'
        write(unit) y
        call uniform_filter_nd(z,[2,3,4],w,status=s,mode=trim(modes(m)), &
            cval=cmplx(-1.75_dp,0.6_dp,dp))
        if (s/=0) error stop 'uniform complex'
        write(unit) w
        call minimum_filter_nd(x,[2,3,4],y,status=s, &
            mode=trim(modes(m)),cval=-1.75_dp)
        if (s/=0) error stop 'minimum filter'
        write(unit) y
        call maximum_filter_nd(x,[2,3,4],y,status=s, &
            mode=trim(modes(m)),cval=-1.75_dp)
        if (s/=0) error stop 'maximum filter'
        write(unit) y
        call sobel_nd(x,2,y,status=s,mode=trim(modes(m)),cval=-1.75_dp)
        if (s/=0) error stop 'sobel'
        write(unit) y
        call laplace_nd(x,y,status=s,mode=trim(modes(m)),cval=-1.75_dp)
        if (s/=0) error stop 'laplace'
        write(unit) y
        call gaussian_filter_nd(x,1.0_dp,y,status=s,mode=trim(modes(m)))
        if (s/=0) error stop 'gaussian scalar'
        write(unit) y
        call uniform_filter_nd(x,3,y,status=s,mode=trim(modes(m)))
        if (s/=0) error stop 'uniform scalar'
        write(unit) y
        call minimum_filter_nd(x,3,y,status=s,mode=trim(modes(m)))
        if (s/=0) error stop 'minimum scalar'
        write(unit) y
        call maximum_filter_nd(x,3,y,status=s,mode=trim(modes(m)))
        if (s/=0) error stop 'maximum scalar'
        write(unit) y
    end do
    close(unit)
    print '(a)', 'PASS: ndimage reference driver'
end program ndimage_reference_driver
