! SPDX-License-Identifier: MIT
! Binary reference output for periodic coefficients and shaped coordinates.
program spline_periodic_grid_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_filter_nd, &
        map_coordinates_grid_nd
    implicit none
    integer, parameter :: lengths(11)=[2,3,4,5,7,11,17,63,128,512,1024]
    character(len=13), parameter :: modes(5)=[character(len=13) :: &
        'grid-wrap','nearest','grid-constant','constant','mirror']
    real(dp) :: line(1024), filtered(1024), sample(7,9), coord(2,5,4), &
        output(5,4)
    complex(dp) :: z(7,9), z_out(5,4), zline(1024), zcoeff(1024)
    integer :: i, j, k, n, degree, status

    open(newunit=k,file='spline_periodic_grid.bin',form='unformatted', &
        access='stream',status='replace')
    do j=1,size(lengths)
        n=lengths(j)
        do i=1,n
            line(i)=sin(0.37_dp*real(i,dp))+.013_dp*i*i
        end do
        zline(:n)=cmplx(line(:n),-.11_dp*line(:n)+.07_dp,dp)
        do degree=2,5
            call spline_filter_nd(line(:n),filtered(:n),order=degree, &
                status=status,mode='grid-wrap')
            if (status/=0) error stop 'periodic filter failed'
            write(k) filtered(:n)
            call spline_filter_nd(zline(:n),zcoeff(:n),order=degree, &
                status=status,mode='grid-wrap')
            if (status/=0) error stop 'complex periodic filter failed'
            write(k) zcoeff(:n)
        end do
    end do
    do j=1,9
        do i=1,7
            sample(i,j)=sin(.21_dp*i+.33_dp*j)+cos(.07_dp*i*j)
        end do
    end do
    z=cmplx(sample,0.7_dp*sample-0.13_dp,dp)
    do j=1,4
        do i=1,5
            coord(:,i,j)=[-.2_dp+.91_dp*i+.43_dp*j, 1.3_dp+.55_dp*i+1.2_dp*j]
        end do
    end do
    do j=1,size(modes)
        do degree=0,5
            call map_coordinates_grid_nd(sample,coord,output,order=degree, &
                status=status,mode=trim(modes(j)),cval=-1.25_dp)
            if (status/=0) error stop 'real grid failed'
            write(k) output
            call map_coordinates_grid_nd(z,coord,z_out,order=degree, &
                status=status,mode=trim(modes(j)),cval=(-1.25_dp,0.42_dp))
            if (status/=0) error stop 'complex grid failed'
            write(k) z_out
        end do
    end do
    ! Explicit no-prefilter path, preserving source coefficients.
    call map_coordinates_grid_nd(sample,coord,output,order=5,mode='grid-wrap', &
        prefilter=.false.,status=status)
    if (status/=0) error stop 'raw grid failed'
    write(k) output
    close(k)
end program spline_periodic_grid_driver
