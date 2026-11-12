! SPDX-License-Identifier: MIT
! Boundary prefiltering with singleton axes and nonzero cval.
program spline_prepad_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : map_coordinates_nd
    implicit none
    real(dp) :: x(3,1,4), coords(3,31),values(31)
    complex(dp) :: z(3,1,4),results(31)
    integer :: i,k,order,status,unit,m
    character(len=16), parameter :: modes(2)=[character(len=16):: &
        'nearest','grid-constant']
    do k=1,4
        do i=1,3
            x(i,1,k)=sin(0.3_dp*i)+cos(0.21_dp*k)+0.15_dp*i*k
        end do
    end do
    z=cmplx(x,0.3_dp*x+0.2_dp,dp)
    do k=1,31
        coords(:,k)=[1.0_dp+0.5_dp*real(modulo(k*7,17)-4,dp), &
                     1.0_dp+0.4_dp*real(modulo(k*3,9)-4,dp), &
                     1.0_dp+0.6_dp*real(modulo(k*5,19)-6,dp)]
    end do
    coords(:,1)=[1.0_dp,1.0_dp,1.0_dp]
    coords(:,2)=[3.0_dp,1.0_dp,4.0_dp]
    coords(:,3)=[3.0_dp,0.5_dp,4.0_dp]
    coords(:,4)=[1.0_dp,1.5_dp,1.0_dp]
    open(newunit=unit,file='spline_prepad_results.bin',access='stream', &
         form='unformatted',status='replace',action='write')
    do m=1,2
        do order=2,5
            call map_coordinates_nd(x,coords,values,order,status,mode=trim(modes(m)), &
                prefilter=.true.,cval=-2.0_dp)
            if (status/=0) error stop 'real prepad'
            write(unit) values
            call map_coordinates_nd(z,coords,results,order,status,mode=trim(modes(m)), &
                prefilter=.true.,cval=cmplx(-2.0_dp,0.75_dp,dp))
            if (status/=0) error stop 'complex prepad'
            write(unit) results
        end do
    end do
    close(unit)
end program spline_prepad_reference_driver
