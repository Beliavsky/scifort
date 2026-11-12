! SPDX-License-Identifier: MIT
! Numeric reference driver for multidimensional spline boundary modes and transforms.
program spline_modes_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : map_coordinates_nd, affine_transform_nd, shift_nd
    implicit none
    character(len=16), parameter :: modes(8) = [character(len=16) :: &
        'mirror','constant','reflect','grid-mirror', 'nearest','wrap','grid-wrap','grid-constant']
    real(dp) :: x(4,5),coords(2,41),v(41),a(3,6),sh(4,5),matrix(2,2),offset(2)
    complex(dp) :: z(4,5),vc(41),ac(3,6),sc(4,5)
    integer :: i,j,k,m,o,p,status,s2,s3,s4,s5,s6
    do j=1,5
        do i=1,4
            x(i,j)=sin(0.31_dp*real(i,dp)+0.42_dp*real(j,dp)) + &
                cos(0.23_dp*real(i,dp)-0.19_dp*real(j,dp))
        end do
    end do
    z=cmplx(x,0.37_dp*x+0.21_dp,dp)
    do k=1,41
        coords(1,k)=0.2_dp+0.37_dp*real(k-11,dp)
        coords(2,k)=1.3_dp+0.49_dp*real(modulo(9*k,41)-12,dp)
    end do
    ! Anchor exact boundaries and fractional coordinates.
    coords(:,1)=[1.0_dp,1.0_dp]
    coords(:,2)=[4.0_dp,5.0_dp]
    coords(:,3)=[0.35_dp,5.4_dp]
    coords(:,4)=[-2.0_dp,6.0_dp]
    coords(:,5)=[8.0_dp,-4.0_dp]
    coords(:,6)=[0.5_dp,5.5_dp]
    coords(:,7)=[0.5_dp,1.0_dp]
    coords(:,8)=[4.5_dp,1.0_dp]
    coords(:,9)=[1.0_dp,5.5_dp]
    coords(:,10)=[4.0_dp,0.5_dp]
    coords(:,11)=[3.5_dp,5.5_dp]
    matrix(1,:)=[0.83_dp,0.12_dp]
    matrix(2,:)=[-0.16_dp,0.93_dp]
    offset=[0.3_dp,-0.5_dp]
    open(newunit=i,file='spline_modes_results.bin',access='stream',status='replace', &
        form='unformatted',action='write')
    do m=1,size(modes)
        do o=0,5
            do p=0,1
                call map_coordinates_nd(x,coords,v,order=o,mode=trim(modes(m)), &
                    prefilter=p==1,cval=-2.25_dp,status=status)
                call map_coordinates_nd(z,coords,vc,order=o,mode=trim(modes(m)), &
                    prefilter=p==1,cval=cmplx(-2.25_dp,0.5_dp,dp),status=s2)
                call affine_transform_nd(x,matrix,offset,a,order=o,mode=trim(modes(m)), &
                    prefilter=p==1,cval=-2.25_dp,status=s3)
                call affine_transform_nd(z,matrix,offset,ac,order=o,mode=trim(modes(m)), &
                    prefilter=p==1,cval=cmplx(-2.25_dp,0.5_dp,dp),status=s4)
                call shift_nd(x,[0.45_dp,-1.3_dp],sh,order=o,mode=trim(modes(m)), &
                    prefilter=p==1,cval=-2.25_dp,status=s5)
                call shift_nd(z,[0.45_dp,-1.3_dp],sc,order=o,mode=trim(modes(m)), &
                    prefilter=p==1,cval=cmplx(-2.25_dp,0.5_dp,dp),status=s6)
                write(i) status,s2,s3,s4,s5,s6
                if (status==0) write(i) v,vc,a,ac,sh,sc
            end do
        end do
    end do
    close(i)
end program spline_modes_reference_driver
