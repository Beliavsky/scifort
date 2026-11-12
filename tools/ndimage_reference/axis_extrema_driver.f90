! SPDX-License-Identifier: MIT
! Binary reference arrays for 1D uniform/minimum/maximum filters.
program axis_extrema_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : uniform_filter1d_nd, minimum_filter1d_nd, maximum_filter1d_nd
    implicit none
    real(dp) :: a(5,7),b(5,7)
    complex(dp) :: z(5,7),w(5,7)
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'reflect','grid-mirror','mirror','nearest','wrap', &
        'grid-wrap','constant','grid-constant']
    integer, parameter :: widths(5)=[1,2,3,4,6]
    integer :: i,j,d,m,k,c,width,origin,s,u
    do j=1,7
        do i=1,5
            a(i,j)=sin(0.4_dp*(i*2+j))+real(i*j,dp)/17.0_dp
            z(i,j)=cmplx(a(i,j),cos(.6_dp*i+.17_dp*j),dp)
        end do
    end do
    open(newunit=u,file='axis_extrema_results.bin',form='unformatted', &
        access='stream',status='replace',action='write')
    do m=1,size(modes)
        do k=1,size(widths)
            width=widths(k)
            do c=1,3
                select case(c)
                case(1)
                    origin=-(width/2)
                case(2)
                    origin=0
                case default
                    origin=(width-1)/2
                end select
                do d=1,2
                    call uniform_filter1d_nd(a,width,d,b,status=s,mode=trim(modes(m)), &
                        cval=-1.75_dp,origin=origin)
                    if (s/=0) error stop 'uniform real 1d'
                    write(u) b
                    call minimum_filter1d_nd(a,width,d,b,status=s,mode=trim(modes(m)), &
                        cval=-1.75_dp,origin=origin)
                    if (s/=0) error stop 'minimum 1d'
                    write(u) b
                    call maximum_filter1d_nd(a,width,d,b,status=s,mode=trim(modes(m)), &
                        cval=-1.75_dp,origin=origin)
                    if (s/=0) error stop 'maximum 1d'
                    write(u) b
                    call uniform_filter1d_nd(z,width,d,w,status=s,mode=trim(modes(m)), &
                        cval=cmplx(-1.75_dp,0.6_dp,dp),origin=origin)
                    if (s/=0) error stop 'uniform complex 1d'
                    write(u) w
                end do
            end do
        end do
    end do
    close(u)
end program axis_extrema_driver
