! SPDX-License-Identifier: MIT
! Reference stream for Fourier ellipsoids, mixed order filters and line callbacks.
program ellipsoid_order_line_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : fourier_ellipsoid_nd, median_filter_nd, &
        rank_filter_nd, percentile_filter_nd, generic_filter1d_nd
    implicit none
    complex(dp) :: f1(7),g1(7),f2(5,6),g2(5,6),f3(4,3,2),g3(4,3,2)
    complex(dp) :: r2(5,4),u2(5,4),r3(2,3,5),u3(2,3,5),a4(3,2,2,2),b4(3,2,2,2)
    complex(dp) :: a15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: b15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: a(5,6),b(5,6),c(3,4,2),d(3,4,2)
    logical :: foot(2,3),foot3(2,3,2)
    character(len=16), parameter :: all_modes(8)=[character(len=16) :: &
        'reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']
    character(len=16) :: modes(2),triple(3)
    integer :: unit,s,i,j,k,l,p,q

    open(newunit=unit,file='ellipsoid_order_line_results.bin', &
        access='stream',form='unformatted',status='replace')
    do i=1,7
        f1(i)=cmplx(sin(0.31_dp*i),cos(0.47_dp*i),dp)
    end do
    call fourier_ellipsoid_nd(f1,2.3_dp,g1,s)
    call check(s,'ellipsoid 1d')
    write(unit) g1
    call fourier_ellipsoid_nd(f1,-3.1_dp,g1,s)
    call check(s,'ellipsoid 1d signed size')
    write(unit) g1
    do j=1,6
        do i=1,5
            f2(i,j)=cmplx(sin(0.4_dp*i+0.23_dp*j),cos(0.27_dp*i-0.49_dp*j),dp)
        end do
    end do
    call fourier_ellipsoid_nd(f2,[1.2_dp,3.1_dp],g2,s)
    call check(s,'ellipsoid 2d vector')
    write(unit) g2
    call fourier_ellipsoid_nd(f2,2.7_dp,g2,s)
    call check(s,'ellipsoid 2d scalar')
    write(unit) g2
    do k=1,2
        do j=1,3
            do i=1,4
                f3(i,j,k)=cmplx(sin(.4_dp*i+.7_dp*j-.2_dp*k),sin(.31_dp*i*j*k),dp)
            end do
        end do
    end do
    call fourier_ellipsoid_nd(f3,[1.8_dp,.2_dp,3.4_dp],g3,s)
    call check(s,'ellipsoid 3d vector')
    write(unit) g3
    call fourier_ellipsoid_nd(f3,1.3_dp,g3,s)
    call check(s,'ellipsoid 3d scalar')
    write(unit) g3
    do j=1,4
        do i=1,5
            r2(i,j)=cmplx(sin(.3_dp*i+j),cos(.44_dp*i-j),dp)
        end do
    end do
    call fourier_ellipsoid_nd(r2,[.9_dp,2.4_dp],u2,s,n=8,axis=1)
    call check(s,'ellipsoid rfft first')
    write(unit) u2
    do k=1,5
        do j=1,3
            do i=1,2
                r3(i,j,k)=cmplx(sin(.5_dp*i+.6_dp*j+.3_dp*k),cos(.41_dp*i-.53_dp*j*k),dp)
            end do
        end do
    end do
    call fourier_ellipsoid_nd(r3,[1.1_dp,1.6_dp,3.8_dp],u3,s,n=9,axis=-1)
    call check(s,'ellipsoid rfft last')
    write(unit) u3
    do l=1,2
        do k=1,2
            do j=1,2
                do i=1,3
                    a4(i,j,k,l)=cmplx(sin(.17_dp*i*j+.11_dp*k*l),cos(.21_dp*i*k-.13_dp*j*l),dp)
                end do
            end do
        end do
    end do
    call fourier_ellipsoid_nd(a4,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],b4,s)
    call check(s,'ellipsoid 4d extension')
    write(unit) b4
    a15=cmplx(2.0_dp,-3.0_dp,dp)
    call fourier_ellipsoid_nd(a15,1.5_dp,b15,s)
    call check(s,'ellipsoid rank15 extension')
    write(unit) b15
    do p=1,4
        select case(p)
        case(1)
            call fourier_ellipsoid_nd(a15,5.2_dp,b15,s)
        case(2)
            call fourier_ellipsoid_nd(a15,12.0_dp,b15,s)
        case(3)
            call fourier_ellipsoid_nd(a15,43.0_dp,b15,s)
        case(4)
            call fourier_ellipsoid_nd(a15,100.0_dp,b15,s)
        end select
        call check(s,'ellipsoid rank15 large frequency')
        write(unit) b15
    end do

    do j=1,6
        do i=1,5
            a(i,j)=sin(.4_dp*i+.7_dp*j)+real(i*i-j,dp)/7.0_dp
        end do
    end do
    foot=reshape([.true.,.false.,.true.,.true.,.false.,.true.],shape(foot))
    do p=1,8
        do q=1,8
            modes=[all_modes(p),all_modes(q)]
            call median_filter_nd(a,[2,3],b,status=s,modes=modes,cval=-.7_dp,origin=[-1,1])
            call check(s,'median mixed')
            write(unit) b
            call median_filter_nd(a,[2,3],b,status=s,modes=modes,cval=-.7_dp, &
                footprint=foot,origin=[-1,1])
            call check(s,'median mixed footprint')
            write(unit) b
            call rank_filter_nd(a,-2,[2,3],b,status=s,modes=modes,cval=-.7_dp, &
                footprint=foot,origin=[-1,1])
            call check(s,'rank mixed footprint')
            write(unit) b
            call percentile_filter_nd(a,35.0_dp,[2,3],b,status=s,modes=modes,cval=-.7_dp, &
                footprint=foot,origin=[-1,1])
            call check(s,'percentile mixed footprint')
            write(unit) b
        end do
    end do
    foot3=reshape([.true.,.true.,.false.,.true.,.true.,.false., &
                   .false.,.true.,.true.,.false.,.true.,.true.],shape(foot3))
    do k=1,2
        do j=1,4
            do i=1,3
                c(i,j,k)=sin(.23_dp*i+.51_dp*j+.37_dp*k)+i*k*.09_dp
            end do
        end do
    end do
    do p=1,8
        triple=[all_modes(p),all_modes(mod(p,8)+1),all_modes(mod(p+3,8)+1)]
        call median_filter_nd(c,[2,3,2],d,status=s,modes=triple,cval=1.3_dp, &
            footprint=foot3,origin=[-1,0,-1])
        call check(s,'3d median footprint')
        write(unit) d
        call rank_filter_nd(c,2,[2,3,2],d,status=s,modes=triple,cval=1.3_dp, &
            footprint=foot3,origin=[-1,0,-1])
        call check(s,'3d rank footprint')
        write(unit) d
        call percentile_filter_nd(c,72.0_dp,[2,3,2],d,status=s,modes=triple,cval=1.3_dp, &
            footprint=foot3,origin=[-1,0,-1])
        call check(s,'3d percentile footprint')
        write(unit) d
    end do
    do p=1,8
        call generic_filter1d_nd(a,line_reduce,4,1,b,status=s, &
            mode=all_modes(p),cval=-.7_dp,origin=1)
        call check(s,'line filter axis1')
        write(unit) b
        call generic_filter1d_nd(a,line_reduce,3,2,b,status=s, &
            mode=all_modes(p),cval=-.7_dp,origin=-1)
        call check(s,'line filter axis2')
        write(unit) b
        call generic_filter1d_nd(c,line_reduce,2,3,d,status=s, &
            mode=all_modes(p),cval=1.3_dp,origin=0)
        call check(s,'line filter axis3')
        write(unit) d
    end do
    close(unit)
    print '(a)', 'PASS: Fourier ellipsoids, mixed order filters and line callbacks'
contains
    subroutine check(status,name)
        integer, intent(in) :: status
        character(len=*), intent(in) :: name
        if (status/=0) then
            print '(3a,i0)', 'FAIL ',name,': status=',status
            error stop 1
        end if
    end subroutine check
    subroutine line_reduce(input_line,output_line)
        real(dp), intent(in) :: input_line(:)
        real(dp), intent(out) :: output_line(:)
        integer :: j,n,filter_width
        n=size(output_line)
        filter_width=size(input_line)-n+1
        do j=1,n
            output_line(j)=input_line(j)+2.0_dp*input_line(j+1)-0.5_dp*input_line(j+filter_width-1)
        end do
    end subroutine line_reduce
end program ellipsoid_order_line_driver
