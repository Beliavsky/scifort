! SPDX-License-Identifier: MIT
! Finite full-dimensional kernels and real neighborhood order-statistics.
program extended_filter_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : correlate_nd, convolve_nd, median_filter_nd, &
        prewitt_nd, &
        rank_filter_nd, percentile_filter_nd
    implicit none
    real(dp) :: x(4,5),y(4,5),w(2,3),v(4,2),a(4,3,2),b(4,3,2),w3(2,3,2)
    complex(dp) :: z(4,5),zz(4,5),cw(2,3)
    character(len=15), parameter :: modes(8)=[character(len=15) :: &
        'reflect','mirror','nearest','wrap','constant','grid-constant', &
        'grid-mirror','grid-wrap']
    integer :: i,j,k,m,s,u

    do j=1,5
        do i=1,4
            x(i,j)=sin(0.3_dp*i+0.7_dp*j)+real(modulo(i*17+j*13,9),dp)
            z(i,j)=cmplx(x(i,j),-0.37_dp*x(i,j)+0.11_dp*i,dp)
        end do
    end do
    do j=1,3
        do i=1,2
            w(i,j)=sin(0.6_dp*i-0.29_dp*j)+0.1_dp*i*j
            cw(i,j)=cmplx(w(i,j),0.16_dp*j-0.23_dp*i,dp)
        end do
    end do
    do j=1,2
        do i=1,4
            v(i,j)=cos(0.21_dp*i+0.44_dp*j)+0.13_dp*i*j
        end do
    end do
    do k=1,2
        do j=1,3
            do i=1,4
                a(i,j,k)=sin(0.52_dp*i-0.19_dp*j+0.11_dp*k)+ &
                    real(modulo(i*13+j*7+k*11,7),dp)
            end do
        end do
    end do
    do k=1,2
        do j=1,3
            do i=1,2
                w3(i,j,k)=sin(0.22_dp*i+0.17_dp*j-0.28_dp*k)
            end do
        end do
    end do
    open(newunit=u,file='extended_filter_results.bin',status='replace', &
        form='unformatted',access='stream',action='write')
    do m=1,size(modes)
        call correlate_nd(x,w,y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '2D correlate'
        write(u) y
        call convolve_nd(x,w,y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '2D convolve'
        write(u) y
        call correlate_nd(x,v,y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '2D correlate even'
        write(u) y
        call convolve_nd(x,v,y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '2D convolve even'
        write(u) y
        call median_filter_nd(x,[2,3],y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop 'median vector'
        write(u) y
        call rank_filter_nd(x,0,[3,2],y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop 'rank 0'
        write(u) y
        call rank_filter_nd(x,-2,[2,3],y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop 'rank negative'
        write(u) y
        call percentile_filter_nd(x,25.0_dp,[2,3],y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop 'percentile 25'
        write(u) y
        call percentile_filter_nd(x,75.0_dp,[3,2],y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop 'percentile 75'
        write(u) y
        call median_filter_nd(x,3,y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop 'median scalar'
        write(u) y
        call prewitt_nd(x,2,y,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '2D prewitt'
        write(u) y
        call correlate_nd(z,cw,zz,s,modes(m),cval=cmplx(-2.3_dp,0.7_dp,dp))
        if (s/=0) error stop 'complex correlate'
        write(u) zz
        call convolve_nd(z,cw,zz,s,modes(m),cval=cmplx(-2.3_dp,0.7_dp,dp))
        if (s/=0) error stop 'complex convolve'
        write(u) zz
        call correlate_nd(a,w3,b,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '3D correlate'
        write(u) b
        call convolve_nd(a,w3,b,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '3D convolve'
        write(u) b
        call rank_filter_nd(a,4,[3,2,2],b,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '3D rank'
        write(u) b
        call percentile_filter_nd(a,95.0_dp,[2,3,2],b,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '3D percentile'
        write(u) b
        call median_filter_nd(a,[3,1,2],b,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '3D median'
        write(u) b
        call prewitt_nd(a,3,b,s,modes(m),cval=-2.3_dp)
        if (s/=0) error stop '3D prewitt'
        write(u) b
    end do
    close(u)
end program extended_filter_driver
