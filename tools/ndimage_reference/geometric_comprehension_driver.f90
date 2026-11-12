! SPDX-License-Identifier: MIT
! Generate deterministic geometric-interpolation and region-callback cases.
program geometric_comprehension_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : geometric_transform_nd, labeled_comprehension_nd
    implicit none
    real(dp) :: a2(6,7),r2(5,8),a3(4,3,5),r3(6,4),a1(6),r1(3,2,2)
    complex(dp) :: z2(6,7),c2(5,8)
    integer :: labs(6,7),requests(8),i,j,k,mi,oi,s,unit
    integer, parameter :: orders(4)=[0,1,3,5]
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'reflect','grid-mirror','mirror','nearest', &
        'wrap','grid-wrap','constant','grid-constant']
    real(dp) :: reductions(8)

    do j=1,7
        do i=1,6
            a2(i,j)=sin(0.17_dp*i+0.31_dp*j)+cos(0.23_dp*i-0.51_dp*j)
            z2(i,j)=cmplx(a2(i,j),0.3_dp*a2(i,j)+sin(0.21_dp*j),dp)
            labs(i,j)=modulo(i+3*j,5)-2
        end do
    end do
    do k=1,5
        do j=1,3
            do i=1,4
                a3(i,j,k)=sin(0.2_dp*i+0.7_dp*j-0.3_dp*k)+0.03_dp*i*k
            end do
        end do
    end do
    do i=1,6
        a1(i)=sin(0.8_dp*i)+0.03_dp*i*i
    end do
    requests=[-2,-1,0,1,2,99,1,-2]
    open(newunit=unit,file='geometric_comprehension.bin',access='stream', &
        form='unformatted',status='replace')
    do mi=1,8
        do oi=1,4
            call geometric_transform_nd(a2,map2,r2,order=orders(oi),status=s, &
                mode=trim(modes(mi)),cval=-0.75_dp)
            if (s/=0) error stop 'geometric real'
            write(unit) r2
            call geometric_transform_nd(z2,map2,c2,order=orders(oi),status=s, &
                mode=trim(modes(mi)),cval=cmplx(-0.75_dp,0.4_dp,dp))
            if (s/=0) error stop 'geometric complex'
            write(unit) c2
        end do
    end do
    call geometric_transform_nd(a2,map2,r2,order=3,status=s,mode='mirror', &
        prefilter=.false.,cval=-0.75_dp)
    if (s/=0) error stop 'no prefilter'
    write(unit) r2
    call geometric_transform_nd(a3,map3,r3,order=3,status=s,mode='nearest')
    if (s/=0) error stop 'rank reducing map'
    write(unit) r3
    call geometric_transform_nd(a1,map1,r1,order=1,status=s,mode='reflect')
    if (s/=0) error stop 'rank increasing map'
    write(unit) r1
    call labeled_comprehension_nd(a2,labs,requests,reduce_symmetric,-19.5_dp, &
        reductions,s)
    if (s/=0) error stop 'batched comprehension'
    write(unit) reductions
    close(unit)
    print '(a)', 'PASS: geometric transform and labeled comprehension driver'
contains
    subroutine map2(out,src)
        integer, intent(in) :: out(:) !! one-based two-dimensional destination
        real(dp), intent(out) :: src(:) !! one-based two-dimensional source
        src(1)=real(out(1),dp)+0.35_dp*real(out(2)-3,dp)
        src(2)=real(out(2),dp)-0.2_dp*real(out(1)-2,dp)
    end subroutine map2

    subroutine map3(out,src)
        integer, intent(in) :: out(:) !! one-based two-dimensional destination
        real(dp), intent(out) :: src(:) !! one-based three-dimensional source
        src(1)=1.0_dp+0.4_dp*real(out(1),dp)
        src(2)=1.5_dp+0.25_dp*real(out(2),dp)
        src(3)=2.0_dp+sin(0.11_dp*real(out(1),dp))
    end subroutine map3

    subroutine map1(out,src)
        integer, intent(in) :: out(:) !! one-based three-dimensional destination
        real(dp), intent(out) :: src(:) !! one-based one-dimensional source
        src(1)=1.0_dp+0.8_dp*real(out(1)-1,dp)+0.33_dp*real(out(2)-1,dp) &
            -0.1_dp*real(out(3)-1,dp)
    end subroutine map1

    function reduce_symmetric(samples) result(value)
        real(dp), intent(in) :: samples(:) !! all pixels of selected label
        real(dp) :: value !! nonlinear symmetric sum
        value=sum(sin(samples)+samples*samples)
    end function reduce_symmetric
end program geometric_comprehension_driver
