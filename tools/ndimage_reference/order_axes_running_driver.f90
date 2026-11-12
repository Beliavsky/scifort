! SPDX-License-Identifier: MIT
! Independent numeric outputs for SciPy order-filter axes and running means.
program order_axes_running_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : median_filter_nd, rank_filter_nd, percentile_filter_nd, &
        uniform_filter1d_nd, uniform_filter_nd
    implicit none
    real(dp) :: x(4,5,3),y(4,5,3)
    complex(dp) :: z(4,5,3),w(4,5,3)
    logical :: fp(2,3)
    integer :: i,j,k,u,status,m,n
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'reflect','mirror','nearest','wrap','grid-wrap','constant','grid-constant','grid-mirror']

    do k=1,3
        do j=1,5
            do i=1,4
                x(i,j,k)=sin(.17_dp*i+.31_dp*j-.19_dp*k)+.05_dp*(i*i-2*j+k)
            end do
        end do
    end do
    z=cmplx(x,.4_dp*x+.2_dp,dp)
    fp=reshape([.true.,.false.,.true.,.true.,.false.,.true.],shape(fp))
    open(newunit=u,file='order_axes_output.bin',status='replace',access='stream',form='unformatted')
    call median_filter_nd(x,[4,3],y,status=status,mode='nearest',axes=[3,1])
    call check_status('median axes 3,1')
    write(u) y
    call median_filter_nd(x,5,y,status=status,mode='mirror',origin=[-1],axes=[1])
    call check_status('median origin')
    write(u) y
    call median_filter_nd(x,3,y,status=status,axes=[integer ::])
    call check_status('median no axes')
    write(u) y
    call rank_filter_nd(x,-2,[2,4],y,status=status,mode='reflect',origin=[0,-1],axes=[1,3])
    call check_status('rank axes/origins')
    write(u) y
    call rank_filter_nd(x,0,y=y,status=status,footprint=fp,origin=[0,0],axes=[3,1], &
        mode='constant',cval=-1.5_dp)
    call check_status('rank footprint')
    write(u) y
    call percentile_filter_nd(x,35.0_dp,[4,3],y,status=status,origin=[1,-1],axes=[2,1],mode='wrap')
    call check_status('percentile axes')
    write(u) y
    call percentile_filter_nd(x,100.0_dp,y=y,status=status,footprint=fp,origin=[0,0], &
        axes=[3,1],mode='mirror')
    call check_status('percentile footprint')
    write(u) y
    call percentile_filter_nd(x,80.0_dp,3,y,status=status,axes=[integer ::])
    call check_status('percentile no axes')
    write(u) y
    call median_filter_nd(x,[2,3,4],y,status=status,mode='constant',cval=-3.5_dp)
    call check_status('median all axes')
    write(u) y
    call rank_filter_nd(x,-1,3,y,status=status,mode='nearest')
    call check_status('rank all axes')
    write(u) y
    call percentile_filter_nd(x,75.0_dp,2,y,status=status,mode='mirror')
    call check_status('percentile all axes')
    write(u) y
    call rank_filter_nd(x,0,3,y,status=status,axes=[2],mode='reflect')
    call check_status('rank single axis')
    write(u) y
    call rank_filter_nd(x,1,y=y,status=status,footprint=fp,origin=[-1,0], &
        axes=[3,1],mode='nearest')
    call check_status('unsorted axes shifted footprint')
    write(u) y
    call median_filter_nd(x,3,y,status=status,footprint=fp,origin=[-1,1], &
        axes=[1,3],mode='reflect')
    call check_status('sorted axes shifted footprint')
    write(u) y
    call percentile_filter_nd(x,60.0_dp,[2,4],y,status=status,axes=[3,1], &
        origin=[-1,0],mode='nearest')
    call check_status('unsorted origin interpretation')
    write(u) y
    close(u)

    open(newunit=u,file='running_uniform_output.bin' ,status='replace',access='stream',form='unformatted')
    call uniform_filter1d_nd(x,6,1,y,status=status,mode='reflect',origin=-1)
    call check_status('real uniform1d reflect')
    write(u) y
    call uniform_filter1d_nd(x,4,3,y,status=status,mode='grid-wrap',origin=1)
    call check_status('real uniform1d wrap')
    write(u) y
    call uniform_filter1d_nd(x,50000,2,y,status=status,mode='constant',cval=3.25_dp)
    call check_status('wide constant uniform')
    write(u) y
    call uniform_filter1d_nd(x,50000,2,y,status=status,mode='nearest',origin=-25000)
    call check_status('wide nearest uniform')
    write(u) y
    call uniform_filter1d_nd(x,3,2,y,status=status,mode='mirror')
    call check_status('real uniform1d mirror')
    write(u) y
    call uniform_filter_nd(x,[4,3],y,status=status,axes=[3,1], &
        modes=[character(len=16) :: 'wrap','constant'],cval=-1.25_dp)
    call check_status('real selected axes uniform')
    write(u) y
    call uniform_filter_nd(x,3,y,status=status,mode='wrap')
    call check_status('real full uniform')
    write(u) y
    call uniform_filter_nd(x,3,y,status=status,axes=[integer ::])
    call check_status('real empty axes')
    write(u) y
    call uniform_filter1d_nd(x,1,2,y,status=status)
    call check_status('real width one')
    write(u) y
    close(u)

    open(newunit=u,file='running_uniform_complex_output.bin',status='replace',access='stream',form='unformatted')
    call uniform_filter1d_nd(z,6,1,w,status=status,mode='constant',cval=cmplx(1.25_dp,-.7_dp,dp),origin=-1)
    call check_status('complex uniform1d constant')
    write(u) w
    call uniform_filter1d_nd(z,5,3,w,status=status,mode='reflect',origin=1)
    call check_status('complex uniform1d reflect')
    write(u) w
    call uniform_filter_nd(z,[3,4],w,status=status,axes=[2,1],mode='mirror')
    call check_status('complex selected uniform')
    write(u) w
    close(u)
    open(newunit=u,file='order_axes_random_output.bin',status='replace',access='stream',form='unformatted')
    do m=1,size(modes)
        call median_filter_nd(x,[2,4],y,status=status,mode=trim(modes(m)),axes=[3,1])
        call check_status('boundary median')
        write(u) y
        call rank_filter_nd(x,2,[2,4],y,status=status,mode=trim(modes(m)),axes=[3,1])
        call check_status('boundary rank')
        write(u) y
        call percentile_filter_nd(x,55.0_dp,[3,2],y,status=status,mode=trim(modes(m)),axes=[1,3])
        call check_status('boundary percentile')
        write(u) y
    end do
    close(u)
    open(newunit=u,file='running_uniform_random_output.bin',status='replace',access='stream',form='unformatted')
    do m=1,size(modes)
        do n=2,9
            call uniform_filter1d_nd(x,n,2,y,status=status,mode=trim(modes(m)))
            call check_status('random uniform1d')
            write(u) y
        end do
    end do
    close(u)
    open(newunit=u,file='running_uniform_random_complex.bin',status='replace',access='stream',form='unformatted')
    do m=1,size(modes)
        do n=2,9
            call uniform_filter1d_nd(z,n,3,w,status=status,mode=trim(modes(m)))
            call check_status('random complex uniform1d')
            write(u) w
        end do
    end do
    close(u)
    print '(a)', 'PASS: order-statistic axes and uniform running-sum reference driver' 
contains
    subroutine check_status(label)
        character(len=*), intent(in) :: label
        if (status/=0) then
            print '(a,1x,a,1x,i0)', 'FAIL:',label,status
            error stop 1
        end if
    end subroutine check_status
end program order_axes_running_driver
