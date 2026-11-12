! SPDX-License-Identifier: MIT
! Fixed image and mixed-boundary reference stream for SciPy comparisons.
program mixed_complex_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : gaussian_filter_nd, uniform_filter_nd, generic_filter_nd, &
        minimum_filter_nd, maximum_filter_nd, correlate_nd, convolve_nd
    implicit none
    real(dp) :: a(5,6),b(5,6),a3(4,3,2),b3(4,3,2)
    complex(dp) :: z(5,6),w(5,6),cw(2,3)
    real(dp) :: weights(2,3)
    character(len=16), parameter :: all_modes(8) = [character(len=16) :: &
        'reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']
    character(len=16) :: modes(2),triple(3)
    integer :: i,j,k,mode1,mode2,s,unit
    logical :: foot(2,3)
    open(newunit=unit,file='mixed_complex_results.bin',access='stream',form='unformatted',status='replace')
    do j=1,6
        do i=1,5
            a(i,j)=sin(0.4_dp*i+0.7_dp*j)+real(i*i-j,dp)/7.0_dp
            z(i,j)=cmplx(a(i,j),cos(real(i*j,dp)*0.3_dp),dp)
        end do
    end do
    weights=reshape([0.4_dp,-0.6_dp,1.1_dp,0.8_dp,-0.3_dp,0.2_dp],[2,3])
    cw=cmplx(weights,0.3_dp*weights+0.15_dp,dp)
    foot=reshape([.true.,.false.,.true.,.true.,.false.,.true.],shape(foot))
    do mode1=1,8
        do mode2=1,8
            modes=[all_modes(mode1),all_modes(mode2)]
            call gaussian_filter_nd(a,[0.85_dp,1.37_dp],b,order=[1,0], &
                status=s,modes=modes,cval=-0.7_dp,truncate=3.0_dp)
            call check(s,'gaussian real')
            write(unit) b
            call gaussian_filter_nd(z,[0.85_dp,1.37_dp],w,order=[0,2], &
                status=s,modes=modes,cval=cmplx(-0.7_dp,0.6_dp,dp),radius=[2,3])
            call check(s,'gaussian complex')
            write(unit) w
            call uniform_filter_nd(a,[2,3],b,status=s,modes=modes,cval=-0.7_dp)
            call check(s,'uniform real')
            write(unit) b
            call uniform_filter_nd(z,[2,4],w,status=s,modes=modes,cval=cmplx(-0.7_dp,0.6_dp,dp))
            call check(s,'uniform complex')
            write(unit) w
            call generic_filter_nd(a,reduce_real,b,s,footprint=foot,origin=[-1,1], &
                modes=modes,cval=-0.7_dp)
            call check(s,'generic real')
            write(unit) b
            call generic_filter_nd(z,reduce_complex,w,s,footprint=foot,origin=[-1,1], &
                modes=modes,cval=cmplx(-0.7_dp,0.6_dp,dp))
            call check(s,'generic complex')
            write(unit) w
            call minimum_filter_nd(a,[2,3],b,status=s,modes=modes,cval=-0.7_dp)
            call check(s,'minimum mixed')
            write(unit) b
            call maximum_filter_nd(a,[3,2],b,status=s,modes=modes,cval=-0.7_dp)
            call check(s,'maximum mixed')
            write(unit) b
            call correlate_nd(a,weights,b,status=s,modes=modes,cval=-0.7_dp)
            call check(s,'correlate mixed')
            write(unit) b
            call convolve_nd(z,cw,w,status=s,modes=modes,cval=cmplx(-0.7_dp,0.6_dp,dp))
            call check(s,'convolve mixed')
            write(unit) w
        end do
    end do
    do k=1,2
        do j=1,3
            do i=1,4
                a3(i,j,k)=sin(real(i+j*2+k*3,dp)*0.3_dp)
            end do
        end do
    end do
    do mode1=1,8
        triple=[all_modes(mode1),all_modes(mod(mode1,8)+1),all_modes(mod(mode1+3,8)+1)]
        call gaussian_filter_nd(a3,[0.8_dp,0.4_dp,1.2_dp],b3,order=[0,1,0], &
            status=s,modes=triple,cval=1.3_dp)
        call check(s,'gaussian 3d')
        write(unit) b3
        call uniform_filter_nd(a3,[2,3,4],b3,status=s,modes=triple,cval=1.3_dp)
        call check(s,'uniform 3d')
        write(unit) b3
    end do
    close(unit)
    print *, 'PASS: mixed modes and complex generic reference driver'
contains
    subroutine check(status,name)
        integer, intent(in) :: status
        character(len=*), intent(in) :: name
        if (status/=0) then
            print *, 'FAIL: ',name,status
            stop 1
        end if
    end subroutine check

    function reduce_real(values) result(value)
        real(dp), intent(in) :: values(:)
        real(dp) :: value
        integer :: n
        value=0.0_dp
        do n=1,size(values)
            value=value+real(n,dp)*values(n)
        end do
    end function reduce_real

    function reduce_complex(values) result(value)
        complex(dp), intent(in) :: values(:)
        complex(dp) :: value
        integer :: n
        value=(0.0_dp,0.0_dp)
        do n=1,size(values)
            value=value+cmplx(real(n,dp),0.0_dp,dp)*values(n)
        end do
    end function reduce_complex
end program mixed_complex_driver
