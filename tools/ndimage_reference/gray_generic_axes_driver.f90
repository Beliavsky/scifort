! SPDX-License-Identifier: MIT
! Numerical output for SciPy-selected-axis grayscale morphology and generic filters.
program gray_generic_axes_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : grey_dilation_nd, grey_erosion_nd, grey_opening_nd, &
        grey_closing_nd, morphological_gradient_nd, morphological_laplace_nd, &
        white_tophat_nd, black_tophat_nd, generic_filter_nd
    implicit none
    real(dp) :: x(4,5,3),y(4,5,3),s(2,3)
    complex(dp) :: z(4,5,3),cy(4,5,3)
    logical :: fp(2,3)
    integer :: i,j,k,u,code,m
    character(len=16), parameter :: mode_names(8)=[character(len=16) :: &
        'reflect','nearest','mirror','wrap','constant','grid-constant', &
        'grid-wrap','grid-mirror']

    do k=1,3
        do j=1,5
            do i=1,4
                x(i,j,k)=sin(.17_dp*i+.31_dp*j-.19_dp*k)+ &
                    .05_dp*(i*i-2*j+k)
            end do
        end do
    end do
    z=cmplx(x,.35_dp*x+.11_dp,dp)
    fp=reshape([.true.,.false.,.true.,.true.,.false.,.true.],[2,3])
    s=reshape([-.3_dp,.8_dp,0.2_dp,-.1_dp,.4_dp,-.2_dp],[2,3])
    open(newunit=u,file='gray_axes_output.bin',access='stream',form='unformatted',status='replace')
    call grey_dilation_nd(x,y,code,sizes=[2,3],axes=[3,1]); call checked_real()
    call grey_erosion_nd(x,y,code,sizes=[3,2],axes=[1,3],origin=[-1,0], &
        mode='constant',cval=-3.0_dp); call checked_real()
    call grey_opening_nd(x,y,code,sizes=[2,3],axes=[3,1],origin=[0,1]); call checked_real()
    call grey_closing_nd(x,y,code,sizes=[2,3],axes=[3,1],mode='nearest'); call checked_real()
    call morphological_gradient_nd(x,y,code,footprint=fp,axes=[3,1],mode='mirror'); call checked_real()
    call morphological_laplace_nd(x,y,code,structure=s,axes=[3,1],mode='grid-wrap'); call checked_real()
    call white_tophat_nd(x,y,code,footprint=fp,structure=s,axes=[3,1],mode='nearest'); call checked_real()
    call black_tophat_nd(x,y,code,footprint=fp,structure=s,axes=[3,1],mode='nearest'); call checked_real()
    call grey_dilation_nd(x,y,code,sizes=[2,3],axes=[3,1],mode='constant',cval=-1.2_dp); call checked_real()
    call grey_erosion_nd(x,y,code,sizes=[integer::],axes=[integer::]); call checked_real()
    call grey_dilation_nd(x,y,code,sizes=[2,3,2],mode='reflect'); call checked_real()
    call grey_dilation_nd(x,y,code,sizes=[2,3],axes=[3,1],modes=[character(len=16):: &
        'constant','nearest'],cval=-1.2_dp); call checked_real()
    call grey_dilation_nd(x,y,code,footprint=fp,structure=s,axes=[3,1], &
        modes=[character(len=16):: 'grid-wrap','constant'],cval=-2.0_dp)
    call checked_real()
    call morphological_gradient_nd(x,y,code,footprint=fp,structure=s,axes=[3,1], &
        modes=[character(len=16):: 'grid-wrap','constant'],cval=-2.0_dp)
    call checked_real()
    do m=1,size(mode_names)
        call grey_erosion_nd(x,y,code,structure=s,axes=[3,1],mode=trim(mode_names(m)))
        call checked_real()
        call grey_dilation_nd(x,y,code,structure=s,axes=[3,1],mode=trim(mode_names(m)))
        call checked_real()
        call grey_erosion_nd(x,y,code,structure=s,axes=[3,1],origin=[-1,0], &
            mode=trim(mode_names(m)))
        call checked_real()
        call grey_dilation_nd(x,y,code,structure=s,axes=[3,1],origin=[-1,0], &
            mode=trim(mode_names(m)))
        call checked_real()
    end do
    close(u)

    open(newunit=u,file='generic_axes_output.bin',access='stream',form='unformatted',status='replace')
    call generic_filter_nd(x,weighted,y,code,sizes=[2,3],axes=[3,1]); call checked_real()
    call generic_filter_nd(x,weighted,y,code,footprint=fp,axes=[3,1],origin=[1,-1], &
        mode='constant',cval=-1.0_dp); call checked_real()
    call generic_filter_nd(x,weighted,y,code,sizes=[3,2],axes=[1,3],origin=[1,-1], &
        mode='nearest'); call checked_real()
    call generic_filter_nd(x,weighted,y,code,sizes=[integer::],axes=[integer::]); call checked_real()
    call generic_filter_nd(x,weighted,y,code,sizes=[2,3,2],mode='mirror'); call checked_real()
    call generic_filter_nd(x,weighted,y,code,sizes=[2,3],axes=[3,1], &
        modes=[character(len=16):: 'wrap','constant'],cval=-2.0_dp); call checked_real()
    do m=1,size(mode_names)
        call generic_filter_nd(x,weighted,y,code,footprint=fp,axes=[3,1], &
            mode=trim(mode_names(m))); call checked_real()
    end do
    close(u)

    open(newunit=u,file='generic_axes_complex_output.bin',access='stream',form='unformatted',status='replace')
    call generic_filter_nd(z,weighted_complex,cy,code,sizes=[2,3],axes=[3,1]); call checked_complex()
    call generic_filter_nd(z,weighted_complex,cy,code,footprint=fp,axes=[3,1], &
        mode='constant',cval=cmplx(2.0_dp,-.7_dp,dp)); call checked_complex()
    call generic_filter_nd(z,weighted_complex,cy,code,sizes=[integer::],axes=[integer::]); call checked_complex()
    call generic_filter_nd(z,weighted_complex,cy,code,sizes=[2,3],axes=[3,1], &
        modes=[character(len=16):: 'grid-wrap','constant'],cval=cmplx(-2.0_dp,.5_dp,dp))
    call checked_complex()
    close(u)
contains
    subroutine checked_real()
        if (code/=0) error stop 'gray/generic axes reference rejected valid input'
        write(u) y
    end subroutine checked_real
    subroutine checked_complex()
        if (code/=0) error stop 'complex generic axes reference rejected valid input'
        write(u) cy
    end subroutine checked_complex
    function weighted(v) result(r)
        real(dp), intent(in) :: v(:)
        real(dp) :: r
        integer :: n
        r=0.0_dp
        do n=1,size(v)
            r=r+real(n,dp)*v(n)
        end do
    end function weighted
    function weighted_complex(v) result(r)
        complex(dp), intent(in) :: v(:)
        complex(dp) :: r
        integer :: n
        r=cmplx(0.0_dp,0.0_dp,dp)
        do n=1,size(v)
            r=r+real(n,dp)*v(n)
        end do
    end function weighted_complex
end program gray_generic_axes_driver
