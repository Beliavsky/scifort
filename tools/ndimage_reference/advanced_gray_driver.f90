! SPDX-License-Identifier: MIT
! Reference output for non-flat grayscale morphology and shifted origins.
program advanced_gray_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : grey_dilation_nd, grey_erosion_nd, &
        grey_opening_nd, grey_closing_nd, morphological_gradient_nd, &
        morphological_laplace_nd, white_tophat_nd, black_tophat_nd
    implicit none
    character(len=16), parameter :: modes(8) = [character(len=16) :: &
        'reflect','mirror','nearest','wrap','constant','grid-constant', &
        'grid-mirror','grid-wrap']
    real(dp) :: x(5,6),y(5,6),structure(2,3)
    logical :: footprint(2,3)
    integer :: i,j,m,op,case_id,code,unit,origin(2)

    do j=1,6
        do i=1,5
            x(i,j)=sin(real(i,dp)*0.31_dp)+real(modulo(i*i+3*j+i*j,11),dp)
        end do
    end do
    structure=reshape([0.2_dp,-0.7_dp,1.3_dp,0.4_dp,-1.6_dp,2.2_dp],[2,3])
    footprint=reshape([.true.,.false.,.true.,.true.,.false.,.true.],[2,3])
    open(newunit=unit,file='advanced_gray_results.bin',form='unformatted', &
         access='stream',status='replace')
    do m=1,8
        do case_id=1,4
            select case(case_id)
            case(1)
                origin=[0,0]
            case(2)
                origin=[-1,1]
            case(3)
                origin=[0,-1]
            case(4)
                origin=[-1,-1]
            end select
            do op=1,8
                select case(op)
                case(1)
                    call grey_dilation_nd(x,y,code,footprint=footprint,structure=structure, &
                        origin=origin,mode=trim(modes(m)),cval=-2.5_dp)
                case(2)
                    call grey_erosion_nd(x,y,code,footprint=footprint,structure=structure, &
                        origin=origin,mode=trim(modes(m)),cval=-2.5_dp)
                case(3)
                    call grey_opening_nd(x,y,code,footprint=footprint,structure=structure, &
                        origin=origin,mode=trim(modes(m)),cval=-2.5_dp)
                case(4)
                    call grey_closing_nd(x,y,code,footprint=footprint,structure=structure, &
                        origin=origin,mode=trim(modes(m)),cval=-2.5_dp)
                case(5)
                    call morphological_gradient_nd(x,y,code,footprint=footprint,structure=structure, &
                        origin=origin,mode=trim(modes(m)),cval=-2.5_dp)
                case(6)
                    call morphological_laplace_nd(x,y,code,footprint=footprint,structure=structure, &
                        origin=origin,mode=trim(modes(m)),cval=-2.5_dp)
                case(7)
                    call white_tophat_nd(x,y,code,footprint=footprint,structure=structure, &
                        origin=origin,mode=trim(modes(m)),cval=-2.5_dp)
                case(8)
                    call black_tophat_nd(x,y,code,footprint=footprint,structure=structure, &
                        origin=origin,mode=trim(modes(m)),cval=-2.5_dp)
                end select
                if (code/=0) error stop 'advanced grayscale status'
                write(unit) y
            end do
        end do
    end do
    close(unit)
    print '(a)', 'PASS: advanced grayscale reference driver'
end program advanced_gray_driver
