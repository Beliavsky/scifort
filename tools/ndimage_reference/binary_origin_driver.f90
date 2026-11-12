! SPDX-License-Identifier: MIT
program binary_origin_driver
    use scifort_ndimage, only : binary_dilation_nd, binary_erosion_nd, &
        binary_opening_nd, binary_closing_nd, binary_hit_or_miss_nd
    implicit none
    logical :: x(6,7),y(6,7),a(3,3),b(2,3)
    integer :: i,j,case_id,op,st,orig(2),unit
    do j=1,7
        do i=1,6
            x(i,j)=modulo(i*i+2*j*i+3*j,7)<3
        end do
    end do
    a=reshape([.true.,.false.,.false.,.true.,.true.,.false.,.false.,.true.,.true.],[3,3])
    b=reshape([.true.,.false.,.false.,.true.,.true.,.false.],[2,3])
    open(newunit=unit,file='binary_origin_results.bin',form='unformatted', &
        access='stream',status='replace')
    do case_id=1,4
        select case(case_id)
        case(1)
            orig=[0,0]
        case(2)
            orig=[-1,1]
        case(3)
            orig=[1,-1]
        case(4)
            orig=[-1,0]
        end select
        do op=1,5
            if (case_id==4) then
                select case(op)
                case(1)
                    call binary_dilation_nd(x,y,st,structure=b,origin=orig)
                case(2)
                    call binary_erosion_nd(x,y,st,structure=b,origin=orig)
                case(3)
                    call binary_opening_nd(x,y,st,structure=b,origin=orig)
                case(4)
                    call binary_closing_nd(x,y,st,structure=b,origin=orig)
                case(5)
                    call binary_hit_or_miss_nd(x,y,st,structure1=b,origin1=orig)
                end select
            else
                select case(op)
                case(1)
                    call binary_dilation_nd(x,y,st,structure=a,origin=orig)
                case(2)
                    call binary_erosion_nd(x,y,st,structure=a,origin=orig)
                case(3)
                    call binary_opening_nd(x,y,st,structure=a,origin=orig)
                case(4)
                    call binary_closing_nd(x,y,st,structure=a,origin=orig)
                case(5)
                    call binary_hit_or_miss_nd(x,y,st,structure1=a,origin1=orig)
                end select
            end if
            if (st/=0) error stop 'binary origin status'
            write(unit) merge(1,0,y)
        end do
    end do
    close(unit)
    print '(a)', 'PASS: binary origin reference driver'
end program binary_origin_driver
