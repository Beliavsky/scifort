! SPDX-License-Identifier: MIT
! SciPy numerical reference stream for complex line filters and axis morphology.
program complex_line_binary_axes_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : generic_filter1d_nd, binary_dilation_nd, &
        binary_erosion_nd, binary_opening_nd, binary_closing_nd
    implicit none
    complex(dp) :: a2(5,6),b2(5,6),a3(3,4,2),b3(3,4,2)
    logical :: m3(4,5,3),result3(4,5,3),mask3(4,5,3)
    logical :: m4(2,3,4,2),result4(2,3,4,2),mask4(2,3,4,2)
    logical :: s23(2,3),s33(3,3),s4(4),s323(3,2,3)
    character(len=16), parameter :: all_modes(8)=[character(len=16) :: &
        'reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']
    integer :: unit,p,q,op,s,i,j,k,l
    complex(dp), parameter :: cval=cmplx(1.5_dp,-0.75_dp,dp)

    open(newunit=unit,file='complex_line_binary_axes_results.bin', &
        access='stream',form='unformatted',status='replace')
    do j=1,6
        do i=1,5
            a2(i,j)=cmplx(sin(0.37_dp*i+0.13_dp*j),cos(0.41_dp*j-0.19_dp*i),dp)
        end do
    end do
    do k=1,2
        do j=1,4
            do i=1,3
                a3(i,j,k)=cmplx(sin(0.14_dp*i*j+0.7_dp*k),cos(0.21_dp*i-0.31_dp*j*k),dp)
            end do
        end do
    end do
    do p=1,8
        call generic_filter1d_nd(a2,weighted_line,4,1,b2,s, &
            mode=trim(all_modes(p)),cval=cval,origin=1)
        call check(s,'complex 2D axis1')
        write(unit) b2
        call generic_filter1d_nd(a2,weighted_line,3,2,b2,s, &
            mode=trim(all_modes(p)),cval=cval,origin=-1)
        call check(s,'complex 2D axis2')
        write(unit) b2
        call generic_filter1d_nd(a3,weighted_line,2,3,b3,s, &
            mode=trim(all_modes(p)),cval=cval)
        call check(s,'complex 3D axis3')
        write(unit) b3
        call generic_filter1d_nd(a3,weighted_line,5,2,b3,s, &
            mode=trim(all_modes(p)),cval=cval,origin=1)
        call check(s,'complex 3D axis2')
        write(unit) b3
    end do

    do k=1,3
        do j=1,5
            do i=1,4
                m3(i,j,k)=mod(3*i+5*j+7*k,11)<4
                mask3(i,j,k)=mod(i+2*j+3*k,5)>0
            end do
        end do
    end do
    do l=1,2
        do k=1,4
            do j=1,3
                do i=1,2
                    m4(i,j,k,l)=mod(i+2*j+3*k+5*l,7)<3
                    mask4(i,j,k,l)=mod(2*i+j+k+3*l,4)/=0
                end do
            end do
        end do
    end do
    s23=reshape([.true.,.false.,.true.,.true.,.true.,.false.],shape(s23))
    s33=reshape([.true.,.true.,.false.,.false.,.true.,.true.,.true.,.false.,.true.],shape(s33))
    s4=[.true.,.true.,.false.,.true.]
    s323=reshape([(mod(i,4)/=0,i=1,size(s323))],shape(s323))
    do op=1,4
        do q=1,9
            select case(q)
            case(1)
                call morph3(op,s,axes=[1,3],structure=s23,origin=[-1,1], &
                    iterations=2,mask=mask3)
            case(2)
                call morph3(op,s,axes=[3,1],structure=s33,origin=[1,-1], &
                    iterations=1,border_value=.true.)
            case(3)
                call morph3(op,s,axes=[2],structure=s4,origin=[1], &
                    iterations=1,mask=mask3)
            case(4)
                call morph3(op,s,axes=[1,2,3],structure=s323,origin=[1,0,-1], &
                    iterations=1,border_value=.true.)
            case(5)
                call morph3(op,s,axes=[integer::])
            case(6)
                call morph3(op,s,axes=[2,3],connectivity=1,origin=[0,1],iterations=2)
            case(7)
                call morph3(op,s,axes=[1,2],connectivity=2,iterations=0)
            case(8)
                call morph3(op,s,axes=[1,2,3],iterations=0,mask=mask3)
            case(9)
                call morph3(op,s,axes=[1,3],origin=[-1,0])
            end select
            call check(s,'3D binary morphology')
            write(unit) merge(1,0,result3)
        end do
        call morph4(op,s,axes=[1,4],structure=s23,origin=[0,0],iterations=2,mask=mask4)
        call check(s,'4D binary morphology')
        write(unit) merge(1,0,result4)
    end do
    close(unit)
    print '(a)', 'PASS: complex line callbacks and selected-axis binary morphology reference driver'
contains
    subroutine check(code,name)
        integer, intent(in) :: code
        character(len=*), intent(in) :: name
        if (code/=0) then
            print '(3a,i0)', 'FAIL ',name,': status=',code
            error stop 1
        end if
    end subroutine check
    subroutine weighted_line(input_line,output_line)
        complex(dp), intent(in) :: input_line(:)
        complex(dp), intent(out) :: output_line(:)
        integer :: j,t,width
        width=size(input_line)-size(output_line)+1
        output_line=(0.0_dp,0.0_dp)
        do j=1,size(output_line)
            do t=1,width
                output_line(j)=output_line(j)+ &
                    cmplx(0.2_dp*t,0.03_dp*(-1)**t*t,dp)*input_line(j+t-1)
            end do
        end do
    end subroutine weighted_line
    subroutine morph3(op,status,axes,structure,origin,iterations,border_value,mask,connectivity)
        integer, intent(in) :: op
        integer, intent(out) :: status
        integer, optional, intent(in) :: axes(:)
        logical, optional, intent(in) :: structure(..)
        integer, optional, intent(in) :: origin(:)
        integer, optional, intent(in) :: iterations
        logical, optional, intent(in) :: border_value
        logical, optional, intent(in) :: mask(..)
        integer, optional, intent(in) :: connectivity
        select case(op)
        case(1)
            call binary_dilation_nd(m3,result3,status,connectivity,iterations, &
                border_value,structure,origin,mask,axes)
        case(2)
            call binary_erosion_nd(m3,result3,status,connectivity,iterations, &
                border_value,structure,origin,mask,axes)
        case(3)
            call binary_opening_nd(m3,result3,status,connectivity,iterations, &
                border_value,structure,origin,mask,axes)
        case(4)
            call binary_closing_nd(m3,result3,status,connectivity,iterations, &
                border_value,structure,origin,mask,axes)
        end select
    end subroutine morph3
    subroutine morph4(op,status,axes,structure,origin,iterations,mask)
        integer, intent(in) :: op
        integer, intent(out) :: status
        integer, intent(in) :: axes(:)
        logical, intent(in) :: structure(..)
        integer, intent(in) :: origin(:)
        integer, intent(in) :: iterations
        logical, intent(in) :: mask(..)
        select case(op)
        case(1)
            call binary_dilation_nd(m4,result4,status,structure=structure, &
                origin=origin,iterations=iterations,mask=mask,axes=axes)
        case(2)
            call binary_erosion_nd(m4,result4,status,structure=structure, &
                origin=origin,iterations=iterations,mask=mask,axes=axes)
        case(3)
            call binary_opening_nd(m4,result4,status,structure=structure, &
                origin=origin,iterations=iterations,mask=mask,axes=axes)
        case(4)
            call binary_closing_nd(m4,result4,status,structure=structure, &
                origin=origin,iterations=iterations,mask=mask,axes=axes)
        end select
    end subroutine morph4
end program complex_line_binary_axes_driver
