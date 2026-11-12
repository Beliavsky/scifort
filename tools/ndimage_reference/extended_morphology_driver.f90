! SPDX-License-Identifier: MIT
! Reference cases for opening, closing, hole filling, and propagation.
program extended_morphology_driver
    use scifort_ndimage, only : binary_opening_nd, binary_closing_nd, &
        binary_fill_holes_nd, binary_propagation_nd
    implicit none
    logical :: m2(5,6),m3(4,5,3),r2(5,6),r3(4,5,3)
    logical :: s2(5,6),s3(4,5,3),a2(5,6),a3(4,5,3)
    logical :: h2(7,8),h3(5,5,4),out2(7,8),out3(5,5,4)
    integer :: i,j,k,conn,it,b,s,u

    do j=1,6
        do i=1,5
            m2(i,j)=modulo(7*i+11*j+i*j,8)<3
            s2(i,j)=m2(i,j).and.(modulo(3*i+2*j,4)==0)
            a2(i,j)=modulo(2*i+5*j,7)<5
        end do
    end do
    do k=1,3
        do j=1,5
            do i=1,4
                m3(i,j,k)=modulo(3*i+7*j+5*k+i*j+k*j,9)<4
                s3(i,j,k)=m3(i,j,k).and.(modulo(i+3*j+7*k,5)==0)
                a3(i,j,k)=modulo(2*i+j+3*k,7)<5
            end do
        end do
    end do
    h2=.true.
    h2(3:5,3:6)=.false.
    h2(4,5)=.true.
    h2(1,4)=.false.
    h2(2,4)=.false.
    h3=.true.
    h3(2:4,2:4,2:3)=.false.
    h3(1,1,1)=.false.
    open(newunit=u,file='extended_morphology_results.bin',form='unformatted', &
         access='stream',status='replace',action='write')
    do conn=1,2
        do it=1,2
            do b=0,1
                call binary_opening_nd(m2,r2,s,connectivity=conn,iterations=it,border_value=b==1)
                if (s/=0) error stop '2d opening'
                write(u) merge(1,0,r2)
                call binary_closing_nd(m2,r2,s,connectivity=conn,iterations=it,border_value=b==1)
                if (s/=0) error stop '2d closing'
                write(u) merge(1,0,r2)
            end do
        end do
        call binary_fill_holes_nd(h2,out2,s,connectivity=conn)
        if (s/=0) error stop '2d fill holes'
        write(u) merge(1,0,out2)
        call binary_propagation_nd(s2,a2,r2,s,connectivity=conn)
        if (s/=0) error stop '2d propagation'
        write(u) merge(1,0,r2)
    end do
    do conn=1,3
        do it=1,2
            do b=0,1
                call binary_opening_nd(m3,r3,s,connectivity=conn,iterations=it,border_value=b==1)
                if (s/=0) error stop '3d opening'
                write(u) merge(1,0,r3)
                call binary_closing_nd(m3,r3,s,connectivity=conn,iterations=it,border_value=b==1)
                if (s/=0) error stop '3d closing'
                write(u) merge(1,0,r3)
            end do
        end do
        call binary_fill_holes_nd(h3,out3,s,connectivity=conn)
        if (s/=0) error stop '3d fill holes'
        write(u) merge(1,0,out3)
        call binary_propagation_nd(s3,a3,r3,s,connectivity=conn)
        if (s/=0) error stop '3d propagation'
        write(u) merge(1,0,r3)
    end do
    close(u)
end program extended_morphology_driver
