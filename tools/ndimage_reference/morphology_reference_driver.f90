! SPDX-License-Identifier: MIT
! Deterministic 2D and 3D morphology reference output.
program morphology_reference_driver
    use scifort_ndimage_morphology, only : binary_dilation_nd, binary_erosion_nd, &
        label_nd
    implicit none
    logical :: m2(5,6), d2(5,6),m3(4,5,3),d3(4,5,3)
    integer :: l2(5,6),l3(4,5,3),i,j,k,c,it,b,n,s,u

    do j=1,6
        do i=1,5
            m2(i,j)=modulo(7*i+11*j+i*j,8)<3
        end do
    end do
    do k=1,3
        do j=1,5
            do i=1,4
                m3(i,j,k)=modulo(3*i+7*j+5*k+i*j+k*j,9)<4
            end do
        end do
    end do
    open(newunit=u,file='morphology_reference_results.bin',access='stream', &
         form='unformatted',status='replace',action='write')
    do c=1,2
        do it=1,2
            do b=0,1
                call binary_dilation_nd(m2,d2,s,connectivity=c,iterations=it,border_value=b==1)
                if (s/=0) error stop '2d dilate'
                write(u) merge(1,0,d2)
                call binary_erosion_nd(m2,d2,s,connectivity=c,iterations=it,border_value=b==1)
                if (s/=0) error stop '2d erode'
                write(u) merge(1,0,d2)
            end do
        end do
        call label_nd(m2,l2,n,s,connectivity=c)
        if (s/=0) error stop '2d label'
        write(u) n,l2
    end do
    do c=1,3
        do it=1,2
            do b=0,1
                call binary_dilation_nd(m3,d3,s,connectivity=c,iterations=it,border_value=b==1)
                if (s/=0) error stop '3d dilate'
                write(u) merge(1,0,d3)
                call binary_erosion_nd(m3,d3,s,connectivity=c,iterations=it,border_value=b==1)
                if (s/=0) error stop '3d erode'
                write(u) merge(1,0,d3)
            end do
        end do
        call label_nd(m3,l3,n,s,connectivity=c)
        if (s/=0) error stop '3d label'
        write(u) n,l3
    end do
    close(u)
end program morphology_reference_driver
