! SPDX-License-Identifier: MIT
! Binary stream of selected-axis propagation and fill-holes results.
program propagation_axes_driver
    use scifort_ndimage, only : binary_propagation_nd, binary_fill_holes_nd
    implicit none
    logical :: seeds(7,4,7),mask(7,4,7),image(7,4,7),y(7,4,7),foot(3,3)
    integer :: i,j,k,u,s

    do k=1,7
        do j=1,4
            do i=1,7
                seeds(i,j,k)=modulo(7*i+11*j+3*k,23)==0
                mask(i,j,k)=modulo(i+3*j+5*k,13)/=0
            end do
        end do
    end do
    foot=reshape([.false.,.true.,.false.,.true.,.true.,.true., &
        .false.,.true.,.false.],shape(foot))
    image=.true.
    do j=1,4
        image(2:6,j,2:6)=.false.
        image(4,j,4)=.true.
        image(1,j,3)=.false.
    end do
    open(newunit=u,file='propagation_axes.bin',status='replace',form='unformatted', &
        access='stream',action='write')
    call binary_propagation_nd(seeds,mask,y,s,axes=[3,1],structure=foot)
    call check(s)
    write(u) y
    call binary_propagation_nd(seeds,mask,y,s,axes=[2])
    call check(s)
    write(u) y
    call binary_propagation_nd(seeds,mask,y,s,axes=[integer::])
    call check(s)
    write(u) y
    call binary_fill_holes_nd(image,y,s,axes=[1,3],structure=foot)
    call check(s)
    write(u) y
    call binary_fill_holes_nd(image,y,s,axes=[2])
    call check(s)
    write(u) y
    call binary_fill_holes_nd(image,y,s,axes=[integer::])
    call check(s)
    write(u) y
    close(u)
    print '(a)', 'PASS: axis-selective segmentation reference stream'
contains
    subroutine check(status)
        integer, intent(in) :: status !! zero for success
        if (status/=0) error stop 'invalid selected-axis segmentation'
    end subroutine check
end program propagation_axes_driver
