! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program lebedev_reference_driver
    use scifort_integrate_lebedev, only : lebedev_rule
    use scifort_kinds, only : dp
    implicit none
    integer, parameter :: orders(*)=[3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35, &
        41,47,53,59,65,71,77,83,89,95,101,107,113,119,125,131]
    real(dp), allocatable :: x(:,:),w(:)
    integer :: i,j,status
    do i=1,size(orders)
        call lebedev_rule(orders(i),x,w,status)
        if (status/=0) error stop 1
        print '(2(i0,1x))',orders(i),size(w)
        do j=1,size(w)
            print '(4es26.17e3)',x(:,j),w(j)
        end do
    end do
end program lebedev_reference_driver
