! SPDX-License-Identifier: MIT
! SciPy value_indices reference fixture in Fortran's one-based CSR format.
program value_indices_driver
    use scifort_ndimage, only : value_indices_nd
    implicit none
    integer :: a1(11),a2(9,7),a3(4,5,3),a4(64,48),empty_labels(2,2)
    integer :: high(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: i,j,k,unit

    do i=1,11
        a1(i)=mod(7*i,11)-5
    end do
    do j=1,7
        do i=1,9
            a2(i,j)=mod(13*i+17*j+3*i*j,31)-13
        end do
    end do
    do k=1,3
        do j=1,5
            do i=1,4
                a3(i,j,k)=mod(17*i+31*j+7*k+i*j*k,13)-6
            end do
        end do
    end do
    do j=1,48
        do i=1,64
            a4(i,j)=mod(39*i+17*j+3*i*j,25)-12
        end do
    end do
    empty_labels=1
    high=17
    high(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=-4
    open(newunit=unit,file='value_indices_results.bin', &
        status='replace',form='unformatted',access='stream')
    call output(a1)
    call output(a2,0)
    call output(a3,-1)
    call output(a3,55)
    call output(a4,0)
    call output(empty_labels,1)
    call output(high,17)
    close(unit)
    print '(a)', 'PASS: value_indices Fortran reference driver'
contains
    subroutine output(arr,ignore)
        integer, intent(in) :: arr(..)
        integer, optional, intent(in) :: ignore
        integer, allocatable :: values(:),starts(:),coords(:,:)
        integer :: s
        call value_indices_nd(arr,values,starts,coords,s,ignore)
        if (s/=0) error stop 'value_indices driver status'
        write(unit) [rank(arr),size(values),size(coords,2)]
        write(unit) values
        write(unit) starts
        write(unit) coords
    end subroutine output
end program value_indices_driver
