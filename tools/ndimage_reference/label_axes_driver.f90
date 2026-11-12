! SPDX-License-Identifier: MIT
! Write deterministic selected-axis component labels for SciPy slice checks.
program label_axes_driver
    use scifort_ndimage, only : label_nd
    implicit none
    logical :: mask3(6,4,7),mask4(2,4,3,5),footprint(3,3)
    integer :: labels3(6,4,7),labels4(2,4,3,5),i,j,k,l,s,n,unit

    do k=1,7
        do j=1,4
            do i=1,6
                mask3(i,j,k)=modulo(2*i+3*j+7*k,9)<4
            end do
        end do
    end do
    do l=1,5
        do k=1,3
            do j=1,4
                do i=1,2
                    mask4(i,j,k,l)=modulo(2*i+3*j+k+5*l,11)<5
                end do
            end do
        end do
    end do
    footprint=reshape([.true.,.false.,.true.,.false.,.true.,.false., &
        .true.,.false.,.true.],[3,3])
    open(newunit=unit,file='selected_axis_labels.bin',access='stream', &
        form='unformatted',status='replace')
    call label_nd(mask3,labels3,n,s,axes=[3,1],structure=footprint)
    if (s/=0) error stop '3D axis labels'
    write(unit) n
    write(unit) labels3
    call label_nd(mask4,labels4,n,s,axes=[4,2],connectivity=1)
    if (s/=0) error stop '4D axis labels'
    write(unit) n
    write(unit) labels4
    call label_nd(mask3,labels3,n,s,axes=[integer::])
    if (s/=0) error stop 'empty axes labels'
    write(unit) n
    write(unit) labels3
    close(unit)
    print '(a)', 'PASS: selected-axis labels driver'
end program label_axes_driver
