! SPDX-License-Identifier: MIT
! Standalone driver used for independent SciPy rotation-group comparisons.
program rotation_group_driver
    use scifort_kinds, only : dp
    use scifort_spatial_rotation_stack, only : rotation_stack_type, rotation_stack_from_quat
    use scifort_spatial_rotation_group, only : rotation_create_group, rotation_reduce
    implicit none
    character(len=24) :: gname, gaxis, lname, laxis, rname, raxis
    type(rotation_stack_type) :: group, left, right, input, output
    real(dp), allocatable :: q(:,:), qout(:,:)
    integer, allocatable :: il(:), ir(:)
    integer :: num, iter, n, k, status

    read(*,*) num
    do iter = 1, num
        read(*,*) gname, gaxis, lname, laxis, rname, raxis, n
        allocate(q(n,4))
        do k=1,n
            read(*,*) q(k,:)
        end do
        group = rotation_create_group(trim(gname),trim(gaxis),status)
        if (status/=0) error stop 'invalid group in test'
        qout = group%as_quat()
        write(*,'(I0,1X,I0)') group%size(),n
        do k=1,group%size()
            write(*,'(4(ES26.17,1X))') qout(k,:)
        end do
        if (lname/='NONE') left=rotation_create_group(trim(lname),trim(laxis),status)
        if (rname/='NONE') right=rotation_create_group(trim(rname),trim(raxis),status)
        if (lname=='NONE') then
            write(*,'(I0)') 0
        else
            qout = left%as_quat()
            write(*,'(I0)') left%size()
            do k=1,left%size()
                write(*,'(4(ES26.17,1X))') qout(k,:)
            end do
        end if
        if (rname=='NONE') then
            write(*,'(I0)') 0
        else
            qout = right%as_quat()
            write(*,'(I0)') right%size()
            do k=1,right%size()
                write(*,'(4(ES26.17,1X))') qout(k,:)
            end do
        end if
        input=rotation_stack_from_quat(q,status=status)
        if (status/=0) error stop 'input rotation stack invalid'
        if (lname=='NONE' .and. rname=='NONE') then
            output=rotation_reduce(input,index_left=il,index_right=ir,status=status)
        else if (lname=='NONE') then
            output=rotation_reduce(input,right=right,index_left=il,index_right=ir,status=status)
        else if (rname=='NONE') then
            output=rotation_reduce(input,left=left,index_left=il,index_right=ir,status=status)
        else
            output=rotation_reduce(input,left=left,right=right,index_left=il,index_right=ir,status=status)
        end if
        if (status/=0) error stop 'reduction failed'
        qout=output%as_quat()
        do k=1,n
            write(*,'(I0,1X,I0,1X,4(ES26.17,1X))') il(k),ir(k),qout(k,:)
        end do
        deallocate(q)
    end do
end program rotation_group_driver
