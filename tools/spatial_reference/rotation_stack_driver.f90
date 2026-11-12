! SPDX-License-Identifier: MIT
! Stdin/stdout oracle bridge for compare_rotation_stack_scipy.py
program rotation_stack_driver
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type
    use scifort_spatial_rotation_stack, only : rotation_stack_type, &
        rotation_stack_from_quat, rotation_stack_from_mrp, &
        rotation_stack_from_matrix, rotation_stack_from_rotvec
    implicit none
    type(rotation_stack_type) :: a, b, combined, powered, inverse, back
    type(rotation_type) :: average
    real(dp), allocatable :: qa(:,:), qb(:,:), vectors(:,:), weights(:)
    real(dp), allocatable :: q(:,:), rotvec(:,:), mrp(:,:), matrices(:,:,:), applied(:,:), invapplied(:,:)
    real(dp), allocatable :: magnitude(:), compq(:,:), powq(:,:)
    integer :: ncases, n, m, k, icase, i, code, code2, j

    read(*,*) ncases
    do icase = 1, ncases
        read(*,*) n, m, k
        allocate(qa(n,4),qb(k,4),vectors(m,3),weights(n))
        do i = 1, n
            read(*,*) qa(i,:)
        end do
        do i = 1, k
            read(*,*) qb(i,:)
        end do
        do i = 1, m
            read(*,*) vectors(i,:)
        end do
        read(*,*) weights
        a = rotation_stack_from_quat(qa,status=code)
        b = rotation_stack_from_quat(qb,status=code2)
        if (code == 0 .and. code2 == 0) then
            q = a%as_quat(canonical=.true.)
            rotvec = a%as_rotvec()
            mrp = a%as_mrp()
            matrices = a%as_matrix()
            magnitude = a%magnitude()
            applied = a%apply(vectors,status=code)
            invapplied = a%apply(vectors,inverse=.true.,status=code2)
            combined = a%compose(b)
            compq = combined%as_quat()
            powered = a%power(0.7_dp)
            powq = powered%as_quat()
            average = a%mean(weights=weights)
            back = rotation_stack_from_rotvec(rotvec,status=code)
            if (code == 0) back = rotation_stack_from_mrp(mrp,status=code)
            if (code == 0) back = rotation_stack_from_matrix(matrices,status=code)
        end if
        write(*,'(i0)') max(code,code2)
        if (code == 0 .and. code2 == 0) then
            do i = 1, n
                write(*,'(*(es26.17,1x))') q(i,:)
            end do
            do i = 1, n
                write(*,'(*(es26.17,1x))') rotvec(i,:)
            end do
            do i = 1, n
                write(*,'(*(es26.17,1x))') mrp(i,:)
            end do
            write(*,'(*(es26.17,1x))') magnitude
            do i = 1, n
                write(*,'(*(es26.17,1x))') reshape(matrices(i,:,:),[9])
            end do
            do i = 1, size(applied,1)
                write(*,'(*(es26.17,1x))') applied(i,:)
            end do
            do i = 1, size(invapplied,1)
                write(*,'(*(es26.17,1x))') invapplied(i,:)
            end do
            do i = 1, size(compq,1)
                write(*,'(*(es26.17,1x))') compq(i,:)
            end do
            do i = 1, n
                write(*,'(*(es26.17,1x))') powq(i,:)
            end do
            write(*,'(*(es26.17,1x))') average%as_quat()
        end if
        deallocate(qa,qb,vectors,weights)
    end do
end program rotation_stack_driver
