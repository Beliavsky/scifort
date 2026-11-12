! SPDX-License-Identifier: MIT
! Fortran bridge for compare_rotation_grid_scipy.py; numbers only on stdout.
program rotation_grid_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_spatial_rotation_stack, only : rotation_rng_type, rotation_rng, &
        rotation_stack_type, rotation_stack_random
    use scifort_spatial_rotation_grid, only : rotation_grid_type, rotation_grid_from_quat, &
        rotation_grid_from_mrp, rotation_grid_from_euler, rotation_grid_from_matrix
    implicit none
    type(rotation_rng_type) :: rng
    type(rotation_stack_type) :: samples, mean0, mean1
    type(rotation_grid_type) :: a, b, composed, roundtrip
    real(dp), allocatable :: qa(:,:,:), qb(:,:,:), v(:,:,:), q(:,:,:), &
        rotvec(:,:,:), matrix(:,:,:,:), applied(:,:,:), inverse_applied(:,:,:), &
        composed_q(:,:,:), magnitudes(:,:), mrps(:,:,:), eulers(:,:,:), &
        back_mrp(:,:,:), back_matrix(:,:,:), back_euler(:,:,:)

    real(dp), allocatable :: weights0(:), weights1(:), mq0(:,:), mq1(:,:)
    integer(int64) :: seed
    integer :: num_random, ncases, icase, i, j, n, m, k, l, vn, vm, code, c2

    read(*,*) num_random, seed
    rng = rotation_rng(seed)
    samples = rotation_stack_random(num_random, rng, code)
    write(*,'(i0)') code
    if (code == 0) then
        mq0 = samples%as_quat()
        do i = 1, num_random
            write(*,'(*(es26.17,1x))') mq0(i,:)
        end do
    end if

    read(*,*) ncases
    do icase = 1, ncases
        read(*,*) n, m, k, l, vn, vm
        allocate(qa(n,m,4), qb(k,l,4), v(vn,vm,3), weights0(n), weights1(m))
        do j = 1, m
            do i = 1, n
                read(*,*) qa(i,j,:)
            end do
        end do
        do j = 1, l
            do i = 1, k
                read(*,*) qb(i,j,:)
            end do
        end do
        do j = 1, vm
            do i = 1, vn
                read(*,*) v(i,j,:)
            end do
        end do
        read(*,*) weights0
        read(*,*) weights1
        a = rotation_grid_from_quat(qa,status=code)
        b = rotation_grid_from_quat(qb,status=c2)
        if (code == 0 .and. c2 == 0) then
            q = a%as_quat(canonical=.true.)
            rotvec = a%as_rotvec()
            matrix = a%as_matrix()
            mrps = a%as_mrp()
            eulers = a%as_euler('ZYX',status=code)
            if (code == 0) then
                roundtrip = rotation_grid_from_mrp(mrps,status=code)
                if (code == 0) back_mrp = roundtrip%as_quat()
            end if
            if (code == 0) then
                roundtrip = rotation_grid_from_matrix(matrix,status=code)
                if (code == 0) back_matrix = roundtrip%as_quat()
            end if
            if (code == 0) then
                roundtrip = rotation_grid_from_euler('ZYX',eulers,status=code)
                if (code == 0) back_euler = roundtrip%as_quat()
            end if
            magnitudes = a%magnitude()
            applied = a%apply(v,status=code)
            inverse_applied = a%apply(v,inverse=.true.,status=c2)
            composed = a%compose(b)
            composed_q = composed%as_quat()
            mean0 = a%mean(1,weights0,status=code)
            mean1 = a%mean(2,weights1,status=c2)
            if (code == 0 .and. c2 == 0) then
                mq0 = mean0%as_quat()
                mq1 = mean1%as_quat()
            end if
        end if
        write(*,'(i0)') max(code,c2)
        if (code == 0 .and. c2 == 0) then
            do j = 1, m
                do i = 1, n
                    write(*,'(*(es26.17,1x))') q(i,j,:)
                end do
            end do
            do j = 1, m
                do i = 1, n
                    write(*,'(*(es26.17,1x))') rotvec(i,j,:)
                end do
            end do
            do j = 1, m
                do i = 1, n
                    write(*,'(*(es26.17,1x))') reshape(matrix(i,j,:,:),[9])
                end do
            end do
            do j = 1, m
                do i = 1, n
                    write(*,'(*(es26.17,1x))') mrps(i,j,:)
                end do
            end do
            do j = 1, m
                do i = 1, n
                    write(*,'(*(es26.17,1x))') eulers(i,j,:)
                end do
            end do
            do j = 1, m
                do i = 1, n
                    write(*,'(*(es26.17,1x))') back_mrp(i,j,:)
                end do
            end do
            do j = 1, m
                do i = 1, n
                    write(*,'(*(es26.17,1x))') back_matrix(i,j,:)
                end do
            end do
            do j = 1, m
                do i = 1, n
                    write(*,'(*(es26.17,1x))') back_euler(i,j,:)
                end do
            end do
            do j = 1, m
                write(*,'(*(es26.17,1x))') magnitudes(:,j)
            end do
            do j = 1, size(applied,2)
                do i = 1, size(applied,1)
                    write(*,'(*(es26.17,1x))') applied(i,j,:)
                end do
            end do
            do j = 1, size(inverse_applied,2)
                do i = 1, size(inverse_applied,1)
                    write(*,'(*(es26.17,1x))') inverse_applied(i,j,:)
                end do
            end do
            do j = 1, size(composed_q,2)
                do i = 1, size(composed_q,1)
                    write(*,'(*(es26.17,1x))') composed_q(i,j,:)
                end do
            end do
            do j = 1, size(mq0,1)
                write(*,'(*(es26.17,1x))') mq0(j,:)
            end do
            do i = 1, size(mq1,1)
                write(*,'(*(es26.17,1x))') mq1(i,:)
            end do
        end if
        deallocate(qa,qb,v,weights0,weights1)
    end do
end program rotation_grid_driver
