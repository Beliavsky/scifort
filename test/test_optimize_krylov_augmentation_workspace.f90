! SPDX-License-Identifier: MIT
! The optional reusable pruning workspace must match allocation-per-call paths.
program test_optimize_krylov_augmentation_workspace
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : &
        lgmres_pruning_workspace, retain_lgmres_direction, retain_lgmres_image_direction
    implicit none

    type(lgmres_pruning_workspace) :: scratch, shared_scratch
    real(dp), allocatable :: fresh(:, :), cached(:, :), shared(:, :), images(:, :)
    real(dp), allocatable :: direction(:), jv(:), arnoldi_basis(:, :)
    integer :: n, capacity, count_a, count_b, count_c, k, mode, index, count_before
    integer :: failures
    logical :: accepted_a, accepted_b, accepted_c

    failures = 0
    capacity = 5
    ! Reuse one workspace over changes in dimension and retained-basis size.
    do index = 1, 4
        select case (index)
        case (1)
            n = 4
        case (2)
            n = 16
        case (3)
            n = 3
        case default
            n = 24
        end select
        allocate(fresh(n, capacity), cached(n, capacity), shared(n, capacity), &
            images(n, capacity), direction(n), jv(n), arnoldi_basis(n, capacity))
        do mode = 1, 3
            fresh = 0.0_dp
            cached = 0.0_dp
            shared = 0.0_dp
            images = 0.0_dp
            count_a = 0
            count_b = 0
            count_c = 0
            do k = 1, 4 * capacity + 3
                direction = 0.0_dp
                jv = 0.0_dp
                ! Vary scales and include exact duplicates and near dependence.
                direction(1 + mod(k - 1, n)) = 1.0_dp
                if (mod(k, 4) == 0) direction(1) = direction(1) + 1.0e-10_dp
                jv = direction
                if (k >= capacity) jv(1) = jv(1) + 1.0_dp
                if (mod(k, 7) == 0) then
                    direction = 0.0_dp
                    jv = 0.0_dp
                end if
                count_before = count_a
                ! Retired Arnoldi columns may contain unrelated old data.
                arnoldi_basis = -22.0_dp
                select case (mode)
                case (1)
                    call retain_lgmres_image_direction(fresh, count_a, direction, &
                        images, jv, 1.0e-8_dp, accepted_a)
                    call retain_lgmres_image_direction(cached, count_b, direction, &
                        images, jv, 1.0e-8_dp, accepted_b, workspace=scratch)
                    call retain_lgmres_image_direction(shared, count_c, direction, &
                        images, jv, 1.0e-8_dp, accepted_c, workspace=shared_scratch, &
                        basis_work=arnoldi_basis)
                case (2)
                    call retain_lgmres_image_direction(fresh, count_a, direction, &
                        images, jv, 0.0_dp, accepted_a, auto_threshold=.true.)
                    call retain_lgmres_image_direction(cached, count_b, direction, &
                        images, jv, 0.0_dp, accepted_b, auto_threshold=.true., workspace=scratch)
                    call retain_lgmres_image_direction(shared, count_c, direction, &
                        images, jv, 0.0_dp, accepted_c, auto_threshold=.true., &
                        workspace=shared_scratch, basis_work=arnoldi_basis)
                case (3)
                    call retain_lgmres_direction(fresh, count_a, direction, &
                        1.0e-8_dp, accepted_a)
                    call retain_lgmres_direction(cached, count_b, direction, &
                        1.0e-8_dp, accepted_b, workspace=scratch)
                    call retain_lgmres_direction(shared, count_c, direction, &
                        1.0e-8_dp, accepted_c, workspace=shared_scratch, &
                        basis_work=arnoldi_basis)
                end select
                call check(accepted_a .eqv. accepted_b, 'matching insertion flag')
                call check(accepted_a .eqv. accepted_c, 'shared basis insertion flag')
                call check(count_a == count_b .and. count_a == count_c, &
                    'matching retained count')
                call check(all(fresh == cached), 'bitwise-identical retained directions')
                call check(all(fresh == shared), 'shared Arnoldi basis matches')
                if (accepted_a .and. mode < 3) then
                    if (count_before == capacity) then
                        images(:, 1:capacity - 1) = images(:, 2:capacity)
                    end if
                    images(:, count_a) = jv
                end if
            end do
        end do
        deallocate(fresh, cached, shared, images, direction, jv, arnoldi_basis)
    end do

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_krylov_augmentation_workspace: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_augmentation_workspace: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! expected condition
        character(len=*), intent(in) :: description !! test case

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check

end program test_optimize_krylov_augmentation_workspace
