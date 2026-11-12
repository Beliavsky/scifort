! SPDX-License-Identifier: MIT
! Time repeated image-space LGMRES pruning with and without reusable scratch.
! gfortran -O2 -std=f2018 kinds.f90 krylov_augmentation.f90 \
!     benchmark_lgmres_pruning_workspace.f90 -o benchmark_lgmres_pruning_workspace
program benchmark_lgmres_pruning_workspace
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : &
        lgmres_pruning_workspace, retain_lgmres_image_direction
    implicit none

    type(lgmres_pruning_workspace) :: scratch
    real(dp), allocatable :: stored(:, :), images(:, :), candidate(:), candidate_image(:)
    real(dp) :: elapsed_alloc, elapsed_reuse, start, finish
    character(len=32) :: argument
    integer :: n, nretain, repetitions, j, iteration, count
    logical :: added

    n = 512
    repetitions = 4000
    nretain = 8
    if (command_argument_count() >= 1) then
        call get_command_argument(1, argument)
        read(argument, *) n
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, argument)
        read(argument, *) repetitions
    end if
    if (n < nretain .or. repetitions < 1) error stop 'invalid benchmark parameters'
    allocate(stored(n, nretain + 1), images(n, nretain + 1), &
        candidate(n), candidate_image(n))
    stored = 0.0_dp
    images = 0.0_dp
    candidate = 0.0_dp
    candidate_image = 0.0_dp
    do j = 1, nretain
        stored(j, j) = 1.0_dp
        images(j, j) = 1.0_dp
    end do
    candidate(1:2) = 1.0_dp
    candidate_image(1:2) = 1.0_dp

    count = nretain
    call retain_lgmres_image_direction(stored, count, candidate, images, &
        candidate_image, 1.0e-8_dp, added, workspace=scratch)
    if (added .or. count /= nretain) error stop 'unexpected benchmark setup'
    call cpu_time(start)
    do iteration = 1, repetitions
        call retain_lgmres_image_direction(stored, count, candidate, images, &
            candidate_image, 1.0e-8_dp, added)
        if (added .or. count /= nretain) error stop 'allocation path changed history'
    end do
    call cpu_time(finish)
    elapsed_alloc = finish - start

    call cpu_time(start)
    do iteration = 1, repetitions
        call retain_lgmres_image_direction(stored, count, candidate, images, &
            candidate_image, 1.0e-8_dp, added, workspace=scratch)
        if (added .or. count /= nretain) error stop 'reuse path changed history'
    end do
    call cpu_time(finish)
    elapsed_reuse = finish - start

    print '(a,i0)', 'n = ', n
    print '(a,i0)', 'retained = ', nretain
    print '(a,i0)', 'repetitions = ', repetitions
    print '(a,es12.4)', 'allocation_per_call_cpu_s = ', elapsed_alloc
    print '(a,es12.4)', 'reused_workspace_cpu_s = ', elapsed_reuse
    if (elapsed_reuse > 0.0_dp) &
        print '(a,f10.3)', 'allocation_over_reuse_ratio = ', elapsed_alloc / elapsed_reuse
end program benchmark_lgmres_pruning_workspace
