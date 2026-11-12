! SPDX-License-Identifier: MIT
! Compare rank-revealing FIFO QR downdates against legacy span reconstruction.
program test_optimize_krylov_augmentation_fifo_qr
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        reset_lgmres_direction_pruning_stream, retain_lgmres_direction, &
        retain_lgmres_direction_stream, retain_lgmres_direction_fifo
    implicit none

    integer, parameter :: n = 96
    integer, parameter :: capacity = 12
    integer, parameter :: ntrial = 240
    real(dp), parameter :: tolerance = 1.0e-8_dp
    real(dp) :: old_vectors(n, capacity)
    real(dp) :: new_vectors(n, capacity)
    real(dp) :: vec(n)
    type(lgmres_pruning_workspace) :: work
    integer :: old_count
    integer :: new_count
    integer :: failures
    integer :: i
    integer :: j
    integer :: accepted_after_fill
    logical :: added_old
    logical :: added_new

    failures = 0
    old_vectors = 0.0_dp
    new_vectors = 0.0_dp
    old_count = 0
    new_count = 0
    accepted_after_fill = 0
    call reset_lgmres_direction_pruning_stream(work)
    do j = 1, ntrial
        do i = 1, n
            vec(i) = sin(real(17*i*i + 19*i*j + 23*j*j, dp)) + &
                0.7_dp * cos(real(7*i*i + 29*i*j + 3*j*j, dp))
        end do
        if (mod(j, 13) == 0) vec = new_vectors(:, max(1, new_count))
        if (mod(j, 29) == 0) vec = 0.0_dp
        if (mod(j, 37) == 0) vec = vec * (huge(1.0_dp) / 8.0_dp)
        if (mod(j, 41) == 0) vec = vec * (tiny(1.0_dp) / 8.0_dp)
        call retain_lgmres_direction(old_vectors, old_count, vec, &
            tolerance, added_old)
        if (new_count < capacity) then
            call retain_lgmres_direction_stream(new_vectors, new_count, &
                vec, tolerance, added_new, work)
        else
            call retain_lgmres_direction_fifo(new_vectors, new_count, &
                vec, tolerance, added_new, work)
            if (added_new) accepted_after_fill = accepted_after_fill + 1
        end if
        call check(added_new .eqv. added_old, 'matching retention decision')
        call check(new_count == old_count, 'matching active size')
        call check(all(new_vectors == old_vectors), 'matching retained directions')
    end do
    call check(accepted_after_fill > 100, 'many FIFO evictions performed')

    ! Reset must discard the old QR history before another solve starts.
    old_vectors = 0.0_dp
    new_vectors = 0.0_dp
    old_count = 0
    new_count = 0
    call reset_lgmres_direction_pruning_stream(work)
    do j = ntrial, 1, -1
        do i = 1, n
            vec(i) = sin(real(11*i*i + 31*i*j + 13*j*j, dp))
        end do
        call retain_lgmres_direction(old_vectors, old_count, vec, &
            tolerance, added_old)
        if (new_count < capacity) then
            call retain_lgmres_direction_stream(new_vectors, new_count, &
                vec, tolerance, added_new, work)
        else
            call retain_lgmres_direction_fifo(new_vectors, new_count, &
                vec, tolerance, added_new, work)
        end if
        call check(added_new .eqv. added_old, 'reset retention decision')
        call check(all(new_vectors == old_vectors), 'reset retained history')
    end do

    ! A nearly singular full history must safely fall back to reconstruction.
    old_vectors = 0.0_dp
    new_vectors = 0.0_dp
    old_count = 0
    new_count = 0
    call reset_lgmres_direction_pruning_stream(work)
    do j = 1, capacity
        vec = 0.0_dp
        vec(1) = 1.0_dp
        if (j > 1) vec(j) = 1.0e-15_dp
        call retain_lgmres_direction(old_vectors, old_count, vec, &
            1.0e-18_dp, added_old)
        call retain_lgmres_direction_stream(new_vectors, new_count, &
            vec, 1.0e-18_dp, added_new, work)
        call check(added_old .eqv. added_new, 'near-rank-deficient initial history')
    end do
    vec = 0.0_dp
    vec(n) = 1.0_dp
    call retain_lgmres_direction(old_vectors, old_count, vec, &
        1.0e-18_dp, added_old)
    call retain_lgmres_direction_fifo(new_vectors, new_count, vec, &
        1.0e-18_dp, added_new, work)
    call check(added_old .eqv. added_new, 'fallback accepts new direction')
    call check(all(new_vectors == old_vectors), 'fallback retains original semantics')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_krylov_augmentation_fifo_qr: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_augmentation_fifo_qr: PASS'
contains
    subroutine check(condition, description)
        logical, intent(in) :: condition !! required result
        character(len=*), intent(in) :: description !! failure context
        if (.not. condition) then
            failures = failures + 1
            if (failures <= 25) print '(a)', 'FAIL: '//description
        end if
    end subroutine check
end program test_optimize_krylov_augmentation_fifo_qr
