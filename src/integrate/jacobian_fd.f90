! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Utilities for finite-difference Jacobian evaluation.

module scifort_integrate_jacobian_fd
    implicit none
    private

    public :: group_fd_columns

contains

    pure subroutine group_fd_columns(sparsity, groups, n_groups)
        logical, intent(in) :: sparsity(:, :) !! nonzero Jacobian pattern, shape (n,n)
        integer, intent(out) :: groups(:) !! greedy perturbation group for each column; 0 if empty
        integer, intent(out) :: n_groups !! number of nonempty perturbation groups

        logical :: group_rows(size(sparsity, 1), size(sparsity, 2))
        integer :: g
        integer :: j

        groups = 0
        group_rows = .false.
        n_groups = 0
        do j = 1, size(sparsity, 2)
            if (.not. any(sparsity(:, j))) cycle
            do g = 1, n_groups
                if (.not. any(group_rows(:, g) .and. sparsity(:, j))) then
                    groups(j) = g
                    group_rows(:, g) = group_rows(:, g) .or. sparsity(:, j)
                    exit
                end if
            end do
            if (groups(j) == 0) then
                n_groups = n_groups + 1
                groups(j) = n_groups
                group_rows(:, n_groups) = sparsity(:, j)
            end if
        end do
    end subroutine group_fd_columns

end module scifort_integrate_jacobian_fd
