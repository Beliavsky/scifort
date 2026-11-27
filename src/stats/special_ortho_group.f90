! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special_ortho_group
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state
    use scifort_random_matrix_helpers, only : real_determinant
    use scifort_ortho_group, only : ortho_group_rvs, ortho_group_status_success
    implicit none
    private

    integer, parameter, public :: special_ortho_group_status_success = 0
    integer, parameter, public :: special_ortho_group_status_invalid_shape = 1
    integer, parameter, public :: special_ortho_group_status_numerical_failure = 2

    public :: special_ortho_group_rvs
    public :: special_ortho_group_rvs_array

contains

    subroutine special_ortho_group_rvs(state, sample, status)
        !! Draw a Haar-distributed special-orthogonal matrix.
        !! @dimension n = size(sample,1)
        !! @require n >= 0
        !! Empty square matrices succeed without advancing state.
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the orthogonal draw
        real(dp), intent(out) :: sample(:, :) !! shape (n,n); Haar-distributed special-orthogonal matrix
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: info
        real(dp) :: det

        if (size(sample,1) /= size(sample,2)) then
            sample = 0.0_dp
            if (present(status)) status = special_ortho_group_status_invalid_shape
            return
        end if
        call ortho_group_rvs(state, sample, info)
        if (info /= ortho_group_status_success) then
            if (present(status)) status = special_ortho_group_status_numerical_failure
            return
        end if
        if (size(sample,1) > 0) then
            det = real_determinant(sample)
            if (det < 0.0_dp) sample(1,:) = -sample(1,:)
        end if
        if (present(status)) status = special_ortho_group_status_success
    end subroutine special_ortho_group_rvs

    subroutine special_ortho_group_rvs_array(state, samples, status)
        !! Draw a batch of special-orthogonal Haar matrices.
        !! @dimension n = size(samples,1)
        !! @require n >= 0
        !! @dimension n_samples = size(samples,3)
        !! @require n_samples >= 0
        !! Empty square matrices or batches succeed without advancing state.
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        real(dp), intent(out) :: samples(:, :, :) !! shape (n,n,n_samples); Haar matrices along the third axis
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: info
        integer :: k

        if (size(samples,1) /= size(samples,2)) then
            samples = 0.0_dp
            if (present(status)) status = special_ortho_group_status_invalid_shape
            return
        end if
        do k = 1, size(samples,3)
            call special_ortho_group_rvs(state, samples(:,:,k), info)
            if (info /= special_ortho_group_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = special_ortho_group_status_success
    end subroutine special_ortho_group_rvs_array

end module scifort_special_ortho_group
