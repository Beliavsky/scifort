! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_unitary_group
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state
    use scifort_random_matrix_helpers, only : haar_unitary, random_matrix_status_success
    implicit none
    private

    integer, parameter, public :: unitary_group_status_success = 0
    integer, parameter, public :: unitary_group_status_invalid_shape = 1
    integer, parameter, public :: unitary_group_status_numerical_failure = 2

    public :: unitary_group_rvs
    public :: unitary_group_rvs_array

contains

    subroutine unitary_group_rvs(state, sample, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by complex Gaussian draws
        complex(dp), intent(out) :: sample(:, :) !! square Haar-distributed unitary matrix
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: info

        if (size(sample,1) /= size(sample,2)) then
            sample = cmplx(0.0_dp, 0.0_dp, kind=dp)
            if (present(status)) status = unitary_group_status_invalid_shape
            return
        end if
        call haar_unitary(state, sample, info)
        if (info == random_matrix_status_success) then
            if (present(status)) status = unitary_group_status_success
        else
            if (present(status)) status = unitary_group_status_numerical_failure
        end if
    end subroutine unitary_group_rvs

    subroutine unitary_group_rvs_array(state, samples, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        complex(dp), intent(out) :: samples(:, :, :) !! square matrices indexed by third dimension
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: info
        integer :: k

        if (size(samples,1) /= size(samples,2)) then
            samples = cmplx(0.0_dp, 0.0_dp, kind=dp)
            if (present(status)) status = unitary_group_status_invalid_shape
            return
        end if
        do k = 1, size(samples,3)
            call unitary_group_rvs(state, samples(:,:,k), info)
            if (info /= unitary_group_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = unitary_group_status_success
    end subroutine unitary_group_rvs_array

end module scifort_unitary_group
