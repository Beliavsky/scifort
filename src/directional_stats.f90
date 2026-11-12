! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_directional_stats
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    type, public :: directional_stats_result
        real(dp), allocatable :: mean_direction(:) !! unit directional mean
        real(dp) :: mean_resultant_length = 0.0_dp !! norm of the sample mean vector
    end type directional_stats_result

    type, public :: directional_stats_batch_result
        real(dp), allocatable :: mean_direction(:, :) !! batched unit directional means
        real(dp), allocatable :: mean_resultant_length(:) !! batched mean-resultant lengths
    end type directional_stats_batch_result

    interface directional_stats
        module procedure directional_stats_rank2
        module procedure directional_stats_rank3
    end interface directional_stats

    public :: directional_stats

contains

    pure function directional_stats_rank2(samples, axis, normalize) result(res)
        real(dp), intent(in) :: samples(:, :) !! vector observations
        integer, intent(in), optional :: axis !! one-based observation axis, default 1
        logical, intent(in), optional :: normalize !! normalize each observation, default true
        type(directional_stats_result) :: res

        integer :: axis_value
        integer :: i
        integer :: n_obs
        integer :: n_vec
        logical :: normalize_value
        real(dp), allocatable :: mean(:)
        real(dp), allocatable :: unit(:)
        real(dp), allocatable :: vector(:)
        real(dp) :: norm_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        normalize_value = .true.
        if (present(normalize)) normalize_value = normalize

        if (axis_value < 1 .or. axis_value > 2) then
            allocate(res%mean_direction(0))
            res%mean_resultant_length = quiet_nan(0.0_dp)
            return
        end if

        if (axis_value == 1) then
            n_obs = size(samples, 1)
            n_vec = size(samples, 2)
        else
            n_obs = size(samples, 2)
            n_vec = size(samples, 1)
        end if

        allocate(res%mean_direction(n_vec), mean(n_vec), unit(n_vec), vector(n_vec))
        if (n_obs == 0 .or. n_vec == 0) then
            res%mean_direction = quiet_nan(0.0_dp)
            res%mean_resultant_length = quiet_nan(0.0_dp)
            return
        end if

        mean = 0.0_dp
        do i = 1, n_obs
            if (axis_value == 1) then
                vector = samples(i, :)
            else
                vector = samples(:, i)
            end if
            if (normalize_value) then
                call normalize_vector(vector, unit, norm_value)
                mean = mean + unit
            else
                mean = mean + vector
            end if
        end do
        mean = mean / real(n_obs, dp)

        call normalize_vector(mean, res%mean_direction, res%mean_resultant_length)
    end function directional_stats_rank2

    pure function directional_stats_rank3(samples, axis, normalize) result(res)
        real(dp), intent(in) :: samples(:, :, :) !! batched vector observations
        integer, intent(in), optional :: axis !! one-based observation axis, default 1
        logical, intent(in), optional :: normalize !! normalize each observation, default true
        type(directional_stats_batch_result) :: res

        integer :: axis_value
        integer :: i
        integer :: j
        integer :: k
        logical :: normalize_value
        real(dp), allocatable :: mean(:)
        real(dp), allocatable :: unit(:)
        real(dp), allocatable :: vector(:)
        real(dp) :: norm_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        normalize_value = .true.
        if (present(normalize)) normalize_value = normalize

        select case (axis_value)
        case (1)
            allocate(res%mean_direction(size(samples, 2), size(samples, 3)))
            allocate(res%mean_resultant_length(size(samples, 2)))
            allocate(mean(size(samples, 3)), unit(size(samples, 3)), vector(size(samples, 3)))
            if (size(samples, 1) == 0 .or. size(samples, 3) == 0) then
                res%mean_direction = quiet_nan(0.0_dp)
                res%mean_resultant_length = quiet_nan(0.0_dp)
                return
            end if
            do j = 1, size(samples, 2)
                mean = 0.0_dp
                do i = 1, size(samples, 1)
                    vector = samples(i, j, :)
                    if (normalize_value) then
                        call normalize_vector(vector, unit, norm_value)
                        mean = mean + unit
                    else
                        mean = mean + vector
                    end if
                end do
                mean = mean / real(size(samples, 1), dp)
                call normalize_vector(mean, res%mean_direction(j, :), &
                    res%mean_resultant_length(j))
            end do

        case (2)
            allocate(res%mean_direction(size(samples, 1), size(samples, 3)))
            allocate(res%mean_resultant_length(size(samples, 1)))
            allocate(mean(size(samples, 3)), unit(size(samples, 3)), vector(size(samples, 3)))
            if (size(samples, 2) == 0 .or. size(samples, 3) == 0) then
                res%mean_direction = quiet_nan(0.0_dp)
                res%mean_resultant_length = quiet_nan(0.0_dp)
                return
            end if
            do i = 1, size(samples, 1)
                mean = 0.0_dp
                do j = 1, size(samples, 2)
                    vector = samples(i, j, :)
                    if (normalize_value) then
                        call normalize_vector(vector, unit, norm_value)
                        mean = mean + unit
                    else
                        mean = mean + vector
                    end if
                end do
                mean = mean / real(size(samples, 2), dp)
                call normalize_vector(mean, res%mean_direction(i, :), &
                    res%mean_resultant_length(i))
            end do

        case (3)
            allocate(res%mean_direction(size(samples, 1), size(samples, 2)))
            allocate(res%mean_resultant_length(size(samples, 1)))
            allocate(mean(size(samples, 2)), unit(size(samples, 2)), vector(size(samples, 2)))
            if (size(samples, 3) == 0 .or. size(samples, 2) == 0) then
                res%mean_direction = quiet_nan(0.0_dp)
                res%mean_resultant_length = quiet_nan(0.0_dp)
                return
            end if
            do i = 1, size(samples, 1)
                mean = 0.0_dp
                do k = 1, size(samples, 3)
                    vector = samples(i, :, k)
                    if (normalize_value) then
                        call normalize_vector(vector, unit, norm_value)
                        mean = mean + unit
                    else
                        mean = mean + vector
                    end if
                end do
                mean = mean / real(size(samples, 3), dp)
                call normalize_vector(mean, res%mean_direction(i, :), &
                    res%mean_resultant_length(i))
            end do

        case default
            allocate(res%mean_direction(0, 0), res%mean_resultant_length(0))
        end select
    end function directional_stats_rank3

    pure subroutine normalize_vector(vector, unit, norm_value)
        real(dp), intent(in) :: vector(:) !! vector to normalize
        real(dp), intent(out) :: unit(:) !! normalized vector or NaNs if undefined
        real(dp), intent(out) :: norm_value !! Euclidean norm

        integer :: i
        real(dp) :: scale
        real(dp) :: sumsq

        if (size(unit) /= size(vector)) then
            norm_value = quiet_nan(0.0_dp)
            return
        end if
        if (size(vector) == 0) then
            norm_value = 0.0_dp
            return
        end if
        if (any(ieee_is_nan(vector))) then
            norm_value = quiet_nan(0.0_dp)
            unit = quiet_nan(0.0_dp)
            return
        end if

        scale = maxval(abs(vector))
        if (.not. ieee_is_finite(scale)) then
            norm_value = positive_infinity(0.0_dp)
            do i = 1, size(vector)
                if (ieee_is_finite(vector(i))) then
                    unit(i) = 0.0_dp
                else
                    unit(i) = quiet_nan(0.0_dp)
                end if
            end do
            return
        end if
        if (scale == 0.0_dp) then
            norm_value = 0.0_dp
            unit = quiet_nan(0.0_dp)
            return
        end if

        sumsq = sum((vector / scale) ** 2)
        norm_value = scale * sqrt(sumsq)
        if (norm_value == 0.0_dp .or. .not. ieee_is_finite(norm_value)) then
            unit = quiet_nan(0.0_dp)
            return
        end if
        unit = vector / norm_value
    end subroutine normalize_vector

end module scifort_directional_stats
