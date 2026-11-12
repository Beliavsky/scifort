! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Directed and symmetric Hausdorff distances for finite point clouds.
! Independent implementation of the mathematical max-min definition,
! with safe scaled Euclidean distances and deterministic early abandonment.
module scifort_spatial_hausdorff
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : spatial_distance_success, spatial_distance_invalid_input
    implicit none
    private
    public :: directed_hausdorff, hausdorff_distance

contains

    pure function euclidean_distance_scaled(a, b) result(d)
        real(dp), intent(in) :: a(:) !! coordinates of the first point
        real(dp), intent(in) :: b(:) !! coordinates of the second point
        real(dp) :: d
        real(dp) :: delta(size(a)), scale
        delta = a-b
        scale = maxval(abs(delta))
        if (scale == 0.0_dp) then
            d = 0.0_dp
        else if (.not. ieee_is_finite(scale)) then
            ! Both points can be finite while a-b exceeds the binary64 range.
            d = ieee_value(0.0_dp, ieee_positive_inf)
        else
            d = scale*sqrt(sum((delta/scale)**2))
        end if
    end function euclidean_distance_scaled

    function directed_hausdorff(u, v, index_u, index_v, status, early_break) result(d)
        real(dp), intent(in) :: u(:,:) !! source points (observations, dimensions)
        real(dp), intent(in) :: v(:,:) !! target points (observations, dimensions)
        integer, intent(out), optional :: index_u !! 1-based source witness index, 0 on error
        integer, intent(out), optional :: index_v !! 1-based target witness index, 0 on error
        integer, intent(out), optional :: status !! zero if valid; one for invalid input
        logical, intent(in), optional :: early_break !! abandon non-maximal points, default true
        real(dp) :: d
        real(dp) :: nearest, candidate
        integer :: i, j, near_index, witness_u, witness_v, code
        logical :: prune

        d = ieee_value(0.0_dp, ieee_quiet_nan)
        code = spatial_distance_invalid_input
        witness_u = 0
        witness_v = 0
        if (size(u,1) == 0 .or. size(v,1) == 0) goto 900
        if (size(u,2) == 0 .or. size(u,2) /= size(v,2)) goto 900
        if (.not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v))) goto 900

        prune = .true.
        if (present(early_break)) prune = early_break
        d = -1.0_dp
        do i = 1, size(u,1)
            nearest = ieee_value(0.0_dp, ieee_positive_inf)
            near_index = 1
            do j = 1, size(v,1)
                candidate = euclidean_distance_scaled(u(i,:), v(j,:))
                if (j == 1 .or. candidate < nearest) then
                    nearest = candidate
                    near_index = j
                end if
                ! The nearest distance of this row cannot increase as j grows.
                ! A row already bounded by d cannot become the first maximizer.
                if (prune .and. nearest <= d) exit
            end do
            if (nearest > d) then
                d = nearest
                witness_u = i
                witness_v = near_index
            end if
        end do
        code = spatial_distance_success
900     if (present(index_u)) index_u = witness_u
        if (present(index_v)) index_v = witness_v
        if (present(status)) status = code
    end function directed_hausdorff

    function hausdorff_distance(u, v, index_u, index_v, from_u, status) result(d)
        real(dp), intent(in) :: u(:,:) !! first point cloud (observations, dimensions)
        real(dp), intent(in) :: v(:,:) !! second point cloud (observations, dimensions)
        integer, intent(out), optional :: index_u !! 1-based witness index in u; 0 on error
        integer, intent(out), optional :: index_v !! 1-based witness index in v; 0 on error
        logical, intent(out), optional :: from_u !! true if directed u-to-v distance maximizes
        integer, intent(out), optional :: status !! zero if valid; one for invalid input
        real(dp) :: d
        real(dp) :: d_uv, d_vu
        integer :: i_uv, j_uv, j_vu, i_vu, code
        logical :: forward

        d_uv = directed_hausdorff(u, v, i_uv, j_uv, code)
        if (code == spatial_distance_success) then
            d_vu = directed_hausdorff(v, u, j_vu, i_vu, code)
        end if
        forward = .true.
        if (code /= spatial_distance_success) then
            d = ieee_value(0.0_dp, ieee_quiet_nan)
            i_uv = 0
            j_uv = 0
            forward = .false.
        else if (d_vu > d_uv) then
            d = d_vu
            i_uv = i_vu
            j_uv = j_vu
            forward = .false.
        else
            d = d_uv
        end if
        if (present(index_u)) index_u = i_uv
        if (present(index_v)) index_v = j_uv
        if (present(from_u)) from_u = forward
        if (present(status)) status = code
    end function hausdorff_distance
end module scifort_spatial_hausdorff
