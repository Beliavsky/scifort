! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Great-circle interpolation on the unit sphere in arbitrary dimensions.
! Independent implementation of the geometric spherical interpolation formula.
module scifort_spatial_slerp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : spatial_distance_success, spatial_distance_invalid_input
    implicit none
    private
    public :: geometric_slerp

    interface geometric_slerp
        module procedure geometric_slerp_scalar
        module procedure geometric_slerp_array
    end interface geometric_slerp

contains

    pure subroutine prepare_arc(start, finish, tol, x, tangent, angle, valid)
        real(dp), intent(in) :: start(:) !! starting unit-sphere vector
        real(dp), intent(in) :: finish(:) !! ending unit-sphere vector
        real(dp), intent(in), optional :: tol !! antipodal chord threshold, default 1e-7
        real(dp), intent(out) :: x(:) !! normalized starting coordinate vector
        real(dp), intent(out) :: tangent(:) !! unit tangent towards finish, if distinct
        real(dp), intent(out) :: angle !! shortest arc angle in radians, from 0 to pi
        logical, intent(out) :: valid !! true for valid unit points and unique shortest arc
        real(dp) :: na, nb, dot_product, tangent_norm, antipodal_tol

        valid = .false.
        angle = 0.0_dp
        x = 0.0_dp
        tangent = 0.0_dp
        if (size(start) < 2 .or. size(finish) /= size(start)) return
        if (.not. all(ieee_is_finite(start)) .or. .not. all(ieee_is_finite(finish))) return
        antipodal_tol = 1.0e-7_dp
        if (present(tol)) antipodal_tol = tol
        if (.not. ieee_is_finite(antipodal_tol) .or. antipodal_tol < 0.0_dp) return
        na = norm2(start)
        nb = norm2(finish)
        if (.not. ieee_is_finite(na) .or. .not. ieee_is_finite(nb)) return
        if (abs(na-1.0_dp) > 1.0e-7_dp .or. abs(nb-1.0_dp) > 1.0e-7_dp) return
        x = start/na
        tangent = finish/nb
        if (norm2(x+tangent) <= antipodal_tol) return
        dot_product = max(-1.0_dp, min(1.0_dp, sum(x*tangent)))
        tangent = tangent - dot_product*x
        tangent_norm = norm2(tangent)
        angle = atan2(tangent_norm, dot_product)
        if (tangent_norm > 0.0_dp) tangent = tangent/tangent_norm
        valid = .true.
    end subroutine prepare_arc

    pure function point_on_arc(x, tangent, angle, fraction, finish) result(point)
        real(dp), intent(in) :: x(:) !! normalized starting vector
        real(dp), intent(in) :: tangent(:) !! tangent defining the great circle
        real(dp), intent(in) :: angle !! arc angle, radians
        real(dp), intent(in) :: fraction !! arc parameter (extrapolation allowed)
        real(dp), intent(in) :: finish(:) !! exact supplied end vector
        real(dp) :: point(size(x))
        if (fraction == 0.0_dp) then
            point = x
        else if (fraction == 1.0_dp) then
            point = finish
        else if (angle == 0.0_dp) then
            point = x
        else
            point = cos(fraction*angle)*x + sin(fraction*angle)*tangent
        end if
    end function point_on_arc

    function geometric_slerp_scalar(start, finish, t, tol, status) result(point)
        real(dp), intent(in) :: start(:) !! starting point on a unit sphere
        real(dp), intent(in) :: finish(:) !! ending point on a unit sphere
        real(dp), intent(in) :: t !! path fraction, any finite value (may extrapolate)
        real(dp), intent(in), optional :: tol !! antipodal chord threshold, default 1e-7
        integer, intent(out), optional :: status !! zero on success; one on invalid input
        real(dp), allocatable :: point(:)
        real(dp) :: x(size(start)), tangent(size(start)), theta
        logical :: valid
        integer :: code

        call prepare_arc(start, finish, tol, x, tangent, theta, valid)
        code = spatial_distance_invalid_input
        if (valid .and. ieee_is_finite(t)) then
            if (ieee_is_finite(t*theta)) then
                allocate(point(size(start)))
                point = point_on_arc(x, tangent, theta, t, finish/norm2(finish))
                code = spatial_distance_success
            end if
        end if
        if (code /= spatial_distance_success) allocate(point(0))
        if (present(status)) status = code
    end function geometric_slerp_scalar

    function geometric_slerp_array(start, finish, t, tol, status) result(points)
        real(dp), intent(in) :: start(:) !! starting point on a unit sphere
        real(dp), intent(in) :: finish(:) !! ending point on a unit sphere
        real(dp), intent(in) :: t(:) !! finite path fractions, may be unsorted/extrapolated
        real(dp), intent(in), optional :: tol !! antipodal chord threshold, default 1e-7
        integer, intent(out), optional :: status !! zero on success; one on invalid input
        real(dp), allocatable :: points(:,:)
        real(dp) :: x(size(start)), tangent(size(start)), end_unit(size(start)), theta
        logical :: valid
        integer :: code, i

        call prepare_arc(start, finish, tol, x, tangent, theta, valid)
        code = spatial_distance_invalid_input
        if (valid .and. all(ieee_is_finite(t))) then
            if (all(ieee_is_finite(t*theta))) then
                allocate(points(size(t),size(start)))
                end_unit = finish/norm2(finish)
                do i=1,size(t)
                    points(i,:) = point_on_arc(x, tangent, theta, t(i), end_unit)
                end do
                code = spatial_distance_success
            end if
        end if
        if (code /= spatial_distance_success) allocate(points(0,0))
        if (present(status)) status = code
    end function geometric_slerp_array
end module scifort_spatial_slerp
