! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_stats, only : normal_cdf, normal_pdf, normal_ppf
    implicit none
    private

    integer(c_int), parameter :: status_ok = 0_c_int
    integer(c_int), parameter :: status_invalid_argument = 1_c_int

    public :: scifort_normal_cdf_f64
    public :: scifort_normal_cdf_vec_f64
    public :: scifort_normal_pdf_f64
    public :: scifort_normal_pdf_vec_f64
    public :: scifort_normal_ppf_f64
    public :: scifort_normal_ppf_vec_f64
    public :: scifort_version_major
    public :: scifort_version_minor
    public :: scifort_version_patch

contains

    function scifort_version_major() result(version) bind(c, name="scifort_version_major")
        integer(c_int) :: version

        version = 0_c_int
    end function scifort_version_major

    function scifort_version_minor() result(version) bind(c, name="scifort_version_minor")
        integer(c_int) :: version

        version = 1_c_int
    end function scifort_version_minor

    function scifort_version_patch() result(version) bind(c, name="scifort_version_patch")
        integer(c_int) :: version

        version = 0_c_int
    end function scifort_version_patch

    function scifort_normal_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_normal_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double) :: y

        y = real(normal_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_normal_pdf_f64

    function scifort_normal_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_normal_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double) :: y

        y = real(normal_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_normal_cdf_f64

    function scifort_normal_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_normal_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double) :: x

        x = real(normal_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_normal_ppf_f64

    subroutine scifort_normal_pdf_vec_f64(n, x, loc, scale, y, status) &
            bind(c, name="scifort_normal_pdf_vec_f64")
        integer(c_size_t), value, intent(in) :: n !! number of elements of x and y
        real(c_double), intent(in) :: x(*) !! points of evaluation
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double), intent(out) :: y(*) !! densities; unchanged if status is nonzero
        integer(c_int), intent(out) :: status !! 0 on success, 1 for an invalid loc or scale

        integer(c_size_t) :: i

        if (.not. valid_c_loc_scale(loc, scale)) then
            status = status_invalid_argument
            return
        end if

        do i = 1_c_size_t, n
            y(i) = real(normal_pdf(real(x(i), dp), real(loc, dp), &
                real(scale, dp)), c_double)
        end do
        status = status_ok
    end subroutine scifort_normal_pdf_vec_f64

    subroutine scifort_normal_cdf_vec_f64(n, x, loc, scale, y, status) &
            bind(c, name="scifort_normal_cdf_vec_f64")
        integer(c_size_t), value, intent(in) :: n !! number of elements of x and y
        real(c_double), intent(in) :: x(*) !! points of evaluation
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double), intent(out) :: y(*) !! lower-tail probabilities; unchanged if status /= 0
        integer(c_int), intent(out) :: status !! 0 on success, 1 for an invalid loc or scale

        integer(c_size_t) :: i

        if (.not. valid_c_loc_scale(loc, scale)) then
            status = status_invalid_argument
            return
        end if

        do i = 1_c_size_t, n
            y(i) = real(normal_cdf(real(x(i), dp), real(loc, dp), &
                real(scale, dp)), c_double)
        end do
        status = status_ok
    end subroutine scifort_normal_cdf_vec_f64

    subroutine scifort_normal_ppf_vec_f64(n, p, loc, scale, x, status) &
            bind(c, name="scifort_normal_ppf_vec_f64")
        integer(c_size_t), value, intent(in) :: n !! number of elements of p and x
        real(c_double), intent(in) :: p(*) !! lower-tail probabilities
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double), intent(out) :: x(*) !! quantiles; unchanged if status is nonzero
        integer(c_int), intent(out) :: status !! 0 on success, 1 for an invalid loc or scale

        integer(c_size_t) :: i

        if (.not. valid_c_loc_scale(loc, scale)) then
            status = status_invalid_argument
            return
        end if

        do i = 1_c_size_t, n
            x(i) = real(normal_ppf(real(p(i), dp), real(loc, dp), &
                real(scale, dp)), c_double)
        end do
        status = status_ok
    end subroutine scifort_normal_ppf_vec_f64

    pure logical function valid_c_loc_scale(loc, scale) result(valid)
        real(c_double), intent(in) :: loc !! mean to check
        real(c_double), intent(in) :: scale !! standard deviation to check

        valid = ieee_is_finite(loc) .and. ieee_is_finite(scale) .and. &
            scale > 0.0_c_double
    end function valid_c_loc_scale

end module scifort_c_api
