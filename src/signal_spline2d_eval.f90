! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Evaluation of tensor-product quadratic and cubic spline coefficients.
! Mirrors scipy.ndimage.map_coordinates(..., mode="mirror", prefilter=.false.)
! on an existing coefficient array. Coordinates use Fortran one-based indexing.
module scifort_signal_spline2d_eval
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spline, only : spline_status_success, spline_status_invalid_input
    implicit none
    private
    public :: cspline2d_eval, qspline2d_eval
    public :: cspline2d_grid, qspline2d_grid

    interface cspline2d_eval
        module procedure c_eval_real_vector
        module procedure c_eval_real_matrix
        module procedure c_eval_real_scalar
        module procedure c_eval_complex_vector
        module procedure c_eval_complex_matrix
        module procedure c_eval_complex_scalar
    end interface

    interface cspline2d_grid
        module procedure c_grid_real_matrix
        module procedure c_grid_complex_matrix
    end interface

    interface qspline2d_eval
        module procedure q_eval_real_vector
        module procedure q_eval_real_matrix
        module procedure q_eval_real_scalar
        module procedure q_eval_complex_vector
        module procedure q_eval_complex_matrix
        module procedure q_eval_complex_scalar
    end interface

    interface qspline2d_grid
        module procedure q_grid_real_matrix
        module procedure q_grid_complex_matrix
    end interface

contains

    pure logical function valid_real(coeff) result(ok)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficients
        ok = size(coeff, 1) > 0 .and. size(coeff, 2) > 0
        if (ok) ok = all(ieee_is_finite(coeff))
    end function valid_real

    pure logical function valid_complex(coeff) result(ok)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficients
        ok = size(coeff, 1) > 0 .and. size(coeff, 2) > 0
        if (ok) ok = all(ieee_is_finite(real(coeff, dp))) .and. &
            all(ieee_is_finite(aimag(coeff)))
    end function valid_complex

    pure subroutine spline_support(x, n, degree, indices, weights, count)
        real(dp), intent(in) :: x !! one-based sample coordinate
        integer, intent(in) :: n !! number of coefficients
        integer, intent(in) :: degree !! quadratic (2) or cubic (3)
        integer, intent(out) :: indices(4) !! reflected support indices
        real(dp), intent(out) :: weights(4) !! nonnegative spline weights
        integer, intent(out) :: count !! three or four spline neighbors
        real(dp) :: u, t, z, period
        integer :: center, k, raw

        indices = 1
        weights = 0.0_dp
        ! Reflect continuous coordinates before calculating compact support.
        ! For n=1, every coordinate samples the sole coefficient.
        if (n == 1) then
            u = 0.0_dp
        else
            period = 2.0_dp * real(n - 1, dp)
            u = modulo(x - 1.0_dp, period)
            if (u > real(n - 1, dp)) u = period - u
        end if
        if (degree == 3) then
            center = floor(u)
            t = u - real(center, dp)
            weights(1) = ((1.0_dp - t)**3) / 6.0_dp
            weights(2) = (4.0_dp - 6.0_dp*t*t + 3.0_dp*t*t*t) / 6.0_dp
            weights(3) = (1.0_dp + 3.0_dp*t + 3.0_dp*t*t - 3.0_dp*t*t*t) / 6.0_dp
            weights(4) = (t*t*t) / 6.0_dp
            count = 4
            do k = 1, count
                raw = center + k - 2
                indices(k) = mirror_index(raw, n)
            end do
        else
            center = floor(u + 0.5_dp)
            count = 3
            do k = 1, count
                raw = center + k - 2
                z = abs(u - real(raw, dp))
                if (z < 0.5_dp) then
                    weights(k) = 0.75_dp - z*z
                else if (z < 1.5_dp) then
                    weights(k) = 0.5_dp*(1.5_dp-z)**2
                end if
                indices(k) = mirror_index(raw, n)
            end do
        end if
    end subroutine spline_support

    pure integer function mirror_index(raw, n) result(index)
        integer, intent(in) :: raw !! zero-based source index
        integer, intent(in) :: n !! axis length
        integer :: t, period
        index = 1
        if (n == 1) return
        period = 2*(n-1)
        t = modulo(raw, period)
        if (t >= n) t = period - t
        index = t + 1
    end function mirror_index

    pure real(dp) function eval_real(coeff, row, col, degree) result(value)
        real(dp), intent(in) :: coeff(:,:) !! real spline coefficient matrix
        real(dp), intent(in) :: row !! one-based row coordinate
        real(dp), intent(in) :: col !! one-based column coordinate
        integer, intent(in) :: degree !! 2 or 3
        integer :: ri(4), ci(4), nr, nc, i, j
        real(dp) :: rw(4), cw(4)
        call spline_support(row, size(coeff,1), degree, ri, rw, nr)
        call spline_support(col, size(coeff,2), degree, ci, cw, nc)
        value = 0.0_dp
        do j = 1, nc
            do i = 1, nr
                value = value + rw(i)*cw(j)*coeff(ri(i),ci(j))
            end do
        end do
    end function eval_real

    pure complex(dp) function eval_complex(coeff, row, col, degree) result(value)
        complex(dp), intent(in) :: coeff(:,:) !! complex spline coefficient matrix
        real(dp), intent(in) :: row !! one-based row coordinate
        real(dp), intent(in) :: col !! one-based column coordinate
        integer, intent(in) :: degree !! 2 or 3
        integer :: ri(4), ci(4), nr, nc, i, j
        real(dp) :: rw(4), cw(4)
        call spline_support(row, size(coeff,1), degree, ri, rw, nr)
        call spline_support(col, size(coeff,2), degree, ci, cw, nc)
        value = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, nc
            do i = 1, nr
                value = value + rw(i)*cw(j)*coeff(ri(i),ci(j))
            end do
        end do
    end function eval_complex

    function c_eval_real_scalar(coeff, rows, cols, status) result(values)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows !! one-based row coordinates
        real(dp), intent(in) :: cols !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        real(dp) :: values
        integer :: code

        code = spline_status_invalid_input
        values = 0.0_dp
        if (.not. valid_real(coeff)) goto 900
        if ((.not. ieee_is_finite(rows))) goto 900
        if ((.not. ieee_is_finite(cols))) goto 900
        values = eval_real(coeff, rows, cols, 3)
        code = spline_status_success
900     if (present(status)) status = code
    end function c_eval_real_scalar

    function c_eval_real_vector(coeff, rows, cols, status) result(values)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: values(:)
        integer :: code, i

        code = spline_status_invalid_input
        allocate(values(0))
        if (.not. valid_real(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        if (any(shape(rows) /= shape(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows)))
        do i = 1, size(rows)
            values(i) = eval_real(coeff, rows(i), cols(i), 3)
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function c_eval_real_vector

    function c_eval_real_matrix(coeff, rows, cols, status) result(values)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:,:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:,:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: values(:,:)
        integer :: code, i, j

        code = spline_status_invalid_input
        allocate(values(0,0))
        if (.not. valid_real(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        if (any(shape(rows) /= shape(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows,1),size(rows,2)))
        do j = 1, size(rows,2)
            do i = 1, size(rows,1)
                values(i,j) = eval_real(coeff, rows(i,j), cols(i,j), 3)
            end do
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function c_eval_real_matrix

    function c_grid_real_matrix(coeff, rows, cols, status) result(values)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: values(:,:)
        integer :: code, i, j

        code = spline_status_invalid_input
        allocate(values(0,0))
        if (.not. valid_real(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows),size(cols)))
        do j = 1, size(cols)
            do i = 1, size(rows)
                values(i,j) = eval_real(coeff, rows(i), cols(j), 3)
            end do
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function c_grid_real_matrix

    function c_eval_complex_scalar(coeff, rows, cols, status) result(values)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows !! one-based row coordinates
        real(dp), intent(in) :: cols !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        complex(dp) :: values
        integer :: code

        code = spline_status_invalid_input
        values = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_complex(coeff)) goto 900
        if ((.not. ieee_is_finite(rows))) goto 900
        if ((.not. ieee_is_finite(cols))) goto 900
        values = eval_complex(coeff, rows, cols, 3)
        code = spline_status_success
900     if (present(status)) status = code
    end function c_eval_complex_scalar

    function c_eval_complex_vector(coeff, rows, cols, status) result(values)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: values(:)
        integer :: code, i

        code = spline_status_invalid_input
        allocate(values(0))
        if (.not. valid_complex(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        if (any(shape(rows) /= shape(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows)))
        do i = 1, size(rows)
            values(i) = eval_complex(coeff, rows(i), cols(i), 3)
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function c_eval_complex_vector

    function c_eval_complex_matrix(coeff, rows, cols, status) result(values)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:,:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:,:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: values(:,:)
        integer :: code, i, j

        code = spline_status_invalid_input
        allocate(values(0,0))
        if (.not. valid_complex(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        if (any(shape(rows) /= shape(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows,1),size(rows,2)))
        do j = 1, size(rows,2)
            do i = 1, size(rows,1)
                values(i,j) = eval_complex(coeff, rows(i,j), cols(i,j), 3)
            end do
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function c_eval_complex_matrix

    function c_grid_complex_matrix(coeff, rows, cols, status) result(values)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: values(:,:)
        integer :: code, i, j

        code = spline_status_invalid_input
        allocate(values(0,0))
        if (.not. valid_complex(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows),size(cols)))
        do j = 1, size(cols)
            do i = 1, size(rows)
                values(i,j) = eval_complex(coeff, rows(i), cols(j), 3)
            end do
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function c_grid_complex_matrix

    function q_eval_real_scalar(coeff, rows, cols, status) result(values)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows !! one-based row coordinates
        real(dp), intent(in) :: cols !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        real(dp) :: values
        integer :: code

        code = spline_status_invalid_input
        values = 0.0_dp
        if (.not. valid_real(coeff)) goto 900
        if ((.not. ieee_is_finite(rows))) goto 900
        if ((.not. ieee_is_finite(cols))) goto 900
        values = eval_real(coeff, rows, cols, 2)
        code = spline_status_success
900     if (present(status)) status = code
    end function q_eval_real_scalar

    function q_eval_real_vector(coeff, rows, cols, status) result(values)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: values(:)
        integer :: code, i

        code = spline_status_invalid_input
        allocate(values(0))
        if (.not. valid_real(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        if (any(shape(rows) /= shape(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows)))
        do i = 1, size(rows)
            values(i) = eval_real(coeff, rows(i), cols(i), 2)
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function q_eval_real_vector

    function q_eval_real_matrix(coeff, rows, cols, status) result(values)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:,:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:,:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: values(:,:)
        integer :: code, i, j

        code = spline_status_invalid_input
        allocate(values(0,0))
        if (.not. valid_real(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        if (any(shape(rows) /= shape(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows,1),size(rows,2)))
        do j = 1, size(rows,2)
            do i = 1, size(rows,1)
                values(i,j) = eval_real(coeff, rows(i,j), cols(i,j), 2)
            end do
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function q_eval_real_matrix

    function q_grid_real_matrix(coeff, rows, cols, status) result(values)
        real(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: values(:,:)
        integer :: code, i, j

        code = spline_status_invalid_input
        allocate(values(0,0))
        if (.not. valid_real(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows),size(cols)))
        do j = 1, size(cols)
            do i = 1, size(rows)
                values(i,j) = eval_real(coeff, rows(i), cols(j), 2)
            end do
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function q_grid_real_matrix

    function q_eval_complex_scalar(coeff, rows, cols, status) result(values)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows !! one-based row coordinates
        real(dp), intent(in) :: cols !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        complex(dp) :: values
        integer :: code

        code = spline_status_invalid_input
        values = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_complex(coeff)) goto 900
        if ((.not. ieee_is_finite(rows))) goto 900
        if ((.not. ieee_is_finite(cols))) goto 900
        values = eval_complex(coeff, rows, cols, 2)
        code = spline_status_success
900     if (present(status)) status = code
    end function q_eval_complex_scalar

    function q_eval_complex_vector(coeff, rows, cols, status) result(values)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: values(:)
        integer :: code, i

        code = spline_status_invalid_input
        allocate(values(0))
        if (.not. valid_complex(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        if (any(shape(rows) /= shape(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows)))
        do i = 1, size(rows)
            values(i) = eval_complex(coeff, rows(i), cols(i), 2)
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function q_eval_complex_vector

    function q_eval_complex_matrix(coeff, rows, cols, status) result(values)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:,:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:,:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: values(:,:)
        integer :: code, i, j

        code = spline_status_invalid_input
        allocate(values(0,0))
        if (.not. valid_complex(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        if (any(shape(rows) /= shape(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows,1),size(rows,2)))
        do j = 1, size(rows,2)
            do i = 1, size(rows,1)
                values(i,j) = eval_complex(coeff, rows(i,j), cols(i,j), 2)
            end do
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function q_eval_complex_matrix

    function q_grid_complex_matrix(coeff, rows, cols, status) result(values)
        complex(dp), intent(in) :: coeff(:,:) !! spline coefficient image
        real(dp), intent(in) :: rows(:) !! one-based row coordinates
        real(dp), intent(in) :: cols(:) !! one-based column coordinates
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: values(:,:)
        integer :: code, i, j

        code = spline_status_invalid_input
        allocate(values(0,0))
        if (.not. valid_complex(coeff)) goto 900
        if (any(.not. ieee_is_finite(rows))) goto 900
        if (any(.not. ieee_is_finite(cols))) goto 900
        deallocate(values)
        allocate(values(size(rows),size(cols)))
        do j = 1, size(cols)
            do i = 1, size(rows)
                values(i,j) = eval_complex(coeff, rows(i), cols(j), 2)
            end do
        end do
        code = spline_status_success
900     if (present(status)) status = code
    end function q_grid_complex_matrix

end module scifort_signal_spline2d_eval
