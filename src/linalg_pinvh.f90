! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Moore--Penrose pseudoinverse of a real symmetric or complex Hermitian matrix.
! Eigenvalues are selected by their absolute magnitude so indefinite matrices
! (including negative eigenvalues) are handled without an SVD.
module scifort_linalg_pinvh
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_dense_eigen, only : eigh
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: pinvh

    interface pinvh
        module procedure pinvh_real
        module procedure pinvh_complex
        module procedure pinvh_real_batch
        module procedure pinvh_complex_batch
    end interface pinvh

contains

    subroutine pinvh_real(a, inverse, status, rank, atol, rtol, lower, check_finite, max_iterations)
        real(dp), intent(in) :: a(:, :) !! square real symmetric input (n,n)
        real(dp), intent(out) :: inverse(:, :) !! Moore--Penrose pseudoinverse (n,n), zero on failure
        integer, intent(out) :: status !! success, invalid shape/option, or numerical failure
        integer, intent(out), optional :: rank !! effective numerical rank of the input
        real(dp), intent(in), optional :: atol !! absolute eigenvalue threshold, default 0; nonnegative
        real(dp), intent(in), optional :: rtol !! relative threshold, default n*epsilon(1.0_dp)
        logical, intent(in), optional :: lower !! use only lower triangle, default true
        logical, intent(in), optional :: check_finite !! reject nonfinite input, default true
        integer, intent(in), optional :: max_iterations !! eigensolver iteration limit, default 120

        real(dp), allocatable :: values(:), vectors(:, :), scaled(:, :)
        real(dp) :: absolute_cutoff, relative_cutoff, threshold
        integer :: n, j, eig_status, kept
        logical :: do_check

        inverse = 0.0_dp
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(inverse, 1) /= n .or. size(inverse, 2) /= n) return
        if (.not. valid_tolerances(n, atol, rtol, max_iterations)) return

        do_check = .true.
        if (present(check_finite)) do_check = check_finite
        if (do_check) then
            if (.not. all(ieee_is_finite(a))) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        allocate(values(n), vectors(n,n), scaled(n,n))
        call eigh(a, values, eig_status, z=vectors, lower=lower, max_iterations=max_iterations)
        status = eig_status
        if (status /= linalg_status_success) return

        absolute_cutoff = 0.0_dp
        if (present(atol)) absolute_cutoff = atol
        relative_cutoff = real(n,dp)*epsilon(1.0_dp)
        if (present(rtol)) relative_cutoff = rtol
        threshold = absolute_cutoff + relative_cutoff*maxval(abs(values))
        scaled = 0.0_dp
        kept = 0
        do j = 1, n
            if (abs(values(j)) <= threshold) cycle
            kept = kept+1
            scaled(:,j) = vectors(:,j)/values(j)
        end do
        inverse = matmul(scaled, transpose(vectors))
        if (.not. all(ieee_is_finite(inverse))) then
            inverse = 0.0_dp
            status = linalg_status_no_convergence
            return
        end if
        if (present(rank)) rank = kept
    end subroutine pinvh_real

    subroutine pinvh_complex(a, inverse, status, rank, atol, rtol, lower, check_finite, max_iterations)
        complex(dp), intent(in) :: a(:, :) !! square complex Hermitian input (n,n)
        complex(dp), intent(out) :: inverse(:, :) !! Moore--Penrose pseudoinverse (n,n), zero on failure
        integer, intent(out) :: status !! success, invalid shape/option, or numerical failure
        integer, intent(out), optional :: rank !! effective numerical rank of the input
        real(dp), intent(in), optional :: atol !! absolute eigenvalue threshold, default 0; nonnegative
        real(dp), intent(in), optional :: rtol !! relative threshold, default n*epsilon(1.0_dp)
        logical, intent(in), optional :: lower !! use only lower triangle, default true
        logical, intent(in), optional :: check_finite !! reject nonfinite input, default true
        integer, intent(in), optional :: max_iterations !! eigensolver iteration limit, default 120

        complex(dp), allocatable :: vectors(:, :), scaled(:, :)
        real(dp), allocatable :: values(:)
        real(dp) :: absolute_cutoff, relative_cutoff, threshold
        integer :: n, j, eig_status, kept
        logical :: do_check

        inverse = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(inverse, 1) /= n .or. size(inverse, 2) /= n) return
        if (.not. valid_tolerances(n, atol, rtol, max_iterations)) return

        do_check = .true.
        if (present(check_finite)) do_check = check_finite
        if (do_check) then
            if (.not. all(ieee_is_finite(real(a,dp))) .or. &
                .not. all(ieee_is_finite(aimag(a)))) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        allocate(values(n), vectors(n,n), scaled(n,n))
        call eigh(a, values, eig_status, z=vectors, lower=lower, max_iterations=max_iterations)
        status = eig_status
        if (status /= linalg_status_success) return

        absolute_cutoff = 0.0_dp
        if (present(atol)) absolute_cutoff = atol
        relative_cutoff = real(n,dp)*epsilon(1.0_dp)
        if (present(rtol)) relative_cutoff = rtol
        threshold = absolute_cutoff + relative_cutoff*maxval(abs(values))
        scaled = cmplx(0.0_dp, 0.0_dp, dp)
        kept = 0
        do j = 1, n
            if (abs(values(j)) <= threshold) cycle
            kept = kept+1
            scaled(:,j) = vectors(:,j)/values(j)
        end do
        inverse = matmul(scaled, transpose(conjg(vectors)))
        if (.not. all(ieee_is_finite(real(inverse,dp))) .or. &
            .not. all(ieee_is_finite(aimag(inverse)))) then
            inverse = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        if (present(rank)) rank = kept
    end subroutine pinvh_complex

    subroutine pinvh_real_batch(a, inverse, status, rank, atol, rtol, lower, check_finite, max_iterations)
        real(dp), intent(in) :: a(:, :, :) !! input matrices (batch,n,n), batch>=1
        real(dp), intent(out) :: inverse(:, :, :) !! pseudoinverses (batch,n,n), zero on failure
        integer, intent(out) :: status !! shared status; first batch failure aborts entire call
        integer, intent(out), optional :: rank(:) !! effective rank for each member (batch)
        real(dp), intent(in), optional :: atol !! shared absolute eigenvalue threshold
        real(dp), intent(in), optional :: rtol !! shared relative eigenvalue threshold
        logical, intent(in), optional :: lower !! use lower input triangle, default true
        logical, intent(in), optional :: check_finite !! check all input entries, default true
        integer, intent(in), optional :: max_iterations !! eigensolver iteration limit

        integer :: b, n, nb, member_rank

        inverse = 0.0_dp
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        nb = size(a,1)
        n = size(a,2)
        if (nb < 1 .or. size(a,3) /= n) return
        if (size(inverse,1) /= nb .or. size(inverse,2) /= n .or. size(inverse,3) /= n) return
        if (present(rank)) then
            if (size(rank) /= nb) return
        end if
        do b = 1, nb
            call pinvh_real(a(b,:,:), inverse(b,:,:), status, member_rank, &
                atol, rtol, lower, check_finite, max_iterations)
            if (status /= linalg_status_success) then
                inverse = 0.0_dp
                if (present(rank)) rank = 0
                return
            end if
            if (present(rank)) rank(b) = member_rank
        end do
        status = linalg_status_success
    end subroutine pinvh_real_batch

    subroutine pinvh_complex_batch(a, inverse, status, rank, atol, rtol, lower, check_finite, max_iterations)
        complex(dp), intent(in) :: a(:, :, :) !! input matrices (batch,n,n), batch>=1
        complex(dp), intent(out) :: inverse(:, :, :) !! pseudoinverses (batch,n,n), zero on failure
        integer, intent(out) :: status !! shared status; first batch failure aborts entire call
        integer, intent(out), optional :: rank(:) !! effective rank for each member (batch)
        real(dp), intent(in), optional :: atol !! shared absolute eigenvalue threshold
        real(dp), intent(in), optional :: rtol !! shared relative eigenvalue threshold
        logical, intent(in), optional :: lower !! use lower input triangle, default true
        logical, intent(in), optional :: check_finite !! check all input entries, default true
        integer, intent(in), optional :: max_iterations !! eigensolver iteration limit

        integer :: b, n, nb, member_rank

        inverse = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(rank)) rank = 0
        status = linalg_status_invalid_shape
        nb = size(a,1)
        n = size(a,2)
        if (nb < 1 .or. size(a,3) /= n) return
        if (size(inverse,1) /= nb .or. size(inverse,2) /= n .or. size(inverse,3) /= n) return
        if (present(rank)) then
            if (size(rank) /= nb) return
        end if
        do b = 1, nb
            call pinvh_complex(a(b,:,:), inverse(b,:,:), status, member_rank, &
                atol, rtol, lower, check_finite, max_iterations)
            if (status /= linalg_status_success) then
                inverse = cmplx(0.0_dp, 0.0_dp, dp)
                if (present(rank)) rank = 0
                return
            end if
            if (present(rank)) rank(b) = member_rank
        end do
        status = linalg_status_success
    end subroutine pinvh_complex_batch

    pure logical function valid_tolerances(n, atol, rtol, max_iterations) result(ok)
        integer, intent(in) :: n !! matrix order, nonnegative
        real(dp), intent(in), optional :: atol !! absolute tolerance, nonnegative and finite
        real(dp), intent(in), optional :: rtol !! relative tolerance, nonnegative and finite
        integer, intent(in), optional :: max_iterations !! eigensolver iterations, positive
        ok = n >= 0
        if (present(atol)) then
            if (.not. ieee_is_finite(atol)) ok = .false.
            if (atol < 0.0_dp) ok = .false.
        end if
        if (present(rtol)) then
            if (.not. ieee_is_finite(rtol)) ok = .false.
            if (rtol < 0.0_dp) ok = .false.
        end if
        if (present(max_iterations)) then
            if (max_iterations < 1) ok = .false.
        end if
    end function valid_tolerances

end module scifort_linalg_pinvh
