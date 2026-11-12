! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! General analytic functions of small dense matrices.  The algorithm is
! blocked Schur-Parlett with unitary reordering of clustered eigenvalues,
! Taylor diagonal blocks (optionally adaptive), and triangular Sylvester solves.  It is not a
! translation of scipy.linalg.funm.
module scifort_linalg_funm
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_matrix_functions, only : matrix_function_max_order, dense_complex_schur
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    real(dp), parameter :: default_separation_tol = 1.0e-10_dp
    real(dp), parameter :: derivative_block_tol = 4096.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: taylor_relative_tol = 256.0_dp*epsilon(1.0_dp)
    integer, parameter :: taylor_order_limit = 128
    integer, parameter :: required_small_terms = 3

    public :: funm, funm_real, scalar_matrix_function, scalar_matrix_derivative

    abstract interface
        function scalar_matrix_function(z) result(value)
            import dp
            complex(dp), intent(in) :: z !! complex scalar argument on the relevant analytic branch
            complex(dp) :: value
        end function scalar_matrix_function

        function scalar_matrix_derivative(z, order) result(value)
            import dp
            complex(dp), intent(in) :: z !! complex point at which a derivative is evaluated
            integer, intent(in) :: order !! k >= 1; beyond block order - 1 only if opted in
            complex(dp) :: value
        end function scalar_matrix_derivative
    end interface

    interface funm
        module procedure funm_real_input
        module procedure funm_complex_input
    end interface funm

contains

    subroutine funm_real_input(a, func, x, status, commutator_residual, &
            separation_tol, derivative, taylor_max_order)
        real(dp), intent(in) :: a(:, :) !! square real input, order 1..matrix_function_max_order
        procedure(scalar_matrix_function) :: func !! scalar analytic function evaluated at complex arguments
        complex(dp), intent(out) :: x(:, :) !! complex f(A); same square shape as input
        integer, intent(out) :: status !! linalg_status_*; nonzero means output is not valid
        real(dp), intent(out), optional :: commutator_residual !! normalized ||AF-FA||_1; not an error bound
        real(dp), intent(in), optional :: separation_tol !! relative eigenvalue clustering threshold, default 1e-10
        procedure(scalar_matrix_derivative), optional :: derivative !! k-th scalar derivative, needed for Jordan blocks
        integer, intent(in), optional :: taylor_max_order !! optional higher derivative limit, 1..128

        call funm_core(cmplx(a, kind=dp), func, x, status, &
            commutator_residual, separation_tol, derivative, taylor_max_order)
    end subroutine funm_real_input

    subroutine funm_complex_input(a, func, x, status, commutator_residual, &
            separation_tol, derivative, taylor_max_order)
        complex(dp), intent(in) :: a(:, :) !! square complex input, order 1..matrix_function_max_order
        procedure(scalar_matrix_function) :: func !! analytic scalar function defined around the spectrum
        complex(dp), intent(out) :: x(:, :) !! complex f(A), same square shape as input
        integer, intent(out) :: status !! linalg_status_*; nonzero means output is not valid
        real(dp), intent(out), optional :: commutator_residual !! normalized ||AF-FA||_1; not an error bound
        real(dp), intent(in), optional :: separation_tol !! relative eigenvalue clustering threshold, default 1e-10
        procedure(scalar_matrix_derivative), optional :: derivative !! k-th scalar derivative for repeated eigenvalues
        integer, intent(in), optional :: taylor_max_order !! optional higher derivative limit, 1..128

        call funm_core(a, func, x, status, commutator_residual, separation_tol, &
            derivative, taylor_max_order)
    end subroutine funm_complex_input

    subroutine funm_real(a, func, x, status, commutator_residual, separation_tol, &
            derivative, taylor_max_order)
        real(dp), intent(in) :: a(:, :) !! square real input, order 1..matrix_function_max_order
        procedure(scalar_matrix_function) :: func !! analytic function, real on input spectrum
        real(dp), intent(out) :: x(:, :) !! real result f(A) when its imaginary part is negligible
        integer, intent(out) :: status !! nonzero when f(A) is non-real or cannot be evaluated
        real(dp), intent(out), optional :: commutator_residual !! normalized commutator; not an accuracy bound
        real(dp), intent(in), optional :: separation_tol !! eigenvalue clustering threshold, default 1e-10
        procedure(scalar_matrix_derivative), optional :: derivative !! optional k-th scalar derivative
        integer, intent(in), optional :: taylor_max_order !! optional higher derivative limit, 1..128
        complex(dp), allocatable :: result_complex(:, :)
        real(dp) :: scale

        x = 0.0_dp
        allocate(result_complex(size(x, 1), size(x, 2)))
        call funm(a, func, result_complex, status, commutator_residual, &
            separation_tol, derivative, taylor_max_order)
        if (status /= linalg_status_success) return
        scale = max(1.0_dp, maxval(abs(result_complex)))
        if (maxval(abs(aimag(result_complex))) > 1.0e-10_dp*scale) then
            status = linalg_status_no_convergence
            return
        end if
        x = real(result_complex, dp)
    end subroutine funm_real

    subroutine funm_core(a, func, x, status, commutator_residual, &
            separation_tol, derivative, taylor_max_order)
        complex(dp), intent(in) :: a(:, :) !! original square matrix
        procedure(scalar_matrix_function) :: func !! scalar analytic function
        complex(dp), intent(out) :: x(:, :) !! resulting matrix or zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(out), optional :: commutator_residual !! normalized commutation defect, not error bound
        real(dp), intent(in), optional :: separation_tol !! relative spectral gap threshold
        procedure(scalar_matrix_derivative), optional :: derivative !! optional derivatives for repeated eigenvalues
        integer, intent(in), optional :: taylor_max_order !! optional higher derivative limit, 1..128
        complex(dp), allocatable :: t(:, :), z(:, :), f(:, :), block(:, :), rhs(:, :), result(:, :)
        complex(dp) :: value
        real(dp) :: magnitude, tolerance, defect, rounding_floor
        integer :: n, i, j, k, width, code, first, last, nblocks, left, right
        integer :: i1, i2, j1, j2, k1, k2
        integer, allocatable :: block_id(:), block_first(:), block_last(:)

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        if (present(commutator_residual)) commutator_residual = huge(1.0_dp)
        n = size(a, 1)
        if (n < 1 .or. n > matrix_function_max_order) return
        if (size(a, 2) /= n .or. size(x, 1) /= n .or. size(x, 2) /= n) return
        tolerance = default_separation_tol
        if (present(separation_tol)) tolerance = separation_tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_dp .or. tolerance >= 0.1_dp) return
        if (present(taylor_max_order)) then
            if (taylor_max_order < 1 .or. taylor_max_order > taylor_order_limit) return
        end if
        if (.not. finite_matrix(a)) return
        magnitude = maxval(abs(a))
        if (.not. ieee_is_finite(magnitude)) return
        if (magnitude == 0.0_dp) then
            value = func(cmplx(0.0_dp, 0.0_dp, dp))
            if (.not. finite_scalar(value)) then
                status = linalg_status_no_convergence
                return
            end if
            do i = 1, n
                x(i, i) = value
            end do
            if (present(commutator_residual)) commutator_residual = 0.0_dp
            status = linalg_status_success
            return
        end if
        allocate(t(n, n), z(n, n), f(n, n), block(n, n), block_id(n))
        ! Positive scaling leaves Schur vectors unchanged and protects QR from
        ! very large and very small entries.  T is A / magnitude.
        call dense_complex_schur(a/magnitude, t, z, code)
        if (code /= linalg_status_success) then
            status = linalg_status_no_convergence
            return
        end if
        ! Form *global* spectral clusters, even when the scalar derivative
        ! is unavailable.  Reordering can make repeated diagonalizable
        ! eigenvalues into diagonal scalar blocks without differentiation.
        call reorder_schur_clusters(t, z, tolerance, block_id, code)
        if (code /= linalg_status_success) then
            status = code
            return
        end if
        ! Detect excessive drift in the reordered Schur factorization.
        defect = matrix_onenorm(matmul(matmul(z, t), transpose(conjg(z)))-a/magnitude)
        if (.not. ieee_is_finite(defect) .or. &
                defect > 1.0e-8_dp*max(1.0_dp, matrix_onenorm(a/magnitude))) then
            status = linalg_status_no_convergence
            return
        end if

        allocate(block_first(n), block_last(n))
        nblocks = 0
        first = 1
        do while (first <= n)
            last = first
            do while (last < n)
                if (block_id(last+1) /= block_id(first)) exit
                last = last+1
            end do
            nblocks = nblocks+1
            block_first(nblocks) = first
            block_last(nblocks) = last
            first = last+1
        end do

        f = cmplx(0.0_dp, 0.0_dp, dp)
        rounding_floor = 128.0_dp*epsilon(1.0_dp)
        do i = 1, n
            f(i, i) = func(magnitude*t(i, i))
        end do
        if (.not. finite_matrix(f)) then
            status = linalg_status_no_convergence
            return
        end if

        ! Evaluate each diagonal Schur block independently.  Derivatives
        ! give a truncated Taylor evaluation for coupled clusters.  Without
        ! derivatives, a *diagonal* repeated-eigenvalue block needs only the
        ! scalar function; a coupled repeated-eigenvalue block is rejected.
        do i = 1, nblocks
            first = block_first(i)
            last = block_last(i)
            width = last-first+1
            if (width == 1) cycle
            if (present(derivative)) then
                call taylor_block(t(first:last, first:last), magnitude, func, &
                    derivative, tolerance, block(first:last, first:last), code, taylor_max_order)
                if (code /= linalg_status_success) then
                    status = code
                    return
                end if
                f(first:last, first:last) = block(first:last, first:last)
            else
                ! Test the upper off-diagonal explicitly.  Do not confuse
                ! a diagonalizable repeated spectrum with a Jordan block.
                ! The scale is in normalized Schur coordinates.
                do j = first+1, last
                    if (maxval(abs(t(first:j-1, j))) > rounding_floor* &
                            max(1.0_dp, matrix_onenorm(t(first:last, first:last)))) then
                        status = linalg_status_no_convergence
                        return
                    end if
                end do
            end if
        end do

        ! The block Parlett relation T*F = F*T yields, for i < j,
        !     Tii*Fij - Fij*Tjj = Fii*Tij - Tij*Fjj
        !                            + sum(Fik*Tkj - Tik*Fkj).
        ! Solve these *block* Sylvester equations in increasing block
        ! distance; this also works when different clusters have different
        ! sizes.  Each Sylvester solve uses triangular substitutions.
        do width = 1, nblocks-1
            do i = 1, nblocks-width
                j = i+width
                i1 = block_first(i)
                i2 = block_last(i)
                j1 = block_first(j)
                j2 = block_last(j)
                left = i2-i1+1
                right = j2-j1+1
                allocate(rhs(left, right), result(left, right))
                rhs = matmul(f(i1:i2, i1:i2), t(i1:i2, j1:j2)) - &
                    matmul(t(i1:i2, j1:j2), f(j1:j2, j1:j2))
                do k = i+1, j-1
                    k1 = block_first(k)
                    k2 = block_last(k)
                    rhs = rhs + matmul(f(i1:i2, k1:k2), t(k1:k2, j1:j2)) - &
                        matmul(t(i1:i2, k1:k2), f(k1:k2, j1:j2))
                end do
                call triangular_sylvester(t(i1:i2, i1:i2), t(j1:j2, j1:j2), &
                    rhs, tolerance, result, code)
                if (code /= linalg_status_success) then
                    status = code
                    return
                end if
                f(i1:i2, j1:j2) = result
                deallocate(rhs, result)
            end do
        end do
        ! The commutator check is necessary, but not sufficient, for an
        ! accurate matrix function.  Do not advertise it as an error bound.
        defect = matrix_onenorm(matmul(t, f)-matmul(f, t))/ &
            max(tiny(1.0_dp), matrix_onenorm(t)*matrix_onenorm(f))
        if (.not. ieee_is_finite(defect) .or. defect > max(1.0e-9_dp, 100.0_dp*tolerance)) then
            status = linalg_status_no_convergence
            return
        end if
        x = matmul(matmul(z, f), transpose(conjg(z)))
        if (.not. finite_matrix(x)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        if (present(commutator_residual)) commutator_residual = defect
        status = linalg_status_success
    end subroutine funm_core


    ! Solve A*X - X*B = C for upper triangular A and B.  The Schur
    ! reordering guarantees that different diagonal blocks have separated
    ! spectra; the threshold also guards against an unresolved divisor.
    subroutine triangular_sylvester(a, b, c, tolerance, x, status)
        complex(dp), intent(in) :: a(:, :) !! upper triangular left coefficient, (m,m)
        complex(dp), intent(in) :: b(:, :) !! upper triangular right coefficient, (n,n)
        complex(dp), intent(in) :: c(:, :) !! right-hand side, (m,n)
        real(dp), intent(in) :: tolerance !! minimum relative spectral separation
        complex(dp), intent(out) :: x(:, :) !! solution (m,n), zero on failure
        integer, intent(out) :: status !! success or numerically unresolved divisor
        complex(dp) :: rhs, delta
        real(dp) :: threshold
        integer :: i, j, k, m, n

        m = size(a,1)
        n = size(b,1)
        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_no_convergence
        do j = 1, n
            do i = m, 1, -1
                rhs = c(i,j)
                do k = i+1, m
                    rhs = rhs-a(i,k)*x(k,j)
                end do
                do k = 1, j-1
                    rhs = rhs+x(i,k)*b(k,j)
                end do
                delta = a(i,i)-b(j,j)
                threshold = tolerance*max(1.0_dp, abs(a(i,i)), abs(b(j,j)))
                if (abs(delta) <= threshold) return
                x(i,j) = rhs/delta
                if (.not. finite_scalar(x(i,j))) return
            end do
        end do
        status = linalg_status_success
    end subroutine triangular_sylvester

    ! Group eigenvalues by connected components under the user-supplied
    ! separation threshold, then stably reorder the Schur diagonal so that
    ! every group is contiguous.  This does not change eigenvalues or use
    ! nonunitary eigenvector transforms; it also handles diagonalizable
    ! repeated eigenvalues when no derivative callback is available.
    subroutine reorder_schur_clusters(t, z, tolerance, block_id, status)
        complex(dp), intent(inout) :: t(:, :) !! scaled upper triangular Schur matrix
        complex(dp), intent(inout) :: z(:, :) !! Schur basis transformed with t
        real(dp), intent(in) :: tolerance !! relative eigenvalue proximity threshold
        integer, intent(out) :: block_id(:) !! contiguous cluster number for each Schur index
        integer, intent(out) :: status !! success or a failed adjacent unitary swap
        integer :: i, j, k, left_id, right_id, n, temp_id

        n = size(t, 1)
        status = linalg_status_success
        do i = 1, n
            block_id(i) = i
        end do
        do i = 1, n-1
            do j = i+1, n
                if (abs(t(i,i)-t(j,j)) > &
                    tolerance*max(1.0_dp,abs(t(i,i)),abs(t(j,j)))) cycle
                left_id = block_id(i)
                right_id = block_id(j)
                where (block_id == left_id .or. block_id == right_id)
                    block_id = min(left_id, right_id)
                end where
            end do
        end do

        ! Insertion-sort the group identifiers by stable adjacent Schur swaps.
        ! A swap is only made across *different* groups, whose eigenvalues
        ! satisfy the requested separation threshold.
        do i = 2, n
            k = i
            do while (k > 1)
                if (block_id(k) >= block_id(k-1)) exit
                call swap_schur_neighbors(t, z, k-1, status)
                if (status /= linalg_status_success) return
                temp_id = block_id(k)
                block_id(k) = block_id(k-1)
                block_id(k-1) = temp_id
                k = k-1
            end do
        end do
    end subroutine reorder_schur_clusters

    ! Swap the diagonal entries of an upper triangular Schur matrix at p and
    ! p+1.  The normalized eigenvector of the lower-right eigenvalue gives
    ! the first column of a 2x2 complex unitary rotation, ensuring the
    ! transformed 2x2 block remains upper triangular (up to roundoff).
    subroutine swap_schur_neighbors(t, z, p, status)
        complex(dp), intent(inout) :: t(:, :) !! scaled triangular Schur factor
        complex(dp), intent(inout) :: z(:, :) !! unitary Schur vectors
        integer, intent(in) :: p !! adjacent swap index, in 1..n-1
        integer, intent(out) :: status !! success, or numerical swap failure
        complex(dp) :: q11, q12, q21, q22, a, b, v1, v2
        real(dp) :: scale, norm_pair, subdiag
        integer :: j, n

        status = linalg_status_no_convergence
        n = size(t, 1)
        b = t(p,p+1)
        a = t(p+1,p+1)-t(p,p)
        norm_pair = hypot(abs(b),abs(a))
        if (.not. ieee_is_finite(norm_pair) .or. norm_pair == 0.0_dp) return
        q11 = b/norm_pair
        q21 = a/norm_pair
        q12 = -conjg(q21)
        q22 = conjg(q11)

        ! Form Q^H T Q using only the rows/columns touched by the rotation.
        do j = p, n
            v1 = t(p,j)
            v2 = t(p+1,j)
            t(p,j) = conjg(q11)*v1+conjg(q21)*v2
            t(p+1,j) = conjg(q12)*v1+conjg(q22)*v2
        end do
        do j = 1, p+1
            v1 = t(j,p)
            v2 = t(j,p+1)
            t(j,p) = v1*q11+v2*q21
            t(j,p+1) = v1*q12+v2*q22
        end do
        do j = 1, n
            v1 = z(j,p)
            v2 = z(j,p+1)
            z(j,p) = v1*q11+v2*q21
            z(j,p+1) = v1*q12+v2*q22
        end do
        subdiag = abs(t(p+1,p))
        scale = max(1.0_dp, abs(t(p,p)), abs(t(p+1,p+1)), abs(t(p,p+1)))
        if (.not. ieee_is_finite(subdiag) .or. &
                subdiag > 1024.0_dp*epsilon(1.0_dp)*scale) return
        t(p+1,p) = cmplx(0.0_dp,0.0_dp,dp)
        if (.not. finite_matrix(t) .or. .not. finite_matrix(z)) return
        status = linalg_status_success
    end subroutine swap_schur_neighbors

    subroutine taylor_block(t, magnitude, func, derivative, tolerance, f, status, max_order)
        complex(dp), intent(in) :: t(:, :) !! contiguous Schur diagonal block in scaled units
        real(dp), intent(in) :: magnitude !! original matrix magnitude, strictly positive
        procedure(scalar_matrix_function) :: func !! scalar f(z)
        procedure(scalar_matrix_derivative) :: derivative !! f^(k)(z) for requested k
        real(dp), intent(in) :: tolerance !! relative eigenvalue clustering tolerance
        complex(dp), intent(out) :: f(:, :) !! result on the block
        integer, intent(out) :: status !! success only when the selected Taylor mode converges
        integer, intent(in), optional :: max_order !! higher order limit, default n-1
        complex(dp) :: center, value, exact_value
        complex(dp), allocatable :: nil(:, :), power(:, :), identity(:, :), remainder(:, :), term(:, :)
        real(dp) :: defect, bound, f_scale
        integer :: n, i, k, limit, small_terms
        logical :: extended, converged, exact_nilpotent

        n = size(t, 1)
        status = linalg_status_no_convergence
        f = cmplx(0.0_dp, 0.0_dp, dp)
        limit = n-1
        if (present(max_order)) limit = max_order
        if (limit < n-1 .or. limit > taylor_order_limit) then
            status = linalg_status_invalid_shape
            return
        end if
        extended = limit > n-1
        center = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            center = center+t(i, i)
        end do
        center = center/real(n, dp)
        allocate(nil(n,n), power(n,n), identity(n,n), remainder(n,n), term(n,n))
        identity = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            identity(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        nil = t
        do i = 1, n
            nil(i, i) = nil(i, i)-center
        end do
        ! A strictly upper triangular block deviation is nilpotent.
        ! Check equality of its *scaled* diagonal to zero rather than
        ! whether a high power underflows: a close but distinct spectrum
        ! could otherwise be mistaken for an exactly repeated eigenvalue.
        exact_nilpotent = .false.
        if (extended) then
            exact_nilpotent = .true.
            do i = 1, n
                if (nil(i,i) /= cmplx(0.0_dp, 0.0_dp, dp)) exact_nilpotent = .false.
            end do
        end if
        if (.not. extended) then
            ! Preserve the historical callback contract: derivatives are
            ! requested at most through block order n-1.  This is exact
            ! for nilpotent deviations and only a screened approximation
            ! for distinct but very close eigenvalues.
            remainder = identity
            do k = 1, n
                remainder = matmul(remainder, nil)
            end do
            defect = matrix_onenorm(remainder)
            bound = max(derivative_block_tol, 0.01_dp*tolerance)* &
                max(1.0_dp, matrix_onenorm(nil))**n
            if (.not. ieee_is_finite(defect) .or. defect > bound) return
        end if
        value = func(magnitude*center)
        if (.not. finite_scalar(value)) return
        f = value*identity
        power = identity
        converged = .not. extended
        small_terms = 0
        do k = 1, limit
            power = matmul(power, magnitude*nil)/real(k, dp)
            ! If the block deviation is exactly nilpotent, the polynomial
            ! already terminates.  Avoid requesting undefined higher
            ! derivatives from user callbacks in that case.
            if (extended .and. k >= n .and. exact_nilpotent) then
                converged = .true.
                exit
            end if
            value = derivative(magnitude*center, k)
            if (.not. finite_scalar(value)) return
            term = value*power
            f = f+term
            if (.not. finite_matrix(f)) return
            if (.not. extended .or. k < n) cycle
            ! Require three consecutive negligible tail terms, starting
            ! no earlier than degree n.  The stopping test is heuristic,
            ! not a proof of convergence for arbitrary analytic callbacks.
            f_scale = max(tiny(1.0_dp), matrix_onenorm(f))
            if (matrix_onenorm(term) <= taylor_relative_tol*f_scale) then
                small_terms = small_terms+1
            else
                small_terms = 0
            end if
            if (small_terms >= required_small_terms) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) return
        ! Check the scalar spectral mapping on every Schur diagonal.  The
        ! commutator alone cannot detect a truncated Taylor series.
        do i = 1, n
            exact_value = func(magnitude*t(i,i))
            if (.not. finite_scalar(exact_value)) return
            if (abs(f(i,i)-exact_value) > 1.0e-11_dp* &
                    max(1.0_dp, abs(exact_value), abs(f(i,i)))) return
        end do
        status = linalg_status_success
    end subroutine taylor_block

    pure real(dp) function matrix_onenorm(a) result(value)
        complex(dp), intent(in) :: a(:, :) !! dense matrix with finite elements
        value = maxval(sum(abs(a), dim=1))
    end function matrix_onenorm

    pure logical function finite_scalar(value)
        complex(dp), intent(in) :: value !! complex scalar
        finite_scalar = ieee_is_finite(real(value, dp)) .and. ieee_is_finite(aimag(value))
    end function finite_scalar

    pure logical function finite_matrix(a)
        complex(dp), intent(in) :: a(:, :) !! matrix with possibly nonfinite entries
        finite_matrix = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_matrix

end module scifort_linalg_funm
