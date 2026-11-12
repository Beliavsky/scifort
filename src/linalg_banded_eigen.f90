! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Eigenvalues of real symmetric tridiagonal matrices and real symmetric /
! complex Hermitian matrices in SciPy-compatible diagonal-ordered band storage.
! The tridiagonal solver uses implicit-shift QL; the general-band path first
! expands the band and reduces it by SciFort's Householder Hessenberg solver.
! Selected banded eigenvalues use O(n*kd) LDL/Sturm bisection.
! Banded eigenvectors use pivoted band LU inverse iteration with
! O(n*kd+n*k) temporary/output storage in selected or opt-in full modes.
! Full eigenvectors default to dense Householder reduction.
! Narrow-band, large, full eigenvalue-only requests can use Givens reduction
! automatically; method="bulge" remains available for explicit selection.
module scifort_linalg_banded_eigen
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_hessenberg, only : hessenberg
    use scifort_linalg_banded_tridiagonal, only : band_to_tridiagonal
    use scifort_linalg_banded, only : lu_factor_banded, lu_solve_banded
    implicit none
    private

    public :: eigh_tridiagonal, eigvalsh_tridiagonal, eig_banded, eigvals_banded

    interface eig_banded
        module procedure eig_banded_real
        module procedure eig_banded_complex
    end interface eig_banded

    interface eigvals_banded
        module procedure eigvals_banded_real
        module procedure eigvals_banded_complex
    end interface eigvals_banded

contains

    pure logical function auto_bulge_full_values(n, kd) result(use_bulge)
        integer, intent(in) :: n !! matrix order
        integer, intent(in) :: kd !! number of stored offdiagonals

        ! Conservative deterministic crossover for eigenvalues-only requests.
        ! Keep full eigenvector requests on the original dense default.
        use_bulge = n >= 192 .and. kd >= 0 .and. kd <= min(8, max(1,n/16))
    end function auto_bulge_full_values

    subroutine eigh_tridiagonal(d, e, w, status, z, max_iterations, &
        subset_by_index, subset_by_value, num_selected, tol, eigenvector_method)
        real(dp), intent(in) :: d(:) !! finite diagonal of symmetric tridiagonal matrix, order n
        real(dp), intent(in) :: e(:) !! finite offdiagonal, length max(0,n-1)
        real(dp), intent(out) :: w(:) !! ascending eigenvalues, capacity in selection modes
        integer, intent(out) :: status !! success, invalid shape or nonconvergence
        real(dp), intent(out), optional :: z(:, :) !! selected eigenvectors, shape (n,capacity)
        integer, intent(in), optional :: max_iterations !! QL steps or bisection steps, default 120
        integer, intent(in), optional :: subset_by_index(2) !! inclusive one-based [first,last]
        real(dp), intent(in), optional :: subset_by_value(2) !! select (low,high], infinities allowed
        integer, intent(out), optional :: num_selected !! count, even if outputs have insufficient space
        real(dp), intent(in), optional :: tol !! absolute eigenvalue tolerance for bisection, >=0
        character(len=*), intent(in), optional :: eigenvector_method !! auto, inverse, or ql

        real(dp), allocatable :: allw(:), allz(:, :)
        integer :: n, first, last, m, i, limit
        character(len=:), allocatable :: vec_method
        logical :: selected

        w = 0.0_dp
        if (present(z)) z = 0.0_dp
        if (present(num_selected)) num_selected = 0
        status = linalg_status_invalid_shape
        n = size(d)
        if (size(e) /= max(0, n-1)) return
        if (present(z)) then
            if (size(z,1) /= n) return
        end if
        if (present(subset_by_index) .and. present(subset_by_value)) return
        if (present(subset_by_index)) then
            if (n < 1) return
            if (subset_by_index(1) < 1 .or. subset_by_index(2) > n) return
            if (subset_by_index(1) > subset_by_index(2)) return
        end if
        if (present(subset_by_value)) then
            if (any(ieee_is_nan(subset_by_value))) return
            if (subset_by_value(1) >= subset_by_value(2)) return
        end if
        limit = 120
        if (present(max_iterations)) limit = max_iterations
        if (limit < 1) return
        if (present(tol)) then
            if (.not. ieee_is_finite(tol)) return
            if (tol < 0.0_dp) return
        end if
        vec_method = 'auto'
        if (present(eigenvector_method)) vec_method = trim(adjustl(eigenvector_method))
        if (vec_method /= 'auto' .and. vec_method /= 'inverse' .and. vec_method /= 'ql') return
        selected = present(subset_by_index) .or. present(subset_by_value)
        if (.not. selected) then
            if (size(w) /= n) return
            if (present(z)) then
                if (size(z,2) /= n) return
            end if
            call eigh_tridiagonal_all(d, e, w, status, z, max_iterations)
            if (status == linalg_status_success) then
                if (present(num_selected)) num_selected = n
            end if
            return
        end if
        if (.not. all(ieee_is_finite(d)) .or. .not. all(ieee_is_finite(e))) then
            status = linalg_status_no_convergence
            return
        end if

        ! Selected eigenvectors use Sturm bisection plus inverse iteration.
        ! They require O(n*m) output storage instead of the full O(n**2)
        ! QL eigenvector matrix (m is the number of selected eigenvectors).
        if (.not. present(z) .or. vec_method /= 'ql') then
            call selected_tridiagonal_bisection(d, e, w, status, limit, &
                subset_by_index, subset_by_value, m, tol)
            if (present(num_selected)) num_selected = m
            if (status /= linalg_status_success) return
            if (.not. present(z)) return
            if (size(z,2) < m) then
                w = 0.0_dp
                status = linalg_status_invalid_shape
                return
            end if
            if (m == 0) return
            call selected_tridiagonal_vectors(d, e, w(1:m), z(:,1:m), status)
            if (status == linalg_status_success) return
            ! Auto mode keeps the existing full-QL path as a safety net for
            ! small, exceptionally difficult clustered spectra.
            if (vec_method == 'inverse' .or. n > 256) then
                w = 0.0_dp
                z = 0.0_dp
                return
            end if
        end if
        allocate(allw(n), allz(n,n))
        call eigh_tridiagonal_all(d, e, allw, status, allz, max_iterations)
        if (status /= linalg_status_success) return
        first = 1
        last = n
        if (present(subset_by_index)) then
            first = subset_by_index(1)
            last = subset_by_index(2)
        end if
        if (present(subset_by_value)) then
            first = n+1
            last = n
            do i = 1, n
                if (allw(i) > subset_by_value(1)) then
                    first = i
                    exit
                end if
            end do
            last = first-1
            do i = first, n
                if (allw(i) > subset_by_value(2)) exit
                last = i
            end do
        end if
        m = last-first+1
        if (present(num_selected)) num_selected = m
        if (size(w) < m .or. size(z,2) < m) then
            status = linalg_status_invalid_shape
            return
        end if
        if (m > 0) then
            w(1:m) = allw(first:last)
            z(:,1:m) = allz(:,first:last)
        end if
    end subroutine eigh_tridiagonal

    subroutine selected_tridiagonal_vectors(d, e, w, z, status)
        ! Inverse iteration after Sturm bisection, with two-pass reorthogonalization.
        ! The factorization is an independent tridiagonal partial-pivot solve.
        ! It stores O(n) work arrays; the output matrix has n*m elements.
        real(dp), intent(in) :: d(:) !! tridiagonal diagonal, length n
        real(dp), intent(in) :: e(:) !! offdiagonal, length n-1
        real(dp), intent(in) :: w(:) !! selected eigenvalues in ascending order
        real(dp), intent(out) :: z(:, :) !! normalized eigenvectors, shape (n,m)
        integer, intent(out) :: status !! success or inverse-iteration failure

        real(dp), allocatable :: dn(:), en(:), diag(:), sub(:), sup(:), sup2(:)
        real(dp), allocatable :: rhs(:), x(:), residual(:)
        integer, allocatable :: piv(:), ordering(:), workspace(:)
        real(dp) :: sigma, shift, shifted, perturb, normx, res, threshold, proj
        integer :: n, m, j, i, k, it, attempt, pass, candidate_index
        logical :: converged

        n = size(d)
        m = size(w)
        z = 0.0_dp
        status = linalg_status_no_convergence
        if (m == 0) then
            status = linalg_status_success
            return
        end if

        ! A diagonal input has exact coordinate eigenvectors. Stable sorting
        ! also produces independent vectors for eigenvalues of multiplicity > 1.
        if (all(e == 0.0_dp)) then
            allocate(ordering(n), workspace(n))
            call stable_diagonal_order(d, ordering, workspace)
            candidate_index = 1
            do j = 1, m
                do while (candidate_index <= n)
                    if (d(ordering(candidate_index)) == w(j)) exit
                    candidate_index = candidate_index+1
                end do
                if (candidate_index > n) return
                z(ordering(candidate_index),j) = 1.0_dp
                candidate_index = candidate_index+1
            end do
            status = linalg_status_success
            return
        end if

        sigma = max(maxval(abs(d)),maxval(abs(e)))
        if (.not. ieee_is_finite(sigma) .or. sigma <= 0.0_dp) return
        allocate(dn(n),en(n-1),diag(n),sub(n-1),sup(n-1),sup2(max(0,n-2)))
        allocate(rhs(n),x(n),residual(n),piv(n-1))
        dn = d/sigma
        en = e/sigma
        threshold = max(3.0e-12_dp,256.0_dp*epsilon(1.0_dp)*sqrt(real(n,dp)))
        perturb = 16.0_dp*epsilon(1.0_dp)
        do j = 1, m
            shift = w(j)/sigma
            if (.not. ieee_is_finite(shift)) return
            converged = .false.
            do attempt = 1, 4
                ! Perturb the shift to avoid singular factorizations. The
                ! inverse-iteration residual is checked against the *unshifted*
                ! selected eigenvalue, not against the perturbed one.
                shifted = shift + (-1.0_dp)**attempt * perturb*real(attempt,dp)
                call shifted_tridiagonal_factor(dn, en, shifted, diag, sub, sup, sup2, piv)
                do i = 1, n
                    rhs(i) = sin(real(i,dp)*(real(j+3*attempt,dp)*0.681313_dp)) &
                        + cos(real(i,dp)*(real(j+attempt,dp)*0.318309_dp))
                end do
                do it = 1, 16
                    call shifted_tridiagonal_solve(diag, sub, sup, sup2, piv, rhs, x)
                    if (.not. all(ieee_is_finite(x))) exit
                    ! MGS twice is essential for clustered or equal eigenvalues.
                    do pass = 1, 2
                        do k = 1, j-1
                            proj = dot_product(z(:,k),x)
                            x = x-proj*z(:,k)
                        end do
                    end do
                    normx = maxval(abs(x))
                    if (normx < 1.0e-200_dp) exit
                    x = x/normx
                    normx = sqrt(sum(x*x))
                    if (normx <= tiny(1.0_dp)) exit
                    x = x/normx
                    residual = (dn-shift)*x
                    residual(1:n-1) = residual(1:n-1)+en*x(2:n)
                    residual(2:n) = residual(2:n)+en*x(1:n-1)
                    res = sqrt(sum(residual*residual))
                    if (res <= threshold) then
                        z(:,j) = x
                        converged = .true.
                        exit
                    end if
                    rhs = x
                end do
                if (converged) exit
            end do
            if (.not. converged) then
                z = 0.0_dp
                return
            end if
        end do
        status = linalg_status_success
    end subroutine selected_tridiagonal_vectors

    subroutine stable_diagonal_order(d, order, work)
        real(dp), intent(in) :: d(:) !! diagonal values to sort
        integer, intent(out) :: order(:) !! positions in nondecreasing stable order
        integer, intent(out) :: work(:) !! scratch space, same length as order
        integer :: width, lo, mid, hi, left, right, out, n, i

        n = size(d)
        do i = 1, n
            order(i) = i
        end do
        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n)
                mid = min(n,lo+width-1)
                hi = min(n,mid+width)
                left = lo
                right = mid+1
                do out = lo, hi
                    if (left > mid) then
                        work(out) = order(right)
                        right = right+1
                    else if (right > hi) then
                        work(out) = order(left)
                        left = left+1
                    else if (d(order(left)) <= d(order(right))) then
                        work(out) = order(left)
                        left = left+1
                    else
                        work(out) = order(right)
                        right = right+1
                    end if
                end do
                lo = hi+1
            end do
            order = work
            if (width > n/2) exit
            width = 2*width
        end do
    end subroutine stable_diagonal_order

    subroutine shifted_tridiagonal_factor(d, e, shift, diag, sub, sup, sup2, piv)
        ! Partial-pivoted LU of T-shift*I, retaining two superdiagonals.
        real(dp), intent(in) :: d(:) !! scaled diagonal, length n
        real(dp), intent(in) :: e(:) !! scaled offdiagonal, length n-1
        real(dp), intent(in) :: shift !! regularized shift
        real(dp), intent(out) :: diag(:) !! pivots of U, length n
        real(dp), intent(out) :: sub(:) !! multipliers of L, length n-1
        real(dp), intent(out) :: sup(:) !! first superdiagonal of U, length n-1
        real(dp), intent(out) :: sup2(:) !! second superdiagonal of U, length n-2
        integer, intent(out) :: piv(:) !! row pivot indices, length n-1

        real(dp) :: factor, swap, floor_pivot
        integer :: n, i

        n = size(d)
        floor_pivot = 16.0_dp*epsilon(1.0_dp)
        diag = d-shift
        sub = e
        sup = e
        sup2 = 0.0_dp
        do i = 1, n-1
            piv(i) = i
            if (max(abs(diag(i)),abs(sub(i))) < floor_pivot) then
                diag(i) = sign(floor_pivot,diag(i))
            end if
            if (abs(diag(i)) >= abs(sub(i))) then
                factor = sub(i)/diag(i)
                sub(i) = factor
                diag(i+1) = diag(i+1)-factor*sup(i)
            else
                factor = diag(i)/sub(i)
                diag(i) = sub(i)
                sub(i) = factor
                swap = sup(i)
                sup(i) = diag(i+1)
                diag(i+1) = swap-factor*diag(i+1)
                if (i < n-1) then
                    sup2(i) = sup(i+1)
                    sup(i+1) = -factor*sup(i+1)
                end if
                piv(i) = i+1
            end if
        end do
        if (abs(diag(n)) < floor_pivot) diag(n) = sign(floor_pivot,diag(n))
    end subroutine shifted_tridiagonal_factor

    subroutine shifted_tridiagonal_solve(diag, sub, sup, sup2, piv, rhs, x)
        real(dp), intent(in) :: diag(:) !! nonzero diagonal of U
        real(dp), intent(in) :: sub(:) !! multipliers of L
        real(dp), intent(in) :: sup(:) !! first superdiagonal of U
        real(dp), intent(in) :: sup2(:) !! second superdiagonal of U
        integer, intent(in) :: piv(:) !! pivot positions i or i+1
        real(dp), intent(in) :: rhs(:) !! solve right-hand side
        real(dp), intent(out) :: x(:) !! approximate inverse-iteration vector

        integer :: n, i
        real(dp) :: temp

        n = size(diag)
        x = rhs
        do i = 1, n-1
            if (piv(i) == i) then
                x(i+1) = x(i+1)-sub(i)*x(i)
            else
                temp = x(i)
                x(i) = x(i+1)
                x(i+1) = temp-sub(i)*x(i)
            end if
        end do
        x(n) = x(n)/diag(n)
        if (n > 1) x(n-1) = (x(n-1)-sup(n-1)*x(n))/diag(n-1)
        do i = n-2, 1, -1
            x(i) = (x(i)-sup(i)*x(i+1)-sup2(i)*x(i+2))/diag(i)
        end do
    end subroutine shifted_tridiagonal_solve

    subroutine selected_tridiagonal_bisection(d, e, w, status, max_steps, &
        subset_by_index, subset_by_value, num_selected, tol)
        real(dp), intent(in) :: d(:) !! finite tridiagonal diagonal, length n
        real(dp), intent(in) :: e(:) !! finite tridiagonal offdiagonal, length n-1
        real(dp), intent(out) :: w(:) !! output capacity for selected ascending eigenvalues
        integer, intent(out) :: status !! success, shape error, or failure to converge
        integer, intent(in) :: max_steps !! upper bound on bisection iterations, >=1
        integer, intent(in), optional :: subset_by_index(2) !! inclusive one-based range
        real(dp), intent(in), optional :: subset_by_value(2) !! selected eigenvalues in (low,high]
        integer, intent(out), optional :: num_selected !! number of selected eigenvalues
        real(dp), intent(in), optional :: tol !! nonnegative absolute eigenvalue tolerance

        real(dp), allocatable :: dn(:), en(:), sorted(:)
        real(dp) :: sigma, lower, upper, radius, left, right, middle, width_tol
        real(dp) :: margin, eigval, pmin
        integer :: n, i, j, k, first, last, steps
        logical :: converged

        w = 0.0_dp
        if (present(num_selected)) num_selected = 0
        status = linalg_status_no_convergence
        n = size(d)
        ! An exactly diagonal matrix has exact, already representable
        ! eigenvalues, even when diagonal magnitudes differ enormously.
        if (all(e == 0.0_dp)) then
            allocate(sorted(n))
            sorted = d
            call sort_ascending(sorted)
            first = 1
            last = n
            if (present(subset_by_index)) then
                first = subset_by_index(1)
                last = subset_by_index(2)
            end if
            if (present(subset_by_value)) then
                first = n+1
                do i = 1, n
                    if (sorted(i) > subset_by_value(1)) then
                        first = i
                        exit
                    end if
                end do
                last = first-1
                do i = first, n
                    if (sorted(i) > subset_by_value(2)) exit
                    last = i
                end do
            end if
            if (present(num_selected)) num_selected = last-first+1
            status = linalg_status_invalid_shape
            if (size(w) < last-first+1) return
            if (last >= first) w(1:last-first+1) = sorted(first:last)
            status = linalg_status_success
            return
        end if
        sigma = 0.0_dp
        if (n > 0) sigma = maxval(abs(d))
        if (n > 1) sigma = max(sigma, maxval(abs(e)))
        allocate(dn(n), en(max(0,n-1)))
        dn = d/sigma
        en = e/sigma
        lower = huge(1.0_dp)
        upper = -huge(1.0_dp)
        do i = 1, n
            radius = 0.0_dp
            if (i > 1) radius = radius + abs(en(i-1))
            if (i < n) radius = radius + abs(en(i))
            lower = min(lower, dn(i)-radius)
            upper = max(upper, dn(i)+radius)
        end do
        margin = 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lower),abs(upper))
        lower = lower-margin
        upper = upper+margin
        pmin = 4.0_dp*tiny(1.0_dp)
        first = 1
        last = n
        if (present(subset_by_index)) then
            first = subset_by_index(1)
            last = subset_by_index(2)
        end if
        if (present(subset_by_value)) then
            first = sturm_count(dn,en,bound_scaled(subset_by_value(1),sigma),pmin)+1
            last = sturm_count(dn,en,bound_scaled(subset_by_value(2),sigma),pmin)
        end if
        if (present(num_selected)) num_selected = last-first+1
        status = linalg_status_invalid_shape
        if (size(w) < last-first+1) return
        width_tol = 8.0_dp*epsilon(1.0_dp)
        if (present(tol)) then
            if (tol > 0.0_dp) width_tol = max(width_tol,abs(bound_scaled(tol,sigma)))
        end if
        status = linalg_status_no_convergence
        do j = 1, last-first+1
            k = first+j-1
            left = lower
            right = upper
            converged = .false.
            do steps = 1, max_steps
                middle = left+0.5_dp*(right-left)
                if (middle <= left .or. middle >= right) then
                    converged = .true.
                    exit
                end if
                if (sturm_count(dn,en,middle,pmin) < k) then
                    left = middle
                else
                    right = middle
                end if
                if (right-left <= width_tol) then
                    converged = .true.
                    exit
                end if
            end do
            if (.not. converged) then
                w = 0.0_dp
                return
            end if
            eigval = left+0.5_dp*(right-left)
            w(j) = eigval*sigma
        end do
        if (.not. all(ieee_is_finite(w))) then
            status = linalg_status_no_convergence
            w = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine selected_tridiagonal_bisection

    subroutine sort_ascending(x)
        real(dp), intent(inout) :: x(:) !! finite array sorted in ascending order
        real(dp), allocatable :: tmp(:)
        integer :: n, width, lo, mid, hi, left, right, out

        n = size(x)
        if (n <= 1) return
        allocate(tmp(n))
        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n)
                mid = min(n, lo+width-1)
                hi = min(n, mid+width)
                left = lo
                right = mid+1
                do out = lo, hi
                    if (left > mid) then
                        tmp(out) = x(right)
                        right = right+1
                    else if (right > hi) then
                        tmp(out) = x(left)
                        left = left+1
                    else if (x(left) <= x(right)) then
                        tmp(out) = x(left)
                        left = left+1
                    else
                        tmp(out) = x(right)
                        right = right+1
                    end if
                end do
                lo = hi+1
            end do
            x = tmp
            if (width > n/2) exit
            width = 2*width
        end do
    end subroutine sort_ascending

    pure integer function sturm_count(d, e, x, pivmin) result(count)
        real(dp), intent(in) :: d(:) !! normalized tridiagonal diagonal
        real(dp), intent(in) :: e(:) !! normalized offdiagonal
        real(dp), intent(in) :: x !! normalized spectral boundary
        real(dp), intent(in) :: pivmin !! positive minimum magnitude of Sturm pivots
        real(dp) :: pivot, term
        integer :: i

        count = 0
        if (size(d) == 0) return
        pivot = d(1)-x
        if (pivot <= 0.0_dp) count = count+1
        if (abs(pivot) <= pivmin) then
            if (pivot > 0.0_dp) count = count+1
            pivot = -pivmin
        end if
        do i = 2, size(d)
            ! The normalized offdiagonal is at most 1. Avoid tiny squares.
            term = 0.0_dp
            if (abs(e(i-1)) >= sqrt(tiny(1.0_dp))) term = e(i-1)*e(i-1)/pivot
            pivot = (d(i)-x)-term
            if (pivot <= 0.0_dp) count = count+1
            if (abs(pivot) <= pivmin) then
                if (pivot > 0.0_dp) count = count+1
                pivot = -pivmin
            end if
        end do
    end function sturm_count

    pure real(dp) function bound_scaled(bound, sigma) result(x)
        real(dp), intent(in) :: bound !! physical eigenvalue interval endpoint
        real(dp), intent(in) :: sigma !! strictly positive input scaling

        ! Clamp to [-4,4] before division to avoid overflow on wide scales.
        if (.not. ieee_is_finite(bound)) then
            x = sign(4.0_dp,bound)
        else if (abs(bound)/4.0_dp >= sigma) then
            x = sign(4.0_dp,bound)
        else
            x = bound/sigma
        end if
    end function bound_scaled

    subroutine eigh_tridiagonal_all(d, e, w, status, z, max_iterations)
        real(dp), intent(in) :: d(:) !! finite diagonal of a symmetric tridiagonal matrix, length n
        real(dp), intent(in) :: e(:) !! finite first subdiagonal, length max(0,n-1)
        real(dp), intent(out) :: w(:) !! eigenvalues sorted ascending, length n
        integer, intent(out) :: status !! zero on success, shape or iteration error otherwise
        real(dp), intent(out), optional :: z(:, :) !! orthonormal eigenvectors, shape (n,n)
        integer, intent(in), optional :: max_iterations !! QL iterations per eigenvalue, default 120

        real(dp), allocatable :: dd(:), ee(:), column(:)
        real(dp) :: amax, b, c, f, g, p, r, s, t, gap
        integer :: i, j, l, m, n, it, limit, row

        w = 0.0_dp
        if (present(z)) z = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(d)
        if (size(e) /= max(0, n - 1) .or. size(w) /= n) return
        if (present(z)) then
            if (size(z, 1) /= n .or. size(z, 2) /= n) return
        end if
        limit = 120
        if (present(max_iterations)) limit = max_iterations
        if (limit < 1) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(d)) .or. .not. all(ieee_is_finite(e))) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        if (present(z)) then
            do i = 1, n
                z(i, i) = 1.0_dp
            end do
        end if
        amax = maxval(abs(d))
        if (size(e) > 0) amax = max(amax, maxval(abs(e)))
        if (amax == 0.0_dp) then
            status = linalg_status_success
            return
        end if
        allocate(dd(n), ee(n))
        dd = d / amax
        ee = 0.0_dp
        ee(1:n-1) = e / amax

        ! Implicit-shift QL with Wilkinson shifts; e(n) is a sentinel.
        ! Scaling before iteration avoids overflowing squared shifts.
        do l = 1, n
            it = 0
            do
                m = l
                do while (m < n)
                    gap = epsilon(1.0_dp) * (abs(dd(m)) + abs(dd(m+1))) + tiny(1.0_dp)
                    if (abs(ee(m)) <= gap) exit
                    m = m + 1
                end do
                if (m == l) exit
                it = it + 1
                if (it > limit) return
                g = (dd(l+1) - dd(l)) / (2.0_dp * ee(l))
                r = hypot(g, 1.0_dp)
                g = dd(m) - dd(l) + ee(l) / (g + sign(r, g))
                s = 1.0_dp
                c = 1.0_dp
                p = 0.0_dp
                do i = m - 1, l, -1
                    f = s * ee(i)
                    b = c * ee(i)
                    if (f == 0.0_dp .and. g == 0.0_dp) then
                        ee(i+1) = 0.0_dp
                        s = 0.0_dp
                        c = 1.0_dp
                    else if (abs(f) >= abs(g)) then
                        c = g / f
                        r = hypot(c, 1.0_dp)
                        ee(i+1) = f * r
                        s = 1.0_dp / r
                        c = c * s
                    else
                        s = f / g
                        r = hypot(s, 1.0_dp)
                        ee(i+1) = g * r
                        c = 1.0_dp / r
                        s = s * c
                    end if
                    g = dd(i+1) - p
                    r = (dd(i) - g) * s + 2.0_dp * c * b
                    p = s * r
                    dd(i+1) = g + p
                    g = c * r - b
                    if (present(z)) then
                        do row = 1, n
                            t = z(row, i+1)
                            z(row, i+1) = s * z(row, i) + c * t
                            z(row, i) = c * z(row, i) - s * t
                        end do
                    end if
                end do
                dd(l) = dd(l) - p
                ee(l) = g
                ee(m) = 0.0_dp
            end do
        end do

        if (.not. all(ieee_is_finite(dd))) return
        w = dd * amax
        if (.not. all(ieee_is_finite(w))) return
        ! Insertion sort moves the corresponding eigenvector columns.
        if (present(z)) allocate(column(n))
        do i = 2, n
            t = w(i)
            if (present(z)) column = z(:, i)
            j = i - 1
            do while (j >= 1)
                if (w(j) <= t) exit
                w(j+1) = w(j)
                if (present(z)) z(:, j+1) = z(:, j)
                j = j - 1
            end do
            w(j+1) = t
            if (present(z)) z(:, j+1) = column
        end do
        status = linalg_status_success
    end subroutine eigh_tridiagonal_all

    subroutine eigvalsh_tridiagonal(d, e, w, status, max_iterations, &
        subset_by_index, subset_by_value, num_selected, tol)
        real(dp), intent(in) :: d(:) !! finite diagonal of tridiagonal matrix, order n
        real(dp), intent(in) :: e(:) !! finite offdiagonal, length max(0,n-1)
        real(dp), intent(out) :: w(:) !! ascending eigenvalues; capacity for value selection
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_iterations !! QL or bisection steps, default 120
        integer, intent(in), optional :: subset_by_index(2) !! inclusive, one-based index limits
        real(dp), intent(in), optional :: subset_by_value(2) !! lower-exclusive, upper-inclusive bounds
        integer, intent(out), optional :: num_selected !! number selected, even if output is too small
        real(dp), intent(in), optional :: tol !! nonnegative absolute bisection tolerance, default automatic

        call eigh_tridiagonal(d, e, w, status, max_iterations=max_iterations, &
            subset_by_index=subset_by_index, subset_by_value=subset_by_value, &
            num_selected=num_selected, tol=tol)
    end subroutine eigvalsh_tridiagonal

    subroutine eig_banded_real(ab, w, status, z, lower, max_iterations, &
        subset_by_index, subset_by_value, num_selected, tol, method)
        real(dp), intent(in) :: ab(:, :) !! SciPy-packed symmetric bands, (kd+1,n)
        real(dp), intent(out) :: w(:) !! selected ascending eigenvalues, capacity >= count
        integer, intent(out) :: status !! success, invalid shape, or failed convergence
        real(dp), intent(out), optional :: z(:, :) !! eigenvectors (n,capacity)
        logical, intent(in), optional :: lower !! lower-band storage, default false
        integer, intent(in), optional :: max_iterations !! QL or bisection limit
        integer, intent(in), optional :: subset_by_index(2) !! one-based inclusive indices
        real(dp), intent(in), optional :: subset_by_value(2) !! half-open interval (low,high]
        integer, intent(out), optional :: num_selected !! number of selected eigenpairs
        real(dp), intent(in), optional :: tol !! absolute bisection tolerance, >=0
        character(len=*), intent(in), optional :: method !! auto, inverse, bulge, or dense

        complex(dp), allocatable :: bc(:, :)
        real(dp), allocatable :: values(:), vectors(:, :)
        character(len=:), allocatable :: strategy
        integer :: n, m, first, last, j
        logical :: selected, full_inverse

        status = linalg_status_invalid_shape
        w = 0.0_dp
        if (present(z)) z = 0.0_dp
        if (present(num_selected)) num_selected = 0
        n = size(ab,2)
        if (size(ab,1) < 1) return
        if (present(z)) then
            if (size(z,1) /= n) return
        end if
        if (present(subset_by_index) .and. present(subset_by_value)) return
        strategy = 'auto'
        if (present(method)) strategy = trim(adjustl(method))
        if (strategy /= 'auto' .and. strategy /= 'inverse' .and. &
            strategy /= 'dense' .and. strategy /= 'bulge') return
        selected = present(subset_by_index) .or. present(subset_by_value)
        full_inverse = strategy == 'inverse' .and. .not. selected .and. n > 0
        if ((selected .or. full_inverse) .and. strategy /= 'dense' .and. &
            strategy /= 'bulge') then
            if (present(z)) then
                if (full_inverse) then
                    call selected_banded_real_vectors(ab,w,status,z,lower=lower, &
                        max_iterations=max_iterations,subset_by_index=[1,n], &
                        num_selected=num_selected,tol=tol)
                else
                    call selected_banded_real_vectors(ab,w,status,z,lower,max_iterations, &
                        subset_by_index,subset_by_value,num_selected,tol)
                end if
            else
                allocate(bc(size(ab,1),n))
                bc=cmplx(ab,0.0_dp,dp)
                if (full_inverse) then
                    call selected_banded_sturm(bc,w,status,lower=lower, &
                        max_iterations=max_iterations,subset_by_index=[1,n], &
                        num_selected=num_selected,tol=tol)
                else
                    call selected_banded_sturm(bc,w,status,lower,max_iterations, &
                        subset_by_index,subset_by_value,num_selected,tol)
                end if
            end if
            if (full_inverse) return
            if (status /= linalg_status_no_convergence .or. &
                strategy == 'inverse' .or. n > 96) return
            ! Auto: the original dense method is a fallback on small hard cases.
            w = 0.0_dp
            if (present(z)) z = 0.0_dp
            if (present(num_selected)) num_selected = 0
        end if
        if (.not. selected) then
            if (strategy == 'bulge' .or. &
                (strategy == 'auto' .and. .not. present(z) .and. &
                auto_bulge_full_values(n,size(ab,1)-1))) then
                call eig_banded_real_bulge(ab,w,status,z,lower,max_iterations)
                if (strategy == 'auto' .and. status == linalg_status_no_convergence) then
                    call eig_banded_real_dense(ab,w,status,z,lower,max_iterations)
                end if
            else
                call eig_banded_real_dense(ab,w,status,z,lower,max_iterations)
            end if
            if (status == linalg_status_success .and. present(num_selected)) num_selected=n
            return
        end if
        allocate(values(n))
        if (present(z)) then
            allocate(vectors(n,n))
            if (strategy == 'bulge') then
                call eig_banded_real_bulge(ab,values,status,vectors,lower,max_iterations)
            else
                call eig_banded_real_dense(ab,values,status,vectors,lower,max_iterations)
            end if
        else
            if (strategy == 'bulge') then
                call eig_banded_real_bulge(ab,values,status,lower=lower, &
                    max_iterations=max_iterations)
            else
                call eig_banded_real_dense(ab,values,status,lower=lower, &
                    max_iterations=max_iterations)
            end if
        end if
        if (status /= linalg_status_success) return
        call band_selection_indices(values,first,last,status,subset_by_index,subset_by_value)
        if (status /= linalg_status_success) return
        m = last-first+1
        if (present(num_selected)) num_selected=m
        status = linalg_status_invalid_shape
        if (size(w) < m) return
        if (present(z)) then
            if (size(z,2) < m) return
        end if
        if (m > 0) then
            w(1:m)=values(first:last)
            if (present(z)) z(:,1:m)=vectors(:,first:last)
        end if
        status = linalg_status_success
    end subroutine eig_banded_real

    subroutine selected_banded_real_vectors(ab,w,status,z,lower,max_iterations, &
        subset_by_index,subset_by_value,num_selected,tol)
        ! Real inverse iteration uses the caller's real output directly;
        ! it never builds a complex n-by-k eigenvector scratch array.
        real(dp), intent(in) :: ab(:, :) !! SciPy-packed real symmetric bands
        real(dp), intent(out) :: w(:) !! selected eigenvalues in ascending order
        integer, intent(out) :: status !! success, invalid shape, or no convergence
        real(dp), intent(out) :: z(:, :) !! eigenvectors, shape (n,capacity)
        logical, intent(in), optional :: lower !! lower-packed if true
        integer, intent(in), optional :: max_iterations !! Sturm bisection limit
        integer, intent(in), optional :: subset_by_index(2) !! inclusive one-based limits
        real(dp), intent(in), optional :: subset_by_value(2) !! half-open value interval
        integer, intent(out), optional :: num_selected !! count including insufficient storage
        real(dp), intent(in), optional :: tol !! absolute eigenvalue tolerance

        complex(dp), allocatable :: bc(:,:)
        real(dp), allocatable :: packed(:,:), diagonal(:), lu(:,:), shifted(:,:)
        real(dp), allocatable :: rhs(:), x(:), ax(:)
        integer, allocatable :: ipiv(:), ordering(:), scratch(:)
        real(dp) :: sigma, value, amplitude, shift, probe, jitter, xnorm, residual
        real(dp) :: bound, coeff, projection
        integer :: n, kd, ldab, m, j, k, i, it, attempt, pass, st, idx
        logical :: low, diagonal_only, good

        w = 0.0_dp
        z = 0.0_dp
        if (present(num_selected)) num_selected = 0
        status = linalg_status_invalid_shape
        n = size(ab,2)
        ldab = size(ab,1)-1
        if (ldab < 0 .or. size(z,1) /= n) return
        allocate(bc(size(ab,1),n))
        bc = cmplx(ab,0.0_dp,dp)
        call selected_banded_sturm(bc,w,status,lower,max_iterations, &
            subset_by_index,subset_by_value,m,tol)
        if (present(num_selected)) num_selected = m
        if (status /= linalg_status_success) return
        status = linalg_status_invalid_shape
        if (size(z,2) < m) then
            w = 0.0_dp
            return
        end if
        if (m == 0) then
            status = linalg_status_success
            return
        end if
        kd = min(ldab,n-1)
        low = .false.
        if (present(lower)) low = lower
        allocate(diagonal(n),packed(2*kd+1,n))
        packed = 0.0_dp
        sigma = 0.0_dp
        diagonal_only = .true.
        do j = 1, n
            if (low) then
                value = ab(1,j)
            else
                value = ab(ldab+1,j)
            end if
            diagonal(j) = value
            packed(kd+1,j) = value
            sigma = max(sigma,abs(value))
            do k = 1, min(kd,n-j)
                if (low) then
                    value = ab(k+1,j)
                else
                    value = ab(ldab+1-k,j+k)
                end if
                packed(kd+1+k,j) = value
                packed(kd+1-k,j+k) = value
                sigma = max(sigma,abs(value))
                if (value /= 0.0_dp) diagonal_only = .false.
            end do
        end do
        if (diagonal_only) then
            allocate(ordering(n),scratch(n))
            call stable_diagonal_order(diagonal,ordering,scratch)
            idx = 1
            do j = 1, m
                do while (idx <= n)
                    if (diagonal(ordering(idx)) == w(j)) exit
                    idx = idx+1
                end do
                if (idx > n) then
                    status = linalg_status_no_convergence
                    w = 0.0_dp
                    z = 0.0_dp
                    return
                end if
                z(ordering(idx),j) = 1.0_dp
                idx = idx+1
            end do
            status = linalg_status_success
            return
        end if
        if (.not. ieee_is_finite(sigma) .or. sigma <= 0.0_dp) then
            status = linalg_status_no_convergence
            w = 0.0_dp
            return
        end if
        packed = packed/sigma
        allocate(shifted(2*kd+1,n),lu(3*kd+1,n),ipiv(n))
        allocate(rhs(n),x(n),ax(n))
        bound = max(1.0_dp,real(2*kd+1,dp))
        status = linalg_status_no_convergence
        do j = 1, m
            shift = w(j)/sigma
            if (.not. ieee_is_finite(shift)) exit
            good = .false.
            do attempt = 1, 8
                jitter = 16.0_dp*epsilon(1.0_dp)*real(attempt,dp)*max(1.0_dp,abs(shift))
                probe = shift+(-1.0_dp)**attempt*jitter
                shifted = packed
                shifted(kd+1,:) = shifted(kd+1,:)-probe
                call lu_factor_banded(shifted,kd,kd,lu,ipiv,st)
                if (st /= linalg_status_success) cycle
                do i = 1, n
                    coeff = real(i,dp)*real(j+3*attempt,dp)
                    rhs(i) = sin(coeff*0.681313_dp)+cos(coeff*0.318309_dp)
                end do
                xnorm = sqrt(sum(rhs*rhs))
                if (xnorm == 0.0_dp) cycle
                rhs = rhs/xnorm
                do it = 1, 20
                    call lu_solve_banded(lu,kd,kd,ipiv,rhs,x,st)
                    if (st /= linalg_status_success) exit
                    if (.not. all(ieee_is_finite(x))) exit
                    amplitude = maxval(abs(x))
                    if (amplitude <= tiny(1.0_dp)) exit
                    x = x/amplitude
                    do pass = 1, 2
                        do k = 1, j-1
                            projection = dot_product(z(:,k),x)
                            x = x-projection*z(:,k)
                        end do
                    end do
                    xnorm = sqrt(sum(x*x))
                    if (xnorm <= 1.0e-100_dp) exit
                    x = x/xnorm
                    call banded_symmetric_matvec(packed,kd,x,ax)
                    ax = ax-shift*x
                    residual = sqrt(sum(ax*ax))
                    if (residual <= max(8.0e-11_dp, &
                        1024.0_dp*epsilon(1.0_dp)*sqrt(real(n,dp)))*bound) then
                        z(:,j) = x
                        good = .true.
                        exit
                    end if
                    rhs = x
                end do
                if (good) exit
            end do
            if (.not. good) exit
        end do
        if (j > m) then
            status = linalg_status_success
        else
            w = 0.0_dp
            z = 0.0_dp
        end if
    end subroutine selected_banded_real_vectors

    subroutine banded_symmetric_matvec(packed,kd,x,y)
        real(dp), intent(in) :: packed(:, :) !! diagonal-centered symmetric band
        integer, intent(in) :: kd !! number of subdiagonals
        real(dp), intent(in) :: x(:) !! matrix-vector operand
        real(dp), intent(out) :: y(:) !! matrix-vector product
        integer :: n, j, k
        real(dp) :: value

        n = size(x)
        y = 0.0_dp
        do j = 1, n
            y(j) = y(j)+packed(kd+1,j)*x(j)
            do k = 1, min(kd,n-j)
                value = packed(kd+1+k,j)
                y(j+k) = y(j+k)+value*x(j)
                y(j) = y(j)+value*x(j+k)
            end do
        end do
    end subroutine banded_symmetric_matvec

    subroutine eig_banded_complex(ab, w, status, z, lower, max_iterations, &
        subset_by_index, subset_by_value, num_selected, tol, method)
        complex(dp), intent(in) :: ab(:, :) !! SciPy-packed Hermitian bands, (kd+1,n)
        real(dp), intent(out) :: w(:) !! selected ascending eigenvalues, capacity >= count
        integer, intent(out) :: status !! success, invalid shape, or failed convergence
        complex(dp), intent(out), optional :: z(:, :) !! eigenvectors (n,capacity)
        logical, intent(in), optional :: lower !! lower-band storage, default false
        integer, intent(in), optional :: max_iterations !! QL or bisection limit
        integer, intent(in), optional :: subset_by_index(2) !! one-based inclusive indices
        real(dp), intent(in), optional :: subset_by_value(2) !! half-open interval (low,high]
        integer, intent(out), optional :: num_selected !! number of selected eigenpairs
        real(dp), intent(in), optional :: tol !! absolute bisection tolerance, >=0
        character(len=*), intent(in), optional :: method !! auto, inverse, bulge, or dense

        complex(dp), allocatable :: vectors(:, :)
        real(dp), allocatable :: values(:)
        character(len=:), allocatable :: strategy
        integer :: n, m, first, last
        logical :: selected, full_inverse

        status = linalg_status_invalid_shape
        w = 0.0_dp
        if (present(z)) z = cmplx(0.0_dp,0.0_dp,dp)
        if (present(num_selected)) num_selected = 0
        n = size(ab,2)
        if (size(ab,1) < 1) return
        if (present(z)) then
            if (size(z,1) /= n) return
        end if
        if (present(subset_by_index) .and. present(subset_by_value)) return
        strategy = 'auto'
        if (present(method)) strategy = trim(adjustl(method))
        if (strategy /= 'auto' .and. strategy /= 'inverse' .and. &
            strategy /= 'dense' .and. strategy /= 'bulge') return
        selected = present(subset_by_index) .or. present(subset_by_value)
        full_inverse = strategy == 'inverse' .and. .not. selected .and. n > 0
        if ((selected .or. full_inverse) .and. strategy /= 'dense' .and. &
            strategy /= 'bulge') then
            if (present(z)) then
                if (full_inverse) then
                    call selected_banded_eigenvectors(ab,w,status,z,lower=lower, &
                        max_iterations=max_iterations,subset_by_index=[1,n], &
                        num_selected=num_selected,tol=tol,real_input=.false.)
                else
                    call selected_banded_eigenvectors(ab,w,status,z,lower,max_iterations, &
                        subset_by_index,subset_by_value,num_selected,tol,.false.)
                end if
            else
                if (full_inverse) then
                    call selected_banded_sturm(ab,w,status,lower=lower, &
                        max_iterations=max_iterations,subset_by_index=[1,n], &
                        num_selected=num_selected,tol=tol)
                else
                    call selected_banded_sturm(ab,w,status,lower,max_iterations, &
                        subset_by_index,subset_by_value,num_selected,tol)
                end if
            end if
            if (full_inverse) return
            if (status /= linalg_status_no_convergence .or. &
                strategy == 'inverse' .or. n > 96) return
            w = 0.0_dp
            if (present(z)) z = cmplx(0.0_dp,0.0_dp,dp)
            if (present(num_selected)) num_selected = 0
        end if
        if (.not. selected) then
            if (strategy == 'bulge' .or. &
                (strategy == 'auto' .and. .not. present(z) .and. &
                auto_bulge_full_values(n,size(ab,1)-1))) then
                call eig_banded_complex_bulge(ab,w,status,z,lower,max_iterations)
                if (strategy == 'auto' .and. status == linalg_status_no_convergence) then
                    call eig_banded_complex_dense(ab,w,status,z,lower,max_iterations)
                end if
            else
                call eig_banded_complex_dense(ab,w,status,z,lower,max_iterations)
            end if
            if (status == linalg_status_success .and. present(num_selected)) num_selected=n
            return
        end if
        allocate(values(n))
        if (present(z)) then
            allocate(vectors(n,n))
            if (strategy == 'bulge') then
                call eig_banded_complex_bulge(ab,values,status,vectors,lower,max_iterations)
            else
                call eig_banded_complex_dense(ab,values,status,vectors,lower,max_iterations)
            end if
        else
            if (strategy == 'bulge') then
                call eig_banded_complex_bulge(ab,values,status,lower=lower, &
                    max_iterations=max_iterations)
            else
                call eig_banded_complex_dense(ab,values,status,lower=lower, &
                    max_iterations=max_iterations)
            end if
        end if
        if (status /= linalg_status_success) return
        call band_selection_indices(values,first,last,status,subset_by_index,subset_by_value)
        if (status /= linalg_status_success) return
        m = last-first+1
        if (present(num_selected)) num_selected=m
        status = linalg_status_invalid_shape
        if (size(w) < m) return
        if (present(z)) then
            if (size(z,2) < m) return
        end if
        if (m > 0) then
            w(1:m)=values(first:last)
            if (present(z)) z(:,1:m)=vectors(:,first:last)
        end if
        status = linalg_status_success
    end subroutine eig_banded_complex

    subroutine band_selection_indices(values,first,last,status,subset_by_index,subset_by_value)
        real(dp), intent(in) :: values(:) !! sorted eigenvalues
        integer, intent(out) :: first !! first selected one-based index
        integer, intent(out) :: last !! last selected one-based index
        integer, intent(out) :: status !! success or invalid interval
        integer, intent(in), optional :: subset_by_index(2) !! inclusive one-based limits
        real(dp), intent(in), optional :: subset_by_value(2) !! interval (low,high]
        integer :: j, n

        n = size(values)
        first = 1
        last = n
        status = linalg_status_invalid_shape
        if (present(subset_by_index) .and. present(subset_by_value)) return
        if (present(subset_by_index)) then
            if (subset_by_index(1) < 1 .or. subset_by_index(2) > n) return
            if (subset_by_index(1) > subset_by_index(2)) return
            first = subset_by_index(1)
            last = subset_by_index(2)
        end if
        if (present(subset_by_value)) then
            if (any(ieee_is_nan(subset_by_value))) return
            if (subset_by_value(1) >= subset_by_value(2)) return
            first = n+1
            do j=1,n
                if (values(j) > subset_by_value(1)) then
                    first=j
                    exit
                end if
            end do
            last=first-1
            do j=first,n
                if (values(j) > subset_by_value(2)) exit
                last=j
            end do
        end if
        status = linalg_status_success
    end subroutine band_selection_indices

    subroutine selected_banded_eigenvectors(ab,w,status,z,lower,max_iterations, &
        subset_by_index,subset_by_value,num_selected,tol,real_input)
        ! Sturm eigenvalues and shifted partial-pivoted band LU inverse iteration.
        ! Workspace: O(n*kd + n*m) for m requested eigenvectors.
        complex(dp), intent(in) :: ab(:, :) !! SciPy-packed Hermitian matrix
        real(dp), intent(out) :: w(:) !! sorted eigenvalues, capacity >= m
        integer, intent(out) :: status !! success, invalid shape, or no convergence
        complex(dp), intent(out) :: z(:, :) !! selected orthonormal eigenvectors, (n,capacity)
        logical, intent(in), optional :: lower !! lower diagonal packing
        integer, intent(in), optional :: max_iterations !! Sturm bisection budget
        integer, intent(in), optional :: subset_by_index(2) !! inclusive one-based limits
        real(dp), intent(in), optional :: subset_by_value(2) !! interval (low,high]
        integer, intent(out), optional :: num_selected !! number of requested eigenvectors
        real(dp), intent(in), optional :: tol !! nonnegative bisection tolerance
        logical, intent(in) :: real_input !! real symmetric rather than complex Hermitian

        complex(dp), allocatable :: packed(:, :), shifted_band(:, :), lu(:, :)
        complex(dp), allocatable :: rhs(:), x(:), ax(:), diag_band(:, :)
        real(dp), allocatable :: diagonal(:)
        real(dp), allocatable :: shifted_real(:, :), lu_real(:, :), rhs_real(:), x_real(:)
        integer, allocatable :: ipiv(:), ordering(:), work(:)
        real(dp) :: sigma, amplitude, shift, probe, jitter, xnorm, residual
        real(dp) :: bound, magnitude, coeff
        complex(dp) :: value, projection
        integer :: n, kd, ldab, m, j, k, i, it, attempt, pass, st, idx
        logical :: low, diagonal_only, good

        w=0.0_dp
        z=cmplx(0.0_dp,0.0_dp,dp)
        if (present(num_selected)) num_selected=0
        status=linalg_status_invalid_shape
        n=size(ab,2)
        ldab=size(ab,1)-1
        if (ldab < 0 .or. size(z,1) /= n) return
        call selected_banded_sturm(ab,w,status,lower,max_iterations, &
            subset_by_index,subset_by_value,m,tol)
        if (present(num_selected)) num_selected=m
        if (status /= linalg_status_success) return
        status=linalg_status_invalid_shape
        if (size(z,2) < m) then
            w=0.0_dp
            return
        end if
        if (m == 0) then
            status=linalg_status_success
            return
        end if
        kd=min(ldab,n-1)
        low=.false.
        if (present(lower)) low=lower
        allocate(diagonal(n),diag_band(kd+1,n))
        diag_band=cmplx(0.0_dp,0.0_dp,dp)
        sigma=0.0_dp
        diagonal_only=.true.
        do j=1,n
            if (low) then
                value=ab(1,j)
            else
                value=ab(ldab+1,j)
            end if
            diagonal(j)=real(value,dp)
            diag_band(1,j)=cmplx(diagonal(j),0.0_dp,dp)
            sigma=max(sigma,abs(diagonal(j)))
            do k=1,min(kd,n-j)
                if (low) then
                    value=ab(k+1,j)
                else
                    value=conjg(ab(ldab+1-k,j+k))
                end if
                diag_band(k+1,j)=value
                sigma=max(sigma,abs(value))
                if (value /= cmplx(0.0_dp,0.0_dp,dp)) diagonal_only=.false.
            end do
        end do
        if (diagonal_only) then
            allocate(ordering(n),work(n))
            call stable_diagonal_order(diagonal,ordering,work)
            idx=1
            do j=1,m
                do while (idx <= n)
                    if (diagonal(ordering(idx)) == w(j)) exit
                    idx=idx+1
                end do
                if (idx > n) then
                    status=linalg_status_no_convergence
                    w=0.0_dp
                    z=cmplx(0.0_dp,0.0_dp,dp)
                    return
                end if
                z(ordering(idx),j)=cmplx(1.0_dp,0.0_dp,dp)
                idx=idx+1
            end do
            status=linalg_status_success
            return
        end if
        if (sigma <= 0.0_dp .or. .not. ieee_is_finite(sigma)) then
            status=linalg_status_no_convergence
            w=0.0_dp
            return
        end if
        allocate(packed(2*kd+1,n))
        if (real_input) then
            allocate(shifted_real(2*kd+1,n),lu_real(3*kd+1,n))
            allocate(rhs_real(n),x_real(n))
        else
            allocate(shifted_band(2*kd+1,n),lu(3*kd+1,n))
        end if
        allocate(ipiv(n),rhs(n),x(n),ax(n))
        packed=cmplx(0.0_dp,0.0_dp,dp)
        do j=1,n
            packed(kd+1,j)=diag_band(1,j)/sigma
            do k=1,min(kd,n-j)
                value=diag_band(k+1,j)/sigma
                packed(kd+1+k,j)=value
                packed(kd+1-k,j+k)=conjg(value)
            end do
        end do
        bound = max(1.0_dp,real(2*kd+1,dp))
        status=linalg_status_no_convergence
        do j=1,m
            shift=w(j)/sigma
            if (.not. ieee_is_finite(shift)) exit
            good=.false.
            do attempt=1,8
                jitter=16.0_dp*epsilon(1.0_dp)*real(attempt,dp)*max(1.0_dp,abs(shift))
                probe=shift + (-1.0_dp)**attempt*jitter
                if (real_input) then
                    ! Real symmetric bands use real factorization and storage.
                    shifted_real=real(packed,dp)
                    shifted_real(kd+1,:)=shifted_real(kd+1,:)-probe
                    call lu_factor_banded(shifted_real,kd,kd,lu_real,ipiv,st)
                else
                    shifted_band=packed
                    shifted_band(kd+1,:)=shifted_band(kd+1,:)-probe
                    call lu_factor_banded(shifted_band,kd,kd,lu,ipiv,st)
                end if
                if (st /= linalg_status_success) cycle
                do i=1,n
                    coeff=real(i,dp)*real(j+3*attempt,dp)
                    rhs(i)=cmplx(sin(coeff*0.681313_dp)+cos(coeff*0.318309_dp), &
                        0.0_dp,dp)
                    if (.not. real_input) then
                        rhs(i)=rhs(i)+cmplx(0.0_dp,sin(coeff*0.4173_dp),dp)
                    end if
                end do
                xnorm=sqrt(sum(abs(rhs)**2))
                if (xnorm==0.0_dp) cycle
                rhs=rhs/xnorm
                do it=1,20
                    if (real_input) then
                        rhs_real=real(rhs,dp)
                        call lu_solve_banded(lu_real,kd,kd,ipiv,rhs_real,x_real,st)
                        if (st /= linalg_status_success) exit
                        x=cmplx(x_real,0.0_dp,dp)
                    else
                        call lu_solve_banded(lu,kd,kd,ipiv,rhs,x,st)
                        if (st /= linalg_status_success) exit
                    end if
                    if (.not. all(ieee_is_finite(real(x,dp))) .or. &
                        .not. all(ieee_is_finite(aimag(x)))) exit
                    amplitude=maxval(abs(x))
                    if (amplitude <= tiny(1.0_dp)) exit
                    x=x/amplitude
                    do pass=1,2
                        do k=1,j-1
                            projection=dot_product(z(:,k),x)
                            x=x-projection*z(:,k)
                        end do
                    end do
                    xnorm=sqrt(sum(abs(x)**2))
                    if (xnorm <= 1.0e-100_dp) exit
                    x=x/xnorm
                    call banded_hermitian_matvec(packed,kd,x,ax)
                    ax=ax-shift*x
                    residual=sqrt(sum(abs(ax)**2))
                    if (residual <= max(8.0e-11_dp,1024.0_dp*epsilon(1.0_dp)*sqrt(real(n,dp))) &
                        *bound) then
                        z(:,j)=x
                        good=.true.
                        exit
                    end if
                    rhs=x
                end do
                if (good) exit
            end do
            if (.not. good) exit
        end do
        if (j > m) then
            status=linalg_status_success
        else
            w=0.0_dp
            z=cmplx(0.0_dp,0.0_dp,dp)
        end if
    end subroutine selected_banded_eigenvectors

    subroutine banded_hermitian_matvec(packed,kd,x,y)
        complex(dp), intent(in) :: packed(:, :) !! compact symmetric (2*kd+1,n) band
        integer, intent(in) :: kd !! half-bandwidth
        complex(dp), intent(in) :: x(:) !! vector length n
        complex(dp), intent(out) :: y(:) !! matrix product of length n
        integer :: n,j,k
        complex(dp) :: value

        n=size(x)
        y=cmplx(0.0_dp,0.0_dp,dp)
        do j=1,n
            y(j)=y(j)+packed(kd+1,j)*x(j)
            do k=1,min(kd,n-j)
                value=packed(kd+1+k,j)
                y(j+k)=y(j+k)+value*x(j)
                y(j)=y(j)+conjg(value)*x(j+k)
            end do
        end do
    end subroutine banded_hermitian_matvec

    subroutine eig_banded_real_bulge(ab, w, status, z, lower, max_iterations)
        real(dp), intent(in) :: ab(:, :) !! symmetric band, SciPy packing
        real(dp), intent(out) :: w(:) !! ascending eigenvalues, length n
        integer, intent(out) :: status !! success or numerical failure
        real(dp), intent(out), optional :: z(:, :) !! eigenvectors, shape (n,n)
        logical, intent(in), optional :: lower !! lower triangle, default false
        integer, intent(in), optional :: max_iterations !! QL step limit

        real(dp), allocatable :: d(:), e(:), q(:, :), u(:, :)
        integer :: n

        n = size(ab,2)
        w = 0.0_dp
        if (present(z)) z = 0.0_dp
        status = linalg_status_invalid_shape
        if (size(w) /= n) return
        if (present(z)) then
            if (size(z,1) /= n .or. size(z,2) /= n) return
        end if
        allocate(d(n),e(max(0,n-1)))
        if (present(z)) then
            allocate(q(n,n),u(n,n))
            call band_to_tridiagonal(ab,d,e,status,q,lower)
            if (status /= linalg_status_success) return
            call eigh_tridiagonal(d,e,w,status,u,max_iterations)
            if (status == linalg_status_success) z = matmul(q,u)
        else
            call band_to_tridiagonal(ab,d,e,status,lower=lower)
            if (status /= linalg_status_success) return
            call eigh_tridiagonal(d,e,w,status,max_iterations=max_iterations)
        end if
    end subroutine eig_banded_real_bulge

    subroutine eig_banded_complex_bulge(ab, w, status, z, lower, max_iterations)
        complex(dp), intent(in) :: ab(:, :) !! Hermitian band, SciPy packing
        real(dp), intent(out) :: w(:) !! ascending eigenvalues, length n
        integer, intent(out) :: status !! success or numerical failure
        complex(dp), intent(out), optional :: z(:, :) !! eigenvectors, shape (n,n)
        logical, intent(in), optional :: lower !! lower triangle, default false
        integer, intent(in), optional :: max_iterations !! QL step limit

        real(dp), allocatable :: d(:), e(:), u(:, :)
        complex(dp), allocatable :: q(:, :)
        integer :: n

        n = size(ab,2)
        w = 0.0_dp
        if (present(z)) z = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (size(w) /= n) return
        if (present(z)) then
            if (size(z,1) /= n .or. size(z,2) /= n) return
        end if
        allocate(d(n),e(max(0,n-1)))
        if (present(z)) then
            allocate(q(n,n),u(n,n))
            call band_to_tridiagonal(ab,d,e,status,q,lower)
            if (status /= linalg_status_success) return
            call eigh_tridiagonal(d,e,w,status,u,max_iterations)
            if (status == linalg_status_success) z = matmul(q,cmplx(u,0.0_dp,dp))
        else
            call band_to_tridiagonal(ab,d,e,status,lower=lower)
            if (status /= linalg_status_success) return
            call eigh_tridiagonal(d,e,w,status,max_iterations=max_iterations)
        end if
    end subroutine eig_banded_complex_bulge

    subroutine eig_banded_real_dense(ab, w, status, z, lower, max_iterations)
        real(dp), intent(in) :: ab(:, :) !! upper/lower diagonal-order bands, shape (kd+1,n)
        real(dp), intent(out) :: w(:) !! sorted ascending eigenvalues, length n
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(out), optional :: z(:, :) !! orthonormal eigenvectors, (n,n)
        logical, intent(in), optional :: lower !! lower packed form, default false (upper)
        integer, intent(in), optional :: max_iterations !! tridiagonal QL limit, default 120

        real(dp), allocatable :: a(:, :), h(:, :), q(:, :), u(:, :), d(:), e(:)
        integer :: n, kd, i, j, k
        logical :: low

        w = 0.0_dp
        if (present(z)) z = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(ab, 2)
        kd = size(ab, 1) - 1
        if (kd < 0 .or. size(w) /= n) return
        if (present(z)) then
            if (size(z,1) /= n .or. size(z,2) /= n) return
        end if
        low = .false.
        if (present(lower)) low = lower
        allocate(a(n,n),h(n,n),d(n),e(max(0,n-1)))
        a = 0.0_dp
        status = linalg_status_no_convergence
        do j = 1, n
            do k = 0, min(kd, merge(n-j,j-1,low))
                if (low) then
                    i = j + k
                    if (.not. ieee_is_finite(ab(1+k,j))) return
                    a(i,j) = ab(1+k,j)
                    a(j,i) = ab(1+k,j)
                else
                    i = j - k
                    if (.not. ieee_is_finite(ab(kd+1-k,j))) return
                    a(i,j) = ab(kd+1-k,j)
                    a(j,i) = ab(kd+1-k,j)
                end if
            end do
        end do
        if (present(z)) then
            allocate(q(n,n),u(n,n))
            call hessenberg(a,h,status,q)
        else
            call hessenberg(a,h,status)
        end if
        if (status /= linalg_status_success) return
        do i = 1, n
            d(i) = h(i,i)
        end do
        do i = 1, n-1
            e(i) = 0.5_dp*h(i+1,i) + 0.5_dp*h(i,i+1)
        end do
        if (present(z)) then
            call eigh_tridiagonal(d,e,w,status,u,max_iterations)
            if (status == linalg_status_success) z = matmul(q,u)
        else
            call eigh_tridiagonal(d,e,w,status,max_iterations=max_iterations)
        end if
    end subroutine eig_banded_real_dense

    subroutine eig_banded_complex_dense(ab, w, status, z, lower, max_iterations)
        complex(dp), intent(in) :: ab(:, :) !! Hermitian bands, shape (kd+1,n)
        real(dp), intent(out) :: w(:) !! sorted ascending eigenvalues, length n
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: z(:, :) !! unitary eigenvectors, shape (n,n)
        logical, intent(in), optional :: lower !! lower packed form, default false
        integer, intent(in), optional :: max_iterations !! tridiagonal QL limit, default 120

        complex(dp), allocatable :: a(:, :), h(:, :), q(:, :), phases(:), off(:)
        real(dp), allocatable :: d(:), e(:), u(:, :)
        real(dp) :: scalar, mag
        integer :: n, kd, i, j, k
        logical :: low

        w = 0.0_dp
        if (present(z)) z = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(ab,2)
        kd = size(ab,1) - 1
        if (kd < 0 .or. size(w) /= n) return
        if (present(z)) then
            if (size(z,1) /= n .or. size(z,2) /= n) return
        end if
        low = .false.
        if (present(lower)) low = lower
        allocate(a(n,n),h(n,n),d(n),e(max(0,n-1)),phases(n),off(max(0,n-1)))
        a = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_no_convergence
        do j = 1, n
            do k = 0, min(kd, merge(n-j,j-1,low))
                if (low) then
                    i = j+k
                    if (.not. finite_complex_scalar(ab(1+k,j))) return
                    a(i,j) = ab(1+k,j)
                else
                    i = j-k
                    if (.not. finite_complex_scalar(ab(kd+1-k,j))) return
                    a(i,j) = ab(kd+1-k,j)
                end if
                if (i == j) then
                    ! The imaginary diagonal is ignored, matching Hermitian storage.
                    a(i,j) = cmplx(real(a(i,j),dp),0.0_dp,dp)
                else
                    a(j,i) = conjg(a(i,j))
                end if
            end do
        end do
        if (present(z)) then
            allocate(q(n,n),u(n,n))
            call hessenberg(a,h,status,q)
        else
            call hessenberg(a,h,status)
        end if
        if (status /= linalg_status_success) return
        do i = 1, n
            d(i) = real(h(i,i),dp)
        end do
        if (n >= 1) phases(1) = cmplx(1.0_dp,0.0_dp,dp)
        do i = 1, n-1
            off(i) = 0.5_dp*h(i+1,i) + 0.5_dp*conjg(h(i,i+1))
            scalar = max(abs(real(off(i),dp)),abs(aimag(off(i))))
            if (scalar == 0.0_dp) then
                e(i) = 0.0_dp
                phases(i+1) = phases(i)
            else
                mag = hypot(real(off(i),dp)/scalar,aimag(off(i))/scalar)
                e(i) = scalar*mag
                if (.not. ieee_is_finite(e(i))) return
                phases(i+1) = phases(i)*(off(i)/scalar)/mag
            end if
        end do
        if (present(z)) then
            call eigh_tridiagonal(d,e,w,status,u,max_iterations)
            if (status /= linalg_status_success) return
            do j = 1, n
                q(:,j) = q(:,j)*phases(j)
            end do
            z = matmul(q,cmplx(u,0.0_dp,dp))
        else
            call eigh_tridiagonal(d,e,w,status,max_iterations=max_iterations)
        end if
    end subroutine eig_banded_complex_dense

    subroutine eigvals_banded_real(ab, w, status, lower, max_iterations, &
        subset_by_index, subset_by_value, num_selected, tol, method)
        real(dp), intent(in) :: ab(:, :) !! diagonal-ordered real symmetric bands
        real(dp), intent(out) :: w(:) !! eigenvalues; selected-capacity array for subsets
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! lower packing, default false
        integer, intent(in), optional :: max_iterations !! iteration limit, default automatic
        integer, intent(in), optional :: subset_by_index(2) !! one-based inclusive index selection
        real(dp), intent(in), optional :: subset_by_value(2) !! half-open (low,high] interval
        integer, intent(out), optional :: num_selected !! selected eigenvalue count
        real(dp), intent(in), optional :: tol !! bisection absolute tolerance, >=0
        character(len=*), intent(in), optional :: method !! auto (default), sturm, bulge, dense

        complex(dp), allocatable :: bc(:, :)
        real(dp), allocatable :: full(:)
        character(len=:), allocatable :: strategy
        integer :: n

        w = 0.0_dp
        if (present(num_selected)) num_selected = 0
        status = linalg_status_invalid_shape
        strategy = 'auto'
        if (present(method)) strategy = trim(adjustl(method))
        if (strategy /= 'auto' .and. strategy /= 'dense' .and. &
            strategy /= 'sturm' .and. strategy /= 'bulge') return
        n = size(ab,2)
        if (strategy == 'sturm' .or. &
            (strategy == 'auto' .and. &
            (present(subset_by_index) .or. present(subset_by_value)))) then
            allocate(bc(size(ab,1),n))
            bc = cmplx(ab,0.0_dp,dp)
            call selected_banded_sturm(bc,w,status,lower,max_iterations, &
                subset_by_index,subset_by_value,num_selected,tol)
        else
            allocate(full(n))
            if (strategy == 'bulge' .or. &
                (strategy == 'auto' .and. .not. present(subset_by_index) .and. &
                .not. present(subset_by_value) .and. &
                auto_bulge_full_values(n,size(ab,1)-1))) then
                call eig_banded_real_bulge(ab,full,status,lower=lower, &
                    max_iterations=max_iterations)
                if (strategy == 'auto' .and. status == linalg_status_no_convergence) then
                    call eig_banded_real_dense(ab,full,status,lower=lower, &
                        max_iterations=max_iterations)
                end if
            else
                call eig_banded_real_dense(ab,full,status,lower=lower,max_iterations=max_iterations)
            end if
            if (status /= linalg_status_success) return
            call select_dense_band_values(full,w,status,subset_by_index, &
                subset_by_value,num_selected)
        end if
    end subroutine eigvals_banded_real

    subroutine eigvals_banded_complex(ab, w, status, lower, max_iterations, &
        subset_by_index, subset_by_value, num_selected, tol, method)
        complex(dp), intent(in) :: ab(:, :) !! diagonal-ordered Hermitian bands
        real(dp), intent(out) :: w(:) !! eigenvalues; selected-capacity array for subsets
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! lower packing, default false
        integer, intent(in), optional :: max_iterations !! iteration limit, default automatic
        integer, intent(in), optional :: subset_by_index(2) !! one-based inclusive index selection
        real(dp), intent(in), optional :: subset_by_value(2) !! half-open (low,high] interval
        integer, intent(out), optional :: num_selected !! selected eigenvalue count
        real(dp), intent(in), optional :: tol !! bisection absolute tolerance, >=0
        character(len=*), intent(in), optional :: method !! auto (default), sturm, bulge, dense

        real(dp), allocatable :: full(:)
        character(len=:), allocatable :: strategy
        integer :: n

        w = 0.0_dp
        if (present(num_selected)) num_selected = 0
        status = linalg_status_invalid_shape
        strategy = 'auto'
        if (present(method)) strategy = trim(adjustl(method))
        if (strategy /= 'auto' .and. strategy /= 'dense' .and. &
            strategy /= 'sturm' .and. strategy /= 'bulge') return
        n = size(ab,2)
        if (strategy == 'sturm' .or. &
            (strategy == 'auto' .and. &
            (present(subset_by_index) .or. present(subset_by_value)))) then
            call selected_banded_sturm(ab,w,status,lower,max_iterations, &
                subset_by_index,subset_by_value,num_selected,tol)
        else
            allocate(full(n))
            if (strategy == 'bulge' .or. &
                (strategy == 'auto' .and. .not. present(subset_by_index) .and. &
                .not. present(subset_by_value) .and. &
                auto_bulge_full_values(n,size(ab,1)-1))) then
                call eig_banded_complex_bulge(ab,full,status,lower=lower, &
                    max_iterations=max_iterations)
                if (strategy == 'auto' .and. status == linalg_status_no_convergence) then
                    call eig_banded_complex_dense(ab,full,status,lower=lower, &
                        max_iterations=max_iterations)
                end if
            else
                call eig_banded_complex_dense(ab,full,status,lower=lower,max_iterations=max_iterations)
            end if
            if (status /= linalg_status_success) return
            call select_dense_band_values(full,w,status,subset_by_index, &
                subset_by_value,num_selected)
        end if
    end subroutine eigvals_banded_complex

    subroutine selected_banded_sturm(ab, w, status, lower, max_iterations, &
        subset_by_index, subset_by_value, num_selected, tol)
        ! Inertia of a Hermitian band via unpivoted LDL^H, followed by bisection.
        ! Only (kd+1)*n storage is needed. No eigenvectors are constructed.
        complex(dp), intent(in) :: ab(:, :) !! SciPy diagonal-ordered Hermitian bands
        real(dp), intent(out) :: w(:) !! sorted selected eigenvalues, capacity >= selection
        integer, intent(out) :: status !! success, invalid shape or nonconvergence
        logical, intent(in), optional :: lower !! lower packed bands if true
        integer, intent(in), optional :: max_iterations !! maximum bisection steps, default 140
        integer, intent(in), optional :: subset_by_index(2) !! inclusive one-based range
        real(dp), intent(in), optional :: subset_by_value(2) !! selected interval (low,high]
        integer, intent(out), optional :: num_selected !! selected multiplicity
        real(dp), intent(in), optional :: tol !! nonnegative absolute bisection tolerance

        complex(dp), allocatable :: band(:, :), scratch(:, :)
        real(dp), allocatable :: sorted(:), radii(:)
        real(dp) :: sigma, lo, hi, left, right, middle, margin, eigval, width_tol, pivmin
        real(dp) :: limit, bound, actual
        integer :: n, ldab, kd, j, k, first, last, selected, max_steps, step, count
        logical :: low, good, converged, diagonal

        status = linalg_status_invalid_shape
        w = 0.0_dp
        if (present(num_selected)) num_selected = 0
        n = size(ab,2)
        ldab = size(ab,1)-1
        if (ldab < 0) return
        if (present(subset_by_index) .and. present(subset_by_value)) return
        if (present(subset_by_index)) then
            if (subset_by_index(1) < 1 .or. subset_by_index(2) > n) return
            if (subset_by_index(1) > subset_by_index(2)) return
        end if
        if (present(subset_by_value)) then
            if (any(ieee_is_nan(subset_by_value))) return
            if (subset_by_value(1) >= subset_by_value(2)) return
        end if
        max_steps = 140
        if (present(max_iterations)) max_steps = max_iterations
        if (max_steps < 1) return
        if (present(tol)) then
            if (.not. ieee_is_finite(tol)) return
            if (tol < 0.0_dp) return
        end if
        if (.not. present(subset_by_index) .and. .not. present(subset_by_value)) then
            if (size(w) /= n) return
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        low = .false.
        if (present(lower)) low = lower
        kd = min(ldab,n-1)
        allocate(band(kd+1,n),scratch(kd+1,n),radii(n))
        band = cmplx(0.0_dp,0.0_dp,dp)
        sigma = 0.0_dp
        diagonal = .true.
        status = linalg_status_no_convergence
        do j = 1, n
            if (low) then
                actual = real(ab(1,j),dp)
                if (.not. finite_complex_scalar(ab(1,j))) return
            else
                actual = real(ab(ldab+1,j),dp)
                if (.not. finite_complex_scalar(ab(ldab+1,j))) return
            end if
            band(1,j) = cmplx(actual,0.0_dp,dp)
            sigma = max(sigma,abs(actual))
            do k = 1, min(kd,n-j)
                if (low) then
                    if (.not. finite_complex_scalar(ab(k+1,j))) return
                    band(k+1,j) = ab(k+1,j)
                else
                    if (.not. finite_complex_scalar(ab(ldab+1-k,j+k))) return
                    band(k+1,j) = conjg(ab(ldab+1-k,j+k))
                end if
                sigma = max(sigma,abs(band(k+1,j)))
                if (band(k+1,j) /= cmplx(0.0_dp,0.0_dp,dp)) diagonal = .false.
            end do
        end do

        ! Preserve subnormal and extreme diagonal values exactly when no
        ! off-diagonal couplings exist: global normalization could lose them.
        if (diagonal) then
            allocate(sorted(n))
            sorted = real(band(1,:),dp)
            call sort_ascending(sorted)
            first = 1
            last = n
            if (present(subset_by_index)) then
                first = subset_by_index(1)
                last = subset_by_index(2)
            end if
            if (present(subset_by_value)) then
                first = n+1
                do j = 1, n
                    if (sorted(j) > subset_by_value(1)) then
                        first = j
                        exit
                    end if
                end do
                last = first-1
                do j = first, n
                    if (sorted(j) > subset_by_value(2)) exit
                    last = j
                end do
            end if
            selected = last-first+1
            if (present(num_selected)) num_selected = selected
            status = linalg_status_invalid_shape
            if (size(w) < selected) return
            if (selected > 0) w(1:selected) = sorted(first:last)
            status = linalg_status_success
            return
        end if
        ! Normalize the entire band to keep LDL pivots near unit scale.
        band = band/sigma
        radii = 0.0_dp
        do j = 1, n
            do k = 1, min(kd,n-j)
                radii(j) = radii(j)+abs(band(k+1,j))
                radii(j+k) = radii(j+k)+abs(band(k+1,j))
            end do
        end do
        lo = huge(1.0_dp)
        hi = -huge(1.0_dp)
        do j = 1, n
            lo = min(lo, real(band(1,j),dp)-radii(j))
            hi = max(hi, real(band(1,j),dp)+radii(j))
        end do
        margin = 128.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lo),abs(hi))
        lo = lo-margin
        hi = hi+margin
        ! Tiny shifted pivots are treated as negative, as in Sturm counting.
        ! This is not a backward-stable indefinite factorization; guard
        ! growth and return a failure status if a finite count is unavailable.
        pivmin = 32.0_dp*epsilon(1.0_dp)
        first = 1
        last = n
        if (present(subset_by_index)) then
            first = subset_by_index(1)
            last = subset_by_index(2)
        end if
        if (present(subset_by_value)) then
            limit = max(4.0_dp,abs(lo),abs(hi))+2.0_dp
            bound = subset_by_value(1)
            if (.not. ieee_is_finite(bound)) then
                actual = sign(limit,bound)
            else if (sigma < 1.0_dp .and. abs(bound) > sigma*limit) then
                actual = sign(limit,bound)
            else
                actual = bound/sigma
            end if
            call banded_inertia(band, scratch, actual, pivmin, count, good)
            if (.not. good) return
            first = count+1
            bound = subset_by_value(2)
            if (.not. ieee_is_finite(bound)) then
                actual = sign(limit,bound)
            else if (sigma < 1.0_dp .and. abs(bound) > sigma*limit) then
                actual = sign(limit,bound)
            else
                actual = bound/sigma
            end if
            call banded_inertia(band, scratch, actual, pivmin, count, good)
            if (.not. good) return
            last = count
        end if
        selected = last-first+1
        if (selected < 0) return
        if (present(num_selected)) num_selected = selected
        status = linalg_status_invalid_shape
        if (size(w) < selected) return
        if (selected <= 0) then
            status = linalg_status_success
            return
        end if

        width_tol = 16.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lo),abs(hi))
        if (present(tol)) then
            if (tol > 0.0_dp) then
                if (sigma < 1.0_dp .and. tol > sigma*max(1.0_dp,hi-lo)) then
                    width_tol = max(width_tol,hi-lo)
                else
                    width_tol = max(width_tol,tol/sigma)
                end if
            end if
        end if
        status = linalg_status_no_convergence
        do j = 1, selected
            k = first+j-1
            left = lo
            right = hi
            converged = .false.
            do step = 1, max_steps
                middle = left+0.5_dp*(right-left)
                if (middle <= left .or. middle >= right) then
                    converged = .true.
                    exit
                end if
                call banded_inertia(band, scratch, middle, pivmin, count, good)
                if (.not. good) then
                    w = 0.0_dp
                    return
                end if
                if (count < k) then
                    left = middle
                else
                    right = middle
                end if
                if (right-left <= width_tol) then
                    converged = .true.
                    exit
                end if
            end do
            if (.not. converged) then
                w = 0.0_dp
                return
            end if
            eigval = (left+0.5_dp*(right-left))*sigma
            if (.not. ieee_is_finite(eigval)) then
                w = 0.0_dp
                return
            end if
            w(j) = eigval
        end do
        status = linalg_status_success
    end subroutine selected_banded_sturm

    subroutine banded_inertia(band, work, x, pivmin, count, good)
        ! Negative inertia of band-x*I via unpivoted Hermitian LDL^H.
        ! Schur updates of a bandwidth-kd matrix retain the same bandwidth.
        complex(dp), intent(in) :: band(:, :) !! normalized lower band, (kd+1,n)
        complex(dp), intent(inout) :: work(:, :) !! LDL scratch of the same shape
        real(dp), intent(in) :: x !! normalized real spectral shift
        real(dp), intent(in) :: pivmin !! regularizer for near-zero 1x1 pivots
        integer, intent(out) :: count !! approximate count of eigenvalues <= x
        logical, intent(out) :: good !! finite factorization and count

        complex(dp) :: ratio, change, old
        real(dp) :: pivot
        integer :: n, kd, k, j, i, last

        n = size(band,2)
        kd = size(band,1)-1
        work = band
        count = 0
        good = .false.
        do k = 1, n
            pivot = real(work(1,k),dp)-x
            if (.not. ieee_is_finite(pivot)) return
            if (abs(pivot) <= pivmin) pivot = -pivmin
            if (pivot < 0.0_dp) count = count+1
            last = min(n,k+kd)
            do j = k+1, last
                ratio = conjg(work(j-k+1,k))/pivot
                if (.not. finite_complex_scalar(ratio)) return
                do i = j, last
                    old = work(i-j+1,j)
                    ! Before multiplication, reject extreme pivot growth
                    ! rather than triggering overflow in a checked build.
                    if (abs(ratio) > 1.0_dp) then
                        if (abs(work(i-k+1,k)) > huge(1.0_dp)/abs(ratio)) return
                    end if
                    change = work(i-k+1,k)*ratio
                    if (.not. finite_complex_scalar(change)) return
                    old = old-change
                    if (.not. finite_complex_scalar(old)) return
                    if (i == j) then
                        work(1,j) = cmplx(real(old,dp),0.0_dp,dp)
                    else
                        work(i-j+1,j) = old
                    end if
                end do
            end do
        end do
        good = .true.
    end subroutine banded_inertia

    subroutine select_dense_band_values(full, w, status, subset_by_index, &
        subset_by_value, num_selected)
        real(dp), intent(in) :: full(:) !! ascending eigenvalues from dense reduction
        real(dp), intent(out) :: w(:) !! selected values; capacity >= selection
        integer, intent(out) :: status !! success or invalid shape
        integer, intent(in), optional :: subset_by_index(2) !! inclusive one-based slice
        real(dp), intent(in), optional :: subset_by_value(2) !! half-open interval
        integer, intent(out), optional :: num_selected !! required capacity
        integer :: n, first, last, j, selected

        n = size(full)
        w = 0.0_dp
        if (present(num_selected)) num_selected = 0
        status = linalg_status_invalid_shape
        if (present(subset_by_index) .and. present(subset_by_value)) return
        first = 1
        last = n
        if (present(subset_by_index)) then
            if (subset_by_index(1) < 1 .or. subset_by_index(2) > n) return
            if (subset_by_index(1) > subset_by_index(2)) return
            first = subset_by_index(1)
            last = subset_by_index(2)
        end if
        if (present(subset_by_value)) then
            if (any(ieee_is_nan(subset_by_value))) return
            if (subset_by_value(1) >= subset_by_value(2)) return
            first = n+1
            do j = 1, n
                if (full(j) > subset_by_value(1)) then
                    first = j
                    exit
                end if
            end do
            last = first-1
            do j = first, n
                if (full(j) > subset_by_value(2)) exit
                last = j
            end do
        end if
        selected = last-first+1
        if (present(num_selected)) num_selected = selected
        if (size(w) < selected) return
        if (selected > 0) w(1:selected) = full(first:last)
        status = linalg_status_success
    end subroutine select_dense_band_values
    pure logical function finite_complex_scalar(x) result(ok)
        complex(dp), intent(in) :: x !! tested value
        ok = ieee_is_finite(real(x,dp)) .and. ieee_is_finite(aimag(x))
    end function finite_complex_scalar

end module scifort_linalg_banded_eigen
