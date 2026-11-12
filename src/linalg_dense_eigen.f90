! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Standard and generalized definite symmetric/Hermitian eigenproblems via Cholesky
! reduction, followed by Householder reduction
! and SciFort's real symmetric tridiagonal implicit-shift QL eigensolver.
! Only the selected triangle of the input matrix is referenced.
! A full spectrum is computed before optional index/value selection.
module scifort_linalg_dense_eigen
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_dense, only : cholesky, solve_triangular
    use scifort_linalg_hessenberg, only : hessenberg
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal
    implicit none
    private

    public :: eigh, eigvalsh

    interface eigh
        module procedure eigh_real
        module procedure eigh_complex
    end interface eigh

    interface eigvalsh
        module procedure eigvalsh_real
        module procedure eigvalsh_complex
    end interface eigvalsh

contains

    recursive subroutine eigh_real(a, w, status, z, lower, subset_by_index, &
        max_iterations, b, problem_type, subset_by_value, num_selected)
        real(dp), intent(in) :: a(:, :) !! square symmetric matrix; unused triangle is ignored
        real(dp), intent(out) :: w(:) !! ascending eigenvalues, capacity in value mode
        integer, intent(out) :: status !! success, invalid shape, or numerical failure
        real(dp), intent(out), optional :: z(:, :) !! eigenvectors, capacity in value mode
        logical, intent(in), optional :: lower !! use lower input triangle; default true
        integer, intent(in), optional :: subset_by_index(2) !! inclusive 1-based [first,last]
        integer, intent(in), optional :: max_iterations !! QL iteration limit; default 120
        real(dp), intent(in), optional :: b(:, :) !! positive-definite mass matrix, default I
        integer, intent(in), optional :: problem_type !! generalized type 1, 2, or 3; default 1
        real(dp), intent(in), optional :: subset_by_value(2) !! (low, high] eigenvalue interval
        integer, intent(out), optional :: num_selected !! number returned; may exceed output capacity

        real(dp), allocatable :: aherm(:, :), h(:, :), q(:, :), u(:, :)
        real(dp), allocatable :: d(:), e(:), allw(:)
        real(dp) :: value
        integer :: n, first, last, m, i, j
        logical :: use_lower

        w = 0.0_dp
        if (present(z)) z = 0.0_dp
        status = linalg_status_invalid_shape
        if (present(num_selected)) num_selected = 0
        if (present(b)) then
            call generalized_real(a, b, w, status, z, lower, &
                subset_by_index, max_iterations, problem_type, subset_by_value, num_selected)
            return
        end if
        n = size(a, 1)
        if (size(a, 2) /= n) return
        if (present(subset_by_index) .and. present(subset_by_value)) return
        call selection_bounds(n, subset_by_index, first, last)
        if (first < 1 .or. last > n .or. last < first - 1) return
        if (present(subset_by_value)) then
            if (.not. valid_value_interval(subset_by_value)) return
        else
            if (size(w) /= last - first + 1) return
        end if
        if (present(z)) then
            if (size(z, 1) /= n) return
            if (.not. present(subset_by_value)) then
                if (size(z, 2) /= size(w)) return
            end if
        end if
        if (present(max_iterations)) then
            if (max_iterations < 1) return
        end if

        use_lower = .true.
        if (present(lower)) use_lower = lower
        allocate(aherm(n, n), h(n, n), d(n), e(max(0, n - 1)), allw(n))
        aherm = 0.0_dp
        status = linalg_status_no_convergence
        do j = 1, n
            if (.not. ieee_is_finite(a(j, j))) return
            aherm(j, j) = a(j, j)
            do i = j + 1, n
                if (use_lower) then
                    value = a(i, j)
                else
                    value = a(j, i)
                end if
                if (.not. ieee_is_finite(value)) return
                aherm(i, j) = value
                aherm(j, i) = value
            end do
        end do

        if (present(z)) then
            allocate(q(n, n), u(n, n))
            call hessenberg(aherm, h, status, q)
        else
            call hessenberg(aherm, h, status)
        end if
        if (status /= linalg_status_success) return
        do j = 1, n
            d(j) = h(j, j)
        end do
        do j = 1, n - 1
            e(j) = 0.5_dp * h(j + 1, j) + 0.5_dp * h(j, j + 1)
        end do
        if (present(z)) then
            call eigh_tridiagonal(d, e, allw, status, u, max_iterations)
            if (status /= linalg_status_success) return
            if (present(subset_by_value)) then
                call value_selection_bounds(allw, subset_by_value, first, last)
            end if
            m = last - first + 1
            if (present(num_selected)) num_selected = m
            if (size(w) < m .or. size(z, 2) < m) then
                status = linalg_status_invalid_shape
                return
            end if
            z(:, 1:m) = matmul(q, u(:, first:last))
        else
            call eigh_tridiagonal(d, e, allw, status, max_iterations=max_iterations)
            if (status /= linalg_status_success) return
            if (present(subset_by_value)) then
                call value_selection_bounds(allw, subset_by_value, first, last)
            end if
            m = last - first + 1
            if (present(num_selected)) num_selected = m
            if (size(w) < m) then
                status = linalg_status_invalid_shape
                return
            end if
        end if
        w(1:m) = allw(first:last)
    end subroutine eigh_real

    recursive subroutine eigh_complex(a, w, status, z, lower, subset_by_index, &
        max_iterations, b, problem_type, subset_by_value, num_selected)
        complex(dp), intent(in) :: a(:, :) !! square Hermitian matrix; unused triangle ignored
        real(dp), intent(out) :: w(:) !! ascending eigenvalues, capacity in value mode
        integer, intent(out) :: status !! success, invalid shape, or numerical failure
        complex(dp), intent(out), optional :: z(:, :) !! unitary eigenvectors, capacity in value mode
        logical, intent(in), optional :: lower !! use lower input triangle; default true
        integer, intent(in), optional :: subset_by_index(2) !! inclusive 1-based [first,last]
        integer, intent(in), optional :: max_iterations !! QL iteration limit; default 120
        complex(dp), intent(in), optional :: b(:, :) !! positive-definite mass matrix, default I
        integer, intent(in), optional :: problem_type !! generalized type 1, 2, or 3; default 1
        real(dp), intent(in), optional :: subset_by_value(2) !! (low, high] eigenvalue interval
        integer, intent(out), optional :: num_selected !! number returned; may exceed output capacity

        complex(dp), allocatable :: aherm(:, :), h(:, :), q(:, :), phases(:)
        real(dp), allocatable :: d(:), e(:), allw(:), u(:, :)
        complex(dp) :: value, off
        real(dp) :: scalar, mag
        integer :: n, first, last, m, i, j
        logical :: use_lower

        w = 0.0_dp
        if (present(z)) z = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        if (present(num_selected)) num_selected = 0
        if (present(b)) then
            call generalized_complex(a, b, w, status, z, lower, &
                subset_by_index, max_iterations, problem_type, subset_by_value, num_selected)
            return
        end if
        n = size(a, 1)
        if (size(a, 2) /= n) return
        if (present(subset_by_index) .and. present(subset_by_value)) return
        call selection_bounds(n, subset_by_index, first, last)
        if (first < 1 .or. last > n .or. last < first - 1) return
        if (present(subset_by_value)) then
            if (.not. valid_value_interval(subset_by_value)) return
        else
            if (size(w) /= last - first + 1) return
        end if
        if (present(z)) then
            if (size(z, 1) /= n) return
            if (.not. present(subset_by_value)) then
                if (size(z, 2) /= size(w)) return
            end if
        end if
        if (present(max_iterations)) then
            if (max_iterations < 1) return
        end if

        use_lower = .true.
        if (present(lower)) use_lower = lower
        allocate(aherm(n, n), h(n, n), d(n), e(max(0, n - 1)), allw(n), phases(n))
        aherm = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_no_convergence
        do j = 1, n
            if (.not. ieee_is_finite(real(a(j, j), dp))) return
            ! As with LAPACK's Hermitian drivers, ignore imaginary diagonal parts.
            aherm(j, j) = cmplx(real(a(j, j), dp), 0.0_dp, dp)
            do i = j + 1, n
                if (use_lower) then
                    value = a(i, j)
                else
                    value = a(j, i)
                end if
                if (.not. finite_complex(value)) return
                if (use_lower) then
                    aherm(i, j) = value
                    aherm(j, i) = conjg(value)
                else
                    aherm(i, j) = conjg(value)
                    aherm(j, i) = value
                end if
            end do
        end do

        if (present(z)) then
            allocate(q(n, n), u(n, n))
            call hessenberg(aherm, h, status, q)
        else
            call hessenberg(aherm, h, status)
        end if
        if (status /= linalg_status_success) return
        do j = 1, n
            d(j) = real(h(j, j), dp)
        end do
        if (n >= 1) phases(1) = cmplx(1.0_dp, 0.0_dp, dp)
        do j = 1, n - 1
            ! D^H H D is real tridiagonal for a suitable unitary diagonal D.
            off = 0.5_dp * h(j + 1, j) + 0.5_dp * conjg(h(j, j + 1))
            scalar = max(abs(real(off, dp)), abs(aimag(off)))
            if (scalar == 0.0_dp) then
                e(j) = 0.0_dp
                phases(j + 1) = phases(j)
            else
                mag = hypot(real(off, dp) / scalar, aimag(off) / scalar)
                e(j) = scalar * mag
                if (.not. ieee_is_finite(e(j))) return
                phases(j + 1) = phases(j) * (off / scalar) / mag
            end if
        end do
        if (present(z)) then
            call eigh_tridiagonal(d, e, allw, status, u, max_iterations)
            if (status /= linalg_status_success) return
            if (present(subset_by_value)) then
                call value_selection_bounds(allw, subset_by_value, first, last)
            end if
            m = last - first + 1
            if (present(num_selected)) num_selected = m
            if (size(w) < m .or. size(z, 2) < m) then
                status = linalg_status_invalid_shape
                return
            end if
            do j = 1, n
                q(:, j) = q(:, j) * phases(j)
            end do
            z(:, 1:m) = matmul(q, cmplx(u(:, first:last), 0.0_dp, dp))
        else
            call eigh_tridiagonal(d, e, allw, status, max_iterations=max_iterations)
            if (status /= linalg_status_success) return
            if (present(subset_by_value)) then
                call value_selection_bounds(allw, subset_by_value, first, last)
            end if
            m = last - first + 1
            if (present(num_selected)) num_selected = m
            if (size(w) < m) then
                status = linalg_status_invalid_shape
                return
            end if
        end if
        w(1:m) = allw(first:last)
    end subroutine eigh_complex

    subroutine eigvalsh_real(a, w, status, lower, subset_by_index, &
        max_iterations, b, problem_type, subset_by_value, num_selected)
        real(dp), intent(in) :: a(:, :) !! symmetric matrix, size (n,n)
        real(dp), intent(out) :: w(:) !! eigenvalues of selected spectrum
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! use lower triangle; default true
        integer, intent(in), optional :: subset_by_index(2) !! inclusive 1-based [first,last]
        integer, intent(in), optional :: max_iterations !! QL iteration limit; default 120
        real(dp), intent(in), optional :: b(:, :) !! positive-definite mass matrix, default I
        integer, intent(in), optional :: problem_type !! generalized type 1, 2, or 3; default 1
        real(dp), intent(in), optional :: subset_by_value(2) !! (low, high] eigenvalue interval
        integer, intent(out), optional :: num_selected !! number returned; may exceed output capacity
        call eigh_real(a, w, status, lower=lower, subset_by_index=subset_by_index, &
            max_iterations=max_iterations, b=b, problem_type=problem_type, &
            subset_by_value=subset_by_value, num_selected=num_selected)
    end subroutine eigvalsh_real

    subroutine eigvalsh_complex(a, w, status, lower, subset_by_index, &
        max_iterations, b, problem_type, subset_by_value, num_selected)
        complex(dp), intent(in) :: a(:, :) !! Hermitian matrix, size (n,n)
        real(dp), intent(out) :: w(:) !! eigenvalues of selected spectrum
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: lower !! use lower triangle; default true
        integer, intent(in), optional :: subset_by_index(2) !! inclusive 1-based [first,last]
        integer, intent(in), optional :: max_iterations !! QL iteration limit; default 120
        complex(dp), intent(in), optional :: b(:, :) !! positive-definite mass matrix, default I
        integer, intent(in), optional :: problem_type !! generalized type 1, 2, or 3; default 1
        real(dp), intent(in), optional :: subset_by_value(2) !! (low, high] eigenvalue interval
        integer, intent(out), optional :: num_selected !! number returned; may exceed output capacity
        call eigh_complex(a, w, status, lower=lower, subset_by_index=subset_by_index, &
            max_iterations=max_iterations, b=b, problem_type=problem_type, &
            subset_by_value=subset_by_value, num_selected=num_selected)
    end subroutine eigvalsh_complex


    ! Cholesky reduction of a real symmetric definite matrix pencil.
    ! Type 1: A*z = lambda*B*z; type 2: A*B*z = lambda*z;
    ! type 3: B*A*z = lambda*z. B = L*transpose(L) in every case.
    subroutine generalized_real(a, b, w, status, z, lower, subset, &
        max_iterations, problem_type, subset_by_value, num_selected)
        real(dp), intent(in) :: a(:, :) !! symmetric left-hand coefficient, shape (n,n)
        real(dp), intent(in) :: b(:, :) !! positive-definite symmetric coefficient, shape (n,n)
        real(dp), intent(out) :: w(:) !! ascending eigenvalues; capacity in value mode
        integer, intent(out) :: status !! success, invalid shape, non-positive B, numerical failure
        real(dp), intent(out), optional :: z(:, :) !! generalized vectors; capacity in value mode
        logical, intent(in), optional :: lower !! triangle to read in A and B; default true
        integer, intent(in), optional :: subset(2) !! inclusive one-based interval
        integer, intent(in), optional :: max_iterations !! tridiagonal QL iteration limit
        integer, intent(in), optional :: problem_type !! eigenproblem type 1, 2 or 3; default 1
        real(dp), intent(in), optional :: subset_by_value(2) !! (low, high] eigenvalue interval
        integer, intent(out), optional :: num_selected !! number returned, including on capacity errors
        real(dp), allocatable :: am(:, :), bm(:, :), l(:, :), tmp(:, :), c(:, :), q(:, :)
        real(dp) :: av, bv
        integer :: n, first, last, m, selected_count, kind, i, j
        logical :: use_lower

        w = 0.0_dp
        if (present(z)) z = 0.0_dp
        status = linalg_status_invalid_shape
        if (present(num_selected)) num_selected = 0
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) return
        if (present(subset) .and. present(subset_by_value)) return
        call selection_bounds(n, subset, first, last)
        if (first < 1 .or. last > n .or. last < first - 1) return
        if (present(subset_by_value)) then
            if (.not. valid_value_interval(subset_by_value)) return
        else
            if (size(w) /= last - first + 1) return
        end if
        if (present(z)) then
            if (size(z, 1) /= n) return
            if (.not. present(subset_by_value)) then
                if (size(z, 2) /= size(w)) return
            end if
        end if
        kind = 1
        if (present(problem_type)) kind = problem_type
        if (kind < 1 .or. kind > 3) return
        if (present(max_iterations)) then
            if (max_iterations < 1) return
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        use_lower = .true.
        if (present(lower)) use_lower = lower
        allocate(am(n,n), bm(n,n), l(n,n), tmp(n,n), c(n,n))
        am = 0.0_dp
        bm = 0.0_dp
        status = linalg_status_no_convergence
        do j = 1, n
            av = a(j,j)
            bv = b(j,j)
            if (.not. ieee_is_finite(av) .or. .not. ieee_is_finite(bv)) return
            am(j,j) = av
            bm(j,j) = bv
            do i = j+1, n
                if (use_lower) then
                    av = a(i,j)
                    bv = b(i,j)
                else
                    av = a(j,i)
                    bv = b(j,i)
                end if
                if (.not. ieee_is_finite(av) .or. .not. ieee_is_finite(bv)) return
                am(i,j) = av
                am(j,i) = av
                bm(i,j) = bv
                bm(j,i) = bv
            end do
        end do
        call cholesky(bm,l,status,lower=.true.)
        if (status /= linalg_status_success) return
        if (kind == 1) then
            ! C = L**(-1) A L**(-T), via two triangular solves.
            call solve_triangular(l,am,tmp,status,lower=.true.)
            if (status /= linalg_status_success) return
            call solve_triangular(l,transpose(tmp),c,status,lower=.true.)
            if (status /= linalg_status_success) return
            c = transpose(c)
        else
            ! The type-2 and type-3 problems share C = L**T A L.
            c = matmul(transpose(l),matmul(am,l))
        end if
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(c))) return
        c = 0.5_dp*c + 0.5_dp*transpose(c)
        if (present(z)) then
            allocate(q(n,size(w)))
            call eigh_real(c,w,status,q,subset_by_index=subset, &
                max_iterations=max_iterations, subset_by_value=subset_by_value, &
                num_selected=selected_count)
            if (present(num_selected)) num_selected = selected_count
            if (status /= linalg_status_success) return
            m = selected_count
            if (size(z, 2) < m) then
                status = linalg_status_invalid_shape
                return
            end if
            if (kind == 3) then
                z(:, 1:m) = matmul(l,q(:, 1:m))
            else
                call solve_triangular(l,q(:, 1:m),z(:, 1:m),status,lower=.true.,trans='T')
                if (status /= linalg_status_success) return
            end if
            if (.not. all(ieee_is_finite(z))) then
                z = 0.0_dp
                status = linalg_status_no_convergence
            end if
        else
            call eigh_real(c,w,status,subset_by_index=subset, &
                max_iterations=max_iterations, subset_by_value=subset_by_value, &
                num_selected=selected_count)
            if (present(num_selected)) num_selected = selected_count
        end if
    end subroutine generalized_real

    ! For Hermitian definite pencils use B=L*L**H and conjugate-transpose
    ! triangular operations. Imaginary parts of diagonal entries are ignored.
    subroutine generalized_complex(a, b, w, status, z, lower, subset, &
        max_iterations, problem_type, subset_by_value, num_selected)
        complex(dp), intent(in) :: a(:, :) !! Hermitian left-hand coefficient, shape (n,n)
        complex(dp), intent(in) :: b(:, :) !! positive-definite Hermitian coefficient, shape (n,n)
        real(dp), intent(out) :: w(:) !! ascending eigenvalues; capacity in value mode
        integer, intent(out) :: status !! success, invalid shape, non-positive B, numerical failure
        complex(dp), intent(out), optional :: z(:, :) !! generalized vectors; capacity in value mode
        logical, intent(in), optional :: lower !! triangle to read in A and B; default true
        integer, intent(in), optional :: subset(2) !! inclusive one-based interval
        integer, intent(in), optional :: max_iterations !! tridiagonal QL iteration limit
        integer, intent(in), optional :: problem_type !! eigenproblem type 1, 2 or 3; default 1
        real(dp), intent(in), optional :: subset_by_value(2) !! (low, high] eigenvalue interval
        integer, intent(out), optional :: num_selected !! number returned, including on capacity errors
        complex(dp), allocatable :: am(:, :), bm(:, :), l(:, :), tmp(:, :), c(:, :), q(:, :)
        complex(dp) :: av, bv
        integer :: n, first, last, m, selected_count, kind, i, j
        logical :: use_lower

        w = 0.0_dp
        if (present(z)) z = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (present(num_selected)) num_selected = 0
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) return
        if (present(subset) .and. present(subset_by_value)) return
        call selection_bounds(n, subset, first, last)
        if (first < 1 .or. last > n .or. last < first - 1) return
        if (present(subset_by_value)) then
            if (.not. valid_value_interval(subset_by_value)) return
        else
            if (size(w) /= last - first + 1) return
        end if
        if (present(z)) then
            if (size(z, 1) /= n) return
            if (.not. present(subset_by_value)) then
                if (size(z, 2) /= size(w)) return
            end if
        end if
        kind = 1
        if (present(problem_type)) kind = problem_type
        if (kind < 1 .or. kind > 3) return
        if (present(max_iterations)) then
            if (max_iterations < 1) return
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        use_lower = .true.
        if (present(lower)) use_lower = lower
        allocate(am(n,n), bm(n,n), l(n,n), tmp(n,n), c(n,n))
        am = cmplx(0.0_dp,0.0_dp,dp)
        bm = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_no_convergence
        do j = 1, n
            if (.not. ieee_is_finite(real(a(j,j),dp))) return
            if (.not. ieee_is_finite(real(b(j,j),dp))) return
            am(j,j) = cmplx(real(a(j,j),dp),0.0_dp,dp)
            bm(j,j) = cmplx(real(b(j,j),dp),0.0_dp,dp)
            do i = j+1, n
                if (use_lower) then
                    av = a(i,j)
                    bv = b(i,j)
                    if (.not. finite_complex(av) .or. .not. finite_complex(bv)) return
                    am(i,j) = av
                    am(j,i) = conjg(av)
                    bm(i,j) = bv
                    bm(j,i) = conjg(bv)
                else
                    av = a(j,i)
                    bv = b(j,i)
                    if (.not. finite_complex(av) .or. .not. finite_complex(bv)) return
                    am(j,i) = av
                    am(i,j) = conjg(av)
                    bm(j,i) = bv
                    bm(i,j) = conjg(bv)
                end if
            end do
        end do
        call cholesky(bm,l,status,lower=.true.)
        if (status /= linalg_status_success) return
        if (kind == 1) then
            ! C = L**(-1) A L**(-H).
            call solve_triangular(l,am,tmp,status,lower=.true.)
            if (status /= linalg_status_success) return
            call solve_triangular(l,conjg(transpose(tmp)),c,status,lower=.true.)
            if (status /= linalg_status_success) return
            c = conjg(transpose(c))
        else
            ! The type-2 and type-3 problems share C = L**H A L.
            c = matmul(conjg(transpose(l)),matmul(am,l))
        end if
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(c,dp)))) return
        if (.not. all(ieee_is_finite(aimag(c)))) return
        c = 0.5_dp*c + 0.5_dp*conjg(transpose(c))
        if (present(z)) then
            allocate(q(n,size(w)))
            call eigh_complex(c,w,status,q,subset_by_index=subset, &
                max_iterations=max_iterations, subset_by_value=subset_by_value, &
                num_selected=selected_count)
            if (present(num_selected)) num_selected = selected_count
            if (status /= linalg_status_success) return
            m = selected_count
            if (size(z, 2) < m) then
                status = linalg_status_invalid_shape
                return
            end if
            if (kind == 3) then
                z(:, 1:m) = matmul(l,q(:, 1:m))
            else
                call solve_triangular(l,q(:, 1:m),z(:, 1:m),status,lower=.true.,trans='C')
                if (status /= linalg_status_success) return
            end if
            if (.not. all(ieee_is_finite(real(z,dp))) .or. &
                .not. all(ieee_is_finite(aimag(z)))) then
                z = cmplx(0.0_dp,0.0_dp,dp)
                status = linalg_status_no_convergence
            end if
        else
            call eigh_complex(c,w,status,subset_by_index=subset, &
                max_iterations=max_iterations, subset_by_value=subset_by_value, &
                num_selected=selected_count)
            if (present(num_selected)) num_selected = selected_count
        end if
    end subroutine generalized_complex

    pure subroutine selection_bounds(n, subset, first, last)
        integer, intent(in) :: n !! matrix order, >=0
        integer, intent(in), optional :: subset(2) !! optional inclusive 1-based range
        integer, intent(out) :: first !! lowest requested 1-based index
        integer, intent(out) :: last !! highest requested 1-based index
        first = 1
        last = n
        if (present(subset)) then
            first = subset(1)
            last = subset(2)
            ! An explicitly requested empty interval is invalid.
            if (first > last) first = 0
        end if
    end subroutine selection_bounds

    pure logical function valid_value_interval(bounds) result(ok)
        use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
        real(dp), intent(in) :: bounds(2) !! lower-exclusive and upper-inclusive eigenvalue bounds
        ! IEEE infinities are allowed; NaNs, reversed, and equal bounds are not.
        ok = .false.
        if (any(ieee_is_nan(bounds))) return
        ok = bounds(1) < bounds(2)
    end function valid_value_interval

    pure subroutine value_selection_bounds(values, bounds, first, last)
        real(dp), intent(in) :: values(:) !! sorted ascending eigenvalues
        real(dp), intent(in) :: bounds(2) !! half-open (low, high] interval
        integer, intent(out) :: first !! 1-based first matching eigenvalue
        integer, intent(out) :: last !! 1-based last matching eigenvalue; first-1 if none
        integer :: i
        first = size(values) + 1
        last = size(values)
        do i = 1, size(values)
            if (values(i) > bounds(1)) then
                first = i
                exit
            end if
        end do
        last = first - 1
        do i = first, size(values)
            if (values(i) > bounds(2)) exit
            last = i
        end do
    end subroutine value_selection_bounds

    pure logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! scalar to validate
        ok = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

end module scifort_linalg_dense_eigen
