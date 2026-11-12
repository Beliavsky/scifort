! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Quadratic assignment: independent Frank-Wolfe (FAQ) relaxation and
! best-improvement 2-opt permutation search. See CODE_PROVENANCE.md.
module scifort_optimize_quadratic_assignment
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_assignment, only : linear_sum_assignment, assignment_status_success
    implicit none
    private

    integer, parameter, public :: qap_success = 0
    integer, parameter, public :: qap_invalid_input = 1
    integer, parameter, public :: qap_assignment_failure = 2
    integer, parameter, public :: qap_numerical_failure = 3

    type, public :: qap_result
        integer, allocatable :: col_ind(:) !! one-based mapping of each A row to a B row
        real(dp) :: fun = 0.0_dp !! objective value for col_ind
        integer :: nit = 0 !! accepted swaps (2opt) or FW steps (FAQ)
        integer :: status = qap_invalid_input !! one of the qap_* constants
        logical :: success = .false. !! true when a valid permutation was returned
    end type qap_result

    public :: quadratic_assignment, quadratic_assignment_score

contains

    pure function quadratic_assignment(a, b, method, maximize, partial_match, &
            partial_guess, p0, maxiter, tol) result(res)
        real(dp), intent(in) :: a(:, :) !! square, finite matrix of interaction weights
        real(dp), intent(in) :: b(:, :) !! equally sized, square, finite matrix
        character(len=*), intent(in), optional :: method !! 'faq' (default) or '2opt'
        logical, intent(in), optional :: maximize !! maximize objective (default false)
        integer, intent(in), optional :: partial_match(:, :) !! fixed one-based (A,B) pairs, (m,2)
        integer, intent(in), optional :: partial_guess(:, :) !! initial pairs for 2opt, shape (g,2)
        real(dp), intent(in), optional :: p0(:, :) !! initial doubly stochastic free block for FAQ
        integer, intent(in), optional :: maxiter !! FAQ steps=30 or 2opt swaps=10000, by default
        real(dp), intent(in), optional :: tol !! FAQ movement tolerance (default 0.03)
        type(qap_result) :: res !! permutation, objective, iterations and status

        character(len=:), allocatable :: chosen
        integer :: i, j, k, n, nf, nseed, row, col, status, limit
        integer :: best_i, best_j, pass, ri, rj
        integer, allocatable :: fixed(:), perm(:), free_row(:), free_col(:)
        integer, allocatable :: guess(:), owner(:), rr(:), cc(:)
        logical :: want_maximum
        real(dp) :: threshold, score, trial_score, delta, best_delta
        real(dp) :: sign_objective, lin, quad, step, change, tiny_step
        real(dp), allocatable :: x(:, :), gradient(:, :), direction(:, :), weights(:, :)

        n = size(a, 1)
        allocate(res%col_ind(0))
        res%fun = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(a, 2) /= n) return
        if (size(b, 1) /= n .or. size(b, 2) /= n) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return

        chosen = 'faq'
        if (present(method)) chosen = lowercase(trim(adjustl(method)))
        if (chosen /= 'faq' .and. chosen /= '2opt') return
        if (present(partial_guess) .and. chosen /= '2opt') return
        if (present(p0) .and. chosen /= 'faq') return

        want_maximum = .false.
        if (present(maximize)) want_maximum = maximize
        sign_objective = 1.0_dp
        if (want_maximum) sign_objective = -1.0_dp
        limit = 30
        if (chosen == '2opt') limit = 10000
        if (present(maxiter)) limit = maxiter
        if (limit < 0) return
        threshold = 0.03_dp
        if (present(tol)) threshold = tol
        if (.not. ieee_is_finite(threshold)) return
        if (threshold < 0.0_dp) return
        if (present(partial_match)) then
            if (size(partial_match, 2) /= 2) return
            if (size(partial_match, 1) > n) return
        end if
        if (present(partial_guess)) then
            if (size(partial_guess, 2) /= 2) return
            if (size(partial_guess, 1) > n) return
        end if

        allocate(fixed(n), owner(n), guess(n), perm(n))
        fixed = 0
        owner = 0
        guess = 0
        if (present(partial_match)) then
            do k = 1, size(partial_match, 1)
                row = partial_match(k, 1)
                col = partial_match(k, 2)
                if (row < 1 .or. row > n .or. col < 1 .or. col > n) return
                if (fixed(row) /= 0 .or. owner(col) /= 0) return
                fixed(row) = col
                owner(col) = row
            end do
        end if
        nseed = count(fixed /= 0)
        nf = n - nseed
        allocate(free_row(nf), free_col(nf))
        j = 0
        k = 0
        do i = 1, n
            if (fixed(i) == 0) then
                j = j + 1
                free_row(j) = i
            end if
            if (owner(i) == 0) then
                k = k + 1
                free_col(k) = i
            end if
        end do

        if (present(partial_guess)) then
            do k = 1, size(partial_guess, 1)
                row = partial_guess(k, 1)
                col = partial_guess(k, 2)
                if (row < 1 .or. row > n .or. col < 1 .or. col > n) return
                if (fixed(row) /= 0) then
                    if (fixed(row) /= col) return
                else
                    if (owner(col) /= 0 .or. guess(row) /= 0) return
                    if (any(guess == col)) return
                    guess(row) = col
                end if
            end do
        end if

        if (present(p0)) then
            if (size(p0, 1) /= nf .or. size(p0, 2) /= nf) return
            if (.not. all(ieee_is_finite(p0))) return
            if (nf > 0) then
                if (minval(p0) < -1.0e-10_dp .or. maxval(p0) > 1.0_dp + 1.0e-10_dp) return
                do i = 1, nf
                    if (abs(sum(p0(i, :)) - 1.0_dp) > 1.0e-8_dp) return
                    if (abs(sum(p0(:, i)) - 1.0_dp) > 1.0e-8_dp) return
                end do
            end if
        end if

        perm = fixed
        if (nf == 0) then
            score = quadratic_assignment_score(a, b, perm)
            if (.not. ieee_is_finite(score)) then
                res%status = qap_numerical_failure
                return
            end if
            call finish(perm, score, res)
            return
        end if

        select case (chosen)
        case ('2opt')
            ! Apply optional partial guesses; fill all remaining slots in index order.
            do i = 1, n
                if (guess(i) /= 0) perm(i) = guess(i)
            end do
            j = 1
            do k = 1, nf
                row = free_row(k)
                if (perm(row) /= 0) cycle
                do while (j <= nf)
                    col = free_col(j)
                    j = j + 1
                    if (any(perm == col)) cycle
                    perm(row) = col
                    exit
                end do
            end do
            if (any(perm == 0)) return
            score = quadratic_assignment_score(a, b, perm)
            if (.not. ieee_is_finite(score)) then
                res%status = qap_numerical_failure
                return
            end if
            tiny_step = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(score))
            do pass = 1, limit
                best_delta = -tiny_step
                best_i = 0
                best_j = 0
                do i = 1, nf - 1
                    do j = i + 1, nf
                        ri = free_row(i)
                        rj = free_row(j)
                        delta = swap_delta(a, b, perm, ri, rj)
                        if (.not. ieee_is_finite(delta)) then
                            res%status = qap_numerical_failure
                            return
                        end if
                        if (sign_objective * delta < best_delta) then
                            best_delta = sign_objective * delta
                            best_i = ri
                            best_j = rj
                        end if
                    end do
                end do
                if (best_i == 0) exit
                col = perm(best_i)
                perm(best_i) = perm(best_j)
                perm(best_j) = col
                trial_score = quadratic_assignment_score(a, b, perm)
                if (.not. ieee_is_finite(trial_score)) then
                    res%status = qap_numerical_failure
                    return
                end if
                ! Re-evaluate any accepted swap to prevent accumulated delta drift.
                if (sign_objective * (trial_score - score) >= -tiny_step) then
                    col = perm(best_i)
                    perm(best_i) = perm(best_j)
                    perm(best_j) = col
                    exit
                end if
                score = trial_score
                res%nit = res%nit + 1
            end do
        case ('faq')
            allocate(x(n, n), gradient(n, n), direction(n, n), weights(nf, nf))
            x = 0.0_dp
            do i = 1, n
                if (fixed(i) /= 0) x(i, fixed(i)) = 1.0_dp
            end do
            do j = 1, nf
                do i = 1, nf
                    x(free_row(i), free_col(j)) = 1.0_dp / real(nf, dp)
                    if (present(p0)) x(free_row(i), free_col(j)) = p0(i, j)
                end do
            end do
            do pass = 1, limit
                gradient = matmul(a, matmul(x, transpose(b))) + &
                    matmul(transpose(a), matmul(x, b))
                if (.not. all(ieee_is_finite(gradient))) then
                    res%status = qap_numerical_failure
                    return
                end if
                do j = 1, nf
                    do i = 1, nf
                        weights(i, j) = sign_objective * gradient(free_row(i), free_col(j))
                    end do
                end do
                call linear_sum_assignment(weights, rr, cc, status)
                if (status /= assignment_status_success) then
                    res%status = qap_assignment_failure
                    return
                end if
                direction = 0.0_dp
                do j = 1, nf
                    do i = 1, nf
                        direction(free_row(i), free_col(j)) = -x(free_row(i), free_col(j))
                    end do
                end do
                do k = 1, nf
                    direction(free_row(rr(k)), free_col(cc(k))) = &
                        direction(free_row(rr(k)), free_col(cc(k))) + 1.0_dp
                end do
                lin = sum(gradient * direction)
                quad = dense_objective(a, b, direction)
                if (.not. ieee_is_finite(lin) .or. .not. ieee_is_finite(quad)) then
                    res%status = qap_numerical_failure
                    return
                end if
                if (sign_objective * quad > 0.0_dp) then
                    step = max(0.0_dp, min(1.0_dp, -lin / (2.0_dp * quad)))
                else if (sign_objective * (lin + quad) < 0.0_dp) then
                    step = 1.0_dp
                else
                    step = 0.0_dp
                end if
                if (step <= 0.0_dp) exit
                change = step * sqrt(sum(direction * direction))
                x = x + step * direction
                res%nit = res%nit + 1
                if (change <= threshold * sqrt(real(nf, dp))) exit
            end do
            do j = 1, nf
                do i = 1, nf
                    weights(i, j) = -x(free_row(i), free_col(j))
                end do
            end do
            call linear_sum_assignment(weights, rr, cc, status)
            if (status /= assignment_status_success) then
                res%status = qap_assignment_failure
                return
            end if
            do k = 1, nf
                perm(free_row(rr(k))) = free_col(cc(k))
            end do
            score = quadratic_assignment_score(a, b, perm)
            if (.not. ieee_is_finite(score)) then
                res%status = qap_numerical_failure
                return
            end if
        end select
        call finish(perm, score, res)
    end function quadratic_assignment

    pure function quadratic_assignment_score(a, b, col_ind) result(value)
        real(dp), intent(in) :: a(:, :) !! square interaction weights
        real(dp), intent(in) :: b(:, :) !! square interaction weights, same dimension
        integer, intent(in) :: col_ind(:) !! complete one-based permutation
        real(dp) :: value !! sum of A(i,j)*B(col_ind(i),col_ind(j)), or NaN if invalid
        integer :: i, j, n
        n = size(col_ind)
        value = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(a, 1) /= n .or. size(a, 2) /= n) return
        if (size(b, 1) /= n .or. size(b, 2) /= n) return
        if (any(col_ind < 1) .or. any(col_ind > n)) return
        do i = 1, n
            if (count(col_ind == col_ind(i)) /= 1) return
        end do
        value = 0.0_dp
        do j = 1, n
            do i = 1, n
                value = value + a(i, j) * b(col_ind(i), col_ind(j))
            end do
        end do
    end function quadratic_assignment_score

    pure function dense_objective(a, b, x) result(value)
        real(dp), intent(in) :: a(:, :) !! square first interaction matrix
        real(dp), intent(in) :: b(:, :) !! square second interaction matrix
        real(dp), intent(in) :: x(:, :) !! relaxed assignment or search direction
        real(dp) :: value !! relaxed objective of x
        value = sum(a * matmul(matmul(x, b), transpose(x)))
    end function dense_objective

    pure function swap_delta(a, b, perm, r, s) result(value)
        real(dp), intent(in) :: a(:, :) !! first interaction matrix
        real(dp), intent(in) :: b(:, :) !! second interaction matrix
        integer, intent(in) :: perm(:) !! current one-based assignment
        integer, intent(in) :: r !! first row to swap
        integer, intent(in) :: s !! second row to swap
        real(dp) :: value !! exact objective change from swapping assignments at r and s
        integer :: k, u, v, w
        u = perm(r)
        v = perm(s)
        value = (a(r, r) - a(s, s)) * (b(v, v) - b(u, u)) + &
            (a(r, s) - a(s, r)) * (b(v, u) - b(u, v))
        do k = 1, size(perm)
            if (k == r .or. k == s) cycle
            w = perm(k)
            value = value + (a(r, k) - a(s, k)) * (b(v, w) - b(u, w)) + &
                (a(k, r) - a(k, s)) * (b(w, v) - b(w, u))
        end do
    end function swap_delta

    pure subroutine finish(perm, value, res)
        integer, intent(in) :: perm(:) !! validated complete permutation
        real(dp), intent(in) :: value !! evaluated objective value
        type(qap_result), intent(inout) :: res !! solver result to finalize
        res%col_ind = perm
        res%fun = value
        res%status = qap_success
        res%success = .true.
    end subroutine finish

    pure function lowercase(source) result(out)
        character(len=*), intent(in) :: source !! method string to normalize
        character(len=len(source)) :: out !! lower-case method string
        integer :: i, c
        out = source
        do i = 1, len(source)
            c = iachar(source(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
        end do
    end function lowercase

end module scifort_optimize_quadratic_assignment
