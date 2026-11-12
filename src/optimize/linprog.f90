! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dense two-phase simplex linear programming (independent implementation).
! Scope: continuous LPs with inequality/equality constraints and variable bounds.
! No external BLAS/LAPACK calls or third-party source are used.
module scifort_optimize_linprog
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: linprog_optimal = 0
    integer, parameter, public :: linprog_iteration_limit = 1
    integer, parameter, public :: linprog_infeasible = 2
    integer, parameter, public :: linprog_unbounded = 3
    integer, parameter, public :: linprog_numerical = 4
    integer, parameter, public :: linprog_invalid = -1

    type, public :: linprog_result
        real(dp), allocatable :: x(:) !! decision vector; allocated after valid setup
        real(dp), allocatable :: slack(:) !! b_ub - A_ub*x
        real(dp), allocatable :: con(:) !! b_eq - A_eq*x
        real(dp) :: fun = 0.0_dp !! objective c^T*x
        integer :: status = linprog_invalid !! SciPy-compatible: 0 optimal, 1 limit, 2 infeasible, 3 unbounded, 4 numerical
        integer :: nit = 0 !! total Phase I and II pivot steps
        logical :: success = .false. !! optimal feasible solution obtained
        character(len=180) :: message = 'invalid linear-programming input' !! termination explanation
    end type linprog_result

    public :: linprog

contains

    function linprog(c, a_ub, b_ub, a_eq, b_eq, lower, upper, tol, maxiter) result(res)
        real(dp), intent(in) :: c(:) !! length-n objective to minimize
        real(dp), intent(in), optional :: a_ub(:, :) !! m_ub-by-n inequality coefficients
        real(dp), intent(in), optional :: b_ub(:) !! length-m_ub upper right-hand sides
        real(dp), intent(in), optional :: a_eq(:, :) !! m_eq-by-n equality coefficients
        real(dp), intent(in), optional :: b_eq(:) !! length-m_eq equality right-hand sides
        real(dp), intent(in), optional :: lower(:) !! length-n lower bounds; default zero; -infinity permitted
        real(dp), intent(in), optional :: upper(:) !! length-n upper bounds; default +infinity
        real(dp), intent(in), optional :: tol !! simplex pivot and feasibility tolerance, default 1e-9
        integer, intent(in), optional :: maxiter !! maximum total simplex pivots across both phases
        type(linprog_result) :: res

        real(dp), allocatable :: lo(:), hi(:), trans(:, :), offset(:), ct(:), &
            row_a(:, :), row_b(:), tableau(:, :), cmax(:), y(:), coeff(:)
        integer, allocatable :: senses(:), basis(:), primary_col(:)
        logical, allocatable :: art(:)
        real(dp) :: inf, tolerance, objective_scale, phase1_rhs, feasibility_tol
        integer :: n, nub, neq, nbound, ndec, m, mtotal, ncols, nslack, nart
        integer :: i, j, k, p, r, info, limit, slack_col, art_col

        n = size(c)
        res%status = linprog_invalid
        if (.not. all(ieee_is_finite(c))) return
        if (present(a_ub) .neqv. present(b_ub)) return
        if (present(a_eq) .neqv. present(b_eq)) return
        nub = 0
        neq = 0
        if (present(a_ub)) then
            nub = size(a_ub, 1)
            if (size(a_ub, 2) /= n .or. size(b_ub) /= nub) return
            if (.not. all(ieee_is_finite(a_ub)) .or. .not. all(ieee_is_finite(b_ub))) return
        end if
        if (present(a_eq)) then
            neq = size(a_eq, 1)
            if (size(a_eq, 2) /= n .or. size(b_eq) /= neq) return
            if (.not. all(ieee_is_finite(a_eq)) .or. .not. all(ieee_is_finite(b_eq))) return
        end if
        if (present(lower)) then
            if (size(lower) /= n) return
            if (any(ieee_is_nan(lower))) return
        end if
        if (present(upper)) then
            if (size(upper) /= n) return
            if (any(ieee_is_nan(upper))) return
        end if
        tolerance = 1.0e-9_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0.0_dp .or. tolerance >= 1.0e-2_dp) return
        limit = max(500, 40*(n + nub + neq + 1))
        if (present(maxiter)) limit = maxiter
        if (limit < 0) return

        inf = ieee_value(1.0_dp, ieee_positive_inf)
        allocate(lo(n), hi(n), offset(n), primary_col(n))
        lo = 0.0_dp
        hi = inf
        if (present(lower)) lo = lower
        if (present(upper)) hi = upper
        if (any(lo > hi) .or. any(lo >= inf) .or. any(hi <= -inf)) then
            res%message = 'inconsistent or invalid variable bounds'
            return
        end if

        ! x = offset + trans*y, y >= 0, including unrestricted variables
        ! represented by the difference of two nonnegative variables.
        ndec = n + count(.not. ieee_is_finite(lo) .and. .not. ieee_is_finite(hi))
        allocate(trans(n, ndec), ct(ndec), coeff(ndec))
        trans = 0.0_dp
        offset = 0.0_dp
        p = 0
        nbound = 0
        do i = 1, n
            p = p + 1
            primary_col(i) = p
            if (ieee_is_finite(lo(i))) then
                offset(i) = lo(i)
                trans(i,p) = 1.0_dp
                if (ieee_is_finite(hi(i))) nbound = nbound + 1
            else if (ieee_is_finite(hi(i))) then
                offset(i) = hi(i)
                trans(i,p) = -1.0_dp
            else
                trans(i,p) = 1.0_dp
                p = p + 1
                trans(i,p) = -1.0_dp
            end if
        end do
        if (ndec > 0) ct = matmul(c, trans)
        mtotal = nub + neq + nbound
        allocate(row_a(mtotal,ndec), row_b(mtotal), senses(mtotal))
        row_a = 0.0_dp
        row_b = 0.0_dp
        senses = 0
        r = 0
        if (present(a_ub)) then
            do i = 1, nub
                coeff = matmul(a_ub(i,:),trans)
                call add_row(coeff, b_ub(i)-dot_product(a_ub(i,:),offset), -1)
            end do
        end if
        if (present(a_eq)) then
            do i = 1, neq
                coeff = matmul(a_eq(i,:),trans)
                call add_row(coeff, b_eq(i)-dot_product(a_eq(i,:),offset), 0)
            end do
        end if
        do i = 1,n
            if (.not. ieee_is_finite(lo(i)) .or. .not. ieee_is_finite(hi(i))) cycle
            coeff = 0.0_dp
            coeff(primary_col(i)) = 1.0_dp
            call add_row(coeff,hi(i)-lo(i),-1)
        end do
        if (.not. all(ieee_is_finite(row_a)) .or. .not. all(ieee_is_finite(row_b)) .or. &
            .not. all(ieee_is_finite(ct))) then
            res%status = linprog_numerical
            res%message = 'transformed LP has nonfinite coefficients'
            return
        end if

        nslack = count(senses /= 0)
        nart = count(senses /= -1)
        ncols = ndec + nslack + nart
        m = mtotal
        allocate(tableau(0:m, ncols+1), basis(m), art(ncols), cmax(ncols))
        tableau = 0.0_dp
        art = .false.
        slack_col = ndec
        art_col = ndec + nslack
        do i = 1,m
            if (ndec > 0) tableau(i,1:ndec) = row_a(i,:)
            tableau(i,ncols+1) = row_b(i)
            select case (senses(i))
            case (-1) ! a*y <= b: slack variable is basic.
                slack_col = slack_col + 1
                tableau(i,slack_col) = 1.0_dp
                basis(i) = slack_col
            case (1) ! a*y >= b: surplus and artificial variables.
                slack_col = slack_col + 1
                tableau(i,slack_col) = -1.0_dp
                art_col = art_col + 1
                tableau(i,art_col) = 1.0_dp
                basis(i) = art_col
                art(art_col) = .true.
            case default ! a*y = b: artificial variable.
                art_col = art_col + 1
                tableau(i,art_col) = 1.0_dp
                basis(i) = art_col
                art(art_col) = .true.
            end select
        end do

        if (nart > 0) then
            ! Phase I: maximize -sum(artificial variables).
            cmax = 0.0_dp
            where (art) cmax = -1.0_dp
            call reset_objective(tableau,basis,m,ncols,cmax)
            feasibility_tol = max(tolerance, 128.0_dp*epsilon(1.0_dp))
            if (m > 0) feasibility_tol = max(feasibility_tol, &
                128.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(row_b))))
            call simplex_iterate(tableau,basis,art,m,ncols,.true.,tolerance,limit,res%nit,info)
            if (info /= linprog_optimal) then
                res%status = linprog_numerical
                if (info == linprog_iteration_limit) res%status = info
                res%message = 'Phase I did not converge'
                return
            end if
            phase1_rhs = tableau(0,ncols+1)
            if (phase1_rhs < -feasibility_tol) then
                res%status = linprog_infeasible
                res%message = 'linear constraints are infeasible (Phase I)'
                return
            end if

            ! Eliminate zero-valued artificial basic variables. Drop linearly
            ! redundant rows if they contain no usable nonartificial pivot.
            i = 1
            do while (i <= m)
                if (.not. art(basis(i))) then
                    i = i + 1
                    cycle
                end if
                if (abs(tableau(i,ncols+1)) > feasibility_tol) then
                    res%status = linprog_infeasible
                    res%message = 'nonzero artificial variable after Phase I'
                    return
                end if
                j = 0
                do k=1,ncols
                    if (art(k)) cycle
                    if (abs(tableau(i,k)) > tolerance) then
                        j = k
                        exit
                    end if
                end do
                if (j /= 0) then
                    call pivot(tableau,basis,m,i,j)
                    i = i + 1
                else
                    if (i < m) then
                        tableau(i:m-1,:) = tableau(i+1:m,:)
                        basis(i:m-1) = basis(i+1:m)
                    end if
                    m = m - 1
                end if
            end do
        end if

        ! Phase II: minimize c^T*x = c^T*offset + ct^T*y.
        objective_scale = 1.0_dp
        if (ndec > 0) objective_scale = maxval(abs(ct))
        if (objective_scale <= 0.0_dp) objective_scale = 1.0_dp
        cmax = 0.0_dp
        if (ndec > 0) cmax(1:ndec) = -ct/objective_scale
        call reset_objective(tableau,basis,m,ncols,cmax)
        call simplex_iterate(tableau,basis,art,m,ncols,.false.,tolerance,limit,res%nit,info)
        res%status = info
        select case (info)
        case (linprog_optimal)
            res%message = 'optimal solution found by two-phase simplex'
            res%success = .true.
        case (linprog_unbounded)
            res%message = 'objective is unbounded below'
        case (linprog_iteration_limit)
            res%message = 'simplex iteration limit reached'
        case default
            res%message = 'numerical difficulties in simplex pivoting'
        end select
        allocate(y(ndec),res%x(n),res%slack(nub),res%con(neq))
        y = 0.0_dp
        do i=1,m
            if (basis(i) <= ndec) y(basis(i)) = max(0.0_dp,tableau(i,ncols+1))
        end do
        res%x = offset + matmul(trans,y)
        res%fun = dot_product(c,res%x)
        if (present(a_ub)) res%slack = b_ub - matmul(a_ub,res%x)
        if (present(a_eq)) res%con = b_eq - matmul(a_eq,res%x)
        if (.not. all(ieee_is_finite(res%x)) .or. .not. ieee_is_finite(res%fun)) then
            res%success = .false.
            res%status = linprog_numerical
            res%message = 'nonfinite returned solution'
        end if

    contains

        subroutine add_row(a, b, s)
            real(dp), intent(in) :: a(:) !! transformed constraint coefficients
            real(dp), intent(in) :: b !! constraint right-hand side
            integer, intent(in) :: s !! -1 <=, 0 =, 1 >=
            real(dp) :: scale
            r = r + 1
            ! Scale using coefficients, not the RHS: a constraint x <= 1e12
            ! still needs an O(1) pivot rather than a tiny one.
            scale = 1.0_dp
            if (size(a)>0) then
                scale = maxval(abs(a))
                if (scale <= tiny(1.0_dp)) scale = 1.0_dp
            end if
            if (b < 0.0_dp) then
                row_a(r,:) = -a/scale
                row_b(r) = -b/scale
                senses(r) = -s
            else
                row_a(r,:) = a/scale
                row_b(r) = b/scale
                senses(r) = s
            end if
        end subroutine add_row

    end function linprog

    subroutine reset_objective(tab,basis,m,ncols,cmax)
        real(dp), intent(inout) :: tab(0:,:) !! simplex tableau, objective at row 0
        integer, intent(in) :: basis(:) !! current basic column indices
        integer, intent(in) :: m !! active constraint row count
        integer, intent(in) :: ncols !! number of variable columns
        real(dp), intent(in) :: cmax(:) !! maximization objective coefficients
        integer :: i
        tab(0,:) = 0.0_dp
        tab(0,1:ncols) = -cmax
        do i=1,m
            tab(0,:) = tab(0,:) + cmax(basis(i))*tab(i,:)
        end do
    end subroutine reset_objective

    subroutine pivot(tab,basis,m,row,col)
        real(dp), intent(inout) :: tab(0:,:) !! tableau modified in-place
        integer, intent(inout) :: basis(:) !! basic indices modified in-place
        integer, intent(in) :: m !! active row count
        integer, intent(in) :: row !! pivot row from 1 to m
        integer, intent(in) :: col !! pivot variable column
        integer :: i
        real(dp) :: factor
        factor = tab(row,col)
        tab(row,:) = tab(row,:)/factor
        do i=0,m
            if (i==row) cycle
            factor = tab(i,col)
            tab(i,:) = tab(i,:) - factor*tab(row,:)
        end do
        basis(row) = col
    end subroutine pivot

    subroutine simplex_iterate(tab,basis,art,m,ncols,phase1,tol,maxiter,nit,status)
        real(dp), intent(inout) :: tab(0:,:) !! simplex tableau (row 0 objective)
        integer, intent(inout) :: basis(:) !! basic variable indices
        logical, intent(in) :: art(:) !! true for artificial columns
        integer, intent(in) :: m !! active rows
        integer, intent(in) :: ncols !! number of variable columns
        logical, intent(in) :: phase1 !! allow artificial columns to enter (Phase I)
        real(dp), intent(in) :: tol !! reduced-cost and pivot tolerance
        integer, intent(in) :: maxiter !! overall pivot budget
        integer, intent(inout) :: nit !! accumulated pivots
        integer, intent(out) :: status !! 0 optimum, 1 iteration limit, 3 unbounded, 4 numerical failure
        integer :: j, i, entering, leaving
        real(dp) :: ratio, smallest, ratio_tie_tol

        status = linprog_optimal
        do
            entering = 0
            ! Bland's smallest-index rule prevents cycling with degenerate LPs.
            do j=1,ncols
                if (.not. phase1 .and. art(j)) cycle
                if (tab(0,j) < -tol) then
                    entering = j
                    exit
                end if
            end do
            if (entering == 0) return
            if (nit >= maxiter) then
                status = linprog_iteration_limit
                return
            end if
            leaving = 0
            smallest = huge(1.0_dp)
            do i=1,m
                if (tab(i,entering) <= tol) cycle
                if (tab(i,ncols+1) < -tol) then
                    status = linprog_numerical
                    return
                end if
                ratio = max(0.0_dp,tab(i,ncols+1))/tab(i,entering)
                if (leaving == 0) then
                    leaving = i
                    smallest = ratio
                else
                    ! Do not use the user simplex tolerance scaled by the
                    ! *size of the ratio* to declare ties. Near 1e12 it can
                    ! treat RHS differences of order 1 as equal and pivot
                    ! onto a nonminimum row, breaking Phase I feasibility.
                    ratio_tie_tol = 16.0_dp*epsilon(1.0_dp)* &
                        max(1.0_dp,abs(ratio),abs(smallest))
                    if (ratio < smallest - ratio_tie_tol) then
                        leaving = i
                        smallest = ratio
                    else if (abs(ratio-smallest) <= ratio_tie_tol) then
                        if (basis(i) < basis(leaving)) then
                            leaving = i
                            smallest = ratio
                        end if
                    end if
                end if
            end do
            if (leaving == 0) then
                status = linprog_unbounded
                return
            end if
            call pivot(tab,basis,m,leaving,entering)
            nit = nit + 1
            if (.not. all(ieee_is_finite(tab(0:m,:)))) then
                status = linprog_numerical
                return
            end if
        end do
    end subroutine simplex_iterate
end module scifort_optimize_linprog
