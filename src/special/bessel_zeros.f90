! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Positive real zeros of integer-order cylindrical Bessel functions.  The
! implementation is original SciFort code based on the interlacing/simple-root
! properties of real Bessel functions and safeguarded Newton refinement.

module scifort_bessel_zeros
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_bessel_derivatives, only : jvp, yvp
    use scifort_bessel_jy, only : jv, yv
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: bessel_zeros_status_success = 0
    integer, parameter, public :: bessel_zeros_status_invalid_argument = 1
    integer, parameter, public :: bessel_zeros_status_numerical_failure = 2

    integer, parameter :: root_j = 1
    integer, parameter :: root_jp = 2
    integer, parameter :: root_y = 3
    integer, parameter :: root_yp = 4
    integer, parameter :: jnjnp_max_zeros = 1200

    public :: jn_zeros
    public :: jnjnp_zeros
    public :: jnp_zeros
    public :: jnyn_zeros
    public :: y0_zeros
    public :: y1_zeros
    public :: y1p_zeros
    public :: yn_zeros
    public :: ynp_zeros

contains

    pure subroutine jn_zeros(n, nt, zeros, status)
        integer, intent(in) :: n !! integer Bessel order; negative orders use abs(n)
        integer, intent(in) :: nt !! number of positive zeros requested, >= 1
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of J_n
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_zeros(n, nt, root_j, zeros, status)
    end subroutine jn_zeros

    pure subroutine jnjnp_zeros(nt, zeros, orders, serials, kinds, status)
        integer, intent(in) :: nt !! number of merged J/J-prime zeros requested, 1..1200
        real(dp), allocatable, intent(out) :: zeros(:) !! zeros ordered by increasing magnitude
        integer, allocatable, intent(out) :: orders(:) !! Bessel order associated with each zero
        integer, allocatable, intent(out) :: serials(:) !! one-based root index; J_0'(0) uses zero
        integer, allocatable, intent(out) :: kinds(:) !! zero for J_n roots; one for J_n-prime roots
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: capacity
        integer :: i
        integer :: info
        integer :: mmax
        integer :: n
        integer :: ncandidate
        integer :: nmax
        real(dp) :: cutoff
        real(dp) :: omitted_min
        real(dp), allocatable :: candidates(:)
        real(dp), allocatable :: roots(:)
        integer, allocatable :: candidate_kinds(:)
        integer, allocatable :: candidate_orders(:)
        integer, allocatable :: candidate_serials(:)

        if (nt < 1 .or. nt > jnjnp_max_zeros) then
            allocate(zeros(0), orders(0), serials(0), kinds(0))
            if (present(status)) status = bessel_zeros_status_invalid_argument
            return
        end if

        allocate(zeros(nt), orders(nt), serials(nt), kinds(nt))
        zeros = quiet_nan(0.0_dp)
        orders = -1
        serials = -1
        kinds = -1

        nmax = max(4, ceiling(1.5_dp * sqrt(real(nt, dp))))
        mmax = nmax

        do
            capacity = 1 + 2 * (nmax + 1) * mmax
            allocate(candidates(capacity))
            allocate(candidate_orders(capacity))
            allocate(candidate_serials(capacity))
            allocate(candidate_kinds(capacity))

            ncandidate = 1
            candidates(1) = 0.0_dp
            candidate_orders(1) = 0
            candidate_serials(1) = 0
            candidate_kinds(1) = 1
            omitted_min = huge(1.0_dp)

            do n = 0, nmax
                call compute_zeros(n, mmax + 1, root_j, roots, info)
                if (info /= bessel_zeros_status_success) then
                    if (present(status)) status = info
                    return
                end if
                do i = 1, mmax
                    ncandidate = ncandidate + 1
                    candidates(ncandidate) = roots(i)
                    candidate_orders(ncandidate) = n
                    candidate_serials(ncandidate) = i
                    candidate_kinds(ncandidate) = 0
                end do
                omitted_min = min(omitted_min, roots(mmax + 1))

                call compute_zeros(n, mmax + 1, root_jp, roots, info)
                if (info /= bessel_zeros_status_success) then
                    if (present(status)) status = info
                    return
                end if
                do i = 1, mmax
                    ncandidate = ncandidate + 1
                    candidates(ncandidate) = roots(i)
                    candidate_orders(ncandidate) = n
                    candidate_serials(ncandidate) = i
                    candidate_kinds(ncandidate) = 1
                end do
                omitted_min = min(omitted_min, roots(mmax + 1))
            end do

            call compute_zeros(nmax + 1, 1, root_j, roots, info)
            if (info /= bessel_zeros_status_success) then
                if (present(status)) status = info
                return
            end if
            omitted_min = min(omitted_min, roots(1))
            call compute_zeros(nmax + 1, 1, root_jp, roots, info)
            if (info /= bessel_zeros_status_success) then
                if (present(status)) status = info
                return
            end if
            omitted_min = min(omitted_min, roots(1))

            call coalesce_j0p_j1(candidates(:ncandidate), &
                candidate_orders(:ncandidate), candidate_serials(:ncandidate), &
                candidate_kinds(:ncandidate))
            call select_first_jnjnp(nt, candidates(:ncandidate), &
                candidate_orders(:ncandidate), candidate_serials(:ncandidate), &
                candidate_kinds(:ncandidate), zeros, orders, serials, kinds)
            cutoff = zeros(nt)
            if (cutoff + 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(cutoff)) < omitted_min) then
                if (present(status)) status = bessel_zeros_status_success
                return
            end if

            deallocate(candidates, candidate_orders, candidate_serials, candidate_kinds)
            nmax = nmax + max(4, nmax / 2)
            mmax = mmax + max(4, mmax / 2)
        end do
    end subroutine jnjnp_zeros

    pure subroutine coalesce_j0p_j1(candidates, candidate_orders, candidate_serials, &
            candidate_kinds)
        real(dp), intent(inout) :: candidates(:) !! roots with J_0-prime/J_1 ties coalesced
        integer, intent(in) :: candidate_orders(:) !! Bessel order for each candidate
        integer, intent(in) :: candidate_serials(:) !! serial root number for each candidate
        integer, intent(in) :: candidate_kinds(:) !! zero for J and one for J-prime

        integer :: i
        integer :: j
        real(dp) :: tie_tol

        do i = 1, size(candidates)
            if (candidate_orders(i) /= 0 .or. candidate_kinds(i) /= 1) cycle
            if (candidate_serials(i) < 1) cycle
            do j = 1, size(candidates)
                if (candidate_orders(j) /= 1 .or. candidate_kinds(j) /= 0) cycle
                if (candidate_serials(j) /= candidate_serials(i)) cycle
                tie_tol = 64.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(candidates(i)), abs(candidates(j)))
                if (abs(candidates(i) - candidates(j)) <= tie_tol) candidates(i) = candidates(j)
                exit
            end do
        end do
    end subroutine coalesce_j0p_j1

    pure subroutine jnp_zeros(n, nt, zeros, status)
        integer, intent(in) :: n !! integer Bessel order; negative orders use abs(n)
        integer, intent(in) :: nt !! number of positive zeros requested, >= 1
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of J_n'
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_zeros(n, nt, root_jp, zeros, status)
    end subroutine jnp_zeros

    pure subroutine y0_zeros(nt, zeros, derivatives, status)
        integer, intent(in) :: nt !! number of positive real zeros requested, >= 1
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of Y_0
        real(dp), allocatable, intent(out) :: derivatives(:) !! Y_0' evaluated at each zero
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_y_zero_pair(0, nt, root_y, .true., zeros, derivatives, status)
    end subroutine y0_zeros

    pure subroutine y1_zeros(nt, zeros, derivatives, status)
        integer, intent(in) :: nt !! number of positive real zeros requested, >= 1
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of Y_1
        real(dp), allocatable, intent(out) :: derivatives(:) !! Y_1' evaluated at each zero
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_y_zero_pair(1, nt, root_y, .true., zeros, derivatives, status)
    end subroutine y1_zeros

    pure subroutine y1p_zeros(nt, zeros, values, status)
        integer, intent(in) :: nt !! number of positive real zeros requested, >= 1
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of Y_1'
        real(dp), allocatable, intent(out) :: values(:) !! Y_1 evaluated at each derivative zero
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_y_zero_pair(1, nt, root_yp, .false., zeros, values, status)
    end subroutine y1p_zeros

    pure subroutine yn_zeros(n, nt, zeros, status)
        integer, intent(in) :: n !! integer Bessel order; negative orders use abs(n)
        integer, intent(in) :: nt !! number of positive zeros requested, >= 1
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of Y_n
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_zeros(n, nt, root_y, zeros, status)
    end subroutine yn_zeros

    pure subroutine ynp_zeros(n, nt, zeros, status)
        integer, intent(in) :: n !! integer Bessel order; negative orders use abs(n)
        integer, intent(in) :: nt !! number of positive zeros requested, >= 1
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of Y_n'
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_zeros(n, nt, root_yp, zeros, status)
    end subroutine ynp_zeros

    pure subroutine jnyn_zeros(n, nt, jzeros, jpzeros, yzeros, ypzeros, status)
        integer, intent(in) :: n !! integer Bessel order; negative orders use abs(n)
        integer, intent(in) :: nt !! number of positive zeros requested for each family, >= 1
        real(dp), allocatable, intent(out) :: jzeros(:) !! first nt positive zeros of J_n
        real(dp), allocatable, intent(out) :: jpzeros(:) !! first nt positive zeros of J_n'
        real(dp), allocatable, intent(out) :: yzeros(:) !! first nt positive zeros of Y_n
        real(dp), allocatable, intent(out) :: ypzeros(:) !! first nt positive zeros of Y_n'
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: info

        call compute_zeros(n, nt, root_j, jzeros, info)
        if (info /= bessel_zeros_status_success) then
            allocate(jpzeros(0), yzeros(0), ypzeros(0))
            if (present(status)) status = info
            return
        end if
        call compute_zeros(n, nt, root_jp, jpzeros, info)
        if (info /= bessel_zeros_status_success) then
            allocate(yzeros(0), ypzeros(0))
            if (present(status)) status = info
            return
        end if
        call compute_zeros(n, nt, root_y, yzeros, info)
        if (info /= bessel_zeros_status_success) then
            allocate(ypzeros(0))
            if (present(status)) status = info
            return
        end if
        call compute_zeros(n, nt, root_yp, ypzeros, info)
        if (present(status)) status = info
    end subroutine jnyn_zeros

    pure subroutine compute_y_zero_pair(n, nt, which, derivative_values, zeros, values, status)
        integer, intent(in) :: n !! nonnegative integer Bessel Y order
        integer, intent(in) :: nt !! number of positive real zeros requested, >= 1
        integer, intent(in) :: which !! internal selector for Y or Y-prime roots
        logical, intent(in) :: derivative_values !! true to return Y_n'; false to return Y_n
        real(dp), allocatable, intent(out) :: zeros(:) !! computed positive zeros
        real(dp), allocatable, intent(out) :: values(:) !! companion values evaluated at the zeros
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: i
        integer :: info

        call compute_zeros(n, nt, which, zeros, info)
        if (info /= bessel_zeros_status_success) then
            allocate(values(0))
            if (present(status)) status = info
            return
        end if

        allocate(values(nt))
        do i = 1, nt
            if (derivative_values) then
                values(i) = yvp(real(n, dp), zeros(i))
            else
                values(i) = yv(real(n, dp), zeros(i))
            end if
            if (.not. ieee_is_finite(values(i))) then
                if (present(status)) status = bessel_zeros_status_numerical_failure
                return
            end if
        end do
        if (present(status)) status = bessel_zeros_status_success
    end subroutine compute_y_zero_pair

    pure subroutine compute_zeros(n, nt, which, zeros, status)
        integer, intent(in) :: n !! integer Bessel order
        integer, intent(in) :: nt !! number of positive zeros requested
        integer, intent(in) :: which !! internal root-family selector
        real(dp), allocatable, intent(out) :: zeros(:) !! computed positive zeros
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: count
        integer :: info
        integer :: max_scans
        integer :: scans
        real(dp) :: a
        real(dp) :: b
        real(dp) :: fa
        real(dp) :: fb
        real(dp) :: nu
        real(dp) :: root
        real(dp) :: step

        if (nt < 1) then
            allocate(zeros(0))
            if (present(status)) status = bessel_zeros_status_invalid_argument
            return
        end if

        allocate(zeros(nt))
        zeros = quiet_nan(0.0_dp)
        nu = abs(real(n, dp))
        step = 0.25_dp * scifort_pi
        if (nu == 0.0_dp) then
            a = sqrt(epsilon(1.0_dp))
        else
            a = nu
        end if
        fa = bessel_root_value(which, nu, a)
        if (ieee_is_nan(fa) .or. .not. ieee_is_finite(fa)) then
            if (present(status)) status = bessel_zeros_status_numerical_failure
            return
        end if

        max_scans = 64 + 8 * nt + ceiling(12.0_dp * max(1.0_dp, nu) ** (1.0_dp / 3.0_dp))
        count = 0
        scans = 0
        do while (count < nt .and. scans < max_scans)
            b = a + step
            fb = bessel_root_value(which, nu, b)
            scans = scans + 1
            if (ieee_is_nan(fb) .or. .not. ieee_is_finite(fb)) then
                if (present(status)) status = bessel_zeros_status_numerical_failure
                return
            end if

            if (fb == 0.0_dp) then
                root = b
                info = bessel_zeros_status_success
            else if (opposite_signs(fa, fb)) then
                call refine_root(which, nu, a, b, fa, fb, root, info)
            else
                info = bessel_zeros_status_invalid_argument
            end if

            if (info == bessel_zeros_status_success) then
                count = count + 1
                zeros(count) = root
            else if (info == bessel_zeros_status_numerical_failure) then
                if (present(status)) status = info
                return
            end if

            a = b
            fa = fb
        end do

        if (count == nt) then
            if (present(status)) status = bessel_zeros_status_success
        else
            if (present(status)) status = bessel_zeros_status_numerical_failure
        end if
    end subroutine compute_zeros

    pure subroutine select_first_jnjnp(nt, candidates, candidate_orders, &
            candidate_serials, candidate_kinds, zeros, orders, serials, kinds)
        integer, intent(in) :: nt !! number of output entries to select
        real(dp), intent(in) :: candidates(:) !! candidate root magnitudes
        integer, intent(in) :: candidate_orders(:) !! Bessel order for each candidate
        integer, intent(in) :: candidate_serials(:) !! serial root number for each candidate
        integer, intent(in) :: candidate_kinds(:) !! zero for J and one for J-prime
        real(dp), intent(out) :: zeros(:) !! selected roots in increasing order
        integer, intent(out) :: orders(:) !! selected Bessel orders
        integer, intent(out) :: serials(:) !! selected serial root numbers
        integer, intent(out) :: kinds(:) !! selected zero-family indicators

        integer :: best
        integer :: i
        integer :: k
        logical, allocatable :: used(:)

        allocate(used(size(candidates)))
        used = .false.
        do k = 1, nt
            best = 0
            do i = 1, size(candidates)
                if (used(i)) cycle
                if (best == 0) then
                    best = i
                else if (jnjnp_candidate_precedes(candidates(i), candidate_orders(i), &
                        candidate_serials(i), candidate_kinds(i), candidates(best), &
                        candidate_orders(best), candidate_serials(best), &
                        candidate_kinds(best))) then
                    best = i
                end if
            end do
            used(best) = .true.
            zeros(k) = candidates(best)
            orders(k) = candidate_orders(best)
            serials(k) = candidate_serials(best)
            kinds(k) = candidate_kinds(best)
        end do
    end subroutine select_first_jnjnp

    pure logical function jnjnp_candidate_precedes(xa, na, ma, ta, xb, nb, mb, tb) result(precedes)
        real(dp), intent(in) :: xa !! first candidate root magnitude
        integer, intent(in) :: na !! first candidate Bessel order
        integer, intent(in) :: ma !! first candidate serial root number
        integer, intent(in) :: ta !! first candidate zero-family indicator
        real(dp), intent(in) :: xb !! second candidate root magnitude
        integer, intent(in) :: nb !! second candidate Bessel order
        integer, intent(in) :: mb !! second candidate serial root number
        integer, intent(in) :: tb !! second candidate zero-family indicator

        if (xa /= xb) then
            precedes = xa < xb
        else if (ta /= tb) then
            precedes = ta < tb
        else if (na /= nb) then
            precedes = na < nb
        else
            precedes = ma < mb
        end if
    end function jnjnp_candidate_precedes

    pure subroutine refine_root(which, nu, left, right, fleft, fright, root, status)
        integer, intent(in) :: which !! internal root-family selector
        real(dp), intent(in) :: nu !! nonnegative integer-valued Bessel order
        real(dp), intent(in) :: left !! lower sign-changing bracket endpoint
        real(dp), intent(in) :: right !! upper sign-changing bracket endpoint
        real(dp), intent(in) :: fleft !! root function at left
        real(dp), intent(in) :: fright !! root function at right
        real(dp), intent(out) :: root !! refined positive root
        integer, intent(out) :: status !! zero on success; numerical-failure code otherwise

        integer :: iter
        real(dp) :: a
        real(dp) :: b
        real(dp) :: candidate
        real(dp) :: fa
        real(dp) :: fb
        real(dp) :: fx
        real(dp) :: slope
        real(dp) :: tol
        real(dp) :: x

        a = left
        b = right
        fa = fleft
        fb = fright
        x = 0.5_dp * (a + b)

        do iter = 1, 64
            fx = bessel_root_value(which, nu, x)
            if (ieee_is_nan(fx) .or. .not. ieee_is_finite(fx)) then
                root = quiet_nan(x)
                status = bessel_zeros_status_numerical_failure
                return
            end if
            if (fx == 0.0_dp) then
                root = x
                status = bessel_zeros_status_success
                return
            end if

            if (opposite_signs(fa, fx)) then
                b = x
                fb = fx
            else
                a = x
                fa = fx
            end if

            tol = 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x))
            if (b - a <= tol) then
                root = 0.5_dp * (a + b)
                status = bessel_zeros_status_success
                return
            end if

            slope = bessel_root_slope(which, nu, x)
            candidate = 0.5_dp * (a + b)
            if (ieee_is_finite(slope)) then
                if (slope /= 0.0_dp) then
                    candidate = x - fx / slope
                    if (.not. ieee_is_finite(candidate)) candidate = 0.5_dp * (a + b)
                    if (candidate <= a .or. candidate >= b) candidate = 0.5_dp * (a + b)
                end if
            end if
            x = candidate
        end do

        root = 0.5_dp * (a + b)
        if (b - a <= 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(root))) then
            status = bessel_zeros_status_success
        else
            status = bessel_zeros_status_numerical_failure
        end if
    end subroutine refine_root

    pure function bessel_root_value(which, nu, x) result(value)
        integer, intent(in) :: which !! internal root-family selector
        real(dp), intent(in) :: nu !! nonnegative integer-valued Bessel order
        real(dp), intent(in) :: x !! positive real argument
        real(dp) :: value

        select case (which)
        case (root_j)
            value = jv(nu, x)
        case (root_jp)
            value = jvp(nu, x)
        case (root_y)
            value = yv(nu, x)
        case (root_yp)
            value = yvp(nu, x)
        case default
            value = quiet_nan(x)
        end select
    end function bessel_root_value

    pure function bessel_root_slope(which, nu, x) result(value)
        integer, intent(in) :: which !! internal root-family selector
        real(dp), intent(in) :: nu !! nonnegative integer-valued Bessel order
        real(dp), intent(in) :: x !! positive real argument
        real(dp) :: value

        select case (which)
        case (root_j)
            value = jvp(nu, x)
        case (root_jp)
            value = jvp(nu, x, 2)
        case (root_y)
            value = yvp(nu, x)
        case (root_yp)
            value = yvp(nu, x, 2)
        case default
            value = quiet_nan(x)
        end select
    end function bessel_root_slope

    pure logical function opposite_signs(a, b) result(different)
        real(dp), intent(in) :: a !! first finite nonzero value
        real(dp), intent(in) :: b !! second finite nonzero value

        different = (a < 0.0_dp .and. b > 0.0_dp) .or. &
            (a > 0.0_dp .and. b < 0.0_dp)
    end function opposite_signs

end module scifort_bessel_zeros
