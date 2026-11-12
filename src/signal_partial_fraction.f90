! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Partial-fraction expansion of analog and digital transfer functions.
! Analog coefficient arrays are in descending powers of s; digital arrays
! are in ascending powers of z**(-1), as in scipy.signal.
module scifort_signal_partial_fraction
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_analog, only : polynomial_roots
    implicit none
    private

    integer, parameter, public :: partial_fraction_success = 0
    integer, parameter, public :: partial_fraction_invalid_input = 1
    integer, parameter, public :: partial_fraction_singular = 2

    type, public :: partial_fraction_result
        complex(dp), allocatable :: residues(:) !! coefficients of ascending inverse powers
        complex(dp), allocatable :: poles(:) !! poles, repeated once for each inverse power
        complex(dp), allocatable :: direct(:) !! direct polynomial coefficients
        integer :: status = partial_fraction_invalid_input
    end type partial_fraction_result

    type, public :: partial_fraction_tf_result
        complex(dp), allocatable :: b(:) !! transfer-function numerator
        complex(dp), allocatable :: a(:) !! monic transfer-function denominator
        integer :: status = partial_fraction_invalid_input
    end type partial_fraction_tf_result

    interface residue
        module procedure residue_real
        module procedure residue_complex
    end interface residue
    interface residuez
        module procedure residuez_real
        module procedure residuez_complex
    end interface residuez
    interface invres
        module procedure invres_complex
        module procedure invres_real_direct
    end interface invres
    interface invresz
        module procedure invresz_complex
        module procedure invresz_real_direct
    end interface invresz
    public :: residue, residuez, invres, invresz

contains

    pure logical function finite_complex(v) result(ok)
        complex(dp), intent(in) :: v(:) !! values to check
        ok = all(ieee_is_finite(real(v,dp))) .and. all(ieee_is_finite(aimag(v)))
    end function finite_complex

    function trim_poly(v, digital) result(out)
        complex(dp), intent(in) :: v(:) !! input coefficient vector
        logical, intent(in) :: digital !! trim trailing zeros if true
        complex(dp), allocatable :: out(:)
        integer :: first, last

        first = 1
        last = size(v)
        if (digital) then
            do while (last >= first)
                if (abs(v(last)) /= 0.0_dp) exit
                last = last - 1
            end do
        else
            do while (first <= last)
                if (abs(v(first)) /= 0.0_dp) exit
                first = first + 1
            end do
        end if
        out = v(first:last)
    end function trim_poly

    pure function multiply_poly(a, b) result(c)
        complex(dp), intent(in) :: a(:) !! first polynomial
        complex(dp), intent(in) :: b(:) !! second polynomial
        complex(dp), allocatable :: c(:)
        integer :: i, j

        if (size(a) == 0 .or. size(b) == 0) then
            allocate(c(0))
            return
        end if
        allocate(c(size(a)+size(b)-1))
        c = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1, size(b)
            do i = 1, size(a)
                c(i+j-1) = c(i+j-1) + a(i)*b(j)
            end do
        end do
    end function multiply_poly

    subroutine divide_desc(b, a, q, rem)
        complex(dp), intent(in) :: b(:) !! dividend, descending powers
        complex(dp), intent(in) :: a(:) !! divisor, descending powers and nonzero a(1)
        complex(dp), allocatable, intent(out) :: q(:) !! polynomial quotient
        complex(dp), allocatable, intent(out) :: rem(:) !! remainder
        complex(dp), allocatable :: temp(:)
        integer :: i, nq

        if (size(b) < size(a)) then
            allocate(q(0))
            rem = b
            return
        end if
        nq = size(b)-size(a)+1
        allocate(q(nq))
        temp = b
        do i = 1, nq
            q(i) = temp(i)/a(1)
            temp(i:i+size(a)-1) = temp(i:i+size(a)-1) - q(i)*a
        end do
        rem = temp(nq+1:)
    end subroutine divide_desc

    pure logical function less_pole(x, y) result(less)
        complex(dp), intent(in) :: x !! first pole
        complex(dp), intent(in) :: y !! second pole
        real(dp) :: ax, ay

        ax = abs(x)
        ay = abs(y)
        if (ax /= ay) then
            less = ax < ay
        else if (real(x,dp) /= real(y,dp)) then
            less = real(x,dp) < real(y,dp)
        else
            less = aimag(x) < aimag(y)
        end if
    end function less_pole

    subroutine group_poles(p, tolerance, rtype, sort_first, unique, counts, ok)
        complex(dp), intent(in) :: p(:) !! poles, with near-duplicates adjacent when not sorting
        real(dp), intent(in) :: tolerance !! absolute grouping tolerance
        character(*), intent(in) :: rtype !! representative choice
        logical, intent(in) :: sort_first !! sort by magnitude for forward transforms
        complex(dp), allocatable, intent(out) :: unique(:) !! representative of each group
        integer, allocatable, intent(out) :: counts(:) !! group multiplicities
        logical, intent(out) :: ok !! true when grouping arguments are valid
        complex(dp), allocatable :: work(:)
        complex(dp) :: item, representative
        integer :: i, j, n, first, ng

        ok = .false.
        allocate(unique(0),counts(0))
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        select case (trim(rtype))
        case ('avg','mean','min','minimum','max','maximum')
        case default
            return
        end select
        if (.not. finite_complex(p)) return
        work = p
        n = size(work)
        if (sort_first) then
            do i = 2, n
                item = work(i)
                j = i-1
                do while (j >= 1)
                    if (.not. less_pole(item,work(j))) exit
                    work(j+1) = work(j)
                    j = j-1
                end do
                work(j+1) = item
            end do
        end if
        first = 1
        do while (first <= n)
            j = first+1
            do while (j <= n)
                if (abs(work(j)-work(first)) > tolerance) exit
                j = j+1
            end do
            representative = work(first)
            select case (trim(rtype))
            case ('avg','mean')
                representative = sum(work(first:j-1))/real(j-first,dp)
            case ('min','minimum')
                do i = first+1, j-1
                    if (real(work(i),dp) < real(representative,dp)) then
                        representative = work(i)
                    else if (real(work(i),dp) == real(representative,dp) .and. &
                        aimag(work(i)) < aimag(representative)) then
                        representative = work(i)
                    end if
                end do
            case ('max','maximum')
                do i = first+1, j-1
                    if (real(work(i),dp) > real(representative,dp)) then
                        representative = work(i)
                    else if (real(work(i),dp) == real(representative,dp) .and. &
                        aimag(work(i)) > aimag(representative)) then
                        representative = work(i)
                    end if
                end do
            end select
            ng = size(unique)
            unique = [unique,representative]
            counts = [counts,j-first]
            first = j
        end do
        ok = .true.
    end subroutine group_poles

    function factor_product(unique, counts, omit_group, omit_powers) result(out)
        complex(dp), intent(in) :: unique(:) !! distinct pole representatives
        integer, intent(in) :: counts(:) !! multiplicities
        integer, intent(in), optional :: omit_group !! group to reduce by omit_powers
        integer, intent(in), optional :: omit_powers !! number of factors to omit
        complex(dp), allocatable :: out(:)
        complex(dp) :: linear(2)
        integer :: i, j, n

        out = [cmplx(1.0_dp,0.0_dp,dp)]
        do i = 1, size(unique)
            linear = [cmplx(1.0_dp,0.0_dp,dp),-unique(i)]
            n = counts(i)
            if (present(omit_group) .and. present(omit_powers)) then
                if (i == omit_group) n = n-omit_powers
            end if
            do j = 1, n
                out = multiply_poly(out,linear)
            end do
        end do
    end function factor_product

    subroutine solve_complex(a, b, x, ok)
        complex(dp), intent(in) :: a(:,:) !! nonsingular coefficient matrix
        complex(dp), intent(in) :: b(:) !! right-hand side
        complex(dp), allocatable, intent(out) :: x(:) !! solution
        logical, intent(out) :: ok !! true for finite nonsingular solution
        complex(dp), allocatable :: work(:,:), rhs(:), temp_row(:)
        complex(dp) :: ratio, temp
        integer :: n, k, i, pivot

        ok = .false.
        n = size(b)
        allocate(x(n))
        if (n == 0) then
            ok = .true.
            return
        end if
        work = a
        rhs = b
        do k = 1, n
            pivot = k-1+maxloc(abs(work(k:n,k)),dim=1)
            if (abs(work(pivot,k)) <= tiny(1.0_dp)) return
            if (pivot /= k) then
                temp_row = work(k,:)
                work(k,:) = work(pivot,:)
                work(pivot,:) = temp_row
                temp = rhs(k)
                rhs(k) = rhs(pivot)
                rhs(pivot) = temp
            end if
            do i = k+1, n
                ratio = work(i,k)/work(k,k)
                work(i,k:n) = work(i,k:n)-ratio*work(k,k:n)
                rhs(i) = rhs(i)-ratio*rhs(k)
            end do
        end do
        do i = n, 1, -1
            x(i) = (rhs(i)-sum(work(i,i+1:n)*x(i+1:n)))/work(i,i)
        end do
        ok = finite_complex(x)
    end subroutine solve_complex

    function residue_real(b,a,tol,rtype) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in), optional :: tol !! pole grouping tolerance, default 1e-3
        character(*), intent(in), optional :: rtype !! avg, min or max representative
        type(partial_fraction_result) :: out
        out = residue_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),.false.,tol,rtype)
    end function residue_real

    function residue_complex(b,a,tol,rtype) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in), optional :: tol !! pole grouping tolerance, default 1e-3
        character(*), intent(in), optional :: rtype !! avg, min or max representative
        type(partial_fraction_result) :: out
        out = residue_core(b,a,.false.,tol,rtype)
    end function residue_complex

    function residuez_real(b,a,tol,rtype) result(out)
        real(dp), intent(in) :: b(:) !! numerator, ascending powers of z**(-1)
        real(dp), intent(in) :: a(:) !! denominator, ascending powers of z**(-1)
        real(dp), intent(in), optional :: tol !! pole grouping tolerance, default 1e-3
        character(*), intent(in), optional :: rtype !! avg, min or max representative
        type(partial_fraction_result) :: out
        out = residue_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),.true.,tol,rtype)
    end function residuez_real

    function residuez_complex(b,a,tol,rtype) result(out)
        complex(dp), intent(in) :: b(:) !! numerator, ascending powers of z**(-1)
        complex(dp), intent(in) :: a(:) !! denominator, ascending powers of z**(-1)
        real(dp), intent(in), optional :: tol !! pole grouping tolerance, default 1e-3
        character(*), intent(in), optional :: rtype !! avg, min or max representative
        type(partial_fraction_result) :: out
        out = residue_core(b,a,.true.,tol,rtype)
    end function residuez_complex

    function residue_core(b,a,digital,tol,rtype) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        logical, intent(in) :: digital !! true for z-inverse coefficients
        real(dp), intent(in), optional :: tol !! absolute grouping tolerance
        character(*), intent(in), optional :: rtype !! grouping representative
        type(partial_fraction_result) :: out
        complex(dp), allocatable :: bc(:), ac(:), poles(:), unique(:), qq(:), rem(:), basis(:)
        complex(dp), allocatable :: mat(:,:), rhs(:), residues(:)
        integer, allocatable :: counts(:)
        integer :: n, i, col, offset, power
        real(dp) :: tolerance
        character(:), allocatable :: kind
        logical :: ok

        if (size(a) < 1 .or. .not. finite_complex(a) .or. .not. finite_complex(b)) return
        ac = trim_poly(a,digital)
        bc = trim_poly(b,digital)
        if (size(ac) == 0) return
        if (digital) then
            if (abs(ac(1)) == 0.0_dp) return
        end if
        tolerance = 1.0e-3_dp
        if (present(tol)) tolerance = tol
        kind = 'avg'
        if (present(rtype)) kind = trim(rtype)
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        call polynomial_roots(ac,poles,ok)
        if (.not. ok) then
            out%status = partial_fraction_singular
            return
        end if
        call group_poles(poles,tolerance,kind,.true.,unique,counts,ok)
        if (.not. ok) return
        ! Polynomial roots at real repeated poles can have small spurious
        ! imaginary parts. Project these representatives onto the real axis
        ! before computing large, cancellation-sensitive residues.
        if (all(aimag(ac) == 0.0_dp)) then
            do i = 1, size(unique)
                if (abs(aimag(unique(i))) <= &
                    min(tolerance,1.0e-7_dp*(1.0_dp+abs(real(unique(i),dp))))) then
                    unique(i) = cmplx(real(unique(i),dp),0.0_dp,dp)
                end if
            end do
        end if
        n = sum(counts)
        allocate(out%poles(n),out%residues(n))
        offset = 0
        do i = 1, size(unique)
            out%poles(offset+1:offset+counts(i)) = unique(i)
            offset = offset+counts(i)
        end do
        out%residues = cmplx(0.0_dp,0.0_dp,dp)
        if (size(bc) == 0) then
            allocate(out%direct(0))
            out%status = partial_fraction_success
            return
        end if
        if (digital) then
            call divide_desc(bc(size(bc):1:-1),ac(size(ac):1:-1),qq,rem)
            out%direct = qq(size(qq):1:-1)
            rem = rem(size(rem):1:-1)
        else
            call divide_desc(bc,ac,qq,rem)
            out%direct = qq
        end if
        if (n == 0) then
            out%status = partial_fraction_success
            return
        end if
        allocate(mat(n,n),rhs(n))
        rhs = cmplx(0.0_dp,0.0_dp,dp)
        if (digital) then
            rhs(1:size(rem)) = rem/ac(1)
        else
            rhs(n-size(rem)+1:n) = rem/ac(1)
        end if
        col = 0
        do i = 1, size(unique)
            do power = 1, counts(i)
                basis = factor_product(unique,counts,i,power)
                col = col+1
                mat(:,col) = cmplx(0.0_dp,0.0_dp,dp)
                if (digital) then
                    mat(1:size(basis),col) = basis
                else
                    mat(n-size(basis)+1:n,col) = basis
                end if
            end do
        end do
        call solve_complex(mat,rhs,residues,ok)
        if (.not. ok) then
            out%status = partial_fraction_singular
            return
        end if
        if (maxval(abs(matmul(mat,residues)-rhs)) > &
            1.0e-8_dp*(1.0_dp+maxval(abs(rhs)))) then
            out%status = partial_fraction_singular
            return
        end if
        out%residues = residues
        out%status = partial_fraction_success
    end function residue_core

    function invres_complex(r,p,k,tol,rtype) result(out)
        complex(dp), intent(in) :: r(:) !! residues, ascending inverse pole powers
        complex(dp), intent(in) :: p(:) !! corresponding poles
        complex(dp), intent(in), optional :: k(:) !! direct polynomial coefficients
        real(dp), intent(in), optional :: tol !! pole grouping tolerance
        character(*), intent(in), optional :: rtype !! avg, min or max representative
        type(partial_fraction_tf_result) :: out
        complex(dp), allocatable :: direct(:)
        allocate(direct(0))
        if (present(k)) direct = k
        out = inverse_core(r,p,direct,.false.,tol,rtype)
    end function invres_complex

    function invres_real_direct(r,p,k,tol,rtype) result(out)
        complex(dp), intent(in) :: r(:) !! residues
        complex(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: k(:) !! real direct polynomial coefficients
        real(dp), intent(in), optional :: tol !! pole grouping tolerance
        character(*), intent(in), optional :: rtype !! avg, min or max representative
        type(partial_fraction_tf_result) :: out
        out = inverse_core(r,p,cmplx(k,0.0_dp,dp),.false.,tol,rtype)
    end function invres_real_direct

    function invresz_complex(r,p,k,tol,rtype) result(out)
        complex(dp), intent(in) :: r(:) !! residues, ascending inverse pole powers
        complex(dp), intent(in) :: p(:) !! corresponding poles
        complex(dp), intent(in), optional :: k(:) !! direct delay polynomial coefficients
        real(dp), intent(in), optional :: tol !! pole grouping tolerance
        character(*), intent(in), optional :: rtype !! avg, min or max representative
        type(partial_fraction_tf_result) :: out
        complex(dp), allocatable :: direct(:)
        allocate(direct(0))
        if (present(k)) direct = k
        out = inverse_core(r,p,direct,.true.,tol,rtype)
    end function invresz_complex

    function invresz_real_direct(r,p,k,tol,rtype) result(out)
        complex(dp), intent(in) :: r(:) !! residues
        complex(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: k(:) !! real direct delay polynomial coefficients
        real(dp), intent(in), optional :: tol !! pole grouping tolerance
        character(*), intent(in), optional :: rtype !! avg, min or max representative
        type(partial_fraction_tf_result) :: out
        out = inverse_core(r,p,cmplx(k,0.0_dp,dp),.true.,tol,rtype)
    end function invresz_real_direct

    function inverse_core(r,p,k,digital,tol,rtype) result(out)
        complex(dp), intent(in) :: r(:) !! residues
        complex(dp), intent(in) :: p(:) !! grouped poles
        complex(dp), intent(in) :: k(:) !! direct polynomial coefficients
        logical, intent(in) :: digital !! true for negative powers of z
        real(dp), intent(in), optional :: tol !! pole grouping tolerance
        character(*), intent(in), optional :: rtype !! grouping representative
        type(partial_fraction_tf_result) :: out
        complex(dp), allocatable :: unique(:), a(:), b(:), direct(:), basis(:), contribution(:)
        integer, allocatable :: counts(:)
        integer :: n, i, power, index, numerator_size
        real(dp) :: tolerance
        character(:), allocatable :: kind
        logical :: ok

        if (size(r) /= size(p)) return
        if (.not. finite_complex(r) .or. .not. finite_complex(p) .or. .not. finite_complex(k)) return
        tolerance = 1.0e-3_dp
        if (present(tol)) tolerance = tol
        kind = 'avg'
        if (present(rtype)) kind = trim(rtype)
        call group_poles(p,tolerance,kind,.false.,unique,counts,ok)
        if (.not. ok) return
        n = sum(counts)
        a = factor_product(unique,counts)
        direct = trim_poly(k,digital)
        numerator_size = n
        if (size(direct) /= 0) numerator_size = max(numerator_size,size(direct)+n)
        numerator_size = max(1,numerator_size)
        allocate(b(numerator_size))
        b = cmplx(0.0_dp,0.0_dp,dp)
        if (size(direct) /= 0) then
            contribution = multiply_poly(direct,a)
            if (digital) then
                b(1:size(contribution)) = contribution
            else
                b(numerator_size-size(contribution)+1:) = contribution
            end if
        end if
        index = 0
        do i = 1, size(unique)
            do power = 1, counts(i)
                index = index+1
                basis = factor_product(unique,counts,i,power)
                if (digital) then
                    b(1:size(basis)) = b(1:size(basis)) + r(index)*basis
                else
                    b(numerator_size-size(basis)+1:) = &
                        b(numerator_size-size(basis)+1:) + r(index)*basis
                end if
            end do
        end do
        if (.not. finite_complex(a) .or. .not. finite_complex(b)) then
            out%status = partial_fraction_singular
            return
        end if
        out%a = a
        out%b = trim_poly(b,digital)
        if (size(out%b) == 0) out%b = [cmplx(0.0_dp,0.0_dp,dp)]
        out%status = partial_fraction_success
    end function inverse_core
end module scifort_signal_partial_fraction
