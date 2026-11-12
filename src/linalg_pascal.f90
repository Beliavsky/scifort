! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Pascal matrices and their inverses, corresponding to scipy.linalg.pascal
! and scipy.linalg.invpascal. The exact decimal interfaces are SciFort extensions.
!
! The formulas follow the binomial identity for Pascal's triangle and the
! factorization P = L L^T, where L(i,j) = C(i-1,j-1) for i >= j.
module scifort_linalg_pascal
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    implicit none
    private

    integer(int64), parameter :: limb_base = 1000000000_int64
    integer, parameter :: symmetric_kind = 1, lower_kind = 2, upper_kind = 3

    type :: positive_integer
        integer(int64), allocatable :: digit(:) ! least-significant limb first; absent means zero
    end type positive_integer

    public :: pascal, invpascal, pascal_exact_decimal, invpascal_exact_decimal

    interface pascal
        module procedure pascal_real
        module procedure pascal_int64
    end interface pascal

    interface invpascal
        module procedure invpascal_real
        module procedure invpascal_int64
    end interface invpascal

contains

    pure integer function parse_kind(kind) result(which)
        character(*), intent(in), optional :: kind !! 'symmetric' (default), 'lower', or 'upper'

        which = symmetric_kind
        if (.not. present(kind)) return
        select case (kind)
        case ('symmetric')
            which = symmetric_kind
        case ('lower')
            which = lower_kind
        case ('upper')
            which = upper_kind
        case default
            which = 0
        end select
    end function parse_kind

    pure subroutine pascal_real(n, p, status, kind)
        integer, intent(in) :: n !! order, nonnegative
        real(dp), intent(out) :: p(:,:) !! output Pascal matrix, shape (n,n)
        integer, intent(out) :: status !! 0 success; 1 invalid shape/kind; 4 numeric overflow
        character(*), intent(in), optional :: kind !! 'symmetric' (default), 'lower', or 'upper'
        real(dp), allocatable :: work(:,:)
        integer :: which, i, j

        p = 0.0_dp
        status = linalg_status_invalid_shape
        which = parse_kind(kind)
        if (n < 0 .or. which == 0) return
        if (size(p,1) /= n .or. size(p,2) /= n) return
        status = linalg_status_success
        if (n == 0) return
        if (which == symmetric_kind) then
            p(1,:) = 1.0_dp
            p(:,1) = 1.0_dp
            do j=2,n
                do i=2,n
                    if (p(i-1,j) > huge(1.0_dp)-p(i,j-1)) then
                        p = 0.0_dp
                        status = linalg_status_no_convergence
                        return
                    end if
                    p(i,j) = p(i-1,j)+p(i,j-1)
                end do
            end do
        else
            allocate(work(n,n))
            call build_lower_real(n,work,status)
            if (status /= linalg_status_success) return
            if (which == lower_kind) then
                p = work
            else
                p = transpose(work)
            end if
        end if
    end subroutine pascal_real

    pure subroutine pascal_int64(n, p, status, kind)
        integer, intent(in) :: n !! order, nonnegative
        integer(int64), intent(out) :: p(:,:) !! exact Pascal entries when representable, shape (n,n)
        integer, intent(out) :: status !! 0 success; 1 invalid shape/kind; 4 integer overflow
        character(*), intent(in), optional :: kind !! 'symmetric' (default), 'lower', or 'upper'
        integer(int64), allocatable :: work(:,:)
        integer :: which, i, j

        p = 0_int64
        status = linalg_status_invalid_shape
        which = parse_kind(kind)
        if (n < 0 .or. which == 0) return
        if (size(p,1) /= n .or. size(p,2) /= n) return
        status = linalg_status_success
        if (n == 0) return
        if (which == symmetric_kind) then
            p(1,:) = 1_int64
            p(:,1) = 1_int64
            do j=2,n
                do i=2,n
                    if (p(i-1,j) > huge(1_int64)-p(i,j-1)) then
                        p = 0_int64
                        status = linalg_status_no_convergence
                        return
                    end if
                    p(i,j) = p(i-1,j)+p(i,j-1)
                end do
            end do
        else
            allocate(work(n,n))
            call build_lower_int64(n,work,status)
            if (status /= linalg_status_success) return
            if (which == lower_kind) then
                p = work
            else
                p = transpose(work)
            end if
        end if
    end subroutine pascal_int64

    pure subroutine invpascal_real(n, p, status, kind)
        integer, intent(in) :: n !! order, nonnegative
        real(dp), intent(out) :: p(:,:) !! inverse Pascal matrix, shape (n,n)
        integer, intent(out) :: status !! 0 success; 1 invalid shape/kind; 4 numeric overflow
        character(*), intent(in), optional :: kind !! 'symmetric' (default), 'lower', or 'upper'
        real(dp), allocatable :: lower(:,:), work(:,:)
        real(dp) :: term, accumulation, limit
        integer :: which, i, j, k

        p = 0.0_dp
        status = linalg_status_invalid_shape
        which = parse_kind(kind)
        if (n < 0 .or. which == 0) return
        if (size(p,1) /= n .or. size(p,2) /= n) return
        status = linalg_status_success
        if (n == 0) return
        allocate(lower(n,n),work(n,n))
        call build_lower_real(n,lower,status)
        if (status /= linalg_status_success) return
        work = 0.0_dp
        if (which == symmetric_kind) then
            limit = huge(1.0_dp)
            do j=1,n
                do i=j,n
                    accumulation = 0.0_dp
                    do k=i,n
                        if (lower(k,i) > 0.0_dp .and. lower(k,j) > limit/lower(k,i)) exit
                        term = lower(k,i)*lower(k,j)
                        if (accumulation > limit-term) exit
                        accumulation = accumulation+term
                    end do
                    if (k <= n) then
                        status = linalg_status_no_convergence
                        return
                    end if
                    if (mod(i+j,2) /= 0) accumulation = -accumulation
                    work(i,j) = accumulation
                    work(j,i) = accumulation
                end do
            end do
        else
            do i=1,n
                do j=1,i
                    work(i,j) = lower(i,j)
                    if (mod(i-j,2) /= 0) work(i,j) = -work(i,j)
                end do
            end do
            if (which == upper_kind) work = transpose(work)
        end if
        p = work
    end subroutine invpascal_real

    pure subroutine invpascal_int64(n, p, status, kind)
        integer, intent(in) :: n !! order, nonnegative
        integer(int64), intent(out) :: p(:,:) !! exact inverse Pascal matrix, shape (n,n)
        integer, intent(out) :: status !! 0 success; 1 invalid shape/kind; 4 integer overflow
        character(*), intent(in), optional :: kind !! 'symmetric' (default), 'lower', or 'upper'
        integer(int64), allocatable :: lower(:,:), work(:,:)
        integer(int64) :: term, accumulation, limit
        integer :: which, i, j, k

        p = 0_int64
        status = linalg_status_invalid_shape
        which = parse_kind(kind)
        if (n < 0 .or. which == 0) return
        if (size(p,1) /= n .or. size(p,2) /= n) return
        status = linalg_status_success
        if (n == 0) return
        allocate(lower(n,n),work(n,n))
        call build_lower_int64(n,lower,status)
        if (status /= linalg_status_success) return
        work = 0_int64
        if (which == symmetric_kind) then
            limit = huge(1_int64)
            do j=1,n
                do i=j,n
                    accumulation = 0_int64
                    do k=i,n
                        if (lower(k,i) > 0_int64 .and. lower(k,j) > limit/lower(k,i)) exit
                        term = lower(k,i)*lower(k,j)
                        if (accumulation > limit-term) exit
                        accumulation = accumulation+term
                    end do
                    if (k <= n) then
                        status = linalg_status_no_convergence
                        return
                    end if
                    if (mod(i+j,2) /= 0) accumulation = -accumulation
                    work(i,j) = accumulation
                    work(j,i) = accumulation
                end do
            end do
        else
            do i=1,n
                do j=1,i
                    work(i,j) = lower(i,j)
                    if (mod(i-j,2) /= 0) work(i,j) = -work(i,j)
                end do
            end do
            if (which == upper_kind) work = transpose(work)
        end if
        p = work
    end subroutine invpascal_int64

    pure subroutine build_lower_real(n, lower, status)
        integer, intent(in) :: n !! nonnegative matrix order
        real(dp), intent(out) :: lower(n,n) !! positive lower-triangular binomial matrix
        integer, intent(out) :: status !! zero, or numeric overflow
        integer :: i, j

        lower = 0.0_dp
        status = linalg_status_success
        if (n == 0) return
        lower(1,1) = 1.0_dp
        do i=2,n
            lower(i,1) = 1.0_dp
            lower(i,i) = 1.0_dp
            do j=2,i-1
                if (lower(i-1,j-1) > huge(1.0_dp)-lower(i-1,j)) then
                    status = linalg_status_no_convergence
                    return
                end if
                lower(i,j) = lower(i-1,j-1)+lower(i-1,j)
            end do
        end do
    end subroutine build_lower_real

    pure subroutine build_lower_int64(n, lower, status)
        integer, intent(in) :: n !! nonnegative matrix order
        integer(int64), intent(out) :: lower(n,n) !! positive exact lower-triangular binomial matrix
        integer, intent(out) :: status !! zero, or integer overflow
        integer :: i, j

        lower = 0_int64
        status = linalg_status_success
        if (n == 0) return
        lower(1,1) = 1_int64
        do i=2,n
            lower(i,1) = 1_int64
            lower(i,i) = 1_int64
            do j=2,i-1
                if (lower(i-1,j-1) > huge(1_int64)-lower(i-1,j)) then
                    status = linalg_status_no_convergence
                    return
                end if
                lower(i,j) = lower(i-1,j-1)+lower(i-1,j)
            end do
        end do
    end subroutine build_lower_int64

    subroutine pascal_exact_decimal(n, p, status, kind)
        integer, intent(in) :: n !! order, nonnegative
        character(*), intent(out) :: p(:,:) !! exact decimal matrix entries, shape (n,n)
        integer, intent(out) :: status !! 0 success; 1 bad shape/kind or insufficient string length
        character(*), intent(in), optional :: kind !! 'symmetric' (default), 'lower', or 'upper'
        type(positive_integer), allocatable :: matrix(:,:)
        integer :: which, i, j
        character(:), allocatable :: entry

        p = ''
        status = linalg_status_invalid_shape
        which = parse_kind(kind)
        if (n < 0 .or. which == 0) return
        if (size(p,1) /= n .or. size(p,2) /= n) return
        status = linalg_status_success
        if (n == 0) return
        allocate(matrix(n,n))
        if (which == symmetric_kind) then
            do i=1,n
                matrix(i,1) = big_one()
                matrix(1,i) = big_one()
            end do
            do j=2,n
                do i=2,n
                    matrix(i,j) = big_add(matrix(i-1,j),matrix(i,j-1))
                end do
            end do
        else
            call build_lower_big(n,matrix)
        end if
        do j=1,n
            do i=1,n
                if (which == upper_kind) then
                    entry = big_text(matrix(j,i))
                else
                    entry = big_text(matrix(i,j))
                end if
                if (len(entry) > len(p)) then
                    p = ''
                    status = linalg_status_invalid_shape
                    return
                end if
                p(i,j) = entry
            end do
        end do
    end subroutine pascal_exact_decimal

    subroutine invpascal_exact_decimal(n, p, status, kind)
        integer, intent(in) :: n !! order, nonnegative
        character(*), intent(out) :: p(:,:) !! exact signed decimal matrix entries, shape (n,n)
        integer, intent(out) :: status !! 0 success; 1 bad shape/kind or insufficient string length
        character(*), intent(in), optional :: kind !! 'symmetric' (default), 'lower', or 'upper'
        type(positive_integer), allocatable :: lower(:,:)
        type(positive_integer) :: accumulator, term
        integer :: which, i, j, k
        character(:), allocatable :: entry

        p = ''
        status = linalg_status_invalid_shape
        which = parse_kind(kind)
        if (n < 0 .or. which == 0) return
        if (size(p,1) /= n .or. size(p,2) /= n) return
        status = linalg_status_success
        if (n == 0) return
        allocate(lower(n,n))
        call build_lower_big(n,lower)
        if (which == symmetric_kind) then
            do j=1,n
                do i=j,n
                    accumulator = positive_integer()
                    do k=i,n
                        term = big_product(lower(k,i),lower(k,j))
                        accumulator = big_add(accumulator,term)
                    end do
                    entry = big_text(accumulator)
                    if (mod(i+j,2) /= 0) entry = '-'//entry
                    if (len(entry) > len(p)) then
                        p = ''
                        status = linalg_status_invalid_shape
                        return
                    end if
                    p(i,j) = entry
                    p(j,i) = entry
                end do
            end do
        else
            do j=1,n
                do i=1,n
                    if (which == upper_kind) then
                        entry = big_text(lower(j,i))
                    else
                        entry = big_text(lower(i,j))
                    end if
                    if (entry /= '0' .and. mod(i+j,2) /= 0) entry = '-'//entry
                    if (len(entry) > len(p)) then
                        p = ''
                        status = linalg_status_invalid_shape
                        return
                    end if
                    p(i,j) = entry
                end do
            end do
        end if
    end subroutine invpascal_exact_decimal

    pure subroutine build_lower_big(n, lower)
        integer, intent(in) :: n !! nonnegative matrix order
        type(positive_integer), intent(out) :: lower(n,n) !! lower binomial coefficients
        integer :: i, j

        if (n == 0) return
        lower(1,1) = big_one()
        do i=2,n
            lower(i,1) = big_one()
            lower(i,i) = big_one()
            do j=2,i-1
                lower(i,j) = big_add(lower(i-1,j-1),lower(i-1,j))
            end do
        end do
    end subroutine build_lower_big

    pure function big_one() result(value)
        type(positive_integer) :: value !! exact one
        allocate(value%digit(1))
        value%digit(1) = 1_int64
    end function big_one

    pure function big_add(a,b) result(c)
        type(positive_integer), intent(in) :: a !! first nonnegative integer
        type(positive_integer), intent(in) :: b !! second nonnegative integer
        type(positive_integer) :: c !! exact sum
        integer :: na, nb, n, i
        integer(int64) :: carry, t

        if (.not. allocated(a%digit)) then
            if (allocated(b%digit)) then
                allocate(c%digit(size(b%digit)))
                c%digit = b%digit
            end if
            return
        end if
        if (.not. allocated(b%digit)) then
            allocate(c%digit(size(a%digit)))
            c%digit = a%digit
            return
        end if
        na = size(a%digit)
        nb = size(b%digit)
        n = max(na,nb)
        allocate(c%digit(n+1))
        c%digit = 0_int64
        carry = 0_int64
        do i=1,n
            t = carry
            if (i <= na) t = t+a%digit(i)
            if (i <= nb) t = t+b%digit(i)
            c%digit(i) = mod(t,limb_base)
            carry = t/limb_base
        end do
        if (carry == 0_int64) then
            c%digit = c%digit(:n)
        else
            c%digit(n+1) = carry
        end if
    end function big_add

    pure function big_product(a,b) result(c)
        type(positive_integer), intent(in) :: a !! first nonnegative integer
        type(positive_integer), intent(in) :: b !! second nonnegative integer
        type(positive_integer) :: c !! exact product
        integer :: na, nb, i, j, last
        integer(int64) :: carry, t

        if (.not. allocated(a%digit) .or. .not. allocated(b%digit)) return
        na = size(a%digit)
        nb = size(b%digit)
        allocate(c%digit(na+nb))
        c%digit = 0_int64
        do i=1,na
            carry = 0_int64
            do j=1,nb
                t = c%digit(i+j-1)+a%digit(i)*b%digit(j)+carry
                c%digit(i+j-1) = mod(t,limb_base)
                carry = t/limb_base
            end do
            c%digit(i+nb) = carry
        end do
        last = size(c%digit)
        do while (last > 0)
            if (c%digit(last) /= 0_int64) exit
            last = last-1
        end do
        if (last == 0) then
            deallocate(c%digit)
        else if (last < size(c%digit)) then
            c%digit = c%digit(:last)
        end if
    end function big_product

    function big_text(a) result(value)
        type(positive_integer), intent(in) :: a !! nonnegative arbitrary-precision integer
        character(:), allocatable :: value !! exact decimal representation
        character(32) :: buffer
        integer :: i, n

        value = '0'
        if (.not. allocated(a%digit)) return
        n = size(a%digit)
        if (n == 0) return
        write(buffer,'(i0)') a%digit(n)
        value = trim(buffer)
        do i=n-1,1,-1
            write(buffer,'(i9.9)') a%digit(i)
            value = value//buffer(:9)
        end do
    end function big_text

end module scifort_linalg_pascal
