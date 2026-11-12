! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Hilbert and inverse Hilbert matrix constructors, with real64, exact int64,
! and arbitrary-precision decimal-string results. Exact decimal strings are
! provided because standard Fortran lacks a portable arbitrary-width integer.
module scifort_linalg_hilbert
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    implicit none
    private

    integer(int64), parameter :: base = 1000000000_int64
    public :: hilbert, invhilbert, invhilbert_exact_decimal

    interface invhilbert
        module procedure invhilbert_real
        module procedure invhilbert_int64
    end interface invhilbert

contains

    pure subroutine hilbert(n, h, status)
        integer, intent(in) :: n !! order, nonnegative
        real(dp), intent(out) :: h(:,:) !! Hilbert matrix, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        h = 0.0_dp
        status = linalg_status_invalid_shape
        if (n < 0) return
        if (size(h,1) /= n .or. size(h,2) /= n) return
        do j=1,n
            do i=1,n
                h(i,j) = 1.0_dp/(real(i,dp)+real(j,dp)-1.0_dp)
            end do
        end do
        status = linalg_status_success
    end subroutine hilbert

    ! Integer Hilbert inverse formula, with one-based indices:
    ! (-1)^(i+j) (i+j-1) C(n+i-1,n-j) C(n+j-1,n-i) C(i+j-2,i-1)^2.
    ! All factors are positive except the checkerboard sign.
    pure subroutine invhilbert_real(n, h, status)
        integer, intent(in) :: n !! order, nonnegative
        real(dp), intent(out) :: h(:,:) !! floating inverse Hilbert matrix, shape (n,n)
        integer, intent(out) :: status !! success, invalid_shape, or no_convergence on overflow
        integer :: i, j, t
        real(dp) :: factors(4), value
        logical :: ok

        h = 0.0_dp
        status = linalg_status_invalid_shape
        if (n < 0 .or. size(h,1) /= n .or. size(h,2) /= n) return
        status = linalg_status_success
        do i=1,n
            do j=1,i
                call choose_real(n+i-1,n-j,factors(1),ok)
                if (.not. ok) exit
                call choose_real(n+j-1,n-i,factors(2),ok)
                if (.not. ok) exit
                call choose_real(i+j-2,i-1,factors(3),ok)
                if (.not. ok) exit
                factors(4) = factors(3)
                value = real(i,dp)+real(j,dp)-1.0_dp
                do t=1,4
                    if (factors(t) > 0.0_dp) then
                        if (value > huge(value)/factors(t)) exit
                    end if
                    value = value*factors(t)
                end do
                if (t <= 4) exit
                if (mod(i+j,2) /= 0) value = -value
                h(i,j) = value
                h(j,i) = value
            end do
            if (j <= i) exit
        end do
        if (i <= n) then
            h = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine invhilbert_real

    pure subroutine choose_real(n, k, v, ok)
        integer, intent(in) :: n !! upper argument
        integer, intent(in) :: k !! lower argument
        real(dp), intent(out) :: v !! real binomial coefficient
        logical, intent(out) :: ok !! false if overflow would occur
        integer :: t, kk
        real(dp) :: factor

        v = 1.0_dp
        ok = .true.
        if (k < 0 .or. k > n) then
            v = 0.0_dp
            return
        end if
        kk = min(k,n-k)
        do t=1,kk
            factor = real(n-kk+t,dp)/real(t,dp)
            if (v > huge(v)/factor) then
                ok = .false.
                return
            end if
            v = v*factor
        end do
    end subroutine choose_real

    pure subroutine invhilbert_int64(n, h, status)
        integer, intent(in) :: n !! order, 0 <= n <= 14
        integer(int64), intent(out) :: h(:,:) !! exact inverse, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape (including n>14)
        integer :: i, j
        integer(int64) :: value, c

        h = 0_int64
        status = linalg_status_invalid_shape
        if (n < 0 .or. n > 14) return
        if (size(h,1) /= n .or. size(h,2) /= n) return
        do i=1,n
            do j=1,i
                value = int(i+j-1,int64)
                value = value*choose_int64(n+i-1,n-j)
                value = value*choose_int64(n+j-1,n-i)
                c = choose_int64(i+j-2,i-1)
                value = value*c*c
                if (mod(i+j,2) /= 0) value = -value
                h(i,j) = value
                h(j,i) = value
            end do
        end do
        status = linalg_status_success
    end subroutine invhilbert_int64

    pure integer(int64) function choose_int64(n, k) result(v)
        integer, intent(in) :: n !! upper argument, n <= 27 for invhilbert exact int64
        integer, intent(in) :: k !! lower argument
        integer :: t, kk

        v = 0_int64
        if (k < 0 .or. k > n) return
        kk = min(k,n-k)
        v = 1_int64
        do t=1,kk
            v = (v*int(n-kk+t,int64))/int(t,int64)
        end do
    end function choose_int64

    subroutine invhilbert_exact_decimal(n, h, status)
        integer, intent(in) :: n !! order, nonnegative; exact results for n>14
        character(*), intent(out) :: h(:,:) !! decimal representations of exact integers, shape (n,n)
        integer, intent(out) :: status !! success or invalid_shape (including too-short strings)
        integer :: i, j, p, m, used, exponent, t, k, capacity
        integer(int64), allocatable :: digits(:)
        logical, allocatable :: sieve(:)
        character(:), allocatable :: value
        logical :: ok

        h = ''
        status = linalg_status_invalid_shape
        if (n < 0 .or. size(h,1) /= n .or. size(h,2) /= n) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        if (len(h) < 1 .or. n > (huge(n)-32)/2) return
        m = 2*n
        ! At most 8*n decimal digits is a loose bound for the products
        ! of three binomial factors. This allocates more than needed.
        capacity = n+16
        allocate(digits(capacity),sieve(0:m))
        sieve = .true.
        sieve(0:1) = .false.
        do p=2,m
            if (p > m/p) exit
            if (.not. sieve(p)) cycle
            do t=p*p,m,p
                sieve(t) = .false.
            end do
        end do

        do i=1,n
            do j=1,i
                digits = 0_int64
                digits(1) = 1_int64
                used = 1
                ok = .true.
                do p=2,m
                    if (.not. sieve(p)) cycle
                    exponent = v_binomial(n+i-1,n-j,p) + &
                        v_binomial(n+j-1,n-i,p) + &
                        2*v_binomial(i+j-2,i-1,p) + v_number(i+j-1,p)
                    do k=1,exponent
                        call multiply_small(digits,used,int(p,int64),ok)
                        if (.not. ok) exit
                    end do
                    if (.not. ok) exit
                end do
                if (.not. ok) then
                    h = ''
                    return
                end if
                value = decimal_string(digits,used)
                if (mod(i+j,2) /= 0) value = '-'//value
                if (len(value) > len(h)) then
                    h = ''
                    return
                end if
                h(i,j) = value
                h(j,i) = value
            end do
        end do
        status = linalg_status_success
    end subroutine invhilbert_exact_decimal

    pure integer function v_factorial(n,p) result(v)
        integer, intent(in) :: n !! nonnegative factorial argument
        integer, intent(in) :: p !! prime
        integer :: q

        v = 0
        q = n
        do while (q >= p)
            q = q/p
            v = v+q
        end do
    end function v_factorial

    pure integer function v_binomial(n,k,p) result(v)
        integer, intent(in) :: n !! upper binomial argument
        integer, intent(in) :: k !! lower binomial argument
        integer, intent(in) :: p !! prime

        v = v_factorial(n,p)-v_factorial(k,p)-v_factorial(n-k,p)
    end function v_binomial

    pure integer function v_number(n,p) result(v)
        integer, intent(in) :: n !! positive integer
        integer, intent(in) :: p !! prime
        integer :: q

        q = n
        v = 0
        do while (mod(q,p) == 0 .and. q > 0)
            v = v+1
            q = q/p
        end do
    end function v_number

    pure subroutine multiply_small(digits,used,multiplier,ok)
        integer(int64), intent(inout) :: digits(:) !! little-endian base-10**9 limbs
        integer, intent(inout) :: used !! count of occupied limbs
        integer(int64), intent(in) :: multiplier !! positive integer factor <= 2*n
        logical, intent(out) :: ok !! false on inadequate limb capacity
        integer :: k
        integer(int64) :: carry, temporary

        ok = .true.
        carry = 0_int64
        do k=1,used
            temporary = digits(k)*multiplier+carry
            digits(k) = mod(temporary,base)
            carry = temporary/base
        end do
        do while (carry > 0_int64)
            if (used >= size(digits)) then
                ok = .false.
                return
            end if
            used = used+1
            digits(used) = mod(carry,base)
            carry = carry/base
        end do
    end subroutine multiply_small

    function decimal_string(digits,used) result(value)
        integer(int64), intent(in) :: digits(:) !! little-endian base-10**9 limbs
        integer, intent(in) :: used !! count of occupied limbs
        character(:), allocatable :: value !! exact unsigned decimal number
        character(32) :: buffer
        integer :: k

        write(buffer,'(i0)') digits(used)
        value = trim(buffer)
        do k=used-1,1,-1
            write(buffer,'(i9.9)') digits(k)
            value = value//trim(buffer)
        end do
    end function decimal_string

end module scifort_linalg_hilbert
