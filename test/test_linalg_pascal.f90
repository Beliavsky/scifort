! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_linalg_pascal
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_pascal, only : pascal, invpascal, pascal_exact_decimal, invpascal_exact_decimal
    implicit none
    character(9), parameter :: kinds(3) = [character(9) :: 'symmetric', 'lower', 'upper']
    integer(int64), allocatable :: p(:,:), q(:,:), identity(:,:)
    real(dp), allocatable :: a(:,:), b(:,:)
    character(180), allocatable :: str_p(:,:), str_q(:,:)
    integer :: n, which, i, j, s, ios
    integer(int64) :: exact_value
    real(dp) :: tolerance

    block
        integer(int64) :: m(4,4), inv(5,5)
        call pascal(4,m,s)
        call check(s == linalg_status_success,'Pascal published example')
        call check(all(m(1,:) == 1_int64),'Pascal symmetric first row')
        call check(all(m(2,:) == [1_int64,2_int64,3_int64,4_int64]),'Pascal symmetric second row')
        call check(all(m(4,:) == [1_int64,4_int64,10_int64,20_int64]),'Pascal symmetric last row')
        call invpascal(5,inv,s)
        call check(s == linalg_status_success,'inverse Pascal published example')
        call check(all(inv(1,:) == [5_int64,-10_int64,10_int64,-5_int64,1_int64]), &
            'inverse symmetric Pascal first row')
        call check(all(inv(3,:) == [10_int64,-35_int64,46_int64,-27_int64,6_int64]), &
            'inverse symmetric Pascal middle row')
    end block

    do which=1,3
        do n=0,30
            allocate(p(n,n),q(n,n),a(n,n),b(n,n),str_p(n,n),str_q(n,n))
            call pascal(n,p,s,kind=trim(kinds(which)))
            call check(s == 0,'Pascal int64 status')
            call invpascal(n,q,s,kind=trim(kinds(which)))
            call check(s == 0,'inverse Pascal int64 status')
            call pascal(n,a,s,kind=trim(kinds(which)))
            call check(s == 0,'Pascal floating status')
            call invpascal(n,b,s,kind=trim(kinds(which)))
            call check(s == 0,'inverse Pascal floating status')
            call pascal_exact_decimal(n,str_p,s,kind=trim(kinds(which)))
            call check(s == 0,'Pascal decimal status')
            call invpascal_exact_decimal(n,str_q,s,kind=trim(kinds(which)))
            call check(s == 0,'inverse Pascal decimal status')
            do j=1,n
                do i=1,n
                    read(str_p(i,j),*,iostat=ios) exact_value
                    call check(ios == 0,'Pascal decimal parses as int64')
                    call check(exact_value == p(i,j),'Pascal decimal agrees exactly with int64')
                    read(str_q(i,j),*,iostat=ios) exact_value
                    call check(ios == 0,'inverse Pascal decimal parses as int64')
                    call check(exact_value == q(i,j),'inverse Pascal decimal agrees exactly with int64')
                    tolerance = 3.0e-13_dp
                    call check(abs(a(i,j)-real(p(i,j),dp)) <= tolerance*max(1.0_dp,abs(a(i,j))), &
                        'Pascal floating agrees with integer')
                    call check(abs(b(i,j)-real(q(i,j),dp)) <= tolerance*max(1.0_dp,abs(b(i,j))), &
                        'inverse Pascal floating agrees with integer')
                end do
            end do
            if (n >= 1 .and. n <= 9) then
                allocate(identity(n,n))
                identity = matmul(p,q)
                do i=1,n
                    identity(i,i) = identity(i,i)-1_int64
                end do
                call check(all(identity == 0_int64),'Pascal times inverse exact identity')
                deallocate(identity)
            end if
            if (n > 0) then
                if (which == 1) then
                    call check(all(p == transpose(p)) .and. all(q == transpose(q)), &
                        'symmetric Pascal entries')
                else if (which == 2) then
                    do j=1,n
                        do i=1,j-1
                            call check(p(i,j) == 0_int64 .and. q(i,j) == 0_int64, 'lower triangle')
                        end do
                    end do
                else
                    do i=1,n
                        do j=1,i-1
                            call check(p(i,j) == 0_int64 .and. q(i,j) == 0_int64, 'upper triangle')
                        end do
                    end do
                end if
            end if
            deallocate(p,q,a,b,str_p,str_q)
        end do
    end do

    block
        integer(int64) :: m34(34,34), m35(35,35), m36(36,36), m67(67,67), m68(68,68)
        character(100) :: big(50,50)
        character(7) :: short(50,50)
        character(3) :: lower(5,5)
        character(40) :: empty(0,0)
        real(dp) :: invalid(2,3), zero(0,0)

        call pascal(34,m34,s)
        call check(s == 0 .and. m34(34,34) == 7219428434016265740_int64, &
            'last int64-representable symmetric Pascal order')
        m35 = -1_int64
        call pascal(35,m35,s)
        call check(s == linalg_status_no_convergence .and. all(m35 == 0_int64), &
            'symmetric Pascal int64 overflow clears output')
        call invpascal(35,m35,s)
        call check(s == 0,'inverse symmetric Pascal supports n=35 in int64')
        call invpascal(36,m36,s)
        call check(s == linalg_status_no_convergence .and. all(m36 == 0_int64), &
            'inverse symmetric Pascal int64 overflow clears output')
        call pascal(67,m67,s,kind='lower')
        call check(s == 0,'lower Pascal order 67 fits int64')
        call pascal(68,m68,s,kind='lower')
        call check(s == linalg_status_no_convergence .and. all(m68 == 0_int64), &
            'lower Pascal order 68 triggers int64 overflow')
        call pascal_exact_decimal(50,big,s)
        call check(s == 0,'big Pascal decimal order 50')
        call check(trim(big(50,50)) == '25477612258980856902730428600', &
            'SciPy documented exact Pascal order 50 entry')
        call invpascal_exact_decimal(50,big,s)
        call check(s == 0,'big inverse Pascal decimal order 50')
        call check(trim(big(1,50)) == '-1','inverse symmetric Pascal upper corner sign')
        call pascal_exact_decimal(50,short,s)
        call check(s == linalg_status_invalid_shape .and. all(short == ''), &
            'short decimal output for Pascal rejected atomically')
        call invpascal_exact_decimal(50,short,s)
        call check(s == linalg_status_invalid_shape .and. all(short == ''), &
            'short decimal output for inverse Pascal rejected atomically')
        call invpascal_exact_decimal(5,lower,s,kind='lower')
        call check(s == 0 .and. lower(5,1) == '1' .and. lower(5,2) == '-4', &
            'signed lower inverse decimal')
        call pascal(0,zero,s)
        call check(s == 0,'zero-order real Pascal')
        call invpascal(0,zero,s)
        call check(s == 0,'zero-order real inverse Pascal')
        call pascal_exact_decimal(0,empty,s)
        call check(s == 0,'zero-order decimal Pascal')
        call invpascal_exact_decimal(0,empty,s)
        call check(s == 0,'zero-order decimal inverse Pascal')
        call pascal(-1,zero,s)
        call check(s == linalg_status_invalid_shape,'negative Pascal order rejected')
        call invpascal(-1,zero,s)
        call check(s == linalg_status_invalid_shape,'negative inverse Pascal order rejected')
        invalid = 5.0_dp
        call pascal(3,invalid,s)
        call check(s == linalg_status_invalid_shape .and. all(invalid == 0.0_dp), &
            'bad Pascal shape rejected')
        invalid = 5.0_dp
        call invpascal(3,invalid,s)
        call check(s == linalg_status_invalid_shape .and. all(invalid == 0.0_dp), &
            'bad inverse Pascal shape rejected')
        m34 = -1_int64
        call pascal(34,m34,s,kind='diagonal')
        call check(s == linalg_status_invalid_shape .and. all(m34 == 0_int64), &
            'invalid Pascal kind rejected')
    end block

    print '(a)', 'test_linalg_pascal: PASS'
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition !! test assertion
        character(*), intent(in) :: label !! diagnostic if false
        if (.not. condition) then
            print '(a,1x,a)', 'FAILED:',trim(label)
            error stop 1
        end if
    end subroutine check
end program test_linalg_pascal
