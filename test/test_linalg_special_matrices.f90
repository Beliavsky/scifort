! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_linalg_special_matrices
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_helmert, only : helmert
    use scifort_linalg_hilbert, only : hilbert, invhilbert, invhilbert_exact_decimal
    implicit none
    real(dp), allocatable :: a(:,:), b(:,:), t(:,:), gram(:,:)
    integer(int64), allocatable :: x(:,:)
    character(96), allocatable :: decimal(:,:)
    integer :: n, i, j, s
    real(dp) :: bound

    ! Direct agreement with SciPy examples, including the row orientation.
    block
        real(dp) :: h(5,5), compact(4,5)
        call helmert(5,h,s,full=.true.)
        call ensure(s == linalg_status_success,'Helmert full, size five')
        call ensure(maxval(abs(h(1,:)-1.0_dp/sqrt(5.0_dp))) < 1.0e-15_dp, &
            'Helmert normalized mean row')
        call ensure(abs(h(2,1)-1.0_dp/sqrt(2.0_dp)) < 1.0e-15_dp, 'Helmert positive contrast')
        call ensure(abs(h(2,2)+1.0_dp/sqrt(2.0_dp)) < 1.0e-15_dp, 'Helmert negative contrast')
        call helmert(5,compact,s)
        call ensure(s == linalg_status_success,'Helmert default compact')
        call ensure(all(compact == h(2:,:)),'Helmert compact equals full rows 2..n')
    end block

    ! Orthonormal rows and zero-mean contrasts over a wide range of orders.
    do n=1,96
        allocate(a(n,n),b(n-1,n),gram(n,n))
        call helmert(n,a,s,full=.true.)
        call ensure(s == linalg_status_success,'Helmert full, variable size')
        call helmert(n,b,s)
        call ensure(s == linalg_status_success,'Helmert compact, variable size')
        call ensure(all(a(2:n,:) == b),'Helmert full/compact identity')
        gram = matmul(a,transpose(a))
        do i=1,n
            gram(i,i) = gram(i,i)-1.0_dp
        end do
        bound = 6.0e-14_dp*sqrt(real(n,dp))
        call ensure(maxval(abs(gram)) < bound,'Helmert orthonormal rows')
        if (n > 1) call ensure(maxval(abs(sum(b,dim=2))) < bound, &
            'Helmert compact contrasts have zero sums')
        deallocate(a,b,gram)
    end do

    block
        real(dp) :: one(1,1), compact(0,1), zero(0,0), wrong(2,3)
        call helmert(1,one,s,full=.true.)
        call ensure(s == 0 .and. one(1,1) == 1.0_dp,'Helmert n=1 full')
        call helmert(1,compact,s)
        call ensure(s == 0,'Helmert n=1 compact')
        call helmert(0,zero,s,full=.true.)
        call ensure(s == linalg_status_invalid_shape,'Helmert zero order rejection')
        call helmert(-1,zero,s)
        call ensure(s == linalg_status_invalid_shape,'Helmert negative order rejection')
        wrong = 9.0_dp
        call helmert(3,wrong,s,full=.true.)
        call ensure(s == linalg_status_invalid_shape .and. all(wrong == 0.0_dp), &
            'Helmert invalid output shape cleared')
    end block

    block
        real(dp) :: h(4,4), h_inv(4,4), id(4,4)
        integer(int64) :: exact(4,4)
        call hilbert(4,h,s)
        call ensure(s == 0,'Hilbert order four')
        call ensure(h(1,1) == 1.0_dp .and. h(1,2) == 0.5_dp, &
            'Hilbert first entries')
        call invhilbert(4,h_inv,s)
        call ensure(s == 0,'real inverse Hilbert, order four')
        call invhilbert(4,exact,s)
        call ensure(s == 0,'exact inverse Hilbert, order four')
        call ensure(all(exact(1,:) == [16_int64,-120_int64,240_int64,-140_int64]), &
            'inverse Hilbert first row')
        call ensure(all(exact(2,:) == [-120_int64,1200_int64,-2700_int64,1680_int64]), &
            'inverse Hilbert second row')
        call ensure(maxval(abs(h_inv-real(exact,dp))) < 1.0e-10_dp, &
            'real and exact inverse Hilbert agree')
        id = matmul(h,h_inv)
        do i=1,4
            id(i,i) = id(i,i)-1.0_dp
        end do
        call ensure(maxval(abs(id)) < 1.0e-10_dp,'Hilbert times inverse is identity')
    end block

    do n=0,14
        allocate(a(n,n),b(n,n),x(n,n))
        call hilbert(n,a,s)
        call ensure(s == linalg_status_success,'Hilbert variable size')
        call invhilbert(n,b,s)
        call ensure(s == linalg_status_success,'inverse Hilbert real variable size')
        call invhilbert(n,x,s)
        call ensure(s == linalg_status_success,'inverse Hilbert int64 variable size')
        if (n > 0) then
            do i=1,n
                do j=1,n
                    call ensure(a(i,j) == 1.0_dp/(real(i,dp)+real(j,dp)-1.0_dp), &
                        'Hilbert defining entries')
                    call ensure(x(i,j) == x(j,i),'exact inverse Hilbert symmetry')
                    call ensure(abs(b(i,j)-real(x(i,j),dp)) <= &
                        8.0e-15_dp*max(1.0_dp,abs(real(x(i,j),dp))), &
                        'exact and floating inverse Hilbert agreement')
                end do
            end do
            call ensure(all(a == transpose(a)),'Hilbert symmetry')
            if (n <= 8) then
                allocate(t(n,n))
                t = matmul(a,b)
                do i=1,n
                    t(i,i) = t(i,i)-1.0_dp
                end do
                call ensure(maxval(abs(t)) < 3.0e-7_dp,'inverse Hilbert reconstruction')
                deallocate(t)
            end if
        end if
        deallocate(a,b,x)
    end do

    block
        integer(int64) :: e14(14,14), too_large(15,15)
        real(dp) :: empty(0,0), wrong(3,2)
        real(dp), allocatable :: overflow(:,:)
        call invhilbert(14,e14,s)
        call ensure(s == 0,'int64 inverse Hilbert 14')
        call ensure(maxval(abs(e14)) == 3521767173114190000_int64,'int64 inverse Hilbert maximum')
        too_large = 7_int64
        call invhilbert(15,too_large,s)
        call ensure(s == linalg_status_invalid_shape .and. all(too_large == 0_int64), &
            'int64 inverse Hilbert 15 rejected instead of overflowing')
        call hilbert(0,empty,s)
        call ensure(s == linalg_status_success,'empty Hilbert matrix')
        call invhilbert(0,empty,s)
        call ensure(s == linalg_status_success,'empty inverse Hilbert')
        call hilbert(-1,empty,s)
        call ensure(s == linalg_status_invalid_shape,'negative Hilbert rejected')
        call invhilbert(-1,empty,s)
        call ensure(s == linalg_status_invalid_shape,'negative inverse Hilbert rejected')
        wrong = 3.0_dp
        call hilbert(3,wrong,s)
        call ensure(s == linalg_status_invalid_shape .and. all(wrong == 0.0_dp), &
            'wrong Hilbert output dimensions rejected')
        wrong = 3.0_dp
        call invhilbert(3,wrong,s)
        call ensure(s == linalg_status_invalid_shape .and. all(wrong == 0.0_dp), &
            'wrong inverse Hilbert output dimensions rejected')
        allocate(overflow(250,250))
        call invhilbert(250,overflow,s)
        call ensure(s == linalg_status_no_convergence .and. all(overflow == 0.0_dp), &
            'inverse Hilbert floating overflow is checked')
    end block

    ! Arbitrary-precision decimal path agrees with small int64 values and
    ! produces known values beyond int64 limits.
    do n=0,30
        allocate(decimal(n,n))
        call invhilbert_exact_decimal(n,decimal,s)
        call ensure(s == linalg_status_success,'decimal inverse Hilbert variable size')
        if (n > 0) then
            call ensure(trim(decimal(1,1)) == int_text(int(n,int64)**2), &
                'exact decimal first entry')
            do i=1,n
                do j=1,n
                    call ensure(decimal(i,j) == decimal(j,i),'exact decimal symmetry')
                end do
            end do
        end if
        if (n <= 14) then
            allocate(x(n,n))
            call invhilbert(n,x,s)
            do i=1,n
                do j=1,n
                    call ensure(trim(decimal(i,j)) == int_text(x(i,j)), &
                        'decimal and int64 exact agreement')
                end do
            end do
            deallocate(x)
        end if
        deallocate(decimal)
    end do

    block
        character(90) :: big(16,16), empty(0,0)
        character(8) :: too_short(16,16)
        call invhilbert_exact_decimal(16,big,s)
        call ensure(s == 0,'large exact decimal inverse')
        call ensure(trim(big(8,8)) == '42475099528537378560', &
            'known SciPy exact inverse Hilbert (16, 8, 8)')
        call invhilbert_exact_decimal(16,too_short,s)
        call ensure(s == linalg_status_invalid_shape .and. all(too_short == ''), &
            'short decimal output is rejected atomically')
        call invhilbert_exact_decimal(-1,empty,s)
        call ensure(s == linalg_status_invalid_shape,'negative decimal inverse rejected')
    end block

    print '(a)', 'test_linalg_special_matrices: PASS'
contains
    subroutine ensure(condition, message)
        logical, intent(in) :: condition !! condition that must hold
        character(*), intent(in) :: message !! diagnostic label
        if (.not. condition) then
            print '(a,1x,a)', 'FAILED:',message
            error stop 1
        end if
    end subroutine ensure

    function int_text(value) result(buffer)
        integer(int64), intent(in) :: value !! signed integer to format
        character(48) :: buffer !! left-justified decimal representation
        write(buffer,'(i0)') value
    end function int_text
end program test_linalg_special_matrices
