! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Simple streaming reference driver used by tools/compare_linalg_special_matrices.py.
! Input: case count, followed by one line per case: mode n bad_shape.
! Mode M=helmert compact, F=helmert full, H=hilbert,
! I=invhilbert real, E=invhilbert int64, D=invhilbert exact decimal.
! Output per case: status on first line, then one row per line.
program linalg_special_matrices_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_helmert, only : helmert
    use scifort_linalg_hilbert, only : hilbert, invhilbert, invhilbert_exact_decimal
    implicit none
    character(1) :: mode
    integer :: count, n, bad, nr, nc, s, t, i, j, ios
    real(dp), allocatable :: real_matrix(:,:)
    integer(int64), allocatable :: integer_matrix(:,:)
    character(256), allocatable :: decimal_matrix(:,:)

    read(*,*,iostat=ios) count
    if (ios /= 0 .or. count < 0) error stop 'bad reference case count'
    do t=1,count
        read(*,*,iostat=ios) mode,n,bad
        if (ios /= 0 .or. bad < 0 .or. bad > 1) error stop 'bad reference case'
        nr = max(0,n)
        nc = max(0,n)
        if (mode == 'M') nr = max(0,n-1)
        nr = nr+bad
        select case(mode)
        case ('M','F','H','I')
            allocate(real_matrix(nr,nc))
            select case(mode)
            case ('M')
                call helmert(n,real_matrix,s)
            case ('F')
                call helmert(n,real_matrix,s,full=.true.)
            case ('H')
                call hilbert(n,real_matrix,s)
            case ('I')
                call invhilbert(n,real_matrix,s)
            end select
            print '(i0)', s
            do i=1,nr
                write(*,'(*(es26.17e3,1x))') real_matrix(i,:)
            end do
            deallocate(real_matrix)
        case ('E')
            allocate(integer_matrix(nr,nc))
            call invhilbert(n,integer_matrix,s)
            print '(i0)', s
            do i=1,nr
                write(*,'(*(i0,1x))') integer_matrix(i,:)
            end do
            deallocate(integer_matrix)
        case ('D')
            allocate(decimal_matrix(nr,nc))
            call invhilbert_exact_decimal(n,decimal_matrix,s)
            print '(i0)', s
            do i=1,nr
                write(*,'(*(a,1x))') (trim(decimal_matrix(i,j)),j=1,nc)
            end do
            deallocate(decimal_matrix)
        case default
            error stop 'invalid matrix mode'
        end select
    end do
end program linalg_special_matrices_driver
