! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Streaming independent-reference driver for tools/compare_linalg_pascal.py.
! Input: number of cases, then (P|I) (R|E|D) (symmetric|lower|upper) order.
! Output: status line, then n rows of whitespace-delimited matrix entries
!         if status is zero. E=int64 exact, D=arbitrary exact decimal.
program linalg_pascal_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_pascal, only : pascal, invpascal, pascal_exact_decimal, invpascal_exact_decimal
    implicit none
    integer :: count, case_id, n, rows, status, i, j
    character(1) :: mode, dtype
    character(20) :: matrix_kind
    character(45) :: part
    character(:), allocatable :: row
    real(dp), allocatable :: r(:,:)
    integer(int64), allocatable :: e(:,:)
    character(250), allocatable :: d(:,:)

    read(*,*) count
    do case_id=1,count
        read(*,*) mode,dtype,matrix_kind,n
        rows = max(0,n)
        select case (dtype)
        case ('R')
            allocate(r(rows,rows))
            if (mode == 'P') then
                call pascal(n,r,status,kind=trim(matrix_kind))
            else
                call invpascal(n,r,status,kind=trim(matrix_kind))
            end if
            write(*,'(i0)') status
            if (status == 0) then
                do i=1,rows
                    write(*,'(*(es26.17e3,1x))') r(i,:)
                end do
            end if
            deallocate(r)
        case ('E')
            allocate(e(rows,rows))
            if (mode == 'P') then
                call pascal(n,e,status,kind=trim(matrix_kind))
            else
                call invpascal(n,e,status,kind=trim(matrix_kind))
            end if
            write(*,'(i0)') status
            if (status == 0) then
                do i=1,rows
                    row = ''
                    do j=1,rows
                        write(part,'(i0)') e(i,j)
                        row = row//trim(part)//' '
                    end do
                    write(*,'(a)') trim(row)
                end do
            end if
            deallocate(e)
        case ('D')
            allocate(d(rows,rows))
            if (mode == 'P') then
                call pascal_exact_decimal(n,d,status,kind=trim(matrix_kind))
            else
                call invpascal_exact_decimal(n,d,status,kind=trim(matrix_kind))
            end if
            write(*,'(i0)') status
            if (status == 0) then
                do i=1,rows
                    row = ''
                    do j=1,rows
                        row = row//trim(d(i,j))//' '
                    end do
                    write(*,'(a)') trim(row)
                end do
            end if
            deallocate(d)
        case default
            error stop 'unknown dtype in comparison input'
        end select
    end do
end program linalg_pascal_driver
