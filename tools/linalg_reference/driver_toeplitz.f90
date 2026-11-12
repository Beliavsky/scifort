! SPDX-License-Identifier: MIT
! Input/output driver used by run_toeplitz_checks.py, not a numerical oracle.
program driver_toeplitz
    use scifort_kinds, only : dp
    use scifort_linalg_toeplitz, only : solve_toeplitz
    implicit none
    real(dp), allocatable :: cr(:), rr(:), br(:, :), xr(:, :)
    complex(dp), allocatable :: cc(:), rc(:), bc(:, :), xc(:, :)
    integer :: n, nrhs, kind_id, has_row, ios, i, j, status
    logical :: dense

    do
        read(*,*,iostat=ios) n, nrhs, kind_id, has_row
        if (ios < 0) exit
        if (ios /= 0 .or. n < 0 .or. nrhs < 0) error stop 'malformed case header'
        select case (kind_id)
        case (1)
            allocate(cr(n),rr(n),br(n,nrhs),xr(n,nrhs))
            read(*,*) cr
            if (has_row == 1) read(*,*) rr
            read(*,*) br
            if (has_row == 1) then
                call solve_toeplitz(cr,br,xr,status,r=rr,used_dense=dense)
            else
                call solve_toeplitz(cr,br,xr,status,used_dense=dense)
            end if
            write(*,'(i0,1x,i0)') status, merge(1,0,dense)
            do j = 1, nrhs
                do i = 1, n
                    write(*,'(es26.17e3)') xr(i,j)
                end do
            end do
            deallocate(cr,rr,br,xr)
        case (2)
            allocate(cc(n),rc(n),bc(n,nrhs),xc(n,nrhs))
            read(*,*) cc
            if (has_row == 1) read(*,*) rc
            read(*,*) bc
            if (has_row == 1) then
                call solve_toeplitz(cc,bc,xc,status,r=rc,used_dense=dense)
            else
                call solve_toeplitz(cc,bc,xc,status,used_dense=dense)
            end if
            write(*,'(i0,1x,i0)') status, merge(1,0,dense)
            do j = 1, nrhs
                do i = 1, n
                    write(*,'(2(es26.17e3,1x))') real(xc(i,j),dp),aimag(xc(i,j))
                end do
            end do
            deallocate(cc,rc,bc,xc)
        case (3)
            allocate(cr(n),rr(n),bc(n,nrhs),xc(n,nrhs))
            read(*,*) cr
            if (has_row == 1) read(*,*) rr
            read(*,*) bc
            if (has_row == 1) then
                call solve_toeplitz(cr,bc,xc,status,r=rr,used_dense=dense)
            else
                call solve_toeplitz(cr,bc,xc,status,used_dense=dense)
            end if
            write(*,'(i0,1x,i0)') status, merge(1,0,dense)
            do j = 1, nrhs
                do i = 1, n
                    write(*,'(2(es26.17e3,1x))') real(xc(i,j),dp),aimag(xc(i,j))
                end do
            end do
            deallocate(cr,rr,bc,xc)
        case (4)
            allocate(cc(n),rc(n),br(n,nrhs),xc(n,nrhs))
            read(*,*) cc
            if (has_row == 1) read(*,*) rc
            read(*,*) br
            if (has_row == 1) then
                call solve_toeplitz(cc,br,xc,status,r=rc,used_dense=dense)
            else
                call solve_toeplitz(cc,br,xc,status,used_dense=dense)
            end if
            write(*,'(i0,1x,i0)') status, merge(1,0,dense)
            do j = 1, nrhs
                do i = 1, n
                    write(*,'(2(es26.17e3,1x))') real(xc(i,j),dp),aimag(xc(i,j))
                end do
            end do
            deallocate(cc,rc,br,xc)
        case default
            error stop 'bad case kind'
        end select
    end do
end program driver_toeplitz
