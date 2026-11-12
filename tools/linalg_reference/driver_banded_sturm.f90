! SPDX-License-Identifier: MIT
! Input/output protocol used by run_banded_sturm_checks.py.
program driver_banded_sturm
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eigvals_banded
    implicit none
    integer :: nc, icase, kind, n, kd, select_mode, first, last, capacity
    integer :: j, i, status, nfound
    real(dp) :: low, high, realpart, imagpart
    real(dp), allocatable :: w(:), br(:, :)
    complex(dp), allocatable :: bc(:, :)

    read (*,*) nc
    do icase = 1, nc
        read (*,*) kind, n, kd, select_mode, first, last, low, high, capacity
        allocate(w(capacity))
        if (kind <= 2) then
            allocate(br(kd+1,n))
            do j = 1, n
                do i = 1, kd+1
                    read (*,*) realpart, imagpart
                    br(i,j) = realpart
                end do
            end do
            select case(select_mode)
            case(0)
                call eigvals_banded(br,w,status,lower=(kind==2),method='sturm',num_selected=nfound)
            case(1)
                call eigvals_banded(br,w,status,lower=(kind==2), &
                    subset_by_index=[first,last],num_selected=nfound)
            case(2)
                call eigvals_banded(br,w,status,lower=(kind==2), &
                    subset_by_value=[low,high],num_selected=nfound)
            case(3)
                call eigvals_banded(br,w,status,lower=(kind==2),method='dense', &
                    subset_by_index=[first,last],num_selected=nfound)
            end select
            deallocate(br)
        else
            allocate(bc(kd+1,n))
            do j = 1, n
                do i = 1, kd+1
                    read (*,*) realpart,imagpart
                    bc(i,j) = cmplx(realpart,imagpart,dp)
                end do
            end do
            select case(select_mode)
            case(0)
                call eigvals_banded(bc,w,status,lower=(kind==4),method='sturm',num_selected=nfound)
            case(1)
                call eigvals_banded(bc,w,status,lower=(kind==4), &
                    subset_by_index=[first,last],num_selected=nfound)
            case(2)
                call eigvals_banded(bc,w,status,lower=(kind==4), &
                    subset_by_value=[low,high],num_selected=nfound)
            case(3)
                call eigvals_banded(bc,w,status,lower=(kind==4),method='dense', &
                    subset_by_index=[first,last],num_selected=nfound)
            end select
            deallocate(bc)
        end if
        write(*,'(I0,1X,I0)') status,nfound
        do j=1,capacity
            write(*,'(ES26.17E3)') w(j)
        end do
        deallocate(w)
    end do
end program driver_banded_sturm
