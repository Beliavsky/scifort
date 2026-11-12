! SPDX-License-Identifier: MIT
program driver_real_schur
    use scifort_kinds, only : dp
    use scifort_linalg_schur, only : schur
    implicit none

    real(dp), allocatable :: a(:, :), t(:, :), z(:, :)
    integer :: cases, testcase, n, rule_code, i, j, status, sdim
    character(len=4) :: rule

    read(*,*) cases
    do testcase = 1, cases
        read(*,*) n, rule_code
        allocate(a(n,n), t(n,n), z(n,n))
        do j = 1,n
            do i = 1,n
                read(*,*) a(i,j)
            end do
        end do
        select case (rule_code)
        case (0)
            call schur(a,t,z,status,sdim=sdim)
        case (1:4)
            select case (rule_code)
            case (1)
                rule = 'lhp'
            case (2)
                rule = 'rhp'
            case (3)
                rule = 'iuc'
            case default
                rule = 'ouc'
            end select
            call schur(a,t,z,status,sort=trim(rule),sdim=sdim)
        case (5)
            call schur(a,t,z,status,selector=positive_imag,sdim=sdim)
        case default
            status = -1
            sdim = 0
        end select
        write(*,'(i0,1x,i0)') status,sdim
        if (status == 0) then
            do j = 1,n
                do i = 1,n
                    write(*,'(2(es26.17e3,1x))') t(i,j),z(i,j)
                end do
            end do
        end if
        deallocate(a,t,z)
    end do
contains
    logical function positive_imag(w)
        complex(dp), intent(in) :: w !! Schur eigenvalue
        positive_imag = aimag(w) > 0.0_dp
    end function positive_imag
end program driver_real_schur
