! SPDX-License-Identifier: MIT
! Standalone batch driver for independent numerical comparisons.
program optimize_milp_driver
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_milp, only : milp, milp_result
    implicit none
    type(milp_result) :: res
    real(dp), allocatable :: c(:), lower(:), upper(:), a(:, :), b_l(:), b_u(:)
    integer, allocatable :: integrality(:)
    integer :: n, m, ncases, case_number, i
    real(dp) :: inf
    inf = ieee_value(1.0_dp, ieee_positive_inf)
    read(*,*) ncases
    do case_number = 1, ncases
        read(*,*) n,m
        allocate(c(n), integrality(n), lower(n), upper(n), a(m,n), b_l(m), b_u(m))
        read(*,*) c
        read(*,*) integrality
        read(*,*) lower
        read(*,*) upper
        if (m>0) then
            do i=1,m
                read(*,*) a(i,:)
            end do
            read(*,*) b_l
            read(*,*) b_u
        end if
        where (lower <= -1.0e99_dp) lower = -inf
        where (upper >= 1.0e99_dp) upper = inf
        where (b_l <= -1.0e99_dp) b_l = -inf
        where (b_u >= 1.0e99_dp) b_u = inf
        res = milp(c, integrality=integrality, a=a, b_l=b_l, b_u=b_u, &
            lower=lower, upper=upper, mip_rel_gap=0.0_dp, node_limit=3000)
        write(*,'(i0,1x,i0,1x,es28.17e3,1x,i0,1x,es28.17e3,1x,es28.17e3)') &
            case_number, res%status, res%fun, res%mip_node_count, &
            res%mip_dual_bound, res%mip_gap
        if (allocated(res%x)) then
            write(*,'(*(es28.17e3,1x))') res%x
        else
            write(*,'(a)') 'NONE'
        end if
        deallocate(c,integrality,lower,upper,a,b_l,b_u)
    end do
end program optimize_milp_driver
