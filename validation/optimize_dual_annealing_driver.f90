! SPDX-License-Identifier: MIT
program optimize_dual_annealing_driver
    use scifort_kinds, only : dp
    use scifort_optimize_dual_annealing, only : dual_annealing_core, dual_annealing_result
    implicit none
    integer, parameter :: ncases = 8, nseeds = 6
    integer, parameter :: iterations = 360
    integer :: icase, iseed, n
    real(dp), allocatable :: bounds(:,:)
    type(dual_annealing_result) :: res

    do icase = 1, ncases
        n = 2
        if (icase == 7) n = 3
        if (icase == 8) n = 1
        allocate(bounds(n,2))
        select case (icase)
        case (1,7)
            bounds(:,1) = -5.0_dp
            bounds(:,2) = 5.0_dp
        case (2)
            bounds(:,1) = -5.0_dp
            bounds(:,2) = 5.0_dp
        case (3)
            bounds(:,1) = -5.12_dp
            bounds(:,2) = 5.12_dp
        case (4)
            bounds(:,1) = -3.0_dp
            bounds(:,2) = 3.0_dp
        case (5)
            bounds(:,1) = -8.0_dp
            bounds(:,2) = 8.0_dp
        case (6)
            bounds(:,1) = -6.0_dp
            bounds(:,2) = 6.0_dp
        case (8)
            bounds(:,1) = -4.0_dp
            bounds(:,2) = 4.0_dp
        end select
        do iseed = 1, nseeds
            res = dual_annealing_core(objective, bounds, maxiter=iterations, seed=101 + 127*iseed + icase, &
                no_local_search=.true.)
            print '(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,es24.16,*(1x,es24.16))', &
                'CASE', icase, iseed, res%status, res%nfev, res%fun, res%x
        end do
        deallocate(bounds)
    end do
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! candidate point, shape (n)
        real(dp) :: f !! objective to minimize
        real(dp) :: t
        integer :: j

        select case (icase)
        case (1)
            f = sum((x - [1.2_dp,-0.7_dp])**2)
        case (2)
            f = -20.0_dp*exp(-0.2_dp*sqrt(sum(x**2)/2.0_dp)) - &
                exp(sum(cos(2.0_dp*acos(-1.0_dp)*x))/2.0_dp) + 20.0_dp + exp(1.0_dp)
        case (3)
            f = 20.0_dp + sum(x**2 - 10.0_dp*cos(2.0_dp*acos(-1.0_dp)*x))
        case (4)
            f = 100.0_dp*(x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
        case (5)
            t = 1.0_dp
            do j = 1, size(x)
                t = t*cos(x(j)/sqrt(real(j,dp)))
            end do
            f = 1.0_dp + sum(x**2)/4000.0_dp - t
        case (6)
            f = (x(1)**2 + x(2) - 11.0_dp)**2 + (x(1) + x(2)**2 - 7.0_dp)**2
        case (7)
            f = sum((x - [0.5_dp,-1.0_dp,2.0_dp])**2)
        case (8)
            f = (x(1) - 0.25_dp)**2 + 0.3_dp*(1.0_dp - cos(8.0_dp*x(1) - 2.0_dp))
        end select
    end function objective
end program optimize_dual_annealing_driver
