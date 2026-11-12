! SPDX-License-Identifier: MIT
! Non-interactive numerical reference driver for test/compare_optimize_differential_evolution_scipy.py.
! Input: case, strategy index, seed, update code (1 immediate; 2 deferred), init code (1 LHS; 2 random).
! Output: same identifiers, fun, x(1:3), nit, nfev, status.
program reference_differential_evolution_driver
    use scifort_kinds, only : dp
    use scifort_optimize_differential_evolution, only : differential_evolution, differential_evolution_result
    implicit none
    character(len=20), parameter :: strategies(12) = [character(len=20) :: &
        'best1bin', 'rand1bin', 'best2bin', 'rand2bin', 'randtobest1bin', 'currenttobest1bin', &
        'best1exp', 'rand1exp', 'best2exp', 'rand2exp', 'randtobest1exp', 'currenttobest1exp']
    type(differential_evolution_result) :: r
    real(dp), allocatable :: bounds(:,:)
    real(dp) :: xout(3)
    character(len=15) :: update_mode, init_mode
    integer :: which, strategy_id, seed, update_id, init_id, io, ndim
    do
        read(*,*,iostat=io) which,strategy_id,seed,update_id,init_id
        if (io /= 0) exit
        if (which < 1 .or. which > 6 .or. strategy_id < 1 .or. strategy_id > 12) error stop 'bad driver input'
        ndim = 2
        if (which == 6) ndim = 3
        allocate(bounds(ndim,2))
        select case(which)
        case(1)
            bounds(:,1) = -5.0_dp
            bounds(:,2) = 5.0_dp
        case(2)
            bounds(:,1) = -2.0_dp
            bounds(:,2) = 2.0_dp
        case(3)
            bounds(:,1) = -5.12_dp
            bounds(:,2) = 5.12_dp
        case(4)
            bounds(:,1) = -5.0_dp
            bounds(:,2) = 5.0_dp
        case(5)
            bounds(:,1) = -1.0_dp
            bounds(:,2) = [2.0_dp,1.0_dp]
        case(6)
            bounds(:,1) = [-5.0_dp,-3.0_dp,2.0_dp]
            bounds(:,2) = [5.0_dp,4.0_dp,2.0_dp]
        end select
        if (update_id == 1) then
            update_mode='immediate'
        else
            update_mode='deferred'
        end if
        if (init_id == 1) then
            init_mode='latinhypercube'
        else
            init_mode='random'
        end if
        r = differential_evolution(objective, bounds, strategy=trim(strategies(strategy_id)), &
            maxiter=170, popsize=9, seed=seed, tol=1.0e-3_dp, atol=1.0e-8_dp, &
            mutation_range=[0.5_dp,1.0_dp], recombination=0.7_dp, &
            init=trim(init_mode), updating=trim(update_mode))
        xout = 0.0_dp
        xout(1:ndim) = r%x
        print '(5(i0,1x),4(es24.16e3,1x),3(i0,1x))', &
            which,strategy_id,seed,update_id,init_id,r%fun,xout,r%nit,r%nfev,r%status
        deallocate(bounds)
    end do
contains
    real(dp) function objective(x) result(fun)
        real(dp), intent(in) :: x(:) !! candidate vector
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: a, b
        select case(which)
        case (1)
            fun = (x(1)-0.75_dp)**2 + (x(2)+1.1_dp)**2
        case (2)
            fun = 100.0_dp*(x(2)-x(1)**2)**2+(1.0_dp-x(1))**2
        case (3)
            fun = 20.0_dp+sum(x*x-10.0_dp*cos(2.0_dp*pi*x))
        case (4)
            a = sqrt(sum(x*x)/real(size(x),dp))
            b = sum(cos(2.0_dp*pi*x))/real(size(x),dp)
            fun = -20.0_dp*exp(-0.2_dp*a)-exp(b)+20.0_dp+exp(1.0_dp)
        case (5)
            fun = (x(1)-2.0_dp)**2+(x(2)+1.0_dp)**2
        case (6)
            fun = (x(1)-0.75_dp)**2+(x(2)+1.1_dp)**2+(x(3)-2.0_dp)**2
        end select
    end function objective
end program reference_differential_evolution_driver
