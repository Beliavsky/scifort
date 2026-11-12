! SPDX-License-Identifier: MIT
! Repeatable adaptive SHGO samples for independent SciPy quality checks.
program optimize_shgo_adaptive_driver
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result
    implicit none
    type(shgo_result) :: res
    real(dp) :: bounds(2,2), bounds3(3,2), one(1,2)
    real(dp), parameter :: pi = acos(-1.0_dp)

    bounds(:,1) = -2.0_dp
    bounds(:,2) = 2.0_dp
    res = shgo_core(quadratic,bounds,iters=3,local_iter=0,adaptive=.true.)
    call emit('quadratic',res)

    bounds(:,1) = -6.0_dp
    bounds(:,2) = 6.0_dp
    res = shgo_core(himmelblau,bounds,iters=5,local_iter=0, &
        adaptive=.true.,adaptive_fraction=0.35_dp)
    call emit('himmelblau',res)

    bounds(:,1) = -5.12_dp
    bounds(:,2) = 5.12_dp
    res = shgo_core(rastrigin,bounds,iters=4,local_iter=0,adaptive=.true.)
    call emit('rastrigin',res)

    bounds(:,1) = -5.0_dp
    bounds(:,2) = 5.0_dp
    res = shgo_core(ackley,bounds,iters=3,local_iter=0,adaptive=.true.)
    call emit('ackley',res)

    one(1,:) = [-2.0_dp,2.0_dp]
    res = shgo_core(double_well,one,iters=3,local_iter=0,adaptive=.true.)
    call emit('double_well',res)

    bounds3(:,1) = -2.0_dp
    bounds3(:,2) = 2.0_dp
    res = shgo_core(quadratic3,bounds3,iters=4,local_iter=0,adaptive=.true.)
    call emit('quadratic3',res)

contains
    subroutine emit(name, res)
        character(len=*), intent(in) :: name !! objective identifier
        type(shgo_result), intent(in) :: res !! converged adaptive run
        integer :: j
        write(*,'(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,es24.16)',advance='no') &
            trim(name),res%nfev,res%nmesh,res%nleaf,res%status,res%fun
        do j = 1, size(res%x)
            write(*,'(1x,es24.16)',advance='no') res%x(j)
        end do
        write(*,*)
    end subroutine emit
    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! 2D position
        real(dp) :: f !! shifted quadratic
        f = sum((x-[0.625_dp,-0.425_dp])**2)
    end function quadratic
    function quadratic3(x) result(f)
        real(dp), intent(in) :: x(:) !! 3D position
        real(dp) :: f !! shifted quadratic
        f = sum((x-[0.625_dp,-0.425_dp,0.775_dp])**2)
    end function quadratic3
    function himmelblau(x) result(f)
        real(dp), intent(in) :: x(:) !! 2D position
        real(dp) :: f !! Himmelblau objective
        f = (x(1)**2+x(2)-11.0_dp)**2+(x(1)+x(2)**2-7.0_dp)**2
    end function himmelblau
    function rastrigin(x) result(f)
        real(dp), intent(in) :: x(:) !! 2D position
        real(dp) :: f !! Rastrigin objective
        f = 20.0_dp + sum(x**2-10.0_dp*cos(2.0_dp*pi*x))
    end function rastrigin
    function ackley(x) result(f)
        real(dp), intent(in) :: x(:) !! 2D position
        real(dp) :: f !! Ackley objective
        f = -20.0_dp*exp(-0.2_dp*sqrt(sum(x**2)/real(size(x),dp))) &
            - exp(sum(cos(2.0_dp*pi*x))/real(size(x),dp)) + 20.0_dp + exp(1.0_dp)
    end function ackley
    function double_well(x) result(f)
        real(dp), intent(in) :: x(:) !! 1D position
        real(dp) :: f !! tilted double well
        f = (x(1)**2-1.0_dp)**2+0.1_dp*x(1)
    end function double_well
end program optimize_shgo_adaptive_driver
