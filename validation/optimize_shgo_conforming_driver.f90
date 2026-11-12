! SPDX-License-Identifier: MIT
! Deterministic 2D conforming SHGO quality cases for SciPy comparison.
program optimize_shgo_conforming_driver
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result
    implicit none
    type(shgo_result) :: res
    real(dp) :: bounds(2,2)
    real(dp), parameter :: pi = acos(-1.0_dp)

    bounds(:,1) = -2.0_dp
    bounds(:,2) = 2.0_dp
    res = shgo_core(quadratic,bounds,iters=3,local_iter=0, &
        adaptive=.true.,conforming=.true.)
    call emit('quadratic',res)

    bounds(:,1) = -6.0_dp
    bounds(:,2) = 6.0_dp
    res = shgo_core(himmelblau,bounds,iters=5,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.35_dp)
    call emit('himmelblau',res)

    bounds(:,1) = -5.12_dp
    bounds(:,2) = 5.12_dp
    res = shgo_core(rastrigin,bounds,iters=4,local_iter=0, &
        adaptive=.true.,conforming=.true.)
    call emit('rastrigin',res)

    bounds(:,1) = -5.0_dp
    bounds(:,2) = 5.0_dp
    res = shgo_core(ackley,bounds,iters=3,local_iter=0, &
        adaptive=.true.,conforming=.true.)
    call emit('ackley',res)

    bounds(1,:) = [-3.0_dp,3.0_dp]
    bounds(2,:) = [-2.0_dp,2.0_dp]
    res = shgo_core(six_hump,bounds,iters=5,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.35_dp)
    call emit('six_hump',res)
contains
    subroutine emit(name,r)
        character(len=*), intent(in) :: name !! numerical benchmark name
        type(shgo_result), intent(in) :: r !! best sampled point and mesh
        integer :: valid
        valid = 0
        if (r%conforming_valid) valid = 1
        write(*,'(a,1x,5(i0,1x),es24.16,2(1x,es24.16))') &
            trim(name),r%status,r%nfev,r%nleaf,r%ntriangles,valid, &
            r%fun,r%x(1),r%x(2)
    end subroutine emit
    function quadratic(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional point
        real(dp) :: value !! convex shifted test objective
        value = sum((x-[0.625_dp,-0.425_dp])**2)
    end function quadratic
    function himmelblau(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional point
        real(dp) :: value !! Himmelblau objective
        value = (x(1)**2+x(2)-11.0_dp)**2 + (x(1)+x(2)**2-7.0_dp)**2
    end function himmelblau
    function rastrigin(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional point
        real(dp) :: value !! Rastrigin objective
        value = 20.0_dp+sum(x*x-10.0_dp*cos(2.0_dp*pi*x))
    end function rastrigin
    function ackley(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional point
        real(dp) :: value !! Ackley objective
        value = -20.0_dp*exp(-0.2_dp*sqrt(sum(x*x)/2.0_dp)) &
            -exp(sum(cos(2.0_dp*pi*x))/2.0_dp)+20.0_dp+exp(1.0_dp)
    end function ackley
    function six_hump(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional point
        real(dp) :: value !! six-hump camel function
        value = (4.0_dp-2.1_dp*x(1)**2+x(1)**4/3.0_dp)*x(1)**2 &
            +x(1)*x(2)+(-4.0_dp+4.0_dp*x(2)**2)*x(2)**2
    end function six_hump
end program optimize_shgo_conforming_driver
