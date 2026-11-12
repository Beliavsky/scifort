! SPDX-License-Identifier: MIT
! Independent numerical-quality checks for conforming adaptive 3D SHGO.
program optimize_shgo_conforming_3d_driver
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result
    implicit none
    type(shgo_result) :: res
    real(dp) :: bounds(3,2)
    real(dp), parameter :: pi = acos(-1.0_dp)

    bounds(:,1) = -1.0_dp
    bounds(:,2) = 1.0_dp
    res = shgo_core(quadratic,bounds,iters=4,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.25_dp)
    call emit('quadratic',res)

    bounds(:,1) = -5.12_dp
    bounds(:,2) = 5.12_dp
    res = shgo_core(rastrigin,bounds,iters=3,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.25_dp)
    call emit('rastrigin',res)

    bounds(:,1) = -5.0_dp
    bounds(:,2) = 5.0_dp
    res = shgo_core(ackley,bounds,iters=3,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.25_dp)
    call emit('ackley',res)

    bounds(:,1) = -2.0_dp
    bounds(:,2) = 2.0_dp
    res = shgo_core(double_well,bounds,iters=4,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.25_dp)
    call emit('double_well',res)

    bounds(:,1) = -3.0_dp
    bounds(:,2) = 3.0_dp
    res = shgo_core(rosenbrock,bounds,iters=4,local_iter=0, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.25_dp)
    call emit('rosenbrock',res)
contains
    subroutine emit(name,r)
        character(len=*), intent(in) :: name !! benchmark name
        type(shgo_result), intent(in) :: r !! sampled best and validated topology
        integer :: validated
        validated = 0
        if (r%conforming_valid) validated = 1
        write(*,'(a,1x,5(i0,1x),es24.16,3(1x,es24.16))') &
            trim(name),r%status,r%nfev,r%nleaf,r%ntetrahedra,validated, &
            r%fun,r%x(1),r%x(2),r%x(3)
    end subroutine emit
    function quadratic(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional point
        real(dp) :: value !! shifted convex quadratic
        value = sum((x-[0.3_dp,-0.2_dp,0.15_dp])**2)
    end function quadratic
    function rastrigin(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional point
        real(dp) :: value !! multimodal Rastrigin function
        value = 30.0_dp+sum(x*x-10.0_dp*cos(2.0_dp*pi*x))
    end function rastrigin
    function ackley(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional point
        real(dp) :: value !! multimodal Ackley function
        value = -20.0_dp*exp(-0.2_dp*sqrt(sum(x*x)/3.0_dp)) - &
            exp(sum(cos(2.0_dp*pi*x))/3.0_dp)+20.0_dp+exp(1.0_dp)
    end function ackley
    function double_well(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional point
        real(dp) :: value !! separated minima along first axis
        value = (x(1)*x(1)-1.0_dp)**2 + (x(2)-0.35_dp)**2 + &
            (x(3)+0.55_dp)**2
    end function double_well
    function rosenbrock(x) result(value)
        real(dp), intent(in) :: x(:) !! three-dimensional point
        real(dp) :: value !! curved Rosenbrock valley
        value = 100.0_dp*(x(2)-x(1)**2)**2+(1.0_dp-x(1))**2 + &
            100.0_dp*(x(3)-x(2)**2)**2+(1.0_dp-x(2))**2
    end function rosenbrock
end program optimize_shgo_conforming_3d_driver
