! SPDX-License-Identifier: MIT
! Deterministic benchmark driver for Python/SciPy comparisons.
program optimize_basinhopping_driver
    use scifort_kinds, only : dp
    use scifort_optimize_basinhopping, only : basinhopping_core, basinhopping_result
    implicit none
    type(basinhopping_result) :: r
    integer :: case_id, seed
    real(dp) :: x1(1), x2(2), step, temp
    do case_id = 1, 6
        select case (case_id)
        case (1)
            x1 = [2.4_dp]
            step = 1.5_dp
            temp = 0.5_dp
        case (2)
            x2 = [1.8_dp, 1.7_dp]
            step = 2.8_dp
            temp = 2.0_dp
        case (3)
            x2 = [2.5_dp, 2.1_dp]
            step = 3.5_dp
            temp = 2.0_dp
        case (4)
            x2 = [2.0_dp, 1.8_dp]
            step = 2.3_dp
            temp = 1.5_dp
        case (5)
            x1 = [2.0_dp]
            step = 1.3_dp
            temp = 0.5_dp
        case (6)
            x2 = [2.0_dp, 1.8_dp]
            step = 1.8_dp
            temp = 0.5_dp
        end select
        do seed = 1, 8
            if (case_id == 1 .or. case_id == 5) then
                r = basinhopping_core(benchmark, x1, coordinate_minimize, niter=120, &
                    stepsize=step, temperature=temp, seed=119*seed + 17)
            else
                r = basinhopping_core(benchmark, x2, coordinate_minimize, niter=120, &
                    stepsize=step, temperature=temp, seed=119*seed + 17)
            end if
            if (.not. r%success) error stop 'basinhopping returned invalid result'
            if (abs(benchmark(r%x) - r%fun) > 1.0e-7_dp) error stop 'inconsistent objective'
            if (r%naccepted + r%nrejected /= 120) error stop 'inconsistent acceptance count'
            write(*,'(2(i0,1x),es24.16,1x,i0,1x,i0)') case_id, seed, r%fun, r%naccepted, r%nfev
        end do
    end do
contains
    function benchmark(x) result(fx)
        real(dp), intent(in) :: x(:)
        real(dp) :: fx, pi, sx, sy
        pi = acos(-1.0_dp)
        select case (case_id)
        case (1)
            fx = (x(1)**2 - 1.0_dp)**2 + 0.25_dp*x(1)
        case (2)
            fx = 20.0_dp + sum(x*x - 10.0_dp*cos(2.0_dp*pi*x))
        case (3)
            fx = (x(1)**2 + x(2) - 11.0_dp)**2 + (x(1) + x(2)**2 - 7.0_dp)**2
        case (4)
            sx = sqrt(sum(x*x)/real(size(x), dp))
            sy = sum(cos(2.0_dp*pi*x))/real(size(x), dp)
            fx = -20.0_dp*exp(-0.2_dp*sx) - exp(sy) + 20.0_dp + exp(1.0_dp)
        case (5)
            fx = sin(4.0_dp*x(1))**2 + 0.03_dp*x(1)**2
        case (6)
            fx = (x(1)**2 - 1.0_dp)**2 + 0.20_dp*x(1) + &
                0.5_dp*(x(2)**2 - 1.0_dp)**2 + 0.06_dp*x(2)
        end select
    end function benchmark

    subroutine coordinate_minimize(fun, x0, x, value, nfev, success)
        interface
            function fun(z) result(fz)
                import dp
                real(dp), intent(in) :: z(:)
                real(dp) :: fz
            end function fun
        end interface
        real(dp), intent(in) :: x0(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: value
        integer, intent(out) :: nfev
        logical, intent(out) :: success
        real(dp) :: delta, fz, z(size(x))
        integer :: i, j, direction
        logical :: improved
        x = x0
        value = fun(x)
        nfev = 1
        delta = 0.5_dp
        do i = 1, 500
            improved = .false.
            do j = 1, size(x)
                do direction = -1, 1, 2
                    z = x
                    z(j) = x(j) + real(direction, dp)*delta
                    fz = fun(z)
                    nfev = nfev + 1
                    if (fz < value) then
                        value = fz
                        x = z
                        improved = .true.
                    end if
                end do
            end do
            if (.not.improved) delta = delta/2.0_dp
            if (delta < 1.0e-8_dp) exit
        end do
        success = .true.
    end subroutine coordinate_minimize
end program optimize_basinhopping_driver
