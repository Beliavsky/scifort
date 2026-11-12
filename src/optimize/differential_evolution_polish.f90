! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Optional bridge from SciFort differential evolution to its existing
! L-BFGS-B implementation.  Pass polisher=de_lbfgsb_polisher when calling
! differential_evolution.  Integer-valued coordinates are fixed at the
! incumbent by the caller; constrained solutions are not polished here.
!
module scifort_optimize_differential_evolution_polish
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_differential_evolution, only : de_objective
    use scifort_optimize_minimize, only : lbfgsb, multivariate_minimize_result
    implicit none
    private
    public :: de_lbfgsb_polisher
contains
    subroutine de_lbfgsb_polisher(func, x0, lower, upper, x, fun, nfev, success)
        procedure(de_objective) :: func !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameters, shape (n)
        real(dp), intent(in) :: lower(:) !! bounded minimizer lower bounds
        real(dp), intent(in) :: upper(:) !! bounded minimizer upper bounds
        real(dp), intent(out) :: x(:) !! polished parameters, shape (n)
        real(dp), intent(out) :: fun !! polished objective value
        integer, intent(out) :: nfev !! polishing function evaluations
        logical, intent(out) :: success !! successful bounded minimization
        type(multivariate_minimize_result) :: local

        success = .false.
        nfev = 0
        x = x0
        fun = huge(1.0_dp)
        if (size(x0) /= size(x) .or. size(lower) /= size(x0) .or. &
            size(upper) /= size(x0)) return
        if (any(lower > upper)) return
        local = lbfgsb(func, x0, lower=lower, upper=upper, maxiter=150, maxfun=5000)
        nfev = local%nfev
        if (.not. local%success .or. .not. allocated(local%x)) return
        if (size(local%x) /= size(x)) return
        if (.not. all(ieee_is_finite(local%x))) return
        if (.not. ieee_is_finite(local%fun)) return
        if (any(local%x < lower) .or. any(local%x > upper)) return
        x = local%x
        fun = local%fun
        success = .true.
    end subroutine de_lbfgsb_polisher
end module scifort_optimize_differential_evolution_polish
