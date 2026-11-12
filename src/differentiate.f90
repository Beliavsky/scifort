! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Copyright (c) 2023, Seyed Ali Ghasemi
!
! Adaptive numerical differentiation with SciPy-compatible public semantics.
! The nested finite-difference stencil, weight construction, error estimate,
! and termination rules are adapted from SciPy 1.17.0
! scipy/differentiate/_differentiate.py (BSD-3-Clause).
!
! The scalar complex-step formula is adapted from ForDiff 0.1.1 by Seyed Ali Ghasemi
! (BSD-3-Clause). The vector complex-step
! routines are SciFort implementations of the same standard method.
!
! Jacobian API organization was cross-checked against NumDiff 1.7.0 by Jacob
! Williams (BSD-style license); no NumDiff implementation code is copied here.

module scifort_differentiate
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: differentiate_status_success = 0
    integer, parameter, public :: differentiate_status_error_increase = -1
    integer, parameter, public :: differentiate_status_maxiter = -2
    integer, parameter, public :: differentiate_status_nonfinite = -3
    integer, parameter, public :: differentiate_status_invalid_input = -5

    type, public :: derivative_result
        logical :: success = .false. !! true when the adaptive iteration converged
        integer :: status = differentiate_status_invalid_input !! numerical termination code
        real(dp) :: df = 0.0_dp !! derivative estimate at x
        real(dp) :: error = 0.0_dp !! absolute difference from the previous estimate
        integer :: nit = 0 !! number of adaptive finite-difference iterations performed
        integer :: nfev = 0 !! number of scalar function evaluations used
        real(dp) :: x = 0.0_dp !! abscissa at which the derivative was evaluated
    end type derivative_result

    type, public :: derivative_array_result
        logical, allocatable :: success(:) !! convergence flag for each abscissa
        integer, allocatable :: status(:) !! termination code for each abscissa
        real(dp), allocatable :: df(:) !! derivative estimate for each abscissa
        real(dp), allocatable :: error(:) !! estimated error for each abscissa
        integer, allocatable :: nit(:) !! adaptive iterations for each abscissa
        integer, allocatable :: nfev(:) !! scalar function evaluations for each abscissa
        real(dp), allocatable :: x(:) !! abscissae at which derivatives were evaluated
    end type derivative_array_result

    type, public :: jacobian_result
        logical, allocatable :: success(:, :) !! convergence flag for each Jacobian element
        integer, allocatable :: status(:, :) !! termination code for each Jacobian element
        real(dp), allocatable :: df(:, :) !! Jacobian, shape (number of outputs, number of inputs)
        real(dp), allocatable :: error(:, :) !! estimated error for each Jacobian element
        integer, allocatable :: nit(:, :) !! adaptive iterations for each Jacobian element
        integer, allocatable :: nfev(:, :) !! function evaluations associated with each element
    end type jacobian_result

    type, public :: hessian_result
        logical, allocatable :: success(:, :) !! convergence flag for each Hessian element
        integer, allocatable :: status(:, :) !! termination code for each Hessian element
        real(dp), allocatable :: ddf(:, :) !! Hessian matrix
        real(dp), allocatable :: error(:, :) !! estimated outer-differentiation error
        integer, allocatable :: nfev(:, :) !! outer gradient evaluations for each Hessian element
    end type hessian_result

    abstract interface
        function scalar_real_function(x) result(fx)
            import dp
            implicit none
            real(dp), intent(in) :: x !! scalar abscissa
            real(dp) :: fx !! scalar function value
        end function scalar_real_function

        function scalar_vector_function(x) result(fx)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! vector argument
            real(dp) :: fx !! scalar function value
        end function scalar_vector_function

        function vector_real_function(x) result(fx)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! vector argument
            real(dp), allocatable :: fx(:) !! vector function value
        end function vector_real_function

        function scalar_complex_function(z) result(fz)
            import dp
            implicit none
            complex(dp), intent(in) :: z !! complex scalar argument
            complex(dp) :: fz !! complex scalar function value
        end function scalar_complex_function

        function scalar_complex_vector_function(z) result(fz)
            import dp
            implicit none
            complex(dp), intent(in) :: z(:) !! complex vector argument
            complex(dp) :: fz !! complex scalar function value
        end function scalar_complex_vector_function

        function vector_complex_function(z) result(fz)
            import dp
            implicit none
            complex(dp), intent(in) :: z(:) !! complex vector argument
            complex(dp), allocatable :: fz(:) !! complex vector function value
        end function vector_complex_function
    end interface

    public :: scalar_real_function
    public :: scalar_vector_function
    public :: vector_real_function
    public :: scalar_complex_function
    public :: scalar_complex_vector_function
    public :: vector_complex_function
    public :: derivative
    public :: jacobian
    public :: hessian
    public :: complex_step_derivative
    public :: complex_step_gradient
    public :: complex_step_jacobian

    interface derivative
        module procedure derivative_scalar
        module procedure derivative_array
    end interface derivative

    interface jacobian
        module procedure jacobian_scalar_output
        module procedure jacobian_vector_output
    end interface jacobian

contains

    recursive function derivative_scalar(f, x, atol, rtol, maxiter, order, initial_step, &
            step_factor, step_direction) result(res)
        procedure(scalar_real_function) :: f !! real scalar function to differentiate
        real(dp), intent(in) :: x !! finite point at which the derivative is desired
        real(dp), intent(in), optional :: atol !! absolute convergence tolerance
        real(dp), intent(in), optional :: rtol !! relative convergence tolerance
        integer, intent(in), optional :: maxiter !! maximum adaptive iterations, default 10
        integer, intent(in), optional :: order !! finite-difference order, default 8
        real(dp), intent(in), optional :: initial_step !! positive initial maximum step, default 0.5
        real(dp), intent(in), optional :: step_factor !! step reduction factor, default 2
        integer, intent(in), optional :: step_direction !! negative, zero, or positive stencil direction
        type(derivative_result) :: res

        integer :: direction
        integer :: effective_order
        integer :: max_iterations
        integer :: n
        real(dp) :: abs_tol
        real(dp) :: current_error
        real(dp) :: df_current
        real(dp) :: df_last
        real(dp) :: error_last
        real(dp) :: fac
        real(dp) :: fx
        real(dp) :: h
        real(dp) :: rel_tol
        real(dp), allocatable :: weights_central(:)
        real(dp), allocatable :: weights_right(:)
        logical :: weights_ok

        call initialize_derivative_result(res, x)

        max_iterations = 10
        if (present(maxiter)) max_iterations = maxiter
        effective_order = 8
        if (present(order)) effective_order = order
        h = 0.5_dp
        if (present(initial_step)) h = initial_step
        fac = 2.0_dp
        if (present(step_factor)) fac = step_factor
        direction = 0
        if (present(step_direction)) direction = sign(1, step_direction)
        if (present(step_direction)) then
            if (step_direction == 0) direction = 0
        end if
        abs_tol = tiny(1.0_dp)
        if (present(atol)) abs_tol = atol
        rel_tol = sqrt(epsilon(1.0_dp))
        if (present(rtol)) rel_tol = rtol

        if (max_iterations <= 0 .or. effective_order <= 0) return
        if (.not. ieee_is_finite(abs_tol) .or. abs_tol < 0.0_dp) return
        if (.not. ieee_is_finite(rel_tol) .or. rel_tol < 0.0_dp) return
        if (.not. ieee_is_finite(fac) .or. fac <= 0.0_dp) return
        if (.not. ieee_is_finite(x)) then
            res%status = differentiate_status_nonfinite
            res%df = quiet_nan()
            res%error = quiet_nan()
            return
        end if

        fx = f(x)
        res%nfev = 1
        if (.not. ieee_is_finite(fx)) then
            res%status = differentiate_status_nonfinite
            res%df = quiet_nan()
            res%error = quiet_nan()
            return
        end if
        if (.not. ieee_is_finite(h) .or. h <= 0.0_dp) then
            res%status = differentiate_status_nonfinite
            res%df = quiet_nan()
            res%error = quiet_nan()
            return
        end if
        if (abs(fac - 1.0_dp) <= epsilon(1.0_dp)) return

        n = (effective_order + 1) / 2
        call derivative_weights(n, fac, weights_central, weights_right, weights_ok)
        if (.not. weights_ok) return

        df_last = quiet_nan()
        error_last = quiet_nan()
        current_error = quiet_nan()

        if (direction == 0) then
            call adaptive_central(f, x, fx, n, h, fac, weights_central, &
                abs_tol, rel_tol, max_iterations, res, df_last, error_last)
        else
            call adaptive_one_sided(f, x, fx, n, h, fac, direction, weights_right, &
                abs_tol, rel_tol, max_iterations, res, df_last, error_last)
        end if

        if (res%status == differentiate_status_success) res%success = .true.
        df_current = res%df
        current_error = res%error
        if (.not. ieee_is_finite(df_current)) then
            res%success = .false.
            res%status = differentiate_status_nonfinite
            res%df = quiet_nan()
            res%error = quiet_nan()
        else if (.not. ieee_is_finite(current_error) .and. res%nit > 1) then
            res%success = .false.
            res%status = differentiate_status_nonfinite
            res%df = quiet_nan()
            res%error = quiet_nan()
        end if
    end function derivative_scalar

    function derivative_array(f, x, atol, rtol, maxiter, order, initial_step, &
            step_factor, step_direction) result(res)
        procedure(scalar_real_function) :: f !! real scalar function to differentiate elementwise
        real(dp), intent(in) :: x(:) !! points at which derivatives are desired
        real(dp), intent(in), optional :: atol !! absolute convergence tolerance
        real(dp), intent(in), optional :: rtol !! relative convergence tolerance
        integer, intent(in), optional :: maxiter !! maximum adaptive iterations, default 10
        integer, intent(in), optional :: order !! finite-difference order, default 8
        real(dp), intent(in), optional :: initial_step !! positive initial maximum step, default 0.5
        real(dp), intent(in), optional :: step_factor !! step reduction factor, default 2
        integer, intent(in), optional :: step_direction !! common stencil direction for all points
        type(derivative_array_result) :: res

        integer :: i
        type(derivative_result) :: one

        allocate(res%success(size(x)), res%status(size(x)), res%df(size(x)), &
            res%error(size(x)), res%nit(size(x)), res%nfev(size(x)), res%x(size(x)))

        do i = 1, size(x)
            one = derivative_scalar(f, x(i), atol=atol, rtol=rtol, maxiter=maxiter, &
                order=order, initial_step=initial_step, step_factor=step_factor, &
                step_direction=step_direction)
            res%success(i) = one%success
            res%status(i) = one%status
            res%df(i) = one%df
            res%error(i) = one%error
            res%nit(i) = one%nit
            res%nfev(i) = one%nfev
            res%x(i) = one%x
        end do
    end function derivative_array

    function jacobian_scalar_output(f, x, atol, rtol, maxiter, order, initial_step, &
            step_factor, step_direction) result(res)
        procedure(scalar_vector_function) :: f !! scalar-valued function of a vector argument
        real(dp), intent(in) :: x(:) !! point at which the gradient is desired
        real(dp), intent(in), optional :: atol !! absolute convergence tolerance
        real(dp), intent(in), optional :: rtol !! relative convergence tolerance
        integer, intent(in), optional :: maxiter !! maximum adaptive iterations, default 10
        integer, intent(in), optional :: order !! finite-difference order, default 8
        real(dp), intent(in), optional :: initial_step !! positive initial maximum step, default 0.5
        real(dp), intent(in), optional :: step_factor !! step reduction factor, default 2
        integer, intent(in), optional :: step_direction(:) !! stencil direction for each input variable
        type(jacobian_result) :: res

        integer :: current_input
        integer :: direction
        integer :: j
        type(derivative_result) :: one

        call allocate_jacobian_result(res, 1, size(x))
        if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) then
            call set_jacobian_nonfinite(res)
            return
        end if
        if (present(step_direction)) then
            if (size(step_direction) /= size(x)) return
        end if

        do j = 1, size(x)
            current_input = j
            direction = 0
            if (present(step_direction)) direction = step_direction(j)
            one = derivative_scalar(component_function, x(j), atol=atol, rtol=rtol, &
                maxiter=maxiter, order=order, initial_step=initial_step, &
                step_factor=step_factor, step_direction=direction)
            call store_jacobian_element(res, 1, j, one)
        end do

    contains

        function component_function(t) result(value)
            real(dp), intent(in) :: t !! perturbed coordinate value
            real(dp) :: value
            real(dp) :: z(size(x))

            z = x
            z(current_input) = t
            value = f(z)
        end function component_function

    end function jacobian_scalar_output

    function jacobian_vector_output(f, x, atol, rtol, maxiter, order, initial_step, &
            step_factor, step_direction) result(res)
        procedure(vector_real_function) :: f !! vector-valued function of a vector argument
        real(dp), intent(in) :: x(:) !! point at which the Jacobian is desired
        real(dp), intent(in), optional :: atol !! absolute convergence tolerance
        real(dp), intent(in), optional :: rtol !! relative convergence tolerance
        integer, intent(in), optional :: maxiter !! maximum adaptive iterations, default 10
        integer, intent(in), optional :: order !! finite-difference order, default 8
        real(dp), intent(in), optional :: initial_step !! positive initial maximum step, default 0.5
        real(dp), intent(in), optional :: step_factor !! step reduction factor, default 2
        integer, intent(in), optional :: step_direction(:) !! stencil direction for each input variable
        type(jacobian_result) :: res

        integer :: current_input
        integer :: current_output
        integer :: direction
        integer :: i
        integer :: j
        real(dp), allocatable :: y0(:)
        type(derivative_result) :: one

        if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) then
            call allocate_jacobian_result(res, 1, max(1, size(x)))
            call set_jacobian_nonfinite(res)
            return
        end if
        if (present(step_direction)) then
            if (size(step_direction) /= size(x)) then
                call allocate_jacobian_result(res, 1, size(x))
                return
            end if
        end if

        y0 = f(x)
        if (.not. allocated(y0) .or. size(y0) == 0) then
            call allocate_jacobian_result(res, 1, size(x))
            return
        end if
        call allocate_jacobian_result(res, size(y0), size(x))
        if (any(.not. ieee_is_finite(y0))) then
            call set_jacobian_nonfinite(res)
            return
        end if

        do j = 1, size(x)
            current_input = j
            direction = 0
            if (present(step_direction)) direction = step_direction(j)
            do i = 1, size(y0)
                current_output = i
                one = derivative_scalar(component_function, x(j), atol=atol, rtol=rtol, &
                    maxiter=maxiter, order=order, initial_step=initial_step, &
                    step_factor=step_factor, step_direction=direction)
                call store_jacobian_element(res, i, j, one)
            end do
        end do

    contains

        function component_function(t) result(value)
            real(dp), intent(in) :: t !! perturbed coordinate value
            real(dp) :: value
            real(dp) :: z(size(x))
            real(dp), allocatable :: yz(:)

            z = x
            z(current_input) = t
            yz = f(z)
            if (.not. allocated(yz)) then
                value = quiet_nan()
            else if (size(yz) < current_output) then
                value = quiet_nan()
            else
                value = yz(current_output)
            end if
        end function component_function

    end function jacobian_vector_output

    function hessian(f, x, atol, rtol, maxiter, order, initial_step, step_factor) result(res)
        procedure(scalar_vector_function) :: f !! scalar-valued function whose Hessian is desired
        real(dp), intent(in) :: x(:) !! point at which the Hessian is desired
        real(dp), intent(in), optional :: atol !! absolute convergence tolerance
        real(dp), intent(in), optional :: rtol !! relative convergence tolerance
        integer, intent(in), optional :: maxiter !! maximum adaptive iterations, default 10
        integer, intent(in), optional :: order !! finite-difference order, default 8
        real(dp), intent(in), optional :: initial_step !! positive initial maximum step, default 0.5
        real(dp), intent(in), optional :: step_factor !! step reduction factor, default 2
        type(hessian_result) :: res

        real(dp) :: abs_tol
        integer :: current_input
        integer :: current_output
        integer :: evaluation_count
        integer :: i
        integer :: j
        real(dp) :: inner_rtol
        real(dp) :: rel_tol
        real(dp), allocatable :: current_z(:)
        type(derivative_result) :: one

        call allocate_hessian_result(res, size(x))
        if (size(x) == 0) return
        if (any(.not. ieee_is_finite(x))) then
            res%status = differentiate_status_nonfinite
            return
        end if

        abs_tol = tiny(1.0_dp)
        if (present(atol)) abs_tol = atol
        rel_tol = sqrt(epsilon(1.0_dp))
        if (present(rtol)) rel_tol = rtol
        if (rel_tol > 0.0_dp .and. rel_tol < 100.0_dp * epsilon(1.0_dp)) then
            rel_tol = 100.0_dp * epsilon(1.0_dp)
        end if
        inner_rtol = rel_tol / 100.0_dp
        allocate(current_z(size(x)))

        do j = 1, size(x)
            current_input = j
            do i = 1, size(x)
                current_output = i
                evaluation_count = 0
                one = derivative_scalar(outer_component, x(j), atol=abs_tol, rtol=rel_tol, &
                    maxiter=maxiter, order=order, initial_step=initial_step, &
                    step_factor=step_factor)
                res%success(i, j) = one%success
                res%status(i, j) = one%status
                res%ddf(i, j) = one%df
                res%error(i, j) = one%error
                res%nfev(i, j) = evaluation_count
            end do
        end do

    contains

        function outer_component(t) result(value)
            real(dp), intent(in) :: t !! perturbation of the outer differentiation coordinate
            real(dp) :: value
            type(derivative_result) :: inner

            current_z = x
            current_z(current_input) = t
            inner = derivative_scalar(inner_component, current_z(current_output), &
                atol=abs_tol, rtol=inner_rtol, maxiter=maxiter, order=order, &
                initial_step=initial_step, step_factor=step_factor)
            evaluation_count = evaluation_count + inner%nfev
            value = inner%df
        end function outer_component

        function inner_component(t) result(value)
            real(dp), intent(in) :: t !! perturbation of the inner differentiation coordinate
            real(dp) :: value
            real(dp) :: z(size(x))

            z = current_z
            z(current_output) = t
            value = f(z)
        end function inner_component

    end function hessian

    function complex_step_derivative(f, x, h, status) result(dfdx)
        procedure(scalar_complex_function) :: f !! complex extension of a scalar analytic function
        real(dp), intent(in) :: x !! real point at which the derivative is desired
        real(dp), intent(in), optional :: h !! positive imaginary perturbation, default 1e-30
        integer, intent(out), optional :: status !! zero on success; invalid-input code otherwise
        real(dp) :: dfdx

        real(dp) :: step
        complex(dp) :: value

        step = 1.0e-30_dp
        if (present(h)) step = h
        if (present(status)) status = differentiate_status_invalid_input
        dfdx = quiet_nan()
        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(step) .or. step <= 0.0_dp) return

        value = f(cmplx(x, step, kind=dp))
        if (.not. ieee_is_finite(real(value, dp)) .or. .not. ieee_is_finite(aimag(value))) return
        dfdx = aimag(value) / step
        if (present(status)) status = differentiate_status_success
    end function complex_step_derivative

    function complex_step_gradient(f, x, h, status) result(grad)
        procedure(scalar_complex_vector_function) :: f !! complex extension of a scalar vector function
        real(dp), intent(in) :: x(:) !! real vector point at which the gradient is desired
        real(dp), intent(in), optional :: h !! positive imaginary perturbation, default 1e-30
        integer, intent(out), optional :: status !! zero on success; invalid-input code otherwise
        real(dp), allocatable :: grad(:)

        integer :: i
        real(dp) :: step
        complex(dp), allocatable :: z(:)
        complex(dp) :: value

        allocate(grad(size(x)))
        grad = quiet_nan()
        step = 1.0e-30_dp
        if (present(h)) step = h
        if (present(status)) status = differentiate_status_invalid_input
        if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) return
        if (.not. ieee_is_finite(step) .or. step <= 0.0_dp) return

        allocate(z(size(x)))
        do i = 1, size(x)
            z = cmplx(x, 0.0_dp, kind=dp)
            z(i) = cmplx(x(i), step, kind=dp)
            value = f(z)
            if (.not. ieee_is_finite(real(value, dp)) .or. .not. ieee_is_finite(aimag(value))) return
            grad(i) = aimag(value) / step
        end do
        if (present(status)) status = differentiate_status_success
    end function complex_step_gradient

    function complex_step_jacobian(f, x, h, status) result(jac)
        procedure(vector_complex_function) :: f !! complex extension of a vector-valued function
        real(dp), intent(in) :: x(:) !! real vector point at which the Jacobian is desired
        real(dp), intent(in), optional :: h !! positive imaginary perturbation, default 1e-30
        integer, intent(out), optional :: status !! zero on success; invalid-input code otherwise
        real(dp), allocatable :: jac(:, :)

        integer :: i
        integer :: m
        real(dp) :: step
        complex(dp), allocatable :: value(:)
        complex(dp), allocatable :: z(:)

        step = 1.0e-30_dp
        if (present(h)) step = h
        if (present(status)) status = differentiate_status_invalid_input
        if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) then
            allocate(jac(0, size(x)))
            return
        end if
        if (.not. ieee_is_finite(step) .or. step <= 0.0_dp) then
            allocate(jac(0, size(x)))
            return
        end if

        z = cmplx(x, 0.0_dp, kind=dp)
        value = f(z)
        if (.not. allocated(value)) then
            allocate(jac(0, size(x)))
            return
        end if
        m = size(value)
        allocate(jac(m, size(x)))
        jac = quiet_nan()
        do i = 1, size(x)
            z = cmplx(x, 0.0_dp, kind=dp)
            z(i) = cmplx(x(i), step, kind=dp)
            value = f(z)
            if (.not. allocated(value) .or. size(value) /= m) return
            if (any(.not. ieee_is_finite(real(value, dp))) .or. &
                    any(.not. ieee_is_finite(aimag(value)))) return
            jac(:, i) = aimag(value) / step
        end do
        if (present(status)) status = differentiate_status_success
    end function complex_step_jacobian

    recursive subroutine adaptive_central(f, x, fx, n, h0, fac, weights, atol, rtol, &
            maxiter, res, df_last, error_last)
        procedure(scalar_real_function) :: f !! scalar function being differentiated
        real(dp), intent(in) :: x !! differentiation point
        real(dp), intent(in) :: fx !! function value at x
        integer, intent(in) :: n !! half the effective even finite-difference order
        real(dp), intent(in) :: h0 !! initial maximum step
        real(dp), intent(in) :: fac !! step reduction factor
        real(dp), intent(in) :: weights(:) !! central finite-difference weights
        real(dp), intent(in) :: atol !! absolute convergence tolerance
        real(dp), intent(in) :: rtol !! relative convergence tolerance
        integer, intent(in) :: maxiter !! maximum adaptive iterations
        type(derivative_result), intent(inout) :: res !! result accumulated by the adaptive loop
        real(dp), intent(inout) :: df_last !! previous derivative estimate
        real(dp), intent(inout) :: error_last !! previous error estimate

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: fnew(:)
        real(dp), allocatable :: fstore(:)
        real(dp) :: df
        real(dp) :: err
        real(dp) :: h
        integer :: i
        integer :: iter
        integer :: mid

        h = h0
        allocate(fstore(1))
        fstore(1) = fx
        df = quiet_nan()
        err = quiet_nan()

        do iter = 1, maxiter
            if (iter == 1) then
                allocate(fnew(2 * n))
                do i = 1, n
                    fnew(i) = f(x - h / fac ** real(n - i, dp))
                    fnew(n + i) = f(x + h / fac ** real(i - 1, dp))
                end do
                res%nfev = res%nfev + 2 * n
                fstore = [fnew(1:n), fx, fnew(n + 1:2 * n)]
                deallocate(fnew)
            else
                allocate(fnew(2))
                fnew(1) = f(x - h / fac ** real(n - 1, dp))
                fnew(2) = f(x + h / fac ** real(n - 1, dp))
                res%nfev = res%nfev + 2
                fstore = [fnew(1), fstore, fnew(2)]
                deallocate(fnew)
            end if

            res%nit = iter
            if (any(.not. ieee_is_finite(fstore))) then
                res%status = differentiate_status_nonfinite
                res%df = quiet_nan()
                res%error = quiet_nan()
                return
            end if

            mid = (size(fstore) + 1) / 2
            allocate(current(2 * n + 1))
            current(1:n) = fstore(1:n)
            current(n + 1) = fstore(mid)
            current(n + 2:2 * n + 1) = fstore(size(fstore) - n + 1:size(fstore))
            df = dot_product(current, weights) / h
            deallocate(current)

            if (.not. ieee_is_finite(df)) then
                res%status = differentiate_status_nonfinite
                res%df = quiet_nan()
                res%error = quiet_nan()
                return
            end if

            if (iter == 1) then
                err = quiet_nan()
            else
                err = abs(df - df_last)
            end if
            res%df = df
            res%error = err

            if (iter > 1) then
                if (err < atol + rtol * abs(df)) then
                    res%status = differentiate_status_success
                    return
                end if
                if (iter > 2 .and. ieee_is_finite(error_last)) then
                    if (err > 10.0_dp * error_last) then
                        res%status = differentiate_status_error_increase
                        return
                    end if
                end if
            end if

            df_last = df
            error_last = err
            h = h / fac
            if (.not. ieee_is_finite(h) .or. h <= 0.0_dp) then
                res%status = differentiate_status_nonfinite
                res%df = quiet_nan()
                res%error = quiet_nan()
                return
            end if
        end do
        res%status = differentiate_status_maxiter
    end subroutine adaptive_central

    recursive subroutine adaptive_one_sided(f, x, fx, n, h0, fac, direction, weights, atol, rtol, &
            maxiter, res, df_last, error_last)
        procedure(scalar_real_function) :: f !! scalar function being differentiated
        real(dp), intent(in) :: x !! differentiation point
        real(dp), intent(in) :: fx !! function value at x
        integer, intent(in) :: n !! half the effective even finite-difference order
        real(dp), intent(in) :: h0 !! initial maximum step
        real(dp), intent(in) :: fac !! step reduction factor
        integer, intent(in) :: direction !! -1 for left or +1 for right one-sided differences
        real(dp), intent(in) :: weights(:) !! right-sided finite-difference weights
        real(dp), intent(in) :: atol !! absolute convergence tolerance
        real(dp), intent(in) :: rtol !! relative convergence tolerance
        integer, intent(in) :: maxiter !! maximum adaptive iterations
        type(derivative_result), intent(inout) :: res !! result accumulated by the adaptive loop
        real(dp), intent(inout) :: df_last !! previous derivative estimate
        real(dp), intent(inout) :: error_last !! previous error estimate

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: fnew(:)
        real(dp), allocatable :: fstore(:)
        real(dp) :: d
        real(dp) :: df
        real(dp) :: err
        real(dp) :: h
        real(dp) :: step
        integer :: i
        integer :: iter

        d = sqrt(fac)
        h = h0
        allocate(fstore(1))
        fstore(1) = fx
        df = quiet_nan()
        err = quiet_nan()

        do iter = 1, maxiter
            if (iter == 1) then
                allocate(fnew(2 * n))
                do i = 1, 2 * n
                    step = h / d ** real(i - 1, dp)
                    fnew(i) = f(x + real(direction, dp) * step)
                end do
                res%nfev = res%nfev + 2 * n
                fstore = [fstore, fnew]
                deallocate(fnew)
            else
                allocate(fnew(2))
                fnew(1) = f(x + real(direction, dp) * h / fac ** real(n - 1, dp))
                fnew(2) = f(x + real(direction, dp) * h / &
                    (fac ** real(n - 1, dp) * d))
                res%nfev = res%nfev + 2
                fstore = [fstore, fnew]
                deallocate(fnew)
            end if

            res%nit = iter
            if (any(.not. ieee_is_finite(fstore))) then
                res%status = differentiate_status_nonfinite
                res%df = quiet_nan()
                res%error = quiet_nan()
                return
            end if

            allocate(current(2 * n + 1))
            current(1) = fstore(1)
            current(2:2 * n + 1) = fstore(size(fstore) - 2 * n + 1:size(fstore))
            df = dot_product(current, weights) / h
            if (direction < 0) df = -df
            deallocate(current)

            if (.not. ieee_is_finite(df)) then
                res%status = differentiate_status_nonfinite
                res%df = quiet_nan()
                res%error = quiet_nan()
                return
            end if

            if (iter == 1) then
                err = quiet_nan()
            else
                err = abs(df - df_last)
            end if
            res%df = df
            res%error = err

            if (iter > 1) then
                if (err < atol + rtol * abs(df)) then
                    res%status = differentiate_status_success
                    return
                end if
                if (iter > 2 .and. ieee_is_finite(error_last)) then
                    if (err > 10.0_dp * error_last) then
                        res%status = differentiate_status_error_increase
                        return
                    end if
                end if
            end if

            df_last = df
            error_last = err
            h = h / fac
            if (.not. ieee_is_finite(h) .or. h <= 0.0_dp) then
                res%status = differentiate_status_nonfinite
                res%df = quiet_nan()
                res%error = quiet_nan()
                return
            end if
        end do
        res%status = differentiate_status_maxiter
    end subroutine adaptive_one_sided

    subroutine derivative_weights(n, fac, central, right, ok)
        integer, intent(in) :: n !! half the effective even finite-difference order
        real(dp), intent(in) :: fac !! step reduction factor
        real(dp), allocatable, intent(out) :: central(:) !! central-difference weights
        real(dp), allocatable, intent(out) :: right(:) !! right one-sided weights
        logical, intent(out) :: ok !! true when both Vandermonde systems were solved

        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: nodes(:)
        integer :: i
        integer :: j
        integer :: k
        integer :: m
        real(dp) :: d

        ok = .false.
        m = 2 * n + 1
        allocate(central(m), right(m), a(m, m), b(m), nodes(m))
        central = 0.0_dp
        right = 0.0_dp
        b = 0.0_dp
        b(2) = 1.0_dp

        do i = 1, m
            k = i - (n + 1)
            if (k == 0) then
                nodes(i) = 0.0_dp
            else
                nodes(i) = real(sign(1, k), dp) / fac ** real(abs(k) - 1, dp)
            end if
        end do
        do j = 1, m
            do i = 1, m
                a(j, i) = nodes(i) ** real(j - 1, dp)
            end do
        end do
        call solve_dense(a, b, central, ok)
        if (.not. ok) return
        central(n + 1) = 0.0_dp
        do i = 1, n
            central(m - i + 1) = -central(i)
        end do

        d = sqrt(fac)
        do i = 1, m
            k = i - 1
            if (k == 0) then
                nodes(i) = 0.0_dp
            else
                nodes(i) = 1.0_dp / d ** real(k - 1, dp)
            end if
        end do
        do j = 1, m
            do i = 1, m
                a(j, i) = nodes(i) ** real(j - 1, dp)
            end do
        end do
        call solve_dense(a, b, right, ok)
    end subroutine derivative_weights

    subroutine solve_dense(a_in, b_in, x, ok)
        real(dp), intent(in) :: a_in(:, :) !! square coefficient matrix
        real(dp), intent(in) :: b_in(:) !! right-hand side vector
        real(dp), intent(out) :: x(:) !! solution vector
        logical, intent(out) :: ok !! true when a finite nonsingular solution was obtained

        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:)
        real(dp) :: factor
        real(dp) :: pivot
        real(dp) :: scale
        real(dp) :: temp
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: p

        ok = .false.
        x = 0.0_dp
        n = size(b_in)
        if (n < 1 .or. size(a_in, 1) /= n .or. size(a_in, 2) /= n .or. size(x) /= n) return
        if (any(.not. ieee_is_finite(a_in)) .or. any(.not. ieee_is_finite(b_in))) return

        allocate(a(n, n), b(n))
        a = a_in
        b = b_in
        do k = 1, n - 1
            p = k
            pivot = abs(a(k, k))
            do i = k + 1, n
                if (abs(a(i, k)) > pivot) then
                    pivot = abs(a(i, k))
                    p = i
                end if
            end do
            scale = max(1.0_dp, maxval(abs(a(k:n, k:n))))
            if (pivot <= epsilon(1.0_dp) * scale) return
            if (p /= k) then
                do j = k, n
                    temp = a(k, j)
                    a(k, j) = a(p, j)
                    a(p, j) = temp
                end do
                temp = b(k)
                b(k) = b(p)
                b(p) = temp
            end if
            do i = k + 1, n
                factor = a(i, k) / a(k, k)
                a(i, k) = 0.0_dp
                a(i, k + 1:n) = a(i, k + 1:n) - factor * a(k, k + 1:n)
                b(i) = b(i) - factor * b(k)
            end do
        end do
        scale = max(1.0_dp, maxval(abs(a)))
        if (abs(a(n, n)) <= epsilon(1.0_dp) * scale) return
        do i = n, 1, -1
            temp = b(i)
            if (i < n) temp = temp - dot_product(a(i, i + 1:n), x(i + 1:n))
            if (abs(a(i, i)) <= epsilon(1.0_dp) * scale) return
            x(i) = temp / a(i, i)
        end do
        if (any(.not. ieee_is_finite(x))) return
        ok = .true.
    end subroutine solve_dense

    subroutine initialize_derivative_result(res, x)
        type(derivative_result), intent(out) :: res !! result initialized to invalid-input state
        real(dp), intent(in) :: x !! requested differentiation point

        res%success = .false.
        res%status = differentiate_status_invalid_input
        res%df = quiet_nan()
        res%error = quiet_nan()
        res%nit = 0
        res%nfev = 0
        res%x = x
    end subroutine initialize_derivative_result

    subroutine allocate_jacobian_result(res, m, n)
        type(jacobian_result), intent(out) :: res !! result arrays initialized to invalid-input state
        integer, intent(in) :: m !! number of function outputs
        integer, intent(in) :: n !! number of function inputs

        allocate(res%success(m, n), res%status(m, n), res%df(m, n), res%error(m, n), &
            res%nit(m, n), res%nfev(m, n))
        res%success = .false.
        res%status = differentiate_status_invalid_input
        res%df = quiet_nan()
        res%error = quiet_nan()
        res%nit = 0
        res%nfev = 0
    end subroutine allocate_jacobian_result

    subroutine allocate_hessian_result(res, n)
        type(hessian_result), intent(out) :: res !! result arrays initialized to invalid-input state
        integer, intent(in) :: n !! number of input variables

        allocate(res%success(n, n), res%status(n, n), res%ddf(n, n), &
            res%error(n, n), res%nfev(n, n))
        res%success = .false.
        res%status = differentiate_status_invalid_input
        res%ddf = quiet_nan()
        res%error = quiet_nan()
        res%nfev = 0
    end subroutine allocate_hessian_result

    subroutine set_jacobian_nonfinite(res)
        type(jacobian_result), intent(inout) :: res !! result arrays set to nonfinite termination

        res%success = .false.
        res%status = differentiate_status_nonfinite
        res%df = quiet_nan()
        res%error = quiet_nan()
    end subroutine set_jacobian_nonfinite

    subroutine store_jacobian_element(res, i, j, one)
        type(jacobian_result), intent(inout) :: res !! Jacobian result receiving one scalar derivative
        integer, intent(in) :: i !! output row index
        integer, intent(in) :: j !! input column index
        type(derivative_result), intent(in) :: one !! scalar derivative result to store

        res%success(i, j) = one%success
        res%status(i, j) = one%status
        res%df(i, j) = one%df
        res%error(i, j) = one%error
        res%nit(i, j) = one%nit
        res%nfev(i, j) = one%nfev
    end subroutine store_jacobian_element

    pure function quiet_nan() result(value)
        real(dp) :: value

        value = ieee_value(0.0_dp, ieee_quiet_nan)
    end function quiet_nan

end module scifort_differentiate
