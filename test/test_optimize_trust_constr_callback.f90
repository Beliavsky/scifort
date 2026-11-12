! SPDX-License-Identifier: MIT
! Iteration callbacks, user termination, and barrier slack hiding for trust-constr.
program test_optimize_trust_constr_callback
    use scifort_kinds, only : dp
    use scifort_optimize_minimize, only : minimize, trust_constr, nonlinear_constraint, &
        multivariate_minimize_result
    implicit none

    type(multivariate_minimize_result) :: r
    type(nonlinear_constraint) :: blocks(1)
    integer :: calls
    integer :: failures
    integer :: previous_nit
    integer :: stop_at
    logical :: bad_shape
    logical :: bad_iteration
    real(dp) :: last_x(2)
    real(dp) :: x0(2)

    failures = 0
    x0 = [3.0_dp, 4.0_dp]

    ! Callback receives the initial state as nit=1; a true stop returns status 3.
    call reset_log(1)
    r = trust_constr(objective, x0, jac=gradient, hess=hessian, &
        callback=stop_callback, maxiter=100)
    call check_stopped('unconstrained initial', r, 1, x0)

    ! The callback is invoked again on accepted/rejected SQP iterations.
    call reset_log(3)
    r = minimize(objective, x0, method='trust-constr', jac=gradient, &
        hess=hessian, callback=stop_callback, maxiter=100)
    call check_stopped('unconstrained minimize', r, 3)

    ! The equality-SQP solver forwards the callback with the raw variables.
    call reset_log(2)
    r = trust_constr(objective, [0.0_dp, 0.0_dp], jac=gradient, &
        hess=hessian, constraints=equality, constraint_jac=equality_jac, &
        meq=1, callback=stop_callback, maxiter=100)
    call check_stopped('equality', r, 2)

    ! Barrier internal variables include slacks; callback must see only two
    ! physical variables, both in the direct and the generic entry point.
    call reset_log(2)
    r = trust_constr(objective, [0.5_dp, 0.4_dp], jac=gradient, hess=hessian, &
        lower=[0.0_dp, 0.0_dp], callback=stop_callback, maxiter=100)
    call check_stopped('bounds interior point', r, 2)
    if (trim(r%method) /= 'tr_interior_point') call fail('barrier method lost')

    blocks(1)%fun => equality
    blocks(1)%jac => equality_jac
    blocks(1)%lower = [1.0_dp]
    blocks(1)%upper = [1.0_dp]
    call reset_log(2)
    r = minimize(objective, [0.5_dp, 0.4_dp], constraint_objects=blocks, &
        jac=gradient, hess=hessian, lower=[0.0_dp, 0.0_dp], &
        callback=stop_callback, maxiter=100)
    call check_stopped('constraint objects with bounds', r, 2)
    if (.not. allocated(r%constr_nfev)) call fail('object counts lost')

    ! Hessian-vector products take a recursive route through trust_constr.
    call reset_log(2)
    r = trust_constr(objective, [3.0_dp, 4.0_dp], jac=gradient, &
        hessp=hessian_product, callback=stop_callback, maxiter=100)
    call check_stopped('recursive hessp dispatch', r, 2)

    ! An absent callback must preserve the original convergence behavior.
    r = trust_constr(objective, [3.0_dp, 4.0_dp], jac=gradient, hess=hessian)
    if (.not. r%success .or. maxval(abs(r%x - [1.0_dp, 2.0_dp])) > 1.0e-7_dp) &
        call fail('no-callback convergence regression')
    call reset_log(100)
    r = trust_constr(objective, [3.0_dp, 4.0_dp], jac=gradient, &
        hess=hessian, callback=stop_callback)
    if (.not. r%success) call fail('nonstopping callback suppressed convergence')
    if (calls /= r%nit .or. previous_nit /= r%nit) &
        call fail('callback called wrong number of times on convergence')

    ! Verify that callback iteration indices do not restart at successive
    ! barrier stages, and that allowing all iterations preserves convergence.
    call reset_log(10000)
    r = trust_constr(objective, [0.5_dp, 0.4_dp], jac=gradient, hess=hessian, &
        lower=[0.0_dp, 0.0_dp], callback=stop_callback, maxiter=200)
    if (.not. r%success) call fail('nonstopping barrier callback did not converge')
    if (calls /= r%nit .or. bad_iteration .or. bad_shape) &
        call fail('barrier callback iteration accounting')
    if (maxval(abs(r%x - last_x)) > 1.0e-14_dp) &
        call fail('nonstopping barrier final callback point differs')

    ! This keyword must not silently do nothing for other minimizers.
    call reset_log(1)
    r = minimize(objective, x0, method='bfgs', jac=gradient, callback=stop_callback)
    if (r%status /= -1 .or. r%success .or. calls /= 0) &
        call fail('non trust-constr accepted callback')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_trust_constr_callback: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_trust_constr_callback: PASS'

contains

    subroutine reset_log(nstop)
        integer, intent(in) :: nstop !! iteration at which to stop
        calls = 0
        previous_nit = 0
        stop_at = nstop
        bad_shape = .false.
        bad_iteration = .false.
        last_x = 0.0_dp
    end subroutine reset_log

    subroutine stop_callback(x, nit, stop)
        real(dp), intent(in) :: x(:) !! physical variable vector
        integer, intent(in) :: nit !! cumulative iteration number
        logical, intent(out) :: stop !! true to request termination

        if (size(x) /= 2) then
            bad_shape = .true.
        else
            last_x = x
        end if
        calls = calls + 1
        if (nit /= previous_nit + 1) bad_iteration = .true.
        previous_nit = nit
        stop = nit >= stop_at
    end subroutine stop_callback

    subroutine check_stopped(label, r, expected_nit, expected_x)
        character(len=*), intent(in) :: label !! diagnostic label
        type(multivariate_minimize_result), intent(in) :: r !! optimization result
        integer, intent(in) :: expected_nit !! callback termination iteration
        real(dp), intent(in), optional :: expected_x(:) !! expected early point

        if (r%status /= 3 .or. r%success) then
            print '(a,1x,a,1x,i0)', 'bad callback status:', label, r%status
            call fail('callback status/termination')
        end if
        if (r%nit /= expected_nit .or. calls /= expected_nit) &
            call fail('callback iteration/count')
        if (bad_shape .or. bad_iteration) call fail('callback shape/sequence')
        if (.not. allocated(r%x)) then
            call fail('missing stopped point')
        else
            if (maxval(abs(r%x - last_x)) > 1.0e-14_dp) &
                call fail('returned point does not match callback')
            if (present(expected_x)) then
                if (maxval(abs(r%x - expected_x)) > 0.0_dp) &
                    call fail('initial callback changed x')
            end if
        end if
    end subroutine check_stopped

    subroutine fail(reason)
        character(len=*), intent(in) :: reason !! description of failure
        failures = failures + 1
        print '(a)', '  FAIL: '//reason
    end subroutine fail

    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp) :: f !! quadratic objective
        f = sum((x - [1.0_dp, 2.0_dp])**2)
    end function objective

    subroutine gradient(x, g)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: g(:) !! objective gradient
        g = 2.0_dp * (x - [1.0_dp, 2.0_dp])
    end subroutine gradient

    subroutine hessian(x, h)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: h(:, :) !! exact Hessian
        h = 0.0_dp
        h(1, 1) = 2.0_dp
        h(2, 2) = 2.0_dp
    end subroutine hessian

    subroutine hessian_product(x, p, hp)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(in) :: p(:) !! direction
        real(dp), intent(out) :: hp(:) !! Hessian times direction
        hp = 2.0_dp * p
    end subroutine hessian_product

    subroutine equality(x, c)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: c(:) !! constraint output
        c(1) = sum(x)
    end subroutine equality

    subroutine equality_jac(x, a)
        real(dp), intent(in) :: x(:) !! optimization variables
        real(dp), intent(out) :: a(:, :) !! constraint Jacobian
        a(1, :) = 1.0_dp
    end subroutine equality_jac
end program test_optimize_trust_constr_callback
