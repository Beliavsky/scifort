! SPDX-License-Identifier: MIT
! Native deterministic DIRECT regression and invalid-input tests.
program optimize_direct_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_direct, only : direct, direct_result, &
        direct_maxfun, direct_maxiter, direct_invalid, direct_target, &
        direct_callback_stop, direct_numerical, direct_vol_tol, direct_len_tol
    implicit none
    type(direct_result) :: result, again
    real(dp) :: b(2,2), invalid_b(2,2), shift(2)
    integer :: i, j, checks, callback_calls
    logical :: stop_callback
    checks = 0
    callback_calls = 0
    stop_callback = .false.
    b(:,1) = [-2.0_dp,-2.0_dp]
    b(:,2) = [ 2.0_dp, 2.0_dp]

    do i = 0, 1
        result = direct(quadratic,b,maxfun=1400,maxiter=120, &
            locally_biased=(i==1),vol_tol=0.0_dp,len_tol=0.0_dp)
        call check(allocated(result%x),'allocated solution')
        call check(result%fun < 1.0e-9_dp,'quadratic minimization')
        call check(result%nfev <= 1400,'strict evaluation cap')
        call check(result%nfev > 100,'meaningful evaluations')
        call check(result%nit > 0 .and. result%nit <= 120,'iteration accounting')
        call check(result%nrect == result%nfev,'one center per function call')
        call check(all(result%x >= b(:,1)) .and. all(result%x <= b(:,2)), 'bounds')
        call check(result%status == direct_maxfun .or. &
            result%status == direct_maxiter, 'budget termination')
        again = direct(quadratic,b,maxfun=1400,maxiter=120, &
            locally_biased=(i==1),vol_tol=0.0_dp,len_tol=0.0_dp)
        call check(all(abs(result%x-again%x) <= 0.0_dp) .and. &
            abs(result%fun-again%fun) <= 0.0_dp, &
            'DIRECT must be reproducible without a random seed')
    end do


    ! Repeated shifted convex bowls ensure that the solver explores the
    ! domain instead of merely returning the initial center or a boundary.
    do i = 1, 4
        do j = 1, 4
            shift = [real(i,dp)*0.37_dp-1.0_dp, real(j,dp)*0.41_dp-1.0_dp]
            result = direct(shifted_quadratic,b,maxfun=800,maxiter=120, &
                locally_biased=mod(i+j,2)==0,vol_tol=0.0_dp,len_tol=0.0_dp)
            call check(allocated(result%x),'shifted objective has solution')
            call check(result%fun < 1.0e-7_dp,'shifted quadratic finds optimum')
            call check(result%nfev <= 800 .and. result%nrect == result%nfev, &
                'shifted quadratic respects budget and rectangle accounting')
            call check(all(result%x >= b(:,1)) .and. all(result%x <= b(:,2)), &
                'shifted quadratic respects bounds')
        end do
    end do

    result = direct(quadratic, b, maxfun=1)
    call check(result%nfev == 1,'maxfun=1')
    call check(result%status == direct_maxfun,'maxfun=1 status')
    result = direct(quadratic, b, maxfun=2)
    call check(result%nfev == 1,'no unpaired samples')
    result = direct(quadratic, b, maxfun=3)
    call check(result%nfev == 3,'one trisection consumes two evaluations')
    result = direct(quadratic, b, maxfun=41)
    call check(result%nfev <= 41,'odd evaluation limit respected')
    result = direct(quadratic, b, maxiter=0)
    call check(result%nit == 0 .and. result%nfev == 1, 'maxiter zero')
    call check(result%status == direct_maxiter,'maxiter status')

    ! Constant zero objective reaches supplied known optimum without search.
    result = direct(constant_zero,b,maxfun=100,f_min=0.0_dp)
    call check(result%status == direct_target .and. result%success,'known target')
    call check(result%nfev == 1,'target at center stops after one evaluation')
    result = direct(constant_zero,b,maxfun=100,vol_tol=1.0_dp)
    call check(result%status == direct_vol_tol .and. result%success,'volume termination')
    result = direct(constant_zero,b,maxfun=100,len_tol=1.0_dp)
    call check(result%status == direct_len_tol .and. result%success,'length termination')

    stop_callback = .true.
    result = direct(quadratic,b,maxfun=101,callback=cb)
    call check(result%status == direct_callback_stop,'callback stop')
    call check(callback_calls == 1,'callback invoked once')
    call check(result%nfev > 1,'callback only after a full iteration')
    stop_callback = .false.

    invalid_b = b
    invalid_b(1,2) = invalid_b(1,1)
    result = direct(quadratic,invalid_b)
    call check(result%status == direct_invalid,'zero-width interval rejected')
    call check(result%nfev == 0,'invalid domain does not call objective')
    result = direct(quadratic,b(:,1:1))
    call check(result%status == direct_invalid,'shape mismatch')
    result = direct(quadratic,b,maxfun=0)
    call check(result%status == direct_invalid,'maxfun invalid')
    result = direct(quadratic,b,maxiter=-1)
    call check(result%status == direct_invalid,'negative maxiter')
    result = direct(quadratic,b,eps=-0.1_dp)
    call check(result%status == direct_invalid,'negative eps')
    result = direct(quadratic,b,len_tol=2.0_dp)
    call check(result%status == direct_invalid,'len_tol out of range')
    result = direct(quadratic,b,vol_tol=2.0_dp)
    call check(result%status == direct_invalid,'vol_tol out of range')
    result = direct(quadratic,b,f_min_rtol=-1.0_dp)
    call check(result%status == direct_invalid,'f_min_rtol out of range')
    invalid_b = b
    invalid_b(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    result = direct(quadratic,invalid_b)
    call check(result%status == direct_invalid,'NaN bound rejected')
    invalid_b = b
    invalid_b(1,2) = ieee_value(0.0_dp,ieee_positive_inf)
    result = direct(quadratic,invalid_b)
    call check(result%status == direct_invalid,'infinite bound rejected')
    result = direct(infinite_objective,b,maxfun=100)
    call check(result%status == direct_numerical,'nonfinite objective rejected')
    call check(result%nfev == 1,'initial failure counted once')

    write (*,'(a,i0)') 'native DIRECT checks passed: ', checks
contains
    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! argument of length two
        real(dp) :: f !! quadratic bowl
        f = (x(1)-0.43_dp)**2 + 3.0_dp*(x(2)+0.61_dp)**2
    end function quadratic
    function shifted_quadratic(x) result(f)
        real(dp), intent(in) :: x(:) !! point to evaluate
        real(dp) :: f !! convex bowl shifted inside the search box
        f = sum((x-shift)**2)
    end function shifted_quadratic
    function constant_zero(x) result(f)
        real(dp), intent(in) :: x(:) !! unused argument
        real(dp) :: f !! constant objective
        f = 0.0_dp
        if (size(x) /= 2) error stop 'unexpected test objective shape'
    end function constant_zero
    function infinite_objective(x) result(f)
        real(dp), intent(in) :: x(:) !! unused argument
        real(dp) :: f !! deliberately infinite objective
        f = ieee_value(0.0_dp,ieee_positive_inf)
        if (size(x) /= 2) error stop 'unexpected test objective shape'
    end function infinite_objective
    subroutine cb(x, value, stop)
        real(dp), intent(in) :: x(:) !! current best point
        real(dp), intent(in) :: value !! current best objective
        logical, intent(out) :: stop !! request early termination
        if (size(x) /= 2 .or. value < 0.0_dp) error stop 'bad callback arguments'
        callback_calls = callback_calls + 1
        stop = stop_callback
    end subroutine cb
    subroutine check(ok, message)
        logical, intent(in) :: ok !! test condition
        character(len=*), intent(in) :: message !! assertion label
        if (.not. ok) then
            print *, 'FAILED: ', message
            error stop 1
        end if
        checks = checks+1
    end subroutine check
end program optimize_direct_reference
