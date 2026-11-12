! SPDX-License-Identifier: MIT
program test_optimize_differential_evolution_reference
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_differential_evolution
    implicit none

    type(differential_evolution_result) :: r, s
    character(len=20), parameter :: strategies(12) = [character(len=20) :: &
        'best1bin', 'rand1bin', 'best2bin', 'rand2bin', &
        'randtobest1bin', 'currenttobest1bin', &
        'best1exp', 'rand1exp', 'best2exp', 'rand2exp', &
        'randtobest1exp', 'currenttobest1exp']
    real(dp) :: bounds(3,2), fixed(3,2), one(1,2), errors(12)
    integer :: i, j, k, count_evals, checks, strata(18)
    checks = 0
    count_evals = 0
    bounds(:,1) = [-5.0_dp,-3.0_dp,-2.0_dp]
    bounds(:,2) = [ 5.0_dp, 4.0_dp, 6.0_dp]
    do i = 1,12
        r = differential_evolution(shifted_sphere, bounds, strategy=trim(strategies(i)), &
                                   maxiter=220, popsize=10, seed=2027+i, updating='immediate')
        write(*,'(a,1x,i0,1x,es12.4,1x,i0,1x,i0)') trim(strategies(i)), r%status, r%fun, r%nit, r%nfev
        call must(r%status /= de_status_invalid, 'valid strategy '//trim(strategies(i)))
        call must(r%nfev <= 30*221, 'evaluation limit')
        call must(r%fun < 0.025_dp, 'optimum quality '//trim(strategies(i)))
        call must(all(r%x >= bounds(:,1)) .and. all(r%x <= bounds(:,2)), 'best in bounds')
        call must(all(r%population >= spread(bounds(:,1),2,size(r%population,2))), 'lower bounds')
        call must(all(r%population <= spread(bounds(:,2),2,size(r%population,2))), 'upper bounds')
        call must(abs(r%fun - minval(r%population_energies)) < 1.0e-13_dp, 'returned best is minimum')
        checks = checks + 6
        errors(i) = r%fun
    end do
    r = differential_evolution(shifted_sphere,bounds,seed=47,maxiter=25,popsize=6, &
                               updating='deferred',init='random',mutation=0.65_dp)
    s = differential_evolution(shifted_sphere,bounds,seed=47,maxiter=25,popsize=6, &
                               updating='deferred',init='random',mutation=0.65_dp)
    call must(all(abs(r%population - s%population) <= 0.0_dp), 'seed reproducibility population')
    call must(all(abs(r%population_energies - s%population_energies) <= 0.0_dp), 'seed reproducibility objective')
    call must(r%nfev == s%nfev .and. r%nit == s%nit, 'seed reproducibility counters')
    checks = checks + 3
    fixed = bounds
    fixed(2,:) = 1.25_dp
    r = differential_evolution(shifted_sphere,fixed,maxiter=180,seed=932,popsize=12)
    call must(all(abs(r%population(2,:) - 1.25_dp) <= 0.0_dp), 'fixed coordinates remain fixed')
    call must(r%fun < 1.0e-4_dp, 'fixed-bound optimum')
    checks = checks + 2
    fixed(:,1) = [1.0_dp,1.25_dp,2.0_dp]
    fixed(:,2) = fixed(:,1)
    r = differential_evolution(shifted_sphere,fixed,maxiter=100)
    call must(r%success.and.r%nit==0.and.r%nfev==1, 'all fixed early return')
    call must(maxval(abs(r%x-fixed(:,1))) < 1.0e-14_dp, 'all fixed solution')
    checks = checks + 2
    r = differential_evolution(shifted_sphere,bounds,maxiter=0,popsize=6,seed=11, &
                               init='latinhypercube')
    do j = 1,3
        strata = 0
        do i = 1,size(r%population,2)
            k = 1 + int(real(size(r%population,2),dp)* &
                (r%population(j,i)-bounds(j,1))/(bounds(j,2)-bounds(j,1)))
            k = min(size(strata),max(1,k))
            strata(k) = strata(k) + 1
        end do
        call must(all(strata == 1), 'Latin hypercube covers each stratum exactly once')
        checks = checks + 1
    end do
    r = differential_evolution(shifted_sphere,bounds,maxiter=0,popsize=6,seed=11, &
                               x0=[1.0_dp,1.25_dp,2.0_dp])
    call must(r%fun <= 0.0_dp.and.all(abs(r%x-[1.0_dp,1.25_dp,2.0_dp]) <= 0.0_dp), &
              'valid x0 replaces initial population member')
    checks = checks + 1
    r = differential_evolution(shifted_sphere,bounds,maxiter=0,popsize=5,seed=1)
    call must(r%nfev==15.and.r%nit==0.and.r%status==de_status_maxiter, 'zero maxiter evaluation count')
    checks = checks + 1
    count_evals = 0
    r = differential_evolution(counted_sphere,bounds,maxiter=100,popsize=5,seed=9, &
                               maxfun=22,updating='deferred')
    call must(r%nfev==22.and.count_evals==22, 'exact maxfun evaluation budget')
    call must(r%status==de_status_maxfun.and.r%nit==1, 'maxfun termination')
    checks = checks + 2
    r = differential_evolution(shifted_sphere,bounds,maxiter=100,popsize=6,seed=1, &
                               callback=stop_at_three)
    call must(r%status==de_status_callback.and.r%nit==3, 'callback stop')
    checks = checks + 1
    r = differential_evolution(shifted_sphere,bounds,strategy='foo')
    call must(r%status==de_status_invalid.and..not.allocated(r%x), 'reject invalid strategy')
    r = differential_evolution(shifted_sphere,bounds,popsize=0)
    call must(r%status==de_status_invalid, 'reject zero popsize')
    r = differential_evolution(shifted_sphere,bounds,maxiter=-1)
    call must(r%status==de_status_invalid, 'reject negative maxiter')
    r = differential_evolution(shifted_sphere,bounds,tol=-1.0_dp)
    call must(r%status==de_status_invalid, 'reject negative tol')
    r = differential_evolution(shifted_sphere,bounds,mutation=2.0_dp)
    call must(r%status==de_status_invalid, 'reject invalid mutation')
    r = differential_evolution(shifted_sphere,bounds,mutation_range=[1.0_dp,0.5_dp])
    call must(r%status==de_status_invalid, 'reject inverted dithering interval')
    r = differential_evolution(shifted_sphere,bounds,recombination=1.2_dp)
    call must(r%status==de_status_invalid, 'reject invalid crossover')
    r = differential_evolution(shifted_sphere,bounds,x0=[1.0_dp])
    call must(r%status==de_status_invalid, 'reject invalid x0 shape')
    r = differential_evolution(shifted_sphere,bounds,x0=[100.0_dp,1.0_dp,2.0_dp])
    call must(r%status==de_status_invalid, 'reject out-of-bounds x0')
    r = differential_evolution(shifted_sphere,bounds,init='sobol')
    call must(r%status==de_status_invalid, 'reject unsupported initialization')
    r = differential_evolution(shifted_sphere,bounds,updating='wrong')
    call must(r%status==de_status_invalid, 'reject bad updating')
    fixed = bounds
    fixed(1,1) = 10.0_dp
    r = differential_evolution(shifted_sphere,fixed)
    call must(r%status==de_status_invalid, 'reject backwards bounds')
    r = differential_evolution(shifted_sphere,bounds,maxiter=5,popsize=6,seed=17,maxfun=18)
    call must(r%nfev==18.and.r%nit==0.and.r%status==de_status_maxfun, &
              'budget equal to initial population')
    checks = checks + 1
    r = differential_evolution(always_nonfinite,bounds,maxiter=5,seed=11)
    call must(.not.r%success.and.r%status==de_status_maxiter, &
              'all-NaN objective must not converge')
    checks = checks + 1
    r = differential_evolution(shifted_sphere,bounds,maxfun=4)
    call must(r%status==de_status_invalid, 'reject insufficient initial evaluation budget')
    checks = checks + 13
    one(1,:) = [-5.0_dp,5.0_dp]
    r = differential_evolution(single_sphere,one,maxiter=180,seed=22)
    call must(abs(r%x(1)-1.75_dp) < 0.015_dp, 'one-dimensional optimum')
    checks = checks + 1
    r = differential_evolution(occasionally_nonfinite,bounds,maxiter=150,seed=100)
    call must(r%fun < 0.01_dp, 'nonfinite objective penalized')
    checks = checks + 1
    print '(a,i0)', 'NATIVE CHECKS PASSED: ',checks
contains
    real(dp) function shifted_sphere(x) result(f)
        real(dp), intent(in) :: x(:) !! input vector, shape (3)
        f = (x(1)-1.0_dp)**2 + (x(2)-1.25_dp)**2 + (x(3)-2.0_dp)**2
    end function shifted_sphere
    real(dp) function counted_sphere(x) result(f)
        real(dp), intent(in) :: x(:) !! input vector, shape (3)
        count_evals = count_evals + 1
        f = shifted_sphere(x)
    end function counted_sphere
    real(dp) function single_sphere(x) result(f)
        real(dp), intent(in) :: x(:) !! input vector, shape (1)
        f = (x(1)-1.75_dp)**2
    end function single_sphere
    real(dp) function always_nonfinite(x) result(f)
        real(dp), intent(in) :: x(:) !! ignored candidate, shape (3)
        f = ieee_value(0.0_dp,ieee_quiet_nan)
        if (size(x) < 1) error stop 'bad candidate length'
    end function always_nonfinite
    real(dp) function occasionally_nonfinite(x) result(f)
        real(dp), intent(in) :: x(:) !! input vector, shape (3)
        if (x(1) < 0.0_dp) then
            f = ieee_value(0.0_dp,ieee_quiet_nan)
        else
            f = shifted_sphere(x)
        end if
    end function occasionally_nonfinite
    subroutine stop_at_three(x,fun,generation,stop)
        real(dp), intent(in) :: x(:) !! current minimizer
        real(dp), intent(in) :: fun !! objective at minimizer
        integer, intent(in) :: generation !! generation counter
        logical, intent(out) :: stop !! stop on generation 3
        if (size(x) /= 3 .or. .not. (fun >= 0.0_dp)) error stop 'bad callback inputs'
        stop = generation >= 3
    end subroutine stop_at_three
    subroutine must(test,msg)
        logical, intent(in) :: test !! assertion to verify
        character(len=*), intent(in) :: msg !! diagnostic message
        if (.not.test) then
            print '(a)', 'FAILED: '//msg
            error stop 1
        end if
    end subroutine must
end program test_optimize_differential_evolution_reference
