! SPDX-License-Identifier: MIT
! Compare bounded and unbounded raw Arnoldi product storage for LGMRES.
!
! Run under /usr/bin/time -f 'peak_kb=%M'.
! Optional arguments: system size, outer_k capacity, and raw cache cap (MiB).
! Avoid a dense Jacobian so memory measures mostly solver workspaces.
program benchmark_lgmres_candidate_cache_budget
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: inner_size = 20
    integer :: outer_size
    integer :: n, j, ios, jvp_calls, budget_mib
    logical :: use_preconditioner
    real(dp), allocatable :: initial(:), target(:), diagonal(:), coupling(:)
    type(vector_root_result) :: res
    character(len=40) :: arg

    n = 200000
    outer_size = 4
    budget_mib = 0
    use_preconditioner = .false.
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *, iostat=ios) n
        if (ios /= 0) error stop 'expected integer number of variables'
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *, iostat=ios) outer_size
        if (ios /= 0) error stop 'expected integer outer_k capacity'
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read(arg, *, iostat=ios) budget_mib
        if (ios /= 0) error stop 'third argument must be nonnegative cache MiB'
    end if
    if (n < 2 .or. n > 1000000) error stop 'n must be between 2 and 1000000'
    if (outer_size < 1 .or. outer_size > 128) error stop 'outer_k must be between 1 and 128'
    if (budget_mib < 0) error stop 'cache cap must be nonnegative'
    allocate(initial(n), target(n), diagonal(n), coupling(n))
    do j = 1, n
        target(j) = 1.0_dp + 0.2_dp * sin(0.1_dp * real(j, dp))
        initial(j) = target(j) + 0.9_dp + 0.1_dp * sin(0.17_dp * real(j, dp))
        diagonal(j) = 0.4_dp + real(mod(j, 97), dp) / 120.0_dp
        coupling(j) = sin(0.21_dp * real(j, dp))
    end do
    jvp_calls = 0
    if (use_preconditioner) then
        res = root(residual, initial, method='krylov', inner_method='lgmres', &
            jacobian_product=jvp, inner_m=preconditioner, inner_maxiter=inner_size, &
            outer_k=outer_size, outer_image_prune_auto=.true., &
            outer_image_reuse=.true., outer_candidate_image_reuse=.true., &
            outer_candidate_cache_max_mib=budget_mib, &
            inner_orthogonalization='selective', line_search='none', &
            fatol=1.0e-9_dp, maxiter=12)
    else
        res = root(residual, initial, method='krylov', inner_method='lgmres', &
            jacobian_product=jvp, inner_maxiter=inner_size, outer_k=outer_size, &
            outer_image_prune_auto=.true., outer_image_reuse=.true., &
            outer_candidate_image_reuse=.true., &
            outer_candidate_cache_max_mib=budget_mib, &
            inner_orthogonalization='selective', line_search='none', &
            fatol=1.0e-9_dp, maxiter=12)
    end if
    if (.not. res%success .or. .not. all(ieee_is_finite(res%x))) then
        print '(a,a)', 'status: ', trim(res%message)
        error stop 'nonlinear solve failed'
    end if
    if (maxval(abs(res%fun)) > 1.0e-9_dp) error stop 'residual too large'
    print '(a,i0,a,i0,a,i0,a,i0,a,es12.4)', &
        'n=', n, ' nit=', res%nit, ' nfev=', res%nfev, &
        ' njvp=', jvp_calls, ' max_residual=', maxval(abs(res%fun))

contains
    subroutine preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to transform
        real(dp), intent(out) :: product(:) !! approximate Jacobian inverse action
        logical, intent(out) :: ok !! validity flag

        product = vec / (2.0_dp * diagonal * initial)
        ok = all(ieee_is_finite(product))
    end subroutine preconditioner

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear variables
        real(dp), intent(out) :: fvec(:) !! nonlinear residual
        real(dp) :: delta(size(x))
        real(dp) :: projection

        delta = x - target
        projection = sum(coupling * delta) / real(n, dp)
        fvec = diagonal * (x*x - target*target) + 0.02_dp * projection * delta
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear variables
        real(dp), intent(in) :: vec(:) !! Jacobian product direction
        real(dp), intent(out) :: product(:) !! J(x)*vec
        logical, intent(out) :: ok !! validity flag
        real(dp) :: delta(size(x))
        real(dp) :: projection, vprojection

        jvp_calls = jvp_calls + 1
        delta = x - target
        projection = sum(coupling * delta) / real(n, dp)
        vprojection = sum(coupling * vec) / real(n, dp)
        product = 2.0_dp * diagonal * x * vec + &
            0.02_dp * (projection * vec + vprojection * delta)
        ok = all(ieee_is_finite(product))
    end subroutine jvp
end program benchmark_lgmres_candidate_cache_budget
