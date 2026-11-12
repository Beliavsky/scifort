! SPDX-License-Identifier: MIT
! End-to-end nonlinear LGMRES image-basis reuse benchmark.
! Run with /usr/bin/time -f 'peak_kb=%M' to measure peak process RSS.
! Arguments: problem size, repetitions, and basis-reuse flag (0 or 1).
program benchmark_lgmres_inplace_image_basis
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer :: n, repetitions, j, trial, ios, jvp_calls, total_jvp, total_fev
    integer :: basis_flag, saved_nit
    logical :: reuse_basis
    real(dp), allocatable :: initial(:), target(:), diagonal(:), coupling(:)
    real(dp) :: before, elapsed, checksum, max_residual
    type(vector_root_result) :: res
    character(len=40) :: argument

    n = 30000
    repetitions = 3
    basis_flag = 1
    if (command_argument_count() >= 1) then
        call get_command_argument(1, argument)
        read(argument, *, iostat=ios) n
        if (ios /= 0) error stop 'invalid problem size'
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, argument)
        read(argument, *, iostat=ios) repetitions
        if (ios /= 0) error stop 'invalid repetitions'
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, argument)
        read(argument, *, iostat=ios) basis_flag
        if (ios /= 0) error stop 'invalid reuse flag'
    end if
    if (n < 2 .or. n > 1000000) error stop 'n out of range'
    if (repetitions < 1 .or. repetitions > 100) error stop 'repetitions out of range'
    if (basis_flag /= 0 .and. basis_flag /= 1) error stop 'reuse flag must be 0 or 1'
    reuse_basis = basis_flag == 1
    allocate(initial(n), target(n), diagonal(n), coupling(n))
    do j = 1, n
        target(j) = 1.0_dp + 0.2_dp * sin(0.1_dp * real(j, dp))
        initial(j) = target(j) + 0.9_dp + 0.1_dp * sin(0.17_dp * real(j, dp))
        diagonal(j) = 0.4_dp + real(mod(j, 97), dp) / 120.0_dp
        coupling(j) = sin(0.21_dp * real(j, dp))
    end do

    checksum = 0.0_dp
    max_residual = 0.0_dp
    total_jvp = 0
    total_fev = 0
    saved_nit = -1
    call cpu_time(before)
    do trial = 1, repetitions
        jvp_calls = 0
        res = root(residual, initial, method='krylov', inner_method='lgmres', &
            jacobian_product=jvp, inner_maxiter=10, outer_k=64, &
            outer_image_prune_auto=.true., outer_image_basis_reuse=reuse_basis, &
            inner_orthogonalization='selective', line_search='none', &
            fatol=1.0e-9_dp, maxiter=12)
        if (.not. res%success .or. .not. all(ieee_is_finite(res%x))) then
            print '(a)', trim(res%message)
            error stop 'nonlinear solve failed'
        end if
        max_residual = max(max_residual, maxval(abs(res%fun)))
        if (max_residual > 1.0e-9_dp) error stop 'residual too large'
        if (saved_nit < 0) saved_nit = res%nit
        if (res%nit /= saved_nit) error stop 'iteration count changed across runs'
        checksum = checksum + sum(res%x) / real(n, dp)
        total_jvp = total_jvp + jvp_calls
        total_fev = total_fev + res%nfev
    end do
    call cpu_time(elapsed)
    elapsed = elapsed - before
    print '(a,i0,a,i0,a,l1)', 'n=', n, ' reps=', repetitions, ' basis_reuse=', reuse_basis
    print '(a,i0,a,i0,a,i0)', 'nit=', saved_nit, ' nfev=', total_fev, ' njvp=', total_jvp
    print '(a,es24.15,a,es12.4)', 'checksum=', checksum, ' max_residual=', max_residual
    print '(a,f12.6)', 'cpu_seconds=', elapsed
contains
    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear variables
        real(dp), intent(out) :: fvec(:) !! nonlinear residual
        real(dp) :: delta(size(x)), projection

        delta = x - target
        projection = sum(coupling * delta) / real(n, dp)
        fvec = diagonal * (x*x - target*target) + 0.02_dp * projection * delta
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! linearization point
        real(dp), intent(in) :: vec(:) !! Jacobian-vector direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector result
        logical, intent(out) :: ok !! validity flag
        real(dp) :: delta(size(x)), projection, vprojection

        jvp_calls = jvp_calls + 1
        delta = x - target
        projection = sum(coupling * delta) / real(n, dp)
        vprojection = sum(coupling * vec) / real(n, dp)
        product = 2.0_dp * diagonal * x * vec + &
            0.02_dp * (projection * vec + vprojection * delta)
        ok = all(ieee_is_finite(product))
    end subroutine jvp
end program benchmark_lgmres_inplace_image_basis
