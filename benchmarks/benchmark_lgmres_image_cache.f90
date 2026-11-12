! SPDX-License-Identifier: MIT
! Deterministic scaling benchmark for LGMRES Jacobian-image reuse.
!
! Build against the SciFort project and run this executable. Results depend
! on compiler, CPU and BLAS/matmul implementation. No benchmark ranking is
! asserted by the regression suite.
program benchmark_lgmres_image_cache
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: nsizes(3) = [32, 64, 128]
    integer, parameter :: nmethods = 4
    integer, parameter :: repeats = 30
    real(dp), allocatable :: matrix(:, :)
    real(dp), allocatable :: target(:)
    real(dp), allocatable :: initial(:)
    real(dp), allocatable :: coupling(:)
    real(dp), allocatable :: matleft(:, :)
    real(dp), allocatable :: matright(:, :)
    real(dp), allocatable :: diagonal(:, :)
    real(dp), allocatable :: vl(:)
    real(dp), allocatable :: vr(:)
    real(dp) :: start_time, end_time, elapsed
    real(dp) :: max_residual
    type(vector_root_result) :: res
    integer :: i, n, which, repeat_index, jvp_calls, tracked_calls
    integer :: failures
    character(len=24) :: name

    failures = 0
    print '(a)', 'n,policy,success,nit,nfev,njvp,max_residual,cpu_seconds'
    do i = 1, size(nsizes)
        n = nsizes(i)
        call setup_problem(n)
        do which = 1, nmethods
            elapsed = 0.0_dp
            do repeat_index = 1, repeats
                jvp_calls = 0
                call cpu_time(start_time)
                select case (which)
                case (1)
                    name = 'auto'
                    res = root(residual, initial, method='krylov', &
                        inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                        jacobian_product=jvp, outer_image_prune_auto=.true., &
                        inner_orthogonalization='selective', line_search='armijo', &
                        fatol=1.0e-8_dp, maxiter=80)
                case (2)
                    name = 'reuse-retained'
                    res = root(residual, initial, method='krylov', &
                        inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                        jacobian_product=jvp, outer_image_prune_auto=.true., &
                        outer_image_reuse=.true., &
                        inner_orthogonalization='selective', line_search='armijo', &
                        fatol=1.0e-8_dp, maxiter=80)
                case (3)
                    name = 'reuse-candidate'
                    res = root(residual, initial, method='krylov', &
                        inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                        jacobian_product=jvp, outer_image_prune_auto=.true., &
                        outer_candidate_image_reuse=.true., &
                        inner_orthogonalization='selective', line_search='armijo', &
                        fatol=1.0e-8_dp, maxiter=80)
                case (4)
                    name = 'reuse-both'
                    res = root(residual, initial, method='krylov', &
                        inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                        jacobian_product=jvp, outer_image_prune_auto=.true., &
                        outer_image_reuse=.true., outer_candidate_image_reuse=.true., &
                        inner_orthogonalization='selective', line_search='armijo', &
                        fatol=1.0e-8_dp, maxiter=80)
                end select
                call cpu_time(end_time)
                elapsed = elapsed + (end_time - start_time)
                tracked_calls = jvp_calls
                if (.not. allocated(res%x) .or. .not. allocated(res%fun)) then
                    failures = failures + 1
                    cycle
                end if
                if (.not. all(ieee_is_finite(res%x))) failures = failures + 1
                max_residual = maxval(abs(res%fun))
                if (res%success .and. max_residual > 1.01e-8_dp) failures = failures + 1
            end do
            write (*, '(i0,a,a,a,l1,a,i0,a,i0,a,i0,a,es13.5,a,f10.6)') &
                n, ',', trim(name), ',', res%success, ',', res%nit, ',', &
                res%nfev, ',', tracked_calls, ',', max_residual, ',', &
                elapsed / real(repeats, dp)
        end do
        deallocate(matrix, target, initial, coupling, matleft, matright, diagonal, vl, vr)
    end do
    if (failures /= 0) then
        print '(a,i0)', 'benchmark_lgmres_image_cache: invalid results ', failures
        stop 1
    end if

contains
    subroutine setup_problem(nsize)
        integer, intent(in) :: nsize !! size of square nonlinear system
        integer :: i, j

        allocate(matrix(nsize, nsize), target(nsize), initial(nsize))
        allocate(coupling(nsize), matleft(nsize, nsize), matright(nsize, nsize))
        allocate(diagonal(nsize, nsize), vl(nsize), vr(nsize))
        do i = 1, nsize
            target(i) = 0.08_dp * sin(0.3_dp * real(i, dp))
            initial(i) = target(i) + 0.2_dp + 0.12_dp * sin(1.1_dp * real(i, dp))
            vl(i) = sin(0.73_dp * real(i, dp))
            vr(i) = cos(0.4_dp * real(i, dp))
            coupling(i) = sin(0.2_dp * real(i, dp) + 0.3_dp)
        end do
        vl = vl / norm2(vl)
        vr = vr / norm2(vr)
        coupling = coupling / norm2(coupling)
        diagonal = 0.0_dp
        do j = 1, nsize
            do i = 1, nsize
                matleft(i, j) = -2.0_dp * vl(i) * vl(j)
                matright(i, j) = -2.0_dp * vr(i) * vr(j)
            end do
            matleft(j, j) = matleft(j, j) + 1.0_dp
            matright(j, j) = matright(j, j) + 1.0_dp
            diagonal(j, j) = 1.0e3_dp ** &
                (-real(j - 1, dp) / real(nsize - 1, dp))
        end do
        matrix = matmul(matleft, matmul(diagonal, transpose(matright)))
    end subroutine setup_problem

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: fvec(:) !! residual
        real(dp) :: delta(size(x))
        real(dp) :: projection

        delta = x - target
        projection = dot_product(coupling, delta)
        fvec = matmul(matrix, delta + 0.3_dp * delta * delta + &
            0.03_dp * coupling * projection * projection)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! analytic Jacobian-vector product
        logical, intent(out) :: ok !! valid and finite
        real(dp) :: delta(size(x))
        real(dp) :: projected

        jvp_calls = jvp_calls + 1
        delta = x - target
        projected = dot_product(coupling, delta)
        product = matmul(matrix, (1.0_dp + 0.6_dp * delta) * vec + &
            0.06_dp * coupling * projected * dot_product(coupling, vec))
        ok = all(ieee_is_finite(product))
    end subroutine jvp
end program benchmark_lgmres_image_cache
