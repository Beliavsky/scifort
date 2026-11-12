! SPDX-License-Identifier: MIT
program test_linalg_interpolative_random
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_interpolative_random, only : randomized_interp_decomp
    use scifort_linalg_interpolative, only : reconstruct_matrix_from_id
    implicit none
    real(dp), allocatable :: a(:, :), ar(:, :), proj(:, :), last_proj(:, :)
    complex(dp), allocatable :: z(:, :), zr(:, :), projz(:, :)
    integer, allocatable :: idx(:), idxz(:), last_idx(:)
    integer :: m, n, variant, mode, rank, k, status, i, j, checks
    real(dp) :: eps, error, scale, xnan
    integer(int64), parameter :: seed = 318771_int64
    real(dp), parameter :: roundoff_tolerance = 1.0e-8_dp

    checks = 0
    do m = 0, 11
        do n = 0, 12
            allocate(a(m,n),z(m,n),ar(m,n),zr(m,n))
            do variant = 1, 5
                do j = 1, n
                    do i = 1, m
                        a(i,j) = sin(real(3*i+7*j+variant,dp)) + 0.17_dp*cos(real(11*i-4*j,dp))
                        z(i,j) = cmplx(a(i,j),cos(real(5*i+13*j+variant,dp)),dp)
                    end do
                end do
                select case (variant)
                case (2)
                    a = 0.0_dp
                    z = cmplx(0.0_dp,0.0_dp,dp)
                case (3)
                    if (n > 1) then
                        do j = 2,n
                            a(:,j) = a(:,1)*real(j,dp)
                            z(:,j) = z(:,1)*cmplx(real(j,dp),0.3_dp*real(j,dp),dp)
                        end do
                    end if
                case (4)
                    a = a*1.0e-80_dp
                    z = z*1.0e-80_dp
                case (5)
                    a = a*1.0e80_dp
                    z = z*1.0e80_dp
                end select
                do mode = 1, 6
                    select case (mode)
                    case (1)
                        k = 0
                    case (2)
                        k = min(m,n)/2
                    case (3)
                        k = min(m,n)
                    case (4)
                        eps = 0.2_dp
                    case (5)
                        eps = 1.0e-5_dp
                    case (6)
                        eps = 0.0_dp
                    end select
                    if (mode <= 3) then
                        call randomized_interp_decomp(a,k,idx,proj,status,seed=seed)
                    else
                        call randomized_interp_decomp(a,eps,k,idx,proj,status,seed=seed)
                    end if
                    if (status /= linalg_status_success) then
                        print *, 'real ID status',m,n,variant,mode,status
                        error stop 'real randomized ID'
                    end if
                    call verify_perm(idx,n)
                    if (k < 0 .or. k > min(m,n)) error stop 'invalid real rank'
                    if (size(proj,1) /= k .or. size(proj,2) /= n-k) error stop 'real proj shape'
                    if (m*n > 0) then
                        scale = maxval(abs(a))
                    else
                        scale = 0.0_dp
                    end if
                    if (scale /= 0.0_dp) then
                        call reconstruct_matrix_from_id(a(:,idx(1:k)),idx,proj,ar,status)
                        if (status /= linalg_status_success) error stop 'real ID reconstruction API'
                        error = sqrt(sum(((a-ar)/scale)**2))/sqrt(sum((a/scale)**2))
                        if (mode == 3 .or. mode >= 4) then
                            if (mode == 3) then
                                if (error > roundoff_tolerance) error stop 'real full rank residual'
                            else
                                if (error > max(roundoff_tolerance,eps*1.00001_dp)) then
                                    print *, 'real adaptive residual',m,n,variant,mode,k,error,eps
                                    error stop 'real adaptive tolerance'
                                end if
                            end if
                        end if
                    end if
                    if (mode <= 3) then
                        call randomized_interp_decomp(z,k,idxz,projz,status,seed=seed)
                    else
                        call randomized_interp_decomp(z,eps,rank,idxz,projz,status,seed=seed)
                        k = rank
                    end if
                    if (status /= linalg_status_success) then
                        print *, 'complex ID status',m,n,variant,mode,status
                        error stop 'complex randomized ID'
                    end if
                    call verify_perm(idxz,n)
                    if (k < 0 .or. k > min(m,n)) error stop 'invalid complex rank'
                    if (size(projz,1) /= k .or. size(projz,2) /= n-k) error stop 'complex proj shape'
                    if (m*n > 0) then
                        scale = max(maxval(abs(real(z,dp))),maxval(abs(aimag(z))))
                    else
                        scale = 0.0_dp
                    end if
                    if (scale /= 0.0_dp) then
                        call reconstruct_matrix_from_id(z(:,idxz(1:k)),idxz,projz,zr,status)
                        if (status /= linalg_status_success) error stop 'complex ID reconstruction API'
                        error = sqrt(sum(abs((z-zr)/scale)**2))/sqrt(sum(abs(z/scale)**2))
                        if (mode == 3 .or. mode >= 4) then
                            if (mode == 3) then
                                if (error > roundoff_tolerance) error stop 'complex full rank residual'
                            else
                                if (error > max(roundoff_tolerance,eps*1.00001_dp)) then
                                    print *, 'complex adaptive residual',m,n,variant,mode,k,error,eps
                                    error stop 'complex adaptive tolerance'
                                end if
                            end if
                        end if
                    end if
                    checks = checks+8
                end do
            end do
            deallocate(a,z,ar,zr)
        end do
    end do

    allocate(a(5,4))
    a = 1.0_dp
    call randomized_interp_decomp(a,5,idx,proj,status,seed=seed)
    if (status /= linalg_status_invalid_shape) error stop 'bad rank'
    call randomized_interp_decomp(a,-1,idx,proj,status,seed=seed)
    if (status /= linalg_status_invalid_shape) error stop 'negative rank'
    call randomized_interp_decomp(a,1,idx,proj,status,oversampling=-1)
    if (status /= linalg_status_invalid_shape) error stop 'negative oversampling'
    call randomized_interp_decomp(a,1.0_dp,k,idx,proj,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad eps'
    call randomized_interp_decomp(a,-0.1_dp,k,idx,proj,status)
    if (status /= linalg_status_invalid_shape) error stop 'negative eps'
    xnan = ieee_value(0.0_dp,ieee_quiet_nan)
    a(1,1) = xnan
    call randomized_interp_decomp(a,2,idx,proj,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite fixed'
    call randomized_interp_decomp(a,0.1_dp,k,idx,proj,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite adaptive'
    a = 1.0_dp
    call randomized_interp_decomp(a,2,idx,proj,status,seed=seed)
    if (status /= linalg_status_success) error stop 'repeatability setup'
    last_idx = idx
    last_proj = proj
    call randomized_interp_decomp(a,2,idx,proj,status,seed=seed)
    if (any(idx /= last_idx)) error stop 'seed repeatability pivot'
    if (any(proj /= last_proj)) error stop 'seed repeatability coefficients'
    checks = checks + 9
    print '(a,i0)', 'randomized interpolative ID checks passed: ',checks
contains
    subroutine verify_perm(p, ncols)
        integer, intent(in) :: p(:) !! one-based column indices
        integer, intent(in) :: ncols !! expected number of columns
        logical :: visited(ncols)
        integer :: t
        if (size(p) /= ncols) error stop 'permutation length'
        visited = .false.
        do t=1,ncols
            if (p(t) < 1 .or. p(t) > ncols) error stop 'permutation range'
            if (visited(p(t))) error stop 'duplicate pivot'
            visited(p(t)) = .true.
        end do
    end subroutine verify_perm
end program test_linalg_interpolative_random
