! SPDX-License-Identifier: MIT
program test_linalg_interpolative_operator
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_interpolative_operator, only : randomized_interp_decomp_operator, &
        randomized_interp_decomp_operator_complex
    use scifort_linalg_interpolative, only : reconstruct_matrix_from_id
    implicit none
    real(dp), allocatable :: a(:, :), p(:, :), skel(:, :), a_rec(:, :), old_p(:, :)
    complex(dp), allocatable :: z(:, :), q(:, :), skel_z(:, :), z_rec(:, :)
    integer, allocatable :: idx(:), old_idx(:), idxz(:)
    integer :: m, n, j, i, variant, mode, k, kz, status, checks, num_fwd, num_adj
    real(dp) :: eps, err, true_error, scale
    logical :: fail_fwd, fail_adj, emit_nan
    integer(int64), parameter :: seed = 811831_int64

    checks = 0
    fail_fwd = .false.
    fail_adj = .false.
    emit_nan = .false.
    do m = 0, 10
        do n = 0, 12
            allocate(a(m,n),z(m,n),a_rec(m,n),z_rec(m,n))
            do variant = 1, 5
                do j = 1, n
                    do i = 1, m
                        a(i,j) = sin(real(3*i+2*j+variant,dp)) + cos(real(i-3*j+variant,dp))
                        z(i,j) = cmplx(a(i,j),sin(real(7*i+5*j+variant,dp)),dp)
                    end do
                end do
                select case(variant)
                case(2)
                    a = 0.0_dp
                    z = cmplx(0.0_dp,0.0_dp,dp)
                case(3)
                    if (n > 1) then
                        do j = 2,n
                            a(:,j) = a(:,1)*real(j,dp)
                            z(:,j) = z(:,1)*cmplx(real(j,dp),real(j,dp)*0.1_dp,dp)
                        end do
                    end if
                case(4)
                    a = a*1.0e-80_dp
                    z = z*1.0e-80_dp
                case(5)
                    a = a*1.0e80_dp
                    z = z*1.0e80_dp
                end select
                do mode = 1, 6
                    eps = 0.0_dp
                    select case(mode)
                    case(1)
                        k = 0
                    case(2)
                        k = min(m,n)/2
                    case(3)
                        k = min(m,n)
                    case(4)
                        eps = 0.25_dp
                    case(5)
                        eps = 1.0e-5_dp
                    case(6)
                        eps = 0.0_dp
                    end select
                    num_fwd = 0
                    num_adj = 0
                    if (mode <= 3) then
                        call randomized_interp_decomp_operator(m,n,k,real_forward,real_adjoint, &
                            idx,p,status,seed=seed,oversampling=3,skeleton=skel)
                    else
                        call randomized_interp_decomp_operator(m,n,eps,real_forward,real_adjoint, &
                            k,idx,p,status,seed=seed,oversampling=3,verification_probes=8, &
                            skeleton=skel,estimated_error=err)
                    end if
                    if (status /= linalg_status_success) then
                        print *, 'real failed',m,n,variant,mode,k,status
                        error stop 'real operator ID status'
                    end if
                    call check_perm(idx,n)
                    if (size(skel,1) /= m .or. size(skel,2) /= k) error stop 'real skeleton shape'
                    if (size(p,1) /= k .or. size(p,2) /= n-k) error stop 'real coefficients shape'
                    if (k > 0) then
                        if (any(skel /= a(:,idx(1:k)))) error stop 'real skeleton not exact'
                    end if
                    if (mode <= 3) then
                        if (num_fwd /= k) error stop 'fixed real used unwanted forward products'
                    end if
                    call reconstruct_matrix_from_id(skel,idx,p,a_rec,status)
                    if (status /= linalg_status_success) error stop 'real reconstruction API'
                    scale = 0.0_dp
                    if (m*n > 0) scale = maxval(abs(a))
                    true_error = 0.0_dp
                    if (scale > 0.0_dp) then
                        true_error = sqrt(sum(((a-a_rec)/scale)**2))/sqrt(sum((a/scale)**2))
                    end if
                    if (mode == 3 .or. mode >= 4) then
                        if (mode == 3 .and. true_error > 1.0e-8_dp) error stop 'real full-rank accuracy'
                        if (mode >= 4 .and. true_error > max(1.0e-8_dp,eps*1.03_dp)) then
                            print *, 'real accuracy',m,n,variant,mode,k,true_error,eps
                            error stop 'real adaptive accuracy'
                        end if
                    end if
                    num_fwd = 0
                    num_adj = 0
                    if (mode <= 3) then
                        call randomized_interp_decomp_operator_complex(m,n,k,complex_forward,complex_adjoint, &
                            idxz,q,status,seed=seed,oversampling=3,skeleton=skel_z)
                    else
                        call randomized_interp_decomp_operator_complex(m,n,eps,complex_forward,complex_adjoint, &
                            kz,idxz,q,status,seed=seed,oversampling=3,verification_probes=8, &
                            skeleton=skel_z,estimated_error=err)
                        k = kz
                    end if
                    if (status /= linalg_status_success) then
                        print *, 'complex failed',m,n,variant,mode,k,status
                        error stop 'complex operator ID status'
                    end if
                    call check_perm(idxz,n)
                    if (size(skel_z,1) /= m .or. size(skel_z,2) /= k) error stop 'complex skeleton shape'
                    if (size(q,1) /= k .or. size(q,2) /= n-k) error stop 'complex coefficients shape'
                    if (k > 0) then
                        if (any(skel_z /= z(:,idxz(1:k)))) error stop 'complex skeleton not exact'
                    end if
                    if (mode <= 3) then
                        if (num_fwd /= k) error stop 'fixed complex used unwanted forward products'
                    end if
                    call reconstruct_matrix_from_id(skel_z,idxz,q,z_rec,status)
                    if (status /= linalg_status_success) error stop 'complex reconstruction API'
                    scale = 0.0_dp
                    if (m*n > 0) scale = max(maxval(abs(real(z,dp))),maxval(abs(aimag(z))))
                    true_error = 0.0_dp
                    if (scale > 0.0_dp) then
                        true_error = sqrt(sum(abs((z-z_rec)/scale)**2))/sqrt(sum(abs(z/scale)**2))
                    end if
                    if (mode == 3 .or. mode >= 4) then
                        if (mode == 3 .and. true_error > 1.0e-8_dp) error stop 'complex full-rank accuracy'
                        if (mode >= 4 .and. true_error > max(1.0e-8_dp,eps*1.03_dp)) then
                            print *, 'complex accuracy',m,n,variant,mode,k,true_error,eps
                            error stop 'complex adaptive accuracy'
                        end if
                    end if
                    checks = checks + 10
                end do
            end do
            deallocate(a,z,a_rec,z_rec)
        end do
    end do

    allocate(a(20,50),z(20,50))
    do j = 1,50
        do i = 1,20
            a(i,j) = sin(real(i+j,dp)) + cos(real(i*j,dp))
            z(i,j) = cmplx(a(i,j),sin(real(2*i+3*j,dp)),dp)
        end do
    end do
    num_fwd = 0
    num_adj = 0
    call randomized_interp_decomp_operator(20,50,3,real_forward,real_adjoint, &
        idx,p,status,seed=seed,oversampling=4)
    if (status /= linalg_status_success) error stop 'large fixed rank'
    if (num_fwd /= 3 .or. num_adj /= 10) error stop 'matrix-free callback count'
    old_idx = idx
    old_p = p
    call randomized_interp_decomp_operator(20,50,3,real_forward,real_adjoint, &
        idx,p,status,seed=seed,oversampling=4)
    if (status /= linalg_status_success) error stop 'repeat status'
    if (any(idx /= old_idx) .or. any(p /= old_p)) error stop 'repeatability'
    if (size(idx) /= 50) error stop 'index size'

    ! Validate empty dimensions, bad parameters, and failures of either callback.
    call randomized_interp_decomp_operator(-1,5,0,real_forward,real_adjoint,idx,p,status)
    if (status /= linalg_status_invalid_shape) error stop 'negative dimensions'
    call randomized_interp_decomp_operator(20,50,21,real_forward,real_adjoint,idx,p,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid rank'
    call randomized_interp_decomp_operator(20,50,2,real_forward,real_adjoint, &
        idx,p,status,oversampling=-1)
    if (status /= linalg_status_invalid_shape) error stop 'invalid oversampling'
    call randomized_interp_decomp_operator(20,50,0.2_dp,real_forward,real_adjoint, &
        k,idx,p,status,verification_probes=0)
    if (status /= linalg_status_invalid_shape) error stop 'invalid probe count'
    call randomized_interp_decomp_operator(20,50,-0.1_dp,real_forward,real_adjoint,k,idx,p,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid tolerance'
    call randomized_interp_decomp_operator(20,50,1.0_dp,real_forward,real_adjoint,k,idx,p,status)
    if (status /= linalg_status_invalid_shape) error stop 'tolerance >=1'
    fail_adj = .true.
    call randomized_interp_decomp_operator(20,50,3,real_forward,real_adjoint,idx,p,status)
    if (status /= linalg_status_no_convergence) error stop 'adjoint error propagation'
    fail_adj = .false.
    fail_fwd = .true.
    call randomized_interp_decomp_operator(20,50,3,real_forward,real_adjoint,idx,p,status)
    if (status /= linalg_status_no_convergence) error stop 'forward error propagation'
    fail_fwd = .false.
    emit_nan = .true.
    call randomized_interp_decomp_operator(20,50,3,real_forward,real_adjoint,idx,p,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite products'
    emit_nan = .false.
    call randomized_interp_decomp_operator_complex(20,50,3,complex_forward,complex_adjoint,idxz,q,status)
    if (status /= linalg_status_success) error stop 'complex recovery'
    checks = checks + 11
    print '(a,i0)', 'matrix-free interpolative ID native checks passed: ',checks
contains
    subroutine check_perm(perm,ncols)
        integer, intent(in) :: perm(:) !! one-based permutation
        integer, intent(in) :: ncols !! expected number of columns
        logical :: seen(ncols)
        integer :: t
        if (size(perm) /= ncols) error stop 'permutation size'
        seen = .false.
        do t = 1,ncols
            if (perm(t) < 1 .or. perm(t) > ncols) error stop 'permutation range'
            if (seen(perm(t))) error stop 'duplicate permutation index'
            seen(perm(t)) = .true.
        end do
    end subroutine check_perm

    subroutine real_forward(x,y,istat)
        real(dp), intent(in) :: x(:) !! vector input
        real(dp), intent(out) :: y(:) !! product
        integer, intent(out) :: istat !! success/failure
        num_fwd = num_fwd + 1
        istat = 9
        if (fail_fwd) return
        y = matmul(a,x)
        if (emit_nan .and. size(y) > 0) y(1) = ieee_value(0.0_dp,ieee_quiet_nan)
        istat = 0
    end subroutine real_forward
    subroutine real_adjoint(x,y,istat)
        real(dp), intent(in) :: x(:) !! vector input
        real(dp), intent(out) :: y(:) !! adjoint product
        integer, intent(out) :: istat !! success/failure
        num_adj = num_adj + 1
        istat = 9
        if (fail_adj) return
        y = matmul(transpose(a),x)
        if (emit_nan .and. size(y) > 0) y(1) = ieee_value(0.0_dp,ieee_quiet_nan)
        istat = 0
    end subroutine real_adjoint
    subroutine complex_forward(x,y,istat)
        complex(dp), intent(in) :: x(:) !! vector input
        complex(dp), intent(out) :: y(:) !! product
        integer, intent(out) :: istat !! success/failure
        num_fwd = num_fwd + 1
        istat = 9
        if (fail_fwd) return
        y = matmul(z,x)
        istat = 0
    end subroutine complex_forward
    subroutine complex_adjoint(x,y,istat)
        complex(dp), intent(in) :: x(:) !! vector input
        complex(dp), intent(out) :: y(:) !! adjoint product
        integer, intent(out) :: istat !! success/failure
        num_adj = num_adj + 1
        istat = 9
        if (fail_adj) return
        y = matmul(transpose(conjg(z)),x)
        istat = 0
    end subroutine complex_adjoint
end program test_linalg_interpolative_operator
