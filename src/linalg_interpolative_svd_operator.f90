! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Randomized interpolative SVD for matrix-free linear operators.
! Analogue of scipy.linalg.interpolative.svd(LinearOperator, eps_or_k).
!
! A is accessed ONLY through matvec (A*x) and rmatvec (A^H*x).
! Randomized ID provides selected columns B and interpolation coefficients;
! a QR and a small-core SVD convert B*P into U*diag(s)*V^H.
!
! The tolerance interface uses an independent random-probe estimate of the
! relative Frobenius error, NOT a certified spectral-norm error bound.
module scifort_linalg_interpolative_svd_operator
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_interpolative_operator, only : randomized_interp_decomp_operator, &
        randomized_interp_decomp_operator_complex, id_real_matvec, id_complex_matvec
    use scifort_linalg_interpolative, only : id_to_svd
    implicit none
    private

    public :: randomized_interpolative_svd_operator
    public :: randomized_interpolative_svd_operator_complex

    ! The generic accepts real and complex matrices. The _complex synonym
    ! mirrors the naming of scifort_linalg_interpolative_operator.
    interface randomized_interpolative_svd_operator
        module procedure real_fixed, real_tolerance, complex_fixed, complex_tolerance
    end interface randomized_interpolative_svd_operator
    interface randomized_interpolative_svd_operator_complex
        module procedure complex_fixed, complex_tolerance
    end interface randomized_interpolative_svd_operator_complex

contains

    subroutine convert_real(b, idx, proj, u, s, v, status)
        real(dp), intent(in) :: b(:, :) !! m-by-k exact selected columns of A
        integer, intent(in) :: idx(:) !! one-based permutation of 1:n
        real(dp), intent(in) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        real(dp), intent(out) :: u(:, :) !! m-by-k left singular vectors
        real(dp), intent(out) :: s(:) !! k singular values in descending order
        real(dp), intent(out) :: v(:, :) !! n-by-k right singular vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp) :: scale

        u = 0.0_dp
        s = 0.0_dp
        v = 0.0_dp
        status = linalg_status_success
        if (size(s) == 0) return
        scale = maxval(abs(b))
        if (.not. ieee_is_finite(scale)) then
            status = linalg_status_no_convergence
            return
        end if
        if (scale <= 0.0_dp) then
            call id_to_svd(b,idx,proj,u,s,v,status)
        else
            call id_to_svd(b/scale,idx,proj,u,s,v,status)
            if (status == linalg_status_success) s = s*scale
        end if
        if (status == linalg_status_success) then
            if (.not. all(ieee_is_finite(s)) .or. .not. all(ieee_is_finite(u)) .or. &
                .not. all(ieee_is_finite(v))) status = linalg_status_no_convergence
        end if
        if (status /= linalg_status_success) then
            u = 0.0_dp
            s = 0.0_dp
            v = 0.0_dp
        end if
    end subroutine convert_real

    subroutine convert_complex(b, idx, proj, u, s, v, status)
        complex(dp), intent(in) :: b(:, :) !! m-by-k exact selected columns
        integer, intent(in) :: idx(:) !! one-based permutation of 1:n
        complex(dp), intent(in) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        complex(dp), intent(out) :: u(:, :) !! m-by-k left singular vectors
        real(dp), intent(out) :: s(:) !! k nonnegative singular values, descending
        complex(dp), intent(out) :: v(:, :) !! n-by-k right singular vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp) :: scale

        u = cmplx(0.0_dp,0.0_dp,dp)
        s = 0.0_dp
        v = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_success
        if (size(s) == 0) return
        scale = maxval(abs(b))
        if (.not. ieee_is_finite(scale)) then
            status = linalg_status_no_convergence
            return
        end if
        if (scale <= 0.0_dp) then
            call id_to_svd(b,idx,proj,u,s,v,status)
        else
            call id_to_svd(b/scale,idx,proj,u,s,v,status)
            if (status == linalg_status_success) s = s*scale
        end if
        if (status == linalg_status_success) then
            if (.not. all(ieee_is_finite(s)) .or. &
                .not. all(ieee_is_finite(real(u,dp))) .or. .not. all(ieee_is_finite(aimag(u))) .or. &
                .not. all(ieee_is_finite(real(v,dp))) .or. .not. all(ieee_is_finite(aimag(v)))) &
                status = linalg_status_no_convergence
        end if
        if (status /= linalg_status_success) then
            u = cmplx(0.0_dp,0.0_dp,dp)
            s = 0.0_dp
            v = cmplx(0.0_dp,0.0_dp,dp)
        end if
    end subroutine convert_complex

    subroutine real_fixed(m, n, k, matvec, rmatvec, u, s, v, status, seed, oversampling)
        integer, intent(in) :: m !! number of operator rows
        integer, intent(in) :: n !! number of operator columns
        integer, intent(in) :: k !! target rank, 0..min(m,n)
        procedure(id_real_matvec) :: matvec !! computes A*x, length n to m
        procedure(id_real_matvec) :: rmatvec !! computes A^T*x, length m to n
        real(dp), intent(out) :: u(:, :) !! caller-allocated (m,k) left vectors
        real(dp), intent(out) :: s(:) !! caller-allocated (k) singular values
        real(dp), intent(out) :: v(:, :) !! caller-allocated (n,k) right vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! local reproducible random seed
        integer, intent(in), optional :: oversampling !! sketch oversampling, default 8
        integer, allocatable :: idx(:)
        real(dp), allocatable :: proj(:, :), b(:, :)

        u = 0.0_dp
        s = 0.0_dp
        v = 0.0_dp
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0 .or. k < 0 .or. k > min(m,n)) return
        if (size(u,1) /= m .or. size(u,2) /= k .or. size(s) /= k) return
        if (size(v,1) /= n .or. size(v,2) /= k) return
        call randomized_interp_decomp_operator(m,n,k,matvec,rmatvec,idx,proj,status, &
            seed=seed,oversampling=oversampling,skeleton=b)
        if (status /= linalg_status_success .or. k == 0) return
        call convert_real(b,idx,proj,u,s,v,status)
    end subroutine real_fixed

    subroutine complex_fixed(m, n, k, matvec, rmatvec, u, s, v, status, seed, oversampling)
        integer, intent(in) :: m !! number of operator rows
        integer, intent(in) :: n !! number of operator columns
        integer, intent(in) :: k !! target rank, 0..min(m,n)
        procedure(id_complex_matvec) :: matvec !! computes A*x, length n to m
        procedure(id_complex_matvec) :: rmatvec !! computes A^H*x, length m to n
        complex(dp), intent(out) :: u(:, :) !! caller-allocated (m,k) left vectors
        real(dp), intent(out) :: s(:) !! caller-allocated (k) singular values
        complex(dp), intent(out) :: v(:, :) !! caller-allocated (n,k) right vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! sketch oversampling, default 8
        integer, allocatable :: idx(:)
        complex(dp), allocatable :: proj(:, :), b(:, :)

        u = cmplx(0.0_dp,0.0_dp,dp)
        s = 0.0_dp
        v = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (m < 0 .or. n < 0 .or. k < 0 .or. k > min(m,n)) return
        if (size(u,1) /= m .or. size(u,2) /= k .or. size(s) /= k) return
        if (size(v,1) /= n .or. size(v,2) /= k) return
        call randomized_interp_decomp_operator_complex(m,n,k,matvec,rmatvec,idx,proj,status, &
            seed=seed,oversampling=oversampling,skeleton=b)
        if (status /= linalg_status_success .or. k == 0) return
        call convert_complex(b,idx,proj,u,s,v,status)
    end subroutine complex_fixed

    subroutine real_tolerance(m, n, eps, matvec, rmatvec, k, u, s, v, status, &
        seed, oversampling, verification_probes, estimated_error)
        integer, intent(in) :: m !! number of operator rows
        integer, intent(in) :: n !! number of operator columns
        real(dp), intent(in) :: eps !! relative Frobenius tolerance, [0,1)
        procedure(id_real_matvec) :: matvec !! computes A*x
        procedure(id_real_matvec) :: rmatvec !! computes A^T*x
        integer, intent(out) :: k !! selected rank
        real(dp), allocatable, intent(out) :: u(:, :) !! allocated (m,k) left vectors
        real(dp), allocatable, intent(out) :: s(:) !! allocated (k) singular values
        real(dp), allocatable, intent(out) :: v(:, :) !! allocated (n,k) right vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! local reproducible random seed
        integer, intent(in), optional :: oversampling !! sketch oversampling, default 8
        integer, intent(in), optional :: verification_probes !! residual probes, default 16
        real(dp), intent(out), optional :: estimated_error !! sampled relative Frobenius residual
        integer, allocatable :: idx(:)
        real(dp), allocatable :: proj(:, :), b(:, :)

        k = 0
        allocate(u(0,0),s(0),v(0,0))
        call randomized_interp_decomp_operator(m,n,eps,matvec,rmatvec,k,idx,proj,status, &
            seed=seed,oversampling=oversampling,verification_probes=verification_probes, &
            skeleton=b,estimated_error=estimated_error)
        if (status /= linalg_status_success) then
            k = 0
            return
        end if
        deallocate(u,s,v)
        allocate(u(m,k),s(k),v(n,k))
        if (k == 0) then
            u = 0.0_dp
            s = 0.0_dp
            v = 0.0_dp
            return
        end if
        call convert_real(b,idx,proj,u,s,v,status)
        if (status /= linalg_status_success) then
            k = 0
            deallocate(u,s,v)
            allocate(u(0,0),s(0),v(0,0))
        end if
    end subroutine real_tolerance

    subroutine complex_tolerance(m, n, eps, matvec, rmatvec, k, u, s, v, status, &
        seed, oversampling, verification_probes, estimated_error)
        integer, intent(in) :: m !! number of operator rows
        integer, intent(in) :: n !! number of operator columns
        real(dp), intent(in) :: eps !! relative Frobenius tolerance, [0,1)
        procedure(id_complex_matvec) :: matvec !! computes A*x
        procedure(id_complex_matvec) :: rmatvec !! computes A^H*x
        integer, intent(out) :: k !! selected rank
        complex(dp), allocatable, intent(out) :: u(:, :) !! allocated (m,k) left vectors
        real(dp), allocatable, intent(out) :: s(:) !! allocated (k) singular values
        complex(dp), allocatable, intent(out) :: v(:, :) !! allocated (n,k) right vectors
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! local reproducible random seed
        integer, intent(in), optional :: oversampling !! sketch oversampling, default 8
        integer, intent(in), optional :: verification_probes !! residual probes, default 16
        real(dp), intent(out), optional :: estimated_error !! sampled relative Frobenius residual
        integer, allocatable :: idx(:)
        complex(dp), allocatable :: proj(:, :), b(:, :)

        k = 0
        allocate(u(0,0),s(0),v(0,0))
        call randomized_interp_decomp_operator_complex(m,n,eps,matvec,rmatvec,k,idx,proj,status, &
            seed=seed,oversampling=oversampling,verification_probes=verification_probes, &
            skeleton=b,estimated_error=estimated_error)
        if (status /= linalg_status_success) then
            k = 0
            return
        end if
        deallocate(u,s,v)
        allocate(u(m,k),s(k),v(n,k))
        if (k == 0) then
            u = cmplx(0.0_dp,0.0_dp,dp)
            s = 0.0_dp
            v = cmplx(0.0_dp,0.0_dp,dp)
            return
        end if
        call convert_complex(b,idx,proj,u,s,v,status)
        if (status /= linalg_status_success) then
            k = 0
            deallocate(u,s,v)
            allocate(u(0,0),s(0),v(0,0))
        end if
    end subroutine complex_tolerance

end module scifort_linalg_interpolative_svd_operator
