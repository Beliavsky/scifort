! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Randomized, ID-based truncated SVD for dense real and complex matrices.
! Analogous to scipy.linalg.interpolative.svd(A, eps_or_k, rand=True).
!
! A small Gaussian row sketch G*A provides pivot candidates.  Unlike an
! approximation based only on the sketch, the interpolation coefficients are
! computed by a least-squares fit to the ORIGINAL matrix.  Existing ID-to-SVD
! conversion yields orthonormal U,V and descending singular values.
!
! The adaptive interface explicitly verifies the relative Frobenius norm of
! the residual on A; it grows the sketch until the requested tolerance is
! achieved.  This is not SciPy's spectral-norm stopping rule.  Neither the
! stochastic stream nor individual pivot choices need agree with NumPy.
module scifort_linalg_interpolative_svd_random
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr_factor
    use scifort_linalg_interpolative, only : interp_decomp, id_to_svd
    use scifort_linalg_svd, only : svd_lstsq
    implicit none
    private

    public :: randomized_interpolative_svd

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64
    integer(int64), parameter :: default_seed = 20261110_int64
    real(dp), parameter :: pi = acos(-1.0_dp)

    interface randomized_interpolative_svd
        module procedure random_svd_real_rank
        module procedure random_svd_complex_rank
        module procedure random_svd_real_tol
        module procedure random_svd_complex_tol
    end interface randomized_interpolative_svd

contains

    pure logical function complex_finite(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! candidate input matrix
        ok = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function complex_finite

    pure integer(int64) function initial_state(seed) result(state)
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        state = modulo(default_seed, rng_modulus - 1_int64) + 1_int64
        if (present(seed)) state = modulo(seed, rng_modulus - 1_int64) + 1_int64
    end function initial_state

    real(dp) function draw_uniform(state) result(x)
        integer(int64), intent(inout) :: state !! Park-Miller RNG state
        state = modulo(state * rng_multiplier, rng_modulus)
        x = real(state, dp) / real(rng_modulus, dp)
    end function draw_uniform

    real(dp) function draw_normal(state) result(x)
        integer(int64), intent(inout) :: state !! Park-Miller RNG state
        real(dp) :: u1, u2
        u1 = draw_uniform(state)
        u2 = draw_uniform(state)
        x = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
    end function draw_normal

    subroutine real_sketch(a, l, state, y)
        real(dp), intent(in) :: a(:, :) !! scaled input m-by-n
        integer, intent(in) :: l !! sketch row count, 1..min(m,n)
        integer(int64), intent(inout) :: state !! local random state
        real(dp), allocatable, intent(out) :: y(:, :) !! l-by-n row sketch
        real(dp), allocatable :: g(:, :)
        integer :: i, j
        allocate(g(l, size(a,1)))
        do j = 1, size(g,2)
            do i = 1, l
                g(i,j) = draw_normal(state) / sqrt(real(l,dp))
            end do
        end do
        y = matmul(g, a)
    end subroutine real_sketch

    subroutine complex_sketch(a, l, state, y)
        complex(dp), intent(in) :: a(:, :) !! scaled complex input
        integer, intent(in) :: l !! sketch row count
        integer(int64), intent(inout) :: state !! local random state
        complex(dp), allocatable, intent(out) :: y(:, :) !! l-by-n sketch
        complex(dp), allocatable :: g(:, :)
        real(dp) :: scale
        integer :: i, j
        scale = sqrt(2.0_dp * real(l,dp))
        allocate(g(l, size(a,1)))
        do j = 1, size(g,2)
            do i = 1, l
                g(i,j) = cmplx(draw_normal(state), draw_normal(state), dp) / scale
            end do
        end do
        y = matmul(g, a)
    end subroutine complex_sketch

    ! Compute a true least-squares interpolative fit, using the original
    ! matrix and the pivot order proposed by the small random sketch.
    subroutine original_id_real(a, k, idx, proj, status)
        real(dp), intent(in) :: a(:, :) !! original scaled matrix
        integer, intent(in) :: k !! skeleton size
        integer, intent(in) :: idx(:) !! one-based column permutation
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) fit
        integer, intent(out) :: status !! success or failure
        integer :: n
        n = size(a,2)
        allocate(proj(k,n-k))
        proj = 0.0_dp
        status = linalg_status_success
        if (k == 0 .or. k == n) return
        call svd_lstsq(a(:,idx(1:k)), a(:,idx(k+1:n)), proj, status)
    end subroutine original_id_real

    subroutine original_id_complex(a, k, idx, proj, status)
        complex(dp), intent(in) :: a(:, :) !! original scaled matrix
        integer, intent(in) :: k !! skeleton size
        integer, intent(in) :: idx(:) !! one-based column permutation
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) fit
        integer, intent(out) :: status !! success or failure
        integer :: n
        n = size(a,2)
        allocate(proj(k,n-k))
        proj = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_success
        if (k == 0 .or. k == n) return
        call svd_lstsq(a(:,idx(1:k)), a(:,idx(k+1:n)), proj, status)
    end subroutine original_id_complex

    pure real(dp) function id_error_real(a, k, idx, proj) result(err)
        real(dp), intent(in) :: a(:, :) !! scaled original matrix
        integer, intent(in) :: k !! selected rank
        integer, intent(in) :: idx(:) !! pivot ordering
        real(dp), intent(in) :: proj(:, :) !! original-matrix coefficients
        real(dp) :: denom
        denom = sqrt(sum(a*a))
        err = 0.0_dp
        if (denom == 0.0_dp) return
        if (k == 0) then
            err = 1.0_dp
        else if (k < size(a,2)) then
            err = sqrt(sum((a(:,idx(k+1:)) - matmul(a(:,idx(1:k)),proj))**2)) / denom
        end if
    end function id_error_real

    pure real(dp) function id_error_complex(a, k, idx, proj) result(err)
        complex(dp), intent(in) :: a(:, :) !! scaled original matrix
        integer, intent(in) :: k !! selected rank
        integer, intent(in) :: idx(:) !! pivot ordering
        complex(dp), intent(in) :: proj(:, :) !! original-matrix coefficients
        real(dp) :: denom
        denom = sqrt(sum(abs(a)**2))
        err = 0.0_dp
        if (denom == 0.0_dp) return
        if (k == 0) then
            err = 1.0_dp
        else if (k < size(a,2)) then
            err = sqrt(sum(abs(a(:,idx(k+1:)) - matmul(a(:,idx(1:k)),proj))**2)) / denom
        end if
    end function id_error_complex

    ! QR of the row sketch ONLY, not QR of all m rows of A.  For a fixed
    ! rank k, the sketch has at most k+oversampling rows.
    subroutine choose_real(a, l, state, idx, status)
        real(dp), intent(in) :: a(:, :) !! scaled matrix
        integer, intent(in) :: l !! row sketch size
        integer(int64), intent(inout) :: state !! seeded state
        integer, allocatable, intent(out) :: idx(:) !! pivot ordering
        integer, intent(out) :: status !! QR status
        real(dp), allocatable :: y(:, :), af(:, :), tau(:)
        call real_sketch(a, l, state, y)
        allocate(af(l,size(a,2)),tau(l),idx(size(a,2)))
        call qr_factor(y,af,tau,status,idx)
    end subroutine choose_real

    subroutine choose_complex(a, l, state, idx, status)
        complex(dp), intent(in) :: a(:, :) !! scaled matrix
        integer, intent(in) :: l !! row sketch size
        integer(int64), intent(inout) :: state !! seeded state
        integer, allocatable, intent(out) :: idx(:) !! pivot ordering
        integer, intent(out) :: status !! QR status
        complex(dp), allocatable :: y(:, :), af(:, :), tau(:)
        call complex_sketch(a, l, state, y)
        allocate(af(l,size(a,2)),tau(l),idx(size(a,2)))
        call qr_factor(y,af,tau,status,idx)
    end subroutine choose_complex

    subroutine real_factors(a, k, idx, proj, u, s, v, status)
        real(dp), intent(in) :: a(:, :) !! scaled matrix
        integer, intent(in) :: k !! selected rank
        integer, intent(in) :: idx(:) !! pivot permutation
        real(dp), intent(in) :: proj(:, :) !! interpolation coefficients
        real(dp), intent(out) :: u(:, :) !! left singular vectors
        real(dp), intent(out) :: s(:) !! scaled singular values
        real(dp), intent(out) :: v(:, :) !! right singular vectors
        integer, intent(out) :: status !! factorization status
        call id_to_svd(a(:,idx(1:k)),idx,proj,u,s,v,status)
    end subroutine real_factors

    subroutine complex_factors(a, k, idx, proj, u, s, v, status)
        complex(dp), intent(in) :: a(:, :) !! scaled matrix
        integer, intent(in) :: k !! selected rank
        integer, intent(in) :: idx(:) !! pivot permutation
        complex(dp), intent(in) :: proj(:, :) !! interpolation coefficients
        complex(dp), intent(out) :: u(:, :) !! left singular vectors
        real(dp), intent(out) :: s(:) !! scaled singular values
        complex(dp), intent(out) :: v(:, :) !! right singular vectors
        integer, intent(out) :: status !! factorization status
        call id_to_svd(a(:,idx(1:k)),idx,proj,u,s,v,status)
    end subroutine complex_factors

    subroutine random_svd_real_rank(a, k, u, s, v, status, seed, oversampling)
        real(dp), intent(in) :: a(:, :) !! m-by-n input
        integer, intent(in) :: k !! target rank, 0..min(m,n)
        real(dp), intent(out) :: u(:, :) !! m-by-k orthonormal columns
        real(dp), intent(out) :: s(:) !! k singular values, descending
        real(dp), intent(out) :: v(:, :) !! n-by-k orthonormal columns
        integer, intent(out) :: status !! linalg_status_* code
        integer(int64), optional, intent(in) :: seed !! reproducible local RNG seed
        integer, optional, intent(in) :: oversampling !! extra sketch rows, default 8
        real(dp), allocatable :: w(:, :), proj(:, :)
        integer, allocatable :: idx(:)
        integer(int64) :: state
        integer :: m, n, over, l
        real(dp) :: scale
        u = 0.0_dp
        s = 0.0_dp
        v = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        over = 8
        if (present(oversampling)) over = oversampling
        if (over < 0 .or. k < 0 .or. k > min(m,n)) return
        if (size(u,1) /= m .or. size(u,2) /= k) return
        if (size(s) /= k .or. size(v,1) /= n .or. size(v,2) /= k) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        status = linalg_status_success
        if (k == 0) return
        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            ! Obtain orthonormal zero-singular-value factors from the old path.
            call interp_decomp(a,k,idx,proj,status)
            if (status == linalg_status_success) call real_factors(a,k,idx,proj,u,s,v,status)
            return
        end if
        w = a/scale
        state = initial_state(seed)
        l = min(min(m,n),k + min(over,max(0,min(m,n)-k)))
        call choose_real(w,l,state,idx,status)
        if (status /= linalg_status_success) return
        call original_id_real(w,k,idx,proj,status)
        if (status /= linalg_status_success) return
        call real_factors(w,k,idx,proj,u,s,v,status)
        if (status /= linalg_status_success) return
        s = scale*s
        if (.not. all(ieee_is_finite(s))) status = linalg_status_no_convergence
        if (status /= linalg_status_success) then
            u = 0.0_dp
            s = 0.0_dp
            v = 0.0_dp
        end if
    end subroutine random_svd_real_rank

    subroutine random_svd_complex_rank(a, k, u, s, v, status, seed, oversampling)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex input
        integer, intent(in) :: k !! target rank
        complex(dp), intent(out) :: u(:, :) !! m-by-k orthonormal columns
        real(dp), intent(out) :: s(:) !! k singular values, descending
        complex(dp), intent(out) :: v(:, :) !! n-by-k orthonormal columns
        integer, intent(out) :: status !! linalg_status_* code
        integer(int64), optional, intent(in) :: seed !! reproducible local RNG seed
        integer, optional, intent(in) :: oversampling !! extra sketch rows, default 8
        complex(dp), allocatable :: w(:, :), proj(:, :)
        integer, allocatable :: idx(:)
        integer(int64) :: state
        integer :: m, n, over, l
        real(dp) :: scale
        u = cmplx(0.0_dp,0.0_dp,dp)
        s = 0.0_dp
        v = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        over = 8
        if (present(oversampling)) over = oversampling
        if (over < 0 .or. k < 0 .or. k > min(m,n)) return
        if (size(u,1) /= m .or. size(u,2) /= k) return
        if (size(s) /= k .or. size(v,1) /= n .or. size(v,2) /= k) return
        status = linalg_status_no_convergence
        if (.not. complex_finite(a)) return
        status = linalg_status_success
        if (k == 0) return
        scale = max(maxval(abs(real(a,dp))), maxval(abs(aimag(a))))
        if (scale == 0.0_dp) then
            call interp_decomp(a,k,idx,proj,status)
            if (status == linalg_status_success) call complex_factors(a,k,idx,proj,u,s,v,status)
            return
        end if
        w = a/scale
        state = initial_state(seed)
        l = min(min(m,n),k + min(over,max(0,min(m,n)-k)))
        call choose_complex(w,l,state,idx,status)
        if (status /= linalg_status_success) return
        call original_id_complex(w,k,idx,proj,status)
        if (status /= linalg_status_success) return
        call complex_factors(w,k,idx,proj,u,s,v,status)
        if (status /= linalg_status_success) return
        s = scale*s
        if (.not. all(ieee_is_finite(s))) status = linalg_status_no_convergence
        if (status /= linalg_status_success) then
            u = cmplx(0.0_dp,0.0_dp,dp)
            s = 0.0_dp
            v = cmplx(0.0_dp,0.0_dp,dp)
        end if
    end subroutine random_svd_complex_rank

    subroutine random_svd_real_tol(a, eps, k, u, s, v, status, seed, oversampling)
        real(dp), intent(in) :: a(:, :) !! m-by-n matrix
        real(dp), intent(in) :: eps !! relative Frobenius tolerance, 0 <= eps < 1
        integer, intent(out) :: k !! selected numerical rank
        real(dp), allocatable, intent(out) :: u(:, :) !! allocated m-by-k vectors
        real(dp), allocatable, intent(out) :: s(:) !! allocated k singular values
        real(dp), allocatable, intent(out) :: v(:, :) !! allocated n-by-k vectors
        integer, intent(out) :: status !! linalg_status_* code
        integer(int64), optional, intent(in) :: seed !! reproducible local seed
        integer, optional, intent(in) :: oversampling !! extra sketch rows, default 8
        real(dp), allocatable :: w(:, :), y(:, :), proj(:, :), proposed(:, :)
        integer, allocatable :: idx(:)
        integer(int64) :: state
        integer :: m, n, over, l, rmax, candidate
        real(dp) :: scale, residual
        k = 0
        allocate(u(0,0),s(0),v(0,0))
        status = linalg_status_invalid_shape
        over = 8
        if (present(oversampling)) over = oversampling
        if (over < 0 .or. .not. ieee_is_finite(eps)) return
        if (eps < 0.0_dp .or. eps >= 1.0_dp) return
        m = size(a,1)
        n = size(a,2)
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        rmax = min(m,n)
        if (rmax == 0) then
            status = linalg_status_success
        else
            scale = maxval(abs(a))
            if (scale == 0.0_dp) then
                status = linalg_status_success
            else
                w = a/scale
                state = initial_state(seed)
                l = min(rmax,max(8,2*min(max(1,over),max(1,rmax/2))))
                do
                    call real_sketch(w,l,state,y)
                    call interp_decomp(y,eps,candidate,idx,proposed,status)
                    if (status /= linalg_status_success) exit
                    call original_id_real(w,candidate,idx,proj,status)
                    if (status /= linalg_status_success) exit
                    residual = id_error_real(w,candidate,idx,proj)
                    if (residual <= eps*(1.0_dp + 1.0e-12_dp)) exit
                    if (l == rmax) then
                        ! Full-rank candidate always gives a best-available fit.
                        candidate = rmax
                        deallocate(proj)
                        call original_id_real(w,candidate,idx,proj,status)
                        exit
                    end if
                    l = l + min(l,rmax-l)
                    deallocate(y)
                end do
                if (status == linalg_status_success) then
                    k = candidate
                    deallocate(u,s,v)
                    allocate(u(m,k),s(k),v(n,k))
                    call real_factors(w,k,idx,proj,u,s,v,status)
                    if (status == linalg_status_success) then
                        s = scale*s
                        if (.not. all(ieee_is_finite(s))) status = linalg_status_no_convergence
                    end if
                end if
            end if
        end if
        if (status == linalg_status_success .and. rmax == 0) then
            deallocate(u,s,v)
            allocate(u(m,0),s(0),v(n,0))
        else if (status == linalg_status_success .and. k == 0) then
            deallocate(u,s,v)
            allocate(u(m,0),s(0),v(n,0))
        else if (status /= linalg_status_success) then
            k = 0
            deallocate(u,s,v)
            allocate(u(0,0),s(0),v(0,0))
        end if
    end subroutine random_svd_real_tol

    subroutine random_svd_complex_tol(a, eps, k, u, s, v, status, seed, oversampling)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex matrix
        real(dp), intent(in) :: eps !! relative Frobenius tolerance, 0 <= eps < 1
        integer, intent(out) :: k !! selected numerical rank
        complex(dp), allocatable, intent(out) :: u(:, :) !! allocated m-by-k vectors
        real(dp), allocatable, intent(out) :: s(:) !! allocated k singular values
        complex(dp), allocatable, intent(out) :: v(:, :) !! allocated n-by-k vectors
        integer, intent(out) :: status !! linalg_status_* code
        integer(int64), optional, intent(in) :: seed !! reproducible local seed
        integer, optional, intent(in) :: oversampling !! extra sketch rows, default 8
        complex(dp), allocatable :: w(:, :), y(:, :), proj(:, :), proposed(:, :)
        integer, allocatable :: idx(:)
        integer(int64) :: state
        integer :: m, n, over, l, rmax, candidate
        real(dp) :: scale, residual
        k = 0
        allocate(u(0,0),s(0),v(0,0))
        status = linalg_status_invalid_shape
        over = 8
        if (present(oversampling)) over = oversampling
        if (over < 0 .or. .not. ieee_is_finite(eps)) return
        if (eps < 0.0_dp .or. eps >= 1.0_dp) return
        m = size(a,1)
        n = size(a,2)
        status = linalg_status_no_convergence
        if (.not. complex_finite(a)) return
        rmax = min(m,n)
        if (rmax == 0) then
            status = linalg_status_success
        else
            scale = max(maxval(abs(real(a,dp))), maxval(abs(aimag(a))))
            if (scale == 0.0_dp) then
                status = linalg_status_success
            else
                w = a/scale
                state = initial_state(seed)
                l = min(rmax,max(8,2*min(max(1,over),max(1,rmax/2))))
                do
                    call complex_sketch(w,l,state,y)
                    call interp_decomp(y,eps,candidate,idx,proposed,status)
                    if (status /= linalg_status_success) exit
                    call original_id_complex(w,candidate,idx,proj,status)
                    if (status /= linalg_status_success) exit
                    residual = id_error_complex(w,candidate,idx,proj)
                    if (residual <= eps*(1.0_dp + 1.0e-12_dp)) exit
                    if (l == rmax) then
                        candidate = rmax
                        deallocate(proj)
                        call original_id_complex(w,candidate,idx,proj,status)
                        exit
                    end if
                    l = l + min(l,rmax-l)
                    deallocate(y)
                end do
                if (status == linalg_status_success) then
                    k = candidate
                    deallocate(u,s,v)
                    allocate(u(m,k),s(k),v(n,k))
                    call complex_factors(w,k,idx,proj,u,s,v,status)
                    if (status == linalg_status_success) then
                        s = scale*s
                        if (.not. all(ieee_is_finite(s))) status = linalg_status_no_convergence
                    end if
                end if
            end if
        end if
        if (status == linalg_status_success .and. rmax == 0) then
            deallocate(u,s,v)
            allocate(u(m,0),s(0),v(n,0))
        else if (status == linalg_status_success .and. k == 0) then
            deallocate(u,s,v)
            allocate(u(m,0),s(0),v(n,0))
        else if (status /= linalg_status_success) then
            k = 0
            deallocate(u,s,v)
            allocate(u(0,0),s(0),v(0,0))
        end if
    end subroutine random_svd_complex_tol

end module scifort_linalg_interpolative_svd_random
