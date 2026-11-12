! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Randomized interpolative decomposition of dense real/complex matrices.
! SciPy analogue: scipy.linalg.interpolative.interp_decomp(..., rand=True).
!
! Sketch A on the left with a small Gaussian matrix; column-pivoted QR of
! the sketch selects original columns. Fit interpolation coefficients on the
! ORIGINAL, scaled matrix by minimum-norm least squares, not on the sketch.
! Fixed-rank calls use k+oversampling sketch rows (capped at min(m,n)).
! Adaptive calls grow the sketch and CERTIFY the relative Frobenius residual
! on A. This criterion differs from SciPy's spectral-norm target.
!
! Return a one-based permutation idx of length n and a k-by-(n-k) proj:
! A(:,idx(k+1:n)) ~= matmul(A(:,idx(1:k)),proj).
module scifort_linalg_interpolative_random
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr_factor
    use scifort_linalg_interpolative, only : interp_decomp
    use scifort_linalg_svd, only : svd_lstsq
    implicit none
    private

    public :: randomized_interp_decomp

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64
    integer(int64), parameter :: default_seed = 20261110_int64
    real(dp), parameter :: pi = acos(-1.0_dp)

    interface randomized_interp_decomp
        module procedure random_id_real_fixed
        module procedure random_id_complex_fixed
        module procedure random_id_real_tol
        module procedure random_id_complex_tol
    end interface randomized_interp_decomp

contains

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! input matrix
        ok = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    pure integer(int64) function initial_state(seed) result(state)
        integer(int64), intent(in), optional :: seed !! reproducible RNG seed
        state = modulo(default_seed, rng_modulus-1_int64) + 1_int64
        if (present(seed)) state = modulo(seed, rng_modulus-1_int64) + 1_int64
    end function initial_state

    real(dp) function draw_normal(state) result(z)
        integer(int64), intent(inout) :: state !! local RNG state
        real(dp) :: u1, u2
        state = modulo(state*rng_multiplier, rng_modulus)
        u1 = real(state,dp)/real(rng_modulus,dp)
        state = modulo(state*rng_multiplier, rng_modulus)
        u2 = real(state,dp)/real(rng_modulus,dp)
        z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
    end function draw_normal

    subroutine sketch_real(a, rows, state, y)
        real(dp), intent(in) :: a(:, :) !! scaled m-by-n source
        integer, intent(in) :: rows !! number of random sketch rows
        integer(int64), intent(inout) :: state !! local RNG state
        real(dp), allocatable, intent(out) :: y(:, :) !! rows-by-n sketch
        real(dp), allocatable :: g(:, :)
        integer :: i, j
        allocate(g(rows,size(a,1)))
        do j = 1, size(g,2)
            do i = 1, rows
                g(i,j) = draw_normal(state)/sqrt(real(rows,dp))
            end do
        end do
        y = matmul(g,a)
    end subroutine sketch_real

    subroutine sketch_complex(a, rows, state, y)
        complex(dp), intent(in) :: a(:, :) !! scaled m-by-n source
        integer, intent(in) :: rows !! number of random sketch rows
        integer(int64), intent(inout) :: state !! local RNG state
        complex(dp), allocatable, intent(out) :: y(:, :) !! rows-by-n sketch
        complex(dp), allocatable :: g(:, :)
        integer :: i, j
        real(dp) :: normalization
        normalization = sqrt(2.0_dp*real(rows,dp))
        allocate(g(rows,size(a,1)))
        do j = 1, size(g,2)
            do i = 1, rows
                g(i,j) = cmplx(draw_normal(state),draw_normal(state),dp)/normalization
            end do
        end do
        y = matmul(g,a)
    end subroutine sketch_complex

    subroutine choose_real(a, rows, state, idx, status)
        real(dp), intent(in) :: a(:, :) !! scaled source matrix
        integer, intent(in) :: rows !! number of sketch rows
        integer(int64), intent(inout) :: state !! local RNG state
        integer, allocatable, intent(out) :: idx(:) !! one-based column ordering
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), allocatable :: y(:, :), af(:, :), tau(:)
        call sketch_real(a,rows,state,y)
        allocate(af(rows,size(a,2)),tau(rows),idx(size(a,2)))
        call qr_factor(y,af,tau,status,idx)
    end subroutine choose_real

    subroutine choose_complex(a, rows, state, idx, status)
        complex(dp), intent(in) :: a(:, :) !! scaled source matrix
        integer, intent(in) :: rows !! number of sketch rows
        integer(int64), intent(inout) :: state !! local RNG state
        integer, allocatable, intent(out) :: idx(:) !! one-based column ordering
        integer, intent(out) :: status !! linalg_status_* result
        complex(dp), allocatable :: y(:, :), af(:, :), tau(:)
        call sketch_complex(a,rows,state,y)
        allocate(af(rows,size(a,2)),tau(rows),idx(size(a,2)))
        call qr_factor(y,af,tau,status,idx)
    end subroutine choose_complex

    subroutine fit_real(a, k, idx, proj, status)
        real(dp), intent(in) :: a(:, :) !! scaled ORIGINAL matrix
        integer, intent(in) :: k !! number of selected columns
        integer, intent(in) :: idx(:) !! one-based column ordering
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) weights
        integer, intent(out) :: status !! linalg_status_* result
        integer :: n
        n = size(a,2)
        allocate(proj(k,n-k))
        proj = 0.0_dp
        status = linalg_status_success
        if (k == 0 .or. k == n) return
        call svd_lstsq(a(:,idx(1:k)),a(:,idx(k+1:n)),proj,status)
    end subroutine fit_real

    subroutine fit_complex(a, k, idx, proj, status)
        complex(dp), intent(in) :: a(:, :) !! scaled ORIGINAL matrix
        integer, intent(in) :: k !! number of selected columns
        integer, intent(in) :: idx(:) !! one-based column ordering
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) weights
        integer, intent(out) :: status !! linalg_status_* result
        integer :: n
        n = size(a,2)
        allocate(proj(k,n-k))
        proj = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_success
        if (k == 0 .or. k == n) return
        call svd_lstsq(a(:,idx(1:k)),a(:,idx(k+1:n)),proj,status)
    end subroutine fit_complex

    pure real(dp) function residual_real(a, k, idx, proj) result(err)
        real(dp), intent(in) :: a(:, :) !! scaled ORIGINAL matrix
        integer, intent(in) :: k !! selected rank
        integer, intent(in) :: idx(:) !! selected column permutation
        real(dp), intent(in) :: proj(:, :) !! fitted weights
        real(dp) :: denominator
        denominator = sqrt(sum(a*a))
        err = 0.0_dp
        if (denominator == 0.0_dp) return
        if (k == 0) then
            err = 1.0_dp
        else if (k < size(a,2)) then
            err = sqrt(sum((a(:,idx(k+1:))-matmul(a(:,idx(1:k)),proj))**2))/denominator
        end if
    end function residual_real

    pure real(dp) function residual_complex(a, k, idx, proj) result(err)
        complex(dp), intent(in) :: a(:, :) !! scaled ORIGINAL matrix
        integer, intent(in) :: k !! selected rank
        integer, intent(in) :: idx(:) !! selected column permutation
        complex(dp), intent(in) :: proj(:, :) !! fitted weights
        real(dp) :: denominator
        denominator = sqrt(sum(abs(a)**2))
        err = 0.0_dp
        if (denominator == 0.0_dp) return
        if (k == 0) then
            err = 1.0_dp
        else if (k < size(a,2)) then
            err = sqrt(sum(abs(a(:,idx(k+1:))-matmul(a(:,idx(1:k)),proj))**2))/denominator
        end if
    end function residual_complex

    subroutine random_id_real_fixed(a, k, idx, proj, status, seed, oversampling)
        real(dp), intent(in) :: a(:, :) !! m-by-n input matrix
        integer, intent(in) :: k !! requested rank 0..min(m,n)
        integer, allocatable, intent(out) :: idx(:) !! one-based n-column permutation
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) weights
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! extra sketch rows, default 8
        real(dp), allocatable :: w(:, :)
        integer(int64) :: state
        integer :: m, n, over, rows, i
        real(dp) :: scale
        allocate(idx(0),proj(0,0))
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        over = 8
        if (present(oversampling)) over = oversampling
        if (k < 0 .or. k > min(m,n) .or. over < 0) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        scale = 0.0_dp
        if (m > 0 .and. n > 0) scale = maxval(abs(a))
        if (k == 0 .or. scale == 0.0_dp) then
            deallocate(idx,proj)
            allocate(idx(n),proj(k,n-k))
            idx = [(i,i=1,n)]
            proj = 0.0_dp
            status = linalg_status_success
            return
        end if
        w = a/scale
        state = initial_state(seed)
        rows = k + min(over,min(m,n)-k)
        deallocate(idx)
        call choose_real(w,rows,state,idx,status)
        if (status == linalg_status_success) call fit_real(w,k,idx,proj,status)
        if (status /= linalg_status_success) then
            if (allocated(idx)) deallocate(idx)
            if (allocated(proj)) deallocate(proj)
            allocate(idx(0),proj(0,0))
        end if
    end subroutine random_id_real_fixed

    subroutine random_id_complex_fixed(a, k, idx, proj, status, seed, oversampling)
        complex(dp), intent(in) :: a(:, :) !! m-by-n input matrix
        integer, intent(in) :: k !! requested rank 0..min(m,n)
        integer, allocatable, intent(out) :: idx(:) !! one-based n-column permutation
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) weights
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! extra sketch rows, default 8
        complex(dp), allocatable :: w(:, :)
        integer(int64) :: state
        integer :: m, n, over, rows, i
        real(dp) :: scale
        allocate(idx(0),proj(0,0))
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        over = 8
        if (present(oversampling)) over = oversampling
        if (k < 0 .or. k > min(m,n) .or. over < 0) return
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        if (k == 0) then
            deallocate(idx,proj)
            allocate(idx(n),proj(0,n))
            idx = [(i,i=1,n)]
            status = linalg_status_success
            return
        end if
        scale = max(maxval(abs(real(a,dp))),maxval(abs(aimag(a))))
        if (scale == 0.0_dp) then
            deallocate(idx,proj)
            allocate(idx(n),proj(k,n-k))
            idx = [(i,i=1,n)]
            proj = cmplx(0.0_dp,0.0_dp,dp)
            status = linalg_status_success
            return
        end if
        w = a/scale
        state = initial_state(seed)
        rows = k + min(over,min(m,n)-k)
        deallocate(idx)
        call choose_complex(w,rows,state,idx,status)
        if (status == linalg_status_success) call fit_complex(w,k,idx,proj,status)
        if (status /= linalg_status_success) then
            if (allocated(idx)) deallocate(idx)
            if (allocated(proj)) deallocate(proj)
            allocate(idx(0),proj(0,0))
        end if
    end subroutine random_id_complex_fixed

    subroutine random_id_real_tol(a, eps, k, idx, proj, status, seed, oversampling)
        real(dp), intent(in) :: a(:, :) !! m-by-n input matrix
        real(dp), intent(in) :: eps !! relative Frobenius tolerance, 0 <= eps < 1
        integer, intent(out) :: k !! selected numerical rank
        integer, allocatable, intent(out) :: idx(:) !! one-based n-column permutation
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) weights
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! extra sketch rows, default 8
        real(dp), allocatable :: w(:, :), y(:, :), proposed(:, :)
        integer, allocatable :: candidate_idx(:)
        integer(int64) :: state
        integer :: m, n, rmax, over, rows, candidate, i
        real(dp) :: scale, residual, floor_tolerance
        k = 0
        allocate(idx(0),proj(0,0))
        status = linalg_status_invalid_shape
        over = 8
        if (present(oversampling)) over = oversampling
        if (.not. ieee_is_finite(eps) .or. eps < 0.0_dp .or. eps >= 1.0_dp .or. over < 0) return
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
                rows = min(rmax,max(8,2*min(max(1,over),max(1,rmax/2))))
                floor_tolerance = 100.0_dp*epsilon(1.0_dp)*real(max(m,n),dp)
                do
                    call sketch_real(w,rows,state,y)
                    call interp_decomp(y,eps,candidate,candidate_idx,proposed,status)
                    if (status /= linalg_status_success) exit
                    call fit_real(w,candidate,candidate_idx,proj,status)
                    if (status /= linalg_status_success) exit
                    residual = residual_real(w,candidate,candidate_idx,proj)
                    if (ieee_is_finite(residual)) then
                        if (residual <= eps*(1.0_dp+1.0e-12_dp)+floor_tolerance) exit
                    end if
                    if (rows == rmax) then
                        ! A stable last resort uses pivoted QR of the original.
                        call interp_decomp(w,rmax,candidate_idx,proj,status)
                        candidate = rmax
                        if (status == linalg_status_success) then
                            residual = residual_real(w,candidate,candidate_idx,proj)
                            if (.not. ieee_is_finite(residual)) status = linalg_status_no_convergence
                            if (residual > eps*(1.0_dp+1.0e-12_dp)+floor_tolerance) &
                                status = linalg_status_no_convergence
                        end if
                        exit
                    end if
                    rows = rows + min(rows,rmax-rows)
                end do
                if (status == linalg_status_success) then
                    k = candidate
                    call move_alloc(candidate_idx,idx)
                end if
            end if
        end if
        if (status == linalg_status_success .and. size(idx) == 0) then
            deallocate(idx,proj)
            allocate(idx(n),proj(0,n))
            idx = [(i,i=1,n)]
        end if
        if (status /= linalg_status_success) then
            k = 0
            if (allocated(idx)) deallocate(idx)
            if (allocated(proj)) deallocate(proj)
            allocate(idx(0),proj(0,0))
        end if
    end subroutine random_id_real_tol

    subroutine random_id_complex_tol(a, eps, k, idx, proj, status, seed, oversampling)
        complex(dp), intent(in) :: a(:, :) !! m-by-n input matrix
        real(dp), intent(in) :: eps !! relative Frobenius tolerance, 0 <= eps < 1
        integer, intent(out) :: k !! selected numerical rank
        integer, allocatable, intent(out) :: idx(:) !! one-based n-column permutation
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) weights
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! extra sketch rows, default 8
        complex(dp), allocatable :: w(:, :), y(:, :), proposed(:, :)
        integer, allocatable :: candidate_idx(:)
        integer(int64) :: state
        integer :: m, n, rmax, over, rows, candidate, i
        real(dp) :: scale, residual, floor_tolerance
        k = 0
        allocate(idx(0),proj(0,0))
        status = linalg_status_invalid_shape
        over = 8
        if (present(oversampling)) over = oversampling
        if (.not. ieee_is_finite(eps) .or. eps < 0.0_dp .or. eps >= 1.0_dp .or. over < 0) return
        m = size(a,1)
        n = size(a,2)
        status = linalg_status_no_convergence
        if (.not. finite_complex(a)) return
        rmax = min(m,n)
        if (rmax == 0) then
            status = linalg_status_success
        else
            scale = max(maxval(abs(real(a,dp))),maxval(abs(aimag(a))))
            if (scale == 0.0_dp) then
                status = linalg_status_success
            else
                w = a/scale
                state = initial_state(seed)
                rows = min(rmax,max(8,2*min(max(1,over),max(1,rmax/2))))
                floor_tolerance = 100.0_dp*epsilon(1.0_dp)*real(max(m,n),dp)
                do
                    call sketch_complex(w,rows,state,y)
                    call interp_decomp(y,eps,candidate,candidate_idx,proposed,status)
                    if (status /= linalg_status_success) exit
                    call fit_complex(w,candidate,candidate_idx,proj,status)
                    if (status /= linalg_status_success) exit
                    residual = residual_complex(w,candidate,candidate_idx,proj)
                    if (ieee_is_finite(residual)) then
                        if (residual <= eps*(1.0_dp+1.0e-12_dp)+floor_tolerance) exit
                    end if
                    if (rows == rmax) then
                        call interp_decomp(w,rmax,candidate_idx,proj,status)
                        candidate = rmax
                        if (status == linalg_status_success) then
                            residual = residual_complex(w,candidate,candidate_idx,proj)
                            if (.not. ieee_is_finite(residual)) status = linalg_status_no_convergence
                            if (residual > eps*(1.0_dp+1.0e-12_dp)+floor_tolerance) &
                                status = linalg_status_no_convergence
                        end if
                        exit
                    end if
                    rows = rows + min(rows,rmax-rows)
                end do
                if (status == linalg_status_success) then
                    k = candidate
                    call move_alloc(candidate_idx,idx)
                end if
            end if
        end if
        if (status == linalg_status_success .and. size(idx) == 0) then
            deallocate(idx,proj)
            allocate(idx(n),proj(0,n))
            idx = [(i,i=1,n)]
        end if
        if (status /= linalg_status_success) then
            k = 0
            if (allocated(idx)) deallocate(idx)
            if (allocated(proj)) deallocate(proj)
            allocate(idx(0),proj(0,0))
        end if
    end subroutine random_id_complex_tol

end module scifort_linalg_interpolative_random
