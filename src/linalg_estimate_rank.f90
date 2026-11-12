! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Adaptive randomized numerical-rank estimator for dense matrices and
! matrix-free operators. Analogous to scipy.linalg.interpolative.estimate_rank.
!
! Each random test vector generates a candidate image A*x. Two-pass modified
! Gram-Schmidt finds its component outside the accumulated column space.
! A candidate is retained when the residual norm exceeds eps times the
! largest observed image norm. After `oversampling` consecutive rejections,
! the rank estimate is complete. This is a randomized estimate, not a
! certified spectral-norm rank. It need not equal SciPy's numerical output.
module scifort_linalg_estimate_rank
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    implicit none
    private

    public :: estimate_rank
    public :: estimate_rank_operator, estimate_rank_operator_complex
    public :: rank_real_matvec, rank_complex_matvec

    abstract interface
        subroutine rank_real_matvec(x, y, status)
            import dp
            real(dp), intent(in) :: x(:) !! input vector, size n
            real(dp), intent(out) :: y(:) !! output vector, size m
            integer, intent(out) :: status !! zero for success
        end subroutine rank_real_matvec
        subroutine rank_complex_matvec(x, y, status)
            import dp
            complex(dp), intent(in) :: x(:) !! input vector, size n
            complex(dp), intent(out) :: y(:) !! output vector, size m
            integer, intent(out) :: status !! zero for success
        end subroutine rank_complex_matvec
    end interface

    interface estimate_rank
        module procedure rank_real_dense
        module procedure rank_complex_dense
    end interface estimate_rank

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64

contains

    real(dp) function rng_draw(state) result(x)
        integer(int64), intent(inout) :: state !! local state, 1..modulus-1
        state = modulo(state * rng_multiplier, rng_modulus)
        x = 2.0_dp * real(state,dp) / real(rng_modulus,dp) - 1.0_dp
    end function rng_draw

    subroutine random_real(x, state, seeded)
        real(dp), intent(out) :: x(:) !! random probe vector
        integer(int64), intent(inout) :: state !! deterministic RNG state
        logical, intent(in) :: seeded !! true to use local state
        integer :: j
        real(dp) :: u
        do j = 1, size(x)
            if (seeded) then
                x(j) = rng_draw(state)
            else
                call random_number(u)
                x(j) = 2.0_dp*u - 1.0_dp
            end if
        end do
    end subroutine random_real

    subroutine random_complex(x, state, seeded)
        complex(dp), intent(out) :: x(:) !! random complex probe
        integer(int64), intent(inout) :: state !! deterministic RNG state
        logical, intent(in) :: seeded !! true to use local state
        integer :: j
        real(dp) :: u, v
        do j = 1, size(x)
            if (seeded) then
                u = rng_draw(state)
                v = rng_draw(state)
            else
                call random_number(u)
                call random_number(v)
                u = 2.0_dp*u - 1.0_dp
                v = 2.0_dp*v - 1.0_dp
            end if
            x(j) = cmplx(u,v,dp)
        end do
    end subroutine random_complex

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! matrix to validate
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    pure logical function finite_complex_vector(a) result(ok)
        complex(dp), intent(in) :: a(:) !! vector to validate
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex_vector

    pure real(dp) function norm_real(x) result(r)
        real(dp), intent(in) :: x(:) !! real vector
        real(dp) :: s
        r = 0.0_dp
        if (size(x) == 0) return
        s = maxval(abs(x))
        if (s == 0.0_dp) return
        r = s * sqrt(sum((x/s)**2))
    end function norm_real

    pure real(dp) function norm_complex(x) result(r)
        complex(dp), intent(in) :: x(:) !! complex vector
        real(dp) :: s
        r = 0.0_dp
        if (size(x) == 0) return
        s = maxval(abs(x))
        if (s == 0.0_dp) return
        r = s * sqrt(sum(abs(x/s)**2))
    end function norm_complex

    subroutine rank_real_dense(a, eps, rank, status, seed, oversampling)
        real(dp), intent(in) :: a(:, :) !! dense real m-by-n input
        real(dp), intent(in) :: eps !! relative residual tolerance, 0 <= eps < 1
        integer, intent(out) :: rank !! estimated rank, 0..min(m,n)
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! consecutive small probes, default 8
        real(dp), allocatable :: scaled(:, :)
        real(dp) :: scale
        rank = 0
        status = linalg_status_invalid_shape
        if (.not. ieee_is_finite(eps)) return
        if (eps < 0.0_dp .or. eps >= 1.0_dp) return
        if (present(oversampling)) then
            if (oversampling < 1) return
        end if
        if (.not. all(ieee_is_finite(a))) then
            status = linalg_status_no_convergence
            return
        end if
        if (size(a,1) == 0 .or. size(a,2) == 0) then
            status = linalg_status_success
            return
        end if
        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            status = linalg_status_success
            return
        end if
        scaled = a/scale
        call rank_real_core(size(a,1),size(a,2),eps,rank,status,seed,oversampling,a=scaled)
    end subroutine rank_real_dense

    subroutine rank_complex_dense(a, eps, rank, status, seed, oversampling)
        complex(dp), intent(in) :: a(:, :) !! dense complex m-by-n input
        real(dp), intent(in) :: eps !! relative residual tolerance, 0 <= eps < 1
        integer, intent(out) :: rank !! estimated rank, 0..min(m,n)
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! consecutive small probes, default 8
        complex(dp), allocatable :: scaled(:, :)
        real(dp) :: scale
        rank = 0
        status = linalg_status_invalid_shape
        if (.not. ieee_is_finite(eps)) return
        if (eps < 0.0_dp .or. eps >= 1.0_dp) return
        if (present(oversampling)) then
            if (oversampling < 1) return
        end if
        if (.not. finite_complex(a)) then
            status = linalg_status_no_convergence
            return
        end if
        if (size(a,1) == 0 .or. size(a,2) == 0) then
            status = linalg_status_success
            return
        end if
        ! Magnitudes of finite complex numbers can overflow when both
        ! components approach HUGE(dp); scale by the components instead.
        scale = max(maxval(abs(real(a,dp))),maxval(abs(aimag(a))))
        if (scale == 0.0_dp) then
            status = linalg_status_success
            return
        end if
        scaled = a/scale
        call rank_complex_core(size(a,1),size(a,2),eps,rank,status,seed,oversampling,a=scaled)
    end subroutine rank_complex_dense

    subroutine estimate_rank_operator(m, n, mv, eps, rank, status, seed, oversampling)
        integer, intent(in) :: m !! number of rows, >= 0
        integer, intent(in) :: n !! number of columns, >= 0
        procedure(rank_real_matvec) :: mv !! callback implementing y = A*x
        real(dp), intent(in) :: eps !! relative residual tolerance, 0 <= eps < 1
        integer, intent(out) :: rank !! estimated rank
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        integer, intent(in), optional :: oversampling !! consecutive small probes, default 8
        call rank_real_core(m,n,eps,rank,status,seed,oversampling,mv=mv)
    end subroutine estimate_rank_operator

    subroutine estimate_rank_operator_complex(m, n, mv, eps, rank, status, seed, oversampling)
        integer, intent(in) :: m !! number of rows, >= 0
        integer, intent(in) :: n !! number of columns, >= 0
        procedure(rank_complex_matvec) :: mv !! callback implementing y = A*x
        real(dp), intent(in) :: eps !! relative residual tolerance, 0 <= eps < 1
        integer, intent(out) :: rank !! estimated rank
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        integer, intent(in), optional :: oversampling !! consecutive small probes, default 8
        call rank_complex_core(m,n,eps,rank,status,seed,oversampling,mv=mv)
    end subroutine estimate_rank_operator_complex

    subroutine rank_real_core(m, n, eps, rank, status, seed, oversampling, a, mv)
        integer, intent(in) :: m !! number of output rows
        integer, intent(in) :: n !! number of input columns
        real(dp), intent(in) :: eps !! relative tolerance, 0 <= eps < 1
        integer, intent(out) :: rank !! rank estimate
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! deterministic seed
        integer, intent(in), optional :: oversampling !! consecutive rejected probes
        real(dp), intent(in), optional :: a(:, :) !! scaled dense matrix
        procedure(rank_real_matvec), optional :: mv !! matrix-free forward action
        real(dp), allocatable :: basis(:, :), x(:), y(:), work(:)
        real(dp) :: ref, r, alpha
        integer(int64) :: state
        integer :: tries, extra, maxrank, misses, i, pass, callback_status
        logical :: seeded

        rank = 0
        status = linalg_status_invalid_shape
        extra = 8
        if (present(oversampling)) extra = oversampling
        if (m < 0 .or. n < 0 .or. extra < 1) return
        if (.not. ieee_is_finite(eps)) return
        if (eps < 0.0_dp .or. eps >= 1.0_dp) return
        if (present(a)) then
            if (size(a,1) /= m .or. size(a,2) /= n) return
        else
            if (.not. present(mv)) return
        end if
        status = linalg_status_success
        if (m == 0 .or. n == 0) return
        maxrank = min(m,n)
        allocate(basis(m,maxrank),x(n),y(m),work(m))
        state = 1_int64
        seeded = present(seed)
        if (seeded) state = modulo(seed, rng_modulus-1_int64) + 1_int64
        ref = 0.0_dp
        misses = 0
        ! No need for more than one test per potential basis direction and
        ! `extra` final rejected probes. Guard arithmetic on loop bounds.
        do tries = 1, maxrank
            ! Probe attempts are managed by the inner routine so we can
            ! continue past rejected samples without reducing the budget.
            call take_probe(mv)
            if (status /= linalg_status_success) return
            if (rank == maxrank) return
            if (misses >= extra) return
        end do
        do tries = 1, extra
            call take_probe(mv)
            if (status /= linalg_status_success) return
            if (rank == maxrank .or. misses >= extra) return
        end do
    contains
        subroutine take_probe(mv)
            procedure(rank_real_matvec), optional :: mv !! matrix-free forward action
            call random_real(x,state,seeded)
            r = norm_real(x)
            if (r == 0.0_dp) then
                misses = misses + 1
                return
            end if
            x = x/r
            if (present(a)) then
                y = matmul(a,x)
            else
                call mv(x,y,callback_status)
                if (callback_status /= 0) then
                    status = linalg_status_no_convergence
                    rank = 0
                    return
                end if
            end if
            if (.not. all(ieee_is_finite(y))) then
                status = linalg_status_no_convergence
                rank = 0
                return
            end if
            r = norm_real(y)
            if (.not. ieee_is_finite(r)) then
                status = linalg_status_no_convergence
                rank = 0
                return
            end if
            ref = max(ref,r)
            work = y
            do pass = 1, 2
                do i = 1, rank
                    alpha = dot_product(basis(:,i),work)
                    work = work - alpha*basis(:,i)
                end do
            end do
            r = norm_real(work)
            if (.not. ieee_is_finite(r)) then
                status = linalg_status_no_convergence
                rank = 0
                return
            end if
            if (r > eps*ref .and. r > 0.0_dp) then
                rank = rank + 1
                basis(:,rank) = work/r
                misses = 0
            else
                misses = misses + 1
            end if
        end subroutine take_probe
    end subroutine rank_real_core

    subroutine rank_complex_core(m, n, eps, rank, status, seed, oversampling, a, mv)
        integer, intent(in) :: m !! number of output rows
        integer, intent(in) :: n !! number of input columns
        real(dp), intent(in) :: eps !! relative tolerance, 0 <= eps < 1
        integer, intent(out) :: rank !! rank estimate
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! deterministic seed
        integer, intent(in), optional :: oversampling !! consecutive rejected probes
        complex(dp), intent(in), optional :: a(:, :) !! scaled dense matrix
        procedure(rank_complex_matvec), optional :: mv !! matrix-free forward action
        complex(dp), allocatable :: basis(:, :), x(:), y(:), work(:)
        complex(dp) :: alpha
        real(dp) :: ref, r
        integer(int64) :: state
        integer :: tries, extra, maxrank, misses, i, pass, callback_status
        logical :: seeded

        rank = 0
        status = linalg_status_invalid_shape
        extra = 8
        if (present(oversampling)) extra = oversampling
        if (m < 0 .or. n < 0 .or. extra < 1) return
        if (.not. ieee_is_finite(eps)) return
        if (eps < 0.0_dp .or. eps >= 1.0_dp) return
        if (present(a)) then
            if (size(a,1) /= m .or. size(a,2) /= n) return
        else
            if (.not. present(mv)) return
        end if
        status = linalg_status_success
        if (m == 0 .or. n == 0) return
        maxrank = min(m,n)
        allocate(basis(m,maxrank),x(n),y(m),work(m))
        state = 1_int64
        seeded = present(seed)
        if (seeded) state = modulo(seed, rng_modulus-1_int64) + 1_int64
        ref = 0.0_dp
        misses = 0
        do tries = 1, maxrank
            call take_probe(mv)
            if (status /= linalg_status_success) return
            if (rank == maxrank .or. misses >= extra) return
        end do
        do tries = 1, extra
            call take_probe(mv)
            if (status /= linalg_status_success) return
            if (rank == maxrank .or. misses >= extra) return
        end do
    contains
        subroutine take_probe(mv)
            procedure(rank_complex_matvec), optional :: mv !! matrix-free forward action
            call random_complex(x,state,seeded)
            r = norm_complex(x)
            if (r == 0.0_dp) then
                misses = misses + 1
                return
            end if
            x = x/r
            if (present(a)) then
                y = matmul(a,x)
            else
                call mv(x,y,callback_status)
                if (callback_status /= 0) then
                    status = linalg_status_no_convergence
                    rank = 0
                    return
                end if
            end if
            if (.not. finite_complex_vector(y)) then
                status = linalg_status_no_convergence
                rank = 0
                return
            end if
            r = norm_complex(y)
            if (.not. ieee_is_finite(r)) then
                status = linalg_status_no_convergence
                rank = 0
                return
            end if
            ref = max(ref,r)
            work = y
            do pass = 1, 2
                do i = 1, rank
                    alpha = dot_product(basis(:,i),work)
                    work = work - alpha*basis(:,i)
                end do
            end do
            r = norm_complex(work)
            if (.not. ieee_is_finite(r)) then
                status = linalg_status_no_convergence
                rank = 0
                return
            end if
            if (r > eps*ref .and. r > 0.0_dp) then
                rank = rank + 1
                basis(:,rank) = work/r
                misses = 0
            else
                misses = misses + 1
            end if
        end subroutine take_probe
    end subroutine rank_complex_core
end module scifort_linalg_estimate_rank
