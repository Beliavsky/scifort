! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Randomized power-method estimates of the matrix spectral norm.
! SciPy analogues: scipy.linalg.interpolative.estimate_spectral_norm and
! scipy.linalg.interpolative.estimate_spectral_norm_diff.
!
! A dense matrix or a pair of matrix/adjoint callbacks may be supplied.
! Callbacks allow the algorithm to work with sparse matrices, fast transforms
! and other implicit operators without forming a dense matrix. This source is
! independently implemented from the standard power-iteration identity.
module scifort_linalg_spectral_norm
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    implicit none
    private

    public :: estimate_spectral_norm, estimate_spectral_norm_diff
    public :: estimate_spectral_norm_operator, estimate_spectral_norm_operator_complex
    public :: estimate_spectral_norm_diff_operator, estimate_spectral_norm_diff_operator_complex
    public :: spectral_real_matvec, spectral_complex_matvec

    ! An operator supplies the matrix-vector product and its *adjoint* action.
    ! The output array is caller-allocated. A nonzero callback status aborts
    ! the estimate with linalg_status_no_convergence.
    abstract interface
        subroutine spectral_real_matvec(x, y, status)
            import dp
            real(dp), intent(in) :: x(:) !! input vector, either n or m elements
            real(dp), intent(out) :: y(:) !! result of op*x, either m or n elements
            integer, intent(out) :: status !! zero indicates successful product
        end subroutine spectral_real_matvec
        subroutine spectral_complex_matvec(x, y, status)
            import dp
            complex(dp), intent(in) :: x(:) !! input vector
            complex(dp), intent(out) :: y(:) !! result of op*x
            integer, intent(out) :: status !! zero indicates successful product
        end subroutine spectral_complex_matvec
    end interface

    interface estimate_spectral_norm
        module procedure spectral_real_dense
        module procedure spectral_complex_dense
    end interface estimate_spectral_norm
    interface estimate_spectral_norm_diff
        module procedure spectral_real_diff_dense
        module procedure spectral_complex_diff_dense
    end interface estimate_spectral_norm_diff
    interface estimate_spectral_norm_operator
        module procedure spectral_real_operator
    end interface estimate_spectral_norm_operator
    interface estimate_spectral_norm_diff_operator
        module procedure spectral_real_diff_operator
    end interface estimate_spectral_norm_diff_operator

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64

contains

    pure logical function finite_complex_matrix(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! matrix to check
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex_matrix

    pure logical function finite_complex_vector(x) result(ok)
        complex(dp), intent(in) :: x(:) !! vector to check
        ok = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function finite_complex_vector

    pure real(dp) function real_vector_norm(x) result(r)
        real(dp), intent(in) :: x(:) !! input vector
        real(dp) :: scale
        r = 0.0_dp
        if (size(x) == 0) return
        scale = maxval(abs(x))
        if (scale == 0.0_dp) return
        if (.not. ieee_is_finite(scale)) then
            r = scale
            return
        end if
        r = scale * sqrt(sum((x / scale)**2))
    end function real_vector_norm

    pure real(dp) function complex_vector_norm(x) result(r)
        complex(dp), intent(in) :: x(:) !! input vector
        real(dp) :: scale
        r = 0.0_dp
        if (size(x) == 0) return
        scale = maxval(abs(x))
        if (scale == 0.0_dp) return
        if (.not. ieee_is_finite(scale)) then
            r = scale
            return
        end if
        r = scale * sqrt(sum(abs(x / scale)**2))
    end function complex_vector_norm

    ! This generator does not change Fortran's global RNG when a seed is given.
    ! It is deliberately not bit-compatible with NumPy's Generator.
    real(dp) function seeded_draw(state) result(draw)
        integer(int64), intent(inout) :: state !! Park-Miller state, 1..modulus-1
        state = modulo(state * rng_multiplier, rng_modulus)
        draw = 2.0_dp * real(state,dp) / real(rng_modulus,dp) - 1.0_dp
    end function seeded_draw

    subroutine init_real_start(x, seed)
        real(dp), intent(out) :: x(:) !! randomized starting vector
        integer(int64), intent(in), optional :: seed !! local deterministic seed
        integer(int64) :: state
        integer :: j
        real(dp) :: u
        if (present(seed)) then
            state = modulo(seed, rng_modulus - 1_int64) + 1_int64
            do j = 1, size(x)
                x(j) = seeded_draw(state)
            end do
        else
            do j = 1, size(x)
                call random_number(u)
                x(j) = 2.0_dp * u - 1.0_dp
            end do
        end if
    end subroutine init_real_start

    subroutine init_complex_start(x, seed)
        complex(dp), intent(out) :: x(:) !! randomized complex starting vector
        integer(int64), intent(in), optional :: seed !! local deterministic seed
        integer(int64) :: state
        integer :: j
        real(dp) :: a, b
        if (present(seed)) then
            state = modulo(seed, rng_modulus - 1_int64) + 1_int64
            do j = 1, size(x)
                a = seeded_draw(state)
                b = seeded_draw(state)
                x(j) = cmplx(a,b,dp)
            end do
        else
            do j = 1, size(x)
                call random_number(a)
                call random_number(b)
                x(j) = cmplx(2.0_dp*a-1.0_dp,2.0_dp*b-1.0_dp,dp)
            end do
        end if
    end subroutine init_complex_start

    ! Real/complex power iterations share the same convention: `its` counts
    ! complete forward/adjoint cycles and the output uses one final forward
    ! multiply. This is an estimator, not a guaranteed error bound.
    subroutine real_power(m, n, estimate, status, its, seed, start, a, b, mv, rmv, mv_b, rmv_b)
        integer, intent(in) :: m !! operator rows, >= 0
        integer, intent(in) :: n !! operator columns, >= 0
        real(dp), intent(out) :: estimate !! estimated operator 2-norm
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: its !! number of iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        real(dp), intent(in), optional :: start(:) !! optional initial n-vector
        real(dp), intent(in), optional :: a(:, :) !! scaled dense first matrix
        real(dp), intent(in), optional :: b(:, :) !! scaled dense second matrix
        procedure(spectral_real_matvec), optional :: mv !! first forward product
        procedure(spectral_real_matvec), optional :: rmv !! first adjoint product
        procedure(spectral_real_matvec), optional :: mv_b !! second forward product
        procedure(spectral_real_matvec), optional :: rmv_b !! second adjoint product
        real(dp), allocatable :: x(:), y(:), z(:), temp(:)
        real(dp) :: r
        integer :: k, j, callback_status

        estimate = 0.0_dp
        status = linalg_status_invalid_shape
        k = 20
        if (present(its)) k = its
        if (m < 0 .or. n < 0 .or. k < 0) return
        if (present(start)) then
            if (size(start) /= n) return
            if (.not. all(ieee_is_finite(start))) return
        end if
        if (m == 0 .or. n == 0) then
            status = linalg_status_success
            return
        end if
        if (.not. present(a) .and. .not. present(mv)) then
            status = linalg_status_success
            return
        end if
        allocate(x(n),y(m),z(n),temp(max(m,n)))
        if (present(start)) then
            x = start
        else
            call init_real_start(x,seed)
        end if
        r = real_vector_norm(x)
        if (.not. ieee_is_finite(r) .or. r == 0.0_dp) return
        x = x / r
        status = linalg_status_no_convergence
        do j = 1, k
            call forward(x,y,callback_status,mv,mv_b)
            if (callback_status /= 0) return
            r = real_vector_norm(y)
            if (.not. ieee_is_finite(r)) return
            if (r == 0.0_dp) then
                status = linalg_status_success
                return
            end if
            y = y / r
            call adjoint(y,z,callback_status,rmv,rmv_b)
            if (callback_status /= 0) return
            r = real_vector_norm(z)
            if (.not. ieee_is_finite(r) .or. r == 0.0_dp) return
            x = z / r
        end do
        call forward(x,y,callback_status,mv,mv_b)
        if (callback_status /= 0) return
        r = real_vector_norm(y)
        if (.not. ieee_is_finite(r)) return
        estimate = r
        status = linalg_status_success

    contains
        subroutine forward(xin, yout, istat, mv, mv_b)
            procedure(spectral_real_matvec), optional :: mv !! first forward product
            procedure(spectral_real_matvec), optional :: mv_b !! second forward product
            real(dp), intent(in) :: xin(:) !! input n-vector
            real(dp), intent(out) :: yout(:) !! output m-vector
            integer, intent(out) :: istat !! zero means success
            istat = 0
            if (present(a)) then
                yout = matmul(a,xin)
                if (present(b)) yout = yout - matmul(b,xin)
            else
                if (.not. present(mv)) then
                    istat = 1
                    return
                end if
                call mv(xin,yout,istat)
                if (istat /= 0) return
                if (present(mv_b)) then
                    call mv_b(xin,temp(1:m),istat)
                    if (istat == 0) yout = yout - temp(1:m)
                end if
            end if
            if (.not. all(ieee_is_finite(yout))) istat = 1
        end subroutine forward

        subroutine adjoint(yin, zout, istat, rmv, rmv_b)
            procedure(spectral_real_matvec), optional :: rmv !! first adjoint product
            procedure(spectral_real_matvec), optional :: rmv_b !! second adjoint product
            real(dp), intent(in) :: yin(:) !! input m-vector
            real(dp), intent(out) :: zout(:) !! output n-vector
            integer, intent(out) :: istat !! zero means success
            istat = 0
            if (present(a)) then
                zout = matmul(transpose(a),yin)
                if (present(b)) zout = zout - matmul(transpose(b),yin)
            else
                if (.not. present(rmv)) then
                    istat = 1
                    return
                end if
                call rmv(yin,zout,istat)
                if (istat /= 0) return
                if (present(rmv_b)) then
                    call rmv_b(yin,temp(1:n),istat)
                    if (istat == 0) zout = zout - temp(1:n)
                end if
            end if
            if (.not. all(ieee_is_finite(zout))) istat = 1
        end subroutine adjoint
    end subroutine real_power

    subroutine complex_power(m, n, estimate, status, its, seed, start, a, b, mv, rmv, mv_b, rmv_b)
        integer, intent(in) :: m !! operator rows, >= 0
        integer, intent(in) :: n !! operator columns, >= 0
        real(dp), intent(out) :: estimate !! estimated operator 2-norm
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: its !! number of iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        complex(dp), intent(in), optional :: start(:) !! optional initial n-vector
        complex(dp), intent(in), optional :: a(:, :) !! scaled dense first matrix
        complex(dp), intent(in), optional :: b(:, :) !! scaled dense second matrix
        procedure(spectral_complex_matvec), optional :: mv !! first forward product
        procedure(spectral_complex_matvec), optional :: rmv !! first adjoint product
        procedure(spectral_complex_matvec), optional :: mv_b !! second forward product
        procedure(spectral_complex_matvec), optional :: rmv_b !! second adjoint product
        complex(dp), allocatable :: x(:), y(:), z(:), temp(:)
        real(dp) :: r
        integer :: k, j, callback_status

        estimate = 0.0_dp
        status = linalg_status_invalid_shape
        k = 20
        if (present(its)) k = its
        if (m < 0 .or. n < 0 .or. k < 0) return
        if (present(start)) then
            if (size(start) /= n) return
            if (.not. finite_complex_vector(start)) return
        end if
        if (m == 0 .or. n == 0) then
            status = linalg_status_success
            return
        end if
        if (.not. present(a) .and. .not. present(mv)) then
            status = linalg_status_success
            return
        end if
        allocate(x(n),y(m),z(n),temp(max(m,n)))
        if (present(start)) then
            x = start
        else
            call init_complex_start(x,seed)
        end if
        r = complex_vector_norm(x)
        if (.not. ieee_is_finite(r) .or. r == 0.0_dp) return
        x = x / r
        status = linalg_status_no_convergence
        do j = 1, k
            call forward(x,y,callback_status,mv,mv_b)
            if (callback_status /= 0) return
            r = complex_vector_norm(y)
            if (.not. ieee_is_finite(r)) return
            if (r == 0.0_dp) then
                status = linalg_status_success
                return
            end if
            y = y / r
            call adjoint(y,z,callback_status,rmv,rmv_b)
            if (callback_status /= 0) return
            r = complex_vector_norm(z)
            if (.not. ieee_is_finite(r) .or. r == 0.0_dp) return
            x = z / r
        end do
        call forward(x,y,callback_status,mv,mv_b)
        if (callback_status /= 0) return
        r = complex_vector_norm(y)
        if (.not. ieee_is_finite(r)) return
        estimate = r
        status = linalg_status_success

    contains
        subroutine forward(xin, yout, istat, mv, mv_b)
            procedure(spectral_complex_matvec), optional :: mv !! first forward product
            procedure(spectral_complex_matvec), optional :: mv_b !! second forward product
            complex(dp), intent(in) :: xin(:) !! input n-vector
            complex(dp), intent(out) :: yout(:) !! output m-vector
            integer, intent(out) :: istat !! zero means success
            istat = 0
            if (present(a)) then
                yout = matmul(a,xin)
                if (present(b)) yout = yout - matmul(b,xin)
            else
                if (.not. present(mv)) then
                    istat = 1
                    return
                end if
                call mv(xin,yout,istat)
                if (istat /= 0) return
                if (present(mv_b)) then
                    call mv_b(xin,temp(1:m),istat)
                    if (istat == 0) yout = yout - temp(1:m)
                end if
            end if
            if (.not. finite_complex_vector(yout)) istat = 1
        end subroutine forward

        subroutine adjoint(yin, zout, istat, rmv, rmv_b)
            procedure(spectral_complex_matvec), optional :: rmv !! first adjoint product
            procedure(spectral_complex_matvec), optional :: rmv_b !! second adjoint product
            complex(dp), intent(in) :: yin(:) !! input m-vector
            complex(dp), intent(out) :: zout(:) !! output n-vector
            integer, intent(out) :: istat !! zero means success
            istat = 0
            if (present(a)) then
                zout = matmul(transpose(conjg(a)),yin)
                if (present(b)) zout = zout - matmul(transpose(conjg(b)),yin)
            else
                if (.not. present(rmv)) then
                    istat = 1
                    return
                end if
                call rmv(yin,zout,istat)
                if (istat /= 0) return
                if (present(rmv_b)) then
                    call rmv_b(yin,temp(1:n),istat)
                    if (istat == 0) zout = zout - temp(1:n)
                end if
            end if
            if (.not. finite_complex_vector(zout)) istat = 1
        end subroutine adjoint
    end subroutine complex_power

    subroutine spectral_real_dense(a, estimate, status, its, seed, start)
        real(dp), intent(in) :: a(:, :) !! m-by-n real matrix
        real(dp), intent(out) :: estimate !! estimated spectral norm
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in), optional :: its !! power iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        real(dp), intent(in), optional :: start(:) !! starting n-vector
        real(dp) :: scale
        estimate = 0.0_dp
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        scale = 0.0_dp
        if (size(a) > 0) scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            call real_power(size(a,1),size(a,2),estimate,status,its,seed,start)
        else
            call real_power(size(a,1),size(a,2),estimate,status,its,seed,start,a=a/scale)
            if (status == linalg_status_success) then
                estimate = scale * estimate
                if (.not. ieee_is_finite(estimate)) then
                    estimate = 0.0_dp
                    status = linalg_status_no_convergence
                end if
            end if
        end if
    end subroutine spectral_real_dense

    subroutine spectral_complex_dense(a, estimate, status, its, seed, start)
        complex(dp), intent(in) :: a(:, :) !! m-by-n complex matrix
        real(dp), intent(out) :: estimate !! estimated spectral norm
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in), optional :: its !! power iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        complex(dp), intent(in), optional :: start(:) !! starting n-vector
        real(dp) :: scale
        estimate = 0.0_dp
        status = linalg_status_no_convergence
        if (.not. finite_complex_matrix(a)) return
        scale = 0.0_dp
        if (size(a) > 0) scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            call complex_power(size(a,1),size(a,2),estimate,status,its,seed,start)
        else
            call complex_power(size(a,1),size(a,2),estimate,status,its,seed,start,a=a/scale)
            if (status == linalg_status_success) then
                estimate = scale * estimate
                if (.not. ieee_is_finite(estimate)) then
                    estimate = 0.0_dp
                    status = linalg_status_no_convergence
                end if
            end if
        end if
    end subroutine spectral_complex_dense

    subroutine spectral_real_diff_dense(a, b, estimate, status, its, seed, start)
        real(dp), intent(in) :: a(:, :) !! first m-by-n matrix
        real(dp), intent(in) :: b(:, :) !! second m-by-n matrix
        real(dp), intent(out) :: estimate !! estimated 2-norm of A-B
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in), optional :: its !! power iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        real(dp), intent(in), optional :: start(:) !! starting n-vector
        real(dp) :: scale, max_input
        real(dp), allocatable :: difference(:, :)
        estimate = 0.0_dp
        status = linalg_status_invalid_shape
        if (any(shape(a) /= shape(b))) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        max_input = 0.0_dp
        if (size(a) > 0) max_input = max(maxval(abs(a)),maxval(abs(b)))
        if (max_input <= huge(1.0_dp)/4.0_dp) then
            ! Subtract before scaling to retain small differences between
            ! large, nearly equal matrix entries whenever safe from overflow.
            difference = a - b
        else
            difference = a / max_input - b / max_input
        end if
        scale = 0.0_dp
        if (size(difference) > 0) scale = maxval(abs(difference))
        if (scale == 0.0_dp) then
            call real_power(size(a,1),size(a,2),estimate,status,its,seed,start)
        else
            call real_power(size(a,1),size(a,2),estimate,status,its,seed,start, &
                a=difference/scale)
            if (status == linalg_status_success) then
                estimate = scale * estimate
                if (max_input > huge(1.0_dp)/4.0_dp) estimate = max_input * estimate
                if (.not. ieee_is_finite(estimate)) then
                    estimate = 0.0_dp
                    status = linalg_status_no_convergence
                end if
            end if
        end if
    end subroutine spectral_real_diff_dense

    subroutine spectral_complex_diff_dense(a, b, estimate, status, its, seed, start)
        complex(dp), intent(in) :: a(:, :) !! first m-by-n matrix
        complex(dp), intent(in) :: b(:, :) !! second m-by-n matrix
        real(dp), intent(out) :: estimate !! estimated 2-norm of A-B
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in), optional :: its !! power iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        complex(dp), intent(in), optional :: start(:) !! starting n-vector
        real(dp) :: scale, max_input
        complex(dp), allocatable :: difference(:, :)
        estimate = 0.0_dp
        status = linalg_status_invalid_shape
        if (any(shape(a) /= shape(b))) return
        status = linalg_status_no_convergence
        if (.not. finite_complex_matrix(a) .or. .not. finite_complex_matrix(b)) return
        max_input = 0.0_dp
        if (size(a) > 0) max_input = max(maxval(abs(a)),maxval(abs(b)))
        if (max_input <= huge(1.0_dp)/4.0_dp) then
            difference = a - b
        else
            difference = a / max_input - b / max_input
        end if
        scale = 0.0_dp
        if (size(difference) > 0) scale = maxval(abs(difference))
        if (scale == 0.0_dp) then
            call complex_power(size(a,1),size(a,2),estimate,status,its,seed,start)
        else
            call complex_power(size(a,1),size(a,2),estimate,status,its,seed,start, &
                a=difference/scale)
            if (status == linalg_status_success) then
                estimate = scale * estimate
                if (max_input > huge(1.0_dp)/4.0_dp) estimate = max_input * estimate
                if (.not. ieee_is_finite(estimate)) then
                    estimate = 0.0_dp
                    status = linalg_status_no_convergence
                end if
            end if
        end if
    end subroutine spectral_complex_diff_dense

    subroutine spectral_real_operator(m, n, matvec, rmatvec, estimate, status, its, seed, start)
        integer, intent(in) :: m !! matrix rows
        integer, intent(in) :: n !! matrix columns
        procedure(spectral_real_matvec) :: matvec !! computes A*x
        procedure(spectral_real_matvec) :: rmatvec !! computes transpose(A)*y
        real(dp), intent(out) :: estimate !! estimated spectral norm
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in), optional :: its !! power iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        real(dp), intent(in), optional :: start(:) !! starting n-vector
        call real_power(m,n,estimate,status,its,seed,start,mv=matvec,rmv=rmatvec)
    end subroutine spectral_real_operator

    subroutine estimate_spectral_norm_operator_complex(m, n, matvec, rmatvec, estimate, status, its, seed, start)
        integer, intent(in) :: m !! matrix rows
        integer, intent(in) :: n !! matrix columns
        procedure(spectral_complex_matvec) :: matvec !! computes A*x
        procedure(spectral_complex_matvec) :: rmatvec !! computes A^H*y
        real(dp), intent(out) :: estimate !! estimated spectral norm
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in), optional :: its !! power iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        complex(dp), intent(in), optional :: start(:) !! starting n-vector
        call complex_power(m,n,estimate,status,its,seed,start,mv=matvec,rmv=rmatvec)
    end subroutine estimate_spectral_norm_operator_complex

    subroutine spectral_real_diff_operator(m, n, matvec_a, rmatvec_a, matvec_b, &
        rmatvec_b, estimate, status, its, seed, start)
        integer, intent(in) :: m !! matrix rows
        integer, intent(in) :: n !! matrix columns
        procedure(spectral_real_matvec) :: matvec_a !! computes A*x
        procedure(spectral_real_matvec) :: rmatvec_a !! computes transpose(A)*y
        procedure(spectral_real_matvec) :: matvec_b !! computes B*x
        procedure(spectral_real_matvec) :: rmatvec_b !! computes transpose(B)*y
        real(dp), intent(out) :: estimate !! estimated norm of A-B
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in), optional :: its !! power iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        real(dp), intent(in), optional :: start(:) !! starting n-vector
        call real_power(m,n,estimate,status,its,seed,start, &
            mv=matvec_a,rmv=rmatvec_a,mv_b=matvec_b,rmv_b=rmatvec_b)
    end subroutine spectral_real_diff_operator

    subroutine estimate_spectral_norm_diff_operator_complex(m, n, matvec_a, rmatvec_a, matvec_b, &
        rmatvec_b, estimate, status, its, seed, start)
        integer, intent(in) :: m !! matrix rows
        integer, intent(in) :: n !! matrix columns
        procedure(spectral_complex_matvec) :: matvec_a !! computes A*x
        procedure(spectral_complex_matvec) :: rmatvec_a !! computes A^H*y
        procedure(spectral_complex_matvec) :: matvec_b !! computes B*x
        procedure(spectral_complex_matvec) :: rmatvec_b !! computes B^H*y
        real(dp), intent(out) :: estimate !! estimated norm of A-B
        integer, intent(out) :: status !! linalg_status_* result
        integer, intent(in), optional :: its !! power iterations, default 20
        integer(int64), intent(in), optional :: seed !! deterministic local seed
        complex(dp), intent(in), optional :: start(:) !! starting n-vector
        call complex_power(m,n,estimate,status,its,seed,start, &
            mv=matvec_a,rmv=rmatvec_a,mv_b=matvec_b,rmv_b=rmatvec_b)
    end subroutine estimate_spectral_norm_diff_operator_complex
end module scifort_linalg_spectral_norm
