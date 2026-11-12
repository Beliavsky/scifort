! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Randomized interpolative decomposition of an implicit real/complex matrix.
! Analogue of scipy.linalg.interpolative.interp_decomp(LinearOperator, ...).
!
! Only forward (A*x) and adjoint (A^H*x) callbacks are needed. A short
! Gaussian row sketch chooses columns; k forward products recover the exact
! skeleton, and k adjoint products refit the interpolation coefficients in
! its orthogonal basis. The adaptive interface uses independent Gaussian
! probes to ESTIMATE, not certify, the relative Frobenius residual.
!
! A(:,idx(k+1:n)) ~= A(:,idx(1:k)) * proj, with one-based idx.
module scifort_linalg_interpolative_operator
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr_factor, qr
    use scifort_linalg_svd, only : svd_lstsq
    use scifort_linalg_interpolative, only : interp_decomp
    implicit none
    private

    public :: randomized_interp_decomp_operator
    public :: randomized_interp_decomp_operator_complex
    public :: id_real_matvec, id_complex_matvec

    abstract interface
        subroutine id_real_matvec(x, y, status)
            import dp
            real(dp), intent(in) :: x(:) !! input of length n (forward) or m (adjoint)
            real(dp), intent(out) :: y(:) !! output of length m (forward) or n (adjoint)
            integer, intent(out) :: status !! 0 on success; nonzero on failure
        end subroutine id_real_matvec
        subroutine id_complex_matvec(x, y, status)
            import dp
            complex(dp), intent(in) :: x(:) !! input of length n (forward) or m (adjoint)
            complex(dp), intent(out) :: y(:) !! output of length m (forward) or n (adjoint)
            integer, intent(out) :: status !! 0 on success; nonzero on failure
        end subroutine id_complex_matvec
    end interface

    interface randomized_interp_decomp_operator
        module procedure real_fixed
        module procedure real_tolerance
    end interface
    interface randomized_interp_decomp_operator_complex
        module procedure complex_fixed
        module procedure complex_tolerance
    end interface

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64
    integer(int64), parameter :: default_seed = 20261110_int64
    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    pure integer(int64) function seed_state(seed) result(state)
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        state = modulo(default_seed, rng_modulus-1_int64) + 1_int64
        if (present(seed)) state = modulo(seed, rng_modulus-1_int64) + 1_int64
    end function seed_state

    real(dp) function normal_draw(state) result(value)
        integer(int64), intent(inout) :: state !! local Park-Miller RNG state
        real(dp) :: u, v
        state = modulo(state*rng_multiplier,rng_modulus)
        u = real(state,dp)/real(rng_modulus,dp)
        state = modulo(state*rng_multiplier,rng_modulus)
        v = real(state,dp)/real(rng_modulus,dp)
        value = sqrt(-2.0_dp*log(u))*cos(2.0_dp*pi*v)
    end function normal_draw

    pure logical function finite_complex(x) result(ok)
        complex(dp), intent(in) :: x(:) !! input vector
        ok = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function finite_complex

    ! Scaled sum of squares; avoids overflow in probing very large matrices.
    pure subroutine add_square(x, scale, ssq)
        real(dp), intent(in) :: x !! one real component to accumulate
        real(dp), intent(inout) :: scale !! largest magnitude so far
        real(dp), intent(inout) :: ssq !! scaled sum of squares, initially 1
        real(dp) :: ax
        ax = abs(x)
        if (ax == 0.0_dp) return
        if (scale < ax) then
            ssq = 1.0_dp + ssq*(scale/ax)**2
            scale = ax
        else
            ssq = ssq + (ax/scale)**2
        end if
    end subroutine add_square

    subroutine build_sketch_real(m, n, rows, adjoint, state, y, status)
        integer, intent(in) :: m !! original number of rows
        integer, intent(in) :: n !! original number of columns
        integer, intent(in) :: rows !! number of sketch rows, positive
        procedure(id_real_matvec) :: adjoint !! A^T vector product
        integer(int64), intent(inout) :: state !! local RNG state
        real(dp), allocatable, intent(out) :: y(:, :) !! random row sketch, rows-by-n
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: omega(:), z(:)
        integer :: i, j, callback_status
        allocate(y(rows,n),omega(m),z(n))
        status = linalg_status_no_convergence
        do i = 1, rows
            do j = 1, m
                omega(j) = normal_draw(state)
            end do
            call adjoint(omega,z,callback_status)
            if (callback_status /= 0) return
            if (.not. all(ieee_is_finite(z))) return
            y(i,:) = z
        end do
        status = linalg_status_success
    end subroutine build_sketch_real

    subroutine build_sketch_complex(m, n, rows, adjoint, state, y, status)
        integer, intent(in) :: m !! original number of rows
        integer, intent(in) :: n !! original number of columns
        integer, intent(in) :: rows !! number of sketch rows, positive
        procedure(id_complex_matvec) :: adjoint !! A^H vector product
        integer(int64), intent(inout) :: state !! local RNG state
        complex(dp), allocatable, intent(out) :: y(:, :) !! random row sketch, rows-by-n
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: omega(:), z(:)
        integer :: i, j, callback_status
        allocate(y(rows,n),omega(m),z(n))
        status = linalg_status_no_convergence
        do i = 1, rows
            do j = 1, m
                omega(j) = cmplx(normal_draw(state),normal_draw(state),dp)/sqrt(2.0_dp)
            end do
            call adjoint(omega,z,callback_status)
            if (callback_status /= 0) return
            if (.not. finite_complex(z)) return
            ! y(i,:) = omega^H A = conjugate(A^H omega)^T.
            y(i,:) = conjg(z)
        end do
        status = linalg_status_success
    end subroutine build_sketch_complex

    subroutine refit_real(m, n, k, idx, forward, adjoint, proj, skeleton, status)
        integer, intent(in) :: m !! operator row count
        integer, intent(in) :: n !! operator column count
        integer, intent(in) :: k !! selected column count
        integer, intent(in) :: idx(:) !! one-based column permutation
        procedure(id_real_matvec) :: forward !! forward product callback
        procedure(id_real_matvec) :: adjoint !! adjoint product callback
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) coefficients
        real(dp), allocatable, intent(out) :: skeleton(:, :) !! selected original columns
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: e(:), z(:), bscaled(:, :), q(:, :), r(:, :), h(:, :)
        real(dp) :: scale
        integer :: i, callback_status

        allocate(proj(k,n-k),skeleton(m,k))
        proj = 0.0_dp
        skeleton = 0.0_dp
        status = linalg_status_success
        if (k == 0) return
        allocate(e(n),z(n))
        status = linalg_status_no_convergence
        do i = 1, k
            e = 0.0_dp
            e(idx(i)) = 1.0_dp
            call forward(e,skeleton(:,i),callback_status)
            if (callback_status /= 0) return
            if (.not. all(ieee_is_finite(skeleton(:,i)))) return
        end do
        if (n == k) then
            status = linalg_status_success
            return
        end if
        scale = maxval(abs(skeleton))
        if (scale == 0.0_dp) then
            status = linalg_status_success
            return
        end if
        bscaled = skeleton/scale
        allocate(q(m,k),r(k,k),h(k,n))
        call qr(bscaled,q,r,status)
        if (status /= linalg_status_success) return
        do i = 1, k
            call adjoint(q(:,i),z,callback_status)
            if (callback_status /= 0) then
                status = linalg_status_no_convergence
                return
            end if
            if (.not. all(ieee_is_finite(z))) then
                status = linalg_status_no_convergence
                return
            end if
            h(i,:) = z/scale
        end do
        if (.not. all(ieee_is_finite(h))) then
            status = linalg_status_no_convergence
            return
        end if
        call svd_lstsq(h(:,idx(1:k)),h(:,idx(k+1:n)),proj,status)
    end subroutine refit_real

    subroutine refit_complex(m, n, k, idx, forward, adjoint, proj, skeleton, status)
        integer, intent(in) :: m !! operator row count
        integer, intent(in) :: n !! operator column count
        integer, intent(in) :: k !! selected column count
        integer, intent(in) :: idx(:) !! one-based column permutation
        procedure(id_complex_matvec) :: forward !! forward product callback
        procedure(id_complex_matvec) :: adjoint !! Hermitian-adjoint callback
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) coefficients
        complex(dp), allocatable, intent(out) :: skeleton(:, :) !! selected original columns
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: e(:), z(:), bscaled(:, :), q(:, :), r(:, :), h(:, :)
        real(dp) :: scale
        integer :: i, callback_status

        allocate(proj(k,n-k),skeleton(m,k))
        proj = cmplx(0.0_dp,0.0_dp,dp)
        skeleton = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_success
        if (k == 0) return
        allocate(e(n),z(n))
        status = linalg_status_no_convergence
        do i = 1, k
            e = cmplx(0.0_dp,0.0_dp,dp)
            e(idx(i)) = cmplx(1.0_dp,0.0_dp,dp)
            call forward(e,skeleton(:,i),callback_status)
            if (callback_status /= 0) return
            if (.not. finite_complex(skeleton(:,i))) return
        end do
        if (n == k) then
            status = linalg_status_success
            return
        end if
        scale = max(maxval(abs(real(skeleton,dp))),maxval(abs(aimag(skeleton))))
        if (scale == 0.0_dp) then
            status = linalg_status_success
            return
        end if
        bscaled = skeleton/scale
        allocate(q(m,k),r(k,k),h(k,n))
        call qr(bscaled,q,r,status)
        if (status /= linalg_status_success) return
        do i = 1, k
            call adjoint(q(:,i),z,callback_status)
            if (callback_status /= 0) then
                status = linalg_status_no_convergence
                return
            end if
            if (.not. finite_complex(z)) then
                status = linalg_status_no_convergence
                return
            end if
            h(i,:) = conjg(z/scale)
        end do
        if (.not. all(ieee_is_finite(real(h,dp))) .or. .not. all(ieee_is_finite(aimag(h)))) then
            status = linalg_status_no_convergence
            return
        end if
        call svd_lstsq(h(:,idx(1:k)),h(:,idx(k+1:n)),proj,status)
    end subroutine refit_complex

    subroutine probe_real(m, n, k, idx, proj, skeleton, forward, state, probes, err, status)
        integer, intent(in) :: m !! operator row count
        integer, intent(in) :: n !! operator column count
        integer, intent(in) :: k !! selected column count
        integer, intent(in) :: idx(:) !! one-based column permutation
        real(dp), intent(in) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        real(dp), intent(in) :: skeleton(:, :) !! m-by-k selected original columns
        procedure(id_real_matvec) :: forward !! A*x callback
        integer(int64), intent(inout) :: state !! local RNG state
        integer, intent(in) :: probes !! independent Gaussian test vectors
        real(dp), intent(out) :: err !! estimated relative Frobenius residual
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: x(:), y(:), residual(:), pvec(:)
        real(dp) :: sa, se, qa, qe
        integer :: i, j, callback_status
        err = huge(1.0_dp)
        status = linalg_status_no_convergence
        allocate(x(n),y(m),residual(m),pvec(k))
        sa = 0.0_dp
        se = 0.0_dp
        qa = 1.0_dp
        qe = 1.0_dp
        do i = 1, probes
            do j = 1, n
                x(j) = normal_draw(state)
            end do
            call forward(x,y,callback_status)
            if (callback_status /= 0) return
            if (.not. all(ieee_is_finite(y))) return
            if (k > 0) pvec = x(idx(1:k)) + matmul(proj,x(idx(k+1:n)))
            residual = y
            if (k > 0) residual = y - matmul(skeleton,pvec)
            if (.not. all(ieee_is_finite(residual))) return
            do j = 1, m
                call add_square(y(j),sa,qa)
                call add_square(residual(j),se,qe)
            end do
        end do
        err = 0.0_dp
        if (sa == 0.0_dp) then
            if (se /= 0.0_dp) return
        else
            err = (se/sa)*sqrt(qe/qa)
        end if
        if (.not. ieee_is_finite(err)) return
        status = linalg_status_success
    end subroutine probe_real

    subroutine probe_complex(m, n, k, idx, proj, skeleton, forward, state, probes, err, status)
        integer, intent(in) :: m !! operator row count
        integer, intent(in) :: n !! operator column count
        integer, intent(in) :: k !! selected column count
        integer, intent(in) :: idx(:) !! one-based column permutation
        complex(dp), intent(in) :: proj(:, :) !! k-by-(n-k) interpolation coefficients
        complex(dp), intent(in) :: skeleton(:, :) !! m-by-k selected original columns
        procedure(id_complex_matvec) :: forward !! A*x callback
        integer(int64), intent(inout) :: state !! local RNG state
        integer, intent(in) :: probes !! independent Gaussian test vectors
        real(dp), intent(out) :: err !! estimated relative Frobenius residual
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: x(:), y(:), residual(:), pvec(:)
        real(dp) :: sa, se, qa, qe
        integer :: i, j, callback_status
        err = huge(1.0_dp)
        status = linalg_status_no_convergence
        allocate(x(n),y(m),residual(m),pvec(k))
        sa = 0.0_dp
        se = 0.0_dp
        qa = 1.0_dp
        qe = 1.0_dp
        do i = 1, probes
            do j = 1, n
                x(j) = cmplx(normal_draw(state),normal_draw(state),dp)/sqrt(2.0_dp)
            end do
            call forward(x,y,callback_status)
            if (callback_status /= 0) return
            if (.not. finite_complex(y)) return
            if (k > 0) pvec = x(idx(1:k)) + matmul(proj,x(idx(k+1:n)))
            residual = y
            if (k > 0) residual = y - matmul(skeleton,pvec)
            if (.not. finite_complex(residual)) return
            do j = 1, m
                call add_square(real(y(j),dp),sa,qa)
                call add_square(aimag(y(j)),sa,qa)
                call add_square(real(residual(j),dp),se,qe)
                call add_square(aimag(residual(j)),se,qe)
            end do
        end do
        err = 0.0_dp
        if (sa == 0.0_dp) then
            if (se /= 0.0_dp) return
        else
            err = (se/sa)*sqrt(qe/qa)
        end if
        if (.not. ieee_is_finite(err)) return
        status = linalg_status_success
    end subroutine probe_complex

    subroutine real_fixed(m, n, k, matvec, rmatvec, idx, proj, status, seed, oversampling, skeleton)
        integer, intent(in) :: m !! number of rows of the implicit A
        integer, intent(in) :: n !! number of columns of the implicit A
        integer, intent(in) :: k !! desired skeleton size, 0..min(m,n)
        procedure(id_real_matvec) :: matvec !! A*x callback
        procedure(id_real_matvec) :: rmatvec !! A^T*x callback
        integer, allocatable, intent(out) :: idx(:) !! one-based n-column permutation
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! extra sketch rows, default 8
        real(dp), allocatable, intent(out), optional :: skeleton(:, :) !! m-by-k exact selected columns
        real(dp), allocatable :: y(:, :), af(:, :), tau(:), tproj(:, :), b(:, :)
        integer, allocatable :: p(:)
        integer(int64) :: state
        integer :: over, rows, j
        allocate(idx(0),proj(0,0))
        if (present(skeleton)) allocate(skeleton(0,0))
        status = linalg_status_invalid_shape
        over = 8
        if (present(oversampling)) over = oversampling
        if (m < 0 .or. n < 0 .or. k < 0 .or. k > min(m,n) .or. over < 0) return
        if (k == 0) then
            deallocate(idx,proj)
            allocate(idx(n),proj(0,n))
            idx = [(j,j=1,n)]
            if (present(skeleton)) then
                deallocate(skeleton)
                allocate(skeleton(m,0))
            end if
            status = linalg_status_success
            return
        end if
        state = seed_state(seed)
        rows = k + min(over,min(m,n)-k)
        call build_sketch_real(m,n,rows,rmatvec,state,y,status)
        if (status /= linalg_status_success) return
        allocate(af(rows,n),tau(rows),p(n))
        call qr_factor(y,af,tau,status,p)
        if (status /= linalg_status_success) return
        call refit_real(m,n,k,p,matvec,rmatvec,tproj,b,status)
        if (status /= linalg_status_success) return
        deallocate(idx,proj)
        call move_alloc(p,idx)
        call move_alloc(tproj,proj)
        if (present(skeleton)) then
            deallocate(skeleton)
            call move_alloc(b,skeleton)
        end if
    end subroutine real_fixed

    subroutine complex_fixed(m, n, k, matvec, rmatvec, idx, proj, status, seed, oversampling, skeleton)
        integer, intent(in) :: m !! number of rows of implicit A
        integer, intent(in) :: n !! number of columns of implicit A
        integer, intent(in) :: k !! desired skeleton size, 0..min(m,n)
        procedure(id_complex_matvec) :: matvec !! A*x callback
        procedure(id_complex_matvec) :: rmatvec !! A^H*x callback
        integer, allocatable, intent(out) :: idx(:) !! one-based n-column permutation
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! extra sketch rows, default 8
        complex(dp), allocatable, intent(out), optional :: skeleton(:, :) !! m-by-k exact selected columns
        complex(dp), allocatable :: y(:, :), af(:, :), tau(:), tproj(:, :), b(:, :)
        integer, allocatable :: p(:)
        integer(int64) :: state
        integer :: over, rows, j
        allocate(idx(0),proj(0,0))
        if (present(skeleton)) allocate(skeleton(0,0))
        status = linalg_status_invalid_shape
        over = 8
        if (present(oversampling)) over = oversampling
        if (m < 0 .or. n < 0 .or. k < 0 .or. k > min(m,n) .or. over < 0) return
        if (k == 0) then
            deallocate(idx,proj)
            allocate(idx(n),proj(0,n))
            idx = [(j,j=1,n)]
            if (present(skeleton)) then
                deallocate(skeleton)
                allocate(skeleton(m,0))
            end if
            status = linalg_status_success
            return
        end if
        state = seed_state(seed)
        rows = k + min(over,min(m,n)-k)
        call build_sketch_complex(m,n,rows,rmatvec,state,y,status)
        if (status /= linalg_status_success) return
        allocate(af(rows,n),tau(rows),p(n))
        call qr_factor(y,af,tau,status,p)
        if (status /= linalg_status_success) return
        call refit_complex(m,n,k,p,matvec,rmatvec,tproj,b,status)
        if (status /= linalg_status_success) return
        deallocate(idx,proj)
        call move_alloc(p,idx)
        call move_alloc(tproj,proj)
        if (present(skeleton)) then
            deallocate(skeleton)
            call move_alloc(b,skeleton)
        end if
    end subroutine complex_fixed

    subroutine real_tolerance(m, n, eps, matvec, rmatvec, k, idx, proj, status, &
        seed, oversampling, verification_probes, skeleton, estimated_error)
        integer, intent(in) :: m !! number of rows of implicit A
        integer, intent(in) :: n !! number of columns of implicit A
        real(dp), intent(in) :: eps !! relative Frobenius error target, 0 <= eps < 1
        procedure(id_real_matvec) :: matvec !! A*x callback
        procedure(id_real_matvec) :: rmatvec !! A^T*x callback
        integer, intent(out) :: k !! selected numerical rank
        integer, allocatable, intent(out) :: idx(:) !! one-based n-column permutation
        real(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! initial sketch size hint, default 8
        integer, intent(in), optional :: verification_probes !! Gaussian residual probes, default 16
        real(dp), allocatable, intent(out), optional :: skeleton(:, :) !! m-by-k exact selected columns
        real(dp), intent(out), optional :: estimated_error !! sampled relative Frobenius residual
        real(dp), allocatable :: y(:, :), unused(:, :), tproj(:, :), b(:, :)
        integer, allocatable :: p(:)
        integer(int64) :: state, probe_state
        integer :: over, tests, rows, rmax, candidate, lo, hi, mid, j
        real(dp) :: err, floor_tolerance, threshold

        k = 0
        allocate(idx(0),proj(0,0))
        if (present(skeleton)) allocate(skeleton(0,0))
        if (present(estimated_error)) estimated_error = huge(1.0_dp)
        status = linalg_status_invalid_shape
        over = 8
        tests = 16
        if (present(oversampling)) over = oversampling
        if (present(verification_probes)) tests = verification_probes
        if (m < 0 .or. n < 0 .or. over < 0 .or. tests < 1) return
        if (.not. ieee_is_finite(eps) .or. eps < 0.0_dp .or. eps >= 1.0_dp) return
        rmax = min(m,n)
        if (rmax == 0) then
            call empty_result()
            return
        end if
        state = seed_state(seed)
        rows = min(rmax,max(8,min(rmax,max(1,over))))
        floor_tolerance = 100.0_dp*epsilon(1.0_dp)*real(max(m,n),dp)
        threshold = max(min(0.85_dp,max(0.4_dp,1.0_dp-1.5_dp/sqrt(real(tests,dp))))*eps, &
            floor_tolerance)
        do
            call build_sketch_real(m,n,rows,rmatvec,state,y,status)
            if (status /= linalg_status_success) return
            call interp_decomp(y,eps,candidate,p,unused,status)
            if (status /= linalg_status_success) return
            ! Keep verification vectors fixed while comparing ranks. This
            ! makes errors comparable for nested skeletons and avoids an
            ! accidental rank jump caused only by fresh random probes.
            probe_state = state
            call refit_real(m,n,candidate,p,matvec,rmatvec,tproj,b,status)
            if (status /= linalg_status_success) return
            call probe_real(m,n,candidate,p,tproj,b,matvec,probe_state,tests,err,status)
            if (status /= linalg_status_success) return
            if (err <= threshold) exit
            if (candidate < rows) then
                ! Try all available sketch rows before drawing a bigger
                ! sketch; then bisect to avoid returning needlessly high rank.
                hi = rows
                probe_state = state
                call refit_real(m,n,hi,p,matvec,rmatvec,tproj,b,status)
                if (status /= linalg_status_success) return
                call probe_real(m,n,hi,p,tproj,b,matvec,probe_state,tests,err,status)
                if (status /= linalg_status_success) return
                if (err <= threshold) then
                    lo = candidate+1
                    do while (lo < hi)
                        mid = lo + (hi-lo)/2
                        probe_state = state
                        call refit_real(m,n,mid,p,matvec,rmatvec,tproj,b,status)
                        if (status /= linalg_status_success) return
                        call probe_real(m,n,mid,p,tproj,b,matvec,probe_state,tests,err,status)
                        if (status /= linalg_status_success) return
                        if (err <= threshold) then
                            hi = mid
                        else
                            lo = mid+1
                        end if
                    end do
                    candidate = lo
                    probe_state = state
                    call refit_real(m,n,candidate,p,matvec,rmatvec,tproj,b,status)
                    if (status /= linalg_status_success) return
                    call probe_real(m,n,candidate,p,tproj,b,matvec,probe_state,tests,err,status)
                    if (status /= linalg_status_success) return
                    if (err <= threshold) exit
                end if
            end if
            if (rows == rmax) then
                status = linalg_status_no_convergence
                return
            end if
            rows = rows + min(rows,rmax-rows)
        end do
        k = candidate
        deallocate(idx,proj)
        call move_alloc(p,idx)
        call move_alloc(tproj,proj)
        if (present(skeleton)) then
            deallocate(skeleton)
            call move_alloc(b,skeleton)
        end if
        if (present(estimated_error)) estimated_error = err

    contains
        subroutine empty_result()
            deallocate(idx,proj)
            allocate(idx(n),proj(0,n))
            idx = [(j,j=1,n)]
            if (present(skeleton)) then
                deallocate(skeleton)
                allocate(skeleton(m,0))
            end if
            if (present(estimated_error)) estimated_error = 0.0_dp
            status = linalg_status_success
        end subroutine empty_result
    end subroutine real_tolerance

    subroutine complex_tolerance(m, n, eps, matvec, rmatvec, k, idx, proj, status, &
        seed, oversampling, verification_probes, skeleton, estimated_error)
        integer, intent(in) :: m !! number of rows of implicit A
        integer, intent(in) :: n !! number of columns of implicit A
        real(dp), intent(in) :: eps !! relative Frobenius error target, 0 <= eps < 1
        procedure(id_complex_matvec) :: matvec !! A*x callback
        procedure(id_complex_matvec) :: rmatvec !! A^H*x callback
        integer, intent(out) :: k !! selected numerical rank
        integer, allocatable, intent(out) :: idx(:) !! one-based n-column permutation
        complex(dp), allocatable, intent(out) :: proj(:, :) !! k-by-(n-k) coefficients
        integer, intent(out) :: status !! linalg_status_* outcome
        integer(int64), intent(in), optional :: seed !! reproducible local seed
        integer, intent(in), optional :: oversampling !! initial sketch size hint, default 8
        integer, intent(in), optional :: verification_probes !! Gaussian residual probes, default 16
        complex(dp), allocatable, intent(out), optional :: skeleton(:, :) !! m-by-k exact selected columns
        real(dp), intent(out), optional :: estimated_error !! sampled relative Frobenius residual
        complex(dp), allocatable :: y(:, :), unused(:, :), tproj(:, :), b(:, :)
        integer, allocatable :: p(:)
        integer(int64) :: state, probe_state
        integer :: over, tests, rows, rmax, candidate, lo, hi, mid, j
        real(dp) :: err, floor_tolerance, threshold

        k = 0
        allocate(idx(0),proj(0,0))
        if (present(skeleton)) allocate(skeleton(0,0))
        if (present(estimated_error)) estimated_error = huge(1.0_dp)
        status = linalg_status_invalid_shape
        over = 8
        tests = 16
        if (present(oversampling)) over = oversampling
        if (present(verification_probes)) tests = verification_probes
        if (m < 0 .or. n < 0 .or. over < 0 .or. tests < 1) return
        if (.not. ieee_is_finite(eps) .or. eps < 0.0_dp .or. eps >= 1.0_dp) return
        rmax = min(m,n)
        if (rmax == 0) then
            call empty_result()
            return
        end if
        state = seed_state(seed)
        rows = min(rmax,max(8,min(rmax,max(1,over))))
        floor_tolerance = 100.0_dp*epsilon(1.0_dp)*real(max(m,n),dp)
        threshold = max(min(0.85_dp,max(0.4_dp,1.0_dp-1.5_dp/sqrt(real(tests,dp))))*eps, &
            floor_tolerance)
        do
            call build_sketch_complex(m,n,rows,rmatvec,state,y,status)
            if (status /= linalg_status_success) return
            call interp_decomp(y,eps,candidate,p,unused,status)
            if (status /= linalg_status_success) return
            probe_state = state
            call refit_complex(m,n,candidate,p,matvec,rmatvec,tproj,b,status)
            if (status /= linalg_status_success) return
            call probe_complex(m,n,candidate,p,tproj,b,matvec,probe_state,tests,err,status)
            if (status /= linalg_status_success) return
            if (err <= threshold) exit
            if (candidate < rows) then
                hi = rows
                probe_state = state
                call refit_complex(m,n,hi,p,matvec,rmatvec,tproj,b,status)
                if (status /= linalg_status_success) return
                call probe_complex(m,n,hi,p,tproj,b,matvec,probe_state,tests,err,status)
                if (status /= linalg_status_success) return
                if (err <= threshold) then
                    lo = candidate+1
                    do while (lo < hi)
                        mid = lo + (hi-lo)/2
                        probe_state = state
                        call refit_complex(m,n,mid,p,matvec,rmatvec,tproj,b,status)
                        if (status /= linalg_status_success) return
                        call probe_complex(m,n,mid,p,tproj,b,matvec,probe_state,tests,err,status)
                        if (status /= linalg_status_success) return
                        if (err <= threshold) then
                            hi = mid
                        else
                            lo = mid+1
                        end if
                    end do
                    candidate = lo
                    probe_state = state
                    call refit_complex(m,n,candidate,p,matvec,rmatvec,tproj,b,status)
                    if (status /= linalg_status_success) return
                    call probe_complex(m,n,candidate,p,tproj,b,matvec,probe_state,tests,err,status)
                    if (status /= linalg_status_success) return
                    if (err <= threshold) exit
                end if
            end if
            if (rows == rmax) then
                status = linalg_status_no_convergence
                return
            end if
            rows = rows + min(rows,rmax-rows)
        end do
        k = candidate
        deallocate(idx,proj)
        call move_alloc(p,idx)
        call move_alloc(tproj,proj)
        if (present(skeleton)) then
            deallocate(skeleton)
            call move_alloc(b,skeleton)
        end if
        if (present(estimated_error)) estimated_error = err

    contains
        subroutine empty_result()
            deallocate(idx,proj)
            allocate(idx(n),proj(0,n))
            idx = [(j,j=1,n)]
            if (present(skeleton)) then
                deallocate(skeleton)
                allocate(skeleton(m,0))
            end if
            if (present(estimated_error)) estimated_error = 0.0_dp
            status = linalg_status_success
        end subroutine empty_result
    end subroutine complex_tolerance

end module scifort_linalg_interpolative_operator
