! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Complex generalized Schur factors for modest regular matrix pairs.
! The new 'iterative' method uses an orthogonal Hessenberg--triangular
! reduction and two-sided implicit single-shift QZ iterations, without
! forming inverses. The previous inverse-based method is retained as a
! separate fallback for cases beyond the present iterative implementation.
module scifort_linalg_qz
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_dense, only : lu_factor, lu_solve
    use scifort_linalg_qr, only : qr
    use scifort_linalg_schur, only : schur
    use scifort_linalg_matrix_functions, only : matrix_function_max_order
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: qz

    interface qz
        module procedure qz_rr
        module procedure qz_rc
        module procedure qz_cr
        module procedure qz_cc
    end interface qz

contains

    subroutine qz_rr(a, b, aa, bb, q, z, status, alpha, beta, method)
        real(dp), intent(in) :: a(:, :) !! first real square coefficient matrix, order 0..40
        real(dp), intent(in) :: b(:, :) !! second real square coefficient matrix, may be singular
        complex(dp), intent(out) :: aa(:, :) !! complex upper triangular factor of A
        complex(dp), intent(out) :: bb(:, :) !! complex upper triangular factor of B
        complex(dp), intent(out) :: q(:, :) !! left unitary Schur basis
        complex(dp), intent(out) :: z(:, :) !! right unitary Schur basis
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: alpha(:) !! generalized numerator eigenvalues
        complex(dp), intent(out), optional :: beta(:) !! generalized denominator values

        character(len=*), intent(in), optional :: method !! auto, iterative, inverse

        call qz_core(cmplx(a,kind=dp),cmplx(b,kind=dp),aa,bb,q,z,status,alpha,beta,method)
    end subroutine qz_rr

    subroutine qz_rc(a, b, aa, bb, q, z, status, alpha, beta, method)
        real(dp), intent(in) :: a(:, :) !! first real square coefficient matrix, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second complex square matrix, may be singular
        complex(dp), intent(out) :: aa(:, :) !! complex upper triangular factor of A
        complex(dp), intent(out) :: bb(:, :) !! complex upper triangular factor of B
        complex(dp), intent(out) :: q(:, :) !! left unitary Schur basis
        complex(dp), intent(out) :: z(:, :) !! right unitary Schur basis
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: alpha(:) !! generalized numerator eigenvalues
        complex(dp), intent(out), optional :: beta(:) !! generalized denominator values

        character(len=*), intent(in), optional :: method !! auto, iterative, inverse

        call qz_core(cmplx(a,kind=dp),b,aa,bb,q,z,status,alpha,beta,method)
    end subroutine qz_rc

    subroutine qz_cr(a, b, aa, bb, q, z, status, alpha, beta, method)
        complex(dp), intent(in) :: a(:, :) !! first complex square coefficient matrix, order 0..40
        real(dp), intent(in) :: b(:, :) !! second real square matrix, may be singular
        complex(dp), intent(out) :: aa(:, :) !! complex upper triangular factor of A
        complex(dp), intent(out) :: bb(:, :) !! complex upper triangular factor of B
        complex(dp), intent(out) :: q(:, :) !! left unitary Schur basis
        complex(dp), intent(out) :: z(:, :) !! right unitary Schur basis
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: alpha(:) !! generalized numerator eigenvalues
        complex(dp), intent(out), optional :: beta(:) !! generalized denominator values

        character(len=*), intent(in), optional :: method !! auto, iterative, inverse

        call qz_core(a,cmplx(b,kind=dp),aa,bb,q,z,status,alpha,beta,method)
    end subroutine qz_cr

    subroutine qz_cc(a, b, aa, bb, q, z, status, alpha, beta, method)
        complex(dp), intent(in) :: a(:, :) !! first complex square coefficient matrix, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second complex square coefficient matrix, may be singular
        complex(dp), intent(out) :: aa(:, :) !! complex upper triangular factor of A
        complex(dp), intent(out) :: bb(:, :) !! complex upper triangular factor of B
        complex(dp), intent(out) :: q(:, :) !! left unitary Schur basis
        complex(dp), intent(out) :: z(:, :) !! right unitary Schur basis
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: alpha(:) !! generalized numerator eigenvalues
        complex(dp), intent(out), optional :: beta(:) !! generalized denominator values

        character(len=*), intent(in), optional :: method !! auto, iterative, inverse

        call qz_core(a,b,aa,bb,q,z,status,alpha,beta,method)
    end subroutine qz_cc

    pure logical function finite_matrix(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! complex matrix to check for finite parts
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_matrix


    ! Preserve the original public call sequence; method is appended.
    ! auto uses the two-sided QZ iteration, with the established inverse
    ! construction as fallback on numerical nonconvergence.
    subroutine qz_core(a, b, aa, bb, q, z, status, alpha, beta, method)
        complex(dp), intent(in) :: a(:, :) !! original first matrix
        complex(dp), intent(in) :: b(:, :) !! original second matrix
        complex(dp), intent(out) :: aa(:, :) !! first upper triangular factor
        complex(dp), intent(out) :: bb(:, :) !! second upper triangular factor
        complex(dp), intent(out) :: q(:, :) !! left unitary factor
        complex(dp), intent(out) :: z(:, :) !! right unitary factor
        integer, intent(out) :: status !! linalg status code
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerator
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominator
        character(len=*), intent(in), optional :: method !! auto, iterative, inverse
        character(len=:), allocatable :: selected
        integer :: n

        selected = 'auto'
        if (present(method)) selected = trim(adjustl(method))
        n = size(a,1)
        if (selected /= 'auto' .and. selected /= 'iterative' .and. selected /= 'inverse') then
            status = linalg_status_invalid_shape
            aa = cmplx(0.0_dp,0.0_dp,dp)
            bb = cmplx(0.0_dp,0.0_dp,dp)
            q = cmplx(0.0_dp,0.0_dp,dp)
            z = cmplx(0.0_dp,0.0_dp,dp)
            if (present(alpha)) alpha = cmplx(0.0_dp,0.0_dp,dp)
            if (present(beta)) beta = cmplx(0.0_dp,0.0_dp,dp)
            return
        end if
        if (selected == 'inverse') then
            call qz_inverse_core(a,b,aa,bb,q,z,status,alpha,beta)
            return
        end if
        call qz_iterative_core(a,b,aa,bb,q,z,status,alpha,beta)
        if (selected == 'iterative') return
        if ((status == linalg_status_no_convergence .or. status == linalg_status_singular) .and. &
            n <= matrix_function_max_order) then
            call qz_inverse_core(a,b,aa,bb,q,z,status,alpha,beta)
        end if
    end subroutine qz_core

    ! Construct a real-cosine, complex-sine unitary Givens rotation:
    ! [c,s; -conjg(s),c]*[x;y] = [r;0]. The scaling prevents overflow
    ! in a naive sum of squared magnitudes.
    pure subroutine qz_givens(x, y, c, s)
        complex(dp), intent(in) :: x !! first component to retain
        complex(dp), intent(in) :: y !! second component to eliminate
        real(dp), intent(out) :: c !! cosine in [0,1]
        complex(dp), intent(out) :: s !! complex sine
        real(dp) :: magx, magy, scale, normxy
        magx = abs(x)
        magy = abs(y)
        if (magy == 0.0_dp) then
            c = 1.0_dp
            s = cmplx(0.0_dp,0.0_dp,dp)
        else if (magx == 0.0_dp) then
            c = 0.0_dp
            s = conjg(y) / magy
        else
            scale = max(magx,magy)
            normxy = sqrt((magx/scale)**2 + (magy/scale)**2)
            c = (magx/scale) / normxy
            s = (x/magx) * conjg(y/scale) / normxy
        end if
    end subroutine qz_givens

    ! Left-apply G and accumulate Q <- Q G^H.
    subroutine qz_left(a, b, q, k, c, s)
        complex(dp), intent(inout) :: a(:, :) !! first transformed matrix
        complex(dp), intent(inout) :: b(:, :) !! second transformed matrix
        complex(dp), intent(inout) :: q(:, :) !! accumulated left unitary
        integer, intent(in) :: k !! first affected row, k+1 is second
        real(dp), intent(in) :: c !! real cosine
        complex(dp), intent(in) :: s !! complex sine
        complex(dp) :: x, y
        integer :: j, n
        n = size(a,1)
        do j=1,n
            x = a(k,j)
            y = a(k+1,j)
            a(k,j) = c*x + s*y
            a(k+1,j) = -conjg(s)*x + c*y
            x = b(k,j)
            y = b(k+1,j)
            b(k,j) = c*x + s*y
            b(k+1,j) = -conjg(s)*x + c*y
        end do
        do j=1,n
            x = q(j,k)
            y = q(j,k+1)
            q(j,k) = c*x + conjg(s)*y
            q(j,k+1) = -s*x + c*y
        end do
    end subroutine qz_left

    ! Right-apply V=[c,-conjg(t);t,c] to annihilate B(k+1,k).
    ! Compute the right transformation from [B(k+1,k+1); B(k+1,k)].
    subroutine qz_right(a, b, z, k)
        complex(dp), intent(inout) :: a(:, :) !! first transformed matrix
        complex(dp), intent(inout) :: b(:, :) !! second transformed matrix
        complex(dp), intent(inout) :: z(:, :) !! accumulated right unitary
        integer, intent(in) :: k !! first affected column
        real(dp) :: c
        complex(dp) :: s, t, x, y
        integer :: j, n
        n = size(a,1)
        call qz_givens(b(k+1,k+1),b(k+1,k),c,s)
        t = -conjg(s)
        do j=1,n
            x = a(j,k)
            y = a(j,k+1)
            a(j,k) = c*x + t*y
            a(j,k+1) = -conjg(t)*x + c*y
            x = b(j,k)
            y = b(j,k+1)
            b(j,k) = c*x + t*y
            b(j,k+1) = -conjg(t)*x + c*y
            x = z(j,k)
            y = z(j,k+1)
            z(j,k) = c*x + t*y
            z(j,k+1) = -conjg(t)*x + c*y
        end do
        b(k+1,k) = cmplx(0.0_dp,0.0_dp,dp)
    end subroutine qz_right

    ! A two-sided Hessenberg--triangular reduction, without either inverse.
    subroutine qz_reduce_ht(a, b, q, z, status)
        complex(dp), intent(inout) :: a(:, :) !! first matrix, reduced to Hessenberg
        complex(dp), intent(inout) :: b(:, :) !! second matrix, reduced to triangular
        complex(dp), intent(out) :: q(:, :) !! left unitary
        complex(dp), intent(out) :: z(:, :) !! right unitary
        integer, intent(out) :: status !! linalg status
        complex(dp), allocatable :: r(:, :), old_a(:, :)
        complex(dp) :: s
        real(dp) :: c
        integer :: n, i, j
        n = size(a,1)
        allocate(r(n,n),old_a(n,n))
        old_a = a
        call qr(b,q,r,status)
        if (status /= linalg_status_success) return
        b = r
        a = matmul(transpose(conjg(q)),old_a)
        z = cmplx(0.0_dp,0.0_dp,dp)
        do i=1,n
            z(i,i) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        do j=1,n-2
            do i=n,j+2,-1
                call qz_givens(a(i-1,j),a(i,j),c,s)
                call qz_left(a,b,q,i-1,c,s)
                a(i,j) = cmplx(0.0_dp,0.0_dp,dp)
                call qz_right(a,b,z,i-1)
            end do
        end do
    end subroutine qz_reduce_ht

    ! One implicit single-shift step on active Hessenberg-triangular block
    ! lo:hi. Homogeneous (alpha,beta) avoids dividing by small B diagonals.
    subroutine qz_step(a, b, q, z, lo, hi, alpha, beta)
        complex(dp), intent(inout) :: a(:, :) !! Hessenberg first matrix
        complex(dp), intent(inout) :: b(:, :) !! triangular second matrix
        complex(dp), intent(inout) :: q(:, :) !! accumulated left unitary
        complex(dp), intent(inout) :: z(:, :) !! accumulated right unitary
        integer, intent(in) :: lo !! first active index
        integer, intent(in) :: hi !! last active index
        complex(dp), intent(in) :: alpha !! homogeneous shift numerator
        complex(dp), intent(in) :: beta !! homogeneous shift denominator
        real(dp) :: c
        complex(dp) :: s, x, y
        integer :: k
        x = beta*a(lo,lo) - alpha*b(lo,lo)
        y = beta*a(lo+1,lo)
        call qz_givens(x,y,c,s)
        call qz_left(a,b,q,lo,c,s)
        call qz_right(a,b,z,lo)
        do k=lo+1,hi-1
            call qz_givens(a(k,k-1),a(k+1,k-1),c,s)
            call qz_left(a,b,q,k,c,s)
            a(k+1,k-1) = cmplx(0.0_dp,0.0_dp,dp)
            call qz_right(a,b,z,k)
        end do
    end subroutine qz_step

    ! Wilkinson-style shift from the trailing generalized 2x2 pencil.
    ! For real pencils with a conjugate pair, the complex root is essential:
    ! a real Rayleigh shift cannot triangularize such a block over the reals.
    pure subroutine qz_choose_shift(a,b,hi,alpha,beta)
        complex(dp), intent(in) :: a(:, :) !! Hessenberg first matrix
        complex(dp), intent(in) :: b(:, :) !! upper triangular second matrix
        integer, intent(in) :: hi !! trailing block endpoint, >=2
        complex(dp), intent(out) :: alpha !! normalized homogeneous shift numerator
        complex(dp), intent(out) :: beta !! normalized homogeneous shift denominator
        complex(dp) :: p2, p1, p0, disc, root1, root2, target, shift
        real(dp) :: scale, magnitude
        integer :: k
        k = hi-1
        p2 = b(k,k)*b(hi,hi)
        p1 = -(a(k,k)*b(hi,hi)+a(hi,hi)*b(k,k)-a(hi,k)*b(k,hi))
        p0 = a(k,k)*a(hi,hi)-a(k,hi)*a(hi,k)
        if (abs(p2) > 64.0_dp*epsilon(1.0_dp)* &
            max(abs(p1),abs(p0),tiny(1.0_dp))) then
            disc = sqrt(p1*p1-4.0_dp*p2*p0)
            root1 = (-p1+disc)/(2.0_dp*p2)
            root2 = (-p1-disc)/(2.0_dp*p2)
            if (abs(b(hi,hi)) > 64.0_dp*epsilon(1.0_dp)) then
                target = a(hi,hi)/b(hi,hi)
                if (abs(root1-target) <= abs(root2-target)) then
                    shift = root1
                else
                    shift = root2
                end if
            else
                if (abs(root1) <= abs(root2)) then
                    shift = root1
                else
                    shift = root2
                end if
            end if
            scale = max(1.0_dp,abs(shift))
            alpha = shift/scale
            beta = cmplx(1.0_dp/scale,0.0_dp,dp)
        else if (abs(p1) > 64.0_dp*epsilon(1.0_dp)* &
            max(abs(p0),tiny(1.0_dp))) then
            shift = -p0/p1
            scale = max(1.0_dp,abs(shift))
            alpha = shift/scale
            beta = cmplx(1.0_dp/scale,0.0_dp,dp)
        else
            alpha = a(hi,hi)
            beta = b(hi,hi)
            magnitude = max(abs(alpha),abs(beta))
            if (magnitude > 0.0_dp) then
                alpha = alpha/magnitude
                beta = beta/magnitude
            else
                alpha = cmplx(0.0_dp,0.0_dp,dp)
                beta = cmplx(1.0_dp,0.0_dp,dp)
            end if
        end if
    end subroutine qz_choose_shift

    subroutine qz_iterative_core(a,b,aa,bb,q,z,status,alpha,beta)
        complex(dp), intent(in) :: a(:, :) !! first coefficient matrix
        complex(dp), intent(in) :: b(:, :) !! second coefficient matrix
        complex(dp), intent(out) :: aa(:, :) !! upper triangular first factor
        complex(dp), intent(out) :: bb(:, :) !! upper triangular second factor
        complex(dp), intent(out) :: q(:, :) !! accumulated left unitary
        complex(dp), intent(out) :: z(:, :) !! accumulated right unitary
        integer, intent(out) :: status !! success, invalid shape, no convergence
        complex(dp), intent(out), optional :: alpha(:) !! generalized numerators
        complex(dp), intent(out), optional :: beta(:) !! generalized denominators
        complex(dp), allocatable :: a0(:, :), b0(:, :), work(:, :)
        complex(dp) :: phase, shift_a, shift_b
        integer :: n, i, j, lo, hi, sweep, sweeps_limit
        real(dp) :: sa, sb, err_a, err_b, tolerance, local_threshold, bound

        aa = cmplx(0.0_dp,0.0_dp,dp)
        bb = cmplx(0.0_dp,0.0_dp,dp)
        q = cmplx(0.0_dp,0.0_dp,dp)
        z = cmplx(0.0_dp,0.0_dp,dp)
        if (present(alpha)) alpha = cmplx(0.0_dp,0.0_dp,dp)
        if (present(beta)) beta = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (n > matrix_function_max_order .or. size(a,2) /= n) return
        if (size(b,1) /= n .or. size(b,2) /= n) return
        if (size(aa,1) /= n .or. size(aa,2) /= n) return
        if (size(bb,1) /= n .or. size(bb,2) /= n) return
        if (size(q,1) /= n .or. size(q,2) /= n) return
        if (size(z,1) /= n .or. size(z,2) /= n) return
        if (present(alpha)) then
            if (size(alpha) /= n) return
        end if
        if (present(beta)) then
            if (size(beta) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. finite_matrix(a) .or. .not. finite_matrix(b)) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        sa = max(maxval(abs(real(a,dp))),maxval(abs(aimag(a))))
        sb = max(maxval(abs(real(b,dp))),maxval(abs(aimag(b))))
        if (sa == 0.0_dp) sa = 1.0_dp
        if (sb == 0.0_dp) sb = 1.0_dp
        allocate(a0(n,n),b0(n,n),work(n,n))
        a0 = a/sa
        b0 = b/sb
        aa = a0
        bb = b0
        call qz_reduce_ht(aa,bb,q,z,status)
        if (status /= linalg_status_success) return
        tolerance = 64.0_dp*epsilon(1.0_dp)
        bound = 2.0e-10_dp
        hi = n
        sweep = 0
        sweeps_limit = max(120,200*n)
        do while (hi > 1)
            if (.not. finite_matrix(aa) .or. .not. finite_matrix(bb)) then
                status = linalg_status_no_convergence
                return
            end if
            local_threshold = tolerance*(abs(aa(hi,hi)) + abs(aa(hi-1,hi-1)) + &
                abs(bb(hi,hi)) + abs(bb(hi-1,hi-1)))
            if (abs(aa(hi,hi-1)) <= max(local_threshold,1.0e-300_dp)) then
                aa(hi,hi-1) = cmplx(0.0_dp,0.0_dp,dp)
                hi = hi-1
                cycle
            end if
            lo = 1
            do j=hi-1,2,-1
                local_threshold = tolerance*(abs(aa(j,j)) + abs(aa(j-1,j-1)) + &
                    abs(bb(j,j)) + abs(bb(j-1,j-1)))
                if (abs(aa(j,j-1)) <= max(local_threshold,1.0e-300_dp)) then
                    aa(j,j-1) = cmplx(0.0_dp,0.0_dp,dp)
                    lo = j
                    exit
                end if
            end do
            call qz_choose_shift(aa,bb,hi,shift_a,shift_b)
            if (sweep > sweeps_limit) then
                status = linalg_status_no_convergence
                return
            end if
            call qz_step(aa,bb,q,z,lo,hi,shift_a,shift_b)
            sweep = sweep+1
        end do
        do j=1,n-1
            aa(j+1:n,j) = cmplx(0.0_dp,0.0_dp,dp)
            bb(j+1:n,j) = cmplx(0.0_dp,0.0_dp,dp)
        end do
        ! Rotate rows by diagonal phases, preserving the two identities.
        do i=1,n
            if (abs(bb(i,i)) > 0.0_dp) then
                phase = bb(i,i)/abs(bb(i,i))
            else if (abs(aa(i,i)) > 0.0_dp) then
                phase = aa(i,i)/abs(aa(i,i))
            else
                phase = cmplx(1.0_dp,0.0_dp,dp)
            end if
            q(:,i) = phase*q(:,i)
            aa(i,:) = conjg(phase)*aa(i,:)
            bb(i,:) = conjg(phase)*bb(i,:)
        end do
        work = matmul(q,matmul(aa,transpose(conjg(z))))
        err_a = maxval(abs(work-a0))/max(1.0_dp,maxval(abs(a0)))
        work = matmul(q,matmul(bb,transpose(conjg(z))))
        err_b = maxval(abs(work-b0))/max(1.0_dp,maxval(abs(b0)))
        if (.not. ieee_is_finite(err_a) .or. .not. ieee_is_finite(err_b)) then
            status = linalg_status_no_convergence
            return
        end if
        if (max(err_a,err_b) > bound) then
            status = linalg_status_no_convergence
            return
        end if
        ! A simultaneous zero of the two triangular diagonals means that
        ! the pencil is singular, rather than an infinite eigenvalue.
        do i=1,n
            if (abs(aa(i,i)) <= tolerance .and. abs(bb(i,i)) <= tolerance) then
                status = linalg_status_singular
                return
            end if
        end do
        aa = aa*sa
        bb = bb*sb
        if (.not. finite_matrix(aa) .or. .not. finite_matrix(bb)) then
            status = linalg_status_no_convergence
            return
        end if
        if (present(alpha)) then
            do i=1,n
                alpha(i) = aa(i,i)
            end do
        end if
        if (present(beta)) then
            do i=1,n
                beta(i) = bb(i,i)
            end do
        end if
        status = linalg_status_success
    end subroutine qz_iterative_core

    subroutine qz_inverse_core(a, b, aa, bb, q, z, status, alpha, beta)
        complex(dp), intent(in) :: a(:, :) !! first square matrix, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second square matrix; may be singular
        complex(dp), intent(out) :: aa(:, :) !! Schur factor satisfying A = Q AA Z**H
        complex(dp), intent(out) :: bb(:, :) !! Schur factor satisfying B = Q BB Z**H
        complex(dp), intent(out) :: q(:, :) !! unitary left matrix, order n
        complex(dp), intent(out) :: z(:, :) !! unitary right matrix, order n
        integer, intent(out) :: status !! success, invalid shape, singular B, or no convergence
        complex(dp), intent(out), optional :: alpha(:) !! diagonal of AA, size n
        complex(dp), intent(out), optional :: beta(:) !! diagonal of BB, size n

        complex(dp), allocatable :: a0(:, :), b0(:, :), lu(:, :), c(:, :)
        complex(dp), allocatable :: u(:, :), r(:, :), product(:, :), work(:, :)
        complex(dp) :: phase, shift
        integer, allocatable :: ipiv(:)
        integer :: n, i, j, route, shift_index
        real(dp) :: scale_a, scale_b, residual_a, residual_b, solve_error, bound
        complex(dp), parameter :: shifts(6) = [ &
            cmplx(1.0_dp,0.0_dp,dp),cmplx(-1.0_dp,0.0_dp,dp), &
            cmplx(0.0_dp,1.0_dp,dp),cmplx(0.0_dp,-1.0_dp,dp), &
            cmplx(2.0_dp,0.0_dp,dp),cmplx(0.5_dp,0.0_dp,dp)]

        aa = cmplx(0.0_dp, 0.0_dp, dp)
        bb = cmplx(0.0_dp, 0.0_dp, dp)
        q = cmplx(0.0_dp, 0.0_dp, dp)
        z = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(alpha)) alpha = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(beta)) beta = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (n > matrix_function_max_order .or. size(a,2) /= n) return
        if (size(b,1) /= n .or. size(b,2) /= n) return
        if (size(aa,1) /= n .or. size(aa,2) /= n) return
        if (size(bb,1) /= n .or. size(bb,2) /= n) return
        if (size(q,1) /= n .or. size(q,2) /= n) return
        if (size(z,1) /= n .or. size(z,2) /= n) return
        if (present(alpha)) then
            if (size(alpha) /= n) return
        end if
        if (present(beta)) then
            if (size(beta) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. finite_matrix(a) .or. .not. finite_matrix(b)) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        scale_a = max(maxval(abs(real(a,dp))),maxval(abs(aimag(a))))
        scale_b = max(maxval(abs(real(b,dp))),maxval(abs(aimag(b))))
        if (scale_b <= 0.0_dp) scale_b = 1.0_dp
        if (scale_a <= 0.0_dp) scale_a = 1.0_dp
        allocate(a0(n,n),b0(n,n),lu(n,n),c(n,n),u(n,n),r(n,n))
        allocate(product(n,n),work(n,n),ipiv(n))
        a0 = a / scale_a
        b0 = b / scale_b

        ! Route 1: B nonsingular, C = B**(-1) A.
        ! Route 2: A nonsingular, C = A**(-1) B (allows beta=0).
        ! Route 3: both singular, but a shifted pencil K = A - shift*B
        !          nonsingular; C = K**(-1) B. This handles many regular
        !          pencils with finite and infinite generalized eigenvalues.
        route = 1
        shift = cmplx(0.0_dp,0.0_dp,dp)
        call lu_factor(b0,lu,ipiv,status)
        if (status /= linalg_status_success) then
            route = 2
            call lu_factor(a0,lu,ipiv,status)
        end if
        if (status /= linalg_status_success) then
            route = 3
            do shift_index=1,size(shifts)
                shift = shifts(shift_index)
                work = a0 - shift*b0
                call lu_factor(work,lu,ipiv,status)
                if (status == linalg_status_success) exit
            end do
            if (status /= linalg_status_success) then
                status = linalg_status_singular
                return
            end if
        end if
        if (route == 1) then
            call lu_solve(lu,ipiv,a0,c,status)
        else
            call lu_solve(lu,ipiv,b0,c,status)
        end if
        if (status /= linalg_status_success) return
        if (.not. finite_matrix(c)) then
            status = linalg_status_no_convergence
            return
        end if
        select case(route)
        case(1)
            solve_error = maxval(abs(matmul(b0,c)-a0)) / &
                max(1.0_dp,maxval(abs(a0)))
        case(2)
            solve_error = maxval(abs(matmul(a0,c)-b0)) / &
                max(1.0_dp,maxval(abs(b0)))
        case(3)
            solve_error = maxval(abs(matmul(a0-shift*b0,c)-b0)) / &
                max(1.0_dp,maxval(abs(b0)))
        end select
        bound = 1.0e-9_dp
        if (.not. ieee_is_finite(solve_error) .or. solve_error > bound) then
            status = linalg_status_no_convergence
            return
        end if
        call schur(c,u,z,status)
        if (status /= linalg_status_success) return
        select case(route)
        case(1)
            product = matmul(b0,z)
        case(2)
            product = matmul(a0,z)
        case(3)
            product = matmul(a0-shift*b0,z)
        end select
        call qr(product,q,r,status)
        if (status /= linalg_status_success) return
        select case(route)
        case(1)
            aa = matmul(r,u)
            bb = r
        case(2)
            aa = r
            bb = matmul(r,u)
        case(3)
            bb = matmul(r,u)
            aa = r + shift*bb
        end select

        ! Choose column phases so the generalized Schur B diagonal is real, nonnegative.
        do i = 1, n
            if (abs(bb(i,i)) > 0.0_dp) then
                phase = bb(i,i)/abs(bb(i,i))
            else if (abs(aa(i,i)) > 0.0_dp) then
                phase = aa(i,i)/abs(aa(i,i))
            else
                phase = cmplx(1.0_dp,0.0_dp,dp)
            end if
            q(:,i) = phase*q(:,i)
            aa(i,:) = conjg(phase)*aa(i,:)
            bb(i,:) = conjg(phase)*bb(i,:)
        end do
        do j = 1, n-1
            aa(j+1:n,j) = cmplx(0.0_dp,0.0_dp,dp)
            bb(j+1:n,j) = cmplx(0.0_dp,0.0_dp,dp)
        end do

        ! Verify both original *scaled* equations before restoring physical units.
        work = matmul(q,matmul(aa,transpose(conjg(z))))
        residual_a = maxval(abs(work-a0))/max(1.0_dp,maxval(abs(a0)))
        work = matmul(q,matmul(bb,transpose(conjg(z))))
        residual_b = maxval(abs(work-b0))/max(1.0_dp,maxval(abs(b0)))
        if (.not. ieee_is_finite(residual_a) .or. &
            .not. ieee_is_finite(residual_b)) then
            status = linalg_status_no_convergence
            return
        end if
        if (max(residual_a,residual_b) > bound) then
            status = linalg_status_no_convergence
            return
        end if

        aa = aa*scale_a
        bb = bb*scale_b
        if (.not. finite_matrix(aa) .or. .not. finite_matrix(bb)) then
            status = linalg_status_no_convergence
            return
        end if
        if (present(alpha)) then
            do i = 1, n
                alpha(i) = aa(i,i)
            end do
        end if
        if (present(beta)) then
            do i = 1, n
                beta(i) = bb(i,i)
            end do
        end if
        status = linalg_status_success
    end subroutine qz_inverse_core

end module scifort_linalg_qz
