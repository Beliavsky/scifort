! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Complete cosine-sine decomposition of a partitioned real orthogonal or
! complex unitary matrix.  Uses SciFort's own one-sided Jacobi SVD, not LAPACK.
! The CS angles are intrinsic; factor columns can differ by phases/signs from
! LAPACK's ORCSD/UNCSD gauge.  The middle matrix matches SciPy's layout.
module scifort_linalg_cossin
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_svd, only : svd
    implicit none
    private

    public :: cossin, cossin_separate

    interface cossin
        module procedure cossin_real
        module procedure cossin_complex
    end interface

    interface cossin_separate
        module procedure cossin_separate_real
        module procedure cossin_separate_complex
    end interface

contains

    pure logical function finite_complex(x) result(ok)
        complex(dp), intent(in) :: x(:, :) !! input matrix
        ok = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function finite_complex

    pure logical function output_shape(a, n, enabled) result(ok)
        complex(dp), intent(in) :: a(:, :) !! candidate output
        integer, intent(in) :: n !! square order when enabled
        logical, intent(in) :: enabled !! whether this factor was requested
        ok = size(a,1) == n .and. size(a,2) == n
        if (.not. enabled) ok = ok .or. (size(a,1) == 0 .and. size(a,2) == 0)
    end function output_shape

    pure real(dp) function vector_norm(x) result(norm)
        complex(dp), intent(in) :: x(:) !! orthonormal-basis candidate
        norm = sqrt(sum(abs(x)**2))
    end function vector_norm

    pure logical function approximately_unitary(x) result(ok)
        complex(dp), intent(in) :: x(:, :) !! m-by-m input to verify
        complex(dp), allocatable :: gram(:, :)
        real(dp) :: tolerance
        integer :: j, n

        ok = .false.
        n = size(x,1)
        if (n < 2 .or. size(x,2) /= n) return
        if (.not. finite_complex(x)) return
        gram = matmul(transpose(conjg(x)),x)
        do j = 1, n
            gram(j,j) = gram(j,j) - cmplx(1.0_dp,0.0_dp,dp)
        end do
        tolerance = max(1.0e-10_dp, 256.0_dp*epsilon(1.0_dp)*real(n,dp))
        ok = maxval(abs(gram)) <= tolerance
    end function approximately_unitary

    ! Fill predetermined columns and complete a unitary square matrix.
    ! Known columns are renormalized by two-pass modified Gram-Schmidt.
    ! Remaining columns are chosen greedily from the coordinate basis.
    subroutine complete_basis(a, present_column, status)
        complex(dp), intent(inout) :: a(:, :) !! n-by-n matrix, known columns on entry
        logical, intent(inout) :: present_column(:) !! true for known columns
        integer, intent(out) :: status !! success or no_convergence
        complex(dp), allocatable :: trial(:), best_trial(:), basis(:)
        real(dp) :: best_norm, trial_norm
        integer :: n, j, k, i, pass

        status = linalg_status_no_convergence
        n = size(a,1)
        if (n < 1 .or. size(a,2) /= n .or. size(present_column) /= n) return
        allocate(trial(n),best_trial(n),basis(n))
        do j = 1, n
            if (present_column(j)) then
                trial = a(:,j)
            else
                best_norm = -1.0_dp
                do i = 1, n
                    basis = cmplx(0.0_dp,0.0_dp,dp)
                    basis(i) = cmplx(1.0_dp,0.0_dp,dp)
                    do pass = 1, 2
                        do k = 1, n
                            if (present_column(k)) basis = basis - dot_product(a(:,k),basis)*a(:,k)
                        end do
                    end do
                    trial_norm = vector_norm(basis)
                    if (trial_norm > best_norm) then
                        best_norm = trial_norm
                        best_trial = basis
                    end if
                end do
                trial = best_trial
            end if
            do pass = 1, 2
                do k = 1, n
                    if (present_column(k) .and. k /= j) then
                        trial = trial - dot_product(a(:,k),trial)*a(:,k)
                    end if
                end do
            end do
            trial_norm = vector_norm(trial)
            if (trial_norm < 128.0_dp*epsilon(1.0_dp)) return
            a(:,j) = trial/trial_norm
            present_column(j) = .true.
        end do
        status = linalg_status_success
    end subroutine complete_basis

    subroutine cossin_separate_complex(x, p, q, u1, u2, theta, v1h, v2h, status, swap_sign)
        complex(dp), intent(in) :: x(:, :) !! m-by-m approximately unitary matrix
        integer, intent(in) :: p !! upper block row count, 0 < p < m
        integer, intent(in) :: q !! left block column count, 0 < q < m
        complex(dp), intent(out) :: u1(:, :) !! upper p-by-p unitary factor
        complex(dp), intent(out) :: u2(:, :) !! lower (m-p)-by-(m-p) unitary factor
        real(dp), intent(out) :: theta(:) !! min(p,m-p,q,m-q) CS angles in [0,pi/2]
        complex(dp), intent(out) :: v1h(:, :) !! q-by-q conjugate right factor
        complex(dp), intent(out) :: v2h(:, :) !! (m-q)-by-(m-q) conjugate right factor
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: swap_sign !! exchange signs of sine / off-diagonal identity blocks
        complex(dp), allocatable :: v1(:, :), v2(:, :), u2_work(:, :), z(:)
        real(dp), allocatable :: singular(:)
        logical, allocatable :: u2_known(:), v2_known(:)
        integer :: m, r, n11, n12, n21, n22, i, j
        real(dp) :: sin_angle, cos_angle
        logical :: flip

        u1 = cmplx(0.0_dp,0.0_dp,dp)
        u2 = cmplx(0.0_dp,0.0_dp,dp)
        v1h = cmplx(0.0_dp,0.0_dp,dp)
        v2h = cmplx(0.0_dp,0.0_dp,dp)
        theta = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(x,1)
        if (m < 2 .or. size(x,2) /= m) return
        if (p < 1 .or. p >= m .or. q < 1 .or. q >= m) return
        r = min(p,m-p,q,m-q)
        if (size(u1,1) /= p .or. size(u1,2) /= p) return
        if (size(u2,1) /= m-p .or. size(u2,2) /= m-p) return
        if (size(v1h,1) /= q .or. size(v1h,2) /= q) return
        if (size(v2h,1) /= m-q .or. size(v2h,2) /= m-q) return
        if (size(theta) /= r) return
        status = linalg_status_no_convergence
        if (.not. approximately_unitary(x)) return

        flip = .false.
        if (present(swap_sign)) flip = swap_sign
        n11 = min(p,q)-r
        n12 = min(p,m-q)-r
        n21 = min(m-p,q)-r
        n22 = min(m-p,m-q)-r
        allocate(singular(min(p,q)),v1(q,q),v2(m-q,m-q),u2_work(m-p,m-p))
        allocate(u2_known(m-p),v2_known(m-q),z(m-p))

        ! X11 = U1*diag(singular)*V1H.  The first n11 modes correspond
        ! to identity entries of CS, followed by r cosine-sine modes.
        call svd(x(1:p,1:q),u1,singular,v1h,status,full_matrices=.true.)
        if (status /= linalg_status_success) go to 900
        v1 = transpose(conjg(v1h))
        u2_work = cmplx(0.0_dp,0.0_dp,dp)
        u2_known = .false.
        do i = 1, r
            j = n11+i
            z = matmul(x(p+1:m,1:q),v1(:,j))
            sin_angle = vector_norm(z)
            cos_angle = vector_norm(matmul(x(1:p,1:q),v1(:,j)))
            theta(i) = atan2(sin_angle,cos_angle)
            if (sin_angle > 64.0_dp*epsilon(1.0_dp)) then
                u2_work(:,n22+i) = z/sin_angle
                if (flip) u2_work(:,n22+i) = -u2_work(:,n22+i)
                u2_known(n22+i) = .true.
            end if
        end do
        do i = 1, n21
            j = n11+r+i
            z = matmul(x(p+1:m,1:q),v1(:,j))
            sin_angle = vector_norm(z)
            if (sin_angle < 64.0_dp*epsilon(1.0_dp)) then
                status = linalg_status_no_convergence
                go to 900
            end if
            u2_work(:,n22+r+i) = z/sin_angle
            if (flip) u2_work(:,n22+r+i) = -u2_work(:,n22+r+i)
            u2_known(n22+r+i) = .true.
        end do
        call complete_basis(u2_work,u2_known,status)
        if (status /= linalg_status_success) go to 900

        ! Determine V2 from X12 and the complemented bottom basis U2.
        ! At theta=0, the sine block vanishes and the lower-right block
        ! supplies the otherwise undetermined right vector instead.
        v2 = cmplx(0.0_dp,0.0_dp,dp)
        v2_known = .true.
        do i = 1, n22
            v2(:,i) = matmul(transpose(conjg(x(p+1:m,q+1:m))),u2_work(:,i))
        end do
        do i = 1, r
            j = n11+i
            sin_angle = sin(theta(i))
            if (sin_angle > 64.0_dp*epsilon(1.0_dp)) then
                v2(:,n22+i) = matmul(transpose(conjg(x(1:p,q+1:m))),u1(:,j))/sin_angle
                if (.not. flip) v2(:,n22+i) = -v2(:,n22+i)
            else
                v2(:,n22+i) = matmul(transpose(conjg(x(p+1:m,q+1:m))),u2_work(:,n22+i))
            end if
        end do
        do i = 1, n12
            j = n11+r+i
            v2(:,n22+r+i) = matmul(transpose(conjg(x(1:p,q+1:m))),u1(:,j))
            if (.not. flip) v2(:,n22+r+i) = -v2(:,n22+r+i)
        end do
        call complete_basis(v2,v2_known,status)
        if (status /= linalg_status_success) go to 900
        u2 = u2_work
        v2h = transpose(conjg(v2))
        return
900     continue
        u1 = cmplx(0.0_dp,0.0_dp,dp)
        u2 = cmplx(0.0_dp,0.0_dp,dp)
        v1h = cmplx(0.0_dp,0.0_dp,dp)
        v2h = cmplx(0.0_dp,0.0_dp,dp)
        theta = 0.0_dp
    end subroutine cossin_separate_complex

    subroutine cossin_separate_real(x, p, q, u1, u2, theta, v1t, v2t, status, swap_sign)
        real(dp), intent(in) :: x(:, :) !! m-by-m approximately orthogonal matrix
        integer, intent(in) :: p !! upper block row count
        integer, intent(in) :: q !! left block column count
        real(dp), intent(out) :: u1(:, :) !! p-by-p orthogonal factor
        real(dp), intent(out) :: u2(:, :) !! (m-p)-by-(m-p) orthogonal factor
        real(dp), intent(out) :: theta(:) !! CS angles in radians
        real(dp), intent(out) :: v1t(:, :) !! q-by-q transposed right factor
        real(dp), intent(out) :: v2t(:, :) !! (m-q)-by-(m-q) transposed right factor
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: swap_sign !! exchange sine block signs
        complex(dp), allocatable :: u1c(:, :), u2c(:, :), v1c(:, :), v2c(:, :)
        integer :: m

        u1 = 0.0_dp
        u2 = 0.0_dp
        theta = 0.0_dp
        v1t = 0.0_dp
        v2t = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(x,1)
        if (m < 2 .or. size(x,2) /= m) return
        if (p < 1 .or. p >= m .or. q < 1 .or. q >= m) return
        allocate(u1c(p,p),u2c(m-p,m-p),v1c(q,q),v2c(m-q,m-q))
        call cossin_separate_complex(cmplx(x,0.0_dp,dp),p,q,u1c,u2c,theta,v1c,v2c,status,swap_sign)
        if (status /= linalg_status_success) return
        if (maxval(abs(aimag(u1c))) > 1.0e-10_dp .or. &
            maxval(abs(aimag(u2c))) > 1.0e-10_dp .or. &
            maxval(abs(aimag(v1c))) > 1.0e-10_dp .or. &
            maxval(abs(aimag(v2c))) > 1.0e-10_dp) then
            status = linalg_status_no_convergence
            theta = 0.0_dp
            return
        end if
        u1 = real(u1c,dp)
        u2 = real(u2c,dp)
        v1t = real(v1c,dp)
        v2t = real(v2c,dp)
    end subroutine cossin_separate_real

    ! Build SciPy's canonical CS middle matrix, including the case p/=q.
    pure subroutine middle_factor(m, p, q, theta, cs, flip)
        integer, intent(in) :: m !! full matrix order
        integer, intent(in) :: p !! top block row count
        integer, intent(in) :: q !! left block column count
        real(dp), intent(in) :: theta(:) !! r CS angles
        real(dp), intent(out) :: cs(:, :) !! m-by-m canonical CS factor
        logical, intent(in) :: flip !! swap signs on off-diagonal terms
        integer :: n11, n12, n21, n22, r, i, sgn
        real(dp) :: c, s

        cs = 0.0_dp
        r = size(theta)
        n11 = min(p,q)-r
        n12 = min(p,m-q)-r
        n21 = min(m-p,q)-r
        n22 = min(m-p,m-q)-r
        sgn = 1
        if (flip) sgn = -1
        do i = 1, n11
            cs(i,i) = 1.0_dp
        end do
        do i = 1, n12
            cs(n11+r+i,q+n22+r+i) = -real(sgn,dp)
        end do
        do i = 1, n21
            cs(p+n22+r+i,n11+r+i) = real(sgn,dp)
        end do
        do i = 1, n22
            cs(p+i,q+i) = 1.0_dp
        end do
        do i = 1, r
            c = cos(theta(i))
            s = sin(theta(i))*real(sgn,dp)
            cs(n11+i,n11+i) = c
            cs(p+n22+i,q+n22+i) = c
            cs(n11+i,q+n22+i) = -s
            cs(p+n22+i,n11+i) = s
        end do
    end subroutine middle_factor

    subroutine cossin_complex(x, p, q, u, cs, vh, status, swap_sign, compute_u, compute_vh)
        complex(dp), intent(in) :: x(:, :) !! m-by-m unitary input
        integer, intent(in) :: p !! partition row count, 0 < p < m
        integer, intent(in) :: q !! partition column count, 0 < q < m
        complex(dp), intent(out) :: u(:, :) !! m-by-m block-unitary, or empty when disabled
        real(dp), intent(out) :: cs(:, :) !! m-by-m real CS middle factor
        complex(dp), intent(out) :: vh(:, :) !! m-by-m block-unitary, or empty when disabled
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: swap_sign !! use alternative CS signs; default false
        logical, intent(in), optional :: compute_u !! return U (default true)
        logical, intent(in), optional :: compute_vh !! return V^H (default true)
        complex(dp), allocatable :: u1(:, :),u2(:, :),v1h(:, :),v2h(:, :)
        real(dp), allocatable :: theta(:)
        logical :: flip, want_u, want_vh
        integer :: m,r

        u = cmplx(0.0_dp,0.0_dp,dp)
        cs = 0.0_dp
        vh = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(x,1)
        if (m < 2 .or. size(x,2) /= m) return
        if (p < 1 .or. p >= m .or. q < 1 .or. q >= m) return
        if (size(cs,1) /= m .or. size(cs,2) /= m) return
        want_u = .true.
        want_vh = .true.
        flip = .false.
        if (present(compute_u)) want_u = compute_u
        if (present(compute_vh)) want_vh = compute_vh
        if (present(swap_sign)) flip = swap_sign
        if (.not. output_shape(u,m,want_u)) return
        if (.not. output_shape(vh,m,want_vh)) return
        r = min(p,m-p,q,m-q)
        allocate(u1(p,p),u2(m-p,m-p),v1h(q,q),v2h(m-q,m-q),theta(r))
        call cossin_separate_complex(x,p,q,u1,u2,theta,v1h,v2h,status,swap_sign=flip)
        if (status /= linalg_status_success) return
        call middle_factor(m,p,q,theta,cs,flip)
        if (want_u) then
            u(1:p,1:p) = u1
            u(p+1:m,p+1:m) = u2
        end if
        if (want_vh) then
            vh(1:q,1:q) = v1h
            vh(q+1:m,q+1:m) = v2h
        end if
    end subroutine cossin_complex

    subroutine cossin_real(x, p, q, u, cs, vt, status, swap_sign, compute_u, compute_vh)
        real(dp), intent(in) :: x(:, :) !! m-by-m orthogonal input
        integer, intent(in) :: p !! partition row count
        integer, intent(in) :: q !! partition column count
        real(dp), intent(out) :: u(:, :) !! m-by-m block-orthogonal, or empty when disabled
        real(dp), intent(out) :: cs(:, :) !! m-by-m real CS middle factor
        real(dp), intent(out) :: vt(:, :) !! m-by-m block-orthogonal, or empty when disabled
        integer, intent(out) :: status !! linalg_status_* outcome
        logical, intent(in), optional :: swap_sign !! use alternative CS signs
        logical, intent(in), optional :: compute_u !! return U (default true)
        logical, intent(in), optional :: compute_vh !! return V^T (default true)
        complex(dp), allocatable :: uc(:, :), vc(:, :)
        integer :: m

        u = 0.0_dp
        cs = 0.0_dp
        vt = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(x,1)
        if (m < 2 .or. size(x,2) /= m) return
        if (p < 1 .or. p >= m .or. q < 1 .or. q >= m) return
        allocate(uc(size(u,1),size(u,2)),vc(size(vt,1),size(vt,2)))
        call cossin_complex(cmplx(x,0.0_dp,dp),p,q,uc,cs,vc,status,swap_sign,compute_u,compute_vh)
        if (status /= linalg_status_success) return
        if (size(uc) > 0) then
            if (maxval(abs(aimag(uc))) > 1.0e-10_dp) then
                status = linalg_status_no_convergence
                cs = 0.0_dp
                return
            end if
        end if
        if (size(vc) > 0) then
            if (maxval(abs(aimag(vc))) > 1.0e-10_dp) then
                status = linalg_status_no_convergence
                cs = 0.0_dp
                return
            end if
        end if
        u = real(uc,dp)
        vt = real(vc,dp)
    end subroutine cossin_real

end module scifort_linalg_cossin
