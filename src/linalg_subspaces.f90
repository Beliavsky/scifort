! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Bases of matrix ranges and null spaces, and principal subspace angles.
! Uses SciFort's economy Jacobi SVD, with two-pass orthogonal completion
! of right singular vectors when an input matrix is wide.
module scifort_linalg_subspaces
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_svd, only : svd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: orth, null_space, subspace_angles

    interface orth
        module procedure orth_real
        module procedure orth_complex
    end interface

    interface null_space
        module procedure null_space_real
        module procedure null_space_complex
    end interface

    interface subspace_angles
        module procedure subspace_angles_real
        module procedure subspace_angles_complex
        module procedure subspace_angles_real_complex
        module procedure subspace_angles_complex_real
    end interface

contains

    pure logical function valid_options(m, n, rcond, max_sweeps) result(ok)
        integer, intent(in) :: m !! row count, nonnegative
        integer, intent(in) :: n !! column count, nonnegative
        real(dp), intent(in), optional :: rcond !! relative singular value cutoff, >= 0
        integer, intent(in), optional :: max_sweeps !! Jacobi sweep limit, >= 1
        ok = .true.
        if (m < 0 .or. n < 0) ok = .false.
        if (present(rcond)) then
            if (.not. ieee_is_finite(rcond)) ok = .false.
            if (rcond < 0.0_dp) ok = .false.
        end if
        if (present(max_sweeps)) then
            if (max_sweeps < 1) ok = .false.
        end if
    end function valid_options

    pure logical function finite_complex_matrix(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! matrix to check for finite components
        ok = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex_matrix

    subroutine orth_real(a, q, status, rcond, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! m-by-n real matrix, including empty matrices
        real(dp), allocatable, intent(out) :: q(:, :) !! m-by-r orthonormal range basis
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: rcond !! relative singular value cutoff, >= 0
        integer, intent(in), optional :: max_sweeps !! positive Jacobi sweep limit
        complex(dp), allocatable :: qc(:, :)

        call orth_complex(cmplx(a, 0.0_dp, dp), qc, status, rcond, max_sweeps)
        allocate(q(size(qc,1),size(qc,2)))
        q = 0.0_dp
        if (status /= linalg_status_success) return
        if (size(qc) > 0) then
            if (maxval(abs(aimag(qc))) > 1.0e-11_dp) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        q = real(qc, dp)
    end subroutine orth_real

    subroutine orth_complex(a, q, status, rcond, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! m-by-n finite complex matrix
        complex(dp), allocatable, intent(out) :: q(:, :) !! m-by-r orthonormal range basis
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: rcond !! relative singular value cutoff, >= 0
        integer, intent(in), optional :: max_sweeps !! positive Jacobi sweep limit
        complex(dp), allocatable :: u(:, :), vh(:, :)
        real(dp), allocatable :: s(:)
        real(dp) :: cutoff
        integer :: m, n, k, r

        m = size(a,1)
        n = size(a,2)
        allocate(q(m,0))
        status = linalg_status_invalid_shape
        if (.not. valid_options(m,n,rcond,max_sweeps)) return
        status = linalg_status_no_convergence
        if (.not. finite_complex_matrix(a)) return
        status = linalg_status_success
        if (min(m,n) == 0) return
        k = min(m,n)
        allocate(u(m,k),s(k),vh(k,n))
        call svd(a,u,s,vh,status,max_sweeps)
        if (status /= linalg_status_success) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        r = count(s > cutoff*s(1))
        deallocate(q)
        allocate(q(m,r))
        q = u(:,1:r)
    end subroutine orth_complex

    subroutine null_space_real(a, z, status, rcond, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! m-by-n real matrix, including empty matrices
        real(dp), allocatable, intent(out) :: z(:, :) !! n-by-(n-r) orthonormal kernel basis
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: rcond !! relative singular value cutoff, >= 0
        integer, intent(in), optional :: max_sweeps !! positive Jacobi sweep limit
        complex(dp), allocatable :: zc(:, :)

        call null_space_complex(cmplx(a, 0.0_dp, dp), zc, status, rcond, max_sweeps)
        allocate(z(size(zc,1),size(zc,2)))
        z = 0.0_dp
        if (status /= linalg_status_success) return
        if (size(zc) > 0) then
            if (maxval(abs(aimag(zc))) > 1.0e-11_dp) then
                status = linalg_status_no_convergence
                return
            end if
        end if
        z = real(zc,dp)
    end subroutine null_space_real

    subroutine null_space_complex(a, z, status, rcond, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! m-by-n finite complex matrix
        complex(dp), allocatable, intent(out) :: z(:, :) !! n-by-(n-r) orthonormal kernel basis
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: rcond !! relative singular value cutoff, >= 0
        integer, intent(in), optional :: max_sweeps !! positive Jacobi sweep limit
        complex(dp), allocatable :: u(:, :), vh(:, :), v(:, :)
        complex(dp), allocatable :: trial(:)
        real(dp), allocatable :: s(:), residual_diagonal(:)
        real(dp) :: cutoff, trial_norm
        integer :: m, n, k, r, j, pivot, pass, i

        m = size(a,1)
        n = size(a,2)
        allocate(z(n,0))
        status = linalg_status_invalid_shape
        if (.not. valid_options(m,n,rcond,max_sweeps)) return
        status = linalg_status_no_convergence
        if (.not. finite_complex_matrix(a)) return
        status = linalg_status_success
        if (n == 0) return
        if (m == 0) then
            deallocate(z)
            allocate(z(n,n))
            z = cmplx(0.0_dp,0.0_dp,dp)
            do j = 1, n
                z(j,j) = cmplx(1.0_dp,0.0_dp,dp)
            end do
            return
        end if

        k = min(m,n)
        allocate(u(m,k),s(k),vh(k,n))
        call svd(a,u,s,vh,status,max_sweeps)
        if (status /= linalg_status_success) return
        cutoff = real(max(m,n),dp)*epsilon(1.0_dp)
        if (present(rcond)) cutoff = rcond
        r = count(s > cutoff*s(1))
        if (r == 0) then
            deallocate(z)
            allocate(z(n,n))
            z = cmplx(0.0_dp,0.0_dp,dp)
            do j = 1, n
                z(j,j) = cmplx(1.0_dp,0.0_dp,dp)
            end do
            return
        end if

        ! The economy SVD contains only k right singular directions. Complete
        ! all k, rather than only the r retained directions, so the final
        ! n-r columns include any discarded singular vectors as well.
        allocate(v(n,n),trial(n),residual_diagonal(n))
        v = cmplx(0.0_dp,0.0_dp,dp)
        v(:,1:k) = transpose(conjg(vh))
        ! Diagonal of I - V_k V_k^H, updated as vectors are completed.
        ! Its largest entry identifies the canonical basis vector with the
        ! largest remaining orthogonal component, in O(n^3) total work.
        residual_diagonal = 1.0_dp
        do i = 1, k
            residual_diagonal = residual_diagonal-abs(v(:,i))**2
        end do
        residual_diagonal = max(0.0_dp,residual_diagonal)
        do j = k+1, n
            pivot = maxloc(residual_diagonal,dim=1)
            trial = cmplx(0.0_dp,0.0_dp,dp)
            trial(pivot) = cmplx(1.0_dp,0.0_dp,dp)
            do pass = 1, 2
                do i = 1, j-1
                    trial = trial-dot_product(v(:,i),trial)*v(:,i)
                end do
            end do
            trial_norm = sqrt(sum(abs(trial)**2))
            if (trial_norm <= 32.0_dp*epsilon(1.0_dp)) then
                status = linalg_status_no_convergence
                return
            end if
            v(:,j) = trial/trial_norm
            residual_diagonal = max(0.0_dp,residual_diagonal-abs(v(:,j))**2)
        end do

        deallocate(z)
        allocate(z(n,n-r))
        if (k > r) z(:,1:k-r) = v(:,r+1:k)
        if (n > k) z(:,k-r+1:n-r) = v(:,k+1:n)
        if (.not. finite_complex_matrix(z)) then
            z = cmplx(0.0_dp,0.0_dp,dp)
            status = linalg_status_no_convergence
        end if
    end subroutine null_space_complex

    subroutine subspace_angles_real(a, b, angles, status, rcond, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! first real column-space matrix, m-by-n
        real(dp), intent(in) :: b(:, :) !! second real column-space matrix, m-by-k
        real(dp), allocatable, intent(out) :: angles(:) !! principal angles, descending, radians
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: rcond !! rank cutoff for both inputs, >= 0
        integer, intent(in), optional :: max_sweeps !! positive Jacobi sweep limit

        call subspace_angles_complex(cmplx(a,0.0_dp,dp), cmplx(b,0.0_dp,dp), &
            angles, status, rcond, max_sweeps)
    end subroutine subspace_angles_real

    subroutine subspace_angles_real_complex(a, b, angles, status, rcond, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! first real column-space matrix, m-by-n
        complex(dp), intent(in) :: b(:, :) !! second complex column-space matrix, m-by-k
        real(dp), allocatable, intent(out) :: angles(:) !! principal angles, descending, radians
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: rcond !! rank cutoff for both inputs, >= 0
        integer, intent(in), optional :: max_sweeps !! positive Jacobi sweep limit

        call subspace_angles_complex(cmplx(a,0.0_dp,dp), b, angles, status, rcond, max_sweeps)
    end subroutine subspace_angles_real_complex

    subroutine subspace_angles_complex_real(a, b, angles, status, rcond, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! first complex column-space matrix, m-by-n
        real(dp), intent(in) :: b(:, :) !! second real column-space matrix, m-by-k
        real(dp), allocatable, intent(out) :: angles(:) !! principal angles, descending, radians
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: rcond !! rank cutoff for both inputs, >= 0
        integer, intent(in), optional :: max_sweeps !! positive Jacobi sweep limit

        call subspace_angles_complex(a, cmplx(b,0.0_dp,dp), angles, status, rcond, max_sweeps)
    end subroutine subspace_angles_complex_real

    subroutine subspace_angles_complex(a, b, angles, status, rcond, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! first column-space matrix, m-by-n
        complex(dp), intent(in) :: b(:, :) !! second column-space matrix, m-by-k
        real(dp), allocatable, intent(out) :: angles(:) !! principal angles, descending, radians
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: rcond !! rank cutoff for both inputs, >= 0
        integer, intent(in), optional :: max_sweeps !! positive Jacobi sweep limit
        complex(dp), allocatable :: qa(:, :), qb(:, :), overlap(:, :), residual(:, :)
        complex(dp), allocatable :: us(:, :), vhs(:, :)
        real(dp), allocatable :: cosines(:), sines(:)
        integer :: ra, rb, h, m, i
        real(dp) :: cosine, sine

        allocate(angles(0))
        status = linalg_status_invalid_shape
        if (size(a,1) /= size(b,1)) return
        if (.not. valid_options(size(a,1),size(a,2),rcond,max_sweeps)) return
        if (.not. valid_options(size(b,1),size(b,2),rcond,max_sweeps)) return
        call orth_complex(a,qa,status,rcond,max_sweeps)
        if (status /= linalg_status_success) return
        call orth_complex(b,qb,status,rcond,max_sweeps)
        if (status /= linalg_status_success) return
        m = size(a,1)
        ra = size(qa,2)
        rb = size(qb,2)
        h = min(ra,rb)
        if (h == 0) return
        allocate(overlap(ra,rb),cosines(h),sines(h))
        overlap = matmul(transpose(conjg(qa)),qb)
        allocate(us(ra,h),vhs(h,rb))
        call svd(overlap,us,cosines,vhs,status,max_sweeps)
        if (status /= linalg_status_success) return
        deallocate(us,vhs)

        ! Singular values of overlap give cosines; near zero angles use
        ! singular values of the projection residual to avoid cancellation.
        ! The smaller column-space basis is projected off the larger space.
        if (ra >= rb) then
            allocate(residual(m,rb))
            residual = qb-matmul(qa,overlap)
        else
            allocate(residual(m,ra))
            residual = qa-matmul(qb,transpose(conjg(overlap)))
        end if
        allocate(us(m,h),vhs(h,h))
        call svd(residual,us,sines,vhs,status,max_sweeps)
        if (status /= linalg_status_success) return
        deallocate(angles)
        allocate(angles(h))
        do i = 1, h
            cosine = min(1.0_dp,max(0.0_dp,cosines(h-i+1)))
            sine = min(1.0_dp,max(0.0_dp,sines(i)))
            if (cosine**2 >= 0.5_dp) then
                angles(i) = asin(sine)
            else
                angles(i) = acos(cosine)
            end if
        end do
        if (.not. all(ieee_is_finite(angles))) then
            angles = 0.0_dp
            status = linalg_status_no_convergence
        end if
    end subroutine subspace_angles_complex

end module scifort_linalg_subspaces
