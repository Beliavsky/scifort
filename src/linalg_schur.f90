! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Public real and complex Schur decomposition for dense matrices.
! Real Schur output uses conjugate-pair reordering and nested real invariant
! subspaces. The Schur engine and adjacent unitary Schur swaps are original
! MIT-licensed SciFort algorithms.
module scifort_linalg_schur
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_matrix_functions, only : dense_complex_schur, matrix_function_max_order
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: schur, schur_selector, rsf2csf

    abstract interface
        logical function schur_selector(eigenvalue)
            import dp
            complex(dp), intent(in) :: eigenvalue !! eigenvalue in the original matrix's units
        end function schur_selector
    end interface

    interface schur
        module procedure schur_real_input
        module procedure schur_complex_input
        module procedure schur_real_output
    end interface schur

contains

    ! Convert an already computed *real* Schur form to a complex Schur form.
    ! Each 2x2 real block is triangularized by a local complex unitary
    ! similarity, and the same rotation is accumulated into the Schur basis.
    ! This costs O(n**2) per two-by-two block and does not rerun Schur QR.
    ! The real input T and Z are never modified.
    subroutine rsf2csf(tr, zr, tc, zc, status)
        real(dp), intent(in) :: tr(:, :) !! real quasi-upper-triangular Schur form (n,n)
        real(dp), intent(in) :: zr(:, :) !! real orthogonal Schur vectors (n,n)
        complex(dp), intent(out) :: tc(:, :) !! complex upper-triangular Schur form (n,n)
        complex(dp), intent(out) :: zc(:, :) !! complex unitary Schur vectors (n,n)
        integer, intent(out) :: status !! success, invalid shape, or numerical failure

        real(dp) :: magnitude, local_tol, delta, discriminant, lower, upper
        real(dp) :: diag_left, diag_right, rotation_norm, sine, result_tol
        complex(dp) :: mu, cosine, root, u, v
        integer :: n, i, j, k
        logical :: active(size(tr,1)), previous_active

        tc = cmplx(0.0_dp, 0.0_dp, dp)
        zc = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(tr,1)
        if (size(tr,2) /= n .or. size(zr,1) /= n .or. size(zr,2) /= n) return
        if (size(tc,1) /= n .or. size(tc,2) /= n) return
        if (size(zc,1) /= n .or. size(zc,2) /= n) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(tr)) .or. .not. all(ieee_is_finite(zr))) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        ! Work on a normalized matrix: rotations and intermediate updates
        ! cannot overflow when the input has very large entries.
        magnitude = maxval(abs(tr))
        zc = cmplx(zr, 0.0_dp, kind=dp)
        if (magnitude <= 0.0_dp) then
            status = linalg_status_success
            return
        end if
        tc = cmplx(tr / magnitude, 0.0_dp, kind=dp)

        ! The accepted real Schur form has disjoint 1x1 or 2x2 blocks.
        ! Tiny roundoff below the first subdiagonal is permitted.
        result_tol = 128.0_dp * epsilon(1.0_dp)
        do j = 1, n-2
            if (any(abs(tc(j+2:n,j)) > result_tol)) return
            tc(j+2:n,j) = cmplx(0.0_dp, 0.0_dp, dp)
        end do
        active = .false.
        previous_active = .false.
        do j = 1, n-1
            local_tol = 16.0_dp * epsilon(1.0_dp) * &
                (abs(tc(j,j)) + abs(tc(j+1,j+1)))
            active(j) = abs(tc(j+1,j)) > local_tol
            if (active(j) .and. previous_active) return
            previous_active = active(j)
        end do

        ! Bottom-up processing is convenient because later rotations only
        ! change columns to the right of earlier real Schur blocks.
        do k = n-1, 1, -1
            if (.not. active(k)) then
                tc(k+1,k) = cmplx(0.0_dp, 0.0_dp, dp)
                cycle
            end if
            diag_left = real(tc(k,k),dp)
            diag_right = real(tc(k+1,k+1),dp)
            upper = real(tc(k,k+1),dp)
            lower = real(tc(k+1,k),dp)
            delta = 0.5_dp*(diag_left-diag_right)
            discriminant = delta*delta + upper*lower
            if (discriminant < 0.0_dp) then
                root = cmplx(0.0_dp, sqrt(-discriminant), kind=dp)
            else
                ! For a noncanonical 2x2 block with real eigenvalues, pick
                ! the root with the larger eigenvector component to avoid
                ! loss of precision in delta +/- sqrt(discriminant).
                root = cmplx(sign(sqrt(discriminant), delta), 0.0_dp, kind=dp)
            end if
            mu = cmplx(delta, 0.0_dp, kind=dp) + root
            rotation_norm = hypot(abs(mu), abs(lower))
            if (.not. ieee_is_finite(rotation_norm) .or. rotation_norm <= 0.0_dp) return
            cosine = mu / rotation_norm
            sine = lower / rotation_norm

            ! Left multiply by G=[conjg(c),s;-s,c], then right multiply
            ! by G^H=[c,-s;s,conjg(c)].  The first column of G^H is an
            ! eigenvector of the local 2x2 Schur block.
            do j = k, n
                u = tc(k,j)
                v = tc(k+1,j)
                tc(k,j) = conjg(cosine)*u + sine*v
                tc(k+1,j) = -sine*u + cosine*v
            end do
            do i = 1, k+1
                u = tc(i,k)
                v = tc(i,k+1)
                tc(i,k) = cosine*u + sine*v
                tc(i,k+1) = -sine*u + conjg(cosine)*v
            end do
            do i = 1, n
                u = zc(i,k)
                v = zc(i,k+1)
                zc(i,k) = cosine*u + sine*v
                zc(i,k+1) = -sine*u + conjg(cosine)*v
            end do
            ! The annihilated subdiagonal should be roundoff-sized.
            if (abs(tc(k+1,k)) > 512.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp,abs(tc(k,k)),abs(tc(k+1,k+1)))) return
            tc(k+1,k) = cmplx(0.0_dp, 0.0_dp, dp)
        end do
        tc = tc * magnitude
        if (.not. finite_complex(tc) .or. .not. finite_complex(zc)) then
            tc = cmplx(0.0_dp,0.0_dp,dp)
            zc = cmplx(0.0_dp,0.0_dp,dp)
            return
        end if
        status = linalg_status_success
    end subroutine rsf2csf


    subroutine schur_real_input(a, t, z, status, sort, sdim, selector)
        real(dp), intent(in) :: a(:, :) !! square real input, order 0..40
        complex(dp), intent(out) :: t(:, :) !! complex upper-triangular Schur factor (n,n)
        complex(dp), intent(out) :: z(:, :) !! complex unitary Schur basis (n,n)
        integer, intent(out) :: status !! success, invalid shape, or nonconvergence
        character(len=*), intent(in), optional :: sort !! 'none', 'lhp', 'rhp', 'iuc', or 'ouc'
        integer, intent(out), optional :: sdim !! number of reordered eigenvalues at the top left
        procedure(schur_selector), optional :: selector !! optional custom complex eigenvalue predicate

        t = cmplx(0.0_dp, 0.0_dp, dp)
        z = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(sdim)) sdim = 0
        status = linalg_status_invalid_shape
        if (.not. valid_shape(size(a,1), size(a,2), size(t,1), size(t,2), size(z,1), size(z,2))) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        call schur_core(cmplx(a, kind=dp), t, z, status, sort, sdim, selector)
    end subroutine schur_real_input

    subroutine schur_complex_input(a, t, z, status, sort, sdim, selector)
        complex(dp), intent(in) :: a(:, :) !! square complex input, order 0..40
        complex(dp), intent(out) :: t(:, :) !! complex upper-triangular Schur factor (n,n)
        complex(dp), intent(out) :: z(:, :) !! complex unitary Schur basis (n,n)
        integer, intent(out) :: status !! success, invalid shape, or nonconvergence
        character(len=*), intent(in), optional :: sort !! 'none', 'lhp', 'rhp', 'iuc', or 'ouc'
        integer, intent(out), optional :: sdim !! number of reordered eigenvalues at the top left
        procedure(schur_selector), optional :: selector !! optional custom complex eigenvalue predicate

        t = cmplx(0.0_dp, 0.0_dp, dp)
        z = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(sdim)) sdim = 0
        status = linalg_status_invalid_shape
        if (.not. valid_shape(size(a,1), size(a,2), size(t,1), size(t,2), size(z,1), size(z,2))) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(a,dp))) .or. .not. all(ieee_is_finite(aimag(a)))) return
        call schur_core(a, t, z, status, sort, sdim, selector)
    end subroutine schur_complex_input

    ! A real Schur form consists of 1x1 real blocks and 2x2 real blocks
    ! for complex-conjugate eigenpairs.  Reuse the complex Schur engine,
    ! reorder conjugates to adjacent positions, and construct nested real
    ! invariant subspaces.  The final similarity is wholly real.
    subroutine schur_real_output(a, t, z, status, sort, sdim, selector)
        real(dp), intent(in) :: a(:, :) !! square real matrix, order 0..40
        real(dp), intent(out) :: t(:, :) !! real quasi-upper-triangular Schur form (n,n)
        real(dp), intent(out) :: z(:, :) !! real orthogonal Schur vectors (n,n)
        integer, intent(out) :: status !! success, invalid shape, or nonconvergence
        character(len=*), intent(in), optional :: sort !! 'none', 'lhp', 'rhp', 'iuc', 'ouc'
        integer, intent(out), optional :: sdim !! number of leading selected eigenvalues
        procedure(schur_selector), optional :: selector !! complex-eigenvalue predicate; pairs stay together

        complex(dp), allocatable :: tc(:, :), zc(:, :)
        real(dp), allocatable :: scaled(:, :), q(:, :), r(:, :)
        integer, allocatable :: block_size(:)
        real(dp) :: magnitude, err, orth_err, structural_err, norm_a
        integer :: n, i, k, block_count, count, code
        character(len=:), allocatable :: rule
        logical, allocatable :: allowed(:, :)

        t = 0.0_dp
        z = 0.0_dp
        if (present(sdim)) sdim = 0
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (.not. valid_shape(n, size(a,2), size(t,1), size(t,2), size(z,1), size(z,2))) return
        rule = 'none'
        if (present(sort)) rule = trim(adjustl(sort))
        select case (rule)
        case ('none', 'lhp', 'rhp', 'iuc', 'ouc')
        case default
            return
        end select
        if (present(selector) .and. rule /= 'none') return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        magnitude = maxval(abs(a))
        allocate(scaled(n,n), tc(n,n), zc(n,n), q(n,n), r(n,n), block_size(n), allowed(n,n))
        scaled = 0.0_dp
        if (magnitude > 0.0_dp) scaled = a / magnitude
        ! Scaling also avoids overflowing the real similarity verification.
        call schur_core(cmplx(scaled,kind=dp), tc, zc, code)
        if (code /= linalg_status_success) return
        call pair_conjugate_eigenvalues(tc, zc, block_size, block_count, code)
        if (code /= linalg_status_success) return
        call sort_real_schur_blocks(tc, zc, block_size, block_count, magnitude, &
            rule, count, code, selector)
        if (code /= linalg_status_success) return
        call real_schur_basis(zc, block_size, block_count, q, code)
        if (code /= linalg_status_success) return

        r = matmul(transpose(q), matmul(scaled, q))
        allowed = .false.
        do i = 1, n
            allowed(i,i:n) = .true.
        end do
        k = 1
        do i = 1, block_count
            if (block_size(i) == 2) allowed(k+1,k) = .true.
            k = k+block_size(i)
        end do
        norm_a = max(1.0_dp, sqrt(sum(scaled**2)))
        structural_err = sqrt(sum(r**2, mask=.not. allowed)) / norm_a
        orth_err = sqrt(sum((matmul(transpose(q),q)-identity_real(n))**2)) / real(n,dp)
        if (.not. ieee_is_finite(structural_err) .or. .not. ieee_is_finite(orth_err)) return
        if (structural_err > 1.0e-8_dp .or. orth_err > 1.0e-8_dp) return
        where (.not. allowed) r = 0.0_dp
        err = sqrt(sum((matmul(matmul(q,r),transpose(q))-scaled)**2))/norm_a
        if (.not. ieee_is_finite(err) .or. err > 1.0e-8_dp) return
        t = r*magnitude
        if (.not. all(ieee_is_finite(t))) then
            t = 0.0_dp
            return
        end if
        z = q
        if (present(sdim)) sdim = count
        status = linalg_status_success
    end subroutine schur_real_output

    ! Form a stable block partition: 1x1 for numerically real eigenvalues,
    ! 2x2 for conjugate pairs.  Complex Schur swaps preserve unitarity.
    subroutine pair_conjugate_eigenvalues(t, z, block_size, block_count, status)
        complex(dp), intent(inout) :: t(:, :) !! normalized complex Schur factor
        complex(dp), intent(inout) :: z(:, :) !! complex Schur vectors
        integer, intent(out) :: block_size(:) !! sizes of successive real Schur blocks (1 or 2)
        integer, intent(out) :: block_count !! number of blocks
        integer, intent(out) :: status !! success or nonconvergence

        integer :: n, k, j, best_j, p, code
        real(dp) :: tol, distance, best_distance, spectral_scale
        complex(dp) :: w

        n = size(t,1)
        block_size = 0
        block_count = 0
        k = 1
        status = linalg_status_no_convergence
        spectral_scale = 1.0_dp
        do j = 1, n
            spectral_scale = max(spectral_scale, abs(t(j,j)))
        end do
        tol = 1024.0_dp*epsilon(1.0_dp)*spectral_scale
        do while (k <= n)
            w = t(k,k)
            block_count = block_count+1
            if (abs(aimag(w)) <= tol) then
                block_size(block_count) = 1
                k = k+1
                cycle
            end if
            best_distance = huge(1.0_dp)
            best_j = 0
            do j = k+1, n
                distance = abs(t(j,j)-conjg(w))
                if (distance < best_distance) then
                    best_distance = distance
                    best_j = j
                end if
            end do
            if (best_j == 0 .or. best_distance > tol) return
            do p = best_j-1, k+1, -1
                call swap_adjacent_schur(t, z, p, code)
                if (code /= linalg_status_success) return
            end do
            block_size(block_count) = 2
            k = k+2
        end do
        status = linalg_status_success
    end subroutine pair_conjugate_eigenvalues

    ! Sorting operates on whole real blocks: a predicate selecting one of
    ! a conjugate pair selects both eigenvalues (as in real-output SciPy).
    subroutine sort_real_schur_blocks(t, z, block_size, block_count, magnitude, &
        rule, selected, status, selector)
        complex(dp), intent(inout) :: t(:, :) !! normalized complex Schur factor
        complex(dp), intent(inout) :: z(:, :) !! complex Schur vectors
        integer, intent(inout) :: block_size(:) !! contiguous 1x1 / 2x2 block sizes
        integer, intent(in) :: block_count !! number of real blocks
        real(dp), intent(in) :: magnitude !! input normalization factor
        character(len=*), intent(in) :: rule !! supported selection rule
        integer, intent(out) :: selected !! selected eigenvalue count
        integer, intent(out) :: status !! success or nonconvergence
        procedure(schur_selector), optional :: selector !! selection predicate

        logical, allocatable :: chosen(:)
        logical :: selected_now, tmp_choice
        integer :: j, k, target, start, offset, p, code, selected_blocks, temp_size

        status = linalg_status_no_convergence
        selected = 0
        if (rule == 'none' .and. .not. present(selector)) then
            status = linalg_status_success
            return
        end if
        allocate(chosen(block_count))
        k = 1
        do j = 1, block_count
            if (block_size(j) == 1) then
                ! A nearly real Schur eigenvalue must be tested as real.
                ! Otherwise a roundoff-size imaginary part can spuriously
                ! trigger a one-sided custom predicate.
                selected_now = selected_eigenvalue( &
                    cmplx(real(t(k,k),dp)*magnitude,0.0_dp,dp), rule, selector)
            else
                selected_now = selected_eigenvalue(t(k,k)*magnitude, rule, selector)
            end if
            if (block_size(j) == 2) then
                selected_now = selected_now .or. &
                    selected_eigenvalue(t(k+1,k+1)*magnitude, rule, selector)
            end if
            chosen(j) = selected_now
            k = k+block_size(j)
        end do
        selected_blocks = 0
        do j = 1, block_count
            if (.not. chosen(j)) cycle
            target = 1+sum(block_size(1:selected_blocks))
            start = 1+sum(block_size(1:j-1))
            do offset = 0, block_size(j)-1
                do p = start+offset-1, target+offset, -1
                    call swap_adjacent_schur(t,z,p,code)
                    if (code /= linalg_status_success) return
                end do
            end do
            temp_size = block_size(j)
            tmp_choice = chosen(j)
            do p = j, selected_blocks+2, -1
                block_size(p) = block_size(p-1)
                chosen(p) = chosen(p-1)
            end do
            block_size(selected_blocks+1) = temp_size
            chosen(selected_blocks+1) = tmp_choice
            selected_blocks = selected_blocks+1
            selected = selected+temp_size
        end do
        ! Recheck membership after adjacent similarities, because rounding
        ! can move a Schur eigenvalue across a sorting boundary.
        k = 1
        do j = 1, block_count
            if (block_size(j) == 1) then
                selected_now = selected_eigenvalue( &
                    cmplx(real(t(k,k),dp)*magnitude,0.0_dp,dp), rule, selector)
            else
                selected_now = selected_eigenvalue(t(k,k)*magnitude, rule, selector) .or. &
                    selected_eigenvalue(t(k+1,k+1)*magnitude, rule, selector)
            end if
            if ((j <= selected_blocks) .neqv. selected_now) return
            k = k+block_size(j)
        end do
        status = linalg_status_success
    end subroutine sort_real_schur_blocks

    subroutine real_schur_basis(zc, block_size, block_count, q, status)
        complex(dp), intent(in) :: zc(:, :) !! ordered complex Schur vectors
        integer, intent(in) :: block_size(:) !! ordered real-block sizes
        integer, intent(in) :: block_count !! number of blocks
        real(dp), intent(out) :: q(:, :) !! constructed real orthogonal basis
        integer, intent(out) :: status !! success or nonconvergence

        real(dp) :: trial(size(zc,1)), best(size(zc,1)), candidate(size(zc,1))
        real(dp) :: strength, trial_norm
        integer :: k, ib, j, part, m, pass, n, offset, c

        q = 0.0_dp
        status = linalg_status_no_convergence
        n = size(zc,1)
        k = 1
        do ib = 1, block_count
            do offset = 0, block_size(ib)-1
                m = k+offset-1
                strength = 0.0_dp
                best = 0.0_dp
                do c = 1, 2*block_size(ib)
                    j = k+(c-1)/2
                    part = mod(c-1,2)
                    if (part == 0) then
                        candidate = real(zc(:,j),dp)
                    else
                        candidate = aimag(zc(:,j))
                    end if
                    trial = candidate
                    do pass = 1, 2
                        if (m > 0) then
                            trial = trial-matmul(q(:,1:m),matmul(transpose(q(:,1:m)),trial))
                        end if
                    end do
                    trial_norm = sqrt(sum(trial**2))
                    if (trial_norm > strength) then
                        best = trial
                        strength = trial_norm
                    end if
                end do
                if (.not. ieee_is_finite(strength) .or. strength <= 1.0e-9_dp) return
                q(:,m+1) = best/strength
            end do
            k = k+block_size(ib)
        end do
        status = linalg_status_success
    end subroutine real_schur_basis

    pure function identity_real(n) result(eye)
        integer, intent(in) :: n !! square matrix order
        real(dp) :: eye(n,n)
        integer :: i
        eye = 0.0_dp
        do i = 1, n
            eye(i,i) = 1.0_dp
        end do
    end function identity_real

    logical function valid_shape(n, m, nt, mt, nz, mz)
        integer, intent(in) :: n !! input rows
        integer, intent(in) :: m !! input columns
        integer, intent(in) :: nt !! Schur factor rows
        integer, intent(in) :: mt !! Schur factor columns
        integer, intent(in) :: nz !! transformation rows
        integer, intent(in) :: mz !! transformation columns

        valid_shape = n == m .and. n <= matrix_function_max_order .and. &
            nt == n .and. mt == n .and. nz == n .and. mz == n
    end function valid_shape

    subroutine schur_core(a, t, z, status, sort, sdim, selector)
        complex(dp), intent(in) :: a(:, :) !! finite square input matrix
        complex(dp), intent(out) :: t(:, :) !! complex upper-triangular Schur factor
        complex(dp), intent(out) :: z(:, :) !! unitary Schur basis
        integer, intent(out) :: status !! success, invalid selection, or nonconvergence
        character(len=*), intent(in), optional :: sort !! supported sorting rule
        integer, intent(out), optional :: sdim !! selected dimension
        procedure(schur_selector), optional :: selector !! custom eigenvalue predicate

        complex(dp), allocatable :: scaled(:, :), orig(:, :)
        real(dp) :: magnitude, err, norm_a
        integer :: n, i, j, selected, code
        logical :: upper, lower, reorder, chosen
        character(len=:), allocatable :: rule

        status = linalg_status_invalid_shape
        if (present(sdim)) sdim = 0
        rule = 'none'
        if (present(sort)) rule = trim(adjustl(sort))
        select case (rule)
        case ('none', 'lhp', 'rhp', 'iuc', 'ouc')
        case default
            return
        end select
        if (present(selector) .and. present(sort)) then
            if (rule /= 'none') return
        end if
        reorder = rule /= 'none' .or. present(selector)

        n = size(a,1)
        t = cmplx(0.0_dp,0.0_dp,dp)
        z = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, n
            z(i,i) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        status = linalg_status_success
        if (n == 0) return

        magnitude = max(maxval(abs(real(a,dp))), maxval(abs(aimag(a))))
        if (magnitude <= 0.0_dp) then
            if (reorder) then
                if (selected_eigenvalue(cmplx(0.0_dp,0.0_dp,dp),rule,selector)) then
                    if (present(sdim)) sdim = n
                end if
            end if
            return
        end if
        allocate(scaled(n,n), orig(n,n))
        orig = a/magnitude
        status = linalg_status_no_convergence
        upper = .true.
        lower = .true.
        do j = 1, n-1
            if (any(abs(a(j+1:n,j)) > 0.0_dp)) upper = .false.
            if (any(abs(a(j,j+1:n)) > 0.0_dp)) lower = .false.
        end do
        if (upper) then
            scaled = orig
        else if (lower) then
            scaled = orig(n:1:-1,n:1:-1)
            z = cmplx(0.0_dp,0.0_dp,dp)
            do i = 1, n
                z(i,n+1-i) = cmplx(1.0_dp,0.0_dp,dp)
            end do
        else
            call dense_complex_schur(orig, scaled, z, code)
            if (code /= linalg_status_success) return
        end if

        selected = 0
        if (reorder) then
            ! Move the selected eigenvalues to the upper-left while retaining
            ! their existing relative order. Use unitary similarity swaps.
            do j = 1, n
                chosen = selected_eigenvalue(scaled(j,j)*magnitude, rule, selector)
                if (.not. chosen) cycle
                selected = selected+1
                do i = j-1, selected, -1
                    call swap_adjacent_schur(scaled, z, i, code)
                    if (code /= linalg_status_success) return
                end do
            end do
            ! Make roundoff-induced changes of region membership explicit:
            ! SciPy likewise reports a failure if a sorted eigenvalue moves
            ! outside the selected region.
            do j = 1, n
                chosen = selected_eigenvalue(scaled(j,j)*magnitude, rule, selector)
                if ((j <= selected) .neqv. chosen) return
            end do
        end if
        err = sqrt(sum(abs(matmul(matmul(z,scaled),transpose(conjg(z)))-orig)**2))
        norm_a = sqrt(sum(abs(orig)**2))
        if (.not. ieee_is_finite(err)) return
        if (err/max(1.0_dp,norm_a) > 1.0e-8_dp) return
        t = scaled*magnitude
        if (.not. finite_complex(t) .or. .not. finite_complex(z)) then
            t = cmplx(0.0_dp,0.0_dp,dp)
            z = cmplx(0.0_dp,0.0_dp,dp)
            return
        end if
        if (present(sdim)) sdim = selected
        status = linalg_status_success
    end subroutine schur_core

    logical function selected_eigenvalue(w, rule, selector)
        complex(dp), intent(in) :: w !! candidate eigenvalue in the matrix's original units
        character(len=*), intent(in) :: rule !! built-in selection code
        procedure(schur_selector), optional :: selector !! user-defined eigenvalue predicate

        selected_eigenvalue = .false.
        if (present(selector)) then
            selected_eigenvalue = selector(w)
            return
        end if
        select case (rule)
        case ('lhp')
            selected_eigenvalue = real(w,dp) < 0.0_dp
        case ('rhp')
            selected_eigenvalue = real(w,dp) >= 0.0_dp
        case ('iuc')
            selected_eigenvalue = abs(w) <= 1.0_dp
        case ('ouc')
            selected_eigenvalue = abs(w) > 1.0_dp
        end select
    end function selected_eigenvalue

    ! Adjacent Schur swap from a unitary eigenvector of the lower eigenvalue.
    ! Unlike a permutation similarity, this preserves upper-triangular form.
    subroutine swap_adjacent_schur(t, z, p, status)
        complex(dp), intent(inout) :: t(:, :) !! scaled upper-triangular Schur factor
        complex(dp), intent(inout) :: z(:, :) !! unitary Schur vectors
        integer, intent(in) :: p !! swap positions p and p+1
        integer, intent(out) :: status !! success or nonconvergence

        complex(dp) :: q11, q12, q21, q22, a, b, u, v
        real(dp) :: pair_norm, local_scale, defect
        integer :: i, n

        status = linalg_status_no_convergence
        n = size(t,1)
        b = t(p,p+1)
        a = t(p+1,p+1)-t(p,p)
        pair_norm = hypot(abs(a),abs(b))
        if (.not. ieee_is_finite(pair_norm) .or. pair_norm <= 0.0_dp) return
        q11 = b/pair_norm
        q21 = a/pair_norm
        q12 = -conjg(q21)
        q22 = conjg(q11)
        do i = p, n
            u = t(p,i)
            v = t(p+1,i)
            t(p,i) = conjg(q11)*u+conjg(q21)*v
            t(p+1,i) = conjg(q12)*u+conjg(q22)*v
        end do
        do i = 1, p+1
            u = t(i,p)
            v = t(i,p+1)
            t(i,p) = u*q11+v*q21
            t(i,p+1) = u*q12+v*q22
        end do
        do i = 1, n
            u = z(i,p)
            v = z(i,p+1)
            z(i,p) = u*q11+v*q21
            z(i,p+1) = u*q12+v*q22
        end do
        defect = abs(t(p+1,p))
        local_scale = max(1.0_dp,abs(t(p,p)),abs(t(p+1,p+1)),abs(t(p,p+1)))
        if (.not. ieee_is_finite(defect) .or. &
            defect > 1024.0_dp*epsilon(1.0_dp)*local_scale) return
        t(p+1,p) = cmplx(0.0_dp,0.0_dp,dp)
        if (.not. finite_complex(t) .or. .not. finite_complex(z)) return
        status = linalg_status_success
    end subroutine swap_adjacent_schur

    pure logical function finite_complex(x)
        complex(dp), intent(in) :: x(:, :) !! matrix to check
        finite_complex = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function finite_complex

end module scifort_linalg_schur
