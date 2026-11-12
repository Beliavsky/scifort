! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Partial SVD of real-shape, complex-data CSR/CSC matrices.  This is a
! dependency-free complex Golub--Kahan/one-sided Jacobi implementation,
! adapted independently from SciFort's real sparse SVD.  No normal matrix
! A^H A is formed.  It uses a bounded (non-restarted) trial subspace.
module scifort_sparse_linalg_complex_svds
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csr_complex_matvec, &
        csr_complex_hermitian, csc_complex_to_csr
    implicit none
    private

    integer, parameter, public :: complex_svds_success = 0
    integer, parameter, public :: complex_svds_no_convergence = 1
    integer, parameter, public :: complex_svds_invalid_input = -1
    integer, parameter, public :: complex_svds_breakdown = -2
    public :: svds_complex

    interface svds_complex
        module procedure svds_complex_csr, svds_complex_csc
    end interface svds_complex
contains
    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate complex value.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    pure real(dp) function complex_norm(z) result(value)
        complex(dp), intent(in) :: z(:) !! Complex vector.
        value = norm2(abs(z))
    end function complex_norm

    subroutine svds_complex_csc(a,k,u,s,vh,info,which,tol,ncv,maxiter,v0,iterations,residual_norms)
        type(csc_complex_matrix), intent(in) :: a !! Canonical complex CSC input (m,n).
        integer, intent(in) :: k !! Requested number of singular triplets, 1..min(m,n).
        complex(dp), allocatable, intent(out) :: u(:,:) !! Left singular vectors (m,k).
        real(dp), allocatable, intent(out) :: s(:) !! Singular values, ascending (k).
        complex(dp), allocatable, intent(out) :: vh(:,:) !! Conjugate-transposed right vectors (k,n).
        integer, intent(out) :: info !! 0 success, 1 nonconvergence, negative invalid/breakdown.
        character(len=*), intent(in), optional :: which !! 'LM' (default) or 'SM'.
        real(dp), intent(in), optional :: tol !! True two-sided residual tolerance, default 1e-8.
        integer, intent(in), optional :: ncv !! Trial subspace width; 1..min(m,n).
        integer, intent(in), optional :: maxiter !! Maximum Lanczos steps, <= ncv.
        complex(dp), intent(in), optional :: v0(:) !! Starting vector on shorter side.
        integer, intent(out), optional :: iterations !! Number of Golub--Kahan steps.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True combined residuals.
        type(csr_complex_matrix) :: csr
        integer :: status

        info = complex_svds_invalid_input
        if (present(iterations)) iterations = 0
        allocate(u(0,0),s(0),vh(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        deallocate(u,s,vh)
        if (present(residual_norms)) deallocate(residual_norms)
        call svds_complex_csr(csr,k,u,s,vh,info,which,tol,ncv,maxiter,v0,iterations,residual_norms)
    end subroutine svds_complex_csc

    subroutine svds_complex_csr(a,k,u,s,vh,info,which,tol,ncv,maxiter,v0,iterations,residual_norms)
        type(csr_complex_matrix), intent(in) :: a !! Canonical complex CSR input (m,n).
        integer, intent(in) :: k !! Requested singular triplets; 1..min(m,n).
        complex(dp), allocatable, intent(out) :: u(:,:) !! Left singular vectors, (m,k).
        real(dp), allocatable, intent(out) :: s(:) !! Singular values, ascending.
        complex(dp), allocatable, intent(out) :: vh(:,:) !! Conjugate-transposed right vectors (k,n).
        integer, intent(out) :: info !! 0 success, 1 nonconvergence, negative invalid/breakdown.
        character(len=*), intent(in), optional :: which !! 'LM' (default) or 'SM'.
        real(dp), intent(in), optional :: tol !! Relative two-sided residual tolerance, default 1e-8.
        integer, intent(in), optional :: ncv !! Maximum basis width, default min(q,max(40,8*k+16)).
        integer, intent(in), optional :: maxiter !! Maximum basis steps; default ncv.
        complex(dp), intent(in), optional :: v0(:) !! Starting vector on shorter side, length q.
        integer, intent(out), optional :: iterations !! Steps performed.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True combined residual norms.
        type(csr_complex_matrix) :: op, oph
        complex(dp), allocatable :: vbasis(:,:), ubasis(:,:), work(:), nextv(:), &
            left(:,:), right(:,:)
        real(dp), allocatable :: singular(:), residuals(:)
        real(dp) :: alpha, beta, prev_beta, threshold, tolerance, scale
        integer :: m, n, q, p, width, steps, j, i, pass, status
        character(len=2) :: target
        logical :: is_transposed, converged, valid, seeded

        info = complex_svds_invalid_input
        if (present(iterations)) iterations = 0
        allocate(u(0,0),s(0),vh(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csr_complex(a)) return
        m = a%n_rows
        n = a%n_cols
        q = min(m,n)
        p = max(m,n)
        if (q < 1 .or. k < 1 .or. k > q) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LM' .and. target /= 'SM') return
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        tolerance = max(tolerance,100.0_dp*epsilon(1.0_dp))
        width = min(q,max(40,8*k+16))
        if (target == 'SM') width = q
        if (present(ncv)) width = ncv
        if (width < k .or. width > q) return
        ! No reliable global smallest-singular triplets without a full trial basis.
        if (target == 'SM' .and. width /= q) return
        steps = width
        if (present(maxiter)) steps = min(width,maxiter)
        if (steps < k) return
        if (target == 'SM' .and. steps /= q) return
        if (present(v0)) then
            if (size(v0) /= q) return
            if (.not. all(finite_complex(v0))) return
            if (complex_norm(v0) == 0.0_dp) return
        end if
        scale = 0.0_dp
        if (size(a%data) > 0) scale = maxval(abs(a%data))
        if (scale == 0.0_dp) then
            deallocate(u,s,vh)
            allocate(u(m,k),s(k),vh(k,n))
            u = cmplx(0.0_dp,0.0_dp,dp)
            s = 0.0_dp
            vh = cmplx(0.0_dp,0.0_dp,dp)
            do j = 1,k
                u(j,j) = cmplx(1.0_dp,0.0_dp,dp)
                vh(j,j) = cmplx(1.0_dp,0.0_dp,dp)
            end do
            if (present(residual_norms)) then
                deallocate(residual_norms)
                allocate(residual_norms(k))
                residual_norms = 0.0_dp
            end if
            info = complex_svds_success
            return
        end if
        threshold = 100.0_dp*epsilon(1.0_dp)*scale
        is_transposed = m < n
        if (is_transposed) then
            op = csr_complex_hermitian(a,status)
            if (status /= 0) then
                info = complex_svds_breakdown
                return
            end if
            oph = a
        else
            op = a
            oph = csr_complex_hermitian(a,status)
            if (status /= 0) then
                info = complex_svds_breakdown
                return
            end if
        end if
        allocate(vbasis(q,steps),ubasis(p,steps),work(p),nextv(q))
        vbasis = cmplx(0.0_dp,0.0_dp,dp)
        ubasis = cmplx(0.0_dp,0.0_dp,dp)
        if (present(v0)) then
            vbasis(:,1) = v0/complex_norm(v0)
        else
            do i = 1,q
                vbasis(i,1) = cmplx(sin(0.6180339887498948_dp*real(i+1,dp)) &
                    +cos(0.3316624790355399_dp*real(i+3,dp)), &
                    sin(0.313_dp*real(i+5,dp)),dp)
            end do
            vbasis(:,1) = vbasis(:,1)/complex_norm(vbasis(:,1))
        end if
        prev_beta = 0.0_dp
        info = complex_svds_no_convergence
        do j = 1,steps
            work = csr_complex_matvec(op,vbasis(:,j),status)
            if (status /= 0) then
                info = complex_svds_breakdown
                exit
            end if
            if (j > 1) work = work-prev_beta*ubasis(:,j-1)
            do pass = 1,2
                do i = 1,j-1
                    work = work-dot_product(ubasis(:,i),work)*ubasis(:,i)
                end do
            end do
            alpha = complex_norm(work)
            beta = 0.0_dp
            if (alpha > threshold) then
                ubasis(:,j) = work/alpha
                nextv = csr_complex_matvec(oph,ubasis(:,j),status)
                if (status /= 0) then
                    info = complex_svds_breakdown
                    exit
                end if
                nextv = nextv-alpha*vbasis(:,j)
                do pass = 1,2
                    do i = 1,j
                        nextv = nextv-dot_product(vbasis(:,i),nextv)*vbasis(:,i)
                    end do
                end do
                beta = complex_norm(nextv)
            end if
            if (.not. ieee_is_finite(alpha) .or. .not. ieee_is_finite(beta)) then
                info = complex_svds_breakdown
                exit
            end if
            if (j >= k .and. (j == steps .or. mod(j,8) == 0)) then
                if (target == 'LM' .or. j == q) then
                    call projected_svds(op,oph,vbasis(:,:j),k,target,tolerance, &
                        left,singular,right,residuals,converged,valid)
                    if (.not. valid) then
                        info = complex_svds_breakdown
                        exit
                    end if
                    deallocate(u,s,vh)
                    allocate(u(m,k),s(k),vh(k,n))
                    s = singular
                    if (is_transposed) then
                        u = right
                        vh = conjg(transpose(left))
                    else
                        u = left
                        vh = conjg(transpose(right))
                    end if
                    if (present(residual_norms)) then
                        deallocate(residual_norms)
                        residual_norms = residuals
                    end if
                    ! When the full shorter-side basis is affordable,
                    ! finish it before declaring success. A single Krylov
                    ! chain can miss multiplicities (for example, the second
                    ! of two identical largest singular values) even when
                    ! the Ritz vectors already found have tiny residuals.
                    if (converged .and. (steps < q .or. j == steps)) then
                        info = complex_svds_success
                        exit
                    end if
                end if
            end if
            if (j == steps) exit
            if (beta > threshold) then
                vbasis(:,j+1) = nextv/beta
                prev_beta = beta
            else
                call independent_seed(vbasis(:,:j),nextv,seeded)
                if (.not. seeded) then
                    info = complex_svds_breakdown
                    exit
                end if
                vbasis(:,j+1) = nextv
                prev_beta = 0.0_dp
            end if
        end do
        if (present(iterations)) iterations = min(j,steps)
    end subroutine svds_complex_csr

    subroutine independent_seed(q,v,ok)
        complex(dp), intent(in) :: q(:,:) !! Previously orthonormal directions (n,j).
        complex(dp), intent(out) :: v(:) !! New unit direction outside their span.
        logical, intent(out) :: ok !! Whether an independent direction exists.
        real(dp) :: available, largest, length
        integer :: i, best, j, pass

        best = 0
        largest = 0.0_dp
        do i = 1,size(q,1)
            available = max(0.0_dp,1.0_dp-sum(abs(q(i,:))**2))
            if (available > largest) then
                largest = available
                best = i
            end if
        end do
        ok = .false.
        v = cmplx(0.0_dp,0.0_dp,dp)
        if (best == 0 .or. largest < 1000.0_dp*epsilon(1.0_dp)) return
        v(best) = cmplx(1.0_dp,0.0_dp,dp)
        do pass = 1,2
            do j = 1,size(q,2)
                v = v-dot_product(q(:,j),v)*q(:,j)
            end do
        end do
        length = complex_norm(v)
        if (length < 100.0_dp*epsilon(1.0_dp)) return
        v = v/length
        ok = .true.
    end subroutine independent_seed

    subroutine projected_svds(a,ah,basis,k,target,tolerance,left,singular,right,residuals,converged,valid)
        type(csr_complex_matrix), intent(in) :: a !! Tall p-by-q sparse operator.
        type(csr_complex_matrix), intent(in) :: ah !! Conjugate transpose of a.
        complex(dp), intent(in) :: basis(:,:) !! Orthonormal right trial basis, (q,t).
        integer, intent(in) :: k !! Requested number of triplets.
        character(len=2), intent(in) :: target !! LM or SM.
        real(dp), intent(in) :: tolerance !! Two-sided true-residual tolerance.
        complex(dp), allocatable, intent(out) :: left(:,:) !! Complex left vectors, (p,k).
        real(dp), allocatable, intent(out) :: singular(:) !! Singular values, (k), ascending.
        complex(dp), allocatable, intent(out) :: right(:,:) !! Complex right vectors, (q,k).
        real(dp), allocatable, intent(out) :: residuals(:) !! Combined true residuals.
        logical, intent(out) :: converged !! True iff every requested triplet converged.
        logical, intent(out) :: valid !! Projected SVD converged and is finite.
        complex(dp), allocatable :: w(:,:), rotation(:,:), ul(:,:), &
            left_candidate(:), probe(:), v(:), av(:), ahu(:)
        real(dp), allocatable :: sl(:)
        integer, allocatable :: order(:)
        real(dp) :: pp, qq, crossabs, tau, tangent, cs, sn, sigmax, cutoff, &
            bound, rleft, rright, probe_norm, best_norm, length
        complex(dp) :: cross, phase, candidate_i, candidate_j
        integer :: p, q, t, i, j, h, sweep, pass, choice, best, old, status
        logical :: changed, small_converged

        valid = .false.
        converged = .false.
        p = a%n_rows
        q = a%n_cols
        t = size(basis,2)
        allocate(w(p,t),rotation(t,t),sl(t),ul(p,t),order(t), &
            left_candidate(p),probe(p),v(q),av(p),ahu(q))
        rotation = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1,t
            rotation(j,j) = cmplx(1.0_dp,0.0_dp,dp)
            w(:,j) = csr_complex_matvec(a,basis(:,j),status)
            if (status /= 0) return
        end do
        small_converged = .false.
        do sweep = 1,80
            changed = .false.
            do i = 1,t-1
                do j = i+1,t
                    pp = real(dot_product(w(:,i),w(:,i)),dp)
                    qq = real(dot_product(w(:,j),w(:,j)),dp)
                    cross = dot_product(w(:,i),w(:,j))
                    crossabs = abs(cross)
                    if (.not. all(ieee_is_finite([pp,qq,crossabs]))) return
                    if (pp <= 0.0_dp .or. qq <= 0.0_dp) cycle
                    if (crossabs <= 8.0_dp*epsilon(1.0_dp)*sqrt(pp)*sqrt(qq)) cycle
                    changed = .true.
                    phase = cross/crossabs
                    tau = (qq-pp)/(2.0_dp*crossabs)
                    tangent = sign(1.0_dp,tau)/(abs(tau)+hypot(1.0_dp,tau))
                    cs = 1.0_dp/sqrt(1.0_dp+tangent*tangent)
                    sn = cs*tangent
                    ! Apply the same complex-unitary Jacobi rotation to W and V.
                    left_candidate = cs*w(:,i)-sn*conjg(phase)*w(:,j)
                    w(:,j) = sn*phase*w(:,i)+cs*w(:,j)
                    w(:,i) = left_candidate
                    do h = 1,t
                        candidate_i = rotation(h,i)
                        candidate_j = rotation(h,j)
                        rotation(h,i) = cs*candidate_i-sn*conjg(phase)*candidate_j
                        rotation(h,j) = sn*phase*candidate_i+cs*candidate_j
                    end do
                end do
            end do
            if (.not. changed) then
                small_converged = .true.
                exit
            end if
        end do
        if (.not. small_converged) return
        do j = 1,t
            sl(j) = complex_norm(w(:,j))
        end do
        if (.not. all(ieee_is_finite(sl))) return
        sigmax = maxval(sl)
        cutoff = 100.0_dp*epsilon(1.0_dp)*sigmax
        ul = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1,t
            if (sl(j) > cutoff) then
                ul(:,j) = w(:,j)/sl(j)
            else
                sl(j) = 0.0_dp
            end if
        end do
        ! Complete left vectors for numerical null singular values.
        do j = 1,t
            if (sl(j) /= 0.0_dp) cycle
            best = 0
            best_norm = 0.0_dp
            do i = 1,p
                probe = cmplx(0.0_dp,0.0_dp,dp)
                probe(i) = cmplx(1.0_dp,0.0_dp,dp)
                do pass = 1,2
                    do h = 1,t
                        if (h == j) cycle
                        if (any(ul(:,h) /= cmplx(0.0_dp,0.0_dp,dp))) &
                            probe = probe-dot_product(ul(:,h),probe)*ul(:,h)
                    end do
                end do
                probe_norm = complex_norm(probe)
                if (probe_norm > best_norm) then
                    best_norm = probe_norm
                    best = i
                end if
            end do
            if (best == 0 .or. best_norm < 100.0_dp*epsilon(1.0_dp)) return
            ul(best,j) = cmplx(1.0_dp,0.0_dp,dp)
            do pass = 1,2
                do h = 1,t
                    if (h == j) cycle
                    ul(:,j) = ul(:,j)-dot_product(ul(:,h),ul(:,j))*ul(:,h)
                end do
            end do
            length = complex_norm(ul(:,j))
            if (length <= 100.0_dp*epsilon(1.0_dp)) return
            ul(:,j) = ul(:,j)/length
        end do
        order = [(i,i=1,t)]
        do i = 1,t-1
            choice = i
            do j = i+1,t
                if (target == 'LM') then
                    if (sl(order(j)) > sl(order(choice))) choice = j
                else
                    if (sl(order(j)) < sl(order(choice))) choice = j
                end if
            end do
            old = order(i)
            order(i) = order(choice)
            order(choice) = old
        end do
        allocate(left(p,k),singular(k),right(q,k),residuals(k))
        do j = 1,k
            h = order(j)
            singular(j) = sl(h)
            left(:,j) = ul(:,h)
            right(:,j) = matmul(basis,rotation(:,h))
        end do
        do i = 1,k-1
            do j = i+1,k
                if (singular(j) >= singular(i)) cycle
                pp = singular(i)
                singular(i) = singular(j)
                singular(j) = pp
                left_candidate = left(:,i)
                left(:,i) = left(:,j)
                left(:,j) = left_candidate
                v = right(:,i)
                right(:,i) = right(:,j)
                right(:,j) = v
            end do
        end do
        converged = .true.
        bound = tolerance*max(sigmax,tiny(1.0_dp))
        do j = 1,k
            av = csr_complex_matvec(a,right(:,j),status)
            if (status /= 0) return
            ahu = csr_complex_matvec(ah,left(:,j),status)
            if (status /= 0) return
            rleft = complex_norm(av-singular(j)*left(:,j))
            rright = complex_norm(ahu-singular(j)*right(:,j))
            residuals(j) = hypot(rleft,rright)
            if (.not. ieee_is_finite(residuals(j))) return
            if (residuals(j) > bound) converged = .false.
        end do
        valid = .true.
    end subroutine projected_svds
end module scifort_sparse_linalg_complex_svds
