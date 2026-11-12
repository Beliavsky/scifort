! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciFort extension: continuity of nearly degenerate canonical-mode subspaces.
! Use principal angles for block matching and a common complex unitary
! Procrustes transform for the two groups' paired canonical coefficients.
! Rotated columns are SUBSPACE BASES, not individual canonical eigenvectors.
! Optional conditioned-spectral-metric comparison is congruence-invariant
! under nonsingular, frequency-independent channel coordinate changes.
module scifort_signal_canonical_subspace_tracking
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, &
        canonical_coherence_tf_result
    implicit none
    private

    integer, parameter, public :: canonical_subspace_success = 0
    integer, parameter, public :: canonical_subspace_invalid_input = 1

    type, public :: canonical_subspace_tracks_result
        real(dp), allocatable :: frequencies(:) !! original frequency coordinates
        real(dp), allocatable :: squared_modes(:,:) !! original descending strengths, not rotated
        complex(dp), allocatable :: basis_a(:,:,:) !! (frequency,A channel,rank), aligned basis
        complex(dp), allocatable :: basis_b(:,:,:) !! (frequency,B channel,rank), same rotation
        integer, allocatable :: track_id(:,:) !! persistent block ID, 0 at invalid bins
        integer, allocatable :: block_first(:,:) !! first original rank of containing block
        integer, allocatable :: block_size(:,:) !! dimension of containing block
        real(dp), allocatable :: principal_cosines(:,:) !! descending cosines within each block
        real(dp), allocatable :: similarity(:,:) !! minimum principal cosine of matched block
        logical, allocatable :: matched(:,:) !! parent block linked to preceding bin
        logical, allocatable :: segment_start(:,:) !! newly labeled block at this bin
        logical, allocatable :: valid(:) !! original spectral validity
        logical :: covariance_metric_used = .false. !! spectral inner product selected
        integer :: source_status = 0 !! canonical source's status
        integer :: status = canonical_subspace_invalid_input !! operation status
    end type canonical_subspace_tracks_result

    type, public :: canonical_subspace_tracks_tf_result
        real(dp), allocatable :: frequencies(:) !! original frequency coordinates
        real(dp), allocatable :: times(:) !! original time coordinates
        real(dp), allocatable :: squared_modes(:,:,:) !! original (frequency,time,rank) strengths
        complex(dp), allocatable :: basis_a(:,:,:,:) !! (frequency,time,A channel,rank)
        complex(dp), allocatable :: basis_b(:,:,:,:) !! (frequency,time,B channel,rank)
        integer, allocatable :: track_id(:,:,:) !! block ID, independently numbered per time
        integer, allocatable :: block_first(:,:,:) !! first rank of containing block
        integer, allocatable :: block_size(:,:,:) !! dimension of containing block
        real(dp), allocatable :: principal_cosines(:,:,:) !! ordered cosines of matched block
        real(dp), allocatable :: similarity(:,:,:) !! weakest principal cosine of block
        logical, allocatable :: matched(:,:,:) !! linked to preceding frequency at same time
        logical, allocatable :: segment_start(:,:,:) !! new block at frequency and time
        logical, allocatable :: valid(:,:) !! original per-cell validity
        logical :: covariance_metric_used = .false. !! spectral inner product selected
        integer :: source_status = 0 !! canonical source's status
        integer :: status = canonical_subspace_invalid_input !! operation status
    end type canonical_subspace_tracks_tf_result

    interface track_canonical_subspaces
        module procedure track_subspaces_static
        module procedure track_subspaces_tf
    end interface track_canonical_subspaces
    public :: track_canonical_subspaces

contains

    ! Local small-matrix complex Hermitian Jacobi eigensolver, adapted from
    ! this repository's canonical_coherence module to admit arbitrary PSD
    ! eigenvalue ranges. Avoids an external BLAS/LAPACK dependency.
    pure subroutine psd_eigen(a,lambda,vectors,ok)
        complex(dp), intent(in) :: a(:,:) !! finite Hermitian positive-semidefinite matrix
        real(dp), intent(out) :: lambda(:) !! descending eigenvalues of a
        complex(dp), intent(out) :: vectors(:,:) !! unitary matrix of eigenvectors
        logical, intent(out) :: ok !! true if Jacobi rotations converged
        complex(dp) :: h(size(a,1),size(a,1)), z, oldp, oldq, vtmp(size(a,1))
        real(dp) :: magnitude, scale, tolerance, tau, tangent, c, s, saved
        integer :: n, i, j, p, q, sweep
        logical :: converged

        n = size(a,1)
        lambda = 0.0_dp
        vectors = (0.0_dp,0.0_dp)
        ok = .false.
        if (n < 1 .or. size(a,2) /= n) return
        h = 0.5_dp*(a+transpose(conjg(a)))
        scale = maxval(abs(h))
        if (.not. ieee_is_finite(scale)) return
        do i = 1, n
            vectors(i,i) = (1.0_dp,0.0_dp)
        end do
        if (scale == 0.0_dp) then
            ok = .true.
            return
        end if
        tolerance = 16.0_dp*epsilon(1.0_dp)*scale
        converged = (n == 1)
        do sweep = 1, 100
            if (n == 1) exit
            converged = .true.
            do p = 1, n-1
                do q = p+1, n
                    magnitude = abs(h(p,q))
                    if (magnitude <= tolerance) cycle
                    converged = .false.
                    tau = 0.5_dp*(real(h(q,q),dp)-real(h(p,p),dp))/magnitude
                    tangent = sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau))
                    c = 1.0_dp/sqrt(1.0_dp+tangent*tangent)
                    s = tangent*c
                    z = h(p,q)/magnitude
                    do i = 1, n
                        if (i == p .or. i == q) cycle
                        oldp = h(i,p)
                        oldq = h(i,q)
                        h(i,p) = c*oldp-s*conjg(z)*oldq
                        h(i,q) = s*z*oldp+c*oldq
                        h(p,i) = conjg(h(i,p))
                        h(q,i) = conjg(h(i,q))
                    end do
                    h(p,p) = cmplx(real(h(p,p),dp)-tangent*magnitude,0.0_dp,dp)
                    h(q,q) = cmplx(real(h(q,q),dp)+tangent*magnitude,0.0_dp,dp)
                    h(p,q) = (0.0_dp,0.0_dp)
                    h(q,p) = (0.0_dp,0.0_dp)
                    do i = 1, n
                        oldp = vectors(i,p)
                        oldq = vectors(i,q)
                        vectors(i,p) = c*oldp-s*conjg(z)*oldq
                        vectors(i,q) = s*z*oldp+c*oldq
                    end do
                end do
            end do
            if (converged) exit
        end do
        if (.not. converged) return
        do i = 1, n
            saved = real(h(i,i),dp)
            if (.not. ieee_is_finite(saved)) return
            if (saved < -128.0_dp*epsilon(1.0_dp)*scale) return
            lambda(i) = max(0.0_dp,saved)
        end do
        do i = 2, n
            saved = lambda(i)
            vtmp = vectors(:,i)
            j = i-1
            do while (j >= 1)
                if (lambda(j) >= saved) exit
                lambda(j+1) = lambda(j)
                vectors(:,j+1) = vectors(:,j)
                j = j-1
            end do
            lambda(j+1) = saved
            vectors(:,j+1) = vtmp
        end do
        ok = .true.
    end subroutine psd_eigen

    ! A single scaling per whole A and B block preserves its span under
    ! arbitrary unitary column mixing, unlike column-by-column scaling.
    pure subroutine stacked_basis(a,b,x,ok)
        complex(dp), intent(in) :: a(:,:) !! paired coefficient block of group A
        complex(dp), intent(in) :: b(:,:) !! paired coefficient block of group B
        complex(dp), intent(out) :: x(:,:) !! concatenated, group-balanced block
        logical, intent(out) :: ok !! finite inputs and nonzero block norms
        real(dp) :: sa, sb, scale_a, scale_b, n

        ok = .false.
        x = (0.0_dp,0.0_dp)
        if (.not. all(ieee_is_finite(real(a,dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        if (.not. all(ieee_is_finite(real(b,dp)))) return
        if (.not. all(ieee_is_finite(aimag(b)))) return
        n = real(size(a,2),dp)
        scale_a = maxval(abs(a))
        scale_b = maxval(abs(b))
        if (.not. all(ieee_is_finite([scale_a,scale_b]))) return
        if (min(scale_a,scale_b) <= tiny(1.0_dp)) return
        sa = sqrt(sum((abs(a)/scale_a)**2)/n)
        sb = sqrt(sum((abs(b)/scale_b)**2)/n)
        if (min(sa,sb) <= tiny(1.0_dp)) return
        x(1:size(a,1),:) = (a/scale_a)/(sqrt(2.0_dp)*sa)
        x(size(a,1)+1:,:) = (b/scale_b)/(sqrt(2.0_dp)*sb)
        ok = .true.
    end subroutine stacked_basis

    ! Two-pass modified Gram-Schmidt; only used to find Euclidean principal
    ! angles of block spans, not to rotate the physical coefficient bases.
    pure subroutine orthonormalize(x,q,ok)
        complex(dp), intent(in) :: x(:,:) !! full-column-rank stacked coefficient block
        complex(dp), intent(out) :: q(:,:) !! Euclidean orthonormal span basis
        logical, intent(out) :: ok !! all block columns are independent
        complex(dp) :: vector(size(x,1))
        real(dp) :: length
        integer :: j, k, pass

        q = (0.0_dp,0.0_dp)
        ok = .false.
        do j = 1, size(x,2)
            vector = x(:,j)
            do pass = 1, 2
                do k = 1, j-1
                    vector = vector-sum(conjg(q(:,k))*vector)*q(:,k)
                end do
            end do
            length = sqrt(sum(abs(vector)**2))
            if (.not. ieee_is_finite(length)) return
            if (length < 1.0e-11_dp*sqrt(sum(abs(x(:,j))**2))) return
            q(:,j) = vector/length
        end do
        ok = .true.
    end subroutine orthonormalize

    ! The principal angles are the singular values of Qprev^H Qnext.
    ! The optimal joint A/B gauge is the unitary Procrustes factor of the
    ! *raw* stacked coefficient matrices, not their QR factorizations.
    pure subroutine compare_blocks(prev_a,prev_b,cur_a,cur_b,cosines,rotation,ok)
        complex(dp), intent(in) :: prev_a(:,:) !! previous aligned A block
        complex(dp), intent(in) :: prev_b(:,:) !! previous aligned B block
        complex(dp), intent(in) :: cur_a(:,:) !! new canonical A block
        complex(dp), intent(in) :: cur_b(:,:) !! new canonical B block
        real(dp), intent(out) :: cosines(:) !! descending principal-angle cosines
        complex(dp), intent(out) :: rotation(:,:) !! Procrustes right unitary factor
        logical, intent(out) :: ok !! comparable full-rank blocks
        complex(dp) :: xp(size(prev_a,1)+size(prev_b,1),size(prev_a,2))
        complex(dp) :: xc(size(cur_a,1)+size(cur_b,1),size(cur_a,2))
        complex(dp) :: qp(size(xp,1),size(xp,2)), qc(size(xc,1),size(xc,2))
        complex(dp) :: h(size(cosines),size(cosines)), gram(size(cosines),size(cosines))
        complex(dp) :: eigvec(size(cosines),size(cosines)), inverse(size(cosines),size(cosines))
        real(dp) :: eig(size(cosines)), floor_value
        logical :: good
        integer :: i, d

        d = size(cosines)
        ok = .false.
        cosines = 0.0_dp
        rotation = (0.0_dp,0.0_dp)
        if (d < 1) return
        call stacked_basis(prev_a,prev_b,xp,good)
        if (.not. good) return
        call stacked_basis(cur_a,cur_b,xc,good)
        if (.not. good) return
        call orthonormalize(xp,qp,good)
        if (.not. good) return
        call orthonormalize(xc,qc,good)
        if (.not. good) return
        h = matmul(transpose(conjg(qp)),qc)
        gram = matmul(transpose(conjg(h)),h)
        call psd_eigen(gram,eig,eigvec,good)
        if (.not. good) return
        cosines = sqrt(max(0.0_dp,min(1.0_dp,eig)))

        ! The polar factor of Hraw^H is V U^H for Hraw=U Sigma V^H.
        ! Full-rank block comparisons above threshold remain well-conditioned.
        h = matmul(transpose(conjg(xp)),xc)
        gram = matmul(transpose(conjg(h)),h)
        call psd_eigen(gram,eig,eigvec,good)
        if (.not. good) return
        floor_value = max(1.0e-12_dp,128.0_dp*epsilon(1.0_dp))
        if (minval(eig) <= floor_value**2) return
        inverse = (0.0_dp,0.0_dp)
        do i = 1, d
            inverse(i,i) = cmplx(1.0_dp/sqrt(eig(i)),0.0_dp,dp)
        end do
        rotation = matmul(matmul(eigvec,matmul(inverse,transpose(conjg(eigvec)))), &
            transpose(conjg(h)))
        ok = .true.
    end subroutine compare_blocks

    ! Strictly Hermitian positive-definite metric with a relative pivot check.
    ! The matrix is used as an inner product on coefficient vectors.
    pure subroutine metric_cholesky(a,l,ok)
        complex(dp), intent(in) :: a(:,:) !! Hermitian positive-definite conditional spectrum
        complex(dp), intent(out) :: l(:,:) !! Cholesky factor, lower triangular
        logical, intent(out) :: ok !! finite, sufficiently Hermitian and positive definite
        real(dp) :: scale, tolerance, pivot
        complex(dp) :: residual
        integer :: i, j, k, n

        n = size(a,1)
        ok = .false.
        l = (0.0_dp,0.0_dp)
        if (n < 1 .or. size(a,2) /= n) return
        if (.not. all(ieee_is_finite(real(a,dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        scale = maxval(abs(a))
        if (scale <= 0.0_dp) return
        tolerance = 64.0_dp*epsilon(1.0_dp)*scale
        if (maxval(abs(a-transpose(conjg(a)))) > tolerance) return
        do i = 1, n
            do j = 1, i-1
                residual = a(i,j)
                do k = 1, j-1
                    residual = residual-l(i,k)*conjg(l(j,k))
                end do
                l(i,j) = residual/real(l(j,j),dp)
            end do
            pivot = real(a(i,i),dp)
            do k = 1, i-1
                pivot = pivot-abs(l(i,k))**2
            end do
            if (.not. ieee_is_finite(pivot)) return
            if (pivot <= tolerance) return
            l(i,i) = cmplx(sqrt(pivot),0.0_dp,dp)
        end do
        ok = .true.
    end subroutine metric_cholesky

    ! Use the midpoint conditional auto-spectrum as a COMMON inner product
    ! at adjacent frequency bins. The same metric must be used for both
    ! subspaces; whitening each frequency separately is NOT invariant to a
    ! physical change in the covariance. Cholesky factors need not be stored.
    ! Under S -> T S T^H and c -> T^-H c, all overlap Gram matrices remain
    ! unchanged (even though individual whitened coordinates may rotate).
    pure subroutine compare_blocks_metric(prev_a,prev_b,cur_a,cur_b, &
            prev_sa,prev_sb,cur_sa,cur_sb,cosines,rotation,ok)
        complex(dp), intent(in) :: prev_a(:,:) !! previous aligned A coefficients
        complex(dp), intent(in) :: prev_b(:,:) !! previous aligned B coefficients
        complex(dp), intent(in) :: cur_a(:,:) !! current unaligned A coefficients
        complex(dp), intent(in) :: cur_b(:,:) !! current unaligned B coefficients
        complex(dp), intent(in) :: prev_sa(:,:) !! previous conditioned A spectral metric
        complex(dp), intent(in) :: prev_sb(:,:) !! previous conditioned B spectral metric
        complex(dp), intent(in) :: cur_sa(:,:) !! current conditioned A spectral metric
        complex(dp), intent(in) :: cur_sb(:,:) !! current conditioned B spectral metric
        real(dp), intent(out) :: cosines(:) !! ordered principal-angle cosines
        complex(dp), intent(out) :: rotation(:,:) !! optimal common right unitary
        logical, intent(out) :: ok !! both midpoint metric blocks and overlap are valid
        complex(dp) :: la(size(prev_a,1),size(prev_a,1))
        complex(dp) :: lb(size(prev_b,1),size(prev_b,1))
        complex(dp) :: pa(size(prev_a,1),size(prev_a,2))
        complex(dp) :: pb(size(prev_b,1),size(prev_b,2))
        complex(dp) :: ca(size(cur_a,1),size(cur_a,2))
        complex(dp) :: cb(size(cur_b,1),size(cur_b,2))
        logical :: good

        ok = .false.
        cosines = 0.0_dp
        rotation = (0.0_dp,0.0_dp)
        call metric_cholesky(0.5_dp*(prev_sa+cur_sa),la,good)
        if (.not. good) return
        call metric_cholesky(0.5_dp*(prev_sb+cur_sb),lb,good)
        if (.not. good) return
        pa = matmul(transpose(conjg(la)),prev_a)
        pb = matmul(transpose(conjg(lb)),prev_b)
        ca = matmul(transpose(conjg(la)),cur_a)
        cb = matmul(transpose(conjg(lb)),cur_b)
        call compare_blocks(pa,pb,ca,cb,cosines,rotation,ok)
    end subroutine compare_blocks_metric

    pure logical function valid_metric_series(sa,sb,valid) result(ok)
        complex(dp), intent(in) :: sa(:,:,:) !! (frequency,A,A), conditioned spectra
        complex(dp), intent(in) :: sb(:,:,:) !! (frequency,B,B), conditioned spectra
        logical, intent(in) :: valid(:) !! cells requiring positive-definite metrics
        complex(dp) :: la(size(sa,2),size(sa,2)), lb(size(sb,2),size(sb,2))
        logical :: good
        integer :: f

        ok = .false.
        if (size(sa,1) /= size(valid) .or. size(sb,1) /= size(valid)) return
        if (size(sa,2) < 1 .or. size(sb,2) < 1) return
        if (size(sa,3) /= size(sa,2) .or. size(sb,3) /= size(sb,2)) return
        do f = 1, size(valid)
            if (.not. valid(f)) cycle
            call metric_cholesky(sa(f,:,:),la,good)
            if (.not. good) return
            call metric_cholesky(sb(f,:,:),lb,good)
            if (.not. good) return
        end do
        ok = .true.
    end function valid_metric_series

    ! Global maximum-weight assignment on a square similarity table.
    ! Candidate scores also include one unit per eligible match so that
    ! maximizing cardinality takes precedence over aggregate similarity.
    pure subroutine assign_blocks(scores,assigned)
        real(dp), intent(in) :: scores(:,:) !! finite nonnegative match rewards
        integer, intent(out) :: assigned(:) !! chosen current block for each prior block
        real(dp) :: u(0:size(scores,1)), v(0:size(scores,1))
        real(dp) :: minima(0:size(scores,1)), change, current, ceiling
        integer :: path(0:size(scores,1)), owner(0:size(scores,1))
        integer :: i, j, j0, j1, row, n
        logical :: used(0:size(scores,1))

        n = size(scores,1)
        assigned = 0
        u = 0.0_dp
        v = 0.0_dp
        path = 0
        owner = 0
        ceiling = max(1.0_dp,maxval(scores))
        do i = 1, n
            owner(0) = i
            j0 = 0
            used = .false.
            minima = huge(1.0_dp)
            do
                used(j0) = .true.
                row = owner(j0)
                change = huge(1.0_dp)
                j1 = 0
                do j = 1, n
                    if (used(j)) cycle
                    current = ceiling-scores(row,j)-u(row)-v(j)
                    if (current < minima(j)) then
                        minima(j) = current
                        path(j) = j0
                    end if
                    if (minima(j) < change) then
                        change = minima(j)
                        j1 = j
                    end if
                end do
                do j = 0, n
                    if (used(j)) then
                        u(owner(j)) = u(owner(j))+change
                        v(j) = v(j)-change
                    else
                        minima(j) = minima(j)-change
                    end if
                end do
                j0 = j1
                if (owner(j0) == 0) exit
            end do
            do
                j1 = path(j0)
                owner(j0) = owner(j1)
                j0 = j1
                if (j0 == 0) exit
            end do
        end do
        do j = 1, n
            assigned(owner(j)) = j
        end do
    end subroutine assign_blocks

    pure subroutine find_blocks(values,gap,starts,lengths,nblocks)
        real(dp), intent(in) :: values(:) !! descending squared canonical strengths
        real(dp), intent(in) :: gap !! relative near-tie threshold, nonnegative
        integer, intent(out) :: starts(:) !! first rank of each consecutive block
        integer, intent(out) :: lengths(:) !! count of ranks in each block
        integer, intent(out) :: nblocks !! number of blocks
        integer :: i, n
        real(dp) :: scale

        n = size(values)
        starts = 0
        lengths = 0
        nblocks = 0
        i = 1
        do while (i <= n)
            nblocks = nblocks+1
            starts(nblocks) = i
            lengths(nblocks) = 1
            do while (i < n)
                if (.not. ieee_is_finite(values(i))) exit
                if (.not. ieee_is_finite(values(i+1))) exit
                if (values(i) < 0.0_dp .or. values(i+1) < 0.0_dp) exit
                scale = max(values(i),values(i+1),0.05_dp)
                if (abs(values(i)-values(i+1)) > gap*scale) exit
                i = i+1
                lengths(nblocks) = lengths(nblocks)+1
            end do
            i = i+1
        end do
    end subroutine find_blocks

    pure logical function nonzero_block(values,a,b) result(ok)
        real(dp), intent(in) :: values(:) !! squared singular strengths in one block
        complex(dp), intent(in) :: a(:,:) !! group A coefficients
        complex(dp), intent(in) :: b(:,:) !! group B coefficients
        complex(dp) :: x(size(a,1)+size(b,1),size(a,2))
        logical :: good

        ok = .false.
        if (.not. all(ieee_is_finite(values))) return
        if (minval(values) <= 256.0_dp*epsilon(1.0_dp)) return
        call stacked_basis(a,b,x,good)
        ok = good
    end function nonzero_block

    pure subroutine track_sequence(freq,values,a,b,valid,gap,minimum, &
            out_a,out_b,track_id,block_first,block_size,cosines,similarity,matched,starts, &
            metric_a,metric_b)
        real(dp), intent(in) :: freq(:) !! ordered real frequencies; wraps interrupt matching
        real(dp), intent(in) :: values(:,:) !! (frequency,mode) sorted strengths
        complex(dp), intent(in) :: a(:,:,:) !! (frequency,A channel,mode) source coefficients
        complex(dp), intent(in) :: b(:,:,:) !! (frequency,B channel,mode) source coefficients
        logical, intent(in) :: valid(:) !! frequency-bin validity
        real(dp), intent(in) :: gap !! within-block relative squared-strength gap
        real(dp), intent(in) :: minimum !! minimum weakest principal cosine in [0,1]
        complex(dp), intent(out) :: out_a(:,:,:) !! aligned output bases
        complex(dp), intent(out) :: out_b(:,:,:) !! aligned output bases
        integer, intent(out) :: track_id(:,:) !! persistent block labels
        integer, intent(out) :: block_first(:,:) !! initial source strength rank per block
        integer, intent(out) :: block_size(:,:) !! block dimensions
        real(dp), intent(out) :: cosines(:,:) !! principal cosines per matched block
        real(dp), intent(out) :: similarity(:,:) !! weakest principal cosine per matched block
        logical, intent(out) :: matched(:,:) !! successfully matched to previous cell
        logical, intent(out) :: starts(:,:) !! new label created at this cell
        complex(dp), intent(in), optional :: metric_a(:,:,:) !! conditioned A spectrum
        complex(dp), intent(in), optional :: metric_b(:,:,:) !! conditioned B spectrum
        integer :: nm, nf, f, nold, nnew, i, j, d, p, q, maxb, next_id
        integer :: old_begin(size(values,2)), old_size(size(values,2))
        integer :: new_begin(size(values,2)), new_size(size(values,2))
        integer :: assignment(size(values,2)), matched_old(size(values,2))
        real(dp) :: rewards(size(values,2),size(values,2))
        real(dp) :: angles(size(values,2)), nan
        complex(dp) :: rotation(size(values,2),size(values,2))
        logical :: ok, can_link, old_good(size(values,2)), new_good(size(values,2))

        nf = size(values,1)
        nm = size(values,2)
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        out_a = a
        out_b = b
        track_id = 0
        block_first = 0
        block_size = 0
        cosines = nan
        similarity = nan
        matched = .false.
        starts = .false.
        next_id = 0
        nold = 0
        do f = 1, nf
            if (.not. valid(f)) then
                nold = 0
                cycle
            end if
            call find_blocks(values(f,:),gap,new_begin,new_size,nnew)
            new_good = .false.
            do j = 1, nnew
                q = new_begin(j)
                d = new_size(j)
                new_good(j) = nonzero_block(values(f,q:q+d-1),a(f,:,q:q+d-1), &
                    b(f,:,q:q+d-1))
                block_first(f,q:q+d-1) = q
                block_size(f,q:q+d-1) = d
            end do
            can_link = (f > 1) .and. nold > 0
            if (f > 1) then
                if (.not. ieee_is_finite(freq(f))) can_link = .false.
                if (.not. ieee_is_finite(freq(f-1))) can_link = .false.
                if (freq(f) <= freq(f-1)) can_link = .false.
            end if
            matched_old = 0
            if (can_link) then
                maxb = max(nold,nnew)
                rewards = 0.0_dp
                do i = 1, nold
                    p = old_begin(i)
                    d = old_size(i)
                    if (.not. old_good(i)) cycle
                    do j = 1, nnew
                        q = new_begin(j)
                        if (d /= new_size(j) .or. .not. new_good(j)) cycle
                        if (present(metric_a)) then
                            call compare_blocks_metric(out_a(f-1,:,p:p+d-1), &
                                out_b(f-1,:,p:p+d-1),a(f,:,q:q+d-1), &
                                b(f,:,q:q+d-1),metric_a(f-1,:,:),metric_b(f-1,:,:), &
                                metric_a(f,:,:),metric_b(f,:,:), &
                                angles(1:d),rotation(1:d,1:d),ok)
                        else
                            call compare_blocks(out_a(f-1,:,p:p+d-1), &
                                out_b(f-1,:,p:p+d-1),a(f,:,q:q+d-1), &
                                b(f,:,q:q+d-1),angles(1:d),rotation(1:d,1:d),ok)
                        end if
                        if (.not. ok) cycle
                        if (minval(angles(1:d)) <= 64.0_dp*epsilon(1.0_dp)) cycle
                        if (minval(angles(1:d)) < minimum) cycle
                        rewards(i,j) = 1.0_dp+minval(angles(1:d))/real(maxb+1,dp)
                    end do
                end do
                call assign_blocks(rewards(1:maxb,1:maxb),assignment(1:maxb))
                do i = 1, nold
                    j = assignment(i)
                    if (j < 1 .or. j > nnew) cycle
                    if (rewards(i,j) == 0.0_dp) cycle
                    matched_old(j) = i
                end do
            end if
            do j = 1, nnew
                q = new_begin(j)
                d = new_size(j)
                i = matched_old(j)
                if (i > 0) then
                    p = old_begin(i)
                    if (present(metric_a)) then
                        call compare_blocks_metric(out_a(f-1,:,p:p+d-1), &
                            out_b(f-1,:,p:p+d-1),a(f,:,q:q+d-1), &
                            b(f,:,q:q+d-1),metric_a(f-1,:,:),metric_b(f-1,:,:), &
                            metric_a(f,:,:),metric_b(f,:,:), &
                            angles(1:d),rotation(1:d,1:d),ok)
                    else
                        call compare_blocks(out_a(f-1,:,p:p+d-1), &
                            out_b(f-1,:,p:p+d-1),a(f,:,q:q+d-1), &
                            b(f,:,q:q+d-1),angles(1:d),rotation(1:d,1:d),ok)
                    end if
                    if (ok) then
                        track_id(f,q:q+d-1) = track_id(f-1,p)
                        out_a(f,:,q:q+d-1) = matmul(a(f,:,q:q+d-1),rotation(1:d,1:d))
                        out_b(f,:,q:q+d-1) = matmul(b(f,:,q:q+d-1),rotation(1:d,1:d))
                        cosines(f,q:q+d-1) = angles(1:d)
                        similarity(f,q:q+d-1) = minval(angles(1:d))
                        matched(f,q:q+d-1) = .true.
                        cycle
                    end if
                end if
                next_id = next_id+1
                track_id(f,q:q+d-1) = next_id
                starts(f,q:q+d-1) = .true.
            end do
            old_begin = new_begin
            old_size = new_size
            old_good = new_good
            nold = nnew
        end do
    end subroutine track_sequence

    function track_subspaces_static(source,gap_tolerance,min_similarity,covariance_metric) result(out)
        type(canonical_coherence_result), intent(in) :: source !! canonical result with vectors
        real(dp), intent(in), optional :: gap_tolerance !! near-tie relative gap, default 0.001
        real(dp), intent(in), optional :: min_similarity !! weakest cosine floor, [0,1], default 0.5
        logical, intent(in), optional :: covariance_metric !! use conditional spectral inner product
        type(canonical_subspace_tracks_result) :: out
        integer :: nf, nm, na, nb
        real(dp) :: gap, minimum
        logical :: use_metric

        gap = 0.001_dp
        if (present(gap_tolerance)) gap = gap_tolerance
        minimum = 0.5_dp
        if (present(min_similarity)) minimum = min_similarity
        if (.not. all(ieee_is_finite([gap,minimum]))) return
        if (gap < 0.0_dp .or. gap > 1.0_dp) return
        if (minimum < 0.0_dp .or. minimum > 1.0_dp) return
        use_metric = .false.
        if (present(covariance_metric)) use_metric = covariance_metric
        if (.not. allocated(source%frequencies)) return
        if (.not. allocated(source%squared_modes)) return
        if (.not. allocated(source%coefficients_a)) return
        if (.not. allocated(source%coefficients_b)) return
        if (.not. allocated(source%valid)) return
        nf = size(source%frequencies)
        nm = size(source%squared_modes,2)
        na = size(source%coefficients_a,2)
        nb = size(source%coefficients_b,2)
        if (min(nf,nm,na,nb) < 1) return
        if (any(shape(source%squared_modes) /= [nf,nm])) return
        if (any(shape(source%coefficients_a) /= [nf,na,nm])) return
        if (any(shape(source%coefficients_b) /= [nf,nb,nm])) return
        if (size(source%valid) /= nf) return
        if (use_metric) then
            if (.not. allocated(source%covariance_a)) return
            if (.not. allocated(source%covariance_b)) return
            if (any(shape(source%covariance_a) /= [nf,na,na])) return
            if (any(shape(source%covariance_b) /= [nf,nb,nb])) return
            if (.not. valid_metric_series(source%covariance_a,source%covariance_b,source%valid)) return
        end if
        allocate(out%frequencies(nf),out%squared_modes(nf,nm), &
            out%basis_a(nf,na,nm),out%basis_b(nf,nb,nm),out%track_id(nf,nm), &
            out%block_first(nf,nm),out%block_size(nf,nm), &
            out%principal_cosines(nf,nm),out%similarity(nf,nm), &
            out%matched(nf,nm),out%segment_start(nf,nm),out%valid(nf))
        out%frequencies = source%frequencies
        out%squared_modes = source%squared_modes
        out%valid = source%valid
        out%source_status = source%status
        if (use_metric) then
            call track_sequence(source%frequencies,source%squared_modes,source%coefficients_a, &
                source%coefficients_b,source%valid,gap,minimum, &
                out%basis_a,out%basis_b,out%track_id,out%block_first,out%block_size, &
                out%principal_cosines,out%similarity,out%matched,out%segment_start, &
                source%covariance_a,source%covariance_b)
        else
            call track_sequence(source%frequencies,source%squared_modes,source%coefficients_a, &
                source%coefficients_b,source%valid,gap,minimum, &
                out%basis_a,out%basis_b,out%track_id,out%block_first,out%block_size, &
                out%principal_cosines,out%similarity,out%matched,out%segment_start)
        end if
        out%covariance_metric_used = use_metric
        out%status = canonical_subspace_success
    end function track_subspaces_static

    function track_subspaces_tf(source,gap_tolerance,min_similarity,covariance_metric) result(out)
        type(canonical_coherence_tf_result), intent(in) :: source !! source with vectors
        real(dp), intent(in), optional :: gap_tolerance !! relative tolerance >= 0, default 0.001
        real(dp), intent(in), optional :: min_similarity !! weakest cosine [0,1], default 0.5
        logical, intent(in), optional :: covariance_metric !! use conditional spectral inner product
        type(canonical_subspace_tracks_tf_result) :: out
        integer :: nf, nt, nm, na, nb, t
        real(dp) :: gap, minimum
        logical :: use_metric

        gap = 0.001_dp
        if (present(gap_tolerance)) gap = gap_tolerance
        minimum = 0.5_dp
        if (present(min_similarity)) minimum = min_similarity
        if (.not. all(ieee_is_finite([gap,minimum]))) return
        if (gap < 0.0_dp .or. gap > 1.0_dp) return
        if (minimum < 0.0_dp .or. minimum > 1.0_dp) return
        use_metric = .false.
        if (present(covariance_metric)) use_metric = covariance_metric
        if (.not. allocated(source%frequencies)) return
        if (.not. allocated(source%times)) return
        if (.not. allocated(source%squared_modes)) return
        if (.not. allocated(source%coefficients_a)) return
        if (.not. allocated(source%coefficients_b)) return
        if (.not. allocated(source%valid)) return
        nf = size(source%frequencies)
        nt = size(source%times)
        nm = size(source%squared_modes,3)
        na = size(source%coefficients_a,3)
        nb = size(source%coefficients_b,3)
        if (min(nf,nt,nm,na,nb) < 1) return
        if (any(shape(source%squared_modes) /= [nf,nt,nm])) return
        if (any(shape(source%coefficients_a) /= [nf,nt,na,nm])) return
        if (any(shape(source%coefficients_b) /= [nf,nt,nb,nm])) return
        if (any(shape(source%valid) /= [nf,nt])) return
        if (use_metric) then
            if (.not. allocated(source%covariance_a)) return
            if (.not. allocated(source%covariance_b)) return
            if (any(shape(source%covariance_a) /= [nf,nt,na,na])) return
            if (any(shape(source%covariance_b) /= [nf,nt,nb,nb])) return
            do t = 1, nt
                if (.not. valid_metric_series(source%covariance_a(:,t,:,:), &
                        source%covariance_b(:,t,:,:),source%valid(:,t))) return
            end do
        end if
        allocate(out%frequencies(nf),out%times(nt),out%squared_modes(nf,nt,nm), &
            out%basis_a(nf,nt,na,nm),out%basis_b(nf,nt,nb,nm), &
            out%track_id(nf,nt,nm),out%block_first(nf,nt,nm), &
            out%block_size(nf,nt,nm),out%principal_cosines(nf,nt,nm), &
            out%similarity(nf,nt,nm),out%matched(nf,nt,nm), &
            out%segment_start(nf,nt,nm),out%valid(nf,nt))
        out%frequencies = source%frequencies
        out%times = source%times
        out%squared_modes = source%squared_modes
        out%valid = source%valid
        out%source_status = source%status
        do t = 1, nt
            if (use_metric) then
                call track_sequence(source%frequencies,source%squared_modes(:,t,:), &
                    source%coefficients_a(:,t,:,:),source%coefficients_b(:,t,:,:), &
                    source%valid(:,t),gap,minimum,out%basis_a(:,t,:,:),out%basis_b(:,t,:,:), &
                    out%track_id(:,t,:),out%block_first(:,t,:),out%block_size(:,t,:), &
                    out%principal_cosines(:,t,:),out%similarity(:,t,:), &
                    out%matched(:,t,:),out%segment_start(:,t,:), &
                    source%covariance_a(:,t,:,:),source%covariance_b(:,t,:,:))
            else
                call track_sequence(source%frequencies,source%squared_modes(:,t,:), &
                    source%coefficients_a(:,t,:,:),source%coefficients_b(:,t,:,:), &
                    source%valid(:,t),gap,minimum,out%basis_a(:,t,:,:),out%basis_b(:,t,:,:), &
                    out%track_id(:,t,:),out%block_first(:,t,:),out%block_size(:,t,:), &
                    out%principal_cosines(:,t,:),out%similarity(:,t,:), &
                    out%matched(:,t,:),out%segment_start(:,t,:))
            end if
        end do
        out%covariance_metric_used = use_metric
        out%status = canonical_subspace_success
    end function track_subspaces_tf
end module scifort_signal_canonical_subspace_tracking
