! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Squared canonical coherences for groups of spectral channels, conditional on
! a disjoint third group. Original SciFort extension, not a scipy.signal API.
! Whitening uses Hermitian Cholesky solves; a cyclic complex Jacobi iteration
! obtains eigenpairs of the smaller whitened cross-correlation Gram matrix.
! Optionally return variance-normalized canonical variate coefficients.
module scifort_signal_canonical_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_invalid_input, &
        partial_coherence_singular
    use scifort_signal_group_coherence, only : group_coherence_result, &
        group_coherence_tf_result, group_coherence_matrix, group_coherence_tf_matrix
    use scifort_signal_band_coherence, only : integrate_spectral_bands, &
        integrate_tf_spectral_bands
    implicit none
    private

    type, public :: canonical_coherence_result
        real(dp), allocatable :: frequencies(:) !! sample frequencies or integrated band centers
        real(dp), allocatable :: squared_modes(:,:) !! descending squared canonical correlations (freq,mode)
        complex(dp), allocatable :: coefficients_a(:,:,:) !! (frequency,group-A channel,mode), optional
        complex(dp), allocatable :: coefficients_b(:,:,:) !! (frequency,group-B channel,mode), optional
        complex(dp), allocatable :: covariance_a(:,:,:) !! conditioned/regularized A autospectra
        complex(dp), allocatable :: covariance_b(:,:,:) !! conditioned/regularized B autospectra
        logical, allocatable :: directions_stable(:,:) !! (frequency,mode), isolated nonzero roots
        real(dp), allocatable :: coherence(:) !! aggregate determinant group coherence
        real(dp), allocatable :: log_information(:) !! -log of Wilks' determinant ratio
        real(dp), allocatable :: shrinkage_used(:) !! diagonal shrinkage used at each frequency
        logical, allocatable :: valid(:) !! positive-definite and converged cells
        integer :: status = partial_coherence_invalid_input !! success, invalid or singular
    end type canonical_coherence_result

    type, public :: canonical_coherence_tf_result
        real(dp), allocatable :: frequencies(:) !! sampled frequency or band coordinates
        real(dp), allocatable :: times(:) !! window midpoint times
        real(dp), allocatable :: squared_modes(:,:,:) !! (frequency,time,mode), descending
        complex(dp), allocatable :: coefficients_a(:,:,:,:) !! (frequency,time,A channel,mode), optional
        complex(dp), allocatable :: coefficients_b(:,:,:,:) !! (frequency,time,B channel,mode), optional
        complex(dp), allocatable :: covariance_a(:,:,:,:) !! conditioned/regularized A autospectra
        complex(dp), allocatable :: covariance_b(:,:,:,:) !! conditioned/regularized B autospectra
        logical, allocatable :: directions_stable(:,:,:) !! (frequency,time,mode), optional
        real(dp), allocatable :: coherence(:,:) !! aggregate determinant coherence (freq,time)
        real(dp), allocatable :: log_information(:,:) !! negative log determinant ratio
        real(dp), allocatable :: shrinkage_used(:,:) !! per time-frequency-cell shrinkage
        logical, allocatable :: valid(:,:) !! positive-definite and converged cells
        integer :: status = partial_coherence_invalid_input !! success, invalid or singular
    end type canonical_coherence_tf_result

    interface multitaper_canonical_coherence
        module procedure multitaper_canonical_coherence_real
        module procedure multitaper_canonical_coherence_complex
    end interface
    interface multitaper_tf_canonical_coherence
        module procedure multitaper_tf_canonical_coherence_real
        module procedure multitaper_tf_canonical_coherence_complex
    end interface
    interface multitaper_band_canonical_coherence
        module procedure multitaper_band_canonical_coherence_real
        module procedure multitaper_band_canonical_coherence_complex
    end interface
    interface multitaper_tf_band_canonical_coherence
        module procedure multitaper_tf_band_canonical_coherence_real
        module procedure multitaper_tf_band_canonical_coherence_complex
    end interface

    public :: canonical_coherence_matrix, canonical_coherence_tf_matrix
    public :: band_canonical_coherence_matrix, band_canonical_coherence_tf_matrix
    public :: multitaper_canonical_coherence, multitaper_tf_canonical_coherence
    public :: multitaper_band_canonical_coherence, multitaper_tf_band_canonical_coherence

contains

    ! Cholesky of a strictly positive-definite Hermitian block; use the
    ! same relative pivot test as the determinant-based group estimator.
    pure subroutine cholesky_block(a,l,ok)
        complex(dp), intent(in) :: a(:,:) !! Hermitian positive-definite block
        complex(dp), intent(out) :: l(:,:) !! lower triangular Cholesky factor
        logical, intent(out) :: ok !! true if the factorization succeeded
        real(dp) :: scale, pivot, tolerance
        complex(dp) :: term
        integer :: n, i, j, k

        n = size(a,1)
        ok = .false.
        l = (0.0_dp,0.0_dp)
        if (n == 0) then
            ok = .true.
            return
        end if
        scale = 0.0_dp
        do i = 1, n
            scale = max(scale,abs(a(i,i)))
        end do
        if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) return
        tolerance = 64.0_dp*epsilon(1.0_dp)*scale
        do i = 1, n
            do j = 1, i-1
                term = a(i,j)
                do k = 1, j-1
                    term = term - l(i,k)*conjg(l(j,k))
                end do
                l(i,j) = term/real(l(j,j),dp)
            end do
            pivot = real(a(i,i),dp)
            do k = 1, i-1
                pivot = pivot - abs(l(i,k))**2
            end do
            if (.not. ieee_is_finite(pivot) .or. pivot <= tolerance) return
            l(i,i) = cmplx(sqrt(pivot),0.0_dp,dp)
        end do
        ok = .true.
    end subroutine cholesky_block

    pure subroutine lower_solve(l,b,x)
        complex(dp), intent(in) :: l(:,:) !! positive-diagonal lower triangular matrix
        complex(dp), intent(in) :: b(:,:) !! one or more right-hand-side columns
        complex(dp), intent(out) :: x(:,:) !! solutions l*x=b, same shape as b
        integer :: i, j, k

        x = b
        do j = 1, size(b,2)
            do i = 1, size(l,1)
                do k = 1, i-1
                    x(i,j) = x(i,j)-l(i,k)*x(k,j)
                end do
                x(i,j) = x(i,j)/real(l(i,i),dp)
            end do
        end do
    end subroutine lower_solve

    ! Convert whitened left/right singular vectors to canonical coefficients:
    ! solve L^H * coefficients = whitened vectors.
    pure subroutine upper_conjugate_solve(l,b,x)
        complex(dp), intent(in) :: l(:,:) !! nonsingular lower Cholesky factor
        complex(dp), intent(in) :: b(:,:) !! whitened canonical direction columns
        complex(dp), intent(out) :: x(:,:) !! covariance-normalized canonical coefficients
        integer :: i, j, k

        x = b
        do j = 1, size(b,2)
            do i = size(l,1), 1, -1
                do k = i+1, size(l,1)
                    x(i,j) = x(i,j)-conjg(l(k,i))*x(k,j)
                end do
                x(i,j) = x(i,j)/real(l(i,i),dp)
            end do
        end do
    end subroutine upper_conjugate_solve

    ! Deterministically complete an orthonormal set when a zero singular
    ! value leaves the corresponding opposite-side direction undetermined.
    pure subroutine orthogonal_complement(existing,vector,ok)
        complex(dp), intent(in) :: existing(:,:) !! previously assigned orthonormal directions
        complex(dp), intent(out) :: vector(:) !! selected unit complement
        logical, intent(out) :: ok !! a direction was available
        complex(dp) :: trial(size(vector)), best(size(vector))
        real(dp) :: length2, best_length2
        integer :: i, j, pass

        best_length2 = 0.0_dp
        best = (0.0_dp,0.0_dp)
        do i = 1, size(vector)
            trial = (0.0_dp,0.0_dp)
            trial(i) = (1.0_dp,0.0_dp)
            do pass = 1, 2
                do j = 1, size(existing,2)
                    trial = trial - sum(conjg(existing(:,j))*trial)*existing(:,j)
                end do
            end do
            length2 = sum(abs(trial)**2)
            if (length2 > best_length2) then
                best_length2 = length2
                best = trial
            end if
        end do
        ok = best_length2 > 256.0_dp*epsilon(1.0_dp)
        vector = (0.0_dp,0.0_dp)
        if (ok) vector = best/sqrt(best_length2)
    end subroutine orthogonal_complement

    ! Diagonalize a small Hermitian Gram matrix with complex unitary Jacobi
    ! rotations. Accumulate eigenvectors only when explicitly requested.
    pure subroutine hermitian_eigenvalues(a,values,ok,vectors)
        complex(dp), intent(in) :: a(:,:) !! square Hermitian positive-semidefinite matrix
        real(dp), intent(out) :: values(:) !! descending nonnegative eigenvalues
        logical, intent(out) :: ok !! converged, physical eigenvalues
        complex(dp), intent(out), optional :: vectors(:,:) !! unitary eigenvectors, descending
        complex(dp) :: h(size(a,1),size(a,1)), phase, hp, hq, evp, evq
        complex(dp) :: column(size(a,1))
        real(dp) :: scale, tol, b, tau, tangent, cosine, sine, eig
        integer :: n, sweep, i, p, q
        logical :: converged

        n = size(a,1)
        ok = .false.
        values = 0.0_dp
        if (n == 0) return
        h = 0.5_dp*(a+transpose(conjg(a)))
        if (present(vectors)) then
            vectors = (0.0_dp,0.0_dp)
            do i = 1, n
                vectors(i,i) = (1.0_dp,0.0_dp)
            end do
        end if
        scale = maxval(abs(h))
        if (.not. ieee_is_finite(scale)) return
        if (scale == 0.0_dp) then
            ok = .true.
            return
        end if
        tol = 16.0_dp*epsilon(1.0_dp)*scale
        converged = (n == 1)
        do sweep = 1, 80
            if (n == 1) exit
            converged = .true.
            do p = 1, n-1
                do q = p+1, n
                    b = abs(h(p,q))
                    if (b <= tol) cycle
                    converged = .false.
                    tau = 0.5_dp*(real(h(q,q),dp)-real(h(p,p),dp))/b
                    tangent = sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau))
                    cosine = 1.0_dp/sqrt(1.0_dp+tangent*tangent)
                    sine = tangent*cosine
                    phase = h(p,q)/b
                    do i = 1, n
                        if (i == p .or. i == q) cycle
                        hp = h(i,p)
                        hq = h(i,q)
                        h(i,p) = cosine*hp-sine*conjg(phase)*hq
                        h(i,q) = sine*phase*hp+cosine*hq
                        h(p,i) = conjg(h(i,p))
                        h(q,i) = conjg(h(i,q))
                    end do
                    h(p,p) = cmplx(real(h(p,p),dp)-tangent*b,0.0_dp,dp)
                    h(q,q) = cmplx(real(h(q,q),dp)+tangent*b,0.0_dp,dp)
                    h(p,q) = (0.0_dp,0.0_dp)
                    h(q,p) = (0.0_dp,0.0_dp)
                    if (present(vectors)) then
                        do i = 1, n
                            evp = vectors(i,p)
                            evq = vectors(i,q)
                            vectors(i,p) = cosine*evp-sine*conjg(phase)*evq
                            vectors(i,q) = sine*phase*evp+cosine*evq
                        end do
                    end if
                end do
            end do
            if (converged) exit
        end do
        if (.not. converged) return
        do i = 1, n
            eig = real(h(i,i),dp)
            if (.not. ieee_is_finite(eig)) return
            if (eig < -1.0e-8_dp .or. eig > 1.0_dp+1.0e-8_dp) return
            values(i) = min(1.0_dp,max(0.0_dp,eig))
        end do
        ! Number of channels per group is typically small; simple descending
        ! insertion sort avoids an additional dependency.
        do i = 2, n
            eig = values(i)
            if (present(vectors)) column = vectors(:,i)
            p = i-1
            do while (p >= 1)
                if (values(p) >= eig) exit
                values(p+1) = values(p)
                if (present(vectors)) vectors(:,p+1) = vectors(:,p)
                p = p-1
            end do
            values(p+1) = eig
            if (present(vectors)) vectors(:,p+1) = column
        end do
        ok = .true.
    end subroutine hermitian_eigenvalues

    ! Apply the same global diagonal shrinkage/trace-scaled ridge as group
    ! coherence, residualize A and B against C, then whiten on both sides.
    pure subroutine canonical_cell(s,a,b,c,ridge,alpha,modes,ok,a_coeff,b_coeff,stable,metric_a,metric_b)
        complex(dp), intent(in) :: s(:,:) !! unregularized CSD matrix (all channels)
        integer, intent(in) :: a(:) !! indices of group A
        integer, intent(in) :: b(:) !! indices of group B, disjoint from A
        integer, intent(in) :: c(:) !! indices to condition on, disjoint from A/B
        real(dp), intent(in) :: ridge !! fractional trace diagonal loading
        real(dp), intent(in) :: alpha !! diagonal-target shrinkage coefficient
        real(dp), intent(out) :: modes(:) !! descending squared canonical correlations
        logical, intent(out) :: ok !! positive definite blocks and converged eigensolver
        complex(dp), intent(out), optional :: a_coeff(:,:) !! (A-channel,mode), variance normalized
        complex(dp), intent(out), optional :: b_coeff(:,:) !! (B-channel,mode), variance normalized
        logical, intent(out), optional :: stable(:) !! isolated nonzero singular directions
        complex(dp), intent(out), optional :: metric_a(:,:) !! conditioned A autospectrum
        complex(dp), intent(out), optional :: metric_b(:,:) !! conditioned B autospectrum
        complex(dp) :: reg(size(s,1),size(s,2))
        complex(dp) :: la(size(a),size(a)), lb(size(b),size(b)), lc(size(c),size(c))
        complex(dp) :: ca(size(a),size(a)), cb(size(b),size(b))
        complex(dp) :: cross(size(a),size(b)), work(size(a),size(b))
        complex(dp) :: sca(size(c),size(a)), scb(size(c),size(b))
        complex(dp) :: xa(size(c),size(a)), xb(size(c),size(b))
        complex(dp) :: rhs(size(b),size(a)), solved(size(b),size(a))
        complex(dp) :: gram(size(modes),size(modes))
        complex(dp) :: eigenvectors(size(modes),size(modes))
        complex(dp) :: left(size(a),size(modes)), right(size(b),size(modes))
        complex(dp) :: ac(size(a),size(modes)), bc(size(b),size(modes)), gauge
        real(dp) :: mean_power, sigma, gap_threshold
        integer :: i, n, k, dominant
        logical :: completed

        ok = .false.
        modes = 0.0_dp
        n = size(s,1)
        mean_power = 0.0_dp
        do i = 1, n
            mean_power = mean_power + max(0.0_dp,real(s(i,i),dp))/real(n,dp)
        end do
        reg = (1.0_dp-alpha)*0.5_dp*(s+transpose(conjg(s)))
        do i = 1, n
            reg(i,i) = cmplx(real(s(i,i),dp)+ridge*mean_power,0.0_dp,dp)
        end do
        ca = reg(a,a)
        cb = reg(b,b)
        cross = reg(a,b)
        if (size(c) > 0) then
            call cholesky_block(reg(c,c),lc,ok)
            if (.not. ok) return
            sca = reg(c,a)
            scb = reg(c,b)
            call lower_solve(lc,sca,xa)
            call lower_solve(lc,scb,xb)
            ca = ca - matmul(transpose(conjg(xa)),xa)
            cb = cb - matmul(transpose(conjg(xb)),xb)
            cross = cross - matmul(transpose(conjg(xa)),xb)
        end if
        call cholesky_block(ca,la,ok)
        if (.not. ok) return
        call cholesky_block(cb,lb,ok)
        if (.not. ok) return
        if (present(metric_a)) metric_a = ca
        if (present(metric_b)) metric_b = cb
        call lower_solve(la,cross,work)
        rhs = transpose(conjg(work))
        call lower_solve(lb,rhs,solved)
        work = transpose(conjg(solved))
        if (size(a) <= size(b)) then
            gram = matmul(work,transpose(conjg(work)))
        else
            gram = matmul(transpose(conjg(work)),work)
        end if
        if (.not. present(a_coeff)) then
            call hermitian_eigenvalues(gram,modes,ok)
            return
        end if
        call hermitian_eigenvalues(gram,modes,ok,eigenvectors)
        if (.not. ok) return
        left = (0.0_dp,0.0_dp)
        right = (0.0_dp,0.0_dp)
        if (size(a) <= size(b)) then
            left = eigenvectors
        else
            right = eigenvectors
        end if
        ! Nonzero singular directions are paired through W, with positive
        ! real cross covariance. For zero modes choose an orthogonal complement.
        do k = 1, size(modes)
            sigma = sqrt(modes(k))
            if (modes(k) > 64.0_dp*epsilon(1.0_dp)) then
                if (size(a) <= size(b)) then
                    right(:,k) = matmul(transpose(conjg(work)),left(:,k))/sigma
                else
                    left(:,k) = matmul(work,right(:,k))/sigma
                end if
            else
                if (size(a) <= size(b)) then
                    call orthogonal_complement(right(:,:k-1),right(:,k),completed)
                else
                    call orthogonal_complement(left(:,:k-1),left(:,k),completed)
                end if
                if (.not. completed) then
                    ok = .false.
                    return
                end if
            end if
        end do
        call upper_conjugate_solve(la,left,ac)
        call upper_conjugate_solve(lb,right,bc)
        ! Fix common phase without changing the positive cross covariance.
        ! Choose the largest-magnitude coefficient of group A as positive real.
        do k = 1, size(modes)
            dominant = maxloc(abs(ac(:,k)),dim=1)
            if (abs(ac(dominant,k)) > 0.0_dp) then
                gauge = conjg(ac(dominant,k))/abs(ac(dominant,k))
                ac(:,k) = gauge*ac(:,k)
                bc(:,k) = gauge*bc(:,k)
            end if
        end do
        a_coeff = ac
        b_coeff = bc
        if (present(stable)) then
            stable = .true.
            gap_threshold = 64.0_dp*sqrt(epsilon(1.0_dp))
            do k = 1, size(modes)
                if (sqrt(modes(k)) <= gap_threshold) stable(k) = .false.
                if (k > 1) then
                    if (abs(modes(k-1)-modes(k)) <= gap_threshold) stable(k) = .false.
                end if
                if (k < size(modes)) then
                    if (abs(modes(k+1)-modes(k)) <= gap_threshold) stable(k) = .false.
                end if
            end do
        end if
    end subroutine canonical_cell

    function canonical_coherence_matrix(spectrum,group_a,group_b,conditioned_on,ridge,shrinkage, &
            auto_shrinkage,effective_tapers,cv_shrinkage,return_vectors,return_covariances) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! CSD(frequency,channel,channel)
        integer, intent(in) :: group_a(:) !! first nonempty channel group, one-based
        integer, intent(in) :: group_b(:) !! second nonempty disjoint group, one-based
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels
        real(dp), intent(in), optional :: ridge !! relative trace loading, >=0
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage [0,1]
        logical, intent(in), optional :: auto_shrinkage !! Gaussian-Wishart plug-in shrinkage
        real(dp), intent(in), optional :: effective_tapers !! fixed-weight effective taper count
        logical, intent(in), optional :: cv_shrinkage !! precomputed predictive shrinkage
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_result) :: out
        type(group_coherence_result) :: reference
        integer, allocatable :: c(:)
        integer :: f, nf, nm
        real(dp) :: lambda, nan
        logical :: ok, want_vectors, want_metrics

        allocate(out%frequencies(0),out%squared_modes(0,0),out%coherence(0), &
            out%log_information(0),out%shrinkage_used(0),out%valid(0))
        allocate(c(0))
        if (present(conditioned_on)) c = conditioned_on
        reference = group_coherence_matrix(spectrum,group_a,group_b,c,ridge,shrinkage, &
            auto_shrinkage,effective_tapers,cv_shrinkage)
        if (reference%status == partial_coherence_invalid_input) return
        nf = size(reference%frequencies)
        nm = min(size(group_a),size(group_b))
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        deallocate(out%frequencies,out%squared_modes,out%coherence,out%log_information, &
            out%shrinkage_used,out%valid)
        allocate(out%frequencies(nf),out%squared_modes(nf,nm),out%coherence(nf), &
            out%log_information(nf),out%shrinkage_used(nf),out%valid(nf))
        out%frequencies = reference%frequencies
        out%coherence = reference%coherence
        out%log_information = reference%log_information
        out%shrinkage_used = reference%shrinkage_used
        out%valid = reference%valid
        out%status = reference%status
        out%squared_modes = nan
        want_vectors = .false.
        if (present(return_vectors)) want_vectors = return_vectors
        if (want_vectors) then
            allocate(out%coefficients_a(nf,size(group_a),nm),out%coefficients_b(nf,size(group_b),nm), &
                out%directions_stable(nf,nm))
            out%coefficients_a = cmplx(nan,nan,dp)
            out%coefficients_b = cmplx(nan,nan,dp)
            out%directions_stable = .false.
        end if
        want_metrics = .false.
        if (present(return_covariances)) want_metrics = return_covariances
        if (want_metrics) then
            allocate(out%covariance_a(nf,size(group_a),size(group_a)), &
                out%covariance_b(nf,size(group_b),size(group_b)))
            out%covariance_a = cmplx(nan,nan,dp)
            out%covariance_b = cmplx(nan,nan,dp)
        end if
        lambda = 0.0_dp
        if (present(ridge)) lambda = ridge
        do f = 1, nf
            if (.not. out%valid(f)) cycle
            if (want_vectors .and. want_metrics) then
                call canonical_cell(spectrum%power(f,:,:),group_a,group_b,c,lambda, &
                    out%shrinkage_used(f),out%squared_modes(f,:),ok, &
                    out%coefficients_a(f,:,:),out%coefficients_b(f,:,:),out%directions_stable(f,:), &
                    out%covariance_a(f,:,:),out%covariance_b(f,:,:))
            else if (want_vectors) then
                call canonical_cell(spectrum%power(f,:,:),group_a,group_b,c,lambda, &
                    out%shrinkage_used(f),out%squared_modes(f,:),ok, &
                    out%coefficients_a(f,:,:),out%coefficients_b(f,:,:),out%directions_stable(f,:))
            else if (want_metrics) then
                call canonical_cell(spectrum%power(f,:,:),group_a,group_b,c,lambda, &
                    out%shrinkage_used(f),out%squared_modes(f,:),ok, &
                    metric_a=out%covariance_a(f,:,:),metric_b=out%covariance_b(f,:,:))
            else
                call canonical_cell(spectrum%power(f,:,:),group_a,group_b,c,lambda, &
                    out%shrinkage_used(f),out%squared_modes(f,:),ok)
            end if
            if (ok) cycle
            out%valid(f) = .false.
            out%squared_modes(f,:) = nan
            if (want_vectors) then
                out%coefficients_a(f,:,:) = cmplx(nan,nan,dp)
                out%coefficients_b(f,:,:) = cmplx(nan,nan,dp)
                out%directions_stable(f,:) = .false.
            end if
            if (want_metrics) then
                out%covariance_a(f,:,:) = cmplx(nan,nan,dp)
                out%covariance_b(f,:,:) = cmplx(nan,nan,dp)
            end if
            out%coherence(f) = nan
            out%log_information(f) = nan
            out%status = partial_coherence_singular
        end do
    end function canonical_coherence_matrix

    function canonical_coherence_tf_matrix(spectrum,group_a,group_b,conditioned_on,ridge,shrinkage, &
            auto_shrinkage,effective_tapers,cv_shrinkage,return_vectors,return_covariances) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum !! CSD(f,t,i,j)
        integer, intent(in) :: group_a(:) !! nonempty group A
        integer, intent(in) :: group_b(:) !! nonempty group B, disjoint from A
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning group
        real(dp), intent(in), optional :: ridge !! relative trace loading, >=0
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage [0,1]
        logical, intent(in), optional :: auto_shrinkage !! Gaussian-Wishart plug-in shrinkage
        real(dp), intent(in), optional :: effective_tapers !! effective fixed-weight taper count
        logical, intent(in), optional :: cv_shrinkage !! precomputed predictive shrinkage
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_tf_result) :: out
        type(group_coherence_tf_result) :: reference
        integer, allocatable :: c(:)
        integer :: f, t, nf, nt, nm
        real(dp) :: lambda, nan
        logical :: ok, want_vectors, want_metrics

        allocate(out%frequencies(0),out%times(0),out%squared_modes(0,0,0), &
            out%coherence(0,0),out%log_information(0,0),out%shrinkage_used(0,0),out%valid(0,0))
        allocate(c(0))
        if (present(conditioned_on)) c = conditioned_on
        reference = group_coherence_tf_matrix(spectrum,group_a,group_b,c,ridge,shrinkage, &
            auto_shrinkage,effective_tapers,cv_shrinkage)
        if (reference%status == partial_coherence_invalid_input) return
        nf = size(reference%frequencies)
        nt = size(reference%times)
        nm = min(size(group_a),size(group_b))
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        deallocate(out%frequencies,out%times,out%squared_modes,out%coherence, &
            out%log_information,out%shrinkage_used,out%valid)
        allocate(out%frequencies(nf),out%times(nt),out%squared_modes(nf,nt,nm), &
            out%coherence(nf,nt),out%log_information(nf,nt),out%shrinkage_used(nf,nt), &
            out%valid(nf,nt))
        out%frequencies = reference%frequencies
        out%times = reference%times
        out%coherence = reference%coherence
        out%log_information = reference%log_information
        out%shrinkage_used = reference%shrinkage_used
        out%valid = reference%valid
        out%status = reference%status
        out%squared_modes = nan
        want_vectors = .false.
        if (present(return_vectors)) want_vectors = return_vectors
        if (want_vectors) then
            allocate(out%coefficients_a(nf,nt,size(group_a),nm), &
                out%coefficients_b(nf,nt,size(group_b),nm),out%directions_stable(nf,nt,nm))
            out%coefficients_a = cmplx(nan,nan,dp)
            out%coefficients_b = cmplx(nan,nan,dp)
            out%directions_stable = .false.
        end if
        want_metrics = .false.
        if (present(return_covariances)) want_metrics = return_covariances
        if (want_metrics) then
            allocate(out%covariance_a(nf,nt,size(group_a),size(group_a)), &
                out%covariance_b(nf,nt,size(group_b),size(group_b)))
            out%covariance_a = cmplx(nan,nan,dp)
            out%covariance_b = cmplx(nan,nan,dp)
        end if
        lambda = 0.0_dp
        if (present(ridge)) lambda = ridge
        do t = 1, nt
            do f = 1, nf
                if (.not. out%valid(f,t)) cycle
                if (want_vectors .and. want_metrics) then
                    call canonical_cell(spectrum%power(f,t,:,:),group_a,group_b,c,lambda, &
                        out%shrinkage_used(f,t),out%squared_modes(f,t,:),ok, &
                        out%coefficients_a(f,t,:,:),out%coefficients_b(f,t,:,:), &
                        out%directions_stable(f,t,:), &
                        out%covariance_a(f,t,:,:),out%covariance_b(f,t,:,:))
                else if (want_vectors) then
                    call canonical_cell(spectrum%power(f,t,:,:),group_a,group_b,c,lambda, &
                        out%shrinkage_used(f,t),out%squared_modes(f,t,:),ok, &
                        out%coefficients_a(f,t,:,:),out%coefficients_b(f,t,:,:), &
                        out%directions_stable(f,t,:))
                else if (want_metrics) then
                    call canonical_cell(spectrum%power(f,t,:,:),group_a,group_b,c,lambda, &
                        out%shrinkage_used(f,t),out%squared_modes(f,t,:),ok, &
                        metric_a=out%covariance_a(f,t,:,:),metric_b=out%covariance_b(f,t,:,:))
                else
                    call canonical_cell(spectrum%power(f,t,:,:),group_a,group_b,c,lambda, &
                        out%shrinkage_used(f,t),out%squared_modes(f,t,:),ok)
                end if
                if (ok) cycle
                out%valid(f,t) = .false.
                out%squared_modes(f,t,:) = nan
                if (want_vectors) then
                    out%coefficients_a(f,t,:,:) = cmplx(nan,nan,dp)
                    out%coefficients_b(f,t,:,:) = cmplx(nan,nan,dp)
                    out%directions_stable(f,t,:) = .false.
                end if
                if (want_metrics) then
                    out%covariance_a(f,t,:,:) = cmplx(nan,nan,dp)
                    out%covariance_b(f,t,:,:) = cmplx(nan,nan,dp)
                end if
                out%coherence(f,t) = nan
                out%log_information(f,t) = nan
                out%status = partial_coherence_singular
            end do
        end do
    end function canonical_coherence_tf_matrix

    function band_canonical_coherence_matrix(spectrum,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage,return_vectors,return_covariances) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! frequency-resolved CSD matrix
        real(dp), intent(in) :: bands(:,:) !! band edges in Hz, shape (2,nband)
        integer, intent(in) :: group_a(:) !! nonempty group A
        integer, intent(in) :: group_b(:) !! disjoint nonempty group B
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning group
        logical, intent(in), optional :: average !! mean density rather than total band power
        real(dp), intent(in), optional :: ridge !! relative trace loading after integration
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage after integration
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_result) :: out
        type(multitaper_matrix_result) :: integrated

        integrated = integrate_spectral_bands(spectrum,bands,average)
        out = canonical_coherence_matrix(integrated,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,return_vectors=return_vectors,return_covariances=return_covariances)
    end function band_canonical_coherence_matrix

    function band_canonical_coherence_tf_matrix(spectrum,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage,return_vectors,return_covariances) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum !! CSD(f,t,i,j)
        real(dp), intent(in) :: bands(:,:) !! band edges in Hz, shape (2,nband)
        integer, intent(in) :: group_a(:) !! nonempty group A
        integer, intent(in) :: group_b(:) !! disjoint nonempty group B
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning group
        logical, intent(in), optional :: average !! mean density rather than band power
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! fixed shrinkage after integration
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: integrated

        integrated = integrate_tf_spectral_bands(spectrum,bands,average)
        out = canonical_coherence_tf_matrix(integrated,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,return_vectors=return_vectors,return_covariances=return_covariances)
    end function band_canonical_coherence_tf_matrix
    function multitaper_canonical_coherence_real(x,nw,group_a,group_b,conditioned_on,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage,auto_shrinkage, &
            cv_shrinkage,return_vectors,return_covariances) result(out)
        real(dp), intent(in) :: x(:,:) !! real observations (sample,channel) by default
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: group_a(:) !! first group indices
        integer, intent(in) :: group_b(:) !! second group indices
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean before tapering
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! one-based sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative trace ridge
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        logical, intent(in), optional :: auto_shrinkage !! fixed-taper automatic shrinkage
        logical, intent(in), optional :: cv_shrinkage !! LOO shrinkage including adaptive weights
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_result) :: out
        type(multitaper_matrix_result) :: spec

        spec = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,cv_shrinkage)
        out = canonical_coherence_matrix(spec,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,auto_shrinkage=auto_shrinkage,cv_shrinkage=cv_shrinkage, &
            return_vectors=return_vectors,return_covariances=return_covariances)
    end function multitaper_canonical_coherence_real

    function multitaper_canonical_coherence_complex(x,nw,group_a,group_b,conditioned_on,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage,auto_shrinkage, &
            cv_shrinkage,return_vectors,return_covariances) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex observations
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: group_a(:) !! first group indices
        integer, intent(in) :: group_b(:) !! second group indices
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean before tapering
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! one-based sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative trace ridge
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        logical, intent(in), optional :: auto_shrinkage !! fixed-taper automatic shrinkage
        logical, intent(in), optional :: cv_shrinkage !! LOO shrinkage including adaptive weights
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_result) :: out
        type(multitaper_matrix_result) :: spec

        spec = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,cv_shrinkage)
        out = canonical_coherence_matrix(spec,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,auto_shrinkage=auto_shrinkage,cv_shrinkage=cv_shrinkage, &
            return_vectors=return_vectors,return_covariances=return_covariances)
    end function multitaper_canonical_coherence_complex

    function multitaper_tf_canonical_coherence_real(x,nw,nperseg,group_a,group_b,conditioned_on, &
            noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            auto_shrinkage,cv_shrinkage,return_vectors,return_covariances) result(out)
        real(dp), intent(in) :: x(:,:) !! real samples by channels
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples per window
        integer, intent(in) :: group_a(:) !! first channel group
        integer, intent(in) :: group_b(:) !! second channel group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! Fourier length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis (one-based)
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative loading
        real(dp), intent(in), optional :: shrinkage !! fixed shrinkage
        logical, intent(in), optional :: auto_shrinkage !! fixed-taper automatic shrinkage
        logical, intent(in), optional :: cv_shrinkage !! adaptive-capable predictive shrinkage
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spec

        spec = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode,cv_shrinkage)
        out = canonical_coherence_tf_matrix(spec,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,auto_shrinkage=auto_shrinkage,cv_shrinkage=cv_shrinkage, &
            return_vectors=return_vectors,return_covariances=return_covariances)
    end function multitaper_tf_canonical_coherence_real

    function multitaper_tf_canonical_coherence_complex(x,nw,nperseg,group_a,group_b,conditioned_on, &
            noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            auto_shrinkage,cv_shrinkage,return_vectors,return_covariances) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex samples by channels
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples per window
        integer, intent(in) :: group_a(:) !! first channel group
        integer, intent(in) :: group_b(:) !! second channel group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! Fourier length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis (one-based)
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative loading
        real(dp), intent(in), optional :: shrinkage !! fixed shrinkage
        logical, intent(in), optional :: auto_shrinkage !! fixed-taper automatic shrinkage
        logical, intent(in), optional :: cv_shrinkage !! predictive shrinkage
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spec

        spec = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode,cv_shrinkage)
        out = canonical_coherence_tf_matrix(spec,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,auto_shrinkage=auto_shrinkage,cv_shrinkage=cv_shrinkage, &
            return_vectors=return_vectors,return_covariances=return_covariances)
    end function multitaper_tf_canonical_coherence_complex

    function multitaper_band_canonical_coherence_real(x,nw,bands,group_a,group_b,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage,return_vectors,return_covariances) result(out)
        real(dp), intent(in) :: x(:,:) !! real signals
        real(dp), intent(in) :: nw !! time-bandwidth product
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband)
        integer, intent(in) :: group_a(:) !! first channel group
        integer, intent(in) :: group_b(:) !! second channel group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! Fourier length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: average !! mean PSD per band
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! shrinkage after integration
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_result) :: out
        type(multitaper_matrix_result) :: spec

        spec = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_canonical_coherence_matrix(spec,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage,return_vectors,return_covariances)
    end function multitaper_band_canonical_coherence_real

    function multitaper_band_canonical_coherence_complex(x,nw,bands,group_a,group_b,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage,return_vectors,return_covariances) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex signals
        real(dp), intent(in) :: nw !! time-bandwidth product
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband)
        integer, intent(in) :: group_a(:) !! first channel group
        integer, intent(in) :: group_b(:) !! second channel group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! Fourier length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: average !! mean PSD per band
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! shrinkage after integration
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_result) :: out
        type(multitaper_matrix_result) :: spec

        spec = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_canonical_coherence_matrix(spec,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage,return_vectors,return_covariances)
    end function multitaper_band_canonical_coherence_complex

    function multitaper_tf_band_canonical_coherence_real(x,nw,nperseg,bands,group_a,group_b,conditioned_on, &
            noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,average,ridge, &
            shrinkage,return_vectors,return_covariances) result(out)
        real(dp), intent(in) :: x(:,:) !! real signals
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! window sample count
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband)
        integer, intent(in) :: group_a(:) !! first group
        integer, intent(in) :: group_b(:) !! second group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: noverlap !! sample overlap
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: average !! mean band spectral density
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! shrinkage after integration
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spec

        spec = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode)
        out = band_canonical_coherence_tf_matrix(spec,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage,return_vectors,return_covariances)
    end function multitaper_tf_band_canonical_coherence_real

    function multitaper_tf_band_canonical_coherence_complex(x,nw,nperseg,bands,group_a,group_b,conditioned_on, &
            noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,average,ridge, &
            shrinkage,return_vectors,return_covariances) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex signals
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! window sample count
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband)
        integer, intent(in) :: group_a(:) !! first group
        integer, intent(in) :: group_b(:) !! second group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: noverlap !! sample overlap
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: average !! mean band spectral density
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! shrinkage after integration
        logical, intent(in), optional :: return_vectors !! allocate variance-normalized mode directions
        logical, intent(in), optional :: return_covariances !! return conditional A/B spectral metrics
        type(canonical_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spec

        spec = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode)
        out = band_canonical_coherence_tf_matrix(spec,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage,return_vectors,return_covariances)
    end function multitaper_tf_band_canonical_coherence_complex

end module scifort_signal_canonical_coherence
