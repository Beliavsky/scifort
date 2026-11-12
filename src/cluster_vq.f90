! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent Euclidean vector quantization and Lloyd clustering algorithms.
! The conventions follow scipy.cluster.vq, with one-based Fortran labels.
! See CODE_PROVENANCE.md; no third-party algorithm source is included.
module scifort_cluster_vq
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: cluster_vq_success = 0
    integer, parameter, public :: cluster_vq_invalid_input = 1
    integer, parameter, public :: cluster_vq_empty_cluster = 2

    public :: whiten, vq, kmeans2, kmeans

    interface whiten
        module procedure whiten_matrix
        module procedure whiten_vector
    end interface whiten

    interface vq
        module procedure vq_matrix
        module procedure vq_vector
    end interface vq

    interface kmeans2
        module procedure kmeans2_count
        module procedure kmeans2_matrix
        module procedure kmeans2_vector_count
        module procedure kmeans2_vector_initial
    end interface kmeans2

    interface kmeans
        module procedure kmeans_count
        module procedure kmeans_matrix
        module procedure kmeans_vector_count
        module procedure kmeans_vector_initial
    end interface kmeans

contains

    function whiten_matrix(obs, status) result(result_matrix)
        real(dp), intent(in) :: obs(:,:) !! Observations by rows; finite, nonempty.
        integer, intent(out), optional :: status !! Zero on success, one for invalid observations.
        real(dp), allocatable :: result_matrix(:,:)
        real(dp) :: scale, center, variance
        integer :: j, n

        if (present(status)) status = cluster_vq_invalid_input
        allocate(result_matrix(0,0))
        n = size(obs,1)
        if (n < 1 .or. size(obs,2) < 1) return
        if (.not. all(ieee_is_finite(obs))) return
        deallocate(result_matrix)
        allocate(result_matrix(n,size(obs,2)))
        do j = 1, size(obs,2)
            scale = maxval(abs(obs(:,j)))
            if (scale == 0.0_dp) then
                result_matrix(:,j) = obs(:,j)
                cycle
            end if
            center = sum(obs(:,j)/scale)/real(n,dp)
            variance = sum((obs(:,j)/scale-center)**2)/real(n,dp)
            if (variance <= 0.0_dp) then
                ! SciPy leaves constant features unchanged instead of dividing by zero.
                result_matrix(:,j) = obs(:,j)
            else
                result_matrix(:,j) = (obs(:,j)/scale)/sqrt(variance)
            end if
        end do
        if (present(status)) status = cluster_vq_success
    end function whiten_matrix

    function whiten_vector(obs, status) result(result_vector)
        real(dp), intent(in) :: obs(:) !! One-dimensional observations, at least one finite value.
        integer, intent(out), optional :: status !! Zero on success, one for invalid observations.
        real(dp), allocatable :: result_vector(:)
        real(dp), allocatable :: values(:,:)
        integer :: code

        values = whiten_matrix(reshape(obs,[size(obs),1]), code)
        allocate(result_vector(size(values,1)))
        if (size(result_vector)>0) result_vector = values(:,1)
        if (present(status)) status = code
    end function whiten_vector

    pure function euclidean_distance(x, y) result(distance)
        real(dp), intent(in) :: x(:) !! First feature vector.
        real(dp), intent(in) :: y(:) !! Second feature vector of the same length.
        real(dp) :: distance, scale

        ! Scaling before subtraction/squaring prevents intermediate overflow.
        scale = max(maxval(abs(x)),maxval(abs(y)))
        if (scale == 0.0_dp) then
            distance = 0.0_dp
        else
            distance = scale*sqrt(sum((x/scale-y/scale)**2))
        end if
    end function euclidean_distance

    subroutine vq_matrix(obs, codebook, labels, distances, status)
        real(dp), intent(in) :: obs(:,:) !! Finite observations, shape (n, features).
        real(dp), intent(in) :: codebook(:,:) !! Finite centroids, shape (k, features), k >= 1.
        integer, allocatable, intent(out) :: labels(:) !! One-based nearest-centroid indices.
        real(dp), allocatable, intent(out) :: distances(:) !! Euclidean distance to the nearest centroid.
        integer, intent(out), optional :: status !! Zero on success, one for invalid input.
        real(dp) :: d, best
        integer :: i, j, nearest

        if (present(status)) status = cluster_vq_invalid_input
        allocate(labels(0),distances(0))
        if (size(obs,1) < 1 .or. size(obs,2) < 1 .or. size(codebook,1) < 1) return
        if (size(obs,2) /= size(codebook,2)) return
        if (.not. all(ieee_is_finite(obs)) .or. .not. all(ieee_is_finite(codebook))) return
        deallocate(labels,distances)
        allocate(labels(size(obs,1)),distances(size(obs,1)))
        do i = 1, size(obs,1)
            best = euclidean_distance(obs(i,:),codebook(1,:))
            nearest = 1
            do j = 2, size(codebook,1)
                d = euclidean_distance(obs(i,:),codebook(j,:))
                if (d < best) then
                    best = d
                    nearest = j
                end if
            end do
            labels(i) = nearest
            distances(i) = best
        end do
        if (present(status)) status = cluster_vq_success
    end subroutine vq_matrix

    subroutine vq_vector(obs, codebook, labels, distances, status)
        real(dp), intent(in) :: obs(:) !! Finite one-dimensional observations.
        real(dp), intent(in) :: codebook(:) !! Nonempty vector of candidate centroids.
        integer, allocatable, intent(out) :: labels(:) !! One-based nearest-centroid indices.
        real(dp), allocatable, intent(out) :: distances(:) !! Absolute distance to assigned centroid.
        integer, intent(out), optional :: status !! Zero on success, one for invalid input.
        integer :: code

        call vq_matrix(reshape(obs,[size(obs),1]), reshape(codebook,[size(codebook),1]), &
            labels, distances, code)
        if (present(status)) status = code
    end subroutine vq_vector

    subroutine kmeans2_vector_count(obs,k,centroids,labels,status,state,iterations,init,missing)
        real(dp), intent(in) :: obs(:) !! Finite one-dimensional observations.
        integer, intent(in) :: k !! Requested number of centers, 1 <= k <= n.
        real(dp), allocatable, intent(out) :: centroids(:) !! Final one-dimensional centers.
        integer, allocatable, intent(out) :: labels(:) !! One-based final assignments.
        integer, intent(out) :: status !! Zero on success or diagnostic status.
        type(rng_state), intent(inout), optional :: state !! Explicit RNG; required for count init.
        integer, intent(in), optional :: iterations !! Lloyd iteration count, default 10.
        character(len=*), intent(in), optional :: init !! 'points', '++', or 'random' (Gaussian).
        character(len=*), intent(in), optional :: missing !! 'warn' or 'raise'.
        real(dp), allocatable :: matrix(:,:)

        call kmeans2_count(reshape(obs,[size(obs),1]),k,matrix,labels,status, &
            state,iterations,init,missing)
        allocate(centroids(size(matrix,1)))
        if (size(centroids)>0) centroids=matrix(:,1)
    end subroutine kmeans2_vector_count

    subroutine kmeans2_vector_initial(obs,initial,centroids,labels,status,iterations,missing)
        real(dp), intent(in) :: obs(:) !! Finite one-dimensional observations.
        real(dp), intent(in) :: initial(:) !! Initial one-dimensional centers.
        real(dp), allocatable, intent(out) :: centroids(:) !! Final one-dimensional centers.
        integer, allocatable, intent(out) :: labels(:) !! One-based final assignments.
        integer, intent(out) :: status !! Zero on success or diagnostic status.
        integer, intent(in), optional :: iterations !! Lloyd iteration count, default 10.
        character(len=*), intent(in), optional :: missing !! 'warn' or 'raise'.
        real(dp), allocatable :: matrix(:,:)

        call kmeans2_matrix(reshape(obs,[size(obs),1]), &
            reshape(initial,[size(initial),1]),matrix,labels,status,iterations,missing)
        allocate(centroids(size(matrix,1)))
        if (size(centroids)>0) centroids=matrix(:,1)
    end subroutine kmeans2_vector_initial

    subroutine kmeans2_count(obs, k, centroids, labels, status, state, iterations, init, missing)
        real(dp), intent(in) :: obs(:,:) !! Finite observation matrix, (n, features).
        integer, intent(in) :: k !! Number of clusters, 1 <= k <= n.
        real(dp), allocatable, intent(out) :: centroids(:,:) !! Final k centers.
        integer, allocatable, intent(out) :: labels(:) !! One-based final nearest-center labels.
        integer, intent(out) :: status !! Zero on success; invalid input or empty cluster otherwise.
        type(rng_state), intent(inout), optional :: state !! Explicit RNG, required for initialization.
        integer, intent(in), optional :: iterations !! Lloyd iterations, default 10, >= 1.
        character(len=*), intent(in), optional :: init !! 'points' (default), '++', or 'random'.
        character(len=*), intent(in), optional :: missing !! 'warn' (default) or 'raise'.
        character(len=:), allocatable :: method
        real(dp), allocatable :: initial(:,:)
        integer :: code

        status = cluster_vq_invalid_input
        allocate(centroids(0,0),labels(0))
        if (.not. valid_observations(obs)) return
        if (k < 1 .or. k > size(obs,1)) return
        if (.not. present(state)) return
        method = 'points'
        if (present(init)) method = trim(init)
        if (method /= 'points' .and. method /= '++' .and. method /= 'random') return
        if (.not. valid_options(iterations,missing)) return
        call initialize_centers(obs,k,state,method,initial,code)
        if (code /= cluster_vq_success) return
        call run_kmeans2(obs,initial,centroids,labels,status,iterations,missing)
    end subroutine kmeans2_count

    subroutine kmeans2_matrix(obs, initial, centroids, labels, status, iterations, missing)
        real(dp), intent(in) :: obs(:,:) !! Finite observation matrix, (n, features).
        real(dp), intent(in) :: initial(:,:) !! Initial centroids, shape (k, features), k >= 1.
        real(dp), allocatable, intent(out) :: centroids(:,:) !! Final k centers.
        integer, allocatable, intent(out) :: labels(:) !! One-based final nearest-center labels.
        integer, intent(out) :: status !! Zero on success; invalid input or empty cluster otherwise.
        integer, intent(in), optional :: iterations !! Lloyd iterations, default 10, >= 1.
        character(len=*), intent(in), optional :: missing !! 'warn' (default) or 'raise'.

        status = cluster_vq_invalid_input
        allocate(centroids(0,0),labels(0))
        if (.not. valid_centers(obs,initial)) return
        if (.not. valid_options(iterations,missing)) return
        deallocate(centroids,labels)
        call run_kmeans2(obs,initial,centroids,labels,status,iterations,missing)
    end subroutine kmeans2_matrix

    pure logical function valid_observations(obs) result(valid)
        real(dp), intent(in) :: obs(:,:) !! Observation matrix to validate.
        valid = size(obs,1)>0 .and. size(obs,2)>0
        if (valid) valid = all(ieee_is_finite(obs))
    end function valid_observations

    pure logical function valid_centers(obs,initial) result(valid)
        real(dp), intent(in) :: obs(:,:) !! Observation matrix to validate.
        real(dp), intent(in) :: initial(:,:) !! Proposed initial centers.
        valid = valid_observations(obs) .and. size(initial,1)>0
        valid = valid .and. size(initial,2)==size(obs,2)
        if (valid) valid = all(ieee_is_finite(initial))
    end function valid_centers

    pure logical function valid_options(iterations, missing) result(valid)
        integer, intent(in), optional :: iterations !! Requested iteration count.
        character(len=*), intent(in), optional :: missing !! Empty-cluster policy name.
        valid = .true.
        if (present(iterations)) valid = iterations>=1
        if (present(missing)) then
            valid = valid .and. (missing=='warn' .or. missing=='raise')
        end if
    end function valid_options

    subroutine run_kmeans2(obs, initial, centroids, labels, status, iterations, missing)
        real(dp), intent(in) :: obs(:,:) !! Nonempty finite observation matrix.
        real(dp), intent(in) :: initial(:,:) !! Nonempty finite codebook of starting centers.
        real(dp), allocatable, intent(out) :: centroids(:,:) !! Updated codebook.
        integer, allocatable, intent(out) :: labels(:) !! Final one-based assignments.
        integer, intent(out) :: status !! Success or empty cluster; valid output for 'warn'.
        integer, intent(in), optional :: iterations !! Number of updates, default 10.
        character(len=*), intent(in), optional :: missing !! 'warn' (default) or 'raise'.
        real(dp), allocatable :: new_centers(:,:), distances(:)
        integer, allocatable :: counts(:)
        integer :: nsteps, step, k, i, j, code
        logical :: raise_on_empty, had_empty

        nsteps = 10
        if (present(iterations)) nsteps = iterations
        raise_on_empty = .false.
        if (present(missing)) raise_on_empty = missing=='raise'
        centroids = initial
        k = size(initial,1)
        allocate(new_centers(k,size(obs,2)),counts(k))
        had_empty = .false.
        do step = 1,nsteps
            call vq_matrix(obs,centroids,labels,distances,code)
            counts = 0
            do i = 1,size(labels)
                counts(labels(i)) = counts(labels(i))+1
            end do
            if (any(counts==0)) then
                had_empty = .true.
                if (raise_on_empty) then
                    deallocate(centroids,labels)
                    allocate(centroids(0,0),labels(0))
                    status = cluster_vq_empty_cluster
                    return
                end if
            end if
            new_centers = 0.0_dp
            do i = 1,size(labels)
                j = labels(i)
                new_centers(j,:) = new_centers(j,:) + obs(i,:)/real(counts(j),dp)
            end do
            do j = 1,k
                if (counts(j)==0) new_centers(j,:) = centroids(j,:)
            end do
            centroids = new_centers
        end do
        call vq_matrix(obs,centroids,labels,distances,code)
        status = cluster_vq_success
        if (had_empty) status = cluster_vq_empty_cluster
    end subroutine run_kmeans2

    subroutine initialize_centers(obs,k,state,method,centroids,status)
        real(dp), intent(in) :: obs(:,:) !! Finite observation matrix, (n, features).
        integer, intent(in) :: k !! Requested number of centers, 1 <= k <= n.
        type(rng_state), intent(inout) :: state !! Generator advanced by initialization draws.
        character(len=*), intent(in) :: method !! 'points', '++', or 'random'.
        real(dp), allocatable, intent(out) :: centroids(:,:) !! Selected observation rows.
        integer, intent(out) :: status !! Zero on success.
        integer, allocatable :: order(:), selected(:)
        real(dp), allocatable :: d2(:)
        real(dp) :: best, maxdist, total, target, cumulative, d
        integer :: n, i, j, next, tmp, chosen, t

        n = size(obs,1)
        allocate(centroids(k,size(obs,2)),selected(k))
        status = cluster_vq_success
        if (method=='random') then
            call initialize_gaussian(obs,state,centroids,status)
            return
        end if
        if (method=='points') then
            allocate(order(n))
            order = [(i,i=1,n)]
            do i = 1,k
                next = i + int(rng_uniform(state)*real(n-i+1,dp))
                tmp = order(i)
                order(i) = order(next)
                order(next) = tmp
                centroids(i,:) = obs(order(i),:)
            end do
            return
        end if

        ! K-means++: sample each new point proportional to its squared
        ! Euclidean distance to the closest previously chosen center.
        selected = 0
        selected(1) = 1 + int(rng_uniform(state)*real(n,dp))
        centroids(1,:) = obs(selected(1),:)
        allocate(d2(n))
        do i = 2,k
            maxdist = 0.0_dp
            do j = 1,n
                best = euclidean_distance(obs(j,:),centroids(1,:))
                do t = 2,i-1
                    d = euclidean_distance(obs(j,:),centroids(t,:))
                    if (d<best) best = d
                end do
                d2(j) = best
                maxdist = max(maxdist,best)
            end do
            if (maxdist==0.0_dp) then
                ! All observations have zero residual distance; avoid duplicate
                ! selected indices while allowing equal-valued observations.
                chosen = 0
                do j = 1,n
                    if (.not. any(selected(:i-1)==j)) then
                        chosen = j
                        exit
                    end if
                end do
            else
                d2 = (d2/maxdist)**2
                total = sum(d2)
                target = rng_uniform(state)*total
                cumulative = 0.0_dp
                chosen = n
                do j = 1,n
                    cumulative = cumulative+d2(j)
                    if (cumulative>target) then
                        chosen = j
                        exit
                    end if
                end do
            end if
            selected(i) = chosen
            centroids(i,:) = obs(chosen,:)
        end do
    end subroutine initialize_centers

    subroutine initialize_gaussian(obs, state, centroids, status)
        ! Draw N(mu,C) samples, where C is the unbiased sample covariance.
        ! The projection g^T (X - mu)/sqrt(n-1), g_i ~ N(0,1), has exactly
        ! this covariance. Unlike a Cholesky factorization it handles rank-
        ! deficient data (including n < features) without LAPACK or SVD.
        ! Scale each feature before accumulating to minimize overflow risk.
        real(dp), intent(in) :: obs(:,:) !! Finite rows of input observations.
        type(rng_state), intent(inout) :: state !! Generator advanced by Gaussian draws.
        real(dp), intent(out) :: centroids(:,:) !! Gaussian-initialized centers (k,features).
        integer, intent(out) :: status !! Success or invalid input (n < 2/overflow).
        real(dp), allocatable :: normalized(:,:), means(:), scales(:), g(:)
        real(dp) :: inverse_sd
        integer :: n, p, i, j

        status = cluster_vq_invalid_input
        n = size(obs,1)
        p = size(obs,2)
        if (n < 2 .or. p < 1 .or. size(centroids,2) /= p) return
        allocate(normalized(n,p),means(p),scales(p),g(n))
        do j = 1,p
            scales(j) = maxval(abs(obs(:,j)))
            if (scales(j) > 0.0_dp) then
                normalized(:,j) = obs(:,j)/scales(j)
            else
                normalized(:,j) = 0.0_dp
            end if
            means(j) = sum(normalized(:,j))/real(n,dp)
            normalized(:,j) = normalized(:,j)-means(j)
        end do
        inverse_sd = 1.0_dp/sqrt(real(n-1,dp))
        do i = 1,size(centroids,1)
            do j = 1,n
                g(j) = gaussian_standard(state)
            end do
            centroids(i,:) = (means + matmul(g,normalized)*inverse_sd)*scales
        end do
        if (.not. all(ieee_is_finite(centroids))) return
        status = cluster_vq_success
    end subroutine initialize_gaussian

    function gaussian_standard(state) result(z)
        ! Polar Box-Muller, using the library's open-interval uniform RNG.
        type(rng_state), intent(inout) :: state !! Explicit RNG state.
        real(dp) :: z, u1, u2
        real(dp), parameter :: TWO_PI = 2.0_dp*acos(-1.0_dp)

        ! Two explicit statements define RNG consumption even when a
        ! compiler reorders expression evaluation under optimization.
        u1 = rng_uniform(state)
        u2 = rng_uniform(state)
        z = sqrt(-2.0_dp*log(u1))*cos(TWO_PI*u2)
    end function gaussian_standard

    subroutine kmeans_vector_count(obs,k,centroids,distortion,status,state, &
        restarts,threshold,max_iterations)
        real(dp), intent(in) :: obs(:) !! Finite one-dimensional observations.
        integer, intent(in) :: k !! Requested number of centers, 1 <= k <= n.
        real(dp), allocatable, intent(out) :: centroids(:) !! Best nonempty one-dimensional centers.
        real(dp), intent(out) :: distortion !! Mean nonsquared distance, or huge on error.
        integer, intent(out) :: status !! Zero on success or invalid-input status.
        type(rng_state), intent(inout) :: state !! Explicit generator for restarts.
        integer, intent(in), optional :: restarts !! Number of starts, default 20.
        real(dp), intent(in), optional :: threshold !! Convergence threshold, default 1e-5.
        integer, intent(in), optional :: max_iterations !! Max updates per start, default 300.
        real(dp), allocatable :: matrix(:,:)

        call kmeans_count(reshape(obs,[size(obs),1]),k,matrix,distortion,status,state, &
            restarts,threshold,max_iterations)
        allocate(centroids(size(matrix,1)))
        if (size(centroids)>0) centroids=matrix(:,1)
    end subroutine kmeans_vector_count

    subroutine kmeans_vector_initial(obs,initial,centroids,distortion,status, &
        threshold,max_iterations)
        real(dp), intent(in) :: obs(:) !! Finite one-dimensional observations.
        real(dp), intent(in) :: initial(:) !! Initial one-dimensional centers.
        real(dp), allocatable, intent(out) :: centroids(:) !! Nonempty one-dimensional centers.
        real(dp), intent(out) :: distortion !! Mean nonsquared distance, or huge on error.
        integer, intent(out) :: status !! Zero on success or invalid-input status.
        real(dp), intent(in), optional :: threshold !! Convergence threshold, default 1e-5.
        integer, intent(in), optional :: max_iterations !! Max iterations, default 300.
        real(dp), allocatable :: matrix(:,:)

        call kmeans_matrix(reshape(obs,[size(obs),1]), &
            reshape(initial,[size(initial),1]),matrix,distortion,status, &
            threshold,max_iterations)
        allocate(centroids(size(matrix,1)))
        if (size(centroids)>0) centroids=matrix(:,1)
    end subroutine kmeans_vector_initial

    subroutine kmeans_count(obs, k, centroids, distortion, status, state, restarts, &
        threshold, max_iterations)
        real(dp), intent(in) :: obs(:,:) !! Finite observations, (n, features).
        integer, intent(in) :: k !! Requested number of centers, 1 <= k <= n.
        real(dp), allocatable, intent(out) :: centroids(:,:) !! Best nonempty codebook.
        real(dp), intent(out) :: distortion !! Mean unsquared Euclidean distance, or huge on error.
        integer, intent(out) :: status !! Zero on success; one for invalid input.
        type(rng_state), intent(inout) :: state !! Explicit generator for 'points' restarts.
        integer, intent(in), optional :: restarts !! Number of initializations, default 20.
        real(dp), intent(in), optional :: threshold !! Convergence tolerance, default 1e-5.
        integer, intent(in), optional :: max_iterations !! Max updates per restart, default 300.
        real(dp), allocatable :: guess(:,:), candidate(:,:)
        real(dp) :: candidate_dist
        integer :: ntries, t, code

        allocate(centroids(0,0))
        status = cluster_vq_invalid_input
        distortion = huge(1.0_dp)
        if (.not. valid_observations(obs)) return
        if (k < 1 .or. k > size(obs,1)) return
        if (.not. valid_kmeans_options(restarts,threshold,max_iterations)) return
        ntries = 20
        if (present(restarts)) ntries = restarts
        do t=1,ntries
            call initialize_centers(obs,k,state,'points',guess,code)
            call run_kmeans(obs,guess,candidate,candidate_dist,code,threshold,max_iterations)
            if (code /= cluster_vq_success) cycle
            if (candidate_dist<distortion) then
                distortion = candidate_dist
                centroids = candidate
                status = cluster_vq_success
            end if
        end do
    end subroutine kmeans_count

    subroutine kmeans_matrix(obs, initial, centroids, distortion, status, &
        threshold, max_iterations)
        real(dp), intent(in) :: obs(:,:) !! Finite observations, (n, features).
        real(dp), intent(in) :: initial(:,:) !! Initial codebook, (k, features).
        real(dp), allocatable, intent(out) :: centroids(:,:) !! Nonempty final centers.
        real(dp), intent(out) :: distortion !! Mean unsquared Euclidean distance, or huge on error.
        integer, intent(out) :: status !! Zero on success, one for invalid input.
        real(dp), intent(in), optional :: threshold !! Convergence tolerance, default 1e-5.
        integer, intent(in), optional :: max_iterations !! Max updates, default 300.

        allocate(centroids(0,0))
        status = cluster_vq_invalid_input
        distortion = huge(1.0_dp)
        if (.not. valid_centers(obs,initial)) return
        if (.not. valid_kmeans_options(threshold=threshold,max_iterations=max_iterations)) return
        deallocate(centroids)
        call run_kmeans(obs,initial,centroids,distortion,status,threshold,max_iterations)
    end subroutine kmeans_matrix

    pure logical function valid_kmeans_options(restarts,threshold,max_iterations) result(valid)
        integer, intent(in), optional :: restarts !! Optional positive restart count.
        real(dp), intent(in), optional :: threshold !! Optional finite nonnegative threshold.
        integer, intent(in), optional :: max_iterations !! Optional positive iteration cap.
        valid = .true.
        if (present(restarts)) valid = restarts>=1
        if (present(threshold)) then
            valid = valid .and. ieee_is_finite(threshold)
            valid = valid .and. threshold>=0.0_dp
        end if
        if (present(max_iterations)) valid = valid .and. max_iterations>=1
    end function valid_kmeans_options

    subroutine run_kmeans(obs, initial, centroids, distortion, status, threshold, max_iterations)
        real(dp), intent(in) :: obs(:,:) !! Nonempty finite observation matrix.
        real(dp), intent(in) :: initial(:,:) !! Nonempty finite initial codebook.
        real(dp), allocatable, intent(out) :: centroids(:,:) !! Final nonempty centers.
        real(dp), intent(out) :: distortion !! Mean non-squared Euclidean distance.
        integer, intent(out) :: status !! Zero on success, one for numerical failure.
        real(dp), intent(in), optional :: threshold !! Mean distance convergence threshold.
        integer, intent(in), optional :: max_iterations !! Iteration cap, default 300.
        real(dp), allocatable :: updated(:,:), distances(:)
        integer, allocatable :: labels(:), counts(:)
        real(dp) :: tolerance, previous, current
        integer :: nsteps, iter, i, j, k, nonempty, code

        nsteps = 300
        if (present(max_iterations)) nsteps = max_iterations
        tolerance = 1.0e-5_dp
        if (present(threshold)) tolerance = threshold
        centroids = initial
        previous = huge(1.0_dp)
        status = cluster_vq_success
        distortion = huge(1.0_dp)
        do iter=1,nsteps
            call vq_matrix(obs,centroids,labels,distances,code)
            k = size(centroids,1)
            allocate(counts(k),updated(k,size(obs,2)))
            counts = 0
            do i = 1,size(labels)
                counts(labels(i)) = counts(labels(i))+1
            end do
            updated = 0.0_dp
            do i=1,size(labels)
                j = labels(i)
                updated(j,:) = updated(j,:) + obs(i,:)/real(counts(j),dp)
            end do
            nonempty = count(counts>0)
            centroids = pack_centers(updated,counts,nonempty)
            deallocate(counts,updated)
            call vq_matrix(obs,centroids,labels,distances,code)
            current = sum(distances/real(size(distances),dp))
            distortion = current
            if (.not. ieee_is_finite(current)) then
                status = cluster_vq_invalid_input
                deallocate(centroids)
                allocate(centroids(0,0))
                return
            end if
            if (abs(current-previous)<=tolerance) exit
            previous = current
        end do
    end subroutine run_kmeans

    pure function pack_centers(centers,counts,nkeep) result(compacted)
        real(dp), intent(in) :: centers(:,:) !! Center estimates, some possibly unused.
        integer, intent(in) :: counts(:) !! Population of each estimated center.
        integer, intent(in) :: nkeep !! Number of nonempty centers.
        real(dp), allocatable :: compacted(:,:)
        integer :: i, j

        allocate(compacted(nkeep,size(centers,2)))
        j = 0
        do i = 1,size(counts)
            if (counts(i)==0) cycle
            j = j+1
            compacted(j,:) = centers(i,:)
        end do
    end function pack_centers

end module scifort_cluster_vq
