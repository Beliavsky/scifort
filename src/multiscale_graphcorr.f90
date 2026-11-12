! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Adapted from SciPy 1.17.0 scipy/stats/_mgc.py and scipy/stats/_stats.pyx.
! Retained BSD terms: THIRD_PARTY_LICENSES.md; details: CODE_PROVENANCE.md.
!
! Multiscale graph correlation (MGC) independence and two-sample tests.

module scifort_multiscale_graphcorr
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_beta, only : beta_ppf
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_random, only : rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: mgc_status_success = 0
    integer, parameter, public :: mgc_status_invalid_shape = 1
    integer, parameter, public :: mgc_status_invalid_parameter = 2
    integer, parameter, public :: mgc_status_nonfinite = 3

    abstract interface
        subroutine mgc_distance_callback(data, distance, status)
            import :: dp
            real(dp), intent(in) :: data(:, :) !! observations by rows
            real(dp), allocatable, intent(out) :: distance(:, :) !! square distance matrix
            integer, intent(out) :: status !! zero on success
        end subroutine mgc_distance_callback
    end interface

    type, public :: mgc_result
        real(dp) :: statistic !! multiscale graph correlation statistic
        real(dp) :: pvalue !! permutation p-value
        real(dp) :: stat !! alias of statistic, matching SciPy's legacy attribute
        real(dp), allocatable :: mgc_map(:, :) !! local correlation map
        integer :: opt_scale(2) !! one-based optimal MGC scale
        real(dp), allocatable :: null_dist(:) !! permutation null distribution
        integer :: status !! zero on success; nonzero for invalid input
    end type mgc_result

    interface multiscale_graphcorr
        module procedure multiscale_graphcorr_vv
        module procedure multiscale_graphcorr_vm
        module procedure multiscale_graphcorr_mv
        module procedure multiscale_graphcorr_mm
    end interface multiscale_graphcorr

    public :: mgc_distance_callback
    public :: multiscale_graphcorr

contains

    function multiscale_graphcorr_vv(x, y, state, reps, is_twosamp, compute_distance) result(res)
        real(dp), intent(in) :: x(:) !! first one-dimensional sample
        real(dp), intent(in) :: y(:) !! second one-dimensional sample
        type(rng_state), intent(inout) :: state !! explicit RNG state used by permutations
        integer, intent(in), optional :: reps !! number of permutation replications (default 1000)
        logical, intent(in), optional :: is_twosamp !! force unpaired two-sample test
        procedure(mgc_distance_callback), optional :: compute_distance !! custom pairwise distance callback
        type(mgc_result) :: res

        real(dp), allocatable :: xm(:, :)
        real(dp), allocatable :: ym(:, :)

        allocate(xm(size(x), 1), ym(size(y), 1))
        xm(:, 1) = x
        ym(:, 1) = y
        res = multiscale_graphcorr_mm(xm, ym, state, reps=reps, is_twosamp=is_twosamp, &
            compute_distance=compute_distance)
    end function multiscale_graphcorr_vv

    function multiscale_graphcorr_vm(x, y, state, reps, is_twosamp, compute_distance) result(res)
        real(dp), intent(in) :: x(:) !! first one-dimensional sample
        real(dp), intent(in) :: y(:, :) !! second sample, observations by rows
        type(rng_state), intent(inout) :: state !! explicit RNG state used by permutations
        integer, intent(in), optional :: reps !! number of permutation replications (default 1000)
        logical, intent(in), optional :: is_twosamp !! force unpaired two-sample test
        procedure(mgc_distance_callback), optional :: compute_distance !! custom pairwise distance callback
        type(mgc_result) :: res

        real(dp), allocatable :: xm(:, :)

        allocate(xm(size(x), 1))
        xm(:, 1) = x
        res = multiscale_graphcorr_mm(xm, y, state, reps=reps, is_twosamp=is_twosamp, &
            compute_distance=compute_distance)
    end function multiscale_graphcorr_vm

    function multiscale_graphcorr_mv(x, y, state, reps, is_twosamp, compute_distance) result(res)
        real(dp), intent(in) :: x(:, :) !! first sample, observations by rows
        real(dp), intent(in) :: y(:) !! second one-dimensional sample
        type(rng_state), intent(inout) :: state !! explicit RNG state used by permutations
        integer, intent(in), optional :: reps !! number of permutation replications (default 1000)
        logical, intent(in), optional :: is_twosamp !! force unpaired two-sample test
        procedure(mgc_distance_callback), optional :: compute_distance !! custom pairwise distance callback
        type(mgc_result) :: res

        real(dp), allocatable :: ym(:, :)

        allocate(ym(size(y), 1))
        ym(:, 1) = y
        res = multiscale_graphcorr_mm(x, ym, state, reps=reps, is_twosamp=is_twosamp, &
            compute_distance=compute_distance)
    end function multiscale_graphcorr_mv

    function multiscale_graphcorr_mm(x, y, state, reps, is_twosamp, precomputed, compute_distance) result(res)
        real(dp), intent(in) :: x(:, :) !! first sample or square distance matrix
        real(dp), intent(in) :: y(:, :) !! second sample or square distance matrix
        type(rng_state), intent(inout) :: state !! explicit RNG state used by permutations
        integer, intent(in), optional :: reps !! number of permutation replications (default 1000)
        logical, intent(in), optional :: is_twosamp !! force unpaired two-sample test
        logical, intent(in), optional :: precomputed !! true when x and y are distance matrices
        procedure(mgc_distance_callback), optional :: compute_distance !! custom pairwise distance callback
        type(mgc_result) :: res

        integer :: distance_status
        integer :: nreps
        logical :: distance_input
        logical :: two_sample
        real(dp), allocatable :: distx(:, :)
        real(dp), allocatable :: disty(:, :)
        real(dp), allocatable :: tx(:, :)
        real(dp), allocatable :: ty(:, :)

        call invalidate_mgc_result(res)
        nreps = 1000
        if (present(reps)) nreps = reps
        if (nreps < 0) then
            res%status = mgc_status_invalid_parameter
            return
        end if
        if (size(x, 1) < 5 .or. size(y, 1) < 5) then
            res%status = mgc_status_invalid_shape
            return
        end if
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) then
            res%status = mgc_status_nonfinite
            return
        end if

        distance_input = .false.
        if (present(precomputed)) distance_input = precomputed
        two_sample = .false.
        if (present(is_twosamp)) two_sample = is_twosamp

        if (distance_input) then
            if (present(compute_distance)) then
                res%status = mgc_status_invalid_parameter
                return
            end if
            if (two_sample) then
                res%status = mgc_status_invalid_parameter
                return
            end if
            if (size(x, 1) /= size(x, 2) .or. size(y, 1) /= size(y, 2) .or. &
                    size(x, 1) /= size(y, 1)) then
                res%status = mgc_status_invalid_shape
                return
            end if
            distx = x
            disty = y
        else
            if (size(x, 1) /= size(y, 1)) then
                if (size(x, 2) /= size(y, 2)) then
                    res%status = mgc_status_invalid_shape
                    return
                end if
                two_sample = .true.
            end if

            if (two_sample) then
                if (size(x, 2) /= size(y, 2)) then
                    res%status = mgc_status_invalid_shape
                    return
                end if
                call two_sample_transform(x, y, tx, ty)
                if (present(compute_distance)) then
                    call compute_distance(tx, distx, distance_status)
                    if (distance_status /= 0 .or. .not. valid_distance_matrix(distx, size(tx, 1))) then
                        res%status = mgc_status_invalid_parameter
                        return
                    end if
                    call compute_distance(ty, disty, distance_status)
                    if (distance_status /= 0 .or. .not. valid_distance_matrix(disty, size(ty, 1))) then
                        res%status = mgc_status_invalid_parameter
                        return
                    end if
                else
                    call euclidean_distance_matrix(tx, distx)
                    call euclidean_distance_matrix(ty, disty)
                end if
            else
                if (size(x, 1) /= size(y, 1)) then
                    res%status = mgc_status_invalid_shape
                    return
                end if
                if (present(compute_distance)) then
                    call compute_distance(x, distx, distance_status)
                    if (distance_status /= 0 .or. .not. valid_distance_matrix(distx, size(x, 1))) then
                        res%status = mgc_status_invalid_parameter
                        return
                    end if
                    call compute_distance(y, disty, distance_status)
                    if (distance_status /= 0 .or. .not. valid_distance_matrix(disty, size(y, 1))) then
                        res%status = mgc_status_invalid_parameter
                        return
                    end if
                else
                    call euclidean_distance_matrix(x, distx)
                    call euclidean_distance_matrix(y, disty)
                end if
            end if
        end if

        call mgc_statistic(distx, disty, res%statistic, res%mgc_map, res%opt_scale)
        if (.not. ieee_is_finite(res%statistic)) then
            res%status = mgc_status_invalid_parameter
            return
        end if
        res%stat = res%statistic
        call mgc_permutation_test(distx, disty, res%statistic, nreps, state, &
            res%pvalue, res%null_dist)
        res%status = mgc_status_success
    end function multiscale_graphcorr_mm

    subroutine invalidate_mgc_result(res)
        type(mgc_result), intent(out) :: res !! result to initialize to invalid values

        res%statistic = quiet_nan(0.0_dp)
        res%pvalue = quiet_nan(0.0_dp)
        res%stat = quiet_nan(0.0_dp)
        res%opt_scale = 0
        res%status = mgc_status_invalid_parameter
        if (allocated(res%mgc_map)) deallocate(res%mgc_map)
        if (allocated(res%null_dist)) deallocate(res%null_dist)
    end subroutine invalidate_mgc_result

    logical function valid_distance_matrix(distance, n) result(valid)
        real(dp), intent(in), allocatable :: distance(:, :) !! callback-produced distance matrix
        integer, intent(in) :: n !! required matrix order

        valid = allocated(distance)
        if (.not. valid) return
        valid = size(distance, 1) == n .and. size(distance, 2) == n
        if (.not. valid) return
        valid = all(ieee_is_finite(distance))
    end function valid_distance_matrix

    subroutine euclidean_distance_matrix(data, distance)
        real(dp), intent(in) :: data(:, :) !! observations by rows
        real(dp), allocatable, intent(out) :: distance(:, :) !! pairwise Euclidean distances

        integer :: i
        integer :: j
        real(dp) :: norm_factor
        real(dp) :: scale_value
        real(dp), allocatable :: difference(:)

        allocate(distance(size(data, 1), size(data, 1)), difference(size(data, 2)))
        distance = 0.0_dp
        do j = 2, size(data, 1)
            do i = 1, j - 1
                difference = data(i, :) - data(j, :)
                scale_value = maxval(abs(difference))
                if (scale_value <= 0.0_dp) then
                    distance(i, j) = 0.0_dp
                else
                    norm_factor = sqrt(sum((difference / scale_value)**2))
                    if (scale_value > huge(1.0_dp) / norm_factor) then
                        distance(i, j) = huge(1.0_dp)
                    else
                        distance(i, j) = scale_value * norm_factor
                    end if
                end if
                distance(j, i) = distance(i, j)
            end do
        end do
    end subroutine euclidean_distance_matrix

    subroutine two_sample_transform(x, y, combined, labels)
        real(dp), intent(in) :: x(:, :) !! first sample
        real(dp), intent(in) :: y(:, :) !! second sample
        real(dp), allocatable, intent(out) :: combined(:, :) !! concatenated observations
        real(dp), allocatable, intent(out) :: labels(:, :) !! zero/one sample labels

        integer :: nx
        integer :: ny

        nx = size(x, 1)
        ny = size(y, 1)
        allocate(combined(nx + ny, size(x, 2)), labels(nx + ny, 1))
        combined(1:nx, :) = x
        combined(nx + 1:nx + ny, :) = y
        labels(1:nx, 1) = 0.0_dp
        labels(nx + 1:nx + ny, 1) = 1.0_dp
    end subroutine two_sample_transform

    subroutine mgc_permutation_test(distx, disty, observed, reps, state, pvalue, null_dist)
        real(dp), intent(in) :: distx(:, :) !! first distance matrix
        real(dp), intent(in) :: disty(:, :) !! second distance matrix
        real(dp), intent(in) :: observed !! observed MGC statistic
        integer, intent(in) :: reps !! number of permutations
        type(rng_state), intent(inout) :: state !! RNG state advanced by shuffles
        real(dp), intent(out) :: pvalue !! plus-one permutation p-value
        real(dp), allocatable, intent(out) :: null_dist(:) !! permutation statistics

        integer :: b
        integer :: i
        integer, allocatable :: order(:)
        real(dp) :: perm_stat
        real(dp), allocatable :: perm_y(:, :)
        real(dp), allocatable :: temp_map(:, :)
        integer :: temp_scale(2)

        allocate(null_dist(reps), order(size(disty, 1)), perm_y(size(disty, 1), size(disty, 2)))
        b = 0
        do i = 1, reps
            call random_permutation(state, order)
            call permute_distance_matrix(disty, order, perm_y)
            call mgc_statistic(distx, perm_y, perm_stat, temp_map, temp_scale)
            null_dist(i) = perm_stat
            if (perm_stat >= observed) b = b + 1
            if (allocated(temp_map)) deallocate(temp_map)
        end do
        pvalue = real(b + 1, dp) / real(reps + 1, dp)
    end subroutine mgc_permutation_test

    subroutine random_permutation(state, order)
        type(rng_state), intent(inout) :: state !! RNG state advanced by Fisher-Yates shuffling
        integer, intent(out) :: order(:) !! one-based random permutation

        integer :: i
        integer :: j
        integer :: temp

        do i = 1, size(order)
            order(i) = i
        end do
        do i = size(order), 2, -1
            j = 1 + int(rng_uniform(state) * real(i, dp))
            if (j > i) j = i
            temp = order(i)
            order(i) = order(j)
            order(j) = temp
        end do
    end subroutine random_permutation

    subroutine permute_distance_matrix(distance, order, permuted)
        real(dp), intent(in) :: distance(:, :) !! square distance matrix
        integer, intent(in) :: order(:) !! row/column permutation
        real(dp), intent(out) :: permuted(:, :) !! jointly permuted distance matrix

        integer :: i
        integer :: j

        do j = 1, size(order)
            do i = 1, size(order)
                permuted(i, j) = distance(order(i), order(j))
            end do
        end do
    end subroutine permute_distance_matrix

    subroutine mgc_statistic(distx, disty, statistic, corr_map, opt_scale)
        real(dp), intent(in) :: distx(:, :) !! first square distance matrix
        real(dp), intent(in) :: disty(:, :) !! second square distance matrix
        real(dp), intent(out) :: statistic !! smoothed MGC statistic
        real(dp), allocatable, intent(out) :: corr_map(:, :) !! local correlation map
        integer, intent(out) :: opt_scale(2) !! optimal one-based scale

        logical, allocatable :: significant(:, :)
        integer :: m
        integer :: n
        integer :: sample_size

        call local_correlations(distx, disty, corr_map)
        m = size(corr_map, 1)
        n = size(corr_map, 2)
        if (m == 1 .or. n == 1) then
            statistic = corr_map(m, n)
            opt_scale = [m, n]
            return
        end if

        sample_size = size(distx, 1) - 1
        call threshold_mgc_map(corr_map, sample_size, significant)
        call smooth_mgc_map(significant, corr_map, statistic, opt_scale)
    end subroutine mgc_statistic

    subroutine local_correlations(distx, disty, corr_map)
        real(dp), intent(in) :: distx(:, :) !! first square distance matrix
        real(dp), intent(in) :: disty(:, :) !! second square distance matrix
        real(dp), allocatable, intent(out) :: corr_map(:, :) !! local correlations

        real(dp), allocatable :: centx(:, :)
        real(dp), allocatable :: centy(:, :)
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: varx_map(:, :)
        real(dp), allocatable :: vary_map(:, :)
        real(dp), allocatable :: varx(:)
        real(dp), allocatable :: vary(:)
        integer, allocatable :: rankx(:, :)
        integer, allocatable :: ranky(:, :)
        integer :: i
        integer :: j
        real(dp) :: denom

        call center_distance_matrix(distx, centx, rankx)
        call center_distance_matrix(disty, centy, ranky)
        call local_covariance(centx, transpose(centy), rankx, transpose(ranky), cov)
        call local_covariance(centx, transpose(centx), rankx, transpose(rankx), varx_map)
        call local_covariance(centy, transpose(centy), ranky, transpose(ranky), vary_map)

        allocate(varx(min(size(varx_map, 1), size(varx_map, 2))))
        allocate(vary(min(size(vary_map, 1), size(vary_map, 2))))
        do i = 1, size(varx)
            varx(i) = varx_map(i, i)
        end do
        do i = 1, size(vary)
            vary(i) = vary_map(i, i)
        end do

        allocate(corr_map(size(cov, 1), size(cov, 2)))
        corr_map = 0.0_dp
        do j = 1, size(corr_map, 2)
            do i = 1, size(corr_map, 1)
                if (i <= size(varx) .and. j <= size(vary)) then
                    if (varx(i) > 0.0_dp .and. vary(j) > 0.0_dp) then
                        denom = sqrt(varx(i) * vary(j))
                        if (denom > 0.0_dp) corr_map(i, j) = cov(i, j) / denom
                        if (corr_map(i, j) > 1.0_dp) corr_map(i, j) = 1.0_dp
                    end if
                end if
            end do
        end do
    end subroutine local_correlations

    subroutine center_distance_matrix(distance, centered, ranks)
        real(dp), intent(in) :: distance(:, :) !! square distance matrix
        real(dp), allocatable, intent(out) :: centered(:, :) !! MGC-centered distances
        integer, allocatable, intent(out) :: ranks(:, :) !! dense column-wise ranks

        integer :: i
        integer :: j
        integer :: n
        real(dp) :: expected

        n = size(distance, 1)
        call dense_rank_columns(distance, ranks)
        allocate(centered(n, n))
        do j = 1, n
            expected = sum(distance(:, j)) / real(n - 1, dp)
            do i = 1, n
                centered(i, j) = distance(i, j) - expected
            end do
        end do
        do i = 1, n
            centered(i, i) = 0.0_dp
        end do
    end subroutine center_distance_matrix

    subroutine dense_rank_columns(values, ranks)
        real(dp), intent(in) :: values(:, :) !! matrix to rank independently by columns
        integer, allocatable, intent(out) :: ranks(:, :) !! one-based dense ranks

        integer :: i
        integer :: j
        integer :: k
        integer :: rank_value
        integer, allocatable :: order(:)

        allocate(ranks(size(values, 1), size(values, 2)), order(size(values, 1)))
        do j = 1, size(values, 2)
            do i = 1, size(order)
                order(i) = i
            end do
            call sort_indices_by_values(values(:, j), order)
            rank_value = 1
            ranks(order(1), j) = rank_value
            do k = 2, size(order)
                if (values(order(k), j) < values(order(k - 1), j) .or. &
                    values(order(k), j) > values(order(k - 1), j)) rank_value = rank_value + 1
                ranks(order(k), j) = rank_value
            end do
        end do
    end subroutine dense_rank_columns

    subroutine sort_indices_by_values(values, order)
        real(dp), intent(in) :: values(:) !! finite values
        integer, intent(inout) :: order(:) !! indices to sort by corresponding values

        integer :: i
        integer :: j
        integer :: key

        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (values(order(j)) <= values(key)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine sort_indices_by_values

    subroutine local_covariance(distx, disty, rankx, ranky, cov)
        real(dp), intent(in) :: distx(:, :) !! centered first distances
        real(dp), intent(in) :: disty(:, :) !! centered/transposed second distances
        integer, intent(in) :: rankx(:, :) !! first dense ranks
        integer, intent(in) :: ranky(:, :) !! second dense ranks
        real(dp), allocatable, intent(out) :: cov(:, :) !! cumulative local covariance map

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: n
        integer :: nx
        integer :: ny
        real(dp), allocatable :: expectx(:)
        real(dp), allocatable :: expecty(:)

        n = size(distx, 1)
        nx = maxval(rankx)
        ny = maxval(ranky)
        allocate(cov(nx, ny), expectx(nx), expecty(ny))
        cov = 0.0_dp
        expectx = 0.0_dp
        expecty = 0.0_dp

        do j = 1, n
            do i = 1, n
                k = rankx(i, j)
                l = ranky(i, j)
                cov(k, l) = cov(k, l) + distx(i, j) * disty(i, j)
                expectx(k) = expectx(k) + distx(i, j)
                expecty(l) = expecty(l) + disty(i, j)
            end do
        end do

        do k = 2, nx
            cov(k, 1) = cov(k, 1) + cov(k - 1, 1)
            expectx(k) = expectx(k) + expectx(k - 1)
        end do
        do l = 2, ny
            cov(1, l) = cov(1, l) + cov(1, l - 1)
            expecty(l) = expecty(l) + expecty(l - 1)
        end do
        do k = 1, nx - 1
            do l = 1, ny - 1
                cov(k + 1, l + 1) = cov(k + 1, l + 1) + cov(k + 1, l) + &
                    cov(k, l + 1) - cov(k, l)
            end do
        end do
        do l = 1, ny
            do k = 1, nx
                cov(k, l) = cov(k, l) - expectx(k) * expecty(l) / (real(n, dp) * real(n, dp))
            end do
        end do
    end subroutine local_covariance

    subroutine threshold_mgc_map(corr_map, sample_size, significant)
        real(dp), intent(in) :: corr_map(:, :) !! local correlation map
        integer, intent(in) :: sample_size !! original n minus one
        logical, allocatable, intent(out) :: significant(:, :) !! largest significant component

        real(dp) :: a
        real(dp) :: percentile
        real(dp) :: threshold
        logical, allocatable :: raw(:, :)

        percentile = 1.0_dp - 0.02_dp / real(sample_size, dp)
        a = real(sample_size, dp) * real(sample_size - 3, dp) / 4.0_dp - 0.5_dp
        threshold = 2.0_dp * beta_ppf(percentile, a, a) - 1.0_dp
        threshold = max(threshold, corr_map(size(corr_map, 1), size(corr_map, 2)))
        allocate(raw(size(corr_map, 1), size(corr_map, 2)))
        raw = corr_map > threshold
        call largest_connected_component(raw, significant)
    end subroutine threshold_mgc_map

    subroutine largest_connected_component(mask, component)
        logical, intent(in) :: mask(:, :) !! candidate significant entries
        logical, allocatable, intent(out) :: component(:, :) !! largest 4-connected component

        integer :: best_count
        integer :: count_component
        integer :: head
        integer :: i
        integer :: j
        integer :: k
        integer :: ni
        integer :: nj
        integer :: tail
        integer, allocatable :: qi(:)
        integer, allocatable :: qj(:)
        logical, allocatable :: candidate(:, :)
        logical, allocatable :: visited(:, :)
        integer, parameter :: di(4) = [-1, 1, 0, 0]
        integer, parameter :: dj(4) = [0, 0, -1, 1]

        allocate(component(size(mask, 1), size(mask, 2)))
        component = .false.
        if (.not. any(mask)) return

        allocate(visited(size(mask, 1), size(mask, 2)))
        allocate(candidate(size(mask, 1), size(mask, 2)))
        allocate(qi(size(mask)), qj(size(mask)))
        visited = .false.
        best_count = 0

        do i = 1, size(mask, 1)
            do j = 1, size(mask, 2)
                if (.not. mask(i, j) .or. visited(i, j)) cycle
                candidate = .false.
                head = 1
                tail = 1
                qi(1) = i
                qj(1) = j
                visited(i, j) = .true.
                candidate(i, j) = .true.
                count_component = 1
                do while (head <= tail)
                    do k = 1, 4
                        ni = qi(head) + di(k)
                        nj = qj(head) + dj(k)
                        if (ni < 1 .or. ni > size(mask, 1)) cycle
                        if (nj < 1 .or. nj > size(mask, 2)) cycle
                        if (.not. mask(ni, nj) .or. visited(ni, nj)) cycle
                        tail = tail + 1
                        qi(tail) = ni
                        qj(tail) = nj
                        visited(ni, nj) = .true.
                        candidate(ni, nj) = .true.
                        count_component = count_component + 1
                    end do
                    head = head + 1
                end do
                if (count_component > best_count) then
                    best_count = count_component
                    component = candidate
                end if
            end do
        end do
    end subroutine largest_connected_component

    subroutine smooth_mgc_map(significant, corr_map, statistic, opt_scale)
        logical, intent(in) :: significant(:, :) !! largest significant connected region
        real(dp), intent(in) :: corr_map(:, :) !! local correlation map
        real(dp), intent(out) :: statistic !! smoothed maximum correlation
        integer, intent(out) :: opt_scale(2) !! selected one-based scale

        integer :: area_limit
        integer :: i
        integer :: j
        integer :: best_linear
        integer :: linear_index
        integer :: m
        integer :: n
        real(dp) :: max_corr

        m = size(corr_map, 1)
        n = size(corr_map, 2)
        statistic = corr_map(m, n)
        opt_scale = [m, n]
        if (.not. any(significant)) return

        area_limit = ceiling(0.02_dp * real(max(m, n), dp)) * min(m, n)
        if (count(significant) < area_limit) return
        max_corr = maxval(corr_map, mask=significant)
        if (max_corr < statistic) return

        statistic = max_corr
        best_linear = -1
        do i = 1, m
            do j = 1, n
                if (.not. significant(i, j)) cycle
                if (corr_map(i, j) < max_corr) cycle
                linear_index = (i - 1) * n + (j - 1)
                if (linear_index > best_linear) then
                    best_linear = linear_index
                    opt_scale = [i, j]
                end if
            end do
        end do
    end subroutine smooth_mgc_map

end module scifort_multiscale_graphcorr
