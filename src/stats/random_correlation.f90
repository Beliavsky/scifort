! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Random correlation matrices with prescribed eigenvalues. The construction
! follows the Davies-Higham similarity-plus-Givens algorithm:
! P. I. Davies and N. J. Higham, BIT 40 (2000), 640-651.

module scifort_random_correlation
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_ortho_group, only : ortho_group_rvs, ortho_group_status_success
    use scifort_random, only : rng_state
    implicit none
    private

    integer, parameter, public :: random_correlation_status_success = 0
    integer, parameter, public :: random_correlation_status_invalid_shape = 1
    integer, parameter, public :: random_correlation_status_invalid_parameter = 2
    integer, parameter, public :: random_correlation_status_numerical_failure = 3

    public :: random_correlation_rvs

contains

    subroutine random_correlation_rvs(state, eigs, sample, status, tol, diag_tol)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the orthogonal draw
        real(dp), intent(in) :: eigs(:) !! prescribed eigenvalues; length > 1, nonnegative, sum = dimension
        real(dp), intent(out) :: sample(:, :) !! output correlation matrix with eigenvalues eigs
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure
        real(dp), intent(in), optional :: tol !! parameter tolerance (default 1e-13)
        real(dp), intent(in), optional :: diag_tol !! final unit-diagonal tolerance (default 1e-7)

        integer :: i
        integer :: info
        integer :: n
        real(dp), allocatable :: q(:, :)
        real(dp), allocatable :: scaled(:, :)
        real(dp) :: parameter_tol
        real(dp) :: diagonal_tol

        n = size(eigs)
        sample = 0.0_dp
        parameter_tol = 1.0e-13_dp
        diagonal_tol = 1.0e-7_dp
        if (present(tol)) parameter_tol = tol
        if (present(diag_tol)) diagonal_tol = diag_tol

        if (n <= 1 .or. size(sample,1) /= n .or. size(sample,2) /= n) then
            if (present(status)) status = random_correlation_status_invalid_shape
            return
        end if
        if (.not. ieee_is_finite(parameter_tol) .or. parameter_tol < 0.0_dp .or. &
                .not. ieee_is_finite(diagonal_tol) .or. diagonal_tol < 0.0_dp .or. &
                .not. all(ieee_is_finite(eigs)) .or. &
                abs(sum(eigs) - real(n,dp)) > parameter_tol .or. &
                any(eigs < -parameter_tol)) then
            if (present(status)) status = random_correlation_status_invalid_parameter
            return
        end if

        allocate(q(n,n), scaled(n,n))
        call ortho_group_rvs(state, q, info)
        if (info /= ortho_group_status_success) then
            if (present(status)) status = random_correlation_status_numerical_failure
            return
        end if
        scaled = q
        do i = 1, n
            scaled(:,i) = scaled(:,i) * eigs(i)
        end do
        sample = matmul(scaled, transpose(q))
        call rotate_to_correlation(sample, info)
        if (info /= random_correlation_status_success) then
            if (present(status)) status = info
            return
        end if
        if (maxval(abs([(sample(i,i) - 1.0_dp, i=1,n)])) > diagonal_tol) then
            if (present(status)) status = random_correlation_status_numerical_failure
            return
        end if
        do i = 1, n
            sample(i,i) = 1.0_dp
        end do
        if (present(status)) status = random_correlation_status_success
    end subroutine random_correlation_rvs

    subroutine rotate_to_correlation(m, status)
        real(dp), intent(inout) :: m(:, :) !! symmetric PSD matrix with trace equal to its dimension
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        integer :: i
        integer :: j
        integer :: n
        real(dp) :: c
        real(dp) :: s

        n = size(m,1)
        if (size(m,2) /= n .or. n <= 1) then
            status = random_correlation_status_invalid_shape
            return
        end if

        do i = 1, n - 1
            if (m(i,i) == 1.0_dp) cycle
            j = find_partner(m, i)
            if (j <= i .or. j > n) then
                status = random_correlation_status_numerical_failure
                return
            end if
            call givens_to_one(m(i,i), m(j,j), m(i,j), c, s)
            call similarity_givens(m, i, j, c, s)
            m(i,i) = 1.0_dp
        end do
        status = random_correlation_status_success
    end subroutine rotate_to_correlation

    pure integer function find_partner(m, i) result(jfound)
        real(dp), intent(in) :: m(:, :) !! symmetric working matrix
        integer, intent(in) :: i !! diagonal index requiring a partner

        integer :: j
        integer :: n

        n = size(m,1)
        jfound = 0
        if (m(i,i) > 1.0_dp) then
            do j = i + 1, n
                if (m(j,j) < 1.0_dp) then
                    jfound = j
                    return
                end if
            end do
        else
            do j = i + 1, n
                if (m(j,j) > 1.0_dp) then
                    jfound = j
                    return
                end if
            end do
        end if

        ! Roundoff can make every remaining diagonal lie on the same side of
        ! one by a few ulps. Select the largest opposite-direction correction.
        if (i < n) then
            jfound = i + 1
            do j = i + 2, n
                if (m(i,i) > 1.0_dp) then
                    if (m(j,j) < m(jfound,jfound)) jfound = j
                else
                    if (m(j,j) > m(jfound,jfound)) jfound = j
                end if
            end do
        end if
    end function find_partner

    pure subroutine givens_to_one(aii, ajj, aij, c, s)
        real(dp), intent(in) :: aii !! first diagonal element of the 2x2 principal block
        real(dp), intent(in) :: ajj !! second diagonal element of the 2x2 principal block
        real(dp), intent(in) :: aij !! off-diagonal element of the 2x2 principal block
        real(dp), intent(out) :: c !! cosine of the similarity rotation
        real(dp), intent(out) :: s !! sine of the similarity rotation

        real(dp) :: aiid
        real(dp) :: ajjd
        real(dp) :: dd
        real(dp) :: t

        aiid = aii - 1.0_dp
        ajjd = ajj - 1.0_dp
        if (ajjd == 0.0_dp) then
            c = 0.0_dp
            s = 1.0_dp
            return
        end if
        dd = sqrt(max(aij * aij - aiid * ajjd, 0.0_dp))
        t = (aij + sign(dd, aij)) / ajjd
        c = 1.0_dp / sqrt(1.0_dp + t * t)
        if (c == 0.0_dp) then
            s = 1.0_dp
        else
            s = c * t
        end if
    end subroutine givens_to_one

    subroutine similarity_givens(m, i, j, c, s)
        real(dp), intent(inout) :: m(:, :) !! square matrix transformed in place as G^T M G
        integer, intent(in) :: i !! first rotation index
        integer, intent(in) :: j !! second rotation index
        real(dp), intent(in) :: c !! rotation cosine
        real(dp), intent(in) :: s !! rotation sine

        integer :: k
        integer :: n
        real(dp), allocatable :: col_i(:)
        real(dp), allocatable :: col_j(:)
        real(dp), allocatable :: row_i(:)
        real(dp), allocatable :: row_j(:)

        n = size(m,1)
        allocate(col_i(n), col_j(n), row_i(n), row_j(n))
        col_i = m(:,i)
        col_j = m(:,j)
        m(:,i) = c * col_i - s * col_j
        m(:,j) = s * col_i + c * col_j

        row_i = m(i,:)
        row_j = m(j,:)
        m(i,:) = c * row_i - s * row_j
        m(j,:) = s * row_i + c * row_j

        ! Restore symmetry explicitly to keep later partner searches stable.
        do k = 1, n
            m(k,i) = 0.5_dp * (m(k,i) + m(i,k))
            m(i,k) = m(k,i)
            m(k,j) = 0.5_dp * (m(k,j) + m(j,k))
            m(j,k) = m(k,j)
        end do
    end subroutine similarity_givens

end module scifort_random_correlation
