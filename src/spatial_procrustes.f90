! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Similarity and orthogonal Procrustes fits for real point configurations.
! Independent one-sided Jacobi polar factorization of a cross-product matrix.
! No LAPACK/BLAS or third-party source is required.
module scifort_spatial_procrustes
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : spatial_distance_success, &
        spatial_distance_invalid_input
    implicit none
    private

    integer, parameter, public :: spatial_procrustes_no_convergence = 2
    public :: procrustes, orthogonal_procrustes

contains

    ! One-sided Jacobi orthogonalization of C's columns supplies the polar
    ! factor U V^T directly. Unlike an eigendecomposition of C^T C, it does
    ! not square C's singular-value condition number.
    pure subroutine polar_rotation(cross, rotation, nuclear_norm, status)
        real(dp), intent(in) :: cross(:,:) !! square real cross-product matrix
        real(dp), intent(out) :: rotation(:,:) !! orthogonal factor, same shape as cross
        real(dp), intent(out) :: nuclear_norm !! sum of singular values of cross
        integer, intent(out) :: status !! 0 success; 1 invalid; 2 no convergence
        integer, parameter :: max_sweeps = 120
        real(dp), parameter :: relative_tolerance = 64.0_dp*epsilon(1.0_dp)
        real(dp) :: columns(size(cross,1),size(cross,2))
        real(dp) :: right(size(cross,1),size(cross,2))
        real(dp) :: left(size(cross,1),size(cross,2))
        real(dp) :: sigma(size(cross,1)), work(size(cross,1)), swap_col(size(cross,1))
        real(dp) :: a, b, g, angle, cosine, sine, rank_tol, magnitude, temp, small_column
        integer :: n, i, j, p, q, k, sweep, pass, best
        logical :: changed, converged

        rotation = 0.0_dp
        nuclear_norm = 0.0_dp
        status = spatial_distance_invalid_input
        n = size(cross,1)
        if (n < 1 .or. size(cross,2) /= n) return
        if (size(rotation,1) /= n .or. size(rotation,2) /= n) return
        if (.not. all(ieee_is_finite(cross))) return

        columns = cross
        small_column = 32.0_dp*epsilon(1.0_dp)*norm2(reshape(cross,[n*n]))
        right = 0.0_dp
        do i=1,n
            right(i,i) = 1.0_dp
        end do

        converged = .false.
        do sweep=1,max_sweeps
            changed = .false.
            do p=1,n-1
                do q=p+1,n
                    a = dot_product(columns(:,p),columns(:,p))
                    b = dot_product(columns(:,q),columns(:,q))
                    ! Do not chase roundoff in numerically null columns.
                    if (sqrt(a) <= small_column .or. sqrt(b) <= small_column) cycle
                    g = dot_product(columns(:,p),columns(:,q))
                    if (abs(g) <= relative_tolerance*sqrt(a)*sqrt(b)) cycle
                    ! Keep the angle in [-pi/4,pi/4] without subtracting
                    ! pi/2 from a nearly-pi/2 angle, which loses tiny turns.
                    if (b >= a) then
                        angle = 0.5_dp*atan2(2.0_dp*g,b-a)
                    else
                        angle = -0.5_dp*atan2(2.0_dp*g,a-b)
                    end if
                    ! A sub-ulp turn cannot improve the computed columns.
                    if (abs(angle) <= 2.0_dp*epsilon(1.0_dp)) cycle
                    cosine = cos(angle)
                    sine = sin(angle)
                    work = cosine*columns(:,p) - sine*columns(:,q)
                    columns(:,q) = sine*columns(:,p) + cosine*columns(:,q)
                    columns(:,p) = work
                    work = cosine*right(:,p) - sine*right(:,q)
                    right(:,q) = sine*right(:,p) + cosine*right(:,q)
                    right(:,p) = work
                    changed = .true.
                end do
            end do
            if (.not. changed) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) then
            status = spatial_procrustes_no_convergence
            return
        end if

        do j=1,n
            sigma(j) = norm2(columns(:,j))
        end do
        ! Put nonzero singular directions first, so completion of a null
        ! space cannot obstruct a later nonzero singular direction.
        do j=1,n-1
            best = maxloc(sigma(j:n),dim=1)+j-1
            if (best == j) cycle
            temp = sigma(j)
            sigma(j) = sigma(best)
            sigma(best) = temp
            swap_col = columns(:,j)
            columns(:,j) = columns(:,best)
            columns(:,best) = swap_col
            swap_col = right(:,j)
            right(:,j) = right(:,best)
            right(:,best) = swap_col
        end do
        nuclear_norm = sum(sigma)
        rank_tol = 32.0_dp*epsilon(1.0_dp)*real(n,dp)*maxval(sigma)
        left = 0.0_dp
        do j=1,n
            if (sigma(j) > rank_tol) then
                work = columns(:,j)/sigma(j)
                do pass=1,2
                    do k=1,j-1
                        work = work-dot_product(left(:,k),work)*left(:,k)
                    end do
                end do
                magnitude = norm2(work)
                if (magnitude < rank_tol) then
                    status = spatial_procrustes_no_convergence
                    return
                end if
                left(:,j) = work/magnitude
            else
                ! Deterministic orthogonal completion for rank-deficient C.
                magnitude = -1.0_dp
                do i=1,n
                    work = 0.0_dp
                    work(i) = 1.0_dp
                    do pass=1,2
                        do k=1,j-1
                            work = work-dot_product(left(:,k),work)*left(:,k)
                        end do
                    end do
                    temp = norm2(work)
                    if (temp > magnitude) then
                        magnitude = temp
                        swap_col = work
                    end if
                end do
                if (magnitude <= 0.0_dp) then
                    status = spatial_procrustes_no_convergence
                    return
                end if
                left(:,j) = swap_col/magnitude
            end if
        end do
        rotation = matmul(left,transpose(right))
        if (.not. all(ieee_is_finite(rotation))) then
            status = spatial_procrustes_no_convergence
            return
        end if
        status = spatial_distance_success
    end subroutine polar_rotation

    pure subroutine orthogonal_procrustes(a, b, rotation, scale, status)
        real(dp), intent(in) :: a(:,:) !! input observations (rows, coordinates)
        real(dp), intent(in) :: b(:,:) !! target observations, same shape as a
        real(dp), allocatable, intent(out) :: rotation(:,:) !! orthogonal mapping a R to b
        real(dp), intent(out) :: scale !! nuclear norm of transpose(a) times b
        integer, intent(out), optional :: status !! 0 success; 1 invalid; 2 no convergence
        real(dp), allocatable :: left(:,:), right(:,:), c(:,:)
        real(dp) :: magnitude_a, magnitude_b, norm_a, norm_b, score
        integer :: n, d, code, j

        code = spatial_distance_invalid_input
        scale = ieee_value(0.0_dp,ieee_quiet_nan)
        n = size(a,1)
        d = size(a,2)
        if (n < 1 .or. d < 1 .or. any(shape(a) /= shape(b))) goto 900
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) goto 900
        allocate(left(n,d),right(n,d),c(d,d),rotation(d,d))
        rotation = 0.0_dp
        do j=1,d
            rotation(j,j) = 1.0_dp
        end do
        magnitude_a = maxval(abs(a))
        magnitude_b = maxval(abs(b))
        if (magnitude_a == 0.0_dp .or. magnitude_b == 0.0_dp) then
            scale = 0.0_dp
            code = spatial_distance_success
            goto 900
        end if
        left = a/magnitude_a
        right = b/magnitude_b
        norm_a = norm2(reshape(left,[n*d]))
        norm_b = norm2(reshape(right,[n*d]))
        left = left/norm_a
        right = right/norm_b
        c = matmul(transpose(left),right)
        call polar_rotation(c,rotation,score,code)
        if (code /= spatial_distance_success) goto 900
        ! Grouping small and large factors helps avoid intermediate overflow.
        scale = ((score*norm_a)*magnitude_a)*(norm_b*magnitude_b)
        ! If the true result exceeds binary64 range, +Inf is valid output.
900     if (code /= spatial_distance_success) then
            if (allocated(rotation)) deallocate(rotation)
            allocate(rotation(0,0))
        end if
        if (present(status)) status = code
    end subroutine orthogonal_procrustes

    pure subroutine standardize_points(data, standardized, valid)
        real(dp), intent(in) :: data(:,:) !! finite rows of point coordinates
        real(dp), intent(out) :: standardized(:,:) !! centered, unit-Frobenius data
        logical, intent(out) :: valid !! true if input has nonzero centered norm
        real(dp) :: maximum, magnitude, anchor
        integer :: j, n, d

        valid = .false.
        n = size(data,1)
        d = size(data,2)
        if (n < 2 .or. d < 1) return
        if (.not. all(ieee_is_finite(data))) return
        maximum = maxval(abs(data))
        if (maximum == 0.0_dp) return
        standardized = data/maximum
        do j=1,d
            ! Center differences from a fixed row rather than summing
            ! potentially large offsets and subtracting nearly equal means.
            anchor = standardized(1,j)
            standardized(:,j) = standardized(:,j)-anchor
            standardized(:,j) = standardized(:,j) &
                - sum(standardized(:,j))/real(n,dp)
        end do
        magnitude = norm2(reshape(standardized,[n*d]))
        if (magnitude == 0.0_dp .or. .not. ieee_is_finite(magnitude)) return
        standardized = standardized/magnitude
        valid = .true.
    end subroutine standardize_points

    pure subroutine procrustes(data1, data2, mtx1, mtx2, disparity, status)
        real(dp), intent(in) :: data1(:,:) !! reference point configuration, rows are points
        real(dp), intent(in) :: data2(:,:) !! configuration to align, same shape as data1
        real(dp), allocatable, intent(out) :: mtx1(:,:) !! centered, unit-Frobenius reference
        real(dp), allocatable, intent(out) :: mtx2(:,:) !! rotated/reflected/scaled data2
        real(dp), intent(out) :: disparity !! sum of squared pointwise differences
        integer, intent(out), optional :: status !! 0 success; 1 invalid; 2 no convergence
        real(dp), allocatable :: rotation(:,:)
        real(dp) :: score
        integer :: n, d, code
        logical :: ok

        code = spatial_distance_invalid_input
        disparity = ieee_value(0.0_dp,ieee_quiet_nan)
        n = size(data1,1)
        d = size(data1,2)
        if (n < 2 .or. d < 1 .or. any(shape(data1)/=shape(data2))) goto 900
        allocate(mtx1(n,d),mtx2(n,d))
        call standardize_points(data1,mtx1,ok)
        if (.not. ok) goto 900
        call standardize_points(data2,mtx2,ok)
        if (.not. ok) goto 900
        call orthogonal_procrustes(mtx2,mtx1,rotation,score,code)
        if (code /= spatial_distance_success) goto 900
        mtx2 = score*matmul(mtx2,rotation)
        disparity = sum((mtx1-mtx2)**2)
900     if (code /= spatial_distance_success) then
            if (allocated(mtx1)) deallocate(mtx1)
            if (allocated(mtx2)) deallocate(mtx2)
            allocate(mtx1(0,0),mtx2(0,0))
        end if
        if (present(status)) status = code
    end subroutine procrustes
end module scifort_spatial_procrustes
