! SPDX-License-Identifier: MIT
program test_spatial_procrustes
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial, only : procrustes, orthogonal_procrustes, &
        spatial_distance_success, spatial_distance_invalid_input
    implicit none
    real(dp) :: a(6,3), b(6,3), c(6,3), r(3,3), u(4,5), v(4,5)
    real(dp) :: s, disparity, original_norm, error_max
    real(dp), allocatable :: mt1(:,:), mt2(:,:), rotation(:,:)
    integer :: code, i, j

    ! General nondegenerate 3-D transformation, including a reflection.
    a(1,:) = [1.0_dp,3.0_dp,5.0_dp]
    a(2,:) = [-2.0_dp,7.0_dp,0.0_dp]
    a(3,:) = [6.0_dp,4.0_dp,-3.0_dp]
    a(4,:) = [0.0_dp,0.0_dp,2.0_dp]
    a(5,:) = [1.0_dp,-5.0_dp,8.0_dp]
    a(6,:) = [4.0_dp,-2.0_dp,1.0_dp]
    r(1,:) = [0.0_dp,0.0_dp,-1.0_dp]
    r(2,:) = [-1.0_dp,0.0_dp,0.0_dp]
    r(3,:) = [0.0_dp,1.0_dp,0.0_dp]
    b = 12.0_dp*matmul(a,r)
    do j=1,3
        b(:,j) = b(:,j)+20.0_dp*real(j,dp)
    end do
    call procrustes(a,b,mt1,mt2,disparity,status=code)
    if (code /= spatial_distance_success) error stop 'procrustes rigid fit status'
    if (any(shape(mt1)/=[6,3]) .or. any(shape(mt2)/=[6,3])) &
        error stop 'procrustes shape'
    if (disparity > 1.0e-25_dp) error stop 'procrustes transformed disparity'
    call close(mt1,mt2,2.0e-12_dp)
    if (abs(sum(mt1**2)-1.0_dp)>1.0e-14_dp) error stop 'reference unit trace'
    if (any(abs(sum(mt1,dim=1))>1.0e-14_dp)) error stop 'reference centering'
    call procrustes(b,a,mt1,mt2,s,status=code)
    if (code /= 0 .or. abs(s-disparity)>1.0e-12_dp) &
        error stop 'symmetric disparity'

    ! Orthogonal-only fit must preserve norm and recover a reflection.
    call orthogonal_procrustes(a,matmul(a,r),rotation,s,status=code)
    if (code /= 0) error stop 'orthogonal procrustes status'
    if (any(shape(rotation)/=[3,3])) error stop 'orthogonal shape'
    call close(rotation,r,1.0e-12_dp)
    call check_identity(matmul(transpose(rotation),rotation),2.0e-12_dp)
    original_norm = sum(a**2)
    if (abs(s-original_norm)>2.0e-12_dp*original_norm) &
        error stop 'orthogonal singular value sum'

    ! Rank-deficient fits: only one occupied geometric direction, and
    ! more features than samples. Singular SVD columns must be completed.
    u = 0.0_dp
    do i=1,4
        u(i,:) = real(i,dp)*[1.0_dp,2.0_dp,0.0_dp,0.0_dp,0.0_dp]
    end do
    v(:,1) = -3.0_dp*u(:,2)+8.0_dp
    v(:,2) = 3.0_dp*u(:,1)-2.0_dp
    v(:,3:5) = 7.0_dp
    call procrustes(u,v,mt1,mt2,disparity,status=code)
    if (code /= 0 .or. disparity > 2.0e-24_dp) &
        error stop 'rank-deficient procrustes'
    call close(mt1,mt2,2.0e-12_dp)
    call orthogonal_procrustes(u, v, rotation, s, status=code)
    if (code /= 0) error stop 'rank-deficient orthogonal procrustes'
    call check_identity(matmul(transpose(rotation),rotation),2.0e-12_dp)

    ! Two non-collinear, non-similar point sets have positive disparity.
    c = a
    c(2,1) = c(2,1)+0.75_dp
    call procrustes(a,c,mt1,mt2,disparity,status=code)
    if (code /= 0 .or. disparity <= 0.0_dp .or. disparity >= 1.0_dp) &
        error stop 'positive disparity'
    error_max = maxval(abs(mt1-mt2))
    if (abs(disparity-sum((mt1-mt2)**2))>1.0e-15_dp .or. error_max<=0.0_dp) &
        error stop 'disparity definition'

    ! A one-dimensional sample is valid and the fit can reflect its sign.
    call procrustes(reshape([1.0_dp,2.0_dp,5.0_dp],[3,1]), &
        reshape([-1.0_dp,-2.0_dp,-5.0_dp],[3,1]),mt1,mt2,disparity,status=code)
    if (code /= 0 .or. disparity>1.0e-25_dp) error stop '1D reflection'

    ! Arbitrary binary64 scaling must not overflow intermediate variances.
    call procrustes(a*1.0e200_dp,b*1.0e-200_dp,mt1,mt2,disparity,status=code)
    if (code /= 0 .or. disparity>1.0e-23_dp) &
        error stop 'extreme-scale procrustes'
    call procrustes(a*1.0e-200_dp,b*1.0e200_dp,mt1,mt2,disparity,status=code)
    if (code /= 0 .or. disparity>1.0e-23_dp) &
        error stop 'inverted-scale procrustes'

    ! The orthogonal-only API permits zero matrices, returning identity.
    call orthogonal_procrustes(0.0_dp*a,b,rotation,s,status=code)
    if (code /= 0 .or. s /= 0.0_dp) error stop 'zero orthogonal input'
    call check_identity(rotation,0.0_dp)

    ! Shape mismatch, all-equal point clouds, empty dimensions, and NaNs.
    call procrustes(a,b(:,1:2),mt1,mt2,disparity,status=code)
    call invalid_result(mt1,mt2,disparity,code)
    call procrustes(a,0.0_dp*b,mt1,mt2,disparity,status=code)
    call invalid_result(mt1,mt2,disparity,code)
    call procrustes(a(:1,:),b(:1,:),mt1,mt2,disparity,status=code)
    call invalid_result(mt1,mt2,disparity,code)
    call procrustes(a(:,1:0),b(:,1:0),mt1,mt2,disparity,status=code)
    call invalid_result(mt1,mt2,disparity,code)
    c = a
    c(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call procrustes(a,c,mt1,mt2,disparity,status=code)
    call invalid_result(mt1,mt2,disparity,code)
    call orthogonal_procrustes(a,b(:,1:2),rotation,s,status=code)
    if (code /= spatial_distance_invalid_input .or. any(shape(rotation)/=[0,0]) &
        .or. .not. ieee_is_nan(s)) error stop 'orthogonal invalid shape'
    call orthogonal_procrustes(a,c,rotation,s,status=code)
    if (code /= spatial_distance_invalid_input .or. any(shape(rotation)/=[0,0])) &
        error stop 'orthogonal nonfinite input'

    print '(a)', 'test_spatial_procrustes PASS'
contains
    subroutine close(x, y, tolerance)
        real(dp), intent(in) :: x(:,:) !! actual matrix
        real(dp), intent(in) :: y(:,:) !! expected matrix
        real(dp), intent(in) :: tolerance !! absolute error tolerance
        if (any(shape(x)/=shape(y))) error stop 'shape mismatch'
        if (.not. all(abs(x-y)<=tolerance)) then
            print *, 'largest difference:',maxval(abs(x-y))
            error stop 'matrix discrepancy'
        end if
    end subroutine close

    subroutine check_identity(x, tolerance)
        real(dp), intent(in) :: x(:,:) !! approximately square identity matrix
        real(dp), intent(in) :: tolerance !! maximum entrywise error
        integer :: k
        real(dp) :: identity(size(x,1),size(x,2))
        if (size(x,1)/=size(x,2)) error stop 'identity not square'
        identity = 0.0_dp
        do k=1,size(x,1)
            identity(k,k)=1.0_dp
        end do
        call close(x,identity,tolerance)
    end subroutine check_identity

    subroutine invalid_result(m1,m2,dis,status)
        real(dp), intent(in) :: m1(:,:) !! empty first result
        real(dp), intent(in) :: m2(:,:) !! empty second result
        real(dp), intent(in) :: dis !! NaN disparity
        integer, intent(in) :: status !! invalid-input status
        if (status /= spatial_distance_invalid_input .or. &
            any(shape(m1)/=[0,0]) .or. any(shape(m2)/=[0,0]) &
            .or. .not. ieee_is_nan(dis)) error stop 'invalid procrustes result'
    end subroutine invalid_result
end program test_spatial_procrustes
