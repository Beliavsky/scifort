! SPDX-License-Identifier: MIT
program test_spatial_rotation_estimation
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, &
        rotation_from_quat, rotation_from_rotvec, rotation_align_vectors, rotation_mean, &
        spatial_rotation_success, spatial_rotation_invalid_input
    implicit none
    type(rotation_type) :: r, fitted, samples(2), mean_rot, invalid
    real(dp) :: a(3,3), b(3,3), weights(3), rssd, sensitivity(3,3), target(3,3)
    real(dp) :: inf, pi, q(4)
    integer :: i, status

    pi = acos(-1.0_dp)
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    r = rotation_from_rotvec([0.0_dp, 0.0_dp, 0.5_dp*pi])
    b = 0.0_dp
    do i=1,3
        b(i,i) = 1.0_dp
        a(i,:) = r%apply(b(i,:))
    end do
    fitted = rotation_align_vectors(a,b,rssd=rssd,sensitivity=sensitivity,status=status)
    if (status /= spatial_rotation_success .or. .not. fitted%valid()) error stop 'Wahba alignment'
    call same_rotation(fitted,r,1.0e-13_dp)
    if (abs(rssd) > 1.0e-12_dp) error stop 'Wahba residual'
    target = 0.0_dp
    do i=1,3
        target(i,i) = 0.5_dp
    end do
    call close_matrix(sensitivity,target,2.0e-14_dp)

    weights = [1.0_dp,2.0_dp,4.0_dp]
    fitted = rotation_align_vectors(a,b,weights=weights,sensitivity=sensitivity,status=status)
    if (status /= spatial_rotation_success) error stop 'weighted Wahba alignment'
    call same_rotation(fitted,r,1.0e-13_dp)
    target = 0.0_dp
    target(1,1) = 7.0_dp/15.0_dp
    target(2,2) = 7.0_dp/18.0_dp
    target(3,3) = 7.0_dp/9.0_dp
    call close_matrix(sensitivity,target,2.0e-14_dp)

    ! A single pair has a shortest alignment, but no meaningful sensitivity.
    fitted = rotation_align_vectors(a(1:1,:), b(1:1,:),rssd=rssd,status=status)
    if (status /= 0) error stop 'single-vector alignment'
    call close_vector(fitted%apply(b(1,:)),a(1,:),1.0e-14_dp)
    fitted = rotation_align_vectors(a(1:1,:), b(1:1,:),sensitivity=sensitivity,status=status)
    if (status /= spatial_rotation_invalid_input .or. fitted%valid()) &
        error stop 'single-vector sensitivity accepted'

    ! Near-antipodal data preserve the tiny transverse component.
    b(1,:) = [1.0_dp,0.0_dp,0.0_dp]
    a(1,:) = [-1.0_dp,1.0e-8_dp,0.0_dp]
    fitted = rotation_align_vectors(a(1:1,:),b(1:1,:),status=status)
    if (status /= 0) error stop 'near-antipodal alignment'
    call close_vector(fitted%apply(b(1,:)),a(1,:)/sqrt(sum(a(1,:)**2)),5.0e-13_dp)

    ! A unique positive infinite weight exactly fixes the primary direction;
    ! secondary pairs determine the remaining twist about that direction.
    a(1,:) = [0.0_dp,2.0_dp,0.0_dp]
    b(1,:) = [1.0_dp,0.0_dp,0.0_dp]
    a(2,:) = [-1.0_dp,0.0_dp,0.0_dp]
    b(2,:) = [0.0_dp,1.0_dp,0.0_dp]
    fitted = rotation_align_vectors(a(1:2,:),b(1:2,:), &
                                    weights=[inf,1.0_dp],rssd=rssd,status=status)
    if (status /= 0) error stop 'infinite primary weight'
    call same_rotation(fitted,r,2.0e-14_dp)
    if (abs(rssd) > 1.0e-14_dp) error stop 'infinite weight residual'
    fitted = rotation_align_vectors(a(1:2,:),b(1:2,:), &
                                    weights=[inf,inf],status=status)
    if (status /= 1 .or. fitted%valid()) error stop 'multiple infinite weights accepted'
    fitted = rotation_align_vectors(a(1:2,:),b(1:2,:), &
                                    weights=[inf,1.0_dp],sensitivity=sensitivity,status=status)
    if (status /= 1 .or. fitted%valid()) error stop 'infinite weight sensitivity accepted'

    samples(1) = rotation_identity()
    samples(2) = r
    mean_rot = rotation_mean(samples,status=status)
    if (status /= 0) error stop 'unweighted chordal mean'
    fitted = rotation_from_rotvec([0.0_dp,0.0_dp,0.25_dp*pi])
    call same_rotation(mean_rot,fitted,2.0e-14_dp)
    mean_rot = rotation_mean(samples,weights=[3.0_dp,1.0_dp],status=status)
    if (status /= 0) error stop 'weighted chordal mean'
    fitted = rotation_from_rotvec([0.0_dp,0.0_dp,atan2(1.0_dp,3.0_dp)])
    call same_rotation(mean_rot,fitted,2.0e-14_dp)

    q = r%as_quat()
    samples(2) = rotation_from_quat(-q)
    mean_rot = rotation_mean(samples,weights=[1.0e-300_dp,1.0e300_dp],status=status)
    if (status /= 0) error stop 'extreme chordal mean weights'
    call same_rotation(mean_rot,r,2.0e-14_dp)
    invalid = rotation_mean(samples,weights=[0.0_dp,0.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'zero weight mean accepted'
    invalid = rotation_mean(samples,weights=[1.0_dp,-1.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'negative mean weight accepted'
    invalid = rotation_align_vectors(a(:,1:2),b,status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'wrong vector shape accepted'
    invalid = rotation_align_vectors(a,b,weights=[1.0_dp,0.0_dp,1.0_dp],status=status)
    if (status /= 1 .or. invalid%valid()) error stop 'nonpositive alignment weight accepted'

    print '(a)', 'test_spatial_rotation_estimation PASS'
contains
    subroutine close_vector(actual,expected,tol)
        real(dp), intent(in) :: actual(:) !! computed vector
        real(dp), intent(in) :: expected(:) !! expected vector
        real(dp), intent(in) :: tol !! absolute maximum discrepancy
        if (size(actual) /= size(expected)) error stop 'vector size mismatch'
        if (.not. all(abs(actual-expected) <= tol)) then
            print *, actual,expected
            error stop 'vector mismatch'
        end if
    end subroutine close_vector
    subroutine close_matrix(actual,expected,tol)
        real(dp), intent(in) :: actual(:,:) !! computed matrix
        real(dp), intent(in) :: expected(:,:) !! expected matrix
        real(dp), intent(in) :: tol !! absolute maximum discrepancy
        if (any(shape(actual) /= shape(expected))) error stop 'matrix size mismatch'
        if (.not. all(abs(actual-expected) <= tol)) then
            print *, actual
            print *, expected
            error stop 'matrix mismatch'
        end if
    end subroutine close_matrix
    subroutine same_rotation(first,second,tol)
        type(rotation_type), intent(in) :: first !! first orientation
        type(rotation_type), intent(in) :: second !! second orientation
        real(dp), intent(in) :: tol !! maximal quaternion component discrepancy modulo sign
        real(dp) :: x(4),y(4)
        if (.not. first%valid() .or. .not. second%valid()) error stop 'invalid rotation'
        x = first%as_quat()
        y = second%as_quat()
        if (min(maxval(abs(x-y)), maxval(abs(x+y))) > tol) then
            print *, x,y
            error stop 'rotation mismatch'
        end if
    end subroutine same_rotation
end program test_spatial_rotation_estimation
