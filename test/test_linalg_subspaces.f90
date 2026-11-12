! SPDX-License-Identifier: MIT
program test_linalg_subspaces
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_subspaces, only : orth, null_space, subspace_angles
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(2,5), a2(3,3), r1(2,1), r2(2,1), projector(5,5)
    real(dp) :: empty(0,4), zero_cols(4,0), bad(3,1)
    real(dp), allocatable :: q(:, :), z(:, :), angles(:)
    complex(dp) :: ac(2,5), c1(2,1)
    complex(dp), allocatable :: qc(:, :), zc(:, :)
    real(dp), parameter :: tol = 2.0e-12_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: status, j

    a = 0.0_dp
    a(1,1) = 2.0_dp
    a(2,2) = 3.0_dp
    a(:,3) = a(:,1)
    call orth(a,q,status)
    if (status /= linalg_status_success) error stop 'real orth status'
    if (any(shape(q) /= [2,2])) error stop 'real orth shape'
    if (maxval(abs(matmul(q,transpose(q))-eye(2))) > tol) error stop 'real range projector'
    call null_space(a,z,status)
    if (status /= linalg_status_success) error stop 'real null status'
    if (any(shape(z) /= [5,3])) error stop 'wide null shape/completion'
    if (maxval(abs(matmul(a,z))) > tol) error stop 'wide null residual'
    if (maxval(abs(matmul(transpose(z),z)-eye(3))) > tol) error stop 'wide null orthogonality'
    projector = eye(5)
    projector(1,1) = 0.5_dp
    projector(3,3) = 0.5_dp
    projector(1,3) = -0.5_dp
    projector(3,1) = -0.5_dp
    projector(2,2) = 0.0_dp
    if (maxval(abs(matmul(z,transpose(z))-projector)) > tol) &
        error stop 'wide null projector'

    a2 = 0.0_dp
    a2(1,1) = 1.0_dp
    a2(2,1) = 1.0_dp
    a2(1,2) = 1.0_dp
    a2(2,2) = 1.0_dp
    a2(3,3) = 5.0_dp
    call null_space(a2,z,status)
    if (status /= linalg_status_success .or. size(z,2) /= 1) error stop 'tall deficient null rank'
    if (maxval(abs(matmul(a2,z))) > tol) error stop 'tall deficient null residual'
    if (abs(sum(z(:,1)**2)-1.0_dp) > tol) error stop 'tall deficient orthogonality'
    call null_space(0.0_dp*a2,z,status,rcond=2.0_dp)
    if (status /= linalg_status_success .or. size(z,2) /= 3) error stop 'zero null dimension'
    if (maxval(abs(z-eye(3))) > tol) error stop 'zero null identity'
    call orth(0.0_dp*a2,q,status)
    if (status /= linalg_status_success .or. size(q,2) /= 0) error stop 'zero orth dimension'
    ! Test rank cutoffs far enough from roundoff to give unambiguous ranks.
    a2 = 0.0_dp
    a2(1,1) = 1.0_dp
    a2(2,2) = 1.0e-9_dp
    call orth(a2,q,status,rcond=1.0e-8_dp)
    if (status /= linalg_status_success .or. size(q,2) /= 1) &
        error stop 'truncated orth rank'
    call null_space(a2,z,status,rcond=1.0e-8_dp)
    if (status /= linalg_status_success .or. size(z,2) /= 2) &
        error stop 'truncated null rank'
    call orth(a2,q,status,rcond=1.0e-10_dp)
    if (status /= linalg_status_success .or. size(q,2) /= 2) &
        error stop 'retained orth rank'

    ac = cmplx(a,0.0_dp,dp)
    ac(1,:) = cmplx(0.0_dp,1.0_dp,dp)*ac(1,:)
    call orth(ac,qc,status)
    if (status /= linalg_status_success .or. size(qc,2) /= 2) error stop 'complex orth status'
    call null_space(ac,zc,status)
    if (status /= linalg_status_success .or. size(zc,2) /= 3) error stop 'complex null status'
    if (maxval(abs(matmul(ac,zc))) > tol) error stop 'complex null residual'
    if (maxval(abs(matmul(transpose(conjg(zc)),zc)-eye(3))) > tol) &
        error stop 'complex null orthogonality'

    r1(:,1) = [1.0_dp,0.0_dp]
    r2(:,1) = [cos(1.0e-9_dp),sin(1.0e-9_dp)]
    call subspace_angles(r1,r2,angles,status)
    if (status /= linalg_status_success .or. size(angles) /= 1) error stop 'small angle shape'
    if (abs(angles(1)-1.0e-9_dp) > 2.0e-13_dp) error stop 'near-zero angle accuracy'
    r2(:,1) = [sin(1.0e-9_dp),cos(1.0e-9_dp)]
    call subspace_angles(r1,r2,angles,status)
    if (status /= linalg_status_success) error stop 'near right angle status'
    if (abs(angles(1)-(pi/2.0_dp-1.0e-9_dp)) > 2.0e-13_dp) &
        error stop 'near right angle accuracy'
    c1(:,1) = cmplx(r2(:,1),0.0_dp,dp)
    call subspace_angles(r1,c1,angles,status)
    if (status /= linalg_status_success) error stop 'mixed real-complex status'
    if (abs(angles(1)-(pi/2.0_dp-1.0e-9_dp)) > 2.0e-13_dp) &
        error stop 'mixed real-complex value'
    call subspace_angles(c1,r1,angles,status)
    if (status /= linalg_status_success) error stop 'mixed complex-real status'
    if (abs(angles(1)-(pi/2.0_dp-1.0e-9_dp)) > 2.0e-13_dp) &
        error stop 'mixed complex-real value'

    call orth(empty,q,status)
    if (status /= linalg_status_success .or. any(shape(q) /= [0,0])) &
        error stop 'empty orth'
    call null_space(empty,z,status)
    if (status /= linalg_status_success .or. any(shape(z) /= [4,4])) &
        error stop 'empty rows null'
    if (maxval(abs(z-eye(4))) > tol) error stop 'empty rows null identity'
    call null_space(zero_cols,z,status)
    if (status /= linalg_status_success .or. any(shape(z) /= [0,0])) &
        error stop 'empty columns null'
    call subspace_angles(zero_cols,zero_cols,angles,status)
    if (status /= linalg_status_success .or. size(angles) /= 0) error stop 'empty angles'
    bad = 0.0_dp
    call subspace_angles(r1,bad,angles,status)
    if (status /= linalg_status_invalid_shape) error stop 'mismatched rows accepted'
    call orth(a,q,status,rcond=-1.0_dp)
    if (status /= linalg_status_invalid_shape) error stop 'negative rcond accepted'
    call null_space(ac,zc,status,max_sweeps=0)
    if (status /= linalg_status_invalid_shape) error stop 'invalid max_sweeps accepted'
    call subspace_angles(r1,r2,angles,status,rcond=-1.0_dp)
    if (status /= linalg_status_invalid_shape) error stop 'invalid angle rcond accepted'
    a(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call null_space(a,z,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite null input accepted'
    call orth(a,q,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite orth input accepted'
    do j = 1, 2
        if (j == 1) then
            call subspace_angles(a,r2,angles,status)
        else
            call subspace_angles(r2,a,angles,status)
        end if
        if (status /= linalg_status_no_convergence) error stop 'nonfinite angles input accepted'
    end do

    print '(A)', 'linalg subspaces native tests passed'
contains
    function eye(n) result(e)
        integer, intent(in) :: n !! order of the returned identity
        real(dp) :: e(n,n)
        integer :: k
        e = 0.0_dp
        do k = 1, n
            e(k,k) = 1.0_dp
        end do
    end function eye
end program test_linalg_subspaces
