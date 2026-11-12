! SPDX-License-Identifier: MIT
program test_linalg_procrustes
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_procrustes, only : orthogonal_procrustes
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(4,3), b(4,3), r(3,3), reference(3,3), factor
    real(dp) :: a_big(2,2), b_small(2,2), r2(2,2)
    complex(dp) :: ac(4,3), bc(4,3), rc(3,3), q(3,3)
    real(dp), allocatable :: empty_rows(:, :), empty_columns(:, :)
    real(dp), allocatable :: empty_r(:, :)
    integer :: i, status

    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 3.0_dp
    reference = 0.0_dp
    reference(1,2) = 1.0_dp
    reference(2,3) = 1.0_dp
    reference(3,1) = 1.0_dp
    b = matmul(a, reference)
    call orthogonal_procrustes(a, b, r, factor, status)
    if (status /= linalg_status_success) error stop 'real procrustes failed'
    if (maxval(abs(r-reference)) > 1.0e-12_dp) error stop 'real mapping'
    if (abs(factor-14.0_dp) > 1.0e-12_dp) error stop 'real scale'
    if (maxval(abs(matmul(transpose(r),r)-identity_real(3))) > 1.0e-12_dp) &
        error stop 'real orthogonality'

    ac = cmplx(a, 0.0_dp, dp)
    q = cmplx(0.0_dp, 0.0_dp, dp)
    q(1,2) = cmplx(0.0_dp, 1.0_dp, dp)
    q(2,3) = cmplx(0.0_dp,-1.0_dp, dp)
    q(3,1) = cmplx(1.0_dp, 0.0_dp, dp)
    bc = matmul(ac, q)
    call orthogonal_procrustes(ac, bc, rc, factor, status)
    if (status /= linalg_status_success) error stop 'complex procrustes failed'
    if (maxval(abs(rc-q)) > 1.0e-12_dp) error stop 'complex mapping'
    if (abs(factor-14.0_dp) > 1.0e-12_dp) error stop 'complex scale'
    call orthogonal_procrustes(a, bc, rc, factor, status)
    if (status /= linalg_status_success) error stop 'real-complex mapping failed'
    if (maxval(abs(rc-q)) > 1.0e-12_dp) error stop 'real-complex mapping'
    call orthogonal_procrustes(ac, b, rc, factor, status)
    if (status /= linalg_status_success) error stop 'complex-real mapping failed'
    if (maxval(abs(rc-cmplx(reference,0.0_dp,dp))) > 1.0e-12_dp) &
        error stop 'complex-real mapping'
    if (maxval(abs(matmul(transpose(conjg(rc)),rc)-identity_complex(3))) > 1.0e-12_dp) &
        error stop 'complex unitarity'

    ! Rank-deficient matrices also have valid (non-unique) minimizers.
    a(:,2:3) = a(:,1:2)
    b = a
    call orthogonal_procrustes(a,b,r,factor,status)
    if (status /= linalg_status_success) error stop 'rank-deficient status'
    if (maxval(abs(matmul(a,r)-b)) > 1.0e-11_dp) error stop 'rank-deficient residual'
    if (maxval(abs(matmul(transpose(r),r)-identity_real(3))) > 1.0e-11_dp) &
        error stop 'rank-deficient orthogonality'

    b = 0.0_dp
    call orthogonal_procrustes(a,b,r,factor,status)
    if (status /= linalg_status_success .or. factor /= 0.0_dp) &
        error stop 'zero target'
    if (maxval(abs(r-identity_real(3))) > 1.0e-14_dp) error stop 'zero target identity'

    allocate(empty_rows(0,3))
    call orthogonal_procrustes(empty_rows,empty_rows,r,factor,status)
    if (status /= linalg_status_success .or. factor /= 0.0_dp) error stop 'empty rows'
    if (maxval(abs(r-identity_real(3))) /= 0.0_dp) error stop 'empty rows identity'
    allocate(empty_columns(4,0),empty_r(0,0))
    call orthogonal_procrustes(empty_columns,empty_columns,empty_r,factor,status)
    if (status /= linalg_status_success .or. factor /= 0.0_dp) error stop 'empty columns'

    a = 0.0_dp
    b = 0.0_dp
    call orthogonal_procrustes(a,b,r(1:2,:),factor,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad R shape'
    call orthogonal_procrustes(a,b(1:3,:),r,factor,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad B shape'
    call orthogonal_procrustes(a,b,r,factor,status,max_sweeps=0)
    if (status /= linalg_status_invalid_shape) error stop 'bad sweeps'
    a(1,1) = ieee_value(1.0_dp,ieee_quiet_nan)
    call orthogonal_procrustes(a,b,r,factor,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite data'
    bc(1,1) = cmplx(ieee_value(1.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call orthogonal_procrustes(ac,bc,rc,factor,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite complex data'
    bc(1,1) = cmplx(0.0_dp,0.0_dp,dp)
    call orthogonal_procrustes(a,bc,rc,factor,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite mixed data'

    ! The raw cross-product over/underflows in intermediate arithmetic;
    ! its true nonzero singular values remain representable.
    a_big = 0.0_dp
    b_small = 0.0_dp
    a_big(1,1) = 1.0e200_dp
    a_big(2,2) = 3.0e200_dp
    b_small(1,1) = 2.0e-200_dp
    b_small(2,2) = 4.0e-200_dp
    call orthogonal_procrustes(a_big,b_small,r2,factor,status)
    if (status /= linalg_status_success) error stop 'opposite scales status'
    if (abs(factor-14.0_dp) > 1.0e-12_dp) error stop 'opposite scales value'

    a_big = 0.0_dp
    b_small = 0.0_dp
    a_big(1,1) = 1.0e200_dp
    b_small(1,1) = 1.0e200_dp
    call orthogonal_procrustes(a_big,b_small,r2,factor,status)
    if (status /= linalg_status_no_convergence) error stop 'overflowed true scale'
    if (any(r2 /= 0.0_dp)) error stop 'overflow status leaves nonzero factor'
    print *, 'linalg orthogonal procrustes tests passed'

contains

    pure function identity_real(n) result(eye)
        integer, intent(in) :: n !! matrix order
        real(dp) :: eye(n,n)
        integer :: k
        eye = 0.0_dp
        do k = 1, n
            eye(k,k) = 1.0_dp
        end do
    end function identity_real

    pure function identity_complex(n) result(eye)
        integer, intent(in) :: n !! matrix order
        complex(dp) :: eye(n,n)
        integer :: k
        eye = cmplx(0.0_dp,0.0_dp,dp)
        do k = 1, n
            eye(k,k) = cmplx(1.0_dp,0.0_dp,dp)
        end do
    end function identity_complex

end program test_linalg_procrustes
