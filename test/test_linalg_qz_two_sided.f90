! SPDX-License-Identifier: MIT
! Test the inverse-free two-sided QZ route and public method dispatch.
program test_linalg_qz_two_sided
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_qz, only : qz
    implicit none
    real(dp) :: ar(5,5), br(5,5), wrong(4,4)
    complex(dp) :: ac(5,5), bc(5,5), aa(5,5), bb(5,5), q(5,5), z(5,5)
    complex(dp) :: alpha(5), beta(5)
    complex(dp) :: eye(5,5)
    integer :: i, j, status

    ar = 0.0_dp
    br = 0.0_dp
    eye = cmplx(0.0_dp,0.0_dp,dp)
    do i=1,5
        eye(i,i) = cmplx(1.0_dp,0.0_dp,dp)
        ar(i,i) = real(i,dp) + 0.15_dp
        br(i,i) = 1.0_dp / real(i,dp)
    end do
    ! Includes a conjugate eigenvalue pair, so complex shifts are necessary.
    ar(1,1) = 0.0_dp
    ar(2,2) = 0.0_dp
    ar(1,2) = -2.0_dp
    ar(2,1) = 2.0_dp
    ar(3,2) = 0.7_dp
    ar(4,3) = 0.5_dp
    ar(4,5) = -0.2_dp
    br(1,3) = 0.1_dp
    br(2,4) = -0.05_dp

    call qz(ar,br,aa,bb,q,z,status,alpha,beta,method='iterative')
    if (status /= linalg_status_success) error stop 'real matrix two-sided QZ'
    call verify(cmplx(ar,kind=dp),cmplx(br,kind=dp),aa,bb,q,z)
    do i=1,5
        if (alpha(i) /= aa(i,i) .or. beta(i) /= bb(i,i)) error stop 'homogeneous diagonal mismatch'
    end do
    call qz(ar,br,aa,bb,q,z,status,method='auto')
    if (status /= linalg_status_success) error stop 'auto method'
    call verify(cmplx(ar,kind=dp),cmplx(br,kind=dp),aa,bb,q,z)
    call qz(ar,br,aa,bb,q,z,status,method='inverse')
    if (status /= linalg_status_success) error stop 'inverse compatibility'
    call verify(cmplx(ar,kind=dp),cmplx(br,kind=dp),aa,bb,q,z)

    ac = cmplx(ar,kind=dp)
    bc = cmplx(br,kind=dp)
    ac(1,5) = cmplx(0.25_dp,-0.3_dp,dp)
    ac(4,1) = cmplx(0.15_dp,0.3_dp,dp)
    bc(3,2) = cmplx(0.1_dp,0.2_dp,dp)
    call qz(ac,bc,aa,bb,q,z,status,method='iterative')
    if (status /= linalg_status_success) error stop 'complex matrix two-sided QZ'
    call verify(ac,bc,aa,bb,q,z)

    ! Extreme scale difference between the two input matrices. Neither
    ! their explicit inverse nor the ratio of their scales is formed.
    ac = ac * 1.0e150_dp
    bc = bc * 1.0e-150_dp
    call qz(ac,bc,aa,bb,q,z,status,method='iterative')
    if (status /= linalg_status_success) error stop 'extreme scales'
    call verify(ac,bc,aa,bb,q,z)

    ! A highly ill-conditioned, but nonsingular, second matrix. The
    ! homogeneous representation avoids overflowing alpha/beta ratios.
    ac = cmplx(0.0_dp,0.0_dp,dp)
    bc = cmplx(0.0_dp,0.0_dp,dp)
    do i=1,5
        ac(i,i) = cmplx(real(i,dp),0.0_dp,dp)
        bc(i,i) = cmplx(10.0_dp**(-12*(i-1)),0.0_dp,dp)
    end do
    call qz(ac,bc,aa,bb,q,z,status,method='iterative')
    if (status /= linalg_status_success) error stop 'ill-conditioned B'
    call verify(ac,bc,aa,bb,q,z)

    ! Exactly singular B can still define a regular pencil.
    ac = cmplx(0.0_dp,0.0_dp,dp)
    bc = cmplx(0.0_dp,0.0_dp,dp)
    do i=1,5
        ac(i,i) = cmplx(real(i,dp),0.0_dp,dp)
        if (i < 5) bc(i,i) = cmplx(1.0_dp,0.0_dp,dp)
    end do
    call qz(ac,bc,aa,bb,q,z,status,method='iterative')
    if (status /= linalg_status_success) error stop 'infinite generalized eigenvalue'
    call verify(ac,bc,aa,bb,q,z)

    call qz(ar,br,aa,bb,q,z,status,method='not_a_method')
    if (status /= linalg_status_invalid_shape) error stop 'invalid method'
    wrong = 0.0_dp
    call qz(wrong,br,aa,bb,q,z,status,method='iterative')
    if (status /= linalg_status_invalid_shape) error stop 'invalid shape'
    print '(a)', 'two-sided generalized QZ native tests passed'

contains
    subroutine verify(a,b,aa,bb,q,z)
        complex(dp), intent(in) :: a(:,:) !! first input matrix
        complex(dp), intent(in) :: b(:,:) !! second input matrix
        complex(dp), intent(in) :: aa(:,:) !! triangular first matrix
        complex(dp), intent(in) :: bb(:,:) !! triangular second matrix
        complex(dp), intent(in) :: q(:,:) !! first unitary basis
        complex(dp), intent(in) :: z(:,:) !! second unitary basis
        real(dp) :: scale_a, scale_b, error_a, error_b
        integer :: n
        n = size(a,1)
        scale_a = maxval(abs(a))
        scale_b = maxval(abs(b))
        if (scale_a <= 0.0_dp) scale_a = 1.0_dp
        if (scale_b <= 0.0_dp) scale_b = 1.0_dp
        error_a = maxval(abs(matmul(q,matmul(aa/scale_a,transpose(conjg(z))))-a/scale_a))
        error_b = maxval(abs(matmul(q,matmul(bb/scale_b,transpose(conjg(z))))-b/scale_b))
        if (max(error_a,error_b) > 2.0e-10_dp) error stop 'reconstruction'
        if (maxval(abs(matmul(conjg(transpose(q)),q)-eye)) > 1.0e-10_dp) &
            error stop 'left unitarity'
        if (maxval(abs(matmul(conjg(transpose(z)),z)-eye)) > 1.0e-10_dp) &
            error stop 'right unitarity'
        do j=1,n-1
            if (maxval(abs(aa(j+1:n,j))) > 1.0e-10_dp*scale_a) &
                error stop 'first factor not triangular'
            if (maxval(abs(bb(j+1:n,j))) > 1.0e-10_dp*scale_b) &
                error stop 'second factor not triangular'
        end do
    end subroutine verify
end program test_linalg_qz_two_sided
