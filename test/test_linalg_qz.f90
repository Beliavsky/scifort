! SPDX-License-Identifier: MIT
! Native checks for inverse-based generalized complex Schur factors.
program test_linalg_qz
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_qz, only : qz
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: ar(3,3), br(3,3), invalid(3,2), singular(3,3)
    complex(dp) :: ac(3,3), bc(3,3), aa(3,3), bb(3,3), q(3,3), z(3,3)
    complex(dp) :: alpha(3), beta(3), wrong(2,2)
    complex(dp) :: zero_matrix(0,0)
    integer :: status, i

    ar = 0.0_dp
    br = 0.0_dp
    do i = 1, 3
        br(i,i) = 1.0_dp + 0.2_dp*real(i,dp)
    end do
    ar(1,1) = 1.0_dp
    ar(2,2) = 3.0_dp
    ar(3,3) = -2.0_dp
    ar(1,2) = 0.4_dp
    ar(2,1) = 1.2_dp
    ar(2,3) = -0.3_dp
    ar(3,1) = 0.2_dp
    br(1,2) = 0.5_dp
    br(2,3) = -0.3_dp
    call qz(ar,br,aa,bb,q,z,status,alpha,beta)
    if (status /= linalg_status_success) error stop 'real pair'
    call verify(cmplx(ar,kind=dp),cmplx(br,kind=dp),aa,bb,q,z)
    do i = 1, 3
        if (abs(alpha(i)-aa(i,i)) > 1e-12_dp) error stop 'alpha'
        if (abs(beta(i)-bb(i,i)) > 1e-12_dp) error stop 'beta'
    end do

    ac = cmplx(ar,kind=dp)
    bc = cmplx(br,kind=dp)
    ac(1,3) = cmplx(0.3_dp, -0.5_dp,dp)
    ac(2,1) = cmplx(0.4_dp, 0.3_dp,dp)
    bc(1,2) = cmplx(0.5_dp, -0.1_dp,dp)
    call qz(ac,bc,aa,bb,q,z,status)
    if (status /= linalg_status_success) error stop 'complex pair'
    call verify(ac,bc,aa,bb,q,z)
    call qz(ar,bc,aa,bb,q,z,status)
    if (status /= linalg_status_success) error stop 'mixed real-complex'
    call verify(cmplx(ar,kind=dp),bc,aa,bb,q,z)
    call qz(ac,br,aa,bb,q,z,status)
    if (status /= linalg_status_success) error stop 'mixed complex-real'
    call verify(ac,cmplx(br,kind=dp),aa,bb,q,z)

    ar = 0.0_dp
    call qz(ar,br,aa,bb,q,z,status)
    if (status /= linalg_status_success) error stop 'zero A'
    call verify(cmplx(ar,kind=dp),cmplx(br,kind=dp),aa,bb,q,z)
    ar(1,1) = 1e-200_dp
    br = br*1e200_dp
    call qz(ar,br,aa,bb,q,z,status)
    if (status /= linalg_status_success) error stop 'extreme input scale'
    call verify(cmplx(ar,kind=dp),cmplx(br,kind=dp),aa,bb,q,z)

    ! A nonsingular / B singular (including B=0) admits a valid complex
    ! generalized Schur decomposition, with zero beta for infinite roots.
    singular = 0.0_dp
    ar = 0.0_dp
    do i=1,3
        ar(i,i) = 1.0_dp + 0.5_dp*real(i,dp)
    end do
    call qz(ar,singular,aa,bb,q,z,status,alpha,beta)
    if (status /= linalg_status_success) error stop 'singular B, invertible A'
    call verify(cmplx(ar,kind=dp),cmplx(singular,kind=dp),aa,bb,q,z)
    if (maxval(abs(beta)) > 0.0_dp) error stop 'zero B beta'
    singular(1,1) = 1.0_dp
    call qz(ar,singular,aa,bb,q,z,status)
    if (status /= linalg_status_success) error stop 'rank-deficient B'
    call verify(cmplx(ar,kind=dp),cmplx(singular,kind=dp),aa,bb,q,z)
    ! A and B can both be singular yet define a regular pencil.
    ar = 0.0_dp
    singular = 0.0_dp
    ar(1,1) = 1.0_dp
    ar(2,2) = 2.0_dp
    singular(2,2) = 1.0_dp
    singular(3,3) = 1.0_dp
    call qz(ar,singular,aa,bb,q,z,status)
    if (status /= linalg_status_success) error stop 'regular singular pair'
    call verify(cmplx(ar,kind=dp),cmplx(singular,kind=dp),aa,bb,q,z)
    ar = 0.0_dp
    call qz(ar,singular,aa,bb,q,z,status)
    if (status /= linalg_status_singular) error stop 'singular pencil accepted'
    invalid = 0.0_dp
    call qz(invalid,br,aa,bb,q,z,status)
    if (status /= linalg_status_invalid_shape) error stop 'malformed A'
    call qz(ar,br,wrong,bb,q,z,status)
    if (status /= linalg_status_invalid_shape) error stop 'malformed AA'
    call qz(ar,br,aa,bb,q,z,status,alpha=alpha(1:2))
    if (status /= linalg_status_invalid_shape) error stop 'malformed alpha'
    ar(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call qz(ar,br,aa,bb,q,z,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite A'
    call qz(zero_matrix,zero_matrix,zero_matrix,zero_matrix,zero_matrix,zero_matrix,status)
    if (status /= linalg_status_success) error stop 'empty pair'

    print '(a)', 'generalized Schur native tests passed'

contains
    subroutine verify(a,b,aa,bb,q,z)
        complex(dp), intent(in) :: a(:,:) !! original coefficient matrix
        complex(dp), intent(in) :: b(:,:) !! original nonsingular coefficient matrix
        complex(dp), intent(in) :: aa(:,:) !! upper triangular first factor
        complex(dp), intent(in) :: bb(:,:) !! upper triangular second factor
        complex(dp), intent(in) :: q(:,:) !! first unitary basis
        complex(dp), intent(in) :: z(:,:) !! second unitary basis
        complex(dp) :: identity(size(a,1),size(a,1))
        integer :: k, n
        real(dp) :: r_a, r_b, norm_a, norm_b

        n = size(a,1)
        identity = cmplx(0.0_dp,0.0_dp,dp)
        do k=1,n
            identity(k,k) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        norm_a = max(1.0_dp,maxval(abs(a)))
        norm_b = max(1.0_dp,maxval(abs(b)))
        r_a = maxval(abs(matmul(q,matmul(aa,transpose(conjg(z))))-a))/norm_a
        r_b = maxval(abs(matmul(q,matmul(bb,transpose(conjg(z))))-b))/norm_b
        if (r_a > 1e-10_dp .or. r_b > 1e-10_dp) error stop 'QZ reconstruction'
        if (maxval(abs(matmul(transpose(conjg(q)),q)-identity)) > 1e-10_dp) &
            error stop 'QZ left unitarity'
        if (maxval(abs(matmul(transpose(conjg(z)),z)-identity)) > 1e-10_dp) &
            error stop 'QZ right unitarity'
        do k=1,n
            if (abs(aimag(bb(k,k))) > 1e-10_dp*norm_b) error stop 'beta phase'
            if (real(bb(k,k),dp) < 0.0_dp) error stop 'negative beta'
            if (k<n) then
                if (maxval(abs(aa(k+1:n,k))) > 1e-10_dp*norm_a) &
                    error stop 'not triangular A'
                if (maxval(abs(bb(k+1:n,k))) > 1e-10_dp*norm_b) &
                    error stop 'not triangular B'
            end if
        end do
    end subroutine verify
end program test_linalg_qz
