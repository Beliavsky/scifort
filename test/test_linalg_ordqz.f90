! SPDX-License-Identifier: MIT
program test_linalg_ordqz
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_ordqz, only : ordqz, ordqz_factors
    implicit none
    integer :: status, sdim, i
    real(dp) :: ar(4,4), br(4,4)
    complex(dp) :: a(4,4), b(4,4), aa(4,4), bb(4,4), q(4,4), z(4,4)
    complex(dp) :: alpha(4), beta(4), abefore(4,4)
    complex(dp) :: empty_a(0,0), empty_b(0,0), empty_q(0,0), empty_z(0,0)
    complex(dp) :: empty_aa(0,0), empty_bb(0,0)

    ar = 0.0_dp
    br = 0.0_dp
    do i=1,4
        br(i,i) = 1.0_dp
    end do
    ar(1,1) = 2.0_dp
    ar(2,2) = -3.0_dp
    ar(3,3) = 4.0_dp
    ar(4,4) = -1.0_dp
    ar(1,2) = 0.6_dp
    ar(1,3) = -0.4_dp
    ar(2,3) = 0.5_dp
    ar(2,4) = 0.8_dp
    call ordqz(ar,br,aa,bb,q,z,status,alpha,beta,sort='lhp',sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 2) error stop 'lhp ordering'
    if (real(alpha(1)/beta(1),dp) >= 0.0_dp .or. &
        real(alpha(2)/beta(2),dp) >= 0.0_dp) error stop 'left eigenvalues'
    if (real(alpha(3)/beta(3),dp) <= 0.0_dp .or. &
        real(alpha(4)/beta(4),dp) <= 0.0_dp) error stop 'right eigenvalues'
    call check_factors(cmplx(ar,kind=dp),cmplx(br,kind=dp),aa,bb,q,z)

    b = cmplx(br,kind=dp)
    a = cmplx(ar,kind=dp)
    a(1,1) = cmplx(0.0_dp,2.0_dp,dp)
    a(3,3) = cmplx(0.0_dp,-4.0_dp,dp)
    b(3,3) = cmplx(0.0_dp,0.0_dp,dp) ! one infinite eigenvalue
    a(4,4) = cmplx(0.5_dp,0.0_dp,dp)
    call ordqz(a,b,aa,bb,q,z,status,alpha,beta,sort='ouc',sdim=sdim)
    if (status /= linalg_status_success) error stop 'infinite eigenvalue ordering'
    if (sdim < 1) error stop 'infinite eigenvalue selection'
    if (minval(abs(beta(1:sdim))) > 1.0e-12_dp) &
        error stop 'infinite root not selected'
    call check_factors(a,b,aa,bb,q,z)

    ! Reordering already computed QZ factors does not rerun QZ and can work
    ! with a matrix supplied by a different generalized Schur solver.
    aa = a
    bb = b
    q = cmplx(0.0_dp,0.0_dp,dp)
    z = q
    do i=1,4
        q(i,i) = cmplx(1.0_dp,0.0_dp,dp)
        z(i,i) = cmplx(1.0_dp,0.0_dp,dp)
    end do
    call ordqz_factors(aa,bb,q,z,status,sort='iuc',sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 1) error stop 'factor sorting'
    call check_factors(a,b,aa,bb,q,z)
    call ordqz_factors(aa,bb,q,z,status,sdim=sdim,selector=select_infinite)
    if (status /= linalg_status_success .or. sdim /= 1) error stop 'callback sorting'
    call check_factors(a,b,aa,bb,q,z)

    abefore = aa
    call ordqz_factors(aa,bb,q,z,status,sort='nonsense')
    if (status /= linalg_status_invalid_shape) error stop 'invalid sort status'
    if (any(aa /= abefore)) error stop 'failure changed factors'
    call ordqz_factors(empty_a,empty_b,empty_q,empty_z,status,sort='lhp',sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 0) error stop 'empty pencil'
    call ordqz(empty_a,empty_b,empty_aa,empty_bb,empty_q,empty_z,status,sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 0) error stop 'empty ordqz'
    call ordqz(a,b,aa(1:3,1:3),bb,q,z,status,sort='lhp')
    if (status /= linalg_status_invalid_shape) error stop 'bad output shape'
    print '(a)', 'ordqz native tests passed'

contains
    logical function select_infinite(x,y)
        complex(dp), intent(in) :: x !! homogeneous numerator
        complex(dp), intent(in) :: y !! homogeneous denominator
        select_infinite = abs(y) < 1.0e-13_dp .and. abs(x) > 0.0_dp
    end function select_infinite

    subroutine check_factors(a,b,aa,bb,q,z)
        complex(dp), intent(in) :: a(:,:) !! original A
        complex(dp), intent(in) :: b(:,:) !! original B
        complex(dp), intent(in) :: aa(:,:) !! triangular factor A
        complex(dp), intent(in) :: bb(:,:) !! triangular factor B
        complex(dp), intent(in) :: q(:,:) !! left basis
        complex(dp), intent(in) :: z(:,:) !! right basis
        complex(dp) :: ident(size(q,1),size(q,1))
        real(dp) :: e
        integer :: j
        ident = cmplx(0.0_dp,0.0_dp,dp)
        do j=1,size(q,1)
            ident(j,j) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        e = maxval(abs(a-matmul(q,matmul(aa,transpose(conjg(z))))))
        e = max(e,maxval(abs(b-matmul(q,matmul(bb,transpose(conjg(z)))))))
        e = max(e,maxval(abs(matmul(transpose(conjg(q)),q)-ident)))
        e = max(e,maxval(abs(matmul(transpose(conjg(z)),z)-ident)))
        if (e > 2.0e-9_dp) then
            print *, 'reconstruction error', e
            error stop 'Schur factors incorrect'
        end if
        do j=1,size(aa,1)-1
            if (maxval(abs(aa(j+1:,j))) > 1.0e-10_dp) error stop 'A not upper triangular'
            if (maxval(abs(bb(j+1:,j))) > 1.0e-10_dp) error stop 'B not upper triangular'
        end do
    end subroutine check_factors
end program test_linalg_ordqz
