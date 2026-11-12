! SPDX-License-Identifier: MIT
program test_linalg_svd
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_svd, only : svd, pinv
    use scifort_linalg_qr, only : lstsq
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(4,3), u(4,3), s(3), vh(3,3), ap(3,4), x(3,2), b(4,2)
    real(dp) :: rerr, identity(3,3), v(3)
    complex(dp) :: ac(3,5), uc(3,3), vhc(3,5), apc(5,3), bc(3,2), xc(5,2)
    real(dp) :: sc(3)
    integer :: status, rank, i

    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(1,2) = 1.0_dp
    a(3,3) = 2.0_dp
    b = 0.0_dp
    b(:,1) = [2.0_dp,5.0_dp,4.0_dp,0.0_dp]
    b(:,2) = [6.0_dp,0.0_dp,8.0_dp,1.0_dp]
    identity = 0.0_dp
    do i = 1, 3
        identity(i,i) = 1.0_dp
    end do

    call svd(a,u,s,vh,status)
    if (status /= linalg_status_success) error stop 'real SVD status'
    rerr = maxval(abs(matmul(u*spread(s,1,4),vh)-a))
    if (rerr > 5.0e-13_dp) error stop 'real SVD reconstruction'
    if (maxval(abs(matmul(transpose(u),u)-identity)) > 5.0e-13_dp) &
        error stop 'real U orthogonality'
    if (maxval(abs(matmul(vh,transpose(vh))-identity)) > 5.0e-13_dp) &
        error stop 'real V orthogonality'
    if (abs(s(1)-2.0_dp) > 5.0e-13_dp .or. &
        abs(s(2)-sqrt(2.0_dp)) > 5.0e-13_dp .or. s(3) > 1.0e-13_dp) &
        error stop 'real singular values'

    call pinv(a,ap,status,rank)
    if (status /= linalg_status_success .or. rank /= 2) error stop 'real pinv rank'
    if (maxval(abs(matmul(a,matmul(ap,a))-a)) > 5.0e-13_dp) &
        error stop 'real Moore-Penrose identity'
    call lstsq(a,b,x,status,rank)
    if (status /= linalg_status_success .or. rank /= 2) error stop 'auto rank fallback'
    if (maxval(abs(x(:,1)-[1.0_dp,1.0_dp,2.0_dp])) > 5.0e-13_dp) &
        error stop 'real minimum norm'
    if (maxval(abs(x(:,2)-[3.0_dp,3.0_dp,4.0_dp])) > 5.0e-13_dp) &
        error stop 'real multiple rhs'
    call lstsq(a,b(:,1),v,status,rank,method='svd')
    if (status /= linalg_status_success .or. rank /= 2) error stop 'svd vector status'
    if (maxval(abs(v-x(:,1))) > 5.0e-13_dp) error stop 'real vector result'
    call lstsq(a,b(:,1),v,status,rank,method='qr')
    if (status /= linalg_status_singular) error stop 'explicit QR rank failure'
    call lstsq(a,b(:,1),v,status,method='other')
    if (status /= linalg_status_invalid_shape) error stop 'invalid method accepted'
    call pinv(a,ap,status,rank,rcond=0.9_dp)
    if (status /= linalg_status_success .or. rank /= 1) error stop 'real custom rcond'

    ! Force a genuinely complex rank-deficient wide matrix.
    ac = cmplx(0.0_dp,0.0_dp,dp)
    ac(1,1) = cmplx(1.0_dp,1.0_dp,dp)
    ac(1,2) = ac(1,1)
    ac(2,3) = cmplx(0.0_dp,2.0_dp,dp)
    ac(3,5) = cmplx(1.0_dp,-2.0_dp,dp)
    bc = cmplx(0.0_dp,0.0_dp,dp)
    bc(:,1) = [cmplx(2.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.0_dp,dp), &
               cmplx(-1.0_dp,3.0_dp,dp)]
    bc(:,2) = conjg(bc(:,1))
    call svd(ac,uc,sc,vhc,status)
    if (status /= linalg_status_success) error stop 'complex SVD status'
    if (maxval(abs(matmul(uc*spread(sc,1,3),vhc)-ac)) > 5.0e-13_dp) &
        error stop 'complex reconstruction'
    call pinv(ac,apc,status,rank)
    if (status /= linalg_status_success .or. rank /= 3) error stop 'complex pinv rank'
    if (maxval(abs(matmul(ac,matmul(apc,ac))-ac)) > 5.0e-13_dp) &
        error stop 'complex pinv identity'
    call lstsq(ac,bc,xc,status,rank,method='svd')
    if (status /= linalg_status_success .or. rank /= 3) error stop 'complex wide lstsq'
    if (maxval(abs(matmul(ac,xc)-bc)) > 5.0e-13_dp) error stop 'complex solution'
    if (maxval(abs(xc(:,1)-matmul(apc,bc(:,1)))) > 5.0e-13_dp) &
        error stop 'complex least norm'

    ! Analytical zero case, including a null-space basis.
    a = 0.0_dp
    call svd(a,u,s,vh,status)
    if (status /= linalg_status_success .or. any(s /= 0.0_dp)) error stop 'zero SVD'
    call pinv(a,ap,status,rank)
    if (status /= linalg_status_success .or. rank /= 0 .or. any(ap /= 0.0_dp)) &
        error stop 'zero pseudoinverse'
    call lstsq(a,b,x,status,rank)
    if (status /= linalg_status_success .or. rank /= 0 .or. any(x /= 0.0_dp)) &
        error stop 'zero minimum norm'

    call svd(a,u(1:3,:),s,vh,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad U shape'
    if (any(s /= 0.0_dp)) error stop 'bad U returned singulars'
    call svd(a,u,s,vh,status,max_sweeps=0)
    if (status /= linalg_status_invalid_shape) error stop 'bad sweep count'
    call pinv(a,ap,status,rcond=-0.1_dp)
    if (status /= linalg_status_invalid_shape) error stop 'bad cutoff'
    a(1,1) = ieee_value(1.0_dp,ieee_quiet_nan)
    call svd(a,u,s,vh,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite SVD'
    if (any(s /= 0.0_dp)) error stop 'SVD output on error'
    call lstsq(a,b,x,status,method='svd')
    if (status /= linalg_status_no_convergence) error stop 'nonfinite lstsq'

    print '(A)', 'real and complex SVD, pseudoinverse, rank deficient lstsq and errors passed'
end program test_linalg_svd
