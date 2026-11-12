! SPDX-License-Identifier: MIT
program test_linalg_qr
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_qr, only : qr, qr_factor, qr_unpack, qr_solve, lstsq
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp), allocatable :: a(:, :), q(:, :), r(:, :), af(:, :), tau(:), x(:, :), b(:, :)
    real(dp), allocatable :: truth(:, :), identity(:, :), xv(:), yv(:)
    complex(dp), allocatable :: ac(:, :), qc(:, :), rc(:, :), afc(:, :), tauc(:)
    complex(dp), allocatable :: xc(:, :), bc(:, :), truthc(:, :), xcv(:), ycv(:)
    integer, allocatable :: piv(:)
    integer :: m, n, kmax, mode, variant, i, j, h, status, rank
    real(dp) :: tol

    tol = 5.0e-11_dp
    do mode = 1, 3
        select case(mode)
        case(1)
            m = 11
            n = 7
        case(2)
            m = 5
            n = 9
        case(3)
            m = 8
            n = 8
        end select
        kmax = min(m,n)
        allocate(a(m,n),q(m,kmax),r(kmax,n),af(m,n),tau(kmax),piv(n))
        allocate(ac(m,n),qc(m,kmax),rc(kmax,n),afc(m,n),tauc(kmax))
        allocate(identity(m,m),truth(n,3),truthc(n,3),b(m,3),bc(m,3))
        allocate(x(n,3),xc(n,3),xv(n),xcv(n),yv(m),ycv(m))
        identity = 0.0_dp
        do i = 1, m
            identity(i,i) = 1.0_dp
        end do
        do j = 1, n
            do i = 1, m
                a(i,j) = sin(real(i*(2*j+1),dp)) + 0.31_dp*cos(real(2*i*j+j,dp))
                ac(i,j) = cmplx(a(i,j),0.3_dp*sin(real(i*(3*j-1),dp)),dp)
            end do
            do h = 1, 3
                truth(j,h) = cos(real(j+3*h,dp))
                truthc(j,h) = cmplx(truth(j,h),sin(real(2*j+h,dp)),dp)
            end do
        end do
        a(:,1) = 0.05_dp*a(:,1)
        ac(:,1) = 0.05_dp*ac(:,1)
        if (m < n) then
            ! Known minimum-norm solution in the row space of A.
            truth = matmul(transpose(a),reshape([(sin(real(i,dp)),i=1,m*3)], [m,3]))
            truthc = matmul(transpose(conjg(ac)), &
                cmplx(reshape([(sin(real(i,dp)),i=1,m*3)],[m,3]), &
                reshape([(cos(real(i,dp)),i=1,m*3)],[m,3]),dp))
        end if
        b = matmul(a,truth)
        bc = matmul(ac,truthc)
        do variant = 1, 2
            if (variant == 1) then
                call qr(a,q,r,status)
                if (status /= linalg_status_success) error stop 'real QR'
                if (maxval(abs(matmul(q,r)-a)) > tol) error stop 'real QR reconstruction'
                if (maxval(abs(matmul(transpose(q),q)-identity(1:kmax,1:kmax))) > tol) &
                    error stop 'real Q orthogonality'
                call qr_factor(a,af,tau,status,piv)
                if (status /= linalg_status_success) error stop 'real pivoted factor'
                call qr_unpack(af,tau,q,r,status)
                if (status /= linalg_status_success) error stop 'real factor unpack'
                if (maxval(abs(matmul(q,r)-a(:,piv))) > tol) error stop 'real pivot QR'
                if (m >= n) then
                    call qr_solve(af,tau,piv,b,x,status,rank)
                    if (status /= linalg_status_success .or. rank /= n) error stop 'real reuse QR'
                    if (maxval(abs(x-truth)) > tol) error stop 'real QR solve values'
                    call qr_solve(af,tau,piv,2.0_dp*b,x,status,rank)
                    if (status /= linalg_status_success .or. maxval(abs(x-2.0_dp*truth)) > 2.0_dp*tol) &
                        error stop 'real reused second rhs'
                end if
                call lstsq(a,b,x,status,rank)
                if (status /= linalg_status_success .or. rank /= kmax) then
                    print *, 'real lstsq failure',mode,status,rank,kmax
                    error stop 'real lstsq'
                end if
                if (maxval(abs(x-truth)) > tol) then
                    print *, 'real lstsq mismatch',mode,maxval(abs(x-truth))
                    error stop 'real lstsq values'
                end if
                call lstsq(a,b(:,1),xv,status)
                if (status /= linalg_status_success .or. maxval(abs(xv-truth(:,1))) > tol) &
                    error stop 'real lstsq vector'
                yv = b(:,1)
            else
                call qr(ac,qc,rc,status)
                if (status /= linalg_status_success) error stop 'complex QR'
                if (maxval(abs(matmul(qc,rc)-ac)) > tol) error stop 'complex QR reconstruction'
                if (maxval(abs(matmul(transpose(conjg(qc)),qc)- &
                    cmplx(identity(1:kmax,1:kmax),0.0_dp,dp))) > tol) &
                    error stop 'complex Q orthogonality'
                call qr_factor(ac,afc,tauc,status,piv)
                if (status /= linalg_status_success) error stop 'complex pivoted factor'
                call qr_unpack(afc,tauc,qc,rc,status)
                if (status /= linalg_status_success) error stop 'complex factor unpack'
                if (maxval(abs(matmul(qc,rc)-ac(:,piv))) > tol) error stop 'complex pivot QR'
                if (m >= n) then
                    call qr_solve(afc,tauc,piv,bc,xc,status,rank)
                    if (status /= linalg_status_success .or. rank /= n) error stop 'complex reuse QR'
                    if (maxval(abs(xc-truthc)) > tol) error stop 'complex QR solve values'
                    call qr_solve(afc,tauc,piv,2.0_dp*bc,xc,status,rank)
                    if (status /= linalg_status_success .or. maxval(abs(xc-2.0_dp*truthc)) > 2.0_dp*tol) &
                        error stop 'complex reused second rhs'
                end if
                call lstsq(ac,bc,xc,status,rank)
                if (status /= linalg_status_success .or. rank /= kmax) error stop 'complex lstsq'
                if (maxval(abs(xc-truthc)) > tol) error stop 'complex lstsq values'
                call lstsq(ac,bc(:,1),xcv,status)
                if (status /= linalg_status_success .or. maxval(abs(xcv-truthc(:,1))) > tol) &
                    error stop 'complex lstsq vector'
                ycv = bc(:,1)
            end if
        end do
        ! Test complete Q and R, including the trailing identity subspace.
        deallocate(q,r,qc,rc)
        allocate(q(m,m),r(m,n),qc(m,m),rc(m,n))
        call qr(a,q,r,status,piv)
        if (status /= linalg_status_success) error stop 'real complete QR'
        if (maxval(abs(matmul(q,r)-a(:,piv))) > tol) error stop 'complete QR reconstruction'
        if (maxval(abs(matmul(transpose(q),q)-identity)) > tol) error stop 'complete Q orthogonality'
        call qr(ac,qc,rc,status,piv)
        if (status /= linalg_status_success) error stop 'complex complete QR'
        if (maxval(abs(matmul(qc,rc)-ac(:,piv))) > tol) error stop 'complex complete QR reconstruction'
        if (maxval(abs(matmul(transpose(conjg(qc)),qc)-cmplx(identity,0.0_dp,dp))) > tol) &
            error stop 'complex complete Q orthogonality'
        deallocate(a,q,r,af,tau,x,b,truth,identity,xv,yv,piv)
        deallocate(ac,qc,rc,afc,tauc,xc,bc,truthc,xcv,ycv)
    end do

    allocate(a(5,3),af(5,3),tau(3),piv(3),b(5,2),x(3,2))
    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 1.0_dp
    a(:,3) = a(:,2)
    b = 0.0_dp
    call qr_factor(a,af,tau,status,piv)
    if (status /= linalg_status_success) error stop 'rank deficient factor'
    call qr_solve(af,tau,piv,b,x,status,rank)
    if (status /= linalg_status_singular .or. rank /= 2) error stop 'reused rank deficiency'
    piv(2) = piv(1)
    call qr_solve(af,tau,piv,b,x,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid pivots rejected'
    call lstsq(a,b,x,status,rank,method='qr')
    if (status /= linalg_status_singular .or. rank /= 2 .or. any(x /= 0.0_dp)) &
        error stop 'rank deficient rejected by QR'
    call lstsq(a,b,x,status,rank)
    if (status /= linalg_status_success .or. rank /= 2 .or. any(x /= 0.0_dp)) &
        error stop 'rank deficient handled by SVD fallback' 
    call lstsq(a,b,x,status,rcond=-0.1_dp)
    if (status /= linalg_status_invalid_shape) error stop 'negative rcond rejected'
    call qr_factor(a,af,tau(1:2),status)
    if (status /= linalg_status_invalid_shape) error stop 'tau shape rejected'
    a(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call qr_factor(a,af,tau,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite QR rejected'
    print '(A)', 'real/complex reduced and complete QR, pivoting, tall/wide lstsq, errors passed'
end program test_linalg_qr
