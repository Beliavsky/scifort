! SPDX-License-Identifier: MIT
program test_linalg_dense_cholesky
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_not_positive_definite, linalg_status_no_convergence, linalg_status_singular
    use scifort_linalg_dense, only : cholesky, cho_factor, cho_solve
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    integer, parameter :: n=9
    real(dp) :: r(n,n), a(n,n), a_input(n,n), f(n,n), f2(n,n), &
        truth(n,3), b(n,3), x(n,3), xv(n), id(n,n), bad(n,n)
    complex(dp) :: rc(n,n), ac(n,n), ac_input(n,n), fc(n,n), fc2(n,n), &
        truthc(n,3), bc(n,3), xc(n,3), xcv(n), badc(n,n)
    integer :: i, j, k, status, side
    logical :: lower

    do j=1,n
        do i=1,n
            r(i,j)=0.1_dp*real(mod(3*i+2*j,11)-5,dp)
            rc(i,j)=cmplx(r(i,j),0.02_dp*real(3*i-j,dp),dp)
        end do
    end do
    a=matmul(transpose(r),r)
    ac=matmul(transpose(conjg(rc)),rc)
    id=0.0_dp
    do i=1,n
        a(i,i)=a(i,i)+1.5_dp
        ac(i,i)=ac(i,i)+cmplx(1.5_dp,0.0_dp,dp)
        id(i,i)=1.0_dp
        do k=1,3
            truth(i,k)=real(i*k,dp)/7.0_dp
            truthc(i,k)=cmplx(truth(i,k),real(k-i,dp)/9.0_dp,dp)
        end do
    end do
    b=matmul(a,truth)
    bc=matmul(ac,truthc)

    do side=0,1
        lower=side==1
        a_input=a
        ac_input=ac
        ! The ignored triangle is intentionally incompatible with the matrix.
        do j=1,n
            do i=1,n
                if ((lower .and. i < j) .or. (.not. lower .and. i > j)) then
                    a_input(i,j)=real(17*i+j,dp)
                    ac_input(i,j)=cmplx(real(17*i+j,dp),real(j-i,dp),dp)
                end if
            end do
        end do
        call cho_factor(a_input,f,status,lower=lower)
        if (status /= linalg_status_success) error stop 'real Cholesky factor'
        call cholesky(a_input,f2,status,lower=lower)
        if (status /= linalg_status_success .or. any(f2 /= f)) error stop 'cholesky alias real'
        call cho_factor(ac_input,fc,status,lower=lower)
        if (status /= linalg_status_success) error stop 'complex Cholesky factor'
        call cholesky(ac_input,fc2,status,lower=lower)
        if (status /= linalg_status_success .or. any(fc2 /= fc)) error stop 'cholesky alias complex'
        if (lower) then
            if (maxval(abs(matmul(f,transpose(f))-a)) > 3.0e-12_dp) error stop 'real lower reconstruction'
            if (maxval(abs(matmul(fc,transpose(conjg(fc)))-ac)) > 3.0e-12_dp) &
                error stop 'complex lower reconstruction'
            do j=1,n
                if (any(f(1:j-1,j) /= 0.0_dp)) error stop 'unused upper not zero'
                if (any(fc(1:j-1,j) /= cmplx(0.0_dp,0.0_dp,dp))) error stop 'unused upper complex'
            end do
        else
            if (maxval(abs(matmul(transpose(f),f)-a)) > 3.0e-12_dp) error stop 'real upper reconstruction'
            if (maxval(abs(matmul(transpose(conjg(fc)),fc)-ac)) > 3.0e-12_dp) &
                error stop 'complex upper reconstruction'
            do j=1,n
                if (any(f(j+1:n,j) /= 0.0_dp)) error stop 'unused lower not zero'
                if (any(fc(j+1:n,j) /= cmplx(0.0_dp,0.0_dp,dp))) error stop 'unused lower complex'
            end do
        end if
        call cho_solve(f,b,x,status,lower=lower)
        if (status /= linalg_status_success .or. maxval(abs(x-truth)) > 1.0e-12_dp) &
            error stop 'real multiple right-hand sides'
        call cho_solve(fc,bc,xc,status,lower=lower)
        if (status /= linalg_status_success .or. maxval(abs(xc-truthc)) > 1.0e-12_dp) &
            error stop 'complex multiple right-hand sides'
        call cho_solve(f,b(:,1),xv,status,lower=lower)
        if (status /= linalg_status_success .or. maxval(abs(xv-truth(:,1))) > 1.0e-12_dp) &
            error stop 'real vector right-hand side'
        call cho_solve(fc,bc(:,1),xcv,status,lower=lower)
        if (status /= linalg_status_success .or. maxval(abs(xcv-truthc(:,1))) > 1.0e-12_dp) &
            error stop 'complex vector right-hand side'
        call cho_solve(f, b(:,1:0), x(:,1:0), status, lower=lower)
        if (status /= linalg_status_success) error stop 'empty rhs column support'
    end do

    call cho_factor(a, f(1:n-1,:), status)
    if (status /= linalg_status_invalid_shape) error stop 'bad factor shape'
    call cho_factor(a(1:n-1,:),f,status)
    if (status /= linalg_status_invalid_shape) error stop 'non-square input'
    call cho_solve(f,b,x(:,1:2),status)
    if (status /= linalg_status_invalid_shape) error stop 'bad solution shape'
    call cho_solve(f,b(:,1),xv(1:3),status)
    if (status /= linalg_status_invalid_shape) error stop 'bad vector length'

    bad=id
    bad(4,4)=0.0_dp
    call cho_factor(bad,f,status)
    if (status /= linalg_status_not_positive_definite .or. any(f /= 0.0_dp)) &
        error stop 'semidefinite rejected'
    bad=id
    bad(4,4)=-0.25_dp
    call cho_factor(bad,f,status)
    if (status /= linalg_status_not_positive_definite) error stop 'indefinite rejected'
    badc=cmplx(bad,0.0_dp,dp)
    call cho_factor(badc,fc,status,lower=.true.)
    if (status /= linalg_status_not_positive_definite) error stop 'indefinite complex rejected'
    bad=a
    bad(1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call cho_factor(bad,f,status)
    if (status /= linalg_status_no_convergence .or. any(f /= 0.0_dp)) &
        error stop 'NaN factor rejected'
    badc=ac
    badc(1,1)=cmplx(1.0_dp,ieee_value(0.0_dp,ieee_quiet_nan),dp)
    call cho_factor(badc,fc,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite complex rejected'

    f=id
    f(1,1)=0.0_dp
    call cho_solve(f,b,x,status)
    if (status /= linalg_status_singular .or. any(x /= 0.0_dp)) &
        error stop 'singular triangular factor'
    fc=cmplx(id,0.0_dp,dp)
    fc(1,1)=cmplx(0.0_dp,0.0_dp,dp)
    call cho_solve(fc,bc,xc,status)
    if (status /= linalg_status_singular .or. any(xc /= cmplx(0.0_dp,0.0_dp,dp))) &
        error stop 'complex singular factor'

    print '(A)', 'dense real/complex Cholesky upper/lower, reuse, failures passed'
end program test_linalg_dense_cholesky
