! SPDX-License-Identifier: MIT
program test_linalg_dense
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_dense, only : solve_triangular, lu_factor, lu_solve
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    integer, parameter :: n=7
    real(dp) :: a(n,n), tri(n,n), op(n,n), lu(n,n), b(n), x(n), exact(n)
    real(dp) :: bm(n,2), xm(n,2), bad(n,n)
    complex(dp) :: ac(n,n), tric(n,n), opc(n,n), luc(n,n), bc(n), xc(n), ec(n)
    complex(dp) :: bcm(n,2), xcm(n,2)
    integer :: piv(n), status, i, j, mode, side, u
    character(len=1), parameter :: ops(3) = ['N','T','C']

    do j=1,n
        exact(j)=real(j,dp)/3.0_dp - 1.0_dp
        ec(j)=cmplx(exact(j),real(j-4,dp)*0.1_dp,dp)
        do i=1,n
            a(i,j)=0.02_dp*real(2*i-j,dp)
            ac(i,j)=cmplx(a(i,j),0.01_dp*real(3*i+j,dp),dp)
        end do
        a(j,j)=3.0_dp+0.1_dp*real(j,dp)
        ac(j,j)=cmplx(a(j,j),0.2_dp*real(j,dp),dp)
    end do

    do side=0,1
        do u=0,1
            tri=a
            tric=ac
            do j=1,n
                do i=1,n
                    if ((side == 1 .and. i < j) .or. (side == 0 .and. i > j)) then
                        tri(i,j)=0.0_dp
                        tric(i,j)=cmplx(0.0_dp,0.0_dp,dp)
                    end if
                end do
                if (u == 1) then
                    tri(j,j)=1.0_dp
                    tric(j,j)=cmplx(1.0_dp,0.0_dp,dp)
                    a(j,j)=0.0_dp
                    ac(j,j)=cmplx(0.0_dp,0.0_dp,dp)
                end if
            end do
            do mode=1,3
                op=tri
                opc=tric
                if (mode == 2 .or. mode == 3) then
                    op=transpose(tri)
                    opc=transpose(tric)
                end if
                if (mode == 3) opc=transpose(conjg(tric))
                b=matmul(op,exact)
                bc=matmul(opc,ec)
                call solve_triangular(a,b,x,status,lower=side==1, &
                    trans=ops(mode),unit_diagonal=u==1)
                if (status /= linalg_status_success .or. maxval(abs(x-exact)) > 1.0e-12_dp) &
                    error stop 'real triangular vector'
                call solve_triangular(ac,bc,xc,status,lower=side==1, &
                    trans=ops(mode),unit_diagonal=u==1)
                if (status /= linalg_status_success .or. maxval(abs(xc-ec)) > 1.0e-12_dp) &
                    error stop 'complex triangular vector'
                bm(:,1)=b
                bm(:,2)=b*2.0_dp
                call solve_triangular(a,bm,xm,status,lower=side==1, &
                    trans=ops(mode),unit_diagonal=u==1)
                if (status /= linalg_status_success .or. maxval(abs(xm(:,2)-2.0_dp*exact)) > 1.0e-12_dp) &
                    error stop 'real triangular multiple RHS'
                bcm(:,1)=bc
                bcm(:,2)=bc*2.0_dp
                call solve_triangular(ac,bcm,xcm,status,lower=side==1, &
                    trans=ops(mode),unit_diagonal=u==1)
                if (status /= linalg_status_success .or. maxval(abs(xcm(:,2)-2.0_dp*ec)) > 1.0e-12_dp) &
                    error stop 'complex triangular multiple RHS'
            end do
            if (u == 1) then
                do j=1,n
                    a(j,j)=3.0_dp+0.1_dp*real(j,dp)
                    ac(j,j)=cmplx(3.0_dp+0.1_dp*real(j,dp),0.2_dp*real(j,dp),dp)
                end do
            end if
        end do
    end do

    ! Force multiple pivots and compare normal and transposed dense LU solves.
    a(1,1)=0.0_dp
    a(3,1)=10.0_dp
    ac=cmplx(a,0.0_dp,dp)
    do j=1,n
        do i=1,n
            ac(i,j)=ac(i,j)+cmplx(0.0_dp,real(2*i-j,dp)*0.08_dp,dp)
        end do
    end do
    call lu_factor(a,lu,piv,status)
    if (status /= linalg_status_success .or. piv(1) /= 3) &
        error stop 'real dense LU row pivoting'
    do mode=1,3
        op=a
        if (mode /= 1) op=transpose(a)
        b=matmul(op,exact)
        call lu_solve(lu,piv,b,x,status,trans=ops(mode))
        if (status /= linalg_status_success .or. maxval(abs(x-exact)) > 1.0e-11_dp) &
            error stop 'real dense LU solve'
        bm(:,1)=b
        bm(:,2)=2.0_dp*b
        call lu_solve(lu,piv,bm,xm,status,trans=ops(mode))
        if (status /= linalg_status_success .or. maxval(abs(xm(:,2)-2.0_dp*exact)) > 2.0e-11_dp) &
            error stop 'real dense LU reuse'
    end do
    call lu_factor(ac,luc,piv,status)
    if (status /= linalg_status_success .or. piv(1) /= 3) error stop 'complex dense LU factor'
    do mode=1,3
        opc=ac
        if (mode == 2) opc=transpose(ac)
        if (mode == 3) opc=transpose(conjg(ac))
        bc=matmul(opc,ec)
        call lu_solve(luc,piv,bc,xc,status,trans=ops(mode))
        if (status /= linalg_status_success .or. maxval(abs(xc-ec)) > 1.0e-11_dp) &
            error stop 'complex dense LU solve'
        bcm(:,1)=bc
        bcm(:,2)=2.0_dp*bc
        call lu_solve(luc,piv,bcm,xcm,status,trans=ops(mode))
        if (status /= linalg_status_success .or. maxval(abs(xcm(:,2)-2.0_dp*ec)) > 2.0e-11_dp) &
            error stop 'complex dense LU multiRHS'
    end do

    call solve_triangular(a,b,x,status,trans='X')
    if (status /= linalg_status_invalid_shape .or. any(x /= 0.0_dp)) &
        error stop 'invalid triangle op'
    call lu_solve(lu,piv,b,x,status,trans='X')
    if (status /= linalg_status_invalid_shape) error stop 'invalid LU op'
    call solve_triangular(a,b,x(1:3),status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid vector shape'
    call solve_triangular(a(1:3,1:2),b,x,status)
    if (status /= linalg_status_invalid_shape) error stop 'non-square triangular input'
    call lu_factor(a,lu,piv,status)
    if (status /= linalg_status_success) error stop 'refactor before pivot failure tests'
    piv(1)=0
    call lu_solve(lu,piv,b,x,status)
    if (status /= linalg_status_invalid_shape) error stop 'malformed LU pivot index'
    piv(1)=3
    lu(n,n)=0.0_dp
    call lu_solve(lu,piv,b,x,status)
    if (status /= linalg_status_singular) error stop 'singular LU solve'
    call lu_factor(a,lu,piv(1:2),status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid pivot shape'
    bad=0.0_dp
    call lu_factor(bad,lu,piv,status)
    if (status /= linalg_status_singular .or. any(piv /= 0)) error stop 'singular factor'
    call solve_triangular(bad,b,x,status)
    if (status /= linalg_status_singular) error stop 'singular triangular'
    bad=a
    bad(1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call lu_factor(bad,lu,piv,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite LU input'
    call solve_triangular(bad,b,x,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite triangular input'

    print '(A)', 'dense triangular and LU N/T/C, unit diagonal, pivots, failures passed'
end program test_linalg_dense
