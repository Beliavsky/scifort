! SPDX-License-Identifier: MIT
program test_sparse_lobpcg
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_from_dense, csr_matvec
    use scifort_sparse_linalg, only : lobpcg, sparse_linalg_invalid_input
    implicit none
    type(csr_matrix) :: a, b, prec
    real(dp), allocatable :: values(:), vectors(:,:), errors(:)
    real(dp), allocatable :: x0(:,:), dense(:,:), y(:,:)
    real(dp) :: nan_value
    integer :: n, i, steps, info, j

    n = 90
    a = csr_from_coo(n,n,[(i,i=1,n)],[(i,i=1,n)], &
        [(2.0_dp+real(i,dp)/100.0_dp,i=1,n)])
    allocate(x0(n,3))
    do j=1,3
        do i=1,n
            x0(i,j)=sin(real(i*j+2*j,dp)*0.17_dp)+cos(real(i*(j+2),dp)*0.22_dp)
        end do
    end do
    call lobpcg(a,x0,values,vectors,info,largest=.true.,tol=1.0e-8_dp, &
        maxiter=150,iterations=steps,residual_norms=errors)
    print *, 'Largest:',info,steps,values,errors
    call check(info == 0,'largest LOBPCG converges')
    call check(maxval(abs(values-[2.9_dp,2.89_dp,2.88_dp]))<1.0e-7_dp,'largest values')
    call check(maxval(errors) < 1.0e-8_dp,'largest residuals')
    call check(maxval(abs(matmul(transpose(vectors),vectors)-eye(3)))<1.0e-10_dp, &
        'orthogonal eigenvectors')
    call lobpcg(a,x0,values,vectors,info,largest=.false.,tol=1.0e-8_dp, &
        maxiter=150,jacobi=.true.,residual_norms=errors)
    print *, 'Smallest:',info,values,errors
    call check(info == 0,'smallest LOBPCG converges')
    call check(maxval(abs(values-[2.01_dp,2.02_dp,2.03_dp]))<1.0e-7_dp,'smallest values')

    ! A and B both SPD; generalized A*x=lambda*B*x with (nontrivial) constraints.
    n = 50
    a = csr_from_coo(n,n,[(i,i=1,n)],[(i,i=1,n)], &
        [(1.5_dp+real(i,dp)*0.02_dp,i=1,n)])
    b = csr_from_coo(n,n,[(i,i=1,n)],[(i,i=1,n)], &
        [(1.0_dp+real(i,dp)*0.01_dp,i=1,n)])
    deallocate(x0)
    allocate(x0(n,2),y(n,1))
    do j=1,2
        do i=1,n
            x0(i,j)=sin(real(i*j,dp)*0.27_dp)+cos(real(i*(j+1),dp)*0.32_dp)
        end do
    end do
    y=0.0_dp
    y(n,1)=1.0_dp
    call lobpcg(a,x0,values,vectors,info,largest=.true.,tol=1.0e-8_dp, &
        maxiter=150,b=b,y=y,residual_norms=errors)
    print *, 'Generalized:',info,values,errors
    call check(info == 0,'generalized LOBPCG converges')
    call check(maxval(abs(values-[2.48_dp/1.49_dp,2.46_dp/1.48_dp]))<1.0e-7_dp, &
        'generalized values')
    call check(maxval(abs(vectors(n,:)))<1.0e-10_dp,'constraint respected')
    call check(maxval(errors)<1.0e-8_dp,'generalized true residuals')

    ! Repeated eigenvalue: block initial guesses must converge to independent directions.
    n = 120
    a = csr_from_coo(n,n,[(i,i=1,n)],[(i,i=1,n)], &
        [(merge(12.0_dp,1.0_dp+0.02_dp*real(i,dp),i<=2),i=1,n)])
    deallocate(x0,y)
    allocate(x0(n,2))
    do j=1,2
        do i=1,n
            x0(i,j)=sin(real(i*j,dp)*0.37_dp)+cos(real(i*(j+3),dp)*0.19_dp)
        end do
    end do
    call lobpcg(a,x0,values,vectors,info,largest=.true.,tol=1.0e-9_dp, &
        maxiter=110,jacobi=.true.,residual_norms=errors)
    print *, 'Repeated:',info,values,errors
    call check(info == 0,'repeated eigenvalue converges')
    call check(maxval(abs(values-12.0_dp))<1.0e-8_dp,'repeated eigenvalues')
    call check(maxval(abs(matmul(transpose(vectors),vectors)-eye(2)))<1.0e-10_dp, &
        'repeated eigenvectors independent')
    call lobpcg(a,x0,values,vectors,info,maxiter=0,iterations=steps,residual_norms=errors)
    call check(info>0 .and. steps==0 .and. maxval(errors)>1.0e-8_dp, &
        'zero budget reports nonconvergence')
    call lobpcg(a,x0(:,1:1),values,vectors,info,tol=-1.0_dp)
    call check(info==sparse_linalg_invalid_input,'negative tolerance invalid')
    nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
    call lobpcg(a,x0(:,1:1),values,vectors,info,tol=nan_value)
    call check(info==sparse_linalg_invalid_input,'nonfinite tolerance invalid')
    call lobpcg(a,x0(:,1:2),values,vectors,info,y=x0(:,1:1))
    call check(info==sparse_linalg_invalid_input,'dependent initial direction invalid')
    call lobpcg(a,x0(:,1:2),values,vectors,info,m=a,jacobi=.true.)
    call check(info==sparse_linalg_invalid_input,'conflicting preconditioners')
    block
        real(dp) :: asymmetric(3,3)
        asymmetric=0.0_dp
        asymmetric(1,2)=2.0_dp
        a=csr_from_dense(asymmetric)
        call lobpcg(a,reshape([1.0_dp,1.0_dp,1.0_dp],[3,1]), &
            values,vectors,info)
        call check(info==sparse_linalg_invalid_input,'asymmetric A rejected')
    end block
    block
        real(dp) :: small(1,1)
        small(1,1)=4.0_dp
        a=csr_from_dense(small)
        call lobpcg(a,reshape([2.0_dp],[1,1]),values,vectors,info)
        call check(info==0 .and. abs(values(1)-4.0_dp)<1.0e-13_dp, &
            'singleton matrix')
    end block

    ! The large test stores only three thousand CSR diagonal values, not
    ! any dense matrix, and exposes defects in the block work-space shape.
    block
        integer, parameter :: nn = 3600
        integer :: j
        real(dp), allocatable :: start(:,:)
        a=csr_from_coo(nn,nn,[(i,i=1,nn)],[(i,i=1,nn)], &
            [(merge(12.0_dp-real(i-1,dp),1.0_dp+real(mod(i,97),dp)*0.01_dp,i<=3),i=1,nn)])
        allocate(start(nn,3))
        do j=1,3
            do i=1,nn
                start(i,j)=sin(real(i*j,dp)*0.14_dp)+cos(real(i*(j+3),dp)*0.31_dp)
            end do
        end do
        call lobpcg(a,start,values,vectors,info,largest=.true.,maxiter=70, &
            tol=1.0e-7_dp,residual_norms=errors,iterations=steps)
        call check(info==0,'large CSR LOBPCG converges')
        call check(maxval(abs(values-[12.0_dp,11.0_dp,10.0_dp]))<1.0e-7_dp, &
            'large eigenvalues')
        call check(maxval(errors)<1.0e-7_dp,'large true residuals')
    end block
    print '(a)', 'test_sparse_lobpcg: PASS'
contains
    function eye(n) result(a)
        integer, intent(in) :: n
        real(dp) :: a(n,n)
        integer :: i
        a=0.0_dp
        do i=1,n
            a(i,i)=1.0_dp
        end do
    end function eye
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine check
end program test_sparse_lobpcg
