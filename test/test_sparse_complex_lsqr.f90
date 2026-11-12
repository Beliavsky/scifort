! SPDX-License-Identifier: MIT
program test_sparse_complex_lsqr
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc, &
        csr_complex_from_coo, csr_complex_matvec
    use scifort_sparse_linalg_complex_lsqr, only : lsqr_complex, &
        complex_lsqr_invalid_input
    implicit none

    type(csr_complex_matrix) :: a, under, degenerate, empty, large
    type(csc_complex_matrix) :: csc
    complex(dp) :: d(6,4), u(2,3), exact(4), center(4), b(6), &
        grad(4), b_under(2), b_bad(6), h(4,4)
    complex(dp), allocatable :: x(:), x_csc(:), z(:), target(:), rhs_large(:), values(:)
    complex(dp), allocatable :: zero(:,:)
    real(dp), allocatable :: variance(:)
    integer, allocatable :: rcoords(:), ccoords(:)
    real(dp) :: r1, r2, ar, an, cond, xn, lambda
    integer :: info, iterations, i, j

    do j = 1,4
        do i = 1,6
            d(i,j) = cmplx(real(mod(3*i+2*j,11)-5,dp)/7.0_dp, &
                real(mod(5*i+j*j,13)-6,dp)/9.0_dp,dp)
        end do
    end do
    d(1:4,:) = d(1:4,:) + 3.0_dp*identity4()
    do i = 1,4
        exact(i) = cmplx(real(i,dp)/3.0_dp,real(i-2,dp)/5.0_dp,dp)
        center(i) = cmplx(-real(i,dp)/8.0_dp,real(i+1,dp)/10.0_dp,dp)
    end do
    b = matmul(d,exact)
    a = csr_complex_from_dense(d,info)
    call require(info == 0,'create rectangular CSR')
    csc = csr_complex_to_csc(a,info)
    call require(info == 0,'create rectangular CSC')
    call lsqr_complex(a,b,x,info,itn=iterations,r1norm=r1,r2norm=r2, &
        atol=1.0e-12_dp,btol=1.0e-12_dp,iter_lim=200,arnorm=ar, &
        anorm=an,acond=cond,xnorm=xn,var=variance,calc_var=.true.)
    call require(info == 1 .or. info == 2,'consistent LSQR stop')
    call require(iterations > 0 .and. iterations <= 200,'iteration count')
    call require(maxval(abs(x-exact)) < 1.0e-10_dp,'complex tall system solution')
    call require(r1 < 1.0e-10_dp .and. r2 < 1.0e-10_dp,'true residuals')
    call require(ar < 1.0e-9_dp,'normal equation residual')
    call require(an >= 0.0_dp .and. cond >= 0.0_dp,'estimated norms')
    call require(abs(xn-norm2(abs(x))) < 1.0e-12_dp,'solution norm')
    call require(size(variance) == 4 .and. all(variance >= 0.0_dp),'variance diagonal')

    call lsqr_complex(csc,b,x_csc,info,atol=1.0e-12_dp,btol=1.0e-12_dp, &
        iter_lim=200)
    call require(info == 1 .or. info == 2,'CSC LSQR stop')
    call require(maxval(abs(x_csc-exact)) < 1.0e-10_dp,'CSC least squares')

    ! The minimizer of ||Ax-b||^2 + lambda^2 ||x-x0||^2 satisfies
    ! A^H(b-Ax) - lambda^2(x-x0) = 0.
    lambda = 0.6_dp
    b(6) = b(6) + cmplx(0.2_dp,-0.5_dp,dp)
    call lsqr_complex(a,b,x,info,damp=lambda,x0=center,atol=1.0e-12_dp, &
        btol=1.0e-12_dp,iter_lim=200,r1norm=r1,r2norm=r2,arnorm=ar)
    call require(info == 1 .or. info == 2,'damped LSQR stop')
    grad = matmul(transpose(conjg(d)),b-matmul(d,x)) - lambda**2*(x-center)
    call require(norm2(abs(grad)) < 1.0e-9_dp,'Tikhonov stationarity')
    call require(abs(r1-norm2(abs(b-matmul(d,x)))) < 1.0e-11_dp,'damped r1')
    call require(abs(r2-hypot(r1,lambda*norm2(abs(x-center)))) < 1.0e-11_dp,'damped r2')
    call require(abs(ar-norm2(abs(grad))) < 1.0e-11_dp,'damped normal residual')
    h = matmul(transpose(conjg(d)),d)
    do i = 1,4
        h(i,i) = h(i,i)+lambda**2
    end do
    call require(norm2(abs(matmul(h,x)- &
        (matmul(transpose(conjg(d)),b)+lambda**2*center))) < 1.0e-9_dp, &
        'independent regularized normal equations')

    u(1,:) = [cmplx(1.0_dp,1.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),cmplx(1.0_dp,-2.0_dp,dp)]
    u(2,:) = [cmplx(0.0_dp,0.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp),cmplx(1.0_dp,1.0_dp,dp)]
    b_under = [cmplx(0.6_dp,1.0_dp,dp),cmplx(-1.0_dp,0.7_dp,dp)]
    under = csr_complex_from_dense(u,info)
    call require(info == 0,'underdetermined CSR')
    call lsqr_complex(under,b_under,x,info,atol=1.0e-12_dp,btol=1.0e-12_dp, &
        iter_lim=40,r1norm=r1)
    call require(info == 1 .or. info == 2,'underdetermined LSQR stop')
    call require(r1 < 1.0e-11_dp,'underdetermined true residual')
    call require(abs(dot_product(u(1,:),x)) > 0.0_dp,'complex output nontrivial')
    ! The minimum-norm LSQR solution is orthogonal to the nullspace of A.
    z = [u(1,2)*u(2,3)-u(1,3)*u(2,2), &
        u(1,3)*u(2,1)-u(1,1)*u(2,3), &
        u(1,1)*u(2,2)-u(1,2)*u(2,1)]
    call require(abs(dot_product(z,x)) < 1.0e-10_dp,'minimum-norm nullspace orthogonality')

    ! Zero columns, zero rows and exact x0 are all valid edge cases.
    allocate(zero(3,0))
    degenerate = csr_complex_from_dense(zero,info)
    call require(info == 0,'zero-column constructor')
    call lsqr_complex(degenerate,[cmplx(1.0_dp,1.0_dp,dp), &
        cmplx(0.0_dp,0.0_dp,dp),cmplx(2.0_dp,0.0_dp,dp)],x,info,r1norm=r1)
    call require(info == 2 .and. size(x) == 0,'zero columns')
    deallocate(zero)
    allocate(zero(0,3))
    empty = csr_complex_from_dense(zero,info)
    call require(info == 0,'zero-row constructor')
    call lsqr_complex(empty,[complex(dp)::],x,info,r1norm=r1)
    call require(info == 0 .and. size(x) == 3 .and. r1 == 0.0_dp,'zero rows')
    call lsqr_complex(a,b-b,x,info,itn=iterations)
    call require(info == 0 .and. iterations == 0,'exact zero-residual initial solution')
    call lsqr_complex(a,b,x,info,iter_lim=0,itn=iterations)
    call require(info == 7 .and. iterations == 0,'zero iteration cap')

    call lsqr_complex(a,b,x,info,damp=-1.0_dp)
    call require(info == complex_lsqr_invalid_input .and. size(x) == 0,'reject negative damping')
    call lsqr_complex(a,b,x,info,atol=-0.1_dp)
    call require(info == complex_lsqr_invalid_input,'reject negative tolerance')
    call lsqr_complex(a,b,x,info,iter_lim=-1)
    call require(info == complex_lsqr_invalid_input,'reject negative iteration cap')
    call lsqr_complex(a,b,x,info,x0=center(:3))
    call require(info == complex_lsqr_invalid_input,'reject wrong center size')
    call lsqr_complex(a,b(:5),x,info)
    call require(info == complex_lsqr_invalid_input,'reject wrong RHS size')
    b_bad = b
    b_bad(2) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call lsqr_complex(a,b_bad,x,info)
    call require(info == complex_lsqr_invalid_input,'reject nonfinite RHS')

    ! A large genuinely sparse tall system; construct only COO triplets and CSR.
    allocate(rcoords(3*6000),ccoords(3*6000),values(3*6000),target(120))
    do i = 1,6000
        rcoords(3*i-2:3*i) = i
        ccoords(3*i-2) = mod(i-1,120)+1
        ccoords(3*i-1) = mod(7*i+5,120)+1
        ccoords(3*i) = mod(13*i+17,120)+1
        values(3*i-2) = cmplx(2.0_dp,0.5_dp,dp)
        values(3*i-1) = cmplx(-0.4_dp,0.25_dp,dp)
        values(3*i) = cmplx(0.2_dp,-0.45_dp,dp)
    end do
    do j = 1,120
        target(j) = cmplx(sin(real(j,dp)),cos(real(j,dp))/2.0_dp,dp)
    end do
    large = csr_complex_from_coo(6000,120,rcoords,ccoords,values,info)
    call require(info == 0,'large COO to CSR')
    rhs_large = csr_complex_matvec(large,target,info)
    call require(info == 0,'large sparse matvec')
    call lsqr_complex(large,rhs_large,x,info,atol=1.0e-11_dp,btol=1.0e-11_dp, &
        iter_lim=720,r1norm=r1,itn=iterations)
    call require(info == 1 .or. info == 2,'large sparse LSQR convergence')
    call require(iterations > 0 .and. iterations < 720,'large sparse iteration bound')
    call require(maxval(abs(x-target)) < 5.0e-8_dp,'large sparse solution accuracy')
    call require(r1 < 1.0e-7_dp,'large sparse residual')

    print *, 'PASS complex CSR/CSC LSQR tests' 
contains
    function identity4() result(mat)
        complex(dp) :: mat(4,4)
        integer :: idx
        mat = cmplx(0.0_dp,0.0_dp,dp)
        do idx = 1,4
            mat(idx,idx) = cmplx(1.0_dp,0.0_dp,dp)
        end do
    end function identity4

    subroutine require(condition, message)
        logical, intent(in) :: condition !! Condition expected true.
        character(len=*), intent(in) :: message !! Diagnostic label.
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program test_sparse_complex_lsqr
