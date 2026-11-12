! SPDX-License-Identifier: MIT
program test_linalg_banded_lu
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, linalg_status_singular
    use scifort_linalg_banded, only : lu_factor_banded, lu_solve_banded, solve_banded
    implicit none
    integer, parameter :: n=5, lower=2, upper=1
    real(dp) :: a(n,n), ab(lower+upper+1,n), lu(upper+2*lower+1,n), b(n), x(n), truth(n)
    real(dp) :: rhs(n,2), sols(n,2), empty(0,n), malformed(lower+upper+1,n-1)
    real(dp), allocatable :: large_ab(:,:), large_lu(:,:), large_rhs(:), large_x(:)
    complex(dp) :: ac(n,n), abc(lower+upper+1,n), luc(upper+2*lower+1,n)
    complex(dp) :: bc(n), xc(n), expected(n), rhsc(n,2), xc2(n,2)
    integer :: ipiv(n), status, i,j
    integer, allocatable :: piv_big(:)

    a=0.0_dp
    a(1,1)=0.01_dp
    a(2,1)=3.0_dp
    a(3,1)=7.0_dp
    a(1,2)=0.2_dp
    a(2,2)=0.02_dp
    a(3,2)=1.0_dp
    a(4,2)=8.0_dp
    a(2,3)=2.0_dp
    a(3,3)=0.05_dp
    a(4,3)=1.0_dp
    a(5,3)=5.0_dp
    a(3,4)=-1.0_dp
    a(4,4)=1.0_dp
    a(5,4)=0.5_dp
    a(4,5)=0.6_dp
    a(5,5)=2.5_dp
    call pack_real(a,ab)
    call lu_factor_banded(ab,lower,upper,lu,ipiv,status)
    if (status /= linalg_status_success) error stop 'real LU factor'
    if (ipiv(1) /= 3 .or. ipiv(2) /= 4) error stop 'LU did not perform expected row pivots'
    truth=[1.0_dp,-2.0_dp,0.5_dp,4.0_dp,-1.5_dp]
    b=matmul(a,truth)
    call lu_solve_banded(lu,lower,upper,ipiv,b,x,status)
    if (status /= linalg_status_success .or. maxval(abs(x-truth)) > 2.0e-12_dp) &
        error stop 'real LU solve'
    rhs(:,1)=b
    rhs(:,2)=matmul(a,[3.0_dp,1.0_dp,2.0_dp,1.0_dp,0.2_dp])
    call lu_solve_banded(lu,lower,upper,ipiv,rhs,sols,status)
    if (status /= linalg_status_success .or. maxval(abs(sols(:,1)-truth)) > 2.0e-12_dp) &
        error stop 'real LU matrix reuse'
    if (maxval(abs(sols(:,2)-[3.0_dp,1.0_dp,2.0_dp,1.0_dp,0.2_dp])) > 2.0e-12_dp) &
        error stop 'real LU second RHS'
    b=matmul(transpose(a),truth)
    call lu_solve_banded(lu,lower,upper,ipiv,b,x,status,trans='T')
    if (status /= linalg_status_success .or. maxval(abs(x-truth)) > 2.0e-12_dp) &
        error stop 'real LU transpose'
    call lu_solve_banded(lu,lower,upper,ipiv,b,x,status,trans='C')
    if (status /= linalg_status_success .or. maxval(abs(x-truth)) > 2.0e-12_dp) &
        error stop 'real LU conjugate transpose'
    call solve_banded(ab,lower,upper,rhs,sols,status)
    if (status /= linalg_status_success) error stop 'original general solve regression'

    ac=cmplx(a,0.0_dp,dp)
    do j=1,n
        do i=max(1,j-upper),min(n,j+lower)
            ac(i,j)=ac(i,j)+cmplx(0.0_dp,0.1_dp*real(2*i-j,dp),dp)
        end do
    end do
    call pack_complex(ac,abc)
    call lu_factor_banded(abc,lower,upper,luc,ipiv,status)
    if (status /= linalg_status_success) error stop 'complex LU factor'
    expected=[(cmplx(real(i,dp)*0.3_dp,real(6-i,dp)*0.2_dp,dp),i=1,n)]
    bc=matmul(ac,expected)
    call lu_solve_banded(luc,lower,upper,ipiv,bc,xc,status)
    if (status /= linalg_status_success .or. maxval(abs(xc-expected)) > 2.0e-12_dp) &
        error stop 'complex LU normal'
    bc=matmul(transpose(ac),expected)
    call lu_solve_banded(luc,lower,upper,ipiv,bc,xc,status,trans='T')
    if (status /= linalg_status_success .or. maxval(abs(xc-expected)) > 2.0e-12_dp) &
        error stop 'complex LU transpose'
    bc=matmul(transpose(conjg(ac)),expected)
    call lu_solve_banded(luc,lower,upper,ipiv,bc,xc,status,trans='C')
    if (status /= linalg_status_success .or. maxval(abs(xc-expected)) > 2.0e-12_dp) &
        error stop 'complex LU conjugate transpose'
    rhsc(:,1)=bc
    rhsc(:,2)=matmul(transpose(conjg(ac)),expected*2.0_dp)
    call lu_solve_banded(luc,lower,upper,ipiv,rhsc,xc2,status,trans='C')
    if (status /= linalg_status_success .or. maxval(abs(xc2(:,2)-2.0_dp*expected)) > 3.0e-12_dp) &
        error stop 'complex LU conjugate multiRHS'

    ! Invalid dimension, pivots, operation, singular factor, and zero-width band.
    call lu_factor_banded(ab,lower,upper,empty,ipiv,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad LU output shape'
    call lu_factor_banded(malformed,lower,upper,lu,ipiv,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad LU input shape'
    call lu_factor_banded(ab,lower,upper,lu,ipiv,status)
    call lu_solve_banded(lu,lower,upper,ipiv,b,x,status,trans='X')
    if (status /= linalg_status_invalid_shape .or. any(x /= 0.0_dp)) &
        error stop 'invalid transpose operation'
    ipiv(1)=0
    call lu_solve_banded(lu,lower,upper,ipiv,b,x,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid pivot'
    ab=0.0_dp
    call lu_factor_banded(ab,lower,upper,lu,ipiv,status)
    if (status /= linalg_status_singular .or. any(ipiv /= 0)) error stop 'singular LU'

    ! Linear storage: 50,000-by-50,000 matrix without n^2 workspace.
    block
        integer, parameter :: big=50000
        allocate(large_ab(3,big),large_lu(4,big),large_rhs(big),large_x(big),piv_big(big))
        large_ab=0.0_dp
        large_ab(2,:)=4.0_dp
        large_ab(1,2:big)=-0.3_dp
        large_ab(3,1:big-1)=-0.2_dp
        large_rhs=3.5_dp
        large_rhs(1)=3.7_dp
        large_rhs(big)=3.8_dp
        call lu_factor_banded(large_ab,1,1,large_lu,piv_big,status)
        if (status /= linalg_status_success) error stop 'large LU factor'
        call lu_solve_banded(large_lu,1,1,piv_big,large_rhs,large_x,status)
        if (status /= linalg_status_success .or. maxval(abs(large_x-1.0_dp)) > 3.0e-13_dp) &
            error stop 'large LU solve'
        call lu_solve_banded(large_lu,1,1,piv_big,large_rhs,large_x,status)
        if (status /= linalg_status_success .or. maxval(abs(large_x-1.0_dp)) > 3.0e-13_dp) &
            error stop 'large LU reuse'
    end block
    print '(A)', 'banded LU factor reuse, N/T/C, pivots, failures, n=50000 passed'
contains
    subroutine pack_real(m, packed)
        real(dp), intent(in) :: m(:,:) !! dense matrix with band lower/upper
        real(dp), intent(out) :: packed(:,:) !! packed band
        integer :: ii,jj
        packed=0.0_dp
        do jj=1,size(m,2)
            do ii=max(1,jj-upper),min(size(m,1),jj+lower)
                packed(upper+1+ii-jj,jj)=m(ii,jj)
            end do
        end do
    end subroutine pack_real
    subroutine pack_complex(m, packed)
        complex(dp), intent(in) :: m(:,:) !! dense matrix with band lower/upper
        complex(dp), intent(out) :: packed(:,:) !! packed band
        integer :: ii,jj
        packed=cmplx(0.0_dp,0.0_dp,dp)
        do jj=1,size(m,2)
            do ii=max(1,jj-upper),min(size(m,1),jj+lower)
                packed(upper+1+ii-jj,jj)=m(ii,jj)
            end do
        end do
    end subroutine pack_complex
end program test_linalg_banded_lu
