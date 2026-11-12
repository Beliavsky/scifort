! SPDX-License-Identifier: MIT
program test_linalg_banded_eigen
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal, eigvalsh_tridiagonal, &
        eig_banded, eigvals_banded
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    integer :: status, i
    real(dp) :: d(3), e(2), w(3), vals(3), q(3,3), a(3,3), error
    real(dp) :: band(2,3), wlarge(2)
    complex(dp) :: cb(2,3), v(3,3), ca(3,3)
    real(dp), allocatable :: empty(:), empty2(:, :)

    d = 2.0_dp
    e = 1.0_dp
    call eigh_tridiagonal(d,e,w,status,q)
    call check(status == linalg_status_success, 'tridiagonal status')
    call check(maxval(abs(w - [2.0_dp-sqrt(2.0_dp), 2.0_dp, &
        2.0_dp+sqrt(2.0_dp)])) < 1.0e-13_dp, 'tridiagonal values')
    a = 0.0_dp
    a(1,1)=2.0_dp
    a(2,2)=2.0_dp
    a(3,3)=2.0_dp
    a(1,2)=1.0_dp
    a(2,1)=1.0_dp
    a(2,3)=1.0_dp
    a(3,2)=1.0_dp
    do i=1,3
        call check(maxval(abs(matmul(a,q(:,i)) - w(i)*q(:,i))) < 1.0e-13_dp, &
            'tridiagonal residual')
    end do
    call check(maxval(abs(matmul(transpose(q),q) - identity(3))) < 1.0e-13_dp, &
        'tridiagonal orthogonality')
    call eigvalsh_tridiagonal(d,e,vals,status)
    call check(status == linalg_status_success, 'tridiagonal vals only status')
    call check(maxval(abs(vals-w)) < 1.0e-13_dp,'tridiagonal vals-only')

    ! Upper packed band: first row unused at column 1 (can contain NaN).
    band = 0.0_dp
    band(2,:) = 2.0_dp
    band(1,2:3) = 1.0_dp
    band(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call eig_banded(band,w,status,q)
    call check(status == linalg_status_success, 'upper band status')
    call check(maxval(abs(w-vals)) < 1.0e-13_dp,'upper band values')
    call check(maxval(abs(matmul(a,q) - matmul(q,diag(w)))) < 1.0e-13_dp, &
        'upper band residual')
    band = 0.0_dp
    band(1,:) = 2.0_dp
    band(2,1:2)=1.0_dp
    band(2,3) = ieee_value(0.0_dp,ieee_quiet_nan)
    call eigvals_banded(band,vals,status,lower=.true.)
    call check(status == linalg_status_success,'lower band vals status')
    call check(maxval(abs(w-vals)) < 1.0e-13_dp,'lower band values')

    cb = cmplx(0.0_dp,0.0_dp,dp)
    cb(1,:) = cmplx(2.0_dp,0.0_dp,dp)
    cb(2,1:2) = cmplx(0.0_dp,1.0_dp,dp)
    ca = cmplx(0.0_dp,0.0_dp,dp)
    do i=1,3
        ca(i,i) = cmplx(2.0_dp,0.0_dp,dp)
    end do
    ca(2,1)=cmplx(0.0_dp,1.0_dp,dp)
    ca(1,2)=conjg(ca(2,1))
    ca(3,2)=ca(2,1)
    ca(2,3)=conjg(ca(3,2))
    call eig_banded(cb,w,status,v,lower=.true.)
    call check(status == linalg_status_success,'complex Hermitian status')
    call check(maxval(abs(matmul(ca,v)-matmul(v,diag_complex(w)))) < 2.0e-13_dp, &
        'complex Hermitian residual')
    call check(maxval(abs(matmul(conjg(transpose(v)),v) - diag_complex([1.0_dp, &
        1.0_dp,1.0_dp]))) < 2.0e-13_dp,'complex eigenvector orthogonality')

    allocate(empty(0),empty2(0,0))
    call eigh_tridiagonal(empty,empty,empty,status,empty2)
    call check(status == linalg_status_success,'empty tridiagonal')
    call eigh_tridiagonal(d,e,w,status,q,max_iterations=0)
    call check(status == linalg_status_invalid_shape,'invalid iteration limit')
    call eigh_tridiagonal(d,e(1:1),w,status)
    call check(status == linalg_status_invalid_shape,'bad off-diagonal shape')
    d(1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call eigh_tridiagonal(d,e,w,status)
    call check(status == linalg_status_no_convergence,'NaN tridiagonal')
    d = 2.0e150_dp
    e = 1.0e150_dp
    call eigvalsh_tridiagonal(d,e,w,status)
    call check(status == linalg_status_success, 'huge tridiagonal')
    call check(maxval(abs(w/1.0e150_dp - [2.0_dp-sqrt(2.0_dp),2.0_dp, &
        2.0_dp+sqrt(2.0_dp)])) < 1.0e-13_dp,'huge eigenvalues')
    d=2.0e-150_dp
    e=1.0e-150_dp
    call eigvalsh_tridiagonal(d,e,w,status)
    call check(status == linalg_status_success,'tiny tridiagonal')
    call check(maxval(abs(w/1.0e-150_dp - [2.0_dp-sqrt(2.0_dp),2.0_dp, &
        2.0_dp+sqrt(2.0_dp)])) < 1.0e-13_dp,'tiny eigenvalues')
    call eigvalsh_tridiagonal([2.0_dp, 2.0_dp],[1.0_dp],wlarge,status)
    call check(status == linalg_status_success,'2x2 tridiagonal')
    call check(maxval(abs(wlarge - [1.0_dp,3.0_dp])) < 1.0e-13_dp,'2x2 eigenvalues')
    print '(a)', 'banded eigen native tests passed'

contains
    subroutine check(test, label)
        logical, intent(in) :: test !! required assertion
        character(len=*), intent(in) :: label !! diagnostic label
        if (.not. test) then
            print '(a)', 'FAILED: '//label
            error stop 1
        end if
    end subroutine check

    pure function identity(n) result(x)
        integer, intent(in) :: n !! matrix order
        real(dp) :: x(n,n)
        integer :: k
        x = 0.0_dp
        do k=1,n
            x(k,k)=1.0_dp
        end do
    end function identity

    pure function diag(w) result(x)
        real(dp), intent(in) :: w(:) !! diagonal values
        real(dp) :: x(size(w),size(w))
        integer :: k
        x=0.0_dp
        do k=1,size(w)
            x(k,k)=w(k)
        end do
    end function diag

    pure function diag_complex(w) result(x)
        real(dp), intent(in) :: w(:) !! real diagonal values
        complex(dp) :: x(size(w),size(w))
        integer :: k
        x=cmplx(0.0_dp,0.0_dp,dp)
        do k=1,size(w)
            x(k,k)=cmplx(w(k),0.0_dp,dp)
        end do
    end function diag_complex
end program test_linalg_banded_eigen
