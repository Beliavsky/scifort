! SPDX-License-Identifier: MIT
program test_linalg_dense_eigen
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_dense_eigen, only : eigh, eigvalsh
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: a(3,3), w(3), u(3,3), part(2), upart(3,2), small(1)
    real(dp) :: empty(0,0), wempty(0), uempty(0,0), nan
    complex(dp) :: c(3,3), v(3,3), vc(3,2)
    real(dp) :: eye(3,3), scale, expect(3)
    integer :: i, status

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    eye = 0.0_dp
    do i = 1, 3
        eye(i,i) = 1.0_dp
    end do

    a = 0.0_dp
    do i = 1, 3
        a(i,i) = 2.0_dp
    end do
    a(2,1) = -1.0_dp
    a(3,2) = -1.0_dp
    ! The upper triangle is ignored by default, including nonfinite entries.
    a(1,2) = nan
    a(2,3) = nan
    a(1,3) = nan
    expect = [2.0_dp-sqrt(2.0_dp), 2.0_dp, 2.0_dp+sqrt(2.0_dp)]
    call eigh(a,w,status,u)
    call check(status == linalg_status_success,'real lower status')
    call check(maxval(abs(w-expect)) < 1.0e-13_dp,'real lower eigenvalues')
    a(1,2) = -1.0_dp
    a(2,3) = -1.0_dp
    a(1,3) = 0.0_dp
    call check(maxval(abs(matmul(a,u)-matmul(u,diag_real(w)))) < 1.0e-13_dp, &
        'real lower residual')
    call check(maxval(abs(matmul(transpose(u),u)-eye)) < 1.0e-13_dp, &
        'real lower orthonormality')

    a(2,1) = nan
    a(3,2) = nan
    call eigvalsh(a,w,status,lower=.false.)
    call check(status == linalg_status_success,'real upper status')
    call check(maxval(abs(w-expect)) < 1.0e-13_dp,'real upper eigenvalues')
    a(2,1) = -1.0_dp
    a(3,2) = -1.0_dp
    call eigh(a,part,status,upart,subset_by_index=[2,3])
    call check(status == linalg_status_success,'real subset status')
    call check(maxval(abs(part-expect(2:3))) < 1.0e-13_dp,'real subset values')
    call check(maxval(abs(matmul(a,upart)-matmul(upart,diag_real(part)))) < 1.0e-13_dp, &
        'real subset residual')
    call eigvalsh(a,part,status,subset_by_index=[1,2])
    call check(status == linalg_status_success,'real values-only subset status')
    call check(maxval(abs(part-expect(1:2))) < 1.0e-13_dp,'real values-only subset')

    c = cmplx(0.0_dp, 0.0_dp, dp)
    do i = 1, 3
        c(i,i) = cmplx(2.0_dp, nan, dp)
    end do
    c(2,1) = cmplx(0.0_dp,1.0_dp,dp)
    c(3,2) = cmplx(0.0_dp,1.0_dp,dp)
    c(1,2) = cmplx(nan,0.0_dp,dp)
    c(2,3) = cmplx(nan,0.0_dp,dp)
    call eigh(c,w,status,v)
    call check(status == linalg_status_success,'complex lower status')
    call check(maxval(abs(w-expect)) < 1.0e-13_dp,'complex eigenvalues')
    c(1,2) = conjg(c(2,1))
    c(2,3) = conjg(c(3,2))
    do i = 1, 3
        c(i,i) = cmplx(2.0_dp,0.0_dp,dp)
    end do
    call check(maxval(abs(matmul(c,v)-matmul(v,diag_complex(w)))) < 2.0e-13_dp, &
        'complex lower residual')
    call check(maxval(abs(matmul(conjg(transpose(v)),v)-cmplx(eye,0.0_dp,dp))) &
        < 2.0e-13_dp,'complex lower unitary')
    call eigh(c,part,status,vc,subset_by_index=[2,3])
    call check(status == linalg_status_success,'complex subset')
    call check(maxval(abs(matmul(c,vc)-matmul(vc,diag_complex(part)))) < 2.0e-13_dp, &
        'complex subset residual')
    call eigvalsh(c,w,status,lower=.false.)
    call check(status == linalg_status_success,'complex upper values')

    scale = 1.0e150_dp
    a = 0.0_dp
    a(1,1) = scale
    a(2,2) = -2.0_dp * scale
    a(3,3) = 4.0_dp * scale
    call eigvalsh(a,w,status)
    call check(status == linalg_status_success,'large scale')
    call check(maxval(abs(w/scale - [-2.0_dp,1.0_dp,4.0_dp])) < 1.0e-13_dp, &
        'large scale eigenvalues')
    call eigh(empty,wempty,status,uempty)
    call check(status == linalg_status_success,'empty matrix')

    call eigvalsh(a,small,status)
    call check(status == linalg_status_invalid_shape,'bad result size')
    call eigh(a,part,status,upart,subset_by_index=[0,1])
    call check(status == linalg_status_invalid_shape,'bad subset low')
    call eigh(a,part,status,upart,subset_by_index=[2,4])
    call check(status == linalg_status_invalid_shape,'bad subset high')
    call eigh(a,part,status,upart,subset_by_index=[3,2])
    call check(status == linalg_status_invalid_shape,'reversed subset')
    call eigh(a,w,status,max_iterations=0)
    call check(status == linalg_status_invalid_shape,'bad iteration limit')
    a(1,1) = nan
    call eigvalsh(a,w,status)
    call check(status == linalg_status_no_convergence,'nonfinite used diagonal')
    print '(a)', 'dense eigh native tests passed'

contains

    subroutine check(condition,message)
        logical, intent(in) :: condition !! expected true
        character(*), intent(in) :: message !! test description
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine check

    function diag_real(x) result(y)
        real(dp), intent(in) :: x(:) !! real diagonal
        real(dp) :: y(size(x),size(x))
        integer :: j
        y = 0.0_dp
        do j=1,size(x)
            y(j,j) = x(j)
        end do
    end function diag_real

    function diag_complex(x) result(y)
        real(dp), intent(in) :: x(:) !! real diagonal
        complex(dp) :: y(size(x),size(x))
        integer :: j
        y = cmplx(0.0_dp,0.0_dp,dp)
        do j=1,size(x)
            y(j,j) = cmplx(x(j),0.0_dp,dp)
        end do
    end function diag_complex
end program test_linalg_dense_eigen
