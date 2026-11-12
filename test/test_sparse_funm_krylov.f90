! SPDX-License-Identifier: MIT
! Krylov matrix-function actions: Arnoldi, Lanczos, restart, matrix-free API.
program test_sparse_funm_krylov
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : funm_multiply_krylov, expm_multiply, &
        sparse_linalg_invalid_input, sparse_linalg_breakdown
    implicit none
    type(csr_matrix) :: a
    real(dp) :: dense(12,12), b(12), value, nan_value
    real(dp), allocatable :: y(:), z(:), yref(:), zz(:)
    integer :: info, k, ncycles, nproducts, i
    real(dp) :: update

    ! The default Arnoldi path is exact for a full 2x2 skew-symmetric space.
    dense = 0.0_dp
    dense(1,2) = -1.0_dp
    dense(2,1) = 1.0_dp
    a = csr_from_dense(dense(:2,:2))
    call funm_multiply_krylov(a,[1.0_dp,0.0_dp],dense_exponential,y,info, &
                              t=acos(-1.0_dp)/2.0_dp,cycles=ncycles,products=nproducts)
    call require(info == 0 .and. ncycles == 1 .and. nproducts == 2,'rotation status')
    call close(y,[0.0_dp,1.0_dp],3.0e-13_dp,'rotation')

    ! A polynomial matrix function with independent dense reference.
    dense = 0.0_dp
    do i = 1, 5
        dense(i,i) = -0.2_dp*real(i,dp)
    end do
    dense(1,2) = 0.5_dp
    dense(2,3) = -0.9_dp
    dense(3,4) = 0.7_dp
    dense(4,5) = -0.8_dp
    a = csr_from_dense(dense(:5,:5))
    b(:5) = [1.0_dp,2.0_dp,3.0_dp,-2.0_dp,1.0_dp]
    call funm_multiply_krylov(a,b(:5),cubic_minus_linear,y,info,restart_every_m=5)
    call require(info == 0,'polynomial full-space status')
    z = matmul(dense(:5,:5),b(:5))
    yref = matmul(dense(:5,:5),matmul(dense(:5,:5),z)) - z
    call close(y,yref,2.0e-12_dp,'polynomial')

    ! Real Hermitian Lanczos, short Krylov restarts, vs expm_multiply.
    dense = 0.0_dp
    do i = 1, 12
        dense(i,i) = -2.5_dp + 0.06_dp*real(i,dp)
        if (i < 12) then
            dense(i+1,i) = 0.8_dp
            dense(i,i+1) = 0.8_dp
        end if
    end do
    a = csr_from_dense(dense)
    b = 1.0_dp
    b(1) = 4.0_dp
    call funm_multiply_krylov(a,b,dense_exponential,y,info,t=0.25_dp, &
        restart_every_m=4,max_restarts=4,assume_hermitian=.true., &
        cycles=ncycles,products=nproducts,update_norm=update)
    call require(info == 0 .and. ncycles > 1 .and. nproducts > 4,'Lanczos restarted')
    call expm_multiply(a,b,yref,info,t=0.25_dp)
    call require(info == 0,'reference expm status')
    call close(y,yref,2.0e-9_dp,'Lanczos vs expm')
    call funm_multiply_krylov(12,b,dense_product,dense_exponential,z,info, &
        t=0.25_dp,restart_every_m=4,max_restarts=4,assume_hermitian=.true.)
    call require(info == 0,'matrix-free status')
    call close(y,z,2.0e-13_dp,'CSR vs matrix-free')

    ! General (nonsymmetric) Arnoldi with short cycles.
    dense = 0.0_dp
    do i = 1,12
        dense(i,i) = -real(i,dp)/4.0_dp
        if (i < 12) dense(i,i+1) = 0.6_dp
        if (i > 1) dense(i,i-1) = -0.25_dp
    end do
    a = csr_from_dense(dense)
    call funm_multiply_krylov(a,b,dense_exponential,y,info,t=0.25_dp, &
                              restart_every_m=4,max_restarts=4,cycles=ncycles)
    call require(info == 0 .and. ncycles > 1,'Arnoldi restarted')
    call expm_multiply(a,b,yref,info,t=0.25_dp)
    call require(info == 0,'Arnoldi reference status')
    call close(y,yref,5.0e-9_dp,'Arnoldi vs expm')

    ! A limited cycle budget must be recognizable, never reported as success.
    call funm_multiply_krylov(a,b,dense_exponential,y,info,restart_every_m=2, &
                              max_restarts=1,cycles=ncycles)
    call require(info == 1 .and. ncycles == 1,'non-convergence status')

    ! True zero vectors never invoke the user's dense function.
    call funm_multiply_krylov(a,0.0_dp*b,broken_callback,y,info)
    call require(info == 0 .and. all(y == 0.0_dp),'zero right hand side')
    call funm_multiply_krylov(a,b,broken_callback,y,info)
    call require(info == sparse_linalg_breakdown,'callback failure')

    ! Shape/parameter validation and empty matrix.
    call funm_multiply_krylov(a,b(:11),dense_exponential,y,info)
    call require(info == sparse_linalg_invalid_input .and. size(y) == 0,'shape rejection')
    call funm_multiply_krylov(a,b,dense_exponential,y,info,restart_every_m=0)
    call require(info == sparse_linalg_invalid_input,'invalid restart size')
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    call funm_multiply_krylov(a,b,dense_exponential,y,info,t=nan_value)
    call require(info == sparse_linalg_invalid_input,'NaN time')
    call funm_multiply_krylov(0,[real(dp)::],dense_product,dense_exponential,y,info)
    call require(info == 0 .and. size(y) == 0,'empty matrix')
    a = csr_from_dense(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3]))
    call funm_multiply_krylov(a,[1.0_dp,1.0_dp],dense_exponential,y,info)
    call require(info == sparse_linalg_invalid_input,'rectangular matrix')

    ! An oversized user basis should be clamped to n, not allocate O(m**2)
    ! memory when m is far larger than the matrix order.
    dense = 0.0_dp
    dense(1,1) = -0.5_dp
    dense(2,2) = -2.0_dp
    a = csr_from_dense(dense(:2,:2))
    call funm_multiply_krylov(a,[1.0_dp,2.0_dp],dense_exponential,y,info, &
                              restart_every_m=huge(0),max_restarts=3, &
                              cycles=ncycles,products=nproducts)
    call require(info == 0 .and. ncycles == 1 .and. nproducts <= 2,'oversized basis')
    call close(y,[exp(-0.5_dp),2.0_dp*exp(-2.0_dp)],1.0e-12_dp,'oversized basis result')

    ! A nilpotent zero operator terminates immediately by happy breakdown.
    a = csr_from_dense(dense(:2,:2)*0.0_dp)
    call funm_multiply_krylov(a,[3.0_dp,2.0_dp],dense_exponential,y,info, &
                              restart_every_m=1,max_restarts=3, &
                              cycles=ncycles,products=nproducts)
    call require(info == 0 .and. ncycles == 1 .and. nproducts == 1,'zero A happy breakdown')
    call close(y,[3.0_dp,2.0_dp],1.0e-14_dp,'zero A result')

    ! At t=0 the matrix exponential is the identity, regardless of A.
    a = csr_from_dense(dense(:2,:2))
    call funm_multiply_krylov(a,[1.0_dp,-2.0_dp],dense_exponential,y,info, &
                              t=0.0_dp,restart_every_m=2)
    call require(info == 0,'zero time status')
    call close(y,[1.0_dp,-2.0_dp],1.0e-13_dp,'zero time')

    ! Invalid tolerances and failed callbacks do not masquerade as success.
    call funm_multiply_krylov(a,[1.0_dp,2.0_dp],dense_exponential,y,info,rtol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input,'negative tolerance')
    call funm_multiply_krylov(a,[1.0_dp,2.0_dp],dense_exponential,y,info,max_restarts=0)
    call require(info == sparse_linalg_invalid_input,'zero cycle budget')
    call funm_multiply_krylov(a,[1.0_dp,2.0_dp],nan_callback,y,info)
    call require(info == sparse_linalg_breakdown,'nonfinite dense callback')
    call funm_multiply_krylov(2,[1.0_dp,2.0_dp],nan_product,dense_exponential,y,info)
    call require(info == sparse_linalg_breakdown,'nonfinite matrix-free product')

    call test_large_matrix_free()

    print '(a)', 'test_sparse_funm_krylov: PASS'
contains
    subroutine test_large_matrix_free()
        ! The 2048-dimensional operator is never allocated as a dense matrix.
        integer, parameter :: nlarge = 2048
        real(dp) :: rhs(nlarge), action(nlarge), expected(nlarge)
        real(dp), allocatable :: result(:)
        integer :: j, ierr, niters, mv_count

        do j = 1, nlarge
            rhs(j) = sin(real(j,dp)*0.021_dp) + 0.03_dp*cos(real(j,dp)*0.13_dp)
        end do
        call tridiagonal_operator(rhs, action)
        expected = rhs + 0.35_dp*action
        call funm_multiply_krylov(nlarge,rhs,tridiagonal_operator,affine_function,result,ierr, &
             t=0.35_dp, restart_every_m=5,max_restarts=5, &
             cycles=niters,products=mv_count)
        call require(ierr == 0 .and. niters >= 2 .and. mv_count <= 15,'large matrix-free status')
        call close(result,expected,2.0e-11_dp,'2048-vector matrix-free polynomial')
    end subroutine test_large_matrix_free

    subroutine tridiagonal_operator(x, y)
        real(dp), intent(in) :: x(:) !! Matrix-free input vector.
        real(dp), intent(out) :: y(:) !! Minus discrete second difference.
        integer :: j
        y = -2.0_dp*x
        do j = 1,size(x)-1
            y(j) = y(j)+x(j+1)
            y(j+1) = y(j+1)+x(j)
        end do
    end subroutine tridiagonal_operator

    subroutine affine_function(h, fh, status)
        real(dp), intent(in) :: h(:,:) !! Small projected matrix.
        real(dp), intent(out) :: fh(:,:) !! f(h)=I+h.
        integer, intent(out) :: status !! Zero on success.
        integer :: j
        fh = h
        do j=1,size(h,1)
            fh(j,j) = fh(j,j)+1.0_dp
        end do
        status = 0
    end subroutine affine_function

    subroutine nan_product(x, z)
        real(dp), intent(in) :: x(:) !! Input vector.
        real(dp), intent(out) :: z(:) !! Deliberately invalid matvec.
        z = ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine nan_product

    subroutine nan_callback(h, fh, status)
        real(dp), intent(in) :: h(:,:) !! Dense projection.
        real(dp), intent(out) :: fh(:,:) !! Deliberately invalid output.
        integer, intent(out) :: status !! Status returned by callback.
        fh = ieee_value(0.0_dp,ieee_quiet_nan)
        status = 0
    end subroutine nan_callback

    subroutine dense_product(x, z)
        real(dp), intent(in) :: x(:) !! Vector of length 12 (or zero for empty test).
        real(dp), intent(out) :: z(:) !! Dense matrix-vector product.
        if (size(z) > 0) z = matmul(dense(:size(z),:size(x)),x)
    end subroutine dense_product

    subroutine cubic_minus_linear(h, fh, status)
        real(dp), intent(in) :: h(:,:) !! Dense Krylov projection.
        real(dp), intent(out) :: fh(:,:) !! Polynomial f(H)=H**3-H.
        integer, intent(out) :: status !! Zero on completion.
        fh = matmul(h,matmul(h,h))-h
        status = 0
    end subroutine cubic_minus_linear

    subroutine broken_callback(h, fh, status)
        real(dp), intent(in) :: h(:,:) !! Ignored Krylov projection.
        real(dp), intent(out) :: fh(:,:) !! Unused result.
        integer, intent(out) :: status !! Intentional failure.
        fh = 0.0_dp
        status = 42
    end subroutine broken_callback

    subroutine dense_exponential(h, fh, status)
        real(dp), intent(in) :: h(:,:) !! Real square Krylov projection.
        real(dp), intent(out) :: fh(:,:) !! exp(H).
        integer, intent(out) :: status !! Zero for success.
        real(dp) :: mat(size(h,1),size(h,2)), term(size(h,1),size(h,2)), anorm
        integer :: j, i, scale_count
        status = 0
        fh = 0.0_dp
        do i = 1,size(h,1)
            fh(i,i) = 1.0_dp
        end do
        if (size(h,1) == 0) return
        anorm = maxval(sum(abs(h),dim=1))
        scale_count = max(0,ceiling(log(max(1.0_dp,anorm))/log(2.0_dp)))
        mat = h/real(2**scale_count,dp)
        term = fh
        do j = 1,100
            term = matmul(term,mat)/real(j,dp)
            fh = fh+term
            if (maxval(abs(term)) < 1.0e-16_dp) exit
        end do
        do j = 1,scale_count
            fh = matmul(fh,fh)
        end do
    end subroutine dense_exponential

    subroutine require(condition, label)
        logical, intent(in) :: condition !! Condition to assert.
        character(*), intent(in) :: label !! Test name.
        if (.not. condition) then
            print '(2a)', 'FAIL: ',label
            error stop 1
        end if
    end subroutine require

    subroutine close(a0,b0,tol,label)
        real(dp), intent(in) :: a0(:) !! Actual result.
        real(dp), intent(in) :: b0(:) !! Reference result.
        real(dp), intent(in) :: tol !! Absolute tolerance.
        character(*), intent(in) :: label !! Label.
        if (size(a0) /= size(b0)) then
            print *, 'wrong shape: ',label
            error stop 1
        end if
        if (maxval(abs(a0-b0)) > tol) then
            print '(a,a,es13.4)', 'Mismatch: ',label,maxval(abs(a0-b0))
            error stop 1
        end if
    end subroutine close
end program test_sparse_funm_krylov
