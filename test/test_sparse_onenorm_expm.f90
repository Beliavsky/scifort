! SPDX-License-Identifier: MIT
! Native regression for the CSR 1-norm estimator and sparse matrix exponential.
program test_sparse_onenorm_expm
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, &
        csr_to_dense, valid_csr, csr_eye
    use scifort_sparse_linalg, only : onenormest, expm, expm_multiply, &
        sparse_linalg_invalid_input, sparse_linalg_breakdown
    implicit none
    type(csr_matrix) :: a, e, bad
    real(dp), allocatable :: v(:), w(:), dense(:,:), y(:)
    real(dp) :: d(3,3), rotation(2,2), expected(2,2), est, nan_value
    integer :: info, i

    d = reshape([1.0_dp,5.0_dp,0.0_dp, &
                 0.0_dp,8.0_dp,-1.0_dp, &
                 0.0_dp,2.0_dp,0.0_dp],[3,3])
    a = csr_from_dense(d)
    est = onenormest(a,t=3,v=v,w=w,info=info)
    call require(info == 0, 'exact 1-norm status')
    call require(abs(est-9.0_dp) < 1.0e-14_dp, 'exact 1-norm')
    call require(abs(sum(abs(w))-est) < 1.0e-14_dp, 'exact witness output')
    call require(abs(sum(abs(v))-1.0_dp) < 1.0e-14_dp, 'exact witness input')
    call require(maxval(abs(matmul(d,v)-w)) < 1.0e-14_dp, 'exact Av')
    est = onenormest(a,t=1,itmax=7,v=v,w=w,info=info)
    call require(info == 0 .and. est <= 9.0_dp*(1.0_dp+1.0e-12_dp), &
        'iterative lower bound')
    call require(est >= 6.0_dp, 'iterative estimate useful')
    call require(abs(sum(abs(w))-est) < 1.0e-12_dp, 'iterative witness output')
    call require(abs(sum(abs(v))-1.0_dp) < 1.0e-14_dp, 'iterative witness input')
    call require(maxval(abs(matmul(d,v)-w)) < 1.0e-14_dp, 'iterative Av')
    est = onenormest(a,t=0,info=info)
    call require(info == sparse_linalg_invalid_input, 'invalid estimator width')
    est = onenormest(a,itmax=0,info=info)
    call require(info == sparse_linalg_invalid_input, 'invalid estimator iterations')
    bad = csr_matrix()
    est = onenormest(bad,info=info)
    call require(info == sparse_linalg_invalid_input, 'unallocated CSR')
    a = csr_from_dense(reshape([real(dp)::],[0,0]))
    est = onenormest(a,v=v,w=w,info=info)
    call require(info == 0 .and. est == 0.0_dp, 'empty norm')
    call require(size(v) == 0 .and. size(w) == 0, 'empty witness')
    e = expm(a,info)
    call require(info == 0 .and. valid_csr(e), 'empty exponential')
    call require(e%n_rows == 0 .and. e%n_cols == 0, 'empty shape')

    rotation = reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2])
    a = csr_from_dense(rotation)
    e = expm(a,info,block_size=1)
    call require(info == 0 .and. valid_csr(e), 'rotation status')
    dense = csr_to_dense(e)
    expected = reshape([cos(1.0_dp),sin(1.0_dp), &
                        -sin(1.0_dp),cos(1.0_dp)],[2,2])
    call require(maxval(abs(dense-expected)) < 2.0e-13_dp, 'rotation values')
    e = expm(a,info,block_size=2)
    call require(info == 0 .and. maxval(abs(csr_to_dense(e)-expected)) < 2.0e-13_dp, &
        'block exponential')
    call expm_multiply(a,[1.0_dp,3.0_dp],y,info)
    call require(info == 0 .and. maxval(abs(matmul(csr_to_dense(e),[1.0_dp,3.0_dp])-y)) &
        < 3.0e-13_dp,'exponential action agrees')

    d = 0.0_dp
    d(1,2) = 2.0_dp
    d(2,3) = -3.0_dp
    a = csr_from_dense(d)
    e = expm(a,info,block_size=2)
    call require(info == 0, 'nilpotent status')
    d = 0.0_dp
    d(1,1) = 1.0_dp
    d(2,2) = 1.0_dp
    d(3,3) = 1.0_dp
    d(1,2) = 2.0_dp
    d(2,3) = -3.0_dp
    d(1,3) = -3.0_dp
    call require(maxval(abs(csr_to_dense(e)-d)) < 5.0e-13_dp,'nilpotent values')

    a = csr_from_dense(reshape([0.0_dp,0.0_dp,1.0e-8_dp,0.0_dp],[2,2]))
    e = expm(a,info,drop_tol=1.0e-7_dp)
    call require(info == 0 .and. size(e%data) == 2, 'exponential cutoff')
    call require(maxval(abs(csr_to_dense(e)-reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], &
        [2,2]))) < 1.0e-13_dp,'cutoff values')
    e = expm(a,info,drop_tol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'bad cutoff')
    e = expm(a,info,block_size=0)
    call require(info == sparse_linalg_invalid_input, 'bad chunk width')
    e = expm(a,info,tol=2.0_dp)
    call require(info == sparse_linalg_invalid_input, 'bad action tolerance')

    block
        integer, parameter :: n = 7000
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        allocate(rows(n),cols(n),vals(n))
        do i = 1,n
            rows(i) = i
            cols(i) = i
            vals(i) = real(mod(i,13)-6,dp)*0.1_dp
        end do
        a = csr_from_coo(n,n,rows,cols,vals)
        e = expm(a,info)
        call require(info == 0 .and. valid_csr(e), 'large diagonal exponential')
        call require(size(e%data) == n, 'large diagonal nonzeros')
        do i = 1,n
            call require(e%indices(i) == i, 'large diagonal coordinates')
            call require(abs(e%data(i)-exp(vals(i))) < 2.0e-15_dp, 'large diagonal values')
        end do
        ! O(nnz), not O(n*n), even for an exact full-width norm request.
        est = onenormest(a,t=n,info=info)
        call require(info == 0, 'large norm status')
        call require(abs(est-maxval(abs(vals))) < 1.0e-15_dp, 'large exact norm')
    end block

    ! A nontrivial upper-bidiagonal sparse matrix has a known exponential.
    ! This exercises the block assembly path beyond tiny dense examples.
    block
        integer, parameter :: n = 160
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        real(dp) :: expected_value, factorial_k
        integer :: k, j
        allocate(rows(2*n-1),cols(2*n-1),vals(2*n-1))
        do i = 1,n
            rows(i) = i
            cols(i) = i
            vals(i) = -0.2_dp
        end do
        do i = 1,n-1
            rows(n+i) = i
            cols(n+i) = i+1
            vals(n+i) = 0.3_dp
        end do
        a = csr_from_coo(n,n,rows,cols,vals)
        e = expm(a,info,block_size=5)
        call require(info == 0 .and. valid_csr(e), 'bidiagonal exponential status')
        dense = csr_to_dense(e)
        call require(size(e%data) < n*n, 'bidiagonal output retains sparsity')
        do k = 0,9
            factorial_k = 1.0_dp
            do j = 2,k
                factorial_k = factorial_k*real(j,dp)
            end do
            expected_value = exp(-0.2_dp)*0.3_dp**k/factorial_k
            call require(abs(dense(10,10+k)-expected_value) < 5.0e-13_dp, &
                'bidiagonal exponential coefficients')
        end do
        call require(maxval(abs(dense(2:n,1))) == 0.0_dp, &
            'upper-triangular exponential structure')
    end block

    a = csr_from_dense(reshape([1000.0_dp,0.0_dp,0.0_dp,1000.0_dp],[2,2]))
    e = expm(a,info)
    call require(info == sparse_linalg_breakdown, 'overflow status')
    a = csr_from_dense(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3]))
    e = expm(a,info)
    call require(info == sparse_linalg_invalid_input, 'rectangular expm rejected')
    est = onenormest(a,info=info)
    call require(info == sparse_linalg_invalid_input, 'rectangular onenormest rejected')
    a = csr_eye(2)
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    e = expm(a,info,drop_tol=nan_value)
    call require(info == sparse_linalg_invalid_input, 'nan cutoff')

    print '(a)', 'test_sparse_onenorm_expm: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Assertion condition.
        character(len=*), intent(in) :: message !! Diagnostic message.
        if (.not. condition) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine require
end program test_sparse_onenorm_expm
