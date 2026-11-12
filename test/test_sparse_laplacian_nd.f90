! SPDX-License-Identifier: MIT
program test_sparse_laplacian_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_from_coo, csr_from_dense, &
        csr_to_csc, csr_to_dense, valid_csr
    use scifort_sparse_laplacian_nd, only : laplacian_nd_operator, laplacian_nd, &
        laplacian_nd_matvec, laplacian_nd_matmat, laplacian_nd_tocsr, &
        laplacian_nd_eigenvalues, laplacian_nd_eigenvectors, &
        laplacian_nd_shape, valid_laplacian_nd, laplacian_nd_invalid_input
    use scifort_sparse_linalg, only : spbandwidth, is_sptriangular, sparse_linalg_invalid_input
    implicit none
    type(laplacian_nd_operator) :: op
    type(csr_matrix) :: a, broken
    type(csc_matrix) :: b
    real(dp), allocatable :: x(:), y(:), v(:,:), eig(:), dense(:,:), matx(:,:), maty(:,:)
    real(dp) :: nan_value
    integer, allocatable :: grid_dims(:)
    integer :: bc, i, j, n, status, nshape, lo, hi
    logical :: lower, upper
    character(len=9), parameter :: bc_names(3) = [character(len=9) :: 'dirichlet', 'neumann', 'periodic']

    do bc = 1,3
        op = laplacian_nd([3,4],bc_names(bc),status)
        call require(status == 0 .and. valid_laplacian_nd(op), 'valid operator')
        call laplacian_nd_shape(op,grid_dims,nshape,status)
        call require(status == 0 .and. all(grid_dims == [3,4]) .and. nshape == 12, 'shape')
        n = 12
        x = [(real(i,dp)*0.25_dp,i=1,n)]
        call laplacian_nd_matvec(op,x,y,status)
        call require(status == 0 .and. size(y) == n, 'matvec shape')
        a = laplacian_nd_tocsr(op,status)
        call require(status == 0 .and. valid_csr(a), 'CSR construction')
        dense = csr_to_dense(a)
        call require(maxval(abs(matmul(dense,x)-y)) < 1.0e-12_dp, 'matrix free and CSR')
        call laplacian_nd_eigenvalues(op,eig,status)
        call require(status == 0 .and. size(eig) == n, 'full eigenvalues')
        call require(all(eig(2:) >= eig(:n-1)-1.0e-14_dp), 'ordered eigenvalues')
        call laplacian_nd_eigenvectors(op,v,status)
        call require(status == 0 .and. all(shape(v) == [n,n]), 'full eigenvectors')
        do j = 1,n
            call require(maxval(abs(matmul(dense,v(:,j))-eig(j)*v(:,j))) < 4.0e-12_dp, &
                'eigenpair residual')
        end do
        matx = matmul(transpose(v),v)
        do i = 1,n
            matx(i,i) = matx(i,i)-1.0_dp
        end do
        call require(maxval(abs(matx)) < 4.0e-12_dp, 'orthonormal eigenbasis')
        call laplacian_nd_eigenvalues(op,eig,status,m=3)
        call laplacian_nd_eigenvectors(op,v,status,m=3)
        call require(status == 0 .and. size(eig) == 3 .and. size(v,2) == 3, 'truncated spectrum')
        do j = 1,3
            call require(maxval(abs(matmul(dense,v(:,j))-eig(j)*v(:,j))) < 4.0e-12_dp, &
                'truncated eigenpair residual')
        end do
        deallocate(matx)
        allocate(matx(n,2))
        matx(:,1) = x
        matx(:,2) = 0.5_dp*x
        call laplacian_nd_matmat(op,matx,maty,status)
        call require(status == 0 .and. maxval(abs(maty(:,1)-y)) < 1.0e-12_dp, 'matmat')
        call require(maxval(abs(maty(:,2)-0.5_dp*y)) < 1.0e-12_dp, 'matmat second column')
        deallocate(matx)
    end do

    ! Duplicated periodic neighbors in a length-two dimension are summed.
    op = laplacian_nd([2,2],'periodic')
    a = laplacian_nd_tocsr(op)
    dense = csr_to_dense(a)
    call require(dense(1,2) == 2.0_dp .and. dense(1,3) == 2.0_dp &
        .and. dense(1,1) == -4.0_dp, 'periodic n=2 multiedges')

    ! Length-one axes: match the mathematical/matrix-free SciPy operator.
    do bc = 1,3
        op = laplacian_nd([1,4],bc_names(bc))
        a = laplacian_nd_tocsr(op)
        dense = csr_to_dense(a)
        x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
        call laplacian_nd_matvec(op,x,y)
        call require(maxval(abs(matmul(dense,x)-y)) < 1.0e-12_dp, 'singleton axis')
        call laplacian_nd_eigenvalues(op,eig)
        call laplacian_nd_eigenvectors(op,v)
        do j = 1,4
            call require(maxval(abs(matmul(dense,v(:,j))-eig(j)*v(:,j))) < 1.0e-12_dp, &
                'singleton eigenvalues')
        end do
    end do

    op = laplacian_nd([2,3,4],'neumann')
    a = laplacian_nd_tocsr(op)
    dense = csr_to_dense(a)
    call require(maxval(abs(sum(dense,dim=2))) < 1.0e-12_dp, 'neumann constant nullspace')
    call laplacian_nd_eigenvalues(op,eig,m=1)
    call require(abs(eig(1)) < 1.0e-12_dp, 'neumann largest eigenvalue zero')

    ! Nonfinite input, invalid dimensions and invalid mode.
    op = laplacian_nd([2,0],status=status)
    call require(status == laplacian_nd_invalid_input .and. .not. valid_laplacian_nd(op), 'zero dim')
    op = laplacian_nd([integer::],status=status)
    call require(status == laplacian_nd_invalid_input, 'empty grid')
    op = laplacian_nd([huge(0),2],status=status)
    call require(status == laplacian_nd_invalid_input, 'overflow dimensions')
    op = laplacian_nd([4],boundary_conditions='unknown',status=status)
    call require(status == laplacian_nd_invalid_input, 'invalid mode')
    op = laplacian_nd([3],'dirichlet')
    call laplacian_nd_matvec(op,[1.0_dp],y,status)
    call require(status == laplacian_nd_invalid_input .and. size(y) == 0, 'invalid vector length')
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    call laplacian_nd_matvec(op,[1.0_dp,nan_value,2.0_dp],y,status)
    call require(status == laplacian_nd_invalid_input, 'NaN input')
    call laplacian_nd_eigenvalues(op,eig,status,m=0)
    call require(status == laplacian_nd_invalid_input .and. size(eig) == 0, 'invalid m')
    call laplacian_nd_eigenvectors(op,v,status,m=4)
    call require(status == laplacian_nd_invalid_input .and. size(v) == 0, 'm too large')
    call laplacian_nd_matmat(op,reshape([real(dp)::],[3,0]),maty,status)
    call require(status == 0 .and. all(shape(maty) == [3,0]), 'zero matrix columns')

    ! Structural queries count explicitly stored zeros. Test both CSR and CSC.
    a = csr_from_coo(4,4,[1,2,3,4],[1,1,2,3],[0.0_dp,3.0_dp,2.0_dp,1.0_dp])
    b = csr_to_csc(a)
    call spbandwidth(a,lo,hi,status)
    call require(status == 0 .and. lo == 1 .and. hi == 0, 'CSR lower bandwidth')
    call is_sptriangular(a,lower,upper,status)
    call require(status == 0 .and. lower .and. .not. upper, 'CSR triangularity')
    call spbandwidth(b,lo,hi,status)
    call require(status == 0 .and. lo == 1 .and. hi == 0, 'CSC lower bandwidth')
    call is_sptriangular(b,lower,upper,status)
    call require(status == 0 .and. lower .and. .not. upper, 'CSC triangularity')
    a = csr_from_coo(3,4,[1,2,3],[4,2,1],[0.0_dp,3.0_dp,1.0_dp])
    b = csr_to_csc(a)
    call spbandwidth(a,lo,hi,status)
    call require(status == 0 .and. lo == 2 .and. hi == 3, 'rectangular CSR bandwidth')
    call spbandwidth(b,lo,hi,status)
    call require(status == 0 .and. lo == 2 .and. hi == 3, 'rectangular CSC bandwidth')
    call is_sptriangular(a,lower,upper,status)
    call require(status == 0 .and. .not. lower .and. .not. upper, 'rectangular mixed structure')
    a = csr_from_dense(reshape([real(dp)::],[0,0]))
    call spbandwidth(a,lo,hi,status)
    call require(status == 0 .and. lo == 0 .and. hi == 0, 'empty bandwidth')
    call is_sptriangular(a,lower,upper,status)
    call require(status == 0 .and. lower .and. upper, 'empty triangular')
    broken = a
    broken%indptr(1) = 2
    call spbandwidth(broken,lo,hi,status)
    call require(status == sparse_linalg_invalid_input, 'invalid CSR')

    print '(A)', 'test_sparse_laplacian_nd: PASS'
contains
    subroutine require(condition,message)
        logical, intent(in) :: condition !! Test result.
        character(len=*), intent(in) :: message !! Failure label.
        if (.not. condition) then
            print '(A)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_laplacian_nd
