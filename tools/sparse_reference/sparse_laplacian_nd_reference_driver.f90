! SPDX-License-Identifier: MIT
! Machine-readable driver for independent scipy.sparse.linalg.LaplacianNd comparisons.
program sparse_laplacian_nd_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_to_csc
    use scifort_sparse_laplacian_nd, only : laplacian_nd_operator, laplacian_nd, &
        laplacian_nd_matvec, laplacian_nd_tocsr, laplacian_nd_eigenvalues, &
        laplacian_nd_eigenvectors
    use scifort_sparse_linalg, only : spbandwidth, is_sptriangular
    implicit none
    type(laplacian_nd_operator) :: op
    type(csr_matrix) :: a
    type(csc_matrix) :: ac
    real(dp), allocatable :: x(:), y(:), eig(:), vectors(:,:)
    integer, allocatable :: shape(:)
    integer :: cases, case_id, d, bc, m, n, info, lo, hi, lo2, hi2
    logical :: lower, upper, lower2, upper2
    character(len=9), parameter :: bcs(3) = [character(len=9) :: 'dirichlet', 'neumann', 'periodic']

    read(*,*) cases
    do case_id = 1,cases
        read(*,*) d, bc, m
        if (d < 1 .or. bc < 1 .or. bc > 3) error stop 'invalid driver input'
        allocate(shape(d))
        read(*,*) shape
        op = laplacian_nd(shape,bcs(bc),info)
        if (info /= 0) error stop 'laplacian construction'
        n = product(shape)
        allocate(x(n))
        read(*,*) x
        call laplacian_nd_matvec(op,x,y,info)
        if (info /= 0) error stop 'matvec'
        a = laplacian_nd_tocsr(op,info)
        if (info /= 0) error stop 'tocsr'
        ac = csr_to_csc(a)
        call laplacian_nd_eigenvalues(op,eig,info,m=m)
        if (info /= 0) error stop 'eigenvalues'
        call laplacian_nd_eigenvectors(op,vectors,info,m=m)
        if (info /= 0) error stop 'eigenvectors'
        call spbandwidth(a,lo,hi)
        call is_sptriangular(a,lower,upper)
        call spbandwidth(ac,lo2,hi2)
        call is_sptriangular(ac,lower2,upper2)
        write(*,'(A,3(1X,I0))') 'CASE',n,size(a%data),m
        write(*,'(*(ES25.16E3,1X))') y
        write(*,'(*(ES25.16E3,1X))') eig
        write(*,'(*(ES25.16E3,1X))') vectors
        write(*,'(*(I0,1X))') a%indptr
        write(*,'(*(I0,1X))') a%indices
        write(*,'(*(ES25.16E3,1X))') a%data
        write(*,'(8(I0,1X))') lo,hi,merge(1,0,lower),merge(1,0,upper), &
            lo2,hi2,merge(1,0,lower2),merge(1,0,upper2)
        deallocate(shape,x)
    end do
end program sparse_laplacian_nd_reference_driver
