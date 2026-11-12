! SPDX-License-Identifier: MIT
! Sparse bandwidth and structural triangularity for randomized rectangular CSR/CSC inputs.
program sparse_structure_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_from_coo, csr_to_csc
    use scifort_sparse_linalg, only : spbandwidth, is_sptriangular
    implicit none
    type(csr_matrix) :: a
    type(csc_matrix) :: b
    integer, allocatable :: rr(:), cc(:)
    real(dp), allocatable :: values(:)
    integer :: cases, ic, n, m, nnz, i, lo, hi, lo2, hi2, info
    logical :: lower, upper, lower2, upper2

    read(*,*) cases
    do ic = 1,cases
        read(*,*) n,m,nnz
        allocate(rr(nnz),cc(nnz),values(nnz))
        do i = 1,nnz
            read(*,*) rr(i),cc(i),values(i)
        end do
        a = csr_from_coo(n,m,rr,cc,values,info)
        if (info /= 0) error stop 'invalid COO input'
        b = csr_to_csc(a)
        call spbandwidth(a,lo,hi,info)
        if (info /= 0) error stop 'CSR bandwidth'
        call is_sptriangular(a,lower,upper,info)
        if (info /= 0) error stop 'CSR triangularity'
        call spbandwidth(b,lo2,hi2,info)
        if (info /= 0) error stop 'CSC bandwidth'
        call is_sptriangular(b,lower2,upper2,info)
        if (info /= 0) error stop 'CSC triangularity'
        write(*,'(8(I0,1X))') lo,hi,merge(1,0,lower),merge(1,0,upper), &
            lo2,hi2,merge(1,0,lower2),merge(1,0,upper2)
        deallocate(rr,cc,values)
    end do
end program sparse_structure_reference_driver
