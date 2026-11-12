! SPDX-License-Identifier: MIT
! Simple stdin/stdout batch driver for independent SciPy N-D COO comparisons.
program sparse_coo_nd_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse_coo, only : coo_matrix
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, &
        coo_nd_to_dense_flat, coo_nd_sum_duplicates, coo_nd_transpose, &
        coo_nd_reshape, coo_nd_from_dense_flat, coo_nd_to_coo, coo_nd_from_coo, &
        coo_nd_get, valid_coo_nd
    implicit none
    type(coo_nd_array) :: a, b, c, t, z
    type(coo_matrix) :: old
    integer, allocatable :: dims(:), perm(:), target(:), coords(:,:)
    real(dp), allocatable :: data(:), vals(:), flat(:)
    real(dp) :: looked_up
    integer :: ncases, trial, ndim, nnz, ndim_target, layout_flag, ierr, k, j
    integer :: roundtrip
    character(len=1) :: order

    read(*,*) ncases
    do trial = 1, ncases
        read(*,*) ndim, nnz, ndim_target, layout_flag
        allocate(dims(ndim),perm(ndim),target(ndim_target),coords(ndim,nnz),data(nnz))
        read(*,*) dims
        read(*,*) target
        read(*,*) perm
        do k = 1, nnz
            read(*,*) coords(:,k),data(k)
        end do
        a = coo_nd_from_triplets(dims,coords,data,ierr)
        if (ierr /= 0 .or. .not. valid_coo_nd(a)) error stop 'failed construction'
        flat = coo_nd_to_dense_flat(a,ierr)
        if (ierr /= 0) error stop 'dense conversion'
        write(*,*) 'CASE',trial,ndim,nnz
        call write_dense('BASE',flat)
        t = coo_nd_transpose(a,perm,ierr)
        if (ierr /= 0) error stop 'transpose'
        vals = coo_nd_to_dense_flat(t,ierr)
        if (ierr /= 0) error stop 'transpose dense'
        call write_dense('TRANSPOSE',vals)
        order = 'C'
        if (layout_flag == 1) order = 'F'
        b = coo_nd_reshape(a,target,order,ierr)
        if (ierr /= 0) error stop 'reshape'
        vals = coo_nd_to_dense_flat(b,ierr)
        if (ierr /= 0) error stop 'reshape dense'
        call write_dense('RESHAPE',vals)
        c = a
        call coo_nd_sum_duplicates(c,ierr)
        if (ierr /= 0) error stop 'sum_duplicates'
        write(*,*) 'CANON',size(c%data)
        do k = 1, size(c%data)
            write(*,*) c%coords(:,k), c%data(k)
        end do
        z = coo_nd_from_dense_flat(dims,flat,ierr)
        if (ierr /= 0) error stop 'from_dense_flat'
        vals = coo_nd_to_dense_flat(z,ierr)
        if (ierr /= 0) error stop 'from_flat roundtrip'
        call write_dense('FROM_FLAT',vals)
        roundtrip = 1
        if (ndim == 2) then
            old = coo_nd_to_coo(a,ierr)
            if (ierr /= 0) error stop 'to 2D COO'
            z = coo_nd_from_coo(old,ierr)
            if (ierr /= 0) error stop 'from 2D COO'
            if (any(z%coords /= a%coords) .or. any(z%data /= a%data)) roundtrip = 0
        end if
        write(*,*) 'ROUNDTRIP',roundtrip
        if (nnz > 0) then
            looked_up = coo_nd_get(a,coords(:,1),ierr)
            if (ierr /= 0) error stop 'get'
            write(*,*) 'GET',looked_up
        else
            write(*,*) 'GET',0.0_dp
        end if
        deallocate(dims,perm,target,coords,data)
    end do
contains
    subroutine write_dense(name,values)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: values(:)
        integer :: k
        write(*,*) name,size(values)
        do k = 1, size(values)
            write(*,'(ES26.17E3)') values(k)
        end do
    end subroutine write_dense
end program sparse_coo_nd_reference_driver
