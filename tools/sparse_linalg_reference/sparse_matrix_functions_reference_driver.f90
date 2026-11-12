! SPDX-License-Identifier: MIT
! Reference transport for matrix_power and sparse inverse randomized tests.
program sparse_matrix_functions_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csc_matrix, csr_from_coo, &
        csr_to_dense, csr_to_csc, csc_to_csr
    use scifort_sparse_linalg, only : matrix_power, inv
    implicit none
    type(csr_matrix) :: a, p, q
    type(csc_matrix) :: ac, pc, qc
    integer, allocatable :: rows(:), cols(:)
    real(dp), allocatable :: data(:), dense(:,:)
    integer :: ncases, case_num, n, nnz, power, i, ierr, sp, sc, si, sic

    read(*,*) ncases
    do case_num = 1, ncases
        read(*,*) n,nnz,power
        allocate(rows(nnz),cols(nnz),data(nnz))
        do i = 1,nnz
            read(*,*) rows(i),cols(i),data(i)
        end do
        a = csr_from_coo(n,n,rows,cols,data,ierr)
        if (ierr /= 0) error stop 'invalid reference input'
        ac = csr_to_csc(a,ierr)
        if (ierr /= 0) error stop 'invalid CSC conversion'

        p = matrix_power(a,power,sp)
        pc = matrix_power(ac,power,sc)
        if (sc == 0) then
            q = csc_to_csr(pc,ierr)
            if (ierr /= 0) error stop 'power CSC conversion'
        end if
        write(*,'(4(i0,1x))') case_num,sp,sc,power
        if (sp == 0) then
            dense = csr_to_dense(p)
        else
            allocate(dense(n,n))
            dense = 0.0_dp
        end if
        call write_dense(dense)
        if (sc == 0) then
            dense = csr_to_dense(q)
        else
            dense = 0.0_dp
        end if
        call write_dense(dense)

        p = inv(a,si,block_size=3)
        qc = inv(ac,sic,block_size=5)
        if (sic == 0) then
            q = csc_to_csr(qc,ierr)
            if (ierr /= 0) error stop 'inverse CSC conversion'
        end if
        write(*,'(3(i0,1x))') case_num,si,sic
        if (si == 0) then
            dense = csr_to_dense(p)
        else
            dense = 0.0_dp
        end if
        call write_dense(dense)
        if (sic == 0) then
            dense = csr_to_dense(q)
        else
            dense = 0.0_dp
        end if
        call write_dense(dense)
        deallocate(rows,cols,data)
    end do
contains
    subroutine write_dense(a)
        real(dp), intent(in) :: a(:,:) !! Matrix in Fortran column-major layout.
        write(*,'(*(es26.17e3,1x))') reshape(a,[size(a)])
    end subroutine write_dense
end program sparse_matrix_functions_reference_driver
