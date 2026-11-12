! SPDX-License-Identifier: MIT
program sparse_formats_norm_reference_driver
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, dia_matrix, bsr_matrix, &
        csr_from_coo, csr_to_dia, csr_to_bsr, dia_to_dense, bsr_to_dense, &
        dia_matvec, bsr_matvec
    use scifort_sparse_linalg, only : norm, norm_axis
    implicit none
    type(csr_matrix) :: a
    type(dia_matrix) :: d
    type(bsr_matrix) :: b
    real(dp), allocatable :: vals(:), x(:), dmat(:,:), bmat(:,:), yd(:), yb(:), vcol(:), vrow(:)
    integer, allocatable :: rows(:), cols(:)
    real(dp) :: norms(6)
    integer :: m,n,nnz,br,bc,i,k,ncases,stat,stat2,stat3,stat4

    read(*,*) ncases
    do i = 1,ncases
        read(*,*) m,n,nnz,br,bc
        allocate(rows(nnz),cols(nnz),vals(nnz),x(n))
        do k = 1,nnz
            read(*,*) rows(k),cols(k),vals(k)
        end do
        do k = 1,n
            x(k) = real(k,dp)/real(max(1,n),dp)
        end do
        a = csr_from_coo(m,n,rows,cols,vals,stat)
        if (stat /= 0) error stop 'CSR construction'
        d = csr_to_dia(a,stat)
        b = csr_to_bsr(a,br,bc,stat2)
        if (stat /= 0 .or. stat2 /= 0) error stop 'format conversion'
        yd = dia_matvec(d,x,stat)
        yb = bsr_matvec(b,x,stat2)
        dmat = dia_to_dense(d,stat3)
        bmat = bsr_to_dense(b,stat4)
        if (stat /= 0 .or. stat2 /= 0 .or. stat3 /= 0 .or. stat4 /= 0) error stop 'format values'
        norms(1) = norm(a)
        norms(2) = norm(a,1)
        norms(3) = norm(a,-1)
        norms(4) = norm(a,ieee_value(0.0_dp,ieee_positive_inf))
        norms(5) = norm(a,ieee_value(0.0_dp,ieee_negative_inf))
        norms(6) = -1.0_dp
        if (mod(i,5) == 0) then
            norms(6) = norm(a,2,stat)
            if (stat /= 0) norms(6) = -1.0_dp
        end if
        vcol = norm_axis(a,1)
        vrow = norm_axis(a,2,ord=1.0_dp)
        write(*,'(a,1x,*(es24.16,1x))') 'N',norms
        write(*,'(a,1x,*(es24.16,1x))') 'C',vcol
        write(*,'(a,1x,*(es24.16,1x))') 'R',vrow
        write(*,'(a,1x,*(es24.16,1x))') 'D',reshape(dmat,[m*n])
        write(*,'(a,1x,*(es24.16,1x))') 'B',reshape(bmat,[m*n])
        write(*,'(a,1x,*(es24.16,1x))') 'Y',yd
        write(*,'(a,1x,*(es24.16,1x))') 'Z',yb
        deallocate(rows,cols,vals,x)
    end do
end program sparse_formats_norm_reference_driver
