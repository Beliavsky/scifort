! SPDX-License-Identifier: MIT
! Streamed Fortran/NumPy/SciPy complex LSQR comparison driver.
program complex_lsqr_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc
    use scifort_sparse_linalg_complex_lsqr, only : lsqr_complex
    implicit none
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    complex(dp), allocatable :: dense(:,:), b(:), x0(:), x(:)
    real(dp) :: re, im, damping, tol_a, tol_b, condlim, r1, r2, arnorm
    integer :: nc, case_id, m, n, fmt, cap, usex0, i, j, status, stop_code, itn

    read(*,*) nc
    do case_id = 1,nc
        read(*,*) m,n,fmt,cap,usex0,damping,tol_a,tol_b,condlim
        allocate(dense(m,n),b(m),x0(n))
        do i = 1,m
            do j = 1,n
                read(*,*) re,im
                dense(i,j) = cmplx(re,im,dp)
            end do
        end do
        do i = 1,m
            read(*,*) re,im
            b(i) = cmplx(re,im,dp)
        end do
        do j = 1,n
            read(*,*) re,im
            x0(j) = cmplx(re,im,dp)
        end do
        csr = csr_complex_from_dense(dense,status)
        if (status /= 0) error stop 'complex LSQR driver: invalid CSR input'
        if (fmt == 2) then
            csc = csr_complex_to_csc(csr,status)
            if (status /= 0) error stop 'complex LSQR driver: CSC conversion'
            if (usex0 == 1) then
                call lsqr_complex(csc,b,x,stop_code,itn=itn,r1norm=r1,r2norm=r2, &
                    arnorm=arnorm,damp=damping,atol=tol_a,btol=tol_b,conlim=condlim, &
                    iter_lim=cap,x0=x0)
            else
                call lsqr_complex(csc,b,x,stop_code,itn=itn,r1norm=r1,r2norm=r2, &
                    arnorm=arnorm,damp=damping,atol=tol_a,btol=tol_b,conlim=condlim, &
                    iter_lim=cap)
            end if
        else
            if (usex0 == 1) then
                call lsqr_complex(csr,b,x,stop_code,itn=itn,r1norm=r1,r2norm=r2, &
                    arnorm=arnorm,damp=damping,atol=tol_a,btol=tol_b,conlim=condlim, &
                    iter_lim=cap,x0=x0)
            else
                call lsqr_complex(csr,b,x,stop_code,itn=itn,r1norm=r1,r2norm=r2, &
                    arnorm=arnorm,damp=damping,atol=tol_a,btol=tol_b,conlim=condlim, &
                    iter_lim=cap)
            end if
        end if
        write(*,'(i0,1x,i0,1x,i0,1x,3(es25.16,1x))') stop_code,itn,size(x),r1,r2,arnorm
        do i = 1,size(x)
            write(*,'(2(es25.16,1x))') real(x(i),dp),aimag(x(i))
        end do
        deallocate(dense,b,x0,x)
    end do
end program complex_lsqr_driver
