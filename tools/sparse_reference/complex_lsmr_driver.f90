! SPDX-License-Identifier: BSD-3-Clause
! Batch comparison driver for complex CSR/CSC LSMR.
program complex_lsmr_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc
    use scifort_sparse_linalg_complex_lsmr, only : lsmr_complex
    implicit none
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    complex(dp), allocatable :: dense(:,:), b(:), x0(:), x(:)
    real(dp) :: re, im, damping, tol_a, tol_b, condlim, rnorm, arnorm, anorm, acond, xnorm
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
        if (status /= 0) error stop 'complex LSMR driver: invalid CSR'
        if (fmt == 2) then
            csc = csr_complex_to_csc(csr,status)
            if (status /= 0) error stop 'complex LSMR driver: invalid CSC'
            if (usex0 == 1) then
                call lsmr_complex(csc,b,x,stop_code,itn=itn,normr=rnorm,normar=arnorm, &
                    norma=anorm,conda=acond,normx=xnorm,damp=damping,atol=tol_a, &
                    btol=tol_b,conlim=condlim,maxiter=cap,x0=x0)
            else
                call lsmr_complex(csc,b,x,stop_code,itn=itn,normr=rnorm,normar=arnorm, &
                    norma=anorm,conda=acond,normx=xnorm,damp=damping,atol=tol_a, &
                    btol=tol_b,conlim=condlim,maxiter=cap)
            end if
        else
            if (usex0 == 1) then
                call lsmr_complex(csr,b,x,stop_code,itn=itn,normr=rnorm,normar=arnorm, &
                    norma=anorm,conda=acond,normx=xnorm,damp=damping,atol=tol_a, &
                    btol=tol_b,conlim=condlim,maxiter=cap,x0=x0)
            else
                call lsmr_complex(csr,b,x,stop_code,itn=itn,normr=rnorm,normar=arnorm, &
                    norma=anorm,conda=acond,normx=xnorm,damp=damping,atol=tol_a, &
                    btol=tol_b,conlim=condlim,maxiter=cap)
            end if
        end if
        write(*,'(i0,1x,i0,1x,i0,1x,5(es25.16,1x))') &
            stop_code,itn,size(x),rnorm,arnorm,anorm,acond,xnorm
        do i = 1,size(x)
            write(*,'(2(es25.16,1x))') real(x(i),dp),aimag(x(i))
        end do
        deallocate(dense,b,x0,x)
    end do
end program complex_lsmr_driver
