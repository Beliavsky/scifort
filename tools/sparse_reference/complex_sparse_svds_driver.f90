! SPDX-License-Identifier: MIT
! Streaming reference driver for complex sparse SVD comparisons.
program complex_sparse_svds_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_coo, csr_complex_to_csc
    use scifort_sparse_linalg_complex_svds, only : svds_complex
    implicit none
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    complex(dp), allocatable :: u(:,:),vh(:,:),values(:),v0(:)
    integer, allocatable :: rows(:),cols(:)
    real(dp), allocatable :: singular(:),residual(:)
    integer :: cases, case_id, m, n, nz, k, width, cap, fmt, use_v0, &
        j, status, iterations
    real(dp) :: re,im,tolerance
    character(len=2) :: which

    read(*,*) cases
    do case_id = 1,cases
        read(*,*) m,n,nz,k,width,cap,which,tolerance,fmt,use_v0
        allocate(rows(nz),cols(nz),values(nz))
        do j = 1,nz
            read(*,*) rows(j),cols(j),re,im
            values(j) = cmplx(re,im,dp)
        end do
        if (use_v0 /= 0) then
            allocate(v0(min(m,n)))
            do j = 1,size(v0)
                read(*,*) re,im
                v0(j) = cmplx(re,im,dp)
            end do
        end if
        csr = csr_complex_from_coo(m,n,rows,cols,values,status)
        if (status /= 0) error stop 'failed to construct CSR matrix'
        csc = csr_complex_to_csc(csr,status)
        if (status /= 0) error stop 'failed to construct CSC matrix'
        if (fmt == 1) then
            if (use_v0 /= 0) then
                call svds_complex(csr,k,u,singular,vh,status,which=which,tol=tolerance, &
                    ncv=width,maxiter=cap,v0=v0,iterations=iterations,residual_norms=residual)
            else
                call svds_complex(csr,k,u,singular,vh,status,which=which,tol=tolerance, &
                    ncv=width,maxiter=cap,iterations=iterations,residual_norms=residual)
            end if
        else
            if (use_v0 /= 0) then
                call svds_complex(csc,k,u,singular,vh,status,which=which,tol=tolerance, &
                    ncv=width,maxiter=cap,v0=v0,iterations=iterations,residual_norms=residual)
            else
                call svds_complex(csc,k,u,singular,vh,status,which=which,tol=tolerance, &
                    ncv=width,maxiter=cap,iterations=iterations,residual_norms=residual)
            end if
        end if
        print '(a,1x,i0)', 'CASE', case_id
        print '(3(i0,1x))', status,iterations,size(singular)
        if (size(singular) > 0) then
            do j = 1,k
                print '(2es26.17e3)', singular(j),residual(j)
                print '( *(2es26.17e3,1x) )', u(:,j)
                print '( *(2es26.17e3,1x) )', vh(j,:)
            end do
        end if
        deallocate(rows,cols,values)
        if (allocated(v0)) deallocate(v0)
        deallocate(u,singular,vh,residual)
    end do
end program complex_sparse_svds_driver
