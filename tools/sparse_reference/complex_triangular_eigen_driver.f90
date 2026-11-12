! SPDX-License-Identifier: MIT
! Streaming numerical reference driver for complex triangular and Hermitian eigen APIs.
program complex_triangular_eigen_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc
    use scifort_sparse_linalg_complex_triangular, only : spsolve_triangular_complex
    use scifort_sparse_linalg_complex_eigen, only : eigsh_complex
    implicit none
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    complex(dp), allocatable :: dense(:,:), rhs(:,:), sol(:,:), vectors(:,:)
    real(dp), allocatable :: values(:), residuals(:)
    integer :: cases, case_no, n, nrhs, lower_flag, unit_flag, fmt, k, ncv, info, i, j, steps
    real(dp) :: re, im, tolerance
    character(len=1) :: task, trans
    character(len=2) :: which

    read(*,*) cases
    do case_no = 1,cases
        read(*,*) task
        if (task == 'T') then
            read(*,*) n,nrhs,lower_flag,unit_flag,trans,fmt
            allocate(dense(n,n),rhs(n,nrhs))
            do i = 1,n
                do j = 1,n
                    read(*,*) re,im
                    dense(i,j) = cmplx(re,im,dp)
                end do
            end do
            do i = 1,n
                do j = 1,nrhs
                    read(*,*) re,im
                    rhs(i,j) = cmplx(re,im,dp)
                end do
            end do
            csr = csr_complex_from_dense(dense,info)
            if (info /= 0) error stop 'input matrix rejected'
            csc = csr_complex_to_csc(csr,info)
            if (info /= 0) error stop 'CSC conversion failed'
            if (fmt == 1) then
                call spsolve_triangular_complex(csr,rhs,sol,info, &
                    lower=lower_flag==1,unit_diagonal=unit_flag==1,trans=trans)
            else
                call spsolve_triangular_complex(csc,rhs,sol,info, &
                    lower=lower_flag==1,unit_diagonal=unit_flag==1,trans=trans)
            end if
            print '(a,1x,i0,1x,i0,1x,i0)', 'T',info,size(sol,1),size(sol,2)
            if (info == 0) then
                do i = 1,n
                    do j = 1,nrhs
                        print '(2es26.17e3)', real(sol(i,j),dp),aimag(sol(i,j))
                    end do
                end do
            end if
            deallocate(dense,rhs,sol)
        else if (task == 'E') then
            read(*,*) n,k,which,ncv,tolerance,fmt
            allocate(dense(n,n))
            do i = 1,n
                do j = 1,n
                    read(*,*) re,im
                    dense(i,j) = cmplx(re,im,dp)
                end do
            end do
            csr = csr_complex_from_dense(dense,info)
            if (info /= 0) error stop 'bad eigen matrix'
            csc = csr_complex_to_csc(csr,info)
            if (info /= 0) error stop 'CSC conversion failed'
            if (fmt == 1) then
                call eigsh_complex(csr,k,values,vectors,info,which=which,tol=tolerance, &
                    ncv=ncv,residual_norms=residuals,iterations=steps)
            else
                call eigsh_complex(csc,k,values,vectors,info,which=which,tol=tolerance, &
                    ncv=ncv,residual_norms=residuals,iterations=steps)
            end if
            print '(a,1x,i0,1x,i0,1x,i0,1x,i0)', 'E',info,size(values),size(vectors,1),steps
            if (info == 0 .or. info == 1) then
                do j = 1,k
                    print '(2es26.17e3)', values(j),residuals(j)
                    do i = 1,n
                        print '(2es26.17e3)', real(vectors(i,j),dp),aimag(vectors(i,j))
                    end do
                end do
            end if
            deallocate(dense,values,vectors,residuals)
        else
            error stop 'unknown task'
        end if
    end do
end program complex_triangular_eigen_driver
