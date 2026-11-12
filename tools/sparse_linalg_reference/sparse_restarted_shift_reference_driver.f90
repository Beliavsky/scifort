! SPDX-License-Identifier: MIT
! Input: case count, then (mode,n,k,ncv,cap,which,tol,sigma_real,sigma_imag), then dense A.
! Output: original-operator eigenpairs/residuals for independent Python comparisons.
program sparse_restarted_shift_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : eigsh_lanczos, eigs_arnoldi
    implicit none
    integer :: count, ci, n, k, p, cap, info, cycles, locks, ops, j
    character(len=1) :: kind
    character(len=2) :: which
    real(dp) :: tol, sr, si
    real(dp), allocatable :: dense(:,:), rval(:), rvec(:,:), errors(:)
    complex(dp), allocatable :: cval(:), cvec(:,:)
    type(csr_matrix) :: a
    read(*,*) count
    do ci=1,count
        read(*,*) kind,n,k,p,cap,which,tol,sr,si
        allocate(dense(n,n))
        read(*,*) dense
        a=csr_from_dense(dense)
        if (kind == 'S') then
            call eigsh_lanczos(a,k,rval,rvec,info,which=which,tol=tol, &
                maxiter=cap,ncv=p,sigma=sr,inner_rtol=1.0e-12_dp, &
                inner_maxiter=max(100,12*n),iterations=cycles, &
                locked_count=locks,matvecs=ops,residual_norms=errors)
        else
            call eigs_arnoldi(a,k,cval,cvec,info,which=which,tol=tol, &
                maxiter=cap,ncv=p,sigma=cmplx(sr,si,dp), &
                inner_rtol=1.0e-11_dp,inner_maxiter=max(100,12*n), &
                inner_restart=min(n,40),iterations=cycles, &
                locked_count=locks,matvecs=ops,residual_norms=errors)
        end if
        write(*,'(a,1x,i0)') 'CASE',ci
        write(*,'(4(i0,1x))') info,cycles,locks,ops
        if (info >= 0) then
            do j=1,k
                if (kind == 'S') then
                    write(*,'(*(es26.17,1x))') rval(j),0.0_dp,errors(j),rvec(:,j)
                else
                    write(*,'(*(es26.17,1x))') real(cval(j),dp),aimag(cval(j)), &
                        errors(j),real(cvec(:,j),dp),aimag(cvec(:,j))
                end if
            end do
        end if
        deallocate(dense)
    end do
end program sparse_restarted_shift_reference_driver
