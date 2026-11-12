! SPDX-License-Identifier: MIT
program sparse_arnoldi_reference_driver
    use scifort_kinds, only: dp
    use scifort_sparse, only: csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only: eigs_arnoldi
    implicit none
    integer :: n,k,p,cap,info,cycles,locks,mvs,i
    character(len=2) :: mode
    real(dp) :: tol
    real(dp), allocatable :: dense(:,:),errors(:)
    complex(dp), allocatable :: values(:),vectors(:,:)
    type(csr_matrix) :: a

    read(*,*) n,k,p,cap,mode,tol
    allocate(dense(n,n))
    read(*,*) dense
    a=csr_from_dense(dense)
    call eigs_arnoldi(a,k,values,vectors,info,which=mode,tol=tol, &
        ncv=p,maxiter=cap,iterations=cycles,locked_count=locks,matvecs=mvs,residual_norms=errors)
    write(*,'(4(i0,1x))') info,cycles,locks,mvs
    do i=1,size(values)
        write(*,'(3(es26.17,1x))') real(values(i),dp),aimag(values(i)),errors(i)
    end do
end program sparse_arnoldi_reference_driver
