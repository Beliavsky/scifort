! SPDX-License-Identifier: MIT
! Text transport for independent SciPy comparisons. Each case contains a
! potentially non-Hermitian unused triangle to test selected-triangle behavior.
program driver_dense_cholesky
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_dense, only : cho_factor, cho_solve
    implicit none
    real(dp), allocatable :: a(:,:), b(:,:), b2(:,:), f(:,:), x(:,:), x2(:,:)
    complex(dp), allocatable :: ac(:,:), bc(:,:), bc2(:,:), fc(:,:), xc(:,:), xc2(:,:)
    real(dp) :: ar, ai
    integer :: ncase, icase, variant, n, nrhs, lower_flag, i, j, status
    logical :: lower

    read(*,*) ncase
    do icase=1,ncase
        read(*,*) variant, n, nrhs, lower_flag
        lower=lower_flag==1
        if (allocated(a)) deallocate(a,b,b2,f,x,x2,ac,bc,bc2,fc,xc,xc2)
        allocate(a(n,n), b(n,nrhs), b2(n,nrhs), f(n,n), x(n,nrhs), x2(n,nrhs))
        allocate(ac(n,n), bc(n,nrhs), bc2(n,nrhs), fc(n,n), xc(n,nrhs), xc2(n,nrhs))
        do j=1,n
            do i=1,n
                read(*,*) ar,ai
                a(i,j)=ar
                ac(i,j)=cmplx(ar,ai,dp)
            end do
        end do
        do j=1,nrhs
            do i=1,n
                read(*,*) ar,ai
                b(i,j)=ar
                bc(i,j)=cmplx(ar,ai,dp)
            end do
        end do
        do j=1,nrhs
            do i=1,n
                read(*,*) ar,ai
                b2(i,j)=ar
                bc2(i,j)=cmplx(ar,ai,dp)
            end do
        end do
        if (variant == 1) then
            call cho_factor(a,f,status,lower=lower)
            if (status == linalg_status_success) call cho_solve(f,b,x,status,lower=lower)
            if (status == linalg_status_success) call cho_solve(f,b2,x2,status,lower=lower)
        else
            call cho_factor(ac,fc,status,lower=lower)
            if (status == linalg_status_success) call cho_solve(fc,bc,xc,status,lower=lower)
            if (status == linalg_status_success) call cho_solve(fc,bc2,xc2,status,lower=lower)
        end if
        write(*,'(I0)') status
        if (status /= linalg_status_success) cycle
        do j=1,n
            do i=1,n
                if (variant == 1) then
                    write(*,'(ES28.18,1X,ES28.18)') f(i,j),0.0_dp
                else
                    write(*,'(ES28.18,1X,ES28.18)') real(fc(i,j),dp),aimag(fc(i,j))
                end if
            end do
        end do
        do j=1,nrhs
            do i=1,n
                if (variant == 1) then
                    write(*,'(ES28.18,1X,ES28.18)') x(i,j),0.0_dp
                else
                    write(*,'(ES28.18,1X,ES28.18)') real(xc(i,j),dp),aimag(xc(i,j))
                end if
            end do
        end do
        do j=1,nrhs
            do i=1,n
                if (variant == 1) then
                    write(*,'(ES28.18,1X,ES28.18)') x2(i,j),0.0_dp
                else
                    write(*,'(ES28.18,1X,ES28.18)') real(xc2(i,j),dp),aimag(xc2(i,j))
                end if
            end do
        end do
    end do
end program driver_dense_cholesky
