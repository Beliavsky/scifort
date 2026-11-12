! SPDX-License-Identifier: MIT
! Stdin/stdout driver for independent SciPy cross-checks.
program driver_circulant
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_circulant, only : solve_circulant
    implicit none
    integer :: n, nrhs, input_kind, mode, code, ios, i, j
    real(dp) :: tol
    real(dp), allocatable :: cr(:), ci(:), br(:,:), bi(:,:), xr(:,:)
    complex(dp), allocatable :: c(:), b(:,:), x(:,:)
    character(:), allocatable :: policy

    do
        read(*,*,iostat=ios) n,nrhs,input_kind,mode,tol
        if (ios /= 0) exit
        if (n < 0) exit
        allocate(cr(n),ci(n),br(n,nrhs),bi(n,nrhs),xr(n,nrhs),c(n),b(n,nrhs),x(n,nrhs))
        if (n > 0) then
            read(*,*) cr
            read(*,*) ci
        end if
        do j=1,nrhs
            if (n > 0) then
                read(*,*) br(:,j)
                read(*,*) bi(:,j)
            end if
        end do
        c = cmplx(cr,ci,dp)
        b = cmplx(br,bi,dp)
        policy = 'raise'
        if (mode == 1) policy = 'lstsq'
        if (tol < 0.0_dp) then
            select case (input_kind)
            case (0)
                call solve_circulant(cr,br,xr,code,singular=policy)
                x = cmplx(xr,kind=dp)
            case (1)
                call solve_circulant(c,b,x,code,singular=policy)
            case (2)
                call solve_circulant(cr,b,x,code,singular=policy)
            case (3)
                call solve_circulant(c,br,x,code,singular=policy)
            end select
        else
            select case (input_kind)
            case (0)
                call solve_circulant(cr,br,xr,code,singular=policy,tol=tol)
                x = cmplx(xr,kind=dp)
            case (1)
                call solve_circulant(c,b,x,code,singular=policy,tol=tol)
            case (2)
                call solve_circulant(cr,b,x,code,singular=policy,tol=tol)
            case (3)
                call solve_circulant(c,br,x,code,singular=policy,tol=tol)
            end select
        end if
        write(*,'(i0)') code
        if (code == linalg_status_success) then
            do j = 1,nrhs
                do i = 1,n
                    write(*,'(es26.17e3,1x,es26.17e3)') real(x(i,j),dp),aimag(x(i,j))
                end do
            end do
        end if
        deallocate(cr,ci,br,bi,xr,c,b,x)
    end do
end program driver_circulant
