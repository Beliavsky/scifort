! SPDX-License-Identifier: MIT
! Text I/O driver for SciPy comparisons of the dense solvers.
program driver_dense
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_no_convergence
    use scifort_linalg_dense, only : solve_triangular, lu_factor, lu_solve
    implicit none
    integer :: ncases, icase, variant, n, nrhs, lower_i, unit_i, i, j, status
    integer, allocatable :: ipiv(:)
    character(len=1) :: op
    real(dp) :: r, im
    real(dp), allocatable :: a_r(:,:), b_r(:,:), x_r(:,:), lu_r(:,:), y_r(:,:)
    complex(dp), allocatable :: a_c(:,:), b_c(:,:), x_c(:,:), lu_c(:,:), y_c(:,:)

    read(*,*) ncases
    do icase=1,ncases
        read(*,*) variant, n, nrhs, op, lower_i, unit_i
        allocate(a_r(n,n),b_r(n,nrhs),x_r(n,nrhs),lu_r(n,n),y_r(n,nrhs))
        allocate(a_c(n,n),b_c(n,nrhs),x_c(n,nrhs),lu_c(n,n),y_c(n,nrhs))
        allocate(ipiv(n))
        do j=1,n
            do i=1,n
                read(*,*) r, im
                a_r(i,j)=r
                a_c(i,j)=cmplx(r,im,dp)
            end do
        end do
        do j=1,nrhs
            do i=1,n
                read(*,*) r, im
                b_r(i,j)=r
                b_c(i,j)=cmplx(r,im,dp)
            end do
        end do
        select case (variant)
        case(1)
            call solve_triangular(a_r,b_r,x_r,status,lower=lower_i==1, &
                trans=op,unit_diagonal=unit_i==1)
            x_c=cmplx(x_r,0.0_dp,dp)
        case(2)
            call solve_triangular(a_c,b_c,x_c,status,lower=lower_i==1, &
                trans=op,unit_diagonal=unit_i==1)
        case(3)
            call lu_factor(a_r,lu_r,ipiv,status)
            if (status == linalg_status_success) then
                call lu_solve(lu_r,ipiv,b_r,x_r,status,trans=op)
                if (status == linalg_status_success) then
                    call lu_solve(lu_r,ipiv,2.0_dp*b_r,y_r,status,trans=op)
                    if (status == linalg_status_success) then
                        if (maxval(abs(y_r-2.0_dp*x_r)) > &
                            2.0e-10_dp*max(1.0_dp,maxval(abs(x_r)))) &
                            status=linalg_status_no_convergence
                    end if
                end if
            end if
            x_c=cmplx(x_r,0.0_dp,dp)
        case(4)
            call lu_factor(a_c,lu_c,ipiv,status)
            if (status == linalg_status_success) then
                call lu_solve(lu_c,ipiv,b_c,x_c,status,trans=op)
                if (status == linalg_status_success) then
                    call lu_solve(lu_c,ipiv,2.0_dp*b_c,y_c,status,trans=op)
                    if (status == linalg_status_success) then
                        if (maxval(abs(y_c-2.0_dp*x_c)) > &
                            2.0e-10_dp*max(1.0_dp,maxval(abs(x_c)))) &
                            status=linalg_status_no_convergence
                    end if
                end if
            end if
        case default
            status=linalg_status_no_convergence
        end select
        write(*,'(I0)') status
        if (status == linalg_status_success) then
            do j=1,nrhs
                do i=1,n
                    write(*,'(ES26.17E3,1X,ES26.17E3)') real(x_c(i,j),dp),aimag(x_c(i,j))
                end do
            end do
        end if
        deallocate(a_r,b_r,x_r,lu_r,y_r,a_c,b_c,x_c,lu_c,y_c,ipiv)
    end do
end program driver_dense
