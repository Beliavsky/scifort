! SPDX-License-Identifier: MIT
! Reference protocol for reusable banded LU in both arithmetic kinds.
program driver_banded_lu
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_banded, only : lu_factor_banded, lu_solve_banded
    implicit none
    integer :: cases, case_id, kind_id, lower, upper, n, nrhs, i, j, op, rep, status
    integer, allocatable :: ipiv(:)
    real(dp) :: re, im
    real(dp), allocatable :: ab_r(:,:), b_r(:,:,:), x_r(:,:), lu_r(:,:)
    complex(dp), allocatable :: ab_c(:,:), b_c(:,:,:), x_c(:,:), lu_c(:,:)
    character(len=1), parameter :: modes(3)=['N','T','C']

    read(*,*) cases
    do case_id=1,cases
        read(*,*) kind_id, lower, upper, n, nrhs
        allocate(ab_r(lower+upper+1,n),b_r(n,nrhs,2),x_r(n,nrhs))
        allocate(ab_c(lower+upper+1,n),b_c(n,nrhs,2),x_c(n,nrhs))
        allocate(lu_r(upper+2*lower+1,n),lu_c(upper+2*lower+1,n),ipiv(n))
        do j=1,n
            do i=1,lower+upper+1
                read(*,*) re,im
                ab_c(i,j)=cmplx(re,im,dp)
            end do
        end do
        do rep=1,2
            do j=1,nrhs
                do i=1,n
                    read(*,*) re,im
                    b_c(i,j,rep)=cmplx(re,im,dp)
                end do
            end do
        end do
        ab_r=real(ab_c,dp)
        b_r=real(b_c,dp)
        if (kind_id == 1) then
            call lu_factor_banded(ab_r,lower,upper,lu_r,ipiv,status)
        else
            call lu_factor_banded(ab_c,lower,upper,lu_c,ipiv,status)
        end if
        write(*,'(I0)') status
        if (status == linalg_status_success) then
            do op=1,3
                do rep=1,2
                    if (kind_id == 1) then
                        call lu_solve_banded(lu_r,lower,upper,ipiv,b_r(:,:,rep),x_r,status,modes(op))
                        x_c=cmplx(x_r,0.0_dp,dp)
                    else
                        call lu_solve_banded(lu_c,lower,upper,ipiv,b_c(:,:,rep),x_c,status,modes(op))
                    end if
                    write(*,'(I0)') status
                    if (status == linalg_status_success) then
                        do j=1,nrhs
                            do i=1,n
                                write(*,'(ES26.17E3,1X,ES26.17E3)') real(x_c(i,j),dp),aimag(x_c(i,j))
                            end do
                        end do
                    end if
                end do
            end do
        end if
        deallocate(ab_r,b_r,x_r,lu_r,ab_c,b_c,x_c,lu_c,ipiv)
    end do
end program driver_banded_lu
