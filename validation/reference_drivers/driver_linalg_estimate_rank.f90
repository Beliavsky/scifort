! SPDX-License-Identifier: MIT
! Numerical comparison driver: matrix entries supplied by Python, in row-major order.
program driver_linalg_estimate_rank
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_estimate_rank, only : estimate_rank, estimate_rank_operator, &
        estimate_rank_operator_complex
    implicit none
    real(dp), allocatable :: a(:, :)
    complex(dp), allocatable :: z(:, :)
    real(dp) :: eps, re, im
    integer :: cases, icase, kind, m, n, i, j, rd, ro, sd, so
    integer(int64) :: seed
    read (*,*) cases
    do icase = 1, cases
        read (*,*) kind,m,n,eps,seed
        allocate(a(m,n),z(m,n))
        do i = 1, m
            do j = 1, n
                read (*,*) re,im
                a(i,j) = re
                z(i,j) = cmplx(re,im,dp)
            end do
        end do
        if (kind == 0) then
            call estimate_rank(a,eps,rd,sd,seed=seed)
            call estimate_rank_operator(m,n,mv_real,eps,ro,so,seed=seed)
        else
            call estimate_rank(z,eps,rd,sd,seed=seed)
            call estimate_rank_operator_complex(m,n,mv_complex,eps,ro,so,seed=seed)
        end if
        print '(a,1x,i0,4(1x,i0))', 'CASE',icase,rd,ro,sd,so
        deallocate(a,z)
    end do
contains
    subroutine mv_real(x,y,status)
        real(dp), intent(in) :: x(:) !! n-vector
        real(dp), intent(out) :: y(:) !! m-vector
        integer, intent(out) :: status !! success
        y = matmul(a,x)
        status = 0
    end subroutine mv_real
    subroutine mv_complex(x,y,status)
        complex(dp), intent(in) :: x(:) !! n-vector
        complex(dp), intent(out) :: y(:) !! m-vector
        integer, intent(out) :: status !! success
        y = matmul(z,x)
        status = 0
    end subroutine mv_complex
end program driver_linalg_estimate_rank
