! SPDX-License-Identifier: MIT
! Streamed real/complex indefinite linear solves for SciPy comparisons.
program driver_indefinite_solve
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_ldl, only : ldl
    use scifort_linalg_indefinite_solve, only : ldl_solve, solve_indefinite
    implicit none

    real(dp), allocatable :: ar(:, :), br(:, :), lr(:, :), dr(:, :)
    real(dp), allocatable :: xr(:, :), xrf(:, :)
    complex(dp), allocatable :: a(:, :), b(:, :), lu(:, :), d(:, :)
    complex(dp), allocatable :: x(:, :), xf(:, :)
    integer, allocatable :: perm(:)
    integer :: n, nrhs, family, orient, sd, sf, i, j, ios
    real(dp) :: re, im

    do
        read(*,*,iostat=ios) n, nrhs, family, orient
        if (ios /= 0) exit
        allocate(a(n,n),b(n,nrhs),lu(n,n),d(n,n),x(n,nrhs),xf(n,nrhs), &
            ar(n,n),br(n,nrhs),lr(n,n),dr(n,n),xr(n,nrhs),xrf(n,nrhs),perm(n))
        do j = 1, n
            do i = 1, n
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
                ar(i,j) = re
            end do
        end do
        do j = 1, nrhs
            do i = 1, n
                read(*,*) re, im
                b(i,j) = cmplx(re,im,dp)
                br(i,j) = re
            end do
        end do

        select case (family)
        case (1)
            call solve_indefinite(ar,br,xr,sd,lower=(orient == 1))
            call ldl(ar,lr,dr,perm,sf,lower=(orient == 1))
            xrf = 0.0_dp
            if (sf == linalg_status_success) &
                call ldl_solve(lr,dr,perm,br,xrf,sf,lower=(orient == 1))
            x = cmplx(xr,kind=dp)
            xf = cmplx(xrf,kind=dp)
        case (2)
            call solve_indefinite(a,b,x,sd,lower=(orient == 1))
            call ldl(a,lu,d,perm,sf,lower=(orient == 1))
            xf = cmplx(0.0_dp,0.0_dp,dp)
            if (sf == linalg_status_success) &
                call ldl_solve(lu,d,perm,b,xf,sf,lower=(orient == 1))
        case (3)
            call solve_indefinite(a,b,x,sd,lower=(orient == 1),hermitian=.false.)
            call ldl(a,lu,d,perm,sf,lower=(orient == 1),hermitian=.false.)
            xf = cmplx(0.0_dp,0.0_dp,dp)
            if (sf == linalg_status_success) &
                call ldl_solve(lu,d,perm,b,xf,sf,lower=(orient == 1),hermitian=.false.)
        case default
            sd = -1
            sf = -1
        end select
        print '(2(i0,1x))', sd,sf
        do j = 1, nrhs
            do i = 1, n
                print '(4(es26.17e3,1x))', real(x(i,j),dp), aimag(x(i,j)), &
                    real(xf(i,j),dp),aimag(xf(i,j))
            end do
        end do
        deallocate(a,b,lu,d,x,xf,ar,br,lr,dr,xr,xrf,perm)
    end do
end program driver_indefinite_solve
