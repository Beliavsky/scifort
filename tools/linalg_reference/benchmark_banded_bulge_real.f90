! SPDX-License-Identifier: MIT
! Deterministic microbenchmark, not a substitute for correctness testing.
program bench_banded_methods
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_banded_eigen, only : eigvals_banded
    use scifort_linalg_banded_tridiagonal, only : band_to_tridiagonal
    implicit none
    integer :: n, kd, i, j, r, st, reps, trial
    real(dp), allocatable :: ab(:,:), d(:), e(:), w(:)
    complex(dp), allocatable :: cab(:,:)
    real(dp) :: start, finish, dt_real, dt_complex, dt_bulge, dt_dense
    integer, parameter :: sizes(*) = [120, 320, 500]
    integer, parameter :: widths(*) = [1, 2, 4, 8, 16]

    print '(a)', 'n kd tridiag_real tridiag_complex eigvals_bulge eigvals_dense'
    do i=1,size(sizes)
        n=sizes(i)
        do j=1,size(widths)
            kd=widths(j)
            allocate(ab(kd+1,n),cab(kd+1,n),d(n),e(n-1),w(n))
            ab=0.0_dp
            do r=1,n
                ab(1,r)=2.0_dp+sin(0.13_dp*real(r,dp))
                do trial=1,min(kd,n-r)
                    ab(trial+1,r)=0.06_dp*cos(real(5*r+7*trial,dp))
                end do
            end do
            cab=cmplx(ab,0.0_dp,dp)
            reps=2
            call cpu_time(start)
            do trial=1,reps
                call band_to_tridiagonal(ab,d,e,st,lower=.true.)
                if (st /= linalg_status_success) error stop 'real reduction failed'
            end do
            call cpu_time(finish)
            dt_real=(finish-start)/reps
            call cpu_time(start)
            do trial=1,reps
                call band_to_tridiagonal(cab,d,e,st,lower=.true.)
                if (st /= linalg_status_success) error stop 'complex reduction failed'
            end do
            call cpu_time(finish)
            dt_complex=(finish-start)/reps
            call cpu_time(start)
            call eigvals_banded(ab,w,st,lower=.true.,method='bulge')
            if (st /= linalg_status_success) error stop 'bulge solve failed'
            call cpu_time(finish)
            dt_bulge=finish-start
            call cpu_time(start)
            call eigvals_banded(ab,w,st,lower=.true.,method='dense')
            if (st /= linalg_status_success) error stop 'dense solve failed'
            call cpu_time(finish)
            dt_dense=finish-start
            print '(i4,1x,i3,4(1x,es12.4))',n,kd,dt_real,dt_complex,dt_bulge,dt_dense
            deallocate(ab,cab,d,e,w)
        end do
    end do
end program bench_banded_methods
