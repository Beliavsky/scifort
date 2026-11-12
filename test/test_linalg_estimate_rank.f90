! SPDX-License-Identifier: MIT
program test_linalg_estimate_rank
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_estimate_rank, only : estimate_rank, estimate_rank_operator, &
        estimate_rank_operator_complex
    implicit none
    real(dp), allocatable :: a(:, :), lr(:, :), rr(:, :)
    complex(dp), allocatable :: ac(:, :)
    integer :: m, n, k, rank, ro, rc, rco, rank2, status, i, j, variant, checks
    real(dp) :: t, scale, nan
    real(dp), parameter :: pi = acos(-1.0_dp)

    checks = 0
    do m = 0, 15
        do n = 0, 17
            allocate(a(m,n),ac(m,n),lr(m,min(m,n)),rr(min(m,n),n))
            do variant = 1, 4
                select case (variant)
                case (1)
                    k = 0
                case (2)
                    k = min(1,min(m,n))
                case (3)
                    k = min(m,n)/2
                case (4)
                    k = min(m,n)
                end select
                lr = 0.0_dp
                rr = 0.0_dp
                do j = 1, k
                    do i = 1, m
                        lr(i,j) = sqrt(2.0_dp/real(m+1,dp)) * &
                            sin(pi*real(i*j,dp)/real(m+1,dp))
                    end do
                    do i = 1, n
                        rr(j,i) = sqrt(2.0_dp/real(n+1,dp)) * &
                            sin(pi*real(j*i,dp)/real(n+1,dp))
                    end do
                end do
                a = matmul(lr,rr)
                do j = 1, n
                    do i = 1, m
                        t = real(i+j,dp)*0.31_dp
                        ac(i,j) = cmplx(cos(t),sin(t),dp)*a(i,j)
                    end do
                end do
                select case (mod(m+n,3))
                case (0)
                    scale = 1.0_dp
                case (1)
                    scale = 1.0e-90_dp
                case (2)
                    scale = 1.0e90_dp
                end select
                a = a*scale
                ac = ac*scale
                call estimate_rank(a,1.0e-9_dp,rank,status,seed=123456_int64)
                if (status /= linalg_status_success .or. rank /= k) then
                    print *, 'real dense', m,n,k,rank,status,scale
                    error stop 'real dense rank'
                end if
                call estimate_rank(ac,1.0e-9_dp,rc,status,seed=123456_int64)
                if (status /= linalg_status_success .or. rc /= k) then
                    print *, 'complex dense', m,n,k,rc,status,scale
                    error stop 'complex dense rank'
                end if
                call estimate_rank_operator(m,n,mv_real,1.0e-9_dp,ro,status,seed=123456_int64)
                if (status /= linalg_status_success .or. ro /= k) then
                    print *, 'real callback',m,n,k,ro,status,scale
                    error stop 'real operator rank'
                end if
                call estimate_rank_operator_complex(m,n,mv_complex,1.0e-9_dp,rco,status, &
                    seed=123456_int64)
                if (status /= linalg_status_success .or. rco /= k) then
                    print *, 'complex callback',m,n,k,rco,status,scale
                    error stop 'complex operator rank'
                end if
                call estimate_rank(a,1.0e-9_dp,rank2,status,seed=123456_int64)
                if (status /= linalg_status_success .or. rank2 /= rank) &
                    error stop 'seed reproducibility'
                checks = checks+5
            end do
            deallocate(a,ac,lr,rr)
        end do
    end do

    allocate(a(6,6),ac(6,6))
    a = 0.0_dp
    do i = 1, 5
        a(i,i) = 1.0_dp
    end do
    a(6,6) = 1.0e-11_dp
    ac = cmplx(a,0.0_dp,dp)
    call estimate_rank(a,1.0e-7_dp,rank,status,seed=71_int64)
    if (rank /= 5 .or. status /= linalg_status_success) error stop 'real thresholded rank'
    call estimate_rank(ac,1.0e-7_dp,rank,status,seed=71_int64)
    if (rank /= 5 .or. status /= linalg_status_success) error stop 'complex thresholded rank'
    call estimate_rank(a,1.0e-13_dp,rank,status,seed=71_int64)
    if (rank /= 6 .or. status /= linalg_status_success) error stop 'small singular retained'
    checks = checks+3

    call estimate_rank(a,-1.0_dp,rank,status)
    if (status /= linalg_status_invalid_shape .or. rank /= 0) error stop 'negative eps'
    call estimate_rank(a,1.0_dp,rank,status)
    if (status /= linalg_status_invalid_shape) error stop 'eps >= 1'
    call estimate_rank(ac,0.1_dp,rank,status,oversampling=0)
    if (status /= linalg_status_invalid_shape) error stop 'invalid oversampling'
    nan = ieee_value(1.0_dp,ieee_quiet_nan)
    a(3,2) = nan
    call estimate_rank(a,0.1_dp,rank,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite real'
    ac(2,3) = cmplx(0.0_dp,nan,dp)
    call estimate_rank(ac,0.1_dp,rank,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite complex'
    call estimate_rank_operator(6,6,mv_fail,1.0e-8_dp,rank,status,seed=1_int64)
    if (status /= linalg_status_no_convergence .or. rank /= 0) error stop 'failed callback'
    call estimate_rank_operator(6,6,mv_real,2.0_dp,rank,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid operator tolerance'
    call estimate_rank_operator_complex(-1,2,mv_complex,0.1_dp,rank,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid operator dimensions'
    checks = checks+8
    print '(a,i0)', 'PASS randomized-rank native checks: ',checks
contains
    subroutine mv_real(x,y,status)
        real(dp), intent(in) :: x(:) !! input vector
        real(dp), intent(out) :: y(:) !! output vector
        integer, intent(out) :: status !! success
        y = matmul(a,x)
        status = 0
    end subroutine mv_real
    subroutine mv_complex(x,y,status)
        complex(dp), intent(in) :: x(:) !! input vector
        complex(dp), intent(out) :: y(:) !! output vector
        integer, intent(out) :: status !! success
        y = matmul(ac,x)
        status = 0
    end subroutine mv_complex
    subroutine mv_fail(x,y,status)
        real(dp), intent(in) :: x(:) !! input vector
        real(dp), intent(out) :: y(:) !! output vector
        integer, intent(out) :: status !! failure
        y = 0.0_dp
        status = 1
    end subroutine mv_fail
end program test_linalg_estimate_rank
