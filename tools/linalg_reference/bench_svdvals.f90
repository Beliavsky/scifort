! SPDX-License-Identifier: MIT
! Reproducible performance comparison: economy svd vs singular values only.
program bench_svdvals
    use scifort_kinds, only : dp
    use scifort_linalg_svd, only : svd, svdvals
    implicit none
    integer, parameter :: m = 120, n = 80, repetitions = 4
    real(dp) :: ar(m,n), br(m,n), ur(m,n), vhr(n,n), s(n), ref(n)
    complex(dp) :: ac(m,n), uc(m,n), vhc(n,n)
    real(dp) :: start, finish, t_svd, t_vals, err
    integer :: status, j, seed_len
    integer, allocatable :: seed_values(:)

    call random_seed(size=seed_len)
    allocate(seed_values(seed_len))
    seed_values = 1701
    call random_seed(put=seed_values)
    call random_number(ar)
    call random_number(br)
    ar = ar-0.5_dp
    ac = cmplx(ar,br-0.5_dp,dp)
    call svd(ar,ur,ref,vhr,status)
    if (status /= 0) error stop 'real svd warmup'
    call svdvals(ar,s,status)
    if (status /= 0) error stop 'real svdvals warmup'
    err = maxval(abs(ref-s))/max(1.0_dp,ref(1))
    if (err > 1.0e-11_dp) error stop 'real spectrum parity'
    call cpu_time(start)
    do j = 1, repetitions
        call svd(ar,ur,ref,vhr,status)
        if (status /= 0) error stop 'real svd timed'
    end do
    call cpu_time(finish)
    t_svd = finish-start
    call cpu_time(start)
    do j = 1, repetitions
        call svdvals(ar,s,status)
        if (status /= 0) error stop 'real svdvals timed'
    end do
    call cpu_time(finish)
    t_vals = finish-start
    print '(A,ES12.4,A,ES12.4,A,F8.3)', 'real: svd=', t_svd, &
        ' svdvals=',t_vals,' speedup=',t_svd/t_vals
    call svd(ac,uc,ref,vhc,status)
    if (status /= 0) error stop 'complex svd warmup'
    call svdvals(ac,s,status)
    if (status /= 0) error stop 'complex svdvals warmup'
    err = maxval(abs(ref-s))/max(1.0_dp,ref(1))
    if (err > 1.0e-11_dp) error stop 'complex spectrum parity'
    call cpu_time(start)
    do j = 1, repetitions
        call svd(ac,uc,ref,vhc,status)
        if (status /= 0) error stop 'complex svd timed'
    end do
    call cpu_time(finish)
    t_svd = finish-start
    call cpu_time(start)
    do j = 1, repetitions
        call svdvals(ac,s,status)
        if (status /= 0) error stop 'complex svdvals timed'
    end do
    call cpu_time(finish)
    t_vals = finish-start
    print '(A,ES12.4,A,ES12.4,A,F8.3)', 'complex: svd=', t_svd, &
        ' svdvals=',t_vals,' speedup=',t_svd/t_vals
    print '(A,ES12.4)', 'singular-value agreement: ',err
end program bench_svdvals
