! SPDX-License-Identifier: MIT
program test_linalg_svd_condition
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_svd, only : svdvals, cond2, matrix_rank
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    implicit none
    real(dp) :: a(3,3), s(3), cond, wide(2,4), tall(4,2), &
        empty(0,3), absent(0)
    complex(dp) :: z(3,3), zz(2,4)
    integer :: status, rank, i

    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 0.5_dp
    a(3,3) = 0.1_dp
    call svdvals(a,s,status)
    if (status /= linalg_status_success) error stop 'svdvals real status'
    if (maxval(abs(s-[1.0_dp,0.5_dp,0.1_dp])) > 1e-14_dp) &
        error stop 'svdvals real diagonal'
    call cond2(a,cond,status)
    if (status /= linalg_status_success .or. abs(cond-10.0_dp) > 1e-12_dp) &
        error stop 'real cond2'
    call matrix_rank(a,rank,status)
    if (status /= linalg_status_success .or. rank /= 3) &
        error stop 'real matrix_rank'
    call matrix_rank(a,rank,status,tol=0.2_dp)
    if (status /= linalg_status_success .or. rank /= 2) &
        error stop 'real matrix_rank threshold'

    z = cmplx(0.0_dp,0.0_dp,dp)
    z(1,1) = cmplx(0.0_dp,2.0_dp,dp)
    z(2,2) = cmplx(-1.0_dp,0.0_dp,dp)
    z(3,3) = cmplx(0.0_dp,0.5_dp,dp)
    call svdvals(z,s,status)
    if (status /= linalg_status_success) error stop 'svdvals complex status'
    if (maxval(abs(s-[2.0_dp,1.0_dp,0.5_dp])) > 1e-14_dp) &
        error stop 'svdvals complex diagonal'
    call cond2(z,cond,status)
    if (status /= linalg_status_success .or. abs(cond-4.0_dp) > 1e-12_dp) &
        error stop 'complex cond2'
    call matrix_rank(z,rank,status,tol=1.1_dp)
    if (status /= linalg_status_success .or. rank /= 1) &
        error stop 'complex matrix_rank custom tol'

    wide = 0.0_dp
    tall = 0.0_dp
    zz = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1, 2
        wide(i,i) = 1.0_dp
        tall(i,i) = 1.0_dp
        zz(i,i) = cmplx(0.0_dp,1.0_dp,dp)
    end do
    call cond2(wide,cond,status)
    if (status /= 0 .or. abs(cond-1.0_dp) > 1e-12_dp) &
        error stop 'wide cond2'
    call matrix_rank(wide,rank,status)
    if (status /= 0 .or. rank /= 2) error stop 'wide rank'
    call cond2(tall,cond,status)
    if (status /= 0 .or. abs(cond-1.0_dp) > 1e-12_dp) &
        error stop 'tall cond2'
    call matrix_rank(tall,rank,status)
    if (status /= 0 .or. rank /= 2) error stop 'tall rank'
    call cond2(zz,cond,status)
    if (status /= 0 .or. abs(cond-1.0_dp) > 1e-12_dp) &
        error stop 'complex wide cond2'

    ! Exact singularity is indicated by +infinity, not by an error status.
    a(3,3) = 0.0_dp
    call cond2(a,cond,status)
    if (status /= 0 .or. ieee_is_finite(cond) .or. cond < 0.0_dp) &
        error stop 'singular cond2 is positive infinity'
    call matrix_rank(a,rank,status)
    if (status /= 0 .or. rank /= 2) error stop 'singular rank'
    a = 0.0_dp
    call cond2(a,cond,status)
    if (status /= 0 .or. ieee_is_finite(cond)) error stop 'zero matrix cond2'
    call matrix_rank(a,rank,status)
    if (status /= 0 .or. rank /= 0) error stop 'zero matrix rank'

    a = 0.0_dp
    a(1,1) = 1.0e-110_dp
    a(2,2) = 1.0e-111_dp
    a(3,3) = 0.5e-111_dp
    call cond2(a,cond,status)
    if (status /= 0 .or. abs(cond-20.0_dp) > 1e-12_dp) &
        error stop 'small-scale cond2'
    a = a*1.0e220_dp
    call cond2(a,cond,status)
    if (status /= 0 .or. abs(cond-20.0_dp) > 1e-12_dp) &
        error stop 'large-scale cond2'

    call svdvals(a,s(:2),status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid svdvals shape'
    call svdvals(a,s,status,max_sweeps=0)
    if (status /= linalg_status_invalid_shape) error stop 'invalid sweeps'
    call cond2(a,cond,status,max_sweeps=0)
    if (status /= linalg_status_invalid_shape .or. cond /= 0.0_dp) &
        error stop 'invalid condition sweeps'
    call matrix_rank(a,rank,status,tol=-1.0_dp)
    if (status /= linalg_status_invalid_shape .or. rank /= 0) &
        error stop 'invalid matrix_rank tolerance'
    call matrix_rank(z,rank,status,tol=-1.0_dp)
    if (status /= linalg_status_invalid_shape .or. rank /= 0) &
        error stop 'invalid complex rank tolerance'
    call matrix_rank(a,rank,status,max_sweeps=0)
    if (status /= linalg_status_invalid_shape .or. rank /= 0) &
        error stop 'invalid rank sweeps'
    call cond2(empty,cond,status)
    if (status /= linalg_status_invalid_shape .or. cond /= 0.0_dp) &
        error stop 'empty cond2'
    call matrix_rank(empty,rank,status)
    if (status /= linalg_status_invalid_shape .or. rank /= 0) &
        error stop 'empty matrix rank'
    call svdvals(empty,absent,status)
    if (status /= linalg_status_invalid_shape) error stop 'empty svdvals'

    a(1,1) = ieee_value(1.0_dp,ieee_quiet_nan)
    call svdvals(a,s,status)
    if (status /= linalg_status_no_convergence .or. any(s /= 0.0_dp)) &
        error stop 'nonfinite svdvals'
    call cond2(a,cond,status)
    if (status /= linalg_status_no_convergence .or. cond /= 0.0_dp) &
        error stop 'nonfinite cond2'
    call matrix_rank(a,rank,status)
    if (status /= linalg_status_no_convergence .or. rank /= 0) &
        error stop 'nonfinite rank'
    z(1,1) = cmplx(1.0_dp,ieee_value(1.0_dp,ieee_quiet_nan),dp)
    call cond2(z,cond,status)
    if (status /= linalg_status_no_convergence .or. cond /= 0.0_dp) &
        error stop 'nonfinite complex cond2'
    call matrix_rank(z,rank,status)
    if (status /= linalg_status_no_convergence .or. rank /= 0) &
        error stop 'nonfinite complex rank'

    print '(A)', 'singular-value-only, cond2 and matrix_rank native checks passed'
end program test_linalg_svd_condition
