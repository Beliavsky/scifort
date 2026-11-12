! SPDX-License-Identifier: MIT
program test_linalg_balance
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_balance, only : matrix_balance
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: a(5,5), b(5,5), t(5,5), fac(5)
    real(dp) :: wrong_b(4,4), wrong_t(3,3), wrong_fac(4), e(0,0), be(0,0), te(0,0)
    real(dp) :: huge_a(5,5), small_a(5,5), a0(5,5), b0(5,5), nan_a(5,5)
    complex(dp) :: ac(5,5), bc(5,5), tc(5,5), anc(5,5)
    integer :: perm(5), wrong_perm(4), stat, i, j

    a = 0.0_dp
    a(1,1) = 7.0_dp
    a(1,5) = 2.0_dp
    a(5,5) = -3.0_dp
    do j = 2, 4
        do i = 2, 4
            a(i,j) = real(i * (3 + j) + 2 * j,dp)
        end do
    end do
    a(2,3) = 2.0_dp**22
    a(3,2) = 2.0_dp**(-22)
    call matrix_balance(a,b,stat,t,scale_factors=fac,permutation=perm)
    if (stat /= linalg_status_success) error stop 'real balance status'
    if (all(perm == [1,2,3,4,5])) error stop 'permutation path not exercised'
    do i = 1, 5
        if (abs(fac(i) * fac(i)) <= 0.0_dp) error stop 'invalid scale factors'
        if (fac(i) /= 2.0_dp**nint(log(fac(i))/log(2.0_dp))) &
            error stop 'nonbinary scale factor'
        do j = 1, 5
            if (t(j,i) /= merge(fac(i),0.0_dp,j==perm(i))) &
                error stop 'invalid transformation matrix'
        end do
    end do
    if (maxval(abs(matmul(a,t)-matmul(t,b))) > &
        1.0e-12_dp * max(1.0_dp,maxval(abs(matmul(a,t))))) &
        error stop 'real similarity identity'

    call matrix_balance(a,b,stat,t,permute=.false.,scale=.false., &
        scale_factors=fac,permutation=perm)
    if (stat /= linalg_status_success) error stop 'identity options status'
    if (any(b /= a)) error stop 'identity options altered input'
    if (any(perm /= [1,2,3,4,5]) .or. any(fac /= 1.0_dp)) &
        error stop 'identity options returned transforms'
    do i = 1, 5
        do j = 1, 5
            if (t(j,i) /= merge(1.0_dp,0.0_dp,j==i)) error stop 'identity t'
        end do
    end do
    call matrix_balance(a,b,stat,permute=.false.,scale_factors=fac,permutation=perm)
    if (stat /= linalg_status_success) error stop 'scaling-only status'
    if (any(perm /= [1,2,3,4,5])) error stop 'scaling-only permutation'
    if (all(fac == 1.0_dp)) error stop 'scaling not exercised'
    call matrix_balance(a,b,stat,scale=.false.,scale_factors=fac,permutation=perm)
    if (stat /= linalg_status_success) error stop 'permutation-only status'
    if (any(fac /= 1.0_dp)) error stop 'permutation-only factors'
    do i = 1, 5
        do j = 1, 5
            if (b(i,j) /= a(perm(i),perm(j))) error stop 'permutation-only transform'
        end do
    end do

    ac = cmplx(a, 0.17_dp*a, dp)
    call matrix_balance(ac,bc,stat,tc,scale_factors=fac,permutation=perm)
    if (stat /= linalg_status_success) error stop 'complex balance status'
    if (maxval(abs(matmul(ac,tc)-matmul(tc,bc))) > &
        1.0e-12_dp * max(1.0_dp,maxval(abs(matmul(ac,tc))))) &
        error stop 'complex similarity identity'

    call matrix_balance(e,be,stat,te)
    if (stat /= linalg_status_success) error stop 'empty status'
    call matrix_balance(a,wrong_b,stat)
    if (stat /= linalg_status_invalid_shape) error stop 'wrong B shape'
    call matrix_balance(a,b,stat,wrong_t)
    if (stat /= linalg_status_invalid_shape) error stop 'wrong T shape'
    call matrix_balance(a,b,stat,scale_factors=wrong_fac)
    if (stat /= linalg_status_invalid_shape) error stop 'wrong factors shape'
    call matrix_balance(a,b,stat,permutation=wrong_perm)
    if (stat /= linalg_status_invalid_shape) error stop 'wrong permutation shape'
    nan_a = a
    nan_a(2,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    call matrix_balance(nan_a,b,stat,t)
    if (stat /= linalg_status_no_convergence) error stop 'NaN real status'
    if (any(b /= 0.0_dp) .or. any(t /= 0.0_dp)) error stop 'NaN outputs'
    anc = ac
    anc(2,2) = cmplx(1.0_dp,ieee_value(0.0_dp,ieee_quiet_nan),dp)
    call matrix_balance(anc,bc,stat,tc)
    if (stat /= linalg_status_no_convergence) error stop 'NaN complex status'
    if (any(bc /= cmplx(0.0_dp,0.0_dp,dp))) error stop 'NaN complex output'

    a0 = a
    a0(1,:) = a0(1,:) * 0.001_dp
    huge_a = a0 * 1.0e280_dp
    small_a = a0 * 1.0e-280_dp
    call matrix_balance(huge_a,b,stat)
    if (stat /= linalg_status_success) error stop 'huge status'
    call matrix_balance(a0,b0,stat)
    if (stat /= linalg_status_success) error stop 'base status'
    if (maxval(abs(b/1.0e280_dp-b0)) > 1.0e-10_dp * maxval(abs(b0))) &
        error stop 'huge scale invariance'
    call matrix_balance(small_a,b,stat)
    if (stat /= linalg_status_success) error stop 'tiny status'
    if (maxval(abs(b/1.0e-280_dp-b0)) > 1.0e-10_dp * maxval(abs(b0))) &
        error stop 'small scale invariance'

    print *, 'matrix_balance native tests passed'
end program test_linalg_balance
