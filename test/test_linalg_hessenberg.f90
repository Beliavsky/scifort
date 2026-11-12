! SPDX-License-Identifier: MIT
program test_linalg_hessenberg
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_hessenberg, only : hessenberg
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(5,5), h(5,5), q(5,5), ident(5,5), nan_a(5,5)
    real(dp) :: bad_h(4,4), bad_q(4,4), bad_a(3,4), only_h(5,5)
    real(dp) :: empty_r(0,0), empty_h(0,0), empty_q(0,0)
    complex(dp) :: ac(5,5), hc(5,5), qc(5,5), idc(5,5), only_hc(5,5)
    integer :: i, j, stat
    a = 0.0_dp
    ident = 0.0_dp
    do j = 1, 5
        ident(j,j) = 1.0_dp
        do i = 1, 5
            a(i,j) = real(3*i-j, dp)/real(1+i+j, dp)
        end do
    end do
    call hessenberg(a,h,stat,q)
    if (stat /= linalg_status_success) error stop 'real hessenberg status'
    if (maxval(abs(matmul(q,transpose(q))-ident)) > 2.0e-13_dp) &
        error stop 'real Q orthogonality'
    if (maxval(abs(matmul(matmul(q,h),transpose(q))-a)) > 2.0e-12_dp) &
        error stop 'real reconstruction'
    do j = 1, 3
        if (maxval(abs(h(j+2:5,j))) /= 0.0_dp) error stop 'real Hessenberg shape'
    end do
    call hessenberg(a,only_h,stat)
    if (stat /= linalg_status_success) error stop 'real no-Q status'
    if (maxval(abs(only_h-h)) > 2.0e-13_dp) error stop 'real no-Q values'

    ac = cmplx(a,transpose(a)*0.15_dp,dp)
    idc = cmplx(ident,0.0_dp,dp)
    call hessenberg(ac,hc,stat,qc)
    if (stat /= linalg_status_success) error stop 'complex status'
    if (maxval(abs(matmul(qc,transpose(conjg(qc)))-idc)) > 2.0e-13_dp) &
        error stop 'complex Q unitarity'
    if (maxval(abs(matmul(matmul(qc,hc),transpose(conjg(qc)))-ac)) > 2.0e-12_dp) &
        error stop 'complex reconstruction'
    do j = 1, 3
        if (maxval(abs(hc(j+2:5,j))) /= 0.0_dp) error stop 'complex Hessenberg shape'
    end do
    call hessenberg(ac,only_hc,stat)
    if (stat /= linalg_status_success) error stop 'complex no-Q status'
    if (maxval(abs(only_hc-hc)) > 2.0e-13_dp) error stop 'complex no-Q values'

    call hessenberg(a,bad_h,stat)
    if (stat /= linalg_status_invalid_shape) error stop 'invalid output shape status'
    call hessenberg(a,h,stat,bad_q)
    if (stat /= linalg_status_invalid_shape) error stop 'invalid q shape status'
    bad_a = 0.0_dp
    call hessenberg(bad_a,h,stat)
    if (stat /= linalg_status_invalid_shape) error stop 'nonsquare input status'
    call hessenberg(empty_r,empty_h,stat,empty_q)
    if (stat /= linalg_status_success) error stop 'empty shape status'
    nan_a = a
    nan_a(2,3) = ieee_value(0.0_dp,ieee_quiet_nan)
    call hessenberg(nan_a,h,stat,q)
    if (stat /= linalg_status_no_convergence) error stop 'nonfinite status'
    if (maxval(abs(q)) /= 0.0_dp) error stop 'nonfinite q not cleared'
    ac(1,3) = cmplx(0.0_dp, ieee_value(0.0_dp,ieee_quiet_nan),dp)
    call hessenberg(ac,hc,stat,qc)
    if (stat /= linalg_status_no_convergence) error stop 'complex nonfinite status'
    print *, 'hessenberg native tests passed'
end program test_linalg_hessenberg
