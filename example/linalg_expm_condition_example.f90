! SPDX-License-Identifier: MIT
program linalg_expm_condition_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_expm_condition, only : expm_cond
    implicit none
    real(dp) :: a(3,3), cond
    complex(dp) :: b(2,2)
    integer :: status, iterations

    a(1,:) = [-0.3_dp, 0.2_dp, 0.6_dp]
    a(2,:) = [0.6_dp, 0.3_dp, -0.1_dp]
    a(3,:) = [-0.7_dp, 1.2_dp, 0.9_dp]
    call expm_cond(a, cond, status, iterations=iterations)
    if (status /= linalg_status_success) error stop 'real expm_cond failed'
    print '(a,f12.8,a,i0,a)', 'relative condition = ', cond, ' (steps = ', iterations, ')'
    print '(a)', 'SciPy example reference: 1.7787805864469866'

    b(1,:) = [cmplx(0.3_dp,0.2_dp,dp), cmplx(0.2_dp,0.5_dp,dp)]
    b(2,:) = [cmplx(0.4_dp,-0.5_dp,dp), cmplx(-0.1_dp,0.4_dp,dp)]
    call expm_cond(b, cond, status)
    if (status /= linalg_status_success) error stop 'complex expm_cond failed'
    print '(a,f12.8)', 'complex relative condition = ', cond
end program linalg_expm_condition_example
