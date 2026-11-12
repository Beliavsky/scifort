! SPDX-License-Identifier: MIT
program linalg_circulant_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_circulant, only : circulant, solve_circulant
    implicit none
    real(dp) :: c(3), b(3), x(3), a(3,3)
    real(dp) :: cs(4), bs(4), xs(4)
    integer :: status

    c = [2.0_dp,2.0_dp,4.0_dp]
    b = [1.0_dp,2.0_dp,3.0_dp]
    call circulant(c,a,status)
    if (status /= linalg_status_success) error stop 'circulant failed'
    call solve_circulant(c,b,x,status)
    if (status /= linalg_status_success) error stop 'solve_circulant failed'
    print '(a)', 'SciPy solve_circulant example:'
    print '(a,*(f10.5,1x))', 'Solution: ', x
    print '(a,es12.4)', 'Residual: ', maxval(abs(matmul(a,x)-b))

    cs = [1.0_dp,1.0_dp,0.0_dp,0.0_dp]
    bs = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    call solve_circulant(cs,bs,xs,status,singular='lstsq')
    if (status /= linalg_status_success) error stop 'least squares failed'
    print '(a,*(f10.5,1x))', 'Minimum-norm least squares: ', xs
end program linalg_circulant_example
