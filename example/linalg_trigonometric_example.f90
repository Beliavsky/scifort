! SPDX-License-Identifier: MIT
program linalg_trigonometric_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_trigonometric, only : sinm, cosm, tanm, sinhm, coshm, tanhm
    implicit none
    real(dp) :: a(2, 2), s(2, 2), c(2, 2), t(2, 2)
    complex(dp) :: z(2, 2), sz(2, 2)
    integer :: status

    a(1, :) = [0.3_dp, 0.4_dp]
    a(2, :) = [0.0_dp, -0.2_dp]
    call sinm(a, s, status)
    if (status /= linalg_status_success) error stop 'sinm failed'
    call cosm(a, c, status)
    if (status /= linalg_status_success) error stop 'cosm failed'
    call tanm(a, t, status)
    if (status /= linalg_status_success) error stop 'tanm failed'
    print '(a, es12.4)', '||tan(A)*cos(A)-sin(A)||_max = ', maxval(abs(matmul(t,c)-s))

    call sinhm(a, s, status)
    if (status /= linalg_status_success) error stop 'sinhm failed'
    call coshm(a, c, status)
    if (status /= linalg_status_success) error stop 'coshm failed'
    call tanhm(a, t, status)
    if (status /= linalg_status_success) error stop 'tanhm failed'
    print '(a, es12.4)', '||tanh(A)*cosh(A)-sinh(A)||_max = ', maxval(abs(matmul(t,c)-s))

    z = cmplx(a, kind=dp)
    z(1, 1) = cmplx(0.3_dp, 0.4_dp, dp)
    call sinm(z, sz, status)
    if (status /= linalg_status_success) error stop 'complex sinm failed'
    print '(a, 2es14.5)', 'complex sinm(A)(1,1) = ', real(sz(1,1),dp), aimag(sz(1,1))
end program linalg_trigonometric_example
