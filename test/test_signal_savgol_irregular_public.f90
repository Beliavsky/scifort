! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_savgol_irregular_public
    use scifort_signal, only : dp, savgol_irregular, savgol_status_success
    implicit none
    real(dp) :: t(9), values(9), u(5)
    complex(dp) :: z(9)
    real(dp), allocatable :: y(:)
    complex(dp), allocatable :: out(:)
    integer :: i, stat

    t = [0.0_dp, 0.13_dp, 0.29_dp, 0.51_dp, 0.68_dp, 0.83_dp, 1.12_dp, 1.3_dp, 1.51_dp]
    values = t**2 - 3.0_dp*t + 0.5_dp
    z = cmplx(values, 2.0_dp*values, dp)
    u = [0.1_dp, 0.2_dp, 0.5_dp, 0.9_dp, 1.1_dp]
    y = savgol_irregular(t, values, 5, 2, deriv=1, query=u, status=stat)
    if (stat /= savgol_status_success) error stop 'public savgol irregular status'
    if (maxval(abs(y - (2.0_dp*u-3.0_dp))) > 1.0e-11_dp) error stop 'public real regression'
    out = savgol_irregular(t, z, 5, 2, deriv=1, query=u, status=stat)
    if (stat /= savgol_status_success) error stop 'public complex status'
    if (maxval(abs(out-cmplx(2.0_dp*u-3.0_dp,4.0_dp*u-6.0_dp,dp))) > 1.0e-11_dp) &
        error stop 'public complex regression'
    print *, 'Irregular Savitzky-Golay public API passed.'
end program test_signal_savgol_irregular_public
