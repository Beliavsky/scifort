! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_envelope_deconvolve_public
    use scifort_signal, only : dp, deconvolve, envelope, &
        deconvolution_status_success, deconvolution_status_invalid_input, &
        envelope_status_success, envelope_status_invalid_input
    implicit none
    real(dp), allocatable :: q(:), r(:), y(:), remainder(:)
    integer :: status

    call deconvolve([1.0_dp,3.0_dp,2.0_dp], [1.0_dp,1.0_dp], q, r, status)
    if (status /= deconvolution_status_success) error stop 'public deconvolution export'
    if (maxval(abs(q-[1.0_dp,2.0_dp])) > 1.0e-12_dp) error stop 'public deconvolution result'
    call deconvolve([1.0_dp], [0.0_dp], q, r, status)
    if (status /= deconvolution_status_invalid_input) error stop 'public deconvolution error'

    y = envelope([1.0_dp,1.0_dp,1.0_dp,1.0_dp], bp_low=0, bp_high=1, &
        remainder=remainder, status=status)
    if (status /= envelope_status_success) error stop 'public envelope export'
    if (maxval(abs(y-1.0_dp)) > 1.0e-12_dp) error stop 'public envelope result'
    y = envelope([1.0_dp,1.0_dp], bp_low=2, status=status)
    if (status /= envelope_status_invalid_input) error stop 'public envelope error'

    print *, 'PASS test_signal_envelope_deconvolve_public'
end program test_signal_envelope_deconvolve_public
