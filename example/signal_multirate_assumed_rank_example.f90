! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program signal_multirate_assumed_rank_example
    use scifort_signal, only : dp, resample_poly_into, resample_status_success
    implicit none
    real(dp) :: x(2,1,9,3), y(2,1,14,3)
    complex(dp) :: complex_y(2,1,14,3), taps(3)
    integer :: i, j, k, status
    do k = 1,3
        do j = 1,9
            do i = 1,2
                x(i,1,j,k) = sin(0.17_dp*real(i+j+k,dp))
            end do
        end do
    end do
    ! Four-dimensional real input, resampled along dimension three.
    call resample_poly_into(x,y,3,2,axis=3,status=status)
    if (status /= resample_status_success) error stop 'real resampling failed'
    taps=cmplx([0.2_dp,0.6_dp,0.2_dp],[0.1_dp,-0.1_dp,0.05_dp],dp)
    ! Complex filter coefficients change the output type even for real input.
    call resample_poly_into(x,complex_y,3,2,taps=taps,axis=3,status=status)
    if (status /= resample_status_success) error stop 'complex FIR resampling failed'
    print '(a,2(1x,es14.6))', 'real and complex result:', y(1,1,1,1), &
        real(complex_y(1,1,1,1),dp)
end program signal_multirate_assumed_rank_example
