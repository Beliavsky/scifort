! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_firls_rank
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_signal, only : dp, firls, fir_advanced_status_success, &
        fir_advanced_status_singular
    implicit none
    real(dp), allocatable :: taps(:), ref(:)
    real(dp) :: e, freq, target
    integer :: status, n, i, j
    real(dp), parameter :: pi = acos(-1.0_dp)

    ! These valid narrow bands used to be rejected by the Cholesky solver.
    do j = 1, 3
        select case (j)
        case (1)
            n = 9
        case (2)
            n = 31
        case default
            n = 51
        end select
        taps = firls(n,[0.4_dp,0.405_dp],[1.0_dp,1.0_dp],status=status)
        call check_status(taps,n,status)
        e = 0.0_dp
        do i = 0, 20
            freq = 0.4_dp + 0.005_dp*real(i,dp)/20.0_dp
            e = max(e,abs(zero_phase_gain(taps,freq)-1.0_dp))
        end do
        if (e > 2.0e-5_dp) error stop 'narrowband least-squares response inaccurate'
    end do

    taps = firls(31,[0.4_dp,0.4001_dp],[1.0_dp,1.0_dp],status=status)
    call check_status(taps,31,status)
    do i = 0, 20
        freq = 0.4_dp + 0.0001_dp*real(i,dp)/20.0_dp
        if (abs(zero_phase_gain(taps,freq)-1.0_dp) > 2.0e-5_dp) &
            error stop 'almost rank-deficient system not solved'
    end do

    ! The linear target slopes across the band; the solution need not have
    ! small coefficients because the response changes rapidly in frequency.
    taps = firls(19,[0.10855911_dp,0.14932278_dp], &
        [0.52381273_dp,1.46214332_dp],status=status)
    call check_status(taps,19,status)
    do i = 0, 30
        freq = 0.10855911_dp + (0.14932278_dp-0.10855911_dp)*real(i,dp)/30.0_dp
        target = 0.52381273_dp + (1.46214332_dp-0.52381273_dp)*real(i,dp)/30.0_dp
        if (abs(zero_phase_gain(taps,freq)-target) > 1.0e-3_dp) &
            error stop 'tilted narrowband response inaccurate'
    end do

    ! Rescaling weights by a common factor must leave the filter unchanged,
    ! even for weights large enough to overflow unscaled Gram arithmetic.
    ref = firls(11,[0.0_dp,0.2_dp,0.3_dp,1.0_dp], &
        [1.0_dp,0.8_dp,0.0_dp,0.0_dp],weight=[1.0_dp,4.0_dp],status=status)
    call check_status(ref,11,status)
    taps = firls(11,[0.0_dp,0.2_dp,0.3_dp,1.0_dp], &
        [1.0_dp,0.8_dp,0.0_dp,0.0_dp],weight=[1.0e280_dp,4.0e280_dp],status=status)
    call check_status(taps,11,status)
    if (maxval(abs(taps-ref)) > 3.0e-11_dp) error stop 'firls scaling invariance'

    ! A zero-weight band must not contribute to the design.
    taps = firls(11,[0.0_dp,0.2_dp,0.3_dp,1.0_dp], &
        [1.0_dp,1.0_dp,0.0_dp,0.0_dp],weight=[0.0_dp,2.0_dp],status=status)
    call check_status(taps,11,status)
    ref = firls(11,[0.3_dp,1.0_dp],[0.0_dp,0.0_dp],status=status)
    call check_status(ref,11,status)
    if (maxval(abs(taps-ref)) > 2.0e-10_dp) error stop 'zero-weight band affected design'

    taps = firls(31,[0.4_dp,0.4_dp],[1.0_dp,1.0_dp],status=status)
    if (size(taps) /= 0) error stop 'zero-width band accepted'
    taps = firls(15,[0.4_dp,0.41_dp],[1.0_dp,1.0_dp], &
        weight=[0.0_dp],status=status)
    if (status /= fir_advanced_status_singular .or. size(taps) /= 0) &
        error stop 'all-zero weights not rejected'

    print '(a)', 'PASS: rank-revealing FIRLS, narrow bands, slopes, scaled weights'
contains
    real(dp) function zero_phase_gain(h, freq) result(gain)
        real(dp), intent(in) :: h(:) !! symmetric FIR taps
        real(dp), intent(in) :: freq !! normalized frequency relative to Nyquist
        integer :: k, center
        center = (size(h)+1)/2
        gain = 0.0_dp
        do k = 1, size(h)
            gain = gain + h(k)*cos(pi*freq*real(k-center,dp))
        end do
    end function zero_phase_gain

    subroutine check_status(h, expected_n, status)
        real(dp), intent(in) :: h(:) !! returned FIR taps
        integer, intent(in) :: expected_n !! required tap count
        integer, intent(in) :: status !! expected success status
        if (status /= fir_advanced_status_success) error stop 'firls failed valid narrow band'
        if (size(h) /= expected_n) error stop 'firls invalid output size'
        if (.not. all(ieee_is_finite(h))) error stop 'firls nonfinite coefficients'
        if (maxval(abs(h-h(expected_n:1:-1))) > 1.0e-12_dp) &
            error stop 'firls symmetry error'
    end subroutine check_status
end program test_signal_firls_rank
