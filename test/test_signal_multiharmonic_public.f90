! SPDX-License-Identifier: MIT
! Ensure both extensions work through the full public signal facade.
program test_signal_multiharmonic_public
    use scifort_signal, only : dp, lombscargle, lombscargle_multiharmonic, irregular_status_success
    implicit none
    real(dp) :: t(24), y(24)
    real(dp), parameter :: f(2) = [0.8_dp, 1.9_dp]
    real(dp), allocatable :: a(:), b(:)
    integer :: i, status

    do i = 1, size(t)
        t(i) = 0.37_dp*real(i,dp) + 0.02_dp*sin(real(i,dp))
    end do
    y = 1.0_dp + sin(1.9_dp*t) + 0.4_dp*cos(3.8_dp*t)
    a = lombscargle(t, y, f, nterms=2, floating_mean=.true., &
        normalize='normalize', status=status)
    if (status /= irregular_status_success) error stop 'public nterms status'
    b = lombscargle_multiharmonic(t, y, f, 2, floating_mean=.true., &
        normalize='normalize', status=status)
    if (status /= irregular_status_success) error stop 'public direct status'
    if (maxval(abs(a-b)) > 1.0e-12_dp) error stop 'public nterms values'
    if (abs(a(2)-1.0_dp) > 1.0e-12_dp) error stop 'public multi-harmonic reconstruction'
    print '(a)', 'test_signal_multiharmonic_public passed'
end program test_signal_multiharmonic_public
