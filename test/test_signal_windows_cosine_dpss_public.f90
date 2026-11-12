! SPDX-License-Identifier: MIT
! Public SciFort signal umbrella exports for flexible windows.
program test_signal_windows_cosine_dpss_public
    use scifort_signal, only : dp, get_window, general_cosine, spectral_status_success
    implicit none
    real(dp), allocatable :: w(:), v(:)
    integer :: code

    w = general_cosine(16,[0.5_dp,0.5_dp],sym=.false.,status=code)
    if (code /= spectral_status_success) error stop 'public general_cosine status'
    v = get_window('general_cosine',16,coefficients=[0.5_dp,0.5_dp],status=code)
    if (code /= spectral_status_success) error stop 'public cosine get_window status'
    if (any(abs(w-v)>1.0e-13_dp)) error stop 'public cosine values'
    w = get_window('dpss',16,parameter=2.5_dp,status=code)
    if (code /= spectral_status_success) error stop 'public dpss status'
    if (size(w) /= 16 .or. minval(w) <= 0.0_dp) error stop 'public dpss values'
    print '(a)', 'test_signal_windows_cosine_dpss_public passed'
end program test_signal_windows_cosine_dpss_public
