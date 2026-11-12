! SPDX-License-Identifier: MIT
program test_signal_var1_taper_spectrum_public
    use scifort_signal, only : dp, var1_expected_taper_spectrum, &
        var1_taper_spectrum_result, var1_spectrum_success
    implicit none
    type(var1_taper_spectrum_result) :: result
    real(dp) :: a(2,2), sigma(2,2)

    a = 0.0_dp
    sigma = reshape([1.0_dp,0.2_dp,0.2_dp,1.0_dp],[2,2])
    result = var1_expected_taper_spectrum(a,sigma,[1.0_dp,0.0_dp,1.0_dp],[0.0_dp,0.25_dp], &
        normalized=.true.)
    if (result%status /= var1_spectrum_success .or. size(result%power,1) /= 2) &
        error stop 'VAR(1) public interface'
    if (abs(result%power(1,1,2)-0.2_dp)>1.0e-12_dp) &
        error stop 'VAR(1) public white-noise spectrum'
    print '(a)', 'PASS public VAR(1) finite-window spectrum interface'
end program test_signal_var1_taper_spectrum_public
