! SPDX-License-Identifier: MIT
program test_signal_var1_taper_spectrum
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_var1_taper_spectrum, only : var1_taper_spectrum_result, &
        var1_expected_taper_spectrum, var1_spectrum_success, &
        var1_spectrum_invalid_input, var1_spectrum_nonstationary
    implicit none
    type(var1_taper_spectrum_result) :: r, scaled, bad
    real(dp) :: a(2,2), sigma(2,2), taper(8), freq(3), expected(2,2)
    real(dp) :: scalar_a(1,1), scalar_sigma(1,1), scalar_taper(3), exact
    integer :: i, j, f

    a = 0.0_dp
    sigma = reshape([2.0_dp,0.4_dp,0.4_dp,1.0_dp],[2,2])
    taper = [(real(i,dp)/8.0_dp,i=1,8)]
    freq = [-0.33_dp,0.0_dp,0.40_dp]
    r = var1_expected_taper_spectrum(a,sigma,taper,freq)
    call check(r%status == var1_spectrum_success,'white noise accepted')
    call check(maxval(abs(r%lag0_covariance-sigma)) < 1.0e-13_dp,'white covariance')
    do f = 1, 3
        expected = sum(taper**2)*sigma
        call check(maxval(abs(r%power(f,:,:)-expected)) < 1.0e-12_dp, &
            'white spectral matrix independent of frequency')
    end do
    scaled = var1_expected_taper_spectrum(a,sigma,taper,freq,normalized=.true.)
    call check(maxval(abs(scaled%power-r%power/sum(taper**2))) < 1.0e-13_dp, &
        'optional taper norm')

    scalar_a = 0.5_dp
    scalar_sigma = 1.2_dp
    scalar_taper = [1.0_dp,2.0_dp,3.0_dp]
    r = var1_expected_taper_spectrum(scalar_a,scalar_sigma,scalar_taper,[0.0_dp])
    exact = scalar_sigma(1,1)/(1.0_dp-scalar_a(1,1)**2)* &
        (sum(scalar_taper**2)+2.0_dp*sum(scalar_taper(1:2)*scalar_taper(2:3))*0.5_dp+ &
        2.0_dp*scalar_taper(1)*scalar_taper(3)*0.5_dp**2)
    call check(abs(real(r%power(1,1,1),dp)-exact)<1.0e-12_dp,'scalar AR(1) exact')

    a = reshape([0.65_dp,0.0_dp,0.28_dp,0.43_dp],[2,2])
    r = var1_expected_taper_spectrum(a,sigma,taper,freq)
    call check(r%status == var1_spectrum_success,'coupled VAR(1) accepted')
    expected = matmul(a,matmul(r%lag0_covariance,transpose(a)))+sigma
    call check(maxval(abs(expected-r%lag0_covariance))<1.0e-11_dp,'Lyapunov identity')
    do f = 1, 3
        do i = 1, 2
            do j = 1, 2
                call check(abs(r%power(f,i,j)-conjg(r%power(f,j,i)))<1.0e-11_dp, &
                    'Hermitian output')
            end do
        end do
    end do
    call check(abs(aimag(r%power(1,1,2)))>1.0e-4_dp,'nontrivial complex cross phase')

    bad = var1_expected_taper_spectrum(a,sigma,[0.0_dp,0.0_dp],[0.0_dp])
    call check(bad%status == var1_spectrum_invalid_input,'zero-energy taper invalid')
    bad = var1_expected_taper_spectrum(a,reshape([1.0_dp,0.4_dp,0.0_dp,1.0_dp], &
        [2,2]),taper,freq)
    call check(bad%status == var1_spectrum_invalid_input,'nonsymmetric noise invalid')
    bad = var1_expected_taper_spectrum(a,reshape([1.0_dp,2.0_dp,2.0_dp,1.0_dp], &
        [2,2]),taper,freq)
    call check(bad%status == var1_spectrum_invalid_input,'indefinite noise invalid')
    bad = var1_expected_taper_spectrum(a,sigma,taper, &
        [ieee_value(0.0_dp,ieee_quiet_nan)])
    call check(bad%status == var1_spectrum_invalid_input,'nonfinite frequencies invalid')
    bad = var1_expected_taper_spectrum(reshape([1.0_dp,0.0_dp,0.0_dp,0.5_dp], &
        [2,2]),sigma,taper,freq)
    call check(bad%status == var1_spectrum_nonstationary,'unit root rejected')
    bad = var1_expected_taper_spectrum(reshape([1.1_dp,0.0_dp,0.0_dp,0.5_dp], &
        [2,2]),sigma,taper,freq)
    call check(bad%status == var1_spectrum_nonstationary,'explosive VAR rejected')
    print '(a)', 'PASS VAR(1) exact finite-window spectral targets'

contains
    subroutine check(condition,message)
        logical, intent(in) :: condition !! property to verify
        character(*), intent(in) :: message !! context on failure
        if (.not. condition) then
            print *, 'FAIL ',message
            error stop 1
        end if
    end subroutine check
end program test_signal_var1_taper_spectrum
