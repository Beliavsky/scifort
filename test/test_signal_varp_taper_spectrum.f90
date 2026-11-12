! SPDX-License-Identifier: MIT
program test_signal_varp_taper_spectrum
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_var1_taper_spectrum, only : var1_expected_taper_spectrum, &
        var1_taper_spectrum_result
    use scifort_signal_varp_taper_spectrum, only : varp_expected_taper_spectrum, &
        varp_expected_band_taper_spectrum, varp_taper_spectrum_result, &
        varp_spectrum_success, varp_spectrum_nonstationary, varp_spectrum_invalid_input
    implicit none
    type(varp_taper_spectrum_result) :: result, bands_result, bad, normalized
    type(var1_taper_spectrum_result) :: legacy
    real(dp) :: a(2,2,3), sigma(2,2), taper(8), freq(4), bands(2,3)
    real(dp) :: scalar_a(1,1,2), scalar_sigma(1,1), scalar_taper(3)
    real(dp) :: variance, expected, phi1, phi2
    integer :: i, j, f

    a = 0.0_dp
    sigma = reshape([2.0_dp,0.4_dp,0.4_dp,1.5_dp],[2,2])
    taper = [(real(i,dp)/8.0_dp,i=1,8)]
    freq = [-0.44_dp,0.0_dp,0.12_dp,0.45_dp]
    result = varp_expected_taper_spectrum(a,sigma,taper,freq)
    call check(result%status == varp_spectrum_success,'white VAR(3) accepted')
    call check(maxval(abs(result%lag0_covariance-sigma))<1.0e-12_dp, &
        'white VAR(3) covariance')
    do f = 1, size(freq)
        call check(maxval(abs(result%power(f,:,:)-sigma*sum(taper**2)))<1.0e-11_dp, &
            'white VAR(p) expected window power')
    end do

    bands(:,1) = [-0.5_dp,0.5_dp]
    bands(:,2) = [0.05_dp,0.45_dp]
    bands(:,3) = [-0.3_dp,0.1_dp]
    bands_result = varp_expected_band_taper_spectrum(a,sigma,taper,bands)
    call check(bands_result%status == varp_spectrum_success,'exact bands accepted')
    call check(allocated(bands_result%bands),'band metadata returned')
    call check(maxval(abs(bands_result%bands-bands))==0.0_dp,'band metadata exact')
    do f = 1, size(bands,2)
        call check(maxval(abs(bands_result%power(f,:,:)- &
            sigma*sum(taper**2)*(bands(2,f)-bands(1,f))))<1.0e-11_dp, &
            'white VAR(p) exact band power')
    end do
    normalized = varp_expected_band_taper_spectrum(a,sigma,taper,bands, &
        average=.true.,normalized=.true.)
    do f = 1, size(bands,2)
        call check(maxval(abs(normalized%power(f,:,:)-sigma))<1.0e-12_dp, &
            'white band mean and taper normalization')
    end do

    a = 0.0_dp
    a(:,:,1) = reshape([0.65_dp,0.12_dp,0.25_dp,0.48_dp],[2,2])
    legacy = var1_expected_taper_spectrum(a(:,:,1),sigma,taper,freq)
    result = varp_expected_taper_spectrum(a,sigma,taper,freq)
    call check(result%status == varp_spectrum_success,'padded VAR(1) accepted')
    call check(maxval(abs(legacy%lag0_covariance-result%lag0_covariance))<1.0e-12_dp, &
        'padded VAR(1) lag-zero compatibility')
    call check(maxval(abs(legacy%power-result%power))<1.0e-11_dp, &
        'padded VAR(1) spectral compatibility')
    call check(abs(aimag(result%power(1,1,2)))>1.0e-5_dp, &
        'complex channel phase preserved')
    do f = 1, size(freq)
        do j = 1, 2
            do i = 1, 2
                call check(abs(result%power(f,i,j)-conjg(result%power(f,j,i)))<1.0e-11_dp, &
                    'Hermitian spectral output')
            end do
        end do
    end do

    bands(:,1) = [-0.5_dp,0.0_dp]
    bands(:,2) = [0.0_dp,0.5_dp]
    bands(:,3) = [-0.5_dp,0.5_dp]
    bands_result = varp_expected_band_taper_spectrum(a,sigma,taper,bands)
    call check(maxval(abs(bands_result%power(1,:,:)+bands_result%power(2,:,:) &
        -bands_result%power(3,:,:)))<1.0e-11_dp,'adjacent band additivity')
    call check(maxval(abs(bands_result%power(3,:,:)- &
        result%lag0_covariance*sum(taper**2)))<1.0e-11_dp, &
        'full-period integral of correlated VAR')
    normalized = varp_expected_taper_spectrum(a,sigma,taper,freq+1.0_dp)
    call check(maxval(abs(normalized%power-result%power))<1.0e-11_dp, &
        'discrete-time frequency periodicity')

    scalar_a(1,1,1) = 0.55_dp
    scalar_a(1,1,2) = -0.2_dp
    scalar_sigma(1,1) = 1.5_dp
    scalar_taper = [1.0_dp,2.0_dp,3.0_dp]
    result = varp_expected_taper_spectrum(scalar_a,scalar_sigma,scalar_taper,[0.0_dp])
    phi1 = scalar_a(1,1,1)
    phi2 = scalar_a(1,1,2)
    variance = scalar_sigma(1,1)*(1.0_dp-phi2)/ &
        ((1.0_dp+phi2)*((1.0_dp-phi2)**2-phi1**2))
    call check(abs(result%lag0_covariance(1,1)-variance)<1.0e-12_dp, &
        'scalar AR(2) Yule-Walker variance')
    expected = variance*(sum(scalar_taper**2)+ &
        2.0_dp*sum(scalar_taper(1:2)*scalar_taper(2:3))*phi1/(1.0_dp-phi2)+ &
        2.0_dp*scalar_taper(1)*scalar_taper(3)* &
        (phi1*phi1/(1.0_dp-phi2)+phi2))
    call check(abs(real(result%power(1,1,1),dp)-expected)<1.0e-11_dp, &
        'scalar AR(2) finite-taper spectrum')
    bands_result = varp_expected_band_taper_spectrum(scalar_a,scalar_sigma,scalar_taper, &
        reshape([-0.5_dp,0.5_dp],[2,1]))
    call check(abs(real(bands_result%power(1,1,1),dp)- &
        variance*sum(scalar_taper**2))<1.0e-12_dp, &
        'full frequency-period band integrates to time-domain variance')

    bad = varp_expected_taper_spectrum(a,sigma,taper,[ieee_value(0.0_dp,ieee_quiet_nan)])
    call check(bad%status == varp_spectrum_invalid_input,'nonfinite frequency rejected')
    bad = varp_expected_taper_spectrum(a,sigma,[0.0_dp,0.0_dp],freq)
    call check(bad%status == varp_spectrum_invalid_input,'zero-energy taper rejected')
    bad = varp_expected_taper_spectrum(a,reshape([1.0_dp,0.0_dp,0.3_dp,1.0_dp],[2,2]), &
        taper,freq)
    call check(bad%status == varp_spectrum_invalid_input,'nonsymmetric innovations rejected')
    bad = varp_expected_taper_spectrum(a,reshape([1.0_dp,2.0_dp,2.0_dp,1.0_dp],[2,2]), &
        taper,freq)
    call check(bad%status == varp_spectrum_invalid_input,'non-PD innovations rejected')
    bad = varp_expected_band_taper_spectrum(a,sigma,taper, &
        reshape([0.4_dp,0.2_dp],[2,1]))
    call check(bad%status == varp_spectrum_invalid_input,'reversed band rejected')
    bad = varp_expected_band_taper_spectrum(a,sigma,taper, &
        reshape([0.2_dp,0.2_dp],[2,1]))
    call check(bad%status == varp_spectrum_invalid_input,'zero-width band rejected')
    a = 0.0_dp
    a(1,1,1) = 1.0_dp
    a(2,2,1) = 0.25_dp
    bad = varp_expected_taper_spectrum(a,sigma,taper,freq)
    call check(bad%status == varp_spectrum_nonstationary,'unit-root VAR rejected')
    a(1,1,1) = 1.1_dp
    bad = varp_expected_taper_spectrum(a,sigma,taper,freq)
    call check(bad%status == varp_spectrum_nonstationary,'explosive VAR rejected')
    print '(a)', 'PASS VAR(p) exact finite-window and band population targets'

contains
    subroutine check(condition,message)
        logical, intent(in) :: condition !! required assertion
        character(*), intent(in) :: message !! test name
        if (.not. condition) then
            print *, 'FAIL ',message
            error stop 1
        end if
    end subroutine check
end program test_signal_varp_taper_spectrum
