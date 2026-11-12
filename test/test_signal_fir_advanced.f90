! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_fir_advanced
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_signal, only : dp, firwin, firwin2, firls, minimum_phase, &
        fir_advanced_status_success, fir_advanced_status_invalid_input, fir_advanced_status_singular
    implicit none
    real(dp), allocatable :: taps(:), h(:)
    integer :: st
    real(dp), parameter :: tol = 4.0e-11_dp

    taps = firwin2(9,[0.0_dp,0.35_dp,0.5_dp,1.0_dp], &
        [1.0_dp,1.0_dp,0.0_dp,0.0_dp],status=st)
    call require_success(st)
    call check(taps,[-0.00456239477846077_dp,-0.01570510184853538_dp, &
        0.03904135394049031_dp,0.26493469306373818_dp,0.421875_dp, &
        0.26493469306373807_dp,0.03904135394049028_dp, &
        -0.01570510184853539_dp,-0.00456239477846077_dp],tol)
    taps = firwin2(8,[0.0_dp,0.4_dp,1.0_dp],[1.0_dp,1.0_dp,0.0_dp], &
        window='hann',status=st)
    call require_success(st)
    call check(taps,[0.0_dp,-0.00510759705953223_dp,-0.01404961800472379_dp, &
        0.51952489694569581_dp,0.51952489694569581_dp,-0.01404961800472382_dp, &
        -0.00510759705953222_dp,0.0_dp],tol)

    taps = firwin2(9,[0.0_dp,0.3_dp,0.6_dp,1.0_dp], &
        [0.0_dp,1.0_dp,1.0_dp,0.0_dp],antisymmetric=.true.,window='boxcar',status=st)
    call require_success(st)
    call check(taps,[0.002143810689331689_dp,-0.00479616515449573_dp, &
        0.0432856698119900_dp,0.5138538385140303_dp,0.0_dp, &
        -0.5138538385140302_dp,-0.0432856698119899_dp, &
        0.00479616515449569_dp,-0.00214381068933166_dp],tol)
    if (taps(5) /= 0.0_dp) error stop 'firwin2 type III center tap'
    taps = firwin2(8,[0.0_dp,0.25_dp,0.75_dp,1.0_dp], &
        [0.0_dp,1.0_dp,1.0_dp,1.0_dp],antisymmetric=.true.,status=st)
    call require_success(st)
    call check(taps,[0.00118859122592374_dp,0.01644814112578356_dp, &
        0.11004375131273506_dp,0.5940279057551432_dp, &
        -0.5940279057551433_dp,-0.11004375131273492_dp, &
        -0.01644814112578361_dp,-0.00118859122592373_dp],tol)
    taps = firwin2(11,[0.0_dp,0.35_dp,0.35_dp,1.0_dp], &
        [1.0_dp,1.0_dp,0.0_dp,0.0_dp],status=st)
    call require_success(st)
    call check(taps,[-0.00409957469459569_dp,-0.01266346923989769_dp, &
        -0.00419807006133749_dp,0.09085297892655222_dp, &
        0.25647381471376013_dp,0.34375_dp,0.25647381471376013_dp, &
        0.09085297892655220_dp,-0.00419807006133752_dp, &
        -0.01266346923989769_dp,-0.00409957469459569_dp],tol)

    taps = firls(9,[0.0_dp,0.3_dp,0.45_dp,1.0_dp], &
        [1.0_dp,1.0_dp,0.0_dp,0.0_dp],status=st)
    call require_success(st)
    call check(taps,[-0.06042018796155106_dp,-0.03783847436999917_dp, &
        0.10312389447412734_dp,0.29066384482104035_dp, &
        0.37837913862168027_dp,0.29066384482104035_dp, &
        0.10312389447412734_dp,-0.03783847436999917_dp, &
        -0.06042018796155106_dp],tol)
    taps = firls(11,[0.0_dp,0.2_dp,0.28_dp,0.45_dp,0.62_dp,0.8_dp,0.87_dp,1.0_dp], &
        [1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.7_dp,0.7_dp,0.0_dp,0.0_dp], &
        weight=[1.0_dp,3.0_dp,2.0_dp,4.0_dp],status=st)
    call require_success(st)
    call check(taps,[-0.03206869278587658_dp,-0.03472568320636387_dp, &
        0.22169831743561452_dp,0.07584355728405495_dp,0.10295779099219110_dp, &
        0.43832595167172334_dp,0.10295779099219110_dp, &
        0.07584355728405495_dp,0.22169831743561452_dp, &
        -0.03472568320636387_dp,-0.03206869278587658_dp],tol)
    taps = firls(1,[0.0_dp,1.0_dp],[0.4_dp,0.4_dp],status=st)
    call require_success(st)
    call check(taps,[0.4_dp],tol)

    h = firwin(9,[0.35_dp])
    taps = minimum_phase(h,n_fft=2048,status=st)
    call require_success(st)
    call check(taps,[0.30520735644857777_dp,0.46194254151965880_dp, &
        0.24296762214145329_dp,0.00353767379938179_dp, &
        -0.00845973161614533_dp],tol)
    taps = minimum_phase(h,n_fft=2048,half=.false.,status=st)
    call require_success(st)
    call check(taps,[0.09315157830838994_dp,0.28197656599209314_dp, &
        0.36170192988272065_dp,0.22663361232390250_dp,0.05713771894503718_dp, &
        -0.01132412045083919_dp,-0.00891485066165325_dp, &
        -0.00078881539921613_dp,0.00042645055903873_dp],tol)
    taps = minimum_phase(h,method='hilbert',n_fft=2048,status=st)
    call require_success(st)
    call check(taps,[0.29691424482051576_dp,0.47954605304088338_dp, &
        0.22299508622574232_dp,0.02160803678193611_dp, &
        -0.02104666546867615_dp],tol)

    ! The default FFT length is 2048 for a 9-tap input.
    taps = minimum_phase(h,status=st)
    call require_success(st)
    call check(taps,[0.30520735644857777_dp,0.46194254151965880_dp, &
        0.24296762214145329_dp,0.00353767379938179_dp, &
        -0.00845973161614533_dp],tol)
    ! Non-power-of-two FFT length exercises the direct-DFT fallback.
    taps = minimum_phase(h,n_fft=257,status=st)
    call require_success(st)
    call check(taps,[0.30315920672023394_dp,0.46243476387356319_dp, &
        0.24453826329097989_dp,0.00330121277526493_dp, &
        -0.00796719708559603_dp],tol)
    taps = firwin2(1,[0.0_dp,1.0_dp],[1.0_dp,0.0_dp],window='none',status=st)
    call require_success(st)
    call check(taps,[0.5_dp],tol)
    taps = firwin2(1,[0.0_dp,1.0_dp],[0.0_dp,0.0_dp], &
        window='none',antisymmetric=.true.,status=st)
    call require_success(st)
    call check(taps,[0.0_dp],tol)

    ! Invalid designs return empty coefficients and the corresponding status.
    taps = firwin2(6,[0.0_dp,1.0_dp],[1.0_dp,1.0_dp],status=st)
    call require_invalid(taps,st)
    taps = firwin2(7,[0.0_dp,1.0_dp],[1.0_dp,0.0_dp],nfreqs=7,status=st)
    call require_invalid(taps,st)
    taps = firwin2(5,[0.0_dp,0.0_dp,1.0_dp], &
        [1.0_dp,1.0_dp,0.0_dp],status=st)
    call require_invalid(taps,st)
    taps = firwin2(5,[0.0_dp,0.4_dp,0.4_dp,0.4_dp,1.0_dp], &
        [1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp],status=st)
    call require_invalid(taps,st)
    taps = firwin2(7,[0.0_dp,1.0_dp],[0.1_dp,0.0_dp], &
        antisymmetric=.true.,status=st)
    call require_invalid(taps,st)
    taps = firwin2(7,[0.0_dp,1.0_dp],[1.0_dp,0.0_dp], &
        window='not-a-window',status=st)
    call require_invalid(taps,st)
    taps = firls(8,[0.0_dp,0.4_dp],[1.0_dp,0.0_dp],status=st)
    call require_invalid(taps,st)
    taps = firls(9,[0.0_dp,0.4_dp],[1.0_dp,0.0_dp], &
        weight=[1.0_dp,2.0_dp],status=st)
    call require_invalid(taps,st)
    taps = firls(9,[0.0_dp,0.4_dp],[1.0_dp,0.0_dp], &
        weight=[0.0_dp],status=st)
    if (st /= fir_advanced_status_singular .or. size(taps) /= 0) &
        error stop 'expected singular FIRLS system'
    taps = minimum_phase([0.0_dp,0.0_dp,0.0_dp],n_fft=8,status=st)
    if (st /= fir_advanced_status_singular .or. size(taps) /= 0) &
        error stop 'expected singular minimum phase filter'
    taps = minimum_phase(h,n_fft=8,status=st)
    call require_invalid(taps,st)
    taps = minimum_phase(h,method='hilbert',half=.false.,status=st)
    call require_invalid(taps,st)
    taps = minimum_phase(h,method='wrong',status=st)
    call require_invalid(taps,st)
    print '(a)', 'PASS: advanced FIR designs and minimum-phase conversion'
contains
    subroutine check(actual,expected,atol)
        real(dp), intent(in) :: actual(:) !! computed coefficients
        real(dp), intent(in) :: expected(:) !! independently computed SciPy reference
        real(dp), intent(in) :: atol !! absolute comparison tolerance
        if (size(actual) /= size(expected)) error stop 'unexpected FIR array size'
        if (.not. all(ieee_is_finite(actual))) error stop 'nonfinite FIR coefficients'
        if (maxval(abs(actual-expected)) > atol) then
            print *, 'actual: ',actual
            print *, 'expected: ',expected
            error stop 'FIR mismatch against SciPy reference'
        end if
    end subroutine check

    subroutine require_success(status)
        integer, intent(in) :: status !! returned status code
        if (status /= fir_advanced_status_success) error stop 'advanced FIR status'
    end subroutine require_success

    subroutine require_invalid(coeffs,status)
        real(dp), intent(in) :: coeffs(:) !! expected empty return
        integer, intent(in) :: status !! expected invalid status
        if (status /= fir_advanced_status_invalid_input .or. size(coeffs) /= 0) &
            error stop 'advanced FIR invalid input not rejected'
    end subroutine require_invalid
end program test_signal_fir_advanced
