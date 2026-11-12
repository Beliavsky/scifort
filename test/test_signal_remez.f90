! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_remez
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_signal, only : dp, remez, remez_status_success, &
        remez_status_invalid_input, remez_status_no_convergence
    implicit none
    real(dp), allocatable :: taps(:)
    integer :: st, niter
    real(dp) :: ripple

    ! Reference design 1: bandpass, 17 taps.
    call check_reference(17,[0.000000000000000e+00_dp,3.000000000000000e-01_dp,4.000000000000000e-01_dp, &
            1.000000000000000e+00_dp],[1.000000000000000e+00_dp,0.000000000000000e+00_dp],[1.000000000000000e+00_dp, &
            1.000000000000000e+00_dp],'bandpass', &
        [1.032264618598786e-02_dp,7.111385831311677e-02_dp, &
            8.898799975735883e-03_dp,-3.987832771227279e-02_dp,-7.086214909093244e-02_dp,-1.590333344425915e-02_dp, &
            1.266677702877821e-01_dp,2.823815497604519e-01_dp,3.499458652828532e-01_dp,2.823815497604519e-01_dp, &
            1.266677702877821e-01_dp,-1.590333344425915e-02_dp,-7.086214909093244e-02_dp,-3.987832771227279e-02_dp, &
            8.898799975735883e-03_dp,7.111385831311677e-02_dp,1.032264618598786e-02_dp],3.0e-4_dp)
    ! Reference design 2: bandpass, 20 taps.
    call check_reference(20,[0.000000000000000e+00_dp,3.000000000000000e-01_dp,4.000000000000000e-01_dp, &
            1.000000000000000e+00_dp],[1.000000000000000e+00_dp,0.000000000000000e+00_dp],[1.000000000000000e+00_dp, &
            1.000000000000000e+00_dp],'bandpass', &
        [-3.251135928419532e-02_dp,1.713735140414507e-02_dp, &
            4.080128172298027e-02_dp,3.788869547482611e-02_dp,-7.302399402976278e-03_dp,-6.015124288055940e-02_dp, &
            -5.431892711925260e-02_dp,4.861172665675077e-02_dp,2.105885650255683e-01_dp,3.334116418726581e-01_dp, &
            3.334116418726581e-01_dp,2.105885650255683e-01_dp,4.861172665675077e-02_dp,-5.431892711925260e-02_dp, &
            -6.015124288055940e-02_dp,-7.302399402976278e-03_dp,3.788869547482611e-02_dp,4.080128172298027e-02_dp, &
            1.713735140414507e-02_dp,-3.251135928419532e-02_dp],3.0e-4_dp)
    ! Reference design 3: bandpass, 21 taps.
    call check_reference(21,[0.000000000000000e+00_dp,1.800000000000000e-01_dp,2.500000000000000e-01_dp, &
            5.500000000000000e-01_dp,6.500000000000000e-01_dp,1.000000000000000e+00_dp],[0.000000000000000e+00_dp, &
            1.000000000000000e+00_dp,0.000000000000000e+00_dp],[1.000000000000000e+00_dp,2.000000000000000e+00_dp, &
            1.000000000000000e+00_dp],'bandpass', &
        [-1.409488510223903e-02_dp,2.746831839823665e-02_dp, &
            3.311247872427668e-02_dp,1.048377386392366e-01_dp,-9.508460874760599e-03_dp,8.615801187409861e-03_dp, &
            4.467545820552646e-02_dp,-1.561368377556695e-01_dp,-2.482257659120359e-01_dp,1.054270101224847e-01_dp, &
            3.880823499184650e-01_dp,1.054270101224847e-01_dp,-2.482257659120359e-01_dp,-1.561368377556695e-01_dp, &
            4.467545820552646e-02_dp,8.615801187409861e-03_dp,-9.508460874760599e-03_dp,1.048377386392366e-01_dp, &
            3.311247872427668e-02_dp,2.746831839823665e-02_dp,-1.409488510223903e-02_dp],3.0e-4_dp)
    ! Reference design 4: bandpass, 17 taps.
    call check_reference(17,[0.000000000000000e+00_dp,3.000000000000000e-01_dp,4.000000000000000e-01_dp, &
            1.000000000000000e+00_dp],[0.000000000000000e+00_dp,1.000000000000000e+00_dp],[1.000000000000000e+00_dp, &
            1.000000000000000e+00_dp],'bandpass', &
        [-1.032264618598808e-02_dp,-7.111385831311694e-02_dp, &
            -8.898799975735899e-03_dp,3.987832771227273e-02_dp,7.086214909093243e-02_dp,1.590333344425934e-02_dp, &
            -1.266677702877821e-01_dp,-2.823815497604518e-01_dp,6.500541347171468e-01_dp,-2.823815497604518e-01_dp, &
            -1.266677702877821e-01_dp,1.590333344425934e-02_dp,7.086214909093243e-02_dp,3.987832771227273e-02_dp, &
            -8.898799975735899e-03_dp,-7.111385831311694e-02_dp,-1.032264618598808e-02_dp],3.0e-4_dp)
    ! Reference design 5: hilbert, 11 taps.
    call check_reference(11,[1.000000000000000e-01_dp,9.000000000000000e-01_dp],[1.000000000000000e+00_dp], &
            [1.000000000000000e+00_dp],'hilbert', &
        [1.075000170973511e-01_dp,8.603991924916471e-17_dp, &
            1.830134183305607e-01_dp,-6.139167567625127e-17_dp,6.263879094245898e-01_dp,0.000000000000000e+00_dp, &
            -6.263879094245898e-01_dp,6.139167567625127e-17_dp,-1.830134183305607e-01_dp,-8.603991924916471e-17_dp, &
            -1.075000170973511e-01_dp],3.0e-4_dp)
    ! Reference design 6: hilbert, 12 taps.
    call check_reference(12,[1.000000000000000e-01_dp,9.000000000000000e-01_dp],[1.000000000000000e+00_dp], &
            [1.000000000000000e+00_dp],'hilbert', &
        [7.360270075141666e-02_dp,3.645292963791466e-02_dp, &
            9.033682840713524e-02_dp,9.857842633044739e-02_dp,2.224714415649244e-01_dp,6.161900431661770e-01_dp, &
            -6.161900431661770e-01_dp,-2.224714415649244e-01_dp,-9.857842633044739e-02_dp,-9.033682840713524e-02_dp, &
            -3.645292963791466e-02_dp,-7.360270075141666e-02_dp],3.0e-4_dp)
    ! Reference design 7: differentiator, 11 taps.
    call check_reference(11,[1.000000000000000e-01_dp,9.000000000000000e-01_dp],[1.000000000000000e+00_dp], &
            [1.000000000000000e+00_dp],'differentiator', &
        [2.290690021866995e-02_dp, &
            -3.868300899492869e-02_dp,4.723955497288590e-02_dp,-7.471793694391624e-02_dp,1.567608294587378e-01_dp, &
            0.000000000000000e+00_dp,-1.567608294587378e-01_dp,7.471793694391624e-02_dp,-4.723955497288590e-02_dp, &
            3.868300899492869e-02_dp,-2.290690021866995e-02_dp],3.0e-4_dp)
    ! Reference design 8: differentiator, 12 taps.
    call check_reference(12,[1.000000000000000e-01_dp,9.000000000000000e-01_dp],[1.000000000000000e+00_dp], &
            [1.000000000000000e+00_dp],'differentiator', &
        [-6.312547847440328e-04_dp, &
            1.534984112616448e-03_dp,-3.029591418719435e-03_dp,6.927567000769218e-03_dp,-2.130791408581981e-02_dp, &
            2.014112191076097e-01_dp,-2.014112191076097e-01_dp,2.130791408581981e-02_dp,-6.927567000769218e-03_dp, &
            3.029591418719435e-03_dp,-1.534984112616448e-03_dp,6.312547847440328e-04_dp],3.0e-4_dp)

    ! High-order designs that previously failed in the double-precision
    ! exchange. The bounds below are independent dense-grid SciPy 1.17
    ! ripple references, with modest allowance for different exchange grids.
    ! Compilers without an extended real kind retain the dp fallback.
    if (selected_real_kind(30,400) > 0) then
        call check_hard_case(86, [55.135560157134343_dp, 332.85949520688717_dp], &
            0.90131079005307935_dp, 1.9460550103868526_dp, 'differentiator', 800.0_dp, 2.0e-12_dp)
        call check_hard_case(93, [1.8267911366373659_dp, 20.261063050362562_dp], &
            0.98787784438475501_dp, 1.7266745840927831_dp, 'hilbert', 48.0_dp, 5.2e-6_dp)
        call check_hard_case(82, [68.220056243703809_dp, 360.65709662922114_dp], &
            1.0360208487266178_dp, 1.1862817144586622_dp, 'differentiator', 800.0_dp, 3.5e-9_dp)
        call check_hard_case(88, [3.519590482101552_dp, 22.975486288373880_dp], &
            1.1154409076018119_dp, 2.4481028756498815_dp, 'differentiator', 48.0_dp, 9.4e-6_dp)
        call check_hard_case(85, [3.7619099997082794_dp, 22.266802637592082_dp], &
            1.0584257568289865_dp, 2.3415720176257500_dp, 'differentiator', 48.0_dp, 4.8e-5_dp)
        call check_hard_case(91, [4.2182192875463773_dp, 21.428004287599016_dp], &
            1.1380849911176247_dp, 1.9806073112159517_dp, 'differentiator', 48.0_dp, 9.5e-8_dp)
        call check_hard_case(90, [0.082540500105575482_dp, 0.44151340122423904_dp], &
            0.75104218149354518_dp, 1.1451377665575944_dp, 'hilbert', 1.0_dp, 8.5e-12_dp)
        call check_hard_case(88, [61.603190639171544_dp, 372.16765104445665_dp], &
            0.91547729638725761_dp, 1.5440788621504540_dp, 'differentiator', 800.0_dp, 9.5e-8_dp)
    end if

    ! Frequency units do not change the impulse response when band edges and fs
    ! are scaled by the same factor.
    taps = remez(17,[0.0_dp,6.0_dp,8.0_dp,20.0_dp], &
        [1.0_dp,0.0_dp],fs=40.0_dp,status=st,iterations=niter,ripple=ripple)
    if (st /= remez_status_success .or. size(taps) /= 17) error stop 'remez scaled fs'
    if (niter <= 0 .or. ripple <= 0.0_dp) error stop 'remez diagnostics'

    ! SciPy-like invalid input conventions; no undefined coefficients on errors.
    taps = remez(0,[0.0_dp,1.0_dp],[1.0_dp],status=st)
    call check_invalid(taps,st)
    taps = remez(11,[0.0_dp,0.4_dp,0.3_dp,1.0_dp],[1.0_dp,0.0_dp],status=st)
    call check_invalid(taps,st)
    taps = remez(11,[0.0_dp,0.5_dp],[1.0_dp,0.0_dp],status=st)
    call check_invalid(taps,st)
    taps = remez(11,[0.0_dp,0.5_dp],[1.0_dp],weight=[0.0_dp],status=st)
    call check_invalid(taps,st)
    taps = remez(11,[0.0_dp,0.5_dp],[1.0_dp],weight=[1.0_dp,2.0_dp],status=st)
    call check_invalid(taps,st)
    taps = remez(11,[0.0_dp,0.5_dp],[1.0_dp],fs=-1.0_dp,status=st)
    call check_invalid(taps,st)
    taps = remez(11,[0.0_dp,0.5_dp],[1.0_dp],filter_type='other',status=st)
    call check_invalid(taps,st)
    taps = remez(11,[0.0_dp,0.5_dp],[1.0_dp],grid_density=1,status=st)
    call check_invalid(taps,st)
    taps = remez(11,[0.0_dp,0.5_dp],[1.0_dp],maxiter=0,status=st)
    call check_invalid(taps,st)

    ! Failure to reach an equiripple design is reported, not silently accepted.
    taps = remez(31,[0.0_dp,0.3_dp,0.4_dp,1.0_dp],[1.0_dp,0.0_dp], &
        fs=2.0_dp,maxiter=1,status=st)
    if (st /= remez_status_no_convergence .or. size(taps) /= 0) &
        error stop 'expected Remez nonconvergence'
    print '(a)', 'PASS: Remez equiripple FIR designs'

contains
    subroutine check_reference(n,bands,desired,weight,kind,expected,tol)
        integer, intent(in) :: n !! number of taps
        real(dp), intent(in) :: bands(:) !! band endpoints
        real(dp), intent(in) :: desired(:) !! target response by band
        real(dp), intent(in) :: weight(:) !! error weights by band
        character(*), intent(in) :: kind !! filter type
        real(dp), intent(in) :: expected(:) !! independent SciPy reference taps
        real(dp), intent(in) :: tol !! tolerance for finite-grid extrema differences
        real(dp), allocatable :: h(:)
        integer :: status, iterations
        real(dp) :: ripple

        h = remez(n,bands,desired,weight=weight,filter_type=kind,fs=2.0_dp, &
            status=status,iterations=iterations,ripple=ripple)
        if (status /= remez_status_success .or. size(h) /= n) then
            print *, 'Remez failure ', n, kind, status, iterations
            error stop 'Remez reference case did not converge'
        end if
        if (.not. all(ieee_is_finite(h)) .or. .not. ieee_is_finite(ripple)) &
            error stop 'Remez returned nonfinite values'
        if (maxval(abs(h-expected)) > tol) then
            print *, 'Remez ',n,kind,' max coefficient error ',maxval(abs(h-expected))
            error stop 'Remez differs from SciPy reference'
        end if
        if (kind == 'bandpass') then
            if (maxval(abs(h-h(n:1:-1))) > 1.0e-13_dp) error stop 'Remez symmetry'
        else
            if (maxval(abs(h+h(n:1:-1))) > 1.0e-13_dp) error stop 'Remez antisymmetry'
        end if
        if (iterations < 1 .or. ripple <= 0.0_dp) error stop 'Remez diagnostics'
    end subroutine check_reference

    subroutine check_hard_case(n, bands, gain, weight, kind, rate, tolerance)
        integer, intent(in) :: n !! FIR tap count
        real(dp), intent(in) :: bands(2) !! positive, nonzero band limits
        real(dp), intent(in) :: gain !! desired Hilbert gain or differentiator slope
        real(dp), intent(in) :: weight !! positive error weight
        character(*), intent(in) :: kind !! Hilbert or differentiator
        real(dp), intent(in) :: rate !! sampling frequency
        real(dp), intent(in) :: tolerance !! maximum weighted in-band error
        integer, parameter :: nfreq = 1025
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), allocatable :: h(:)
        real(dp) :: f, omega, desired, scale, error_max
        complex(dp) :: response
        integer :: st, j, k, niter

        h = remez(n, bands, [gain], weight=[weight], filter_type=kind, &
            fs=rate, maxiter=30, status=st, iterations=niter)
        if (st /= remez_status_success .or. size(h) /= n) then
            print *, 'High-order Remez failed:', n, kind, st, niter
            error stop 'high-order Remez nonconvergence'
        end if
        if (.not. all(ieee_is_finite(h))) error stop 'high-order Remez nonfinite taps'
        error_max = 0.0_dp
        do j = 1, nfreq
            f = bands(1) + (bands(2)-bands(1))*real(j-1,dp)/real(nfreq-1,dp)
            omega = 2.0_dp*pi*f/rate
            response = cmplx(0.0_dp,0.0_dp,kind=dp)
            do k = 1, n
                response = response + h(k)*exp(cmplx(0.0_dp, &
                    (0.5_dp*real(n-1,dp)-real(k-1,dp))*omega, kind=dp))
            end do
            desired = gain
            scale = weight
            if (kind == 'differentiator') then
                desired = gain*f/rate
                scale = weight*rate/f
            end if
            error_max = max(error_max,abs(aimag(response)-desired)*scale)
        end do
        if (error_max > tolerance) then
            print *, 'High-order Remez ripple:', n, kind, error_max, tolerance
            error stop 'high-order Remez response mismatch'
        end if
    end subroutine check_hard_case

    subroutine check_invalid(h,status)
        real(dp), intent(in) :: h(:) !! expected empty array
        integer, intent(in) :: status !! expected invalid status
        if (size(h) /= 0 .or. status /= remez_status_invalid_input) &
            error stop 'Remez failed to reject invalid input'
    end subroutine check_invalid
end program test_signal_remez
