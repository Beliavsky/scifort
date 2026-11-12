! SPDX-License-Identifier: MIT
program test_signal_band_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_multiple_coherence, only : multiple_coherence_result, &
        multiple_coherence_tf_result, multiple_coherence_matrix, multiple_coherence_tf_matrix
    use scifort_signal_partial_coherence, only : partial_coherence_result, &
        partial_coherence_tf_result, partial_coherence_matrix, partial_coherence_tf_matrix, &
        partial_coherence_success, partial_coherence_invalid_input
    use scifort_signal_band_coherence, only : integrate_spectral_bands, &
        integrate_tf_spectral_bands, band_multiple_coherence_matrix, &
        band_multiple_coherence_tf_matrix, band_partial_coherence_matrix, &
        band_partial_coherence_tf_matrix, multitaper_band_multiple_coherence, &
        multitaper_tf_band_multiple_coherence, multitaper_band_partial_coherence, &
        multitaper_tf_band_partial_coherence
    implicit none

    type(multitaper_matrix_result) :: s, integrated, averaged, bad, s_direct
    type(multitaper_spectral_matrix_spectrogram_result) :: st, bt, st_direct
    type(multiple_coherence_result) :: r, re, r_direct
    type(multiple_coherence_tf_result) :: rt, rt_direct
    type(partial_coherence_result) :: p, pe
    type(partial_coherence_tf_result) :: pt, pte
    real(dp) :: bands(2,3), x(128,3), expected, width, f, theta
    complex(dp) :: z(128,3), a(3,3), b(3,3), want(3,3)
    integer :: i, j, k, m, t

    ! A positive-definite Hermitian CSD linear in frequency is integrated
    ! *exactly*, even across FFT-wrapped negative frequencies and partial bins.
    allocate(s%frequencies(5),s%power(5,3,3),s%cv_shrinkage(5))
    s%status = 0
    s%effective_tapers = 5.0_dp
    s%cv_shrinkage = 0.2_dp
    s%frequencies = [0.0_dp,0.7_dp,1.5_dp,-1.8_dp,-0.4_dp]
    a = (0.0_dp,0.0_dp)
    a(1,1) = 2.0_dp
    a(2,2) = 1.5_dp
    a(3,3) = 1.2_dp
    a(1,2) = cmplx(0.4_dp,0.2_dp,dp)
    a(1,3) = cmplx(0.2_dp,-0.1_dp,dp)
    a(2,3) = cmplx(-0.3_dp,0.1_dp,dp)
    b = (0.0_dp,0.0_dp)
    b(1,1) = 0.15_dp
    b(2,2) = 0.1_dp
    b(3,3) = -0.12_dp
    b(1,2) = cmplx(0.025_dp,-0.01_dp,dp)
    b(1,3) = cmplx(-0.02_dp,0.015_dp,dp)
    b(2,3) = cmplx(0.01_dp,0.02_dp,dp)
    do j = 1, 3
        do i = j+1, 3
            a(i,j) = conjg(a(j,i))
            b(i,j) = conjg(b(j,i))
        end do
    end do
    do k = 1, 5
        s%power(k,:,:) = a+s%frequencies(k)*b
    end do
    bands(:,1) = [-1.4_dp,-0.1_dp]
    bands(:,2) = [-0.5_dp,1.1_dp]
    bands(:,3) = [-1.8_dp,1.5_dp]
    integrated = integrate_spectral_bands(s,bands)
    call require(integrated%status == 0,'static integration success')
    call require(.not. allocated(integrated%cv_shrinkage),'do not reuse frequency CV shrinkage')
    call near(integrated%effective_tapers,0.0_dp,0.0_dp,'do not reuse taper count')
    call require(size(integrated%frequencies) == 3,'three integrated bands')
    do m = 1, 3
        width = bands(2,m)-bands(1,m)
        f = 0.5_dp*(bands(1,m)+bands(2,m))
        want = width*(a+f*b)
        call near(maxval(abs(integrated%power(m,:,:)-want)),0.0_dp,1.0e-13_dp, &
            'exact integration of linearly interpolated CSD')
        call near(integrated%frequencies(m),f,1.0e-14_dp,'band center metadata')
    end do
    averaged = integrate_spectral_bands(s,bands,average=.true.)
    do m = 1, 3
        width = bands(2,m)-bands(1,m)
        call near(maxval(abs(averaged%power(m,:,:)-integrated%power(m,:,:)/width)), &
            0.0_dp,8.0e-15_dp,'mean vs integrated power')
    end do
    r = band_multiple_coherence_matrix(s,bands,ridge=0.03_dp,shrinkage=0.15_dp)
    re = multiple_coherence_matrix(integrated,ridge=0.03_dp,shrinkage=0.15_dp)
    p = band_partial_coherence_matrix(s,bands,ridge=0.03_dp,shrinkage=0.15_dp, &
        return_coherency=.true.)
    pe = partial_coherence_matrix(integrated,ridge=0.03_dp,shrinkage=0.15_dp, &
        return_coherency=.true.)
    call require(r%status == partial_coherence_success,'band multiple coherence succeeded')
    call require(p%status == partial_coherence_success,'band partial coherence succeeded')
    call near(maxval(abs(r%coherence-re%coherence)),0.0_dp,1.0e-13_dp, &
        'band multiple vs pre-integrated matrix')
    call near(maxval(abs(p%coherency-pe%coherency)),0.0_dp,1.0e-13_dp, &
        'band partial complex coherency vs pre-integrated matrix')
    call require(all(r%coherence >= 0.0_dp) .and. all(r%coherence <= 1.0_dp), &
        'multiple coherence stays in unit interval')
    ! Coherence of the integrated cross spectrum is generally NOT the
    ! arithmetic mean of coherence at different frequencies.
    call require(abs(r%coherence(1,1)-r%coherence(2,1)) > 1.0e-7_dp, &
        'frequency-dependent band multiple coherence')

    allocate(st%frequencies(5),st%times(3),st%power(5,3,3,3),st%cv_shrinkage(5,3))
    st%status = 0
    st%frequencies = s%frequencies
    st%times = [0.125_dp,0.25_dp,0.375_dp]
    st%effective_tapers = 6.0_dp
    st%cv_shrinkage = 0.1_dp
    do t = 1, 3
        st%power(:,t,:,:) = s%power*real(t,dp)
    end do
    bt = integrate_tf_spectral_bands(st,bands)
    call require(bt%status == 0,'TF band integration succeeded')
    call require(.not. allocated(bt%cv_shrinkage),'TF CV metadata discarded')
    call require(size(bt%power,2) == 3,'time axis preserved')
    do t = 1, 3
        call near(maxval(abs(bt%power(:,t,:,:)-integrated%power*real(t,dp))), &
            0.0_dp,3.0e-14_dp,'TF integrated band power')
    end do
    rt = band_multiple_coherence_tf_matrix(st,bands,ridge=0.015_dp,shrinkage=0.25_dp)
    rt_direct = multiple_coherence_tf_matrix(bt,ridge=0.015_dp,shrinkage=0.25_dp)
    pt = band_partial_coherence_tf_matrix(st,bands,ridge=0.015_dp,shrinkage=0.25_dp, &
        return_coherency=.true.)
    pte = partial_coherence_tf_matrix(bt,ridge=0.015_dp,shrinkage=0.25_dp, &
        return_coherency=.true.)
    call require(rt%status == 0 .and. pt%status == 0,'time-frequency band coherence')
    call near(maxval(abs(rt%coherence-rt_direct%coherence)),0.0_dp,1.0e-13_dp, &
        'TF multiple coherence equivalence')
    call near(maxval(abs(pt%coherency-pte%coherency)),0.0_dp,1.0e-13_dp, &
        'TF partial coherency equivalence')

    ! Invalid geometry/frequencies/CSD never produces deceptively valid bands.
    bad = integrate_spectral_bands(s,reshape([-2.0_dp,0.1_dp],[2,1]))
    call require(bad%status /= 0 .and. size(bad%power) == 0,'reject low out-of-range band')
    bad = integrate_spectral_bands(s,reshape([-0.1_dp,2.0_dp],[2,1]))
    call require(bad%status /= 0,'reject high out-of-range band')
    bad = integrate_spectral_bands(s,reshape([1.0_dp,0.5_dp],[2,1]))
    call require(bad%status /= 0,'reject reversed band')
    bad = integrate_spectral_bands(s,reshape([1.0_dp,1.0_dp],[2,1]))
    call require(bad%status /= 0,'reject zero-width band')
    bad = integrate_spectral_bands(s,reshape([0.1_dp,0.2_dp,0.3_dp],[3,1]))
    call require(bad%status /= 0,'reject invalid band shape')
    bad = integrate_spectral_bands(s,reshape([0.1_dp,0.2_dp],[2,0]))
    call require(bad%status /= 0,'reject zero bands')
    bad = integrate_spectral_bands(s,reshape([ieee_value(0.0_dp,ieee_quiet_nan),1.0_dp],[2,1]))
    call require(bad%status /= 0,'reject nonfinite band')
    s%frequencies(1) = s%frequencies(2)
    bad = integrate_spectral_bands(s,bands)
    call require(bad%status /= 0,'reject repeated frequency')
    s%frequencies(1) = 0.0_dp
    s%power(1,1,1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    bad = integrate_spectral_bands(s,bands)
    call require(bad%status /= 0,'reject nonfinite spectrum')
    s%power(1,1,1) = a(1,1)
    s%status = 1
    bad = integrate_spectral_bands(s,bands)
    call require(bad%status /= 0,'reject failed upstream spectral estimate')
    r = band_multiple_coherence_matrix(s,bands)
    call require(r%status == partial_coherence_invalid_input,'reject bad matrix in band coherence')
    s%status = 0
    st%times(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    bt = integrate_tf_spectral_bands(st,bands)
    call require(bt%status /= 0,'reject nonfinite segment time')
    st%times(1) = 0.125_dp
    bt = integrate_tf_spectral_bands(st,reshape([1.0_dp,0.2_dp],[2,1]))
    call require(bt%status /= 0,'reject invalid TF bands')

    ! Signal-level convenience with real/complex inputs and wrapped two-sided FFT.
    do k = 1, 128
        theta = real(k,dp)
        x(k,1) = sin(0.07_dp*theta)+0.2_dp*cos(0.31_dp*theta)
        x(k,2) = 0.5_dp*x(k,1)+sin(0.23_dp*theta)
        x(k,3) = 0.4_dp*x(k,2)+0.3_dp*cos(0.19_dp*theta)
        do j = 1, 3
            z(k,j) = cmplx(x(k,j),0.3_dp*cos(0.043_dp*theta*real(j,dp)),dp)
        end do
    end do
    block
        real(dp) :: positive(2,1), negative(2,2)
        positive(:,1) = [0.08_dp,0.35_dp]
        negative(:,1) = [-0.30_dp,-0.05_dp]
        negative(:,2) = [-0.1_dp,0.35_dp]
        s_direct = multitaper_spectral_matrix(x,3.5_dp,kmax=6,nfft=128,weighting='eigen')
        r = multitaper_band_multiple_coherence(x,3.5_dp,positive, &
            kmax=6,nfft=128,weighting='eigen',ridge=0.02_dp)
        re = band_multiple_coherence_matrix(s_direct,positive,ridge=0.02_dp)
        call require(r%status == 0,'real band multitaper result')
        call near(maxval(abs(r%coherence-re%coherence)),0.0_dp,5.0e-13_dp, &
            'direct real band vs spectral matrix')
        s_direct = multitaper_spectral_matrix(z,3.5_dp,kmax=6,nfft=128,weighting='adaptive')
        r = multitaper_band_multiple_coherence(z,3.5_dp,negative, &
            kmax=6,nfft=128,weighting='adaptive',ridge=0.03_dp)
        re = band_multiple_coherence_matrix(s_direct,negative,ridge=0.03_dp)
        call require(r%status == 0,'complex adaptive band multitaper result')
        call near(maxval(abs(r%coherence-re%coherence)),0.0_dp,5.0e-13_dp, &
            'direct complex adaptive band vs spectral matrix')
        st_direct = multitaper_spectral_matrix_spectrogram(z,3.5_dp,64, &
            noverlap=32,kmax=5,nfft=64,weighting='equal')
        rt = multitaper_tf_band_multiple_coherence(z,3.5_dp,64,negative, &
            noverlap=32,kmax=5,nfft=64,weighting='equal',ridge=0.02_dp)
        rt_direct = band_multiple_coherence_tf_matrix(st_direct,negative,ridge=0.02_dp)
        call require(rt%status == 0,'complex TF band multitaper result')
        call near(maxval(abs(rt%coherence-rt_direct%coherence)),0.0_dp,5.0e-13_dp, &
            'direct complex TF band vs spectral matrix')
        st_direct = multitaper_spectral_matrix_spectrogram(x,3.5_dp,64, &
            noverlap=32,kmax=5,nfft=64,weighting='eigen')
        rt = multitaper_tf_band_multiple_coherence(x,3.5_dp,64,positive, &
            noverlap=32,kmax=5,nfft=64,weighting='eigen',ridge=0.02_dp)
        rt_direct = band_multiple_coherence_tf_matrix(st_direct,positive,ridge=0.02_dp)
        call require(rt%status == 0,'real TF band multitaper result')
        p = multitaper_band_partial_coherence(x,3.5_dp,positive, &
            kmax=6,nfft=128,weighting='eigen',ridge=0.02_dp,return_coherency=.true.)
        pe = band_partial_coherence_matrix(multitaper_spectral_matrix(x,3.5_dp, &
            kmax=6,nfft=128,weighting='eigen'),positive, &
            ridge=0.02_dp,return_coherency=.true.)
        call require(p%status == 0,'real direct band partial coherence')
        call near(maxval(abs(p%coherency-pe%coherency)),0.0_dp,5.0e-13_dp, &
            'real direct band partial phase')
        p = multitaper_band_partial_coherence(z,3.5_dp,negative, &
            kmax=6,nfft=128,weighting='adaptive',ridge=0.02_dp)
        call require(p%status == 0,'complex direct adaptive band partial coherence')
        pt = multitaper_tf_band_partial_coherence(z,3.5_dp,64,negative, &
            noverlap=32,kmax=5,nfft=64,weighting='equal',ridge=0.02_dp, &
            return_coherency=.true.)
        pte = band_partial_coherence_tf_matrix(multitaper_spectral_matrix_spectrogram( &
            z,3.5_dp,64,noverlap=32,kmax=5,nfft=64,weighting='equal'), &
            negative,ridge=0.02_dp,return_coherency=.true.)
        call require(pt%status == 0,'complex direct TF band partial coherence')
        call near(maxval(abs(pt%coherency-pte%coherency)),0.0_dp,5.0e-13_dp, &
            'complex direct time-frequency conditional phase')
        pt = multitaper_tf_band_partial_coherence(x,3.5_dp,64,positive, &
            noverlap=32,kmax=5,nfft=64,weighting='eigen',ridge=0.02_dp)
        call require(pt%status == 0,'real direct TF band partial coherence')
        call near(maxval(abs(rt%coherence-rt_direct%coherence)),0.0_dp,5.0e-13_dp, &
            'direct real TF band vs spectral matrix')
    end block

    print *, 'PASS spectral band integration / partial / multiple coherence regressions'

contains
    subroutine require(condition,description)
        logical, intent(in) :: condition !! test predicate
        character(*), intent(in) :: description !! failure label
        if (.not. condition) then
            print *, 'FAIL: ',description
            error stop 1
        end if
    end subroutine require

    subroutine near(actual,reference,tolerance,description)
        real(dp), intent(in) :: actual !! measured value
        real(dp), intent(in) :: reference !! target
        real(dp), intent(in) :: tolerance !! acceptable absolute error
        character(*), intent(in) :: description !! failure label
        call require(abs(actual-reference) <= tolerance,description)
    end subroutine near
end program test_signal_band_coherence
