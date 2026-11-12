! SPDX-License-Identifier: MIT
program test_signal_canonical_block_jackknife
    use scifort_kinds, only : dp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result, multitaper_matrix_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, canonical_coherence_matrix
    use scifort_signal_canonical_jackknife, only : canonical_window_jackknife_result, &
        canonical_window_jackknife, canonical_block_jackknife, multitaper_canonical_block_jackknife
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra, bad
    type(multitaper_matrix_result) :: matrix
    type(canonical_window_jackknife_result) :: block, ordinary, singleton, band, direct
    type(canonical_coherence_result) :: point, leave
    complex(dp) :: factor(5,5), z(256,5)
    real(dp) :: x(256,5), roots(4), groups(4), bands(2,2), avg, expected, zvalues(4)
    integer :: f, t, i, j, k

    allocate(spectra%power(4,12,5,5),spectra%frequencies(4),spectra%times(12))
    spectra%frequencies = [0.0_dp,0.25_dp,0.5_dp,0.75_dp]
    spectra%times = [(real(t,dp),t=1,12)]
    spectra%status = 0
    do f = 1, 4
        do t = 1, 12
            factor = (0.0_dp,0.0_dp)
            do i = 1, 5
                factor(i,i) = cmplx(1.0_dp+0.2_dp*i+0.03_dp*t,0.0_dp,dp)
                do j = 1, i-1
                    factor(i,j) = cmplx(0.21_dp*cos(0.19_dp*i*j+0.18_dp*t), &
                        0.14_dp*sin(0.21_dp*j*f+0.10_dp*t),dp)
                end do
            end do
            spectra%power(f,t,:,:) = matmul(factor,transpose(conjg(factor)))*(0.8_dp+0.06_dp*f)
        end do
    end do

    block = canonical_block_jackknife(spectra,[1,2],[3,4],3,conditioned_on=[5], &
        ridge=0.04_dp,shrinkage=0.08_dp)
    call require(block%status == 0,'valid block result')
    call require(block%n_windows == 12 .and. block%n_blocks == 4 .and. block%block_length == 3, &
        'block metadata distinguishes windows from deletion groups')
    call require(all(block%coherence_valid),'aggregate block SE valid')
    call require(all(block%valid),'isolated modes valid')
    call require(all(block%stderr >= 0.0_dp),'nonnegative uncertainty')
    call require(all(block%lower95 <= block%squared_modes) .and. &
        all(block%upper95 >= block%squared_modes),'Fisher intervals contain center')

    allocate(matrix%frequencies(4),matrix%power(4,5,5))
    matrix%frequencies = spectra%frequencies
    matrix%power = sum(spectra%power,dim=2)/12.0_dp
    matrix%status = 0
    point = canonical_coherence_matrix(matrix,[1,2],[3,4],conditioned_on=[5], &
        ridge=0.04_dp,shrinkage=0.08_dp)
    call near(maxval(abs(point%squared_modes-block%squared_modes)),0.0_dp,1.0e-14_dp, &
        'full estimate independent of block length')
    do k = 1, 4
        matrix%power = (sum(spectra%power,dim=2)- &
            sum(spectra%power(:,(k-1)*3+1:k*3,:,:),dim=2))/9.0_dp
        leave = canonical_coherence_matrix(matrix,[1,2],[3,4],conditioned_on=[5], &
            ridge=0.04_dp,shrinkage=0.08_dp)
        roots(k) = leave%squared_modes(2,1)
        groups(k) = leave%coherence(2)
    end do
    avg = sum(roots)/4.0_dp
    expected = sqrt(3.0_dp/4.0_dp*sum((roots-avg)**2))
    call near(block%stderr(2,1),expected,2.0e-14_dp, &
        'delete-equal-block SE uses block count')
    call near(block%bias(2,1),3.0_dp*(avg-point%squared_modes(2,1)),2.0e-14_dp, &
        'delete-block bias uses block count')
    zvalues = atanh(sqrt(roots))
    avg = sum(zvalues)/4.0_dp
    expected = sqrt(3.0_dp/4.0_dp*sum((zvalues-avg)**2))
    call near(block%stderr_fisher(2,1),expected,2.0e-14_dp,'Fisher-scale block SE')
    avg = sum(groups)/4.0_dp
    expected = sqrt(3.0_dp/4.0_dp*sum((groups-avg)**2))
    call near(block%coherence_stderr(2),expected,2.0e-14_dp,'aggregate block SE')

    ordinary = canonical_window_jackknife(spectra,[1,2],[3,4],conditioned_on=[5], &
        ridge=0.04_dp,shrinkage=0.08_dp)
    singleton = canonical_block_jackknife(spectra,[1,2],[3,4],1,conditioned_on=[5], &
        ridge=0.04_dp,shrinkage=0.08_dp)
    call require(singleton%status == 0 .and. singleton%n_blocks == 12,'block length 1 alias')
    call near(maxval(abs(ordinary%stderr-singleton%stderr)),0.0_dp,1.0e-14_dp, &
        'one-window blocks reproduce old jackknife')
    call near(maxval(abs(ordinary%bias-singleton%bias)),0.0_dp,1.0e-14_dp, &
        'one-window blocks reproduce old bias')
    call near(maxval(abs(ordinary%stderr_fisher-singleton%stderr_fisher)),0.0_dp,1.0e-14_dp, &
        'one-window blocks reproduce Fisher SE')
    call require(all(ordinary%valid .eqv. singleton%valid),'one-window block validity')
    call require(ordinary%block_length == 1 .and. ordinary%n_blocks == 12,'legacy metadata defaults')

    bands(:,1) = [0.04_dp,0.37_dp]
    bands(:,2) = [0.42_dp,0.71_dp]
    band = canonical_block_jackknife(spectra,[1,2],[3,4],3,conditioned_on=[5], &
        ridge=0.04_dp,shrinkage=0.08_dp,bands=bands)
    call require(band%status == 0 .and. size(band%frequencies) == 2,'band block inference')
    call require(all(band%coherence_valid),'band aggregate SE')
    singleton = canonical_block_jackknife(spectra,[1,2],[3,4],3,conditioned_on=[5], &
        ridge=0.04_dp,shrinkage=0.08_dp,bands=bands,average=.true.)
    call near(maxval(abs(band%squared_modes-singleton%squared_modes)), &
        0.0_dp,1.0e-13_dp,'band integration versus band average')

    bad = spectra
    direct = canonical_block_jackknife(bad,[1],[2],0)
    call require(direct%status == 1,'reject zero block length')
    direct = canonical_block_jackknife(bad,[1],[2],5)
    call require(direct%status == 1,'reject unequal block sizes')
    direct = canonical_block_jackknife(bad,[1],[2],4)
    call require(direct%status == 1,'require at least four blocks')
    bad%power = (0.0_dp,0.0_dp)
    direct = canonical_block_jackknife(bad,[1],[2],3)
    call require(direct%status == 2 .and. .not. any(direct%coherence_valid), &
        'handle singular blocks')
    bad = spectra
    bad%power = (0.0_dp,0.0_dp)
    do t = 1, 12
        do f = 1, 4
            do i = 1, 5
                bad%power(f,t,i,i) = (1.0_dp,0.0_dp)
            end do
        end do
    end do
    direct = canonical_block_jackknife(bad,[1],[2],3)
    call require(direct%status == 0 .and. .not. any(direct%valid), &
        'zero modes should not receive intervals')
    call require(all(ieee_is_nan(direct%stderr)),'invalid mode SE NaN')

    do i = 1, 256
        do j = 1, 5
            x(i,j) = sin(0.095_dp*i*j)+0.23_dp*cos(0.16_dp*i+0.09_dp*j) &
                +0.04_dp*sin(0.17_dp*i*i+0.17_dp*j)
            z(i,j) = cmplx(x(i,j),0.16_dp*cos(0.08_dp*i*j+0.12_dp*j),dp)
        end do
    end do
    direct = multitaper_canonical_block_jackknife(x,3.0_dp,16,[1,2],[3,4],4, &
        conditioned_on=[5],kmax=5,weighting='eigen',ridge=0.10_dp,shrinkage=0.14_dp)
    call require(direct%status == 0 .and. direct%n_windows == 16 .and. direct%n_blocks == 4, &
        'direct real DPSS block jackknife')
    direct = multitaper_canonical_block_jackknife(z,3.0_dp,16,[1,2],[3,4],2, &
        conditioned_on=[5],kmax=5,weighting='adaptive',ridge=0.10_dp,shrinkage=0.14_dp, &
        bands=reshape([0.04_dp,0.18_dp,0.24_dp,0.42_dp],[2,2]))
    call require(direct%status == 0 .and. direct%n_blocks == 8 .and. &
        size(direct%frequencies) == 2,'direct complex adaptive band block jackknife')
    call require(all(direct%coherence_valid),'complex adaptive grouped estimates')
    print '(a)', 'PASS canonical block jackknife, analytic variance, bands and real/complex adaptive cases'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! assertion
        character(*), intent(in) :: message !! failed assertion
        if (.not. ok) then
            print '(a)', 'FAIL '//message
            error stop 1
        end if
    end subroutine require

    subroutine near(value,target,tolerance,message)
        real(dp), intent(in) :: value !! actual numerical result
        real(dp), intent(in) :: target !! independent expected value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        character(*), intent(in) :: message !! failed assertion
        if (.not. abs(value-target) <= tolerance) then
            print '(a,2es24.15)', 'FAIL '//message//' got/expected ',value,target
            error stop 1
        end if
    end subroutine near
end program test_signal_canonical_block_jackknife
