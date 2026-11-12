! SPDX-License-Identifier: MIT
program test_signal_canonical_jackknife
    use scifort_kinds, only : dp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result, multitaper_matrix_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, canonical_coherence_matrix
    use scifort_signal_canonical_jackknife, only : canonical_window_jackknife, &
        multitaper_canonical_window_jackknife, canonical_window_jackknife_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra, broken
    type(multitaper_matrix_result) :: matrix
    type(canonical_window_jackknife_result) :: uncertainty, band, direct
    type(canonical_coherence_result) :: point, leave
    complex(dp) :: tri(5,5), cov(5,5), z(192,5)
    real(dp) :: x(192,5), values(8), avg, expected, bands(2,2), phase
    integer :: f, t, i, j

    allocate(spectra%power(4,8,5,5),spectra%frequencies(4),spectra%times(8))
    spectra%frequencies = [0.0_dp,0.5_dp,1.0_dp,1.5_dp]
    spectra%times = [(real(t,dp),t=1,8)]
    spectra%status = 0
    do f = 1, 4
        do t = 1, 8
            tri = (0.0_dp,0.0_dp)
            do i = 1, 5
                tri(i,i) = cmplx(1.0_dp+0.15_dp*i+0.03_dp*t,0.0_dp,dp)
                do j = 1, i-1
                    tri(i,j) = cmplx(0.17_dp*cos(0.2_dp*i*j+0.3_dp*t), &
                        0.12_dp*sin(0.3_dp*j*f+0.1_dp*t),dp)
                end do
            end do
            cov = matmul(tri,transpose(conjg(tri)))
            spectra%power(f,t,:,:) = cov*(0.8_dp+0.06_dp*f)
        end do
    end do
    uncertainty = canonical_window_jackknife(spectra,[1,2],[3,4], &
        conditioned_on=[5],ridge=0.02_dp,shrinkage=0.08_dp)
    call require(uncertainty%status == 0,'valid window jackknife')
    call require(uncertainty%n_windows == 8,'eight windows')
    call require(all(shape(uncertainty%squared_modes)==[4,2]),'shape')
    call require(all(uncertainty%valid),'two isolated roots')
    call require(all(uncertainty%coherence_valid),'group aggregate SE')
    call require(all(uncertainty%stderr >= 0.0_dp),'nonnegative SE')
    call require(all(uncertainty%lower95 <= uncertainty%squared_modes) .and. &
        all(uncertainty%upper95 >= uncertainty%squared_modes),'bounds cover estimate')
    allocate(matrix%frequencies(4),matrix%power(4,5,5))
    matrix%frequencies = spectra%frequencies
    matrix%power = sum(spectra%power,dim=2)/8.0_dp
    matrix%status = 0
    point = canonical_coherence_matrix(matrix,[1,2],[3,4],conditioned_on=[5], &
        ridge=0.02_dp,shrinkage=0.08_dp)
    call near(maxval(abs(point%squared_modes-uncertainty%squared_modes)),0.0_dp,1.0e-14_dp, &
        'jackknife uses correct full estimator')
    do t = 1, 8
        matrix%power = (sum(spectra%power,dim=2)-spectra%power(:,t,:,:))/7.0_dp
        leave = canonical_coherence_matrix(matrix,[1,2],[3,4],conditioned_on=[5], &
            ridge=0.02_dp,shrinkage=0.08_dp)
        values(t) = leave%squared_modes(2,1)
    end do
    avg = sum(values)/8.0_dp
    expected = sqrt(7.0_dp/8.0_dp*sum((values-avg)**2))
    call near(uncertainty%stderr(2,1),expected,1.0e-14_dp,'delete-one jackknife variance')
    call near(uncertainty%bias(2,1),7.0_dp*(avg-point%squared_modes(2,1)), &
        1.0e-14_dp,'delete-one bias')

    bands(:,1) = [0.1_dp,0.7_dp]
    bands(:,2) = [0.9_dp,1.4_dp]
    band = canonical_window_jackknife(spectra,[1,2],[3,4],conditioned_on=[5], &
        bands=bands,ridge=0.03_dp,shrinkage=0.1_dp)
    call require(band%status == 0 .and. all(band%coherence_valid),'band windows')
    call require(all(shape(band%squared_modes)==[2,2]),'integrated band shape')
    call require(all(abs(band%frequencies-[0.4_dp,1.15_dp]) < 1.0e-14_dp),'band centers')
    block
        type(canonical_window_jackknife_result) :: density
        density = canonical_window_jackknife(spectra,[1,2],[3,4],conditioned_on=[5], &
            bands=bands,average=.true.,ridge=0.03_dp,shrinkage=0.1_dp)
        call near(maxval(abs(density%squared_modes-band%squared_modes)),0.0_dp,1.0e-13_dp, &
            'band density and integral have same scale-invariant roots')
    end block

    broken = spectra
    broken%times(4) = broken%times(3)
    direct = canonical_window_jackknife(broken,[1],[2])
    call require(direct%status == 1,'reject nonmonotone time stamps')
    broken = spectra
    broken%power = (0.0_dp,0.0_dp)
    direct = canonical_window_jackknife(broken,[1],[2],ridge=0.0_dp)
    call require(direct%status == 2 .and. .not. any(direct%coherence_valid), &
        'singular windows yield missing uncertainties')
    direct = canonical_window_jackknife(spectra,[1,1],[2])
    call require(direct%status == 1,'reject overlapping channel indices')
    direct = canonical_window_jackknife(spectra,[1],[2],gap_tolerance=0.0_dp)
    call require(direct%status == 1,'reject invalid tolerance')
    direct = canonical_window_jackknife(spectra,[1],[2],bands=reshape([0.1_dp,4.0_dp],[2,1]))
    call require(direct%status == 1,'reject out-of-range bands')
    broken = spectra
    broken%power = (0.0_dp,0.0_dp)
    do t = 1, 8
        do f = 1, 4
            do i = 1, 5
                broken%power(f,t,i,i) = (1.0_dp,0.0_dp)
            end do
        end do
    end do
    direct = canonical_window_jackknife(broken,[1],[2])
    call require(direct%status == 0 .and. .not. any(direct%valid), &
        'zero canonical roots have no regular jackknife interval')
    call require(all(ieee_is_nan(direct%stderr)),'invalid mode SE is NaN')
    ! Repeated canonical eigenvalues identify only a subspace, not separate
    ! reliably ordered roots. They have valid aggregate coherence but no
    ! ordinary rank-specific jackknife uncertainty.
    broken = spectra
    broken%power = (0.0_dp,0.0_dp)
    do t = 1, 8
        do f = 1, 4
            do i = 1, 5
                broken%power(f,t,i,i) = (1.0_dp,0.0_dp)
            end do
            broken%power(f,t,1,3) = (0.7_dp,0.0_dp)
            broken%power(f,t,3,1) = (0.7_dp,0.0_dp)
            broken%power(f,t,2,4) = (0.7_dp,0.0_dp)
            broken%power(f,t,4,2) = (0.7_dp,0.0_dp)
        end do
    end do
    direct = canonical_window_jackknife(broken,[1,2],[3,4])
    call require(direct%status == 0 .and. all(direct%coherence_valid), &
        'repeated roots preserve aggregate inference')
    call require(.not. any(direct%valid),'repeated roots suppress mode-specific SE')
    broken%power(:,:,1,3) = (0.9999999_dp,0.0_dp)
    broken%power(:,:,3,1) = (0.9999999_dp,0.0_dp)
    broken%power(:,:,2,4) = (0.0_dp,0.0_dp)
    broken%power(:,:,4,2) = (0.0_dp,0.0_dp)
    direct = canonical_window_jackknife(broken,[1,2],[3,4])
    call require(direct%status == 0 .and. .not. any(direct%valid), &
        'boundary roots suppress uncalibrated jackknife intervals')

    do t = 1, 192
        do i = 1, 5
            phase = 0.13_dp*real(t*i,dp)
            x(t,i) = sin(phase)+0.4_dp*cos(0.08_dp*t+0.23_dp*i)+ &
                0.05_dp*sin(0.21_dp*t*t+0.17_dp*i)
            z(t,i) = cmplx(x(t,i),0.3_dp*cos(0.06_dp*t*i+0.12_dp*i),dp)
        end do
    end do
    direct = multitaper_canonical_window_jackknife(x,3.0_dp,24,[1,2],[3,4], &
        conditioned_on=[5],kmax=5,weighting='eigen',ridge=0.08_dp,shrinkage=0.2_dp)
    call require(direct%status == 0 .and. direct%n_windows == 8,'direct real eigen multitaper')
    call require(all(direct%coherence_valid),'real multitaper complete fits')
    direct = multitaper_canonical_window_jackknife(z,3.0_dp,24,[1,2],[3,4], &
        conditioned_on=[5],kmax=5,weighting='adaptive',ridge=0.08_dp,shrinkage=0.2_dp, &
        bands=reshape([0.02_dp,0.12_dp,0.14_dp,0.30_dp],[2,2]))
    call require(direct%status == 0 .and. direct%n_windows == 8, &
        'direct complex adaptive multitaper band')
    call require(all(direct%coherence_valid),'complex adaptive windows valid')
    print '(a)', 'PASS canonical window jackknife and band, real/complex, adaptive tests'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! assertion
        character(*), intent(in) :: message !! failed assertion description
        if (.not. ok) then
            print '(a)', 'FAIL '//message
            error stop 1
        end if
    end subroutine require
    subroutine near(value,target,tolerance,message)
        real(dp), intent(in) :: value !! actual result
        real(dp), intent(in) :: target !! expected result
        real(dp), intent(in) :: tolerance !! accepted absolute difference
        character(*), intent(in) :: message !! description
        if (.not. abs(value-target) <= tolerance) then
            print '(a,2es24.15)', 'FAIL '//message//' got/expected ',value,target
            error stop 1
        end if
    end subroutine near
end program test_signal_canonical_jackknife
