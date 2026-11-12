! SPDX-License-Identifier: MIT
! Stdin: number of cases, then n,p,length,ntaper,nfreq,nband;
! A, Sigma, taper bank, frequencies, bands, taper weights;
! normalized, average. Arrays are written in Fortran column-major order.
! Stdout: two status codes, full complex spectra (real,imag pairs), and
! conditional canonical (aggregate, squared mode 1, squared mode 2).
program varp_multitaper_target_driver
    use scifort_kinds, only : dp
    use scifort_signal_varp_taper_spectrum, only : varp_taper_spectrum_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_result
    use scifort_signal_varp_multitaper_target, only : varp_expected_multitaper_spectrum, &
        varp_expected_band_multitaper_spectrum, &
        varp_expected_multitaper_canonical_coherence, &
        varp_expected_band_multitaper_canonical_coherence
    implicit none
    real(dp), allocatable :: a(:,:,:), sigma(:,:), tapers(:,:), frequencies(:), &
        bands(:,:), weights(:)
    type(varp_taper_spectrum_result) :: point, integrated
    type(canonical_coherence_result) :: canonical_point, canonical_band
    integer :: count, case_number, n, p, length, ntaper, nf, nb, norm, av
    integer :: i, j, f

    read(*,*) count
    do case_number = 1, count
        read(*,*) n,p,length,ntaper,nf,nb
        allocate(a(n,n,p),sigma(n,n),tapers(length,ntaper),frequencies(nf), &
            bands(2,nb),weights(ntaper))
        read(*,*) a
        read(*,*) sigma
        read(*,*) tapers
        read(*,*) frequencies
        read(*,*) bands
        read(*,*) weights
        read(*,*) norm,av
        point = varp_expected_multitaper_spectrum(a,sigma,tapers,frequencies, &
            taper_weights=weights,normalized=(norm==1))
        integrated = varp_expected_band_multitaper_spectrum(a,sigma,tapers,bands, &
            taper_weights=weights,normalized=(norm==1),average=(av==1))
        write(*,'(2(i0,1x))',advance='no') point%status,integrated%status
        if (point%status == 0 .and. integrated%status == 0) then
            canonical_point = varp_expected_multitaper_canonical_coherence( &
                a,sigma,tapers,frequencies,[1,2],[3,4],conditioned_on=[5], &
                taper_weights=weights,normalized=(norm==1))
            canonical_band = varp_expected_band_multitaper_canonical_coherence( &
                a,sigma,tapers,bands,[1,2],[3,4],conditioned_on=[5], &
                taper_weights=weights,normalized=(norm==1),average=(av==1))
            do f = 1, nf
                do j = 1, n
                    do i = 1, n
                        write(*,'(2(es25.16,1x))',advance='no') &
                            real(point%power(f,i,j),dp),aimag(point%power(f,i,j))
                    end do
                end do
            end do
            do f = 1, nb
                do j = 1, n
                    do i = 1, n
                        write(*,'(2(es25.16,1x))',advance='no') &
                            real(integrated%power(f,i,j),dp),aimag(integrated%power(f,i,j))
                    end do
                end do
            end do
            do f = 1, nf
                write(*,'(3(es25.16,1x))',advance='no') &
                    canonical_point%coherence(f),canonical_point%squared_modes(f,:)
            end do
            do f = 1, nb
                write(*,'(3(es25.16,1x))',advance='no') &
                    canonical_band%coherence(f),canonical_band%squared_modes(f,:)
            end do
        end if
        print *
        deallocate(a,sigma,tapers,frequencies,bands,weights)
    end do
end program varp_multitaper_target_driver
