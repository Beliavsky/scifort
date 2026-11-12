! SPDX-License-Identifier: MIT
! Stdin: cases; p,length,nf,nb; A(5,5,p), Sigma(5,5), taper,
! frequencies, bands(2,nb). Groups A=[1,2], B=[3,4], conditioned=[5].
! Stdout: status of each, then (coherence,mode1,mode2) per point and band.
program varp_canonical_target_driver
    use scifort_kinds, only : dp
    use scifort_signal_varp_canonical_target, only : varp_expected_canonical_coherence, &
        varp_expected_band_canonical_coherence
    use scifort_signal_canonical_coherence, only : canonical_coherence_result
    implicit none
    integer :: count, k, p, l, nf, nb, f
    real(dp), allocatable :: a(:,:,:), sigma(:,:), taper(:), freq(:), bands(:,:)
    type(canonical_coherence_result) :: point, band

    read(*,*) count
    do k = 1, count
        read(*,*) p,l,nf,nb
        allocate(a(5,5,p),sigma(5,5),taper(l),freq(nf),bands(2,nb))
        read(*,*) a
        read(*,*) sigma
        read(*,*) taper
        read(*,*) freq
        read(*,*) bands
        point = varp_expected_canonical_coherence(a,sigma,taper,freq,[1,2],[3,4], &
            conditioned_on=[5])
        band = varp_expected_band_canonical_coherence(a,sigma,taper,bands,[1,2],[3,4], &
            conditioned_on=[5])
        write(*,'(2(i0,1x))',advance='no') point%status,band%status
        if (point%status == 0 .and. band%status == 0) then
            do f = 1, nf
                write(*,'(3(es25.16,1x))',advance='no') &
                    point%coherence(f),point%squared_modes(f,1),point%squared_modes(f,2)
            end do
            do f = 1, nb
                write(*,'(3(es25.16,1x))',advance='no') &
                    band%coherence(f),band%squared_modes(f,1),band%squared_modes(f,2)
            end do
        end if
        print *
        deallocate(a,sigma,taper,freq,bands)
    end do
end program varp_canonical_target_driver
