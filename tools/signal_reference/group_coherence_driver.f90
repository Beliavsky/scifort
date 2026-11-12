! SPDX-License-Identifier: MIT
! Text protocol driver for independent NumPy canonical-correlation comparisons.
program group_coherence_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_group_coherence, only : group_coherence_tf_result, group_coherence_tf_matrix
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(group_coherence_tf_result) :: result
    integer :: ncase, icase, nf, nt, nc, na, nb, nd, f, t, i, j
    integer, allocatable :: a(:), b(:), c(:)
    real(dp) :: ridge, shrinkage, re, im

    read(*,*) ncase
    do icase = 1, ncase
        read(*,*) nc,nf,nt,na,nb,nd,ridge,shrinkage
        allocate(a(na),b(nb),c(nd))
        read(*,*) a
        read(*,*) b
        if (nd > 0) read(*,*) c
        allocate(spectra%frequencies(nf),spectra%times(nt),spectra%power(nf,nt,nc,nc))
        spectra%frequencies = [(real(f,dp),f=1,nf)]
        spectra%times = [(real(t,dp),t=1,nt)]
        spectra%status = 0
        do t = 1, nt
            do f = 1, nf
                do j = 1, nc
                    do i = 1, nc
                        read(*,*) re,im
                        spectra%power(f,t,i,j) = cmplx(re,im,dp)
                    end do
                end do
            end do
        end do
        result = group_coherence_tf_matrix(spectra,a,b,conditioned_on=c,ridge=ridge,shrinkage=shrinkage)
        write(*,'(i0,1x,i0,1x,i0)') result%status,nf,nt
        do t = 1, nt
            do f = 1, nf
                write(*,'(es25.16e3,1x,es25.16e3,1x,l1)') &
                    result%coherence(f,t),result%log_information(f,t),result%valid(f,t)
            end do
        end do
        deallocate(a,b,c,spectra%frequencies,spectra%times,spectra%power)
    end do
end program group_coherence_driver
