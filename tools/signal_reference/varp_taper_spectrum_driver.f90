! SPDX-License-Identifier: MIT
! Stdin: cases; n,p,length,nf,nb; A; Sigma; taper; frequencies;
! bands; average (0/1), normalized (0/1). All multidimensional arrays
! are in Fortran column-major order.
! Stdout: status, stationary lag-zero covariance, then complex frequency
! spectral matrices and exact band spectral matrices as (real,imag) pairs.
program varp_taper_spectrum_driver
    use scifort_kinds, only : dp
    use scifort_signal_varp_taper_spectrum, only : varp_expected_taper_spectrum, &
        varp_expected_band_taper_spectrum, varp_taper_spectrum_result
    implicit none
    integer :: nc, k, n, p, l, nf, nb, i, j, f, av, norm
    real(dp), allocatable :: a(:,:,:), sigma(:,:), taper(:), freq(:), bands(:,:)
    type(varp_taper_spectrum_result) :: point, integrated

    read(*,*) nc
    do k = 1, nc
        read(*,*) n,p,l,nf,nb
        allocate(a(n,n,p),sigma(n,n),taper(l),freq(nf),bands(2,nb))
        read(*,*) a
        read(*,*) sigma
        read(*,*) taper
        read(*,*) freq
        read(*,*) bands
        read(*,*) av,norm
        point = varp_expected_taper_spectrum(a,sigma,taper,freq,normalized=(norm==1))
        integrated = varp_expected_band_taper_spectrum(a,sigma,taper,bands, &
            average=(av==1),normalized=(norm==1))
        write(*,'(2(i0,1x))',advance='no') point%status,integrated%status
        if (point%status == 0 .and. integrated%status == 0) then
            do j = 1, n
                do i = 1, n
                    write(*,'(es25.16,1x)',advance='no') point%lag0_covariance(i,j)
                end do
            end do
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
        end if
        print *
        deallocate(a,sigma,taper,freq,bands)
    end do
end program varp_taper_spectrum_driver
