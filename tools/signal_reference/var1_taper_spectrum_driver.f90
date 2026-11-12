! SPDX-License-Identifier: MIT
! Stdin: count, then n length nf, A(n,n), Sigma(n,n), taper, frequencies.
! Stdout: one line per scenario with status and, on success,
! Gamma(0) and (real,imag) cross-spectral matrices in Fortran column order.
program var1_taper_spectrum_driver
    use scifort_kinds, only : dp
    use scifort_signal_var1_taper_spectrum, only : var1_expected_taper_spectrum, &
        var1_taper_spectrum_result
    implicit none
    integer :: count, k, n, l, nf, i, j, f
    real(dp), allocatable :: a(:,:), sigma(:,:), taper(:), frequencies(:)
    type(var1_taper_spectrum_result) :: result

    read(*,*) count
    do k = 1, count
        read(*,*) n,l,nf
        allocate(a(n,n),sigma(n,n),taper(l),frequencies(nf))
        read(*,*) a
        read(*,*) sigma
        read(*,*) taper
        read(*,*) frequencies
        result = var1_expected_taper_spectrum(a,sigma,taper,frequencies)
        write(*,'(i0)',advance='no') result%status
        if (result%status == 0) then
            do j = 1, n
                do i = 1, n
                    write(*,'(1x,es25.16)',advance='no') result%lag0_covariance(i,j)
                end do
            end do
            do f = 1, nf
                do j = 1, n
                    do i = 1, n
                        write(*,'(2(1x,es25.16))',advance='no') &
                            real(result%power(f,i,j),dp),aimag(result%power(f,i,j))
                    end do
                end do
            end do
        end if
        print *
        deallocate(a,sigma,taper,frequencies)
    end do
end program var1_taper_spectrum_driver
