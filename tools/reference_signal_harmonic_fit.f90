! SPDX-License-Identifier: MIT
! Text reference driver for independently generated NumPy harmonic fits.
program reference_signal_harmonic_fit
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : lombscargle_harmonic_fit, &
        lombscargle_harmonic_predict, harmonic_fit_result
    implicit none
    integer :: ncases, casenum, n, m, nnew, nterms, center, floating, weighted
    integer :: k, j, h, status
    real(dp), allocatable :: x(:), y(:), freqs(:), w(:), xnew(:), prediction(:)
    type(harmonic_fit_result) :: fit

    read (*, *) ncases
    do casenum = 1, ncases
        read (*, *) n, m, nnew, nterms, center, floating, weighted
        allocate(x(n),y(n),freqs(m),w(n),xnew(nnew))
        read (*, *) x
        read (*, *) y
        read (*, *) freqs
        read (*, *) w
        read (*, *) xnew
        if (weighted == 1) then
            fit = lombscargle_harmonic_fit(x,y,freqs,nterms, &
                precenter=center == 1, weights=w, floating_mean=floating == 1)
        else
            fit = lombscargle_harmonic_fit(x,y,freqs,nterms, &
                precenter=center == 1, floating_mean=floating == 1)
        end if
        write (*, '(i0,1x,i0,1x,i0)') fit%status, nterms, m
        do k = 1, m
            write (*, '(i0)') merge(1,0,fit%valid(k))
            write (*, '(es26.17e3)') fit%offset(k)
            do h = 1, nterms
                write (*, '(es26.17e3)') fit%cosine(h,k)
                write (*, '(es26.17e3)') fit%sine(h,k)
            end do
            write (*, '(es26.17e3)') fit%power(k)
            write (*, '(es26.17e3)') fit%normalized_power(k)
            write (*, '(es26.17e3)') fit%rms(k)
            prediction = lombscargle_harmonic_predict(fit,xnew,k,status)
            write (*, '(i0)') status
            do j = 1, nnew
                write (*, '(es26.17e3)') prediction(j)
            end do
        end do
        deallocate(x,y,freqs,w,xnew)
    end do
end program reference_signal_harmonic_fit
