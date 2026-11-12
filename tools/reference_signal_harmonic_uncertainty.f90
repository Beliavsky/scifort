! SPDX-License-Identifier: MIT
! Text driver for NumPy comparisons of weighted harmonic parameter covariance.
program reference_signal_harmonic_uncertainty
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : harmonic_fit_result, &
        lombscargle_harmonic_fit, lombscargle_harmonic_predict_std
    implicit none

    integer :: ncases, casenum, n, m, nnew, nterms, center, floating, weighted
    integer :: k, j, h, col, row, status, p
    real(dp), allocatable :: x(:), y(:), freqs(:), w(:), xnew(:), se(:)
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
                precenter=center == 1, weights=w, floating_mean=floating == 1, &
                estimate_covariance=.true.)
        else
            fit = lombscargle_harmonic_fit(x,y,freqs,nterms, &
                precenter=center == 1, floating_mean=floating == 1, &
                estimate_covariance=.true.)
        end if
        write (*, '(i0,1x,i0)') fit%status,fit%degrees_of_freedom
        p = 1+2*nterms
        do k = 1, m
            write (*, '(i0)') merge(1,0,fit%valid(k))
            write (*, '(es26.17e3)') fit%residual_variance(k)
            do j = 1, p
                write (*, '(es26.17e3)') fit%standard_error(j,k)
            end do
            do col = 1, p
                do row = 1, p
                    write (*, '(es26.17e3)') fit%covariance(row,col,k)
                end do
            end do
            se = lombscargle_harmonic_predict_std(fit,xnew,k,status)
            write (*, '(i0)') status
            do h = 1, nnew
                write (*, '(es26.17e3)') se(h)
            end do
        end do
        deallocate(x,y,freqs,w,xnew)
    end do
end program reference_signal_harmonic_uncertainty
