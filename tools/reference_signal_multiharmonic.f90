! SPDX-License-Identifier: MIT
! Text-based reference driver for multi-harmonic Lomb-Scargle comparisons.
program reference_signal_multiharmonic
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : lombscargle_multiharmonic
    implicit none
    integer :: ncases, casenum, n, m, nterms, cent, floating, weighted, normalized, status, i
    real(dp), allocatable :: x(:), y(:), f(:), w(:), power(:)
    character(len=9) :: mode

    read (*, *) ncases
    do casenum = 1, ncases
        read (*, *) n, m, nterms, cent, floating, weighted, normalized
        allocate(x(n), y(n), f(m), w(n))
        read (*, *) x
        read (*, *) y
        read (*, *) f
        read (*, *) w
        mode = 'power'
        if (normalized == 1) mode = 'normalize'
        if (weighted == 1) then
            power = lombscargle_multiharmonic(x, y, f, nterms, &
                precenter=cent == 1, floating_mean=floating == 1, &
                normalize=trim(mode), weights=w, status=status)
        else
            power = lombscargle_multiharmonic(x, y, f, nterms, &
                precenter=cent == 1, floating_mean=floating == 1, &
                normalize=trim(mode), status=status)
        end if
        write (*, '(i0,1x,i0)') status, size(power)
        do i = 1, size(power)
            write (*, '(es26.17e3)') power(i)
        end do
        deallocate(x, y, f, w, power)
    end do
end program reference_signal_multiharmonic
