! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Standalone batch interface for the optional scipy.signal.remez comparison.
program remez_driver
    use scifort_kinds, only : dp
    use scifort_signal_remez, only : remez
    implicit none
    integer :: ncases, n, nb, ftype, gd, st, it, k
    real(dp) :: rate, ripple
    real(dp), allocatable :: h(:), bands(:), desired(:), weight(:)
    character(16) :: kind

    read(*,*) ncases
    do k = 1, ncases
        read(*,*) n, nb, ftype, gd, rate
        allocate(bands(2*nb),desired(nb),weight(nb))
        read(*,*) bands
        read(*,*) desired
        read(*,*) weight
        select case (ftype)
        case (1)
            kind = 'bandpass'
        case (2)
            kind = 'hilbert'
        case (3)
            kind = 'differentiator'
        end select
        h = remez(n,bands,desired,weight=weight,filter_type=kind,fs=rate, &
            grid_density=gd,maxiter=30,status=st,iterations=it,ripple=ripple)
        print '(3(i0,1x),es25.16)', st, it, size(h), ripple
        if (size(h) > 0) print '(*(es25.16,1x))', h
        deallocate(h,bands,desired,weight)
    end do
end program remez_driver
