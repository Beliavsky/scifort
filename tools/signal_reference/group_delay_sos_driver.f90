! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Streamed group-delay reference driver: structured input is supplied by the
! independent SciPy comparison script; output contains frequencies and delays.
program group_delay_sos_driver
    use scifort_kinds, only : dp
    use scifort_signal_design, only : group_delay, group_delay_sos, group_delay_result
    implicit none
    integer :: nc, ns, nf, is_complex, i, j, k, ignored
    real(dp) :: rate, components(12)
    real(dp), allocatable :: w(:)
    complex(dp), allocatable :: sos(:,:)
    type(group_delay_result) :: out

    read(*,*) nc
    do i = 1, nc
        read(*,*) ns, nf, is_complex, rate
        allocate(sos(ns,6),w(nf))
        do j = 1, ns
            read(*,*) components
            do k = 1, 6
                sos(j,k) = cmplx(components(2*k-1),components(2*k),dp)
            end do
        end do
        read(*,*) w
        if (is_complex == 0) then
            out = group_delay_sos(real(sos,dp),fs=rate,w=w)
        else
            out = group_delay(sos,fs=rate,w=w)
        end if
        ignored = 0
        if (allocated(out%delay)) ignored = size(out%delay)
        write(*,'(i0,1x,i0)') out%status, ignored
        if (ignored > 0) then
            do j = 1, ignored
                write(*,'(es26.17,1x,es26.17)') out%frequencies(j),out%delay(j)
            end do
        end if
        deallocate(w,sos)
    end do
end program group_delay_sos_driver
