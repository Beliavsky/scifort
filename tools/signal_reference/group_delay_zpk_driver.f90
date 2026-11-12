! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Independent Python reference sends ZPK arrays and sampling grids.
program group_delay_zpk_driver
    use scifort_kinds, only : dp
    use scifort_signal_design, only : group_delay_zpk, group_delay_result
    implicit none
    integer :: nc, nz, np, nf, i, j, kind, observed
    real(dp) :: rate, gr, gi, zr, zi
    real(dp), allocatable :: w(:)
    complex(dp), allocatable :: z(:), p(:)
    type(group_delay_result) :: out

    read(*,*) nc
    do i = 1, nc
        read(*,*) nz, np, nf, kind, rate, gr, gi
        allocate(z(nz),p(np),w(nf))
        do j = 1, nz
            read(*,*) zr, zi
            z(j) = cmplx(zr,zi,dp)
        end do
        do j = 1, np
            read(*,*) zr, zi
            p(j) = cmplx(zr,zi,dp)
        end do
        read(*,*) w
        if (kind == 0) then
            out = group_delay_zpk(real(z,dp),real(p,dp),gr,fs=rate,w=w)
        else
            out = group_delay_zpk(z,p,cmplx(gr,gi,dp),fs=rate,w=w)
        end if
        observed = 0
        if (allocated(out%delay)) observed = size(out%delay)
        write(*,'(i0,1x,i0)') out%status, observed
        do j = 1, observed
            write(*,'(es26.17,1x,es26.17)') out%frequencies(j),out%delay(j)
        end do
        deallocate(z,p,w)
    end do
end program group_delay_zpk_driver
