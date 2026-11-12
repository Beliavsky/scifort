! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program signal_detrend_multiaxis_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_detrend, only : detrend
    implicit none
    integer :: cases, case_id, rank, isc, trend, axis, nb, dims(3), total, status, i
    integer, allocatable :: bp(:)
    real(dp), allocatable :: xr(:), xi(:)
    character(8) :: mode

    read(*,*) cases
    do case_id=1,cases
        read(*,*) rank,isc,trend,axis,nb,dims
        allocate(bp(nb))
        if (nb > 0) read(*,*) bp
        total = product(dims(:rank))
        allocate(xr(total),xi(total))
        if (total > 0) then
            read(*,*) xr
            if (isc == 1) read(*,*) xi
        end if
        if (isc == 0) xi = 0.0_dp
        mode = 'linear'
        if (trend == 1) mode = 'constant'
        select case (rank)
        case (1)
            if (isc == 0) then
                call report_real(detrend(xr,kind=trim(mode),breaks=bp,status=status))
            else
                call report_complex(detrend(cmplx(xr,xi,dp),kind=trim(mode),breaks=bp,status=status))
            end if
        case (2)
            if (isc == 0) then
                call report_real(detrend(reshape(xr,dims(:2)),kind=trim(mode),breaks=bp, &
                    axis=axis,status=status))
            else
                call report_complex(detrend(reshape(cmplx(xr,xi,dp),dims(:2)),kind=trim(mode), &
                    breaks=bp,axis=axis,status=status))
            end if
        case (3)
            if (isc == 0) then
                call report_real(detrend(reshape(xr,dims),kind=trim(mode),breaks=bp, &
                    axis=axis,status=status))
            else
                call report_complex(detrend(reshape(cmplx(xr,xi,dp),dims),kind=trim(mode), &
                    breaks=bp,axis=axis,status=status))
            end if
        end select
        deallocate(bp,xr,xi)
    end do
contains
    subroutine report_real(y)
        real(dp), intent(in) :: y(..) !! arbitrary rank 1 through 3
        real(dp), allocatable :: flat(:)
        select rank(y)
        rank(1)
            flat = y
        rank(2)
            flat = reshape(y,[size(y)])
        rank(3)
            flat = reshape(y,[size(y)])
        end select
        write(*,'(i0,1x,i0)') status,size(flat)
        do i=1,size(flat)
            write(*,'(es26.17e3,1x,es26.17e3)') flat(i),0.0_dp
        end do
    end subroutine report_real

    subroutine report_complex(y)
        complex(dp), intent(in) :: y(..) !! arbitrary rank 1 through 3
        complex(dp), allocatable :: flat(:)
        select rank(y)
        rank(1)
            flat = y
        rank(2)
            flat = reshape(y,[size(y)])
        rank(3)
            flat = reshape(y,[size(y)])
        end select
        write(*,'(i0,1x,i0)') status,size(flat)
        do i=1,size(flat)
            write(*,'(es26.17e3,1x,es26.17e3)') real(flat(i),dp),aimag(flat(i))
        end do
    end subroutine report_complex
end program signal_detrend_multiaxis_reference_driver
