! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Text-stream interface for independent NumPy/SciPy reference comparisons.
program reference_signal_savgol_irregular
    use scifort_kinds, only : dp
    use scifort_signal_savgol, only : savgol_irregular
    implicit none
    real(dp), allocatable :: time(:), r(:), imag(:), weights(:), queries(:), output(:)
    complex(dp), allocatable :: complex_output(:)
    integer :: ncases, k, n, m, wl, order, derivative, cmode, weighted, stat, i

    read(*,*) ncases
    do k = 1, ncases
        read(*,*) n, m, wl, order, derivative, cmode, weighted
        allocate(time(n),r(n),imag(n),weights(n),queries(m))
        read(*,*) (time(i),i=1,n)
        read(*,*) (r(i),i=1,n)
        read(*,*) (imag(i),i=1,n)
        read(*,*) (weights(i),i=1,n)
        read(*,*) (queries(i),i=1,m)
        if (cmode == 0) then
            if (weighted == 1) then
                output = savgol_irregular(time,r,wl,order,deriv=derivative, &
                    weights=weights,query=queries,status=stat)
            else
                output = savgol_irregular(time,r,wl,order,deriv=derivative,query=queries,status=stat)
            end if
            write(*,'(I0,1X,I0)') stat, size(output)
            do i = 1, size(output)
                write(*,'(ES25.16E3,1X,ES25.16E3)') output(i),0.0_dp
            end do
        else
            if (weighted == 1) then
                complex_output = savgol_irregular(time,cmplx(r,imag,dp),wl,order,deriv=derivative, &
                    weights=weights,query=queries,status=stat)
            else
                complex_output = savgol_irregular(time,cmplx(r,imag,dp),wl,order,deriv=derivative, &
                    query=queries,status=stat)
            end if
            write(*,'(I0,1X,I0)') stat, size(complex_output)
            do i = 1, size(complex_output)
                write(*,'(ES25.16E3,1X,ES25.16E3)') real(complex_output(i),dp),aimag(complex_output(i))
            end do
        end if
        deallocate(time,r,imag,weights,queries)
    end do
end program reference_signal_savgol_irregular
