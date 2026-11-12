! SPDX-License-Identifier: MIT
! Compare generalized FIR windows in resample_poly against independent SciPy cases.
program reference_signal_poly_windows
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : resample_poly
    implicit none

    integer :: count, trial, kind, family, n, up, down, pad, st, i
    real(dp), allocatable :: xr(:), xi(:), yr(:)
    complex(dp), allocatable :: z(:), yz(:)
    character(len=8), parameter :: padding(5) = [character(len=8) :: &
        'constant', 'mean', 'line', 'reflect', 'median']

    read(*,*) count
    do trial = 1, count
        read(*,*) kind, family, n, up, down, pad
        allocate(xr(n), xi(n), z(n))
        read(*,*) xr
        read(*,*) xi
        z = cmplx(xr, xi, dp)
        if (kind == 0) then
            select case (family)
            case (1)
                yr = resample_poly(xr, up, down, window='gaussian', &
                    window_parameter=2.75_dp, padtype=trim(padding(pad)), status=st)
            case (2)
                yr = resample_poly(xr, up, down, window='kaiser', &
                    window_parameter=4.5_dp, padtype=trim(padding(pad)), status=st)
            case (3)
                yr = resample_poly(xr, up, down, window='tukey', &
                    window_parameter=0.65_dp, padtype=trim(padding(pad)), status=st)
            case (4)
                yr = resample_poly(xr, up, down, window='general_gaussian', &
                    window_parameter=1.3_dp, window_parameter2=4.4_dp, &
                    padtype=trim(padding(pad)), status=st)
            case (5)
                yr = resample_poly(xr, up, down, window='general_cosine', &
                    window_coefficients=[0.42_dp,0.50_dp,0.08_dp], &
                    padtype=trim(padding(pad)), status=st)
            case (6)
                yr = resample_poly(xr, up, down, window='taylor', &
                    window_parameter=4.0_dp, window_parameter2=32.0_dp, &
                    padtype=trim(padding(pad)), status=st)
            case (7)
                yr = resample_poly(xr, up, down, window='hann', &
                    padtype=trim(padding(pad)), status=st)
            case (8)
                yr = resample_poly(xr, up, down, window='kaiser', beta=5.0_dp, &
                    padtype=trim(padding(pad)), status=st)
            end select
            write(*,'(i0,1x,i0)') st,size(yr)
            do i = 1, size(yr)
                write(*,'(es26.17e3,1x,es26.17e3)') yr(i), 0.0_dp
            end do
        else
            select case (family)
            case (1)
                yz = resample_poly(z, up, down, window='gaussian', &
                    window_parameter=2.75_dp, padtype=trim(padding(pad)), status=st)
            case (2)
                yz = resample_poly(z, up, down, window='kaiser', &
                    window_parameter=4.5_dp, padtype=trim(padding(pad)), status=st)
            case (3)
                yz = resample_poly(z, up, down, window='tukey', &
                    window_parameter=0.65_dp, padtype=trim(padding(pad)), status=st)
            case (4)
                yz = resample_poly(z, up, down, window='general_gaussian', &
                    window_parameter=1.3_dp, window_parameter2=4.4_dp, &
                    padtype=trim(padding(pad)), status=st)
            case (5)
                yz = resample_poly(z, up, down, window='general_cosine', &
                    window_coefficients=[0.42_dp,0.50_dp,0.08_dp], &
                    padtype=trim(padding(pad)), status=st)
            case (6)
                yz = resample_poly(z, up, down, window='taylor', &
                    window_parameter=4.0_dp, window_parameter2=32.0_dp, &
                    padtype=trim(padding(pad)), status=st)
            case (7)
                yz = resample_poly(z, up, down, window='hann', &
                    padtype=trim(padding(pad)), status=st)
            case (8)
                yz = resample_poly(z, up, down, window='kaiser', beta=5.0_dp, &
                    padtype=trim(padding(pad)), status=st)
            end select
            write(*,'(i0,1x,i0)') st,size(yz)
            do i = 1, size(yz)
                write(*,'(es26.17e3,1x,es26.17e3)') real(yz(i),dp),aimag(yz(i))
            end do
        end if
        deallocate(xr, xi, z)
    end do
end program reference_signal_poly_windows
