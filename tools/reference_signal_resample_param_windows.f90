! SPDX-License-Identifier: MIT
! Numeric reference output for SciPy comparison of named Fourier-resampling windows.
program reference_signal_resample_param_windows
    use scifort_kinds, only : dp
    use scifort_signal_fourier_resample, only : resample
    implicit none
    integer :: cases, trial, mode, kind, nx, ny, st, i
    real(dp) :: param, param2
    real(dp), allocatable :: xr(:), xi(:), yr(:)
    complex(dp), allocatable :: x(:), y(:)
    character(len=32) :: family
    real(dp), parameter :: coefficients(3) = [0.6_dp, 0.3_dp, 0.1_dp]

    read(*,*) cases
    do trial = 1, cases
        read(*,*) mode, kind, nx, ny, param, param2
        allocate(xr(nx), xi(nx), x(nx))
        read(*,*) xr
        read(*,*) xi
        x = cmplx(xr,xi,dp)
        select case(kind)
        case(1)
            family = 'kaiser'
        case(2)
            family = 'gaussian'
        case(3)
            family = 'tukey'
        case(4)
            family = 'general_gaussian'
        case(5)
            family = 'taylor'
        case(6)
            family = 'general_cosine'
        case(7)
            family = 'dpss'
        case default
            error stop 'invalid test family'
        end select
        if (mode == 0) then
            if (kind == 6) then
                yr = resample(xr,ny,window=trim(family), &
                    window_coefficients=coefficients,status=st)
            else if (kind == 4 .or. kind == 5) then
                yr = resample(xr,ny,window=trim(family), &
                    window_parameter=param, window_parameter2=param2,status=st)
            else
                yr = resample(xr,ny,window=trim(family), &
                    window_parameter=param,status=st)
            end if
            y = cmplx(yr,0.0_dp,dp)
        else
            if (kind == 6) then
                y = resample(x,ny,window=trim(family), &
                    window_coefficients=coefficients,status=st, &
                    domain=merge('freq','time',mode==2))
            else if (kind == 4 .or. kind == 5) then
                y = resample(x,ny,window=trim(family), &
                    window_parameter=param,window_parameter2=param2,status=st, &
                    domain=merge('freq','time',mode==2))
            else
                y = resample(x,ny,window=trim(family), &
                    window_parameter=param,status=st, &
                    domain=merge('freq','time',mode==2))
            end if
        end if
        print '(i0,1x,i0)',st,size(y)
        do i = 1, size(y)
            print '(es26.17e3,1x,es26.17e3)',real(y(i),dp),aimag(y(i))
        end do
        deallocate(xr,xi,x,y)
    end do
end program reference_signal_resample_param_windows
