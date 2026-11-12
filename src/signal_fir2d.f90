! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Two-dimensional window-method FIR design: separable outer products and
! the radial-window construction used by scipy.signal.firwin_2d.
! Independent axis cutoffs and elliptical radial profiles are SciFort extensions.
module scifort_signal_fir2d
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_design, only : firwin, design_status_success, design_status_invalid_input
    implicit none
    private
    interface firwin_2d
        module procedure firwin_2d_scalar
        module procedure firwin_2d_vector
    end interface
    public :: firwin_2d
contains

    function firwin_2d_scalar(hsize, window, fc, window2, circular, fs, &
            pass_zero, scale, status, window_parameter, window_parameter2, &
            fc2, fs2, pass_zero2, scale2, ellipse_axes, ellipse_angle) result(taps)
        integer, intent(in) :: hsize(2) !! number of taps along the two filter axes
        character(*), intent(in) :: window !! first-axis or radial window name
        real(dp), intent(in) :: fc !! first-axis scalar cutoff in the units of fs
        character(*), intent(in), optional :: window2 !! second-axis window, default window
        logical, intent(in), optional :: circular !! use a radial profile, default false
        real(dp), intent(in), optional :: fs !! first-axis sampling frequency, default 2
        character(*), intent(in), optional :: pass_zero !! first-axis passband type
        logical, intent(in), optional :: scale !! normalize first-axis passband, default true
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), intent(in), optional :: window_parameter !! first or radial window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second-axis window shape parameter
        real(dp), intent(in), optional :: fc2(:) !! second-axis cutoff edges, default first-axis fc
        real(dp), intent(in), optional :: fs2 !! second-axis sample rate, default fs
        character(*), intent(in), optional :: pass_zero2 !! second-axis passband, default pass_zero
        logical, intent(in), optional :: scale2 !! second-axis scaling flag, default scale
        real(dp), intent(in), optional :: ellipse_axes(2) !! ellipse scales (>0), default [1,1]
        real(dp), intent(in), optional :: ellipse_angle !! ellipse rotation in radians, default 0
        real(dp), allocatable :: taps(:,:)

        taps = firwin_2d_vector(hsize,window,[fc],window2,circular,fs,pass_zero, &
            scale,status,window_parameter,window_parameter2,fc2,fs2,pass_zero2, &
            scale2,ellipse_axes,ellipse_angle)
    end function firwin_2d_scalar

    function firwin_2d_vector(hsize, window, fc, window2, circular, fs, &
            pass_zero, scale, status, window_parameter, window_parameter2, &
            fc2, fs2, pass_zero2, scale2, ellipse_axes, ellipse_angle) result(taps)
        integer, intent(in) :: hsize(2) !! row and column FIR lengths, both positive
        character(*), intent(in) :: window !! first-axis or radial window name
        real(dp), intent(in) :: fc(:) !! first-axis ordered cutoff edges in (0,fs/2)
        character(*), intent(in), optional :: window2 !! second-axis window, default window
        logical, intent(in), optional :: circular !! use a radial profile, default false
        real(dp), intent(in), optional :: fs !! first-axis sampling frequency, default 2
        character(*), intent(in), optional :: pass_zero !! first-axis passband type
        logical, intent(in), optional :: scale !! normalize first-axis passband, default true
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), intent(in), optional :: window_parameter !! first or radial window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second-axis window shape parameter
        real(dp), intent(in), optional :: fc2(:) !! second-axis ordered cutoff edges, default fc
        real(dp), intent(in), optional :: fs2 !! second-axis sample rate, default fs
        character(*), intent(in), optional :: pass_zero2 !! second-axis passband, default pass_zero
        logical, intent(in), optional :: scale2 !! second-axis scaling flag, default scale
        real(dp), intent(in), optional :: ellipse_axes(2) !! radial axis scales (>0), default [1,1]
        real(dp), intent(in), optional :: ellipse_angle !! radial rotation in radians, default 0
        real(dp), allocatable :: taps(:,:), row_taps(:), col_taps(:), radial_taps(:)
        real(dp), allocatable :: second_cutoffs(:)
        real(dp) :: xr, xc, r, t, position, rate2, axes(2), angle, cosa, sina, u, v
        character(:), allocatable :: second_window, second_mode
        logical :: radial, use_scale2
        integer :: i, j, nr, st, idx

        allocate(taps(0,0))
        if (present(status)) status = design_status_invalid_input
        if (any(hsize < 1) .or. size(fc) < 1) return
        if (any(hsize > 100000)) return
        if (real(hsize(1),dp)*real(hsize(2),dp) > 20000000.0_dp) return
        radial = .false.
        if (present(circular)) radial = circular
        if (radial) then
            ! Independent sampling rates and cutoffs describe separable filters,
            ! so do not silently ignore them in the radial design.
            if (present(fc2) .or. present(fs2) .or. present(pass_zero2) .or. present(scale2)) return
        else
            if (present(ellipse_axes) .or. present(ellipse_angle)) return
        end if
        if (.not. radial) then
            second_window = window
            if (present(window2)) second_window = window2
            second_cutoffs = fc
            if (present(fc2)) second_cutoffs = fc2
            rate2 = 2.0_dp
            if (present(fs)) rate2 = fs
            if (present(fs2)) rate2 = fs2
            second_mode = 'lowpass'
            if (size(second_cutoffs) > 1) second_mode = 'bandstop'
            if (present(pass_zero)) second_mode = pass_zero
            if (present(pass_zero2)) second_mode = pass_zero2
            use_scale2 = .true.
            if (present(scale)) use_scale2 = scale
            if (present(scale2)) use_scale2 = scale2
            row_taps = firwin(hsize(1),fc,window=window,fs=fs,pass_zero=pass_zero, &
                scale=scale,status=st,window_parameter=window_parameter)
            if (st /= design_status_success) return
            if (present(window_parameter2)) then
                col_taps = firwin(hsize(2),second_cutoffs,window=second_window,fs=rate2, &
                    pass_zero=second_mode,scale=use_scale2,status=st, &
                    window_parameter=window_parameter2)
            else
                col_taps = firwin(hsize(2),second_cutoffs,window=second_window,fs=rate2, &
                    pass_zero=second_mode,scale=use_scale2,status=st, &
                    window_parameter=window_parameter)
            end if
            if (st /= design_status_success) return
            deallocate(taps)
            allocate(taps(hsize(1),hsize(2)))
            do j = 1,hsize(2)
                taps(:,j) = row_taps*col_taps(j)
            end do
        else
            ! For non-square shapes SciPy's meshgrid convention transposes the
            ! requested axes. Preserve this convention for compatibility.
            axes = [1.0_dp,1.0_dp]
            if (present(ellipse_axes)) axes = ellipse_axes
            if (.not. all(ieee_is_finite(axes))) return
            if (any(axes <= 0.0_dp)) return
            angle = 0.0_dp
            if (present(ellipse_angle)) angle = ellipse_angle
            if (.not. ieee_is_finite(angle)) return
            cosa = cos(angle)
            sina = sin(angle)
            nr = 8*maxval(hsize)
            radial_taps = firwin(nr,fc,window=window,fs=fs,pass_zero=pass_zero, &
                scale=scale,status=st,window_parameter=window_parameter)
            if (st /= design_status_success) return
            deallocate(taps)
            allocate(taps(hsize(2),hsize(1)))
            do j = 1,hsize(1)
                xr = -1.0_dp
                if (hsize(1) > 1) xr = -1.0_dp+2.0_dp*real(j-1,dp)/real(hsize(1)-1,dp)
                do i = 1,hsize(2)
                    xc = -1.0_dp
                    if (hsize(2) > 1) xc = -1.0_dp+2.0_dp*real(i-1,dp)/real(hsize(2)-1,dp)
                    u = cosa*xr+sina*xc
                    v = -sina*xr+cosa*xc
                    ! Check the axis limits before dividing; even a tiny positive
                    ! ellipse scale cannot cause an intermediate overflow.
                    if (abs(u) >= axes(1) .or. abs(v) >= axes(2)) then
                        taps(i,j) = radial_taps(nr)
                        cycle
                    end if
                    r = sqrt((u/axes(1))**2+(v/axes(2))**2)
                    position = r*real(nr-1,dp)
                    if (position <= 0.0_dp) then
                        taps(i,j) = radial_taps(1)
                    else if (position >= real(nr-1,dp)) then
                        taps(i,j) = radial_taps(nr)
                    else
                        idx = floor(position)+1
                        t = position-real(idx-1,dp)
                        taps(i,j) = (1.0_dp-t)*radial_taps(idx)+t*radial_taps(idx+1)
                    end if
                end do
            end do
        end if
        if (.not. all(ieee_is_finite(taps))) then
            deallocate(taps)
            allocate(taps(0,0))
            return
        end if
        if (present(status)) status = design_status_success
    end function firwin_2d_vector
end module scifort_signal_fir2d
