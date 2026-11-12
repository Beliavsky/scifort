! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Deterministic FIR2D designs for independent NumPy/SciPy comparisons.
program fir2d_anisotropic_driver
    use scifort_kinds, only : dp
    use scifort_signal_fir2d, only : firwin_2d
    implicit none
    real(dp), allocatable :: h(:,:)
    integer :: family, mode, ir, ic, nr, nc, st

    do mode = 1,7
        do ir = 1,5
            do ic = 1,5
                nr = 3+2*ir
                nc = 3+2*ic
                select case (mode)
                case (1)
                    h = firwin_2d([nr,nc],'hann',0.17_dp, &
                        window2='blackman',fc2=[0.36_dp],fs2=2.8_dp,status=st)
                case (2)
                    h = firwin_2d([nr,nc],'blackman',0.39_dp, &
                        window2='hamming',fc2=[0.72_dp],fs2=3.1_dp, &
                        pass_zero='highpass',status=st)
                case (3)
                    h = firwin_2d([nr,nc],'hamming',[0.12_dp,0.29_dp], &
                        window2='hann',fc2=[0.19_dp,0.37_dp], &
                        pass_zero='bandpass',pass_zero2='bandstop',status=st)
                case (4)
                    h = firwin_2d([nr,nc],'kaiser',[0.13_dp,0.35_dp], &
                        window_parameter=6.0_dp,window2='kaiser', &
                        window_parameter2=3.0_dp,fc2=[0.20_dp,0.48_dp], &
                        pass_zero='bandstop',pass_zero2='bandpass',scale2=.false.,status=st)
                case (5)
                    h = firwin_2d([nr,nc],'hann',0.16_dp, &
                        fs=2.5_dp,fs2=3.1_dp,fc2=[0.21_dp,0.61_dp], &
                        pass_zero2='bandpass',status=st)
                case (6)
                    h = firwin_2d([nr,nc],'hamming',[0.08_dp,0.33_dp], &
                        window2='blackman',fc2=[0.43_dp], &
                        pass_zero='bandpass',pass_zero2='highpass', &
                        scale=.false.,scale2=.true.,status=st)
                case (7)
                    nr = 4+2*ir
                    nc = 2+2*ic
                    h = firwin_2d([nr,nc],'hamming',0.23_dp, &
                        window2='hann',fc2=[0.37_dp],status=st)
                end select
                call print_design('S',mode,nr,nc,st,h)
            end do
        end do
    end do

    do mode = 1,4
        do ir = 1,4
            do ic = 1,4
                nr = 3+2*ir
                nc = 4+2*ic
                select case (mode)
                case (1)
                    h = firwin_2d([nr,nc],'hamming',0.18_dp,circular=.true., &
                        ellipse_axes=[1.4_dp,0.65_dp],status=st)
                case (2)
                    h = firwin_2d([nr,nc],'blackman',0.27_dp,circular=.true., &
                        ellipse_axes=[0.65_dp,1.8_dp],ellipse_angle=0.37_dp,status=st)
                case (3)
                    h = firwin_2d([nr,nc],'hamming',[0.16_dp,0.29_dp], &
                        circular=.true.,pass_zero='bandpass', &
                        ellipse_axes=[0.9_dp,1.2_dp],ellipse_angle=-0.58_dp,status=st)
                case (4)
                    h = firwin_2d([nr,nc],'hann',0.26_dp,fs=3.0_dp,circular=.true., &
                        ellipse_axes=[1.0_dp,1.0_dp],ellipse_angle=1.2_dp,status=st)
                end select
                call print_design('E',mode,nr,nc,st,h)
            end do
        end do
    end do
contains
    subroutine print_design(family, mode, nr, nc, status, h)
        character(*), intent(in) :: family !! 'S' for separable or 'E' for elliptical
        integer, intent(in) :: mode !! reference case number
        integer, intent(in) :: nr !! requested first-axis length
        integer, intent(in) :: nc !! requested second-axis length
        integer, intent(in) :: status !! status returned by design routine
        real(dp), intent(in) :: h(:,:) !! resulting FIR coefficients in Fortran array order

        write(*,'(a,6(1x,i0))',advance='no') family,mode,nr,nc,status,size(h,1),size(h,2)
        write(*,'(*(1x,es24.16))') reshape(h,[size(h)])
    end subroutine print_design
end program fir2d_anisotropic_driver
