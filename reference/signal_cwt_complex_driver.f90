! SPDX-License-Identifier: MIT
! Deterministic complex CWT output for independent SciPy verification.
program signal_cwt_complex_driver
    use scifort_kinds, only : dp
    use scifort_signal_cwt_peaks, only : cwt_morlet2, cwt_complex_wavelet, &
        cwt_peak_success, morlet2_wavelet
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: lengths(8) = [7,8,16,17,31,32,61,64]
    real(dp), allocatable :: x(:), widths(:)
    complex(dp), allocatable :: z(:), y(:,:)
    real(dp) :: omega
    integer :: ti, mode, kind, trial, n, j, k, state, icase

    icase = 0
    do kind = 1, 2
        do mode = 1, 3
            do ti = 1, size(lengths)
                n = lengths(ti)
                allocate(x(n), z(n), widths(4))
                do trial = 1, 5
                    icase = icase+1
                    omega = 1.7_dp + 0.37_dp*trial
                    do j = 1, n
                        x(j) = sin(0.14_dp*j+0.31_dp*trial) + &
                            0.3_dp*cos(0.65_dp*j-0.17_dp*trial)
                        z(j) = cmplx(x(j), &
                            0.6_dp*cos(0.37_dp*j+0.11_dp*trial)-sin(0.07_dp*j),dp)
                    end do
                    widths = [0.43_dp + 0.08_dp*trial, 1.1_dp + 0.23_dp*trial, &
                        2.2_dp + 0.51_dp*trial, 4.3_dp + 0.9_dp*trial]
                    if (kind == 1) then
                        select case(mode)
                        case(1)
                            y = cwt_morlet2(x,widths,omega,status=state)
                        case(2)
                            y = cwt_complex_wavelet(x,widths,morlet2_wavelet,status=state)
                        case(3)
                            y = cwt_complex_wavelet(x,widths,asymmetric_wavelet,status=state)
                        end select
                    else
                        select case(mode)
                        case(1)
                            y = cwt_morlet2(z,widths,omega,status=state)
                        case(2)
                            y = cwt_complex_wavelet(z,widths,morlet2_wavelet,status=state)
                        case(3)
                            y = cwt_complex_wavelet(z,widths,asymmetric_wavelet,status=state)
                        end select
                    end if
                    if (state /= cwt_peak_success) error stop 'CWT reference failed'
                    write(*,'(a,5(1x,i0),1x,es24.16e3)') 'CASE',icase,kind,mode,n,size(widths),omega
                    do j = 1, n
                        if (kind == 1) then
                            write(*,'(2(1x,es24.16e3))') x(j),0.0_dp
                        else
                            write(*,'(2(1x,es24.16e3))') real(z(j),dp),aimag(z(j))
                        end if
                    end do
                    do j = 1, size(widths)
                        write(*,'(1x,es24.16e3)') widths(j)
                    end do
                    do j = 1, size(widths)
                        do k = 1, n
                            write(*,'(2(1x,es24.16e3))') real(y(j,k),dp),aimag(y(j,k))
                        end do
                    end do
                end do
                deallocate(x,z,widths)
            end do
        end do
    end do
    write(*,'(a,1x,i0)') 'TOTAL',icase
contains
    function asymmetric_wavelet(points,width) result(w)
        real(dp), intent(in) :: points !! possibly noninteger sample count
        real(dp), intent(in) :: width !! positive wavelet scale
        complex(dp), allocatable :: w(:)
        real(dp) :: r
        integer :: j

        allocate(w(ceiling(points)))
        do j = 1, size(w)
            r = (real(j-1,dp)-0.5_dp*(points-1.0_dp))/width
            w(j) = exp(-0.5_dp*r*r)*cmplx(1.0_dp+0.1_dp*r,0.7_dp*sin(r)+0.2_dp,dp)
        end do
    end function asymmetric_wavelet
end program signal_cwt_complex_driver
