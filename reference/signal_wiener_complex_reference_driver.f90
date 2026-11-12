! SPDX-License-Identifier: MIT
! SciPy comparison input generator for complex 1D, 2D and 3D local Wiener filters.
program signal_wiener_complex_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_denoise, only : wiener, denoise_status_success
    implicit none
    complex(dp) :: x1(13), x2(8, 7), x3(5, 4, 3)
    complex(dp), allocatable :: z1(:), z2(:,:), z3(:,:,:)
    integer :: i, j, k, mode, st

    do i = 1, size(x1)
        x1(i) = cmplx(real(mod(17*i, 29)-14,dp)/7.0_dp, &
            real(mod(19*i+3,23)-11,dp)/9.0_dp,dp)
    end do
    do j = 1, size(x2,2)
        do i = 1, size(x2,1)
            x2(i,j) = cmplx(real(mod(7*i+11*j,37)-18,dp)/8.0_dp, &
                real(mod(13*i+26*j,31)-15,dp)/10.0_dp,dp)
        end do
    end do
    do k = 1, size(x3,3)
        do j = 1, size(x3,2)
            do i = 1, size(x3,1)
                x3(i,j,k) = cmplx(real(mod(11*i+7*j+5*k,41)-20,dp)/11.0_dp, &
                    real(mod(17*i+13*j+19*k,43)-21,dp)/13.0_dp,dp)
            end do
        end do
    end do

    do mode = 1, 9
        select case (mode)
        case (1)
            z1 = wiener(x1, status=st)
        case (2)
            z1 = wiener(x1, kernel_size=2, status=st)
        case (3)
            z1 = wiener(x1, kernel_size=5, noise=0.12_dp, status=st)
        case (4)
            z1 = wiener(x1, kernel_size=3, noise=cmplx(0.3_dp,-0.14_dp,dp), status=st)
        case (5)
            z1 = wiener(x1, kernel_shape=[4], noise=cmplx(-0.2_dp,0.13_dp,dp),status=st)
        case (6)
            z1 = wiener(x1, kernel_shape=[1], status=st)
        case (7)
            z1 = wiener(x1, kernel_size=7, noise=0.0_dp, status=st)
        case (8)
            z1 = wiener(x1, kernel_size=4, noise=cmplx(0.0_dp,0.0_dp,dp),status=st)
        case (9)
            z1 = wiener(x1, kernel_size=9, noise=0.65_dp, status=st)
        end select
        if (st /= denoise_status_success) error stop 'complex Wiener 1D reference invalid'
        call emit_case(mode,reshape(z1,[size(z1)]))
    end do

    do mode = 1, 9
        select case (mode)
        case (1)
            z2 = wiener(x2, status=st)
        case (2)
            z2 = wiener(x2, kernel_size=2, status=st)
        case (3)
            z2 = wiener(x2, kernel_shape=[3,2], noise=0.12_dp, status=st)
        case (4)
            z2 = wiener(x2, kernel_shape=[4,3], noise=cmplx(0.3_dp,-0.14_dp,dp), status=st)
        case (5)
            z2 = wiener(x2, kernel_shape=[1,5], status=st)
        case (6)
            z2 = wiener(x2, kernel_shape=[7,1], noise=cmplx(-0.2_dp,0.13_dp,dp), status=st)
        case (7)
            z2 = wiener(x2, kernel_shape=[2,4], noise=0.0_dp, status=st)
        case (8)
            z2 = wiener(x2, kernel_size=5, noise=cmplx(0.0_dp,0.0_dp,dp), status=st)
        case (9)
            z2 = wiener(x2, kernel_shape=[5,5], noise=0.65_dp, status=st)
        end select
        if (st /= denoise_status_success) error stop 'complex Wiener 2D reference invalid'
        call emit_case(100+mode,reshape(z2,[size(z2)]))
    end do

    do mode = 1, 9
        select case (mode)
        case (1)
            z3 = wiener(x3, status=st)
        case (2)
            z3 = wiener(x3, kernel_size=2, status=st)
        case (3)
            z3 = wiener(x3, kernel_shape=[3,2,5], noise=0.12_dp, status=st)
        case (4)
            z3 = wiener(x3, kernel_shape=[4,3,2], noise=cmplx(0.3_dp,-0.14_dp,dp), status=st)
        case (5)
            z3 = wiener(x3, kernel_shape=[1,5,2], status=st)
        case (6)
            z3 = wiener(x3, kernel_shape=[5,1,3], noise=cmplx(-0.2_dp,0.13_dp,dp), status=st)
        case (7)
            z3 = wiener(x3, kernel_shape=[2,4,1], noise=0.0_dp, status=st)
        case (8)
            z3 = wiener(x3, kernel_size=3, noise=cmplx(0.0_dp,0.0_dp,dp), status=st)
        case (9)
            z3 = wiener(x3, kernel_shape=[3,3,3], noise=0.65_dp, status=st)
        end select
        if (st /= denoise_status_success) error stop 'complex Wiener 3D reference invalid'
        call emit_case(200+mode,reshape(z3,[size(z3)]))
    end do

contains

    subroutine emit_case(label, values)
        integer, intent(in) :: label !! case identity; 1D, 2D+100, or 3D+200
        complex(dp), intent(in) :: values(:) !! Fortran column-major flattened result
        integer :: n
        print '(a,1x,i0,1x,i0)', 'CASE', label, size(values)
        do n = 1, size(values)
            print '(es26.17,1x,es26.17)', real(values(n),dp), aimag(values(n))
        end do
    end subroutine emit_case

end program signal_wiener_complex_reference_driver
