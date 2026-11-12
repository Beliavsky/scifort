! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_fresnel_zeros
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_constants, only : scifort_pi
    use scifort_special, only : dp, fresnel, fresnel_zeros, fresnelc_zeros, &
        fresnels_zeros, fresnel_zeros_status_invalid_argument, &
        fresnel_zeros_status_success
    implicit none

    complex(dp), allocatable :: c_zeros(:)
    complex(dp), allocatable :: s_zeros(:)
    complex(dp), allocatable :: pair_c(:)
    complex(dp), allocatable :: pair_s(:)
    complex(dp) :: c_value
    complex(dp) :: s_value
    integer :: k
    integer :: status

    call fresnelc_zeros(5, c_zeros, status)
    call require(status == fresnel_zeros_status_success, 'fresnelc status')
    call check_complex('fresnelc zero 1', c_zeros(1), &
        cmplx(1.7436674861820027_dp, 0.3057350636296493_dp, kind=dp), 2.0e-13_dp)
    call check_complex('fresnelc zero 5', c_zeros(5), &
        cmplx(4.361063517028432_dp, 0.19066973234885032_dp, kind=dp), 2.0e-13_dp)

    call fresnels_zeros(5, s_zeros, status)
    call require(status == fresnel_zeros_status_success, 'fresnels status')
    call check_complex('fresnels zero 1', s_zeros(1), &
        cmplx(2.0092570118220077_dp, 0.2885478973334392_dp, kind=dp), 2.0e-13_dp)
    call check_complex('fresnels zero 5', s_zeros(5), &
        cmplx(4.474189295216029_dp, 0.18768858914403957_dp, kind=dp), 2.0e-13_dp)

    call fresnel_zeros(12, pair_s, pair_c, status)
    call require(status == fresnel_zeros_status_success, 'fresnel_zeros status')
    do k = 1, 12
        call fresnel(pair_s(k), s_value, c_value)
        call require(abs(s_value) <= 2.0e-12_dp * max(1.0_dp, abs(pair_s(k))), &
            'fresnels residual')
        call fresnel(pair_c(k), s_value, c_value)
        call require(abs(c_value) <= 2.0e-12_dp * max(1.0_dp, abs(pair_c(k))), &
            'fresnelc residual')
        call require(real(pair_s(k), dp) > 0.0_dp .and. aimag(pair_s(k)) > 0.0_dp, &
            'fresnels first quadrant')
        call require(real(pair_c(k), dp) > 0.0_dp .and. aimag(pair_c(k)) > 0.0_dp, &
            'fresnelc first quadrant')
        if (k > 1) then
            call require(abs(pair_s(k)) > abs(pair_s(k - 1)), 'fresnels ordering')
            call require(abs(pair_c(k)) > abs(pair_c(k - 1)), 'fresnelc ordering')
        end if
    end do

    call fresnelc_zeros(332, c_zeros, status)
    call require(status == fresnel_zeros_status_success, 'fresnelc 332 status')
    call require(ieee_is_finite(real(c_zeros(332), dp)) .and. &
        ieee_is_finite(aimag(c_zeros(332))), 'fresnelc 332 finite')
    call fresnel(c_zeros(332), s_value, c_value)
    call require(abs(c_value) <= 2.0e-11_dp * abs(c_zeros(332)), 'fresnelc 332 residual')

    call fresnels_zeros(0, s_zeros, status)
    call require(status == fresnel_zeros_status_invalid_argument, 'fresnels invalid status')
    call require(size(s_zeros) == 0, 'fresnels invalid output')

    print *, 'special Fresnel-zero tests passed'

contains

    subroutine check_complex(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex

    subroutine require(condition, name)
        logical, intent(in) :: condition !! assertion condition expected to be true
        character(len=*), intent(in) :: name !! check label printed on failure

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine require

end program test_special_fresnel_zeros
