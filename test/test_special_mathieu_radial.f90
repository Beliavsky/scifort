! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_special_mathieu_radial
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : dp, mathieu_modcem1, mathieu_modcem2, &
        mathieu_modsem1, mathieu_modsem2
    implicit none
    real(dp) :: derivative, value
    real(dp) :: values(3), derivatives(3)

    call mathieu_modcem1(0, 0.1_dp, 0.5_dp, value, derivative)
    call check_close(value, 0.92582719629762533_dp, 5.0e-12_dp, 'modcem1 m0 value')
    call check_close(derivative, -0.11312798006378917_dp, 5.0e-12_dp, 'modcem1 m0 derivative')
    call mathieu_modcem1(2, 10.0_dp, 1.0_dp, value, derivative)
    call check_close(value, 0.099390623888148719_dp, 8.0e-12_dp, 'modcem1 m2 value')
    call check_close(derivative, 2.0854684665454077_dp, 8.0e-12_dp, 'modcem1 m2 derivative')
    call mathieu_modcem1(3, 1.0_dp, -0.5_dp, value, derivative)
    call check_close(value, 0.088368463436639053_dp, 5.0e-12_dp, 'modcem1 negative x value')
    call check_close(derivative, -0.19440766180656968_dp, 5.0e-12_dp, 'modcem1 negative x derivative')

    call mathieu_modcem2(0, 1.0_dp, 0.0_dp, value, derivative)
    call check_close(value, -0.094272874162047604_dp, 8.0e-12_dp, 'modcem2 origin value')
    call check_close(derivative, 0.94449953478604609_dp, 8.0e-12_dp, 'modcem2 origin derivative')
    call mathieu_modcem2(2, 10.0_dp, 1.0_dp, value, derivative)
    call check_close(value, -0.25963201902894473_dp, 8.0e-12_dp, 'modcem2 m2 value')
    call check_close(derivative, 0.95748854423426266_dp, 8.0e-12_dp, 'modcem2 m2 derivative')

    call mathieu_modsem1(1, 0.1_dp, 0.5_dp, value, derivative)
    call check_close(value, 0.16055823639367958_dp, 5.0e-12_dp, 'modsem1 m1 value')
    call check_close(derivative, 0.33792378774581339_dp, 5.0e-12_dp, 'modsem1 m1 derivative')
    call mathieu_modsem1(3, 10.0_dp, 1.0_dp, value, derivative)
    call check_close(value, -0.25828405941416616_dp, 8.0e-12_dp, 'modsem1 m3 value')
    call check_close(derivative, 0.99005338105276763_dp, 8.0e-12_dp, 'modsem1 m3 derivative')

    call mathieu_modsem2(1, 1.0_dp, 0.0_dp, value, derivative)
    call check_close(value, -0.70966143924481495_dp, 8.0e-12_dp, 'modsem2 origin value')
    call check_close(derivative, 0.17190795458198405_dp, 8.0e-12_dp, 'modsem2 origin derivative')
    call mathieu_modsem2(3, 10.0_dp, 1.0_dp, value, derivative)
    call check_close(value, -0.10356844959906460_dp, 8.0e-12_dp, 'modsem2 m3 value')
    call check_close(derivative, -2.0678065842817559_dp, 8.0e-12_dp, 'modsem2 m3 derivative')

    call mathieu_modcem1(0, 1.0_dp, 1.0_dp, value, derivative)
    call check_close(value, -0.10932055554864945_dp, 8.0e-12_dp, 'modcem1 radian argument')

    call mathieu_modcem1([0, 1, 2], [1.0_dp, 1.0_dp, 1.0_dp], &
        [0.0_dp, 0.5_dp, 1.0_dp], values, derivatives)
    call check_close(values(1), 0.6740286775385139_dp, 8.0e-12_dp, 'modcem1 elemental 1')
    call check_close(values(2), 0.66102716577587128_dp, 8.0e-12_dp, 'modcem1 elemental 2')
    call check_close(values(3), 0.4930389783192972_dp, 8.0e-12_dp, 'modcem1 elemental 3')

    call mathieu_modcem1(-1, 1.0_dp, 0.5_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'modcem1 negative order')
    call mathieu_modsem1(0, 1.0_dp, 0.5_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'modsem1 zero order')
    call mathieu_modcem1(0, -1.0_dp, 0.5_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'modcem1 negative q')
    call mathieu_modcem2(0, 0.0_dp, 0.5_dp, value, derivative)
    call require(ieee_is_nan(value) .and. ieee_is_nan(derivative), 'modcem2 degenerate q')

    print *, 'test_special_mathieu_radial: PASS'
contains
    subroutine check_close(actual, expected, relative_tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! independent reference value
        real(dp), intent(in) :: relative_tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label
        if (abs(actual - expected) > relative_tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close
    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label
        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require
end program test_special_mathieu_radial
