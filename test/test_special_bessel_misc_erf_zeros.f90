! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_bessel_misc_erf_zeros
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_special, only : bessel_misc_status_invalid_argument, &
        bessel_misc_status_success, besselpoly, dp, erf, erf_zeros, &
        erf_zeros_status_invalid_argument, erf_zeros_status_success, lmbda
    implicit none

    complex(dp), allocatable :: zeros(:)
    complex(dp) :: residual
    real(dp), allocatable :: dl(:)
    real(dp), allocatable :: vl(:)
    real(dp) :: array_values(3)
    integer :: k
    integer :: status

    call lmbda(2.3_dp, 2.0_dp, vl, dl, status)
    call require(status == bessel_misc_status_success, 'lmbda status')
    call require(size(vl) == 3 .and. size(dl) == 3, 'lmbda sequence size')
    call check_close(vl(1), 0.382047946219906815358362760704_dp, 3.0e-14_dp, 'lmbda v0')
    call check_close(vl(2), 0.626220271291433545134813197307_dp, 3.0e-14_dp, 'lmbda v1')
    call check_close(vl(3), 0.730075251963864761124539150069_dp, 3.0e-14_dp, 'lmbda v2')
    call check_close(dl(1), -0.481707900993410493443379535376_dp, 4.0e-14_dp, 'lmbda d0')
    call check_close(dl(2), -0.317424022592984699874920553278_dp, 4.0e-14_dp, 'lmbda d1')
    call check_close(dl(3), -0.238866455546591730022251523287_dp, 4.0e-14_dp, 'lmbda d2')

    call lmbda(0.3_dp, -3.0_dp, vl, dl, status)
    call check_close(vl(1), -0.0534462732734492499582162502494_dp, 5.0e-14_dp, &
        'lmbda negative x value')
    call check_close(dl(1), 0.349411439734674703007167314330_dp, 5.0e-14_dp, &
        'lmbda negative x derivative')

    call lmbda(5.0_dp, 0.0_dp, vl, dl, status)
    call require(maxval(abs(vl - 1.0_dp)) <= 8.0_dp * epsilon(1.0_dp), 'lmbda x0 values')
    call require(maxval(abs(dl)) <= tiny(1.0_dp), 'lmbda x0 derivatives')

    call lmbda(-0.5_dp, 1.0_dp, vl, dl, status)
    call require(status == bessel_misc_status_invalid_argument, 'lmbda invalid status')
    call require(size(vl) == 0 .and. size(dl) == 0, 'lmbda invalid outputs')

    call check_close(besselpoly(1.0_dp, 1.0_dp, 1.0_dp), &
        0.244497183728638808217787634807_dp, 3.0e-14_dp, 'besselpoly ordinary')
    call check_close(besselpoly(1.0_dp, -0.5_dp, 1.0_dp), &
        0.538070676400706851112487004252_dp, 3.0e-13_dp, 'besselpoly endpoint power')
    call check_close(besselpoly(1.0_dp, 1.0_dp, -0.5_dp), &
        0.115795841723399756008325647029_dp, 4.0e-13_dp, 'besselpoly negative order')
    call check_close(besselpoly(10.0_dp, 1.0_dp, 2.0_dp), &
        0.000823220469504581947413748567_dp, 2.0e-12_dp, 'besselpoly oscillatory')
    call check_close(besselpoly(100.0_dp, 1.0_dp, 1.0_dp), &
        0.000100831549666290103844739779_dp, 3.0e-11_dp, 'besselpoly large scale')
    call check_close(besselpoly(100.0_dp, 0.5_dp, 2.3_dp), &
        0.000751238049640820468851117306_dp, 3.0e-10_dp, 'besselpoly large noninteger')
    call check_close(besselpoly(1.0_dp, -1.49_dp, 0.5_dp), &
        112.498892856513809650003621933_dp, 5.0e-13_dp, 'besselpoly near convergence boundary')
    call check_close(besselpoly(10.0_dp, -1.69_dp, 0.7_dp), &
        538.031263278571008587335560486_dp, 2.0e-12_dp, 'besselpoly singular oscillatory')
    call require(.not. ieee_is_finite(besselpoly(1.0_dp, -2.0_dp, 1.0_dp)), &
        'besselpoly divergent endpoint')
    call require(ieee_is_nan(besselpoly(-1.0_dp, 1.0_dp, 1.0_dp)), &
        'besselpoly negative scale')
    array_values = besselpoly([1.0_dp, 1.0_dp, 10.0_dp], &
        [1.0_dp, -1.0_dp, 1.0_dp], [1.0_dp, 1.0_dp, 2.0_dp])
    call check_close(array_values(1), 0.244497183728638808217787634807_dp, &
        3.0e-14_dp, 'besselpoly elemental 1')
    call check_close(array_values(2), 0.849045485440153181214346863956_dp, &
        3.0e-14_dp, 'besselpoly elemental 2')
    call check_close(array_values(3), 0.000823220469504581947413748567_dp, &
        2.0e-12_dp, 'besselpoly elemental 3')

    call erf_zeros(10, zeros, status)
    call require(status == erf_zeros_status_success, 'erf_zeros status')
    call check_complex(zeros(1), cmplx(1.4506161632436827_dp, 1.8809430001533076_dp, kind=dp), &
        3.0e-13_dp, 'erf zero 1')
    call check_complex(zeros(5), cmplx(3.7690055670141998_dp, 4.060697233933303_dp, kind=dp), &
        4.0e-13_dp, 'erf zero 5')
    call check_complex(zeros(10), cmplx(5.452192201109879_dp, 5.688837437036479_dp, kind=dp), &
        8.0e-13_dp, 'erf zero 10')
    do k = 1, 10
        residual = erf(zeros(k))
        call require(abs(residual) <= 3.0e-12_dp * max(1.0_dp, abs(zeros(k))), &
            'erf zero residual')
        call require(real(zeros(k), dp) > 0.0_dp .and. aimag(zeros(k)) > 0.0_dp, &
            'erf zero first quadrant')
        if (k > 1) call require(abs(zeros(k)) > abs(zeros(k - 1)), 'erf zero ordering')
    end do

    call erf_zeros(1000, zeros, status)
    call require(status == erf_zeros_status_success, 'erf_zeros 1000 status')
    call require(ieee_is_finite(real(zeros(1000), dp)) .and. &
        ieee_is_finite(aimag(zeros(1000))), 'erf_zeros 1000 finite')
    residual = erf(zeros(1000))
    call require(abs(residual) <= 2.0e-10_dp * abs(zeros(1000)), 'erf_zeros 1000 residual')

    call erf_zeros(0, zeros, status)
    call require(status == erf_zeros_status_invalid_argument, 'erf_zeros invalid status')
    call require(size(zeros) == 0, 'erf_zeros invalid size')

    print *, 'test_special_bessel_misc_erf_zeros: PASS'

contains

    subroutine check_close(actual, expected, relative_tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: relative_tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > relative_tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_complex(actual, expected, relative_tolerance, label)
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: relative_tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > relative_tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_complex

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_special_bessel_misc_erf_zeros
