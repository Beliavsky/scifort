! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_legendre_q
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_special, only : associated_legendre_status_invalid_argument, &
        associated_legendre_status_success, dp, lqmn, lqn
    implicit none

    complex(dp), allocatable :: cq(:)
    complex(dp), allocatable :: cd(:)
    complex(dp), allocatable :: cqm(:, :)
    complex(dp), allocatable :: cdm(:, :)
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: d(:)
    real(dp), allocatable :: qm(:, :)
    real(dp), allocatable :: dm(:, :)
    integer :: status

    call lqn(4, 0.2_dp, q, d, status)
    call check_true('lqn success', status == associated_legendre_status_success)
    call check_true('lqn zero lower bound', lbound(q, 1) == 0 .and. ubound(q, 1) == 4)
    call check_close('lqn q0 interior', q(0), 0.2027325540540821_dp, 3.0e-15_dp)
    call check_close('lqn q1 interior', q(1), -0.9594534891891836_dp, 4.0e-15_dp)
    call check_close('lqn q4 interior', q(4), 0.4703672858738804_dp, 8.0e-15_dp)
    call check_close('lqn q2 derivative', d(2), -1.836693800900885_dp, 2.0e-14_dp)

    call lqmn(2, 4, 2.0_dp, qm, dm, status)
    call check_true('lqmn success', status == associated_legendre_status_success)
    call check_true('lqmn zero lower bounds', all(lbound(qm) == [0, 0]))
    call check_close('lqmn exterior q00', qm(0, 0), 0.5493061443340549_dp, 3.0e-15_dp)
    call check_close('lqmn exterior q10', qm(1, 0), -0.5773502691896258_dp, 3.0e-15_dp)
    call check_close('lqmn exterior q11', qm(1, 1), -0.20327438748290558_dp, 3.0e-15_dp)
    call check_close('lqmn exterior q24', qm(2, 4), 0.0368160162716521_dp, 3.0e-15_dp)
    call check_close('lqmn exterior derivative q24', dm(2, 4), -0.1102632700721333_dp, 8.0e-15_dp)

    call lqmn(1, 1, -2.0_dp, qm, dm, status)
    call check_close('lqmn negative exterior q00', qm(0, 0), -0.5493061443340549_dp, 3.0e-15_dp)
    call check_close('lqmn negative exterior q01', qm(0, 1), 0.09861228866810969_dp, 3.0e-15_dp)
    call check_close('lqmn negative exterior q10', qm(1, 0), -0.5773502691896258_dp, 3.0e-15_dp)
    call check_close('lqmn negative exterior q11', qm(1, 1), 0.20327438748290558_dp, 3.0e-15_dp)

    call lqmn(2, 4, cmplx(0.2_dp, 0.1_dp, dp), cqm, cdm, status)
    call check_complex_close('lqmn complex interior q00', cqm(0, 0), &
        cmplx(0.2005866181312343_dp, 0.10374811321760134_dp, dp), 4.0e-15_dp)
    call check_complex_close('lqmn complex interior q10', cqm(1, 0), &
        cmplx(-1.0146994934664402_dp, -0.020912752285606087_dp, dp), 5.0e-15_dp)
    call check_complex_close('lqmn complex interior q24', cqm(2, 4), &
        cmplx(-9.363712350715861_dp, -3.5811058614616207_dp, dp), 4.0e-14_dp)
    call check_complex_close('lqmn complex interior derivative q24', cdm(2, 4), &
        cmplx(-37.80421981391196_dp, 13.054613242573865_dp, dp), 2.0e-13_dp)

    call lqmn(2, 4, cmplx(0.8_dp, 0.8_dp, dp), cqm, cdm, status)
    call check_complex_close('lqmn complex exterior q00', cqm(0, 0), &
        cmplx(0.4353744086117919_dp, -0.6987753301712658_dp, dp), 5.0e-15_dp)
    call check_complex_close('lqmn complex exterior q24', cqm(2, 4), &
        cmplx(0.07999745499772007_dp, 0.28563308709590346_dp, dp), 8.0e-15_dp)
    call check_complex_close('lqmn complex exterior derivative q24', cdm(2, 4), &
        cmplx(-1.1117698790290809_dp, -0.12668702867800924_dp, dp), 2.0e-14_dp)

    call lqn(2, cmplx(0.0_dp, 1.0_dp, dp), cq, cd, status)
    call check_complex_close('lqn unit circle canonical q0', cq(0), &
        cmplx(0.0_dp, -0.25_dp * scifort_pi, dp), 4.0e-15_dp)
    call check_complex_close('lqn unit circle canonical q1', cq(1), &
        cmplx(0.25_dp * scifort_pi - 1.0_dp, 0.0_dp, dp), 4.0e-15_dp)
    call check_complex_close('lqn unit circle q0 derivative', cd(0), &
        cmplx(0.5_dp, 0.0_dp, dp), 3.0e-15_dp)

    call lqmn(0, 4, cmplx(0.2_dp, 0.1_dp, dp), cqm, cdm, status)
    call lqn(4, cmplx(0.2_dp, 0.1_dp, dp), cq, cd, status)
    call check_complex_close('lqn agrees lqmn value', cq(4), cqm(0, 4), 2.0e-15_dp)
    call check_complex_close('lqn agrees lqmn derivative', cd(4), cdm(0, 4), 2.0e-14_dp)

    call lqmn(-1, 2, 0.2_dp, qm, dm, status)
    call check_true('lqmn invalid order status', status == associated_legendre_status_invalid_argument)
    call check_true('lqmn invalid order empty', size(qm) == 0 .and. size(dm) == 0)
    call lqn(-1, 0.2_dp, q, d, status)
    call check_true('lqn invalid degree status', status == associated_legendre_status_invalid_argument)
    call check_true('lqn invalid degree empty', size(q) == 0 .and. size(d) == 0)

    call lqn(2, 1.0_dp, q, d, status)
    call check_true('lqn singular endpoint nan', all(ieee_is_nan(q)) .and. all(ieee_is_nan(d)))

    print *, 'special Legendre-Q tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_complex_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum absolute complex error

        if (actual == expected) return
        if (ieee_is_nan(real(actual, dp)) .or. ieee_is_nan(aimag(actual)) .or. &
                abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_legendre_q
