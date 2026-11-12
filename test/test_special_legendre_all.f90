program test_special_legendre_all
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, &
        ieee_positive_inf, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_special, only : assoc_legendre_p_all, assoc_legendre_p_derivatives, &
        associated_legendre_status_invalid_argument, associated_legendre_status_success, &
        legendre_p, legendre_p_all, legendre_p_derivatives, sph_harm_y_all, &
        sph_legendre_p_all, sph_legendre_p_derivatives
    implicit none

    real(dp), allocatable :: assoc_values(:, :, :)
    complex(dp), allocatable :: cassoc_values(:, :, :)
    complex(dp), allocatable :: harmonic_values(:, :, :, :)
    real(dp), allocatable :: legendre_values(:, :)
    complex(dp), allocatable :: clegendre_values(:, :)
    real(dp), allocatable :: sph_values(:, :, :)
    real(dp) :: assoc_derivatives(0:1)
    real(dp) :: derivatives(0:4)
    real(dp) :: nan_value
    real(dp) :: inf_value
    real(dp) :: sph_derivatives(0:2)
    integer :: degree
    integer :: order
    integer :: status

    call legendre_p_all(4, 0.3_dp, legendre_values, status, diff_n=2)
    call check_true('legendre all status', status == associated_legendre_status_success)
    call check_true('legendre all lower derivative', lbound(legendre_values, 1) == 0)
    call check_true('legendre all lower degree', lbound(legendre_values, 2) == 0)
    call check_true('legendre all shape d', ubound(legendre_values, 1) == 2)
    call check_true('legendre all shape n', ubound(legendre_values, 2) == 4)
    call check_close('legendre P4', legendre_values(0, 4), &
        0.07293750000000002_dp, 2.0e-15_dp)
    call check_close('legendre P4 prime', legendre_values(1, 4), &
        -1.7774999999999999_dp, 3.0e-15_dp)
    call check_close('legendre P4 second', legendre_values(2, 4), &
        -2.775_dp, 4.0e-15_dp)
    call check_close('legendre scalar alias', legendre_p(4, 0.3_dp), &
        legendre_values(0, 4), 1.0e-15_dp)

    call legendre_p_derivatives(4, 0.3_dp, 4, derivatives)
    call check_close('legendre derivative 3', derivatives(3), 31.5_dp, 2.0e-14_dp)
    call check_close('legendre derivative 4', derivatives(4), 105.0_dp, 5.0e-14_dp)

    call legendre_p_all(4, cmplx(0.3_dp, 0.2_dp, dp), clegendre_values, status, &
        diff_n=2)
    call check_complex_close('complex legendre P4', clegendre_values(0, 4), &
        cmplx(0.13543750000000004_dp, -0.3975000000000001_dp, dp), 4.0e-15_dp)
    call check_complex_close('complex legendre P4 prime', clegendre_values(1, 4), &
        cmplx(-2.4075_dp, -0.6950000000000001_dp, dp), 5.0e-15_dp)
    call check_complex_close('complex legendre P4 second', clegendre_values(2, 4), &
        cmplx(-4.875000000000001_dp, 6.300000000000001_dp, dp), 8.0e-15_dp)

    call assoc_legendre_p_all(3, 2, 0.3_dp, assoc_values, status, diff_n=1)
    call check_true('assoc all status', status == associated_legendre_status_success)
    call check_true('assoc all order lower bound', lbound(assoc_values, 3) == -2)
    call check_true('assoc all order upper bound', ubound(assoc_values, 3) == 2)
    do degree = 0, 3
        do order = -2, 2
            if (abs(order) <= degree) then
                call assoc_legendre_p_derivatives(degree, order, 0.3_dp, 1, &
                    assoc_derivatives)
                call check_close('assoc all value consistency', &
                    assoc_values(0, degree, order), assoc_derivatives(0), 2.0e-14_dp)
                call check_close('assoc all derivative consistency', &
                    assoc_values(1, degree, order), assoc_derivatives(1), 5.0e-14_dp)
            else
                call check_close('assoc all out-of-range order value', &
                    assoc_values(0, degree, order), 0.0_dp, 0.0_dp)
                call check_close('assoc all out-of-range order derivative', &
                    assoc_values(1, degree, order), 0.0_dp, 0.0_dp)
            end if
        end do
    end do
    call check_close('assoc P3 -2', assoc_values(0, 3, -2), &
        0.034125_dp, 2.0e-15_dp)
    call check_close('assoc P3 -1 derivative', assoc_values(1, 3, -1), &
        0.3793480752887442_dp, 8.0e-15_dp)
    call check_close('assoc P3 +2', assoc_values(0, 3, 2), &
        4.095_dp, 5.0e-14_dp)

    call assoc_legendre_p_all(3, 2, 2.0_dp, assoc_values, status, &
        branch_cut=3, diff_n=1)
    call check_close('assoc branch3 -2 value', assoc_values(0, 3, -2), &
        0.75_dp, 2.0e-14_dp)
    call check_close('assoc branch3 -1 derivative', assoc_values(1, 3, -1), &
        7.072540797572914_dp, 8.0e-14_dp)
    call check_close('assoc branch3 order0 value', assoc_values(0, 3, 0), &
        17.000000000000004_dp, 8.0e-14_dp)
    call check_close('assoc branch3 +1 value', assoc_values(0, 3, 1), &
        49.363448015713004_dp, 2.0e-13_dp)
    call check_close('assoc branch3 +2 derivative', assoc_values(1, 3, 2), &
        165.0_dp, 4.0e-13_dp)

    call assoc_legendre_p_all(3, 2, cmplx(0.3_dp, 0.2_dp, dp), &
        cassoc_values, status, norm=.true., diff_n=1)
    call check_complex_close('normalized assoc -1', cassoc_values(0, 3, -1), &
        cmplx(-0.5634990034149943_dp, 0.5120156570845191_dp, dp), 2.0e-14_dp)
    call check_complex_close('normalized assoc order0 derivative', &
        cassoc_values(1, 3, 0), &
        cmplx(-2.1046822800603415_dp, 1.6837458240482737_dp, dp), 4.0e-14_dp)
    call check_complex_close('normalized assoc +2 derivative', &
        cassoc_values(1, 3, 2), &
        cmplx(2.177477037766415_dp, -0.9222255689363639_dp, dp), 5.0e-14_dp)

    call sph_legendre_p_all(3, 2, 0.7_dp, sph_values, status, diff_n=2)
    call check_true('sph all status', status == associated_legendre_status_success)
    call check_true('sph all order lower bound', lbound(sph_values, 3) == -2)
    do order = -2, 2
        call sph_legendre_p_derivatives(3, order, 0.7_dp, 2, sph_derivatives)
        call check_close('sph all value consistency', sph_values(0, 3, order), &
            sph_derivatives(0), 8.0e-15_dp)
        call check_close('sph all derivative consistency', sph_values(1, 3, order), &
            sph_derivatives(1), 2.0e-14_dp)
        call check_close('sph all second consistency', sph_values(2, 3, order), &
            sph_derivatives(2), 5.0e-14_dp)
    end do
    call check_close('sph P3 -2', sph_values(0, 3, -2), &
        0.324400748475796_dp, 8.0e-15_dp)
    call check_close('sph P3 -1 derivative', sph_values(1, 3, -1), &
        -0.5500405041508836_dp, 1.5e-14_dp)
    call check_close('sph P3 order0 second', sph_values(2, 3, 0), &
        1.90539619882026_dp, 4.0e-14_dp)

    call sph_harm_y_all(3, 2, 0.7_dp, 1.1_dp, harmonic_values, status, diff_n=2)
    call check_true('harmonic all status', status == associated_legendre_status_success)
    call check_true('harmonic order lower bound', lbound(harmonic_values, 4) == -2)
    call check_complex_close('harmonic +1 value', harmonic_values(0, 0, 3, 1), &
        cmplx(-0.1817853589955138_dp, -0.35716453963284894_dp, dp), 8.0e-15_dp)
    call check_complex_close('harmonic +1 theta', harmonic_values(1, 0, 3, 1), &
        cmplx(0.24949623930980994_dp, 0.49020014563117_dp, dp), 1.5e-14_dp)
    call check_complex_close('harmonic +1 phi', harmonic_values(0, 1, 3, 1), &
        cmplx(0.35716453963284894_dp, -0.1817853589955138_dp, dp), 1.5e-14_dp)
    call check_complex_close('harmonic +1 theta theta', harmonic_values(2, 0, 3, 1), &
        cmplx(1.4471922658780163_dp, 2.8433849802793922_dp, dp), 5.0e-14_dp)
    call check_complex_close('harmonic +1 theta phi', harmonic_values(1, 1, 3, 1), &
        cmplx(-0.49020014563117_dp, 0.24949623930980994_dp, dp), 2.0e-14_dp)
    call check_complex_close('harmonic +1 phi phi', harmonic_values(0, 2, 3, 1), &
        cmplx(0.1817853589955138_dp, 0.35716453963284894_dp, dp), 1.5e-14_dp)
    call check_complex_close('harmonic -2 theta', harmonic_values(1, 0, 3, -2), &
        cmplx(-0.2925117125360703_dp, -0.401859335056979_dp, dp), 2.0e-14_dp)
    call check_complex_close('harmonic -2 phi phi', harmonic_values(0, 2, 3, -2), &
        cmplx(0.7636408116659054_dp, 1.0491073541562577_dp, dp), 3.0e-14_dp)
    call check_complex_close('harmonic unused mixed order', harmonic_values(2, 1, 3, 1), &
        cmplx(0.0_dp, 0.0_dp, dp), 0.0_dp)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    inf_value = ieee_value(0.0_dp, ieee_positive_inf)
    call legendre_p_all(2, nan_value, legendre_values, status, diff_n=2)
    call check_close('nan P0', legendre_values(0, 0), 1.0_dp, 0.0_dp)
    call check_true('nan P1', ieee_is_nan(legendre_values(0, 1)))
    call check_close('nan P1 prime', legendre_values(1, 1), 1.0_dp, 0.0_dp)
    call check_close('nan P2 second', legendre_values(2, 2), 3.0_dp, 0.0_dp)
    call legendre_p_all(2, inf_value, legendre_values, status, diff_n=2)
    call check_true('infinity P2', .not. ieee_is_nan(legendre_values(0, 2)))
    call check_close('infinity P2 second', legendre_values(2, 2), 3.0_dp, 0.0_dp)

    call assoc_legendre_p_all(2, 1, 2.0_dp, assoc_values, status, diff_n=0)
    call check_close('branch2 outside order0 retained', assoc_values(0, 2, 0), &
        5.5_dp, 2.0e-15_dp)
    call check_true('branch2 outside nonzero order nan', &
        ieee_is_nan(assoc_values(0, 2, 1)))

    call sph_legendre_p_all(1, 1, nan_value, sph_values, status, diff_n=1)
    call check_close('nan theta spherical constant', sph_values(0, 0, 0), &
        0.28209479177387814_dp, 3.0e-16_dp)
    call check_close('nan theta constant derivative', sph_values(1, 0, 0), 0.0_dp, 0.0_dp)
    call check_true('nan theta degree1', ieee_is_nan(sph_values(0, 1, 0)))

    call legendre_p_all(-1, 0.3_dp, legendre_values, status)
    call check_true('legendre invalid status', &
        status == associated_legendre_status_invalid_argument)
    call check_true('legendre invalid empty', size(legendre_values) == 0)
    call assoc_legendre_p_all(2, -1, 0.3_dp, assoc_values, status)
    call check_true('assoc invalid status', &
        status == associated_legendre_status_invalid_argument)
    call check_true('assoc invalid empty', size(assoc_values) == 0)
    call sph_harm_y_all(2, 1, 0.7_dp, 1.1_dp, harmonic_values, status, diff_n=-1)
    call check_true('harmonic invalid status', &
        status == associated_legendre_status_invalid_argument)
    call check_true('harmonic invalid empty', size(harmonic_values) == 0)

    print *, 'special Legendre all-family tests passed'

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
        if (complex_is_nan(actual) .or. abs(actual - expected) > tolerance) then
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

    pure function complex_is_nan(z) result(is_nan)
        complex(dp), intent(in) :: z !! complex value to inspect
        logical :: is_nan

        is_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_is_nan

end program test_special_legendre_all
