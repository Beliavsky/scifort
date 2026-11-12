! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_noncentral_distribution_helpers
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_special, only : chndtr, chndtridf, chndtrinc, chndtrix, dp, &
        ncfdtr, ncfdtri, ncfdtridfd, ncfdtridfn, ncfdtrinc, nctdtr, &
        nctdtridf, nctdtrinc, nctdtrit
    implicit none

    real(dp) :: p
    real(dp) :: x

    call check_close('chndtr ordinary', chndtr(5.0_dp, 4.0_dp, 2.0_dp), &
        0.48196384244277046_dp, 3.0e-12_dp)
    call check_close('chndtrix ordinary', chndtrix(0.7_dp, 4.0_dp, 2.0_dp), &
        7.358409860000526_dp, 3.0e-9_dp)
    p = chndtr(5.0_dp, 4.0_dp, 2.0_dp)
    call check_close('chndtridf round trip', chndtridf(5.0_dp, p, 2.0_dp), &
        4.0_dp, 2.0e-10_dp)
    call check_close('chndtrinc round trip', chndtrinc(5.0_dp, 4.0_dp, p), &
        2.0_dp, 2.0e-10_dp)

    call check_close('ncfdtr ordinary', ncfdtr(5.0_dp, 10.0_dp, 2.0_dp, 1.5_dp), &
        0.5640340505917344_dp, 5.0e-12_dp)
    call check_close('ncfdtri ordinary', ncfdtri(5.0_dp, 10.0_dp, 2.0_dp, 0.7_dp), &
        1.9820381269020202_dp, 8.0e-9_dp)
    p = ncfdtr(5.0_dp, 10.0_dp, 2.0_dp, 1.5_dp)
    call check_close('ncfdtridfd round trip', ncfdtridfd(5.0_dp, p, 2.0_dp, 1.5_dp), &
        10.0_dp, 3.0e-9_dp)
    call check_close('ncfdtridfn round trip', ncfdtridfn(p, 10.0_dp, 2.0_dp, 1.5_dp), &
        5.0_dp, 3.0e-9_dp)
    call check_close('ncfdtrinc round trip', ncfdtrinc(5.0_dp, 10.0_dp, p, 1.5_dp), &
        2.0_dp, 3.0e-9_dp)

    p = ncfdtr(20.0_dp, 5.0_dp, 3.0_dp, 0.8_dp)
    x = ncfdtridfd(20.0_dp, p, 3.0_dp, 0.8_dp)
    call check_true('ncfdtridfd smaller-root branch', x < 5.0_dp)
    call check_close('ncfdtridfd smaller-root residual', &
        ncfdtr(20.0_dp, x, 3.0_dp, 0.8_dp), p, 2.0e-13_dp)
    p = ncfdtr(1.5_dp, 30.0_dp, 1.0_dp, 1.2_dp)
    x = ncfdtridfn(p, 30.0_dp, 1.0_dp, 1.2_dp)
    call check_true('ncfdtridfn smaller-root branch', x < 1.5_dp)
    call check_close('ncfdtridfn smaller-root residual', &
        ncfdtr(x, 30.0_dp, 1.0_dp, 1.2_dp), p, 2.0e-13_dp)

    call check_close('nctdtr ordinary', nctdtr(8.0_dp, 1.25_dp, 0.5_dp), &
        0.22372827804528184_dp, 5.0e-10_dp)
    call check_close('nctdtrit ordinary', nctdtrit(8.0_dp, 1.25_dp, 0.7_dp), &
        1.888834779485227_dp, 5.0e-8_dp)
    p = nctdtr(8.0_dp, 1.25_dp, 0.5_dp)
    call check_close('nctdtridf round trip', nctdtridf(p, 1.25_dp, 0.5_dp), &
        8.0_dp, 2.0e-7_dp)
    call check_close('nctdtrinc round trip', nctdtrinc(8.0_dp, p, 0.5_dp), &
        1.25_dp, 2.0e-8_dp)

    call check_close('chndtr lower support', chndtr(0.0_dp, 4.0_dp, 2.0_dp), &
        0.0_dp, 0.0_dp)
    call check_close('ncfdtr lower support', ncfdtr(5.0_dp, 10.0_dp, 2.0_dp, 0.0_dp), &
        0.0_dp, 0.0_dp)
    x = nctdtrit(8.0_dp, 1.25_dp, 0.0_dp)
    call check_true('nctdtrit p=0 negative infinity', .not. ieee_is_finite(x) .and. x < 0.0_dp)
    x = nctdtrit(8.0_dp, 1.25_dp, 1.0_dp)
    call check_true('nctdtrit p=1 positive infinity', .not. ieee_is_finite(x) .and. x > 0.0_dp)
    call check_true('chndtr negative x invalid', ieee_is_nan(chndtr(-1.0_dp, 4.0_dp, 2.0_dp)))
    call check_true('ncfdtr negative f invalid', ieee_is_nan(ncfdtr(5.0_dp, 10.0_dp, 2.0_dp, -1.0_dp)))
    call check_true('nctdtr zero df invalid', ieee_is_nan(nctdtr(0.0_dp, 1.0_dp, 0.0_dp)))

    print *, 'special noncentral-distribution helper tests passed'

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

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_noncentral_distribution_helpers
