program test_distribution_infrastructure
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_random, only : rng_seed, rng_state
    use scifort_stats, only : continuous_distribution, distribution_status_invalid_bounds, &
        distribution_status_invalid_support, distribution_status_success, &
        distribution_status_zero_mass, make_distribution, normal_cdf, normal_isf, &
        normal_logcdf, normal_logpdf, normal_logsf, normal_pdf, normal_ppf, normal_sf, &
        rayleigh_cdf, rayleigh_isf, rayleigh_logcdf, rayleigh_logpdf, rayleigh_logsf, &
        rayleigh_pdf, rayleigh_ppf, rayleigh_sf, truncate, truncated_distribution, &
        uniform_cdf, uniform_pdf, uniform_ppf
    use test_distribution_infrastructure_reference
    implicit none

    type(continuous_distribution) :: base
    type(continuous_distribution) :: uniform_base
    type(truncated_distribution) :: dist
    integer :: failures
    integer :: status
    type(rng_state) :: state
    type(rng_state) :: state2
    real(dp), allocatable :: sample1(:)
    real(dp), allocatable :: sample2(:)
    real(dp) :: bounds(2)
    real(dp) :: qnan

    failures = 0
    qnan = quiet_nan(0.0_dp)

    base = make_distribution(normal_pdf_callback, normal_cdf_callback, normal_ppf_callback, &
        [negative_infinity(0.0_dp), positive_infinity(0.0_dp)], &
        parameters=[0.0_dp, 1.0_dp], logpdf=normal_logpdf_callback, &
        logcdf=normal_logcdf_callback, sf=normal_sf_callback, &
        logsf=normal_logsf_callback, isf=normal_isf_callback, status=status)
    call check_int('make normal status', status, distribution_status_success, failures)
    bounds = base%support()
    call check_true('make normal lower support', bounds(1) < -huge(1.0_dp), failures)
    call check_true('make normal upper support', bounds(2) > huge(1.0_dp), failures)
    call check_close('make normal median', base%median(), 0.0_dp, 0.0_dp, failures)
    call check_close('make normal pdf', base%pdf(0.3_dp), normal_pdf(0.3_dp), 0.0_dp, failures)
    call check_close('make normal cdf', base%cdf(0.3_dp), normal_cdf(0.3_dp), 0.0_dp, failures)
    call check_close('make normal sf', base%ccdf(0.3_dp), normal_sf(0.3_dp), 0.0_dp, failures)
    call check_close('make normal icdf', base%icdf(0.7_dp), normal_ppf(0.7_dp), 0.0_dp, failures)
    call check_close('make normal iccdf', base%iccdf(0.7_dp), normal_isf(0.7_dp), 0.0_dp, failures)

    dist = truncate(base, -1.0_dp, 2.0_dp, status)
    call check_int('truncate normal status', status, distribution_status_success, failures)
    call check_array('truncate normal support', dist%support(), [-1.0_dp, 2.0_dp], &
        0.0_dp, failures)
    call check_array('truncate normal pdf', eval_pdf(dist, truncate_x_ref), &
        truncate_normal_pdf_ref, 4.0e-15_dp, failures)
    call check_array('truncate normal logpdf', eval_logpdf(dist, truncate_x_ref(3:5)), &
        truncate_normal_logpdf_ref, 4.0e-15_dp, failures)
    call check_array('truncate normal cdf', eval_cdf(dist, truncate_x_ref), &
        truncate_normal_cdf_ref, 5.0e-15_dp, failures)
    call check_array('truncate normal logcdf', eval_logcdf(dist, truncate_x_ref(3:5)), &
        truncate_normal_logcdf_ref, 5.0e-15_dp, failures)
    call check_array('truncate normal sf', eval_sf(dist, truncate_x_ref), &
        truncate_normal_sf_ref, 5.0e-15_dp, failures)
    call check_array('truncate normal logsf', eval_logsf(dist, truncate_x_ref(3:5)), &
        truncate_normal_logsf_ref, 5.0e-15_dp, failures)
    call check_true('truncate logpdf lower boundary', .not. ieee_is_finite(dist%logpdf(-1.0_dp)), failures)
    call check_true('truncate logcdf lower boundary', .not. ieee_is_finite(dist%logcdf(-1.0_dp)), failures)
    call check_true('truncate logsf upper boundary', .not. ieee_is_finite(dist%logsf(2.0_dp)), failures)
    call check_array('truncate normal ppf', eval_ppf(dist, truncate_p_ref), &
        truncate_normal_ppf_ref, 7.0e-15_dp, failures)
    call check_array('truncate normal isf', eval_isf(dist, truncate_p_ref), &
        truncate_normal_isf_ref, 7.0e-15_dp, failures)
    call check_close('truncate normal ilogcdf', dist%ilogcdf(log(0.1_dp)), &
        truncate_normal_ppf_ref(2), 7.0e-15_dp, failures)
    call check_close('truncate normal ilogccdf', dist%ilogccdf(log(0.1_dp)), &
        truncate_normal_isf_ref(2), 7.0e-15_dp, failures)

    base = make_distribution(normal_pdf_callback, normal_cdf_callback, normal_ppf_callback, &
        [negative_infinity(0.0_dp), positive_infinity(0.0_dp)], &
        parameters=[1.0_dp, 2.0_dp], logpdf=normal_logpdf_callback, &
        logcdf=normal_logcdf_callback, sf=normal_sf_callback, &
        logsf=normal_logsf_callback, isf=normal_isf_callback)
    dist = truncate(base, -2.0_dp, 2.0_dp)
    call check_array('truncate shifted pdf', eval_pdf(dist, truncate_x_ref), &
        truncate_shifted_pdf_ref, 5.0e-15_dp, failures)
    call check_array('truncate shifted cdf', eval_cdf(dist, truncate_x_ref), &
        truncate_shifted_cdf_ref, 3.0e-14_dp, failures)
    call check_array('truncate shifted ppf', eval_ppf(dist, truncate_p_ref), &
        truncate_shifted_ppf_ref, 1.0e-14_dp, failures)

    base = make_distribution(rayleigh_pdf_callback, rayleigh_cdf_callback, &
        rayleigh_ppf_callback, [0.0_dp, positive_infinity(0.0_dp)], &
        parameters=[0.0_dp, 1.0_dp], inclusive=[.true., .true.], &
        logpdf=rayleigh_logpdf_callback, logcdf=rayleigh_logcdf_callback, &
        sf=rayleigh_sf_callback, logsf=rayleigh_logsf_callback, &
        isf=rayleigh_isf_callback)
    dist = truncate(base, 0.5_dp, 3.0_dp)
    call check_array('truncate rayleigh pdf', eval_pdf(dist, truncate_rayleigh_x_ref), &
        truncate_rayleigh_pdf_ref, 6.0e-15_dp, failures)
    call check_array('truncate rayleigh cdf', eval_cdf(dist, truncate_rayleigh_x_ref), &
        truncate_rayleigh_cdf_ref, 6.0e-15_dp, failures)
    call check_array('truncate rayleigh ppf', eval_ppf(dist, truncate_rayleigh_p_ref), &
        truncate_rayleigh_ppf_ref, 8.0e-15_dp, failures)

    ! Exercise make_distribution fallbacks when only PDF/CDF/PPF are supplied.
    uniform_base = make_distribution(uniform_pdf_callback, uniform_cdf_callback, &
        uniform_ppf_callback, [0.0_dp, 1.0_dp], parameters=[0.0_dp, 1.0_dp], &
        inclusive=[.true., .true.])
    dist = truncate(uniform_base, -2.0_dp, 0.8_dp)
    call check_array('truncate uniform support', dist%support(), [0.0_dp, 0.8_dp], &
        0.0_dp, failures)
    call check_array('truncate uniform pdf', eval_pdf(dist, truncate_uniform_x_ref), &
        truncate_uniform_pdf_ref, 3.0e-15_dp, failures)
    call check_array('truncate uniform cdf', eval_cdf(dist, truncate_uniform_x_ref), &
        truncate_uniform_cdf_ref, 3.0e-15_dp, failures)
    call check_close('make fallback logpdf', uniform_base%logpdf(0.4_dp), 0.0_dp, &
        0.0_dp, failures)
    call check_close('make fallback logcdf', uniform_base%logcdf(0.4_dp), log(0.4_dp), &
        2.0e-15_dp, failures)
    call check_close('make fallback sf', uniform_base%sf(0.4_dp), 0.6_dp, &
        2.0e-16_dp, failures)

    call rng_seed(state, 24680)
    call rng_seed(state2, 24680)
    sample1 = dist%sample(12, state)
    sample2 = dist%sample(12, state2)
    call check_array('truncate sample replay', sample1, sample2, 0.0_dp, failures)
    call check_true('truncate sample support', all(sample1 >= 0.0_dp .and. sample1 <= 0.8_dp), &
        failures)

    dist = truncate(uniform_base, 1.0_dp, 2.0_dp, status)
    call check_int('truncate zero mass status', status, distribution_status_zero_mass, failures)
    call check_true('truncate zero mass invalid', ieee_is_nan(dist%pdf(1.5_dp)), failures)
    dist = truncate(uniform_base, 0.8_dp, 0.8_dp, status)
    call check_int('truncate bad bounds status', status, distribution_status_invalid_bounds, failures)
    dist = truncate(uniform_base, qnan, 0.8_dp, status)
    call check_int('truncate nan bounds status', status, distribution_status_invalid_bounds, failures)
    base = make_distribution(uniform_pdf_callback, uniform_cdf_callback, uniform_ppf_callback, &
        [1.0_dp, 1.0_dp], status=status)
    call check_int('make bad support status', status, distribution_status_invalid_support, failures)
    call check_true('make bad support invalid', ieee_is_nan(base%pdf(1.0_dp)), failures)

    if (failures /= 0) then
        print *, 'distribution infrastructure failures:', failures
        error stop 1
    end if
    print *, 'distribution infrastructure tests passed'

contains

    pure function normal_pdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = normal_pdf(x, par(1), par(2))
    end function normal_pdf_callback

    pure function normal_logpdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = normal_logpdf(x, par(1), par(2))
    end function normal_logpdf_callback

    pure function normal_cdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = normal_cdf(x, par(1), par(2))
    end function normal_cdf_callback

    pure function normal_logcdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = normal_logcdf(x, par(1), par(2))
    end function normal_logcdf_callback

    pure function normal_sf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = normal_sf(x, par(1), par(2))
    end function normal_sf_callback

    pure function normal_logsf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = normal_logsf(x, par(1), par(2))
    end function normal_logsf_callback

    pure function normal_ppf_callback(p, par) result(x)
        real(dp), intent(in) :: p !! lower-tail probability
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: x
        x = normal_ppf(p, par(1), par(2))
    end function normal_ppf_callback

    pure function normal_isf_callback(p, par) result(x)
        real(dp), intent(in) :: p !! upper-tail probability
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: x
        x = normal_isf(p, par(1), par(2))
    end function normal_isf_callback

    pure function rayleigh_pdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = rayleigh_pdf(x, par(1), par(2))
    end function rayleigh_pdf_callback

    pure function rayleigh_logpdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = rayleigh_logpdf(x, par(1), par(2))
    end function rayleigh_logpdf_callback

    pure function rayleigh_cdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = rayleigh_cdf(x, par(1), par(2))
    end function rayleigh_cdf_callback

    pure function rayleigh_logcdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = rayleigh_logcdf(x, par(1), par(2))
    end function rayleigh_logcdf_callback

    pure function rayleigh_sf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = rayleigh_sf(x, par(1), par(2))
    end function rayleigh_sf_callback

    pure function rayleigh_logsf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = rayleigh_logsf(x, par(1), par(2))
    end function rayleigh_logsf_callback

    pure function rayleigh_ppf_callback(p, par) result(x)
        real(dp), intent(in) :: p !! lower-tail probability
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: x
        x = rayleigh_ppf(p, par(1), par(2))
    end function rayleigh_ppf_callback

    pure function rayleigh_isf_callback(p, par) result(x)
        real(dp), intent(in) :: p !! upper-tail probability
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: x
        x = rayleigh_isf(p, par(1), par(2))
    end function rayleigh_isf_callback

    pure function uniform_pdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = uniform_pdf(x, par(1), par(2))
    end function uniform_pdf_callback

    pure function uniform_cdf_callback(x, par) result(y)
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: y
        y = uniform_cdf(x, par(1), par(2))
    end function uniform_cdf_callback

    pure function uniform_ppf_callback(p, par) result(x)
        real(dp), intent(in) :: p !! lower-tail probability
        real(dp), intent(in) :: par(:) !! [loc, scale]
        real(dp) :: x
        x = uniform_ppf(p, par(1), par(2))
    end function uniform_ppf_callback

    function eval_pdf(d, values) result(out)
        type(truncated_distribution), intent(in) :: d !! distribution
        real(dp), intent(in) :: values(:) !! arguments
        real(dp) :: out(size(values))
        integer :: i
        do i = 1, size(values)
            out(i) = d%pdf(values(i))
        end do
    end function eval_pdf

    function eval_logpdf(d, values) result(out)
        type(truncated_distribution), intent(in) :: d !! distribution
        real(dp), intent(in) :: values(:) !! arguments
        real(dp) :: out(size(values))
        integer :: i
        do i = 1, size(values)
            out(i) = d%logpdf(values(i))
        end do
    end function eval_logpdf

    function eval_cdf(d, values) result(out)
        type(truncated_distribution), intent(in) :: d !! distribution
        real(dp), intent(in) :: values(:) !! arguments
        real(dp) :: out(size(values))
        integer :: i
        do i = 1, size(values)
            out(i) = d%cdf(values(i))
        end do
    end function eval_cdf

    function eval_logcdf(d, values) result(out)
        type(truncated_distribution), intent(in) :: d !! distribution
        real(dp), intent(in) :: values(:) !! arguments
        real(dp) :: out(size(values))
        integer :: i
        do i = 1, size(values)
            out(i) = d%logcdf(values(i))
        end do
    end function eval_logcdf

    function eval_sf(d, values) result(out)
        type(truncated_distribution), intent(in) :: d !! distribution
        real(dp), intent(in) :: values(:) !! arguments
        real(dp) :: out(size(values))
        integer :: i
        do i = 1, size(values)
            out(i) = d%sf(values(i))
        end do
    end function eval_sf

    function eval_logsf(d, values) result(out)
        type(truncated_distribution), intent(in) :: d !! distribution
        real(dp), intent(in) :: values(:) !! arguments
        real(dp) :: out(size(values))
        integer :: i
        do i = 1, size(values)
            out(i) = d%logsf(values(i))
        end do
    end function eval_logsf

    function eval_ppf(d, values) result(out)
        type(truncated_distribution), intent(in) :: d !! distribution
        real(dp), intent(in) :: values(:) !! probabilities
        real(dp) :: out(size(values))
        integer :: i
        do i = 1, size(values)
            out(i) = d%ppf(values(i))
        end do
    end function eval_ppf

    function eval_isf(d, values) result(out)
        type(truncated_distribution), intent(in) :: d !! distribution
        real(dp), intent(in) :: values(:) !! probabilities
        real(dp) :: out(size(values))
        integer :: i
        do i = 1, size(values)
            out(i) = d%isf(values(i))
        end do
    end function eval_isf

    subroutine check_int(name, actual, expected, count)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual !! computed value
        integer, intent(in) :: expected !! expected value
        integer, intent(inout) :: count !! accumulated failures
        if (actual /= expected) then
            print *, 'FAIL ', trim(name), actual, expected
            count = count + 1
        end if
    end subroutine check_int

    subroutine check_true(name, condition, count)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition expected to be true
        integer, intent(inout) :: count !! accumulated failures
        if (.not. condition) then
            print *, 'FAIL ', trim(name)
            count = count + 1
        end if
    end subroutine check_true

    subroutine check_close(name, actual, expected, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! SciPy reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failures

        if (.not. same_real(actual, expected, tolerance)) then
            print *, 'FAIL ', trim(name), actual, expected
            count = count + 1
        end if
    end subroutine check_close

    subroutine check_array(name, actual, expected, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! SciPy reference values
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failures

        integer :: i

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(name), ' size'
            count = count + 1
            return
        end if
        do i = 1, size(actual)
            if (.not. same_real(actual(i), expected(i), tolerance)) then
                print *, 'FAIL ', trim(name), i, actual(i), expected(i)
                count = count + 1
                return
            end if
        end do
    end subroutine check_array

    pure logical function same_real(actual, expected, tolerance) result(ok)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! allowed absolute difference

        if (ieee_is_nan(expected)) then
            ok = ieee_is_nan(actual)
        else if (.not. ieee_is_finite(expected)) then
            ok = actual == expected
        else
            ok = ieee_is_finite(actual) .and. abs(actual - expected) <= tolerance
        end if
    end function same_real

end program test_distribution_infrastructure
