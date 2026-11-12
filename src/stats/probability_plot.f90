! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_probability_plot
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_normal, only : normal_ppf
    use scifort_tukeylambda, only : tukeylambda_ppf
    use scifort_weibull, only : weibull_ppf
    use scifort_student_t, only : t_ppf
    use scifort_loggamma, only : loggamma_ppf
    use scifort_gamma, only : gamma_ppf
    use scifort_logistic, only : logistic_ppf
    use scifort_exponential, only : exponential_ppf
    use scifort_laplace, only : laplace_ppf
    use scifort_cauchy, only : cauchy_ppf
    use scifort_power_transform, only : boxcox, yeojohnson
    implicit none
    private

    public :: probplot
    public :: ppcc_max
    public :: ppcc_plot
    public :: boxcox_normplot
    public :: yeojohnson_normplot
    public :: probability_plot_ppf

    abstract interface
        function probability_plot_ppf(p, params) result(x)
            import :: dp
            real(dp), intent(in) :: p
            real(dp), intent(in) :: params(:)
            real(dp) :: x
        end function probability_plot_ppf
    end interface

contains

    subroutine probplot(x, osm, osr, slope, intercept, r, dist, sparams, fit, ppf)
        real(dp), intent(in) :: x(:) !! sample values
        real(dp), allocatable, intent(out) :: osm(:) !! theoretical quantiles
        real(dp), allocatable, intent(out) :: osr(:) !! ordered responses
        real(dp), intent(out), optional :: slope !! least-squares slope
        real(dp), intent(out), optional :: intercept !! least-squares intercept
        real(dp), intent(out), optional :: r !! probability-plot correlation
        character(len=*), intent(in), optional :: dist !! distribution name; default norm
        real(dp), intent(in), optional :: sparams(:) !! shape(s), then loc and scale
        logical, intent(in), optional :: fit !! compute least-squares fit; default true
        procedure(probability_plot_ppf), optional :: ppf !! custom PPF callback

        real(dp), allocatable :: probs(:)
        real(dp) :: slope_value
        real(dp) :: intercept_value
        real(dp) :: r_value
        logical :: fit_value
        integer :: i

        fit_value = .true.
        if (present(fit)) fit_value = fit

        if (size(x) == 0) then
            allocate(osm(0), osr(0))
            slope_value = quiet_nan(0.0_dp)
            intercept_value = quiet_nan(0.0_dp)
            r_value = 0.0_dp
            call assign_fit_outputs(slope, intercept, r, slope_value, intercept_value, r_value)
            return
        end if

        call filliben_probabilities(size(x), probs)
        allocate(osm(size(x)), osr(size(x)))
        do i = 1, size(x)
            if (present(ppf)) then
                if (present(sparams)) then
                    osm(i) = ppf(probs(i), sparams)
                else
                    osm(i) = ppf(probs(i), [real(dp) ::])
                end if
            else
                osm(i) = named_ppf(probs(i), dist, sparams)
            end if
        end do
        osr = x
        call sort_real_nan_last(osr)

        if (fit_value) then
            call linear_fit(osm, osr, slope_value, intercept_value, r_value)
        else
            slope_value = quiet_nan(0.0_dp)
            intercept_value = quiet_nan(0.0_dp)
            r_value = quiet_nan(0.0_dp)
        end if
        call assign_fit_outputs(slope, intercept, r, slope_value, intercept_value, r_value)
    end subroutine probplot

    function ppcc_max(x, brack, dist, ppf) result(shape_value)
        real(dp), intent(in) :: x(:) !! one-dimensional sample
        real(dp), intent(in), optional :: brack(:) !! two- or three-point starting bracket
        character(len=*), intent(in), optional :: dist !! one-shape distribution name
        procedure(probability_plot_ppf), optional :: ppf !! custom PPF callback
        real(dp) :: shape_value

        real(dp) :: interval(2)
        logical :: fixed_bracket

        if (size(x) == 0) then
            error stop 'scifort: ppcc_max requires a nonempty sample', quiet=.true.
        end if
        interval = [0.0_dp, 1.0_dp]
        fixed_bracket = .false.
        if (present(brack)) then
            select case (size(brack))
            case (2)
                if (brack(1) >= brack(2)) then
                    error stop 'scifort: ppcc_max two-point bracket must be increasing', quiet=.true.
                end if
                interval = brack
            case (3)
                if (.not. (brack(1) < brack(2) .and. brack(2) < brack(3))) then
                    error stop 'scifort: ppcc_max three-point bracket must be increasing', quiet=.true.
                end if
                interval = [brack(1), brack(3)]
                fixed_bracket = .true.
            case default
                error stop 'scifort: ppcc_max bracket must contain two or three values', quiet=.true.
            end select
        end if
        call validate_shape_distribution(dist, ppf)
        shape_value = maximize_ppcc(x, interval, dist, ppf, fixed_bracket)
    end function ppcc_max

    subroutine ppcc_plot(x, a, b, svals, ppcc, dist, n, ppf)
        real(dp), intent(in) :: x(:) !! one-dimensional sample
        real(dp), intent(in) :: a !! lower shape value
        real(dp), intent(in) :: b !! upper shape value
        real(dp), allocatable, intent(out) :: svals(:) !! shape grid
        real(dp), allocatable, intent(out) :: ppcc(:) !! probability-plot correlations
        character(len=*), intent(in), optional :: dist !! one-shape distribution name
        integer, intent(in), optional :: n !! number of grid points; default 80
        procedure(probability_plot_ppf), optional :: ppf !! custom PPF callback

        integer :: i
        integer :: n_value

        if (b <= a) error stop 'scifort: ppcc_plot requires b > a', quiet=.true.
        n_value = 80
        if (present(n)) n_value = n
        if (n_value < 0) error stop 'scifort: ppcc_plot n must be nonnegative', quiet=.true.
        call validate_shape_distribution(dist, ppf)
        call linspace(a, b, n_value, svals)
        allocate(ppcc(n_value))
        if (size(x) == 0) then
            ppcc = 0.0_dp
            return
        end if
        do i = 1, n_value
            ppcc(i) = shape_correlation(x, svals(i), dist, ppf)
        end do
    end subroutine ppcc_plot

    subroutine boxcox_normplot(x, la, lb, lmbdas, ppcc, n)
        real(dp), intent(in) :: x(:) !! positive sample
        real(dp), intent(in) :: la !! lower lambda
        real(dp), intent(in) :: lb !! upper lambda
        real(dp), allocatable, intent(out) :: lmbdas(:) !! lambda grid
        real(dp), allocatable, intent(out) :: ppcc(:) !! normal PPCC values
        integer, intent(in), optional :: n !! number of grid points; default 80

        call normplot_impl(x, la, lb, lmbdas, ppcc, n, .false.)
    end subroutine boxcox_normplot

    subroutine yeojohnson_normplot(x, la, lb, lmbdas, ppcc, n)
        real(dp), intent(in) :: x(:) !! sample
        real(dp), intent(in) :: la !! lower lambda
        real(dp), intent(in) :: lb !! upper lambda
        real(dp), allocatable, intent(out) :: lmbdas(:) !! lambda grid
        real(dp), allocatable, intent(out) :: ppcc(:) !! normal PPCC values
        integer, intent(in), optional :: n !! number of grid points; default 80

        call normplot_impl(x, la, lb, lmbdas, ppcc, n, .true.)
    end subroutine yeojohnson_normplot

    subroutine normplot_impl(x, la, lb, lmbdas, ppcc, n, is_yj)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: la
        real(dp), intent(in) :: lb
        real(dp), allocatable, intent(out) :: lmbdas(:)
        real(dp), allocatable, intent(out) :: ppcc(:)
        integer, intent(in), optional :: n
        logical, intent(in) :: is_yj

        real(dp), allocatable :: transformed(:)
        integer :: i
        integer :: n_value

        if (size(x) == 0) then
            allocate(lmbdas(0), ppcc(0))
            return
        end if
        if (lb <= la) error stop 'scifort: normality-plot upper bound must exceed lower bound', quiet=.true.
        if (.not. is_yj .and. any(x <= 0.0_dp)) then
            error stop 'scifort: boxcox_normplot data must be positive', quiet=.true.
        end if
        n_value = 80
        if (present(n)) n_value = n
        if (n_value < 0) error stop 'scifort: normality-plot n must be nonnegative', quiet=.true.
        call linspace(la, lb, n_value, lmbdas)
        allocate(ppcc(n_value), transformed(size(x)))
        do i = 1, n_value
            if (is_yj) then
                transformed = yeojohnson(x, lmbdas(i))
            else
                transformed = boxcox(x, lmbdas(i))
            end if
            ppcc(i) = normal_probability_correlation(transformed)
        end do
    end subroutine normplot_impl

    function shape_correlation(x, shape, dist, ppf) result(value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: shape
        character(len=*), intent(in), optional :: dist
        procedure(probability_plot_ppf), optional :: ppf
        real(dp) :: value

        real(dp), allocatable :: expected(:)
        real(dp), allocatable :: ordered(:)
        real(dp), allocatable :: probs(:)
        integer :: i

        if (size(x) < 2) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call filliben_probabilities(size(x), probs)
        allocate(expected(size(x)), ordered(size(x)))
        do i = 1, size(x)
            if (present(ppf)) then
                expected(i) = ppf(probs(i), [shape])
            else if (present(dist)) then
                expected(i) = named_ppf(probs(i), dist, [shape])
            else
                expected(i) = named_ppf(probs(i), 'tukeylambda', [shape])
            end if
        end do
        ordered = x
        call sort_real_nan_last(ordered)
        value = correlation(expected, ordered)
    end function shape_correlation

    function normal_probability_correlation(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value

        real(dp), allocatable :: expected(:)
        real(dp), allocatable :: ordered(:)
        real(dp), allocatable :: probs(:)
        integer :: i

        if (size(x) < 2) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call filliben_probabilities(size(x), probs)
        allocate(expected(size(x)), ordered(size(x)))
        do i = 1, size(x)
            expected(i) = normal_ppf(probs(i))
        end do
        ordered = x
        call sort_real_nan_last(ordered)
        value = correlation(expected, ordered)
    end function normal_probability_correlation

    function maximize_ppcc(x, initial, dist, ppf, fixed_bracket) result(best_x)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: initial(2)
        character(len=*), intent(in), optional :: dist
        procedure(probability_plot_ppf), optional :: ppf
        logical, intent(in) :: fixed_bracket
        real(dp) :: best_x

        real(dp) :: a
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        real(dp) :: fc
        real(dp) :: fd
        real(dp) :: left
        real(dp) :: right
        real(dp) :: ra
        real(dp) :: rb
        real(dp) :: rc
        real(dp) :: temp
        real(dp) :: xa
        real(dp) :: xb
        real(dp) :: xc
        integer :: i
        real(dp), parameter :: bracket_gold = 1.618034_dp
        real(dp), parameter :: phi = 0.6180339887498948482_dp
        real(dp), parameter :: tol = 1.0e-11_dp

        if (fixed_bracket) then
            left = initial(1)
            right = initial(2)
        else
            xa = initial(1)
            xb = initial(2)
            ra = safe_shape_correlation(x, xa, dist, ppf)
            rb = safe_shape_correlation(x, xb, dist, ppf)
            if (rb < ra) then
                temp = xa
                xa = xb
                xb = temp
                temp = ra
                ra = rb
                rb = temp
            end if
            xc = xb + bracket_gold * (xb - xa)
            rc = safe_shape_correlation(x, xc, dist, ppf)
            do i = 1, 100
                if (rc <= rb) exit
                xa = xb
                xb = xc
                rb = rc
                xc = xb + bracket_gold * (xb - xa)
                rc = safe_shape_correlation(x, xc, dist, ppf)
            end do
            left = min(xa, xc)
            right = max(xa, xc)
        end if

        a = left
        b = right
        c = b - phi * (b - a)
        d = a + phi * (b - a)
        fc = safe_shape_correlation(x, c, dist, ppf)
        fd = safe_shape_correlation(x, d, dist, ppf)
        do i = 1, 240
            if (abs(b - a) <= tol * (1.0_dp + abs(c) + abs(d))) exit
            if (fc > fd) then
                b = d
                d = c
                fd = fc
                c = b - phi * (b - a)
                fc = safe_shape_correlation(x, c, dist, ppf)
            else
                a = c
                c = d
                fc = fd
                d = a + phi * (b - a)
                fd = safe_shape_correlation(x, d, dist, ppf)
            end if
        end do
        if (fc >= fd) then
            best_x = c
        else
            best_x = d
        end if
    end function maximize_ppcc

    function safe_shape_correlation(x, shape, dist, ppf) result(value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: shape
        character(len=*), intent(in), optional :: dist
        procedure(probability_plot_ppf), optional :: ppf
        real(dp) :: value

        value = shape_correlation(x, shape, dist, ppf)
        if (.not. ieee_is_finite(value)) value = -huge(1.0_dp)
    end function safe_shape_correlation

    function named_ppf(p, dist, sparams) result(x)
        real(dp), intent(in) :: p
        character(len=*), intent(in), optional :: dist
        real(dp), intent(in), optional :: sparams(:)
        real(dp) :: x

        character(len=32) :: name
        real(dp), allocatable :: params(:)
        real(dp) :: loc
        real(dp) :: scale
        real(dp) :: shape

        name = 'norm'
        if (present(dist)) name = lowercase(trim(adjustl(dist)))
        if (present(sparams)) then
            params = sparams
        else
            allocate(params(0))
        end if

        select case (trim(name))
        case ('norm', 'normal')
            call parse_loc_scale(params, 0, loc, scale)
            x = normal_ppf(p, loc, scale)
        case ('logistic')
            call parse_loc_scale(params, 0, loc, scale)
            x = logistic_ppf(p, loc, scale)
        case ('expon', 'exponential')
            call parse_loc_scale(params, 0, loc, scale)
            x = exponential_ppf(p, loc, scale)
        case ('laplace')
            call parse_loc_scale(params, 0, loc, scale)
            x = laplace_ppf(p, loc, scale)
        case ('cauchy')
            call parse_loc_scale(params, 0, loc, scale)
            x = cauchy_ppf(p, loc, scale)
        case ('tukeylambda')
            call parse_one_shape(params, shape, loc, scale)
            x = tukeylambda_ppf(p, shape, loc, scale)
        case ('weibull_min', 'weibull')
            call parse_one_shape(params, shape, loc, scale)
            x = weibull_ppf(p, shape, loc, scale)
        case ('t')
            call parse_one_shape(params, shape, loc, scale)
            x = t_ppf(p, shape, loc, scale)
        case ('loggamma')
            call parse_one_shape(params, shape, loc, scale)
            x = loggamma_ppf(p, shape, loc, scale)
        case ('gamma')
            call parse_one_shape(params, shape, loc, scale)
            x = gamma_ppf(p, shape, loc, scale)
        case default
            error stop 'scifort: unsupported probability-plot distribution', quiet=.true.
        end select
    end function named_ppf

    subroutine validate_shape_distribution(dist, ppf)
        character(len=*), intent(in), optional :: dist
        procedure(probability_plot_ppf), optional :: ppf
        character(len=32) :: name

        if (present(ppf)) return
        name = 'tukeylambda'
        if (present(dist)) name = lowercase(trim(adjustl(dist)))
        select case (trim(name))
        case ('tukeylambda', 'weibull_min', 'weibull', 't', 'loggamma', 'gamma')
        case default
            error stop 'scifort: PPCC distribution must have one supported shape parameter', quiet=.true.
        end select
    end subroutine validate_shape_distribution

    subroutine parse_loc_scale(params, nshapes, loc, scale)
        real(dp), intent(in) :: params(:)
        integer, intent(in) :: nshapes
        real(dp), intent(out) :: loc
        real(dp), intent(out) :: scale

        if (size(params) < nshapes .or. size(params) > nshapes + 2) then
            error stop 'scifort: invalid probability-plot distribution parameter count', quiet=.true.
        end if
        loc = 0.0_dp
        scale = 1.0_dp
        if (size(params) >= nshapes + 1) loc = params(nshapes + 1)
        if (size(params) >= nshapes + 2) scale = params(nshapes + 2)
    end subroutine parse_loc_scale

    subroutine parse_one_shape(params, shape, loc, scale)
        real(dp), intent(in) :: params(:)
        real(dp), intent(out) :: shape
        real(dp), intent(out) :: loc
        real(dp), intent(out) :: scale

        if (size(params) < 1 .or. size(params) > 3) then
            error stop 'scifort: one-shape probability distribution requires shape[,loc,scale]', quiet=.true.
        end if
        shape = params(1)
        call parse_loc_scale(params, 1, loc, scale)
    end subroutine parse_one_shape

    subroutine filliben_probabilities(n, probabilities)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: probabilities(:)
        integer :: i

        if (n < 0) error stop 'scifort: sample size must be nonnegative', quiet=.true.
        allocate(probabilities(n))
        if (n == 0) return
        probabilities(n) = 0.5_dp ** (1.0_dp / real(n, dp))
        probabilities(1) = 1.0_dp - probabilities(n)
        do i = 2, n - 1
            probabilities(i) = (real(i, dp) - 0.3175_dp) / (real(n, dp) + 0.365_dp)
        end do
    end subroutine filliben_probabilities

    subroutine linear_fit(x, y, slope, intercept, r)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: slope
        real(dp), intent(out) :: intercept
        real(dp), intent(out) :: r

        real(dp) :: mx
        real(dp) :: my
        real(dp) :: sxx
        real(dp) :: sxy
        real(dp) :: syy

        if (size(x) /= size(y) .or. size(x) < 2 .or. &
            any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) then
            slope = quiet_nan(0.0_dp)
            intercept = quiet_nan(0.0_dp)
            r = quiet_nan(0.0_dp)
            return
        end if
        mx = sum(x) / real(size(x), dp)
        my = sum(y) / real(size(y), dp)
        sxx = sum((x - mx) * (x - mx))
        syy = sum((y - my) * (y - my))
        sxy = sum((x - mx) * (y - my))
        if (sxx <= 0.0_dp) then
            slope = quiet_nan(0.0_dp)
            intercept = quiet_nan(0.0_dp)
            r = quiet_nan(0.0_dp)
            return
        end if
        slope = sxy / sxx
        intercept = my - slope * mx
        if (syy <= 0.0_dp) then
            r = quiet_nan(0.0_dp)
        else
            r = sxy / sqrt(sxx * syy)
        end if
    end subroutine linear_fit

    function correlation(x, y) result(value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)
        real(dp) :: value
        real(dp) :: slope
        real(dp) :: intercept

        call linear_fit(x, y, slope, intercept, value)
    end function correlation

    subroutine assign_fit_outputs(slope, intercept, r, slope_value, intercept_value, r_value)
        real(dp), intent(out), optional :: slope
        real(dp), intent(out), optional :: intercept
        real(dp), intent(out), optional :: r
        real(dp), intent(in) :: slope_value
        real(dp), intent(in) :: intercept_value
        real(dp), intent(in) :: r_value

        if (present(slope)) slope = slope_value
        if (present(intercept)) intercept = intercept_value
        if (present(r)) r = r_value
    end subroutine assign_fit_outputs

    subroutine linspace(a, b, n, values)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: values(:)
        integer :: i

        allocate(values(n))
        if (n == 0) return
        if (n == 1) then
            values(1) = a
            return
        end if
        do i = 1, n
            values(i) = a + real(i - 1, dp) * (b - a) / real(n - 1, dp)
        end do
    end subroutine linspace

    subroutine sort_real_nan_last(values)
        real(dp), intent(inout) :: values(:)
        integer :: i
        integer :: j
        real(dp) :: key

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (.not. comes_after(values(j), key)) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_real_nan_last

    function comes_after(a, b) result(answer)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: answer

        if (ieee_is_nan(a)) then
            answer = .not. ieee_is_nan(b)
        else if (ieee_is_nan(b)) then
            answer = .false.
        else
            answer = a > b
        end if
    end function comes_after

    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: lower
        integer :: c
        integer :: i

        lower = text
        do i = 1, len(text)
            c = iachar(text(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) lower(i:i) = achar(c + 32)
        end do
    end function lowercase

end module scifort_probability_plot
