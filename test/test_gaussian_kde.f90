program test_gaussian_kde
    use scifort_stats, only : gaussian_kde, gaussian_kde_init, &
        gaussian_kde_status_success, gaussian_kde_status_invalid_parameter, &
        gaussian_kde_status_invalid_shape, gaussian_kde_status_singular
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state
    use gaussian_kde_reference
    implicit none

    type(gaussian_kde) :: k1
    type(gaussian_kde) :: k1b
    type(gaussian_kde) :: k2
    type(gaussian_kde) :: km
    type(rng_state) :: state1
    type(rng_state) :: state2
    integer :: info
    integer :: failures
    real(dp) :: x1(5)
    real(dp) :: x1b(4)
    real(dp) :: p1(1,3)
    real(dp) :: x2(2,6)
    real(dp) :: w2(6)
    real(dp) :: p2(2,3)
    real(dp) :: mean2(2)
    real(dp) :: cov2(2,2)
    real(dp) :: inv2(2,2)
    real(dp) :: low2(2)
    real(dp) :: high2(2)
    real(dp), allocatable :: values(:)
    real(dp) :: samples1(2,8)
    real(dp) :: samples2(2,8)

    failures = 0
    x1 = [-2.0_dp, -0.5_dp, 0.25_dp, 1.5_dp, 3.0_dp]
    x1b = [-1.2_dp, 0.1_dp, 0.8_dp, 2.2_dp]
    p1(1,:) = [-1.0_dp, 0.0_dp, 2.0_dp]
    call gaussian_kde_init(k1, x1, info)
    call check_int('1d init', info, gaussian_kde_status_success, failures)
    call check_close('1d factor', k1%factor, ref_1d_factor, 2.0e-14_dp, failures)
    call check_close('1d covariance', k1%covariance(1,1), ref_1d_cov, 2.0e-14_dp, failures)
    call k1%set_bandwidth('silverman', info)
    call check_close('1d silverman factor', k1%factor, ref_1d_silverman_factor, 2.0e-14_dp, failures)
    call check_close('1d silverman cov', k1%covariance(1,1), ref_1d_silverman_cov, 3.0e-14_dp, failures)
    call k1%set_bandwidth('scott', info)
    values = k1%pdf(p1, info)
    call check_int('1d pdf status', info, gaussian_kde_status_success, failures)
    call check_vec('1d pdf', values, ref_1d_pdf, 3.0e-14_dp, failures)
    values = k1%logpdf(p1, info)
    call check_vec('1d logpdf', values, ref_1d_logpdf, 3.0e-14_dp, failures)
    call check_close('1d point pdf', k1%pdf([0.0_dp]), ref_1d_pdf(2), 3.0e-14_dp, failures)
    call check_close('1d box', k1%integrate_box_1d(-0.75_dp, 1.25_dp, info), &
        ref_1d_box, 3.0e-14_dp, failures)
    call check_close('1d gaussian', k1%integrate_gaussian([0.4_dp], reshape([0.7_dp],[1,1]), info), &
        ref_1d_gaussian, 3.0e-14_dp, failures)
    call gaussian_kde_init(k1b, x1b, info)
    call check_close('1d kde product', k1%integrate_kde(k1b, info), ref_1d_kde, 5.0e-14_dp, failures)

    x2(1,:) = [-1.0_dp, 0.0_dp, 0.5_dp, 2.0_dp, 1.0_dp, 2.5_dp]
    x2(2,:) = [0.5_dp, -0.2_dp, 1.5_dp, 1.2_dp, 2.4_dp, 0.7_dp]
    w2 = [0.10_dp, 0.20_dp, 0.25_dp, 0.15_dp, 0.20_dp, 0.10_dp]
    p2(:,1) = [0.0_dp, 0.2_dp]
    p2(:,2) = [1.2_dp, 1.7_dp]
    p2(:,3) = [2.0_dp, 1.0_dp]
    call gaussian_kde_init(k2, x2, info, weights=w2)
    call check_int('2d init', info, gaussian_kde_status_success, failures)
    call check_close('2d neff', k2%neff, ref_2d_neff, 2.0e-14_dp, failures)
    call check_close('2d factor', k2%factor, ref_2d_factor, 2.0e-14_dp, failures)
    call check_mat('2d covariance', k2%covariance, ref_2d_cov, 3.0e-14_dp, failures)
    call k2%inv_cov(inv2, info)
    call check_mat('2d inv cov', inv2, ref_2d_inv_cov, 5.0e-13_dp, failures)
    values = k2%pdf(p2, info)
    call check_vec('2d pdf', values, ref_2d_pdf, 5.0e-14_dp, failures)
    values = k2%logpdf(p2, info)
    call check_vec('2d logpdf', values, ref_2d_logpdf, 5.0e-14_dp, failures)
    mean2 = [0.3_dp, 0.8_dp]
    cov2 = reshape([0.8_dp,0.1_dp,0.1_dp,0.5_dp], [2,2])
    call check_close('2d gaussian', k2%integrate_gaussian(mean2, cov2, info), &
        ref_2d_gaussian, 5.0e-14_dp, failures)
    low2 = [-0.5_dp, 0.0_dp]
    high2 = [1.5_dp, 2.0_dp]
    call rng_seed(state1, 20260922)
    call check_close('2d box', k2%integrate_box(low2, high2, maxpts=400000, state=state1, status=info), &
        ref_2d_box, 3.0e-5_dp, failures)

    call k2%set_bandwidth(0.42_dp, info)
    call check_int('constant bandwidth status', info, gaussian_kde_status_success, failures)
    call check_close('constant factor', k2%factor, ref_constant_factor, 2.0e-14_dp, failures)
    call check_mat('constant covariance', k2%covariance, ref_constant_cov, 3.0e-14_dp, failures)
    call k2%set_bandwidth('scott', info)
    km = k2%marginal([2], info)
    call check_int('marginal status', info, gaussian_kde_status_success, failures)
    call check_close('marginal factor', km%factor, ref_marginal_factor, 2.0e-14_dp, failures)
    values = km%pdf(reshape([0.0_dp,1.0_dp,2.0_dp],[1,3]), info)
    call check_vec('marginal pdf', values, ref_marginal_pdf, 6.0e-14_dp, failures)

    call rng_seed(state1, 918273)
    call rng_seed(state2, 918273)
    call k2%resample(state1, samples1, info)
    call check_int('resample status', info, gaussian_kde_status_success, failures)
    call k2%resample(state2, samples2, info)
    call check_mat('resample deterministic', samples1, samples2, 0.0_dp, failures)

    call invalid_cases(failures)

    if (failures /= 0) then
        write(*,'(a,i0)') 'test_gaussian_kde: FAIL ', failures
        error stop 1
    end if
    write(*,'(a)') 'test_gaussian_kde: PASS'

contains

    subroutine invalid_cases(failures)
        integer, intent(inout) :: failures !! failure counter
        type(gaussian_kde) :: bad
        integer :: info
        real(dp) :: singular_data(2,3)

        singular_data(1,:) = [0.0_dp, 1.0_dp, 2.0_dp]
        singular_data(2,:) = 2.0_dp * singular_data(1,:)
        call gaussian_kde_init(bad, singular_data, info)
        call check_int('singular data', info, gaussian_kde_status_singular, failures)
        call gaussian_kde_init(bad, [1.0_dp], info)
        call check_int('single observation', info, gaussian_kde_status_invalid_shape, failures)
        call gaussian_kde_init(bad, [0.0_dp,1.0_dp,2.0_dp], info, bw_factor=-1.0_dp)
        call check_int('negative bandwidth', info, gaussian_kde_status_invalid_parameter, failures)
    end subroutine invalid_cases

    subroutine check_close(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: got !! observed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: failures !! failure counter

        if (abs(got - expected) > tol) then
            write(*,'(a,2es24.15,es12.3)') trim(name)//' mismatch: ', got, expected, abs(got-expected)
            failures = failures + 1
        end if
    end subroutine check_close

    subroutine check_vec(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: got(:) !! observed vector
        real(dp), intent(in) :: expected(:) !! reference vector
        real(dp), intent(in) :: tol !! maximum absolute tolerance
        integer, intent(inout) :: failures !! failure counter

        if (size(got) /= size(expected) .or. maxval(abs(got-expected)) > tol) then
            write(*,'(a,es12.3)') trim(name)//' max error: ', maxval(abs(got-expected))
            failures = failures + 1
        end if
    end subroutine check_vec

    subroutine check_mat(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: got(:, :) !! observed matrix
        real(dp), intent(in) :: expected(:, :) !! reference matrix
        real(dp), intent(in) :: tol !! maximum absolute tolerance
        integer, intent(inout) :: failures !! failure counter

        if (any(shape(got) /= shape(expected)) .or. maxval(abs(got-expected)) > tol) then
            write(*,'(a,es12.3)') trim(name)//' max error: ', maxval(abs(got-expected))
            failures = failures + 1
        end if
    end subroutine check_mat

    subroutine check_int(name, got, expected, failures)
        character(len=*), intent(in) :: name !! test label
        integer, intent(in) :: got !! observed integer
        integer, intent(in) :: expected !! expected integer
        integer, intent(inout) :: failures !! failure counter

        if (got /= expected) then
            write(*,'(a,2i8)') trim(name)//' mismatch: ', got, expected
            failures = failures + 1
        end if
    end subroutine check_int

end program test_gaussian_kde
