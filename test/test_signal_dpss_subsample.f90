! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Check the half-sample DPSS normalization without external dependencies.
program test_signal_dpss_subsample
    use scifort_kinds, only : dp
    use scifort_signal, only : dpss, dpss_result, dpss_status_success, &
        dpss_status_invalid_input
    implicit none
    type(dpss_result) :: sub, base, approx, default_window, invalid
    integer, parameter :: lengths(*) = [2, 3, 4, 5, 6, 7, 8, 15, 16, 25, 32, 33, 64, 65]
    integer :: j, n, k, symmetric, cases
    real(dp) :: nw, ratio, max_error, l2_norm

    cases = 0
    do j = 1, size(lengths)
        n = lengths(j)
        nw = min(2.5_dp, 0.22_dp*real(n,dp))
        k = min(3, n)
        do symmetric = 0, 1
            sub = dpss(n,nw,kmax=k,sym=symmetric == 1,norm='subsample')
            base = dpss(n,nw,kmax=k,sym=symmetric == 1,norm='l2')
            approx = dpss(n,nw,kmax=k,sym=symmetric == 1,norm='approximate')
            if (sub%status /= dpss_status_success) error stop 'subsample status'
            if (base%status /= dpss_status_success) error stop 'baseline status'
            if (approx%status /= dpss_status_success) error stop 'approximate status'
            if (any(shape(sub%tapers) /= [n,k])) error stop 'subsample shape'
            if (maxval(abs(sub%ratios-base%ratios)) > 2.0e-13_dp) &
                error stop 'subsample changed concentration ratios'
            ! Normalization is a single real scale factor shared by all tapers.
            ratio = dot_product(sub%tapers(:,1),base%tapers(:,1)) / &
                dot_product(base%tapers(:,1),base%tapers(:,1))
            if (ratio <= 0.0_dp) error stop 'negative subsample scaling'
            max_error = maxval(abs(sub%tapers-ratio*base%tapers))
            if (max_error > 2.0e-12_dp) error stop 'nonuniform subsample scaling'
            if (symmetric == 1 .and. mod(n,2) == 0) then
                ! An independent direct DFT must reconstruct unit amplitude
                ! halfway between the two central time samples.
                if (abs(midpoint_from_dft(sub%tapers(:,1))-1.0_dp) > 2.0e-12_dp) &
                    error stop 'subsample interpolation peak not unity'
            end if
            if (symmetric == 1 .and. mod(n,2) == 1) then
                ! Odd-length normalization is simply the common maximum.
                l2_norm = maxval(base%tapers)
                if (maxval(abs(sub%tapers-base%tapers/l2_norm)) > 1.0e-12_dp) &
                    error stop 'odd-length maximum normalization'
                if (maxval(abs(sub%tapers-approx%tapers)) > 1.0e-12_dp) &
                    error stop 'odd-length corrections must agree'
            end if
            cases = cases + 1
        end do
    end do
    sub = dpss(64,2.5_dp,norm='subsample')
    base = dpss(64,2.5_dp,kmax=1,norm='subsample')
    if (sub%status /= dpss_status_success) error stop 'singleton norm status'
    if (maxval(abs(sub%tapers-base%tapers)) > 1.0e-13_dp) &
        error stop 'implicit single taper norm'
    default_window = dpss(64,2.5_dp)
    approx = dpss(64,2.5_dp,norm='approximate')
    if (maxval(abs(default_window%tapers-approx%tapers)) > 1.0e-13_dp) &
        error stop 'default normalization changed'
    invalid = dpss(8,2.0_dp,kmax=2,norm='bogus')
    if (invalid%status /= dpss_status_invalid_input) error stop 'invalid norm accepted'
    sub = dpss(1,0.4_dp,norm='subsample')
    if (sub%status /= dpss_status_success .or. sub%tapers(1,1) /= 1.0_dp) &
        error stop 'singleton length'
    print '(a,i0,a)', 'test_signal_dpss_subsample passed (',cases,' cases)'
contains
    real(dp) function midpoint_from_dft(x) result(value)
        real(dp), intent(in) :: x(:) !! symmetric even-length taper to interpolate
        real(dp) :: angle, shift, pi
        complex(dp) :: coefficient, rotated
        integer :: n, j, k

        pi = acos(-1.0_dp)
        n = size(x)
        value = 0.0_dp
        do k = 0, n/2
            coefficient = (0.0_dp,0.0_dp)
            do j = 0, n-1
                angle = -2.0_dp*pi*real(j*k,dp)/real(n,dp)
                coefficient = coefficient + x(j+1)*cmplx(cos(angle),sin(angle),dp)
            end do
            shift = pi*real(k,dp)*(1.0_dp-1.0_dp/real(n,dp))
            rotated = coefficient*cmplx(cos(shift),sin(shift),dp)
            if (k == 0) then
                value = value + real(rotated,dp)
            else
                value = value + 2.0_dp*real(rotated,dp)
            end if
        end do
        value = value / real(n,dp)
    end function midpoint_from_dft
end program test_signal_dpss_subsample
