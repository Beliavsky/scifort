! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_aaa_large
    use scifort_kinds, only : dp
    use scifort_constants, only : pi => scifort_pi
    use scifort_interpolate, only : aaa,aaa_type,complex_aaa_type
    use scifort_interpolate_aaa_spectrum, only : aaa_arrowhead_eigenvalues, &
        aaa_spectral_max_terms
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_value,ieee_quiet_nan
    implicit none
    type(aaa_type) :: fit
    type(complex_aaa_type) :: complex_fit
    real(dp) :: x(65),f(65),before,after,nan,t
    complex(dp) :: cx(81),cf(81),q,baseline,value
    complex(dp), allocatable :: poles(:),zeros(:),residues(:),support(:),weights(:)
    complex(dp), allocatable :: alpha(:),beta(:)
    integer :: i,status,removed,checks,old_terms
    checks = 0
    x = [(-1.0_dp+real(i-1,dp)/32.0_dp,i=1,65)]
    f = 1.1_dp+sin(1000*x)+0.3_dp*cos(747*x)
    fit = aaa(x,f,rtol=0.0_dp,max_terms=65,clean_up=.false.,status=status)
    call require(status == 0 .and. fit%valid() .and. fit%converged(),'large real fit')
    call require(size(fit%support_points()) > 39,'large real support count')
    poles = fit%poles(status)
    call require(status == 0 .and. size(poles) > 39,'large real poles')
    zeros = fit%roots(status)
    call require(status == 0 .and. size(zeros) > 39,'large real zeros')
    residues = fit%residues(status)
    call require(status == 0 .and. size(residues) == size(poles),'large real residues')
    call require(all(ieee_is_finite(real(residues,dp))), 'finite real-fit residues')
    before = fit%value(0.12345_dp)
    removed = fit%clean_up(1e-18_dp,status)
    call require(status == 0 .and. removed == 0,'large real cleanup')
    after = fit%value(0.12345_dp)
    call require(abs(after-before) < 3e-12_dp*(1+abs(before)),'cleanup preserves real fit')
    old_terms = size(fit%support_points())
    removed = fit%clean_up(0.1_dp,status)
    call require(status == 0 .and. removed > 0 .and. fit%valid(),'large cleanup removes poles')
    call require(size(fit%support_points()) == old_terms-removed,'cleanup support accounting')
    do i = 1, size(cx)
        t = -1.0_dp+2.0_dp*real(i-1,dp)/real(size(cx)-1,dp)
        cx(i) = cmplx(t,0.2_dp*sin(3*t),dp)
        cf(i) = cmplx(sin(239.0_dp*real(i,dp)),cos(137.0_dp*real(i,dp)),dp)
    end do
    complex_fit = aaa(cx,cf,rtol=0.0_dp,max_terms=81,clean_up=.false.,status=status)
    call require(status == 0 .and. complex_fit%converged(),'large complex fit')
    call require(size(complex_fit%support_points()) > 39,'large complex support count')
    poles = complex_fit%poles(status)
    call require(status == 0 .and. size(poles) > 39,'large complex poles')
    zeros = complex_fit%roots(status)
    call require(status == 0 .and. size(zeros) > 39,'large complex zeros')
    residues = complex_fit%residues(status)
    call require(status == 0 .and. size(residues) == size(poles),'large complex residues')
    q = cmplx(0.13_dp,0.17_dp,dp)
    baseline = complex_fit%value(q)
    removed = complex_fit%clean_up(1e-18_dp,status)
    call require(status == 0 .and. removed == 0,'large complex cleanup')
    value = complex_fit%value(q)
    call require(abs(value-baseline) < 3e-12_dp*(1+abs(baseline)), &
        'cleanup preserves complex fit')
    complex_fit = aaa(cx,cf,rtol=0.0_dp,max_terms=81,status=status)
    call require(status == 0 .and. complex_fit%converged(),'large automatic cleanup')
    call require(complex_fit%cleanup_status() == 0,'stored successful cleanup')

    ! Analytic secular numerator z**63-0.4 gives 63 simple known roots.
    allocate(support(64),weights(64),alpha(65),beta(65))
    do i = 1, 64
        support(i) = exp(cmplx(0.0_dp,2*pi*real(i-1,dp)/64.0_dp,dp))
    end do
    weights = (1-0.4_dp*support)/1.4_dp
    call aaa_arrowhead_eigenvalues(support,weights,alpha,beta,status)
    call require(status == 0,'analytic large pencil')
    call require(count(abs(beta) > 64*epsilon(1.0_dp)*abs(alpha)) == 63, &
        'two homogeneous infinite eigenvalues')
    do i = 1, 65
        if (abs(beta(i)) <= 64*epsilon(1.0_dp)*abs(alpha(i))) cycle
        q = alpha(i)/beta(i)
        call require(abs(q**63-0.4_dp) < 3e-11_dp,'analytic secular polynomial')
    end do
    call aaa_arrowhead_eigenvalues(support,weights(:63),alpha,beta,status)
    call require(status == 1,'pencil shape validation')
    call aaa_arrowhead_eigenvalues(support,weights,alpha(:64),beta,status)
    call require(status == 1,'homogeneous output shape validation')
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    weights(1) = cmplx(nan,0.0_dp,dp)
    call aaa_arrowhead_eigenvalues(support,weights,alpha,beta,status)
    call require(status == 1,'nonfinite weights rejected')
    weights = 3.0_dp
    call aaa_arrowhead_eigenvalues(support,weights,alpha,beta,status)
    call require(status == 1,'normalization bound enforced')
    deallocate(support,weights,alpha,beta)
    allocate(support(aaa_spectral_max_terms+1),weights(aaa_spectral_max_terms+1), &
        alpha(aaa_spectral_max_terms+2),beta(aaa_spectral_max_terms+2))
    support = 0.0_dp
    weights = 1.0_dp
    call aaa_arrowhead_eigenvalues(support,weights,alpha,beta,status)
    call require(status == 1,'larger resource guard before solver allocation')
    print '(a,i0,a)', 'Large AAA spectra: ',checks,' checks passed.'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Assertion that must hold.
        character(len=*), intent(in) :: label !! Failure description.
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',label
            error stop 1
        end if
        checks = checks+1
    end subroutine require
end program test_interpolate_aaa_large
