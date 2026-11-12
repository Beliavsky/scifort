! SPDX-License-Identifier: MIT
program test_signal_partial_fraction
    use scifort_signal, only : dp, residue, residuez, invres, invresz, &
        partial_fraction_result, partial_fraction_tf_result, &
        partial_fraction_success, partial_fraction_invalid_input
    implicit none
    type(partial_fraction_result) :: pf
    type(partial_fraction_tf_result) :: tf
    complex(dp), parameter :: zero = cmplx(0.0_dp,0.0_dp,dp)

    ! Proper, distinct analog poles, ordered by increasing magnitude.
    pf = residue([1.0_dp,3.0_dp],[1.0_dp,3.0_dp,2.0_dp])
    call check(pf%status == partial_fraction_success,'analog distinct status')
    call check(size(pf%residues) == 2,'analog pole count')
    call check(maxval(abs(pf%residues-[cmplx(2.0_dp,0.0_dp,dp), &
        cmplx(-1.0_dp,0.0_dp,dp)])) < 1.0e-9_dp,'analog distinct residues')
    call check(maxval(abs(pf%poles-[cmplx(-1.0_dp,0.0_dp,dp), &
        cmplx(-2.0_dp,0.0_dp,dp)])) < 1.0e-9_dp,'analog distinct poles')
    call check(size(pf%direct) == 0,'analog proper has no polynomial')
    tf = invres(pf%residues,pf%poles,pf%direct)
    call check(tf%status == partial_fraction_success,'analog inverse status')
    call close_poly(tf%b,cmplx([1.0_dp,3.0_dp],0.0_dp,dp),'analog inverse b')
    call close_poly(tf%a,cmplx([1.0_dp,3.0_dp,2.0_dp],0.0_dp,dp),'analog inverse a')

    ! Repeated pole and nontrivial direct polynomial.
    pf = residue([2.0_dp,7.0_dp,9.0_dp],[1.0_dp,-4.0_dp,4.0_dp])
    call check(pf%status == partial_fraction_success,'analog repeated status')
    call check(maxval(abs(pf%residues-cmplx([15.0_dp,31.0_dp],0.0_dp,dp))) < &
        1.0e-5_dp,'analog repeated residues')
    call check(maxval(abs(pf%poles-cmplx([2.0_dp,2.0_dp],0.0_dp,dp))) < &
        1.0e-5_dp,'analog repeated poles')
    call close_poly(pf%direct,[cmplx(2.0_dp,0.0_dp,dp)],'analog direct')
    tf = invres(pf%residues,pf%poles,pf%direct)
    call close_poly(tf%b,cmplx([2.0_dp,7.0_dp,9.0_dp],0.0_dp,dp),'analog repeat inverse b')

    ! Discrete filters use ascending powers of z**(-1).
    pf = residuez([1.0_dp,0.25_dp],[1.0_dp,-1.2_dp,0.32_dp])
    call check(pf%status == partial_fraction_success,'digital distinct status')
    call check(maxval(abs(pf%residues-cmplx([-1.625_dp,2.625_dp],0.0_dp,dp))) < &
        1.0e-8_dp,'digital distinct residues')
    call check(maxval(abs(pf%poles-cmplx([0.4_dp,0.8_dp],0.0_dp,dp))) < &
        1.0e-8_dp,'digital distinct poles')
    tf = invresz(pf%residues,pf%poles,pf%direct)
    call close_poly(tf%b,cmplx([1.0_dp,0.25_dp],0.0_dp,dp),'digital inverse b')
    call close_poly(tf%a,cmplx([1.0_dp,-1.2_dp,0.32_dp],0.0_dp,dp),'digital inverse a')

    pf = residuez([1.0_dp,2.0_dp],[1.0_dp,-1.0_dp,0.25_dp])
    call check(pf%status == partial_fraction_success,'digital repeated status')
    call check(maxval(abs(pf%residues-cmplx([-4.0_dp,5.0_dp],0.0_dp,dp))) < &
        1.0e-5_dp,'digital repeated residues')
    tf = invresz(pf%residues,pf%poles,pf%direct)
    call close_poly(tf%b,cmplx([1.0_dp,2.0_dp],0.0_dp,dp),'digital repeated inverse b')

    ! A direct digital term and a nonzero pole.
    pf = residuez([1.0_dp,2.0_dp],[1.0_dp,-0.5_dp])
    call close_poly(pf%residues,[cmplx(5.0_dp,0.0_dp,dp)],'digital residue')
    call close_poly(pf%direct,[cmplx(-4.0_dp,0.0_dp,dp)],'digital direct')
    tf = invresz(pf%residues,pf%poles,pf%direct)
    call close_poly(tf%b,cmplx([1.0_dp,2.0_dp],0.0_dp,dp),'digital direct inverse')

    ! Complex-valued transfer function with an explicit real direct term.
    tf = invres([cmplx(1.0_dp,-1.0_dp,dp)], &
        [cmplx(-0.5_dp,0.25_dp,dp)],[3.0_dp])
    call check(tf%status == partial_fraction_success,'complex inverse status')
    pf = residue(tf%b,tf%a)
    call check(pf%status == partial_fraction_success,'complex residue status')
    call check(abs(pf%residues(1)-cmplx(1.0_dp,-1.0_dp,dp)) < 1.0e-9_dp, &
        'complex residue value')
    call check(abs(pf%poles(1)-cmplx(-0.5_dp,0.25_dp,dp)) < 1.0e-9_dp, &
        'complex pole value')

    ! Empty/zero numerator and validation.
    pf = residue([0.0_dp,0.0_dp],[1.0_dp,2.0_dp,1.0_dp])
    call check(pf%status == partial_fraction_success,'zero numerator status')
    call check(all(abs(pf%residues) == 0.0_dp),'zero numerator residues')
    pf = residue([1.0_dp],[0.0_dp,0.0_dp])
    call check(pf%status == partial_fraction_invalid_input,'zero denominator rejected')
    pf = residuez([1.0_dp],[0.0_dp,1.0_dp])
    call check(pf%status == partial_fraction_invalid_input,'zero leading digital coefficient rejected')
    pf = residue([1.0_dp],[1.0_dp,2.0_dp],rtype='invalid')
    call check(pf%status == partial_fraction_invalid_input,'invalid representative rejected')
    tf = invres([zero],[cmplx(1.0_dp,0.0_dp,dp)],tol=-1.0_dp)
    call check(tf%status == partial_fraction_invalid_input,'negative tolerance rejected')
    print '(a)', 'PASS: signal partial-fraction decomposition and reconstruction'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion condition
        character(*), intent(in) :: label !! assertion description
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check

    subroutine close_poly(actual,expected,label)
        complex(dp), intent(in) :: actual(:) !! computed coefficients
        complex(dp), intent(in) :: expected(:) !! expected coefficients
        character(*), intent(in) :: label !! assertion description
        call check(size(actual) == size(expected),label//' length')
        if (size(actual) == size(expected)) then
            call check(maxval(abs(actual-expected)) < 1.0e-5_dp,label)
        end if
    end subroutine close_poly
end program test_signal_partial_fraction
