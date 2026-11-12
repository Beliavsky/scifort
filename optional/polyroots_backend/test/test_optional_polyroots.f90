! SPDX-License-Identifier: MIT
program test_optional_polyroots
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly, ppoly_type
    use scifort_optional_polyroots, only : polyroots_segment_solver
    implicit none
    type(ppoly_type) :: pp
    real(dp), allocatable :: roots(:), rnative(:)
    real(dp) :: coeff(6)
    ! An explicit quintic with known real roots 0,1,2.
    ! p(t)=t*(t-1)*(t-2)*(t*t+1), with +/-i as nonreal roots.
    coeff=[1.0_dp,-3.0_dp,3.0_dp,-3.0_dp,2.0_dp,0.0_dp]
    pp=ppoly(reshape(coeff,[6,1]),[0.0_dp,3.0_dp],extrapolate=.false.)
    roots=pp%roots(root_solver=polyroots_segment_solver)
    rnative=pp%roots()
    if (size(roots)/=3) error stop 'polyroots callback expected 3 roots'
    if (size(rnative)/=size(roots)) error stop 'native and callback root counts'
    if (maxval(abs(roots-[0.0_dp,1.0_dp,2.0_dp]))>1.0e-8_dp) &
        error stop 'polyroots locations'
    if (maxval(abs(roots-rnative))>1.0e-8_dp) error stop 'polyroots/native disagreement'
    print '(a)', 'optional polyroots-fortran callback PASS'
end program test_optional_polyroots
