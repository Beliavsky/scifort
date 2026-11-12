! SPDX-License-Identifier: MIT
! Native tests for transactional PPoly/BPoly extension and Hermite degrees.
program test_interpolate_piecewise_extend_orders
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly, bpoly, ppoly_type, bpoly_type, &
        complex_ppoly_type, complex_bpoly_type, bpoly_from_derivatives, &
        piecewise_status_invalid_input
    implicit none
    type(ppoly_type) :: pp
    type(bpoly_type) :: bp
    type(complex_ppoly_type) :: cp
    type(complex_bpoly_type) :: cb
    real(dp), allocatable :: x(:), c(:,:,:)
    complex(dp) :: z(2,1)
    real(dp) :: deriv(2,3)
    integer :: status
    real(dp), parameter :: tol=1.0e-11_dp

    ! Basic append, with different power-basis degrees and scalar wrapper.
    pp=ppoly(reshape([2.0_dp,1.0_dp],[2,1]),[0.0_dp,1.0_dp])
    call pp%extend(reshape([0.0_dp,2.0_dp,3.0_dp],[3,1]),[2.0_dp],status)
    if (status/=0) error stop 'PPoly append rejected'
    if (abs(pp%value(0.5_dp)-2.0_dp)>tol) error stop 'PPoly old interval'
    if (abs(pp%value(1.5_dp)-4.0_dp)>tol) error stop 'PPoly new interval'
    c=pp%coefficients()
    if (size(c,1)/=3 .or. size(c,2)/=2) error stop 'PPoly degree promotion'
    if (any(abs(c(:,1,1)-[0.0_dp,2.0_dp,1.0_dp])>tol)) error stop 'PPoly zero padding'
    ! Invalid append must not modify the original, even with 2D input.
    call pp%extend(reshape([1.0_dp],[1,1]),[1.25_dp],status)
    if (status/=piecewise_status_invalid_input) error stop 'PPoly invalid overlap'
    x=pp%knots()
    if (size(x)/=3 .or. x(3)/=2.0_dp) error stop 'PPoly invalid altered knots'
    if (abs(pp%value(1.5_dp)-4.0_dp)>tol) error stop 'PPoly invalid altered coefficients'
    ! Prepend with new higher order; the local origin changes only for the new segment.
    call pp%extend(reshape([0.0_dp,0.0_dp,2.0_dp,-1.0_dp],[4,1]),[-1.0_dp],status)
    if (status/=0 .or. abs(pp%value(-0.5_dp)-0.0_dp)>tol) &
        error stop 'PPoly prepend higher order'

    ! BPoly degree elevation preserves the old linear polynomial.
    bp=bpoly(reshape([1.0_dp,3.0_dp],[2,1]),[0.0_dp,1.0_dp])
    call bp%extend(reshape([3.0_dp,11.0_dp/3.0_dp,13.0_dp/3.0_dp,5.0_dp],[4,1]),[2.0_dp],status)
    if (status/=0) error stop 'BPoly append rejected'
    if (abs(bp%value(0.5_dp)-2.0_dp)>tol .or. abs(bp%value(1.5_dp)-4.0_dp)>tol) &
        error stop 'BPoly extension changed values'
    c=bp%coefficients()
    if (any(abs(c(:,1,1)-[1.0_dp,5.0_dp/3.0_dp,7.0_dp/3.0_dp,3.0_dp])>tol)) &
        error stop 'Bernstein degree elevation'
    call bp%extend(reshape([4.0_dp],[1,1]),[0.5_dp],status)
    if (status/=piecewise_status_invalid_input .or. abs(bp%value(0.5_dp)-2.0_dp)>tol) &
        error stop 'BPoly invalid extension mutated object'

    ! Descending knots and complex-valued power basis.
    z(:,1)=[cmplx(1.0_dp,2.0_dp,dp),cmplx(3.0_dp,-1.0_dp,dp)]
    cp=ppoly(z,[2.0_dp,1.0_dp],status=status)
    call cp%extend(reshape([cmplx(1.0_dp,0.0_dp,dp),cmplx(2.0_dp,0.0_dp,dp), &
        cmplx(3.0_dp,-2.0_dp,dp)],[3,1]),[0.0_dp],status)
    if (status/=0) error stop 'descending complex PPoly extension'
    if (abs(cp%value(1.5_dp)-cmplx(2.5_dp,-2.0_dp,dp))>tol) &
        error stop 'descending complex PPoly original data'
    if (abs(cp%value(0.5_dp)-cmplx(2.25_dp,-2.0_dp,dp))>tol) &
        error stop 'descending complex PPoly appended data'
    cb=bpoly(z,[2.0_dp,1.0_dp],status=status)
    call cb%extend(reshape([cmplx(3.0_dp,-1.0_dp,dp)],[1,1]),[0.0_dp],status)
    if (status/=0 .or. abs(cb%value(1.5_dp)-cmplx(2.0_dp,0.5_dp,dp))>tol) &
        error stop 'complex BPoly extension'

    ! Asymmetric Hermite degree selection: three conditions on the left,
    ! only a function value on the right. SciPy 1.17 rejects some such
    ! mathematically valid cases; verify cubic polynomial reproduction.
    deriv=0.0_dp
    deriv(1,:)=[1.0_dp,2.0_dp,6.0_dp]
    deriv(2,1)=10.0_dp
    bp=bpoly_from_derivatives([0.0_dp,1.0_dp],deriv,[3,1],order=3,status=status)
    if (status/=0) error stop 'asymmetric Hermite degree rejected'
    if (abs(bp%value(0.5_dp)-3.25_dp)>tol) error stop 'asymmetric Hermite interpolation'
    bp=bpoly_from_derivatives([0.0_dp,1.0_dp],deriv,[3,1],orders=[1],status=status)
    if (status/=0 .or. abs(bp%value(0.5_dp)-5.5_dp)>tol) &
        error stop 'Hermite degree one selection'
    bp=bpoly_from_derivatives([0.0_dp,1.0_dp],deriv,[3,1],order=4,status=status)
    if (status/=piecewise_status_invalid_input .or. bp%valid()) &
        error stop 'Hermite impossible polynomial degree accepted'
    bp=bpoly_from_derivatives([0.0_dp,1.0_dp],deriv,[3,1],orders=[3],order=3,status=status)
    if (status/=piecewise_status_invalid_input .or. bp%valid()) &
        error stop 'Hermite conflicting orders accepted'
    print '(a)', 'piecewise extension and Hermite orders native regression PASS'
end program test_interpolate_piecewise_extend_orders
